! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Crystal geometry and elasticity: crystal structures (fcc, bcc, hcp, bct), their slip systems
!> and Schmid tensors in the crystal frame, slip-system families (the active families, in the
!> order given, ROADMAP D48), and the crystal-frame Mandel stiffness matrices.
!>
!> Theory : recoupled_formulation.md §2 (crystal-frame Mandel stiffness [C_c]: cubic; hexagonal
!>          with C33 independent and C66 = (C11 - C12)/2; tetragonal with C66 independent;
!>          Schmid tensor p = sym(s (x) m) as a Mandel 6-vector); the skew part
!>          q = skw(s (x) m) (FEPX arXiv, Dawson & Boyce 2015, Eq. 27) as a skew 3-vector
!>          (conventions §1.4). Crystal frames: cubic axes along the cube edges; hexagonal and
!>          tetragonal with Z_c along c and, for hcp, X_c along a1 = [2 -1 -1 0] (Neper >= 4.10
!>          documentation figure `directionconvention.png`, Neper commit 817c0d41; the same as
!>          v2's Miller-Bravais conversion below).
!> Source : slip systems: port of v2-fixed crys_type_mod2.f90::schmidtensors (tables, order,
!>          Miller and Miller-Bravais to Cartesian conversions, d (x) n), including the bcc
!>          {112}<111> pairing of CHANGES.md C12; families: v2 read_input_cfg.f90 (parse_g_0)
!>          and solve_evp_rstar_mod2.f90 (per-family m) index ranges; family tokens and the
!>          selection and order of the active families: dev/cfg_keywords.md §2.1, ROADMAP D48.
!>
!> Full slip-system tables (v2 order, the "canonical" order; d slip direction, n plane normal,
!> unit vectors, crystal frame), with the canonical families and their tokens:
!>   fcc  12  family 1 {111}<110> 1-12
!>   bcc  48  family 1 {110}<111> 1-12; family 2 {112}<111> 13-24; family 3 {123}<111> 25-48
!>            (family 3: not in v2; added in step 7.3 for Barton et al. (2013), who use all three
!>            families, J. Mech. Phys. Solids 61, 341, §3.1; its table: the 24 {123} planes, each
!>            with the one <111> direction it contains, grouped by direction)
!>   hcp  18  family 1 {0001}<11-20> (basal) 1-3; 2 {10-10}<11-20> (prismatic) 4-6;
!>            3 {10-11}<11-23> (pyramidal <c+a>) 7-18
!>   bct  32  10 families: 1-2, 3-4, 5-6, 7-10, 11-12, 13-16, 17-18, 19-20, 21-24, 25-32
!>            (tokens in bct_tokens)
!> Active families (crystal_init): the slip systems of the given families only, family after
!> family in the order given, each in its table order; default: every family, except bcc
!> {110}<111> only. With the default (or every family in canonical order) the systems are v2's.
!> Miller-Bravais to Cartesian (hcp, c/a = c_over_a):
!>   plane (h k i l):     n = (h, (2k + h)/sqrt3, l/(c/a))
!>   direction [u v t w]: d = (3u/2, (v + u/2) sqrt3, w (c/a))
!> Miller to Cartesian (bct): n = (h, k, l/(c/a)), d = (u, v, w (c/a)).
!>
!> Public parameters:
!>   crystal_fcc, crystal_bcc, crystal_hcp, crystal_bct   integer  structure codes 1..4
!>
!> Public types:
!>   crystal_t   slip geometry of one crystal structure (components documented at the type)
!>
!> Public routines:
!>   crystal_structure_from_name  structure code of "fcc", "bcc", "hcp", "bct" (0 if unknown)
!>   crystal_structure_name       name of a structure code
!>   crystal_init                 build the slip geometry of a structure for the given active
!>                                families (validated input)
!>   crystal_num_families         number of slip families of a structure
!>   crystal_family_token         token of a family ("{111}<110>", ...)
!>   crystal_family_from_token    family of a token (0 if not a token of the structure)
!>   crystal_family_size          number of slip systems of a family
!>   crystal_default_families     default active families of a structure
!>   crystal_loop_normals         loop normals of a loop-normal family token of a cubic
!>                                structure (irradiation, step 7.3)
!>   crystal_system_label         Miller (hcp: Miller-Bravais) label "(plane)[direction]" of a
!>                                system of the full table
!>   stiffness_cubic              crystal-frame Mandel stiffness, cubic (C11, C12, C44)
!>   stiffness_hexagonal          ... hexagonal (C11, C12, C13, C33, C44)
!>   stiffness_tetragonal         ... tetragonal (C11, C12, C13, C33, C44, C66)
module crystal_mod

  use kinds_mod, only: rk
  use constants_mod, only: sqrt3
  use tensor_mod, only: mandel_from_matrix, skew_from_matrix
  use linalg_mod, only: ldlt_factor, ldlt_substitute
  use error_mod, only: fatal_error, internal_error

  implicit none

  private

  public :: crystal_fcc, crystal_bcc, crystal_hcp, crystal_bct
  public :: crystal_t
  public :: crystal_structure_from_name, crystal_structure_name, crystal_init
  public :: crystal_num_families, crystal_family_token, crystal_family_from_token
  public :: crystal_family_size, crystal_default_families, crystal_system_label
  public :: crystal_loop_normals
  public :: stiffness_cubic, stiffness_hexagonal, stiffness_tetragonal, shear_modulus_hill

  ! Structure codes.
  integer, parameter :: crystal_fcc = 1
  integer, parameter :: crystal_bcc = 2
  integer, parameter :: crystal_hcp = 3
  integer, parameter :: crystal_bct = 4

  character(len=3), parameter :: structure_names(4) = ['fcc', 'bcc', 'hcp', 'bct']

  !> Slip geometry of one crystal structure (crystal frame).
  !>
  !> Components:
  !>   structure      integer              structure code
  !>   c_over_a       real(rk)             lattice ratio c/a (hcp, bct; 0 for cubic) [-]
  !>   num_slip       integer              number of slip systems
  !>   num_families   integer              number of active slip-system families
  !>   families       integer(num_families) canonical family (module header) of each active
  !>                                       family, in the order given
  !>   slip_dir_c     real(rk)(3,num_slip) unit slip directions d, crystal frame
  !>   slip_normal_c  real(rk)(3,num_slip) unit slip-plane normals n, crystal frame
  !>   schmid_sym_c   real(rk)(6,num_slip) Mandel 6-vectors of p = sym(d (x) n), crystal frame
  !>   schmid_skew_c  real(rk)(3,num_slip) skew 3-vectors of q = skw(d (x) n), crystal frame
  !>   slip_family    integer(num_slip)    active-family index (1..num_families, the position
  !>                                       in families) of each system
  !>   slip_table     integer(num_slip)    index of each system in the full table of the
  !>                                       structure (v2's numbering)
  type :: crystal_t
    integer :: structure = 0
    real(rk) :: c_over_a = 0.0_rk
    integer :: num_slip = 0
    integer :: num_families = 0
    integer, allocatable :: families(:)
    real(rk), allocatable :: slip_dir_c(:, :)
    real(rk), allocatable :: slip_normal_c(:, :)
    real(rk), allocatable :: schmid_sym_c(:, :)
    real(rk), allocatable :: schmid_skew_c(:, :)
    integer, allocatable :: slip_family(:)
    integer, allocatable :: slip_table(:)
  end type crystal_t

  ! --- v2 slip-system tables (crys_type_mod2.f90::schmidtensors), one vector per column ---

  ! Cubic <110> (fcc directions, bcc {110} normals).
  real(rk), parameter :: cub_110(3, 12) = reshape([ &
    & 0.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, -1.0_rk, 1.0_rk, -1.0_rk, 0.0_rk, &
    & 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 0.0_rk, &
    & 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, &
    & 0.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 0.0_rk], [3, 12])

  ! Cubic <111> (fcc normals, bcc directions).
  real(rk), parameter :: cub_111(3, 12) = reshape([ &
    & 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, &
    & 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, &
    & 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, &
    & 1.0_rk, -1.0_rk, -1.0_rk, 1.0_rk, -1.0_rk, -1.0_rk, 1.0_rk, -1.0_rk, -1.0_rk], [3, 12])

  ! Cubic {112} (bcc {112}<111> normals; column i pairs with cub_111(:, i), CHANGES.md C12).
  real(rk), parameter :: cub_211(3, 12) = reshape([ &
    & -1.0_rk, -1.0_rk, 2.0_rk, -1.0_rk, 2.0_rk, -1.0_rk, 2.0_rk, -1.0_rk, -1.0_rk, &
    & -1.0_rk, -1.0_rk, -2.0_rk, -1.0_rk, 2.0_rk, 1.0_rk, 2.0_rk, -1.0_rk, 1.0_rk, &
    & -1.0_rk, 1.0_rk, 2.0_rk, -1.0_rk, -2.0_rk, -1.0_rk, 2.0_rk, 1.0_rk, -1.0_rk, &
    & -1.0_rk, 1.0_rk, -2.0_rk, -1.0_rk, -2.0_rk, 1.0_rk, 2.0_rk, 1.0_rk, 1.0_rk], [3, 12])

  ! Cubic {123} (bcc {123}<111> normals; column i pairs with cub_123_dir(:, i); not in v2).
  real(rk), parameter :: cub_123(3, 24) = reshape([ &
    & 1.0_rk, 2.0_rk, -3.0_rk, 1.0_rk, -3.0_rk, 2.0_rk, 2.0_rk, 1.0_rk, -3.0_rk, 2.0_rk, -3.0_rk, &
    & 1.0_rk, 3.0_rk, -1.0_rk, -2.0_rk, 3.0_rk, -2.0_rk, -1.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 1.0_rk, &
    & -3.0_rk, -2.0_rk, 2.0_rk, 1.0_rk, 3.0_rk, 2.0_rk, -3.0_rk, -1.0_rk, 3.0_rk, -1.0_rk, 2.0_rk, &
    & 3.0_rk, -2.0_rk, 1.0_rk, 1.0_rk, -2.0_rk, -3.0_rk, 1.0_rk, 3.0_rk, 2.0_rk, 2.0_rk, -1.0_rk, &
    & -3.0_rk, 2.0_rk, 3.0_rk, 1.0_rk, 3.0_rk, 1.0_rk, -2.0_rk, 3.0_rk, 2.0_rk, -1.0_rk, 1.0_rk, &
    & -2.0_rk, 3.0_rk, 1.0_rk, 3.0_rk, -2.0_rk, 2.0_rk, -1.0_rk, 3.0_rk, 2.0_rk, 3.0_rk, -1.0_rk, &
    & 3.0_rk, 1.0_rk, 2.0_rk, 3.0_rk, 2.0_rk, 1.0_rk], [3, 24])
  ! <111> direction of each {123} plane of cub_123 (the one <111> in the plane).
  real(rk), parameter :: cub_123_dir(3, 24) = reshape([ &
    & 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, &
    & 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, &
    & 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, &
    & 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, &
    & 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, &
    & -1.0_rk, -1.0_rk, 1.0_rk, -1.0_rk, -1.0_rk, 1.0_rk, -1.0_rk, -1.0_rk, 1.0_rk, -1.0_rk, &
    & -1.0_rk, 1.0_rk, -1.0_rk, -1.0_rk, 1.0_rk, -1.0_rk, -1.0_rk], [3, 24])

  ! hcp plane normals (h k i l), Miller-Bravais.
  real(rk), parameter :: hcp_n(4, 18) = reshape([ &
    & 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, &
    & 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, -1.0_rk, 0.0_rk, &
    & -1.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, -1.0_rk, 0.0_rk, 0.0_rk, &
    & 1.0_rk, 0.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, -1.0_rk, 1.0_rk, &
    & 0.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, &
    & -1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, &
    & -1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, &
    & 0.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, &
    & 1.0_rk, -1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 0.0_rk, 1.0_rk], [4, 18])

  ! hcp slip directions [u v t w], Miller-Bravais.
  real(rk), parameter :: hcp_d(4, 18) = reshape([ &
    & 2.0_rk, -1.0_rk, -1.0_rk, 0.0_rk, -1.0_rk, 2.0_rk, -1.0_rk, 0.0_rk, &
    & -1.0_rk, -1.0_rk, 2.0_rk, 0.0_rk, 2.0_rk, -1.0_rk, -1.0_rk, 0.0_rk, &
    & -1.0_rk, 2.0_rk, -1.0_rk, 0.0_rk, -1.0_rk, -1.0_rk, 2.0_rk, 0.0_rk, &
    & -2.0_rk, 1.0_rk, 1.0_rk, 3.0_rk, -1.0_rk, -1.0_rk, 2.0_rk, 3.0_rk, &
    & -1.0_rk, -1.0_rk, 2.0_rk, 3.0_rk, 1.0_rk, -2.0_rk, 1.0_rk, 3.0_rk, &
    & 1.0_rk, -2.0_rk, 1.0_rk, 3.0_rk, 2.0_rk, -1.0_rk, -1.0_rk, 3.0_rk, &
    & 2.0_rk, -1.0_rk, -1.0_rk, 3.0_rk, 1.0_rk, 1.0_rk, -2.0_rk, 3.0_rk, &
    & 1.0_rk, 1.0_rk, -2.0_rk, 3.0_rk, -1.0_rk, 2.0_rk, -1.0_rk, 3.0_rk, &
    & -1.0_rk, 2.0_rk, -1.0_rk, 3.0_rk, -2.0_rk, 1.0_rk, 1.0_rk, 3.0_rk], [4, 18])

  ! bct plane normals (h k l), Miller.
  real(rk), parameter :: bct_n(3, 32) = reshape([ &
    & 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, &
    & 1.0_rk, -1.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, &
    & 1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, -1.0_rk, 0.0_rk, &
    & 1.0_rk, -1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, -1.0_rk, 0.0_rk, &
    & 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, &
    & 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, &
    & 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, &
    & 1.0_rk, 0.0_rk, -1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, -1.0_rk, &
    & 1.0_rk, 2.0_rk, 1.0_rk, -1.0_rk, 2.0_rk, 1.0_rk, -1.0_rk, -2.0_rk, 1.0_rk, &
    & 1.0_rk, -2.0_rk, 1.0_rk, 2.0_rk, 1.0_rk, 1.0_rk, -2.0_rk, 1.0_rk, 1.0_rk, &
    & -2.0_rk, -1.0_rk, 1.0_rk, 2.0_rk, -1.0_rk, 1.0_rk], [3, 32])

  ! bct slip directions [u v w], Miller.
  real(rk), parameter :: bct_d(3, 32) = reshape([ &
    & 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, &
    & 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, &
    & 1.0_rk, -1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, 1.0_rk, &
    & -1.0_rk, -1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, &
    & 0.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, &
    & 1.0_rk, 0.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, &
    & 1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, -1.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, -1.0_rk, &
    & 1.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, -1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, &
    & -1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, &
    & -1.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, -1.0_rk, 1.0_rk, &
    & 0.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk], [3, 32])

  ! bct family of each slip system (v2 parse_g_0: 1-2, 3-4, 5-6, 7-10, 11-12, 13-16, 17-18,
  ! 19-20, 21-24, 25-32).
  integer, parameter :: bct_family(32) = [1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, &
    & 7, 7, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10]

  ! Family tokens (dev/cfg_keywords.md §2.1), canonical order; each is the family's name in
  ! Miller (hcp: Miller-Bravais) notation, {plane}<direction>.
  integer, parameter :: token_len = 16
  character(len=token_len), parameter :: fcc_tokens(1) = [character(len=token_len) :: &
    & '{111}<110>']
  character(len=token_len), parameter :: bcc_tokens(3) = [character(len=token_len) :: &
    & '{110}<111>', '{112}<111>', '{123}<111>']
  character(len=token_len), parameter :: hcp_tokens(3) = [character(len=token_len) :: &
    & '{0001}<11-20>', '{10-10}<11-20>', '{10-11}<11-23>']
  character(len=token_len), parameter :: bct_tokens(10) = [character(len=token_len) :: &
    & '{100}<001>', '{110}<001>', '{100}<010>', '{110}<1-11>', '{110}<1-10>', '{100}<011>', &
    & '{001}<100>', '{001}<110>', '{101}<10-1>', '{121}<-101>']

contains

  !> Structure code of a crystal structure name ("fcc", "bcc", "hcp", "bct"), 0 if unknown.
  !>
  !> Arguments:
  !>   name  (in)  character  structure name
  !>   code  (out, result) integer  structure code
  !> Locals:
  !>   i     integer  loop index over structures
  pure function crystal_structure_from_name(name) result(code)

    character(len=*), intent(in) :: name
    integer :: code

    integer :: i

    ! Table lookup.
    code = 0
    do i = 1, size(structure_names)
      if (trim(name) == structure_names(i)) code = i
    end do

  end function crystal_structure_from_name

  !> Name of a structure code ("unknown" if invalid).
  !>
  !> Arguments:
  !>   code  (in)  integer  structure code
  !>   name  (out, result) character(:), allocatable  structure name
  pure function crystal_structure_name(code) result(name)

    integer, intent(in) :: code
    character(len=:), allocatable :: name

    ! Table lookup with a guard.
    if (code >= 1 .and. code <= size(structure_names)) then
      name = structure_names(code)
    else
      name = 'unknown'
    end if

  end function crystal_structure_name

  !> Build the slip geometry of a crystal structure for its active families: unit slip
  !> directions and normals, Schmid tensors p = sym(d (x) n) (Mandel) and q = skw(d (x) n) (skew
  !> 3-vector), family indices. The systems are those of the given families, family after family
  !> in the order given (ROADMAP D48), each family in its table order (module header). Invalid
  !> input (unknown structure, c/a <= 0 for hcp or bct, an empty family list, a family outside
  !> the structure's, a family given twice) is a fatal error.
  !>
  !> Source : port of v2-fixed crys_type_mod2.f90::schmidtensors (module header); family
  !>          selection and order: dev/cfg_keywords.md §2.1, ROADMAP D48.
  !>
  !> Arguments:
  !>   crys           (out) type(crystal_t)  slip geometry
  !>   structure      (in)  integer          structure code
  !>   c_over_a       (in, optional) real(rk)  c/a, required for hcp and bct [-]
  !>   families       (in, optional) integer(:)  active canonical families, in the order wanted
  !>                                           (default: crystal_default_families)
  !> Locals:
  !>   ratio          real(rk)               c/a used for hcp and bct
  !>   active         integer(:), allocatable  the active families
  !>   all_dir        real(rk)(3,:)          full-table slip directions (not normalized)
  !>   all_normal     real(rk)(3,:)          full-table plane normals (not normalized)
  !>   all_family     integer(:)             canonical family of each full-table system
  !>   i, j, k        integer                active family, full-table system, built system
  !>   detail         character(80)          error message
  subroutine crystal_init(crys, structure, c_over_a, families)

    type(crystal_t), intent(out) :: crys
    integer, intent(in) :: structure
    real(rk), intent(in), optional :: c_over_a
    integer, intent(in), optional :: families(:)

    character(len=*), parameter :: routine = 'crystal_init'
    real(rk) :: ratio
    integer, allocatable :: active(:)
    real(rk), allocatable :: all_dir(:, :), all_normal(:, :)
    integer, allocatable :: all_family(:)
    integer :: i, j, k
    character(len=80) :: detail

    ! Structure code.
    if (structure < 1 .or. structure > size(structure_names)) then
      write (detail, '(a,i0)') 'invalid crystal structure code ', structure
      call internal_error(routine, trim(detail))
    end if

    ! c/a is required and must be positive for the non-cubic structures.
    ratio = 0.0_rk
    if (present(c_over_a)) ratio = c_over_a
    if (structure == crystal_hcp .or. structure == crystal_bct) then
      if (.not. (ratio > 0.0_rk)) then
        call internal_error(routine, 'c_over_a must be given and positive for ' // &
          & crystal_structure_name(structure))
      end if
    else
      ratio = 0.0_rk
    end if
    crys%structure = structure
    crys%c_over_a = ratio

    ! Active families: given (checked: not empty, within the structure's, each once) or the
    ! default.
    if (present(families)) then
      if (size(families) == 0) call internal_error(routine, 'no active slip family')
      do i = 1, size(families)
        if (families(i) < 1 .or. families(i) > crystal_num_families(structure)) then
          write (detail, '(a,i0,a,a)') 'invalid slip family ', families(i), ' for ', &
            & crystal_structure_name(structure)
          call internal_error(routine, trim(detail))
        end if
        if (count(families(1:i - 1) == families(i)) > 0) then
          write (detail, '(a,i0,a)') 'slip family ', families(i), ' given twice'
          call internal_error(routine, trim(detail))
        end if
      end do
      allocate (active(size(families)))
      active(:) = families
    else
      ! Allocated first, then assigned (gfortran -Wmaybe-uninitialized, ROADMAP section G).
      allocate (active(size(crystal_default_families(structure))))
      active(:) = crystal_default_families(structure)
    end if

    ! The full table of the structure.
    call full_table(structure, ratio, all_dir, all_normal, all_family)

    ! The systems of the active families, family after family in the given order.
    call allocate_systems(crys, sum([(crystal_family_size(structure, active(i)), &
      & i = 1, size(active))]), size(active))
    crys%families = active
    k = 0
    do i = 1, size(active)
      do j = 1, size(all_family)
        if (all_family(j) /= active(i)) cycle
        k = k + 1
        crys%slip_dir_c(:, k) = all_dir(:, j)
        crys%slip_normal_c(:, k) = all_normal(:, j)
        crys%slip_family(k) = i
        crys%slip_table(k) = j
      end do
    end do

    ! Normalize, then build the symmetric and skew parts of d (x) n.
    do k = 1, crys%num_slip
      crys%slip_dir_c(:, k) = crys%slip_dir_c(:, k)/norm2(crys%slip_dir_c(:, k))
      crys%slip_normal_c(:, k) = crys%slip_normal_c(:, k)/norm2(crys%slip_normal_c(:, k))
      crys%schmid_sym_c(:, k) = mandel_from_matrix(dyad(crys%slip_dir_c(:, k), &
        & crys%slip_normal_c(:, k)))
      crys%schmid_skew_c(:, k) = skew_from_matrix(dyad(crys%slip_dir_c(:, k), &
        & crys%slip_normal_c(:, k)))
    end do

  end subroutine crystal_init

  !> Full slip-system table of a structure, in v2's order, with the canonical family of each
  !> system: Cartesian directions and normals in the crystal frame, not normalized.
  !>
  !> Source : port of v2-fixed crys_type_mod2.f90::schmidtensors (tables, Miller and
  !>          Miller-Bravais to Cartesian conversions; module header).
  !>
  !> Arguments:
  !>   structure  (in)  integer                     structure code (valid)
  !>   ratio      (in)  real(rk)                    c/a (hcp, bct; > 0) [-]
  !>   dir        (out) real(rk)(3,:), allocatable  slip directions d
  !>   normal     (out) real(rk)(3,:), allocatable  plane normals n
  !>   family     (out) integer(:), allocatable     canonical family of each system
  subroutine full_table(structure, ratio, dir, normal, family)

    integer, intent(in) :: structure
    real(rk), intent(in) :: ratio
    real(rk), allocatable, intent(out) :: dir(:, :)
    real(rk), allocatable, intent(out) :: normal(:, :)
    integer, allocatable, intent(out) :: family(:)

    select case (structure)

    case (crystal_fcc)
      ! {111}<110>: d from <110>, n from {111}.
      allocate (dir(3, 12), normal(3, 12), family(12))
      dir = cub_110
      normal = cub_111
      family = 1

    case (crystal_bcc)
      ! {110}<111> (1-12): d from <111>, n from {110}; {112}<111> (13-24): d from <111>, n
      ! from {112}, column i paired with column i (CHANGES.md C12); {123}<111> (25-48): the
      ! paired tables cub_123, cub_123_dir.
      allocate (dir(3, 48), normal(3, 48), family(48))
      dir(:, 1:12) = cub_111
      normal(:, 1:12) = cub_110
      family(1:12) = 1
      dir(:, 13:24) = cub_111
      normal(:, 13:24) = cub_211
      family(13:24) = 2
      dir(:, 25:48) = cub_123_dir
      normal(:, 25:48) = cub_123
      family(25:48) = 3

    case (crystal_hcp)
      ! Miller-Bravais to Cartesian (module header); basal, prismatic, pyramidal <c+a>.
      allocate (dir(3, 18), normal(3, 18), family(18))
      normal(1, :) = hcp_n(1, :)
      normal(2, :) = (2.0_rk*hcp_n(2, :) + hcp_n(1, :))/sqrt3
      normal(3, :) = hcp_n(4, :)/ratio
      dir(1, :) = 1.5_rk*hcp_d(1, :)
      dir(2, :) = (hcp_d(2, :) + 0.5_rk*hcp_d(1, :))*sqrt3
      dir(3, :) = hcp_d(4, :)*ratio
      family(1:3) = 1
      family(4:6) = 2
      family(7:18) = 3

    case default
      ! bct: Miller to Cartesian (module header); 10 families.
      allocate (dir(3, 32), normal(3, 32), family(32))
      normal(1:2, :) = bct_n(1:2, :)
      normal(3, :) = bct_n(3, :)/ratio
      dir(1:2, :) = bct_d(1:2, :)
      dir(3, :) = bct_d(3, :)*ratio
      family = bct_family

    end select

  end subroutine full_table

  !> Allocate the per-system arrays of a crystal_t and set its counts.
  !>
  !> Arguments:
  !>   crys          (inout) type(crystal_t)  crystal being built
  !>   num_slip      (in)    integer          number of slip systems
  !>   num_families  (in)    integer          number of active families
  subroutine allocate_systems(crys, num_slip, num_families)

    type(crystal_t), intent(inout) :: crys
    integer, intent(in) :: num_slip
    integer, intent(in) :: num_families

    crys%num_slip = num_slip
    crys%num_families = num_families
    allocate (crys%families(num_families))
    allocate (crys%slip_dir_c(3, num_slip), crys%slip_normal_c(3, num_slip), &
      & crys%schmid_sym_c(6, num_slip), crys%schmid_skew_c(3, num_slip), &
      & crys%slip_family(num_slip), crys%slip_table(num_slip))

  end subroutine allocate_systems

  !> Number of slip families of a structure (0 for an invalid code).
  !>
  !> Arguments:
  !>   structure  (in)  integer  structure code
  !>   num        (out, result) integer  number of families
  pure function crystal_num_families(structure) result(num)

    integer, intent(in) :: structure
    integer :: num

    select case (structure)
    case (crystal_fcc)
      num = size(fcc_tokens)
    case (crystal_bcc)
      num = size(bcc_tokens)
    case (crystal_hcp)
      num = size(hcp_tokens)
    case (crystal_bct)
      num = size(bct_tokens)
    case default
      num = 0
    end select

  end function crystal_num_families

  !> Token of a family of a structure, e.g. "{111}<110>" (empty for an invalid structure or
  !> family).
  !>
  !> Arguments:
  !>   structure  (in)  integer  structure code
  !>   family     (in)  integer  canonical family
  !>   token      (out, result) character(:), allocatable  the token
  pure function crystal_family_token(structure, family) result(token)

    integer, intent(in) :: structure
    integer, intent(in) :: family
    character(len=:), allocatable :: token

    ! Guard, then table lookup.
    token = ''
    if (family < 1 .or. family > crystal_num_families(structure)) return
    select case (structure)
    case (crystal_fcc)
      token = trim(fcc_tokens(family))
    case (crystal_bcc)
      token = trim(bcc_tokens(family))
    case (crystal_hcp)
      token = trim(hcp_tokens(family))
    case default
      token = trim(bct_tokens(family))
    end select

  end function crystal_family_token

  !> Family of a token for a structure, matched exactly (ROADMAP D48); 0 if the token is not
  !> one of the structure's.
  !>
  !> Arguments:
  !>   structure  (in)  integer    structure code
  !>   token      (in)  character  token, e.g. "{112}<111>"
  !>   family     (out, result) integer  canonical family, or 0
  !> Locals:
  !>   i          integer          family index
  pure function crystal_family_from_token(structure, token) result(family)

    integer, intent(in) :: structure
    character(len=*), intent(in) :: token
    integer :: family

    integer :: i

    family = 0
    do i = 1, crystal_num_families(structure)
      if (crystal_family_token(structure, i) == token) family = i
    end do

  end function crystal_family_from_token

  !> Loop normals of a loop-normal family (irradiation softening, step 7.3, D61,
  !> addendum_irradiation.md §5): the token "<111>", "<100>" or "<110>" of an fcc or bcc
  !> structure gives the 4, 3 or 6 cubic directions of that family, the +/- pairs counted once
  !> (unit vectors, crystal frame). Any other token or structure gives no normals (hcp and bct
  !> loops are not defined yet).
  !>
  !> Arguments:
  !>   structure  (in)  integer                 structure code
  !>   token      (in)  character(*)            the family token
  !>   normals_c  (out) real(rk)(3, :), allocatable  the normals (size 0 when not a family)
  !>
  !> Locals:
  !>   i  integer  column
  pure subroutine crystal_loop_normals(structure, token, normals_c)

    integer, intent(in) :: structure
    character(len=*), intent(in) :: token
    real(rk), allocatable, intent(out) :: normals_c(:, :)

    integer :: i

    ! Cubic structures only; the directions of each family, one per +/- pair.
    allocate (normals_c(3, 0))
    if (structure /= crystal_fcc .and. structure /= crystal_bcc) return
    select case (token)
    case ('<111>')
      normals_c = cub_111(:, [1, 4, 7, 10])
    case ('<100>')
      normals_c = reshape([1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
        & 1.0_rk], [3, 3])
    case ('<110>')
      normals_c = reshape([1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, -1.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, &
        & 1.0_rk, 1.0_rk, 0.0_rk, -1.0_rk, 0.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 1.0_rk, -1.0_rk], &
        & [3, 6])
    end select
    do i = 1, size(normals_c, 2)
      normals_c(:, i) = normals_c(:, i)/norm2(normals_c(:, i))
    end do

  end subroutine crystal_loop_normals

  !> Number of slip systems of a family of a structure (0 for an invalid structure or family).
  !>
  !> Arguments:
  !>   structure  (in)  integer  structure code
  !>   family     (in)  integer  canonical family
  !>   num        (out, result) integer  number of systems
  pure function crystal_family_size(structure, family) result(num)

    integer, intent(in) :: structure
    integer, intent(in) :: family
    integer :: num

    ! Guard, then the family ranges of the module header.
    num = 0
    if (family < 1 .or. family > crystal_num_families(structure)) return
    select case (structure)
    case (crystal_fcc)
      num = 12
    case (crystal_bcc)
      num = merge(24, 12, family == 3)
    case (crystal_hcp)
      num = merge(12, 3, family == 3)
    case default
      num = count(bct_family == family)
    end select

  end function crystal_family_size

  !> Default active families of a structure (dev/cfg_keywords.md §2.1): every family in
  !> canonical order, except bcc: {110}<111> only. Empty for an invalid code.
  !>
  !> Arguments:
  !>   structure  (in)  integer  structure code
  !>   families   (out, result) integer(:), allocatable  canonical families
  !> Locals:
  !>   i          integer  family index
  pure function crystal_default_families(structure) result(families)

    integer, intent(in) :: structure
    integer, allocatable :: families(:)

    integer :: i

    if (structure == crystal_bcc) then
      families = [1]
    else
      families = [(i, i = 1, crystal_num_families(structure))]
    end if

  end function crystal_default_families

  !> Label of a system of the full table of a structure in Miller (hcp: Miller-Bravais)
  !> indices, "(plane)[direction]", e.g. "(111)[01-1]" (fcc system 1); empty for an invalid
  !> structure or system. The indices are those of v2's tables (module header).
  !>
  !> Arguments:
  !>   structure  (in)  integer  structure code
  !>   system     (in)  integer  full-table system index (v2's numbering)
  !>   label      (out, result) character(:), allocatable  the label
  !> Locals:
  !>   plane, dir  real(rk)(:), allocatable  Miller or Miller-Bravais indices
  pure function crystal_system_label(structure, system) result(label)

    integer, intent(in) :: structure
    integer, intent(in) :: system
    character(len=:), allocatable :: label

    real(rk), allocatable :: plane(:), dir(:)

    ! The indices of the system in its table.
    label = ''
    select case (structure)
    case (crystal_fcc)
      if (system < 1 .or. system > 12) return
      plane = cub_111(:, system)
      dir = cub_110(:, system)
    case (crystal_bcc)
      if (system < 1 .or. system > 48) return
      if (system <= 12) then
        plane = cub_110(:, system)
        dir = cub_111(:, system)
      else if (system <= 24) then
        plane = cub_211(:, system - 12)
        dir = cub_111(:, system - 12)
      else
        plane = cub_123(:, system - 24)
        dir = cub_123_dir(:, system - 24)
      end if
    case (crystal_hcp)
      if (system < 1 .or. system > 18) return
      plane = hcp_n(:, system)
      dir = hcp_d(:, system)
    case (crystal_bct)
      if (system < 1 .or. system > 32) return
      plane = bct_n(:, system)
      dir = bct_d(:, system)
    case default
      return
    end select

    ! "(" plane ")[" direction "]", indices written one after the other.
    label = '(' // indices_text(plane) // ')[' // indices_text(dir) // ']'

  end function crystal_system_label

  !> Miller indices written one after the other ("1-10" for 1, -1, 0).
  !>
  !> Arguments:
  !>   values  (in)  real(rk)(:)  the indices (whole numbers)
  !>   text    (out, result) character(:), allocatable  the text
  !> Locals:
  !>   i       integer        index
  !>   buf     character(8)   formatting buffer
  pure function indices_text(values) result(text)

    real(rk), intent(in) :: values(:)
    character(len=:), allocatable :: text

    integer :: i
    character(len=8) :: buf

    text = ''
    do i = 1, size(values)
      write (buf, '(i0)') nint(values(i))
      text = text // trim(buf)
    end do

  end function indices_text

  !> Dyadic product a (x) b (3x3 matrix with entries a_i b_j).
  !>
  !> Arguments:
  !>   vec_a, vec_b  (in)  real(rk)(3)    vectors
  !>   mat           (out, result) real(rk)(3,3)  a (x) b
  pure function dyad(vec_a, vec_b) result(mat)

    real(rk), intent(in) :: vec_a(3)
    real(rk), intent(in) :: vec_b(3)
    real(rk) :: mat(3, 3)

    mat = spread(vec_a, 2, 3)*spread(vec_b, 1, 3)

  end function dyad

  !> Crystal-frame Mandel stiffness of a cubic crystal.
  !>
  !> Theory : recoupled_formulation.md §2: normal block [[C11,C12,C12],[C12,C11,C12],
  !>          [C12,C12,C11]], shear diagonal (2 C44, 2 C44, 2 C44) (order 23 13 12).
  !>
  !> Arguments:
  !>   c11, c12, c44  (in)  real(rk)  Voigt elastic constants [MPa]
  !>   stiff_el_c     (out, result) real(rk)(6,6)  Mandel stiffness, crystal frame [MPa]
  pure function stiffness_cubic(c11, c12, c44) result(stiff_el_c)

    real(rk), intent(in) :: c11, c12, c44
    real(rk) :: stiff_el_c(6, 6)

    ! Cubic is tetragonal with C13 = C12, C33 = C11, C66 = C44.
    stiff_el_c = stiffness_tetragonal(c11, c12, c12, c11, c44, c44)

  end function stiffness_cubic

  !> Crystal-frame Mandel stiffness of a hexagonal crystal (c-axis = crystal 3), with C33 an
  !> independent constant (ROADMAP D3) and C66 = (C11 - C12)/2.
  !>
  !> Theory : recoupled_formulation.md §2.
  !>
  !> Arguments:
  !>   c11, c12, c13, c33, c44  (in)  real(rk)  Voigt elastic constants [MPa]
  !>   stiff_el_c     (out, result) real(rk)(6,6)  Mandel stiffness, crystal frame [MPa]
  pure function stiffness_hexagonal(c11, c12, c13, c33, c44) result(stiff_el_c)

    real(rk), intent(in) :: c11, c12, c13, c33, c44
    real(rk) :: stiff_el_c(6, 6)

    ! Tetragonal form with the transversely isotropic C66.
    stiff_el_c = stiffness_tetragonal(c11, c12, c13, c33, c44, 0.5_rk*(c11 - c12))

  end function stiffness_hexagonal

  !> Crystal-frame Mandel stiffness of a tetragonal crystal (c-axis = crystal 3), with C33 and
  !> C66 independent constants.
  !>
  !> Theory : recoupled_formulation.md §2: normal block [[C11,C12,C13],[C12,C11,C13],
  !>          [C13,C13,C33]], shear diagonal (2 C44, 2 C44, 2 C66) (order 23 13 12).
  !>
  !> Arguments:
  !>   c11, c12, c13, c33, c44, c66  (in)  real(rk)  Voigt elastic constants [MPa]
  !>   stiff_el_c     (out, result) real(rk)(6,6)  Mandel stiffness, crystal frame [MPa]
  pure function stiffness_tetragonal(c11, c12, c13, c33, c44, c66) result(stiff_el_c)

    real(rk), intent(in) :: c11, c12, c13, c33, c44, c66
    real(rk) :: stiff_el_c(6, 6)

    ! Normal block.
    stiff_el_c = 0.0_rk
    stiff_el_c(1, 1:3) = [c11, c12, c13]
    stiff_el_c(2, 1:3) = [c12, c11, c13]
    stiff_el_c(3, 1:3) = [c13, c13, c33]

    ! Shear block: Mandel factor 2 on the Voigt shear constants.
    stiff_el_c(4, 4) = 2.0_rk*c44
    stiff_el_c(5, 5) = 2.0_rk*c44
    stiff_el_c(6, 6) = 2.0_rk*c66

  end function stiffness_tetragonal

  !> Voigt-Reuss-Hill average shear modulus of a crystal's stiffness (the matrix shear modulus
  !> G_m of precipitate bowing, step 7.1).
  !>
  !> Theory : addendum_hardening_modifiers.md §5 (D59): Hill (1952), Proc. Phys. Soc. A 65:349,
  !>          general formulas in Voigt notation (engineering shear strains):
  !>            G_V = [(C11+C22+C33) - (C12+C13+C23) + 3(C44+C55+C66)]/15,
  !>            G_R = 15/[4(S11+S22+S33) - 4(S12+S13+S23) + 3(S44+S55+S66)], S = C^-1,
  !>            G = (G_V + G_R)/2.
  !>          With Mandel matrices (factor sqrt(2) per shear index) the Voigt entries are
  !>          C_V(4..6, 4..6) = C_M/2 on the shear diagonal and S_V(4..6, 4..6) = 2 S_M.
  !> Source : replaces v2-fixed crys_type_mod1.f90::crys_initparams' per-structure formulas
  !>          (equal for cubic and bct, C13; equal for hcp with v2's C33 = C11 + C12 - C13).
  !>
  !> Arguments:
  !>   stiff_el_c  (in)  real(rk)(6, 6)  Mandel stiffness, crystal frame [MPa]
  !>   g_hill      (out, result) real(rk)  Voigt-Reuss-Hill shear modulus [MPa]; 0 if the
  !>                                       stiffness is singular
  !> Locals:
  !>   compl       real(rk)(6, 6)  Mandel compliance S_M [1/MPa]
  !>   fact        real(rk)(6, 6)  LDL^T factors of the stiffness
  !>   g_voigt, g_reuss  real(rk)  Voigt and Reuss bounds [MPa]
  !>   i, status   integer         column, solve status
  pure function shear_modulus_hill(stiff_el_c) result(g_hill)

    real(rk), intent(in) :: stiff_el_c(6, 6)
    real(rk) :: g_hill

    real(rk) :: compl(6, 6), fact(6, 6), g_voigt, g_reuss
    integer :: i, status

    ! Voigt bound from the stiffness (Mandel shear diagonal = 2 C44 etc.).
    g_voigt = ((stiff_el_c(1, 1) + stiff_el_c(2, 2) + stiff_el_c(3, 3)) - (stiff_el_c(1, 2) + &
      & stiff_el_c(1, 3) + stiff_el_c(2, 3)) + 1.5_rk*(stiff_el_c(4, 4) + stiff_el_c(5, 5) + &
      & stiff_el_c(6, 6)))/15.0_rk

    ! Compliance: one factorization, then column by column.
    g_hill = 0.0_rk
    fact = stiff_el_c
    call ldlt_factor(fact, status)
    if (status /= 0) return
    do i = 1, 6
      compl(:, i) = 0.0_rk
      compl(i, i) = 1.0_rk
      call ldlt_substitute(fact, compl(:, i))
    end do

    ! Reuss bound from the compliance (Voigt S44 = 2 S_M44 etc.), and the average.
    g_reuss = 15.0_rk/(4.0_rk*(compl(1, 1) + compl(2, 2) + compl(3, 3)) - 4.0_rk*(compl(1, 2) + &
      & compl(1, 3) + compl(2, 3)) + 6.0_rk*(compl(4, 4) + compl(5, 5) + compl(6, 6)))
    g_hill = 0.5_rk*(g_voigt + g_reuss)

  end function shear_modulus_hill

end module crystal_mod
