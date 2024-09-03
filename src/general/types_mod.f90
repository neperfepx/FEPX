! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module types_mod

! Module to define types

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use crys_type_mod3
  use parallel_mod, only: par_max

  implicit none

  public

  type printing_type

    logical :: print_coo
    logical :: print_disp
    logical :: print_vel
    logical :: print_ori
    logical :: print_crss
    logical :: print_stress_eq
    logical :: print_stress
    logical :: print_strain
    logical :: print_strain_el
    logical :: print_strain_pl
    logical :: print_sliprate
    logical :: print_rss
    logical :: print_defrate_eq
    logical :: print_defrate_pl_eq
    logical :: print_strain_eq
    logical :: print_strain_el_eq
    logical :: print_strain_pl_eq
    logical :: print_velgrad
    logical :: print_defrate_pl
    logical :: print_spinrate
    logical :: print_slip
    logical :: print_restart
    logical :: print_forces
    logical :: print_conv
    logical :: print_work
    logical :: print_work_pl
    logical :: print_defrate
    logical :: print_workrate
    logical :: print_workrate_pl
    logical :: print_rotrate
    logical :: print_rotrate_spin
    logical :: print_rotrate_slip

    logical :: restart
    character(len=50) :: restart_file_handling
    integer :: restart_initial_step

    logical :: legacy_printing = .false.

    integer, public :: &
      & ori_u = 1, &
      & coo_u = 2, &
      & vel_u = 3, &
      & stress_u = 4, &
      & strain_el_u = 5, &
      & strain_eq_u = 7,&
      & defrate_eq_u = 8, &
      & sliprate_u = 9, &
      & rss_u = 10, &
      & stats_u = 17, &
      & log_u = 18, &
      & debug_u = 20, &    !  debugging file
      & defrate_pl_eq_u = 21, &
      & strain_pl_eq_u = 22, &
      & conv_u = 23, &
      & bcs_iter_1_u = 24, &
      & bcs_iter_2_u = 25, &
      & bcs_iter_log_u = 26, &
      & velgrad_u = 33, &
      & defrate_pl_u = 34, &
      & spinrate_u = 35, &
      & slip_u = 36, &
      & crss_u = 37, &
      & stress_eq_u = 38, &
      & dot_sim_u = 39, &
      & work_u = 40, &
      & work_pl_u = 41, &
      & defrate_u = 42, &
      & workrate_u = 44, &
      & workrate_pl_u = 45, &
      & strain_el_eq_u = 46, &
      & strain_pl_u = 47, &
      & strain_u = 48, &
      & disp_u = 49, &
      & rotrate_u = 50, &
      & rotrate_spin_u = 51, &
      & rotrate_slip_u = 52, &
      & force_u = 60

  end type printing_type

  type trace

    ! Buffer : data to be received

    real(rk), allocatable :: buffer(:)
    integer, allocatable :: copies(:)
    integer, allocatable :: copy_address_1(:)
    integer, allocatable :: copy_address_2(:)
    integer, allocatable :: node_address(:)
    integer, allocatable :: locals(:)
    integer, allocatable :: n_to_move(:)
    integer, allocatable :: nodes_to_move(:)
    integer, allocatable :: min_n_to_move(:)
    integer, allocatable :: node_trace(:)
    integer, allocatable :: elt_trace(:)
    integer, allocatable :: min_ptr(:)
    integer, allocatable :: copy_ptr(:)
    integer, allocatable :: node_ptr(:)
    integer, allocatable :: req(:)

    integer :: dim1_sub
    integer :: dim1_sup
    integer :: dim2_sub
    integer :: dim2_sup
    integer :: elt_sub
    integer :: elt_sup
    integer :: nlocals
    integer :: nreceived
    integer :: nelt_traces
    integer :: nnode_traces
    integer :: ntotal
    integer :: setup

    ! num_procs: Number of processes
    ! myid: Rank of this process
    ! gs_comm: Communicator for gather/scatter

    integer :: num_procs
    integer :: myid
    integer :: gs_comm

  end type trace

  type userdef_mpc_type

    character(len=256) :: nset
    integer :: nbdof
    character(len=256), allocatable :: dir(:)
    real(rk), allocatable :: coeff(:)
    real(rk), allocatable :: offset(:)

  end type userdef_mpc_type

  type imposed_strain_rate_state_type

    integer :: nb_comp = 0
    character(len=256) :: comp(6) = "00"
    real(rk) :: val(6) = 0.0d0

  end type imposed_strain_rate_state_type

  !> Table for periodic relation direction
  type pv_table_type

    logical :: flag_label(12) = .false.
    integer :: label(12) = 0
    integer :: primary_drive(12) = 0
    integer :: secondary_drive(12) = 0

  end type pv_table_type

  type loading_options_type

    !> can be uniaxial_strain_target, uniaxial_load_target,
    !> triaxial_constant_strain_rate triaxial_constant_load_rate
    character(len = 50) :: def_control_by

    !> number of steps
    integer :: num_steps

    character(len = 50) :: target

    !> target strain
    real(rk), allocatable :: target_time(:)

    !> target strain
    real(rk), allocatable :: target_strain(:)

    !> target load along up to 3 directions
    real(rk), allocatable :: target_load(:, :)

    !> number of increments of a step
    integer,  allocatable :: step_num_incr(:)

    !> print at end of step
    logical,  allocatable :: step_print(:)

    !> can be general uniaxial_grip uniaxial_symmetry uniaxial_minimal triaxial
    character(len = 50) :: boundary_conditions

    !> strain rate
    real(rk) :: strain_rate

    integer :: num_bcs

    character(len=256) :: bc_type(1024)
    character(len=256) :: bc_var(1024)
    character(len=256) :: bc_nset(1024)
    character(len=2)   :: bc_dir(1024)
    real(rk)           :: bc_vel(1024)

    ! should be merged with strain_rate (or velocity) somehow
    integer :: number_of_strain_rate_jumps
    real(rk), allocatable :: strain_rate_jump(:, :)

    ! load/triaxial-related ----------------------------------------------------

    !> min time step
    real(rk), allocatable :: step_dt_min(:)
    !> max time step
    real(rk), allocatable :: step_dt_max(:)

    !> min strain step
    real(rk), allocatable :: step_dstrain_min(:)
    !> max strain step
    real(rk), allocatable :: step_dstrain_max(:)

    real(rk) :: dwell_max_strain_incr                         ! copied to loading
    real(rk), allocatable :: step_target_time_incr(:)   ! ?

    real(rk) :: load_rate                               ! not used in sim
    integer :: number_of_load_rate_jumps                ! not used in sim
    integer :: number_of_dwell_episodes                 ! copied to loading
    real(rk), allocatable :: load_rate_jump(:, :)       ! not used in sim
    real(rk), allocatable :: dwell_episode(:, :)        ! not used in sim

    !> loading direction: 1, 2 or 3
    integer  :: loading_direction
    ! --------------------------------------------------------------------------

    ! parsing helpers
    integer :: curr_strain_jump = 0
    integer :: curr_target = 0
    integer :: curr_load_rate_jump = 0
    integer :: curr_dwell_episode = 0

    integer :: nbsdof = 0
    integer :: nb_tot_dof = 0       ! total number of dof (probably redundant with an other variable)
    integer :: nb_primary_dof = 0   ! number of primary node dofs
    integer :: nb_secondary_dof = 0 ! number of secondary node dofs
    logical :: mpc_status = .false.
    integer :: num_mpcs
    type(userdef_mpc_type) :: general_mpc (1024)

    character(len = 50) :: periodicity = "xyz" ! full periodicity as default input
    type(imposed_strain_rate_state_type) :: imposed_strain_rate_state

  end type loading_options_type

  type loading_type

    ! General loading properties -----------------------------------------------

    !> can be uniaxial_strain_target, uniaxial_load_target,
    !> triaxial_constant_strain_rate triaxial_constant_load_rate
    character(len = 50) :: def_control_by

    !> index of the loading face, rarely used, somewhat redundant with loading_direction
    integer :: loading_face
    !> loading direction: 1, 2 or 3
    integer  :: loading_direction
    
    !> gage length
    real(rk) :: gage_length
    real(rk) :: gage_length_alldir(3)
    real(rk), allocatable :: gage_length_array(:)
    ! --------------------------------------------------------------------------

    ! Step definition ----------------------------------------------------------

    !> number of steps
    integer :: num_steps

    ! All variables below defined the steps, are defined 1:num_steps and won't
    ! move during simulation

    !> target displacement
    real(rk), allocatable :: target_disp(:)

    !> target loads
    real(rk), allocatable :: target_load(:, :)

    ! strain_rate does not exist (but should) as it is defined as one value in
    ! loading_options (same for all steps) and it is used to initialize bcs_vel

    !> load rate
    real(rk), allocatable :: load_rate(:)

    !> target (cumulated) time of a step
    real(rk), allocatable :: target_time(:)

    !> step velocity; used only as a proxy, or in multiaxial
    real(rk), allocatable :: step_velocity(:)

    ! control direction, signed, used in triaxial
    integer, allocatable :: control_dir(:)

    !> min time step allowed at a step
    real(rk), allocatable :: step_dt_min(:)
    !> max time step allowed at a step
    real(rk), allocatable :: step_dt_max(:)

    !> print at end of step
    logical,  allocatable :: step_print(:)

    ! --------------------------------------------------------------------------

    ! velocity arrays, apply to all steps, indexed dof_sub:dof_sup
    ! the bcs_vel_defined values of bcs_vel are pushed to results%vel
    logical, allocatable :: bcs_vel_defined(:)
    real(rk), allocatable :: bcs_vel(:)

    ! triaxial-related ---------------------------------------------------------

    real(rk), allocatable :: target_sign(:)
    real(rk), allocatable :: vel_factor(:)
    real(rk), allocatable :: step_target_time_incr(:)   ! ?

    ! --------------------------------------------------------------------------

    ! mpc variables
    logical :: mpc_status = .false.
    integer :: nb_primary_dof
    integer :: nb_secondary_dof
    real(rk), allocatable :: coeff_ps(:)
    real(rk), allocatable :: offset_ps(:)
    !> connectivity matrix to primary dofs
    real(rk), allocatable :: conn_mpc(:,:)
    !> multiplicity of conn_mpc
    integer, allocatable :: count_mpc(:)
    !> mask if a dof is primary (=1) or secondary (=0)
    real(rk), allocatable :: mask_mpc(:)
    integer, allocatable :: ps_dof(:)
    integer, allocatable :: primary_dof(:)
    type(trace) :: mpc_trace

    logical :: pbc_status = .false.

    !> periodic boundary conditions variables
    integer, allocatable :: primary_drive(:)    ! primary node of the driven relation
    integer, allocatable :: secondary_drive(:)  ! secondary node of the driven relation
    integer, allocatable :: label_pv(:)         ! label of the periodic vector
    integer, allocatable :: label_pv_drive(:)   ! label of the periodic vector of the driven relation
    integer, allocatable :: label_sg(:)         ! sign of label_pv*label_pv_drive
    integer, allocatable :: imposed_state(:)    ! flag for imposed state (1 = true, 0 = false), easier with int

    !> ebe PBC variables (as output of part_gather, easier with real)
    real(rk), allocatable :: primary_drive_ebe(:,:)    ! primary node of the driven relation connectivity matrix
    real(rk), allocatable :: secondary_drive_ebe(:,:)  ! secondary node of the driven relation conectivity matrix
    real(rk), allocatable :: label_pv_ebe(:,:)    ! label pv connectivity matrix
    real(rk), allocatable :: label_pv_drive_ebe(:,:)    ! label pv connectivity matrix
    real(rk), allocatable :: label_sg_ebe(:,:)    ! label sg sonnectivity matrix
    real(rk), allocatable :: imposed_state_ebe(:,:)  ! imposed state conectivity matrix
    real(rk), allocatable :: offset_ps_ebe(:,:)  ! imposed state conectivity matrix
    real(rk), allocatable :: mask_drive(:,:)  ! mask for the drive node term in periodic relations
   
    !> temp velocity array, used only for visco-plastic solver 
    logical, allocatable :: bcs_vel_defined_isovp(:)
    real(rk), allocatable :: bcs_vel_isovp(:)
    
    type(trace) :: primary_drive_trace
    type(trace) :: secondary_drive_trace
    type(trace) :: label_pv_trace
    type(trace) :: imposed_state_trace

    type(pv_table_type) :: pv_table

    ! --------------------------------------------------------------------------

    ! helpers
    real(rk) :: prev_dt

    !> used for necking test and dtime control
    real(rk), dimension(3) :: curr_load
    real(rk), dimension(3) :: prev_load

    logical :: step_complete = .true. ! helps initial printing
    logical :: all_steps_complete

    integer :: curr_step

    ! Variables that should typically vanish -----------------------------------

    ! dwell time, used only in multiaxial
    ! should be removed and replaced by a zero-velocity step
    real(rk), allocatable :: dwell_time(:)
    real(rk) :: dwell_max_strain_incr

    real(rk) :: prev_disp
    real(rk) :: curr_disp

    ! --------------------------------------------------------------------------

  end type loading_type

  type crys_type

    ! Required crystal parameters for all phases.
    character(len=4) :: structure
    real(rk) :: m
    real(rk) :: gammadot_0
    real(rk) :: h_0
    real(rk), allocatable :: g_0(:,:)
    real(rk), allocatable :: g_0_tmp(:)
    logical :: g_0_112 = .false.
    real(rk) :: g_s0
    real(rk) :: n
    real(rk) :: c11
    real(rk) :: c12
    real(rk) :: c13 ! hcp and bct only
    real(rk) :: c44
    real(rk) :: c66 ! bct only

    real(rk) :: kelas(5)
    real(rk) :: keinv(5)
    real(rk) :: elas_coeffs(5)

    real(rk) :: c_over_a ! hcp and bct only

    real(rk) :: m_prime
    real(rk) :: gammadot_s0

    real(rk) :: cyclic_a
    real(rk) :: cyclic_c

    character(len=255) :: hardening

    logical :: anisotropic
    logical :: cyclic
    logical :: saturation_evolution
    ! Precipitation hardening parameters
    logical :: precipitation
    logical :: precipitation_cutting
    logical :: precipitation_bowing

    real(rk), allocatable :: interaction_matrix_parameters (:)
    integer(rk) :: interaction_matrix_parameters_num

    real(rk) :: a_p
    real(rk) :: f_p
    real(rk) :: b_p
    real(rk) :: r_p
    real(rk) :: b_m
    real(rk) :: c_p
    real(rk) :: precipitation_strength

    logical :: use_aniso_m
    real(rk) :: aniso_m (10) ! hardwiring as a temporary patch

    integer :: num_vals ! Number of input strength values
    integer :: num_g_0_vals_112 ! Number of input strength values

    real(rk) :: bulk_shear_mod
    real(rk) :: bulk_mod


    real(rk) :: t_min

    real(rk), allocatable :: fcc_h1(:, :)
    real(rk), allocatable :: fcc_h2(:, :)
    real(rk), allocatable :: fcc_h3(:, :)
    real(rk), allocatable :: fcc_h4(:, :)

    real(rk), allocatable :: bcc_h1(:, :)
    real(rk), allocatable :: bcc_112_h1(:, :)
    real(rk) :: bcc_112_h2(2, 2), bcc_112_h3(2, 2), bcc_112_h4(2, 2), bcc_112_h5(2, 2), bcc_112_h6(2, 2)
    real(rk) :: bcc_h2(2, 2), bcc_h3(2, 2), bcc_h4(2, 2), bcc_h5(2, 2), bcc_h6(2, 2)
    real(rk), allocatable :: hcp_h1(:, :)
    real(rk) :: hcp_h2(2, 2), hcp_h3(2, 2), hcp_h4(2, 2), hcp_h5(2, 2), &
        & hcp_h6(2, 2), hcp_h7(2, 2), hcp_vert(3)
    real(rk), allocatable :: bct_h1(:, :)
    real(rk) :: bct_h2(2, 2), bct_h3(2, 2), bct_h4(4, 4), bct_h5(2, 2), &
        & bct_h6(4, 4), bct_h7(2, 2), bct_h8(2, 2), bct_h9(4, 4), bct_h10(8, 8)

    ! Hardening-related parameters.
    integer :: numslip
    integer :: numvertices

    real(rk), pointer :: schmid_3x3(:, :, :)
    real(rk), pointer :: vertices(:, :)
    real(rk), pointer :: vertices3x3(:, :, :)

  end type crys_type

  type exec_type

    type(trace) :: dof_trace
    type(trace) :: node_trace

    logical :: restart
    real(rk) :: clock_start
    integer :: auto_time = 0

    ! Conjugate Gradient (cg)
    integer :: cg_max_iters = 16000
    real(rk) :: cg_tol = 1.0d-8

    ! Global Nonlinear Iteration (NL)
    integer :: nl_max_iters = 50
    real(rk) :: nl_tol_strict = 5.0d-4
    real(rk) :: nl_tol_loose = 5.0d-4
    real(rk) :: nl_tol_min = 1.0d-10

    ! Newton-Raphson
    real(rk) :: nr_tol_switch_ref = 1.0d-2
    real(rk) :: nr_tol_conv = 0.2d0
    integer :: nr_slope_start = 2

    ! Single Crystal/ State Calculations (sx)
    integer :: sx_max_iters_state = 100
    integer :: sx_max_iters_newton = 100
    real(rk) :: sx_tol = 1.0d-4

    integer  :: max_bc_iter = 10
    integer :: max_incr = 50000
    real(rk) :: max_total_time = 12000.0d0
    real(rk) :: max_strain = 0.2d0
    real(rk) :: max_eqstrain = 0.2d0
    real(rk) :: load_tol_rel = 0.0d0
    real(rk) :: load_tol_abs = 0.1d0
    integer :: max_iter_hard_limit = 10
    real(rk) :: dtime_factor = 1.001d0
    logical :: check_necking = .false.

    real(rk) :: toler_state = 1.0d-5
    real(rk) :: toler_hard = 1.0d-4, ls_cutoff = 0.001d0

    ! helpers during execution

    !> number of NR iterations
    integer :: num_nr_iters

    !> NR attempted during increment
    logical :: nr_attempted

    !> NR slowly converging
    logical :: nr_isslow

    !> number of SA iterations
    integer :: num_sa_iters

    real(rk) :: nr_tol_switch
    real(rk) :: nr_conv
    real(rk) :: nr_conv_prev

    logical :: iter_converged

    character(len=2) :: itmethod

    !> minimum fraction of the control velocity by which the secondary and
    ! tertiary surface velocities are perturbed during boundary condition iterations
    real(rk) :: min_pert_frac

  end type exec_type

  type nset_type

    character(len=50) :: nset_label
    integer :: num_nset_nodes
    integer, allocatable :: nset_nodes(:)

  end type nset_type
  
  type periodicity_type

    integer :: primary              ! primary node
    integer :: secondary            ! secondary node
    integer, dimension(3) :: pvect  ! coordinates of the periodicity vector
    integer :: pvect_label          ! label of the periodicity vector

  end type periodicity_type

  type surface_section

    ! type: Currently a number indicating nodes per element
    ! nsel: Number of surface elements in this section
    ! semin to semax: Range of surface element numbers
    ! conn: Connectivity
    ! conn3d: dof connectivity
    ! econn: Elemental connectivity to be used for gather/scatter stress
    ! elt:
    ! tr: The trace data structure for gather/scatter operations
    ! etr: Elemental trace for stress

    integer :: type, semin, semax
    integer, pointer :: conn(:, :)
    integer, pointer :: conn3d(:, :)
    integer, pointer :: econn(:, :)
    real(rk), pointer :: coo(:, :)
    integer, pointer :: elt(:)
    type(trace) :: tr
    type(trace) :: etr

  end type surface_section

  !> Input mesh
  type mesh_type

    !> global partitioning information
    integer, allocatable :: global_info(:, :)

    !> global number of nodes
    integer:: num_nodes = 0
    !> global number of cells
    integer:: num_cell = 0
    !> global number of element sets
    integer:: num_elsets = 0
    !> global number of elements
    integer:: num_elts = 0

    integer :: num_phases = 1

    !> flag for periodicity
    logical :: periodic_mesh = .false.

    integer, allocatable :: elt_phase(:)

    !> element nodes
    integer, allocatable :: elt_nodes(:,:)

    !> element dofs
    integer, allocatable :: elt_dofs(:,:)

    !> multiplicity
    real(rk), allocatable :: g_ones(:)

    !> initial node coordinates (one row), indexed dof_sub:suf_sup
    real(rk), allocatable :: coo(:)

    !> initial element orientations, indexed elt_sub:elt_sup
    real(rk), allocatable :: ori(:,:,:)

    !> initial element crss, indexed elt_sub:elt_sup  (optional)
    real(rk), allocatable :: crss(:,:,:)
    !> initial element saturation strength, indexed elt_sub:elt_sup  (optional)
    real(rk), allocatable :: sat_str(:,:)

    !> orientation parameterization (rodrigues, etc.)
    character(len=50) :: orientation_parameterization

    !> orientation convention (active or passive)
    character(len=50) :: orientation_convention

    
    ! Optional input flags

    logical :: g_s_from_file = .false.
    logical :: g_0_from_file = .false.

    !> crss file hardening
    logical :: crss_defined = .false.
    !> saturation strength file hardening
    logical :: sat_str_defined = .false.

    integer :: num_fasets

    character(len=50), allocatable :: faset_labels(:)

    type(surface_section), allocatable :: fasets(:)

    !> start id of the elements in the mesh file, only used during parsing
    integer:: elt_startid = -1

    ! numslip: Max number of slip systems
    integer :: maxnumslip = 0

    type(nset_type), allocatable :: nsets(:)
    
    type(periodicity_type), allocatable :: periodicity(:)

    integer :: num_periodicity = 0

  end type mesh_type

  type results_type

    !> Results defined on dofs -------------------------------------------------

    !> dof velocities
    real(rk), allocatable :: vel(:)

    !> dof coordinates
    real(rk), allocatable :: coo(:)

    !> dof forces
    real(rk), allocatable :: force(:)

    !> Results defined on quad points ------------------------------------------

    !> orientations (rotation matrix), in reference frame
    real(rk), allocatable :: ori(:, :, :, :)

    !> lattice rotation (rotation matrix) from undeformed state (mesh%ori)
    real(rk), allocatable :: rstar(:, :, :, :)

    !> lattice rotation (rotation matrix) during time step (usual)
    real(rk), allocatable :: d_rstar(:, :, :, :)

    real(rk), allocatable :: d_rstar_spin(:, :, :, :)
    real(rk), allocatable :: d_rstar_slip(:, :, :, :)

    !> slip rates
    real(rk), allocatable :: sliprate(:, :, :)

    !> resolved shear stresses (dual of sliprate)
    real(rk), allocatable :: rss(:, :, :)

    !> resolved shear stresses
    real(rk), allocatable :: crss(:, :, :)

    !> accumulated slips
    real(rk), allocatable :: acmslip(:, :, :)

    !> volumetric lattice strain, scalar
    real(rk), allocatable :: e_elas_kk_bar(:, :)

    !> deviatoric kirchhoff stress - 5-vec, in crystal frame
    real(rk), allocatable :: sig_vec(:, :, :)

    !> keinv * sig_vec
    real(rk), allocatable :: e_vec(:, :, :)

    !>
    real(rk), allocatable :: sig_kk(:, :)

    !> velocity gradient, in reference frame
    real(rk), allocatable :: velgrad(:, :, :, :)

    !> equivalent deformation rate
    !> has some usage in-code (not only post-processing)
    real(rk), allocatable :: defrate_eq(:, :)

    !> mean/volumetric part of the symmetric part of the vel gradient
    real(rk), allocatable :: d_kk(:, :)

    !> deviatoric part of the symmetric part of the vel gradient, in reference frame
    !> matrix form
    real(rk), allocatable :: d(:, :, :, :)
    !> vector form
    real(rk), allocatable :: d_vec(:, :, :)

    !> skew part of the vel gradient, in reference frame
    !> matrix form
    real(rk), allocatable :: w(:, :, :, :)
    !> vector form
    real(rk), allocatable :: w_vec(:, :, :)
    !> vector form, crystal frame
    real(rk), allocatable :: w_vec_lat(:, :, :)

    !> elemental plastic deformation rate tensor (5-vector), in crystal frame
    !> rename as defrate_pl?
    real(rk), allocatable :: dp_hat(:, :, :)

    !> elemental plastic spin rate tensor (3-vector), in crystal frame
    !> same as w/w_vec, but in crystal frame
    !> rename as spinrate?
    real(rk), allocatable :: wp_hat(:, :, :)

    ! Results defined only on the central quad point ---------------------------
    ! (only used for post-processing) ------------------------------------------
    ! They can be changed into arrays on quadrature points as necessary --------

    !> elastic strain tensors, in crystal frame, 6-vector, i.e. upper triangle
    real(rk), allocatable :: elas_tot6(:, :, :)

    !> elemental deformation rate, in reference frame
    real(rk), allocatable :: defrate(:, :, :)

    !> elemental Cauchy stress, in reference frame
    real(rk), allocatable :: stress(:, :, :)

    !> elemental strain, in reference frame
    real(rk), allocatable :: strain(:, :, :)

    !> elemental plastic strain, in crystal frame, time integration of dp_hat
    real(rk), allocatable :: strain_pl(:, :)

    !> elemental work
    real(rk), allocatable :: work(:)

    !> elemental plastic work
    real(rk), allocatable :: work_pl(:)

    !> elemental work rate
    real(rk), allocatable :: workrate(:)

    !> elemental plastic work rate
    real(rk), allocatable :: workrate_pl(:)

    !> elemental slips
    real(rk), allocatable :: slip(:, :)

    !> elemental rotation associated to slip, in crystal frame, instantaneous
    real(rk), allocatable :: wp_ss(:, :, :)

    real(rk), allocatable :: e_bar_vec(:, :, :)

    real(rk), allocatable :: detv(:, :)

  end type results_type

  end module types_mod
