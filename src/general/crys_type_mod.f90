! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module crys_type_mod

! Module to handle processing of parsed crystal data.

! Contains subroutines:
! process_material_parameters: Primary routine for processing crystal parameters
! set_tmin: Set minimum rss to prevent overflow
! set_fcc_block_matrices: Define hardening interaction matrix for fcc crystal
! set_bcc_block_matrices: Define hardening interaction matrix for bcc crystal
! set_hcp_block_matrices: Define hardening interaction matrix for hcp crystal
! set_bct_block_matrices: Define hardening interaction matrix for bct crystal

use, intrinsic :: iso_fortran_env, only: rk => real64
use general_mod
use types_mod
use crys_type_mod2
use crys_type_mod3
use matrix_operations_mod
use parallel_mod
use orientation_conversion_mod
use units_mod
use utils_mod

implicit none

public

contains

subroutine crys_initparams(crys)

  ! Routine to process all crystal parameter info from data read in from the
  ! configuration file.

  !---------------------------------------------------------------------------

  ! Arguments:
  ! mesh: Mesh type
  type(crys_type) :: crys (:)
  ! Locals:
  ! iphase: Looping index over number of phases
  ! c11,c12,c13,c33,c44: Crystal elastic constants for phase.
  ! cutting_precip: Temporary variable for cutting based strengthening
  ! bowing_precip: Temporary variable for bowing based strengthening
  ! l_precip: Temporary variable for bowing based strengthening
  ! g_voigt: Temporay variable for Voigt shear modulus calculation
  ! g_reuss: Temporay variable for Reuss shear modulus calculation
  ! c_sqr: Temporary variable for Voigt and Reuss calculations
  integer :: iphase
  real(rk) :: c11, c12, c13, c33, c44, c66
  real(rk) :: cutting_precip, bowing_precip, l_precip
  real(rk) :: g_voigt, g_reuss, c_sqr
  real(rk):: pi = 4*atan(1.0_rk)

  ! Notes:
  ! The elastic constants for each phase (Cij) are written in the Strength of
  ! Materials (som) convention and c44 is multiplied by 2 after being read in
  ! from the *.cfg file.

  ! The scaling terms for the hardening interaction
  ! matrix (diag, h1-h7) follow the implementation of anisotropic
  ! hardening into FEpX by Carson et. al. in:
  ! Characterizing heterogeneous intragranular deformations in poly-
  ! -crystalline solids using diffraction-based and mechanics-based metrics
  ! https://doi.org/10.1088/1361-651x/aa6dc5

  ! Similarly, the cyclic hardening parameters follow the implementation of
  ! Turkmen's model into FEpX by Turkmen et. al in:
  ! A formulation for the pseudo‐saturation behavior observed during
  ! variable amplitude multiaxial cyclic plasticity
  ! https://doi.org/10.1063/1.1766796

  !---------------------------------------------------------------------------

  ! Loop over each phase and push the material data into the proper arrays.

  do iphase = 1, size(crys)

    select case (crys(iphase)%structure)

      case("fcc", "bcc")

        ! Internal assignment for brevity
        c11 = crys(iphase)%c11
        c12 = crys(iphase)%c12
        c13 = crys(iphase)%c13 ! Auto assigned from c12 in read_cfg_mod.f90
        c44 = crys(iphase)%c44

        ! Elastic coefficents and stiffness tensor
        crys(iphase)%elas_coeffs = (/c11, c12, c12, c44, 0.0d0/)
        crys(iphase)%kelas = (/c11 - c12, c11 - c12, 2.0d0*c44, 2.0d0*c44, 2.0d0*c44/)
        crys(iphase)%keinv = 1.0d0 / crys(iphase)%kelas

        ! Bulk modulus
        crys(iphase)%bulk_mod = (c11 + 2.0d0*c12)/3.0d0

        ! Bulk shear modulus
        ! Calculated using an average of Voigt and Reuss values, from:
        ! doi.org/10.1016/j.jmrt.2020.01.086, Eqs. 6 and 7
        g_voigt = (c11 - c12 + (3*c44)) / 5
        g_reuss = (5*(c11 - c12)*c44) / ((3*c11) - (3*c12) + (4*c44))
        crys(iphase)%bulk_shear_mod = (g_voigt + g_reuss)/2

      case ("hcp")

        ! Internal assignment for brevity
        c11 = crys(iphase)%c11
        c12 = crys(iphase)%c12
        c13 = crys(iphase)%c13
        c44 = crys(iphase)%c44
        ! Compute c33 and c66 (from deviatioric/hydrostatic decoupling)
        c33 = c11 + c12 - c13
        c66 = (c11 - c12)/2

        ! Elastic coefficents and stiffness tensor
        crys(iphase)%elas_coeffs = (/c11, c12, c13, c44, 0.0d0/)
        crys(iphase)%kelas = (/c11 - c12, c11 + c12 - 2.0d0*c13, 2.0d0*c44, 2.0d0*c44, c11 - c12/)
        crys(iphase)%keinv = 1.0d0 / crys(iphase)%kelas

        ! Bulk modulus
        crys(iphase)%bulk_mod = (c11 + c12 + c13)/3.0d0

        ! Bulk shear modulus
        ! Calculated using an average of Voigt and Reuss values, from:
        ! doi.org/10.1016/j.jmrt.2020.01.086, Eqs. 9 and 11
        g_voigt = ((7*c11) - (5*c12) + (12*c44) + (2*c33) - (4*c13))/30
        c_sqr = ((c11 + c12)*c33) - (2*(c13**2))
        g_reuss = (5*c_sqr*c44*c66)/ &
          & ((6*crys(iphase)%bulk_mod*c44*c66) + (2*c_sqr*(c44 + c66)))
        crys(iphase)%bulk_shear_mod = (g_voigt + g_reuss)/2

      case ("bct")

        ! Internal assignment for brevity
        c11 = crys(iphase)%c11
        c12 = crys(iphase)%c12
        c13 = crys(iphase)%c13
        c44 = crys(iphase)%c44
        c66 = crys(iphase)%c66
        ! Compute c33 (from deviatioric/hydrostatic decoupling)
        c33 = c11 + c12 - c13

        ! Elastic coefficents and stiffness tensor
        crys(iphase)%elas_coeffs = (/c11, c12, c13, c44, c66/)
        crys(iphase)%kelas = (/c11 - c12, c11 + c12 - 2.0d0*c13, 2.0d0*c44, 2.0d0*c44, 2.0d0*c66/)
        crys(iphase)%keinv = 1.0d0 / crys(iphase)%kelas

        ! Bulk modulus
        crys(iphase)%bulk_mod = (c11 + c12 + c13)/3.0d0

        ! Bulk shear modulus
        ! Calculated using an average of Voigt and Reuss values, from:
        ! doi.org/10.1016/j.jmrt.2020.01.086, Eqs. 9 and 11 (assumed hcp)
        g_voigt = ((7*c11) - (5*c12) + (12*c44) + (2*c33) - (4*c13))/30
        c_sqr = ((c11 + c12)*c33) - (2*(c13**2))
        g_reuss = (5*c_sqr*c44*c66)/ &
          & ((6*crys(iphase)%bulk_mod*c44*c66) + (2*c_sqr*(c44 + c66)))
        crys(iphase)%bulk_shear_mod = (g_voigt + g_reuss)/2

      case default

        call par_quit("Error  : crys_initparams: invalid crystal type")

    end select

    ! Consider precipitation hardening model if available
    if (crys(iphase)%precipitation) then
      l_precip = crys(iphase)%r_p*sqrt(pi/crys(iphase)%f_p)
      bowing_precip = (crys(iphase)%c_p*crys(iphase)%bulk_shear_mod*crys(iphase)%b_p)/(l_precip - (2*crys(iphase)%r_p))
      cutting_precip = crys(iphase)%a_p*(sqrt((crys(iphase)%f_p*crys(iphase)%r_p)/crys(iphase)%b_p))
      ! Lesser of the two values will determine the strengthening. If cutting
      ! is greater than bowing, we must be in the bowing regime and vice versa
      if (bowing_precip .ge. cutting_precip) then
        crys(iphase)%g_0 = crys(iphase)%g_0 + cutting_precip
      else if (bowing_precip .lt. cutting_precip) then
        crys(iphase)%g_0 = crys(iphase)%g_0 + bowing_precip
      end if
    else if (crys(iphase)%precipitation_cutting) then
      crys(iphase)%g_0 = crys(iphase)%g_0 + &
        & crys(iphase)%a_p*(sqrt((crys(iphase)%f_p*crys(iphase)%r_p)/crys(iphase)%b_p))
    else if (crys(iphase)%precipitation_bowing) then
      l_precip = crys(iphase)%r_p*sqrt(pi/crys(iphase)%f_p)
      crys(iphase)%g_0 = crys(iphase)%g_0 + &
        & (crys(iphase)%c_p*crys(iphase)%bulk_shear_mod*crys(iphase)%b_p)/(l_precip - (2*crys(iphase)%r_p))
    end if

    ! Make the CrystalType object (crys).
    select case (crys(iphase)%structure)

      case ("fcc")
        call crys_create(crys(iphase), crys(iphase)%structure, &
          & hratio_cubic = crys(iphase)%hratio_cubic(1:12), &
          & self = crys(iphase))

      case ("bcc")
        if (crys(iphase)%g_0_bcc_112 .le. 0.0d0) then
          call crys_create(crys(iphase), crys(iphase)%structure, &
            & hratio_cubic = crys(iphase)%hratio_cubic(1:12), &
            & self = crys(iphase))
        else if (crys(iphase)%g_0_bcc_112 .gt. 0.0d0) then
          call crys_create(crys(iphase), crys(iphase)%structure, &
            & hratio_cubic = crys(iphase)%hratio_cubic(1:12), &
            & hratio_cubic_112 = crys(iphase)%hratio_cubic_112(1:12), &
            & self = crys(iphase))
        end if

      case ("hcp")
        call crys_create(crys(iphase), crys(iphase)%structure, &
          & c_over_a = crys(iphase)%c_over_a, &
          & hratio_hcp = crys(iphase)%hratio_hcp(1:18), &
          & self = crys(iphase))

      case ("bct")
        call crys_create(crys(iphase), crys(iphase)%structure, &
          & c_over_a = crys(iphase)%c_over_a, &
          & hratio_bct = crys(iphase)%hratio_bct(1:32), &
          & self = crys(iphase))

    end select

    ! Set minimum value for tau to prevent overflow.
    call set_tmin(crys(iphase)%t_min, crys(iphase)%m)

    ! Construct anisotropic hardening interaction matrices
    if (ut_list_testelt (crys(iphase)%hardening, ',', "anisotropic")) then
      select case (crys(iphase)%structure)

        case ("fcc")
          call set_fcc_block_matrices(&
            & crys(iphase)%interaction_matrix_parameters(1), &
            & crys(iphase)%interaction_matrix_parameters(2), &
            & crys(iphase)%interaction_matrix_parameters(3), &
            & crys(iphase)%interaction_matrix_parameters(4), &
            & crys(iphase)%interaction_matrix_parameters(5), crys(iphase))

        case ("bcc")
          call set_bcc_block_matrices(&
            & crys(iphase)%interaction_matrix_parameters(1), &
            & crys(iphase)%interaction_matrix_parameters(2), &
            & crys(iphase)%interaction_matrix_parameters(3), &
            & crys(iphase)%interaction_matrix_parameters(4), &
            & crys(iphase)%interaction_matrix_parameters(5), &
            & crys(iphase)%interaction_matrix_parameters(6), &
            & crys(iphase)%interaction_matrix_parameters(7), crys(iphase))

        case ("hcp")
          call set_hcp_block_matrices(&
            & crys(iphase)%interaction_matrix_parameters(1), &
            & crys(iphase)%interaction_matrix_parameters(2), &
            & crys(iphase)%interaction_matrix_parameters(3), &
            & crys(iphase)%interaction_matrix_parameters(4), &
            & crys(iphase)%interaction_matrix_parameters(5), &
            & crys(iphase)%interaction_matrix_parameters(6), &
            & crys(iphase)%interaction_matrix_parameters(7), &
            & crys(iphase)%interaction_matrix_parameters(8), crys(iphase))

        case ("bct")
          call set_bct_block_matrices(&
            & crys(iphase)%interaction_matrix_parameters(1), &
            & crys(iphase)%interaction_matrix_parameters(2), &
            & crys(iphase)%interaction_matrix_parameters(3), &
            & crys(iphase)%interaction_matrix_parameters(4), &
            & crys(iphase)%interaction_matrix_parameters(5), &
            & crys(iphase)%interaction_matrix_parameters(6), &
            & crys(iphase)%interaction_matrix_parameters(7), &
            & crys(iphase)%interaction_matrix_parameters(8), &
            & crys(iphase)%interaction_matrix_parameters(9), &
            & crys(iphase)%interaction_matrix_parameters(10), &
            & crys(iphase)%interaction_matrix_parameters(11), crys(iphase))

        case default
          call par_quit("Error  : crys_initparams: invalid crystal type")

      end select
    end if
  end do

end subroutine crys_initparams

  subroutine crys_get(self, dev, skw, pptrans, vertices, numslip, &
      & numvertices)

    ! Return deviatoric parts of Schmid tensors

    !---------------------------------------------------------------------------

    ! Arguments:

    ! self: The CrystalType object
    ! dev: Deviatoric part of Schmid tensors
    ! skw: Skew part of Schmid tensors
    ! pptrans: Matrices of Schmid diads
    ! vertices: Vertices in 5-vector form
    ! numslip: Number of slip systems
    ! numvertices: Number of vertices

    type(crys_type) :: self
    real(rk), pointer, optional :: dev(:, :)
    real(rk), pointer, optional :: skw(:, :)
    real(rk), pointer, optional :: pptrans(:, :, :)
    real(rk), pointer, optional :: vertices(:, :)
    integer, optional :: numslip
    integer, optional :: numvertices

    ! Locals:

    integer :: dshape(2)
    integer :: sshape(2)
    integer :: argpptshape(3)
    integer :: mypptshape(3)
    integer :: i

    real(rk), pointer :: mydev(:, :)

    !---------------------------------------------------------------------------

    if (present(dev)) then
      sshape = (/5, self%numslip/)
      !shape(self%schmid_sym)

      if (associated(dev)) then
        ! Could check dimensions here ...

        dshape = shape(dev)

        if ((dshape(1) .ne. sshape(1)) .or. (dshape(2) .ne. sshape(2))) &
            & then
          deallocate (dev)
          allocate (dev(sshape(1), sshape(2)))
        end if

      else
        allocate (dev(sshape(1), sshape(2)))
      end if

      call tensor3ddecompose(self%schmid_3x3, dev=dev)
    end if

    if (present(skw)) then
      sshape = (/3, self%numslip/)

      if (associated(skw)) then
        ! Could check dimensions here ...

        dshape = shape(skw)

        if ((dshape(1) .ne. sshape(1)) .or. (dshape(2) .ne. sshape(2))) &
            & then
          deallocate (skw)
          allocate (skw(sshape(1), sshape(2)))
        end if

      else
        allocate (skw(sshape(1), sshape(2)))
      end if

      call tensor3ddecompose(self%schmid_3x3, skw=skw)
    end if

    if (present(pptrans)) then
      mypptshape = (/5, 5, self%numslip/)

      if (associated(pptrans)) then
        !  Could check dimensions here ...

        argpptshape = shape(pptrans)

        if ((argpptshape(1) .ne. mypptshape(1)) .or. &
            & (argpptshape(2) .ne. mypptshape(2)) .or. &
            & (argpptshape(3) .ne. mypptshape(3))) then
          deallocate (pptrans)
          allocate (pptrans(mypptshape(1), mypptshape(2), mypptshape(3)))
        end if

      else
        allocate (pptrans(mypptshape(1), mypptshape(2), mypptshape(3)))
      end if

      allocate (mydev(5, self%numslip))

      call tensor3ddecompose(self%schmid_3x3, dev=mydev)

      do i = 1, self%numslip
        pptrans(:, :, i) = matmul(&
            & reshape(mydev(:, i), shape=(/5, 1/)), &
            & reshape(mydev(:, i), shape=(/1, 5/)))
      end do

      deallocate (mydev)
    end if

    if (present(vertices)) then
      sshape = (/5, self%numvertices/)

      if (associated(vertices)) then
        !  Could check dimensions here ...

        dshape = shape(vertices)

        if ((dshape(1) .ne. sshape(1)) .or. (dshape(2) .ne. sshape(2))) then
          deallocate (vertices)
          allocate (vertices(sshape(1), sshape(2)))
        end if

      else
        allocate (vertices(sshape(1), sshape(2)))
      end if

      vertices = self%vertices
    end if

    if (present(numslip)) then
      numslip = self%numslip
    end if

    if (present(numvertices)) then
      numvertices = self%numvertices
    end if

    ! Need call to get vertices as 3x3.

  end subroutine crys_get

  subroutine crys_maxnumvert (crys, maxnumvert)

    type(crys_type), intent(in) :: crys (:)
    integer, intent(out) :: maxnumvert

    integer :: i

    maxnumvert = 0
    do i = 1, size(crys)
      maxnumvert = max(maxnumvert, crys(i)%numvertices)
    end do

  end subroutine crys_maxnumvert

  subroutine crys_set_default (crys)

    ! Allocates and initializes all variables that are dependent on the
    ! number of phases present in the configuration file.
    ! Arguments:

    type(crys_type), intent(out) :: crys

    !-----------------------------------------------------------------------------

    allocate(crys%interaction_matrix_parameters(11))
    allocate(crys%hratio_cubic(12)) ! FCC and BCC
    allocate(crys%hratio_hcp(18)) ! HCP
    allocate(crys%hratio_bct(32)) ! BCT
    allocate(crys%hratio_cubic_112(12)) ! FCC and BCC

    ! Initialize all variables as -1 or -1.0d0
    crys%structure = ""
    crys%m = -1.0d0
    crys%gammadot_0 = -1.0d0
    crys%h_0 = -1.0d0
    crys%g_0 = -1.0d0
    crys%g_0_bcc_112 = -1.0d0
    crys%g_s0 = -1.0d0
    crys%n = -1.0d0
    crys%c11 = -1.0d0
    crys%c12 = -1.0d0
    crys%c13 = -1.0d0
    crys%c44 = -1.0d0
    crys%c66 = -1.0d0
    crys%c_over_a = -1.0d0
    crys%m_prime = -1.0d0
    crys%gammadot_s0 = -1.0d0

    crys%hardening = 'saturation,isotropic'
    crys%anisotropic = .false.
    crys%cyclic = .false.
    crys%saturation_evolution = .false.
    crys%precipitation = .false.
    crys%precipitation_cutting = .false.
    crys%precipitation_bowing = .false.

    crys%cyclic_a = -1.0d0
    crys%cyclic_c = -1.0d0
    crys%interaction_matrix_parameters = -1.0d0
    crys%interaction_matrix_parameters_num = -1.0d0
    crys%a_p = -1.0d0
    crys%f_p = -1.0d0
    crys%b_p = -1.0d0
    crys%r_p = -1.0d0
    crys%b_m = -1.0d0
    crys%c_p = -1.0d0
    crys%use_aniso_m = .false.
    crys%aniso_m = -1.0d0
    crys%bulk_mod = -1.0d0
    crys%bulk_shear_mod = -1.0d0
    ! Assume isotropic g_0 unless input specifies otherwise
    crys%hratio_num = 1.0d0
    crys%hratio_num_112 = 1.0d0
    crys%hratio_cubic = 1.0d0
    crys%hratio_cubic_112 = 1.0d0
    crys%hratio_hcp = 1.0d0
    crys%hratio_bct = 1.0d0
    crys%numslip = 0.0d0
    crys%numvertices = 0.0d0

  end subroutine crys_set_default

end module crys_type_mod
