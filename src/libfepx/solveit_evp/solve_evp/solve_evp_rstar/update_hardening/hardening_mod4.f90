! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module hardening_mod4

! Module containing definition of hardening evolution equations

! Contains subroutines:
! hard_law: Evaluate hardening rate or derivative of hardening rate
! calculate_shrate: Calculate the shear rate

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use crys_type_mod2
  use matrix_operations_mod
  use parallel_mod, only: par_quit
  use gather_scatter_mod

  implicit none

  public

contains

  subroutine calculate_shrate(mesh, sliprate, shrate, structure, ind, my_crys)

    ! SHRATEcalculations are based upon the assumption that slip systems that
    !   share the same slip plane interact instead of the same slip direction.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! shrate: Shear rate
    ! structure: Crystal type
    ! ind:
    ! phase: Phase id
    type(mesh_type), intent(in) :: mesh
    real(rk), intent(in) :: sliprate(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(out) :: shrate(mesh%maxnumslip)
    character(len=*), intent(in) :: structure
    integer, intent(in) :: ind
    type(crys_type), intent(inout) :: my_crys

    ! Locals:

    real(rk) :: ones(mesh%maxnumslip, mesh%maxnumslip)

    !---------------------------------------------------------------------------

    ones = 1.0d0

    select case (structure) ! Values are hard wired

    case ("fcc")
      ! Fully anisotropic

      if (my_crys%interaction_matrix_parameters_num .eq. 2) then
        call matrix_vec_mult(my_crys%fcc_h1, &
            & abs(sliprate(1:12, ind)), shrate(1:12), 12)

        ! Co-planar

      else if (my_crys%interaction_matrix_parameters_num .eq. 5) then
        call matrix_vec_mult(my_crys%fcc_h1, &
            & abs(sliprate(1:3, ind)), shrate(1:3), 3)
        call matrix_vec_mult(my_crys%fcc_h2, &
            & abs(sliprate(4:6, ind)), shrate(4:6), 3)
        call matrix_vec_mult(my_crys%fcc_h3, &
            & abs(sliprate(7:9, ind)), shrate(7:9), 3)
        call matrix_vec_mult(my_crys%fcc_h4, &
            & abs(sliprate(10:12, ind)), shrate(10:12), 3)
      end if

    case ("bcc")
      ! Fully anisotropic

      if (my_crys%interaction_matrix_parameters_num .eq. 2) then
        call matrix_vec_mult(my_crys%bcc_h1, &
            & abs(sliprate(1:12, ind)), shrate(1:12), 12)

        ! Co-planar

      else if (my_crys%interaction_matrix_parameters_num .eq. 7) then
        call matrix_vec_mult(my_crys%bcc_h1, &
            & abs(sliprate(1:2, ind)), shrate(1:2), 2)
        call matrix_vec_mult(my_crys%bcc_h2, &
            & abs(sliprate(3:4, ind)), shrate(3:4), 2)
        call matrix_vec_mult(my_crys%bcc_h3, &
            & abs(sliprate(5:6, ind)), shrate(5:6), 2)
        call matrix_vec_mult(my_crys%bcc_h4, &
            & abs(sliprate(7:8, ind)), shrate(7:8), 2)
        call matrix_vec_mult(my_crys%bcc_h5, &
            & abs(sliprate(9:10, ind)), shrate(9:10), 2)
        call matrix_vec_mult(my_crys%bcc_h6, &
            & abs(sliprate(11:12, ind)), shrate(11:12), 2)
      end if

      if (my_crys%g_0_bcc_112 .gt. 0.0d0) then
        ! Fully anisotropic

        if (my_crys%interaction_matrix_parameters_112_num .eq. 2) then
          call matrix_vec_mult(my_crys%bcc_112_h1, &
              & abs(sliprate(13:24, ind)), shrate(13:24), 12)

          ! Co-planar

        else if (my_crys%interaction_matrix_parameters_112_num .eq. 7) then
          call matrix_vec_mult(my_crys%bcc_112_h1, &
              & abs(sliprate(13:14, ind)), shrate(13:14), 2)
          call matrix_vec_mult(my_crys%bcc_112_h2, &
              & abs(sliprate(15:16, ind)), shrate(15:16), 2)
          call matrix_vec_mult(my_crys%bcc_112_h3, &
              & abs(sliprate(17:18, ind)), shrate(17:18), 2)
          call matrix_vec_mult(my_crys%bcc_112_h4, &
              & abs(sliprate(19:20, ind)), shrate(19:20), 2)
          call matrix_vec_mult(my_crys%bcc_112_h5, &
              & abs(sliprate(21:22, ind)), shrate(21:22), 2)
          call matrix_vec_mult(my_crys%bcc_112_h6, &
              & abs(sliprate(23:24, ind)), shrate(23:24), 2)
        end if
      end if

    case ("hcp")
      ! Fully anisotropic

      if (my_crys%interaction_matrix_parameters_num .eq. 2) then
        call matrix_vec_mult(my_crys%hcp_h1, &
            & abs(sliprate(1:18, ind)), shrate(1:18), 18)

        ! Co-planar

      else if (my_crys%interaction_matrix_parameters_num .eq. 8) then
        call matrix_vec_mult(my_crys%hcp_h1, &
            & abs(sliprate(1:3, ind)), shrate(1:3), 3)

        shrate(4:6) = my_crys%hcp_vert*abs(sliprate(4:6, ind))

        call matrix_vec_mult(my_crys%hcp_h2, &
            & abs(sliprate(7:8, ind)), shrate(7:8), 2)
        call matrix_vec_mult(my_crys%hcp_h3, &
            & abs(sliprate(9:10, ind)), shrate(9:10), 2)
        call matrix_vec_mult(my_crys%hcp_h4, &
            & abs(sliprate(11:12, ind)), shrate(11:12), 2)
        call matrix_vec_mult(my_crys%hcp_h5, &
            & abs(sliprate(13:14, ind)), shrate(13:14), 2)
        call matrix_vec_mult(my_crys%hcp_h6, &
            & abs(sliprate(15:16, ind)), shrate(15:16), 2)
        call matrix_vec_mult(my_crys%hcp_h7, &
            & abs(sliprate(17:18, ind)), shrate(17:18), 2)
      end if

    case ("bct")
      ! Fully anisotropic

      if (my_crys%interaction_matrix_parameters_num .eq. 2) then
        call matrix_vec_mult(my_crys%bct_h1, &
            & abs(sliprate(1:32, ind)), shrate(1:32), 32)

        ! Co-planar

      else if (my_crys%interaction_matrix_parameters_num .eq. 11) then
        call matrix_vec_mult(my_crys%bct_h1, &
            & abs(sliprate(1:2, ind)), shrate(1:2), 2)
        call matrix_vec_mult(my_crys%bct_h2, &
            & abs(sliprate(3:4, ind)), shrate(3:4), 2)
        call matrix_vec_mult(my_crys%bct_h3, &
            & abs(sliprate(5:6, ind)), shrate(5:6), 2)
        call matrix_vec_mult(my_crys%bct_h4, &
            & abs(sliprate(7:10, ind)), shrate(7:10), 4)
        call matrix_vec_mult(my_crys%bct_h5, &
            & abs(sliprate(11:12, ind)), shrate(11:12), 2)
        call matrix_vec_mult(my_crys%bct_h6, &
            & abs(sliprate(13:16, ind)), shrate(13:16), 4)
        call matrix_vec_mult(my_crys%bct_h7, &
            & abs(sliprate(17:18, ind)), shrate(17:18), 2)
        call matrix_vec_mult(my_crys%bct_h8, &
            & abs(sliprate(19:20, ind)), shrate(19:20), 2)
        call matrix_vec_mult(my_crys%bct_h9, &
            & abs(sliprate(21:24, ind)), shrate(21:24), 4)
        call matrix_vec_mult(my_crys%bct_h10, &
            & abs(sliprate(25:32, ind)), shrate(25:32), 8)
      end if

    case default
      call matrix_vec_mult(ones, abs(sliprate(:, ind)), &
          & shrate, mesh%maxnumslip)

    end select

  end subroutine calculate_shrate

end module hardening_mod4
