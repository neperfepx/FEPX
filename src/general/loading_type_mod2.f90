! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module loading_type_mod2

! Module to define types

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use parallel_mod
  use types_mod, only: loading_type, exec_type

  implicit none

  public

  contains

  !> Check if load is in specified range.
  function loading_load_isinrange(loading, exec, load) result(status)

    implicit none

    type(loading_type), intent(in) :: loading
    type(exec_type), intent(in) :: exec
    real(rk), dimension(3), intent(in) :: load

    logical :: status

    if ((load(loading%loading_direction) - loading%target_load(loading%curr_step, 1)) &
         & * loading%target_sign(loading%curr_step) .gt. -exec%load_tol_rel) then
      status = .true.
    else
      status = .false.
    end if

    return

    !---------------------------------------------------------------------------

  end function loading_load_isinrange

end module loading_type_mod2
