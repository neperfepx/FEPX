! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module driver_uniaxial_control_mod2

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use kinematics_mod
  use matrix_operations_mod
  use solve_evp_mod
  use quadrature_mod
  use shape_3d_mod
  use surface_mod
  use units_mod
  use gather_scatter_mod
  use parallel_mod

  implicit none

  public

contains

  subroutine driver_uniaxial_control_printtoterminal (loading, printing, step, &
             & time, dtime, restart_incr, incr)

    type(loading_type), intent(in) :: loading
    type(printing_type), intent(in) :: printing
    real(rk), intent(in) :: time, dtime
    integer, intent(in) :: restart_incr, incr, step

    character(len=16) :: field, time_string, dtime_string

    if (myid .eq. 0) then
      if (loading%step_complete .or. (incr .eq. 1)) then
        write (*, '(a,i0,a)') 'Info   : Running step ', step, '...'
      end if

      if (loading%step_complete .or. (incr .eq. 1)) then
        write (field, '(f0.4)') time
        if (field(1:1) .eq. '.') field = '0'//field(1:15)
        time_string = field
        write (field, '(f0.4)') dtime
        if (field(1:1) .eq. '.') field = '0'//field(1:15)
        dtime_string = field

        write (*, '(a,i0,a,a,a,a,a)') 'Info   :   - &
            &Increment ', incr, ': t = ', trim(time_string), '&
            & secs, dt = ', trim(dtime_string), ' secs'

      else if ((printing%restart .eqv. .true.) .and. &
              & restart_incr .eq. (incr - 1)) then
        write (field, '(f0.4)') time
        if (field(1:1) .eq. '.') field = '0'//field(1:15)
        time_string = field
        write (field, '(f0.4)') dtime
        if (field(1:1) .eq. '.') field = '0'//field(1:15)
        dtime_string = field

        write (*, '(a,i0,a,a,a,a,a)') 'Info   :   - &
            &Increment ', incr, ': t = ', trim(time_string), '&
            & secs, dt = ', trim(dtime_string), ' secs'

      else
        write (field, '(f0.4)') time
        if (field(1:1) .eq. '.') field = '0'//field(1:15)
        time_string = field
        write (field, '(f0.4)') dtime
        if (field(1:1) .eq. '.') field = '0'//field(1:15)
        dtime_string = field

        write (*, '(/,a,i0,a,a,a,a,a)') 'Info   :   - &
            &Increment ', incr, ': t = ', trim(time_string), '&
            & secs, dt = ', trim(dtime_string), ' secs'
      end if
    end if

  end subroutine driver_uniaxial_control_printtoterminal

end module driver_uniaxial_control_mod2
