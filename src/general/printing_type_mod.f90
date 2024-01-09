! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module printing_type_mod

! Module to define types

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use types_mod, only: printing_type

  implicit none

  public

  contains

  subroutine printing_set_default (printing)

    type(printing_type), intent(out) :: printing

    ! Initialize printing options
    printing%print_conv = .false.
    printing%print_coo = .false.
    printing%print_crss = .false.
    printing%print_defrate = .false.
    printing%print_defrate_eq = .false.
    printing%print_defrate_pl = .false.
    printing%print_defrate_pl_eq = .false.
    printing%print_disp = .false.
    printing%print_forces = .false.
    printing%print_ori = .false.
    printing%print_restart = .false.
    printing%print_slip = .false.
    printing%print_sliprate = .false.
    printing%print_spinrate = .false.
    printing%print_rotrate = .false.
    printing%print_rotrate_slip = .false.
    printing%print_rotrate_spin = .false.
    printing%print_strain = .false.
    printing%print_strain_el = .false.
    printing%print_strain_el_eq = .false.
    printing%print_strain_eq = .false.
    printing%print_strain_pl = .false.
    printing%print_strain_pl_eq = .false.
    printing%print_stress = .false.
    printing%print_stress_eq = .false.
    printing%print_vel = .false.
    printing%print_velgrad = .false.
    printing%print_work = .false.
    printing%print_work_pl = .false.
    printing%print_workrate = .false.
    printing%print_workrate_pl = .false.
    printing%print_rss = .false.
    printing%restart = .false.
    printing%restart_file_handling = "normal"
    printing%restart_initial_step = -1 ! could probably be removed

  end subroutine printing_set_default

  end module printing_type_mod
