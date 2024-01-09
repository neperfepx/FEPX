! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module loading_options_type_mod

! Module to define types

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use types_mod, only: loading_options_type

  implicit none

  public

  contains

  subroutine loading_options_set_default(loading_options)

    type(loading_options_type), intent(out) :: loading_options

      loading_options%def_control_by = ""
      loading_options%target = ""

      loading_options%num_bcs = 0

      loading_options%bc_var = ""
      loading_options%bc_nset = ""
      loading_options%bc_dir = ""
      loading_options%bc_vel = 0.0d0

      loading_options%num_mpcs = 0

      loading_options%num_steps = 0
      loading_options%number_of_strain_rate_jumps = 0
      loading_options%number_of_load_rate_jumps = 0
      loading_options%number_of_dwell_episodes = 0
      loading_options%dwell_max_strain_incr = 0.001d0
      loading_options%boundary_conditions = ""
      loading_options%loading_direction = 0
      loading_options%load_rate = 0.0d0
      loading_options%strain_rate = 0.0d0

  end subroutine loading_options_set_default

  end module loading_options_type_mod
