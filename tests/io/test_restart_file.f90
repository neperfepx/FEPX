! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of restart_file_mod (step 8.4, dev/theory/addendum_restart.md §5): the round trip
!> of every record kind (real ranks 0-4, integer ranks 0-1, logical ranks 0-1, character), of a
!> whole material state (every component filled with distinct values) and of a schedule step
!> (step_record_t), bit for bit; the header of a file read back. Error paths: test
!> restart_file_errors (tests/CMakeLists.txt).
!>
!> Source : FEPX (ROADMAP step 8.4).
!>
!> Locals:
!>   name      character(*)            scratch file name (working directory)
!>   file      type(restart_file_t)    the file
!>   r0, r1, r2, r3, r4  real(rk)      values written (ranks 0-4)
!>   s0, s1, s2, s3, s4  real(rk)      values read back
!>   i0, i1, j0, j1      integer       integers written, read back
!>   l0, l1, m0, m1      logical       logicals written, read back
!>   text      character(:), allocatable  a string read back
!>   state, back  type(material_state_t)  a state written, read back
!>   rec, rec_back  type(step_record_t)   a schedule step written, read back
!>   k         integer                 index
program test_restart_file

  use kinds_mod, only: rk, ik
  use restart_file_mod, only: restart_file_t, restart_create, restart_open, restart_close, &
    & restart_put, restart_get, restart_put_state, restart_get_state, restart_put_step, &
    & restart_get_step, restart_rename
  use material_state_mod, only: material_state_t, material_state_alloc
  use step_schedule_mod, only: step_record_t
  use assert_mod, only: assert_true, test_summary

  implicit none

  character(len=*), parameter :: name = 'test_restart_file.rst'
  type(restart_file_t) :: file
  real(rk) :: r0, r1(3), r2(2, 3), r3(2, 2, 2), r4(2, 1, 2, 2)
  real(rk) :: s0, s1(3), s2(2, 3), s3(2, 2, 2), s4(2, 1, 2, 2)
  integer :: i0, i1(4), j0, j1(4)
  logical :: l0, l1(3), m0, m1(3)
  character(len=:), allocatable :: text
  type(material_state_t) :: state, back
  type(step_record_t) :: rec, rec_back
  integer :: k

  ! Values with every bit used (thirds, negative, tiny, huge).
  r0 = 1.0_rk/3.0_rk
  r1 = [-2.0_rk/7.0_rk, tiny(1.0_rk), huge(1.0_rk)]
  r2 = reshape([(real(k, rk)/11.0_rk, k = 1, 6)], [2, 3])
  r3 = reshape([(-real(k, rk)/13.0_rk, k = 1, 8)], [2, 2, 2])
  r4 = reshape([(real(k, rk)**1.5_rk, k = 1, 8)], [2, 1, 2, 2])
  i0 = -42
  i1 = [1, -2, huge(1), 0]
  l0 = .true.
  l1 = [.false., .true., .false.]

  ! A state of 2 points, 3 elements, 4 systems, every component different.
  call material_state_alloc(state, 2_ik, 3_ik, 4_ik)
  call fill_state(state)
  call material_state_alloc(back, 2_ik, 3_ik, 4_ik)

  ! A schedule step with every component set.
  rec%target_kind = 3
  rec%is_load = .true.
  rec%target_value = -12.5_rk
  rec%target_tol = 0.1_rk
  rec%load_sign = -1.0_rk
  rec%dt_min = 1.0e-3_rk
  rec%strain_rate = [1.0e-3_rk, -2.0e-3_rk, 0.0_rk]
  rec%vel_factor = -1.0_rk
  rec%time_start = 1.0_rk/3.0_rk
  rec%time_end = 2.0_rk/3.0_rk
  rec%strain_start = [0.1_rk, 0.2_rk, 0.3_rk]
  rec%strain_end = [0.4_rk, 0.5_rk, 0.6_rk]
  rec%num_incrs = 7
  rec%dt = 0.25_rk
  rec%dt_last = 0.125_rk
  rec%print_results = .true.
  rec%bc_kind = [1, 0, 2]
  rec%bc_target = [1.0_rk, 2.0_rk, 3.0_rk]
  rec%bc_tol = [0.5_rk, 0.6_rk, 0.7_rk]
  rec%bc_mode = [3, 0, 1]
  rec%bc_rate = [-1.0_rk, 0.0_rk, 1.0_rk]
  rec%rate_kind = 2
  rec%is_rate = .true.
  rec%is_dwell = .true.
  rec%is_open = .true.
  rec%ctrl_rate = 100.0_rk
  rec%per_rate = [(real(k, rk)*1.0e-4_rk, k = 1, 9)]

  ! Write, then rename (the file's final name is set by restart_rename).
  call restart_create(file, name // '.new', 2, 1)
  call restart_put(file, 'r0', r0)
  call restart_put(file, 'r1', r1)
  call restart_put(file, 'r2', r2)
  call restart_put(file, 'r3', r3)
  call restart_put(file, 'r4', r4)
  call restart_put(file, 'i0', i0)
  call restart_put(file, 'i1', i1)
  call restart_put(file, 'l0', l0)
  call restart_put(file, 'l1', l1)
  call restart_put(file, 'text', 'a line: 1 2 3')
  call restart_put(file, 'empty', '')
  call restart_put_state(file, state)
  call restart_put_step(file, rec)
  call restart_close(file)
  call restart_rename(name // '.new', name)

  ! Read back in the same order (header checked: 2 ranks, rank 1, even).
  call restart_open(file, name, 2, 1)
  call restart_get(file, 'r0', s0)
  call restart_get(file, 'r1', s1)
  call restart_get(file, 'r2', s2)
  call restart_get(file, 'r3', s3)
  call restart_get(file, 'r4', s4)
  call restart_get(file, 'i0', j0)
  call restart_get(file, 'i1', j1)
  call restart_get(file, 'l0', m0)
  call restart_get(file, 'l1', m1)
  call assert_true(same(s0, r0) .and. all(same1(s1, r1)) .and. &
    & all(same1(reshape(s2, [6]), reshape(r2, [6]))) .and. &
    & all(same1(reshape(s3, [8]), reshape(r3, [8]))) .and. &
    & all(same1(reshape(s4, [8]), reshape(r4, [8]))), 'reals bit for bit')
  call assert_true(j0 == i0 .and. all(j1 == i1), 'integers')
  call assert_true((m0 .eqv. l0) .and. all(m1 .eqv. l1), 'logicals')
  call restart_get(file, 'text', text)
  call assert_true(text == 'a line: 1 2 3', 'character')
  call restart_get(file, 'empty', text)
  call assert_true(len(text) == 0, 'empty character')
  call restart_get_state(file, back)
  call assert_true(states_equal(state, back), 'material state, every component')
  call restart_get_step(file, rec_back)
  call assert_true(records_equal(rec, rec_back), 'schedule step, every component')
  call restart_close(file)

  call test_summary()

contains

  !> Whether two reals have the same bits (equal, as an ordered comparison; NaN not used).
  !> Arguments: a, b (in) real(rk); t (out, result) logical.
  elemental function same(a, b) result(t)
    real(rk), intent(in) :: a, b
    logical :: t
    t = .not. (a < b .or. a > b)
  end function same

  !> same for rank-1 arrays. Arguments: a, b (in) real(rk)(:); t (out, result) logical(:).
  pure function same1(a, b) result(t)
    real(rk), intent(in) :: a(:), b(:)
    logical :: t(size(a))
    t = same(a, b)
  end function same1

  !> Fill every component of a state with distinct values (component c: c + i/1000 for its
  !> i-th value). Arguments: st (inout) type(material_state_t) (allocated).
  subroutine fill_state(st)
    type(material_state_t), intent(inout) :: st
    st%start_from_vertex = .true.
    call fill3(st%stress_s, 1)
    call fill3(st%strain_el_s, 2)
    st%rot_c2s = reshape([(3.0_rk + real(k, rk)/1000.0_rk, k = 1, size(st%rot_c2s))], &
      & shape(st%rot_c2s))
    call fill3(st%spin_lat_s, 4)
    call fill3(st%crss, 5)
    call fill3(st%sliprate, 6)
    call fill3(st%rss, 7)
    call fill3(st%slip, 8)
    call fill3(st%slip_accum, 9)
    call fill3(st%strain_pl_c, 10)
    call fill3(st%defrate_pl_s, 11)
    call fill3(st%spin_pl_s, 12)
    call fill3(st%rotrate_s, 13)
    call fill2(st%workrate, 14)
    call fill2(st%work, 15)
    call fill3(st%defrate_s, 16)
    call fill3(st%spin_s, 17)
    call fill3(st%strain_s, 18)
    call fill3(st%strain_pl_s, 19)
    call fill3(st%spin_hat_update_s, 20)
    call fill3(st%spin_slip_update_s, 21)
    call fill2(st%workrate_pl, 22)
    call fill2(st%work_pl, 23)
    call fill2(st%km_rho, 24)
    call fill2(st%km_k2, 25)
    call fill2(st%km_rhorate, 26)
    call fill3(st%dmg_c, 27)
    call fill3(st%dmgrate_c, 28)
  end subroutine fill_state

  !> Fill a rank-3 array: c + i/1000. Arguments: a (inout) real(rk)(:, :, :); c (in) integer.
  subroutine fill3(a, c)
    real(rk), intent(inout) :: a(:, :, :)
    integer, intent(in) :: c
    a = reshape([(real(c, rk) + real(k, rk)/1000.0_rk, k = 1, size(a))], shape(a))
  end subroutine fill3

  !> Fill a rank-2 array: c + i/1000. Arguments: a (inout) real(rk)(:, :); c (in) integer.
  subroutine fill2(a, c)
    real(rk), intent(inout) :: a(:, :)
    integer, intent(in) :: c
    a = reshape([(real(c, rk) + real(k, rk)/1000.0_rk, k = 1, size(a))], shape(a))
  end subroutine fill2

  !> Whether every component of two states is the same. Arguments: a, b (in)
  !> type(material_state_t); t (out, result) logical.
  function states_equal(a, b) result(t)
    type(material_state_t), intent(in) :: a, b
    logical :: t
    t = (a%start_from_vertex .eqv. b%start_from_vertex) .and. &
      & all(same1(pack(a%stress_s, .true.), pack(b%stress_s, .true.))) .and. &
      & all(same1(pack(a%strain_el_s, .true.), pack(b%strain_el_s, .true.))) .and. &
      & all(same1(pack(a%rot_c2s, .true.), pack(b%rot_c2s, .true.))) .and. &
      & all(same1(pack(a%spin_lat_s, .true.), pack(b%spin_lat_s, .true.))) .and. &
      & all(same1(pack(a%crss, .true.), pack(b%crss, .true.))) .and. &
      & all(same1(pack(a%sliprate, .true.), pack(b%sliprate, .true.))) .and. &
      & all(same1(pack(a%rss, .true.), pack(b%rss, .true.))) .and. &
      & all(same1(pack(a%slip, .true.), pack(b%slip, .true.))) .and. &
      & all(same1(pack(a%slip_accum, .true.), pack(b%slip_accum, .true.))) .and. &
      & all(same1(pack(a%strain_pl_c, .true.), pack(b%strain_pl_c, .true.))) .and. &
      & all(same1(pack(a%defrate_pl_s, .true.), pack(b%defrate_pl_s, .true.))) .and. &
      & all(same1(pack(a%spin_pl_s, .true.), pack(b%spin_pl_s, .true.))) .and. &
      & all(same1(pack(a%rotrate_s, .true.), pack(b%rotrate_s, .true.))) .and. &
      & all(same1(pack(a%workrate, .true.), pack(b%workrate, .true.))) .and. &
      & all(same1(pack(a%work, .true.), pack(b%work, .true.))) .and. &
      & all(same1(pack(a%defrate_s, .true.), pack(b%defrate_s, .true.))) .and. &
      & all(same1(pack(a%spin_s, .true.), pack(b%spin_s, .true.))) .and. &
      & all(same1(pack(a%strain_s, .true.), pack(b%strain_s, .true.))) .and. &
      & all(same1(pack(a%strain_pl_s, .true.), pack(b%strain_pl_s, .true.))) .and. &
      & all(same1(pack(a%spin_hat_update_s, .true.), pack(b%spin_hat_update_s, .true.))) .and. &
      & all(same1(pack(a%spin_slip_update_s, .true.), pack(b%spin_slip_update_s, .true.))) &
      & .and. all(same1(pack(a%workrate_pl, .true.), pack(b%workrate_pl, .true.))) .and. &
      & all(same1(pack(a%work_pl, .true.), pack(b%work_pl, .true.))) .and. &
      & all(same1(pack(a%km_rho, .true.), pack(b%km_rho, .true.))) .and. &
      & all(same1(pack(a%km_k2, .true.), pack(b%km_k2, .true.))) .and. &
      & all(same1(pack(a%km_rhorate, .true.), pack(b%km_rhorate, .true.))) .and. &
      & all(same1(pack(a%dmg_c, .true.), pack(b%dmg_c, .true.))) .and. &
      & all(same1(pack(a%dmgrate_c, .true.), pack(b%dmgrate_c, .true.)))
  end function states_equal

  !> Whether every component of two schedule steps is the same. Arguments: a, b (in)
  !> type(step_record_t); t (out, result) logical.
  function records_equal(a, b) result(t)
    type(step_record_t), intent(in) :: a, b
    logical :: t
    t = a%target_kind == b%target_kind .and. (a%is_load .eqv. b%is_load) .and. &
      & same(a%target_value, b%target_value) .and. same(a%target_tol, b%target_tol) .and. &
      & same(a%load_sign, b%load_sign) .and. same(a%dt_min, b%dt_min) .and. &
      & all(same1(a%strain_rate, b%strain_rate)) .and. same(a%vel_factor, b%vel_factor) .and. &
      & same(a%time_start, b%time_start) .and. same(a%time_end, b%time_end) .and. &
      & all(same1(a%strain_start, b%strain_start)) .and. &
      & all(same1(a%strain_end, b%strain_end)) .and. a%num_incrs == b%num_incrs .and. &
      & same(a%dt, b%dt) .and. same(a%dt_last, b%dt_last) .and. &
      & (a%print_results .eqv. b%print_results) .and. all(a%bc_kind == b%bc_kind) .and. &
      & all(same1(a%bc_target, b%bc_target)) .and. all(same1(a%bc_tol, b%bc_tol)) .and. &
      & all(a%bc_mode == b%bc_mode) .and. all(same1(a%bc_rate, b%bc_rate)) .and. &
      & a%rate_kind == b%rate_kind .and. (a%is_rate .eqv. b%is_rate) .and. &
      & (a%is_dwell .eqv. b%is_dwell) .and. (a%is_open .eqv. b%is_open) .and. &
      & same(a%ctrl_rate, b%ctrl_rate) .and. all(same1(a%per_rate, b%per_rate))
  end function records_equal

end program test_restart_file
