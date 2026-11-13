! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Writes FEPX's Schmid tensors d (x) n (src/material/crystal_mod.f90::crystal_init, crystal
!> frame) for every crystal type to standard output, as the fixture schmid_reference.txt of
!> tests/material/test_crystal.f90.
!> Not part of FEPX; built against the objects of an FEPX build (README.md in this directory).
!>
!> Source : FEPX (ROADMAP step 9.13, D91: replaces v2_dump_schmid.f90, which wrote v2-fixed's
!>          tensors; FEPX's systems are v2's, in the same order).
!>
!> Output: for each case a header line "<structure> <num_slip> <c_over_a>", then one line per slip
!> system with the 9 entries of d (x) n row by row (es25.17). Cases: fcc; bcc with the {112}
!> systems (24); hcp, c/a = 1.587; bct, c/a = 0.5456.
!>
!> Locals:
!>   crys  type(crystal_t)  crystal of the current case
program dump_schmid

  use kinds_mod, only: rk
  use crystal_mod, only: crystal_t, crystal_init, crystal_fcc, crystal_bcc, crystal_hcp, &
    & crystal_bct

  implicit none

  type(crystal_t) :: crys

  ! fcc.
  call crystal_init(crys, crystal_fcc)
  call dump('fcc')
  ! bcc with {112}<111> (systems 13-24).
  call crystal_init(crys, crystal_bcc, families=[1, 2])
  call dump('bcc')
  ! hcp.
  call crystal_init(crys, crystal_hcp, c_over_a=1.587_rk)
  call dump('hcp')
  ! bct.
  call crystal_init(crys, crystal_bct, c_over_a=0.5456_rk)
  call dump('bct')

contains

  !> Writes the header and the Schmid tensors of the current crystal.
  !>
  !> Arguments:
  !>   name  (in)  character       structure name of the case
  !> Locals:
  !>   dn    real(rk)(3, 3)        d (x) n of one system (crystal frame) [-]
  !>   i     integer               slip system index
  subroutine dump(name)
    character(len=*), intent(in) :: name
    real(rk) :: dn(3, 3)
    integer :: i
    write (*, '(a,1x,i0,1x,es25.17)') name, crys%num_slip, crys%c_over_a
    do i = 1, crys%num_slip
      dn = spread(crys%slip_dir_c(:, i), 2, 3)*spread(crys%slip_normal_c(:, i), 1, 3)
      write (*, '(9es25.17)') transpose(dn)
    end do
  end subroutine dump

end program dump_schmid
