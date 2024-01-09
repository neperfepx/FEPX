! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module res_init_mod2

  use general_mod
  use types_mod
  use crys_type_mod2
  use matrix_operations_mod

  implicit none

  public

contains

  subroutine res_init_crss(mesh, crys, results, my_phase)

    ! Initialize the crystal slip system hardness

    !---------------------------------------------------------------------------

    ! Arguments:
    ! my_phase: Crystal phase currently being initialized on this processor

    type(mesh_type) :: mesh
    type(crys_type) :: crys (:)
    type(results_type) :: results
    integer, intent(in)  :: my_phase(elt_sub:elt_sup)

    ! Locals:
    ! indices: (??)
    ! islip: Loop index over number of slip systems
    ! iphase: Loop index over number of crystal phases

    integer, pointer :: indices(:) => null()
    integer :: islip, iphase, num_ind

    !---------------------------------------------------------------------------

    ! Initialize
    results%crss = 0.0d0

    ! Assign the crss based on the crystal type
    do iphase = 1, mesh%num_phases
      ! Finds the indices corresponding to the current phase the loop is on
      call find_indices(my_phase, iphase, indices, num_ind, elt_sub - 1)

      ! Initialize the crystal slip strengths corresponding to the slip number
      do islip = 1, crys(iphase)%numslip
        results%crss(islip, indices, :) = crys(iphase)%g_0
      end do
    end do

  end subroutine res_init_crss

end module res_init_mod2
