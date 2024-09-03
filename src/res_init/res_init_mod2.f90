! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module res_init_mod2

  use general_mod
  use types_mod
  use crys_type_mod2
  use matrix_operations_mod
  use parallel_mod

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
    !
    if (mesh%crss_defined .eqv. .false.) then
      ! Assign the crss based on the crystal type
      do iphase = 1, mesh%num_phases
        ! Finds the indices corresponding to the current phase the loop is on
        call find_indices(my_phase, iphase, indices, num_ind, elt_sub - 1)
        allocate(crys(iphase)%g_0(crys(iphase)%numslip,elt_sub:elt_sup))
        ! if (crys(iphase)%precipitation .eqv. .false.) 
        ! Initialize the crystal slip strengths corresponding to the slip number
        do islip = 1, crys(iphase)%numslip
          crys(iphase)%g_0(islip,:) = crys(iphase)%g_0_tmp(islip) + &
            & crys(iphase)%precipitation_strength
          results%crss(islip, indices, : ) = crys(iphase)%g_0_tmp(islip) + &
          & crys(iphase)%precipitation_strength
        end do
      end do 
    else if (mesh%crss_defined .eqv. .true.) then
      ! Assign the crss based on the crystal type
      do iphase = 1, mesh%num_phases
        ! Finds the indices corresponding to the current phase the loop is on
        call find_indices(my_phase, iphase, indices, num_ind, elt_sub - 1)
        ! Initialize the crystal slip strengths corresponding to the slip number
        allocate(crys(iphase)%g_0(crys(iphase)%numslip, elt_sub:elt_sup))
        do islip = 1, crys(iphase)%numslip
          crys(iphase)%g_0(islip,:) = mesh%crss(islip, indices, 5)
          results%crss(islip, indices, :) = mesh%crss(islip, indices,:)
        end do 
      end do
      deallocate(mesh%crss)
    else 
      call par_quit("Error  :     > g_0 not properly initialized")
    end if

  end subroutine res_init_crss

end module res_init_mod2
