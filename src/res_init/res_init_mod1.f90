! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module res_init_mod

  use general_mod
  use types_mod
  use utils_mod
  use res_init_mod2

  implicit none

  public

contains

  subroutine res_init (mesh, crys, loading, printing, results)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys (:)
    type(loading_type), intent(in) :: loading
    type(printing_type), intent(in) :: printing
    type(results_type), intent(out) :: results
    integer :: i

    call res_qpt_alloc (mesh, results)

    results%ori = spread(mesh%ori, 4, nqpt)

    do i = 1, 3
      results%d_rstar(i, i, :, :) = 1.0d0
      results%rstar(i, i, :, :) = 1.0d0
    end do

    call res_init_crss(mesh, crys, results, mesh%elt_phase(elt_sub:elt_sup))

    if (printing%print_rotrate_slip) then
      call alloc_4d (results%d_rstar_slip, 1, 3, 1, 3, elt_sub, elt_sup, 1, nqpt)
    end if

    if (printing%print_rotrate_spin) then
      call alloc_4d (results%d_rstar_spin, 1, 3, 1, 3, elt_sub, elt_sup, 1, nqpt)
    end if

    if (printing%print_slip) then
      call alloc_2d (results%slip, 1, mesh%maxnumslip, elt_sub, elt_sup)
    end if

    if (printing%print_strain) then
      if (.not. allocated(results%defrate)) then
        call alloc_3d (results%defrate, 1, 3, 1, 3, elt_sub, elt_sup)
      end if

      call alloc_3d (results%strain, 1, 3, 1, 3, elt_sub, elt_sup)
    end if

    if (printing%print_strain_pl) then
      call alloc_2d (results%strain_pl, 1, 5, elt_sub, elt_sup)
    end if

    if (printing%print_defrate) then
      call alloc_3d (results%defrate, 1, 3, 1, 3, elt_sub, elt_sup)
    end if

    if ((printing%print_work) .or. (printing%print_workrate) &
        & .or. printing%print_restart) then
      call alloc_3d (results%defrate, 1, 3, 1, 3, elt_sub, elt_sup, "soft")
      call alloc_1d (results%workrate, elt_sub, elt_sup)
      call alloc_1d (results%work, elt_sub, elt_sup)
    end if

    if ((printing%print_work_pl) .or. (printing%print_workrate_pl) &
        & .or. printing%print_restart) then
      call alloc_1d (results%workrate_pl, elt_sub, elt_sup)
      call alloc_1d (results%work_pl, elt_sub, elt_sup)
    end if

    ! handling dependencies

    if (printing%print_strain_eq) then
      call alloc_3d (results%strain, 1, 3, 1, 3, elt_sub, elt_sup, "soft")
    end if

    if (printing%print_strain_el_eq) then
      call alloc_3d (results%elas_tot6, 1, 6, elt_sub, elt_sup, 1, nqpt, "soft")
    end if

    if (printing%print_strain_pl_eq) then
      call alloc_2d (results%strain_pl, 1, 5, elt_sub, elt_sup, "soft")
    end if

    if (allocated(results%work) .or. allocated(results%workrate) &
        & .or. printing%print_restart) then
      call alloc_3d (results%defrate, 1, 3, 1, 3, elt_sub, elt_sup, "soft")
      call alloc_1d (results%workrate, elt_sub, elt_sup, "soft")
    end if

    if (allocated(results%work_pl) .or. (allocated(results%workrate_pl) &
        & .or. printing%print_restart)) then
      call alloc_1d (results%workrate_pl, elt_sub, elt_sup, "soft")
      call alloc_1d (results%work_pl, elt_sub, elt_sup)
    end if

    call res_qpt_init (mesh, loading, results)

  end subroutine res_init

  subroutine res_qpt_alloc (mesh, results)

    type(mesh_type), intent(in) :: mesh
    type(results_type), intent(out) :: results

    call alloc_1d (results%coo, dof_sub, dof_sup)
    call alloc_1d (results%vel, dof_sub, dof_sup)
    call alloc_1d (results%force, dof_sub, dof_sup)

    call alloc_4d (results%ori, 1, 3, 1, 3, elt_sub, elt_sup, 1, nqpt)
    call alloc_4d (results%rstar, 1, 3, 1, 3, elt_sub, elt_sup, 1, nqpt)
    call alloc_4d (results%d_rstar, 1, 3, 1, 3, elt_sub, elt_sup, 1, nqpt)
    call alloc_3d (results%sliprate, 1, mesh%maxnumslip, elt_sub, elt_sup, 1, nqpt)
    call alloc_3d (results%rss, 1, mesh%maxnumslip, elt_sub, elt_sup, 1, nqpt)
    call alloc_3d (results%crss, 1, mesh%maxnumslip, elt_sub, elt_sup, 1, nqpt)
    call alloc_3d (results%acmslip, 1, mesh%maxnumslip, elt_sub, elt_sup, 1, nqpt)
    call alloc_3d (results%sig_vec, 1, 5, elt_sub, elt_sup, 1, nqpt)
    call alloc_3d (results%e_vec, 1, 5, elt_sub, elt_sup, 1, nqpt)
    call alloc_2d (results%e_elas_kk_bar, elt_sub, elt_sup, 1, nqpt)
    call alloc_2d (results%d_kk, elt_sub, elt_sup, 1, nqpt)
    call alloc_2d (results%sig_kk, elt_sub, elt_sup, 1, nqpt)
    call alloc_2d (results%defrate_eq, elt_sub, elt_sup, 1, nqpt)

    call alloc_4d (results%velgrad, 1, 3, 1, 3, elt_sub, elt_sup, 1, nqpt)
    call alloc_4d (results%d, 1, 3, 1, 3, elt_sub, elt_sup, 1, nqpt)
    call alloc_4d (results%w, 1, 3, 1, 3, elt_sub, elt_sup, 1, nqpt)

    ! central qpt only
    call alloc_3d (results%elas_tot6, 1, 6, elt_sub, elt_sup, 1, nqpt)
    call alloc_3d (results%stress, 1, 3, 1, 3, elt_sub, elt_sup)

    call alloc_1d (results%work, elt_sub, elt_sup)
    if (allocated(results%work)) then
      call alloc_1d (results%work, elt_sub, elt_sup)
    end if
    if (allocated(results%work_pl)) then
      call alloc_1d (results%work_pl, elt_sub, elt_sup)
    end if

    call alloc_3d (results%dp_hat, 1, 5, elt_sub, elt_sup, 1, nqpt)
    call alloc_3d (results%wp_hat, 1, 3, elt_sub, elt_sup, 1, nqpt)
    call alloc_3d (results%wp_ss, 1, 3, elt_sub, elt_sup, 1, nqpt)

    call alloc_3d (results%e_bar_vec, 1, 5, elt_sub, elt_sup, 1, nqpt)

    call alloc_3d (results%d_vec, 1, 5, elt_sub, elt_sup, 1, nqpt)
    call alloc_3d (results%w_vec, 1, 3, elt_sub, elt_sup, 1, nqpt)
    call alloc_3d (results%w_vec_lat, 1, 3, elt_sub, elt_sup, 1, nqpt)
    call alloc_2d (results%detv, elt_sub, elt_sup, 1, nqpt)

  end subroutine res_qpt_alloc

  subroutine res_qpt_init_prev (results, results_prev)

    type(results_type), intent(in) :: results
    type(results_type), intent(inout) :: results_prev

    results_prev%coo = results%coo
    results_prev%vel = results%vel

    results_prev%ori = results%ori
    results_prev%rstar = results%rstar
    results_prev%d_rstar = results%d_rstar
    if (allocated(results%d_rstar_spin)) then
      results_prev%d_rstar_spin = results%d_rstar_spin
    end if
    if (allocated(results%d_rstar_slip)) then
      results_prev%d_rstar_slip = results%d_rstar_slip
    end if
    results_prev%crss = results%crss
    results_prev%acmslip = results%acmslip
    results_prev%sig_vec = results%sig_vec
    results_prev%e_vec = results%e_vec
    results_prev%e_elas_kk_bar = results%e_elas_kk_bar
    results_prev%d_kk = results%d_kk
    results_prev%sig_kk = results%sig_kk
    results_prev%defrate_eq = results%defrate_eq

    results_prev%velgrad = results%velgrad
    results_prev%d = results%d
    results_prev%w = results%w

    results_prev%work = results%work

    if (allocated(results%work_pl)) then
      results_prev%work_pl = results%work_pl
    end if

    if (allocated(results%slip)) then
      results_prev%slip = results%slip
    end if

  end subroutine res_qpt_init_prev

  subroutine res_qpt_init (mesh, loading, results)

    type(mesh_type), intent(in) :: mesh
    type(loading_type), intent(in) :: loading
    type(results_type), intent(inout) :: results

    results%coo = mesh%coo
    if (mesh%num_periodicity .eq. 0) then
      results%vel = loading%bcs_vel
    else
      results%vel = loading%offset_ps
    end if
    results%force = 0.0d0

  end subroutine res_qpt_init

end module res_init_mod
