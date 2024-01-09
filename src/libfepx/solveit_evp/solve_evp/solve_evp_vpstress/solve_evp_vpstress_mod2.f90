! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module solve_evp_vpstress_mod2

! Routines for solving viscoplastic crystal stress equations

! Contains subroutines:
! solve_evp_vpstress: Routine which takes care of scaling and initial guesses
! solve_newton_vp: Nonlinear solution of viscoplastic crystal stress equations
! scale_up_sigm: Rescale the stress after solution is found

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use matrix_operations_mod
  use units_mod
  use parallel_mod
  use crys_type_mod2
  use stiffness_mod
  use solve_evp_vpstress_mod3

  implicit none

  public

contains

  subroutine scale_down_defr(d_vec, epseff)

    ! Rescale the deformation rate to unit size.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! d_vec: Array of deformation rates as 5-vectors (input/output)
    ! epseff: Effective deformation of these tensors (input)

    real(rk), intent(inout) :: d_vec(5, elt_sub:elt_sup)
    real(rk), intent(in) :: epseff(elt_sub:elt_sup)

    ! Locals:

    integer :: i

    !---------------------------------------------------------------------------

    !-tm from donald's version:

    do i = 1, 5
      where (epseff > 0.0d0)
        d_vec(i, :) = d_vec(i, :)/epseff
      end where
    end do

  end subroutine scale_down_defr

  subroutine find_vertex(mesh, crys, vertex, direction, plwork)

    ! Select the vertex which maximizes the plastic work.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! vertex: List of optimal vertex numbers for each grain (output)
    ! direction: Sign of the vertex (1 or -1) (output)
    ! plwork: Array of virtual plastic work for each grain and all vertices (in)

    type(mesh_type) :: mesh
    type(crys_type) :: crys (:)
    integer, intent(out) :: vertex(elt_sub:elt_sup)
    real(rk), intent(out) :: direction(elt_sub:elt_sup)
    real(rk), allocatable, intent(inout) :: plwork(:,:)

    ! Locals:

    integer, pointer :: indices(:) => null()
    real(rk) :: pa(elt_sub:elt_sup)
    integer :: i
    integer :: j
    integer :: iphase
    integer :: num_ind
    integer :: n_edge
    integer  :: my_phase(elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    vertex = 0
    direction = -1.0d0
    do iphase = 1, mesh%num_phases
      call find_indices(my_phase, iphase, indices, num_ind, elt_sub - 1)

      n_edge = crys(iphase)%numvertices
      do i = 1, n_edge
        pa(indices) = abs(plwork(i, indices))

        do j = elt_sub, elt_sup
          if (pa(j) .gt. direction(j)) then
            direction(j) = pa(j)
            vertex(j) = i - 1
          end if
        end do
      end do

      deallocate (indices)
    end do !num_phases

    ! rc 6/24/2016: Reordered loops for better memory striding

    do j = elt_sub, elt_sup
      direction(j) = plwork(vertex(j) + 1, j)/direction(j)
    end do

    ! To keep from the vertex being reused

    !  rc 6/24/2016: Reordered loops for better memory striding
    do j = elt_sub, elt_sup
      plwork(vertex(j) + 1, j) = 0.0d0
    end do

  end subroutine find_vertex

  subroutine compute_work(mesh, crys, plwork, d_vec)

    ! Compute virtual plastic work for array of deformation rates applied to
    !   vertex stresses.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! plwork:
    ! d_vec:

    type(mesh_type) :: mesh
    type(crys_type) :: crys (:)
    real(rk), allocatable, intent(out) :: plwork(:,:)
    real(rk), intent(in) :: d_vec(5, elt_sub:elt_sup)

    ! Locals:

    real(rk), pointer :: sig_fs(:, :) => null()
    integer :: my_phase(elt_sub:elt_sup)
    integer :: n_edge
    integer :: i
    integer :: j
    integer :: iphase
    integer :: maxnumvert

    !---------------------------------------------------------------------------
    call crys_maxnumvert (crys, maxnumvert)

    allocate(plwork(maxnumvert, elt_sub:elt_sup))

    ! tsh, 1/26/03
    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    plwork = 0.0d0

    do iphase = 1, mesh%num_phases
      call crys_get(crys(iphase), vertices=sig_fs)

      n_edge = crys(iphase)%numvertices
      sig_fs = sig_fs*crys(iphase)%g_0

      do i = 1, n_edge
        do j = 1, 5
          where (my_phase .eq. iphase)
            plwork(i, :) = plwork(i, :) + sig_fs(j, i)*d_vec(j, :)
          end where
        end do
      end do

      deallocate (sig_fs)
    end do !num_phases

  end subroutine compute_work

  subroutine vertex_stress(mesh, crys, sig, vertex, direction)

    !     Set initial guess (vertex stress) for nonlinear solver.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! sig: Initial guesses for stresses (output)
    ! sig_fs: Vertex stresses (input)
    ! vertex: Optimal vertex numbers
    ! direction: Sign to multiply vertex stress by

    type(mesh_type) :: mesh
    type(crys_type) :: crys (:)
    real(rk), intent(out) :: sig(5, elt_sub:elt_sup)
    integer, intent(in) :: vertex(elt_sub:elt_sup)
    real(rk), intent(in) :: direction(elt_sub:elt_sup)

    ! Locals:

    integer :: my_phase(elt_sub:elt_sup)
    real(rk), pointer :: sig_fs(:, :) => null()
    integer :: n_edge
    integer :: i
    integer :: j
    integer :: iphase

    !---------------------------------------------------------------------------

    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    do iphase = 1, mesh%num_phases
      call crys_get(crys(iphase), vertices=sig_fs)

      n_edge = crys(iphase)%numvertices
      ! rc 6/24/2016 reordered the where loops as well so it only runs the
      !   parts needed

      sig_fs = sig_fs*crys(iphase)%g_0

      ! rc 6/24/2016: Reordered loops for better memory striding
      do j = 1, n_edge
        do i = 1, 5
          where (vertex .eq. j - 1)
            where (my_phase .eq. iphase)
              sig(i, :) = sig_fs(i, j)*direction(:)
            end where
          end where
        end do
      end do

      deallocate (sig_fs)
    end do !num_phases

  end subroutine vertex_stress

  subroutine scale_stress(mesh, crys, sig, crss)

    ! Rescale stress initial guess according to hardness.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! sig: Initial guess for stress (input/output)
    ! crss: Crystal hardnesses

    type(mesh_type) :: mesh
    type(crys_type) :: crys (:)
    real(rk), intent(inout) :: sig(5, elt_sub:elt_sup)
    real(rk), intent(in) :: crss(elt_sub:elt_sup)

    ! Locals:

    integer :: my_phase(elt_sub:elt_sup)
    integer :: islip
    integer :: i
    integer :: k
    integer :: iphase
    integer :: num_ind
    integer, pointer :: indices(:) => null()
    real(rk), pointer  :: p_hat_vec(:, :) => null()
    real(rk) :: taumax(elt_sub:elt_sup)
    real(rk) :: tau(elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    !-tm  from donald's verion:

    taumax = 1.0d-8  ! deb/ prevent div by 0
    !taumax = 0.0d0

    do iphase = 1, mesh%num_phases
      call crys_get(crys(iphase), dev=p_hat_vec)

      call find_indices(my_phase, iphase, indices, num_ind, elt_sub - 1)

      do islip = 1, crys(iphase)%numslip
        call ss_project(p_hat_vec(:, islip), sig, num_ind, indices, tau)

        tau(indices) = dabs(tau(indices))

        ! maybe: where ((tau(:,indices) .gt. taumax) taumax=tau
        where ((my_phase .eq. iphase) .and. (tau .gt. taumax))
          taumax = tau
        end where
      end do

      deallocate (p_hat_vec)
      deallocate (indices)
    end do !num_phases

    tau = crss/taumax

    do k = elt_sub, elt_sup
      do i = 1, 5
        sig(i, k) = sig(i, k)*tau(k)
      end do
    end do

  end subroutine scale_stress

  subroutine solve_newton_vp(mesh, crys, exec, sig, d_vec, crss, irc, eps, converged, &
      &vp_log)

    ! Nonlinear solution of viscoplastic crystal stress equations.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! sig: Initial guess for stresses and final solution (input/output)
    ! d_vec: Deformation rate for which to solve
    ! crss: Crystal hardnesses
    ! irc: Return flag
    ! eps: Error tolerance for nonlinear solution
    ! converged: Array telling what crystals have already converged
    ! vp_log: Write viscoplastic convergence output to log files

    type(mesh_type) :: mesh
    type(crys_type) :: crys (:)
    type(exec_type), intent(in) :: exec
    real(rk), intent(inout) :: sig(5, elt_sub:elt_sup)
    real(rk), intent(in) :: d_vec(5, elt_sub:elt_sup)
    real(rk), intent(in) :: crss(elt_sub:elt_sup)
    integer, intent(out) :: irc
    real(rk), intent(in) :: eps
    logical, intent(inout) :: converged(elt_sub:elt_sup)
    logical, intent(in) :: vp_log

    ! Locals:

    logical :: newton_ok(elt_sub:elt_sup)

    integer :: iter
    integer :: inewton

    real(rk) :: res(elt_sub:elt_sup)
    real(rk) :: res_n(elt_sub:elt_sup)
    real(rk) :: fact(elt_sub:elt_sup)
    real(rk) :: ratio_res(elt_sub:elt_sup)
    real(rk) :: res_aux(elt_sub:elt_sup)
    real(rk) :: sig_0(5, elt_sub:elt_sup)
    real(rk) :: rhs(5, elt_sub:elt_sup)
    real(rk) :: rss(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: shear(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: stif(5, 5, elt_sub:elt_sup)
    real(rk) :: del_s(5, elt_sub:elt_sup)
    real(rk) :: xlambda(5, elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    irc = 0
    sig_0 = 0.0d0

    newton_ok = .true.

    do iter = 1, exec%sx_max_iters_newton
      sig_0 = sig

      call get_res(mesh, crys, res_n, rhs, rss, shear, sig, d_vec, crss)

      call form_crystif(mesh, crys, stif, rss, shear, crss)

      del_s = rhs

      call solvit(stif, del_s)

      call check_diagonals(stif, newton_ok)

      xlambda = sig_0 + del_s

      call get_res(mesh, crys, res, rhs, rss, shear, xlambda, d_vec, crss)

      fact = 1.0d0
      ratio_res = res/res_n

      do while (any(ratio_res .gt. 1.0d0 .and. newton_ok .and. .not. &
          &converged))

        where (ratio_res .gt. 1.0d0 .and. newton_ok .and. .not. converged)
          fact = fact*0.5d0
        end where

        if (any(fact .lt. 0.001d0)) then
          if (vp_log .and. (myid .eq. 0)) then
            write (*, '(a)') 'Warning:       . solve_newton_vp: &
                &Line search failure for ', count(fact .lt. 0.001d0), &
                & ' grains.'
          end if

          where (fact .lt. 0.001d0)
            newton_ok = .false.
          end where
        end if

        xlambda(1, :) = sig_0(1, :) + fact*del_s(1, :)
        xlambda(2, :) = sig_0(2, :) + fact*del_s(2, :)
        xlambda(3, :) = sig_0(3, :) + fact*del_s(3, :)
        xlambda(4, :) = sig_0(4, :) + fact*del_s(4, :)
        xlambda(5, :) = sig_0(5, :) + fact*del_s(5, :)

        call get_res(mesh, crys, res_aux, rhs, rss, shear, xlambda, d_vec, crss)

        where (ratio_res .gt. 1.0d0 .and. newton_ok .and. .not. converged)
          res = res_aux
          ratio_res = res/res_n
        end where
      end do

      where (newton_ok .and. .not. converged)
        sig(1, :) = sig_0(1, :) + fact*del_s(1, :)
        sig(2, :) = sig_0(2, :) + fact*del_s(2, :)
        sig(3, :) = sig_0(3, :) + fact*del_s(3, :)
        sig(4, :) = sig_0(4, :) + fact*del_s(4, :)
        sig(5, :) = sig_0(5, :) + fact*del_s(5, :)
      end where

      where (res .le. eps .and. newton_ok)
        converged = .true.
      end where

      ! Return if all grains have converged or solution is only an estimate.

      inewton = count(.not. newton_ok)

      if (((count(converged) + inewton) .eq. elt_sup - elt_sub + 1) .and. (myid .eq. 0)) &
          & then
        if ((inewton .gt. 0) .and. vp_log) then
          write (*, '(a)') 'Info   :     > solve_newton_vp: &
              &Converged = ', count(converged), ' remaining = ', &
              & inewton, ' minval res = ', &
              & minval(res, mask=converged), ' maxval res = ', &
              & maxval(res, mask=converged)
        end if

        irc = inewton

        return
      end if
    end do !newton_iterations

    irc = -2

  end subroutine solve_newton_vp

  !> Computes the average strength of the crystal slip systems
  subroutine compute_avg_crss(mesh, crys, crss, crss_avg)

    type(mesh_type) :: mesh
    type(crys_type) :: crys (:)
    real(rk), intent(in) :: crss(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(inout) :: crss_avg(elt_sub:elt_sup)

    integer, pointer :: indices(:) => null()
    integer :: islip
    integer :: n_slip
    integer :: iphase
    integer :: num_ind

    !---------------------------------------------------------------------------

    crss_avg = 0.0d0

    ! Goes through the total number of phases in the material
    do iphase = 1, mesh%num_phases
      n_slip = crys(iphase)%numslip

      ! Finds the indices corresponding to the current phase the loop is on

      call find_indices(mesh%elt_phase(elt_sub:elt_sup), iphase, indices, num_ind, elt_sub - 1)

      ! Sums up the crystal slip strengths corresponding to the given indices

      do islip = 1, n_slip
        crss_avg(indices) = crss_avg(indices) + crss(islip, indices)
      end do

      ! Calculates the average of these slip systems

      crss_avg(indices) = crss_avg(indices)/n_slip

      deallocate (indices)
    end do

  end subroutine compute_avg_crss

end module solve_evp_vpstress_mod2
