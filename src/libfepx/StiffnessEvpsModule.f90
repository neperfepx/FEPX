! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE StiffnessEvpsModule

    USE parallel_mod
    USE parallel_matrix_mod

    USE IntrinsicTypesModule, RK=>REAL_KIND
    USE DimsModule
    USE units_mod

    USE quadrature_mod
    USE READ_INPUT_MOD
    USE shape_3d_mod
    use microstructure_mod
    USE MATRIX_OPERATIONS_MOD
    USE MaterialMatrixEvpsModule
    use StiffnessVpModule

    IMPLICIT NONE

CONTAINS
  !
  SUBROUTINE element_stif_evps(&
       &    gstiff, gtanstiff, f_vec, gcoords, gvel,&
       &    c0_angs, c_angs, sig_vec_n, sig_vec, crss_n,&
       &    crss, rstar_n, rstar,  &
       &    e_bar_vec, e_elas_kk_bar, e_elas_kk, sig_kk, jiter_state, &
       &    keinv, wts, epseff, &
       &    dtime, incr, &
       &    converged_solution, auto_time, NR)
    !
    !----------------------------------------------------------------------
    !
    !     Arguments:
    !
    LOGICAL, INTENT(INOUT) :: converged_solution
    !
    REAL(RK) ::  dtime
    REAL(RK) ::  bulk_fac1(el_sub1:el_sup1), bulk_fac2(el_sub1:el_sup1)
    !
    INTEGER  :: incr, auto_time
    INTEGER, INTENT(OUT)  :: jiter_state(0:ngrain1, el_sub1:el_sup1)    
    !
    REAL(RK), INTENT(IN)    :: keinv(0:TVEC1,1:numphases)
    REAL(RK), INTENT(INOUT) :: gstiff(0:kdim1, 0:kdim1, el_sub1:el_sup1)  !enters=0
    REAL(RK), INTENT(INOUT) :: gtanstiff(0:kdim1, 0:kdim1, el_sub1:el_sup1)  !enters=0
    REAL(RK), INTENT(OUT)   :: e_elas_kk_bar(el_sub1:el_sup1, 0:nqpt1)
    REAL(RK), INTENT(OUT)   :: e_elas_kk(el_sub1:el_sup1, 0:nqpt1)
    REAL(RK), INTENT(OUT)   :: sig_kk(el_sub1:el_sup1, 0:nqpt1)
    REAL(RK), INTENT(OUT)   :: epseff(el_sub1:el_sup1, 0:nqpt1)
    REAL(RK), INTENT(IN)    :: gcoords(0:kdim1, el_sub1:el_sup1)
    REAL(RK), INTENT(IN)    :: gvel   (0:kdim1, el_sub1:el_sup1)
    REAL(RK), INTENT(OUT)   :: f_vec  (0:kdim1, el_sub1:el_sup1)
    REAL(RK), INTENT(IN)    :: c0_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(INOUT) :: c_angs (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
    REAL(RK), INTENT(IN)    :: rstar_n(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(INOUT) :: rstar  (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
    REAL(RK), INTENT(OUT)   :: sig_vec_n(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
    REAL(RK), INTENT(OUT)   :: sig_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
    REAL(RK), INTENT(OUT)   :: e_bar_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
    REAL(RK), INTENT(IN)    :: crss_n(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(OUT)   :: crss  (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1,0:nqpt1)
    REAL(RK), INTENT(IN)    :: wts(0:ngrain1, el_sub1:el_sup1)
    LOGICAL,  INTENT(IN)    :: NR
    !
    !
    !     Locals:
    !
    INTEGER  ::  iqpt, jqpt, kqpt, m_el, iphase, numind
    INTEGER, POINTER :: indices(:) => NULL()
    REAL(RK), pointer :: e_bar_vec_tmp(:,:,:) => NULL()
    INTEGER  ::  i, j, k, l, i1, i2, i3, j1, j2, j3, ier
    !
    REAL(RK) :: wt
    !
    REAL(RK) :: dndx(0:nnpe, el_sub1:el_sup1, 0:nqpt1)
    REAL(RK) :: dndy(0:nnpe, el_sub1:el_sup1, 0:nqpt1)
    REAL(RK) :: dndz(0:nnpe, el_sub1:el_sup1, 0:nqpt1)
    !
    REAL(RK) :: det (el_sub1:el_sup1, 0:nqpt1)
    REAL(RK) :: detv(el_sub1:el_sup1, 0:nqpt1)
    REAL(RK) :: loc0
    REAL(RK) :: loc1
    REAL(RK) :: loc2    
    !
    REAL(RK) :: s11(el_sub1:el_sup1), s12(el_sub1:el_sup1), s13(el_sub1:el_sup1)
    REAL(RK) :: s21(el_sub1:el_sup1), s22(el_sub1:el_sup1), s23(el_sub1:el_sup1)
    REAL(RK) :: s31(el_sub1:el_sup1), s32(el_sub1:el_sup1), s33(el_sub1:el_sup1)
    !
    REAL(RK) :: c(TVEC, TVEC, el_sub1:el_sup1, 0:nqpt1)
    REAL(RK) :: c_tan(TVEC, TVEC, el_sub1:el_sup1, 0:nqpt1)
    REAL(RK) :: f(TVEC, el_sub1:el_sup1, 0:nqpt1)
    !
    REAL(RK) :: xni  (3, 5, el_sub1:el_sup1)
    REAL(RK) :: xnj  (5, 3, el_sub1:el_sup1)
    REAL(RK) :: temp1(3, 5, el_sub1:el_sup1)
    REAL(RK) :: temp2(3, 3, el_sub1:el_sup1)
    REAL(RK) :: temp3(3, 5, el_sub1:el_sup1)
    REAL(RK) :: temp4(3, 3, el_sub1:el_sup1)
    !
    REAL(RK) :: ftemp(3, el_sub1:el_sup1)
    !
    INTEGER :: my_phase(el_sub1:el_sup1)
    !
    !----------------------------------------------------------------------
    !
    m_el = el_sup1 - el_sub1 + 1
    !
    my_phase(:) = phase(el_sub1:el_sup1)
    !
    !
    f_vec = 0.0_RK
    !
    !----- quadrature points -----

    !
    !prd
    ! update internal variables at QPs.
    sig_vec_n = gsig_vec_n
    e_elas_kk_bar = gela_kk_bar(0,:,:)


    do iqpt = 0, nqpt1

       ! coordinates in the parent element
       loc0 = qploc(0, iqpt)
       loc1 = qploc(1, iqpt)
       loc2 = qploc(2, iqpt)
       ! weigth
       wt = wtqp(0, iqpt)

       ! This subroutine is moved here (tsh version)
       ! compute elastic strain: sig_vec_n --> e_bar_vec
       do iphase=1,numphases
          call find_indices(numind, iphase, my_phase, indices)
          indices=indices+el_sub1
          if (associated(e_bar_vec_tmp)) then
             deallocate(e_bar_vec_tmp)
          endif
          allocate(e_bar_vec_tmp(0:TVEC1,0:ngrain1,0:(numind-1)))
          call vec_d_vec5(keinv(:,iphase), sig_vec_n(:,:,indices,iqpt), e_bar_vec_tmp, ngrain, numind)
          e_bar_vec(:,:,indices,iqpt)=e_bar_vec_tmp
          deallocate(e_bar_vec_tmp)
          deallocate(indices)
       enddo !numphases

       ! Compute quadrature quantities given a set of local coordinates.
       call sfder_hpar(loc0, loc1, loc2, gcoords, dndx(:,:,iqpt), dndy(:,:,iqpt),  &
            &     dndz(:,:,iqpt), det(:,iqpt), s11, s12, s13, s21, s22, s23, s31, s32,  &
            &     s33)

    enddo

   ! IN: e_bar_vec
   ! tm268_M13: removed eqplas
   call material_matrix_evps(c, c_tan, f, detv, dndx, dndy, dndz, gvel, &
        &     c0_angs, c_angs, sig_vec_n, sig_vec, crss_n,&
        &     crss, rstar_n, rstar, keinv,&
        &     e_bar_vec, e_elas_kk_bar, e_elas_kk, sig_kk, jiter_state, &
        &     wts, epseff, &
        &     dtime, incr,&
        &     converged_solution, auto_time, NR)

   if (.not. converged_solution .and. auto_time .eq. 1) RETURN

    do iqpt = 0, nqpt1

       do iphase = 1,numphases
          where (my_phase == iphase)
             bulk_fac1 = crystal_parm(8,iphase) / detv(:, iqpt) * dtime
             bulk_fac2 = crystal_parm(8,iphase) / detv(:, iqpt) * e_elas_kk_bar(:, iqpt)
          endwhere
       enddo

       if (minval(det) .lt. 0.0d0 ) then
          !
          do i = el_sub1, el_sup1
             if ( det(i,iqpt) .lt. 0.0 ) then
                WRITE(DFLT_U, *) 'Error  :       . Element: ',i,', Determinant: ',det(i,iqpt)
             endif
          end do
          !
          CALL par_quit('Error  :       . ELEMENT_STIF_EVPS: Negative Jacobian(s)', ABORT=.TRUE.)
       endif

       det(:, iqpt) = det(:, iqpt) * wtqp(0, iqpt)

       do i = 0, nnpe
          i1 = 3 * i
          i2 = i1 + 1
          i3 = i2 + 1

          xni(1, 1, :) = dndx(i, :,iqpt)
          xni(2, 1, :) = -dndy(i, :,iqpt)
          xni(3, 1, :) = 0.0
          xni(1, 2, :) = -dndx(i, :,iqpt) / 3.0
          xni(2, 2, :) = -dndy(i, :,iqpt) / 3.0
          xni(3, 2, :) = 2.0 * dndz(i, :,iqpt) / 3.0
          xni(1, 3, :) = 0.5 * dndy(i, :,iqpt)
          xni(2, 3, :) = 0.5 * dndx(i, :,iqpt)
          xni(3, 3, :) = 0.0
          xni(1, 4, :) = 0.5 * dndz(i, :,iqpt)
          xni(2, 4, :) = 0.0
          xni(3, 4, :) = 0.5 * dndx(i, :,iqpt)
          xni(1, 5, :) = 0.0
          xni(2, 5, :) = 0.5 * dndz(i, :,iqpt)
          xni(3, 5, :) = 0.5 * dndy(i, :,iqpt)

          ftemp = 0.0_RK
!  RC 6/24/2016: Reordered for better memory striding
          do l = 1, 5
             do k = 1, 3
                ftemp(k, :) = ftemp(k, :) + xni(k, l, :) * f(l, :, iqpt)
             enddo
          enddo

          ftemp(1, :) = ftemp(1, :) - dndx(i, :, iqpt) * bulk_fac2
          ftemp(2, :) = ftemp(2, :) - dndy(i, :, iqpt) * bulk_fac2
          ftemp(3, :) = ftemp(3, :) - dndz(i, :, iqpt) * bulk_fac2

          f_vec(i1, :) = f_vec(i1, :) + ftemp(1, :) * det(:,iqpt)
          f_vec(i2, :) = f_vec(i2, :) + ftemp(2, :) * det(:,iqpt)
          f_vec(i3, :) = f_vec(i3, :) + ftemp(3, :) * det(:,iqpt)

          temp1 = 0.0_RK
          temp3 = 0.0_RK
          call gen_matrix_mult(temp1, xni, c(:,:,:,iqpt), 1, 2, ier)
          if (NR) call gen_matrix_mult(temp3, xni, c_tan(:,:,:,iqpt), 1, 2, ier)

          do j = 0, i
             j1 = 3 * j
             j2 = j1 + 1
             j3 = j2 + 1

             xnj(1, 1, :) = dndx(j, :, iqpt)
             xnj(1, 2, :) = -dndy(j, :, iqpt)
             xnj(1, 3, :) = 0.0
             xnj(2, 1, :) = -dndx(j, :, iqpt) / 3.0
             xnj(2, 2, :) = -dndy(j, :, iqpt) / 3.0
             xnj(2, 3, :) = 2.0 * dndz(j, :, iqpt)  /3.0
             xnj(3, 1, :) = 0.5 * dndy(j, :, iqpt)
             xnj(3, 2, :) = 0.5 * dndx(j, :, iqpt)
             xnj(3, 3, :) = 0.0
             xnj(4, 1, :) = 0.5 * dndz(j, :, iqpt)
             xnj(4, 2, :) = 0.0
             xnj(4, 3, :) = 0.5 * dndx(j, :, iqpt)
             xnj(5, 1, :) = 0.0
             xnj(5, 2, :) = 0.5 * dndz(j, :, iqpt)
             xnj(5, 3, :) = 0.5 * dndy(j, :, iqpt)

             temp2 = 0.0_RK
             temp4 = 0.0_RK
             call gen_matrix_mult(temp2, temp1, xnj, 1, 2, ier)
             if (NR) call gen_matrix_mult(temp4, temp3, xnj, 1, 2, ier)

             ! Assemble secant + volumentric stiffness

             s11 = temp2(1, 1, :) * det(:, iqpt)&
                  &         + dndx(i, :, iqpt) * dndx(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
             s12 = temp2(1, 2, :) * det(:, iqpt)&
                  &         + dndx(i, :, iqpt) * dndy(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
             s13 = temp2(1, 3, :) * det(:, iqpt)&
                  &         + dndx(i, :, iqpt) * dndz(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
             s21 = temp2(2, 1, :) * det(:, iqpt)&
                  &         + dndy(i, :, iqpt) * dndx(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
             s22 = temp2(2, 2, :) * det(:, iqpt)&
                  &         + dndy(i, :, iqpt) * dndy(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
             s23 = temp2(2, 3, :) * det(:, iqpt)&
                  &         + dndy(i, :, iqpt) * dndz(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
             s31 = temp2(3, 1, :) * det(:, iqpt)&
                  &         + dndz(i, :, iqpt) * dndx(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
             s32 = temp2(3, 2, :) * det(:, iqpt)&
                  &         + dndz(i, :, iqpt) * dndy(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
             s33 = temp2(3, 3, :) * det(:, iqpt)&
                  &         + dndz(i, :, iqpt) * dndz(j, :, iqpt) * bulk_fac1 * det(:, iqpt)

             call add_to_stiffness(gstiff, s11, i1, j1)
             call add_to_stiffness(gstiff, s22, i2, j2)
             call add_to_stiffness(gstiff, s33, i3, j3)
             call add_to_stiffness(gstiff, s12, i1, j2)
             call add_to_stiffness(gstiff, s13, i1, j3)
             call add_to_stiffness(gstiff, s21, i2, j1)
             call add_to_stiffness(gstiff, s23, i2, j3)
             call add_to_stiffness(gstiff, s31, i3, j1)
             call add_to_stiffness(gstiff, s32, i3, j2)

             ! Assemble tangent + volumentric stiffness
             if (NR) then
                s11 = temp4(1, 1, :) * det(:, iqpt)&
                     &         + dndx(i, :, iqpt) * dndx(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
                s12 = temp4(1, 2, :) * det(:, iqpt)&
                     &         + dndx(i, :, iqpt) * dndy(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
                s13 = temp4(1, 3, :) * det(:, iqpt)&
                     &         + dndx(i, :, iqpt) * dndz(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
                s21 = temp4(2, 1, :) * det(:, iqpt)&
                     &         + dndy(i, :, iqpt) * dndx(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
                s22 = temp4(2, 2, :) * det(:, iqpt)&
                     &         + dndy(i, :, iqpt) * dndy(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
                s23 = temp4(2, 3, :) * det(:, iqpt)&
                     &         + dndy(i, :, iqpt) * dndz(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
                s31 = temp4(3, 1, :) * det(:, iqpt)&
                     &         + dndz(i, :, iqpt) * dndx(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
                s32 = temp4(3, 2, :) * det(:, iqpt)&
                     &         + dndz(i, :, iqpt) * dndy(j, :, iqpt) * bulk_fac1 * det(:, iqpt)
                s33 = temp4(3, 3, :) * det(:, iqpt)&
                     &         + dndz(i, :, iqpt) * dndz(j, :, iqpt) * bulk_fac1 * det(:, iqpt)

                call add_to_stiffness(gtanstiff, s11, i1, j1)
                call add_to_stiffness(gtanstiff, s22, i2, j2)
                call add_to_stiffness(gtanstiff, s33, i3, j3)
                call add_to_stiffness(gtanstiff, s12, i1, j2)
                call add_to_stiffness(gtanstiff, s13, i1, j3)
                call add_to_stiffness(gtanstiff, s21, i2, j1)
                call add_to_stiffness(gtanstiff, s23, i2, j3)
                call add_to_stiffness(gtanstiff, s31, i3, j1)
                call add_to_stiffness(gtanstiff, s32, i3, j2)
             end if ! NR

          enddo
       enddo
    enddo! nqpt1

    END SUBROUTINE element_stif_evps

END MODULE StiffnessEvpsModule
