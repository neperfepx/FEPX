! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
!**********************************************************************
!
      SUBROUTINE check_cstress_cstrain_evps(&
  &   stif, fe, d, d_vec, elpress,&
  &   sig_vec, e_vec, c_angs, e_elas_kk, sig_kk, wts, sig_avg,&
  &   sig_avg_kk, keinv, dtime, incr, strain )
!
!----------------------------------------------------------------------
!
      USE parallel_mod
      USE READ_INPUT_MOD
      USE DimsModule
      USE IntrinsicTypesModule, RK=>REAL_KIND
      USE UtilsCrystalModule
      use KinematicsModule
      USE microstructure_mod
!
      IMPLICIT  NONE
!
!     Arguments:
!
      INTEGER incr
!
      REAL(RK), INTENT(OUT)   ::  sig_avg(0:TVEC1, el_sub1:el_sup1)
      REAL(RK), INTENT(OUT)   ::  sig_avg_kk(el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) ::  sig_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) ::  strain(0:DIMS1, 0:DIMS1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    ::  d(0:DIMS1, 0:DIMS1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    ::  d_vec(0:TVEC1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    ::  dtime
      REAL(RK), INTENT(IN)    ::  keinv(0:TVEC1,1:numphases)
      REAL(RK), INTENT(IN)    ::  stif(TVEC, TVEC, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    ::  fe(TVEC, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    ::  elpress(el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    ::  c_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    ::  e_elas_kk(el_sub1:el_sup1), sig_kk(el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    ::  wts(0:ngrain1, el_sub1:el_sup1)

      ! tm268_M11:
      REAL(RK), INTENT(OUT)   :: e_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
!
!     Locals:
!
      INTEGER  :: i, j, ij, m_el
!
      REAL(RK) :: sig_avg_chk(0:TVEC1, el_sub1:el_sup1)
      REAL(RK) :: sig_avg_kk_chk      (el_sub1:el_sup1)
      REAL(RK) :: d_eff
      REAL(RK) :: sav_chk(0:DIMS1, 0:DIMS1, el_sub1:el_sup1)
      REAL(RK) :: sav    (0:DIMS1, 0:DIMS1, el_sub1:el_sup1)
      REAL(RK) :: s_diff(6, el_sub1:el_sup1)
      REAL(RK) :: s_diff_norm(6)
      REAL(RK) :: d_mean(6)
      REAL(RK) :: d_std(6)
      REAL(RK) :: d_var(6, el_sub1:el_sup1)
      REAL(RK) :: part_sdiff_norm, whole_sdiff_norm
      REAL(RK) :: part_d_mean, whole_d_mean
!
!----------------------------------------------------------------------
!
      m_el = el_sup1 - el_sub1 + 1
!
!     -------------------
!
!     Average stresses using : [stiff]{d_vec} - {fe}
!
      sig_avg_chk = 0.0_RK
      do i  = 0, TVEC1
         do j = 0, TVEC1
            sig_avg_chk(i, :) = sig_avg_chk(i, :) + stif(i + 1, j + 1, :) * d_vec(j, :)
         enddo
         sig_avg_chk(i, :) = sig_avg_chk(i, :) - fe(i + 1, :)
      enddo
      sig_avg_kk_chk = - 3.0d0 * elpress
!
!     Average stresses using averaging procedure.
!
      call compute_cstress_evps(sig_vec, e_vec, c_angs, sig_kk, sig_avg,&
  &        sig_avg_kk, e_elas_kk, wts, keinv, ngrain, m_el)
!
!     Print stress values at element '0'.
!
      ! {5} --> [3x3]sym
      call vec_mat_symm(sig_avg_chk, sav_chk, m_el)
      call vec_mat_symm(sig_avg, sav, m_el)

!
!     Develop norm for stress differences.
!
      ij = 0
      do i = 0, DIMS1
         do j = i, DIMS1
            ij = ij + 1
            s_diff(ij, :) = (sav(i, j, :) - sav_chk(i, j, :))**2
         enddo
      enddo

      do ij = 1, 6
         part_sdiff_norm = sum(s_diff(ij, :))
         call par_sum(part_sdiff_norm, whole_sdiff_norm)
         s_diff_norm(ij) =  dsqrt(whole_sdiff_norm) / maxel
      enddo
!
!     Average deformation rate & strain history at the center.
!
      ij = 0
      do i = 0, DIMS1
         do j = i, DIMS1
            ij = ij + 1
            part_d_mean = sum( d(i, j, :) )
            call par_sum(part_d_mean, whole_d_mean)
            d_mean(ij) =  whole_d_mean/maxel
         enddo
      enddo
!
      d_eff = d_mean(1)**2 + d_mean(4)**2 + d_mean(6)**2 +&
  &        2.d0 * (d_mean(2)**2 + d_mean(3)**2 + d_mean(5)**2)
      d_eff = dsqrt(2.d0 / 3.d0 * d_eff)
!
      d_var(1, :) = d(0, 0, :) - d_mean(1)
      d_var(2, :) = d(0, 1, :) - d_mean(2)
      d_var(3, :) = d(0, 2, :) - d_mean(3)
      d_var(4, :) = d(1, 1, :) - d_mean(4)
      d_var(5, :) = d(1, 2, :) - d_mean(5)
      d_var(6, :) = d(2, 2, :) - d_mean(6)

      do ij = 1, 6
         d_std(ij) = sum( d_var(ij, :) * d_var(ij, :) )
         if (maxel .gt. 1) d_std(ij) = dsqrt(d_std(ij) / (maxel - 1))
      enddo

!     update strain

      do i = 0, DIMS1
         do j = 0, DIMS1
            strain(i, j, :) = strain(i, j, :) + dtime * d(i, j, :)
         enddo
      enddo

      RETURN
      END

!**********************************************************************
!
      SUBROUTINE compute_cstress_evps(&
  &   sig_vec, e_vec, c_angs, s_kk, s_avg, s_avg_kk, e_elas_kk, wts,&
  &   keinv, n, m )
!
!----------------------------------------------------------------------
!
      USE READ_INPUT_MOD
      use microstructure_mod
      USE UtilsCrystalModule
      USE IntrinsicTypesModule, RK=>REAL_KIND
      !
      IMPLICIT  NONE
!
!
!     Arguments:
!
      INTEGER n, m
!
      REAL(RK), INTENT(OUT)   ::  s_avg(0:TVEC1, 0:(m - 1))
      REAL(RK), INTENT(OUT)   ::  s_avg_kk(0:(m - 1))
      REAL(RK), INTENT(INOUT) ::  sig_vec(0:TVEC1, 0:(n - 1), 0:(m - 1))
                                  ! enters as kirchhoff stress (deviatoric part)
                                  ! exits as cauchy stress (deviatoric part)
                                  ! (both in crystal coords) 
      REAL(RK), INTENT(IN)    ::  c_angs(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    ::  s_kk(0:(m - 1)), e_elas_kk(0:(m - 1))
      REAL(RK), INTENT(IN)    ::  keinv(0:TVEC1,1:numphases)
      REAL(RK), INTENT(IN)    ::  wts(0:(n - 1), 0:(m - 1))
      !
      ! tm268_M11:
      REAL(RK), INTENT(OUT)   :: e_vec(0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER  :: i, j, k, iphase, numind
      INTEGER  :: my_phase(0:(m-1))
      INTEGER, POINTER :: indices(:) => NULL()
!
      REAL(RK)  ::  e_elas_kk_grn(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  s_kk_grn(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  e_elas_vec_dev(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), pointer ::    e_elas_vec_dev_tmp(:,:,:) => NULL()
      REAL(RK)  ::  v_tensor(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  determ_v(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  qr5x5(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  sig_sm(0:TVEC1, 0:(n - 1), 0:(m - 1))

!
!----------------------------------------------------------------------
!
      s_avg = 0.0d0
      s_avg_kk = 0.0d0
      my_phase(:) = phase(el_sub1:el_sup1)
!
      do iphase=1,numphases
         call find_indices(numind, iphase, my_phase, indices)
           if (associated(e_elas_vec_dev_tmp)) then
              deallocate(e_elas_vec_dev_tmp)
           endif
           allocate(e_elas_vec_dev_tmp(0:TVEC1,0:(n-1),0:(numind-1)))
         call vec_d_vec5(keinv(:,iphase), sig_vec(:,:,indices), e_elas_vec_dev_tmp, n, numind)
         e_elas_vec_dev(:,:,indices)=e_elas_vec_dev_tmp
         !
         deallocate(e_elas_vec_dev_tmp)
         deallocate(indices)
         !
      enddo !numphases


      ! tm268_M11:
      ! Deviatoric elastic strains
      e_vec = e_elas_vec_dev

      ! spread over grains
      e_elas_kk_grn = spread(e_elas_kk, dim = 1, ncopies = n)
      s_kk_grn = spread(s_kk, dim = 1, ncopies = n)
!
!     Determinant of tensor V* in lattice axes.
!     -----------------------------------------
      ! V*=I+e*
      ! 
      ! e_elas_vec_dev --> v_tensor
      ! {5} --> [3x3]sym
      call vec_mat_symm_grn(e_elas_vec_dev, v_tensor, n, m)
      do i = 0, DIMS1
         v_tensor(i, i, :, :) = v_tensor(i, i, :, :) +  e_elas_kk_grn / 3. + 1.
      enddo
      ! det(V*)
      call determinant_grn(v_tensor, determ_v, n, m)
!
!     Cauchy Stress at current configuration.
!     -----------------------------------------
      ! deviatoric part
      do i = 0,TVEC1
         do j = 0,(n - 1)
            do k = 0,(m - 1)
               sig_vec(i, j, k) = sig_vec(i, j, k) / determ_v(j, k)
            enddo
         enddo
      enddo
      ! volumetric part
      s_kk_grn = s_kk_grn / determ_v
!
!     Average stresses in sample coords
!     -----------------------------------------
      ! c_angs [3x3] --> qr5x5 [5x5] 
      call rot_mat_symm(c_angs, qr5x5, n, m)
      ! deviatoric cauhcy stress: crystal --> sample
      call mat_x_vec5(qr5x5, sig_vec, sig_sm, n, m)
      ! weighted average
      do i = 0, TVEC1
         s_avg(i, :) = sum(sig_sm(i, :, :) * wts, dim = 1)
      enddo
      s_avg_kk = sum(s_kk_grn * wts, dim = 1)

      RETURN
      END
!
!**********************************************************************
!
      SUBROUTINE qpstress_deform_evps(&
  &   c, fe, velocity, dtrace, sigv_def_qp, elem_volume )
!
!----------------------------------------------------------------------
!
      USE gather_scatter
      USE parallel_mod
      USE parallel_matrix_mod
      USE units_mod
      USE quadrature_mod
      USE READ_INPUT_MOD
      USE shape_3d_mod
      USE DimsModule
      USE IntrinsicTypesModule, RK=>REAL_KIND
      USE UtilsCrystalModule
!
      IMPLICIT  NONE
!
!     Arguments:
!
      TYPE(trace) dtrace
!
      INTEGER :: idummy
!
      REAL(RK), INTENT(IN)   ::  c(TVEC, TVEC, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)   ::  fe(TVEC, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)   ::  velocity(0:maxdof1)
      REAL(RK), INTENT(OUT)  ::  sigv_def_qp(0:TVEC1, 0:nqpt1, el_sub1:el_sup1)
      REAL(RK), INTENT(OUT)  ::  elem_volume(el_sub1:el_sup1)
!
!     Locals:
!
      INTEGER   i, j, k, l, iqpt, ier, idug, m_el
!
      REAL(RK)&
  &   gcoords(0:kdim1, el_sub1:el_sup1),&
  &   gvel   (0:kdim1, el_sub1:el_sup1)
      REAL(RK)&
  &   dndx(0:nnpe, el_sub1:el_sup1), dndy(0:nnpe, el_sub1:el_sup1),&
  &   dndz(0:nnpe, el_sub1:el_sup1), det(el_sub1:el_sup1)
      REAL(RK)&
  &   loc0, loc1,&
  &   loc2
      REAL(RK)&
  &   s11(el_sub1:el_sup1), s12(el_sub1:el_sup1), s13(el_sub1:el_sup1),&
  &   s21(el_sub1:el_sup1), s22(el_sub1:el_sup1), s23(el_sub1:el_sup1),&
  &   s31(el_sub1:el_sup1), s32(el_sub1:el_sup1), s33(el_sub1:el_sup1)
      REAL(RK)&
  &   d(0:DIMS1, 0:DIMS1, el_sub1:el_sup1),&
  &   d_vec     (0:TVEC1, el_sub1:el_sup1)
      REAL(RK)  sigv_def(0:TVEC1, el_sub1:el_sup1)
      REAL(RK)  wt
!
!----------------------------------------------------------------------
!
      m_el = el_sup1 - el_sub1 + 1
!
      elem_volume = 0.0d0
!
!     Evaluate the stresses at each quadrature point using constitutive law.
!
      call part_gather(gcoords, coords, nodes, dtrace)
      call part_gather(gvel, velocity, nodes, dtrace)
!
      iqpt = 0
      do i = 0, nqpt1
         loc0 = qploc(0, i)
         loc1 = qploc(1, i)
         loc2 = qploc(2, i)

         call sfder_hpar(loc0, loc1, loc2, gcoords, dndx, dndy,&
  &      dndz, det, s11, s12, s13, s21, s22, s23, s31, s32, s33)

         call defrate(d, dndx, dndy, dndz, gvel)

         call mat_vec_symm(d, d_vec, m_el)

         sigv_def = 0.0d0

         call gen_matrix_vector_mult(sigv_def, c, d_vec,&
  &      idummy, idummy, idummy, idummy, ier)

         if ( ier .ne. 0 )&
  &      call par_quit('qpstress_deform_evps: matrix/vector multipy error.')

         do l = 0, TVEC1
           sigv_def_qp(l, iqpt, :) = sigv_def(l, :) - fe(l + 1, :)
         enddo
!
!        Contribution to averaged stress over domain
!

         wt = wtqp(0, i)
         elem_volume = elem_volume + det * wt

         iqpt = iqpt + 1

      end do !nqpt1

      RETURN
      END
!
