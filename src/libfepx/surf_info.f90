! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE surf_info_mod
  !
  ! Maintain surface information.
  !
  USE parallel_mod
  USE shape_2_mod
  USE quadrature_mod
  USE surface_mod
  USE IntrinsicTypesModule, RK=>REAL_KIND
  !
  IMPLICIT NONE
  !
  ! Globals.
  !
  INTEGER, PARAMETER:: MAXSURFACES=6
  TYPE(surface_section) :: surfaces(MAXSURFACES)
  INTEGER :: nsurfaces
  CHARACTER(LEN=12) :: faset(6)

  ! 3 node triangular element
  !INTEGER, PARAMETER :: SFTYPE=3
  ! 6 node triangular element
  INTEGER, PARAMETER :: SFTYPE=6
  ! 4 node  quadralateral element
  !INTEGER, PARAMETER :: SFTYPE=4
  REAL(RK) :: sfqp2d(SFTYPE, MAXQP2D), sfgqp2d(NDIM2, SFTYPE, MAXQP2D)
  !
CONTAINS
  !
  SUBROUTINE upd_surf(unit, coords, sig_all, load_out, area_out)
  !
  !-----------------------------------------------------------------------  
    !
    ! Update surface information.
    !
    IMPLICIT NONE
    !
    !   Arguments:
    !
    INTEGER :: unit !--> DEBUG_U
    REAL(RK), INTENT(IN) :: coords(:), sig_all (:)
    REAL(RK), INTENT(INOUT) :: load_out(nsurfaces,3), area_out(nsurfaces) ! these variables enter as 0
    !
    !   Locals:
    !
    INTEGER :: i, iel, jq, jn, jdof, iload, j
    REAL(RK) :: tangent(3, 2, MAXQP2D), normal(3, MAXQP2D), nmag, sjac(MAXQP2D)
    REAL(RK), ALLOCATABLE:: sig_se(:,:)
    REAL(RK):: load(3), p_total_load(3), total_load(3), p_area, area
    !
    !----------------------------------------------------------------------
    !
    !write(unit,*) 'Updating surface info.'
    !
    !
    load_out=0.0_RK
    area_out=0.0_RK
    !    
    !
    do i=1, nsurfaces
       call part_gather(surfaces(i)%crds, coords,surfaces(i)%conn3d, surfaces(i)%tr)
       ALLOCATE(sig_se(0:5,surfaces(i)%semin:surfaces(i)%semax))
       call part_gather(sig_se, sig_all, surfaces(i)%econn, surfaces(i)%etr)
       !
       !  Now compute normal vectors at quadrature points.
       !
       p_total_load = 0.0d0
       p_area = 0.0d0
       !
       do iel=surfaces(i)%semin, surfaces(i)%semax
          !
          tangent = 0.0d0
          do jq=1, nqp2d
             do jn=1, SFTYPE !number of nodes in the surface element
                jdof = 3*(jn-1)
                !
                !  First tangent vector.
                !  
                tangent(1, 1, jq) = tangent(1, 1, jq) +&
                     & surfaces(i)%crds(jdof,iel)  *sfgqp2d(1,jn, jq)
                tangent(2, 1, jq) = tangent(2, 1, jq) +&
                     & surfaces(i)%crds(jdof+1,iel)*sfgqp2d(1,jn, jq)
                tangent(3, 1, jq) = tangent(3, 1, jq) +&
                     & surfaces(i)%crds(jdof+2,iel)*sfgqp2d(1,jn, jq)
                !
                !  Second tangent vector.
                !  
                tangent(1, 2, jq) = tangent(1, 2, jq) +&
                     & surfaces(i)%crds(jdof,iel)  *sfgqp2d(2,jn, jq)
                tangent(2, 2, jq) = tangent(2, 2, jq) +&
                     & surfaces(i)%crds(jdof+1,iel)*sfgqp2d(2,jn, jq)
                tangent(3, 2, jq) = tangent(3, 2, jq) +&
                     & surfaces(i)%crds(jdof+2,iel)*sfgqp2d(2,jn, jq)
             end do
             !
             !  Now compute normals.
             !
             normal(1, jq) =&
                  & tangent(2,1,jq)*tangent(3,2,jq) -&
                  & tangent(3,1,jq)*tangent(2,2,jq)
             normal(2, jq) =&
                  & tangent(3,1,jq)*tangent(1,2,jq) -&
                  & tangent(1,1,jq)*tangent(3,2,jq)
             normal(3, jq) =&
                  & tangent(1,1,jq)*tangent(2,2,jq) -&
                  & tangent(2,1,jq)*tangent(1,2,jq)
             nmag = dsqrt(&
                  & normal(1, jq)*normal(1, jq) +&
                  & normal(2, jq)*normal(2, jq) +&
                  & normal(3, jq)*normal(3, jq)  )
             if (nmag > 0) then
                normal(:,jq) = normal(:,jq)/nmag
                sjac(jq) = nmag
             else
                call par_quit('Error  :     > Surface normal zero magnitude.')
             end if
          end do
          !write(unit, *) 'element: ', iel
          !write(unit, '(3e14.4)') normal
!
          load = 0.0d0
          do jq=1, nqp2d
             load(1) = load(1) + wt2d(jq)*sjac(jq)*&
                  & (sig_se(0,iel)*normal(1,jq) +&
                  & sig_se(1,iel)*normal(2,jq) +&
                  & sig_se(2,iel)*normal(3,jq) )
             load(2) = load(2) + wt2d(jq)*sjac(jq)*&
                  & (sig_se(1,iel)*normal(1,jq) +&
                  & sig_se(3,iel)*normal(2,jq) +&
                  & sig_se(4,iel)*normal(3,jq) )
             load(3) = load(3) + wt2d(jq)*sjac(jq)*&
                  & (sig_se(2,iel)*normal(1,jq) +&
                  & sig_se(4,iel)*normal(2,jq) +&
                  & sig_se(5,iel)*normal(3,jq) )
             p_area = p_area + wt2d(jq)*sjac(jq)
          end do
          p_total_load = p_total_load + load
       end do !iel

       do iload=1,3
          call par_sum(p_total_load(iload), total_load(iload))
       end do
       call par_sum(p_area, area)
       !if (myid .eq. 0) then
       !  write(unit, '(a9,i3,a1)', ADVANCE='NO') 'load surf', i, ':'
       !  write(unit, '(5x,3e14.4)') total_load
       !  write(unit, *) 'area = ', area
       !endif
       !
       load_out(i,:)=total_load
       area_out(i)=area
       !
       DEALLOCATE(sig_se)
       !
    end do !nsurfaces
    !
    RETURN
    !
  END SUBROUTINE upd_surf
  !
  !
  !**************************************************************************
  !  
  SUBROUTINE init_surf()
  !
  !--------------------------------------------------------------------------  
    !
    ! Initialize surface arrays -- shape functions at qp.
    !
    IMPLICIT NONE
    !
    !   Arguments:
    !

    !
    !   Locals:
    !
    INTEGER :: istat
    !
    !----------------------------------------------------------------------
    !
    call sf2d(SFTYPE, nqp2d, qp2d, sfqp2d, SFTYPE, istat)
    call sf2dg(SFTYPE, nqp2d, qp2d, sfgqp2d, SFTYPE, istat)
    !
    RETURN
    !
  END SUBROUTINE init_surf
  !
  !
  !***********************************************************************
  !
  SUBROUTINE compute_area(coords, area0)
  !
  !-----------------------------------------------------------------------  
    !
    ! compute surface area
    !
    IMPLICIT NONE
    !
    !   Arguments:
    !
    REAL(RK), INTENT(IN) :: coords(:)
    REAL(RK), INTENT(INOUT) :: area0(nsurfaces)  ! enters =0
    !
    !   Locals:
    !
    INTEGER :: i, iel, jq, jn, jdof
    REAL(RK) :: tangent(3, 2, MAXQP2D), normal(3, MAXQP2D), nmag, sjac(MAXQP2D)
    REAL(RK) :: p_area, area
    !
    !----------------------------------------------------------------------
    !
    area0=0.0_RK
    !
    !   
    do i=1, nsurfaces
       !
       call part_gather(surfaces(i)%crds, coords,surfaces(i)%conn3d, surfaces(i)%tr)
       !
       !  Now compute normal vectors at quadrature points.
       !
       p_area = 0.0_RK
       !
       do iel=surfaces(i)%semin, surfaces(i)%semax
          tangent = 0.0_RK
          do jq=1, nqp2d
             do jn=1, SFTYPE
                jdof = 3*(jn-1)
                !
                !  First tangent vector.
                !  
                tangent(1, 1, jq) = tangent(1, 1, jq) +&
                     & surfaces(i)%crds(jdof,iel)  *sfgqp2d(1,jn, jq)
                tangent(2, 1, jq) = tangent(2, 1, jq) +&
                     & surfaces(i)%crds(jdof+1,iel)*sfgqp2d(1,jn, jq)
                tangent(3, 1, jq) = tangent(3, 1, jq) +&
                     & surfaces(i)%crds(jdof+2,iel)*sfgqp2d(1,jn, jq)
                !
                !  Second tangent vector.
                !  
                tangent(1, 2, jq) = tangent(1, 2, jq) +&
                     & surfaces(i)%crds(jdof,iel)  *sfgqp2d(2,jn, jq)
                tangent(2, 2, jq) = tangent(2, 2, jq) +&
                     & surfaces(i)%crds(jdof+1,iel)*sfgqp2d(2,jn, jq)
                tangent(3, 2, jq) = tangent(3, 2, jq) +&
                     & surfaces(i)%crds(jdof+2,iel)*sfgqp2d(2,jn, jq)
             enddo
             !
             !  Now compute normals.
             !
             normal(1, jq) =&
                  & tangent(2,1,jq)*tangent(3,2,jq) -&
                  & tangent(3,1,jq)*tangent(2,2,jq)
             normal(2, jq) =&
                  & tangent(3,1,jq)*tangent(1,2,jq) -&
                  & tangent(1,1,jq)*tangent(3,2,jq)
             normal(3, jq) =&
                  & tangent(1,1,jq)*tangent(2,2,jq) -&
                  & tangent(2,1,jq)*tangent(1,2,jq)
             nmag = dsqrt(&
                  & normal(1, jq)*normal(1, jq) +&
                  & normal(2, jq)*normal(2, jq) +&
                  & normal(3, jq)*normal(3, jq)  )
             !     
             if (nmag > 0) then
                normal(:,jq) = normal(:,jq)/nmag
                sjac(jq) = nmag
             else
                call par_quit('Error  :     > Surface normal zero magnitude.')
             endif
             !
          enddo
          do jq=1, nqp2d
             p_area = p_area + wt2d(jq)*sjac(jq)
          enddo
       enddo !iel
       ! 
       call par_sum(p_area, area)
       !
       area0(i)=area
       !
    enddo !nsurfaces
    !
    RETURN
    !
  END SUBROUTINE compute_area
!    
!
END MODULE surf_info_mod
!
