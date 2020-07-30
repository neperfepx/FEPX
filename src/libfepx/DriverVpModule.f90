! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE DriverVpModule
  !
  !  Driver for viscoplastic problem
  !
  !  ==================== Other Modules (USE statements)
  !
  USE LibF95, ONLY: RK => REAL_KIND

  USE gather_scatter
  USE parallel_mod
  USE units_mod
  USE READ_INPUT_MOD
  USE DimsModule
  use ItMethodVpModule
  USE WRITE_OUTPUT_MOD, ONLY: PRINT_STEP

  IMPLICIT NONE
  !
  !  ==================== Public Entities
  !
  !  variables, procedures, constants, derived types and namelist groups
  !
  PRIVATE   ! all objects are private unless declared otherwise
  PUBLIC :: driver_vp_solve
  !
CONTAINS ! ============================================= MODULE PROCEDURES
  !
  !
  !
  !**********************************************************************
  !
  SUBROUTINE driver_vp_solve(&
       &   itype, bcs, velocity, pforce,&
       &   dof_trace,&
       &   np_trace, crss,&
       &   c0_angs, rstar, wts)
    !
    !     Driver for viscoplastic solution.
    !
    !----------------------------------------------------------------------
    !
    !     Modules:
    !
    !     
    !     Included files:
    !
    IMPLICIT  NONE
    !
    !
    !     Arguments:
    !
    TYPE(trace) dof_trace, np_trace
    !
    LOGICAL :: bcs(dof_sub1:dof_sup1)
    !
    !     itype: flag, 0=isotropic viscoplastic,1=anisotropic viscoplastic
    !
    INTEGER :: itype
    !
    REAL(RK) :: velocity(dof_sub1:dof_sup1)
    REAL(RK) :: pforce  (dof_sub1:dof_sup1)
    REAL(RK) :: c0_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK) :: rstar  (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK) :: crss(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK) :: wts (0:ngrain1, el_sub1:el_sup1)
    !
    !     Locals:
    !
    INTEGER :: incr, ier, m_el
    !
    REAL(RK) :: dtime
    REAL(RK) :: elpress(el_sub1:el_sup1)
    REAL(RK) :: evel(0:kdim1, el_sub1:el_sup1)
    REAL(RK) :: c_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK) :: qr5x5 (0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK) :: sig_vec        (0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK) :: epseff(el_sub1:el_sup1)
    ! tm268_M13:
    !      REAL(RK) :: eqplas(el_sub1:el_sup1)

    !
    !     sig_dummy is currently just for passing to print_incr, but
    !     in the future one could store and print the grain stresses.
    !
    REAL(RK) :: sig_avg (0:TVEC1, el_sub1:el_sup1)
    REAL(RK) :: sig_avg_kk       (el_sub1:el_sup1)
    REAL(RK) :: sig_dummy(0:TVEC, el_sub1:el_sup1)
    REAL(RK) :: sdummy3x3(0:DIMS1,0:DIMS1,el_sub1:el_sup1)
    !
    !----------------------------------------------------------------------
    ! 
    m_el = el_sup1 - el_sub1 + 1
    !
    GO TO (100, 200) itype + 1
    !
100 CONTINUE
    !----------------------------------------------------------------------
    !     Isotropic Solution : Viscoplastic
    !----------------------------------------------------------------------
    if (myid .eq. 0)  write(DFLT_U,'(A)') 'Info   :   - Initializing fields &
        &from isotropic viscoplastic solution'

    elpress = 0.0_RK
    evel    = 0.0_RK
    dtime   = 0.0_RK
    incr    = 0
    ! tm268_M13:
    !      eqplas  = 0.0_RK
    !
    !
    ! tm268_M13: removed eqplas
    ier =  itmethod_vp(itype, bcs, pforce, velocity,&
         &          elpress, evel,&
         &          dof_trace, np_trace,&
         &          qr5x5, wts, epseff,&
         &          dtime, incr)
    !      
    if (ier .lt. 0) then
      call par_quit('Error  :     > ITMETHOD_VP failed to converge.')
    endif
    !
    !Output computed quantities
    !
    call part_gather(element_crds, coords, nodes, dof_trace)
    !
    call print_step(incr, itype, coords, velocity, c0_angs, wts, crss)
    ! 
    !  
    RETURN
    !
    ! 
200 CONTINUE
    !----------------------------------------------------------------------
    !     Anisotropic Solution : Viscoplastic
    !----------------------------------------------------------------------
    call par_quit('Error  :     > ANISOTROPIC_VP is no longer implemented.')
    !
    !
    RETURN
  END SUBROUTINE driver_vp_solve
  !      
END MODULE DriverVpModule
