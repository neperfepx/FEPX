! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE RestartModule
  !
  !  *** Program Unit:  Module
  !  ***    Unit Name:  RestartModule
  !  ***  Description:
  !
  !  This module allows for restart capability.
  !
  !  *** Use Statements:
  !
  USE IntrinsicTypesModule, &
       &  RK=>REAL_KIND, IK=>INTEGER_KIND, LK=>LOGICAL_KIND
  USE parallel_mod
  USE units_mod
  USE READ_INPUT_MOD
  USE post_update_n_mod
  USE stress_strain_mod
  USE SIMULATION_CONFIGURATION_MOD
  USE DimsModule
  !
  !  *** End:
  !
  IMPLICIT NONE
  !
  !----------------------------------------------------------------------
  !
  CONTAINS
  !
  SUBROUTINE RestartWriteField(velocity, c0_angs, c_angs, &
            & rstar, rstar_n, wts, crss, crss_n, &
            & e_elas_kk_bar, sig_vec_n, EqStrain, EqPlStrain, gamma)
    !
    !  *** Program Unit:  subroutine
    !  ***    Unit Name:  RestartWriteField
    !
    !  *** Unit Declaration: 
    !
    !  ***  Description:
    !
    !  Writes field data for restarting simulation
    !
    !  *** Argument Declarations:
    !
    REAL(RK), INTENT(IN) :: velocity(dof_sub1:dof_sup1)
    REAL(RK), INTENT(IN) :: c0_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(IN) :: c_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(IN) :: rstar(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(IN) :: rstar_n(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(IN) :: wts(0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(IN) :: crss(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(IN) :: crss_n(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(IN) :: e_elas_kk_bar(el_sub1:el_sup1)
    REAL(RK), INTENT(IN) :: sig_vec_n(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(IN) :: EqStrain(el_sub1:el_sup1)
    REAL(RK), INTENT(IN) :: EqPlStrain(el_sub1:el_sup1)
    REAL(RK), INTENT(IN) :: gamma(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
    !
    !  *** End:
    !
    INTEGER :: myunit
    INTEGER :: IOSTATUS
    !
    INTRINSIC :: TRIM
    !
    !----------------------------------------------------------------------
    !
    myunit = NewUnitNumber()
    open(UNIT=myunit, FILE=TRIM(options%rsfield_base_out)//'.'//TRIM(myidstr), &
         &   FORM='UNFORMATTED', ACTION='WRITE', IOSTAT = IOSTATUS)
    !
    IF (IOSTATUS .NE. 0) THEN
        !
        CALL PAR_QUIT("Error  :     > Failure to write restart field data.")
        !
    END IF
    !
    !  Velocity and coordinates.
    !
    write(myunit) coords
    write(myunit) velocity  
    !
    !  Orientations, weights and hardnesses.
    !
    write(myunit) c0_angs
    write(myunit) c_angs
    write(myunit) rstar
    write(myunit) rstar_n
    write(myunit) wts
    write(myunit) crss
    write(myunit) crss_n
    !
    !  Elastic Strains.
    !
    write(myunit) gela_kk_bar
    write(myunit) gsig_vec_n
    write(myunit) pela_kk_bar
    write(myunit) psig_vec_n
    write(myunit) e_elas_kk_bar
    write(myunit) sig_vec_n
    !
    !  Equivalent Strains
    !
    write(myunit) EqStrain
    write(myunit) EqPlStrain
    write(myunit) gamma
    !
    close(myunit)
    !
  END SUBROUTINE RestartWriteField
  !
  !
  ! *********************************************************************
  !
  !
  SUBROUTINE RestartReadField(velocity, c0_angs, c_angs, &
            & rstar, rstar_n, wts, crss, crss_n, &
            & e_elas_kk_bar, sig_vec_n, EqStrain, EqPlStrain, gamma)
    !
    !  *** Program Unit:  subroutine
    !  ***    Unit Name:  RestartRead
    !
    !  *** Unit Declaration: 
    !
    !  ***  Description:  
    !
    !  Reads field data for restarting simulation
    !
    !  *** Argument Declarations:
    !
    REAL(RK), INTENT(OUT) :: velocity(dof_sub1:dof_sup1)
    REAL(RK), INTENT(OUT) :: c0_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(OUT) :: c_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(OUT) :: rstar(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(OUT) :: rstar_n(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(OUT) :: wts(0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(OUT) :: crss(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(OUT) :: crss_n(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(OUT) :: e_elas_kk_bar(el_sub1:el_sup1)
    REAL(RK), INTENT(OUT) :: sig_vec_n(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
    REAL(RK), INTENT(OUT) :: EqStrain(el_sub1:el_sup1)
    REAL(RK), INTENT(OUT) :: EqPlStrain(el_sub1:el_sup1)
    REAL(RK), INTENT(OUT) :: gamma(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
    !
    !  *** Locals:
    !
    INTEGER :: myunit
    INTEGER :: IOSTATUS
    !
    INTRINSIC :: TRIM
    !
    !----------------------------------------------------------------------
    !
    myunit = NewUnitNumber()
    open(UNIT=myunit, FILE=TRIM(options%rsfield_base_in)//'.'//TRIM(myidstr), &
         &   FORM='UNFORMATTED', ACTION='READ', IOSTAT = IOSTATUS)
    !
    IF (IOSTATUS .NE. 0) THEN
        !
        CALL PAR_QUIT("Error  :     > Failure to read restart field data.")
        !
    END IF
    !
    !  Velocity and coordinates.
    !
    read(myunit) coords
    read(myunit) velocity  
    !
    !  Orientations, weights and hardnesses.
    !
    read(myunit) c0_angs
    read(myunit) c_angs
    read(myunit) rstar
    read(myunit) rstar_n
    read(myunit) wts
    read(myunit) crss
    read(myunit) crss_n
    !
    !  Elastic Strains.
    !
    read(myunit) gela_kk_bar
    read(myunit) gsig_vec_n
    read(myunit) pela_kk_bar
    read(myunit) psig_vec_n
    read(myunit) e_elas_kk_bar
    read(myunit) sig_vec_n
    !
    !  Equivalent Strains
    !
    read(myunit) EqStrain
    read(myunit) EqPlStrain
    read(myunit) gamma
    !
    close(myunit)
    !
  END SUBROUTINE RestartReadField
  !
END MODULE RestartModule
