! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE HardeningModule
  !
  ! Hardening laws.
  !
  ! ==================== Other Modules (USE statements)
  !
  USE LibF95, ONLY: RK => REAL_KIND, RK_ONE, RK_ZERO
  USE parallel_mod, ONLY: Quit => par_quit
  USE DimsModule
  USE READ_INPUT_MOD
  USE microstructure_mod
  USE MATRIX_OPERATIONS_MOD

  IMPLICIT NONE
  !
  !  ==================== Public Entities
  !
  !  variables, procedures, constants, derived types and namelist groups
  !
  PRIVATE   ! all objects are private unless declared otherwise
  !
  PUBLIC :: HARD_LAW_VOCE, EVAL_RATE, EVAL_RATE_DER
  PUBLIC :: hard_law
  !
  !  ==================== Module Data
  !
  INTEGER, PARAMETER :: HARD_LAW_VOCE=2
  INTEGER, PARAMETER :: EVAL_RATE=1, EVAL_RATE_DER=2
  INTEGER            :: n_slip
  !
CONTAINS ! ============================================= MODULE PROCEDURES
  !
  SUBROUTINE hard_law(func, dfunc, crss, crss_sat, shear, shrate, epseff, icode, ihard, n, m)
    !
    ! Evaluate hardening rate or derivative of hardening rate
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: icode, ihard, n, m
    REAL(RK), INTENT(OUT) :: func(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
    REAL(RK), INTENT(OUT) :: dfunc(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
    REAL(RK), INTENT(IN) :: crss(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
    REAL(RK), INTENT(IN) :: crss_sat(0:(n - 1), 0:(m - 1))
    REAL(RK), INTENT(IN) :: shrate(0:(n - 1), 0:(m - 1))
    REAL(RK), INTENT(IN) :: shear(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
    REAL(RK), INTENT(IN) :: epseff(0:(n - 1), 0:(m - 1))
    !
    ! *icode* -  the hardening law
    ! *ihard* - flag indicating return value of rate or derivative of rate
    ! *n*, *m* - array dimensions
    ! *crss* - hardness array
    ! *crss_sat* - saturation hardness array
    ! *shrate* - shear rates
    ! *func* - array of computed hardening rates
    ! *dfunc* - array of computed hardening rate derivatives
    ! 
    !  ========== Locals
    ! 
    REAL(RK) :: mysign(0:MAXSLIP1, 0:(n - 1), 0:(m - 1)), ratio(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
    INTEGER :: numind, n_slip
    INTEGER, POINTER :: indices(:) => NULL()
    INTEGER :: iphase, islip, my_phase(0:(m - 1))
    REAL(RK) :: c3, c2, c10
    
    !
    !----------------------------------------------------------------------
    !

    IF (ihard /= HARD_LAW_VOCE) THEN
      !write(*,*) 'ihard: ', ihard
      CALL PAR_QUIT('Error  :     > Voce hardening must be used.')
    END IF
    !RC
    !Initialize these to zero or else the terms that are not used in a dual phase material get assigned random values.
    func = RK_ZERO
    dfunc = RK_ZERO
    !RC
    !In future updates this will need to be set up such that it makes a call for either a isotropic case or latent hardening case
    my_phase(:) = phase(el_sub1:el_sup1)

    mysign = RK_ONE

    do islip=0,MAXSLIP1

        WHERE( (crss_sat - crss(islip,:,:)) .LE. RK_ZERO) mysign(islip,:,:) = RK_ZERO

    enddo
    !
    !Options for which hardening model to use
    IF(options%hard_type.EQ.'isotropic') THEN
        call isotropic_hardening(func, dfunc, crss, crss_sat, shrate, icode, ihard, n, m, mysign, my_phase)
    ELSEIF(options%hard_type .EQ. 'cyclic_isotropic') THEN
        call cyclic_hardening(func, dfunc, crss, crss_sat, shear, shrate, epseff, icode, ihard, n, m, mysign, my_phase)
    ELSEIF(options%hard_type .EQ. 'cyclic_anisotropic') THEN
        CALL PAR_QUIT('Error  :     > Anisotropic cyclic hardening is not currently supported.')   
    ELSEIF(options%hard_type .EQ. 'latent') THEN
        call anisotropic_hardening(func, dfunc, crss, crss_sat,icode, ihard, n, m, mysign, my_phase)
    ELSE
        CALL PAR_QUIT('Error  :     > Invalid hardening model option input.')   
    ENDIF

    END SUBROUTINE hard_law

    SUBROUTINE isotropic_hardening(func, dfunc, crss, crss_sat, shrate, icode, ihard, n, m, mysign, my_phase)


        IMPLICIT NONE

        INTEGER, INTENT(IN) :: icode, ihard, n, m, my_phase(0:(m - 1))
        REAL(RK), INTENT(OUT) :: func(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(OUT) :: dfunc(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(IN) :: crss(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(IN) :: crss_sat(0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(IN) :: shrate(0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(IN) :: mysign(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))

        REAL(RK) :: ratio(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
        INTEGER :: numind, n_slip
        INTEGER, POINTER :: indices(:) => NULL()
        INTEGER :: iphase, islip
        REAL(RK) :: c3, c2, c10

        DO iphase=1,numphases !start of loop

          c2 = crystal_parm(2,iphase)
          c3 = crystal_parm(3,iphase)

          call CrystalTypeGet(ctype(iphase))
          n_slip=ctype(iphase)%numslip
          
          CALL find_indices(numind, iphase, my_phase, indices)

            !RC
            !All of the slip systems are just looped over instead of just the one.
            !Since only the isotropic case is being examined here it should be possible to do only
            !   one computation of the slip system and then assign all the other slip systems that value
            !Future update should look into that implementation
          do islip=0,n_slip-1
             if (islip==0) then
                ratio(islip,:,indices) = (crss_sat(:,indices) - crss(islip,:,indices)) / &
                     &             (crss_sat(:,indices) - c3)

                SELECT CASE(icode)
                  !                                                                                   
                CASE (EVAL_RATE)
                  !                                                                                   
                  func(islip,:,indices) = shrate(:,indices) * mysign(islip,:,indices) * c2 * &
                       &  ratio(islip,:,indices)**n_voce(iphase)

                CASE (EVAL_RATE_DER)
                  !                                                                                   
                  dfunc(islip,:,indices) = - c2 /(crss_sat(:,indices) - crystal_parm(3,iphase))&
                       & * shrate(:,indices) * mysign(islip,:,indices) * n_voce(iphase)
                CASE DEFAULT
                  !                                                                                   
                END SELECT
             else
                select case(icode)

                case (EVAL_RATE)
                   func(islip,:,indices)=func(0,:,indices)

                case (EVAL_RATE_DER)
                   dfunc(islip,:,indices)=dfunc(0,:,indices)
                case DEFAULT

                end select
             endif

           enddo !n_slip
        
          !
          DEALLOCATE(indices)
          !
        ENDDO !numphases

    END SUBROUTINE isotropic_hardening


    SUBROUTINE cyclic_hardening(func, dfunc, crss, crss_sat, shear, shrate, epseff, icode, ihard, n, m, mysign, my_phase)


        IMPLICIT NONE

        INTEGER, INTENT(IN) :: icode, ihard, n, m, my_phase(0:(m - 1))
        REAL(RK), INTENT(OUT) :: func(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(OUT) :: dfunc(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(IN) :: crss(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(IN) :: crss_sat(0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(IN) :: shrate(0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(IN) :: mysign(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(IN) :: shear(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(IN) :: epseff(0:(n - 1), 0:(m - 1))

        REAL(RK) :: ratio(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
        REAL(RK) :: ratio_sh(0:(n - 1), 0:(m - 1))
        REAL(RK) :: active_shrate(0:(n - 1), 0:(m - 1))
        REAL(RK) :: shear_crit(0:(n - 1), 0:(m - 1))
        INTEGER  :: numind, n_slip, n_gt1, n_lt0
        INTEGER, POINTER :: indices(:) => NULL()
        INTEGER  :: iphase, islip
        REAL(RK) :: c3, c2, c10
        REAL(RK) :: shr_max(0:(n - 1), 0:(m - 1)), shr_min(0:(n - 1), 0:(m - 1))
        REAL(RK), PARAMETER :: mach_eps = 2.22e-16

        active_shrate(:,:) = 0.0_RK
        n_lt0 = 0

        DO iphase = 1, numphases !start of loop

          c2 = crystal_parm(2,iphase)
          c3 = crystal_parm(3,iphase)

          call CrystalTypeGet(ctype(iphase))
          n_slip=ctype(iphase)%numslip
          
          CALL find_indices(numind, iphase, my_phase, indices)

            !RC
            !All of the slip systems are just looped over instead of just the one.
            !Since only the isotropic case is being examined here it should be possible to do only
            !   one computation of the slip system and then assign all the other slip systems that value
            !Future update should look into that implementation

          ratio_sh(:, indices) = dabs(crss(0, :, indices)/crss_sat(:, indices))

          shear_crit(:,indices) = cyclic_parm(0, iphase) * &
          &     (ratio_sh(:, indices)**cyclic_parm(1, iphase))

          where(shear_crit(:, indices) .LT. mach_eps)

            shear_crit(:, indices) = 0.0_RK
            !n_lt0 = n_lt0 + 1

          endwhere

          if (any(dabs(shear_crit) .GE. 1.0)) then
            n_gt1 = count((dabs(shear_crit) .GE. 1.0))
            !print *, 'Number greater than 1: ', n_gt1
            call par_quit('Error  :     > `SHEAR_CRIT` is greater than 1.')
          endif

          DO islip = 0, n_slip-1

            where (accumshear(islip,0,indices+el_sub1) .GE. shear_crit(0,indices))
                active_shrate(0,indices) = active_shrate(0,indices) + &
                    &     dabs(shear(islip, 0, indices))
            endwhere

          ENDDO

          shr_min = 1.0d-6 * epseff
          shr_max = 1.0d1 * epseff

          where (active_shrate .LE. shr_min)   active_shrate = shr_min
          where (active_shrate .GE. shr_max)   active_shrate = shr_max

          DO islip = 0, n_slip-1
             if (islip==0) then

                ratio(islip,:,indices) = (crss_sat(:,indices) - crss(islip,:,indices)) / &
                     &             (crss_sat(:,indices) - c3)

                SELECT CASE(icode)
                  !                                                                                   
                CASE (EVAL_RATE)
                  !                                                                                   
                  func(islip,:,indices) = active_shrate(:,indices) * mysign(islip,:,indices) * c2 * &
                       &  ratio(islip,:,indices)**n_voce(iphase)

                CASE (EVAL_RATE_DER)
                  !                                                                                   
                  dfunc(islip,:,indices) = - c2 /(crss_sat(:,indices) - crystal_parm(3,iphase))&
                       & * active_shrate(:,indices) * mysign(islip,:,indices) * n_voce(iphase)
                CASE DEFAULT
                  !                                                                                   
                END SELECT
             else
                select case(icode)

                case (EVAL_RATE)
                   func(islip,:,indices)=func(0,:,indices)

                case (EVAL_RATE_DER)
                   dfunc(islip,:,indices)=dfunc(0,:,indices)
                case DEFAULT

                end select
             endif

           ENDDO !n_slip
        
          !
          DEALLOCATE(indices)
          !
        ENDDO !numphases

    END SUBROUTINE cyclic_hardening


    SUBROUTINE anisotropic_hardening (func, dfunc, crss, crss_sat,icode, ihard, n, m, mysign, my_phase)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: icode, ihard, n, m, my_phase(0:(m - 1))
        REAL(RK), INTENT(OUT) :: func(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(OUT) :: dfunc(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(IN) :: crss(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(IN) :: crss_sat(0:(n - 1), 0:(m - 1))
        REAL(RK), INTENT(IN) :: mysign(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))

        REAL(RK) :: ratio(0:MAXSLIP1, 0:(n - 1), 0:(m - 1)), shrate(0:MAXSLIP1)
        INTEGER :: numind, nslip, xtype
        INTEGER, ALLOCATABLE :: indices(:)
        INTEGER :: iphase, ielms,dim, index
        REAL(RK) :: c3, c2, c10
        LOGICAL :: logVec(0:(m-1))

        shrate=RK_ZERO

        DO iphase=1,numphases !start of loop

          c2 = crystal_parm(2,iphase)
          c3 = crystal_parm(3,iphase)

          call CrystalTypeGet(ctype(iphase))
          nslip=ctype(iphase)%numslip-1
          if (options%hard_type.eq."test")then
            xtype=4
            !Set it outside the available crystal types to trigger the default calculation of the isotropic case
          else
            xtype=ctype(iphase)%class
          endif
          logVec = (my_phase .eq. iphase)
!          CALL indexFinder(logVec, indices)
          dim=size(logVec)
!          write(*,*) logVec

           do ielms=0,dim-1
              if(logVec(ielms)) then !Find where the element is the same phase iphase
                  ratio(0:nslip,0,ielms) = (crss_sat(0,ielms) - crss(0:nslip,0,ielms)) / &
                       &             (crss_sat(0,ielms) - c3)
                  call calculate_shrate(shrate(:),xtype,ielms)
                                  
                  SELECT CASE(icode)
                    !
                  CASE (EVAL_RATE)! all the slip system func and dfunc values for the element are calculated in one go
                    !
                    func(0:nslip,0,ielms) = shrate(0:nslip) * mysign(0:nslip,0,ielms) * c2 * &
                         &  ratio(0:nslip,0,ielms)**n_voce(iphase)

                  CASE (EVAL_RATE_DER)
                    !
                    dfunc(0:nslip,0,ielms) = - c2 /(crss_sat(0,ielms) - crystal_parm(3,iphase))&
                         & * shrate(0:nslip) * mysign(0:nslip,0,ielms) * n_voce(iphase)
                  CASE DEFAULT 
                    !
                  END SELECT
                ENDIF

           enddo !n_slip
        
          !
!          DEALLOCATE(indices)
          !
        ENDDO !numphases

    END SUBROUTINE anisotropic_hardening

    SUBROUTINE matrix_vec_mult(m,v,mv,dim)
        !A simple matrix vector multiplication routine that takes advantage of the compiler's vectorization optimizations
        !It ends up being faster than using the built in Fortran matmul method
        IMPLICIT NONE
        integer, intent(in) :: dim
        REAL(RK), intent(in) :: m(0:dim-1,0:dim-1)
        REAL(RK), intent(in) :: v(0:dim-1)
        REAL(RK), intent(out) :: mv(0:dim-1)
        integer :: i
        mv=0.0_RK
        do i=0,dim-1
            mv(:)=mv(:)+m(:,i)*v(i)
        enddo
!        call v_print(v)

    END SUBROUTINE matrix_vec_mult

    SUBROUTINE calculate_shrate(shrate,xtype,ind)
        !shrate calculations are based upon the assumption that slip systems that share the same slip plane interact instead of the same slip direction

        IMPLICIT NONE
        REAL(RK), INTENT(out) :: shrate(0:MAXSLIP1)
        INTEGER, INTENT(in) :: xtype, ind

        REAL(RK) :: tmp(0:1), tmpO(0:1), ones(0:MAXSLIP1,0:MAXSLIP1)

        ones=RK_ONE

            SELECT CASE(xtype) !Values are hard wired

            CASE (1) !FCC xtype

                call matrix_vec_mult(fcc_h1,abs(gammadot(0:2,0,ind+el_sub1)),shrate(0:2),3)
                call matrix_vec_mult(fcc_h2,abs(gammadot(3:5,0,ind+el_sub1)),shrate(3:5),3)
                call matrix_vec_mult(fcc_h3,abs(gammadot(6:8,0,ind+el_sub1)),shrate(6:8),3)
                call matrix_vec_mult(fcc_h4,abs(gammadot(9:11,0,ind+el_sub1)),shrate(9:11),3)

            CASE (3) !HCP xtype

                call matrix_vec_mult(hcp_h1,abs(gammadot(0:2,0,ind+el_sub1)),shrate(0:2),3)
                shrate(3:5)=hcp_vert*abs(gammadot(3:5,0,ind+el_sub1))
                call matrix_vec_mult(hcp_h2,abs(gammadot(6:7,0,ind+el_sub1)),shrate(6:7),2)
                call matrix_vec_mult(hcp_h3,abs(gammadot(8:9,0,ind+el_sub1)),shrate(8:9),2)
                call matrix_vec_mult(hcp_h4,abs(gammadot(10:11,0,ind+el_sub1)),shrate(10:11),2)
                call matrix_vec_mult(hcp_h5,abs(gammadot(12:13,0,ind+el_sub1)),shrate(12:13),2)
                call matrix_vec_mult(hcp_h6,abs(gammadot(14:15,0,ind+el_sub1)),shrate(14:15),2)
                call matrix_vec_mult(hcp_h7,abs(gammadot(16:17,0,ind+el_sub1)),shrate(16:17),2)

            CASE (2) !BCC xtype

                tmp=abs(gammadot((/ 0,9/),0,ind+el_sub1))
                call matrix_vec_mult(bcc_h1,tmp,tmpO,2)
                shrate((/ 0,9/))=tmpO

                tmp=abs(gammadot((/ 1,7/),0,ind+el_sub1))
                call matrix_vec_mult(bcc_h2,tmp,tmpO,2)
                shrate((/ 1,7/))=tmpO

                tmp=abs(gammadot((/ 2,5/),0,ind+el_sub1))
                call matrix_vec_mult(bcc_h3,tmp,tmpO,2)
                shrate((/ 2,5/))=tmpO

                tmp=abs(gammadot((/ 3,6/),0,ind+el_sub1))
                call matrix_vec_mult(bcc_h4,tmp,tmpO,2)
                shrate((/ 3,6/))=tmpO

                tmp=abs(gammadot((/ 4,10/),0,ind+el_sub1))
                call matrix_vec_mult(bcc_h5,tmp,tmpO,2)
                shrate((/ 4,10/))=tmpO

                tmp=abs(gammadot((/ 8,11/),0,ind+el_sub1))
                call matrix_vec_mult(bcc_h6,tmp,tmpO,2)
                shrate((/ 8,11/))=tmpO

            CASE DEFAULT

                call matrix_vec_mult(ones,abs(gammadot(:,0,ind+el_sub1)),shrate,MAXSLIP1+1)
!                print *, "Assuming Isotropic Case since crystal type is not a default"
!                CALL Quit('no such crystal type implemented: exiting', ABORT=.TRUE.)
!                call par_quit('Invalid crystal type was choosen for the hardening model')

            END SELECT

    END SUBROUTINE


  !
END MODULE HardeningModule
