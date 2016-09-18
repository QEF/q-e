!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
MODULE charg_resp
  !----------------------------------------------------------------------------
  ! This module contains charge response calculation related variables & 
  ! subroutines.
  ! 
  ! Created by Osman Baris Malcioglu 2009
  ! Cleaned up and refactored Simon Binnie, 2013
  !----------------------------------------------------------------------------

  USE kinds,                ONLY : dp
  USE lr_variables,         ONLY : lr_verbosity
  USE io_global,            ONLY : ionode, stdout,ionode_id


  REAL(kind=dp), ALLOCATABLE    ::  rho_1_tot(:,:)  !response charge density
                                                    !(mode 2)
  COMPLEX(kind=dp), ALLOCATABLE ::  rho_1_tot_im(:,:)!response charge density,
                                  ! imaginary part used in resonance condition

  CHARACTER(len=256)         :: w_T_prefix ! prefix for storage of previous calculation
  INTEGER                    :: w_T_npol   ! number of polarization directions
                                           ! considered in previous run
  REAL(kind=dp), ALLOCATABLE :: &  ! the required parts of the lanczos matrix 
                                w_T_beta_store(:),&         ! for w_T (iter) 
                                w_T_gamma_store(:)

  COMPLEX(kind=dp), ALLOCATABLE :: w_T (:) ! The solution to (omega-T) (iter)
  REAL(kind=dp)                 :: omeg    ! frequencies for calculating charge
                                           ! response
  REAL(kind=dp)                 :: epsil   ! Broadening
  REAL(kind=dp)                 :: w_T_norm0_store       ! The norm for this step
  COMPLEX(kind=dp), ALLOCATABLE :: w_T_zeta_store(:,:)   ! The zeta coefficients
                                                         ! from file
  COMPLEX(kind=dp), ALLOCATABLE :: chi(:,:)              ! The susceptibility
                                                         ! tensor for the given
                                                         ! frequency
  LOGICAL :: resonance_condition

CONTAINS
!-----------------------------------------------------------------------
  SUBROUTINE lr_calc_project(ip)
    !-----------------------------------------------------------------------
    ! Handles the calculation and output of the results of projection analysis
    ! for the polarisation ip, or if ip=4, for the entire strength function.
    !
    USE kinds,         ONLY : DP
    USE lr_variables,  ONLY : F, R, nbnd_total, LR_polarization, n_ipol
    USE wvfct,         ONLY : nbnd

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ip
    !
    COMPLEX(kind=DP) :: sum_f, sum_c
    INTEGER          :: ibnd_occ, ibnd_virt, ic

    IF(ip==4) THEN
       !
       WRITE(stdout,'(/,/5x,"Participation of virtual states to absorbtion &
            &coefficent")')
       WRITE(stdout,'(5x,"occ",1x,"vir",5x,"Re(Tr(F.R))",6x,"Im(TR(F.R))",5x,&
            &"fraction in alpha")')
       !
       DO ibnd_occ = 1, nbnd
          !
          DO ibnd_virt = 1, (nbnd_total-nbnd)
             !
             sum_f = CMPLX(0.0d0,0.0d0,dp)
             sum_c = CMPLX(0.0d0,0.0d0,dp)
             !
             DO ic = 1, n_ipol
                !
                sum_f = sum_f + F(ibnd_occ, ibnd_virt, ic)&
                     &* CONJG( R(ibnd_occ, ibnd_virt, ic) )
                sum_c = sum_c + chi(ic, ic)
                !
             ENDDO
             !
             WRITE(stdout,'(5x,i3,1x,i3,3x,E16.8,2X,E16.8,2X,F8.5)') &
                  ibnd_occ, ibnd_virt, DBLE( sum_F ),AIMAG( sum_F ),&
                  &( AIMAG( sum_f ) / AIMAG( sum_c ) )
          ENDDO
       ENDDO
       !
    ELSE
       !
       WRITE(stdout,'(/,/5x,"Projection of virtual states for polarization &
            &direction",1x,i8)') LR_polarization
       WRITE(stdout,'(2x,"occ",1x,"vir",8x,"Re(F)",14x,"Im(F)",8x, " Frac. &
            &pres. in Re(chi_",I1,"_",I1,") and Im(chi_",I1,"_",I1,")")') &
            &  ip,ip,ip,ip
       !
       sum_f=CMPLX(0.0d0,0.0d0,dp)
       !
       DO ibnd_occ=1,nbnd
          !
          DO ibnd_virt=1,(nbnd_total-nbnd)
             !
             F(ibnd_occ, ibnd_virt, ip) = F(ibnd_occ, ibnd_virt, ip)&
                                         &* CMPLX(w_T_norm0_store, 0.0d0, dp)
             !
             sum_f = F(ibnd_occ, ibnd_virt, ip) &
                    &* CONJG( R(ibnd_occ, ibnd_virt, ip) )
             !
             WRITE(stdout,'(2x,i3,1x,i3,3x,E16.8,2X,E16.8,17X,F8.5,2x,F8.5)') &
                  ibnd_occ, ibnd_virt, DBLE( F(ibnd_occ, ibnd_virt, ip) ),&
                  & AIMAG( F(ibnd_occ, ibnd_virt, ip) ),&
                  & ( DBLE(sum_f) / DBLE( chi(ip ,ip) )),&
                  & ( AIMAG(sum_f)/ AIMAG( chi(ip, ip) ))
          ENDDO
          !
       ENDDO
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE lr_calc_project
  !---------------------------------------------------------------------
  SUBROUTINE lr_project_init()
    !-----------------------------------------------------------------------
    ! Handles the initial calculation of the oscilator strengths for 
    ! projection analysis and writes them to stdout.
    !
    USE kinds,         ONLY : DP
    USE lr_variables,  ONLY : R, nbnd_total, n_ipol, project
    USE wvfct,         ONLY : nbnd
    USE control_lr,    ONLY : nbnd_occ

    IMPLICIT NONE

    INTEGER :: ip, ibnd_occ, ibnd_virt

    !
    IF(nbnd_total > nbnd_occ(1)) THEN
       WRITE(stdout,'(/,5X,"Virtual states in ground state run will be used in&
            & projection analysis")')
    ELSE
       !
       WRITE(stdout,'(/,5X,"No virtual states for projection found")')
       project=.FALSE.
       !
       RETURN
       !
    ENDIF
    !
    CALL sd0psi() 
    !
    CALL lr_calc_R()
    !
    DO ip=1, n_ipol
       !
       WRITE(stdout,'(/,/5x,"Oscillator strengths for polarization &
            &direction",1x,i8)') ip
       WRITE(stdout,'(5x,"occ",1x,"con",8x,"Re(R)",14x,"Im(R)")')
       !
       DO ibnd_occ = 1, nbnd
          !
          DO ibnd_virt = 1, (nbnd_total-nbnd)
             !
             WRITE(stdout,'(5x,i3,1x,i3,3x,E16.8,2X,E16.8)') &
                  ibnd_occ, ibnd_virt, DBLE( R(ibnd_occ, ibnd_virt, ip) ),&
                  &AIMAG( R(ibnd_occ, ibnd_virt, ip) )
          ENDDO
          !
       ENDDO
       !
    ENDDO
    !
    RETURN
    !
  END SUBROUTINE lr_project_init
  !---------------------------------------------------------------------
  SUBROUTINE read_wT_beta_gamma_z()
    !-----------------------------------------------------------------------
    ! Reads beta_gamma_z from a previous calculation for the given polarization
    ! direction pol.
    !---------------------------------------------------------------------
    !
    USE mp,                   ONLY : mp_bcast, mp_barrier
    USE lr_variables,         ONLY : LR_polarization, itermax
    USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
    USE mp_world,             ONLY : world_comm
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), EXTERNAL :: trimcheck
    !
    CHARACTER(len=6), EXTERNAL :: int_to_char
    ! local
    LOGICAL :: exst
    INTEGER :: iter_restart,i,j
    CHARACTER(len=256) :: filename
    REAL(KIND=DP) :: discard
    

    CALL start_clock( 'post-processing' )
    IF (lr_verbosity > 5) WRITE(stdout,'("<read_wT_beta_gamma_z>")')
#if defined(__MPI)
    IF (ionode) THEN
#endif
       !
       !
       !if (.not. allocated(w_T_beta_store)) print *, "aaaaaaaaaaaaaa"
       filename = trim(w_T_prefix) // trim(int_to_char(LR_polarization))
       !
       WRITE(stdout,'(/,/5x,"Reading Pre-calculated lanczos coefficents from ",A50)') filename
       !
       INQUIRE (file = filename, exist = exst)
       !
       IF (.not.exst) CALL errore(' read_beta_gamma_z ','Stage 1 Lanczos coefficents not found ',1)
       !
       !
       OPEN (158, file = filename, form = 'formatted', status = 'old')
       INQUIRE (file = filename, opened = exst)
       IF (.not.exst) CALL errore(' read_beta_gamma_z ','Stage 1 Lanczos coefficents can not be opened ',1)
       !
       !
       READ(158,*,end=301,err=302) iter_restart
       !print *,iter_restart
       !write(stdout,'(/,5X,"Reading of precalculated Lanczos Matrix")')
       IF (iter_restart < itermax) CALL errore ('read_beta_gamma_z', 'Lanczos iteration mismatch', 1 )
       !
       IF (.not. allocated(w_T_beta_store))  ALLOCATE(w_T_beta_store(iter_restart))
       IF (.not. allocated(w_T_gamma_store))  ALLOCATE(w_T_gamma_store(iter_restart))
       READ(158,*,end=301,err=303) w_T_norm0_store
       READ(158,*,end=301,err=303) discard
       READ(158,*,end=301,err=303) discard
       READ(158,*,end=301,err=303) discard
       READ(158,*,end=301,err=303) discard
       READ(158,*,end=301,err=303) discard
       READ(158,*,end=301,err=303) discard
       READ(158,*,end=301,err=303) discard
       !print *, discard
       !
       !write(stdout,'("--------------Lanczos Matrix-------------------")')
       DO i=1,itermax
          !
          !print *, "Iter=",i
          READ(158,*,end=301,err=303) w_T_beta_store(i)
          !print *, w_T_beta_store(i)
          READ(158,*,end=301,err=303) w_T_gamma_store(i)
          !print *, w_T_gamma_store(i)
          DO j=1,w_T_npol
             READ(158,*,end=301,err=303) w_T_zeta_store(j,i)
          ENDDO
          !print *, discard2(:)
          !
       ENDDO
       !print *, "closing file"
       !
       CLOSE(158)
       !
       !print *, "starting broadcast"
#if defined(__MPI)
    ENDIF
    CALL mp_barrier(world_comm)
    CALL mp_bcast (w_T_beta_store(:), ionode_id, world_comm)
    CALL mp_bcast (w_T_gamma_store(:), ionode_id, world_comm)
    CALL mp_bcast (w_T_zeta_store(:,:), ionode_id, world_comm)
    CALL mp_bcast (w_T_norm0_store, ionode_id, world_comm)
#endif
    !print *, "broadcast complete"
    WRITE(stdout,'(5x,I8,1x,"steps succesfully read for polarization index",1x,I3)') itermax,LR_polarization
    CALL stop_clock( 'post-processing' )
    RETURN
301 CALL errore ('read_beta_gamma_z', 'File is corrupted, no data', i )
302 CALL errore ('read_beta_gamma_z', 'File is corrupted, itermax not found', 1 )
303 CALL errore ('read_beta_gamma_z', 'File is corrupted, data number follows:', i )
  END SUBROUTINE read_wT_beta_gamma_z
  !-----------------------------------------------------------------------
  SUBROUTINE lr_calc_w_T()
    !---------------------------------------------------------------------
    ! ... calculates the w_T from equation (\freq - L ) e_1
    ! ... by solving tridiagonal problem for each value of freq
    !---------------------------------------------------------------------
    !
    USE lr_variables,         ONLY : itermax,beta_store,gamma_store, &
         LR_polarization,charge_response, n_ipol, &
         itermax_int,project,sum_rule
    USE fft_base,             ONLY : dfftp
    USE noncollin_module,     ONLY : nspin_mag

    !
    IMPLICIT NONE
    !
    !integer, intent(in) :: freq ! Input : The frequency identifier (1 o 5) for w_T
    COMPLEX(kind=dp), ALLOCATABLE :: a(:), b(:), c(:),r(:)
    real(kind=dp) :: average,av_amplitude
    COMPLEX(kind=dp) :: norm
    !
    INTEGER :: i, info,ip,ip2 !used for error reporting
    INTEGER :: counter
    LOGICAL :: skip
    !Solver:
    real(kind=dp), EXTERNAL :: ddot
    COMPLEX(kind=dp), EXTERNAL :: zdotc
    !
    CALL start_clock( 'post-processing' )
    IF (lr_verbosity > 5) THEN
       WRITE(stdout,'("<lr_calc_w_T>")')
    ENDIF
    !IF (omeg == 0.D0) THEN
    !   CALL stop_clock( 'post-processing' )
    !   RETURN
    !END IF
    !
    ALLOCATE(a(itermax_int))
    ALLOCATE(b(itermax_int))
    ALLOCATE(c(itermax_int))
    ALLOCATE(r(itermax_int))
    a(:) = (0.0d0,0.0d0)
    b(:) = (0.0d0,0.0d0)
    c(:) = (0.0d0,0.0d0)
    w_T(:) = (0.0d0,0.0d0)
    !
    WRITE(stdout,'(/,5X,"Calculating response coefficients")')
     
    if (sum_rule == -2 ) THEN
    WRITE(stdout,'(/,5X,"S - 2 will be attempted")')
    resonance_condition=.true.
    IF (allocated(rho_1_tot)) DEALLOCATE (rho_1_tot)
    IF (.not. allocated(rho_1_tot_im)) ALLOCATE(rho_1_tot_im(dfftp%nnr,nspin_mag))
    rho_1_tot_im(:,:)=cmplx(0.0d0,0.0d0,dp)
    !
    !
    !
    ! prepare tridiagonal (w-L) for the given polarization
    !
    omeg=0
    w_t(:)=(0.0d0,0.0d0)
    do while (omeg < 100) 
    omeg=omeg+0.01
    a(:) = cmplx(omeg,epsil,dp)
    !
       !Read the actual iterations
       DO i=1,itermax
          !
          !b(i)=-w_T_beta_store(i)
          !c(i)=-w_T_gamma_store(i)
          b(i)=cmplx(-w_T_beta_store(i),0.0d0,dp)
          c(i)=cmplx(-w_T_gamma_store(i),0.0d0,dp)
          !
       ENDDO
       IF (itermax_int>itermax .and. itermax > 151) THEN
          !calculation of the average
          !OBM: (I am using the code from tddfpt_pp, I am not very confortable
          ! with the mechanism discarding the "bad" points.
          average=0.d0
          av_amplitude=0.d0
          counter=0
          skip=.false.
          !
          DO i=151,itermax
             !
             IF (skip .eqv. .true.) THEN
                skip=.false.
                CYCLE
             ENDIF
             !
             IF (mod(i,2)==1) THEN
                !
                IF ( i/=151 .and. abs( w_T_beta_store(i)-average/counter ) > 2.d0 ) THEN
                   !
                   !if ( i.ne.151 .and. counter == 0) counter = 1
                   skip=.true.
                   !
                ELSE
                   !
                   average=average+w_T_beta_store(i)
                   av_amplitude=av_amplitude+w_T_beta_store(i)
                   counter=counter+1
                   !print *, "t1 ipol",ip,"av_amp",av_amplitude(ip)
                   !
                ENDIF
                !
             ELSE
                !
                IF ( i/=151 .and. abs( w_T_beta_store(i)-average/counter ) > 2.d0 ) THEN
                   !
                   !if ( i.ne.151 .and. counter == 0) counter = 1
                   skip=.true.
                   !
                ELSE
                   !
                   average=average+w_T_beta_store(i)
                   av_amplitude=av_amplitude-w_T_beta_store(i)
                   counter=counter+1
                   !print *, "t2 ipol",ip,"av_amp",av_amplitude(ip)
                   !
                ENDIF
                !
             ENDIF
             !
             !
          ENDDO
          average=average/counter
          av_amplitude=av_amplitude/counter
          !
          !
          !extrapolated part of b and c
          DO i=itermax,itermax_int
             !
             IF (mod(i,2)==1) THEN
                !
                b(i)=cmplx((-average-av_amplitude),0.0d0,dp)
                c(i)=b(i)
                !
             ELSE
                !
                b(i)=cmplx((-average+av_amplitude),0.0d0,dp)
                c(i)=b(i)
                !
             ENDIF
             !
          ENDDO

    ENDIF
    !
    r(:) =(0.0d0,0.0d0)
    r(1)=(1.0d0,0.0d0)
    !
    ! solve the equation
    !
    CALL zgtsv(itermax_int,1,b,a,c,r(:),itermax_int,info)
    IF(info /= 0) CALL errore ('calc_w_T', 'unable to solve tridiagonal system', 1 )
    !
    !Check if we are close to a resonance
    !
    norm=sum(abs(aimag(r(:))/dble(r(:))))
    norm=norm/(1.0d0*itermax_int)
    ! Is this a correct way to find imaginary part?
    IF (abs(norm) > 0.1) THEN
         w_t(:)=w_t(:)+r(:)/omeg
    END IF
    !print *,"norm",norm
    !
    !
    END DO
    ELSE IF (sum_rule == -99 ) THEN
    !
    !
    !
    ! prepare tridiagonal (w-L) for the given polarization
    !
    a(:) = cmplx(omeg,epsil,dp)
    !
    IF (charge_response == 2) THEN
       DO i=1,itermax-1
          !
          !Memory mapping in case of selected polarization direction
          info=1
          IF ( n_ipol /= 1 ) info=LR_polarization
          !
          b(i)=-beta_store(info,i)
          c(i)=-gamma_store(info,i)
          !
       ENDDO
    ENDIF
    IF (charge_response == 1) THEN
       !Read the actual iterations
       DO i=1,itermax
          !
          !b(i)=-w_T_beta_store(i)
          !c(i)=-w_T_gamma_store(i)
          b(i)=cmplx(-w_T_beta_store(i),0.0d0,dp)
          c(i)=cmplx(-w_T_gamma_store(i),0.0d0,dp)
          !
       ENDDO
       IF (itermax_int>itermax .and. itermax > 151) THEN
          !calculation of the average
          !OBM: (I am using the code from tddfpt_pp, I am not very confortable
          ! with the mechanism discarding the "bad" points.
          average=0.d0
          av_amplitude=0.d0
          counter=0
          skip=.false.
          !
          DO i=151,itermax
             !
             IF (skip .eqv. .true.) THEN
                skip=.false.
                CYCLE
             ENDIF
             !
             IF (mod(i,2)==1) THEN
                !
                IF ( i/=151 .and. abs( w_T_beta_store(i)-average/counter ) > 2.d0 ) THEN
                   !
                   !if ( i.ne.151 .and. counter == 0) counter = 1
                   skip=.true.
                   !
                ELSE
                   !
                   average=average+w_T_beta_store(i)
                   av_amplitude=av_amplitude+w_T_beta_store(i)
                   counter=counter+1
                   !print *, "t1 ipol",ip,"av_amp",av_amplitude(ip)
                   !
                ENDIF
                !
             ELSE
                !
                IF ( i/=151 .and. abs( w_T_beta_store(i)-average/counter ) > 2.d0 ) THEN
                   !
                   !if ( i.ne.151 .and. counter == 0) counter = 1
                   skip=.true.
                   !
                ELSE
                   !
                   average=average+w_T_beta_store(i)
                   av_amplitude=av_amplitude-w_T_beta_store(i)
                   counter=counter+1
                   !print *, "t2 ipol",ip,"av_amp",av_amplitude(ip)
                   !
                ENDIF
                !
             ENDIF
             !
             !
          ENDDO
          average=average/counter
          av_amplitude=av_amplitude/counter
          !
          !
          WRITE(stdout,'(/,5X,"Charge Response extrapolation average: ",E15.5)') average
          WRITE(stdout,'(5X,"Charge Response extrapolation oscillation amplitude: ",E15.5)') av_amplitude

          !extrapolated part of b and c
          DO i=itermax,itermax_int
             !
             IF (mod(i,2)==1) THEN
                !
                b(i)=cmplx((-average-av_amplitude),0.0d0,dp)
                c(i)=b(i)
                !
             ELSE
                !
                b(i)=cmplx((-average+av_amplitude),0.0d0,dp)
                c(i)=b(i)
                !
             ENDIF
             !
          ENDDO

       ENDIF
    ENDIF
    !
    r(:) =(0.0d0,0.0d0)
    r(1)=(1.0d0,0.0d0)
    !
    ! solve the equation
    !
    CALL zgtsv(itermax_int,1,b,a,c,r(:),itermax_int,info)
    IF(info /= 0) CALL errore ('calc_w_T', 'unable to solve tridiagonal system', 1 )
    w_t(:)=r(:)
    !
    !Check if we are close to a resonance
    !
    norm=sum(abs(aimag(w_T(:))/dble(w_T(:))))
    norm=norm/(1.0d0*itermax_int)
    !print *,"norm",norm
    IF (abs(norm) > 0.1) THEN
       resonance_condition=.true.
       IF (allocated(rho_1_tot)) DEALLOCATE (rho_1_tot)
       IF (.not. allocated(rho_1_tot_im)) ALLOCATE(rho_1_tot_im(dfftp%nnr,nspin_mag))
       rho_1_tot_im(:,:)=cmplx(0.0d0,0.0d0,dp)
    ELSE
       resonance_condition=.false.
       IF (allocated(rho_1_tot_im)) DEALLOCATE (rho_1_tot_im)
       IF (.not. allocated(rho_1_tot)) ALLOCATE(rho_1_tot(dfftp%nnr,nspin_mag))
       rho_1_tot(:,:)=0.0d0
    ENDIF
    IF (resonance_condition)  THEN
       WRITE(stdout,'(5X,"Resonance frequency mode enabled")')
       !write(stdout,'(5X,"Response charge density multiplication factor=",E15.8)') 1.0d0/epsil**2
    ENDIF
    !
    ! normalize so that the final charge densities are normalized
    !
    norm=zdotc(itermax_int,w_T(:),1,w_T(:),1)
    WRITE(stdout,'(5X,"Charge Response renormalization factor: ",2(E15.5,1x))') norm
    !w_T(:)=w_T(:)/norm
    !norm=sum(w_T(:))
    !write(stdout,'(3X,"Initial sum of lanczos vectors",F8.5)') norm
    !w_T(:)=w_T(:)/norm
    !
    !
    !Calculate polarizability tensor in the case of projection
    IF (project) THEN
       DO ip=1,w_T_npol
          !
          chi(LR_polarization,ip)=ZDOTC(itermax,w_T_zeta_store(ip,:),1,w_T(:),1)
          chi(LR_polarization,ip)=chi(LR_polarization,ip)*cmplx(w_T_norm0_store,0.0d0,dp)
          !
          WRITE(stdout,'(5X,"Chi_",I1,"_",I1,"=",2(E15.5,1x))') LR_polarization,ip,chi(LR_polarization,ip)
       ENDDO
    ENDIF
    ENDIF
    !
    !
    DEALLOCATE(a)
    DEALLOCATE(b)
    DEALLOCATE(c)
    DEALLOCATE(r)
    !
    IF ( lr_verbosity > 3 ) THEN
       WRITE(stdout,'("--------Lanczos weight coefficients in the direction ", &
            &  I1," for freq=",D15.8," Ry ----------")') LR_polarization, omeg
       DO i=1,itermax
          WRITE(stdout,'(I5,3X,2D15.8)') i, w_T(i)
       ENDDO
       WRITE(stdout,'("------------------------------------------------------------------------")')
       WRITE(stdout,'("NR1=",I15," NR2=",I15," NR3=",I15)') dfftp%nr1, dfftp%nr2, dfftp%nr3
       WRITE(stdout,'("------------------------------------------------------------------------")')
    ENDIF
    CALL stop_clock( 'post-processing' )
    RETURN
    !
  END SUBROUTINE lr_calc_w_T
  !-----------------------------------------------------------------------
  SUBROUTINE lr_calc_F(evc1)
    !-------------------------------------------------------------------------------
    ! Calculates the projection of empty states to response orbitals
    !
    USE lsda_mod,                 ONLY : nspin
    USE mp,                       ONLY : mp_sum
    USE mp_global,                ONLY : inter_pool_comm, intra_bgrp_comm
    USE uspp,                     ONLY : okvan,qq,vkb
    USE wvfct,                    ONLY : wg,nbnd,npwx
    USE uspp_param,               ONLY : upf, nh
    USE becmod,                   ONLY : becp,calbec
    USE ions_base,                ONLY : ityp,nat,ntyp=>nsp
    USE realus,                   ONLY : real_space_debug,invfft_orbital_gamma,calbec_rs_gamma
    USE gvect,                    ONLY : gstart
    USE klist,                    ONLY : nks, ngk
    USE lr_variables,             ONLY : lr_verbosity, itermax, LR_iteration, LR_polarization, &
                                         project,evc0_virt,F,nbnd_total,n_ipol, becp1_virt

    IMPLICIT NONE
    !
    !input
    COMPLEX(kind=dp), INTENT(in) :: evc1(npwx,nbnd,nks)
    !
    !internal variables
    INTEGER :: ibnd_occ,ibnd_virt,ipol
    real(kind=dp) :: w1,w2,scal
    INTEGER :: ir,ik,ibnd,jbnd,ig,ijkb0,np,na,ijh,ih,jh,ikb,jkb,ispin
    !complex(kind=dp) :: SSUM
    real(kind=dp)     :: SSUM
    !
    !functions
    real(kind=dp), EXTERNAL    :: DDOT
    !complex(kind=dp), external    :: ZDOTC
    !
    scal=0.0d0
    !
    ! I calculate the projection <virtual|\rho^\prime|occupied> from
    ! F=2(<evc0_virt|evc1>+\sum Q <evc0_virt|beta><beta|evc1>
    IF ( .not. project) RETURN
    IF (n_ipol>1) THEN
       ipol=LR_polarization
    ELSE
       ipol=1
    ENDIF
    IF (okvan) THEN
       !BECP initialisation for evc1
       IF (real_space_debug >6) THEN
          DO ibnd=1,nbnd,2
             CALL invfft_orbital_gamma(evc1(:,:,1),ibnd,nbnd)
             CALL calbec_rs_gamma(ibnd,nbnd,becp%r)
          ENDDO
       ELSE
          CALL calbec(ngk(1), vkb, evc1(:,:,1), becp)
       ENDIF
    ENDIF
    !
!!! Actual projection starts here
    !
    DO ibnd_occ=1,nbnd
       DO ibnd_virt=1,(nbnd_total-nbnd)
          !
          !ultrasoft part
          !
          IF (okvan) THEN
             !initialization
             scal = 0.0d0
             !
             !Calculation of  qq<evc0|beta><beta|evc1>
             !
             w1 = wg(ibnd,1)
             ijkb0 = 0
             !
             DO np = 1, ntyp
                !
                IF ( upf(np)%tvanp ) THEN
                   !
                   DO na = 1, nat
                      !
                      IF ( ityp(na) == np ) THEN
                         !
                         ijh = 1
                         !
                         DO ih = 1, nh(np)
                            !
                            ikb = ijkb0 + ih
                            !
                            !  <beta_i|beta_i> terms
                            !
                            scal = scal + qq(ih,ih,np) *1.d0 *  becp%r(ikb,ibnd_occ) * becp1_virt(ikb,ibnd_virt)
                            !
                            ijh = ijh + 1
                            !
                            ! <beta_i|beta_j> terms
                            !
                            DO jh = ( ih + 1 ), nh(np)
                               !
                               jkb = ijkb0 + jh
                               !
                               scal = scal + qq(ih,jh,np) *1.d0  * (becp%r(ikb,ibnd_occ) * becp1_virt(jkb,ibnd_virt)+&
                                    becp%r(jkb,ibnd_occ) * becp1_virt(ikb,ibnd_virt))
                               !
                               ijh = ijh + 1
                               !
                            ENDDO
                            !
                         ENDDO
                         !
                         ijkb0 = ijkb0 + nh(np)
                         !
                      ENDIF
                      !
                   ENDDO
                   !
                ELSE
                   !
                   DO na = 1, nat
                      !
                      IF ( ityp(na) == np ) ijkb0 = ijkb0 + nh(np)
                      !
                   ENDDO
                   !
                ENDIF
                !
             ENDDO
             !
          ENDIF
          ! US part finished
          !first part
          ! the dot  product <evc1|evc0> taken from lr_dot
          SSUM=(2.D0*wg(ibnd_occ,1)*DDOT(2*ngk(1),evc0_virt(:,ibnd_virt,1),1,evc1(:,ibnd_occ,1),1))
          IF (gstart==2) SSUM = SSUM - (wg(ibnd_occ,1)*dble(evc1(1,ibnd_occ,1))*dble(evc0_virt(1,ibnd_virt,1)))
          !US contribution
          SSUM=SSUM+scal
#if defined(__MPI)
          CALL mp_sum(SSUM, intra_bgrp_comm)
#endif
          IF(nspin/=2) SSUM=SSUM/2.0D0
          !
          !
          !and finally (note:parellization handled in dot product, each node has the copy of F)
          !
          F(ibnd_occ,ibnd_virt,ipol)=F(ibnd_occ,ibnd_virt,ipol)+cmplx(SSUM,0.0d0,dp)*w_T(LR_iteration)
       ENDDO
    ENDDO
  END SUBROUTINE lr_calc_F
  !-------------------------------------------------------------------------------
  SUBROUTINE lr_calc_R()
    !-------------------------------------------------------------------------------
    ! Calculates the oscillator strengths
    !
    USE lsda_mod,                 ONLY : nspin
    USE mp,                       ONLY : mp_sum
    USE mp_global,                ONLY : inter_pool_comm, intra_bgrp_comm
    USE uspp,                     ONLY : okvan,qq,vkb
    USE wvfct,                    ONLY : wg,nbnd,npwx
    USE uspp_param,               ONLY : upf, nh
    USE becmod,                   ONLY : becp,calbec
    USE ions_base,                ONLY : ityp,nat,ntyp=>nsp
    USE realus,                   ONLY : real_space_debug,invfft_orbital_gamma,calbec_rs_gamma
    USE gvect,                    ONLY : gstart
    USE klist,                    ONLY : nks, ngk
    USE lr_variables,             ONLY : lr_verbosity, itermax, LR_iteration, LR_polarization, &
         project,evc0_virt,R,nbnd_total,n_ipol, becp1_virt,d0psi

    IMPLICIT NONE
    !
    !
    !internal variables
    INTEGER :: ibnd_occ,ibnd_virt,ipol
    real(kind=dp)     :: SSUM
    !
    !functions
    real(kind=dp), EXTERNAL    :: DDOT
    !
    DO ipol=1,n_ipol
       DO ibnd_occ=1,nbnd
          DO ibnd_virt=1,(nbnd_total-nbnd)
             ! the dot  product <evc0|sd0psi> taken from lr_dot
             SSUM=(2.D0*wg(ibnd_occ,1)*DDOT(2*ngk(1),evc0_virt(:,ibnd_virt,1),1,d0psi(:,ibnd_occ,1,ipol),1))
             IF (gstart==2) SSUM = SSUM - (wg(ibnd_occ,1)*dble(d0psi(1,ibnd_occ,1,ipol))*dble(evc0_virt(1,ibnd_virt,1)))
#if defined(__MPI)
             CALL mp_sum(SSUM, intra_bgrp_comm)
#endif
             IF(nspin/=2) SSUM=SSUM/2.0D0
             !
             R(ibnd_occ,ibnd_virt,ipol)=cmplx(SSUM,0.0d0,dp)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE lr_calc_R
  !-----------------------------------------------------------------------
  SUBROUTINE lr_dump_rho(plot_type)
    !-----------------------------------------------------------------------
    !
    ! This routine acts as an interface to all of the various lr_dump_* routines
    ! and calls the apropriate ones, based on the plot type specified in input
    ! and whether the resonance condition is satisfied.
    !
    !-----------------------------------------------------------------------
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)             :: plot_type
    !
    IF (resonance_condition) THEN
       !response charge density, absorbtive
       IF (plot_type == 1 .or. plot_type == 5) &
            CALL lr_dump_rho_tot_xyzd(aimag(rho_1_tot_im(:,1)),"absorbtive")
       IF (plot_type == 2 .or. plot_type == 5) &
            CALL lr_dump_rho_tot_xcrys(aimag(rho_1_tot_im(:,1)),"absorbtive")
       IF (plot_type == 3 .or. plot_type == 5) &
            CALL lr_dump_rho_tot_cube(aimag(rho_1_tot_im(:,1)),"absorbtive")
       !response charge density, dispersive
       IF (plot_type == 1 .or. plot_type == 5) &
            CALL lr_dump_rho_tot_xyzd(dble(rho_1_tot_im(:,1)),"dispersive")
       IF (plot_type == 2 .or. plot_type == 5) &
            CALL lr_dump_rho_tot_xcrys(dble(rho_1_tot_im(:,1)),"dispersive")
       IF (plot_type == 3 .or. plot_type == 5) &
            CALL lr_dump_rho_tot_cube(dble(rho_1_tot_im(:,1)),"dispersive")
    ELSE
       IF (plot_type == 1 .OR. plot_type == 5) &
            & CALL lr_dump_rho_tot_xyzd(rho_1_tot(:,1),"summed-rho")
       IF (plot_type == 2 .OR. plot_type == 5) &
            & CALL lr_dump_rho_tot_xcrys(rho_1_tot(:,1),"summed-rho")
       IF (plot_type == 3 .OR. plot_type == 5) &
            & CALL lr_dump_rho_tot_cube(rho_1_tot(:,1),"summed-rho")
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE lr_dump_rho
  !-----------------------------------------------------------------------
  SUBROUTINE lr_dump_rho_tot_compat1()
    !-----------------------------------------------------------------------
    ! dump a density file in a format compatible to type 1 charge response calculation
    !-----------------------------------------------------------------------
    USE io_files,              ONLY : prefix
    USE lr_variables,          ONLY : LR_polarization, LR_iteration, cube_save
    USE fft_base,              ONLY : dfftp
    USE mp_global,             ONLY : inter_pool_comm, intra_bgrp_comm

    IMPLICIT NONE
    CHARACTER(len=6), EXTERNAL :: int_to_char
    !
    !Local
    CHARACTER (len=80):: filename
    real(kind=dp), ALLOCATABLE :: rho_sum_resp_x(:),rho_sum_resp_y(:),rho_sum_resp_z(:)
    INTEGER ir,i,j,k
    CALL start_clock( 'post-processing' )
    IF (lr_verbosity > 5) WRITE(stdout,'("<lr_dump_rho_tot_compat1>")')
#if defined(__MPI)
    IF (ionode) THEN
#endif
       !
       IF ( .not. allocated(cube_save) ) CALL lr_set_boxes_density()
       ALLOCATE( rho_sum_resp_x( dfftp%nr1 ) )
       ALLOCATE( rho_sum_resp_y( dfftp%nr2 ) )
       ALLOCATE( rho_sum_resp_z( dfftp%nr3 ) )
       !
       rho_sum_resp_x = 0.D0
       rho_sum_resp_y = 0.D0
       rho_sum_resp_z = 0.D0
       !
       DO ir=1,dfftp%nnr
          !
          i=cube_save(ir,1)+1
          j=cube_save(ir,2)+1
          k=cube_save(ir,3)+1
          !
          rho_sum_resp_x(i)=rho_sum_resp_x(i)+rho_1_tot(ir,1)
          rho_sum_resp_y(j)=rho_sum_resp_y(j)+rho_1_tot(ir,1)
          rho_sum_resp_z(k)=rho_sum_resp_z(k)+rho_1_tot(ir,1)
          !
       ENDDO
       !
       !
       filename = trim(prefix) // "-summed-density-pol" //trim(int_to_char(LR_polarization))// "_x"
       !
       OPEN (158, file = filename, form = 'formatted', status = 'unknown', position = 'append')
       !
       DO i=1,dfftp%nr1
          WRITE(158,*) rho_sum_resp_x(i)
       ENDDO
       !
       CLOSE(158)
       !
       filename = trim(prefix) // "-summed-density-pol" //trim(int_to_char(LR_polarization))// "_y"
       !
       OPEN (158, file = filename, form = 'formatted', status = 'unknown', position = 'append')
       !
       DO i=1,dfftp%nr2
          WRITE(158,*) rho_sum_resp_y(i)
       ENDDO
       !
       CLOSE(158)
       !
       filename = trim(prefix) // "-summed-density-pol" //trim(int_to_char(LR_polarization))// "_z"
       !
       OPEN (158, file = filename, form = 'formatted', status = 'unknown', position = 'append')
       !
       DO i=1,dfftp%nr3
          WRITE(158,*) rho_sum_resp_z(i)
       ENDDO
       !
       CLOSE(158)
       DEALLOCATE( rho_sum_resp_x )
       DEALLOCATE( rho_sum_resp_y )
       DEALLOCATE( rho_sum_resp_z )
       !
#if defined(__MPI)
    ENDIF
#endif
    CALL stop_clock( 'post-processing' )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE lr_dump_rho_tot_compat1
  !-----------------------------------------------------------------------
  SUBROUTINE lr_dump_rho_tot_cube(rho,identifier)
    !-----------------------------------------------------------------------
    ! dump a density file in the gaussian cube format. "Inspired" by
    ! Modules/cube.f90 :)
    !-----------------------------------------------------------------------
    USE io_files,              ONLY : prefix
    USE lr_variables,          ONLY : LR_polarization, LR_iteration, cube_save
    USE cell_base
    USE ions_base,                ONLY : nat, ityp, atm, ntyp => nsp, tau
    USE mp,                   ONLY : mp_barrier, mp_sum, mp_bcast, mp_get
    USE mp_world,             ONLY : world_comm
    USE mp_global,            ONLY : me_image, intra_image_comm, me_bgrp, &
                                     nproc_bgrp, intra_bgrp_comm, my_bgrp_id

    USE constants,            ONLY : BOHR_RADIUS_ANGS
    USE fft_base,             ONLY : dfftp !this contains dfftp%npp (number of z planes per processor
    ! and dfftp%ipp (offset of the first z plane of the processor

    !
    IMPLICIT NONE
    !
    real (kind=dp), INTENT(in)   :: rho(:)
    CHARACTER(len=10), INTENT(in) :: identifier
    !
    CHARACTER(len=80) :: filename
    !
    CHARACTER(len=6), EXTERNAL :: int_to_char
    !
    !Local
    INTEGER          :: i, nt, i1, i2, i3, at_num, iopool_id, ldr, kk, ionode_pool, six_count
    real(DP)    :: at_chrg, tpos(3), inpos(3)
    REAL(DP), ALLOCATABLE :: rho_plane(:),rho_temp(:)
    INTEGER,  ALLOCATABLE :: kowner(:)
    !
    INTEGER, EXTERNAL:: atomic_number
    !

    CALL start_clock( 'post-processing' )
    IF (lr_verbosity > 5) WRITE(stdout,'("<lr_dump_rho_tot_cube>")')
    !
    six_count=0
#if defined(__MPI)
    ALLOCATE( rho_temp(dfftp%npp(1)+1) )
    IF (ionode) THEN
       filename = trim(prefix) // "-" // identifier // "-pol" //trim(int_to_char(LR_polarization))// ".cube"
       WRITE(stdout,'(/5X,"Writing Cube file for response charge density")')
       !write(stdout, *) filename
       !write(stdout,'(5X,"|rho|=",D15.8)') rho_sum
       OPEN (158, file = filename, form = 'formatted', status = 'replace', err=501)

       !C     WRITE A FORMATTED 'DENSITY-STYLE' CUBEFILE VERY SIMILAR
       !C     TO THOSE CREATED BY THE GAUSSIAN PROGRAM OR THE CUBEGEN UTILITY.
       !C     THE FORMAT IS AS FOLLOWS (LAST CHECKED AGAINST GAUSSIAN 98):
       !C
       !C     LINE   FORMAT      CONTENTS
       !C     ===============================================================
       !C      1     A           TITLE
       !C      2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
       !C      3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
       !C      4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
       !C      #ATOMS LINES OF ATOM COORDINATES:
       !C      ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
       !C      REST: 6E13.5      CUBE DATA
       !C
       !C     ALL COORDINATES ARE GIVEN IN ATOMIC UNITS.

       WRITE(158,*) 'Cubfile created from TDDFPT calculation'
       WRITE(158,*) identifier
       !                        origin is forced to (0.0,0.0,0.0)
       WRITE(158,'(I5,3F12.6)') nat, 0.0d0, 0.0d0, 0.0d0
       WRITE(158,'(I5,3F12.6)') dfftp%nr1, (alat*at(i,1)/dble(dfftp%nr1),i=1,3)
       WRITE(158,'(I5,3F12.6)') dfftp%nr2, (alat*at(i,2)/dble(dfftp%nr2),i=1,3)
       WRITE(158,'(I5,3F12.6)') dfftp%nr3, (alat*at(i,3)/dble(dfftp%nr3),i=1,3)

       DO i=1,nat
          nt = ityp(i)
          ! find atomic number for this atom.
          at_num = atomic_number(trim(atm(nt)))
          at_chrg= dble(at_num)
          ! at_chrg could be alternatively set to valence charge
          ! positions are in cartesian coordinates and a.u.
          !
          ! wrap coordinates back into cell.
          tpos = matmul( transpose(bg), tau(:,i) )
          tpos = tpos - nint(tpos - 0.5d0)
          inpos = alat * matmul( at, tpos )
          WRITE(158,'(I5,5F12.6)') at_num, at_chrg, inpos
       ENDDO
    ENDIF
    ! Header is complete, now dump the charge density, as derived from xyzd subroutine
    ALLOCATE( rho_plane( dfftp%nr3x ) )
    !ALLOCATE( kowner( dfftp%nr3 ) )
    !
    ! ... find the index of the pool that will write rho
    !
    IF ( ionode ) iopool_id = my_bgrp_id
    !
    CALL mp_bcast( iopool_id, ionode_id, intra_image_comm )
    !
    ! ... find the index of the ionode within its own pool
    !
    IF ( ionode ) ionode_pool = me_bgrp
    !
    CALL mp_bcast( ionode_pool, ionode_id, intra_image_comm )
    !
    ! ... find out the owner of each "z" plane
    !
    !
    !IF (nproc_bgrp > 1) THEN
    ! DO i = 1, nproc_bgrp
    !    !
    !    kowner( (dfftp%ipp(i)+1):(dfftp%ipp(i)+dfftp%npp(i)) ) = i - 1
    !    !
    ! END DO
    !ELSE
    ! kowner = ionode_id
    !ENDIF
    ldr = dfftp%nr1x*dfftp%nr2x
    !
    !
    ! Each processor is on standby to send its plane to ionode
    DO i1 = 1, dfftp%nr1
       !
       DO i2 = 1, dfftp%nr2
          !
          !Parallel gather of Z plane
          rho_plane(:)=0
          DO i = 1, nproc_bgrp
             rho_temp(:)=0
             IF( (i-1) == me_bgrp ) THEN
                !
                !
                DO  i3=1, dfftp%npp(i)
                   !
                   rho_temp(i3) = rho(i1+(i2-1)*dfftp%nr1x+(i3-1)*ldr)
                   !
                   !
                ENDDO
                !print *, "get 1=",rho_plane(1)," 2=",rho_plane(2)," ",dfftp%npp(i),"=",rho_plane(dfftp%npp(i))
             ENDIF
             !call mp_barrier( world_comm )
             IF ( my_bgrp_id == iopool_id ) &
                                !Send plane to ionode
                                ! Send and recieve rho_plane,
                  CALL mp_get( rho_temp, rho_temp, &
                  me_bgrp, ionode_pool, (i-1), i-1, intra_bgrp_comm )

             !
             !call mp_barrier( world_comm )
             IF(ionode) THEN
                rho_plane( (dfftp%ipp(i)+1):(dfftp%ipp(i)+dfftp%npp(i)) ) = rho_temp(1:dfftp%npp(i))
                !print *, "get (",dfftp%ipp(i)+1,")=",rho_plane(dfftp%ipp(i)+1)," (",dfftp%ipp(i)+dfftp%npp(i),")=",rho_plane(dfftp%ipp(i)+dfftp%npp(i))
                !print *, "data of proc ",i," written I2=",i2,"I1=",i1
             ENDIF
          ENDDO
          ! End of parallel send
          IF (ionode) THEN
             DO  i3=1, dfftp%nr3
                six_count=six_count+1
                WRITE(158,'(E13.5)',advance='no') rho_plane(i3)
                IF (six_count == 6 ) THEN
                   WRITE(158,'("")')
                   six_count=0
                ENDIF
                !print *, rho_plane(i3)
             ENDDO
          ENDIF
          CALL mp_barrier( world_comm )
       ENDDO
    ENDDO
    !
    DEALLOCATE( rho_plane )
    DEALLOCATE( rho_temp )

    IF (ionode) CLOSE(158)

    CALL mp_barrier( world_comm )

#else
    !
    !

    filename = trim(prefix) // "-" // identifier // "-pol" //trim(int_to_char(LR_polarization))// ".cube"
    WRITE(stdout,'(/5X,"Writing Cube file for response charge density")')
    !write(stdout, *) filename
    !write(stdout,'(5X,"|rho|=",D15.8)') rho_sum
    OPEN (158, file = filename, form = 'formatted', status = 'replace', err=501)

    !C     WRITE A FORMATTED 'DENSITY-STYLE' CUBEFILE VERY SIMILAR
    !C     TO THOSE CREATED BY THE GAUSSIAN PROGRAM OR THE CUBEGEN UTILITY.
    !C     THE FORMAT IS AS FOLLOWS (LAST CHECKED AGAINST GAUSSIAN 98):
    !C
    !C     LINE   FORMAT      CONTENTS
    !C     ===============================================================
    !C      1     A           TITLE
    !C      2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
    !C      3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
    !C      4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
    !C      #ATOMS LINES OF ATOM COORDINATES:
    !C      ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
    !C      REST: 6E13.5      CUBE DATA
    !C
    !C     ALL COORDINATES ARE GIVEN IN ATOMIC UNITS.

    WRITE(158,*) 'Cubefile created from TDDFPT calculation'
    WRITE(158,*) identifier
    !                        origin is forced to (0.0,0.0,0.0)
    WRITE(158,'(I5,3F12.6)') nat, 0.0d0, 0.0d0, 0.0d0
    WRITE(158,'(I5,3F12.6)') dfftp%nr1, (alat*at(i,1)/dble(dfftp%nr1),i=1,3)
    WRITE(158,'(I5,3F12.6)') dfftp%nr2, (alat*at(i,2)/dble(dfftp%nr2),i=1,3)
    WRITE(158,'(I5,3F12.6)') dfftp%nr3, (alat*at(i,3)/dble(dfftp%nr3),i=1,3)

    DO i=1,nat
       nt = ityp(i)
       ! find atomic number for this atom.
       at_num = atomic_number(trim(atm(nt)))
       at_chrg= dble(at_num)
       ! at_chrg could be alternatively set to valence charge
       ! positions are in cartesian coordinates and a.u.
       !
       ! wrap coordinates back into cell.
       tpos = matmul( transpose(bg), tau(:,i) )
       tpos = tpos - nint(tpos - 0.5d0)
       inpos = alat * matmul( at, tpos )
       WRITE(158,'(I5,5F12.6)') at_num, at_chrg, inpos
    ENDDO
    i=0
    DO i1=1,dfftp%nr1
       DO i2=1,dfftp%nr2
          DO i3=1,dfftp%nr3
             !i(i3-1)*dfftp%nr1x*dfftp%nr2x+(i2-1)*dfftp%nr1x+(i1-1)+1
             i=i+1
             WRITE(158,'(E13.5)',advance='no') (rho((i3-1)*dfftp%nr1x*dfftp%nr2x+(i2-1)*dfftp%nr1x+i1))
             IF (i == 6 ) THEN
                WRITE(158,'("")')
                i=0
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    CLOSE(158)
    !
    !
#endif
    CALL stop_clock( 'post-processing' )
    RETURN
    !
501 CALL errore ('lr_dump_rho_tot_cube', 'Unable to open file for writing', 1 )
    !-----------------------------------------------------------------------
  END SUBROUTINE lr_dump_rho_tot_cube
  !-----------------------------------------------------------------------
  SUBROUTINE lr_dump_rho_tot_xyzd(rho,identifier)

    ! dump a density file in the x y z density format.
    !-----------------------------------------------------------------------
    USE io_files,             ONLY : prefix
    USE lr_variables,         ONLY : LR_polarization, LR_iteration, cube_save
    USE cell_base
    USE ions_base,            ONLY : nat, ityp, atm, ntyp => nsp, tau
    USE mp,                   ONLY : mp_barrier, mp_sum, mp_bcast, mp_get
    USE mp_world,             ONLY : world_comm
    USE mp_global,            ONLY : me_image, intra_image_comm, me_bgrp, &
                                     nproc_bgrp, intra_bgrp_comm, my_bgrp_id

    USE constants,            ONLY : BOHR_RADIUS_ANGS
    USE fft_base,             ONLY : dfftp !this contains dfftp%npp (number of z planes per processor
    ! and dfftp%ipp (offset of the first z plane of the processor
    !
    IMPLICIT NONE
    !
    real (kind=dp), INTENT(in)   :: rho(:)
    CHARACTER(len=10), INTENT(in) :: identifier
    !
    CHARACTER(len=80) :: filename
    !
    CHARACTER(len=6), EXTERNAL :: int_to_char
    !
    !Local
    INTEGER          :: i, nt, i1, i2, i3, at_num, iopool_id,ldr,kk,ionode_pool
    REAL(DP), ALLOCATABLE :: rho_plane(:)
    INTEGER,  ALLOCATABLE :: kowner(:)
    REAL(DP)              :: tpos(3), inpos(3)
    INTEGER, EXTERNAL:: atomic_number
    !

    CALL start_clock( 'post-processing' )
    IF (lr_verbosity > 5) WRITE(stdout,'("<lr_dump_rho_tot_xyzd>")')
    !

#if defined(__MPI)
    !Derived From Modules/xml_io_base.f90
    ALLOCATE( rho_plane( dfftp%nr1*dfftp%nr2 ) )
    ALLOCATE( kowner( dfftp%nr3 ) )
    IF (ionode) THEN
       filename = trim(prefix) // "-" // identifier // "-pol" //trim(int_to_char(LR_polarization))// ".xyzd"
       WRITE(stdout,'(/5X,"Writing xyzd file for response charge density")')
       !write(stdout, *) filename
       !write(stdout,'(5X,"|rho|=",D15.8)') rho_sum
       OPEN (158, file = filename, form = 'formatted', status = 'replace', err=501)
       WRITE(158,'("#NAT=",I5)') nat
       WRITE(158,'("#NR1=",I5,"at1=",3F12.6)') dfftp%nr1, (alat*at(i,1)/dble(dfftp%nr1),i=1,3)
       WRITE(158,'("#NR2=",I5,"at2=",3F12.6)') dfftp%nr2, (alat*at(i,2)/dble(dfftp%nr2),i=1,3)
       WRITE(158,'("#NR3=",I5,"at3=",3F12.6)') dfftp%nr3, (alat*at(i,3)/dble(dfftp%nr3),i=1,3)

       DO i=1,nat
          ! wrap coordinates back into cell.
          tpos = matmul( transpose(bg), tau(:,i) )
          tpos = tpos - nint(tpos - 0.5d0)
          inpos = alat * matmul( at, tpos )
          WRITE(158,'("#",A3,1X,I5,1X,3F12.6)') &
               atm(ityp(i)), atomic_number(trim(atm(ityp(i)))), inpos
       ENDDO
    ENDIF
    !
    ! ... find the index of the pool that will write rho
    !
    IF ( ionode ) iopool_id = my_bgrp_id
    !
    CALL mp_bcast( iopool_id, ionode_id, intra_image_comm )
    !
    ! ... find the index of the ionode within its own pool
    !
    IF ( ionode ) ionode_pool = me_bgrp
    !
    CALL mp_bcast( ionode_pool, ionode_id, intra_image_comm )
    !
    ! ... find out the owner of each "z" plane
    !
    !
    IF (nproc_bgrp > 1) THEN
       DO i = 1, nproc_bgrp
          !
          kowner( (dfftp%ipp(i)+1):(dfftp%ipp(i)+dfftp%npp(i)) ) = i - 1
          !
       ENDDO
    ELSE
       kowner = ionode_id
    ENDIF
    ldr = dfftp%nr1x*dfftp%nr2x
    !
    !
    ! Each processor is on standby to send its plane to ionode
    !
    DO i3 = 1, dfftp%nr3
       !
       IF( kowner(i3) == me_bgrp ) THEN
          !
          kk = i3
          !
          IF ( nproc_bgrp > 1 ) kk = i3 - dfftp%ipp(me_bgrp+1)
          !
          DO i2 = 1, dfftp%nr2
             !
             DO i1 = 1, dfftp%nr1
                !
                rho_plane(i1+(i2-1)*dfftp%nr1) = rho(i1+(i2-1)&
                     &*dfftp%nr1x+(kk-1)*ldr) 
                !
             ENDDO
             !
          ENDDO
          !
       ENDIF
       !Send plane to ionode
       IF ( kowner(i3) /= ionode_pool .and. my_bgrp_id == iopool_id ) &
            CALL mp_get( rho_plane, rho_plane, me_bgrp, ionode_pool,&
            & kowner(i3), i3, intra_bgrp_comm ) 
       !
       ! write
       IF ( ionode ) THEN
          DO i2 = 1, dfftp%nr2
             !
             DO i1 = 1, dfftp%nr1
                !
                WRITE(158,'(f15.8,3X)', advance='no') (DBLE(i1-1)&
                     &*(alat*BOHR_RADIUS_ANGS*(at(1,1)+at(2,1)+at(3&
                     &,1))/DBLE(dfftp%nr1-1))) 
                WRITE(158,'(f15.8,3X)', advance='no') (DBLE(i2-1)&
                     &*(alat*BOHR_RADIUS_ANGS*(at(1,2)+at(2,2)+at(3&
                     &,2))/DBLE(dfftp%nr2-1))) 
                WRITE(158,'(f15.8,3X)', advance='no') (DBLE(i3-1)&
                     &*(alat*BOHR_RADIUS_ANGS*(at(1,3)+at(2,3)+at(3&
                     &,3))/DBLE(dfftp%nr3-1))) 
                WRITE(158,'(e13.5)') rho_plane((i2-1)*dfftp%nr1+i1)
             ENDDO
          ENDDO
       ENDIF
       !
    ENDDO
    !
    DEALLOCATE( rho_plane )
    DEALLOCATE( kowner )

    IF (ionode) CLOSE(158)
#else
    !
    !

    filename = TRIM(prefix) // "-" // identifier // "-pol" &
         &//TRIM(int_to_char(LR_polarization))// ".xyzd" 
    WRITE(stdout,'(/5X,"Writing xyzd file for response charge density")')
    !write(stdout, *) filename
    !write(stdout,'(5X,"|rho|=",D15.8)') rho_sum
    OPEN (158, file = filename, form = 'formatted', status = 'replace', err=501)

    WRITE(158,*) "# x         y          z        density"
    DO i3=0,(dfftp%nr3-1)
       DO i2=0,(dfftp%nr2-1)
          DO i1=0,(dfftp%nr1-1)
             WRITE(158,'(f15.8,3X)', advance='no') (DBLE(i1)*(alat&
                  &*BOHR_RADIUS_ANGS*(at(1,1)+at(2,1)+at(3,1))& 
                  &/DBLE(dfftp%nr1-1))) 
             WRITE(158,'(f15.8,3X)', advance='no') (DBLE(i2)*(alat&
                  &*BOHR_RADIUS_ANGS*(at(1,2)+at(2,2)+at(3,2))&
                  &/DBLE(dfftp%nr2-1))) 
             WRITE(158,'(f15.8,3X)', advance='no') (DBLE(i3)*(alat&
                  &*BOHR_RADIUS_ANGS*(at(1,3)+at(2,3)+at(3,3))&
                  &/DBLE(dfftp%nr3-1))) 
             WRITE(158,'(e13.5)') &
                  &rho(i3*dfftp%nr1*dfftp%nr2+i2*dfftp%nr1+i1+1)
          ENDDO
       ENDDO
    ENDDO
    CLOSE(158)
    !
    !
#endif
    CALL stop_clock( 'post-processing' )
    RETURN
    !
501 CALL errore ('lr_dump_rho_tot_xyzd', 'Unable to open file for writing', 1 )
    !-----------------------------------------------------------------------
  END SUBROUTINE lr_dump_rho_tot_xyzd
  !-----------------------------------------------------------------------

  SUBROUTINE lr_dump_rho_tot_xcrys(rho, identifier)
    !---------------------------------------------------------------------------
    ! This routine dumps the charge density in xcrysden format, copyright information from
    ! the derived routines follows
    !---------------------------------------------------------------------------
    ! Copyright (C) 2003 Tone Kokalj
    ! This file is distributed under the terms of the
    ! GNU General Public License. See the file `License'
    ! in the root directory of the present distribution,
    ! or http://www.gnu.org/copyleft/gpl.txt .
    !
    ! This file holds XSF (=Xcrysden Structure File) utilities.
    ! Routines written by Tone Kokalj on Mon Jan 27 18:51:17 CET 2003
    !
    ! -------------------------------------------------------------------
    !   this routine writes the crystal structure in XSF format
    ! -------------------------------------------------------------------
    ! -------------------------------------------------------------------
    !   this routine writes the 3D scalar field (i.e. uniform mesh of points)
    !   in XSF format using the FFT mesh (i.e. fast write)
    ! -------------------------------------------------------------------
    USE constants,             ONLY : BOHR_RADIUS_ANGS
    USE io_files,              ONLY : prefix
    USE lr_variables,          ONLY : LR_polarization, LR_iteration, cube_save
    USE cell_base
    USE ions_base,             ONLY : nat, ityp, atm, ntyp => nsp, tau
    USE mp,                   ONLY : mp_barrier, mp_sum, mp_bcast, mp_get
    USE mp_world,             ONLY : world_comm
    USE mp_global,            ONLY : me_image, intra_image_comm, me_bgrp, nproc_bgrp, &
         intra_bgrp_comm, my_bgrp_id

    USE constants,            ONLY : BOHR_RADIUS_ANGS
    USE fft_base,             ONLY : dfftp !this contains dfftp%npp (number of z planes per processor
    ! and dfftp%ipp (offset of the first z plane of the processor



    IMPLICIT NONE
    !
    real (kind=dp), INTENT(in)   :: rho(:)
    CHARACTER(len=10), INTENT(in) :: identifier
    ! INTERNAL
    CHARACTER(len=80) :: filename
    ! --
    INTEGER          :: i, j, n
    INTEGER       :: i1, i2, i3, ix, iy, iz, count, &
         ind_x(10), ind_y(10),ind_z(10)

    real(DP)    :: at1 (3, 3)
    CHARACTER(len=6), EXTERNAL :: int_to_char
    !Local
    INTEGER          :: iopool_id,ldr,kk,ionode_pool,six_count
    REAL(DP), ALLOCATABLE :: rho_plane(:)
    INTEGER,  ALLOCATABLE :: kowner(:)
    !
    six_count=0
    CALL start_clock( 'post-processing' )
    IF (lr_verbosity > 5) WRITE(stdout,'("<lr_dump_rho_tot_xsf>")')
#if defined(__MPI)
    IF (ionode) THEN
       !
       !
       filename = trim(prefix) // "-" // identifier // "-pol" //trim(int_to_char(LR_polarization))// ".xsf"
       WRITE(stdout,'(/5X,"Writing xsf file for response charge density")')
       !write(stdout, *) filename
       OPEN (158, file = filename, form = 'formatted', status = 'replace', err=501)


       ! convert lattice vectors to ANGSTROM units ...
       DO i=1,3
          DO j=1,3
             at1(j,i) = at(j,i)*alat*BOHR_RADIUS_ANGS
          ENDDO
       ENDDO

       WRITE(158,*) 'CRYSTAL'
       WRITE(158,*) 'PRIMVEC'
       WRITE(158,'(2(3F15.9/),3f15.9)') at1
       WRITE(158,*) 'PRIMCOORD'
       WRITE(158,*) nat, 1

       DO n=1,nat
          ! positions are in Angstroms
          WRITE(158,'(a3,3x,3f15.9)') atm(ityp(n)), &
               tau(1,n)*alat*BOHR_RADIUS_ANGS, &
               tau(2,n)*alat*BOHR_RADIUS_ANGS, &
               tau(3,n)*alat*BOHR_RADIUS_ANGS
       ENDDO

       ! --
       ! XSF scalar-field header
       WRITE(158,'(a)') 'BEGIN_BLOCK_DATAGRID_3D'
       WRITE(158,'(a)') '3D_PWSCF'
       WRITE(158,'(a)') 'DATAGRID_3D_UNKNOWN'

       ! number of points in each direction
       WRITE(158,*) dfftp%nr1+1, dfftp%nr2+1, dfftp%nr3+1
       ! origin
       WRITE(158,'(3f10.6)') 0.0d0, 0.0d0, 0.0d0
       ! 1st spanning (=lattice) vector
       WRITE(158,'(3f10.6)') (BOHR_RADIUS_ANGS*alat*at(i,1),i=1,3) ! in ANSTROMS
       ! 2nd spanning (=lattice) vector
       WRITE(158,'(3f10.6)') (BOHR_RADIUS_ANGS*alat*at(i,2),i=1,3)
       ! 3rd spanning (=lattice) vector
       WRITE(158,'(3f10.6)') (BOHR_RADIUS_ANGS*alat*at(i,3),i=1,3)
    ENDIF
    ALLOCATE( rho_plane( dfftp%nr1*dfftp%nr2 ) )
    ALLOCATE( kowner( dfftp%nr3 ) )
    !
    ! ... find the index of the pool that will write rho
    !
    IF ( ionode ) iopool_id = my_bgrp_id
    !
    CALL mp_bcast( iopool_id, ionode_id, intra_image_comm )
    !
    ! ... find the index of the ionode within its own pool
    !
    IF ( ionode ) ionode_pool = me_bgrp
    !
    CALL mp_bcast( ionode_pool, ionode_id, intra_image_comm )
    !
    ! ... find out the owner of each "z" plane
    !
    !
    IF (nproc_bgrp > 1) THEN
       DO i = 1, nproc_bgrp
          !
          kowner( (dfftp%ipp(i)+1):(dfftp%ipp(i)+dfftp%npp(i)) ) = i - 1
          !
       ENDDO
    ELSE
       kowner = ionode_id
    ENDIF
    ldr = dfftp%nr1x*dfftp%nr2x
    !
    !
    ! Each processor is on standby to send its plane to ionode
    !
    DO i3 = 1, dfftp%nr3
       !
       IF( kowner(i3) == me_bgrp ) THEN
          !
          kk = i3
          !
          IF ( nproc_bgrp > 1 ) kk = i3 - dfftp%ipp(me_bgrp+1)
          !
          DO i2 = 1, dfftp%nr2
             !
             DO i1 = 1, dfftp%nr1
                !
                rho_plane(i1+(i2-1)*dfftp%nr1) = rho(i1+(i2-1)*dfftp%nr1x+(kk-1)*ldr)
                !
             ENDDO
             !
          ENDDO
          !
       ENDIF
       !Send plane to ionode
       IF ( kowner(i3) /= ionode_pool .and. my_bgrp_id == iopool_id ) &
            CALL mp_get( rho_plane, rho_plane, &
            me_bgrp, ionode_pool, kowner(i3), i3, intra_bgrp_comm )
       !
       ! write
       IF ( ionode ) THEN
          DO i2 = 1, dfftp%nr2
             !
             DO i1 = 1, dfftp%nr1
                six_count=six_count+1
                WRITE(158,'(e13.5)',advance='no') rho_plane((i2-1)*dfftp%nr1+i1)
                IF (six_count == 6 ) THEN
                   WRITE(158,'("")')
                   six_count=0
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       !
    ENDDO
    !
    DEALLOCATE( rho_plane )
    DEALLOCATE( kowner )

    IF (ionode) CLOSE(158)

#else
    !
    !
    filename = trim(prefix) // "-" // identifier // "-pol" //trim(int_to_char(LR_polarization))// ".xsf"
    WRITE(stdout,'(/5X,"Writing xsf file for response charge density")')
    !write(stdout, *) filename
    OPEN (158, file = filename, form = 'formatted', status = 'replace', err=501)


    ! convert lattice vectors to ANGSTROM units ...
    DO i=1,3
       DO j=1,3
          at1(j,i) = at(j,i)*alat*BOHR_RADIUS_ANGS
       ENDDO
    ENDDO

    WRITE(158,*) 'CRYSTAL'
    WRITE(158,*) 'PRIMVEC'
    WRITE(158,'(2(3F15.9/),3f15.9)') at1
    WRITE(158,*) 'PRIMCOORD'
    WRITE(158,*) nat, 1

    DO n=1,nat
       ! positions are in Angstroms
       WRITE(158,'(a3,3x,3f15.9)') atm(ityp(n)), &
            tau(1,n)*alat*BOHR_RADIUS_ANGS, &
            tau(2,n)*alat*BOHR_RADIUS_ANGS, &
            tau(3,n)*alat*BOHR_RADIUS_ANGS
    ENDDO

    ! --
    ! XSF scalar-field header
    WRITE(158,'(a)') 'BEGIN_BLOCK_DATAGRID_3D'
    WRITE(158,'(a)') '3D_PWSCF'
    WRITE(158,'(a)') 'DATAGRID_3D_UNKNOWN'

    ! number of points in each direction
    WRITE(158,*) dfftp%nr1+1, dfftp%nr2+1, dfftp%nr3+1
    ! origin
    WRITE(158,'(3f10.6)') 0.0d0, 0.0d0, 0.0d0
    ! 1st spanning (=lattice) vector
    WRITE(158,'(3f10.6)') (BOHR_RADIUS_ANGS*alat*at(i,1),i=1,3) ! in ANSTROMS
    ! 2nd spanning (=lattice) vector
    WRITE(158,'(3f10.6)') (BOHR_RADIUS_ANGS*alat*at(i,2),i=1,3)
    ! 3rd spanning (=lattice) vector
    WRITE(158,'(3f10.6)') (BOHR_RADIUS_ANGS*alat*at(i,3),i=1,3)

    count=0
    DO i3=0,dfftp%nr3
       iz = mod(i3,dfftp%nr3)
       !iz = mod(i3,dfftp%nr3) + 1

       DO i2=0,dfftp%nr2
          iy = mod(i2,dfftp%nr2)
          !iy = mod(i2,dfftp%nr2) + 1

          DO i1=0,dfftp%nr1
             ix = mod(i1,dfftp%nr1)
             !ix = mod(i1,dfftp%nr1) + 1

             !ii = (1+ix) + iy*dfftp%nr1x + iz*dfftp%nr1x*dfftp%nr2x
             IF (count<6) THEN
                count = count + 1
                !ind(count) = ii
             ELSE
                WRITE(158,'(6e13.5)') &
                     (rho(ind_x(i)+1+dfftp%nr1*ind_y(i)+dfftp%nr1*dfftp%nr2*ind_z(i)),i=1,6)
                count=1
                !ind(count) = ii
             ENDIF
             ind_x(count) = ix
             ind_y(count) = iy
             ind_z(count) = iz
          ENDDO
       ENDDO
    ENDDO
    WRITE(158,'(6e13.5:)') (rho(ind_x(i)+1+dfftp%nr1*ind_y(i)+dfftp%nr1*dfftp%nr2*ind_z(i)),i=1,count)
    WRITE(158,'(a)') 'END_DATAGRID_3D'
    WRITE(158,'(a)') 'END_BLOCK_DATAGRID_3D'
#endif
    CALL stop_clock( 'post-processing' )
    RETURN
    !
501 CALL errore ('lr_dump_rho_tot_xyzd', 'Unable to open file for writing', 1 )
    !-----------------------------------------------------------------------
  END SUBROUTINE lr_dump_rho_tot_xcrys
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE lr_dump_rho_tot_pxyd(rho,identifier)
    !-----------------------------------------------------------------------
    ! dump a density file in the x y plane density summed over z planes format.
    !-----------------------------------------------------------------------
    USE io_files,              ONLY : prefix
    USE lr_variables,          ONLY : LR_polarization, LR_iteration, cube_save
    USE fft_base,              ONLY : dfftp
    USE cell_base
    USE ions_base,             ONLY : nat, ityp, atm, ntyp => nsp, tau
    USE mp,                    ONLY : mp_barrier, mp_sum
    USE mp_world,              ONLY : world_comm
    USE mp_global,             ONLY : intra_bgrp_comm
    USE constants,             ONLY : BOHR_RADIUS_ANGS
    !
    IMPLICIT NONE
    !
    real (kind=dp), INTENT(in)   :: rho(:)
    CHARACTER(len=10), INTENT(in) :: identifier
    !
    CHARACTER(len=80) :: filename
    !
    CHARACTER(len=6), EXTERNAL :: int_to_char
    !
    !Local
    INTEGER          :: i, nt, i1, i2, i3, at_num
    INTEGER, EXTERNAL:: atomic_number
    real(DP)    :: at_chrg, tpos(3), inpos(3),rho_sum
    !

    CALL start_clock( 'post-processing' )
    IF (lr_verbosity > 5) WRITE(stdout,'("<lr_dump_rho_tot_pxyd>")')
    rho_sum=0.0d0
    DO i=1,dfftp%nnr
       rho_sum=rho_sum+rho(i)
    ENDDO
    !

#if defined(__MPI)
    IF (ionode) THEN
#endif
       !
       !

       filename = trim(prefix) // "-" // identifier // "-pol" //trim(int_to_char(LR_polarization))// ".pxyd"
       WRITE(stdout,'(/5X,"Writing z plane averaged pxyd file for response charge density")')
       !write(stdout, *) filename
       WRITE(stdout,'(5X,"|rho|=",D15.8)') rho_sum
       OPEN (158, file = filename, form = 'formatted', status = 'replace', err=501)

       WRITE(158,*) "# x         y          z        density"
       DO i1=0,(dfftp%nr1-1)
          DO i2=0,(dfftp%nr2-1)
             rho_sum=0
             DO i3=0,(dfftp%nr3-1)
                rho_sum=rho_sum+rho(i3*dfftp%nr1*dfftp%nr2+i2*dfftp%nr1+i1+1)
             ENDDO
             WRITE(158,'(f15.8,3X)', advance='no') (dble(i1)*(alat*BOHR_RADIUS_ANGS*(at(1,1)+at(2,1)+at(3,1))/dble(dfftp%nr1-1)))
             WRITE(158,'(f15.8,3X)', advance='no') (dble(i2)*(alat*BOHR_RADIUS_ANGS*(at(1,2)+at(2,2)+at(3,2))/dble(dfftp%nr2-1)))
             WRITE(158,'(e13.5)') rho_sum
          ENDDO
       ENDDO
       CLOSE(158)
       !
       !
#if defined(__MPI)
    ENDIF
    CALL mp_barrier(world_comm)
#endif
    CALL stop_clock( 'post-processing' )
    RETURN
    !
501 CALL errore ('lr_dump_rho_tot_pxyd', 'Unable to open file for writing', 1 )
    !-----------------------------------------------------------------------
  END SUBROUTINE lr_dump_rho_tot_pxyd

END MODULE charg_resp
