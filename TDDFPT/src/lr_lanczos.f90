!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE lr_lanczos

CONTAINS

!-------------------------------------------------------------------------------------------!
!
! There are 2 flavors of the Lanczos algorithm implemented in TDDFPT:
! - non-Hermitian Lanczos biorthogonalization algorithm
! - pseudo-Hermitian Lanczos algorithm (twice faster)
! For details see:
! - O. Malcioglu, R. Gebauer, D. Rocca, S. Baroni,
!   Comput. Phys. Commun. 182, 1744 (2011). 
! - X. Ge, S. Binnie, D. Rocca, R. Gebauer, S. Baroni, 
!   Comput. Phys. Commun. 185, 2080 (2014). 
!
! Non-Hermitian algorithm:
! This subroutine handles two interleaved non-Hermitian chains for x and y at the same time,
! in the beginning
! for even steps evc1(:,:,:,1) corresponds to q of y and evc1(:,:,:,2) corresponds to p of x
! for odd steps evc1(:,:,:,1) corresponds to q of x and evc1(:,:,:,2) corresponds to p of y
!
! using the terminology of Lanczos chains by Y.Saad, this is effectively equal to
! altering the A matrix correspondingly for even and odd steps correspondingly.
! This change is controlled by the interaction parameter in lr_apply_liouvillian
!
! For further reference please refer to eq. (32) and (33) in
! Ralph Gebauer, Brent Walker J. Chem. Phys., 127, 164106 (2007)
!
! Modified by Osman Baris Malcioglu (2009)
! Modified by Xiaochuan Ge to introduce pseudo-Hermitian algorithm (2013)
! Modified by Iurii Timrov to introduce EELS (2013)
!
!--------------------------------------------------------------------------------------------!
SUBROUTINE one_lanczos_step()
    
    USE io_global,                ONLY : ionode, stdout
    USE kinds,                    ONLY : dp
    USE klist,                    ONLY : nks,xk
    USE lr_variables,             ONLY : n_ipol, ltammd, pseudo_hermitian, itermax,  &
                                         evc1, evc1_new, evc1_old, sevc1_new, sevc1, &
                                         evc0, sevc0, d0psi, d0psi2, eels, test_case_no, lr_verbosity, &
                                         alpha_store, beta_store, gamma_store, zeta_store,     &
                                         charge_response, size_evc, LR_polarization, LR_iteration
    USE uspp,                     ONLY : vkb, nkb, okvan
    USE wvfct,                    ONLY : nbnd, npwx
    USE control_flags,            ONLY : gamma_only, tqr
    USE becmod,                   ONLY : bec_type, becp, calbec
    USE realus,                   ONLY : real_space, invfft_orbital_gamma, initialisation_level,    &
                                         fwfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                         v_loc_psir, s_psir_gamma, real_space_debug
    USE charg_resp,               ONLY : w_T_beta_store, w_T, lr_calc_F
    USE lr_us,                    ONLY : lr_apply_s
    USE noncollin_module,         ONLY : npol
    ! Debugging
    USE iso_c_binding,            ONLY : c_int

    USE qpoint,                   ONLY : nksq
    !
    IMPLICIT NONE
    !
    COMPLEX(kind=dp),EXTERNAL :: lr_dot
    INTEGER :: ik, ip, ibnd,ig, pol_index
    CHARACTER(len=6), EXTERNAL :: int_to_char
    !
    ! Local variables
    !
    REAL(kind=dp) :: alpha, beta, gamma, angle
    COMPLEX(kind=dp) :: zeta
    INTEGER(kind=c_int) :: kilobytes
    !
    IF (lr_verbosity > 5) THEN
       WRITE(stdout,'("<lr_lanczos_one_step>")')
    ENDIF
    !
    IF (lr_verbosity > 10) THEN
       PRINT *, "Real space = ", real_space
       PRINT *, "Real space debug ", real_space_debug
       PRINT *, "TQR = ", tqr
    ENDIF
    !
    ! Memory usage
    !
    !IF (eels) THEN
    !   CALL memstat( kilobytes )
    !   IF ( kilobytes > 0 ) WRITE(stdout,'(5X,"lr_lanczos, & 
    !                      & per-process dynamical memory:",f7.1,"Mb")' ) kilobytes/1000.0
    !ENDIF
    !
    CALL start_clock('one_step')
    !
    pol_index = 1
    !
    IF ( n_ipol /= 1 ) pol_index = LR_polarization
    !
    ! Application of the Liouvillian superoperator
    !
    IF (eels) THEN
       !
       ! EELS case
       !
       IF (mod(LR_iteration,2)==0) THEN
          !
          CALL lr_apply_liouvillian_eels(evc1(1,1,1,1),evc1_new(1,1,1,1),.true.)
          IF (.NOT.pseudo_hermitian) &
          CALL lr_apply_liouvillian_eels(evc1(1,1,1,2),evc1_new(1,1,1,2),.false.)
          !
       ELSE
          !
          CALL lr_apply_liouvillian_eels(evc1(1,1,1,1),evc1_new(1,1,1,1),.false.)
          IF (.not. pseudo_hermitian) &
          CALL lr_apply_liouvillian_eels(evc1(1,1,1,2),evc1_new(1,1,1,2),.true.)
          !
       ENDIF
       !
    ELSE
       !
       ! Optical case
       !
       IF (.NOT.ltammd) THEN
          !
          ! Normal regime
          !
          IF (mod(LR_iteration,2)==0) THEN
             !
             CALL lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),.true.)
             IF (.NOT.pseudo_hermitian) &
             CALL lr_apply_liouvillian(evc1(1,1,1,2),evc1_new(1,1,1,2),.false.)
             !
          ELSE
             !
             CALL lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),.false.)
             IF (.not. pseudo_hermitian) &
             CALL lr_apply_liouvillian(evc1(1,1,1,2),evc1_new(1,1,1,2),.true.)
             !
          ENDIF
          !
       ELSE
          !
          ! Tamm-Dancoff approximation
          !
          CALL lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),.true.)
          !
          IF (.not. pseudo_hermitian) THEN
             CALL zcopy(size_evc,evc1(1,1,1,1),1,evc1(1,1,1,2),1) 
             CALL zcopy(size_evc,evc1_new(1,1,1,1),1,evc1_new(1,1,1,2),1) 
          ENDIF
          !
       ENDIF
       !
    ENDIF 
    !
    ! By construction <p|Lq>=0 should be 0, forcing this both conserves 
    ! resources and increases stability.
    !
    alpha = 0.0d0
    alpha_store(pol_index,LR_iteration) = alpha
    WRITE(stdout,'(5X,"alpha(",i8.8,")=",f10.6)') LR_iteration, alpha
    !
    ! Apply S operator if USPP, otherwise just copy 
    ! one array into another.
    !
    IF (pseudo_hermitian) THEN
       CALL lr_apply_s(evc1_new(:,:,:,1), sevc1_new(:,:,:))
    ELSE
       CALL lr_apply_s(evc1(:,:,:,2), sevc1(:,:,:))
    ENDIF
    !
    ! Orthogonality requirement: <v|\bar{L}|v> = 1
    !
    IF (pseudo_hermitian) THEN
       beta = dble(lr_dot(evc1(:,:,:,1), sevc1_new(:,:,:)))
    ELSE
       beta = dble(lr_dot(evc1(:,:,:,1), sevc1(:,:,:)))
    ENDIF
    !
    ! beta<0 is a serious error for the pseudo-Hermitian algorithm
    !
    IF ( beta < 0.0d0 .and. pseudo_hermitian) CALL error_beta()
    !
    IF ( abs(beta)<1.0d-12 ) THEN
       !
       WRITE(stdout,'(5x,"lr_lanczos: Left and right Lanczos vectors are orthogonal, &
                      & this is a violation of oblique projection")')
       !
    ELSEIF ( beta<0.0d0 ) THEN
       !
       beta = sqrt(-beta)
       gamma = -beta
       !
    ELSEIF ( beta>0.0d0 ) THEN 
       !
       ! X. Ge: Actually, this is the only case in the pseudo-Hermitian algorithm.
       !
       beta = sqrt(beta)
       gamma = beta
       !
    ENDIF
    !
    beta_store (pol_index,LR_iteration) = beta
    gamma_store(pol_index,LR_iteration) = gamma
    WRITE(stdout,'(5X,"beta (",i8.8,")=",f10.6)') LR_iteration, beta
    WRITE(stdout,'(5X,"gamma(",i8.8,")=",f10.6)') LR_iteration, gamma
    !
    !IF (LR_iteration==1) WRITE(stdout,'(5X,"Norm of initial Lanczos vectors=",1x,f21.15)') beta
    !
    ! Analysis of the evolution of the beta coefficients.
    !
    IF (lr_verbosity > 3) THEN
       !
       IF ( LR_iteration > 6 ) THEN
          !
          WRITE(stdout,'(5X,"(oscillatory variation) : ",f6.2,"%")') &
              & abs(100*(beta_store(pol_index,LR_iteration)- &
              & (beta_store(pol_index,LR_iteration-6)+beta_store(pol_index,LR_iteration-4)+&
              & beta_store(pol_index,LR_iteration-2))/3.0)/beta_store(pol_index,LR_iteration))
          !  
          WRITE(stdout,'(5X,"(linear variation) : ",f6.2,"%")') &
              & abs(100*(beta_store(pol_index,LR_iteration)- &
              & (beta_store(pol_index,LR_iteration-3)+beta_store(pol_index,LR_iteration-2)+&
              & beta_store(pol_index,LR_iteration-1))/3.0)/beta_store(pol_index,LR_iteration))
          !
       ENDIF
       !
    ENDIF
    !
    ! Renormalize q(i) and Lq(i), also p(i) and Lp(i) in the non-Hermitian case
    !
    CALL zscal(size_evc,cmplx(1.0d0/beta,0.0d0,kind=dp),evc1(1,1,1,1),1)
    CALL zscal(size_evc,cmplx(1.0d0/beta,0.0d0,kind=dp),evc1_new(1,1,1,1),1)
    !
    IF (.not.pseudo_hermitian) THEN
       CALL zscal(size_evc,cmplx(1.0d0/gamma,0.0d0,kind=dp),evc1(1,1,1,2),1)
       CALL zscal(size_evc,cmplx(1.0d0/gamma,0.0d0,kind=dp),evc1_new(1,1,1,2),1)
    ENDIF
    !
    ! Calculation of zeta coefficients.
    ! See Eq.(35) in Malcioglu et al., Comput. Phys. Commun. 182, 1744 (2011).
    !
    IF (mod(LR_iteration,2)==0) THEN
       !
       DO ip = 1, n_ipol
          !
          ! Optics: In the ultrasoft case, the S operator was already
          ! applied to d0psi, so we have <S*d0psi|evc1>.
          !
          IF (eels) THEN
             zeta = lr_dot(d0psi2(:,:,:,ip),evc1(:,:,:,1))
          ELSE
             zeta = lr_dot(d0psi(:,:,:,ip),evc1(:,:,:,1))
          ENDIF
          !
          zeta_store (pol_index,ip,LR_iteration) = zeta
          WRITE(stdout,'(5x,"z1= ",1x,i6,2(1x,e22.15))') ip,real(zeta),aimag(zeta)
          !
       ENDDO
       !
       ! OBM: evc1(:,:,:,1) contains the q of x for even steps, 
       ! lets calculate the response related observables.
       ! 
       IF (charge_response == 1 .and. .not.eels) THEN
          CALL lr_calc_dens(evc1(:,:,:,1), .true.)
          CALL lr_calc_F(evc1(:,:,:,1))
       ENDIF
       !
    ELSE
       !
       DO ip = 1, n_ipol
          !
          zeta = (0.0d0,0.0d0)
          zeta_store (pol_index,ip,LR_iteration) = zeta
          WRITE(stdout,'(5x,"z1= ",1x,i6,2(1x,e22.15))') ip,real(zeta),aimag(zeta)
          !
       ENDDO
       !
    ENDIF
    !
    ! X. Ge: q(i+1) = Lq(i) - beta(i)*q(i-1); 
    ! Renormalization will be done in the begining of the next iteration.
    ! In the non-Hermitian case, similar operation needs to be done also for p(i).
    !
    CALL zaxpy(size_evc,-cmplx(gamma,0.0d0,kind=dp),evc1_old(1,1,1,1),1,evc1_new(1,1,1,1),1)
    IF (.not. pseudo_hermitian) &
     CALL zaxpy(size_evc,-cmplx(beta,0.0d0,kind=dp),evc1_old(1,1,1,2),1,evc1_new(1,1,1,2),1)
    !
    ! X. Ge: To increase the stability, apply lr_ortho.
    ! I.Timrov: Actually, without this trick, it turns out that 
    ! the Lanczos chain is not stable when pseudo_hermitian=.false., 
    ! because there is a warning from lr_calc_dens that the integral of 
    ! the response charge density does not some to zero, i.e. the non-zero 
    ! value of the integral increases during the Lanczos chain.
    !
    IF (.not.eels) THEN
       !
       DO ik=1, nks
          CALL lr_ortho(evc1_new(:,:,ik,1), evc0(:,:,ik), ik, ik, sevc0(:,:,ik),.true.)
          IF (.not. pseudo_hermitian) &
          CALL lr_ortho(evc1_new(:,:,ik,2), evc0(:,:,ik), ik, ik, sevc0(:,:,ik),.true.)
       ENDDO
       !
    ENDIF
    !
    ! X. Ge: Throw away q(i-1), and make q(i+1) to be the current vector,
    ! be ready for the next iteration. evc1_new will be free again after this step
    !
    CALL zcopy(size_evc,evc1(1,1,1,1),1,evc1_old(1,1,1,1),1)    ! evc1_old = evc1
    CALL zcopy(size_evc,evc1_new(1,1,1,1),1,evc1(1,1,1,1),1)    ! evc1 = evc1_new
    !
    IF (.not.pseudo_hermitian) THEN
       CALL zcopy(size_evc,evc1(1,1,1,2),1,evc1_old(1,1,1,2),1) ! evc1_old = evc1
       CALL zcopy(size_evc,evc1_new(1,1,1,2),1,evc1(1,1,1,2),1) ! evc1 = evc1_new
    ENDIF
    !
    IF (ionode) THEN
       IF ( charge_response == 1 .and. lr_verbosity > 0) THEN
           WRITE (stdout,'(5x,"(calc=",e22.15," read=",e22.15,")")') &
                                & beta, w_T_beta_store(LR_iteration)
           WRITE (stdout,'(5x,"Weight for this step=",2(e12.5,1x))') w_T(LR_iteration)
       ENDIF
    ENDIF
    !
    CALL stop_clock('one_step')
    !
    RETURN
    !
CONTAINS
    !
 SUBROUTINE error_beta()
    !
    ! IT: This subroutine is called when beta<0 in the pseudo-Hermitian
    ! algorithm. This subroutine prints the response orbitals to
    ! the file for the analysis of the error. Note, in the parallel
    ! execution only a part of the information about the response
    ! orbitals will be written to the file. 
    !
    ! Written by X. Ge (2013)
    !
    IMPLICIT NONE
    !
    INTEGER :: ibnd, ig
    !
    WRITE(stdout,'(5x,"Error: Beta is negative:",5x,e22.15)') beta
    !
    OPEN(301, file="evc1_1.dat")  
    DO ibnd = 1, nbnd
       DO ig = 1, npwx
          WRITE(301,*) ibnd, ig, dble(evc1(ig,ibnd,1,1)), aimag(evc1(ig,ibnd,1,1))
       ENDDO
    ENDDO
    CLOSE(301)
    !
    OPEN(302, file="evc1_2.dat")
    DO ibnd = 1, nbnd
       DO ig = 1, npwx
          WRITE(302,*) ibnd, ig, dble(evc1(ig,ibnd,1,2)), aimag(evc1(ig,ibnd,1,2))
       ENDDO
    ENDDO
    CLOSE(302)
    !
    WRITE (stdout,'(/5x,"!!!WARNING!!! The pseudo-Hermitian Lanczos algorithm is stopping...")') 
    WRITE (stdout,'(/5x,"Try to use the non-Hermitian Lanczos algorithm.")') 
    !
    itermax = LR_iteration-1
    CALL clean_pw( .FALSE. )
    CALL stop_lr( .TRUE. )
    !
 END SUBROUTINE error_beta

END SUBROUTINE one_lanczos_step

END MODULE lr_lanczos
