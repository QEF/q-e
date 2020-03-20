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
                                         charge_response, size_evc, LR_polarization, LR_iteration, &
                                         magnons, evc1_rgt, evc1_lft, evc1_rgt_new, evc1_lft_new, &
                                         alpha_magnons_store, gamma_magnons_store, evc1_rgt_old, &
                                         evc1_lft_old, O_psi, n_op
    USE uspp,                     ONLY : vkb, nkb, okvan
    USE wvfct,                    ONLY : nbnd, npwx
    USE control_flags,            ONLY : gamma_only
    USE becmod,                   ONLY : bec_type, becp, calbec
    USE realus,                   ONLY : invfft_orbital_gamma, initialisation_level,    &
                                         fwfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                         v_loc_psir, s_psir_gamma
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
    REAL(kind=dp) :: alpha, beta, gamma
    COMPLEX(kind=dp) :: alphac, gammac
    COMPLEX(kind=dp) ,allocatable :: zeta(:)
    INTEGER(kind=c_int) :: kilobytes
    !
    IF (lr_verbosity > 5) THEN
       WRITE(stdout,'("<lr_lanczos_one_step>")')
    ENDIF
    !
    IF (magnons) THEN 
       ALLOCATE (zeta(n_op))
    ELSE 
       ALLOCATE (zeta(n_ipol))
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
    ELSEIF (magnons) THEN
       !
       CALL lr_apply_liouvillian_magnons(evc1_rgt(:,:,:,:), evc1_rgt_new(:,:,:,:),.false.)
       ! 
       IF (.not. pseudo_hermitian) &
       CALL lr_apply_liouvillian_magnons(evc1_lft(:,:,:,:), evc1_lft_new(:,:,:,:),.true.)
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
    ! Apply S operator if USPP, otherwise just copy 
    ! one array into another.
    !
    IF (pseudo_hermitian) THEN
       CALL lr_apply_s(evc1_new(:,:,:,1), sevc1_new(:,:,:))
    ELSE
       CALL lr_apply_s(evc1(:,:,:,2), sevc1(:,:,:))
    ENDIF
    !
    ! call general lanczos iteration routines
    ! O. Baseggio (2019)
    !
    IF (pseudo_hermitian) THEN
       IF (eels) THEN
          CALL lanczos_pseudohermitian(LR_iteration,size(evc1,1), size(evc1,2), size(evc1,3),&
                                      &evc1(:,:,:,1), evc1_new(:,:,:,1), sevc1_new(:,:,:), &
                                      &evc1_old(:,:,:,1), n_ipol, d0psi2, alpha, beta, &
                                      &gamma, zeta)
       ELSEIF (magnons) THEN
          call lanczos_pseudohermitian_c(LR_iteration,size(evc1_rgt,1), size(evc1_rgt,2), &
                                        &size(evc1_rgt,3), evc1_rgt, evc1_rgt_new, evc1_rgt_old, &
                                        &n_op, O_psi, alphac, beta, gammac, zeta)          
       ELSE
          CALL lanczos_pseudohermitian(LR_iteration,size(evc1,1), size(evc1,2), size(evc1,3),&
                                      &evc1(:,:,:,1), evc1_new(:,:,:,1), sevc1_new(:,:,:), &
                                      &evc1_old(:,:,:,1), n_ipol, d0psi(:,:,:,:), alpha, beta, &
                                      &gamma, zeta)
       ENDIF
    ELSE
       IF (eels) THEN
          CALL lanczos_nonhermitian(LR_iteration,size(evc1,1), size(evc1,2), size(evc1,3),&
                                   &evc1(:,:,:,:), evc1_new(:,:,:,:), sevc1(:,:,:), &
                                   &evc1_old(:,:,:,1), n_ipol, d0psi2, alpha, beta, &
                                   &gamma, zeta)
       ELSEIF (magnons) THEN
          call lanczos_nonhermitian_c(LR_iteration,size(evc1_rgt,1), size(evc1_rgt,2), & 
                                     &size(evc1_rgt,3), evc1_rgt, evc1_rgt_new, evc1_rgt_old, & 
                                     &evc1_lft, evc1_lft_new, evc1_lft_old, n_op, O_psi, alphac, &
                                     &beta, gammac, zeta)
       ELSE
          CALL lanczos_nonhermitian(LR_iteration,size(evc1,1), size(evc1,2), size(evc1,3),&
                                   &evc1(:,:,:,1), evc1_new(:,:,:,1), sevc1(:,:,:), &
                                   &evc1_old(:,:,:,1), n_ipol, d0psi, alpha, beta, &
                                   &gamma, zeta)
       ENDIF
    ENDIF 
    !
    ! beta<0 is a serious error for the pseudo-Hermitian algorithm
    !
    IF (magnons) THEN
       alpha_magnons_store(pol_index,LR_iteration) = alphac
       WRITE(stdout,'(5X,"alpha(",i8.8,")=",2(f21.14))') LR_iteration,dble(alphac), aimag(alphac)
    ELSE
       alpha_store(pol_index,LR_iteration) = alpha
       WRITE(stdout,'(5X,"alpha(",i8.8,")=",f10.6)') LR_iteration, alpha
    ENDIF
    !
    IF ( beta < 0.0d0 .and. pseudo_hermitian) CALL error_beta()
    !
    IF ( abs(beta)<1.0d-12 ) THEN
       !
       WRITE(stdout,'(5x,"lr_lanczos: Left and right Lanczos vectors are orthogonal, &
                      & this is a violation of oblique projection")')
       !
    ENDIF
    !
    IF (magnons) THEN
       beta_store (pol_index,LR_iteration) = beta
       WRITE(stdout,'(5X,"beta (",i8.8,")=",f21.14)') LR_iteration, beta
       gamma_magnons_store(pol_index,LR_iteration) = gammac
       WRITE(stdout,'(5X,"gamma(",i8.8,")=",2(f21.14))') LR_iteration, dble(gammac), aimag(gammac) 
    ELSE
        beta_store (pol_index,LR_iteration) = beta
       WRITE(stdout,'(5X,"beta (",i8.8,")=",f10.6)') LR_iteration, beta
       gamma_store(pol_index,LR_iteration) = gamma
       WRITE(stdout,'(5X,"gamma(",i8.8,")=",f10.6)') LR_iteration, gamma
    ENDIF
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
    IF (magnons) THEN
       !
       DO ip = 1, n_op
          zeta_store (pol_index,ip,LR_iteration) = zeta(ip)
          WRITE(stdout,'(5x,"z1= ",1x,i6,2(1x,e22.15))') ip,real(zeta(ip)),aimag(zeta(ip))
       ENDDO
       !
    ELSEIF (mod(LR_iteration,2)==0) THEN
       !
       DO ip = 1, n_ipol
          !
          zeta_store (pol_index,ip,LR_iteration) = zeta(ip)
          WRITE(stdout,'(5x,"z1= ",1x,i6,2(1x,e22.15))') ip,real(zeta(ip)),aimag(zeta(ip))
          !
       ENDDO
       !
       ! OBM: evc1(:,:,:,1) contains the q of x for even steps, 
       ! lets calculate the response related observables.
       ! 
       IF (charge_response == 1 .and. .not.eels .and. .not. magnons) THEN
          CALL lr_calc_dens(evc1_old(:,:,:,1), .true.)
          CALL lr_calc_F(evc1_old(:,:,:,1))
       ENDIF
       !
    ELSE
       !
       DO ip = 1, n_ipol
          !
          zeta_store (pol_index,ip,LR_iteration) = zeta(ip)
          WRITE(stdout,'(5x,"z1= ",1x,i6,2(1x,e22.15))') ip,real(zeta(ip)),aimag(zeta(ip))
          !
       ENDDO
       !
    ENDIF
    !
    ! X. Ge: To increase the stability, apply lr_ortho.
    ! I.Timrov: Actually, without this trick, it turns out that 
    ! the Lanczos chain is not stable when pseudo_hermitian=.false., 
    ! because there is a warning from lr_calc_dens that the integral of 
    ! the response charge density does not some to zero, i.e. the non-zero 
    ! value of the integral increases during the Lanczos chain.
    !
    IF (.not.eels .and. .not. magnons) THEN
       !
       DO ik=1, nks
          CALL lr_ortho(evc1(:,:,ik,1), evc0(:,:,ik), ik, ik, sevc0(:,:,ik),.true.)
          IF (.not. pseudo_hermitian) &
          CALL lr_ortho(evc1(:,:,ik,2), evc0(:,:,ik), ik, ik, sevc0(:,:,ik),.true.)
       ENDDO
       !
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
    DEALLOCATE(zeta)
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
