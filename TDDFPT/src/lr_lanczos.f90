MODULE lr_lanczos
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OBM :
!
! This subroutine handles two interleaved non-hermitian chains for x and y at the same time,
!
! in the beginning
! for odd steps evc1(:,:,:,1) corresponds to q of y and evc1(:,:,:,2) corresponds to p of x
! for even steps evc1(:,:,:,1) corresponds to q of x and evc1(:,:,:,2) corresponds to p of y
!
! using the terminology of Lanczos chains by Y.Saad, this is effectively equal to
! altering the A matrix correspondingly for even and odd steps correspondingly.
! This change is controlled by the interaction parameter in lr_apply_liouvillian
!
! For further reference please refer to eq. (32) and (33) in
! Ralph Gebauer, Brent Walker J. Chem. Phys., 127, 164106 (2007)
!
! Modified by Osman Baris Malcioglu (2009)
  SUBROUTINE one_lanczos_step()
    !
    !   Non-Hermitian Lanczos
    !
    USE io_global,                ONLY : ionode, stdout
    USE kinds,                ONLY : dp
    USE klist,                ONLY : nks,xk
    USE lr_variables,         ONLY : n_ipol, ltammd, itermax,&
                                     evc1, evc1_new, sevc1_new, evc1_old, &
                                     evc0, sevc0, d0psi, &
                                     alpha_store, beta_store, gamma_store, zeta_store,&
                                    charge_response, size_evc, LR_polarization, LR_iteration,&!,real_space
                                     test_case_no
    USE uspp,                 ONLY : vkb, nkb, okvan
    USE wvfct,                ONLY : nbnd, npwx, npw
    USE control_flags,         ONLY : gamma_only,tqr
    USE becmod,               ONLY :  bec_type, becp, calbec
    !use real_beta,           only : ccalbecr_gamma,s_psir,fft_orbital_gamma,bfft_orbital_gamma
    USE realus,               ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                                    bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                    v_loc_psir, s_psir_gamma, igk_k,npw_k, real_space_debug
    USE lr_variables,   ONLY : lr_verbosity
    USE charg_resp,               ONLY : w_T_beta_store,w_T,lr_calc_F
    !Debugging
    USE lr_variables, ONLY: check_all_bands_gamma, check_density_gamma,check_vector_gamma
    !
    IMPLICIT NONE
    !
    COMPLEX(kind=dp),EXTERNAL :: lr_dot
    !
    !integer,intent(in) :: iter, pol
    INTEGER :: ik, ip, ibnd,ig, pol_index
    !
    CHARACTER(len=6), EXTERNAL :: int_to_char
    !
    !   Local variables
    !
    real(kind=dp) :: alpha, beta, gamma, temp
    !
    COMPLEX(kind=dp) :: zeta
    !
    !integer :: n
    !

    !
    IF (lr_verbosity > 5) THEN
      WRITE(stdout,'("<lr_lanczos_one_step>")')
    ENDIF
    IF (lr_verbosity > 10) THEN
      PRINT *, "Real space = ", real_space
      PRINT *, "Real space debug ", real_space_debug
      PRINT *, "TQR = ", tqr
    ENDIF
    !
    CALL start_clock('one_step')
    pol_index=1
    IF ( n_ipol /= 1 ) pol_index=LR_polarization
    !
    ! Calculation of zeta coefficients
    !

    IF (mod(LR_iteration,2)==0) THEN
       !
       DO ip=1,n_ipol
          !
          zeta = lr_dot(d0psi(:,:,:,ip),evc1(:,:,:,1)) !Why gamma point dot?
          zeta_store (pol_index,ip,LR_iteration) = zeta
          WRITE(stdout,'(5x,"z1= ",1x,i6,2(1x,e21.15))') ip,real(zeta),aimag(zeta)
          !
       ENDDO
       !evc1(:,:,:,1) contains the q of x for even steps, lets calculate the response related observables
       !
       IF (charge_response == 1) THEN
        CALL lr_calc_dens(evc1(:,:,:,1), .true.)
        CALL lr_calc_F(evc1(:,:,:,1))
       ENDIF
       !
       !
    ELSE
       !
       DO ip=1,n_ipol
          !
          zeta = (0.0d0,0.0d0)
          zeta_store (pol_index,ip,LR_iteration) = zeta
          WRITE(stdout,'(5x,"z1= ",1x,i6,2(1x,e21.15))') ip,real(zeta),aimag(zeta)
          !
       ENDDO
       !
    ENDIF
    !
    ! Application of the Liouvillian superoperator for the first iteration
    !
    IF(LR_iteration==1) THEN
       LR_iteration=0 !Charge response dump related part is disabled by setting this to zero
       IF(.not.ltammd) THEN
          CALL lr_apply_liouvillian(evc1(:,:,:,1),evc1_new(:,:,:,1),sevc1_new(:,:,:,1),.false.)
          CALL lr_apply_liouvillian(evc1(:,:,:,2),evc1_new(:,:,:,2),sevc1_new(:,:,:,2),.true.)
       ELSE
          CALL lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),sevc1_new(1,1,1,1),.true.)
          evc1(:,:,:,2)=evc1(:,:,:,1)
          evc1_new(:,:,:,2)=evc1_new(:,:,:,1)
          sevc1_new(:,:,:,2)=sevc1_new(:,:,:,1)
       ENDIF
       LR_iteration=1
    ENDIF
    !
    ! The lanczos algorithm starts here
    ! http://www.cs.utk.edu/~dongarra/etemplates/node245.html
    !
    ! Left and right vectors are orthogonalised wrto ground state wf

    !OBM: Notice that here "orthogonalization" is not strictly the true word, as the norm of the vectors change
    !This is due to how the uspp scheme is implemented, the beta are evc1_new(left).sevc1_new(right), that is,
    !a mixing of two vectors, thus the resultant vector from belov should be devoid from S, which affects the norm
    !the modification in lr_ortho subroutine handles this reversal (the last flag).
    DO ik=1, nks
     CALL lr_ortho(evc1_new(:,:,ik,1), evc0(:,:,ik), ik, ik, sevc0(:,:,ik),.true.)
     CALL lr_ortho(evc1_new(:,:,ik,2), evc0(:,:,ik), ik, ik, sevc0(:,:,ik),.true.)
    ENDDO

    !
    ! By construction <p|Lq>=0 should be 0, forcing this both conserves resources and increases stability
    alpha=0.0d0
    alpha_store(pol_index,LR_iteration)   = alpha
    WRITE(stdout,'(5X,"alpha(",i8.8,")=",e21.15)') LR_iteration,alpha
    !
    IF ( gamma_only ) THEN
       !
       IF ( nkb >0 ) THEN
        !
        IF (real_space_debug>5) THEN
        ! real space & nkb > 0
         !
          !
          ! The following part converts evc1_new(:,:,1,2) to real space
          ! then performs ccalbecr (rbecp calculation in real space)
          !
          DO ibnd=1,nbnd,2
          !
           !
            !
             CALL fft_orbital_gamma(evc1_new(:,:,1,2),ibnd,nbnd)
             CALL calbec_rs_gamma(ibnd,nbnd,becp%r)
             CALL s_psir_gamma(ibnd,nbnd)
             CALL bfft_orbital_gamma(sevc1_new(:,:,1,2),ibnd,nbnd)
            !
           !
          !
          ENDDO
         !
        !
        !
        ELSE
        ! nkb > 0 & not real space
         !
          !
          CALL calbec(npw,vkb,evc1_new(:,:,1,2),becp)
          CALL s_psi(npwx,npw,nbnd,evc1_new(1,1,1,2),sevc1_new(1,1,1,2))
          !
         !
        !
        ENDIF
        !
       !
       ELSE
       ! nkb = 0, (not an us pp)
        !
         ! This line just copies the array, leave it alone
         CALL s_psi(npwx,npw,nbnd,evc1_new(1,1,1,2),sevc1_new(1,1,1,2))
         !
        !
        !
       !
       ENDIF
    ELSE
       !This is the generalised K point part
       !
       DO ik=1,nks
          !
          IF (  nkb > 0 .and. okvan ) THEN
             CALL init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
             CALL calbec(npw_k(ik),vkb,evc1_new(:,:,ik,2),becp)
          ENDIF
          !
          CALL s_psi(npwx,npw_k(ik),nbnd,evc1_new(1,1,ik,2),sevc1_new(1,1,ik,2))
          !
       ENDDO
         !
    ENDIF
    !print *, "norm of sevc1,1 after spsi", lr_dot(sevc1_new(1,1,1,1),sevc1_new(1,1,1,1))
    !print *, "norm of sevc1,2 after spsi", lr_dot(sevc1_new(1,1,1,1),sevc1_new(1,1,1,2))
    !Resume the LR
    !
    ! Orthogonality requirement as proposed by Y. Saad beta=sqrt(|qdash.pdash|) gamma=sign(qdash.pdash)*beta
    beta=dble(lr_dot(evc1_new(1,1,1,1),sevc1_new(1,1,1,2)))
    !
    !
    IF ( abs(beta)<1.0d-12 ) THEN
       !
       WRITE(stdout,'(5x,"lr_lanczos: Left and right Lanczos vectors are orthogonal, this is a violation of oblique projection")')
       !
    END IF
    IF ( beta<0.0d0 ) THEN
       !
       beta=sqrt(-beta)
       gamma=-beta
       !
    ELSEIF ( beta>0.0d0 ) THEN
       !
       beta=sqrt(beta)
       gamma=beta
       !
    ENDIF
    !
    !
    IF (ionode) THEN
     IF ( charge_response == 1 .and. lr_verbosity > 0) THEN
       !print *, "beta=",beta,"w_T_beta_store", w_T_beta_store(LR_iteration)
        WRITE (stdout,'(5x,"(calc=",e21.15," read=",e21.15,")")') beta, w_T_beta_store(LR_iteration)
        WRITE (stdout,'(5x,"Weight for this step=",2(e11.5,1x))')  w_T(LR_iteration)
     ENDIF
    ENDIF
    beta_store (pol_index,LR_iteration) = beta
    gamma_store(pol_index,LR_iteration) = gamma
    WRITE(stdout,'(5X,"beta (",i8.8,")=",e21.15)') LR_iteration+1,beta
    WRITE(stdout,'(5X,"gamma(",i8.8,")=",e21.15)') LR_iteration+1,gamma
    IF (lr_verbosity > 3) THEN
    IF ( LR_iteration > 6 ) THEN
      WRITE(stdout,'(5X,"(oscillatory variation) : ",f6.2,"%")') abs(100*(beta_store(pol_index,LR_iteration)- &
      (beta_store(pol_index,LR_iteration-6)+beta_store(pol_index,LR_iteration-4)+ &
      beta_store(pol_index,LR_iteration-2))/3.0)/beta_store(pol_index,LR_iteration))
      WRITE(stdout,'(5X,"(linear variation) : ",f6.2,"%")') abs(100*(beta_store(pol_index,LR_iteration)- &
      (beta_store(pol_index,LR_iteration-3)+beta_store(pol_index,LR_iteration-2)+ &
      beta_store(pol_index,LR_iteration-1))/3.0)/beta_store(pol_index,LR_iteration))
    ENDIF
    ENDIF
    !
    !Since beta and gamma are known now, we can calculate the proper q from qdash
    ! V matrix is reset for a new step. Notice that evc1_old and evc1 are scaled so that they are q and p
    ! however evc1_new is qdash and pdash

    !OBM, lets try BLAS

    CALL zcopy(size_evc*2,evc1(1,1,1,1),1,evc1_old(1,1,1,1),1) !evc1_old = evc1
    CALL zcopy(size_evc*2,evc1_new(1,1,1,1),1,evc1(1,1,1,1),1) !evc1 = evc1_new
    !
    CALL zscal(size_evc,cmplx(1.0d0/beta,0.0d0,kind=dp),evc1(1,1,1,1),1)
    CALL zscal(size_evc,cmplx(1.0d0/gamma,0.0d0,kind=dp),evc1(1,1,1,2),1)
    !
    evc1_new(:,:,:,:)=(0.0d0,0.0d0)
    sevc1_new(:,:,:,:)=(0.0d0,0.0d0)
    !
    !
    !
    IF(.not.ltammd) THEN
    !
    IF ( mod(LR_iteration,2)==0 ) THEN
       CALL lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),sevc1_new(1,1,1,1),.false.)
       CALL lr_apply_liouvillian(evc1(1,1,1,2),evc1_new(1,1,1,2),sevc1_new(1,1,1,2),.true.)
    ELSE
       CALL lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),sevc1_new(1,1,1,1),.true.)
       CALL lr_apply_liouvillian(evc1(1,1,1,2),evc1_new(1,1,1,2),sevc1_new(1,1,1,2),.false.)
    ENDIF
    !
    ELSE
       CALL lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),sevc1_new(1,1,1,1),.true.)
       CALL zcopy(size_evc,evc1(1,1,1,1),1,evc1(1,1,1,2),1) !evc1(,1) = evc1(,2)
       CALL zcopy(size_evc,evc1_new(1,1,1,1),1,evc1_new(1,1,1,2),1) !evc1_new(,1) = evc1_new(,2)

    ENDIF
    !
    ! qdash(i+1)=f(q(i))-gamma*q(i-1)
    ! pdash(i+1)=f(p(i))-beta*p(i-1)
    ! where f(p(i)) or f(q(i)) are calculated by lr_apply_liovillian
    !
    !OBM BLAS
    CALL zaxpy(size_evc,-cmplx(gamma,0.0d0,kind=dp),evc1_old(1,1,1,1),1,evc1_new(1,1,1,1),1)
    CALL zaxpy(size_evc,-CMPLX(beta,0.0d0,kind=dp),evc1_old(1,1,1,2),1,evc1_new(1,1,1,2),1)
    !
    CALL stop_clock('one_step')
    !
    RETURN
    !
  END SUBROUTINE one_lanczos_step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE lr_lanczos
