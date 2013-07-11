MODULE lr_lanczos
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OBM :
!
! This subroutine handles two interleaved non-hermitian chains for x and y at the same time,
!
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
! Modified by Xiaochuan Ge to introduce Pseudo Hermitian Algorithm (2013)

  SUBROUTINE one_lanczos_step()
    !
    !   Non-Hermitian Lanczos
    !
    USE io_global,                ONLY : ionode, stdout
    USE kinds,                ONLY : dp
    USE klist,                ONLY : nks,xk
    USE lr_variables,         ONLY : n_ipol, ltammd, pseudo_hermitian, itermax,&
                                     evc1, evc1_new, sevc1_new, evc1_old, sevc1, &
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
    use lr_us,          only : lr_apply_s
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
    real(kind=dp) :: alpha, beta, gamma, angle
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

    ! Apply Liouvillian
    if(.not. ltammd) then
      if(mod(LR_iteration,2)==0) then
        CALL lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),sevc1_new(1,1,1,1),.true.)
        if(.not. pseudo_hermitian) &
          CALL lr_apply_liouvillian(evc1(1,1,1,2),evc1_new(1,1,1,2),sevc1_new(1,1,1,2),.false.)
      ELSE
        CALL lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),sevc1_new(1,1,1,1),.false.)
        if(.not. pseudo_hermitian) &
          CALL lr_apply_liouvillian(evc1(1,1,1,2),evc1_new(1,1,1,2),sevc1_new(1,1,1,2),.true.)
      ENDIF
    ELSE
      CALL lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),sevc1_new(1,1,1,1),.true.)
      if(.not. pseudo_hermitian) then
        CALL zcopy(size_evc,evc1(1,1,1,1),1,evc1(1,1,1,2),1) !evc1(,1) = evc1(,2)
        CALL zcopy(size_evc,evc1_new(1,1,1,1),1,evc1_new(1,1,1,2),1) !evc1_new(,1) = evc1_new(,2)
      endif
    ENDIF

    DO ik=1, nks
      CALL lr_ortho(evc1_new(:,:,ik,1), evc0(:,:,ik), ik, ik, sevc0(:,:,ik),.true.)
      if(.not. pseudo_hermitian) &
        CALL lr_ortho(evc1_new(:,:,ik,2), evc0(:,:,ik), ik, ik, sevc0(:,:,ik),.true.)
    ENDDO

    ! By construction <p|Lq>=0 should be 0, forcing this both conserves resources and increases stability
    alpha=0.0d0
    alpha_store(pol_index,LR_iteration)=alpha
    WRITE(stdout,'(5X,"alpha(",i8.8,")=",e22.15)') LR_iteration,alpha

    ! Apply S* operation
    if(pseudo_hermitian) then
      call lr_apply_s(evc1_new(:,:,1,1),sevc1_new(:,:,1,1))
    else
      call lr_apply_s(evc1(:,:,1,2),sevc1(:,:,1,2))
    endif

    ! Orthogonality requirement: <v|\bar{L}|v>=1
    if(pseudo_hermitian) then
      beta=dble(lr_dot(evc1(1,1,1,1),sevc1_new(1,1,1,1)))
    else
      beta=dble(lr_dot(evc1(1,1,1,1),sevc1(1,1,1,2)))
      !call lr_apply_s(evc1(:,:,1,1),sevc1(:,:,1,1))
      !angle=dble(lr_dot(evc1(1,1,1,1),sevc1(1,1,1,1)))*dble(lr_dot(evc1(1,1,1,2),sevc1(1,1,1,2)))
      !angle=beta/sqrt(angle)
      !angle = acos(angle)
      !print *, "Pol: ", LR_polarization, "Iteration: ",LR_iteration, "Angle:", angle
    endif

    ! Beta is less than 0 is a serious error for pseudo hermitian algorithm
    if( beta < 0.0d0 .and. pseudo_hermitian) call error_beta()
    
    IF ( abs(beta)<1.0d-12 )   WRITE(stdout,'(5x,"lr_lanczos: Left and right Lanczos vectors are orthogonal,&
                    &this is a violation of oblique projection")')
    IF ( beta<0.0d0 ) THEN
       beta=sqrt(-beta)
       gamma=-beta
    ELSEIF ( beta>0.0d0 ) THEN ! Actually this is the only case in pseudo hermitian algorithm
       beta=sqrt(beta)
       gamma=beta
    ENDIF

    beta_store (pol_index,LR_iteration) = beta
    gamma_store(pol_index,LR_iteration) = gamma
    WRITE(stdout,'(5X,"beta (",i8.8,")=",e22.15)') LR_iteration,beta
    WRITE(stdout,'(5X,"gamma(",i8.8,")=",e22.15)') LR_iteration,gamma

    ! renormalize q(i) and Lq(i), also p(i) and Lp(i) in non pseudo hermitian case
    CALL zscal(size_evc,cmplx(1.0d0/beta,0.0d0,kind=dp),evc1(1,1,1,1),1)
    CALL zscal(size_evc,cmplx(1.0d0/beta,0.0d0,kind=dp),evc1_new(1,1,1,1),1)
    if(.not. pseudo_hermitian) then
      CALL zscal(size_evc,cmplx(1.0d0/gamma,0.0d0,kind=dp),evc1(1,1,1,2),1)
      CALL zscal(size_evc,cmplx(1.0d0/gamma,0.0d0,kind=dp),evc1_new(1,1,1,2),1)
    endif

    ! Calculation of zeta coefficients
    if (mod(LR_iteration,2)==0) then
      do ip=1,n_ipol
        zeta = lr_dot(d0psi(:,:,:,ip),evc1(:,:,:,1)) !Why gamma point dot?
        zeta_store (pol_index,ip,LR_iteration) = zeta
        write(stdout,'(5x,"z1= ",1x,i6,2(1x,e22.15))') ip,real(zeta),aimag(zeta)
      enddo
      
       IF (charge_response == 1) THEN
        CALL lr_calc_dens(evc1(:,:,:,1), .true.)
        CALL lr_calc_F(evc1(:,:,:,1))
      ENDIF
    else
      do ip = 1, n_ipol
        zeta = (0.0d0,0.0d0)
        zeta_store (pol_index,ip,LR_iteration) = zeta
        write(stdout,'(5x,"z1= ",1x,i6,2(1x,e22.15))') ip,real(zeta),aimag(zeta)
      enddo
    endif

    ! XC. Ge : q(i+1)=Lq(i)-beta(i)*q(i-1); renormalization will be done in the
    ! begining of the next iteration 
    ! In non pseudo hermitian case, similar operation needs to be done also for p(i)
    CALL zaxpy(size_evc,-cmplx(gamma,0.0d0,kind=dp),evc1_old(1,1,1,1),1,evc1_new(1,1,1,1),1)
    if(.not. pseudo_hermitian) &
      CALL zaxpy(size_evc,-cmplx(beta,0.0d0,kind=dp),evc1_old(1,1,1,2),1,evc1_new(1,1,1,2),1)

    ! To increase stability, apply lr_ortho once more
    DO ik=1, nks
      CALL lr_ortho(evc1_new(:,:,ik,1), evc0(:,:,ik), ik, ik, sevc0(:,:,ik),.true.)
      if(.not. pseudo_hermitian) &
        CALL lr_ortho(evc1_new(:,:,ik,2), evc0(:,:,ik), ik, ik, sevc0(:,:,ik),.true.)
    ENDDO
    
    ! XC. Ge : throw away q(i-1), and make q(i+1) to be the current vector,
    ! be ready for the next iteration. evc1_new will be free again after this step
    CALL zcopy(size_evc,evc1(1,1,1,1),1,evc1_old(1,1,1,1),1) !evc1_old = evc1
    CALL zcopy(size_evc,evc1_new(1,1,1,1),1,evc1(1,1,1,1),1) !evc1 = evc1_new
    if(.not. pseudo_hermitian) then
      CALL zcopy(size_evc,evc1(1,1,1,2),1,evc1_old(1,1,1,2),1) !evc1_old = evc1
      CALL zcopy(size_evc,evc1_new(1,1,1,2),1,evc1(1,1,1,2),1) !evc1 = evc1_new
    endif

    IF (ionode) THEN
      IF ( charge_response == 1 .and. lr_verbosity > 0) THEN
       !print *, "beta=",beta,"w_T_beta_store", w_T_beta_store(LR_iteration)
        WRITE (stdout,'(5x,"(calc=",e22.15," read=",e22.15,")")') beta, w_T_beta_store(LR_iteration)
        WRITE (stdout,'(5x,"Weight for this step=",2(e12.5,1x))')  w_T(LR_iteration)
      ENDIF
    ENDIF

    CALL stop_clock('one_step')
    RETURN
!----------------------------------------------------------------------------------------------
contains

  subroutine error_beta()
    implicit none
    integer :: ib,ik
    WRITE(stdout,'(5x,"Error: Beta is negative:",5x,e22.15)') beta
      do ib = 1, nbnd
        do ik = 1, npwx
          write(301,*) ib,ik,dble(evc1(ik,ib,1,1)),aimag(evc1(ik,ib,1,1))
        enddo
      enddo
      close(301)

      do ib = 1, nbnd
        do ik = 1, npwx
          write(302,*) ib,ik,dble(evc1(ik,ib,1,2)),aimag(evc1(ik,ib,1,2))
        enddo
      enddo
      close(302)

      do ib = 1, nbnd
        do ik = 1, npwx
          write(303,*) ib,ik,dble(sevc1_new(ik,ib,1,2)),aimag(evc1(ik,ib,1,2))
        enddo
      enddo
      close(303)

      stop "Beta is negative"
  end subroutine error_beta

  END SUBROUTINE one_lanczos_step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE lr_lanczos
