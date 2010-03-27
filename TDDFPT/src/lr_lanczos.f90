module lr_lanczos
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OBM :
!
! This subroutine handles two interleaved non-hermitian chains for x and y at the same time, 
!
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
  subroutine one_lanczos_step()
    !
    !   Non-Hermitian Lanczos
    !
    USE io_global,                ONLY : ionode, stdout
    use kinds,                only : dp
    use klist,                only : nks,xk
    use lr_variables,         only : n_ipol, ltammd, itermax,&
                                     evc1, evc1_new, sevc1_new, evc1_old, &
                                     evc0, sevc0, d0psi, &
                                     alpha_store, beta_store, gamma_store, zeta_store,&
                                    charge_response, size_evc, LR_polarization, LR_iteration,&!,real_space
                                     test_case_no
    use uspp,                 only : vkb, nkb, okvan
    use wvfct,                only : nbnd, npwx, npw
    use control_flags,         only : gamma_only,tqr
    use becmod,               only :  bec_type, becp, calbec  
    !use real_beta,           only : ccalbecr_gamma,s_psir,fft_orbital_gamma,bfft_orbital_gamma
    USE realus,               ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                                    bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                    v_loc_psir, s_psir_gamma, igk_k,npw_k, real_space_debug
    USE lr_variables,   ONLY : lr_verbosity, charge_response
    use charg_resp,               only : w_T_beta_store,w_T,lr_calc_F

    !
    implicit none
    !
    complex(kind=dp),external :: lr_dot
    !
    !integer,intent(in) :: iter, pol
    integer :: ik, ip, ibnd,ig, pol_index
    !
    character(len=6), external :: int_to_char
    !
    !   Local variables
    !
    real(kind=dp) :: alpha, beta, gamma
    !
    complex(kind=dp) :: zeta
    !
    !integer :: n
    !
    
    !
    If (lr_verbosity > 5) THEN
      WRITE(stdout,'("<lr_lanczos_one_step>")')
      print *, "Real space = ", real_space
      print *, "Real space debug ", real_space_debug
      print *, "TQR = ", tqr
    endif
    !
    call start_clock('one_step')
    
    !print *, "norm of d0psi", lr_dot(d0psi(1,1,1,1),d0psi(1,1,1,1))
    pol_index=1
    if ( n_ipol /= 1 ) pol_index=LR_polarization
    !
    ! Calculation of zeta coefficients
    !
    
    if (mod(LR_iteration,2)==0) then
       !
       do ip=1,n_ipol
          !
          zeta = lr_dot(d0psi(:,:,:,ip),evc1(:,:,:,1))
          zeta_store (pol_index,ip,LR_iteration) = zeta
          write(stdout,'(5x,"z1= ",1x,i6,2(1x,e21.15))') ip,real(zeta),aimag(zeta)
          !
       end do
       !evc1(:,:,:,1) contains the q of x for even steps, lets calculate the response related observables
       !
       if (charge_response == 2) then
        call lr_calc_dens(evc1(:,:,:,1), .true.)
        call lr_calc_F(evc1(:,:,:,1))
       endif
       !
       !
    else
       !
       do ip=1,n_ipol
          !
          zeta = (0.0d0,0.0d0)
          zeta_store (pol_index,ip,LR_iteration) = zeta
          write(stdout,'(5x,"z1= ",1x,i6,2(1x,e21.15))') ip,real(zeta),aimag(zeta)
          !
       end do
       !
    endif
    !    
    ! Application of the Liouvillian superoperator for the first iteration
    !
    if(LR_iteration==1) then
       LR_iteration=0 !Charge response dump related part is disabled by setting this to zero
       if(.not.ltammd) then
          call lr_apply_liouvillian(evc1(:,:,:,1),evc1_new(:,:,:,1),sevc1_new(:,:,:,1),.false.)
          call lr_apply_liouvillian(evc1(:,:,:,2),evc1_new(:,:,:,2),sevc1_new(:,:,:,2),.true.)
       else
          call lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),sevc1_new(1,1,1,1),.true.)
          evc1(:,:,:,2)=evc1(:,:,:,1)
          evc1_new(:,:,:,2)=evc1_new(:,:,:,1)
          sevc1_new(:,:,:,2)=sevc1_new(:,:,:,1)
       endif
       LR_iteration=1
    endif
    !
    ! The lanczos algorithm starts here
    ! http://www.cs.utk.edu/~dongarra/etemplates/node245.html
    !
    ! Left and right vectors are orthogonalised wrto ground state wf

    !print *, "norm of evc1,1 before lr_ortho", lr_dot(evc1_new(1,1,1,1),evc1_new(1,1,1,1))
    !print *, "norm of evc1,2 before lr_ortho", lr_dot(evc1_new(1,1,1,2),evc1_new(1,1,1,2))
    !print *, "norm of evc0 before lr_ortho", lr_dot(evc0(1,1,1),evc0(1,1,1))
    !print *, "norm of sevc0 before lr_ortho", lr_dot(sevc0(1,1,1),sevc0(1,1,1))
    !print *,"nks=", nks
     
    !OBM: Notice that here "orthogonalization" is not strictly the true word, as the norm of the vectors change
    !This is due to how the uspp scheme is implemented, the beta are evc1_new(left).sevc1_new(right), that is,
    !a mixing of two vectors, thus the resultant vector from belov should be devoid from S, which affects the norm
    !the modification in lr_ortho subroutine handles this reversal (the last flag). 
    do ik=1, nks
     call lr_ortho(evc1_new(:,:,ik,1), evc0(:,:,ik), ik, ik, sevc0(:,:,ik),.true.) 
     call lr_ortho(evc1_new(:,:,ik,2), evc0(:,:,ik), ik, ik, sevc0(:,:,ik),.true.)
    enddo

    !print *, "norm of evc1,1 after lr_ortho", lr_dot(evc1_new(1,1,1,1),evc1_new(1,1,1,1))
    !print *, "norm of evc1,2 after lr_ortho", lr_dot(evc1_new(1,1,1,1),evc1_new(1,1,1,2))


    !
    ! By construction <p|Lq>=0 should be 0, forcing this both conserves resources and increases stability
    alpha=0.0d0
    alpha_store(pol_index,LR_iteration)   = alpha
    write(stdout,'(5X,"alpha(",i8.8,")=",e21.15)') LR_iteration,alpha
    !
    if ( gamma_only ) then
       ! 
       if ( nkb >0 ) then
        !
        if (real_space_debug>5) then
        ! real space & nkb > 0 
         !
          !
          ! The following part converts evc1_new(:,:,1,2) to real space
          ! then performs ccalbecr (rbecp calculation in real space)
          !
          do ibnd=1,nbnd,2
          ! 
           !
            !
             call fft_orbital_gamma(evc1_new(:,:,1,2),ibnd,nbnd)
             call calbec_rs_gamma(ibnd,nbnd,becp%r)
             call s_psir_gamma(ibnd,nbnd)
             call bfft_orbital_gamma(sevc1_new(:,:,1,2),ibnd,nbnd)
            !
           !
          !
          enddo
         !
        !
        !
        else
        ! nkb > 0 & not real space
         ! 
          !
          call calbec(npw,vkb,evc1_new(:,:,1,2),becp)
          call s_psi(npwx,npw,nbnd,evc1_new(1,1,1,2),sevc1_new(1,1,1,2))
          !
         !
        !
        endif
        !
       !
       else
       ! nkb = 0, (not an us pp) 
        !
         ! This line just copies the array, leave it alone
         call s_psi(npwx,npw,nbnd,evc1_new(1,1,1,2),sevc1_new(1,1,1,2))
         !
        !
        !
       !
       endif
    else
       !This is the generalised K point part
       !
       do ik=1,nks
          !
          if (  nkb > 0 .and. okvan ) then
             call init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
             call calbec(npw_k(ik),vkb,evc1_new(:,:,ik,2),becp)
          end if
          !
          call s_psi(npwx,npw_k(ik),nbnd,evc1_new(1,1,ik,2),sevc1_new(1,1,ik,2))
          !
       end do
         !
    end if
    !print *, "norm of sevc1,1 after spsi", lr_dot(sevc1_new(1,1,1,1),sevc1_new(1,1,1,1))
    !print *, "norm of sevc1,2 after spsi", lr_dot(sevc1_new(1,1,1,1),sevc1_new(1,1,1,2))

    !Resume the LR
    !
    ! Orthogonality requirement as proposed by Y. Saad beta=sqrt(|qdash.pdash|) gamma=sign(qdash.pdash)*beta
    beta=dble(lr_dot(evc1_new(1,1,1,1),sevc1_new(1,1,1,2)))
    !
    !
    if ( abs(beta)<1.0d-12 ) then
       !
       write(stdout,'(5x,"lr_lanczos: Left and right Lanczos vectors are orthogonal, this is a violation of oblique projection")')
       !
    else if ( beta<0.0d0 ) then
       !
       beta=sqrt(-beta)
       gamma=-beta
       !
    else if ( beta>0.0d0 ) then
       !
       beta=sqrt(beta)
       gamma=beta
       !
    endif
    !
    ! 
    if (ionode) then
     if ( charge_response == 2 .and. lr_verbosity > 0) then
       !print *, "beta=",beta,"w_T_beta_store", w_T_beta_store(LR_iteration)
        write (stdout,'(5x,"(calc=",e21.15," read=",e21.15,")")') beta, w_T_beta_store(LR_iteration)
        write (stdout,'(5x,"Weight for this step=",e21.15)')  w_T(LR_iteration)
     endif
    endif
    beta_store (pol_index,LR_iteration) = beta
    gamma_store(pol_index,LR_iteration) = gamma
    write(stdout,'(5X,"beta (",i8.8,")=",e21.15)') LR_iteration+1,beta
    write(stdout,'(5X,"gamma(",i8.8,")=",e21.15)') LR_iteration+1,gamma
    if (lr_verbosity > 3) then
    if ( LR_iteration > 6 ) then
      write(stdout,'(5X,"(oscillatory variation) : ",f6.2,"%")') abs(100*(beta_store(pol_index,LR_iteration)- &
      (beta_store(pol_index,LR_iteration-6)+beta_store(pol_index,LR_iteration-4)+ &
      beta_store(pol_index,LR_iteration-2))/3.0)/beta_store(pol_index,LR_iteration))
      write(stdout,'(5X,"(linear variation) : ",f6.2,"%")') abs(100*(beta_store(pol_index,LR_iteration)- &
      (beta_store(pol_index,LR_iteration-3)+beta_store(pol_index,LR_iteration-2)+ &
      beta_store(pol_index,LR_iteration-1))/3.0)/beta_store(pol_index,LR_iteration))
    endif
    endif
    ! 
    !Since beta and gamma are known now, we can calculate the proper q from qdash
    ! V matrix is reset for a new step. Notice that evc1_old and evc1 are scaled so that they are q and p
    ! however evc1_new is qdash and pdash

    !OBM, lets try BLAS
    
    call zcopy(size_evc*2,evc1(:,:,:,:),1,evc1_old(:,:,:,:),1) !evc1_old = evc1
    call zcopy(size_evc*2,evc1_new(:,:,:,:),1,evc1(:,:,:,:),1) !evc1 = evc1_new
    !
    call zscal(size_evc,cmplx(1.0d0/beta,0.0d0,kind=dp),evc1(1,1,1,1),1)
    call zscal(size_evc,cmplx(1.0d0/gamma,0.0d0,kind=dp),evc1(1,1,1,2),1)
    !
    evc1_new(:,:,:,:)=(0.0d0,0.0d0)
    sevc1_new(:,:,:,:)=(0.0d0,0.0d0)
    !
    ! 
    !
    if(.not.ltammd) then
    !
    if ( mod(LR_iteration,2)==0 ) then
       call lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),sevc1_new(1,1,1,1),.false.)
       call lr_apply_liouvillian(evc1(1,1,1,2),evc1_new(1,1,1,2),sevc1_new(1,1,1,2),.true.)
    else
       call lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),sevc1_new(1,1,1,1),.true.) 
       call lr_apply_liouvillian(evc1(1,1,1,2),evc1_new(1,1,1,2),sevc1_new(1,1,1,2),.false.)
    end if
    !
    else
       call lr_apply_liouvillian(evc1(1,1,1,1),evc1_new(1,1,1,1),sevc1_new(1,1,1,1),.true.)
       call zcopy(size_evc,evc1(:,:,:,1),1,evc1(:,:,:,2),1) !evc1(,1) = evc1(,2)
       call zcopy(size_evc,evc1_new(:,:,:,1),1,evc1_new(:,:,:,2),1) !evc1_new(,1) = evc1_new(,2)
       call zcopy(size_evc,evc1_new(:,:,:,1),1,evc1_new(:,:,:,2),1) !evc1_new(,1) = evc1_new(,2)
       
    end if
    !
    ! qdash(i+1)=f(q(i))-gamma*q(i-1)
    ! pdash(i+1)=f(p(i))-beta*p(i-1) 
    ! where f(p(i)) or f(q(i)) are calculated by lr_apply_liovillian 
    !
    !OBM BLAS
    call zaxpy(size_evc,-cmplx(gamma,0.0d0,kind=dp),evc1_old(:,:,:,1),1,evc1_new(:,:,:,1),1)
    call zaxpy(size_evc,-cmplx(beta,0.0d0,kind=dp),evc1_old(:,:,:,2),1,evc1_new(:,:,:,2),1)
    !
    ! Writing files for restart
    !
    !if ( mod(LR_iteration,restart_step)==0 .or. LR_iteration==itermax ) then
    !   !
    !   nwordrestart = 2 * nbnd * npwx * nks
    !   !
    !   call diropn ( iunrestart, 'restart_lanczos.'//trim(int_to_char(LR_polarization)), nwordrestart, exst)
       !
    !   call davcio(evc1(:,:,:,1),nwordrestart,iunrestart,1,1)
    !   call davcio(evc1(:,:,:,2),nwordrestart,iunrestart,2,1)
    !   call davcio(evc1_new(:,:,:,1),nwordrestart,iunrestart,3,1)
    !   call davcio(evc1_new(:,:,:,2),nwordrestart,iunrestart,4,1)
    !   !
    !   close( unit = iunrestart)
    !   if (charge_response == 2 ) then 
    !    call diropn ( iunrestart, 'restart_lanczos-rho_tot.'//trim(int_to_char(LR_polarization)), 2*nrxx, exst)
    !     call davcio(rho_1_tot(:,:),2*nrxx*nspin_mag,iunrestart,1,1)
    !     close( unit = iunrestart)
    !   endif
       !
    !end if
    !
    call stop_clock('one_step')
    !
    return
    !
  end subroutine one_lanczos_step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module lr_lanczos
