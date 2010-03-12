!-----------------------------------------------------------------------
subroutine lr_apply_liouvillian( evc1, evc1_new, sevc1_new, interaction )
  !---------------------------------------------------------------------
  ! ... applies linear response operator to response wavefunctions
  ! OBM: or to be more exact this function is responsible for calculating L.q(i) and (L^T).p(i)
  ! where q is evc1 and return partial qdash(i+1) and pdash(i+1) in evc1_new . Ultrasoft additions are
  ! handled here...
  ! interaction=.true. corresponds to eq.(32) 
  ! interaction=.false. corresponds to eq.(33)
  ! in Ralph Gebauer, Brent Walker  J. Chem. Phys., 127, 164106 (2007) 
  !---------------------------------------------------------------------   
  !
  ! OBM
  ! 050608 : gamma_only correction
  
#include "f_defs.h"
  !
  use ions_base,            only : ityp, nat, ntyp=>nsp
  use cell_base,            only : tpiba2
  use gsmooth,              only : nr1s, nr2s, nr3s,&
       nrx1s, nrx2s, nrx3s, nrxxs, nls, nlsm
  use gvect,                only : nr1, nr2, nr3, nrx1, nrx2, nrx3,&
       nrxx, nl, ngm, gstart, g, gg
  use io_global,            only : stdout
  use kinds,                only : dp
  use klist,                only : nks, xk
  use lr_variables,         only : evc0, revc0, rho_1, lr_verbosity, ltammd, size_evc, no_hxc
  use realus,               only : igk_k,npw_k
  use lsda_mod,             only : nspin
  use uspp,                 only : vkb, nkb, okvan
  use uspp_param,           only : nhm, nh
  use wavefunctions_module, only : psic
  use wvfct,                only : nbnd, npwx, igk, g2kin, et 
  use control_flags,        only : gamma_only
  USE realus,               ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                                   bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                   v_loc_psir, s_psir_gamma, real_space_debug, &
                                   betasave, box_beta, maxbox_beta,newq_r
  USE lr_variables,   ONLY : lr_verbosity, charge_response
  USE io_global,      ONLY : stdout
  USE DFUNCT,         ONLY : newq,add_paw_to_deeq
  USE control_flags,  ONLY : tqr
  

  !
  implicit none
  !
  complex(kind=dp) :: evc1(npwx,nbnd,nks), evc1_new(npwx,nbnd,nks), sevc1_new(npwx,nbnd,nks)
  logical, intent(in) :: interaction
  !
  !   Local variables
  !
  integer :: ir, ibnd, ik, ig, ia, mbia
  integer :: ijkb0, na, nt, ih, jh, ikb, jkb, iqs,jqs
  real(kind=dp), allocatable :: dvrs(:,:), dvrss(:)
  complex(kind=dp), allocatable :: dvrs_temp(:,:) !OBM This waste of memory was already there in lr_dv_of_drho
  real(kind=dp), allocatable :: d_deeq(:,:,:,:)
  complex(kind=dp), allocatable :: spsi1(:,:)
  complex(kind=dp) :: fp, fm
  REAL(DP), allocatable, dimension(:) :: w1, w2
  !
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_apply_liouvillian>")')
  endif
  !

  call start_clock('lr_apply')
  if (interaction) call start_clock('lr_apply_int')
  if (.not.interaction) call start_clock('lr_apply_no')
  !
  allocate( dvrs(nrxx, nspin) )
  allocate( dvrss(nrxxs) )
  dvrs(:,:)=0.0d0
  dvrss(:)=0.0d0
  allocate( d_deeq(nhm, nhm, nat, nspin) )
  d_deeq(:,:,:,:)=0.0d0
  allocate( spsi1(npwx, nbnd) )
  spsi1(:,:)=(0.0d0,0.0d0)
  !
  if( interaction ) then
     !
     call lr_calc_dens( evc1, .false. )
     !
     if (no_hxc) then
     !OBM no_hxc controls the hartree excange correlation addition, if true, they are not added
      dvrs(:,1)=0.0d0
      call interpolate (dvrs(:,1),dvrss,-1)
     else
     dvrs(:,1)=rho_1(:,1)
     !
     !call lr_dv_of_drho(dvrs)
     allocate( dvrs_temp(nrxx, nspin) ) 
       dvrs_temp=CMPLX(dvrs,0.0d0)         !OBM: This memory copy was hidden in lr_dv_of_drho, can it be avoided?
       call dv_of_drho(0,dvrs_temp,.false.)
       dvrs=DBLE(dvrs_temp)
     deallocate(dvrs_temp)
       !
       if ( okvan )  then 
        if ( tqr ) then
         call newq_r(dvrs,d_deeq,.true.)
        else 
         call newq(dvrs,d_deeq,.true.)
        endif
       endif
       call add_paw_to_deeq(d_deeq)
       !
     call interpolate (dvrs(:,1),dvrss,-1)
     endif
     !
  endif
  !
  if( gamma_only ) then
     !
     call lr_apply_liouvillian_gamma()
     !
  else
     !
     call lr_apply_liouvillian_k()
     !
  endif
  !
  if ( interaction .and. (.not.ltammd) ) then
     !
     !   Normal interaction
     !
     write(stdout,'(5X,"lr_apply_liouvillian: applying interaction: normal")')
     !
     !   Here evc1_new contains the interaction
     !
     !OBM, blas
     !sevc1_new=sevc1_new+(1.0d0,0.0d0)*evc1_new
     call zaxpy(size_evc,cmplx(1.0d0,0.0d0,kind=dp),evc1_new(:,:,:),1,sevc1_new(:,:,:),1)
     !
     !
  else if ( interaction .and. ltammd ) then
     !
     !   Tamm-dancoff interaction
     !
     write(stdout,'(5X,"lr_apply_liouvillian: applying interaction: tamm-dancoff")') 
     !
     !   Here evc1_new contains the interaction
     !
     !OBM, blas
     !sevc1_new=sevc1_new+(0.50d0,0.0d0)*evc1_new
     call zaxpy(size_evc,cmplx(0.5d0,0.0d0,kind=dp),evc1_new(:,:,:),1,sevc1_new(:,:,:),1)
     !
     !
  else
     !
     !   Non interacting
     !
     write(stdout,'(5X,"lr_apply_liouvillian: not applying interaction")')
     !
  end if
  !
  if (gstart == 2 .and. gamma_only ) sevc1_new(1,:,:)=cmplx(real(sevc1_new(1,:,:),dp),0.0d0,dp)
  ! (OBM: Why there is this check?)
  if(gstart==2 .and. gamma_only) then
     !
     do ik=1,nks
        !
        do ibnd=1,nbnd
           !
           if(lr_verbosity>6) write(stdout,9000) ibnd,1,sevc1_new(1,ibnd,ik)
           !
           if (abs(aimag(sevc1_new(1,ibnd,ik)))>1.0d-12) then
              !
              call errore(' lr_apply_liouvillian ',&
                   'Imaginary part of G=0 '// &
                   'component does not equal zero',1)
              !
           endif
           !
        enddo
        !
     enddo
     !
  end if
  !
  do ik=1,nks
     !
     call sm1_psi(.false.,ik,npwx,npw_k(ik),nbnd,sevc1_new(1,1,ik),evc1_new(1,1,ik))
     !
  enddo
  !
  deallocate(dvrs)
  deallocate(dvrss)
  deallocate(d_deeq)
  deallocate(spsi1)
  !
9000 FORMAT(/5x,'lr_apply_liouvillian: ibnd=',1X,i2,1X,'sevc1_new(G=0)[',i1,']',2(1x,e12.5)/)
  !
  call stop_clock('lr_apply')
  if (interaction) call stop_clock('lr_apply_int')
  if (.not.interaction) call stop_clock('lr_apply_no')
  !
  return
  !
contains
  !
  subroutine lr_apply_liouvillian_gamma()
    !
    !use becmod,              only : bec_type,becp
    use lr_variables,        only : becp1
    !
    real(kind=dp), allocatable :: becp2(:,:)
    !
    if ( nkb > 0 .and. okvan ) then
       !
       allocate(becp2(nkb,nbnd))
       becp2(:,:)=0.0d0
       !
    end if
    !
    !    Now apply to the ground state wavefunctions
    !    and convert to real space
    !
    if ( interaction ) then
       !
       call start_clock('interaction')
       
       if (nkb > 0 .and. okvan) then
          ! calculation of becp2
          becp2(:,:) = 0.0d0
          !
          ijkb0 = 0
          !
          do nt = 1, ntyp
             !
             do na = 1, nat
                !
                if ( ityp(na) == nt ) then
                   !
                   do ibnd = 1, nbnd
                      !
                      do jh = 1, nh(nt)
                         !
                         jkb = ijkb0 + jh
                         !
                         do ih = 1, nh(nt)
                            !
                            ikb = ijkb0 + ih
                            becp2(ikb, ibnd) = becp2(ikb, ibnd) + &
                                 d_deeq(ih,jh,na,1) * becp1(jkb,ibnd)
                            !
                         enddo
                         !
                      enddo
                      ! 
                   enddo
                   !
                   ijkb0 = ijkb0 + nh(nt)
                   !
                endif
                !
             enddo
             !
          enddo
          !end: calculation of becp2
       endif
       !
       !   evc1_new is used as a container for the interaction
       !
       evc1_new(:,:,:)=(0.0d0,0.0d0)
       !
       do ibnd=1,nbnd,2
          !
          !   Product with the potential vrs = (vltot+vr)
          ! revc0 is on smooth grid. psic is used upto smooth grid
          do ir=1,nrxxs
             !
             psic(ir)=revc0(ir,ibnd,1)*cmplx(dvrss(ir),0.0d0,dp)
             !
          enddo

          !
          !print *,"1"
          if (real_space_debug > 7 .and. okvan .and. nkb > 0) then
          !THE REAL SPACE PART (modified from s_psi)
                  !print *, "lr_apply_liouvillian:Experimental interaction part not using vkb"
                  !fac = sqrt(omega)
                  !
                  ijkb0 = 0
                  iqs = 0
                  jqs = 0
                  !
                  DO nt = 1, ntyp
                  !
                    DO ia = 1, nat
                      !
                      IF ( ityp(ia) == nt ) THEN
                        !
                        mbia = maxbox_beta(ia)
                        !print *, "mbia=",mbia
                        ALLOCATE( w1(nh(nt)),  w2(nh(nt)) )
                        w1 = 0.D0
                        w2 = 0.D0
                        !
                        DO ih = 1, nh(nt)
                        !
                          DO jh = 1, nh(nt)
                          !
                          jkb = ijkb0 + jh
                          w1(ih) = w1(ih) + becp2(jkb, ibnd) 
                          IF ( ibnd+1 .le. nbnd ) w2(ih) = w2(ih) + becp2(jkb, ibnd+1) 
                          !
                          END DO
                        !
                        END DO
                        !
                        !w1 = w1 * fac
                        !w2 = w2 * fac
                        ijkb0 = ijkb0 + nh(nt)
                        !
                        DO ih = 1, nh(nt)
                          !
                          DO ir = 1, mbia
                          !
                           iqs = jqs + ir
                           psic( box_beta(ir,ia) ) = psic(  box_beta(ir,ia) ) + betasave(ia,ih,ir)*CMPLX( w1(ih), w2(ih) )
                          !
                          END DO
                        !
                        jqs = iqs
                        !
                        END DO
                        !
                      DEALLOCATE( w1, w2 )
                        !
                      END IF
                    !
                    END DO
                  !
                  END DO

          endif
          !print *,"2"
          !
          !   Back to reciprocal space This part is equivalent to bfft_orbital_gamma
          !
          call bfft_orbital_gamma (evc1_new(:,:,1), ibnd, nbnd,.false.)
          !print *,"3"
       enddo
       !
       !
       if( nkb > 0 .and. okvan .and. real_space_debug <= 7) then
         !The non real_space part
          call dgemm( 'N', 'N', 2*npw_k(1), nbnd, nkb, 1.d0, vkb, &
               2*npwx, becp2, nkb, 1.d0, evc1_new, 2*npwx )
          !print *, "lr_apply_liouvillian:interaction part using vkb"
          !
       end if
       call stop_clock('interaction')
       !
    end if
    !
    !   Call h_psi on evc1 such that h.evc1 = sevc1_new
    !
    call h_psi(npwx,npw_k(1),nbnd,evc1(1,1,1),sevc1_new(1,1,1))
    !
    ! spsi1 = s*evc1 
    !
    if (real_space_debug > 9 ) then 
        do ibnd=1,nbnd,2
         call fft_orbital_gamma(evc1(:,:,1),ibnd,nbnd)
         call s_psir_gamma(ibnd,nbnd)
         call bfft_orbital_gamma(spsi1,ibnd,nbnd)
        enddo
    else
       call s_psi(npwx,npw_k(1),nbnd,evc1(1,1,1),spsi1)
    endif
    !
    !   Subtract the eigenvalues
    !
    do ibnd=1,nbnd
       !
       call zaxpy(npw_k(1), cmplx(-et(ibnd,1),0.0d0,dp), spsi1(:,ibnd), 1, sevc1_new(:,ibnd,1), 1)
       !
    enddo
    !
    if( nkb > 0 .and. okvan ) deallocate(becp2)
    !
    return
    !
  end subroutine lr_apply_liouvillian_gamma
  !
  subroutine lr_apply_liouvillian_k()
    !
    !use becmod,              only : becp
    use lr_variables,        only : becp1_c
    !
    complex(kind=dp), allocatable :: becp2(:,:)
    !
    if( nkb > 0 .and. okvan ) then
       !
       allocate(becp2(nkb,nbnd))
       becp2(:,:)=(0.0d0,0.0d0)
       !
    endif
    !
    !   Now apply to the ground state wavefunctions
    !   and  convert to real space
    !
    if ( interaction ) then
       !
       call start_clock('interaction')
       !
       !   evc1_new is used as a container for the interaction
       !        
       evc1_new(:,:,:)=(0.0d0,0.0d0)
       !
       do ik=1,nks
          !
          do ibnd=1,nbnd
             !
             !   Product with the potential vrs = (vltot+vr)
             !
             do ir=1,nrxxs
                !
                psic(ir)=revc0(ir,ibnd,ik)*cmplx(dvrss(ir),0.0d0,dp)
                !
             enddo
             !
             !   Back to reciprocal space
             !
             call cft3s(psic,nr1s,nr2s,nr3s,nrx1s,nrx2s,nrx3s,-2)
             !
             do ig=1,npw_k(ik)
                !
                evc1_new(ig,ibnd,ik)=psic(nls(igk_k(ig,ik)))
                !
             enddo
             !
          enddo
          !
       enddo
       !
       call stop_clock('interaction')
       !
       if ( nkb > 0 .and. okvan ) then
          !
          do ik=1,nks
             !
             call init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
             !
             becp2(:,:) = 0.0d0
             !
             ijkb0 = 0
             !
             do nt = 1, ntyp
                !
                do na = 1, nat
                   !
                   if ( ityp(na) == nt ) then
                      !
                      do ibnd = 1, nbnd
                         !
                         do jh = 1, nh(nt)
                            !
                            jkb = ijkb0 + jh
                            !
                            do ih = 1, nh(nt)
                               !
                               ikb = ijkb0 + ih
                               becp2(ikb, ibnd) = becp2(ikb, ibnd) + &
                                    d_deeq(ih,jh,na,1) * becp1_c(jkb,ibnd,ik)
                                !
                            enddo
                            !
                         enddo
                         !
                      enddo
                      !
                      ijkb0 = ijkb0 + nh(nt)
                      !
                   endif
                   !
                enddo
                !
             enddo
             !
             !evc1_new(ik) = evc1_new(ik) + vkb*becp2(ik) 
             call zgemm( 'N', 'N', npw_k(ik), nbnd, nkb, (1.d0,0.d0), vkb, &
                  npwx, becp2, nkb, (1.d0,0.d0), evc1_new(:,:,ik), npwx )
             !
          enddo
          !
       endif
       !
    endif
    !
    !   Call h_psi on evc1
    !   h_psi uses arrays igk and npw, so restore those
    !
    do ik=1,nks
       !
       call init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
       !
       do ig=1,npw_k(ik)
          !
          g2kin(ig)=((xk(1,ik)+g(1,igk_k(ig,ik)))**2 &
               +(xk(2,ik)+g(2,igk_k(ig,ik)))**2 &
               +(xk(3,ik)+g(3,igk_k(ig,ik)))**2)*tpiba2
          !
       enddo
       !
       igk(:)=igk_k(:,ik)
       !
       call h_psi(npwx,npw_k(ik),nbnd,evc1(1,1,ik),sevc1_new(1,1,ik))
       !
       call s_psi(npwx,npw_k(ik),nbnd,evc1(1,1,ik),spsi1)
       !
       !   Subtract the eigenvalues
       !
       do ibnd=1,nbnd
          !
          do ig=1,npw_k(ik)
             !
             sevc1_new(ig,ibnd,ik)=sevc1_new(ig,ibnd,ik) &
                  -cmplx(et(ibnd,ik),0.0d0,dp)*spsi1(ig,ibnd)
             !
          enddo
          !
       enddo
       !
    enddo ! end k loop
    !
    if( nkb > 0 .and. okvan ) deallocate(becp2)
    !
    return
  end subroutine lr_apply_liouvillian_k
  !
end subroutine lr_apply_liouvillian
!-----------------------------------------------------------------------
