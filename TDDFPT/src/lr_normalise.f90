!-----------------------------------------------------------------------
! ... normalises the two components of a supervector so that they
! ... have an inner product of 1
!-----------------------------------------------------------------------
! OBM :
! 050608 Modified for calbec interface in v4.0 (evcx(1,1,ik)->evcx(:,:,ik)
!        gamma_only correction
subroutine lr_normalise(evc1,norm)
  !
#include "f_defs.h"
  !
  use gvect,                only : gstart
  use cell_base,            only : omega
  use io_global,            only : stdout
  use kinds,                only : dp
  use klist,                only : nks,xk
  use lsda_mod,             only : nspin
  use lr_variables,         only : lanc_norm
  use realus,               only : igk_k,npw_k
  use uspp,                 only : vkb,nkb,okvan
  use wvfct,                only : nbnd,npwx,npw,wg
  use control_flags,        only : gamma_only
  USE lr_variables,   ONLY : lr_verbosity
  !
  implicit none
  !
  real(kind=dp), intent(out) :: norm
  !
  ! local variables
  integer :: ik
  complex(kind=dp) :: evc1(npwx,nbnd,nks)
  complex(kind=dp), allocatable :: spsi(:,:,:)
  integer :: ibnd,ig
  !
  allocate(spsi(npwx,nbnd,nks)) 
  !
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_normalise>")')
  endif
  if(gamma_only) then
     call lr_normalise_gamma()
  else
     call lr_normalise_k()
  endif
  !
  deallocate(spsi)
  !
  return
  !
contains
  !
  subroutine lr_normalise_gamma()
    !
    use becmod,                   only : bec_type, becp,calbec
    !use lr_variables,             only : real_space
    !use real_beta,                only : ccalbecr_gamma,s_psir,fft_orbital_gamma,bfft_orbital_gamma
      USE realus,               ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                                    bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                    v_loc_psir, s_psir_gamma, real_space_debug

    !
    !
    !
    implicit none
    !
    real(kind=dp) :: prod
    complex(kind=dp), external :: lr_dot
    integer :: ibnd,ig
    !
    prod=0.0d0
    !
    if ( nkb > 0 ) then
     !
     if (real_space_debug>6) then
     ! real space & nkb > 0 
      !
      do ibnd=1,nbnd,2
          call fft_orbital_gamma(evc1(:,:,1),ibnd,nbnd)
          call calbec_rs_gamma(ibnd,nbnd,becp%r)
          call s_psir_gamma(ibnd,nbnd)
          call bfft_orbital_gamma(spsi(:,:,1),ibnd,nbnd)
      enddo
      !
     !
     else
     !
      !the non real_space & nkb > 0 case 
       !
       call calbec(npw_k(1),vkb,evc1(:,:,1),becp)
       !call pw_gemm('Y',nkb,nbnd,npw_k(1),vkb,npwx,evc1(1,1,1),npwx,rbecp,nkb) 
       !
       call s_psi(npwx,npw_k(1),nbnd,evc1(1,1,1),spsi)
      !
     !
     endif
    else
    ! The nkb == 0 part
      ! JUST array copying
       call s_psi(npwx,npw_k(1),nbnd,evc1(1,1,1),spsi)
      !
     !
    !
    endif
    !The below two lines are the replicated part in real space implementation 
    !call calbec(npw_k(1),vkb,evc1(:,:,1),rbecp)
    !call s_psi(npwx,npw_k(1),nbnd,evc1(1,1,1),spsi)
    !
    prod=dble( lr_dot( evc1(1,1,1),spsi(1,1,1) ) )
    prod=1.0d0/sqrt(abs(prod))
    !
    evc1(:,:,1)=cmplx(prod,0.0d0,dp)*evc1(:,:,1)
    !
    write(stdout,'(5X,"Norm of initial Lanczos vectors=",1x,f21.15)') 1.0d0/prod
    lanc_norm=1.d0/prod**2/omega
    norm=1.0d0/prod
    !
    return
  end subroutine lr_normalise_gamma
  !
  subroutine lr_normalise_k()
    !
    use becmod,              only : becp,calbec
    !
    real(kind=dp) :: prod
    complex(kind=dp), external :: lr_dot
    !
    prod=0.0d0
    !
    do ik=1,nks
       !
       if ( nkb > 0 .and. okvan) then
          !
          call init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
          !
          !call ccalbec(nkb,npwx,npw_k(ik),nbnd,becp,vkb,evc1(1,1,ik)) 
          call calbec(npw_k(ik),vkb,evc1(:,:,ik),becp)
          !
       endif
          !
          call s_psi(npwx,npw_k(ik),nbnd,evc1(:,:,ik),spsi(:,:,ik)) 
          !
    end do
    !
    prod=dble( lr_dot( evc1(1,1,1),spsi(1,1,1) ) )
    prod=1.0d0/sqrt(abs(prod))
    !
    evc1(:,:,:)=cmplx(prod,0.0d0,dp)*evc1(:,:,:)
    !
    write(stdout,'(5X,"Norm of initial Lanczos vectors=",1x,f21.15)') 1.0d0/prod
    lanc_norm=1.d0/prod**2/omega
    norm=1.0d0/prod
    !
    return
    !
  end subroutine lr_normalise_k
  !
end subroutine lr_normalise
!-----------------------------------------------------------------------
