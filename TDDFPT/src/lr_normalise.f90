!-----------------------------------------------------------------------
! ... normalises the two components of a supervector so that they
! ... have an inner product of 1
!-----------------------------------------------------------------------
! Modified by Osman Baris Malcioglu (2009)
SUBROUTINE lr_normalise(evc1,norm)
  !
  !
  USE kinds,                ONLY : dp
  USE gvect,                ONLY : gstart
  USE cell_base,            ONLY : omega
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : nks,xk
  USE lsda_mod,             ONLY : nspin
  USE lr_variables,         ONLY : lanc_norm
  USE realus,               ONLY : igk_k, npw_k
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE wvfct,                ONLY : nbnd, npwx, npw, wg
  USE control_flags,        ONLY : gamma_only
  USE lr_variables,         ONLY : lr_verbosity
  !
  IMPLICIT NONE
  !
  real(kind=dp), INTENT(out) :: norm
  !
  ! local variables
  INTEGER :: ik
  COMPLEX(kind=dp) :: evc1(npwx,nbnd,nks)
  COMPLEX(kind=dp), ALLOCATABLE :: spsi(:,:,:)
  INTEGER :: ibnd,ig
  !
  ALLOCATE(spsi(npwx,nbnd,nks))
  !
  IF (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_normalise>")')
  ENDIF
  IF(gamma_only) THEN
     CALL lr_normalise_gamma()
  ELSE
     CALL lr_normalise_k()
  ENDIF
  !
  DEALLOCATE(spsi)
  !
  RETURN
  !
CONTAINS
  !
!--------------------------------------------------------------------
  SUBROUTINE lr_normalise_gamma()
    !
    USE becmod,         ONLY : bec_type, becp,calbec
    USE realus,         ONLY : real_space, fft_orbital_gamma,&
         & initialisation_level,&
         & bfft_orbital_gamma, calbec_rs_gamma,&
         & add_vuspsir_gamma, v_loc_psir,&
         & s_psir_gamma, real_space_debug 
    !
    !
    !
    IMPLICIT NONE
    !
    REAL(kind=dp) :: prod
    COMPLEX(kind=dp), EXTERNAL :: lr_dot
    INTEGER :: ibnd,ig
    !
    prod=0.0d0
    !
    IF ( nkb > 0 ) THEN
       !
       IF (real_space_debug>6) THEN
          ! real space & nkb > 0
          !
          DO ibnd=1,nbnd,2
             CALL fft_orbital_gamma(evc1(:,:,1),ibnd,nbnd)
             CALL calbec_rs_gamma(ibnd,nbnd,becp%r)
             CALL s_psir_gamma(ibnd,nbnd)
             CALL bfft_orbital_gamma(spsi(:,:,1),ibnd,nbnd)
          ENDDO
          !
          !
       ELSE
          !
          !the non real_space & nkb > 0 case
          !
          CALL calbec(npw_k(1),vkb,evc1(:,:,1),becp)
          !
          CALL s_psi(npwx,npw_k(1),nbnd,evc1(1,1,1),spsi)
          !
          !
       ENDIF
    ELSE
       ! The nkb == 0 part
       ! JUST array copying
       CALL s_psi(npwx,npw_k(1),nbnd,evc1(1,1,1),spsi)
       !
       !
       !
    ENDIF
    !
    prod=dble( lr_dot( evc1(1,1,1),spsi(1,1,1) ) )
    prod=1.0d0/sqrt(abs(prod))
    !
    evc1(:,:,1)=cmplx(prod,0.0d0,dp)*evc1(:,:,1)
    !
    WRITE(stdout,'(5X,"Norm of initial Lanczos vectors=",1x,f21.15)')&
         & 1.0d0/prod
    !
    lanc_norm=1.d0/prod**2/omega
    norm=1.0d0/prod
    !
    RETURN
  END SUBROUTINE lr_normalise_gamma
!--------------------------------------------------------------------
  SUBROUTINE lr_normalise_k()
    !
    USE becmod,              ONLY : becp,calbec
    !
    real(kind=dp) :: prod
    COMPLEX(kind=dp), EXTERNAL :: lr_dot
    !
    prod=0.0d0
    !
    DO ik=1,nks
       !
       IF ( nkb > 0 .and. okvan) THEN
          !
          CALL init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
          !
          CALL calbec(npw_k(ik),vkb,evc1(:,:,ik),becp)
          !
       ENDIF
          !
          CALL s_psi(npwx,npw_k(ik),nbnd,evc1(:,:,ik),spsi(:,:,ik))
          !
    ENDDO
    !
    prod=dble( lr_dot( evc1(1,1,1),spsi(1,1,1) ) )
    prod=1.0d0/sqrt(abs(prod))
    !
    evc1(:,:,:)=cmplx(prod,0.0d0,dp)*evc1(:,:,:)
    !
    WRITE(stdout,'(5X,"Norm of initial Lanczos vectors=",1x,f21.15)')&
         & 1.0d0/prod
    !
    lanc_norm=1.d0/prod**2/omega
    norm=1.0d0/prod
    !
    RETURN
    !
  END SUBROUTINE lr_normalise_k
  !
END SUBROUTINE lr_normalise
!-----------------------------------------------------------------------
