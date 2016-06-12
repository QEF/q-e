!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE lr_normalise (evc1, norm)
  !--------------------------------------------------------------------
  !
  ! This subroutine normalises the two components of a super-vector 
  ! so that they have an inner product equal to 1.
  !
  ! I. Timrov's note: This subroutine is not used any longer
  ! after a modification of the Lanczos algorithm by X. Ge.
  !
  USE kinds,                ONLY : dp
  USE gvect,                ONLY : gstart
  USE cell_base,            ONLY : omega
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE lsda_mod,             ONLY : nspin
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE wvfct,                ONLY : nbnd, npwx, wg
  USE control_flags,        ONLY : gamma_only
  USE lr_variables,         ONLY : lr_verbosity, eels
  USE noncollin_module,     ONLY : npol
  USE qpoint,               ONLY : nksq
  !
  IMPLICIT NONE
  !
  COMPLEX(kind=dp), INTENT(inout) :: evc1(npwx*npol,nbnd,nksq)
  REAL(kind=dp), INTENT(out) :: norm
  !
  ! local variables
  !
  INTEGER :: ik, ibnd
  COMPLEX(kind=dp), ALLOCATABLE :: spsi(:,:,:)
  !
  ALLOCATE(spsi(npwx*npol,nbnd,nksq))
  !spsi(:,:,:) = (0.0d0, 0.0d0)
  !
  IF (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_normalise>")')
  ENDIF
  !
  IF (eels) THEN
     !
     call lr_normalise_k_eels()
     !
  ELSE
     !
     IF (gamma_only) THEN
        CALL lr_normalise_gamma()
     ELSE
        CALL lr_normalise_k()
     ENDIF
     !
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
    ! Optical case: gamma_only case
    !
    USE becmod,         ONLY : bec_type, becp,calbec
    USE realus,         ONLY : real_space, invfft_orbital_gamma, initialisation_level,   &
                             & fwfft_orbital_gamma, calbec_rs_gamma,add_vuspsir_gamma, &
                             & v_loc_psir, s_psir_gamma, real_space_debug 
    
    IMPLICIT NONE
    REAL(kind=dp) :: prod
    COMPLEX(kind=dp), EXTERNAL :: lr_dot
    !
    ! Calculation of spsi : spsi = S * evc1
    !
    IF ( nkb > 0 ) THEN
       !
       IF (real_space_debug>6) THEN
          !
          ! real space & nkb > 0
          !
          DO ibnd=1,nbnd,2
             CALL invfft_orbital_gamma(evc1(:,:,1),ibnd,nbnd)
             CALL calbec_rs_gamma(ibnd,nbnd,becp%r)
             CALL s_psir_gamma(ibnd,nbnd)
             CALL fwfft_orbital_gamma(spsi(:,:,1),ibnd,nbnd)
          ENDDO
          !
       ELSE
          !
          ! Non real_space & nkb > 0 case
          !
          CALL calbec(ngk(1),vkb,evc1(:,:,1),becp)
          CALL s_psi(npwx,ngk(1),nbnd,evc1(1,1,1),spsi)
          !
       ENDIF
       !
    ELSE
       !
       ! nkb = 0 (just array copying)
       !
       CALL s_psi(npwx,ngk(1),nbnd,evc1(1,1,1),spsi)
       !
    ENDIF
    !
    prod = 0.0d0
    prod = dble( lr_dot( evc1(1,1,1),spsi(1,1,1) ) )
    norm = sqrt(abs(prod))
    prod = 1.0d0/norm
    !
    evc1(:,:,1)=cmplx(prod,0.0d0,dp)*evc1(:,:,1)
    !
    WRITE(stdout,'(5X,"Norm of initial Lanczos vectors=",1x,f21.15)') norm
    !
    RETURN
    !
  END SUBROUTINE lr_normalise_gamma
!--------------------------------------------------------------------
  SUBROUTINE lr_normalise_k()
    !
    ! Optical case: generalized k-point case
    !
    USE becmod,              ONLY : becp,calbec
    !
    IMPLICIT NONE
    REAL(kind=dp) :: prod
    COMPLEX(kind=dp), EXTERNAL :: lr_dot
    !
    DO ik = 1, nks  
       !
       IF ( nkb > 0 .and. okvan) THEN
          !
          CALL init_us_2(ngk(ik),igk_k(1,ik),xk(1,ik),vkb)
          CALL calbec(ngk(ik),vkb,evc1(:,:,ik),becp)
          !
       ENDIF
       !
       CALL s_psi(npwx,ngk(ik),nbnd,evc1(:,:,ik),spsi(:,:,ik))
       !
    ENDDO
    !
    prod = 0.0d0
    prod = dble( lr_dot( evc1(1,1,1),spsi(1,1,1) ) )
    norm = sqrt(abs(prod))
    prod = 1.0d0/norm
    !
    evc1(:,:,:) = cmplx(prod,0.0d0,dp) * evc1(:,:,:)
    !
    WRITE(stdout,'(5X,"Norm of initial Lanczos vectors=",1x,f21.15)') norm
    !
    RETURN
    !
  END SUBROUTINE lr_normalise_k
!-------------------------------------------------------------------------
  SUBROUTINE lr_normalise_k_eels()
    !
    ! EELS: generalized k-point case
    !
    use becmod,              only : becp, calbec
    use qpoint,              only : ikks, ikqs
    use control_lr,          only : nbnd_occ
    !
    IMPLICIT NONE
    !
    INTEGER :: ik,  &
               ikk, & ! index of the point k
               ikq, & ! index of the point k+q
               npwq   ! number of the plane-waves at point k+q
    REAL(kind=dp) :: prod
    COMPLEX(kind=dp), EXTERNAL :: lr_dot
    !
    ! Calculation of spsi : spsi = S * evc1 
    !
    DO ik = 1, nksq
       !
       ikk  = ikks(ik)
       ikq  = ikqs(ik)
       npwq = ngk(ikq)
       !
       IF ( okvan .and. nkb > 0 ) THEN
          !
          ! Calculate beta-functions vkb at point k+q
          !
          CALL init_us_2(npwq, igk_k(1,ikq), xk(1,ikq), vkb)
          !
          ! Calculate the product of beta-functions vkb with
          ! the response orbitals evc1 : becp%k = <vkb|evc1>
          !
          CALL calbec(npwq, vkb, evc1(:,:,ik), becp, nbnd_occ(ikk))
          !
       ENDIF
       !
       CALL s_psi(npwx, npwq, nbnd_occ(ikk), evc1(:,:,ik), spsi(:,:,ik))
       !
    ENDDO
    !
    prod = 0.0d0
    prod = dble( lr_dot( evc1(1,1,1), spsi(1,1,1) ) )
    norm = sqrt(abs(prod))
    prod = 1.0d0/norm
    !
    evc1(:,:,:) = cmplx(prod,0.0d0,dp) * evc1(:,:,:)
    !
    WRITE(stdout,'(5X,"Norm of initial Lanczos vectors=",1x,f21.15)') norm
    !
    RETURN
    !
END SUBROUTINE lr_normalise_k_eels

END SUBROUTINE lr_normalise
!-----------------------------------------------------------------------
