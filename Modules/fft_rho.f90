!
! Copyright (C) 2017 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
MODULE fft_rho
  !-----------------------------------------------------------------------------
  !
  ! ... FFT and inverse FFT of rho on the dense grid
  !
  USE kinds,     ONLY : DP
  USE fft_base,       ONLY: dfftp
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvect,          ONLY: ngm,  nl, nlm
  USE control_flags,  ONLY: gamma_only
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rho_r2g, rho_g2r
  !
CONTAINS
  !
  SUBROUTINE rho_r2g ( rhor, rhog )
    !
    REAL(dp),    INTENT(in) :: rhor(:,:)
    COMPLEX(dp), INTENT(OUT):: rhog(:,:)
    !
    INTEGER :: ir, ig, iss, isup, isdw
    INTEGER :: nspin
    COMPLEX(dp):: fp, fm
    COMPLEX(dp), ALLOCATABLE :: psi(:)

    nspin= SIZE (rhor, 2)

    ALLOCATE( psi( dfftp%nnr ) )
    IF( nspin == 1 ) THEN
       iss=1
       DO ir=1,dfftp%nnr
          psi(ir)=CMPLX(rhor(ir,iss),0.0_dp,kind=dp)
       END DO
       CALL fwfft('Dense', psi, dfftp )
       DO ig=1,ngm
          rhog(ig,iss)=psi(nl(ig))
       END DO
    ELSE
       isup=1
       isdw=2
       DO ir=1,dfftp%nnr
          psi(ir)=CMPLX(rhor(ir,isup),rhor(ir,isdw),kind=dp)
       END DO
       CALL fwfft('Dense', psi, dfftp )
       DO ig=1,ngm
          fp=psi(nl(ig))+psi(nlm(ig))
          fm=psi(nl(ig))-psi(nlm(ig))
          rhog(ig,isup)=0.5_dp*CMPLX( DBLE(fp),AIMAG(fm),kind=DP)
          rhog(ig,isdw)=0.5_dp*CMPLX(AIMAG(fp),-DBLE(fm),kind=DP)
       END DO
    ENDIF
    
    DEALLOCATE( psi )

  END SUBROUTINE rho_r2g
  !
  SUBROUTINE rho_g2r ( rhog, rhor )
    !
    COMPLEX(dp), INTENT(in ):: rhog(:,:)
    REAL(dp),    INTENT(out):: rhor(:,:)
    !
    INTEGER :: ir, ig, iss, isup, isdw
    INTEGER :: nspin
    COMPLEX(dp), PARAMETER :: ci=(0.0_dp, 1.0_dp)
    COMPLEX(dp), ALLOCATABLE :: psi(:)

    nspin= SIZE (rhog, 2)

    ALLOCATE( psi( dfftp%nnr ) )
    IF ( gamma_only ) THEN
       IF( nspin == 1 ) THEN
          iss=1
          psi (:) = (0.0_dp, 0.0_dp)
!$omp parallel do
          DO ig=1,ngm
             psi(nlm(ig))=CONJG(rhog(ig,iss))
             psi(nl (ig))=      rhog(ig,iss)
          END DO
!$omp end parallel do
          CALL invfft('Dense',psi, dfftp )
!$omp parallel do
          DO ir=1,dfftp%nnr
             rhor(ir,iss)=DBLE(psi(ir))
          END DO
!$omp end parallel do
       ELSE
          isup=1
          isdw=2
          psi (:) = (0.0_dp, 0.0_dp)
!$omp parallel do
          DO ig=1,ngm
             psi(nlm(ig))=CONJG(rhog(ig,isup))+ci*CONJG(rhog(ig,isdw))
             psi(nl(ig))=rhog(ig,isup)+ci*rhog(ig,isdw)
          END DO
!$omp end parallel do
          CALL invfft('Dense',psi, dfftp )
!$omp parallel do
          DO ir=1,dfftp%nnr
             rhor(ir,isup)= DBLE(psi(ir))
             rhor(ir,isdw)=AIMAG(psi(ir))
          END DO
       !$omp end parallel do
       ENDIF
       !
    ELSE
       !
       DO iss=1, nspin
          psi (:) = (0.0_dp, 0.0_dp)
!$omp parallel do
          DO ig=1,ngm
             psi(nl (ig))=      rhog(ig,iss)
          END DO
!$omp end parallel do
          CALL invfft('Dense',psi, dfftp )
!$omp parallel do
          DO ir=1,dfftp%nnr
             rhor(ir,iss)=DBLE(psi(ir))
          END DO
!$omp end parallel do
       END DO
    END IF
    
    DEALLOCATE( psi )

  END SUBROUTINE rho_g2r

END MODULE fft_rho
