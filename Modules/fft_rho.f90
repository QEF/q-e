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
  USE fft_interfaces, ONLY: fwfft, invfft
  USE control_flags,  ONLY: gamma_only
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rho_r2g, rho_g2r, smooth_rho_g2r, smooth_rho_r2g
  !
  INTERFACE rho_g2r
    MODULE PROCEDURE rho_g2r_x, rho_g2r_sum_components
  END INTERFACE
  !
CONTAINS
  !
  SUBROUTINE rho_r2g ( rhor, rhog, v )
    USE gvect,          ONLY: ngm,  nl, nlm
    USE fft_base,       ONLY: dfftp
    !
    REAL(dp),    INTENT(in) :: rhor(:,:)
    COMPLEX(dp), INTENT(OUT):: rhog(:,:)
    REAL(dp),    OPTIONAL, INTENT(in) :: v(:)
    !
    INTEGER :: ir, ig, iss, isup, isdw
    INTEGER :: nspin
    COMPLEX(dp):: fp, fm
    COMPLEX(dp), ALLOCATABLE :: psi(:)

    nspin= SIZE (rhor, 2)

    ALLOCATE( psi( dfftp%nnr ) )
    IF( nspin == 1 ) THEN
       iss=1
       IF( PRESENT( v ) ) THEN
          DO ir=1,dfftp%nnr
             psi(ir)=CMPLX(rhor(ir,iss)+v(ir),0.0_dp,kind=dp)
          END DO
       ELSE
          DO ir=1,dfftp%nnr
             psi(ir)=CMPLX(rhor(ir,iss),0.0_dp,kind=dp)
          END DO
       END IF
       CALL fwfft('Dense', psi, dfftp )
       DO ig=1,ngm
          rhog(ig,iss)=psi(nl(ig))
       END DO
    ELSE
       isup=1
       isdw=2
       IF( PRESENT( v ) ) THEN
          DO ir=1,dfftp%nnr
             psi(ir)=CMPLX(rhor(ir,isup)+v(ir),rhor(ir,isdw)+v(ir),kind=dp)
          END DO
       ELSE
          DO ir=1,dfftp%nnr
             psi(ir)=CMPLX(rhor(ir,isup),rhor(ir,isdw),kind=dp)
          END DO
       END IF
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
  SUBROUTINE smooth_rho_r2g ( rhor, rhog, v )
    USE gvecs,          ONLY: ngms,  nls, nlsm
    USE fft_base,       ONLY: dffts
    !
    REAL(dp),    INTENT(in) :: rhor(:,:)
    COMPLEX(dp), INTENT(OUT):: rhog(:,:)
    REAL(dp),    OPTIONAL, INTENT(in) :: v(:)
    !
    INTEGER :: ir, ig, iss, isup, isdw
    INTEGER :: nspin
    COMPLEX(dp):: fp, fm
    COMPLEX(dp), ALLOCATABLE :: psi(:)

    nspin= SIZE (rhor, 2)

    ALLOCATE( psi( dffts%nnr ) )
    IF( nspin == 1 ) THEN
       iss=1
       IF( PRESENT( v ) ) THEN
          DO ir=1,dffts%nnr
             psi(ir)=CMPLX(rhor(ir,iss)+v(ir),0.0_dp,kind=dp)
          END DO
       ELSE
          DO ir=1,dffts%nnr
             psi(ir)=CMPLX(rhor(ir,iss),0.0_dp,kind=dp)
          END DO
       END IF
       CALL fwfft('Smooth', psi, dffts )
       DO ig=1,ngms
          rhog(ig,iss)=psi(nls(ig))
       END DO
    ELSE
       isup=1
       isdw=2
       IF( PRESENT( v ) ) THEN
          DO ir=1,dffts%nnr
             psi(ir)=CMPLX(rhor(ir,isup)+v(ir),rhor(ir,isdw)+v(ir),kind=dp)
          END DO
       ELSE
          DO ir=1,dffts%nnr
             psi(ir)=CMPLX(rhor(ir,isup),rhor(ir,isdw),kind=dp)
          END DO
       END IF
       CALL fwfft('Smooth', psi, dffts )
       DO ig=1,ngms
          fp=psi(nls(ig))+psi(nlsm(ig))
          fm=psi(nls(ig))-psi(nlsm(ig))
          rhog(ig,isup)=0.5_dp*CMPLX( DBLE(fp),AIMAG(fm),kind=DP)
          rhog(ig,isdw)=0.5_dp*CMPLX(AIMAG(fp),-DBLE(fm),kind=DP)
       END DO
    ENDIF
    
    DEALLOCATE( psi )

  END SUBROUTINE smooth_rho_r2g
  !
  SUBROUTINE rho_g2r_x ( rhog, rhor )
    USE gvect,          ONLY: ngm,  nl, nlm
    USE fft_base,       ONLY: dfftp
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

  END SUBROUTINE rho_g2r_x
  !
  SUBROUTINE rho_g2r_sum_components ( rhog, rhor )
    USE gvect,          ONLY: ngm,  nl, nlm
    USE fft_base,       ONLY: dfftp
    !
    COMPLEX(dp), INTENT(in ):: rhog(:,:)
    REAL(dp),    INTENT(out):: rhor(:)
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
             rhor(ir)=DBLE(psi(ir))
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
             rhor(ir)= DBLE(psi(ir))+AIMAG(psi(ir))
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
          IF( iss == 1 ) THEN
!$omp parallel do
             DO ir=1,dfftp%nnr
                rhor(ir)=DBLE(psi(ir))
             END DO
!$omp end parallel do
          ELSE
!$omp parallel do
             DO ir=1,dfftp%nnr
                rhor(ir)=rhor(ir) + DBLE(psi(ir))
             END DO
!$omp end parallel do
          END IF
       END DO
    END IF
    
    DEALLOCATE( psi )

  END SUBROUTINE rho_g2r_sum_components

  SUBROUTINE smooth_rho_g2r ( rhog, rhor )
    USE gvecs,          ONLY: ngms,  nls, nlsm
    USE fft_base,       ONLY: dffts
    !
    COMPLEX(dp), INTENT(in ):: rhog(:,:)
    REAL(dp),    INTENT(out):: rhor(:,:)
    !
    INTEGER :: ir, ig, iss, isup, isdw
    INTEGER :: nspin
    COMPLEX(dp), PARAMETER :: ci=(0.0_dp, 1.0_dp)
    COMPLEX(dp), ALLOCATABLE :: psi(:)

    nspin= SIZE (rhog, 2)

    ALLOCATE( psi( dffts%nnr ) )
    IF ( gamma_only ) THEN
       IF( nspin == 1 ) THEN
          iss=1
          psi (:) = (0.0_dp, 0.0_dp)
!$omp parallel do
          DO ig=1,ngms
             psi(nlsm(ig))=CONJG(rhog(ig,iss))
             psi(nls (ig))=      rhog(ig,iss)
          END DO
!$omp end parallel do
          CALL invfft('Smooth',psi, dffts )
!$omp parallel do
          DO ir=1,dffts%nnr
             rhor(ir,iss)=DBLE(psi(ir))
          END DO
!$omp end parallel do
       ELSE
          isup=1
          isdw=2
          psi (:) = (0.0_dp, 0.0_dp)
!$omp parallel do
          DO ig=1,ngms
             psi(nlsm(ig))=CONJG(rhog(ig,isup))+ci*CONJG(rhog(ig,isdw))
             psi(nls(ig))=rhog(ig,isup)+ci*rhog(ig,isdw)
          END DO
!$omp end parallel do
          CALL invfft('Smooth',psi, dffts )
!$omp parallel do
          DO ir=1,dffts%nnr
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
          DO ig=1,ngms
             psi(nls (ig))=      rhog(ig,iss)
          END DO
!$omp end parallel do
          CALL invfft('Smooth',psi, dffts )
!$omp parallel do
          DO ir=1,dffts%nnr
             rhor(ir,iss)=DBLE(psi(ir))
          END DO
!$omp end parallel do
       END DO
    END IF
    
    DEALLOCATE( psi )

  END SUBROUTINE smooth_rho_g2r


END MODULE fft_rho
