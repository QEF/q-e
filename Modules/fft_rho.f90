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
  PUBLIC :: rho_r2g, rho_g2r
  !
  INTERFACE rho_g2r
    MODULE PROCEDURE rho_g2r_x, rho_g2r_sum_components
  END INTERFACE
  !
CONTAINS
  !
  SUBROUTINE rho_r2g ( desc, rhor, rhog, v )
    USE fft_types,              ONLY: fft_type_descriptor
    USE fft_helper_subroutines, ONLY: fftx_threed2oned
    !
    TYPE(fft_type_descriptor), INTENT(in) :: desc
    REAL(dp),    INTENT(in) :: rhor(:,:)
    COMPLEX(dp), INTENT(OUT):: rhog(:,:)
    REAL(dp),    OPTIONAL, INTENT(in) :: v(:)
    !
    INTEGER :: ir, ig, iss, isup, isdw
    INTEGER :: nspin
    COMPLEX(dp):: fp, fm
    COMPLEX(dp), ALLOCATABLE :: psi(:)

    nspin= SIZE (rhor, 2)

    ALLOCATE( psi( desc%nnr ) )
    IF( nspin == 1 ) THEN
       iss=1
       IF( PRESENT( v ) ) THEN
          DO ir=1,desc%nnr
             psi(ir)=CMPLX(rhor(ir,iss)+v(ir),0.0_dp,kind=dp)
          END DO
       ELSE
          DO ir=1,desc%nnr
             psi(ir)=CMPLX(rhor(ir,iss),0.0_dp,kind=dp)
          END DO
       END IF
       CALL fwfft('Rho', psi, desc )
       CALL fftx_threed2oned( desc, psi, rhog(:,iss) )
    ELSE
       isup=1
       isdw=2
       IF( PRESENT( v ) ) THEN
          DO ir=1,desc%nnr
             psi(ir)=CMPLX(rhor(ir,isup)+v(ir),rhor(ir,isdw)+v(ir),kind=dp)
          END DO
       ELSE
          DO ir=1,desc%nnr
             psi(ir)=CMPLX(rhor(ir,isup),rhor(ir,isdw),kind=dp)
          END DO
       END IF
       CALL fwfft('Rho', psi, desc )
       CALL fftx_threed2oned( desc, psi, rhog(:,isup), rhog(:,isdw) )
    ENDIF
    
    DEALLOCATE( psi )

  END SUBROUTINE rho_r2g
  !
  SUBROUTINE rho_g2r_x ( desc, rhog, rhor )
    USE fft_types,              ONLY: fft_type_descriptor
    USE fft_helper_subroutines, ONLY: fftx_threed2oned, fftx_oned2threed
    !
    TYPE(fft_type_descriptor), INTENT(in) :: desc
    COMPLEX(dp), INTENT(in ):: rhog(:,:)
    REAL(dp),    INTENT(out):: rhor(:,:)
    !
    INTEGER :: ir, ig, iss, isup, isdw
    INTEGER :: nspin
    COMPLEX(dp), ALLOCATABLE :: psi(:)

    nspin= SIZE (rhog, 2)

    ALLOCATE( psi( desc%nnr ) )
    IF ( gamma_only ) THEN
       IF( nspin == 1 ) THEN
          iss=1
          CALL fftx_oned2threed( desc, psi, rhog(:,iss) )
          CALL invfft('Rho',psi, desc )
!$omp parallel do
          DO ir=1,desc%nnr
             rhor(ir,iss)=DBLE(psi(ir))
          END DO
!$omp end parallel do
       ELSE
          isup=1
          isdw=2
          CALL fftx_oned2threed( desc, psi, rhog(:,isup), rhog(:,isdw) )
          CALL invfft('Rho',psi, desc )
!$omp parallel do
          DO ir=1,desc%nnr
             rhor(ir,isup)= DBLE(psi(ir))
             rhor(ir,isdw)=AIMAG(psi(ir))
          END DO
!$omp end parallel do
       ENDIF
       !
    ELSE
       !
       DO iss=1, nspin
          CALL fftx_oned2threed( desc, psi, rhog(:,iss) )
          CALL invfft('Rho',psi, desc )
!$omp parallel do
          DO ir=1,desc%nnr
             rhor(ir,iss)=DBLE(psi(ir))
          END DO
!$omp end parallel do
       END DO
    END IF
    
    DEALLOCATE( psi )

  END SUBROUTINE rho_g2r_x
  !
  SUBROUTINE rho_g2r_sum_components ( desc, rhog, rhor )
    USE fft_types,              ONLY: fft_type_descriptor
    USE fft_helper_subroutines, ONLY: fftx_threed2oned, fftx_oned2threed
    !
    TYPE(fft_type_descriptor), INTENT(in) :: desc
    COMPLEX(dp), INTENT(in ):: rhog(:,:)
    REAL(dp),    INTENT(out):: rhor(:)
    !
    INTEGER :: ir, ig, iss, isup, isdw
    INTEGER :: nspin
    COMPLEX(dp), ALLOCATABLE :: psi(:)

    nspin= SIZE (rhog, 2)

    ALLOCATE( psi( desc%nnr ) )
    IF ( gamma_only ) THEN
       IF( nspin == 1 ) THEN
          iss=1
          CALL fftx_oned2threed( desc, psi, rhog(:,iss) )
          CALL invfft('Rho',psi, desc )
!$omp parallel do
          DO ir=1,desc%nnr
             rhor(ir)=DBLE(psi(ir))
          END DO
!$omp end parallel do
       ELSE
          isup=1
          isdw=2
          CALL fftx_oned2threed( desc, psi, rhog(:,isup), rhog(:,isdw) )
          CALL invfft('Rho',psi, desc )
!$omp parallel do
          DO ir=1,desc%nnr
             rhor(ir)= DBLE(psi(ir))+AIMAG(psi(ir))
          END DO
!$omp end parallel do
       ENDIF
       !
    ELSE
       !
       DO iss=1, nspin
          CALL fftx_oned2threed( desc, psi, rhog(:,iss) )
          CALL invfft('Rho',psi, desc )
          IF( iss == 1 ) THEN
!$omp parallel do
             DO ir=1,desc%nnr
                rhor(ir)=DBLE(psi(ir))
             END DO
!$omp end parallel do
          ELSE
!$omp parallel do
             DO ir=1,desc%nnr
                rhor(ir)=rhor(ir) + DBLE(psi(ir))
             END DO
!$omp end parallel do
          END IF
       END DO
    END IF
    
    DEALLOCATE( psi )

  END SUBROUTINE rho_g2r_sum_components

END MODULE fft_rho
