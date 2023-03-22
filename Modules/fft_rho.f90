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
  !! FFT and inverse FFT of rho on the dense grid.
  !
  USE kinds,          ONLY: DP
  USE fft_interfaces, ONLY: fwfft, invfft
  USE control_flags,  ONLY: gamma_only
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: rho_r2g, rho_g2r
  !
  INTERFACE rho_r2g
    MODULE PROCEDURE rho_r2g_1spin, rho_r2g_Nspin
  END INTERFACE
  !
  INTERFACE rho_g2r
    MODULE PROCEDURE rho_g2r_1spin, rho_g2r_Nspin, rho_g2r_sum_spin
  END INTERFACE
  !
CONTAINS
  !
  !-----------------------------------------------------------------
  SUBROUTINE rho_r2g_1spin( desc, rhor, rhog, v, igs )
    !---------------------------------------------------------------
    !! Bring charge density rho from real to G- space - 1-dimensional
    !! input (so 1 spin component only).
    !
    USE fft_types,              ONLY: fft_type_descriptor
    USE fft_helper_subroutines, ONLY: fftx_threed2oned
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: desc
    REAL(DP), INTENT(IN) :: rhor(:)
    !! rho in real space
    COMPLEX(DP), INTENT(OUT) :: rhog(:,:)
    !! rho in G-space
    REAL(DP), INTENT(IN), OPTIONAL :: v(:)
    INTEGER,  INTENT(IN), OPTIONAL :: igs
    !! index of first G-vector for this processor
    !
    ! ... local variables
    !
    INTEGER :: ir
    COMPLEX(DP) :: fp, fm
    COMPLEX(DP), ALLOCATABLE :: psi(:)
    !
    !$acc data present_or_copyin( rhor, v ) present_or_copyout( rhog )
    !
    ALLOCATE( psi(desc%nnr) )
    !$acc data create( psi )
    !
    IF( PRESENT(v) ) THEN
       !$acc parallel loop
       DO ir = 1, desc%nnr
          psi(ir) = CMPLX(rhor(ir)+v(ir),0.0_DP,kind=DP)
       ENDDO
    ELSE
       !$acc parallel loop
       DO ir = 1, desc%nnr
          psi(ir) = CMPLX(rhor(ir),0.0_DP,kind=DP)
       ENDDO
    ENDIF
    !$acc host_data use_device( psi )
    CALL fwfft( 'Rho', psi, desc )
    !$acc end host_data
    IF ( PRESENT(igs) ) THEN
      CALL fftx_threed2oned( desc, psi, rhog(:,1), iigs=igs )
    ELSE
      CALL fftx_threed2oned( desc, psi, rhog(:,1) )
    ENDIF
    !
    !$acc end data
    DEALLOCATE( psi )
    !
    IF (.NOT.PRESENT(igs) .AND. SIZE(rhog,1)>desc%ngm) THEN
      !$acc kernels
      rhog(desc%ngm+1:,1) = (0.d0,0.d0)
      !$acc end kernels
    ENDIF
    !
    !$acc end data
    !
  END SUBROUTINE rho_r2g_1spin
  !
  !-----------------------------------------------------------------
  SUBROUTINE rho_r2g_Nspin( desc, rhor, rhog, v, igs )
    !---------------------------------------------------------------
    !! Bring charge density rho from real to G- space. N-dimensional
    !! input (unpolarized, polarized, etc.).
    !
    USE fft_types,              ONLY: fft_type_descriptor
    USE fft_helper_subroutines, ONLY: fftx_threed2oned
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: desc
    REAL(DP), INTENT(IN) :: rhor(:,:)
    !! rho in real space
    COMPLEX(DP), INTENT(OUT) :: rhog(:,:)
    !! rho in G-space
    REAL(DP), INTENT(IN), OPTIONAL :: v(:)
    INTEGER,  INTENT(IN), OPTIONAL :: igs
    !! index of first G-vector for this processor
    !
    ! ... local variables
    !
    INTEGER :: ir, ig, iss, isup, isdw
    INTEGER :: nspin
    COMPLEX(DP) :: fp, fm
    COMPLEX(DP), ALLOCATABLE :: psi(:)
    !
    !$acc data present_or_copyin( rhor, v ) present_or_copyout( rhog )
    !
    nspin = SIZE(rhor,2)
    !
    ALLOCATE( psi(desc%nnr) )
    !$acc data create( psi )
    !
    IF( nspin == 1 ) THEN
       !
       iss = 1
       IF( PRESENT(v) ) THEN
          !$acc parallel loop
          DO ir = 1, desc%nnr
             psi(ir) = CMPLX(rhor(ir,iss)+v(ir),0.0_DP,kind=DP)
          ENDDO
       ELSE
          !$acc parallel loop
          DO ir = 1, desc%nnr
             psi(ir) = CMPLX(rhor(ir,iss),0.0_DP,kind=DP)
          ENDDO
       ENDIF
       !$acc host_data use_device( psi )
       CALL fwfft( 'Rho', psi, desc )
       !$acc end host_data
       IF ( PRESENT(igs) ) THEN
         CALL fftx_threed2oned( desc, psi, rhog(:,iss), iigs=igs )
       ELSE
         CALL fftx_threed2oned( desc, psi, rhog(:,iss) )
       ENDIF
       !
    ELSE
       !
       IF ( gamma_only ) THEN
          ! nspin/2 = 1 for LSDA, = 2 for noncolinear
          DO iss = 1, nspin/2
             isup = 1 + (iss-1)*nspin/2 ! 1 for LSDA, 1 and 3 for noncolinear
             isdw = 2 + (iss-1)*nspin/2 ! 2 for LSDA, 2 and 4 for noncolinear
             IF( PRESENT( v ) ) THEN
                !$acc parallel loop
                DO ir = 1, desc%nnr
                    psi(ir) = CMPLX(rhor(ir,isup)+v(ir),rhor(ir,isdw)+v(ir),kind=DP)
                ENDDO
             ELSE
                !$acc parallel loop
                DO ir = 1, desc%nnr
                   psi(ir) = CMPLX(rhor(ir,isup),rhor(ir,isdw),kind=DP)
                ENDDO
             ENDIF
             !$acc host_data use_device( psi )
             CALL fwfft( 'Rho', psi, desc )
             !$acc end host_data
             IF ( PRESENT(igs) ) THEN
               CALL fftx_threed2oned( desc, psi, rhog(:,isup), rhog(:,isdw), iigs=igs )
             ELSE
               CALL fftx_threed2oned( desc, psi, rhog(:,isup), rhog(:,isdw) )
             ENDIF
          ENDDO
       ELSE
          DO iss = 1, nspin
             IF( PRESENT( v ) ) THEN
                !$acc parallel loop
                DO ir = 1, desc%nnr
                    psi(ir) = CMPLX(rhor(ir,iss)+v(ir),0.0_DP,kind=DP)
                ENDDO
             ELSE
                !$acc parallel loop
                DO ir = 1, desc%nnr
                   psi(ir) = CMPLX(rhor(ir,iss),0.0_DP,kind=DP)
                ENDDO
             ENDIF
             !$acc host_data use_device( psi )
             CALL fwfft( 'Rho', psi, desc )
             !$acc end host_data
             IF ( PRESENT(igs) ) THEN
               CALL fftx_threed2oned( desc, psi, rhog(:,iss), iigs=igs )
             ELSE
               CALL fftx_threed2oned( desc, psi, rhog(:,iss) )
             ENDIF
          ENDDO
       ENDIF
    ENDIF
    !$acc end data
    DEALLOCATE( psi )
    !
    IF (.NOT.PRESENT(igs) .AND. SIZE(rhog,1)>desc%ngm) THEN
      !$acc kernels
      rhog(desc%ngm+1:,:) = (0.d0,0.d0)
      !$acc end kernels
    ENDIF
    !
    !$acc end data
    !
  END SUBROUTINE rho_r2g_Nspin
  !
  !------------------------------------------------------------------
  SUBROUTINE rho_g2r_1spin( desc, rhog, rhor )
    !----------------------------------------------------------------
    !! Bring charge density rho from G-space to real space. 1-dimensional
    !! input (1 spin component only).
    !
    USE fft_types,              ONLY: fft_type_descriptor
    USE fft_helper_subroutines, ONLY: fftx_oned2threed
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: desc
    COMPLEX(DP), INTENT(IN)  :: rhog(:)
    REAL(DP),    INTENT(OUT) :: rhor(:)
    !
    INTEGER :: ir
    COMPLEX(DP), ALLOCATABLE :: psi(:)
    !
    ALLOCATE( psi(desc%nnr) )
    !
    !$acc data present_or_copyin(rhog) present_or_copyout(rhor) create(psi)
    !
    CALL fftx_oned2threed( desc, psi, rhog )
    !
    !$acc host_data use_device( psi )
    CALL invfft( 'Rho', psi, desc )
    !$acc end host_data
    !
#if defined(_OPENACC)
!$acc parallel loop
#else
!$omp parallel do
#endif
    DO ir = 1, desc%nnr
       rhor(ir) = DBLE(psi(ir))
    ENDDO
#if !defined(_OPENACC)
!$omp end parallel do
#endif
    !
    !$acc end data
    DEALLOCATE( psi )
    !
  END SUBROUTINE rho_g2r_1spin
  !
  !----------------------------------------------------------------------
  SUBROUTINE rho_g2r_Nspin( desc, rhog, rhor )
    !---------------------------------------------------------------------
    !! Bring charge density rho from G- to real space. N-dimensional
    !! input (unpolarized, polarized, etc.).
    !
    USE fft_types,              ONLY: fft_type_descriptor
    USE fft_helper_subroutines, ONLY: fftx_threed2oned, fftx_oned2threed
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: desc
    COMPLEX(DP), INTENT(IN)  :: rhog(:,:)
    REAL(DP),    INTENT(OUT) :: rhor(:,:)
    !
    INTEGER :: ir, ig, iss, isup, isdw
    INTEGER :: nspin
    COMPLEX(DP), ALLOCATABLE :: psi(:)
    !
    nspin = SIZE(rhog,2)
    !
    ALLOCATE( psi( desc%nnr ) )
    !
    !$acc data present_or_copyin(rhog) present_or_copyout(rhor) create(psi)
    !
    IF ( gamma_only ) THEN
       IF( nspin == 1 ) THEN
          iss=1
          !
          CALL fftx_oned2threed( desc, psi, rhog(:,iss) )
          !
          !$acc host_data use_device( psi )
          CALL invfft( 'Rho', psi, desc )
          !$acc end host_data
          !
#if defined(_OPENACC)
!$acc parallel loop
#else
!$omp parallel do
#endif
          DO ir = 1, desc%nnr
             rhor(ir,iss) = DBLE(psi(ir))
          ENDDO
#if !defined(_OPENACC)
!$omp end parallel do
#endif
       ELSE
          ! nspin/2 = 1 for LSDA, = 2 for noncolinear
          DO iss = 1, nspin/2
             isup = 1+(iss-1)*nspin/2 ! 1 for LSDA, 1 and 3 for noncolinear
             isdw = 2+(iss-1)*nspin/2 ! 2 for LSDA, 2 and 4 for noncolinear
             !
             CALL fftx_oned2threed( desc, psi, rhog(:,isup), rhog(:,isdw) )
             !
             !$acc host_data use_device( psi )
             CALL invfft( 'Rho', psi, desc )
             !$acc end host_data
             !
#if defined(_OPENACC)
!$acc parallel loop
#else
!$omp parallel do
#endif
             DO ir = 1, desc%nnr
                rhor(ir,isup) = DBLE(psi(ir))
                rhor(ir,isdw) = AIMAG(psi(ir))
             ENDDO
#if !defined(_OPENACC)
!$omp end parallel do
#endif
          ENDDO
       ENDIF
       !
    ELSE
       !
       DO iss = 1, nspin
          !
          CALL fftx_oned2threed( desc, psi, rhog(:,iss) )
          !
          !$acc host_data use_device( psi )
          CALL invfft( 'Rho', psi, desc )
          !$acc end host_data
          !
#if defined(_OPENACC)
!$acc parallel loop
#else
!$omp parallel do
#endif
          DO ir = 1, desc%nnr
             rhor(ir,iss) = DBLE(psi(ir))
          ENDDO
#if !defined(_OPENACC)
!$omp end parallel do
#endif
       ENDDO
    ENDIF
    !
    !$acc end data
    !
    DEALLOCATE( psi )
    !
  END SUBROUTINE rho_g2r_Nspin
  !
  !-------------------------------------------------------------------------
  SUBROUTINE rho_g2r_sum_spin( desc, rhog, rhor )
    !-----------------------------------------------------------------------
    !! Bring charge density rho from G- to real space and sum over spin
    !! components.
    !
    USE fft_types,              ONLY: fft_type_descriptor
    USE fft_helper_subroutines, ONLY: fftx_threed2oned, fftx_oned2threed
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor), INTENT(IN) :: desc
    COMPLEX(DP), INTENT(IN)  :: rhog(:,:)
    REAL(DP),    INTENT(OUT) :: rhor(:)
    !
    INTEGER :: ir, ig, iss, isup, isdw
    INTEGER :: nspin
    COMPLEX(DP), ALLOCATABLE :: psi(:)
    !
    nspin = SIZE(rhog,2)
    !
    ALLOCATE( psi(desc%nnr) )
    !
    !$acc data present_or_copyin(rhog) present_or_copyout(rhor) create(psi)
    !
    IF ( gamma_only ) THEN
       IF( nspin == 1 ) THEN
          iss = 1
          !
          CALL fftx_oned2threed( desc, psi, rhog(:,iss) )
          !
          !$acc host_data use_device( psi )
          CALL invfft( 'Rho', psi, desc )
          !$acc end host_data
          !
#if defined(_OPENACC)
!$acc parallel loop
#else
!$omp parallel do
#endif
          DO ir = 1, desc%nnr
             rhor(ir) = DBLE(psi(ir))
          ENDDO
#if !defined(_OPENACC)
!$omp end parallel do
#endif
          !
       ELSEIF ( nspin == 2) THEN
          !
          isup = 1
          isdw = 2
          !
          CALL fftx_oned2threed( desc, psi, rhog(:,isup), rhog(:,isdw) )
          !
          !$acc host_data use_device( psi )
          CALL invfft( 'Rho', psi, desc )
          !$acc end host_data
          !
#if defined(_OPENACC)
!$acc parallel loop
#else
!$omp parallel do
#endif
          DO ir = 1, desc%nnr
             rhor(ir) = DBLE(psi(ir)) + AIMAG(psi(ir))
          ENDDO
#if !defined(_OPENACC)
!$omp end parallel do
#endif
          !
       ELSE
          CALL errore( 'rho_g2r_sum_components', 'noncolinear case?', nspin )
       ENDIF
       !
    ELSE
       !
       DO iss = 1, nspin
          !
          CALL fftx_oned2threed( desc, psi, rhog(:,iss) )
          !
          !$acc host_data use_device( psi )
          CALL invfft( 'Rho', psi, desc )
          !$acc end host_data
          !
          IF( iss == 1 ) THEN
#if defined(_OPENACC)
!$acc parallel loop
#else
!$omp parallel do
#endif
             DO ir = 1, desc%nnr
                rhor(ir) = DBLE(psi(ir))
             ENDDO
#if !defined(_OPENACC)
!$omp end parallel do
#endif
          ELSE
#if defined(_OPENACC)
!$acc parallel loop
#else
!$omp parallel do
#endif
             DO ir = 1, desc%nnr
                rhor(ir) = rhor(ir) + DBLE(psi(ir))
             ENDDO
#if !defined(_OPENACC)
!$omp end parallel do
#endif
          ENDIF
       ENDDO
    ENDIF
    !
    !$acc end data
    !
    DEALLOCATE( psi )
    !
  END SUBROUTINE rho_g2r_sum_spin
  !
  !
END MODULE fft_rho
