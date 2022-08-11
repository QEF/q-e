!
! Copyright (C) Quantum ESPRESSO Foundation
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------------------
MODULE fft_helper_subroutines
  !--------------------------------------------------------------------------------------
  !! Contains subroutines aimed at rearrange wave and density arrays to/from Fourier 
  !! space grid.
  !
  USE fft_param
  USE fft_types,  ONLY: fft_type_descriptor
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTERFACE tg_reduce_rho
    MODULE PROCEDURE tg_reduce_rho_1,tg_reduce_rho_2,tg_reduce_rho_3,tg_reduce_rho_4, &
                     tg_reduce_rho_5
  END INTERFACE
  !
  PRIVATE
  !
  PUBLIC :: fftx_threed2oned, fftx_oned2threed
  PUBLIC :: tg_reduce_rho
  PUBLIC :: tg_get_nnr, tg_get_recip_inc, fftx_ntgrp, fftx_tgpe, &
            tg_get_group_nr3
  ! ... Used only in CP
  PUBLIC :: fftx_add_threed2oned_gamma, fftx_psi2c_gamma, c2psi_gamma, &
            fftx_add_field, c2psi_gamma_tg, c2psi_k, c2psi_k_tg
  PUBLIC :: fft_dist_info
  ! ... Used only in CP+EXX
  PUBLIC :: fftx_tgcomm
  !
#if defined(__CUDA)
  PUBLIC :: fftx_psi2c_gamma_gpu, c2psi_gamma_gpu
  !
  ! ... nlm and nl array: hold conversion indices from 3D to
  !     1-D vectors. Columns along the z-direction are stored
  !     contigiously.
  INTEGER, POINTER, DEVICE :: nl_d(:), nlm_d(:)
#else
  INTEGER, ALLOCATABLE :: nl_d(:), nlm_d(:)
#endif
  !
CONTAINS
  !
  !------------------------------------------------------------------------------------------
  SUBROUTINE tg_reduce_rho_1( rhos, tg_rho_nc, tg_rho, ispin, noncolin, domag, desc )
     !---------------------------------------------------------------------------------------
     !! Get charge density by summation over task groups. Non-collinear case enabled.
     !! Real-space.
     !
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(IN) :: ispin
     LOGICAL, INTENT(IN) :: noncolin, domag
     REAL(DP), INTENT(INOUT)  :: tg_rho(:)
     REAL(DP), INTENT(INOUT)  :: tg_rho_nc(:,:)
     REAL(DP), INTENT(OUT) :: rhos(:,:)
     !
     INTEGER :: ierr, ioff, ir3, ir, ipol, ioff_tg, nxyp, npol_
     !     write (*,*) ' enter tg_reduce_rho_1'
     !
#if defined(__MPI)
     IF( noncolin) THEN
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tg_rho_nc, SIZE(tg_rho_nc), &
                            MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
     ELSE
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tg_rho, SIZE(tg_rho), &
                            MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
     ENDIF
#endif
     !
     ! copy the charge back to the proper processor location
     !
     nxyp = desc%nr1x * desc%my_nr2p
     IF (noncolin) THEN
         npol_ = 1 ; IF (domag) npol_ = 4
     ENDIF
     DO ir3=1,desc%my_nr3p
        ioff    = desc%nr1x * desc%my_nr2p * (ir3-1)
        ioff_tg = desc%nr1x * desc%nr2x    * (ir3-1) + desc%nr1x * desc%my_i0r2p
        IF (noncolin) THEN
!$omp parallel do
           DO ipol=1, npol_
              DO ir = 1, nxyp
                 rhos(ir+ioff,ipol) = rhos(ir+ioff,ipol) + tg_rho_nc(ir+ioff_tg,ipol)
              END DO
           END DO
!$omp end parallel do
        ELSE
!$omp parallel do
           DO ir = 1, nxyp
              rhos(ir+ioff,ispin) = rhos(ir+ioff,ispin) + tg_rho(ir+ioff_tg)
           END DO
!$omp end parallel do
        END IF
     END DO
     !
  END SUBROUTINE
  !
  !----------------------------------------------------------
  SUBROUTINE tg_reduce_rho_2( rhos, tmp_rhos, ispin, desc )
     !--------------------------------------------------------
     !! Get charge density by summation over task groups - 1 selected
     !! spin component only. Real-space.
     !
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(IN) :: ispin
     REAL(DP), INTENT(INOUT)  :: tmp_rhos(:)
     REAL(DP), INTENT(OUT) :: rhos(:,:)
     !
     INTEGER :: ierr, ioff, ir3, nxyp, ioff_tg
     !     write (*,*) ' enter tg_reduce_rho_2'
     !
     IF ( desc%nproc2 > 1 ) THEN
#if defined(__MPI)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_rhos, SIZE(tmp_rhos), MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
#endif
     ENDIF
     !
     !BRING CHARGE DENSITY BACK TO ITS ORIGINAL POSITION
     !
     nxyp = desc%nr1x * desc%my_nr2p
     DO ir3 = 1, desc%my_nr3p
        ioff    = desc%nr1x * desc%my_nr2p * (ir3-1)
        ioff_tg = desc%nr1x * desc%nr2x    * (ir3-1) + desc%nr1x * desc%my_i0r2p
        rhos(ioff+1:ioff+nxyp,ispin) = rhos(ioff+1:ioff+nxyp,ispin) + tmp_rhos(ioff_tg+1:ioff_tg+nxyp)
     END DO
     !
  END SUBROUTINE
  !
  !--------------------------------------------------------------
  SUBROUTINE tg_reduce_rho_3( rhos, tmp_rhos, desc )
     !-----------------------------------------------------------
     !! Get charge density by summation over task groups - multiple
     !! spin components. Real space.
     !
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     REAL(DP), INTENT(INOUT)  :: tmp_rhos(:,:)
     REAL(DP), INTENT(OUT) :: rhos(:,:)
     !
     INTEGER :: ierr, from, ir3, ioff, nxyp, ioff_tg
     !     write (*,*) ' enter tg_reduce_rho_3'
     !
     IF ( desc%nproc2 > 1 ) THEN
#if defined(__MPI)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_rhos, SIZE(tmp_rhos), MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
#endif
     ENDIF
     !
     !BRING CHARGE DENSITY BACK TO ITS ORIGINAL POSITION
     !
     !If the current processor is not the "first" processor in its
     !orbital group then does a local copy (reshuffling) of its data
     !
     nxyp = desc%nr1x * desc%my_nr2p
     DO ir3 = 1, desc%my_nr3p
        ioff    = desc%nr1x * desc%my_nr2p * (ir3-1)
        ioff_tg = desc%nr1x * desc%nr2x    * (ir3-1) + desc%nr1x * desc%my_i0r2p
        rhos(ioff+1:ioff+nxyp,:) = rhos(ioff+1:ioff+nxyp,:) + tmp_rhos(ioff_tg+1:ioff_tg+nxyp,:)
     END DO
     !
  END SUBROUTINE
  !
  !-----------------------------------------------------------
  SUBROUTINE tg_reduce_rho_4( rhos, tmp_rhos, desc )
     !--------------------------------------------------------
     !! Get charge density by summation over task groups. 1-d
     !! input arrays only. G-space
     !
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     COMPLEX(DP), INTENT(INOUT)  :: tmp_rhos(:)
     COMPLEX(DP), INTENT(OUT) :: rhos(:)
     !
     INTEGER :: ierr, from, ir3, ioff, nxyp, ioff_tg
     !     write (*,*) ' enter tg_reduce_rho_4'
     !
     IF ( desc%nproc2 > 1 ) THEN
#if defined(__MPI)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_rhos, 2*SIZE(tmp_rhos), MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
#endif
     ENDIF
     !
     !BRING CHARGE DENSITY BACK TO ITS ORIGINAL POSITION
     !
     !If the current processor is not the "first" processor in its
     !orbital group then does a local copy (reshuffling) of its data
     !
     nxyp = desc%nr1x * desc%my_nr2p
     DO ir3 = 1, desc%my_nr3p
        ioff    = desc%nr1x * desc%my_nr2p * (ir3-1)
        ioff_tg = desc%nr1x * desc%nr2x    * (ir3-1) + desc%nr1x * desc%my_i0r2p
        rhos(ioff+1:ioff+nxyp) = rhos(ioff+1:ioff+nxyp) + tmp_rhos(ioff_tg+1:ioff_tg+nxyp)
     END DO
     !
  END SUBROUTINE
  !
  !-------------------------------------------------------
  SUBROUTINE tg_reduce_rho_5( rhos, tmp_rhos, desc )
     !-----------------------------------------------------
     !! Get charge density by summation over task groups - multiple
     !! spin components. G-space.
     !
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     COMPLEX(DP), INTENT(INOUT)  :: tmp_rhos(:,:)
     COMPLEX(DP), INTENT(OUT) :: rhos(:,:)
     !
     INTEGER :: ierr, from, ir3, ioff, nxyp, ioff_tg
     !     write (*,*) ' enter tg_reduce_rho_5'
     !
     IF ( desc%nproc2 > 1 ) THEN
#if defined(__MPI)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_rhos, 2*SIZE(tmp_rhos), MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
#endif
     ENDIF
     !
     !BRING CHARGE DENSITY BACK TO ITS ORIGINAL POSITION
     !
     !If the current processor is not the "first" processor in its
     !orbital group then does a local copy (reshuffling) of its data
     !
     nxyp = desc%nr1x * desc%my_nr2p
     DO ir3 = 1, desc%my_nr3p
        ioff    = desc%nr1x * desc%my_nr2p * (ir3-1)
        ioff_tg = desc%nr1x * desc%nr2x    * (ir3-1) + desc%nr1x * desc%my_i0r2p
        rhos(ioff+1:ioff+nxyp,:) = rhos(ioff+1:ioff+nxyp,:) + tmp_rhos(ioff_tg+1:ioff_tg+nxyp,:)
     END DO
     !
  END SUBROUTINE
  !
  !-------------------------------------------------
  SUBROUTINE tg_get_nnr( desc, right_nnr )
     !----------------------------------------------
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: right_nnr
     right_nnr = desc%nnr
  END SUBROUTINE
  !
  !-------------------------------------------------
  SUBROUTINE tg_get_local_nr3( desc, val )
     !----------------------------------------------
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: val
     val = desc%my_nr3p
  END SUBROUTINE
  !
  !------------------------------------------------
  SUBROUTINE tg_get_group_nr3( desc, val )
     !---------------------------------------------
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: val
     val = desc%my_nr3p
  END SUBROUTINE
  !
  !------------------------------------------------
  SUBROUTINE tg_get_recip_inc( desc, val )
     !---------------------------------------------
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: val
     val = desc%nnr
  END SUBROUTINE
  !
  !------------------------------------------------
  PURE FUNCTION fftx_ntgrp( desc )
     !---------------------------------------------
     INTEGER :: fftx_ntgrp
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     fftx_ntgrp = desc%nproc2
  END FUNCTION
  !
  !------------------------------------------------
  PURE FUNCTION fftx_tgpe( desc )
     !---------------------------------------------
     INTEGER :: fftx_tgpe
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     fftx_tgpe = desc%mype2
  END FUNCTION
  !
  !------------------------------------------------
  PURE FUNCTION fftx_tgcomm( desc )
     !---------------------------------------------
     INTEGER :: fftx_tgcomm
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     fftx_tgcomm = desc%comm2
  END FUNCTION
  !
  !------------------------------------------------
  SUBROUTINE fftx_add_field( r, f, desc )
     !---------------------------------------------
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     REAL(DP), INTENT(INOUT)  :: r(:,:)
     REAL(DP), INTENT(IN) :: f(:)
     INTEGER :: nspin, ir
     !
     nspin = SIZE( r, 2 )
     IF (nspin==1) THEN
!$omp parallel do
        DO ir=1,desc%nr1*desc%nr2*desc%my_nr3p
           r(ir,1)=r(ir,1)+f(ir)
        END DO
!$omp end parallel do
     ELSE IF (nspin==2) THEN
!$omp parallel do
        DO ir=1,desc%nr1*desc%nr2*desc%my_nr3p
           r(ir,1)=r(ir,1)+f(ir)
        END DO
!$omp end parallel do
!$omp parallel do
        DO ir=1,desc%nr1*desc%nr2*desc%my_nr3p
           r(ir,2)=r(ir,2)+f(ir)
        END DO
!$omp end parallel do
     END IF
  END SUBROUTINE
  !
  !---------------------------------------------------------------------
  SUBROUTINE c2psi_gamma( desc, psi, c, ca )
     !------------------------------------------------------------------
     !! Copy wave-functions from 1D array (c_bgrp) to 3D array (psi) in 
     !! Fourier space - gamma case.
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     !! fft descriptor
     COMPLEX(DP), INTENT(OUT) :: psi(:)
     !! w.f. 3D array in Fourier space
     COMPLEX(DP), INTENT(IN) :: c(:)
     !! stores the Fourier expansion coefficients
     COMPLEX(DP), OPTIONAL, INTENT(IN) :: ca(:)
     COMPLEX(DP), PARAMETER :: ci=(0.0d0,1.0d0)
     INTEGER :: ig
     !
     CALL alloc_nl_pntrs( desc )
     !
     !$acc data present_or_copyin(c) present_or_copyout(psi)
     !
     !$acc kernels
     psi = 0.0d0
     !$acc end kernels
     !
     IF( PRESENT(ca) ) THEN
        !$acc parallel loop present_or_copyin(ca)
        DO ig = 1, desc%ngw
          psi(nlm_d(ig)) = CONJG(c(ig)) + ci * CONJG(ca(ig))
          psi(nl_d(ig)) = c(ig) + ci * ca(ig)
        ENDDO
     ELSE
        !$acc parallel loop
        DO ig = 1, desc%ngw
          psi(nlm_d(ig)) = CONJG(c(ig))
          psi(nl_d(ig)) = c(ig)
        ENDDO
     ENDIF
     !
     !$acc end data
     !
     CALL dealloc_nl_pntrs( desc )
     !
  END SUBROUTINE
  !
  !
#ifdef __CUDA
  !---------------------------------------------------------------------
  SUBROUTINE c2psi_gamma_gpu( desc, psi, c, ca )
     !------------------------------------------------------------------
     !! Provisional gpu double of c2psi_gamma for CPV calls (CPV/src/exx_psi.f90).
     !! To be removed after porting exx_psi to openacc.
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     complex(DP), DEVICE, INTENT(OUT) :: psi(:)
     complex(DP), DEVICE, INTENT(IN) :: c(:)
     complex(DP), DEVICE, OPTIONAL, INTENT(IN) :: ca(:)
     complex(DP), parameter :: ci=(0.0d0,1.0d0)
     integer :: ig
     integer, device, pointer :: nlm_d(:), nl_d(:)
     nlm_d => desc%nlm_d
     nl_d => desc%nl_d
     psi = 0.0d0
     IF( PRESENT(ca) ) THEN
        !$cuf kernel do (1)
        do ig = 1, desc%ngw
           psi( nlm_d( ig ) ) = CONJG( c( ig ) ) + ci * conjg( ca( ig ))
           psi( nl_d( ig ) ) = c( ig ) + ci * ca( ig )
        end do
     ELSE
        !$cuf kernel do (1)
        do ig = 1, desc%ngw
           psi( nlm_d( ig ) ) = CONJG( c( ig ) )
           psi( nl_d( ig ) ) = c( ig )
        end do
     END IF
  END SUBROUTINE
#endif
  !
  !--------------------------------------------------------------------------------
  SUBROUTINE c2psi_k( desc, psi, c, igk, ngk, howmany_set )
     !-----------------------------------------------------------------------------
     !! Copy wave-functions from 1D array (c/evc) ordered according (k+G) index igk 
     !! to 3D array (psi) in Fourier space.
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(IN) :: desc
     !! fft descriptor
     COMPLEX(DP), INTENT(OUT) :: psi(:)
     !! w.f. 3D array in Fourier space
     COMPLEX(DP), INTENT(IN) :: c(:,:)
     !! stores the Fourier expansion coefficients of the wave function
     INTEGER, INTENT(IN) :: igk(:), ngk
     INTEGER, OPTIONAL, INTENT(IN) :: howmany_set(3)
     ! hm_set(1)=howmany ; hm_set(2)=ibnd ; hm_set(3)=ibnd
     !
     INTEGER :: nnr, i, j, ig, howmany
     !
     CALL alloc_nl_pntrs( desc )
     !
     !$acc data present_or_copyin(c,igk) present_or_copyout(psi)
     !
     IF (PRESENT(howmany_set)) THEN
        !
        nnr = desc%nnr
        howmany = howmany_set(1)
        ! == OPTIMIZE HERE == (setting to 0 and setting elements!)
        !$acc kernels
        psi(1:nnr*howmany) = (0.d0,0.d0)
        !$acc end kernels
        !
        !$acc parallel loop collapse(2)
        DO i = 0, howmany-1
          DO j = 1, ngk
            psi(nl_d(igk(j))+i*nnr) = c(j,howmany_set(2)+i)
          ENDDO
        ENDDO
        !
     ELSE
        !$acc kernels
        psi = (0.d0,0.d0)
        !$acc end kernels
        !
#if defined(_OPENACC)
        !$acc parallel loop
#else
        !$omp parallel
        !$omp do
#endif
        DO ig = 1, ngk
          psi(nl_d(igk(ig))) = c(ig,1)
        ENDDO
#if !defined(_OPENACC)
        !$omp end do nowait
        !$omp end parallel
#endif
        !
     ENDIF
     !
     !$acc end data
     !
     CALL dealloc_nl_pntrs( desc )
     !
  END SUBROUTINE c2psi_k
  !
  !
  !-------------------------------------------------------------------------
  SUBROUTINE fftx_oned2threed( desc, psi, c, ca )
     !----------------------------------------------------------------------
     !! Copy charge density from 1D array (c) to 3D array (psi) in Fourier
     !! space. If gpu_args_ is true it assumes the dummy arrays are present
     !! on device (openACC porting).
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(IN) :: desc
     !! fft descriptor
     COMPLEX(DP), INTENT(OUT) :: psi(:)
     !! w.f. 3D array in Fourier space
     COMPLEX(DP), INTENT(IN) :: c(:)
     !! stores the Fourier expansion coefficients
     COMPLEX(DP), OPTIONAL, INTENT(IN) :: ca(:)
     COMPLEX(DP), PARAMETER :: ci=(0.0d0,1.0d0)
     INTEGER :: ig
     !
     CALL alloc_nl_pntrs( desc )
     !
     !$acc data present_or_copyin(c) present_or_copyout(psi)
     !$acc kernels
     psi = (0.0d0,0.d0)
     !$acc end kernels
     IF ( PRESENT(ca) ) THEN
        IF( desc%lgamma ) THEN
           !$acc parallel loop present_or_copyin(ca)
           DO ig = 1, desc%ngm
              psi(nlm_d(ig)) = CONJG(c(ig)) + ci * CONJG(ca(ig))
              psi(nl_d(ig)) = c(ig) + ci * ca(ig)
           ENDDO
        ELSE
           !$acc parallel loop present_or_copyin(ca)
           DO ig = 1, desc%ngm
              psi(nl_d(ig)) = c(ig) + ci * ca(ig)
           ENDDO
        ENDIF
     ELSE
        IF( desc%lgamma ) THEN
           !$acc parallel loop
           DO ig = 1, desc%ngm
              psi(nlm_d(ig)) = CONJG(c(ig))
              psi(nl_d(ig)) = c(ig)
           ENDDO
        ELSE
           !$acc parallel loop
           DO ig = 1, desc%ngm
              psi(nl_d(ig)) = c(ig)
           ENDDO
        ENDIF
     ENDIF
     !$acc end data
     !
     CALL dealloc_nl_pntrs( desc )
     !
  END SUBROUTINE fftx_oned2threed
  !
  !-------------------------------------------------------------------
  SUBROUTINE fftx_add_threed2oned_gamma( desc, vin, vout1, vout2 )
     !----------------------------------------------------------------
     !
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     complex(DP), INTENT(INOUT) :: vout1(:)
     complex(DP), OPTIONAL, INTENT(INOUT) :: vout2(:)
     complex(DP), INTENT(IN) :: vin(:)
     COMPLEX(DP) :: fp, fm
     INTEGER :: ig
     !
     IF( PRESENT( vout2 ) ) THEN
        DO ig=1,desc%ngm
           fp=vin(desc%nl(ig))+vin(desc%nlm(ig))
           fm=vin(desc%nl(ig))-vin(desc%nlm(ig))
           vout1(ig) = vout1(ig) + 0.5d0*CMPLX( DBLE(fp),AIMAG(fm),kind=DP)
           vout2(ig) = vout2(ig) + 0.5d0*CMPLX(AIMAG(fp),-DBLE(fm),kind=DP)
        END DO
     ELSE
        DO ig=1,desc%ngm
           vout1(ig) = vout1(ig) + vin(desc%nl(ig))
        END DO
     END IF
  END SUBROUTINE fftx_add_threed2oned_gamma
  !
  !-----------------------------------------------------------------------------
  SUBROUTINE fftx_threed2oned( desc, vin, vout1, vout2 )
     !--------------------------------------------------------------------------
     !! Copy charge density from 3D array to 1D array in Fourier space.
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(IN) :: desc
     COMPLEX(DP), INTENT(OUT) :: vout1(:)
     COMPLEX(DP), OPTIONAL, INTENT(OUT) :: vout2(:)
     COMPLEX(DP), INTENT(IN) :: vin(:)
     !
     COMPLEX(DP) :: fp, fm
     INTEGER :: ig
     !
     CALL alloc_nl_pntrs( desc )
     !
     !$acc data present_or_copyin(vin) present_or_copyout(vout1)
     !
     IF( PRESENT( vout2 ) ) THEN
        !$acc parallel loop present_or_copyout(vout2)
        DO ig = 1, desc%ngm
           fp = vin(nl_d(ig))+vin(nlm_d(ig))
           fm = vin(nl_d(ig))-vin(nlm_d(ig))
           vout1(ig) = CMPLX(0.5d0,0.d0,kind=DP)*CMPLX( DBLE(fp),AIMAG(fm),kind=DP)
           vout2(ig) = CMPLX(0.5d0,0.d0,kind=DP)*CMPLX(AIMAG(fp),-DBLE(fm),kind=DP)
        ENDDO
     ELSE
        !$acc parallel loop
        DO ig = 1, desc%ngm
           vout1(ig) = vin(nl_d(ig))
        ENDDO
     ENDIF
     !$acc end data
     !
     CALL dealloc_nl_pntrs( desc )
     !
  END SUBROUTINE fftx_threed2oned
  !
  !------------------------------------------------------------
  SUBROUTINE fftx_psi2c_gamma( desc, vin, vout1, vout2 )
     !---------------------------------------------------------
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     complex(DP), INTENT(OUT) :: vout1(:)
     complex(DP), OPTIONAL, INTENT(OUT) :: vout2(:)
     complex(DP), INTENT(IN) :: vin(:)
     COMPLEX(DP) :: fp, fm
     INTEGER :: ig
     !
     IF( PRESENT( vout2 ) ) THEN
        DO ig=1,desc%ngw
           fp=vin(desc%nl(ig))+vin(desc%nlm(ig))
           fm=vin(desc%nl(ig))-vin(desc%nlm(ig))
           vout1(ig) = CMPLX( DBLE(fp),AIMAG(fm),kind=DP)
           vout2(ig) = CMPLX(AIMAG(fp),-DBLE(fm),kind=DP)
        END DO
     ELSE
        DO ig=1,desc%ngw
           vout1(ig) = vin(desc%nl(ig))
        END DO
     END IF
     !
  END SUBROUTINE fftx_psi2c_gamma
  !
  !--------------------------------------------------------------------
  SUBROUTINE fftx_psi2c_gamma_gpu( desc, vin, vout1, vout2 )
     !-----------------------------------------------------------------
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     complex(DP), INTENT(OUT) :: vout1(:)
     complex(DP), OPTIONAL, INTENT(OUT) :: vout2(:)
     complex(DP), INTENT(IN) :: vin(:)
     INTEGER,     POINTER     :: nl(:), nlm(:)
#if defined (__CUDA)
     attributes(DEVICE) :: vout1, vout2, vin, nl, nlm
#endif
     INTEGER :: ig
     !
     nl  => desc%nl_d
     nlm => desc%nlm_d
     IF( PRESENT( vout2 ) ) THEN
!$cuf kernel do(1)
        DO ig=1,desc%ngw
           vout1(ig) = CMPLX( DBLE(vin(nl(ig))+vin(nlm(ig))),AIMAG(vin(nl(ig))-vin(nlm(ig))),kind=DP)
           vout2(ig) = CMPLX(AIMAG(vin(nl(ig))+vin(nlm(ig))),-DBLE(vin(nl(ig))-vin(nlm(ig))),kind=DP)
        END DO
     ELSE
!$cuf kernel do(1)
        DO ig=1,desc%ngw
           vout1(ig) = vin(nl(ig))
        END DO
     END IF
  END SUBROUTINE fftx_psi2c_gamma_gpu
  !
  !--------------------------------------------------------------------
  SUBROUTINE c2psi_gamma_tg( desc, psis, c_bgrp, i, nbsp_bgrp )
     !-----------------------------------------------------------------
     !! Copy all wave-functions of an orbital group from 1D array (c_bgrp)
     !! to 3D array (psi) in Fourier space.
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(IN) :: desc
     COMPLEX(DP), INTENT(OUT) :: psis(:)
     COMPLEX(DP), INTENT(INOUT) :: c_bgrp(:,:)
     INTEGER, INTENT(IN) :: i, nbsp_bgrp
     INTEGER :: eig_offset, eig_index, right_nnr
     !
     !  the i-th column of c_bgrp corresponds to the i-th state (in this band group)
     !
     !  The outer loop goes through i : i + 2*NOGRP to cover
     !  2*NOGRP eigenstates at each iteration
     !
     CALL tg_get_nnr( desc, right_nnr )
     !
!$omp  parallel
!$omp  single
     !
     DO eig_index = 1, 2*fftx_ntgrp(desc), 2
        !
!$omp task default(none) &
!$omp          firstprivate( eig_index, i, nbsp_bgrp, right_nnr ) &
!$omp          private( eig_offset ) &
!$omp          shared( c_bgrp, desc, psis )
        !
        !  here we pack 2*nogrp electronic states in the psis array
        !  note that if nogrp == nproc_bgrp each proc perform a full 3D
        !  fft and the scatter phase is local (without communication)
        !
        !  important: if n is odd => c(*,n+1)=0.
        !
        IF ( (eig_index+i-1)==nbsp_bgrp ) c_bgrp(:,eig_index+i) = 0._DP
        !
        eig_offset = (eig_index-1)/2
        !
        IF ( (i+eig_index-1) <= nbsp_bgrp ) THEN
           !
           !  The  eig_index loop is executed only ONCE when NOGRP=1.
           !
           CALL c2psi_gamma( desc, psis(eig_offset*right_nnr+1: &
                             eig_offset*right_nnr+ right_nnr), &
                             c_bgrp(:,i+eig_index-1), c_bgrp(:,i+eig_index) )
           !
        ENDIF
!$omp end task
        !
     ENDDO
!$omp  end single
!$omp  end parallel
     RETURN
  END SUBROUTINE c2psi_gamma_tg
  !
  !
  !--------------------------------------------------------------------
  SUBROUTINE c2psi_k_tg( desc, psis, c_bgrp, igk, ngk, i, nbsp_bgrp )
     !-----------------------------------------------------------------
     !! Copy all wave-functions of an orbital group from 1D array (c_bgrp)
     !! to 3D array (psi) in Fourier space.
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(IN) :: desc
     COMPLEX(DP), INTENT(OUT) :: psis(:)
     COMPLEX(DP), INTENT(INOUT) :: c_bgrp(:,:)
     INTEGER, INTENT(IN) :: igk(:), ngk, i, nbsp_bgrp
     INTEGER :: right_nnr, idx, j, js, je, numblock, ntgrp
     INTEGER, PARAMETER :: blocksize = 256
     !
     CALL alloc_nl_pntrs( desc )
     !
     CALL tg_get_nnr( desc, right_nnr )
     !
     ntgrp = fftx_ntgrp( desc )
     !
     ! ... compute the number of chuncks
     numblock = (ngk+blocksize-1)/blocksize
     !
     ! ... ntgrp ffts at the same time
     !
     !$omp parallel
     !$omp do collapse(2)
     DO idx = 0, MIN(ntgrp-1,nbsp_bgrp-i)
       DO j = 1, numblock
          js = (j-1)*blocksize+1
          je = MIN(j*blocksize,ngk)
          psis(nl_d(igk(js:je))+right_nnr*idx) = c_bgrp(js:je,idx+i)
       ENDDO
     ENDDO
     !$omp end do nowait
     !$omp end parallel
     !
     CALL dealloc_nl_pntrs( desc )
     !
     RETURN
     !
  END SUBROUTINE c2psi_k_tg
  !
  !
  !----------------------------------------------------------
  SUBROUTINE fft_dist_info( desc, unit )
     !-------------------------------------------------------
     !! Prints infos on fft parallel distribution.
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: unit
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER :: i, j, nr3l
     !
     CALL tg_get_local_nr3( desc, nr3l )
     !
     WRITE( stdout,1000) desc%nr1, desc%nr2, desc%nr3, &
                         desc%nr1, desc%my_nr2p, desc%my_nr3p, &
                         1, desc%nproc2, desc%nproc3
     WRITE( stdout,1010) desc%nr1x, desc%nr2x, desc%nr3x
     WRITE( stdout,1020) desc%nnr
     WRITE( stdout,*) '  Number of x-y planes for each processors: '
     WRITE( stdout, fmt = '( 5("  |",I4,",",I4) )' ) ( ( desc%nr2p(j), &
             desc%nr3p(i), i = 1, desc%nproc3 ), j=1,desc%nproc2 )

     IF ( .not. desc%use_pencil_decomposition ) WRITE( stdout,*) '  Using Slab Decomposition'
     IF (       desc%use_pencil_decomposition ) WRITE( stdout,*) '  Using Pencil Decomposition'
1000  FORMAT(3X, &
         'Global Dimensions   Local  Dimensions   Processor Grid',/,3X, &
         '.X.   .Y.   .Z.     .X.   .Y.   .Z.     .X.   .Y.   .Z.',/, &
         3(1X,I5),2X,3(1X,I5),2X,3(1X,I5) )
1010  FORMAT(3X, 'Array leading dimensions ( nr1x, nr2x, nr3x )   = ', 3(1X,I5))
1020  FORMAT(3X, 'Local number of cell to store the grid ( nrxx ) = ', 1X, I9 )
     RETURN
  END SUBROUTINE fft_dist_info
  !
  !----------------------------------------------------------------
  SUBROUTINE alloc_nl_pntrs( desc )
    !------------------------------------------------------------
    !! Workaround to use nl and nlm arrays with or without openacc+cuda
    !! - allocation.
    IMPLICIT NONE
    TYPE(fft_type_descriptor), INTENT(IN) :: desc
#if defined(__CUDA) && defined(_OPENACC)
    nl_d  => desc%nl_d
    nlm_d => desc%nlm_d
#else
    ALLOCATE( nl_d(desc%ngm) )
    nl_d = desc%nl
    IF ( desc%lgamma ) THEN
      ALLOCATE( nlm_d(desc%ngm) )
      nlm_d = desc%nlm
    ENDIF
    !$acc enter data copyin( nl_d, nlm_d )
#endif
    !
  END SUBROUTINE alloc_nl_pntrs
  !
  !----------------------------------------------------------------
  SUBROUTINE dealloc_nl_pntrs( desc )
    !------------------------------------------------------------
    !! Workaround to use nl and nlm arrays with or without openacc+cuda
    !! - deallocation.
    IMPLICIT NONE
    TYPE(fft_type_descriptor), INTENT(IN) :: desc
#if !defined(__CUDA) || !defined(_OPENACC)
    DEALLOCATE( nl_d )
    IF ( desc%lgamma ) DEALLOCATE( nlm_d )
    !$acc exit data delete( nl_d, nlm_d )
#endif
    !
  END SUBROUTINE dealloc_nl_pntrs
  !
  !
END MODULE fft_helper_subroutines
