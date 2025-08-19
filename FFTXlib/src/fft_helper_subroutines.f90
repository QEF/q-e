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
  ! ... Used in CP too
  PUBLIC :: fftx_add_threed2oned_gamma, fftx_psi2c_gamma, fftx_c2psi_gamma,   &
            fftx_c2psi_gamma_tg, fftx_c2psi_k, fftx_c2psi_k_tg, fftx_psi2c_k, &
            fftx_psi2c_gamma_tg, fftx_psi2c_k_tg, fftx_add_field
  PUBLIC :: fft_dist_info
  ! ... Used only in CP+EXX
  PUBLIC :: fftx_tgcomm
  !
#if defined(__CUDA)
  PUBLIC :: fftx_psi2c_gamma_gpu, fftx_c2psi_gamma_gpu
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
  SUBROUTINE fftx_c2psi_gamma( desc, psi, c, ca, howmany_set )
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
     COMPLEX(DP), INTENT(IN) :: c(:,:)
     !! stores the Fourier expansion coefficients
     COMPLEX(DP), OPTIONAL, INTENT(IN) :: ca(:)
     INTEGER, OPTIONAL, INTENT(IN) :: howmany_set(2)
     ! howmany_set(1)=group_size ; howmany_set(2)=npw
     !
     COMPLEX(DP), PARAMETER :: ci=(0.0d0,1.0d0)
     INTEGER :: ig, idx, n, v_siz, pack_size, remainder, howmany, &
                group_size
     !
     !$acc data present_or_copyin(c,desc) present_or_copyout(psi)
     !$acc data present_or_copyin(desc%nl,desc%nlm)
     !
     IF (PRESENT(howmany_set)) THEN
       !
       ! ... FFT batching strategy is defined here:
       ! ...  the buffer in psi_c can contain 2*many_fft bands.
       ! ...  * group_size: is the number of bands that will be transformed.
       ! ...  * pack_size:  is the number of slots in psi_c used to store
       ! ...                the (couples) of bands
       ! ...  * remainder:  can be 1 or 0, if 1, a spare band should be added
       ! ...                in the the first slot of psi_c not occupied by
       ! ...                the couples of bands.
       !
       group_size= howmany_set(1)
       n = howmany_set(2)
       v_siz = desc%nnr
       pack_size = (group_size/2) ! This is FLOOR(group_size/2)
       remainder = group_size - 2*pack_size
       howmany   = pack_size + remainder
       !
       !$acc kernels
       psi(1:desc%nnr*howmany) = (0.d0,0.d0)
       !$acc end kernels
       !
       ! ... two ffts at the same time (remember, v_siz = dffts%nnr)
       IF ( pack_size > 0 ) THEN
          !$acc parallel loop collapse(2)
          DO idx = 0, pack_size-1
             DO ig = 1, n
                psi(desc%nl(ig) + idx*v_siz) = c(ig,2*idx+1) + (0.d0,1.d0)*c(ig,2*idx+2)
                psi(desc%nlm(ig) + idx*v_siz) = CONJG(c(ig,2*idx+1) - (0.d0,1.d0)*c(ig,2*idx+2))
             ENDDO
          ENDDO
       ENDIF
       !
       IF (remainder > 0) THEN
          !$acc parallel loop
          DO ig = 1, n
             psi(desc%nl(ig) + pack_size*v_siz) = c(ig,group_size)
             psi(desc%nlm(ig) + pack_size*v_siz) = CONJG(c(ig,group_size))
          ENDDO
       ENDIF
       !
     ELSE
       !
       n = desc%ngw
       !
       !$acc kernels
       psi = 0.0d0
       !$acc end kernels
       !
       IF( PRESENT(ca) ) THEN
          !$acc parallel loop present_or_copyin(ca)
          DO ig = 1, n
            psi(desc%nlm(ig)) = CONJG(c(ig,1)) + ci * CONJG(ca(ig))
            psi(desc%nl(ig)) = c(ig,1) + ci * ca(ig)
          ENDDO
       ELSE
          !$acc parallel loop
          DO ig = 1, n
            psi(desc%nlm(ig)) = CONJG(c(ig,1))
            psi(desc%nl(ig)) = c(ig,1)
          ENDDO
       ENDIF
       !
     ENDIF
     !
     !$acc end data
     !$acc end data
     !
  END SUBROUTINE fftx_c2psi_gamma
  !
  !
#ifdef __CUDA
  !---------------------------------------------------------------------
  SUBROUTINE fftx_c2psi_gamma_gpu( desc, psi, c, ca )
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
  END SUBROUTINE fftx_c2psi_gamma_gpu
#endif
  !
  !--------------------------------------------------------------------------------
  SUBROUTINE fftx_c2psi_k( desc, psi, c, igk, ngk, howmany )
     !-----------------------------------------------------------------------------
     !! Copy wave-functions from 1D array (c/evc) ordered according (k+G) index igk 
     !! to 3D array (psi) in Fourier space.
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(IN) :: desc
     !! FFT descriptor
     COMPLEX(DP), INTENT(OUT) :: psi(:)
     !! w.f. 3D array in Fourier space
     COMPLEX(DP), INTENT(IN) :: c(:,:)
     !! stores the Fourier expansion coefficients of the wave function
     INTEGER, INTENT(IN) :: igk(:)
     !! index of G corresponding to a given index of k+G
     INTEGER, INTENT(IN) :: ngk
     !! size of c(:,1) or 
     INTEGER, OPTIONAL, INTENT(IN) :: howmany
     !! 
     !
     INTEGER :: nnr, i, j, ig
     !
     !$acc data present_or_copyin(c,igk,desc) present_or_copyout(psi)
     !$acc data present_or_copyin(desc%nl)
     !
     IF (PRESENT(howmany)) THEN
        !
        nnr = desc%nnr
        ! == OPTIMIZE HERE == (setting to 0 and setting elements!)
        !$acc kernels
        psi(1:nnr*howmany) = (0.d0,0.d0)
        !$acc end kernels
        !
        !$acc parallel loop collapse(2)
        DO i = 0, howmany-1
          DO j = 1, ngk
            psi(desc%nl(igk(j))+i*nnr) = c(j,i+1)
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
          psi(desc%nl(igk(ig))) = c(ig,1)
        ENDDO
#if !defined(_OPENACC)
        !$omp end do nowait
        !$omp end parallel
#endif
        !
     ENDIF
     !
     !$acc end data
     !$acc end data
     !
  END SUBROUTINE fftx_c2psi_k
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
     INTEGER :: ig, desc_ngm
     !
     desc_ngm = desc%ngm
     !
     !$acc data present_or_copyin(c,desc) present_or_copyout(psi)
     !$acc data present_or_copyin(desc%nl,desc%nlm)
     !$acc kernels
     psi = (0.0d0,0.d0)
     !$acc end kernels
     IF ( PRESENT(ca) ) THEN
        IF( desc%lgamma ) THEN
           !$acc parallel loop present_or_copyin(ca)
           DO ig = 1, desc_ngm
              psi(desc%nlm(ig)) = CONJG(c(ig)) + ci * CONJG(ca(ig))
              psi(desc%nl(ig)) = c(ig) + ci * ca(ig)
           ENDDO
        ELSE
           !$acc parallel loop present_or_copyin(ca)
           DO ig = 1, desc_ngm
              psi(desc%nl(ig)) = c(ig) + ci * ca(ig)
           ENDDO
        ENDIF
     ELSE
        IF( desc%lgamma ) THEN
           !$acc parallel loop
           DO ig = 1, desc_ngm
              psi(desc%nlm(ig)) = CONJG(c(ig))
              psi(desc%nl(ig)) = c(ig)
           ENDDO
        ELSE
           !$acc parallel loop
           DO ig = 1, desc_ngm
              psi(desc%nl(ig)) = c(ig)
           ENDDO
        ENDIF
     ENDIF
     !$acc end data
     !$acc end data
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
  SUBROUTINE fftx_threed2oned( desc, vin, vout1, vout2, iigs )
     !--------------------------------------------------------------------------
     !! Copy charge density from 3D array to 1D array in Fourier space.
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(IN) :: desc
     COMPLEX(DP), INTENT(OUT) :: vout1(:)
     COMPLEX(DP), OPTIONAL, INTENT(OUT) :: vout2(:)
     COMPLEX(DP), INTENT(IN) :: vin(:)
     INTEGER, INTENT(IN), OPTIONAL :: iigs
     !
     COMPLEX(DP) :: fp, fm
     INTEGER :: ig, iigs_, nng
     !
     iigs_ = 0
     nng = desc%ngm
     IF (PRESENT(iigs)) THEN
        iigs_ = iigs-1
        nng = SIZE(vout1)
     ENDIF
     !
     !$acc data present_or_copyin(vin,desc) present_or_copyout(vout1)
     !$acc data present_or_copyin(desc%nl,desc%nlm)
     !
     IF( PRESENT( vout2 ) ) THEN
        !$acc parallel loop present_or_copyout(vout2)
        DO ig = 1, nng
           fp = vin(desc%nl(iigs_+ig))+vin(desc%nlm(iigs_+ig))
           fm = vin(desc%nl(iigs_+ig))-vin(desc%nlm(iigs_+ig))
           vout1(ig) = CMPLX(0.5d0,0.d0,kind=DP)*CMPLX( DBLE(fp),AIMAG(fm),kind=DP)
           vout2(ig) = CMPLX(0.5d0,0.d0,kind=DP)*CMPLX(AIMAG(fp),-DBLE(fm),kind=DP)
        ENDDO
     ELSE
        !$acc parallel loop
        DO ig = 1, nng
           vout1(ig) = vin(desc%nl(iigs_+ig))
        ENDDO
     ENDIF
     !
     !$acc end data
     !$acc end data
     !
  END SUBROUTINE fftx_threed2oned
  !
  !------------------------------------------------------------
  SUBROUTINE fftx_psi2c_gamma( desc, vin, vout1, vout2, howmany_set )
     !---------------------------------------------------------
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(IN) :: desc
     COMPLEX(DP), INTENT(OUT) :: vout1(:,:)
     COMPLEX(DP), OPTIONAL, INTENT(OUT) :: vout2(:)
     COMPLEX(DP), INTENT(IN) :: vin(:)
     INTEGER, OPTIONAL, INTENT(IN) :: howmany_set(2)
     !
     COMPLEX(DP) :: fp, fm
     INTEGER :: ig, idx, n, v_siz, pack_size, remainder, howmany, &
                group_size, ioff
     !
     !$acc data present_or_copyin(vin,desc) present_or_copyout(vout1)
     !$acc data present_or_copyin(desc%nl,desc%nlm)
     !
     IF (PRESENT(howmany_set)) THEN
       !
       group_size = howmany_set(1)
       n = howmany_set(2)
       v_siz = desc%nnr
       pack_size = (group_size/2)
       remainder = group_size - 2*pack_size
       howmany = pack_size + remainder
       !
       IF ( pack_size > 0 ) THEN
         ! ... two ffts at the same time
         !$acc parallel loop collapse(2)
         DO idx = 0, pack_size-1
           DO ig = 1, n
             ioff = idx*v_siz
             fp = (vin(ioff+desc%nl(ig)) + vin(ioff+desc%nlm(ig)))*0.5d0
             fm = (vin(ioff+desc%nl(ig)) - vin(ioff+desc%nlm(ig)))*0.5d0
             vout1(ig,idx*2+1)   = CMPLX(DBLE(fp),AIMAG(fm),KIND=DP)
             vout1(ig,idx*2+2) = CMPLX(AIMAG(fp),-DBLE(fm),KIND=DP)
           ENDDO
         ENDDO
       ENDIF
       IF (remainder > 0) THEN
         !$acc parallel loop
         DO ig = 1, n
           vout1(ig,group_size) = vin(pack_size*v_siz+desc%nl(ig))
         ENDDO
       ENDIF
       !
     ELSE
       !
       n = desc%ngw
       !
       IF( PRESENT(vout2) ) THEN
          !$acc parallel loop present_or_copyout(vout2)
          DO ig = 1, n
             fp = vin(desc%nl(ig))+vin(desc%nlm(ig))
             fm = vin(desc%nl(ig))-vin(desc%nlm(ig))
             vout1(ig,1) = CMPLX( DBLE(fp),AIMAG(fm),kind=DP)
             vout2(ig) = CMPLX(AIMAG(fp),-DBLE(fm),kind=DP)
          ENDDO
       ELSE
          !$acc parallel loop
          DO ig = 1, n
             vout1(ig,1) = vin(desc%nl(ig))
          ENDDO
       ENDIF
       !
     ENDIF
     !
     !$acc end data
     !$acc end data
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
     !
     INTEGER, POINTER     :: nl(:), nlm(:)
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
  !------------------------------------------------------------
  SUBROUTINE fftx_psi2c_k( desc, vin, vout, igk, howmany_set )
     !---------------------------------------------------------
     !
     USE fft_types,      ONLY : fft_type_descriptor
     !
     TYPE(fft_type_descriptor), INTENT(IN) :: desc
     COMPLEX(DP), INTENT(IN) :: vin(:)
     COMPLEX(DP), INTENT(OUT) :: vout(:,:)
     INTEGER, INTENT(IN) :: igk(:)
     INTEGER, OPTIONAL, INTENT(IN) :: howmany_set(2)
     !
     INTEGER :: ig, igmax, idx, n, group_size, v_siz
     !
     !$acc data present_or_copyin(vin,igk,desc) present_or_copyout(vout)
     !$acc data present_or_copyin(desc%nl)
     !
     IF (PRESENT(howmany_set)) THEN
        !
        group_size = howmany_set(1)
        n = howmany_set(2)
        v_siz = desc%nnr
        !
        !$acc parallel loop collapse(2)
        DO idx = 0, group_size-1
           DO ig = 1, n
              vout(ig,idx+1) = vin(idx*v_siz+desc%nl(igk(ig)))
           ENDDO
        ENDDO
        !
     ELSE
        !
        igmax = MIN(desc%ngw,SIZE(vout(:,1)))
        !$acc parallel loop
        DO ig = 1, igmax
          vout(ig,1) = vin(desc%nl(igk(ig)))
        ENDDO
        !
     ENDIF
     !
     !$acc end data
     !$acc end data
     !
     RETURN
     !
  END SUBROUTINE fftx_psi2c_k
  !
  !--------------------------------------------------------------------
  SUBROUTINE fftx_c2psi_gamma_tg( desc, psis, c_bgrp, n, dbnd )
     !-----------------------------------------------------------------
     !! Copy all wave-functions of an orbital group from 1D array (c_bgrp)
     !! to 3D array (psi) in Fourier space.
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(IN) :: desc
     COMPLEX(DP), INTENT(OUT) :: psis(:)
     COMPLEX(DP), INTENT(IN) :: c_bgrp(:,:)
     INTEGER, INTENT(IN) :: n, dbnd
     !
     INTEGER :: eig_offset, eig_index, right_nnr, ig, ib, ieg
     COMPLEX(DP), PARAMETER :: ci=(0.0d0,1.0d0)
     !
     ! ... the i-th column of c_bgrp corresponds to the i-th state (in this band group)
     !
     ! ... The outer loop goes through i : i + 2*NOGRP to cover
     ! ... 2*NOGRP eigenstates at each iteration
     !
     right_nnr = desc%nnr
     !
#if defined(_OPENACC)
     !$acc data present_or_copyin(c_bgrp,desc) present_or_copyout(psis)
     !$acc data present_or_copyin(desc%nl,desc%nlm)
#else
     !$omp  parallel
     !$omp  single
#endif
     !
     DO eig_index = 1, 2*fftx_ntgrp(desc), 2
        !
#if !defined(_OPENACC)
        !$omp task default(none) &
        !$omp          firstprivate( eig_index, dbnd, right_nnr ) &
        !$omp          private( eig_offset, ib, ieg, ig ) &
        !$omp          shared( c_bgrp, desc, psis, n )
#endif
        !
        ! ... here we pack 2*nogrp electronic states in the psis array
        ! ... note that if nogrp == nproc_bgrp each proc perform a full 3D
        ! ... fft and the scatter phase is local (without communication)
        !
        ! ... important: if n is odd => c(*,n+1)=0.
        !
        eig_offset = (eig_index-1)/2
        !
        ib = eig_offset*right_nnr
        ieg = eig_index
        !
        ! ... The eig_index loop is executed only ONCE when NOGRP=1.
        IF ( ieg < dbnd ) THEN
           !$acc parallel loop
           DO ig = 1, n !desc%ngw
             psis(ib+desc%nlm(ig)) = CONJG(c_bgrp(ig,ieg)) + ci * CONJG(c_bgrp(ig,ieg+1))
             psis(ib+desc%nl(ig)) = c_bgrp(ig,ieg) + ci * c_bgrp(ig,ieg+1)
           ENDDO
        ELSEIF ( ieg == dbnd ) THEN
           ! ... important: if n is odd => c(*,n+1)=0.
           !$acc parallel loop
           DO ig = 1, n !desc%ngw
             psis(ib+desc%nlm(ig)) = CONJG(c_bgrp(ig,ieg))
             psis(ib+desc%nl(ig)) = c_bgrp(ig,ieg)
           ENDDO
        ENDIF
#if !defined(_OPENACC)
        !$omp end task
#endif
        !
     ENDDO
#if defined(_OPENACC)
     !$acc end data
     !$acc end data
#else
     !$omp end single
     !$omp end parallel
#endif
     !
     RETURN
     !
  END SUBROUTINE fftx_c2psi_gamma_tg
  !
  !
  !--------------------------------------------------------------------
  SUBROUTINE fftx_c2psi_k_tg( desc, psis, c_bgrp, igk, ngk, dbnd )
     !-----------------------------------------------------------------
     !! Copy all wave-functions of an orbital group from 1D array (c_bgrp)
     !! to 3D array (psi) in Fourier space.
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(IN) :: desc
     COMPLEX(DP), INTENT(OUT) :: psis(:)
     COMPLEX(DP), INTENT(INOUT) :: c_bgrp(:,:)
     INTEGER, INTENT(IN) :: igk(:), ngk, dbnd
     !
     INTEGER :: right_nnr, idx, j, js, je, numblock, ntgrp
     INTEGER, PARAMETER :: blocksize = 256
     !
     right_nnr = desc%nnr
     !
#if defined(_OPENACC)
     !$acc data present_or_copyin(c_bgrp,igk,desc) present_or_copyout(psis)
     !$acc data present_or_copyin(desc%nl)
#endif
     !
     ntgrp = fftx_ntgrp( desc )
     !
     ! ... compute the number of chuncks
     numblock = (ngk+blocksize-1)/blocksize
     !
     ! ... ntgrp ffts at the same time
     !
#if !defined(_OPENACC)
     !$omp parallel do collapse(2) private(js,je)
#else
     !$acc parallel loop collapse(2)
#endif
     DO idx = 0, MIN(ntgrp-1,dbnd-1)
       DO j = 1, numblock
          js = (j-1)*blocksize+1
          je = MIN(j*blocksize,ngk)
          psis(desc%nl(igk(js:je))+right_nnr*idx) = c_bgrp(js:je,idx+1)
       ENDDO
     ENDDO
#if !defined(_OPENACC)
     !$omp end parallel do
#else
     !$acc end data
     !$acc end data
#endif
     !
     RETURN
     !
  END SUBROUTINE fftx_c2psi_k_tg
  !
  !
  !--------------------------------------------------------------------
  SUBROUTINE fftx_psi2c_gamma_tg( desc, vin, vout, n, dbnd )
     !-----------------------------------------------------------------
     !! Copy all wave-functions of an orbital group from 3D array (psi)
     !! in Fourier space to 1D array (c_bgrp). Gamma case.
     !
     USE fft_types,   ONLY : fft_type_descriptor
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(IN) :: desc
     COMPLEX(DP), INTENT(IN) :: vin(:)
     COMPLEX(DP), INTENT(OUT) :: vout(:,:)
     INTEGER, INTENT(IN) :: n, dbnd
     !
     INTEGER :: right_inc, idx, j, ioff
     COMPLEX(DP) :: fp, fm
     !
     ioff = 0
     CALL tg_get_recip_inc( desc, right_inc )
     !
     !$acc data present_or_copyin(vin,desc) present_or_copyout(vout)
     !$acc data present_or_copyin(desc%nl,desc%nlm)
     !
     DO idx = 1, 2*fftx_ntgrp(desc), 2
       !
       IF ( idx<dbnd ) THEN
         !$acc parallel loop
         DO j = 1, n
           fp = ( vin(desc%nl(j)+ioff) + vin(desc%nlm(j)+ioff) )
           fm = ( vin(desc%nl(j)+ioff) - vin(desc%nlm(j)+ioff) )
           vout(j,idx)   = CMPLX( DBLE(fp),AIMAG(fm),KIND=DP)
           vout(j,idx+1) = CMPLX(AIMAG(fp),-DBLE(fm),KIND=DP)
         ENDDO
       ELSEIF ( idx==dbnd ) THEN
         !$acc parallel loop
         DO j = 1, n
            vout(j,idx) = vin(desc%nl(j)+ioff)
         ENDDO
       ENDIF
       !
       ioff = ioff + right_inc
       !
     ENDDO
     !
     !$acc end data
     !$acc end data
     !
     RETURN
     !
  END SUBROUTINE fftx_psi2c_gamma_tg
  !
  !
  !--------------------------------------------------------------------
  SUBROUTINE fftx_psi2c_k_tg( desc, vin, vout, igk, n, dbnd )
     !-----------------------------------------------------------------
     !! Copy all wave-functions of an orbital group from 3D array (psi)
     !! in Fourier space to 1D array (c_bgrp).
     !
     USE fft_types,   ONLY : fft_type_descriptor
     !
     IMPLICIT NONE
     !
     TYPE(fft_type_descriptor), INTENT(IN) :: desc
     COMPLEX(DP), INTENT(IN) :: vin(:)
     COMPLEX(DP), INTENT(OUT) :: vout(:,:)
     INTEGER, INTENT(IN) :: igk(:), n, dbnd
     !
     INTEGER :: right_inc, idx, j, iin, numblock
     INTEGER, PARAMETER :: blocksize = 256
     !
     CALL tg_get_recip_inc( desc, right_inc )
     !
     numblock = (n+blocksize-1)/blocksize
     !
     !$omp parallel do collapse(2)
     DO idx = 0, MIN(fftx_ntgrp(desc)-1, dbnd-1)
       DO j = 1, numblock
          DO iin = (j-1)*blocksize+1, MIN(j*blocksize,n)
            vout(iin,1+idx) = vin(desc%nl(igk(iin))+right_inc*idx)
          ENDDO
       ENDDO
     ENDDO
     !$omp end parallel do
     !
     RETURN
     !
  END SUBROUTINE fftx_psi2c_k_tg
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
  !
END MODULE fft_helper_subroutines
