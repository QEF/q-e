MODULE fft_helper_subroutines

  IMPLICIT NONE
  SAVE

  INTERFACE tg_reduce_rho
    MODULE PROCEDURE tg_reduce_rho_1,tg_reduce_rho_2,tg_reduce_rho_3,tg_reduce_rho_4, &
&                    tg_reduce_rho_5
  END INTERFACE


CONTAINS

  SUBROUTINE tg_reduce_rho_1( rhos, tg_rho_nc, tg_rho, ispin, noncolin, domag, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor

     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(IN) :: ispin
     LOGICAL, INTENT(IN) :: noncolin, domag
     REAL(DP), INTENT(INOUT)  :: tg_rho(:)
     REAL(DP), INTENT(INOUT)  :: tg_rho_nc(:,:)
     REAL(DP), INTENT(OUT) :: rhos(:,:)

     INTEGER :: ierr, ioff, idx, ir3, ir, ipol, ioff_tg, nxyp, npol_

#ifdef __MPI
     IF( noncolin) THEN
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tg_rho_nc, SIZE(tg_rho_nc), MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
     ELSE
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tg_rho, SIZE(tg_rho), MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
     END IF
#endif
     !
     ! copy the charge back to the proper processor location
     !
     nxyp = desc%nr1x * desc%my_nr2p
     IF (noncolin) THEN
         npol_ = 1 ; if (domag) npol_ = 4
     endif
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

  END SUBROUTINE



  SUBROUTINE tg_reduce_rho_2( rhos, tmp_rhos, ispin, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor

     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(IN) :: ispin
     REAL(DP), INTENT(INOUT)  :: tmp_rhos(:)
     REAL(DP), INTENT(OUT) :: rhos(:,:)

     INTEGER :: ierr, ioff, idx, ir3, nxyp, ioff_tg

     IF ( desc%nproc2 > 1 ) THEN
#ifdef __MPI
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
 
  END SUBROUTINE

  SUBROUTINE tg_reduce_rho_3( rhos, tmp_rhos, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor

     TYPE(fft_type_descriptor), INTENT(in) :: desc
     REAL(DP), INTENT(INOUT)  :: tmp_rhos(:,:)
     REAL(DP), INTENT(OUT) :: rhos(:,:)

     INTEGER :: ierr, from, ir3, ioff, nxyp, ioff_tg

     IF ( desc%nproc2 > 1 ) THEN
#ifdef __MPI
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
  END SUBROUTINE


  SUBROUTINE tg_reduce_rho_4( rhos, tmp_rhos, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor

     TYPE(fft_type_descriptor), INTENT(in) :: desc
     COMPLEX(DP), INTENT(INOUT)  :: tmp_rhos(:)
     COMPLEX(DP), INTENT(OUT) :: rhos(:)

     INTEGER :: ierr, from, ir3, ioff, nxyp, ioff_tg

     IF ( desc%nproc2 > 1 ) THEN
#ifdef __MPI
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
        rhos(ioff+1:ioff+nxyp) = rhos(ioff+1:ioff+nxyp) + tmp_rhos(ioff_tg+1:ioff_tg+nxyp)
     END DO
  END SUBROUTINE


  SUBROUTINE tg_reduce_rho_5( rhos, tmp_rhos, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor

     TYPE(fft_type_descriptor), INTENT(in) :: desc
     COMPLEX(DP), INTENT(INOUT)  :: tmp_rhos(:,:)
     COMPLEX(DP), INTENT(OUT) :: rhos(:,:)

     INTEGER :: ierr, from, ir3, ioff, nxyp, ioff_tg

     IF ( desc%nproc2 > 1 ) THEN
#ifdef __MPI
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
  END SUBROUTINE



  SUBROUTINE tg_get_nnr( desc, right_nnr )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: right_nnr
     right_nnr = desc%nnr
  END SUBROUTINE

  SUBROUTINE tg_get_local_nr3( desc, val )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: val
     val = desc%my_nr3p
  END SUBROUTINE

  SUBROUTINE tg_get_group_nr3( desc, val )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: val
     val = desc%my_nr3p
  END SUBROUTINE

  SUBROUTINE tg_get_recip_inc( desc, val )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: val
     val = desc%nnr
  END SUBROUTINE

  PURE FUNCTION fftx_ntgrp( desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     INTEGER :: fftx_ntgrp
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     fftx_ntgrp = desc%nproc2
  END FUNCTION

  PURE FUNCTION fftx_tgpe( desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     INTEGER :: fftx_tgpe
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     fftx_tgpe = desc%mype2
  END FUNCTION

  SUBROUTINE fftx_add_field( r, f, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     REAL(DP), INTENT(INOUT)  :: r(:,:)
     REAL(DP), INTENT(IN) :: f(:)
     INTEGER :: nspin, ir
     nspin = SIZE( r, 2 )
     IF (nspin.EQ.1) THEN
!$omp parallel do
        DO ir=1,desc%nr1*desc%nr2*desc%my_nr3p
           r(ir,1)=r(ir,1)+f(ir)
        END DO
!$omp end parallel do
     ELSE IF (nspin.EQ.2) THEN
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


END MODULE fft_helper_subroutines
