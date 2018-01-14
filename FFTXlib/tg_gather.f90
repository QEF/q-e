
SUBROUTINE tg_gather( dffts, v, tg_v )
   !
   USE fft_param
   USE fft_types,      ONLY : fft_type_descriptor

   ! T.G.
   ! NOGRP:      Number of processors per orbital task group

   IMPLICIT NONE

   TYPE(fft_type_descriptor), INTENT(in) :: dffts
   REAL(DP), INTENT(IN)  :: v(dffts%nnr)
   REAL(DP), INTENT(OUT) :: tg_v(dffts%nnr_tg)


   INTEGER :: nsiz, i, ierr, nsiz_tg
   INTEGER :: recv_cnt( dffts%nogrp ), recv_displ( dffts%nogrp )

   nsiz_tg = dffts%tg_nnr * dffts%nogrp

   IF( size( tg_v ) < nsiz_tg ) &
      CALL fftx_error__( ' tg_gather ', ' tg_v too small ', ( nsiz_tg - size( tg_v ) ) )

   nsiz = dffts%npp( dffts%mype+1 ) * dffts%nr1x * dffts%nr2x

   IF( size( v ) < nsiz ) &
      CALL fftx_error__( ' tg_gather ', ' v too small ',  ( nsiz - size( v ) ) )

   !
   !  The potential in v is distributed across all processors
   !  We need to redistribute it so that it is completely contained in the
   !  processors of an orbital TASK-GROUP
   !
   recv_cnt(1)   = dffts%npp( dffts%nolist(1) + 1 ) * dffts%nr1x * dffts%nr2x
   recv_displ(1) = 0
   DO i = 2, dffts%nogrp
      recv_cnt(i) = dffts%npp( dffts%nolist(i) + 1 ) * dffts%nr1x * dffts%nr2x
      recv_displ(i) = recv_displ(i-1) + recv_cnt(i-1)
   ENDDO

   ! clean only elements that will not be overwritten
   !
   DO i = recv_displ(dffts%nogrp) + recv_cnt( dffts%nogrp ) + 1, size( tg_v )
      tg_v( i ) = 0.0d0
   ENDDO

#if defined(__MPI)

   CALL MPI_Allgatherv( v(1), nsiz, MPI_DOUBLE_PRECISION, &
        tg_v(1), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, dffts%ogrp_comm, IERR)

   IF( ierr /= 0 ) &
      CALL fftx_error__( ' tg_gather ', ' MPI_Allgatherv ', abs( ierr ) )

#endif

END SUBROUTINE tg_gather


SUBROUTINE tg_cgather( dffts, v, tg_v )
   !
   USE fft_param
   USE fft_types,      ONLY : fft_type_descriptor

   ! T.G.
   ! NOGRP:      Number of processors per orbital task group

   IMPLICIT NONE

   TYPE(fft_type_descriptor), INTENT(in) :: dffts
   COMPLEX(DP), INTENT(IN)  :: v(dffts%nnr)
   COMPLEX(DP), INTENT(OUT) :: tg_v(dffts%nnr_tg)

   INTEGER :: nsiz, i, ierr, nsiz_tg
   INTEGER :: recv_cnt( dffts%nogrp ), recv_displ( dffts%nogrp )

   nsiz_tg = dffts%tg_nnr * dffts%nogrp

   IF( size( tg_v ) < nsiz_tg ) &
      CALL fftx_error__( ' tg_gather ', ' tg_v too small ', ( nsiz_tg - size( tg_v ) ) )

   nsiz = dffts%npp( dffts%mype+1 ) * dffts%nr1x * dffts%nr2x

   IF( size( v ) < nsiz ) &
      CALL fftx_error__( ' tg_gather ', ' v too small ',  ( nsiz - size( v ) ) )

   !
   !  The potential in v is distributed across all processors
   !  We need to redistribute it so that it is completely contained in the
   !  processors of an orbital TASK-GROUP
   !
   recv_cnt(1)   = dffts%npp( dffts%nolist(1) + 1 ) * dffts%nr1x * dffts%nr2x
   recv_displ(1) = 0
   DO i = 2, dffts%nogrp
      recv_cnt(i) = dffts%npp( dffts%nolist(i) + 1 ) * dffts%nr1x * dffts%nr2x
      recv_displ(i) = recv_displ(i-1) + recv_cnt(i-1)
   ENDDO

   ! clean only elements that will not be overwritten
   !
   DO i = recv_displ(dffts%nogrp) + recv_cnt( dffts%nogrp ) + 1, size( tg_v )
      tg_v( i ) = (0.0d0,0.0d0)
   ENDDO
   !
   ! The quantities are complex, multiply the cunters by 2 and gather
   ! real numbers
   !
   nsiz = 2 * nsiz
   recv_cnt = 2 * recv_cnt
   recv_displ = 2 * recv_displ
#if defined(__MPI)

   CALL MPI_Allgatherv( v(1), nsiz, MPI_DOUBLE_PRECISION, &
        tg_v(1), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, dffts%ogrp_comm, IERR)

   IF( ierr /= 0 ) &
      CALL fftx_error__( ' tg_cgather ', ' MPI_Allgatherv ', abs( ierr ) )

#endif

END SUBROUTINE tg_cgather
