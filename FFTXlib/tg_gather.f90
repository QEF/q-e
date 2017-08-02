!----------------------------------------------------------------------------------------------------------------
!-real version 
SUBROUTINE tg_gather( dffts, v, tg_v )
  !
  USE fft_param
  USE fft_types,      ONLY : fft_type_descriptor

  IMPLICIT NONE

  TYPE(fft_type_descriptor), INTENT(in) :: dffts
  REAL(DP), INTENT(IN)  :: v(dffts%nnr)
  REAL(DP), INTENT(OUT) :: tg_v(dffts%nnr_tg)

  INTEGER :: nxyp, ir3, off, tg_off
  INTEGER :: i, nsiz, ierr

  nxyp   = dffts%nr1x*dffts%my_nr2p
  !
  !  The potential in v is distributed so that each Z-plane is shared among nproc2 processors.
  !  We collect the data of whole planes in tg_v to be used with task group distributed wfcs.
  !
  tg_v(:) = (0.d0,0.d0)
  do ir3 =1, dffts%my_nr3p
     off    = dffts%nr1x*dffts%my_nr2p*(ir3-1)
     tg_off = dffts%nr1x*dffts%nr2x   *(ir3-1) + dffts%nr1x*dffts%my_i0r2p
     tg_v(tg_off+1:tg_off+nxyp) = v(off+1:off+nxyp)
  end do
  !write (6,*) ' tg_v ', dffts%my_i0r2p, dffts%my_nr2p
  !write (6,'(20f12.7)') (v(dffts%my_i0r2p+i+dffts%nr1x*(i-1)), i=1,dffts%my_nr2p)
  !write (6,'(20f12.7)') (tg_v(i+dffts%nr1x*(i-1)), i=1,dffts%nr2x)
#if defined(__MPI)
!used to be   CALL mp_sum(tg_v, dffts%comm2 )
  nsiz =dffts%nnr_tg
  CALL MPI_ALLREDUCE( MPI_IN_PLACE, tg_v, nsiz, MPI_DOUBLE_PRECISION, MPI_SUM, dffts%comm2, ierr )
  IF( ierr /= 0 ) CALL fftx_error__( ' tg_gather ', ' MPI_ALLREDUCE ', abs( ierr ) )
!- could be done (more efficintly?) with an ALLgatherv but the loigc of the ALLREDUCE is simpler
!  CALL MPI_Allgatherv( v(1), nsiz, MPI_DOUBLE_PRECISION, &
!        tg_v(1), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, dffts%comm2, IERR)
!  IF( ierr /= 0 ) CALL fftx_error__( ' tg_gather ', ' MPI_Allgatherv ', abs( ierr ) )
#endif
  !write (6,'(20f12.7)') (tg_v(i+dffts%nr1x*(i-1)), i=1,dffts%nr1x)
  RETURN
END SUBROUTINE tg_gather

!-complex version of previous routine
SUBROUTINE tg_cgather( dffts, v, tg_v )
  !
  USE fft_param
  USE fft_types,      ONLY : fft_type_descriptor

  IMPLICIT NONE

  TYPE(fft_type_descriptor), INTENT(in) :: dffts
  COMPLEX(DP), INTENT(IN)  :: v(dffts%nnr)
  COMPLEX(DP), INTENT(OUT) :: tg_v(dffts%nnr_tg)

  INTEGER :: nxyp, ir3, off, tg_off
  INTEGER :: i, nsiz, ierr

  nxyp   = dffts%nr1x*dffts%my_nr2p
  !
  !  The potential in v is distributed so that each Z-plane is shared among nproc2 processors.
  !  We collect the data of whole planes in tg_v to be used with task group distributed wfcs.
  !
  tg_v(:) = (0.d0,0.d0)
  do ir3 =1, dffts%my_nr3p
     off    = dffts%nr1x*dffts%my_nr2p*(ir3-1)
     tg_off = dffts%nr1x*dffts%nr2x   *(ir3-1) + dffts%nr1x*dffts%my_i0r2p
     tg_v(tg_off+1:tg_off+nxyp) = v(off+1:off+nxyp)
  end do
  !write (6,*) ' tg_v ', dffts%my_i0r2p, dffts%my_nr2p
  !write (6,'(20f12.7)') (v(dffts%my_i0r2p+i+dffts%nr1x*(i-1)), i=1,dffts%my_nr2p)
  !write (6,'(20f12.7)') (tg_v(i+dffts%nr1x*(i-1)), i=1,dffts%nr2x)
#if defined(__MPI)
!used to be   CALL mp_sum(tg_v, dffts%comm2 )
  nsiz =2 * dffts%nnr_tg
  CALL MPI_ALLREDUCE( MPI_IN_PLACE, tg_v, nsiz, MPI_DOUBLE_PRECISION, MPI_SUM, dffts%comm2, ierr )
  IF( ierr /= 0 ) CALL fftx_error__( ' tg_gather ', ' MPI_ALLREDUCE ', abs( ierr ) )
!- could be done (more efficintly?) with an ALLgatherv but the loigc of the ALLREDUCE is simpler
!  CALL MPI_Allgatherv( v(1), nsiz, MPI_DOUBLE_PRECISION, &
!        tg_v(1), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, dffts%comm2, IERR)
!  IF( ierr /= 0 ) CALL fftx_error__( ' tg_gather ', ' MPI_Allgatherv ', abs( ierr ) )
#endif
  !write (6,'(20f12.7)') (tg_v(i+dffts%nr1x*(i-1)), i=1,dffts%nr1x)
  RETURN
END SUBROUTINE tg_cgather

