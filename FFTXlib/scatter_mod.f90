!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
! FFT base Module.
! Written by Carlo Cavazzoni, modified by Paolo Giannozzi
! Rewritten by Stefano de Gironcoli
!----------------------------------------------------------------------
!
!=----------------------------------------------------------------------=!
   MODULE scatter_mod
!=----------------------------------------------------------------------=!

        USE fft_types, ONLY: fft_type_descriptor
        USE fft_param

        IMPLICIT NONE

        INTERFACE gather_grid
           MODULE PROCEDURE gather_real_grid, gather_complex_grid
        END INTERFACE

        INTERFACE scatter_grid
           MODULE PROCEDURE scatter_real_grid, scatter_complex_grid
        END INTERFACE

        SAVE

        PRIVATE

        PUBLIC :: gather_grid, scatter_grid
        PUBLIC :: cgather_sym, cgather_sym_many, cscatter_sym_many

!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!
!
!
!----------------------------------------------------------------------------
SUBROUTINE gather_real_grid ( dfft, f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers a distributed real-space FFT grid to dfft%root, that is,
  ! ... the first processor of input descriptor dfft - version for real arrays
  !
  ! ... REAL*8  f_in  = distributed variable (dfft%nnr)
  ! ... REAL*8  f_out = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(in) :: f_in (:)
  REAL(DP), INTENT(inout):: f_out(:)
  TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info, offset_in, offset_aux, ir3
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfft%nproc-1), recvcount(0:dfft%nproc-1)
  REAL(DP), ALLOCATABLE ::  f_aux(:)
  !
  IF( size( f_in ) < dfft%nnr ) &
     CALL fftx_error__( ' gather_real_grid ', ' f_in too small ', dfft%nnr-size( f_in ) )
  !
  CALL start_clock( 'rgather_grid' )

  !write (6,*) 'gather grid ok 0 ', dfft%nproc, dfft%nproc2, dfft%nproc3
  ALLOCATE ( f_aux(dfft%nr1x * dfft%nr2x * dfft%my_nr3p) )
  ! 1) gather within the comm2 communicator
  displs = 0
  DO proc =0, ( dfft%nproc2 -1 )
     recvcount(proc) = dfft%nr1x * dfft%nr2p(proc+1)
     if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
  end do
  offset_in = 1; offset_aux = 1
  do ir3 = 1, dfft%my_nr3p
    !write (6,*) 'gather grid ok 1 ir3=', ir3
     info = 0
     CALL MPI_GATHERV( f_in(offset_in) , recvcount(dfft%mype2), MPI_DOUBLE_PRECISION, &
                       f_aux(offset_aux), recvcount, displs, MPI_DOUBLE_PRECISION, dfft%root,      &
                       dfft%comm2, info )
     CALL fftx_error__( 'gather_real_grid', 'info<>0', info )
    !write (6,*) 'gather grid ok 2 ir3=', ir3
     offset_in  = offset_in + dfft%nr1x * dfft%my_nr2p
     offset_aux = offset_aux + dfft%nr1x * dfft%nr2
  end do

  ! 2) gather within the comm3 communicator
  displs = 0
  DO proc = 0, ( dfft%nproc3 - 1 )
     recvcount(proc) = dfft%nr1x * dfft%nr2x * dfft%nr3p(proc+1)
     if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
  ENDDO
  info = 0
  !write (6,*) 'gather grid ok 3'
  CALL MPI_GATHERV( f_aux, recvcount(dfft%mype3), MPI_DOUBLE_PRECISION, &
                    f_out, recvcount, displs, MPI_DOUBLE_PRECISION, dfft%root,      &
                    dfft%comm3, info )
  !write (6,*) 'gather grid ok 4'
  CALL fftx_error__( ' gather_real_grid', 'info<>0', info )
  !
  ! ... the following check should be performed only on processor dfft%root
  ! ... otherwise f_out must be allocated on all processors even if not used
  !
  info = size( f_out ) - displs( dfft%nproc3-1 ) - recvcount( dfft%nproc3-1 )
  IF( info < 0 ) &
     CALL fftx_error__( ' gather_real_grid ', ' f_out too small ', -info )
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'rgather_grid' )
  !
#else
  CALL fftx_error__(' gather_real_grid', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE gather_real_grid

!----------------------------------------------------------------------------
SUBROUTINE gather_complex_grid ( dfft, f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers a distributed real-space FFT grid to dfft%root, that is,
  ! ... the first processor of input descriptor dfft - complex arrays
  !
  ! ... COMPLEX*16  f_in  = distributed variable (dfft%nnr)
  ! ... COMPLEX*16  f_out = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(in) :: f_in (:)
  COMPLEX(DP), INTENT(inout):: f_out(:)
  TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:)
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info, offset_in, offset_aux, ir3
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfft%nproc-1), recvcount(0:dfft%nproc-1)
  !

  CALL start_clock( 'cgather_grid' )
  !write (*,*) 'gcgather_grid size(f_in),dfft%nnr',size(f_in), dfft%nnr ; FLUSH(6)
  IF( 2*size( f_in ) < dfft%nnr ) &
     CALL fftx_error__( ' gather_complex_grid ', ' f_in too small ', dfft%nnr-size( f_in ) )
  !
  !write (6,*) 'gather grid ok 0 ', dfft%nproc, dfft%nproc2, dfft%nproc3
  ALLOCATE ( f_aux(dfft%nr1x * dfft%nr2x * dfft%my_nr3p ) )
  ! 1) gather within the comm2 communicator
  displs = 0
  DO proc =0, ( dfft%nproc2 -1 )
     recvcount(proc) = 2 * dfft%nr1x * dfft%nr2p(proc+1)
     if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
  end do
  offset_in = 1; offset_aux = 1
  do ir3 = 1, dfft%my_nr3p
     info = 0
     !write (6,*) 'gather grid ok 1 ir3=', ir3
     CALL MPI_GATHERV( f_in(offset_in) , recvcount(dfft%mype2), MPI_DOUBLE_PRECISION, &
                       f_aux(offset_aux), recvcount, displs, MPI_DOUBLE_PRECISION, dfft%root,      &
                       dfft%comm2, info )
     CALL fftx_error__( 'gather_complex_grid', 'info<>0', info )
     !write (6,*) 'gather grid ok 2 ir3=', ir3
     offset_in  = offset_in + dfft%nr1x * dfft%my_nr2p
     offset_aux = offset_aux + dfft%nr1x * dfft%nr2
  end do

  ! 2) gather within the comm3 communicator
  displs = 0
  DO proc = 0, ( dfft%nproc3 - 1 )
     recvcount(proc) = 2 * dfft%nr1x * dfft%nr2x * dfft%nr3p(proc+1)
     if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
  ENDDO
  !
  ! ... the following check should be performed only on processor dfft%root
  ! ... otherwise f_out must be allocated on all processors even if not used
  !
  !write (*,*) 'gcgather_grid 2*size(f_out)',2*size(f_out) ; FLUSH(6)
  !write (*,*) 'gcgather_grid displ+recv',dfft%nproc3, displs(dfft%nproc3-1) + recvcount(dfft%nproc3-1); FLUSH(6)
  info = 2*size( f_out ) - displs( dfft%nproc3 - 1 ) - recvcount( dfft%nproc3-1 ) ; FLUSH(6)
  IF( info < 0 ) CALL fftx_error__( ' gather_complex_grid ', ' f_out too small ', -info )

  info = 0
  !write (6,*) 'gather grid ok 3'
  CALL MPI_GATHERV( f_aux, recvcount(dfft%mype3), MPI_DOUBLE_PRECISION, &
                    f_out, recvcount, displs, MPI_DOUBLE_PRECISION, dfft%root,      &
                    dfft%comm3, info )
  !write (6,*) 'gather grid ok 4'
  CALL fftx_error__( 'gather_complex_grid', 'info<>0', info )
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'cgather_grid' )
  !
#else
  CALL fftx_error__('gather_complex_grid', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE gather_complex_grid

!----------------------------------------------------------------------------
SUBROUTINE scatter_real_grid ( dfft, f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters a real-space FFT grid from dfft%root, first processor of
  ! ... input descriptor dfft, to all others - opposite of "gather_grid"
  !
  ! ... REAL*8  f_in  = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
  ! ... REAL*8  f_out = distributed variable (dfft%nnr)
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(in) :: f_in (:)
  REAL(DP), INTENT(inout):: f_out(:)
  TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
  REAL(DP), ALLOCATABLE ::  f_aux(:)
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info, offset_in, offset_aux, ir3
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfft%nproc-1), sendcount(0:dfft%nproc-1)
  !
  !
  CALL start_clock( 'rscatter_grid' )
  !
  !write (6,*) 'scatter grid ok 0'
  ALLOCATE ( f_aux(dfft%nr1x * dfft%nr2x * dfft%my_nr3p) )
  ! 1) scatter within the comm3 communicator
  displs = 0
  DO proc = 0, ( dfft%nproc3 - 1 )
     sendcount(proc) = dfft%nr1x * dfft%nr2x * dfft%nr3p(proc+1)
     if (proc > 0) displs(proc) = displs(proc-1) + sendcount(proc-1)
  ENDDO
  info = size( f_in ) - displs( dfft%nproc3 - 1 ) - sendcount( dfft%nproc3 - 1 )
  IF( info < 0 ) CALL fftx_error__( ' scatter_real_grid ', ' f_in too small ', -info )
  info = 0
  !write (6,*) 'scatter grid ok 1'
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_aux, sendcount(dfft%mype3), MPI_DOUBLE_PRECISION, &
                     dfft%root, dfft%comm3, info )
  !write (6,*) 'scatter grid ok 2'
  CALL fftx_error__( 'scatter_real_grid', 'info<>0', info )

  ! 2) scatter within the comm2 communicator
  IF( size( f_out ) < dfft%nnr ) &
     CALL fftx_error__( ' scatter_real_grid ', ' f_out too small ', dfft%nnr-size( f_out ) )
  !
  displs = 0 ; f_out = 0.0D0
  DO proc =0, ( dfft%nproc2 -1 )
     sendcount(proc) = dfft%nr1x * dfft%nr2p(proc+1)
     if (proc > 0) displs(proc) = displs(proc-1) + sendcount(proc-1)
  end do
  offset_in = 1 ; offset_aux = 1
  do ir3 = 1, dfft%my_nr3p
     info = 0
  !write (6,*) 'scatter grid ok 3, ir3=', ir3
     CALL MPI_SCATTERV( f_aux(offset_aux), sendcount, displs, MPI_DOUBLE_PRECISION,   &
                        f_out(offset_in), sendcount(dfft%mype2), MPI_DOUBLE_PRECISION, &
                        dfft%root, dfft%comm2, info )
  !write (6,*) 'scatter grid ok 4, ir3=', ir3
     CALL fftx_error__( 'scatter_real_grid', 'info<>0', info )
     offset_in  = offset_in + dfft%nr1x * dfft%my_nr2p
     offset_aux = offset_aux + dfft%nr1x * dfft%nr2
  end do
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'rscatter_grid' )
  !
#else
  CALL fftx_error__('scatter_real_grid', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE scatter_real_grid
!----------------------------------------------------------------------------
SUBROUTINE scatter_complex_grid ( dfft, f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters a real-space FFT grid from dfft%root, first processor of
  ! ... input descriptor dfft, to all others - opposite of "gather_grid"
  !
  ! ... COMPLEX*16  f_in  = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
  ! ... COMPLEX*16  f_out = distributed variable (dfft%nnr)
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(in) :: f_in (:)
  COMPLEX(DP), INTENT(inout):: f_out(:)
  TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:)
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info, offset_in, offset_aux, ir3
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfft%nproc-1), sendcount(0:dfft%nproc-1)
  !
  CALL start_clock( 'cscatter_grid' )
  !
  !write (6,*) 'scatter grid ok 0'
  ALLOCATE ( f_aux(dfft%nr1x * dfft%nr2x * dfft%my_nr3p ) )
  ! 1) scatter within the comm3 communicator
  displs = 0
  DO proc = 0, ( dfft%nproc3 - 1 )
     sendcount(proc) = 2 * dfft%nr1x * dfft%nr2x * dfft%nr3p(proc+1)
     if (proc > 0) displs(proc) = displs(proc-1) + sendcount(proc-1)
  ENDDO
  !
  !write(*,*) 'cscatter_grid 2*size(f_in) ', 2*size(f_in); FLUSH(6)
  !write(*,*) 'cscatter_grid displ+send ', dfft%nproc3, displs(dfft%nproc3-1) + sendcount(dfft%nproc3-1); FLUSH(6)
  info = 2*size( f_in ) - displs( dfft%nproc3 - 1 ) - sendcount( dfft%nproc3 - 1 )
  IF( info < 0 ) &
     CALL fftx_error__( ' scatter_complex_grid ', ' f_in too small ', -info )
  !
  info = 0
  !write (6,*) 'scatter grid ok 1'
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_aux, sendcount(dfft%mype3), MPI_DOUBLE_PRECISION, &
                     dfft%root, dfft%comm3, info )
  !write (6,*) 'scatter grid ok 2'
  CALL fftx_error__( ' scatter_complex_grid', 'info<>0', info )

  ! 2) scatter within the comm2 communicator
  !write(*,*) 'cscatter_grid size(f_out), dfft%nnr ', size(f_out),dfft%nnr; FLUSH(6)
  IF( size( f_out ) < dfft%nnr ) &
     CALL fftx_error__( ' scatter_complex_grid ', ' f_out too small ', dfft%nnr-size( f_out ) )
  !
  displs = 0 ; f_out = 0.0D0
  DO proc =0, ( dfft%nproc2 -1 )
     sendcount(proc) = 2 * dfft%nr1x * dfft%nr2p(proc+1)
     if (proc > 0) displs(proc) = displs(proc-1) + sendcount(proc-1)
  end do
  offset_in = 1 ; offset_aux = 1
  do ir3 = 1, dfft%my_nr3p
     info = 0
  !write (6,*) 'scatter grid ok 3, ir3=', ir3
     CALL MPI_SCATTERV( f_aux(offset_aux), sendcount, displs, MPI_DOUBLE_PRECISION,   &
                        f_out(offset_in), sendcount(dfft%mype2), MPI_DOUBLE_PRECISION, &
                        dfft%root, dfft%comm2, info )
  !write (6,*) 'scatter grid ok 4, ir3=', ir3
     CALL fftx_error__( 'scatter_complex_grid', 'info<>0', info )
     offset_in  = offset_in + dfft%nr1x * dfft%my_nr2p
     offset_aux = offset_aux + dfft%nr1x * dfft%nr2x
  end do
  !
  ! ... the following check should be performed only on processor dfft%root
  ! ... otherwise f_in must be allocated on all processors even if not used
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'cscatter_grid' )
  !
#else
  CALL fftx_error__('scatter_complex_grid', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE scatter_complex_grid
!
! ... "gather"-like subroutines
!
!-----------------------------------------------------------------------
SUBROUTINE cgather_sym( dfftp, f_in, f_out )
  !-----------------------------------------------------------------------
  !
  ! ... gather complex data for symmetrization (used in phonon code)
  ! ... Differs from gather_grid because mpi_allgatherv is used instead
  ! ... of mpi_gatherv - all data is gathered on ALL processors
  ! ... COMPLEX*16  f_in  = distributed variable (nrxx)
  ! ... COMPLEX*16  f_out = gathered variable (nr1x*nr2x*nr3x)
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfftp
  COMPLEX(DP) :: f_in( : ), f_out(:)
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:)
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info, offset_in, offset_aux, ir3
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfftp%nproc-1), recvcount(0:dfftp%nproc-1)
  !
  CALL start_clock( 'cgather' )

  ALLOCATE ( f_aux(dfftp%nr1x * dfftp%nr2x * dfftp%my_nr3p ) )
  !
  CALL MPI_BARRIER( dfftp%comm, info )
  !
  ! 1) gather within the comm2 communicator
  displs = 0
  DO proc =0, ( dfftp%nproc2 -1 )
     recvcount(proc) = 2 * dfftp%nr1x * dfftp%nr2p(proc+1)
     if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
  end do
  offset_in = 1; offset_aux = 1
  do ir3 = 1, dfftp%my_nr3p
     info = 0
     CALL MPI_ALLGATHERV( f_in(offset_in), recvcount(dfftp%mype2), MPI_DOUBLE_PRECISION, &
                          f_aux(offset_aux),recvcount, displs, MPI_DOUBLE_PRECISION, &
                          dfftp%comm2, info )
     CALL fftx_error__( 'cgather_sym', 'info<>0', info )
     offset_in  = offset_in  + dfftp%nr1x * dfftp%my_nr2p
     offset_aux = offset_aux + dfftp%nr1x * dfftp%nr2x
  end do

  ! 2) gather within the comm3 communicator
  displs = 0
  DO proc = 0, ( dfftp%nproc3 - 1 )
     recvcount(proc) = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3p(proc+1)
     if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
  ENDDO
  info = 0
  CALL MPI_ALLGATHERV( f_aux, recvcount(dfftp%mype3), MPI_DOUBLE_PRECISION, &
                       f_out, recvcount, displs, MPI_DOUBLE_PRECISION, &
                       dfftp%comm3, info )
  CALL fftx_error__( 'cgather_sym', 'info<>0', info )
  !
  ! ... the following check should be performed only on processor dfft%root
  ! ... otherwise f_out must be allocated on all processors even if not used
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'cgather' )
  !
#else
  CALL fftx_error__('cgather_sym', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE cgather_sym
!
!
!-----------------------------------------------------------------------
SUBROUTINE cgather_sym_many( dfftp, f_in, f_out, nbnd, nbnd_proc, start_nbnd_proc )
  !-----------------------------------------------------------------------
  !
  ! ... Written by A. Dal Corso
  !
  ! ... This routine generalizes cgather_sym, receiveng nbnd complex
  ! ... distributed functions and collecting nbnd_proc(dfftp%mype+1)
  ! ... functions in each processor.
  ! ... start_nbnd_proc(dfftp%mype+1), says where the data for each processor
  ! ... start in the distributed variable
  ! ... COMPLEX*16  f_in  = distributed variable (nrxx,nbnd)
  ! ... COMPLEX*16  f_out = gathered variable (nr1x*nr2x*nr3x,nbnd_proc(dfftp%mype+1))
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfftp
  INTEGER :: nbnd, nbnd_proc(dfftp%nproc), start_nbnd_proc(dfftp%nproc)
  COMPLEX(DP) :: f_in(dfftp%nnr,nbnd)
  COMPLEX(DP) :: f_out(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,nbnd_proc(dfftp%mype+1))
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:)
  !
#if defined(__MPI)
  !
  INTEGER :: proc_, proc2_, proc3_, proc, info, offset_in, offset_aux, ir3, nr2px
  INTEGER :: ibnd, jbnd, iii
  INTEGER :: displs(0:dfftp%nproc-1), recvcount(0:dfftp%nproc-1)
  !
  !
  CALL start_clock( 'cgather' )
  !
  ALLOCATE ( f_aux(dfftp%nr1x * dfftp%nr2x * dfftp%my_nr3p ) )
  !
  CALL MPI_BARRIER( dfftp%comm, info )
  !
  f_out = (0.d0,0.d0)

  !write (*,*) 'enter cgather ', nbnd
  !write (*,*) nbnd_proc
  !write (*,*) start_nbnd_proc
  !write (*,*) 'dfftp%nproc',dfftp%nproc
  !do ibnd =1,3
  !write (*,*) 'evc = ',ibnd
  !write (*,*) f_in(1:3,ibnd)
  !end do

  nr2px = MAXVAL ( dfftp%nr2p )  ! maximum number of Y values to be disributed

  DO proc_ = 0, dfftp%nproc - 1
     ! define the processor index of the two sub-communicators
     proc2_ = dfftp%iproc2(proc_ + 1) -1 ; proc3_ = dfftp%iproc3(proc_ + 1) -1
     !write (*,*) ' proc_, proc2_, proc3_, dfftp%mype2', proc_, proc2_, proc3_, dfftp%mype2
     DO ibnd = 1, nbnd_proc(proc_ +1)
        jbnd = start_nbnd_proc(proc_ +1) + ibnd - 1
        !write (*,*) ' proc_, ibnd, jbnd ', proc_, ibnd, jbnd
        f_aux(:) = (0.d0,0.d0)
        ! 1) gather within the comm2 communicator
        displs = 0
        DO proc =0, ( dfftp%nproc2 -1 )
           recvcount(proc) = 2 * dfftp%nr1x * dfftp%nr2p(proc+1)
           if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
        end do
        !write (*,*) 'dfftp%nr1x ', dfftp%nr1x
        !write (*,*) 'dfftp%nr2x ', dfftp%nr2x
        !write (*,*) 'dfftp%nr3x ', dfftp%nr3x
        !write (*,*) 'dfftp%nr2p ', dfftp%nr2p, ' nnr3px ', nr2px
        !write (*,*) 'recvcount ', recvcount(0:dfftp%nproc2 -1)
        !write (*,*) 'displs ', displs(0:dfftp%nproc2 -1)
        offset_in = 1; offset_aux = 1
        do ir3 = 1, dfftp%my_nr3p
           info = 0
           CALL MPI_GATHERV( f_in(offset_in,jbnd), recvcount(dfftp%mype2), MPI_DOUBLE_PRECISION, &
                             f_aux(offset_aux), recvcount, displs, MPI_DOUBLE_PRECISION, &
                             proc2_, dfftp%comm2, info )
           CALL fftx_error__( 'cgather_sym_many', 'info<>0', info )
           offset_in  = offset_in  + dfftp%nr1x * dfftp%my_nr2p
           offset_aux = offset_aux + dfftp%nr1x * dfftp%nr2x
        end do
        !write(*,*) ' -> f_in(...+1:...+3,jbnd) ',jbnd
        !offset_in = 0
        !do iii =1,2*dfftp%nr2x
        !   write(*,'(i4,3("(",2f10.7,") "))') iii, f_in(offset_in+1:offset_in+3, jbnd)
        !   offset_in  = offset_in  + dfftp%nr1x
        !enddo
        IF (dfftp%mype2==proc2_) THEN
        !   write(*,*) ' -> f_aux(...+1:...+3)'
        !   offset_aux = 0
        !   do iii =1,4*dfftp%nr2x
        !      write(*,'(i4,3("(",2f10.7,") "))') iii, f_aux(offset_aux+1:offset_aux+3)
        !      offset_aux = offset_aux + dfftp%nr1x
        !   enddo
        !   write(*,*) ' -> F_AUX(...+1:...+3)'
        !   offset_aux = 0
        !   do iii =1,dfftp%my_nr3p
        !      write(*,'(i4,3("(",2f10.7,") "))') iii, f_aux(offset_aux+1:offset_aux+3)
        !      offset_aux  = offset_aux  + dfftp%nr1x*dfftp%nr2x
        !   enddo
           ! 2) gather within the comm3 communicator
           displs = 0
           DO proc = 0, ( dfftp%nproc3 - 1 )
              recvcount(proc) = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3p(proc+1)
              if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
           ENDDO
           info = 0
           CALL MPI_GATHERV( f_aux, recvcount(dfftp%mype3), MPI_DOUBLE_PRECISION, &
                             f_out(1,ibnd), recvcount, displs, MPI_DOUBLE_PRECISION, &
                             proc3_, dfftp%comm3, info )
           CALL fftx_error__( 'cgather_sym_many', 'info<>0', info )

        !   IF (dfftp%mype3==proc3_) THEN
        !      write(*,*) ' -> f_out(...+1:...+3,ibnd) ',ibnd
        !      offset_in  = 0
        !      do iii =1,dfftp%nr3x
        !         write(*,'(i4,3("(",2f10.7,") "))') iii, f_out(offset_in+1:offset_in+3, ibnd)
        !         offset_in  = offset_in  + dfftp%nr1x*dfftp%nr2x
        !      enddo
        !   END IF
        END IF

     END DO
  END DO
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'cgather' )
  !
#else
  CALL fftx_error__('cgather_sym_many', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE cgather_sym_many
!
!----------------------------------------------------------------------------
SUBROUTINE cscatter_sym_many( dfftp, f_in, f_out, target_ibnd, nbnd, nbnd_proc, &
                              start_nbnd_proc   )
  !----------------------------------------------------------------------------
  !
  ! ... Written by A. Dal Corso
  !
  ! ... generalizes cscatter_sym. It assumes that each processor has
  ! ... a certain number of bands (nbnd_proc(dfftp%mype+1)). The processor
  ! ... that has target_ibnd scatters it to all the other processors
  ! ... that receive a distributed part of the target function.
  ! ... start_nbnd_proc(dfftp%mype+1) is used to identify the processor
  ! ... that has the required band
  !
  ! ... COMPLEX*16  f_in  = gathered variable (nr1x*nr2x*nr3x, nbnd_proc(dfftp%mype+1) )
  ! ... COMPLEX*16  f_out = distributed variable (nrxx)
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfftp
  INTEGER :: nbnd, nbnd_proc(dfftp%nproc), start_nbnd_proc(dfftp%nproc)
  COMPLEX(DP) :: f_in(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,nbnd_proc(dfftp%mype+1))
  COMPLEX(DP) :: f_out(dfftp%nnr)
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:)
  INTEGER :: target_ibnd
  !
#if defined(__MPI)
  !
  INTEGER :: proc_, proc2_, proc3_, proc, info, offset_out, offset_aux, ir3
  INTEGER :: displs(0:dfftp%nproc-1), sendcount(0:dfftp%nproc-1)
  INTEGER :: ibnd, jbnd
  !
  CALL start_clock( 'cscatter_sym' )
  !
  ALLOCATE ( f_aux(dfftp%nr1x * dfftp%nr2x * dfftp%my_nr3p ) )
  !
  f_out = (0.0_DP, 0.0_DP)
  !
  CALL MPI_BARRIER( dfftp%comm, info )
  !
  DO proc_ = 0, dfftp%nproc - 1
     ! define the processor index of the two sub-communicators
     proc2_ = dfftp%iproc2(proc_ + 1) -1 ; proc3_ = dfftp%iproc3(proc_ + 1) -1
     DO ibnd = 1, nbnd_proc(proc_+1)
        jbnd = start_nbnd_proc(proc_+1) + ibnd - 1
        IF (jbnd/=target_ibnd) CYCLE
        IF (dfftp%mype2==proc2_) THEN
           ! 1) scatter within the comm3 communicator
           displs = 0
           DO proc = 0, ( dfftp%nproc3 - 1 )
              sendcount(proc) = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3p(proc+1)
              if (proc > 0) displs(proc) = displs(proc-1) + sendcount(proc-1)
           ENDDO
           info = 0
           CALL MPI_SCATTERV( f_in(1,ibnd), sendcount, displs, MPI_DOUBLE_PRECISION,   &
                              f_aux, sendcount(dfftp%mype3), MPI_DOUBLE_PRECISION, &
                              proc3_, dfftp%comm3, info )
        ELSE
           f_aux=(0.d0,0.d0)
        END IF
        ! 2) scatter within the comm2 communicator
        displs = 0 ; f_out = 0.0D0
        DO proc =0, ( dfftp%nproc2 -1 )
           sendcount(proc) = 2 * dfftp%nr1x * dfftp%nr2p(proc+1)
           if (proc > 0) displs(proc) = displs(proc-1) + sendcount(proc-1)
        end do
        offset_out = 1 ; offset_aux = 1
        do ir3 = 1, dfftp%my_nr3p
           info = 0
           CALL MPI_SCATTERV( f_aux(offset_aux), sendcount, displs, MPI_DOUBLE_PRECISION,   &
                              f_out(offset_out), sendcount(dfftp%mype2), MPI_DOUBLE_PRECISION, &
                              proc2_, dfftp%comm2, info )
           CALL fftx_error__( 'gather_grid', 'info<>0', info )
           offset_out = offset_out + dfftp%nr1x * dfftp%my_nr2p
           offset_aux = offset_aux + dfftp%nr1x * dfftp%nr2x
        end do
     ENDDO
  ENDDO
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'cscatter_sym' )
  !
#else
  CALL fftx_error__('cscatter_sym_many', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE cscatter_sym_many

!=----------------------------------------------------------------------=!
   END MODULE scatter_mod
!=----------------------------------------------------------------------=!
!
