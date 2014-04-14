!
! Copyright (C) 2006-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!
!----------------------------------------------------------------------
! FFT base Module.
! Written by Carlo Cavazzoni
!----------------------------------------------------------------------
!
!=----------------------------------------------------------------------=!
   MODULE fft_base
!=----------------------------------------------------------------------=!

        USE kinds, ONLY: DP
        USE parallel_include

        USE fft_types, ONLY: fft_dlay_descriptor

        IMPLICIT NONE

        ! ... data structure containing all information
        ! ... about fft data distribution for a given
        ! ... potential grid, and its wave functions sub-grid.

        TYPE ( fft_dlay_descriptor ) :: dfftp ! descriptor for dense grid
             !  Dimensions of the 3D real and reciprocal space FFT grid
             !  relative to the charge density and potential ("dense" grid)
        TYPE ( fft_dlay_descriptor ) :: dffts ! descriptor for smooth grid
             !  Dimensions of the 3D real and reciprocal space
             !  FFT grid relative to the smooth part of the charge density
             !  (may differ from the full charge density grid for USPP )
        TYPE ( fft_dlay_descriptor ) :: dfftb ! descriptor for box grids
             !  Dimensions of the 3D real and reciprocal space
             !  FFT grid relative to the "small box" computation
             !  of the atomic augmentation part of the 
             !  charge density used in USPP (to speed up CPV iterations)

        SAVE

        PRIVATE

        PUBLIC :: fft_scatter, grid_gather, grid_scatter
        PUBLIC :: dfftp, dffts, dfftb, fft_dlay_descriptor
        PUBLIC :: cgather_sym, cgather_smooth, cgather_custom
        PUBLIC :: cscatter_sym, cscatter_smooth, cscatter_custom
        PUBLIC :: gather_smooth, scatter_smooth
        PUBLIC :: tg_gather, tg_cgather



!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!
!
!
!
!   ALLTOALL based SCATTER, should be better on network
!   with a defined topology, like on bluegene and cray machine
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter ( dfft, f_in, nr3x, nxx_, f_aux, ncp_, npp_, isgn, use_tg )
  !-----------------------------------------------------------------------
  !
  ! transpose the fft grid across nodes
  ! a) From columns to planes (isgn > 0)
  !
  !    "columns" (or "pencil") representation:
  !    processor "me" has ncp_(me) contiguous columns along z
  !    Each column has nr3x elements for a fft of order nr3
  !    nr3x can be =nr3+1 in order to reduce memory conflicts.
  !
  !    The transpose take places in two steps:
  !    1) on each processor the columns are divided into slices along z
  !       that are stored contiguously. On processor "me", slices for
  !       processor "proc" are npp_(proc)*ncp_(me) big
  !    2) all processors communicate to exchange slices
  !       (all columns with z in the slice belonging to "me"
  !        must be received, all the others must be sent to "proc")
  !    Finally one gets the "planes" representation:
  !    processor "me" has npp_(me) complete xy planes
  !
  !  b) From planes to columns (isgn < 0)
  !
  !  Quite the same in the opposite direction
  !
  !  The output is overwritten on f_in ; f_aux is used as work space
  !
  !  If optional argument "use_tg" is true the subroutines performs
  !  the trasposition using the Task Groups distribution
  !
#ifdef __MPI
  USE parallel_include
#endif
  USE kinds,       ONLY : DP

  IMPLICIT NONE

  TYPE (fft_dlay_descriptor), INTENT(in) :: dfft
  INTEGER, INTENT(in)           :: nr3x, nxx_, isgn, ncp_ (:), npp_ (:)
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)
  LOGICAL, OPTIONAL, INTENT(in) :: use_tg

#ifdef __MPI

  INTEGER :: dest, from, k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
  INTEGER :: me_p, nppx, mc, j, npp, nnp, ii, it, ip, ioff, sendsiz, ncpx
  !
  LOGICAL :: use_tg_

#if defined __HPM
     !       CALL f_hpmstart( 10, 'scatter' )
#endif

  !
  !  Task Groups

  use_tg_ = .false.

  IF( present( use_tg ) ) use_tg_ = use_tg

  me     = dfft%mype + 1
  !
  IF( use_tg_ ) THEN
    !  This is the number of procs. in the plane-wave group
     nprocp = dfft%npgrp
  ELSE
     nprocp = dfft%nproc
  ENDIF
  !
  CALL start_clock ('fft_scatter')
  !
  ncpx = 0
  nppx = 0
  IF( use_tg_ ) THEN
     DO proc = 1, nprocp
        gproc = dfft%nplist( proc ) + 1
        ncpx = max( ncpx, ncp_ ( gproc ) )
        nppx = max( nppx, npp_ ( gproc ) )
     ENDDO
  ELSE
     DO proc = 1, nprocp
        ncpx = max( ncpx, ncp_ ( proc ) )
        nppx = max( nppx, npp_ ( proc ) )
     ENDDO
     IF ( dfft%nproc == 1 ) THEN
        nppx = dfft%nr3x
     END IF
  ENDIF
  sendsiz = ncpx * nppx
  !

  ierr = 0
  IF (isgn.gt.0) THEN

     IF (nprocp==1) GO TO 10
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
     offset = 1

     DO proc = 1, nprocp
        from = offset
        IF( use_tg_ ) THEN
           gproc = dfft%nplist(proc)+1
        ELSE
           gproc = proc
        ENDIF
        dest = 1 + ( proc - 1 ) * sendsiz
        !
        DO k = 1, ncp_ (me)
           kdest = dest + (k - 1) * nppx - 1
           kfrom = from + (k - 1) * nr3x - 1
           DO i = 1, npp_ ( gproc )
              f_aux ( kdest + i ) =  f_in ( kfrom + i )
           ENDDO
        ENDDO
        offset = offset + npp_ ( gproc )
     ENDDO

     !
     ! maybe useless; ensures that no garbage is present in the output
     !
     f_in = 0.0_DP
     !
     ! step two: communication
     !
     IF( use_tg_ ) THEN
        gcomm = dfft%pgrp_comm
     ELSE
        gcomm = dfft%comm
     ENDIF

     ! CALL mpi_barrier (gcomm, ierr)  ! why barrier? for buggy openmpi over ib

     CALL mpi_alltoall (f_aux(1), sendsiz, MPI_DOUBLE_COMPLEX, f_in(1), sendsiz, MPI_DOUBLE_COMPLEX, gcomm, ierr)

     IF( abs(ierr) /= 0 ) CALL errore ('fft_scatter', 'info<>0', abs(ierr) )
     !
10   CONTINUE
     !
     f_aux = (0.d0, 0.d0)
     !
     IF( isgn == 1 ) THEN

!!$omp parallel default(none) private(ip,ioff,i,mc,it,j) shared(dfft,nppx,sendsiz,me,f_in,f_aux)
!!$omp do

        DO ip = 1, dfft%nproc
           ioff = dfft%iss( ip )
           DO i = 1, dfft%nsp( ip )
              mc = dfft%ismap( i + ioff )
              it = ( i - 1 ) * nppx + ( ip - 1 ) * sendsiz
              DO j = 1, dfft%npp( me )
                 f_aux( mc + ( j - 1 ) * dfft%nnp ) = f_in( j + it )
              ENDDO
           ENDDO
        ENDDO

!!$omp end do
!!$omp end parallel

     ELSE

        IF( use_tg_ ) THEN
           npp  = dfft%tg_npp( me )
           nnp  = dfft%nr1x * dfft%nr2x
        ELSE
           npp  = dfft%npp( me )
           nnp  = dfft%nnp
        ENDIF
        !
!!$omp parallel default(none) private(ip,ioff,i,mc,it,j,gproc,ii) shared(dfft,nppx,npp,nnp,sendsiz,use_tg_,f_in,f_aux)
!!$omp do
        DO ip = 1, dfft%nproc

           IF( use_tg_ ) THEN
              gproc = ( ip - 1 ) / dfft%nogrp + 1
              IF( MOD( ip - 1, dfft%nogrp ) == 0 ) ii = 0
           ELSE
              gproc = ip
              ii = 0
           ENDIF
           !
           ioff = dfft%iss( ip )
           !
           DO i = 1, dfft%nsw( ip )
              !
              mc = dfft%ismap( i + ioff )
              !
              it = ii * nppx + ( gproc - 1 ) * sendsiz
              !
              DO j = 1, npp
                 f_aux( mc + ( j - 1 ) * nnp ) = f_in( j + it )
              ENDDO
              !
              ii = ii + 1
              !
           ENDDO
           !
        ENDDO
!!$omp end do
!!$omp end parallel

     END IF

  ELSE
     !
     !  "backward" scatter from planes to columns
     !
     IF( isgn == -1 ) THEN

        DO ip = 1, dfft%nproc
           ioff = dfft%iss( ip )
           DO i = 1, dfft%nsp( ip )
              mc = dfft%ismap( i + ioff )
              it = ( i - 1 ) * nppx + ( ip - 1 ) * sendsiz
              DO j = 1, dfft%npp( me )
                 f_in( j + it ) = f_aux( mc + ( j - 1 ) * dfft%nnp )
              ENDDO
           ENDDO
        ENDDO

     ELSE

        IF( use_tg_ ) THEN
           npp  = dfft%tg_npp( me )
           nnp  = dfft%nr1x * dfft%nr2x
        ELSE
           npp  = dfft%npp( me )
           nnp  = dfft%nnp
        ENDIF

        DO ip = 1, dfft%nproc

           IF( use_tg_ ) THEN
              gproc = ( ip - 1 ) / dfft%nogrp + 1
              IF( MOD( ip - 1, dfft%nogrp ) == 0 ) ii = 0
           ELSE
              gproc = ip
              ii = 0
           ENDIF
           !
           ioff = dfft%iss( ip )
           !
           DO i = 1, dfft%nsw( ip )
              !
              mc = dfft%ismap( i + ioff )
              !
              it = ii * nppx + ( gproc - 1 ) * sendsiz
              !
              DO j = 1, npp
                 f_in( j + it ) = f_aux( mc + ( j - 1 ) * nnp )
              ENDDO
              !
              ii = ii + 1
              !
           ENDDO

        ENDDO

     END IF

     IF( nprocp == 1 ) GO TO 20
     !
     !  step two: communication
     !
     IF( use_tg_ ) THEN
        gcomm = dfft%pgrp_comm
     ELSE
        gcomm = dfft%comm
     ENDIF

     ! CALL mpi_barrier (gcomm, ierr)  ! why barrier? for buggy openmpi over ib

     CALL mpi_alltoall (f_in(1), sendsiz, MPI_DOUBLE_COMPLEX, f_aux(1), sendsiz, MPI_DOUBLE_COMPLEX, gcomm, ierr)

     IF( abs(ierr) /= 0 ) CALL errore ('fft_scatter', 'info<>0', abs(ierr) )
     !
     !  step one: store contiguously the columns
     !
     f_in = 0.0_DP
     !
     offset = 1

     DO proc = 1, nprocp
        from = offset
        IF( use_tg_ ) THEN
           gproc = dfft%nplist(proc)+1
        ELSE
           gproc = proc
        ENDIF
        dest = 1 + ( proc - 1 ) * sendsiz
        !
        DO k = 1, ncp_ (me)
           kdest = dest + (k - 1) * nppx - 1
           kfrom = from + (k - 1) * nr3x - 1
           DO i = 1, npp_ ( gproc )  
              f_in ( kfrom + i ) = f_aux ( kdest + i )
           ENDDO
        ENDDO
        offset = offset + npp_ ( gproc )
     ENDDO

20   CONTINUE

  ENDIF

  CALL stop_clock ('fft_scatter')

#endif

#if defined __HPM
     !       CALL f_hpmstop( 10 )
#endif

  RETURN

END SUBROUTINE fft_scatter


!----------------------------------------------------------------------------
SUBROUTINE grid_gather( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers nproc distributed data on the first processor of every pool
  !
  ! ... REAL*8  f_in  = distributed variable (nxx)
  ! ... REAL*8  f_out = gathered variable (nr1x*nr2x*nr3x)
  !
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  REAL(DP) :: f_in( : ), f_out( : )
  !
#if defined (__MPI)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:dfftp%nproc-1), recvcount(0:dfftp%nproc-1)
  !
  IF( size( f_in ) < dfftp%nnr ) &
     CALL errore( ' grid_gather ', ' f_in too small ', dfftp%nnr - size( f_in ) )
  !
  CALL start_clock( 'gather' )
  !
  DO proc = 0, ( dfftp%nproc - 1 )
     !
     recvcount(proc) = dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + recvcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  info = size( f_out ) - displs( dfftp%nproc - 1 ) - recvcount( dfftp%nproc - 1 )
  !
  IF( info < 0 ) &
     CALL errore( ' grid_gather ', ' f_out too small ', -info )
  !
  info = 0
  !
  CALL MPI_GATHERV( f_in, recvcount(dfftp%mype), MPI_DOUBLE_PRECISION, f_out, &
                    recvcount, displs, MPI_DOUBLE_PRECISION, dfftp%root,    &
                    dfftp%comm, info )
  !
  CALL errore( 'grid_gather', 'info<>0', info )
  !
  CALL stop_clock( 'gather' )
  !
#else
  CALL errore('grid_gather', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE grid_gather


!----------------------------------------------------------------------------
SUBROUTINE grid_scatter( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters data from the first processor of every pool
  !
  ! ... REAL*8  f_in  = gathered variable (nr1x*nr2x*nr3x)
  ! ... REAL*8  f_out = distributed variable (nxx)
  !
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  REAL(DP) :: f_in( : ), f_out( : )
  !
#if defined (__MPI)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:dfftp%nproc-1), sendcount(0:dfftp%nproc-1)
  !
  IF( size( f_out ) < dfftp%nnr ) &
     CALL errore( ' grid_scatter ', ' f_out too small ', dfftp%nnr - size( f_in ) )
  !
  CALL start_clock( 'scatter' )
  !
  DO proc = 0, ( dfftp%nproc - 1 )
     !
     sendcount(proc) = dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + sendcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  info = size( f_in ) - displs( dfftp%nproc - 1 ) - sendcount( dfftp%nproc - 1 )
  !
  IF( info < 0 ) &
     CALL errore( ' grid_scatter ', ' f_in too small ', -info )
  !
  info = 0
  !
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_out, sendcount(dfftp%mype), MPI_DOUBLE_PRECISION, &
                     dfftp%root, dfftp%comm, info )
  !
  CALL errore( 'grid_scatter', 'info<>0', info )
  !
  IF ( sendcount(dfftp%mype) /= dfftp%nnr ) f_out(sendcount(dfftp%mype)+1:dfftp%nnr) = 0.D0
  !
  CALL stop_clock( 'scatter' )
  !
#else
  CALL errore('grid_scatter', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE grid_scatter
!
! ... "gather"-like subroutines
!
!-----------------------------------------------------------------------
SUBROUTINE cgather_sym( f_in, f_out )
  !-----------------------------------------------------------------------
  !
  ! ... gather complex data for symmetrization (in phonon code)
  ! ... COMPLEX*16  f_in  = distributed variable (nrxx)
  ! ... COMPLEX*16  f_out = gathered variable (nr1x*nr2x*nr3x)
  !
  USE mp,        ONLY : mp_barrier
  USE parallel_include
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in( : ), f_out(:)
  !
#if defined (__MPI)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:dfftp%nproc-1), recvcount(0:dfftp%nproc-1)
  !
  !
  CALL start_clock( 'cgather' )
  !
  DO proc = 0, ( dfftp%nproc - 1 )
     !
     recvcount(proc) = 2 * dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + recvcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_barrier( dfftp%comm )
  !
  CALL MPI_ALLGATHERV( f_in, recvcount(dfftp%mype), MPI_DOUBLE_PRECISION, &
                       f_out, recvcount, displs, MPI_DOUBLE_PRECISION, &
                       dfftp%comm, info )
  !
  CALL errore( 'cgather_sym', 'info<>0', info )
  !
!  CALL mp_barrier( dfftp%comm )
  !
  CALL stop_clock( 'cgather' )
  !
#else
  CALL errore('cgather_sym', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE cgather_sym
!
!----------------------------------------------------------------------------
SUBROUTINE cgather_smooth ( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers data on the smooth AND complex fft grid
  !
  ! ... gathers nproc distributed data on the first processor of every pool
  !
  ! ... COMPLEX*16  f_in  = distributed variable ( dffts%nnr )
  ! ... COMPLEX*16  f_out = gathered variable (nr1sx*nr2sx*nr3sx)
  !
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in(:), f_out(:)
  !
#if defined (__MPI)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:dfftp%nproc-1), recvcount(0:dfftp%nproc-1)
  !
  !
  CALL start_clock( 'gather' )
  !
  DO proc = 0, ( dfftp%nproc - 1 )
     !
     recvcount(proc) = 2 * dffts%nnp * dffts%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + recvcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_barrier( dfftp%comm )
  !
  CALL MPI_GATHERV( f_in, recvcount(dfftp%mype), MPI_DOUBLE_PRECISION, f_out, &
                    recvcount, displs, MPI_DOUBLE_PRECISION, dfftp%root,    &
                    dfftp%comm, info )
  !
  CALL errore( 'cgather_smooth', 'info<>0', info )
  !
  CALL stop_clock( 'gather' )
  !
#else
  CALL errore('cgather_smooth', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE cgather_smooth
!
!----------------------------------------------------------------------------
SUBROUTINE cgather_custom ( f_in, f_out, dfftt )
  !----------------------------------------------------------------------------
  !
  ! ... gathers data on the custom AND complex fft grid
  !
  ! ... gathers nproc distributed data on the first processor of every pool
  !
  ! ... COMPLEX*16  f_in  = distributed variable ( dfftt%nnr )
  ! ... COMPLEX*16  f_out = gathered variable (nr1sx*nr2sx*nr3sx)
  !
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in(:), f_out(:)
  TYPE ( fft_dlay_descriptor ), INTENT(IN) :: dfftt 
  !
#if defined (__MPI)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:dfftp%nproc-1), recvcount(0:dfftp%nproc-1)
  !
  !
  CALL start_clock( 'gather' )
  !
  DO proc = 0, ( dfftp%nproc - 1 )
     !
     recvcount(proc) = 2 * dfftt%nnp * dfftt%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + recvcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_barrier( dfftp%comm )
  !
  CALL MPI_GATHERV( f_in, recvcount(dfftp%mype), MPI_DOUBLE_PRECISION, f_out, &
                    recvcount, displs, MPI_DOUBLE_PRECISION, dfftp%root,    &
                    dfftp%comm, info )
  !
  CALL errore( 'cgather_custom', 'info<>0', info )
  !
  CALL stop_clock( 'gather' )
  !
#else
  CALL errore('cgather_custom', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE cgather_custom
!
! ... "scatter"-like subroutines
!
!----------------------------------------------------------------------------
SUBROUTINE cscatter_sym( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters data from the first processor of every pool
  !
  ! ... COMPLEX*16  f_in  = gathered variable (nr1x*nr2x*nr3x)
  ! ... COMPLEX*16  f_out = distributed variable (nxx)
  !
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in(:), f_out(:)
  !
#if defined (__MPI)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:dfftp%nproc-1), sendcount(0:dfftp%nproc-1)
  !
  !
  CALL start_clock( 'cscatter_sym' )
  !
  DO proc = 0, ( dfftp%nproc - 1 )
     !
     sendcount(proc) = 2 * dfftp%nnp * dfftp%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + sendcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_barrier( dfftp%comm )
  !
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_out, sendcount(dfftp%mype), MPI_DOUBLE_PRECISION, &
                     dfftp%root, dfftp%comm, info )
  !
  CALL errore( 'cscatter_sym', 'info<>0', info )
  !
  IF ( sendcount(dfftp%mype) /=  dfftp%nnr  ) f_out(sendcount(dfftp%mype)+1: dfftp%nnr ) = 0.D0
  !
  CALL stop_clock( 'cscatter_sym' )
  !
#else
  CALL errore('cscatter_sym', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE cscatter_sym
!
!----------------------------------------------------------------------------
SUBROUTINE cscatter_smooth( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters data on the smooth AND complex fft grid
  ! ... scatters data from the first processor of every pool
  !
  ! ... COMPLEX*16  f_in  = gathered variable (nr1sx*nr2sx*nr3sx)
  ! ... COMPLEX*16  f_out = distributed variable ( dffts%nnr)
  !
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in(:), f_out(:)
  !
#if defined (__MPI)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:dfftp%nproc-1), sendcount(0:dfftp%nproc-1)
  !
  !
  CALL start_clock( 'scatter' )
  !
  DO proc = 0, ( dfftp%nproc - 1 )
     !
     sendcount(proc) = 2 * dffts%nnp * dffts%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + sendcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_barrier( dfftp%comm )
  !
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_out, sendcount(dfftp%mype), MPI_DOUBLE_PRECISION, &
                     dfftp%root, dfftp%comm, info )
  !
  CALL errore( 'scatter', 'info<>0', info )
  !
  IF ( sendcount(dfftp%mype) /=  dffts%nnr  ) f_out(sendcount(dfftp%mype)+1: dffts%nnr ) = 0.D0
  !
  CALL stop_clock( 'scatter' )
  !
#else
  CALL errore('cscatter_smooth', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE cscatter_smooth
!
!----------------------------------------------------------------------------
SUBROUTINE cscatter_custom( f_in, f_out, dfftt )
  !----------------------------------------------------------------------------
  !
  ! ... scatters data on the custom AND complex fft grid
  ! ... scatters data from the first processor of every pool
  !
  ! ... COMPLEX*16  f_in  = gathered variable (nr1sx*nr2sx*nr3sx)
  ! ... COMPLEX*16  f_out = distributed variable ( dfftt%nnr)
  !
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: f_in(:), f_out(:)
  TYPE ( fft_dlay_descriptor ), INTENT(IN) :: dfftt 
  !
#if defined (__MPI)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:dfftp%nproc-1), sendcount(0:dfftp%nproc-1)
  !
  !
  CALL start_clock( 'scatter' )
  !
  DO proc = 0, ( dfftp%nproc - 1 )
     !
     sendcount(proc) = 2 * dfftt%nnp * dfftt%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + sendcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_barrier( dfftp%comm )
  !
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_out, sendcount(dfftp%mype), MPI_DOUBLE_PRECISION, &
                     dfftp%root, dfftp%comm, info )
  !
  CALL errore( 'scatter', 'info<>0', info )
  !
  IF ( sendcount(dfftp%mype) /=  dfftt%nnr  ) f_out(sendcount(dfftp%mype)+1: dfftt%nnr ) = 0.D0
  !
  CALL stop_clock( 'scatter' )
  !
#else
  CALL errore('cscatter_custom', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE cscatter_custom
!
!----------------------------------------------------------------------------
SUBROUTINE gather_smooth ( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers data on the smooth AND real fft grid
  !
  ! ... gathers nproc distributed data on the first processor of every pool
  !
  ! ... REAL*8      f_in  = distributed variable ( dffts%nnr )
  ! ... REAL*8      f_out = gathered variable (nr1sx*nr2sx*nr3sx)
  !
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  REAL(DP) :: f_in(:), f_out(:)
  !
#if defined (__MPI)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:dffts%nproc-1), recvcount(0:dffts%nproc-1)
  !
  !
  CALL start_clock( 'gather' )
  !
  DO proc = 0, ( dffts%nproc - 1 )
     !
     recvcount(proc) = dffts%nnp * dffts%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + recvcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_barrier( dffts%comm )
  !
  CALL MPI_GATHERV( f_in, recvcount(dffts%mype), MPI_DOUBLE_PRECISION, f_out, &
                    recvcount, displs, MPI_DOUBLE_PRECISION, dffts%root,    &
                    dffts%comm, info )
  !
  CALL errore( 'gather', 'info<>0', info )
  !
  CALL stop_clock( 'gather' )
  !
#else
  CALL errore('gather_smooth', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE gather_smooth
!
!----------------------------------------------------------------------------
SUBROUTINE scatter_smooth( f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters data on the smooth AND real fft grid
  ! ... scatters data from the first processor of every pool
  !
  ! ... REAL*8      f_in  = gathered variable (nr1sx*nr2sx*nr3sx)
  ! ... REAL*8      f_out = distributed variable ( dffts%nnr)
  !
  USE mp,        ONLY : mp_barrier
  USE kinds,     ONLY : DP
  USE parallel_include
  !
  IMPLICIT NONE
  !
  REAL(DP) :: f_in(:), f_out(:)
  !
#if defined (__MPI)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:dffts%nproc-1), sendcount(0:dffts%nproc-1)
  !
  !
  CALL start_clock( 'scatter' )
  !
  DO proc = 0, ( dffts%nproc - 1 )
     !
     sendcount(proc) = dffts%nnp * dffts%npp(proc+1)
     !
     IF ( proc == 0 ) THEN
        !
        displs(proc) = 0
        !
     ELSE
        !
        displs(proc) = displs(proc-1) + sendcount(proc-1)
        !
     ENDIF
     !
  ENDDO
  !
  CALL mp_barrier( dffts%comm )
  !
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_out, sendcount(dffts%mype), MPI_DOUBLE_PRECISION, &
                     dffts%root, dffts%comm, info )
  !
  CALL errore( 'scatter', 'info<>0', info )
  !
  IF ( sendcount(dffts%mype) /=  dffts%nnr  ) f_out(sendcount(dffts%mype)+1: dffts%nnr ) = 0.D0
  !
  CALL stop_clock( 'scatter' )
  !
#else
  CALL errore('scatter_smooth', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE scatter_smooth


!
SUBROUTINE tg_gather( dffts, v, tg_v )
   !
   USE parallel_include
   !
   USE fft_types,      ONLY : fft_dlay_descriptor

   ! T.G.
   ! NOGRP:      Number of processors per orbital task group

   IMPLICIT NONE

   TYPE(fft_dlay_descriptor), INTENT(in) :: dffts

   REAL(DP) :: v(:)
   REAL(DP) :: tg_v(:)

   INTEGER :: nsiz, i, ierr, nsiz_tg
   INTEGER :: recv_cnt( dffts%nogrp ), recv_displ( dffts%nogrp )

   nsiz_tg = dffts%tg_nnr * dffts%nogrp

   IF( size( tg_v ) < nsiz_tg ) &
      CALL errore( ' tg_gather ', ' tg_v too small ', ( nsiz_tg - size( tg_v ) ) )

   nsiz = dffts%npp( dffts%mype+1 ) * dffts%nr1x * dffts%nr2x

   IF( size( v ) < nsiz ) &
      CALL errore( ' tg_gather ', ' v too small ',  ( nsiz - size( v ) ) )

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

#if defined (__MPI)

   CALL MPI_Allgatherv( v(1), nsiz, MPI_DOUBLE_PRECISION, &
        tg_v(1), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, dffts%ogrp_comm, IERR)

   IF( ierr /= 0 ) &
      CALL errore( ' tg_gather ', ' MPI_Allgatherv ', abs( ierr ) )

#endif

END SUBROUTINE tg_gather
!
!  Complex version of previous routine
!
SUBROUTINE tg_cgather( dffts, v, tg_v )
   !
   USE parallel_include
   !
   USE fft_types,      ONLY : fft_dlay_descriptor

   ! T.G.
   ! NOGRP:      Number of processors per orbital task group

   IMPLICIT NONE

   TYPE(fft_dlay_descriptor), INTENT(in) :: dffts

   COMPLEX(DP) :: v(:)
   COMPLEX(DP) :: tg_v(:)

   INTEGER :: nsiz, i, ierr, nsiz_tg
   INTEGER :: recv_cnt( dffts%nogrp ), recv_displ( dffts%nogrp )

   nsiz_tg = dffts%tg_nnr * dffts%nogrp

   IF( size( tg_v ) < nsiz_tg ) &
      CALL errore( ' tg_gather ', ' tg_v too small ', ( nsiz_tg - size( tg_v ) ) )

   nsiz = dffts%npp( dffts%mype+1 ) * dffts%nr1x * dffts%nr2x

   IF( size( v ) < nsiz ) &
      CALL errore( ' tg_gather ', ' v too small ',  ( nsiz - size( v ) ) )

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
#if defined (__MPI)

   CALL MPI_Allgatherv( v(1), nsiz, MPI_DOUBLE_PRECISION, &
        tg_v(1), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, dffts%ogrp_comm, IERR)

   IF( ierr /= 0 ) &
      CALL errore( ' tg_cgather ', ' MPI_Allgatherv ', abs( ierr ) )

#endif

END SUBROUTINE tg_cgather

!=----------------------------------------------------------------------=!
   END MODULE fft_base
!=----------------------------------------------------------------------=!
