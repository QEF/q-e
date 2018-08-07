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
!----------------------------------------------------------------------
!
!=----------------------------------------------------------------------=!
   MODULE scatter_mod_2d
!=----------------------------------------------------------------------=!

        USE fft_types, ONLY: fft_type_descriptor
        USE fft_param

        IMPLICIT NONE



#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
#if __INTEL_COMPILER  < 9999
!dir$ attributes align: 4096 :: scra
#endif
#endif
#endif
        COMPLEX(DP), ALLOCATABLE :: scra(:)

        SAVE

        PRIVATE

        PUBLIC :: fft_type_descriptor
        PUBLIC :: fft_scatter

!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!
!
!
#if ! defined __NON_BLOCKING_SCATTER
!
!   ALLTOALL based SCATTER, should be better on network
!   with a defined topology, like on bluegene and cray machine
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter ( dfft, f_in, nr3x, nxx_, f_aux, ncp_, npp_, isgn )
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
  !    f_in  contains input columns, is destroyed on output
  !    f_aux contains output planes
  !
  !  b) From planes to columns (isgn < 0)
  !
  !    Quite the same in the opposite direction
  !    f_aux contains input planes, is destroyed on output
  !    f_in  contains output columns
  !
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in) :: dfft
  INTEGER, INTENT(in)           :: nr3x, nxx_, isgn, ncp_ (:), npp_ (:)
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)

#if defined(__MPI)

  INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
  INTEGER :: me_p, nppx, mc, j, npp, nnp, ii, it, ip, ioff, sendsiz, ncpx, ipp, nblk, nsiz
  !
  !  Task Groups


  me     = dfft%mype + 1
  !
  nprocp = dfft%nproc
  !
  CALL start_clock ('fft_scatter')
  !
  ncpx = 0
  nppx = 0
  DO proc = 1, nprocp
     ncpx = max( ncpx, ncp_ ( proc ) )
     nppx = max( nppx, npp_ ( proc ) )
  ENDDO
  IF ( dfft%nproc == 1 ) THEN
     nppx = dfft%nr3x
  END IF
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
     offset = 0

     DO proc = 1, nprocp
        gproc = proc
        !
        kdest = ( proc - 1 ) * sendsiz
        kfrom = offset
        !
        DO k = 1, ncp_ (me)
           DO i = 1, npp_ ( gproc )
              f_aux ( kdest + i ) =  f_in ( kfrom + i )
           ENDDO
           kdest = kdest + nppx
           kfrom = kfrom + nr3x
        ENDDO
        offset = offset + npp_ ( gproc )
     ENDDO
     !
     ! step two: communication
     !

     gcomm = dfft%comm

     CALL mpi_alltoall (f_aux(1), sendsiz, MPI_DOUBLE_COMPLEX, f_in(1), sendsiz, MPI_DOUBLE_COMPLEX, gcomm, ierr)

     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
     !
10   CONTINUE
     !
     f_aux = (0.d0, 0.d0)
     !
     IF( isgn == 1 ) THEN

        DO ip = 1, dfft%nproc
           ioff = dfft%iss( ip )
           it = ( ip - 1 ) * sendsiz
           DO i = 1, dfft%nsp( ip )
              mc = dfft%ismap( i + ioff )
              DO j = 1, dfft%nr3p( me )
                 f_aux( mc + ( j - 1 ) * dfft%nnp ) = f_in( j + it )
              ENDDO
              it = it + nppx
           ENDDO
        ENDDO

     ELSE

        npp  = dfft%nr3p( me )
        nnp  = dfft%nnp

        nblk = dfft%nproc 
        nsiz = 1
        !
        ip = 1
        !
        DO gproc = 1, nblk
           !
           ii = 0
           !
           DO ipp = 1, nsiz
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
              ip = ip + 1
              !
           ENDDO
           !
        ENDDO

     END IF

  ELSE
     !
     !  "backward" scatter from planes to columns
     !
     IF( isgn == -1 ) THEN

        DO ip = 1, dfft%nproc
           ioff = dfft%iss( ip )
           it = ( ip - 1 ) * sendsiz
           DO i = 1, dfft%nsp( ip )
              mc = dfft%ismap( i + ioff )
              DO j = 1, dfft%nr3p( me )
                 f_in( j + it ) = f_aux( mc + ( j - 1 ) * dfft%nnp )
              ENDDO
              it = it + nppx
           ENDDO
        ENDDO

     ELSE

        npp  = dfft%nr3p( me )
        nnp  = dfft%nnp

        nblk = dfft%nproc 
        nsiz = 1
        !
        ip = 1
        !
        DO gproc = 1, nblk
           !
           ii = 0
           !
           DO ipp = 1, nsiz
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
              !
              ip = ip + 1
              !
           ENDDO
           !
        ENDDO

     END IF

     IF( nprocp == 1 ) GO TO 20
     !
     !  step two: communication
     !
     gcomm = dfft%comm

     CALL mpi_alltoall (f_in(1), sendsiz, MPI_DOUBLE_COMPLEX, f_aux(1), sendsiz, MPI_DOUBLE_COMPLEX, gcomm, ierr)

     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
     !
     !  step one: store contiguously the columns
     !
     !! f_in = 0.0_DP
     !
     offset = 0

     DO proc = 1, nprocp
        gproc = proc
        !
        kdest = ( proc - 1 ) * sendsiz
        kfrom = offset 
        !
        DO k = 1, ncp_ (me)
           DO i = 1, npp_ ( gproc )  
              f_in ( kfrom + i ) = f_aux ( kdest + i )
           ENDDO
           kdest = kdest + nppx
           kfrom = kfrom + nr3x
        ENDDO
        offset = offset + npp_ ( gproc )
     ENDDO

20   CONTINUE

  ENDIF

  CALL stop_clock ('fft_scatter')

#endif

  RETURN

END SUBROUTINE fft_scatter
!
#else
!
!   NON BLOCKING SCATTER, NOT IMPLEMENTED YET
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter ( dfft, f_in, nr3x, nxx_, f_aux, ncp_, npp_, isgn )
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
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in) :: dfft
  INTEGER, INTENT(in)           :: nr3x, nxx_, isgn, ncp_ (:), npp_ (:)
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)

#if defined(__MPI)

  INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
  INTEGER :: me_p, nppx, mc, j, npp, nnp, ii, it, ip, ioff, sendsiz, ncpx, ipp, nblk, nsiz, ijp
  INTEGER :: sh(dfft%nproc), rh(dfft%nproc)
  INTEGER :: istat( MPI_STATUS_SIZE )
  !
  LOGICAL :: use_tg_
  !
  CALL start_clock ('fft_scatter')

  use_tg_ = .false.

  IF( dfft%has_task_groups ) use_tg_ = .true.
  IF( abs(isgn) == 1 ) use_tg_ = .false.

  IF( .NOT. ALLOCATED( scra ) ) THEN
     ALLOCATE( scra( nxx_ ) )
  END IF
  IF( SIZE( scra ) .LT. nxx_ ) THEN
     DEALLOCATE( scra )
     ALLOCATE( scra( nxx_ ) )
  END IF

  me     = dfft%mype + 1
  !
  ncpx = 0
  nppx = 0
  IF( use_tg_ ) THEN
     !  This is the number of procs. in the plane-wave group
     nprocp = dfft%npgrp
     ncpx   = dfft%tg_ncpx
     nppx   = dfft%tg_nppx
     gcomm  = dfft%pgrp_comm
  ELSE
     nprocp = dfft%nproc
     DO proc = 1, nprocp
        ncpx = max( ncpx, ncp_ ( proc ) )
        nppx = max( nppx, npp_ ( proc ) )
     ENDDO
     IF ( dfft%nproc == 1 ) THEN
        nppx = dfft%nr3x
     END IF
     gcomm = dfft%comm
  ENDIF
  ! 
  sendsiz = ncpx * nppx
  !
  IF ( isgn .gt. 0 ) THEN
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
     offset = 0

     IF( use_tg_ ) THEN
        DO proc = 1, nprocp
           gproc = dfft%nplist(proc)+1
           kdest = ( proc - 1 ) * sendsiz
           kfrom = offset 
           DO k = 1, ncp_ (me)
              DO i = 1, npp_ ( gproc )
                 f_aux ( kdest + i ) =  f_in ( kfrom + i )
              ENDDO
              kdest = kdest + nppx
              kfrom = kfrom + nr3x
           ENDDO
           offset = offset + npp_ ( gproc )
           ! post the non-blocking send, f_aux can't be overwritten until operation has completed
           CALL mpi_isend( f_aux( (proc-1)*sendsiz + 1 ), sendsiz, MPI_DOUBLE_COMPLEX, proc-1, me, gcomm, sh( proc ), ierr )
           ! IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', ' forward send info<>0', abs(ierr) )
        ENDDO
     ELSE
        DO proc = 1, nprocp
           kdest = ( proc - 1 ) * sendsiz
           kfrom = offset 
           DO k = 1, ncp_ (me)
              DO i = 1, npp_ ( proc )
                 f_aux ( kdest + i ) =  f_in ( kfrom + i )
              ENDDO
              kdest = kdest + nppx
              kfrom = kfrom + nr3x
           ENDDO
           offset = offset + npp_ ( proc )
           ! post the non-blocking send, f_aux can't be overwritten until operation has completed
           CALL mpi_isend( f_aux( (proc-1)*sendsiz + 1 ), sendsiz, MPI_DOUBLE_COMPLEX, proc-1, me, gcomm, sh( proc ), ierr )
           ! IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', ' forward send info<>0', abs(ierr) )
        ENDDO
     ENDIF
     !
     ! step two: receive
     !
     DO proc = 1, nprocp
        !
        ! now post the receive
        !
        CALL mpi_irecv( f_in( (proc-1)*sendsiz + 1 ), sendsiz, MPI_DOUBLE_COMPLEX, proc-1, MPI_ANY_TAG, gcomm, rh( proc ), ierr )
        !IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', ' forward receive info<>0', abs(ierr) )
        !
        !
     ENDDO
     !
     scra = (0.d0, 0.d0)
     !
     call mpi_waitall( nprocp, rh, MPI_STATUSES_IGNORE, ierr )
     !
     IF( isgn == 1 ) THEN

        DO ip = 1, dfft%nproc
           ioff = dfft%iss( ip )
           it = ( ip - 1 ) * sendsiz
           DO i = 1, dfft%nsp( ip )
              mc = dfft%ismap( i + ioff )
              DO j = 1, dfft%npp( me )
                 scra( mc + ( j - 1 ) * dfft%nnp ) = f_in( j + it )
              ENDDO
              it = it + nppx
           ENDDO
        ENDDO

     ELSE

        IF( use_tg_ ) THEN
           npp  = dfft%tg_npp( me )
           nnp  = dfft%nr1x * dfft%nr2x
           nblk = dfft%nproc / dfft%nogrp
           nsiz = dfft%nogrp
        ELSE
           npp  = dfft%npp( me )
           nnp  = dfft%nnp
           nblk = dfft%nproc 
           nsiz = 1
        ENDIF
        !
        ip = 1
        !
        DO gproc = 1, nblk
           !
           ii = 0
           !
           DO ipp = 1, nsiz
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
                    scra( mc + ( j - 1 ) * nnp ) = f_in( j + it )
                 ENDDO
                 !
                 ii = ii + 1
                 !
              ENDDO
              !
              ip = ip + 1
              !
           ENDDO
           !
        ENDDO
        !
     END IF

     call mpi_waitall( nprocp, sh, MPI_STATUSES_IGNORE, ierr )

     f_aux = scra( 1:SIZE(f_aux) )

  ELSE
     !
     !  "backward" scatter from planes to columns
     !
     IF( isgn == -1 ) THEN

        nblk = dfft%nproc 

        DO ip = 1, dfft%nproc
           ioff = dfft%iss( ip )
           it = ( ip - 1 ) * sendsiz
           DO i = 1, dfft%nsp( ip )
              mc = dfft%ismap( i + ioff )
              DO j = 1, dfft%npp( me )
                 f_in( j + it ) = f_aux( mc + ( j - 1 ) * dfft%nnp )
              ENDDO
              it = it + nppx
           ENDDO

           CALL mpi_isend( f_in( (ip-1)*sendsiz + 1 ), sendsiz, MPI_DOUBLE_COMPLEX, ip-1, me, gcomm, sh( ip ), ierr )
           ! IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', ' backward send info<>0', abs(ierr) )

        ENDDO

        DO ip = 1, dfft%nproc
           CALL mpi_irecv( f_aux( (ip-1)*sendsiz + 1 ), sendsiz, MPI_DOUBLE_COMPLEX, ip-1, MPI_ANY_TAG, gcomm, rh(ip), ierr )
           ! IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', ' backward receive info<>0', abs(ierr) )
        ENDDO

     ELSE

        IF( use_tg_ ) THEN
           npp  = dfft%tg_npp( me )
           nnp  = dfft%nr1x * dfft%nr2x
           nblk = dfft%nproc / dfft%nogrp
           nsiz = dfft%nogrp
        ELSE
           npp  = dfft%npp( me )
           nnp  = dfft%nnp
           nblk = dfft%nproc 
           nsiz = 1
        ENDIF
        !
        DO gproc = 0, nblk-1

           ii = 0
           DO ipp = 1, nsiz
              ioff = dfft%iss(  gproc*nsiz + ipp  )
              DO i = 1, dfft%nsw(  gproc*nsiz + ipp  )
                  mc = dfft%ismap( i + ioff )
                  it = ii * nppx + gproc * sendsiz
                  DO j = 1, npp
                     f_in( j + it ) = f_aux( mc + ( j - 1 ) * nnp )
                  ENDDO
                  ii = ii + 1
              ENDDO
           ENDDO

           CALL mpi_isend( f_in( gproc*sendsiz + 1 ), sendsiz, MPI_DOUBLE_COMPLEX, gproc, me, gcomm, sh( gproc+1 ), ierr )
        ENDDO

        DO gproc = 1, nblk
           CALL mpi_irecv( f_aux( (gproc-1)*sendsiz + 1 ), sendsiz, MPI_DOUBLE_COMPLEX, gproc-1, MPI_ANY_TAG, gcomm, rh(gproc), ierr )
        ENDDO

     END IF
     !
     call mpi_waitall( nblk, rh, MPI_STATUSES_IGNORE, ierr )
     !
     offset = 0

     IF( use_tg_ ) THEN
        DO proc = 1, nprocp
           gproc = dfft%nplist(proc) + 1
           kdest = ( proc - 1 ) * sendsiz
           kfrom = offset 
           DO k = 1, ncp_ (me)
              DO i = 1, npp_ ( gproc )  
                 scra ( kfrom + i ) = f_aux ( kdest + i )
              ENDDO
              kdest = kdest + nppx
              kfrom = kfrom + nr3x
           ENDDO
           offset = offset + npp_ ( gproc )
        ENDDO
     ELSE
        DO proc = 1, nprocp
           kdest = ( proc - 1 ) * sendsiz 
           kfrom = offset 
           DO k = 1, ncp_ (me)
              DO i = 1, npp_ ( proc )  
                 scra ( kfrom + i ) = f_aux ( kdest + i )
              ENDDO
              kdest = kdest + nppx
              kfrom = kfrom + nr3x
           ENDDO
           offset = offset + npp_ ( proc )
        ENDDO
     ENDIF

     call mpi_waitall( nblk, sh, MPI_STATUSES_IGNORE, ierr )

     f_in = scra( 1:SIZE(f_in) )

  ENDIF

  CALL stop_clock ('fft_scatter')

#endif

  RETURN

END SUBROUTINE fft_scatter
!
#endif
!


!=----------------------------------------------------------------------=!
   END MODULE scatter_mod_2d
!=----------------------------------------------------------------------=!

