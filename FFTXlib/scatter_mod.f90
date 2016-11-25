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
   MODULE scatter_mod
!=----------------------------------------------------------------------=!

        USE fft_types, ONLY: fft_type_descriptor
        USE task_groups, ONLY: task_groups_descriptor
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

        PUBLIC :: fft_type_descriptor
        PUBLIC :: fft_scatter, gather_grid, scatter_grid
        PUBLIC :: cgather_sym, cgather_sym_many, cscatter_sym_many
        PUBLIC :: maps_sticks_to_3d

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
SUBROUTINE fft_scatter ( dfft, f_in, nr3x, nxx_, f_aux, ncp_, npp_, isgn, dtgs )
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
  !  If optional argument "dtgs" is present the subroutines performs
  !  the trasposition using the Task Groups distribution
  !
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in) :: dfft
  INTEGER, INTENT(in)           :: nr3x, nxx_, isgn, ncp_ (:), npp_ (:)
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)
  TYPE (task_groups_descriptor), OPTIONAL, INTENT(in) :: dtgs

#if defined(__MPI)

  INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
  INTEGER :: me_p, nppx, mc, j, npp, nnp, ii, it, ip, ioff, sendsiz, ncpx, ipp, nblk, nsiz
  !
  LOGICAL :: use_tg_

  !
  !  Task Groups

  use_tg_ = .false.

  IF( present( dtgs ) ) use_tg_ = .true.

  me     = dfft%mype + 1
  !
  IF( use_tg_ ) THEN
    !  This is the number of procs. in the plane-wave group
     nprocp = dtgs%npgrp
  ELSE
     nprocp = dfft%nproc
  ENDIF
  !
  CALL start_clock ('fft_scatter')
  !
  ncpx = 0
  nppx = 0
  IF( use_tg_ ) THEN
     ncpx   = dtgs%tg_ncpx
     nppx   = dtgs%tg_nppx
     gcomm  = dtgs%pgrp_comm
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
     offset = 0

     DO proc = 1, nprocp
        IF( use_tg_ ) THEN
           gproc = dtgs%nplist(proc)+1
        ELSE
           gproc = proc
        ENDIF
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
     ! maybe useless; ensures that no garbage is present in the output
     !
     !! f_in = 0.0_DP
     !
     ! step two: communication
     !
     IF( use_tg_ ) THEN
        gcomm = dtgs%pgrp_comm
     ELSE
        gcomm = dfft%comm
     ENDIF

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
              DO j = 1, dfft%npp( me )
                 f_aux( mc + ( j - 1 ) * dfft%nnp ) = f_in( j + it )
              ENDDO
              it = it + nppx
           ENDDO
        ENDDO

     ELSE

        IF( use_tg_ ) THEN
           npp  = dtgs%tg_npp( me )
           nnp  = dfft%nr1x * dfft%nr2x
        ELSE
           npp  = dfft%npp( me )
           nnp  = dfft%nnp
        ENDIF

        IF( use_tg_ ) THEN
           nblk = dfft%nproc / dtgs%nogrp
           nsiz = dtgs%nogrp
        ELSE
           nblk = dfft%nproc 
           nsiz = 1
        END IF
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
              DO j = 1, dfft%npp( me )
                 f_in( j + it ) = f_aux( mc + ( j - 1 ) * dfft%nnp )
              ENDDO
              it = it + nppx
           ENDDO
        ENDDO

     ELSE

        IF( use_tg_ ) THEN
           npp  = dtgs%tg_npp( me )
           nnp  = dfft%nr1x * dfft%nr2x
        ELSE
           npp  = dfft%npp( me )
           nnp  = dfft%nnp
        ENDIF

        IF( use_tg_ ) THEN
           nblk = dtgs%nproc / dtgs%nogrp
           nsiz = dtgs%nogrp
        ELSE
           nblk = dfft%nproc 
           nsiz = 1
        END IF
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
     IF( use_tg_ ) THEN
        gcomm = dtgs%pgrp_comm
     ELSE
        gcomm = dfft%comm
     ENDIF

     CALL mpi_alltoall (f_in(1), sendsiz, MPI_DOUBLE_COMPLEX, f_aux(1), sendsiz, MPI_DOUBLE_COMPLEX, gcomm, ierr)

     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
     !
     !  step one: store contiguously the columns
     !
     !! f_in = 0.0_DP
     !
     offset = 0

     DO proc = 1, nprocp
        IF( use_tg_ ) THEN
           gproc = dtgs%nplist(proc)+1
        ELSE
           gproc = proc
        ENDIF
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
!   NON BLOCKING SCATTER, should be better on switched network
!   like infiniband, ethernet, myrinet
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter ( dfft, f_in, nr3x, nxx_, f_aux, ncp_, npp_, isgn, dtgs )
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
  !  If optional argument "dtgs" is present the subroutines performs
  !  the trasposition using the Task Groups distribution
  !
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in) :: dfft
  INTEGER, INTENT(in)           :: nr3x, nxx_, isgn, ncp_ (:), npp_ (:)
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)
  TYPE (task_groups_descriptor), OPTIONAL, INTENT(in) :: dtgs

#if defined(__MPI)

  INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
  INTEGER :: me_p, nppx, mc, j, npp, nnp, ii, it, ip, ioff, sendsiz, ncpx, ipp, nblk, nsiz, ijp
  INTEGER :: sh(dfft%nproc), rh(dfft%nproc)
  INTEGER :: istat( MPI_STATUS_SIZE )
  !
  INTEGER, SAVE, ALLOCATABLE :: indmap(:,:)
  INTEGER, SAVE, ALLOCATABLE :: indmap_bw(:)
  INTEGER, SAVE  :: nijp
  INTEGER, SAVE  :: dimref(4) = 0
  INTEGER, SAVE  :: dimref_bw(4) = 0
  !
  LOGICAL :: use_tg_
  !
  CALL start_clock ('fft_scatter')

  use_tg_ = .false.

  IF( present( dtgs ) ) use_tg_ = .true.

  me     = dfft%mype + 1
  !
  ncpx = 0
  nppx = 0
  IF( use_tg_ ) THEN
     !  This is the number of procs. in the plane-wave group
     nprocp = dtgs%npgrp
     ncpx   = dtgs%tg_ncpx
     nppx   = dtgs%tg_nppx
     gcomm  = dtgs%pgrp_comm
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
           gproc = dtgs%nplist(proc)+1
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
     ! maybe useless; ensures that no garbage is present in the output
     !
     !f_in( nprocp*sendsiz + 1 : size( f_in )  ) = 0.0_DP
     !
     call mpi_waitall( nprocp, sh, MPI_STATUSES_IGNORE, ierr )
     !
     f_aux = (0.d0, 0.d0)
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
                 f_aux( mc + ( j - 1 ) * dfft%nnp ) = f_in( j + it )
              ENDDO
              it = it + nppx
           ENDDO
        ENDDO

     ELSE

        IF( use_tg_ ) THEN
           npp  = dtgs%tg_npp( me )
           nnp  = dfft%nr1x * dfft%nr2x
           nblk = dtgs%nproc / dtgs%nogrp
           nsiz = dtgs%nogrp
        ELSE
           npp  = dfft%npp( me )
           nnp  = dfft%nnp
           nblk = dfft%nproc 
           nsiz = 1
        ENDIF
        !
        IF( ( dimref(1) .ne. npp ) .or. ( dimref(2) .ne. nnp ) .or. &
            ( dimref(3) .ne. nblk ) .or. ( dimref(4) .ne. nsiz ) ) THEN
           !
           IF( ALLOCATED( indmap ) )  &
              DEALLOCATE( indmap )
           ALLOCATE( indmap(2,SIZE(f_aux)) )
           !
           ijp = 0
           !
           DO gproc = 1, nblk
              ii = 0
              DO ipp = 1, nsiz
                 ioff = dfft%iss( (gproc-1)*nsiz + ipp )
                 DO i = 1, dfft%nsw( (gproc-1)*nsiz + ipp )
                    mc = dfft%ismap( i + ioff )
                    it = ii * nppx + (gproc-1) * sendsiz
                    DO j = 1, npp
                       ijp = ijp + 1
                       indmap(1,ijp) = mc + ( j - 1 ) * nnp
                       indmap(2,ijp) = j + it 
                    ENDDO
                    ii = ii + 1
                 ENDDO
              ENDDO
           ENDDO
           !
           nijp = ijp
           CALL fftsort( nijp, indmap )
           dimref(1) = npp
           dimref(2) = nnp
           dimref(3) = nblk
           dimref(4) = nsiz
           !
        END IF
        !
        DO ijp = 1, nijp
           f_aux( indmap(1,ijp) ) = f_in( indmap(2,ijp) )
        END DO
        !
     END IF

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
           npp  = dtgs%tg_npp( me )
           nnp  = dfft%nr1x * dfft%nr2x
           nblk = dtgs%nproc / dtgs%nogrp
           nsiz = dtgs%nogrp
        ELSE
           npp  = dfft%npp( me )
           nnp  = dfft%nnp
           nblk = dfft%nproc 
           nsiz = 1
        ENDIF
        !
        IF( ( dimref_bw(1) .ne. npp ) .or. ( dimref_bw(2) .ne. nnp ) .or. &
            ( dimref_bw(3) .ne. nblk ) .or. ( dimref_bw(4) .ne. nsiz ) ) THEN
           !
           IF( ALLOCATED( indmap_bw ) )  &
              DEALLOCATE( indmap_bw )
           !
           ALLOCATE( indmap_bw(SIZE(f_aux)) )
           !
           ijp = 0
           DO gproc = 1, nblk
              ii = 0
              DO ipp = 1, nsiz
                 ioff = dfft%iss(  (gproc-1)*nsiz + ipp  )
                 DO i = 1, dfft%nsw(  (gproc-1)*nsiz + ipp  )
                    mc = dfft%ismap( i + ioff )
                    it = ii * nppx + ( gproc - 1 ) * sendsiz
                    DO j = 1, npp
                       indmap_bw( j + ijp ) = mc + ( j - 1 ) * nnp
                    ENDDO
                    ijp = ijp + npp 
                    ii = ii + 1
                 ENDDO
              ENDDO
           ENDDO

           dimref_bw(1) = npp
           dimref_bw(2) = nnp
           dimref_bw(3) = nblk
           dimref_bw(4) = nsiz

        END IF


        ijp = 0
        DO gproc = 0, nblk-1
           ii = gproc * sendsiz
           DO ipp = gproc*nsiz+1, gproc*nsiz+nsiz
              DO i = 1, dfft%nsw(  ipp  )
                 DO j = 1, npp
                    f_in( j + ii ) = f_aux( indmap_bw( j + ijp ) )
                 ENDDO
                 ijp = ijp + npp
                 ii = ii + nppx 
              ENDDO
           ENDDO
           CALL mpi_isend( f_in( gproc*sendsiz + 1 ), sendsiz, MPI_DOUBLE_COMPLEX, gproc, me, gcomm, sh( gproc+1 ), ierr )
           ! IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', ' backward send info<>0', abs(ierr) )
        ENDDO

        DO gproc = 1, nblk
           CALL mpi_irecv( f_aux( (gproc-1)*sendsiz + 1 ), sendsiz, MPI_DOUBLE_COMPLEX, gproc-1, MPI_ANY_TAG, gcomm, rh(gproc), ierr )
           ! IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', ' backward receive info<>0', abs(ierr) )
        ENDDO

     END IF
     !
     call mpi_waitall( nblk, rh, MPI_STATUSES_IGNORE, ierr )
     call mpi_waitall( nblk, sh, MPI_STATUSES_IGNORE, ierr )
     !
     offset = 0

     IF( use_tg_ ) THEN
        DO proc = 1, nprocp
           gproc = dtgs%nplist(proc) + 1
           kdest = ( proc - 1 ) * sendsiz
           kfrom = offset 
           DO k = 1, ncp_ (me)
              DO i = 1, npp_ ( gproc )  
                 f_in ( kfrom + i ) = f_aux ( kdest + i )
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
                 f_in ( kfrom + i ) = f_aux ( kdest + i )
              ENDDO
              kdest = kdest + nppx
              kfrom = kfrom + nr3x
           ENDDO
           offset = offset + npp_ ( proc )
        ENDDO
     ENDIF

  ENDIF

  CALL stop_clock ('fft_scatter')

#endif

  RETURN

END SUBROUTINE fft_scatter
!
#endif
!
!
SUBROUTINE maps_sticks_to_3d( dffts, dtgs, f_in, nxx_, f_aux, isgn )
  !
  ! this subroutine copy sticks stored in 1D array into the 3D array
  ! to be used with 3D FFT. 
  ! This is meant for the use of 3D scalar FFT in parallel build 
  ! once the data have been "rotated" to have a single band in a single task 
  !
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in) :: dffts
  TYPE (task_groups_descriptor), INTENT(in) :: dtgs
  INTEGER, INTENT(in)           :: nxx_, isgn
  COMPLEX (DP), INTENT(in)      :: f_in (nxx_)
  COMPLEX (DP), INTENT(out)     :: f_aux (nxx_)

  INTEGER :: ijp, ii, i, j, it, ioff, ipp, mc, jj, ip, gproc, nr12x
  !
  f_aux = 0.0d0
  !
  IF( isgn == 2 ) THEN
     ip = 1
     nr12x = dffts%nr1x * dffts%nr2x
     DO gproc = 1, dtgs%nproc / dtgs%nogrp
        ii = 0
        DO ipp = 1, dtgs%nogrp
           ioff = dffts%iss( ip )
           DO i = 1, dffts%nsw( ip )
              mc = dffts%ismap( i + ioff )
              it = ( ii + ( gproc - 1 ) * dtgs%tg_ncpx ) * dtgs%tg_nppx
              DO j = 1, dtgs%tg_npp( dffts%mype + 1 )
                 f_aux( mc + ( j - 1 ) * nr12x ) = f_in( j + it )
              ENDDO
              ii = ii + 1
           ENDDO
           ip = ip + 1
        ENDDO
     ENDDO
  ELSE
     CALL fftx_error__ (' maps_sticks_to_3d ', ' isgn .ne. 2  not implemented ', 999 )
  END IF
  RETURN
END SUBROUTINE maps_sticks_to_3d
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
  INTEGER :: proc, info
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfft%nproc-1), recvcount(0:dfft%nproc-1)
  !
  IF( size( f_in ) < dfft%nnr ) &
     CALL fftx_error__( ' gather_grid ', ' f_in too small ', dfft%nnr-size( f_in ) )
  !
  CALL start_clock( 'gather_grid' )
  !
  DO proc = 0, ( dfft%nproc - 1 )
     !
     recvcount(proc) = dfft%nnp * dfft%npp(proc+1)
     IF ( proc == 0 ) THEN
        displs(proc) = 0
     ELSE
        displs(proc) = displs(proc-1) + recvcount(proc-1)
     ENDIF
     !
  ENDDO
  !
  ! ... the following check should be performed only on processor dfft%root
  ! ... otherwise f_out must be allocated on all processors even if not used
  !
  info = size( f_out ) - displs( dfft%nproc-1 ) - recvcount( dfft%nproc-1 )
  IF( info < 0 ) &
     CALL fftx_error__( ' gather_grid ', ' f_out too small ', -info )
  !
  info = 0
  !
  CALL MPI_GATHERV( f_in, recvcount(dfft%mype), MPI_DOUBLE_PRECISION, f_out, &
                    recvcount, displs, MPI_DOUBLE_PRECISION, dfft%root,      &
                    dfft%comm, info )
  !
  CALL fftx_error__( 'gather_grid', 'info<>0', info )
  !
  CALL stop_clock( 'gather_grid' )
  !
#else
  CALL fftx_error__('gather_grid', 'do not use in serial execution', 1)
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
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfft%nproc-1), recvcount(0:dfft%nproc-1)
  !
  IF( 2*size( f_in ) < dfft%nnr ) &
     CALL fftx_error__( ' gather_grid ', ' f_in too small ', dfft%nnr-size( f_in ) )
  !
  CALL start_clock( 'gather_grid' )
  !
  DO proc = 0, ( dfft%nproc - 1 )
     !
     recvcount(proc) = 2*dfft%nnp * dfft%npp(proc+1)
     IF ( proc == 0 ) THEN
        displs(proc) = 0
     ELSE
        displs(proc) = displs(proc-1) + recvcount(proc-1)
     ENDIF
     !
  ENDDO
  !
  ! ... the following check should be performed only on processor dfft%root
  ! ... otherwise f_out must be allocated on all processors even if not used
  !
  info = 2*size( f_out ) - displs( dfft%nproc - 1 ) - recvcount( dfft%nproc-1 )
  IF( info < 0 ) &
     CALL fftx_error__( ' gather_grid ', ' f_out too small ', -info )
  !
  info = 0
  !
  CALL MPI_GATHERV( f_in, recvcount(dfft%mype), MPI_DOUBLE_PRECISION, f_out, &
                    recvcount, displs, MPI_DOUBLE_PRECISION, dfft%root,      &
                    dfft%comm, info )
  !
  CALL fftx_error__( 'gather_grid', 'info<>0', info )
  !
  CALL stop_clock( 'gather_grid' )
  !
#else
  CALL fftx_error__('gather_grid', 'do not use in serial execution', 1)
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
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfft%nproc-1), sendcount(0:dfft%nproc-1)
  !
  IF( size( f_out ) < dfft%nnr ) &
     CALL fftx_error__( ' scatter_grid ', ' f_out too small ', dfft%nnr-size( f_in ) )
  !
  CALL start_clock( 'scatter_grid' )
  !
  DO proc = 0, ( dfft%nproc - 1 )
     !
     sendcount(proc) = dfft%nnp * dfft%npp(proc+1)
     IF ( proc == 0 ) THEN
        displs(proc) = 0
     ELSE
        displs(proc) = displs(proc-1) + sendcount(proc-1)
     ENDIF
     !
  ENDDO
  !
  ! ... the following check should be performed only on processor dfft%root
  ! ... otherwise f_in must be allocated on all processors even if not used
  !
  info = size( f_in ) - displs( dfft%nproc - 1 ) - sendcount( dfft%nproc - 1 )
  IF( info < 0 ) &
     CALL fftx_error__( ' scatter_grid ', ' f_in too small ', -info )
  !
  info = 0
  !
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_out, sendcount(dfft%mype), MPI_DOUBLE_PRECISION, &
                     dfft%root, dfft%comm, info )
  !
  CALL fftx_error__( 'scatter_grid', 'info<>0', info )
  !
  IF ( sendcount(dfft%mype) /= dfft%nnr ) &
     f_out(sendcount(dfft%mype)+1:dfft%nnr) = 0.D0
  !
  CALL stop_clock( 'scatter_grid' )
  !
#else
  CALL fftx_error__('scatter_grid', 'do not use in serial execution', 1)
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
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfft%nproc-1), sendcount(0:dfft%nproc-1)
  !
  IF( 2*size( f_out ) < dfft%nnr ) &
     CALL fftx_error__( ' scatter_grid ', ' f_out too small ', dfft%nnr-size( f_in ) )
  !
  CALL start_clock( 'scatter_grid' )
  !
  DO proc = 0, ( dfft%nproc - 1 )
     !
     sendcount(proc) = 2*dfft%nnp * dfft%npp(proc+1)
     IF ( proc == 0 ) THEN
        displs(proc) = 0
     ELSE
        displs(proc) = displs(proc-1) + sendcount(proc-1)
     ENDIF
     !
  ENDDO
  !
  ! ... the following check should be performed only on processor dfft%root
  ! ... otherwise f_in must be allocated on all processors even if not used
  !
  info = 2*size( f_in ) - displs( dfft%nproc - 1 ) - sendcount( dfft%nproc - 1 )
  IF( info < 0 ) &
     CALL fftx_error__( ' scatter_grid ', ' f_in too small ', -info )
  !
  info = 0
  !
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_out, sendcount(dfft%mype), MPI_DOUBLE_PRECISION, &
                     dfft%root, dfft%comm, info )
  !
  CALL fftx_error__( 'scatter_grid', 'info<>0', info )
  !
  IF ( sendcount(dfft%mype) /= dfft%nnr ) &
     f_out(sendcount(dfft%mype)+1:dfft%nnr) = 0.D0
  !
  CALL stop_clock( 'scatter_grid' )
  !
#else
  CALL fftx_error__('scatter_grid', 'do not use in serial execution', 1)
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
  !
#if defined(__MPI)
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
     IF ( proc == 0 ) THEN
        displs(proc) = 0
     ELSE
        displs(proc) = displs(proc-1) + recvcount(proc-1)
     ENDIF
     !
  ENDDO
  !
  CALL MPI_BARRIER( dfftp%comm, info )
  !
  CALL MPI_ALLGATHERV( f_in, recvcount(dfftp%mype), MPI_DOUBLE_PRECISION, &
                       f_out, recvcount, displs, MPI_DOUBLE_PRECISION, &
                       dfftp%comm, info )
  !
  CALL fftx_error__( 'cgather_sym', 'info<>0', info )
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
  ! ... COMPLEX*16  f_out = gathered variable (nr1x*nr2x*nr3x, 
  !                                             nbnd_proc(dfftp%mype+1))
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfftp
  INTEGER :: nbnd, nbnd_proc(dfftp%nproc), start_nbnd_proc(dfftp%nproc)
  COMPLEX(DP) :: f_in(dfftp%nnr,nbnd)
  COMPLEX(DP) :: f_out(dfftp%nnp*dfftp%nr3x,nbnd_proc(dfftp%mype+1))
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info
  INTEGER :: ibnd, jbnd
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
  CALL MPI_BARRIER( dfftp%comm, info )
  !
  DO proc = 0, dfftp%nproc - 1
     DO ibnd = 1, nbnd_proc(proc+1)
        jbnd = start_nbnd_proc(proc+1) + ibnd - 1
        CALL MPI_GATHERV( f_in(1,jbnd), recvcount(dfftp%mype), &
                        MPI_DOUBLE_PRECISION, f_out(1,ibnd), recvcount, &
                        displs, MPI_DOUBLE_PRECISION, proc, dfftp%comm, info )
     END DO
  END DO
  !
  CALL fftx_error__( 'cgather_sym_many', 'info<>0', info )
  !
!  CALL mp_barrier( dfftp%comm )
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
  ! ... COMPLEX*16  f_in  = gathered variable (nr1x*nr2x*nr3x,
  !                                                nbnd_proc(dfftp%mype+1) )
  ! ... COMPLEX*16  f_out = distributed variable (nrxx)
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfftp
  INTEGER :: nbnd, nbnd_proc(dfftp%nproc), start_nbnd_proc(dfftp%nproc)
  COMPLEX(DP) :: f_in(dfftp%nnp*dfftp%nr3x,nbnd_proc(dfftp%mype+1))
  COMPLEX(DP) :: f_out(dfftp%nnr)
  INTEGER :: target_ibnd
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info
  INTEGER :: displs(0:dfftp%nproc-1), sendcount(0:dfftp%nproc-1)
  INTEGER :: ibnd, jbnd
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
  f_out = (0.0_DP, 0.0_DP)
  !
  CALL MPI_BARRIER( dfftp%comm, info )
  !
  DO proc = 0, dfftp%nproc - 1
     DO ibnd = 1, nbnd_proc(proc+1)
        jbnd = start_nbnd_proc(proc+1) + ibnd - 1
        IF (jbnd==target_ibnd) &
        CALL MPI_SCATTERV( f_in(1,ibnd), sendcount, displs, &
               MPI_DOUBLE_PRECISION, f_out, sendcount(dfftp%mype), &
               MPI_DOUBLE_PRECISION, proc, dfftp%comm, info )
     ENDDO
  ENDDO
  !
  CALL fftx_error__( 'cscatter_sym_many', 'info<>0', info )
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
!
!---------------------------------------------------------------------
subroutine fftsort (n, ia)  
  !---------------------------------------------------------------------
  ! sort an integer array ia(1:n) into ascending order using heapsort algorithm.
  ! n is input, ia is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ia).
  ! in case of equal values in the data array (ia) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  implicit none  
  !-input/output variables
  integer :: n  
  integer :: ia (2,n)  
  !-local variables
  integer :: i, ir, j, l
  integer :: iia(2)  
  ! nothing to order
  if (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  
  ir = n  
10 continue  
  ! still in hiring phase
  if (l.gt.1) then  
     l = l - 1  
     iia(:) = ia (:,l)  
     ! in retirement-promotion phase.
  else  
     ! clear a space at the end of the array
     iia(:) = ia (:,ir)  
     !
     ! retire the top of the heap into it
     ia (:,ir) = ia (:,1)  
     !
     ! decrease the size of the corporation
     ir = ir - 1  
     ! done with the last promotion
     if (ir.eq.1) then  
        ! the least competent worker at all !
        ia (:,1) = iia(:)  
        !
        return  
     endif
  endif
  ! wheter in hiring or promotion phase, we
  i = l  
  ! set up to place iia in its proper level
  j = l + l  
  !
  do while (j.le.ir)  
     if (j.lt.ir) then  
        if (ia (1,j) .lt. ia (1,j + 1) ) then  
           j = j + 1  
        endif
     endif
     ! demote iia
     if (iia(1).lt.ia (1,j) ) then  
        ia (:,i) = ia (:,j)  
        i = j  
        j = j + j  
     else  
        ! set j to terminate do-while loop
        j = ir + 1  
     endif
  enddo
  ia (:,i) = iia(:)  
  goto 10  
  !
end subroutine fftsort

