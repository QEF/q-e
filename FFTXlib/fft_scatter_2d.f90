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
   MODULE fft_scatter_2d
!=----------------------------------------------------------------------=!

        USE fft_types, ONLY: fft_type_descriptor
        USE fft_param

        IMPLICIT NONE

        SAVE

        PRIVATE

        PUBLIC :: fft_scatter

!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!
!
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

  INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, i, kdest, kfrom
  INTEGER :: me_p, nppx, mc, j, npp, nnp, ii, it, ip, ioff, sendsiz, ncpx, ipp, nsiz
  !
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
     !
     DO proc = 1, nprocp
        !
        kdest = ( proc - 1 ) * sendsiz
        kfrom = offset
        !
        DO k = 1, ncp_ (me)
           DO i = 1, npp_ ( proc )
              f_aux ( kdest + i ) =  f_in ( kfrom + i )
           ENDDO
           kdest = kdest + nppx
           kfrom = kfrom + nr3x
        ENDDO
        offset = offset + npp_ ( proc )
     ENDDO
     !
     ! step two: communication
     !
     CALL mpi_alltoall (f_aux(1), sendsiz, MPI_DOUBLE_COMPLEX, f_in(1), sendsiz, MPI_DOUBLE_COMPLEX, dfft%comm, ierr)
     !
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
        !
        DO ip = 1, dfft%nproc
           !
           ii = 0
           !
           ioff = dfft%iss( ip )
           !
           DO i = 1, dfft%nsw( ip )
              !
              mc = dfft%ismap( i + ioff )
              !
              it = ii * nppx + ( ip - 1 ) * sendsiz
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
        !
        npp  = dfft%nr3p( me )
        nnp  = dfft%nnp
        !
        !
        DO gproc = 1, dfft%nproc
           !
           ii = 0
           !
           ioff = dfft%iss( gproc )
           !
           DO i = 1, dfft%nsw( gproc )
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
        ENDDO

     END IF
     !
     IF( nprocp == 1 ) GO TO 20
     !
     !  step two: communication
     !
     CALL mpi_alltoall (f_in(1), sendsiz, MPI_DOUBLE_COMPLEX, f_aux(1), sendsiz, MPI_DOUBLE_COMPLEX, dfft%comm, ierr)
     !
     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
     !
     !  step one: store contiguously the columns
     !
     !! f_in = 0.0_DP
     !
     offset = 0
     !
     DO proc = 1, nprocp
        !
        kdest = ( proc - 1 ) * sendsiz
        kfrom = offset
        !
        DO k = 1, ncp_ (me)
           DO i = 1, npp_ ( proc )
              f_in ( kfrom + i ) = f_aux ( kdest + i )
           ENDDO
           kdest = kdest + nppx
           kfrom = kfrom + nr3x
        ENDDO
        offset = offset + npp_ ( proc )
     ENDDO

20   CONTINUE

  ENDIF

  CALL stop_clock ('fft_scatter')

#endif

  RETURN

END SUBROUTINE fft_scatter
!
!=----------------------------------------------------------------------=!
   END MODULE fft_scatter_2d
!=----------------------------------------------------------------------=!

