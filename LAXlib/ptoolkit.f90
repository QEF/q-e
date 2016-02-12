!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!==----------------------------------------------==!
    MODULE parallel_toolkit
!==----------------------------------------------==!

    IMPLICIT NONE
    SAVE

#include "la_param.f90"

    PRIVATE

    PUBLIC :: rep_matmul_drv
    PUBLIC :: zrep_matmul_drv
    PUBLIC :: dsqmdst, dsqmcll, dsqmred, dsqmsym
    PUBLIC :: zsqmdst, zsqmcll, zsqmred, zsqmher

    CONTAINS

! ---------------------------------------------------------------------------------

SUBROUTINE dsqmdst( n, ar, ldar, a, lda, desc )
  !
  !  Double precision SQuare Matrix DiSTribution
  !  This sub. take a replicated square matrix "ar" and distribute it
  !  across processors as described by descriptor "desc"
  !
  USE descriptors
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: ldar
  REAL(DP)            :: ar(ldar,*)  !  matrix to be splitted, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  REAL(DP)            :: a(lda,*)
  TYPE(la_descriptor), INTENT(IN) :: desc
  !
  REAL(DP), PARAMETER :: zero = 0_DP
  !
  INTEGER :: i, j, nr, nc, ic, ir, nx
  !
  IF( desc%active_node <= 0 ) THEN
     RETURN
  END IF

  nx  = desc%nrcx
  ir  = desc%ir
  ic  = desc%ic
  nr  = desc%nr
  nc  = desc%nc

  IF( lda < nx ) &
     CALL lax_error__( " dsqmdst ", " inconsistent dimension lda ", lda )
  IF( n /= desc%n ) &
     CALL lax_error__( " dsqmdst ", " inconsistent dimension n ", n )

  DO j = 1, nc
     DO i = 1, nr
        a( i, j ) = ar( i + ir - 1, j + ic - 1 )
     END DO
     DO i = nr+1, nx
        a( i, j ) = zero
     END DO
  END DO
  DO j = nc + 1, nx
     DO i = 1, nx
        a( i, j ) = zero
     END DO
  END DO

  RETURN

END SUBROUTINE dsqmdst


SUBROUTINE zsqmdst( n, ar, ldar, a, lda, desc )
  !
  !  double complex (Z) SQuare Matrix DiSTribution
  !  This sub. take a replicated square matrix "ar" and distribute it
  !  across processors as described by descriptor "desc"
  !
  USE descriptors
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: ldar
  COMPLEX(DP)         :: ar(ldar,*)  !  matrix to be splitted, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  COMPLEX(DP)         :: a(lda,*)
  TYPE(la_descriptor), INTENT(IN) :: desc
  !
  COMPLEX(DP), PARAMETER :: zero = ( 0_DP , 0_DP )
  !
  INTEGER :: i, j, nr, nc, ic, ir, nx
  !
  IF( desc%active_node <= 0 ) THEN
     RETURN
  END IF

  nx  = desc%nrcx
  ir  = desc%ir
  ic  = desc%ic
  nr  = desc%nr
  nc  = desc%nc

  IF( lda < nx ) &
     CALL lax_error__( " zsqmdst ", " inconsistent dimension lda ", lda )
  IF( n /= desc%n ) &
     CALL lax_error__( " zsqmdst ", " inconsistent dimension n ", n )

  DO j = 1, nc
     DO i = 1, nr
        a( i, j ) = ar( i + ir - 1, j + ic - 1 )
     END DO
     DO i = nr+1, nx
        a( i, j ) = zero
     END DO
  END DO
  DO j = nc + 1, nx
     DO i = 1, nx
        a( i, j ) = zero
     END DO
  END DO

  RETURN

END SUBROUTINE zsqmdst

! ---------------------------------------------------------------------------------

SUBROUTINE dsqmcll( n, a, lda, ar, ldar, desc, comm )
  !
  !  Double precision SQuare Matrix CoLLect
  !  This sub. take a distributed square matrix "a" and collect 
  !  the block assigned to processors into a replicated matrix "ar",
  !  matrix is distributed as described by descriptor desc
  !
  USE descriptors
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: ldar
  REAL(DP)            :: ar(ldar,*)  !  matrix to be merged, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  REAL(DP)            :: a(lda,*)
  TYPE(la_descriptor), INTENT(IN) :: desc
  INTEGER, INTENT(IN) :: comm
  !
  INTEGER :: i, j

#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif

#if defined __MPI
  !
  INTEGER :: np, nx, ipc, ipr, npr, npc, noff
  INTEGER :: ierr, ir, ic, nr, nc

  REAL(DP), ALLOCATABLE :: buf(:,:)
  !
  IF( desc%active_node > 0 ) THEN
     !
     np = desc%npr * desc%npc
     nx = desc%nrcx
     npr = desc%npr
     npc = desc%npc
     !
     IF( desc%myr == 0 .AND. desc%myc == 0 ) THEN
        ALLOCATE( buf( nx, nx * np ) )
     ELSE
        ALLOCATE( buf( 1, 1 ) )
     END IF
     !
     IF( lda /= nx ) &
        CALL lax_error__( " dsqmcll ", " inconsistent dimension lda ", lda )
     !
     IF( desc%n /= n ) &
        CALL lax_error__( " dsqmcll ", " inconsistent dimension n ", n )
     !
     CALL mpi_gather( a, nx*nx, mpi_double_precision, &
                      buf, nx*nx, mpi_double_precision, 0, desc%comm , ierr )
     !
     IF( ierr /= 0 ) &
        CALL lax_error__( " dsqmcll ", " in gather ", ABS( ierr ) )
     !
     IF( desc%myr == 0 .AND. desc%myc == 0 ) THEN
        DO ipc = 1, npc
           CALL descla_local_dims( ic, nc, n, desc%nx, npc, ipc-1 )
           DO ipr = 1, npr
              CALL descla_local_dims( ir, nr, n, desc%nx, npr, ipr-1 )
              noff = ( ipc - 1 + npc * ( ipr - 1 ) ) * nx
              DO j = 1, nc
                 DO i = 1, nr
                    ar( i + ir - 1, j + ic - 1 ) = buf( i, j + noff )
                 END DO
              END DO
           END DO
        END DO
     END IF
     !
     DEALLOCATE( buf )
     !
  END IF
  !
  CALL mpi_bcast( ar,  ldar * n, mpi_double_precision, 0, comm, ierr )   
  !
  IF( ierr /= 0 ) &
     CALL lax_error__( " dsqmcll ", " in bcast ", ABS( ierr ) )

#else

  DO j = 1, n
     DO i = 1, n
        ar( i, j ) = a( i, j )
     END DO
  END DO

#endif

  RETURN
END SUBROUTINE dsqmcll


SUBROUTINE zsqmcll( n, a, lda, ar, ldar, desc, comm )
  !
  !  double complex (Z) SQuare Matrix CoLLect
  !  This sub. take a distributed square matrix "a" and collect 
  !  the block assigned to processors into a replicated matrix "ar",
  !  matrix is distributed as described by descriptor desc
  !
  USE descriptors
  !
  implicit none
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  !
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: ldar
  COMPLEX(DP)         :: ar(ldar,*)  !  matrix to be merged, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  COMPLEX(DP)         :: a(lda,*)
  TYPE(la_descriptor), INTENT(IN) :: desc
  INTEGER, INTENT(IN) :: comm
  !
  INTEGER :: i, j

#if defined __MPI
  !
  INTEGER :: np, nx, ipc, ipr, npr, npc, noff
  INTEGER :: ierr, ir, ic, nr, nc

  COMPLEX(DP), ALLOCATABLE :: buf(:,:)
  !
  IF( desc%active_node > 0 ) THEN
     !
     np = desc%npr * desc%npc 
     nx = desc%nrcx
     npr = desc%npr
     npc = desc%npc
     !
     IF( desc%myr == 0 .AND. desc%myc == 0 ) THEN
        ALLOCATE( buf( nx, nx * np ) )
     ELSE
        ALLOCATE( buf( 1, 1 ) )
     END IF
     !
     IF( lda /= nx ) &
        CALL lax_error__( " zsqmcll ", " inconsistent dimension lda ", lda )
     !
     IF( desc%n /= n ) &
        CALL lax_error__( " zsqmcll ", " inconsistent dimension n ", n )
     !
     CALL mpi_gather( a, nx*nx, mpi_double_complex, &
                      buf, nx*nx, mpi_double_complex, 0, desc%comm , ierr )
     !
     IF( ierr /= 0 ) &
        CALL lax_error__( " zsqmcll ", " in gather ", ABS( ierr ) )
     !
     IF( desc%myr == 0 .AND. desc%myc == 0 ) THEN
        DO ipc = 1, npc
           CALL descla_local_dims( ic, nc, n, desc%nx, npc, ipc-1 )
           DO ipr = 1, npr
              CALL descla_local_dims( ir, nr, n, desc%nx, npr, ipr-1 )
              noff = ( ipc - 1 + npc * ( ipr - 1 ) ) * nx
              DO j = 1, nc
                 DO i = 1, nr
                    ar( i + ir - 1, j + ic - 1 ) = buf( i, j + noff )
                 END DO
              END DO
           END DO
        END DO
     END IF
     !
     DEALLOCATE( buf )
     !
  END IF
  !
  CALL mpi_bcast( ar,  ldar * n, mpi_double_complex, 0, comm, ierr )   
  !
  IF( ierr /= 0 ) &
     CALL lax_error__( " zsqmcll ", " in bcast ", ABS( ierr ) )

#else

  DO j = 1, n
     DO i = 1, n
        ar( i, j ) = a( i, j )
     END DO
  END DO

#endif

  RETURN
END SUBROUTINE zsqmcll


! ---------------------------------------------------------------------------------

SUBROUTINE dsqmwpb( n, a, lda, desc )
   !
   ! Double precision SQuare Matrix WiPe Border subroutine
   !
   USE descriptors
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda
   REAL(DP)            :: a(lda,*)  !  matrix to be redistributed into b
   TYPE(la_descriptor), INTENT(IN) :: desc
   !
   INTEGER :: i, j
   !
   DO j = 1, desc%nc
      DO i = desc%nr + 1, desc%nrcx
         a( i, j ) = 0_DP
      END DO
   END DO
   DO j = desc%nc + 1, desc%nrcx
      DO i = 1, desc%nrcx
         a( i, j ) = 0_DP
      END DO
   END DO
   !
   RETURN
END SUBROUTINE dsqmwpb

! ---------------------------------------------------------------------------------

SUBROUTINE dsqmsym( n, a, lda, desc )
   !
   ! Double precision SQuare Matrix SYMmetrization
   !
   USE descriptors
   !
   IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda
   REAL(DP)            :: a(lda,*) 
   TYPE(la_descriptor), INTENT(IN) :: desc
#if defined __MPI
   INTEGER :: istatus( MPI_STATUS_SIZE )
#endif
   INTEGER :: i, j
   INTEGER :: comm 
   INTEGER :: nr, nc, dest, sreq, ierr, sour
   REAL(DP) :: atmp

#if defined __MPI

   IF( desc%active_node <= 0 ) THEN
      RETURN
   END IF

   IF( n /= desc%n ) &
      CALL lax_error__( " dsqmsym ", " wrong global dim n ", n )
   IF( lda /= desc%nrcx ) &
      CALL lax_error__( " dsqmsym ", " wrong leading dim lda ", lda )

   comm = desc%comm

   nr = desc%nr 
   nc = desc%nc 
   IF( desc%myc == desc%myr ) THEN
      !
      !  diagonal block, procs work locally
      !
      DO j = 1, nc
         DO i = j + 1, nr
            a(i,j) = a(j,i)
         END DO
      END DO
      !
   ELSE IF( desc%myc > desc%myr ) THEN
      !
      !  super diagonal block, procs send the block to sub diag.
      !
      CALL GRID2D_RANK( 'R', desc%npr, desc%npc, &
                             desc%myc, desc%myr, dest )
      CALL mpi_isend( a, lda*lda, MPI_DOUBLE_PRECISION, dest, 1, comm, sreq, ierr )
      !
      IF( ierr /= 0 ) &
         CALL lax_error__( " dsqmsym ", " in isend ", ABS( ierr ) )
      !
   ELSE IF( desc%myc < desc%myr ) THEN
      !
      !  sub diagonal block, procs receive the block from super diag,
      !  then transpose locally
      !
      CALL GRID2D_RANK( 'R', desc%npr, desc%npc, &
                             desc%myc, desc%myr, sour )
      CALL mpi_recv( a, lda*lda, MPI_DOUBLE_PRECISION, sour, 1, comm, istatus, ierr )
      !
      IF( ierr /= 0 ) &
         CALL lax_error__( " dsqmsym ", " in recv ", ABS( ierr ) )
      !
      DO j = 1, lda
         DO i = j + 1, lda
            atmp = a(i,j)
            a(i,j) = a(j,i)
            a(j,i) = atmp
         END DO
      END DO
      !
   END IF

   IF( desc%myc > desc%myr ) THEN
      !
      CALL MPI_Wait( sreq, istatus, ierr )
      !
      IF( ierr /= 0 ) &
         CALL lax_error__( " dsqmsym ", " in wait ", ABS( ierr ) )
      !
   END IF

#else

   DO j = 1, n
      !
      DO i = j + 1, n
         !
         a(i,j) = a(j,i)
         !
      END DO
      !
   END DO

#endif

   RETURN
END SUBROUTINE dsqmsym


SUBROUTINE zsqmher( n, a, lda, desc )
   !
   ! double complex (Z) SQuare Matrix HERmitianize
   !
   USE descriptors
   !
   IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda
   COMPLEX(DP)         :: a(lda,lda) 
   TYPE(la_descriptor), INTENT(IN) :: desc
#if defined __MPI
   INTEGER :: istatus( MPI_STATUS_SIZE )
#endif
   INTEGER :: i, j
   INTEGER :: comm, myid
   INTEGER :: nr, nc, dest, sreq, ierr, sour
   COMPLEX(DP) :: atmp
   COMPLEX(DP), ALLOCATABLE :: tst1(:,:)
   COMPLEX(DP), ALLOCATABLE :: tst2(:,:)

#if defined __MPI

   IF( desc%active_node <= 0 ) THEN
      RETURN
   END IF

   IF( n /= desc%n ) &
      CALL lax_error__( " zsqmsym ", " wrong global dim n ", n )
   IF( lda /= desc%nrcx ) &
      CALL lax_error__( " zsqmsym ", " wrong leading dim lda ", lda )

   comm = desc%comm

   nr = desc%nr 
   nc = desc%nc 
   IF( desc%myc == desc%myr ) THEN
      !
      !  diagonal block, procs work locally
      !
      DO j = 1, nc
         a(j,j) = CMPLX( DBLE( a(j,j) ), 0_DP, KIND=DP )
         DO i = j + 1, nr
            a(i,j) = CONJG( a(j,i) )
         END DO
      END DO
      !
   ELSE IF( desc%myc > desc%myr ) THEN
      !
      !  super diagonal block, procs send the block to sub diag.
      !
      CALL GRID2D_RANK( 'R', desc%npr, desc%npc, &
                             desc%myc, desc%myr, dest )
      CALL mpi_isend( a, lda*lda, MPI_DOUBLE_COMPLEX, dest, 1, comm, sreq, ierr )
      !
      IF( ierr /= 0 ) &
         CALL lax_error__( " zsqmher ", " in mpi_isend ", ABS( ierr ) )
      !
   ELSE IF( desc%myc < desc%myr ) THEN
      !
      !  sub diagonal block, procs receive the block from super diag,
      !  then transpose locally
      !
      CALL GRID2D_RANK( 'R', desc%npr, desc%npc, &
                             desc%myc, desc%myr, sour )
      CALL mpi_recv( a, lda*lda, MPI_DOUBLE_COMPLEX, sour, 1, comm, istatus, ierr )
      !
      IF( ierr /= 0 ) &
         CALL lax_error__( " zsqmher ", " in mpi_recv ", ABS( ierr ) )
      !
      DO j = 1, lda
         DO i = j + 1, lda
            atmp   = a(i,j)
            a(i,j) = a(j,i)
            a(j,i) = atmp
         END DO
      END DO
      DO j = 1, nc
         DO i = 1, nr
            a(i,j)  = CONJG( a(i,j) )
         END DO
      END DO
      !
   END IF

   IF( desc%myc > desc%myr ) THEN
      !
      CALL MPI_Wait( sreq, istatus, ierr )
      !
      IF( ierr /= 0 ) &
         CALL lax_error__( " zsqmher ", " in MPI_Wait ", ABS( ierr ) )
      !
   END IF

#if defined __PIPPO
   CALL MPI_Comm_rank( comm, myid, ierr )
   ALLOCATE( tst1( n, n ) )
   ALLOCATE( tst2( n, n ) )
   tst1 = 0.0d0
   tst2 = 0.0d0
   do j = 1, desc%nc
   do i = 1, desc%nr
      tst1( i + desc%ir - 1, j + desc%ic - 1 ) = a( i , j )
   end do
   end do
   CALL MPI_REDUCE( tst1, tst2, n*n, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, comm, ierr )
   IF( myid == 0 ) THEN
   DO j = 1, n
      !
      IF( tst2(j,j) /=  CMPLX( DBLE( tst2(j,j) ), 0_DP, KIND=DP ) ) &
        WRITE( 4000, * ) j, tst2(j,j)
      !
      DO i = j + 1, n
         !
         IF( tst2(i,j) /= CONJG( tst2(j,i) ) )  WRITE( 4000, * ) i,j, tst2(i,j)
         !
      END DO
      !
   END DO
   END IF
   
   DEALLOCATE( tst1 )
   DEALLOCATE( tst2 )
#endif

#else

   DO j = 1, n
      !
      a(j,j) = CMPLX( DBLE( a(j,j) ), 0_DP, KIND=DP )
      !
      DO i = j + 1, n
         !
         a(i,j) = CONJG( a(j,i) )
         !
      END DO
      !
   END DO

#endif

   RETURN
END SUBROUTINE zsqmher


! ---------------------------------------------------------------------------------


SUBROUTINE dsqmred( na, a, lda, desca, nb, b, ldb, descb )
   !
   ! Double precision SQuare Matrix REDistribution
   ! 
   ! Copy a global "na * na" matrix locally stored in "a",
   !  and distributed as described by "desca", into a larger
   !  global "nb * nb" matrix stored in "b" and distributed
   !  as described in "descb".
   ! 
   ! If you want to read, get prepared for an headache!
   ! Written struggling by Carlo Cavazzoni.
   !
   USE descriptors
   !
   IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
   !
   INTEGER, INTENT(IN) :: na
   INTEGER, INTENT(IN) :: lda
   REAL(DP)            :: a(lda,lda)  !  matrix to be redistributed into b
   TYPE(la_descriptor), INTENT(IN) :: desca
   INTEGER, INTENT(IN) :: nb
   INTEGER, INTENT(IN) :: ldb
   REAL(DP)            :: b(ldb,ldb)
   TYPE(la_descriptor), INTENT(IN) :: descb

   INTEGER :: ipc, ipr, npc, npr
   INTEGER :: ipr_old, ir_old, nr_old, irx_old
   INTEGER :: ipc_old, ic_old, nc_old, icx_old
   INTEGER :: myrow, mycol, ierr, rank
   INTEGER :: col_comm, row_comm, comm, sreq
   INTEGER :: nr_new, ir_new, irx_new, ir, nr, nrtot, irb, ire
   INTEGER :: nc_new, ic_new, icx_new, ic, nc, nctot, icb, ice
   INTEGER :: ib, i, j, myid
   INTEGER :: nrsnd( desca%npr )
   INTEGER :: ncsnd( desca%npr )
   INTEGER :: displ( desca%npr )
   INTEGER :: irb_new( desca%npr )
   INTEGER :: ire_new( desca%npr )
   INTEGER :: icb_new( desca%npr )
   INTEGER :: ice_new( desca%npr )
   REAL(DP), ALLOCATABLE :: buf(:)
   REAL(DP), ALLOCATABLE :: ab(:,:)
   REAL(DP), ALLOCATABLE :: tst1(:,:)
   REAL(DP), ALLOCATABLE :: tst2(:,:)
#if defined __MPI
    INTEGER :: istatus( MPI_STATUS_SIZE )
#endif

   IF( desca%active_node <= 0 ) THEN
      RETURN
   END IF

   ! preliminary consistency checks

   IF( nb < na ) &
      CALL lax_error__( " dsqmred ", " nb < na, this sub. work only with nb >= na ", nb )
   IF( nb /= descb%n ) &
      CALL lax_error__( " dsqmred ", " wrong global dim nb ", nb )
   IF( na /= desca%n ) &
      CALL lax_error__( " dsqmred ", " wrong global dim na ", na )
   IF( ldb /= descb%nrcx ) &
      CALL lax_error__( " dsqmred ", " wrong leading dim ldb ", ldb )
   IF( lda /= desca%nrcx ) &
      CALL lax_error__( " dsqmred ", " wrong leading dim lda ", lda )

   npr   = desca%npr
   myrow = desca%myr
   npc   = desca%npc
   mycol = desca%myc
   comm  = desca%comm

#if defined __MPI

   ! split communicator into row and col communicators

   CALL MPI_Comm_rank( comm, myid, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " dsqmred ", " in MPI_Comm_rank 1 ", ABS( ierr ) )

   CALL MPI_Comm_split( comm, mycol, myrow, col_comm, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " dsqmred ", " in MPI_Comm_split 1 ", ABS( ierr ) )
    
   CALL MPI_Comm_split( comm, myrow, mycol, row_comm, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " dsqmred ", " in MPI_Comm_split 2 ", ABS( ierr ) )

   CALL MPI_Comm_rank( col_comm, rank, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " dsqmred ", " in MPI_Comm_rank 2 ", ABS( ierr ) )
   IF( rank /= myrow ) &
      CALL lax_error__( " dsqmred ", " building col_comm ", rank )

   CALL MPI_Comm_rank( row_comm, rank, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " dsqmred ", " in MPI_Comm_rank 3 ", ABS( ierr ) )
   IF( rank /= mycol ) &
      CALL lax_error__( " dsqmred ", " building row_comm ", rank )

   ALLOCATE( buf( descb%nrcx * descb%nrcx ) )
   ALLOCATE( ab( descb%nrcx, desca%nrcx ) )

   ! write( 3000 + myid, * ) 'na, nb = ', na, nb

   DO j = 1, descb%nc
      DO i = 1, descb%nr
         b( i, j ) = 0.0d0
      END DO
   END DO

   ab = 0.0d0

   ! first redistribute rows, column groups work in parallel

   DO ipr = 1, npr
      !
      CALL descla_local_dims( ir_new, nr_new, nb, descb%nx, npr, ipr-1 )
      !
      irx_new = ir_new + nr_new - 1

      ! write( 3000 + myid, * ) 'ir_new, nr_new, irx_new = ', ir_new, nr_new, irx_new
      !
      DO ipr_old = 1, npr
         !
         CALL descla_local_dims( ir_old, nr_old, na, desca%nx, npr, ipr_old-1 )
         !
         irx_old = ir_old + nr_old - 1
         !
         ! write( 3000 + myid, * ) 'ir_old, nr_old, irx_old = ', ir_old, nr_old, irx_old
         !
         IF( ir_old >= ir_new .AND. ir_old <= irx_new ) THEN
            !
            nrsnd( ipr_old ) = MIN( nr_old, irx_new - ir_old + 1 )
            irb = 1
            ire = nrsnd( ipr_old )
            irb_new( ipr_old ) = ir_old - ir_new + 1
            ire_new( ipr_old ) = irb_new( ipr_old ) + nrsnd( ipr_old ) - 1
            !
         ELSE IF( ir_new >= ir_old .AND. ir_new <= irx_old ) THEN
            !
            nrsnd( ipr_old ) = irx_old - ir_new + 1
            irb = ir_new - ir_old + 1 
            ire = nr_old
            irb_new( ipr_old ) = 1
            ire_new( ipr_old ) = nrsnd( ipr_old )
            !
         ELSE
            nrsnd( ipr_old ) = 0
            irb = 0
            ire = 0
            irb_new( ipr_old ) = 0
            ire_new( ipr_old ) = 0
         END IF
         !
         ! write( 3000 + myid, * ) 'ipr_old, nrsnd            = ', ipr_old, nrsnd( ipr_old )
         ! write( 3000 + myid, * ) 'ipr_old, irb, ire         = ', ipr_old, irb, ire
         ! write( 3000 + myid, * ) 'ipr_old, irb_new, ire_new = ', ipr_old, irb_new( ipr_old ), ire_new( ipr_old )
         !
         IF( ( myrow == ipr_old - 1 ) .AND. ( nrsnd( ipr_old ) > 0 ) ) THEN
            IF(  myrow /= ipr - 1 ) THEN
               ib = 0
               DO j = 1, desca%nc
                  DO i = irb, ire
                     ib = ib + 1
                     buf( ib ) = a( i, j )
                  END DO
               END DO
               CALL mpi_isend( buf, ib, MPI_DOUBLE_PRECISION, ipr-1, ipr, col_comm, sreq, ierr )
               IF( ierr /= 0 ) &
                  CALL lax_error__( " dsqmred ", " in mpi_isend ", ABS( ierr ) )
            ELSE
               DO j = 1, desca%nc
                  ib = irb
                  DO i = irb_new( ipr_old ), ire_new( ipr_old )
                     ab( i, j ) = a( ib, j )
                     ib = ib + 1
                  END DO
               END DO
            END IF
         END IF
         !
         IF( nrsnd( ipr_old ) /= ire - irb + 1 ) &
            CALL lax_error__( " dsqmred ", " somthing wrong with row 1 ", nrsnd( ipr_old ) )
         IF( nrsnd( ipr_old ) /= ire_new( ipr_old ) - irb_new( ipr_old ) + 1 ) &
            CALL lax_error__( " dsqmred ", " somthing wrong with row 2 ", nrsnd( ipr_old ) )
         !
         nrsnd( ipr_old ) = nrsnd( ipr_old ) * desca%nc
         !
      END DO
      !
      IF( myrow == ipr - 1 ) THEN
         DO ipr_old = 1, npr
            IF( nrsnd( ipr_old ) > 0 ) THEN
               IF(  myrow /= ipr_old - 1 ) THEN
                  CALL mpi_recv( buf, nrsnd(ipr_old), MPI_DOUBLE_PRECISION, ipr_old-1, ipr, col_comm, istatus, ierr )
                  IF( ierr /= 0 ) &
                     CALL lax_error__( " dsqmred ", " in mpi_recv ", ABS( ierr ) )
                  CALL mpi_get_count( istatus, MPI_DOUBLE_PRECISION, ib, ierr) 
                  IF( ierr /= 0 ) &
                     CALL lax_error__( " dsqmred ", " in mpi_get_count ", ABS( ierr ) )
                  IF( ib /= nrsnd(ipr_old) ) &
                     CALL lax_error__( " dsqmred ", " somthing wrong with row 3 ", ib )
                  ib = 0
                  DO j = 1, desca%nc
                     DO i = irb_new( ipr_old ), ire_new( ipr_old )
                        ib = ib + 1
                        ab( i, j ) = buf( ib )
                     END DO
                  END DO
               END IF
            END IF
         END DO
      ELSE
         DO ipr_old = 1, npr
            IF( myrow == ipr_old - 1 .AND. nrsnd( ipr_old ) > 0 ) THEN
               CALL MPI_Wait( sreq, istatus, ierr )
               IF( ierr /= 0 ) &
                  CALL lax_error__( " dsqmred ", " in MPI_Wait ", ABS( ierr ) )
            END IF
         END DO
      END IF
      !
   END DO

   ! then redistribute cols, row groups work in parallel

   DO ipc = 1, npc
      !
      CALL descla_local_dims( ic_new, nc_new, nb, descb%nx, npc, ipc-1 )
      !
      icx_new = ic_new + nc_new - 1
      !
      ! write( 3000 + myid, * ) 'ic_new, nc_new, icx_new = ', ic_new, nc_new, icx_new
      !
      DO ipc_old = 1, npc
         !
         CALL descla_local_dims( ic_old, nc_old, na, desca%nx, npc, ipc_old-1 )
         !
         icx_old = ic_old + nc_old - 1
         !
         ! write( 3000 + myid, * ) 'ic_old, nc_old, icx_old = ', ic_old, nc_old, icx_old
         !
         IF( ic_old >= ic_new .AND. ic_old <= icx_new ) THEN
            !
            ncsnd( ipc_old ) = MIN( nc_old, icx_new - ic_old + 1 )
            icb = 1
            ice = ncsnd( ipc_old )
            icb_new( ipc_old ) = ic_old - ic_new + 1
            ice_new( ipc_old ) = icb_new( ipc_old ) + ncsnd( ipc_old ) - 1
            !
         ELSE IF( ic_new >= ic_old .AND. ic_new <= icx_old ) THEN
            !
            ncsnd( ipc_old ) = icx_old - ic_new + 1
            icb = ic_new - ic_old + 1 
            ice = nc_old
            icb_new( ipc_old ) = 1
            ice_new( ipc_old ) = ncsnd( ipc_old )
            !
         ELSE
            ncsnd( ipc_old ) = 0
            icb = 0
            ice = 0
            icb_new( ipc_old ) = 0
            ice_new( ipc_old ) = 0
         END IF
         !
         ! write( 3000 + myid, * ) 'ipc_old, ncsnd            = ', ipc_old, ncsnd( ipc_old )
         ! write( 3000 + myid, * ) 'ipc_old, icb, ice         = ', ipc_old, icb, ice
         ! write( 3000 + myid, * ) 'ipc_old, icb_new, ice_new = ', ipc_old, icb_new( ipc_old ), ice_new( ipc_old )

         IF( ( mycol == ipc_old - 1 ) .AND. ( ncsnd( ipc_old ) > 0 ) ) THEN
            IF(  mycol /= ipc - 1 ) THEN
               ib = 0
               DO j = icb, ice
                  DO i = 1, descb%nrcx
                     ib = ib + 1
                     buf( ib ) = ab( i, j )
                  END DO
               END DO
               CALL mpi_isend( buf, ib, MPI_DOUBLE_PRECISION, ipc-1, ipc, row_comm, sreq, ierr )
               IF( ierr /= 0 ) &
                  CALL lax_error__( " dsqmred ", " in mpi_isend 2 ", ABS( ierr ) )
            ELSE
               ib = icb
               DO j = icb_new( ipc_old ), ice_new( ipc_old )
                  DO i = 1, descb%nrcx
                        b( i, j ) = ab( i, ib )
                  END DO
                  ib = ib + 1
               END DO
            END IF
         END IF

         IF( ncsnd( ipc_old ) /= ice-icb+1 ) &
            CALL lax_error__( " dsqmred ", " somthing wrong with col 1 ", ncsnd( ipc_old ) )
         IF( ncsnd( ipc_old ) /= ice_new( ipc_old ) - icb_new( ipc_old ) + 1 ) &
            CALL lax_error__( " dsqmred ", " somthing wrong with col 2 ", ncsnd( ipc_old ) )
         !
         ncsnd( ipc_old ) = ncsnd( ipc_old ) * descb%nrcx
         !
      END DO
      !
      IF( mycol == ipc - 1 ) THEN
         DO ipc_old = 1, npc
            IF( ncsnd( ipc_old ) > 0 ) THEN
               IF(  mycol /= ipc_old - 1 ) THEN
                  ib = icb_new( ipc_old )
                  CALL mpi_recv( b( 1, ib ), ncsnd(ipc_old), MPI_DOUBLE_PRECISION, ipc_old-1, ipc, row_comm, istatus, ierr )
                  IF( ierr /= 0 ) &
                     CALL lax_error__( " dsqmred ", " in mpi_recv 2 ", ABS( ierr ) )
                  CALL MPI_GET_COUNT( istatus, MPI_DOUBLE_PRECISION, ib, ierr ) 
                  IF( ierr /= 0 ) &
                     CALL lax_error__( " dsqmred ", " in MPI_GET_COUNT 2 ", ABS( ierr ) )
                  IF( ib /= ncsnd(ipc_old) ) &
                     CALL lax_error__( " dsqmred ", " somthing wrong with col 3 ", ib )
               END IF
            END IF
         END DO
      ELSE
         DO ipc_old = 1, npc
            IF( mycol == ipc_old - 1 .AND. ncsnd( ipc_old ) > 0 ) THEN
               CALL MPI_Wait( sreq, istatus, ierr )
               IF( ierr /= 0 ) &
                  CALL lax_error__( " dsqmred ", " in MPI_Wait 2 ", ABS( ierr ) )
            END IF
         END DO
      END IF
      !
   END DO

   DEALLOCATE( ab )
   DEALLOCATE( buf )

   CALL mpi_comm_free( col_comm, ierr )

   IF( ierr /= 0 ) &
      CALL lax_error__( " dsqmred ", " in mpi_comm_free 1 ", ABS( ierr ) )

   CALL mpi_comm_free( row_comm, ierr )

   IF( ierr /= 0 ) &
      CALL lax_error__( " dsqmred ", " in mpi_comm_free 2 ", ABS( ierr ) )

#if defined __PIPPO

   !  this is for debugging, tests through global matrix, if
   !  the two matrix (pre and before the redistribution) coincide.

   ALLOCATE( tst1( nb, nb ) )
   ALLOCATE( tst2( nb, nb ) )
   ALLOCATE( ab( nb, nb ) )

   ab = 0.0d0

   do j = 1, desca%nc
   do i = 1, desca%nr
      ab( i + desca%ir - 1, j + desca%ic - 1 ) = a( i , j )
   end do
   end do

   CALL MPI_REDUCE( ab, tst1, nb*nb, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr )

   ab = 0.0d0

   do j = 1, descb%nc
   do i = 1, descb%nr
      ab( i + descb%ir - 1, j + descb%ic - 1 ) = b( i , j )
   end do
   end do

   CALL MPI_REDUCE( ab, tst2, nb*nb, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr )

   IF( myid == 0 ) THEN
      write( 1000, * ) na, nb, SUM( ABS( tst2 - tst1 ) ) 
   END IF

   DEALLOCATE( ab )
   DEALLOCATE( tst2 )
   DEALLOCATE( tst1 )

#endif

#endif

   RETURN
END SUBROUTINE dsqmred



SUBROUTINE zsqmred( na, a, lda, desca, nb, b, ldb, descb )
   !
   ! double complex (Z) SQuare Matrix REDistribution
   ! 
   ! Copy a global "na * na" matrix locally stored in "a",
   !  and distributed as described by "desca", into a larger
   !  global "nb * nb" matrix stored in "b" and distributed
   !  as described in "descb".
   ! 
   ! If you want to read, get prepared for an headache!
   ! Written struggling by Carlo Cavazzoni.
   !
   USE descriptors
   !
   IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
   !
   INTEGER, INTENT(IN) :: na
   INTEGER, INTENT(IN) :: lda
   COMPLEX(DP)         :: a(lda,lda)  !  matrix to be redistributed into b
   TYPE(la_descriptor), INTENT(IN) :: desca
   INTEGER, INTENT(IN) :: nb
   INTEGER, INTENT(IN) :: ldb
   COMPLEX(DP)         :: b(ldb,ldb)
   TYPE(la_descriptor), INTENT(IN) :: descb

   INTEGER :: ipc, ipr, npc, npr
   INTEGER :: ipr_old, ir_old, nr_old, irx_old
   INTEGER :: ipc_old, ic_old, nc_old, icx_old
   INTEGER :: myrow, mycol, ierr, rank
   INTEGER :: col_comm, row_comm, comm, sreq
   INTEGER :: nr_new, ir_new, irx_new, ir, nr, nrtot, irb, ire
   INTEGER :: nc_new, ic_new, icx_new, ic, nc, nctot, icb, ice
   INTEGER :: ib, i, j, myid
   INTEGER :: nrsnd( desca%npr )
   INTEGER :: ncsnd( desca%npr )
   INTEGER :: displ( desca%npr )
   INTEGER :: irb_new( desca%npr )
   INTEGER :: ire_new( desca%npr )
   INTEGER :: icb_new( desca%npr )
   INTEGER :: ice_new( desca%npr )
   COMPLEX(DP), ALLOCATABLE :: buf(:)
   COMPLEX(DP), ALLOCATABLE :: ab(:,:)
   COMPLEX(DP), ALLOCATABLE :: tst1(:,:)
   COMPLEX(DP), ALLOCATABLE :: tst2(:,:)
#if defined __MPI
    INTEGER :: istatus( MPI_STATUS_SIZE )
#endif

   IF( desca%active_node <= 0 ) THEN
      RETURN
   END IF

   ! preliminary consistency checks

   IF( nb < na ) &
      CALL lax_error__( " zsqmred ", " nb < na, this sub. work only with nb >= na ", nb )
   IF( nb /= descb%n ) &
      CALL lax_error__( " zsqmred ", " wrong global dim nb ", nb )
   IF( na /= desca%n ) &
      CALL lax_error__( " zsqmred ", " wrong global dim na ", na )
   IF( ldb /= descb%nrcx ) &
      CALL lax_error__( " zsqmred ", " wrong leading dim ldb ", ldb )
   IF( lda /= desca%nrcx ) &
      CALL lax_error__( " zsqmred ", " wrong leading dim lda ", lda )

   npr   = desca%npr
   myrow = desca%myr
   npc   = desca%npc
   mycol = desca%myc
   comm  = desca%comm

#if defined __MPI

   ! split communicator into row and col communicators

   CALL MPI_Comm_rank( comm, myid, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " zsqmred ", " in MPI_Comm_rank 1 ", ABS( ierr ) )

   CALL MPI_Comm_split( comm, mycol, myrow, col_comm, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " zsqmred ", " in MPI_Comm_split 1 ", ABS( ierr ) )

   CALL MPI_Comm_split( comm, myrow, mycol, row_comm, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " zsqmred ", " in MPI_Comm_split 2 ", ABS( ierr ) )

   CALL MPI_Comm_rank( col_comm, rank, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " zsqmred ", " in MPI_Comm_rank 2 ", ABS( ierr ) )
   IF( rank /= myrow ) &
      CALL lax_error__( " zsqmred ", " building col_comm ", rank )

   CALL MPI_Comm_rank( row_comm, rank, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " zsqmred ", " in MPI_Comm_rank 3 ", ABS( ierr ) )
   IF( rank /= mycol ) &
      CALL lax_error__( " zsqmred ", " building row_comm ", rank )

   ALLOCATE( buf( descb%nrcx * descb%nrcx ) )
   ALLOCATE( ab( descb%nrcx, desca%nrcx ) )

   DO j = 1, descb%nc
      DO i = 1, descb%nr
         b( i, j ) = ( 0_DP , 0_DP )
      END DO
   END DO

   ab = ( 0_DP , 0_DP )

   ! first redistribute rows, column groups work in parallel

   DO ipr = 1, npr
      !
      CALL descla_local_dims( ir_new, nr_new, nb, descb%nx, npr, ipr-1 )
      !
      irx_new = ir_new + nr_new - 1
      !
      DO ipr_old = 1, npr
         !
         CALL descla_local_dims( ir_old, nr_old, na, desca%nx, npr, ipr_old-1 )
         !
         irx_old = ir_old + nr_old - 1
         !
         IF( ir_old >= ir_new .AND. ir_old <= irx_new ) THEN
            !
            nrsnd( ipr_old ) = MIN( nr_old, irx_new - ir_old + 1 )
            irb = 1
            ire = nrsnd( ipr_old )
            irb_new( ipr_old ) = ir_old - ir_new + 1
            ire_new( ipr_old ) = irb_new( ipr_old ) + nrsnd( ipr_old ) - 1
            !
         ELSE IF( ir_new >= ir_old .AND. ir_new <= irx_old ) THEN
            !
            nrsnd( ipr_old ) = irx_old - ir_new + 1
            irb = ir_new - ir_old + 1 
            ire = nr_old
            irb_new( ipr_old ) = 1
            ire_new( ipr_old ) = nrsnd( ipr_old )
            !
         ELSE
            nrsnd( ipr_old ) = 0
            irb = 0
            ire = 0
            irb_new( ipr_old ) = 0
            ire_new( ipr_old ) = 0
         END IF
         !
         IF( ( myrow == ipr_old - 1 ) .AND. ( nrsnd( ipr_old ) > 0 ) ) THEN
            IF(  myrow /= ipr - 1 ) THEN
               ib = 0
               DO j = 1, desca%nc
                  DO i = irb, ire
                     ib = ib + 1
                     buf( ib ) = a( i, j )
                  END DO
               END DO
               CALL mpi_isend( buf, ib, MPI_DOUBLE_COMPLEX, ipr-1, ipr, col_comm, sreq, ierr )
               IF( ierr /= 0 ) &
                  CALL lax_error__( " zsqmred ", " in mpi_isend 1 ", ABS( ierr ) )
            ELSE
               DO j = 1, desca%nc
                  ib = irb
                  DO i = irb_new( ipr_old ), ire_new( ipr_old )
                     ab( i, j ) = a( ib, j )
                     ib = ib + 1
                  END DO
               END DO
            END IF
         END IF
         !
         IF( nrsnd( ipr_old ) /= ire - irb + 1 ) &
            CALL lax_error__( " zsqmred ", " somthing wrong with row 1 ", nrsnd( ipr_old ) )
         IF( nrsnd( ipr_old ) /= ire_new( ipr_old ) - irb_new( ipr_old ) + 1 ) &
            CALL lax_error__( " zsqmred ", " somthing wrong with row 2 ", nrsnd( ipr_old ) )
         !
         nrsnd( ipr_old ) = nrsnd( ipr_old ) * desca%nc
         !
      END DO
      !
      IF( myrow == ipr - 1 ) THEN
         DO ipr_old = 1, npr
            IF( nrsnd( ipr_old ) > 0 ) THEN
               IF(  myrow /= ipr_old - 1 ) THEN
                  CALL mpi_recv( buf, nrsnd(ipr_old), MPI_DOUBLE_COMPLEX, ipr_old-1, ipr, col_comm, istatus, ierr )
                  IF( ierr /= 0 ) &
                     CALL lax_error__( " zsqmred ", " in mpi_recv 1 ", ABS( ierr ) )
                  CALL MPI_GET_COUNT( istatus, MPI_DOUBLE_COMPLEX, ib, ierr) 
                  IF( ierr /= 0 ) &
                     CALL lax_error__( " zsqmred ", " in MPI_GET_COUNT 1 ", ABS( ierr ) )
                  IF( ib /= nrsnd(ipr_old) ) &
                     CALL lax_error__( " zsqmred ", " somthing wrong with row 3 ", ib )
                  ib = 0
                  DO j = 1, desca%nc
                     DO i = irb_new( ipr_old ), ire_new( ipr_old )
                        ib = ib + 1
                        ab( i, j ) = buf( ib )
                     END DO
                  END DO
               END IF
            END IF
         END DO
      ELSE
         DO ipr_old = 1, npr
            IF( myrow == ipr_old - 1 .AND. nrsnd( ipr_old ) > 0 ) THEN
               CALL MPI_Wait( sreq, istatus, ierr )
               IF( ierr /= 0 ) &
                  CALL lax_error__( " zsqmred ", " in MPI_Wait 1 ", ABS( ierr ) )
            END IF
         END DO
      END IF
      !
   END DO

   ! then redistribute cols, row groups work in parallel

   DO ipc = 1, npc
      !
      CALL descla_local_dims( ic_new, nc_new, nb, descb%nx, npc, ipc-1 )
      !
      icx_new = ic_new + nc_new - 1
      !
      DO ipc_old = 1, npc
         !
         CALL descla_local_dims( ic_old, nc_old, na, desca%nx, npc, ipc_old-1 )
         !
         icx_old = ic_old + nc_old - 1
         !
         IF( ic_old >= ic_new .AND. ic_old <= icx_new ) THEN
            !
            ncsnd( ipc_old ) = MIN( nc_old, icx_new - ic_old + 1 )
            icb = 1
            ice = ncsnd( ipc_old )
            icb_new( ipc_old ) = ic_old - ic_new + 1
            ice_new( ipc_old ) = icb_new( ipc_old ) + ncsnd( ipc_old ) - 1
            !
         ELSE IF( ic_new >= ic_old .AND. ic_new <= icx_old ) THEN
            !
            ncsnd( ipc_old ) = icx_old - ic_new + 1
            icb = ic_new - ic_old + 1 
            ice = nc_old
            icb_new( ipc_old ) = 1
            ice_new( ipc_old ) = ncsnd( ipc_old )
            !
         ELSE
            ncsnd( ipc_old ) = 0
            icb = 0
            ice = 0
            icb_new( ipc_old ) = 0
            ice_new( ipc_old ) = 0
         END IF
         !
         IF( ( mycol == ipc_old - 1 ) .AND. ( ncsnd( ipc_old ) > 0 ) ) THEN
            IF(  mycol /= ipc - 1 ) THEN
               ib = 0
               DO j = icb, ice
                  DO i = 1, descb%nrcx
                     ib = ib + 1
                     buf( ib ) = ab( i, j )
                  END DO
               END DO
               CALL mpi_isend( buf, ib, MPI_DOUBLE_COMPLEX, ipc-1, ipc, row_comm, sreq, ierr )
               IF( ierr /= 0 ) &
                  CALL lax_error__( " zsqmred ", " in mpi_isend 2 ", ABS( ierr ) )
            ELSE
               ib = icb
               DO j = icb_new( ipc_old ), ice_new( ipc_old )
                  DO i = 1, descb%nrcx
                        b( i, j ) = ab( i, ib )
                  END DO
                  ib = ib + 1
               END DO
            END IF
         END IF

         IF( ncsnd( ipc_old ) /= ice-icb+1 ) &
            CALL lax_error__( " zsqmred ", " somthing wrong with col 1 ", ncsnd( ipc_old ) )
         IF( ncsnd( ipc_old ) /= ice_new( ipc_old ) - icb_new( ipc_old ) + 1 ) &
            CALL lax_error__( " zsqmred ", " somthing wrong with col 2 ", ncsnd( ipc_old ) )
         !
         ncsnd( ipc_old ) = ncsnd( ipc_old ) * descb%nrcx
         !
      END DO
      !
      IF( mycol == ipc - 1 ) THEN
         DO ipc_old = 1, npc
            IF( ncsnd( ipc_old ) > 0 ) THEN
               IF(  mycol /= ipc_old - 1 ) THEN
                  ib = icb_new( ipc_old )
                  CALL mpi_recv( b( 1, ib ), ncsnd(ipc_old), MPI_DOUBLE_COMPLEX, ipc_old-1, ipc, row_comm, istatus, ierr )
                  IF( ierr /= 0 ) &
                     CALL lax_error__( " zsqmred ", " in mpi_recv 2 ", ABS( ierr ) )
                  CALL MPI_GET_COUNT( istatus, MPI_DOUBLE_COMPLEX, ib, ierr ) 
                  IF( ierr /= 0 ) &
                     CALL lax_error__( " zsqmred ", " in MPI_GET_COUNT 2 ", ABS( ierr ) )
                  IF( ib /= ncsnd(ipc_old) ) &
                     CALL lax_error__( " zsqmred ", " somthing wrong with col 3 ", ib )
               END IF
            END IF
         END DO
      ELSE
         DO ipc_old = 1, npc
            IF( mycol == ipc_old - 1 .AND. ncsnd( ipc_old ) > 0 ) THEN
               CALL MPI_Wait( sreq, istatus, ierr )
               IF( ierr /= 0 ) &
                  CALL lax_error__( " zsqmred ", " in MPI_Wait 2 ", ABS( ierr ) )
            END IF
         END DO
      END IF
      !
   END DO

   DEALLOCATE( ab )
   DEALLOCATE( buf )

   CALL mpi_comm_free( col_comm, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " zsqmred ", " in mpi_comm_free 1 ", ABS( ierr ) )

   CALL mpi_comm_free( row_comm, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " zsqmred ", " in mpi_comm_free 2 ", ABS( ierr ) )

#if defined __PIPPO

   !  this is for debugging, tests through global matrix, if
   !  the two matrix (pre and before the redistribution) coincide.

   ALLOCATE( tst1( nb, nb ) )
   ALLOCATE( tst2( nb, nb ) )
   ALLOCATE( ab( nb, nb ) )

   ab = 0.0d0

   do j = 1, desca%nc
   do i = 1, desca%nr
      ab( i + desca%ir - 1, j + desca%ic - 1 ) = a( i , j )
   end do
   end do

   CALL MPI_REDUCE( ab, tst1, nb*nb, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, comm, ierr )

   ab = 0.0d0

   do j = 1, descb%nc
   do i = 1, descb%nr
      ab( i + descb%ir - 1, j + descb%ic - 1 ) = b( i , j )
   end do
   end do

   CALL MPI_REDUCE( ab, tst2, nb*nb, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, comm, ierr )

   IF( myid == 0 ) THEN
      write( 4000, * ) na, nb, SUM( ABS( tst2 - tst1 ) ) 
   END IF

   DEALLOCATE( ab )
   DEALLOCATE( tst2 )
   DEALLOCATE( tst1 )

#endif

#endif

   RETURN
END SUBROUTINE zsqmred



! ---------------------------------------------------------------------------------


SUBROUTINE rep_matmul_drv( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, comm )
  !
  !  Parallel matrix multiplication with replicated matrix
  !  written by Carlo Cavazzoni
  !
  implicit none
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  !
  CHARACTER(LEN=1), INTENT(IN) :: transa, transb
  INTEGER, INTENT(IN) :: m, n, k
  REAL(DP), INTENT(IN) :: alpha, beta
  INTEGER, INTENT(IN) :: lda, ldb, ldc
  REAL(DP) :: a(lda,*), b(ldb,*), c(ldc,*)
  INTEGER, INTENT(IN) :: comm
  !
  !  DGEMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
  !
  !     C := ALPHA*OP( A )*OP( B ) + BETA*C,
  !
  !  WHERE  OP( X ) IS ONE OF
  !
  !     OP( X ) = X   OR   OP( X ) = X',
  !
  !  ALPHA AND BETA ARE SCALARS, AND A, B AND C ARE MATRICES, WITH OP( A )
  !  AN M BY K MATRIX,  OP( B )  A  K BY N MATRIX AND  C AN M BY N MATRIX.
  !
  !
  !

#if defined __MPI

  !

  INTEGER :: ME, I, II, J, JJ, IP, SOUR, DEST, INFO, IERR, ioff, ldx
  INTEGER :: NB, IB_S, NB_SOUR, IB_SOUR, IBUF
  INTEGER :: nproc, mpime, q, r

  REAL(DP), ALLOCATABLE :: auxa( : )
  REAL(DP), ALLOCATABLE :: auxc( : )

  !
  ! ... BODY
  !

  CALL MPI_COMM_SIZE(comm, NPROC, IERR)
  CALL MPI_COMM_RANK(comm, MPIME, IERR)

  IF ( NPROC == 1 ) THEN

     !  if there is only one proc no need of using parallel alg.

     CALL dgemm(TRANSA, TRANSB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

     RETURN

  END IF

  ME = MPIME + 1
  Q = INT( m / NPROC )
  R = MOD( m , NPROC )

  ! ... Find out the number of elements in the local block
  !     along "M" first dimension os matrix A

  NB = Q
  IF( ME <= R ) NB = NB + 1

  ! ... Find out the global index of the local first row

  IF( ME <= R ) THEN
     ib_s = (Q+1)*(ME-1) + 1
  ELSE
     ib_s = Q*(ME-1) + R + 1
  END IF

  ldx = m / nproc + 1

  ALLOCATE( auxa( MAX( n, k ) * ldx ) )
  ALLOCATE( auxc( MAX( n, m ) * ldx ) )

  IF( TRANSA == 'N' .OR. TRANSA == 'n' ) THEN
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, k
        DO I = 1, NB
           auxa( ibuf + I ) = A( I + ioff, J )
        END DO
        ibuf = ibuf + ldx
     END DO
  ELSE
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, k
        DO I = 1, NB
           auxa( ibuf + I ) = A( J, I + ioff )
        END DO
        ibuf = ibuf + ldx
     END DO
     !ioff = ib_s - 1
     !call mytranspose( A( 1, ioff + 1 ), lda, auxa(1), ldx, m, nb)
  END IF

  IF( beta /= 0.0_DP ) THEN
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, n
        DO I = 1, NB
           auxc( ibuf + I ) = C( I + ioff, J )
        END DO
        ibuf = ibuf + ldx
     END DO
  END IF

  CALL dgemm( 'N', transb, nb, n, k, alpha, auxa(1), ldx, B, ldb, beta, auxc(1), ldx )

  ! ... Here processors exchange blocks

  DO IP = 0, NPROC-1

     ! ...    Find out the number of elements in the block of processor SOUR

     NB_SOUR = q
     IF( (IP+1) .LE. r ) NB_SOUR = NB_SOUR+1

     ! ...    Find out the global index of the first row owned by SOUR

     IF( (IP+1) .LE. r ) THEN
        ib_sour = (Q+1)*IP + 1
     ELSE
        ib_sour = Q*IP + R + 1
     END IF

     IF( mpime == ip ) auxa(1:n*ldx) = auxc(1:n*ldx)

     CALL MPI_BCAST( auxa(1), ldx*n, mpi_double_precision, ip, comm, IERR)

     IF( ierr /= 0 ) &
        CALL lax_error__( " rep_matmul_drv ", " in MPI_BCAST ", ABS( ierr ) )

     IBUF = 0
     ioff = IB_SOUR - 1
     DO J = 1, N
        DO I = 1, NB_SOUR
           C( I + ioff, J ) = AUXA( IBUF + I )
        END DO
        IBUF = IBUF + ldx
     END DO

  END DO

  DEALLOCATE( auxa, auxc )

#else

     !  if we are not compiling with __MPI this is equivalent to a blas call

     CALL dgemm(TRANSA, TRANSB, m, N, k, alpha, A, lda, B, ldb, beta, C, ldc)

#endif

  RETURN

END SUBROUTINE rep_matmul_drv


SUBROUTINE zrep_matmul_drv( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, comm )
  !
  !  Parallel matrix multiplication with replicated matrix
  !  written by Carlo Cavazzoni
  !
  implicit none
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  !
  CHARACTER(LEN=1), INTENT(IN) :: transa, transb
  INTEGER, INTENT(IN) :: m, n, k
  COMPLEX(DP), INTENT(IN) :: alpha, beta
  INTEGER, INTENT(IN) :: lda, ldb, ldc
  COMPLEX(DP) :: a(lda,*), b(ldb,*), c(ldc,*)
  INTEGER, INTENT(IN) :: comm
  !
  !  DGEMM  PERFORMS ONE OF THE MATRIX-MATRIX OPERATIONS
  !
  !     C := ALPHA*OP( A )*OP( B ) + BETA*C,
  !
  !  WHERE  OP( X ) IS ONE OF
  !
  !     OP( X ) = X   OR   OP( X ) = X',
  !
  !  ALPHA AND BETA ARE SCALARS, AND A, B AND C ARE MATRICES, WITH OP( A )
  !  AN M BY K MATRIX,  OP( B )  A  K BY N MATRIX AND  C AN M BY N MATRIX.
  !
  !
  !

#if defined __MPI

  !

  INTEGER :: ME, I, II, J, JJ, IP, SOUR, DEST, INFO, IERR, ioff, ldx
  INTEGER :: NB, IB_S, NB_SOUR, IB_SOUR, IBUF
  INTEGER :: nproc, mpime, q, r

  COMPLEX(DP), ALLOCATABLE :: auxa( : )
  COMPLEX(DP), ALLOCATABLE :: auxc( : )

  !
  ! ... BODY
  !

  CALL MPI_COMM_SIZE(comm, NPROC, IERR)
  CALL MPI_COMM_RANK(comm, MPIME, IERR)

  IF ( NPROC == 1 ) THEN

     !  if there is only one proc no need of using parallel alg.

     CALL zgemm(TRANSA, TRANSB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

     RETURN

  END IF

  ME = MPIME + 1
  Q = INT( m / NPROC )
  R = MOD( m , NPROC )

  ! ... Find out the number of elements in the local block
  !     along "M" first dimension os matrix A

  NB = Q
  IF( ME <= R ) NB = NB + 1

  ! ... Find out the global index of the local first row

  IF( ME <= R ) THEN
     ib_s = (Q+1)*(ME-1) + 1
  ELSE
     ib_s = Q*(ME-1) + R + 1
  END IF

  ldx = m / nproc + 1

  ALLOCATE( auxa( MAX( n, k ) * ldx ) )
  ALLOCATE( auxc( MAX( n, m ) * ldx ) )

  IF( TRANSA == 'N' .OR. TRANSA == 'n' ) THEN
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, k
        DO I = 1, NB
           auxa( ibuf + I ) = A( I + ioff, J )
        END DO
        ibuf = ibuf + ldx
     END DO
  ELSE
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, k
        DO I = 1, NB
           auxa( ibuf + I ) = CONJG( A( J, I + ioff ) )
        END DO
        ibuf = ibuf + ldx
     END DO
     !ioff = ib_s - 1
     !call mytranspose( A( 1, ioff + 1 ), lda, auxa(1), ldx, m, nb)
  END IF

  IF( beta /= 0.0_DP ) THEN
     ibuf = 0
     ioff = ib_s - 1
     DO J = 1, n
        DO I = 1, NB
           auxc( ibuf + I ) = C( I + ioff, J )
        END DO
        ibuf = ibuf + ldx
     END DO
  END IF

  CALL zgemm( 'N', transb, nb, n, k, alpha, auxa(1), ldx, B, ldb, beta, auxc(1), ldx )

  ! ... Here processors exchange blocks

  DO IP = 0, NPROC-1

     ! ...    Find out the number of elements in the block of processor SOUR

     NB_SOUR = q
     IF( (IP+1) .LE. r ) NB_SOUR = NB_SOUR+1

     ! ...    Find out the global index of the first row owned by SOUR

     IF( (IP+1) .LE. r ) THEN
        ib_sour = (Q+1)*IP + 1
     ELSE
        ib_sour = Q*IP + R + 1
     END IF

     IF( mpime == ip ) auxa(1:n*ldx) = auxc(1:n*ldx)

     CALL MPI_BCAST( auxa(1), ldx*n, mpi_double_complex, ip, comm, IERR)

     IF( ierr /= 0 ) &
        CALL lax_error__( " zrep_matmul_drv ", " in MPI_BCAST ", ABS( ierr ) )

     IBUF = 0
     ioff = IB_SOUR - 1
     DO J = 1, N
        DO I = 1, NB_SOUR
           C( I + ioff, J ) = AUXA( IBUF + I )
        END DO
        IBUF = IBUF + ldx
     END DO

  END DO

  DEALLOCATE( auxa, auxc )

#else

     !  if we are not compiling with __MPI this is equivalent to a blas call

     CALL zgemm(TRANSA, TRANSB, m, N, k, alpha, A, lda, B, ldb, beta, C, ldc)

#endif

  RETURN

END SUBROUTINE zrep_matmul_drv


!==----------------------------------------------==!
END MODULE parallel_toolkit
!==----------------------------------------------==!

!
!
!=----------------------------------------------------------------------------=!
!
!
!  Cannon's algorithms for parallel matrix multiplication
!  written by Carlo Cavazzoni
!  
!
!

SUBROUTINE sqr_mm_cannon( transa, transb, n, alpha, a, lda, b, ldb, beta, c, ldc, desc )
   !
   !  Parallel square matrix multiplication with Cannon's algorithm
   !
   USE descriptors
   !
   IMPLICIT NONE

#include "la_param.f90"
   !
   CHARACTER(LEN=1), INTENT(IN) :: transa, transb
   INTEGER, INTENT(IN) :: n
   REAL(DP), INTENT(IN) :: alpha, beta
   INTEGER, INTENT(IN) :: lda, ldb, ldc
   REAL(DP) :: a(lda,*), b(ldb,*), c(ldc,*)
   TYPE(la_descriptor), INTENT(IN) :: desc
   !
   !  performs one of the matrix-matrix operations
   !
   !     C := ALPHA*OP( A )*OP( B ) + BETA*C,
   !
   !  where  op( x ) is one of
   !
   !     OP( X ) = X   OR   OP( X ) = X',
   !
   !  alpha and beta are scalars, and a, b and c are square matrices
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
   integer :: ierr
   integer :: np
   integer :: i, j, nr, nc, nb, iter, rowid, colid
   logical :: ta, tb
   INTEGER :: comm
   !
   !
   real(DP), allocatable :: bblk(:,:), ablk(:,:)
   !
#if defined (__MPI)
   !
   integer :: istatus( MPI_STATUS_SIZE )
   !
#endif
   !
   IF( desc%active_node < 0 ) THEN
      !
      !  processors not interested in this computation return quickly
      !
      RETURN
      !
   END IF

   IF( n < 1 ) THEN
      RETURN
   END IF

   IF( desc%npr == 1 ) THEN 
      !
      !  quick return if only one processor is used 
      !
      CALL dgemm( TRANSA, TRANSB, n, n, n, alpha, a, lda, b, ldb, beta, c, ldc)
      !
      RETURN
      !
   END IF

   IF( desc%npr /= desc%npc ) &
      CALL lax_error__( ' sqr_mm_cannon ', ' works only with square processor mesh ', 1 )
   !
   !  Retrieve communicator and mesh geometry
   !
   np    = desc%npr
   comm  = desc%comm
   rowid = desc%myr
   colid = desc%myc
   !
   !  Retrieve the size of the local block
   !
   nr    = desc%nr 
   nc    = desc%nc 
   nb    = desc%nrcx
   !
#if defined (__MPI)
   CALL MPI_BARRIER( comm, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " sqr_mm_cannon ", " in MPI_BARRIER ", ABS( ierr ) )
#endif
   !
   allocate( ablk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         ablk( i, j ) = a( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         ablk( i, j ) = 0.0_DP
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         ablk( i, j ) = 0.0_DP
      END DO
   END DO
   !
   !
   allocate( bblk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         bblk( i, j ) = b( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         bblk( i, j ) = 0.0_DP
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         bblk( i, j ) = 0.0_DP
      END DO
   END DO
   !
   !
   ta = ( TRANSA == 'T' .OR. TRANSA == 't' )
   tb = ( TRANSB == 'T' .OR. TRANSB == 't' )
   !
   !  Shift A rowid+1 places to the west
   ! 
   IF( ta ) THEN
      CALL shift_exch_block( ablk, 'W', 1 )
   ELSE
      CALL shift_block( ablk, 'W', rowid+1, 1 )
   END IF
   !
   !  Shift B colid+1 places to the north
   ! 
   IF( tb ) THEN
      CALL shift_exch_block( bblk, 'N', np+1 )
   ELSE
      CALL shift_block( bblk, 'N', colid+1, np+1 )
   END IF
   !
   !  Accumulate on C
   !
   CALL dgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, beta, c, ldc)
   !
   DO iter = 2, np
      !
      !  Shift A 1 places to the east
      ! 
      CALL shift_block( ablk, 'E', 1, iter )
      !
      !  Shift B 1 places to the south
      ! 
      CALL shift_block( bblk, 'S', 1, np+iter )
      !
      !  Accumulate on C
      !
      CALL dgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, 1.0_DP, c, ldc)
      !
   END DO

   deallocate( ablk, bblk )
   
   RETURN

CONTAINS

   SUBROUTINE shift_block( blk, dir, ln, tag )
      !
      !   Block shift 
      !
      IMPLICIT NONE
      REAL(DP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir      ! shift direction
      INTEGER,          INTENT(IN) :: ln       ! shift length
      INTEGER,          INTENT(IN) :: tag      ! communication tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      !
      IF( dir == 'W' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid - ln + np, np )
         icsrc = MOD( colid + ln + np, np )
         !
      ELSE IF( dir == 'E' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid + ln + np, np )
         icsrc = MOD( colid - ln + np, np )
         !
      ELSE IF( dir == 'N' ) THEN

         irdst = MOD( rowid - ln + np, np )
         irsrc = MOD( rowid + ln + np, np )
         icdst = colid
         icsrc = colid

      ELSE IF( dir == 'S' ) THEN

         irdst = MOD( rowid + ln + np, np )
         irsrc = MOD( rowid - ln + np, np )
         icdst = colid
         icsrc = colid

      ELSE

         CALL lax_error__( ' sqr_mm_cannon ', ' unknown shift direction ', 1 )

      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_PRECISION, &
           idest, tag, isour, tag, comm, istatus, ierr)
      IF( ierr /= 0 ) &
         CALL lax_error__( " sqr_mm_cannon ", " in MPI_SENDRECV_REPLACE ", ABS( ierr ) )
      !
#endif
      RETURN
   END SUBROUTINE shift_block

   SUBROUTINE shift_exch_block( blk, dir, tag )
      !
      !   Combined block shift and exchange
      !   only used for the first step
      !
      IMPLICIT NONE
      REAL(DP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir
      INTEGER,          INTENT(IN) :: tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      INTEGER :: icol, irow
      !
      IF( dir == 'W' ) THEN
         !
         icol = rowid
         irow = colid
         !
         irdst = irow
         icdst = MOD( icol - irow-1 + np, np )
         !
         irow = rowid
         icol = MOD( colid + rowid+1 + np, np )
         !
         irsrc = icol
         icsrc = irow
         !
      ELSE IF( dir == 'N' ) THEN
         !
         icol = rowid
         irow = colid
         !
         icdst = icol
         irdst = MOD( irow - icol-1 + np, np )
         !
         irow = MOD( rowid + colid+1 + np, np )
         icol = colid
         !
         irsrc = icol
         icsrc = irow

      ELSE

         CALL lax_error__( ' sqr_mm_cannon ', ' unknown shift_exch direction ', 1 )

      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_PRECISION, &
           idest, tag, isour, tag, comm, istatus, ierr)
      IF( ierr /= 0 ) &
         CALL lax_error__( " sqr_mm_cannon ", " in MPI_SENDRECV_REPLACE 2 ", ABS( ierr ) )
      !
#endif
      RETURN
   END SUBROUTINE shift_exch_block

END SUBROUTINE sqr_mm_cannon


!=----------------------------------------------------------------------------=!

SUBROUTINE sqr_zmm_cannon( transa, transb, n, alpha, a, lda, b, ldb, beta, c, ldc, desc )
   !
   !  Parallel square matrix multiplication with Cannon's algorithm
   !
   USE descriptors
   !
   IMPLICIT NONE
   !
#include "la_param.f90"
   !
   CHARACTER(LEN=1), INTENT(IN) :: transa, transb
   INTEGER, INTENT(IN) :: n
   COMPLEX(DP), INTENT(IN) :: alpha, beta
   INTEGER, INTENT(IN) :: lda, ldb, ldc
   COMPLEX(DP) :: a(lda,*), b(ldb,*), c(ldc,*)
   TYPE(la_descriptor), INTENT(IN) :: desc
   !
   !  performs one of the matrix-matrix operations
   !
   !     C := ALPHA*OP( A )*OP( B ) + BETA*C,
   !
   !  where  op( x ) is one of
   !
   !     OP( X ) = X   OR   OP( X ) = X',
   !
   !  alpha and beta are scalars, and a, b and c are square matrices
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
   INTEGER :: ierr
   INTEGER :: np
   INTEGER :: i, j, nr, nc, nb, iter, rowid, colid
   LOGICAL :: ta, tb
   INTEGER :: comm
   !
   !
   COMPLEX(DP), ALLOCATABLE :: bblk(:,:), ablk(:,:)
   COMPLEX(DP) :: zone = ( 1.0_DP, 0.0_DP )
   COMPLEX(DP) :: zzero = ( 0.0_DP, 0.0_DP )
   !
#if defined (__MPI)
   !
   integer :: istatus( MPI_STATUS_SIZE )
   !
#endif
   !
   IF( desc%active_node < 0 ) THEN
      !
      !  processors not interested in this computation return quickly
      !
      RETURN
      !
   END IF

   IF( n < 1 ) THEN
      RETURN
   END IF

   IF( desc%npr == 1 ) THEN 
      !
      !  quick return if only one processor is used 
      !
      CALL zgemm( TRANSA, TRANSB, n, n, n, alpha, a, lda, b, ldb, beta, c, ldc)
      !
      RETURN
      !
   END IF

   IF( desc%npr /= desc%npc ) &
      CALL lax_error__( ' sqr_zmm_cannon ', ' works only with square processor mesh ', 1 )
   !
   !  Retrieve communicator and mesh geometry
   !
   np    = desc%npr
   comm  = desc%comm
   rowid = desc%myr
   colid = desc%myc
   !
   !  Retrieve the size of the local block
   !
   nr    = desc%nr 
   nc    = desc%nc 
   nb    = desc%nrcx
   !
#if defined (__MPI)
   CALL MPI_BARRIER( comm, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " sqr_zmm_cannon ", " in MPI_BARRIER ", ABS( ierr ) )
#endif
   !
   allocate( ablk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         ablk( i, j ) = a( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         ablk( i, j ) = zzero
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         ablk( i, j ) = zzero
      END DO
   END DO
   !
   !
   allocate( bblk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         bblk( i, j ) = b( i, j )
      END DO
   END DO
   !
   !  Clear memory outside the matrix block
   !
   DO j = nc+1, nb
      DO i = 1, nb
         bblk( i, j ) = zzero
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         bblk( i, j ) = zzero
      END DO
   END DO
   !
   !
   ta = ( TRANSA == 'C' .OR. TRANSA == 'c' )
   tb = ( TRANSB == 'C' .OR. TRANSB == 'c' )
   !
   !  Shift A rowid+1 places to the west
   ! 
   IF( ta ) THEN
      CALL shift_exch_block( ablk, 'W', 1 )
   ELSE
      CALL shift_block( ablk, 'W', rowid+1, 1 )
   END IF
   !
   !  Shift B colid+1 places to the north
   ! 
   IF( tb ) THEN
      CALL shift_exch_block( bblk, 'N', np+1 )
   ELSE
      CALL shift_block( bblk, 'N', colid+1, np+1 )
   END IF
   !
   !  Accumulate on C
   !
   CALL zgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, beta, c, ldc)
   !
   DO iter = 2, np
      !
      !  Shift A 1 places to the east
      ! 
      CALL shift_block( ablk, 'E', 1, iter )
      !
      !  Shift B 1 places to the south
      ! 
      CALL shift_block( bblk, 'S', 1, np+iter )
      !
      !  Accumulate on C
      !
      CALL zgemm( TRANSA, TRANSB, nr, nc, nb, alpha, ablk, nb, bblk, nb, zone, c, ldc)
      !
   END DO

   deallocate( ablk, bblk )
   
   RETURN

CONTAINS

   SUBROUTINE shift_block( blk, dir, ln, tag )
      !
      !   Block shift 
      !
      IMPLICIT NONE
      COMPLEX(DP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir      ! shift direction
      INTEGER,          INTENT(IN) :: ln       ! shift length
      INTEGER,          INTENT(IN) :: tag      ! communication tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      !
      IF( dir == 'W' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid - ln + np, np )
         icsrc = MOD( colid + ln + np, np )
         !
      ELSE IF( dir == 'E' ) THEN
         !
         irdst = rowid
         irsrc = rowid
         icdst = MOD( colid + ln + np, np )
         icsrc = MOD( colid - ln + np, np )
         !
      ELSE IF( dir == 'N' ) THEN

         irdst = MOD( rowid - ln + np, np )
         irsrc = MOD( rowid + ln + np, np )
         icdst = colid
         icsrc = colid

      ELSE IF( dir == 'S' ) THEN

         irdst = MOD( rowid + ln + np, np )
         irsrc = MOD( rowid - ln + np, np )
         icdst = colid
         icsrc = colid

      ELSE

         CALL lax_error__( ' sqr_zmm_cannon ', ' unknown shift direction ', 1 )

      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_COMPLEX, &
           idest, tag, isour, tag, comm, istatus, ierr)
      IF( ierr /= 0 ) &
         CALL lax_error__( " sqr_zmm_cannon ", " in MPI_SENDRECV_REPLACE 1 ", ABS( ierr ) )
      !
#endif
      RETURN
   END SUBROUTINE shift_block
   !
   SUBROUTINE shift_exch_block( blk, dir, tag )
      !
      !   Combined block shift and exchange
      !   only used for the first step
      !
      IMPLICIT NONE
      COMPLEX(DP) :: blk( :, : )
      CHARACTER(LEN=1), INTENT(IN) :: dir
      INTEGER,          INTENT(IN) :: tag
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      INTEGER :: icol, irow
      !
      IF( dir == 'W' ) THEN
         !
         icol = rowid
         irow = colid
         !
         irdst = irow
         icdst = MOD( icol - irow-1 + np, np )
         !
         irow = rowid
         icol = MOD( colid + rowid+1 + np, np )
         !
         irsrc = icol
         icsrc = irow
         !
      ELSE IF( dir == 'N' ) THEN
         !
         icol = rowid
         irow = colid
         !
         icdst = icol
         irdst = MOD( irow - icol-1 + np, np )
         !
         irow = MOD( rowid + colid+1 + np, np )
         icol = colid
         !
         irsrc = icol
         icsrc = irow

      ELSE

         CALL lax_error__( ' sqr_zmm_cannon ', ' unknown shift_exch direction ', 1 )

      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_COMPLEX, &
           idest, tag, isour, tag, comm, istatus, ierr)
      IF( ierr /= 0 ) &
         CALL lax_error__( " sqr_zmm_cannon ", " in MPI_SENDRECV_REPLACE 2 ", ABS( ierr ) )
      !
#endif
      RETURN
   END SUBROUTINE shift_exch_block

END SUBROUTINE sqr_zmm_cannon

!
!
!
!

SUBROUTINE sqr_tr_cannon( n, a, lda, b, ldb, desc )
   !
   !  Parallel square matrix transposition with Cannon's algorithm
   !
   USE descriptors
   !
   IMPLICIT NONE
   !
#include "la_param.f90"
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   REAL(DP)            :: a(lda,*), b(ldb,*)
   TYPE(la_descriptor), INTENT(IN) :: desc
   !
#if defined (__MPI)
   !
   INCLUDE 'mpif.h'
   !
#endif
   !
   INTEGER :: ierr
   INTEGER :: np, rowid, colid
   INTEGER :: i, j, nr, nc, nb
   INTEGER :: comm
   !
   REAL(DP), ALLOCATABLE :: ablk(:,:)
   !
#if defined (__MPI)
   !
   INTEGER :: istatus( MPI_STATUS_SIZE )
   !
#endif
   !
   IF( desc%active_node < 0 ) THEN
      RETURN
   END IF

   IF( n < 1 ) THEN
     RETURN
   END IF

   IF( desc%npr == 1 ) THEN
      CALL mytranspose( a, lda, b, ldb, n, n )
      RETURN
   END IF

   IF( desc%npr /= desc%npc ) &
      CALL lax_error__( ' sqr_tr_cannon ', ' works only with square processor mesh ', 1 )
   IF( n /= desc%n ) &
      CALL lax_error__( ' sqr_tr_cannon ', ' inconsistent size n  ', 1 )
   IF( lda /= desc%nrcx ) &
      CALL lax_error__( ' sqr_tr_cannon ', ' inconsistent size lda  ', 1 )
   IF( ldb /= desc%nrcx ) &
      CALL lax_error__( ' sqr_tr_cannon ', ' inconsistent size ldb  ', 1 )

   comm = desc%comm

   rowid = desc%myr
   colid = desc%myc
   np    = desc%npr
   !
   !  Compute the size of the local block
   !
   nr = desc%nr 
   nc = desc%nc 
   nb = desc%nrcx
   !
   allocate( ablk( nb, nb ) )
   DO j = 1, nc
      DO i = 1, nr
         ablk( i, j ) = a( i, j )
      END DO
   END DO
   DO j = nc+1, nb
      DO i = 1, nb
         ablk( i, j ) = 0.0_DP
      END DO
   END DO
   DO j = 1, nb
      DO i = nr+1, nb
         ablk( i, j ) = 0.0_DP
      END DO
   END DO
   !
   CALL exchange_block( ablk )
   !
#if defined (__MPI)
   CALL MPI_BARRIER( comm, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " sqr_tr_cannon ", " in MPI_BARRIER ", ABS( ierr ) )
#endif
   !
   DO j = 1, nr
      DO i = 1, nc
         b( j, i ) = ablk( i, j )
      END DO
   END DO
   !
   deallocate( ablk )
   
   RETURN

CONTAINS

   SUBROUTINE exchange_block( blk )
      !
      !   Block exchange ( transpose )
      !
      IMPLICIT NONE
      REAL(DP) :: blk( :, : )
      !
      INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
      !
      irdst = colid
      icdst = rowid
      irsrc = colid
      icsrc = rowid
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_PRECISION, &
           idest, np+np+1, isour, np+np+1, comm, istatus, ierr)
      IF( ierr /= 0 ) &
         CALL lax_error__( " sqr_tr_cannon ", " in MPI_SENDRECV_REPLACE ", ABS( ierr ) )
      !
#endif

      RETURN
   END SUBROUTINE


END SUBROUTINE

!

SUBROUTINE redist_row2col( n, a, b, ldx, nx, desc )
   !
   !  redistribute a, array whose second dimension is distributed over processor row,
   !  to obtain b, with the second dim. distributed over processor clolumn 
   !
   USE descriptors
   !
   IMPLICIT NONE
   !
#include "la_param.f90"
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: ldx, nx
   REAL(DP)            :: a(ldx,nx), b(ldx,nx)
   TYPE(la_descriptor), INTENT(IN) :: desc
   !
#if defined (__MPI)
   !
   INCLUDE 'mpif.h'
   !
#endif
   !
   INTEGER :: ierr
   INTEGER :: np, rowid, colid
   INTEGER :: comm
   INTEGER :: icdst, irdst, icsrc, irsrc, idest, isour
   !
#if defined (__MPI)
   !
   INTEGER :: istatus( MPI_STATUS_SIZE )
   !
#endif
   !
   IF( desc%active_node < 0 ) THEN
      RETURN
   END IF

   IF( n < 1 ) THEN
     RETURN
   END IF

   IF( desc%npr == 1 ) THEN
      b = a
      RETURN
   END IF

   IF( desc%npr /= desc%npc ) &
      CALL lax_error__( ' redist_row2col ', ' works only with square processor mesh ', 1 )
   IF( n /= desc%n ) &
      CALL lax_error__( ' redist_row2col ', ' inconsistent size n  ', 1 )
   IF( nx /= desc%nrcx ) &
      CALL lax_error__( ' redist_row2col ', ' inconsistent size lda  ', 1 )

   comm = desc%comm

   rowid = desc%myr
   colid = desc%myc
   np    = desc%npr
   !
   irdst = colid
   icdst = rowid
   irsrc = colid
   icsrc = rowid
   !
   CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
   CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
   !
#if defined (__MPI)
   !
   CALL MPI_BARRIER( comm, ierr )
   IF( ierr /= 0 ) &
      CALL lax_error__( " redist_row2col ", " in MPI_BARRIER ", ABS( ierr ) )
   !
   CALL MPI_SENDRECV(a, ldx*nx, MPI_DOUBLE_PRECISION, idest, np+np+1, &
                     b, ldx*nx, MPI_DOUBLE_PRECISION, isour, np+np+1, comm, istatus, ierr)
   IF( ierr /= 0 ) &
      CALL lax_error__( " redist_row2col ", " in MPI_SENDRECV ", ABS( ierr ) )
   !
#else
   b = a
#endif
   !
   RETURN

END SUBROUTINE redist_row2col

!
!
!

SUBROUTINE cyc2blk_redist( n, a, lda, nca, b, ldb, ncb, desc )
   !
   !  Parallel square matrix redistribution.
   !  A (input) is cyclically distributed by rows across processors
   !  B (output) is distributed by block across 2D processors grid
   !
   USE descriptors
   !
   IMPLICIT NONE
   !
#include "la_param.f90"
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, nca, ldb, ncb
   REAL(DP) :: a( lda, nca ), b( ldb, ncb )
   TYPE(la_descriptor), INTENT(IN) :: desc
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
   integer :: ierr, itag
   integer :: np, ip, me, nproc, comm_a
   integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
   integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
   INTEGER :: me_ortho(2), np_ortho(2)
   !
   real(DP), allocatable :: rcvbuf(:,:,:)
   real(DP), allocatable :: sndbuf(:,:)
   TYPE(la_descriptor) :: ip_desc
   !
   character(len=256) :: msg
   !
#if defined (__MPI)

   IF( desc%active_node < 0 ) THEN
      RETURN
   END IF

   np = desc%npr  !  dimension of the processor mesh
   nb = desc%nrcx    !  leading dimension of the local matrix block
   me = desc%mype   !  my processor id (starting from 0)
   comm_a = desc%comm
   nproc  = desc%npr * desc%npc

   IF( np /= desc%npc ) &
      CALL lax_error__( ' cyc2blk_redist ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL lax_error__( ' cyc2blk_redist ', ' incorrect first argument (n <= 0)', 1 )
   IF( desc%n < nproc ) &
      CALL lax_error__( ' cyc2blk_redist ', ' number of bands < number of proc ', 1 )

   nbuf = (nb/nproc+2) * nb
   !
   ALLOCATE( sndbuf( nb/nproc+2, nb ) )
   ALLOCATE( rcvbuf( nb/nproc+2, nb, nproc ) )
   !
   ! loop over all processors
   ! 
   DO ip = 0, nproc - 1
      !
      ! 2D proc ortho grid sizes
      !
      np_ortho(1) = desc%npr
      np_ortho(2) = desc%npc
      !
      ! compute other processor coordinates me_ortho
      !
      CALL GRID2D_COORDS( 'R', ip, np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2) )
      !
      ! initialize other processor descriptor
      !
      CALL descla_init( ip_desc, desc%n, desc%nx, np_ortho, me_ortho, desc%comm, desc%cntx, 1 )

      IF( ip_desc%nrcx /= nb ) &
         CALL lax_error__( ' cyc2blk_redist ', ' inconsistent block dim nb ', 1 )
      !
      IF( ip_desc%active_node > 0 ) THEN

         ip_nr = ip_desc%nr
         ip_nc = ip_desc%nc
         ip_ir = ip_desc%ir
         ip_ic = ip_desc%ic
         !
         DO j = 1, ip_nc
            jj = j + ip_ic - 1
            il = 1
            DO i = 1, ip_nr
               ii = i + ip_ir - 1
               IF( MOD( ii - 1, nproc ) == me ) THEN
                  CALL check_sndbuf_index()
                  sndbuf( il, j ) = a( ( ii - 1 )/nproc + 1, jj )
                  il = il + 1
               END IF 
            END DO
         END DO
   
      END IF

      CALL mpi_barrier( comm_a, ierr )
     
      CALL mpi_gather( sndbuf, nbuf, mpi_double_precision, &
                       rcvbuf, nbuf, mpi_double_precision, ip, comm_a, ierr )
      IF( ierr /= 0 ) &
         CALL lax_error__( " cyc2blk_redist ", " in mpi_gather ", ABS( ierr ) )

   END DO

   !
   nr  = desc%nr 
   nc  = desc%nc
   ir  = desc%ir
   ic  = desc%ic
   !
   DO ip = 0, nproc - 1
      DO j = 1, nc
         il = 1
         DO i = 1, nr
            ii = i + ir - 1
            IF( MOD( ii - 1, nproc ) == ip ) THEN
               CALL check_rcvbuf_index()
               b( i, j ) = rcvbuf( il, j, ip+1 )
               il = il + 1
            END IF 
         END DO
      END DO
   END DO
   !
   DEALLOCATE( rcvbuf )
   DEALLOCATE( sndbuf )
   
#else

   b( 1:n, 1:n ) = a( 1:n, 1:n )   

#endif

   RETURN

CONTAINS

   SUBROUTINE check_sndbuf_index()
      CHARACTER(LEN=38), SAVE :: msg = ' check_sndbuf_index in cyc2blk_redist '
      IF( j  > SIZE(sndbuf,2) ) CALL lax_error__( msg, ' j > SIZE(sndbuf,2) ', ip+1 )
      IF( il > SIZE(sndbuf,1) ) CALL lax_error__( msg, ' il > SIZE(sndbuf,1) ', ip+1 )
      IF( ( ii - 1 )/nproc + 1 < 1 ) CALL lax_error__( msg, ' ( ii - 1 )/nproc + 1 < 1 ', ip+1 )
      IF( ( ii - 1 )/nproc + 1 > lda ) CALL lax_error__( msg, ' ( ii - 1 )/nproc + 1 > SIZE(a,1) ', ip+1 )
      IF( jj < 1 ) CALL lax_error__( msg, ' jj < 1 ', ip+1 )
      IF( jj > n ) CALL lax_error__( msg, ' jj > n ', ip+1 )
      RETURN
   END SUBROUTINE check_sndbuf_index

   SUBROUTINE check_rcvbuf_index()
      CHARACTER(LEN=38), SAVE :: msg = ' check_rcvbuf_index in cyc2blk_redist '
      IF( i > ldb ) CALL lax_error__( msg, ' i > ldb ', ip+1 )
      IF( j > ldb ) CALL lax_error__( msg, ' j > ldb ', ip+1 )
      IF( j > nb  ) CALL lax_error__( msg, ' j > nb  ', ip+1 )
      IF( il > SIZE( rcvbuf, 1 ) ) CALL lax_error__( msg, ' il too large ', ip+1 )
      RETURN
   END SUBROUTINE check_rcvbuf_index

END SUBROUTINE cyc2blk_redist


SUBROUTINE cyc2blk_zredist( n, a, lda, nca, b, ldb, ncb, desc )
   !
   !  Parallel square matrix redistribution.
   !  A (input) is cyclically distributed by rows across processors
   !  B (output) is distributed by block across 2D processors grid
   !
   USE descriptors
   !
   IMPLICIT NONE
   !
#include "la_param.f90"
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, nca, ldb, ncb
   COMPLEX(DP) :: a( lda, nca ), b( ldb, ncb )
   TYPE(la_descriptor), INTENT(IN) :: desc
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
   integer :: ierr, itag
   integer :: np, ip, me, nproc, comm_a
   integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
   integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
   INTEGER :: me_ortho(2), np_ortho(2)
   !
   COMPLEX(DP), allocatable :: rcvbuf(:,:,:)
   COMPLEX(DP), allocatable :: sndbuf(:,:)
   TYPE(la_descriptor) :: ip_desc
   !
   character(len=256) :: msg
   !
#if defined (__MPI)

   IF( desc%active_node < 0 ) THEN
      RETURN
   END IF

   np = desc%npr  !  dimension of the processor mesh
   nb = desc%nrcx    !  leading dimension of the local matrix block
   me = desc%mype   !  my processor id (starting from 0)
   comm_a = desc%comm
   nproc  = desc%npr * desc%npc

   IF( np /= desc%npc ) &
      CALL lax_error__( ' cyc2blk_zredist ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL lax_error__( ' cyc2blk_zredist ', ' n less or equal zero ', 1 )
   IF( desc%n < nproc ) &
      CALL lax_error__( ' cyc2blk_zredist ', ' nb less than the number of proc ', 1 )
   !
   nbuf = (nb/nproc+2) * nb
   !
   ALLOCATE( sndbuf( nb/nproc+2, nb ) )
   ALLOCATE( rcvbuf( nb/nproc+2, nb, nproc ) )
   
   DO ip = 0, nproc - 1
      !
      ! 2D proc ortho grid sizes
      !
      np_ortho(1) = desc%npr
      np_ortho(2) = desc%npc
      !
      ! compute other processor coordinates me_ortho
      !
      CALL GRID2D_COORDS( 'R', ip, np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2) )
      !
      ! initialize other processor descriptor
      !
      CALL descla_init( ip_desc, desc%n, desc%nx, np_ortho, me_ortho, desc%comm, desc%cntx, 1 )

      ip_nr = ip_desc%nr
      ip_nc = ip_desc%nc
      ip_ir = ip_desc%ir
      ip_ic = ip_desc%ic
      !
      DO j = 1, ip_nc
         jj = j + ip_ic - 1
         il = 1
         DO i = 1, ip_nr
            ii = i + ip_ir - 1
            IF( MOD( ii - 1, nproc ) == me ) THEN
               CALL check_sndbuf_index()
               sndbuf( il, j ) = a( ( ii - 1 )/nproc + 1, jj )
               il = il + 1
            END IF 
         END DO
      END DO
   
      CALL mpi_barrier( comm_a, ierr )

      CALL mpi_gather( sndbuf, nbuf, mpi_double_complex, &
                       rcvbuf, nbuf, mpi_double_complex, ip, comm_a, ierr )
      IF( ierr /= 0 ) &
         CALL lax_error__( " cyc2blk_zredist ", " in mpi_gather ", ABS( ierr ) )
     
   END DO

   !
   nr  = desc%nr 
   nc  = desc%nc
   ir  = desc%ir
   ic  = desc%ic
   !
   DO ip = 0, nproc - 1
      DO j = 1, nc
         il = 1
         DO i = 1, nr
            ii = i + ir - 1
            IF( MOD( ii - 1, nproc ) == ip ) THEN
               CALL check_rcvbuf_index()
               b( i, j ) = rcvbuf( il, j, ip+1 )
               il = il + 1
            END IF 
         END DO
      END DO
   END DO
   !
   !
   DEALLOCATE( rcvbuf )
   DEALLOCATE( sndbuf )
   
#else

   b( 1:n, 1:n ) = a( 1:n, 1:n )   

#endif

   RETURN

CONTAINS

   SUBROUTINE check_sndbuf_index()
      CHARACTER(LEN=40), SAVE :: msg = ' check_sndbuf_index in cyc2blk_zredist  '
      IF( j  > SIZE(sndbuf,2) ) CALL lax_error__( msg, ' j > SIZE(sndbuf,2) ', ip+1 )
      IF( il > SIZE(sndbuf,1) ) CALL lax_error__( msg, ' il > SIZE(sndbuf,1) ', ip+1 )
      IF( ( ii - 1 )/nproc + 1 < 1 ) CALL lax_error__( msg, ' ( ii - 1 )/nproc + 1 < 1 ', ip+1 )
      IF( ( ii - 1 )/nproc + 1 > SIZE(a,1) ) CALL lax_error__( msg, ' ( ii - 1 )/nproc + 1 > SIZE(a,1) ', ip+1 )
      IF( jj < 1 ) CALL lax_error__( msg, ' jj < 1 ', ip+1 )
      IF( jj > n ) CALL lax_error__( msg, ' jj > n ', ip+1 )
      RETURN
   END SUBROUTINE check_sndbuf_index

   SUBROUTINE check_rcvbuf_index()
      CHARACTER(LEN=40), SAVE :: msg = ' check_rcvbuf_index in cyc2blk_zredist  '
      IF( i > ldb ) CALL lax_error__( msg, ' i > ldb ', ip+1 )
      IF( j > ldb ) CALL lax_error__( msg, ' j > ldb ', ip+1 )
      IF( j > nb  ) CALL lax_error__( msg, ' j > nb  ', ip+1 )
      IF( il > SIZE( rcvbuf, 1 ) ) CALL lax_error__( msg, ' il too large ', ip+1 )
      RETURN
   END SUBROUTINE check_rcvbuf_index

END SUBROUTINE cyc2blk_zredist




SUBROUTINE blk2cyc_redist( n, a, lda, nca, b, ldb, ncb, desc )
   !
   !  Parallel square matrix redistribution.
   !  A (output) is cyclically distributed by rows across processors
   !  B (input) is distributed by block across 2D processors grid
   !
   USE descriptors
   !
   IMPLICIT NONE
   !
#include "la_param.f90"
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, nca, ldb, ncb
   REAL(DP) :: a( lda, nca ), b( ldb, ncb )
   TYPE(la_descriptor), INTENT(IN) :: desc
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
   integer :: ierr, itag
   integer :: np, ip, me, comm_a, nproc
   integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
   integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
   INTEGER :: me_ortho(2), np_ortho(2)
   !
   REAL(DP), allocatable :: rcvbuf(:,:,:)
   REAL(DP), allocatable :: sndbuf(:,:)
   TYPE(la_descriptor) :: ip_desc
   !
   character(len=256) :: msg
   !
#if defined (__MPI)

   IF( desc%active_node < 0 ) THEN
      RETURN
   END IF

   np = desc%npr  !  dimension of the processor mesh
   nb = desc%nrcx    !  leading dimension of the local matrix block
   me = desc%mype   !  my processor id (starting from 0)
   comm_a = desc%comm
   nproc  = desc%npr * desc%npc

   IF( np /= desc%npc ) &
      CALL lax_error__( ' blk2cyc_redist ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL lax_error__( ' blk2cyc_redist ', ' n less or equal zero ', 1 )
   IF( desc%n < nproc ) &
      CALL lax_error__( ' blk2cyc_redist ', ' nb less than the number of proc ', 1 )
   !
   nbuf = (nb/nproc+2) * nb
   !
   ALLOCATE( sndbuf( nb/nproc+2, nb ) )
   ALLOCATE( rcvbuf( nb/nproc+2, nb, nproc ) )
   !
   nr  = desc%nr 
   nc  = desc%nc
   ir  = desc%ir
   ic  = desc%ic
   !
   DO ip = 0, nproc - 1
      DO j = 1, nc
         il = 1
         DO i = 1, nr
            ii = i + ir - 1
            IF( MOD( ii - 1, nproc ) == ip ) THEN
               sndbuf( il, j ) = b( i, j )
               il = il + 1
            END IF 
         END DO
      END DO
      CALL mpi_barrier( comm_a, ierr )
      CALL mpi_gather( sndbuf, nbuf, mpi_double_precision, &
                       rcvbuf, nbuf, mpi_double_precision, ip, comm_a, ierr )
      IF( ierr /= 0 ) &
         CALL lax_error__( " blk2cyc_redist ", " in mpi_gather ", ABS( ierr ) )
   END DO
   !
   
   DO ip = 0, nproc - 1
      !
      ! 2D proc ortho grid sizes
      !
      np_ortho(1) = desc%npr
      np_ortho(2) = desc%npc
      !
      ! compute other processor coordinates me_ortho
      !
      CALL GRID2D_COORDS( 'R', ip, np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2) )
      !
      ! initialize other processor descriptor
      !
      CALL descla_init( ip_desc, desc%n, desc%nx, np_ortho, me_ortho, desc%comm, desc%cntx, 1 )
      !
      ip_nr = ip_desc%nr
      ip_nc = ip_desc%nc
      ip_ir = ip_desc%ir
      ip_ic = ip_desc%ic
      !
      DO j = 1, ip_nc
         jj = j + ip_ic - 1
         il = 1
         DO i = 1, ip_nr
            ii = i + ip_ir - 1
            IF( MOD( ii - 1, nproc ) == me ) THEN
               a( ( ii - 1 )/nproc + 1, jj ) = rcvbuf( il, j, ip+1 )
               il = il + 1
            END IF 
         END DO
      END DO
   
   END DO
   !
   DEALLOCATE( rcvbuf )
   DEALLOCATE( sndbuf )
   
#else

   a( 1:n, 1:n ) = b( 1:n, 1:n )   

#endif

   RETURN

END SUBROUTINE blk2cyc_redist


SUBROUTINE blk2cyc_zredist( n, a, lda, nca, b, ldb, ncb, desc )
   !
   !  Parallel square matrix redistribution.
   !  A (output) is cyclically distributed by rows across processors
   !  B (input) is distributed by block across 2D processors grid
   !
   USE descriptors
   !
   IMPLICIT NONE
   !
#include "la_param.f90"
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, nca, ldb, ncb
   COMPLEX(DP) :: a( lda, nca ), b( ldb, ncb )
   TYPE(la_descriptor), INTENT(IN) :: desc
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
   integer :: ierr, itag
   integer :: np, ip, me, comm_a, nproc
   integer :: ip_ir, ip_ic, ip_nr, ip_nc, il, nbuf, ip_irl
   integer :: i, ii, j, jj, nr, nc, nb, nrl, irl, ir, ic
   INTEGER :: me_ortho(2), np_ortho(2)
   !
   COMPLEX(DP), allocatable :: rcvbuf(:,:,:)
   COMPLEX(DP), allocatable :: sndbuf(:,:)
   TYPE(la_descriptor) :: ip_desc
   !
   character(len=256) :: msg
   !
#if defined (__MPI)

   IF( desc%active_node < 0 ) THEN
      RETURN
   END IF

   np = desc%npr  !  dimension of the processor mesh
   nb = desc%nrcx    !  leading dimension of the local matrix block
   me = desc%mype   !  my processor id (starting from 0)
   comm_a = desc%comm
   nproc  = desc%npr * desc%npc

   IF( np /= desc%npc ) &
      CALL lax_error__( ' blk2cyc_zredist ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL lax_error__( ' blk2cyc_zredist ', ' n less or equal zero ', 1 )
   IF( desc%n < nproc ) &
      CALL lax_error__( ' blk2cyc_zredist ', ' nb less than the number of proc ', 1 )
   !
   nbuf = (nb/nproc+2) * nb
   !
   ALLOCATE( sndbuf( nb/nproc+2, nb ) )
   ALLOCATE( rcvbuf( nb/nproc+2, nb, nproc ) )
   !
   nr  = desc%nr 
   nc  = desc%nc
   ir  = desc%ir
   ic  = desc%ic
   !
   DO ip = 0, nproc - 1
      DO j = 1, nc
         il = 1
         DO i = 1, nr
            ii = i + ir - 1
            IF( MOD( ii - 1, nproc ) == ip ) THEN
               sndbuf( il, j ) = b( i, j )
               il = il + 1
            END IF 
         END DO
      END DO
      CALL mpi_barrier( comm_a, ierr )
      CALL mpi_gather( sndbuf, nbuf, mpi_double_complex, &
                       rcvbuf, nbuf, mpi_double_complex, ip, comm_a, ierr )
      IF( ierr /= 0 ) &
         CALL lax_error__( " blk2cyc_zredist ", " in mpi_gather ", ABS( ierr ) )
   END DO
   !
   
   DO ip = 0, nproc - 1
      !
      ! 2D proc ortho grid sizes
      !
      np_ortho(1) = desc%npr
      np_ortho(2) = desc%npc
      !
      ! compute other processor coordinates me_ortho
      !
      CALL GRID2D_COORDS( 'R', ip, np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2) )
      !
      ! initialize other processor descriptor
      !
      CALL descla_init( ip_desc, desc%n, desc%nx, np_ortho, me_ortho, desc%comm, desc%cntx, 1 )
      !
      ip_nr = ip_desc%nr
      ip_nc = ip_desc%nc
      ip_ir = ip_desc%ir
      ip_ic = ip_desc%ic
      !
      DO j = 1, ip_nc
         jj = j + ip_ic - 1
         il = 1
         DO i = 1, ip_nr
            ii = i + ip_ir - 1
            IF( MOD( ii - 1, nproc ) == me ) THEN
               a( ( ii - 1 )/nproc + 1, jj ) = rcvbuf( il, j, ip+1 )
               il = il + 1
            END IF 
         END DO
      END DO
   
   END DO
   !
   DEALLOCATE( rcvbuf )
   DEALLOCATE( sndbuf )
   
#else

   a( 1:n, 1:n ) = b( 1:n, 1:n )   

#endif

   RETURN

END SUBROUTINE blk2cyc_zredist
!
!
!
!  Double Complex and Double Precision Cholesky Factorization of
!  an Hermitan/Symmetric block distributed matrix
!  written by Carlo Cavazzoni
!
!

SUBROUTINE qe_pzpotrf( sll, ldx, n, desc )
   !
   use descriptors
   !
   implicit none
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
#include "la_param.f90"
   !
   integer :: n, ldx
   TYPE(la_descriptor), INTENT(IN) :: desc
   real(DP)  :: one, zero
   complex(DP) :: sll( ldx, ldx ), cone, czero
   integer :: myrow, mycol, ierr
   integer :: jb, info, ib, kb
   integer :: jnr, jir, jic, jnc
   integer :: inr, iir, iic, inc
   integer :: knr, kir, kic, knc
   integer :: nr, nc
   integer :: rcomm, ccomm, color, key, myid, np
   complex(DP), allocatable :: ssnd( :, : ), srcv( :, : )

   one   = 1.0_DP
   cone  = 1.0_DP
   zero  = 0.0_DP
   czero = 0.0_DP

#if defined __MPI

   myrow = desc%myr
   mycol = desc%myc
   myid  = desc%mype
   np    = desc%npr

   IF( desc%npr /= desc%npc ) THEN
      CALL lax_error__( ' pzpotrf ', ' only square grid are allowed ', 1 ) 
   END IF

   IF( ldx /= desc%nrcx ) THEN
      CALL lax_error__( ' pzpotrf ', ' wrong leading dimension ldx ', ldx ) 
   END IF

   nr = desc%nr
   nc = desc%nc

   ALLOCATE( ssnd( ldx, ldx ) )
   ALLOCATE( srcv( ldx, ldx ) )

   DO jb = 1, np
      !
      !    Update and factorize the current diagonal block and test
      !    for non-positive-definiteness.
      !
      CALL descla_local_dims( jir, jnr, n, desc%nx, np, jb-1 )
      !
      !    since we loop on diagonal blocks/procs we have jnc == jnr
      !
      jnc = jnr   
      !
      !    prepare row and colum communicators
      IF( ( myrow >= ( jb-1 ) ) .AND. ( mycol <= ( jb-1 ) ) ) THEN 
          color = mycol
          key   = myrow
      ELSE
          color = np
          key   = myid
      END IF
      !
      CALL mpi_comm_split( desc%comm , color, key, ccomm, ierr )
      IF( ierr /= 0 ) &
         CALL lax_error__( " pzpotrf ", " in mpi_comm_split 1 ", ABS( ierr ) )
      !  
      IF( myrow >= jb-1 .and. mycol <= jb-1 ) THEN
          color = myrow
          key   = mycol
      ELSE
          color = np
          key   = myid
      END IF
      !
      CALL mpi_comm_split( desc%comm, color, key, rcomm, ierr )
      IF( ierr /= 0 ) &
         CALL lax_error__( " pzpotrf ", " in mpi_comm_split 2 ", ABS( ierr ) )
      !
      !    here every process can work independently, then we need a reduce.
      !
      IF( jb > 1 ) THEN
         !
         DO ib = 1, jb - 1
            IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( ib - 1 ) ) ) THEN
               !
               !  this is because only the lover triangle of ssnd will be set to 0 by ZHERK
               !
               ssnd = 0.0_DP  
               !
               !  remember: matrix ssnd is nr*nr, and procs on the diagonale have nr == nc 
               !
               CALL ZHERK( 'L', 'N', nr, nc, -ONE, sll, ldx, zero, ssnd, ldx )
               !
            END IF
         END DO
         IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
            ssnd = sll
         END IF
         !
         IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol <= ( jb - 1 ) ) ) THEN
            !
            !  accumulate on the diagonal block/proc
            !
            CALL mpi_barrier( rcomm, ierr )

            CALL MPI_REDUCE( ssnd, sll, ldx*ldx, MPI_DOUBLE_COMPLEX, MPI_SUM, jb-1, rcomm, ierr )
            IF( ierr /= 0 ) &
               CALL lax_error__( " pzpotrf ", " in MPI_REDUCE 1 ", ABS( ierr ) )
            !
         END IF
         !
      END IF
      !
      ! Only proj ( jb-1, jb-1 ) operates this
      !
      info = 0
      !
      IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
         CALL ZPOTF2( 'L', jnr, sll, ldx, INFO )
         IF( info /= 0 ) &
            CALL lax_error__( " pzpotrf ", " problems computing cholesky decomposition ", ABS( info ) )
      END IF
      !
      IF( ( jb > 1 ) .AND. ( jb < np ) ) THEN
         !
         !           Compute the current block column.
         !
         ! processors ( 1 : jb - 1, jb ) should bcast their blocs
         ! along column to processor ( 1 : jb - 1, jb + 1 : nb )
         !
         IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol < ( jb - 1 ) ) ) THEN
            CALL mpi_barrier( ccomm, ierr )
            CALL mpi_bcast( sll,  ldx*ldx, MPI_DOUBLE_COMPLEX, 0, ccomm, ierr )   
            IF( ierr /= 0 ) &
               CALL lax_error__( " pzpotrf ", " in mpi_bcast 1 ", ABS( ierr ) )
         ELSE IF( ( myrow > ( jb - 1 ) ) .AND. ( mycol < ( jb - 1 ) ) ) THEN
            CALL mpi_barrier( ccomm, ierr )
            CALL mpi_bcast( srcv, ldx*ldx, MPI_DOUBLE_COMPLEX, 0, ccomm, ierr )   
            IF( ierr /= 0 ) &
               CALL lax_error__( " pzpotrf ", " in mpi_bcast 2 ", ABS( ierr ) )
         END IF
         !
         DO ib = jb + 1, np
            CALL descla_local_dims( iir, inr, n, desc%nx, np, ib-1 )
            DO kb = 1, jb - 1
               CALL descla_local_dims( kic, knc, n, desc%nx, np, kb-1 )
               IF( ( myrow == ( ib - 1 ) ) .AND. ( mycol == ( kb - 1 ) ) ) THEN
                  CALL zgemm( 'N', 'C', inr, jnr, knc, -CONE, sll, ldx, srcv, ldx, czero, ssnd, ldx )
               END IF
            END DO
            IF( ( myrow == ( ib - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
               ssnd = sll
            END IF
         END DO
         !
         ! processors ( jb, jb + 1 : nb ) should collect block along row,
         ! from processors  ( 1 : jb - 1, jb + 1 : nb )
         !
         DO kb = jb + 1, np
            IF( ( myrow == ( kb - 1 ) ) .AND. ( mycol <= ( jb - 1 ) ) ) THEN
               IF( ( jb == 1 ) ) THEN
                  IF( mycol == ( jb - 1 ) ) THEN
                     sll = ssnd
                  END IF
               ELSE
                  CALL mpi_barrier( rcomm, ierr )
                  CALL MPI_REDUCE( ssnd, sll, ldx*ldx, MPI_DOUBLE_COMPLEX, MPI_SUM, jb-1, rcomm, ierr )
                  IF( ierr /= 0 ) &
                     CALL lax_error__( " pzpotrf ", " in mpi_reduce 2 ", ABS( ierr ) )
               END IF
            END IF
         END DO
         !
      END IF
      !
      IF( jb < np ) THEN
         !
         ! processor "jb,jb" should broadcast his block to procs ( jb+1 : nb, jb )
         !
         IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
            CALL mpi_barrier( ccomm, ierr )
            CALL mpi_bcast( sll,  ldx*ldx, MPI_DOUBLE_COMPLEX, 0, ccomm, ierr )   
            IF( ierr /= 0 ) &
               CALL lax_error__( " pzpotrf ", " in mpi_bcast 3 ", ABS( ierr ) )
         ELSE IF( ( myrow > ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
            CALL mpi_barrier( ccomm, ierr )
            CALL mpi_bcast( srcv,  ldx*ldx, MPI_DOUBLE_COMPLEX, 0, ccomm, ierr )   
            IF( ierr /= 0 ) &
               CALL lax_error__( " pzpotrf ", " in mpi_bcast 4 ", ABS( ierr ) )
         END IF
         !
         DO ib = jb + 1, np
            IF( ( myrow == ( ib - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
               CALL ZTRSM( 'R', 'L', 'C', 'N', nr, nc, CONE, srcv, ldx, sll, ldx )
            END IF
         END DO
         !
      END IF
      !
      CALL mpi_comm_free( rcomm, ierr )
      IF( ierr /= 0 ) &
         CALL lax_error__( " pzpotrf ", " in mpi_comm_free 1 ", ABS( ierr ) )
      !
      CALL mpi_comm_free( ccomm, ierr )
      IF( ierr /= 0 ) &
         CALL lax_error__( " pzpotrf ", " in mpi_comm_free 2 ", ABS( ierr ) )
      !
   END DO

   DEALLOCATE( srcv, ssnd )

#else

   CALL ZPOTRF( 'L', n, sll, ldx, info )

   IF( info /= 0 ) &
      CALL lax_error__( " pzpotrf ", " problems computing cholesky decomposition ", ABS( info ) )

#endif

   return
END SUBROUTINE qe_pzpotrf

!  now the Double Precision subroutine

SUBROUTINE qe_pdpotrf( sll, ldx, n, desc )
   !
   use descriptors
   !
   implicit none
   !
#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   !
#include "la_param.f90"
   !
   integer  :: n, ldx
   TYPE(la_descriptor), INTENT(IN) :: desc
   REAL(DP) :: one, zero
   REAL(DP) :: sll( ldx, ldx )
   integer  :: myrow, mycol, ierr
   integer  :: jb, info, ib, kb
   integer  :: jnr, jir, jic, jnc
   integer  :: inr, iir, iic, inc
   integer  :: knr, kir, kic, knc
   integer  :: nr, nc
   integer  :: rcomm, ccomm, color, key, myid, np
   REAL(DP), ALLOCATABLE :: ssnd( :, : ), srcv( :, : )

   one   = 1.0_DP
   zero  = 0.0_DP

#if defined __MPI

   myrow = desc%myr
   mycol = desc%myc
   myid  = desc%mype
   np    = desc%npr

   IF( desc%npr /= desc%npc ) THEN
      CALL lax_error__( ' pdpotrf ', ' only square grid are allowed ', 1 ) 
   END IF

   IF( ldx /= desc%nrcx ) THEN
      CALL lax_error__( ' pdpotrf ', ' wrong leading dimension ldx ', ldx ) 
   END IF

   nr = desc%nr
   nc = desc%nc

   ALLOCATE( ssnd( ldx, ldx ) )
   ALLOCATE( srcv( ldx, ldx ) )

   DO jb = 1, np
      !
      !    Update and factorize the current diagonal block and test
      !    for non-positive-definiteness.
      !
      CALL descla_local_dims( jir, jnr, n, desc%nx, np, jb-1 )
      !
      !    since we loop on diagonal blocks/procs we have jnc == jnr
      !
      jnc = jnr   
      !
      !    prepare row and colum communicators
      IF( ( myrow >= ( jb-1 ) ) .AND. ( mycol <= ( jb-1 ) ) ) THEN 
          color = mycol
          key   = myrow
      ELSE
          color = np
          key   = myid
      END IF
      !
      CALL mpi_comm_split( desc%comm , color, key, ccomm, ierr )
      IF( ierr /= 0 ) &
         CALL lax_error__( " pdpotrf ", " in mpi_comm_split 1 ", ABS( ierr ) )
      !  
      IF( myrow >= jb-1 .and. mycol <= jb-1 ) THEN
          color = myrow
          key   = mycol
      ELSE
          color = np
          key   = myid
      END IF
      !
      CALL mpi_comm_split( desc%comm, color, key, rcomm, ierr )
      IF( ierr /= 0 ) &
         CALL lax_error__( " pdpotrf ", " in mpi_comm_split 2 ", ABS( ierr ) )
      !
      !    here every process can work independently, then we need a reduce.
      !
      IF( jb > 1 ) THEN
         !
         DO ib = 1, jb - 1
            IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( ib - 1 ) ) ) THEN
               !
               !  this is because only the lover triangle of ssnd will be set to 0 by ZHERK
               !
               ssnd = 0_DP
               !
               !  remember: matrix ssnd is nr*nr, and procs on the diagonale have nr == nc 
               !
               CALL DSYRK( 'L', 'N', nr, nc, -ONE, sll, ldx, zero, ssnd, ldx )
               !
            END IF
         END DO
         IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
            ssnd = sll
         END IF
         !
         IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol <= ( jb - 1 ) ) ) THEN
            !
            !  accumulate on the diagonal block/proc
            !
            CALL MPI_REDUCE( ssnd, sll, ldx*ldx, MPI_DOUBLE_PRECISION, MPI_SUM, jb-1, rcomm, ierr )
            IF( ierr /= 0 ) &
               CALL lax_error__( " pdpotrf ", " in MPI_REDUCE 1 ", ABS( ierr ) )
            !
         END IF
         !
      END IF
      !
      ! Only proj ( jb-1, jb-1 ) operates this
      !
      info = 0
      !
      IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
         CALL DPOTRF( 'L', jnr, sll, ldx, INFO )
         IF( info /= 0 ) &
            CALL lax_error__( " pdpotrf ", " problems computing cholesky decomposition ", ABS( info ) )
      END IF
      !
      IF( ( jb > 1 ) .AND. ( jb < np ) ) THEN
         !
         !           Compute the current block column.
         !
         ! processors ( 1 : jb - 1, jb ) should bcast their blocs
         ! along column to processor ( 1 : jb - 1, jb + 1 : nb )
         !
         IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol < ( jb - 1 ) ) ) THEN
            CALL mpi_barrier( ccomm, ierr )
            CALL mpi_bcast( sll,  ldx*ldx, MPI_DOUBLE_PRECISION, 0, ccomm, ierr )   
            IF( ierr /= 0 ) &
               CALL lax_error__( " pdpotrf ", " in mpi_bcast 1 ", ABS( ierr ) )
         ELSE IF( ( myrow > ( jb - 1 ) ) .AND. ( mycol < ( jb - 1 ) ) ) THEN
            CALL mpi_barrier( ccomm, ierr )
            CALL mpi_bcast( srcv, ldx*ldx, MPI_DOUBLE_PRECISION, 0, ccomm, ierr )   
            IF( ierr /= 0 ) &
               CALL lax_error__( " pdpotrf ", " in mpi_bcast 2 ", ABS( ierr ) )
         END IF
         !
         DO ib = jb + 1, np
            CALL descla_local_dims( iir, inr, n, desc%nx, np, ib-1 )
            DO kb = 1, jb - 1
               CALL descla_local_dims( kic, knc, n, desc%nx, np, kb-1 )
               IF( ( myrow == ( ib - 1 ) ) .AND. ( mycol == ( kb - 1 ) ) ) THEN
                  CALL dgemm( 'N', 'T', inr, jnr, knc, -ONE, sll, ldx, srcv, ldx, zero, ssnd, ldx )
               END IF
            END DO
            IF( ( myrow == ( ib - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
               ssnd = sll
            END IF
         END DO
         !
         ! processors ( jb, jb + 1 : nb ) should collect block along row,
         ! from processors  ( 1 : jb - 1, jb + 1 : nb )
         !
         DO kb = jb + 1, np
            IF( ( myrow == ( kb - 1 ) ) .AND. ( mycol <= ( jb - 1 ) ) ) THEN
               IF( ( jb == 1 ) ) THEN
                  IF( mycol == ( jb - 1 ) ) THEN
                     sll = ssnd
                  END IF
               ELSE
                  CALL MPI_REDUCE( ssnd, sll, ldx*ldx, MPI_DOUBLE_PRECISION, MPI_SUM, jb-1, rcomm, ierr )
                  IF( ierr /= 0 ) &
                     CALL lax_error__( " pdpotrf ", " in mpi_reduce 2 ", ABS( ierr ) )
               END IF
            END IF
         END DO
         !
      END IF
      !
      IF( jb < np ) THEN
         !
         ! processor "jb,jb" should broadcast his block to procs ( jb+1 : nb, jb )
         !
         IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
            CALL mpi_barrier( ccomm, ierr )
            CALL mpi_bcast( sll,  ldx*ldx, MPI_DOUBLE_PRECISION, 0, ccomm, ierr )   
            IF( ierr /= 0 ) &
               CALL lax_error__( " pdpotrf ", " in mpi_bcast 3 ", ABS( ierr ) )
         ELSE IF( ( myrow > ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
            CALL mpi_barrier( ccomm, ierr )
            CALL mpi_bcast( srcv,  ldx*ldx, MPI_DOUBLE_PRECISION, 0, ccomm, ierr )   
            IF( ierr /= 0 ) &
               CALL lax_error__( " pdpotrf ", " in mpi_bcast 4 ", ABS( ierr ) )
         END IF
         !
         DO ib = jb + 1, np
            IF( ( myrow == ( ib - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
               CALL DTRSM( 'R', 'L', 'T', 'N', nr, nc, ONE, srcv, ldx, sll, ldx )
            END IF
         END DO
         !
      END IF
      !
      CALL mpi_comm_free( rcomm, ierr )
      IF( ierr /= 0 ) &
         CALL lax_error__( " pdpotrf ", " in mpi_comm_free 1 ", ABS( ierr ) )

      CALL mpi_comm_free( ccomm, ierr )
      IF( ierr /= 0 ) &
         CALL lax_error__( " pdpotrf ", " in mpi_comm_free 2 ", ABS( ierr ) )
      !
   END DO

   DEALLOCATE( srcv, ssnd )

#else

   CALL DPOTRF( 'L', n, sll, ldx, info )

   IF( info /= 0 ) &
      CALL lax_error__( " pzpotrf ", " problems computing cholesky decomposition ", ABS( info ) )

#endif

   return
END SUBROUTINE qe_pdpotrf

!
!
!
!

SUBROUTINE qe_pztrtri ( sll, ldx, n, desc )
    
    ! pztrtri computes the parallel inversion of a lower triangular matrix 
    ! distribuited among the processes using a 2-D block partitioning. 
    ! The algorithm is based on the schema below and executes the model 
    ! recursively to each column C2 under the diagonal.     
    !
    !     |-------|-------|      |--------------------|--------------------|
    !     |   A1  |   0   |      |   C1 = trtri(A1)   |          0         |
    ! A = |-------|-------|  C = |--------------------|--------------------|
    !     |   A2  |   A3  |      | C2 = -C3 * A2 * C1 |   C3 = trtri(A3)   | 
    !     |-------|-------|      |--------------------|--------------------|
    !
    ! The recursive steps of multiplication (C2 = -C3 * A2 * C1) is based on the Cannon's algorithms 
    ! for parallel matrix multiplication and is done with BLACS(dgemm)
    !
    !
    ! Arguments
    ! ============
    !
    ! sll   = local block of data
    ! ldx   = leading dimension of one block
    ! n     = size of the global array diributed among the blocks
    ! desc  = descriptor of the matrix distribution
    !
    !
    !  written by Ivan Girotto
    !

    USE descriptors

    IMPLICIT NONE

#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   
#include "la_param.f90"

    INTEGER, INTENT( IN ) :: n, ldx
    TYPE(la_descriptor), INTENT(IN) :: desc
    COMPLEX(DP), INTENT( INOUT ) :: sll( ldx, ldx )

    COMPLEX(DP), PARAMETER :: ONE = (1.0_DP, 0.0_DP)
    COMPLEX(DP), PARAMETER :: ZERO = (0.0_DP, 0.0_DP)

#if defined __MPI
    INTEGER :: status(MPI_STATUS_SIZE)
#endif
    INTEGER :: req(2), ierr, col_comm
    INTEGER :: send, recv, group_rank, group_size
    INTEGER :: myrow, mycol, np, myid, comm

    ! counters
    INTEGER :: k, i, j, count, step_count, shiftcount, cicle 
    INTEGER :: C3dim   ! Dimension of submatrix B
    INTEGER :: nc, nr ! Local dimension of block
    INTEGER :: info, sup_recv
    INTEGER :: idrowref, idcolref, idref, idrecv 

    ! B and BUF_RECV are used to overload the computation of matrix multiplication and the shift of the blocks
    COMPLEX(DP), ALLOCATABLE, DIMENSION( :, : ) :: B, C, BUF_RECV 
    COMPLEX(DP) :: first

    myrow = desc%myr
    mycol = desc%myc
    myid  = desc%mype
    np    = desc%npr
    comm  = desc%comm

    IF( desc%npr /= desc%npc ) THEN
       CALL lax_error__( ' pztrtri ', ' only square grid are allowed ', 1 ) 
    END IF
    IF( ldx /= desc%nrcx ) THEN
       CALL lax_error__( ' pztrtri ', ' wrong leading dimension ldx ', ldx ) 
    END IF

    nr = desc%nr
    nc = desc%nc

    !  clear elements outside local meaningful block nr*nc

    DO j = nc+1, ldx
       DO i = 1, ldx
          sll( i, j ) = zero
       END DO
    END DO
    DO j = 1, ldx
       DO i = nr+1, ldx
          sll( i, j ) = zero
       END DO
    END DO

#if defined __MPI

    ALLOCATE( B( ldx, ldx ) )
    ALLOCATE( C( ldx, ldx ) )
    ALLOCATE( BUF_RECV ( ldx, ldx ) )

    IF( np == 2 ) THEN
       !
       !  special case with 4 proc, 2x2 grid
       !
       IF( myrow == mycol ) THEN
          CALL compute_ztrtri()
       END IF
       !
       CALL GRID2D_RANK( 'R', np, np, 1, 0, idref )
       !
       IF( myrow == 0 .AND. mycol == 0 ) THEN
          CALL MPI_Send(sll, ldx*ldx, MPI_DOUBLE_COMPLEX, idref, 0, comm, ierr)
          IF( ierr /= 0 ) &
             CALL lax_error__( " pztrtri ", " in mpi_send 1 ", ABS( ierr ) )
       END IF
       !
       IF( myrow == 1 .AND. mycol == 1 ) THEN
          CALL MPI_Send(sll, ldx*ldx, MPI_DOUBLE_COMPLEX, idref, 1, comm, ierr)
          IF( ierr /= 0 ) &
             CALL lax_error__( " pztrtri ", " in mpi_send 2 ", ABS( ierr ) )
       END IF
       !
       IF( myrow == 1 .AND. mycol == 0 ) THEN
          !
          CALL GRID2D_RANK( 'R', np, np, 0, 0, i )
          CALL GRID2D_RANK( 'R', np, np, 1, 1, j )
          !
          CALL MPI_Irecv( B, ldx*ldx, MPI_DOUBLE_COMPLEX, i, 0, comm, req(1), ierr)
          IF( ierr /= 0 ) &
             CALL lax_error__( " pztrtri ", " in mpi_irecv 3 ", ABS( ierr ) )
          !
          CALL MPI_Irecv( C, ldx*ldx, MPI_DOUBLE_COMPLEX, j, 1, comm, req(2), ierr)
          IF( ierr /= 0 ) &
             CALL lax_error__( " pztrtri ", " in mpi_irecv 4 ", ABS( ierr ) )
          !
          CALL MPI_Wait(req(1), status, ierr)
          IF( ierr /= 0 ) &
             CALL lax_error__( " pztrtri ", " in MPI_Wait 5 ", ABS( ierr ) )
          !
          CALL zgemm('N', 'N', ldx, ldx, ldx, ONE, sll, ldx, b, ldx, ZERO, buf_recv, ldx)
          !
          CALL MPI_Wait(req(2), status, ierr)
          IF( ierr /= 0 ) &
             CALL lax_error__( " pztrtri ", " in MPI_Wait 6 ", ABS( ierr ) )
          !
          CALL zgemm('N', 'N', ldx, ldx, ldx, -ONE, c, ldx, buf_recv, ldx, ZERO, sll, ldx)
          !
       END IF
       !
       IF( myrow == 0 .AND. mycol == 1 ) THEN
          !
          sll = zero
          !
       END IF
       !
       DEALLOCATE( b, c, buf_recv )
       !
       RETURN
       !
    END IF

  
    IF( myrow >= mycol ) THEN
       !
       !  only procs on lower triangle partecipates
       !
       CALL MPI_Comm_split( comm, mycol, myrow, col_comm, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pztrtri ", " in MPI_Comm_split 9 ", ABS( ierr ) )

       CALL MPI_Comm_size( col_comm, group_size, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pztrtri ", " in MPI_Comm_size 10 ", ABS( ierr ) )
       !
       CALL MPI_Comm_rank( col_comm, group_rank, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pztrtri ", " in MPI_Comm_rank 11 ", ABS( ierr ) )
       !
    ELSE
       !
       !  other procs stay at the window!
       !
       CALL MPI_Comm_split( comm, MPI_UNDEFINED, MPI_UNDEFINED, col_comm, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pztrtri ", " in MPI_Comm_split 12 ", ABS( ierr ) )
       !
       sll = zero
       !
    END IF
    ! 

    ! Compute the inverse of a lower triangular 
    ! along the diagonal of the global array with BLAS(ztrtri) 
    !
    IF( mycol == myrow ) THEN
       !
       CALL compute_ztrtri()
       !
    ELSE IF( myrow > mycol ) THEN
       !
       buf_recv = sll
       !
    END IF

    IF( myrow >= mycol ) THEN
       !
       ! Broadcast the diagonal blocks to the processors under the diagonal
       !
       CALL MPI_Bcast( sll, ldx*ldx, MPI_DOUBLE_COMPLEX, 0, col_comm, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pztrtri ", " in MPI_Bcast 13 ", ABS( ierr ) )
       !
    END IF

    ! Compute A2 * C1 and start the Cannon's algorithm shifting the blocks of column one place to the North
    !
    IF( myrow > mycol ) THEN
       !
       CALL zgemm( 'N', 'N', ldx, ldx, ldx, ONE, buf_recv, ldx, sll, ldx, ZERO, c, ldx )
       !
       send = shift( 1, group_rank, 1, ( group_size - 1 ), 'N' )
       recv = shift( 1, group_rank, 1, ( group_size - 1 ), 'S' )
       !
       CALL MPI_Sendrecv( c, ldx*ldx, MPI_DOUBLE_COMPLEX, send, 0, buf_recv, &
            ldx*ldx, MPI_DOUBLE_COMPLEX, recv, 0, col_comm, status, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pztrtri ", " in MPI_Sendrecv 14 ", ABS( ierr ) )
       !
    END IF

    ! Execute the Cannon's algorithm to compute ricorsively the multiplication of C2 = -C3 * A2 * C1
    !
    DO count = ( np - 2 ), 0, -1
       C3dim = (np-1) - count ! Dimension of the submatrix C3
       first = ZERO
       cicle = 0
       IF( ( myrow > count ) .AND. ( mycol >= count ) ) THEN
          idcolref = count + 1
          idrowref = myrow
          CALL GRID2D_RANK( 'R', np, np, idrowref, idcolref, idref )
          idrecv = idref - 1
          ! Compute C2 = -C3 * A2 * C1
          DO shiftcount = count, np-2
             IF(mycol>count)THEN
                ! Execute the virtual shift of the matrix C3 along the row in order to know which processor 
                ! have to send the block to C2 
                IF( cicle == 0)THEN
                   ! virtual shift of the block i,j of the submatrix C3 i place to West
                   send = shift(idref, myid, myrow-count, C3dim, 'W')
                ELSE
                   ! virtual shift of the block i,j of the submatrix C3 i place to West
                   send = shift(idref, send, 1, C3dim, 'E')
                END IF
                IF(send==idref)THEN
                   CALL MPI_Send(sll, ldx*ldx, MPI_DOUBLE_COMPLEX, idrecv, myid, comm, ierr)
                   IF( ierr /= 0 ) &
                      CALL lax_error__( " pztrtri ", " in MPI_Send 15 ", ABS( ierr ) )
                END IF
             ELSE
                IF( cicle == 0)THEN
                   ! virtual shift of the block i,j of the submatrix C3 i place to West
                   sup_recv = shift(idref, myid+1, myrow-count, C3dim, 'E')
                ELSE
                   ! virtual shift of the block i,j of the submatrix C3 i place to West
                   sup_recv = shift(idref, sup_recv, 1, C3dim, 'W')
                END IF
                CALL MPI_Recv(C, ldx*ldx, MPI_DOUBLE_COMPLEX, sup_recv, sup_recv, comm, status, ierr)
                IF( ierr /= 0 ) &
                   CALL lax_error__( " pztrtri ", " in MPI_Recv 16 ", ABS( ierr ) )
                send = shift(1, group_rank, 1, (group_size-1), 'S')
                recv = shift(1, group_rank, 1, (group_size-1), 'N')
                ! with the no-blocking communication the computation and the shift of the column block are overapped  
                !
                IF( MOD( cicle, 2 ) == 0 ) THEN
                   CALL MPI_Isend(BUF_RECV, ldx*ldx, MPI_DOUBLE_COMPLEX, send, group_rank+cicle, col_comm, req(1), ierr)
                   IF( ierr /= 0 ) &
                      CALL lax_error__( " pztrtri ", " in MPI_Isend 17 ", ABS( ierr ) )
                   CALL MPI_Irecv(B, ldx*ldx, MPI_DOUBLE_COMPLEX, recv, recv+cicle, col_comm, req(2), ierr)
                   IF( ierr /= 0 ) &
                      CALL lax_error__( " pztrtri ", " in MPI_Irecv 18 ", ABS( ierr ) )
                   CALL zgemm('N', 'N', ldx, ldx, ldx, -ONE, C, ldx, BUF_RECV, ldx, first, sll, ldx)
                ELSE
                   CALL MPI_Isend(B, ldx*ldx, MPI_DOUBLE_COMPLEX, send, group_rank+cicle, col_comm, req(1), ierr)
                   IF( ierr /= 0 ) &
                      CALL lax_error__( " pztrtri ", " in MPI_Isend 19 ", ABS( ierr ) )
                   CALL MPI_Irecv(BUF_RECV, ldx*ldx, MPI_DOUBLE_COMPLEX, recv, recv+cicle, col_comm, req(2), ierr)
                   IF( ierr /= 0 ) &
                      CALL lax_error__( " pztrtri ", " in MPI_Irecv 20 ", ABS( ierr ) )
                   CALL zgemm('N', 'N', ldx, ldx, ldx, -ONE, C, ldx, B, ldx, ONE, sll, ldx)
                END IF
                !
                CALL MPI_Wait(req(1), status, ierr)
                IF( ierr /= 0 ) &
                   CALL lax_error__( " pztrtri ", " in MPI_Wait 21 ", ABS( ierr ) )
                !
                CALL MPI_Wait(req(2), status, ierr)
                IF( ierr /= 0 ) &
                   CALL lax_error__( " pztrtri ", " in MPI_Wait 22 ", ABS( ierr ) )
                !
             END IF
             cicle = cicle + 1
             first = ONE 
          END DO
       END IF
    END DO

    IF( myrow >= mycol ) THEN
       CALL mpi_comm_free( col_comm, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pztrtri ", " in mpi_comm_free 25 ", ABS( ierr ) )
    END IF

    DEALLOCATE(B)
    DEALLOCATE(C)
    DEALLOCATE(BUF_RECV)

#else

    CALL compute_ztrtri()

#endif

  CONTAINS

     SUBROUTINE compute_ztrtri()
       !
       !  clear the upper triangle (excluding diagonal terms) and
       !
       DO j = 1, ldx
          DO i = 1, j-1
             sll ( i, j ) = zero
          END DO
       END DO
       !
       CALL ztrtri( 'L', 'N', nr, sll, ldx, info )
       !
       IF( info /= 0 ) THEN
          CALL lax_error__( ' pztrtri ', ' problem in the local inversion ', info )
       END IF
       !
     END SUBROUTINE compute_ztrtri


     INTEGER FUNCTION shift ( idref, id, pos, size, dir )

       IMPLICIT NONE
   
       INTEGER :: idref, id, pos, size
       CHARACTER ( LEN = 1 ) :: dir
   
       IF( ( dir == 'E' ) .OR. ( dir == 'S' ) ) THEN
          shift = idref + MOD ( ( id - idref ) + pos, size )
       ELSE IF( ( dir == 'W' ) .OR. ( dir == 'N' ) ) THEN
          shift = idref + MOD ( ( id - idref ) - pos + size, size )
       ELSE
          shift = -1
       END IF
   
       RETURN

     END FUNCTION shift

END SUBROUTINE qe_pztrtri

!  now the Double Precision subroutine

SUBROUTINE qe_pdtrtri ( sll, ldx, n, desc )
    
    ! pztrtri computes the parallel inversion of a lower triangular matrix 
    ! distribuited among the processes using a 2-D block partitioning. 
    ! The algorithm is based on the schema below and executes the model 
    ! recursively to each column C2 under the diagonal.     
    !
    !     |-------|-------|      |--------------------|--------------------|
    !     |   A1  |   0   |      |   C1 = trtri(A1)   |          0         |
    ! A = |-------|-------|  C = |--------------------|--------------------|
    !     |   A2  |   A3  |      | C2 = -C3 * A2 * C1 |   C3 = trtri(A3)   | 
    !     |-------|-------|      |--------------------|--------------------|
    !
    ! The recursive steps of multiplication (C2 = -C3 * A2 * C1) is based on the Cannon's algorithms 
    ! for parallel matrix multiplication and is done with BLACS(dgemm)
    !
    !
    ! Arguments
    ! ============
    !
    ! sll   = local block of data
    ! ldx   = leading dimension of one block
    ! n     = size of the global array diributed among the blocks
    ! desc  = descriptor of the matrix distribution
    !
    !
    !  written by Ivan Girotto
    !

    USE descriptors

    IMPLICIT NONE

#if defined (__MPI)
   !
   include 'mpif.h'
   !
#endif
   
#include "la_param.f90"

    INTEGER, INTENT( IN ) :: n, ldx
    TYPE(la_descriptor), INTENT(IN) :: desc
    REAL(DP), INTENT( INOUT ) :: sll( ldx, ldx )

    REAL(DP), PARAMETER :: ONE = 1.0_DP
    REAL(DP), PARAMETER :: ZERO = 0.0_DP

#if defined __MPI
    INTEGER :: status(MPI_STATUS_SIZE)
#endif
    INTEGER :: req(2), ierr, col_comm
    INTEGER :: send, recv, group_rank, group_size
    INTEGER :: myrow, mycol, np, myid, comm

    ! counters
    INTEGER :: k, i, j, count, step_count, shiftcount, cicle 
    INTEGER :: C3dim   ! Dimension of submatrix B
    INTEGER :: nc, nr ! Local dimension of block
    INTEGER :: info, sup_recv
    INTEGER :: idrowref, idcolref, idref, idrecv 

    ! B and BUF_RECV are used to overload the computation of matrix multiplication and the shift of the blocks
    REAL(DP), ALLOCATABLE, DIMENSION( :, : ) :: B, C, BUF_RECV 
    REAL(DP) :: first

    myrow = desc%myr
    mycol = desc%myc
    myid  = desc%mype
    np    = desc%npr
    comm  = desc%comm

    IF( desc%npr /= desc%npc ) THEN
       CALL lax_error__( ' pdtrtri ', ' only square grid are allowed ', 1 ) 
    END IF
    IF( ldx /= desc%nrcx ) THEN
       CALL lax_error__( ' pdtrtri ', ' wrong leading dimension ldx ', ldx ) 
    END IF

    nr = desc%nr
    nc = desc%nc

    !  clear elements outside local meaningful block nr*nc

    DO j = nc+1, ldx
       DO i = 1, ldx
          sll( i, j ) = zero
       END DO
    END DO
    DO j = 1, ldx
       DO i = nr+1, ldx
          sll( i, j ) = zero
       END DO
    END DO

#if defined __MPI

    ALLOCATE( B( ldx, ldx ) )
    ALLOCATE( C( ldx, ldx ) )
    ALLOCATE( BUF_RECV ( ldx, ldx ) )

    IF( np == 2 ) THEN
       !
       !  special case with 4 proc, 2x2 grid
       !
       IF( myrow == mycol ) THEN
          CALL compute_dtrtri()
       END IF
       !
       CALL GRID2D_RANK( 'R', np, np, 1, 0, idref )
       !
       IF( myrow == 0 .AND. mycol == 0 ) THEN
          CALL MPI_Send(sll, ldx*ldx, MPI_DOUBLE_PRECISION, idref, 0, comm, ierr)
          IF( ierr /= 0 ) &
             CALL lax_error__( " pdtrtri ", " in MPI_Send 1 ", ABS( ierr ) )
       END IF
       !
       IF( myrow == 1 .AND. mycol == 1 ) THEN
          CALL MPI_Send(sll, ldx*ldx, MPI_DOUBLE_PRECISION, idref, 1, comm, ierr)
          IF( ierr /= 0 ) &
             CALL lax_error__( " pdtrtri ", " in MPI_Send 2 ", ABS( ierr ) )
       END IF
       !
       IF( myrow == 1 .AND. mycol == 0 ) THEN
          !
          CALL GRID2D_RANK( 'R', np, np, 0, 0, i )
          CALL GRID2D_RANK( 'R', np, np, 1, 1, j )
          !
          CALL MPI_Irecv( B, ldx*ldx, MPI_DOUBLE_PRECISION, i, 0, comm, req(1), ierr)
          IF( ierr /= 0 ) &
             CALL lax_error__( " pdtrtri ", " in MPI_Irecv 3 ", ABS( ierr ) )
          !
          CALL MPI_Irecv( C, ldx*ldx, MPI_DOUBLE_PRECISION, j, 1, comm, req(2), ierr)
          IF( ierr /= 0 ) &
             CALL lax_error__( " pdtrtri ", " in MPI_Irecv 4 ", ABS( ierr ) )
          !
          CALL MPI_Wait(req(1), status, ierr)
          IF( ierr /= 0 ) &
             CALL lax_error__( " pdtrtri ", " in MPI_Wait 5 ", ABS( ierr ) )
          !
          CALL dgemm('N', 'N', ldx, ldx, ldx, ONE, sll, ldx, b, ldx, ZERO, buf_recv, ldx)
          !
          CALL MPI_Wait(req(2), status, ierr)
          IF( ierr /= 0 ) &
             CALL lax_error__( " pdtrtri ", " in MPI_Wait 6 ", ABS( ierr ) )
          !
          CALL dgemm('N', 'N', ldx, ldx, ldx, -ONE, c, ldx, buf_recv, ldx, ZERO, sll, ldx)
          !
       END IF
       !
       IF( myrow == 0 .AND. mycol == 1 ) THEN
          !
          sll = zero
          !
       END IF
       !
       DEALLOCATE( b, c, buf_recv )
       !
       RETURN
       !
    END IF

  
    IF( myrow >= mycol ) THEN
       !
       !  only procs on lower triangle partecipates
       !
       CALL MPI_Comm_split( comm, mycol, myrow, col_comm, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pdtrtri ", " in MPI_Comm_split 9 ", ABS( ierr ) )

       CALL MPI_Comm_size( col_comm, group_size, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pdtrtri ", " in MPI_Comm_size 10 ", ABS( ierr ) )

       CALL MPI_Comm_rank( col_comm, group_rank, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pdtrtri ", " in MPI_Comm_rank 11 ", ABS( ierr ) )
       !
    ELSE
       !
       !  other procs stay at the window!
       !
       CALL MPI_Comm_split( comm, MPI_UNDEFINED, MPI_UNDEFINED, col_comm, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pdtrtri ", " in MPI_Comm_split 12 ", ABS( ierr ) )
       !
       sll = zero
       !
    END IF
    ! 

    ! Compute the inverse of a lower triangular 
    ! along the diagonal of the global array with BLAS(ztrtri) 
    !
    IF( mycol == myrow ) THEN
       !
       CALL compute_dtrtri()
       !
    ELSE IF( myrow > mycol ) THEN
       !
       buf_recv = sll
       !
    END IF

    IF( myrow >= mycol ) THEN
       !
       ! Broadcast the diagonal blocks to the processors under the diagonal
       !
       CALL MPI_Bcast( sll, ldx*ldx, MPI_DOUBLE_PRECISION, 0, col_comm, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pdtrtri ", " in MPI_Bcast 13 ", ABS( ierr ) )
       !
    END IF

    ! Compute A2 * C1 and start the Cannon's algorithm shifting the blocks of column one place to the North
    !
    IF( myrow > mycol ) THEN
       !
       CALL dgemm( 'N', 'N', ldx, ldx, ldx, ONE, buf_recv, ldx, sll, ldx, ZERO, c, ldx )
       !
       send = shift( 1, group_rank, 1, ( group_size - 1 ), 'N' )
       recv = shift( 1, group_rank, 1, ( group_size - 1 ), 'S' )
       !
       CALL MPI_Sendrecv( c, ldx*ldx, MPI_DOUBLE_PRECISION, send, 0, buf_recv, &
            ldx*ldx, MPI_DOUBLE_PRECISION, recv, 0, col_comm, status, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pdtrtri ", " in MPI_Sendrecv 14 ", ABS( ierr ) )
       !
    END IF

    ! Execute the Cannon's algorithm to compute ricorsively the multiplication of C2 = -C3 * A2 * C1
    !
    DO count = ( np - 2 ), 0, -1
       C3dim = (np-1) - count ! Dimension of the submatrix C3
       first = ZERO
       cicle = 0
       IF( ( myrow > count ) .AND. ( mycol >= count ) ) THEN
          idcolref = count + 1
          idrowref = myrow
          CALL GRID2D_RANK( 'R', np, np, idrowref, idcolref, idref )
          idrecv = idref - 1
          ! Compute C2 = -C3 * A2 * C1
          DO shiftcount = count, np-2
             IF(mycol>count)THEN
                ! Execute the virtual shift of the matrix C3 along the row in order to know which processor 
                ! have to send the block to C2 
                IF( cicle == 0)THEN
                   ! virtual shift of the block i,j of the submatrix C3 i place to West
                   send = shift(idref, myid, myrow-count, C3dim, 'W')
                ELSE
                   ! virtual shift of the block i,j of the submatrix C3 i place to West
                   send = shift(idref, send, 1, C3dim, 'E')
                END IF
                IF(send==idref)THEN
                   CALL MPI_Send(sll, ldx*ldx, MPI_DOUBLE_PRECISION, idrecv, myid, comm, ierr)
                   IF( ierr /= 0 ) &
                      CALL lax_error__( " pdtrtri ", " in MPI_Send 15 ", ABS( ierr ) )
                END IF
             ELSE
                IF( cicle == 0)THEN
                   ! virtual shift of the block i,j of the submatrix C3 i place to West
                   sup_recv = shift(idref, myid+1, myrow-count, C3dim, 'E')
                ELSE
                   ! virtual shift of the block i,j of the submatrix C3 i place to West
                   sup_recv = shift(idref, sup_recv, 1, C3dim, 'W')
                END IF
                CALL MPI_Recv(C, ldx*ldx, MPI_DOUBLE_PRECISION, sup_recv, sup_recv, comm, status, ierr)
                IF( ierr /= 0 ) &
                   CALL lax_error__( " pdtrtri ", " in MPI_Recv 16 ", ABS( ierr ) )
                send = shift(1, group_rank, 1, (group_size-1), 'S')
                recv = shift(1, group_rank, 1, (group_size-1), 'N')
                ! with the no-blocking communication the computation and the shift of the column block are overapped  
                IF( MOD( cicle, 2 ) == 0 ) THEN
                   !
                   CALL MPI_Isend(BUF_RECV, ldx*ldx, MPI_DOUBLE_PRECISION, send, group_rank+cicle, col_comm, req(1), ierr)
                   IF( ierr /= 0 ) &
                      CALL lax_error__( " pdtrtri ", " in MPI_Isend 17 ", ABS( ierr ) )
                   CALL MPI_Irecv(B, ldx*ldx, MPI_DOUBLE_PRECISION, recv, recv+cicle, col_comm, req(2), ierr)
                   IF( ierr /= 0 ) &
                      CALL lax_error__( " pdtrtri ", " in MPI_Irecv 18 ", ABS( ierr ) )
                   !
                   CALL dgemm('N', 'N', ldx, ldx, ldx, -ONE, C, ldx, BUF_RECV, ldx, first, sll, ldx)
                   !
                ELSE
                   !
                   CALL MPI_Isend(B, ldx*ldx, MPI_DOUBLE_PRECISION, send, group_rank+cicle, col_comm, req(1), ierr)
                   IF( ierr /= 0 ) &
                      CALL lax_error__( " pdtrtri ", " in MPI_Isend 19 ", ABS( ierr ) )
                   CALL MPI_Irecv(BUF_RECV, ldx*ldx, MPI_DOUBLE_PRECISION, recv, recv+cicle, col_comm, req(2), ierr)
                   IF( ierr /= 0 ) &
                      CALL lax_error__( " pdtrtri ", " in MPI_Irecv 20 ", ABS( ierr ) )
                   !
                   CALL dgemm('N', 'N', ldx, ldx, ldx, -ONE, C, ldx, B, ldx, ONE, sll, ldx)
                   !
                END IF
                !
                CALL MPI_Wait(req(1), status, ierr)
                IF( ierr /= 0 ) &
                   CALL lax_error__( " pdtrtri ", " in MPI_Wait 21 ", ABS( ierr ) )
                CALL MPI_Wait(req(2), status, ierr)
                IF( ierr /= 0 ) &
                   CALL lax_error__( " pdtrtri ", " in MPI_Wait 22 ", ABS( ierr ) )
                !
             END IF
             cicle = cicle + 1
             first = ONE 
          END DO
       END IF
    END DO

    IF( myrow >= mycol ) THEN
       CALL mpi_comm_free( col_comm, ierr )
       IF( ierr /= 0 ) &
          CALL lax_error__( " pdtrtri ", " in mpi_comm_free 25 ", ABS( ierr ) )
    END IF

    DEALLOCATE(B)
    DEALLOCATE(C)
    DEALLOCATE(BUF_RECV)

#else

    CALL compute_dtrtri()

#endif

  CONTAINS

     SUBROUTINE compute_dtrtri()
       !
       !  clear the upper triangle (excluding diagonal terms) and
       !
       DO j = 1, ldx
          DO i = 1, j-1
             sll ( i, j ) = zero
          END DO
       END DO
       !
       CALL dtrtri( 'L', 'N', nr, sll, ldx, info )
       !
       IF( info /= 0 ) THEN
          CALL lax_error__( ' pdtrtri ', ' problem in the local inversion ', info )
       END IF
       !
     END SUBROUTINE compute_dtrtri


     INTEGER FUNCTION shift ( idref, id, pos, size, dir )

       IMPLICIT NONE
   
       INTEGER :: idref, id, pos, size
       CHARACTER ( LEN = 1 ) :: dir
   
       IF( ( dir == 'E' ) .OR. ( dir == 'S' ) ) THEN
          shift = idref + MOD ( ( id - idref ) + pos, size )
       ELSE IF( ( dir == 'W' ) .OR. ( dir == 'N' ) ) THEN
          shift = idref + MOD ( ( id - idref ) - pos + size, size )
       ELSE
          shift = -1
       END IF
   
       RETURN

     END FUNCTION shift

END SUBROUTINE qe_pdtrtri



SUBROUTINE qe_pdsyevd( tv, n, desc, hh, ldh, e )
   USE descriptors
   USE dspev_module,     ONLY : pdspev_drv
   IMPLICIT NONE
#include "la_param.f90"
   LOGICAL, INTENT(IN) :: tv
       ! if tv is true compute eigenvalues and eigenvectors (not used)
   INTEGER, INTENT(IN) :: n, ldh
       ! n = matrix size, ldh = leading dimension of hh
   TYPE(la_descriptor), INTENT(IN) :: desc
       ! desc = descrittore della matrice 
   REAL(DP) :: hh( ldh, ldh )
       ! input:  hh = matrix to be diagonalized
   REAL(DP) :: e( n )
       ! output: hh = eigenvectors, e = eigenvalues

   INTEGER :: nrlx, nrl
   REAL(DP), ALLOCATABLE :: diag(:,:), vv(:,:)
   CHARACTER :: jobv

   nrl  = desc%nrl
   nrlx = desc%nrlx

   ALLOCATE( diag( nrlx, n ) )
   ALLOCATE( vv( nrlx, n ) )

   jobv = 'N'
   IF( tv ) jobv = 'V'
   !
   !  Redistribute matrix "hh" into "diag",  
   !  matrix "hh" is block distributed, matrix diag is cyclic distributed
   CALL blk2cyc_redist( n, diag, nrlx, n, hh, ldh, ldh, desc )
   !
   CALL pdspev_drv( jobv, diag, nrlx, e, vv, nrlx, nrl, n, &
        desc%npc * desc%npr, desc%mype, desc%comm )
   !
   IF( tv ) CALL cyc2blk_redist( n, vv, nrlx, n, hh, ldh, ldh, desc )
   !
   DEALLOCATE( vv )
   DEALLOCATE( diag )

   RETURN
END SUBROUTINE



SUBROUTINE qe_pzheevd( tv, n, desc, hh, ldh, e )
   USE descriptors
   USE zhpev_module,     ONLY : pzhpev_drv
   IMPLICIT NONE
#include "la_param.f90"
   LOGICAL, INTENT(IN) :: tv
       ! if tv is true compute eigenvalues and eigenvectors (not used)
   INTEGER, INTENT(IN) :: n, ldh
       ! n = matrix size, ldh = leading dimension of hh
   TYPE(la_descriptor), INTENT(IN) :: desc
       ! desc = descrittore della matrice 
   COMPLEX(DP) :: hh( ldh, ldh )
       ! input:  hh = matrix to be diagonalized
   REAL(DP) :: e( n )
       ! output: hh = eigenvectors, e = eigenvalues

   INTEGER :: nrlx, nrl
   COMPLEX(DP), ALLOCATABLE :: diag(:,:), vv(:,:)
   CHARACTER :: jobv

   nrl  = desc%nrl
   nrlx = desc%nrlx
   !
   ALLOCATE( diag( nrlx, n ) )
   ALLOCATE( vv( nrlx, n ) )
   !
   jobv = 'N'
   IF( tv ) jobv = 'V'

   CALL blk2cyc_zredist( n, diag, nrlx, n, hh, ldh, ldh, desc )
   !
   CALL pzhpev_drv( jobv, diag, nrlx, e, vv, nrlx, nrl, n, &
        desc%npc * desc%npr, desc%mype, desc%comm )
   !
   if( tv ) CALL cyc2blk_zredist( n, vv, nrlx, n, hh, ldh, ldh, desc )
   !
   DEALLOCATE( vv ) 
   DEALLOCATE( diag )

   RETURN
END SUBROUTINE



SUBROUTINE sqr_dsetmat( what, n, alpha, a, lda, desc )
   !
   !  Set the values of a square distributed matrix 
   !
   USE descriptors
   !
   IMPLICIT NONE

#include "la_param.f90"
   !
   CHARACTER(LEN=1), INTENT(IN) :: what
     ! what = 'A' set all the values of "a" equal to alpha
     ! what = 'U' set the values in the upper triangle of "a" equal to alpha
     ! what = 'L' set the values in the lower triangle of "a" equal to alpha
     ! what = 'D' set the values in the diagonal of "a" equal to alpha
   INTEGER, INTENT(IN) :: n
     ! dimension of the matrix
   REAL(DP), INTENT(IN) :: alpha
     ! value to be assigned to elements of "a"
   INTEGER, INTENT(IN) :: lda
     ! leading dimension of a
   REAL(DP) :: a(lda,*)
     ! matrix whose values have to be set
   TYPE(la_descriptor), INTENT(IN) :: desc
     ! descriptor of matrix a

   INTEGER :: i, j

   IF( desc%active_node < 0 ) THEN
      !
      !  processors not interested in this computation return quickly
      !
      RETURN
      !
   END IF

   SELECT CASE( what )
     CASE( 'U', 'u' )
        IF( desc%myc > desc%myr ) THEN
           DO j = 1, desc%nc
              DO i = 1, desc%nr
                 a( i, j ) = alpha
              END DO
           END DO
        ELSE IF( desc%myc == desc%myr ) THEN
           DO j = 1, desc%nc
              DO i = 1, j - 1
                 a( i, j ) = alpha
              END DO
           END DO
        END IF
     CASE( 'L', 'l' )
        IF( desc%myc < desc%myr ) THEN
           DO j = 1, desc%nc
              DO i = 1, desc%nr
                 a( i, j ) = alpha
              END DO
           END DO
        ELSE IF( desc%myc == desc%myr ) THEN
           DO j = 1, desc%nc
              DO i = j + 1, desc%nr
                 a( i, j ) = alpha
              END DO
           END DO
        END IF
     CASE( 'D', 'd' )
        IF( desc%myc == desc%myr ) THEN
           DO i = 1, desc%nr
              a( i, i ) = alpha
           END DO
        END IF
     CASE DEFAULT
        DO j = 1, desc%nc
           DO i = 1, desc%nr
              a( i, j ) = alpha
           END DO
        END DO
   END SELECT
   !
   RETURN
END SUBROUTINE sqr_dsetmat


SUBROUTINE sqr_zsetmat( what, n, alpha, a, lda, desc )
   !
   !  Set the values of a square distributed matrix 
   !
   USE descriptors
   !
   IMPLICIT NONE

#include "la_param.f90"
   !
   CHARACTER(LEN=1), INTENT(IN) :: what
     ! what = 'A' set all the values of "a" equal to alpha
     ! what = 'U' set the values in the upper triangle of "a" equal to alpha
     ! what = 'L' set the values in the lower triangle of "a" equal to alpha
     ! what = 'D' set the values in the diagonal of "a" equal to alpha
     ! what = 'H' clear the imaginary part of the diagonal of "a" 
   INTEGER, INTENT(IN) :: n
     ! dimension of the matrix
   COMPLEX(DP), INTENT(IN) :: alpha
     ! value to be assigned to elements of "a"
   INTEGER, INTENT(IN) :: lda
     ! leading dimension of a
   COMPLEX(DP) :: a(lda,*)
     ! matrix whose values have to be set
   TYPE(la_descriptor), INTENT(IN) :: desc
     ! descriptor of matrix a

   INTEGER :: i, j

   IF( desc%active_node < 0 ) THEN
      !
      !  processors not interested in this computation return quickly
      !
      RETURN
      !
   END IF

   SELECT CASE( what )
     CASE( 'U', 'u' )
        IF( desc%myc > desc%myr ) THEN
           DO j = 1, desc%nc
              DO i = 1, desc%nr
                 a( i, j ) = alpha
              END DO
           END DO
        ELSE IF( desc%myc == desc%myr ) THEN
           DO j = 1, desc%nc
              DO i = 1, j - 1
                 a( i, j ) = alpha
              END DO
           END DO
        END IF
     CASE( 'L', 'l' )
        IF( desc%myc < desc%myr ) THEN
           DO j = 1, desc%nc
              DO i = 1, desc%nr
                 a( i, j ) = alpha
              END DO
           END DO
        ELSE IF( desc%myc == desc%myr ) THEN
           DO j = 1, desc%nc
              DO i = j + 1, desc%nr
                 a( i, j ) = alpha
              END DO
           END DO
        END IF
     CASE( 'D', 'd' )
        IF( desc%myc == desc%myr ) THEN
           DO i = 1, desc%nr
              a( i, i ) = alpha
           END DO
        END IF
     CASE( 'H', 'h' )
        IF( desc%myc == desc%myr ) THEN
           DO i = 1, desc%nr
              a( i, i ) = CMPLX( DBLE( a(i,i) ), 0_DP, KIND=DP )
           END DO
        END IF
     CASE DEFAULT
        DO j = 1, desc%nc
           DO i = 1, desc%nr
              a( i, j ) = alpha
           END DO
        END DO
   END SELECT
   !
   RETURN
END SUBROUTINE sqr_zsetmat
