!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!==----------------------------------------------==!
    MODULE parallel_toolkit
!==----------------------------------------------==!

    USE kinds,     ONLY : DP
    USE io_global, ONLY : stdout
    USE parallel_include

    IMPLICIT NONE
    SAVE
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
  USE kinds
  USE descriptors
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: ldar
  REAL(DP)            :: ar(ldar,*)  !  matrix to be splitted, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  REAL(DP)            :: a(lda,*)
  INTEGER, INTENT(IN) :: desc( descla_siz_ )
  !
  REAL(DP), PARAMETER :: zero = 0_DP
  !
  INTEGER :: i, j, nr, nc, ic, ir, nx
  !
  IF( desc( lambda_node_ ) <= 0 ) THEN
     RETURN
  END IF

  nx  = desc( nlax_ )
  ir  = desc( ilar_ )
  ic  = desc( ilac_ )
  nr  = desc( nlar_ )
  nc  = desc( nlac_ )

  IF( lda < nx ) &
     CALL errore( " dsqmdst ", " inconsistent dimension lda ", lda )
  IF( n /= desc( la_n_ ) ) &
     CALL errore( " dsqmdst ", " inconsistent dimension n ", n )

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
  USE kinds
  USE descriptors
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: ldar
  COMPLEX(DP)         :: ar(ldar,*)  !  matrix to be splitted, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  COMPLEX(DP)         :: a(lda,*)
  INTEGER, INTENT(IN) :: desc( descla_siz_ )
  !
  COMPLEX(DP), PARAMETER :: zero = ( 0_DP , 0_DP )
  !
  INTEGER :: i, j, nr, nc, ic, ir, nx
  !
  IF( desc( lambda_node_ ) <= 0 ) THEN
     RETURN
  END IF

  nx  = desc( nlax_ )
  ir  = desc( ilar_ )
  ic  = desc( ilac_ )
  nr  = desc( nlar_ )
  nc  = desc( nlac_ )

  IF( lda < nx ) &
     CALL errore( " zsqmdst ", " inconsistent dimension lda ", lda )
  IF( n /= desc( la_n_ ) ) &
     CALL errore( " zsqmdst ", " inconsistent dimension n ", n )

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
  USE kinds
  USE descriptors
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: ldar
  REAL(DP)            :: ar(ldar,*)  !  matrix to be merged, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  REAL(DP)            :: a(lda,*)
  INTEGER, INTENT(IN) :: desc( descla_siz_ )
  INTEGER, INTENT(IN) :: comm
  !
  INTEGER :: i, j

#if defined __MPI
  !
  INTEGER :: np, nx, ipc, ipr, npr, npc, noff
  INTEGER :: ierr, ir, ic, nr, nc

  REAL(DP), ALLOCATABLE :: buf(:,:)
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     np = desc( la_npr_ ) * desc( la_npc_ ) 
     nx = desc( nlax_ )
     npr = desc( la_npr_ )
     npc = desc( la_npc_ )
     !
     IF( desc( la_myr_ ) == 0 .AND. desc( la_myc_ ) == 0 ) THEN
        ALLOCATE( buf( nx, nx * np ) )
     ELSE
        ALLOCATE( buf( 1, 1 ) )
     END IF
     !
     IF( lda /= nx ) &
        CALL errore( " dsqmcll ", " inconsistent dimension lda ", lda )
     !
     IF( desc( la_n_ ) /= n ) &
        CALL errore( " dsqmcll ", " inconsistent dimension n ", n )
     !
     CALL mpi_gather( a, nx*nx, mpi_double_precision, &
                      buf, nx*nx, mpi_double_precision, 0, desc( la_comm_ ) , ierr )
     !
     IF( desc( la_myr_ ) == 0 .AND. desc( la_myc_ ) == 0 ) THEN
        DO ipc = 1, npc
           CALL descla_local_dims( ic, nc, n, desc( la_nx_ ), npc, ipc-1 )
           DO ipr = 1, npr
              CALL descla_local_dims( ir, nr, n, desc( la_nx_ ), npr, ipr-1 )
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
  USE kinds
  USE descriptors
  !
  implicit none
  !
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: ldar
  COMPLEX(DP)         :: ar(ldar,*)  !  matrix to be merged, replicated on all proc
  INTEGER, INTENT(IN) :: lda
  COMPLEX(DP)         :: a(lda,*)
  INTEGER, INTENT(IN) :: desc( descla_siz_ )
  INTEGER, INTENT(IN) :: comm
  !
  INTEGER :: i, j

#if defined __MPI
  !
  INTEGER :: np, nx, ipc, ipr, npr, npc, noff
  INTEGER :: ierr, ir, ic, nr, nc

  COMPLEX(DP), ALLOCATABLE :: buf(:,:)
  !
  IF( desc( lambda_node_ ) > 0 ) THEN
     !
     np = desc( la_npr_ ) * desc( la_npc_ ) 
     nx = desc( nlax_ )
     npr = desc( la_npr_ )
     npc = desc( la_npc_ )
     !
     IF( desc( la_myr_ ) == 0 .AND. desc( la_myc_ ) == 0 ) THEN
        ALLOCATE( buf( nx, nx * np ) )
     ELSE
        ALLOCATE( buf( 1, 1 ) )
     END IF
     !
     IF( lda /= nx ) &
        CALL errore( " zsqmcll ", " inconsistent dimension lda ", lda )
     !
     IF( desc( la_n_ ) /= n ) &
        CALL errore( " zsqmcll ", " inconsistent dimension n ", n )
     !
     CALL mpi_gather( a, nx*nx, mpi_double_complex, &
                      buf, nx*nx, mpi_double_complex, 0, desc( la_comm_ ) , ierr )
     !
     IF( desc( la_myr_ ) == 0 .AND. desc( la_myc_ ) == 0 ) THEN
        DO ipc = 1, npc
           CALL descla_local_dims( ic, nc, n, desc( la_nx_ ), npc, ipc-1 )
           DO ipr = 1, npr
              CALL descla_local_dims( ir, nr, n, desc( la_nx_ ), npr, ipr-1 )
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
   USE kinds
   USE descriptors
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda
   REAL(DP)            :: a(lda,*)  !  matrix to be redistributed into b
   INTEGER, INTENT(IN) :: desc( descla_siz_ )
   !
   INTEGER :: i, j
   !
   DO j = 1, desc( nlac_ )
      DO i = desc( nlar_ ) + 1, desc( nlax_ )
         a( i, j ) = 0_DP
      END DO
   END DO
   DO j = desc( nlac_ ) + 1, desc( nlax_ )
      DO i = 1, desc( nlax_ )
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
   USE kinds
   USE descriptors
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda
   REAL(DP)            :: a(lda,*) 
   INTEGER, INTENT(IN) :: desc( descla_siz_ )
#if defined __MPI
   INTEGER :: istatus( MPI_STATUS_SIZE )
#endif
   INTEGER :: i, j
   INTEGER :: comm 
   INTEGER :: nr, nc, dest, sreq, ierr, sour
   REAL(DP) :: atmp

#if defined __MPI

   IF( desc( lambda_node_ ) <= 0 ) THEN
      RETURN
   END IF

   IF( n /= desc( la_n_ ) ) &
      CALL errore( " dsqmsym ", " wrong global dim n ", n )
   IF( lda /= desc( nlax_ ) ) &
      CALL errore( " dsqmsym ", " wrong leading dim lda ", lda )

   comm = desc( la_comm_ )

   nr = desc( nlar_ ) 
   nc = desc( nlac_ ) 
   IF( desc( la_myc_ ) == desc( la_myr_ ) ) THEN
      !
      !  diagonal block, procs work locally
      !
      DO j = 1, nc
         DO i = j + 1, nr
            a(i,j) = a(j,i)
         END DO
      END DO
      !
   ELSE IF( desc( la_myc_ ) > desc( la_myr_ ) ) THEN
      !
      !  super diagonal block, procs send the block to sub diag.
      !
      CALL GRID2D_RANK( 'R', desc( la_npr_ ), desc( la_npc_ ), &
                             desc( la_myc_ ), desc( la_myr_ ), dest )
      CALL mpi_isend( a, lda*lda, MPI_DOUBLE_PRECISION, dest, 1, comm, sreq, ierr )
      !
   ELSE IF( desc( la_myc_ ) < desc( la_myr_ ) ) THEN
      !
      !  sub diagonal block, procs receive the block from super diag,
      !  then transpose locally
      !
      CALL GRID2D_RANK( 'R', desc( la_npr_ ), desc( la_npc_ ), &
                             desc( la_myc_ ), desc( la_myr_ ), sour )
      CALL mpi_recv( a, lda*lda, MPI_DOUBLE_PRECISION, sour, 1, comm, istatus, ierr )
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

   IF( desc( la_myc_ ) > desc( la_myr_ ) ) THEN
      CALL MPI_Wait( sreq, istatus, ierr )
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
   USE kinds
   USE descriptors
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda
   COMPLEX(DP)         :: a(lda,*) 
   INTEGER, INTENT(IN) :: desc( descla_siz_ )
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

   IF( desc( lambda_node_ ) <= 0 ) THEN
      RETURN
   END IF

   IF( n /= desc( la_n_ ) ) &
      CALL errore( " zsqmsym ", " wrong global dim n ", n )
   IF( lda /= desc( nlax_ ) ) &
      CALL errore( " zsqmsym ", " wrong leading dim lda ", lda )

   comm = desc( la_comm_ )

   nr = desc( nlar_ ) 
   nc = desc( nlac_ ) 
   IF( desc( la_myc_ ) == desc( la_myr_ ) ) THEN
      !
      !  diagonal block, procs work locally
      !
      DO j = 1, nc
         a(j,j) = CMPLX( REAL( a(j,j) ), 0_DP )
         DO i = j + 1, nr
            a(i,j) = CONJG( a(j,i) )
         END DO
      END DO
      !
   ELSE IF( desc( la_myc_ ) > desc( la_myr_ ) ) THEN
      !
      !  super diagonal block, procs send the block to sub diag.
      !
      CALL GRID2D_RANK( 'R', desc( la_npr_ ), desc( la_npc_ ), &
                             desc( la_myc_ ), desc( la_myr_ ), dest )
      CALL mpi_isend( a, lda*lda, MPI_DOUBLE_COMPLEX, dest, 1, comm, sreq, ierr )
      !
   ELSE IF( desc( la_myc_ ) < desc( la_myr_ ) ) THEN
      !
      !  sub diagonal block, procs receive the block from super diag,
      !  then transpose locally
      !
      CALL GRID2D_RANK( 'R', desc( la_npr_ ), desc( la_npc_ ), &
                             desc( la_myc_ ), desc( la_myr_ ), sour )
      CALL mpi_recv( a, lda*lda, MPI_DOUBLE_COMPLEX, sour, 1, comm, istatus, ierr )
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

   IF( desc( la_myc_ ) > desc( la_myr_ ) ) THEN
      CALL MPI_Wait( sreq, istatus, ierr )
   END IF

#if defined __PIPPO
   CALL MPI_Comm_rank( comm, myid, ierr )
   ALLOCATE( tst1( n, n ) )
   ALLOCATE( tst2( n, n ) )
   tst1 = 0.0d0
   tst2 = 0.0d0
   do j = 1, desc( nlac_ )
   do i = 1, desc( nlar_ )
      tst1( i + desc( ilar_ ) - 1, j + desc( ilac_ ) - 1 ) = a( i , j )
   end do
   end do
   CALL MPI_REDUCE( tst1, tst2, n*n, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, comm, ierr )
   IF( myid == 0 ) THEN
   DO j = 1, n
      !
      IF( tst2(j,j) /=  CMPLX( REAL( tst2(j,j) ), 0_DP ) ) WRITE( 4000, * ) j, tst2(j,j)
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
      a(j,j) = CMPLX( REAL( a(j,j) ), 0_DP )
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
   USE kinds
   USE descriptors
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: na
   INTEGER, INTENT(IN) :: lda
   REAL(DP)            :: a(lda,*)  !  matrix to be redistributed into b
   INTEGER, INTENT(IN) :: desca( descla_siz_ )
   INTEGER, INTENT(IN) :: nb
   INTEGER, INTENT(IN) :: ldb
   REAL(DP)            :: b(ldb,*)
   INTEGER, INTENT(IN) :: descb( descla_siz_ )

   INTEGER :: ipc, ipr, npc, npr
   INTEGER :: ipr_old, ir_old, nr_old, irx_old
   INTEGER :: ipc_old, ic_old, nc_old, icx_old
   INTEGER :: myrow, mycol, ierr, rank
   INTEGER :: col_comm, row_comm, comm, sreq
   INTEGER :: nr_new, ir_new, irx_new, ir, nr, nrtot, irb, ire
   INTEGER :: nc_new, ic_new, icx_new, ic, nc, nctot, icb, ice
   INTEGER :: ib, i, j, myid
   INTEGER :: nrsnd( desca( la_npr_ ) )
   INTEGER :: ncsnd( desca( la_npr_ ) )
   INTEGER :: displ( desca( la_npr_ ) )
   INTEGER :: irb_new( desca( la_npr_ ) )
   INTEGER :: ire_new( desca( la_npr_ ) )
   INTEGER :: icb_new( desca( la_npr_ ) )
   INTEGER :: ice_new( desca( la_npr_ ) )
   REAL(DP), ALLOCATABLE :: buf(:)
   REAL(DP), ALLOCATABLE :: ab(:,:)
   REAL(DP), ALLOCATABLE :: tst1(:,:)
   REAL(DP), ALLOCATABLE :: tst2(:,:)
#if defined __MPI
    INTEGER :: istatus( MPI_STATUS_SIZE )
#endif

   IF( desca( lambda_node_ ) <= 0 ) THEN
      RETURN
   END IF

   ! preliminary consistency checks

   IF( nb < na ) &
      CALL errore( " dsqmred ", " nb < na, this sub. work only with nb >= na ", nb )
   IF( nb /= descb( la_n_ ) ) &
      CALL errore( " dsqmred ", " wrong global dim nb ", nb )
   IF( na /= desca( la_n_ ) ) &
      CALL errore( " dsqmred ", " wrong global dim na ", na )
   IF( ldb /= descb( nlax_ ) ) &
      CALL errore( " dsqmred ", " wrong leading dim ldb ", ldb )
   IF( lda /= desca( nlax_ ) ) &
      CALL errore( " dsqmred ", " wrong leading dim lda ", lda )

   npr   = desca( la_npr_ )
   myrow = desca( la_myr_ )
   npc   = desca( la_npc_ )
   mycol = desca( la_myc_ )
   comm  = desca( la_comm_ )

#if defined __MPI

   ! split communicator into row and col communicators

   CALL MPI_Comm_rank( comm, myid, ierr )

   CALL MPI_Comm_split( comm, mycol, myrow, col_comm, ierr )
   CALL MPI_Comm_split( comm, myrow, mycol, row_comm, ierr )

   CALL MPI_Comm_rank( col_comm, rank, ierr )
   IF( rank /= myrow ) &
      CALL errore( " dsqmred ", " building col_comm ", rank )

   CALL MPI_Comm_rank( row_comm, rank, ierr )
   IF( rank /= mycol ) &
      CALL errore( " dsqmred ", " building row_comm ", rank )

   ALLOCATE( buf( descb( nlax_ ) * descb( nlax_ ) ) )
   ALLOCATE( ab( descb( nlax_ ), desca( nlax_ ) ) )

   ! write( 3000 + myid, * ) 'na, nb = ', na, nb

   DO j = 1, descb( nlac_ )
      DO i = 1, descb( nlar_ )
         b( i, j ) = 0.0d0
      END DO
   END DO

   ab = 0.0d0

   ! first redistribute rows, column groups work in parallel

   DO ipr = 1, npr
      !
      CALL descla_local_dims( ir_new, nr_new, nb, descb( la_nx_ ), npr, ipr-1 )
      !
      irx_new = ir_new + nr_new - 1

      ! write( 3000 + myid, * ) 'ir_new, nr_new, irx_new = ', ir_new, nr_new, irx_new
      !
      DO ipr_old = 1, npr
         !
         CALL descla_local_dims( ir_old, nr_old, na, desca( la_nx_ ), npr, ipr_old-1 )
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
               DO j = 1, desca( nlac_ )
                  DO i = irb, ire
                     ib = ib + 1
                     buf( ib ) = a( i, j )
                  END DO
               END DO
               CALL mpi_isend( buf, ib, MPI_DOUBLE_PRECISION, ipr-1, ipr, col_comm, sreq, ierr )
            ELSE
               DO j = 1, desca( nlac_ )
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
            CALL errore( " dsqmred ", " somthing wrong with row 1 ", nrsnd( ipr_old ) )
         IF( nrsnd( ipr_old ) /= ire_new( ipr_old ) - irb_new( ipr_old ) + 1 ) &
            CALL errore( " dsqmred ", " somthing wrong with row 2 ", nrsnd( ipr_old ) )
         !
         nrsnd( ipr_old ) = nrsnd( ipr_old ) * desca( nlac_ )
         !
      END DO
      !
      IF( myrow == ipr - 1 ) THEN
         DO ipr_old = 1, npr
            IF( nrsnd( ipr_old ) > 0 ) THEN
               IF(  myrow /= ipr_old - 1 ) THEN
                  CALL mpi_recv( buf, nrsnd(ipr_old), MPI_DOUBLE_PRECISION, ipr_old-1, ipr, col_comm, istatus, ierr )
                  CALL MPI_GET_COUNT( istatus, MPI_DOUBLE_PRECISION, ib, ierr) 
                  IF( ib /= nrsnd(ipr_old) ) &
                     CALL errore( " dsqmred ", " somthing wrong with row 3 ", ib )
                  ib = 0
                  DO j = 1, desca( nlac_ )
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
            END IF
         END DO
      END IF
      !
   END DO

   ! then redistribute cols, row groups work in parallel

   DO ipc = 1, npc
      !
      CALL descla_local_dims( ic_new, nc_new, nb, descb( la_nx_ ), npc, ipc-1 )
      !
      icx_new = ic_new + nc_new - 1
      !
      ! write( 3000 + myid, * ) 'ic_new, nc_new, icx_new = ', ic_new, nc_new, icx_new
      !
      DO ipc_old = 1, npc
         !
         CALL descla_local_dims( ic_old, nc_old, na, desca( la_nx_ ), npc, ipc_old-1 )
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
                  DO i = 1, descb( nlax_ )
                     ib = ib + 1
                     buf( ib ) = ab( i, j )
                  END DO
               END DO
               CALL mpi_isend( buf, ib, MPI_DOUBLE_PRECISION, ipc-1, ipc, row_comm, sreq, ierr )
            ELSE
               ib = icb
               DO j = icb_new( ipc_old ), ice_new( ipc_old )
                  DO i = 1, descb( nlax_ )
                        b( i, j ) = ab( i, ib )
                  END DO
                  ib = ib + 1
               END DO
            END IF
         END IF

         IF( ncsnd( ipc_old ) /= ice-icb+1 ) &
            CALL errore( " dsqmred ", " somthing wrong with col 1 ", ncsnd( ipc_old ) )
         IF( ncsnd( ipc_old ) /= ice_new( ipc_old ) - icb_new( ipc_old ) + 1 ) &
            CALL errore( " dsqmred ", " somthing wrong with col 2 ", ncsnd( ipc_old ) )
         !
         ncsnd( ipc_old ) = ncsnd( ipc_old ) * descb( nlax_ )
         !
      END DO
      !
      IF( mycol == ipc - 1 ) THEN
         DO ipc_old = 1, npc
            IF( ncsnd( ipc_old ) > 0 ) THEN
               IF(  mycol /= ipc_old - 1 ) THEN
                  ib = icb_new( ipc_old )
                  CALL mpi_recv( b( 1, ib ), ncsnd(ipc_old), MPI_DOUBLE_PRECISION, ipc_old-1, ipc, row_comm, istatus, ierr )
                  CALL MPI_GET_COUNT( istatus, MPI_DOUBLE_PRECISION, ib, ierr ) 
                  IF( ib /= ncsnd(ipc_old) ) &
                     CALL errore( " dsqmred ", " somthing wrong with col 3 ", ib )
               END IF
            END IF
         END DO
      ELSE
         DO ipc_old = 1, npc
            IF( mycol == ipc_old - 1 .AND. ncsnd( ipc_old ) > 0 ) THEN
                CALL MPI_Wait( sreq, istatus, ierr )
            END IF
         END DO
      END IF
      !
   END DO

   DEALLOCATE( ab )
   DEALLOCATE( buf )

   CALL mpi_comm_free( col_comm, ierr )
   CALL mpi_comm_free( row_comm, ierr )

#if defined __PIPPO

   !  this is for debugging, tests through global matrix, if
   !  the two matrix (pre and before the redistribution) coincide.

   ALLOCATE( tst1( nb, nb ) )
   ALLOCATE( tst2( nb, nb ) )
   ALLOCATE( ab( nb, nb ) )

   ab = 0.0d0

   do j = 1, desca( nlac_ )
   do i = 1, desca( nlar_ )
      ab( i + desca( ilar_ ) - 1, j + desca( ilac_ ) - 1 ) = a( i , j )
   end do
   end do

   CALL MPI_REDUCE( ab, tst1, nb*nb, MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr )

   ab = 0.0d0

   do j = 1, descb( nlac_ )
   do i = 1, descb( nlar_ )
      ab( i + descb( ilar_ ) - 1, j + descb( ilac_ ) - 1 ) = b( i , j )
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
   USE kinds
   USE descriptors
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: na
   INTEGER, INTENT(IN) :: lda
   COMPLEX(DP)         :: a(lda,*)  !  matrix to be redistributed into b
   INTEGER, INTENT(IN) :: desca( descla_siz_ )
   INTEGER, INTENT(IN) :: nb
   INTEGER, INTENT(IN) :: ldb
   COMPLEX(DP)         :: b(ldb,*)
   INTEGER, INTENT(IN) :: descb( descla_siz_ )

   INTEGER :: ipc, ipr, npc, npr
   INTEGER :: ipr_old, ir_old, nr_old, irx_old
   INTEGER :: ipc_old, ic_old, nc_old, icx_old
   INTEGER :: myrow, mycol, ierr, rank
   INTEGER :: col_comm, row_comm, comm, sreq
   INTEGER :: nr_new, ir_new, irx_new, ir, nr, nrtot, irb, ire
   INTEGER :: nc_new, ic_new, icx_new, ic, nc, nctot, icb, ice
   INTEGER :: ib, i, j, myid
   INTEGER :: nrsnd( desca( la_npr_ ) )
   INTEGER :: ncsnd( desca( la_npr_ ) )
   INTEGER :: displ( desca( la_npr_ ) )
   INTEGER :: irb_new( desca( la_npr_ ) )
   INTEGER :: ire_new( desca( la_npr_ ) )
   INTEGER :: icb_new( desca( la_npr_ ) )
   INTEGER :: ice_new( desca( la_npr_ ) )
   COMPLEX(DP), ALLOCATABLE :: buf(:)
   COMPLEX(DP), ALLOCATABLE :: ab(:,:)
   COMPLEX(DP), ALLOCATABLE :: tst1(:,:)
   COMPLEX(DP), ALLOCATABLE :: tst2(:,:)
#if defined __MPI
    INTEGER :: istatus( MPI_STATUS_SIZE )
#endif

   IF( desca( lambda_node_ ) <= 0 ) THEN
      RETURN
   END IF

   ! preliminary consistency checks

   IF( nb < na ) &
      CALL errore( " zsqmred ", " nb < na, this sub. work only with nb >= na ", nb )
   IF( nb /= descb( la_n_ ) ) &
      CALL errore( " zsqmred ", " wrong global dim nb ", nb )
   IF( na /= desca( la_n_ ) ) &
      CALL errore( " zsqmred ", " wrong global dim na ", na )
   IF( ldb /= descb( nlax_ ) ) &
      CALL errore( " zsqmred ", " wrong leading dim ldb ", ldb )
   IF( lda /= desca( nlax_ ) ) &
      CALL errore( " zsqmred ", " wrong leading dim lda ", lda )

   npr   = desca( la_npr_ )
   myrow = desca( la_myr_ )
   npc   = desca( la_npc_ )
   mycol = desca( la_myc_ )
   comm  = desca( la_comm_ )

#if defined __MPI

   ! split communicator into row and col communicators

   CALL MPI_Comm_rank( comm, myid, ierr )

   CALL MPI_Comm_split( comm, mycol, myrow, col_comm, ierr )
   CALL MPI_Comm_split( comm, myrow, mycol, row_comm, ierr )

   CALL MPI_Comm_rank( col_comm, rank, ierr )
   IF( rank /= myrow ) &
      CALL errore( " zsqmred ", " building col_comm ", rank )

   CALL MPI_Comm_rank( row_comm, rank, ierr )
   IF( rank /= mycol ) &
      CALL errore( " zsqmred ", " building row_comm ", rank )

   ALLOCATE( buf( descb( nlax_ ) * descb( nlax_ ) ) )
   ALLOCATE( ab( descb( nlax_ ), desca( nlax_ ) ) )

   DO j = 1, descb( nlac_ )
      DO i = 1, descb( nlar_ )
         b( i, j ) = ( 0_DP , 0_DP )
      END DO
   END DO

   ab = ( 0_DP , 0_DP )

   ! first redistribute rows, column groups work in parallel

   DO ipr = 1, npr
      !
      CALL descla_local_dims( ir_new, nr_new, nb, descb( la_nx_ ), npr, ipr-1 )
      !
      irx_new = ir_new + nr_new - 1
      !
      DO ipr_old = 1, npr
         !
         CALL descla_local_dims( ir_old, nr_old, na, desca( la_nx_ ), npr, ipr_old-1 )
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
               DO j = 1, desca( nlac_ )
                  DO i = irb, ire
                     ib = ib + 1
                     buf( ib ) = a( i, j )
                  END DO
               END DO
               CALL mpi_isend( buf, ib, MPI_DOUBLE_COMPLEX, ipr-1, ipr, col_comm, sreq, ierr )
            ELSE
               DO j = 1, desca( nlac_ )
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
            CALL errore( " zsqmred ", " somthing wrong with row 1 ", nrsnd( ipr_old ) )
         IF( nrsnd( ipr_old ) /= ire_new( ipr_old ) - irb_new( ipr_old ) + 1 ) &
            CALL errore( " zsqmred ", " somthing wrong with row 2 ", nrsnd( ipr_old ) )
         !
         nrsnd( ipr_old ) = nrsnd( ipr_old ) * desca( nlac_ )
         !
      END DO
      !
      IF( myrow == ipr - 1 ) THEN
         DO ipr_old = 1, npr
            IF( nrsnd( ipr_old ) > 0 ) THEN
               IF(  myrow /= ipr_old - 1 ) THEN
                  CALL mpi_recv( buf, nrsnd(ipr_old), MPI_DOUBLE_COMPLEX, ipr_old-1, ipr, col_comm, istatus, ierr )
                  CALL MPI_GET_COUNT( istatus, MPI_DOUBLE_COMPLEX, ib, ierr) 
                  IF( ib /= nrsnd(ipr_old) ) &
                     CALL errore( " zsqmred ", " somthing wrong with row 3 ", ib )
                  ib = 0
                  DO j = 1, desca( nlac_ )
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
            END IF
         END DO
      END IF
      !
   END DO

   ! then redistribute cols, row groups work in parallel

   DO ipc = 1, npc
      !
      CALL descla_local_dims( ic_new, nc_new, nb, descb( la_nx_ ), npc, ipc-1 )
      !
      icx_new = ic_new + nc_new - 1
      !
      DO ipc_old = 1, npc
         !
         CALL descla_local_dims( ic_old, nc_old, na, desca( la_nx_ ), npc, ipc_old-1 )
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
                  DO i = 1, descb( nlax_ )
                     ib = ib + 1
                     buf( ib ) = ab( i, j )
                  END DO
               END DO
               CALL mpi_isend( buf, ib, MPI_DOUBLE_COMPLEX, ipc-1, ipc, row_comm, sreq, ierr )
            ELSE
               ib = icb
               DO j = icb_new( ipc_old ), ice_new( ipc_old )
                  DO i = 1, descb( nlax_ )
                        b( i, j ) = ab( i, ib )
                  END DO
                  ib = ib + 1
               END DO
            END IF
         END IF

         IF( ncsnd( ipc_old ) /= ice-icb+1 ) &
            CALL errore( " zsqmred ", " somthing wrong with col 1 ", ncsnd( ipc_old ) )
         IF( ncsnd( ipc_old ) /= ice_new( ipc_old ) - icb_new( ipc_old ) + 1 ) &
            CALL errore( " zsqmred ", " somthing wrong with col 2 ", ncsnd( ipc_old ) )
         !
         ncsnd( ipc_old ) = ncsnd( ipc_old ) * descb( nlax_ )
         !
      END DO
      !
      IF( mycol == ipc - 1 ) THEN
         DO ipc_old = 1, npc
            IF( ncsnd( ipc_old ) > 0 ) THEN
               IF(  mycol /= ipc_old - 1 ) THEN
                  ib = icb_new( ipc_old )
                  CALL mpi_recv( b( 1, ib ), ncsnd(ipc_old), MPI_DOUBLE_COMPLEX, ipc_old-1, ipc, row_comm, istatus, ierr )
                  CALL MPI_GET_COUNT( istatus, MPI_DOUBLE_COMPLEX, ib, ierr ) 
                  IF( ib /= ncsnd(ipc_old) ) &
                     CALL errore( " zsqmred ", " somthing wrong with col 3 ", ib )
               END IF
            END IF
         END DO
      ELSE
         DO ipc_old = 1, npc
            IF( mycol == ipc_old - 1 .AND. ncsnd( ipc_old ) > 0 ) THEN
                CALL MPI_Wait( sreq, istatus, ierr )
            END IF
         END DO
      END IF
      !
   END DO

   DEALLOCATE( ab )
   DEALLOCATE( buf )

   CALL mpi_comm_free( col_comm, ierr )
   CALL mpi_comm_free( row_comm, ierr )

#if defined __PIPPO

   !  this is for debugging, tests through global matrix, if
   !  the two matrix (pre and before the redistribution) coincide.

   ALLOCATE( tst1( nb, nb ) )
   ALLOCATE( tst2( nb, nb ) )
   ALLOCATE( ab( nb, nb ) )

   ab = 0.0d0

   do j = 1, desca( nlac_ )
   do i = 1, desca( nlar_ )
      ab( i + desca( ilar_ ) - 1, j + desca( ilac_ ) - 1 ) = a( i , j )
   end do
   end do

   CALL MPI_REDUCE( ab, tst1, nb*nb, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, comm, ierr )

   ab = 0.0d0

   do j = 1, descb( nlac_ )
   do i = 1, descb( nlar_ )
      ab( i + descb( ilar_ ) - 1, j + descb( ilac_ ) - 1 ) = b( i , j )
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

  REAL(8), ALLOCATABLE :: auxa( : )
  REAL(8), ALLOCATABLE :: auxc( : )

  !
  ! ... BODY
  !

  CALL MPI_COMM_SIZE(comm, NPROC, IERR)
  CALL MPI_COMM_RANK(comm, MPIME, IERR)

  IF ( NPROC == 1 ) THEN

     !  if there is only one proc no need of using parallel alg.

     CALL DGEMM(TRANSA, TRANSB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

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

  ALLOCATE( auxa( MAX( n, m ) * ldx ) )
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

  CALL DGEMM( 'N', transb, nb, n, k, alpha, auxa(1), ldx, B, ldb, beta, auxc(1), ldx )

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

     CALL DGEMM(TRANSA, TRANSB, m, N, k, alpha, A, lda, B, ldb, beta, C, ldc)

#endif

  RETURN

END SUBROUTINE rep_matmul_drv


SUBROUTINE zrep_matmul_drv( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC, comm )
  !
  !  Parallel matrix multiplication with replicated matrix
  !  written by Carlo Cavazzoni
  !
  implicit none
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

     CALL ZGEMM(TRANSA, TRANSB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

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

  ALLOCATE( auxa( MAX( n, m ) * ldx ) )
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

  CALL ZGEMM( 'N', transb, nb, n, k, alpha, auxa(1), ldx, B, ldb, beta, auxc(1), ldx )

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

     CALL ZGEMM(TRANSA, TRANSB, m, N, k, alpha, A, lda, B, ldb, beta, C, ldc)

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
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , &
                           la_comm_ , lambda_node_ , la_npr_ , la_npc_ , la_myr_ , la_myc_
   !
   IMPLICIT NONE
   !
   CHARACTER(LEN=1), INTENT(IN) :: transa, transb
   INTEGER, INTENT(IN) :: n
   REAL(DP), INTENT(IN) :: alpha, beta
   INTEGER, INTENT(IN) :: lda, ldb, ldc
   REAL(DP) :: a(lda,*), b(ldb,*), c(ldc,*)
   INTEGER, INTENT(IN) :: desc(*)
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
   IF( desc( lambda_node_ ) < 0 ) THEN
      !
      !  processors not interested in this computation return quickly
      !
      RETURN
      !
   END IF

   IF( desc( la_npr_ ) == 1 ) THEN 
      !
      !  quick return if only one processor is used 
      !
      CALL dgemm( TRANSA, TRANSB, n, n, n, alpha, a, lda, b, ldb, beta, c, ldc)
      !
      RETURN
      !
   END IF

   IF( desc( la_npr_ ) /= desc( la_npc_ ) ) &
      CALL errore( ' sqr_mm_cannon ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' sqr_mm_cannon ', ' n less or equal zero ', 1 )
   IF( n < desc( la_npr_ ) ) &
      CALL errore( ' sqr_mm_cannon ', ' n less than the mesh size ', 1 )
   IF( n < desc( nlax_ ) ) &
      CALL errore( ' sqr_mm_cannon ', ' n less than the block size ', 1 )

   !
   !  Retrieve communicator and mesh geometry
   !
   np    = desc( la_npr_ )
   comm  = desc( la_comm_ )
   rowid = desc( la_myr_  )
   colid = desc( la_myc_  )
   !
   !  Retrieve the size of the local block
   !
   nr    = desc( nlar_ ) 
   nc    = desc( nlac_ ) 
   nb    = desc( nlax_ )
   !
#if defined (__MPI)
   CALL MPI_BARRIER( comm, ierr )
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
      INTEGER,          INTENT(IN) :: ln       ! shift lenght
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

         CALL errore( ' sqr_mm_cannon ', ' unknown shift direction ', 1 )

      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_PRECISION, &
           idest, tag, isour, tag, comm, istatus, ierr)
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

         CALL errore( ' sqr_mm_cannon ', ' unknown shift_exch direction ', 1 )

      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_PRECISION, &
           idest, tag, isour, tag, comm, istatus, ierr)
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
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , &
                           la_comm_ , lambda_node_ , la_npr_ , la_npc_ , la_myr_ , la_myc_
   !
   IMPLICIT NONE
   !
   CHARACTER(LEN=1), INTENT(IN) :: transa, transb
   INTEGER, INTENT(IN) :: n
   COMPLEX(DP), INTENT(IN) :: alpha, beta
   INTEGER, INTENT(IN) :: lda, ldb, ldc
   COMPLEX(DP) :: a(lda,*), b(ldb,*), c(ldc,*)
   INTEGER, INTENT(IN) :: desc(*)
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
   IF( desc( lambda_node_ ) < 0 ) THEN
      !
      !  processors not interested in this computation return quickly
      !
      RETURN
      !
   END IF

   IF( desc( la_npr_ ) == 1 ) THEN 
      !
      !  quick return if only one processor is used 
      !
      CALL zgemm( TRANSA, TRANSB, n, n, n, alpha, a, lda, b, ldb, beta, c, ldc)
      !
      RETURN
      !
   END IF

   IF( desc( la_npr_ ) /= desc( la_npc_ ) ) &
      CALL errore( ' sqr_zmm_cannon ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' sqr_zmm_cannon ', ' n less or equal zero ', 1 )
   IF( n < desc( la_npr_ ) ) &
      CALL errore( ' sqr_zmm_cannon ', ' n less than the mesh size ', 1 )
   IF( n < desc( nlax_ ) ) &
      CALL errore( ' sqr_zmm_cannon ', ' n less than the block size ', 1 )

   !
   !  Retrieve communicator and mesh geometry
   !
   np    = desc( la_npr_ )
   comm  = desc( la_comm_ )
   rowid = desc( la_myr_  )
   colid = desc( la_myc_  )
   !
   !  Retrieve the size of the local block
   !
   nr    = desc( nlar_ ) 
   nc    = desc( nlac_ ) 
   nb    = desc( nlax_ )
   !
#if defined (__MPI)
   CALL MPI_BARRIER( comm, ierr )
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
      INTEGER,          INTENT(IN) :: ln       ! shift lenght
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

         CALL errore( ' sqr_zmm_cannon ', ' unknown shift direction ', 1 )

      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_COMPLEX, &
           idest, tag, isour, tag, comm, istatus, ierr)
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

         CALL errore( ' sqr_zmm_cannon ', ' unknown shift_exch direction ', 1 )

      END IF
      !
      CALL GRID2D_RANK( 'R', np, np, irdst, icdst, idest )
      CALL GRID2D_RANK( 'R', np, np, irsrc, icsrc, isour )
      !
#if defined (__MPI)
      !
      CALL MPI_SENDRECV_REPLACE(blk, nb*nb, MPI_DOUBLE_COMPLEX, &
           idest, tag, isour, tag, comm, istatus, ierr)
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
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , la_npc_ , la_n_ , &
                           la_comm_ , lambda_node_ , la_npr_ , la_myr_ , la_myc_
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   REAL(DP)            :: a(lda,*), b(ldb,*)
   INTEGER, INTENT(IN) :: desc(*)
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
   IF( desc( lambda_node_ ) < 0 ) THEN
      RETURN
   END IF

   IF( desc( la_npr_ ) == 1 ) THEN
      CALL mytranspose( a, lda, b, ldb, n, n )
      RETURN
   END IF

   IF( desc( la_npr_ ) /= desc( la_npc_ ) ) &
      CALL errore( ' sqr_tr_cannon ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' sqr_tr_cannon ', ' n less or equal zero ', 1 )
   IF( n /= desc( la_n_ ) ) &
      CALL errore( ' sqr_tr_cannon ', ' inconsistent size n  ', 1 )
   IF( lda /= desc( nlax_ ) ) &
      CALL errore( ' sqr_tr_cannon ', ' inconsistent size lda  ', 1 )
   IF( ldb /= desc( nlax_ ) ) &
      CALL errore( ' sqr_tr_cannon ', ' inconsistent size ldb  ', 1 )

   comm = desc( la_comm_ )


   rowid = desc( la_myr_ )
   colid = desc( la_myc_ )
   np    = desc( la_npr_ )
   !
   !  Compute the size of the local block
   !
   nr = desc( nlar_ ) 
   nc = desc( nlac_ ) 
   nb = desc( nlax_ )
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
      !
#endif

      RETURN
   END SUBROUTINE


END SUBROUTINE

!
!
!

SUBROUTINE cyc2blk_redist( n, a, lda, b, ldb, desc )
   !
   !  Parallel square matrix redistribution.
   !  A (input) is cyclically distributed by rows across processors
   !  B (output) is distributed by block across 2D processors grid
   !
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , lambda_node_ , la_npr_ , &
                           descla_siz_ , la_npc_ , la_n_ , la_me_ , la_comm_
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   REAL(DP) :: a(lda,*), b(ldb,*)
   INTEGER  :: desc(*)
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
   !
   real(DP), allocatable :: rcvbuf(:,:,:)
   real(DP), allocatable :: sndbuf(:,:)
   integer, allocatable :: ip_desc(:,:)
   !
   character(len=256) :: msg
   !
#if defined (__MPI)

   IF( desc( lambda_node_ ) < 0 ) THEN
      RETURN
   END IF

   np = desc( la_npr_ )  !  dimension of the processor mesh
   nb = desc( nlax_ )    !  leading dimension of the local matrix block
   me = desc( la_me_ )   !  my processor id (starting from 0)
   comm_a = desc( la_comm_ )
   nproc  = desc( la_npr_ ) * desc( la_npc_ )

   IF( np /= desc( la_npc_ ) ) &
      CALL errore( ' cyc2blk_redist ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' cyc2blk_redist ', ' n less or equal zero ', 1 )
   IF( desc( la_n_ ) < nproc ) &
      CALL errore( ' cyc2blk_redist ', ' nb less than the number of proc ', 1 )

   ALLOCATE( ip_desc( descla_siz_ , nproc ) )

   CALL mpi_allgather( desc, descla_siz_ , mpi_integer, ip_desc, descla_siz_ , mpi_integer, comm_a, ierr )
   !
   nbuf = (nb/nproc+2) * nb
   !
   ALLOCATE( sndbuf( nb/nproc+2, nb ) )
   ALLOCATE( rcvbuf( nb/nproc+2, nb, nproc ) )
   
   DO ip = 0, nproc - 1
      !
      IF( ip_desc( lambda_node_ , ip + 1 ) > 0 ) THEN

         ip_nr = ip_desc( nlar_ , ip + 1) 
         ip_nc = ip_desc( nlac_ , ip + 1) 
         ip_ir = ip_desc( ilar_ , ip + 1) 
         ip_ic = ip_desc( ilac_ , ip + 1) 
         !
         DO j = 1, ip_nc
            jj = j + ip_ic - 1
            il = 1
            DO i = 1, ip_nr
               ii = i + ip_ir - 1
               IF( MOD( ii - 1, nproc ) == me ) THEN
                  sndbuf( il, j ) = a( ( ii - 1 )/nproc + 1, jj )
                  il = il + 1
               END IF 
            END DO
         END DO
   
      END IF

      CALL mpi_gather( sndbuf, nbuf, mpi_double_precision, &
                       rcvbuf, nbuf, mpi_double_precision, ip, comm_a, ierr )
     
   END DO

   !
   nr  = desc( nlar_ ) 
   nc  = desc( nlac_ )
   ir  = desc( ilar_ )
   ic  = desc( ilac_ )
   !
   DO ip = 0, nproc - 1
      DO j = 1, nc
         il = 1
         DO i = 1, nr
            ii = i + ir - 1
            IF( MOD( ii - 1, nproc ) == ip ) THEN
               b( i, j ) = rcvbuf( il, j, ip+1 )
               il = il + 1
            END IF 
         END DO
      END DO
   END DO
   !
   !
   DEALLOCATE( ip_desc )
   DEALLOCATE( rcvbuf )
   DEALLOCATE( sndbuf )
   
#else

   b( 1:n, 1:n ) = a( 1:n, 1:n )   

#endif

   RETURN

END SUBROUTINE cyc2blk_redist


SUBROUTINE cyc2blk_zredist( n, a, lda, b, ldb, desc )
   !
   !  Parallel square matrix redistribution.
   !  A (input) is cyclically distributed by rows across processors
   !  B (output) is distributed by block across 2D processors grid
   !
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , lambda_node_ , la_npr_ , &
                           descla_siz_ , la_npc_ , la_n_ , la_me_ , la_comm_
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   COMPLEX(DP) :: a(lda,*), b(ldb,*)
   INTEGER  :: desc(*)
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
   !
   COMPLEX(DP), allocatable :: rcvbuf(:,:,:)
   COMPLEX(DP), allocatable :: sndbuf(:,:)
   integer, allocatable :: ip_desc(:,:)
   !
   character(len=256) :: msg
   !
#if defined (__MPI)

   IF( desc( lambda_node_ ) < 0 ) THEN
      RETURN
   END IF

   np = desc( la_npr_ )  !  dimension of the processor mesh
   nb = desc( nlax_ )    !  leading dimension of the local matrix block
   me = desc( la_me_ )   !  my processor id (starting from 0)
   comm_a = desc( la_comm_ )
   nproc  = desc( la_npr_ ) * desc( la_npc_ )

   IF( np /= desc( la_npc_ ) ) &
      CALL errore( ' cyc2blk_zredist ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' cyc2blk_zredist ', ' n less or equal zero ', 1 )
   IF( desc( la_n_ ) < nproc ) &
      CALL errore( ' cyc2blk_zredist ', ' nb less than the number of proc ', 1 )

   ALLOCATE( ip_desc( descla_siz_ , nproc ) )

   CALL mpi_allgather( desc, descla_siz_ , mpi_integer, ip_desc, descla_siz_ , mpi_integer, comm_a, ierr )
   !
   nbuf = (nb/nproc+2) * nb
   !
   ALLOCATE( sndbuf( nb/nproc+2, nb ) )
   ALLOCATE( rcvbuf( nb/nproc+2, nb, nproc ) )
   
   DO ip = 0, nproc - 1
      !
      IF( ip_desc( lambda_node_ , ip + 1 ) > 0 ) THEN

         ip_nr = ip_desc( nlar_ , ip + 1) 
         ip_nc = ip_desc( nlac_ , ip + 1) 
         ip_ir = ip_desc( ilar_ , ip + 1) 
         ip_ic = ip_desc( ilac_ , ip + 1) 
         !
         DO j = 1, ip_nc
            jj = j + ip_ic - 1
            il = 1
            DO i = 1, ip_nr
               ii = i + ip_ir - 1
               IF( MOD( ii - 1, nproc ) == me ) THEN
                  sndbuf( il, j ) = a( ( ii - 1 )/nproc + 1, jj )
                  il = il + 1
               END IF 
            END DO
         END DO
   
      END IF

      CALL mpi_gather( sndbuf, nbuf, mpi_double_complex, &
                       rcvbuf, nbuf, mpi_double_complex, ip, comm_a, ierr )
     
   END DO

   !
   nr  = desc( nlar_ ) 
   nc  = desc( nlac_ )
   ir  = desc( ilar_ )
   ic  = desc( ilac_ )
   !
   DO ip = 0, nproc - 1
      DO j = 1, nc
         il = 1
         DO i = 1, nr
            ii = i + ir - 1
            IF( MOD( ii - 1, nproc ) == ip ) THEN
               b( i, j ) = rcvbuf( il, j, ip+1 )
               il = il + 1
            END IF 
         END DO
      END DO
   END DO
   !
   !
   DEALLOCATE( ip_desc )
   DEALLOCATE( rcvbuf )
   DEALLOCATE( sndbuf )
   
#else

   b( 1:n, 1:n ) = a( 1:n, 1:n )   

#endif

   RETURN

END SUBROUTINE cyc2blk_zredist


SUBROUTINE blk2cyc_redist( n, a, lda, b, ldb, desc )
   !
   !  Parallel square matrix redistribution.
   !  A (output) is cyclically distributed by rows across processors
   !  B (input) is distributed by block across 2D processors grid
   !
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , lambda_node_ , la_npr_ , &
                           descla_siz_ , la_npc_ , la_n_ , la_me_ , la_comm_
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   REAL(DP) :: a(lda,*), b(ldb,*)
   INTEGER  :: desc(*)
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
   !
   real(DP), allocatable :: rcvbuf(:,:,:)
   real(DP), allocatable :: sndbuf(:,:)
   integer, allocatable :: ip_desc(:,:)
   !
   character(len=256) :: msg
   !
#if defined (__MPI)

   IF( desc( lambda_node_ ) < 0 ) THEN
      RETURN
   END IF

   np = desc( la_npr_ )  !  dimension of the processor mesh
   nb = desc( nlax_ )    !  leading dimension of the local matrix block
   me = desc( la_me_ )   !  my processor id (starting from 0)
   comm_a = desc( la_comm_ )
   nproc  = desc( la_npr_ ) * desc( la_npc_ )

   IF( np /= desc( la_npc_ ) ) &
      CALL errore( ' blk2cyc_redist ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' blk2cyc_redist ', ' n less or equal zero ', 1 )
   IF( desc( la_n_ ) < nproc ) &
      CALL errore( ' blk2cyc_redist ', ' nb less than the number of proc ', 1 )

   ALLOCATE( ip_desc( descla_siz_ , nproc ) )

   CALL mpi_allgather( desc, descla_siz_ , mpi_integer, ip_desc, descla_siz_ , mpi_integer, comm_a, ierr )
   !
   nbuf = (nb/nproc+2) * nb
   !
   ALLOCATE( sndbuf( nb/nproc+2, nb ) )
   ALLOCATE( rcvbuf( nb/nproc+2, nb, nproc ) )
   !
   nr  = desc( nlar_ ) 
   nc  = desc( nlac_ )
   ir  = desc( ilar_ )
   ic  = desc( ilac_ )
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
      CALL mpi_gather( sndbuf, nbuf, mpi_double_precision, &
                       rcvbuf, nbuf, mpi_double_precision, ip, comm_a, ierr )
   END DO
   !
   
   DO ip = 0, nproc - 1
      !
      IF( ip_desc( lambda_node_ , ip + 1 ) > 0 ) THEN

         ip_nr = ip_desc( nlar_ , ip + 1) 
         ip_nc = ip_desc( nlac_ , ip + 1) 
         ip_ir = ip_desc( ilar_ , ip + 1) 
         ip_ic = ip_desc( ilac_ , ip + 1) 
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
   
      END IF
     
   END DO
   !
   DEALLOCATE( ip_desc )
   DEALLOCATE( rcvbuf )
   DEALLOCATE( sndbuf )
   
#else

   a( 1:n, 1:n ) = b( 1:n, 1:n )   

#endif

   RETURN

END SUBROUTINE blk2cyc_redist


SUBROUTINE blk2cyc_zredist( n, a, lda, b, ldb, desc )
   !
   !  Parallel square matrix redistribution.
   !  A (output) is cyclically distributed by rows across processors
   !  B (input) is distributed by block across 2D processors grid
   !
   USE kinds,       ONLY : DP
   USE descriptors, ONLY : ilar_ , nlar_ , ilac_ , nlac_ , nlax_ , lambda_node_ , la_npr_ , &
                           descla_siz_ , la_npc_ , la_n_ , la_me_ , la_comm_
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: n
   INTEGER, INTENT(IN) :: lda, ldb
   COMPLEX(DP) :: a(lda,*), b(ldb,*)
   INTEGER  :: desc(*)
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
   !
   COMPLEX(DP), allocatable :: rcvbuf(:,:,:)
   COMPLEX(DP), allocatable :: sndbuf(:,:)
   integer, allocatable :: ip_desc(:,:)
   !
   character(len=256) :: msg
   !
#if defined (__MPI)

   IF( desc( lambda_node_ ) < 0 ) THEN
      RETURN
   END IF

   np = desc( la_npr_ )  !  dimension of the processor mesh
   nb = desc( nlax_ )    !  leading dimension of the local matrix block
   me = desc( la_me_ )   !  my processor id (starting from 0)
   comm_a = desc( la_comm_ )
   nproc  = desc( la_npr_ ) * desc( la_npc_ )

   IF( np /= desc( la_npc_ ) ) &
      CALL errore( ' blk2cyc_zredist ', ' works only with square processor mesh ', 1 )
   IF( n < 1 ) &
      CALL errore( ' blk2cyc_zredist ', ' n less or equal zero ', 1 )
   IF( desc( la_n_ ) < nproc ) &
      CALL errore( ' blk2cyc_zredist ', ' nb less than the number of proc ', 1 )

   ALLOCATE( ip_desc( descla_siz_ , nproc ) )

   CALL mpi_allgather( desc, descla_siz_ , mpi_integer, ip_desc, descla_siz_ , mpi_integer, comm_a, ierr )
   !
   nbuf = (nb/nproc+2) * nb
   !
   ALLOCATE( sndbuf( nb/nproc+2, nb ) )
   ALLOCATE( rcvbuf( nb/nproc+2, nb, nproc ) )
   !
   nr  = desc( nlar_ ) 
   nc  = desc( nlac_ )
   ir  = desc( ilar_ )
   ic  = desc( ilac_ )
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
      CALL mpi_gather( sndbuf, nbuf, mpi_double_complex, &
                       rcvbuf, nbuf, mpi_double_complex, ip, comm_a, ierr )
   END DO
   !
   
   DO ip = 0, nproc - 1
      !
      IF( ip_desc( lambda_node_ , ip + 1 ) > 0 ) THEN

         ip_nr = ip_desc( nlar_ , ip + 1) 
         ip_nc = ip_desc( nlac_ , ip + 1) 
         ip_ir = ip_desc( ilar_ , ip + 1) 
         ip_ic = ip_desc( ilac_ , ip + 1) 
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
   
      END IF
     
   END DO
   !
   DEALLOCATE( ip_desc )
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

SUBROUTINE pzpotf( sll, ldx, n, desc )
   !
   use descriptors, ONLY: descla_local_dims, descla_siz_ , la_myr_ , la_myc_ , la_me_ , nlax_ ,&
                          nlar_ , nlac_ , ilar_ , ilac_ , la_comm_ , la_nx_ , la_npr_ , la_npc_
   use parallel_include
   use kinds
   !
   implicit none
   !
   integer :: n, ldx
   integer :: desc( descla_siz_ )
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

   myrow = desc( la_myr_ )
   mycol = desc( la_myc_ )
   myid  = desc( la_me_ )
   np    = desc( la_npr_ )

   IF( desc( la_npr_ ) /= desc( la_npc_ ) ) THEN
      CALL errore( ' pzpotf ', ' only square grid are allowed ', 1 ) 
   END IF

   IF( ldx /= desc( nlax_ ) ) THEN
      CALL errore( ' pzpotf ', ' wrong leading dimension ldx ', ldx ) 
   END IF

   nr = desc( nlar_ )
   nc = desc( nlac_ )

   ALLOCATE( ssnd( ldx, ldx ) )
   ALLOCATE( srcv( ldx, ldx ) )

   DO jb = 1, np
      !
      !    Update and factorize the current diagonal block and test
      !    for non-positive-definiteness.
      !
      CALL descla_local_dims( jir, jnr, n, desc( la_nx_ ), np, jb-1 )
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
      CALL mpi_comm_split( desc( la_comm_ ) , color, key, ccomm, ierr )
      !  
      IF( myrow >= jb-1 .and. mycol <= jb-1 ) THEN
          color = myrow
          key   = mycol
      ELSE
          color = np
          key   = myid
      END IF
      !
      CALL mpi_comm_split( desc( la_comm_ ), color, key, rcomm, ierr )
      !
      !    here every process can work independently, then we need a reduce.
      !
      IF( jb > 1 ) THEN
         !
         DO ib = 1, jb - 1
            IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( ib - 1 ) ) ) THEN
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
            CALL MPI_REDUCE( ssnd, sll, ldx*ldx, MPI_DOUBLE_COMPLEX, MPI_SUM, jb-1, rcomm, ierr )
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
            CALL errore( " pzpotrf ", " problems computing cholesky decomposition ", ABS( info ) )
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
            CALL mpi_bcast( sll,  ldx*ldx, MPI_DOUBLE_COMPLEX, 0, ccomm, ierr )   
         ELSE IF( ( myrow > ( jb - 1 ) ) .AND. ( mycol < ( jb - 1 ) ) ) THEN
            CALL mpi_bcast( srcv, ldx*ldx, MPI_DOUBLE_COMPLEX, 0, ccomm, ierr )   
         END IF
         !
         DO ib = jb + 1, np
            CALL descla_local_dims( iir, inr, n, desc( la_nx_ ), np, ib-1 )
            DO kb = 1, jb - 1
               CALL descla_local_dims( kic, knc, n, desc( la_nx_ ), np, kb-1 )
               IF( ( myrow == ( ib - 1 ) ) .AND. ( mycol == ( kb - 1 ) ) ) THEN
                  CALL ZGEMM( 'N', 'C', inr, jnr, knc, -CONE, sll, ldx, srcv, ldx, czero, ssnd, ldx )
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
                  CALL MPI_REDUCE( ssnd, sll, ldx*ldx, MPI_DOUBLE_COMPLEX, MPI_SUM, jb-1, rcomm, ierr )
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
            CALL mpi_bcast( sll,  ldx*ldx, MPI_DOUBLE_COMPLEX, 0, ccomm, ierr )   
         ELSE IF( ( myrow > ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
            CALL mpi_bcast( srcv,  ldx*ldx, MPI_DOUBLE_COMPLEX, 0, ccomm, ierr )   
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
      CALL mpi_comm_free( ccomm, ierr )
      !
   END DO

   DEALLOCATE( srcv, ssnd )

#else

   CALL ZPOTRF( 'L', n, sll, ldx, info )

   IF( info /= 0 ) &
      CALL errore( " pzpotrf ", " problems computing cholesky decomposition ", ABS( info ) )

#endif

   return
END SUBROUTINE pzpotf

!  now the Double Precision subroutine

SUBROUTINE pdpotf( sll, ldx, n, desc )
   !
   use descriptors, ONLY: descla_local_dims, descla_siz_ , la_myr_ , la_myc_ , la_me_ , nlax_ , &
                          nlar_ , nlac_ , ilar_ , ilac_ , la_comm_ , la_nx_ , la_npr_ , la_npc_
   use parallel_include
   use kinds
   !
   implicit none
   !
   integer  :: n, ldx
   integer  :: desc( descla_siz_ )
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

   myrow = desc( la_myr_ )
   mycol = desc( la_myc_ )
   myid  = desc( la_me_ )
   np    = desc( la_npr_ )

   IF( desc( la_npr_ ) /= desc( la_npc_ ) ) THEN
      CALL errore( ' pdpotf ', ' only square grid are allowed ', 1 ) 
   END IF

   IF( ldx /= desc( nlax_ ) ) THEN
      CALL errore( ' pdpotf ', ' wrong leading dimension ldx ', ldx ) 
   END IF

   nr = desc( nlar_ )
   nc = desc( nlac_ )

   ALLOCATE( ssnd( ldx, ldx ) )
   ALLOCATE( srcv( ldx, ldx ) )

   DO jb = 1, np
      !
      !    Update and factorize the current diagonal block and test
      !    for non-positive-definiteness.
      !
      CALL descla_local_dims( jir, jnr, n, desc( la_nx_ ), np, jb-1 )
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
      CALL mpi_comm_split( desc( la_comm_ ) , color, key, ccomm, ierr )
      !  
      IF( myrow >= jb-1 .and. mycol <= jb-1 ) THEN
          color = myrow
          key   = mycol
      ELSE
          color = np
          key   = myid
      END IF
      !
      CALL mpi_comm_split( desc( la_comm_ ), color, key, rcomm, ierr )
      !
      !    here every process can work independently, then we need a reduce.
      !
      IF( jb > 1 ) THEN
         !
         DO ib = 1, jb - 1
            IF( ( myrow == ( jb - 1 ) ) .AND. ( mycol == ( ib - 1 ) ) ) THEN
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
            CALL errore( " pzpotrf ", " problems computing cholesky decomposition ", ABS( info ) )
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
            CALL mpi_bcast( sll,  ldx*ldx, MPI_DOUBLE_PRECISION, 0, ccomm, ierr )   
         ELSE IF( ( myrow > ( jb - 1 ) ) .AND. ( mycol < ( jb - 1 ) ) ) THEN
            CALL mpi_bcast( srcv, ldx*ldx, MPI_DOUBLE_PRECISION, 0, ccomm, ierr )   
         END IF
         !
         DO ib = jb + 1, np
            CALL descla_local_dims( iir, inr, n, desc( la_nx_ ), np, ib-1 )
            DO kb = 1, jb - 1
               CALL descla_local_dims( kic, knc, n, desc( la_nx_ ), np, kb-1 )
               IF( ( myrow == ( ib - 1 ) ) .AND. ( mycol == ( kb - 1 ) ) ) THEN
                  CALL DGEMM( 'N', 'T', inr, jnr, knc, -ONE, sll, ldx, srcv, ldx, zero, ssnd, ldx )
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
            CALL mpi_bcast( sll,  ldx*ldx, MPI_DOUBLE_PRECISION, 0, ccomm, ierr )   
         ELSE IF( ( myrow > ( jb - 1 ) ) .AND. ( mycol == ( jb - 1 ) ) ) THEN
            CALL mpi_bcast( srcv,  ldx*ldx, MPI_DOUBLE_PRECISION, 0, ccomm, ierr )   
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
      CALL mpi_comm_free( ccomm, ierr )
      !
   END DO

   DEALLOCATE( srcv, ssnd )

#else

   CALL DPOTRF( 'L', n, sll, ldx, info )

   IF( info /= 0 ) &
      CALL errore( " pzpotrf ", " problems computing cholesky decomposition ", ABS( info ) )

#endif

   return
END SUBROUTINE pdpotf

!
!
!
!

SUBROUTINE pztrtri ( sll, ldx, n, desc )
    
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

    USE kinds
    USE parallel_include
    USE descriptors, ONLY: descla_local_dims, descla_siz_ , la_myr_ , la_myc_ , la_me_ , nlax_ , &
                          nlar_ , nlac_ , ilar_ , ilac_ , la_comm_ , la_nx_ , la_npr_ , la_npc_

    IMPLICIT NONE

    INTEGER, INTENT( IN ) :: n, ldx
    INTEGER, INTENT( IN ) :: desc( descla_siz_ )
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

    myrow = desc( la_myr_ )
    mycol = desc( la_myc_ )
    myid  = desc( la_me_ )
    np    = desc( la_npr_ )
    comm  = desc( la_comm_ )

    IF( desc( la_npr_ ) /= desc( la_npc_ ) ) THEN
       CALL errore( ' pztrtri ', ' only square grid are allowed ', 1 ) 
    END IF
    IF( ldx /= desc( nlax_ ) ) THEN
       CALL errore( ' pztrtri ', ' wrong leading dimension ldx ', ldx ) 
    END IF

    nr = desc( nlar_ )
    nc = desc( nlac_ )

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
       END IF
       !
       IF( myrow == 1 .AND. mycol == 1 ) THEN
          CALL MPI_Send(sll, ldx*ldx, MPI_DOUBLE_COMPLEX, idref, 1, comm, ierr)
       END IF
       !
       IF( myrow == 1 .AND. mycol == 0 ) THEN
          !
          CALL GRID2D_RANK( 'R', np, np, 0, 0, i )
          CALL GRID2D_RANK( 'R', np, np, 1, 1, j )
          !
          CALL MPI_Irecv( B, ldx*ldx, MPI_DOUBLE_COMPLEX, i, 0, comm, req(1), ierr)
          CALL MPI_Irecv( C, ldx*ldx, MPI_DOUBLE_COMPLEX, j, 1, comm, req(2), ierr)
          !
          CALL MPI_Wait(req(1), status, ierr)
          !
          CALL zgemm('N', 'N', ldx, ldx, ldx, ONE, sll, ldx, b, ldx, ZERO, buf_recv, ldx)
          !
          CALL MPI_Wait(req(2), status, ierr)
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
       CALL MPI_Comm_size( col_comm, group_size, ierr )
       CALL MPI_Comm_rank( col_comm, group_rank, ierr )
       !
    ELSE
       !
       !  other procs stay at the window!
       !
       CALL MPI_Comm_split( comm, MPI_UNDEFINED, MPI_UNDEFINED, col_comm, ierr )
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
                send = shift(1, group_rank, 1, (group_size-1), 'S')
                recv = shift(1, group_rank, 1, (group_size-1), 'N')
                ! with the no-blocking communication the computation and the shift of the column block are overapped  
                IF( MOD( cicle, 2 ) == 0 ) THEN
                   CALL MPI_Isend(BUF_RECV, ldx*ldx, MPI_DOUBLE_COMPLEX, send, group_rank+cicle, col_comm, req(1), ierr)
                   CALL MPI_Irecv(B, ldx*ldx, MPI_DOUBLE_COMPLEX, recv, recv+cicle, col_comm, req(2), ierr)
                   CALL zgemm('N', 'N', ldx, ldx, ldx, -ONE, C, ldx, BUF_RECV, ldx, first, sll, ldx)
                ELSE
                   CALL MPI_Isend(B, ldx*ldx, MPI_DOUBLE_COMPLEX, send, group_rank+cicle, col_comm, req(1), ierr)
                   CALL MPI_Irecv(BUF_RECV, ldx*ldx, MPI_DOUBLE_COMPLEX, recv, recv+cicle, col_comm, req(2), ierr)
                   CALL zgemm('N', 'N', ldx, ldx, ldx, -ONE, C, ldx, B, ldx, ONE, sll, ldx)
                END IF
                CALL MPI_Wait(req(1), status, ierr)
                CALL MPI_Wait(req(2), status, ierr)
             END IF
             cicle = cicle + 1
             first = ONE 
          END DO
       END IF
    END DO

    IF( myrow >= mycol ) THEN
       CALL mpi_comm_free( col_comm, ierr )
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
          CALL errore( ' pztrtri ', ' problem in the local inversion ', info )
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

END SUBROUTINE pztrtri

!  now the Double Precision subroutine

SUBROUTINE pdtrtri ( sll, ldx, n, desc )
    
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

    USE kinds
    USE parallel_include
    USE descriptors, ONLY: descla_local_dims, descla_siz_ , la_myr_ , la_myc_ , la_me_ , nlax_ , &
                          nlar_ , nlac_ , ilar_ , ilac_ , la_comm_ , la_nx_ , la_npr_ , la_npc_

    IMPLICIT NONE

    INTEGER, INTENT( IN ) :: n, ldx
    INTEGER, INTENT( IN ) :: desc( descla_siz_ )
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

    myrow = desc( la_myr_ )
    mycol = desc( la_myc_ )
    myid  = desc( la_me_ )
    np    = desc( la_npr_ )
    comm  = desc( la_comm_ )

    IF( desc( la_npr_ ) /= desc( la_npc_ ) ) THEN
       CALL errore( ' pdtrtri ', ' only square grid are allowed ', 1 ) 
    END IF
    IF( ldx /= desc( nlax_ ) ) THEN
       CALL errore( ' pdtrtri ', ' wrong leading dimension ldx ', ldx ) 
    END IF

    nr = desc( nlar_ )
    nc = desc( nlac_ )

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
       END IF
       !
       IF( myrow == 1 .AND. mycol == 1 ) THEN
          CALL MPI_Send(sll, ldx*ldx, MPI_DOUBLE_PRECISION, idref, 1, comm, ierr)
       END IF
       !
       IF( myrow == 1 .AND. mycol == 0 ) THEN
          !
          CALL GRID2D_RANK( 'R', np, np, 0, 0, i )
          CALL GRID2D_RANK( 'R', np, np, 1, 1, j )
          !
          CALL MPI_Irecv( B, ldx*ldx, MPI_DOUBLE_PRECISION, i, 0, comm, req(1), ierr)
          CALL MPI_Irecv( C, ldx*ldx, MPI_DOUBLE_PRECISION, j, 1, comm, req(2), ierr)
          !
          CALL MPI_Wait(req(1), status, ierr)
          !
          CALL dgemm('N', 'N', ldx, ldx, ldx, ONE, sll, ldx, b, ldx, ZERO, buf_recv, ldx)
          !
          CALL MPI_Wait(req(2), status, ierr)
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
       CALL MPI_Comm_size( col_comm, group_size, ierr )
       CALL MPI_Comm_rank( col_comm, group_rank, ierr )
       !
    ELSE
       !
       !  other procs stay at the window!
       !
       CALL MPI_Comm_split( comm, MPI_UNDEFINED, MPI_UNDEFINED, col_comm, ierr )
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
                send = shift(1, group_rank, 1, (group_size-1), 'S')
                recv = shift(1, group_rank, 1, (group_size-1), 'N')
                ! with the no-blocking communication the computation and the shift of the column block are overapped  
                IF( MOD( cicle, 2 ) == 0 ) THEN
                   CALL MPI_Isend(BUF_RECV, ldx*ldx, MPI_DOUBLE_PRECISION, send, group_rank+cicle, col_comm, req(1), ierr)
                   CALL MPI_Irecv(B, ldx*ldx, MPI_DOUBLE_PRECISION, recv, recv+cicle, col_comm, req(2), ierr)
                   CALL dgemm('N', 'N', ldx, ldx, ldx, -ONE, C, ldx, BUF_RECV, ldx, first, sll, ldx)
                ELSE
                   CALL MPI_Isend(B, ldx*ldx, MPI_DOUBLE_PRECISION, send, group_rank+cicle, col_comm, req(1), ierr)
                   CALL MPI_Irecv(BUF_RECV, ldx*ldx, MPI_DOUBLE_PRECISION, recv, recv+cicle, col_comm, req(2), ierr)
                   CALL dgemm('N', 'N', ldx, ldx, ldx, -ONE, C, ldx, B, ldx, ONE, sll, ldx)
                END IF
                CALL MPI_Wait(req(1), status, ierr)
                CALL MPI_Wait(req(2), status, ierr)
             END IF
             cicle = cicle + 1
             first = ONE 
          END DO
       END IF
    END DO

    IF( myrow >= mycol ) THEN
       CALL mpi_comm_free( col_comm, ierr )
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
          CALL errore( ' pdtrtri ', ' problem in the local inversion ', info )
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

END SUBROUTINE pdtrtri
