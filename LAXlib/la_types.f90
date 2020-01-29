!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

   MODULE laxlib_descriptor
      !
      IMPLICIT NONE
      SAVE

      INTEGER, EXTERNAL ::  ldim_block, ldim_cyclic, ldim_block_sca
      INTEGER, EXTERNAL ::  gind_block, gind_cyclic, gind_block_sca

      !  Descriptor for linear algebra data distribution (like in Cannon's algorithm)
      !
      !  Remember here we use square matrixes block distributed on a square grid of processors
      !
      TYPE la_descriptor  
         INTEGER :: ir        = 0 !  global index of the first row in the local block of the distributed matrix
         INTEGER :: nr        = 0 !  number of row in the local block of the distributed matrix
         INTEGER :: ic        = 0 !  global index of the first column in the local block of the distributed matrix
         INTEGER :: nc        = 0 !  number of column in the local block of the distributed matrix
         INTEGER :: nrcx        = 0 !  leading dimension of the distribute matrix (greather than nr and nc)
         INTEGER :: active_node = 0 !  if > 0 the proc holds a block of the lambda matrix
         INTEGER :: n        = 0 !  global dimension of the matrix
         INTEGER :: nx       = 0 !  global leading dimension ( >= n )
         INTEGER :: npr      = 0 !  number of row processors 
         INTEGER :: npc      = 0 !  number of column processors 
         INTEGER :: myr      = 0 !  processor row index
         INTEGER :: myc      = 0 !  processor column index
         INTEGER :: comm     = 0 !  communicator
         INTEGER :: cntx     =-1 !  scalapack context
         INTEGER :: mype     = 0 !  processor index ( from 0 to desc( la_npr_ ) * desc( la_npc_ ) - 1 )
         INTEGER :: nrl      = 0 !  number of local rows, when the matrix rows are cyclically distributed across proc
         INTEGER :: nrlx     = 0 !  leading dimension, when the matrix is distributed by row
      END TYPE
      !
   CONTAINS
   !
   !
   SUBROUTINE descla_init( descla, n, nx, np, me, comm, cntx, includeme )
      !
      IMPLICIT NONE  
      TYPE(la_descriptor), INTENT(OUT) :: descla
      INTEGER, INTENT(IN)  :: n   !  the size of this matrix
      INTEGER, INTENT(IN)  :: nx  !  the max among different matrixes sharing 
                                  !  this descriptor or the same data distribution
      INTEGER, INTENT(IN)  :: np(2), me(2), comm, cntx
      INTEGER, INTENT(IN)  :: includeme
      INTEGER  :: ir, nr, ic, nc, lnode, nrcx, nrl, nrlx
      INTEGER  :: ip, npp
      
      IF( np(1) /= np(2) ) &
         CALL lax_error__( ' descla_init ', ' only square grid of proc are allowed ', 2 )
      IF( n < 0 ) &
         CALL lax_error__( ' descla_init ', ' dummy argument n less than 1 ', 3 )
      IF( nx < n ) &
         CALL lax_error__( ' descla_init ', ' dummy argument nx less than n ', 4 )
      IF( np(1) < 1 ) &
         CALL lax_error__( ' descla_init ', ' dummy argument np less than 1 ', 5 )

      ! find the block maximum dimensions

#if __SCALAPACK
      nrcx = ldim_block_sca( nx, np(1), 0 )
      descla%cntx = cntx
#else
      nrcx = ldim_block( nx, np(1), 0 )
      DO ip = 1, np(1) - 1
         nrcx = MAX( nrcx, ldim_block( nx, np(1), ip ) )
      END DO
      descla%cntx = -1
#endif
      !
      ! find local dimensions, if appropriate
      !
      IF( includeme == 1 ) THEN
         !
         CALL descla_local_dims( ir, nr, n, nx, np(1), me(1) )
         CALL descla_local_dims( ic, nc, n, nx, np(2), me(2) )
         !
         lnode = 1
         !
      ELSE
         !
         nr = 0
         nc = 0
         !  
         ir = 0
         ic = 0
         !
         lnode = -1
         !
      END IF

      descla%ir = ir    ! global index of the first row in the local block of lambda
      descla%nr = nr    ! number of row in the local block of lambda ( the "2" accounts for spin)
      descla%ic = ic    ! global index of the first column in the local block of lambda
      descla%nc = nc    ! number of column in the local block of lambda
      descla%nrcx = nrcx  ! leading dimension of the distribute lambda matrix
      descla%active_node = lnode 
                          !  if > 0 the proc holds a block of the lambda matrix
      descla%n = n        ! global dimension of the matrix
      descla%nx = nx      ! global leading dimension
      descla%npr = np(1)  ! number of row processors 
      descla%npc = np(2)  ! number of column processors 
      descla%myr = me(1)  ! processor row index
      descla%myc = me(2)  ! processor column index
      descla%comm = comm  ! communicator
      descla%mype = descla%myc + descla%myr * descla%npr 
                          ! processor index ( from 0 to desc( la_npr_ ) * desc( la_npc_ ) - 1 )
     
      npp = np(1) * np(2)

      !  Compute local dimension of the cyclically distributed matrix
      !
      IF( includeme == 1 ) THEN
         nrl  = ldim_cyclic( n, npp, descla%mype )
      ELSE
         nrl = 0
      END IF
      nrlx = n / npp + 1

      descla%nrl  = nrl  ! number of local rows, when the matrix rows are cyclically distributed across procs
      descla%nrlx = nrlx ! leading dimension 

      IF( nr < 0 .OR. nc < 0 ) &
         CALL lax_error__( ' descla_init ', ' wrong valune for computed nr and nc ', 1 )
      IF( nrcx < 1 ) &
         CALL lax_error__( ' descla_init ', ' wrong value for computed nrcx ', 2 )
      IF( nrcx < nr ) &
         CALL lax_error__( ' descla_init ', ' nrcx < nr ', ( nr - nrcx ) )
      IF( nrcx < nc ) &
         CALL lax_error__( ' descla_init ', ' nrcx < nc ', ( nc - nrcx ) )
      IF( nrlx < nrl ) &
         CALL lax_error__( ' descla_init ', ' nrlx < nrl ', ( nrl - nrlx ) )
      IF( nrl < 0 ) &
         CALL lax_error__( ' descla_init ', ' nrl < 0 ', ABS( nrl ) )


      RETURN
   END SUBROUTINE descla_init


   SUBROUTINE laxlib_desc_to_intarray( idesc, descla )
      IMPLICIT NONE  
      include 'laxlib_param.fh'
      TYPE(la_descriptor), INTENT(IN) :: descla
      INTEGER, INTENT(OUT)  :: idesc(LAX_DESC_SIZE)
      idesc(LAX_DESC_IR) = descla%ir        !  global index of the first row in the local block of the distributed matrix
      idesc(LAX_DESC_NR) = descla%nr        !  number of row in the local block of the distributed matrix
      idesc(LAX_DESC_IC) = descla%ic        !  global index of the first column in the local block of the distributed matrix
      idesc(LAX_DESC_NC) = descla%nc        !  number of column in the local block of the distributed matrix
      idesc(LAX_DESC_NRCX) = descla%nrcx      !  leading dimension of the distribute matrix (greather than nr and nc)
      idesc(LAX_DESC_ACTIVE_NODE) = descla%active_node !  if > 0 the proc holds a block of the lambda matrix
      idesc(LAX_DESC_N) = descla%n        !  global dimension of the matrix
      idesc(LAX_DESC_NX) = descla%nx       !  global leading dimension ( >= n )
      idesc(LAX_DESC_NPR) = descla%npr      !  number of row processors 
      idesc(LAX_DESC_NPC) = descla%npc      !  number of column processors 
      idesc(LAX_DESC_MYR) = descla%myr      !  processor row index
      idesc(LAX_DESC_MYC) = descla%myc      !  processor column index
      idesc(LAX_DESC_COMM) = descla%comm     !  communicator
      idesc(LAX_DESC_CNTX) = descla%cntx     !  scalapack context
      idesc(LAX_DESC_MYPE) = descla%mype     !  processor index ( from 0 to desc( la_npr_ ) * desc( la_npc_ ) - 1 )
      idesc(LAX_DESC_NRL) = descla%nrl      !  number of local rows, when the matrix rows are cyclically distributed across proc
      idesc(LAX_DESC_NRLX) = descla%nrlx     !  leading dimension, when the matrix is distributed by row
      RETURN
   END SUBROUTINE 
   SUBROUTINE laxlib_intarray_to_desc( descla, idesc )
      IMPLICIT NONE  
      include 'laxlib_param.fh'
      TYPE(la_descriptor), INTENT(OUT) :: descla
      INTEGER, INTENT(IN)  :: idesc(LAX_DESC_SIZE)
      descla%ir = idesc(LAX_DESC_IR) !  global index of the first row in the local block of the distributed matrix
      descla%nr = idesc(LAX_DESC_NR) !  number of row in the local block of the distributed matrix
      descla%ic = idesc(LAX_DESC_IC) !  global index of the first column in the local block of the distributed matrix
      descla%nc = idesc(LAX_DESC_NC) !  number of column in the local block of the distributed matrix
      descla%nrcx = idesc(LAX_DESC_NRCX) !  leading dimension of the distribute matrix (greather than nr and nc)
      descla%active_node = idesc(LAX_DESC_ACTIVE_NODE) !  if > 0 the proc holds a block of the lambda matrix
      descla%n = idesc(LAX_DESC_N) !  global dimension of the matrix
      descla%nx = idesc(LAX_DESC_NX) !  global leading dimension ( >= n )
      descla%npr = idesc(LAX_DESC_NPR) !  number of row processors 
      descla%npc = idesc(LAX_DESC_NPC) !  number of column processors 
      descla%myr = idesc(LAX_DESC_MYR) !  processor row index
      descla%myc = idesc(LAX_DESC_MYC) !  processor column index
      descla%comm = idesc(LAX_DESC_COMM) !  communicator
      descla%cntx = idesc(LAX_DESC_CNTX) !  scalapack context
      descla%mype = idesc(LAX_DESC_MYPE) !  processor index ( from 0 to desc( la_npr_ ) * desc( la_npc_ ) - 1 )
      descla%nrl = idesc(LAX_DESC_NRL) !  number of local rows, when the matrix rows are cyclically distributed across proc
      descla%nrlx = idesc(LAX_DESC_NRLX) !  leading dimension, when the matrix is distributed by row
      RETURN
   END SUBROUTINE 

   END MODULE laxlib_descriptor
