!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

   MODULE parallel_types
      USE kinds
      IMPLICIT NONE
      PRIVATE
      SAVE


      TYPE processors_grid
         INTEGER :: context  !  Communication handle, grid identification
         INTEGER :: nproc    !  number of processors in the grid
         INTEGER :: my_pe    !  process index (0 ... nproc -1)
         INTEGER :: npx      !  Grid dimensions :  
         INTEGER :: npy      !  (nprows, npcolumns, npplanes)
         INTEGER :: npz      !  
         INTEGER :: mex      !  Processor coordinates:
         INTEGER :: mey      !  (mex, mey, mez)
         INTEGER :: mez      !  0 <= mex < npx-1
                             !  0 <= mey < npy-1
                             !  0 <= mez < npz-1
      END TYPE

        ! ...   Valid values for data distribution
        !
        INTEGER, PARAMETER :: BLOCK_CYCLIC_DIST    = 1
        INTEGER, PARAMETER :: BLOCK_PARTITION_DIST = 2
        INTEGER, PARAMETER :: FREE_PATTERN_DIST    = 3
        INTEGER, PARAMETER :: REPLICATED_DATA_DIST = 4
        INTEGER, PARAMETER :: CYCLIC_DIST          = 5

!  ----------------------------------------------
!  BEGIN manual
!
!  Given the Array  |a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11|
!  and three processors P0, P1, P2 
!
!  in the BLOCK_PARTITION_DIST scheme, the Array is partitioned 
!  as follow
!       P0            P1            P2
!  |a1 a2 a3 a4| |a5 a6 a7 a8| |a9 a10 a11|
!
!  in the BLOCK_CYCLIC_DIST scheme the Array is first partitioned 
!  into blocks (i.e. of size 2)  |a1 a2|a3 a4|a5 a6|a7 a8|a9 a10|a11|
!  Then the block are distributed cyclically among P0, P1 and P2
!       P0             P1              P2
!  |a1 a2|a7 a8|  |a3 a4|a9 a10|  |a5 a6|a11|
!
!  in the CYCLIC_DIST scheme the Array elements are distributed round robin
!  among P0, P1 and P2
!       P0             P1              P2
!  |a1 a4 a7 a10|  |a2 a5 a8 a11|  |a3 a6 a9|
!
!  ----------------------------------------------
!  END manual



        TYPE descriptor
          INTEGER :: matrix_type     ! = 1, for dense matrices
          TYPE (processors_grid) :: grid ! Communication handle
          INTEGER :: nx     ! rows, number of rows in the global array
          INTEGER :: ny     ! columns, number of columns in the global array
          INTEGER :: nz     ! planes, number of planes in the global array
          INTEGER :: nxblk  ! row_block, if DIST = BLOCK_CICLYC_DIST,
                            ! this value represent the blocking factor
                            ! used to distribute the rows of the array,
                            ! otherwise this is the size of local block of rows
          INTEGER :: nyblk  ! column_block, same as row_block but for columns
          INTEGER :: nzblk  ! plane_block, same as row_block but for planes
          INTEGER :: nxl    ! local_rows, number of rows in the local array
          INTEGER :: nyl    ! local_columns, number of columns in the local array
          INTEGER :: nzl    ! local_planes, number of planes in the local array
          INTEGER :: ixl    ! irow
          INTEGER :: iyl    ! icolumn
          INTEGER :: izl    ! iplane
          INTEGER :: ipexs  ! row_src_pe, process row over which the first row 
                            ! of the array is distributed
          INTEGER :: ipeys  ! column_src_pe, process column over which the first column
                            ! of the array is distributed
          INTEGER :: ipezs  ! plane_src_pe, process plane over which the first plane
                            ! of the array is distributed
          INTEGER :: ldx    ! local_ld, leading dimension of the local sub-block of the array
          INTEGER :: ldy    ! local_sub_ld, sub-leading dimension of the local sub-block
                            ! of the array
          INTEGER :: ldz    ! 

          INTEGER :: xdist ! row_dist
          INTEGER :: ydist ! column_dist
          INTEGER :: zdist ! plane_dist
        END TYPE
        
        
        PUBLIC :: processors_grid

        PUBLIC ::  BLOCK_CYCLIC_DIST, BLOCK_PARTITION_DIST, &
          FREE_PATTERN_DIST, REPLICATED_DATA_DIST, CYCLIC_DIST

        INTEGER NUMROC
        EXTERNAL NUMROC


      END MODULE parallel_types
