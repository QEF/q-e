!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      MODULE descriptors_module
        USE parallel_types
        USE io_global,  ONLY : stdout
        IMPLICIT NONE
        SAVE

        INTEGER  ldim_block, ldim_cyclic, ldim_block_cyclic
        INTEGER  lind_block, lind_cyclic, lind_block_cyclic
        EXTERNAL ldim_block, ldim_cyclic, ldim_block_cyclic
        EXTERNAL lind_block, lind_cyclic, lind_block_cyclic
       
        CONTAINS


          SUBROUTINE desc_init_blacs(desc, matrix_type, rows, columns, &
            row_block, column_block, row_src_pe, column_src_pe, grid, local_ld)

              TYPE (descriptor) :: desc
              INTEGER, INTENT(IN) :: matrix_type    
              TYPE (processors_grid), INTENT(IN) :: grid 
              INTEGER, INTENT(IN) :: rows   
              INTEGER, INTENT(IN) :: columns
              INTEGER, INTENT(IN) :: row_block
              INTEGER, INTENT(IN) :: column_block
              INTEGER, INTENT(IN) :: row_src_pe   
              INTEGER, INTENT(IN) :: column_src_pe
              INTEGER, INTENT(IN), OPTIONAL :: local_ld

              desc%matrix_type = matrix_type
              desc%grid = grid

              CALL desc_init_x(desc%nx, desc%xdist, desc%nxl, &
                desc%nxblk, desc%ixl, desc%ipexs, rows, &
                BLOCK_CYCLIC_DIST, row_block, row_src_pe, grid%mex, &
                grid%npx)
              CALL desc_init_x(desc%ny, desc%ydist, &
                desc%nyl, desc%nyblk, desc%iyl, &
                desc%ipeys, columns, BLOCK_CYCLIC_DIST, column_block, &
                column_src_pe, grid%mey, grid%npy)

              desc%nz = 1
              desc%nzl = 1
              desc%nzblk = 1
              desc%ipezs = 0
              desc%zdist = REPLICATED_DATA_DIST

              IF(PRESENT(local_ld)) THEN
                desc%ldx = local_ld
              ELSE
                desc%ldx = ldim_block_cyclic( rows, row_block, grid%npx, grid%mex )
              END IF
              desc%ldy = 1

            RETURN
          END SUBROUTINE desc_init_blacs


          SUBROUTINE desc_init_x(desc_nxs, desc_nx_dist, desc_local_nxs, &
            desc_nx_block, desc_ix, desc_nx_src_pe, nxs, nx_dist, nx_block, &
            nx_src_pe, mype, npes)
!

              IMPLICIT NONE
              INTEGER, INTENT(OUT) ::  desc_nxs
              INTEGER, INTENT(OUT) ::  desc_nx_dist
              INTEGER, INTENT(OUT) ::  desc_local_nxs
              INTEGER, INTENT(OUT) ::  desc_nx_block
              INTEGER, INTENT(OUT) ::  desc_ix
              INTEGER, INTENT(OUT) ::  desc_nx_src_pe
              INTEGER, INTENT(IN)  ::  nxs
              INTEGER, INTENT(IN)  ::  nx_dist
              INTEGER, INTENT(IN)  ::  nx_block
              INTEGER, INTENT(IN)  ::  nx_src_pe
              INTEGER, INTENT(IN)  ::  mype
              INTEGER, INTENT(IN)  ::  npes
                
              desc_nxs      = nxs
              desc_nx_dist  = nx_dist

              SELECT CASE (nx_dist)
              CASE ( BLOCK_CYCLIC_DIST )
                desc_local_nxs = ldim_block_cyclic( nxs, nx_block, npes, mype )
                desc_ix        = lind_block_cyclic( 1, nxs, nx_block, npes, mype)
                desc_nx_block  = nx_block
                desc_nx_src_pe = nx_src_pe
              CASE ( BLOCK_PARTITION_DIST )
                desc_local_nxs = ldim_block( nxs, npes, mype )
                desc_ix        = lind_block( 1, nxs, npes, mype)
                desc_nx_block  = desc_local_nxs 
                desc_nx_src_pe = 0             
              CASE ( CYCLIC_DIST )
                desc_local_nxs = ldim_cyclic( nxs, npes, mype )
                desc_ix        = lind_cyclic( 1, nxs, npes, mype)
                desc_nx_block  = 1 
                desc_nx_src_pe = 0             
              CASE ( REPLICATED_DATA_DIST )
                desc_local_nxs = nxs
                desc_ix        = 1
                desc_nx_block  = nxs
                desc_nx_src_pe = mype
              END SELECT

            RETURN
          END SUBROUTINE desc_init_x

!=----------------------------------------------------------------------------=!
!  BEGIN manual

          SUBROUTINE pblas_descriptor(pb_desc, desc)
!
!  END manual
!=----------------------------------------------------------------------------=!

            INTEGER :: pb_desc(:)
            TYPE (descriptor) :: desc
                pb_desc(1) = desc%matrix_type
                pb_desc(2) = desc%grid%context
                pb_desc(3) = desc%nx
                pb_desc(4) = desc%ny
                pb_desc(5) = desc%nxblk
                pb_desc(6) = desc%nyblk
                pb_desc(7) = desc%ipexs
                pb_desc(8) = desc%ipeys
                pb_desc(9) = desc%ldx
            RETURN
          END SUBROUTINE pblas_descriptor


      END MODULE descriptors_module
