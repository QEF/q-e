!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      MODULE descriptors_module
        USE parallel_types
        IMPLICIT NONE
        SAVE

        INTERFACE desc_init
          MODULE PROCEDURE desc_init_1d, desc_init_2d, desc_init_3d
        END INTERFACE
        INTERFACE global_index
          MODULE PROCEDURE globalindex_desc, globalindex_shape
        END INTERFACE
        INTERFACE local_index
          MODULE PROCEDURE localindex_desc, localindex_shape
        END INTERFACE
        INTERFACE local_dimension
          MODULE PROCEDURE localdim_desc, localdim_shape
        END INTERFACE
        INTERFACE owner_of
          MODULE PROCEDURE ownerof_desc, ownerof_shape
        END INTERFACE
        INTERFACE get_local_dims
          MODULE PROCEDURE desc_ldims
        END INTERFACE
        INTERFACE get_global_dims
          MODULE PROCEDURE desc_gdims
        END INTERFACE

        INTEGER NUMROC
        EXTERNAL NUMROC
       
        CONTAINS


!=----------------------------------------------------------------------------=!
!  BEGIN manual
!
          SUBROUTINE desc_init_1d(desc, matrix_type, rows, &
            row_block, row_src_pe, grid, row_shape)
!
!  END manual
!=----------------------------------------------------------------------------=!

              TYPE (descriptor) :: desc
              INTEGER, INTENT(IN) :: matrix_type    
              TYPE (processors_grid), INTENT(IN) :: grid 
              INTEGER, INTENT(IN) :: rows   
              INTEGER, INTENT(IN) :: row_block
              INTEGER, INTENT(IN) :: row_src_pe   
              INTEGER, INTENT(IN) :: row_shape   

              desc%matrix_type = matrix_type
              desc%grid = grid

              CALL desc_init_x(desc%nx, desc%xshape, desc%nxl, &
                desc%nxblk, desc%ixl, desc%ipexs, rows, row_shape, &
                row_block, row_src_pe, grid%mex, grid%npx)

              desc%ny = 1
              desc%nyl = 1
              desc%nyblk = 1
              desc%ipeys = 0
              desc%yshape = REPLICATED_DATA_SHAPE

              desc%nz = 1
              desc%nzl = 1
              desc%nzblk = 1
              desc%ipezs = 0
              desc%zshape = REPLICATED_DATA_SHAPE

              desc%ldx = 1
              desc%ldy = 1
            RETURN
          END SUBROUTINE desc_init_1d

!=----------------------------------------------------------------------------=!
!  BEGIN manual

          SUBROUTINE desc_init_blacs(desc, matrix_type, rows, columns, &
            row_block, column_block, row_src_pe, column_src_pe, grid, local_ld)
!
!
!  END manual
!=----------------------------------------------------------------------------=!

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

              CALL desc_init_x(desc%nx, desc%xshape, desc%nxl, &
                desc%nxblk, desc%ixl, desc%ipexs, rows, &
                BLOCK_CYCLIC_SHAPE, row_block, row_src_pe, grid%mex, &
                grid%npx)
              CALL desc_init_x(desc%ny, desc%yshape, &
                desc%nyl, desc%nyblk, desc%iyl, &
                desc%ipeys, columns, BLOCK_CYCLIC_SHAPE, column_block, &
                column_src_pe, grid%mey, grid%npy)

              desc%nz = 1
              desc%nzl = 1
              desc%nzblk = 1
              desc%ipezs = 0
              desc%zshape = REPLICATED_DATA_SHAPE

              IF(PRESENT(local_ld)) THEN
                desc%ldx = local_ld
              ELSE
                desc%ldx = localdim_shape( rows, row_block, grid%mex, &
                  row_src_pe, grid%npx, desc%xshape)
              END IF
              desc%ldy = 1

            RETURN
          END SUBROUTINE desc_init_blacs

!=----------------------------------------------------------------------------=!
!  BEGIN manual

          SUBROUTINE desc_init_2d(desc, matrix_type, rows, columns, &
            row_block, column_block, row_src_pe, &
            column_src_pe, grid, row_shape, column_shape, local_ld)
!
!  END manual
!=----------------------------------------------------------------------------=!

              TYPE (descriptor) :: desc
              INTEGER, INTENT(IN) :: matrix_type    
              TYPE (processors_grid), INTENT(IN) :: grid 
              INTEGER, INTENT(IN) :: rows   
              INTEGER, INTENT(IN) :: columns
              INTEGER, INTENT(IN) :: row_block
              INTEGER, INTENT(IN) :: column_block
              INTEGER, INTENT(IN) :: row_src_pe   
              INTEGER, INTENT(IN) :: column_src_pe
              INTEGER, INTENT(IN) :: row_shape   
              INTEGER, INTENT(IN) :: column_shape
              INTEGER, INTENT(IN), OPTIONAL :: local_ld
 
              LOGICAL :: debug = .FALSE.

              desc%matrix_type = matrix_type
              desc%grid = grid

              CALL desc_init_x(desc%nx, desc%xshape, desc%nxl, &
                desc%nxblk, desc%ixl, desc%ipexs, rows, row_shape, &
                row_block, row_src_pe, grid%mex, grid%npx)
              CALL desc_init_x(desc%ny, desc%yshape, &
                desc%nyl, desc%nyblk, desc%iyl, &
                desc%ipeys, columns, column_shape, column_block, &
                column_src_pe, grid%mey, grid%npy)

              IF( debug ) THEN
                WRITE(6,fmt="(' desc%nx      = ', I6 )") desc%nx
                WRITE(6,fmt="(' desc%xshape  = ', I6 )") desc%xshape
                WRITE(6,fmt="(' desc%nxl     = ', I6 )") desc%nxl
                WRITE(6,fmt="(' desc%nxblk   = ', I6 )") desc%nxblk
                WRITE(6,fmt="(' desc%ixl     = ', I6 )") desc%ixl
                WRITE(6,fmt="(' desc%ipexs   = ', I6 )") desc%ipexs

                WRITE(6,fmt="(' desc%ny      = ', I6 )") desc%ny
                WRITE(6,fmt="(' desc%yshape  = ', I6 )") desc%yshape
                WRITE(6,fmt="(' desc%nyl     = ', I6 )") desc%nyl
                WRITE(6,fmt="(' desc%nyblk   = ', I6 )") desc%nyblk
                WRITE(6,fmt="(' desc%iyl     = ', I6 )") desc%iyl
                WRITE(6,fmt="(' desc%ipeys   = ', I6 )") desc%ipeys
              END IF

              desc%nz = 1
              desc%nzl = 1
              desc%nzblk = 1
              desc%ipezs = 0
              desc%zshape = REPLICATED_DATA_SHAPE

              IF(PRESENT(local_ld)) THEN
                desc%ldx = local_ld
              ELSE
                desc%ldx = localdim_shape( rows, row_block, grid%mex, &
                  row_src_pe, grid%npx, desc%xshape)
              END IF
              desc%ldy = 1
            RETURN
          END SUBROUTINE desc_init_2d


!=----------------------------------------------------------------------------=!
!  BEGIN manual

          SUBROUTINE desc_init_3d(desc, matrix_type, rows, columns, &
            planes, row_block, column_block, plane_block, row_src_pe, &
            column_src_pe, plane_src_pe, grid, row_shape, column_shape, &
            plane_shape, local_ld, local_sub_ld)
!
!  END manual
!=----------------------------------------------------------------------------=!

              TYPE (descriptor) :: desc
              INTEGER, INTENT(IN) :: matrix_type    
              TYPE (processors_grid), INTENT(IN) :: grid 
              INTEGER, INTENT(IN) :: rows   
              INTEGER, INTENT(IN) :: columns
              INTEGER, INTENT(IN) :: planes
              INTEGER, INTENT(IN) :: row_block
              INTEGER, INTENT(IN) :: column_block
              INTEGER, INTENT(IN) :: plane_block
              INTEGER, INTENT(IN) :: row_src_pe   
              INTEGER, INTENT(IN) :: column_src_pe
              INTEGER, INTENT(IN) :: plane_src_pe
              INTEGER, INTENT(IN) :: row_shape   
              INTEGER, INTENT(IN) :: column_shape
              INTEGER, INTENT(IN) :: plane_shape
              INTEGER, INTENT(IN), OPTIONAL :: local_ld
              INTEGER, INTENT(IN), OPTIONAL :: local_sub_ld

              desc%matrix_type = matrix_type
              desc%grid = grid

              CALL desc_init_x(desc%nx, desc%xshape, desc%nxl, &
                desc%nxblk, desc%ixl, desc%ipexs, rows, row_shape, &
                row_block, row_src_pe, grid%mex, grid%npx)
              CALL desc_init_x(desc%ny, desc%yshape, &
                desc%nyl, desc%nyblk, desc%iyl, &
                desc%ipeys, columns, column_shape, column_block, &
                column_src_pe, grid%mey, grid%npy)
              CALL desc_init_x(desc%nz, desc%zshape, &
                desc%nzl, desc%nzblk, desc%izl, &
                desc%ipezs, planes, plane_shape, plane_block, &
                plane_src_pe, grid%mez, grid%npz)

              IF(PRESENT(local_ld)) THEN
                desc%ldx = local_ld
              ELSE
                desc%ldx = localdim_shape( rows, row_block, grid%mex, &
                  row_src_pe, grid%npx, desc%xshape)
              END IF
              IF(PRESENT(local_sub_ld)) THEN
                desc%ldy = local_sub_ld
              ELSE
                desc%ldy =  localdim_shape( columns, column_block, &
                  grid%mey, column_src_pe, grid%npy, &
                  desc%yshape)
              END IF
            RETURN
          END SUBROUTINE desc_init_3d

!=----------------------------------------------------------------------------=!
!  BEGIN manual

          SUBROUTINE desc_init_x(desc_nxs, desc_nx_shape, desc_local_nxs, &
            desc_nx_block, desc_ix, desc_nx_src_pe, nxs, nx_shape, nx_block, &
            nx_src_pe, mype, npes)
!
!  END manual
!=----------------------------------------------------------------------------=!

              IMPLICIT NONE
              INTEGER, INTENT(OUT) ::  desc_nxs
              INTEGER, INTENT(OUT) ::  desc_nx_shape
              INTEGER, INTENT(OUT) ::  desc_local_nxs
              INTEGER, INTENT(OUT) ::  desc_nx_block
              INTEGER, INTENT(OUT) ::  desc_ix
              INTEGER, INTENT(OUT) ::  desc_nx_src_pe
              INTEGER, INTENT(IN)  ::  nxs
              INTEGER, INTENT(IN)  ::  nx_shape
              INTEGER, INTENT(IN)  ::  nx_block
              INTEGER, INTENT(IN)  ::  nx_src_pe
              INTEGER, INTENT(IN)  ::  mype
              INTEGER, INTENT(IN)  ::  npes
                
              desc_nxs       = nxs
              desc_nx_shape  = nx_shape
              desc_local_nxs = localdim_shape( nxs, nx_block, mype, nx_src_pe, npes, nx_shape)
              desc_ix        = localindex_shape( 1, nxs, nx_block, mype, npes, nx_shape)

              SELECT CASE (nx_shape)
              CASE ( BLOCK_CYCLIC_SHAPE )
                desc_nx_block = nx_block
                desc_nx_src_pe = nx_src_pe
              CASE ( BLOCK_PARTITION_SHAPE )
                desc_nx_block  = desc_local_nxs 
                desc_nx_src_pe = 0             
              CASE ( CYCLIC_SHAPE )
                desc_nx_block  = 1 
                desc_nx_src_pe = 0             
              CASE ( REPLICATED_DATA_SHAPE )
                desc_nx_block = nxs
                desc_nx_src_pe = mype
              END SELECT

            RETURN
          END SUBROUTINE

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

!=----------------------------------------------------------------------------=!
!  BEGIN manual

          INTEGER FUNCTION globalindex_shape( lind, n, nb, me, isrc, np, pshape )

!  This function computes the global index of a distributed array entry
!  pointed to by the local index lind of the process indicated by me.
!  lind      local index of the distributed matrix entry.
!  N         is the size of the global array.
!  NB        size of the blocks the distributed matrix is split into.
!  me        The coordinate of the process whose local array row or
!            column is to be determined.
!  isrc      The coordinate of the process that possesses the first
!            row/column of the distributed matrix.
!  np        The total number processes over which the distributed
!            matrix is distributed.
!
!  END manual
!=----------------------------------------------------------------------------=!


            INTEGER, INTENT(IN) :: lind, n, nb, me, isrc, np, pshape
            INTEGER r, q

            IF( pshape .EQ. BLOCK_PARTITION_SHAPE ) THEN

              q = INT(n/np)
              r = MOD(n,np)
              IF( me < r ) THEN
                GLOBALINDEX_SHAPE = (Q+1)*me + lind
              ELSE
                GLOBALINDEX_SHAPE = Q*me + R + lind
              END IF

            ELSE IF ( pshape .EQ. BLOCK_CYCLIC_SHAPE ) THEN

              GLOBALINDEX_SHAPE = np*NB*((lind-1)/NB) + &
                MOD(lind-1,NB) + MOD(np+me-isrc, np)*NB + 1

            ELSE IF ( pshape .EQ. CYCLIC_SHAPE ) THEN

                GLOBALINDEX_SHAPE = (lind-1) * np + me + 1

            ELSE

              GLOBALINDEX_SHAPE = lind

            END IF
            RETURN
          END FUNCTION

!=----------------------------------------------------------------------------=!
!  BEGIN manual

          INTEGER FUNCTION globalindex_desc( lind, desc, what )

!  END manual
!=----------------------------------------------------------------------------=!

            INTEGER, INTENT(IN) :: lind
            TYPE (descriptor) :: desc
            CHARACTER(LEN=*) :: what
            INTEGER N, nb, src_pe, my_pe, np, pshape
            IF ( what(1:1) .EQ. 'R' .OR. what(1:1) .EQ. 'r' ) THEN
              NB     = desc%nxblk;   N      = desc%nx; 
              np     = desc%grid%npx; src_pe = desc%ipexs; 
              my_pe  = desc%grid%mex; pshape = desc%xshape
            ELSE IF ( what(1:1) .EQ. 'C' .OR. what(1:1) .EQ. 'c' ) THEN
              NB     = desc%nyblk;   N      = desc%ny; 
              np     = desc%grid%npy; src_pe = desc%ipeys; 
              my_pe  = desc%grid%mey; pshape = desc%yshape
            ELSE IF ( what(1:1) .EQ. 'P' .OR. what(1:1) .EQ. 'p' ) THEN
              NB     = desc%nzblk;   N      = desc%nz; 
              np     = desc%grid%npz; src_pe = desc%ipezs; 
              my_pe  = desc%grid%mez; pshape = desc%zshape
            END IF
            globalindex_desc = globalindex_shape(lind, n, nb, my_pe, src_pe, np, pshape )
            RETURN
          END FUNCTION

!=----------------------------------------------------------------------------=!
!  BEGIN manual

          INTEGER FUNCTION localdim_shape( n, nb, me, isrc, np, pshape)
            
!  N     = Global dimension of the array
!  NB    = Size of the blocks ( meaningful only for BLOCK_CYCLIC_SHAPE )
!  me = Index of the callig processor
!  isrc = Index of the processor owning the first element of the array
!  np = Number of processors among which the array is subdivided
!  pshape  = Shape of the distributed data
!
!  This function return the number of array elements owned
!  by the callig processor
!
!  END manual
!=----------------------------------------------------------------------------=!

           IMPLICIT NONE
           INTEGER, INTENT(IN) :: n, nb, me, isrc, np, pshape

           IF( pshape .EQ. BLOCK_PARTITION_SHAPE ) THEN

             LOCALDIM_SHAPE = INT(N/np)
             IF( me < MOD(N,np) ) LOCALDIM_SHAPE = LOCALDIM_SHAPE + 1

           ELSE IF( pshape .EQ. BLOCK_CYCLIC_SHAPE ) THEN

             LOCALDIM_SHAPE = NUMROC( N, NB, me, isrc, np )

           ELSE IF( pshape .EQ. CYCLIC_SHAPE ) THEN

             LOCALDIM_SHAPE = INT(N/np)
             IF( me < MOD(N,np) ) LOCALDIM_SHAPE = LOCALDIM_SHAPE + 1

           ELSE

             LOCALDIM_SHAPE = n

           END IF
           RETURN
         END FUNCTION

!=----------------------------------------------------------------------------=!
!  BEGIN manual

          INTEGER FUNCTION localdim_desc( desc, what )

!  END manual
!=----------------------------------------------------------------------------=!

            TYPE (descriptor) :: desc
            CHARACTER(LEN=*) :: what
            INTEGER n, nb, src_pe, my_pe, np, pshape
            IF ( what(1:1) .EQ. 'R' .OR. what(1:1) .EQ. 'r' ) THEN
              NB     = desc%nxblk;   N      = desc%nx; 
              np     = desc%grid%npx; src_pe = desc%ipexs; 
              my_pe  = desc%grid%mex; pshape = desc%xshape
            ELSE IF ( what(1:1) .EQ. 'C' .OR. what(1:1) .EQ. 'c' ) THEN
              NB     = desc%nyblk;   N      = desc%ny; 
              np     = desc%grid%npy; src_pe = desc%ipeys; 
              my_pe  = desc%grid%mey; pshape = desc%yshape
            ELSE IF ( what(1:1) .EQ. 'P' .OR. what(1:1) .EQ. 'p' ) THEN
              NB     = desc%nzblk;   N      = desc%nz; 
              np     = desc%grid%npz; src_pe = desc%ipezs; 
              my_pe  = desc%grid%mez; pshape = desc%zshape
            END IF
            localdim_desc = localdim_shape( N, NB, my_pe, src_pe, np, pshape)
            RETURN
          END FUNCTION

!=----------------------------------------------------------------------------=!
!  BEGIN manual

         INTEGER FUNCTION  localindex_shape(ig, n, nb, me, np, pshape)

!  ig  global index of the x dimension of array element
!  n   dimension of the global array
!  nb  dimension of the block the global array is split into.
!  np  number of processors onto which the array is distributed
!
!  This function return the index of the element in the local block
!
!  END manual
!=----------------------------------------------------------------------------=!

           INTEGER ig, n, np, pshape, nb, me, q, r

           IF( pshape .EQ. BLOCK_PARTITION_SHAPE ) THEN

             q = INT(n/np)
             r = MOD(n,np)
             IF( me < r ) THEN
               LOCALINDEX_SHAPE = ig - (q+1) * me
             ELSE
               LOCALINDEX_SHAPE = ig - (q+1) * r - q * (me - r)
             END IF

           ELSE IF ( pshape .EQ. BLOCK_CYCLIC_SHAPE ) THEN

             LOCALINDEX_SHAPE = NB*((IG-1)/(NB*NP))+MOD(IG-1,NB)+1

           ELSE IF ( pshape .EQ. CYCLIC_SHAPE ) THEN
 
             LOCALINDEX_SHAPE = (ig-1)/np + 1

           ELSE
             LOCALINDEX_SHAPE = ig
           END IF
           RETURN
         END FUNCTION

!=----------------------------------------------------------------------------=!
!  BEGIN manual

          INTEGER FUNCTION localindex_desc(ig, desc, what )

!  END manual
!=----------------------------------------------------------------------------=!

            TYPE (descriptor) :: desc
            CHARACTER(LEN=*) :: what
            INTEGER ig, n, nb, np, pshape, me
            IF ( what(1:1) .EQ. 'R' .OR. what(1:1) .EQ. 'r' ) THEN
              NB = desc%nxblk;      N      = desc%nx; 
              np = desc%grid%npx;    pshape = desc%xshape
              me = desc%grid%mex
            ELSE IF ( what(1:1) .EQ. 'C' .OR. what(1:1) .EQ. 'c' ) THEN
              NB = desc%nyblk;   N      = desc%ny; 
              np = desc%grid%npy; pshape = desc%yshape
              me = desc%grid%mey
            ELSE IF ( what(1:1) .EQ. 'P' .OR. what(1:1) .EQ. 'p' ) THEN
              NB = desc%nzblk;    N      = desc%nz; 
              np = desc%grid%npz;  pshape = desc%zshape
              me = desc%grid%mez
            END IF
            localindex_desc = localindex_shape(ig,n,nb,me,np,pshape)
            RETURN
          END FUNCTION

!=----------------------------------------------------------------------------=!
!  BEGIN manual

         INTEGER FUNCTION  ownerof_shape(ig,n,nb,src_pe,np,pshape)
!
!  ig      global index of the x dimension of array element
!  n       dimension of the global array
!  nb      dimension of the block 
!  src_pe  index of the processor owning the first element of the array
!             at the moment meaningfull only for pshape = BLOCK_CYCLIC_SHAPE
!  np      number of processors
!
!  This function return the index of the processor owning the array element
!  whose global index is "ig"
!
!  END manual
!=----------------------------------------------------------------------------=!
 
       IMPLICIT NONE
       INTEGER ig, n, nb, np, pshape, src_pe, r, q
       IF( pshape .EQ. BLOCK_PARTITION_SHAPE ) THEN
         q = INT(n/np);  r = MOD(n,np)
         IF ( ig <= ((q+1)*r) ) THEN
           ownerof_shape = INT((ig-1)/(q+1))
         ELSE
           ownerof_shape = INT((ig-1-r*(q+1))/q)+r
         END IF
       ELSE IF( pshape .EQ. BLOCK_CYCLIC_SHAPE ) THEN
         ownerof_shape = MOD( src_pe + (ig - 1) / NB, NP )
       ELSE IF( pshape .EQ. CYCLIC_SHAPE ) THEN
         ownerof_shape = MOD( ig-1, np )
       END IF
       RETURN
       END FUNCTION


!=----------------------------------------------------------------------------=!
!  BEGIN manual

          INTEGER FUNCTION ownerof_desc(ig, desc, what )

!  END manual
!=----------------------------------------------------------------------------=!

            TYPE (descriptor) :: desc
            CHARACTER(LEN=*) :: what
            INTEGER ig, n, nb, src_pe, np, pshape
            IF ( what(1:1) .EQ. 'R' .OR. what(1:1) .EQ. 'r' ) THEN
              NB = desc%nxblk;      N      = desc%nx; 
              np = desc%grid%npx;    pshape = desc%xshape
              src_pe = desc%ipexs
            ELSE IF ( what(1:1) .EQ. 'C' .OR. what(1:1) .EQ. 'c' ) THEN
              NB = desc%nyblk;   N      = desc%ny; 
              np = desc%grid%npy; pshape = desc%yshape
              src_pe = desc%ipeys
            ELSE IF ( what(1:1) .EQ. 'P' .OR. what(1:1) .EQ. 'p' ) THEN
              NB = desc%nzblk;    N      = desc%nz; 
              np = desc%grid%npz;  pshape = desc%zshape
              src_pe = desc%ipezs
            END IF
            ownerof_desc = ownerof_shape(ig, n, nb, src_pe, np, pshape)
            RETURN
          END FUNCTION

!=----------------------------------------------------------------------------=!
!  BEGIN manual

          SUBROUTINE desc_gdims(d, nx, ny, nz )

!  END manual
!=----------------------------------------------------------------------------=!

            TYPE (descriptor), INTENT(IN) :: d
            INTEGER, INTENT(OUT) :: nx, ny, nz
            nx = d%nx
            ny = d%ny
            nz = d%nz
            RETURN
          END SUBROUTINE

!=----------------------------------------------------------------------------=!
!  BEGIN manual

          SUBROUTINE desc_ldims(d, nxl, nyl, nzl )

!  END manual
!=----------------------------------------------------------------------------=!

            TYPE (descriptor), INTENT(IN) :: d
            INTEGER, INTENT(OUT) :: nxl
            INTEGER, OPTIONAL, INTENT(OUT) :: nyl, nzl
            nxl = d%nxl
            IF( PRESENT( nyl ) ) nyl = d%nyl
            IF( PRESENT( nzl ) ) nzl = d%nzl

            RETURN
          END SUBROUTINE



      END MODULE descriptors_module
