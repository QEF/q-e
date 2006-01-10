!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
      MODULE processors_grid_module

        USE io_global, ONLY: stdout
        USE parallel_types, ONLY: processors_grid
        USE blacs

        IMPLICIT NONE
        SAVE

        PRIVATE

        PUBLIC :: get_blacs_grid, free_blacs_grid, grid_barrier, grid_init, &
          get_grid_dims, calculate_grid_dims, get_grid_info, get_grid_coor

        CONTAINS

          SUBROUTINE calculate_grid_dims(nproc,nprow,npcol)
! ...     This subroutine factorizes the number of processors (NPROC)
! ...     into NPROW and NPCOL,  that are the sizes of the 2D processors mesh.
            integer, intent(in) :: nproc
            integer, intent(out) :: nprow, npcol
            integer sqrtnp,i
            sqrtnp = INT( SQRT( DBLE(nproc) ) + 1 )
            DO i=1,sqrtnp
              IF(MOD(nproc,i).EQ.0) nprow = i
            END DO
            npcol = nproc/nprow
            RETURN
          END SUBROUTINE calculate_grid_dims

          SUBROUTINE get_grid_dims(gr, npx, npy, npz)
            INTEGER, INTENT(OUT) :: npx, npy, npz
            TYPE (processors_grid), INTENT(IN) :: gr
            npx = gr%npx
            npy = gr%npy
            npz = gr%npz
            RETURN
          END SUBROUTINE get_grid_dims

          SUBROUTINE get_grid_coor(gr, mex, mey, mez)
            INTEGER, INTENT(OUT) :: mex, mey, mez
            TYPE (processors_grid), INTENT(IN) :: gr
            mex = gr%mex
            mey = gr%mey
            mez = gr%mez
            RETURN
          END SUBROUTINE get_grid_coor


          SUBROUTINE get_grid_info(gr, nproc, my_pe, npx, mex, npy, mey, npz, mez)
            INTEGER, INTENT(OUT) :: nproc, my_pe
            INTEGER, INTENT(OUT) :: npx, npy, npz
            INTEGER, INTENT(OUT) :: mex, mey, mez
            TYPE (processors_grid), INTENT(IN) :: gr
            nproc = gr%nproc
            my_pe = gr%my_pe
            npx = gr%npx
            npy = gr%npy
            npz = gr%npz
            mex = gr%mex
            mey = gr%mey
            mez = gr%mez
            RETURN
          END SUBROUTINE get_grid_info



          SUBROUTINE free_blacs_grid(grid)
            TYPE (processors_grid), INTENT(OUT) :: grid
#if defined __SCALAPACK
            CALL BLACS_GRIDEXIT( grid%context )
#else
            grid%context = -1
#endif
            RETURN
          END SUBROUTINE free_blacs_grid



          SUBROUTINE get_blacs_grid(grid, debug)

            TYPE (processors_grid), INTENT(OUT) :: grid
            INTEGER, INTENT(IN), OPTIONAL :: debug
            INTEGER :: iam, nproc , nprow, npcol, context, myrow, mycol 

            INTEGER :: ndims, dims(2), coor(2)
            LOGICAL :: periods(2), reorder
            INTEGER :: comm_cart
            INTEGER :: ierr

#if defined __SCALAPACK
            CALL BLACS_PINFO( iam, nproc  )
#else
            ndims = 2

#endif


#if defined __SCALAPACK
            CALL BLACS_GET( -1, 0, context )
            CALL BLACS_GRIDINIT( context, 'R', nprow, npcol )
            CALL BLACS_GRIDINFO( context, nprow, npcol, myrow, mycol )
#else
            context = -1
            nproc  = 0
            iam = -1
            nprow = 0
            npcol = 0
            myrow = 0
            mycol = 0
#endif

            grid%context = context
            grid%nproc   = nproc 
            grid%my_pe = iam
            grid%npx  = nprow
            grid%npy  = npcol
            grid%npz = 1
            grid%mex  = myrow
            grid%mey  = mycol
            grid%mez = 0

            IF(PRESENT(debug)) THEN
              WRITE(stdout,100) iam, 'context', grid%context 
              WRITE(stdout,100) iam, 'nproc ', grid%nproc  
              WRITE(stdout,100) iam, 'my_pe', grid%my_pe 
              WRITE(stdout,100) iam, 'nprows', grid%npx 
              WRITE(stdout,100) iam, 'npcolumns', grid%npy 
              WRITE(stdout,100) iam, 'npplanes', grid%npz 
              WRITE(stdout,100) iam, 'my_row', grid%mex 
              WRITE(stdout,100) iam, 'my_column', grid%mey 
              WRITE(stdout,100) iam, 'my_plane', grid%mez 
  100         FORMAT(I4,A,I4)
            END IF

            RETURN
          END SUBROUTINE get_blacs_grid

          SUBROUTINE grid_init(grid, context, nproc , iam, &
            rows, columns, planes, my_row, my_column, my_plane, debug)
            TYPE (processors_grid), INTENT(OUT) :: grid
            INTEGER, INTENT(IN) :: rows, columns, planes
            INTEGER, INTENT(IN), OPTIONAL :: my_row, my_column,  my_plane
            INTEGER, INTENT(IN), OPTIONAL :: debug
            INTEGER, INTENT(IN) :: context, nproc , iam
            LOGICAL :: tand, tor

            IF((rows * columns * planes).NE.nproc ) THEN
              WRITE(stdout,10) rows * columns * planes, nproc  
 10           FORMAT('## grid_init Inconsistent processor grid : ',2I4) 
              STOP 'grid_init' 
            END IF

            tand = PRESENT(my_row).AND.PRESENT(my_column).AND.PRESENT(my_plane)
            tor  = PRESENT(my_row).OR.PRESENT(my_column).OR.PRESENT(my_plane)

            IF(.NOT.tand .AND. tor) THEN
              WRITE(stdout,20)
 20           FORMAT('## grid_init, optional arguments my_row, my_column, ', &
                     'my_plane',/,'## must be all present or all absent')
              STOP 'grid_init' 
            END IF

            grid%context = context
            grid%nproc   = nproc 
            grid%my_pe = iam
            grid%npx  = rows
            grid%npy  = columns
            grid%npz = planes
            IF(tand) THEN
              grid%mex  = my_row
              grid%mey  = my_column
              grid%mez = my_plane
            ELSE
              grid%mex  = MOD(MOD(iam,rows * columns), rows)
              grid%mey  = MOD(iam,rows * columns) / rows
              grid%mez = iam / (rows * columns)
            END IF


            IF(PRESENT(debug)) THEN
              WRITE(stdout,100) iam, 'context', grid%context
              WRITE(stdout,100) iam, 'nproc ', grid%nproc 
              WRITE(stdout,100) iam, 'my_pe', grid%my_pe
              WRITE(stdout,100) iam, 'nprows', grid%npx
              WRITE(stdout,100) iam, 'npcolumns', grid%npy
              WRITE(stdout,100) iam, 'npplanes', grid%npz
              WRITE(stdout,100) iam, 'my_row', grid%mex
              WRITE(stdout,100) iam, 'my_column', grid%mey
              WRITE(stdout,100) iam, 'my_plane', grid%mez
  100         FORMAT(I4,A,I4)
            END IF


            RETURN
          END SUBROUTINE grid_init
 

          SUBROUTINE grid_barrier(grid, scope)
            TYPE (processors_grid), INTENT(IN) :: grid
            CHARACTER(len=*), INTENT(IN), OPTIONAL :: scope

#if defined __SCALAPACK

            IF(PRESENT(scope)) THEN
              CALL BLACS_BARRIER (grid%context, scope)
            ELSE
              CALL BLACS_BARRIER (grid%context, 'A')
            END IF

#endif

            RETURN
          END SUBROUTINE grid_barrier


      END MODULE processors_grid_module
