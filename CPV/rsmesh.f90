!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Last modified: Tue Nov 20 10:50:18 CET 2001
!  by Carlo Cavazzoni
!  ----------------------------------------------

!=----------------------------------------------------------------------------=!
      MODULE real_space_mesh
!=----------------------------------------------------------------------------=!

        USE kinds
        USE parallel_types, ONLY: processors_grid, descriptor, BLOCK_PARTITION_DIST
        USE processors_grid_module, ONLY: grid_init, get_grid_info, calculate_grid_dims
        USE grid_dimensions, ONLY: nr1,  nr2,  nr3, nr1x, nr2x, nr3x
        USE grid_dimensions, ONLY: nr1l, nr2l, nr3l, nnrx
        USE smooth_grid_dimensions, ONLY: nr1s,  nr2s,  nr3s, nr1sx, nr2sx, nr3sx
        USE smooth_grid_dimensions, ONLY: nr1sl, nr2sl, nr3sl, nnrsx
        USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx, nnrbx
        USE smallbox_grid_dimensions, ONLY: nr1bl, nr2bl, nr3bl


        IMPLICIT NONE
        SAVE

        PRIVATE

        TYPE (processors_grid) :: real_space_grid  

        LOGICAL :: depend = .FALSE.

        PUBLIC :: real_space_grid, realspace_procgrid_init

!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!


!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE realspace_procgrid_init( twod_grid )

!  This subroutines defines the processors grid for the real space data
!  distribution.
!  END manual
!  ----------------------------------------------


        USE mp_global, ONLY: nproc, group, mpime
        LOGICAL, INTENT(IN), OPTIONAL :: twod_grid
        INTEGER :: npx, npy, npz
! ..
        npx = 1
        npy = 1
        npz = nproc
        IF ( PRESENT( twod_grid ) ) THEN
          IF( twod_grid ) CALL calculate_grid_dims(nproc, npy, npz) 
        END IF
        CALL grid_init( real_space_grid, group, nproc, mpime, npx, npy, npz )
        IF ( ( npy * npz * npx ) /= nproc ) THEN
          CALL errore(" real_space_mesh_setup "," npx*npy*npz .ne. nproc ", npx*npy*npz)
        END IF
        depend = .TRUE.

        RETURN
      END SUBROUTINE realspace_procgrid_init

!=----------------------------------------------------------------------------=!
      END MODULE real_space_mesh
!=----------------------------------------------------------------------------=!
