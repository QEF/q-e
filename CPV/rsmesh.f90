!
! Copyright (C) 2002 FPMD group
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
        USE parallel_types, ONLY: processors_grid, descriptor, BLOCK_PARTITION_SHAPE
        USE processors_grid_module, ONLY: grid_init, get_grid_info, calculate_grid_dims
        USE descriptors_module, ONLY: desc_init, get_local_dims
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

        PUBLIC :: real_space_grid, real_space_mesh_setup

!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!


!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE real_space_mesh_setup(t2dpegrid_inp)

!  This subroutines defines the processors grid for the real space data
!  distribution.
!  END manual
!  ----------------------------------------------


        USE mp, ONLY: mp_env
        LOGICAL, INTENT(IN) :: t2dpegrid_inp
        INTEGER :: npx, npy, npz
        INTEGER :: nproc, mpime, group
! ..
        CALL mp_env(nproc, mpime, group)
        npx = 1
        IF ( t2dpegrid_inp ) THEN
          CALL calculate_grid_dims(nproc, npy, npz) 
        ELSE
          npy = 1
          npz = nproc
        END IF
        CALL grid_init( real_space_grid, group, nproc, mpime, npx, npy, npz )
        IF ( ( npy * npz * npx ) /= nproc ) THEN
          CALL errore(" real_space_mesh_setup "," npx*npy*npz .ne. nproc ", npx*npy*npz)
        END IF
        depend = .TRUE.

        RETURN
      END SUBROUTINE real_space_mesh_setup

!=----------------------------------------------------------------------------=!
      END MODULE real_space_mesh
!=----------------------------------------------------------------------------=!
