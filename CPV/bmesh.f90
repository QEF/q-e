!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE bands_mesh
!=----------------------------------------------------------------------------=!

        USE kinds
        USE parallel_types, ONLY: processors_grid

        IMPLICIT NONE
        SAVE

        PRIVATE

        TYPE (processors_grid) :: bands_grid  

        PUBLIC :: bands_grid, bands_procgrid_init

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!

   SUBROUTINE bands_procgrid_init( )

!
!  This subroutine initialize the processor grid, among which the electronic 
!  bands are distributed. Note that bands are not distributed all across the
!  execution, but only in some subroutines, and only if it is convenient to
!  do so.
!  

     USE mp_global, ONLY: mpime, nproc, group
     USE processors_grid_module, ONLY: grid_init

       INTEGER :: npx, npy, npz

       npx = nproc
       npy = 1
       npz = 1

       CALL grid_init( bands_grid, group, nproc, mpime, npx, npy, npz )

       IF ( ( npy * npz * npx ) /= nproc ) THEN
          CALL errore(" bands_mesh_setup "," npx * npy * npz .ne. nproc ",  1 )
       END IF

     RETURN

   END SUBROUTINE bands_procgrid_init

!=----------------------------------------------------------------------------=!
   END MODULE bands_mesh
!=----------------------------------------------------------------------------=!
