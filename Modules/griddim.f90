!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE grid_dimensions
!=----------------------------------------------------------------------------=!

     !  This module contains the dimensions of the 3D real and reciprocal space
     !  grid relative to the charde density and potential

     IMPLICIT NONE
     SAVE

     INTEGER :: nr1  = 0   ! local first dimension of the 3D grid 
     INTEGER :: nr2  = 0   ! local second  "           "
     INTEGER :: nr3  = 0   ! local third   "           "
     INTEGER :: nr1x = 0   ! leading dimension
     INTEGER :: nr2x = 0
     INTEGER :: nr3x = 0
     INTEGER :: nr1t = 0   ! global first dimension (sum over all proc of nr1)
     INTEGER :: nr2t = 0   ! 
     INTEGER :: nr3t = 0   !
     INTEGER :: nnrx  = 0  ! nr1x * nr2x * nr3
     INTEGER :: nnrxt = 0  ! nr1x * nr2x * nr3t

!=----------------------------------------------------------------------------=!
   END MODULE grid_dimensions
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   MODULE smooth_grid_dimensions
!=----------------------------------------------------------------------------=!

     !  This module contains the dimensions of the 3D real and reciprocal space
     !  grid relative to the smooth charge density ( see Vanderbilt Pseudopot )

     IMPLICIT NONE
     SAVE

     !  parameter description: same as above but for smooth grid

     INTEGER :: nr1s  = 0
     INTEGER :: nr2s  = 0
     INTEGER :: nr3s  = 0
     INTEGER :: nr1sx = 0
     INTEGER :: nr2sx = 0
     INTEGER :: nr3sx = 0
     INTEGER :: nr1st = 0
     INTEGER :: nr2st = 0
     INTEGER :: nr3st = 0
     INTEGER :: nnrsx  = 0
     INTEGER :: nnrsxt = 0

!=----------------------------------------------------------------------------=!
   END MODULE smooth_grid_dimensions
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   MODULE smallbox_grid_dimensions
!=----------------------------------------------------------------------------=!

     !  This module contains the dimensions of the 3D real and reciprocal space
     !  sub grid relative to the atomic augmentation charge density 
     !  ( see Vanderbilt Pseudopot )

     IMPLICIT NONE
     SAVE

     !  parameter description: same as above but for small box grid

     INTEGER :: nr1b  = 0
     INTEGER :: nr2b  = 0
     INTEGER :: nr3b  = 0
     INTEGER :: nr1bx = 0
     INTEGER :: nr2bx = 0
     INTEGER :: nr3bx = 0
     INTEGER :: nr1bt = 0
     INTEGER :: nr2bt = 0
     INTEGER :: nr3bt = 0
     INTEGER :: nnrbx  = 0
     INTEGER :: nnrbxt = 0

!=----------------------------------------------------------------------------=!
   END MODULE smallbox_grid_dimensions
!=----------------------------------------------------------------------------=!

