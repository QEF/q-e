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

     INTEGER :: nr1  = 0   ! global first dimension of the 3D grid 
     INTEGER :: nr2  = 0   ! global second  "           "
     INTEGER :: nr3  = 0   ! global third   "           "
     INTEGER :: nr1x = 0   ! global leading dimension
     INTEGER :: nr2x = 0
     INTEGER :: nr3x = 0
     INTEGER :: nr1l = 0   ! local first dimension 
     INTEGER :: nr2l = 0   ! 
     INTEGER :: nr3l = 0   !
     INTEGER :: nnrxl = 0  ! nr1x * nr2x * nr3l
     INTEGER :: nnrx  = 0  ! nr1x * nr2x * nr3

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
     INTEGER :: nr1sl = 0
     INTEGER :: nr2sl = 0
     INTEGER :: nr3sl = 0
     INTEGER :: nnrsxl = 0
     INTEGER :: nnrsx  = 0

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
     INTEGER :: nr1bl = 0
     INTEGER :: nr2bl = 0
     INTEGER :: nr3bl = 0
     INTEGER :: nnrbxl = 0
     INTEGER :: nnrbx  = 0

!=----------------------------------------------------------------------------=!
   END MODULE smallbox_grid_dimensions
!=----------------------------------------------------------------------------=!




!=----------------------------------------------------------------------------=!
   MODULE grid_subroutines
!=----------------------------------------------------------------------------=!

     ! This module contains subroutines that are related to grids 
     ! parameters

     USE kinds, ONLY: dbl

     IMPLICIT NONE
     SAVE

   CONTAINS


     SUBROUTINE realspace_grids_init( alat, a1, a2, a3, gcutd, gcuts, ng, ngs )
       !
       USE grid_dimensions, ONLY: nr1, nr2, nr3, nr1x, nr2x, nr3x
       USE smooth_grid_dimensions, ONLY: nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx
       USE smallbox_grid_dimensions, ONLY: nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx, &
           nnrbx, nr1bl, nr2bl, nr3bl
       USE fft_scalar, only: good_fft_dimension, good_fft_order
       !
       IMPLICIT NONE
       !
       REAL(dbl), INTENT(IN) :: alat
       REAL(dbl), INTENT(IN) :: a1(3), a2(3), a3(3)
       REAL(dbl), INTENT(IN) :: gcutd, gcuts
       INTEGER, INTENT(OUT) :: ng, ngs
       !
       REAL(dbl) :: qk(3) = 0.0d0

       ! ... This subroutines calculates the size of the real and reciprocal dense grids
       CALL ngnr_set( alat, a1, a2, a3, gcutd, qk, ng, nr1, nr2, nr3 )

       nr1 = good_fft_order( nr1 )
       nr2 = good_fft_order( nr2 )
       nr3 = good_fft_order( nr3 )

       nr1x  = good_fft_dimension(nr1 )
       nr2x  = nr2
       nr3x  = good_fft_dimension(nr3 )

       ! ... This subroutines calculates the size of the real and reciprocal smoth grids
       CALL ngnr_set( alat, a1, a2, a3, gcuts, qk, ngs, nr1s, nr2s, nr3s )

       nr1s = good_fft_order( nr1s )
       nr2s = good_fft_order( nr2s )
       nr3s = good_fft_order( nr3s )

       nr1sx = good_fft_dimension(nr1s)
       nr2sx = nr2s
       nr3sx = good_fft_dimension(nr3s)

       IF ( nr1s > nr1 .or. nr2s > nr2 .or. nr3s > nr3 ) &
         CALL errore(' realspace_grids_init ', ' smooth grid larger than dense grid?',1)

       IF( nr1b < 1 ) nr1b = 1
       IF( nr2b < 1 ) nr2b = 1
       IF( nr3b < 1 ) nr3b = 1

       nr1b = good_fft_order( nr1b ) ! small box is not parallelized
       nr2b = good_fft_order( nr2b ) ! small box is not parallelized
       nr3b = good_fft_order( nr3b ) ! small box is not parallelized
     
       nr1bx = good_fft_dimension( nr1b )
       nr2bx = nr2b
       nr3bx = nr3b
       nnrbx = nr1bx * nr2bx * nr3bx

       nr1bl = nr1b
       nr2bl = nr2b
       nr3bl = nr3b

       IF ( nr1b > nr1 .or. nr2b > nr2 .or. nr3b > nr3 ) THEN
         CALL errore(' realspace_grids_init ', ' box grid larger than dense grid?',1)
       END IF

       RETURN

     END SUBROUTINE realspace_grids_init

!=----------------------------------------------------------------------------=!
   END MODULE grid_subroutines
!=----------------------------------------------------------------------------=!
