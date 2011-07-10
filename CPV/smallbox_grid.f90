!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
   MODULE smallbox_grid_dim
!=----------------------------------------------------------------------------=!

     !  Dimensions of the 3D real and reciprocal space FFT subgrids
     !  used for atomic augmentation charge density (USPP)
     !  Dependencies:
     !     fft_scalar       good_fft_dimension, good_fft_order
     !     grid_dimensions  nr1, nr2, nr3
     !     io_global        stdout, ionode
     !

     IMPLICIT NONE
     SAVE

     !  dimensions of the "small box" 3D grid (global)
     INTEGER :: nr1b  = 0, nr2b  = 0, nr3b  = 0

     !  dimensions of the arrays for the "small box" 3D grid (global)
     !  may differ from nr1b,nr2b,nr3b in order to boost performances
     INTEGER :: nr1bx = 0, nr2bx = 0, nr3bx = 0

     !  dimensions of the "small box" 3D grid (local on each processor)
     INTEGER :: nr1bl = 0, nr2bl = 0, nr3bl = 0

     ! size of the arrays allocated for the FFT, local to each processor:
     ! in parallel execution may differ from nr1bx*nr2bx*nr3bx
     INTEGER :: nnrbx  = 0

     PRIVATE
     PUBLIC :: nr1b, nr2b, nr3b, nr1bx, nr2bx, nr3bx, nnrbx
     PUBLIC :: smallbox_grid_init, smallbox_grid_info

   CONTAINS

     SUBROUTINE smallbox_grid_init( dense )
       !
       USE fft_scalar, only: good_fft_dimension, good_fft_order
       USE grid_types, only: grid_dim
       !
       IMPLICIT NONE
       !
       TYPE(grid_dim), INTENT(IN) :: dense
       !
       ! no default values for grid box: if nr*b=0, ignore

       IF( nr1b > 0 .AND. nr2b > 0 .AND. nr3b > 0 ) THEN

          nr1b = good_fft_order( nr1b )
          nr2b = good_fft_order( nr2b )
          nr3b = good_fft_order( nr3b )
          nr1bx = good_fft_dimension( nr1b )

       ELSE
 
          nr1bx = nr1b

       END IF

       nr2bx = nr2b
       nr3bx = nr3b
       nnrbx = nr1bx * nr2bx * nr3bx

       ! small box grid is not distributed

       nr1bl = nr1b
       nr2bl = nr2b
       nr3bl = nr3b

       IF ( nr1b > dense%nr1 .or. nr2b > dense%nr2 .or. nr3b > dense%nr3 ) &
          CALL errore(' smallbox_grid_init ', ' box grid larger than dense grid?',1)
       RETURN

     END SUBROUTINE smallbox_grid_init

     SUBROUTINE smallbox_grid_info( )
       !
       USE io_global, ONLY: stdout, ionode
       !
       IF ( ionode ) THEN
         IF ( nr1b > 0 .AND. nr2b > 0 .AND. nr3b > 0 ) THEN
           WRITE( stdout,*)
           WRITE( stdout,*) '  Small Box Real Mesh'
           WRITE( stdout,*) '  -------------------'
           WRITE( stdout,1000) nr1b, nr2b, nr3b, nr1bl, nr2bl, nr3bl, 1, 1, 1
           WRITE( stdout,1010) nr1bx, nr2bx, nr3bx
           WRITE( stdout,1020) nnrbx
         END IF
       END IF

1000  FORMAT(3X, &
         'Global Dimensions   Local  Dimensions   Processor Grid',/,3X, &
         '.X.   .Y.   .Z.     .X.   .Y.   .Z.     .X.   .Y.   .Z.',/, &
         3(1X,I5),2X,3(1X,I5),2X,3(1X,I5) )
1010  FORMAT(3X, 'Array leading dimensions ( nr1x, nr2x, nr3x )   = ', 3(1X,I5) )
1020  FORMAT(3X, 'Local number of cell to store the grid ( nrxx ) = ', 1X, I9 )


     END SUBROUTINE smallbox_grid_info

!=----------------------------------------------------------------------------=!
   END MODULE smallbox_grid_dim
!=----------------------------------------------------------------------------=!
