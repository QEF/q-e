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
     !     io_global        stdout, ionode
     !

     IMPLICIT NONE
     SAVE

     PRIVATE
     PUBLIC :: smallbox_grid_init, smallbox_grid_info

   CONTAINS

     SUBROUTINE smallbox_grid_init( dfftp, dfftb )
       !
       USE fft_support, only: good_fft_dimension, good_fft_order
       USE fft_types,  only: fft_type_descriptor
       USE fft_smallbox_type,  only: fft_box_descriptor
       !
       IMPLICIT NONE
       !
       TYPE(fft_type_descriptor), INTENT(IN)    :: dfftp
       TYPE(fft_box_descriptor), INTENT(INOUT) :: dfftb
       !
       ! no default values for grid box: if nr*b=0, ignore

       IF( dfftb%nr1 > 0 .AND. dfftb%nr2 > 0 .AND. dfftb%nr3 > 0 ) THEN

          dfftb%nr1 = good_fft_order( dfftb%nr1 )
          dfftb%nr2 = good_fft_order( dfftb%nr2 )
          dfftb%nr3 = good_fft_order( dfftb%nr3 )
          dfftb%nr1x = good_fft_dimension( dfftb%nr1 )

       ELSE
 
          dfftb%nr1x = dfftb%nr1

       END IF

       dfftb%nr2x = dfftb%nr2
       dfftb%nr3x = dfftb%nr3
       dfftb%nnr  = dfftb%nr1x * dfftb%nr2x * dfftb%nr3x

       IF ( dfftb%nr1 > dfftp%nr1 .or. dfftb%nr2 > dfftp%nr2 .or. dfftb%nr3 > dfftp%nr3 ) &
          CALL errore(' smallbox_grid_init ', ' box grid larger than dense grid?',1)
       RETURN

     END SUBROUTINE smallbox_grid_init

     SUBROUTINE smallbox_grid_info( dfftb )
       !
       USE io_global,  ONLY: stdout, ionode
       USE fft_smallbox_type,  only: fft_box_descriptor
       !
       TYPE(fft_box_descriptor), INTENT(IN) :: dfftb
       !
       IF ( ionode ) THEN
         IF ( dfftb%nr1 > 0 .AND. dfftb%nr2 > 0 .AND. dfftb%nr3 > 0 ) THEN
           WRITE( stdout,*)
           WRITE( stdout,*) '  Small Box Real Mesh'
           WRITE( stdout,*) '  -------------------'
           WRITE( stdout,1000) dfftb%nr1, dfftb%nr2, dfftb%nr3, dfftb%nr1, dfftb%nr2, dfftb%nr3, 1, 1, 1
           WRITE( stdout,1010) dfftb%nr1x, dfftb%nr2x, dfftb%nr3x
           WRITE( stdout,1020) dfftb%nnr
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
