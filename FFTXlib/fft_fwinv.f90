!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------------=!
SUBROUTINE invfft_x( grid_type, f, dfft, dtgs, howmany )
  !! Compute G-space to R-space for a specific grid type
  !! 
  !! **grid_type = 'Dense'** : 
  !!   inverse fourier transform of potentials and charge density f
  !!   on the dense grid . On output, f is overwritten
  !! 
  !! **grid_type = 'Smooth'** :
  !!   inverse fourier transform of  potentials and charge density f
  !!   on the smooth grid . On output, f is overwritten
  !! 
  !! **grid_type = 'Wave'** :
  !!   inverse fourier transform of  wave functions f
  !!   on the smooth grid . On output, f is overwritten
  !! 
  !! **grid_type = 'Custom'** : 
  !!   inverse fourier transform of potentials and charge density f
  !!   on a custom grid. On output, f is overwritten
  !! 
  !! **grid_type = 'CustomWave'** :
  !!   inverse fourier transform of  wave functions f
  !!   on a custom grid. On output, f is overwritten
  !! 
  !! **dfft = FFT descriptor**, IMPORTANT NOTICE: grid is specified only by dfft.
  !!   No check is performed on the correspondence between dfft and grid_type.
  !!   grid_type is now used only to distinguish cases 'Wave' / 'CustomWave' 
  !!   from all other cases
  
  USE fft_scalar,    ONLY: cfft3d, cfft3ds
  USE fft_smallbox,  ONLY: cft_b, cft_b_omp
  USE fft_parallel,  ONLY: tg_cft3s
  USE fft_types,     ONLY: fft_type_descriptor
  USE task_groups,   ONLY: task_groups_descriptor
  USE fft_param,     ONLY: DP

  IMPLICIT NONE

  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: grid_type
  COMPLEX(DP) :: f(:)
  TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
  INTEGER, OPTIONAL, INTENT(IN) :: howmany
  INTEGER :: howmany_ = 1

  IF(PRESENT(howmany) ) THEN
     howmany_ = howmany
  END IF
  !
  IF( grid_type == 'Dense' ) THEN
     CALL start_clock( 'fft' )
  ELSE IF( grid_type == 'Smooth' ) THEN
     CALL start_clock( 'ffts' )
  ELSE IF( grid_type == 'Wave' ) THEN
     CALL start_clock('fftw')
  ELSE IF( grid_type == 'Custom' ) THEN
     CALL start_clock('fftc')
  ELSE IF( grid_type == 'CustomWave' ) THEN
     CALL start_clock('fftcw')
  ELSE 
     CALL fftx_error__( ' invfft ', ' unknown grid: '//grid_type , 1 )
  END IF

  IF( dfft%lpara ) THEN

     IF( howmany_ /= 1 ) THEN
        CALL fftx_error__( ' invfft ', ' howmany not yet implemented for parallel driver ', 1 )
     END IF
     
     IF( grid_type == 'Dense' .OR. grid_type == 'Smooth' .OR. &
          grid_type == 'Custom' ) THEN
        CALL tg_cft3s( f, dfft, 1 )
     ELSE IF( grid_type == 'Wave' .OR. grid_type == 'CustomWave' ) THEN
        CALL tg_cft3s( f, dfft, 2, dtgs )
     END IF

  ELSE

     IF( grid_type == 'Dense' .OR. grid_type == 'Smooth' .OR. &
         grid_type == 'Custom' ) THEN
        CALL cfft3d( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , 1)
     ELSE IF( grid_type == 'Wave' .OR. grid_type == 'CustomWave' ) THEN
        CALL cfft3ds( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                         dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , 1, &
                         dfft%isind, dfft%iplw )
     END IF

  END IF

  IF( grid_type == 'Dense' ) THEN
     CALL stop_clock( 'fft' )
  ELSE IF( grid_type == 'Smooth' ) THEN
     CALL stop_clock( 'ffts' )
  ELSE IF( grid_type == 'Wave' ) THEN
     CALL stop_clock('fftw')
  ELSE IF( grid_type == 'Custom' ) THEN
     CALL stop_clock('fftc')
  ELSE IF( grid_type == 'CustomWave' ) THEN
     CALL stop_clock('fftcw')
  END IF

  RETURN

END SUBROUTINE invfft_x
!=---------------------------------------------------------------------------=!
!
!=---------------------------------------------------------------------------=!
SUBROUTINE fwfft_x( grid_type, f, dfft, dtgs, howmany )
  !! Compute R-space to G-space for a specific grid type
  !! 
  !! **grid_type = 'Dense'**
  !!   forward fourier transform of potentials and charge density f
  !!   on the dense grid . On output, f is overwritten
  !! 
  !! **grid_type = 'Smooth'**
  !!   forward fourier transform of potentials and charge density f
  !!   on the smooth grid . On output, f is overwritten
  !! 
  !! **grid_type = 'Wave'**
  !!   forward fourier transform of  wave functions f
  !!   on the smooth grid . On output, f is overwritten
  !! 
  !! **grid_type = 'Custom'**
  !!   forward fourier transform of potentials and charge density f
  !!   on a custom grid . On output, f is overwritten
  !! 
  !! **grid_type = 'CustomWave'**
  !!   forward fourier transform of  wave functions
  !!   on a custom grid . On output, f is overwritten
  !! 
  !! **dfft = FFT descriptor**, IMPORTANT NOTICE: grid is specified only by dfft.
  !!   No check is performed on the correspondence between dfft and grid_type.
  !!   grid_type is now used only to distinguish cases 'Wave' / 'CustomWave' 
  !!   from all other cases
  
  USE fft_scalar,    ONLY: cfft3d, cfft3ds
  USE fft_parallel,  ONLY: tg_cft3s
  USE fft_types,     ONLY: fft_type_descriptor
  USE task_groups,   ONLY: task_groups_descriptor
  USE fft_param,     ONLY: DP

  IMPLICIT NONE

  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: grid_type
  COMPLEX(DP) :: f(:)
  TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
  INTEGER, OPTIONAL, INTENT(IN) :: howmany
  INTEGER :: howmany_ = 1

  IF(PRESENT(howmany) ) THEN
     howmany_ = howmany
  END IF

  IF( grid_type == 'Dense' ) THEN
     CALL start_clock( 'fft' )
  ELSE IF( grid_type == 'Smooth' ) THEN
     CALL start_clock( 'ffts' )
  ELSE IF( grid_type == 'Wave' ) THEN
     CALL start_clock( 'fftw' )
  ELSE IF( grid_type == 'Custom' ) THEN
     CALL start_clock('fftc')
  ELSE IF( grid_type == 'CustomWave' ) THEN
     CALL start_clock('fftcw')
  ELSE
     CALL fftx_error__( ' fwfft ', ' unknown grid: '//grid_type , 1 )
  END IF

  IF( dfft%lpara ) THEN

     IF( howmany_ /= 1 ) THEN
        CALL fftx_error__( ' fwfft ', ' howmany not yet implemented for parallel driver ', 1 )
     END IF
     
     IF( grid_type == 'Dense' .OR. grid_type == 'Smooth' .OR. &
         grid_type == 'Custom' ) THEN
        CALL tg_cft3s(f,dfft,-1)
     ELSE IF( grid_type == 'Wave' .OR. grid_type == 'CustomWave' ) THEN
        CALL tg_cft3s(f,dfft,-2, dtgs )
     END IF

  ELSE

     IF( grid_type == 'Dense' .OR. grid_type == 'Smooth' .OR. &
         grid_type == 'Custom' ) THEN
        CALL cfft3d( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , -1)
     ELSE IF( grid_type == 'Wave' .OR. grid_type == 'CustomWave' ) THEN
        CALL cfft3ds( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                         dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , -1, &
                         dfft%isind, dfft%iplw )
     END IF

  END IF

  IF( grid_type == 'Dense' ) THEN
     CALL stop_clock( 'fft' )
  ELSE IF( grid_type == 'Smooth' ) THEN
     CALL stop_clock( 'ffts' )
  ELSE IF( grid_type == 'Wave' ) THEN
     CALL stop_clock( 'fftw' )
  ELSE IF( grid_type == 'Custom' ) THEN
     CALL stop_clock('fftc')
  ELSE IF( grid_type == 'CustomWave' ) THEN
     CALL stop_clock('fftcw')
  END IF
  
  RETURN
  !
END SUBROUTINE fwfft_x
!=---------------------------------------------------------------------------=!
!
!=---------------------------------------------------------------------------=!
SUBROUTINE invfft_b( f, dfft, ia )
  !! Not-so-parallel 3d fft for box grid, implemented ONLY for sign=1
  !!
  !! ComputeG-space to R-space, $$ output = \sum_G f(G)exp(+iG*R) $$
  !! The array f (overwritten on output) is NOT distributed:
  !! a copy is present on each processor.
  !! The fft along z  is done on the entire grid.
  !! The fft along xy is done ONLY on planes that have components on the
  !! dense grid for each processor. Note that the final array will no
  !! longer be the same on all processors.
  !! 
  !! **grid_type** = 'Box' (only allowed value!)
  !! 
  !! **dfft** = fft descriptor for the box grid
  !! 
  !! **ia**   = index of the atom with a box grid. Used to find the number
  !!         of planes on this processors, contained in dfft%np3(ia)
  
  USE fft_scalar,    ONLY: cfft3d, cfft3ds
  USE fft_smallbox,  ONLY: cft_b, cft_b_omp
  USE fft_parallel,  ONLY: tg_cft3s
  USE fft_smallbox_type, ONLY: fft_box_descriptor
  USE fft_param,     ONLY: DP
  
  IMPLICIT NONE

  TYPE(fft_box_descriptor), INTENT(IN) :: dfft
! Removed the 'OPTIONAL' attribute. When present, the specific interfaces
! 'invfft_x' and 'invfft_b' cannot be disambiguated when the generic interface
! call is made. This is a violation the Fortran standard. The Cray compiler
! errors out on this, while the Intel only issues a warning. --rbw
! INTEGER, OPTIONAL, INTENT(IN) :: ia
  INTEGER, INTENT(IN) :: ia
  COMPLEX(DP) :: f(:)
  !
  INTEGER :: imin3, imax3, np3

  ! clocks called inside a parallel region do not work properly!
  ! in the future we probably need a thread safe version of the clock

!$omp master
  CALL start_clock( 'fftb' )
!$omp end master 

#if defined(__MPI) && !defined(__USE_3D_FFT)
     
  IF( dfft%np3( ia ) > 0 ) THEN

#if defined(__OPENMP)

     CALL cft_b_omp( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, &
                        dfft%imin3( ia ), dfft%imax3( ia ), 1 )

#else
     CALL cft_b( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                    dfft%nr1x,dfft%nr2x,dfft%nr3x, &
                    dfft%imin3( ia ), dfft%imax3( ia ), 1 )
#endif

  END IF

#else

#if defined(__OPENMP)
  CALL cft_b_omp( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                     dfft%nr1x,dfft%nr2x,dfft%nr3x, &
                     dfft%imin3( ia ), dfft%imax3( ia ), 1 )
#else
  CALL cfft3d( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                  dfft%nr1x,dfft%nr2x,dfft%nr3x, 1, 1)
#endif

#endif

!$omp master
  CALL stop_clock( 'fftb' )
!$omp end master

  RETURN
END SUBROUTINE invfft_b
!=---------------------------------------------------------------------------=!
