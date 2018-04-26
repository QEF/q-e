!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------------=!
SUBROUTINE invfft_y( fft_kind, f, dfft, howmany )
  !! Compute G-space to R-space for a specific grid type
  !! 
  !! **fft_kind = 'Rho'** : 
  !!   inverse fourier transform of potentials and charge density f
  !!   On output, f is overwritten
  !! 
  !! **fft_kind = 'Wave'** :
  !!   inverse fourier transform of  wave functions f
  !!   On output, f is overwritten
  !!
  !! **fft_kind = 'tgWave'** :
  !!   inverse fourier transform of  wave functions f with task group
  !!   On output, f is overwritten
  !!
  !! **dfft = FFT grid descriptor**, IMPORTANT NOTICE: grid is specified only by dfft.
  !!   No check is performed on the correspondence between dfft and fft_kind. 
  !!   from all other cases
  
  USE fft_scalar,    ONLY: cfft3d, cfft3ds
  USE fft_smallbox,  ONLY: cft_b, cft_b_omp
  USE fft_parallel,  ONLY: tg_cft3s
  USE fft_types,     ONLY: fft_type_descriptor
  USE fft_param,     ONLY: DP

  IMPLICIT NONE

  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: fft_kind
  COMPLEX(DP) :: f(:)
  INTEGER, OPTIONAL, INTENT(IN) :: howmany
  INTEGER :: howmany_ = 1
  CHARACTER(LEN=12) :: clock_label

  IF(PRESENT(howmany) ) THEN
     howmany_ = howmany
  END IF
  !
  IF( fft_kind == 'Rho' ) THEN
     clock_label = dfft%rho_clock_label
  ELSE IF( fft_kind == 'Wave' .OR. fft_kind == 'tgWave' ) THEN
     clock_label = dfft%wave_clock_label
  ELSE
     CALL fftx_error__( ' invfft ', ' unknown fft kind : '//fft_kind , 1 )
  END IF
  IF (clock_label == ' ') CALL fftx_error__( ' invfft ', ' uninitialized fft kind : '//fft_kind , 1 )

  CALL start_clock(clock_label)

  IF( dfft%lpara ) THEN

     IF( howmany_ /= 1 ) THEN
        CALL fftx_error__( ' invfft ', ' howmany not yet implemented for parallel driver ', 1 )
     END IF
     
     IF( fft_kind == 'Rho' ) THEN
        CALL tg_cft3s( f, dfft, 1 )
     ELSE IF( fft_kind == 'Wave' ) THEN
        CALL tg_cft3s( f, dfft, 2 )
     ELSE IF( fft_kind == 'tgWave' ) THEN
        CALL tg_cft3s( f, dfft, 3 )
     END IF

  ELSE

     IF( fft_kind == 'Rho' ) THEN
        CALL cfft3d( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x, dfft%nr2x, dfft%nr3x, howmany_ , 1)
     ELSE 
        CALL cfft3ds( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , 1, &
                        dfft%isind, dfft%iplw )
     END IF

  END IF

  CALL stop_clock( clock_label )

  RETURN

END SUBROUTINE invfft_y
!
!=---------------------------------------------------------------------------=!
!
SUBROUTINE fwfft_y( fft_kind, f, dfft, howmany )
  !! Compute R-space to G-space for a specific grid type
  !! 
  !! **fft_kind = 'Rho'**
  !!   forward fourier transform of potentials and charge density f
  !!   On output, f is overwritten
  !! 
  !! **fft_kind = 'Wave'**
  !!   forward fourier transform of  wave functions f
  !!   On output, f is overwritten
  !!
  !! **fft_kind = 'tgWave'**
  !!   forward fourier transform of wave functions f with task group
  !!   On output, f is overwritten
  !! 
  
  USE fft_scalar,    ONLY: cfft3d, cfft3ds
  USE fft_parallel,  ONLY: tg_cft3s
  USE fft_types,     ONLY: fft_type_descriptor
  USE fft_param,     ONLY: DP

  IMPLICIT NONE

  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: fft_kind
  COMPLEX(DP) :: f(:)
  INTEGER, OPTIONAL, INTENT(IN) :: howmany
  INTEGER :: howmany_ = 1
  CHARACTER(LEN=12) :: clock_label

  IF(PRESENT(howmany) ) THEN
     howmany_ = howmany
  END IF

  IF( fft_kind == 'Rho' ) THEN
     clock_label = dfft%rho_clock_label
  ELSE IF( fft_kind == 'Wave' .OR. fft_kind == 'tgWave' ) THEN
     clock_label = dfft%wave_clock_label
  ELSE
     CALL fftx_error__( ' fwfft ', ' unknown fft kind: '//fft_kind , 1 )
  END IF
  IF (clock_label == ' ') CALL fftx_error__( ' fwfft ', ' uninitialized fft kind : '//fft_kind , 1 )

  CALL start_clock(clock_label)

  IF( dfft%lpara ) THEN

     IF( howmany_ /= 1 ) THEN
        CALL fftx_error__( ' fwfft ', ' howmany not yet implemented for parallel driver ', 1 )
     END IF
     
     IF( fft_kind == 'Rho' ) THEN
        CALL tg_cft3s(f,dfft,-1)
     ELSE IF( fft_kind == 'Wave' ) THEN
        CALL tg_cft3s(f,dfft,-2 )
     ELSE IF( fft_kind == 'tgWave' ) THEN
        CALL tg_cft3s(f,dfft,-3 )
     END IF

  ELSE

     IF( fft_kind == 'Rho' ) THEN
        CALL cfft3d( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , -1)
     ELSE 
        CALL cfft3ds( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                         dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , -1, &
                         dfft%isind, dfft%iplw )
     END IF

  END IF

  CALL stop_clock( clock_label )
  
  RETURN
  !
END SUBROUTINE fwfft_y
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
  !! The fft along y  is done ONLY on planes that have components on the dense 
  !! dense grid for each processor. In addition the fft along x is done ONLY on
  !! the y-sections that have components on the dense grid for each processor. 
  !! Note that the final array will no longer be the same on all processors.
  !! 
  !! **fft_kind** = 'Box' (only allowed value!)
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
!  INTEGER :: imin2, imax2, np2, imin3, imax2, np2, imax3, np3

  ! clocks called inside a parallel region do not work properly!
  ! in the future we probably need a thread safe version of the clock

!$omp master
  CALL start_clock( 'fftb' )
!$omp end master 

#if defined(__MPI) && !defined(__USE_3D_FFT)
     
  IF( (dfft%np3( ia ) > 0) .AND. (dfft%np2( ia ) > 0) ) THEN

#if defined(_OPENMP)

     CALL cft_b_omp( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, &
                        dfft%imin2( ia ), dfft%imax2( ia ), &
                        dfft%imin3( ia ), dfft%imax3( ia ), 1 )

#else
     CALL cft_b( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                    dfft%nr1x,dfft%nr2x,dfft%nr3x, &
                    dfft%imin2( ia ), dfft%imax2( ia ), &
                    dfft%imin3( ia ), dfft%imax3( ia ), 1 )
#endif

  END IF

#else

#if defined(_OPENMP)
  CALL cft_b_omp( f, dfft%nr1, dfft%nr2, dfft%nr3, &
                     dfft%nr1x,dfft%nr2x,dfft%nr3x, &
                     dfft%imin2( ia ), dfft%imax2( ia ), &
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

#if defined(__CUDA)

SUBROUTINE invfft_y_gpu( fft_kind, f_d, dfft, howmany, stream )
  !! Compute G-space to R-space for a specific grid type
  !! 
  !! **fft_kind = 'Rho'** : 
  !!   inverse fourier transform of potentials and charge density f
  !!   On output, f is overwritten
  !! 
  !! **fft_kind = 'Wave'** :
  !!   inverse fourier transform of  wave functions f
  !!   On output, f is overwritten
  !!
  !! **fft_kind = 'tgWave'** :
  !!   inverse fourier transform of  wave functions f with task group
  !!   On output, f is overwritten
  !!
  !! **dfft = FFT grid descriptor**, IMPORTANT NOTICE: grid is specified only by dfft.
  !!   No check is performed on the correspondence between dfft and fft_kind. 
  !!   from all other cases
  USE cudafor
  USE fft_scalar,    ONLY: cfft3d_gpu, cfft3ds_gpu
  USE fft_smallbox,  ONLY: cft_b, cft_b_omp
  USE fft_parallel,  ONLY: tg_cft3s_gpu
  USE fft_types,     ONLY: fft_type_descriptor
  USE fft_param,     ONLY: DP

  IMPLICIT NONE

  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: fft_kind
  COMPLEX(DP), DEVICE :: f_d(:)
  INTEGER, OPTIONAL, INTENT(IN) :: howmany
  INTEGER(kind = cuda_stream_kind), OPTIONAL, INTENT(IN) :: stream
  !
  INTEGER                          :: howmany_ = 1
  INTEGER(kind = cuda_stream_kind) :: stream_  = 0

  CHARACTER(LEN=12) :: clock_label

  IF(PRESENT(howmany) ) THEN
     howmany_ = howmany
  END IF
  !
  IF( present( stream ) ) THEN
    stream_ = stream
  ELSE
    stream_ = 0
  END IF
  !
  IF( fft_kind == 'Rho' ) THEN
     clock_label = dfft%rho_clock_label
  ELSE IF( fft_kind == 'Wave' .OR. fft_kind == 'tgWave' ) THEN
     clock_label = dfft%wave_clock_label
  ELSE
     CALL fftx_error__( ' invfft ', ' unknown fft kind : '//fft_kind , 1 )
  END IF
  IF (clock_label == ' ') CALL fftx_error__( ' invfft ', ' uninitialized fft kind : '//fft_kind , 1 )

  CALL start_clock(clock_label)

  IF( dfft%lpara ) THEN

     IF( howmany_ /= 1 ) THEN
        CALL fftx_error__( ' invfft ', ' howmany not yet implemented for parallel driver ', 1 )
     END IF
     IF( stream_ /= 0 ) THEN
        CALL fftx_error__( ' invfft ', ' stream support not implemented for parallel driver ', 1 )
     END IF
     
     IF( fft_kind == 'Rho' ) THEN
        CALL tg_cft3s_gpu( f_d, dfft, 1, 1 )
     ELSE IF( fft_kind == 'Wave' ) THEN
        CALL tg_cft3s_gpu( f_d, dfft, 2, 1 )
     ELSE IF( fft_kind == 'tgWave' ) THEN
        CALL tg_cft3s_gpu( f_d, dfft, 3, 1 )
     END IF

  ELSE

     IF( fft_kind == 'Rho' ) THEN
        CALL cfft3d_gpu( f_d, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x, dfft%nr2x, dfft%nr3x, howmany_ , 1, stream_)
     ELSE 
        CALL cfft3ds_gpu( f_d, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , 1, &
                        dfft%isind, dfft%iplw, stream_ )
     END IF

  END IF

  CALL stop_clock( clock_label )

  RETURN

END SUBROUTINE invfft_y_gpu
!
!=---------------------------------------------------------------------------=!
!
SUBROUTINE fwfft_y_gpu( fft_kind, f_d, dfft, howmany, stream )
  !! Compute R-space to G-space for a specific grid type
  !! 
  !! **fft_kind = 'Rho'**
  !!   forward fourier transform of potentials and charge density f
  !!   On output, f is overwritten
  !! 
  !! **fft_kind = 'Wave'**
  !!   forward fourier transform of  wave functions f
  !!   On output, f is overwritten
  !!
  !! **fft_kind = 'tgWave'**
  !!   forward fourier transform of wave functions f with task group
  !!   On output, f is overwritten
  !! 
  USE cudafor
  USE fft_scalar,    ONLY: cfft3d_gpu, cfft3ds_gpu
  USE fft_parallel,  ONLY: tg_cft3s_gpu
  USE fft_types,     ONLY: fft_type_descriptor
  USE fft_param,     ONLY: DP

  IMPLICIT NONE

  TYPE(fft_type_descriptor), INTENT(IN) :: dfft
  CHARACTER(LEN=*), INTENT(IN) :: fft_kind
  COMPLEX(DP), DEVICE :: f_d(:)
  INTEGER, OPTIONAL, INTENT(IN) :: howmany
  INTEGER(kind = cuda_stream_kind), OPTIONAL, INTENT(IN) :: stream

  INTEGER                          :: howmany_ = 1
  INTEGER(kind = cuda_stream_kind) :: stream_  = 0

  CHARACTER(LEN=12) :: clock_label
  !
  IF(PRESENT(howmany) ) THEN
     howmany_ = howmany
  END IF
  !
  IF( present( stream ) ) THEN
    stream_ = stream
  ELSE
    stream_ = 0
  END IF
  !
  IF( fft_kind == 'Rho' ) THEN
     clock_label = dfft%rho_clock_label
  ELSE IF( fft_kind == 'Wave' .OR. fft_kind == 'tgWave' ) THEN
     clock_label = dfft%wave_clock_label
  ELSE
     CALL fftx_error__( ' fwfft ', ' unknown fft kind: '//fft_kind , 1 )
  END IF
  IF (clock_label == ' ') CALL fftx_error__( ' fwfft ', ' uninitialized fft kind : '//fft_kind , 1 )

  CALL start_clock(clock_label)

  IF( dfft%lpara ) THEN

     IF( howmany_ /= 1 ) THEN
        CALL fftx_error__( ' fwfft ', ' howmany not yet implemented for parallel driver ', 1 )
     END IF
     IF( stream_ /= 0 ) THEN
        CALL fftx_error__( ' fwfft ', ' stream support not implemented for parallel driver ', 1 )
     END IF
     
     IF( fft_kind == 'Rho' ) THEN
        CALL tg_cft3s_gpu(f_d,dfft,-1, 1)
     ELSE IF( fft_kind == 'Wave' ) THEN
        CALL tg_cft3s_gpu(f_d,dfft,-2, 1)
     ELSE IF( fft_kind == 'tgWave' ) THEN
        CALL tg_cft3s_gpu(f_d,dfft,-3, 1)
     END IF

  ELSE

     IF( fft_kind == 'Rho' ) THEN
        CALL cfft3d_gpu( f_d, dfft%nr1, dfft%nr2, dfft%nr3, &
                        dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , -1, stream_ )
     ELSE 
        CALL cfft3ds_gpu( f_d, dfft%nr1, dfft%nr2, dfft%nr3, &
                         dfft%nr1x,dfft%nr2x,dfft%nr3x, howmany_ , -1, &
                         dfft%isind, dfft%iplw, stream_ )
     END IF

  END IF

  CALL stop_clock( clock_label )
  
  RETURN
  !
END SUBROUTINE fwfft_y_gpu
!=---------------------------------------------------------------------------=!

#endif
