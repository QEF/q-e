!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------------=!
MODULE fft_interfaces

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: fwfft, invfft, fft_interpolate

  
  INTERFACE invfft
     !! invfft is the interface to both the standard fft **invfft_x**,
     !! and to the "box-grid" version **invfft_b**, used only in CP 
     !! (the latter has an additional argument)
     
     SUBROUTINE invfft_y( fft_kind, f, dfft, howmany )
       USE fft_types,  ONLY: fft_type_descriptor
       USE fft_param,  ONLY :DP
       IMPLICIT NONE
       CHARACTER(LEN=*),  INTENT(IN) :: fft_kind
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft
       INTEGER, OPTIONAL, INTENT(IN) :: howmany
       COMPLEX(DP) :: f(:)
     END SUBROUTINE invfft_y
     !
     SUBROUTINE invfft_b( f, dfft, ia )
       USE fft_smallbox_type,  ONLY: fft_box_descriptor
       USE fft_param,  ONLY :DP
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ia
       TYPE(fft_box_descriptor), INTENT(IN) :: dfft
       COMPLEX(DP) :: f(:)
     END SUBROUTINE invfft_b
#if defined(__CUDA)
     SUBROUTINE invfft_y_gpu( grid_type, f_d, dfft, howmany, stream )
       USE fft_types,  ONLY: fft_type_descriptor
       USE fft_param,  ONLY: DP
       USE cudafor
       IMPLICIT NONE
       CHARACTER(LEN=*),  INTENT(IN) :: grid_type
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft
       INTEGER, OPTIONAL, INTENT(IN) :: howmany
       INTEGER(kind = cuda_stream_kind), OPTIONAL, INTENT(IN) :: stream
       COMPLEX(DP), DEVICE :: f_d(:)
     END SUBROUTINE invfft_y_gpu
#endif
  END INTERFACE

  INTERFACE fwfft
     SUBROUTINE fwfft_y( fft_kind, f, dfft, howmany )
       USE fft_types,  ONLY: fft_type_descriptor
       USE fft_param,  ONLY :DP
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(IN) :: fft_kind
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft
       INTEGER, OPTIONAL, INTENT(IN) :: howmany
       COMPLEX(DP) :: f(:)
     END SUBROUTINE fwfft_y
#if defined(__CUDA)
     SUBROUTINE fwfft_y_gpu( grid_type, f_d, dfft, howmany, stream )
       USE fft_types,  ONLY: fft_type_descriptor
       USE fft_param,  ONLY: DP
       USE cudafor
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(IN) :: grid_type
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft
       INTEGER, OPTIONAL, INTENT(IN) :: howmany
       INTEGER(kind = cuda_stream_kind), OPTIONAL, INTENT(IN) :: stream
       COMPLEX(DP), DEVICE :: f_d(:)
     END SUBROUTINE fwfft_y_gpu
#endif
  END INTERFACE

  INTERFACE fft_interpolate
     !! fft_interpolate  is the interface to utility that fourier interpolate
     !! real/complex arrays between two grids 
     
     SUBROUTINE fft_interpolate_real( dfft_in, v_in, dfft_out, v_out )
       USE fft_param,  ONLY :DP
       USE fft_types,  ONLY: fft_type_descriptor
       IMPLICIT NONE
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft_in, dfft_out
       REAL(DP), INTENT(IN)  :: v_in(:)
       REAL(DP), INTENT(OUT) :: v_out(:)
     END SUBROUTINE fft_interpolate_real
     !
     SUBROUTINE fft_interpolate_complex( dfft_in, v_in, dfft_out, v_out )
       USE fft_param,  ONLY :DP
       USE fft_types,  ONLY: fft_type_descriptor
       IMPLICIT NONE
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft_in, dfft_out
       COMPLEX(DP), INTENT(IN)  :: v_in(:)
       COMPLEX(DP), INTENT(OUT) :: v_out(:)
     END SUBROUTINE fft_interpolate_complex
  END INTERFACE

END MODULE fft_interfaces
!=---------------------------------------------------------------------------=!
