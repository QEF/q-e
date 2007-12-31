!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!
!----------------------------------------------------------------------------
SUBROUTINE cft3s( f, n1, n2, n3, nx1, nx2, nx3, sign )
  !----------------------------------------------------------------------------
  !
  ! ... sign = +-1 : parallel 3d fft for rho and for the potential
  ! ... sign = +-2 : parallel 3d fft for wavefunctions
  !
  ! ... sign = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  ! ...              fft along z using pencils        (cft_1z)
  ! ...              transpose across nodes           (fft_scatter)
  ! ...                 and reorder
  ! ...              fft along y (using planes) and x (cft_2xy)
  ! ... sign = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  ! ...              fft along x and y(using planes)  (cft_2xy)
  ! ...              transpose across nodes           (fft_scatter)
  ! ...                 and reorder
  ! ...              fft along z using pencils        (cft_1z)
  !
  ! ...  The array "planes" signals whether a fft is needed along y :
  ! ...    planes(i)=0 : column f(i,*,*) empty , don't do fft along y
  ! ...    planes(i)=1 : column f(i,*,*) filled, fft along y needed
  ! ...  "empty" = no active components are present in f(i,*,*)
  ! ...            after (sign>0) or before (sign<0) the fft on z direction
  !
  ! ...  Note that if sign=+/-1 (fft on rho and pot.) all fft's are needed
  ! ...  and all planes(i) are set to 1
  !
  USE fft_base,     ONLY : dffts
  USE kinds,        ONLY : DP
  USE fft_parallel, ONLY : tg_cft3s
  USE fft_scalar,   ONLY : cfft3ds, cfft3d  !  common scalar fft driver
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)    :: n1, n2, n3, nx1, nx2, nx3, sign

#if defined (__PARA) && !defined(__USE_3D_FFT)
!
  COMPLEX(DP), INTENT(INOUT) :: f( dffts%nnr )
  !
  ! ... call the general purpose parallel driver
  !
  CALL tg_cft3s( f, dffts, sign )
  !
#else
  !
  ! ... serial case
  !
  COMPLEX(DP), INTENT(INOUT) :: f(nx1*nx2*nx3)
  !
  !
  CALL start_clock( 'cft3s' )
  !
  ! ... sign = +-1 : complete 3d fft (for rho and for the potential)
  !
  IF ( sign == 1 ) THEN
     !
     CALL cfft3d( f, n1, n2, n3, nx1, nx2, nx3, 1 )
     !
  ELSE IF ( sign == -1 ) THEN
     !
     CALL cfft3d( f, n1, n2, n3, nx1, nx2, nx3, -1 )
     !
     ! ... sign = +-2 : if available, call the "short" fft (for psi's)
     !
  ELSE IF ( sign == 2 ) THEN
     !
#if (defined __ESSL || defined __LINUX_ESSL || defined __FFTW || defined __FFTW3 || defined __FFTMKL8) && !defined(__USE_3D_FFT)
     !
     CALL cfft3ds( f, n1, n2, n3, nx1, nx2, nx3, 1, dffts%isind, dffts%iplw )
     !
#else
     !
     CALL cfft3d( f, n1, n2, n3, nx1, nx2, nx3, 1 )
     !
#endif
     !
  ELSE IF ( sign == -2 ) THEN
     !
#if (defined __ESSL || defined __LINUX_ESSL || defined __FFTW || defined __FFTW3 || defined __FFTMKL8) && !defined(__USE_3D_FFT)
     !
     CALL cfft3ds( f, n1, n2, n3, nx1, nx2, nx3, -1, dffts%isind, dffts%iplw )
     !
#else
     !
     CALL cfft3d( f, n1, n2, n3, nx1, nx2, nx3, -1 )
     !
#endif
     !
  ELSE
     !
     CALL errore( 'cft3', 'what should i do?', 1 )
     !
  END IF
  !
  CALL stop_clock ('cft3s')
  !
#endif
  !
  RETURN
  !
END SUBROUTINE cft3s
