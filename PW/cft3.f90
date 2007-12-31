!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE cft3( f, n1, n2, n3, nx1, nx2, nx3, sign )
  !----------------------------------------------------------------------------
  !
  ! ...  sign = +-1 : parallel 3d fft for rho and for the potential
  !
  ! ...  sign = +1  : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  ! ...               fft along z using pencils (cft_1z)
  ! ...               transpose across nodes    (fft_scatter)
  ! ...                  and reorder
  ! ...               fft along y and x         (cft_2xy)
  !
  ! ...  sign = -1  : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  ! ...               fft along x and y         (cft_2xy)
  ! ...               transpose across nodes    (fft_scatter)
  ! ...                  and reorder
  ! ...               fft along z using pencils (cft_1z)
  !
  USE fft_scalar,   ONLY : cfft3d
  USE fft_base,     ONLY : dfftp
  USE fft_parallel, ONLY : tg_cft3s
  USE kinds,        ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)    :: n1, n2, n3, nx1, nx2, nx3, sign

#if defined (__PARA) && ! defined (__USE_3D_FFT)

  COMPLEX(DP), INTENT(INOUT) :: f( dfftp%nnr )
  !
  CALL tg_cft3s( f, dfftp, sign )

#else
!
  COMPLEX(DP), INTENT(INOUT) :: f(nx1*nx2*nx3)
  !
  !
  CALL start_clock( 'cft3' )
  !
  ! ... sign = +-1 : complete 3d fft (for rho and for the potential)
  !
  IF ( sign == 1 ) THEN
     !
     CALL cfft3d( f, n1, n2, n3, nx1, nx2, nx3, 1 )
     !
  ELSE IF ( sign == - 1 ) THEN
     !
     CALL cfft3d( f, n1, n2, n3, nx1, nx2, nx3, -1 )
     !
  ELSE
     !
     CALL errore( 'cft3', 'what should i do?', 1 )
     !
  ENDIF
  !
  CALL stop_clock( 'cft3' )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE cft3
