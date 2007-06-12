!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#if defined (__PARA) && !defined(__USE_3D_FFT)
!
! ... parallel case
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
  USE fft_scalar, ONLY : cft_1z, cft_2xy
  USE fft_base,   ONLY : fft_scatter
  USE kinds,      ONLY : DP
  USE mp_global,  ONLY : me_pool, nproc_pool
  USE pffts,      ONLY : ncts, ncplanes, ncp0s, nkcp, nxxs, npps, ncps
  USE sticks,     ONLY : dffts
  !
  IMPLICIT NONE
  !
  INTEGER,          INTENT(IN)    :: n1, n2, n3, nx1, nx2, nx3, sign
  COMPLEX(DP), INTENT(INOUT) :: f(nxxs)
  !
  INTEGER                       :: mc, i, j, ii, iproc, k, nppx
  INTEGER                       :: me_p
  COMPLEX(DP), ALLOCATABLE :: aux (:)
  INTEGER                       :: planes(nx1)
  !
  !
  CALL start_clock( 'cft3s' )
  !
  ALLOCATE( aux( nxxs ) )
  !
  me_p = me_pool + 1
  !
  ! ... see comments in cft3.f90 for the logic (or lack of it) of the following
  !
  IF ( nproc_pool == 1 ) THEN
     !
     nppx = nx3
     !
  ELSE
     !
     nppx = npps(me_p)
     !
  END IF
  !
  IF ( sign > 0 ) THEN
     !
     IF ( sign /= 2 ) THEN
        !
        CALL cft_1z( f, ncps(me_p), n3, nx3, sign, aux )
        !
        CALL fft_scatter( aux, nx3, nxxs, f, ncps, npps, sign )
        !
        f(:) = ( 0.D0, 0.D0 )
        !
        DO i = 1, ncts
           !
           mc = dffts%ismap(i)
           !
           DO j = 1, npps(me_p)
              !
              f(mc+(j-1)*ncplanes) = aux(j+(i-1)*nppx)
              !
           END DO
           !
        END DO
        !
        planes(:) = 1
        !
     ELSE
        !
        CALL cft_1z( f, nkcp(me_p), n3, nx3, sign, aux )
        !
        CALL fft_scatter( aux, nx3, nxxs, f, nkcp, npps, sign )
        !
        f(:)      = ( 0.D0 , 0.D0 )
        ii        = 0
        planes(:) = 0
        !
        DO iproc = 1, nproc_pool
           !
           DO i = 1, nkcp(iproc)
              !
              mc = dffts%ismap(i+ncp0s(iproc))
              !
              ii = ii + 1
              !
              k = MOD( mc - 1, nx1 ) + 1
              !
              planes(k) = 1
              !
              DO j = 1, npps (me_p)
                 !
                 f(mc+(j-1)*ncplanes) = aux(j+(ii-1)*nppx)
                 !
              END DO
              !
           END DO
           !
        END DO
        !
     END IF
     !
     CALL cft_2xy( f, npps(me_p), n1, n2, nx1, nx2, sign, planes )
     !
  ELSE
     !
     IF ( sign /= -2 ) THEN
        !
        planes(:) = 1
        !
     ELSE
        !
        planes(:)  = 0
        !
        DO iproc = 1, nproc_pool
           !
           DO i = 1, nkcp(iproc)
              !
              mc = dffts%ismap(i+ncp0s(iproc))
              !
              k = MOD( mc - 1, nx1 ) + 1
              !
              planes(k) = 1
              !
           END DO
           !
        END DO
        !
     END IF
     !
     CALL cft_2xy( f, npps(me_p), n1, n2, nx1, nx2, sign, planes )
     !
     IF ( sign /= -2 ) THEN
        !
        DO i = 1, ncts
           !
           mc = dffts%ismap(i)
           !
           DO j = 1, npps(me_p)
              !
              aux(j+(i-1)*nppx) = f(mc+(j-1)*ncplanes)
              !
           END DO
           !
        END DO
        !
        CALL fft_scatter( aux, nx3, nxxs, f, ncps, npps, sign )
        !
        CALL cft_1z( aux, ncps(me_p), n3, nx3, sign, f )
        !
     ELSE
        !
        ii = 0
        !
        DO iproc = 1, nproc_pool
           !
           DO i = 1, nkcp(iproc)
              !
              mc = dffts%ismap(i+ncp0s(iproc))
              !
              ii = ii + 1
              !
              DO j = 1, npps(me_p)
                 !
                 aux(j+(ii-1)*nppx) = f(mc+(j-1)*ncplanes)
                 !
              END DO
              !
           END DO
           !
        END DO
        !
        CALL fft_scatter( aux, nx3, nxxs, f, nkcp, npps, sign )
        !
        CALL cft_1z( aux, nkcp(me_p), n3, nx3, sign, f )
        !
     END IF
     !
  END IF
  !
  DEALLOCATE( aux )
  !
  CALL stop_clock( 'cft3s' )
  !
  RETURN
  !
END SUBROUTINE cft3s
!
#else
!
! ... serial case
!
#if defined (__HPM)
#  include "/cineca/prod/hpm/include/f_hpm.h"
#endif
!
!----------------------------------------------------------------------------
SUBROUTINE cft3s( f, n1, n2, n3, nx1, nx2, nx3, sign )
  !----------------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE fft_scalar, ONLY : cfft3ds, cfft3d  !  common scalar fft driver
  USE sticks,     ONLY : dffts            !  data structure for fft data layout
  !
  IMPLICIT NONE
  !
  INTEGER,          INTENT(IN)    :: n1, n2, n3, nx1, nx2, nx3, sign
  COMPLEX(DP), INTENT(INOUT) :: f(nx1*nx2*nx3)
  !
  !
  CALL start_clock( 'cft3s' )

  !
#if defined (__HPM)
  !
  CALL f_hpmstart( 20, 'cft3s' )
  !
#endif
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
#if defined (__HPM)
  CALL f_hpmstop( 20 )
#endif
  !
  CALL stop_clock ('cft3s')
  !
  RETURN
  !
END SUBROUTINE cft3s
!
#endif
