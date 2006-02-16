!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#if defined (__PARA)
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
  USE fft_scalar, ONLY : cft_1z, cft_2xy
  USE sticks,     ONLY : dfftp
  USE fft_base,   ONLY : fft_scatter
  USE kinds,      ONLY : DP
  USE mp_global,  ONLY : nproc_pool, me_pool
  USE pfft,       ONLY : nct, ncp, ncplane, nxx, npp
  !
  IMPLICIT NONE
  !
  INTEGER,          INTENT(IN)    :: n1, n2, n3, nx1, nx2, nx3, sign
  COMPLEX(DP), INTENT(INOUT) :: f(nxx)
  !
  INTEGER                        :: nxx_save, mc, i, j, ii, iproc, nppx
  INTEGER                        :: me_p
  COMPLEX(DP), ALLOCATABLE  :: aux(:)
  !
  !
  CALL start_clock( 'cft3' )
  !
  ALLOCATE( aux( nxx ) )
  !
  me_p = me_pool + 1
  !
  ! ... the following is needed for the parallel fft running on one processor
  ! ... for the special case nx3.ne.n3. Not an elegant solution, but simple,
  ! ... fast, better than the preceding one that did not work in some cases.
  ! ... Note that fft_scatter does nothing if nproc_pool=1. PG
  !
  IF ( nproc_pool == 1 ) THEN
     !
     nppx = nx3
     !
  ELSE
     !
     nppx = npp(me_p)
     !
  END IF
  !
  IF ( sign == 1 ) THEN
     !
     CALL cft_1z( f, ncp(me_p), n3, nx3, sign, aux )
     !
     CALL fft_scatter( aux, nx3, nxx, f, ncp, npp, sign )
     !
     f(:) = ( 0.D0, 0.D0 )
     !
     DO i = 1, nct
        !
        mc = dfftp%ismap(i)
        !
        DO j = 1, npp(me_p)
           !
           f(mc+(j-1)*ncplane) = aux(j+(i-1)*nppx)
           !
        END DO
        !
     END DO
     !
     CALL cft_2xy( f, npp(me_p), n1, n2, nx1, nx2, sign )
     !
  ELSE IF ( sign == - 1) THEN
     !
     CALL cft_2xy( f, npp(me_p), n1, n2, nx1, nx2, sign )
     !
     DO i = 1, nct
        !
        mc = dfftp%ismap(i)
        !
        DO j = 1, npp (me_p)
           !
           aux(j+(i-1)*nppx) = f(mc+(j-1)*ncplane)
           !
        END DO
        !
     END DO
     !
     CALL fft_scatter( aux, nx3, nxx, f, ncp, npp, sign )
     !
     CALL cft_1z( aux, ncp(me_p), n3, nx3, sign, f )
     !
  ELSE
     !
     CALL errore( 'cft3', 'not allowed', ABS( sign ) )
     !
  END IF
  !
  DEALLOCATE( aux )
  !
  CALL stop_clock( 'cft3' )
  !
  RETURN
  !
END SUBROUTINE cft3
!
#else
!
!----------------------------------------------------------------------------
SUBROUTINE cft3( f, n1, n2, n3, nx1, nx2, nx3, sign )
  !----------------------------------------------------------------------------
  !
  USE fft_scalar, ONLY : cfft3d
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,          INTENT(IN)    :: n1, n2, n3, nx1, nx2, nx3, sign
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
  RETURN
  !
END SUBROUTINE cft3
!
#endif
