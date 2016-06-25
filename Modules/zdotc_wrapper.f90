!
! Copyright (C) 2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
! Wrapper for nonstandard implementation of complex BLAS function zdotc
! (e.g. some versions of optimized BLAS for Mac)
! In order to activate it, add -Dzdotc=zdotc_wrapper to DFLAGS in make.inc
!--------------------------------------------------------------------------
FUNCTION zdotc_wrapper(n,a,ia,b,ib) RESULT(c)
!--------------------------------------------------------------------------
   USE kinds, ONLY: dp
   IMPLICIT NONE
   COMPLEX(dp), INTENT(in) :: a(*), b(*)
   INTEGER, INTENT(in):: n,ia,ib
   COMPLEX(dp) :: c
#ifdef zdotc
#undef zdotc
   CALL zdotc(c,n,a,ia,b,ib)
#else
   c=0.0_dp
#endif
   RETURN
END FUNCTION zdotc_wrapper
