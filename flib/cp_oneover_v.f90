!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "../include/fpmd.h"

      subroutine cp_alphaover_v(n,alpha,v,alpha_over_v)

      USE kinds
      integer, intent(in) :: n
      REAL(dbl), intent(in) :: alpha 
      REAL(dbl), intent(in) :: v(*)
      REAL(dbl), intent(out) :: alpha_over_v(*)

#if defined __BENCHLIB

      call oneover_v(n,v,alpha_over_v)
      if(alpha.ne.1.0d0) then
        call FPMD_DSCAL(n,alpha,alpha_over_v,1)
      end if

#elif defined __LIBMASS

      call vrec(alpha_over_v, v, n)
      if(alpha.ne.1.0d0) then
        call FPMD_DSCAL(n,alpha,alpha_over_v,1)
      end if

#else

      do i=1,n
        alpha_over_v(i) = alpha/v(i)
      end do 

#endif

      return
      end
