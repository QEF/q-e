!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

      subroutine force_conv(FNEW,FOLD,DF,NA,NAX,NSP)
      
      USE kinds
      implicit none

      integer nax,nsp,na(nsp)
      REAL(dbl) fnew(3,nax,nsp)
      REAL(dbl) fold(3,nax,nsp)

      REAL(dbl) df

      integer is,ia,j,nat
      REAL(dbl) diffnorm, normold, normnew


      nat = 0
      diffnorm  = 0.0d0
      normold   = 0.0d0
      normnew   = 0.0d0
      do is = 1,nsp
        do ia = 1,na(is)
          do j = 1,3
            diffnorm = diffnorm + (fnew(j,ia,is)-fold(j,ia,is))**2
            normold = normold + fold(j,ia,is)**2
            normnew = normnew + fnew(j,ia,is)**2
          end do
          nat = nat + na(is)
        end do
      end do

      df = (normold**0.25d0 * normnew**0.25d0)
      if(df.gt.0.0d0) then
        df = sqrt(diffnorm)/ df
      else
        df = sqrt(diffnorm)/ 1.d-10
      end if

      return
      end
