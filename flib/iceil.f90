!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      integer function iceil(i,j)

      implicit none

      integer i,j
      real*8 a

      a = dble(i)/dble(j)
   
      iceil = CEILING(a)

      return
      end function iceil

