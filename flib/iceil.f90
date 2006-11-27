!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if 0
! CURRENTLY UNUSED. PLEASE REMOVE ON NEXT OCCASION.
! AXEL KOHLMEYER. 2006-11-25
      integer function iceil(i,j)

      implicit none

      integer i,j
      real(8) a

      a = DBLE(i)/DBLE(j)
   
      iceil = CEILING(a)

      return
      end function iceil
#endif
