!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
     subroutine matinv(hm, hi, dh)
!-----------------------------------------------------------------------
      implicit none
      real(kind=8) hm(3,3), hi(3,3), dh
      real(kind=8) d11, d12, d13, d22, d23, d33, d21, d31, d32
!
!
      d11=hm(2,2)*hm(3,3)-hm(2,3)*hm(3,2)
      d12=hm(2,3)*hm(3,1)-hm(2,1)*hm(3,3)
      d13=hm(2,1)*hm(3,2)-hm(3,1)*hm(2,2)
      d22=hm(1,1)*hm(3,3)-hm(1,3)*hm(3,1)
      d23=hm(3,1)*hm(1,2)-hm(1,1)*hm(3,2)
      d33=hm(1,1)*hm(2,2)-hm(1,2)*hm(2,1)
      d21=hm(3,2)*hm(1,3)-hm(1,2)*hm(3,3)
      d31=hm(1,2)*hm(2,3)-hm(2,2)*hm(1,3)
      d32=hm(1,3)*hm(2,1)-hm(1,1)*hm(2,3)
!
      dh=hm(1,1)*d11+hm(1,2)*d12+hm(1,3)*d13

! ... check for singular matrices
      IF(ABS(dh).LT.1.d-20) CALL error(' matinv ',' singular matrix ', 1)
!
      hi(1,1)=d11/dh
      hi(2,2)=d22/dh
      hi(3,3)=d33/dh
      hi(1,2)=d21/dh
      hi(1,3)=d31/dh
      hi(2,3)=d32/dh
      hi(2,1)=d12/dh
      hi(3,1)=d13/dh
      hi(3,2)=d23/dh
!
      return
      end
