!
! Copyright (C) 2006 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! GNU License
!
! Eyvaz Isaev

! Theoretical Physics Department,
! Moscow State Institute of Steel and Alloys
! (Technological University)
!
! Condensed Matter Theory Group,
! Uppsala University, Sweden
!
! Eyvaz.Isaev@fysik.uu.se, eyvaz_isaev@yahoo.com
!
	real*8 function det3(x,y,z)
	real*8 x,y,z, s1, s2, s3
	dimension x(3),y(3),z(3)
	s1=x(1)*(y(2)*z(3)-y(3)*z(2))
	s2=y(1)*(x(2)*z(3)-x(3)*z(2))
	s3=z(1)*(x(2)*y(3)-x(3)*y(2))
	det3=s1-s2+s3
	return
	end 
