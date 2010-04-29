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

        function det4(pnt1)
        real*8 pnt1,x,y,z,s, s1, s2, s3, s4, det3, det4
	dimension pnt1(4,4),x(3),y(3),z(3)
	do 1 i=1,3
	x(i)=pnt1(i,2)
	y(i)=pnt1(i,3)
1	z(i)=pnt1(i,4)
	s1=det3(x,y,z)
	do 2 i=1,3
	x(i)=pnt1(i,1)
	y(i)=pnt1(i,3)
2	z(i)=pnt1(i,4)
	s2=det3(x,y,z)
	do 3 i=1,3
	x(i)=pnt1(i,1)
	y(i)=pnt1(i,2)
3	z(i)=pnt1(i,4)
	s3=det3(x,y,z)
	do 4 i=1,3
	x(i)=pnt1(i,1)
	y(i)=pnt1(i,2)
4	z(i)=pnt1(i,3)
	s4=det3(x,y,z)
	s=s1-s2+s3-s4

	if(dabs(s).le.1.d-12) then
	write(6,5)
5	format(8('*'),' WARNING: error in det4 due wrong input in ttrinp 
     *  ',8('*') /
     *  6('*'),' volume of a tetrahedron must not equal zero ',6('*'))
	stop
	endif
	det4=s
	return
 	end 
