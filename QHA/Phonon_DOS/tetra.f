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
! Eyvaz.Isaev@fysik.uu.se, 
! eyvaz_isaev@yahoo.com
!
C	subroutine tetra
C	Program tetra to generate k-points for integration
C       over the BZ
	program tetra 
	include 'parameters.h'
	open(9,file='kpts_out',form='formatted')

	wgt=1.0
	call k_brillouin
	call generate_tetra(npnt)

	write(9,71) npnt
	do j=1,npnt
	write(9,788) (pnt(i,j),i=1,3), wgt
	enddo
788	format(3x,4f10.5)
71      format(i4)
	close(9)
	stop
	
C	return
	end
