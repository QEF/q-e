!
! Copyright (C) 2006, PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Eyvaz Isaev, 
! Theoretical Physics Department, Moscow State Institute of Steel and Alloys, 
! Department of Physics, Chemistry and Biology (IFM), Linkoping University, Sweden 
! Condensed Matter Theory Group, Uppsala University, Sweden
! evyaz_isaev@yahoo.com
! isaev@ifm.liu.se
! Eyvaz.Isaev@fysik.uu.se
!
! An auxiliary program to calculate atom projected phonon density of states
! In this program all quantities to integrate over the BZ are calculated 
! and collected in partial_DOS file
!
	implicit real*8(a-h,o-z)
!	
	dimension e_xr(3),e_xi(3), e_yr(3),e_yi(3),e_zr(3),e_zi(3), &
     & 	 Ank(100,25,3), Bnk(225,225),xx(225),t(3)
	complex*16 xmode(3), tt(3)
	open(unit=9, file='matdyn.modes',form='formatted')
	
	read(5,*) kpnts, nmodes
	natoms=nmodes/3
	
	lrecl=14*natoms*3*nmodes
	
	open(unit=18, file='partial_DOS',access='direct',recl=lrecl,form='unformatted')
	
	do kpt=1,kpnts
!
        read(9,*)
	read(9,*)
!	
1	read(9,'(6X,3f12.4 )', end=999) qx,qy,qz
	read(9,*)
	
	do modes=1,nmodes

	read(9,'(11X, I5, 3X,F15.6,8X,f15.6,8X)') nmode,omega_THz, omega_cm
!
	    do iatom=1,natoms
	    read(9, '(2x,f10.6,f11.6,f13.6,f11.6,f13.6,f11.6,5x)') &
     &	    	  e_xr(1), e_xi(1), e_yr(2), e_yi(2), e_zr(3), e_zi(3)

    	    xmode(1)=dcmplx(e_xr(1),e_xi(1))
	    xmode(2)=dcmplx(e_yr(2),e_yi(2))
    	    xmode(3)=dcmplx(e_zr(3),e_zi(3))
!
    	    Ank(modes,iatom,1)=dconjg(xmode(1))*xmode(1)
	    Ank(modes,iatom,2)=dconjg(xmode(2))*xmode(2)
    	    Ank(modes,iatom,3)=dconjg(xmode(3))*xmode(3)
!
            enddo
	enddo

	read(9,*)
!			
!  Write to file
	do modes=1,nmodes
!	
	ij=0
!	
	    do i=1,natoms
	    do j=1,3
	    ij=ij+1
	    Bnk(modes,ij)=Ank(modes,i,j)
	    enddo
	    enddo
!
	enddo
	write(18,rec=kpt) ((Bnk(i,j),j=1,ij),i=1,nmodes)

	enddo
999     continue
	close(9)
		
	stop 'Partial_DOS finished'
	end
	
	
