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
! Eyvaz.Isaev@fysik.uu.se
! evyaz_isaev@yahoo.com
!
! This program calculates phonon density of states using the tetrahedra method
!	frequences are in cm^{-1} 
!       a frequency step delta_e (cm^{-1}) is read from temporary phdos.in file
! Suitable for Espresso 3.0 and  higher versions

	include 'parameters.h'
        
!  initial parameters
	parameter (nmax=10000)

c	reading frequences 

	read(5,'(12X,I4,6X,I4)') nzone, Nkpt

!!! NB nzone=3*natoms
	irec=14*nzone
	nstar=nzone
	natoms=nzone/3
	
	print*, 'natoms==', natoms
	print*,'irec====',irec
	
	open(2,file='eigenv',access='direct',form='unformatted',
     *  recl=irec)

	do 4 ki=1,nkpt
	read(5,*) 
        read(5,*) (e(l),l=1,nzone)
        write(2,rec=ki) (e(l),l=1,nzone)

	if(ki.eq.1) then
	E_min=e(1)
	E_max=e(1)
	endif
	
	do i=1,nzone
	if(e(i).le.E_min) E_min=e(i)
	if(e(i).ge.E_max) E_max=e(i)
	enddo

4	continue
12	format(8f10.4)
	close(2)
	
!	print*,'Testing eigenv'

        open(9,file='phdos.in',access='sequential')
	read(9,*) delta_e
	print*, delta_e
        read(9,*) (atom(i),i=1,natoms)
	print*,atom(1),atom(2)
        close(9)

! give some delta ~0.0001 then E_min=0 

	if(E_min.gt.0.0) then
	Emin=0.0
	else
	Emin=E_min
	write(6,'("It seems you have imaginary frequences.\
     *        Hopefully you know what you are doing")')
        endif
! We take somewhat larger frequency range
	if(E_max.le.100) then
	    Emax=E_max*(1+0.06)
	else if(E_max.gt.100.and.E_max.lt.500) then
	    Emax=E_max*(1+0.05)
	else if(E_max.gt.500.and.E_max.lt.1000) then
	    Emax=E_max*(1+0.04)
	else if(E_max.gt.1000) then
	    Emax=E_max*(1+0.03)
	endif

	nstep=aint((Emax-Emin)/delta_e) - 1

	if(nstep.gt.nmax) then 
      	write(6,
     * 	'("Be sure your phonon spectrum is correct: Emax >2500 cm-1")')
      	write(6,'("Or you do not use very small delta_e")')
	stop 
	endif

	print*,'nstep====', nstep
	
	open(17,file='eigenv',access='direct',form='unformatted',
     *  recl=irec)

	do 41 ki=1,11
	read(17,rec=ki) (e(i),i=1,3)
	write(6,121) (e(i),i=1,3)
41	continue
121	format(8f10.4)
	close(17)

	print*,'E_min=',E_min,'   E_max=', E_max
	print*,'nstep====', nstep
	
	call k_brillouin

	call generate_tetra(npnt)
!
	open(unit=8,file='PHDOS.out',access='sequential',form='formatted')
	rewind 8

	print*,'before integration:  E_min=', E_min,'   E_max=', E_max
	write(8,'("#", I4)') natoms
	write(8,'("#", I6,"  ",2f14.8)') nstep, emax, delta_e
	do 1 i=1,nstep
	e1=emin+delta_e*(i-1)
	call Integration(e1,tdos,dos,1)
	write(8,2)  e1,dos 
1	continue
2	format(2f14.8)
	close(8)
	stop
	end
