!  Find atomic symbol and atomic mass for current atom from the list
!	
	character*4 Atom(100), dummy, dummy1,dummy2, dummy3, dummy4
	real*8 amass(100)
	
	read(5,*) N
!        read(5,'(a,a,a,a,a,a,a,a,a,a)') (Atom(i), i=1,N)
        read(5,'(10a)') (Atom(i), i=1,N)
	read(5,*) (amass(i), i=1,N)
		
	open(1,file='atom_name')
	read(1,'(a4)') dummy
	
	dummy1='  '//dummy
	dummy2=' '//dummy
	dummy3='  '//dummy
	dummy4=' '//dummy//' '
	
	do i=1,N
	
	if(Atom(i).eq.dummy.or.Atom(i).eq.dummy1.or. &
     &     Atom(i).eq.dummy2.or.Atom(i).eq.dummy3.or.Atom(i).eq.dummy4) then
	write(6,'(f10.4)') amass(i)
        stop
	endif
	
	enddo
	
	close(1)
	
	stop 'equivalent atomic symbols are not found '
	end
