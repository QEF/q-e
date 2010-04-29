! To correct NaN problem in Phonon DOS calculatons
! To avoid this problem try delta_e \sim 1 cm^{-1} in Phonon_DOS
! These NaN are replaced by (dos(i-1)+dos(i+1))/2
!
	implicit real*8(a-h,o-z)
  
!        parameter (max_kpoints=100000, max_bands=500)

    	dimension ef(10000),dos(10000),pdos(10000),gx(10000),gy(10000),gz(10000) 
        dimension iline(10000)
	
	character*3 dummy, dummy1
	character*100 line 
	character*95  line1
	character*30  line2
	character*9   dos_file_name
	character*9   phdos

	dummy='NaN'
	dummy1='nan'
	phdos='PHDOS.out'

	read(5,'(a)') line
	dos_file_name=line(1:9)
	if(dos_file_name.eq.phdos) then 
	open(9,file='PHDOS.out')
	else
	open(9,file=line)
	endif
	
	j=0
	do
	read(9,'(a)',end=99) line
	j=j+1
	enddo
99	continue	
	
	nlines=j
	
	if(nlines.gt.10000) stop 'Tooo many lines, nlines > 10000'
	
	rewind(9)

        j=0
	do i=1,nlines
	read(9,'(a)',end=100) line
	if(line(27:29).eq.dummy.or.line(27:29).eq.dummy1.or. &
     &	   line(26:28).eq.dummy.or.line(26:28).eq.dummy1.or. &
     &	   line(23:25).eq.dummy.or.line(23:25).eq.dummy1.or. &
     &	   line(17:19).eq.dummy.or.line(17:19).eq.dummy1.or. &
     &	   line(16:18).eq.dummy.or.line(16:18).eq.dummy1)  then 
	j=j+1
	iline(i)=i
	else
	iline(i)=0
	endif
	enddo

100     continue

	REWIND(9)
	
	j=0
	do i=1,nlines 
	if(i.le.2) then 
	read(9,'(a)', end=100) line
	write(6,'(a)') line
	else 
	if(iline(i).eq.0) then 
		if(dos_file_name.ne.phdos) then 
		read(9, *,end=101) ef(i),dos(i),pdos(i),gx(i),gy(i),gz(i)
		    else
		read(9, *,end=101) ef(i),dos(i)
		endif
	else
	read(9,'(f12.4)') ef(i)
	endif
	endif
	enddo
101	continue	

	rewind(9)
		
	do i=1,nlines
	if(i.le.2) then
!	
	read(9,'(a)') line
	else 
	if(iline(i).eq.0) then
	if(dos_file_name.ne.phdos) then 
	read(9,*) ef(i), dos(i), pdos(i), gx(i),gy(i),gz(i)
	write(6,'(f12.4,5(3XG14.7))') ef(i), dos(i), pdos(i), gx(i),gy(i),gz(i)
	else 
	read(9,*) ef(i), dos(i)
	write(6,'(f12.4,(3XG14.7))') ef(i), dos(i)
	endif
        else
	read(9,'(a)') line
	if(dos_file_name.ne.phdos) then
	dos(i)=(dos(i-1)+dos(i+1))/2
	pdos(i)=(pdos(i-1)+pdos(i+1))/2
	gx(i)=(gx(i-1)+gx(i+1))/2
	gy(i)=(gy(i-1)+gy(i+1))/2
	gz(i)=(gz(i-1)+gz(i+1))/2
	else
	dos(i)=(dos(i-1)+dos(i+1))/2
	write(6,'(f12.4,(3XG14.7))') ef(i), dos(i)
	endif

	endif
	endif
	enddo 
	close(9)	
	stop 
	end	 	
	
