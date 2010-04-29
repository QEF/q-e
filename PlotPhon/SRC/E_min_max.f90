! Copyright Eyvaz Isaev, PWSCF group, 2009-2010
!
! Department of Physics, Chemistry and Biology (IFM)
! Linkoping Univrsity, Sweden
!
! Theoretical Physics Department
! Moscow State Institute of Steel and Alloys, Russia
!
! Program prepares a part of gnuplot file
! 1-defines axis labels
! 2-defines axis
!
!	implicit real*8(a-h,o-z)
	implicit none

!        real*8  kpoint(3,1000)
        integer maxlin, idiv, integs(100)
        integer i, j, k, nbnd, kpnt
        real*8 ax(3),by(3),cz(3)
	real*8 dl,kx,ky,kz
	real*8 E_min, E_max, THz, e_min_low
        character*10 dummy
        character*1 dummy1
        character*11 dummy3
        character*3 dummy2
	character*65 plotchar 
	
	real*8, allocatable :: e(:)
	real*8, allocatable :: kpoint(:,:)
!
	read(5,'(12x,i4,6X,i5,1X)') nbnd, kpnt
	allocate (e(nbnd))
	allocate (kpoint(3,kpnt))
	
	do i=1,kpnt
	read(5,*)
	read(5,*) (e(j),j=1,nbnd)
	
        if(i.eq.1) then
        E_min=e(1)
        E_max=e(1)
        endif
	
	do k=1,nbnd
        if(e(k).le.E_min) E_min=e(k)
        if(e(k).ge.E_max) E_max=e(k)
        enddo
	
	enddo
	
	e_min=dint(e_min*1.2)
	e_max=dint(e_max*1.2)
        e_min_low=e_min-33.0/2
	
	open(12,file='Freq_plot_unit')
	open(15,file='plot.GNU.tmp')

	read(12,*) THz

	write(15,*)
	write(15,'("THz=", f8.4)') THz
	write(15,*)
	write(15,'("E_min=",f8.1)') E_min/Thz
	write(15,'("E_max=",f8.1)') E_max/Thz
	write(15,'("Label_position=",f8.1)') E_min_low/Thz
	write(15,*)

        close(12)

! Writing labels
!	
	idiv=1
	maxlin=0

	open(9,file='kpts')
	open(10,file='K_points')

        read(10,*) dummy
        print*, dummy

        read(9,*)
	read(9,*)
	
	read(10,*) ax(1),by(1),cz(1)
	read(10,*) ax(2),by(2),cz(2)
	read(10,*) ax(3),by(3),cz(3)

	read(10,*) 
11	read(10,*,end=99) integs(idiv), (kpoint(j,idiv),j=1,3), dummy1
	print*, dummy1
	if(dummy1.eq.'G') dummy3="{/Symbol G}"
	dummy2="//dummy1//"
!	write(6,'(i6,2x,3f10.6)')  integs(idiv), (kpoint(j,idiv),j=1,3)
	
	read(9,*) dl, kx,ky,kz 
	
	dl=dl-0.03
	if(dummy1.eq.'G') then
	write(15,'("set label  ", 1h", a, 1h", "  at  ", f6.2, " , ", "Label_position")') dummy3, dl 
	else
	write(15,'("set label  ", 1h",  a, 1h", "  at  ", f6.2, " , ", "Label_position")') dummy1, dl 
	endif

	idiv=idiv+1
	maxlin=maxlin+1
	
	goto 11
99	continue	        

	write(15,*)
        close(9)

        open(9,file='kpts')
	open(15,file='plot.GNU.tmp')

	read(9,'(1x,I5)') maxlin
	read(9,*)
	
	do i=1,maxlin
	read(9,*) dl,kx,ky,kz
!	
        if(i.ne.1.and.i.ne.maxlin) then
	 write(15,'("set arrow nohead from  ",f8.5, ", E_min", "   to  ", f8.5,", E_max", "  lw 3")') dl,  dl 
        endif
!
	enddo
	
	close(9)
	close(10)
	close(14)

	write(15,*)

        plotchar="plot [:] [E_min:E_max] 'Frequency' u 1:($2/THz) w lines lt 1 lw 3"
        write(15,'(a)') plotchar
	write(15,*)

!	deallocate (e(nbnd))
!	deallocate (kpoint(3,kpnt))

	stop
	end
	
