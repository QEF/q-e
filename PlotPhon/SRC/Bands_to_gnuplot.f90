!  For Quantum Espresso 3.0 and higher
!
!  Copyright Eyvaz Isaev (2006-2007)
!  Condensed Matter Theory Group, Uppsala University, Sweden, 
!
!  Theoretical Physics Department, 
!  Moscow State Institute of Steel and Alloys, Russia
!
!  Copyright Eyvaz Isaev (2007-2010): 
!  Department of Physics, Chemistry and Biology (IFM), 
!  Linkoping University, Linkoping, Sweden  
!
!  isaev@ifm.liu.se, eyvaz_isaev@yahoo.com
!  
!  Program makes it possible plotting phonon dispersion curves along high symmetry directions
!  using Gnuplot (www.gnuplot.info).
!  Phonon frequences  are calculated by means of matdyn.x for a given number of k-vectors.
!  Then the program yields a Frequency file to be used in a gnuplot script.  

	implicit none
	integer nbnd, nks, i, j
	real*8, allocatable :: e(:,:), qx(:), qy(:), qz(:), ql(:)
!        dimension e(1000,500),qx(1000),qy(1000),qz(1000)
!        dimension ql(1000)

        open(10,file='Frequency')
        read(5,'(12X, I4, 6X,I4,2x)') nbnd, nks
        print*, 'nbands ===', nbnd

	allocate (e(nks,nbnd))
	allocate  (qx(nks),qy(nks),qz(nks),ql(nks))

	if(nbnd.gt.500) print*, 'nbnd > 500, You have too many bands!'
	if(nks.gt.1000)  stop   'nks > 1000, Too many k-points, usually 300-500 is enough.'

        do i=1,nks
        read(5,'(10X, 3f10.6)') qx(i),qy(i),qz(i)
        if(i.eq.1) then
        ql(i)=0.
        else
        ql(i)=ql(i-1)+ sqrt((qx(i)-qx(i-1))**2+(qy(i)-qy(i-1))**2+(qz(i-1)-qz(i))**2)
        endif
        read(5,'(6f10.4)') (e(i,j),j=1,nbnd)
!        write(6,'(f8.4,6f14.6)') ql(i),(e(i,j),j=1,6)
        enddo

        do i=1, nbnd
        do j=1,nks
        write(10, '(2f10.4)') ql(j), e(j,i)
        enddo
        write(10,*)
        enddo

!	deallocate (e(nks,nbnd))
!	deallocate  (qx(nks),qy(nks),qz(nks),ql(nks))

        stop
        end


