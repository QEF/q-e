!C *********************************************************************
!C k-points generation for band structure and phonon dispersion relations calculation.
!C Written by E.I. Isaev, 2003-2010
!C GNU General Public License 
!  
!  Department of Physics, Chemistry and Biology (IFM),
!  Linkoping University, Sweden
!
!  Theoretical Physics Department, 
!  Moscow State Institute of Steel and Alloys, Russia
!
!C *********************************************************************
!C *********************************************************************
!C
!C Input paremeters:
!C maxlin - maximal number of symmetric directions
!C ax(3), by(3), cz(3) - basis vectors for the reciprocal lattice 
!C                       (taken from scf.out)
!C idiv - number of dividing for a given direction
!C kpoint(3,i) - Irreducible Brillouin zone vertexes??? for a given crystal symmetry
!C               and a successive pair forms a symmetric direction 
!C 
!C Output parameters
!C kpts(3,i) - k-vectors to be used for band structure calculations
!C 
!C Output files
!C q.grid   - to be used for band sturcture and phonon dispersion 
!C            relation calculations (band.in, matdyn.in) 
!C q.points - to be used in conjunction with Bands.f file 
!C            which generates an output file for Gnuplot    
!C
!
      implicit real*8(a-h,o-z)
      real*8  kpoint(3,2000), kpts(3,2000)
      integer maxlin, integs(100)
      real*8 ax(3),by(3),cz(3)
      character*10 dummy
      character*1 dummy1

	idiv=1
	maxlin=0

        read(5,*) dummy

	read(5,*) ax(1),by(1),cz(1)
	read(5,*) ax(2),by(2),cz(2)
	read(5,*) ax(3),by(3),cz(3)

	read(5,*) 
11	read(5,*,end=99) integs(idiv), (kpoint(j,idiv),j=1,3), dummy1
	print*, dummy1
	idiv=idiv+1
	maxlin=maxlin+1
	goto 11
99	continue	

        do i=1, maxlin
	x=kpoint(1,i)
	y=kpoint(2,i)
	z=kpoint(3,i)

        x1=ax(1)*x +ax(2)*y + ax(3)*z
        x2=by(1)*x +by(2)*y + by(3)*z
        x3=cz(1)*x +cz(2)*y + cz(3)*z
	
	kpoint(1,i)=x1
	kpoint(2,i)=x2
	kpoint(3,i)=x3
        write(6,'(3f10.6)') (kpoint(j,i), j=1,3) 
	enddo
	
	open(9,file='kpts')
	
	do i=1,maxlin
	if(i.eq.1) then
	write(9,'("#", i5)') maxlin
	write(9,'("#")') 
	dl=0.0
        write(9,'(4f10.5)')  dl, kpoint(1,i),kpoint(2,i),kpoint(3,i)
	else
	dl=dl+sqrt((kpoint(1,i)-kpoint(1,i-1))**2 +  &
     &            (kpoint(2,i)-kpoint(2,i-1))**2 +  &
     &            (kpoint(3,i)-kpoint(3,i-1))**2)
        write(9,'(4f10.5)')  dl, kpoint(1,i),kpoint(2,i),kpoint(3,i)
        
	endif
        enddo

        close(9)
	
          nk = 0
	  do i=1,3
	  kpts(i,1)=kpoint(i,1)
	  enddo

	  nk1=2
 
         do il = 2,maxlin
             nkl = integs(il)
             nk = nk + nkl
            do ik = 1,nkl
            do ix = 1,3
            kpts(ix,nk1) = &
     &       kpoint(ix,il-1)+ & 
     &       (kpoint(ix,il) - kpoint(ix,il-1))*dfloat(ik)/integs(il)
         enddo

	write(6,'(3f10.6)') (kpts(i2,nk1),i2=1,3)
        nk1=nk1+1

	if(nk1.gt.500)  print*, 'You have more than 500 k-points. Be sure it is OK.'
	if(nk1.gt.1000) stop 'You have lots of k-points nk1>1000, this is unreasanoble.'

        enddo
	enddo

        open(11, file='ph.grid', status='unknown')

	nk1=nk1-1
	wgt=1.
	path=0.d0

	write(11,'(I6)') nk1 

	do i1=1,nk1
	if(i1.eq.1) then

	path=0.d0
	write(11,'(3f10.5,f6.2)') (kpts(i2,i1),i2=1,3),wgt
	    else
     	path=path + dsqrt((kpts(1,i1)-kpts(1,i1-1))**2+ &
     & 	(kpts(2,i1)-kpts(2,i1-1))**2+ &
     &  (kpts(3,i1)-kpts(3,i1-1))**2)

	write(11,'(3f10.5,f6.2)') (kpts(i2,i1),i2=1,3),wgt
	endif	
	enddo

      stop 
      end


