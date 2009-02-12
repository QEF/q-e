SUBROUTINE wannier_u_matrix(U,hJ)

  use io_global, only: stdout, ionode, ionode_id
  use io_files
  use kinds, only: DP 
  use wannier_new, only: nwan, pp, wannier_occ, wannier_energy,wan_in
  
  implicit none
  integer i,j, k,c, iwan, l
  real(DP) :: U, hJ, u2(10,10)
  
  integer :: atoms(10)
  real(DP) :: rotm(10,10), unew(10,10), tmp
  
  write(stdout,'(5x,a34)') 'Generation of interaction matrix U'
  write(stdout,'(5x,a29)') '(works only for nspin=1 case)'
  write(stdout,*)
  u2 = 0.d0
  call mk_u(2,5,U,hJ,u2)
  
  !rotation from TB-LMTO basis to our new
  
  rotm = 0.d0
  c = 0
  do iwan=1, nwan
  	do j=1,wan_in(iwan,1)%ning
		if(wan_in(iwan,1)%ing(j)%l.eq.2) then
			c = c+1
			SELECT CASE(wan_in(iwan,1)%ing(j)%m)
			 	CASE(1) 
			 		rotm(c,3) = wan_in(iwan,1)%ing(j)%c
			 	CASE(2) 
			 		rotm(c,4) = wan_in(iwan,1)%ing(j)%c
			 	CASE(3) 
			 		rotm(c,2) = wan_in(iwan,1)%ing(j)%c
			 	CASE(4) 
			 		rotm(c,5) = wan_in(iwan,1)%ing(j)%c
			 	CASE(5) 
			 		rotm(c,1) = wan_in(iwan,1)%ing(j)%c
			END SELECT
		end if
  	end do
  end do
  
  if(c.gt.5) call errore('Too many interactiong atoms - cant construct U matrix',c)
  
  do i=1,5
  	do j=1,5
  		rotm(i+5,j+5) = rotm(i,j)
  	end do
  end do
    
  do i = 1,10
	do j = 1, 10
		tmp = 0.d0
		do k=1,10
			do l=1,10
		  		tmp=tmp+rotm(i,k)*u2(k,l)*rotm(j,l)
			enddo
		enddo
		unew(i,j)=tmp 
	enddo 
  enddo
 
!output  
	do i=1,c
		write(stdout,'(5x,10f5.2)') (unew(i,j),j=1,c), (unew(i,j+5),j=1,c)
	end do  
	do i=6,5+c
		write(stdout,'(5x,10f5.2)') (unew(i,j),j=1,c), (unew(i,j+5),j=1,c)
	end do  
	write(stdout,*)
	
	open(70,file='umatrix',status='unknown',form='formatted')
		do i=1,c
			write(70,'(10f5.2)') (unew(i,j),j=1,c), (unew(i,j+5),j=1,c)
		end do  
		do i=6,5+c
			write(70,'(10f5.2)') (unew(i,j),j=1,c), (unew(i,j+5),j=1,c)
		end do  
		write(70,*)	
	close(70)


END SUBROUTINE wannier_u_matrix
