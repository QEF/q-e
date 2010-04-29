!
! Copyright (C) 2004-2008 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Eyvaz Isaev, 
! Department of Physics, Chemistry, and Biology (IFM), Linkoping University, Sweden
! Condensed Matter Theory Group, Uppsala University, Sweden 
! Theoretical Physics Department, Moscow State Institute of Steel and Alloys, Russia
! E-mail: isaev@ifm.liu.se, eyvaz_isaev@yahoo.com
!
! Debye temperature calculation at various T
!
! Some definitions of former fqha.f90 were used
!
Program Debye_Temperature
!
! As the Debye function, D(x)=1/x^3(\int_{0}^{x} y^4*exp(y)/(exp(y)-1)^2 dy, 
! is an universal function and depends on only x = \Theta/T,  
! x (and then T_D) can be found comparing calcuated C_v(T) with  the 9*N_at*D(x) with a given precision, 
! therefore, \Theta, using the given temperature.
!
!  At low temperatures (usually, 1 - 25K) the C_v=234k_B(T/T_D)^3 is used to find T_D.
!  For T > w2 (the second moment of phonon frequencies) we put T_D = \Theta_{infty}, 
!  this mostly happens around 300K, that is why calculations are restricted by T < 300K
!
!  implicit none
  implicit real*8(a-h,o-z)
  integer, parameter:: ndivx=10000
  real(kind=8) :: dos(ndivx),nu(ndivx) ,T,a1,a2,a3,a4
  real(kind=8) :: de, emax 
  real(kind=8) :: CV, CV_tot, q1,q2,q3, cv1,pi,x0
  integer :: i,ndiv, n1
  integer :: natoms
  character(len=80) :: filename
  real(kind=8) :: pdos(ndivx),gx(ndivx),gy(ndivx),gz(ndivx)
  real(kind=8) :: x(10),D(10)

  character(len=1) :: dummy
  character(len=18) :: PHDOS_file
  
  external Debye_T
!  
!
  a1=0.5d0/13.6058d0/8065.5d0   ! 0.5\hbar cm^{-1} to Ry
  a2=8.617d-5/13.6058d0         ! K_B to Ry
  a3=1.0d0/8065.5d0/8.617d-5    ! \hbar/K_B

  pi=3.141592653
  
! K_B
   a4=1.d0/13.6058d0/8065.5d0
!
        open(9, file='T_Debye.in')
        read(9,'(a)') PHDOS_file
	read(9,*)  accuracy
	read(9,*)  T_low_start, T_low, T_low_delta
	read(9,*)  T_high, T_high_delta

	if(T_low.lt.1) then
	write(6,'("Please choose T_low around 10-20K")')
	stop
	endif
	
	if(accuracy.lt.0.00005) then
	stop " Accuracy should not be better than 5d-5, this one is quite enough"
	endif 
	
	close(9)

	open(unit=10,file=PHDOS_file)
	
        read(10,*)  dummy, natoms
        read (10,*) dummy, ndiv, emax, de
        if (ndiv.gt.ndivx) stop ' ndiv too large'

	if(PHDOS_file.eq.'PHDOS.out') then

! More safe reading a file PHDOS.out
        i=1
100     read(10,*,end=99) nu(i),dos(i) 
        i=i+1
        goto 100
99      continue
!	
	else
!
! More safe reading a file projected_DOS
        i=1
101     read(10,*,end=98) nu(i),dos(i),pdos(i),gx(i),gy(i),gz(i)
        i=i+1
        goto 101
98      continue

	 endif	
!
        close(10)
	
     w1=0.             
     w2=0.
     tdos=0.d0

     do i=2,ndiv-3,3

     tdos = tdos + 3*dos(i) + 3*dos(i+1) + 2*dos(i+2)

     w1 = w1 + 3*dos(i)*nu(i) + &
   &           3*dos(i+1)*nu(i+1) + &
   &           2*dos(i+2)*nu(i+2)


     w2 = w2 + 3*dos(i)*nu(i)**2 + &
   &           3*dos(i+1)*nu(i+1)**2 + &
   &           2*dos(i+2)*nu(i+2)**2

     enddo 

     tdos=tdos*de*3/8
     
     if((tdos/(3*natoms)).ge.2) then 
       factor=2
     else
       factor=1
      endif 
     
     write(6,'("# factor ===", f8.4)') factor

     w1=w1*de*3/8/(3*natoms)/factor
     w2=w2*de*3/8/(3*natoms)/factor

     write(6,'("# The first   moment of phonon frequencies <w>  : ",2x,f12.2," [cm^-1]", f8.2," [THz]")') w1, w1/33
     write(6,'("# The second  moment of phonon frequencies <w^2>: ",2x,f12.2, " [cm^-2]", f8.2,"[THz^2]")') w2, w2/33**2
     
     w1=4*w1/3*a3
     w2=a3*sqrt(5./3*w2)

     write(6,'("# Debye temperature via the first   moment of phonon frequencies: ",2x,f8.2, "[ K]")') w1
     write(6,'("# Debye_{\infty}    via the second  moment of phonon frequencies: ",2x,f8.2, "[ K]")') w2

! T_D, low temperatures: T < T_low
     
     do T = T_low_start, T_low, T_low_delta

! The heat capacity in a standard way
!
     CV=0.

     do i=2,ndiv-3,3

     q1=0.5*a3*nu(i)/T
     q2=0.5*a3*nu(i+1)/T
     q3=0.5*a3*nu(i+2)/T

    
     CV = CV + 3*dos(i)*q1*q1/sinh(q1)**2 + &
   &           3*dos(i+1)*q2*q2/sinh(q2)**2 + &   	
   &           2*dos(i+2)*q3*q3/sinh(q3)**2 

     enddo 

! factor is due old program which gives 2 times larger DOS, and  should be removed in the next future

        CV_tot=CV*de*3/8/factor
	
	t3=dexp(dlog(CV_tot)/3)
	coefficient=dexp(dlog(234.d0*natoms)/3)
	TD=T*coefficient/t3

        write(6,'(f8.2, 4x,f8.2,4x,f12.8,4x,f12.8)')  T, TD, CV_tot, 234*natoms*dexp(dlog(T/TD)*3)
!      
       enddo

! for intermediate and higher T , T_high usually up to room temperature (300K)
     
	do T = T_low+5, T_high, T_high_delta

! The heat capacity in a standard way
!
     CV=0.

     do i=2,ndiv-3,3

     q1=0.5*a3*nu(i)/T
     q2=0.5*a3*nu(i+1)/T
     q3=0.5*a3*nu(i+2)/T

    
     CV = CV + 3*dos(i)*q1*q1/sinh(q1)**2 + &
   &           3*dos(i+1)*q2*q2/sinh(q2)**2 + &   	
   &           2*dos(i+2)*q3*q3/sinh(q3)**2 

     enddo 

! factor is due old program which gives 2 times larger DOS, and  should be removed in the next future

        CV_tot=CV*de*3/8/factor

	x0=0.00
!	
!      do while (.true.) 
	TD=0

       do while  (x0.le.500.0)
!      
       if(dabs(CV_tot-9*natoms*Debye_T(x0)).le.accuracy &
     &    .and.(dabs(CV_tot-3*natoms).gt.0.01)) then    
         TD=T*x0	
         write(6,'(f8.2, 4x,f8.2,4x,f12.8,4x,f12.8)') T, TD, CV_tot, 9*natoms*Debye_T(x0)
         if(dabs(TD-w2).lt.1) then
	 stop "# There is no need to calculate the Debye temperature at higher temperatures"
	 endif
         goto 102
!	 else
!	 goto 102
       endif	 
         x0=x0+0.0001 
!         if(dabs(CV_tot-3*natoms).le.0.001) x0=x0-0.0001
       enddo
 102     continue  

      enddo
     
  stop 
end program Debye_Temperature



