!
! Copyright (C) 2006 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Eyvaz Isaev, Copyright 2005-2008
!
! Department of Physics, Chemistry, and Biophysics (IFM), Linkoping University, Sweden
!
! Theoretical Physics Department, Moscow State Institute of Steel and Alloys, Russia
! 
! Materials Theory Group, Institute of Physics and Materials Science, Uppsala University, Sweden
!
! isaev@ifm.liu.se
! eyvaz_isaev@yahoo.com
!
! Given projected phonon DOS, calculates: 
! 
! u^2   - mean square displacement, in Ang^2, 
!         as well <(u_x^)2>, <(u_y)^2>, and <(u_z)^2>
!
! Parameters used 
! T in K, 
! phonon DOS has the norm of 3*N
!
! Output file is suitable for visualization using Gnuplot (www.gnuplot.info) or Xmgace (Xmgr).
!
! Keep in mind: preliminary u^2 estimations are for simple metals or semiconductors (one atom type only!)
! If Theta_D is unknown, then set up Theta_D=0 or negatve, so that preliminary u^2 estimations are avoided.
! High temperature estimaton for u^2 is given at room temperature
!

program mean_square_displacement
  !
  implicit none
  integer, parameter:: ndivx=10000
  real(kind=8) :: dos(ndivx),nu(ndivx),a1,a2,a3
  real(kind=8) :: pdos(ndivx),gx(ndivx),gy(ndivx),gz(ndivx)
  real(kind=8) :: de, emax, norm,norm_tot
!
  real(kind=8) :: q1,q2,q3
  integer      :: T_start, T_end, T_delta,  T300
  real(kind=8) :: x, coth, sN,s, norm1, norm_partial

  real(kind=8) :: Theta_D,prefactor,amass

  real(kind=8) :: u2, u2_x, u2_y, u2_z, T, factor2

  integer :: i,ndiv,n
  character(len=80) :: filename, outfile
  character(len=1)  :: dummy
!
  coth(x)=(dexp(x)+dexp(-x))/(dexp(x) - dexp(-x))
!
!
! Used constants
  a1=0.5d0/13.6058d0/8065.5d0
  a2=8.617d-5/13.6058d0
  a3=1.0d0/8065.5d0/8.617d-5
  
  read(5,'(a)') filename
!  read(5,'(a)') outfile
  read(5, *)    T_start, T_end, T_delta
  read(5, *)    amass

!  read(5, *)    amass, Theta_D

  open(unit=1,file=filename)
  read (1,*)
  read (1,*) dummy, ndiv, emax, de
  if (ndiv.gt.ndivx) stop ' ndiv  too large'
    print*, 'ndiv from file ===', ndiv
!
! Safe reading of input dos file  
      i=1
2     read(1,*,end=98) nu(i),dos(i),pdos(i),gx(i),gy(i),gz(i)
      i=i+1
      goto 2
!
98 continue

    close(1)
!
    ndiv = i - 1
!    
    print*, 'ndiv===', ndiv

    norm1=0.
    do i=2,ndiv-3,3
    norm1 = norm1 + 3*pdos(i)+3*pdos(i+1)+2*pdos(i+2)
    enddo
    norm_partial=norm1*de*3/8

    print*, 'norm_partial==', norm_partial
    
!    if(Theta_D.gt.0.d0) then
!
!      if(Theta_D.lt.100.d0) then 
!      write(6,'("Unlikely, the Debye temperature is lower than 100K")')
!      write(6,'("So, be careful about low and high temperature  estimations for MSD")')
!      endif
!      
! Some estimations before theoretical calculations
! This part can be helpful for estimation purposes, so, I leave this part.
!
! Low-Temperature limit, for test case
! See Reissland's textbook "The physics of phonons" 
!
! Al
!    amass=26.98
!    Theta_D=400
!
! Si
!    amass=14.01
!    Theta_D=2000
!
! Cu
!    amass=14.01
!    Theta_D=315
!
! Pd
!    amass=106.42
!    Theta_D=275
!
! Pa
!    amass=231
!    Theta_D=175
!
!    prefactor=9*1.0545**2/(4*1.66053886*1.3804)
!    u2=prefactor*100/amass/Theta_D
!
!    print*, '#  Low temperature u^2===', u2
!
! High Temperature limit (at room temperature), for test case
!
!    T300=300
!    prefactor=9*1.0545**2/(1.66053886*1.3804)
!    u2=prefactor*100*T300/amass/Theta_D**2
!
!    print*, '#  High temperature u^2(300K)===', u2
!
!    else
!    
!    write(6,'("# Preliminary Low and High temperature u^2 estimations are not calculated")')    
!    
!    endif
    
!    open(unit=1,file=outfile,status='unknown')

    open(unit=9,file='Displacements',status='unknown')

    write(9,'(56("#"))') 
    write(9,'("#   T in K, u^2, u^2_x, u^2_y, u^2_z are  in Ang^2")')
    write(9,'("#   u^2 = u^2_x + u^2_y + u^2_z ")')
    write(9,'("#   T         u^2        u^2_x       u^2_y      u^2_z")')  
    write(9,'(56("#"))')  

!

    do T=T_start,T_end,T_delta
!
    norm=0.d0
!
    u2=0    
    u2_x=0
    u2_y=0
    u2_z=0

    do i=2,ndiv-3,3

     norm = norm + 3*dos(i)+3*dos(i+1)+2*dos(i+2)

     q1=0.5*a3*nu(i)/T
     q2=0.5*a3*nu(i+1)/T
     q3=0.5*a3*nu(i+2)/T

! Mean square displacements
!
     u2 = u2 + 3*pdos(i)*coth(q1)/q1  + &
   &           3*pdos(i+1)*coth(q2)/q2 + &   	
   &           2*pdos(i+2)*coth(q3)/q3 
!!

     u2_x = u2_x + 3*gx(i)*coth(q1)/q1  + &
   &           3*gx(i+1)*coth(q2)/q2 + &   	
   &           2*gx(i+2)*coth(q3)/q3 
!
     u2_y = u2_y + 3*gy(i)*coth(q1)/q1  + &
   &           3*gy(i+1)*coth(q2)/q2 + &   	
   &           2*gy(i+2)*coth(q3)/q3 
!
     u2_z = u2_z + 3*gz(i)*coth(q1)/q1  + &
   &           3*gz(i+1)*coth(q2)/q2 + &   	
   &           2*gz(i+2)*coth(q3)/q3 

   enddo 

    norm_tot=norm*de*3/8
    
    prefactor=1.0545**2*100/(1.66053886*1.3804)/4

    factor2=1.
!    if(dabs(norm_par/3 - 1).gt.0.1) factor2=2.
!    print*, 'factor2===', factor2
    
    u2=prefactor*u2*de*3/8/amass/T/factor2
!    u2=prefactor*u2*de*3/8/amass/T

    u2_x=prefactor*u2_x*de*3/8/amass/T/factor2
    u2_y=prefactor*u2_y*de*3/8/amass/T/factor2
    u2_z=prefactor*u2_z*de*3/8/amass/T/factor2

!
    write(9,'(f8.2,4f12.6)') T, u2, u2_x,u2_y,u2_z

    enddo
!
!10 close(1)
   close(9)
!
  stop 
end program mean_square_displacement
!
