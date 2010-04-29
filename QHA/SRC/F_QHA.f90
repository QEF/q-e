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
! Given phonon DOS, calculates: 
! 
! F0    - zero vibration energy, ZVE (Ry/cell)
! E_int - the internal energy   (Ry/cell)
! F_vib - the phonon contribution to the free energy (Ry/cell) including the ZVE
! C_v   - the heat capacity at a constant volume (in units of R, the universal gas constant)
! S     - entropy (in units of R, the universal gas constant )
! u^2   - mean square displacement, in Ang^2
!
! Parameters used 
! T in K, 
! F_vib in Ry/cell, 
! C_v in R (the universal gas constant, so that C_v(\infty) \to 3NR)
!        where N is the number of atoms considered, 
! S in R
! Zero vibration energy in Ry/cell, 
! phonon DOS has the norm of 3*N
!
! from R to J/(mol K) multiply C_v by 8.314 
! from R to cal/(mol K) multiply C_v by 1.985   
! from Web: R equal to 8.314 joules per Kelvin or 1.985 calories 
!
! PHDOS.out - name of the phonon DOS file obtained by means of either matdyn.x
!             applying dos=.true. option (then copy your phonon DOS file to PHDOS.out)
!             or the program with projected phonon DOS.
!           
! Avoid using another DOS filename otherwise the program will try to read projected DOS.
!
! Output file is suitable for visualization using Gnuplot (www.gnuplot.info) or Xmgace (Xmgr).
!
! MSD calculations moved to MSD.f90 
! Keep in mind: preliminary u^2 estimations are for simple metals or semiconductors (one atom type only!)
! If Theta_D is unknown, then set up Theta_D=0 or negatve, so that preliminary u^2 estimations are avoided.
! High temperature estimaton for u^2 is given at room temperature
!

program fqha
  !
  implicit none
  integer, parameter:: ndivx=10000
  real(kind=8) :: dos(ndivx),nu(ndivx),a1,a2,a3
  real(kind=8) :: pdos(ndivx),gx(ndivx),gy(ndivx),gz(ndivx)
  real(kind=8) :: de, emax, norm,norm_tot
  real(kind=8) :: F0, F0_tot, Ftot, Free, Ftot_sum, Free_sum, S0, E_internal
!
  real(kind=8) :: CV, CV_tot, E_int, q1,q2,q3,Entropy
  integer      :: T_start, T_end, T_delta
  real(kind=8) :: coth, norm_par, x, T

!  real(kind=8) :: Theta_D,prefactor,amass, factor2

  integer :: i,ndiv,n
  character(len=80) :: filename, outfile
  character(len=9)  :: dos_file, dosfile
  character(len=1)  :: dummy
!
  coth(x)=(dexp(x)+dexp(-x))/(dexp(x) - dexp(-x))
!
!
! Used constants
  a1=0.5d0/13.6058d0/8065.5d0
  a2=8.617d-5/13.6058d0
  a3=1.0d0/8065.5d0/8.617d-5
  
  dos_file='PHDOS.out'
 
  read(5,'(a)') filename
  read(5,'(a)') outfile
  read(5, *)    T_start, T_end, T_delta

  dosfile=filename(1:9)
  
  open(unit=1,file=filename)
  read (1,*)
  read (1,*) dummy, ndiv, emax, de
  if (ndiv.gt.ndivx) stop ' ndivx too small'
    print*, 'ndiv from file ===', ndiv
  

  if(dosfile.eq.dos_file) then 
      i=1
1     read(1,*, end=98) nu(i),dos(i)
      i=i+1
      goto 1
!
   else 
      i=1
2     read(1,*,end=98) nu(i),dos(i),pdos(i),gx(i),gy(i),gz(i)
      i=i+1
      goto 2
!
   endif
98 continue

    close(1)
!

    ndiv = i - 1
!    
    print*, 'ndiv===', ndiv

    if(dosfile.ne.dos_file) then 
    
    norm_par=0.
    do i=2,ndiv-3,3
    norm_par = norm_par + 3*pdos(i)+3*pdos(i+1)+2*pdos(i+2)
    enddo
    norm_par=norm_par*de*3/8

    print*, 'norm_partial==', norm_par
    
    endif

    open(unit=1,file=outfile,status='unknown')

    do T=T_start,T_end,T_delta
!
    norm=0.d0
    F0=0.d0
    Free=0.d0
    Ftot=0.d0
    CV=0.d0
    S0=0.d0
    E_int=0.d0
!
    do i=2,ndiv-3,3

     norm = norm + 3*dos(i)+3*dos(i+1)+2*dos(i+2)

     F0= F0 + 3*dos(i)*a1*nu(i)+3*dos(i+1)*a1*nu(i+1)+ & 
   &           2*dos(i+2)*a1*nu(i+2)

     q1=0.5*a3*nu(i)/T
     q2=0.5*a3*nu(i+1)/T
     q3=0.5*a3*nu(i+2)/T

! The next is the phonon contribution to the free energy without ZVE 
!
!     Ftot=Ftot + 3*dos(i)*a2*T*dlog(1.d0-exp(-2*q1))   + &
!   &             3*dos(i+1)*a2*T*dlog(1.d0-exp(-2*q2)) + &      
!   &	         2*dos(i+2)*a2*T*dlog(1.d0-exp(-2*q3)) 
!    
    E_int = E_int + 3*dos(i)*a1*nu(i)*coth(q1)   + &
   &             3*dos(i+1)*a1*nu(i+1)*coth(q2) + &
   &             2*dos(i+2)*a1*nu(i+2)*coth(q3) 
!
!
    Free = Free + 3*dos(i)*a2*T*dlog(2.*sinh(q1))   + &
   &             3*dos(i+1)*a2*T*dlog(2.*sinh(q2)) + &
   &             2*dos(i+2)*a2*T*dlog(2.*sinh(q3)) 
!
!    
     CV = CV + 3*dos(i)*q1*q1/sinh(q1)**2   + &
   &           3*dos(i+1)*q2*q2/sinh(q2)**2 + &   	
   &           2*dos(i+2)*q3*q3/sinh(q3)**2 
!
! S: a2*dos(i)*[ ]
!
     S0 = S0 + 3*dos(i)*(q1*coth(q1)   - dlog(2*sinh(q1))) + &
   &           3*dos(i+1)*(q2*coth(q2) - dlog(2*sinh(q2))) + &   	
   &           2*dos(i+2)*(q3*coth(q3) - dlog(2*sinh(q3)))


   enddo 

    norm_tot=norm*de*3/8
    F0_tot=F0*de*3/8
!    Ftot_sum=Ftot*de*3/8 + F0_tot
    E_internal=E_int*de*3/8
    Free_sum=Free*de*3/8 
    CV_tot=CV*de*3/8
    Entropy=S0*de*3/8
    
    if(T.eq.T_start) then
    write(1,'("# Zero vibration energy:",f18.10,2x, "(Ry/cell)")') F0_tot 
    write(1,'("# Phonon DOS norm      :",2x, f12.6, 6x, "! 3N for check purpose, N number of atoms in the unit cell")') norm_tot
    write(1,'("# T in K, F_vib in Ry/cell, C_v in R (the universal gas constant by 3N modes), S in k_B ")')
    write(1,'("#")') 
!    write(1,'("#   T         E_internal        F_vibration           Specific heat (C_v)      Entropy             u2")')  
    write(1,'("#   T", 9X,"E_internal",8X, "F_vibration",10X, "Specific heat (C_v)",7X,"Entropy")')  
    write(1,'(108("#"))') 

!!   &    "Specific heat (C_v)",7X,"Entropy",13X, "u2")')  
!
    endif
!
    write(1,'(f8.2,2f18.10,5X,f18.10,4X,f18.10,4X,f12.6)') T, E_internal, Free_sum,CV_tot,Entropy 
    write(6,'(f8.2,5f14.8)') T,norm_tot,F0_tot,Free_sum, Entropy

    enddo
!
10 close(1)
!
  stop 
end program fqha
!
