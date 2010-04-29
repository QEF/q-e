!
! Copyright (C) 2004-2009 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! GNU License
!
! Eyvaz Isaev
!
! Department of Physics, Chemistry and Biology,
! Institute of Physics, Linkoping University, Sweden
!
! Theoretical Physics Department,
! Moscow State Institute of Steel and Alloys, Russia 
!
! isaev@ifm.liu.se
! evyaz_isaev@yahoo.com
!
!
!   Integrate phonon DOS, Free_Energy, Entropy, MSD
!
program Atom_projected_properties
  !
  implicit none
  integer, parameter:: ndivx=10000
  real(kind=8) :: dos(ndivx),nu(ndivx) ,T,a1,a2,a3
  real(kind=8) :: de, emax, norm,norm_tot
  real(kind=8) :: F0, F0_tot, Ftot, Free, Ftot_sum, Free_sum, S0
  real(kind=8) :: T_start, T_end, T_delta, CV, CV_tot, q1,q2,q3,Entropy
  real(kind=8) :: x, coth, sinh, prefactor, amass
  integer :: i, ndiv
  integer :: i1, i2, i3
  real(kind=8) :: t_dos, p_dos
  character(len=80) :: filename
  character(len=2)  :: dummy
  real(kind=8) :: pdos(ndivx),gx(ndivx),gy(ndivx),gz(ndivx)
  real(kind=8) :: t_dos_A, p_dos_A, t_dos_O, p_dos_O, F0_A, F0_O, &
 &                Free_sum_A, Free_sum_O, Entropy_A, Entropy_O, &
 &                Heat_capacity_A, Heat_Capacity_O, E_int, E_int_A, E_int_O
  real(kind=8) :: rest
!  
!
  sinh(x)=(dexp(x)-dexp(-x))/2
  coth(x)=(dexp(x)+dexp(-x))/(dexp(x) - dexp(-x))
!
  !
  a1=0.5d0/13.6058d0/8065.5d0
  a2=8.617d-5/13.6058d0
  a3=1.0d0/8065.5d0/8.617d-5
!  amass=911.3684d0

    open(unit=9, file='projected.DOS')
    open(unit=12,file='Temperature')

    read(12,*) T_start, T_end, T_delta
        
  read(9,*)
  read (9,*) dummy, ndiv, emax, de
  if (ndiv.gt.ndivx) stop ' ndivx too small'
!
! do i=1,ndiv
!     read(9,*) nu(i),dos(i),pdos(i),gx(i),gy(i),gz(i)
!  enddo

        i=1
100     read(9,*,end=99) nu(i),dos(i),pdos(i),gx(i),gy(i),gz(i)
	i=i+1
        goto 100
99  	continue

    ndiv=i
    
! "Acoustic" branches
!
    do i=2,ndiv
    if(dos(i).le.1.d-9) then
    i1=i
    goto 1
    endif
    enddo
1   continue

! Optical modes
!
    do i=i1, ndiv
    if(dos(i).gt.1.d-8) then
    i2=i
    goto 2
    endif
    enddo
2   continue    

    i3=i1+3
    
! If i1 and i2 differ no more than 5% we suggest there is no band gap      
    rest=float(ndiv-i1)/ndiv*100
    if(rest.lt.5.0) then
    write(6,'("# Presumably, the phonon spectrum has no band gap")')
    endif

! Contribution from Acoustic part

    t_dos=0.0
    p_dos=0.0
    F0=0.0
!    
! "Acoustic part", Low frequencies     
!
    do i=2,i3-3,3
     t_dos = t_dos + 3*dos(i)+3*dos(i+1)+2*dos(i+2)
     p_dos = p_dos + 3*pdos(i)+3*pdos(i+1)+2*pdos(i+2)
     F0= F0 + 3*pdos(i)*a1*nu(i)+3*pdos(i+1)*a1*nu(i+1)+ &
   &           2*pdos(i+2)*a1*nu(i+2)
    enddo

    t_dos_A=t_dos*de*3/8
    p_dos_A=p_dos*de*3/8
    F0_A=F0*de*3/8
!
! Optical part
!
    if(rest.gt.5.0) then
    t_dos=0.0
    p_dos=0.0
    F0=0.0

    do i=i2,ndiv-3,3
     t_dos = t_dos + 3*dos(i)+3*dos(i+1)+2*dos(i+2)
     p_dos = p_dos + 3*pdos(i)+3*pdos(i+1)+2*pdos(i+2)
     F0= F0 + 3*pdos(i)*a1*nu(i)+3*pdos(i+1)*a1*nu(i+1)+ &
   &           2*pdos(i+2)*a1*nu(i+2)
    enddo
!
    t_dos_O=t_dos*de*3/8
    p_dos_O=p_dos*de*3/8
    F0_O=F0*de*3/8
    endif
!    
   if(rest.gt.5.0) then
    write(6, &
  &  '("#Contribution From:", / &
  &  "#", 25X,"Acoustic part", 5X, "Optical part", / &
  &  "#Integrated DOS        :", f14.8, 4x, f14.8, / &
  &  "#Integrated Partial DOS:", f14.8, 4x, f14.8, / &
  &  "#Zero vibration energy :", f14.8, 4x, f14.8)') &
  &  t_dos_A, T_dos_O, p_dos_A, p_dos_O, F0_A, F0_O

    write(6,  &
  &  '("#Contribution From:", / "#",25X,"Acoustic part", 36X, "Optical part")')
 
    write(6, &
  & '(105("#"),/"#","   T,K",5X,"F_vib", 9X, " C_v", 9x," S", &
  &   9X, "E_int",8X, "F_vib", 7X," C_v", 7X, " S",  & 
  &   10x, "E_int",/ 105("#"))')
    
!  & '(96("#"),/"#","   T,K",5X,"Vibr.energy", 3X, "Heat capacity", 2x,"Entropy", &
!  &   3X, "Int.energy",8X, "Vibr.energy", 3X,"Heat capacity", 2X, "Entropy",  & 
!  &   3x, "Int.energy",/ 96("#"))')

   else
   
    write(6, &
  &  '("#Contribution From:", / &
  &  "#", 25X,"Acoustic part" / &
  &  "#Integrated Phonon  DOS:", f14.8, / &
  &  "#Integrated Partial DOS:", f14.8, / &
  &  "#Zero vibration energy :", f14.8)') &
  &  t_dos_A, p_dos_A,  F0_A

    write(6,  &
  &  '("#Contribution From:", / "#",15X,"Acoustic part")')
 
    write(6, &
  & '(66("#"),/"#","   T,K",5X,"Vibr.energy", 3X, "Heat capacity", 2x,"Entropy", &
  &   7X, "Int.energy", / 66("#"))')
  
   endif

    do T=T_start, T_end, T_delta

    Free=0.0  
    S0=0.0
    CV=0.0
    E_int=0.0

    do i=2,i3-3,3
   
     q1=0.5*a3*nu(i)/T
     q2=0.5*a3*nu(i+1)/T
     q3=0.5*a3*nu(i+2)/T

!
    E_int = E_int + 3*nu(i)*pdos(i)*a1*coth(q1) + &
   &                3*nu(i+1)*pdos(i+1)*a1*coth(q2) + &
   &                2*nu(i+2)*pdos(i+2)*a1*coth(q3) 

    Free = Free + 3*pdos(i)*a2*T*dlog(2.*sinh(q1)) + &
   &              3*pdos(i+1)*a2*T*dlog(2.*sinh(q2)) + &
   &              2*pdos(i+2)*a2*T*dlog(2.*sinh(q3))
!
     S0 = S0 + 3*pdos(i)*(q1*coth(q1) - dlog(2*sinh(q1))) + &
   &           3*pdos(i+1)*(q2*coth(q2) - dlog(2*sinh(q2))) + &
   &           2*pdos(i+2)*(q3*coth(q3) - dlog(2*sinh(q3)))

     CV = CV + 3*pdos(i)*q1*q1/sinh(0.5d0*a3*nu(i)/T)**2 + &
   &           3*pdos(i+1)*q2*q2/sinh(q2)**2 + &
   &           2*pdos(i+2)*q3*q3/sinh(q3)**2

    enddo 

    Free_sum_A=Free*de*3/8
    Entropy_A=S0*de*3/8
    Heat_capacity_A=CV*de*3/8
    E_int_A=E_int*de*3/8

! Contribution from Optical part
!
    if(rest.gt.5.0) then 
!    
    Free=0.0
    S0=0.0
    CV=0.0
    E_int=0.0
    
    do i=i2,ndiv-3,3
!
     q1=0.5*a3*nu(i)/T
     q2=0.5*a3*nu(i+1)/T
     q3=0.5*a3*nu(i+2)/T
!

    E_int = E_int + 3*nu(i)*pdos(i)*a1*coth(q1) + &
   &                3*nu(i+1)*pdos(i+1)*a1*coth(q2) + &
   &                2*nu(i+2)*pdos(i+2)*a1*coth(q3) 

    Free =Free + 3*pdos(i)*a2*T*dlog(2.*sinh(q1)) + &
   &             3*pdos(i+1)*a2*T*dlog(2.*sinh(q2)) + &
   &             2*pdos(i+2)*a2*T*dlog(2.*sinh(q3))
!
     S0 = S0 + 3*pdos(i)*(q1*coth(q1) - dlog(2*sinh(q1))) + &
   &           3*pdos(i+1)*(q2*coth(q2) - dlog(2*sinh(q2))) + &
   &           2*pdos(i+2)*(q3*coth(q3) - dlog(2*sinh(q3)))

     CV = CV + 3*pdos(i)*q1*q1/sinh(q1)**2 + &
   &           3*pdos(i+1)*q2*q2/sinh(q2)**2 + &
   &           2*pdos(i+2)*q3*q3/sinh(q3)**2

    enddo 

    Free_sum_O=Free*de*3/8
    Entropy_O=S0*de*3/8
    Heat_capacity_O=CV*de*3/8
    E_int_O=E_int*de*3/8

    endif
!            
    if(rest.gt.5.0) then
    
    write(6,'(f8.2, 4f12.6, 1x, 4f12.6)') &
  & T, Free_sum_A,  Heat_capacity_A, Entropy_A, E_int_A, Free_sum_O, & 
  & Heat_capacity_O, Entropy_O, E_int_O

   else
   
    write(6,'(f8.2, 4f14.8)') &
  & T, Free_sum_A,  Heat_capacity_A, Entropy_A, E_int_A

    endif
   enddo

  close(9)
  close(12)
  
  stop 
end program Atom_projected_properties
!
