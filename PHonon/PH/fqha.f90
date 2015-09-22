!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!     Calculate Free Energy F
!     Given phonon DOS, calculate F at various temperatures
!
program fqha
  !
  implicit none
  integer, parameter:: ndivx=10000
  real(8) :: dos(ndivx),nu(ndivx), T, a2,a3,Ftot,norm,F0
  real(8) :: de, de_, nu_,dos_
  integer :: i,ndiv
  character(len=256) :: filename
  !
  !
  write (*,"('File containing the dos >>> ')",advance="no")
  read(*,'(a)') filename
  open(unit=1,file=filename,status='old')
  !
  de = 0d0
  do i=1,ndivx
     ! nu(i) = frequencies (cm^{-1}), dos(i) in states/cm^{-1}
     read(1,*,end=10,err=20) nu(i),dos(i)
      if ( nu(i) < -1.d0 ) then
         stop ' wrong grid: omega < 0'
      else if ( nu(i) < 0.d0 ) then
         nu(i) = 0.d0
      end if
      if ( i > 1 ) then
         de = nu(i) - nu(i-1)
         if ( i > 2 ) then
            de_ = nu(i) - nu(i-1)
            if ( abs(de - de_) > 1.0d-4 ) stop ' wrong grid: not uniform'
         end if
      end if
      ndiv=i
  enddo
  read(1,*,end=10,err=20) nu_,dos_
  write(*,"('File read only up to line # ',i5)") ndivx
10 close(1)
  write(*,"('Read ',i5,' lines; Delta e (cm^-1) =',f10.6)") ndiv,de
  ! zero point energy : \sum (\hbar\omega/2) g(omega) d\omega
  F0 = 0.5 * de * dot_product ( dos(1:ndiv), nu(1:ndiv) )
  ! result is in cm^{-1}, bring it to Ry
  F0 = F0 / 8065.5d0 / 13.6058d0
  ! normalization check: \sum g(omega) d\omega = 3*Nat
  norm = sum (dos(1:ndiv)) * de
  write(*,"('Check: 3*Nat = ',f8.4,5x,'zero-point energy (Ry)=',f15.8)") norm,F0
  write (*,"('Output file for the Free energy >>> ')",advance="no")
  read(*,'(a)') filename
  if ( filename == ' ') then
     filename = 'fqha.out'
     write(*,"(' output to file ',a)") trim(filename)
  end if
  open(unit=1,file=filename,status='unknown')
  !
1 continue
  write (*,"('Temperature (K) >>> ')",advance="no")
  read (*,*,end=20,err=20) T
  if ( T < 0d0 ) then
     write(*,"('Incorrect T < 0, stopping')")
     go to 20
  end if
  ! this is Kb in Ry/K
  a2=8.617d-5/13.6058d0
  ! this is 1/Kb in cm^{-1}/K
  a3=1.0d0/8065.5d0/8.617d-5
  Ftot=0.0d0
  do i=1,ndiv
     if (T > 0.d0 .and. nu(i) > 0.d0) Ftot=Ftot+dos(i)*a2*T*log(1.0d0-exp(-a3*nu(i)/T))
  enddo
  Ftot=F0+Ftot*de
  write(*,"('T=',f8.2,'K,  F(T)= ',f15.8,' Ry')") T,Ftot
  write(1,*) T,Ftot
  !
  go to 1
20 close(1)
  !
  stop
end program fqha
!
