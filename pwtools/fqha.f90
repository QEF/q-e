!
! Copyright (C) 2001-2003 PWSCF group
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
  real(kind=8) :: dos(ndivx),nu(ndivx),T,a1,a2,a3,Ftot,norm,F0
  real(kind=8) :: de, emax
  integer :: i,ndiv
  character(len=80) :: filename
  !
  a1=0.5d0/13.6058d0/8065.5d0
  a2=8.617d-5/13.6058d0
  a3=1.0d0/8065.5d0/8.617d-5
  !
  write (*,'(a,$)') ' file containing the dos >>>'
  read(*,'(a)') filename
  open(unit=1,file=filename,status='old')
  read (1,*) ndiv, emax, de
  if (ndiv.gt.ndivx) stop ' ndivx too small'
  do i=1,ndiv
     read(1,*) nu(i),dos(i)
     if (abs(nu(i)-de*(i-0.5d0)).gt.1.0d-6) stop ' wrong grid'
  enddo
  close(1)
  write (*,'(a,$)') ' file for the Free energy >>>'
  read(*,'(a)') filename
  open(unit=1,file=filename,status='new')
  !
1 continue
  write (*,'(a,$)') ' temperature >>>'
  read (*,*,err=10) T
  Ftot=0.0d0
  F0=0.0d0
  norm=0.d0
  do i=1,ndiv
     F0=F0+dos(i)*a1*nu(i)
     if (T.gt.0.d0) Ftot=Ftot+dos(i)*a2*T*log(1.0d0-exp(-a3*nu(i)/T))
     norm=norm+dos(i)
  enddo
  Ftot=(F0+Ftot)*de
  norm=norm*de
  write(1,*) T,Ftot
  write(*,*) T,Ftot
  write(*,*) norm,F0
  !
  go to 1
10 close(1)
  !
  stop 
end program fqha
!
