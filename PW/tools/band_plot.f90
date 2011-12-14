!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
program prog
  implicit none
  real, allocatable :: e(:,:)
  real, allocatable :: k(:,:)
  real, dimension(3) ::k0,a
  integer nbnd, nbnd2, nks, i, n, j
  real ef, dk
  character(len=32):: input, output
  namelist/plot/ nbnd, nks

  write(6,*) 'Number of bands to be plotted:'
  read(5,*) nbnd2
  write(6,*) 'Fermi level (eV):'
  read(5,*) ef
  write(6,*) 'Name of the bands file (produced by band.x):'
  read(5,*) input
  write(6,*) 'Name of the output file:'
  read(5,*) output

  open(10,file=input, status='old')
  open(11,file=output, status='unknown')

  read(10,plot)

  if (nbnd2.gt.nbnd) nbnd2=nbnd

  allocate(e(nks,nbnd))
  allocate(k(nks,3))

  do i=1,nks
     read(10,*)(k(i,j),j=1,3)
     write(6,9020)(k(i,j),j=1,3)

     read(10,*) (e(i,n),n=1,nbnd)
     write(6,9030) (e(i,n),n=1,nbnd)

9020 format (14x,3f7.4)
9030 format (8f9.4)  
  enddo

  do j=1,nbnd2
     dk=0.0
     do i=1,nks
        if (i.eq.1) then
           k0=k(i,:)
        endif
        a=k(i,:)-k0
        dk=dk+sqrt(dot_product(a,a))

        write(11,*)dk,e(i,j)-ef
        k0=k(i,:)
     enddo
     write(11,*)
  enddo
  stop
end program prog

