!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
program prog
  real, allocatable :: e(:,:)
  real, allocatable :: k(:,:)
  real, dimension(3) ::k0,a
  character(len=32):: input, output

  write(6,*) 'number of calculated bands' 
  read(5,*) nbands
  write(6,*) 'number of bands to be plotted'
  read(5,*) nbands2
  write(6,*) 'number of k-points'
  read(5,*) nk
  write(6,*) 'fermi level (eV)'
  read(5,*) ef
  write(6,*) 'nome bande.in e bande.out'
  read(5,*) input,output

  allocate(e(nk,nbands))
  allocate(k(nk,3))
  open(10,file=input, status='old')
  open(11,file=output, status='new')
  do i=1,nk
     !  read(10,*) 
     read(10,9020)(k(i,j),j=1,3)
     write(6,*)(k(i,j),j=1,3)
     !  read(10,*)
     read(10,9030) (e(i,n),n=1,nbands)
     write(13,9030) (e(i,n),n=1,nbands)
9020 format (14x,3f7.4)
9030 format ('   ',8f9.4)  
  enddo
  do j=1,nbands2
     dk=0
     do i=1,nk
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

