!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
program wannier
  !-----------------------------------------------------------------------
  !
  character (len=3) :: nodenumber
  !
  call start_postproc (nodenumber)  
  call do_wannier (nodenumber)  
  !
  call stop_pp
  stop
end program wannier

!
!-----------------------------------------------------------------------
subroutine do_wannier(nodenumber)
  !-----------------------------------------------------------------------
  !
  use pwcom  
  use io  

  implicit none
  character(len=3) :: nodenumber
  integer , dimension(3):: nk
  real(kind=8), dimension(3):: s0
  integer :: ik, i

  namelist / inputpp / tmp_dir, filpun, nk, s0
  !
  nd_nmbr = nodenumber  
  !
  !   set default values for variables in namelist
  !
  filpun = ' '  
  tmp_dir = './'
  nk = 0
  s0 = 0.d0
  !
  !     reading the namelist inputpp
  !
  read (5, inputpp)
  !
  !     Check of namelist variables
  !
  if (filpun.eq.' ') &
       call errore ('wannier', 'Missing input file name', 1)
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  call read_file  
  call openfil  

  do ik = 1, nks
    write(6,fmt="(3F9.6)" ) xk(1,ik), xk(2,ik), xk(3,ik)
  end do
  !
  if (nk(1)==0 .and. nk(2)==0 .and. nk(3) == 0) then
     nk(1) = nk1
     nk(2) = nk2
     nk(3) = nk3
     s0(1) = k1/2d0
     s0(2) = k2/2d0
     s0(3) = k3/2d0
  end if
  call write_wannier (nk, s0)

  call cryst_to_cart (nks,xk,at,-1)
  do ik = 1, nks
    write(6,fmt="(' ik = ',I3,3F10.6)" ) ik,xk(1,ik), xk(2,ik), xk(3,ik)
  end do
  write(6,*)
  do i = 1, 3
    write(6,fmt="(' a(',I1,')',3F10.6)" ) i,at(1,i), at(2,i), at(3,i)
  end do

  !
  return  
end subroutine do_wannier
!
!-----------------------------------------------------------------------
subroutine write_wannier (nk, s0)
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  use pwcom  
  use io  

  implicit none
  integer , dimension(3):: nk
  real(kind=8), dimension(3):: s0

  integer :: i,j,k, ig, ik, ibnd
  integer, allocatable :: kisort(:,:)
  real(kind=8), allocatable :: ei_k(:,:)

  open (unit=4, file='KG_GRIDS.DAT', &
       form='formatted', status='unknown') 
  write(4,'(3i4,3f12.6,i8)') nk, s0, ngm
  write(4,*) (ig1(ig), ig2(ig), ig3(ig), ig=1,ngm)
  close(unit=4)

  open (unit=40, file='vec.in', &
       form='unformatted', status='unknown')

  do ik = 1, nkstot
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)  
     write (40) ((DREAL(evc(ig,ibnd)), ig=1,npwx),ibnd=1,nbnd)
     write (40) ((1.0*DIMAG(evc(ig,ibnd)), ig=1,npwx),ibnd=1,nbnd)
  end do
  close(unit=40)

  allocate (kisort(npwx,nkstot))
  kisort = 0
  open (unit=41, file='val.in', &
       form='unformatted', status='unknown')
  do ik = 1, nkstot
     npw = npwx
     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, &
          kisort(1,ik), g2kin)
     ngk (ik) = npw
     write (6,fmt="(' k ',I4,3F12.6)") ik, xk(1,ik), xk(2,ik), xk(3,ik)
  end do
  write (41) npwx, nbnd, nkstot
  write (6,*) ' *** DEBUG 1 = ', npwx, nbnd, nkstot
  write (41)((kisort(ig,ik),ig=1,npwx),ik=1,nkstot)

  do i=1,npwx
     write(*,'(1(3x,f12.6),4i7)') dreal(evc(i,1)*conjg(evc(i,1))),kisort(i,1), &
          ig1(kisort(i,1)), ig2(kisort(i,1)), ig3(kisort(i,1))
  end do
 
  deallocate (kisort)
  allocate (ei_k(npwx,nkstot))
  ei_k = 0.0
  ei_k(1:nbnd,:) = et(1:nbnd,:) / 2.0d0  !  Rydberg to Hartree conversion
  write (41) ((ei_k(ig,ik),ig=1,npwx),ik=1,nkstot)
  deallocate (ei_k)
  write (41) (ngk(ik), ik=1,nkstot)
  write (41) (nbnd, ik=1,nkstot)
  write (41) nr1, nr2, nr3, ngm
  close (unit=41)

end subroutine write_wannier
