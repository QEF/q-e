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

  namelist / inputpp / tmp_dir, nk, s0, prefix
  !
  nd_nmbr = nodenumber  
  !
  !   set default values for variables in namelist
  !
  tmp_dir = './'
  prefix = ' '
  nk = 0
  s0 = 0.d0
  !
  !     reading the namelist inputpp
  !
  read (5, inputpp)
  !
  !     Check of namelist variables
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  call read_file  
  call openfil  

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

  write( 6, * ) 'K-points (reciprocal lattice coordinates) '
  do ik = 1, nks
    write(6,fmt="(' ik = ',I3,3F10.6)" ) ik, xk(1,ik), xk(2,ik), xk(3,ik)
  end do
  !
  call cryst_to_cart ( nks, xk, at, -1)
  !
  write( 6, * ) 'K-points in unit of 2PI/alat (cartesian coordinates) '
  do ik = 1, nks
    write(6,fmt="(' ik = ',I3,3F10.6)" ) ik, xk(1,ik), xk(2,ik), xk(3,ik)
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

  integer :: i,j,k, ig, ik, ibnd, na
  integer, allocatable :: kisort(:,:)
  real(kind=8), allocatable :: ei_k(:,:)
  real(kind=8), allocatable :: rat(:,:,:)
  integer, allocatable :: natom(:)

  allocate ( kisort( npwx, nkstot ) )
  allocate (ei_k(npwx,nkstot))
  allocate ( rat( 3, nat, ntyp ) )
  allocate ( natom( ntyp ) )
  kisort = 0
  ei_k = 0.0
  natom = 0 

  do i = 1, nat
    if( ityp(i) > ntyp ) &
      call errore( ' write_wannier ', ' type index out of range ', 1 )
    natom( ityp( i ) ) = natom( ityp( i ) ) + 1
    rat( :, natom( ityp( i ) ), ityp( i ) ) = tau( :, i )
    call cryst_to_cart ( 1, rat( :, natom( ityp( i ) ), ityp( i ) ), at, -1)
  end do

  !  Open file launch.dat

  open (unit=40, file='launch.dat', form='unformatted', status='unknown')

  write(6,*) 'Crystal basis vectors, in unit of alat = ', alat
  do i = 1, 3
    write(6,fmt="(' a(',I1,')',3F10.6)" ) i,at(1,i), at(2,i), at(3,i)
  end do
  write( 40 ) alat
  write( 40 ) ( at( i, 1 ), i = 1, 3 )  ! save A1
  write( 40 ) ( at( i, 2 ), i = 1, 3 )  ! save A2
  write( 40 ) ( at( i, 3 ), i = 1, 3 )  ! save A3

  write(6,*) 'Atomic positions ( lattice coordinates )'
  write( 40 ) ntyp
  do i = 1, ntyp
    write( 40 ) natom(i), atm(i)
    write( 40 ) ( ( rat( k, j, i ), k = 1, 3 ), j = 1, natom(i) )
    write(6,*) 'Specie ', atm(i), ' atoms = ', natom(i)
    do j = 1, natom( i )
      write(6,fmt="(' tau(',I1,')',3F10.6)" ) j, ( rat( k, j, i ), k = 1, 3 )
    end do
  end do

  write (40) ecutwfc, nbnd
  write (40) ( nk(i), i = 1, 3 ), ( s0(j), j = 1, 3 ), ngm
  write (40)( ig1(ig), ig2(ig), ig3(ig), ig = 1, ngm )

  do ik = 1, nkstot
     npw = npwx
     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, kisort(1,ik), g2kin)
     ngk (ik) = npw
     write (6,fmt="(' k ',2I4,4F12.6)") ik, npw, ecutwfc / tpiba2,  xk(1,ik), xk(2,ik), xk(3,ik)
  end do
  write (40) npwx, nbnd, nkstot
  write (6,*) 'Wave function dimensions ( npwx, nbnd, nkstot) = ', npwx, nbnd, nkstot
  write (40)( ( kisort(ig,ik), ig = 1, npwx ), ik = 1, nkstot )

  do i = 1, npwx
     ! write(*,'(1(3x,f12.6),4i7)') dreal(evc(i,1)*conjg(evc(i,1))),kisort(i,1), &
     !     ig1(kisort(i,1)), ig2(kisort(i,1)), ig3(kisort(i,1))
  end do
 
  ei_k(1:nbnd,:) = et(1:nbnd,:) / 2.0d0  !  Rydberg to Hartree conversion

  write (40) ( ( ei_k( ig, ik ), ig = 1, npwx ), ik = 1, nkstot )
  write (40) ( ngk( ik ), ik = 1, nkstot )
  write (40) ( nbnd, ik = 1, nkstot )
  write (40) nr1, nr2, nr3, ngm

  do ik = 1, nkstot
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)  
     write (40) ((DREAL(evc(ig,ibnd)), ig=1,npwx),ibnd=1,nbnd)
     write (40) ((1.0*DIMAG(evc(ig,ibnd)), ig=1,npwx),ibnd=1,nbnd)
  end do

  deallocate (kisort)
  deallocate (ei_k)
  deallocate ( rat )
  deallocate ( natom )

  close ( unit = 40 )

end subroutine write_wannier
