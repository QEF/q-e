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
  use para, only: kunit
  use io_files

  implicit none
  character(len=3) :: nodenumber
  integer , dimension(3):: nk
  real(kind=8), dimension(3):: s0
  integer :: ik, i, kunittmp

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

#if defined __PARA
  kunittmp = kunit
#else
  kunittmp = 1
#endif

  call write_wannier (nk, s0, kunittmp)

  write( 6, * ) 'K-points (reciprocal lattice coordinates) '
  do ik = 1, nkstot
    write(6,fmt="(' ik = ',I3,3F10.6)" ) ik, xk(1,ik), xk(2,ik), xk(3,ik)
  end do
  !
  call cryst_to_cart ( nkstot, xk, at, -1)
  !
  write( 6, * ) 'K-points in unit of 2PI/alat (cartesian coordinates) '
  do ik = 1, nkstot
    write(6,fmt="(' ik = ',I3,3F10.6)" ) ik, xk(1,ik), xk(2,ik), xk(3,ik)
  end do
  !
  return  
end subroutine do_wannier
!
!-----------------------------------------------------------------------
subroutine write_wannier (nk, s0, kunit)
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  use pwcom  
  use io_files, only: nd_nmbr, tmp_dir, prefix
  use io_base, only: write_restart_wfc
  use io_global, only: ionode
  use mp_global, only: nproc, nproc_pool
  use mp_global, only: my_pool_id, intra_pool_comm, inter_pool_comm
  use mp, only: mp_sum, mp_max


  implicit none
  integer , dimension(3):: nk
  real(kind=8), dimension(3):: s0
  integer :: kunit

  integer :: i, j, k, ig, ik, ibnd, na, ngg
  integer, allocatable :: kisort(:)
  real(kind=8), allocatable :: ei_k(:,:)
  real(kind=8), allocatable :: rat(:,:,:)
  integer, allocatable :: natom(:)
  integer :: npool, nkbl, nkl, nkr, npwx_g
  integer :: ike, iks, npw_g, ispin
  integer, allocatable :: ngk_g( : )
  integer, allocatable :: itmp( :, : )
  integer, allocatable :: igwk( : )

  real(kind=8) :: wfc_scal 
  logical :: twf0, twfm

  IF( nkstot > 0 ) THEN

     IF( ( kunit < 1 ) .OR. ( MOD( nkstot, kunit ) /= 0 ) ) &
       CALL errore( ' write_wannier ',' wrong kunit ', 1 )

     IF( ( nproc_pool > nproc ) .OR. ( MOD( nproc, nproc_pool ) /= 0 ) ) &
       CALL errore( ' write_wannier ',' nproc_pool ', 1 )

     !  find out the number of pools
     npool = nproc / nproc_pool

     !  find out number of k points blocks
     nkbl = nkstot / kunit

     !  k points per pool
     nkl = kunit * ( nkbl / npool )

     !  find out the reminder
     nkr = ( nkstot - nkl * npool ) / kunit

     !  Assign the reminder to the first nkr pools
     IF( my_pool_id < nkr ) nkl = nkl + kunit

     !  find out the index of the first k point in this pool
     iks = nkl * my_pool_id + 1
     IF( my_pool_id >= nkr ) iks = iks + nkr * kunit

     !  find out the index of the last k point in this pool
     ike = iks + nkl - 1

  END IF

  ngm_g = ngm
  call mp_sum( ngm_g , intra_pool_comm )

  allocate ( ei_k( nbnd, nkstot ) )
  allocate ( rat( 3, nat, ntyp ) )
  allocate ( natom( ntyp ) )
  ei_k   = 0.0d0
  natom  = 0 

  do i = 1, nat
    if( ityp(i) > ntyp ) &
      call errore( ' write_wannier ', ' type index out of range ', 1 )
    natom( ityp( i ) ) = natom( ityp( i ) ) + 1
    rat( :, natom( ityp( i ) ), ityp( i ) ) = tau( :, i )
    call cryst_to_cart ( 1, rat( :, natom( ityp( i ) ), ityp( i ) ), at, -1)
  end do

  !  Open file launch.dat

  if( ionode ) then
    open (unit=40, file='launch.dat', form='unformatted', status='unknown')
  end if

  if( ionode ) then

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
  
  end if

  allocate( itmp( 3, ngm_g ) )
  itmp = 0
  do  ig = 1, ngm
    itmp( 1, ig_l2g( ig ) ) = ig1( ig )
    itmp( 2, ig_l2g( ig ) ) = ig2( ig )
    itmp( 3, ig_l2g( ig ) ) = ig3( ig )
  end do
  call mp_sum( itmp , intra_pool_comm )

  if( ionode ) then
    write (40) ecutwfc, nbnd
    write (40) ( nk(i), i = 1, 3 ), ( s0(j), j = 1, 3 ), ngm_g
    write (40) ( itmp( 1, ig ), itmp( 2, ig ), itmp( 3, ig ), ig = 1, ngm_g )
  end if

  deallocate( itmp )

  allocate ( kisort( npwx ) )
  do ik = 1, nks
     kisort = 0
     npw = npwx
     call gk_sort (xk (1, ik+iks-1), ngm, g, ecutwfc / tpiba2, npw, kisort(1), g2kin)
     call gk_l2gmap (ngm, ig_l2g(1), npw, kisort(1), igk_l2g(1,ik) )
     ngk (ik) = npw
  end do
  deallocate (kisort)

  allocate( ngk_g( nkstot ) )
  ngk_g = 0
  ngk_g( iks:ike ) = ngk( 1:nks )
  CALL mp_sum( ngk_g )

  do ik = 1, nkstot
    if( ionode ) then
      write (6,fmt="(' k ',2I4,4F12.6)")  &
        ik, ngk_g( ik), ecutwfc / tpiba2,  xk(1,ik), xk(2,ik), xk(3,ik)
    end if
  end do

  npw_g = MAXVAL( igk_l2g(:,:) )
  CALL mp_max( npw_g )

  npwx_g = MAXVAL( ngk_g( 1:nkstot ) )

  if( ionode ) then
    write (40) npwx_g, nbnd, nkstot
    write (6,*) 'Wave function dimensions ( npwx, nbnd, nkstot) = ', npwx_g, nbnd, nkstot
  end if

  write (6,*) 'combining indexes'

  allocate( igwk( npwx_g ) )

  do ik = 1, nkstot
    igwk = 0
    allocate( itmp( npw_g, 1 ) )
    itmp = 0
    if( ik >= iks .AND. ik <= ike ) then 
      do  ig = 1, ngk( ik-iks+1 )
        itmp( igk_l2g( ig, ik-iks+1 ), 1 ) = igk_l2g( ig, ik-iks+1 ) 
      end do
    end if
    call mp_sum( itmp )
    ngg = 0
    do  ig = 1, npw_g
      if( itmp( ig, 1 ) == ig ) then
        ngg = ngg + 1
        igwk( ngg ) = ig
      end if
    end do
    if( ngg /= ngk_g( ik ) ) then
      write(6,*) ' ik, ngg, ngk_g = ', ik, ngg, ngk_g( ik )
    end if
    deallocate( itmp )
    if( ionode ) then
      write (40)( igwk(ig), ig = 1, npwx_g )
    end if
  end do

  deallocate( igwk )


#ifdef __PARA
  call poolrecover (et, nbnd, nkstot, nks)
#endif
 
  ei_k(1:nbnd,:) = et(1:nbnd,:) / 2.0d0  !  Rydberg to Hartree conversion

  if( ionode ) then
    write (40) ( ( ei_k( i, ik ), i = 1, nbnd ), ik = 1, nkstot )
    write (40) ( ngk_g( ik ), ik = 1, nkstot )
    write (40) ( nbnd, ik = 1, nkstot )
    write (40) nr1, nr2, nr3, ngm_g, npw_g
    write (6,*) 'Grid  ( nr1, nr2, nr3, ngm_g, npw_g ) = ', nr1, nr2, nr3, ngm_g, npw_g
  end if

  

  wfc_scal = 1.0d0
  twf0 = .true.
  twfm = .false.

  do ik = 1, nkstot
     IF( (ik >= iks) .AND. (ik <= ike) ) THEN
       call davcio (evc, nwordwfc, iunwfc, (ik-iks+1), - 1)
     END IF
     ispin = isk( ik )
     !  ! write(6,*) ' ### ', ik,nkstot,iks,ike,kunit,nproc,nproc_pool ! DEBUG
     CALL write_restart_wfc(40, ik, nkstot, kunit, ispin, nspin, &
         wfc_scal, evc, twf0, evc, twfm, npw_g, nbnd, igk_l2g(:,ik-iks+1), ngk(ik-iks+1) )
  end do

  deallocate (ei_k)
  deallocate ( rat )
  deallocate ( natom )
  deallocate ( ngk_g )

  if( ionode ) then
    close ( unit = 40 )
  end if

end subroutine write_wannier
