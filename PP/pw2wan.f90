!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
program wannier
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout, ionode
  USE mp_global,  ONLY : mpime
  USE mp,         ONLY : mp_bcast
  use pwcom  
  use para,       ONLY : kunit
  use io_files
  !
  implicit none
  integer , dimension(3):: nk
  real(kind=8), dimension(3):: s0
  integer :: ik, i, kunittmp
  CHARACTER(LEN=4) :: spin_component
  integer :: ispinw

  namelist / inputpp / tmp_dir, nk, s0, prefix, spin_component
  !
  call start_postproc (nd_nmbr)
  !
  !   set default values for variables in namelist
  !
  tmp_dir = './'
  prefix = ' '
  spin_component = 'none'
  nk = 0
  s0 = 0.d0
  !
  !     reading the namelist inputpp
  !
  read (5, inputpp)
  !
  !     Check of namelist variables
  !
  SELECT CASE ( TRIM( spin_component ) )
    CASE ( 'up' )
      ispinw = 1
    CASE ( 'down' )
      ispinw = 2
    CASE DEFAULT
      ispinw = 0
  END SELECT
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  call read_file  
  call openfil_pp
  !
  call mp_bcast( isk, 0 )
  !
  if (nk(1)==0 .and. nk(2)==0 .and. nk(3) == 0) then
     nk(1) = nk1
     nk(2) = nk2
     nk(3) = nk3
     s0(1) = k1/2d0
     s0(2) = k2/2d0
     s0(3) = k3/2d0
  end if

  IF( nspin /= 1 .AND. nspin /= 2 ) THEN
    CALL errore(' pw2wan ', ' nspin not allowed ', 1 )
  END IF

#if defined __PARA
  kunittmp = kunit
#else
  kunittmp = 1
#endif

  call write_wannier (nk, s0, kunittmp, ispinw)

  IF( ionode ) THEN
    WRITE( stdout, * ) 'K-points (reciprocal lattice coordinates) '
    do ik = 1, nkstot
      WRITE( stdout,fmt="(' ik = ',I3,3F10.6,I2)" ) &
        ik, xk(1,ik), xk(2,ik), xk(3,ik), isk(ik)
    end do
    !
    call cryst_to_cart ( nkstot, xk, at, -1)
    !
    WRITE( stdout, * ) 'K-points in unit of 2PI/alat (cartesian coordinates) '
    do ik = 1, nkstot
      WRITE( stdout,fmt="(' ik = ',I3,3F10.6,I2)" ) &
        ik, xk(1,ik), xk(2,ik), xk(3,ik), isk(ik)
    end do
  END IF
  !
  call stop_pp
  stop
end program wannier
!
!-----------------------------------------------------------------------
subroutine write_wannier (nk, s0, kunit, ispinw)
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  USE io_global,      ONLY : stdout
  use pwcom  
  USE wavefunctions_module,  ONLY : evc
  use io_files,       only : nd_nmbr, tmp_dir, prefix, iunwfc, nwordwfc
  use io_base,        only : write_restart_wfc
  use io_global,      only : ionode
  use mp_global,      only : nproc, nproc_pool, mpime
  use mp_global,      only : my_pool_id, my_image_id, intra_pool_comm
  use mp,             only : mp_sum, mp_max


  implicit none
  integer , dimension(3):: nk
  real(kind=8), dimension(3):: s0
  integer :: kunit
  integer :: ispinw

  integer :: i, j, k, ig, ik, ibnd, na, ngg, ikw
  integer, allocatable :: kisort(:)
  real(kind=8), allocatable :: ei_k(:,:)
  real(kind=8), allocatable :: ei_kw(:,:)
  real(kind=8), allocatable :: rat(:,:,:)
  real(kind=8) :: hmat(3,3), rr(3)
  integer, allocatable :: natom(:)
  integer :: npool, nkbl, nkl, nkr, npwx_g
  integer :: ike, iks, npw_g, ispin
  integer, allocatable :: ngk_g( : )
  integer, allocatable :: ngk_gw( : )
  integer, allocatable :: itmp( :, : )
  integer, allocatable :: igwk( : )

  real(kind=8) :: wfc_scal 
  logical :: twf0, twfm, twrite_wfc

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

  ! find out the global number of G vectors: ngm_g  
  ngm_g = ngm
  call mp_sum( ngm_g, intra_pool_comm )

  allocate ( ei_k ( nbnd, nkstot ) )  ! eigenvectors
  allocate ( ei_kw( nbnd, nkstot/nspin ) )  ! eigenvectors
  allocate ( rat( 3, nat, ntyp ) )   ! atomic positions
  allocate ( natom( ntyp ) )         ! number of atoms
  ei_k   = 0.0d0
  ei_kw  = 0.0d0
  rat    = 0.0d0
  natom  = 0 

  hmat = alat * at

  do i = 1, nat
    if( ityp(i) > ntyp ) &
      call errore( ' write_wannier ', ' type index out of range ', 1 )
    natom( ityp( i ) ) = natom( ityp( i ) ) + 1

    rat( :, natom( ityp( i ) ), ityp( i ) ) = tau( :, i ) 

    call cryst_to_cart ( 1, rat( :, natom( ityp( i ) ), ityp( i ) ), bg, -1)
  end do

  !  Open file launch.dat

  if( ionode ) then
    open (unit=40, file='launch.dat', form='unformatted', status='unknown')
  end if

  if( ionode ) then

    ! First write to launch.dat crystal base vectors
    WRITE( stdout,*) 'Crystal basis vectors, in unit of alat = ', alat
    do i = 1, 3
      WRITE( stdout,fmt="(' a(',I1,')',3F10.6)" ) i,at(1,i), at(2,i), at(3,i)
    end do
    write( 40 ) alat
    write( 40 ) ( at( i, 1 ), i = 1, 3 )  ! save A1
    write( 40 ) ( at( i, 2 ), i = 1, 3 )  ! save A2
    write( 40 ) ( at( i, 3 ), i = 1, 3 )  ! save A3

    ! write to launch.dat atomic positions 
    WRITE( stdout,*) 'Atomic positions ( lattice coordinates )'
    write( 40 ) ntyp
    do i = 1, ntyp
      write( 40 ) natom(i), atm(i)
      write( 40 ) ( ( rat( k, j, i ), k = 1, 3 ), j = 1, natom(i) )
      WRITE( stdout,*) 'Specie ', atm(i), ' atoms = ', natom(i)
      do j = 1, natom( i )
        WRITE( stdout,fmt="(' tau(',I1,')',3F10.6)" ) j, ( rat( k, j, i ), k = 1, 3 )
      end do
    end do
    write(6,*) 'Atomic positions ( From PW )'
    do i = 1, nat
      write(6,fmt="(' tau(',I1,')',3F10.6)" ) i, ( tau( k, i ), k = 1, 3 )
    end do
  
  end if

  ! collect all G vectors across processors within the pools
  allocate( itmp( 3, ngm_g ) )
  itmp = 0
  do  ig = 1, ngm
    itmp( 1, ig_l2g( ig ) ) = ig1( ig )
    itmp( 2, ig_l2g( ig ) ) = ig2( ig )
    itmp( 3, ig_l2g( ig ) ) = ig3( ig )
  end do
  call mp_sum( itmp, intra_pool_comm )

  ! write G space parameters and vectors
  if( ionode ) then
    write (40) ecutwfc, nbnd
    write (40) ( nk(i), i = 1, 3 ), ( s0(j), j = 1, 3 ), ngm_g
    write (40) ( itmp( 1, ig ), itmp( 2, ig ), itmp( 3, ig ), ig = 1, ngm_g )
  end if

  deallocate( itmp )

  ! build the G+k array indexes
  allocate ( kisort( npwx ) )
  do ik = 1, nks
     kisort = 0
     call gk_sort (xk (1, ik+iks-1), ngm, g, ecutwfc / tpiba2, npw, kisort(1), g2kin)
     call gk_l2gmap (ngm, ig_l2g(1), npw, kisort(1), igk_l2g(1,ik) )
     ngk (ik) = npw
  end do
  deallocate (kisort)

  ! compute the global number of G+k vectors for each k point
  allocate( ngk_g( nkstot ) )
  allocate( ngk_gw( nkstot/nspin ) )
  ngk_g = 0
  ngk_g( iks:ike ) = ngk( 1:nks )
  CALL mp_sum( ngk_g )

  do ik = 1, nkstot
    if( ionode ) then
      WRITE( stdout,fmt="(' k ',2I4,4F12.6,I2)")  &
        ik, ngk_g( ik), ecutwfc / tpiba2,  xk(1,ik), xk(2,ik), xk(3,ik), isk(ik)
    end if
  end do

  ! compute the Maximum G vector index among all G+k an processors
  npw_g = MAXVAL( igk_l2g(:,:) )
  CALL mp_max( npw_g )

  ! compute the Maximum number of G vector among all k points
  npwx_g = MAXVAL( ngk_g( 1:nkstot ) )

  if( ionode ) then
    write (40) npwx_g, nbnd, nkstot/nspin
    WRITE( stdout,*) 'Wave function dimensions ( npwx, nbnd, nkstot, nspin) = ', &
      npwx_g, nbnd, nkstot, nspin
  end if

  ! for each k point build and write the global G+k indexes array
  WRITE( stdout,*) 'combining indexes'

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
      WRITE( stdout,*) ' ik, ngg, ngk_g = ', ik, ngg, ngk_g( ik )
    end if
    deallocate( itmp )
    if( ionode ) then
      IF( ( ispinw == 0 ) .OR. ( isk(ik) == ispinw ) ) THEN
        write (40)( igwk(ig), ig = 1, npwx_g )
      END IF 
    end if
  end do

  deallocate( igwk )


#ifdef __PARA
  call poolrecover (et, nbnd, nkstot, nks)
#endif
 
  ei_k(1:nbnd,:) = et(1:nbnd,:) / 2.0d0  !  Rydberg to Hartree conversion

  !WRITE( stdout, * ) isk

  if( ionode ) then
    ikw = 0
    DO ik = 1, nkstot
      IF( ( ispinw == 0 ) .OR. ( isk(ik) == ispinw ) ) THEN
        ikw = ikw + 1
        ei_kw( :, ikw ) = ei_k( :, ik )
        ngk_gw( ikw ) = ngk_g( ik )
      END IF
    END DO
    write (40) ( ( ei_kw( i, ik ), i = 1, nbnd ), ik = 1, ikw )
    write (40) ( ngk_gw( ik ), ik = 1, ikw )
    write (40) ( nbnd, ik = 1, ikw )
    write (40) nr1, nr2, nr3, ngm_g, npw_g
    WRITE( stdout,*) 'Grid  ( nr1, nr2, nr3, ngm_g, npw_g, ikw ) = ', nr1, nr2, nr3, ngm_g, npw_g, ikw
  end if

  wfc_scal = 1.0d0
  twf0 = .true.
  twfm = .false.

  ! WRITE( stdout, * ) ispinw, iks, ike, nkstot

  do ik = 1, nkstot
     IF( (ik >= iks) .AND. (ik <= ike) ) THEN
       call davcio (evc, nwordwfc, iunwfc, (ik-iks+1), - 1)
     END IF
     ispin = isk( ik )
     IF( ( ispinw == 0 ) .OR. ( isk(ik) == ispinw ) ) THEN
       ! WRITE( stdout,*) ' ### ', ik,nkstot,iks,ike,kunit,nproc,nproc_pool ! DEBUG
       CALL write_restart_wfc(40, ik, nkstot, kunit, ispin, nspin, &
         wfc_scal, evc, twf0, evc, twfm, npw_g, nbnd, igk_l2g(:,ik-iks+1), ngk(ik-iks+1) )
     END IF
  end do

  deallocate (ei_k)
  deallocate (ei_kw)
  deallocate ( rat )
  deallocate ( natom )
  deallocate ( ngk_g )
  deallocate ( ngk_gw )

  if( ionode ) then
    close ( unit = 40 )
  end if

end subroutine write_wannier
