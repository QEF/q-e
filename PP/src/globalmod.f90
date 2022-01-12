!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! Authors: Ivan Carnimeo, Pietro Delugas (September 2021)
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
MODULE globalmod
  USE kinds,                ONLY : dp
implicit none
  !INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  ! 
  ! method
  character(len=50) :: method
  !
  ! band indexes
  integer :: Nb
  !
  ! uniform grid of q-points in the IBZ (cart. coord. in units 2pi/alat)
  integer :: Nq, Nlines
  real(dp), allocatable :: q(:,:), eq(:,:)
  ! 
  ! path of k-points in the IBZ
  integer :: Nk
  real(dp), allocatable :: k(:,:), ek(:,:), t(:)
  !
  ! list of k-points for building the path in the BZ (cart. coord. in units 2pi/alat)
  integer :: Nkl
  integer,  allocatable :: kln(:)
  real(dp), allocatable :: kl(:,:)
  ! 
  ! Crystal data 
  real(dp) :: at(3,3)   ! real-space lattice translation vectors (cart. coord. in units of alat) a_i(:) = at(:,i)/alat
  real(dp) :: bg(3,3)   ! reciprocal-lattice (cart. coord. in units 2 pi/alat)                   b_i(:) = bg(:,i)/tpiba
  integer :: Nsym                     ! number of symmetry operations of the point group of the crystal  
  real(dp), allocatable :: Op(:,:,:)      ! 3x3 symmetry matrices: Op(Nsym,3,3) (cart. units)
  real(dp), allocatable :: Op_tmp(:,:,:)  ! this is just a buffer to convert Op in cartesian units. 
                                          ! quite redundant here, but useful to use s_axis_to_cart without modifications 
  !
CONTAINS
!----------------------------------------------------------------------------
subroutine read_input ()
!
! read the input file and make all allocations 
! 
use fouriermod,  ONLY : NMax, allocate_fourier, NC, C, NUser, VecUser
use shepardmod,  ONLY : PMetric, ScaleSphere 
use qes_read_module, only: qes_read 
use qes_types_module, only: band_structure_type, atomic_structure_type, symmetries_type, basis_set_type
use fox_dom 
implicit none
  integer :: ik, iq, ib, ikl, ivec
  integer :: ilines
  integer :: isym
  character(len=50) :: string
  ! 
  type (node),pointer      :: doc, pnode1, pnode2 
  type (band_structure_type) :: bandstr
  type (atomic_structure_type) :: atstr 
  type (symmetries_type)       :: symstr 
  type (basis_set_type)        :: basisstr 
  !
  ! read data from the xml_file 
  doc => parsefile('pwscf.xml') 
  pnode1 => item(getElementsByTagname(doc, 'output'),0) 
  pnode2 => item(getElementsByTagname(pnode1, 'atomic_structure'),0) 
  call qes_read (pnode2, atstr) 
  pnode2 => item(getElementsByTagname(pnode1, 'symmetries'),0) 
  call qes_read (pnode2, symstr) 
  pnode2 => item(getElementsByTagname(pnode1, 'band_structure'),0) 
  call qes_read(pnode2, bandstr) 
  pnode2 => item(getElementsByTagname(pnode1, 'basis_set'),0) 
  call qes_read(pnode2, basisstr)
  call destroy (doc) 
  nullify (doc, pnode1, pnode2) 

  read(*, *) 
  read(*, *) method
  write(*,*) 'Interpolation method: ', method
  if( TRIM(method).ne.'shepard'.and.TRIM(method).ne.'shepard-sphere'&
        .and.TRIM(method).ne.'fourier'.and.TRIM(method).ne.'fourier-diff' ) then
    write(*,*) 'ERROR: Wrong method ', method
    stop
  end if 
  !
  ! optionally add user-given star functions to the basis set
  ! 
  NUser = 0 
  read(*,*)
  read(*,*) NUser
  if(NUser.gt.0) then 
    write(*,*) NUser, ' user-given star functions found'
    allocate( VecUser(3,NUser) ) 
    VecUser = 0.0d0
    do ivec = 1, NUser
      read(*,*) VecUser(1:3,ivec) 
      write(*,'(3f12.6)') VecUser(1:3,ivec)
    end do 
  elseif(NUser.eq.0) then 
    write(*,*) 'No user-given star functions provided'
  else
    write(*,*) 'ERROR: Wrong NUser'
    write(*,*) '       Please provide non-negative NUser'
    stop
  end if 
  !
  ! read specific parameters for the interpolation methods
  ! 
  if( TRIM(method).eq.'shepard'.or.TRIM(method).eq.'shepard-sphere' ) THEN 
    !
    ! read parameters for Shepard interpolation 
    ! 
    read(*,*)
    read(*,*) PMetric, ScaleSphere
    write(*,*) 'PMetric: ', PMetric, 'ScaleSphere ', ScaleSphere
  elseif( TRIM(method).eq.'fourier'.or.TRIM(method).eq.'fourier-diff' ) 
    !
    ! read parameters for Fourier interpolation 
    ! 
    read(*, *) 
    read(*, *) string, NMax
    if( NMax.le.0) then 
      write(*,*) 'Wrong NMax: ', NMax
      write(*,*) 'NMax must be greater than 0 '
      stop
    end if
    read(*, *) string, NC
    if( NC.le.0) then 
      write(*,*) 'Wrong NC: ', NC  
      write(*,*) 'NC must be greater than 0 '
      stop
    end if
    Call allocate_fourier( )
    read(*, *) string, C(1:NC)
    write(*,*) NC, ' coefficients read for rho expansion: ', C(:)
  end if 
  !
  ! read the list of Nkl special points
  !
  read(*, *) 
  read(*, *) Nkl
  !
  write(*,*) Nkl, ' special points read'
  !
  ! create the abscissa values for bands plotting
  !  
  allocate( kl(3,Nkl), kln(Nkl ) )
  !
  do ikl = 1, Nkl
    read(*, *) kl(:,ikl), kln(ikl) 
    write(*,'(3f12.6,I5)') kl(:,ikl), kln(ikl) 
  end do 
  !
  Nlines = Nkl - 1
  Nk = sum(kln(1:Nlines)) + 1
  write(*,*) 'Creating ', Nlines, ' lines connecting ', Nkl, ' special points with ', Nk, ' k-points' 
  !
  allocate( k(3,Nk), t(Nk) )
  !
  Call build_kpath()
  !
  deallocate( kl, kln )
  !
  ! read the uniform grid of q-points from xml
  !
  Nq = bandstr%nks  
  Nb = bandstr%nbnd 
  !
  write(*,*) Nq, ' points on the uniform grid, ', Nb, ' bands'
  !write(*,'(A)') 'iq, q(iq, :), e(iq, :)'
  !
  allocate( q(3, Nq), eq(Nq, Nb), ek(Nk,Nb)  )
  !
  do iq = 1, Nq
    q(:,iq) = bandstr%ks_energies(iq)%k_point%k_point(:) 
  end do 
  do iq = 1, Nq
    eq(iq,:) = bandstr%ks_energies(iq)%eigenvalues%vector(:)*27.211386245988 
    !write(*,'(I5,11f12.6)') iq, q(iq, :), eq(iq, :)
  end do 
  !
  ! read from xml crystalline group specifications 
  ! (direct and reciprocal vectors, symmetry operations)
  !
  at(1:3,1) = atstr%cell%a1 / atstr%alat
  at(1:3,2) = atstr%cell%a2 / atstr%alat 
  at(1:3,3) = atstr%cell%a3 / atstr%alat
  write(*,*) ' Crystal lattice vectors found  '
  write(*,*) 'Ra: ' , at(:,1)
  write(*,*) 'Rb: ' , at(:,2)
  write(*,*) 'Rc: ' , at(:,3)
  bg(1:3,1) = basisstr%reciprocal_lattice%b1
  bg(1:3,2) = basisstr%reciprocal_lattice%b2
  bg(1:3,3) = basisstr%reciprocal_lattice%b3  
  write(*,*) ' Reciprocal lattice vectors found  '
  write(*,*) 'Ga: ' , bg(:,1)
  write(*,*) 'Gb: ' , bg(:,2)
  write(*,*) 'Gc: ' , bg(:,3)
  Nsym = symstr%nsym   
  write(*,*) Nsym, ' symmetry operations found '
  !
  allocate( Op(3,3,Nsym), Op_tmp(3,3,Nsym) )
  !
  do isym = 1, Nsym
    Op_tmp(:,:,isym) = reshape(symstr%symmetry(isym)%rotation%matrix,[3,3])  
  end do 
  Call s_axis_to_cart()
  deallocate(Op_tmp)
  !
  return
  !
end subroutine read_input 
!----------------------------------------------------------------------------
subroutine build_kpath ()
!
! build the path of k-points connecting the Nkl special points
! 
implicit none
  integer :: ik, iik, iline
  real(dp) :: tk, dt, modt
  !
  k = 0.0d0
  t = 0.0d0
  ! do the first point (t = 0) 
  ik = 1  
  t(ik) = 0.0d0
  k(1, ik) = kl(1, 1) 
  k(2, ik) = kl(2, 1) 
  k(3, ik) = kl(3, 1) 
  ! loop over lines
  do iline = 1, Nlines   
    tk = 0.0d0
    do iik = 1, kln(iline)
      ik = ik + 1  
      modt = sqrt( (kl(1,iline+1)-kl(1,iline))**2 + (kl(2,iline+1)-kl(2,iline))**2 + (kl(3,iline+1)-kl(3,iline))**2 ) 
      tk = tk + 1.0d0/dble(kln(iline))
      dt = modt/dble(kln(iline))
      t(ik) = t(ik-1) + dt
      k(1,ik) = kl(1,iline) + (kl(1,iline+1) - kl(1,iline)) * tk   
      k(2,ik) = kl(2,iline) + (kl(2,iline+1) - kl(2,iline)) * tk 
      k(3,ik) = kl(3,iline) + (kl(3,iline+1) - kl(3,iline)) * tk 
    end do 
  end do 
  !
  return
  !
end subroutine build_kpath 
!----------------------------------------------------------------------------
subroutine print_bands (label)
!
! print band structure
!
implicit none
  integer :: ik
  character(len=*), intent(in) :: label 
  character(len=100) :: formt, filename
  !
  write(formt,'(A,I5,A)') '(', Nb+1 ,'f24.6)'  
  write(filename, '(A,A)')  TRIM(label),'.dat'
  !
  write(*,*) 'writing band structure on ', filename
  !
  open(2, file=filename, status='unknown')
  !
  do ik = 1, Nk
    write(2,formt) t(ik), ek(ik,:)
    !write(*,formt) t(ik), ek(ik,:)
  end do 
  !
  close(2)
  !
  return
  !
end subroutine print_bands
!----------------------------------------------------------------------------
subroutine deallocate_global ()
implicit none
  deallocate(q, eq, k, ek, t, Op)
end subroutine deallocate_global
!----------------------------------------------------------------------
SUBROUTINE s_axis_to_cart()
  !----------------------------------------------------------------------
  !! This routine transforms symmetry matrices expressed in the
  !! basis of the crystal axis into rotations in cartesian axis.
!
!civn 2FIX: better remove this one and use PW/src/symm_base.f90 instead 
!           (change Op_tmp --> sr and   Op --> s) 
  !
  IMPLICIT NONE
  !
  INTEGER :: isym
  REAL(DP) :: sa(3,3), sb(3,3)
  !
  DO isym = 1,nsym
     sa(:,:) = DBLE( Op_tmp(:,:,isym) )
     sb = MATMUL( bg, sa )
     Op(:,:,isym) = MATMUL( at, TRANSPOSE(sb) )
  ENDDO
  !
 END SUBROUTINE s_axis_to_cart
!----------------------------------------------------------------------------
END MODULE
!----------------------------------------------------------------------------
