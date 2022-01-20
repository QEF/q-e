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
  ! 
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
  logical :: trough = .false.    
  logical :: tuser  = .false.    
  !
CONTAINS
!----------------------------------------------------------------------------
subroutine card_user_stars( input_line )
use fouriermod ,  only : NUser, VecUser
use parser,       ONLY : read_line
implicit none
  integer  :: ivec
  LOGICAL            :: tend,terr
  CHARACTER(len=256) :: input_line
  !
  IF ( tuser  ) THEN
     CALL errore( ' card_user_stars  ', ' two occurrences', 2 )
  ENDIF
  !
  ! ... input Star vectors 
  !
  CALL read_line( input_line, end_of_file = tend, error = terr )
  IF (tend) GOTO 10
  IF (terr) GOTO 20
  READ(input_line, *, END=10, ERR=20) NUser 
  !
  if(NUser.gt.0) then 
    allocate( VecUser(3,NUser) ) 
    !
    do ivec = 1, NUser
      CALL read_line( input_line, end_of_file = tend, error = terr )
      IF (tend) GOTO 10
      IF (terr) GOTO 20
      READ(input_line,*, END=10, ERR=20) VecUser(1:3,ivec)
    end do 
  end if 
  !
  tuser  = .true. 
  !
  RETURN
  !
10 CALL errore ('card_user_stars', ' end of file while reading roughness function ', 1) 
20 CALL errore ('card_user_stars', ' error while reading roughness function', 1) 
  !
end subroutine card_user_stars
!----------------------------------------------------------------------------
subroutine card_roughness( input_line )
use fouriermod ,  only : RoughN, RoughC
use parser,       ONLY : read_line
implicit none
  LOGICAL            :: tend,terr
  CHARACTER(len=256) :: input_line
  !
  IF ( trough ) THEN
     CALL errore( ' card_roughness  ', ' two occurrences', 2 )
  ENDIF
  !
  ! ... input coefficients for the roughness function
  !
  CALL read_line( input_line, end_of_file = tend, error = terr )
  IF (tend) GOTO 10
  IF (terr) GOTO 20
  READ(input_line, *, END=10, ERR=20) RoughN 
  !
  if(RoughN.gt.1) then 
    deallocate( RoughC ) 
    allocate( RoughC(RoughN) ) 
  end if 
  !
  if(RoughN.gt.0) then 
    CALL read_line( input_line, end_of_file = tend, error = terr )
    IF (tend) GOTO 10
    IF (terr) GOTO 20
    READ(input_line,*, END=10, ERR=20) RoughC(1:RoughN) 
    !
  else
    Call errore( ' card_roughness  ', ' RoughN must be a positive integer ', 2 )
  end if
  !
  trough = .true. 
  !
  RETURN
  !
10 CALL errore ('card_roughness', ' end of file while reading roughness function ', 1) 
20 CALL errore ('card_roughness', ' error while reading roughness function', 1) 
  !
end subroutine card_roughness
!----------------------------------------------------------------------------
subroutine read_xml_input ()
!
! read the xml input file and make all allocations 
! 
USE io_global, ONLY : stdout
use qes_read_module,  ONLY : qes_read 
use qes_types_module, ONLY : band_structure_type, atomic_structure_type, symmetries_type, basis_set_type
use fox_dom 
implicit none
  integer :: iq
  integer :: isym
  ! 
  type (node),pointer          :: doc, pnode1, pnode2 
  type (band_structure_type)   :: bandstr
  type (atomic_structure_type) :: atstr 
  type (symmetries_type)       :: symstr 
  type (basis_set_type)        :: basisstr 
  !
  ! read data from the xml_file 
  !
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
  !
  ! read the uniform grid of q-points from xml
  !
  Nq = bandstr%nks  
  Nb = bandstr%nbnd 
  !
  write(stdout,'(2(I5,A))') Nq, ' points on the uniform grid, ', Nb, ' bands'
  !write(stdout,'(A)') 'iq, q(iq, :), e(iq, :)'
  !
  allocate( q(3, Nq), eq(Nq, Nb), ek(Nk,Nb)  )
  !
  do iq = 1, Nq
    q(:,iq) = bandstr%ks_energies(iq)%k_point%k_point(:) 
  end do 
  do iq = 1, Nq
    eq(iq,:) = bandstr%ks_energies(iq)%eigenvalues%vector(:)*27.211386245988 
    !write(stdout, '(I5,11f12.6)') iq, q(iq, :), eq(iq, :)
  end do 
  !
  ! read from xml crystalline group specifications 
  ! (direct and reciprocal vectors, symmetry operations)
  !
  at(1:3,1) = atstr%cell%a1 / atstr%alat
  at(1:3,2) = atstr%cell%a2 / atstr%alat 
  at(1:3,3) = atstr%cell%a3 / atstr%alat
  write(stdout,'(A)')        ' Crystal lattice vectors found  '
  write(stdout,'(A,3f12.6)') 'Ra: ' , at(:,1)
  write(stdout,'(A,3f12.6)') 'Rb: ' , at(:,2)
  write(stdout,'(A,3f12.6)') 'Rc: ' , at(:,3)
  bg(1:3,1) = basisstr%reciprocal_lattice%b1
  bg(1:3,2) = basisstr%reciprocal_lattice%b2
  bg(1:3,3) = basisstr%reciprocal_lattice%b3  
  write(stdout,'(A)')        ' Reciprocal lattice vectors found  '
  write(stdout,'(A,3f12.6)') 'Ga: ' , bg(:,1)
  write(stdout,'(A,3f12.6)') 'Gb: ' , bg(:,2)
  write(stdout,'(A,3f12.6)') 'Gc: ' , bg(:,3)
  Nsym = symstr%nsym   
  write(stdout,'(I5,A)') Nsym, ' symmetry operations found '
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
end subroutine read_xml_input 
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
