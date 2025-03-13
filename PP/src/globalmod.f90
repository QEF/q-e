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
  USE kinds,            ONLY : dp
  USE io_global,        ONLY : stdout
implicit none
  ! 
  ! a string describing the method used for interpolation  
  CHARACTER(len=80) :: method = ' '
  !
  ! band indexes
  integer :: Nb
  !
  ! uniform grid of q-points in the IBZ (cart. coord. in units 2pi/alat)
  integer :: Nq
  real(dp), allocatable :: q(:,:), eq(:,:)
  ! 
  ! band energies of the path of k-points in the IBZ
  real(dp), allocatable :: ek(:,:)
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
subroutine read_xml_input ()
!
! read the xml input file and make all allocations 
! 
use qes_read_module,  ONLY : qes_read 
use qes_types_module, ONLY : band_structure_type, atomic_structure_type, symmetries_type, basis_set_type
#if defined (__fox)
  USE FoX_dom,   ONLY : parseFile, node, item, getElementsByTagname, &
                        destroy
#else
  USE     dom,   ONLY : parseFile, node, item, getElementsByTagname, &
                        destroy
#endif
use input_parameters, ONLY : nkstot 
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
  write(stdout,'(A,2(I5,A))') 'The uniform grid contains ', Nq, ' q-points and ', Nb, ' bands'
  !write(stdout,'(A)') 'iq, q(iq, :), e(iq, :)'
  !
  allocate( q(3, Nq), eq(Nq, Nb), ek(nkstot,Nb)  )
  !
  do iq = 1, Nq
    q(:,iq) = bandstr%ks_energies(iq)%k_point%k_point(:) 
  end do 
  do iq = 1, Nq
    eq(iq,:) = bandstr%ks_energies(iq)%eigenvalues%vector(:)*27.2113862459880d0 
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
  write(stdout,'(A,3f12.4)') '      Ra: ' , at(:,1)
  write(stdout,'(A,3f12.4)') '      Rb: ' , at(:,2)
  write(stdout,'(A,3f12.4)') '      Rc: ' , at(:,3)
  bg(1:3,1) = basisstr%reciprocal_lattice%b1
  bg(1:3,2) = basisstr%reciprocal_lattice%b2
  bg(1:3,3) = basisstr%reciprocal_lattice%b3  
  write(stdout,'(A)')        ' Reciprocal lattice vectors found  '
  write(stdout,'(A,3f12.4)') '      Ga: ' , bg(:,1)
  write(stdout,'(A,3f12.4)') '      Gb: ' , bg(:,2)
  write(stdout,'(A,3f12.4)') '      Gc: ' , bg(:,3)
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
subroutine deallocate_global ()
USE input_parameters, ONLY : xk
implicit none
  deallocate(q, eq, ek, Op, xk)
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
subroutine print_bands (label)
!
! print band structure
!
USE kinds, ONLY : dp
USE input_parameters, ONLY : nkstot, xk
implicit none
  character(len=*), intent(in) :: label 
  !
  real(dp) :: x  
  integer :: ik
  character(len=100) :: formt, filename
  !
  write(formt,'(A,I5,A)') '(', Nb+1 ,'f24.6)'  
  write(filename, '(A,A)')  TRIM(label),'.dat'
  !
  write(stdout,'(A)') ' '
  write(stdout,'(A,A)') 'writing band structure on ', filename
  !
  open(2, file=filename, status='unknown')
  !
  x = 0.0d0 
  !
  write(2,formt) x, ek(1,:) 
  !
  do ik = 2, nkstot
    x = x + sqrt(dot_product(xk(:,ik)-xk(:,ik-1),xk(:,ik)-xk(:,ik-1)))
    write(2,formt) x, ek(ik,:)
  end do 
  !
  close(2)
  !
  return
  !
end subroutine print_bands
!----------------------------------------------------------------------------
END MODULE
!----------------------------------------------------------------------------
