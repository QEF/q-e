!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine check_atoms (nvec, vec, trmat)
  !-----------------------------------------------------------------------
  !
  !     This routine tests that the atomic coordinates (or k-points)
  !     are different and not related by a lattice translation
  !
  !
  USE kinds
  implicit none
  !
  real(DP), parameter :: accep=1.d-5
  !
  integer, intent(in) :: nvec
  ! nvec : number of atomic positions (or k-points)
  real(DP), intent(in) :: vec (3, nvec), trmat (3, 3)
  ! vec  : cartesian coordinates of atomic positions (or k-points)
  ! trmat: transformation matrix to crystal axis
  !        ( = bg , basis of the real-space lattice, for atoms
  !          = at , basis of the rec.-space lattice, for k-points )
  !
  integer :: nv1, nv2
  real(DP), allocatable :: vaux(:,:)
  real(DP) :: zero (3) = 0.0_dp
  character(len=80) :: message
  logical, external :: eqvect
  !
  !   Copy input positions and transform them to crystal units
  !
  allocate ( vaux(3,nvec) )
  vaux = vec
  call cryst_to_cart ( nvec, vaux, trmat, -1)
  !
  !   Test that all the atomic positions (or k-points) are different
  !
  do nv1 = 1, nvec-1
     do nv2 = nv1+1, nvec
        if ( eqvect ( vaux (1,nv1), vaux (1,nv2), zero, accep ) ) then
           zero(:) = vaux (:,nv1)-vaux(:,nv2)
           IF ( ALL ( abs(zero) < accep ) ) THEN
              WRITE (message,'("atoms #",i4," and #",i4," overlap!")') nv1, nv2
           ELSE
              WRITE (message,'("atoms #",i4," and #",i4," differ by lattice &
                    &vector (",i2,",",i2,",",i2,") in crystal axis")') &
                    nv1, nv2, NINT(zero)
           ENDIF
           CALL errore ( 'check_atoms', TRIM(message), 1)
        end if
     enddo
  enddo
  !
  deallocate(vaux)
  return
end subroutine check_atoms

