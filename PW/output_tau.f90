!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine output_tau (plot_lattice)
  !-----------------------------------------------------------------------
  !
  USE io_global,  ONLY : stdout
  USE kinds, only: DP
  use constants, only: bohr_radius_angs
  use brilz, only: alat, at, bg
  use basis, only: nat, tau, atomic_positions, ityp, atm
  implicit none
  LOGICAL, INTENT(IN):: plot_lattice
  real (kind=DP), allocatable:: tau_out (:, :)
  integer :: na, i, k
  !
  ! tau in output format
  !
  allocate (tau_out(3,nat))
  tau_out(:,:) = tau(:,:)
  !
  ! print cell parameters if required
  !
  IF (plot_lattice) THEN
     WRITE( stdout, '("CELL_PARAMETERS (alat)")') 
     WRITE( stdout, '(3f14.9)') ( (at (i, k) , i = 1, 3) , k = 1, 3)
  END IF
  !
  SELECT CASE ( atomic_positions )
     !
     !  convert output atomic positions from internally used format
     !  (a0 units) to the same format used in input
     !
  CASE ('alat')
     !
     WRITE( stdout, '("ATOMIC_POSITIONS (alat)")')
  CASE ('bohr')
     !
     WRITE( stdout, '("ATOMIC_POSITIONS (bohr)")')
     tau_out(:,:) = tau_out(:,:)*alat
  CASE ('crystal')
     !
     WRITE( stdout, '("ATOMIC_POSITIONS (crystal)")')
     call cryst_to_cart (nat, tau_out, bg,-1)
  CASE ('angstrom')
     !
     WRITE( stdout, '("ATOMIC_POSITIONS (angstrom)")')
     !
     tau_out(:,:) = tau_out(:,:)*alat*bohr_radius_angs
  CASE DEFAULT
     WRITE( stdout, '("ATOMIC_POSITIONS")')
  END SELECT
  !
  do na = 1, nat
     WRITE( stdout,'(a3,3x,3f14.9)') atm(ityp(na)), (tau_out(i,na), i=1,3)
  enddo
  WRITE( stdout, '(/)')
  deallocate(tau_out)
  return
end subroutine output_tau
