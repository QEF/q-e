!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine output_tau
  !-----------------------------------------------------------------------
  !
  use parameters, only: DP
  use constants, only: bohr_radius_angs
  use brilz, only: alat, at, bg
  use basis, only: nat, tau, atomic_positions, ityp, atm
  real (kind=DP), allocatable:: tau_out (:, :)
  !
  ! tau in output format
  !
  allocate (tau_out(3,nat))
  tau_out(:,:) = tau(:,:)
  SELECT CASE ( atomic_positions )
     !
     !  convert output atomic positions from internally used format
     !  (a0 units) to the same format used in input
     !
  CASE ('alat')
     !
     write (6, '("ATOMIC_POSITIONS (alat)")')
  CASE ('bohr')
     !
     write (6, '("ATOMIC_POSITIONS (bohr)")')
     tau_out(:,:) = tau_out(:,:)*alat
  CASE ('crystal')
     !
     write (6, '("ATOMIC_POSITIONS (crystal)")')
     call cryst_to_cart (nat, tau_out, at, 1)
  CASE ('angstrom')
     !
     write (6, '("ATOMIC_POSITIONS (angstrom)")')
     !
     tau_out(:,:) = tau_out(:,:)*alat*bohr_radius_angs
  CASE DEFAULT
     write (6, '("ATOMIC_POSITIONS")')
  END SELECT
  do na = 1, nat
     write (6,'(a3,3x,3f14.9)') atm(ityp(na)), (tau_out(i,na), i=1,3)
  enddo
  write (6, '(/)')
  return
end subroutine output_tau
