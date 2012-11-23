!
! Copyright (C) 2003-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE output_tau( print_lattice, print_final  )
  !----------------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE constants, ONLY : bohr_radius_angs
  USE cell_base, ONLY : alat, at, bg, omega, cell_units
  USE ions_base, ONLY : nat, tau, ityp, atm, if_pos, tau_format
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN)         :: print_lattice, print_final
  REAL (DP), ALLOCATABLE :: tau_out(:,:)
  INTEGER                     :: na, i, k
  !
  !
  ! ... tau in output format
  !
  ALLOCATE( tau_out(3,nat) )
  !
  tau_out(:,:) = tau(:,:)
  !
  ! ... print cell parameters if required
  !
  IF ( print_final  ) WRITE( stdout, '("Begin final coordinates")') 
  IF ( print_lattice ) THEN
     !
     WRITE( stdout, '(5x,a,1F12.5," a.u.^3 ( ",1F11.5," Ang^3 )")') &
                    "new unit-cell volume = ",omega, omega*bohr_radius_angs**3 

     SELECT CASE (cell_units)
     !
     ! ... convert output cell from internally used format
     ! ... (alat units) to the same format used in input
     !
     CASE( 'alat' )
        WRITE( stdout, '(/"CELL_PARAMETERS (alat=",f12.8,")")') alat 
        WRITE( stdout, '(3F14.9)') ( ( at(i,k), i = 1, 3), k = 1, 3 )
     CASE( 'bohr' )
        WRITE( stdout, '(/"CELL_PARAMETERS (bohr)")')  
        WRITE( stdout, '(3F14.9)') ( ( at(i,k) * alat, i = 1, 3), k = 1, 3 )
     CASE( 'angstrom' )
        WRITE( stdout, '(/"CELL_PARAMETERS (angstrom)")') 
        WRITE( stdout, '(3F14.9)') &
             ( ( at(i,k) * alat * bohr_radius_angs, i = 1, 3), k = 1, 3 )
     CASE DEFAULT
        WRITE( stdout, '(/"CELL_PARAMETERS (alat=",f12.8,")")') alat 
        WRITE( stdout, '(3F14.9)') ( ( at(i,k), i = 1, 3), k = 1, 3 )
     END SELECT
     !
  END IF
  !
  SELECT CASE( tau_format )
     !
     ! ... convert output atomic positions from internally used format
     ! ... (a0 units) to the same format used in input
     !
  CASE( 'alat' )
     !
     WRITE( stdout, '(/"ATOMIC_POSITIONS (alat)")' )
     !
  CASE( 'bohr' )
     !
     WRITE( stdout, '(/"ATOMIC_POSITIONS (bohr)")' )
     tau_out(:,:) = tau_out(:,:) * alat
     !
  CASE( 'crystal' )
     !
     WRITE( stdout, '(/"ATOMIC_POSITIONS (crystal)")' )
     !
     call cryst_to_cart( nat, tau_out, bg, -1 )
     !
  CASE( 'angstrom' )
     !
     WRITE( stdout, '(/"ATOMIC_POSITIONS (angstrom)")' )
     !
     tau_out(:,:) = tau_out(:,:) * alat * bohr_radius_angs
     !
  CASE DEFAULT
     !
     WRITE( stdout, '(/"ATOMIC_POSITIONS")' )
     !
  END SELECT
  !
  DO na = 1, nat
     !
     IF ( ANY( if_pos(:,na) == 0 ) ) THEN
        WRITE( stdout,'(A3,3X,3F14.9,1X,3i4)') &
                        atm(ityp(na)), tau_out(:,na), if_pos(:,na)
     ELSE
        WRITE( stdout,'(A3,3X,3F14.9)') &
                        atm(ityp(na)), tau_out(:,na)
     END IF
     !
  END DO
  !
  IF ( print_final  ) WRITE( stdout, '("End final coordinates")') 
  WRITE( stdout, '(/)' )
  !
  DEALLOCATE( tau_out )
  !
  RETURN
  !
END SUBROUTINE output_tau
