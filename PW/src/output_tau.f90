!
! Copyright (C) 2003-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE output_tau( print_lattice, print_final )
  !----------------------------------------------------------------------------
  !! Print cell parameters and atomic positions.
  !
  USE io_global,   ONLY : stdout
  USE kinds,       ONLY : DP
  USE constants,   ONLY : bohr_radius_angs, AVOGADRO
  USE cell_base,   ONLY : alat, at, bg, omega, cell_units
  USE ions_base,   ONLY : nat, tau, ityp, atm, if_pos, tau_format, amass
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: print_lattice
  !! if .TRUE. prints lattice parameters
  LOGICAL, INTENT(IN) :: print_final
  !! if .TRUE. prints the final coordinates
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: tau_out(:,:)
  INTEGER :: na, i, k
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
     WRITE( stdout, '(5x,a,1F12.5," g/cm^3")') &
                    "density = ", SUM( amass(ityp(1:nat)) )&
                                  /(omega*bohr_radius_angs**3 * 1.d-24)/AVOGADRO
     !
     ! ... convert output cell from internally used format
     ! ... (alat units) to the same format used in input
     SELECT CASE( cell_units )
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
  ENDIF
  !
  ! ... convert output atomic positions from internally used format
  ! ... (a0 units) to the same format used in input
  SELECT CASE( tau_format )
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
     CALL cryst_to_cart( nat, tau_out, bg, -1 )
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
     IF ( ALLOCATED( if_pos ) ) THEN
        IF ( ANY( if_pos(:,na) == 0 ) ) THEN
           WRITE( stdout,'(A3,3X,3F20.10,1X,3i4)') &
                        atm(ityp(na)), tau_out(:,na), if_pos(:,na)
        ELSE
           WRITE( stdout,'(A3,3X,3F20.10)') atm(ityp(na)), tau_out(:,na)
        END IF
     ELSE
        WRITE( stdout,'(A3,3X,3F20.10)') atm(ityp(na)), tau_out(:,na)
     ENDIF
     !
  ENDDO
  !
  IF ( print_final  ) WRITE( stdout, '("End final coordinates")') 
  WRITE( stdout, '(/)' )
  !
  DEALLOCATE( tau_out )
  !
  RETURN
  !
END SUBROUTINE output_tau


SUBROUTINE output_tau_rescaled(rescale)
  USE io_global,   ONLY : stdout
  USE kinds,       ONLY : DP
  USE ions_base,   ONLY : nat, tau, ityp, atm, if_pos, tau_format
  IMPLICIT NONE
  REAL(DP),INTENT(in) :: rescale
  INTEGER :: na
 
  IF(rescale==1._dp) RETURN
  IF(tau_format/="alat") RETURN

  WRITE( stdout, '(/"Atomic positions rescaled with new alat:")' )

  DO na = 1, nat
     !
     IF ( ALLOCATED( if_pos ) ) THEN
        IF ( ANY( if_pos(:,na) == 0 ) ) THEN
           WRITE( stdout,'(A3,3X,3F20.10,1X,3i4)') &
                        atm(ityp(na)), tau(:,na)*rescale, if_pos(:,na)
        ELSE
           WRITE( stdout,'(A3,3X,3F20.10)') atm(ityp(na)), tau(:,na)*rescale
        END IF
     ELSE
        WRITE( stdout,'(A3,3X,3F20.10)') atm(ityp(na)), tau(:,na)*rescale
     ENDIF
     !
  ENDDO

END SUBROUTINE
