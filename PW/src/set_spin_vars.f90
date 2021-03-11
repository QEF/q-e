!
! Copyright (C) 2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE set_spin_vars( lsda, noncolin, domag, &
     npol, nspin, nspin_lsda, nspin_mag, nspin_gga, current_spin )
  !------------------------------------------------------------------------
  !
  !  Set various spin-related variables
  !
  LOGICAL, INTENT(IN)  :: lsda, noncolin, domag
  INTEGER, INTENT(OUT) :: npol, nspin, nspin_lsda, nspin_mag, nspin_gga
  INTEGER, INTENT(OUT) :: current_spin
  !
  IF ( lsda ) THEN
     !
     ! ... wavefunctions have up and down spin
     !
     npol  = 1
     nspin = 2
     nspin_mag  = 2
     nspin_lsda = 2
     nspin_gga  = 2
     current_spin = -1
  ELSE IF ( noncolin ) THEN
     !
     ! ... wavefunctions are spinors with 2 components
     !
     npol  = 2
     nspin = 4
     nspin_lsda = 1
     IF (domag) THEN
        nspin_gga = 2
        nspin_mag = 4
     ELSE
        nspin_gga = 1
        nspin_mag = 1
     ENDIF
     current_spin = 1
  ELSE
     !
     ! ... wavefunctions are scalars
     !
     npol  = 1
     nspin = 1
     nspin_mag  = 1
     nspin_lsda = 1
     nspin_gga  = 1
     current_spin = 1
  END IF
  !
END SUBROUTINE set_spin_vars
