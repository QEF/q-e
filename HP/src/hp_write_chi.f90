!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------
SUBROUTINE hp_write_chi()
!------------------------------------------------------------
  !
  ! Write the nah_pert-th column (corresponding to the currently 
  ! perturbed atom) of chi0 and chi on file for restart.
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : ionode
  USE ions_base,     ONLY : nat
  USE io_files,      ONLY : prefix, tmp_dir
  USE ldaU_hp,       ONLY : nath_sc, nath, nah_pert, chi0, chi
  !
  IMPLICIT NONE
  INTEGER :: iunitchi ! unit number
  CHARACTER(len=50) :: filenamechi
  CHARACTER(len=256) :: tempfile
  CHARACTER(len=6), EXTERNAL :: int_to_char
  INTEGER, EXTERNAL :: find_free_unit
  !
  IF (ionode) THEN
     !
     iunitchi = find_free_unit()
     filenamechi = TRIM(prefix) // TRIM(".chi.pert_") // TRIM(int_to_char(nah_pert)) // TRIM(".dat")
     tempfile = TRIM(tmp_dir) // TRIM(filenamechi)
     !
     OPEN(iunitchi, file = tempfile, form = 'formatted', status = 'unknown')
     !
     CALL write_chi(chi0, 'chi0')
     CALL write_chi(chi,  'chi')
     !
     CLOSE(iunitchi)
     !
  ENDIF   
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE write_chi (chi_, name_)
  !  
  IMPLICIT NONE
  !
  CHARACTER(len=4), INTENT(IN) :: name_
  REAL(DP),         INTENT(IN) :: chi_(nath_sc, nat)
  INTEGER :: na
  !
  WRITE(iunitchi,'(6x,"row",2x,"column",2x,a4," matrix elements")') TRIM(name_)
  DO na = 1, nath_sc
     WRITE(iunitchi,'(1x,i7,2x,i4,3x,5f19.15)') na, nah_pert, chi_(na, nah_pert)
  ENDDO
  WRITE(iunitchi,*)
  !
  RETURN
  !
END SUBROUTINE write_chi
  !
END SUBROUTINE hp_write_chi
