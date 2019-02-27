!
! Copyright (C) 2001-2018 Quantum_ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_write_chi_full
!-----------------------------------------------------------------------
  !
  ! This routine writes the response functions chi0 and chi to file.
  ! Number of columns corresponds to the number of Hubbard atoms in
  ! the primitive unit cell (nath), and the number of rows is equal to the 
  ! product of nath with the number of q points (i.e. nath_sc).
  !
  USE ions_base,    ONLY : nat
  USE io_files,     ONLY : prefix, tmp_dir
  USE ldaU_hp,      ONLY : nath_sc, nath, chi0, chi
  !
  IMPLICIT NONE
  !
  INTEGER :: na, nb
  INTEGER :: iunitchi ! unit number
  CHARACTER(len=50)  :: filenamechi
  CHARACTER(len=256) :: tempfile
  INTEGER, EXTERNAL  :: find_free_unit
  !
  ! Find and open unit
  !  
  iunitchi = find_free_unit()
  filenamechi = trim(prefix) // ".chi.dat"
  tempfile = trim(tmp_dir) // trim(filenamechi)
  !
  OPEN(iunitchi, file = tempfile, form = 'formatted', status = 'unknown')
  !
  WRITE(iunitchi,'(9x,"chi0 :")')
  DO na = 1, nath_sc
     WRITE(iunitchi,'(1x,5f19.15)') (chi0(na,nb), nb=1,nath) 
  ENDDO
  !
  WRITE(iunitchi,'(/9x,"chi :")')
  DO na = 1, nath_sc
     WRITE(iunitchi,'(1x,5f19.15)') (chi(na,nb), nb=1,nath) 
  ENDDO
  !
  CLOSE(iunitchi)
  !
  DEALLOCATE(chi0)
  DEALLOCATE(chi)
  !
  RETURN
  !
END SUBROUTINE hp_write_chi_full
