  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Roxana Margine
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE eliashberg_eqs()
  !-----------------------------------------------------------------------
  !!
  !! This is the main driver for solving the Eliashberg equations 
  !!
  USE io_global,         ONLY : stdout 
  USE epwcom,            ONLY : liso, fila2f, gap_edge, lreal, limag, laniso 
  USE eliashbergcom,     ONLY : gap0
  USE supercond, ONLY : eliashberg_init, evaluate_a2f_lambda, & 
                                estimate_tc_gap, deallocate_eliashberg_aniso, & 
                                deallocate_eliashberg_elphon
  USE io_eliashberg,     ONLY : read_a2f, read_frequencies, read_eigenvalues, & 
                                read_ephmat, read_kqmap  
  USE supercond_iso,     ONLY : eliashberg_iso_iaxis, eliashberg_iso_raxis
  USE supercond_aniso,   ONLY : eliashberg_aniso_iaxis
  !
  IMPLICIT NONE
  !
  CALL start_clock('ELIASHBERG')
  !
  IF (liso) THEN
    WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
    WRITE(stdout, '(5x, "Solve isotropic Eliashberg equations")')
    WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
    CALL eliashberg_init()
    IF (fila2f == ' ') THEN
      CALL read_frequencies()
      CALL read_eigenvalues()
      CALL read_kqmap()
      CALL read_ephmat()
      CALL evaluate_a2f_lambda()
      CALL deallocate_eliashberg_elphon()
    ENDIF
    !
    CALL read_a2f()
    CALL estimate_tc_gap()
    IF (gap_edge > 0.d0) THEN
      gap0 = gap_edge
    ENDIF
    IF (lreal) CALL eliashberg_iso_raxis()
    IF (limag) CALL eliashberg_iso_iaxis()
  ENDIF
  !
  IF (laniso) THEN
    WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
    WRITE(stdout, '(5x, "Solve anisotropic Eliashberg equations")')
    WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
    CALL eliashberg_init()
    CALL read_frequencies()
    CALL read_eigenvalues()
    CALL read_kqmap()
    CALL read_ephmat()
    CALL evaluate_a2f_lambda()
    CALL estimate_tc_gap()
    IF (gap_edge > 0.d0) THEN 
      gap0 = gap_edge
    ENDIF
    CALL eliashberg_aniso_iaxis()
  ENDIF
  !
  IF (.NOT. liso .AND. .NOT. laniso) THEN 
    WRITE(stdout, '(/5x, a)') REPEAT('=', 67)
    WRITE(stdout, '(5x, "Calculate Eliashberg spectral function")')
    WRITE(stdout, '(5x, a/)') REPEAT('=', 67)
    CALL eliashberg_init()
    CALL read_frequencies()
    CALL read_eigenvalues()
    CALL read_kqmap()
    CALL read_ephmat()
    CALL evaluate_a2f_lambda()
    CALL estimate_tc_gap() 
    CALL deallocate_eliashberg_aniso()
  ENDIF
  !
  CALL stop_clock('ELIASHBERG')
  !
  RETURN
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE eliashberg_eqs
  !-----------------------------------------------------------------------
