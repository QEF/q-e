  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  !
  ! Copyright (C) 2001 PWSCF group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  !
  ! Modified from stop_ph
  !
  !-----------------------------------------------------------------------
  SUBROUTINE stop_epw()
  !-----------------------------------------------------------------------
  !!
  !! Close all files and synchronize processes before stopping.
  !! Called at the end of the run
  !!
  USE mp,        ONLY : mp_end, mp_barrier
  USE mp_global, ONLY : inter_pool_comm, mp_global_end
  USE io_global, ONLY : stdout
  USE printing,  ONLY : print_clock_epw
  USE epwcom,    ONLY : eliashberg, plselfen, specfun_pl, scattering, iterative_bte
  USE elph2,     ONLY : adapt_smearing
  ! 
  IMPLICIT NONE
  !
  CALL print_clock_epw()
  ! 
  WRITE(stdout, '(a)') "                                                                                          "
  WRITE(stdout, '(a)') " Please consider citing:                                                                  "
  WRITE(stdout, '(a)') " S. Ponce, E. R. Margine, C. Verdi and F. Giustino, Comput. Phys. Commun. 209, 116 (2016) " 
  WRITE(stdout, '(a)') "                                                                                          "
  ! 
  ! Specific functionalities
  IF (eliashberg .OR. plselfen .OR. specfun_pl .OR. scattering .OR. iterative_bte .OR. adapt_smearing) THEN
    WRITE(stdout, '(a)') " In addition, since you have used the following functionalities, please cite:             "
  ENDIF
  ! 
  ! Eliashberg superconductivity
  IF (eliashberg) THEN
    WRITE(stdout, '(a)') "   eliashberg :: E. R. Margine and F. Giustino, Phys. Rev. B 87, 024505 (2013)              "
  ENDIF  
  ! 
  ! Plasmons
  IF (plselfen .OR. specfun_pl) THEN
    WRITE(stdout, '(a)') &
  "   plselfen or specfun_pl :: F. Caruso, C. Verdi, S. Ponce and F. Giustino, Phys. Rev. B 97, 165113 (2018) "
  ENDIF 
  ! 
  ! Transport module 
  IF (scattering) THEN 
    WRITE(stdout, '(a)') "   scattering :: S. Ponce, E. R. Margine and F. Giustino, Phys. Rev. B 97, 121201 (2018)     "
  ENDIF
  IF (iterative_bte) THEN
    WRITE(stdout, '(a)') "   iterative_bte :: S. Ponce, E. R. Margine and F. Giustino, Phys. Rev. B 97, 121201 (2018)  "
    WRITE(stdout, '(a)') "   iterative_bte :: F. Macheda and N. Bonini, Phys. Rev. B 98, 201201 (2018)                 "
  ENDIF  
  ! 
  ! Improvements
  IF (adapt_smearing) THEN
    WRITE(stdout, '(a)') "   adapt_smearing :: F. Macheda and N. Bonini, Phys. Rev. B 98, 201201 (2018)                " 
  ENDIF
  !  
  CALL mp_end(inter_pool_comm)
  CALL mp_global_end()
  ! 
  STOP
  ! 
  RETURN
  !-----------------------------------------------------------------------
  END SUBROUTINE stop_epw
  !-----------------------------------------------------------------------
