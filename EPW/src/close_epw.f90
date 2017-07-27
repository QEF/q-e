  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Original code adapted from PH/close_phq - Quantum-ESPRESSO group
  ! 09/2009 This subroutine is functional and probably complete
  ! a few more files may be deleted to clean the working directory
  ! 
  !------------------------------------------------------------------
  SUBROUTINE close_epw
  !------------------------------------------------------------------
  !
  USE phcom,     ONLY : iuwfc, iudwf, fildrho, iudrho
  USE mp_global, ONLY : me_pool,root_pool
  USE io_epw,    ONLY : iunepmatwe
  USE epwcom,    ONLY : etf_mem
  !
  implicit none
  !
  IF (etf_mem == 1 .OR. etf_mem == 2) THEN
    CLOSE (unit = iunepmatwe, status = 'delete')
  ENDIF
  !
  CLOSE (unit = iuwfc, status = 'keep')
  CLOSE (unit = iudwf, status = 'keep')
  IF (me_pool == root_pool ) THEN
    IF (fildrho.ne.' ') CLOSE (unit = iudrho, status = 'keep')
  ENDIF
  !
  END SUBROUTINE close_epw
