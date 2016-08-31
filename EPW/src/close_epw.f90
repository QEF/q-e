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
  USE phcom,     ONLY : iuwfc, iudwf, iudrhous, iudvkb3, fildrho, iudrho
  USE uspp,      ONLY : okvan      
  USE mp_global, ONLY : me_pool,root_pool
  !
  implicit none
  !
  CLOSE (unit = iuwfc, status = 'keep')
  CLOSE (unit = iudwf, status = 'keep')
  IF(okvan) CLOSE(unit = iudrhous, status = 'delete')
  IF(okvan) CLOSE (unit = iudvkb3, status = 'delete')
  IF (me_pool == root_pool ) THEN
    IF (fildrho.ne.' ') CLOSE (unit = iudrho, status = 'keep')
  ENDIF
  !
  END SUBROUTINE close_epw
