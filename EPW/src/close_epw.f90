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
  USE io_files, ONLY: iunigk
  USE io_epw,   ONLY: iunepmatf, iuetf
  USE phcom,    ONLY: iuwfc, iudwf, iudrhous, iudvkb3, fildrho, iudrho
  USE epwcom,   ONLY: elinterp, iuncuf
  USE uspp,     ONLY: okvan      
#ifdef __PARA
  USE mp_global, ONLY : me_pool,root_pool
#endif
  !
  implicit none
  !
  CLOSE (unit = iuwfc, status = 'keep')
  CLOSE (unit = iudwf, status = 'keep')
  IF(okvan) CLOSE(unit = iudrhous, status = 'delete')
  IF(okvan) CLOSE (unit = iudvkb3, status = 'delete')
#ifdef __PARA
   IF (me_pool  /= root_pool ) go to 100  
#endif
  IF (fildrho.ne.' ') CLOSE (unit = iudrho, status = 'keep')
#ifdef __PARA
100 continue
#endif
  
  IF (elinterp) then
     !
     !  the temporary storage for Wannier interpolation
     !
     CLOSE (unit = iuncuf,     status = 'delete')
     CLOSE (unit = iunepmatf,  status = 'delete')
     CLOSE (unit = iuetf,      status = 'delete')
!     CLOSE (unit = iunepmatwe, status = 'delete')
!     CLOSE (unit = iunepmatwp, status = 'delete')
     !
  ENDIF
     !
  CLOSE (unit = iunigk, status = 'delete')
  !
  END SUBROUTINE close_epw
