  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------------
  SUBROUTINE loadbm()
  !----------------------------------------------------------------------------
  !!
  !! Load the information on the band manifold determined in Wannierization step
  !! 
  !-----------------------------------------------------------------------
  USE wvfct,         ONLY : nbnd
  USE epwcom,        ONLY : filukk, nbndskip
  USE io_var,        ONLY : iunukk
  USE elph2,         ONLY : ibndstart, ibndend, nbndep
  USE io_global,     ONLY : ionode_id, meta_ionode, stdout
  USE mp_global,     ONLY : inter_pool_comm
  USE mp,            ONLY : mp_bcast
  !
  IMPLICIT NONE
  ! 
  ! Local variables 
  LOGICAL :: exbands(nbnd)
  !! Excluded bands in Wannierization step
  INTEGER :: ibnd
  !! Counter on band index
  INTEGER :: ios
  !! INTEGER variable for I/O control
  INTEGER :: ierr
  !! Error status
  CHARACTER*100 dummy
  LOGICAL lll
  !
  IF (meta_ionode) THEN
    !
    OPEN(iunukk, FILE = filukk, STATUS = 'old', FORM = 'formatted', POSITION='REWIND', IOSTAT = ios)
    IF (ios /=0) CALL errore('loadexb', 'error opening ukk file', iunukk)
    !
    READ(iunukk, *) ibndstart, ibndend
    nbndep = ibndend - ibndstart + 1
    nbndskip = ibndstart - 1
    !
    CLOSE(iunukk)
  ENDIF ! meta_ionode
  !
  CALL mp_bcast(ibndstart, ionode_id, inter_pool_comm)
  CALL mp_bcast(ibndend, ionode_id, inter_pool_comm)
  CALL mp_bcast(nbndep, ionode_id, inter_pool_comm)
  CALL mp_bcast(nbndskip, ionode_id, inter_pool_comm)
  !
  RETURN
  !-----------------------------------------------------------------------
  END SUBROUTINE loadbm
  !-----------------------------------------------------------------------
