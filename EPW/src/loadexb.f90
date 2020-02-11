  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------------
  SUBROUTINE loadexb()
  !----------------------------------------------------------------------------
  !!
  !! Load the information on the excluded bands in Wannierization
  !! 
  !-----------------------------------------------------------------------
  USE wvfct,         ONLY : nbnd
  USE epwcom,        ONLY : filukk
  USE io_var,        ONLY : iunukk
  USE elph2,         ONLY : ibndstart, ibndend, nbndep
  USE io_global,     ONLY : ionode_id, meta_ionode
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
  !
  IF (meta_ionode) THEN
    !
    OPEN(iunukk, FILE = filukk, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
    IF (ios /=0) CALL errore('loadexb', 'error opening ukk file', iunukk)
    !
    DO ibnd = 1, nbnd
      READ(iunukk, *) exbands(ibnd)
    ENDDO
    !
    ibndstart = 1
    DO ibnd = 1, nbnd
       IF (exbands(ibnd) == .FALSE.) THEN
          ibndstart = ibnd
          EXIT
       ENDIF
    ENDDO
    !
    ibndend= nbnd
    DO ibnd = nbnd, 1, -1
       IF (exbands(ibnd) == .FALSE.) THEN
          ibndend = ibnd
          EXIT
       ENDIF
    ENDDO
    nbndep = ibndend - ibndstart + 1
    !
    CLOSE(iunukk)
  ENDIF ! meta_ionode
  !
  CALL mp_bcast(ibndstart, ionode_id, inter_pool_comm)
  CALL mp_bcast(ibndend, ionode_id, inter_pool_comm)
  CALL mp_bcast(nbndep, ionode_id, inter_pool_comm)
  !
  RETURN
  !-----------------------------------------------------------------------
  END SUBROUTINE loadexb
  !-----------------------------------------------------------------------
