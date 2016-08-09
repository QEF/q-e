  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Roxana Margine
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------
  SUBROUTINE deallocate_eliashberg
  !----------------------------------------------------------------------
  !!
  !!  deallocates the variables allocated by allocate_eliashberg
  !!
  USE epwcom,        ONLY : liso, laniso, limag
  USE eliashbergcom, ONLY : wsph, estemp, gap, wsph, agap
  !
  IMPLICIT NONE
  !
  IF( ALLOCATED(wsph) )  DEALLOCATE(wsph)
  IF( ALLOCATED(estemp)) DEALLOCATE(estemp)
  !
  IF ( liso ) THEN 
     IF ( limag ) THEN
        CALL deallocate_eliashberg_iso_iaxis
     ENDIF
     CALL deallocate_eliashberg_iso_raxis
  ENDIF
  !
  IF ( laniso ) THEN
     IF ( limag ) THEN
        CALL deallocate_eliashberg_aniso_iaxis
     ENDIF
     CALL deallocate_eliashberg_aniso_raxis
     IF( ALLOCATED(gap))  DEALLOCATE(gap)
     IF( ALLOCATED(Agap)) DEALLOCATE(Agap)
  ENDIF
  !
  CALL deallocate_elphon
  !
  RETURN
  !
  END SUBROUTINE deallocate_eliashberg
  !
  !----------------------------------------------------------------------
  SUBROUTINE deallocate_eliashberg_iso_iaxis
  !----------------------------------------------------------------------
  !!
  !!  deallocates the variables allocated by allocate_eliashberg_iso_iaxis
  !!
  !----------------------------------------------------------------------
  !
  USE eliashbergcom
  !
  IMPLICIT NONE
  !
  IF( ALLOCATED(wsi) )     DEALLOCATE(wsi)
  IF( ALLOCATED(Deltai) )  DEALLOCATE(Deltai)
  IF( ALLOCATED(Deltaip) ) DEALLOCATE(Deltaip)
  IF( ALLOCATED(Znormi) )  DEALLOCATE(Znormi)
  IF( ALLOCATED(NZnormi) ) DEALLOCATE(NZnormi)
  IF( ALLOCATED(Keri) )    DEALLOCATE(Keri)
  !
  RETURN
  !
  END SUBROUTINE deallocate_eliashberg_iso_iaxis
  !                         
  !----------------------------------------------------------------------
  SUBROUTINE deallocate_eliashberg_iso_raxis
  !----------------------------------------------------------------------
  !!
  !!  deallocates the variables allocated by allocate_eliashberg_iso_raxis
  !!
  USE epwcom, ONLY : lreal, limag, lacon
  USE eliashbergcom
  !
  IMPLICIT NONE
  !
  IF( ALLOCATED(ws) )    DEALLOCATE(ws)
  !
  IF( ALLOCATED(Delta))  DEALLOCATE(Delta)
  IF( ALLOCATED(Deltap)) DEALLOCATE(Deltap)
  IF( ALLOCATED(Znorm))  DEALLOCATE(Znorm)
  IF( ALLOCATED(Znormp)) DEALLOCATE(Znormp)
  IF( ALLOCATED(gap))    DEALLOCATE(gap)
  !
  IF ( lreal ) THEN
     IF( ALLOCATED(dws) )   DEALLOCATE(dws)
     IF( ALLOCATED(fdwp) )  DEALLOCATE(fdwp)
     IF( ALLOCATED(bewph) ) DEALLOCATE(bewph)
     IF( ALLOCATED(Kp))     DEALLOCATE(Kp)    
     IF( ALLOCATED(Km))     DEALLOCATE(Km) 
  ENDIF
  !
  IF ( limag .AND. lacon ) THEN
     IF( ALLOCATED(Gp))     DEALLOCATE(Gp)
     IF( ALLOCATED(Gm))     DEALLOCATE(Gm)
     IF( ALLOCATED(Dsumi) ) DEALLOCATE(Dsumi)
     IF( ALLOCATED(Zsumi) ) DEALLOCATE(Zsumi)
  ENDIF
  !
  RETURN
  !
  END SUBROUTINE deallocate_eliashberg_iso_raxis
  !
  !----------------------------------------------------------------------
  SUBROUTINE deallocate_eliashberg_aniso_iaxis
  !----------------------------------------------------------------------
  !!
  !!  deallocates the variables allocated by allocate_eliashberg_aniso_iaxis
  !!
  USE eliashbergcom
  !
  IMPLICIT NONE
  !
  IF( ALLOCATED(wsi) )      DEALLOCATE(wsi)
  !
  IF( ALLOCATED(Deltai) )   DEALLOCATE(Deltai)
  IF( ALLOCATED(Znormi) )   DEALLOCATE(Znormi)
  ! 
  IF( ALLOCATED(ADeltai) )  DEALLOCATE(ADeltai)
  IF( ALLOCATED(ADeltaip) ) DEALLOCATE(ADeltaip)
  IF( ALLOCATED(AZnormi) )  DEALLOCATE(AZnormi)
  IF( ALLOCATED(NAZnormi) ) DEALLOCATE(NAZnormi)
  !
  RETURN
  !
  END SUBROUTINE deallocate_eliashberg_aniso_iaxis
  !                                
  !----------------------------------------------------------------------
  SUBROUTINE deallocate_eliashberg_aniso_raxis
  !----------------------------------------------------------------------
  !!
  !!  deallocates the variables allocated by allocate_eliashberg_aniso_raxis
  !!
  USE eliashbergcom
  !
  IMPLICIT NONE
  !
  IF( ALLOCATED(ws))       DEALLOCATE(ws) 
  !
  IF( ALLOCATED(Delta))    DEALLOCATE(Delta)
  IF( ALLOCATED(Znorm))    DEALLOCATE(Znorm)
  !
  IF( ALLOCATED(ADelta) )  DEALLOCATE(ADelta)
  IF( ALLOCATED(ADeltap) ) DEALLOCATE(ADeltap)
  IF( ALLOCATED(AZnorm) )  DEALLOCATE(AZnorm)
  IF( ALLOCATED(AZnormp) ) DEALLOCATE(AZnormp)
  !
  RETURN
  !
  END SUBROUTINE deallocate_eliashberg_aniso_raxis
  !
  !----------------------------------------------------------------------
  SUBROUTINE deallocate_elphon
  !----------------------------------------------------------------------
  !!
  !!  deallocates the variables allocated by electron-phonon
  !!
  USE elph2,         ONLY : wf, wqf
  USE eliashbergcom, ONLY : ekfs, xkfs, wkfs, xkff, g2, a2f_iso, w0g, ixkff, ixkqf, ixqfs, nqfs
  !
  IMPLICIT NONE
  !
  IF( ALLOCATED(wf) )      DEALLOCATE(wf)
  IF( ALLOCATED(wqf) )     DEALLOCATE(wqf)
  IF( ALLOCATED(ekfs) )    DEALLOCATE(ekfs)
  IF( ALLOCATED(xkfs) )    DEALLOCATE(xkfs)
  IF( ALLOCATED(wkfs) )    DEALLOCATE(wkfs)
  IF( ALLOCATED(xkff) )    DEALLOCATE(xkff)
  IF( ALLOCATED(g2) )      DEALLOCATE(g2)
  IF( ALLOCATED(a2f_iso) ) DEALLOCATE(a2f_iso)
  IF( ALLOCATED(w0g) )     DEALLOCATE(w0g)
  IF( ALLOCATED(ixkff) )   DEALLOCATE(ixkff)
  IF( ALLOCATED(ixkqf) )   DEALLOCATE(ixkqf)
  IF( ALLOCATED(ixqfs) )   DEALLOCATE(ixqfs)
  IF( ALLOCATED(nqfs) )    DEALLOCATE(nqfs)
  !
  RETURN
  !
  END SUBROUTINE deallocate_elphon
  !                                                                            
  !----------------------------------------------------------------------
