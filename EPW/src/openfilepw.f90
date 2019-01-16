  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/openfilq - Quantum-ESPRESSO group                
  !-----------------------------------------------------------------------
  subroutine openfilepw
  !-----------------------------------------------------------------------
  !!
  !!     This subroutine opens all the files necessary for the EPW
  !!     calculation.
  !!
  !! RM - Nov/Dec 2014
  !! Imported the noncolinear case implemented by xlzhang
  !!
  !-----------------------------------------------------------------------
  USE io_files,         ONLY : prefix, diropn, seqopn
  USE units_lr,         ONLY : iuwfc, lrwfc
  USE wvfct,            ONLY : nbnd, npwx
  USE noncollin_module, ONLY : npol, nspin_mag
  USE units_ph,         ONLY : lrdrho, lrdrhous
  USE fft_base,         ONLY : dfftp
  USE uspp,             ONLY : okvan
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: find_free_unit
  !! integer variable for I/O control
  ! used for extracting the fildvscf0 directory
  !
  LOGICAL :: exst
  !! logical variable to check file existe
  !
  IF (len_trim(prefix) == 0) CALL errore('openfilepw', 'wrong prefix', 1)
  !
  !     The file with the wavefunctions
  !
  iuwfc = 20 
  lrwfc = 2 * nbnd * npwx * npol 
  CALL diropn(iuwfc, 'wfc', lrwfc, exst) 
  IF (.not. exst) CALL errore ('openfilepw','file '//TRIM( prefix )//'.wfc'//' not found',1)
  !
  !   file for setting unitary gauges of eigenstates
  !
  lrdrho = 2 * dfftp%nr1x *dfftp%nr2x *dfftp%nr3x * nspin_mag
  !
  ! RM - The file should be opened in the same place as fildvscf?
  ! RM - I need to see if we need this file.?
  !
  !   open a file with the static change of the charge
  !
  IF (okvan) THEN
    lrdrhous =  dfftp%nnr * nspin_mag
  ENDIF
  !
  RETURN
  !
  END SUBROUTINE openfilepw
