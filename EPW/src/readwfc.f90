  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !------------------------------------------------------------
  SUBROUTINE readwfc(ipool, recn, evc0)
  !------------------------------------------------------------
  !!
  !! Open wfc files as direct access, read, and close again
  !!
  !! RM - Nov/Dec 2014
  !! Imported the noncolinear case implemented by xlzhang
  !!
  !
  USE kinds,    ONLY : DP
  USE io_files, ONLY : prefix, tmp_dir
  USE units_lr, ONLY : lrwfc, iuwfc
  USE wvfct,    ONLY : npwx
  USE pwcom,    ONLY : nbnd
  USE low_lvl,  ONLY : set_ndnmbr
  USE noncollin_module, ONLY : npol
  USE mp_global,        ONLY : nproc_pool, me_pool, npool
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: recn
  !! kpoint number
  INTEGER, INTENT(in) :: ipool
  !! poolfile number to be read (not used in serial case)
  COMPLEX(KIND = DP), INTENT(out) :: evc0(npwx * npol, nbnd)
  !! wavefunction is read from file
  !
  ! Local variables
  CHARACTER(LEN = 256) :: tempfile
  !! Temp file 
  CHARACTER(LEN = 3) :: nd_nmbr0
  !! File number
  INTEGER :: unf_recl
  !! Rcl unit
  INTEGER :: ios
  !! Error number
  REAL(KIND = DP) :: dummy
  !! Dummy variable 
  !
  ! Open the wfc file, read and close
  CALL set_ndnmbr(ipool, me_pool, nproc_pool, npool, nd_nmbr0)
  !
#if defined(__MPI)
  tempfile = TRIM(tmp_dir) // TRIM(prefix) // '.wfc' // nd_nmbr0
# else
  tempfile = TRIM(tmp_dir) // TRIM(prefix) // '.wfc'
#endif
  INQUIRE(IOLENGTH = unf_recl) dummy 
  unf_recl = unf_recl * lrwfc
  !
  OPEN(iuwfc, FILE = tempfile, FORM = 'unformatted', ACCESS = 'direct', IOSTAT = ios, RECL = unf_recl)
  IF (ios /= 0) CALL errore('readwfc', 'error opening wfc file', iuwfc)
  READ(iuwfc, REC = recn) evc0
  CLOSE(iuwfc, STATUS = 'keep')
  !
  RETURN
  !
  !------------------------------------------------------------
  END SUBROUTINE readwfc
  !------------------------------------------------------------
