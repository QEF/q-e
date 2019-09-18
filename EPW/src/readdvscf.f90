  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-------------------------------------------------------------
  SUBROUTINE readdvscf(dvscf, recn, iq, nqc)
  !-------------------------------------------------------------
  !!
  !! Open dvscf files as direct access, read, and close again
  !!
  !! SP - Nov 2017
  !! Replaced fstat by Fortran instric function inquire. 
  !! 
  !! RM - Nov/Dec 2014
  !! Imported the noncolinear case implemented by xlzhang
  !!
  !-------------------------------------------------------------
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : prefix
  USE units_ph,  ONLY : lrdrho
  USE fft_base,  ONLY : dfftp
  !USE pwcom
  USE epwcom,    ONLY : dvscf_dir
  USE io_epw,    ONLY : iudvscf
  USE low_lvl,   ONLY : set_ndnmbr 
  USE noncollin_module, ONLY : nspin_mag
  !
  IMPLICIT NONE
  ! 
  INTEGER, INTENT(in) :: recn
  !! perturbation number
  INTEGER, INTENT(in) :: iq
  !! the current q-point
  INTEGER, INTENT(in) :: nqc
  !! the total number of q-points in the list
  COMPLEX(KIND = DP), INTENT(out) :: dvscf(dfftp%nnr, nspin_mag) 
  !! dVscf potential is read from file
  !
  ! Local variables
  !
  CHARACTER(LEN = 256) :: tempfile
  !! Temp file 
  CHARACTER(LEN = 3) :: filelab
  !! File number
  INTEGER :: unf_recl
  !! Rcl unit
  INTEGER :: ios
  !! Error number
  INTEGER(KIND = 8) :: mult_unit
  !! Record length
  INTEGER(KIND = 8) :: file_size
  !! File size
  REAL(KIND = DP) :: dummy
  !! Dummy variable 
  !
  !  the call to set_ndnmbr is just a trick to get quickly
  !  a file label by exploiting an existing subroutine
  !  (if you look at the sub you will find that the original 
  !  purpose was for pools and nodes)
  !   
  CALL set_ndnmbr(0, iq, 1, nqc, filelab)
  tempfile = TRIM(dvscf_dir) // TRIM(prefix) // '.dvscf_q' // filelab
  INQUIRE(IOLENGTH = unf_recl) dummy 
  unf_recl = unf_recl  * lrdrho
  mult_unit = unf_recl
  mult_unit = recn * mult_unit
  !
  !  open the dvscf file, read and close
  !
  OPEN(iudvscf, FILE = tempfile, FORM = 'unformatted', &
       ACCESS = 'direct', IOSTAT = ios, RECL = unf_recl, STATUS = 'old')
  IF (ios /= 0) CALL errore('readdvscf', 'error opening ' // tempfile, iudvscf)
  !
  ! check that the binary file is long enough
  INQUIRE(FILE = tempfile, SIZE = file_size)
  IF (mult_unit > file_size) CALL errore('readdvscf', &
       TRIM(tempfile) //' too short, check ecut', iudvscf)
  !
  READ(iudvscf, REC = recn) dvscf
  CLOSE(iudvscf, STATUS = 'keep')
  !
  RETURN
  !
  !-------------------------------------------------------------
  END SUBROUTINE readdvscf
  !-------------------------------------------------------------
