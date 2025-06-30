!
! Copyright (C) 2008-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE read_ps_new ( psfile, upf, printout, ierr )
  !
  !! Read PP file "psfile" into structure "upf"
  !! Print some information if "printout" is true
  !! Return an error code in "ierr" as follows:
  !! ierr = 81  file cannot be opened (not found or not accessible)
  !! ierr >  0  error reading file
  !! ierr <= 0  file successfully read, file format is:
  !!    ierr =  0  UPF xml (unstable, experimental)
  !!    ierr = -1  UPF v.1 (deprecated)
  !!    ierr = -2  UPF v.2 (default choice)
  !!    ierr = -3  psml (experimental, Norm-Conserving only)
  !!    ierr = -4  old Vanderbilt formatted USPP (deprecated)
  !!    ierr = -5  old RRKJ3 USPP format (deprecated)
  !!    ierr = -6  old PWscf NCPP format (deorecated)
  !!    ierr = -7  Goedecker-Teter-Hutter NCPP
  !! Should be executed on a single processor
  !
  USE upf_io,             ONLY: stdout
  USE pseudo_types,       ONLY: pseudo_upf, reset_upf
  USE read_upf_v1_module, ONLY: read_upf_v1
  USE read_upf_new_module,ONLY: read_upf_new
  USE read_uspp_module,   ONLY: readvan, readrrkj
  USE read_psml_module,   ONLY: read_psml
  !USE upf_to_internal,    ONLY: set_upf_q
  !
  IMPLICIT NONE
  CHARACTER(LEN=256), INTENT(in) :: psfile
  LOGICAL, INTENT(IN)            :: printout
  TYPE(pseudo_upf), INTENT(out) :: upf
  INTEGER, INTENT(out) :: ierr
  !
  INTEGER :: iunps, l, lm3, lm4, lm5
  !
  !
  CALL reset_upf( upf )
  CALL read_upf_new( psfile, upf, ierr )
  !
  !! start reading - check  first if file is readable as xml file
  !! (ierr=0) or as UPF v.2 file (ierr=-2)
  !
  IF (ierr == 81 ) THEN
     WRITE (stdout, '("read_ps_new: file ",A," could not be opened")') trim(psfile)
     RETURN
  ELSE IF (ierr > 0 ) THEN
     !! file is not xml or UPF v.2 
     CALL  read_upf_v1 ( psfile, upf, ierr )
     !! try to read UPF v.1 file
     IF ( ierr == 0 ) ierr = -1
  END IF
  !
  IF ( ierr > 0 ) THEN
     !! file is not in any UPF format, try other formats 
     OPEN ( NEWUNIT=iunps, FILE=psfile, STATUS = 'old', &
          FORM = 'formatted', IOSTAT = ierr )
     IF (ierr > 0 ) GO TO 10
     !
     l = len_trim(psfile)
     lm3 = max(l-3,1)
     lm4 = max(l-4,1)
     lm5 = max(l-5,1)
     !! For unlikely short file names, avoid OOB error
     IF (psfile (lm4:l) == '.psml') THEN
        CALL read_psml (psfile, upf, ierr)
        IF ( ierr == 0 ) ierr = -3
     ELSE IF (psfile (lm3:l) =='.vdb' .OR. psfile (lm3:l) =='.van') THEN
        CALL readvan (iunps, upf, ierr)
        IF ( ierr == 0 ) ierr = -4
     ELSE IF (psfile (lm5:l) =='.RRKJ3') THEN
        CALL readrrkj (iunps, upf, ierr)
        IF ( ierr == 0 ) ierr = -5
     ELSE IF (psfile (lm3:l) =='.gth' .OR. psfile(lm3:l) == '.GTH' ) THEN
        !! FIXME: should follow the same logic of the other cases
        ierr = -7
     ELSE
        CALL read_ncpp (iunps, upf, ierr)
        IF ( ierr == 0 ) ierr = -6
     END IF
     !
     CLOSE (iunps)
10   IF ( ierr > 0 ) THEN
        WRITE (stdout, '("readpp: file ",A," could not be read")') trim(psfile)
        RETURN
     END IF
  END IF
  !
  !!! CALL set_upf_q (upf)
  !! reconstruct Q(r) if needed
   
  IF (printout) THEN
     !
     SELECT CASE(ierr)
        CASE(0)
           WRITE( stdout, "('file format is UPF xml (experimental)')") 
        CASE(-1)
           WRITE( stdout, "('file format is UPF v.1')")
        CASE(-2)
           WRITE( stdout, "('file format is UPF v.2')")
        CASE(-3)
           WRITE( stdout, "('file format is PSML (experimental)')")
        CASE(-4)
           WRITE( stdout, "('file format is Vanderbilt US PP')")
        CASE(-5)
           WRITE( stdout, "('file format is RRKJ3')")
        CASE(-6)
           WRITE( stdout, "('file format is old PWscf NC format')")
        CASE(-7)
           WRITE( stdout, "('file format is GTH (Goedecker-Teter-Hutter)')")
        CASE DEFAULT
           WRITE( stdout, "('file format could not be determined')")
     END SELECT
     !
  END IF
  !
END SUBROUTINE read_ps_new
