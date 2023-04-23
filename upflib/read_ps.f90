!
! Copyright (C) 2008-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE read_ps ( psfile, upf, ierr )
  !
  ! stripped-down version of readpp in Modules/read_pseudo.f90
  ! for serial execution only and for usage in conversion codes
  !
  USE read_upf_v1_module, ONLY: read_upf_v1
  USE read_upf_new_module,ONLY: read_upf_new
  USE pseudo_types,     ONLY: pseudo_upf
  USE upf_io ,          ONLY: stdout
  USE upf_to_internal,  ONLY: set_upf_q
  USE read_uspp_module, ONLY: readvan, readrrkj
  USE cpmd_module,      ONLY: read_cpmd
  USE m_gth,            ONLY: readgth
  USE fhi,              ONLY: readfhi
  !
  IMPLICIT NONE
  CHARACTER(LEN=256), INTENT(in) :: psfile
  TYPE(pseudo_upf), INTENT(out) :: upf
  INTEGER, INTENT(out) :: ierr
  !
  INTEGER :: iunps, l, lm3, lm4, lm5
  !
  !
  CALL read_upf_new( psfile, upf, ierr )
  !
  IF (ierr == 81 ) THEN
     WRITE (stdout, '("readpp: file ",A1,"  not readable")') trim(psfile)
     return
  ELSE IF (ierr > 0 ) THEN
     CALL  read_upf_v1 ( psfile, upf, ierr )
     !! try to read UPF v.1 file
     IF ( ierr == 0 ) WRITE( stdout, "('file type is UPF v.1')")
  ELSE IF ( ierr == 0 ) THEN
     WRITE( stdout, "('file type is xml')") 
  ELSE IF ( ierr ==-2 ) THEN
     ierr = 0
     WRITE( stdout, "('file type is UPF v.2')")
  END IF
  !
  IF ( ierr /= 0 ) THEN
     !
     l = LEN_TRIM (psfile)
     lm3 = max(l-3,1)
     lm4 = max(l-4,1)
     lm5 = max(l-5,1)
     !
     !     The type of the pseudopotential is determined by the file name
     !
     IF (psfile (lm4:l) =='.psml') THEN
        !
        ! Unlike following cases, file must be opened with xml tools
        !
        WRITE( stdout, "('file type is PSML (experimental)')")
        CALL read_psml (psfile, upf)
        !
     ELSE
        !
        OPEN ( NEWUNIT=iunps, FILE=psfile, STATUS = 'old', &
             FORM = 'formatted', IOSTAT = ierr )
        IF (ierr /= 0 ) RETURN
        !
        IF (psfile (lm3:l) =='.vdb' .OR. psfile (lm3:l) =='.van') THEN
           !
           WRITE( stdout, "('file type is Vanderbilt US PP')")
           CALL readvan (iunps, upf, ierr)
           !
        ELSE IF (psfile (lm3:l) =='.gth') THEN
           !
           WRITE( stdout, "('file type is GTH (analytical)')")
           CALL readgth (iunps, 1, upf, ierr)
           !
        ELSE IF (psfile (lm5:l) =='.RRKJ3') THEN
           !
           WRITE( stdout, "('file type is RRKJ3')")
           CALL readrrkj (iunps, upf, ierr)
           !
        ELSE IF (psfile (lm3:l) =='.cpi' .OR. psfile (l-3:l) =='.fhi') THEN
           !
           WRITE( stdout, "('file type is FHI .cpi or .fhi format')")
           CALL readfhi (iunps, upf)
           ierr = 0 
           !
        ELSE IF (psfile (lm4:l) =='.cpmd') THEN
           !
           WRITE( stdout, "('file type is CPMD NC format')")
           CALL read_cpmd (iunps, upf, ierr)
           !
        ELSE
           !
           WRITE( stdout, "('file type is old PWscf NC format')")
           CALL read_ncpp (iunps, upf, ierr)
           !
        ENDIF
        !
        ! end of reading
        CLOSE (iunps)
        ! Error return
        IF ( ierr /= 0 ) RETURN
        !
     END IF
     !
  ENDIF
  !
  ! reconstruct Q(r) if needed
  !
  CALL set_upf_q (upf)
  ! Normal return
  ierr = 0 
  !
END SUBROUTINE read_ps
