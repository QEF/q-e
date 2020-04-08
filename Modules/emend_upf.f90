! Copyright (C) 2017 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!=----------------------------------------------------------------------------=!
MODULE emend_upf_module 
!=----------------------------------------------------------------------------=!
  !! author: Pietro Delugas
  !! Contains utility to make the old UPF format readable by FoX

PRIVATE 
PUBLIC make_emended_upf_copy

CONTAINS 
FUNCTION  make_emended_upf_copy( filename, tempname)  RESULT(xml_check)
  !! author: Pietro Delugas
  !! Utility to make the old UPF format readable by FoX
  !! Replaces "&" with "&amp;" in file "filename", writes to file "tempname"
  !
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)      :: filename, tempname
  LOGICAL                          :: xml_check
  !
  INTEGER                          :: iun_source, iun_dest, ierr 
  INTEGER,EXTERNAL                 :: find_free_unit 
  LOGICAL                          :: icopy = .FALSE.
  CHARACTER(LEN=1024)              :: line 
  ! 
  iun_source = find_free_unit()
  OPEN (UNIT = iun_source, FILE = TRIM(filename), STATUS = 'old', &
       ACTION = 'read', FORM='formatted', iostat=ierr)
  IF ( ierr /= 0 ) CALL errore ("make_emended_upf", &
          "error opening file " // TRIM (filename),abs(ierr))
  READ ( iun_source, "(a)", IOSTAT = ierr ) line
     IF ( ierr < 0 ) CALL errore ("make_emended_upf", &
             TRIM (filename) // " is empty",abs(ierr))
     IF (INDEX(line, '<?xml') == 0 .AND. INDEX(line,'<UPF') == 0) THEN
        xml_check = .FALSE. 
        CLOSE ( iun_source )
        RETURN 
     ELSE 
        xml_check = .TRUE. 
        REWIND( iun_source )
     END IF
  iun_dest = find_free_unit()
  OPEN (UNIT = iun_dest, FILE = TRIM(tempname), STATUS = 'unknown', &
       ACTION = 'write', FORM = 'formatted')
  copy_loop: DO
     ! 
     READ(iun_source, "(a)", IOSTAT = ierr ) line 
     IF (ierr < 0 ) EXIT copy_loop
     !  
     IF ( INDEX(line,"<UPF") /= 0 ) icopy = .TRUE. 
     IF ( .NOT. icopy ) CYCLE copy_loop
     ! 
     WRITE ( iun_dest,"(a)") TRIM(check( line ))
     ! 
     IF ( INDEX( line, "</UPF") /= 0 ) EXIT copy_loop
  END DO copy_loop
  ! 
  CLOSE ( iun_source) 
  CLOSE ( iun_dest )
  !
END FUNCTION make_emended_upf_copy
!

FUNCTION check(in) RESULT (out) 
      CHARACTER (LEN = *)     :: in
#if defined(__PGI)
      INTEGER, PARAMETER      :: length = 255 
      CHARACTER(LEN=length )  :: out 
#else
      CHARACTER(LEN = LEN(in) )  :: out 
#endif 
      INTEGER                :: i, o, disp
      ! 
      disp = 0
      DO i = 1, LEN(in) 
         o = i + disp
         IF ( o > LEN (in) ) EXIT 
         IF (in(i:i) == '&') THEN 
            out(o:o+4) = '&amp;'
            disp = disp+4
         ELSE 
           out(o:o) = in (i:i) 
         END IF
      END DO
END FUNCTION check
END MODULE emend_upf_module 
