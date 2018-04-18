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
PUBLIC make_emended_upf_copy, check_upf_file 

CONTAINS 
SUBROUTINE make_emended_upf_copy( filename, tempname, xml_check) 
  !! author: Pietro Delugas
  !! Utility to make the old UPF format readable by FoX
  !! Replaces "&" with "&amp;" in file "filename", writes to file "tempname"
  !
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)      :: filename, tempname
  LOGICAL,INTENT(OUT)              :: xml_check
  !
  INTEGER                          :: iun_source, iun_dest, ierr 
  INTEGER,EXTERNAL                 :: find_free_unit 
  LOGICAL                          :: icopy = .FALSE.
  CHARACTER(LEN=256)               :: line 
  ! 
  iun_source = find_free_unit()
  OPEN (UNIT = iun_source, FILE = TRIM(filename), STATUS = 'old', &
       ACTION = 'read', FORM='formatted', iostat=ierr)
  IF ( ierr /= 0 ) CALL errore ("make_emended_upf", &
          "error opening file " // TRIM (filename),abs(ierr))
  READ ( iun_source, "(a256)", IOSTAT = ierr ) line
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
     READ(iun_source, "(a256)", IOSTAT = ierr ) line 
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
END SUBROUTINE  make_emended_upf_copy 
!
FUNCTION check_upf_file(filename, errcode) RESULT(ok)
   !! checks whether the upf file filename is complian to xml syntax
   !! the errorcode returned by the checking routine may optionally be 
   !! written in the errorcode argument
   USE FoX_dom,   ONLY: Node, DOMException, parseFile, getExceptionCode                  
   IMPLICIT NONE
   CHARACTER(LEN=*),INTENT(IN)   :: filename
   !! name of the upf file being checked 
   INTEGER,OPTIONAL,INTENT(OUT)  :: errcode
   !! if present contains the error code returnd by the upf check
   LOGICAL                       :: ok 
   !!  if true the upf file is compliant to xml syntax
   !
   TYPE(Node),POINTER    :: doc 
   TYPE(DOMException)    :: dom_ex
   INTEGER               :: ierr 
   doc => parseFile(TRIM(filename), EX = dom_ex) 
   ierr = getExceptionCode(dom_ex) 
   IF (PRESENT(errcode)) errcode=ierr 
   IF (ierr /= 0 ) THEN
      ok = .FALSE.
   ELSE
      ok =.TRUE.
   ENDIF 
   !
END FUNCTION check_upf_file


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
