MODULE emend_upf_module 

PRIVATE 
PUBLIC make_emended_upf_copy 

CONTAINS 
SUBROUTINE make_emended_upf_copy( filename, tempname) 
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN)      :: filename, tempname
  !
  INTEGER                          :: iun_source, iun_dest, ierr 
  INTEGER,EXTERNAL                 :: find_free_unit 
  LOGICAL                          :: icopy = .FALSE.
  CHARACTER(LEN=256)               :: line 
  ! 
  iun_source = find_free_unit()
  OPEN (UNIT = iun_source, FILE = TRIM(filename), STATUS = 'old', ACTION = 'read', FORM='formatted')
  iun_dest = find_free_unit()
  OPEN (UNIT = iun_dest, FILE = TRIM(tempname), STATUS = 'unknown', ACTION = 'write', FORM = 'formatted')
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
