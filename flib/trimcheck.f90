!
!-----------------------------------------------------------------------
FUNCTION trimcheck ( outdir )
  !-----------------------------------------------------------------------
  !
  ! ... verify if outdir ends with /, add one if needed; 
  ! ... trim white spaces and put the result in trimcheck
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: outdir
  CHARACTER (LEN=256) :: trimcheck
  INTEGER  :: l
  !
  l = LEN_TRIM( outdir )
  IF ( l == 0 ) CALL errore( 'trimcheck', ' input name empty', 1)
  !
  IF ( outdir(l:l) == '/' ) THEN
     trimcheck = TRIM ( outdir)
  ELSE
     IF ( l < LEN( trimcheck ) ) THEN
        trimcheck = TRIM ( outdir ) // '/'
     ELSE
        CALL errore(  'trimcheck', ' input name too long', l )
     END IF
  END IF
  !
  RETURN
  !
END FUNCTION trimcheck

