!
! Copyright (C) 2018 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
PROGRAM upf_fixer
   !---------------------------------------------------------------------
   !!  fix upf files which contain prohibited escape characters 
   !!  usage: fix_upf.x [--inplace] <filename> 
   !!  if the --inplace option is specified is rewritten otherwise the 
   !!  emended file written as <prefix>_sano.UPF, where prefix is the part 
   !!  of the filename preceding the .UPF extension
   !
   USE wrappers, ONLY: f_copy, f_remove
   USE emend_upf_module, ONLY: make_emended_upf_copy, check_upf_file
   IMPLICIT NONE
   CHARACTER(LEN=256)   :: filein, fileout, line 
   INTEGER              :: ios, argc, prefix_len, iarg
   LOGICAL              :: exst, sano, is_xml, inplace = .FALSE.
   argc = command_argument_count() 
   IF ( .NOT. argc > 0 ) THEN 
      PRINT *, 'usage: fix_upf [--inplace ] filename' 
      STOP
   END IF 
   filein=''
   DO iarg = 1, argc 
      CALL get_command_argument(iarg, line)
      IF ( TRIM(line)  == '--inplace') THEN 
         inplace = .TRUE. 
      ELSE IF ( INDEX(line,'--') == 1 ) THEN 
         PRINT '("unrecognized option ",A)', TRIM(line) 
      ELSE IF ( TRIM (filein) == '') THEN 
         filein=TRIM(line)
      ELSE
         PRINT '("unrecognized argument ",A)', TRIM(line) 
      ENDIF 
   ENDDO
   INQUIRE(FILE = TRIM(filein), EXIST = exst)  
   IF ( exst) THEN
      OPEN (unit = 2 , file=TRIM(filein) ,ACTION = 'read' ) 
      READ (2,*) line 
      CLOSE(2)
      IF (INDEX(line, '<?xml') == 0 .AND. INDEX(line,'<UPF') == 0) THEN
         PRINT '(A," is not an xml file, stopping")', TRIM(filein)
         STOP
      ENDIF
   ELSE
      PRINT '("File ",A," not found")', TRIM(filein) 
      STOP
   ENDIF
   PRINT '("Checking file ", A)', TRIM(filein)
   IF ( inplace) THEN
      fileout = TRIM(filein)
      ios = f_copy(TRIM(filein),TRIM(filein)//'_bak')  
      IF ( ios /= 0 ) THEN 
         PRINT '("error while making backup copy of ",A)' , TRIM(filein) 
         STOP
      ENDIF
   ELSE
      IF ( INDEX(TRIM(filein),'.UPF' ) > 0) THEN 
         prefix_len = INDEX(TRIM(filein),'.UPF' )
      ELSE IF (INDEX(TRIM(filein),'.upf') > 0 ) THEN
         prefix_len = INDEX(TRIM(filein),'.upf') 
      ELSE 
         prefix_len = LEN(TRIM(filein)) 
      ENDIF
      fileout = filein(1:prefix_len) //'_sano.UPF'
   ENDIF
   sano =  check_upf_file(TRIM(filein)) 
   IF (.NOT. sano ) THEN 
      CALL make_emended_upf_copy( filein, './temp.UPF', is_xml)
      IF ( .NOT. is_xml ) THEN 
         PRINT *, "This file is not in xml format !!!! Stopping"
         STOP
      END IF
      ios = f_copy('./temp.UPF', TRIM(fileout)) 
      IF ( ios /= 0 ) THEN 
         PRINT '("error while writing the emended copy of ",A)', TRIM(filein) 
         STOP
      ENDIF
      ios = f_remove('./temp.UPF') 
      PRINT '("Error detected in ",A,A,"Writing emended copy in ",A)', &
                                                         TRIM(filein), '.'//new_line('a'),TRIM(fileout)
   ELSE
      PRINT '("No error detected in file ",A)' , TRIM(filein) 
   ENDIF  
END PROGRAM
