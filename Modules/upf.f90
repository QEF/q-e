! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
      MODULE upf_module
!=----------------------------------------------------------------------------=!
!! author: Unknown 
!! this module handles reading of unified pseudopotential format (UPF)
!! in either v1 or v2 or schema format.
!! @Note
!! 14/11/17 Pietro Delugas: new revision passed from iotk to FoX lib,
!!                          added support for new schema format
      !
      USE kinds,               ONLY: DP
      USE pseudo_types,        ONLY: pseudo_upf, deallocate_pseudo_upf
      USE read_upf_v1_module,  ONLY: scan_begin, scan_end
      !
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: read_upf, scan_begin, scan_end
      !
      CONTAINS

!------------------------------------------------+
SUBROUTINE read_upf(upf, grid, ierr, unit,  filename, xml_only) !
   !---------------------------------------------+
   !! Reads pseudopotential in UPF format (either v.1 or v.2 or upf_schema).
   !! Derived-type variable *upf* and optionally *grid* store in output the 
   !! data read from file. 
   !! If unit number is provided with the *unit* argument, only UPF v1 format
   !! is chhecked; the PP file must be opened and closed outside the routine.  
   !! Otherwise the *filename* argument must be given, file is opened and closed
   !! inside the routine, all formats will be  checked. The optional *xml_only* 
   !! argument may be set to true value to prevent the routine from checking
   !! v1 format. 
   !! @Note last revision: 14-11-2017
   !
   USE radial_grids, ONLY: radial_grid_type, deallocate_radial_grid
   USE read_upf_v1_module,ONLY: read_upf_v1
   USE read_upf_v2_module,ONLY: read_upf_v2
   USE read_upf_schema_module ,ONLY: read_upf_schema
   USE mp,           ONLY: mp_bcast, mp_sum
   USE mp_images,    ONLY: intra_image_comm, my_image_id
   USE io_global,    ONLY: ionode, ionode_id, stdout
   USE io_files,     ONLY: tmp_dir
   USE FoX_DOM,      ONLY: Node, domException, parseFile, getFirstChild, getExceptionCode,&
                              getTagName    
   USE wrappers,     ONLY: f_remove
   USE emend_upf_module, ONLY: make_emended_upf_copy
   IMPLICIT NONE
   INTEGER,INTENT(IN), OPTIONAL            :: unit
   !! i/o unit:    
   CHARACTER(len=*),INTENT(IN),OPTIONAL    :: filename  
   !! i/o filename
   LOGICAL,INTENT(IN), OPTIONAL            :: xml_only
   !! if present and true the program will parse only xml documents neglecting version 1 upf format
   TYPE(pseudo_upf),INTENT(INOUT) :: upf       
   !! the derived type storing the pseudo data
   TYPE(radial_grid_type),OPTIONAL,INTENT(INOUT),TARGET :: grid
   !! derived type where is possible to store data on the radial mesh
   INTEGER,INTENT(OUT) :: ierr
   !! ierr=0: xml schema, ierr=-1: UPF v.1,  ierr=-2: UPF v.2
   !! ierr>0: error reading PP file
   !
   LOGICAL            :: xml_only_ = .FALSE. 
   TYPE(Node),POINTER :: u,doc     
   INTEGER            :: u_temp,&    ! i/o unit in case of upf v1
                         iun, ferr  
   TYPE(DOMException) :: ex 
   INTEGER, EXTERNAL  :: find_free_unit
   CHARACTER(LEN=256) :: temp_upf_file
   CHARACTER(LEN=1024) :: msg
   LOGICAL             :: should_be_xml
   IF (PRESENT(xml_only) ) xml_only_ = xml_only
   ierr = 0

   IF ( present ( unit ) ) THEN 
      REWIND (unit) 
      CALL deallocate_pseudo_upf(upf) 
      CALL deallocate_radial_grid( grid ) 
      CALL read_upf_v1 (unit, upf, grid, ierr ) 
      IF (ierr == 0 ) ierr = -1     
      !
      RETURN
      ! 
   ELSE IF (PRESENT(filename) ) THEN
      doc => parseFile(TRIM(filename), EX = ex )
      ierr = getExceptionCode( ex )
      IF ( ierr ==  81 ) THEN 
         WRITE(temp_upf_file, '("tmp_",I0,".UPF")') my_image_id  
         IF ( ionode ) THEN
            CALL make_emended_upf_copy( TRIM(filename), TRIM(tmp_dir)//trim(temp_upf_file), should_be_xml)  
         END IF   
         CALL mp_bcast ( should_be_xml, ionode_id, intra_image_comm)     
         IF ( should_be_xml) THEN 
            doc => parseFile(TRIM(tmp_dir)//trim(temp_upf_file), EX = ex, IOSTAT = ferr )
            ierr = getExceptionCode( ex ) 
            CALL mp_sum(ferr,intra_image_comm) 
            IF ( ferr /= 0 ) THEN 
               WRITE (msg, '(A)')  'Failure while trying to fix '//trim(filename) // '.'// new_line('a') // &
                                    'For fixing manually UPF files see: '// new_line('a') // &
                                    'https://gitlab.com/QEF/q-e/blob/master/upftools/how_to_fix_upf.md'
               CALL errore('read_upf: ', TRIM(msg), ferr ) 
            ELSE 
               WRITE ( msg, '(A)') 'Pseudo file '// trim(filename) // ' has been successfully fixed on the fly.' &
                              // new_line('a') // 'To avoid this message in the future you can permanently fix ' &
                              // new_line('a') // ' your pseudo files following instructions given in: ' &
                              // new_line('a') // 'https://gitlab.com/QEF/q-e/blob/master/upftools/how_to_fix_upf.md'
               CALL infomsg('read_upf:', trim(msg) )    
            END IF
         END IF
            ! 
            IF (ionode) ferr = f_remove(TRIM(tmp_dir)//TRIM(temp_upf_file) )
            temp_upf_file=""
         END IF 
         IF ( ierr == 0 ) THEN 
            u => getFirstChild(doc) 
            SELECT CASE (TRIM(getTagname(u))) 
               CASE ('UPF') 
                  CALL read_upf_v2( u, upf, grid, ierr )
                  IF ( ierr == 0 ) ierr = -2
               CASE ('qe_pp:pseudo') 
                  CALL read_upf_schema( u, upf, grid, ierr)
               CASE default 
                  ierr = 1
                  CALL errore('read_upf', 'xml format '//TRIM(getTagName(u))//' not implemented', ierr) 
               END SELECT 
            IF ( ierr > 0 ) CALL errore( 'read_upf', 'File is Incomplete or wrong: '//TRIM(filename), ierr)
            RETURN
            !  
         ELSE IF ( ierr > 0 ) THEN
         ! 
         IF ( .NOT. xml_only_ ) THEN
            u_temp = find_free_unit()
            OPEN (UNIT = u_temp, FILE = TRIM(filename), STATUS = 'old', FORM = 'formatted', IOSTAT = ierr)
            CALL errore ("upf_module:read_upf", "error while opening file " // TRIM(filename), ierr) 
            CALL deallocate_pseudo_upf( upf )
            CALL deallocate_radial_grid( grid )
            CALL read_upf_v1( u_temp, upf, grid, ierr )
            IF ( ierr == 0 ) ierr = -1
            CLOSE ( u_temp)  
         END IF
          !
         RETURN
          !
       END IF
   ELSE 
       CALL errore('read_upf',&
         'Nothing to read !!! you must provide one of filename or unit optional arguments',1)
   END IF
END SUBROUTINE read_upf
!=----------------------------------------------------------------------------=!
      END MODULE upf_module
!=----------------------------------------------------------------------------=!

