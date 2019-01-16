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
      PUBLIC :: read_upf, check_upf_file, scan_begin, scan_end
      !
      CONTAINS

!------------------------------------------------+
SUBROUTINE read_upf(upf, grid, ierr, unit,  filename) !
   !---------------------------------------------+
   !! Reads pseudopotential in UPF format (either v.1 or v.2 or upf_schema).
   !! Derived-type variable *upf* and optionally *grid* store in output the 
   !! data read from file. 
   !! If unit number is provided with the *unit* argument, only UPF v1 format
   !! is checked; the PP file must be opened and closed outside the routine.  
   !! Otherwise the *filename* argument must be given, file is opened and closed
   !! inside the routine, all formats will be  checked. 
   !! @Note last revision: 01-01-2019 PG - upf fix moved out from here
   !! @Note last revision: 11-05-2018 PG - removed xml_only
   !
   USE radial_grids, ONLY: radial_grid_type, deallocate_radial_grid
   USE read_upf_v1_module,ONLY: read_upf_v1
   USE read_upf_v2_module,ONLY: read_upf_v2
   USE read_upf_schema_module ,ONLY: read_upf_schema
   USE FoX_DOM,      ONLY: Node, domException, parseFile, getFirstChild, &
        getExceptionCode, getTagName    
   IMPLICIT NONE
   INTEGER,INTENT(IN), OPTIONAL            :: unit
   !! i/o unit:    
   CHARACTER(len=*),INTENT(IN),OPTIONAL    :: filename  
   !! i/o filename
   TYPE(pseudo_upf),INTENT(INOUT) :: upf       
   !! the derived type storing the pseudo data
   TYPE(radial_grid_type),OPTIONAL,INTENT(INOUT),TARGET :: grid
   !! derived type where is possible to store data on the radial mesh
   INTEGER,INTENT(INOUT) :: ierr
   !! On input:
   !! ierr =0:   return if not a valid xml schema or UPF v.2 file
   !! ierr/=0: continue if not a valid xml schema or UPF v.2 file
   !! On output:
   !! ierr=0: xml schema, ierr=-1: UPF v.1,  ierr=-2: UPF v.2
   !! ierr>0: error reading PP file
   !! ierr=-81: error reading PP file, possibly UPF fix needed
   !
   TYPE(Node),POINTER :: u,doc     
   INTEGER            :: u_temp,&    ! i/o unit in case of upf v1
                         iun, ferr  
   TYPE(DOMException) :: ex 
   INTEGER, EXTERNAL  :: find_free_unit

   ferr = ierr
   ierr = 0
   IF ( present ( unit ) ) THEN 
      REWIND (unit) 
      CALL deallocate_pseudo_upf(upf) 
      CALL deallocate_radial_grid( grid ) 
      CALL read_upf_v1 (unit, upf, grid, ierr ) 
      IF (ierr == 0 ) ierr = -1     
      !
   ELSE IF (PRESENT(filename) ) THEN
      doc => parseFile(TRIM(filename), EX = ex )
      ierr = getExceptionCode( ex )
      IF ( ferr == 0 .AND. ierr ==  81 ) THEN
         ierr = -81
         RETURN
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
         !  
      ELSE IF ( ierr > 0 ) THEN
         ! 
         u_temp = find_free_unit()
         OPEN (UNIT = u_temp, FILE = TRIM(filename), STATUS = 'old', FORM = 'formatted', IOSTAT = ierr)
         CALL errore ("upf_module:read_upf", "error while opening file " // TRIM(filename), ierr) 
         CALL deallocate_pseudo_upf( upf )
         CALL deallocate_radial_grid( grid )
         CALL read_upf_v1( u_temp, upf, grid, ierr )
         IF ( ierr == 0 ) ierr = -1
         CLOSE ( u_temp)  
         !
      END IF
   ELSE 
      CALL errore('read_upf', 'Nothing to read !!! Provide either filename or unit optional arguments',1)
   END IF
   !
 END SUBROUTINE read_upf

FUNCTION check_upf_file(filename, errcode) RESULT(ok)
   !! checks whether the upf file filename is compliant with xml syntax
   !! the error code returned by the checking routine may optionally be 
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

!=----------------------------------------------------------------------------=!
      END MODULE upf_module
!=----------------------------------------------------------------------------=!

