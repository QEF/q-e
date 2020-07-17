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
      USE upf_kinds,           ONLY: DP
      !
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: read_ps, read_upf_new
      !
    CONTAINS
      !
SUBROUTINE read_ps ( filein, upf_in )
  !
  ! stripped-down version of readpp in Modules/read_pseudo.f90:
  ! for serial execution only
  !
  !USE read_upf_new_module,ONLY: read_upf_new
  USE read_upf_v1_module, ONLY: read_upf_v1
  USE emend_upf_module, ONLY: make_emended_upf_copy
  USE pseudo_types,     ONLY: pseudo_upf
  USE upf_auxtools,     ONLY: upf_get_pp_format 
  USE upf_to_internal,  ONLY: set_upf_q
  USE read_uspp_module, ONLY: readvan, readrrkj
  USE cpmd_module,      ONLY: read_cpmd
  USE m_gth,            ONLY: readgth
  USE fhi,              ONLY: readfhi
  !
  IMPLICIT NONE
  CHARACTER(LEN=256), INTENT(in) :: filein
  TYPE(pseudo_upf), INTENT(out) :: upf_in
  !
  LOGICAL :: is_xml
  INTEGER :: isupf, iunps = 4, stdout = 6
  !
  is_xml = .false.
  isupf = 0
  CALL read_upf_new( filein, upf_in, isupf )
  IF (isupf ==-81 ) THEN
     is_xml = make_emended_upf_copy( filein, 'tmp.upf' )
     IF (is_xml) THEN
        CALL  read_upf_new( 'tmp.upf', upf_in, isupf )
        !! correction succeeded, try to read corrected file
        OPEN ( unit=iunps, iostat=isupf, file='tmp.upf', status='old')
        CLOSE( unit=iunps, status='delete' )
     ELSE
        CALL  read_upf_v1 ( filein, upf_in, isupf )
        !! try to read UPF v.1 file
        IF ( isupf == 0 ) isupf = -1
     END IF
     !
  END IF
  !
  IF (isupf == -2 .OR. isupf == -1 .OR. isupf == 0) THEN
     !
     IF ( isupf == 0 ) THEN
        WRITE( stdout, "('file type is xml')") 
     ELSE
        IF ( is_xml ) THEN
           WRITE( stdout, "('file type is UPF v.',I1,' with invalid characters &')") ABS(isupf)
        ELSE
           WRITE( stdout, "('file type is UPF v.',I1)") ABS(isupf)
        END IF
     END IF
     !
  ELSE
     !
     OPEN ( UNIT=iunps, FILE=filein, STATUS = 'old', FORM = 'formatted' ) 
     !
     !     The type of the pseudopotential is determined by the file name:
     !    *.xml or *.XML  UPF format with schema              pp_format=0
     !    *.upf or *.UPF  UPF format                          pp_format=1
     !    *.vdb or *.van  Vanderbilt US pseudopotential code  pp_format=2
     !    *.gth           Goedecker-Teter-Hutter NC pseudo    pp_format=3
     !    *.RRKJ3         PWSCF US PP format ("atomic" code)  pp_format=4
     !    *.fhi or *cpi   FHI/Abinit numerical NC PP          pp_format=6
     !    *.cpmd          CPMD NC pseudopot. file format      pp_format=7
     !    none of the above: PWSCF norm-conserving format     pp_format=5
     !
     IF ( upf_get_pp_format( filein ) == 2  ) THEN
        !
        WRITE( stdout, "('file type is Vanderbilt US PP')")
        CALL readvan (iunps, 1, upf_in)
        !
     ELSE IF ( upf_get_pp_format( filein ) == 3 ) THEN
        !
        WRITE( stdout, "('file type is GTH (analytical)')")
        CALL readgth (iunps, 1, upf_in)
        !
     ELSE IF ( upf_get_pp_format( filein ) == 4 ) THEN
        !
        WRITE( stdout, "('file type is RRKJ3')")
        CALL readrrkj (iunps, 1, upf_in)
        !
     ELSE IF ( upf_get_pp_format( filein ) == 5 ) THEN
        !
        WRITE( stdout, "('file type is old PWscf NC format')")
        CALL read_ncpp (iunps, 1, upf_in)
        !
     ELSE IF ( upf_get_pp_format( filein ) == 6 ) THEN
        !
        WRITE( stdout, "('file type is FHI .cpi or .fhi format')")
        CALL readfhi (iunps, upf_in)
        !
     ELSE IF ( upf_get_pp_format( filein ) == 7 ) THEN
        !
        WRITE( stdout, "('file type is CPMD NC format')")
        CALL read_cpmd (iunps, upf_in)
        !
     ELSE
        !
        CALL upf_error('readpp', 'file '//TRIM(filein)//' not readable',1)
        !
     ENDIF
     !
     ! end of reading
     !
     CLOSE (iunps)
     !
  ENDIF
  !
  ! reconstruct Q(r) if needed
  !
  CALL set_upf_q (upf_in)
  !
  RETURN
END SUBROUTINE read_ps

!------------------------------------------------+
SUBROUTINE read_upf_new(filename, upf, ierr)         !
   !---------------------------------------------+
   !! Reads pseudopotential in UPF format (either v.2 or upf_schema).
   !! Derived-type variable *upf* store in output the data read from file. 
   !! File *filename* is opened and closed inside the routine
   !! @Note last revision: 29-04-2020 PG - UPF v.1 no longer read here
   !! @Note last revision: 27-04-2020 PG - "grid" variable moved out
   !! @Note last revision: 01-01-2019 PG - upf fix moved out from here
   !! @Note last revision: 11-05-2018 PG - removed xml_only
   !
   USE read_upf_v2_module, ONLY: read_upf_v2
   USE read_upf_schema_module, ONLY: read_upf_schema
   USE pseudo_types,        ONLY: pseudo_upf, deallocate_pseudo_upf
   USE FoX_dom, ONLY: Node, domException, parseFile, getFirstChild, &
        getExceptionCode, getTagName    
   IMPLICIT NONE
   CHARACTER(len=*),INTENT(IN),OPTIONAL    :: filename  
   !! i/o filename
   TYPE(pseudo_upf),INTENT(INOUT) :: upf       
   !! the derived type storing the pseudo data
   INTEGER,INTENT(INOUT) :: ierr
   !! On input:
   !! ierr =0:   return if not a valid xml schema or UPF v.2 file
   !! ierr/=0: continue if not a valid xml schema or UPF v.2 file
   !! On output:
   !! ierr=0: xml schema, ierr=-2: UPF v.2
   !! ierr>0: error reading PP file
   !! ierr=-81: error reading PP file, possibly UPF fix needed
   !
   TYPE(Node),POINTER :: u,doc     
   INTEGER            :: iun, ferr  
   TYPE(DOMException) :: ex 

   ferr = ierr
   ierr = 0
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
         CALL read_upf_v2( u, upf, ierr )
         IF ( ierr == 0 ) ierr = -2
      CASE ('qe_pp:pseudo') 
         CALL read_upf_schema( u, upf, ierr)
      CASE default 
         ierr = 1
         CALL upf_error('read_upf_new', 'format '//TRIM(getTagName(u))//' not implemented', ierr) 
      END SELECT
   END IF
   IF ( ierr > 0 ) CALL upf_error( 'read_upf_new', 'File is incomplete or wrong: '//TRIM(filename), ierr)
   !
 END SUBROUTINE read_upf_new
 
!=----------------------------------------------------------------------------=!
      END MODULE upf_module
!=----------------------------------------------------------------------------=!

