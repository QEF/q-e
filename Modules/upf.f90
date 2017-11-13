! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
      MODULE upf_module
!=----------------------------------------------------------------------------=!
!  this module handles reading of unified pseudopotential format (UPF)
!  in either v1 or v2 format
!
! A macro to trim both from left and right
#define TRIM(a) trim(adjustl(a))
      !
      USE kinds,        ONLY: DP
      USE pseudo_types, ONLY: pseudo_upf, deallocate_pseudo_upf
      !
      IMPLICIT NONE
      PUBLIC
      !PRIVATE
      !PUBLIC :: read_upf, pseudo_upf, deallocate_pseudo_upf
      !
      CONTAINS

!------------------------------------------------+
SUBROUTINE read_upf(upf, grid, ierr, unit,  filename, xml_only) !
   !---------------------------------------------+
   ! Read pseudopotential in UPF format (either v.1 or v.2 or upf_schema)
   ! ierr = -2 : read upf_schema
   ! ierr = -1 : read UPF v.1 
   ! ierr =  0 : read UPF v.2 
   ! ierr =  1 : not an UPF file, or error while reading
   !
   USE radial_grids, ONLY: radial_grid_type, deallocate_radial_grid
   USE read_upf_v1_module,ONLY: read_upf_v1
   USE read_upf_v2_module,ONLY: read_upf_v2
   USE read_upf_schema_module,ONLY: read_upf_schema
   USE mp,           ONLY: mp_barrier
   USE mp_world,     ONLY: world_comm
   USE io_global,    ONLY: ionode
   USE io_files,     ONLY: tmp_dir
   USE FoX_DOM,      ONLY: Node, domException, parseFile, getFirstChild, getExceptionCode,&
                              getTagName    
   USE wrappers,     ONLY: f_remove
   IMPLICIT NONE
   INTEGER,INTENT(IN), OPTIONAL            :: unit
   !! i/o unit: used only to read upf version 1 files !!!!        
   CHARACTER(len=*),INTENT(IN),OPTIONAL    :: filename  
   !! i/o filename
   LOGICAL,INTENT(IN), OPTIONAL            :: xml_only
   !! if present and true the program will parse only xml documents neglecting version 1 upf format
   TYPE(pseudo_upf),INTENT(INOUT) :: upf       
!! the pseudo data
   TYPE(radial_grid_type),OPTIONAL,INTENT(INOUT),TARGET :: grid
   INTEGER,INTENT(OUT) :: ierr
   !
   LOGICAL            :: xml_only_ = .FALSE. 
   TYPE(Node),POINTER :: u,doc     
   INTEGER            :: u_temp,&    ! i/o unit in case of upf v1
                         iun, ferr  
   TYPE(DOMException) :: ex 
   INTEGER, EXTERNAL  :: find_free_unit

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
          IF ( ionode ) CALL make_emended_upf_copy( TRIM(filename), TRIM(tmp_dir)//'tmp.UPF')     
          CALL mp_barrier ( world_comm) 
          doc => parseFile(TRIM(tmp_dir)//'tmp.UPF', EX = ex )
          ierr = getExceptionCode( ex ) 
          CALL mp_barrier(world_comm) 
          IF (ionode) ferr = f_remove(TRIM(tmp_dir)//'tmp.UPF')
       END IF 
       IF ( ierr == 0 ) THEN 
           u => getFirstChild(doc) 
           SELECT CASE (TRIM(getTagname(u))) 
              CASE ('UPF') 
                 CALL read_upf_v2( u, upf, grid, ierr )
              CASE ('pseudo') 
                 CALL read_upf_schema( u, upf, grid, ierr)
                 IF ( ierr == 0 ) ierr = -2
              CASE default 
                 ierr = 1
                 CALL errore('read_upf', 'unrecognized xml format '//TRIM(getTagName(u)),ierr) 
           END SELECT 
           IF ( ierr > 0 ) CALL errore( 'read_upf', 'File is Incomplete or wrong: '//TRIM(filename), ierr)
           !
           RETURN
           !  
       ELSE IF ( ierr > 0 ) THEN
          ! 
          IF ( .NOT. xml_only ) THEN 
             u_temp = find_free_unit()
             open (unit = u_temp, file = TRIM(filename), status = 'old', form = 'formatted', iostat = ierr)
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
  CONTAINS 
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
END SUBROUTINE  make_emended_upf_copy 
!=----------------------------------------------------------------------------=!
      END MODULE upf_module
!=----------------------------------------------------------------------------=!
#undef TRIM

