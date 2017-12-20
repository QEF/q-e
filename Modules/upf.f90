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
!! in either v1 or v2 of schema format.
!! @Note
!! 14/11/17 Pietro Delugas: new revision passed from iotk to FoX lib, added support 
!!  new schema for UPF files 
!
!
! A macro to trim both from left and right
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
   !! Reads pseudopotential in UPF format (either v.1 or v.2 or upf_schema).
   !! Derived type variable *upf* and optionally *grid* store in output the data read
   !! from file. 
   !! If the unit number is provided with the *unit* argument  only UPF v1 format 
   !! is checked, the pseudo file must be opened and closed outside the routine.  
   !! Otherwise the *filename* argument must be given, file is opened and closed inside 
   !! the routine   and all formats will be  checked. The logical *xml_only* optional 
   !! argument may be given with true value to prevent the routine to check v1 format. 
   !! @Note last revision: 14-11-2017
   !
   USE radial_grids, ONLY: radial_grid_type, deallocate_radial_grid
   USE read_upf_v1_module,ONLY: read_upf_v1
   USE read_upf_v2_module,ONLY: read_upf_v2
   USE read_upf_schema_module ,ONLY: read_upf_schema
   USE mp,           ONLY: mp_barrier
   USE mp_world,     ONLY: world_comm
   USE io_global,    ONLY: ionode
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
                 !CALL errore('read_upf', 'pseudo format reader under maintence' ,ierr)
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
!=----------------------------------------------------------------------------=!
      END MODULE upf_module
!=----------------------------------------------------------------------------=!

