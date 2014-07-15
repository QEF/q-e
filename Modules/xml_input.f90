!
! Copyright (C) 2002-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!=----------------------------------------------------------------------------=!
!
MODULE xml_input

   USE xml_io_base, ONLY : attr
   USE qexml_module, ONLY: qexml_init, qexml_write_header, qexml_openfile, &
                           qexml_closefile
   USE iotk_module
   USE kinds

   IMPLICIT NONE
   PRIVATE
   
   PUBLIC :: xml_input_dump

   INTERFACE dump_keyword
      MODULE PROCEDURE dump_keyword_str, dump_keyword_i
   END INTERFACE

   CONTAINS

   SUBROUTINE xml_input_dump
      
      USE io_global,        ONLY : ionode, stdout
      USE io_files,         ONLY : iunpun
      USE global_version,   ONLY : version_number
      USE input_parameters

      CHARACTER(LEN=256) :: filename
      INTEGER            :: ierr

      filename = 'qe_input.xml'
      
      IF ( ionode ) THEN
         !
         ! ... Open XML descriptor
         !
         WRITE( stdout, '(/,3X,"Dumping input parameters",/)' )
         !
         CALL qexml_init( iunpun )
         CALL qexml_openfile( filename, 'write', .FALSE., ierr )
         !
      END IF

      IF ( ionode ) THEN

         CALL iotk_write_attr( attr, "targetNamespace", "http://www.deisa.org/pwscf/3_2", FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "elementFormDefault", "qualified" )
         CALL iotk_write_attr( attr, "xmlns", "http://www.w3.org/2001/XMLSchema" )
         CALL iotk_write_attr( attr, "xmlns:tns", "http://www.deisa.org/pwscf/3_2" )
         CALL iotk_write_begin( iunpun, "schema", attr )

         CALL qexml_write_header( "Quantum ESPRESSO", TRIM(version_number) )

         CALL iotk_write_attr( attr, "section_type", "namelist", FIRST = .TRUE. )
         CALL iotk_write_begin( iunpun, "CONTROLS", attr )
           CALL dump_keyword( "title", title, "namelist", " " )
           CALL dump_keyword( "calculation", calculation, "namelist", " ", calculation_allowed )
           CALL dump_keyword( "verbosity", verbosity, "namelist", " ", verbosity_allowed )
           CALL dump_keyword( "restart_mode", restart_mode, "namelist", " ", restart_mode_allowed )
           CALL dump_keyword( "nstep", nstep, "namelist", " ", min_value = 1 )
           CALL dump_keyword( "iprint", iprint, "namelist", " ", min_value = 1 )
         CALL iotk_write_end( iunpun, "CONTROLS" )

         CALL iotk_write_attr( attr, "section_type", "namelist", FIRST = .TRUE. )
         CALL iotk_write_begin( iunpun, "SYSTEM", attr )
         CALL iotk_write_end( iunpun, "SYSTEM" )

         CALL iotk_write_attr( attr, "section_type", "namelist", FIRST = .TRUE. )
         CALL iotk_write_begin( iunpun, "ELECTRONS", attr )
         CALL iotk_write_end( iunpun, "ELECTRONS" )

         CALL iotk_write_attr( attr, "section_type", "namelist", FIRST = .TRUE. )
         CALL iotk_write_begin( iunpun, "IONS", attr )
         CALL iotk_write_end( iunpun, "IONS" )

         CALL iotk_write_attr( attr, "section_type", "namelist", FIRST = .TRUE. )
         CALL iotk_write_begin( iunpun, "CELL", attr )
         CALL iotk_write_end( iunpun, "CELL" )

         CALL iotk_write_attr( attr, "section_type", "card", FIRST = .TRUE. )
         CALL iotk_write_begin( iunpun, "ATOMIC_SPECIES", attr )
         CALL iotk_write_end( iunpun, "ATOMIC_SPECIES" )

         CALL iotk_write_attr( attr, "section_type", "card", FIRST = .TRUE. )
         CALL iotk_write_begin( iunpun, "ATOMIC_POSITIONS", attr )
         CALL iotk_write_end( iunpun, "ATOMIC_POSITIONS" )

         CALL iotk_write_attr( attr, "section_type", "card", FIRST = .TRUE. )
         CALL iotk_write_begin( iunpun, "K_POINTS", attr )
         CALL iotk_write_end( iunpun, "K_POINTS" )

         CALL iotk_write_end( iunpun, "schema" )

      END IF

      IF ( ionode ) CALL qexml_closefile( 'write', IERR=ierr )

      RETURN
   END SUBROUTINE


   SUBROUTINE dump_keyword_str( kname, defval, usage, descr, allowed )
      USE io_files,         ONLY : iunpun
      CHARACTER(LEN=*) :: kname
      CHARACTER(LEN=*) :: defval
      CHARACTER(LEN=*) :: usage
      CHARACTER(LEN=*) :: descr
      CHARACTER(LEN=*), OPTIONAL :: allowed(:)
         CALL iotk_write_attr( attr, "required", "no", FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "repeat", "no")
         CALL iotk_write_begin( iunpun, "KEYWORD", ATTR = attr )
         CALL iotk_write_attr( attr, "type", "default", FIRST = .TRUE. )
         CALL iotk_write_dat( iunpun, "NAME", kname, ATTR = attr )
         CALL iotk_write_attr( attr, "kind", "STRING", FIRST = .TRUE. )  ! type
         CALL iotk_write_begin( iunpun, "DATA_TYPE", ATTR = attr )
         CALL iotk_write_dat( iunpun, "N_VAR", 1 )
         CALL iotk_write_end( iunpun, "DATA_TYPE" )
         IF( usage == "namelist" ) THEN
            CALL iotk_write_dat( iunpun, "USAGE", kname//" = value" )
         ELSE
            CALL iotk_write_dat( iunpun, "USAGE", usage )
         END IF
         IF( PRESENT( allowed ) ) THEN
            CALL iotk_write_dat( iunpun, "ALLOWED_VALUES", allowed )
         END IF
         CALL iotk_write_dat( iunpun, "DESCRIPTION", descr )
         CALL iotk_write_dat( iunpun, "DEFAULT_VALUE", defval )
         CALL iotk_write_end( iunpun, "KEYWORD" )
      RETURN
   END SUBROUTINE

   SUBROUTINE dump_keyword_i( kname, defval, usage, descr, min_value, max_value )
      USE io_files,         ONLY : iunpun
      CHARACTER(LEN=*) :: kname
      INTEGER          :: defval                                         ! type
      CHARACTER(LEN=*) :: usage
      CHARACTER(LEN=*) :: descr
      INTEGER, OPTIONAL :: min_value                                  ! type
      INTEGER, OPTIONAL :: max_value                                  ! type
         CALL iotk_write_attr( attr, "required", "no", FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "repeat", "no")
         CALL iotk_write_begin( iunpun, "KEYWORD", ATTR = attr )
         CALL iotk_write_attr( attr, "type", "default", FIRST = .TRUE. )
         CALL iotk_write_dat( iunpun, "NAME", kname, ATTR = attr )
         CALL iotk_write_attr( attr, "kind", "INTEGER", FIRST = .TRUE. )  ! type
         CALL iotk_write_begin( iunpun, "DATA_TYPE", ATTR = attr )
         CALL iotk_write_dat( iunpun, "N_VAR", 1 )
         CALL iotk_write_end( iunpun, "DATA_TYPE" )
         IF( usage == "namelist" ) THEN
            CALL iotk_write_dat( iunpun, "USAGE", kname//" = value" )
         ELSE
            CALL iotk_write_dat( iunpun, "USAGE", usage )
         END IF
         IF( PRESENT( min_value ) ) THEN
            CALL iotk_write_dat( iunpun, "MIN_VALUE", min_value )
         END IF
         IF( PRESENT( max_value ) ) THEN
            CALL iotk_write_dat( iunpun, "MAX_VALUE", max_value )
         END IF
         CALL iotk_write_dat( iunpun, "DESCRIPTION", descr )
         CALL iotk_write_dat( iunpun, "DEFAULT_VALUE", defval )
         CALL iotk_write_end( iunpun, "KEYWORD" )
      RETURN
   END SUBROUTINE


END MODULE
