!
! Copyright (C) 2008-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!--------------------------------------------------------
MODULE  write_upf_module
   !-----------------------------------------------------
   !! this module collects the interfaces for the subroutines 
   !! writing pseudopotential usin XML. 
   !! write_upf_v2 (xf, upf, conf, u_input) writes using UPF v2.1 format 
   !! write_upf_schema(xf, upf, conf, u_input) writes according to qe_pseudo schema
   !! *xf* is a xmlf_f derived type from FoX library: is the xml file descriptor, it must be initialized with
   !! xml_OpenFile before calling the routine, the file is closed before inside write_upf routines. 
   !! *upf* is  pseudo_upf derived type from qe Module pseudo_types: stores the whole pseudo data set
   !! *conf* (optional) is a pseudo_conf derived type from qe Module pseudo_types: stores the  info about 
   !!        the configuration used to generate the pseudo. 
   !! *u_input* is the Fortran unit number pointing to the input file of the pseudo creation program
   
   USE write_upf_v2_module, ONLY: write_upf_v2
   USE write_upf_schema_module, ONLY: write_upf_schema
   PUBLIC 
CONTAINS 
   SUBROUTINE write_upf(filename, upf, schema , conf, u_input) 
      !! this routine writes an upf file with name = filename using 
      !! data contained in pseudo_upf structure upf 
      !! if schema = 'qe_pp' the new qe_pp schema the old upf v2.1 is used 
      !! otherwise
      USE FoX_wxml 
      USE pseudo_types, ONLY: pseudo_upf, pseudo_config
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN)   :: filename
      !! name of the output file
      TYPE(pseudo_upf),INTENT(IN)   :: upf
      !! pseudo_upf structure containing all the pseudo data
      CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: schema
      !!  optional, character flag which selects what schema will be used on writing 
      TYPE(pseudo_config),OPTIONAL,INTENT(IN)   :: conf 
      !!  optional, pseudo_conf data structure containing the atomic configuration used
      !!  to generate the pseudo
      INTEGER,OPTIONAL                        :: u_input
      !!  optional: unit of stdin for  the generation program, used to write the 
      !!  generation input in the upf file
      TYPE(xmlf_t)                   :: xf_desc
      CHARACTER(LEN=10)              :: schema_
      !
      CALL  xml_OpenFile (filename = TRIM(filename), XF = xf_desc, UNIT = 2, PRETTY_PRINT =.true., &
                                 REPLACE = .true., NAMESPACE = .true. )
      IF ( PRESENT(schema )) schema_ = TRIM(schema) 
      SELECT CASE (TRIM(schema_))
         CASE ('qe_pp', 'QE_PP') 
            CALL write_upf_schema(xf_desc, upf = upf, CONF = conf, U_INPUT = u_input) 
         CASE ('V2', 'v2' ,'upf', 'UPF') 
            CALL write_upf_v2(xf_desc,upf=upf, CONF = conf, U_INPUT = u_input)
      END SELECT
   END SUBROUTINE write_upf

END MODULE write_upf_module

