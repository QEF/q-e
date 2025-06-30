! Copyright (C) 2003-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE qexsd_module
  !----------------------------------------------------------------------------
  !! This module contains subroutines used to read and write in XML format,
  !! according to the "schema", the data produced by Quantum ESPRESSO.
  !
  !! Based on initial work by Carlo Sbraccia (2003)
  !! and on the qexml.f90 routines written by Andrea Ferretti (2006)
  !! Modified by Simone Ziraldo (2013).
  !! Rewritten by Giovanni Borghi, A. Ferretti, et al. (2015).
  !! Heavily modified by Pietro Delugas and Paolo Giannozzi (2016 on)
  !
  !
  USE kinds,            ONLY : DP
  USE input_parameters, ONLY : input_xml_schema_file
  USE mp_world,         ONLY : nproc
  USE mp_images,        ONLY : nimage,nproc_image
  USE mp_pools,         ONLY : npool
  USE mp_bands,         ONLY : ntask_groups, nproc_bgrp, nbgrp
  USE global_version,   ONLY : version_number
  !
  USE qes_types_module
  USE qes_write_module, ONLY : qes_write
  USE qes_reset_module, ONLY : qes_reset
  USE qes_init_module,  ONLY : qes_init
  !
#if defined (__fox)
  USE FoX_wxml,  ONLY : xmlf_t, xml_OpenFile, xml_DeclareNamespace, &
                        xml_NewElement, xml_addAttribute, xml_addComment,&
                        xml_AddCharacters, xml_EndElement, xml_Close
  USE FoX_dom,   ONLY : parseFile, item, getElementsByTagname, &
                        destroy, nodeList, Node
#else
  USE     wxml,  ONLY : xmlf_t, xml_OpenFile, xml_DeclareNamespace, &
                        xml_NewElement, xml_addAttribute, xml_addComment,&
                        xml_AddCharacters, xml_EndElement, xml_Close
  USE     dom,   ONLY : parseFile, item, getElementsByTagname, &
                        destroy, nodeList, Node
#endif
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  ! definitions for the fmt
  !
  CHARACTER(5),  PARAMETER :: fmt_name = "QEXSD"
  CHARACTER(8),  PARAMETER :: fmt_version = "25.05.21"
  !
  ! internal data to be set
  !
  TYPE(xmlf_t)     :: qexsd_xf
  !
  ! vars to manage back compatibility
  !
  CHARACTER(10)    :: qexsd_current_version = " "
  CHARACTER(10)    :: qexsd_default_version = trim( fmt_version  )
  LOGICAL          :: qexsd_current_version_init = .FALSE.
  !
  TYPE (step_type), ALLOCATABLE                :: steps(:)
  INTEGER                                      :: exit_status
  TYPE ( closed_type )                         :: qexsd_closed_element
  INTEGER                                      :: step_counter
  !
  ! private vars used to set the list of clocks to be saved in the xml file
  CHARACTER(:),dimension(:), ALLOCATABLE     :: clock_list
  INTEGER                                    :: clock_list_dim =0
  INTEGER                                    :: clock_list_last=0 
  !
  ! end of declarations
  !
  PUBLIC :: qexsd_xf  
  PUBLIC :: qexsd_openschema, qexsd_closeschema
  PUBLIC :: qexsd_readschema
  PUBLIC :: qexsd_step_addstep, qexsd_reset_steps
  PUBLIC :: qexsd_current_version, qexsd_default_version, qexsd_current_version_init
  PUBLIC :: qexsd_set_status
  PUBLIC :: qexsd_allocate_clock_list
  PUBLIC :: qexsd_add_label
  PUBLIC :: qexsd_add_all_clocks 
  ! 
CONTAINS
!
!-------------------------------------------
! ... basic subroutines
!-------------------------------------------
!
    !
    !-------------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_set_status(status_int)
    !-------------------------------------------------------------------------------------------------
    INTEGER, INTENT(IN)      :: status_int
    !
    exit_status = status_int
    !
    END SUBROUTINE qexsd_set_status
    !
!
!-------------------------------------------
! ... subroutine writing header, general, parallel info to file
!-------------------------------------------
!
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_openschema(filename, ounit, prog, title)
      !------------------------------------------------------------------------
      !
      USE qexsd_input, ONLY: qexsd_input_obj
      !
      CHARACTER(len=*), INTENT(IN) :: filename, prog, title
      INTEGER, INTENT(IN)          :: ounit
      TYPE (general_info_type)  :: general_info
      TYPE (parallel_info_type) :: parallel_info
      CHARACTER(len=16) :: subname = 'qexsd_openschema'
      INTEGER :: ierr, len_steps, i_step
      !
      ! we need a qes-version number here
      CALL xml_OpenFile(FILENAME = TRIM(filename), XF = qexsd_xf, UNIT = ounit,&
              PRETTY_PRINT = .TRUE., REPLACE  = .TRUE., NAMESPACE = .TRUE., &
              IOSTAT = ierr ) 
      !
      CALL xml_DeclareNamespace (XF=qexsd_xf, PREFIX = "xsi", nsURI ="http://www.w3.org/2001/XMLSchema-instance")
      CALL xml_DeclareNamespace (XF=qexsd_xf, PREFIX = "qes", nsURI ="http://www.quantum-espresso.org/ns/qes/qes-1.0")
      CALL xml_NewElement (XF=qexsd_xf, NAME = "qes:espresso")
      CALL xml_addAttribute(XF=qexsd_xf, NAME = "xsi:schemaLocation", &
                            VALUE = "http://www.quantum-espresso.org/ns/qes/qes-1.0 "//&
                                    "http://www.quantum-espresso.org/ns/qes/qes_250521.xsd" )
      CALL xml_addAttribute(XF=qexsd_xf, NAME="Units", VALUE="Hartree atomic units")
      CALL xml_addComment(XF = qexsd_xf, &
              COMMENT = "All quantities are in Hartree atomic units unless otherwise specified" ) 
      !
      IF (ierr /= 0) call errore(subname, 'opening xml output file', ierr)
      ! the input file is mandatory to have a validating schema 
      ! here an error should be issued, instead
      !
      CALL qexsd_init_general_info(general_info, prog(1:2), title )
      CALL qes_write (qexsd_xf,general_info)
      CALL qes_reset (general_info)
      !
      CALL qexsd_init_parallel_info(parallel_info)
      CALL qes_write (qexsd_xf,parallel_info)
      CALL qes_reset (parallel_info) 
      IF ( check_file_exst(input_xml_schema_file) )  THEN
         CALL xml_addComment( XF = qexsd_xf, COMMENT= "")
#if defined(__fox)
         CALL qexsd_cp_line_by_line(ounit ,input_xml_schema_file, spec_tag="input")
#else
         CALL qexsd_cp_line_by_line(qexsd_xf%unit,input_xml_schema_file, spec_tag="input")
#endif 
      ELSE IF ( TRIM(qexsd_input_obj%tagname) == "input") THEN 
         CALL qes_write (qexsd_xf, qexsd_input_obj)
      END IF
      ! 
      IF (ALLOCATED(steps) ) THEN 
         len_steps= step_counter 
         IF (TRIM (steps(1)%tagname ) .EQ. 'step') THEN
            DO i_step = 1, len_steps
               CALL qes_write (qexsd_xf, steps(i_step) )
            END DO 
         END IF
      END IF
      ! 
    END SUBROUTINE qexsd_openschema
    !
    !
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_general_info(obj, prog, title )
    !---------------------------------------------------------------------------------------
      TYPE( general_info_type )         ::  obj
      CHARACTER(LEN=*),INTENT(IN)       ::  prog
      CHARACTER(LEN=*),INTENT(IN)       ::  title
      CHARACTER(LEN=*),PARAMETER        ::  TAGNAME="general_info"
      TYPE( creator_type )              ::  creator_obj
      TYPE( created_type )              ::  created_obj
      TYPE( xml_format_type)            ::  xml_fmt_obj
      CHARACTER(LEN=256)                ::  version
      CHARACTER(9)                      ::  cdate, ctime
      CHARACTER(60)                     ::  timestamp
      !
      version=TRIM(version_number)
      SELECT CASE( prog(1:2))
         CASE ('pw','PW') 
            CALL qes_init (creator_obj, "creator", "PWSCF", version, "XML file generated by PWSCF")
         CASE ('cp', 'CP') 
            CALL qes_init (creator_obj, "creator", "CP", version, "XML file generated by CP") 
      END SELECT
      !
      CALL date_and_tim(cdate, ctime) 
      timestamp = 'This run was terminated on:  ' // ctime // ' ' // cdate(1:2) // & 
                  ' '//cdate(3:5) // ' '// cdate (6:9)
     
      CALL qes_init (created_obj, "created", cdate, ctime, timestamp ) 
      !
      CALL qes_init (xml_fmt_obj, "xml_format", fmt_name, fmt_version, fmt_name//"_"//fmt_version)
      !
      CALL qes_init ( obj, TAGNAME, XML_FORMAT = xml_fmt_obj, CREATOR = creator_obj, CREATED = created_obj, &
                      JOB=title)
      !
      CALL qes_reset (creator_obj)
      CALL qes_reset (created_obj)
      CALL qes_reset (xml_fmt_obj) 
    END SUBROUTINE qexsd_init_general_info    
    !
    !---------------------------------------------------------------------------------------------
    SUBROUTINE   qexsd_init_parallel_info(obj)
    !---------------------------------------------------------------------------------------------
      TYPE ( parallel_info_type )           :: obj
      !
      INTEGER                               :: nthreads=1
#if defined(__OMP) 
      INTEGER,EXTERNAL                      :: omp_get_max
      !     
      nthreads = omp_get_max()
#endif      
      CALL qes_init (obj, "parallel_info", nproc, nthreads, ntask_groups, &
                                  nbgrp, npool, nproc_bgrp)
    END SUBROUTINE qexsd_init_parallel_info
    !
!
!-------------------------------------------
! ... subroutine writing status and timing info to file and closing it
!-------------------------------------------
!
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_closeschema()
      !------------------------------------------------------------------------
      USE mytime,    ONLY: nclock, clock_label
      REAL(DP),EXTERNAL     :: get_clock
      TYPE(timing_type) :: qexsd_timing_  
      !
      CHARACTER(len=17) :: subname = 'qexsd_closeschema'
      INTEGER :: ierr
      !
      IF (exit_status .ge. 0 ) THEN 
         CALL xml_NewElement(qexsd_xf, "exit_status")
         CALL xml_AddCharacters(qexsd_xf, exit_status)
         CALL xml_EndElement(qexsd_xf, "exit_status")          
         CALL qexsd_set_closed()
         IF (get_clock('PWSCF') > get_clock('CP'))  THEN 
            CALL qexsd_init_clocks (qexsd_timing_, 'PWSCF       ' , clock_list)
         ELSE 
            CALL qexsd_init_clocks (qexsd_timing_, 'CP          ', clock_list) 
         END IF 
         CALL qes_write ( qexsd_xf, qexsd_timing_) 
         CALL qes_reset(qexsd_timing_) 
         !CALL xml_NewElement (qexsd_xf, "cputime")
         !CALL xml_addCharacters(qexsd_xf, MAX(nint(get_clock('PWSCF')),nint(get_clock('CP'))) )
         !CALL xml_EndElement ( qexsd_xf, "cputime")
         CALL qes_write (qexsd_xf, qexsd_closed_element)
      END IF
         CALL xml_Close(qexsd_xf) 
      !
    END SUBROUTINE qexsd_closeschema
    !
    !
!-------------------------------------------
! ... function reading xml and storing inofo into objects
!-------------------------------------------
!
!------------------------------------------------------------------------
    SUBROUTINE qexsd_readschema (filename, ierr, output_obj, parinfo_obj, &
         geninfo_obj, input_obj)
!------------------------------------------------------------------------
      !
      USE qes_read_module, ONLY : qes_read
      !
      CHARACTER(LEN=*), INTENT(IN) :: filename
      INTEGER, INTENT(OUT)         :: ierr
      TYPE( output_type ), OPTIONAL,       INTENT(OUT)   :: output_obj
      TYPE(parallel_info_type), OPTIONAL,  INTENT(OUT)   :: parinfo_obj
      TYPE(general_info_type ), OPTIONAL,  INTENT(OUT)   :: geninfo_obj
      TYPE(input_type), OPTIONAL,          INTENT(OUT)   :: input_obj
      ! 
      TYPE(Node), POINTER     :: root, nodePointer
      TYPE(nodeList),POINTER  :: listPointer
      LOGICAL                 :: found
      CHARACTER(LEN=80)       :: errmsg = ' '
      CHARACTER(len=17)       :: subname = 'qexsd_readschema'
      ! 
      ierr = 0
      ! 
      INQUIRE ( file=filename, exist=found )
      IF (.NOT. found ) THEN
         ierr = 1
         errmsg='xml data file ' // TRIM(filename) // ' not found'
         GOTO 100
      END IF
      !
      ! read XML file into "root" object
      !
      root => parseFile(filename)
      !
      ! copy from "root" object into geninfo, parinfo, output objs
      !
      IF ( PRESENT ( geninfo_obj ) ) THEN 
         nodePointer => item ( getElementsByTagname(root, "general_info"),0)
         IF (ASSOCIATED(nodePointer)) THEN
            CALL qes_read( nodePointer, geninfo_obj, ierr)
         ELSE
            ierr = 2
         END IF
         IF ( ierr /= 0 ) THEN
            errmsg='error reading header of xml data file'
            ierr = 2
            GOTO 100
         END IF
      END IF 
      ! 
      IF ( PRESENT ( parinfo_obj ) ) THEN 
         nodePointer => item ( getElementsByTagname(root,"parallel_info"),0)
         IF (ASSOCIATED(nodePointer)) THEN
            CALL qes_read(nodePointer, parinfo_obj, ierr)
         ELSE
            ierr = 3
         END IF
         IF ( ierr /= 0 ) THEN  
            errmsg='error in parallel_info  of xsd data file' 
            ierr = 3
            GOTO 100
         END IF
      END IF  
      ! 
      IF ( PRESENT ( output_obj ) ) THEN
         nodePointer => item ( getElementsByTagname(root, "output"),0)
         IF (ASSOCIATED(nodePointer)) THEN
            CALL qes_read ( nodePointer, output_obj, ierr ) 
         ELSE
            ierr = 4
         END IF
         IF ( ierr /= 0 ) THEN  
            errmsg = 'error reading output_obj of xsd data file' 
            ierr = 4
            GOTO 100 
         END IF 
      END IF 
      !
      IF (PRESENT (input_obj)) THEN
         nodePointer => item( getElementsByTagname(root, "input"),0)
         IF ( ASSOCIATED(nodePointer) ) THEN
            CALL qes_read (nodePointer, input_obj, ierr ) 
         ELSE 
            ierr =-1
         END IF
         IF ( ierr /= 0 ) THEN
            errmsg = 'input info not found or not readable in xml file'
            IF ( TRIM(input_obj%tagname) == 'input' )  CALL qes_reset (input_obj)
            ierr =-1
         END IF
      END IF
      ! 
      CALL destroy(root)       
      !
 100  IF ( ierr /= 0 ) CALL infomsg(subname,TRIM(errmsg))
      !
    END SUBROUTINE qexsd_readschema
!
!-------------------------------------------
! ... utilities
!-------------------------------------------
!
    !
    !------------------------------------------------------------------------
    FUNCTION check_file_exst( filename )
      !------------------------------------------------------------------------
      !
      LOGICAL          :: check_file_exst
      CHARACTER(len=*) :: filename
      !
      LOGICAL :: lexists
      !
      INQUIRE( FILE = trim( filename ), EXIST = lexists )
      !
      check_file_exst = lexists
      RETURN
      !
    END FUNCTION check_file_exst
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_cp_line_by_line(iun_out,filename,spec_tag)
      !------------------------------------------------------------------------
      implicit none
      !
      integer,      intent(in) :: iun_out
      character(*), intent(in) :: filename
      character(*), optional, intent(in) :: spec_tag
      !
      integer :: iun, ierr
      character(256) :: str
      logical :: icopy, exists

      INQUIRE(FILE=trim(filename), EXIST=exists)
      !
      IF(.not.exists) THEN
         CALL errore('qexsd_cp_line_by_line', 'input xml file "' // & 
        &             TRIM(filename) // '" not found', 1)
      ENDIF
      !
      open(NEWUNIT=iun,FILE=trim(filename),status="old", IOSTAT=ierr)
      !
      icopy=.false.
      copy_loop: do
         !
         read(iun,"(a256)",iostat=ierr) str
         if (ierr<0) exit copy_loop
         if (present(spec_tag)) then
            !
            if (index(str,"<"//trim(adjustl(spec_tag))//">")/=0) then
               !
               icopy=.true.
               !
            endif
            !
         else
            !
            icopy=.true.
            !
         endif
         ! 
         ! filtering
         ! 
         if ( index(str,"<Root>")/=0 .or. index(str,"<Root>")/=0 .or. &
              index(str,"<?")/=0     .or. .not.icopy) cycle copy_loop
         !
         write(iun_out,"(a)") trim(str)
         !
         if (present(spec_tag)) then
            if (index(str,"</input>")/=0) icopy=.false.
         endif
         ! 
      enddo copy_loop
      !
      close(iun)
      ! 
    END SUBROUTINE qexsd_cp_line_by_line
    !
!
!-------------------------------------------
! ... subroutine related to MD steps
!-------------------------------------------
!
    ! 
    !----------------------------------------------------------------------------------------
    SUBROUTINE qexsd_step_addstep(i_step, max_steps, ntyp, atm, ityp, nat, tau, alat, a1, a2, a3, &
                                  etot, eband, ehart, vtxc, etxc, ewald, degauss, demet, forces,  &
                                  stress, scf_has_converged, n_scf_steps, scf_error, efieldcorr, potstat_contr,      &
                                  fcp_force, fcp_tot_charge, gatefield_en)
    !-----------------------------------------------------------------------------------------
    !! This routing initializes le steps array containing up to max_steps elements of the step_type
    !! data structure. Each element contains structural and energetic info for m.d. trajectories and 
    !! structural minimization paths. All quantities must be provided directly in Hartree atomic units. 
    !! @Note updated on April 10th 2018 by Pietro Delugas
    USE qexsd_init, ONLY : qexsd_init_atomic_structure, qexsd_init_total_energy
    USE control_flags, ONLY : tstress
    ! 
    INTEGER ,INTENT(IN)             :: i_step, max_steps, ntyp, nat, n_scf_steps, ityp(:)
    REAL(DP),INTENT(IN)             :: tau(3,nat), alat, a1(3), a2(3), a3(3), etot, eband, ehart, vtxc, &
                                       etxc, ewald, scf_error, forces(3,nat), stress(3,3) 
    LOGICAL,INTENT(IN)              :: scf_has_converged 
    REAL(DP),OPTIONAL,INTENT(IN)    :: degauss, demet, gatefield_en, efieldcorr
    REAL(DP),OPTIONAL,INTENT (IN)   :: potstat_contr, fcp_force, fcp_tot_charge       
    CHARACTER(LEN=*),INTENT(IN)     :: atm(:)
    TYPE (step_type)                :: step_obj
    TYPE ( scf_conv_type )          :: scf_conv_obj
    TYPE ( atomic_structure_type )  :: atomic_struct_obj
    TYPE ( total_energy_type )      :: tot_en_obj
    TYPE ( matrix_type )            :: mat_forces, mat_stress  
    !    
    IF ( i_step .EQ. 1 ) THEN 
       ALLOCATE (steps(max_steps))
       step_counter = 0
    END IF 
    step_counter = step_counter+1
    !
    step_obj%tagname="step"
    step_obj%n_step = i_step 
    !
    CALL qes_init( scf_conv_obj,"scf_conv", scf_has_converged, n_scf_steps, scf_error )
    !
    step_obj%scf_conv = scf_conv_obj 
    CALL qes_reset(scf_conv_obj)
    ! 
    CALL qexsd_init_atomic_structure(atomic_struct_obj, ntyp, atm, ityp, nat, tau, &
                                     alat, a1, a2, a3, 0)
    step_obj%atomic_structure=atomic_struct_obj
    CALL qes_reset( atomic_struct_obj )
    ! 
    CALL qexsd_init_total_energy (tot_en_obj, etot, eband, ehart, &
          vtxc, etxc, ewald, degauss, demet, efieldcorr, potstat_contr, gatefield_en)  
    step_obj%total_energy=tot_en_obj
    CALL qes_reset( tot_en_obj )
    ! 
    CALL  qes_init( mat_forces, "forces", [3, nat], forces ) 
    step_obj%forces=mat_forces
    CALL qes_reset ( mat_forces )
    ! 
    CALL qes_init( mat_stress, "stress", [3, 3], stress ) 
    step_obj%stress = mat_stress
    IF( tstress ) THEN
       step_obj%stress_ispresent = .TRUE.
    END IF
    CALL qes_reset ( mat_stress ) 
    IF ( PRESENT ( fcp_force ) ) THEN 
       step_obj%FCP_force = fcp_force
       step_obj%FCP_force_ispresent = .TRUE.
    END IF 
    IF (PRESENT( fcp_tot_charge)) THEN 
       step_obj%FCP_tot_charge = fcp_tot_charge
       step_obj%FCP_tot_charge_ispresent = .TRUE. 
    END IF 
    !  
    ! 
    steps(step_counter) = step_obj
    steps(step_counter)%lwrite  = .TRUE.
    steps(step_counter)%lread   = .TRUE. 
    call qes_reset (step_obj)
    END SUBROUTINE qexsd_step_addstep 
    !
    !------------------------------------------------------------------------------------
    SUBROUTINE qexsd_reset_steps()
       INTEGER  :: i_step
       IF (ALLOCATED(steps)) THEN
          DO i_step =1, SIZE(steps) 
            CALL qes_reset(steps(i_step))
          END DO
          DEALLOCATE (steps)
      END IF
   END SUBROUTINE qexsd_reset_steps
    !
    !--------------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_set_closed() 
    ! 
    CHARACTER(LEN=9)                  :: cdate, time_string
    CHARACTER(LEN=12)                 :: date_string
    !
    CALL date_and_tim( cdate, time_string ) 
    date_string = cdate(1:2) // ' ' // cdate(3:5) // ' ' // cdate (6:9)
    CALL qes_init (qexsd_closed_element, "closed", date_string, time_string,&
                          "")
    END SUBROUTINE qexsd_set_closed 
    
!
!-------------------------------------------
! ... subroutine related to timing information
!-------------------------------------------
!
    ! 

SUBROUTINE qexsd_init_clocks (timing_, total_clock, partial_clocks)
      !
      USE mytime,  ONLY: nclock, clock_label, cputime, walltime, called
      !
      TYPE(timing_type),INTENT(INOUT)          :: timing_ 
      CHARACTER(LEN=*),INTENT(IN)             :: total_clock 
      CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: partial_clocks(:) 
      ! 
      TYPE (clock_type)                 :: total_
      TYPE(clock_type),ALLOCATABLE      :: partial_(:)  
      LOGICAL,ALLOCATABLE               :: match(:)
      INTEGER                           :: partial_ndim = 0, ic, ipar, nc 
      REAL (DP)                         :: t(2)
      INTERFACE
         FUNCTION get_cpu_and_wall(n_) result(t_)
            IMPORT :: DP
            IMPLICIT NONE
            INTEGER    :: n_ 
            REAL(DP)   t_(2)
         END FUNCTION get_cpu_and_wall 
      END INTERFACE
      ! 
      IF (PRESENT(partial_clocks)) partial_ndim = clock_list_last
      DO ic = 1, nclock
         IF ( TRIM(total_clock) == clock_label(ic) ) EXIT 
      END DO 
      t = get_cpu_and_wall(ic) 
      CALL qes_init ( total_, "total", TRIM(clock_label(ic)), CPU = t(1), WALL = t(2) ) 
      IF ( partial_ndim .GT.  0 ) THEN  
         ALLOCATE(partial_(partial_ndim), match(nclock) ) 
         DO ipar = 1, partial_ndim 
            match = clock_label(1:nclock) == TRIM(partial_clocks(ipar)) 
            IF ( ANY (match))  THEN
               nc = get_index(.TRUE., match)
               IF (nc == ic .OR. called(nc) == 0 ) CYCLE
               t = get_cpu_and_wall(nc) 
               CALL qes_init(partial_(ipar), "partial", TRIM(clock_label(nc)), CPU = t(1), WALL = t(2), & 
                              CALLS = called(nc))
            ELSE 
               CALL qes_init (partial_(ipar), "partial", "not_found",  CPU = -1.d0, WALL = -1.d0, CALLS = 0)  
               partial_(ipar)%lwrite=.FALSE. 
            END IF 
         END DO
      END IF 
      CALL qes_init( timing_, "timing_info", total_, partial_)
      CALL qes_reset ( total_) 
      DO ipar =1, partial_ndim 
         CALL qes_reset(partial_(ipar)) 
      END DO
      CONTAINS 
         FUNCTION get_index(val, array)  result(n) 
            IMPLICIT NONE
            LOGICAL                     :: val 
            LOGICAL                     :: array(:)
            INTEGER                     :: n 
            ! 
            INTEGER                     :: i 
            !
            n = - 1
            DO i =1, SIZE(array) 
               IF (array(i) .EQV. val) EXIT 
            END DO
            IF ( array(i) .EQV. val )  n = i 
         END FUNCTION get_index 
   END SUBROUTINE qexsd_init_clocks  

   SUBROUTINE qexsd_allocate_clock_list(prog)
     !! allocates the list of clock labels    
     CHARACTER(*), INTENT(IN) :: prog
     !! name of the program
     IF (ALLOCATED(clock_list)) DEALLOCATE (clock_list) 
     IF (prog == 'PW') THEN 
      ALLOCATE(character(len=32) :: clock_list(100))
      clock_list_dim = 100 
     ELSE IF (prog == 'CPV') THEN 
      ALLOCATE (character(len=32) :: clock_list(100))
      clock_list_dim = 100 
     END IF 
   END SUBROUTINE qexsd_allocate_clock_list

   SUBROUTINE qexsd_add_all_clocks()
     !! allocates the list of clock labels copying all active clocks
     USE mytime,  ONLY: nclock, clock_label
     IF (ALLOCATED(clock_list))  DEALLOCATE (clock_list) 
     ALLOCATE (clock_list, SOURCE=clock_label(1:nclock)) 
     clock_list_dim = nclock
     clock_list_last = nclock
    
   END SUBROUTINE qexsd_add_all_clocks

   SUBROUTINE qexsd_add_label (label)
      !! adds a clock label to the clock list that will be reported in the xml file
      CHARACTER(LEN=*),INTENT(IN)  :: label
      !! clock label to be added to the list 
      !
      IF (clock_list_dim == 0) THEN 
         CALL infomsg("qexsd_add_label:", "trying to add label before allocation FIXME")
         RETURN
      END IF
      IF ( clock_list_last .GE. clock_list_dim) THEN 
        CALL infomsg("qexsd_add_label:", "too many clocks FIXME")
        RETURN
      END IF 
      clock_list(clock_list_last+1) = label
      clock_list_last = clock_list_last + 1 
   END SUBROUTINE qexsd_add_label
   !
END MODULE qexsd_module
