
! Copyright (C) 2003-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module contains routines for reading XML files written according to the
! Quantum Espresso XML input/output schema (QEXSD). The routines use the iotk
! library and all assume that the iotk library has been initialized in the calling wrapper. 
! This module should be, in the future, replaced by a script generated library. 
! 
! Written by Pietro Delugas -- May-June 2016. 
!--------------------------------------------------------------------------------------
MODULE qexsd_reader_module
!--------------------------------------------------------------------------------------
USE   kinds,    ONLY: DP
USE   io_global,ONLY: stderr =>stdout
USE   iotk_module
USE   qes_module 
!
IMPLICIT NONE
!
CONTAINS 
!-------------------------------------------------------------------
SUBROUTINE qexsd_get_general_info(iunit, obj, general_info_ispresent)
!-------------------------------------------------------------------
IMPLICIT NONE 
!
INTEGER, INTENT(IN)                      :: iunit
TYPE(general_info_type),INTENT(OUT)      :: obj
LOGICAL,INTENT(OUT)                      :: general_info_ispresent     
!
!
TYPE (xml_format_type)                   :: xml_fmt_obj
TYPE ( creator_type )                    :: creator_obj
TYPE ( created_type )                    :: created_obj
INTEGER                                  :: ierr
CHARACTER(256)                           :: dummy, name_, version_, creator_, xml_format_, date_, &
                                            time_, job_, created_
LOGICAL                                  :: found, xml_format_ispresent, &
                                            creator_ispresent, created_ispresent, job_ispresent 

CHARACTER(iotk_attlenx)                  :: attr


CALL iotk_scan_begin(iunit, "general_info" ,IERR = ierr , FOUND = general_info_ispresent )
IF(ierr /= 0) RETURN 
   !
   CALL iotk_scan_dat  (iunit, "xml_format", dat=xml_format_, IERR = ierr ) 
   IF (ierr /=0 ) RETURN
   CALL iotk_scan_begin(iunit, "xml_format", ATTR = attr, IERR = ierr , FOUND = xml_format_ispresent )
   IF ( (ierr /=0) ) RETURN 
      CALL iotk_scan_attr(attr, "NAME", name_, IERR = ierr )
      IF (ierr /= 0) RETURN 
      CALL iotk_scan_attr(attr, "VERSION", version_, IERR = ierr) 
      IF (ierr /= 0) RETURN 
   CALL iotk_scan_end(iunit,"xml_format", IERR = ierr) 
   IF (ierr /=0) RETURN 
   CALL qes_init_xml_format(xml_fmt_obj, "xml_format", name_, version_, xml_format_)
   ! 
   CALL iotk_scan_dat(iunit, "creator", dat=creator_, IERR = ierr )
   IF (ierr /= 0) RETURN 
   CALL iotk_scan_begin(iunit, "creator",ATTR=attr, FOUND = creator_ispresent)
   IF (ierr /=0 ) RETURN 
      CALL iotk_scan_attr( attr, "NAME", name_, IERR=ierr ) 
      CALL iotk_scan_attr( attr, "VERSION", version_, IERR=ierr )
   CALL iotk_scan_end(iunit, "creator", IERR = ierr )
   CALL qes_init_creator(creator_obj, "creator", name_, version_, creator_)  
   ! 
   CALL iotk_scan_dat(iunit, "created", created_, IERR = ierr )
   IF (ierr /=0) RETURN
   CALL iotk_scan_begin(iunit, "created", ATTR = attr, FOUND = created_ispresent ,IERR = ierr)
   IF (ierr /=0) RETURN 
      CALL iotk_scan_attr(attr,"DATE",date_,IERR = ierr)
      CALL iotk_scan_attr(attr,"TIME", time_, IERR = ierr )
   CALL iotk_scan_end(iunit,"created", IERR = ierr) 
   IF (ierr /= 0) RETURN  
   CALL qes_init_created( created_obj, "created", date_, time_, created_)
   CALL iotk_scan_dat(iunit, "job", job_, IERR = ierr ) 
CALL iotk_scan_end(iunit, "general_info", IERR = ierr ) 
IF(ierr /=0) RETURN
CALL qes_init_general_info(obj,"general_info", xml_fmt_obj, creator_obj, created_obj, job_) 
CALL qes_reset_xml_format(xml_fmt_obj)
CALL qes_reset_creator(creator_obj)
CALL qes_reset_created(created_obj)
! 
END SUBROUTINE qexsd_get_general_info
!
!---------------------------------------------------------------------------------------------- 
SUBROUTINE qexsd_get_parallel_info(iunit, obj, parallel_info_ispresent)
!
IMPLICIT NONE 
INTEGER,INTENT(IN)                  :: iunit
TYPE (parallel_info_type),INTENT(out)    :: obj 
LOGICAL,INTENT(OUT)                 :: parallel_info_ispresent
! 

INTEGER                             :: nprocs_, nthreads_, ntasks_, nbgrp_, npool_, ndiag_
INTEGER                             :: ierr, sum_err
CHARACTER(iotk_attlenx)             :: attr
!
parallel_info_ispresent = .FALSE.
CALL iotk_scan_begin(iunit, "parallel_info", IERR = ierr, FOUND = parallel_info_ispresent ) 
IF (ierr /=0 .OR. .NOT. parallel_info_ispresent) RETURN 
! 
sum_err = 0 
CALL iotk_scan_dat(iunit, "nprocs", nprocs_, IERR = ierr ) 
sum_err = sum_err + ierr
CALL iotk_scan_dat(iunit, "nthreads", nthreads_, IERR = ierr )
sum_err = sum_err + ierr
CALL iotk_scan_dat(iunit, "ntasks", ntasks_, IERR = ierr )
sum_err = sum_err + ierr 
CALL iotk_scan_dat(iunit, "nbgrp", nbgrp_, IERR = ierr )
sum_err = sum_err + ierr
CALL iotk_scan_dat(iunit, "npool", npool_, IERR = ierr )
sum_err = sum_err + ierr
CALL iotk_scan_dat(iunit, "ndiag", ndiag_, IERR = ierr )
sum_err = sum_err + ierr
IF (sum_err/=0) RETURN
CALL iotk_scan_end(iunit, "parallel_info", IERR = ierr ) 
IF (ierr /=0) RETURN

CALL qes_init_parallel_info(obj, "parallel_info", nprocs_, nthreads_, ntasks_, nbgrp_, npool_, ndiag_)
END SUBROUTINE qexsd_get_parallel_info
!
!---------------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_control_variables(iunit, obj, control_variables_ispresent)
!
IMPLICIT NONE 
!   
INTEGER,INTENT(IN)                            :: iunit 
TYPE (control_variables_type),INTENT(OUT)     :: obj
LOGICAL,INTENT(OUT)                           :: control_variables_ispresent
! 
INTEGER          :: ierr 
CHARACTER(256)   :: buffer, title_, calculation_, prefix_, pseudo_dir_, &
                    disk_io_, verbosity_, outdir_, restart_
LOGICAL          :: found, stress_, forces_, wf_collect_, nstep_ispresent
INTEGER          :: max_seconds_, print_every_, nstep
REAL(DP)         :: etot_conv_thr_, forc_conv_thr_, press_conv_thr_   
CHARACTER(iotk_attlenx)         :: attr
! 
control_variables_ispresent = .FALSE.
CALL iotk_scan_begin(iunit, "control_variables", IERR = ierr, FOUND = control_variables_ispresent )
IF (ierr /= 0 .OR. (.NOT. control_variables_ispresent ) ) RETURN 
!
CALL iotk_scan_dat(iunit, "title", dat = title_, IERR = ierr) 
IF (ierr /= 0) RETURN
!
CALL iotk_scan_dat(iunit, "calculation", dat = calculation_, FOUND = found ) 
IF ( (.NOT. found) .OR. ( TRIM(calculation_)== "")) calculation_ = "scf"
!
CALL iotk_scan_dat (iunit, "restart_mode", dat = restart_, FOUND = found ) 
IF ( (.NOT. found) .OR. ( TRIM(restart_)== "")) restart_ = "from_scratch" 
! 
CALL iotk_scan_dat (iunit, "prefix", dat = prefix_, FOUND = found ) 
IF ( (.NOT. found) .OR. ( TRIM(prefix_)== "")) prefix_ = "pwscf"
! 
CALL iotk_scan_dat ( iunit, "pseudo_dir", dat = pseudo_dir_, FOUND = found)
IF ( (.NOT. found) .OR. ( TRIM(pseudo_dir_)== "")) pseudo_dir_ = "./"
! 
CALL iotk_scan_dat ( iunit, "outdir", dat = outdir_, FOUND = found)
IF ( (.NOT. found) .OR. ( TRIM(outdir_)== "")) outdir_ = "./"
!
stress_ = .FALSE.  
CALL iotk_scan_dat ( iunit, "stress", dat = buffer, FOUND = found)
IF ( found ) THEN
   SELECT CASE (TRIM(buffer))
      CASE ("true") 
         stress_ = .TRUE.
      CASE ("false") 
          stress_ = .FALSE.
   END SELECT
END IF 
! 
forces_ = .FALSE. 
CALL iotk_scan_dat ( iunit, "forces", dat = buffer, FOUND = found ) 
IF ( found ) THEN
   SELECT CASE (TRIM(buffer))
      CASE ("true")
         forces_ = .TRUE.
      CASE ("false")
          forces_ = .FALSE.
   END SELECT
END IF
! 
wf_collect_ = .FALSE. 
CALL iotk_scan_dat ( iunit, "wf_collect", dat = buffer, FOUND = found )
IF ( found ) THEN
   SELECT CASE (TRIM(buffer))
      CASE ("true")
         wf_collect_ = .TRUE.
      CASE ("false")
          wf_collect_ = .FALSE.
   END SELECT
END IF 
! 
CALL iotk_scan_dat (iunit, "disk_io", dat = disk_io_, FOUND = found ) 
IF ( (.NOT. found ) .OR. TRIM(disk_io_) == "") disk_io_ = "low"
! 
CALL iotk_scan_dat (iunit, "max_seconds", dat = max_seconds_, FOUND = found, IERR = ierr ) 
IF ((.NOT. found ) .OR. ierr /= 0) RETURN 
! 
CALL iotk_scan_dat (iunit, "etot_conv_thr", dat = etot_conv_thr_, FOUND = found, IERR = ierr  )
IF (.NOT. found)  etot_conv_thr_=1.d-5
! 
CALL iotk_scan_dat (iunit, "forc_conv_thr", dat = forc_conv_thr_, FOUND = found, IERR = ierr  )
IF (.NOT. found ) THEN 
      IF ( forces_ ) forc_conv_thr_=1.d-3
END IF 
! 
CALL iotk_scan_dat ( iunit, "press_conv_thr", dat = press_conv_thr_, FOUND = found, IERR = ierr ) 
IF (.NOT. found ) THEN 
      IF ( stress_ ) press_conv_thr_=forc_conv_thr_
END IF
! 
CALL iotk_scan_dat ( iunit, "verbosity", dat = verbosity_, FOUND = found, IERR = ierr )
IF ((.NOT. found) .OR. TRIM(verbosity_) == "" ) verbosity_="low"
! 
CALL iotk_scan_dat ( iunit, "print_every", dat = print_every_, FOUND = found, IERR = ierr )
IF ( (.NOT. found ) .OR. (ierr /=0) ) RETURN 
!
CALL iotk_scan_dat ( iunit, "nstep", dat = nstep, FOUND = nstep_ispresent, IERR = ierr ) 
CALL iotk_scan_end( iunit, "control_variables", IERR = ierr ) 
IF (ierr /=0) RETURN    
CALL qes_init_control_variables(obj, "control_variables", title_, calculation_, restart_, &
                                prefix_, pseudo_dir_, outdir_, stress_, forces_, wf_collect_,  &
                                disk_io_, max_seconds_, nstep_ispresent, nstep, etot_conv_thr_, &
                                forc_conv_thr_, press_conv_thr_, verbosity_, print_every_) 

END SUBROUTINE  qexsd_get_control_variables
!
!---------------------------------------------------------------------------------------------------------- 
SUBROUTINE qexsd_get_atomic_species(iunit, obj, atomic_species_ispresent)
! 
IMPLICIT NONE 
! 
INTEGER,INTENT(IN)                        :: iunit 
TYPE (atomic_species_type),INTENT(OUT)    :: obj
LOGICAL,INTENT(OUT)                       :: atomic_species_ispresent
! 
INTEGER               :: ntyp 
TYPE (species_type),ALLOCATABLE           :: the_species(:)
LOGICAL                                   :: pippo, mass_ispresent, spin_theta_ispresent, &
                                             spin_phi_ispresent, start_mag_ispresent
INTEGER                                   :: ierr, ctyp, t_typ, ntyp_
REAL(DP)                                  :: mass_, start_mag_, spin_phi_, spin_theta_
CHARACTER(LEN=256)                        :: buffer, pseudo_file_, name_spec_
CHARACTER(iotk_attlenx)                   :: attr
!
atomic_species_ispresent =.FALSE.
CALL iotk_scan_begin(iunit, "atomic_species", ATTR = attr,  FOUND = atomic_species_ispresent, &
                                                                                     IERR = ierr ) 
IF (ierr /= 0 ) RETURN 
! 
CALL iotk_scan_attr(attr, "ntyp", ntyp_, IERR = ierr ) 
IF (ierr /=0 ) RETURN 
! 
ALLOCATE (the_species(ntyp_)) 

speciescyc:DO ctyp =1 , ntyp_
   CALL iotk_scan_begin(iunit,"species", ATTR = attr, IERR = ierr ) 
   IF ( ierr /= 0 ) RETURN
 
   CALL iotk_scan_attr(attr,"name",  name_spec_, IERR = ierr ) 
   IF (ierr /= 0 ) RETURN 
   CALL iotk_scan_dat(iunit,"mass", mass_, FOUND = mass_ispresent, IERR = ierr)
   CALL iotk_scan_dat(iunit,"pseudo_file",dat = pseudo_file_, IERR =ierr ) 
   IF (ierr /= 0 ) RETURN
   CALL iotk_scan_dat(iunit, "starting_magnetization", dat = start_mag_, &
                      FOUND = start_mag_ispresent, IERR = ierr)
   CALL iotk_scan_dat(iunit, "spin_teta", dat = spin_theta_,            &
                      FOUND = spin_theta_ispresent, IERR = ierr ) 
   CALL iotk_scan_dat(iunit, "spin_phi", dat = spin_phi_, FOUND = spin_phi_ispresent, & 
                      IERR = ierr ) 
   CALL iotk_scan_end(iunit, "species", IERR = ierr ) 
   CALL qes_init_species(the_species(ctyp), "species", name_spec_, mass_ispresent, mass_, &
                         pseudo_file_, start_mag_ispresent, start_mag_, spin_theta_ispresent, &
                         spin_theta_, spin_phi_ispresent, spin_phi_)
END DO speciescyc     
! 
CALL iotk_scan_end(iunit, "atomic_species", IERR = ierr ) 
IF (ierr /= 0 ) RETURN    
! 
CALL qes_init_atomic_species(obj, "atomic_species", ntyp_, ntyp_, the_species)
DO t_typ=1, ntyp_
   CALL qes_reset_species(the_species(t_typ))
END DO
DEALLOCATE (the_species)
!
END SUBROUTINE qexsd_get_atomic_species
! 
!-------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_atomic_structure(iunit,obj,atomic_structure_ispresent) 
! 
IMPLICIT NONE 
! 
INTEGER, INTENT(IN)                        :: iunit
TYPE(atomic_structure_type),INTENT(OUT)    :: obj
LOGICAL,INTENT(OUT)                         :: atomic_structure_ispresent
! 
INTEGER                                    :: ierr, sum_err
LOGICAL                                    :: nat_ispresent, alat_ispresent, atomic_positions_ispresent, &
                                              wyckoff_positions_ispresent, index_ispresent, ibrav_ispresent
TYPE ( atomic_positions_type )             :: at_pos_obj
TYPE ( wyckoff_positions_type)             :: wyck_pos_obj
TYPE (cell_type)                           :: cell_obj
TYPE ( atom_type),ALLOCATABLE              :: atoms_obj(:)
REAL(DP)                                   :: alat_, the_coordinates(3), a1_(3), a2_(3), a3_(3)
INTEGER                                    :: iat, nat_,spc_grp_, index_, ibrav_, direction_
CHARACTER(LEN=256)                         :: label_, wyck_label_
CHARACTER(iotk_attlenx)                    :: attr

!  
atomic_structure_ispresent = .FALSE. 
CALL  iotk_scan_begin(iunit, "atomic_structure", ATTR = attr, IERR = ierr, &
                               FOUND = atomic_structure_ispresent ) 
IF ( ierr /= 0 ) CALL errore (  "qexsd_get_atomic_structure", "problem opening atomic_structure",ierr) 
CALL iotk_scan_attr(attr, "nat", nat_, IERR = ierr ) 
IF (ierr /= 0) CALL errore ( "qexsd_get_atomic_structure", "mandatory nat attribure is missing",ierr) 
CALL iotk_scan_attr(attr, "alat", alat_, FOUND = alat_ispresent, IERR = ierr ) 
! 
CALL iotk_scan_attr ( attr,"bravais_index", ibrav_, FOUND = ibrav_ispresent, IERR = ierr )
IF ( ierr /= 0 ) CALL errore ( "qexsd_get_atomic_structure", "weird error reading atomic_structure attributes", ierr) 
! 
ALLOCATE (atoms_obj(nat_))
CALL iotk_scan_begin(iunit, "atomic_positions", ATTR=attr, IERR = ierr , FOUND = atomic_positions_ispresent ) 
IF (atomic_positions_ispresent ) THEN 
  DO iat = 1, nat_ 
     CALL iotk_scan_dat(iunit, "atom", ATTR = attr, dat = the_coordinates, IERR = ierr )
     IF (ierr /= 0) RETURN
     CALL iotk_scan_attr(attr, "name", label_, IERR = ierr ) 
     IF (ierr /= 0) RETURN 
     CALL iotk_scan_attr( attr, "index", index_, IERR = ierr, FOUND = index_ispresent ) 
     IF ( ierr /= 0 ) RETURN
     ! 
     CALL qes_init_atom(atoms_obj(iat), "atom", label_, "", .FALSE. , index_, index_ispresent, &
                         the_coordinates)
  END DO         
  CALL iotk_scan_end(iunit, "atomic_positions", IERR = ierr)
  IF (ierr /=0) RETURN 
  CALL qes_init_atomic_positions( at_pos_obj, "atomic_positions", nat_, atoms_obj)
END IF 
!  
CALL iotk_scan_begin(iunit, "wyckoff_positions", ATTR = attr, FOUND = wyckoff_positions_ispresent, IERR = ierr )
IF (atomic_positions_ispresent .AND. wyckoff_positions_ispresent) RETURN 
!
IF (wyckoff_positions_ispresent) THEN 
   CALL iotk_scan_attr(attr,"space_group", spc_grp_, IERR = ierr) 
   IF (ierr /= 0 ) RETURN  
   DO iat = 1, nat_ 
      CALL iotk_scan_begin(iunit, "atom", ATTR = attr, IERR = ierr) 
      IF (ierr /= 0) RETURN 
      CALL iotk_scan_attr(attr, "name", label_, IERR = ierr ) 
      IF (ierr /= 0 ) RETURN 
      CALL iotk_scan_attr(attr, "position", wyck_label_, IERR = ierr ) 
      IF (ierr /=0) RETURN 
      CALL iotk_scan_dat(iunit,"atom", the_coordinates, IERR = ierr ) 
      IF (ierr /= 0) RETURN 
      CALL iotk_scan_end(iunit, "wyckoff_positions", IERR = ierr)
      IF (ierr /=0) RETURN
      CALL qes_init_atom(atoms_obj(iat), "atom", label_, wyck_label_, .TRUE. , 0, .FALSE., the_coordinates)
   END DO
   CALL iotk_scan_end(iunit, "wyckoff_positions", IERR = ierr)
   IF (ierr /= 0 ) RETURN
   CALL qes_init_wyckoff_positions(wyck_pos_obj , "wyckoff_positions", spc_grp_, "",.FALSE., nat_, atoms_obj)
END IF 
!
IF (.NOT. (atomic_positions_ispresent .OR. wyckoff_positions_ispresent ) ) RETURN 
DO iat = 1, nat_
   CALL qes_reset_atom(atoms_obj(iat)) 
END DO
DEALLOCATE ( atoms_obj) 
CALL iotk_scan_begin(iunit, "cell", ATTR = attr, IERR = ierr ) 
IF (ierr /= 0 ) RETURN 
sum_err = 0 
CALL iotk_scan_dat( iunit, "a1", dat = a1_, IERR = ierr ) 
sum_err = sum_err + ierr 
CALL iotk_scan_dat( iunit, "a2", dat = a2_, IERR = ierr )
sum_err = sum_err + ierr  
CALL iotk_scan_dat( iunit, "a3", dat = a3_, IERR = ierr )
sum_err = sum_err + ierr
!
IF (sum_err /= 0 ) CALL errore ( "qexsd_get_atomic_structure", "some cell vector is missing",1) 
! 
CALL iotk_scan_end(iunit, "cell", IERR = ierr )
IF (ierr /= 0 ) CALL errore ( "qexsd_get_atomic_structure", "problems closing cell element",ierr)
!
CALL qes_init_cell ( cell_obj, "cell", a1_, a2_, a3_)

CALL iotk_scan_end (iunit, "atomic_structure", IERR = ierr ) 
IF ( ierr /=0) CALL errore ( "qexsd_get_atomic_structure", "problems closing atomic_structure",ierr)
! 
CALL qes_init_atomic_structure( obj, "atomic_structure", nat_, alat_, alat_ispresent, ibrav_, ibrav_ispresent,&
                                 atomic_positions_ispresent, at_pos_obj, wyckoff_positions_ispresent, &
                                 wyck_pos_obj, cell_obj) 
IF (atomic_positions_ispresent)  CALL qes_reset_atomic_positions(at_pos_obj)
IF (wyckoff_positions_ispresent) CALL qes_reset_wyckoff_positions(wyck_pos_obj) 
CALL qes_reset_cell (cell_obj) 

END SUBROUTINE qexsd_get_atomic_structure
!
!---------------------------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_hybrid_dft( iunit, obj, ispresent ) 
!--------------------------------------------------------------------------------------------------------------------
!
IMPLICIT NONE 
! 
INTEGER,INTENT(IN)             :: iunit
TYPE (hybrid_type),INTENT(OUT) :: obj
LOGICAL,INTENT(OUT)            :: ispresent
! 
INTEGER                        :: ierr, sum_err  
INTEGER                        :: nqx1_, nqx2_, nqx3_ 
TYPE (qpoint_grid_type)        :: qgrid_obj
REAL(DP)                       :: ecutfock_, exx_fraction_, screening_parameter_, &
                                  ecutvcut_
LOGICAL                        :: x_gamma_extrapolation_
CHARACTER(LEN=256)             :: exxdiv_treatment_
CHARACTER(iotk_attlenx)         :: attr
!
!
!
ispresent = .FALSE.
CALL iotk_scan_begin(iunit,"hybrid", ATTR = attr, IERR = ierr, FOUND = ispresent)  
IF (ierr /= 0 ) RETURN 
IF (.NOT. ispresent ) RETURN 
CALL iotk_scan_begin(iunit, "qpoint_grid", ATTR = attr, IERR = ierr ) 
IF (ierr /=0) RETURN
sum_err = 0 
CALL iotk_scan_attr(attr, "nqx1", nqx1_, IERR = ierr ) 
sum_err = sum_err + ierr 
CALL iotk_scan_attr(attr, "nqx2", nqx2_, IERR = ierr )
sum_err = sum_err + ierr
CALL iotk_scan_attr(attr, "nqx3", nqx3_, IERR = ierr )
sum_err = sum_err + ierr
IF (sum_err /= 0 ) RETURN
CALL iotk_scan_end(iunit, "qpoint_grid", IERR = ierr ) 
CALL qes_init_qpoint_grid(qgrid_obj, "qpoint_grid", nqx1_, nqx2_, nqx3_, "")
IF (ierr /=0) RETURN
sum_err = 0 
CALL iotk_scan_dat(iunit, "ecutfock", ecutfock_, IERR = ierr ) 
sum_err = sum_err + ierr 
CALL iotk_scan_dat(iunit, "exx_fraction", exx_fraction_, IERR = ierr )
sum_err = sum_err + ierr
CALL iotk_scan_dat(iunit, "screening_parameter", screening_parameter_, IERR = ierr )
sum_err = sum_err + ierr
CALL iotk_scan_dat(iunit, "exxdiv_treatment", exxdiv_treatment_, IERR = ierr )
sum_err = sum_err + ierr
CALL iotk_scan_dat(iunit, "x_gamma_extrapolation",  x_gamma_extrapolation_ , IERR = ierr )
sum_err = sum_err + ierr      
CALL iotk_scan_dat(iunit, "ecutvcut",  ecutvcut_ , IERR = ierr )
sum_err = sum_err + ierr
IF (sum_err /= 0 ) RETURN
CALL iotk_scan_end(iunit, "hybrid", IERR = ierr ) 
IF (ierr /=0) RETURN
!
CALL qes_init_hybrid(obj, "hybrid", qgrid_obj, ecutfock_, exx_fraction_, screening_parameter_, &
                     exxdiv_treatment_, x_gamma_extrapolation_, ecutvcut_) 
!
CALL qes_reset_qpoint_grid( qgrid_obj ) 

END SUBROUTINE qexsd_get_hybrid_dft
!
!-----------------------------------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_hubbard_common(iunit, tag, obj, ispresent) 
!--------------------------------------------------------------------------------------------------------------------------
! 
IMPLICIT NONE 
! 
INTEGER,INTENT(IN)                            :: iunit 
TYPE (hubbardCommon_type),INTENT(OUT)         :: obj
CHARACTER(LEN = *), INTENT(IN)                      :: tag
LOGICAL,INTENT(OUT)                           :: ispresent 
!
REAL(DP)                                      :: my_hubbard_val 
INTEGER                                       :: ierr
CHARACTER(LEN=256)                            :: tag_, specie_, label_ 
CHARACTER(iotk_attlenx)         :: attr
! 
ispresent = .FALSE. 
!
IF (len (TRIM ( tag)) .GT. 256  ) CALL errore ( 'qexsd_get_hubbard_u:', 'tag too long!!!',1)  
tag_=TRIM(tag)  
CALL iotk_scan_dat(iunit,  tag_ , my_hubbard_val, ATTR = attr, IERR = ierr, FOUND = ispresent ) 
IF (ierr /=0 ) THEN 
   ispresent = .FALSE.
   RETURN 
END IF 
IF ( .NOT. ispresent ) RETURN 
! 
CALL iotk_scan_attr(attr, "specie",  specie_, IERR = ierr )
IF ( ierr /=0 ) RETURN  
CALL iotk_scan_attr(attr, "label", label_, IERR = ierr ) 
IF (ierr /=0 ) RETURN 
! 
CALL qes_init_hubbardCommon(obj, TRIM(tag_) , TRIM (specie_), TRIM (label_), my_hubbard_val )
! 

END SUBROUTINE qexsd_get_hubbard_common 
! 
!--------------------------------------------------------------------------
SUBROUTINE qexsd_get_hubbard_J(iunit, tag, obj, ispresent)
!-------------------------------------------------------------------------- 
IMPLICIT NONE
! 
INTEGER,INTENT(IN)                    :: iunit
TYPE (hubbardJ_type),INTENT(OUT)      :: obj
CHARACTER(LEN=*),INTENT(IN)           :: tag
LOGICAL,INTENT(OUT)                   :: ispresent
!
REAL(DP)                              :: my_hubbard_J(3) 
INTEGER                               :: ierr
CHARACTER(LEN=256)                    :: tag_, specie_, label_
CHARACTER(iotk_attlenx)               :: attr
!
IF (len (TRIM(tag)) .GT. 256 )   CALL errore ( 'qexsd_get_hubbard_J','tag too long' , 1)   
ispresent = .FALSE.
!
tag_ = TRIM ( tag) 
CALL iotk_scan_dat(iunit, TRIM(tag_), my_hubbard_J, ATTR = attr, IERR = ierr , FOUND = ispresent )
IF (ierr /=0 ) THEN 
   ispresent = .FALSE. 
   RETURN
END IF
IF (.NOT. ispresent ) RETURN 
! 
CALL iotk_scan_attr(attr, "specie",  specie_, IERR = ierr )
IF ( ierr /=0 ) RETURN
CALL iotk_scan_attr(attr, "label", label_, IERR = ierr )
IF (ierr /=0 ) RETURN
! 
CALL qes_init_hubbardJ(obj, TRIM(tag_), TRIM (specie_), TRIM (label_), my_hubbard_J )
! 
END SUBROUTINE qexsd_get_hubbard_J     
! 
!--------------------------------------------------------------------------------------------- 
SUBROUTINE qexsd_get_starting_ns(iunit, obj, lmax, ispresent )  
!----------------------------------------------------------------------------------------------
! 
IMPLICIT NONE 
! 
INTEGER, INTENT(IN)                 :: iunit, lmax
TYPE (starting_ns_type),INTENT(OUT) :: obj
LOGICAL,INTENT(OUT)                 :: ispresent
! 
CHARACTER(256)                      :: label_, specie_
REAL(DP),ALLOCATABLE                :: the_vec(:) 
INTEGER                             :: spin_  
INTEGER                             :: ierr, ndim_vec
CHARACTER(iotk_attlenx)             :: attr
!

ndim_vec = 2*lmax+1
ALLOCATE (the_vec(ndim_vec)) 

ispresent = .FALSE. 
CALL iotk_scan_dat(iunit, "starting_ns", DAT= the_vec, ATTR = attr, IERR = ierr , FOUND = ispresent)  
IF (ierr /=0) THEN 
   ispresent = .FALSE. 
   RETURN 
END IF 
IF ( .NOT. ispresent  ) RETURN 
CALL iotk_scan_attr(attr, "specie", specie_, IERR = ierr ) 
IF (ierr /= 0 ) RETURN 
! 
CALL iotk_scan_attr(attr,"label", label_, IERR = ierr ) 
IF (ierr /= 0 ) RETURN  
!
CALL iotk_scan_attr(attr, "spin", spin_, IERR = ierr ) 
IF (ierr /= 0 ) RETURN 
CALL qes_init_starting_ns(obj, "starting_ns", TRIM(specie_), TRIM(label_), spin_, ndim_vec, the_vec ) 
! 
DEALLOCATE (the_vec) 
!
END SUBROUTINE qexsd_get_starting_ns
!
!------------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_hubbard_ns( iunit, obj, lmax, ispresent )
!----------------------------------------------------------------------------------------------------
! 
IMPLICIT NONE 
! 
INTEGER,INTENT ( IN )                     :: iunit, lmax 
TYPE ( hubbard_ns_type ), INTENT (out)    :: obj
LOGICAL, INTENT(OUT)                      :: ispresent 
! 
INTEGER                                   :: ierr
INTEGER                                   :: ndim_matrix_ns
REAL(DP),ALLOCATABLE                      :: the_matrix(:,:)
CHARACTER(256)                            :: specie_, label_
INTEGER                                   :: spin_, index_
CHARACTER(iotk_attlenx)                   :: attr
!
ndim_matrix_ns = 2*lmax+1 
ALLOCATE (the_matrix(ndim_matrix_ns, ndim_matrix_ns))
ispresent = .FALSE. 
CALL iotk_scan_dat(iunit, "Hubbard_ns", dat = the_matrix,  ATTR = attr, FOUND = ispresent, IERR = ierr )   
IF ( ierr /= 0) RETURN
IF ( .NOT. ispresent ) RETURN
! 
CALL iotk_scan_attr(attr, "specie", specie_, IERR = ierr ) 
IF (ierr /= 0 ) RETURN 
! 
CALL iotk_scan_attr(attr, "label", label_, IERR = ierr )
IF (ierr /= 0) RETURN 
!
CALL iotk_scan_attr(attr, "spin", spin_, IERR = ierr )
IF (ierr /= 0) RETURN
! 
CALL iotk_scan_attr(attr, "index", index_, IERR = ierr )
IF (ierr /= 0) RETURN
! 
CALL qes_init_hubbard_ns( obj, "Hubbard_ns", specie_, label_, spin_, index_, ndim_matrix_ns, ndim_matrix_ns, &
                          the_matrix ) 
! 
DEALLOCATE (the_matrix) 
! 
END SUBROUTINE qexsd_get_hubbard_ns
!-----------------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_dftU( iunit, obj, ispresent ) 
!-----------------------------------------------------------------------------------------------------------
! 
USE PARAMETERS,                      ONLY: ntypx
IMPLICIT NONE 
! 
INTEGER, INTENT (IN)                     :: iunit
TYPE (dftU_type),INTENT(OUT)             :: obj
LOGICAL,INTENT(OUT)                      :: ispresent
! 
INTEGER                                  :: ierr, iobj, l, lmax, lda_plus_u_kind_, &
                                            ndim_hubbard_alpha, ndim_hubbard_beta, ndim_hubbard_J, &
                                            ndim_hubbard_J0, ndim_hubbard_ns, ndim_hubbard_U, & 
                                            ndim_starting_ns 
CHARACTER(LEN=256)                       :: label_, specie_, u_projection_string
LOGICAL                                  :: lda_plus_u_kind_ispresent, hub_alpha_ispresent, &
                                            hub_beta_ispresent, hub_J0_ispresent, hub_Jvec_ispresent, &
                                            hub_U_ispresent, hub_ns_ispresent, starting_ns_ispresent, &
                                            u_projection_type_ispresent
TYPE ( hubbardCommon_type )             ::  hub_aux0, hub_U_obj(ntypx), hub_J0_obj(ntypx),& 
                                            hub_alpha_obj(ntypx),&
                                            hub_beta_obj(ntypx)
! 
TYPE (hubbardJ_type )                   ::  hub_Jvec_obj(ntypx), hubJ_aux0
!
TYPE (  starting_ns_type)               ::  ns_aux0, starting_ns_obj(ntypx)
!
TYPE   hubbard_ns_list 
   TYPE (hubbard_ns_type)              :: obj
   TYPE (hubbard_ns_list),POINTER      :: next
END TYPE hubbard_ns_list            
TYPE ( hubbard_ns_type )                :: hub_ns_aux0
TYPE ( hubbard_ns_type ),ALLOCATABLE    :: hub_ns_obj(:) 
TYPE ( hubbard_ns_list ),TARGET         :: hub_ns_list
TYPE ( hubbard_ns_list ),POINTER        :: hub_ns_last, hub_ns_first, hub_ns_ptr
CHARACTER(iotk_attlenx)                 :: attr
! 
ispresent = .FALSE. 
CALL iotk_scan_begin(iunit, "dftU", IERR = ierr, FOUND = ispresent) 
IF (ierr /=0) RETURN  
IF (.NOT. ispresent ) RETURN 
lda_plus_u_kind_ispresent = .FALSE. 
CALL iotk_scan_dat(iunit, "lda_plus_u_kind", lda_plus_u_kind_, FOUND = lda_plus_u_kind_ispresent ) 
!
ndim_hubbard_U=0
hub_U_ispresent= .FALSE.
CALL qexsd_get_hubbard_common(iunit, "Hubbard_U", hub_aux0, hub_U_ispresent)
IF (hub_U_ispresent) THEN  
   ndim_hubbard_U=1 
   count_hubbard_U:DO
      hub_U_obj(ndim_hubbard_U) = hub_aux0
      CALL qexsd_get_hubbard_common(iunit, "Hubbard_U", hub_aux0, hub_U_ispresent )
      IF (hub_U_obj(ndim_hubbard_U)%specie == hub_aux0%specie ) EXIT count_hubbard_U
      ndim_hubbard_U=ndim_hubbard_U+1
   END DO count_hubbard_U
END IF 
!
! 
ndim_hubbard_J0=0
hub_J0_ispresent= .FALSE.
CALL qexsd_get_hubbard_common(iunit, "Hubbard_J0", hub_aux0, hub_J0_ispresent)
IF (hub_J0_ispresent) THEN
   ndim_hubbard_J0=1
   count_hubbard_J0:DO
      hub_J0_obj(ndim_hubbard_J0) = hub_aux0
      CALL qexsd_get_hubbard_common(iunit, "Hubbard_J0", hub_aux0, hub_J0_ispresent )
      IF (hub_J0_obj(ndim_hubbard_J0)%specie == hub_aux0%specie) EXIT count_hubbard_J0
      ndim_hubbard_J0=ndim_hubbard_J0+1
   END DO count_hubbard_J0
END IF 
!
! 
ndim_hubbard_alpha=0
hub_alpha_ispresent= .FALSE.
CALL qexsd_get_hubbard_common(iunit, "Hubbard_alpha", hub_aux0,  hub_alpha_ispresent)
IF (hub_alpha_ispresent) THEN
   ndim_hubbard_alpha=1
   count_hubbard_alpha:DO
     hub_alpha_obj(ndim_hubbard_alpha) = hub_aux0
     CALL qexsd_get_hubbard_common(iunit, "Hubbard_alpha", hub_aux0,  hub_alpha_ispresent )
     IF (hub_alpha_obj(ndim_hubbard_alpha)%specie == hub_aux0%specie) EXIT count_hubbard_alpha
     ndim_hubbard_alpha=ndim_hubbard_alpha+1
   END DO count_hubbard_alpha
END IF
!
!
ndim_hubbard_beta=0
hub_beta_ispresent= .FALSE.
CALL qexsd_get_hubbard_common(iunit, "Hubbard_beta", hub_aux0, hub_beta_ispresent)
IF (hub_beta_ispresent) THEN
   ndim_hubbard_beta=1
   count_hubbard_beta:DO
      hub_beta_obj(ndim_hubbard_beta)  = hub_aux0
      CALL qexsd_get_hubbard_common(iunit, "Hubbard_beta", hub_aux0, hub_beta_ispresent )
      IF (hub_beta_obj(ndim_hubbard_beta)%specie == hub_aux0%specie)  EXIT count_hubbard_beta
      ndim_hubbard_beta=ndim_hubbard_beta+1
   END DO count_hubbard_beta
END IF
! 
ndim_hubbard_J=0
hub_Jvec_ispresent= .FALSE.
CALL qexsd_get_hubbard_J(iunit, "Hubbard_J", hubJ_aux0, hub_Jvec_ispresent)
IF (hub_Jvec_ispresent) THEN
   ndim_hubbard_J=1
   count_hubbard_J:DO
      hub_Jvec_obj(ndim_hubbard_J) = hubJ_aux0
      CALL qexsd_get_hubbard_J(iunit, "Hubbard_J", hubJ_aux0, hub_Jvec_ispresent )
      IF (hub_Jvec_obj(ndim_hubbard_J)%specie == hubJ_aux0%specie) EXIT count_hubbard_J
      ndim_hubbard_J=ndim_hubbard_J+1
   END DO count_hubbard_J
END IF
!
lmax=0
DO iobj =1, ndim_hubbard_U
   label_ = TRIM ( hub_U_obj(iobj)%label )
   SELECT CASE ( TRIM(label_))
      CASE ("1s", "2s", "3s", "4s", "5s", "6s", "7s") 
         l = 0
      CASE ("2p", "3p", "4p", "5p", "6p", "7p") 
         l = 1 
      CASE ("3d", "4d", "5d", "6d" )
         l = 2 
      CASE ("4f", "5f" )
         l = 3 
      CASE DEFAULT
         l = 0
   END SELECT
   lmax=max(l,lmax)
END DO
! 
ndim_starting_ns=0
starting_ns_ispresent= .FALSE.
CALL qexsd_get_starting_ns(iunit, ns_aux0, lmax, starting_ns_ispresent)
IF (starting_ns_ispresent) THEN
   ndim_starting_ns=1  
   count_starting_ns:DO
      starting_ns_obj(ndim_starting_ns) = ns_aux0
      CALL qexsd_get_starting_ns(iunit, ns_aux0, lmax, starting_ns_ispresent )
      IF (starting_ns_obj(ndim_starting_ns)%specie == ns_aux0%specie   .AND. & 
          starting_ns_obj(ndim_starting_ns)%label == ns_aux0%label   .AND.   &
          starting_ns_obj(ndim_starting_ns)%spin   == ns_aux0%spin   ) EXIT count_starting_ns
          ndim_starting_ns=ndim_starting_ns+1
   END DO count_starting_ns
END IF

ndim_Hubbard_ns=0
Hub_ns_ispresent= .FALSE.
CALL qexsd_get_Hubbard_ns(iunit, hub_ns_aux0, lmax, hub_ns_ispresent)
IF (hub_ns_ispresent) THEN 
   ndim_hubbard_ns=1
   hub_ns_last => hub_ns_list
   count_hubbard_ns:DO
     hub_ns_last%obj = hub_ns_aux0
     CALL qexsd_get_hubbard_ns(iunit, hub_ns_aux0, lmax, hub_ns_ispresent )
     IF ( hub_ns_last%obj%index == hub_ns_aux0%index .AND. &
  hub_ns_last%obj%spin  == hub_ns_aux0%spin ) EXIT count_hubbard_ns
     ! 
     ndim_hubbard_ns=ndim_hubbard_ns+1
     ALLOCATE (hub_ns_last%next) 
     hub_ns_last => hub_ns_last%next
     hub_ns_last%next =>null(hub_ns_last)
   END DO count_hubbard_ns
   ALLOCATE ( hub_ns_obj(ndim_hubbard_ns) ) 
   iobj = 1
   hub_ns_obj(1) = hub_ns_list%obj
   IF (associated( hub_ns_list%next ) ) THEN 
     hub_ns_first => hub_ns_list%next    
     copy_hubbard_ns:DO
        iobj = iobj +1 
        hub_ns_obj(iobj) = hub_ns_first%obj
        IF ( ASSOCIATED( hub_ns_first%next) ) THEN 
           hub_ns_ptr => hub_ns_first
           hub_ns_first => hub_ns_ptr%next 
           DEALLOCATE ( hub_ns_ptr)
        ELSE 
           DEALLOCATE ( hub_ns_first) 
           EXIT copy_hubbard_ns
        END IF 

     END DO copy_hubbard_ns
   END IF             
END IF           
      
!
CALL iotk_scan_dat(iunit, "U_projection_type", u_projection_string, FOUND = U_projection_type_ispresent, &
                   IERR = ierr ) 
! 
CALL iotk_scan_end( iunit, "dftU", IERR = ierr ) 
IF (ierr /=0) RETURN 
!
CALL qes_init_dftU(obj, "dftU", lda_plus_u_kind_ispresent, lda_plus_u_kind_, hub_U_ispresent, & 
                   ndim_hubbard_U, hub_U_obj, hub_J0_ispresent, ndim_hubbard_j0, hub_J0_obj, & 
                   hub_alpha_ispresent, ndim_hubbard_alpha, hub_alpha_obj,           &
                   hub_beta_ispresent, ndim_hubbard_beta, hub_beta_obj, hub_Jvec_ispresent,&   
                   ndim_hubbard_j, hub_Jvec_obj, starting_ns_ispresent, ndim_starting_ns, & 
                   starting_ns_obj, hub_ns_ispresent, ndim_hubbard_ns, hub_ns_obj, &
                   u_projection_type_ispresent, u_projection_string ) 
!
IF ( hub_U_ispresent ) THEN 
   DO iobj=1, ndim_hubbard_U
      CALL qes_reset_hubbardcommon(hub_U_obj(iobj))
   END DO
END IF 
! 
IF ( hub_J0_ispresent ) THEN
   DO iobj=1, ndim_hubbard_j0
      CALL qes_reset_hubbardcommon(hub_j0_obj(iobj))
   END DO
END IF
! 
IF ( hub_alpha_ispresent ) THEN
   DO iobj=1, ndim_hubbard_alpha
      CALL qes_reset_hubbardcommon(hub_alpha_obj(iobj))
   END DO
END IF
! 
IF ( hub_beta_ispresent ) THEN
   DO iobj=1, ndim_hubbard_beta
      CALL qes_reset_hubbardcommon(hub_beta_obj(iobj))
   END DO
END IF
! 
IF ( hub_Jvec_ispresent   ) THEN
   DO iobj=1, ndim_hubbard_j
      CALL qes_reset_hubbardJ(hub_Jvec_obj(iobj))
   END DO
END IF
! 
IF ( starting_ns_ispresent) THEN
   DO iobj=1, ndim_starting_ns
      CALL qes_reset_starting_ns(starting_ns_obj(iobj))
   END DO
END IF
! 
IF (hub_ns_ispresent )  THEN 
  DO iobj=1, ndim_hubbard_ns
     CALL qes_reset_hubbard_ns( hub_ns_obj(iobj)) 
  END DO 
  DEALLOCATE (hub_ns_obj)
END IF
! 
END SUBROUTINE qexsd_get_dftU
! 
!-------------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_vdW(iunit, obj, ispresent ) 
! 
USE parameters,            ONLY:  ntypx
IMPLICIT NONE 
! 
INTEGER, INTENT(IN)            :: iunit 
TYPE ( vdw_type ),INTENT(OUT)  :: obj
LOGICAL, INTENT(OUT)           :: ispresent
! 
INTEGER                        :: ierr 
CHARACTER(LEN=256)             :: vdw_corr_, non_locc_
REAL(DP)                       :: london_s6_, london_rcut_, xdm_a1_, xdm_a2_, ts_vdw_econv_thr_, london_rcut  
LOGICAL                        :: non_locterm_ispresent, lond_s6_ispresent, ts_vdw_econv_thr_ispresent,&
                                  ts_vdw_isolated_, ts_vdw_isolated_ispresent, london_s6_ispresent, &
                                  london_rcut_ispresent, xdm_a1_ispresent, xdm_a2_ispresent, london_c6_ispresent
INTEGER                        :: ndim_london_c6, iobj
TYPE ( hubbardCommon_type )   :: hub_aux0, london_c6_obj(ntypx)
CHARACTER(iotk_attlenx)         :: attr

! 
ispresent = .FALSE. 
CALL iotk_scan_begin(iunit, "vdW", ATTR = attr, IERR = ierr, FOUND = ispresent) 
IF (ierr /= 0 ) RETURN 
IF ( .NOT. ispresent ) RETURN 
! 
CALL iotk_scan_dat( iunit, "vdw_corr", vdw_corr_, IERR = ierr ) 
IF (ierr /=0 ) THEN 
   ispresent = .FALSE. 
   RETURN
END IF
CALL iotk_scan_dat ( iunit, "non_local_term", non_locc_, FOUND = non_locterm_ispresent, IERR = ierr ) 
IF ( ierr /= 0 ) THEN 
   ispresent = .FALSE. 
   RETURN
END IF 
CALL iotk_scan_dat( iunit, "london_s6", london_s6_, FOUND = lond_s6_ispresent, IERR = ierr )
IF (ierr /=0 ) THEN
   ispresent = .FALSE.
   RETURN
END IF
!
CALL iotk_scan_dat( iunit, "ts_vdw_econv_thr", ts_vdw_econv_thr_, FOUND = ts_vdw_econv_thr_ispresent, &
                    IERR = ierr ) 
IF (ierr /=0 ) THEN
   ispresent = .FALSE.
   RETURN
END IF
! 
CALL iotk_scan_dat( iunit, "ts_vdw_isolated", ts_vdw_isolated_, FOUND = ts_vdw_isolated_ispresent, &
                    IERR = ierr ) 
IF (ierr /=0 ) THEN
   ispresent = .FALSE.
   RETURN
END IF
! 
CALL iotk_scan_dat(iunit, "london_rcut", london_rcut_, FOUND = london_rcut_ispresent, IERR = ierr ) 
IF (ierr /=0 ) THEN
   ispresent = .FALSE.
   RETURN
END IF
!
CALL iotk_scan_dat(iunit, "xdm_a1", xdm_a1_, FOUND = xdm_a1_ispresent, IERR = ierr ) 
IF (ierr /= 0 ) THEN 
   ispresent = .FALSE. 
   RETURN
END IF
!
CALL iotk_scan_dat(iunit, "xdm_a2", xdm_a2_, FOUND = xdm_a2_ispresent, IERR = ierr )
IF (ierr /= 0 ) THEN
   ispresent = .FALSE.
   RETURN
END IF
!
ndim_london_c6=0
london_c6_ispresent= .FALSE.
CALL qexsd_get_hubbard_common(iunit, "london_c6", hub_aux0, london_c6_ispresent )
IF ( london_c6_ispresent ) THEN  
   ndim_london_c6=1 
   count_london_c6:DO
      london_c6_obj(ndim_london_c6) = hub_aux0
      CALL qexsd_get_hubbard_common(iunit, "london_c6", hub_aux0, london_c6_ispresent )
      IF (london_c6_obj(ndim_london_c6)%specie == hub_aux0%specie ) EXIT count_london_c6
      ndim_london_c6=ndim_london_c6 + 1
   END DO count_london_c6
END IF 
!

CALL iotk_scan_end (iunit, "vdW", IERR = ierr ) 
IF (ierr /= 0 ) RETURN 
!


CALL qes_init_vdw(obj, "vdW", vdw_corr_, non_locterm_ispresent, non_locc_, london_s6_ispresent, london_s6_,& 
                  ts_vdw_econv_thr_ispresent, ts_vdw_econv_thr_, ts_vdw_isolated_ispresent, ts_vdw_isolated_, &
                  london_rcut_ispresent, london_rcut, xdm_a1_ispresent, xdm_a1_, xdm_a2_ispresent, xdm_a2_, &
                  london_c6_ispresent, ndim_london_c6, london_c6_obj(1:ndim_london_c6) )
!  
IF ( london_c6_ispresent ) THEN 
   CALL qes_reset_hubbardCommon ( hub_aux0 )
   DO iobj = 1, ndim_london_c6
      CALL qes_reset_hubbardCommon ( london_c6_obj ( iobj ))
   END DO 
END IF
END SUBROUTINE qexsd_get_vdW
!---------------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_dft(iunit, obj, dft_ispresent) 
!
IMPLICIT NONE 
! 
INTEGER,INTENT(IN)                       :: iunit
TYPE( dft_type ),INTENT(OUT)             :: obj
LOGICAL,INTENT(OUT)                      :: dft_ispresent
! 
INTEGER                                  :: ierr, sum_err  
CHARACTER(LEN=256)                       :: functional_
LOGICAL                                  :: hybrid_ispresent, dftU_ispresent, vdW_ispresent
TYPE ( hybrid_type )                     :: hybrid_obj
TYPE ( dftU_type )                       :: dftU_obj
TYPE ( vdw_type )                         :: vdw_obj
CHARACTER(iotk_attlenx)                   :: attr
!  
dft_ispresent = .FALSE.
CALL iotk_scan_begin(iunit, "dft", ATTR = attr, IERR = ierr, FOUND = dft_ispresent)
IF (ierr /= 0 ) RETURN
CALL iotk_scan_dat( iunit, "functional", dat = functional_, ATTR = attr, IERR = ierr )
IF (ierr /=0) RETURN 
!
CALL qexsd_get_hybrid_dft(iunit, hybrid_obj, hybrid_ispresent) 
! 
CALL qexsd_get_dftU(iunit, dftU_obj, dftU_ispresent) 
! 
CALL qexsd_get_vdW( iunit, vdw_obj, vdW_ispresent ) 

! 
CALL iotk_scan_end ( iunit, "dft", IERR = ierr ) 
IF (ierr /= 0 ) RETURN 
!
CALL qes_init_dft ( obj, "dft", functional_, hybrid_ispresent, hybrid_obj, dftU_ispresent, dftU_obj, &
                    vdW_ispresent, vdw_obj ) 
! 
IF ( hybrid_ispresent )  CALL qes_reset_hybrid(hybrid_obj) 
IF ( dftU_ispresent )    CALL qes_reset_dftU( dftU_obj ) 
IF ( vdW_ispresent )     CALL qes_reset_vdw ( vdw_obj )  
!
END SUBROUTINE qexsd_get_dft   
! 
!---------------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_spin(iunit, obj, ispresent ) 
!---------------------------------------------------------------------------------------------------------
! 
IMPLICIT NONE 
! 
INTEGER, INTENT ( IN )                :: iunit
TYPE ( spin_type ), INTENT(OUT)       :: obj
LOGICAL, INTENT(OUT)                  :: ispresent 
! 
LOGICAL                               :: lsda_, noncolin_, spinorbit_
INTEGER                               :: ierr 
CHARACTER(iotk_attlenx)               :: attr
!
ispresent = .FALSE. 
CALL iotk_scan_begin(iunit, "spin", ATTR = attr, IERR = ierr , FOUND = ispresent)
IF ( ierr /= 0 ) RETURN 
!  
CALL iotk_scan_dat ( iunit, "lsda", lsda_, IERR = ierr )
IF (ierr /=0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "noncolin", noncolin_, IERR = ierr )
IF (ierr /= 0) RETURN 
! 
CALL iotk_scan_dat ( iunit, "spinorbit", spinorbit_, IERR = ierr )
IF (ierr /= 0) RETURN 
! 
CALL iotk_scan_end( iunit, "spin" , IERR = ierr ) 
IF (ierr /= 0 ) RETURN 
! 
CALL qes_init_spin( obj, "spin", lsda_, noncolin_, spinorbit_ ) 
!  
END SUBROUTINE qexsd_get_spin
! 
!------------------------------------------------------------------------------------
SUBROUTINE  qexsd_get_smearing( iunit, obj, ispresent ) 
!-------------------------------------------------------------------------------------
! 
IMPLICIT NONE 
! 
INTEGER,INTENT(IN)                   :: iunit
TYPE ( smearing_type ),INTENT(OUT)   :: obj
LOGICAL,INTENT(OUT)                  :: ispresent 
!
INTEGER                              :: ierr  
REAL(DP)                             :: degauss_
CHARACTER(256)                       :: smearing_choice_
CHARACTER(iotk_attlenx)              :: attr
! 
ispresent = .FALSE. 
CALL iotk_scan_dat( iunit, "smearing", smearing_choice_, ATTR = attr, IERR = ierr , FOUND = ispresent )
IF ( ierr /= 0 ) CALL errore ( "qexsd_get_smearing", "error reading smearing element", 1)
IF ( ispresent ) THEN   
   CALL iotk_scan_attr ( attr, "degauss", degauss_, IERR = ierr ) 
   IF (ierr /= 0 ) CALL errore ( "qexsd_get_smearing", "error reading degauss value in smearing element", 1) 
   CALL qes_init_smearing( obj, "smearing", degauss_, TRIM( smearing_choice_ ) )
END IF 
END SUBROUTINE qexsd_get_smearing
!-------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_bands( iunit, obj, ispresent  ) 
!-------------------------------------------------------------------------------------------------
!
! 
IMPLICIT NONE 
! 
INTEGER, INTENT ( IN )                    :: iunit
TYPE ( bands_type ), INTENT(OUT)          :: obj
LOGICAL, INTENT(OUT)                      :: ispresent 
! 
INTEGER                                   :: ierr
INTEGER                                   :: nbnd_, occ_spin_, ndim_inputOcc, &
                                             ispin_,  ndim_occ_vec
REAL(DP)                                  :: spin_factor_
LOGICAL                                   :: nbnd_ispresent, smearing_ispresent, tot_chg_ispresent, &
                                             inputOcc_is, spin_occ_is, tot_magnetization_ispresent, found
TYPE (smearing_type)                      :: smearing_obj
TYPE ( occupations_type)                  :: occ_obj
TYPE ( inputOccupations_type)             :: inpOcc_obj(2) 
REAL(DP)                                  :: tot_chg_, tot_mag_
REAL(DP), ALLOCATABLE                     :: input_occ_vec(:)
CHARACTER(LEN=256)                        :: occ_string_
INTEGER                                   :: iobj 
CHARACTER(iotk_attlenx)                   :: attr
 
ispresent = .FALSE. 
CALL iotk_scan_begin( iunit, "bands", ATTR = attr, IERR = ierr , FOUND = ispresent ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat( iunit, "nbnd", nbnd_, IERR = ierr , FOUND = nbnd_ispresent )  
IF ( ierr /=0 ) RETURN  
CALL qexsd_get_smearing ( iunit, smearing_obj, smearing_ispresent ) 
! 
CALL iotk_scan_dat ( iunit, "tot_charge", tot_chg_, IERR = ierr, FOUND = tot_chg_ispresent )
IF ( ierr /= 0) RETURN 
! 
CALL iotk_scan_dat ( iunit, "tot_magnetization", tot_mag_, IERR = ierr, FOUND = tot_magnetization_ispresent) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "occupations",  occ_string_, ATTR = attr, IERR = ierr, FOUND = inputOcc_is )
IF ( ierr /= 0 ) RETURN 
!
CALL iotk_scan_attr (attr, "spin", occ_spin_, IERR = ierr , FOUND = spin_occ_is) 
IF  ( ierr /= 0 ) RETURN
!
CALL qes_init_occupations( occ_obj, "occupations", occ_spin_, spin_occ_is, occ_string_) 
!
inputOcc_is = .FALSE.   
IF ( nbnd_ispresent ) THEN
   ndim_occ_vec = nbnd_
   ALLOCATE (input_occ_vec(ndim_occ_vec)) 
   CALL iotk_scan_dat(iunit, "inputOccupations", input_occ_vec, ATTR = attr, IERR = ierr, &
                      FOUND = inputOcc_is )   
   IF ( ierr /= 0 ) RETURN 
   IF (inputOcc_is ) THEN 
      CALL iotk_scan_attr(attr, "ispin",      ispin_, IERR = ierr ) 
      CALL iotk_scan_attr(attr, "spin_factor",spin_factor_, IERR = ierr )
      CALL qes_init_inputOccupations ( inpOcc_obj(1), "inputOccupations", ispin_, spin_factor_, &
                                    ndim_occ_vec, input_occ_vec) 
   ! 
      CALL iotk_scan_dat(iunit, "inputOccupations", input_occ_vec, ATTR = attr, IERR = ierr ,   &
                         FOUND = found )
      IF ( ierr /= 0 ) RETURN 
      CALL iotk_scan_attr(attr, "ispin",      ispin_, IERR = ierr )
      CALL iotk_scan_attr(attr, "spin_factor",spin_factor_, IERR = ierr ) 
      CALL qes_init_inputOccupations ( inpOcc_obj(2), "inputOccupations", ispin_, spin_factor_, &
                                    ndim_occ_vec, input_occ_vec)
   ! 
      IF (inpOcc_obj(1)%ispin == inpOcc_obj(2)%ispin ) THEN 
         ndim_inputOcc = 1 
         CALL qes_reset_inputOccupations(inpOcc_obj(2))
      ELSE 
         ndim_inputOcc = 2 
      END IF   
   END IF
   DEALLOCATE ( input_occ_vec ) 
END IF          
! 
      
CALL iotk_scan_end ( iunit, "bands", IERR = ierr ) 
IF (ierr /=0) RETURN 
!
CALL qes_init_bands(obj, "bands", nbnd_ispresent, nbnd_, smearing_ispresent, smearing_obj, &  
     tot_chg_ispresent, tot_chg_, tot_magnetization_ispresent, tot_mag_, occ_obj, inputOcc_is, &
     ndim_inputOcc, inpOcc_obj( 1:ndim_inputOcc ) )
!  
CALL qes_reset_smearing ( smearing_obj ) 
CALL qes_reset_occupations ( occ_obj ) 
IF ( inputOcc_is ) THEN 
   DO iobj = 1, ndim_inputOcc 
      CALL qes_reset_inputOccupations ( inpOcc_obj ( iobj ))
   END DO 
END IF
END SUBROUTINE qexsd_get_bands
! 
!-------------------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_basis(iunit, obj, ispresent ) 
!-------------------------------------------------------------------------------------------------------------
! 
IMPLICIT NONE 
! 
! 
INTEGER, INTENT ( IN )                    :: iunit
TYPE ( basis_type ), INTENT(OUT)          :: obj
LOGICAL, INTENT(OUT)                      :: ispresent
! 
INTEGER                                   :: ierr 
LOGICAL                                   :: gamma_only_, gamma_only_ispresent, ecutrho_ispresent,&
                                             ngms_ispresent, box_ispresent, smooth_ispresent, grid_ispresent
CHARACTER(LEN=256)                        :: empty_str
REAL(DP),DIMENSION(3)                     :: b1_, b2_, b3_
REAL(DP)                                  :: ecutrho_, ecutwfc_  
INTEGER                                   :: fft_nr1_,    fft_nr2_,    fft_nr3_,&
                                             smooth_nr1_, smooth_nr2_, smooth_nr3_,&
                                             box_nr1_,    box_nr2_,    box_nr3_,&
                                             ngm_, ngms_, npwx_
TYPE ( basisSetItem_type)                :: grid_obj, smooth_obj, box_obj
CHARACTER(iotk_attlenx)                  :: attr
! 
ispresent = .FALSE. 
CALL iotk_scan_begin ( iunit, "basis", ATTR = attr, IERR = ierr, FOUND = ispresent ) 
IF (ierr /= 0) RETURN   
! 
gamma_only_ispresent = .FALSE. 
CALL iotk_scan_dat( iunit, "gamma_only", gamma_only_, FOUND = gamma_only_ispresent, IERR = ierr )
! 
CALL iotk_scan_dat ( iunit, "ecutwfc", ecutwfc_, IERR = ierr ) 
IF (ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "ecutrho", ecutrho_, IERR = ierr , FOUND = ecutrho_ispresent )
! 
CALL iotk_scan_dat (iunit, "fft_grid", dat = empty_str, ATTR = attr, IERR = ierr , FOUND = grid_ispresent )  
IF ( grid_ispresent ) THEN 
   CALL iotk_scan_attr(attr, "nr1", fft_nr1_, IERR = ierr ) 
   IF ( ierr /= 0) RETURN 
   CALL iotk_scan_attr( attr, "nr2", fft_nr2_, IERR = ierr ) 
   IF ( ierr /= 0) RETURN
   CALL iotk_scan_attr( attr, "nr3", fft_nr3_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL qes_init_basisSetItem ( grid_obj, "fft_grid", fft_nr1_, fft_nr2_, fft_nr3_ ,"")  
END IF
! 
CALL iotk_scan_dat (iunit, "fft_smooth", dat = empty_str, ATTR = attr, &
    FOUND = smooth_ispresent, IERR = ierr )
IF ( smooth_ispresent ) THEN  
   CALL iotk_scan_attr(attr, "nr1", smooth_nr1_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL iotk_scan_attr( attr, "nr2", smooth_nr2_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL iotk_scan_attr( attr, "nr3", smooth_nr3_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL qes_init_basisSetItem ( smooth_obj, "smooth_grid", smooth_nr1_, smooth_nr2_, smooth_nr3_ ,"")
END IF
!
CALL iotk_scan_dat (iunit, "fft_box", dat = empty_str, ATTR = attr, &
    FOUND = box_ispresent, IERR = ierr )
IF ( box_ispresent ) THEN  
   CALL iotk_scan_attr(attr, "nr1", box_nr1_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL iotk_scan_attr( attr, "nr2", box_nr2_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL iotk_scan_attr( attr, "nr3", box_nr3_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL qes_init_basisSetItem ( box_obj, "fft_box", box_nr1_, box_nr2_, box_nr3_ ,"")
END IF
CALL iotk_scan_end( iunit, "basis", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
CALL qes_init_basis ( obj, "basis", gamma_only_ispresent, gamma_only_, ecutwfc_, ecutrho_ispresent, &
      ecutrho_, grid_ispresent, grid_obj, smooth_ispresent, smooth_obj, box_ispresent, &
      box_obj )  
!
IF ( grid_ispresent) CALL qes_reset_basisSetItem(grid_obj)
IF ( smooth_ispresent) CALL qes_reset_basisSetItem(smooth_obj)
IF ( box_ispresent) CALL qes_reset_basisSetItem(box_obj)
END SUBROUTINE qexsd_get_basis 
! 
!-------------------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_electron_control( iunit, obj, ispresent ) 
!-------------------------------------------------------------------------------------------------------
! 
IMPLICIT NONE 
! 
INTEGER,INTENT (IN)                         :: iunit
TYPE (electron_control_type),INTENT(OUT)    :: obj
LOGICAL, INTENT(OUT)                        :: ispresent 
! 
INTEGER                                     :: ierr
CHARACTER(LEN=256)                          :: diago_str_, mixing_str_ 
REAL(DP)                                    :: mixing_beta_, conv_thr_, diago_thr_init_
INTEGER                                     :: mixing_ndim_, max_nstep_, diago_cg_maxiter_
LOGICAL                                     :: real_space_q_, diago_full_acc_, found_, &
                                               tq_smoothing_ = .FALSE., tbeta_smoothing_ =.FALSE.
CHARACTER(iotk_attlenx)                     :: attr
! 
ispresent = .FALSE.
CALL iotk_scan_begin( iunit, "electron_control", FOUND = ispresent, IERR = ierr )
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat( iunit, "diagonalization", diago_str_, IERR = ierr) 
IF ( ierr/= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "mixing_mode", mixing_str_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL  iotk_scan_dat( iunit, "mixing_beta", mixing_beta_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!  
CALL  iotk_scan_dat( iunit, "conv_thr", conv_thr_, IERR = ierr )
IF ( ierr /= 0 ) RETURN
!
CALL  iotk_scan_dat( iunit, "mixing_ndim" , mixing_ndim_, IERR = ierr )
IF ( ierr /= 0 ) RETURN
! 
CALL  iotk_scan_dat( iunit, "max_nstep" , max_nstep_, IERR = ierr )
IF ( ierr /= 0 ) RETURN
!
CALL  iotk_scan_dat( iunit, "real_space_q" , real_space_q_, IERR = ierr )
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "tq_smoothing", tq_smoothing_, IERR = ierr , FOUND = found_) 
IF ( ierr /= 0 ) RETURN 
!
CALL iotk_scan_dat ( iunit, "tbeta_smoothing", tbeta_smoothing_, IERR = ierr , FOUND = found_)
IF ( ierr /= 0 ) RETURN
! 
CALL  iotk_scan_dat( iunit, "diago_thr_init" , diago_thr_init_, IERR = ierr )
IF ( ierr /= 0 ) RETURN
!
CALL  iotk_scan_dat( iunit, "diago_full_acc" , diago_full_acc_, IERR = ierr )
IF ( ierr /= 0 ) RETURN
!
CALL  iotk_scan_dat( iunit, "diago_cg_maxiter" , diago_cg_maxiter_, IERR = ierr )
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_end( iunit, "electron_control", IERR = ierr ) 
IF  ( ierr /= 0 ) THEN 
   ispresent = .FALSE.
   RETURN 
END IF
CALL qes_init_electron_control ( obj, "electron_control", diago_str_, mixing_str_, mixing_beta_, & 
                                 conv_thr_, mixing_ndim_, max_nstep_, real_space_q_, &
                                 tq_smoothing_, tbeta_smoothing_, diago_thr_init_, &
                                 diago_full_acc_, diago_cg_maxiter_) 
! 
END SUBROUTINE qexsd_get_electron_control
! 
!-----------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_k_points_IBZ( iunit, obj, ispresent , tagname ) 
!-----------------------------------------------------------------------------------------
! 
IMPLICIT NONE 
! 
INTEGER, INTENT (IN)                   :: iunit
TYPE (k_points_IBZ_type),INTENT(OUT)   :: obj
LOGICAL, INTENT (OUT)                  :: ispresent
CHARACTER(LEN=*),OPTIONAL,INTENT(IN)   :: tagname
! 
INTEGER                                :: ierr, nks_, ik, & 
                                          nk1_, nk2_, nk3_, k1_, k2_, k3_  
LOGICAL                                :: monkh_pack_ispresent, nks_ispresent, weight_ispresent, &
                                          label_ispresent 
CHARACTER(LEN=256)                     :: label, empty_str, tagname_
TYPE (k_point_type),ALLOCATABLE        :: kp_obj(:)
TYPE ( monkhorst_pack_type )           :: mp_obj
REAL (DP)                              :: wk_, xk_(3)
CHARACTER(iotk_attlenx)                :: attr
! 
ispresent = .FALSE.
IF (PRESENT(tagname)) THEN   
   tagname_ = TRIM(tagname(1:MIN(LEN(tagname),256)))
ELSE 
   tagname_ = "k_points_IBS"
END IF
CALL iotk_scan_begin( iunit, TRIM(tagname_) , FOUND = ispresent, IERR = ierr ) 
IF ( ierr /= 0 ) CALL errore ( "qexsd_get_k_points_IBZ", "error reading element from xml file", ierr) 
IF (.NOT. ispresent ) CALL errore ( "qexsd_get_k_points_IBZ", TRIM(tagname_)//" not found" , ierr ) 
! 
CALL iotk_scan_dat( iunit, "monkhorst_pack", empty_str, ATTR = attr, FOUND = monkh_pack_ispresent, IERR = ierr)
IF ( monkh_pack_ispresent ) THEN 
   CALL iotk_scan_attr( attr, "nk1", nk1_, IERR = ierr ) 
   IF (ierr /= 0 ) RETURN    
   CALL iotk_scan_attr( attr, "nk2", nk2_, IERR = ierr )
   IF (ierr /= 0 ) RETURN
   CALL iotk_scan_attr( attr,  "nk3", nk3_, IERR = ierr )
   IF (ierr /= 0 ) RETURN
   ! 
   CALL iotk_scan_attr( attr, "k1", k1_, IERR = ierr )
   IF (ierr /= 0 ) RETURN
   CALL iotk_scan_attr( attr, "k2", k2_, IERR = ierr )
   IF (ierr /= 0 ) RETURN
   CALL iotk_scan_attr( attr, "k3", k3_, IERR = ierr )
   IF (ierr /= 0 ) RETURN
   ! 
   CALL qes_init_monkhorst_pack(mp_obj, "monkhorst_pack", nk1_, nk2_, nk3_, k1_, k2_, k3_, "") 
END IF 
CALL iotk_scan_dat(iunit, "nk", nks_, FOUND = nks_ispresent , IERR = ierr ) 
IF ( nks_ispresent .AND. monkh_pack_ispresent ) THEN
   ispresent = .FALSE. 
   RETURN 
END IF 
IF ( .NOT. ( nks_ispresent .OR. monkh_pack_ispresent ) ) THEN 
   ispresent = .FALSE. 
   RETURN
END IF 
IF ( nks_ispresent) THEN 
   ALLOCATE ( kp_obj(nks_)) 
   DO ik= 1, nks_
      CALL iotk_scan_dat(iunit, "k_point", xk_, ATTR = attr, IERR = ierr ) 
      IF ( ierr /= 0 ) RETURN 
      wk_ = 0.d0
      CALL iotk_scan_attr(attr, "weight", wk_, FOUND = weight_ispresent ) 
      CALL iotk_scan_attr(attr, "label", label, FOUND = label_ispresent )
      IF (.NOT. label_ispresent ) label=" "
      CALL qes_init_k_point (kp_obj(ik),"k_point", wk_, weight_ispresent, TRIM(label), label_ispresent, xk_) 
   END DO
END IF    
!  
CALL iotk_scan_end(iunit, TRIM(tagname_) , IERR = ierr ) 
IF (ierr /= 0 ) CALL errore ( "qexsd_get_k_points_IBZ", "k_points_IBZ element is not correctly closed",1) 
CALL qes_init_k_points_IBZ( obj, "k_points_IBZ" , monkh_pack_ispresent, mp_obj, nks_ispresent, & 
                            nks_, nks_ispresent, nks_, kp_obj )
!
IF ( monkh_pack_ispresent ) CALL qes_reset_monkhorst_pack(mp_obj) 
IF ( nks_ispresent ) THEN 
   DO ik = 1, nks_
      CALL qes_reset_k_point(kp_obj(ik))
   END DO
   DEALLOCATE (kp_obj)
END IF 
! 
END SUBROUTINE qexsd_get_k_points_IBZ 
!
!-----------------------------------------------------------------------------
SUBROUTINE qexsd_get_bfgs( iunit, obj, ispresent ) 
! 
INTEGER, INTENT(IN)                   :: iunit
TYPE (bfgs_type),INTENT(OUT)          :: obj
LOGICAL,INTENT(OUT)                   :: ispresent
! 
INTEGER                               :: ierr, ndim_ 
REAL(DP)                              :: trust_radius_min_, trust_radius_max_, trust_radius_init_, &
                                         w1_, w2_
CHARACTER(iotk_attlenx)               :: attr
! 
ispresent = .FALSE. 
CALL iotk_scan_begin( iunit, "bfgs", IERR = ierr , FOUND = ispresent ) 
IF ( ierr /= 0 ) RETURN
IF ( .NOT. ispresent ) RETURN 
! 

CALL iotk_scan_dat( iunit, "ndim", ndim_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat( iunit, "trust_radius_min", trust_radius_min_, IERR = ierr ) 
IF (ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat(iunit, "trust_radius_max", trust_radius_max_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat(iunit, "trust_radius_init", trust_radius_init_, IERR =ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat(iunit, "w1", w1_, IERR =ierr )
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat(iunit, "w2", w2_, IERR =ierr )
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_end ( iunit, "bfgs", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 

CALL qes_init_bfgs( obj, "bfgs", ndim_, trust_radius_min_, trust_radius_max_, trust_radius_init_, & 
    w1_, w2_ ) 
END SUBROUTINE qexsd_get_bfgs
!
!----------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_md( iunit, obj, ispresent )
! 
IMPLICIT NONE 
! 
INTEGER, INTENT(IN)                   :: iunit
TYPE (md_type),INTENT(OUT)            :: obj
LOGICAL,INTENT(OUT)                   :: ispresent
! 
INTEGER                               :: ierr, nraise_
CHARACTER(LEN =256)                   :: pot_extrapolation_, wfc_extrapolation_, ion_temperature_
REAL(DP)                              :: timestep_, tempw_, tolp_, deltaT_ 
CHARACTER(iotk_attlenx)               :: attr
! 
ispresent = .FALSE.
CALL iotk_scan_begin(iunit, "md", IERR = ierr, FOUND = ispresent)
IF ( ierr /= 0 ) RETURN 
IF (.NOT. ispresent ) RETURN 
! 
CALL iotk_scan_dat( iunit, "pot_extrapolation", pot_extrapolation_, IERR = ierr) 
IF (ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat(iunit, "wfc_extrapolation", wfc_extrapolation_, IERR = ierr ) 
IF (ierr /= 0 ) RETURN 
! 
 CALL iotk_scan_dat(iunit, "ion_temperature", ion_temperature_, IERR = ierr )
IF (ierr /= 0 ) RETURN
!  
CALL iotk_scan_dat(iunit, "timestep", timestep_, DEFAULT = 20.d0, IERR = ierr )
IF (ierr /= 0 ) RETURN
!
CALL iotk_scan_dat(iunit, "timestep", tempw_, IERR = ierr )
IF (ierr /= 0 ) RETURN
!
CALL iotk_scan_dat(iunit, "tolp", tolp_, IERR = ierr )
IF (ierr /= 0 ) RETURN
!
CALL iotk_scan_dat(iunit, "deltaT", deltaT_, IERR = ierr )
IF (ierr /= 0 ) RETURN
!
CALL iotk_scan_dat(iunit, "nraise", nraise_, IERR = ierr )
IF (ierr /= 0 ) RETURN
!
CALL iotk_scan_end(iunit, "md", IERR = ierr )
IF ( ierr /= 0 ) RETURN 
CALL qes_init_md(obj, "md", pot_extrapolation_, wfc_extrapolation_, ion_temperature_, timestep_, &
                 tempw_, tolp_, deltaT_, nraise_)   
!
END SUBROUTINE qexsd_get_md
!----------------------------------------------------------------------------------
SUBROUTINE  qexsd_get_ion_control(iunit, obj, ispresent) 
!----------------------------------------------------------------------------------
! 
INTEGER, INTENT(IN)                   :: iunit
TYPE (ion_control_type),INTENT(OUT)   :: obj 
LOGICAL,INTENT(OUT)                   :: ispresent
! 
INTEGER                               :: ierr 
CHARACTER(LEN=256)                    :: ion_dynamics_
REAL(DP)                              :: upscale_
LOGICAL                               :: upscale_ispresent, refold_pos_ispresent, refold_pos_, &
                                         bfgs_ispresent, md_ispresent, remove_rig_rot_is, remove_rig_rot_ 
TYPE ( bfgs_type )                    :: bfgs_obj
TYPE ( md_type )                      :: md_obj 
CHARACTER(iotk_attlenx)               :: attr
! 
ispresent = .FALSE. 
CALL iotk_scan_begin(iunit, "ion_control", IERR= ierr, FOUND = ispresent)
IF (  .NOT. ispresent ) RETURN 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat(iunit, "ion_dynamics", ion_dynamics_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat(iunit, "upscale", upscale_, IERR = ierr , FOUND = upscale_ispresent)
! 
CALL iotk_scan_dat(iunit, "remove_rigid_rot", remove_rig_rot_, FOUND = remove_rig_rot_is)
CALL iotk_scan_dat(iunit, "refold_pos", refold_pos_, IERR = ierr, FOUND = refold_pos_ispresent )     
! 
CALL qexsd_get_bfgs( iunit, bfgs_obj, bfgs_ispresent )
! 
CALL qexsd_get_md( iunit, md_obj, md_ispresent ) 
!  
CALL iotk_scan_end ( iunit, "ion_control", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
!
CALL qes_init_ion_control( obj, "ion_control", ion_dynamics_, upscale_ispresent, upscale_, &
                           remove_rig_rot_is, remove_rig_rot_, refold_pos_ispresent, refold_pos_, &
                           bfgs_ispresent, bfgs_obj, md_ispresent, md_obj)
END SUBROUTINE qexsd_get_ion_control                                
!
!------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_cell_control( iunit, obj, ispresent) 
! 
IMPLICIT NONE 
! 
INTEGER, INTENT(IN)                   :: iunit
TYPE (cell_control_type),INTENT(OUT)  :: obj
LOGICAL,INTENT(OUT)                   :: ispresent
!                                
INTEGER                               :: ierr 
CHARACTER ( LEN = 256 )               :: cell_dynamics_
REAL(DP)                              :: wmass_, cell_factor_, pressure_
LOGICAL                               :: cell_fac_ispresent, wmass_ispresent, fix_vol_ispresent, & 
                                         fix_area_ispresent, free_cell_ispresent, isotropic_ispresent
LOGICAL                               :: fix_volume_, fix_area_, isotropic_, found_pressure
INTEGER                               :: free_cell_mat_(3,3)
TYPE (integerMatrix_type)             :: mat_obj
CHARACTER(iotk_attlenx)               :: attr
 
!
ispresent = .FALSE. 
CALL iotk_scan_begin( iunit, "cell_control", IERR = ierr , FOUND = ispresent) 
IF ( ierr /= 0 ) RETURN 
IF ( .NOT. ispresent ) RETURN 
! 
CALL iotk_scan_dat( iunit, "cell_dynamics", cell_dynamics_, IERR = ierr ) 
IF (ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat( iunit, "wmass", wmass_, IERR = ierr , FOUND = wmass_ispresent)
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat( iunit, "pressure", pressure_, IERR = ierr, FOUND = found_pressure) 
IF ( ierr /= 0 ) RETURN 
IF (.NOT. found_pressure) pressure_=0.d0
!
CALL iotk_scan_dat( iunit, "cell_factor", cell_factor_, IERR = ierr , FOUND = cell_fac_ispresent)
IF ( ierr /= 0 ) RETURN
!  
CALL iotk_scan_dat( iunit, "fix_volume", fix_volume_, IERR = ierr , FOUND = fix_vol_ispresent)
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat( iunit, "fix_area", fix_area_, IERR = ierr, FOUND = fix_area_ispresent )
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat( iunit, "isotropic", isotropic_, IERR = ierr, FOUND = isotropic_ispresent )
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat ( iunit, "free_cell", free_cell_mat_, IERR = ierr, FOUND = free_cell_ispresent ) 
IF ( ierr /= 0 ) RETURN 
CALL qes_init_integerMatrix(mat_obj, "free_cell", 3, 3, free_cell_mat_)
! 
CALL iotk_scan_end ( iunit, "cell_control", IERR = ierr )
IF ( ierr /= 0)  RETURN 
! 
CALL qes_init_cell_control( obj, "cell_control", cell_dynamics_, pressure_, wmass_ispresent, wmass_, & 
                            cell_fac_ispresent, cell_factor_, fix_vol_ispresent, fix_volume_, &
                            fix_area_ispresent, fix_area_, isotropic_ispresent, isotropic_, & 
                            free_cell_ispresent, mat_obj)
END SUBROUTINE qexsd_get_cell_control 
! 
!-------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_symmetry_flags( iunit, obj, ispresent )
! 
IMPLICIT NONE
!
INTEGER, INTENT(IN)                    :: iunit
TYPE (symmetry_flags_type),INTENT(OUT) :: obj
LOGICAL,INTENT(OUT)                    :: ispresent
! 
INTEGER                                :: ierr 
LOGICAL                                :: nosym_= .FALSE.,&
                                          nosym_evc_ = .FALSE., &
                                          noinv_ = .FALSE.,&
                                          no_t_rev_ = .FALSE.,&
                                          force_symmorphic_ = .FALSE., &
                                          use_all_frac_ = .FALSE. , found_
CHARACTER(iotk_attlenx)                :: attr
! 
ispresent = .FALSE.
CALL iotk_scan_begin( iunit, "symmetry_flags", IERR = ierr, FOUND = ispresent ) 
IF ( ierr /= 0) RETURN 
IF ( .NOT. ispresent ) RETURN 
! 
CALL iotk_scan_dat( iunit, "nosym", nosym_, IERR = ierr , FOUND = found_) 
IF (ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat( iunit, "nosym_evc", nosym_evc_, IERR = ierr , FOUND = found_)
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat( iunit, "noinv", noinv_, IERR = ierr , FOUND = found_)
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat( iunit, "no_t_rev", no_t_rev_, IERR = ierr , FOUND = found_)
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat( iunit, "force_symmorphic", force_symmorphic_, IERR = ierr , FOUND = found_)
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat( iunit, "use_all_frac", use_all_frac_, IERR = ierr , FOUND = found_)
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_end( iunit, "symmetry_flags", IERR = ierr)
IF ( ierr /= 0 ) RETURN 
! 
CALL qes_init_symmetry_flags( obj, "symmetry_flags", nosym_, nosym_evc_, noinv_, no_t_rev_, & 
                              force_symmorphic_, use_all_frac_ ) 
END SUBROUTINE qexsd_get_symmetry_flags  
!
!--------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_esm(iunit, obj, ispresent )
! 
IMPLICIT NONE
! 
INTEGER, INTENT(IN)               :: iunit
TYPE (esm_type),INTENT(OUT)       :: obj
LOGICAL,INTENT(OUT)               :: ispresent
!
INTEGER                           :: ierr, nfit_
CHARACTER ( LEN =256 )            :: bc_
REAL(DP)                          :: efield_, w_
CHARACTER(iotk_attlenx)           :: attr
!
ispresent = .FALSE. 
CALL iotk_scan_begin(iunit, "esm", IERR = ierr, FOUND = ispresent)
IF ( ierr /= 0 ) RETURN 
IF ( .NOT.  ispresent ) RETURN 
! 
CALL iotk_scan_dat(iunit, "bc", bc_, IERR = ierr )
IF ( ierr /= 0) RETURN 
!
CALL iotk_scan_dat ( iunit, "nfit", nfit_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "w", w_, IERR = ierr )
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat ( iunit, "efield", efield_, IERR = ierr )
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_end( iunit, "esm", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL qes_init_esm( obj, "esm", bc_, nfit_, w_, efield_)
END SUBROUTINE qexsd_get_esm
! 
!---------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_boundary_conditions(iunit, obj, ispresent ) 
!---------------------------------------------------------------------------------------------------
! 
IMPLICIT NONE 
! 
INTEGER, INTENT(IN)                         :: iunit
TYPE (boundary_conditions_type),INTENT(OUT) :: obj
LOGICAL,INTENT(OUT)                         :: ispresent
! 
INTEGER                                     :: ierr 
CHARACTER( LEN = 256 )                      :: assume_isolated_
LOGICAL                                     :: fcp_opt_ = .FALSE. 
REAL(DP)                                    :: fcp_mu_ 
TYPE ( esm_type )                           :: esm_obj
LOGICAL                                     :: esm_ispresent, fcp_opt_ispresent, fcp_mu_ispresent
CHARACTER(iotk_attlenx)                     :: attr
! 
ispresent = .FALSE. 
CALL iotk_scan_begin( iunit, "boundary_conditions", IERR = ierr, FOUND = ispresent) 
IF ( ierr /= 0 ) RETURN
IF (.NOT. ispresent ) RETURN 
! 
CALL iotk_scan_dat( iunit, "assume_isolated", assume_isolated_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL qexsd_get_esm( iunit, esm_obj, esm_ispresent )  
! 
CALL iotk_scan_dat(iunit, "fcp_opt", fcp_opt_, FOUND = fcp_opt_ispresent, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
!
CALL iotk_scan_dat(iunit, "fcp_mu", fcp_mu_, FOUND = fcp_mu_ispresent, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
IF (  fcp_opt_ .AND. ( .NOT. fcp_mu_ispresent ) ) &
   CALL errore( "qexsd_get_boundary_conditions","found fcp_opt true but no value for fcp_mu has been found", 10) 
CALL iotk_scan_end(iunit, "boundary_conditions", IERR = ierr)
IF ( ierr /= 0 ) RETURN
!
CALL qes_init_boundary_conditions(obj, "boundary_conditions", assume_isolated_, esm_ispresent, esm_obj, &
                                  fcp_opt_ispresent, fcp_opt_, fcp_mu_ispresent, fcp_mu_ ) 
!
END SUBROUTINE qexsd_get_boundary_conditions
!
!-----------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_ekin_functional( iunit, obj, ispresent) 
! 
IMPLICIT NONE
! 
INTEGER, INTENT(IN)                         :: iunit
TYPE (ekin_functional_type),INTENT(OUT) :: obj
LOGICAL,INTENT(OUT)                         :: ispresent
! 
INTEGER                                     :: ierr
REAL(DP)                                    :: ecfixed_, qcutz_, q2sigma_
CHARACTER(iotk_attlenx)                     :: attr

!
ispresent = .FALSE. 
CALL iotk_scan_begin ( iunit, "ekin_functional", FOUND = ispresent, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
IF ( .NOT. ispresent ) RETURN 
! 
CALL iotk_scan_dat( iunit, "ecfixed", ecfixed_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
!  
CALL iotk_scan_dat( iunit, "qcutz",  qcutz_, IERR = ierr )
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat( iunit, "q2sigma", q2sigma_, IERR = ierr )
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_end( iunit, "ekin_functional", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
! 
CALL qes_init_ekin_functional( obj, "ekin_functional", ecfixed_, qcutz_, q2sigma_)
END SUBROUTINE qexsd_get_ekin_functional
!
!-----------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_electric_field( iunit, obj, ispresent )
! 
IMPLICIT NONE
!  
INTEGER, INTENT(IN)                      :: iunit
TYPE (electric_field_type),INTENT(OUT)   :: obj
LOGICAL,INTENT(OUT)                      :: ispresent
!
INTEGER                                  :: ierr 
CHARACTER(LEN=256)                       :: electric_pot_
LOGICAL                                  :: dipole_correction_, edir_ispresent, pomax_pos_ispresent,&
                                            dwn_width_ispresent, eamp_ispresent, nppstr_ispresent, & 
                                            nberry_cyc_ispresent, efield_vec_ispresent
INTEGER                                  :: efield_direction_
REAL(DP)                                 :: pot_max_pos_, down_width_, eamp_, efield_vec_(3)
INTEGER                                  :: nppstr_, nberry_cyc_
CHARACTER(iotk_attlenx)                  :: attr
! 
ispresent = .FALSE. 
CALL iotk_scan_begin( iunit, "electric_field", FOUND = ispresent, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
IF ( .NOT. ispresent ) RETURN 
! 
CALL iotk_scan_dat( iunit, "electric_potential", electric_pot_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat( iunit, "dipole_correction", dipole_correction_, DEFAULT = .FALSE., IERR = ierr )
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat( iunit, "electric_field_direction", efield_direction_, IERR = ierr , FOUND = edir_ispresent )
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat( iunit, "potential_max_position", pot_max_pos_, IERR = ierr , FOUND = pomax_pos_ispresent )
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat( iunit, "potential_decrease_width", down_width_, IERR = ierr, FOUND = dwn_width_ispresent )
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat( iunit, "electric_field_amplitude", eamp_, IERR = ierr, FOUND = eamp_ispresent )
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "electric_field_vector", efield_vec_, IERR = ierr, FOUND = efield_vec_ispresent)
IF ( ierr /= 0) RETURN  
CALL iotk_scan_dat( iunit, "nk_per_string", nppstr_, IERR = ierr, FOUND = nppstr_ispresent )
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat( iunit, "n_berry_cycles", nberry_cyc_, IERR = ierr, FOUND = nberry_cyc_ispresent) 
IF ( ierr /= 0) RETURN 
! 
CALL iotk_scan_end( iunit, "electric_field" , IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
!
CALL qes_init_electric_field( obj, "electric_field", electric_pot_, .TRUE., dipole_correction_, & 
                              edir_ispresent, efield_direction_, pomax_pos_ispresent, pot_max_pos_,& 
                              dwn_width_ispresent, down_width_, eamp_ispresent, eamp_, efield_vec_ispresent,&
                              efield_vec_, nppstr_ispresent, nppstr_, nberry_cyc_ispresent, nberry_cyc_)  
END SUBROUTINE qexsd_get_electric_field
! 
!------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_atomic_constraints(iunit, obj, ispresent)
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                          :: iunit
TYPE (atomic_constraints_type),INTENT(OUT)   :: obj
LOGICAL,INTENT(OUT)                          :: ispresent
!
INTEGER                                      :: ierr, num_of_constraints_, iconstr
REAL(DP)                                     :: tolerance_
TYPE ( atomic_constraint_type),ALLOCATABLE   :: at_constr_obj(:)
CHARACTER(iotk_attlenx)                      :: attr
!
ispresent = .FALSE.
CALL iotk_scan_begin(iunit, "atomic_constraints", IERR = ierr, FOUND = ispresent)
IF ( ierr /= 0 ) RETURN
IF ( .NOT. ispresent ) RETURN
!
CALL iotk_scan_dat( iunit, "num_of_constraints", num_of_constraints_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat( iunit, "tolerance", tolerance_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
DO iconstr =1, num_of_constraints_
   CALL qexsd_get_atomic_constraint( iunit, at_constr_obj(iconstr), ispresent ) 
   if ( .NOT. ispresent )  RETURN
END DO
CALL iotk_scan_end ( iunit, "atomic_constraints", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
CALL qes_init_atomic_constraints(obj, "atomic_constraints", num_of_constraints_, tolerance_, &
                                 num_of_constraints_, at_constr_obj)
DO iconstr =1 , num_of_constraints_
   CALL qes_reset_atomic_constraint(at_constr_obj(iconstr))
END DO
DEALLOCATE ( at_constr_obj )
END SUBROUTINE qexsd_get_atomic_constraints
! 
!------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_atomic_constraint( iunit, obj, ispresent)
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                          :: iunit
TYPE (atomic_constraint_type),INTENT(OUT)    :: obj
LOGICAL,INTENT(OUT)                          :: ispresent
!
INTEGER                                      :: ierr, iconstr
REAL(DP)                                     :: constr_target_, constr_parms_(4)
CHARACTER ( LEN = 256 )                      :: constr_type_
CHARACTER(iotk_attlenx)                      :: attr
ispresent = .FALSE.
CALL iotk_scan_begin(iunit, "atomic_constraint", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat(iunit, "constr_parms", constr_parms_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "constr_type", constr_type_, IERR = ierr )
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat( iunit, "constr_target", constr_target_, IERR = ierr )
IF ( ierr /= 0 ) RETURN 
!
CALL iotk_scan_end( iunit, "atomic_constraint", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
! 
CALL qes_init_atomic_constraint (obj, "atomic_constraint", constr_parms_,  constr_type_, constr_target_)
END SUBROUTINE qexsd_get_atomic_constraint
! 
!-------------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_spin_constraints( iunit, obj, ispresent)
!-----------------------------------------------------------------------------------------------------
 ! 
IMPLICIT NONE
! 
INTEGER, INTENT(IN)                          :: iunit
TYPE (spin_constraints_type),INTENT(OUT)     :: obj
LOGICAL,INTENT(OUT)                          :: ispresent
!
INTEGER                                      :: ierr
CHARACTER(LEN=256)                           :: spin_constraints_
REAL(DP)                                     :: lambda_, targ_mag_(3)
LOGICAL                                      :: targ_mag_ispresent
CHARACTER(iotk_attlenx)                      :: attr
! 
 
!
ispresent = .FALSE.
CALL iotk_scan_begin( iunit, "spin_constraints", IERR = ierr, FOUND = ispresent )
IF ( ierr /= 0 ) RETURN
IF (.NOT. ispresent ) RETURN
!
CALL iotk_scan_dat( iunit, "spin_constraints", spin_constraints_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat( iunit, "lagrange_multiplier", lambda_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "target_magnetization", targ_mag_, FOUND = targ_mag_ispresent, & 
     IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_end( iunit, "spin_constraints", IERR = ierr )
IF ( ierr /= 0) RETURN
!
CALL qes_init_spin_constraints( obj, "spin_constraints", spin_constraints_, lambda_, targ_mag_ispresent,&
                                targ_mag_) 
END SUBROUTINE qexsd_get_spin_constraints
!
!--------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_total_energy( iunit, obj, ispresent ) 
!--------------------------------------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                       :: iunit 
TYPE ( total_energy_type ),INTENT(OUT)    :: obj
LOGICAL,INTENT(OUT)                       :: ispresent
! 
INTEGER                                   :: ierr
REAL(DP)                                  :: etot_, eband_, ehart_, vtxc_, etxc_, ewald_, demet_, & 
                                             efield_corr_, potstat_contr_ 
LOGICAL                                   :: efield_corr_ispresent, eband_ispresent, ehart_ispresent, & 
                                             vtxc_ispresent, etxc_ispresent, ewald_ispresent, & 
                                             demet_ispresent, potstat_ispresent
CHARACTER(iotk_attlenx)                   :: attr
! 

! 
ispresent = .FALSE. 
CALL iotk_scan_begin ( iunit, "total_energy", IERR = ierr, FOUND = ispresent ) 
IF ( ierr /= 0 ) THEN 
   ispresent = .FALSE.
   RETURN 
END IF 
IF ( .NOT. ispresent ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "etot", etot_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "eband", eband_, IERR = ierr, FOUND = eband_ispresent ) 
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "ehart", ehart_, IERR = ierr, FOUND = ehart_ispresent ) 
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "vtxc", vtxc_, IERR = ierr, FOUND = vtxc_ispresent ) 
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "etxc", etxc_, IERR = ierr, FOUND = etxc_ispresent ) 
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "ewald", ewald_, IERR = ierr, FOUND = ewald_ispresent ) 
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "demet", demet_, IERR = ierr, FOUND = demet_ispresent ) 
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "efield_corr", efield_corr_, IERR = ierr, FOUND = efield_corr_ispresent) 
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "potentiostat_contr", potstat_contr_, IERR = ierr, FOUND = potstat_ispresent ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_end( iunit, "total_energy", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL qes_init_total_energy ( obj, "total_energy", etot_, eband_ispresent, eband_, ehart_ispresent, & 
                             ehart_, vtxc_ispresent, vtxc_, etxc_ispresent, etxc_, ewald_ispresent, & 
                            ewald_, demet_ispresent, demet_, efield_corr_ispresent, efield_corr_, potstat_ispresent,& 
                            potstat_contr_)  
END SUBROUTINE qexsd_get_total_energy
!
!------------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_scf_conv( iunit, obj, ispresent ) 
!----------------------------------------------------------------------------------------------------
!
IMPLICIT NONE
! 
INTEGER, INTENT(IN)                  :: iunit 
TYPE ( scf_conv_type ),INTENT(OUT)         :: obj
LOGICAL,INTENT(OUT)                  :: ispresent
! 
INTEGER                              :: ierr, nscf_steps_
REAL(DP)                             :: scf_error_
CHARACTER(iotk_attlenx)              :: attr
!
ispresent = .FALSE. 
CALL iotk_scan_begin ( iunit, "scf_conv", IERR = ierr, FOUND = ispresent ) 
IF ( ierr /= 0 ) RETURN
! 
IF ( .NOT. ispresent ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "n_scf_steps", nscf_steps_, IERR = ierr ) 
IF  ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "scf_error", scf_error_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
CALL iotk_scan_end ( iunit, "scf_conv", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
CALL  qes_init_scf_conv( obj, "scf_conv", nscf_steps_, scf_error_) 
! 
END SUBROUTINE qexsd_get_scf_conv
!-----------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_step( iunit, obj, ispresent) 
!-----------------------------------------------------------------------------------------------
! 
IMPLICIT NONE 
! 
INTEGER, INTENT(IN)                  :: iunit 
TYPE (step_type ),INTENT(OUT)        :: obj
LOGICAL,INTENT(OUT)                  :: ispresent
! 
INTEGER                              :: ierr, n_step_, nat_ 
LOGICAL                              :: found, stress_ispresent, fcp_force_ispresent, fcp_charge_ispresent
TYPE( scf_conv_type )                :: scf_conv_obj
TYPE( atomic_structure_type )        :: at_struct_obj
TYPE( total_energy_type )            :: tot_en_obj
TYPE ( matrix_type )                 :: for_mat_obj, stress_mat_obj
REAL(DP),ALLOCATABLE                 :: forces_mat(:,:)
REAL(DP)                             :: stress_mat(3,3), fcp_force_, fcp_tot_charge_
CHARACTER(iotk_attlenx)              :: attr
! 

! 
ispresent  = .FALSE.
CALL  iotk_scan_begin( iunit, "step", ATTR = attr, IERR = ierr , FOUND = ispresent )
IF ( ierr /= 0 ) RETURN 
IF (.NOT. ispresent ) RETURN 
!
CALL iotk_scan_attr( attr, "n_step", n_step_, IERR = ierr )
IF ( ierr /= 0) RETURN
!
CALL qexsd_get_scf_conv( iunit, scf_conv_obj, found ) 
IF ( .NOT. found ) RETURN 
!
CALL qexsd_get_atomic_structure( iunit, at_struct_obj, found) 
IF ( .NOT. found ) RETURN 
nat_ = at_struct_obj%nat
!
CALL qexsd_get_total_energy ( iunit, tot_en_obj, found ) 
IF ( .NOT. found ) RETURN 
! 
ALLOCATE ( forces_mat(3,nat_) )
CALL iotk_scan_dat( iunit, "forces", forces_mat, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
CALL qes_init_matrix( for_mat_obj, "forces", 3, nat_, forces_mat) 
!
CALL iotk_scan_dat( iunit, "stress", stress_mat, FOUND = stress_ispresent, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
IF ( stress_ispresent ) CALL qes_init_matrix( stress_mat_obj, "stress", 3, 3, stress_mat )
! 
CALL iotk_scan_dat( iunit, "FCP_force", fcp_force_, FOUND = fcp_force_ispresent, IERR = ierr ) 
IF (ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "FCP_tot_charge", fcp_tot_charge_, FOUND = fcp_charge_ispresent, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_end( iunit, "step", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
!
CALL qes_init_step( obj, "step", n_step_, scf_conv_obj, at_struct_obj, tot_en_obj, for_mat_obj, stress_ispresent,&
                    stress_mat_obj, fcp_force_ispresent, fcp_force_, fcp_charge_ispresent, fcp_tot_charge_ ) 
DEALLOCATE ( forces_mat )  
!
END SUBROUTINE qexsd_get_step
!
!---------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_opt_conv( iunit, obj, ispresent ) 
! 
IMPLICIT NONE
! 
INTEGER, INTENT(IN)                          :: iunit
TYPE( opt_conv_type ),INTENT(OUT)    :: obj
LOGICAL, INTENT ( OUT )                      :: ispresent
!
INTEGER                                      :: ierr, nopt_steps_
REAL(DP)                                     :: grad_norm_
CHARACTER(iotk_attlenx)                      :: attr
!
ispresent = .FALSE. 
CALL iotk_scan_begin ( iunit, "opt_conv", IERR = ierr, FOUND = ispresent ) 
IF ( ierr /= 0 ) RETURN 
IF ( .NOT. ispresent ) RETURN
!
CALL iotk_scan_dat ( iunit, "n_opt_steps", nopt_steps_, IERR = ierr) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "grad_norm", grad_norm_, IERR = ierr) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_end( iunit, "opt_conv", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL qes_init_opt_conv( obj, "opt_conv", nopt_steps_, grad_norm_)
END SUBROUTINE qexsd_get_opt_conv
!---------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_convergence_info( iunit, obj, ispresent ) 
!----------------------------------------------------------------------------------------------
! 
IMPLICIT NONE
! 
INTEGER, INTENT(IN)                          :: iunit
TYPE( convergence_info_type ),INTENT(OUT)    :: obj
LOGICAL, INTENT ( OUT )                      :: ispresent
!
INTEGER                                      :: ierr
LOGICAL                                      :: found, opt_conv_ispresent  
TYPE ( scf_conv_type )                       :: scf_conv_obj
TYPE ( opt_conv_type )                       :: opt_conv_obj
CHARACTER(iotk_attlenx)                      :: attr

!
ispresent = .FALSE. 
CALL iotk_scan_begin( iunit, "convergence_info", IERR = ierr , FOUND = ispresent) 
IF ( ierr /= 0 ) RETURN  
! 
CALL qexsd_get_scf_conv( iunit, scf_conv_obj, found )
IF ( .NOT. found ) RETURN
! 
CALL qexsd_get_opt_conv( iunit, opt_conv_obj, opt_conv_ispresent ) 
CALL iotk_scan_end( iunit, "convergence_info", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
! 
CALL  qes_init_convergence_info( obj, "convergence_info", scf_conv_obj, opt_conv_ispresent, &
                                 opt_conv_obj)
CALL qes_reset_scf_conv( scf_conv_obj) 
CALL qes_reset_opt_conv ( opt_conv_obj ) 
!  
END SUBROUTINE qexsd_get_convergence_info
!
!--------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_algorithmic_info( iunit, obj, ispresent ) 
! 
IMPLICIT NONE
! 
INTEGER, INTENT(IN)                          :: iunit
TYPE( algorithmic_info_type ),INTENT(OUT)    :: obj
LOGICAL, INTENT ( OUT )                      :: ispresent
!
INTEGER                                      :: ierr
LOGICAL                                      :: real_space_q_, uspp_, paw_
CHARACTER(iotk_attlenx)                      :: attr
! 
ispresent = .FALSE. 
CALL iotk_scan_begin( iunit, "algorithmic_info", IERR = ierr , FOUND = ispresent ) 
IF ( ierr /= 0 ) RETURN 
! 
IF ( .NOT. ispresent ) RETURN
! 
CALL iotk_scan_dat ( iunit, "real_space_q", real_space_q_, IERR = ierr)
IF ( ierr /= 0 ) RETURN 
!  
CALL iotk_scan_dat ( iunit, "uspp", uspp_, IERR = ierr )
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "paw", paw_, IERR = ierr )
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_end( iunit , "algorithmic_info", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
! 
CALL qes_init_algorithmic_info ( obj, "algorithmic_info", real_space_q_, uspp_, paw_ )
!
END SUBROUTINE qexsd_get_algorithmic_info
!
!-----------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_symmetry( iunit, obj, ispresent) 
!-------------------------------------------------------------------------------
   ! 
   IMPLICIT NONE
   ! 
   INTEGER, INTENT(IN)                          :: iunit
   TYPE( symmetry_type ),INTENT(OUT)          :: obj
   LOGICAL, INTENT ( OUT )                      :: ispresent
   !
   INTEGER                                      :: ierr, natoms_
   INTEGER,ALLOCATABLE                          :: equiv_atoms_(:)
   LOGICAL                                      :: found, frac_ispresent, equiv_at_ispresent,&
                                                   t_rev_, t_rev_is, symm_name_is, class_is
   REAL(DP)                                     :: mat_rot_(3,3), frac_tras_(3) 
   CHARACTER(LEN=256)                           :: symm_name_, info_str_, class_
   TYPE ( info_type )                           :: info_obj
   TYPE ( matrix_type )                         :: mat_rot_obj
   TYPE (equivalent_atoms_type )                :: eqv_at_obj
 CHARACTER(iotk_attlenx)                       :: attr
   !
   
   ispresent = .FALSE.
   CALL iotk_scan_begin ( iunit, "symmetry", IERR = ierr, FOUND = ispresent ) 
   IF ( ierr /= 0 ) RETURN 
   IF ( .NOT. ispresent ) RETURN
   ! 
   ! CALL qexsd_get_info( iunit, info_obj, found ) 
   !IF ( .NOT. found ) RETURN 
   CALL iotk_scan_dat(iunit, "info", DAT = info_str_ , ATTR = attr, IERR = ierr ) 
   IF ( ierr /= 0 ) RETURN 
   CALL iotk_scan_attr( attr, "name", symm_name_, IERR = ierr, FOUND = symm_name_is)
   IF ( ierr /= 0 ) RETURN
   CALL iotk_scan_attr( attr, "class", class_, IERR = ierr, FOUND = class_is)
   IF ( ierr /= 0 ) RETURN
   CALL iotk_scan_attr( attr, "time_reversal", t_rev_, IERR = ierr, FOUND = t_rev_is)
   IF ( ierr /= 0 ) RETURN
   CALL qes_init_info( info_obj, "info", symm_name_, symm_name_is, class_, class_is, t_rev_,&
       t_rev_is, info_str_ ) 
   ! 
   CALL iotk_scan_dat( iunit, "rotation", mat_rot_, IERR = ierr ) 
   IF ( ierr /= 0 ) RETURN 
   ! 
   CALL  qes_init_matrix( mat_rot_obj, "rotation", 3, 3, mat_rot_) 
   ! 
   CALL iotk_scan_dat( iunit, "fractional_translation", frac_tras_, IERR = ierr , &
       FOUND = frac_ispresent) 
   IF ( ierr /= 0 ) RETURN
   ! 
   CALL iotk_scan_begin( iunit, "equivalent_atoms", ATTR = attr, IERR = ierr , FOUND = equiv_at_ispresent)
   IF ( ierr /= 0 ) RETURN 
   IF ( equiv_at_ispresent ) THEN 
      CALL iotk_scan_end(iunit,"equivalent_atoms",IERR = ierr )
      IF ( ierr /= 0) RETURN
      CALL iotk_scan_attr( attr, "nat", natoms_, IERR = ierr ) 
      ALLOCATE ( equiv_atoms_(natoms_)) 
      CALL iotk_scan_dat(iunit, "equivalent_atoms", equiv_atoms_, IERR = ierr ) 
      IF ( ierr /= 0 ) RETURN
      CALL qes_init_equivalent_atoms( eqv_at_obj, "equivalent_atoms", natoms_, natoms_, equiv_atoms_) 
      DEALLOCATE ( equiv_atoms_)
   END IF 
   
   !
   CALL iotk_scan_end( iunit, "symmetry", IERR = ierr ) 
   IF ( ierr /= 0) RETURN
   !
   CALL qes_init_symmetry(obj, "symmetry", info_obj, mat_rot_obj, frac_ispresent, frac_tras_, &
          equiv_at_ispresent, eqv_at_obj) 
   CALL qes_reset_info(info_obj)
   CALL qes_reset_matrix(mat_rot_obj)
   IF (equiv_at_ispresent) CALL qes_reset_equivalent_atoms(eqv_at_obj)
END SUBROUTINE qexsd_get_symmetry
!---------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_symmetries( iunit, obj, ispresent ) 
!---------------------------------------------------------------------------------------
! 
IMPLICIT NONE 
! 
INTEGER, INTENT(IN)                          :: iunit
TYPE( symmetries_type ),INTENT(OUT)          :: obj
LOGICAL, INTENT ( OUT )                      :: ispresent
!
INTEGER                                      :: ierr, nsym_, nrot_, space_group_, isym
LOGICAL                                      :: found 
TYPE ( symmetry_type ),ALLOCATABLE           :: symm_objs(:)
CHARACTER(iotk_attlenx)                      :: attr
! 
ispresent = .FALSE. 
CALL iotk_scan_begin ( iunit, "symmetries", IERR = ierr , FOUND = ispresent) 
IF ( ierr /= 0 ) RETURN
IF ( .NOT. ispresent) RETURN
! 
CALL iotk_scan_dat ( iunit, "nsym", nsym_, IERR = ierr ) 
IF ( ierr /= 0   )  RETURN 
IF ( nsym_ .LT. 1 )  RETURN 
! 
CALL iotk_scan_dat ( iunit, "nrot", nrot_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "space_group", space_group_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
ALLOCATE (symm_objs(nrot_))
DO isym = 1, nrot_ 
   CALL qexsd_get_symmetry( iunit, symm_objs(isym), found) 
   IF ( .NOT. found ) RETURN
END DO
!
CALL iotk_scan_end( iunit, "symmetries", IERR = ierr ) 
IF  ( ierr /= 0 ) RETURN 
! 
CALL qes_init_symmetries ( obj, "symmetries", nsym_, nrot_, space_group_, nrot_, symm_objs ) 
DO isym =1, nrot_
   CALL qes_reset_symmetry( symm_objs(isym))
END DO
DEALLOCATE ( symm_objs )
END SUBROUTINE qexsd_get_symmetries
! 
!-----------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_basis_set( iunit, obj, ispresent )
!----------------------------------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                          :: iunit
TYPE( basis_set_type ),INTENT(OUT)           :: obj
LOGICAL, INTENT ( OUT )                      :: ispresent
!
INTEGER                                      :: ierr

!
LOGICAL                                      :: gamma_only_, gamma_only_ispresent, ecutrho_ispresent, found,&
                                                ngms_ispresent, box_ispresent, smooth_ispresent, grid_ispresent 
CHARACTER(LEN=256)                           :: empty_str
REAL(DP)                                     :: ecutrho_, ecutwfc_  
INTEGER                                      ::    fft_nr1_,    fft_nr2_,    fft_nr3_,&
                                                   box_nr1_,    box_nr2_,    box_nr3_,&
                                                   smooth_nr1_, smooth_nr2_, smooth_nr3_,&
                                                   ngm_, ngms_, npwx_ 
TYPE ( basisSetItem_type)                    :: grid_obj, smooth_obj, box_obj
TYPE ( reciprocal_lattice_type )             :: recip_obj 
CHARACTER(iotk_attlenx)                      :: attr
! 
ispresent = .FALSE. 
CALL iotk_scan_begin ( iunit, "basis_set", ATTR = attr, IERR = ierr, FOUND = ispresent ) 
IF (ierr /= 0) RETURN   
! 
IF ( .NOT. ispresent ) RETURN
gamma_only_ispresent = .FALSE. 
CALL iotk_scan_dat( iunit, "gamma_only", gamma_only_, FOUND = gamma_only_ispresent, IERR = ierr )
! 
CALL iotk_scan_dat ( iunit, "ecutwfc", ecutwfc_, IERR = ierr ) 
IF (ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "ecutrho", ecutrho_, IERR = ierr , FOUND = ecutrho_ispresent ) 
!
CALL iotk_scan_dat (iunit, "fft_grid", dat = empty_str, ATTR = attr, IERR = ierr , FOUND = grid_ispresent )  
IF ( ierr /= 0 ) RETURN
IF ( grid_ispresent ) THEN 
   CALL iotk_scan_attr(attr, "nr1", fft_nr1_, IERR = ierr ) 
   IF ( ierr /= 0) RETURN 
   CALL iotk_scan_attr( attr, "nr2", fft_nr2_, IERR = ierr ) 
   IF ( ierr /= 0) RETURN
   CALL iotk_scan_attr( attr, "nr3", fft_nr3_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL qes_init_basisSetItem ( grid_obj, "fft_grid", fft_nr1_, fft_nr2_, fft_nr3_ ,"")
ELSE
   RETURN  
END IF
! 

CALL iotk_scan_dat (iunit, "fft_smooth", dat = empty_str, ATTR = attr, &
                    FOUND = smooth_ispresent, IERR = ierr )
IF ( smooth_ispresent ) THEN  
   CALL iotk_scan_attr(attr, "nr1", smooth_nr1_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL iotk_scan_attr( attr, "nr2", smooth_nr2_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL iotk_scan_attr( attr, "nr3", smooth_nr3_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL qes_init_basisSetItem ( smooth_obj, "smooth_grid", smooth_nr1_, smooth_nr2_, smooth_nr3_ ,"")
END IF
!

CALL iotk_scan_dat (iunit, "fft_box", dat = empty_str, ATTR = attr, &
    FOUND = box_ispresent, IERR = ierr )
IF ( box_ispresent ) THEN  
   CALL iotk_scan_attr(attr, "nr1", box_nr1_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL iotk_scan_attr( attr, "nr2", box_nr2_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL iotk_scan_attr( attr, "nr3", box_nr3_, IERR = ierr )
   IF ( ierr /= 0) RETURN
   CALL qes_init_basisSetItem ( box_obj, "fft_box", box_nr1_, box_nr2_, box_nr3_ ,"")
END IF
!

CALL iotk_scan_dat(iunit, "ngm", ngm_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat(iunit, "ngms", ngms_, IERR = ierr, FOUND = ngms_ispresent ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat(iunit, "npwx", npwx_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!

CALL qexsd_get_reciprocal_lattice ( iunit, recip_obj, found) 
IF ( .NOT. found ) RETURN 
! 

CALL iotk_scan_end( iunit, "basis_set", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!

CALL  qes_init_basis_set ( obj, "basis_set", gamma_only_ispresent, gamma_only_, ecutwfc_, &
                         ecutrho_ispresent, ecutrho_, grid_obj, smooth_ispresent, & 
                         smooth_obj, box_ispresent, box_obj, ngm_, ngms_ispresent, ngms_, npwx_, recip_obj )
!
CALL qes_reset_basisSetItem(grid_obj)
IF (box_ispresent ) CALL qes_reset_basisSetItem(box_obj)
IF (smooth_ispresent ) CALL qes_reset_basisSetItem(smooth_obj)
CALL qes_reset_reciprocal_lattice ( recip_obj )
END SUBROUTINE qexsd_get_basis_set
!
!---------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_reciprocal_lattice( iunit, obj, ispresent ) 
!---------------------------------------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                          :: iunit
TYPE( reciprocal_lattice_type ),INTENT(OUT)           :: obj
LOGICAL, INTENT ( OUT )                      :: ispresent
!
INTEGER                                      :: ierr
REAL(DP)                                     :: b1_(3), b2_(3), b3_(3)
CHARACTER(iotk_attlenx)                      :: attr
!
ispresent = .FALSE.
CALL iotk_scan_begin ( iunit, "reciprocal_lattice", IERR = ierr , FOUND = ispresent)
IF ( ierr /= 0 ) RETURN
IF ( .NOT. ispresent ) RETURN
!
CALL iotk_scan_dat ( iunit, "b1", b1_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat ( iunit, "b2", b2_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "b3", b3_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
! 

CALL iotk_scan_end(iunit, "reciprocal_lattice", IERR = ierr ) 
IF ( ierr /= 0  ) RETURN
!
CALL qes_init_reciprocal_lattice ( obj, "reciprocal_lattice", b1_, b2_, b3_)
!
END SUBROUTINE qexsd_get_reciprocal_lattice
!
!--------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_magnetization ( iunit, obj, ispresent ) 
!--------------------------------------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                          :: iunit
TYPE( magnetization_type ),INTENT(OUT)       :: obj
LOGICAL, INTENT ( OUT )                      :: ispresent
!
INTEGER                                      :: ierr
LOGICAL                                      :: lsda_, noncolin_, spinorbit_, do_magnetization_
REAL(DP)                                     :: total_, absolute_ 
CHARACTER(iotk_attlenx)                      :: attr
!
ispresent = .FALSE.
CALL iotk_scan_begin( iunit, "magnetization", IERR = ierr, FOUND = ispresent) 
IF ( ierr /= 0 ) RETURN 
IF ( .NOT. ispresent ) RETURN 
!
CALL iotk_scan_dat ( iunit, "lsda", lsda_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "noncolin", noncolin_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "spinorbit", spinorbit_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "total", total_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "absolute", absolute_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "do_magnetization", do_magnetization_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_end( iunit, "magnetization", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL qes_init_magnetization ( obj, "magnetization", lsda_, noncolin_, spinorbit_, total_, absolute_, &
                             do_magnetization_)
!
END SUBROUTINE qexsd_get_magnetization
!
!--------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_ks_energies ( iunit,  obj, ispresent, lsda, noncolin , nbnd, nbnd_up, nbnd_dw)
!--------------------------------------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER,INTENT(IN)                          :: iunit, nbnd
INTEGER,INTENT(IN),OPTIONAL                 :: nbnd_up, nbnd_dw
TYPE (ks_energies_type),INTENT(OUT)         :: obj
LOGICAL,INTENT(IN)                          :: lsda, noncolin
LOGICAL,INTENT(OUT)                         :: ispresent
!
INTEGER                                     :: ierr, nbnd_aux_up, nbnd_aux_dw, ndim_eigenval, npw_
REAL(DP),ALLOCATABLE                        :: eiv_(:), occ_(:)
REAL(DP)                                    :: xk_(3), wk_
TYPE ( k_point_type )                       :: kp_obj
LOGICAL                                     :: wk_ispresent
CHARACTER(iotk_attlenx)                     :: attr
!
ispresent = .FALSE.
CALL iotk_scan_begin( iunit, "ks_energies", IERR = ierr, FOUND = ispresent )
IF  ( ierr /= 0 ) RETURN
IF ( .NOT. ispresent ) RETURN
! 
CALL iotk_scan_dat( iunit, "k_point", xk_, ATTR = attr, IERR = ierr )
IF ( ierr /= 0 ) RETURN
!F
CALL iotk_scan_attr(attr, "weight", wk_, IERR = ierr, FOUND = wk_ispresent ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL qes_init_k_point ( kp_obj, "k_point", wk_, wk_ispresent, "", .FALSE., xk_)
! 
CALL iotk_scan_dat ( iunit, "npw", npw_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
!
IF ( lsda ) THEN 
   IF ( PRESENT ( nbnd_up ) ) THEN 
      nbnd_aux_up = nbnd_up
   ELSE 
      nbnd_aux_up = nbnd
   END IF 
   IF ( PRESENT ( nbnd_dw ) ) THEN 
      nbnd_aux_dw = nbnd_dw
   ELSE 
      nbnd_aux_dw = nbnd
   END IF  
   ndim_eigenval = nbnd_aux_up+nbnd_aux_dw
ELSE
   ndim_eigenval = nbnd
END IF 
ALLOCATE (eiv_(ndim_eigenval), occ_(ndim_eigenval))
!
! 
CALL iotk_scan_dat( iunit, "eigenvalues", eiv_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
!
CALL iotk_scan_dat ( iunit, "occupations", occ_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_end ( iunit, "ks_energies", IERR = ierr ) 
IF ( ierr /= 0) RETURN
! 
CALL qes_init_ks_energies( obj, "ks_energies", kp_obj, npw_, ndim_eigenval, eiv_, ndim_eigenval, occ_)
! 
DEALLOCATE ( eiv_, occ_)
CALL qes_reset_k_point( kp_obj)
! 
END SUBROUTINE qexsd_get_ks_energies  
!--------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_band_structure( iunit, obj, ispresent ) 
!--------------------------------------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                          :: iunit
TYPE( band_structure_type ),INTENT(OUT)      :: obj
LOGICAL, INTENT ( OUT )                      :: ispresent
!
INTEGER                                      :: ierr, ik, nbnd_, nbnd_up_, nbnd_dw_, nks_, n_wfc_at_, occ_spin_
LOGICAL                                      :: found, lsda_, noncolin_, spinorbit_, nbdup_ispresent, &
                                                nbddw_ispresent, fermi_energy_ispresent, ks_eng_found,&
                                                two_fermi_energies_ispresent, HOL_ispresent, n_wfc_at_ispresent,&
                                                spin_occ_is, smearing_ispresent, wf_collected_
REAL(DP)                                     :: nelec_, fermi_energy_,ef_updw_(2), HOL_energy_
TYPE ( ks_energies_type ),ALLOCATABLE        :: ks_energies_obj(:)                    
TYPE ( k_points_IBZ_type)                    :: starting_k_points_ 
TYPE ( occupations_type )                    :: occupations_obj
TYPE ( smearing_type )                       :: smearing_obj 
CHARACTER(LEN = 256)                         :: occupations_string_
CHARACTER(iotk_attlenx)                      :: attr
!
ispresent = .FALSE.
CALL iotk_scan_begin ( iunit , "band_structure", IERR = ierr, FOUND = ispresent ) 
IF ( ierr /= 0 ) RETURN 
IF (.NOT. ispresent ) RETURN
! 
CALL iotk_scan_dat ( iunit, "lsda", lsda_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat ( iunit, "noncolin", noncolin_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat ( iunit, "spinorbit", spinorbit_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat ( iunit, "nbnd", nbnd_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat ( iunit, "nbnd_up", nbnd_up_, IERR = ierr , FOUND = nbdup_ispresent) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat ( iunit, "nbnd_dw", nbnd_dw_, IERR = ierr , FOUND = nbddw_ispresent ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat ( iunit, "nelec", nelec_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat ( iunit, "num_of_atomic_wfc", n_wfc_at_, IERR = ierr, FOUND = n_wfc_at_ispresent ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL iotk_scan_dat ( iunit, "wf_collected", wf_collected_, IERR = ierr ) 
IF ( ierr/=0) CALL errore ( "qexsd_get_band_structure", "wf_collected not found", 1) 
!
CALL iotk_scan_dat ( iunit, "fermi_energy", fermi_energy_, IERR = ierr, FOUND = fermi_energy_ispresent ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_dat ( iunit, "highestOccupiedLevel", hol_energy_, IERR = ierr , FOUND = HOL_ispresent ) 
IF ( ierr /= 0 ) RETURN 
!
CALL iotk_scan_dat( iunit, "two_fermi_energies", ef_updw_, IERR = ierr, FOUND = two_fermi_energies_ispresent )
IF ( ierr /= 0 ) RETURN 
!
CALL qexsd_get_k_points_IBZ(iunit, starting_k_points_, ISPRESENT = found, TAGNAME = "starting_k_points")
IF ( .NOT. found ) CALL errore ( "qexsd_get_band_structure", "starting_k_points element not found", 1)    
!   
CALL iotk_scan_dat ( iunit, "nks", nks_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
!
CALL iotk_scan_dat ( iunit, "occupations_kind",  occupations_string_, ATTR = attr, IERR = ierr )
IF ( ierr /= 0 ) CALL errore ( "qexsd_get_band_structure", "error reading occupations_kind element", 1) 
!
CALL iotk_scan_attr (attr, "spin", occ_spin_, IERR = ierr , FOUND = spin_occ_is) 
IF  ( ierr /= 0 ) CALL errore ( "qexsd_get_band_structure", "error reading occupations_kind spin element", 1 )
!
CALL qes_init_occupations( occupations_obj, "occupations_kind", occ_spin_, spin_occ_is, occupations_string_) 
!
CALL qexsd_get_smearing ( iunit, smearing_obj, ISPRESENT = smearing_ispresent ) 
ALLOCATE ( ks_energies_obj(nks_) )
DO ik = 1, nks_
   IF ((.NOT. nbdup_ispresent) .AND. (.NOT. nbddw_ispresent )) THEN  
     CALL qexsd_get_ks_energies ( iunit, ks_energies_obj(ik), ks_eng_found , lsda_, noncolin_, nbnd_)
   ELSE IF ( nbdup_ispresent .AND. nbddw_ispresent ) THEN 
     CALL qexsd_get_ks_energies ( iunit, ks_energies_obj(ik), ks_eng_found , lsda_, noncolin_, nbnd_, &
                                  nbnd_up_, nbnd_dw_)
   ELSE IF (nbdup_ispresent ) THEN 
     CALL qexsd_get_ks_energies ( iunit, ks_energies_obj(ik), ks_eng_found , lsda_, noncolin_, nbnd_, &
     NBND_UP = nbnd_up_)
   ELSE IF (nbddw_ispresent ) THEN 
     CALL qexsd_get_ks_energies ( iunit, ks_energies_obj(ik), ks_eng_found , lsda_, noncolin_, nbnd_, &
     NBND_DW = nbnd_dw_)
   END IF
   IF ( .NOT. ks_eng_found ) RETURN 
END DO
! 
CALL iotk_scan_end( iunit, "band_structure", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
!
CALL qes_init_band_structure ( obj, "band_structure", lsda_, noncolin_, spinorbit_, nbnd_, nbdup_ispresent,          &
       nbnd_up_, nbddw_ispresent, nbnd_dw_, nelec_, n_wfc_at_ispresent, n_wfc_at_, wf_collected_,                    &
       fermi_energy_ispresent, fermi_energy_,  HOL_ispresent, HOL_energy_, two_fermi_energies_ispresent, 2, ef_updw_,& 
       starting_k_points_, SMEARING_ISPRESENT = smearing_ispresent, SMEARING = smearing_obj, NKS = nks_,             &
       OCCUPATIONS_KIND = occupations_obj, NDIM_KS_ENERGIES = nks_,   KS_ENERGIES = ks_energies_obj)
!  
CALL qes_reset_k_points_IBZ(starting_k_points_) 
CALL qes_reset_occupations(occupations_obj) 
DO ik = 1, nks_
   CALL qes_reset_ks_energies(ks_energies_obj(ik))
END DO
DEALLOCATE ( ks_energies_obj )

END SUBROUTINE qexsd_get_band_structure
!--------------------------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_polarization ( iunit, obj, ispresent )
!--------------------------------------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER,INTENT(IN)                      :: iunit
TYPE( polarization_type ),INTENT(OUT)   :: obj
LOGICAL,INTENT(OUT)                     :: ispresent
!
INTEGER                                 :: ierr 
REAL(DP)                                :: polarization_, modulus_, direction_(3)
CHARACTER(LEN=256)                      :: units_
TYPE ( scalarQuantity_type)             :: polarization_obj
CHARACTER(iotk_attlenx)                 :: attr
!
ispresent = .FALSE.
CALL iotk_scan_begin(iunit, "polarization", IERR = ierr, FOUND = ispresent ) 
IF ( ierr /= 0 ) RETURN
IF (.NOT. ispresent ) RETURN
! 
CALL iotk_scan_dat(iunit, "polarization", polarization_, ATTR = attr, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
CALL iotk_scan_attr(attr, "Units", units_, IERR = ierr ) 
IF (ierr /= 0 ) RETURN
CALL qes_init_scalarQuantity(polarization_obj, "polarization", TRIM(units_), polarization_)
!
CALL iotk_scan_dat ( iunit, "modulus", modulus_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
! 
CALL iotk_scan_dat ( iunit, "direction", direction_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_end(iunit, "polarization", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL  qes_init_polarization(obj, "polarization", polarization_obj, modulus_, direction_)
CALL qes_reset_scalarQuantity(polarization_obj)
!
END SUBROUTINE qexsd_get_polarization
!
!-------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_phase(iunit, tag, obj, ispresent ) 
!-------------------------------------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                              :: iunit
CHARACTER(LEN=*),INTENT(IN)                      :: tag
TYPE (phase_type ),INTENT (OUT)                  :: obj
LOGICAL,INTENT (OUT)                              :: ispresent
!
INTEGER                                          :: ierr 
REAL(DP)                                         :: phase_, ionic_, electronic_
CHARACTER(LEN=256)                               :: modulus_
LOGICAL                                          :: ionic_ispresent, electronic_ispresent, &
                                                    modulus_ispresent
CHARACTER(iotk_attlenx)                          :: attr
!
ispresent =.FALSE.
CALL iotk_scan_dat(iunit, TRIM(tag), phase_, ATTR = attr, IERR = ierr , FOUND = ispresent)
IF ( ierr /= 0 ) RETURN
IF (.NOT. ispresent ) RETURN 
!
CALL iotk_scan_attr ( attr, "ionic", ionic_, IERR = ierr , FOUND = ionic_ispresent)
IF ( ierr /= 0 ) RETURN 
CALL iotk_scan_attr ( attr, "electronic", electronic_, IERR = ierr , FOUND = electronic_ispresent)
IF ( ierr /= 0 ) RETURN 
CALL iotk_scan_attr( attr, "modulus", modulus_, IERR = ierr, FOUND = modulus_ispresent ) 
IF ( ierr /= 0 ) RETURN 
!
CALL qes_init_phase( obj, TRIM(tag), ionic_, ionic_ispresent, electronic_, electronic_ispresent, &
                     TRIM(modulus_), modulus_ispresent, phase_)
END SUBROUTINE qexsd_get_phase 
!
!--------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_ionic_polarization ( iunit, obj, ispresent) 
!--------------------------------------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                         :: iunit
TYPE ( ionicPolarization_type ),INTENT(OUT) :: obj 
LOGICAL,INTENT(OUT)                         :: ispresent
!
INTEGER                                     :: ierr 
LOGICAL                                     :: found
REAL(DP)                                    :: charge_
TYPE ( atom_type )                          :: ion_obj
TYPE ( phase_type )                         :: phase_obj
CHARACTER(iotk_attlenx)                     :: attr
! 
ispresent = .FALSE.
CALL iotk_scan_begin(iunit, "ionicPolarization", IERR = ierr, FOUND = ispresent ) 
IF ( ierr /= 0 ) RETURN 
IF ( .NOT. ispresent ) RETURN
! 
CALL qexsd_get_atom( iunit, "ion", ion_obj, found)
IF (.NOT. found ) RETURN 
CALL iotk_scan_dat ( iunit, "charge", charge_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
CALL qexsd_get_phase( iunit, "phase", phase_obj, found ) 
IF (.NOT. found ) RETURN 
! 
CALL iotk_scan_end(iunit, "ionicPolarization", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
CALL qes_init_ionicPolarization(obj, "ionicPolarization", ion_obj, charge_, phase_obj)
CALL qes_reset_atom(ion_obj)
CALL qes_reset_phase(phase_obj)
END SUBROUTINE qexsd_get_ionic_polarization
!
!-------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_atom(iunit, tag, obj, ispresent ) 
!-------------------------------------------------------------------------------------------------
!
IMPLICIT NONE
! 
INTEGER,INTENT (IN)                             :: iunit
CHARACTER(LEN=*),INTENT(IN)                     :: tag
TYPE ( atom_type ),INTENT(OUT)                  :: obj
LOGICAL,INTENT(OUT)                             :: ispresent
! 
INTEGER                                         :: ierr, index_ 
REAL(DP)                                        :: the_coordinates(3)
CHARACTER(LEN=256)                              :: label_
LOGICAL                                         :: index_ispresent
CHARACTER(iotk_attlenx)                         :: attr
! 
ispresent = .FALSE.
CALL iotk_scan_dat(iunit, TRIM(tag), ATTR = attr, dat = the_coordinates, IERR = ierr, FOUND = ispresent ) 
IF (ierr /= 0) RETURN
IF (.NOT. ispresent ) RETURN 
CALL iotk_scan_attr(attr, "name", label_, IERR = ierr ) 
IF (ierr /= 0) RETURN 
CALL iotk_scan_attr( attr, "index", index_, IERR = ierr, FOUND = index_ispresent ) 
IF ( ierr /= 0 ) RETURN
! 
CALL qes_init_atom(obj, TRIM(tag), label_, "", .FALSE. , index_, index_ispresent,  the_coordinates)
END SUBROUTINE qexsd_get_atom
!-------------------------------------------------------------------------------------------------- 
SUBROUTINE qexsd_get_berryPhaseOutput( iunit, obj, ispresent)
!-------------------------------------------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                          :: iunit
TYPE ( berryPhaseOutput_type ),INTENT(OUT)    :: obj
LOGICAL, INTENT(OUT)                         :: ispresent 
! 
INTEGER                                      :: ierr, iobj, ndim_el_pol, ndim_ion_pol, pippo 
LOGICAL                                      :: found
TYPE ( polarization_type )                   :: polarization_obj
TYPE ( phase_type )                          :: total_phase_obj
!
TYPE(ionicPolarization_type)              :: ion_pol_aux
TYPE(ionicPolarization_type),ALLOCATABLE  :: ion_pol_obj(:)
!
TYPE  ionicPolarization_list 
  TYPE ( ionicPolarization_type)          :: obj
  TYPE (ionicPolarization_list),POINTER   :: next 
END TYPE ionicPolarization_list
!
TYPE ( ionicPolarization_list),TARGET     :: ion_pol_list
TYPE (ionicPolarization_list),POINTER     :: ion_pol_first, ion_pol_last, ion_pol_ptr    
!
TYPE (electronicPolarization_type)              :: el_pol_aux
TYPE ( electronicPolarization_type),ALLOCATABLE  :: el_pol_obj(:)
TYPE electronicPolarization_list
   TYPE ( electronicPolarization_type )        :: obj 
   TYPE ( electronicPolarization_list), POINTER :: next 
END TYPE electronicPolarization_list
! 
TYPE ( electronicPolarization_list ),TARGET     :: el_pol_list
TYPE ( electronicPolarization_list ),POINTER    :: el_pol_first, el_pol_last, el_pol_ptr
CHARACTER(iotk_attlenx)                      :: attr
!
ispresent = .FALSE.
CALL iotk_scan_begin(iunit, "BerryPhase", IERR = ierr, FOUND = ispresent ) 
IF ( ierr /= 0 ) RETURN
IF ( .NOT. ispresent ) RETURN
! 
CALL qexsd_get_polarization( iunit, polarization_obj, found ) 
IF ( .NOT. found ) CALL errore ( "qexsd_get_berry_phase", "failed to get polarization field", 1)
!
CALL qexsd_get_phase( iunit, "totalPhase", total_phase_obj,  found  )
IF ( .NOT. found ) CALL errore ( "qexsd_get_berry_phase", "failed to get phase field", 1)
!
CALL qexsd_get_ionic_polarization ( iunit, ion_pol_aux, found ) 
IF ( .NOT. found ) CALL errore ( "qexsd_get_berry_phase", "failed to get ionic polarization field", 1)
!
ndim_ion_pol = 1
ion_pol_last => ion_pol_list 
conta_ion_polarization:DO
  ion_pol_last%obj = ion_pol_aux
  CALL qexsd_get_ionic_polarization ( iunit, ion_pol_aux, found ) 
  IF (ALL(ion_pol_last%obj%ion%atom == ion_pol_aux%ion%atom) ) EXIT conta_ion_polarization
  ndim_ion_pol = ndim_ion_pol+1
  ALLOCATE (ion_pol_last%next)
  ion_pol_last => ion_pol_last%next
  ion_pol_last%next => null (ion_pol_last)
END DO conta_ion_polarization
ALLOCATE ( ion_pol_obj(ndim_ion_pol))

iobj =1 
ion_pol_obj( 1 ) = ion_pol_list%obj
IF (ASSOCIATED( ion_pol_list%next)) THEN 
   ion_pol_first => ion_pol_list%next
   copy_ion_pol:DO
      iobj = iobj + 1 
      ion_pol_obj( iobj ) = ion_pol_first%obj
      IF ( ASSOCIATED( ion_pol_first%next )) THEN 
         ion_pol_ptr => ion_pol_first
         ion_pol_first => ion_pol_ptr%next
         DEALLOCATE( ion_pol_ptr) 
      ELSE
         DEALLOCATE ( ion_pol_first ) 
         EXIT copy_ion_pol
      END IF
      
   END DO copy_ion_pol
END IF 
CALL qexsd_get_electronicPolarization ( iunit, el_pol_aux, found ) 
IF (.NOT. found ) CALL errore ( "qexsd_get_berry_phase", "failed to get elec. polarization field", 1)
!
ndim_el_pol = 1
el_pol_last => el_pol_list 
count_el_pol:DO
  el_pol_last%obj = el_pol_aux
  CALL qexsd_get_electronicPolarization ( iunit, el_pol_aux, found ) 
  IF (ALL(el_pol_last%obj%firstKeyPoint%k_point == el_pol_aux%firstKeyPoint%k_point) ) THEN 
     IF (el_pol_last%obj%spin_ispresent ) THEN
        IF (el_pol_last%obj%spin == el_pol_aux%spin  ) EXIT count_el_pol
     ELSE 
        EXIT count_el_pol
     END IF
  END IF 
  ndim_el_pol = ndim_el_pol+1
  ALLOCATE (el_pol_last%next)
  el_pol_last => el_pol_last%next
  el_pol_last%next => null (el_pol_last)
END DO count_el_pol
ALLOCATE ( el_pol_obj(ndim_el_pol))
iobj =1 
el_pol_obj( 1 ) = el_pol_list%obj
IF (ASSOCIATED( el_pol_list%next)) THEN 
   el_pol_first => el_pol_list%next
   copy_el_pol:DO
      iobj = iobj + 1 
      el_pol_obj( iobj ) = el_pol_first%obj
      IF ( ASSOCIATED( el_pol_first%next )) THEN 
         el_pol_ptr => el_pol_first
         el_pol_first => el_pol_ptr%next
         DEALLOCATE( el_pol_ptr) 
      ELSE
         DEALLOCATE ( el_pol_first ) 
         EXIT copy_el_pol
      END IF
   END DO copy_el_pol
END IF 
!
CALL iotk_scan_end(iunit, "BerryPhase", IERR = ierr ) 
IF ( ierr /= 0 ) CALL errore ( "qexsd_get_BerryPhase", "iotk_scan_end failed." , 2 )
! 
CALL qes_init_BerryPhaseOutput(obj, "BerryPhase", polarization_obj, total_phase_obj, ndim_ion_pol, &
                         ion_pol_obj, ndim_el_pol, el_pol_obj)
CALL qes_reset_polarization(polarization_obj) 
CALL qes_reset_phase( total_phase_obj)
DO iobj = 1, ndim_ion_pol
   CALL qes_reset_ionicPolarization(ion_pol_obj(iobj))
END DO
DEALLOCATE ( ion_pol_obj)
DO iobj =1, ndim_el_pol
   CALL qes_reset_electronicPolarization(el_pol_obj(iobj))
END DO
DEALLOCATE ( el_pol_obj)
!
END SUBROUTINE qexsd_get_berryPhaseOutput
!
!-------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_electronicPolarization( iunit, obj, ispresent )
!-------------------------------------------------------------------------------------------------
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                               :: iunit
TYPE ( electronicPolarization_type ),INTENT(OUT)  :: obj
LOGICAL,INTENT(OUT)                               :: ispresent
!
INTEGER                                           :: ierr, ispin_ 
REAL(DP)                                          :: xk_(3), wk_
TYPE ( k_point_type)                              :: kp_obj
TYPE ( phase_type )                               :: phase_obj
LOGICAL                                           :: found, weight_ispresent, spin_ispresent
CHARACTER(iotk_attlenx)                           :: attr
! 
ispresent = .FALSE.
CALL iotk_scan_begin ( iunit, "electronicPolarization", IERR = ierr , FOUND = ispresent ) 
IF ( ierr  /= 0 ) RETURN
IF (.NOT. ispresent ) RETURN 
!
CALL iotk_scan_dat(iunit, "firstKeyPoint", xk_, ATTR = attr, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
wk_ = 0.d0
CALL iotk_scan_attr(attr, "weight", wk_, FOUND = weight_ispresent ) 
CALL qes_init_k_point (kp_obj, "firstKeyPoint", wk_, weight_ispresent, LABEL = "", LABEL_ISPRESENT = .FALSE., &
                       K_POINT =  xk_) 
!
CALL iotk_scan_dat ( iunit, "spin", ispin_, IERR = ierr, FOUND = spin_ispresent )
!
CALL qexsd_get_phase(iunit, "phase", phase_obj, found ) 
IF (.NOT. found ) RETURN
CALL iotk_scan_end(iunit, "electronicPolarization", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
CALL qes_init_electronicPolarization ( obj, "electronicPolarization", kp_obj, spin_ispresent, ispin_, &
                                       phase_obj )  
CALL qes_reset_k_point( kp_obj) 
CALL qes_reset_phase ( phase_obj ) 
END SUBROUTINE qexsd_get_electronicPolarization 
!
!------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_finiteFieldOut (iunit, obj, ispresent ) 
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)                          :: iunit
TYPE ( finiteFieldOut_type),INTENT(OUT) :: obj
LOGICAL, INTENT(OUT)                         :: ispresent 
! 
INTEGER                                      :: ierr 
REAL(DP)                                     :: el_dip_(3), ion_dip_(3)
CHARACTER(iotk_attlenx)                      :: attr
!
ispresent = .FALSE.
CALL iotk_scan_begin( iunit, "finiteElectricFieldInfo", IERR = ierr,  FOUND = ispresent ) 
IF (  ierr /= 0 ) RETURN
IF (.NOT. ispresent ) RETURN
!
CALL iotk_scan_dat ( iunit, "electronicDipole", el_dip_, IERR =  ierr )
IF ( ierr /= 0 ) RETURN
! 
CALL  iotk_scan_dat ( iunit, "ionicDipole", ion_dip_,  IERR = ierr )
IF ( ierr /= 0 ) RETURN
!
CALL iotk_scan_end(iunit, "finiteElectricFieldInfo", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
CALL qes_init_finiteFieldOut(obj, "finiteElectricFieldInfo", el_dip_, ion_dip_ )
END SUBROUTINE qexsd_get_finiteFieldOut
!-------------------------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_scalarQuantity ( iunit, tag, obj, ispresent ) 
!-------------------------------------------------------------------------------------------------
! 
IMPLICIT NONE
!
INTEGER, INTENT(IN)                       :: iunit
CHARACTER(LEN=*),INTENT(IN)               :: tag
TYPE ( scalarQuantity_type ),INTENT(OUT)  :: obj
LOGICAL, INTENT(OUT)                      :: ispresent 
! 
INTEGER                                   :: ierr 
REAL(DP)                                  :: quantity_
CHARACTER(LEN=256)                        :: units_
CHARACTER(iotk_attlenx)                   :: attr
! 
ispresent = .FALSE.
CALL iotk_scan_dat ( iunit, TRIM(tag), quantity_, ATTR = attr, IERR = ierr, FOUND = ispresent )
IF ( ierr /= 0 ) RETURN
IF (.NOT. ispresent ) RETURN 
!
CALL iotk_scan_attr(attr, "Units", units_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
CALL qes_init_scalarQuantity( obj, TRIM(tag), units_, quantity_)
END SUBROUTINE qexsd_get_scalarQuantity  
!-------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_dipoleOutput( iunit, obj, ispresent ) 
!
IMPLICIT NONE
!
INTEGER,INTENT(IN)                        :: iunit
TYPE ( dipoleOutput_type ),INTENT(OUT)    :: obj
LOGICAL,INTENT(OUT)                       :: ispresent 
!
INTEGER                                   :: ierr 
TYPE (scalarQuantity_type)                :: dipole_obj, ion_dipole_obj, elec_dipole_obj, & 
                                             dipoleField_obj, potentialAmp_obj, totalLength_obj
INTEGER                                   :: idir_
LOGICAL                                   :: found
CHARACTER(iotk_attlenx)                   :: attr
! 
ispresent = .FALSE.
CALL iotk_scan_begin( iunit, "dipoleInfo", IERR = ierr, FOUND = ispresent ) 
IF ( ierr /= 0 ) RETURN
IF (.NOT. ispresent ) RETURN 
!
CALL iotk_scan_dat(iunit, "idir", idir_, IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
!
CALL qexsd_get_scalarQuantity(iunit, "dipole", dipole_obj, found)
IF (.NOT. found ) RETURN 
! 
CALL qexsd_get_scalarQuantity(iunit, "ion_dipole", ion_dipole_obj, found)
IF (.NOT. found ) RETURN 
! 
CALL qexsd_get_scalarQuantity(iunit, "elec_dipole",elec_dipole_obj, found)
IF (.NOT. found ) RETURN 
! 
CALL qexsd_get_scalarQuantity(iunit, "dipoleField", dipoleField_obj, found)
IF (.NOT. found ) RETURN 
! 
CALL qexsd_get_scalarQuantity(iunit, "potentialAmp", potentialAmp_obj, found)
IF (.NOT. found ) RETURN 
! 
CALL qexsd_get_scalarQuantity(iunit, "totalLength", totalLength_obj, found)
IF (.NOT. found ) RETURN 
! 
CALL iotk_scan_end ( iunit, "dipoleInfo", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN
!
CALL qes_init_dipoleOutput (obj, "dipoleInfo", idir_, dipole_obj, ion_dipole_obj, elec_dipole_obj,&
                             dipoleField_obj, potentialAmp_obj, totalLength_obj )
CALL qes_reset_scalarQuantity ( dipole_obj)  
CALL qes_reset_scalarQuantity ( ion_dipole_obj)  
CALL qes_reset_scalarQuantity ( elec_dipole_obj)  
CALL qes_reset_scalarQuantity ( dipoleField_obj)  
CALL qes_reset_scalarQuantity ( potentialAmp_obj)  
CALL qes_reset_scalarQuantity ( totalLength_obj)
END SUBROUTINE qexsd_get_dipoleOutput
!  
!-------------------------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_outputElectricField( iunit, obj, ispresent ) 
!-------------------------------------------------------------------------------------------------
IMPLICIT NONE
!
INTEGER, INTENT(IN)                          :: iunit
TYPE ( outputElectricField_type),INTENT(OUT) :: obj
LOGICAL, INTENT(OUT)                         :: ispresent 
! 
INTEGER                                      :: ierr 
LOGICAL                                      :: bp_out_ispresent, finite_field_ispresent, &
                                                dipole_out_ispresent
TYPE ( BerryPhaseOutput_type)                :: bp_out_obj
TYPE ( finiteFieldOut_type )                 :: finite_field_obj
TYPE ( dipoleOutput_type )                   :: dipole_out_obj
CHARACTER(iotk_attlenx)                      :: attr
! 
ispresent = .FALSE. 
CALL iotk_scan_begin ( iunit, "electric_field", IERR = ierr, FOUND = ispresent ) 
IF ( ierr /= 0 ) RETURN
IF (.NOT. ispresent ) RETURN 
! 
CALL qexsd_get_BerryPhaseOutput( iunit, bp_out_obj, bp_out_ispresent) 
! 
CALL qexsd_get_finiteFieldOut ( iunit, finite_field_obj, finite_field_ispresent ) 
! 
CALL qexsd_get_dipoleOutput ( iunit, dipole_out_obj, dipole_out_ispresent ) 

CALL iotk_scan_end ( iunit, "electric_field", IERR = ierr ) 
IF ( ierr /= 0 ) RETURN 
! 
CALL qes_init_outputElectricField( obj, "electric_field", bp_out_ispresent, bp_out_obj, &
                                   finite_field_ispresent, finite_field_obj, dipole_out_ispresent, &
                                   dipole_out_obj ) 
! 
IF ( bp_out_ispresent ) CALL qes_reset_BerryPhaseOutput(bp_out_obj)
IF ( finite_field_ispresent ) CALL qes_reset_finiteFieldOut(finite_field_obj)
IF ( dipole_out_ispresent )  CALL qes_reset_dipoleOutput( dipole_out_obj ) 
END SUBROUTINE qexsd_get_outputElectricField
!
!--------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_forces( iunit, obj, natoms, ispresent ) 
! 
IMPLICIT NONE 
! 
INTEGER,INTENT (IN)                   :: iunit, natoms
TYPE ( matrix_type ), INTENT (OUT)    :: obj
LOGICAL,INTENT (OUT)                  :: ispresent
! 
INTEGER                               :: ierr 
REAL(DP), ALLOCATABLE                 :: force_mat_(:,:) 
CHARACTER(iotk_attlenx)               :: attr
! 
ispresent = .FALSE.
ALLOCATE (force_mat_(3,natoms)) 
CALL iotk_scan_dat ( iunit, "forces", force_mat_, IERR = ierr, FOUND = ispresent ) 
IF ( ierr /= 0 ) THEN
  DEALLOCATE ( force_mat_) 
  RETURN 
END IF
IF (.NOT. ispresent ) THEN 
   DEALLOCATE ( force_mat_)
   RETURN
END IF 
!
CALL qes_init_matrix( obj, "forces", 3, natoms, force_mat_) 
DEALLOCATE ( force_mat_)
END SUBROUTINE qexsd_get_forces
!
!-------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_stress ( iunit, obj, ispresent ) 
! 
IMPLICIT NONE 
! 
INTEGER, INTENT (IN)                    :: iunit
TYPE ( matrix_type ),INTENT(OUT)        :: obj
LOGICAL,INTENT ( OUT)                   :: ispresent
! 
INTEGER                                 :: ierr 
REAL(DP)                                :: stress_mat_(3,3)
CHARACTER(iotk_attlenx)                 :: attr
!
ispresent = .FALSE.
CALL iotk_scan_dat( iunit, "stress", stress_mat_, IERR = ierr, FOUND = ispresent ) 
IF ( ierr /= 0 ) RETURN 
IF ( .NOT. ispresent ) RETURN
! 
CALL qes_init_matrix( obj, "stress", 3, 3, stress_mat_) 
END SUBROUTINE qexsd_get_stress 
!--------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_output( iunit, obj, ispresent, obj_tagname) 
!--------------------------------------------------------------------------------------------------
! 
IMPLICIT NONE
! 
INTEGER,INTENT(IN)                            :: iunit
CHARACTER (LEN = 256 ),OPTIONAL,INTENT(IN)      :: obj_tagname 
TYPE (output_type ),INTENT(OUT)               :: obj 
LOGICAL,INTENT(OUT)                           :: ispresent 
! 
INTEGER                                       :: ierr 
LOGICAL                                       :: found 
CHARACTER(iotk_attlenx)                      :: attr
! 
IF ( PRESENT ( obj_tagname ) ) THEN
   obj%tagname = TRIM( obj_tagname )
ELSE 
   obj%tagname = "output"
END IF
ispresent = .FALSE.
CALL iotk_scan_begin(iunit, "output", IERR = ierr, FOUND = ispresent ) 
IF ( ierr /= 0) RETURN 
IF (.NOT. ispresent ) RETURN
!
CALL qexsd_get_convergence_info(iunit, obj%convergence_info, found )
IF ( .NOT. found ) CALL  errore ("qexsd_get_output", "reading xml-output: convergence_info not found", 1)
! 
CALL qexsd_get_algorithmic_info ( iunit, obj%algorithmic_info, found ) 
IF ( .NOT. found ) CALL errore ("qexsd_get_output", "reading xml-output: algorithmic_info not found", 1)
! 
CALL qexsd_get_atomic_species( iunit, obj%atomic_species, found ) 
IF ( .NOT. found ) CALL errore ("qexsd_get_output", "reading xml-output: atomic_species not found", 1) 
! 
CALL qexsd_get_atomic_structure ( iunit, obj%atomic_structure, found ) 
IF (.NOT. found ) CALL errore ( "qexsd_get_output", "reading xml-output: atomic_species not found", 1 ) 
! 
CALL qexsd_get_symmetries( iunit, obj%symmetries, found ) 
IF (.NOT. found ) THEN 
   obj%symmetries_ispresent = .FALSE. 
   CALL infomsg ( "qexsd_get_output", "reading xml-output: symmetries not found") 
END IF 
CALL qexsd_get_basis_set (iunit, obj%basis_set, found ) 
IF ( .NOT. found ) CALL errore ( "qexsd_get_output", "reading xml-output: basis_set not found", 1) 
! 

CALL qexsd_get_dft ( iunit, obj%dft, found ) 
IF ( .NOT. found ) CALL errore ( "qexsd_get_output", "reading xml-output: dft not found", 1)
! 

CALL qexsd_get_magnetization ( iunit, obj%magnetization, found) 
IF ( .NOT. found ) CALL errore ( "qexsd_get_output", "reading xml-output: magnetization not found", 1) 
! 

CALL qexsd_get_total_energy ( iunit, obj%total_energy, found ) 
IF ( .NOT. found ) CALL errore ( "qexsd_get_output", "reading xml-output: total_energy not found", 1)
! 

CALL qexsd_get_band_structure ( iunit, obj%band_structure, found ) 
IF (.NOT. found ) CALL errore ( "qexsd_get_output", "reading xml-output: total_energy not found", 1) 
! 


CALL qexsd_get_forces(iunit, obj%forces, obj%atomic_structure%nat, found ) 
obj%forces_ispresent = found 
! 
CALL qexsd_get_stress( iunit, obj%stress, found ) 
obj%stress_ispresent = found 
! 
CALL qexsd_get_outputElectricField ( iunit, obj%electric_field, found ) 
obj%electric_field_ispresent = found 

CALL iotk_scan_end( iunit, "output",IERR = ierr ) 
END SUBROUTINE qexsd_get_output
!
!

!---------------------------------------------------------------------------------------------------------
SUBROUTINE qexsd_get_input ( iunit, obj, ispresent ) 
! 
IMPLICIT NONE
! 
INTEGER, INTENT (IN)                      :: iunit
TYPE ( input_type ),INTENT (OUT)          :: obj 
LOGICAL, INTENT (OUT)                     :: ispresent
! 
INTEGER                                   :: ierr 
LOGICAL                                   :: found, symm_flag_ispresent, boundary_ispresent, &
                                             ekin_fun_ispresent, ext_forc_ispresent, free_pos_ispresent,&
                                             starting_vel_ispresent, electric_field_ispresent, & 
                                             atomic_constr_ispresent, spin_constr_ispresent  
TYPE ( control_variables_type )           :: controls_obj
TYPE ( atomic_species_type )              :: at_specs_obj
TYPE ( atomic_structure_type )            :: at_struct_obj 
TYPE ( dft_type )                         :: dft_in_obj
TYPE ( spin_type )                        :: spin_in_obj
TYPE ( bands_type )                       :: bands_in_obj
TYPE ( basis_type )                       :: basis_obj
TYPE (k_points_IBZ_type )                 :: kipo_obj
TYPE ( electron_control_type )            :: electrons_obj
TYPE ( ion_control_type )                 :: ions_obj
TYPE ( cell_control_type )                :: cell_obj
TYPE ( symmetry_flags_type )              :: symm_flags_obj
TYPE ( boundary_conditions_type)          :: boundary_obj
TYPE ( ekin_functional_type  )            :: ekin_funct_obj
REAL(DP),ALLOCATABLE                      :: ext_forc_(:,:), starting_vel_(:,:) 
TYPE ( matrix_type )                      :: ext_forc_obj, starting_vel_obj
INTEGER, ALLOCATABLE                      :: free_pos_(:,:)
TYPE ( integerMatrix_type )               :: free_pos_obj
TYPE ( electric_field_type)               :: elecfield_obj
TYPE ( atomic_constraints_type )          :: atomic_constr_obj
TYPE ( spin_constraints_type )            :: spin_constr_obj 
CHARACTER(iotk_attlenx)                   :: attr
! 
ispresent = .FALSE. 
CALL iotk_scan_begin( iunit, "input", IERR = ierr , FOUND = ispresent ) 
IF ( ierr /= 0 ) RETURN 
IF ( .NOT. ispresent ) RETURN 
! 
CALL qexsd_get_control_variables(iunit, controls_obj, found ) 
IF (.NOT. found ) CALL errore ( 'qexsd_get_input',' missing mandatory element "control_variables"' , 1) 
! 
CALL qexsd_get_atomic_species( iunit, at_specs_obj, found )  
IF ( .NOT. found ) CALL errore ( 'qexsd_get_input',' missing mandatory element "atomic_species"' , 1)
! 
CALL qexsd_get_atomic_structure ( iunit, at_struct_obj, found )   
IF ( .NOT. found ) CALL errore ( 'qexsd_get_input',' missing mandatory element "atomic_structure"' , 1)
! 
CALL qexsd_get_dft ( iunit, dft_in_obj, found ) 
IF ( .NOT. found ) CALL errore ( 'qexsd_get_input',' missing mandatory element "dft"' , 1)
! 
CALL qexsd_get_spin( iunit, spin_in_obj, found ) 
IF ( .NOT. found ) CALL errore ( 'qexsd_get_input',' missing mandatory element "spin"' , 1)
! 
CALL qexsd_get_bands( iunit, bands_in_obj, found ) 
IF ( .NOT. found ) CALL errore ( 'qexsd_get_input',' missing mandatory element "bands"' , 1)
! 
CALL qexsd_get_basis ( iunit, basis_obj, found ) 
IF ( .NOT. found ) CALL errore ( 'qexsd_get_input',' missing mandatory element "basis"' , 1)
! 
CALL qexsd_get_electron_control( iunit, electrons_obj, found ) 
IF ( .NOT. found ) CALL errore ( 'qexsd_get_input',' missing mandatory element "electron_control"' , 1)
! 
CALL qexsd_get_k_points_IBZ ( iunit, kipo_obj, found ) 
IF ( .NOT. found ) CALL errore ( 'qexsd_get_input',' missing mandatory element "k_points_IBZ"' , 1)
! 
CALL qexsd_get_ion_control( iunit, ions_obj, found ) 
IF ( .NOT. found ) CALL errore ( 'qexsd_get_input',' missing mandatory element "ion_controls"' , 1)
! 
CALL qexsd_get_cell_control( iunit, cell_obj, found ) 
IF ( .NOT. found ) CALL errore ( 'qexsd_get_input',' missing mandatory element "cell_controls"' , 1)
!                  
CALL qexsd_get_symmetry_flags ( iunit, symm_flags_obj, symm_flag_ispresent ) 
! 
CALL qexsd_get_boundary_conditions ( iunit, boundary_obj, boundary_ispresent ) 
!
CALL qexsd_get_ekin_functional ( iunit, ekin_funct_obj, ekin_fun_ispresent ) 
!
ALLOCATE ( ext_forc_(3, at_struct_obj%nat))

CALL iotk_scan_dat( iunit, "external_forces", ext_forc_, IERR = ierr, FOUND = ext_forc_ispresent )
IF (ierr /= 0 ) RETURN 
IF ( ext_forc_ispresent ) CALL qes_init_matrix( ext_forc_obj, "external_forces", 3, at_struct_obj%nat, &
                                               ext_forc_ ) 
DEALLOCATE (ext_forc_)
! 
ALLOCATE (free_pos_(3, at_struct_obj%nat)) 
CALL iotk_scan_dat( iunit, "free_positions", free_pos_, IERR = ierr, FOUND = free_pos_ispresent )
IF ( ierr /= 0 ) RETURN 
IF ( free_pos_ispresent ) CALL qes_init_integerMatrix( free_pos_obj, "free_positions", 3, &
                                                       at_struct_obj%nat, free_pos_)   
DEALLOCATE (free_pos_) 
! 
ALLOCATE ( starting_vel_(3, at_struct_obj%nat) )
CALL iotk_scan_dat ( iunit, "starting_atomic_velocities", starting_vel_, IERR = ierr, &
                     FOUND = starting_vel_ispresent ) 
IF (ierr /= 0 ) RETURN
IF  ( starting_vel_ispresent ) CALL  qes_init_matrix( starting_vel_obj, "starting_atomic_velocities",3, & 
                                                      at_struct_obj%nat, starting_vel_)
DEALLOCATE ( starting_vel_)
! 
CALL qexsd_get_electric_Field (iunit, elecfield_obj, electric_field_ispresent ) 
!  
CALL qexsd_get_atomic_constraints( iunit, atomic_constr_obj, atomic_constr_ispresent ) 
! 
CALL qexsd_get_spin_constraints ( iunit, spin_constr_obj, spin_constr_ispresent ) 
CALL  iotk_scan_end( iunit, "input", IERR = ierr ) 
IF ( ierr /= 0 ) CALL errore ( "qexsd_get_input", "failed to close input element on reading: check file",1)
! 
CALL qes_init_input(obj, "input", controls_obj, at_specs_obj, at_struct_obj, dft_in_obj, spin_in_obj,& 
                    bands_in_obj, basis_obj, electrons_obj, kipo_obj, ions_obj, cell_obj, & 
                    symm_flag_ispresent, symm_flags_obj, boundary_ispresent, boundary_obj,    & 
                    ekin_fun_ispresent, ekin_funct_obj, ext_forc_ispresent, ext_forc_obj, free_pos_ispresent,&
                    free_pos_obj, starting_vel_ispresent, starting_vel_obj, electric_field_ispresent,& 
                    elecfield_obj, atomic_constr_ispresent, atomic_constr_obj, spin_constr_ispresent, &
                    spin_constr_obj )

CALL qes_reset_control_variables( controls_obj) 
CALL qes_reset_atomic_species ( at_specs_obj ) 
CALL qes_reset_atomic_structure ( at_struct_obj ) 
CALL qes_reset_dft ( dft_in_obj ) 
CALL qes_reset_spin ( spin_in_obj ) 
CALL qes_reset_bands ( bands_in_obj ) 
CALL qes_reset_basis ( basis_obj) 
CALL qes_reset_electron_control ( electrons_obj) 
CALL qes_reset_k_points_IBZ ( kipo_obj ) 
CALL qes_reset_ion_control ( ions_obj ) 
CALL qes_reset_cell_control ( cell_obj )
IF ( symm_flag_ispresent) CALL qes_reset_symmetry_flags( symm_flags_obj ) 
IF (  boundary_ispresent ) CALL qes_reset_boundary_conditions ( boundary_obj ) 
IF ( ekin_fun_ispresent ) CALL qes_reset_ekin_functional ( ekin_funct_obj ) 
IF ( ext_forc_ispresent ) CALL qes_reset_matrix ( ext_forc_obj ) 
IF ( free_pos_ispresent ) CALL qes_reset_integerMatrix( free_pos_obj ) 
IF ( starting_vel_ispresent )   CALL qes_reset_matrix( starting_vel_obj ) 
IF ( electric_field_ispresent ) CALL qes_reset_electric_Field( elecfield_obj ) 
IF ( atomic_constr_ispresent )  CALL qes_reset_atomic_constraints ( atomic_constr_obj ) 
IF ( spin_constr_ispresent )   CALL qes_reset_spin_constraints ( spin_constr_obj) 

END SUBROUTINE qexsd_get_input  
!---------------------------------------------------------------------------------------------------
END MODULE
