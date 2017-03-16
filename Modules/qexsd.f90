! Copyright (C) 2003-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE qexsd_module
  !----------------------------------------------------------------------------
  !
  ! This module contains some common subroutines used to read and write
  ! in XML format the data produced by Quantum ESPRESSO package.
  !
  ! Written by Giovanni Borghi, A. Ferretti, ... (2015).
  !
  ! Based on the qexml.f90 routine:
  ! Written by Andrea Ferretti (2006).
  ! Initial work by Carlo Sbraccia (xml_io_base.f90)
  ! Modified by Simone Ziraldo (2013).
  !
  USE kinds,            ONLY : DP
  USE input_parameters, ONLY : input_xml_schema_file, title
  USE mp_world,         ONLY : nproc
  USE mp_images,         ONLY : nimage,nproc_image
  USE mp_pools,         ONLY : npool
  USE mp_bands,         ONLY : ntask_groups, nproc_bgrp, nbgrp
  USE global_version,   ONLY:  version_number, svn_revision
  !
  USE iotk_base,        ONLY : iotk_indent, iotk_maxindent
  USE constants,        ONLY : e2
  USE iotk_module
  USE qes_module
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  ! definitions for the fmt
  !
  CHARACTER(5), PARAMETER :: fmt_name = "QEXSD"
  CHARACTER(5), PARAMETER :: fmt_version = "0.1.0"
  !
  ! some default for kinds
  !
  !INTEGER,   PARAMETER :: DP = selected_real_kind( 14, 200 )
  !
  ! internal data to be set
  !
  CHARACTER(256)   :: datadir_in, datadir_out
  INTEGER          :: iunit, ounit
  !
  ! vars to manage back compatibility
  !
  CHARACTER(10)    :: qexsd_current_version = " "
  CHARACTER(10)    :: qexsd_default_version = trim( fmt_version  )
  LOGICAL          :: qexsd_current_version_init = .FALSE.
  !
  LOGICAL          :: qexsd_use_large_indent = .FALSE.
  !
  CHARACTER(iotk_attlenx) :: attr
  ! 
  TYPE (input_type)                :: qexsd_input_obj
  TYPE (general_info_type)         :: general_info
  TYPE (parallel_info_type)        :: parallel_info
  TYPE (berryPhaseOutput_type)     :: qexsd_bp_obj
  TYPE (dipoleOutput_type )        :: qexsd_dipol_obj
  TYPE ( step_type),ALLOCATABLE    :: steps(:)
  TYPE ( status_type )             :: exit_status
  TYPE ( closed_type )             :: qexsd_closed_element
  INTEGER                          :: step_counter
  !
  ! end of declarations
  !
  PUBLIC :: qexsd_current_version, qexsd_default_version
  PUBLIC :: qexsd_current_version_init
  !
  PUBLIC :: qexsd_input_obj
  ! 
  PUBLIC :: qexsd_init_schema,  qexsd_openschema, qexsd_closeschema
  !
  PUBLIC :: qexsd_init_convergence_info, qexsd_init_algorithmic_info, &
            qexsd_init_atomic_species, qexsd_init_atomic_structure, &
            qexsd_init_symmetries, qexsd_init_basis_set, qexsd_init_dft, &
            qexsd_init_magnetization, qexsd_init_band_structure, & 
            qexsd_init_total_energy, qexsd_init_forces, qexsd_init_stress, &
            qexsd_init_dipole_info, qexsd_init_outputElectricField
  !
  PUBLIC :: qexsd_step_addstep, qexsd_set_status    
  ! 
  PUBLIC :: qexsd_init_berryPhaseOutput, qexsd_bp_obj, qexsd_dipol_obj
CONTAINS
!
!-------------------------------------------
! ... basic (public) subroutines
!-------------------------------------------
!
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_schema( unit_in, unit_out )
      !------------------------------------------------------------------------
      !
      ! just init module data
      !
      IMPLICIT NONE
      INTEGER,                INTENT(in) :: unit_in
      INTEGER,      OPTIONAL, INTENT(in) :: unit_out
      !
      iunit       = unit_in
      ounit       = unit_in
      IF ( present( unit_out ) ) ounit  = unit_out
      !
      ! increase the xml indentation for better reading
      IF ( qexsd_use_large_indent ) THEN
          iotk_indent = 8
          iotk_maxindent = 64
      ENDIF
      !
    END SUBROUTINE qexsd_init_schema
    !
    !
    !------------------------------------------------------------------------
    FUNCTION check_file_exst( filename )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
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
    SUBROUTINE qexsd_openschema( filename )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      CHARACTER(iotk_attlenx)  :: attr
      CHARACTER(len=*), INTENT(IN) :: filename
      CHARACTER(len=16) :: subname = 'qexsd_openschema'
      INTEGER :: ierr, len_steps, i_step
      !
      ! we need a qes-version number here
      CALL iotk_write_attr(attr,"xmlns:qes","http://www.quantum-espresso.org/ns/qes/qes-1.0",FIRST=.true.) 
      !
      ! we need a schema location number (put a string in fortran-input)?
      CALL iotk_write_attr(attr,"xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance",newline=.true.) 
      !
      CALL iotk_write_attr(attr,"xsi:schemaLocation","http://www.quantum-espresso.org/ns/qes/qes-1.0 espresso.xsd",&
                           newline=.true.)
      !
      CALL iotk_open_write(ounit,FILE=filename , root="qes:espresso",attr=attr,binary=.false., &
                           skip_head=.true.,IERR=ierr)
      !
      IF (ierr /= 0) call errore(subname, 'opening xml output file', ierr)
      ! the input file is mandatory to have a validating schema 
      ! here an error should be issued, instead
      !
      CALL qexsd_init_general_info(general_info)
      CALL qes_write_general_info(ounit,general_info)
      CALL qes_reset_general_info(general_info)
      !
      CALL qexsd_init_parallel_info(parallel_info)
      CALL qes_write_parallel_info(ounit,parallel_info)
      CALL qes_reset_parallel_info(parallel_info) 
      IF ( check_file_exst(input_xml_schema_file) )  THEN
         CALL qexsd_cp_line_by_line(ounit,input_xml_schema_file, spec_tag="input")
      ELSE IF ( TRIM(qexsd_input_obj%tagname) == "input") THEN 
         CALL qes_write_input(ounit,qexsd_input_obj)
      END IF
      ! 
      !CALL qes_reset_input(qexsd_input_obj)     
      IF (ALLOCATED(steps) ) THEN 
         len_steps= step_counter 
         DO i_step = 1, len_steps
            CALL qes_write_step(ounit, steps(i_step) )
         END DO 
      END IF
      !
      
    END SUBROUTINE qexsd_openschema
    !
    !
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_general_info(obj)
    !---------------------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE( general_info_type )         ::  obj
      CHARACTER(LEN=*),PARAMETER        ::  TAGNAME="general_info"
      TYPE( creator_type )              ::  creator_obj
      TYPE( created_type )              ::  created_obj
      TYPE( xml_format_type)            ::  xml_fmt_obj
      CHARACTER(LEN=256)                ::  version
      CHARACTER(9)                      ::  cdate, ctime
      CHARACTER(60)                     ::  timestamp
      !
      version=TRIM(version_number)
      IF (svn_revision .NE. "unknown") version=TRIM(version)//&
                                       " (svn rev. " // TRIM (svn_revision) // ")"
      CALL qes_init_creator(creator_obj, "creator", "PWSCF", version, "XML file generated by PWSCF")
      !
      CALL date_and_tim(cdate, ctime) 
      timestamp = 'This run was terminated on:  ' // ctime // ' ' // cdate(1:2) // & 
                  ' '//cdate(3:5) // ' '// cdate (6:9)
     
      CALL qes_init_created (created_obj, "created", cdate, ctime, timestamp ) 
      !
      CALL qes_init_xml_format(xml_fmt_obj, "xml_format", fmt_name, fmt_version, fmt_name//"_"//fmt_version)
      !
      CALL qes_init_general_info( obj, TAGNAME, xml_fmt_obj, creator = creator_obj, created = created_obj,&
                                  job=title)
      !
      CALL qes_reset_creator(creator_obj)
      CALL qes_reset_created(created_obj)
      CALL qes_reset_xml_format(xml_fmt_obj) 
    END SUBROUTINE qexsd_init_general_info    
    !
    !---------------------------------------------------------------------------------------------
    SUBROUTINE   qexsd_init_parallel_info(obj)
    !---------------------------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE ( parallel_info_type )           :: obj
      !
      INTEGER                               :: nthreads=1
#ifdef  __OMP 
      INTEGER,EXTERNAL                      :: omp_get_max
      !     
      nthreads = omp_get_max()
#endif      
      CALL qes_init_parallel_info(obj, "parallel_info", nproc, nthreads, ntask_groups, &
                                  nbgrp, npool, nproc_bgrp)
    END SUBROUTINE qexsd_init_parallel_info
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_closeschema()
      !------------------------------------------------------------------------
      USE mytime,    ONLY: nclock, clock_label
      IMPLICIT NONE
      REAL(DP),EXTERNAL    :: get_clock
      !
      CHARACTER(len=17) :: subname = 'qexsd_closeschema'
      INTEGER :: ierr
      !
      IF (TRIM(exit_status%tagname) == "status") THEN 
         CALL qes_write_status(ounit, exit_status) 
         CALL qexsd_set_closed()
         CALL iotk_write_begin(ounit, "cputime", attr="",new_line=.FALSE.)
      !
         WRITE(ounit, '(I0)', advance = 'no' )  nint(get_clock('PWSCF'))
         CALL iotk_write_end(ounit, "cputime",indentation=.FALSE.)
         CALL qes_write_closed(ounit, qexsd_closed_element)
      END IF
         CALL iotk_close_write(ounit, IERR=ierr)
      !
      CALL errore(subname, 'closing xml input file', ierr)
      !
    END SUBROUTINE qexsd_closeschema
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

      call iotk_free_unit(iun)
      !
      INQUIRE(FILE=trim(filename), EXIST=exists)
      !
      IF(.not.exists) THEN
         CALL errore('qexsd_cp_line_by_line', 'input xml file "' // & 
        &             TRIM(filename) // '" not found', 1)
      ENDIF
      !
      open(iun,FILE=trim(filename),status="old", IOSTAT=ierr)
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
! ... write subroutines
!-------------------------------------------
!
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_convergence_info(obj, n_scf_steps, scf_error, &
                                           opt_conv_ispresent, n_opt_steps, grad_norm )
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(convergence_info_type)   :: obj
      INTEGER,           INTENT(IN) :: n_scf_steps
      REAL(DP),          INTENT(IN) :: scf_error
      LOGICAL,           INTENT(IN) :: opt_conv_ispresent
      INTEGER, OPTIONAL, INTENT(in) :: n_opt_steps
      REAL(DP),OPTIONAL, INTENT(IN) :: grad_norm
      !
      CHARACTER(27)       :: subname="qexsd_init_convergence_info"
      TYPE(scf_conv_type) :: scf_conv
      TYPE(opt_conv_type) :: opt_conv
      !
      call qes_init_scf_conv(scf_conv, "scf_conv", n_scf_steps, scf_error)
      !
      IF ( opt_conv_ispresent ) THEN
          !
          IF ( .NOT. PRESENT(n_opt_steps) ) CALL errore(subname,"n_opt_steps not present",10)
          IF ( .NOT. PRESENT(grad_norm) )   CALL errore(subname,"grad_norm not present",10)
          !
          call qes_init_opt_conv(opt_conv, "opt_conv", n_opt_steps, grad_norm)
      ENDIF
      !
      call qes_init_convergence_info(obj, "convergence_info", scf_conv, opt_conv_ispresent, opt_conv)
      !
      call qes_reset_scf_conv(scf_conv)
      call qes_reset_opt_conv(opt_conv)
      !
    END SUBROUTINE qexsd_init_convergence_info
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_algorithmic_info(obj, real_space_q, uspp, paw )
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(algorithmic_info_type)   :: obj
      LOGICAL,           INTENT(IN) :: real_space_q, uspp, paw
      !
      CALL qes_init_algorithmic_info(obj, "algorithmic_info", real_space_q, uspp, paw)
      !
    END SUBROUTINE qexsd_init_algorithmic_info
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_atomic_species(obj, nsp, atm, psfile, amass, starting_magnetization,&
                                         angle1,angle2) 
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(atomic_species_type)    :: obj
      INTEGER,          INTENT(IN) :: nsp
      CHARACTER(len=*), INTENT(IN) :: atm(:)
      CHARACTER(len=*), INTENT(IN) :: psfile(:)
      REAL(DP), OPTIONAL, INTENT(IN) :: amass(:)
      REAL(DP), OPTIONAL, INTENT(IN) :: starting_magnetization(:)
      REAL(DP), OPTIONAL, INTENT(IN) :: angle1(:),angle2(:)
      !
      TYPE(species_type), ALLOCATABLE :: species(:)
      REAL(DP)  :: amass_ = 0.0d0
      REAL(DP)  :: start_mag_ = 0.0d0
      REAL(DP)  :: spin_teta = 0.0d0
      REAL(DP)  :: spin_phi  = 0.0d0
      INTEGER   :: i
      
      ALLOCATE(species(nsp))
      !
      DO i = 1, nsp
          !
          IF ( PRESENT(amass) ) amass_=amass(i)
          IF ( PRESENT(starting_magnetization) ) start_mag_=starting_magnetization(i)
          IF ( PRESENT( angle1 ) )  spin_teta =angle1(i)
          IF ( PRESENT( angle2 ) )  spin_phi = angle2(i)
          !
          CALL qes_init_species( species(i), "species", TRIM(atm(i)),PRESENT(amass),amass_, &
                               TRIM(psfile(i)), PRESENT(starting_magnetization), start_mag_,&
                               PRESENT(angle1),spin_teta,PRESENT(angle2),spin_phi)
      ENDDO
      !
      CALL qes_init_atomic_species(obj, "atomic_species", nsp, SIZE(species), species)
      !
      DO i = 1, nsp
          CALL qes_reset_species(species(i))
      ENDDO
      DEALLOCATE(species)
      !
    END SUBROUTINE qexsd_init_atomic_species
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_atomic_structure(obj, nsp, atm, ityp, nat, tau, &
                                           alat, a1, a2, a3, ibrav)
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(atomic_structure_type)  :: obj
      INTEGER,          INTENT(IN) :: nsp, nat
      INTEGER,          INTENT(in) :: ityp(:)
      CHARACTER(LEN=*), INTENT(in) :: atm(:)
      REAL(DP),         INTENT(IN) :: tau(3,*)! cartesian atomic positions, a.u.
      REAL(DP),         INTENT(IN) :: alat
      REAL(DP),         INTENT(IN) :: a1(:), a2(:), a3(:)
      INTEGER,          INTENT(IN) :: ibrav
      !
      INTEGER         :: ia
      TYPE(atom_type), ALLOCATABLE :: atom(:)
      TYPE(cell_type) :: cell
      TYPE(atomic_positions_type)  :: atomic_pos
      TYPE(wyckoff_positions_type) :: wyckoff_pos
      REAL(DP)                     :: new_alat
      LOGICAL                      :: ibrav_ispresent
      !
      ! atomic positions
      !
      IF ( ibrav .gt. 0 ) THEN 
         ibrav_ispresent = .TRUE.
      ELSE
         ibrav_ispresent = .FALSE.
      END IF
      !
      ALLOCATE(atom(nat))
      DO ia = 1, nat
          CALL qes_init_atom( atom(ia), "atom", name=trim(atm(ityp(ia))), &
                             position="", position_ispresent=.FALSE., &
                             atom=tau(1:3,ia), index_ispresent = .TRUE.,&
                             index = ia )
      ENDDO
      !
      CALL qes_init_atomic_positions(atomic_pos, "atomic_positions", SIZE(atom), atom)
      !
      DO ia = 1, nat
          CALL qes_reset_atom( atom(ia) )
      ENDDO
      DEALLOCATE(atom)
      !
      ! cell
      !
      CALL qes_init_cell(cell, "cell", a1, a2, a3)
      !
      ! global init
      !
      CALL qes_init_atomic_structure(obj, "atomic_structure", nat=nat, &
                     alat=alat, alat_ispresent=.TRUE., atomic_positions_ispresent=.TRUE., &
                     atomic_positions=atomic_pos, wyckoff_positions_ispresent=.FALSE., &
                     wyckoff_positions=wyckoff_pos, cell=cell ,& 
                     bravais_index_ispresent = ibrav_ispresent, bravais_index=ibrav)
      ! 
      ! cleanup 
      ! 
      CALL qes_reset_atomic_positions(atomic_pos)
      CALL qes_reset_cell(cell)
      !
    END SUBROUTINE qexsd_init_atomic_structure
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_symmetries(obj, nsym, nrot, space_group, s, ft, sname, t_rev, nat, irt, &
                                     class_names, verbosity, noncolin)
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(symmetries_type)    :: obj
      INTEGER,          INTENT(IN) :: nsym, nrot, nat
      INTEGER,          INTENT(IN) :: space_group
      INTEGER,          INTENT(IN) :: s(:,:,:), irt(:,:)
      REAL(DP),         INTENT(IN) :: ft(:,:)
      INTEGER,          INTENT(IN) :: t_rev(:)
      CHARACTER(LEN=*), INTENT(IN) :: sname(:), verbosity
      CHARACTER(LEN=15),INTENT(IN) :: class_names(:)
      LOGICAL,INTENT(IN)           :: noncolin
      !
      TYPE(symmetry_type), ALLOCATABLE  :: symm(:)
      TYPE(equivalent_atoms_type)  :: equiv_atm
      TYPE(info_type)              :: info
      TYPE(matrix_type)            :: matrix
      CHARACTER(LEN=15)            :: classname
      CHARACTER(LEN=256)            :: la_info
      LOGICAL                      :: class_ispresent = .FALSE., time_reversal_ispresent = .FALSE.
      INTEGER                      :: i
      
      ALLOCATE(symm(nrot))
      !
      IF ( TRIM(verbosity) .EQ. 'high' .OR. TRIM(verbosity) .EQ. 'medium')  class_ispresent= .TRUE.
      IF ( noncolin  ) time_reversal_ispresent = .TRUE.
      DO i = 1, nrot
          !
          classname = class_names(i)
          IF ( i .LE. nsym ) THEN 
             la_info = "crystal_symmetry"
          ELSE 
             la_info = "lattice_symmetry"
          END IF
          CALL qes_init_info(info, "info", name=sname(i), name_ispresent=.TRUE., &
                             class=classname, class_ispresent = class_ispresent,   &
                             time_reversal=(t_rev(i)==1), time_reversal_ispresent = time_reversal_ispresent, &
                             INFO= TRIM(la_info) )
          !
          CALL qes_init_matrix(matrix, "rotation", ndim1_mat=3, ndim2_mat=3, mat=real(s(:,:,i),DP))
          !
          IF ( i .LE. nsym ) THEN 
             CALL qes_init_equivalent_atoms(equiv_atm, "equivalent_atoms", nat=nat, ndim_index_list=nat, &
                                         index_list=irt(i,1:nat)  )
          !
             CALL qes_init_symmetry(symm(i),"symmetry", info=info, rotation=matrix, &
                                 fractional_translation_ispresent=.TRUE., fractional_translation=ft(:,i), &
                                 equivalent_atoms_ispresent=.TRUE., equivalent_atoms=equiv_atm)
          ELSE 
             CALL qes_init_symmetry ( symm(i), "symmetry", INFO = info, ROTATION = matrix, &
                                      FRACTIONAL_TRANSLATION_ISPRESENT = .FALSE., FRACTIONAL_TRANSLATION=ft(:,i), &
                                      EQUIVALENT_ATOMS_ISPRESENT = .FALSE.,  EQUIVALENT_ATOMS=equiv_atm) 
          END IF
          !
          CALL qes_reset_info(info)
          CALL qes_reset_matrix(matrix)
          IF ( i .LT. nsym ) THEN 
             CALL qes_reset_equivalent_atoms( equiv_atm )
          ELSE IF ( i .EQ. nrot ) THEN  
            CALL qes_reset_equivalent_atoms( equiv_atm )
          END IF
          !
      ENDDO
      !
      CALL qes_init_symmetries(obj,"symmetries",NSYM = nsym, NROT=nrot, SPACE_GROUP = space_group, &
                               NDIM_SYMMETRY=SIZE(symm), SYMMETRY=symm )
      !
      DO i = 1, nsym
         CALL qes_reset_symmetry(symm(i))
      ENDDO
      DEALLOCATE(symm)
      !
    END SUBROUTINE qexsd_init_symmetries
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_basis_set(obj, gamma_only, ecutwfc, ecutrho, &
                                    nr1, nr2, nr3, nr1s, nr2s, nr3s, &
                                    fft_box_ispresent, nr1b, nr2b, nr3b, &
                                    ngm, ngms, npwx, b1, b2, b3 )
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(basis_set_type)    :: obj
      LOGICAL,          INTENT(IN) :: gamma_only
      INTEGER,          INTENT(IN) :: nr1, nr2, nr3
      INTEGER,          INTENT(IN) :: nr1s, nr2s, nr3s
      LOGICAL,          INTENT(IN) :: fft_box_ispresent
      INTEGER,          INTENT(IN) :: nr1b, nr2b, nr3b
      INTEGER,          INTENT(IN) :: ngm, ngms, npwx
      REAL(DP),         INTENT(IN) :: ecutwfc, ecutrho
      REAL(DP),         INTENT(IN) :: b1(3), b2(3), b3(3)
      !
      TYPE(basisSetItem_type) :: fft_grid
      TYPE(basisSetItem_type) :: fft_smooth
      TYPE(basisSetItem_type) :: fft_box
      TYPE(reciprocal_lattice_type) :: recipr_latt

      CALL qes_init_basisSetItem(fft_grid, "fft_grid", nr1, nr2, nr3, "")
      CALL qes_init_basisSetItem(fft_smooth, "fft_smooth", nr1s, nr2s, nr3s, "")
      CALL qes_init_basisSetItem(fft_box, "fft_box", nr1b, nr2b, nr3b, "" )
      CALL qes_init_reciprocal_lattice(recipr_latt, "reciprocal_lattice", b1, b2, b3)

      CALL qes_init_basis_set(obj, "basis_set", GAMMA_ONLY_ISPRESENT=.TRUE., GAMMA_ONLY=gamma_only, &
                              ECUTWFC=ecutwfc, ECUTRHO_ISPRESENT=.TRUE., ECUTRHO=ecutrho, FFT_GRID=fft_grid, &
                              FFT_SMOOTH_ISPRESENT=.TRUE., FFT_SMOOTH=fft_smooth, &
                              FFT_BOX_ISPRESENT=fft_box_ispresent, FFT_BOX=fft_box, NGM=ngm, &
                              NGMS_ISPRESENT=.TRUE., NGMS=ngms, NPWX=npwx, RECIPROCAL_LATTICE=recipr_latt )
      !
      CALL qes_reset_basisSetItem(fft_grid)
      CALL qes_reset_basisSetItem(fft_smooth)
      CALL qes_reset_basisSetItem(fft_box)
      CALL qes_reset_reciprocal_lattice(recipr_latt)
      !
    END SUBROUTINE qexsd_init_basis_set
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_dft(obj, functional, root_is_output, dft_is_hybrid, nqx1, nqx2, nqx3, ecutfock,       &
                   exx_fraction, screening_parameter, exxdiv_treatment, x_gamma_extrapolation, ecutvcut,        &
                   dft_is_lda_plus_U, lda_plus_U_kind, llmax, noncolin, nspin, nsp, ldim, nat, species, ityp,   &
                   Hubbard_U, Hubbard_J0, Hubbard_alpha, Hubbard_beta, Hubbard_J, starting_ns, Hubbard_ns,      &
                   Hubbard_ns_nc, U_projection_type, dft_is_vdW, vdw_corr, nonlocal_term, london_s6, london_c6, &
                   london_rcut, xdm_a1, xdm_a2 ,ts_vdw_econv_thr, ts_vdw_isolated, is_hubbard, psd)
      !------------------------------------------------------------------------
      USE  parameters,           ONLY:  lqmax
      USE  input_parameters,     ONLY:  nspinx
      IMPLICIT NONE
      !
      TYPE(dft_type)    :: obj
      CHARACTER(len=*), INTENT(IN) :: functional, nonlocal_term 
      LOGICAL,          INTENT(IN) :: dft_is_hybrid
      LOGICAL,          INTENT(IN) :: root_is_output
      INTEGER,          INTENT(IN) :: nqx1, nqx2, nqx3
      REAL(DP),         INTENT(IN) :: ecutfock
      REAL(DP),         INTENT(IN) :: exx_fraction
      REAL(DP),         INTENT(IN) :: screening_parameter
      CHARACTER(len=*), INTENT(IN) :: exxdiv_treatment
      LOGICAL,          INTENT(IN) :: x_gamma_extrapolation
      REAL(DP),         INTENT(IN) :: ecutvcut
      !
      LOGICAL,          INTENT(IN) :: dft_is_lda_plus_U, noncolin 
      INTEGER,          INTENT(IN) :: lda_plus_U_kind
      INTEGER,          INTENT(IN) :: llmax, nspin, nsp, ldim, nat
      CHARACTER(len=*), INTENT(IN) :: species(nsp)
      INTEGER,          INTENT(IN) :: ityp(nat)
      REAL(DP),         INTENT(IN) :: Hubbard_U(nsp)
      REAL(DP),         INTENT(IN) :: Hubbard_J0(nsp)
      REAL(DP),         INTENT(IN) :: Hubbard_alpha(nsp)
      REAL(DP),         INTENT(IN) :: Hubbard_beta(nsp)
      REAL(DP),         INTENT(IN) :: Hubbard_J(3,nsp)
      REAL(DP),         INTENT(IN) :: starting_ns(lqmax,nspinx,nsp)
      REAL(DP),         INTENT(IN) :: Hubbard_ns(ldim,ldim,nspin,nat)
      COMPLEX(DP),      INTENT(IN) :: Hubbard_ns_nc(ldim,ldim,nspin,nat)
      CHARACTER(len=*), INTENT(IN) :: U_projection_type
      LOGICAL,INTENT(IN)           :: is_hubbard(nsp)
      CHARACTER(LEN=2),INTENT(IN)  :: psd(nsp)
      !
      LOGICAL,          INTENT(IN) :: dft_is_vdW, ts_vdw_isolated
      CHARACTER(len=*), INTENT(IN) :: vdw_corr
      REAL(DP),         INTENT(IN) :: london_s6
      REAL(DP),         INTENT(IN) :: london_rcut
      REAL(DP),         INTENT(IN) :: xdm_a1
      REAL(DP),         INTENT(IN) :: xdm_a2
      REAL(DP),         INTENT(IN) :: london_c6(nsp), ts_vdw_econv_thr
      !
      INTEGER  :: i, is, isp, ind,hubb_l,hubb_n
      TYPE(hybrid_type) :: hybrid
      TYPE(qpoint_grid_type) :: qpoint_grid
      TYPE(dftU_type) :: dftU
      TYPE(vdW_type) :: vdW
      TYPE(HubbardCommon_type), ALLOCATABLE :: Hubbard_U_(:)
      TYPE(HubbardCommon_type), ALLOCATABLE :: Hubbard_J0_(:)
      TYPE(HubbardCommon_type), ALLOCATABLE :: Hubbard_alpha_(:)
      TYPE(HubbardCommon_type), ALLOCATABLE :: Hubbard_beta_(:)
      TYPE(HubbardJ_type),      ALLOCATABLE :: Hubbard_J_(:)
      TYPE(starting_ns_type),   ALLOCATABLE :: starting_ns_(:)
      TYPE(Hubbard_ns_type),    ALLOCATABLE :: Hubbard_ns_(:)
      TYPE(HubbardCommon_type), ALLOCATABLE :: london_c6_obj(:)
      REAL(DP),                 ALLOCATABLE :: Hubb_occ_aux(:,:) 
      INTEGER                               :: m1, m2
      LOGICAL  :: Hubbard_U_ispresent
      LOGICAL  :: Hubbard_J0_ispresent
      LOGICAL  :: Hubbard_alpha_ispresent
      LOGICAL  :: Hubbard_beta_ispresent
      LOGICAL  :: Hubbard_J_ispresent
      LOGICAL  :: starting_ns_ispresent
      LOGICAL  :: Hubbard_ns_ispresent
      LOGICAL  :: london_c6_ispresent, london_s6_ispresent, london_rvdw_ispresent, ts_vdw_econv_thr_ispresent, & 
                  london_rcut_ispresent, ts_vdw_isolated_ispresent, xdm_a1_ispresent, xdm_a2_ispresent, &
                  empirical_vdw = .FALSE. 
      INTEGER  :: ndim_london_c6                   
      CHARACTER(10), ALLOCATABLE :: label(:)
      CHARACTER                  :: hubbard_shell 
      INTEGER,EXTERNAL           :: set_hubbard_l,set_hubbard_n
      !
      !
      IF ( dft_is_hybrid ) THEN
          !
          CALL qes_init_qpoint_grid(qpoint_grid, "qpoint_grid", nqx1, nqx2, nqx3, "")
          !
          CALL qes_init_hybrid(hybrid, "hybrid", qpoint_grid, ecutfock, exx_fraction, &
                               screening_parameter, exxdiv_treatment, x_gamma_extrapolation, ecutvcut)
          !
          CALL qes_reset_qpoint_grid(qpoint_grid)
          !
      ENDIF
      !
      IF ( dft_is_lda_plus_U ) THEN
          !
          ALLOCATE(label(nsp))
          DO i = 1, nsp
             IF (is_hubbard(i)) THEN
                hubb_l=set_hubbard_l(psd(i))
                hubb_n=set_hubbard_n(psd(i))
                SELECT CASE ( hubb_l ) 
                   CASE ( 0)  
                      hubbard_shell='s'
                   CASE ( 1 ) 
                      hubbard_shell='p'
                   CASE( 2 ) 
                      hubbard_shell='d'
                   CASE( 3 ) 
                      hubbard_shell='f'
                END SELECT
                WRITE (label(i),'(I0,A)') hubb_n,hubbard_shell
             ELSE
                label(i)="no Hubbard"
             END IF
          END DO
          !
          Hubbard_U_ispresent = (SIZE(Hubbard_U)>0)
          Hubbard_J0_ispresent = (SIZE(Hubbard_J0)>0)
          Hubbard_alpha_ispresent = (SIZE(Hubbard_alpha)>0)
          Hubbard_beta_ispresent = (SIZE(Hubbard_beta)>0)
          Hubbard_J_ispresent = (SIZE(Hubbard_J)>0)
          Hubbard_ns_ispresent = (SIZE(Hubbard_ns)>0)
          starting_ns_ispresent = (SIZE(starting_ns)>0)
          !
          ALLOCATE( Hubbard_U_(nsp) )
          ALLOCATE( Hubbard_J0_(nsp) )
          ALLOCATE( Hubbard_alpha_(nsp) )
          ALLOCATE( Hubbard_beta_(nsp) )
          ALLOCATE( Hubbard_J_(nsp) )
          !
          IF (noncolin ) THEN 
             ALLOCATE (starting_ns_(nsp))
             ALLOCATE (Hubbard_ns_(nat))
          ELSE 
            ALLOCATE( starting_ns_(min(nspin,nspinx)*nsp) )
            ALLOCATE( Hubbard_ns_(nspin*nat) )
          END IF
          !
          DO i = 1, nsp
              CALL qes_init_HubbardCommon(Hubbard_U_(i),"Hubbard_U",TRIM(species(i)),TRIM(label(i)),Hubbard_U(i))
              CALL qes_init_HubbardCommon(Hubbard_J0_(i),"Hubbard_J0",TRIM(species(i)),TRIM(label(i)),Hubbard_J0(i))
              CALL qes_init_HubbardCommon(Hubbard_alpha_(i),"Hubbard_alpha",TRIM(species(i)),TRIM(label(i)),&
                                          Hubbard_alpha(i))
              CALL qes_init_HubbardCommon(Hubbard_beta_(i),"Hubbard_beta",TRIM(species(i)),TRIM(label(i)),&
                                          Hubbard_beta(i))
              CALL qes_init_HubbardJ(Hubbard_J_(i),"Hubbard_J",TRIM(species(i)),TRIM(label(i)),Hubbard_J(1:3,i))
          ENDDO
          !
          ind = 0
          IF (starting_ns_ispresent) THEN 
             IF (noncolin) THEN 
                DO i = 1, nsp 
                   ind = ind + 1 
                   CALL qes_init_starting_ns(starting_ns_(ind), "starting_ns", TRIM (species(i)),TRIM (label(i)),&
                                             1,2*llmax, starting_ns(1:2*llmax, 1, i))
                END DO
             ELSE 
                DO is = 1, MIN(nspin,nspinx) 
                   DO i  = 1, nsp
                      ind = ind+1
                      CALL qes_init_starting_ns(starting_ns_(ind),"starting_ns",TRIM(species(i)),TRIM(label(i)), &
                                                is, llmax, starting_ns(1:llmax,is,i) )
                  ENDDO
                ENDDO
             END IF
          END IF
          !
          ind = 0
          IF (noncolin) THEN 
             ALLOCATE (Hubb_occ_aux(2*ldim,2*ldim))
             DO i = 1, nat 
                Hubb_occ_aux = 0.d0
                DO m1 =1, ldim 
                   DO m2 = 1, ldim 
                      Hubb_occ_aux(     m1,     m2) = SQRT(DCONJG(Hubbard_ns_nc(m1,m2,1,i))*Hubbard_ns_nc(m1,m2,1,i))
                      Hubb_occ_aux(     m1,ldim+m2) = SQRT(DCONJG(Hubbard_ns_nc(m1,m2,2,i))*Hubbard_ns_nc(m1,m2,2,i)) 
                      Hubb_occ_aux(ldim+m1,     m2) = SQRT(DCONJG(Hubbard_ns_nc(m1,m2,3,i))*Hubbard_ns_nc(m1,m2,3,i))
                      Hubb_occ_aux(ldim+m1,ldim+m2) = SQRT(DCONJG(Hubbard_ns_nc(m1,m2,4,i))*Hubbard_ns_nc(m1,m2,4,i))
                   END DO
                END DO
                CALL qes_init_Hubbard_ns(Hubbard_ns_(i),"Hubbard_ns_mod", TRIM(species(ityp(i))),TRIM(label(ityp(i))), &
                                         1, i, 2*ldim, 2*ldim, Hubb_occ_aux(:,:))
             END DO 
             DEALLOCATE ( Hubb_occ_aux) 
          ELSE 
             DO i = 1, nat
                DO is = 1, nspin
                   ind = ind+1
                   CALL qes_init_Hubbard_ns(Hubbard_ns_(ind),"Hubbard_ns", TRIM(species(ityp(i))),TRIM(label(ityp(i))), &
                                       is, i, ldim, ldim, Hubbard_ns(:,:,is,i) )
                ENDDO
             ENDDO
          END IF
          !
          ! main init
          CALL qes_init_dftU(dftU, "dftU", .TRUE., lda_plus_u_kind, &
                              Hubbard_U_ispresent, SIZE(Hubbard_U_), Hubbard_U_, &
                              Hubbard_J0_ispresent, SIZE(Hubbard_J0_), Hubbard_J0_, &
                              Hubbard_alpha_ispresent, SIZE(Hubbard_alpha_), Hubbard_alpha_, &
                              Hubbard_beta_ispresent, SIZE(Hubbard_beta_), Hubbard_beta_, &
                              Hubbard_J_ispresent, SIZE(Hubbard_J_), Hubbard_J_, &
                              starting_ns_ispresent, SIZE(starting_ns_), starting_ns_, &
                              Hubbard_ns_ispresent, SIZE(Hubbard_ns_), Hubbard_ns_, &
                              .TRUE., U_projection_type)
          !
          DO i = 1, nsp
              CALL qes_reset_HubbardCommon(Hubbard_U_(i))
              CALL qes_reset_HubbardCommon(Hubbard_J0_(i))
              CALL qes_reset_HubbardCommon(Hubbard_alpha_(i))
              CALL qes_reset_HubbardCommon(Hubbard_beta_(i))
              CALL qes_reset_HubbardJ(Hubbard_J_(i))
          ENDDO
          !
          DEALLOCATE(Hubbard_U_)
          DEALLOCATE(Hubbard_J0_)
          DEALLOCATE(Hubbard_alpha_)
          DEALLOCATE(Hubbard_beta_)
          DEALLOCATE(Hubbard_J_)
          !
          DO i = 1, SIZE(starting_ns_)
              CALL qes_reset_starting_ns(starting_ns_(i))
          ENDDO
          DEALLOCATE(starting_ns_)
          !
          DO i = 1, SIZE(Hubbard_ns_)
              CALL qes_reset_Hubbard_ns(Hubbard_ns_(i))
          ENDDO
          DEALLOCATE(Hubbard_ns_)
          !
          DEALLOCATE(label)
          !
      ENDIF
      !
      SELECT CASE ( TRIM (vdw_corr )) 
      CASE ( 'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d') 
           empirical_vdw = .TRUE.
           london_s6_ispresent = .TRUE. 
           london_rcut_ispresent = .TRUE. 
           xdm_a1_ispresent = .TRUE. 
           xdm_a2_ispresent = .TRUE. 
           IF ( ANY(london_c6 .GT.  0.d0 )) THEN 
              london_c6_ispresent = .TRUE.
              ndim_london_c6 = 0 
              DO isp = 1, nsp 
                 IF ( london_c6(isp) .GT. 0.d0 ) THEN 
                    ndim_london_c6 = ndim_london_c6 + 1
                 END IF 
              END DO
              ALLOCATE (london_c6_obj(ndim_london_c6))
              ndim_london_c6 = 0 
              DO isp = 1, nsp
                 IF ( london_c6(isp) .GT. 0.d0) THEN
                    ndim_london_c6 = ndim_london_c6 + 1  
                    CALL qes_init_hubbardcommon(london_c6_obj(ndim_london_c6), "london_c6", TRIM(species(isp)),"",&
                                                london_c6(isp))
                 END IF 
              END DO                        
           ELSE 
              london_c6_ispresent = .FALSE. 
              ALLOCATE ( london_c6_obj(1))
           END IF
           ts_vdw_econv_thr_ispresent = .FALSE. 
           ts_vdw_isolated_ispresent = .FALSE. 
      CASE ( 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler')
           empirical_vdw = .TRUE.
           london_s6_ispresent =   .FALSE.
           london_c6_ispresent = .FALSE.
           ALLOCATE ( london_c6_obj(1)) 
           london_rcut_ispresent = .FALSE.
           xdm_a1_ispresent =      .FALSE.
           xdm_a2_ispresent =      .FALSE.
           ts_vdw_econv_thr_ispresent = .TRUE. 
           ts_vdw_isolated_ispresent  = .TRUE. 
      CASE default 
           empirical_vdw = .FALSE.
           ts_vdw_econv_thr_ispresent = .FALSE.
           ts_vdw_isolated_ispresent = .FALSE.
           london_s6_ispresent =   .FALSE.
           london_c6_ispresent = .FALSE.
           ALLOCATE (london_c6_obj(1))
           london_rcut_ispresent = .FALSE.
           xdm_a1_ispresent =      .FALSE.
           xdm_a2_ispresent =      .FALSE.
           london_c6_ispresent   = .FALSE.
      END SELECT 

      IF ( dft_is_vdW .OR. empirical_vdw ) THEN
          !
          CALL qes_init_vdW(vdW, "vdW", TRIM(vdw_corr), root_is_output,  TRIM(nonlocal_term), london_s6_ispresent, london_s6, &
                            ts_vdw_econv_thr_ispresent, ts_vdw_econv_thr, ts_vdw_isolated_ispresent, ts_vdw_isolated,& 
                            london_rcut_ispresent, london_rcut, xdm_a1_ispresent, xdm_a1, xdm_a2_ispresent, xdm_a2, &
                            london_c6_ispresent, ndim_london_c6, london_c6_obj )
          !
          IF (london_c6_ispresent )   THEN
             DO isp=1, ndim_london_c6
                CALL qes_reset_hubbardcommon(london_c6_obj(isp))
             END DO 
          END IF
          DEALLOCATE ( london_c6_obj) 
      ENDIF
        
      CALL qes_init_dft(obj, "dft", functional, dft_is_hybrid, hybrid, &
                             dft_is_lda_plus_U, dftU, (dft_is_vdW .OR. empirical_vdw) , vdW)
      !
      IF (dft_is_hybrid)      CALL qes_reset_hybrid(hybrid)
      IF (dft_is_lda_plus_U)  CALL qes_reset_dftU(dftU)
      IF (dft_is_vdW .OR. empirical_vdw )  CALL qes_reset_vdW(vdW)
      !
    END SUBROUTINE qexsd_init_dft
    !
    !
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_magnetization(obj, lsda, noncolin, spinorbit, total_mag, total_mag_nc, &
                                        absolute_mag, do_magnetization)
      !------------------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(magnetization_type)    :: obj
      LOGICAL,         INTENT(IN) :: lsda, noncolin, spinorbit
      REAL(DP),        INTENT(IN) :: total_mag, absolute_mag
      REAL(DP),        INTENT(IN) :: total_mag_nc(3)
      LOGICAL,         INTENT(IN) :: do_magnetization
      !
      CALL qes_init_magnetization(obj, "magnetization", lsda, noncolin, spinorbit, total_mag, absolute_mag, &
                                 do_magnetization)
      !
    END SUBROUTINE qexsd_init_magnetization 
    !
    ! 
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_band_structure(obj, lsda, noncolin, lspinorb, nbnd, nelec, n_wfc_at, occupations_are_fixed, & 
                                         fermi_energy, two_fermi_energies, ef_updw, et, wg, nks, xk, ngk, wk, & 
                                         starting_kpoints, occupation_kind, smearing, wf_collected)
    !----------------------------------------------------------------------------------------
    IMPLICIT NONE
    !
    TYPE(band_structure_type)               :: obj
    CHARACTER(LEN=*), PARAMETER             :: TAGNAME="band_structure"
    LOGICAL,INTENT(IN)                      :: lsda, noncolin, lspinorb, occupations_are_fixed
    INTEGER,INTENT(IN)                      :: nbnd, nks, n_wfc_at
    REAL(DP),INTENT(IN)                     :: nelec, fermi_energy
    REAL(DP),DIMENSION(:,:),INTENT(IN)      :: et, wg, xk
    REAL(DP),DIMENSION(:),INTENT(IN)        :: wk
    INTEGER,DIMENSION(:),INTENT(IN)         :: ngk      
    REAL(DP),DIMENSION(2),INTENT(IN)        :: ef_updw 
    LOGICAL,INTENT(IN)                      :: two_fermi_energies
    TYPE(k_points_IBZ_type),INTENT(IN)      :: starting_kpoints
    TYPE(occupations_type), INTENT(IN)      :: occupation_kind
    TYPE(smearing_type),OPTIONAL,INTENT(IN) :: smearing
    LOGICAL,INTENT(IN)                      :: wf_collected                    
    ! 
    LOGICAL                                 :: nbnd_up_ispresent, nbnd_dw_ispresent, &
                                               fermi_energy_ispresent, HOL_ispresent, & 
                                               n_wfc_at_ispresent = .TRUE.  
    INTEGER                                 :: nbnd_up,nbnd_dw
    INTEGER                                 :: ndim_ks_energies, nbnd_tot, ik
    TYPE(k_point_type)                      :: kp_obj
    TYPE(ks_energies_type),ALLOCATABLE      :: ks_objs(:)
    TYPE (k_points_IBZ_type)                :: starting_k_points_
    TYPE ( occupations_type)                :: occupations_kind_ 
    REAL(DP),DIMENSION(:),ALLOCATABLE       :: eigenvalues, occupations
    TYPE (smearing_type)                    :: smearing_ 
    !
    !
    ndim_ks_energies=nks   
    nbnd_tot=nbnd
    !
    IF ( lsda ) THEN 
       ndim_ks_energies=ndim_ks_energies/2
       nbnd_up=nbnd
       nbnd_dw=nbnd
       nbnd_tot=nbnd_up+nbnd_dw
       nbnd_up_ispresent=.true.
       nbnd_dw_ispresent=.true.
    ELSE 
       nbnd_up=0
       nbnd_dw=0
       nbnd_up_ispresent=.false.
       nbnd_dw_ispresent=.false. 
    END IF 
    IF (fermi_energy.GT.-1.D6 .AND. ( .NOT. two_fermi_energies ) ) THEN
      IF ( occupations_are_fixed ) THEN 
         fermi_energy_ispresent = .FALSE. 
         HOL_ispresent = .TRUE. 
      ELSE 
         fermi_energy_ispresent = .TRUE.
         HOL_ispresent = .FALSE. 
      END IF 
    ELSE 
      fermi_energy_ispresent=.FALSE.
      HOL_ispresent = .FALSE.
    END IF  
    !
    !   
    ALLOCATE(eigenvalues(nbnd_tot),occupations(nbnd_tot))
    ALLOCATE(ks_objs(ndim_ks_energies))
    !  
    DO ik=1,ndim_ks_energies
       CALL qes_init_k_point(kp_obj,"k_point",wk(ik),.true.,"",.FALSE., xk(:,ik))
       IF ( lsda ) THEN 
          eigenvalues(1:nbnd_up)=et(1:nbnd_up,ik)/e2
          eigenvalues(nbnd_up+1:nbnd_tot)=et(1:nbnd_dw,ndim_ks_energies+ik)/e2
       ELSE 
          eigenvalues(1:nbnd_tot)= et(1:nbnd_tot,ik)/e2
       END IF
       !
       !
       IF (lsda) THEN 
          IF ( ABS(wk(ik)).GT.1.d-10) THEN 
             occupations(1:nbnd_up)=wg(1:nbnd_up,ik)/wk(ik)
             occupations(nbnd_up+1:nbnd_tot)=wg(1:nbnd_dw,ndim_ks_energies+ik)/wk(ndim_ks_energies+ik)
          ELSE 
             occupations(1:nbnd_up)=wg(1:nbnd_up,ik)
             occupations(nbnd_up+1:nbnd_tot)=wg(1:nbnd_dw,ik) 
          END IF            
       ELSE 
          IF (ABS(wk(ik)).GT.1.d-10) THEN
              occupations(1:nbnd_tot)=wg(1:nbnd_tot,ik)/wk(ik)
          ELSE
              occupations(1:nbnd_tot)=wg(1:nbnd_tot,ik)
          END IF
       END IF
       !
       !
       CALL  qes_init_ks_energies(ks_objs(ik),"ks_energies",kp_obj,ngk(ik),nbnd_tot,eigenvalues,& 
                      nbnd_tot,occupations)
       !
       eigenvalues=0.d0
       occupations=0.d0
       CALL qes_reset_k_point(kp_obj)  
    END DO 
    !
    IF ( PRESENT(smearing) ) smearing_ = smearing
!
    starting_k_points_ = starting_kpoints
    starting_k_points_%tagname = "starting_k_points"
!
    occupations_kind_  = occupation_kind
    occupations_kind_%tagname = "occupations_kind"
! 
    CALL qes_init_band_structure( obj,TAGNAME,lsda,noncolin,lspinorb, nbnd , nbnd_up_ispresent,&
                  nbnd_up,nbnd_dw_ispresent,nbnd_dw,nelec, n_wfc_at_ispresent, n_wfc_at, wf_collected, & 
                  fermi_energy_ispresent, fermi_energy/e2, HOL_ispresent, fermi_energy/e2,     &
                  two_fermi_energies, 2, ef_updw/e2, starting_k_points_, ndim_ks_energies,      &
                  occupations_kind_, PRESENT(smearing), smearing_, ndim_ks_energies, ks_objs )
    DO ik=1,ndim_ks_energies
       CALL qes_reset_ks_energies(ks_objs(ik))
    END DO
    CALL qes_reset_k_points_IBZ ( starting_k_points_ ) 
    DEALLOCATE (ks_objs,eigenvalues,occupations)
    END SUBROUTINE qexsd_init_band_structure 
    !
    ! 
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_total_energy(obj,etot,eband,ehart,vtxc,etxc,ewald,degauss,demet,electric_field_corr,&
                                       potentiostat_contr)
    !----------------------------------------------------------------------------------------
    !
    ! 
    IMPLICIT NONE
    ! 
    TYPE (total_energy_type)        :: obj
    REAL(DP),INTENT(IN)             :: etot,eband,ehart,vtxc,etxc,ewald,demet
    REAL(DP),INTENT(IN)             :: degauss
    REAL(DP),OPTIONAL,INTENT(IN)    :: electric_field_corr
    REAL(DP),OPTIONAL,INTENT(IN)    :: potentiostat_contr
    !
    LOGICAL                         :: demet_ispresent
    CHARACTER(LEN=*),PARAMETER      :: TAGNAME="total_energy"
    REAL(DP)                        :: etot_har,eband_har,vtxc_har,etxc_har,ewald_har,&
                                       demet_har,ehart_har,efield_corr

    etot_har  = etot/e2
    eband_har = etot/e2
    vtxc_har  = vtxc/e2
    etxc_har  = etxc/e2
    ewald_har = ewald/e2
    ehart_har = ehart/e2
    IF (PRESENT(electric_field_corr)) THEN
       efield_corr=electric_field_corr/e2
    ELSE
       efield_corr=0.d0
    END IF

    IF (degauss .GT. 0.D0) THEN 
       demet_ispresent=.TRUE.
       demet_har=demet/e2
    ELSE 
       demet_ispresent=.FALSE.
       demet_har=0.d0
    ENDIF
    
    CALL  qes_init_total_energy(obj,TAGNAME,etot_har,eband_ispresent=.TRUE.,eband=eband_har,&
                               ehart_ispresent=.TRUE., ehart=ehart_har,vtxc_ispresent=.TRUE.,& 
                               vtxc=vtxc_har,etxc_ispresent=.TRUE., etxc=etxc_har,ewald_ispresent=.TRUE.,&
                               ewald=ewald_har, demet_ispresent=demet_ispresent,demet=demet_har, &
                               efieldcorr_ispresent=PRESENT(electric_field_corr), efieldcorr=efield_corr,&
                               POTENTIOSTAT_CONTR_ISPRESENT = PRESENT(potentiostat_contr), & 
                               POTENTIOSTAT_CONTR = potentiostat_contr)

    END SUBROUTINE qexsd_init_total_energy
    ! 
    !
    !-------------------------------------------------------------------------------------------------------- 
    SUBROUTINE qexsd_init_forces(obj,nat,forces,tprnfor)
    !-------------------------------------------------------------------------------------------------------- 
    !
    IMPLICIT NONE
    !
    TYPE(matrix_type)                            :: obj
    INTEGER,INTENT(IN)                           :: nat 
    REAL(DP),DIMENSION(:,:),INTENT(IN)           :: forces
    LOGICAL,INTENT(IN)                           :: tprnfor
    !
    CHARACTER(LEN=*),PARAMETER                   :: TAGNAME="forces"
    REAL(DP),DIMENSION(:,:),ALLOCATABLE          :: forces_aux
    ! 
    IF (.NOT. tprnfor) THEN
       obj%lwrite=.FALSE.
       obj%lread =.FALSE.
       RETURN
    END IF 
    !
    ALLOCATE (forces_aux(3,nat))
    forces_aux(1:3,1:nat)=forces(1:3,1:nat)/e2
    !
    CALL qes_init_matrix(obj,TAGNAME,3,nat,forces_aux)
    !
    DEALLOCATE (forces_aux)
    !
    END SUBROUTINE qexsd_init_forces
    ! 
    ! 
    !---------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_stress(obj,stress,tstress) 
    !---------------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE( matrix_type)                           :: obj
    REAL(DP),DIMENSION(3,3),INTENT(IN)           :: stress
    LOGICAL,INTENT(IN)                           :: tstress
    ! 
    CHARACTER(LEN=*),PARAMETER                   :: TAGNAME="stress"
    REAL(DP),DIMENSION(3,3)                      :: stress_aux
    
    IF ( .NOT. tstress ) THEN 
       obj%lwrite = .FALSE.
       obj%lread  = .FALSE.
       stress_aux = 0.d0
       RETURN
    END IF
    ! 
    stress_aux=stress/e2
    CALL qes_init_matrix(obj,TAGNAME,3,3,stress_aux)
    ! 
    END SUBROUTINE qexsd_init_stress
    !
    !
    !------------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_dipole_info (dipole_info, el_dipole, ion_dipole, edir, eamp, emaxpos, eopreg) 
       !------------------------------------------------------------------------------------------------
       ! 
       USE kinds,           ONLY : DP
       USE constants,       ONLY : e2, fpi 
       USE qes_types_module,ONLY : dipoleOutput_type, scalarQuantity_type
       USE qes_libs_module, ONLY : qes_init_scalarQuantity, qes_reset_scalarQuantity
       USE cell_base,       ONLY : alat, at, omega
       ! 
       IMPLICIT NONE  
       ! 
       TYPE ( dipoleOutput_type ), INTENT(OUT)  :: dipole_info
       REAL(DP),INTENT(IN)                      :: el_dipole, ion_dipole, eamp, emaxpos, eopreg
       INTEGER , INTENT(IN)                     :: edir
       ! 
       REAL(DP)                                 :: tot_dipole, length, vamp, fac
       TYPE ( scalarQuantity_type)              :: temp_qobj
       ! 
       tot_dipole = -el_dipole+ion_dipole
       ! 
       dipole_info%idir = edir  
       fac=omega/fpi
       CALL qes_init_scalarQuantity(dipole_info%ion_dipole,"ion_dipole" , units="Atomic Units", &
                                    scalarQuantity= ion_dipole*fac)
       CALL qes_init_scalarQuantity(dipole_info%elec_dipole,"elec_dipole" , units="Atomic Units",&
                                     scalarQuantity= el_dipole*fac)
       CALL qes_init_scalarQuantity(dipole_info%dipole,"dipole" , units="Atomic Units", &
                                    scalarQuantity= tot_dipole*fac)
       CALL qes_init_scalarQuantity(dipole_info%dipoleField,"dipoleField" , units="Atomic Units", &
                                    scalarQuantity= tot_dipole)
       ! 
       length=(1._DP-eopreg)*(alat*SQRT(at(1,edir)**2+at(2,edir)**2+at(3,edir)**2))
       vamp=e2*(eamp-tot_dipole)*length
       !
       CALL qes_init_scalarQuantity(dipole_info%potentialAmp,"potentialAmp" , units="Atomic Units",&
                                     scalarQuantity= vamp)
       CALL qes_init_scalarQuantity(dipole_info%totalLength, "totalLength", units = "Bohr",&
                                     scalarQuantity = length ) 
  
    END SUBROUTINE qexsd_init_dipole_info
    !---------------------------------------------------------------------------------------------
    SUBROUTINE  qexsd_init_outputElectricField(obj, lelfield, tefield, ldipole, lberry, bp_obj, el_pol, &
                                               ion_pol, dipole_obj )
    !---------------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    ! 
    TYPE(outputElectricField_type)                    :: obj 
    ! 
    LOGICAL,INTENT(IN)                                :: lberry, lelfield, tefield, ldipole
    REAL(DP),OPTIONAL,DIMENSION(3),INTENT(IN)         :: el_pol, ion_pol
    TYPE(berryPhaseOutput_type),OPTIONAL,INTENT(IN)   :: bp_obj
    TYPE ( dipoleOutput_type ),OPTIONAL, INTENT(IN)   :: dipole_obj 
    ! 
    CHARACTER(LEN=*),PARAMETER                        :: TAGNAME="electric_field" 
    TYPE ( berryPhaseOutput_type )                    :: bp_loc_obj
    TYPE ( dipoleOutput_type )                        :: dip_loc_obj
    TYPE ( finiteFieldOut_type )                      :: finiteField_obj
    LOGICAL                                           :: bp_is = .FALSE. , finfield_is = .FALSE. , &
                                                         dipo_is = .FALSE.
    ! 
    

    IF (lberry .AND. PRESENT ( bp_obj))  THEN
       bp_is = .TRUE. 
       bp_loc_obj = bp_obj
    END IF 
    IF ( lelfield .AND. PRESENT(el_pol) .AND. PRESENT (ion_pol ) ) THEN 
       finfield_is=.TRUE.
       CALL qes_init_finiteFieldOut (finiteField_obj, "finiteElectricFieldInfo", el_pol, ion_pol)
    END IF 
    IF ( ldipole .AND. PRESENT( dipole_obj ) ) THEN
       dipo_is = .TRUE.
       dip_loc_obj=dipole_obj
    END IF 
    CALL  qes_init_outputElectricField(obj, TAGNAME, BerryPhase_ispresent = bp_is, &
                                       BerryPhase = bp_loc_obj, &
                                       finiteElectricFieldInfo_ispresent = finfield_is, & 
                                       finiteElectricFieldInfo = finiteField_obj, &
                                        dipoleInfo_ispresent = dipo_is, dipoleInfo = dip_loc_obj)
    IF (dipo_is) CALL qes_reset_dipoleOutput( dip_loc_obj )
    IF ( bp_is ) CALL qes_reset_berryPhaseOutput( bp_loc_obj ) 
    !
    END SUBROUTINE qexsd_init_outputElectricField
    ! 
    !----------------------------------------------------------------------------------------
    SUBROUTINE   qexsd_step_addstep( i_step, max_steps, ntyp, atm, ityp, nat,& 
                                   tau, alat, a1, a2, a3, etot, eband, ehart, vtxc, etxc,&
                                   ewald, degauss, demet, forces, stress, n_scf_steps, scf_error, potstat_contr, &
                                   fcp_force, fcp_tot_charge )
    !-----------------------------------------------------------------------------------------
    IMPLICIT NONE 
    ! 
    INTEGER ,INTENT(IN)             :: i_step, max_steps, ntyp, nat, n_scf_steps, ityp(:)
    REAL(DP),INTENT(IN)             :: tau(3,nat), alat, a1(3), a2(3), a3(3), etot, eband, ehart, vtxc, &
                                       etxc, ewald, degauss, demet, scf_error, forces(3,nat), stress(3,3) 
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
    step_obj%n_step = step_counter
    !
    CALL qes_init_scf_conv( scf_conv_obj,"scf_conv", n_scf_steps, scf_error )
    !
    step_obj%scf_conv = scf_conv_obj 
    CALL qes_reset_scf_conv(scf_conv_obj)
    ! 
    CALL qexsd_init_atomic_structure(atomic_struct_obj, ntyp, atm, ityp, nat, tau, &
                                     alat, a1, a2, a3, 0)
    step_obj%atomic_structure=atomic_struct_obj
    CALL qes_reset_atomic_structure( atomic_struct_obj )
    ! 
    CALL qexsd_init_total_energy ( tot_en_obj, etot/e2, eband/e2, ehart/e2, vtxc/e2, etxc/e2, ewald/e2, degauss/e2, &
                                   demet/e2 )
    IF ( PRESENT ( potstat_contr )) THEN  
       tot_en_obj%potentiostat_contr_ispresent = .TRUE. 
       tot_en_obj%potentiostat_contr = potstat_contr/e2 
    END IF  
    step_obj%total_energy=tot_en_obj
    CALL qes_reset_total_energy( tot_en_obj )
    ! 
    CALL  qes_init_matrix( mat_forces, "forces", 3, nat, forces ) 
    step_obj%forces=mat_forces
    CALL qes_reset_matrix ( mat_forces )
    ! 
    CALL qes_init_matrix( mat_stress, "stress", 3, 3, stress ) 
    step_obj%stress = mat_stress
    CALL qes_reset_matrix ( mat_stress ) 
    IF ( PRESENT ( fcp_force ) ) THEN 
       step_obj%FCP_force = fcp_force
       step_obj%FCP_force_ispresent = .TRUE.
       step_obj%FCP_tot_charge = fcp_tot_charge
       step_obj%FCP_tot_charge_ispresent = .TRUE. 
    END IF 
    !  
    ! 
    steps(step_counter) = step_obj
    call qes_reset_step(step_obj)
    END SUBROUTINE qexsd_step_addstep 
    !-------------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_berryPhaseOutput( obj, gpar, gvec, nppstr, nkort, xk, pdl_ion,  &    
                                            mod_ion, pdl_ion_tot, mod_ion_tot, nstring, pdl_elec,  &
                                          mod_elec, wstring, pdl_elec_up, mod_elec_up, pdl_elec_dw,& 
                                          mod_elec_dw, pdl_elec_tot,mod_elec_tot, pdl_tot, mod_tot,&
                                          upol, rmod)
    !---------------------------------------------------------------------------------------------------
    !
    USE ions_base,            ONLY: nat, tau, atm, zv, ityp
    USE cell_base,            ONLY: omega
    USE noncollin_module,     ONLY : noncolin, nspin_lsda
    IMPLICIT NONE 
    ! 
    TYPE (berryPhaseOutput_type)                      :: obj
    REAL(DP),INTENT(IN)                               :: gpar(3), gvec, pdl_ion(nat), pdl_ion_tot, xk(3,*) 

    REAL(DP),INTENT(IN)                               :: pdl_elec(:), pdl_elec_up, pdl_elec_dw, pdl_elec_tot,    & 
                                                         pdl_tot, upol(3), rmod 
    !  
    INTEGER,INTENT(IN)                                :: mod_ion(nat), mod_ion_tot, mod_elec(:), mod_elec_up,    &
                                                         mod_elec_dw, mod_elec_tot, mod_tot, nppstr, nkort, nstring  
    !  
    REAL(DP),INTENT(IN)                               :: wstring(nstring)      
    ! 
#if !defined (__OLDXLM)
    CHARACTER(LEN=*),PARAMETER                        :: TAGNAME = "BerryPhase"
    TYPE ( polarization_type)                         :: tot_pol_obj
    ! 
    TYPE ( electronicPolarization_type),ALLOCATABLE   :: str_pol_obj(:)
    TYPE ( ionicPolarization_type ),    ALLOCATABLE   :: ion_pol_obj(:)
    TYPE ( k_point_type )                             :: kp_obj
    TYPE ( phase_type)                                :: el_phase, ion_phase, tot_phase
    TYPE ( atom_type )                                :: atom_obj
    TYPE ( scalarQuantity_type )                      :: pol_val
    INTEGER                                           :: iat, istring, indstring, ispin 
    CHARACTER(10)                                     :: mod_string
    LOGICAL                                           :: spin_is = .FALSE. 
    !
    ALLOCATE (ion_pol_obj(nat), str_pol_obj(nat)) 
    DO iat =1, nat 
       WRITE(mod_string,'("(mod" ,I1,")")') mod_ion(iat) 
       CALL qes_init_phase(ion_phase,"phase", 0.d0,.FALSE.,0.d0,.FALSE.,TRIM(mod_string),.TRUE., pdl_ion(iat) )
       CALL qes_init_atom(atom_obj,"ion",name=TRIM(atm(ityp(iat))),position_ispresent=.FALSE.,atom = tau(:,iat), &
                          index_ispresent = .FALSE.)
       CALL qes_init_ionicPolarization(ion_pol_obj(iat), "ionicPolarization", atom_obj, zv(ityp(iat)), ion_phase )       
       CALL qes_reset_phase(ion_phase)
       CALL qes_reset_atom(atom_obj)
    END DO
    ! 
    IF ( nspin_lsda .EQ. 2 ) spin_is  = .TRUE.
    DO  istring= 1, nstring
        indstring = 1+(istring-1)*nppstr
        WRITE(mod_string,'("(mod ",I1,")")') mod_elec(istring)
        CALL qes_init_phase(el_phase, "phase", 0.d0, .FALSE., 0.d0, .FALSE., TRIM (mod_string), .TRUE., &
                            pdl_elec(istring) )
        IF (istring .LE. nstring/nspin_lsda) THEN 
           ispin = 1 
        ELSE 
           ispin = 2 
        END IF
        CALL qes_init_k_point(kp_obj, "firstKeyPoint", wstring(istring), .TRUE., "",.FALSE., xk(:,indstring))
        CALL qes_init_electronicPolarization(str_pol_obj(istring),"electronicPolarization", kp_obj, spin_is, ispin, &
                                             el_phase )
        CALL qes_reset_phase ( el_phase ) 
        CALL qes_reset_k_point(kp_obj)
    END DO
    ! 
    WRITE(mod_string,'("(mod ",I1,")")') mod_tot
    CALL qes_init_phase(tot_phase, "totalPhase", pdl_ion_tot, .TRUE. , pdl_elec_tot, .TRUE., TRIM(mod_string), &
                        .TRUE., pdl_tot)
    ! 
    CALL qes_init_scalarQuantity ( pol_val, "polarization", Units="e/bohr^2", scalarQuantity=(rmod/omega)*pdl_tot )
    !
    CALL qes_init_polarization(tot_pol_obj, "polarization", pol_val, modulus = (rmod/omega)*dble(mod_tot), &
                               direction = upol )  
    ! 
    CALL qes_init_berryPhaseOutput( obj, TAGNAME, tot_pol_obj, tot_phase, nat, ion_pol_obj, nstring, str_pol_obj )
    ! 
    DO istring=1,nstring     
       CALL  qes_reset_electronicPolarization(str_pol_obj(istring))
    END DO 
    DEALLOCATE (str_pol_obj)
    DO iat=1, nat
       CALL qes_reset_ionicPolarization(ion_pol_obj(iat))
    END DO
    DEALLOCATE (ion_pol_obj)
    CALL qes_reset_polarization(tot_pol_obj)
    CALL qes_reset_scalarQuantity(pol_val)
    CALL qes_reset_phase(tot_phase) 
#endif 
    !
    END SUBROUTINE qexsd_init_berryPhaseOutput
    !
    !-------------------------------------------------------------------------------------------------         
    SUBROUTINE qexsd_set_status(status_int)
    !-------------------------------------------------------------------------------------------------
    IMPLICIT NONE 
    !
    INTEGER      :: status_int
#if !defined(__OLDXML)  
    CALL qes_init_status( exit_status, "status", status_int)
#endif
    END SUBROUTINE qexsd_set_status 
    !
    !--------------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_set_closed() 
    ! 
    IMPLICIT NONE 
    CHARACTER(LEN=9)                  :: cdate, time_string
    CHARACTER(LEN=12)                 :: date_string
    !
    CALL date_and_tim( cdate, time_string ) 
    date_string = cdate(1:2) // ' ' // cdate(3:5) // ' ' // cdate (6:9)
    CALL qes_init_closed (qexsd_closed_element, "closed", date_string, time_string,&
                          "")
    END SUBROUTINE qexsd_set_closed 
     
!-------------------------------------------------------------------------
!
!-------------------------------------------
! ... read subroutines
!-------------------------------------------
! 





    !
END MODULE qexsd_module

!
!----------------
! ... dummy defs
!----------------
!

! 
! 

