! Copyright (C) 2019 Quantum ESPRESSO foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE qexsd_init
  !----------------------------------------------------------------------------
  !
  ! This module contains some common subroutines used to copy data used by
  ! the Quantum ESPRESSO package into XML format
  !
  ! Written by Paolo Giannozzi, building upon pre-existing code qexsd.f90
  !
  !
  USE kinds,            ONLY : DP
  ! 
  USE qes_types_module
  USE qes_reset_module, ONLY:  qes_reset 
  USE qes_init_module,  ONLY:  qes_init 
  ! FIXME: none of the following modules should be used here
  USE constants,        ONLY : e2
  USE mp_world,         ONLY : nproc
  USE mp_images,        ONLY : nimage,nproc_image
  USE mp_pools,         ONLY : npool
  USE mp_bands,         ONLY : ntask_groups, nproc_bgrp, nbgrp
  USE global_version,   ONLY:  version_number
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  ! type of objects
  ! 
  TYPE (berryPhaseOutput_type), TARGET         :: qexsd_bp_obj
  TYPE (k_points_IBZ_type)                     :: qexsd_start_k_obj 
  TYPE (occupations_type)                      :: qexsd_occ_obj
  !
  PUBLIC :: qexsd_bp_obj, qexsd_start_k_obj, qexsd_occ_obj
  !
  ! public subroutines. They all work in the same way:
  !    call qexsd_init_*( xml object, list of QE variables)
  ! copies QE variables into the xml object
  !
  PUBLIC :: qexsd_init_convergence_info, qexsd_init_algorithmic_info, &
            qexsd_init_atomic_species, qexsd_init_atomic_structure, &
            qexsd_init_symmetries, qexsd_init_basis_set, qexsd_init_dft, &
            qexsd_init_magnetization, qexsd_init_band_structure, & 
            qexsd_init_total_energy, qexsd_init_forces, qexsd_init_stress, &
            qexsd_init_dipole_info, qexsd_init_outputElectricField,   &
            qexsd_init_outputPBC, qexsd_init_gate_info, qexsd_init_hybrid, &
            qexsd_init_dftU, qexsd_init_vdw, qexsd_init_berryPhaseOutput
  !
CONTAINS
  !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_convergence_info(obj, n_scf_steps, scf_has_converged, scf_error, &
                                           optimization_has_converged, n_opt_steps, grad_norm )
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(convergence_info_type)   :: obj
      INTEGER,           INTENT(IN) :: n_scf_steps
      LOGICAL,           INTENT(IN) :: scf_has_converged   
      REAL(DP),          INTENT(IN) :: scf_error
      LOGICAL, OPTIONAL, INTENT(IN) :: optimization_has_converged
      INTEGER, OPTIONAL, INTENT(in) :: n_opt_steps
      REAL(DP),OPTIONAL, INTENT(IN) :: grad_norm
      !
      CHARACTER(27)       :: subname="qexsd_init_convergence_info"
      TYPE(scf_conv_type) :: scf_conv
      TYPE(opt_conv_type),POINTER :: opt_conv  
      !
      NULLIFY(opt_conv) 
      call qes_init (scf_conv, "scf_conv", scf_has_converged, n_scf_steps, scf_error)
      !
      IF ( PRESENT(optimization_has_converged ))  THEN
          !
          IF ( .NOT. PRESENT(n_opt_steps) ) CALL errore(subname,"n_opt_steps not present",10)
          IF ( .NOT. PRESENT(grad_norm) )   CALL errore(subname,"grad_norm not present",10)
          ALLOCATE ( opt_conv) 
          !
          call qes_init (opt_conv, "opt_conv", optimization_has_converged, n_opt_steps, grad_norm)
      ENDIF
      !
      call qes_init (obj, "convergence_info", scf_conv, opt_conv)
      !
      call qes_reset (scf_conv)
      IF (ASSOCIATED(opt_conv)) THEN
         CALL qes_reset (opt_conv)
         NULLIFY ( opt_conv) 
      END IF 
      !
    END SUBROUTINE qexsd_init_convergence_info
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_algorithmic_info(obj, real_space_beta, real_space_q, uspp, paw )
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(algorithmic_info_type)   :: obj
      LOGICAL,           INTENT(IN) :: real_space_beta, real_space_q, uspp, paw
      !
      CALL qes_init (obj, "algorithmic_info", REAL_SPACE_Q = real_space_q, &
                     REAL_SPACE_BETA = real_space_beta, USPP = uspp, PAW = paw)
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
      REAL(DP), OPTIONAL,TARGET, INTENT(IN) :: amass(:)
      REAL(DP), OPTIONAL,TARGET, INTENT(IN) :: starting_magnetization(:)
      REAL(DP), OPTIONAL,TARGET, INTENT(IN) :: angle1(:),angle2(:)
      !
      TYPE(species_type), ALLOCATABLE :: species(:)
      REAL(DP),POINTER  :: amass_    
      REAL(DP),POINTER  :: start_mag_ 
      REAL(DP),POINTER  :: spin_teta  
      REAL(DP),POINTER  :: spin_phi   
      INTEGER   :: i
      
      ALLOCATE(species(nsp))
      NULLIFY ( amass_, start_mag_, spin_teta, spin_phi) 
      !
      DO i = 1, nsp
          !
          IF ( PRESENT(amass) ) THEN 
             IF (amass(i) .GT. 0._DP) amass_=>amass(i)
          END IF 
          IF ( PRESENT(starting_magnetization) ) THEN
             IF (ANY( starting_magnetization(1:nsp) /= 0.0_DP)) start_mag_ => starting_magnetization(i)
          END IF 
          IF ( PRESENT( angle1 ) ) THEN 
             IF (ANY ( angle1(1:nsp) /= 0.0_DP))    spin_teta => angle1(i)
          END IF 
          IF ( PRESENT( angle2 ) )  THEN 
             IF (ANY(angle2(1:nsp) /= 0.0_DP)) spin_phi  => angle2(i)
          END IF 
          !
          CALL qes_init ( species(i), "species", NAME = TRIM(atm(i)), PSEUDO_FILE = TRIM(psfile(i)), MASS = amass_, &
                               STARTING_MAGNETIZATION = start_mag_, SPIN_TETA = spin_teta, SPIN_PHI = spin_phi )
      ENDDO
      !
      CALL qes_init (obj, "atomic_species", nsp, species)
      !
      DO i = 1, nsp
          CALL qes_reset (species(i))
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
      INTEGER,TARGET              :: ibrav_tgt
      INTEGER,POINTER             :: ibrav_ptr 
      CHARACTER(LEN=256),POINTER  :: use_alt_axes_ 
      CHARACTER(LEN=256),TARGET   :: use_alt_axes 
      !
      ! atomic positions
      !
      NULLIFY(use_alt_axes_, ibrav_ptr) 
      IF ( ibrav .ne. 0 ) THEN 
         ibrav_tgt =  abs(ibrav) 
         ibrav_ptr => ibrav_tgt
         use_alt_axes_ => use_alt_axes 
         SELECT CASE(abs(ibrav))
            CASE(-3)
               use_alt_axes="b:a-b+c:-c"
            CASE(-5)
               use_alt_axes="3fold-111"
            CASE(-9)
               use_alt_axes="-b:a:c"
            CASE (91)
               ibrav_tgt = 9 
               use_alt_axes ="bcoA-type"
            CASE(-12,-13)
               use_alt_axes="unique-axis-b" 
            CASE default 
               NULLIFY (use_alt_axes_) 
         END SELECT 
      END IF 
      
      !
      ALLOCATE(atom(nat))
      DO ia = 1, nat
          CALL qes_init ( atom(ia), "atom", name=trim(atm(ityp(ia))), atom=tau(1:3,ia), index = ia )
      ENDDO
      !
      CALL qes_init (atomic_pos, "atomic_positions", atom)
      !
      DO ia = 1, nat
          CALL qes_reset ( atom(ia) )
      ENDDO
      DEALLOCATE(atom)
      !
      ! cell
      !
      CALL qes_init (cell, "cell", a1, a2, a3)
      !
      ! global init
      !
      CALL qes_init (obj, "atomic_structure", NAT=nat, ALAT=alat, &
              ATOMIC_POSITIONS=atomic_pos, CELL=cell , &
              BRAVAIS_INDEX=ibrav_ptr, ALTERNATIVE_AXES = use_alt_axes_ )
      ! 
      ! cleanup 
      ! 
      CALL qes_reset (atomic_pos)
      CALL qes_reset (cell)
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
      CHARACTER(LEN=15),POINTER    :: classname 
      CHARACTER(LEN=256)           :: la_info
      LOGICAL                      :: class_ispresent = .FALSE., time_reversal_ispresent = .FALSE.
      INTEGER                      :: i
      REAL(DP)                     :: mat_(3,3)
      LOGICAL                      :: true_=.TRUE., false_ = .FALSE. 
      LOGICAL,POINTER              :: trev                       
      TARGET                       :: class_names, true_, false_  
      ALLOCATE(symm(nrot))
      NULLIFY( classname, trev) 
      !
      IF ( TRIM(verbosity) .EQ. 'high' .OR. TRIM(verbosity) .EQ. 'medium')  class_ispresent= .TRUE.
      IF ( noncolin  ) time_reversal_ispresent = .TRUE.
      DO i = 1, nrot
          !
          IF  (class_ispresent ) classname => class_names(i)
          IF (time_reversal_ispresent) THEN 
             SELECT CASE (t_rev(i)) 
               CASE (1) 
                  trev => true_ 
               CASE default
                  trev => false_ 
            END SELECT 
         END IF 
          IF ( i .LE. nsym ) THEN 
             la_info = "crystal_symmetry"
          ELSE 
             la_info = "lattice_symmetry"
          END IF
          CALL qes_init (info, "info", name=sname(i), class=classname, time_reversal= trev, INFO= TRIM(la_info) )
          !
          mat_ = real(s(:,:,i),DP)
          CALL qes_init (matrix, "rotation", DIMS=[3,3], mat=mat_ )
          !
          IF ( i .LE. nsym ) THEN 
             CALL qes_init (equiv_atm, "equivalent_atoms", nat=nat, equivalent_atoms = irt(i,1:nat)  )
          !
             CALL qes_init (symm(i),"symmetry", info=info, rotation=matrix, fractional_translation=ft(:,i),  & 
                            equivalent_atoms=equiv_atm)
          ELSE 
             CALL qes_init ( symm(i), "symmetry", INFO = info, ROTATION = matrix ) 
          END IF
          !
          CALL qes_reset (info)
          CALL qes_reset (matrix)
          IF ( i .LT. nsym ) THEN 
             CALL qes_reset ( equiv_atm )
          ELSE IF ( i .EQ. nrot ) THEN  
            CALL qes_reset ( equiv_atm )
          END IF
          !
      ENDDO
      !
      CALL qes_init (obj,"symmetries",NSYM = nsym, NROT=nrot, SPACE_GROUP = space_group, SYMMETRY=symm )
      !
      DO i = 1, nsym
         CALL qes_reset (symm(i))
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

      CALL qes_init (fft_grid, "fft_grid", nr1, nr2, nr3, "")
      CALL qes_init (fft_smooth, "fft_smooth", nr1s, nr2s, nr3s, "")
      CALL qes_init (fft_box, "fft_box", nr1b, nr2b, nr3b, "" )
      CALL qes_init (recipr_latt, "reciprocal_lattice", b1, b2, b3)

      CALL qes_init (obj, "basis_set", GAMMA_ONLY=gamma_only, ECUTWFC=ecutwfc, ECUTRHO=ecutrho, FFT_GRID=fft_grid, &
                     FFT_SMOOTH=fft_smooth, FFT_BOX=fft_box, NGM=ngm, NGMS=ngms, NPWX=npwx,  &
                     RECIPROCAL_LATTICE=recipr_latt )
      !
      CALL qes_reset(fft_grid)
      CALL qes_reset(fft_smooth)
      CALL qes_reset(fft_box)
      CALL qes_reset(recipr_latt)
      !
    END SUBROUTINE qexsd_init_basis_set
    !
    !
    SUBROUTINE  qexsd_init_dft (obj, functional, hybrid_, vdW_, dftU_) 
       IMPLICIT NONE 
       TYPE (dft_type),INTENT(INOUT)           :: obj 
       CHARACTER(LEN=*),INTENT(IN)             :: functional
       TYPE(hybrid_type),OPTIONAL,INTENT(IN)   :: hybrid_ 
       TYPE(vdW_type),OPTIONAL,INTENT(IN)      :: vdW_
       TYPE(dftU_type),OPTIONAL,INTENT(IN)     :: dftU_ 
       ! 
       CALL qes_init(obj, 'dft', functional, hybrid_, dftU_, vdW_)
    END SUBROUTINE qexsd_init_dft 

    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_hybrid ( obj, dft_is_hybrid, nq1, nq2, nq3, ecutfock, exx_fraction, screening_parameter,&
                                   exxdiv_treatment, x_gamma_extrapolation, ecutvcut, local_thr ) 
         IMPLICIT NONE 
         TYPE (hybrid_type),INTENT(INOUT)        :: obj 
         LOGICAL,INTENT(IN)                      :: dft_is_hybrid 
         INTEGER,OPTIONAL, INTENT(IN)            :: nq1, nq2, nq3 
         REAL(DP),OPTIONAL,INTENT(IN)            :: ecutfock, exx_fraction, screening_parameter, ecutvcut,&
                                                    local_thr
         CHARACTER(LEN=*), INTENT(IN)            :: exxdiv_treatment 
         LOGICAL,OPTIONAL,INTENT(IN)             :: x_gamma_extrapolation 
         ! 
         TYPE (qpoint_grid_type),TARGET          :: qpoint_grid 
         TYPE (qpoint_grid_type),POINTER         :: qpoint_grid_opt 
         !
         NULLIFY ( qpoint_grid_opt) 
         IF (.NOT. dft_is_hybrid) RETURN 
         IF (PRESENT(nq1) .AND. PRESENT(nq2) .AND. PRESENT(nq3) ) THEN
            qpoint_grid_opt => qpoint_grid           
            CALL qes_init (qpoint_grid, "qpoint_grid", nq1, nq2, nq3, "")
         END IF 
         !
         CALL qes_init ( obj, "hybrid", qpoint_grid_opt, ecutfock, exx_fraction, &
                        screening_parameter, exxdiv_treatment, x_gamma_extrapolation, ecutvcut,&
                        local_thr )
         !
         IF (ASSOCIATED (qpoint_grid_opt)) CALL qes_reset (qpoint_grid_opt)
         !
      END SUBROUTINE qexsd_init_hybrid 
      !
      SUBROUTINE qexsd_init_dftU (obj, nsp, psd, species, ityp, is_hubbard, lda_plus_u_kind, U_projection_type, &
                                   U, J0, alpha, beta, J, noncolin, starting_ns, Hub_ns, Hub_ns_nc )
         IMPLICIT NONE 
         TYPE(dftU_type),INTENT(INOUT)  :: obj 
         INTEGER,INTENT(IN)             :: nsp
         CHARACTER(LEN=*),INTENT(IN)    :: psd(nsp)
         CHARACTER(LEN=*),INTENT(IN)    :: species(nsp)
         INTEGER,INTENT(IN)             :: ityp(:)
         LOGICAL,INTENT(IN)             :: is_hubbard(nsp)
         INTEGER,INTENT(IN)             :: lda_plus_u_kind
         CHARACTER(LEN=*),INTENT(IN)    :: U_projection_type
         LOGICAL,OPTIONAL,INTENT(IN)    :: noncolin 
         REAL(DP),OPTIONAL,INTENT(IN)   :: U(:), J0(:), alpha(:), beta(:), J(:,:)
         REAL(DP),OPTIONAL,INTENT(IN)   :: starting_ns(:,:,:), Hub_ns(:,:,:,:)
         COMPLEX(DP),OPTIONAL,INTENT(IN) :: Hub_ns_nc(:,:,:,:)
         !
         CHARACTER(10), ALLOCATABLE            :: label(:)
         TYPE(HubbardCommon_type),ALLOCATABLE  :: U_(:), J0_(:), alpha_(:), beta_(:) 
         TYPE(HubbardJ_type),ALLOCATABLE       :: J_(:) 
         TYPE(starting_ns_type),ALLOCATABLE    :: starting_ns_(:) 
         TYPE(Hubbard_ns_type),ALLOCATABLE     :: Hubbard_ns_(:)
         LOGICAL                               :: noncolin_ =.FALSE.
         !
         CALL set_labels ()
         IF ( PRESENT(noncolin)) noncolin_ = noncolin 
         !
         IF (PRESENT(U))   CALL init_hubbard_commons(U, U_, label, "Hubbard_U") 
         IF (PRESENT(J0))  CALL init_hubbard_commons(J0, J0_, label, "Hubbard_J0" ) 
         IF (PRESENT(alpha)) CALL init_hubbard_commons(alpha, alpha_,label, "Hubbard_alpha") 
         IF (PRESENT(beta))  CALL init_hubbard_commons(beta, beta_, label, "Hubbard_beta")
         IF (PRESENT(J))     CALL init_hubbard_J (J, J_, label, "Hubbard_J" )
         IF (PRESENT(starting_ns)) CALL init_starting_ns(starting_ns_ , label)
         IF (PRESENT(Hub_ns))  CALL init_Hubbard_ns(Hubbard_ns_ , label)
         !
         CALL qes_init (obj, "dftU", lda_plus_u_kind, U_, J0_, alpha_, beta_, J_, starting_ns_, Hubbard_ns_, &
                           U_projection_type)
         ! 
         CALL reset_hubbard_commons(U_)
         CALL reset_hubbard_commons(beta_) 
         CALL reset_hubbard_commons(J0_)
         CALL reset_hubbard_commons(alpha_) 
         CALL reset_hubbard_J(J_)
         CALL reset_starting_ns(starting_ns_) 
         CALL reset_Hubbard_ns(Hubbard_ns_) 
      CONTAINS 
         SUBROUTINE set_labels() 
            IMPLICIT NONE 
            CHARACTER                     :: hubbard_shell(4)=['s','p','d','f']
            INTEGER,EXTERNAL              :: set_hubbard_l,set_hubbard_n
            INTEGER                       :: i, hubb_l, hubb_n 
            ! 
            ALLOCATE(label(nsp))
            DO i = 1, nsp
               IF (is_hubbard(i)) THEN
                  hubb_l=set_hubbard_l(psd(i))
                  hubb_n=set_hubbard_n(psd(i))
                  WRITE (label(i),'(I0,A)') hubb_n,hubbard_shell(hubb_l+1) 
               ELSE
                  label(i)="no Hubbard"
               END IF
            END DO
         END SUBROUTINE set_labels 

         SUBROUTINE init_hubbard_commons(dati, objs, labs, tag)
            IMPLICIT NONE
            REAL(DP)   :: dati(:) 
            TYPE(HubbardCommon_type),ALLOCATABLE  :: objs(:)
            CHARACTER(LEN=*) :: labs(:), tag
            INTEGER          :: i
            !

            ALLOCATE (objs(nsp)) 
            DO i = 1, nsp 
               CALL qes_init( objs(i), TRIM(tag), TRIM(species(i)), dati(i), TRIM(labs(i)))
               IF (TRIM(labs(i)) =='no Hubbard') objs(i)%lwrite = .FALSE. 
            END DO 
         END SUBROUTINE init_hubbard_commons 
         !
         SUBROUTINE init_hubbard_J(dati, objs, labs, tag)
            IMPLICIT NONE
            REAL(DP)  :: dati(:,:)
            TYPE(HubbardJ_type),ALLOCATABLE :: objs(:)
            CHARACTER(LEN=*)  :: labs(:), tag 
            INTEGER           :: i
            !
            IF ( SIZE(dati,2) .LE. 0 ) RETURN 
            ALLOCATE (objs(nsp)) 
            DO i = 1, nsp 
               CALL qes_init( objs(i), TRIM(tag), TRIM(species(i)), HubbardJ = dati(1:3,i), LABEL = TRIM(labs(i)))
               IF (TRIM(labs(i)) =='no Hubbard') objs(i)%lwrite = .FALSE. 
            END DO 
         END SUBROUTINE init_hubbard_J
         !
         SUBROUTINE reset_hubbard_commons(objs) 
            IMPLICIT NONE 
            TYPE(HubbardCommon_type),ALLOCATABLE :: objs(:)
            INTEGER  ::   i 
            IF (.NOT. ALLOCATED(objs)) RETURN 
            DO i =1, SIZE(objs) 
               CALL qes_reset(objs(i)) 
            END DO 
            DEALLOCATE(objs)
         END SUBROUTINE reset_hubbard_commons 
         ! 
         SUBROUTINE reset_hubbard_J(objs) 
            IMPLICIT NONE
            TYPE(HubbardJ_type),ALLOCATABLE :: objs(:)
            INTEGER  ::   i
            IF (.NOT. ALLOCATED(objs)) RETURN
            DO i =1, SIZE(objs)
               CALL qes_reset(objs(i))
            END DO
            DEALLOCATE(objs)
         END SUBROUTINE reset_hubbard_J 
         !
         SUBROUTINE init_starting_ns(objs, labs )
            IMPLICIT NONE
            TYPE(starting_ns_type), ALLOCATABLE   :: objs(:)
            CHARACTER(len=*)                      :: labs(nsp)
            INTEGER                               :: i, is, ind, llmax, nspin
            !  
            IF ( .NOT. PRESENT(starting_ns)) RETURN
            
            IF (noncolin_) THEN
               llmax = SIZE(starting_ns,1)
               nspin = 1
               ALLOCATE(objs(nsp))
               DO i = 1, nsp
                  IF (.NOT. ANY(starting_ns(1:2*llmax,1,i)>0.d0)) CYCLE
                  ind = ind + 1 
                  CALL qes_init(objs(ind),"starting_ns", TRIM(species(i)), TRIM(labs(i)), 1, &
                                MAX(starting_ns(1:2*llmax,1,i),0._DP)) 
               END DO 
               RETURN 
            ELSE
               llmax = SIZE (starting_ns, 1)
               nspin = SIZE(starting_ns, 2)
               ALLOCATE(objs(nspin*nsp))
               ind = 0 
               DO is = 1, nspin
                  DO i = 1, nsp 
                     IF (.NOT. ANY (starting_ns(1:llmax,is,i) > 0.d0)) CYCLE 
                     ind = ind + 1 
                     CALL qes_init(objs(ind), "starting_ns", TRIM(species(i)), TRIM (labs(i)), & 
                                    is,MAX(starting_ns(1:llmax,is,i),0._DP))
                  END DO 
               END DO
               RETURN 
            END IF 
         END SUBROUTINE init_starting_ns 
         ! 
         SUBROUTINE init_Hubbard_ns(objs, labs  )
            IMPLICIT NONE 
            TYPE (Hubbard_ns_type),ALLOCATABLE  :: objs(:) 
            CHARACTER(LEN=*)                    :: labs(nsp) 
            !
            REAL(DP), ALLOCATABLE               :: Hubb_occ_aux(:,:) 
            INTEGER                             :: i, is,ind, ldim, m1, m2, llmax, nat, nspin
            ! 
            ! 
            IF (PRESENT(Hub_ns_nc )) THEN
               llmax = SIZE ( Hub_ns_nc, 1) 
               nat = size(Hub_ns_nc,4)
               ALLOCATE (objs(nat))
               ldim = SIZE(Hub_ns_nc,1)
               ALLOCATE (Hubb_occ_aux(2*ldim, 2*ldim)) 
               DO i =1, nat 
                  Hubb_occ_aux = 0._DP 
                  DO m2 = 1, ldim
                     DO m1 =1, ldim 
                        Hubb_occ_aux(m1,m2)=SQRT(DCONJG(Hub_ns_nc(m1,m2,1,i))*Hub_ns_nc(m1,m2,1,i))
                        Hubb_occ_aux(m1,ldim+m2)=SQRT(DCONJG(Hub_ns_nc(m1,m2,2,i))*Hub_ns_nc(m1,m2,2,i))
                        Hubb_occ_aux(ldim+m1,m2)=SQRT(DCONJG(Hub_ns_nc(m1,m2,3,i))*Hub_ns_nc(m1,m2,3,i))
                        Hubb_occ_aux(ldim+m1,ldim+m2)=SQRT(DCONJG(Hub_ns_nc(m1,m2,4,i))*Hub_ns_nc(m1,m2,4,i))
                     END DO
                  END DO
                  CALL qes_init (objs(i), TAGNAME = "Hubbard_ns_mod", SPECIE = TRIM(species(ityp(i))), &
                               LABEL = TRIM(labs(ityp(i))), SPIN =1, INDEX = i,ORDER ='F',Hubbard_NS = Hubb_occ_aux) 
                  IF (TRIM(labs(ityp(i))) == 'no Hubbard') objs(i)%lwrite = .FALSE. 
               END DO
               RETURN 
            ELSE IF (PRESENT (Hub_ns)) THEN
               llmax = SIZE ( Hub_ns,1) 
               nat = size(Hub_ns,4)
               nspin = size(Hub_ns,3)
               ALLOCATE( objs(nspin*nat) )
               ind = 0 
               DO i = 1, nat
                  DO is = 1, nspin
                     ind = ind+1
                     CALL qes_init(objs(ind),"Hubbard_ns", SPECIE = TRIM(species(ityp(i))), SPIN = is, &
                        ORDER = 'F', INDEX = ind, LABEL = TRIM(labs(ityp(i))), Hubbard_NS = Hub_ns(:,:,is,i))
                        IF (TRIM(labs(ityp(i))) =='no Hubbard' )  objs(ind)%lwrite=.FALSE. 
                  END DO
               END DO
            END IF
            RETURN
            !
         END SUBROUTINE init_Hubbard_ns 
         
         SUBROUTINE reset_Hubbard_ns(objs) 
            IMPLICIT NONE 
            ! 
            TYPE(hubbard_ns_type),OPTIONAL    :: objs(:) 
            INTEGER   :: i_ 

            IF ( .NOT. PRESENT (objs)) RETURN 
            DO i_ = 1, SIZE(objs) 
               CALL qes_reset(objs(i_)) 
            END DO 
         END SUBROUTINE reset_Hubbard_ns

         SUBROUTINE reset_starting_ns(obj) 
            IMPLICIT NONE 
            TYPE (starting_ns_type), OPTIONAL  :: obj(:)  
            INTEGER  :: i 
            IF ( .NOT. PRESENT(obj) ) RETURN 
            DO i = 1, SIZE(obj)
               CALL qes_reset(obj(i)) 
            END DO 
         END SUBROUTINE reset_starting_ns 
         !
             
      END SUBROUTINE qexsd_init_dftU 
         ! 
         !
         SUBROUTINE qexsd_init_vdw(obj, non_local_term, vdw_corr, vdw_term, ts_thr, ts_isol,& 
                                   london_s6, london_c6, london_rcut, species, xdm_a1, xdm_a2,&
                                   dftd3_version, dftd3_threebody )
            IMPLICIT NONE 
            TYPE(vdW_type)  :: obj 
            CHARACTER(LEN=*),OPTIONAL,INTENT(IN)            :: non_local_term, vdw_corr
            REAL(DP),OPTIONAL,INTENT(IN)           :: vdw_term, london_c6(:), london_rcut,  xdm_a1, xdm_a2, ts_thr,&
                                                      london_s6
            INTEGER,OPTIONAL,INTENT(IN)                     :: dftd3_version           
            CHARACTER(LEN=*),OPTIONAL              :: species(:)
            LOGICAL,OPTIONAL,INTENT(IN)            :: ts_isol, dftd3_threebody 
            !
            LOGICAL         :: empirical_vdw = .FALSE. , dft_is_vdw  = .FALSE. 
            TYPE(HubbardCommon_type),ALLOCATABLE :: london_c6_obj(:)  
            INTEGER                              :: isp
            ! 
            empirical_vdw = PRESENT(vdw_corr) 
            dft_is_vdw = PRESENT(non_local_term) 
            IF ( .NOT. (dft_is_vdW .OR. empirical_vdw)) RETURN
            IF ( PRESENT (london_c6)) CALL init_londonc6(london_c6, london_c6_obj) 
            CALL qes_init (obj, "vdW", VDW_CORR = vdw_corr, NON_LOCAL_TERM = non_local_term,&
                           TOTAL_ENERGY_TERM = vdw_term, LONDON_S6  = london_s6,& 
                            TS_VDW_ECONV_THR = ts_thr,  TS_VDW_ISOLATED  = ts_isol, LONDON_RCUT = london_rcut, &
                            XDM_A1 = xdm_a1, XDM_A2  = xdm_a2, LONDON_C6 = london_c6_obj, &
                            DFTD3_VERSION = dftd3_version, DFTD3_THREEBODY = dftd3_threebody)
          !
          IF (ALLOCATED(london_c6_obj))   THEN
             DO isp=1, SIZE(london_c6_obj,1) 
                CALL qes_reset(london_c6_obj(isp))
             END DO 
          END IF
          CONTAINS
          ! 
          SUBROUTINE init_londonc6(c6data, c6objs )
            USE constants, ONLY: eps16
            IMPLICIT NONE 
            REAL(DP),INTENT(IN)  :: c6data(:)
            TYPE(HubbardCommon_type),ALLOCATABLE,INTENT(INOUT) :: c6objs(:) 
            ! 
            INTEGER :: ndim_london_c6, isp, ind, nsp
            !
            IF (.NOT. PRESENT ( species)) RETURN
            nsp = SIZE(c6data)
            ndim_london_c6 = COUNT ( c6data .GT. -eps16) 
            IF ( ndim_london_c6 .GT. 0 ) THEN 
               ALLOCATE (c6objs(ndim_london_c6))
               DO isp = 1, nsp
                  IF ( c6data(isp) .GT. -eps16 ) THEN
                     ind  = ind + 1  
                     CALL qes_init(c6objs(ind ), "london_c6", SPECIE = TRIM(species(isp)), HUBBARDCOMMON = c6data(isp))
                  END IF 
               END DO                        
            END IF 
         END SUBROUTINE init_londonc6 
         !   
   END SUBROUTINE qexsd_init_vdw  
    !--------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_outputPBC(obj,assume_isolated)
    !--------------------------------------------------------------------------------------------
    ! 
    IMPLICIT NONE
    ! 
    TYPE (outputPBC_type)                       :: obj
    CHARACTER(LEN=*),INTENT(IN)                  :: assume_isolated
    CHARACTER(LEN=*),PARAMETER                   :: TAGNAME="boundary_conditions"
    !
    CALL qes_init (obj,TAGNAME,ASSUME_ISOLATED =assume_isolated)
    END SUBROUTINE qexsd_init_outputPBC
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
      CALL qes_init(obj, "magnetization", lsda, noncolin, spinorbit, total_mag, absolute_mag, &
                                 do_magnetization)
      !
    END SUBROUTINE qexsd_init_magnetization 
    !
    ! 
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_band_structure(obj, lsda, noncolin, lspinorb, nelec, n_wfc_at,  et, wg, nks, xk, ngk, wk, & 
                                         starting_kpoints, occupations_kind, wf_collected, & 
                                         smearing, nbnd, nbnd_up, nbnd_dw, fermi_energy, ef_updw, homo, lumo)
    !----------------------------------------------------------------------------------------
    IMPLICIT NONE
    !
    TYPE(band_structure_type)               :: obj
    CHARACTER(LEN=*), PARAMETER             :: TAGNAME="band_structure"
    LOGICAL,INTENT(IN)                      :: lsda, noncolin, lspinorb 
    INTEGER,INTENT(IN)                      :: nks, n_wfc_at
    INTEGER,OPTIONAL,INTENT(IN)             :: nbnd, nbnd_up, nbnd_dw 
    REAL(DP),INTENT(IN)                     :: nelec 
    REAL(DP),OPTIONAL,INTENT(IN)            :: fermi_energy, ef_updw(2), homo, lumo
    REAL(DP),DIMENSION(:,:),INTENT(IN)      :: et, wg, xk
    REAL(DP),DIMENSION(:),INTENT(IN)        :: wk
    INTEGER,DIMENSION(:),INTENT(IN)         :: ngk      
    TYPE(k_points_IBZ_type),INTENT(IN)      :: starting_kpoints
    TYPE(occupations_type), INTENT(IN)      :: occupations_kind
    TYPE(smearing_type),OPTIONAL,INTENT(IN) :: smearing
    LOGICAL,INTENT(IN)                      :: wf_collected                    
    ! 
    LOGICAL                                 :: n_wfc_at_ispresent = .TRUE.  
    INTEGER                                 :: ndim_ks_energies, ik
    INTEGER,TARGET                          :: nbnd_, nbnd_up_, nbnd_dw_
    INTEGER,POINTER                         :: nbnd_opt, nbnd_up_opt, nbnd_dw_opt 
    TYPE(k_point_type)                      :: kp_obj
    TYPE(ks_energies_type),ALLOCATABLE      :: ks_objs(:)
    TYPE (k_points_IBZ_type)                :: starting_k_points_ 
    REAL(DP),DIMENSION(:),ALLOCATABLE       :: eigenvalues, occupations
    TYPE (smearing_type)                    :: smearing_ 
    !
    !
    ndim_ks_energies=nks   
    !
    NULLIFY( nbnd_opt, nbnd_up_opt, nbnd_dw_opt)  
    IF ( lsda ) THEN 
       ndim_ks_energies=ndim_ks_energies/2
       nbnd_up_opt => nbnd_up_
       nbnd_dw_opt => nbnd_dw_ 
       IF ( PRESENT(nbnd_up) .AND. PRESENT(nbnd_dw) ) THEN
            nbnd_ = nbnd_up+nbnd_dw
            nbnd_up_ = nbnd_up 
            nbnd_dw_ = nbnd_dw 
       ELSE IF ( PRESENT (nbnd) ) THEN
            nbnd_ = 2*nbnd
            nbnd_up_ = nbnd
            nbnd_dw_  = nbnd
       ELSE
            CALL errore ( "qexsd:qexsd_init_band_structure: ", &
                          "in case of lsda nbnd_up+nbnd_dw or nbnd must be givens as arguments", 10)
       END IF
    ELSE 
       IF (.NOT. PRESENT(nbnd) ) &
          CALL errore ("qexsd:qexsd_init_band_structure:", "lsda is false but needed nbnd argument is missing", 10)
       nbnd_=nbnd
       nbnd_opt => nbnd_ 
    END IF  
    !
    !   
    ALLOCATE(eigenvalues(nbnd_),occupations(nbnd_))
    ALLOCATE(ks_objs(ndim_ks_energies))
    !  
    ks_objs%tagname="ks_energies"
    DO ik=1,ndim_ks_energies
       CALL qes_init(kp_obj,"k_point",WEIGHT = wk(ik), K_POINT = xk(:,ik))
       IF ( lsda ) THEN 
          eigenvalues(1:nbnd_up_)=et(1:nbnd_up_,ik)/e2
          eigenvalues(nbnd_up_+1:nbnd_)=et(1:nbnd_dw_,ndim_ks_energies+ik)/e2
       ELSE 
          eigenvalues(1:nbnd_)= et(1:nbnd_,ik)/e2
       END IF
       !
       !
       IF (lsda) THEN 
          IF ( ABS(wk(ik)).GT.1.d-10) THEN 
             occupations(1:nbnd_up_)=wg(1:nbnd_up_,ik)/wk(ik)
             occupations(nbnd_up_+1:nbnd_)=wg(1:nbnd_dw_,ndim_ks_energies+ik)/wk(ndim_ks_energies+ik)
          ELSE 
             occupations(1:nbnd_up_)=wg(1:nbnd_up_,ik)
             occupations(nbnd_up_+1:nbnd_)=wg(1:nbnd_dw_,ik) 
          END IF            
       ELSE 
          IF (ABS(wk(ik)).GT.1.d-10) THEN
              occupations(1:nbnd_)=wg(1:nbnd_,ik)/wk(ik)
          ELSE
              occupations(1:nbnd_)=wg(1:nbnd_,ik)
          END IF
       END IF
       !
       !
       ks_objs(ik)%k_point = kp_obj
       ks_objs(ik)%npw = ngk(ik)
       CALL  qes_init(ks_objs(ik)%eigenvalues, "eigenvalues",eigenvalues)
       CALL  qes_init(ks_objs(ik)%occupations, "occupations",occupations)
       !
       eigenvalues=0.d0
       occupations=0.d0
       CALL qes_reset(kp_obj)  
    END DO 
    ks_objs%lwrite = .TRUE.
    ks_objs%lread  = .TRUE.
    !
    IF ( PRESENT(smearing) ) smearing_ = smearing
!
    starting_k_points_ = starting_kpoints
    starting_k_points_%tagname = "starting_k_points"
!
! 
   CALL qes_init  (obj, TAGNAME, LSDA = lsda, NONCOLIN = noncolin, SPINORBIT = lspinorb, NBND = nbnd_opt,   &
                   NELEC = nelec, WF_COLLECTED = wf_collected, STARTING_K_POINTS = starting_k_points_,      & 
                   NKS = ndim_ks_energies, OCCUPATIONS_KIND = occupations_kind, KS_ENERGIES = ks_objs,      &
                   NBND_UP = nbnd_up_opt, NBND_DW = nbnd_dw_opt, NUM_OF_ATOMIC_WFC = n_wfc_at,              &
                   FERMI_ENERGY = fermi_energy, HIGHESTOCCUPIEDLEVEL = homo,  TWO_FERMI_ENERGIES = ef_updw, & 
                   SMEARING = smearing, LOWESTUNOCCUPIEDLEVEL = lumo)
    DO ik=1,ndim_ks_energies
       CALL qes_reset(ks_objs(ik))
    END DO
    CALL qes_reset( starting_k_points_ ) 
    DEALLOCATE (ks_objs,eigenvalues, occupations)
    END SUBROUTINE qexsd_init_band_structure 
    !
    ! 
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_total_energy(obj, etot, eband, ehart, vtxc, etxc, ewald, degauss, demet, &
                       electric_field_corr, potentiostat_contr, gate_contribution, dispersion_contribution )
    !----------------------------------------------------------------------------------------
    !
    ! 
    IMPLICIT NONE
    ! 
    TYPE (total_energy_type)        :: obj
    REAL(DP),INTENT(IN)             :: etot, ehart,vtxc,etxc
    REAL(DP),OPTIONAL,INTENT(IN)    :: ewald,demet, eband, degauss 
    REAL(DP),OPTIONAL               :: electric_field_corr
    REAL(DP),OPTIONAL               :: potentiostat_contr
    REAL(DP),OPTIONAL               :: gate_contribution
    REAL(DP),OPTIONAL               :: dispersion_contribution  
    !
    CHARACTER(LEN=*),PARAMETER      :: TAGNAME="total_energy"
    ! 
    CALL  qes_init (obj, TAGNAME, ETOT = etot, EBAND = eband, EHART = ehart, VTXC = vtxc, ETXC = etxc, & 
                    EWALD = ewald, DEMET = demet, EFIELDCORR = electric_field_corr, POTENTIOSTAT_CONTR = potentiostat_contr,  &
                                  GATEFIELD_CONTR = gate_contribution, vdW_term = dispersion_contribution ) 
    
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
    CALL qes_init(obj,TAGNAME,[3,nat],forces_aux )
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
    CALL qes_init(obj,TAGNAME,[3,3],stress_aux )
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
       USE qes_libs_module, ONLY : qes_init
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
       dipole_info%tagname = "dipoleInfo"
       dipole_info%lwrite  = .TRUE.
       dipole_info%lread   = .TRUE.
       CALL qes_init (dipole_info%ion_dipole, "ion_dipole" , units="Atomic Units", &
                                    scalarQuantity= ion_dipole*fac)
       CALL qes_init (dipole_info%elec_dipole,"elec_dipole" , units="Atomic Units",&
                                     scalarQuantity= el_dipole*fac)
       CALL qes_init (dipole_info%dipole,"dipole" , units="Atomic Units", &
                                    scalarQuantity= tot_dipole*fac)
       CALL qes_init (dipole_info%dipoleField,"dipoleField" , units="Atomic Units", &
                                    scalarQuantity= tot_dipole)
       ! 
       length=(1._DP-eopreg)*(alat*SQRT(at(1,edir)**2+at(2,edir)**2+at(3,edir)**2))
       vamp=e2*(eamp-tot_dipole)*length
       !
       CALL qes_init (dipole_info%potentialAmp,"potentialAmp" , units="Atomic Units",&
                                     scalarQuantity= vamp)
       CALL qes_init (dipole_info%totalLength, "totalLength", units = "Bohr",&
                                     scalarQuantity = length )
  
    END SUBROUTINE qexsd_init_dipole_info
    !---------------------------------------------------------------------------------------------
    SUBROUTINE  qexsd_init_outputElectricField(obj, lelfield, tefield, ldipole, lberry, bp_obj, el_pol, &
                                               ion_pol, dipole_obj , gateInfo)
    !---------------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    ! 
    TYPE(outputElectricField_type)                    :: obj 
    ! 
    LOGICAL,INTENT(IN)                                :: lberry, lelfield, tefield, ldipole
    REAL(DP),OPTIONAL,INTENT(IN)                      :: el_pol(:), ion_pol(:)
    TYPE(berryPhaseOutput_type),OPTIONAL,INTENT(IN)   :: bp_obj
    TYPE ( dipoleOutput_type ),OPTIONAL, INTENT(IN)   :: dipole_obj 
    TYPE ( gateInfo_type),OPTIONAL,INTENT(IN)         :: gateInfo
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
       CALL qes_init (finiteField_obj, "finiteElectricFieldInfo", el_pol, ion_pol)
    END IF 
    IF ( ldipole .AND. PRESENT( dipole_obj ) ) THEN
       dipo_is = .TRUE.
       dip_loc_obj=dipole_obj
    END IF 
    CALL  qes_init (obj, TAGNAME,   BerryPhase = bp_obj, &
                                    finiteElectricFieldInfo  = finiteField_obj, &
                                    dipoleInfo = dipole_obj, &
                                    GATEINFO =  gateInfo   )
    IF ( finfield_is) CALL qes_reset ( finiteField_obj) 
    !
    END SUBROUTINE qexsd_init_outputElectricField
    ! 
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
    USE noncollin_module,     ONLY: nspin_lsda
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
    CHARACTER(LEN=*),PARAMETER                        :: TAGNAME = "BerryPhase"
    TYPE ( polarization_type)                         :: tot_pol_obj
    ! 
    TYPE ( electronicPolarization_type),ALLOCATABLE   :: str_pol_obj(:)
    TYPE ( ionicPolarization_type ),    ALLOCATABLE   :: ion_pol_obj(:)
    TYPE ( k_point_type )                             :: kp_obj
    TYPE ( phase_type)                                :: el_phase, ion_phase, tot_phase
    TYPE ( atom_type )                                :: atom_obj
    TYPE ( scalarQuantity_type )                      :: pol_val
    INTEGER                                           :: iat, istring, indstring
    INTEGER,POINTER                                   :: ispin  
    INTEGER, TARGET                                   :: spin_val 
    CHARACTER(10)                                     :: mod_string
    LOGICAL                                           :: spin_is = .FALSE. 
    !
    ALLOCATE (ion_pol_obj(nat))
    ALLOCATE (str_pol_obj(nstring))
    NULLIFY(ispin) 
    DO iat =1, nat 
       WRITE(mod_string,'("(mod" ,I1,")")') mod_ion(iat) 
       CALL qes_init (ion_phase,"phase", modulus = TRIM(mod_string), phase = pdl_ion(iat) )
       CALL qes_init (atom_obj,"ion",name=TRIM(atm(ityp(iat))),atom = tau(:,iat))
       CALL qes_init (ion_pol_obj(iat), "ionicPolarization", atom_obj, zv(ityp(iat)), ion_phase )       
       CALL qes_reset (ion_phase)
       CALL qes_reset (atom_obj)
    END DO
    ! 
    IF ( nspin_lsda .EQ. 2 ) ispin => spin_val 
    DO  istring= 1, nstring
        indstring = 1+(istring-1)*nppstr
        WRITE(mod_string,'("(mod ",I1,")")') mod_elec(istring)
        CALL qes_init(el_phase, "phase", modulus = TRIM (mod_string), phase = pdl_elec(istring) )
        IF ( istring .LE. nstring/nspin_lsda ) THEN 
           spin_val  = 1 
        ELSE  
           spin_val  = 2 
        END IF
        CALL qes_init(kp_obj, "firstKeyPoint", WEIGHT = wstring(istring), K_POINT = xk(:,indstring))
        CALL qes_init(str_pol_obj(istring),"electronicPolarization", kp_obj, el_phase, ispin  )
        CALL qes_reset( el_phase ) 
        CALL qes_reset(kp_obj)
    END DO
    ! 
    WRITE(mod_string,'("(mod ",I1,")")') mod_tot
    CALL qes_init (tot_phase, "totalPhase", IONIC = pdl_ion_tot, ELECTRONIC = pdl_elec_tot,  &
                                                       MODULUS = TRIM(mod_string), PHASE = pdl_tot)
    ! 
    CALL qes_init ( pol_val, "polarization", Units="e/bohr^2", scalarQuantity=(rmod/omega)*pdl_tot )
    !
    CALL qes_init (tot_pol_obj, "totalPolarization", pol_val, modulus = (rmod/omega)*dble(mod_tot), &
                               direction = upol )  
    ! 
    CALL qes_init ( obj, TAGNAME, totalPolarization = tot_pol_obj, totalPhase = tot_phase, & 
                         ionicPolarization = ion_pol_obj, electronicPolarization = str_pol_obj )
    ! 
    DO istring=1,nstring     
       CALL  qes_reset (str_pol_obj(istring))
    END DO 
    DEALLOCATE (str_pol_obj)
    DO iat=1, nat
       CALL qes_reset (ion_pol_obj(iat))
    END DO
    DEALLOCATE (ion_pol_obj)
    CALL qes_reset (tot_pol_obj)
    CALL qes_reset (pol_val)
    CALL qes_reset (tot_phase) 
    !
    END SUBROUTINE qexsd_init_berryPhaseOutput
    !
!-----------------------------------------------------------------------------------
SUBROUTINE qexsd_init_gate_info(obj, tagname, gatefield_en, zgate_, nelec_, alat_, at_, bg_, zv_, ityp_) 
   !--------------------------------------------------------------------------------
   USE kinds,         ONLY : DP
   USE constants,     ONLY : tpi
   !
   IMPLICIT NONE
   TYPE (gateInfo_type),INTENT(INOUT)      :: obj;
   CHARACTER(LEN=*)                        :: tagname
   REAL(DP), INTENT(IN)                    :: gatefield_en, zgate_, alat_, at_(3,3), bg_(3,3), zv_(:), nelec_ 
   INTEGER,INTENT(IN)                      :: ityp_(:) 
   ! 
   REAL(DP)                                :: bmod, area, ionic_charge, gateamp, gate_gate_term
   ! 
   bmod=SQRT(bg_(1,3)**2+bg_(2,3)**2+bg_(3,3)**2)
   ionic_charge = SUM( zv_(ityp_(:)) )
   area = ABS((at_(1,1)*at_(2,2)-at_(2,1)*at_(1,2))*alat_**2) 
   gateamp = (-(nelec_-ionic_charge)/area*tpi)
   gate_gate_term =  (- (nelec_-ionic_charge) * gateamp * (alat_/bmod) / 6.0)
   obj = gateInfo_type( TAGNAME = TRIM(tagname), lwrite = .TRUE., lread = .FALSE., POT_PREFACTOR = gateamp, &
                        GATE_ZPOS = zgate_,  GATE_GATE_TERM = gate_gate_term, GATEFIELDENERGY = gatefield_en) 
   ! 
END SUBROUTINE qexsd_init_gate_info 



 END MODULE qexsd_init
