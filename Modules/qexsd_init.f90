! Copyright (C) 2019 Quantum ESPRESSO foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE qexsd_init
  !----------------------------------------------------------------------------
  !! This module contains some common subroutines used to copy data used by
  !! the Quantum ESPRESSO package into XML format.
  !
  !! Written by Paolo Giannozzi, building upon pre-existing code qexsd.f90.
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
            qexsd_init_dftU, qexsd_init_vdw,     &
            qexsd_init_berryPhaseOutput, qexsd_init_rism3d, qexsd_init_rismlaue,& 
            qexsd_init_esm, qexsd_init_sawtooth_info
 !
CONTAINS
  !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_convergence_info(obj, n_scf_steps, scf_has_converged, scf_error, &
                                           optimization_has_converged, n_opt_steps, grad_norm, wf_collected)
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
      LOGICAL,OPTIONAL,  INTENT(IN) :: wf_collected
      !
      CHARACTER(27)       :: subname="qexsd_init_convergence_info"
      TYPE(scf_conv_type) :: scf_conv
      TYPE(opt_conv_type) :: opt_conv  
      !
      call qes_init (scf_conv, "scf_conv", scf_has_converged, n_scf_steps, scf_error)
      !
      IF ( PRESENT(optimization_has_converged ))  THEN
          !
          IF ( .NOT. PRESENT(n_opt_steps) ) CALL errore(subname,"n_opt_steps not present",10)
          IF ( .NOT. PRESENT(grad_norm) )   CALL errore(subname,"grad_norm not present",10)
          !
          call qes_init (opt_conv, "opt_conv", optimization_has_converged, n_opt_steps, grad_norm)
          call qes_init (obj, "convergence_info", scf_conv, opt_conv, wf_collected)
          call qes_reset (scf_conv)
          CALL qes_reset (opt_conv)
          !
      ELSE 
          !
          call qes_init (obj, "convergence_info", scf_conv, WF_COLLECTED = wf_collected)
          call qes_reset (scf_conv)
          !
      ENDIF
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
      CALL qes_init (obj, "atomic_species", NTYP = nsp, SPECIES = species)
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
                                           alat, a1, a2, a3, ibrav, natomwfc)
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
      INTEGER,OPTIONAL, INTENT(IN) :: natomwfc
      !
      INTEGER         :: ia 
      TYPE(atom_type), ALLOCATABLE :: atom(:)
      TYPE(cell_type) :: cell
      TYPE(atomic_positions_type)  :: atomic_pos
      TYPE(wyckoff_positions_type) :: wyckoff_pos
      REAL(DP)                     :: new_alat
      INTEGER,TARGET              :: ibrav_tgt
      INTEGER,POINTER             :: ibrav_ptr 
      CHARACTER(LEN=256),TARGET   :: use_alt_axes 
      CHARACTER(LEN=256),POINTER  :: use_alt_axes_ 
      !
      ! atomic positions
      !
      NULLIFY(use_alt_axes_, ibrav_ptr) 
      IF ( ibrav .ne. 0 ) THEN 
         ibrav_tgt =  abs(ibrav) 
         ibrav_ptr => ibrav_tgt
         SELECT CASE(ibrav)
            CASE(-3)
               use_alt_axes="b:a-b+c:-c"
               use_alt_axes_ => use_alt_axes 
            CASE(-5)
               use_alt_axes="3fold-111"
               use_alt_axes_ => use_alt_axes 
            CASE(-9)
               use_alt_axes="-b:a:c"
               use_alt_axes_ => use_alt_axes 
            CASE (91)
               ibrav_tgt = 9 
               use_alt_axes ="bcoA-type"
               use_alt_axes_ => use_alt_axes 
            CASE(-12,-13)
               use_alt_axes="unique-axis-b" 
               use_alt_axes_ => use_alt_axes 
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
              BRAVAIS_INDEX=ibrav_ptr, ALTERNATIVE_AXES = use_alt_axes_ , NUM_OF_ATOMIC_WFC= natomwfc)
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
    SUBROUTINE qexsd_init_symmetries(obj, space_group, nsym, nrot, s, ft, &
         sname, t_rev, nat, irt, class_names, verbosity, noncolin, colin_mag_)
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
      INTEGER, OPTIONAL, INTENT(IN) :: colin_mag_
      !
      INTEGER                      :: colin_mag
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

      IF ( PRESENT(colin_mag_) ) THEN
         colin_mag = colin_mag_
      ELSE
         colin_mag = -1
      END IF
      
      IF ( noncolin .OR. (colin_mag == 2) ) time_reversal_ispresent = .TRUE.

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
      CALL qes_init (obj,"symmetries",NSYM = nsym, NROT=nrot, SPACE_GROUP = space_group, & 
        SYMMETRY=symm, COLIN_MAG=colin_mag)
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
            CALL qes_init (qpoint_grid, "qpoint_grid", max(nq1,1), max(nq2,1), max(nq3,1), "")
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
      SUBROUTINE qexsd_init_dftU (obj, nsp, psd, species, ityp, is_hubbard, &
                                  is_hubbard_back, backall, hubb_n2, hubb_l2, hubb_n3, hubb_l3,    & 
                                  noncolin, lda_plus_u_kind, U_projection_type, hubb_occ, U, U2,  J0, J, Um,  &
                                  n, l, alpha, beta, alpha_back, starting_ns, Hub_ns, order_um,Hub_ns_nc, Hub_nsg, Hubbard_V )
         IMPLICIT NONE 
         TYPE(dftU_type),INTENT(INOUT)  :: obj 
         INTEGER,INTENT(IN)             :: nsp
         CHARACTER(LEN=*),INTENT(IN)    :: psd(nsp)
         CHARACTER(LEN=*),INTENT(IN)    :: species(nsp)
         INTEGER,INTENT(IN)             :: ityp(:)
         LOGICAL,INTENT(IN)             :: is_hubbard(nsp)
         LOGICAL,OPTIONAL,INTENT(IN)    :: is_hubbard_back(nsp)
         LOGICAL,OPTIONAL,INTENT(IN)    :: backall(nsp)
         INTEGER,OPTIONAL,INTENT(IN)    :: hubb_n2(nsp)
         INTEGER,OPTIONAL,INTENT(IN)    :: hubb_n3(nsp)
         INTEGER,OPTIONAL,INTENT(IN)    :: hubb_l2(nsp)
         INTEGER,OPTIONAL,INTENT(IN)    :: hubb_l3(nsp) 
         INTEGER,OPTIONAL,INTENT(IN)    :: n(nsp), l(nsp)
         INTEGER,INTENT(IN)             :: lda_plus_u_kind
         CHARACTER(LEN=*),INTENT(IN)    :: U_projection_type
         LOGICAL,OPTIONAL,INTENT(IN)    :: noncolin 
         REAL(DP),OPTIONAL,INTENT(IN)   :: U(:), U2(:), J0(:), alpha(:), alpha_back(:), &
                                           beta(:), J(:,:), hubb_occ(:,:)                           
         REAL(DP),OPTIONAL,INTENT(IN)   :: hubbard_v(:,:,:)
         REAL(DP),OPTIONAL,INTENT(IN)   :: Um(:,:,:)
         INTEGER, OPTIONAL, INTENT(IN)  :: order_um(:,:,:)
         REAL(DP),OPTIONAL,INTENT(IN)   :: starting_ns(:,:,:), Hub_ns(:,:,:,:), Hub_nsg(:,:,:,:)
         COMPLEX(DP),OPTIONAL,INTENT(IN) :: Hub_ns_nc(:,:,:,:)
         !
         CHARACTER(10), ALLOCATABLE            :: label(:)
         TYPE(HubbardCommon_type),ALLOCATABLE  :: U_(:), U2_(:), J0_(:), alpha_(:), &
                                                  alpha_back_(:), beta_(:)
         TYPE(HubbardM_type),ALLOCATABLE       :: Um_(:)
         TYPE(HubbardOcc_type),ALLOCATABLE     :: hubb_occ_(:)
         TYPE(HubbardJ_type),ALLOCATABLE       :: J_(:) 
         TYPE(starting_ns_type),ALLOCATABLE    :: starting_ns_(:) 
         TYPE(Hubbard_ns_type),ALLOCATABLE     :: Hubbard_ns_(:), Hubbard_ns_nc_(:)
         TYPE(HubbardBack_type),ALLOCATABLE    :: Hub_back_(:)
         TYPE(HubbardInterSpecieV_type),ALLOCATABLE :: Hub_V_(:) 
         TYPE(orderUm_type), ALLOCATABLE       :: order_Um_(:)
         LOGICAL                               :: noncolin_ =.FALSE.
         INTEGER                               :: icheck 
         !
         IF (PRESENT(n) .AND. PRESENT(l)) THEN
             CALL set_labels (nsp, n, l)
         ELSE
             ALLOCATE(label(nsp))
             label(:)="no Hubbard"
         ENDIF
         IF (PRESENT(hubb_occ)) CALL init_hubbard_occs(hubb_occ, hubb_occ_, label, n2 = hubb_n2, n3=hubb_n3,  &
                                                        l2 = hubb_l2, l3 = hubb_l3) 
         IF ( PRESENT(noncolin)) noncolin_ = noncolin 
         !
         IF (lda_plus_u_kind == 2 ) THEN 
           IF (PRESENT(hubbard_v))   icheck =  check_and_init_Hubbard_V (hub_v_, hubbard_v, species, label)
         ELSE  
           IF (PRESENT(U))           CALL init_hubbard_commons(U, U_, label, "Hubbard_U")
           IF (PRESENT(is_hubbard_back)) THEN
            IF (ANY(is_hubbard_back)) THEN
              IF (.NOT.PRESENT(hubb_l2) .OR. .NOT.PRESENT(hubb_n2) .OR. .NOT. PRESENT(U2)) &
                 CALL errore('qexsd_init_dft:',&
                        'Internal error: second Hubbard channel is present but hubb_n2 or hubb_l2 or U2 is not present',1)
              CALL init_Hubbard_back(is_hubbard_back, Hub_back_, U2, hubb_n2, hubb_l2, backall, n3 = hubb_n3, l3 = hubb_l3) 
            ENDIF
         END IF
         END IF 
         IF (PRESENT(J0))          CALL init_hubbard_commons(J0, J0_, label, "Hubbard_J0" ) 
         IF (PRESENT(alpha))       CALL init_hubbard_commons(alpha, alpha_,label, "Hubbard_alpha") 
         IF (PRESENT(alpha_back))  CALL init_hubbard_commons(alpha_back, alpha_back_,label, "Hubbard_alpha_back") 
         IF (PRESENT(beta))        CALL init_hubbard_commons(beta, beta_, label, "Hubbard_beta")
         IF (PRESENT(J))           CALL init_hubbard_J (J, J_, label, "Hubbard_J" )
         IF (PRESENT(Um))  THEN 
           CALL init_hubbardM (Um, Um_, label, "Hubbard_Um") 
           IF (present(order_um)) CALL init_orderUm  (order_um, order_Um_, label, "Hub_m_order") 
         END IF 
         IF (PRESENT(starting_ns)) CALL init_starting_ns(starting_ns_ , label)
         IF (PRESENT(Hub_ns))   THEN 
                                   CALL init_Hubbard_ns(Hubbard_ns_ , label, Hub_ns)
         ELSE IF (PRESENT(Hub_ns_nc))  THEN 
                                   CALL init_Hubbard_ns(Hubbard_ns_nc_ , label) 
         ELSE IF (PRESENT(Hub_nsg)) THEN 
                                   CALL init_Hubbard_ns(Hubbard_ns_, label, Hub_nsg) 
         END IF 
         
         !
         CALL qes_init (obj, "dftU", .true., lda_plus_u_kind, hubb_occ_,  U_, Um_, J0_, alpha_, beta_,  J_, & 
                       starting_ns_, Hub_V_, Hubbard_ns_, order_Um_, U_projection_type, Hub_back_, alpha_back_, Hubbard_ns_nc_)
         ! 
         CALL reset_hubbard_occs(hubb_occ_)
         CALL reset_hubbard_commons(U_)
         CALL reset_hubbardM(Um_) 
         CALL reset_order_Um(order_Um_) 
         CALL reset_hubbard_commons(U2_)
         CALL reset_hubbard_commons(beta_) 
         CALL reset_hubbard_commons(J0_)
         CALL reset_hubbard_commons(alpha_) 
         CALL reset_hubbard_commons(alpha_back_) 
         CALL reset_hubbard_J(J_)
         CALL reset_starting_ns(starting_ns_) 
         CALL reset_Hubbard_ns(Hubbard_ns_) 
         !
      CONTAINS 
         SUBROUTINE set_labels(ldim, n_, l_) 
            IMPLICIT NONE 
            CHARACTER                     :: hubbard_shell(4)=['s','p','d','f']
            INTEGER                       :: i, hubb_l, hubb_n 
            INTEGER                       :: ldim
            INTEGER                       :: n_(ldim), l_(ldim)
            ! 
            ALLOCATE(label(nsp))
            DO i = 1, nsp
               IF (is_hubbard(i)) THEN
                  hubb_n=n_(i)
                  hubb_l=l_(i) 
                  WRITE (label(i),'(I0,A)') hubb_n,hubbard_shell(hubb_l+1) 
               ELSE
                  label(i)="no Hubbard"
               ENDIF
            ENDDO
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
               CALL qes_init( objs(i), TRIM(tag), TRIM(species(i)), TRIM(labs(i)), dati(i))
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
         SUBROUTINE init_hubbardM(dati, objs, labs, tag) 
           IMPLICIT NONE 
           real(dp) :: dati(:,:,:)
           type(hubbardM_type),allocatable :: objs(:) 
           character(len=*) :: labs(:), tag
           !
           integer :: nhubm, ihubm,nt, nspin, msize, it, ispin, nhubmtot, iobj, ldim
           integer, allocatable :: packdati(:),channels_per_specimen(:),hubm(:) 
           real(dp)  :: uvalue
           !
           nspin = size(dati,2) 
           nhubm = count([(any(dati(:,:,it)/=0.0_dp),it=1,nsp)])
           msize = size(dati,1) 
           if (nhubm .EQ.  0) RETURN
           allocate(channels_per_specimen(nhubm), &
                   packdati(nhubm),               &
                   hubm(size(dati,1))             &
                   )  
           packdati = pack([(it,it=1,nsp)], mask=[(any(dati(:,:,it)/=0.0),it=1,nsp)])
           do ihubm =1, nhubm
               channels_per_specimen(ihubm)  = count([(any(dati(:,ispin, packdati(ihubm))/=0),ispin=1, nspin)])
           end do 
           nhubmtot = sum(channels_per_specimen(:))  
           allocate(objs(nhubmtot)) 
           !
           do ihubm = 1, nhubm 
             it = packdati(ihubm)
             ldim = 2 * l(it) + 1 
             if (noncolin ) ldim = 2 * ldim
             iobj = 1  
             if (ihubm .gt. 1)  iobj = iobj + sum(channels_per_specimen(1:ihubm-1)) 
             do ispin = 1, nspin
               if (ANY(dati(:,ispin,ihubm)/=0.0_DP)) THEN 
                  call qes_init(objs(iobj), TRIM(tag), TRIM(species(it)),TRIM(labs(it)), ispin, &
                                HubbardM = dati(1:ldim,ispin, ihubm)) 
                  if (nspin == 1) objs(iobj)%spin_ispresent = .FALSE. 
                  iobj = iobj + 1
               end if
             end do 
           end do
         END SUBROUTINE init_hubbardM
         !
         SUBROUTINE init_orderUm(dati, objs, labs, tag) 
           IMPLICIT NONE 
           integer  :: dati(:,:,:) 
           type(orderUm_type),allocatable :: objs(:) 
           character(len=*)  :: labs(:), tag 
           ! 
           integer :: iat, nat, nspin, ldim, ispin, iobj
           character(len=10)  :: label_
           character(len=4)  :: lab2ldim = 'spdf' 
           !  
           nat = size(dati,3)
           nspin = size(dati,2) 
           allocate (objs(nspin * nat)) 

           do iat = 1, nat
              label_ = labs(ityp(iat)) 
              ldim  = 2 * index(lab2ldim, label_(2:2)) - 1
              if (noncolin_) ldim = ldim * 2 
              do ispin = 1, nspin 
                 iobj = nspin*(iat-1) + ispin  
                 if (ldim >= 1) &
                   call qes_init(objs(iobj), tagname=tag, specie=TRIM(species(ityp(iat))), label=TRIM(label_),&
                                 spin=ispin, atomidx=iat, OrderUm = dati(1:ldim,ispin,iat)) 
                 if (label_ == "no hubbard") objs(iobj)%lwrite = .FALSE. 
                 if (nspin == 1) objs(iobj)%spin_ispresent = .FALSE. 
              end do 
           end do    
         END SUBROUTINE init_orderUm  
         !
         FUNCTION check_and_init_Hubbard_V (objs, hubbard_v_, specs, labs) result( ndim ) 
           IMPLICIT NONE 
           TYPE(HubbardInterSpecieV_type), ALLOCATABLE   :: objs(:)
           REAL(DP)                                      :: hubbard_v_(:,:,:) 
           CHARACTER(len=*)                              :: labs(:), specs(:)
           INTEGER                                       :: ndim 
           !
           INTEGER    :: nat_, nbt_, na, nb, idim, nb2isp
           CHARACTER(LEN=4)  :: lab1, spec1, lab2, spec2, lab1_bkg = "no_H", lab2_bkg='no_H' 
           CHARACTER         :: hubbard_shell(4) = ['s','p','d','f'] 
           !
           nat_ = SIZE(ityp) 
           nbt_  = SIZE(hubbard_v_, 2) / SIZE( hubbard_v_, 1) * nat_ 
           !
           ndim = COUNT( hubbard_v_(:,:,1:4) /= 0._DP ) 
           IF (ndim == 0 ) RETURN 
           ALLOCATE (objs(ndim)) 

           idim = 0 
           DO na =1, nat_ 
             spec1 = TRIM(species(ityp(na)))
             lab1  = TRIM(label(ityp(na))) 
             IF (PRESENT(hubb_n2) .AND.PRESENT(hubb_l2)) THEN  
               IF (hubb_n2(ityp(na)) /= -1 .AND. hubb_l2(ityp(na)) /= -1) & 
                  WRITE(lab1_bkg,'(I0,A)') hubb_n2(ityp(na)), hubbard_shell( hubb_l2(ityp(na)) + 1 ) 
             END IF      
             DO nb = 1, nbt_ 
              IF (ALL(hubbard_v(na,nb,:) == 0._DP)) CYCLE 
              nb2isp = ityp ( mod(nb -1, nat_) +1)  
              spec2 = TRIM(species(nb2isp))
              lab2  = TRIM(label (nb2isp)) 
              IF (PRESENT(hubb_n2) .AND. PRESENT(hubb_l2)) THEN 
               IF ( hubb_n2(nb2isp) /= -1 .AND. hubb_l2(nb2isp) /= -1 ) & 
                 WRITE (lab2_bkg, '(I0,A)') hubb_n2(nb2isp), hubbard_shell( hubb_l2(nb2isp)  + 1 )
              END IF 
              IF (hubbard_v(na, nb, 1) /= 0._DP) THEN  
                idim = idim + 1 
                CALL qes_init(objs(idim), "Hubbard_V", spec1, na, lab1, spec2, nb, lab2, Hubbard_V_(na,nb,1))
              END IF 
              IF ( hubbard_v(na, nb, 2) /= 0._DP) THEN 
                idim = idim + 1 
                CALL qes_init (objs(idim), "Hubbard_V", spec1, na, lab1, spec2, nb, lab2_bkg, Hubbard_V_(na,nb,2)) 
              END IF 
              IF (hubbard_v(na, nb, 3) /= 0._DP  ) THEN 
               idim = idim + 1 
               CALL qes_init(objs(idim), "Hubbard_V", spec1, na, lab1_bkg, spec2, nb, lab2_bkg, Hubbard_V_(na,nb,3)) 
              END IF 
              IF (hubbard_v(na,nb, 4) /= 0._DP ) THEN 
               idim = idim + 1 
               CALL qes_init(objs(idim), "Hubbard_V", spec1, na, lab1_bkg, spec2, nb, lab2, Hubbard_V_(na,nb,4)) 
              END IF   
             END DO
           END DO 
         END FUNCTION check_and_init_Hubbard_V 
         !
         !
         SUBROUTINE reset_hubbard_occs(objs) 
            IMPLICIT NONE 
            TYPE (HubbardOcc_type),INTENT(INOUT),ALLOCATABLE  :: objs(:)
            !
            INTEGER :: i
            IF (.NOT. ALLOCATED(objs)) RETURN 
            DO i =1, SIZE (objs)
               CALL qes_reset(objs(i))
            END DO 
            DEALLOCATE (objs)
         END SUBROUTINE reset_hubbard_occs 
         !
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
         SUBROUTINE reset_hubbardM(objs) 
           IMPLICIT NONE 
           TYPE(HubbardM_type), ALLOCATABLE :: objs(:) 
           INTEGER :: i 
           IF (.NOT. ALLOCATED(objs)) RETURN 
           DO i =1, SIZE(objs) 
             CALL qes_reset(objs(i)) 
           END DO 
           DEALLOCATE(objs)
         END SUBROUTINE reset_hubbardM
         !
         SUBROUTINE reset_order_Um(objs) 
           IMPLICIT NONE 
           TYPE(orderUm_type), ALLOCATABLE :: objs(:)
           INTEGER :: i
           IF (.NOT. ALLOCATED(objs)) RETURN 
           DO i =1, SIZE(objs)
              CALL qes_reset(objs(i)) 
           END DO 
           DEALLOCATE(objs) 
         END SUBROUTINE reset_order_Um
         !
         !
         SUBROUTINE init_starting_ns(objs, labs )
            IMPLICIT NONE
            TYPE(starting_ns_type), ALLOCATABLE   :: objs(:)
            REAL(DP), ALLOCATABLE                 :: dati(:)  
            CHARACTER(len=*)                      :: labs(nsp)
            INTEGER                               :: i, is, ind, llmax, nspin
            !  
            IF ( .NOT. PRESENT(starting_ns)) RETURN
            
            IF (noncolin_) THEN
               llmax = SIZE(starting_ns,1) 
               nspin = 1
               ALLOCATE(objs(nsp))
               ind = 0
               DO i = 1, nsp
                  IF (.NOT. ANY(starting_ns(1:llmax,1:2,i)>0.d0)) CYCLE
                  ind = ind + 1 
                  ALLOCATE (dati(2*llmax)) 
                  dati(1:llmax) =         MAX(starting_ns(1:llmax,1,i),0._DP)  
                  dati(llmax+1:2*llmax) = MAX(starting_ns(1:llmax,2,i),0._DP)  
                  CALL qes_init(objs(ind),"starting_ns", TRIM(species(i)), TRIM(labs(i)), 1, dati) 
                  DEALLOCATE(dati) 
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
         SUBROUTINE init_Hubbard_ns(objs, labs, hub_ns_)
            IMPLICIT NONE 
            TYPE (Hubbard_ns_type),ALLOCATABLE  :: objs(:) 
            CHARACTER(LEN=*)                    :: labs(nsp) 
            REAL(DP),OPTIONAL,INTENT(IN)        :: hub_ns_(:,:,:,:) 
            !
            REAL(DP), ALLOCATABLE               :: Hubb_occ_aux(:,:) 
            INTEGER                             :: i, is,ind, ldim, m1, m2, llmax, nat, nspin
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
                  CALL qes_init (objs(i), TAGNAME = "Hubbard_ns_nc", SPECIE = TRIM(species(ityp(i))), &
                               LABEL = TRIM(labs(ityp(i))), SPIN =1, INDEX = i,ORDER ='F',Hubbard_NS = Hubb_occ_aux) 
                  IF (TRIM(labs(ityp(i))) == 'no Hubbard') objs(i)%lwrite = .FALSE. 
               END DO
               RETURN 
            ELSE IF (PRESENT (hub_ns_)) THEN
               llmax = SIZE ( hub_ns_,1) 
               nat = size(hub_ns_,4)
               nspin = size(hub_ns_,3)
               ALLOCATE( objs(nspin*nat) )
               ind = 0 
               DO i = 1, nat
                  DO is = 1, nspin
                     ind = ind+1
                     CALL qes_init(objs(ind),"Hubbard_ns", SPECIE = TRIM(species(ityp(i))), SPIN = is, &
                        ORDER = 'F', INDEX = i, LABEL = TRIM(labs(ityp(i))), Hubbard_NS = hub_ns_(:,:,is,i))
                        IF (TRIM(labs(ityp(i))) =='no Hubbard' )  objs(ind)%lwrite=.FALSE. 
                        IF (nspin ==1 ) objs(ind)%spin_ispresent = .FALSE.
                  END DO
               END DO
            END IF
            RETURN
            !
         END SUBROUTINE init_Hubbard_ns 
         
         SUBROUTINE init_Hubbard_back(is_back, objs, U2, n2, l2, backall_, n3, l3) 
            IMPLICIT NONE
            LOGICAL, INTENT(IN)                                  :: is_back(nsp) 
            REAL(DP)                                             :: U2 (nsp)
            INTEGER, INTENT(IN)                                  :: n2(nsp)
            INTEGER, INTENT(IN)                                  :: l2(nsp)
            TYPE(HubbardBack_type),ALLOCATABLE,INTENT(INOUT)     :: objs(:) 
            LOGICAL,OPTIONAL,INTENT(IN)                          :: backall_(nsp) 
            INTEGER,OPTIONAL,INTENT(IN)                          :: l3(nsp)  
            INTEGER,OPTIONAL,INTENT(IN)                          :: n3(nsp)
            !
            INTEGER  :: isp, il, ndimbackL, l2_, l3_, n2_, n3_
            REAL(DP) :: u2_   
            LOGICAL,ALLOCATABLE  :: temp(:) 
            CHARACTER(LEN=16)    :: backchar
            !
            l2_ = -1; l3_ = -1; n2_ = -1; n3_ = -1; u2_ = 0._DP  
            ALLOCATE(objs(nsp), temp(nsp)) 
            IF (PRESENT(backall_)) THEN 
               temp(1:nsp)  = backall_(1:nsp)
            ELSE 
               temp(1:nsp) = .FALSE.
            END IF 
            DO isp =1, nsp 
               u2_ = U2(isp)
               l2_ = l2(isp)
               n2_ = n2(isp) 
               ndimbackL = 1 
               IF (temp(isp) .AND. PRESENT(l3) ) THEN 
                  IF (l3(isp) >=0) THEN
                     ndimbackL=2 
                     l3_ = l3(isp) 
                     IF (PRESENT(n3)) n3_ = n3(isp)
                  ELSE
                     CALL errore ('qexsd_init_dftU:', 'internal error: l3 < 0',1)
                  END IF
               ELSEIF (temp(isp) .AND. .NOT.PRESENT(l3) ) THEN
                     CALL errore ('qexsd_init_dftU:', 'internal error: backall is true but l3 is not present',1) 
               END IF 
               IF (temp(isp)) THEN
                  backchar = 'two_orbitals'
               ELSE 
                  backchar = 'one_orbital'
               END IF 

               CALL qes_init(objs(isp), "Hubbard_back", SPECIES = TRIM(species(isp)), Hubbard_U2= u2_ , &
                             background=TRIM(backchar),  l2_number = l2_, l3_number=l3_, n2_number = n2_, n3_number = n3_)    
               IF (.NOT. is_back(isp)) THEN 
                 objs(isp)%lwrite = .FALSE.
               ELSE IF (.NOT. temp(isp)) THEN 
                 objs(isp)%l3_number_ispresent = .FALSE. 
                 objs(isp)%n3_number_ispresent = .FALSE. 
               END IF 
            END DO 
         END SUBROUTINE init_Hubbard_back 

         

         SUBROUTINE init_hubbard_occs(data, objs, labels_, n2, n3, l2, l3)
            IMPLICIT NONE 
            ! 
            REAL(DP),INTENT(IN)  :: data(:,:)
            TYPE(HubbardOcc_type),ALLOCATABLE, INTENT(INOUT)  :: objs(:)
            CHARACTER(LEN=10),INTENT(IN) :: labels_(:)
            INTEGER,OPTIONAL,INTENT(IN)  :: n2(:), n3(:), l2(:), l3(:) 
            !
            CHARACTER(LEN=10) :: lbl_(3) 
            ! 
            INTEGER  :: i, ich, idx_i, ndim, nchannels
            INTEGER, ALLOCATABLE :: idx(:)
            LOGICAL :: n2_ispresent, n3_ispresent
            TYPE(ChannelOcc_type)   :: channels(3)
            CHARACTER   :: hubbard_shell(4) = ['s','p','d','f']
            n2_ispresent  = PRESENT(n2) .AND. PRESENT(l2)
            n3_ispresent  = PRESENT(n3) .AND. PRESENT(l3) 
            allocate(idx(COUNT(labels_(1:nsp) /= "no Hubbard")))
            idx = PACK( [(i,i=1, nsp)] , [(INDEX(labels_(i),"no")==0,i=1,nsp)])
            ndim = SIZE(idx)
            ALLOCATE (objs(ndim))
            DO i = 1, ndim
               idx_i = idx(i)
               nchannels = 1 
               lbl_(1) = TRIM(labels_(idx_i))
               IF (n2_ispresent ) THEN  
                  IF ( n2(idx_i) > 0  ) THEN 
                    nchannels = nchannels + 1
                    WRITE (lbl_(nchannels),'(I0,A)') n2(idx_i),hubbard_shell(l2(idx_i)+1) 
                  END IF
               END IF 
               IF (n3_ispresent ) THEN 
                  IF ( n3(idx_i) > 0  ) THEN 
                    nchannels = nchannels + 1
                    WRITE (lbl_(nchannels),'(I0,A)') n3(idx_i),hubbard_shell(l3(idx_i)+1) 
                  END IF
               END IF  
               DO ich =1, nchannels 
                  CALL qes_init(channels(ich), 'channel_occ', species(idx_i), lbl_(ich), ich, data(idx_i, ich))
               END DO  
               CALL qes_init(objs(i),"Hubbard_Occ", nchannels, species(idx_i), channels(1:nchannels))
               DO ich = 1, nchannels
                  CALL qes_reset(channels(ich))
               END DO 
            END DO 
         END SUBROUTINE init_hubbard_occs 
         !
         !
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
               ind = 0
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
    SUBROUTINE qexsd_init_esm(esm_obj, bc, nfit, w, efield, a, zb, debug, debug_gpmax) 
      IMPLICIT NONE 
      TYPE(esm_type),INTENT(INOUT)  :: esm_obj 
      CHARACTER(LEN=*),INTENT(IN)  :: bc
      INTEGER,INTENT(IN)            :: nfit
      REAL(DP),INTENT(IN)           :: w, efield, a 
      REAL(DP),INTENT(IN),OPTIONAL  :: zb 
      LOGICAL,INTENT(IN),OPTIONAL   :: debug, debug_gpmax
      !
      CALL qes_init (esm_obj, "esm", bc=TRIM(bc), nfit=nfit, w=w, efield=efield, a=a )
    END SUBROUTINE qexsd_init_esm
    !
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_magnetization(obj, lsda, noncolin, spinorbit, total_mag, total_mag_nc, absolute_mag, &
                                       atm, ityp, site_mag_pol, site_mag, site_charges, do_magnetization)
      !------------------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(magnetization_type)     :: obj
      !! magnetization object to initialize
      LOGICAL,          INTENT(IN) :: lsda, noncolin, spinorbit
      !! flag: true for spin-polarized calculations 
      !! flag: true for noncolinear calculations 
      !! flag: true for fully relativistic calculations 
      REAL(DP),OPTIONAL,INTENT(IN) :: total_mag
      !! total scalar magnetization present for  spin-polarized calculation 
      REAL(DP),         INTENT(IN) :: absolute_mag
      !! sum  of the magnetization's module  (scalar or colinear)
      CHARACTER(LEN=*),INTENT(IN)  :: atm(:)
      !! species labels 
      INTEGER,INTENT(IN)           :: ityp(:)
      !! species index for each atom 
      REAL(DP),OPTIONAL,INTENT(IN) :: total_mag_nc(3)
      !! total magnetization (vector) present for noncolinear magnetic calculations
      REAL(DP),OPTIONAL,INTENT(IN) :: site_charges(:) 
      !! array with the estimate charge of the site
      LOGICAL,OPTIONAL, INTENT(IN) :: do_magnetization 
      !! flag: true if the noncolinear calculation has a finite magnetization 
      REAL(DP),DIMENSION(:,:),OPTIONAL,INTENT(IN) :: site_mag 
      !! array with the estimated magnetization per site (noncollinear case)
      REAL(DP),DIMENSION(:,:)  ,OPTIONAL,INTENT(IN) :: site_mag_pol
      !! array with the magnetic polarization per site (nspin=2 case)
      !
      INTEGER  :: iobj
      TYPE (scalmags_type),TARGET  :: smag_obj 
      TYPE (scalmags_type),POINTER :: smag_ptr 
      TYPE (d3mags_type),  TARGET  :: vmag_obj 
      TYPE (d3mags_type),  POINTER :: vmag_ptr 
      NULLIFY(smag_ptr, vmag_ptr) 
      IF (PRESENT(site_mag_pol)) THEN 
         CALL qexsd_init_scalmags(smag_obj, SIZE(site_mag_pol,2), site_mag_pol(1,:), ityp, atm, site_charges) 
         smag_ptr => smag_obj 
      ELSE IF (PRESENT(site_mag)) THEN 
         CALL qexsd_init_d3mags (vmag_obj, SIZE(site_mag, 2), site_mag, ityp, atm, site_charges ) 
         vmag_ptr => vmag_obj
      END IF       
      CALL qes_init(obj, "magnetization", lsda, noncolin, spinorbit, absolute_mag, total_mag, total_mag_nc,&
                                 smag_ptr, vmag_ptr, do_magnetization)
      !
    END SUBROUTINE qexsd_init_magnetization 
    !
    !----------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_scalmags(obj, nat_, data,  ityp, atm, charges)
      !! stores site scalar magnetization into a structure for XML output
    !-------------------------------------------------------------------------------------------
    IMPLICIT NONE 
    TYPE (scalmags_type),INTENT(INOUT) :: obj 
    !! structure where to store data
    INTEGER                            :: nat_ 
    !! number of atomic sites 
    REAL(DP), INTENT(IN)               :: data (:)
    !! site magnetizations 
    REAL(DP), OPTIONAL,INTENT(IN)      :: charges(:)
    !! site charges 
    INTEGER,INTENT(IN)                 :: ityp(:)
    !! species index for each atom 
    CHARACTER(LEN=*),INTENT(IN)        :: atm(:)
    !! species labels 
    ! 
    INTEGER   :: ia 
    TYPE(SiteMoment_type),ALLOCATABLE  :: site_obj(:)
    ALLOCATE (site_obj(nat_))
    DO ia = 1, nat_
      CALL qes_init(site_obj(ia), "SiteMagnetization", SPECIES=atm(ityp(ia)), ATOM=ia, CHARGE = charges(ia), &
                    SiteMoment= data(ia))
    END DO 
    CALL qes_init(obj, "Scalar_Site_Magnetic_Moments", NAT = nat_, SiteMagnetization=site_obj )
    DEALLOCATE (site_obj)
    END SUBROUTINE qexsd_init_scalmags
    !
    !
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_d3mags (obj, nat_, data, ityp, atm, charges)
    !---------------------------------------------------------------------------------------
    IMPLICIT NONE 
    TYPE (d3mags_type),INTENT(INOUT) :: obj 
    !! structure where to store data
    INTEGER                            :: nat_ 
    !! number of atomic sites 
    REAL(DP), INTENT(IN)               :: data (:,:)
    !! site magnetizations (vectors)
    REAL(DP), OPTIONAL,INTENT(IN)      :: charges(:)
    !! site charges 
    INTEGER,INTENT(IN)                 :: ityp(:)
    !! species index for each atom 
    CHARACTER(LEN=*),INTENT(IN)        :: atm(:)
    !! species labels 
    ! 
    INTEGER   :: ia 
    TYPE(SitMag_type),ALLOCATABLE  :: site_obj(:)
    ALLOCATE (site_obj(nat_))
    DO ia = 1, nat_
      CALL qes_init(site_obj(ia), "SiteMagnetization", SPECIES=atm(ityp(ia)), ATOM=ia, CHARGE = charges(ia), &
                    SitMag = data(1:3,ia))
    END DO 
    CALL qes_init(obj, "Site_Magnetizations", NAT = nat_, SiteMagnetization = site_obj )
    DEALLOCATE (site_obj)
    END SUBROUTINE qexsd_init_d3mags 
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
    REAL(DP),OPTIONAL,INTENT(IN)            :: fermi_energy, ef_updw(:), homo, lumo
    REAL(DP),DIMENSION(:,:),INTENT(IN)      :: et, wg, xk
    REAL(DP),DIMENSION(:),INTENT(IN)        :: wk
    INTEGER,DIMENSION(:),INTENT(IN)         :: ngk      
    TYPE(k_points_IBZ_type),INTENT(IN)      :: starting_kpoints
    TYPE(occupations_type), INTENT(IN)      :: occupations_kind
    TYPE(smearing_type),OPTIONAL,INTENT(IN) :: smearing
    LOGICAL,INTENT(IN)                      :: wf_collected                    
    ! 
    INTEGER                                 :: ndim_ks_energies, ik
    INTEGER,TARGET                          :: nbnd_, nbnd_up_, nbnd_dw_
    INTEGER,POINTER                         :: nbnd_opt, nbnd_up_opt, nbnd_dw_opt 
    TYPE(k_point_type)                      :: kp_obj
    TYPE(ks_energies_type),ALLOCATABLE      :: ks_objs(:)
    TYPE (k_points_IBZ_type)                :: starting_k_points_ 
    REAL(DP),DIMENSION(:),ALLOCATABLE       :: eigenvalues, occupations
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
    starting_k_points_ = starting_kpoints
    starting_k_points_%tagname = "starting_k_points"
    !
    CALL qes_init  (obj, TAGNAME, LSDA = lsda, NONCOLIN = noncolin, SPINORBIT = lspinorb, NBND = nbnd_opt,  &
                   NELEC = nelec, STARTING_K_POINTS = starting_k_points_,      & 
                   NKS = ndim_ks_energies, OCCUPATIONS_KIND = occupations_kind, KS_ENERGIES = ks_objs,      &
                   NBND_UP = nbnd_up_opt, NBND_DW = nbnd_dw_opt,              &
                   FERMI_ENERGY = fermi_energy, HIGHESTOCCUPIEDLEVEL = homo,  TWO_FERMI_ENERGIES = ef_updw, & 
                   SMEARING = smearing, LOWESTUNOCCUPIEDLEVEL = lumo)
    DO ik=1,ndim_ks_energies
       CALL qes_reset(ks_objs(ik))
    END DO
    CALL qes_reset( starting_k_points_ ) 
    DEALLOCATE (ks_objs,eigenvalues, occupations)
    !
    END SUBROUTINE qexsd_init_band_structure 
    ! 
    !---------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_total_energy(obj, etot, eband, ehart, vtxc, etxc, ewald, degauss, demet, &
                       electric_field_corr, potentiostat_contr, gate_contribution, dispersion_contribution, &
                       esol, vsol)
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
    REAL(DP),OPTIONAL               :: esol, vsol
    !
    CHARACTER(LEN=*),PARAMETER      :: TAGNAME="total_energy"
    ! 
    CALL  qes_init (obj, TAGNAME, ETOT = etot, EBAND = eband, EHART = ehart, VTXC = vtxc, ETXC = etxc, & 
                    EWALD = ewald, DEMET = demet, EFIELDCORR = electric_field_corr, POTENTIOSTAT_CONTR = potentiostat_contr,  &
                    GATEFIELD_CONTR = gate_contribution, vdW_term = dispersion_contribution, ESOL = esol, &
                    levelshift_contr = vsol)
    
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
    !-----------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_sawtooth_info(sawtooth_energy, efield_corr, edir, eamp, emaxpos, eopreg) 
    !------------------------------------------------------------------------------------------------
      !
      IMPLICIT NONE
      TYPE(sawtoothEnergy_type), INTENT(OUT)  :: sawtooth_energy 
      REAL(DP),INTENT(IN)                     :: efield_corr, eamp, emaxpos, eopreg 
      INTEGER                                 :: edir 
      call qes_init(sawtooth_energy,TAGNAME="sawtoothEnergy", EDIR=edir, EAMP=eamp, EMAXPOS=emaxpos, & 
                    EOPREG=eopreg, sawtoothEnergy=efield_corr) 
    END SUBROUTINE qexsd_init_sawtooth_info 
    !------------------------------------------------------------------------------------------------
    SUBROUTINE qexsd_init_dipole_info (dipole_info, el_dipole, ion_dipole, edir, eamp, emaxpos, eopreg) 
       !------------------------------------------------------------------------------------------------
       ! 
       USE constants,       ONLY : fpi 
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
                                               ion_pol, sawtooth_obj, dipole_obj , gateInfo)
    !---------------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    ! 
    TYPE(outputElectricField_type)                    :: obj 
    ! 
    LOGICAL,INTENT(IN)                                :: lberry, lelfield, tefield, ldipole
    REAL(DP),OPTIONAL,INTENT(IN)                      :: el_pol(:), ion_pol(:)
    TYPE ( berryPhaseOutput_type),OPTIONAL,INTENT(IN) :: bp_obj
    TYPE ( sawtoothEnergy_type),OPTIONAL,INTENT(IN)   :: sawtooth_obj 
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
                                    sawtoothEnergy = sawtooth_obj, & 
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
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_rism3d(obj, nmol, slabel, molfile, dens1, dens2, ecutsolv)
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(rism3d_type)            :: obj
      INTEGER,          INTENT(IN) :: nmol
      CHARACTER(LEN=*), INTENT(IN) :: slabel(:)
      CHARACTER(LEN=*), INTENT(IN) :: molfile(:)
      REAL(DP),         INTENT(IN) :: dens1(:)
      REAL(DP),         INTENT(IN) :: dens2(:)
      REAL(DP),         INTENT(IN) :: ecutsolv
      !
      TYPE(solvent_type), ALLOCATABLE :: solvents(:)
      INTEGER :: i
      !
      ALLOCATE(solvents(nmol))
      !
      DO i = 1, nmol
          CALL qes_init (solvents(i), "solvent", TRIM(slabel(i)), TRIM(molfile(i)), dens1(i), dens2(i))
      ENDDO
      !
      CALL qes_init (obj, "rism3d", nmol, solvents, ecutsolv)
      !
      DO i = 1, nmol
          CALL qes_reset (solvents(i))
      ENDDO
      !
      DEALLOCATE(solvents)
      !
    END SUBROUTINE qexsd_init_rism3d
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_init_rismlaue(obj, both_hands, nfit, pot_ref, charge, &
                                   right_start,  right_expand, &
                                   right_buffer, right_buffer_u, right_buffer_v, &
                                   left_start,   left_expand,  &
                                   left_buffer,  left_buffer_u,  left_buffer_v)
      !------------------------------------------------------------------------
      IMPLICIT NONE
      !
      TYPE(rismlaue_type)  :: obj
      LOGICAL,  INTENT(IN) :: both_hands
      INTEGER,  INTENT(IN) :: nfit
      INTEGER,  INTENT(IN) :: pot_ref
      REAL(DP), INTENT(IN) :: charge
      REAL(DP), INTENT(IN) :: right_start
      REAL(DP), INTENT(IN) :: right_expand
      REAL(DP), INTENT(IN) :: right_buffer
      REAL(DP), INTENT(IN) :: right_buffer_u
      REAL(DP), INTENT(IN) :: right_buffer_v
      REAL(DP), INTENT(IN) :: left_start
      REAL(DP), INTENT(IN) :: left_expand
      REAL(DP), INTENT(IN) :: left_buffer
      REAL(DP), INTENT(IN) :: left_buffer_u
      REAL(DP), INTENT(IN) :: left_buffer_v
      !
      CALL qes_init (obj, "rismlaue", both_hands, nfit, pot_ref, charge, &
                     right_start,  right_expand, &
                     right_buffer, right_buffer_u, right_buffer_v, &
                     left_start,   left_expand,  &
                     left_buffer,  left_buffer_u,  left_buffer_v)
      !
    END SUBROUTINE qexsd_init_rismlaue

 END MODULE qexsd_init
