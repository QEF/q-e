! Copyright (C) 2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE qexsd_copy
  !----------------------------------------------------------------------------
  !! This module contains some common subroutines used to copy data read from
  !! XML format into data used by the Quantum ESPRESSO package.
  !
  !! Written by Paolo Giannozzi, building upon pre-existing code.
  !
  USE kinds, ONLY : dp
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  PUBLIC:: qexsd_copy_geninfo, qexsd_copy_parallel_info, &
       qexsd_copy_atomic_species, qexsd_copy_atomic_structure, &
       qexsd_copy_symmetry, qexsd_copy_algorithmic_info, &
       qexsd_copy_basis_set, qexsd_copy_dft, qexsd_copy_band_structure, &
       qexsd_copy_efield, qexsd_copy_magnetization, qexsd_copy_kpoints, copy_order_um_from_xml, &
       qexsd_copy_efermi, qexsd_copy_rism3d, qexsd_copy_rismlaue, qexsd_copy_esm, qexsd_copy_twochem 
 !
 INTERFACE 
  SUBROUTINE copy_order_um_from_xml(dftU_obj, nsp, nat, Hubbard_lmax, order_um)
    !! routine that copies order_um values from xml file to allocatable order_um, used 
    !! for the initialization of the order_um variable in ldaU module. 
    !FIXME to be collected with other related copying from dftU object or better read from the scf file
     USE qes_types_module, only: dftU_type
     IMPLICIT NONE
     TYPE(dftU_type),INTENT(IN)  :: dftU_obj
     !! object containing the dftU infortmation read from the XML file
     INTEGER,INTENT(IN)          :: nsp, nat, Hubbard_lmax
     !! number of species defined for this calculation 
     !! number of atoms for this calculation
     !! max l value for species in this calculatio 
     INTEGER, ALLOCATABLE,INTENT(INOUT) :: order_um(:,:,:)  
     !! allocatable 3-rank array to be written dimensions will be 2*Hubbarl_lmax+1, nspin, nat  
  END SUBROUTINE copy_order_um_from_xml 
END INTERFACE
CONTAINS
  !-------------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_geninfo (geninfo_obj, qexsd_fmt, qexsd_version)
    !-------------------------------------------------------------------------------
    ! 
    USE io_files,            ONLY: qexsd_init
    USE qes_types_module,    ONLY: general_info_type
    IMPLICIT NONE 
    !
    TYPE (general_info_type ),INTENT(IN)  :: geninfo_obj   
    CHARACTER(LEN=*), INTENT(OUT) :: qexsd_fmt, qexsd_version
    ! 
    IF ( qexsd_init ) RETURN 
    qexsd_fmt = TRIM (geninfo_obj%xml_format%NAME)
    qexsd_version = TRIM ( geninfo_obj%xml_format%VERSION)
    qexsd_init = .TRUE. 
    !
  END SUBROUTINE qexsd_copy_geninfo
  ! 
  !
  !---------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_parallel_info (parinfo_obj, nproc_file, &
       nproc_pool_file, nproc_image_file, ntask_groups_file, &
       nproc_bgrp_file, nproc_ortho_file)
    !---------------------------------------------------------------------------    !
    USE qes_types_module, ONLY : parallel_info_type
    !
    IMPLICIT NONE 
    !
    TYPE ( parallel_info_type ),INTENT(IN)     :: parinfo_obj
    INTEGER, INTENT(OUT) :: nproc_file, nproc_pool_file, &
                            nproc_image_file, ntask_groups_file, &
                            nproc_bgrp_file, nproc_ortho_file
    ! 
    nproc_file = parinfo_obj%nprocs
    nproc_pool_file = nproc_file/parinfo_obj%npool
    nproc_image_file = nproc_file 
    ntask_groups_file = parinfo_obj%ntasks
    nproc_bgrp_file = nproc_image_file / parinfo_obj%npool / parinfo_obj%nbgrp 
    nproc_ortho_file = parinfo_obj%ndiag
    !
  END SUBROUTINE qexsd_copy_parallel_info
  !
  !--------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_atomic_species (atomic_species, nsp, atm, amass, &
       starting_magnetization, angle1, angle2, psfile, pseudo_dir)
    !---------------------------------------------------------------------------    !
    USE qes_types_module, ONLY : atomic_species_type
    !
    IMPLICIT NONE 
    !
    TYPE ( atomic_species_type ),INTENT(IN)    :: atomic_species
    INTEGER, INTENT(out) :: nsp
    CHARACTER(LEN=*), INTENT(out) :: atm(:)
    REAL(dp), INTENT(out) :: amass(:)
    REAL(dp), OPTIONAL, INTENT(out) :: starting_magnetization(:), &
         angle1(:), angle2(:)
    CHARACTER(LEN=*), OPTIONAL, INTENT(out) :: psfile(:), pseudo_dir
    !
    INTEGER :: isp
    !
    nsp = atomic_species%ntyp
    DO isp = 1, nsp 
       amass(isp) = 0.d0 
       IF (atomic_species%species(isp)%mass_ispresent) &
            amass(isp) = atomic_species%species(isp)%mass
       atm(isp) = TRIM ( atomic_species%species(isp)%name )
       IF ( PRESENT (psfile) ) THEN
          psfile(isp) = TRIM ( atomic_species%species(isp)%pseudo_file) 
       END IF
       IF ( PRESENT (starting_magnetization) ) THEN
          IF ( atomic_species%species(isp)%starting_magnetization_ispresent) THEN
             starting_magnetization(isp) = atomic_species%species(isp)%starting_magnetization
          END IF
       END IF
       IF ( PRESENT (angle1) ) THEN
          IF ( atomic_species%species(isp)%spin_teta_ispresent ) THEN 
             angle1(isp) =  atomic_species%species(isp)%spin_teta 
          END IF
       END IF
       IF ( PRESENT (angle2) ) THEN
          IF ( atomic_species%species(isp)%spin_phi_ispresent ) THEN
             angle2(isp) = atomic_species%species(isp)%spin_phi
          END IF
       END IF
    END DO
    ! 
    ! ... this is where PP files were originally found (if available)
    !
    IF ( PRESENT (pseudo_dir) ) THEN
       IF ( atomic_species%pseudo_dir_ispresent) THEN 
          pseudo_dir = TRIM(atomic_species%pseudo_dir)
       ELSE 
          pseudo_dir = ' '
       END IF
    END IF
    !
  END SUBROUTINE qexsd_copy_atomic_species

  !--------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_atomic_structure (atomic_structure, nsp, atm, &
       nat, tau, ityp, alat, a1, a2, a3, ibrav, natomwfc )
  !--------------------------------------------------------------------------
    
    USE qes_types_module, ONLY : atomic_structure_type
    !
    IMPLICIT NONE 
    !
    TYPE ( atomic_structure_type ),INTENT(IN)  :: atomic_structure
    INTEGER, INTENT(in) :: nsp 
    CHARACTER(LEN = 6), INTENT(in) :: atm(:)
    !
    INTEGER, INTENT(out)  :: nat, ibrav, natomwfc
    REAL(dp), INTENT(out) :: alat, a1(:), a2(:), a3(:)
    INTEGER, INTENT(inout),  ALLOCATABLE :: ityp(:)
    REAL(dp), INTENT(inout), ALLOCATABLE :: tau(:,:)
    !
    CHARACTER(LEN=3), ALLOCATABLE :: symbols(:)
    INTEGER :: iat, idx, isp
    !
    nat = atomic_structure%nat 
    IF ( atomic_structure%num_of_atomic_wfc_ispresent ) THEN 
      natomwfc = atomic_structure%num_of_atomic_wfc
    ELSE 
      natomwfc = 0 
    END IF   
    alat = atomic_structure%alat 
    IF ( atomic_structure%bravais_index_ispresent ) THEN 
       ibrav = atomic_structure%bravais_index 
       IF (atomic_structure%alternative_axes_ispresent ) THEN 
         SELECT CASE(ibrav) 
            CASE(3)
               IF (TRIM(atomic_structure%alternative_axes)=="b:a-b+c:-c") THEN 
                  ibrav = -ibrav
               ELSE 
                  CALL errore("qexsd_copy_atomic_structure:","alternative axes not recognised", 1) 
               END IF
            CASE(5) 
               IF (TRIM(atomic_structure%alternative_axes)=="3fold-111") THEN
                    ibrav = -ibrav
               ELSE
                    CALL errore("qexsd_copy_atomic_structure:","alternative axes not recognised", 1)
               END IF
            CASE(9)
                IF (TRIM(atomic_structure%alternative_axes)=="-b:a:c") THEN
                      ibrav = -ibrav
                ELSE IF( TRIM(atomic_structure%alternative_axes)=="bcoA-type") THEN 
                     ibrav = 91
                ELSE
                      CALL errore("qexsd_copy_atomic_structure:","alternative axes not recognised", 1)
                END IF
            CASE(13,12) 
                IF (TRIM(atomic_structure%alternative_axes)=="unique-axis-b") THEN
                      ibrav = -ibrav
                 ELSE
                      CALL errore("qexsd_copy_atomic_structure:","alternativ axes not recognised", 1)
                 END IF
            END SELECT
       END IF 
    ELSE 
       ibrav = 0
    END IF
    IF ( .NOT. ALLOCATED(tau) ) ALLOCATE(tau(3,nat))
    IF ( .NOT. ALLOCATED(ityp)) ALLOCATE(ityp(nat))
    ALLOCATE ( symbols(nat) )
    loop_on_atoms:DO iat = 1, nat
       idx = atomic_structure%atomic_positions%atom(iat)%index
       tau(:,idx) = atomic_structure%atomic_positions%atom(iat)%atom 
       symbols(idx)  = TRIM ( atomic_structure%atomic_positions%atom(idx)%name)
       loop_on_species:DO isp = 1, nsp
          IF ( TRIM(symbols(idx)) == TRIM (atm(isp))) THEN 
             ityp(iat) = isp 
             exit loop_on_species
          END IF 
       END  DO loop_on_species
    END DO loop_on_atoms
    DEALLOCATE (symbols)
    IF ( atomic_structure%alat_ispresent ) alat = atomic_structure%alat 
    a1(:) = atomic_structure%cell%a1
    a2(:) = atomic_structure%cell%a2
    a3(:) = atomic_structure%cell%a3

  END SUBROUTINE qexsd_copy_atomic_structure
  !
  !------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_symmetry ( symms_obj, spacegroup, &
       nsym, nrot, s, ft, sname, t_rev, invsym, irt, &
       noinv, nosym, no_t_rev, colin_mag, flags_obj)
    !------------------------------------------------------------------------
    ! 
    USE qes_types_module,ONLY : symmetries_type, symmetry_flags_type
    ! 
    IMPLICIT NONE   
    ! 
    TYPE ( symmetries_type )             :: symms_obj 
    TYPE (symmetry_flags_type), OPTIONAL :: flags_obj
    INTEGER, INTENT(OUT) :: spacegroup
    INTEGER, INTENT(OUT) :: nrot
    INTEGER, INTENT(OUT) :: nsym
    INTEGER, INTENT(OUT) :: s(:,:,:)
    LOGICAL, INTENT(OUT) :: invsym
    REAL(dp), INTENT(OUT):: ft(:,:)
    INTEGER, INTENT(OUT) :: irt(:,:)
    INTEGER, INTENT(OUT) :: t_rev(:)
    INTEGER, INTENT(OUT) :: colin_mag
    CHARACTER(len=45) ::  sname(:)
    !
    LOGICAL, INTENT(OUT) :: noinv, nosym, no_t_rev
    !
    INTEGER :: isym 
    ! 
    IF ( PRESENT(flags_obj) ) THEN 
       noinv = flags_obj%noinv
       nosym = flags_obj%nosym
       no_t_rev = flags_obj%no_t_rev
    ELSE
       noinv = .FALSE.
       nosym = .FALSE.
       no_t_rev=.FALSE.
    ENDIF
    !
    spacegroup = symms_obj%space_group
    nrot = symms_obj%nrot 
    nsym = symms_obj%nsym
    IF (symms_obj%colin_mag_ispresent) THEN 
      colin_mag = symms_obj%colin_mag 
    ELSE 
      colin_mag = -1 
    END IF 
    !  
    invsym = .FALSE. 
    DO isym = 1, nrot
       s(:,:,isym) = reshape(symms_obj%symmetry(isym)%rotation%matrix, [3,3]) 
       sname(isym) = TRIM ( symms_obj%symmetry(isym)%info%name )  
       IF ( (TRIM(sname(isym)) == "inversion") .AND. (isym .LE. nsym) ) invsym = .TRUE.
       IF ( symms_obj%symmetry(isym)%fractional_translation_ispresent .AND. (isym .LE. nsym) ) THEN
          ft(1:3,isym)  =  symms_obj%symmetry(isym)%fractional_translation(1:3) 
       END IF
       IF ( symms_obj%symmetry(isym)%info%time_reversal_ispresent ) THEN  
          IF (symms_obj%symmetry(isym)%info%time_reversal) THEN 
             t_rev( isym ) = 1
          ELSE
             t_rev( isym ) = 0 
          END IF
       END IF
       IF ( symms_obj%symmetry(isym)%equivalent_atoms_ispresent .AND. (isym .LE. nsym) )   &
            irt(isym,:) = symms_obj%symmetry(isym)%equivalent_atoms%equivalent_atoms(:)
    END DO
    !
  END SUBROUTINE qexsd_copy_symmetry
  !

  !--------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_basis_set ( basis_set, gamma_only, ecutwfc, ecutrho, &
       nr1s, nr2s, nr3s, nr1, nr2, nr3, nr1b, nr2b, nr3b, &
       ngm_g, ngms_g, npw_g, b1, b2, b3 )
  !--------------------------------------------------------------------------   
    !
    USE qes_types_module, ONLY : basis_set_type
    !
    IMPLICIT NONE 
    TYPE ( basis_set_type ),INTENT(IN)         :: basis_set
    LOGICAL, INTENT(out)  :: gamma_only
    INTEGER, INTENT(out)  :: ngm_g, ngms_g, npw_g
    INTEGER, INTENT(out)  :: nr1s, nr2s, nr3s, nr1, nr2, nr3
    INTEGER, INTENT(inout)  :: nr1b, nr2b, nr3b
    REAL(dp), INTENT(out) :: ecutwfc, ecutrho, b1(:), b2(:), b3(:)
    ! 
    ecutwfc = basis_set%ecutwfc
    ecutrho = basis_set%ecutrho
    gamma_only= basis_set%gamma_only
    nr1 = basis_set%fft_grid%nr1
    nr2 = basis_set%fft_grid%nr2          
    nr3 = basis_set%fft_grid%nr3
    nr1s= basis_set%fft_smooth%nr1
    nr2s= basis_set%fft_smooth%nr2
    nr3s= basis_set%fft_smooth%nr3
    IF ( basis_set%fft_box_ispresent ) THEN
       nr1b = basis_set%fft_box%nr1
       nr2b = basis_set%fft_box%nr2
       nr3b = basis_set%fft_box%nr3
    END IF
    ngm_g     = basis_set%ngm
    ngms_g    = basis_set%ngms
    npw_g     = basis_set%npwx
    !
    b1 =  basis_set%reciprocal_lattice%b1
    b2 =  basis_set%reciprocal_lattice%b2
    b3 =  basis_set%reciprocal_lattice%b3
    !
  END SUBROUTINE qexsd_copy_basis_set
  !
  !-----------------------------------------------------------------------
  SUBROUTINE qexsd_copy_dft ( dft_obj, nsp, atm, &
       dft_name, nq1, nq2, nq3, ecutfock, exx_fraction, screening_parameter, &
       exxdiv_treatment, x_gamma_extrapolation, ecutvcut, local_thr, &
       lda_plus_U, apply_u, lda_plus_U_kind, U_projection, Hubbard_n, Hubbard_l, Hubbard_lmax, Hubbard_occ, &
       Hubbard_n2, Hubbard_l2, Hubbard_n3, Hubbard_l3, backall, Hubbard_lmax_back, Hubbard_alpha_back, &
       Hubbard_U, Hubbard_Um, Hubbard_U2, Hubbard_J0, Hubbard_alpha, Hubbard_alpha_m, Hubbard_beta, Hubbard_J, Hubbard_V, &
       vdw_corr, dftd3_version, dftd3_3body, scal6, lon_rcut, vdw_isolated )
    !-------------------------------------------------------------------
    ! 
    USE qes_types_module, ONLY : dft_type
    USE upf_utils,        ONLY : spdf_to_l
    !
    IMPLICIT NONE 
    TYPE ( dft_type ),INTENT(in) :: dft_obj
    INTEGER, INTENT(in)          :: nsp 
    CHARACTER(LEN=*), INTENT(in) :: atm(nsp)
    ! 
    CHARACTER(LEN=*), INTENT(out) :: dft_name
    ! Variables that may or may not be present should be intent(inout)
    ! so that they do not forget their default value (if any)
    CHARACTER(LEN=*), INTENT(inout) :: exxdiv_treatment
    REAL(dp), INTENT(inout) :: ecutfock, exx_fraction, screening_parameter, &
         ecutvcut, local_thr
    INTEGER, INTENT(inout) :: nq1, nq2, nq3
    LOGICAL, INTENT(inout) :: x_gamma_extrapolation
    !
    LOGICAL, INTENT(out) :: lda_plus_U, apply_u
    INTEGER, INTENT(inout) :: lda_plus_U_kind, Hubbard_lmax, Hubbard_lmax_back
    CHARACTER(LEN=*), INTENT(inout) :: U_projection
    INTEGER, INTENT(inout) :: Hubbard_n(:), Hubbard_l(:), Hubbard_n2(:), Hubbard_l2(:), Hubbard_n3(:), Hubbard_l3(:) 
    REAL(dp), INTENT(inout) :: Hubbard_U(:), Hubbard_U2(:), Hubbard_J0(:), Hubbard_J(:,:), Hubbard_V(:,:,:), &
                               Hubbard_alpha(:), Hubbard_alpha_back(:), Hubbard_beta(:), Hubbard_occ(:,:),   & 
                               Hubbard_Um(:,:,:), Hubbard_alpha_m(:,:,:) 
    LOGICAL, INTENT(inout) :: backall(:)
    OPTIONAL    :: Hubbard_U2, Hubbard_n2, Hubbard_l2, Hubbard_lmax_back, Hubbard_alpha_back, &
                   Hubbard_l3, Hubbard_Um
    !
    CHARACTER(LEN=*), INTENT(out) :: vdw_corr
    LOGICAL,INTENT(inout)   :: dftd3_3body 
    INTEGER,INTENT(inout)   :: dftd3_version 
    REAL(dp), INTENT(inout) :: scal6, lon_rcut
    LOGICAL, INTENT(inout) :: vdw_isolated
    !
    CHARACTER(LEN=256 ) :: label
    CHARACTER(LEN=3 )   :: symbol
    INTEGER :: ihub, isp, hu_n, hu_l, idx1, idx2, idx3, ich, im, ispin, objldim
    !
    dft_name = TRIM(dft_obj%functional)
    IF ( dft_obj%hybrid_ispresent ) THEN
       nq1 = dft_obj%hybrid%qpoint_grid%nqx1
       nq2 = dft_obj%hybrid%qpoint_grid%nqx2
       nq3 = dft_obj%hybrid%qpoint_grid%nqx3
       ecutfock = dft_obj%hybrid%ecutfock
       exx_fraction = dft_obj%hybrid%exx_fraction
       screening_parameter = dft_obj%hybrid%screening_parameter
       exxdiv_treatment = dft_obj%hybrid%exxdiv_treatment
       x_gamma_extrapolation = dft_obj%hybrid%x_gamma_extrapolation
       ecutvcut = dft_obj%hybrid%ecutvcut
       IF (dft_obj%hybrid%localization_threshold_ispresent) THEN
          local_thr = dft_obj%hybrid%localization_threshold  
       ELSE 
          local_thr = 0._DP 
      END IF 
    END IF
    !
    lda_plus_u = dft_obj%dftU_ispresent 
    IF ( lda_plus_u ) THEN 
       Hubbard_U  = 0.0_DP
       Hubbard_U2 = 0.0_DP
       Hubbard_Um = 0.0_DP
       Hubbard_alpha = 0.0_DP
       Hubbard_alpha_back = 0.0_DP
       Hubbard_J = 0.0_DP
       Hubbard_J0 = 0.0_DP
       Hubbard_beta = 0.0_DP
       Hubbard_V  = 0.0_DP 
       lda_plus_u_kind = dft_obj%dftU%lda_plus_u_kind
       U_projection = TRIM ( dft_obj%dftU%U_projection_type )
       Hubbard_n  =-1 
       Hubbard_l  =-1 
       Hubbard_n2 =-1 
       Hubbard_l2 =-1 
       Hubbard_l3 =-1 
       backall = .false.
       !FIXME read and write of Hubbard_alpha_m from XML not yet implemented, temporarily 
       ! we just set it to 0 not to break other cases 
       Hubbard_alpha_m = 0.0_DP
       !
       IF ( dft_obj%dftU%Hubbard_U_ispresent) THEN 
          loop_on_hubbardU:DO ihub =1, dft_obj%dftU%ndim_Hubbard_U
             symbol = TRIM(dft_obj%dftU%Hubbard_U(ihub)%specie)
             label  = TRIM(dft_obj%dftU%Hubbard_U(ihub)%label ) 
             loop_on_speciesU:DO isp = 1, nsp
                IF ( TRIM(symbol) == TRIM ( atm(isp) ) ) THEN 
                     Hubbard_U(isp) = dft_obj%dftU%Hubbard_U(ihub)%HubbardCommon
                     READ (label(1:1),'(i1)', END=14, ERR=15) hu_n
                     hu_l = spdf_to_l( label(2:2) )
                     Hubbard_n(isp) = hu_n
                     Hubbard_l(isp) = hu_l
                     IF (Hubbard_n(isp)<0 .OR. Hubbard_l(isp)<0) &
                        CALL errore ("qexsd_copy_dft:", &
                            &"Problem while reading Hubbard_n and/or Hubbard_l", 1 )
                     EXIT loop_on_speciesU
                END IF 
             END DO loop_on_speciesU
          END DO loop_on_hubbardU
       END IF 
       ! 
       Hubbard_Um = 0.0_dp
       IF ( dft_obj%dftU%Hubbard_Um_ispresent) THEN 
          apply_u = .TRUE. 
          loop_on_hubbardUm:DO ihub =1, dft_obj%dftU%ndim_Hubbard_Um
             symbol = TRIM(dft_obj%dftU%Hubbard_Um(ihub)%specie)
             label  = TRIM(dft_obj%dftU%Hubbard_Um(ihub)%label )  
             IF (dft_obj%dftU%Hubbard_Um(ihub)%spin_ispresent) THEN
                     ispin = dft_obj%dftU%Hubbard_Um(ihub)%spin
             ELSE 
                     ispin = 1
             END IF
             loop_on_speciesUm:DO isp = 1, nsp
                IF ( TRIM(symbol) == TRIM ( atm(isp) ) ) THEN 
                  READ (label(1:1),'(i1)', END=14, ERR=15) hu_n
                  hu_l = spdf_to_l( label(2:2) )
                  Hubbard_n(isp) = hu_n
                  Hubbard_l(isp) = hu_l
                  objldim = dft_obj%dftU%Hubbard_Um(ihub)%size 
                  IF (objldim == 2*hu_l + 1)  THEN 
                    Hubbard_Um(1:2*hu_l+1,ispin,isp) = dft_obj%dftU%Hubbard_Um(ihub)%HubbardM 
                  ELSE IF (objldim == 2*(2*hu_l +1)) THEN 
                    Hubbard_Um(1:2*hu_l+1,1,isp) = dft_obj%dftU%Hubbard_Um(ihub)%HubbardM(1:2*hu_l+1)
                    Hubbard_Um(1:2*hu_l+1,2,isp) = dft_obj%dftU%Hubbard_Um(ihub)%HubbardM(2*hu_l+2:4*hu_l+2)
                  ELSE 
                    call errore("qexsd_copy_dft:", & 
                       "size of Hubbard_Um element not compatible with label",1) 
                  END IF  
                  IF (Hubbard_n(isp)<0 .OR. Hubbard_l(isp)<0) &
                        CALL errore ("qexsd_copy_dft:", &
                            &"Problem while reading Hubbard_n and/or Hubbard_l", 1 )
                  EXIT loop_on_speciesUm
                END IF 
             END DO loop_on_speciesUm
          END DO loop_on_hubbardUm
       END IF 

       !
       IF ( dft_obj%dftU%Hubbard_back_ispresent) THEN
          loop_hubbardBack:DO ihub =1, dft_obj%dftU%ndim_Hubbard_back
            symbol = TRIM(dft_obj%dftU%Hubbard_back(ihub)%species)
            loop_on_species_back:DO isp = 1, nsp
              IF ( TRIM(symbol) == TRIM ( atm(isp) ) ) THEN
                Hubbard_U2(isp) = dft_obj%dftU%Hubbard_back(ihub)%Hubbard_U2
                Hubbard_l2(isp) = dft_obj%dftU%Hubbard_back(ihub)%l2_number 
                Hubbard_n2 =      dft_obj%dftU%Hubbard_back(ihub)%n2_number
                SELECT  CASE (TRIM (dft_obj%dftU%Hubbard_back(ihub)%background))
                  CASE ('one_orbital')
                    backall(isp) = .FALSE. 
                  CASE ('two_orbitals')
                    backall(isp) = .TRUE.
                    IF ( .NOT. dft_obj%dftU%Hubbard_back(ihub)%l3_number_ispresent) &
                      CALL errore ("qexsd_copy_dft", "Only 1 l number found for 2 orbitals Hubbard Background", 1)
                    Hubbard_l3(isp) = dft_obj%dftU%Hubbard_back(ihub)%l3_number
                    IF (dft_obj%dftU%Hubbard_back(ihub)%n3_number_ispresent) & 
                      Hubbard_n3 = dft_obj%dftU%Hubbard_back(ihub)%n3_number 
                  CASE DEFAULT
                    CALL errore ("qexsd_copy_dft", "Unrecognized Background type", 1)
                END SELECT 
                EXIT loop_on_species_back
              END IF
              END DO loop_on_species_back
            END DO loop_hubbardBack
       END IF
       !
      IF (dft_obj%dftU%Hubbard_Occ_ispresent) THEN 
         loop_on_hubbard_occ: DO ihub =1, dft_obj%dftU%ndim_Hubbard_Occ 
            symbol = TRIM(dft_obj%dftU%Hubbard_Occ(ihub)%specie) 
            loop_on_species: DO isp = 1, nsp
               IF (TRIM(symbol) == TRIM(atm(isp))) THEN 
                  DO ich = 1, dft_obj%dftU%Hubbard_Occ(ihub)%channels 
                     Hubbard_occ(isp,ich) = dft_obj%dftU%Hubbard_Occ(ihub)%channel_occ(ich)%ChannelOcc 
                  END DO 
               END IF 
            END DO loop_on_species 
         END DO loop_on_hubbard_occ
      END IF 
      ! 
       IF ( dft_obj%dftU%Hubbard_J0_ispresent ) THEN 
            loop_on_hubbardj0:DO ihub =1, dft_obj%dftU%ndim_Hubbard_J0
               symbol = TRIM(dft_obj%dftU%Hubbard_J0(ihub)%specie)
               loop_on_speciesj0:DO isp = 1, nsp
                  IF ( TRIM(symbol) == TRIM (atm(isp)) ) THEN
                     Hubbard_J0(isp) = dft_obj%dftU%Hubbard_J0(ihub)%HubbardCommon
                     EXIT loop_on_speciesj0
                  END IF
               END DO loop_on_speciesj0
            END DO loop_on_hubbardj0
       END IF
       !
       IF ( dft_obj%dftU%Hubbard_alpha_ispresent) THEN 
            loop_on_hubbardAlpha:DO ihub =1, dft_obj%dftU%ndim_Hubbard_alpha
               symbol = TRIM(dft_obj%dftU%Hubbard_alpha(ihub)%specie)
               loop_on_speciesAlpha:DO isp = 1, nsp
                  IF ( TRIM(symbol) == TRIM (atm(isp)) ) THEN 
                     Hubbard_alpha(isp) = dft_obj%dftU%Hubbard_alpha(ihub)%HubbardCommon
                     EXIT loop_on_speciesAlpha
                  END IF
               END DO loop_on_speciesAlpha
            END DO loop_on_hubbardAlpha
       END IF
       !
       IF ( dft_obj%dftU%Hubbard_alpha_back_ispresent) THEN    
            loop_on_hubbardAlphaBack:DO ihub =1, dft_obj%dftU%ndim_Hubbard_alpha_back
               symbol = TRIM(dft_obj%dftU%Hubbard_alpha_back(ihub)%specie)
               loop_on_speciesAlphaBack:DO isp = 1, nsp
                  IF ( TRIM(symbol) == TRIM (atm(isp)) ) THEN
                     Hubbard_alpha_back(isp) = dft_obj%dftU%Hubbard_alpha_back(ihub)%HubbardCommon
                     EXIT loop_on_speciesAlphaBack
                  END IF
               END DO loop_on_speciesAlphaBack
            END DO loop_on_hubbardAlphaBack
       END IF
       !
       IF ( dft_obj%dftU%Hubbard_beta_ispresent) THEN 
            loop_on_hubbardBeta:DO ihub =1, dft_obj%dftU%ndim_Hubbard_beta
               symbol = TRIM(dft_obj%dftU%Hubbard_beta(ihub)%specie)
               loop_on_speciesBeta:DO isp = 1, nsp
                  IF ( TRIM(symbol) == TRIM (atm(isp)) ) THEN 
                     Hubbard_beta(isp) = dft_obj%dftU%Hubbard_beta(ihub)%HubbardCommon
                     EXIT loop_on_speciesBeta
                  END IF
               END DO loop_on_speciesBeta
            END DO loop_on_hubbardBeta
       END IF
       !
       IF ( dft_obj%dftU%Hubbard_J_ispresent) THEN 
            loop_on_hubbardJ:DO ihub =1, dft_obj%dftU%ndim_Hubbard_J
               symbol = TRIM(dft_obj%dftU%Hubbard_J(ihub)%specie)
               loop_on_speciesJ:DO isp = 1, nsp
                  IF ( TRIM(symbol) == TRIM (atm(isp)) ) THEN 
                     Hubbard_J(:,isp) = dft_obj%dftU%Hubbard_J(ihub)%HubbardJ
                     EXIT loop_on_speciesJ
                  END IF
               END DO loop_on_speciesJ
            END DO loop_on_hubbardJ
       END IF
       !
       IF (dft_obj%dftU%Hubbard_V_ispresent .AND. SIZE (hubbard_v,1) > 1) THEN 
         DO ihub = 1, dft_obj%dftU%ndim_Hubbard_V 
           idx1 = dft_obj%dftU%Hubbard_V(ihub)%index1
           idx2 = dft_obj%dftU%Hubbard_V(ihub)%index2
           IF (Hubbard_V(idx1, idx2,1 ) == 0._DP ) THEN 
             idx3 = 1 
           ELSE IF (Hubbard_V(idx1, idx2, 2) == 0._DP) THEN 
             idx3 = 2 
           ELSE IF (Hubbard_V(idx1, idx2, 3) == 0._DP) THEN 
             idx3 = 3 
           ELSE IF (Hubbard_V(idx1, idx2, 4) == 0._DP) THEN
             idx3 = 4  
           END IF 
           Hubbard_V(idx1, idx2, idx3 ) = dft_obj%dftU%Hubbard_V(ihub)%HubbardInterSpecieV
           symbol = TRIM(dft_obj%dftU%Hubbard_V(ihub)%specie1) 
           label  = TRIM(dft_obj%dftU%hubbard_V(ihub)%label1) 
           DO isp = 1, nsp
             IF (TRIM(symbol) == TRIM(atm(isp)) .AND. & 
                  ( Hubbard_n(isp) == -1 .OR. Hubbard_n2(isp) == -1 ))  THEN 
               READ (label(1:1),'(i1)', END=14, ERR=15) hu_n
               hu_l = spdf_to_l( label(2:2) )
               IF ( idx3 == 1 .OR. idx3 == 2 ) THEN 
                 Hubbard_n(isp) = hu_n
                 Hubbard_l(isp) = hu_l
                 IF (Hubbard_n(isp)<0 .OR. Hubbard_l(isp)<0) &
                    CALL errore ("qexsd_copy_dft:", "Problem while reading Hubbard_n and/or Hubbard_l", 1)
               ELSE IF ( idx3 == 3 .OR. idx3 == 4 ) THEN 
                 Hubbard_n2 = hu_n 
                 Hubbard_l2 = hu_l 
               END IF 
             END IF 
           END DO
         END DO     
       END IF  
       !
       Hubbard_lmax      = MAXVAL( Hubbard_l(1:nsp) )
       Hubbard_lmax_back = MAXVAL( Hubbard_l2(1:nsp) ) 
       ! IT: What about Hubbard_l3?
       !  
    END IF

      IF ( dft_obj%vdW_ispresent ) THEN 
         vdw_corr = TRIM( dft_obj%vdW%vdw_corr ) 
      ELSE
         vdw_corr = ''
      END IF
      
      IF ( dft_obj%vdW_ispresent ) THEN
         IF (dft_obj%vdW%dftd3_threebody_ispresent) THEN 
            dftd3_3body = dft_obj%vdw%dftd3_threebody 
         END IF 
         IF (dft_obj%vdw%dftd3_version_ispresent) THEN 
            dftd3_version = dft_obj%vdW%dftd3_version 
         END IF  
         IF (dft_obj%vdW%london_s6_ispresent ) THEN 
            scal6 = dft_obj%vdW%london_s6
         END IF
         IF ( dft_obj%vdW%london_rcut_ispresent ) THEN 
            lon_rcut = dft_obj%vdW%london_rcut
         END IF
         IF (dft_obj%vdW%ts_vdW_isolated_ispresent ) THEN 
            vdW_isolated = dft_obj%vdW%ts_vdW_isolated
         END IF
      END IF

      RETURN

14    CALL errore ('qexsd_copy_dft:', ' End of file while parsing Hubbard manifolds', 1)
15    CALL errore ('qexsd_copy_dft:', ' Error while parsing Hubbard manifolds', 1)

    END SUBROUTINE qexsd_copy_dft
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_copy_band_structure( band_struct_obj, lsda, nkstot, &
         isk, nbnd, nbnd_up, nbnd_dw, nelec, xk, wk, wg, &
         ef, ef_up, ef_dw, et )
      !------------------------------------------------------------------------
      !
      ! IMPORTANT NOTICE: IN LSDA CASE CONVERTS TO "PWSCF" LOGIC for k-points
      !
      USE qes_types_module, ONLY : band_structure_type
      !
      IMPLICIT NONE
      TYPE ( band_structure_type)         :: band_struct_obj
      LOGICAL, INTENT(out) :: lsda
      INTEGER, INTENT(out) :: nkstot, nbnd, nbnd_up, nbnd_dw, &
              isk(:)
      REAL(dp), INTENT(out):: nelec, ef, ef_up, ef_dw, xk(:,:), wk(:)
      REAL(dp), INTENT(inout), ALLOCATABLE ::  wg(:,:), et(:,:)
      !
      LOGICAL :: two_fermi_energies
      INTEGER :: ik
      ! 
      lsda = band_struct_obj%lsda
      nkstot = band_struct_obj%nks  
      !
      IF ( lsda) THEN
         !
         IF (band_struct_obj%nbnd_ispresent) THEN 
            nbnd  = band_struct_obj%nbnd / 2
         ELSE IF ( band_struct_obj%nbnd_up_ispresent .AND. band_struct_obj%nbnd_dw_ispresent ) THEN 
            nbnd = (band_struct_obj%nbnd_up + band_struct_obj%nbnd_dw)/2 
         ELSE 
            CALL errore ('qexsd_copy_band_structure: ','both nbnd and nbnd_up+nbnd_dw missing', 1)  
         END IF
         !
         ! the quantity used below (not sure about the logic ...):
         !    band_struct_obj%ks_energies(ik)%eigenvalues%size
         ! should be the same for all k-points so ik=1 does the job
         ik = 1
         IF ( band_struct_obj%nbnd_up_ispresent .AND. &
              band_struct_obj%nbnd_dw_ispresent ) THEN
            nbnd_up = band_struct_obj%nbnd_up
            nbnd_dw = band_struct_obj%nbnd_dw 
         ELSE IF ( band_struct_obj%nbnd_up_ispresent ) THEN 
            nbnd_up = band_struct_obj%nbnd_up
            nbnd_dw = band_struct_obj%ks_energies(ik)%eigenvalues%size - nbnd_up
         ELSE IF ( band_struct_obj%nbnd_dw_ispresent ) THEN 
            nbnd_dw = band_struct_obj%nbnd_dw
            nbnd_up = band_struct_obj%ks_energies(ik)%eigenvalues%size - nbnd_dw
         ELSE 
            nbnd_up = band_struct_obj%ks_energies(ik)%eigenvalues%size/2  
            nbnd_dw = band_struct_obj%ks_energies(ik)%eigenvalues%size/2
         END IF
         !
         nkstot = nkstot * 2 
         isk(1:nkstot/2) = 1
         isk(nkstot/2+1:nkstot) = 2
      ELSE
         IF (band_struct_obj%nbnd_ispresent) THEN 
            nbnd  = band_struct_obj%nbnd
         ELSE 
            CALL errore ('qexsd_copy_band_structure: ','nbnd missing', 1)
         END IF  
         nbnd_up = nbnd
         nbnd_dw = nbnd
         isk(1:nkstot)   = 1 
      END IF
      !
      CALL qexsd_copy_efermi ( band_struct_obj, &
           nelec, ef, two_fermi_energies, ef_up, ef_dw )
      !
      IF ( .NOT. ALLOCATED(et) ) ALLOCATE( et(nbnd,nkstot) )
      !$acc enter data create(et)
      IF ( .NOT. ALLOCATED(wg) ) ALLOCATE( wg(nbnd,nkstot) )
      !
      DO ik =1, band_struct_obj%ndim_ks_energies
         IF ( band_struct_obj%lsda) THEN
            xk(:,ik) = band_struct_obj%ks_energies(ik)%k_point%k_point(:) 
            xk(:,ik + band_struct_obj%ndim_ks_energies) = xk(:,ik)
            wk(ik) = band_struct_obj%ks_energies(ik)%k_point%weight
            wk(ik + band_struct_obj%ndim_ks_energies ) = wk(ik) 
            et(1:nbnd_up,ik) = band_struct_obj%ks_energies(ik)%eigenvalues%vector(1:nbnd_up)
            et(1:nbnd_dw,ik+band_struct_obj%ndim_ks_energies) =  &
                 band_struct_obj%ks_energies(ik)%eigenvalues%vector(nbnd_up+1:nbnd_up+nbnd_dw)
            wg(1:nbnd_up,ik) = &
                 band_struct_obj%ks_energies(ik)%occupations%vector(1:nbnd_up)*wk(ik)
            wg(1:nbnd_dw,ik+band_struct_obj%ndim_ks_energies) =  &
                 band_struct_obj%ks_energies(ik)%occupations%vector(nbnd_up+1:nbnd_up+nbnd_dw)*wk(ik)
         ELSE 
            xk(:,ik) = band_struct_obj%ks_energies(ik)%k_point%k_point(:) 
            wk(ik) = band_struct_obj%ks_energies(ik)%k_point%weight
            et (1:nbnd,ik) = band_struct_obj%ks_energies(ik)%eigenvalues%vector(1:nbnd)
            wg (1:nbnd,ik) = band_struct_obj%ks_energies(ik)%occupations%vector(1:nbnd)*wk(ik)
         END IF
         !
      END DO
      !$acc update device(et)
      !
    END SUBROUTINE qexsd_copy_band_structure
    !
    SUBROUTINE qexsd_copy_efermi ( band_struct_obj, &
         nelec, ef, two_fermi_energies, ef_up, ef_dw, nbnd )
      !------------------------------------------------------------------------
      !
      USE qes_types_module, ONLY : band_structure_type
      !
      IMPLICIT NONE
      TYPE ( band_structure_type) :: band_struct_obj
      LOGICAL, INTENT(out) :: two_fermi_energies
      REAL(dp), INTENT(out):: nelec, ef, ef_up, ef_dw
      INTEGER, OPTIONAL, INTENT(out) :: nbnd
      !
      nelec = band_struct_obj%nelec
      two_fermi_energies = band_struct_obj%two_fermi_energies_ispresent 
      IF ( band_struct_obj%fermi_energy_ispresent) THEN 
         ef = band_struct_obj%fermi_energy
         ef_up = 0.d0
         ef_dw = 0.d0
      ELSE IF ( two_fermi_energies ) THEN 
         ef = 0.d0 
         ef_up = band_struct_obj%two_fermi_energies(1)
         ef_dw = band_struct_obj%two_fermi_energies(2)
      ELSE 
         ef = 0.d0
         ef_up = 0.d0
         ef_dw = 0.d0
      END IF      
      !
      IF ( PRESENT(nbnd) ) THEN
         !
         IF ( band_struct_obj%lsda ) THEN
            !
            IF (band_struct_obj%nbnd_ispresent) THEN
               nbnd  = band_struct_obj%nbnd / 2
            ELSE IF ( band_struct_obj%nbnd_up_ispresent .AND. band_struct_obj%nbnd_dw_ispresent ) THEN
               nbnd = (band_struct_obj%nbnd_up + band_struct_obj%nbnd_dw)/2
            ELSE
               CALL errore ('qexsd_copy_efermi: ','both nbnd and nbnd_up+nbnd_dw missing', 1)
            END IF
            !
         ELSE
            !
            IF (band_struct_obj%nbnd_ispresent) THEN
               nbnd  = band_struct_obj%nbnd
            ELSE
               CALL errore ('qexsd_copy_efermi: ','nbnd missing', 1)
            END IF
            !
         END IF
         !
      END IF
      !
    END SUBROUTINE qexsd_copy_efermi
    !-----------------------------------------------------------------------
    SUBROUTINE qexsd_copy_algorithmic_info ( algo_obj, &
         real_space, tqr, okvan, okpaw )
      USE qes_types_module, ONLY: algorithmic_info_type
      IMPLICIT NONE 
      TYPE(algorithmic_info_type),INTENT(IN)   :: algo_obj
      LOGICAL, INTENT(OUT) :: real_space, tqr, okvan, okpaw
      !
      tqr = algo_obj%real_space_q 
      real_space = algo_obj%real_space_beta
      okvan = algo_obj%uspp
      okpaw = algo_obj%paw
      !
    END SUBROUTINE qexsd_copy_algorithmic_info
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    SUBROUTINE qexsd_copy_twochem ( two_chem_obj, &
         twochem, nbnd_cond, nelec_cond, degauss_cond, ef_cond)
      USE qes_types_module, ONLY: two_chem_type
      IMPLICIT NONE 
      TYPE(two_chem_type),INTENT(IN)     ::  two_chem_obj
      LOGICAL,INTENT(OUT)               ::  twochem
      REAL(DP), INTENT(OUT)             ::  degauss_cond
      REAL(DP), INTENT(OUT)             ::  nelec_cond
      INTEGER, INTENT(OUT)              ::  nbnd_cond
      REAL(DP),OPTIONAL, INTENT(OUT)    :: ef_cond
      !
      twochem = two_chem_obj%twochem
      degauss_cond = two_chem_obj%degauss_cond
      nelec_cond = two_chem_obj%nelec_cond
      nbnd_cond = two_chem_obj%nbnd_cond
      IF (PRESENT(ef_cond)) ef_cond = two_chem_obj%ef_cond 
      !
    END SUBROUTINE qexsd_copy_twochem
    !-----------------------------------------------------------------------

    !
    !---------------------------------------------------------------------------
    SUBROUTINE qexsd_copy_efield ( efield_obj, tefield, dipfield, edir, &
         emaxpos, eopreg, eamp, gate, zgate, &
         block_, block_1, block_2, block_height, relaxz )
      !---------------------------------------------------------------------------
      USE qes_types_module,    ONLY: electric_field_type
      IMPLICIT NONE 
      ! 
      TYPE ( electric_field_type),OPTIONAL, INTENT(IN)    :: efield_obj
      LOGICAL, INTENT(OUT) :: tefield, dipfield
      INTEGER, INTENT(INOUT) :: edir
      LOGICAL, INTENT(INOUT) :: gate, block_, relaxz 
      REAL(dp), INTENT(INOUT) :: emaxpos, eopreg, eamp, &
              zgate, block_1, block_2, block_height
      !
      !
      tefield = .FALSE. 
      dipfield = .FALSE. 
      IF ( .NOT. PRESENT( efield_obj) ) RETURN 
      IF (TRIM(efield_obj%electric_potential) == 'sawtooth_potential') THEN 
         tefield = .TRUE. 
         IF ( efield_obj%dipole_correction_ispresent ) THEN 
            dipfield = efield_obj%dipole_correction
         ELSE 
            dipfield = .FALSE. 
         END IF
         IF ( efield_obj%electric_field_direction_ispresent ) THEN 
            edir = efield_obj%electric_field_direction
         ELSE 
            edir = 3 
         END IF
         IF ( efield_obj%potential_max_position_ispresent ) THEN 
            emaxpos = efield_obj%potential_max_position
         ELSE 
            emaxpos = 5d-1
         END IF
         IF ( efield_obj%potential_decrease_width_ispresent ) THEN 
            eopreg = efield_obj%potential_decrease_width
         ELSE 
            eopreg = 1.d-1
         END IF
         IF ( efield_obj%electric_field_amplitude_ispresent ) THEN 
            eamp = efield_obj%electric_field_amplitude
         ELSE 
            eamp = 1.d-3
         END IF
         IF (efield_obj%gate_settings_ispresent) THEN 
            gate = efield_obj%gate_settings%use_gate
            IF (efield_obj%gate_settings%zgate_ispresent) &
                 zgate     = efield_obj%gate_settings%zgate
            IF (efield_obj%gate_settings%relaxz_ispresent) &
                 relaxz   = efield_obj%gate_settings%relaxz
            IF (efield_obj%gate_settings%block_ispresent) &
                 block_    = efield_obj%gate_settings%block
            IF (efield_obj%gate_settings%block_1_ispresent) &
                 block_1 = efield_obj%gate_settings%block_1
            IF (efield_obj%gate_settings%block_2_ispresent) &
                 block_2 = efield_obj%gate_settings%block_2
            IF (efield_obj%gate_settings%block_height_ispresent) &
                 block_height = efield_obj%gate_settings%block_height
         END IF
      END IF
      !
    END SUBROUTINE qexsd_copy_efield
    !
    !--------------------------------------------------------------------------
    SUBROUTINE qexsd_copy_magnetization ( magnetization_obj, &
         lsda, noncolin, lspinorb, domag, tot_magnetization )
      !------------------------------------------------------------------------
      ! 
      USE qes_types_module, ONLY : magnetization_type
      !
      IMPLICIT NONE 
      !
      TYPE ( magnetization_type ) ,INTENT(IN)    :: magnetization_obj
      LOGICAL, INTENT(OUT)  :: lsda, noncolin, lspinorb, domag
      REAL(dp), INTENT(OUT) :: tot_magnetization
      ! 
      lsda  =   magnetization_obj%lsda
      noncolin = magnetization_obj%noncolin  
      lspinorb = magnetization_obj%spinorbit 
      IF (magnetization_obj%do_magnetization_ispresent) THEN 
        domag =   magnetization_obj%do_magnetization
      ELSE 
        domag = .FALSE.
      END IF
      IF (magnetization_obj%total_ispresent) THEN 
        tot_magnetization = magnetization_obj%total
      ELSE IF (magnetization_obj%total_vec_ispresent) THEN 
        tot_magnetization = SQRT(dot_product(magnetization_obj%total_vec, magnetization_obj%total_vec))
      ELSE 
        tot_magnetization = 0._DP 
      END IF 
      !
    END SUBROUTINE qexsd_copy_magnetization
    !-----------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    SUBROUTINE qexsd_copy_kpoints ( band_struct_obj, nks_start, xk_start,&
         wk_start, nk1, nk2, nk3, k1, k2, k3, occupations, smearing, degauss )
    !---------------------------------------------------------------------------
       !
       USE qes_types_module, ONLY : band_structure_type
       !
       IMPLICIT NONE
       !
       TYPE ( band_structure_type ),INTENT(IN)    :: band_struct_obj
       INTEGER,  INTENT(out) :: nks_start, nk1, nk2, nk3, k1, k2, k3 
       REAL(dp), ALLOCATABLE, INTENT(inout) :: xk_start(:,:), wk_start(:)
       REAL(dp), INTENT(out) :: degauss
       CHARACTER(LEN=*), intent(out) :: smearing, occupations
       !
       INTEGER :: ik
       !
       occupations = TRIM ( band_struct_obj%occupations_kind%occupations ) 
       smearing    = TRIM ( band_struct_obj%smearing%smearing ) 
       IF (band_struct_obj%smearing%degauss_ispresent) THEN 
         degauss     = band_struct_obj%smearing%degauss
       ELSE 
         degauss = 0._DP
       END IF 
       !   
       IF ( band_struct_obj%starting_k_points%monkhorst_pack_ispresent ) THEN 
          nks_start = 0 
          nk1 = band_struct_obj%starting_k_points%monkhorst_pack%nk1 
          nk2 = band_struct_obj%starting_k_points%monkhorst_pack%nk2
          nk3 = band_struct_obj%starting_k_points%monkhorst_pack%nk3 
           k1 = band_struct_obj%starting_k_points%monkhorst_pack%k1
           k2 = band_struct_obj%starting_k_points%monkhorst_pack%k2
           k3 = band_struct_obj%starting_k_points%monkhorst_pack%k3
       ELSE IF (band_struct_obj%starting_k_points%nk_ispresent ) THEN 
           nks_start = band_struct_obj%starting_k_points%nk
           IF ( nks_start > 0 ) THEN 
              IF ( .NOT. ALLOCATED(xk_start) ) ALLOCATE (xk_start(3,nks_start))
              IF ( .NOT. ALLOCATED(wk_start) ) ALLOCATE (wk_start(nks_start))
              IF ( nks_start == size( band_struct_obj%starting_k_points%k_point ) ) THEN 
                 DO ik =1, nks_start
                    xk_start(:,ik) = band_struct_obj%starting_k_points%k_point(ik)%k_point(:) 
                    IF ( band_struct_obj%starting_k_points%k_point(ik)%weight_ispresent) THEN 
                        wk_start(ik) = band_struct_obj%starting_k_points%k_point(ik)%weight 
                    ELSE 
                        wk_start(ik) = 0.d0
                    END IF 
                 END DO
              ELSE
                 CALL infomsg ( "qexsd_copy_kp: ", &
                      "actual number of start kpoint not equal to nks_start, set nks_start=0")  
                 nks_start = 0 
              END IF
           END IF
       ELSE 
          CALL errore ("qexsd_copy_kp: ", &
               " no information found for initializing brillouin zone information", 1)
       END IF  
       ! 
     END SUBROUTINE qexsd_copy_kpoints
     !
  !
  !---------------------------------------------------------------
  SUBROUTINE qexsd_copy_esm( pbc_obj, bc, nfit, w, efield, a) 
    !------------------------------------------------------------
    USE qes_types_module, ONLY: outputPBC_type 
    IMPLICIT NONE 
    TYPE(outputPBC_type),INTENT(IN) :: pbc_obj 
    CHARACTER(LEN=3),INTENT(OUT)    :: bc 
    INTEGER,INTENT(OUT)             :: nfit 
    REAL(DP),INTENT(OUT)             :: w 
    REAL(DP),INTENT(OUT)             :: efield
    REAL(DP),INTENT(OUT)             :: a 
    ! 
    IF (pbc_obj%esm_ispresent) THEN 
      bc = TRIM(pbc_obj%esm%bc) 
      nfit = pbc_obj%esm%nfit 
      w = pbc_obj%esm%w 
      efield = pbc_obj%esm%efield
      a = pbc_obj%esm%a 
    ELSE 
      CALL errore("qexsd_copy_esm","esm object not present in input", 1) 
    END IF 
  END SUBROUTINE qexsd_copy_esm 
  !
  !---------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_rism3d( rism3d_obj, pseudo_dir, nsolV, solVs, molfile, ecutsolv )
    !---------------------------------------------------------------------------
    !
    USE qes_types_module, ONLY : rism3d_type
    USE molecule_types,   ONLY : molecule, nullify_molecule
    !
    IMPLICIT NONE
    !
    TYPE(rism3d_type), INTENT(IN)  :: rism3d_obj
    CHARACTER(LEN=*),  INTENT(IN)  :: pseudo_dir
    INTEGER,           INTENT(OUT) :: nsolV
    TYPE(molecule),    INTENT(INOUT), ALLOCATABLE :: solVs(:)
    CHARACTER(LEN=*),  INTENT(OUT) :: molfile(:)
    REAL(DP),          INTENT(OUT) :: ecutsolv
    !
    INTEGER :: isolV
    !
    IF ( rism3d_obj%molec_dir_ispresent ) THEN
       IF ( TRIM(pseudo_dir) /= TRIM(rism3d_obj%molec_dir) ) THEN
          CALL errore ("qexsd_copy_rism3d:", "pseudo_dir /= molec_dir", 1)
       END IF
    END IF
    !
    nsolV = rism3d_obj%nmol
    !
    IF ( .NOT. ALLOCATED(solVs) ) ALLOCATE(solVs(nsolV))
    !
    DO isolV = 1, nsolV
       !
       CALL nullify_molecule(solVs(isolV))
       solVs(isolV)%name       = TRIM(rism3d_obj%solvent(isolV)%label)
       solVs(isolV)%density    = rism3d_obj%solvent(isolV)%density1
       solVs(isolV)%subdensity = rism3d_obj%solvent(isolV)%density2
       !
       molfile(isolV) = TRIM(rism3d_obj%solvent(isolV)%molec_file)
       !
    END DO
    !
    ecutsolv = rism3d_obj%ecutsolv
    !
  END SUBROUTINE qexsd_copy_rism3d
  !
  !---------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_rismlaue( rismlaue_obj, both_hands, laue_nfit, ireference, qsol, &
                                  starting_r, expand_r, buffer_r, buffer_ru, buffer_rv,   &
                                  starting_l, expand_l, buffer_l, buffer_lu, buffer_lv )
    !---------------------------------------------------------------------------
    !
    USE qes_types_module, ONLY : rismlaue_type
    !
    IMPLICIT NONE
    !
    TYPE(rismlaue_type), INTENT(IN)  :: rismlaue_obj
    LOGICAL,             INTENT(OUT) :: both_hands
    INTEGER,             INTENT(OUT) :: laue_nfit
    INTEGER,             INTENT(OUT) :: ireference
    REAL(DP),            INTENT(OUT) :: qsol
    REAL(DP),            INTENT(OUT) :: starting_r, starting_l
    REAL(DP),            INTENT(OUT) :: expand_r,   expand_l
    REAL(DP),            INTENT(OUT) :: buffer_r,   buffer_l
    REAL(DP),            INTENT(OUT) :: buffer_ru,  buffer_lu
    REAL(DP),            INTENT(OUT) :: buffer_rv,  buffer_lv
    !
    both_hands = rismlaue_obj%both_hands
    laue_nfit  = rismlaue_obj%nfit
    ireference = rismlaue_obj%pot_ref
    qsol       = rismlaue_obj%charge
    !
    starting_r = rismlaue_obj%right_start
    expand_r   = rismlaue_obj%right_expand
    buffer_r   = rismlaue_obj%right_buffer
    buffer_ru  = rismlaue_obj%right_buffer_u
    buffer_rv  = rismlaue_obj%right_buffer_v
    !
    starting_l = rismlaue_obj%left_start
    expand_l   = rismlaue_obj%left_expand
    buffer_l   = rismlaue_obj%left_buffer
    buffer_lu  = rismlaue_obj%left_buffer_u
    buffer_lv  = rismlaue_obj%left_buffer_v
    !
  END SUBROUTINE qexsd_copy_rismlaue
  !
END MODULE qexsd_copy
