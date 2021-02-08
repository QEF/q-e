! Copyright (C) 2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE qexsd_copy
  !----------------------------------------------------------------------------
  !
  ! This module contains some common subroutines used to copy data read from
  ! XML format into data used by the Quantum ESPRESSO package.
  !
  ! Written by Paolo Giannozzi, building upon pre-existing code
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
       qexsd_copy_efield, qexsd_copy_magnetization, qexsd_copy_kpoints, &
       qexsd_copy_efermi
  !
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
       nat, tau, ityp, alat, a1, a2, a3, ibrav )
  !--------------------------------------------------------------------------
    
    USE qes_types_module, ONLY : atomic_structure_type
    !
    IMPLICIT NONE 
    !
    TYPE ( atomic_structure_type ),INTENT(IN)  :: atomic_structure
    INTEGER, INTENT(in) :: nsp 
    CHARACTER(LEN = 3), INTENT(in) :: atm(:)
    !
    INTEGER, INTENT(out)  :: nat, ibrav
    REAL(dp), INTENT(out) :: alat, a1(:), a2(:), a3(:)
    INTEGER, INTENT(inout),  ALLOCATABLE :: ityp(:)
    REAL(dp), INTENT(inout), ALLOCATABLE :: tau(:,:)
    !
    CHARACTER(LEN=3), ALLOCATABLE :: symbols(:)
    INTEGER :: iat, idx, isp
    !
    nat = atomic_structure%nat 
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
       noinv, nosym, no_t_rev, flags_obj )
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
       lda_plus_U, lda_plus_U_kind, U_projection, Hubbard_l, Hubbard_lmax, &
       Hubbard_l_back, Hubbard_l1_back, backall, Hubbard_lmax_back, Hubbard_alpha_back, &
       Hubbard_U, Hubbard_U_back, Hubbard_J0, Hubbard_alpha, Hubbard_beta, Hubbard_J, &
       vdw_corr, scal6, lon_rcut, vdw_isolated )
    !-------------------------------------------------------------------
    ! 
    USE qes_types_module, ONLY : dft_type
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
    LOGICAL, INTENT(out) :: lda_plus_U
    INTEGER, INTENT(inout) :: lda_plus_U_kind, Hubbard_lmax, Hubbard_lmax_back
    CHARACTER(LEN=*), INTENT(inout) :: U_projection
    INTEGER, INTENT(inout) :: Hubbard_l(:), Hubbard_l_back(:), Hubbard_l1_back(:) 
    REAL(dp), INTENT(inout) :: Hubbard_U(:), Hubbard_U_back(:), Hubbard_J0(:), Hubbard_J(:,:), &
                               Hubbard_alpha(:), Hubbard_alpha_back(:), Hubbard_beta(:)
    LOGICAL, INTENT(inout) :: backall(:)
    OPTIONAL    :: Hubbard_U_back, Hubbard_l_back, Hubbard_lmax_back, Hubbard_alpha_back, &
                   Hubbard_l1_back 
    !
    CHARACTER(LEN=*), INTENT(out) :: vdw_corr
    REAL(dp), INTENT(inout) :: scal6, lon_rcut
    LOGICAL, INTENT(inout) :: vdw_isolated
    !
    CHARACTER(LEN=256 ) :: label
    CHARACTER(LEN=3 )   :: symbol
    INTEGER :: ihub, isp
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
       Hubbard_U = 0.0_DP
       Hubbard_U_back =0.0_DP
       Hubbard_alpha = 0.0_DP
       Hubbard_alpha_back = 0.0_DP
       Hubbard_J = 0.0_DP
       Hubbard_J0 = 0.0_DP
       Hubbard_beta = 0.0_DP
       lda_plus_u_kind = dft_obj%dftU%lda_plus_u_kind
       U_projection = TRIM ( dft_obj%dftU%U_projection_type )
       Hubbard_l =-1 
       Hubbard_l_back =-1 
       backall = .false.
       !
       IF ( dft_obj%dftU%Hubbard_U_ispresent) THEN 
          loop_on_hubbardU:DO ihub =1, dft_obj%dftU%ndim_Hubbard_U
             symbol = TRIM(dft_obj%dftU%Hubbard_U(ihub)%specie)
             label  = TRIM(dft_obj%dftU%Hubbard_U(ihub)%label ) 
             loop_on_speciesU:DO isp = 1, nsp
                IF ( TRIM(symbol) == TRIM ( atm(isp) ) ) THEN 
                     Hubbard_U(isp) = dft_obj%dftU%Hubbard_U(ihub)%HubbardCommon
                     SELECT CASE ( TRIM (label))
                     CASE ( '1s', '2s', '3s', '4s', '5s', '6s', '7s' ) 
                        Hubbard_l(isp) = 0 
                     CASE ( '2p', '3p', '4p', '5p', '6p' ) 
                        Hubbard_l(isp) = 1 
                     CASE ( '3d', '4d', '5d' ) 
                        Hubbard_l( isp ) = 2 
                     CASE ( '4f', '5f' )  
                        Hubbard_l(isp ) = 3
                     CASE  default 
                        IF (Hubbard_U(isp)/=0) &
                             CALL errore ("qexsd_copy_dft:", "unrecognized label for Hubbard "//label, 1 ) 
                     END SELECT
                     EXIT loop_on_speciesU
                END IF 
             END DO loop_on_speciesU
          END DO loop_on_hubbardU
       END IF 
       ! 
       IF ( dft_obj%dftU%Hubbard_U_back_ispresent) THEN
          loop_on_hubbardUback:DO ihub =1, dft_obj%dftU%ndim_Hubbard_U_back
             symbol = TRIM(dft_obj%dftU%Hubbard_U_back(ihub)%specie)
             label  = TRIM(dft_obj%dftU%Hubbard_U_back(ihub)%label )
             loop_on_speciesU_back:DO isp = 1, nsp
                IF ( TRIM(symbol) == TRIM ( atm(isp) ) ) THEN
                     Hubbard_U_back(isp) = dft_obj%dftU%Hubbard_U_back(ihub)%HubbardCommon
                     EXIT loop_on_speciesU_back
                END IF
             END DO loop_on_speciesU_back
          END DO loop_on_hubbardUback
          IF (.NOT. dft_obj%dftU%Hubbard_back_ispresent) CALL errore("qexsd_copy:", &
                                                        "internal error: U_back is present but not Hub_back",1) 
          loop_hubbardBack: DO ihub =1, dft_obj%dftU%ndim_Hubbard_back
              symbol = TRIM(dft_obj%dftU%Hubbard_back(ihub)%species) 
              loop_on_species_2:DO isp = 1, nsp
                 IF ( TRIM(symbol) == TRIM(atm(isp))) THEN 
                    Hubbard_l_back(isp) = dft_obj%dftU%Hubbard_back(ihub)%l_number(1)%backL
                    SELECT CASE ( TRIM (dft_obj%dftU%Hubbard_back(ihub)%background)) 
                       CASE ('one_orbital') 
                         backall(isp) = .FALSE. 
                       CASE ('two_orbitals') 
                         backall(isp)  = .TRUE. 
                         Hubbard_l1_back(isp) = dft_obj%dftU%Hubbard_back(ihub)%l_number(2)%backL 
                    END SELECT
                    EXIT loop_on_species_2 
                 END IF
              END DO loop_on_species_2
         END DO loop_hubbardBack

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
       Hubbard_lmax      = MAXVAL( Hubbard_l(1:nsp) )
       Hubbard_lmax_back = MAXVAL( Hubbard_l_back(1:nsp) ) 
       !  
    END IF

      IF ( dft_obj%vdW_ispresent ) THEN 
         vdw_corr = TRIM( dft_obj%vdW%vdw_corr ) 
      ELSE
         vdw_corr = ''
      END IF
      
      IF ( dft_obj%vdW_ispresent ) THEN 
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
   
    END SUBROUTINE qexsd_copy_dft
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_copy_band_structure( band_struct_obj, lsda, nkstot, &
         isk, natomwfc, nbnd, nbnd_up, nbnd_dw, nelec, xk, wk, wg, &
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
      INTEGER, INTENT(out) :: nkstot, natomwfc, nbnd, nbnd_up, nbnd_dw, &
              isk(:)
      REAL(dp), INTENT(out):: nelec, ef, ef_up, ef_dw, xk(:,:), wk(:)
      REAL(dp), INTENT(inout), ALLOCATABLE ::  wg(:,:), et(:,:)
      !
      LOGICAL :: two_fermi_energies
      INTEGER :: ik
      ! 
      lsda = band_struct_obj%lsda
      nkstot = band_struct_obj%nks 
      natomwfc = band_struct_obj%num_of_atomic_wfc
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
      !
    END SUBROUTINE qexsd_copy_band_structure
    !
    SUBROUTINE qexsd_copy_efermi ( band_struct_obj, &
         nelec, ef, two_fermi_energies, ef_up, ef_dw )
      !------------------------------------------------------------------------
      !
      USE qes_types_module, ONLY : band_structure_type
      !
      IMPLICIT NONE
      TYPE ( band_structure_type) :: band_struct_obj
      LOGICAL, INTENT(out) :: two_fermi_energies
      REAL(dp), INTENT(out):: nelec, ef, ef_up, ef_dw
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
      domag =   magnetization_obj%do_magnetization 
      tot_magnetization = magnetization_obj%total
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
       degauss     = band_struct_obj%smearing%degauss
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
   END MODULE qexsd_copy
