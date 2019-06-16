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
  PUBLIC:: qexsd_copy_geninfo, qexsd_copy_parallel_info, qexsd_copy_dim, &
       qexsd_copy_atomic_species, qexsd_copy_atomic_structure, &
       qexsd_copy_symmetry, &
       qexsd_copy_basis_set, qexsd_copy_dft, qexsd_copy_band_structure
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
  SUBROUTINE qexsd_copy_dim (atomic_structure, band_structure, &
         nat, nkstot, nbnd ) 
      !
    USE qes_types_module, ONLY : atomic_structure_type, band_structure_type
    IMPLICIT NONE 
    !
    TYPE ( atomic_structure_type ),INTENT(IN)  :: atomic_structure
    TYPE ( band_structure_type ),INTENT(IN)    :: band_structure 
    INTEGER, INTENT(OUT) :: nat, nkstot, nbnd
    !
    LOGICAL :: lsda
    !
    nat = atomic_structure%nat 
    nkstot =   band_structure%nks  
    IF (band_structure%nbnd_ispresent) THEN
       nbnd = band_structure%nbnd
    ELSE IF ( band_structure%nbnd_up_ispresent .AND. band_structure%nbnd_dw_ispresent) THEN
       nbnd = ( band_structure%nbnd_up + band_structure%nbnd_dw )
    ELSE 
       CALL errore('init_vars_from_schema: check xml file !!', &
                   'nbnd or nbnd_up+nbnd_dw are missing in band_structure element', 1)
    END IF     
    lsda  =    band_structure%lsda
    IF ( lsda ) THEN
       nkstot = nkstot * 2 
       nbnd   = nbnd / 2
    END IF

  END SUBROUTINE qexsd_copy_dim
  !
  !--------------------------------------------------------------------------
  SUBROUTINE qexsd_copy_atomic_species (atomic_species, nsp, atm, amass, &
       psfile, pseudo_dir)
    !---------------------------------------------------------------------------    !
    USE qes_types_module, ONLY : atomic_species_type
    !
    IMPLICIT NONE 
    !
    TYPE ( atomic_species_type ),INTENT(IN)    :: atomic_species
    INTEGER, INTENT(out) :: nsp
    CHARACTER(LEN=*), INTENT(out) :: atm(:)
    CHARACTER(LEN=*), OPTIONAL, INTENT(out) :: psfile(:), pseudo_dir
    REAL(dp), INTENT(out) :: amass(:)
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
    USE constants,        ONLY : pi
    !
    IMPLICIT NONE 
    !
    TYPE ( atomic_structure_type ),INTENT(IN)  :: atomic_structure
    INTEGER, INTENT(in) :: nsp 
    CHARACTER(LEN = 3), INTENT(in) :: atm(:)
    !
    INTEGER, INTENT(out)  :: nat, ibrav, ityp(:)
    REAL(dp), INTENT(out) :: alat, a1(:), a2(:), a3(:), tau(:,:)
    !
    CHARACTER(LEN=3), ALLOCATABLE :: symbols(:)
    INTEGER :: iat, idx, isp
    !
    nat = atomic_structure%nat 
    alat = atomic_structure%alat 
    IF ( atomic_structure%bravais_index_ispresent ) THEN 
       ibrav = atomic_structure%bravais_index 
    ELSE 
       ibrav = 0
    END IF
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
  SUBROUTINE qexsd_copy_symmetry ( symms_obj, &
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
       exxdiv_treatment, x_gamma_extrapolation, ecutvcut, &
       lda_plus_U, lda_plus_U_kind, U_projection, Hubbard_l, Hubbard_lmax, &
       Hubbard_U, Hubbard_J0, Hubbard_alpha, Hubbard_beta, Hubbard_J, &
       vdw_corr,  llondon, ts_vdw, lxdm, inlc, vdw_table_name, scal6, &
       lon_rcut, vdw_isolated)
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
         ecutvcut
    INTEGER, INTENT(inout) :: nq1, nq2, nq3
    LOGICAL, INTENT(inout) :: x_gamma_extrapolation
    !
    LOGICAL, INTENT(out) :: lda_plus_U
    INTEGER, INTENT(inout) :: lda_plus_U_kind, Hubbard_lmax
    CHARACTER(LEN=*), INTENT(inout) :: U_projection
    INTEGER, INTENT(inout) :: Hubbard_l(:)
    REAL(dp), INTENT(inout) :: Hubbard_U(:), Hubbard_J0(:), Hubbard_J(:,:), &
         Hubbard_alpha(:), Hubbard_beta(:)
    !
    CHARACTER(LEN=256), INTENT(out) :: vdw_corr
    CHARACTER(LEN=256), INTENT(inout) :: vdw_table_name
    LOGICAL, INTENT(out) :: llondon, ts_vdw, lxdm
    INTEGER, INTENT(inout):: inlc
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
    END IF
    !
    lda_plus_u = dft_obj%dftU_ispresent 
    IF ( lda_plus_u ) THEN 
       lda_plus_u_kind = dft_obj%dftU%lda_plus_u_kind
       U_projection = TRIM ( dft_obj%dftU%U_projection_type )
       Hubbard_l =-1 
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
         Hubbard_lmax = MAXVAL( Hubbard_l(1:nsp) )
      END IF

      IF ( dft_obj%vdW_ispresent ) THEN 
         vdw_corr = TRIM( dft_obj%vdW%vdw_corr ) 
      ELSE
         vdw_corr = ''
      END IF
      SELECT CASE( TRIM( dft_obj%vdW%vdw_corr ) )
         !
      CASE( 'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d' )
         !
         llondon= .TRUE.
         ts_vdw= .FALSE.
         lxdm   = .FALSE.
         !
      CASE( 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler' )
         !
         llondon= .FALSE.
         ts_vdw= .TRUE.
         lxdm   = .FALSE.
         !
      CASE( 'XDM', 'xdm' )
         !
         llondon= .FALSE.
         ts_vdw= .FALSE.
         lxdm   = .TRUE.
         !
      CASE DEFAULT
         !
         llondon= .FALSE.
         ts_vdw = .FALSE.
         lxdm   = .FALSE.
         !
      END SELECT
      IF ( dft_obj%vdW_ispresent ) THEN 
         SELECT CASE ( TRIM (dft_obj%vdW%non_local_term))
         CASE ('vdw1')  
            inlc = 1
         CASE ('vdw2') 
            inlc = 2
         CASE ('vv10' ) 
            inlc = 3 
         CASE ( 'vdW-DF-x') 
            inlc = 4
         CASE ( 'vdW-DF-y')
            inlc = 5
         CASE ( 'vdW-DF-z')
            inlc = 6
         CASE default 
            inlc = 0 
         END SELECT
         IF (inlc == 0 ) THEN 
            vdw_table_name = ' '
         ELSE IF ( inlc == 3 ) THEN 
            vdw_table_name = 'rVV10_kernel_table'
         ELSE
            vdw_table_name = 'vdW_kernel_table'
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
   
    END SUBROUTINE qexsd_copy_dft
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexsd_copy_band_structure( band_struct_obj, lsda, nkstot, &
         isk, natomwfc, nbnd_up, nbnd_dw, nelec, wk, wg, ef, ef_up, ef_dw, et )
      !------------------------------------------------------------------------
      !
      USE qes_types_module, ONLY : band_structure_type
      !
      IMPLICIT NONE
      TYPE ( band_structure_type)         :: band_struct_obj
      LOGICAL, INTENT(out) :: lsda
      INTEGER, INTENT(out) :: nkstot, natomwfc, nbnd_up, nbnd_dw, isk(:)
      REAL(dp), INTENT(out):: nelec, wk(:), wg(:,:)
      REAL(dp), INTENT(out):: ef, ef_up, ef_dw, et(:,:)
      !
      INTEGER :: ik, nbnd
      ! 
      lsda = band_struct_obj%lsda
      nkstot = band_struct_obj%nks 
      IF ( lsda) THEN 
         nkstot = nkstot * 2 
         isk(1:nkstot/2) = 1
         isk(nkstot/2+1:nkstot) = 2 
      ELSE 
         isk(1:nkstot)   = 1 
      END IF
      ! 
      nelec = band_struct_obj%nelec
      nbnd  = band_struct_obj%nbnd 
      natomwfc = band_struct_obj%num_of_atomic_wfc
      IF ( band_struct_obj%fermi_energy_ispresent) THEN 
         ef = band_struct_obj%fermi_energy
         ef_up = 0.d0
         ef_dw = 0.d0
      ELSE IF ( band_struct_obj%two_fermi_energies_ispresent ) THEN 
         ef = 0.d0 
         ef_up = band_struct_obj%two_fermi_energies(1)
         ef_dw = band_struct_obj%two_fermi_energies(2)
      ELSE 
         ef = 0.d0
         ef_up = 0.d0
         ef_dw = 0.d0
      END IF
      DO ik =1, band_struct_obj%ndim_ks_energies
         IF ( band_struct_obj%lsda) THEN
            IF ( band_struct_obj%nbnd_up_ispresent .AND. band_struct_obj%nbnd_dw_ispresent) THEN
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
            wk(ik) = band_struct_obj%ks_energies(ik)%k_point%weight
            wk( ik + band_struct_obj%ndim_ks_energies ) = wk(ik) 
            et(1:nbnd_up,ik) = band_struct_obj%ks_energies(ik)%eigenvalues%vector(1:nbnd_up)
            et(1:nbnd_dw,ik+band_struct_obj%ndim_ks_energies) =  &
                 band_struct_obj%ks_energies(ik)%eigenvalues%vector(nbnd_up+1:nbnd_up+nbnd_dw)
            wg(1:nbnd_up,ik) = band_struct_obj%ks_energies(ik)%occupations%vector(1:nbnd_up)*wk(ik)
            wg(1:nbnd_dw,ik+band_struct_obj%ndim_ks_energies) =  &
                 band_struct_obj%ks_energies(ik)%occupations%vector(nbnd_up+1:nbnd_up+nbnd_dw)*wk(ik)
         ELSE 
            wk(ik) = band_struct_obj%ks_energies(ik)%k_point%weight
            nbnd = band_struct_obj%ks_energies(ik)%eigenvalues%size
            et (1:nbnd,ik) = band_struct_obj%ks_energies(ik)%eigenvalues%vector(1:nbnd)
            wg (1:nbnd,ik) = band_struct_obj%ks_energies(ik)%occupations%vector(1:nbnd)*wk(ik)
            nbnd_up = nbnd
            nbnd_dw = nbnd
         END IF
      END DO
    END SUBROUTINE qexsd_copy_band_structure
    !
  END MODULE qexsd_copy
