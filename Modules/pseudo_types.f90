!
! Copyright (C) 2002-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      MODULE pseudo_types

!  this module contains the definitions of several TYPE structures,
!  together with their allocation/deallocation routines

        USE kinds, ONLY: DP
        USE parameters, ONLY: cp_lmax, lmaxx
        use radial_grids, ONLY: ndmx, radial_grid_type

        IMPLICIT NONE
        SAVE
!
TYPE :: paw_t
   !
   ! Type describing a PAW dataset (temporary).
   ! Functions are defined on a logarithmic radial mesh.
   !
   CHARACTER(LEN=2) :: symbol
   REAL (DP) :: zval
   REAL (DP) :: z
   CHARACTER(LEN=80) :: dft
   TYPE(radial_grid_type) :: grid
   REAL (DP) :: rmatch_augfun  ! the matching radius for augmentation charges
   LOGICAL :: nlcc ! nonlinear core correction
   INTEGER :: nwfc ! number of wavefunctions/projectors
   INTEGER :: lmax ! maximum angular momentum of projectors
   INTEGER, POINTER :: l(:) !l(nwfsx) ! angular momentum of projectors
   INTEGER, POINTER :: ikk(:) !ikk(nwfsx) ! cutoff radius for the projectors
   INTEGER :: irc ! r(irc) = radius of the augmentation sphere
   CHARACTER(LEN=2),POINTER ::  els (:) ! the name of the wavefunction
   REAL (DP), POINTER :: &
        oc(:), &          !(nwfsx) the occupations
        enl(:), &         !(nwfsx) the energy of the wavefunctions
        jj (:), &         ! the total angular momentum
        rcutus (:), &     ! the cutoff
        aewfc(:,:), &     !(ndmx,nwfsx) all-electron wavefunctions
        pswfc(:,:), &     !(ndmx,nwfsx) pseudo wavefunctions
        proj(:,:), &      !(ndmx,nwfsx) projectors
        augfun(:,:,:,:), &!(ndmx,nwfsx,nwfsx,0:2*lmaxx+1),
        augmom(:,:,:), &  !(nwfsx,nwfsx,0:2*lmaxx) moments of the augmentation functions
        aeccharge(:), &   !(ndmx) AE core charge * 4PI r^2
        psccharge(:), &   !(ndmx) PS core charge * 4PI r^2
        pscharge(:), &    !(ndmx) PS charge * 4PI r^2
        aeloc(:), &       !(ndmx) descreened AE potential: v_AE-v_H[n1]-v_XC[n1+nc]
        psloc(:), &       !(ndmx) descreened local PS potential: v_PS-v_H[n~+n^]-v_XC[n~+n^+n~c]
        kdiff(:,:), &     !(nwfsx,nwfsx) kinetic energy differences
        dion(:,:)         !(nwfsx,nwfsx) descreened D coeffs
!!!  Notes about screening:
!!!       Without nlcc, the local PSpotential is descreened with n~+n^ only.
!!!       The local AEpotential is descreened ALWAYS with n1+nc. This improves
!!!       the accuracy, and will not cost in the plane wave code (atomic
!!!       contribution only).
END TYPE paw_t

!

        TYPE pseudo_upf
          CHARACTER(LEN=80):: generated   ! 
          CHARACTER(LEN=80):: date_author ! Misc info
          CHARACTER(LEN=80):: comment     !
          CHARACTER(LEN=2) :: psd       ! Element label
          CHARACTER(LEN=20) :: typ      ! Pseudo type ( NC or US )
          LOGICAL  :: tvanp             ! .true. if Ultrasoft
          LOGICAL :: nlcc               ! Non linear core corrections
          CHARACTER(LEN=20) :: dft      ! Exch-Corr type
          REAL(DP) :: zp                ! z valence
          REAL(DP) :: etotps            ! total energy
          REAL(DP) :: ecutwfc           ! suggested cut-off for wfc
          REAL(DP) :: ecutrho           ! suggested cut-off for rho

          LOGICAL :: has_so             ! if .true. includes spin-orbit
          REAL(DP) :: xmin              ! the minimum x of the linear mesh
          REAL(DP) :: rmax              ! the maximum radius of the mesh
          REAL(DP) :: zmesh             ! the nuclear charge used for mesh
          REAL(DP) :: dx                ! the deltax of the linear mesh
          INTEGER, POINTER :: nn(:)     ! nn(nwfc)
          REAL(DP), POINTER :: rcut(:)  ! cut-off radius(nbeta)
          REAL(DP), POINTER :: rcutus(:)! cut-off ultrasoft radius (nbeta)
          REAL(DP), POINTER :: epseu(:) ! energy (nwfc)
          REAL(DP), POINTER :: jchi(:)  ! jchi(nwfc)
          REAL(DP), POINTER :: jjj(:)   ! jjj(nbeta)

          INTEGER :: nv                 ! UPF file version number
          INTEGER :: lmax               ! maximum angular momentum component
          INTEGER :: mesh               ! number of point in the radial mesh
          INTEGER :: nwfc               ! number of wavefunctions
          INTEGER :: nbeta              ! number of projectors
          INTEGER :: kkbeta             ! kkbeta=max(kbeta(:))
          !  kbeta<=mesh is the number of grid points for each beta function
          !              beta(r,nb) = 0 for r > r(kbeta(nb))
          ! kkbeta<=mesh is the largest of such number so that for all beta
          !              beta(r,nb) = 0 for r > r(kkbeta)
          CHARACTER(LEN=2), POINTER :: els(:)  ! els(nwfc)
          CHARACTER(LEN=2), POINTER :: els_beta(:)  ! els(nbeta)
          INTEGER, POINTER :: lchi(:)     ! lchi(nwfc)
          REAL(DP), POINTER :: oc(:)      ! oc(nwfc)
          REAL(DP), POINTER :: r(:)       ! r(mesh)
          REAL(DP), POINTER :: rab(:)     ! rab(mesh)
          REAL(DP), POINTER :: rho_atc(:) ! rho_atc(mesh)
          REAL(DP), POINTER :: vloc(:)    ! vloc(mesh)
          INTEGER, POINTER :: lll(:)      ! lll(nbeta)
          INTEGER, POINTER :: kbeta(:)    ! kbeta(nbeta) 
          REAL(DP), POINTER :: beta(:,:)  ! beta(mesh,nbeta)
          INTEGER :: nd
          REAL(DP), POINTER :: dion(:,:)  ! dion(nbeta,nbeta)
          INTEGER :: nqf
          INTEGER :: nqlc
          REAL(DP), POINTER :: rinner(:)  ! rinner(0:2*lmax)
          REAL(DP), POINTER :: qqq(:,:)   ! qqq(nbeta,nbeta)
          REAL(DP), POINTER :: qfunc(:,:,:) ! qfunc(mesh,nbeta,nbeta)
          REAL(DP), POINTER :: qfcoef(:,:,:,:) ! qfcoef(nqf,0:2*lmax,nbeta,nbeta)
          REAL(DP), POINTER :: chi(:,:)   !  chi(mesh,nwfc)
          REAL(DP), POINTER :: rho_at(:)  !  rho_at(mesh)
          
          LOGICAL :: has_paw              ! Whether PAW data is included
          REAL(DP) :: paw_data_format     ! The version of the format
          LOGICAL :: has_gipaw            ! Whether GIPAW data is included
          REAL(DP) :: gipaw_data_format   ! The version of the format
          INTEGER :: gipaw_ncore_orbitals
          REAL(DP), POINTER :: gipaw_core_orbital_n(:)
          REAL(DP), POINTER :: gipaw_core_orbital_l(:)
          CHARACTER(LEN=2), POINTER :: gipaw_core_orbital_el(:)
          REAL(DP), POINTER :: gipaw_core_orbital(:,:)
          REAL(DP), POINTER :: gipaw_vlocal_ae(:)
          REAL(DP), POINTER :: gipaw_vlocal_ps(:)
          INTEGER :: gipaw_wfs_nchannels
          CHARACTER(LEN=2), POINTER :: gipaw_wfs_el(:)
          INTEGER, POINTER :: gipaw_wfs_ll(:)
          REAL(DP), POINTER :: gipaw_wfs_ae(:,:)
          REAL(DP), POINTER :: gipaw_wfs_rcut(:)
          REAL(DP), POINTER :: gipaw_wfs_rcutus(:)
          REAL(DP), POINTER :: gipaw_wfs_ps(:,:)
        END TYPE

      CONTAINS

        SUBROUTINE nullify_pseudo_upf( upf )
          TYPE( pseudo_upf ), INTENT(INOUT) :: upf
          NULLIFY( upf%els, upf%lchi, upf%jchi, upf%oc )
          NULLIFY( upf%r, upf%rab )  
          NULLIFY( upf%rho_atc, upf%vloc )  
          NULLIFY( upf%nn, upf%rcut)
          NULLIFY( upf%els_beta)
          NULLIFY( upf%rcutus, upf%epseu)
          NULLIFY( upf%lll, upf%jjj, upf%kbeta, upf%beta, upf%dion )  
          NULLIFY( upf%rinner, upf%qqq, upf%qfunc, upf%qfcoef )  
          NULLIFY( upf%chi )  
          NULLIFY( upf%rho_at )  
          NULLIFY ( upf%gipaw_core_orbital_n )
          NULLIFY ( upf%gipaw_core_orbital_l )
          NULLIFY ( upf%gipaw_core_orbital_el )
          NULLIFY ( upf%gipaw_core_orbital )
          NULLIFY ( upf%gipaw_vlocal_ae )
          NULLIFY ( upf%gipaw_vlocal_ps )
          NULLIFY ( upf%gipaw_wfs_el )
          NULLIFY ( upf%gipaw_wfs_ll )
          NULLIFY ( upf%gipaw_wfs_ae )
          NULLIFY ( upf%gipaw_wfs_rcut )
          NULLIFY ( upf%gipaw_wfs_rcutus )
          NULLIFY ( upf%gipaw_wfs_ps )
          RETURN
        END SUBROUTINE nullify_pseudo_upf

        SUBROUTINE deallocate_pseudo_upf( upf )
          TYPE( pseudo_upf ), INTENT(INOUT) :: upf
          IF( ASSOCIATED( upf%els ) ) DEALLOCATE( upf%els )
          IF( ASSOCIATED( upf%lchi ) ) DEALLOCATE( upf%lchi )
          IF( ASSOCIATED( upf%jchi ) ) DEALLOCATE( upf%jchi )
          IF( ASSOCIATED( upf%oc ) ) DEALLOCATE( upf%oc )
          IF( ASSOCIATED( upf%r ) ) DEALLOCATE( upf%r )
          IF( ASSOCIATED( upf%rab ) ) DEALLOCATE( upf%rab )
          IF( ASSOCIATED( upf%nn ) ) DEALLOCATE( upf%nn )
          IF( ASSOCIATED( upf%els_beta ) ) DEALLOCATE( upf%els_beta )
          IF( ASSOCIATED( upf%rcut ) ) DEALLOCATE( upf%rcut )
          IF( ASSOCIATED( upf%rcutus ) ) DEALLOCATE( upf%rcutus )
          IF( ASSOCIATED( upf%epseu ) ) DEALLOCATE( upf%epseu )
          IF( ASSOCIATED( upf%rho_atc ) ) DEALLOCATE( upf%rho_atc )
          IF( ASSOCIATED( upf%vloc ) ) DEALLOCATE( upf%vloc )
          IF( ASSOCIATED( upf%lll ) ) DEALLOCATE( upf%lll )
          IF( ASSOCIATED( upf%jjj ) ) DEALLOCATE( upf%jjj )
          IF( ASSOCIATED( upf%kbeta ) ) DEALLOCATE( upf%kbeta )
          IF( ASSOCIATED( upf%beta ) ) DEALLOCATE( upf%beta )
          IF( ASSOCIATED( upf%dion ) ) DEALLOCATE( upf%dion )
          IF( ASSOCIATED( upf%rinner ) ) DEALLOCATE( upf%rinner )
          IF( ASSOCIATED( upf%qqq ) ) DEALLOCATE( upf%qqq )
          IF( ASSOCIATED( upf%qfunc ) ) DEALLOCATE( upf%qfunc )
          IF( ASSOCIATED( upf%qfcoef ) ) DEALLOCATE( upf%qfcoef )
          IF( ASSOCIATED( upf%chi ) ) DEALLOCATE( upf%chi )
          IF( ASSOCIATED( upf%rho_at ) ) DEALLOCATE( upf%rho_at )
          IF ( ASSOCIATED ( upf%gipaw_core_orbital_n ) ) &
               DEALLOCATE ( upf%gipaw_core_orbital_n )
          IF ( ASSOCIATED ( upf%gipaw_core_orbital_l ) ) &
               DEALLOCATE ( upf%gipaw_core_orbital_l )
          IF ( ASSOCIATED ( upf%gipaw_core_orbital_el ) ) &
               DEALLOCATE ( upf%gipaw_core_orbital_el )
          IF ( ASSOCIATED ( upf%gipaw_core_orbital ) ) &
               DEALLOCATE ( upf%gipaw_core_orbital )
          IF ( ASSOCIATED ( upf%gipaw_vlocal_ae ) ) &
               DEALLOCATE ( upf%gipaw_vlocal_ae )
          IF ( ASSOCIATED ( upf%gipaw_vlocal_ps ) ) &
               DEALLOCATE ( upf%gipaw_vlocal_ps )
          IF ( ASSOCIATED ( upf%gipaw_wfs_el ) ) &
               DEALLOCATE ( upf%gipaw_wfs_el )
          IF ( ASSOCIATED ( upf%gipaw_wfs_ll ) ) &
               DEALLOCATE ( upf%gipaw_wfs_ll )
          IF ( ASSOCIATED ( upf%gipaw_wfs_ae ) ) &
               DEALLOCATE ( upf%gipaw_wfs_ae )
          IF ( ASSOCIATED ( upf%gipaw_wfs_rcut ) ) &
               DEALLOCATE ( upf%gipaw_wfs_rcut )
          IF ( ASSOCIATED ( upf%gipaw_wfs_rcutus ) ) &
               DEALLOCATE ( upf%gipaw_wfs_rcutus )
          IF ( ASSOCIATED ( upf%gipaw_wfs_ps ) ) &
               DEALLOCATE ( upf%gipaw_wfs_ps )
          RETURN
        END SUBROUTINE deallocate_pseudo_upf

      END MODULE pseudo_types
