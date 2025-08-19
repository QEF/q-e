!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE ldaU
  !--------------------------------------------------------------------------
  !! The quantities needed in DFT+U and extended DFT+U calculations.
  !
  USE kinds,         ONLY : DP
  USE upf_params,    ONLY : lqmax
  ! FIXME: lqmax should not be used (see starting_ns* below)
  USE parameters,    ONLY : ntypx, natx, sc_size
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE control_flags, ONLY : dfpt_hub
#if defined (__OSCDFT)
  USE plugin_flags,      ONLY : use_oscdft
  USE oscdft_base,       ONLY : oscdft_ctx
  USE oscdft_functions,  ONLY : oscdft_init_constraint
#endif 
  USE io_global,     ONLY : stdout
  !
  SAVE
  !
  COMPLEX(DP), ALLOCATABLE :: wfcU(:,:)
  !! atomic wfcs with U term - FIXME: should be PINNED
  COMPLEX(DP), ALLOCATABLE :: d_spin_ldau(:,:,:)
  !! the rotations in spin space for all symmetries
  REAL(DP) :: eth
  !! the Hubbard contribution to the energy
  REAL(DP) :: Hubbard_U(ntypx)
  !! the Hubbard U (main Hubbard channel)
  REAL(DP) :: Hubbard_U2(ntypx)
  !! the Hubbard U (second (and third) Hubbard channel)
  REAL(DP) :: Hubbard_Um(lqmax,2,ntypx)
  !! the (spin-)orbital-resolved Hubbard U
  REAL(DP) :: Hubbard_Um_nc(2*lqmax,ntypx)
  !! the (spin-)orbital-resolved Hubbard U (noncolinear case)
  REAL(DP) :: Hubbard_J0(ntypx)
  !! the Hubbard J, in simplified DFT+U
  REAL(DP) :: Hubbard_J(3,ntypx)
  !! extra Hubbard parameters for full DFT+U:  
  !! * p: J(1)=J  
  !! * d: J(1)=J, J(2)=B   
  !! * f: J(1)=J, J(2)=E2, J(3)=E3   
  REAL(DP) :: Hubbard_alpha(ntypx)
  !! the Hubbard alpha (used to calculate U)
  REAL(DP) :: Hubbard_alpha_back(ntypx)
  !! the Hubbard alpha (used to calculate U on background states)
  REAL(DP) :: Hubbard_alpha_m(lqmax,2,ntypx)
  !! the Hubbard alpha used to calculate orbital-resolved U parameters
  REAL(DP) :: Hubbard_alpha_m_nc(2*lqmax,ntypx)
  !! the Hubbard alpha used to calculate orbital-resolved U parameters (noncolinear case)
  REAL(DP) :: Hubbard_beta(ntypx)
  !! the Hubbard beta (used to calculate J0)
  REAL(DP) :: Hubbard_occ(ntypx,3)
  !! the Hubbard occupations
  REAL(DP) :: starting_ns(lqmax,2,ntypx)
  !! starting ns
  !! FIXME: allocate dynamically
  REAL(DP) :: starting_ns_back(lqmax,2,ntypx)
  !! starting ns on background states
  !! FIXME: allocate dynamically, or better, remove
  INTEGER :: nwfcU
  !! total no. of atomic wavefunctions having U term
  INTEGER :: niter_with_fixed_ns
  !! no. of iterations with fixed ns
  LOGICAL :: lda_plus_u
  !! .TRUE. if DFT+Hubbard is used
  INTEGER :: lda_plus_u_kind
  !! 0 --> Simplified rotationally-invariant formulation of DFT+U
  !! 1 --> Full formulation of DFT+U
  !! 2 --> Simplified rotationally-invariant formulation of DFT+U+V
  INTEGER :: Hubbard_n(ntypx)
  !! the principal quantum number of the first Hubbard manifold
  INTEGER :: Hubbard_l(ntypx)
  !! the angular momentum of the first Hubbard manifold
  INTEGER :: Hubbard_n2(ntypx)
  !! the principal quantum number of the second Hubbard manifold
  INTEGER :: Hubbard_l2(ntypx)
  !! the angular momentum of the second Hubbard manifold
  INTEGER :: Hubbard_n3(ntypx)
  !! the principal quantum number of the third Hubbard manifold
  INTEGER :: Hubbard_l3(ntypx)
  !! the angular momentum of the third Hubbard manifold
  INTEGER :: Hubbard_lmax = 0 
  !! maximum angular momentum of Hubbard states
  INTEGER :: Hubbard_lmax_back = 0
  !! maximum angular momentum of Hubbard background states
  INTEGER :: ldmx = -1
  !! max dimension of the manifold where the Hubbard correction
  !! is applied (max of 2*Hubbard_l+1 over all atoms)
  INTEGER :: ldmx_b = -1
  !! like above for background states
  INTEGER :: ldmx_tot
  !! max value of ldim_u(nt) = ldim_u(nt) + ldim_back(nt) over all ntyp
  INTEGER :: l0
  !! index in the array of Hubbard states. 
  !! for every atom one has 2*Hubbard_l+1 + 2*Hubbard_l2 +1 = ldim_u states for example.
  INTEGER :: l0b
  !! like above for background states
  LOGICAL :: is_hubbard(ntypx)
  !! .TRUE. if this atom species has U correction
  LOGICAL :: is_hubbard_back(ntypx)
  !! .TRUE. if this atom species has U correction for background states
  LOGICAL :: conv_ns
  !! .TRUE. if ns are converged
  LOGICAL :: reserv(ntypx), reserv_back(ntypx)
  !! reservoir states
  LOGICAL :: backall(ntypx)
  !! if .TRUE. two Hubbard manifolds can be used in the background (Hubbard_l2 and Hubbard_l3)
  LOGICAL :: hub_back
  !! .TRUE. if at least one species has Hubbard_U in the background states
  LOGICAL :: hub_pot_fix
  !! if .TRUE. do not include into account the change of the Hubbard potential
  !! during the SCF cycle (needed to compute U self-consistently with supercells)
  LOGICAL :: orbital_resolved
  !! if .TRUE., the orbital-resolved formulation of Hubbard corrections is used
  !! currently only affects Hubbard U and Hubbard alpha
  LOGICAL :: apply_U
  !!  In one-step orbital-resolved DFT+U calculations, once set to .TRUE.,
  !!  the Hubbard potential and energy are calculated based on the diagonalized 
  !!  occupations. Before, Hubbard corrections are not applied in order
  !!  to stabilize the eigenstates before storing reference eigenvectors.
  LOGICAL :: iso_sys
  !! .TRUE. if the system is isolated (the code diagonalizes
  !! and prints the full occupation matrix)
  CHARACTER(len=30) :: Hubbard_projectors
  !! 'atomic', 'ortho-atomic', 'file'
  INTEGER, ALLOCATABLE :: oatwfc(:)
  !! specifies how input coordinates are given
  INTEGER, ALLOCATABLE :: oatwfc_back(:), oatwfc_back1(:)
  !! specifies how input coordinates are given for background states
  INTEGER, ALLOCATABLE :: offsetU(:) 
  !! offset of atomic wfcs used for projections
  INTEGER, ALLOCATABLE :: offsetU_back(:), offsetU_back1(:)
  !! offset of atomic wfcs used for projections (background states)
  INTEGER, ALLOCATABLE :: ldim_u(:)
  !! number of U states for each atom
  INTEGER, ALLOCATABLE :: ldim_back(:)
  !! number of U background states for each atom
  INTEGER, ALLOCATABLE :: ll(:,:)
  !! ll(i=1:ldim_u) = Hubbard_l   if 0    <= i <= l0
  !!                = Hubbard_l2  if l0+1 <  i <= ldim_u
  REAL(DP), ALLOCATABLE :: q_ae(:,:,:)
  !! coefficients for projecting onto beta functions
  REAL(DP), ALLOCATABLE :: q_ps(:,:,:)
  !! (matrix elements on AE and PS atomic wfcs)
  REAL(DP), ALLOCATABLE :: lambda_ns(:,:,:)
  !! Array to store the eigenvalues of the occupation matrix,
  !! needed for orbital-resolved DFT+U
  COMPLEX(DP), ALLOCATABLE :: eigenvecs_ref(:,:,:,:)
  !! Eigenvectors used for tracking
  !! orbitals in orbital-resolved DFT+U
  INTEGER, ALLOCATABLE :: order_um(:,:,:) 
  !!
  !***************************************************
  !                Hubbard V part                     !
  !****************************************************
  !
  ! Inter atomic interaction should be cut off at some distance 
  !
  ! Inter atomic interaction should be cut off at some distance 
  ! that is the reason of having so many unitcell information. 
  !
  REAL(DP) :: Hubbard_V(natx,natx*(2*sc_size+1)**3,4) 
  !! The Hubbard_V(I,J,int_type) gives the interaction between atom I (in the unit cell)
  !! with atom J (in the supercell).
  !! If int_type=1, the interaction is between standard orbitals,
  !! If int_type=2, the interaction is between standard (on I) and background (on J) orbitals,
  !! If int_type=3, the interaction is between background orbitals,
  !! If int_type=4, the interaction is between background (on I) and standard (on J) orbitals.
  !! Hubbard_V(I,J,4) is equal to Hubbard_V(J,I,2). It is useful
  !! in cases where Hubbard_V(I,J,2) /= 0 but I is outside the unit cell, J inside.
  INTEGER :: num_uc
  !! Number of unit cells in the supercell = (2*sc_size+1)**3
  INTEGER :: max_num_neighbors
  !! the maximum number of neighbors
  REAL(DP), ALLOCATABLE :: atom_pos(:,:)
  !! matrix with dimensions nat x 3 with the atomic coordinates in 
  !! the primitive basis.
  REAL(DP), ALLOCATABLE :: dist_s(:,:)
  !! Distance between atoms in the 3x3x3 supercell (if sc_size = 1)
  INTEGER,  ALLOCATABLE :: ityp_s(:)
  !! Type of atoms in the 3x3x3 supercell (if sc_size = 1)
  REAL(DP), ALLOCATABLE :: nsnew(:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: nsg(:,:,:,:,:), nsgnew(:,:,:,:,:)
  !! Generalized occupation matrices, which depend on two atomic sites.
  !! These matrices nsg(at1,m1,viz,m2,sp) store the expectation value:
  !! <C^\dagger_{at1,m1,sp}C_{viz,m2,sp}>, where sp = spin and
  !! viz identifies the atom in the neighborhood of at1.
  COMPLEX(DP), ALLOCATABLE :: v_nsg(:,:,:,:,:)
  !! The kernel of the Hubbard potential (see above for the meaning of the
  !! size of the array)
  COMPLEX(DP), ALLOCATABLE :: phase_fac(:)
  !! Phase factor (it is 1 if we have only Hubbard U)
  INTEGER, ALLOCATABLE :: sc_at(:,:,:,:)
  !! Matrix with ranges [1:nat], gives the corresponding atom in the supercell ordering 
  REAL(DP), PARAMETER :: eps_dist = 6.d-4     
  !! Threshold for comparing inter-atomic distances
  !
  TYPE position
     INTEGER :: at, n(3)
     !! Identifies the atom: it is equivalent to atom 'at' in the unit cell,
     !! and it is located in the unit cell (n(1),n(2),n(3))    
  ENDTYPE position                       
  !
  TYPE at_center
     integer :: num_neigh                  
     !! Number of neighbours (i.e., other V-interacting atoms) 
     !! of the considered atom
     integer, ALLOCATABLE :: neigh(:)      
     !! Vector indicating which are the neigbours
  ENDTYPE at_center
  !
  TYPE(position), ALLOCATABLE :: at_sc(:)         
  !! Vector with all the atoms in the supercell
  !
  TYPE(at_center), ALLOCATABLE :: neighood(:)     
  !! Vector with the information about the neighbours 
  !! for all the atoms in the unit cell
  !
CONTAINS
  !
  SUBROUTINE init_hubbard ( psd, nspin, noncolin )
    !
    !! NOTE: average_pp must be called before init_hubbard
    !
    IMPLICIT NONE
    !
    CHARACTER (LEN=2), INTENT(IN) :: psd(:)
    INTEGER, INTENT(IN) :: nspin
    LOGICAL, INTENT(IN) :: noncolin
    !
    INTEGER :: na, nt
    LOGICAL :: lba, lb
    !
    lba = .FALSE.
    lb  = .FALSE.
    hub_back = .FALSE.
    hub_pot_fix = .FALSE.
    orbital_resolved = .FALSE.
    apply_U = .FALSE.
    !
    is_hubbard(:) = .FALSE.
    is_hubbard_back(:) = .FALSE.
    !
    IF ( .NOT. lda_plus_u ) THEN
       Hubbard_lmax = 0
       RETURN
    ENDIF
    !
    Hubbard_lmax = -1
    Hubbard_lmax_back  = -1
    ldmx_tot = -1
    !
    IF (.NOT.ALLOCATED (ldim_u) )    ALLOCATE(ldim_u(ntyp))
    ldim_u(:)=-1
    IF (.NOT.ALLOCATED (ldim_back) ) ALLOCATE(ldim_back(ntyp))
    ldim_back(:)=-1
    !
    IF ( lda_plus_u_kind == 0 ) THEN
       !
       ! DFT+U(+J0) : Dudarev's formulation
       !
       DO nt = 1, ntyp
          !
          is_hubbard(nt) = Hubbard_U(nt) /= 0.0_DP                   .OR. &
                           Hubbard_U2(nt)/= 0.0_dp                   .OR. &
                           ANY(Hubbard_Um(:,:,nt) /= 0.0_DP)         .OR. &
                           ANY(Hubbard_Um_nc(:,nt) /= 0.0_DP)        .OR. &
                           Hubbard_alpha(nt) /= 0.0_DP               .OR. &
                           Hubbard_alpha_back(nt) /= 0.0_dp          .OR. &
                           ANY(Hubbard_alpha_m(:,:,nt) /= 0.0_DP)    .OR. &
                           ANY(Hubbard_alpha_m_nc(:,nt) /= 0.0_DP)   .OR. &
                           Hubbard_J0(nt) /= 0.0_DP                  .OR. &
                           Hubbard_beta(nt) /= 0.0_DP  
          !
          is_hubbard_back(nt) = Hubbard_U2(nt)/= 0.0_dp     .OR. &
                                Hubbard_alpha_back(nt) /= 0.0_dp                 
          !
          orbital_resolved =  orbital_resolved .OR. &
                              ANY(Hubbard_Um(:,:,nt) /= 0.0_DP)       .OR. &
                              ANY(Hubbard_Um_nc(:,nt) /= 0.0_DP)      .OR. &
                              ANY(Hubbard_alpha_m(:,:,nt) /= 0.0_DP)  .OR. &
                              ANY(Hubbard_alpha_m_nc(:,nt) /= 0.0_DP)
          !
          IF ( is_hubbard(nt) ) THEN
             ! Hubbard_l is read from the input file (HUBBARD card)
             Hubbard_lmax = MAX( Hubbard_lmax, Hubbard_l(nt) )
             ldmx = MAX( ldmx, 2*Hubbard_l(nt)+1 )
             ldim_u(nt) = 2*Hubbard_l(nt)+1
             IF (hubbard_occ(nt,1)<0.0d0) CALL determine_hubbard_occ(nt,1) 
          ENDIF
          !
          IF ( is_hubbard_back(nt) ) THEN
             lb = .TRUE.
             hub_back = .TRUE.
             !
             IF (.NOT.backall(nt)) THEN
                ! In this case there is one Hubbard channel for background states.
                ! Hubbard_l2 is read from the input file (HUBBARD card)
                Hubbard_lmax_back = MAX( Hubbard_lmax_back, Hubbard_l2(nt) )
                ldmx_b = MAX( ldmx_b, 2*Hubbard_l2(nt)+1)
                ldim_back(nt) = 2 * Hubbard_l2(nt) + 1
                IF (hubbard_occ(nt,2)<0.0d0) CALL determine_hubbard_occ(nt,2)
             ELSE
                ! In this case there are two Hubbard channels for background states.
                ! Hubbard_l2 and Hubbard_l3 are read from the input file (HUBBARD card)
                ! Note: Hubbard_U2 is used for these two background channels.
                lba = .TRUE.
                Hubbard_lmax_back = MAX( Hubbard_lmax_back, Hubbard_l3(nt) )
                ldmx_b = MAX( ldmx_b, 2*Hubbard_l2(nt)+2*Hubbard_l3(nt)+2 )
                ldim_back(nt) = 2 * (Hubbard_l2(nt) + Hubbard_l3(nt) + 1)
                IF (hubbard_occ(nt,2)<0.0d0) CALL determine_hubbard_occ(nt,2)
                IF (hubbard_occ(nt,3)<0.0d0) CALL determine_hubbard_occ(nt,3)
             ENDIF
             ldim_u(nt) = ldim_u(nt) + ldim_back(nt) 
             Hubbard_lmax_back = MAX( Hubbard_lmax_back, Hubbard_l2(nt) )
          ENDIF
          !
          ldmx_tot = MAX( ldmx_tot, ldim_u(nt) )
          !
       ENDDO !nt
       !
       IF (orbital_resolved) THEN
          IF (noncolin) THEN
             ! need to store eigenvectors and eigenvalues in a 2*ldim array
             ! retain extra-spin dimension for compatibility
             ALLOCATE(lambda_ns(2*ldmx,1,nat))
             ALLOCATE(eigenvecs_ref(2*ldmx,2*ldmx,1,nat))
          ELSE
             ALLOCATE(lambda_ns(ldmx,nspin,nat))
             ALLOCATE(eigenvecs_ref(ldmx,ldmx,nspin,nat))
          ENDIF
          !
          lambda_ns(:,:,:) = 0.0_DP
          eigenvecs_ref(:,:,:,:) = CMPLX(0.d0,0.d0, kind=DP)
       ENDIF
       !
       IF ( ANY(Hubbard_alpha(1:ntyp) /= 0.0_DP)       .OR. &
            ANY(Hubbard_alpha_m(:,:,1:ntyp) /= 0.0_DP) .OR. & 
            ANY(Hubbard_alpha_m_nc(:,1:ntyp) /= 0.0_DP)       ) THEN
          ! To apply LR-cDFT to calculate Hubbard parameters,
          ! we fix the Hubbard potential. Also, if this is an orbital-
          ! resolved calculation, turn on corrections before the first iteration.
          ! This supposes we are restarting from a converged charge density, which
          ! is checked in potinit.f90.
          hub_pot_fix = .TRUE.
          IF (orbital_resolved) apply_U = .TRUE.
          WRITE(stdout,'(/5x,"NONZERO HUBBARD_ALPHA DETECTED:")')
          WRITE(stdout,'(/5x,"FIXING HUBBARD POTENTIAL TO THE &
                              &GROUND STATE ONE (PRB 98, 085127)")')
       ENDIF
       !
       IF (ALLOCATED(ll)) DEALLOCATE (ll) 
       ALLOCATE(ll(ldmx_tot,ntyp))
       !
       ! ll is a label of all the Hubbard states telling the l of that states. 
       ! It is equal to Hubbard_l for the first 2*Hubbard_l+1 states, 
       ! to Hubbard_l2 for the next 2*Hubbard_l2+1, 
       ! to Hubbard_l3 for the next 2*Hubbard_l3+1
       ! (if there are two channels in the background).
       ! 
       ll(:,:) = -1
       DO nt = 1,ntyp
          IF (Hubbard_l(nt).GE.0) THEN
             ll(1:2*Hubbard_l(nt)+1,nt) = Hubbard_l(nt)
             l0 = 2*Hubbard_l(nt)+1
          ELSE
             l0 = 0
          ENDIF
          IF (Hubbard_l2(nt).GE.0) THEN
             ll(l0+1:l0+2*Hubbard_l2(nt)+1,nt) = &
             Hubbard_l2(nt)
             l0b = l0 + 2*Hubbard_l2(nt)+1
          ELSE
             l0b = 0
          ENDIF
          IF (backall(nt) .AND. Hubbard_l3(nt).GE.0) THEN
             ll(l0b+1:l0b+2*Hubbard_l3(nt)+1,nt) = &
             Hubbard_l3(nt)
          ENDIF
       ENDDO  
       IF (noncolin) THEN
          IF ( .NOT. ALLOCATED (d_spin_ldau) ) ALLOCATE( d_spin_ldau(2,2,48) )
          CALL comp_dspinldau()
       ENDIF
    ELSEIF ( lda_plus_u_kind == 1 ) THEN
       !
       ! DFT+U(+J) : Liechtenstein's formulation
       ! 
       IF ( Hubbard_projectors == 'pseudo' ) CALL errore( 'init_hubbard', &
            & 'full DFT+U not implemented with pseudo Hubbard projectors', 1 )
       !
       IF (noncolin) THEN
          IF ( .NOT. ALLOCATED (d_spin_ldau) ) ALLOCATE( d_spin_ldau(2,2,48) )
          CALL comp_dspinldau()
       ENDIF
       !
       DO nt = 1, ntyp
          IF (Hubbard_alpha(nt)/=0.d0 ) CALL errore( 'init_hubbard', &
               'full DFT+U does not support Hubbard_alpha calculation', 1 )

          is_hubbard(nt) = Hubbard_U(nt)/= 0.0_dp .OR. &
                           ANY( Hubbard_J(:,nt) /= 0.0_dp )

          IF ( is_hubbard(nt) ) THEN
             !
             ! Hubbard_l is read from the input file (HUBBARD card)
             Hubbard_lmax = MAX( Hubbard_lmax, Hubbard_l(nt) )
             ldmx = MAX( ldmx, 2*Hubbard_l(nt)+1 )
             ldim_u(nt) = 2*Hubbard_l(nt)+1
             IF (hubbard_occ(nt,1)<0.0d0) CALL determine_hubbard_occ(nt,1)
             !
             IF (Hubbard_U(nt) == 0.0_dp) Hubbard_U(nt) = 1.d-14
             !
             IF ( Hubbard_l(nt) == 2 ) THEN
                IF ( Hubbard_J(2,nt) == 0.d0 ) &
                     Hubbard_J(2,nt) = 0.114774114774d0 * Hubbard_J(1,nt)
             ELSEIF ( Hubbard_l(nt) == 3 ) THEN
                IF ( Hubbard_J(2,nt) == 0.d0 ) &
                     Hubbard_J(2,nt) = 0.002268d0 * Hubbard_J(1,nt)
                IF ( Hubbard_J(3,nt)==0.d0 )   &
                     Hubbard_J(3,nt) = 0.0438d0 * Hubbard_J(1,nt)
             ENDIF
          ENDIF
          !
       ENDDO
       !
    ELSEIF ( lda_plus_u_kind == 2 ) THEN
       !
       ! DFT+U+V(+J0) : Dudarev's formulation
       !
       ! Number of cells in the supercell
       num_uc = (2*sc_size+1)**3
       !
       ! Setup atomic positions in the primitive basis coordinates
       !
       CALL alloc_atom_pos()
       !
       ! In this case is_hubbard is set inside alloc_neighborhood
       !
       CALL alloc_neighborhood()
       !
       ldmx = 0
       ldmx_b = 0
       lba = .FALSE.
       !
       DO nt = 1, ntyp
          !
          IF (noncolin) THEN
             IF ( .NOT. ALLOCATED (d_spin_ldau) ) ALLOCATE( d_spin_ldau(2,2,48) )
             CALL comp_dspinldau()
          ENDIF
          !
          ! Here we account for the remaining cases when we need to 
          ! setup is_hubbard
          !
          is_hubbard(nt) = is_hubbard(nt)              .OR. &
                           Hubbard_alpha(nt) /= 0.0_dp .OR. &
                           Hubbard_beta(nt) /= 0.0_dp  .OR. &
                           Hubbard_J0(nt) /= 0.0_dp
          !
          is_hubbard_back(nt) = is_hubbard_back(nt)    .OR. &
                                Hubbard_alpha_back(nt) /= 0.0_dp
          !
          IF (  is_hubbard(nt) ) THEN
             ! Hubbard_l is read from the input file (HUBBARD card)
             Hubbard_lmax = MAX( Hubbard_lmax, Hubbard_l(nt) )
             ldmx = MAX( ldmx, 2*Hubbard_l(nt)+1 )
             ldim_u(nt) = 2*Hubbard_l(nt)+1
             IF (hubbard_occ(nt,1)<0.0d0) CALL determine_hubbard_occ(nt,1)
          ENDIF
          !
          IF ( is_hubbard_back(nt) ) THEN
             !
             lb = .TRUE.
             hub_back = .TRUE.
             !
             IF (.NOT.backall(nt)) THEN
                ! In this case there is only one Hubbard channel for background states
                ! Hubbard_l2 is read from the input file (HUBBARD card)
                ldim_back(nt) = 2 * Hubbard_l2(nt) + 1
                IF (hubbard_occ(nt,2)<0.0d0) CALL determine_hubbard_occ(nt,2)
             ELSE
                ! In this case there are two Hubbard channels for background states.
                ! Hubbard_l2 and Hubbard_l3 are read from the input file (HUBBARD card)
                ! Note: Hubbard_U2 is used for these two background channels.
                lba = .TRUE.
                ldim_back(nt) = 2 * Hubbard_l2(nt) + 2 * Hubbard_l3(nt) + 2
                Hubbard_lmax_back = MAX( Hubbard_lmax_back, Hubbard_l3(nt) )
                IF (hubbard_occ(nt,2)<0.0d0) CALL determine_hubbard_occ(nt,2)
                IF (hubbard_occ(nt,3)<0.0d0) CALL determine_hubbard_occ(nt,3)
             ENDIF
             !
             ldim_u(nt) = ldim_u(nt) + ldim_back(nt) 
             ldmx_b = MAX( ldmx_b, ldim_back(nt) )
             Hubbard_lmax_back = MAX( Hubbard_lmax_back, Hubbard_l2(nt) )
             ! 
          ENDIF
          !
          ldmx_tot = MAX( ldmx_tot, ldim_u(nt) )
          !
       ENDDO
       !
       ! The allocation should be moved into scf_mod ?
       !
       IF (.NOT.ALLOCATED(v_nsg))     ALLOCATE ( v_nsg ( ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin ) )
       IF (.NOT.ALLOCATED(nsg))       ALLOCATE ( nsg   ( ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin ) )
       IF (.NOT.ALLOCATED(nsgnew))    ALLOCATE ( nsgnew( ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin ) )
       IF (.NOT.ALLOCATED(phase_fac)) ALLOCATE ( phase_fac(nat*num_uc))
       IF (.NOT.ALLOCATED(ll))        ALLOCATE ( ll(ldmx_tot, ntyp))
       !
       ! ll is a label of all the Hubbard states telling the l of that states. 
       ! It is equal to Hubbard_l for the first 2*Hubbard_l+1 states, 
       ! to Hubbard_l2 for the next 2*Hubbard_l2+1, 
       ! to Hubbard_l3 for the next 2*Hubbard_l3+1
       ! (if there are two channels in the background).
       !
       ll(:,:) = -1
       !
       DO nt = 1, ntyp
          !
          IF (Hubbard_l(nt).GE.0) THEN
             ll(1:2*Hubbard_l(nt)+1,nt) = Hubbard_l(nt)
             l0 = 2*Hubbard_l(nt)+1
          ELSE
             l0 = 0
          ENDIF
          !
          IF (Hubbard_l2(nt).GE.0) THEN
             ll(l0+1:l0+2*Hubbard_l2(nt)+1,nt) = Hubbard_l2(nt)
             l0b = l0 + 2*Hubbard_l2(nt)+1
          ELSE
             l0b = 0
          ENDIF
          !
          IF (backall(nt) .AND. Hubbard_l3(nt).GE.0) &
             & ll(l0b+1:l0b+2*Hubbard_l3(nt)+1,nt) = Hubbard_l3(nt)
          ! 
       ENDDO
       !
    ELSE
       !
       CALL errore( 'init_hubbard', 'Not allowed value of lda_plus_u_kind', 1 )
       !
    ENDIF
    !
#if defined (__OSCDFT)
     IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
        CALL oscdft_init_constraint (oscdft_ctx, ldmx, starting_ns, lda_plus_u, lda_plus_u_kind, is_hubbard)
     ENDIF
#endif
    !
    IF ( Hubbard_lmax == -1 ) CALL errore( 'init_hubbard', &
         'DFT+Hubbard calculation but Hubbard_l not set', 1 )
    !
    IF ( Hubbard_lmax > 3 ) &
         CALL errore( 'init_hubbard', 'Hubbard_l should not be > 3 ', 1 )
    !
    ! Compute the index of atomic wfcs used as projectors
    !
    IF ( .NOT.ALLOCATED(oatwfc)) ALLOCATE( oatwfc(nat) )
    CALL offset_atom_wfc ( .FALSE., 1, oatwfc, nwfcU )
    !
    IF ( lb .AND. .NOT.ALLOCATED(oatwfc_back)) THEN
       ALLOCATE ( oatwfc_back(nat) )
       CALL offset_atom_wfc ( .FALSE., 2, oatwfc_back, nwfcU )
    ENDIF
    !
    IF ( lba .AND. .NOT.ALLOCATED(oatwfc_back1)) THEN
       ALLOCATE ( oatwfc_back1(nat) )
       CALL offset_atom_wfc ( .FALSE., 3, oatwfc_back1, nwfcU )
    ENDIF
    !
    ! nwfcU is set to natomwfc by the routine above
    ! check below disabled because it introduces a dependency upon natomwfc
    ! IF ( nwfcU /= natomwfc ) &
    !     CALL errore( 'offset_atom_wfc', 'wrong number of wavefunctions', 1 )
    !
    ! For each atom, compute the index of its projectors (among projectors only)
    !
    IF ( .NOT.ALLOCATED(offsetU)) ALLOCATE( offsetU(nat) )
    CALL offset_atom_wfc( .TRUE., 1, offsetU, nwfcU )
    !
    IF ( lb .AND. .NOT.ALLOCATED(offsetU_back)) THEN
       ALLOCATE ( offsetU_back(nat) )
       CALL offset_atom_wfc ( .TRUE., 2, offsetU_back, nwfcU )
    ENDIF
    !
    IF ( lba .AND. .NOT.ALLOCATED(offsetU_back1)) THEN
       ALLOCATE ( offsetU_back1(nat) )
       CALL offset_atom_wfc ( .TRUE., 3, offsetU_back1, nwfcU )
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE init_hubbard
  !
  SUBROUTINE deallocate_hubbard( flag )
  !
  LOGICAL, INTENT(IN) :: flag
  INTEGER :: na
  !
  IF ( flag ) THEN
     IF ( ALLOCATED( oatwfc ) )        DEALLOCATE( oatwfc )
     IF ( ALLOCATED( oatwfc_back ) )   DEALLOCATE( oatwfc_back )
     IF ( ALLOCATED( oatwfc_back1 ) )  DEALLOCATE( oatwfc_back1 )
     IF ( ALLOCATED( offsetU ) )       DEALLOCATE( offsetU )
     IF ( ALLOCATED( offsetU_back ) )  DEALLOCATE( offsetU_back )
     IF ( ALLOCATED( offsetU_back1 ) ) DEALLOCATE( offsetU_back1 )
     IF ( ALLOCATED( q_ae ) )          DEALLOCATE( q_ae )
     IF ( ALLOCATED( q_ps ) )          DEALLOCATE( q_ps )
     IF ( ALLOCATED( d_spin_ldau ))    DEALLOCATE( d_spin_ldau )
     IF ( ALLOCATED( ll ) )            DEALLOCATE( ll )
     IF ( ALLOCATED( v_nsg ) )         DEALLOCATE( v_nsg )
     IF ( ALLOCATED( nsg ) )           DEALLOCATE( nsg )
     IF ( ALLOCATED( nsgnew ) )        DEALLOCATE( nsgnew )
     IF ( ALLOCATED( phase_fac ) )     DEALLOCATE( phase_fac )
     IF ( ALLOCATED( atom_pos ) )      DEALLOCATE( atom_pos )
     IF ( ALLOCATED( at_sc ) )         DEALLOCATE( at_sc )
     IF ( ALLOCATED( sc_at ) )         DEALLOCATE( sc_at )
     IF ( ALLOCATED( lambda_ns ) )     DEALLOCATE( lambda_ns )
     IF ( ALLOCATED( eigenvecs_ref))   DEALLOCATE( eigenvecs_ref )
     IF ( ALLOCATED( neighood ) ) THEN
        DO na = 1, nat
           CALL deallocate_at_center_type ( neighood(na) )
        ENDDO
        DEALLOCATE( neighood )
     ENDIF
     IF ( ALLOCATED( ldim_u ) )        DEALLOCATE( ldim_u )
     IF ( ALLOCATED( ldim_back ) )     DEALLOCATE( ldim_back )
  END IF
  !
  IF ( ALLOCATED( wfcU ) ) THEN
     !$acc exit data delete(wfcU)
     DEALLOCATE( wfcU )
  END IF
  !
  IF (.NOT.dfpt_hub) THEN
     IF ( ALLOCATED( dist_s ) )        DEALLOCATE( dist_s )
     IF ( ALLOCATED( ityp_s ) )        DEALLOCATE( ityp_s )
  ENDIF
  !
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     IF ( ALLOCATED( oscdft_ctx%constraining ) ) DEALLOCATE( oscdft_ctx%constraining )
     IF ( ALLOCATED( oscdft_ctx%constraint ) ) DEALLOCATE( oscdft_ctx%constraint )
  ENDIF
#endif
  !
  RETURN
  !
  END SUBROUTINE deallocate_hubbard
  !
  SUBROUTINE deallocate_at_center_type (neighood_)
  !
  IMPLICIT NONE
  TYPE (at_center) :: neighood_
  !
  neighood_%num_neigh = 0
  !
  IF ( ALLOCATED( neighood_%neigh ) )  DEALLOCATE( neighood_%neigh )
  !
  RETURN
  !
  END SUBROUTINE deallocate_at_center_type
  !
  SUBROUTINE copy_U_wfc( swfcatom, noncolin )
  !
  !  Copy (orthogonalized) atomic wavefunctions "swfcatom"
  !  having a Hubbard U correction to array "wfcU"
  !
  IMPLICIT NONE
  COMPLEX(KIND=DP), INTENT(IN) :: swfcatom(:,:)
  LOGICAL, INTENT(IN), OPTIONAL :: noncolin
  LOGICAL :: twice
  INTEGER :: na, nt, m1, m2, offU, offa

  IF ( PRESENT(noncolin) ) THEN
     twice = noncolin
  ELSE
     twice = .FALSE.
  ENDIF
  !$acc data present_or_copyin(swfcatom) present_or_copyout(wfcU)
  DO na = 1, nat
     nt = ityp(na)
     IF ( is_hubbard(nt) ) THEN
        m1 = 1
        m2 = 2*hubbard_l(nt)+1
        IF ( twice ) m2 = 2*m2
        offU = offsetU(na)
        offa = oatwfc(na)
        !$acc kernels
        wfcU(:,offU+m1:offU+m2) = swfcatom(:,offa+m1:offa+m2)
        !$acc end kernels
     ENDIF
     IF (is_hubbard_back(nt)) THEN
        m1 = 1
        m2 = 2*Hubbard_l2(nt)+1
        offU = offsetU_back(na)
        offa = oatwfc_back(na)
        !$acc kernels
        wfcU(:,offU+m1:offU+m2) = swfcatom(:,offa+m1:offa+m2)
        !$acc end kernels
        IF (backall(nt)) THEN
           m1 = 1
           m2 = 2*Hubbard_l3(nt)+1
           offU = offsetU_back1(na)
           offa = oatwfc_back1(na)
           !$acc kernels
           wfcU(:,offU+m1:offU+m2) = swfcatom(:,offa+m1:offa+m2)
           !$acc end kernels
        ENDIF
     ENDIF
  ENDDO
  !$acc end data
  RETURN
  !
  END SUBROUTINE copy_U_wfc
  !
END MODULE ldaU
!
