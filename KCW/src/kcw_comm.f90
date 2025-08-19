!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
!! Common variables for the KCW program
!
!
MODULE control_kcw
  !
  USE kinds,        ONLY : DP
  !
  SAVE
  CHARACTER(LEN=256)       :: tmp_dir_kcw, tmp_dir_kcwq, tmp_dir_save
  !
  COMPLEX (DP), ALLOCATABLE:: unimatrx(:,:,:)           ! Unitary matrix Canonical-variational OCC orbitals (n,n,kpoint)
  COMPLEX (DP), ALLOCATABLE:: unimatrx_opt(:,:,:)       ! Optimal Unitary matrix from EMP disentanglement
  COMPLEX (DP), ALLOCATABLE:: evc0(:,:)                 ! Variational orbitals at agiven kpoint
  REAL(DP), ALLOCATABLE    :: alpha_final(:)            ! Value of the screening parameter
  REAL(DP), ALLOCATABLE    :: alpha_final_full(:)       ! Value of the screening parameter in case the full hamiltonian is computed
  LOGICAL                  :: l_unique_manifold         ! IF true no occ/empty distinction from Wannier
  !
  LOGICAL :: kcw_at_ks                                   ! if TRUE calculate alphas for KS orbitals
  LOGICAL :: fix_orb                                    ! if TRUE freeze the response of the orbital we are looking at
  LOGICAL :: homo_only                                  ! if TRUE only screening parameter for HOMO eigenvalue is computed
  ! 
  REAL(DP) :: spread_thr                                ! the tollerance within which two orbital are considered to have the same spread
  !
  LOGICAL  :: read_unitary_matrix                       ! if TRUE read the unitary matrix that defines the minimizing orbital
  LOGICAL  :: has_disentangle                           ! if TRUE read the unitary matrix that defines the minimizing orbital
  LOGICAL  :: have_empty                                ! if TRUE read the U matrices fore empty states 
  LOGICAL  :: write_pot                                 ! if TRUE write bare and screened potential to files
  LOGICAL  :: qp_symm                                   ! if TRUE build and hermitin operator for the KI hamiltonian. This is in hte spirit of scQP GW.
  LOGICAL  :: kipz_corr                                 ! if TRUE add the PZ term of KIPZ. This is in the spirit of a perturbative approach 
  CHARACTER(len=256) :: seedname                        ! the seedname of the Wannier calculation
  INTEGER  :: num_wann_occ                              ! number of wannier function for the OCC manifold
  INTEGER  :: num_wann_emp                              ! number of wannier function for the EMP manifold
  INTEGER  :: num_wann                                  ! total number of wannier function num_wann = num_wann_occ + num_wann_emp
  INTEGER  :: iuwfc_wann, iurho_wann, iuwfc_wann_allk   ! units for the wannier WFC and densities
  LOGICAL  :: check_ks                                  ! compute the KS Hamiltonian on top of the localized representation
  LOGICAL  :: do_bands                                  ! if TRUE computes electronic bands
  LOGICAL  :: use_ws_distance                           ! if TRUE uses Wannier centers in the interpolation
  LOGICAL  :: write_hr                                  ! if TRUE prints KC H(R) to file
  INTEGER,  ALLOCATABLE   :: map_ikq(:)                 ! map k+q -> p+G 
  INTEGER,  ALLOCATABLE   :: map_ikq_minus(:)           ! map k-q -> p+G 
  REAL(DP), ALLOCATABLE :: shift_1bz(:,:)               ! the G to brings k+q inside the 1BZ
  REAL(DP), ALLOCATABLE :: shift_1bz_minus(:,:)         ! the G to brings k-q inside the 1BZ
  INTEGER :: kcw_iverbosity                              ! set the verbosity 
  INTEGER :: spin_component                             ! which spin component to compute 
  INTEGER, ALLOCATABLE :: isq(:)                        ! the spin of the q point
  INTEGER, PARAMETER :: maxter = 100                    ! maximum number of iterations
  REAL (DP), ALLOCATABLE :: dmuxc_pol(:,:,:) 
  INTEGER :: i_orb                                      ! orbital to compute given from input 
  INTEGER :: iorb_start, iorb_end                       
  REAL (DP), ALLOCATABLE :: wq(:)                       ! weight of the point q
  INTEGER :: nqstot                                     ! total number of q points 
  INTEGER :: nrho
  INTEGER :: nkstot_eff
  !
  CHARACTER(len=256) :: calculation                     ! set which calculation (screen or ham)
  !
  !! The Hamiltonian
  LOGICAL :: on_site_only               ! On-site only: only R=0 and i=j
  REAL(DP), ALLOCATABLE :: rvect(:,:)   ! The R vector in the (virtual) supercell
  INTEGER, ALLOCATABLE :: irvect(:,:)   ! The R vector in the (virtual) supercell
  INTEGER :: mp1, mp2, mp3              ! the MP grid 
  !
  COMPLEX(DP), ALLOCATABLE :: Hamlt(:,:,:)   ! the KC hamiltonian 
  COMPLEX(DP), ALLOCATABLE :: Hamlt_R(:,:,:) ! the KC hamiltonian in WF representation
  !
  REAL(DP), ALLOCATABLE :: xk_bands(:,:), wk_bands(:)
  INTEGER :: nks_bands
  REAL (DP) :: eps_inf    ! The macroscopic dielectric funtion
  LOGICAL :: l_vcut       ! IF true use a regularization for the coulomb potential. set (q+G)0=0 otherwise
  LOGICAL :: l_alpha_corr ! add a correction to the second order approximation of the energy based 
                          ! on the un-relaxed energy at N \pm 1
  LOGICAL, ALLOCATABLE :: alpha_corr_done(:) ! keep track wheter the correction to alpha was already conputed
  !
  REAL(DP), ALLOCATABLE :: delta_alpha(:)
  !
  REAL(DP), ALLOCATABLE :: centers(:,:)    ! Wannier centers
  !
  REAL(DP), ALLOCATABLE :: occ_mat(:,:,:)
  !
  LOGICAL, ALLOCATABLE :: l_do_alpha(:)
  INTEGER, ALLOCATABLE :: group_alpha(:)
  LOGICAL :: check_spread 
  !
  REAL(DP), ALLOCATABLE :: x_q(:,:)
  INTEGER               :: niter
  INTEGER               :: nmix
  REAL(DP)              :: tr2
  REAL(DP)              :: alpha_mix(maxter)
  INTEGER               :: iudvwfc, lrdvwfc
  !
  LOGICAL, ALLOCATABLE :: lgamma_iq(:)
  !! if TRUE this q is gamma.
  !
  LOGICAL :: io_sp, io_real_space
  !
  !Giovanni Cistaro
  !
  REAL(DP),    ALLOCATABLE :: r(:,:)! position in the real grid
                                    ! defined by q-e
  INTEGER                  :: ir_end !last index of r array

  INTEGER(DP), ALLOCATABLE    :: s_w(:,:,:,:)
  !symmetries respected by each wannier function
  REAL(DP), ALLOCATABLE    :: ft_w(:,:,:)
  !fractional translation
  INTEGER, ALLOCATABLE     :: nsym_w_k(:)
  !number of symmetries respected by each wannier function, with the constraint 
  !of the center. Will be used for reducing k points
  INTEGER, ALLOCATABLE     :: nsym_w_q(:)
  !number of symmetries respected by each wannier function, without the contraint
  !of the center.
  !will be used only for reducing q points (are in general more than nsym_w_k)
  INTEGER, ALLOCATABLE     :: nqstot_ibz(:)
  !total number of q points in irreducible bz
  REAL(DP), ALLOCATABLE    :: xq_ibz(:,:,:)
  !q coordinates (cartesian, in units of 2pi/alat)
  REAL(DP), ALLOCATABLE    :: wq_ibz(:,:)
  !weights of q points
  INTEGER, ALLOCATABLE     :: ibz2fbz(:,:)
  !index of q point in xq_ibz in array xk
  INTEGER, ALLOCATABLE     :: fbz2ibz(:, :)
  !index of q point xk in array xq_ibz
  LOGICAL                  :: irr_bz
  INTEGER                  :: nsym_old
  !to keep track of number of symmetries of crystal, without restrictions
  LOGICAL                  :: setup_pw
!  !variable to decide wether we want to move the wannier centers in the origin when checking symmetries
  LOGICAL                  :: use_wct
  !variable to decide wether we want to verify if symmetries move a Wannier in an other unitcell in the supercell
  LOGICAL, ALLOCATABLE     :: sym_only_for_q(:,:)
END MODULE control_kcw

