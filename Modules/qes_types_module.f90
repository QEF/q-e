!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE qes_types_module
  !
  ! Auto-generated code: don't edit this file
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes-1.0
  !
  USE kinds, only: DP
  !
  IMPLICIT NONE
  !
  TYPE :: xml_format_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: NAME
    LOGICAL :: NAME_ispresent = .FALSE.
    CHARACTER(len=256) :: VERSION
    LOGICAL :: VERSION_ispresent = .FALSE.
    !
    CHARACTER(len=256) :: xml_format
    !
  END TYPE xml_format_type
  !
  TYPE :: creator_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: NAME
    LOGICAL :: NAME_ispresent = .FALSE.
    CHARACTER(len=256) :: VERSION
    LOGICAL :: VERSION_ispresent = .FALSE.
    !
    CHARACTER(len=256) :: creator
    !
  END TYPE creator_type
  !
  TYPE :: created_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: DATE
    LOGICAL :: DATE_ispresent = .FALSE.
    CHARACTER(len=256) :: TIME
    LOGICAL :: TIME_ispresent = .FALSE.
    !
    CHARACTER(len=256) :: created
    !
  END TYPE created_type
  !
  TYPE :: atom_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: name
    LOGICAL :: name_ispresent = .FALSE.
    CHARACTER(len=256) :: position
    LOGICAL :: position_ispresent = .FALSE.
    INTEGER :: index
    LOGICAL :: index_ispresent = .FALSE.
    !
    REAL(DP), DIMENSION(3) :: atom
    !
  END TYPE atom_type
  !
  TYPE :: qpoint_grid_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: nqx1
    LOGICAL :: nqx1_ispresent = .FALSE.
    INTEGER :: nqx2
    LOGICAL :: nqx2_ispresent = .FALSE.
    INTEGER :: nqx3
    LOGICAL :: nqx3_ispresent = .FALSE.
    !
    CHARACTER(len=256) :: qpoint_grid
    !
  END TYPE qpoint_grid_type
  !
  TYPE :: HubbardCommon_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: specie
    LOGICAL :: specie_ispresent = .FALSE.
    CHARACTER(len=256) :: label
    LOGICAL :: label_ispresent = .FALSE.
    !
    REAL(DP) :: HubbardCommon
    !
  END TYPE HubbardCommon_type
  !
  TYPE :: HubbardInterSpecieV_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: specie1
    INTEGER :: index1
    CHARACTER(len=256) :: label1
    LOGICAL :: label1_ispresent = .FALSE.
    CHARACTER(len=256) :: specie2
    INTEGER :: index2
    CHARACTER(len=256) :: label2
    LOGICAL :: label2_ispresent = .FALSE.
    !
    REAL(DP) :: HubbardInterSpecieV
    !
  END TYPE HubbardInterSpecieV_type
  !
  TYPE :: SiteMoment_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: species
    LOGICAL :: species_ispresent = .FALSE.
    INTEGER :: atom
    LOGICAL :: atom_ispresent = .FALSE.
    REAL(DP) :: charge
    LOGICAL :: charge_ispresent = .FALSE.
    !
    REAL(DP) :: SiteMoment
    !
  END TYPE SiteMoment_type
  !
  TYPE :: HubbardJ_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: specie
    LOGICAL :: specie_ispresent = .FALSE.
    CHARACTER(len=256) :: label
    LOGICAL :: label_ispresent = .FALSE.
    !
    REAL(DP), DIMENSION(3) :: HubbardJ
    !
  END TYPE HubbardJ_type
  !
  TYPE :: vector_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: size
    !
    REAL(DP), DIMENSION(:), ALLOCATABLE :: vector
    !
  END TYPE vector_type
  !
  TYPE :: HubbardM_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: size
    CHARACTER(len=256) :: specie
    LOGICAL :: specie_ispresent = .FALSE.
    CHARACTER(len=256) :: label
    LOGICAL :: label_ispresent = .FALSE.
    INTEGER :: spin
    LOGICAL :: spin_ispresent = .FALSE.
    REAL(DP) :: jjj
    LOGICAL :: jjj_ispresent = .FALSE.
    !
    REAL(DP), DIMENSION(:), ALLOCATABLE :: HubbardM
    !
  END TYPE HubbardM_type
  !
  TYPE :: ChannelOcc_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: specie
    LOGICAL :: specie_ispresent = .FALSE.
    CHARACTER(len=256) :: label
    LOGICAL :: label_ispresent = .FALSE.
    INTEGER :: index
    !
    REAL(DP) :: ChannelOcc
    !
  END TYPE ChannelOcc_type
  !
  TYPE :: SitMag_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: species
    LOGICAL :: species_ispresent = .FALSE.
    INTEGER :: atom
    LOGICAL :: atom_ispresent = .FALSE.
    REAL(DP) :: charge
    LOGICAL :: charge_ispresent = .FALSE.
    !
    REAL(DP), DIMENSION(3) :: SitMag
    !
  END TYPE SitMag_type
  !
  TYPE :: starting_ns_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: size
    CHARACTER(len=256) :: specie
    LOGICAL :: specie_ispresent = .FALSE.
    CHARACTER(len=256) :: label
    LOGICAL :: label_ispresent = .FALSE.
    INTEGER :: spin
    LOGICAL :: spin_ispresent = .FALSE.
    !
    REAL(DP), DIMENSION(:), ALLOCATABLE :: starting_ns
    !
  END TYPE starting_ns_type
  !
  TYPE :: integerVector_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: size
    !
    INTEGER, DIMENSION(:), ALLOCATABLE :: integerVector
    !
  END TYPE integerVector_type
  !
  TYPE :: orderUm_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: size
    CHARACTER(len=256) :: specie
    LOGICAL :: specie_ispresent = .FALSE.
    CHARACTER(len=256) :: label
    LOGICAL :: label_ispresent = .FALSE.
    INTEGER :: spin
    LOGICAL :: spin_ispresent = .FALSE.
    INTEGER :: atomidx
    LOGICAL :: atomidx_ispresent = .FALSE.
    !
    INTEGER, DIMENSION(:), ALLOCATABLE :: orderUm
    !
  END TYPE orderUm_type
  !
  TYPE :: matrix_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: rank
    INTEGER, DIMENSION(:), ALLOCATABLE :: dims
    CHARACTER(len=256) :: order
    LOGICAL :: order_ispresent = .FALSE.
    !
    REAL(DP), DIMENSION(:), ALLOCATABLE :: matrix
    !
  END TYPE matrix_type
  !
  TYPE :: Hubbard_ns_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: rank
    INTEGER, DIMENSION(:), ALLOCATABLE :: dims
    CHARACTER(len=256) :: order
    LOGICAL :: order_ispresent = .FALSE.
    CHARACTER(len=256) :: specie
    LOGICAL :: specie_ispresent = .FALSE.
    CHARACTER(len=256) :: label
    LOGICAL :: label_ispresent = .FALSE.
    INTEGER :: spin
    LOGICAL :: spin_ispresent = .FALSE.
    INTEGER :: index
    LOGICAL :: index_ispresent = .FALSE.
    !
    REAL(DP), DIMENSION(:), ALLOCATABLE :: Hubbard_ns
    !
  END TYPE Hubbard_ns_type
  !
  TYPE :: smearing_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP) :: degauss
    LOGICAL :: degauss_ispresent = .FALSE.
    !
    CHARACTER(len=256) :: smearing
    !
  END TYPE smearing_type
  !
  TYPE :: occupations_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: spin
    LOGICAL :: spin_ispresent = .FALSE.
    !
    CHARACTER(len=256) :: occupations
    !
  END TYPE occupations_type
  !
  TYPE :: basisSetItem_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: nr1
    LOGICAL :: nr1_ispresent = .FALSE.
    INTEGER :: nr2
    LOGICAL :: nr2_ispresent = .FALSE.
    INTEGER :: nr3
    LOGICAL :: nr3_ispresent = .FALSE.
    !
    CHARACTER(len=256) :: basisSetItem
    !
  END TYPE basisSetItem_type
  !
  TYPE :: monkhorst_pack_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: nk1
    LOGICAL :: nk1_ispresent = .FALSE.
    INTEGER :: nk2
    LOGICAL :: nk2_ispresent = .FALSE.
    INTEGER :: nk3
    LOGICAL :: nk3_ispresent = .FALSE.
    INTEGER :: k1
    LOGICAL :: k1_ispresent = .FALSE.
    INTEGER :: k2
    LOGICAL :: k2_ispresent = .FALSE.
    INTEGER :: k3
    LOGICAL :: k3_ispresent = .FALSE.
    !
    CHARACTER(len=256) :: monkhorst_pack
    !
  END TYPE monkhorst_pack_type
  !
  TYPE :: k_point_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP) :: weight
    LOGICAL :: weight_ispresent = .FALSE.
    CHARACTER(len=256) :: label
    LOGICAL :: label_ispresent = .FALSE.
    !
    REAL(DP), DIMENSION(3) :: k_point
    !
  END TYPE k_point_type
  !
  TYPE :: inputOccupations_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: size
    INTEGER :: ispin
    LOGICAL :: ispin_ispresent = .FALSE.
    REAL(DP) :: spin_factor
    LOGICAL :: spin_factor_ispresent = .FALSE.
    !
    REAL(DP), DIMENSION(:), ALLOCATABLE :: inputOccupations
    !
  END TYPE inputOccupations_type
  !
  TYPE :: sawtoothEnergy_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP) :: eamp
    LOGICAL :: eamp_ispresent = .FALSE.
    REAL(DP) :: eopreg
    LOGICAL :: eopreg_ispresent = .FALSE.
    REAL(DP) :: emaxpos
    LOGICAL :: emaxpos_ispresent = .FALSE.
    INTEGER :: edir
    LOGICAL :: edir_ispresent = .FALSE.
    !
    REAL(DP) :: sawtoothEnergy
    !
  END TYPE sawtoothEnergy_type
  !
  TYPE :: phase_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP) :: ionic
    LOGICAL :: ionic_ispresent = .FALSE.
    REAL(DP) :: electronic
    LOGICAL :: electronic_ispresent = .FALSE.
    CHARACTER(len=256) :: modulus
    LOGICAL :: modulus_ispresent = .FALSE.
    !
    REAL(DP) :: phase
    !
  END TYPE phase_type
  !
  TYPE :: equivalent_atoms_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: size
    INTEGER :: nat
    LOGICAL :: nat_ispresent = .FALSE.
    !
    INTEGER, DIMENSION(:), ALLOCATABLE :: equivalent_atoms
    !
  END TYPE equivalent_atoms_type
  !
  TYPE :: info_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: name
    LOGICAL :: name_ispresent = .FALSE.
    CHARACTER(len=256) :: class
    LOGICAL :: class_ispresent = .FALSE.
    LOGICAL :: time_reversal
    LOGICAL :: time_reversal_ispresent = .FALSE.
    !
    CHARACTER(len=256) :: info
    !
  END TYPE info_type
  !
  TYPE :: closed_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: DATE
    LOGICAL :: DATE_ispresent = .FALSE.
    CHARACTER(len=256) :: TIME
    LOGICAL :: TIME_ispresent = .FALSE.
    !
    CHARACTER(len=256) :: closed
    !
  END TYPE closed_type
  !
  TYPE :: cpnumstep_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: ITERATION
    LOGICAL :: ITERATION_ispresent = .FALSE.
    !
    CHARACTER(len=256) :: cpnumstep
    !
  END TYPE cpnumstep_type
  !
  TYPE :: integerMatrix_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: rank
    INTEGER, DIMENSION(:), ALLOCATABLE :: dims
    CHARACTER(len=256) :: order
    LOGICAL :: order_ispresent = .FALSE.
    !
    INTEGER, DIMENSION(:), ALLOCATABLE :: integerMatrix
    !
  END TYPE integerMatrix_type
  !
  TYPE :: scalarQuantity_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: Units
    LOGICAL :: Units_ispresent = .FALSE.
    !
    REAL(DP) :: scalarQuantity
    !
  END TYPE scalarQuantity_type
  !
  TYPE :: general_info_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(xml_format_type) :: xml_format
    TYPE(creator_type) :: creator
    TYPE(created_type) :: created
    CHARACTER(len=256) :: job
    !
  END TYPE general_info_type
  !
  TYPE :: parallel_info_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: nprocs
    INTEGER :: nthreads
    INTEGER :: ntasks
    INTEGER :: nbgrp
    INTEGER :: npool
    INTEGER :: ndiag
    !
  END TYPE parallel_info_type
  !
  TYPE :: clock_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: label
    INTEGER :: calls
    LOGICAL :: calls_ispresent = .FALSE.
    REAL(DP) :: cpu
    REAL(DP) :: wall
    !
  END TYPE clock_type
  !
  TYPE :: control_variables_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: title
    CHARACTER(len=256) :: calculation
    CHARACTER(len=256) :: restart_mode
    CHARACTER(len=256) :: prefix
    CHARACTER(len=256) :: pseudo_dir
    CHARACTER(len=256) :: outdir
    LOGICAL :: stress
    LOGICAL :: forces
    LOGICAL :: wf_collect
    CHARACTER(len=256) :: disk_io
    INTEGER :: max_seconds
    LOGICAL  :: nstep_ispresent = .FALSE.
    INTEGER :: nstep
    REAL(DP) :: etot_conv_thr
    REAL(DP) :: forc_conv_thr
    REAL(DP) :: press_conv_thr
    CHARACTER(len=256) :: verbosity
    INTEGER :: print_every
    LOGICAL :: fcp
    LOGICAL :: rism
    !
  END TYPE control_variables_type
  !
  TYPE :: species_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: name
    LOGICAL :: name_ispresent = .FALSE.
    LOGICAL  :: mass_ispresent = .FALSE.
    REAL(DP) :: mass
    CHARACTER(len=256) :: pseudo_file
    LOGICAL  :: starting_magnetization_ispresent = .FALSE.
    REAL(DP) :: starting_magnetization
    LOGICAL  :: spin_teta_ispresent = .FALSE.
    REAL(DP) :: spin_teta
    LOGICAL  :: spin_phi_ispresent = .FALSE.
    REAL(DP) :: spin_phi
    !
  END TYPE species_type
  !
  TYPE :: atomic_positions_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(atom_type), DIMENSION(:), ALLOCATABLE :: atom
    INTEGER   :: ndim_atom
    !
  END TYPE atomic_positions_type
  !
  TYPE :: wyckoff_positions_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: space_group
    LOGICAL :: space_group_ispresent = .FALSE.
    CHARACTER(len=256) :: more_options
    LOGICAL :: more_options_ispresent = .FALSE.
    TYPE(atom_type), DIMENSION(:), ALLOCATABLE :: atom
    INTEGER   :: ndim_atom
    !
  END TYPE wyckoff_positions_type
  !
  TYPE :: cell_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP), DIMENSION(3) :: a1
    REAL(DP), DIMENSION(3) :: a2
    REAL(DP), DIMENSION(3) :: a3
    !
  END TYPE cell_type
  !
  TYPE :: hybrid_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: qpoint_grid_ispresent = .FALSE.
    TYPE(qpoint_grid_type) :: qpoint_grid
    LOGICAL  :: ecutfock_ispresent = .FALSE.
    REAL(DP) :: ecutfock
    LOGICAL  :: exx_fraction_ispresent = .FALSE.
    REAL(DP) :: exx_fraction
    LOGICAL  :: screening_parameter_ispresent = .FALSE.
    REAL(DP) :: screening_parameter
    LOGICAL  :: exxdiv_treatment_ispresent = .FALSE.
    CHARACTER(len=256) :: exxdiv_treatment
    LOGICAL  :: x_gamma_extrapolation_ispresent = .FALSE.
    LOGICAL :: x_gamma_extrapolation
    LOGICAL  :: ecutvcut_ispresent = .FALSE.
    REAL(DP) :: ecutvcut
    LOGICAL  :: localization_threshold_ispresent = .FALSE.
    REAL(DP) :: localization_threshold
    !
  END TYPE hybrid_type
  !
  TYPE :: HubbardOcc_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: channels
    CHARACTER(len=256) :: specie
    TYPE(ChannelOcc_type), DIMENSION(:), ALLOCATABLE :: channel_occ
    INTEGER   :: ndim_channel_occ
    !
  END TYPE HubbardOcc_type
  !
  TYPE :: HubbardBack_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: background
    CHARACTER(len=256) :: label
    LOGICAL :: label_ispresent = .FALSE.
    CHARACTER(len=256) :: species
    LOGICAL :: species_ispresent = .FALSE.
    REAL(DP) :: Hubbard_U2
    INTEGER :: n2_number
    INTEGER :: l2_number
    LOGICAL  :: n3_number_ispresent = .FALSE.
    INTEGER :: n3_number
    LOGICAL  :: l3_number_ispresent = .FALSE.
    INTEGER :: l3_number
    !
  END TYPE HubbardBack_type
  !
  TYPE :: vdW_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: vdw_corr_ispresent = .FALSE.
    CHARACTER(len=256) :: vdw_corr
    LOGICAL  :: dftd3_version_ispresent = .FALSE.
    INTEGER :: dftd3_version
    LOGICAL  :: dftd3_threebody_ispresent = .FALSE.
    LOGICAL :: dftd3_threebody
    LOGICAL  :: non_local_term_ispresent = .FALSE.
    CHARACTER(len=256) :: non_local_term
    LOGICAL  :: functional_ispresent = .FALSE.
    CHARACTER(len=256) :: functional
    LOGICAL  :: total_energy_term_ispresent = .FALSE.
    REAL(DP) :: total_energy_term
    LOGICAL  :: london_s6_ispresent = .FALSE.
    REAL(DP) :: london_s6
    LOGICAL  :: ts_vdw_econv_thr_ispresent = .FALSE.
    REAL(DP) :: ts_vdw_econv_thr
    LOGICAL  :: ts_vdw_isolated_ispresent = .FALSE.
    LOGICAL :: ts_vdw_isolated
    LOGICAL  :: london_rcut_ispresent = .FALSE.
    REAL(DP) :: london_rcut
    LOGICAL  :: xdm_a1_ispresent = .FALSE.
    REAL(DP) :: xdm_a1
    LOGICAL  :: xdm_a2_ispresent = .FALSE.
    REAL(DP) :: xdm_a2
    LOGICAL  :: london_c6_ispresent = .FALSE.
    TYPE(HubbardCommon_type), DIMENSION(:), ALLOCATABLE :: london_c6
    INTEGER   :: ndim_london_c6
    !
  END TYPE vdW_type
  !
  TYPE :: spin_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL :: lsda
    LOGICAL :: noncolin
    LOGICAL :: spinorbit
    !
  END TYPE spin_type
  !
  TYPE :: bands_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: nbnd_ispresent = .FALSE.
    INTEGER :: nbnd
    LOGICAL  :: smearing_ispresent = .FALSE.
    TYPE(smearing_type) :: smearing
    LOGICAL  :: tot_charge_ispresent = .FALSE.
    REAL(DP) :: tot_charge
    LOGICAL  :: tot_magnetization_ispresent = .FALSE.
    REAL(DP) :: tot_magnetization
    TYPE(occupations_type) :: occupations
    LOGICAL  :: inputOccupations_ispresent = .FALSE.
    TYPE(inputOccupations_type), DIMENSION(:), ALLOCATABLE :: inputOccupations
    INTEGER   :: ndim_inputOccupations
    !
  END TYPE bands_type
  !
  TYPE :: basis_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: gamma_only_ispresent = .FALSE.
    LOGICAL :: gamma_only
    REAL(DP) :: ecutwfc
    LOGICAL  :: ecutrho_ispresent = .FALSE.
    REAL(DP) :: ecutrho
    LOGICAL  :: fft_grid_ispresent = .FALSE.
    TYPE(basisSetItem_type) :: fft_grid
    LOGICAL  :: fft_smooth_ispresent = .FALSE.
    TYPE(basisSetItem_type) :: fft_smooth
    LOGICAL  :: fft_box_ispresent = .FALSE.
    TYPE(basisSetItem_type) :: fft_box
    !
  END TYPE basis_type
  !
  TYPE :: reciprocal_lattice_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP), DIMENSION(3) :: b1
    REAL(DP), DIMENSION(3) :: b2
    REAL(DP), DIMENSION(3) :: b3
    !
  END TYPE reciprocal_lattice_type
  !
  TYPE :: electron_control_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: diagonalization
    CHARACTER(len=256) :: mixing_mode
    REAL(DP) :: mixing_beta
    REAL(DP) :: conv_thr
    INTEGER :: mixing_ndim
    INTEGER :: max_nstep
    LOGICAL  :: exx_nstep_ispresent = .FALSE.
    INTEGER :: exx_nstep
    LOGICAL  :: real_space_q_ispresent = .FALSE.
    LOGICAL :: real_space_q
    LOGICAL  :: real_space_beta_ispresent = .FALSE.
    LOGICAL :: real_space_beta
    LOGICAL :: tq_smoothing
    LOGICAL :: tbeta_smoothing
    REAL(DP) :: diago_thr_init
    LOGICAL :: diago_full_acc
    LOGICAL  :: diago_cg_maxiter_ispresent = .FALSE.
    INTEGER :: diago_cg_maxiter
    LOGICAL  :: diago_ppcg_maxiter_ispresent = .FALSE.
    INTEGER :: diago_ppcg_maxiter
    LOGICAL  :: diago_david_ndim_ispresent = .FALSE.
    INTEGER :: diago_david_ndim
    LOGICAL  :: diago_rmm_ndim_ispresent = .FALSE.
    INTEGER :: diago_rmm_ndim
    LOGICAL  :: diago_gs_nblock_ispresent = .FALSE.
    INTEGER :: diago_gs_nblock
    LOGICAL  :: diago_rmm_conv_ispresent = .FALSE.
    LOGICAL :: diago_rmm_conv
    !
  END TYPE electron_control_type
  !
  TYPE :: fcp_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: fcp_mu_ispresent = .FALSE.
    REAL(DP) :: fcp_mu
    LOGICAL  :: fcp_dynamics_ispresent = .FALSE.
    CHARACTER(len=256) :: fcp_dynamics
    LOGICAL  :: fcp_conv_thr_ispresent = .FALSE.
    REAL(DP) :: fcp_conv_thr
    LOGICAL  :: fcp_ndiis_ispresent = .FALSE.
    INTEGER :: fcp_ndiis
    LOGICAL  :: fcp_rdiis_ispresent = .FALSE.
    REAL(DP) :: fcp_rdiis
    LOGICAL  :: fcp_mass_ispresent = .FALSE.
    REAL(DP) :: fcp_mass
    LOGICAL  :: fcp_velocity_ispresent = .FALSE.
    REAL(DP) :: fcp_velocity
    LOGICAL  :: fcp_temperature_ispresent = .FALSE.
    CHARACTER(len=256) :: fcp_temperature
    LOGICAL  :: fcp_tempw_ispresent = .FALSE.
    REAL(DP) :: fcp_tempw
    LOGICAL  :: fcp_tolp_ispresent = .FALSE.
    REAL(DP) :: fcp_tolp
    LOGICAL  :: fcp_delta_t_ispresent = .FALSE.
    REAL(DP) :: fcp_delta_t
    LOGICAL  :: fcp_nraise_ispresent = .FALSE.
    INTEGER :: fcp_nraise
    LOGICAL  :: freeze_all_atoms_ispresent = .FALSE.
    LOGICAL :: freeze_all_atoms
    !
  END TYPE fcp_type
  !
  TYPE :: solute_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: solute_lj
    REAL(DP) :: epsilon
    REAL(DP) :: sigma
    !
  END TYPE solute_type
  !
  TYPE :: solvent_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: label
    CHARACTER(len=256) :: molec_file
    REAL(DP) :: density1
    LOGICAL  :: density2_ispresent = .FALSE.
    REAL(DP) :: density2
    LOGICAL  :: unit_ispresent = .FALSE.
    CHARACTER(len=256) :: unit
    !
  END TYPE solvent_type
  !
  TYPE :: k_points_IBZ_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: monkhorst_pack_ispresent = .FALSE.
    TYPE(monkhorst_pack_type) :: monkhorst_pack
    LOGICAL  :: nk_ispresent = .FALSE.
    INTEGER :: nk
    LOGICAL  :: k_point_ispresent = .FALSE.
    TYPE(k_point_type), DIMENSION(:), ALLOCATABLE :: k_point
    INTEGER   :: ndim_k_point
    !
  END TYPE k_points_IBZ_type
  !
  TYPE :: bfgs_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: ndim
    REAL(DP) :: trust_radius_min
    REAL(DP) :: trust_radius_max
    REAL(DP) :: trust_radius_init
    REAL(DP) :: w1
    REAL(DP) :: w2
    !
  END TYPE bfgs_type
  !
  TYPE :: md_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: pot_extrapolation
    CHARACTER(len=256) :: wfc_extrapolation
    CHARACTER(len=256) :: ion_temperature
    REAL(DP) :: timestep
    REAL(DP) :: tempw
    REAL(DP) :: tolp
    REAL(DP) :: deltaT
    INTEGER :: nraise
    !
  END TYPE md_type
  !
  TYPE :: cell_control_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: cell_dynamics
    REAL(DP) :: pressure
    LOGICAL  :: wmass_ispresent = .FALSE.
    REAL(DP) :: wmass
    LOGICAL  :: cell_factor_ispresent = .FALSE.
    REAL(DP) :: cell_factor
    LOGICAL  :: cell_do_free_ispresent = .FALSE.
    CHARACTER(len=256) :: cell_do_free
    LOGICAL  :: fix_volume_ispresent = .FALSE.
    LOGICAL :: fix_volume
    LOGICAL  :: fix_area_ispresent = .FALSE.
    LOGICAL :: fix_area
    LOGICAL  :: isotropic_ispresent = .FALSE.
    LOGICAL :: isotropic
    LOGICAL  :: free_cell_ispresent = .FALSE.
    TYPE(integerMatrix_type) :: free_cell
    !
  END TYPE cell_control_type
  !
  TYPE :: symmetry_flags_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL :: nosym
    LOGICAL :: nosym_evc
    LOGICAL :: noinv
    LOGICAL :: no_t_rev
    LOGICAL :: force_symmorphic
    LOGICAL :: use_all_frac
    !
  END TYPE symmetry_flags_type
  !
  TYPE :: esm_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: bc
    LOGICAL  :: nfit_ispresent = .FALSE.
    INTEGER :: nfit
    LOGICAL  :: w_ispresent = .FALSE.
    REAL(DP) :: w
    LOGICAL  :: efield_ispresent = .FALSE.
    REAL(DP) :: efield
    LOGICAL  :: a_ispresent = .FALSE.
    REAL(DP) :: a
    LOGICAL  :: zb_ispresent = .FALSE.
    REAL(DP) :: zb
    LOGICAL  :: debug_ispresent = .FALSE.
    LOGICAL :: debug
    LOGICAL  :: debug_gpmax_ispresent = .FALSE.
    INTEGER :: debug_gpmax
    !
  END TYPE esm_type
  !
  TYPE :: gcscf_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: ignore_mun_ispresent = .FALSE.
    LOGICAL :: ignore_mun
    LOGICAL  :: mu_ispresent = .FALSE.
    REAL(DP) :: mu
    LOGICAL  :: conv_thr_ispresent = .FALSE.
    REAL(DP) :: conv_thr
    LOGICAL  :: gk_ispresent = .FALSE.
    REAL(DP) :: gk
    LOGICAL  :: gh_ispresent = .FALSE.
    REAL(DP) :: gh
    LOGICAL  :: beta_ispresent = .FALSE.
    REAL(DP) :: beta
    !
  END TYPE gcscf_type
  !
  TYPE :: ekin_functional_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP) :: ecfixed
    REAL(DP) :: qcutz
    REAL(DP) :: q2sigma
    !
  END TYPE ekin_functional_type
  !
  TYPE :: spin_constraints_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: spin_constraints
    REAL(DP) :: lagrange_multiplier
    LOGICAL  :: target_magnetization_ispresent = .FALSE.
    REAL(DP), DIMENSION(3) :: target_magnetization
    !
  END TYPE spin_constraints_type
  !
  TYPE :: gate_settings_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL :: use_gate
    LOGICAL  :: zgate_ispresent = .FALSE.
    REAL(DP) :: zgate
    LOGICAL  :: relaxz_ispresent = .FALSE.
    LOGICAL :: relaxz
    LOGICAL  :: block_ispresent = .FALSE.
    LOGICAL :: block
    LOGICAL  :: block_1_ispresent = .FALSE.
    REAL(DP) :: block_1
    LOGICAL  :: block_2_ispresent = .FALSE.
    REAL(DP) :: block_2
    LOGICAL  :: block_height_ispresent = .FALSE.
    REAL(DP) :: block_height
    !
  END TYPE gate_settings_type
  !
  TYPE :: atomic_constraint_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP), DIMENSION(4) :: constr_parms
    CHARACTER(len=256) :: constr_type
    LOGICAL  :: constr_target_ispresent = .FALSE.
    REAL(DP) :: constr_target
    !
  END TYPE atomic_constraint_type
  !
  TYPE :: dipoleOutput_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: idir
    TYPE(scalarQuantity_type) :: dipole
    TYPE(scalarQuantity_type) :: ion_dipole
    TYPE(scalarQuantity_type) :: elec_dipole
    TYPE(scalarQuantity_type) :: dipoleField
    TYPE(scalarQuantity_type) :: potentialAmp
    TYPE(scalarQuantity_type) :: totalLength
    !
  END TYPE dipoleOutput_type
  !
  TYPE :: finiteFieldOut_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP), DIMENSION(3) :: electronicDipole
    REAL(DP), DIMENSION(3) :: ionicDipole
    !
  END TYPE finiteFieldOut_type
  !
  TYPE :: polarization_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(scalarQuantity_type) :: polarization
    REAL(DP) :: modulus
    REAL(DP), DIMENSION(3) :: direction
    !
  END TYPE polarization_type
  !
  TYPE :: ionicPolarization_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(atom_type) :: ion
    REAL(DP) :: charge
    TYPE(phase_type) :: phase
    !
  END TYPE ionicPolarization_type
  !
  TYPE :: electronicPolarization_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(k_point_type) :: firstKeyPoint
    LOGICAL  :: spin_ispresent = .FALSE.
    INTEGER :: spin
    TYPE(phase_type) :: phase
    !
  END TYPE electronicPolarization_type
  !
  TYPE :: gateInfo_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP) :: pot_prefactor
    REAL(DP) :: gate_zpos
    REAL(DP) :: gate_gate_term
    REAL(DP) :: gatefieldEnergy
    !
  END TYPE gateInfo_type
  !
  TYPE :: scf_conv_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL :: convergence_achieved
    INTEGER :: n_scf_steps
    REAL(DP) :: scf_error
    !
  END TYPE scf_conv_type
  !
  TYPE :: opt_conv_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL :: convergence_achieved
    INTEGER :: n_opt_steps
    REAL(DP) :: grad_norm
    !
  END TYPE opt_conv_type
  !
  TYPE :: algorithmic_info_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL :: real_space_q
    LOGICAL  :: real_space_beta_ispresent = .FALSE.
    LOGICAL :: real_space_beta
    LOGICAL :: uspp
    LOGICAL :: paw
    !
  END TYPE algorithmic_info_type
  !
  TYPE :: symmetry_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(info_type) :: info
    TYPE(matrix_type) :: rotation
    LOGICAL  :: fractional_translation_ispresent = .FALSE.
    REAL(DP), DIMENSION(3) :: fractional_translation
    LOGICAL  :: equivalent_atoms_ispresent = .FALSE.
    TYPE(equivalent_atoms_type) :: equivalent_atoms
    !
  END TYPE symmetry_type
  !
  TYPE :: total_energy_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP) :: etot
    LOGICAL  :: eband_ispresent = .FALSE.
    REAL(DP) :: eband
    LOGICAL  :: ehart_ispresent = .FALSE.
    REAL(DP) :: ehart
    LOGICAL  :: vtxc_ispresent = .FALSE.
    REAL(DP) :: vtxc
    LOGICAL  :: etxc_ispresent = .FALSE.
    REAL(DP) :: etxc
    LOGICAL  :: ewald_ispresent = .FALSE.
    REAL(DP) :: ewald
    LOGICAL  :: demet_ispresent = .FALSE.
    REAL(DP) :: demet
    LOGICAL  :: efieldcorr_ispresent = .FALSE.
    REAL(DP) :: efieldcorr
    LOGICAL  :: potentiostat_contr_ispresent = .FALSE.
    REAL(DP) :: potentiostat_contr
    LOGICAL  :: gatefield_contr_ispresent = .FALSE.
    REAL(DP) :: gatefield_contr
    LOGICAL  :: vdW_term_ispresent = .FALSE.
    REAL(DP) :: vdW_term
    LOGICAL  :: esol_ispresent = .FALSE.
    REAL(DP) :: esol
    LOGICAL  :: levelshift_contr_ispresent = .FALSE.
    REAL(DP) :: levelshift_contr
    !
  END TYPE total_energy_type
  !
  TYPE :: ks_energies_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(k_point_type) :: k_point
    INTEGER :: npw
    TYPE(vector_type) :: eigenvalues
    TYPE(vector_type) :: occupations
    !
  END TYPE ks_energies_type
  !
  TYPE :: cpstatus_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(cpnumstep_type) :: STEP
    TYPE(scalarQuantity_type) :: TIME
    CHARACTER(len=256) :: TITLE
    TYPE(scalarQuantity_type) :: KINETIC_ENERGY
    TYPE(scalarQuantity_type) :: HARTREE_ENERGY
    TYPE(scalarQuantity_type) :: EWALD_TERM
    TYPE(scalarQuantity_type) :: GAUSS_SELFINT
    TYPE(scalarQuantity_type) :: LPSP_ENERGY
    TYPE(scalarQuantity_type) :: NLPSP_ENERGY
    TYPE(scalarQuantity_type) :: EXC_ENERGY
    TYPE(scalarQuantity_type) :: AVERAGE_POT
    TYPE(scalarQuantity_type) :: ENTHALPY
    !
  END TYPE cpstatus_type
  !
  TYPE :: cp_ionPos_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP), DIMENSION(:), ALLOCATABLE :: stau
    REAL(DP), DIMENSION(:), ALLOCATABLE :: svel
    LOGICAL  :: taui_ispresent = .FALSE.
    REAL(DP), DIMENSION(:), ALLOCATABLE :: taui
    LOGICAL  :: cdmi_ispresent = .FALSE.
    REAL(DP), DIMENSION(3) :: cdmi
    LOGICAL  :: force_ispresent = .FALSE.
    REAL(DP), DIMENSION(:), ALLOCATABLE :: force
    !
  END TYPE cp_ionPos_type
  !
  TYPE :: cp_ionsNose_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: nhpcl
    INTEGER :: nhpdim
    REAL(DP) :: xnhp
    LOGICAL  :: vnhp_ispresent = .FALSE.
    REAL(DP) :: vnhp
    !
  END TYPE cp_ionsNose_type
  !
  TYPE :: cp_elecNose_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP) :: xnhe
    LOGICAL  :: vnhe_ispresent = .FALSE.
    REAL(DP) :: vnhe
    !
  END TYPE cp_elecNose_type
  !
  TYPE :: cp_cell_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP), DIMENSION(:), ALLOCATABLE :: ht
    LOGICAL  :: htvel_ispresent = .FALSE.
    REAL(DP), DIMENSION(:), ALLOCATABLE :: htvel
    LOGICAL  :: gvel_ispresent = .FALSE.
    REAL(DP), DIMENSION(:), ALLOCATABLE :: gvel
    !
  END TYPE cp_cell_type
  !
  TYPE :: cp_cellNose_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    REAL(DP), DIMENSION(:), ALLOCATABLE :: xnhh
    LOGICAL  :: vnhh_ispresent = .FALSE.
    REAL(DP), DIMENSION(:), ALLOCATABLE :: vnhh
    !
  END TYPE cp_cellNose_type
  !
  TYPE :: scalmags_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: nat
    LOGICAL :: nat_ispresent = .FALSE.
    TYPE(SiteMoment_type), DIMENSION(:), ALLOCATABLE :: SiteMagnetization
    INTEGER   :: ndim_SiteMagnetization
    !
  END TYPE scalmags_type
  !
  TYPE :: d3mags_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: nat
    LOGICAL :: nat_ispresent = .FALSE.
    TYPE(SitMag_type), DIMENSION(:), ALLOCATABLE :: SiteMagnetization
    INTEGER   :: ndim_SiteMagnetization
    !
  END TYPE d3mags_type
  !
  TYPE :: rismlaue_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: both_hands_ispresent = .FALSE.
    LOGICAL :: both_hands
    LOGICAL  :: nfit_ispresent = .FALSE.
    INTEGER :: nfit
    LOGICAL  :: pot_ref_ispresent = .FALSE.
    INTEGER :: pot_ref
    LOGICAL  :: charge_ispresent = .FALSE.
    REAL(DP) :: charge
    LOGICAL  :: right_start_ispresent = .FALSE.
    REAL(DP) :: right_start
    LOGICAL  :: right_expand_ispresent = .FALSE.
    REAL(DP) :: right_expand
    LOGICAL  :: right_buffer_ispresent = .FALSE.
    REAL(DP) :: right_buffer
    LOGICAL  :: right_buffer_u_ispresent = .FALSE.
    REAL(DP) :: right_buffer_u
    LOGICAL  :: right_buffer_v_ispresent = .FALSE.
    REAL(DP) :: right_buffer_v
    LOGICAL  :: left_start_ispresent = .FALSE.
    REAL(DP) :: left_start
    LOGICAL  :: left_expand_ispresent = .FALSE.
    REAL(DP) :: left_expand
    LOGICAL  :: left_buffer_ispresent = .FALSE.
    REAL(DP) :: left_buffer
    LOGICAL  :: left_buffer_u_ispresent = .FALSE.
    REAL(DP) :: left_buffer_u
    LOGICAL  :: left_buffer_v_ispresent = .FALSE.
    REAL(DP) :: left_buffer_v
    !
  END TYPE rismlaue_type
  !
  TYPE :: two_chem_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL :: twochem
    INTEGER :: nbnd_cond
    REAL(DP) :: degauss_cond
    REAL(DP) :: nelec_cond
    LOGICAL  :: ef_cond_ispresent = .FALSE.
    REAL(DP) :: ef_cond
    !
  END TYPE two_chem_type
  !
  TYPE :: timing_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(clock_type) :: total
    LOGICAL  :: partial_ispresent = .FALSE.
    TYPE(clock_type), DIMENSION(:), ALLOCATABLE :: partial
    INTEGER   :: ndim_partial
    !
  END TYPE timing_type
  !
  TYPE :: atomic_species_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: ntyp
    LOGICAL :: ntyp_ispresent = .FALSE.
    CHARACTER(len=256) :: pseudo_dir
    LOGICAL :: pseudo_dir_ispresent = .FALSE.
    TYPE(species_type), DIMENSION(:), ALLOCATABLE :: species
    INTEGER   :: ndim_species
    !
  END TYPE atomic_species_type
  !
  TYPE :: atomic_structure_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: nat
    LOGICAL :: nat_ispresent = .FALSE.
    INTEGER :: num_of_atomic_wfc
    LOGICAL :: num_of_atomic_wfc_ispresent = .FALSE.
    REAL(DP) :: alat
    LOGICAL :: alat_ispresent = .FALSE.
    INTEGER :: bravais_index
    LOGICAL :: bravais_index_ispresent = .FALSE.
    CHARACTER(len=256) :: alternative_axes
    LOGICAL :: alternative_axes_ispresent = .FALSE.
    LOGICAL  :: atomic_positions_ispresent = .FALSE.
    TYPE(atomic_positions_type) :: atomic_positions
    LOGICAL  :: wyckoff_positions_ispresent = .FALSE.
    TYPE(wyckoff_positions_type) :: wyckoff_positions
    LOGICAL  :: crystal_positions_ispresent = .FALSE.
    TYPE(atomic_positions_type) :: crystal_positions
    TYPE(cell_type) :: cell
    !
  END TYPE atomic_structure_type
  !
  TYPE :: dftU_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL :: new_format
    LOGICAL :: new_format_ispresent = .FALSE.
    LOGICAL  :: lda_plus_u_kind_ispresent = .FALSE.
    INTEGER :: lda_plus_u_kind
    LOGICAL  :: Hubbard_Occ_ispresent = .FALSE.
    TYPE(HubbardOcc_type), DIMENSION(:), ALLOCATABLE :: Hubbard_Occ
    INTEGER   :: ndim_Hubbard_Occ
    LOGICAL  :: Hubbard_U_ispresent = .FALSE.
    TYPE(HubbardCommon_type), DIMENSION(:), ALLOCATABLE :: Hubbard_U
    INTEGER   :: ndim_Hubbard_U
    LOGICAL  :: Hubbard_Um_ispresent = .FALSE.
    TYPE(HubbardM_type), DIMENSION(:), ALLOCATABLE :: Hubbard_Um
    INTEGER   :: ndim_Hubbard_Um
    LOGICAL  :: Hubbard_J0_ispresent = .FALSE.
    TYPE(HubbardCommon_type), DIMENSION(:), ALLOCATABLE :: Hubbard_J0
    INTEGER   :: ndim_Hubbard_J0
    LOGICAL  :: Hubbard_alpha_ispresent = .FALSE.
    TYPE(HubbardCommon_type), DIMENSION(:), ALLOCATABLE :: Hubbard_alpha
    INTEGER   :: ndim_Hubbard_alpha
    LOGICAL  :: Hubbard_beta_ispresent = .FALSE.
    TYPE(HubbardCommon_type), DIMENSION(:), ALLOCATABLE :: Hubbard_beta
    INTEGER   :: ndim_Hubbard_beta
    LOGICAL  :: Hubbard_J_ispresent = .FALSE.
    TYPE(HubbardJ_type), DIMENSION(:), ALLOCATABLE :: Hubbard_J
    INTEGER   :: ndim_Hubbard_J
    LOGICAL  :: starting_ns_ispresent = .FALSE.
    TYPE(starting_ns_type), DIMENSION(:), ALLOCATABLE :: starting_ns
    INTEGER   :: ndim_starting_ns
    LOGICAL  :: Hubbard_V_ispresent = .FALSE.
    TYPE(HubbardInterSpecieV_type), DIMENSION(:), ALLOCATABLE :: Hubbard_V
    INTEGER   :: ndim_Hubbard_V
    LOGICAL  :: Hubbard_ns_ispresent = .FALSE.
    TYPE(Hubbard_ns_type), DIMENSION(:), ALLOCATABLE :: Hubbard_ns
    INTEGER   :: ndim_Hubbard_ns
    LOGICAL  :: Hub_m_order_ispresent = .FALSE.
    TYPE(orderUm_type), DIMENSION(:), ALLOCATABLE :: Hub_m_order
    INTEGER   :: ndim_Hub_m_order
    LOGICAL  :: U_projection_type_ispresent = .FALSE.
    CHARACTER(len=256) :: U_projection_type
    LOGICAL  :: Hubbard_back_ispresent = .FALSE.
    TYPE(HubbardBack_type), DIMENSION(:), ALLOCATABLE :: Hubbard_back
    INTEGER   :: ndim_Hubbard_back
    LOGICAL  :: Hubbard_alpha_back_ispresent = .FALSE.
    TYPE(HubbardCommon_type), DIMENSION(:), ALLOCATABLE :: Hubbard_alpha_back
    INTEGER   :: ndim_Hubbard_alpha_back
    LOGICAL  :: Hubbard_ns_nc_ispresent = .FALSE.
    TYPE(Hubbard_ns_type), DIMENSION(:), ALLOCATABLE :: Hubbard_ns_nc
    INTEGER   :: ndim_Hubbard_ns_nc
    !
  END TYPE dftU_type
  !
  TYPE :: basis_set_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: gamma_only_ispresent = .FALSE.
    LOGICAL :: gamma_only
    REAL(DP) :: ecutwfc
    LOGICAL  :: ecutrho_ispresent = .FALSE.
    REAL(DP) :: ecutrho
    TYPE(basisSetItem_type) :: fft_grid
    LOGICAL  :: fft_smooth_ispresent = .FALSE.
    TYPE(basisSetItem_type) :: fft_smooth
    LOGICAL  :: fft_box_ispresent = .FALSE.
    TYPE(basisSetItem_type) :: fft_box
    INTEGER :: ngm
    LOGICAL  :: ngms_ispresent = .FALSE.
    INTEGER :: ngms
    INTEGER :: npwx
    TYPE(reciprocal_lattice_type) :: reciprocal_lattice
    !
  END TYPE basis_set_type
  !
  TYPE :: rism_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: nsolv
    TYPE(solute_type), DIMENSION(:), ALLOCATABLE :: solute
    INTEGER   :: ndim_solute
    LOGICAL  :: closure_ispresent = .FALSE.
    CHARACTER(len=256) :: closure
    LOGICAL  :: tempv_ispresent = .FALSE.
    REAL(DP) :: tempv
    LOGICAL  :: ecutsolv_ispresent = .FALSE.
    REAL(DP) :: ecutsolv
    LOGICAL  :: rmax_lj_ispresent = .FALSE.
    REAL(DP) :: rmax_lj
    LOGICAL  :: rmax1d_ispresent = .FALSE.
    REAL(DP) :: rmax1d
    LOGICAL  :: starting1d_ispresent = .FALSE.
    CHARACTER(len=256) :: starting1d
    LOGICAL  :: starting3d_ispresent = .FALSE.
    CHARACTER(len=256) :: starting3d
    LOGICAL  :: smear1d_ispresent = .FALSE.
    REAL(DP) :: smear1d
    LOGICAL  :: smear3d_ispresent = .FALSE.
    REAL(DP) :: smear3d
    LOGICAL  :: rism1d_maxstep_ispresent = .FALSE.
    INTEGER :: rism1d_maxstep
    LOGICAL  :: rism3d_maxstep_ispresent = .FALSE.
    INTEGER :: rism3d_maxstep
    LOGICAL  :: rism1d_conv_thr_ispresent = .FALSE.
    REAL(DP) :: rism1d_conv_thr
    LOGICAL  :: rism3d_conv_thr_ispresent = .FALSE.
    REAL(DP) :: rism3d_conv_thr
    LOGICAL  :: mdiis1d_size_ispresent = .FALSE.
    INTEGER :: mdiis1d_size
    LOGICAL  :: mdiis3d_size_ispresent = .FALSE.
    INTEGER :: mdiis3d_size
    LOGICAL  :: mdiis1d_step_ispresent = .FALSE.
    REAL(DP) :: mdiis1d_step
    LOGICAL  :: mdiis3d_step_ispresent = .FALSE.
    REAL(DP) :: mdiis3d_step
    LOGICAL  :: rism1d_bond_width_ispresent = .FALSE.
    REAL(DP) :: rism1d_bond_width
    LOGICAL  :: rism1d_dielectric_ispresent = .FALSE.
    REAL(DP) :: rism1d_dielectric
    LOGICAL  :: rism1d_molesize_ispresent = .FALSE.
    REAL(DP) :: rism1d_molesize
    LOGICAL  :: rism1d_nproc_ispresent = .FALSE.
    INTEGER :: rism1d_nproc
    LOGICAL  :: rism1d_nproc_switch_ispresent = .FALSE.
    INTEGER :: rism1d_nproc_switch
    LOGICAL  :: rism3d_conv_level_ispresent = .FALSE.
    REAL(DP) :: rism3d_conv_level
    LOGICAL  :: rism3d_planar_average_ispresent = .FALSE.
    LOGICAL :: rism3d_planar_average
    LOGICAL  :: laue_nfit_ispresent = .FALSE.
    INTEGER :: laue_nfit
    LOGICAL  :: laue_expand_right_ispresent = .FALSE.
    REAL(DP) :: laue_expand_right
    LOGICAL  :: laue_expand_left_ispresent = .FALSE.
    REAL(DP) :: laue_expand_left
    LOGICAL  :: laue_starting_right_ispresent = .FALSE.
    REAL(DP) :: laue_starting_right
    LOGICAL  :: laue_starting_left_ispresent = .FALSE.
    REAL(DP) :: laue_starting_left
    LOGICAL  :: laue_buffer_right_ispresent = .FALSE.
    REAL(DP) :: laue_buffer_right
    LOGICAL  :: laue_buffer_right_solu_ispresent = .FALSE.
    REAL(DP) :: laue_buffer_right_solu
    LOGICAL  :: laue_buffer_right_solv_ispresent = .FALSE.
    REAL(DP) :: laue_buffer_right_solv
    LOGICAL  :: laue_buffer_left_ispresent = .FALSE.
    REAL(DP) :: laue_buffer_left
    LOGICAL  :: laue_buffer_left_solu_ispresent = .FALSE.
    REAL(DP) :: laue_buffer_left_solu
    LOGICAL  :: laue_buffer_left_solv_ispresent = .FALSE.
    REAL(DP) :: laue_buffer_left_solv
    LOGICAL  :: laue_both_hands_ispresent = .FALSE.
    LOGICAL :: laue_both_hands
    LOGICAL  :: laue_reference_ispresent = .FALSE.
    CHARACTER(len=256) :: laue_reference
    LOGICAL  :: laue_wall_ispresent = .FALSE.
    CHARACTER(len=256) :: laue_wall
    LOGICAL  :: laue_wall_z_ispresent = .FALSE.
    REAL(DP) :: laue_wall_z
    LOGICAL  :: laue_wall_rho_ispresent = .FALSE.
    REAL(DP) :: laue_wall_rho
    LOGICAL  :: laue_wall_epsilon_ispresent = .FALSE.
    REAL(DP) :: laue_wall_epsilon
    LOGICAL  :: laue_wall_sigma_ispresent = .FALSE.
    REAL(DP) :: laue_wall_sigma
    LOGICAL  :: laue_wall_lj6_ispresent = .FALSE.
    LOGICAL :: laue_wall_lj6
    !
  END TYPE rism_type
  !
  TYPE :: ion_control_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: ion_dynamics
    LOGICAL  :: upscale_ispresent = .FALSE.
    REAL(DP) :: upscale
    LOGICAL  :: remove_rigid_rot_ispresent = .FALSE.
    LOGICAL :: remove_rigid_rot
    LOGICAL  :: refold_pos_ispresent = .FALSE.
    LOGICAL :: refold_pos
    LOGICAL  :: bfgs_ispresent = .FALSE.
    TYPE(bfgs_type) :: bfgs
    LOGICAL  :: md_ispresent = .FALSE.
    TYPE(md_type) :: md
    !
  END TYPE ion_control_type
  !
  TYPE :: boundary_conditions_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: assume_isolated
    LOGICAL  :: esm_ispresent = .FALSE.
    TYPE(esm_type) :: esm
    LOGICAL  :: gcscf_ispresent = .FALSE.
    TYPE(gcscf_type) :: gcscf
    !
  END TYPE boundary_conditions_type
  !
  TYPE :: solvents_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(solvent_type), DIMENSION(:), ALLOCATABLE :: solvent
    INTEGER   :: ndim_solvent
    !
  END TYPE solvents_type
  !
  TYPE :: electric_field_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: electric_potential
    LOGICAL  :: dipole_correction_ispresent = .FALSE.
    LOGICAL :: dipole_correction
    LOGICAL  :: gate_settings_ispresent = .FALSE.
    TYPE(gate_settings_type) :: gate_settings
    LOGICAL  :: electric_field_direction_ispresent = .FALSE.
    INTEGER :: electric_field_direction
    LOGICAL  :: potential_max_position_ispresent = .FALSE.
    REAL(DP) :: potential_max_position
    LOGICAL  :: potential_decrease_width_ispresent = .FALSE.
    REAL(DP) :: potential_decrease_width
    LOGICAL  :: electric_field_amplitude_ispresent = .FALSE.
    REAL(DP) :: electric_field_amplitude
    LOGICAL  :: electric_field_vector_ispresent = .FALSE.
    REAL(DP), DIMENSION(3) :: electric_field_vector
    LOGICAL  :: nk_per_string_ispresent = .FALSE.
    INTEGER :: nk_per_string
    LOGICAL  :: n_berry_cycles_ispresent = .FALSE.
    INTEGER :: n_berry_cycles
    !
  END TYPE electric_field_type
  !
  TYPE :: atomic_constraints_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: num_of_constraints
    REAL(DP) :: tolerance
    TYPE(atomic_constraint_type), DIMENSION(:), ALLOCATABLE :: atomic_constraint
    INTEGER   :: ndim_atomic_constraint
    !
  END TYPE atomic_constraints_type
  !
  TYPE :: BerryPhaseOutput_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(polarization_type) :: totalPolarization
    TYPE(phase_type) :: totalPhase
    TYPE(ionicPolarization_type), DIMENSION(:), ALLOCATABLE :: ionicPolarization
    INTEGER   :: ndim_ionicPolarization
    TYPE(electronicPolarization_type), DIMENSION(:), ALLOCATABLE :: electronicPolarization
    INTEGER   :: ndim_electronicPolarization
    !
  END TYPE BerryPhaseOutput_type
  !
  TYPE :: convergence_info_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(scf_conv_type) :: scf_conv
    LOGICAL  :: opt_conv_ispresent = .FALSE.
    TYPE(opt_conv_type) :: opt_conv
    LOGICAL  :: wf_collected_ispresent = .FALSE.
    LOGICAL :: wf_collected
    !
  END TYPE convergence_info_type
  !
  TYPE :: symmetries_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: nsym
    LOGICAL  :: colin_mag_ispresent = .FALSE.
    INTEGER :: colin_mag
    INTEGER :: nrot
    INTEGER :: space_group
    TYPE(symmetry_type), DIMENSION(:), ALLOCATABLE :: symmetry
    INTEGER   :: ndim_symmetry
    !
  END TYPE symmetries_type
  !
  TYPE :: outputPBC_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: assume_isolated
    LOGICAL  :: esm_ispresent = .FALSE.
    TYPE(esm_type) :: esm
    !
  END TYPE outputPBC_type
  !
  TYPE :: magnetization_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL :: lsda
    LOGICAL :: noncolin
    LOGICAL :: spinorbit
    LOGICAL  :: total_ispresent = .FALSE.
    REAL(DP) :: total
    LOGICAL  :: total_vec_ispresent = .FALSE.
    REAL(DP), DIMENSION(3) :: total_vec
    REAL(DP) :: absolute
    LOGICAL  :: Scalar_Site_Magnetic_Moments_ispresent = .FALSE.
    TYPE(scalmags_type) :: Scalar_Site_Magnetic_Moments
    LOGICAL  :: Site_Magnetizations_ispresent = .FALSE.
    TYPE(d3mags_type) :: Site_Magnetizations
    LOGICAL  :: do_magnetization_ispresent = .FALSE.
    LOGICAL :: do_magnetization
    !
  END TYPE magnetization_type
  !
  TYPE :: band_structure_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL :: lsda
    LOGICAL :: noncolin
    LOGICAL :: spinorbit
    LOGICAL  :: nbnd_ispresent = .FALSE.
    INTEGER :: nbnd
    LOGICAL  :: nbnd_up_ispresent = .FALSE.
    INTEGER :: nbnd_up
    LOGICAL  :: nbnd_dw_ispresent = .FALSE.
    INTEGER :: nbnd_dw
    REAL(DP) :: nelec
    LOGICAL  :: fermi_energy_ispresent = .FALSE.
    REAL(DP) :: fermi_energy
    LOGICAL  :: highestOccupiedLevel_ispresent = .FALSE.
    REAL(DP) :: highestOccupiedLevel
    LOGICAL  :: lowestUnoccupiedLevel_ispresent = .FALSE.
    REAL(DP) :: lowestUnoccupiedLevel
    LOGICAL  :: two_fermi_energies_ispresent = .FALSE.
    REAL(DP), DIMENSION(2) :: two_fermi_energies
    TYPE(k_points_IBZ_type) :: starting_k_points
    INTEGER :: nks
    TYPE(occupations_type) :: occupations_kind
    LOGICAL  :: smearing_ispresent = .FALSE.
    TYPE(smearing_type) :: smearing
    TYPE(ks_energies_type), DIMENSION(:), ALLOCATABLE :: ks_energies
    INTEGER   :: ndim_ks_energies
    !
  END TYPE band_structure_type
  !
  TYPE :: cpstep_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: ACCUMULATORS_ispresent = .FALSE.
    REAL(DP), DIMENSION(:), ALLOCATABLE :: ACCUMULATORS
    TYPE(cp_ionPos_type) :: IONS_POSITIONS
    TYPE(cp_ionsNose_type) :: IONS_NOSE
    LOGICAL  :: ekincm_ispresent = .FALSE.
    REAL(DP) :: ekincm
    TYPE(cp_elecNose_type) :: ELECTRONS_NOSE
    TYPE(cp_cell_type) :: CELL_PARAMETERS
    TYPE(cp_cellNose_type) :: CELL_NOSE
    !
  END TYPE cpstep_type
  !
  TYPE :: rism3d_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: nmol
    LOGICAL  :: molec_dir_ispresent = .FALSE.
    CHARACTER(len=256) :: molec_dir
    TYPE(solvent_type), DIMENSION(:), ALLOCATABLE :: solvent
    INTEGER   :: ndim_solvent
    REAL(DP) :: ecutsolv
    !
  END TYPE rism3d_type
  !
  TYPE :: step_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: n_step
    LOGICAL :: n_step_ispresent = .FALSE.
    TYPE(scf_conv_type) :: scf_conv
    TYPE(atomic_structure_type) :: atomic_structure
    TYPE(total_energy_type) :: total_energy
    TYPE(matrix_type) :: forces
    LOGICAL  :: stress_ispresent = .FALSE.
    TYPE(matrix_type) :: stress
    LOGICAL  :: fcp_force_ispresent = .FALSE.
    REAL(DP) :: fcp_force
    LOGICAL  :: fcp_tot_charge_ispresent = .FALSE.
    REAL(DP) :: fcp_tot_charge
    !
  END TYPE step_type
  !
  TYPE :: dft_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: functional
    LOGICAL  :: hybrid_ispresent = .FALSE.
    TYPE(hybrid_type) :: hybrid
    LOGICAL  :: dftU_ispresent = .FALSE.
    TYPE(dftU_type) :: dftU
    LOGICAL  :: vdW_ispresent = .FALSE.
    TYPE(vdW_type) :: vdW
    !
  END TYPE dft_type
  !
  TYPE :: outputElectricField_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: BerryPhase_ispresent = .FALSE.
    TYPE(BerryPhaseOutput_type) :: BerryPhase
    LOGICAL  :: finiteElectricFieldInfo_ispresent = .FALSE.
    TYPE(finiteFieldOut_type) :: finiteElectricFieldInfo
    LOGICAL  :: sawtoothEnergy_ispresent = .FALSE.
    TYPE(sawtoothEnergy_type) :: sawtoothEnergy
    LOGICAL  :: dipoleInfo_ispresent = .FALSE.
    TYPE(dipoleOutput_type) :: dipoleInfo
    LOGICAL  :: gateInfo_ispresent = .FALSE.
    TYPE(gateInfo_type) :: gateInfo
    !
  END TYPE outputElectricField_type
  !
  TYPE :: cptimesteps_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    INTEGER :: nt
    LOGICAL :: nt_ispresent = .FALSE.
    TYPE(cpstep_type) :: STEP0
    TYPE(cpstep_type) :: STEPM
    !
  END TYPE cptimesteps_type
  !
  TYPE :: input_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    TYPE(control_variables_type) :: control_variables
    TYPE(atomic_species_type) :: atomic_species
    TYPE(atomic_structure_type) :: atomic_structure
    TYPE(dft_type) :: dft
    TYPE(spin_type) :: spin
    TYPE(bands_type) :: bands
    TYPE(basis_type) :: basis
    TYPE(electron_control_type) :: electron_control
    TYPE(k_points_IBZ_type) :: k_points_IBZ
    TYPE(ion_control_type) :: ion_control
    TYPE(cell_control_type) :: cell_control
    LOGICAL  :: symmetry_flags_ispresent = .FALSE.
    TYPE(symmetry_flags_type) :: symmetry_flags
    LOGICAL  :: boundary_conditions_ispresent = .FALSE.
    TYPE(boundary_conditions_type) :: boundary_conditions
    LOGICAL  :: fcp_settings_ispresent = .FALSE.
    TYPE(fcp_type) :: fcp_settings
    LOGICAL  :: rism_settings_ispresent = .FALSE.
    TYPE(rism_type) :: rism_settings
    LOGICAL  :: solvents_ispresent = .FALSE.
    TYPE(solvents_type) :: solvents
    LOGICAL  :: ekin_functional_ispresent = .FALSE.
    TYPE(ekin_functional_type) :: ekin_functional
    LOGICAL  :: external_atomic_forces_ispresent = .FALSE.
    TYPE(matrix_type) :: external_atomic_forces
    LOGICAL  :: free_positions_ispresent = .FALSE.
    TYPE(integerMatrix_type) :: free_positions
    LOGICAL  :: starting_atomic_velocities_ispresent = .FALSE.
    TYPE(matrix_type) :: starting_atomic_velocities
    LOGICAL  :: electric_field_ispresent = .FALSE.
    TYPE(electric_field_type) :: electric_field
    LOGICAL  :: atomic_constraints_ispresent = .FALSE.
    TYPE(atomic_constraints_type) :: atomic_constraints
    LOGICAL  :: spin_constraints_ispresent = .FALSE.
    TYPE(spin_constraints_type) :: spin_constraints
    LOGICAL  :: twoch__ispresent = .FALSE.
    TYPE(two_chem_type) :: twoch_
    !
  END TYPE input_type
  !
  TYPE :: output_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    LOGICAL  :: convergence_info_ispresent = .FALSE.
    TYPE(convergence_info_type) :: convergence_info
    TYPE(algorithmic_info_type) :: algorithmic_info
    TYPE(atomic_species_type) :: atomic_species
    TYPE(atomic_structure_type) :: atomic_structure
    LOGICAL  :: symmetries_ispresent = .FALSE.
    TYPE(symmetries_type) :: symmetries
    TYPE(basis_set_type) :: basis_set
    TYPE(dft_type) :: dft
    LOGICAL  :: boundary_conditions_ispresent = .FALSE.
    TYPE(outputPBC_type) :: boundary_conditions
    LOGICAL  :: magnetization_ispresent = .FALSE.
    TYPE(magnetization_type) :: magnetization
    TYPE(total_energy_type) :: total_energy
    TYPE(band_structure_type) :: band_structure
    LOGICAL  :: forces_ispresent = .FALSE.
    TYPE(matrix_type) :: forces
    LOGICAL  :: stress_ispresent = .FALSE.
    TYPE(matrix_type) :: stress
    LOGICAL  :: electric_field_ispresent = .FALSE.
    TYPE(outputElectricField_type) :: electric_field
    LOGICAL  :: fcp_force_ispresent = .FALSE.
    REAL(DP) :: fcp_force
    LOGICAL  :: fcp_tot_charge_ispresent = .FALSE.
    REAL(DP) :: fcp_tot_charge
    LOGICAL  :: rism3d_ispresent = .FALSE.
    TYPE(rism3d_type) :: rism3d
    LOGICAL  :: rismlaue_ispresent = .FALSE.
    TYPE(rismlaue_type) :: rismlaue
    LOGICAL  :: two_chem_ispresent = .FALSE.
    TYPE(two_chem_type) :: two_chem
    !
  END TYPE output_type
  !
  TYPE :: espresso_type
    !
    CHARACTER(len=100) :: tagname
    LOGICAL  :: lwrite = .FALSE.
    LOGICAL  :: lread  = .FALSE.
    !
    CHARACTER(len=256) :: Units
    LOGICAL :: Units_ispresent = .FALSE.
    LOGICAL  :: general_info_ispresent = .FALSE.
    TYPE(general_info_type) :: general_info
    LOGICAL  :: parallel_info_ispresent = .FALSE.
    TYPE(parallel_info_type) :: parallel_info
    LOGICAL  :: input_ispresent = .FALSE.
    TYPE(input_type) :: input
    LOGICAL  :: step_ispresent = .FALSE.
    TYPE(step_type), DIMENSION(:), ALLOCATABLE :: step
    INTEGER   :: ndim_step
    LOGICAL  :: output_ispresent = .FALSE.
    TYPE(output_type) :: output
    LOGICAL  :: STATUS_ispresent = .FALSE.
    TYPE(cpstatus_type) :: STATUS
    LOGICAL  :: TIMESTEPS_ispresent = .FALSE.
    TYPE(cptimesteps_type) :: TIMESTEPS
    LOGICAL  :: exit_status_ispresent = .FALSE.
    INTEGER :: exit_status
    LOGICAL  :: cputime_ispresent = .FALSE.
    INTEGER :: cputime
    LOGICAL  :: timing_info_ispresent = .FALSE.
    TYPE(timing_type) :: timing_info
    LOGICAL  :: closed_ispresent = .FALSE.
    TYPE(closed_type) :: closed
    !
  END TYPE espresso_type
  !
  !
END MODULE qes_types_module