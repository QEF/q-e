!
! Copyright (C) 2003-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE qes_types_module

! This module contains the data_types used for describing
! the XML files produced by Quantum ESPRESSO package.
!
! Written by Giovanni Borghi, A. Ferretti, ... (2015).
!
USE kinds, only: DP

TYPE :: closed_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: DATE
   CHARACTER(len=256) :: TIME
   CHARACTER(len=256) :: closed
   !
END TYPE closed_type

TYPE :: status_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
   INTEGER  :: status
   !
END TYPE status_type

TYPE :: scalarQuantity_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: Units
   REAL(DP) :: scalarQuantity
   !
END TYPE scalarQuantity_type

TYPE :: finiteFieldOut_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   REAL(DP), DIMENSION(3) :: electronicDipole
   REAL(DP), DIMENSION(3) :: ionicDipole
   !
END TYPE finiteFieldOut_type

TYPE :: k_point_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: weight_ispresent
   REAL(DP) :: weight
   LOGICAL  :: label_ispresent
   CHARACTER(len=256) :: label
   REAL(DP), DIMENSION(3) :: k_point
   !
END TYPE k_point_type

TYPE :: atom_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: name
   LOGICAL  :: position_ispresent
   CHARACTER(len=256) :: position
   LOGICAL  :: index_ispresent
   INTEGER  :: index
   REAL(DP), DIMENSION(3) :: atom
   !
END TYPE atom_type

TYPE :: phase_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: ionic_ispresent
   REAL(DP) :: ionic
   LOGICAL  :: electronic_ispresent
   REAL(DP) :: electronic
   LOGICAL  :: modulus_ispresent
   CHARACTER(len=256) :: modulus
   REAL(DP) :: phase
   !
END TYPE phase_type

TYPE :: dipoleOutput_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: idir
   TYPE(scalarQuantity_type) :: dipole
   TYPE(scalarQuantity_type) :: ion_dipole
   TYPE(scalarQuantity_type) :: elec_dipole
   TYPE(scalarQuantity_type) :: dipoleField
   TYPE(scalarQuantity_type) :: potentialAmp
   TYPE(scalarQuantity_type) :: totalLength
   !
END TYPE dipoleOutput_type

TYPE :: polarization_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   TYPE(scalarQuantity_type) :: polarization
   REAL(DP) :: modulus
   REAL(DP), DIMENSION(3) :: direction
   !
END TYPE polarization_type

TYPE :: ionicPolarization_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   TYPE(atom_type) :: ion
   REAL(DP) :: charge
   TYPE(phase_type) :: phase
   !
END TYPE ionicPolarization_type

TYPE :: electronicPolarization_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   TYPE(k_point_type) :: firstKeyPoint
   LOGICAL  :: spin_ispresent
   INTEGER  :: spin
   TYPE(phase_type) :: phase
   !
END TYPE electronicPolarization_type

TYPE :: BerryPhaseOutput_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   TYPE(polarization_type) :: polarization
   TYPE(phase_type) :: totalPhase
   TYPE(ionicPolarization_type), DIMENSION(:), ALLOCATABLE :: ionicPolarization
   TYPE(electronicPolarization_type), DIMENSION(:), ALLOCATABLE :: electronicPolarization
   !
   INTEGER  :: ndim_ionicPolarization
   INTEGER  :: ndim_electronicPolarization

   !
END TYPE BerryPhaseOutput_type

TYPE :: vector_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
   INTEGER  :: ndim_vec
   REAL(DP), DIMENSION(:), ALLOCATABLE :: vec
   !
END TYPE vector_type

TYPE :: ks_energies_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   TYPE(k_point_type) :: k_point
   INTEGER  :: npw
   !
   INTEGER  :: ndim_eigenvalues
   REAL(DP), DIMENSION(:), ALLOCATABLE :: eigenvalues
   !
   INTEGER  :: ndim_occupations
   REAL(DP), DIMENSION(:), ALLOCATABLE :: occupations
   !
END TYPE ks_energies_type

TYPE :: magnetization_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: lsda
   LOGICAL  :: noncolin
   LOGICAL  :: spinorbit
   REAL(DP) :: total
   REAL(DP) :: absolute
   LOGICAL  :: do_magnetization
   !
END TYPE magnetization_type

TYPE :: reciprocal_lattice_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   REAL(DP), DIMENSION(3) :: b1
   REAL(DP), DIMENSION(3) :: b2
   REAL(DP), DIMENSION(3) :: b3
   !
END TYPE reciprocal_lattice_type

TYPE :: basisSetItem_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: nr1
   INTEGER  :: nr2
   INTEGER  :: nr3
   CHARACTER(len=256) :: basisSetItem
   !
END TYPE basisSetItem_type

TYPE :: equivalent_atoms_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: nat
   !
   INTEGER  :: ndim_index_list
   INTEGER, DIMENSION(:), ALLOCATABLE :: index_list
   !
END TYPE equivalent_atoms_type

TYPE :: info_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: name_ispresent
   CHARACTER(len=256) :: name
   LOGICAL  :: class_ispresent
   CHARACTER(len=256) :: class
   LOGICAL  :: time_reversal_ispresent
   LOGICAL  :: time_reversal
   CHARACTER(len=256) :: info
   !
END TYPE info_type

TYPE :: matrix_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
   INTEGER  :: ndim1_mat
   INTEGER  :: ndim2_mat
   REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mat
   !
END TYPE matrix_type

TYPE :: symmetry_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   TYPE(info_type) :: info
   TYPE(matrix_type) :: rotation
   LOGICAL  :: fractional_translation_ispresent
   REAL(DP), DIMENSION(3) :: fractional_translation
   LOGICAL  :: equivalent_atoms_ispresent
   TYPE(equivalent_atoms_type) :: equivalent_atoms
   !
END TYPE symmetry_type

TYPE :: algorithmic_info_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: real_space_q
   LOGICAL  :: uspp
   LOGICAL  :: paw
   !
END TYPE algorithmic_info_type

TYPE :: opt_conv_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: n_opt_steps
   REAL(DP) :: grad_norm
   !
END TYPE opt_conv_type

TYPE :: scf_conv_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: n_scf_steps
   REAL(DP) :: scf_error
   !
END TYPE scf_conv_type

TYPE :: species_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: name
   LOGICAL  :: mass_ispresent
   REAL(DP) :: mass
   CHARACTER(len=256) :: pseudo_file
   LOGICAL  :: starting_magnetization_ispresent
   REAL(DP) :: starting_magnetization
   LOGICAL  :: spin_teta_ispresent
   REAL(DP) :: spin_teta
   LOGICAL  :: spin_phi_ispresent
   REAL(DP) :: spin_phi
   !
END TYPE species_type

TYPE :: total_energy_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   REAL(DP) :: etot
   LOGICAL  :: eband_ispresent
   REAL(DP) :: eband
   LOGICAL  :: ehart_ispresent
   REAL(DP) :: ehart
   LOGICAL  :: vtxc_ispresent
   REAL(DP) :: vtxc
   LOGICAL  :: etxc_ispresent
   REAL(DP) :: etxc
   LOGICAL  :: ewald_ispresent
   REAL(DP) :: ewald
   LOGICAL  :: demet_ispresent
   REAL(DP) :: demet
   LOGICAL  :: efieldcorr_ispresent
   REAL(DP) :: efieldcorr
   LOGICAL  :: potentiostat_contr_ispresent
   REAL(DP) :: potentiostat_contr
   !
END TYPE total_energy_type

TYPE :: convergence_info_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   TYPE(scf_conv_type) :: scf_conv
   LOGICAL  :: opt_conv_ispresent
   TYPE(opt_conv_type) :: opt_conv
   !
END TYPE convergence_info_type

TYPE :: outputElectricField_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: BerryPhase_ispresent
   TYPE(BerryPhaseOutput_type) :: BerryPhase
   LOGICAL  :: finiteElectricFieldInfo_ispresent
   TYPE(finiteFieldOut_type) :: finiteElectricFieldInfo
   LOGICAL  :: dipoleInfo_ispresent
   TYPE(dipoleOutput_type) :: dipoleInfo
   !
END TYPE outputElectricField_type

TYPE :: spin_constraints_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: spin_constraints
   REAL(DP) :: lagrange_multiplier
   LOGICAL  :: target_magnetization_ispresent
   REAL(DP), DIMENSION(3) :: target_magnetization
   !
END TYPE spin_constraints_type

TYPE :: constr_type_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
END TYPE constr_type_type

TYPE :: constr_parms_list_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
END TYPE constr_parms_list_type

TYPE :: atomic_constraint_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   REAL(DP), DIMENSION(4) :: constr_parms
   CHARACTER(len=256) :: constr_type
   REAL(DP) :: constr_target
   !
END TYPE atomic_constraint_type

TYPE :: atomic_constraints_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: num_of_constraints
   REAL(DP) :: tolerance
   TYPE(atomic_constraint_type), DIMENSION(:), ALLOCATABLE :: atomic_constraint
   !
   INTEGER  :: ndim_atomic_constraint

   !
END TYPE atomic_constraints_type

TYPE :: electric_potential_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
END TYPE electric_potential_type

TYPE :: electric_field_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: electric_potential
   LOGICAL  :: dipole_correction_ispresent
   LOGICAL  :: dipole_correction
   LOGICAL  :: electric_field_direction_ispresent
   INTEGER  :: electric_field_direction
   LOGICAL  :: potential_max_position_ispresent
   REAL(DP) :: potential_max_position
   LOGICAL  :: potential_decrease_width_ispresent
   REAL(DP) :: potential_decrease_width
   LOGICAL  :: electric_field_amplitude_ispresent
   REAL(DP) :: electric_field_amplitude
   LOGICAL  :: electric_field_vector_ispresent
   REAL(DP), DIMENSION(3) :: electric_field_vector
   LOGICAL  :: nk_per_string_ispresent
   INTEGER  :: nk_per_string
   LOGICAL  :: n_berry_cycles_ispresent
   INTEGER  :: n_berry_cycles
   !
END TYPE electric_field_type

TYPE :: symmetries_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER :: nsym
   INTEGER :: nrot
   INTEGER :: space_group
   TYPE(symmetry_type), DIMENSION(:), ALLOCATABLE :: symmetry
   !
   INTEGER  :: ndim_symmetry

   !
END TYPE symmetries_type

TYPE :: ekin_functional_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   REAL(DP) :: ecfixed
   REAL(DP) :: qcutz
   REAL(DP) :: q2sigma
   !
END TYPE ekin_functional_type

TYPE :: esm_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: bc
   INTEGER  :: nfit
   REAL(DP) :: w
   REAL(DP) :: efield
   !
END TYPE esm_type

TYPE :: boundary_conditions_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: assume_isolated
   LOGICAL  :: esm_ispresent
   TYPE(esm_type) :: esm
   LOGICAL  :: fcp_opt_ispresent
   LOGICAL  :: fcp_opt
   LOGICAL  :: fcp_mu_ispresent
   REAL(DP) :: fcp_mu
   !
END TYPE boundary_conditions_type

TYPE :: symmetry_flags_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: nosym
   LOGICAL  :: nosym_evc
   LOGICAL  :: noinv
   LOGICAL  :: no_t_rev
   LOGICAL  :: force_symmorphic
   LOGICAL  :: use_all_frac
   !
END TYPE symmetry_flags_type

TYPE :: integerMatrix_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
   INTEGER  :: ndim1_int_mat
   INTEGER  :: ndim2_int_mat
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: int_mat
   !
END TYPE integerMatrix_type

TYPE :: cell_control_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: cell_dynamics
   REAL(DP) :: pressure
   LOGICAL  :: wmass_ispresent
   REAL(DP) :: wmass
   LOGICAL  :: cell_factor_ispresent
   REAL(DP) :: cell_factor
   LOGICAL  :: fix_volume_ispresent
   LOGICAL  :: fix_volume
   LOGICAL  :: fix_area_ispresent
   LOGICAL  :: fix_area
   LOGICAL  :: isotropic_ispresent
   LOGICAL  :: isotropic
   LOGICAL  :: free_cell_ispresent
   TYPE(integerMatrix_type) :: free_cell
   !
END TYPE cell_control_type

TYPE :: md_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: pot_extrapolation
   CHARACTER(len=256) :: wfc_extrapolation
   CHARACTER(len=256) :: ion_temperature
   REAL(DP) :: timestep
   REAL(DP) :: tempw
   REAL(DP) :: tolp
   REAL(DP) :: deltaT
   INTEGER  :: nraise
   !
END TYPE md_type

TYPE :: bfgs_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: ndim
   REAL(DP) :: trust_radius_min
   REAL(DP) :: trust_radius_max
   REAL(DP) :: trust_radius_init
   REAL(DP) :: w1
   REAL(DP) :: w2
   !
END TYPE bfgs_type

TYPE :: ion_control_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: ion_dynamics
   LOGICAL  :: upscale_ispresent
   REAL(DP) :: upscale
   LOGICAL  :: remove_rigid_rot_ispresent
   LOGICAL  :: remove_rigid_rot
   LOGICAL  :: refold_pos_ispresent
   LOGICAL  :: refold_pos
   LOGICAL  :: bfgs_ispresent
   TYPE(bfgs_type) :: bfgs
   LOGICAL  :: md_ispresent
   TYPE(md_type) :: md
   !
END TYPE ion_control_type

TYPE :: band_structure_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: lsda
   LOGICAL  :: noncolin
   LOGICAL  :: spinorbit
   INTEGER  :: nbnd
   LOGICAL  :: nbnd_up_ispresent
   INTEGER  :: nbnd_up
   LOGICAL  :: nbnd_dw_ispresent
   INTEGER  :: nbnd_dw
   REAL(DP) :: nelec
   LOGICAL  :: num_of_atomic_wfc_ispresent
   INTEGER  :: num_of_atomic_wfc
   LOGICAL  :: fermi_energy_ispresent
   REAL(DP) :: fermi_energy
   LOGICAL  :: highestOccupiedLevel_ispresent
   REAL(DP) :: highestOccupiedLevel
   LOGICAL  :: two_fermi_energies_ispresent
   !
   INTEGER  :: ndim_two_fermi_energies
   REAL(DP), DIMENSION(:), ALLOCATABLE :: two_fermi_energies
   INTEGER  :: nks
   TYPE(ks_energies_type), DIMENSION(:), ALLOCATABLE :: ks_energies
   !
   INTEGER  :: ndim_ks_energies

   !
END TYPE band_structure_type

TYPE :: monkhorst_pack_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: nk1
   INTEGER  :: nk2
   INTEGER  :: nk3
   INTEGER  :: k1
   INTEGER  :: k2
   INTEGER  :: k3
   CHARACTER(len=256) :: monkhorst_pack
   !
END TYPE monkhorst_pack_type

TYPE :: k_points_IBZ_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: monkhorst_pack_ispresent
   TYPE(monkhorst_pack_type) :: monkhorst_pack
   LOGICAL  :: nk_ispresent
   INTEGER  :: nk
   LOGICAL  :: k_point_ispresent
   TYPE(k_point_type), DIMENSION(:), ALLOCATABLE :: k_point
   !
   INTEGER  :: ndim_k_point

   !
END TYPE k_points_IBZ_type

TYPE :: mixingMode_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
END TYPE mixingMode_type

TYPE :: diago_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
END TYPE diago_type

TYPE :: electron_control_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: diagonalization
   CHARACTER(len=256) :: mixing_mode
   REAL(DP) :: mixing_beta
   REAL(DP) :: conv_thr
   INTEGER  :: mixing_ndim
   INTEGER  :: max_nstep
   LOGICAL  :: real_space_q
   LOGICAL  :: tq_smoothing
   LOGICAL  :: tbeta_smoothing
   REAL(DP) :: diago_thr_init
   LOGICAL  :: diago_full_acc
   INTEGER  :: diago_cg_maxiter
   !
END TYPE electron_control_type

TYPE :: basis_set_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: gamma_only_ispresent
   LOGICAL  :: gamma_only
   REAL(DP) :: ecutwfc
   LOGICAL  :: ecutrho_ispresent
   REAL(DP) :: ecutrho
   TYPE(basisSetItem_type) :: fft_grid
   LOGICAL  :: fft_smooth_ispresent
   TYPE(basisSetItem_type) :: fft_smooth
   LOGICAL  :: fft_box_ispresent
   TYPE(basisSetItem_type) :: fft_box
   INTEGER  :: ngm
   LOGICAL  :: ngms_ispresent
   INTEGER  :: ngms
   INTEGER  :: npwx
   TYPE(reciprocal_lattice_type) :: reciprocal_lattice
   !
END TYPE basis_set_type

TYPE :: basis_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: gamma_only_ispresent
   LOGICAL  :: gamma_only
   REAL(DP) :: ecutwfc
   LOGICAL  :: ecutrho_ispresent
   REAL(DP) :: ecutrho
   LOGICAL  :: fft_grid_ispresent
   TYPE(basisSetItem_type) :: fft_grid
   LOGICAL  :: fft_smooth_ispresent
   TYPE(basisSetItem_type) :: fft_smooth
   LOGICAL  :: fft_box_ispresent
   TYPE(basisSetItem_type) :: fft_box
   !
END TYPE basis_type

TYPE :: inputOccupations_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: ispin
   REAL(DP) :: spin_factor
   !
   INTEGER  :: ndim_vec
   REAL(DP), DIMENSION(:), ALLOCATABLE :: vec
   !
END TYPE inputOccupations_type

TYPE :: occupations_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: spin_ispresent
   INTEGER  :: spin
   CHARACTER(len=256) :: occupations
   !
END TYPE occupations_type

TYPE :: smearing_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   REAL(DP) :: degauss
   CHARACTER(len=256) :: smearing
   !
END TYPE smearing_type

TYPE :: bands_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: nbnd_ispresent
   INTEGER  :: nbnd
   LOGICAL  :: smearing_ispresent
   TYPE(smearing_type) :: smearing
   LOGICAL  :: tot_charge_ispresent
   REAL(DP) :: tot_charge
   LOGICAL  :: tot_magnetization_ispresent
   REAL(DP) :: tot_magnetization
   TYPE(occupations_type) :: occupations
   LOGICAL  :: inputOccupations_ispresent
   TYPE(inputOccupations_type), DIMENSION(:), ALLOCATABLE :: inputOccupations
   !
   INTEGER  :: ndim_inputOccupations

   !
END TYPE bands_type

TYPE :: spin_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: lsda
   LOGICAL  :: noncolin
   LOGICAL  :: spinorbit
   !
END TYPE spin_type

TYPE :: HubbardCommon_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: specie
   CHARACTER(len=256) :: label
   REAL(DP) :: HubbardCommon
   !
END TYPE HubbardCommon_type

TYPE :: HubbardProj_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
END TYPE HubbardProj_type

TYPE :: Hubbard_ns_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: specie
   CHARACTER(len=256) :: label
   INTEGER  :: spin
   INTEGER  :: index
   !
   INTEGER  :: ndim1_mat
   INTEGER  :: ndim2_mat
   REAL(DP), DIMENSION(:,:), ALLOCATABLE :: mat
   !
END TYPE Hubbard_ns_type

TYPE :: starting_ns_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: specie
   CHARACTER(len=256) :: label
   INTEGER  :: spin
   !
   INTEGER  :: ndim_vec
   REAL(DP), DIMENSION(:), ALLOCATABLE :: vec
   !
END TYPE starting_ns_type

TYPE :: HubbardJ_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: specie
   CHARACTER(len=256) :: label
   REAL(DP), DIMENSION(3) :: HubbardJ
   !
END TYPE HubbardJ_type

TYPE :: vdW_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: vdw_corr
   LOGICAL  :: non_local_term_ispresent
   CHARACTER(len=256) :: non_local_term
   LOGICAL  :: london_s6_ispresent
   REAL(DP) :: london_s6
   LOGICAL  :: ts_vdw_econv_thr_ispresent
   REAL(DP) :: ts_vdw_econv_thr
   LOGICAL  :: ts_vdw_isolated_ispresent
   LOGICAL  :: ts_vdw_isolated
   LOGICAL  :: london_rcut_ispresent
   REAL(DP) :: london_rcut
   LOGICAL  :: xdm_a1_ispresent
   REAL(DP) :: xdm_a1
   LOGICAL  :: xdm_a2_ispresent
   REAL(DP) :: xdm_a2
   LOGICAL  :: london_c6_ispresent
   TYPE(HubbardCommon_type), DIMENSION(:), ALLOCATABLE :: london_c6
   !
   INTEGER  :: ndim_london_c6

   !
END TYPE vdW_type

TYPE :: dftU_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   LOGICAL  :: lda_plus_u_kind_ispresent
   INTEGER :: lda_plus_u_kind
   LOGICAL  :: Hubbard_U_ispresent
   TYPE(HubbardCommon_type), DIMENSION(:), ALLOCATABLE :: Hubbard_U
   LOGICAL  :: Hubbard_J0_ispresent
   TYPE(HubbardCommon_type), DIMENSION(:), ALLOCATABLE :: Hubbard_J0
   LOGICAL  :: Hubbard_alpha_ispresent
   TYPE(HubbardCommon_type), DIMENSION(:), ALLOCATABLE :: Hubbard_alpha
   LOGICAL  :: Hubbard_beta_ispresent
   TYPE(HubbardCommon_type), DIMENSION(:), ALLOCATABLE :: Hubbard_beta
   LOGICAL  :: Hubbard_J_ispresent
   TYPE(HubbardJ_type), DIMENSION(:), ALLOCATABLE :: Hubbard_J
   LOGICAL  :: starting_ns_ispresent
   TYPE(starting_ns_type), DIMENSION(:), ALLOCATABLE :: starting_ns
   LOGICAL  :: Hubbard_ns_ispresent
   TYPE(Hubbard_ns_type), DIMENSION(:), ALLOCATABLE :: Hubbard_ns
   LOGICAL  :: U_projection_type_ispresent
   CHARACTER(len=256) :: U_projection_type
   !
   INTEGER  :: ndim_Hubbard_U
   INTEGER  :: ndim_Hubbard_J0
   INTEGER  :: ndim_Hubbard_alpha
   INTEGER  :: ndim_Hubbard_beta
   INTEGER  :: ndim_Hubbard_J
   INTEGER  :: ndim_starting_ns
   INTEGER  :: ndim_Hubbard_ns

   !
END TYPE dftU_type

TYPE :: qpoint_grid_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: nqx1
   INTEGER  :: nqx2
   INTEGER  :: nqx3
   CHARACTER(len=256) :: qpoint_grid
   !
END TYPE qpoint_grid_type

TYPE :: hybrid_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   TYPE(qpoint_grid_type) :: qpoint_grid
   REAL(DP) :: ecutfock
   REAL(DP) :: exx_fraction
   REAL(DP) :: screening_parameter
   CHARACTER(len=256) :: exxdiv_treatment
   LOGICAL  :: x_gamma_extrapolation
   REAL(DP) :: ecutvcut
   !
END TYPE hybrid_type

TYPE :: functional_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
END TYPE functional_type

TYPE :: dft_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: functional
   LOGICAL  :: hybrid_ispresent
   TYPE(hybrid_type) :: hybrid
   LOGICAL  :: dftU_ispresent
   TYPE(dftU_type) :: dftU
   LOGICAL  :: vdW_ispresent
   TYPE(vdW_type) :: vdW
   !
END TYPE dft_type

TYPE :: d3vector_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
END TYPE d3vector_type

TYPE :: cell_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   REAL(DP), DIMENSION(3) :: a1
   REAL(DP), DIMENSION(3) :: a2
   REAL(DP), DIMENSION(3) :: a3
   !
END TYPE cell_type

TYPE :: wyckoff_positions_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER :: space_group
   LOGICAL  :: more_options_ispresent
   CHARACTER(len=256) :: more_options
   TYPE(atom_type), DIMENSION(:), ALLOCATABLE :: atom
   !
   INTEGER  :: ndim_atom

   !
END TYPE wyckoff_positions_type

TYPE :: atomic_positions_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   TYPE(atom_type), DIMENSION(:), ALLOCATABLE :: atom
   !
   INTEGER  :: ndim_atom

   !
END TYPE atomic_positions_type

TYPE :: atomic_structure_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: nat
   LOGICAL  :: alat_ispresent
   REAL(DP) :: alat
   LOGICAL  :: bravais_index_ispresent
   INTEGER  :: bravais_index
   LOGICAL  :: atomic_positions_ispresent
   TYPE(atomic_positions_type) :: atomic_positions
   LOGICAL  :: wyckoff_positions_ispresent
   TYPE(wyckoff_positions_type) :: wyckoff_positions
   TYPE(cell_type) :: cell
   !
END TYPE atomic_structure_type

TYPE :: step_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: n_step
   TYPE(scf_conv_type) :: scf_conv
   TYPE(atomic_structure_type) :: atomic_structure
   TYPE(total_energy_type) :: total_energy
   TYPE(matrix_type) :: forces
   LOGICAL  :: stress_ispresent
   TYPE(matrix_type) :: stress
   LOGICAL  :: FCP_force_ispresent
   REAL(DP) :: FCP_force
   LOGICAL  :: FCP_tot_charge_ispresent
   REAL(DP) :: FCP_tot_charge
   !
END TYPE step_type

TYPE :: atomic_species_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: ntyp
   TYPE(species_type), DIMENSION(:), ALLOCATABLE :: species
   !
   INTEGER  :: ndim_species

   !
END TYPE atomic_species_type

TYPE :: output_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   TYPE(convergence_info_type) :: convergence_info
   TYPE(algorithmic_info_type) :: algorithmic_info
   TYPE(atomic_species_type) :: atomic_species
   TYPE(atomic_structure_type) :: atomic_structure
   TYPE(symmetries_type) :: symmetries
   TYPE(basis_set_type) :: basis_set
   TYPE(dft_type) :: dft
   TYPE(magnetization_type) :: magnetization
   TYPE(total_energy_type) :: total_energy
   TYPE(band_structure_type) :: band_structure
   LOGICAL  :: forces_ispresent
   TYPE(matrix_type) :: forces
   LOGICAL  :: stress_ispresent
   TYPE(matrix_type) :: stress
   LOGICAL  :: electric_field_ispresent
   TYPE(outputElectricField_type) :: electric_field
   LOGICAL  :: FCP_force_ispresent
   REAL(DP) :: FCP_force
   LOGICAL  :: FCP_tot_charge_ispresent
   REAL(DP) :: FCP_tot_charge
   !
END TYPE output_type

TYPE :: lowhigh_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
END TYPE lowhigh_type

TYPE :: controlRestartMode_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
END TYPE controlRestartMode_type

TYPE :: calculation_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   !
END TYPE calculation_type

TYPE :: control_variables_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: title
   CHARACTER(len=256) :: calculation
   CHARACTER(len=256) :: restart_mode
   CHARACTER(len=256) :: prefix
   CHARACTER(len=256) :: pseudo_dir
   CHARACTER(len=256) :: outdir
   LOGICAL  :: stress
   LOGICAL  :: forces
   LOGICAL  :: wf_collect
   CHARACTER(len=256) :: disk_io
   INTEGER  :: max_seconds
   LOGICAL  :: nstep_ispresent
   INTEGER  :: nstep
   REAL(DP) :: etot_conv_thr
   REAL(DP) :: forc_conv_thr
   REAL(DP) :: press_conv_thr
   CHARACTER(len=256) :: verbosity
   INTEGER  :: print_every
   !
END TYPE control_variables_type

TYPE :: input_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
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
   LOGICAL  :: symmetry_flags_ispresent
   TYPE(symmetry_flags_type) :: symmetry_flags
   LOGICAL  :: boundary_conditions_ispresent
   TYPE(boundary_conditions_type) :: boundary_conditions
   LOGICAL  :: ekin_functional_ispresent
   TYPE(ekin_functional_type) :: ekin_functional
   LOGICAL  :: external_atomic_forces_ispresent
   TYPE(matrix_type) :: external_atomic_forces
   LOGICAL  :: free_positions_ispresent
   TYPE(integerMatrix_type) :: free_positions
   LOGICAL  :: starting_atomic_velocities_ispresent
   TYPE(matrix_type) :: starting_atomic_velocities
   LOGICAL  :: electric_field_ispresent
   TYPE(electric_field_type) :: electric_field
   LOGICAL  :: atomic_constraints_ispresent
   TYPE(atomic_constraints_type) :: atomic_constraints
   LOGICAL  :: spin_constraints_ispresent
   TYPE(spin_constraints_type) :: spin_constraints
   !
END TYPE input_type

TYPE :: parallel_info_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   INTEGER  :: nprocs
   INTEGER  :: nthreads
   INTEGER  :: ntasks
   INTEGER  :: nbgrp
   INTEGER  :: npool
   INTEGER  :: ndiag
   !
END TYPE parallel_info_type

TYPE :: created_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: DATE
   CHARACTER(len=256) :: TIME
   CHARACTER(len=256) :: created
   !
END TYPE created_type

TYPE :: creator_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: NAME
   CHARACTER(len=256) :: VERSION
   CHARACTER(len=256) :: creator
   !
END TYPE creator_type

TYPE :: xml_format_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   CHARACTER(len=256) :: NAME
   CHARACTER(len=256) :: VERSION
   CHARACTER(len=256) :: xml_format
   !
END TYPE xml_format_type

TYPE :: general_info_type
   !
   CHARACTER(len=100) :: tagname
   LOGICAL  :: lread = .true.
   LOGICAL  :: lwrite = .true.
   !
   TYPE(xml_format_type) :: xml_format
   TYPE(creator_type) :: creator
   TYPE(created_type) :: created
   CHARACTER(len=256) :: job
   !
END TYPE general_info_type

END MODULE qes_types_module
