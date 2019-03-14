!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE qes_reset_module
  !
  ! Auto-generated code: don't edit or at least don't commit changes
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes-1.0
  !
  USE qes_types_module
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC qes_reset
  !
  INTERFACE qes_reset
    MODULE PROCEDURE qes_reset_espresso
    MODULE PROCEDURE qes_reset_general_info
    MODULE PROCEDURE qes_reset_parallel_info
    MODULE PROCEDURE qes_reset_input
    MODULE PROCEDURE qes_reset_step
    MODULE PROCEDURE qes_reset_output
    MODULE PROCEDURE qes_reset_timing
    MODULE PROCEDURE qes_reset_clock
    MODULE PROCEDURE qes_reset_control_variables
    MODULE PROCEDURE qes_reset_xml_format
    MODULE PROCEDURE qes_reset_creator
    MODULE PROCEDURE qes_reset_created
    MODULE PROCEDURE qes_reset_atomic_species
    MODULE PROCEDURE qes_reset_species
    MODULE PROCEDURE qes_reset_atomic_structure
    MODULE PROCEDURE qes_reset_atomic_positions
    MODULE PROCEDURE qes_reset_atom
    MODULE PROCEDURE qes_reset_wyckoff_positions
    MODULE PROCEDURE qes_reset_cell
    MODULE PROCEDURE qes_reset_dft
    MODULE PROCEDURE qes_reset_hybrid
    MODULE PROCEDURE qes_reset_qpoint_grid
    MODULE PROCEDURE qes_reset_dftU
    MODULE PROCEDURE qes_reset_HubbardCommon
    MODULE PROCEDURE qes_reset_HubbardJ
    MODULE PROCEDURE qes_reset_starting_ns
    MODULE PROCEDURE qes_reset_Hubbard_ns
    MODULE PROCEDURE qes_reset_vdW
    MODULE PROCEDURE qes_reset_spin
    MODULE PROCEDURE qes_reset_bands
    MODULE PROCEDURE qes_reset_smearing
    MODULE PROCEDURE qes_reset_occupations
    MODULE PROCEDURE qes_reset_basis
    MODULE PROCEDURE qes_reset_basis_set
    MODULE PROCEDURE qes_reset_basisSetItem
    MODULE PROCEDURE qes_reset_reciprocal_lattice
    MODULE PROCEDURE qes_reset_electron_control
    MODULE PROCEDURE qes_reset_k_points_IBZ
    MODULE PROCEDURE qes_reset_monkhorst_pack
    MODULE PROCEDURE qes_reset_k_point
    MODULE PROCEDURE qes_reset_ion_control
    MODULE PROCEDURE qes_reset_bfgs
    MODULE PROCEDURE qes_reset_md
    MODULE PROCEDURE qes_reset_cell_control
    MODULE PROCEDURE qes_reset_symmetry_flags
    MODULE PROCEDURE qes_reset_boundary_conditions
    MODULE PROCEDURE qes_reset_esm
    MODULE PROCEDURE qes_reset_ekin_functional
    MODULE PROCEDURE qes_reset_spin_constraints
    MODULE PROCEDURE qes_reset_electric_field
    MODULE PROCEDURE qes_reset_gate_settings
    MODULE PROCEDURE qes_reset_atomic_constraints
    MODULE PROCEDURE qes_reset_atomic_constraint
    MODULE PROCEDURE qes_reset_inputOccupations
    MODULE PROCEDURE qes_reset_outputElectricField
    MODULE PROCEDURE qes_reset_BerryPhaseOutput
    MODULE PROCEDURE qes_reset_dipoleOutput
    MODULE PROCEDURE qes_reset_finiteFieldOut
    MODULE PROCEDURE qes_reset_polarization
    MODULE PROCEDURE qes_reset_ionicPolarization
    MODULE PROCEDURE qes_reset_electronicPolarization
    MODULE PROCEDURE qes_reset_phase
    MODULE PROCEDURE qes_reset_gateInfo
    MODULE PROCEDURE qes_reset_convergence_info
    MODULE PROCEDURE qes_reset_scf_conv
    MODULE PROCEDURE qes_reset_opt_conv
    MODULE PROCEDURE qes_reset_algorithmic_info
    MODULE PROCEDURE qes_reset_symmetries
    MODULE PROCEDURE qes_reset_symmetry
    MODULE PROCEDURE qes_reset_equivalent_atoms
    MODULE PROCEDURE qes_reset_info
    MODULE PROCEDURE qes_reset_outputPBC
    MODULE PROCEDURE qes_reset_magnetization
    MODULE PROCEDURE qes_reset_total_energy
    MODULE PROCEDURE qes_reset_band_structure
    MODULE PROCEDURE qes_reset_ks_energies
    MODULE PROCEDURE qes_reset_closed
    MODULE PROCEDURE qes_reset_vector
    MODULE PROCEDURE qes_reset_integerVector
    MODULE PROCEDURE qes_reset_matrix
    MODULE PROCEDURE qes_reset_integerMatrix
    MODULE PROCEDURE qes_reset_scalarQuantity
  END INTERFACE qes_reset
  !
  CONTAINS
  !
  !
  SUBROUTINE qes_reset_espresso(obj)
    !
    IMPLICIT NONE
    TYPE(espresso_type),INTENT(INOUT)    :: obj
    INTEGER :: i
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (obj%general_info_ispresent) &
      CALL qes_reset_general_info(obj%general_info)
    obj%general_info_ispresent = .FALSE.
    IF (obj%parallel_info_ispresent) &
      CALL qes_reset_parallel_info(obj%parallel_info)
    obj%parallel_info_ispresent = .FALSE.
    CALL qes_reset_input(obj%input)
    IF (obj%step_ispresent) THEN
      IF (ALLOCATED(obj%step)) THEN
        DO i=1, SIZE(obj%step)
          CALL qes_reset_step(obj%step(i))
        ENDDO
        DEALLOCATE(obj%step)
      ENDIF
      obj%ndim_step = 0
      obj%step_ispresent = .FALSE.
    ENDIF
    IF (obj%output_ispresent) &
      CALL qes_reset_output(obj%output)
    obj%output_ispresent = .FALSE.
    obj%status_ispresent = .FALSE.
    obj%cputime_ispresent = .FALSE.
    IF (obj%timing_info_ispresent) &
      CALL qes_reset_timing(obj%timing_info)
    obj%timing_info_ispresent = .FALSE.
    IF (obj%closed_ispresent) &
      CALL qes_reset_closed(obj%closed)
    obj%closed_ispresent = .FALSE.
    obj%Units_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_espresso
  !
  !
  SUBROUTINE qes_reset_general_info(obj)
    !
    IMPLICIT NONE
    TYPE(general_info_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    CALL qes_reset_xml_format(obj%xml_format)
    CALL qes_reset_creator(obj%creator)
    CALL qes_reset_created(obj%created)
    !
  END SUBROUTINE qes_reset_general_info
  !
  !
  SUBROUTINE qes_reset_parallel_info(obj)
    !
    IMPLICIT NONE
    TYPE(parallel_info_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_parallel_info
  !
  !
  SUBROUTINE qes_reset_input(obj)
    !
    IMPLICIT NONE
    TYPE(input_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    CALL qes_reset_control_variables(obj%control_variables)
    CALL qes_reset_atomic_species(obj%atomic_species)
    CALL qes_reset_atomic_structure(obj%atomic_structure)
    CALL qes_reset_dft(obj%dft)
    CALL qes_reset_spin(obj%spin)
    CALL qes_reset_bands(obj%bands)
    CALL qes_reset_basis(obj%basis)
    CALL qes_reset_electron_control(obj%electron_control)
    CALL qes_reset_k_points_IBZ(obj%k_points_IBZ)
    CALL qes_reset_ion_control(obj%ion_control)
    CALL qes_reset_cell_control(obj%cell_control)
    IF (obj%symmetry_flags_ispresent) &
      CALL qes_reset_symmetry_flags(obj%symmetry_flags)
    obj%symmetry_flags_ispresent = .FALSE.
    IF (obj%boundary_conditions_ispresent) &
      CALL qes_reset_boundary_conditions(obj%boundary_conditions)
    obj%boundary_conditions_ispresent = .FALSE.
    IF (obj%ekin_functional_ispresent) &
      CALL qes_reset_ekin_functional(obj%ekin_functional)
    obj%ekin_functional_ispresent = .FALSE.
    IF (obj%external_atomic_forces_ispresent) &
      CALL qes_reset_matrix(obj%external_atomic_forces)
    obj%external_atomic_forces_ispresent = .FALSE.
    IF (obj%free_positions_ispresent) &
      CALL qes_reset_integerMatrix(obj%free_positions)
    obj%free_positions_ispresent = .FALSE.
    IF (obj%starting_atomic_velocities_ispresent) &
      CALL qes_reset_matrix(obj%starting_atomic_velocities)
    obj%starting_atomic_velocities_ispresent = .FALSE.
    IF (obj%electric_field_ispresent) &
      CALL qes_reset_electric_field(obj%electric_field)
    obj%electric_field_ispresent = .FALSE.
    IF (obj%atomic_constraints_ispresent) &
      CALL qes_reset_atomic_constraints(obj%atomic_constraints)
    obj%atomic_constraints_ispresent = .FALSE.
    IF (obj%spin_constraints_ispresent) &
      CALL qes_reset_spin_constraints(obj%spin_constraints)
    obj%spin_constraints_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_input
  !
  !
  SUBROUTINE qes_reset_step(obj)
    !
    IMPLICIT NONE
    TYPE(step_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    CALL qes_reset_scf_conv(obj%scf_conv)
    CALL qes_reset_atomic_structure(obj%atomic_structure)
    CALL qes_reset_total_energy(obj%total_energy)
    CALL qes_reset_matrix(obj%forces)
    IF (obj%stress_ispresent) &
      CALL qes_reset_matrix(obj%stress)
    obj%stress_ispresent = .FALSE.
    obj%FCP_force_ispresent = .FALSE.
    obj%FCP_tot_charge_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_step
  !
  !
  SUBROUTINE qes_reset_output(obj)
    !
    IMPLICIT NONE
    TYPE(output_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (obj%convergence_info_ispresent) &
      CALL qes_reset_convergence_info(obj%convergence_info)
    obj%convergence_info_ispresent = .FALSE.
    CALL qes_reset_algorithmic_info(obj%algorithmic_info)
    CALL qes_reset_atomic_species(obj%atomic_species)
    CALL qes_reset_atomic_structure(obj%atomic_structure)
    IF (obj%symmetries_ispresent) &
      CALL qes_reset_symmetries(obj%symmetries)
    obj%symmetries_ispresent = .FALSE.
    CALL qes_reset_basis_set(obj%basis_set)
    CALL qes_reset_dft(obj%dft)
    IF (obj%boundary_conditions_ispresent) &
      CALL qes_reset_outputPBC(obj%boundary_conditions)
    obj%boundary_conditions_ispresent = .FALSE.
    CALL qes_reset_magnetization(obj%magnetization)
    CALL qes_reset_total_energy(obj%total_energy)
    CALL qes_reset_band_structure(obj%band_structure)
    IF (obj%forces_ispresent) &
      CALL qes_reset_matrix(obj%forces)
    obj%forces_ispresent = .FALSE.
    IF (obj%stress_ispresent) &
      CALL qes_reset_matrix(obj%stress)
    obj%stress_ispresent = .FALSE.
    IF (obj%electric_field_ispresent) &
      CALL qes_reset_outputElectricField(obj%electric_field)
    obj%electric_field_ispresent = .FALSE.
    obj%FCP_force_ispresent = .FALSE.
    obj%FCP_tot_charge_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_output
  !
  !
  SUBROUTINE qes_reset_timing(obj)
    !
    IMPLICIT NONE
    TYPE(timing_type),INTENT(INOUT)    :: obj
    INTEGER :: i
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    CALL qes_reset_clock(obj%total)
    IF (obj%partial_ispresent) THEN
      IF (ALLOCATED(obj%partial)) THEN
        DO i=1, SIZE(obj%partial)
          CALL qes_reset_clock(obj%partial(i))
        ENDDO
        DEALLOCATE(obj%partial)
      ENDIF
      obj%ndim_partial = 0
      obj%partial_ispresent = .FALSE.
    ENDIF
    !
  END SUBROUTINE qes_reset_timing
  !
  !
  SUBROUTINE qes_reset_clock(obj)
    !
    IMPLICIT NONE
    TYPE(clock_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%calls_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_clock
  !
  !
  SUBROUTINE qes_reset_control_variables(obj)
    !
    IMPLICIT NONE
    TYPE(control_variables_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%nstep_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_control_variables
  !
  !
  SUBROUTINE qes_reset_xml_format(obj)
    !
    IMPLICIT NONE
    TYPE(xml_format_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_xml_format
  !
  !
  SUBROUTINE qes_reset_creator(obj)
    !
    IMPLICIT NONE
    TYPE(creator_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_creator
  !
  !
  SUBROUTINE qes_reset_created(obj)
    !
    IMPLICIT NONE
    TYPE(created_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_created
  !
  !
  SUBROUTINE qes_reset_atomic_species(obj)
    !
    IMPLICIT NONE
    TYPE(atomic_species_type),INTENT(INOUT)    :: obj
    INTEGER :: i
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (ALLOCATED(obj%species)) THEN
      DO i=1, SIZE(obj%species)
        CALL qes_reset_species(obj%species(i))
      ENDDO
      DEALLOCATE(obj%species)
    ENDIF
    obj%ndim_species = 0
    obj%pseudo_dir_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_atomic_species
  !
  !
  SUBROUTINE qes_reset_species(obj)
    !
    IMPLICIT NONE
    TYPE(species_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%mass_ispresent = .FALSE.
    obj%starting_magnetization_ispresent = .FALSE.
    obj%spin_teta_ispresent = .FALSE.
    obj%spin_phi_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_species
  !
  !
  SUBROUTINE qes_reset_atomic_structure(obj)
    !
    IMPLICIT NONE
    TYPE(atomic_structure_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (obj%atomic_positions_ispresent) &
      CALL qes_reset_atomic_positions(obj%atomic_positions)
    obj%atomic_positions_ispresent = .FALSE.
    IF (obj%wyckoff_positions_ispresent) &
      CALL qes_reset_wyckoff_positions(obj%wyckoff_positions)
    obj%wyckoff_positions_ispresent = .FALSE.
    IF (obj%crystal_positions_ispresent) &
      CALL qes_reset_atomic_positions(obj%crystal_positions)
    obj%crystal_positions_ispresent = .FALSE.
    CALL qes_reset_cell(obj%cell)
    obj%alat_ispresent = .FALSE.
    obj%bravais_index_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_atomic_structure
  !
  !
  SUBROUTINE qes_reset_atomic_positions(obj)
    !
    IMPLICIT NONE
    TYPE(atomic_positions_type),INTENT(INOUT)    :: obj
    INTEGER :: i
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (ALLOCATED(obj%atom)) THEN
      DO i=1, SIZE(obj%atom)
        CALL qes_reset_atom(obj%atom(i))
      ENDDO
      DEALLOCATE(obj%atom)
    ENDIF
    obj%ndim_atom = 0
    !
  END SUBROUTINE qes_reset_atomic_positions
  !
  !
  SUBROUTINE qes_reset_atom(obj)
    !
    IMPLICIT NONE
    TYPE(atom_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%position_ispresent = .FALSE.
    obj%index_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_atom
  !
  !
  SUBROUTINE qes_reset_wyckoff_positions(obj)
    !
    IMPLICIT NONE
    TYPE(wyckoff_positions_type),INTENT(INOUT)    :: obj
    INTEGER :: i
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (ALLOCATED(obj%atom)) THEN
      DO i=1, SIZE(obj%atom)
        CALL qes_reset_atom(obj%atom(i))
      ENDDO
      DEALLOCATE(obj%atom)
    ENDIF
    obj%ndim_atom = 0
    obj%more_options_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_wyckoff_positions
  !
  !
  SUBROUTINE qes_reset_cell(obj)
    !
    IMPLICIT NONE
    TYPE(cell_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_cell
  !
  !
  SUBROUTINE qes_reset_dft(obj)
    !
    IMPLICIT NONE
    TYPE(dft_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (obj%hybrid_ispresent) &
      CALL qes_reset_hybrid(obj%hybrid)
    obj%hybrid_ispresent = .FALSE.
    IF (obj%dftU_ispresent) &
      CALL qes_reset_dftU(obj%dftU)
    obj%dftU_ispresent = .FALSE.
    IF (obj%vdW_ispresent) &
      CALL qes_reset_vdW(obj%vdW)
    obj%vdW_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_dft
  !
  !
  SUBROUTINE qes_reset_hybrid(obj)
    !
    IMPLICIT NONE
    TYPE(hybrid_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (obj%qpoint_grid_ispresent) &
      CALL qes_reset_qpoint_grid(obj%qpoint_grid)
    obj%qpoint_grid_ispresent = .FALSE.
    obj%ecutfock_ispresent = .FALSE.
    obj%exx_fraction_ispresent = .FALSE.
    obj%screening_parameter_ispresent = .FALSE.
    obj%exxdiv_treatment_ispresent = .FALSE.
    obj%x_gamma_extrapolation_ispresent = .FALSE.
    obj%ecutvcut_ispresent = .FALSE.
    obj%localization_threshold_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_hybrid
  !
  !
  SUBROUTINE qes_reset_qpoint_grid(obj)
    !
    IMPLICIT NONE
    TYPE(qpoint_grid_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_qpoint_grid
  !
  !
  SUBROUTINE qes_reset_dftU(obj)
    !
    IMPLICIT NONE
    TYPE(dftU_type),INTENT(INOUT)    :: obj
    INTEGER :: i
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%lda_plus_u_kind_ispresent = .FALSE.
    IF (obj%Hubbard_U_ispresent) THEN
      IF (ALLOCATED(obj%Hubbard_U)) THEN
        DO i=1, SIZE(obj%Hubbard_U)
          CALL qes_reset_HubbardCommon(obj%Hubbard_U(i))
        ENDDO
        DEALLOCATE(obj%Hubbard_U)
      ENDIF
      obj%ndim_Hubbard_U = 0
      obj%Hubbard_U_ispresent = .FALSE.
    ENDIF
    IF (obj%Hubbard_J0_ispresent) THEN
      IF (ALLOCATED(obj%Hubbard_J0)) THEN
        DO i=1, SIZE(obj%Hubbard_J0)
          CALL qes_reset_HubbardCommon(obj%Hubbard_J0(i))
        ENDDO
        DEALLOCATE(obj%Hubbard_J0)
      ENDIF
      obj%ndim_Hubbard_J0 = 0
      obj%Hubbard_J0_ispresent = .FALSE.
    ENDIF
    IF (obj%Hubbard_alpha_ispresent) THEN
      IF (ALLOCATED(obj%Hubbard_alpha)) THEN
        DO i=1, SIZE(obj%Hubbard_alpha)
          CALL qes_reset_HubbardCommon(obj%Hubbard_alpha(i))
        ENDDO
        DEALLOCATE(obj%Hubbard_alpha)
      ENDIF
      obj%ndim_Hubbard_alpha = 0
      obj%Hubbard_alpha_ispresent = .FALSE.
    ENDIF
    IF (obj%Hubbard_beta_ispresent) THEN
      IF (ALLOCATED(obj%Hubbard_beta)) THEN
        DO i=1, SIZE(obj%Hubbard_beta)
          CALL qes_reset_HubbardCommon(obj%Hubbard_beta(i))
        ENDDO
        DEALLOCATE(obj%Hubbard_beta)
      ENDIF
      obj%ndim_Hubbard_beta = 0
      obj%Hubbard_beta_ispresent = .FALSE.
    ENDIF
    IF (obj%Hubbard_J_ispresent) THEN
      IF (ALLOCATED(obj%Hubbard_J)) THEN
        DO i=1, SIZE(obj%Hubbard_J)
          CALL qes_reset_HubbardJ(obj%Hubbard_J(i))
        ENDDO
        DEALLOCATE(obj%Hubbard_J)
      ENDIF
      obj%ndim_Hubbard_J = 0
      obj%Hubbard_J_ispresent = .FALSE.
    ENDIF
    IF (obj%starting_ns_ispresent) THEN
      IF (ALLOCATED(obj%starting_ns)) THEN
        DO i=1, SIZE(obj%starting_ns)
          CALL qes_reset_starting_ns(obj%starting_ns(i))
        ENDDO
        DEALLOCATE(obj%starting_ns)
      ENDIF
      obj%ndim_starting_ns = 0
      obj%starting_ns_ispresent = .FALSE.
    ENDIF
    IF (obj%Hubbard_ns_ispresent) THEN
      IF (ALLOCATED(obj%Hubbard_ns)) THEN
        DO i=1, SIZE(obj%Hubbard_ns)
          CALL qes_reset_Hubbard_ns(obj%Hubbard_ns(i))
        ENDDO
        DEALLOCATE(obj%Hubbard_ns)
      ENDIF
      obj%ndim_Hubbard_ns = 0
      obj%Hubbard_ns_ispresent = .FALSE.
    ENDIF
    obj%U_projection_type_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_dftU
  !
  !
  SUBROUTINE qes_reset_HubbardCommon(obj)
    !
    IMPLICIT NONE
    TYPE(HubbardCommon_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%label_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_HubbardCommon
  !
  !
  SUBROUTINE qes_reset_HubbardJ(obj)
    !
    IMPLICIT NONE
    TYPE(HubbardJ_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_HubbardJ
  !
  !
  SUBROUTINE qes_reset_starting_ns(obj)
    !
    IMPLICIT NONE
    TYPE(starting_ns_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%size = 0
    IF (ALLOCATED(obj%starting_ns)) THEN
      DEALLOCATE(obj%starting_ns)
    ENDIF
    !
  END SUBROUTINE qes_reset_starting_ns
  !
  !
  SUBROUTINE qes_reset_Hubbard_ns(obj)
    !
    IMPLICIT NONE
    TYPE(Hubbard_ns_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (ALLOCATED(obj%dims)) THEN
      DEALLOCATE(obj%dims)
    ENDIF
    obj%rank = 0
    obj%order = 'F'
    IF (ALLOCATED(obj%Hubbard_ns)) THEN
      DEALLOCATE(obj%Hubbard_ns)
    ENDIF
    !
  END SUBROUTINE qes_reset_Hubbard_ns
  !
  !
  SUBROUTINE qes_reset_vdW(obj)
    !
    IMPLICIT NONE
    TYPE(vdW_type),INTENT(INOUT)    :: obj
    INTEGER :: i
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%vdw_corr_ispresent = .FALSE.
    obj%dftd3_version_ispresent = .FALSE.
    obj%dftd3_threebody_ispresent = .FALSE.
    obj%non_local_term_ispresent = .FALSE.
    obj%functional_ispresent = .FALSE.
    obj%total_energy_term_ispresent = .FALSE.
    obj%london_s6_ispresent = .FALSE.
    obj%ts_vdw_econv_thr_ispresent = .FALSE.
    obj%ts_vdw_isolated_ispresent = .FALSE.
    obj%london_rcut_ispresent = .FALSE.
    obj%xdm_a1_ispresent = .FALSE.
    obj%xdm_a2_ispresent = .FALSE.
    IF (obj%london_c6_ispresent) THEN
      IF (ALLOCATED(obj%london_c6)) THEN
        DO i=1, SIZE(obj%london_c6)
          CALL qes_reset_HubbardCommon(obj%london_c6(i))
        ENDDO
        DEALLOCATE(obj%london_c6)
      ENDIF
      obj%ndim_london_c6 = 0
      obj%london_c6_ispresent = .FALSE.
    ENDIF
    !
  END SUBROUTINE qes_reset_vdW
  !
  !
  SUBROUTINE qes_reset_spin(obj)
    !
    IMPLICIT NONE
    TYPE(spin_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_spin
  !
  !
  SUBROUTINE qes_reset_bands(obj)
    !
    IMPLICIT NONE
    TYPE(bands_type),INTENT(INOUT)    :: obj
    INTEGER :: i
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%nbnd_ispresent = .FALSE.
    IF (obj%smearing_ispresent) &
      CALL qes_reset_smearing(obj%smearing)
    obj%smearing_ispresent = .FALSE.
    obj%tot_charge_ispresent = .FALSE.
    obj%tot_magnetization_ispresent = .FALSE.
    CALL qes_reset_occupations(obj%occupations)
    IF (obj%inputOccupations_ispresent) THEN
      IF (ALLOCATED(obj%inputOccupations)) THEN
        DO i=1, SIZE(obj%inputOccupations)
          CALL qes_reset_inputOccupations(obj%inputOccupations(i))
        ENDDO
        DEALLOCATE(obj%inputOccupations)
      ENDIF
      obj%ndim_inputOccupations = 0
      obj%inputOccupations_ispresent = .FALSE.
    ENDIF
    !
  END SUBROUTINE qes_reset_bands
  !
  !
  SUBROUTINE qes_reset_smearing(obj)
    !
    IMPLICIT NONE
    TYPE(smearing_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_smearing
  !
  !
  SUBROUTINE qes_reset_occupations(obj)
    !
    IMPLICIT NONE
    TYPE(occupations_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%spin_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_occupations
  !
  !
  SUBROUTINE qes_reset_basis(obj)
    !
    IMPLICIT NONE
    TYPE(basis_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%gamma_only_ispresent = .FALSE.
    obj%ecutrho_ispresent = .FALSE.
    IF (obj%fft_grid_ispresent) &
      CALL qes_reset_basisSetItem(obj%fft_grid)
    obj%fft_grid_ispresent = .FALSE.
    IF (obj%fft_smooth_ispresent) &
      CALL qes_reset_basisSetItem(obj%fft_smooth)
    obj%fft_smooth_ispresent = .FALSE.
    IF (obj%fft_box_ispresent) &
      CALL qes_reset_basisSetItem(obj%fft_box)
    obj%fft_box_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_basis
  !
  !
  SUBROUTINE qes_reset_basis_set(obj)
    !
    IMPLICIT NONE
    TYPE(basis_set_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%gamma_only_ispresent = .FALSE.
    obj%ecutrho_ispresent = .FALSE.
    CALL qes_reset_basisSetItem(obj%fft_grid)
    IF (obj%fft_smooth_ispresent) &
      CALL qes_reset_basisSetItem(obj%fft_smooth)
    obj%fft_smooth_ispresent = .FALSE.
    IF (obj%fft_box_ispresent) &
      CALL qes_reset_basisSetItem(obj%fft_box)
    obj%fft_box_ispresent = .FALSE.
    obj%ngms_ispresent = .FALSE.
    CALL qes_reset_reciprocal_lattice(obj%reciprocal_lattice)
    !
  END SUBROUTINE qes_reset_basis_set
  !
  !
  SUBROUTINE qes_reset_basisSetItem(obj)
    !
    IMPLICIT NONE
    TYPE(basisSetItem_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_basisSetItem
  !
  !
  SUBROUTINE qes_reset_reciprocal_lattice(obj)
    !
    IMPLICIT NONE
    TYPE(reciprocal_lattice_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_reciprocal_lattice
  !
  !
  SUBROUTINE qes_reset_electron_control(obj)
    !
    IMPLICIT NONE
    TYPE(electron_control_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%real_space_q_ispresent = .FALSE.
    obj%real_space_beta_ispresent = .FALSE.
    obj%diago_cg_maxiter_ispresent = .FALSE.
    obj%diago_ppcg_maxiter_ispresent = .FALSE.
    obj%diago_david_ndim_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_electron_control
  !
  !
  SUBROUTINE qes_reset_k_points_IBZ(obj)
    !
    IMPLICIT NONE
    TYPE(k_points_IBZ_type),INTENT(INOUT)    :: obj
    INTEGER :: i
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (obj%monkhorst_pack_ispresent) &
      CALL qes_reset_monkhorst_pack(obj%monkhorst_pack)
    obj%monkhorst_pack_ispresent = .FALSE.
    obj%nk_ispresent = .FALSE.
    IF (obj%k_point_ispresent) THEN
      IF (ALLOCATED(obj%k_point)) THEN
        DO i=1, SIZE(obj%k_point)
          CALL qes_reset_k_point(obj%k_point(i))
        ENDDO
        DEALLOCATE(obj%k_point)
      ENDIF
      obj%ndim_k_point = 0
      obj%k_point_ispresent = .FALSE.
    ENDIF
    !
  END SUBROUTINE qes_reset_k_points_IBZ
  !
  !
  SUBROUTINE qes_reset_monkhorst_pack(obj)
    !
    IMPLICIT NONE
    TYPE(monkhorst_pack_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_monkhorst_pack
  !
  !
  SUBROUTINE qes_reset_k_point(obj)
    !
    IMPLICIT NONE
    TYPE(k_point_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%weight_ispresent = .FALSE.
    obj%label_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_k_point
  !
  !
  SUBROUTINE qes_reset_ion_control(obj)
    !
    IMPLICIT NONE
    TYPE(ion_control_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%upscale_ispresent = .FALSE.
    obj%remove_rigid_rot_ispresent = .FALSE.
    obj%refold_pos_ispresent = .FALSE.
    IF (obj%bfgs_ispresent) &
      CALL qes_reset_bfgs(obj%bfgs)
    obj%bfgs_ispresent = .FALSE.
    IF (obj%md_ispresent) &
      CALL qes_reset_md(obj%md)
    obj%md_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_ion_control
  !
  !
  SUBROUTINE qes_reset_bfgs(obj)
    !
    IMPLICIT NONE
    TYPE(bfgs_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_bfgs
  !
  !
  SUBROUTINE qes_reset_md(obj)
    !
    IMPLICIT NONE
    TYPE(md_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_md
  !
  !
  SUBROUTINE qes_reset_cell_control(obj)
    !
    IMPLICIT NONE
    TYPE(cell_control_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%wmass_ispresent = .FALSE.
    obj%cell_factor_ispresent = .FALSE.
    obj%fix_volume_ispresent = .FALSE.
    obj%fix_area_ispresent = .FALSE.
    obj%isotropic_ispresent = .FALSE.
    IF (obj%free_cell_ispresent) &
      CALL qes_reset_integerMatrix(obj%free_cell)
    obj%free_cell_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_cell_control
  !
  !
  SUBROUTINE qes_reset_symmetry_flags(obj)
    !
    IMPLICIT NONE
    TYPE(symmetry_flags_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_symmetry_flags
  !
  !
  SUBROUTINE qes_reset_boundary_conditions(obj)
    !
    IMPLICIT NONE
    TYPE(boundary_conditions_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (obj%esm_ispresent) &
      CALL qes_reset_esm(obj%esm)
    obj%esm_ispresent = .FALSE.
    obj%fcp_opt_ispresent = .FALSE.
    obj%fcp_mu_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_boundary_conditions
  !
  !
  SUBROUTINE qes_reset_esm(obj)
    !
    IMPLICIT NONE
    TYPE(esm_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_esm
  !
  !
  SUBROUTINE qes_reset_ekin_functional(obj)
    !
    IMPLICIT NONE
    TYPE(ekin_functional_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_ekin_functional
  !
  !
  SUBROUTINE qes_reset_spin_constraints(obj)
    !
    IMPLICIT NONE
    TYPE(spin_constraints_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%target_magnetization_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_spin_constraints
  !
  !
  SUBROUTINE qes_reset_electric_field(obj)
    !
    IMPLICIT NONE
    TYPE(electric_field_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%dipole_correction_ispresent = .FALSE.
    IF (obj%gate_settings_ispresent) &
      CALL qes_reset_gate_settings(obj%gate_settings)
    obj%gate_settings_ispresent = .FALSE.
    obj%electric_field_direction_ispresent = .FALSE.
    obj%potential_max_position_ispresent = .FALSE.
    obj%potential_decrease_width_ispresent = .FALSE.
    obj%electric_field_amplitude_ispresent = .FALSE.
    obj%electric_field_vector_ispresent = .FALSE.
    obj%nk_per_string_ispresent = .FALSE.
    obj%n_berry_cycles_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_electric_field
  !
  !
  SUBROUTINE qes_reset_gate_settings(obj)
    !
    IMPLICIT NONE
    TYPE(gate_settings_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%zgate_ispresent = .FALSE.
    obj%relaxz_ispresent = .FALSE.
    obj%block_ispresent = .FALSE.
    obj%block_1_ispresent = .FALSE.
    obj%block_2_ispresent = .FALSE.
    obj%block_height_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_gate_settings
  !
  !
  SUBROUTINE qes_reset_atomic_constraints(obj)
    !
    IMPLICIT NONE
    TYPE(atomic_constraints_type),INTENT(INOUT)    :: obj
    INTEGER :: i
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (ALLOCATED(obj%atomic_constraint)) THEN
      DO i=1, SIZE(obj%atomic_constraint)
        CALL qes_reset_atomic_constraint(obj%atomic_constraint(i))
      ENDDO
      DEALLOCATE(obj%atomic_constraint)
    ENDIF
    obj%ndim_atomic_constraint = 0
    !
  END SUBROUTINE qes_reset_atomic_constraints
  !
  !
  SUBROUTINE qes_reset_atomic_constraint(obj)
    !
    IMPLICIT NONE
    TYPE(atomic_constraint_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_atomic_constraint
  !
  !
  SUBROUTINE qes_reset_inputOccupations(obj)
    !
    IMPLICIT NONE
    TYPE(inputOccupations_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%size = 0
    IF (ALLOCATED(obj%inputOccupations)) THEN
      DEALLOCATE(obj%inputOccupations)
    ENDIF
    !
  END SUBROUTINE qes_reset_inputOccupations
  !
  !
  SUBROUTINE qes_reset_outputElectricField(obj)
    !
    IMPLICIT NONE
    TYPE(outputElectricField_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (obj%BerryPhase_ispresent) &
      CALL qes_reset_BerryPhaseOutput(obj%BerryPhase)
    obj%BerryPhase_ispresent = .FALSE.
    IF (obj%finiteElectricFieldInfo_ispresent) &
      CALL qes_reset_finiteFieldOut(obj%finiteElectricFieldInfo)
    obj%finiteElectricFieldInfo_ispresent = .FALSE.
    IF (obj%dipoleInfo_ispresent) &
      CALL qes_reset_dipoleOutput(obj%dipoleInfo)
    obj%dipoleInfo_ispresent = .FALSE.
    IF (obj%gateInfo_ispresent) &
      CALL qes_reset_gateInfo(obj%gateInfo)
    obj%gateInfo_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_outputElectricField
  !
  !
  SUBROUTINE qes_reset_BerryPhaseOutput(obj)
    !
    IMPLICIT NONE
    TYPE(BerryPhaseOutput_type),INTENT(INOUT)    :: obj
    INTEGER :: i
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    CALL qes_reset_polarization(obj%totalPolarization)
    CALL qes_reset_phase(obj%totalPhase)
    IF (ALLOCATED(obj%ionicPolarization)) THEN
      DO i=1, SIZE(obj%ionicPolarization)
        CALL qes_reset_ionicPolarization(obj%ionicPolarization(i))
      ENDDO
      DEALLOCATE(obj%ionicPolarization)
    ENDIF
    obj%ndim_ionicPolarization = 0
    IF (ALLOCATED(obj%electronicPolarization)) THEN
      DO i=1, SIZE(obj%electronicPolarization)
        CALL qes_reset_electronicPolarization(obj%electronicPolarization(i))
      ENDDO
      DEALLOCATE(obj%electronicPolarization)
    ENDIF
    obj%ndim_electronicPolarization = 0
    !
  END SUBROUTINE qes_reset_BerryPhaseOutput
  !
  !
  SUBROUTINE qes_reset_dipoleOutput(obj)
    !
    IMPLICIT NONE
    TYPE(dipoleOutput_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    CALL qes_reset_scalarQuantity(obj%dipole)
    CALL qes_reset_scalarQuantity(obj%ion_dipole)
    CALL qes_reset_scalarQuantity(obj%elec_dipole)
    CALL qes_reset_scalarQuantity(obj%dipoleField)
    CALL qes_reset_scalarQuantity(obj%potentialAmp)
    CALL qes_reset_scalarQuantity(obj%totalLength)
    !
  END SUBROUTINE qes_reset_dipoleOutput
  !
  !
  SUBROUTINE qes_reset_finiteFieldOut(obj)
    !
    IMPLICIT NONE
    TYPE(finiteFieldOut_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_finiteFieldOut
  !
  !
  SUBROUTINE qes_reset_polarization(obj)
    !
    IMPLICIT NONE
    TYPE(polarization_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    CALL qes_reset_scalarQuantity(obj%polarization)
    !
  END SUBROUTINE qes_reset_polarization
  !
  !
  SUBROUTINE qes_reset_ionicPolarization(obj)
    !
    IMPLICIT NONE
    TYPE(ionicPolarization_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    CALL qes_reset_atom(obj%ion)
    CALL qes_reset_phase(obj%phase)
    !
  END SUBROUTINE qes_reset_ionicPolarization
  !
  !
  SUBROUTINE qes_reset_electronicPolarization(obj)
    !
    IMPLICIT NONE
    TYPE(electronicPolarization_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    CALL qes_reset_k_point(obj%firstKeyPoint)
    obj%spin_ispresent = .FALSE.
    CALL qes_reset_phase(obj%phase)
    !
  END SUBROUTINE qes_reset_electronicPolarization
  !
  !
  SUBROUTINE qes_reset_phase(obj)
    !
    IMPLICIT NONE
    TYPE(phase_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%ionic_ispresent = .FALSE.
    obj%electronic_ispresent = .FALSE.
    obj%modulus_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_phase
  !
  !
  SUBROUTINE qes_reset_gateInfo(obj)
    !
    IMPLICIT NONE
    TYPE(gateInfo_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_gateInfo
  !
  !
  SUBROUTINE qes_reset_convergence_info(obj)
    !
    IMPLICIT NONE
    TYPE(convergence_info_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    CALL qes_reset_scf_conv(obj%scf_conv)
    IF (obj%opt_conv_ispresent) &
      CALL qes_reset_opt_conv(obj%opt_conv)
    obj%opt_conv_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_convergence_info
  !
  !
  SUBROUTINE qes_reset_scf_conv(obj)
    !
    IMPLICIT NONE
    TYPE(scf_conv_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_scf_conv
  !
  !
  SUBROUTINE qes_reset_opt_conv(obj)
    !
    IMPLICIT NONE
    TYPE(opt_conv_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_opt_conv
  !
  !
  SUBROUTINE qes_reset_algorithmic_info(obj)
    !
    IMPLICIT NONE
    TYPE(algorithmic_info_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%real_space_beta_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_algorithmic_info
  !
  !
  SUBROUTINE qes_reset_symmetries(obj)
    !
    IMPLICIT NONE
    TYPE(symmetries_type),INTENT(INOUT)    :: obj
    INTEGER :: i
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (ALLOCATED(obj%symmetry)) THEN
      DO i=1, SIZE(obj%symmetry)
        CALL qes_reset_symmetry(obj%symmetry(i))
      ENDDO
      DEALLOCATE(obj%symmetry)
    ENDIF
    obj%ndim_symmetry = 0
    !
  END SUBROUTINE qes_reset_symmetries
  !
  !
  SUBROUTINE qes_reset_symmetry(obj)
    !
    IMPLICIT NONE
    TYPE(symmetry_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    CALL qes_reset_info(obj%info)
    CALL qes_reset_matrix(obj%rotation)
    obj%fractional_translation_ispresent = .FALSE.
    IF (obj%equivalent_atoms_ispresent) &
      CALL qes_reset_equivalent_atoms(obj%equivalent_atoms)
    obj%equivalent_atoms_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_symmetry
  !
  !
  SUBROUTINE qes_reset_equivalent_atoms(obj)
    !
    IMPLICIT NONE
    TYPE(equivalent_atoms_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%size = 0
    IF (ALLOCATED(obj%equivalent_atoms)) THEN
      DEALLOCATE(obj%equivalent_atoms)
    ENDIF
    !
  END SUBROUTINE qes_reset_equivalent_atoms
  !
  !
  SUBROUTINE qes_reset_info(obj)
    !
    IMPLICIT NONE
    TYPE(info_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%name_ispresent = .FALSE.
    obj%class_ispresent = .FALSE.
    obj%time_reversal_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_info
  !
  !
  SUBROUTINE qes_reset_outputPBC(obj)
    !
    IMPLICIT NONE
    TYPE(outputPBC_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_outputPBC
  !
  !
  SUBROUTINE qes_reset_magnetization(obj)
    !
    IMPLICIT NONE
    TYPE(magnetization_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_magnetization
  !
  !
  SUBROUTINE qes_reset_total_energy(obj)
    !
    IMPLICIT NONE
    TYPE(total_energy_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%eband_ispresent = .FALSE.
    obj%ehart_ispresent = .FALSE.
    obj%vtxc_ispresent = .FALSE.
    obj%etxc_ispresent = .FALSE.
    obj%ewald_ispresent = .FALSE.
    obj%demet_ispresent = .FALSE.
    obj%efieldcorr_ispresent = .FALSE.
    obj%potentiostat_contr_ispresent = .FALSE.
    obj%gatefield_contr_ispresent = .FALSE.
    obj%vdW_term_ispresent = .FALSE.
    !
  END SUBROUTINE qes_reset_total_energy
  !
  !
  SUBROUTINE qes_reset_band_structure(obj)
    !
    IMPLICIT NONE
    TYPE(band_structure_type),INTENT(INOUT)    :: obj
    INTEGER :: i
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    obj%nbnd_ispresent = .FALSE.
    obj%nbnd_up_ispresent = .FALSE.
    obj%nbnd_dw_ispresent = .FALSE.
    obj%num_of_atomic_wfc_ispresent = .FALSE.
    obj%fermi_energy_ispresent = .FALSE.
    obj%highestOccupiedLevel_ispresent = .FALSE.
    obj%lowestUnoccupiedLevel_ispresent = .FALSE.
    obj%two_fermi_energies_ispresent = .FALSE.
    CALL qes_reset_k_points_IBZ(obj%starting_k_points)
    CALL qes_reset_occupations(obj%occupations_kind)
    IF (obj%smearing_ispresent) &
      CALL qes_reset_smearing(obj%smearing)
    obj%smearing_ispresent = .FALSE.
    IF (ALLOCATED(obj%ks_energies)) THEN
      DO i=1, SIZE(obj%ks_energies)
        CALL qes_reset_ks_energies(obj%ks_energies(i))
      ENDDO
      DEALLOCATE(obj%ks_energies)
    ENDIF
    obj%ndim_ks_energies = 0
    !
  END SUBROUTINE qes_reset_band_structure
  !
  !
  SUBROUTINE qes_reset_ks_energies(obj)
    !
    IMPLICIT NONE
    TYPE(ks_energies_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    CALL qes_reset_k_point(obj%k_point)
    CALL qes_reset_vector(obj%eigenvalues)
    CALL qes_reset_vector(obj%occupations)
    !
  END SUBROUTINE qes_reset_ks_energies
  !
  !
  SUBROUTINE qes_reset_closed(obj)
    !
    IMPLICIT NONE
    TYPE(closed_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_closed
  !
  !
  SUBROUTINE qes_reset_vector(obj)
    !
    IMPLICIT NONE
    TYPE(vector_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (ALLOCATED(obj%vector)) THEN
      DEALLOCATE(obj%vector)
    ENDIF
    obj%size = 0
    !
  END SUBROUTINE qes_reset_vector
  !
  !
  SUBROUTINE qes_reset_integerVector(obj)
    !
    IMPLICIT NONE
    TYPE(integerVector_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (ALLOCATED(obj%integerVector)) THEN
      DEALLOCATE(obj%integerVector)
    ENDIF
    obj%size = 0
    !
  END SUBROUTINE qes_reset_integerVector
  !
  !
  SUBROUTINE qes_reset_matrix(obj)
    !
    IMPLICIT NONE
    TYPE(matrix_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (ALLOCATED(obj%matrix)) THEN
      DEALLOCATE(obj%matrix)
    ENDIF
    IF (ALLOCATED(obj%dims)) THEN
      DEALLOCATE(obj%dims)
    ENDIF
    obj%rank = 0
    obj%order = 'F'
    !
  END SUBROUTINE qes_reset_matrix
  !
  !
  SUBROUTINE qes_reset_integerMatrix(obj)
    !
    IMPLICIT NONE
    TYPE(integerMatrix_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    IF (ALLOCATED(obj%integerMatrix)) THEN
      DEALLOCATE(obj%integerMatrix)
    ENDIF
    IF (ALLOCATED(obj%dims)) THEN
      DEALLOCATE(obj%dims)
    ENDIF
    obj%rank = 0
    obj%order = 'F'
    !
  END SUBROUTINE qes_reset_integerMatrix
  !
  !
  SUBROUTINE qes_reset_scalarQuantity(obj)
    !
    IMPLICIT NONE
    TYPE(scalarQuantity_type),INTENT(INOUT)    :: obj
    !
    obj%tagname = ""
    obj%lwrite  = .FALSE.
    obj%lread  = .FALSE.
    !
    !
  END SUBROUTINE qes_reset_scalarQuantity
  !
  !
END MODULE qes_reset_module