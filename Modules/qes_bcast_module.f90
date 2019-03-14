!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE qes_bcast_module
  !
  ! Auto-generated code: don't edit this file
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes-1.0
  !
  USE qes_types_module
  USE io_global, ONLY : ionode
  USE mp, ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  PUBLIC qes_bcast
  !
  INTERFACE qes_bcast
    MODULE PROCEDURE qes_bcast_espresso
    MODULE PROCEDURE qes_bcast_general_info
    MODULE PROCEDURE qes_bcast_parallel_info
    MODULE PROCEDURE qes_bcast_input
    MODULE PROCEDURE qes_bcast_step
    MODULE PROCEDURE qes_bcast_output
    MODULE PROCEDURE qes_bcast_timing
    MODULE PROCEDURE qes_bcast_clock
    MODULE PROCEDURE qes_bcast_control_variables
    MODULE PROCEDURE qes_bcast_xml_format
    MODULE PROCEDURE qes_bcast_creator
    MODULE PROCEDURE qes_bcast_created
    MODULE PROCEDURE qes_bcast_atomic_species
    MODULE PROCEDURE qes_bcast_species
    MODULE PROCEDURE qes_bcast_atomic_structure
    MODULE PROCEDURE qes_bcast_atomic_positions
    MODULE PROCEDURE qes_bcast_atom
    MODULE PROCEDURE qes_bcast_wyckoff_positions
    MODULE PROCEDURE qes_bcast_cell
    MODULE PROCEDURE qes_bcast_dft
    MODULE PROCEDURE qes_bcast_hybrid
    MODULE PROCEDURE qes_bcast_qpoint_grid
    MODULE PROCEDURE qes_bcast_dftU
    MODULE PROCEDURE qes_bcast_HubbardCommon
    MODULE PROCEDURE qes_bcast_HubbardJ
    MODULE PROCEDURE qes_bcast_starting_ns
    MODULE PROCEDURE qes_bcast_Hubbard_ns
    MODULE PROCEDURE qes_bcast_vdW
    MODULE PROCEDURE qes_bcast_spin
    MODULE PROCEDURE qes_bcast_bands
    MODULE PROCEDURE qes_bcast_smearing
    MODULE PROCEDURE qes_bcast_occupations
    MODULE PROCEDURE qes_bcast_basis
    MODULE PROCEDURE qes_bcast_basis_set
    MODULE PROCEDURE qes_bcast_basisSetItem
    MODULE PROCEDURE qes_bcast_reciprocal_lattice
    MODULE PROCEDURE qes_bcast_electron_control
    MODULE PROCEDURE qes_bcast_k_points_IBZ
    MODULE PROCEDURE qes_bcast_monkhorst_pack
    MODULE PROCEDURE qes_bcast_k_point
    MODULE PROCEDURE qes_bcast_ion_control
    MODULE PROCEDURE qes_bcast_bfgs
    MODULE PROCEDURE qes_bcast_md
    MODULE PROCEDURE qes_bcast_cell_control
    MODULE PROCEDURE qes_bcast_symmetry_flags
    MODULE PROCEDURE qes_bcast_boundary_conditions
    MODULE PROCEDURE qes_bcast_esm
    MODULE PROCEDURE qes_bcast_ekin_functional
    MODULE PROCEDURE qes_bcast_spin_constraints
    MODULE PROCEDURE qes_bcast_electric_field
    MODULE PROCEDURE qes_bcast_gate_settings
    MODULE PROCEDURE qes_bcast_atomic_constraints
    MODULE PROCEDURE qes_bcast_atomic_constraint
    MODULE PROCEDURE qes_bcast_inputOccupations
    MODULE PROCEDURE qes_bcast_outputElectricField
    MODULE PROCEDURE qes_bcast_BerryPhaseOutput
    MODULE PROCEDURE qes_bcast_dipoleOutput
    MODULE PROCEDURE qes_bcast_finiteFieldOut
    MODULE PROCEDURE qes_bcast_polarization
    MODULE PROCEDURE qes_bcast_ionicPolarization
    MODULE PROCEDURE qes_bcast_electronicPolarization
    MODULE PROCEDURE qes_bcast_phase
    MODULE PROCEDURE qes_bcast_gateInfo
    MODULE PROCEDURE qes_bcast_convergence_info
    MODULE PROCEDURE qes_bcast_scf_conv
    MODULE PROCEDURE qes_bcast_opt_conv
    MODULE PROCEDURE qes_bcast_algorithmic_info
    MODULE PROCEDURE qes_bcast_symmetries
    MODULE PROCEDURE qes_bcast_symmetry
    MODULE PROCEDURE qes_bcast_equivalent_atoms
    MODULE PROCEDURE qes_bcast_info
    MODULE PROCEDURE qes_bcast_outputPBC
    MODULE PROCEDURE qes_bcast_magnetization
    MODULE PROCEDURE qes_bcast_total_energy
    MODULE PROCEDURE qes_bcast_band_structure
    MODULE PROCEDURE qes_bcast_ks_energies
    MODULE PROCEDURE qes_bcast_closed
    MODULE PROCEDURE qes_bcast_vector
    MODULE PROCEDURE qes_bcast_integerVector
    MODULE PROCEDURE qes_bcast_matrix
    MODULE PROCEDURE qes_bcast_integerMatrix
    MODULE PROCEDURE qes_bcast_scalarQuantity
  END INTERFACE qes_bcast
  !
  CONTAINS
  !
  !
  SUBROUTINE qes_bcast_espresso(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(espresso_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%Units_ispresent, ionode_id, comm)
    IF (obj%Units_ispresent) &
      CALL mp_bcast(obj%Units, ionode_id, comm)
    CALL mp_bcast(obj%general_info_ispresent, ionode_id, comm)
    IF (obj%general_info_ispresent) &
      CALL qes_bcast_general_info(obj%general_info, ionode_id, comm)
    CALL mp_bcast(obj%parallel_info_ispresent, ionode_id, comm)
    IF (obj%parallel_info_ispresent) &
      CALL qes_bcast_parallel_info(obj%parallel_info, ionode_id, comm)
    CALL qes_bcast_input(obj%input, ionode_id, comm)
    CALL mp_bcast(obj%step_ispresent, ionode_id, comm)
    IF (obj%step_ispresent) THEN
      CALL mp_bcast(obj%ndim_step, ionode_id, comm)
      IF (.NOT.ionode) ALLOCATE(obj%step(obj%ndim_step))
      DO i=1, obj%ndim_step
        CALL qes_bcast_step(obj%step(i), ionode_id, comm)
      ENDDO
    ENDIF
    CALL mp_bcast(obj%output_ispresent, ionode_id, comm)
    IF (obj%output_ispresent) &
      CALL qes_bcast_output(obj%output, ionode_id, comm)
    CALL mp_bcast(obj%status_ispresent, ionode_id, comm)
    IF (obj%status_ispresent) &
      CALL mp_bcast(obj%status, ionode_id, comm)
    CALL mp_bcast(obj%cputime_ispresent, ionode_id, comm)
    IF (obj%cputime_ispresent) &
      CALL mp_bcast(obj%cputime, ionode_id, comm)
    CALL mp_bcast(obj%timing_info_ispresent, ionode_id, comm)
    IF (obj%timing_info_ispresent) &
      CALL qes_bcast_timing(obj%timing_info, ionode_id, comm)
    CALL mp_bcast(obj%closed_ispresent, ionode_id, comm)
    IF (obj%closed_ispresent) &
      CALL qes_bcast_closed(obj%closed, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_espresso
  !
  !
  SUBROUTINE qes_bcast_general_info(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(general_info_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL qes_bcast_xml_format(obj%xml_format, ionode_id, comm)
    CALL qes_bcast_creator(obj%creator, ionode_id, comm)
    CALL qes_bcast_created(obj%created, ionode_id, comm)
    CALL mp_bcast(obj%job, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_general_info
  !
  !
  SUBROUTINE qes_bcast_parallel_info(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(parallel_info_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%nprocs, ionode_id, comm)
    CALL mp_bcast(obj%nthreads, ionode_id, comm)
    CALL mp_bcast(obj%ntasks, ionode_id, comm)
    CALL mp_bcast(obj%nbgrp, ionode_id, comm)
    CALL mp_bcast(obj%npool, ionode_id, comm)
    CALL mp_bcast(obj%ndiag, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_parallel_info
  !
  !
  SUBROUTINE qes_bcast_input(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(input_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL qes_bcast_control_variables(obj%control_variables, ionode_id, comm)
    CALL qes_bcast_atomic_species(obj%atomic_species, ionode_id, comm)
    CALL qes_bcast_atomic_structure(obj%atomic_structure, ionode_id, comm)
    CALL qes_bcast_dft(obj%dft, ionode_id, comm)
    CALL qes_bcast_spin(obj%spin, ionode_id, comm)
    CALL qes_bcast_bands(obj%bands, ionode_id, comm)
    CALL qes_bcast_basis(obj%basis, ionode_id, comm)
    CALL qes_bcast_electron_control(obj%electron_control, ionode_id, comm)
    CALL qes_bcast_k_points_IBZ(obj%k_points_IBZ, ionode_id, comm)
    CALL qes_bcast_ion_control(obj%ion_control, ionode_id, comm)
    CALL qes_bcast_cell_control(obj%cell_control, ionode_id, comm)
    CALL mp_bcast(obj%symmetry_flags_ispresent, ionode_id, comm)
    IF (obj%symmetry_flags_ispresent) &
      CALL qes_bcast_symmetry_flags(obj%symmetry_flags, ionode_id, comm)
    CALL mp_bcast(obj%boundary_conditions_ispresent, ionode_id, comm)
    IF (obj%boundary_conditions_ispresent) &
      CALL qes_bcast_boundary_conditions(obj%boundary_conditions, ionode_id, comm)
    CALL mp_bcast(obj%ekin_functional_ispresent, ionode_id, comm)
    IF (obj%ekin_functional_ispresent) &
      CALL qes_bcast_ekin_functional(obj%ekin_functional, ionode_id, comm)
    CALL mp_bcast(obj%external_atomic_forces_ispresent, ionode_id, comm)
    IF (obj%external_atomic_forces_ispresent) &
      CALL qes_bcast_matrix(obj%external_atomic_forces, ionode_id, comm)
    CALL mp_bcast(obj%free_positions_ispresent, ionode_id, comm)
    IF (obj%free_positions_ispresent) &
      CALL qes_bcast_integerMatrix(obj%free_positions, ionode_id, comm)
    CALL mp_bcast(obj%starting_atomic_velocities_ispresent, ionode_id, comm)
    IF (obj%starting_atomic_velocities_ispresent) &
      CALL qes_bcast_matrix(obj%starting_atomic_velocities, ionode_id, comm)
    CALL mp_bcast(obj%electric_field_ispresent, ionode_id, comm)
    IF (obj%electric_field_ispresent) &
      CALL qes_bcast_electric_field(obj%electric_field, ionode_id, comm)
    CALL mp_bcast(obj%atomic_constraints_ispresent, ionode_id, comm)
    IF (obj%atomic_constraints_ispresent) &
      CALL qes_bcast_atomic_constraints(obj%atomic_constraints, ionode_id, comm)
    CALL mp_bcast(obj%spin_constraints_ispresent, ionode_id, comm)
    IF (obj%spin_constraints_ispresent) &
      CALL qes_bcast_spin_constraints(obj%spin_constraints, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_input
  !
  !
  SUBROUTINE qes_bcast_step(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(step_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%n_step, ionode_id, comm)
    CALL qes_bcast_scf_conv(obj%scf_conv, ionode_id, comm)
    CALL qes_bcast_atomic_structure(obj%atomic_structure, ionode_id, comm)
    CALL qes_bcast_total_energy(obj%total_energy, ionode_id, comm)
    CALL qes_bcast_matrix(obj%forces, ionode_id, comm)
    CALL mp_bcast(obj%stress_ispresent, ionode_id, comm)
    IF (obj%stress_ispresent) &
      CALL qes_bcast_matrix(obj%stress, ionode_id, comm)
    CALL mp_bcast(obj%FCP_force_ispresent, ionode_id, comm)
    IF (obj%FCP_force_ispresent) &
      CALL mp_bcast(obj%FCP_force, ionode_id, comm)
    CALL mp_bcast(obj%FCP_tot_charge_ispresent, ionode_id, comm)
    IF (obj%FCP_tot_charge_ispresent) &
      CALL mp_bcast(obj%FCP_tot_charge, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_step
  !
  !
  SUBROUTINE qes_bcast_output(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(output_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%convergence_info_ispresent, ionode_id, comm)
    IF (obj%convergence_info_ispresent) &
      CALL qes_bcast_convergence_info(obj%convergence_info, ionode_id, comm)
    CALL qes_bcast_algorithmic_info(obj%algorithmic_info, ionode_id, comm)
    CALL qes_bcast_atomic_species(obj%atomic_species, ionode_id, comm)
    CALL qes_bcast_atomic_structure(obj%atomic_structure, ionode_id, comm)
    CALL mp_bcast(obj%symmetries_ispresent, ionode_id, comm)
    IF (obj%symmetries_ispresent) &
      CALL qes_bcast_symmetries(obj%symmetries, ionode_id, comm)
    CALL qes_bcast_basis_set(obj%basis_set, ionode_id, comm)
    CALL qes_bcast_dft(obj%dft, ionode_id, comm)
    CALL mp_bcast(obj%boundary_conditions_ispresent, ionode_id, comm)
    IF (obj%boundary_conditions_ispresent) &
      CALL qes_bcast_outputPBC(obj%boundary_conditions, ionode_id, comm)
    CALL qes_bcast_magnetization(obj%magnetization, ionode_id, comm)
    CALL qes_bcast_total_energy(obj%total_energy, ionode_id, comm)
    CALL qes_bcast_band_structure(obj%band_structure, ionode_id, comm)
    CALL mp_bcast(obj%forces_ispresent, ionode_id, comm)
    IF (obj%forces_ispresent) &
      CALL qes_bcast_matrix(obj%forces, ionode_id, comm)
    CALL mp_bcast(obj%stress_ispresent, ionode_id, comm)
    IF (obj%stress_ispresent) &
      CALL qes_bcast_matrix(obj%stress, ionode_id, comm)
    CALL mp_bcast(obj%electric_field_ispresent, ionode_id, comm)
    IF (obj%electric_field_ispresent) &
      CALL qes_bcast_outputElectricField(obj%electric_field, ionode_id, comm)
    CALL mp_bcast(obj%FCP_force_ispresent, ionode_id, comm)
    IF (obj%FCP_force_ispresent) &
      CALL mp_bcast(obj%FCP_force, ionode_id, comm)
    CALL mp_bcast(obj%FCP_tot_charge_ispresent, ionode_id, comm)
    IF (obj%FCP_tot_charge_ispresent) &
      CALL mp_bcast(obj%FCP_tot_charge, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_output
  !
  !
  SUBROUTINE qes_bcast_timing(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(timing_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL qes_bcast_clock(obj%total, ionode_id, comm)
    CALL mp_bcast(obj%partial_ispresent, ionode_id, comm)
    IF (obj%partial_ispresent) THEN
      CALL mp_bcast(obj%ndim_partial, ionode_id, comm)
      IF (.NOT.ionode) ALLOCATE(obj%partial(obj%ndim_partial))
      DO i=1, obj%ndim_partial
        CALL qes_bcast_clock(obj%partial(i), ionode_id, comm)
      ENDDO
    ENDIF
    !
  END SUBROUTINE qes_bcast_timing
  !
  !
  SUBROUTINE qes_bcast_clock(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(clock_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%label, ionode_id, comm)
    CALL mp_bcast(obj%calls_ispresent, ionode_id, comm)
    IF (obj%calls_ispresent) &
      CALL mp_bcast(obj%calls, ionode_id, comm)
    CALL mp_bcast(obj%cpu, ionode_id, comm)
    CALL mp_bcast(obj%wall, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_clock
  !
  !
  SUBROUTINE qes_bcast_control_variables(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(control_variables_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%title, ionode_id, comm)
    CALL mp_bcast(obj%calculation, ionode_id, comm)
    CALL mp_bcast(obj%restart_mode, ionode_id, comm)
    CALL mp_bcast(obj%prefix, ionode_id, comm)
    CALL mp_bcast(obj%pseudo_dir, ionode_id, comm)
    CALL mp_bcast(obj%outdir, ionode_id, comm)
    CALL mp_bcast(obj%stress, ionode_id, comm)
    CALL mp_bcast(obj%forces, ionode_id, comm)
    CALL mp_bcast(obj%wf_collect, ionode_id, comm)
    CALL mp_bcast(obj%disk_io, ionode_id, comm)
    CALL mp_bcast(obj%max_seconds, ionode_id, comm)
    CALL mp_bcast(obj%nstep_ispresent, ionode_id, comm)
    IF (obj%nstep_ispresent) &
      CALL mp_bcast(obj%nstep, ionode_id, comm)
    CALL mp_bcast(obj%etot_conv_thr, ionode_id, comm)
    CALL mp_bcast(obj%forc_conv_thr, ionode_id, comm)
    CALL mp_bcast(obj%press_conv_thr, ionode_id, comm)
    CALL mp_bcast(obj%verbosity, ionode_id, comm)
    CALL mp_bcast(obj%print_every, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_control_variables
  !
  !
  SUBROUTINE qes_bcast_xml_format(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(xml_format_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%NAME, ionode_id, comm)
    CALL mp_bcast(obj%VERSION, ionode_id, comm)
    CALL mp_bcast(obj%xml_format, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_xml_format
  !
  !
  SUBROUTINE qes_bcast_creator(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(creator_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%NAME, ionode_id, comm)
    CALL mp_bcast(obj%VERSION, ionode_id, comm)
    CALL mp_bcast(obj%creator, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_creator
  !
  !
  SUBROUTINE qes_bcast_created(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(created_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%DATE, ionode_id, comm)
    CALL mp_bcast(obj%TIME, ionode_id, comm)
    CALL mp_bcast(obj%created, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_created
  !
  !
  SUBROUTINE qes_bcast_atomic_species(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(atomic_species_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%ntyp, ionode_id, comm)
    CALL mp_bcast(obj%pseudo_dir_ispresent, ionode_id, comm)
    IF (obj%pseudo_dir_ispresent) &
      CALL mp_bcast(obj%pseudo_dir, ionode_id, comm)
    CALL mp_bcast(obj%ndim_species, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%species(obj%ndim_species))
    DO i=1, obj%ndim_species
      CALL qes_bcast_species(obj%species(i), ionode_id, comm)
    ENDDO
    !
  END SUBROUTINE qes_bcast_atomic_species
  !
  !
  SUBROUTINE qes_bcast_species(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(species_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%name, ionode_id, comm)
    CALL mp_bcast(obj%mass_ispresent, ionode_id, comm)
    IF (obj%mass_ispresent) &
      CALL mp_bcast(obj%mass, ionode_id, comm)
    CALL mp_bcast(obj%pseudo_file, ionode_id, comm)
    CALL mp_bcast(obj%starting_magnetization_ispresent, ionode_id, comm)
    IF (obj%starting_magnetization_ispresent) &
      CALL mp_bcast(obj%starting_magnetization, ionode_id, comm)
    CALL mp_bcast(obj%spin_teta_ispresent, ionode_id, comm)
    IF (obj%spin_teta_ispresent) &
      CALL mp_bcast(obj%spin_teta, ionode_id, comm)
    CALL mp_bcast(obj%spin_phi_ispresent, ionode_id, comm)
    IF (obj%spin_phi_ispresent) &
      CALL mp_bcast(obj%spin_phi, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_species
  !
  !
  SUBROUTINE qes_bcast_atomic_structure(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(atomic_structure_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%nat, ionode_id, comm)
    CALL mp_bcast(obj%alat_ispresent, ionode_id, comm)
    IF (obj%alat_ispresent) &
      CALL mp_bcast(obj%alat, ionode_id, comm)
    CALL mp_bcast(obj%bravais_index_ispresent, ionode_id, comm)
    IF (obj%bravais_index_ispresent) &
      CALL mp_bcast(obj%bravais_index, ionode_id, comm)
    CALL mp_bcast(obj%atomic_positions_ispresent, ionode_id, comm)
    IF (obj%atomic_positions_ispresent) &
      CALL qes_bcast_atomic_positions(obj%atomic_positions, ionode_id, comm)
    CALL mp_bcast(obj%wyckoff_positions_ispresent, ionode_id, comm)
    IF (obj%wyckoff_positions_ispresent) &
      CALL qes_bcast_wyckoff_positions(obj%wyckoff_positions, ionode_id, comm)
    CALL mp_bcast(obj%crystal_positions_ispresent, ionode_id, comm)
    IF (obj%crystal_positions_ispresent) &
      CALL qes_bcast_atomic_positions(obj%crystal_positions, ionode_id, comm)
    CALL qes_bcast_cell(obj%cell, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_atomic_structure
  !
  !
  SUBROUTINE qes_bcast_atomic_positions(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(atomic_positions_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%ndim_atom, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%atom(obj%ndim_atom))
    DO i=1, obj%ndim_atom
      CALL qes_bcast_atom(obj%atom(i), ionode_id, comm)
    ENDDO
    !
  END SUBROUTINE qes_bcast_atomic_positions
  !
  !
  SUBROUTINE qes_bcast_atom(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(atom_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%name, ionode_id, comm)
    CALL mp_bcast(obj%position_ispresent, ionode_id, comm)
    IF (obj%position_ispresent) &
      CALL mp_bcast(obj%position, ionode_id, comm)
    CALL mp_bcast(obj%index_ispresent, ionode_id, comm)
    IF (obj%index_ispresent) &
      CALL mp_bcast(obj%index, ionode_id, comm)
    CALL mp_bcast(obj%atom, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_atom
  !
  !
  SUBROUTINE qes_bcast_wyckoff_positions(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(wyckoff_positions_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%space_group, ionode_id, comm)
    CALL mp_bcast(obj%more_options_ispresent, ionode_id, comm)
    IF (obj%more_options_ispresent) &
      CALL mp_bcast(obj%more_options, ionode_id, comm)
    CALL mp_bcast(obj%ndim_atom, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%atom(obj%ndim_atom))
    DO i=1, obj%ndim_atom
      CALL qes_bcast_atom(obj%atom(i), ionode_id, comm)
    ENDDO
    !
  END SUBROUTINE qes_bcast_wyckoff_positions
  !
  !
  SUBROUTINE qes_bcast_cell(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(cell_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%a1, ionode_id, comm)
    CALL mp_bcast(obj%a2, ionode_id, comm)
    CALL mp_bcast(obj%a3, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_cell
  !
  !
  SUBROUTINE qes_bcast_dft(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(dft_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%functional, ionode_id, comm)
    CALL mp_bcast(obj%hybrid_ispresent, ionode_id, comm)
    IF (obj%hybrid_ispresent) &
      CALL qes_bcast_hybrid(obj%hybrid, ionode_id, comm)
    CALL mp_bcast(obj%dftU_ispresent, ionode_id, comm)
    IF (obj%dftU_ispresent) &
      CALL qes_bcast_dftU(obj%dftU, ionode_id, comm)
    CALL mp_bcast(obj%vdW_ispresent, ionode_id, comm)
    IF (obj%vdW_ispresent) &
      CALL qes_bcast_vdW(obj%vdW, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_dft
  !
  !
  SUBROUTINE qes_bcast_hybrid(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(hybrid_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%qpoint_grid_ispresent, ionode_id, comm)
    IF (obj%qpoint_grid_ispresent) &
      CALL qes_bcast_qpoint_grid(obj%qpoint_grid, ionode_id, comm)
    CALL mp_bcast(obj%ecutfock_ispresent, ionode_id, comm)
    IF (obj%ecutfock_ispresent) &
      CALL mp_bcast(obj%ecutfock, ionode_id, comm)
    CALL mp_bcast(obj%exx_fraction_ispresent, ionode_id, comm)
    IF (obj%exx_fraction_ispresent) &
      CALL mp_bcast(obj%exx_fraction, ionode_id, comm)
    CALL mp_bcast(obj%screening_parameter_ispresent, ionode_id, comm)
    IF (obj%screening_parameter_ispresent) &
      CALL mp_bcast(obj%screening_parameter, ionode_id, comm)
    CALL mp_bcast(obj%exxdiv_treatment_ispresent, ionode_id, comm)
    IF (obj%exxdiv_treatment_ispresent) &
      CALL mp_bcast(obj%exxdiv_treatment, ionode_id, comm)
    CALL mp_bcast(obj%x_gamma_extrapolation_ispresent, ionode_id, comm)
    IF (obj%x_gamma_extrapolation_ispresent) &
      CALL mp_bcast(obj%x_gamma_extrapolation, ionode_id, comm)
    CALL mp_bcast(obj%ecutvcut_ispresent, ionode_id, comm)
    IF (obj%ecutvcut_ispresent) &
      CALL mp_bcast(obj%ecutvcut, ionode_id, comm)
    CALL mp_bcast(obj%localization_threshold_ispresent, ionode_id, comm)
    IF (obj%localization_threshold_ispresent) &
      CALL mp_bcast(obj%localization_threshold, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_hybrid
  !
  !
  SUBROUTINE qes_bcast_qpoint_grid(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(qpoint_grid_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%nqx1, ionode_id, comm)
    CALL mp_bcast(obj%nqx2, ionode_id, comm)
    CALL mp_bcast(obj%nqx3, ionode_id, comm)
    CALL mp_bcast(obj%qpoint_grid, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_qpoint_grid
  !
  !
  SUBROUTINE qes_bcast_dftU(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(dftU_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%lda_plus_u_kind_ispresent, ionode_id, comm)
    IF (obj%lda_plus_u_kind_ispresent) &
      CALL mp_bcast(obj%lda_plus_u_kind, ionode_id, comm)
    CALL mp_bcast(obj%Hubbard_U_ispresent, ionode_id, comm)
    IF (obj%Hubbard_U_ispresent) THEN
      CALL mp_bcast(obj%ndim_Hubbard_U, ionode_id, comm)
      IF (.NOT.ionode) ALLOCATE(obj%Hubbard_U(obj%ndim_Hubbard_U))
      DO i=1, obj%ndim_Hubbard_U
        CALL qes_bcast_HubbardCommon(obj%Hubbard_U(i), ionode_id, comm)
      ENDDO
    ENDIF
    CALL mp_bcast(obj%Hubbard_J0_ispresent, ionode_id, comm)
    IF (obj%Hubbard_J0_ispresent) THEN
      CALL mp_bcast(obj%ndim_Hubbard_J0, ionode_id, comm)
      IF (.NOT.ionode) ALLOCATE(obj%Hubbard_J0(obj%ndim_Hubbard_J0))
      DO i=1, obj%ndim_Hubbard_J0
        CALL qes_bcast_HubbardCommon(obj%Hubbard_J0(i), ionode_id, comm)
      ENDDO
    ENDIF
    CALL mp_bcast(obj%Hubbard_alpha_ispresent, ionode_id, comm)
    IF (obj%Hubbard_alpha_ispresent) THEN
      CALL mp_bcast(obj%ndim_Hubbard_alpha, ionode_id, comm)
      IF (.NOT.ionode) ALLOCATE(obj%Hubbard_alpha(obj%ndim_Hubbard_alpha))
      DO i=1, obj%ndim_Hubbard_alpha
        CALL qes_bcast_HubbardCommon(obj%Hubbard_alpha(i), ionode_id, comm)
      ENDDO
    ENDIF
    CALL mp_bcast(obj%Hubbard_beta_ispresent, ionode_id, comm)
    IF (obj%Hubbard_beta_ispresent) THEN
      CALL mp_bcast(obj%ndim_Hubbard_beta, ionode_id, comm)
      IF (.NOT.ionode) ALLOCATE(obj%Hubbard_beta(obj%ndim_Hubbard_beta))
      DO i=1, obj%ndim_Hubbard_beta
        CALL qes_bcast_HubbardCommon(obj%Hubbard_beta(i), ionode_id, comm)
      ENDDO
    ENDIF
    CALL mp_bcast(obj%Hubbard_J_ispresent, ionode_id, comm)
    IF (obj%Hubbard_J_ispresent) THEN
      CALL mp_bcast(obj%ndim_Hubbard_J, ionode_id, comm)
      IF (.NOT.ionode) ALLOCATE(obj%Hubbard_J(obj%ndim_Hubbard_J))
      DO i=1, obj%ndim_Hubbard_J
        CALL qes_bcast_HubbardJ(obj%Hubbard_J(i), ionode_id, comm)
      ENDDO
    ENDIF
    CALL mp_bcast(obj%starting_ns_ispresent, ionode_id, comm)
    IF (obj%starting_ns_ispresent) THEN
      CALL mp_bcast(obj%ndim_starting_ns, ionode_id, comm)
      IF (.NOT.ionode) ALLOCATE(obj%starting_ns(obj%ndim_starting_ns))
      DO i=1, obj%ndim_starting_ns
        CALL qes_bcast_starting_ns(obj%starting_ns(i), ionode_id, comm)
      ENDDO
    ENDIF
    CALL mp_bcast(obj%Hubbard_ns_ispresent, ionode_id, comm)
    IF (obj%Hubbard_ns_ispresent) THEN
      CALL mp_bcast(obj%ndim_Hubbard_ns, ionode_id, comm)
      IF (.NOT.ionode) ALLOCATE(obj%Hubbard_ns(obj%ndim_Hubbard_ns))
      DO i=1, obj%ndim_Hubbard_ns
        CALL qes_bcast_Hubbard_ns(obj%Hubbard_ns(i), ionode_id, comm)
      ENDDO
    ENDIF
    CALL mp_bcast(obj%U_projection_type_ispresent, ionode_id, comm)
    IF (obj%U_projection_type_ispresent) &
      CALL mp_bcast(obj%U_projection_type, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_dftU
  !
  !
  SUBROUTINE qes_bcast_HubbardCommon(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(HubbardCommon_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%specie, ionode_id, comm)
    CALL mp_bcast(obj%label_ispresent, ionode_id, comm)
    IF (obj%label_ispresent) &
      CALL mp_bcast(obj%label, ionode_id, comm)
    CALL mp_bcast(obj%HubbardCommon, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_HubbardCommon
  !
  !
  SUBROUTINE qes_bcast_HubbardJ(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(HubbardJ_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%specie, ionode_id, comm)
    CALL mp_bcast(obj%label, ionode_id, comm)
    CALL mp_bcast(obj%HubbardJ, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_HubbardJ
  !
  !
  SUBROUTINE qes_bcast_starting_ns(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(starting_ns_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%specie, ionode_id, comm)
    CALL mp_bcast(obj%label, ionode_id, comm)
    CALL mp_bcast(obj%spin, ionode_id, comm)
    CALL mp_bcast(obj%size, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%starting_ns(obj%size))
    CALL mp_bcast(obj%starting_ns, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_starting_ns
  !
  !
  SUBROUTINE qes_bcast_Hubbard_ns(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(Hubbard_ns_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: length
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%specie, ionode_id, comm)
    CALL mp_bcast(obj%label, ionode_id, comm)
    CALL mp_bcast(obj%spin, ionode_id, comm)
    CALL mp_bcast(obj%index, ionode_id, comm)
    CALL mp_bcast(obj%rank, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%dims(obj%rank))
    CALL mp_bcast(obj%dims, ionode_id, comm)
    CALL mp_bcast(obj%order, ionode_id, comm)
    IF (.NOT. ionode) THEN
      length = 1
      DO i=1, obj%rank
        length = length * obj%dims(i)
      END DO
      ALLOCATE (obj%Hubbard_ns(length) )
    ENDIF
    CALL mp_bcast(obj%Hubbard_ns, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_Hubbard_ns
  !
  !
  SUBROUTINE qes_bcast_vdW(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(vdW_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%vdw_corr_ispresent, ionode_id, comm)
    IF (obj%vdw_corr_ispresent) &
      CALL mp_bcast(obj%vdw_corr, ionode_id, comm)
    CALL mp_bcast(obj%dftd3_version_ispresent, ionode_id, comm)
    IF (obj%dftd3_version_ispresent) &
      CALL mp_bcast(obj%dftd3_version, ionode_id, comm)
    CALL mp_bcast(obj%dftd3_threebody_ispresent, ionode_id, comm)
    IF (obj%dftd3_threebody_ispresent) &
      CALL mp_bcast(obj%dftd3_threebody, ionode_id, comm)
    CALL mp_bcast(obj%non_local_term_ispresent, ionode_id, comm)
    IF (obj%non_local_term_ispresent) &
      CALL mp_bcast(obj%non_local_term, ionode_id, comm)
    CALL mp_bcast(obj%functional_ispresent, ionode_id, comm)
    IF (obj%functional_ispresent) &
      CALL mp_bcast(obj%functional, ionode_id, comm)
    CALL mp_bcast(obj%total_energy_term_ispresent, ionode_id, comm)
    IF (obj%total_energy_term_ispresent) &
      CALL mp_bcast(obj%total_energy_term, ionode_id, comm)
    CALL mp_bcast(obj%london_s6_ispresent, ionode_id, comm)
    IF (obj%london_s6_ispresent) &
      CALL mp_bcast(obj%london_s6, ionode_id, comm)
    CALL mp_bcast(obj%ts_vdw_econv_thr_ispresent, ionode_id, comm)
    IF (obj%ts_vdw_econv_thr_ispresent) &
      CALL mp_bcast(obj%ts_vdw_econv_thr, ionode_id, comm)
    CALL mp_bcast(obj%ts_vdw_isolated_ispresent, ionode_id, comm)
    IF (obj%ts_vdw_isolated_ispresent) &
      CALL mp_bcast(obj%ts_vdw_isolated, ionode_id, comm)
    CALL mp_bcast(obj%london_rcut_ispresent, ionode_id, comm)
    IF (obj%london_rcut_ispresent) &
      CALL mp_bcast(obj%london_rcut, ionode_id, comm)
    CALL mp_bcast(obj%xdm_a1_ispresent, ionode_id, comm)
    IF (obj%xdm_a1_ispresent) &
      CALL mp_bcast(obj%xdm_a1, ionode_id, comm)
    CALL mp_bcast(obj%xdm_a2_ispresent, ionode_id, comm)
    IF (obj%xdm_a2_ispresent) &
      CALL mp_bcast(obj%xdm_a2, ionode_id, comm)
    CALL mp_bcast(obj%london_c6_ispresent, ionode_id, comm)
    IF (obj%london_c6_ispresent) THEN
      CALL mp_bcast(obj%ndim_london_c6, ionode_id, comm)
      IF (.NOT.ionode) ALLOCATE(obj%london_c6(obj%ndim_london_c6))
      DO i=1, obj%ndim_london_c6
        CALL qes_bcast_HubbardCommon(obj%london_c6(i), ionode_id, comm)
      ENDDO
    ENDIF
    !
  END SUBROUTINE qes_bcast_vdW
  !
  !
  SUBROUTINE qes_bcast_spin(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(spin_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%lsda, ionode_id, comm)
    CALL mp_bcast(obj%noncolin, ionode_id, comm)
    CALL mp_bcast(obj%spinorbit, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_spin
  !
  !
  SUBROUTINE qes_bcast_bands(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(bands_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%nbnd_ispresent, ionode_id, comm)
    IF (obj%nbnd_ispresent) &
      CALL mp_bcast(obj%nbnd, ionode_id, comm)
    CALL mp_bcast(obj%smearing_ispresent, ionode_id, comm)
    IF (obj%smearing_ispresent) &
      CALL qes_bcast_smearing(obj%smearing, ionode_id, comm)
    CALL mp_bcast(obj%tot_charge_ispresent, ionode_id, comm)
    IF (obj%tot_charge_ispresent) &
      CALL mp_bcast(obj%tot_charge, ionode_id, comm)
    CALL mp_bcast(obj%tot_magnetization_ispresent, ionode_id, comm)
    IF (obj%tot_magnetization_ispresent) &
      CALL mp_bcast(obj%tot_magnetization, ionode_id, comm)
    CALL qes_bcast_occupations(obj%occupations, ionode_id, comm)
    CALL mp_bcast(obj%inputOccupations_ispresent, ionode_id, comm)
    IF (obj%inputOccupations_ispresent) THEN
      CALL mp_bcast(obj%ndim_inputOccupations, ionode_id, comm)
      IF (.NOT.ionode) ALLOCATE(obj%inputOccupations(obj%ndim_inputOccupations))
      DO i=1, obj%ndim_inputOccupations
        CALL qes_bcast_inputOccupations(obj%inputOccupations(i), ionode_id, comm)
      ENDDO
    ENDIF
    !
  END SUBROUTINE qes_bcast_bands
  !
  !
  SUBROUTINE qes_bcast_smearing(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(smearing_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%degauss, ionode_id, comm)
    CALL mp_bcast(obj%smearing, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_smearing
  !
  !
  SUBROUTINE qes_bcast_occupations(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(occupations_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%spin_ispresent, ionode_id, comm)
    IF (obj%spin_ispresent) &
      CALL mp_bcast(obj%spin, ionode_id, comm)
    CALL mp_bcast(obj%occupations, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_occupations
  !
  !
  SUBROUTINE qes_bcast_basis(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(basis_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%gamma_only_ispresent, ionode_id, comm)
    IF (obj%gamma_only_ispresent) &
      CALL mp_bcast(obj%gamma_only, ionode_id, comm)
    CALL mp_bcast(obj%ecutwfc, ionode_id, comm)
    CALL mp_bcast(obj%ecutrho_ispresent, ionode_id, comm)
    IF (obj%ecutrho_ispresent) &
      CALL mp_bcast(obj%ecutrho, ionode_id, comm)
    CALL mp_bcast(obj%fft_grid_ispresent, ionode_id, comm)
    IF (obj%fft_grid_ispresent) &
      CALL qes_bcast_basisSetItem(obj%fft_grid, ionode_id, comm)
    CALL mp_bcast(obj%fft_smooth_ispresent, ionode_id, comm)
    IF (obj%fft_smooth_ispresent) &
      CALL qes_bcast_basisSetItem(obj%fft_smooth, ionode_id, comm)
    CALL mp_bcast(obj%fft_box_ispresent, ionode_id, comm)
    IF (obj%fft_box_ispresent) &
      CALL qes_bcast_basisSetItem(obj%fft_box, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_basis
  !
  !
  SUBROUTINE qes_bcast_basis_set(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(basis_set_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%gamma_only_ispresent, ionode_id, comm)
    IF (obj%gamma_only_ispresent) &
      CALL mp_bcast(obj%gamma_only, ionode_id, comm)
    CALL mp_bcast(obj%ecutwfc, ionode_id, comm)
    CALL mp_bcast(obj%ecutrho_ispresent, ionode_id, comm)
    IF (obj%ecutrho_ispresent) &
      CALL mp_bcast(obj%ecutrho, ionode_id, comm)
    CALL qes_bcast_basisSetItem(obj%fft_grid, ionode_id, comm)
    CALL mp_bcast(obj%fft_smooth_ispresent, ionode_id, comm)
    IF (obj%fft_smooth_ispresent) &
      CALL qes_bcast_basisSetItem(obj%fft_smooth, ionode_id, comm)
    CALL mp_bcast(obj%fft_box_ispresent, ionode_id, comm)
    IF (obj%fft_box_ispresent) &
      CALL qes_bcast_basisSetItem(obj%fft_box, ionode_id, comm)
    CALL mp_bcast(obj%ngm, ionode_id, comm)
    CALL mp_bcast(obj%ngms_ispresent, ionode_id, comm)
    IF (obj%ngms_ispresent) &
      CALL mp_bcast(obj%ngms, ionode_id, comm)
    CALL mp_bcast(obj%npwx, ionode_id, comm)
    CALL qes_bcast_reciprocal_lattice(obj%reciprocal_lattice, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_basis_set
  !
  !
  SUBROUTINE qes_bcast_basisSetItem(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(basisSetItem_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%nr1, ionode_id, comm)
    CALL mp_bcast(obj%nr2, ionode_id, comm)
    CALL mp_bcast(obj%nr3, ionode_id, comm)
    CALL mp_bcast(obj%basisSetItem, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_basisSetItem
  !
  !
  SUBROUTINE qes_bcast_reciprocal_lattice(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(reciprocal_lattice_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%b1, ionode_id, comm)
    CALL mp_bcast(obj%b2, ionode_id, comm)
    CALL mp_bcast(obj%b3, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_reciprocal_lattice
  !
  !
  SUBROUTINE qes_bcast_electron_control(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(electron_control_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%diagonalization, ionode_id, comm)
    CALL mp_bcast(obj%mixing_mode, ionode_id, comm)
    CALL mp_bcast(obj%mixing_beta, ionode_id, comm)
    CALL mp_bcast(obj%conv_thr, ionode_id, comm)
    CALL mp_bcast(obj%mixing_ndim, ionode_id, comm)
    CALL mp_bcast(obj%max_nstep, ionode_id, comm)
    CALL mp_bcast(obj%real_space_q_ispresent, ionode_id, comm)
    IF (obj%real_space_q_ispresent) &
      CALL mp_bcast(obj%real_space_q, ionode_id, comm)
    CALL mp_bcast(obj%real_space_beta_ispresent, ionode_id, comm)
    IF (obj%real_space_beta_ispresent) &
      CALL mp_bcast(obj%real_space_beta, ionode_id, comm)
    CALL mp_bcast(obj%tq_smoothing, ionode_id, comm)
    CALL mp_bcast(obj%tbeta_smoothing, ionode_id, comm)
    CALL mp_bcast(obj%diago_thr_init, ionode_id, comm)
    CALL mp_bcast(obj%diago_full_acc, ionode_id, comm)
    CALL mp_bcast(obj%diago_cg_maxiter_ispresent, ionode_id, comm)
    IF (obj%diago_cg_maxiter_ispresent) &
      CALL mp_bcast(obj%diago_cg_maxiter, ionode_id, comm)
    CALL mp_bcast(obj%diago_ppcg_maxiter_ispresent, ionode_id, comm)
    IF (obj%diago_ppcg_maxiter_ispresent) &
      CALL mp_bcast(obj%diago_ppcg_maxiter, ionode_id, comm)
    CALL mp_bcast(obj%diago_david_ndim_ispresent, ionode_id, comm)
    IF (obj%diago_david_ndim_ispresent) &
      CALL mp_bcast(obj%diago_david_ndim, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_electron_control
  !
  !
  SUBROUTINE qes_bcast_k_points_IBZ(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(k_points_IBZ_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%monkhorst_pack_ispresent, ionode_id, comm)
    IF (obj%monkhorst_pack_ispresent) &
      CALL qes_bcast_monkhorst_pack(obj%monkhorst_pack, ionode_id, comm)
    CALL mp_bcast(obj%nk_ispresent, ionode_id, comm)
    IF (obj%nk_ispresent) &
      CALL mp_bcast(obj%nk, ionode_id, comm)
    CALL mp_bcast(obj%k_point_ispresent, ionode_id, comm)
    IF (obj%k_point_ispresent) THEN
      CALL mp_bcast(obj%ndim_k_point, ionode_id, comm)
      IF (.NOT.ionode) ALLOCATE(obj%k_point(obj%ndim_k_point))
      DO i=1, obj%ndim_k_point
        CALL qes_bcast_k_point(obj%k_point(i), ionode_id, comm)
      ENDDO
    ENDIF
    !
  END SUBROUTINE qes_bcast_k_points_IBZ
  !
  !
  SUBROUTINE qes_bcast_monkhorst_pack(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(monkhorst_pack_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%nk1, ionode_id, comm)
    CALL mp_bcast(obj%nk2, ionode_id, comm)
    CALL mp_bcast(obj%nk3, ionode_id, comm)
    CALL mp_bcast(obj%k1, ionode_id, comm)
    CALL mp_bcast(obj%k2, ionode_id, comm)
    CALL mp_bcast(obj%k3, ionode_id, comm)
    CALL mp_bcast(obj%monkhorst_pack, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_monkhorst_pack
  !
  !
  SUBROUTINE qes_bcast_k_point(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(k_point_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%weight_ispresent, ionode_id, comm)
    IF (obj%weight_ispresent) &
      CALL mp_bcast(obj%weight, ionode_id, comm)
    CALL mp_bcast(obj%label_ispresent, ionode_id, comm)
    IF (obj%label_ispresent) &
      CALL mp_bcast(obj%label, ionode_id, comm)
    CALL mp_bcast(obj%k_point, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_k_point
  !
  !
  SUBROUTINE qes_bcast_ion_control(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(ion_control_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%ion_dynamics, ionode_id, comm)
    CALL mp_bcast(obj%upscale_ispresent, ionode_id, comm)
    IF (obj%upscale_ispresent) &
      CALL mp_bcast(obj%upscale, ionode_id, comm)
    CALL mp_bcast(obj%remove_rigid_rot_ispresent, ionode_id, comm)
    IF (obj%remove_rigid_rot_ispresent) &
      CALL mp_bcast(obj%remove_rigid_rot, ionode_id, comm)
    CALL mp_bcast(obj%refold_pos_ispresent, ionode_id, comm)
    IF (obj%refold_pos_ispresent) &
      CALL mp_bcast(obj%refold_pos, ionode_id, comm)
    CALL mp_bcast(obj%bfgs_ispresent, ionode_id, comm)
    IF (obj%bfgs_ispresent) &
      CALL qes_bcast_bfgs(obj%bfgs, ionode_id, comm)
    CALL mp_bcast(obj%md_ispresent, ionode_id, comm)
    IF (obj%md_ispresent) &
      CALL qes_bcast_md(obj%md, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_ion_control
  !
  !
  SUBROUTINE qes_bcast_bfgs(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(bfgs_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%ndim, ionode_id, comm)
    CALL mp_bcast(obj%trust_radius_min, ionode_id, comm)
    CALL mp_bcast(obj%trust_radius_max, ionode_id, comm)
    CALL mp_bcast(obj%trust_radius_init, ionode_id, comm)
    CALL mp_bcast(obj%w1, ionode_id, comm)
    CALL mp_bcast(obj%w2, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_bfgs
  !
  !
  SUBROUTINE qes_bcast_md(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(md_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%pot_extrapolation, ionode_id, comm)
    CALL mp_bcast(obj%wfc_extrapolation, ionode_id, comm)
    CALL mp_bcast(obj%ion_temperature, ionode_id, comm)
    CALL mp_bcast(obj%timestep, ionode_id, comm)
    CALL mp_bcast(obj%tempw, ionode_id, comm)
    CALL mp_bcast(obj%tolp, ionode_id, comm)
    CALL mp_bcast(obj%deltaT, ionode_id, comm)
    CALL mp_bcast(obj%nraise, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_md
  !
  !
  SUBROUTINE qes_bcast_cell_control(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(cell_control_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%cell_dynamics, ionode_id, comm)
    CALL mp_bcast(obj%pressure, ionode_id, comm)
    CALL mp_bcast(obj%wmass_ispresent, ionode_id, comm)
    IF (obj%wmass_ispresent) &
      CALL mp_bcast(obj%wmass, ionode_id, comm)
    CALL mp_bcast(obj%cell_factor_ispresent, ionode_id, comm)
    IF (obj%cell_factor_ispresent) &
      CALL mp_bcast(obj%cell_factor, ionode_id, comm)
    CALL mp_bcast(obj%fix_volume_ispresent, ionode_id, comm)
    IF (obj%fix_volume_ispresent) &
      CALL mp_bcast(obj%fix_volume, ionode_id, comm)
    CALL mp_bcast(obj%fix_area_ispresent, ionode_id, comm)
    IF (obj%fix_area_ispresent) &
      CALL mp_bcast(obj%fix_area, ionode_id, comm)
    CALL mp_bcast(obj%isotropic_ispresent, ionode_id, comm)
    IF (obj%isotropic_ispresent) &
      CALL mp_bcast(obj%isotropic, ionode_id, comm)
    CALL mp_bcast(obj%free_cell_ispresent, ionode_id, comm)
    IF (obj%free_cell_ispresent) &
      CALL qes_bcast_integerMatrix(obj%free_cell, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_cell_control
  !
  !
  SUBROUTINE qes_bcast_symmetry_flags(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(symmetry_flags_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%nosym, ionode_id, comm)
    CALL mp_bcast(obj%nosym_evc, ionode_id, comm)
    CALL mp_bcast(obj%noinv, ionode_id, comm)
    CALL mp_bcast(obj%no_t_rev, ionode_id, comm)
    CALL mp_bcast(obj%force_symmorphic, ionode_id, comm)
    CALL mp_bcast(obj%use_all_frac, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_symmetry_flags
  !
  !
  SUBROUTINE qes_bcast_boundary_conditions(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(boundary_conditions_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%assume_isolated, ionode_id, comm)
    CALL mp_bcast(obj%esm_ispresent, ionode_id, comm)
    IF (obj%esm_ispresent) &
      CALL qes_bcast_esm(obj%esm, ionode_id, comm)
    CALL mp_bcast(obj%fcp_opt_ispresent, ionode_id, comm)
    IF (obj%fcp_opt_ispresent) &
      CALL mp_bcast(obj%fcp_opt, ionode_id, comm)
    CALL mp_bcast(obj%fcp_mu_ispresent, ionode_id, comm)
    IF (obj%fcp_mu_ispresent) &
      CALL mp_bcast(obj%fcp_mu, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_boundary_conditions
  !
  !
  SUBROUTINE qes_bcast_esm(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(esm_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%bc, ionode_id, comm)
    CALL mp_bcast(obj%nfit, ionode_id, comm)
    CALL mp_bcast(obj%w, ionode_id, comm)
    CALL mp_bcast(obj%efield, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_esm
  !
  !
  SUBROUTINE qes_bcast_ekin_functional(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(ekin_functional_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%ecfixed, ionode_id, comm)
    CALL mp_bcast(obj%qcutz, ionode_id, comm)
    CALL mp_bcast(obj%q2sigma, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_ekin_functional
  !
  !
  SUBROUTINE qes_bcast_spin_constraints(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(spin_constraints_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%spin_constraints, ionode_id, comm)
    CALL mp_bcast(obj%lagrange_multiplier, ionode_id, comm)
    CALL mp_bcast(obj%target_magnetization_ispresent, ionode_id, comm)
    IF (obj%target_magnetization_ispresent) &
      CALL mp_bcast(obj%target_magnetization, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_spin_constraints
  !
  !
  SUBROUTINE qes_bcast_electric_field(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(electric_field_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%electric_potential, ionode_id, comm)
    CALL mp_bcast(obj%dipole_correction_ispresent, ionode_id, comm)
    IF (obj%dipole_correction_ispresent) &
      CALL mp_bcast(obj%dipole_correction, ionode_id, comm)
    CALL mp_bcast(obj%gate_settings_ispresent, ionode_id, comm)
    IF (obj%gate_settings_ispresent) &
      CALL qes_bcast_gate_settings(obj%gate_settings, ionode_id, comm)
    CALL mp_bcast(obj%electric_field_direction_ispresent, ionode_id, comm)
    IF (obj%electric_field_direction_ispresent) &
      CALL mp_bcast(obj%electric_field_direction, ionode_id, comm)
    CALL mp_bcast(obj%potential_max_position_ispresent, ionode_id, comm)
    IF (obj%potential_max_position_ispresent) &
      CALL mp_bcast(obj%potential_max_position, ionode_id, comm)
    CALL mp_bcast(obj%potential_decrease_width_ispresent, ionode_id, comm)
    IF (obj%potential_decrease_width_ispresent) &
      CALL mp_bcast(obj%potential_decrease_width, ionode_id, comm)
    CALL mp_bcast(obj%electric_field_amplitude_ispresent, ionode_id, comm)
    IF (obj%electric_field_amplitude_ispresent) &
      CALL mp_bcast(obj%electric_field_amplitude, ionode_id, comm)
    CALL mp_bcast(obj%electric_field_vector_ispresent, ionode_id, comm)
    IF (obj%electric_field_vector_ispresent) &
      CALL mp_bcast(obj%electric_field_vector, ionode_id, comm)
    CALL mp_bcast(obj%nk_per_string_ispresent, ionode_id, comm)
    IF (obj%nk_per_string_ispresent) &
      CALL mp_bcast(obj%nk_per_string, ionode_id, comm)
    CALL mp_bcast(obj%n_berry_cycles_ispresent, ionode_id, comm)
    IF (obj%n_berry_cycles_ispresent) &
      CALL mp_bcast(obj%n_berry_cycles, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_electric_field
  !
  !
  SUBROUTINE qes_bcast_gate_settings(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(gate_settings_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%use_gate, ionode_id, comm)
    CALL mp_bcast(obj%zgate_ispresent, ionode_id, comm)
    IF (obj%zgate_ispresent) &
      CALL mp_bcast(obj%zgate, ionode_id, comm)
    CALL mp_bcast(obj%relaxz_ispresent, ionode_id, comm)
    IF (obj%relaxz_ispresent) &
      CALL mp_bcast(obj%relaxz, ionode_id, comm)
    CALL mp_bcast(obj%block_ispresent, ionode_id, comm)
    IF (obj%block_ispresent) &
      CALL mp_bcast(obj%block, ionode_id, comm)
    CALL mp_bcast(obj%block_1_ispresent, ionode_id, comm)
    IF (obj%block_1_ispresent) &
      CALL mp_bcast(obj%block_1, ionode_id, comm)
    CALL mp_bcast(obj%block_2_ispresent, ionode_id, comm)
    IF (obj%block_2_ispresent) &
      CALL mp_bcast(obj%block_2, ionode_id, comm)
    CALL mp_bcast(obj%block_height_ispresent, ionode_id, comm)
    IF (obj%block_height_ispresent) &
      CALL mp_bcast(obj%block_height, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_gate_settings
  !
  !
  SUBROUTINE qes_bcast_atomic_constraints(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(atomic_constraints_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%num_of_constraints, ionode_id, comm)
    CALL mp_bcast(obj%tolerance, ionode_id, comm)
    CALL mp_bcast(obj%ndim_atomic_constraint, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%atomic_constraint(obj%ndim_atomic_constraint))
    DO i=1, obj%ndim_atomic_constraint
      CALL qes_bcast_atomic_constraint(obj%atomic_constraint(i), ionode_id, comm)
    ENDDO
    !
  END SUBROUTINE qes_bcast_atomic_constraints
  !
  !
  SUBROUTINE qes_bcast_atomic_constraint(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(atomic_constraint_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%constr_parms, ionode_id, comm)
    CALL mp_bcast(obj%constr_type, ionode_id, comm)
    CALL mp_bcast(obj%constr_target, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_atomic_constraint
  !
  !
  SUBROUTINE qes_bcast_inputOccupations(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(inputOccupations_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%ispin, ionode_id, comm)
    CALL mp_bcast(obj%spin_factor, ionode_id, comm)
    CALL mp_bcast(obj%size, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%inputOccupations(obj%size))
    CALL mp_bcast(obj%inputOccupations, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_inputOccupations
  !
  !
  SUBROUTINE qes_bcast_outputElectricField(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(outputElectricField_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%BerryPhase_ispresent, ionode_id, comm)
    IF (obj%BerryPhase_ispresent) &
      CALL qes_bcast_BerryPhaseOutput(obj%BerryPhase, ionode_id, comm)
    CALL mp_bcast(obj%finiteElectricFieldInfo_ispresent, ionode_id, comm)
    IF (obj%finiteElectricFieldInfo_ispresent) &
      CALL qes_bcast_finiteFieldOut(obj%finiteElectricFieldInfo, ionode_id, comm)
    CALL mp_bcast(obj%dipoleInfo_ispresent, ionode_id, comm)
    IF (obj%dipoleInfo_ispresent) &
      CALL qes_bcast_dipoleOutput(obj%dipoleInfo, ionode_id, comm)
    CALL mp_bcast(obj%gateInfo_ispresent, ionode_id, comm)
    IF (obj%gateInfo_ispresent) &
      CALL qes_bcast_gateInfo(obj%gateInfo, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_outputElectricField
  !
  !
  SUBROUTINE qes_bcast_BerryPhaseOutput(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(BerryPhaseOutput_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL qes_bcast_polarization(obj%totalPolarization, ionode_id, comm)
    CALL qes_bcast_phase(obj%totalPhase, ionode_id, comm)
    CALL mp_bcast(obj%ndim_ionicPolarization, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%ionicPolarization(obj%ndim_ionicPolarization))
    DO i=1, obj%ndim_ionicPolarization
      CALL qes_bcast_ionicPolarization(obj%ionicPolarization(i), ionode_id, comm)
    ENDDO
    CALL mp_bcast(obj%ndim_electronicPolarization, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%electronicPolarization(obj%ndim_electronicPolarization))
    DO i=1, obj%ndim_electronicPolarization
      CALL qes_bcast_electronicPolarization(obj%electronicPolarization(i), ionode_id, comm)
    ENDDO
    !
  END SUBROUTINE qes_bcast_BerryPhaseOutput
  !
  !
  SUBROUTINE qes_bcast_dipoleOutput(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(dipoleOutput_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%idir, ionode_id, comm)
    CALL qes_bcast_scalarQuantity(obj%dipole, ionode_id, comm)
    CALL qes_bcast_scalarQuantity(obj%ion_dipole, ionode_id, comm)
    CALL qes_bcast_scalarQuantity(obj%elec_dipole, ionode_id, comm)
    CALL qes_bcast_scalarQuantity(obj%dipoleField, ionode_id, comm)
    CALL qes_bcast_scalarQuantity(obj%potentialAmp, ionode_id, comm)
    CALL qes_bcast_scalarQuantity(obj%totalLength, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_dipoleOutput
  !
  !
  SUBROUTINE qes_bcast_finiteFieldOut(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(finiteFieldOut_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%electronicDipole, ionode_id, comm)
    CALL mp_bcast(obj%ionicDipole, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_finiteFieldOut
  !
  !
  SUBROUTINE qes_bcast_polarization(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(polarization_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL qes_bcast_scalarQuantity(obj%polarization, ionode_id, comm)
    CALL mp_bcast(obj%modulus, ionode_id, comm)
    CALL mp_bcast(obj%direction, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_polarization
  !
  !
  SUBROUTINE qes_bcast_ionicPolarization(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(ionicPolarization_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL qes_bcast_atom(obj%ion, ionode_id, comm)
    CALL mp_bcast(obj%charge, ionode_id, comm)
    CALL qes_bcast_phase(obj%phase, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_ionicPolarization
  !
  !
  SUBROUTINE qes_bcast_electronicPolarization(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(electronicPolarization_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL qes_bcast_k_point(obj%firstKeyPoint, ionode_id, comm)
    CALL mp_bcast(obj%spin_ispresent, ionode_id, comm)
    IF (obj%spin_ispresent) &
      CALL mp_bcast(obj%spin, ionode_id, comm)
    CALL qes_bcast_phase(obj%phase, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_electronicPolarization
  !
  !
  SUBROUTINE qes_bcast_phase(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(phase_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%ionic_ispresent, ionode_id, comm)
    IF (obj%ionic_ispresent) &
      CALL mp_bcast(obj%ionic, ionode_id, comm)
    CALL mp_bcast(obj%electronic_ispresent, ionode_id, comm)
    IF (obj%electronic_ispresent) &
      CALL mp_bcast(obj%electronic, ionode_id, comm)
    CALL mp_bcast(obj%modulus_ispresent, ionode_id, comm)
    IF (obj%modulus_ispresent) &
      CALL mp_bcast(obj%modulus, ionode_id, comm)
    CALL mp_bcast(obj%phase, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_phase
  !
  !
  SUBROUTINE qes_bcast_gateInfo(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(gateInfo_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%pot_prefactor, ionode_id, comm)
    CALL mp_bcast(obj%gate_zpos, ionode_id, comm)
    CALL mp_bcast(obj%gate_gate_term, ionode_id, comm)
    CALL mp_bcast(obj%gatefieldEnergy, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_gateInfo
  !
  !
  SUBROUTINE qes_bcast_convergence_info(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(convergence_info_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL qes_bcast_scf_conv(obj%scf_conv, ionode_id, comm)
    CALL mp_bcast(obj%opt_conv_ispresent, ionode_id, comm)
    IF (obj%opt_conv_ispresent) &
      CALL qes_bcast_opt_conv(obj%opt_conv, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_convergence_info
  !
  !
  SUBROUTINE qes_bcast_scf_conv(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(scf_conv_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%convergence_achieved, ionode_id, comm)
    CALL mp_bcast(obj%n_scf_steps, ionode_id, comm)
    CALL mp_bcast(obj%scf_error, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_scf_conv
  !
  !
  SUBROUTINE qes_bcast_opt_conv(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(opt_conv_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%convergence_achieved, ionode_id, comm)
    CALL mp_bcast(obj%n_opt_steps, ionode_id, comm)
    CALL mp_bcast(obj%grad_norm, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_opt_conv
  !
  !
  SUBROUTINE qes_bcast_algorithmic_info(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(algorithmic_info_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%real_space_q, ionode_id, comm)
    CALL mp_bcast(obj%real_space_beta_ispresent, ionode_id, comm)
    IF (obj%real_space_beta_ispresent) &
      CALL mp_bcast(obj%real_space_beta, ionode_id, comm)
    CALL mp_bcast(obj%uspp, ionode_id, comm)
    CALL mp_bcast(obj%paw, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_algorithmic_info
  !
  !
  SUBROUTINE qes_bcast_symmetries(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(symmetries_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%nsym, ionode_id, comm)
    CALL mp_bcast(obj%nrot, ionode_id, comm)
    CALL mp_bcast(obj%space_group, ionode_id, comm)
    CALL mp_bcast(obj%ndim_symmetry, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%symmetry(obj%ndim_symmetry))
    DO i=1, obj%ndim_symmetry
      CALL qes_bcast_symmetry(obj%symmetry(i), ionode_id, comm)
    ENDDO
    !
  END SUBROUTINE qes_bcast_symmetries
  !
  !
  SUBROUTINE qes_bcast_symmetry(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(symmetry_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL qes_bcast_info(obj%info, ionode_id, comm)
    CALL qes_bcast_matrix(obj%rotation, ionode_id, comm)
    CALL mp_bcast(obj%fractional_translation_ispresent, ionode_id, comm)
    IF (obj%fractional_translation_ispresent) &
      CALL mp_bcast(obj%fractional_translation, ionode_id, comm)
    CALL mp_bcast(obj%equivalent_atoms_ispresent, ionode_id, comm)
    IF (obj%equivalent_atoms_ispresent) &
      CALL qes_bcast_equivalent_atoms(obj%equivalent_atoms, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_symmetry
  !
  !
  SUBROUTINE qes_bcast_equivalent_atoms(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(equivalent_atoms_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%nat, ionode_id, comm)
    CALL mp_bcast(obj%size, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%equivalent_atoms(obj%size))
    CALL mp_bcast(obj%equivalent_atoms, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_equivalent_atoms
  !
  !
  SUBROUTINE qes_bcast_info(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(info_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%name_ispresent, ionode_id, comm)
    IF (obj%name_ispresent) &
      CALL mp_bcast(obj%name, ionode_id, comm)
    CALL mp_bcast(obj%class_ispresent, ionode_id, comm)
    IF (obj%class_ispresent) &
      CALL mp_bcast(obj%class, ionode_id, comm)
    CALL mp_bcast(obj%time_reversal_ispresent, ionode_id, comm)
    IF (obj%time_reversal_ispresent) &
      CALL mp_bcast(obj%time_reversal, ionode_id, comm)
    CALL mp_bcast(obj%info, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_info
  !
  !
  SUBROUTINE qes_bcast_outputPBC(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(outputPBC_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%assume_isolated, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_outputPBC
  !
  !
  SUBROUTINE qes_bcast_magnetization(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(magnetization_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%lsda, ionode_id, comm)
    CALL mp_bcast(obj%noncolin, ionode_id, comm)
    CALL mp_bcast(obj%spinorbit, ionode_id, comm)
    CALL mp_bcast(obj%total, ionode_id, comm)
    CALL mp_bcast(obj%absolute, ionode_id, comm)
    CALL mp_bcast(obj%do_magnetization, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_magnetization
  !
  !
  SUBROUTINE qes_bcast_total_energy(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(total_energy_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%etot, ionode_id, comm)
    CALL mp_bcast(obj%eband_ispresent, ionode_id, comm)
    IF (obj%eband_ispresent) &
      CALL mp_bcast(obj%eband, ionode_id, comm)
    CALL mp_bcast(obj%ehart_ispresent, ionode_id, comm)
    IF (obj%ehart_ispresent) &
      CALL mp_bcast(obj%ehart, ionode_id, comm)
    CALL mp_bcast(obj%vtxc_ispresent, ionode_id, comm)
    IF (obj%vtxc_ispresent) &
      CALL mp_bcast(obj%vtxc, ionode_id, comm)
    CALL mp_bcast(obj%etxc_ispresent, ionode_id, comm)
    IF (obj%etxc_ispresent) &
      CALL mp_bcast(obj%etxc, ionode_id, comm)
    CALL mp_bcast(obj%ewald_ispresent, ionode_id, comm)
    IF (obj%ewald_ispresent) &
      CALL mp_bcast(obj%ewald, ionode_id, comm)
    CALL mp_bcast(obj%demet_ispresent, ionode_id, comm)
    IF (obj%demet_ispresent) &
      CALL mp_bcast(obj%demet, ionode_id, comm)
    CALL mp_bcast(obj%efieldcorr_ispresent, ionode_id, comm)
    IF (obj%efieldcorr_ispresent) &
      CALL mp_bcast(obj%efieldcorr, ionode_id, comm)
    CALL mp_bcast(obj%potentiostat_contr_ispresent, ionode_id, comm)
    IF (obj%potentiostat_contr_ispresent) &
      CALL mp_bcast(obj%potentiostat_contr, ionode_id, comm)
    CALL mp_bcast(obj%gatefield_contr_ispresent, ionode_id, comm)
    IF (obj%gatefield_contr_ispresent) &
      CALL mp_bcast(obj%gatefield_contr, ionode_id, comm)
    CALL mp_bcast(obj%vdW_term_ispresent, ionode_id, comm)
    IF (obj%vdW_term_ispresent) &
      CALL mp_bcast(obj%vdW_term, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_total_energy
  !
  !
  SUBROUTINE qes_bcast_band_structure(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(band_structure_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%lsda, ionode_id, comm)
    CALL mp_bcast(obj%noncolin, ionode_id, comm)
    CALL mp_bcast(obj%spinorbit, ionode_id, comm)
    CALL mp_bcast(obj%nbnd_ispresent, ionode_id, comm)
    IF (obj%nbnd_ispresent) &
      CALL mp_bcast(obj%nbnd, ionode_id, comm)
    CALL mp_bcast(obj%nbnd_up_ispresent, ionode_id, comm)
    IF (obj%nbnd_up_ispresent) &
      CALL mp_bcast(obj%nbnd_up, ionode_id, comm)
    CALL mp_bcast(obj%nbnd_dw_ispresent, ionode_id, comm)
    IF (obj%nbnd_dw_ispresent) &
      CALL mp_bcast(obj%nbnd_dw, ionode_id, comm)
    CALL mp_bcast(obj%nelec, ionode_id, comm)
    CALL mp_bcast(obj%num_of_atomic_wfc_ispresent, ionode_id, comm)
    IF (obj%num_of_atomic_wfc_ispresent) &
      CALL mp_bcast(obj%num_of_atomic_wfc, ionode_id, comm)
    CALL mp_bcast(obj%wf_collected, ionode_id, comm)
    CALL mp_bcast(obj%fermi_energy_ispresent, ionode_id, comm)
    IF (obj%fermi_energy_ispresent) &
      CALL mp_bcast(obj%fermi_energy, ionode_id, comm)
    CALL mp_bcast(obj%highestOccupiedLevel_ispresent, ionode_id, comm)
    IF (obj%highestOccupiedLevel_ispresent) &
      CALL mp_bcast(obj%highestOccupiedLevel, ionode_id, comm)
    CALL mp_bcast(obj%lowestUnoccupiedLevel_ispresent, ionode_id, comm)
    IF (obj%lowestUnoccupiedLevel_ispresent) &
      CALL mp_bcast(obj%lowestUnoccupiedLevel, ionode_id, comm)
    CALL mp_bcast(obj%two_fermi_energies_ispresent, ionode_id, comm)
    IF (obj%two_fermi_energies_ispresent) &
      CALL mp_bcast(obj%two_fermi_energies, ionode_id, comm)
    CALL qes_bcast_k_points_IBZ(obj%starting_k_points, ionode_id, comm)
    CALL mp_bcast(obj%nks, ionode_id, comm)
    CALL qes_bcast_occupations(obj%occupations_kind, ionode_id, comm)
    CALL mp_bcast(obj%smearing_ispresent, ionode_id, comm)
    IF (obj%smearing_ispresent) &
      CALL qes_bcast_smearing(obj%smearing, ionode_id, comm)
    CALL mp_bcast(obj%ndim_ks_energies, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%ks_energies(obj%ndim_ks_energies))
    DO i=1, obj%ndim_ks_energies
      CALL qes_bcast_ks_energies(obj%ks_energies(i), ionode_id, comm)
    ENDDO
    !
  END SUBROUTINE qes_bcast_band_structure
  !
  !
  SUBROUTINE qes_bcast_ks_energies(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(ks_energies_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL qes_bcast_k_point(obj%k_point, ionode_id, comm)
    CALL mp_bcast(obj%npw, ionode_id, comm)
    CALL qes_bcast_vector(obj%eigenvalues, ionode_id, comm)
    CALL qes_bcast_vector(obj%occupations, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_ks_energies
  !
  !
  SUBROUTINE qes_bcast_closed(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(closed_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%DATE, ionode_id, comm)
    CALL mp_bcast(obj%TIME, ionode_id, comm)
    CALL mp_bcast(obj%closed, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_closed
  !
  !
  SUBROUTINE qes_bcast_vector(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(vector_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%size, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%vector(obj%size))
    CALL mp_bcast(obj%vector, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_vector
  !
  !
  SUBROUTINE qes_bcast_integerVector(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(integerVector_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%size, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%integerVector(obj%size))
    CALL mp_bcast(obj%integerVector, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_integerVector
  !
  !
  SUBROUTINE qes_bcast_matrix(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(matrix_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: length
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%rank, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%dims(obj%rank))
    CALL mp_bcast(obj%dims, ionode_id, comm)
    CALL mp_bcast(obj%order, ionode_id, comm)
    IF (.NOT. ionode) THEN
      length = 1
      DO i=1, obj%rank
        length = length * obj%dims(i)
      END DO
      ALLOCATE (obj%matrix(length) )
    ENDIF
    CALL mp_bcast(obj%matrix, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_matrix
  !
  !
  SUBROUTINE qes_bcast_integerMatrix(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(integerMatrix_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    INTEGER :: length
    INTEGER :: i
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%rank, ionode_id, comm)
    IF (.NOT.ionode) ALLOCATE(obj%dims(obj%rank))
    CALL mp_bcast(obj%dims, ionode_id, comm)
    CALL mp_bcast(obj%order, ionode_id, comm)
    IF (.NOT. ionode) THEN
      length = 1
      DO i=1, obj%rank
        length = length * obj%dims(i)
      END DO
      ALLOCATE (obj%integerMatrix(length) )
    ENDIF
    CALL mp_bcast(obj%integerMatrix, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_integerMatrix
  !
  !
  SUBROUTINE qes_bcast_scalarQuantity(obj, ionode_id, comm )
    !
    IMPLICIT NONE
    !
    TYPE(scalarQuantity_type), INTENT(INOUT) :: obj
    INTEGER, INTENT(IN) :: ionode_id, comm
    !
    CALL mp_bcast(obj%tagname, ionode_id, comm)
    CALL mp_bcast(obj%lwrite, ionode_id, comm)
    CALL mp_bcast(obj%lread, ionode_id, comm)
    !
    CALL mp_bcast(obj%Units, ionode_id, comm)
    CALL mp_bcast(obj%scalarQuantity, ionode_id, comm)
    !
  END SUBROUTINE qes_bcast_scalarQuantity
  !
  !
END MODULE qes_bcast_module