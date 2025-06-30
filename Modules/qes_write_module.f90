!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE qes_write_module
  !
  ! Auto-generated code: don't edit this file
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes-1.0
  !
#if defined (__fox) 
  USE  FoX_wxml 
#else 
  USE wxml 
#endif 
  USE qes_types_module
  !
  IMPLICIT NONE
  !
  PUBLIC
  !
  INTERFACE qes_write
    MODULE PROCEDURE qes_write_espresso
    MODULE PROCEDURE qes_write_general_info
    MODULE PROCEDURE qes_write_parallel_info
    MODULE PROCEDURE qes_write_input
    MODULE PROCEDURE qes_write_step
    MODULE PROCEDURE qes_write_output
    MODULE PROCEDURE qes_write_timing
    MODULE PROCEDURE qes_write_clock
    MODULE PROCEDURE qes_write_control_variables
    MODULE PROCEDURE qes_write_xml_format
    MODULE PROCEDURE qes_write_creator
    MODULE PROCEDURE qes_write_created
    MODULE PROCEDURE qes_write_atomic_species
    MODULE PROCEDURE qes_write_species
    MODULE PROCEDURE qes_write_atomic_structure
    MODULE PROCEDURE qes_write_atomic_positions
    MODULE PROCEDURE qes_write_atom
    MODULE PROCEDURE qes_write_wyckoff_positions
    MODULE PROCEDURE qes_write_cell
    MODULE PROCEDURE qes_write_dft
    MODULE PROCEDURE qes_write_hybrid
    MODULE PROCEDURE qes_write_qpoint_grid
    MODULE PROCEDURE qes_write_dftU
    MODULE PROCEDURE qes_write_HubbardCommon
    MODULE PROCEDURE qes_write_HubbardInterSpecieV
    MODULE PROCEDURE qes_write_SiteMoment
    MODULE PROCEDURE qes_write_HubbardJ
    MODULE PROCEDURE qes_write_vector
    MODULE PROCEDURE qes_write_HubbardM
    MODULE PROCEDURE qes_write_ChannelOcc
    MODULE PROCEDURE qes_write_HubbardOcc
    MODULE PROCEDURE qes_write_SitMag
    MODULE PROCEDURE qes_write_starting_ns
    MODULE PROCEDURE qes_write_integerVector
    MODULE PROCEDURE qes_write_orderUm
    MODULE PROCEDURE qes_write_matrix
    MODULE PROCEDURE qes_write_Hubbard_ns
    MODULE PROCEDURE qes_write_HubbardBack
    MODULE PROCEDURE qes_write_vdW
    MODULE PROCEDURE qes_write_spin
    MODULE PROCEDURE qes_write_bands
    MODULE PROCEDURE qes_write_smearing
    MODULE PROCEDURE qes_write_occupations
    MODULE PROCEDURE qes_write_basis
    MODULE PROCEDURE qes_write_basis_set
    MODULE PROCEDURE qes_write_basisSetItem
    MODULE PROCEDURE qes_write_reciprocal_lattice
    MODULE PROCEDURE qes_write_electron_control
    MODULE PROCEDURE qes_write_fcp
    MODULE PROCEDURE qes_write_rism
    MODULE PROCEDURE qes_write_solute
    MODULE PROCEDURE qes_write_solvent
    MODULE PROCEDURE qes_write_k_points_IBZ
    MODULE PROCEDURE qes_write_monkhorst_pack
    MODULE PROCEDURE qes_write_k_point
    MODULE PROCEDURE qes_write_ion_control
    MODULE PROCEDURE qes_write_bfgs
    MODULE PROCEDURE qes_write_md
    MODULE PROCEDURE qes_write_cell_control
    MODULE PROCEDURE qes_write_symmetry_flags
    MODULE PROCEDURE qes_write_boundary_conditions
    MODULE PROCEDURE qes_write_esm
    MODULE PROCEDURE qes_write_gcscf
    MODULE PROCEDURE qes_write_solvents
    MODULE PROCEDURE qes_write_ekin_functional
    MODULE PROCEDURE qes_write_spin_constraints
    MODULE PROCEDURE qes_write_electric_field
    MODULE PROCEDURE qes_write_gate_settings
    MODULE PROCEDURE qes_write_atomic_constraints
    MODULE PROCEDURE qes_write_atomic_constraint
    MODULE PROCEDURE qes_write_inputOccupations
    MODULE PROCEDURE qes_write_outputElectricField
    MODULE PROCEDURE qes_write_BerryPhaseOutput
    MODULE PROCEDURE qes_write_sawtoothEnergy
    MODULE PROCEDURE qes_write_dipoleOutput
    MODULE PROCEDURE qes_write_finiteFieldOut
    MODULE PROCEDURE qes_write_polarization
    MODULE PROCEDURE qes_write_ionicPolarization
    MODULE PROCEDURE qes_write_electronicPolarization
    MODULE PROCEDURE qes_write_phase
    MODULE PROCEDURE qes_write_gateInfo
    MODULE PROCEDURE qes_write_convergence_info
    MODULE PROCEDURE qes_write_scf_conv
    MODULE PROCEDURE qes_write_opt_conv
    MODULE PROCEDURE qes_write_algorithmic_info
    MODULE PROCEDURE qes_write_symmetries
    MODULE PROCEDURE qes_write_symmetry
    MODULE PROCEDURE qes_write_equivalent_atoms
    MODULE PROCEDURE qes_write_info
    MODULE PROCEDURE qes_write_outputPBC
    MODULE PROCEDURE qes_write_magnetization
    MODULE PROCEDURE qes_write_total_energy
    MODULE PROCEDURE qes_write_band_structure
    MODULE PROCEDURE qes_write_ks_energies
    MODULE PROCEDURE qes_write_closed
    MODULE PROCEDURE qes_write_cpstatus
    MODULE PROCEDURE qes_write_cpnumstep
    MODULE PROCEDURE qes_write_cptimesteps
    MODULE PROCEDURE qes_write_cpstep
    MODULE PROCEDURE qes_write_cp_ionPos
    MODULE PROCEDURE qes_write_cp_ionsNose
    MODULE PROCEDURE qes_write_cp_elecNose
    MODULE PROCEDURE qes_write_cp_cell
    MODULE PROCEDURE qes_write_cp_cellNose
    MODULE PROCEDURE qes_write_scalmags
    MODULE PROCEDURE qes_write_d3mags
    MODULE PROCEDURE qes_write_integerMatrix
    MODULE PROCEDURE qes_write_scalarQuantity
    MODULE PROCEDURE qes_write_rism3d
    MODULE PROCEDURE qes_write_rismlaue
    MODULE PROCEDURE qes_write_two_chem
  END INTERFACE qes_write
  !
  CONTAINS
  !
  
   SUBROUTINE qes_write_espresso(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(espresso_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%Units_ispresent) CALL xml_addAttribute(xp, 'Units', TRIM(obj%Units) )
     IF (obj%general_info_ispresent) THEN
        CALL qes_write_general_info (xp, obj%general_info)
     END IF
     IF (obj%parallel_info_ispresent) THEN
        CALL qes_write_parallel_info (xp, obj%parallel_info)
     END IF
     IF (obj%input_ispresent) THEN
        CALL qes_write_input (xp, obj%input)
     END IF
     IF (obj%step_ispresent) THEN
        DO i = 1, obj%ndim_step
           CALL qes_write_step(xp, obj%step(i) )
        END DO
     END IF
     IF (obj%output_ispresent) THEN
        CALL qes_write_output (xp, obj%output)
     END IF
     IF (obj%STATUS_ispresent) THEN
        CALL qes_write_cpstatus (xp, obj%STATUS)
     END IF
     IF (obj%TIMESTEPS_ispresent) THEN
        CALL qes_write_cptimesteps (xp, obj%TIMESTEPS)
     END IF
     IF (obj%exit_status_ispresent) THEN
        CALL xml_NewElement(xp, "exit_status")
           CALL xml_addCharacters(xp, obj%exit_status)
        CALL xml_EndElement(xp, "exit_status")
     END IF
     IF (obj%cputime_ispresent) THEN
        CALL xml_NewElement(xp, "cputime")
           CALL xml_addCharacters(xp, obj%cputime)
        CALL xml_EndElement(xp, "cputime")
     END IF
     IF (obj%timing_info_ispresent) THEN
        CALL qes_write_timing (xp, obj%timing_info)
     END IF
     IF (obj%closed_ispresent) THEN
        CALL qes_write_closed (xp, obj%closed)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_espresso

   SUBROUTINE qes_write_general_info(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(general_info_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_xml_format (xp, obj%xml_format)
     CALL qes_write_creator (xp, obj%creator)
     CALL qes_write_created (xp, obj%created)
     CALL xml_NewElement(xp, 'job')
        CALL xml_addCharacters(xp, TRIM(obj%job))
     CALL xml_EndElement(xp, 'job')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_general_info

   SUBROUTINE qes_write_parallel_info(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(parallel_info_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'nprocs')
        CALL xml_addCharacters(xp, obj%nprocs)
     CALL xml_EndElement(xp, 'nprocs')
     CALL xml_NewElement(xp, 'nthreads')
        CALL xml_addCharacters(xp, obj%nthreads)
     CALL xml_EndElement(xp, 'nthreads')
     CALL xml_NewElement(xp, 'ntasks')
        CALL xml_addCharacters(xp, obj%ntasks)
     CALL xml_EndElement(xp, 'ntasks')
     CALL xml_NewElement(xp, 'nbgrp')
        CALL xml_addCharacters(xp, obj%nbgrp)
     CALL xml_EndElement(xp, 'nbgrp')
     CALL xml_NewElement(xp, 'npool')
        CALL xml_addCharacters(xp, obj%npool)
     CALL xml_EndElement(xp, 'npool')
     CALL xml_NewElement(xp, 'ndiag')
        CALL xml_addCharacters(xp, obj%ndiag)
     CALL xml_EndElement(xp, 'ndiag')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_parallel_info

   SUBROUTINE qes_write_input(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(input_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_control_variables (xp, obj%control_variables)
     CALL qes_write_atomic_species (xp, obj%atomic_species)
     CALL qes_write_atomic_structure (xp, obj%atomic_structure)
     CALL qes_write_dft (xp, obj%dft)
     CALL qes_write_spin (xp, obj%spin)
     CALL qes_write_bands (xp, obj%bands)
     CALL qes_write_basis (xp, obj%basis)
     CALL qes_write_electron_control (xp, obj%electron_control)
     CALL qes_write_k_points_IBZ (xp, obj%k_points_IBZ)
     CALL qes_write_ion_control (xp, obj%ion_control)
     CALL qes_write_cell_control (xp, obj%cell_control)
     IF (obj%symmetry_flags_ispresent) THEN
        CALL qes_write_symmetry_flags (xp, obj%symmetry_flags)
     END IF
     IF (obj%boundary_conditions_ispresent) THEN
        CALL qes_write_boundary_conditions (xp, obj%boundary_conditions)
     END IF
     IF (obj%fcp_settings_ispresent) THEN
        CALL qes_write_fcp (xp, obj%fcp_settings)
     END IF
     IF (obj%rism_settings_ispresent) THEN
        CALL qes_write_rism (xp, obj%rism_settings)
     END IF
     IF (obj%solvents_ispresent) THEN
        CALL qes_write_solvents (xp, obj%solvents)
     END IF
     IF (obj%ekin_functional_ispresent) THEN
        CALL qes_write_ekin_functional (xp, obj%ekin_functional)
     END IF
     IF (obj%external_atomic_forces_ispresent) THEN
        CALL qes_write_matrix (xp, obj%external_atomic_forces)
     END IF
     IF (obj%free_positions_ispresent) THEN
        CALL qes_write_integerMatrix (xp, obj%free_positions)
     END IF
     IF (obj%starting_atomic_velocities_ispresent) THEN
        CALL qes_write_matrix (xp, obj%starting_atomic_velocities)
     END IF
     IF (obj%electric_field_ispresent) THEN
        CALL qes_write_electric_field (xp, obj%electric_field)
     END IF
     IF (obj%atomic_constraints_ispresent) THEN
        CALL qes_write_atomic_constraints (xp, obj%atomic_constraints)
     END IF
     IF (obj%spin_constraints_ispresent) THEN
        CALL qes_write_spin_constraints (xp, obj%spin_constraints)
     END IF
     IF (obj%twoch__ispresent) THEN
        CALL qes_write_two_chem (xp, obj%twoch_)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_input

   SUBROUTINE qes_write_step(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(step_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%n_step_ispresent) CALL xml_addAttribute(xp, 'n_step', obj%n_step )
     CALL qes_write_scf_conv (xp, obj%scf_conv)
     CALL qes_write_atomic_structure (xp, obj%atomic_structure)
     CALL qes_write_total_energy (xp, obj%total_energy)
     CALL qes_write_matrix (xp, obj%forces)
     IF (obj%stress_ispresent) THEN
        CALL qes_write_matrix (xp, obj%stress)
     END IF
     IF (obj%fcp_force_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_force")
           CALL xml_addCharacters(xp, obj%fcp_force, fmt='s16')
        CALL xml_EndElement(xp, "fcp_force")
     END IF
     IF (obj%fcp_tot_charge_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_tot_charge")
           CALL xml_addCharacters(xp, obj%fcp_tot_charge, fmt='s16')
        CALL xml_EndElement(xp, "fcp_tot_charge")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_step

   SUBROUTINE qes_write_output(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(output_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%convergence_info_ispresent) THEN
        CALL qes_write_convergence_info (xp, obj%convergence_info)
     END IF
     CALL qes_write_algorithmic_info (xp, obj%algorithmic_info)
     CALL qes_write_atomic_species (xp, obj%atomic_species)
     CALL qes_write_atomic_structure (xp, obj%atomic_structure)
     IF (obj%symmetries_ispresent) THEN
        CALL qes_write_symmetries (xp, obj%symmetries)
     END IF
     CALL qes_write_basis_set (xp, obj%basis_set)
     CALL qes_write_dft (xp, obj%dft)
     IF (obj%boundary_conditions_ispresent) THEN
        CALL qes_write_outputPBC (xp, obj%boundary_conditions)
     END IF
     IF (obj%magnetization_ispresent) THEN
        CALL qes_write_magnetization (xp, obj%magnetization)
     END IF
     CALL qes_write_total_energy (xp, obj%total_energy)
     CALL qes_write_band_structure (xp, obj%band_structure)
     IF (obj%forces_ispresent) THEN
        CALL qes_write_matrix (xp, obj%forces)
     END IF
     IF (obj%stress_ispresent) THEN
        CALL qes_write_matrix (xp, obj%stress)
     END IF
     IF (obj%electric_field_ispresent) THEN
        CALL qes_write_outputElectricField (xp, obj%electric_field)
     END IF
     IF (obj%fcp_force_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_force")
           CALL xml_addCharacters(xp, obj%fcp_force, fmt='s16')
        CALL xml_EndElement(xp, "fcp_force")
     END IF
     IF (obj%fcp_tot_charge_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_tot_charge")
           CALL xml_addCharacters(xp, obj%fcp_tot_charge, fmt='s16')
        CALL xml_EndElement(xp, "fcp_tot_charge")
     END IF
     IF (obj%rism3d_ispresent) THEN
        CALL qes_write_rism3d (xp, obj%rism3d)
     END IF
     IF (obj%rismlaue_ispresent) THEN
        CALL qes_write_rismlaue (xp, obj%rismlaue)
     END IF
     IF (obj%two_chem_ispresent) THEN
        CALL qes_write_two_chem (xp, obj%two_chem)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_output

   SUBROUTINE qes_write_timing(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(timing_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_clock (xp, obj%total)
     IF (obj%partial_ispresent) THEN
        DO i = 1, obj%ndim_partial
           CALL qes_write_clock(xp, obj%partial(i) )
        END DO
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_timing

   SUBROUTINE qes_write_clock(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(clock_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
     IF (obj%calls_ispresent) CALL xml_addAttribute(xp, 'calls', obj%calls )
     CALL xml_NewElement(xp, 'cpu')
        CALL xml_addCharacters(xp, obj%cpu, fmt='s16')
     CALL xml_EndElement(xp, 'cpu')
     CALL xml_NewElement(xp, 'wall')
        CALL xml_addCharacters(xp, obj%wall, fmt='s16')
     CALL xml_EndElement(xp, 'wall')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_clock

   SUBROUTINE qes_write_control_variables(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(control_variables_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'title')
        CALL xml_addCharacters(xp, TRIM(obj%title))
     CALL xml_EndElement(xp, 'title')
     CALL xml_NewElement(xp, 'calculation')
        CALL xml_addCharacters(xp, TRIM(obj%calculation))
     CALL xml_EndElement(xp, 'calculation')
     CALL xml_NewElement(xp, 'restart_mode')
        CALL xml_addCharacters(xp, TRIM(obj%restart_mode))
     CALL xml_EndElement(xp, 'restart_mode')
     CALL xml_NewElement(xp, 'prefix')
        CALL xml_addCharacters(xp, TRIM(obj%prefix))
     CALL xml_EndElement(xp, 'prefix')
     CALL xml_NewElement(xp, 'pseudo_dir')
        CALL xml_addCharacters(xp, TRIM(obj%pseudo_dir))
     CALL xml_EndElement(xp, 'pseudo_dir')
     CALL xml_NewElement(xp, 'outdir')
        CALL xml_addCharacters(xp, TRIM(obj%outdir))
     CALL xml_EndElement(xp, 'outdir')
     CALL xml_NewElement(xp, 'stress')
        CALL xml_addCharacters(xp, obj%stress)
     CALL xml_EndElement(xp, 'stress')
     CALL xml_NewElement(xp, 'forces')
        CALL xml_addCharacters(xp, obj%forces)
     CALL xml_EndElement(xp, 'forces')
     CALL xml_NewElement(xp, 'wf_collect')
        CALL xml_addCharacters(xp, obj%wf_collect)
     CALL xml_EndElement(xp, 'wf_collect')
     CALL xml_NewElement(xp, 'disk_io')
        CALL xml_addCharacters(xp, TRIM(obj%disk_io))
     CALL xml_EndElement(xp, 'disk_io')
     CALL xml_NewElement(xp, 'max_seconds')
        CALL xml_addCharacters(xp, obj%max_seconds)
     CALL xml_EndElement(xp, 'max_seconds')
     IF (obj%nstep_ispresent) THEN
        CALL xml_NewElement(xp, "nstep")
           CALL xml_addCharacters(xp, obj%nstep)
        CALL xml_EndElement(xp, "nstep")
     END IF
     CALL xml_NewElement(xp, 'etot_conv_thr')
        CALL xml_addCharacters(xp, obj%etot_conv_thr, fmt='s16')
     CALL xml_EndElement(xp, 'etot_conv_thr')
     CALL xml_NewElement(xp, 'forc_conv_thr')
        CALL xml_addCharacters(xp, obj%forc_conv_thr, fmt='s16')
     CALL xml_EndElement(xp, 'forc_conv_thr')
     CALL xml_NewElement(xp, 'press_conv_thr')
        CALL xml_addCharacters(xp, obj%press_conv_thr, fmt='s16')
     CALL xml_EndElement(xp, 'press_conv_thr')
     CALL xml_NewElement(xp, 'verbosity')
        CALL xml_addCharacters(xp, TRIM(obj%verbosity))
     CALL xml_EndElement(xp, 'verbosity')
     CALL xml_NewElement(xp, 'print_every')
        CALL xml_addCharacters(xp, obj%print_every)
     CALL xml_EndElement(xp, 'print_every')
     CALL xml_NewElement(xp, 'fcp')
        CALL xml_addCharacters(xp, obj%fcp)
     CALL xml_EndElement(xp, 'fcp')
     CALL xml_NewElement(xp, 'rism')
        CALL xml_addCharacters(xp, obj%rism)
     CALL xml_EndElement(xp, 'rism')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_control_variables

   SUBROUTINE qes_write_xml_format(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(xml_format_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%NAME_ispresent) CALL xml_addAttribute(xp, 'NAME', TRIM(obj%NAME) )
     IF (obj%VERSION_ispresent) CALL xml_addAttribute(xp, 'VERSION', TRIM(obj%VERSION) )
        CALL xml_AddCharacters(xp, TRIM(obj%xml_format))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_xml_format

   SUBROUTINE qes_write_creator(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(creator_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%NAME_ispresent) CALL xml_addAttribute(xp, 'NAME', TRIM(obj%NAME) )
     IF (obj%VERSION_ispresent) CALL xml_addAttribute(xp, 'VERSION', TRIM(obj%VERSION) )
        CALL xml_AddCharacters(xp, TRIM(obj%creator))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_creator

   SUBROUTINE qes_write_created(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(created_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%DATE_ispresent) CALL xml_addAttribute(xp, 'DATE', TRIM(obj%DATE) )
     IF (obj%TIME_ispresent) CALL xml_addAttribute(xp, 'TIME', TRIM(obj%TIME) )
        CALL xml_AddCharacters(xp, TRIM(obj%created))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_created

   SUBROUTINE qes_write_atomic_species(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(atomic_species_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%ntyp_ispresent) CALL xml_addAttribute(xp, 'ntyp', obj%ntyp )
     IF (obj%pseudo_dir_ispresent) CALL xml_addAttribute(xp, 'pseudo_dir', TRIM(obj%pseudo_dir) )
     DO i = 1, obj%ndim_species
        CALL qes_write_species(xp, obj%species(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_atomic_species

   SUBROUTINE qes_write_species(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(species_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%name_ispresent) CALL xml_addAttribute(xp, 'name', TRIM(obj%name) )
     IF (obj%mass_ispresent) THEN
        CALL xml_NewElement(xp, "mass")
           CALL xml_addCharacters(xp, obj%mass, fmt='s16')
        CALL xml_EndElement(xp, "mass")
     END IF
     CALL xml_NewElement(xp, 'pseudo_file')
        CALL xml_addCharacters(xp, TRIM(obj%pseudo_file))
     CALL xml_EndElement(xp, 'pseudo_file')
     IF (obj%starting_magnetization_ispresent) THEN
        CALL xml_NewElement(xp, "starting_magnetization")
           CALL xml_addCharacters(xp, obj%starting_magnetization, fmt='s16')
        CALL xml_EndElement(xp, "starting_magnetization")
     END IF
     IF (obj%spin_teta_ispresent) THEN
        CALL xml_NewElement(xp, "spin_teta")
           CALL xml_addCharacters(xp, obj%spin_teta, fmt='s16')
        CALL xml_EndElement(xp, "spin_teta")
     END IF
     IF (obj%spin_phi_ispresent) THEN
        CALL xml_NewElement(xp, "spin_phi")
           CALL xml_addCharacters(xp, obj%spin_phi, fmt='s16')
        CALL xml_EndElement(xp, "spin_phi")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_species

   SUBROUTINE qes_write_atomic_structure(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(atomic_structure_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%nat_ispresent) CALL xml_addAttribute(xp, 'nat', obj%nat )
     IF (obj%num_of_atomic_wfc_ispresent) CALL xml_addAttribute(xp, 'num_of_atomic_wfc', obj%num_of_atomic_wfc )
     IF (obj%alat_ispresent) CALL xml_addAttribute(xp, 'alat', obj%alat )
     IF (obj%bravais_index_ispresent) CALL xml_addAttribute(xp, 'bravais_index', obj%bravais_index )
     IF (obj%alternative_axes_ispresent) CALL xml_addAttribute(xp, 'alternative_axes', TRIM(obj%alternative_axes) )
     IF (obj%atomic_positions_ispresent) THEN
        CALL qes_write_atomic_positions (xp, obj%atomic_positions)
     END IF
     IF (obj%wyckoff_positions_ispresent) THEN
        CALL qes_write_wyckoff_positions (xp, obj%wyckoff_positions)
     END IF
     IF (obj%crystal_positions_ispresent) THEN
        CALL qes_write_atomic_positions (xp, obj%crystal_positions)
     END IF
     CALL qes_write_cell (xp, obj%cell)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_atomic_structure

   SUBROUTINE qes_write_atomic_positions(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(atomic_positions_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     DO i = 1, obj%ndim_atom
        CALL qes_write_atom(xp, obj%atom(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_atomic_positions

   SUBROUTINE qes_write_atom(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(atom_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%name_ispresent) CALL xml_addAttribute(xp, 'name', TRIM(obj%name) )
     IF (obj%position_ispresent) CALL xml_addAttribute(xp, 'position', TRIM(obj%position) )
     IF (obj%index_ispresent) CALL xml_addAttribute(xp, 'index', obj%index )
        CALL xml_AddCharacters(xp, obj%atom, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_atom

   SUBROUTINE qes_write_wyckoff_positions(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(wyckoff_positions_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%space_group_ispresent) CALL xml_addAttribute(xp, 'space_group', obj%space_group )
     IF (obj%more_options_ispresent) CALL xml_addAttribute(xp, 'more_options', TRIM(obj%more_options) )
     DO i = 1, obj%ndim_atom
        CALL qes_write_atom(xp, obj%atom(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_wyckoff_positions

   SUBROUTINE qes_write_cell(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(cell_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'a1')
        CALL xml_addCharacters(xp, obj%a1, fmt='s16')
     CALL xml_EndElement(xp, 'a1')
     CALL xml_NewElement(xp, 'a2')
        CALL xml_addCharacters(xp, obj%a2, fmt='s16')
     CALL xml_EndElement(xp, 'a2')
     CALL xml_NewElement(xp, 'a3')
        CALL xml_addCharacters(xp, obj%a3, fmt='s16')
     CALL xml_EndElement(xp, 'a3')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_cell

   SUBROUTINE qes_write_dft(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(dft_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'functional')
        CALL xml_addCharacters(xp, TRIM(obj%functional))
     CALL xml_EndElement(xp, 'functional')
     IF (obj%hybrid_ispresent) THEN
        CALL qes_write_hybrid (xp, obj%hybrid)
     END IF
     IF (obj%dftU_ispresent) THEN
        CALL qes_write_dftU (xp, obj%dftU)
     END IF
     IF (obj%vdW_ispresent) THEN
        CALL qes_write_vdW (xp, obj%vdW)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_dft

   SUBROUTINE qes_write_hybrid(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(hybrid_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%qpoint_grid_ispresent) THEN
        CALL qes_write_qpoint_grid (xp, obj%qpoint_grid)
     END IF
     IF (obj%ecutfock_ispresent) THEN
        CALL xml_NewElement(xp, "ecutfock")
           CALL xml_addCharacters(xp, obj%ecutfock, fmt='s16')
        CALL xml_EndElement(xp, "ecutfock")
     END IF
     IF (obj%exx_fraction_ispresent) THEN
        CALL xml_NewElement(xp, "exx_fraction")
           CALL xml_addCharacters(xp, obj%exx_fraction, fmt='s16')
        CALL xml_EndElement(xp, "exx_fraction")
     END IF
     IF (obj%screening_parameter_ispresent) THEN
        CALL xml_NewElement(xp, "screening_parameter")
           CALL xml_addCharacters(xp, obj%screening_parameter, fmt='s16')
        CALL xml_EndElement(xp, "screening_parameter")
     END IF
     IF (obj%exxdiv_treatment_ispresent) THEN
        CALL xml_NewElement(xp, "exxdiv_treatment")
           CALL xml_addCharacters(xp, TRIM(obj%exxdiv_treatment))
        CALL xml_EndElement(xp, "exxdiv_treatment")
     END IF
     IF (obj%x_gamma_extrapolation_ispresent) THEN
        CALL xml_NewElement(xp, "x_gamma_extrapolation")
           CALL xml_addCharacters(xp, obj%x_gamma_extrapolation)
        CALL xml_EndElement(xp, "x_gamma_extrapolation")
     END IF
     IF (obj%ecutvcut_ispresent) THEN
        CALL xml_NewElement(xp, "ecutvcut")
           CALL xml_addCharacters(xp, obj%ecutvcut, fmt='s16')
        CALL xml_EndElement(xp, "ecutvcut")
     END IF
     IF (obj%localization_threshold_ispresent) THEN
        CALL xml_NewElement(xp, "localization_threshold")
           CALL xml_addCharacters(xp, obj%localization_threshold, fmt='s16')
        CALL xml_EndElement(xp, "localization_threshold")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_hybrid

   SUBROUTINE qes_write_qpoint_grid(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(qpoint_grid_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%nqx1_ispresent) CALL xml_addAttribute(xp, 'nqx1', obj%nqx1 )
     IF (obj%nqx2_ispresent) CALL xml_addAttribute(xp, 'nqx2', obj%nqx2 )
     IF (obj%nqx3_ispresent) CALL xml_addAttribute(xp, 'nqx3', obj%nqx3 )
        CALL xml_AddCharacters(xp, TRIM(obj%qpoint_grid))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_qpoint_grid

   SUBROUTINE qes_write_dftU(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(dftU_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%new_format_ispresent) CALL xml_addAttribute(xp, 'new_format', obj%new_format )
     IF (obj%lda_plus_u_kind_ispresent) THEN
        CALL xml_NewElement(xp, "lda_plus_u_kind")
           CALL xml_addCharacters(xp, obj%lda_plus_u_kind)
        CALL xml_EndElement(xp, "lda_plus_u_kind")
     END IF
     IF (obj%Hubbard_Occ_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_Occ
           CALL qes_write_HubbardOcc(xp, obj%Hubbard_Occ(i) )
        END DO
     END IF
     IF (obj%Hubbard_U_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_U
           CALL qes_write_HubbardCommon(xp, obj%Hubbard_U(i) )
        END DO
     END IF
     IF (obj%Hubbard_Um_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_Um
           CALL qes_write_HubbardM(xp, obj%Hubbard_Um(i) )
        END DO
     END IF
     IF (obj%Hubbard_J0_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_J0
           CALL qes_write_HubbardCommon(xp, obj%Hubbard_J0(i) )
        END DO
     END IF
     IF (obj%Hubbard_alpha_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_alpha
           CALL qes_write_HubbardCommon(xp, obj%Hubbard_alpha(i) )
        END DO
     END IF
     IF (obj%Hubbard_beta_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_beta
           CALL qes_write_HubbardCommon(xp, obj%Hubbard_beta(i) )
        END DO
     END IF
     IF (obj%Hubbard_J_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_J
           CALL qes_write_HubbardJ(xp, obj%Hubbard_J(i) )
        END DO
     END IF
     IF (obj%starting_ns_ispresent) THEN
        DO i = 1, obj%ndim_starting_ns
           CALL qes_write_starting_ns(xp, obj%starting_ns(i) )
        END DO
     END IF
     IF (obj%Hubbard_V_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_V
           CALL qes_write_HubbardInterSpecieV(xp, obj%Hubbard_V(i) )
        END DO
     END IF
     IF (obj%Hubbard_ns_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_ns
           CALL qes_write_Hubbard_ns(xp, obj%Hubbard_ns(i) )
        END DO
     END IF
     IF (obj%Hub_m_order_ispresent) THEN
        DO i = 1, obj%ndim_Hub_m_order
           CALL qes_write_orderUm(xp, obj%Hub_m_order(i) )
        END DO
     END IF
     IF (obj%U_projection_type_ispresent) THEN
        CALL xml_NewElement(xp, "U_projection_type")
           CALL xml_addCharacters(xp, TRIM(obj%U_projection_type))
        CALL xml_EndElement(xp, "U_projection_type")
     END IF
     IF (obj%Hubbard_back_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_back
           CALL qes_write_HubbardBack(xp, obj%Hubbard_back(i) )
        END DO
     END IF
     IF (obj%Hubbard_alpha_back_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_alpha_back
           CALL qes_write_HubbardCommon(xp, obj%Hubbard_alpha_back(i) )
        END DO
     END IF
     IF (obj%Hubbard_ns_nc_ispresent) THEN
        DO i = 1, obj%ndim_Hubbard_ns_nc
           CALL qes_write_Hubbard_ns(xp, obj%Hubbard_ns_nc(i) )
        END DO
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_dftU

   SUBROUTINE qes_write_HubbardCommon(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(HubbardCommon_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%specie_ispresent) CALL xml_addAttribute(xp, 'specie', TRIM(obj%specie) )
     IF (obj%label_ispresent) CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
        CALL xml_AddCharacters(xp, obj%HubbardCommon, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_HubbardCommon

   SUBROUTINE qes_write_HubbardInterSpecieV(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(HubbardInterSpecieV_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'specie1', TRIM(obj%specie1) )
     CALL xml_addAttribute(xp, 'index1', obj%index1 )
     IF (obj%label1_ispresent) CALL xml_addAttribute(xp, 'label1', TRIM(obj%label1) )
     CALL xml_addAttribute(xp, 'specie2', TRIM(obj%specie2) )
     CALL xml_addAttribute(xp, 'index2', obj%index2 )
     IF (obj%label2_ispresent) CALL xml_addAttribute(xp, 'label2', TRIM(obj%label2) )
        CALL xml_AddCharacters(xp, obj%HubbardInterSpecieV, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_HubbardInterSpecieV

   SUBROUTINE qes_write_SiteMoment(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(SiteMoment_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%species_ispresent) CALL xml_addAttribute(xp, 'species', TRIM(obj%species) )
     IF (obj%atom_ispresent) CALL xml_addAttribute(xp, 'atom', obj%atom )
     IF (obj%charge_ispresent) CALL xml_addAttribute(xp, 'charge', obj%charge )
        CALL xml_AddCharacters(xp, obj%SiteMoment, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_SiteMoment

   SUBROUTINE qes_write_HubbardJ(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(HubbardJ_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%specie_ispresent) CALL xml_addAttribute(xp, 'specie', TRIM(obj%specie) )
     IF (obj%label_ispresent) CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
        CALL xml_AddCharacters(xp, obj%HubbardJ, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_HubbardJ

   SUBROUTINE qes_write_vector(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(vector_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'size', obj%size )
     CALL xml_addNewLine(xp)
     DO i = 1, obj%size, 5
        CALL xml_AddCharacters(xp, obj%vector(i:MIN(i+5-1,obj%size)), fmt='s16')
        CALL xml_AddNewLine(xp)
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_vector

   SUBROUTINE qes_write_HubbardM(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(HubbardM_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'size', obj%size )
     IF (obj%specie_ispresent) CALL xml_addAttribute(xp, 'specie', TRIM(obj%specie) )
     IF (obj%label_ispresent) CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
     IF (obj%spin_ispresent) CALL xml_addAttribute(xp, 'spin', obj%spin )
     IF (obj%jjj_ispresent) CALL xml_addAttribute(xp, 'jjj', obj%jjj )
     CALL xml_addNewLine(xp)
     DO i = 1, obj%size, 5
        CALL xml_AddCharacters(xp, obj%HubbardM(i:MIN(i+5-1,obj%size)), fmt='s16')
        CALL xml_AddNewLine(xp)
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_HubbardM

   SUBROUTINE qes_write_ChannelOcc(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(ChannelOcc_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%specie_ispresent) CALL xml_addAttribute(xp, 'specie', TRIM(obj%specie) )
     IF (obj%label_ispresent) CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
     CALL xml_addAttribute(xp, 'index', obj%index )
        CALL xml_AddCharacters(xp, obj%ChannelOcc, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_ChannelOcc

   SUBROUTINE qes_write_HubbardOcc(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(HubbardOcc_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'channels', obj%channels )
     CALL xml_addAttribute(xp, 'specie', TRIM(obj%specie) )
     DO i = 1, obj%ndim_channel_occ
        CALL qes_write_ChannelOcc(xp, obj%channel_occ(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_HubbardOcc

   SUBROUTINE qes_write_SitMag(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(SitMag_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%species_ispresent) CALL xml_addAttribute(xp, 'species', TRIM(obj%species) )
     IF (obj%atom_ispresent) CALL xml_addAttribute(xp, 'atom', obj%atom )
     IF (obj%charge_ispresent) CALL xml_addAttribute(xp, 'charge', obj%charge )
        CALL xml_AddCharacters(xp, obj%SitMag, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_SitMag

   SUBROUTINE qes_write_starting_ns(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(starting_ns_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'size', obj%size )
     IF (obj%specie_ispresent) CALL xml_addAttribute(xp, 'specie', TRIM(obj%specie) )
     IF (obj%label_ispresent) CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
     IF (obj%spin_ispresent) CALL xml_addAttribute(xp, 'spin', obj%spin )
     CALL xml_addNewLine(xp)
     DO i = 1, obj%size, 5
        CALL xml_AddCharacters(xp, obj%starting_ns(i:MIN(i+5-1,obj%size)), fmt='s16')
        CALL xml_AddNewLine(xp)
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_starting_ns

   SUBROUTINE qes_write_integerVector(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(integerVector_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'size', obj%size )
     CALL xml_addNewLine(xp)
     DO i = 1, obj%size, 8
        CALL xml_AddCharacters(xp, obj%integerVector(i:MIN(i+8-1, obj%size)))
        CALL xml_AddNewLine(xp)
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_integerVector

   SUBROUTINE qes_write_orderUm(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(orderUm_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'size', obj%size )
     IF (obj%specie_ispresent) CALL xml_addAttribute(xp, 'specie', TRIM(obj%specie) )
     IF (obj%label_ispresent) CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
     IF (obj%spin_ispresent) CALL xml_addAttribute(xp, 'spin', obj%spin )
     IF (obj%atomidx_ispresent) CALL xml_addAttribute(xp, 'atomidx', obj%atomidx )
     CALL xml_addNewLine(xp)
     DO i = 1, obj%size, 8
        CALL xml_AddCharacters(xp, obj%orderUm(i:MIN(i+8-1, obj%size)))
        CALL xml_AddNewLine(xp)
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_orderUm

   SUBROUTINE qes_write_matrix(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(matrix_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'rank', obj%rank )
     CALL xml_addAttribute(xp, 'dims', obj%dims )
     IF (obj%order_ispresent) CALL xml_addAttribute(xp, 'order', TRIM(obj%order) )
       CALL xml_addNewLine(xp)
        DO i = 1, obj%dims(2)
           CALL xml_AddCharacters(xp, obj%matrix((i-1)*obj%dims(1)+1: i*obj%dims(1)), fmt ='s16')
           CALL xml_addNewLine(xp)
        END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_matrix

   SUBROUTINE qes_write_Hubbard_ns(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(Hubbard_ns_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'rank', obj%rank )
     CALL xml_addAttribute(xp, 'dims', obj%dims )
     IF (obj%order_ispresent) CALL xml_addAttribute(xp, 'order', TRIM(obj%order) )
     IF (obj%specie_ispresent) CALL xml_addAttribute(xp, 'specie', TRIM(obj%specie) )
     IF (obj%label_ispresent) CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
     IF (obj%spin_ispresent) CALL xml_addAttribute(xp, 'spin', obj%spin )
     IF (obj%index_ispresent) CALL xml_addAttribute(xp, 'index', obj%index )
       CALL xml_addNewLine(xp)
        DO i = 1, obj%dims(2)
           CALL xml_AddCharacters(xp, obj%Hubbard_ns((i-1)*obj%dims(1)+1: i*obj%dims(1)), fmt ='s16')
           CALL xml_addNewLine(xp)
        END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_Hubbard_ns

   SUBROUTINE qes_write_HubbardBack(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(HubbardBack_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'background', TRIM(obj%background) )
     IF (obj%label_ispresent) CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
     IF (obj%species_ispresent) CALL xml_addAttribute(xp, 'species', TRIM(obj%species) )
     CALL xml_NewElement(xp, 'Hubbard_U2')
        CALL xml_addCharacters(xp, obj%Hubbard_U2, fmt='s16')
     CALL xml_EndElement(xp, 'Hubbard_U2')
     CALL xml_NewElement(xp, 'n2_number')
        CALL xml_addCharacters(xp, obj%n2_number)
     CALL xml_EndElement(xp, 'n2_number')
     CALL xml_NewElement(xp, 'l2_number')
        CALL xml_addCharacters(xp, obj%l2_number)
     CALL xml_EndElement(xp, 'l2_number')
     IF (obj%n3_number_ispresent) THEN
        CALL xml_NewElement(xp, "n3_number")
           CALL xml_addCharacters(xp, obj%n3_number)
        CALL xml_EndElement(xp, "n3_number")
     END IF
     IF (obj%l3_number_ispresent) THEN
        CALL xml_NewElement(xp, "l3_number")
           CALL xml_addCharacters(xp, obj%l3_number)
        CALL xml_EndElement(xp, "l3_number")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_HubbardBack

   SUBROUTINE qes_write_vdW(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(vdW_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%vdw_corr_ispresent) THEN
        CALL xml_NewElement(xp, "vdw_corr")
           CALL xml_addCharacters(xp, TRIM(obj%vdw_corr))
        CALL xml_EndElement(xp, "vdw_corr")
     END IF
     IF (obj%dftd3_version_ispresent) THEN
        CALL xml_NewElement(xp, "dftd3_version")
           CALL xml_addCharacters(xp, obj%dftd3_version)
        CALL xml_EndElement(xp, "dftd3_version")
     END IF
     IF (obj%dftd3_threebody_ispresent) THEN
        CALL xml_NewElement(xp, "dftd3_threebody")
           CALL xml_addCharacters(xp, obj%dftd3_threebody)
        CALL xml_EndElement(xp, "dftd3_threebody")
     END IF
     IF (obj%non_local_term_ispresent) THEN
        CALL xml_NewElement(xp, "non_local_term")
           CALL xml_addCharacters(xp, TRIM(obj%non_local_term))
        CALL xml_EndElement(xp, "non_local_term")
     END IF
     IF (obj%functional_ispresent) THEN
        CALL xml_NewElement(xp, "functional")
           CALL xml_addCharacters(xp, TRIM(obj%functional))
        CALL xml_EndElement(xp, "functional")
     END IF
     IF (obj%total_energy_term_ispresent) THEN
        CALL xml_NewElement(xp, "total_energy_term")
           CALL xml_addCharacters(xp, obj%total_energy_term, fmt='s16')
        CALL xml_EndElement(xp, "total_energy_term")
     END IF
     IF (obj%london_s6_ispresent) THEN
        CALL xml_NewElement(xp, "london_s6")
           CALL xml_addCharacters(xp, obj%london_s6, fmt='s16')
        CALL xml_EndElement(xp, "london_s6")
     END IF
     IF (obj%ts_vdw_econv_thr_ispresent) THEN
        CALL xml_NewElement(xp, "ts_vdw_econv_thr")
           CALL xml_addCharacters(xp, obj%ts_vdw_econv_thr, fmt='s16')
        CALL xml_EndElement(xp, "ts_vdw_econv_thr")
     END IF
     IF (obj%ts_vdw_isolated_ispresent) THEN
        CALL xml_NewElement(xp, "ts_vdw_isolated")
           CALL xml_addCharacters(xp, obj%ts_vdw_isolated)
        CALL xml_EndElement(xp, "ts_vdw_isolated")
     END IF
     IF (obj%london_rcut_ispresent) THEN
        CALL xml_NewElement(xp, "london_rcut")
           CALL xml_addCharacters(xp, obj%london_rcut, fmt='s16')
        CALL xml_EndElement(xp, "london_rcut")
     END IF
     IF (obj%xdm_a1_ispresent) THEN
        CALL xml_NewElement(xp, "xdm_a1")
           CALL xml_addCharacters(xp, obj%xdm_a1, fmt='s16')
        CALL xml_EndElement(xp, "xdm_a1")
     END IF
     IF (obj%xdm_a2_ispresent) THEN
        CALL xml_NewElement(xp, "xdm_a2")
           CALL xml_addCharacters(xp, obj%xdm_a2, fmt='s16')
        CALL xml_EndElement(xp, "xdm_a2")
     END IF
     IF (obj%london_c6_ispresent) THEN
        DO i = 1, obj%ndim_london_c6
           CALL qes_write_HubbardCommon(xp, obj%london_c6(i) )
        END DO
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_vdW

   SUBROUTINE qes_write_spin(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(spin_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'lsda')
        CALL xml_addCharacters(xp, obj%lsda)
     CALL xml_EndElement(xp, 'lsda')
     CALL xml_NewElement(xp, 'noncolin')
        CALL xml_addCharacters(xp, obj%noncolin)
     CALL xml_EndElement(xp, 'noncolin')
     CALL xml_NewElement(xp, 'spinorbit')
        CALL xml_addCharacters(xp, obj%spinorbit)
     CALL xml_EndElement(xp, 'spinorbit')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_spin

   SUBROUTINE qes_write_bands(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(bands_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%nbnd_ispresent) THEN
        CALL xml_NewElement(xp, "nbnd")
           CALL xml_addCharacters(xp, obj%nbnd)
        CALL xml_EndElement(xp, "nbnd")
     END IF
     IF (obj%smearing_ispresent) THEN
        CALL qes_write_smearing (xp, obj%smearing)
     END IF
     IF (obj%tot_charge_ispresent) THEN
        CALL xml_NewElement(xp, "tot_charge")
           CALL xml_addCharacters(xp, obj%tot_charge, fmt='s16')
        CALL xml_EndElement(xp, "tot_charge")
     END IF
     IF (obj%tot_magnetization_ispresent) THEN
        CALL xml_NewElement(xp, "tot_magnetization")
           CALL xml_addCharacters(xp, obj%tot_magnetization, fmt='s16')
        CALL xml_EndElement(xp, "tot_magnetization")
     END IF
     CALL qes_write_occupations (xp, obj%occupations)
     IF (obj%inputOccupations_ispresent) THEN
        DO i = 1, obj%ndim_inputOccupations
           CALL qes_write_inputOccupations(xp, obj%inputOccupations(i) )
        END DO
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_bands

   SUBROUTINE qes_write_smearing(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(smearing_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%degauss_ispresent) CALL xml_addAttribute(xp, 'degauss', obj%degauss )
        CALL xml_AddCharacters(xp, TRIM(obj%smearing))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_smearing

   SUBROUTINE qes_write_occupations(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(occupations_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%spin_ispresent) CALL xml_addAttribute(xp, 'spin', obj%spin )
        CALL xml_AddCharacters(xp, TRIM(obj%occupations))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_occupations

   SUBROUTINE qes_write_basis(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(basis_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%gamma_only_ispresent) THEN
        CALL xml_NewElement(xp, "gamma_only")
           CALL xml_addCharacters(xp, obj%gamma_only)
        CALL xml_EndElement(xp, "gamma_only")
     END IF
     CALL xml_NewElement(xp, 'ecutwfc')
        CALL xml_addCharacters(xp, obj%ecutwfc, fmt='s16')
     CALL xml_EndElement(xp, 'ecutwfc')
     IF (obj%ecutrho_ispresent) THEN
        CALL xml_NewElement(xp, "ecutrho")
           CALL xml_addCharacters(xp, obj%ecutrho, fmt='s16')
        CALL xml_EndElement(xp, "ecutrho")
     END IF
     IF (obj%fft_grid_ispresent) THEN
        CALL qes_write_basisSetItem (xp, obj%fft_grid)
     END IF
     IF (obj%fft_smooth_ispresent) THEN
        CALL qes_write_basisSetItem (xp, obj%fft_smooth)
     END IF
     IF (obj%fft_box_ispresent) THEN
        CALL qes_write_basisSetItem (xp, obj%fft_box)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_basis

   SUBROUTINE qes_write_basis_set(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(basis_set_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%gamma_only_ispresent) THEN
        CALL xml_NewElement(xp, "gamma_only")
           CALL xml_addCharacters(xp, obj%gamma_only)
        CALL xml_EndElement(xp, "gamma_only")
     END IF
     CALL xml_NewElement(xp, 'ecutwfc')
        CALL xml_addCharacters(xp, obj%ecutwfc, fmt='s16')
     CALL xml_EndElement(xp, 'ecutwfc')
     IF (obj%ecutrho_ispresent) THEN
        CALL xml_NewElement(xp, "ecutrho")
           CALL xml_addCharacters(xp, obj%ecutrho, fmt='s16')
        CALL xml_EndElement(xp, "ecutrho")
     END IF
     CALL qes_write_basisSetItem (xp, obj%fft_grid)
     IF (obj%fft_smooth_ispresent) THEN
        CALL qes_write_basisSetItem (xp, obj%fft_smooth)
     END IF
     IF (obj%fft_box_ispresent) THEN
        CALL qes_write_basisSetItem (xp, obj%fft_box)
     END IF
     CALL xml_NewElement(xp, 'ngm')
        CALL xml_addCharacters(xp, obj%ngm)
     CALL xml_EndElement(xp, 'ngm')
     IF (obj%ngms_ispresent) THEN
        CALL xml_NewElement(xp, "ngms")
           CALL xml_addCharacters(xp, obj%ngms)
        CALL xml_EndElement(xp, "ngms")
     END IF
     CALL xml_NewElement(xp, 'npwx')
        CALL xml_addCharacters(xp, obj%npwx)
     CALL xml_EndElement(xp, 'npwx')
     CALL qes_write_reciprocal_lattice (xp, obj%reciprocal_lattice)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_basis_set

   SUBROUTINE qes_write_basisSetItem(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(basisSetItem_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%nr1_ispresent) CALL xml_addAttribute(xp, 'nr1', obj%nr1 )
     IF (obj%nr2_ispresent) CALL xml_addAttribute(xp, 'nr2', obj%nr2 )
     IF (obj%nr3_ispresent) CALL xml_addAttribute(xp, 'nr3', obj%nr3 )
        CALL xml_AddCharacters(xp, TRIM(obj%basisSetItem))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_basisSetItem

   SUBROUTINE qes_write_reciprocal_lattice(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(reciprocal_lattice_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'b1')
        CALL xml_addCharacters(xp, obj%b1, fmt='s16')
     CALL xml_EndElement(xp, 'b1')
     CALL xml_NewElement(xp, 'b2')
        CALL xml_addCharacters(xp, obj%b2, fmt='s16')
     CALL xml_EndElement(xp, 'b2')
     CALL xml_NewElement(xp, 'b3')
        CALL xml_addCharacters(xp, obj%b3, fmt='s16')
     CALL xml_EndElement(xp, 'b3')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_reciprocal_lattice

   SUBROUTINE qes_write_electron_control(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(electron_control_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'diagonalization')
        CALL xml_addCharacters(xp, TRIM(obj%diagonalization))
     CALL xml_EndElement(xp, 'diagonalization')
     CALL xml_NewElement(xp, 'mixing_mode')
        CALL xml_addCharacters(xp, TRIM(obj%mixing_mode))
     CALL xml_EndElement(xp, 'mixing_mode')
     CALL xml_NewElement(xp, 'mixing_beta')
        CALL xml_addCharacters(xp, obj%mixing_beta, fmt='s16')
     CALL xml_EndElement(xp, 'mixing_beta')
     CALL xml_NewElement(xp, 'conv_thr')
        CALL xml_addCharacters(xp, obj%conv_thr, fmt='s16')
     CALL xml_EndElement(xp, 'conv_thr')
     CALL xml_NewElement(xp, 'mixing_ndim')
        CALL xml_addCharacters(xp, obj%mixing_ndim)
     CALL xml_EndElement(xp, 'mixing_ndim')
     CALL xml_NewElement(xp, 'max_nstep')
        CALL xml_addCharacters(xp, obj%max_nstep)
     CALL xml_EndElement(xp, 'max_nstep')
     IF (obj%exx_nstep_ispresent) THEN
        CALL xml_NewElement(xp, "exx_nstep")
           CALL xml_addCharacters(xp, obj%exx_nstep)
        CALL xml_EndElement(xp, "exx_nstep")
     END IF
     IF (obj%real_space_q_ispresent) THEN
        CALL xml_NewElement(xp, "real_space_q")
           CALL xml_addCharacters(xp, obj%real_space_q)
        CALL xml_EndElement(xp, "real_space_q")
     END IF
     IF (obj%real_space_beta_ispresent) THEN
        CALL xml_NewElement(xp, "real_space_beta")
           CALL xml_addCharacters(xp, obj%real_space_beta)
        CALL xml_EndElement(xp, "real_space_beta")
     END IF
     CALL xml_NewElement(xp, 'tq_smoothing')
        CALL xml_addCharacters(xp, obj%tq_smoothing)
     CALL xml_EndElement(xp, 'tq_smoothing')
     CALL xml_NewElement(xp, 'tbeta_smoothing')
        CALL xml_addCharacters(xp, obj%tbeta_smoothing)
     CALL xml_EndElement(xp, 'tbeta_smoothing')
     CALL xml_NewElement(xp, 'diago_thr_init')
        CALL xml_addCharacters(xp, obj%diago_thr_init, fmt='s16')
     CALL xml_EndElement(xp, 'diago_thr_init')
     CALL xml_NewElement(xp, 'diago_full_acc')
        CALL xml_addCharacters(xp, obj%diago_full_acc)
     CALL xml_EndElement(xp, 'diago_full_acc')
     IF (obj%diago_cg_maxiter_ispresent) THEN
        CALL xml_NewElement(xp, "diago_cg_maxiter")
           CALL xml_addCharacters(xp, obj%diago_cg_maxiter)
        CALL xml_EndElement(xp, "diago_cg_maxiter")
     END IF
     IF (obj%diago_ppcg_maxiter_ispresent) THEN
        CALL xml_NewElement(xp, "diago_ppcg_maxiter")
           CALL xml_addCharacters(xp, obj%diago_ppcg_maxiter)
        CALL xml_EndElement(xp, "diago_ppcg_maxiter")
     END IF
     IF (obj%diago_david_ndim_ispresent) THEN
        CALL xml_NewElement(xp, "diago_david_ndim")
           CALL xml_addCharacters(xp, obj%diago_david_ndim)
        CALL xml_EndElement(xp, "diago_david_ndim")
     END IF
     IF (obj%diago_rmm_ndim_ispresent) THEN
        CALL xml_NewElement(xp, "diago_rmm_ndim")
           CALL xml_addCharacters(xp, obj%diago_rmm_ndim)
        CALL xml_EndElement(xp, "diago_rmm_ndim")
     END IF
     IF (obj%diago_gs_nblock_ispresent) THEN
        CALL xml_NewElement(xp, "diago_gs_nblock")
           CALL xml_addCharacters(xp, obj%diago_gs_nblock)
        CALL xml_EndElement(xp, "diago_gs_nblock")
     END IF
     IF (obj%diago_rmm_conv_ispresent) THEN
        CALL xml_NewElement(xp, "diago_rmm_conv")
           CALL xml_addCharacters(xp, obj%diago_rmm_conv)
        CALL xml_EndElement(xp, "diago_rmm_conv")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_electron_control

   SUBROUTINE qes_write_fcp(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(fcp_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%fcp_mu_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_mu")
           CALL xml_addCharacters(xp, obj%fcp_mu, fmt='s16')
        CALL xml_EndElement(xp, "fcp_mu")
     END IF
     IF (obj%fcp_dynamics_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_dynamics")
           CALL xml_addCharacters(xp, TRIM(obj%fcp_dynamics))
        CALL xml_EndElement(xp, "fcp_dynamics")
     END IF
     IF (obj%fcp_conv_thr_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_conv_thr")
           CALL xml_addCharacters(xp, obj%fcp_conv_thr, fmt='s16')
        CALL xml_EndElement(xp, "fcp_conv_thr")
     END IF
     IF (obj%fcp_ndiis_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_ndiis")
           CALL xml_addCharacters(xp, obj%fcp_ndiis)
        CALL xml_EndElement(xp, "fcp_ndiis")
     END IF
     IF (obj%fcp_rdiis_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_rdiis")
           CALL xml_addCharacters(xp, obj%fcp_rdiis, fmt='s16')
        CALL xml_EndElement(xp, "fcp_rdiis")
     END IF
     IF (obj%fcp_mass_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_mass")
           CALL xml_addCharacters(xp, obj%fcp_mass, fmt='s16')
        CALL xml_EndElement(xp, "fcp_mass")
     END IF
     IF (obj%fcp_velocity_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_velocity")
           CALL xml_addCharacters(xp, obj%fcp_velocity, fmt='s16')
        CALL xml_EndElement(xp, "fcp_velocity")
     END IF
     IF (obj%fcp_temperature_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_temperature")
           CALL xml_addCharacters(xp, TRIM(obj%fcp_temperature))
        CALL xml_EndElement(xp, "fcp_temperature")
     END IF
     IF (obj%fcp_tempw_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_tempw")
           CALL xml_addCharacters(xp, obj%fcp_tempw, fmt='s16')
        CALL xml_EndElement(xp, "fcp_tempw")
     END IF
     IF (obj%fcp_tolp_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_tolp")
           CALL xml_addCharacters(xp, obj%fcp_tolp, fmt='s16')
        CALL xml_EndElement(xp, "fcp_tolp")
     END IF
     IF (obj%fcp_delta_t_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_delta_t")
           CALL xml_addCharacters(xp, obj%fcp_delta_t, fmt='s16')
        CALL xml_EndElement(xp, "fcp_delta_t")
     END IF
     IF (obj%fcp_nraise_ispresent) THEN
        CALL xml_NewElement(xp, "fcp_nraise")
           CALL xml_addCharacters(xp, obj%fcp_nraise)
        CALL xml_EndElement(xp, "fcp_nraise")
     END IF
     IF (obj%freeze_all_atoms_ispresent) THEN
        CALL xml_NewElement(xp, "freeze_all_atoms")
           CALL xml_addCharacters(xp, obj%freeze_all_atoms)
        CALL xml_EndElement(xp, "freeze_all_atoms")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_fcp

   SUBROUTINE qes_write_rism(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(rism_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'nsolv')
        CALL xml_addCharacters(xp, obj%nsolv)
     CALL xml_EndElement(xp, 'nsolv')
     DO i = 1, obj%ndim_solute
        CALL qes_write_solute(xp, obj%solute(i) )
     END DO
     IF (obj%closure_ispresent) THEN
        CALL xml_NewElement(xp, "closure")
           CALL xml_addCharacters(xp, TRIM(obj%closure))
        CALL xml_EndElement(xp, "closure")
     END IF
     IF (obj%tempv_ispresent) THEN
        CALL xml_NewElement(xp, "tempv")
           CALL xml_addCharacters(xp, obj%tempv, fmt='s16')
        CALL xml_EndElement(xp, "tempv")
     END IF
     IF (obj%ecutsolv_ispresent) THEN
        CALL xml_NewElement(xp, "ecutsolv")
           CALL xml_addCharacters(xp, obj%ecutsolv, fmt='s16')
        CALL xml_EndElement(xp, "ecutsolv")
     END IF
     IF (obj%rmax_lj_ispresent) THEN
        CALL xml_NewElement(xp, "rmax_lj")
           CALL xml_addCharacters(xp, obj%rmax_lj, fmt='s16')
        CALL xml_EndElement(xp, "rmax_lj")
     END IF
     IF (obj%rmax1d_ispresent) THEN
        CALL xml_NewElement(xp, "rmax1d")
           CALL xml_addCharacters(xp, obj%rmax1d, fmt='s16')
        CALL xml_EndElement(xp, "rmax1d")
     END IF
     IF (obj%starting1d_ispresent) THEN
        CALL xml_NewElement(xp, "starting1d")
           CALL xml_addCharacters(xp, TRIM(obj%starting1d))
        CALL xml_EndElement(xp, "starting1d")
     END IF
     IF (obj%starting3d_ispresent) THEN
        CALL xml_NewElement(xp, "starting3d")
           CALL xml_addCharacters(xp, TRIM(obj%starting3d))
        CALL xml_EndElement(xp, "starting3d")
     END IF
     IF (obj%smear1d_ispresent) THEN
        CALL xml_NewElement(xp, "smear1d")
           CALL xml_addCharacters(xp, obj%smear1d, fmt='s16')
        CALL xml_EndElement(xp, "smear1d")
     END IF
     IF (obj%smear3d_ispresent) THEN
        CALL xml_NewElement(xp, "smear3d")
           CALL xml_addCharacters(xp, obj%smear3d, fmt='s16')
        CALL xml_EndElement(xp, "smear3d")
     END IF
     IF (obj%rism1d_maxstep_ispresent) THEN
        CALL xml_NewElement(xp, "rism1d_maxstep")
           CALL xml_addCharacters(xp, obj%rism1d_maxstep)
        CALL xml_EndElement(xp, "rism1d_maxstep")
     END IF
     IF (obj%rism3d_maxstep_ispresent) THEN
        CALL xml_NewElement(xp, "rism3d_maxstep")
           CALL xml_addCharacters(xp, obj%rism3d_maxstep)
        CALL xml_EndElement(xp, "rism3d_maxstep")
     END IF
     IF (obj%rism1d_conv_thr_ispresent) THEN
        CALL xml_NewElement(xp, "rism1d_conv_thr")
           CALL xml_addCharacters(xp, obj%rism1d_conv_thr, fmt='s16')
        CALL xml_EndElement(xp, "rism1d_conv_thr")
     END IF
     IF (obj%rism3d_conv_thr_ispresent) THEN
        CALL xml_NewElement(xp, "rism3d_conv_thr")
           CALL xml_addCharacters(xp, obj%rism3d_conv_thr, fmt='s16')
        CALL xml_EndElement(xp, "rism3d_conv_thr")
     END IF
     IF (obj%mdiis1d_size_ispresent) THEN
        CALL xml_NewElement(xp, "mdiis1d_size")
           CALL xml_addCharacters(xp, obj%mdiis1d_size)
        CALL xml_EndElement(xp, "mdiis1d_size")
     END IF
     IF (obj%mdiis3d_size_ispresent) THEN
        CALL xml_NewElement(xp, "mdiis3d_size")
           CALL xml_addCharacters(xp, obj%mdiis3d_size)
        CALL xml_EndElement(xp, "mdiis3d_size")
     END IF
     IF (obj%mdiis1d_step_ispresent) THEN
        CALL xml_NewElement(xp, "mdiis1d_step")
           CALL xml_addCharacters(xp, obj%mdiis1d_step, fmt='s16')
        CALL xml_EndElement(xp, "mdiis1d_step")
     END IF
     IF (obj%mdiis3d_step_ispresent) THEN
        CALL xml_NewElement(xp, "mdiis3d_step")
           CALL xml_addCharacters(xp, obj%mdiis3d_step, fmt='s16')
        CALL xml_EndElement(xp, "mdiis3d_step")
     END IF
     IF (obj%rism1d_bond_width_ispresent) THEN
        CALL xml_NewElement(xp, "rism1d_bond_width")
           CALL xml_addCharacters(xp, obj%rism1d_bond_width, fmt='s16')
        CALL xml_EndElement(xp, "rism1d_bond_width")
     END IF
     IF (obj%rism1d_dielectric_ispresent) THEN
        CALL xml_NewElement(xp, "rism1d_dielectric")
           CALL xml_addCharacters(xp, obj%rism1d_dielectric, fmt='s16')
        CALL xml_EndElement(xp, "rism1d_dielectric")
     END IF
     IF (obj%rism1d_molesize_ispresent) THEN
        CALL xml_NewElement(xp, "rism1d_molesize")
           CALL xml_addCharacters(xp, obj%rism1d_molesize, fmt='s16')
        CALL xml_EndElement(xp, "rism1d_molesize")
     END IF
     IF (obj%rism1d_nproc_ispresent) THEN
        CALL xml_NewElement(xp, "rism1d_nproc")
           CALL xml_addCharacters(xp, obj%rism1d_nproc)
        CALL xml_EndElement(xp, "rism1d_nproc")
     END IF
     IF (obj%rism1d_nproc_switch_ispresent) THEN
        CALL xml_NewElement(xp, "rism1d_nproc_switch")
           CALL xml_addCharacters(xp, obj%rism1d_nproc_switch)
        CALL xml_EndElement(xp, "rism1d_nproc_switch")
     END IF
     IF (obj%rism3d_conv_level_ispresent) THEN
        CALL xml_NewElement(xp, "rism3d_conv_level")
           CALL xml_addCharacters(xp, obj%rism3d_conv_level, fmt='s16')
        CALL xml_EndElement(xp, "rism3d_conv_level")
     END IF
     IF (obj%rism3d_planar_average_ispresent) THEN
        CALL xml_NewElement(xp, "rism3d_planar_average")
           CALL xml_addCharacters(xp, obj%rism3d_planar_average)
        CALL xml_EndElement(xp, "rism3d_planar_average")
     END IF
     IF (obj%laue_nfit_ispresent) THEN
        CALL xml_NewElement(xp, "laue_nfit")
           CALL xml_addCharacters(xp, obj%laue_nfit)
        CALL xml_EndElement(xp, "laue_nfit")
     END IF
     IF (obj%laue_expand_right_ispresent) THEN
        CALL xml_NewElement(xp, "laue_expand_right")
           CALL xml_addCharacters(xp, obj%laue_expand_right, fmt='s16')
        CALL xml_EndElement(xp, "laue_expand_right")
     END IF
     IF (obj%laue_expand_left_ispresent) THEN
        CALL xml_NewElement(xp, "laue_expand_left")
           CALL xml_addCharacters(xp, obj%laue_expand_left, fmt='s16')
        CALL xml_EndElement(xp, "laue_expand_left")
     END IF
     IF (obj%laue_starting_right_ispresent) THEN
        CALL xml_NewElement(xp, "laue_starting_right")
           CALL xml_addCharacters(xp, obj%laue_starting_right, fmt='s16')
        CALL xml_EndElement(xp, "laue_starting_right")
     END IF
     IF (obj%laue_starting_left_ispresent) THEN
        CALL xml_NewElement(xp, "laue_starting_left")
           CALL xml_addCharacters(xp, obj%laue_starting_left, fmt='s16')
        CALL xml_EndElement(xp, "laue_starting_left")
     END IF
     IF (obj%laue_buffer_right_ispresent) THEN
        CALL xml_NewElement(xp, "laue_buffer_right")
           CALL xml_addCharacters(xp, obj%laue_buffer_right, fmt='s16')
        CALL xml_EndElement(xp, "laue_buffer_right")
     END IF
     IF (obj%laue_buffer_right_solu_ispresent) THEN
        CALL xml_NewElement(xp, "laue_buffer_right_solu")
           CALL xml_addCharacters(xp, obj%laue_buffer_right_solu, fmt='s16')
        CALL xml_EndElement(xp, "laue_buffer_right_solu")
     END IF
     IF (obj%laue_buffer_right_solv_ispresent) THEN
        CALL xml_NewElement(xp, "laue_buffer_right_solv")
           CALL xml_addCharacters(xp, obj%laue_buffer_right_solv, fmt='s16')
        CALL xml_EndElement(xp, "laue_buffer_right_solv")
     END IF
     IF (obj%laue_buffer_left_ispresent) THEN
        CALL xml_NewElement(xp, "laue_buffer_left")
           CALL xml_addCharacters(xp, obj%laue_buffer_left, fmt='s16')
        CALL xml_EndElement(xp, "laue_buffer_left")
     END IF
     IF (obj%laue_buffer_left_solu_ispresent) THEN
        CALL xml_NewElement(xp, "laue_buffer_left_solu")
           CALL xml_addCharacters(xp, obj%laue_buffer_left_solu, fmt='s16')
        CALL xml_EndElement(xp, "laue_buffer_left_solu")
     END IF
     IF (obj%laue_buffer_left_solv_ispresent) THEN
        CALL xml_NewElement(xp, "laue_buffer_left_solv")
           CALL xml_addCharacters(xp, obj%laue_buffer_left_solv, fmt='s16')
        CALL xml_EndElement(xp, "laue_buffer_left_solv")
     END IF
     IF (obj%laue_both_hands_ispresent) THEN
        CALL xml_NewElement(xp, "laue_both_hands")
           CALL xml_addCharacters(xp, obj%laue_both_hands)
        CALL xml_EndElement(xp, "laue_both_hands")
     END IF
     IF (obj%laue_reference_ispresent) THEN
        CALL xml_NewElement(xp, "laue_reference")
           CALL xml_addCharacters(xp, TRIM(obj%laue_reference))
        CALL xml_EndElement(xp, "laue_reference")
     END IF
     IF (obj%laue_wall_ispresent) THEN
        CALL xml_NewElement(xp, "laue_wall")
           CALL xml_addCharacters(xp, TRIM(obj%laue_wall))
        CALL xml_EndElement(xp, "laue_wall")
     END IF
     IF (obj%laue_wall_z_ispresent) THEN
        CALL xml_NewElement(xp, "laue_wall_z")
           CALL xml_addCharacters(xp, obj%laue_wall_z, fmt='s16')
        CALL xml_EndElement(xp, "laue_wall_z")
     END IF
     IF (obj%laue_wall_rho_ispresent) THEN
        CALL xml_NewElement(xp, "laue_wall_rho")
           CALL xml_addCharacters(xp, obj%laue_wall_rho, fmt='s16')
        CALL xml_EndElement(xp, "laue_wall_rho")
     END IF
     IF (obj%laue_wall_epsilon_ispresent) THEN
        CALL xml_NewElement(xp, "laue_wall_epsilon")
           CALL xml_addCharacters(xp, obj%laue_wall_epsilon, fmt='s16')
        CALL xml_EndElement(xp, "laue_wall_epsilon")
     END IF
     IF (obj%laue_wall_sigma_ispresent) THEN
        CALL xml_NewElement(xp, "laue_wall_sigma")
           CALL xml_addCharacters(xp, obj%laue_wall_sigma, fmt='s16')
        CALL xml_EndElement(xp, "laue_wall_sigma")
     END IF
     IF (obj%laue_wall_lj6_ispresent) THEN
        CALL xml_NewElement(xp, "laue_wall_lj6")
           CALL xml_addCharacters(xp, obj%laue_wall_lj6)
        CALL xml_EndElement(xp, "laue_wall_lj6")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_rism

   SUBROUTINE qes_write_solute(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(solute_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'solute_lj')
        CALL xml_addCharacters(xp, TRIM(obj%solute_lj))
     CALL xml_EndElement(xp, 'solute_lj')
     CALL xml_NewElement(xp, 'epsilon')
        CALL xml_addCharacters(xp, obj%epsilon, fmt='s16')
     CALL xml_EndElement(xp, 'epsilon')
     CALL xml_NewElement(xp, 'sigma')
        CALL xml_addCharacters(xp, obj%sigma, fmt='s16')
     CALL xml_EndElement(xp, 'sigma')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_solute

   SUBROUTINE qes_write_solvent(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(solvent_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'label')
        CALL xml_addCharacters(xp, TRIM(obj%label))
     CALL xml_EndElement(xp, 'label')
     CALL xml_NewElement(xp, 'molec_file')
        CALL xml_addCharacters(xp, TRIM(obj%molec_file))
     CALL xml_EndElement(xp, 'molec_file')
     CALL xml_NewElement(xp, 'density1')
        CALL xml_addCharacters(xp, obj%density1, fmt='s16')
     CALL xml_EndElement(xp, 'density1')
     IF (obj%density2_ispresent) THEN
        CALL xml_NewElement(xp, "density2")
           CALL xml_addCharacters(xp, obj%density2, fmt='s16')
        CALL xml_EndElement(xp, "density2")
     END IF
     IF (obj%unit_ispresent) THEN
        CALL xml_NewElement(xp, "unit")
           CALL xml_addCharacters(xp, TRIM(obj%unit))
        CALL xml_EndElement(xp, "unit")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_solvent

   SUBROUTINE qes_write_k_points_IBZ(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(k_points_IBZ_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%monkhorst_pack_ispresent) THEN
        CALL qes_write_monkhorst_pack (xp, obj%monkhorst_pack)
     END IF
     IF (obj%nk_ispresent) THEN
        CALL xml_NewElement(xp, "nk")
           CALL xml_addCharacters(xp, obj%nk)
        CALL xml_EndElement(xp, "nk")
     END IF
     IF (obj%k_point_ispresent) THEN
        DO i = 1, obj%ndim_k_point
           CALL qes_write_k_point(xp, obj%k_point(i) )
        END DO
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_k_points_IBZ

   SUBROUTINE qes_write_monkhorst_pack(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(monkhorst_pack_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%nk1_ispresent) CALL xml_addAttribute(xp, 'nk1', obj%nk1 )
     IF (obj%nk2_ispresent) CALL xml_addAttribute(xp, 'nk2', obj%nk2 )
     IF (obj%nk3_ispresent) CALL xml_addAttribute(xp, 'nk3', obj%nk3 )
     IF (obj%k1_ispresent) CALL xml_addAttribute(xp, 'k1', obj%k1 )
     IF (obj%k2_ispresent) CALL xml_addAttribute(xp, 'k2', obj%k2 )
     IF (obj%k3_ispresent) CALL xml_addAttribute(xp, 'k3', obj%k3 )
        CALL xml_AddCharacters(xp, TRIM(obj%monkhorst_pack))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_monkhorst_pack

   SUBROUTINE qes_write_k_point(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(k_point_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%weight_ispresent) CALL xml_addAttribute(xp, 'weight', obj%weight )
     IF (obj%label_ispresent) CALL xml_addAttribute(xp, 'label', TRIM(obj%label) )
        CALL xml_AddCharacters(xp, obj%k_point, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_k_point

   SUBROUTINE qes_write_ion_control(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(ion_control_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'ion_dynamics')
        CALL xml_addCharacters(xp, TRIM(obj%ion_dynamics))
     CALL xml_EndElement(xp, 'ion_dynamics')
     IF (obj%upscale_ispresent) THEN
        CALL xml_NewElement(xp, "upscale")
           CALL xml_addCharacters(xp, obj%upscale, fmt='s16')
        CALL xml_EndElement(xp, "upscale")
     END IF
     IF (obj%remove_rigid_rot_ispresent) THEN
        CALL xml_NewElement(xp, "remove_rigid_rot")
           CALL xml_addCharacters(xp, obj%remove_rigid_rot)
        CALL xml_EndElement(xp, "remove_rigid_rot")
     END IF
     IF (obj%refold_pos_ispresent) THEN
        CALL xml_NewElement(xp, "refold_pos")
           CALL xml_addCharacters(xp, obj%refold_pos)
        CALL xml_EndElement(xp, "refold_pos")
     END IF
     IF (obj%bfgs_ispresent) THEN
        CALL qes_write_bfgs (xp, obj%bfgs)
     END IF
     IF (obj%md_ispresent) THEN
        CALL qes_write_md (xp, obj%md)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_ion_control

   SUBROUTINE qes_write_bfgs(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(bfgs_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'ndim')
        CALL xml_addCharacters(xp, obj%ndim)
     CALL xml_EndElement(xp, 'ndim')
     CALL xml_NewElement(xp, 'trust_radius_min')
        CALL xml_addCharacters(xp, obj%trust_radius_min, fmt='s16')
     CALL xml_EndElement(xp, 'trust_radius_min')
     CALL xml_NewElement(xp, 'trust_radius_max')
        CALL xml_addCharacters(xp, obj%trust_radius_max, fmt='s16')
     CALL xml_EndElement(xp, 'trust_radius_max')
     CALL xml_NewElement(xp, 'trust_radius_init')
        CALL xml_addCharacters(xp, obj%trust_radius_init, fmt='s16')
     CALL xml_EndElement(xp, 'trust_radius_init')
     CALL xml_NewElement(xp, 'w1')
        CALL xml_addCharacters(xp, obj%w1, fmt='s16')
     CALL xml_EndElement(xp, 'w1')
     CALL xml_NewElement(xp, 'w2')
        CALL xml_addCharacters(xp, obj%w2, fmt='s16')
     CALL xml_EndElement(xp, 'w2')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_bfgs

   SUBROUTINE qes_write_md(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(md_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'pot_extrapolation')
        CALL xml_addCharacters(xp, TRIM(obj%pot_extrapolation))
     CALL xml_EndElement(xp, 'pot_extrapolation')
     CALL xml_NewElement(xp, 'wfc_extrapolation')
        CALL xml_addCharacters(xp, TRIM(obj%wfc_extrapolation))
     CALL xml_EndElement(xp, 'wfc_extrapolation')
     CALL xml_NewElement(xp, 'ion_temperature')
        CALL xml_addCharacters(xp, TRIM(obj%ion_temperature))
     CALL xml_EndElement(xp, 'ion_temperature')
     CALL xml_NewElement(xp, 'timestep')
        CALL xml_addCharacters(xp, obj%timestep, fmt='s16')
     CALL xml_EndElement(xp, 'timestep')
     CALL xml_NewElement(xp, 'tempw')
        CALL xml_addCharacters(xp, obj%tempw, fmt='s16')
     CALL xml_EndElement(xp, 'tempw')
     CALL xml_NewElement(xp, 'tolp')
        CALL xml_addCharacters(xp, obj%tolp, fmt='s16')
     CALL xml_EndElement(xp, 'tolp')
     CALL xml_NewElement(xp, 'deltaT')
        CALL xml_addCharacters(xp, obj%deltaT, fmt='s16')
     CALL xml_EndElement(xp, 'deltaT')
     CALL xml_NewElement(xp, 'nraise')
        CALL xml_addCharacters(xp, obj%nraise)
     CALL xml_EndElement(xp, 'nraise')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_md

   SUBROUTINE qes_write_cell_control(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(cell_control_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'cell_dynamics')
        CALL xml_addCharacters(xp, TRIM(obj%cell_dynamics))
     CALL xml_EndElement(xp, 'cell_dynamics')
     CALL xml_NewElement(xp, 'pressure')
        CALL xml_addCharacters(xp, obj%pressure, fmt='s16')
     CALL xml_EndElement(xp, 'pressure')
     IF (obj%wmass_ispresent) THEN
        CALL xml_NewElement(xp, "wmass")
           CALL xml_addCharacters(xp, obj%wmass, fmt='s16')
        CALL xml_EndElement(xp, "wmass")
     END IF
     IF (obj%cell_factor_ispresent) THEN
        CALL xml_NewElement(xp, "cell_factor")
           CALL xml_addCharacters(xp, obj%cell_factor, fmt='s16')
        CALL xml_EndElement(xp, "cell_factor")
     END IF
     IF (obj%cell_do_free_ispresent) THEN
        CALL xml_NewElement(xp, "cell_do_free")
           CALL xml_addCharacters(xp, TRIM(obj%cell_do_free))
        CALL xml_EndElement(xp, "cell_do_free")
     END IF
     IF (obj%fix_volume_ispresent) THEN
        CALL xml_NewElement(xp, "fix_volume")
           CALL xml_addCharacters(xp, obj%fix_volume)
        CALL xml_EndElement(xp, "fix_volume")
     END IF
     IF (obj%fix_area_ispresent) THEN
        CALL xml_NewElement(xp, "fix_area")
           CALL xml_addCharacters(xp, obj%fix_area)
        CALL xml_EndElement(xp, "fix_area")
     END IF
     IF (obj%isotropic_ispresent) THEN
        CALL xml_NewElement(xp, "isotropic")
           CALL xml_addCharacters(xp, obj%isotropic)
        CALL xml_EndElement(xp, "isotropic")
     END IF
     IF (obj%free_cell_ispresent) THEN
        CALL qes_write_integerMatrix (xp, obj%free_cell)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_cell_control

   SUBROUTINE qes_write_symmetry_flags(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(symmetry_flags_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'nosym')
        CALL xml_addCharacters(xp, obj%nosym)
     CALL xml_EndElement(xp, 'nosym')
     CALL xml_NewElement(xp, 'nosym_evc')
        CALL xml_addCharacters(xp, obj%nosym_evc)
     CALL xml_EndElement(xp, 'nosym_evc')
     CALL xml_NewElement(xp, 'noinv')
        CALL xml_addCharacters(xp, obj%noinv)
     CALL xml_EndElement(xp, 'noinv')
     CALL xml_NewElement(xp, 'no_t_rev')
        CALL xml_addCharacters(xp, obj%no_t_rev)
     CALL xml_EndElement(xp, 'no_t_rev')
     CALL xml_NewElement(xp, 'force_symmorphic')
        CALL xml_addCharacters(xp, obj%force_symmorphic)
     CALL xml_EndElement(xp, 'force_symmorphic')
     CALL xml_NewElement(xp, 'use_all_frac')
        CALL xml_addCharacters(xp, obj%use_all_frac)
     CALL xml_EndElement(xp, 'use_all_frac')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_symmetry_flags

   SUBROUTINE qes_write_boundary_conditions(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(boundary_conditions_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'assume_isolated')
        CALL xml_addCharacters(xp, TRIM(obj%assume_isolated))
     CALL xml_EndElement(xp, 'assume_isolated')
     IF (obj%esm_ispresent) THEN
        CALL qes_write_esm (xp, obj%esm)
     END IF
     IF (obj%gcscf_ispresent) THEN
        CALL qes_write_gcscf (xp, obj%gcscf)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_boundary_conditions

   SUBROUTINE qes_write_esm(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(esm_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'bc')
        CALL xml_addCharacters(xp, TRIM(obj%bc))
     CALL xml_EndElement(xp, 'bc')
     IF (obj%nfit_ispresent) THEN
        CALL xml_NewElement(xp, "nfit")
           CALL xml_addCharacters(xp, obj%nfit)
        CALL xml_EndElement(xp, "nfit")
     END IF
     IF (obj%w_ispresent) THEN
        CALL xml_NewElement(xp, "w")
           CALL xml_addCharacters(xp, obj%w, fmt='s16')
        CALL xml_EndElement(xp, "w")
     END IF
     IF (obj%efield_ispresent) THEN
        CALL xml_NewElement(xp, "efield")
           CALL xml_addCharacters(xp, obj%efield, fmt='s16')
        CALL xml_EndElement(xp, "efield")
     END IF
     IF (obj%a_ispresent) THEN
        CALL xml_NewElement(xp, "a")
           CALL xml_addCharacters(xp, obj%a, fmt='s16')
        CALL xml_EndElement(xp, "a")
     END IF
     IF (obj%zb_ispresent) THEN
        CALL xml_NewElement(xp, "zb")
           CALL xml_addCharacters(xp, obj%zb, fmt='s16')
        CALL xml_EndElement(xp, "zb")
     END IF
     IF (obj%debug_ispresent) THEN
        CALL xml_NewElement(xp, "debug")
           CALL xml_addCharacters(xp, obj%debug)
        CALL xml_EndElement(xp, "debug")
     END IF
     IF (obj%debug_gpmax_ispresent) THEN
        CALL xml_NewElement(xp, "debug_gpmax")
           CALL xml_addCharacters(xp, obj%debug_gpmax)
        CALL xml_EndElement(xp, "debug_gpmax")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_esm

   SUBROUTINE qes_write_gcscf(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(gcscf_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%ignore_mun_ispresent) THEN
        CALL xml_NewElement(xp, "ignore_mun")
           CALL xml_addCharacters(xp, obj%ignore_mun)
        CALL xml_EndElement(xp, "ignore_mun")
     END IF
     IF (obj%mu_ispresent) THEN
        CALL xml_NewElement(xp, "mu")
           CALL xml_addCharacters(xp, obj%mu, fmt='s16')
        CALL xml_EndElement(xp, "mu")
     END IF
     IF (obj%conv_thr_ispresent) THEN
        CALL xml_NewElement(xp, "conv_thr")
           CALL xml_addCharacters(xp, obj%conv_thr, fmt='s16')
        CALL xml_EndElement(xp, "conv_thr")
     END IF
     IF (obj%gk_ispresent) THEN
        CALL xml_NewElement(xp, "gk")
           CALL xml_addCharacters(xp, obj%gk, fmt='s16')
        CALL xml_EndElement(xp, "gk")
     END IF
     IF (obj%gh_ispresent) THEN
        CALL xml_NewElement(xp, "gh")
           CALL xml_addCharacters(xp, obj%gh, fmt='s16')
        CALL xml_EndElement(xp, "gh")
     END IF
     IF (obj%beta_ispresent) THEN
        CALL xml_NewElement(xp, "beta")
           CALL xml_addCharacters(xp, obj%beta, fmt='s16')
        CALL xml_EndElement(xp, "beta")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_gcscf

   SUBROUTINE qes_write_solvents(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(solvents_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     DO i = 1, obj%ndim_solvent
        CALL qes_write_solvent(xp, obj%solvent(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_solvents

   SUBROUTINE qes_write_ekin_functional(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(ekin_functional_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'ecfixed')
        CALL xml_addCharacters(xp, obj%ecfixed, fmt='s16')
     CALL xml_EndElement(xp, 'ecfixed')
     CALL xml_NewElement(xp, 'qcutz')
        CALL xml_addCharacters(xp, obj%qcutz, fmt='s16')
     CALL xml_EndElement(xp, 'qcutz')
     CALL xml_NewElement(xp, 'q2sigma')
        CALL xml_addCharacters(xp, obj%q2sigma, fmt='s16')
     CALL xml_EndElement(xp, 'q2sigma')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_ekin_functional

   SUBROUTINE qes_write_spin_constraints(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(spin_constraints_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'spin_constraints')
        CALL xml_addCharacters(xp, TRIM(obj%spin_constraints))
     CALL xml_EndElement(xp, 'spin_constraints')
     CALL xml_NewElement(xp, 'lagrange_multiplier')
        CALL xml_addCharacters(xp, obj%lagrange_multiplier, fmt='s16')
     CALL xml_EndElement(xp, 'lagrange_multiplier')
     IF (obj%target_magnetization_ispresent) THEN
        CALL xml_NewElement(xp, "target_magnetization")
           CALL xml_addCharacters(xp, obj%target_magnetization, fmt='s16')
        CALL xml_EndElement(xp, "target_magnetization")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_spin_constraints

   SUBROUTINE qes_write_electric_field(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(electric_field_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'electric_potential')
        CALL xml_addCharacters(xp, TRIM(obj%electric_potential))
     CALL xml_EndElement(xp, 'electric_potential')
     IF (obj%dipole_correction_ispresent) THEN
        CALL xml_NewElement(xp, "dipole_correction")
           CALL xml_addCharacters(xp, obj%dipole_correction)
        CALL xml_EndElement(xp, "dipole_correction")
     END IF
     IF (obj%gate_settings_ispresent) THEN
        CALL qes_write_gate_settings (xp, obj%gate_settings)
     END IF
     IF (obj%electric_field_direction_ispresent) THEN
        CALL xml_NewElement(xp, "electric_field_direction")
           CALL xml_addCharacters(xp, obj%electric_field_direction)
        CALL xml_EndElement(xp, "electric_field_direction")
     END IF
     IF (obj%potential_max_position_ispresent) THEN
        CALL xml_NewElement(xp, "potential_max_position")
           CALL xml_addCharacters(xp, obj%potential_max_position, fmt='s16')
        CALL xml_EndElement(xp, "potential_max_position")
     END IF
     IF (obj%potential_decrease_width_ispresent) THEN
        CALL xml_NewElement(xp, "potential_decrease_width")
           CALL xml_addCharacters(xp, obj%potential_decrease_width, fmt='s16')
        CALL xml_EndElement(xp, "potential_decrease_width")
     END IF
     IF (obj%electric_field_amplitude_ispresent) THEN
        CALL xml_NewElement(xp, "electric_field_amplitude")
           CALL xml_addCharacters(xp, obj%electric_field_amplitude, fmt='s16')
        CALL xml_EndElement(xp, "electric_field_amplitude")
     END IF
     IF (obj%electric_field_vector_ispresent) THEN
        CALL xml_NewElement(xp, "electric_field_vector")
           CALL xml_addCharacters(xp, obj%electric_field_vector, fmt='s16')
        CALL xml_EndElement(xp, "electric_field_vector")
     END IF
     IF (obj%nk_per_string_ispresent) THEN
        CALL xml_NewElement(xp, "nk_per_string")
           CALL xml_addCharacters(xp, obj%nk_per_string)
        CALL xml_EndElement(xp, "nk_per_string")
     END IF
     IF (obj%n_berry_cycles_ispresent) THEN
        CALL xml_NewElement(xp, "n_berry_cycles")
           CALL xml_addCharacters(xp, obj%n_berry_cycles)
        CALL xml_EndElement(xp, "n_berry_cycles")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_electric_field

   SUBROUTINE qes_write_gate_settings(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(gate_settings_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'use_gate')
        CALL xml_addCharacters(xp, obj%use_gate)
     CALL xml_EndElement(xp, 'use_gate')
     IF (obj%zgate_ispresent) THEN
        CALL xml_NewElement(xp, "zgate")
           CALL xml_addCharacters(xp, obj%zgate, fmt='s16')
        CALL xml_EndElement(xp, "zgate")
     END IF
     IF (obj%relaxz_ispresent) THEN
        CALL xml_NewElement(xp, "relaxz")
           CALL xml_addCharacters(xp, obj%relaxz)
        CALL xml_EndElement(xp, "relaxz")
     END IF
     IF (obj%block_ispresent) THEN
        CALL xml_NewElement(xp, "block")
           CALL xml_addCharacters(xp, obj%block)
        CALL xml_EndElement(xp, "block")
     END IF
     IF (obj%block_1_ispresent) THEN
        CALL xml_NewElement(xp, "block_1")
           CALL xml_addCharacters(xp, obj%block_1, fmt='s16')
        CALL xml_EndElement(xp, "block_1")
     END IF
     IF (obj%block_2_ispresent) THEN
        CALL xml_NewElement(xp, "block_2")
           CALL xml_addCharacters(xp, obj%block_2, fmt='s16')
        CALL xml_EndElement(xp, "block_2")
     END IF
     IF (obj%block_height_ispresent) THEN
        CALL xml_NewElement(xp, "block_height")
           CALL xml_addCharacters(xp, obj%block_height, fmt='s16')
        CALL xml_EndElement(xp, "block_height")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_gate_settings

   SUBROUTINE qes_write_atomic_constraints(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(atomic_constraints_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'num_of_constraints')
        CALL xml_addCharacters(xp, obj%num_of_constraints)
     CALL xml_EndElement(xp, 'num_of_constraints')
     CALL xml_NewElement(xp, 'tolerance')
        CALL xml_addCharacters(xp, obj%tolerance, fmt='s16')
     CALL xml_EndElement(xp, 'tolerance')
     DO i = 1, obj%ndim_atomic_constraint
        CALL qes_write_atomic_constraint(xp, obj%atomic_constraint(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_atomic_constraints

   SUBROUTINE qes_write_atomic_constraint(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(atomic_constraint_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'constr_parms')
        CALL xml_addCharacters(xp, obj%constr_parms, fmt='s16')
     CALL xml_EndElement(xp, 'constr_parms')
     CALL xml_NewElement(xp, 'constr_type')
        CALL xml_addCharacters(xp, TRIM(obj%constr_type))
     CALL xml_EndElement(xp, 'constr_type')
     IF (obj%constr_target_ispresent) THEN
        CALL xml_NewElement(xp, "constr_target")
           CALL xml_addCharacters(xp, obj%constr_target, fmt='s16')
        CALL xml_EndElement(xp, "constr_target")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_atomic_constraint

   SUBROUTINE qes_write_inputOccupations(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(inputOccupations_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'size', obj%size )
     IF (obj%ispin_ispresent) CALL xml_addAttribute(xp, 'ispin', obj%ispin )
     IF (obj%spin_factor_ispresent) CALL xml_addAttribute(xp, 'spin_factor', obj%spin_factor )
     CALL xml_addNewLine(xp)
     DO i = 1, obj%size, 5
        CALL xml_AddCharacters(xp, obj%inputOccupations(i:MIN(i+5-1,obj%size)), fmt='s16')
        CALL xml_AddNewLine(xp)
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_inputOccupations

   SUBROUTINE qes_write_outputElectricField(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(outputElectricField_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%BerryPhase_ispresent) THEN
        CALL qes_write_BerryPhaseOutput (xp, obj%BerryPhase)
     END IF
     IF (obj%finiteElectricFieldInfo_ispresent) THEN
        CALL qes_write_finiteFieldOut (xp, obj%finiteElectricFieldInfo)
     END IF
     IF (obj%sawtoothEnergy_ispresent) THEN
        CALL qes_write_sawtoothEnergy (xp, obj%sawtoothEnergy)
     END IF
     IF (obj%dipoleInfo_ispresent) THEN
        CALL qes_write_dipoleOutput (xp, obj%dipoleInfo)
     END IF
     IF (obj%gateInfo_ispresent) THEN
        CALL qes_write_gateInfo (xp, obj%gateInfo)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_outputElectricField

   SUBROUTINE qes_write_BerryPhaseOutput(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(BerryPhaseOutput_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_polarization (xp, obj%totalPolarization)
     CALL qes_write_phase (xp, obj%totalPhase)
     DO i = 1, obj%ndim_ionicPolarization
        CALL qes_write_ionicPolarization(xp, obj%ionicPolarization(i) )
     END DO
     DO i = 1, obj%ndim_electronicPolarization
        CALL qes_write_electronicPolarization(xp, obj%electronicPolarization(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_BerryPhaseOutput

   SUBROUTINE qes_write_sawtoothEnergy(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(sawtoothEnergy_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%eamp_ispresent) CALL xml_addAttribute(xp, 'eamp', obj%eamp )
     IF (obj%eopreg_ispresent) CALL xml_addAttribute(xp, 'eopreg', obj%eopreg )
     IF (obj%emaxpos_ispresent) CALL xml_addAttribute(xp, 'emaxpos', obj%emaxpos )
     IF (obj%edir_ispresent) CALL xml_addAttribute(xp, 'edir', obj%edir )
        CALL xml_AddCharacters(xp, obj%sawtoothEnergy, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_sawtoothEnergy

   SUBROUTINE qes_write_dipoleOutput(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(dipoleOutput_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'idir')
        CALL xml_addCharacters(xp, obj%idir)
     CALL xml_EndElement(xp, 'idir')
     CALL qes_write_scalarQuantity (xp, obj%dipole)
     CALL qes_write_scalarQuantity (xp, obj%ion_dipole)
     CALL qes_write_scalarQuantity (xp, obj%elec_dipole)
     CALL qes_write_scalarQuantity (xp, obj%dipoleField)
     CALL qes_write_scalarQuantity (xp, obj%potentialAmp)
     CALL qes_write_scalarQuantity (xp, obj%totalLength)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_dipoleOutput

   SUBROUTINE qes_write_finiteFieldOut(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(finiteFieldOut_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'electronicDipole')
        CALL xml_addCharacters(xp, obj%electronicDipole, fmt='s16')
     CALL xml_EndElement(xp, 'electronicDipole')
     CALL xml_NewElement(xp, 'ionicDipole')
        CALL xml_addCharacters(xp, obj%ionicDipole, fmt='s16')
     CALL xml_EndElement(xp, 'ionicDipole')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_finiteFieldOut

   SUBROUTINE qes_write_polarization(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(polarization_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_scalarQuantity (xp, obj%polarization)
     CALL xml_NewElement(xp, 'modulus')
        CALL xml_addCharacters(xp, obj%modulus, fmt='s16')
     CALL xml_EndElement(xp, 'modulus')
     CALL xml_NewElement(xp, 'direction')
        CALL xml_addCharacters(xp, obj%direction, fmt='s16')
     CALL xml_EndElement(xp, 'direction')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_polarization

   SUBROUTINE qes_write_ionicPolarization(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(ionicPolarization_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_atom (xp, obj%ion)
     CALL xml_NewElement(xp, 'charge')
        CALL xml_addCharacters(xp, obj%charge, fmt='s16')
     CALL xml_EndElement(xp, 'charge')
     CALL qes_write_phase (xp, obj%phase)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_ionicPolarization

   SUBROUTINE qes_write_electronicPolarization(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(electronicPolarization_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_k_point (xp, obj%firstKeyPoint)
     IF (obj%spin_ispresent) THEN
        CALL xml_NewElement(xp, "spin")
           CALL xml_addCharacters(xp, obj%spin)
        CALL xml_EndElement(xp, "spin")
     END IF
     CALL qes_write_phase (xp, obj%phase)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_electronicPolarization

   SUBROUTINE qes_write_phase(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(phase_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%ionic_ispresent) CALL xml_addAttribute(xp, 'ionic', obj%ionic )
     IF (obj%electronic_ispresent) CALL xml_addAttribute(xp, 'electronic', obj%electronic )
     IF (obj%modulus_ispresent) CALL xml_addAttribute(xp, 'modulus', TRIM(obj%modulus) )
        CALL xml_AddCharacters(xp, obj%phase, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_phase

   SUBROUTINE qes_write_gateInfo(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(gateInfo_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'pot_prefactor')
        CALL xml_addCharacters(xp, obj%pot_prefactor, fmt='s16')
     CALL xml_EndElement(xp, 'pot_prefactor')
     CALL xml_NewElement(xp, 'gate_zpos')
        CALL xml_addCharacters(xp, obj%gate_zpos, fmt='s16')
     CALL xml_EndElement(xp, 'gate_zpos')
     CALL xml_NewElement(xp, 'gate_gate_term')
        CALL xml_addCharacters(xp, obj%gate_gate_term, fmt='s16')
     CALL xml_EndElement(xp, 'gate_gate_term')
     CALL xml_NewElement(xp, 'gatefieldEnergy')
        CALL xml_addCharacters(xp, obj%gatefieldEnergy, fmt='s16')
     CALL xml_EndElement(xp, 'gatefieldEnergy')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_gateInfo

   SUBROUTINE qes_write_convergence_info(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(convergence_info_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_scf_conv (xp, obj%scf_conv)
     IF (obj%opt_conv_ispresent) THEN
        CALL qes_write_opt_conv (xp, obj%opt_conv)
     END IF
     IF (obj%wf_collected_ispresent) THEN
        CALL xml_NewElement(xp, "wf_collected")
           CALL xml_addCharacters(xp, obj%wf_collected)
        CALL xml_EndElement(xp, "wf_collected")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_convergence_info

   SUBROUTINE qes_write_scf_conv(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(scf_conv_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'convergence_achieved')
        CALL xml_addCharacters(xp, obj%convergence_achieved)
     CALL xml_EndElement(xp, 'convergence_achieved')
     CALL xml_NewElement(xp, 'n_scf_steps')
        CALL xml_addCharacters(xp, obj%n_scf_steps)
     CALL xml_EndElement(xp, 'n_scf_steps')
     CALL xml_NewElement(xp, 'scf_error')
        CALL xml_addCharacters(xp, obj%scf_error, fmt='s16')
     CALL xml_EndElement(xp, 'scf_error')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_scf_conv

   SUBROUTINE qes_write_opt_conv(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(opt_conv_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'convergence_achieved')
        CALL xml_addCharacters(xp, obj%convergence_achieved)
     CALL xml_EndElement(xp, 'convergence_achieved')
     CALL xml_NewElement(xp, 'n_opt_steps')
        CALL xml_addCharacters(xp, obj%n_opt_steps)
     CALL xml_EndElement(xp, 'n_opt_steps')
     CALL xml_NewElement(xp, 'grad_norm')
        CALL xml_addCharacters(xp, obj%grad_norm, fmt='s16')
     CALL xml_EndElement(xp, 'grad_norm')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_opt_conv

   SUBROUTINE qes_write_algorithmic_info(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(algorithmic_info_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'real_space_q')
        CALL xml_addCharacters(xp, obj%real_space_q)
     CALL xml_EndElement(xp, 'real_space_q')
     IF (obj%real_space_beta_ispresent) THEN
        CALL xml_NewElement(xp, "real_space_beta")
           CALL xml_addCharacters(xp, obj%real_space_beta)
        CALL xml_EndElement(xp, "real_space_beta")
     END IF
     CALL xml_NewElement(xp, 'uspp')
        CALL xml_addCharacters(xp, obj%uspp)
     CALL xml_EndElement(xp, 'uspp')
     CALL xml_NewElement(xp, 'paw')
        CALL xml_addCharacters(xp, obj%paw)
     CALL xml_EndElement(xp, 'paw')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_algorithmic_info

   SUBROUTINE qes_write_symmetries(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(symmetries_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'nsym')
        CALL xml_addCharacters(xp, obj%nsym)
     CALL xml_EndElement(xp, 'nsym')
     IF (obj%colin_mag_ispresent) THEN
        CALL xml_NewElement(xp, "colin_mag")
           CALL xml_addCharacters(xp, obj%colin_mag)
        CALL xml_EndElement(xp, "colin_mag")
     END IF
     CALL xml_NewElement(xp, 'nrot')
        CALL xml_addCharacters(xp, obj%nrot)
     CALL xml_EndElement(xp, 'nrot')
     CALL xml_NewElement(xp, 'space_group')
        CALL xml_addCharacters(xp, obj%space_group)
     CALL xml_EndElement(xp, 'space_group')
     DO i = 1, obj%ndim_symmetry
        CALL qes_write_symmetry(xp, obj%symmetry(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_symmetries

   SUBROUTINE qes_write_symmetry(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(symmetry_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_info (xp, obj%info)
     CALL qes_write_matrix (xp, obj%rotation)
     IF (obj%fractional_translation_ispresent) THEN
        CALL xml_NewElement(xp, "fractional_translation")
           CALL xml_addCharacters(xp, obj%fractional_translation, fmt='s16')
        CALL xml_EndElement(xp, "fractional_translation")
     END IF
     IF (obj%equivalent_atoms_ispresent) THEN
        CALL qes_write_equivalent_atoms (xp, obj%equivalent_atoms)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_symmetry

   SUBROUTINE qes_write_equivalent_atoms(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(equivalent_atoms_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'size', obj%size )
     IF (obj%nat_ispresent) CALL xml_addAttribute(xp, 'nat', obj%nat )
     CALL xml_addNewLine(xp)
     DO i = 1, obj%size, 8
        CALL xml_AddCharacters(xp, obj%equivalent_atoms(i:MIN(i+8-1, obj%size)))
        CALL xml_AddNewLine(xp)
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_equivalent_atoms

   SUBROUTINE qes_write_info(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(info_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%name_ispresent) CALL xml_addAttribute(xp, 'name', TRIM(obj%name) )
     IF (obj%class_ispresent) CALL xml_addAttribute(xp, 'class', TRIM(obj%class) )
     IF (obj%time_reversal_ispresent) CALL xml_addAttribute(xp, 'time_reversal', obj%time_reversal )
        CALL xml_AddCharacters(xp, TRIM(obj%info))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_info

   SUBROUTINE qes_write_outputPBC(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(outputPBC_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'assume_isolated')
        CALL xml_addCharacters(xp, TRIM(obj%assume_isolated))
     CALL xml_EndElement(xp, 'assume_isolated')
     IF (obj%esm_ispresent) THEN
        CALL qes_write_esm (xp, obj%esm)
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_outputPBC

   SUBROUTINE qes_write_magnetization(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(magnetization_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'lsda')
        CALL xml_addCharacters(xp, obj%lsda)
     CALL xml_EndElement(xp, 'lsda')
     CALL xml_NewElement(xp, 'noncolin')
        CALL xml_addCharacters(xp, obj%noncolin)
     CALL xml_EndElement(xp, 'noncolin')
     CALL xml_NewElement(xp, 'spinorbit')
        CALL xml_addCharacters(xp, obj%spinorbit)
     CALL xml_EndElement(xp, 'spinorbit')
     IF (obj%total_ispresent) THEN
        CALL xml_NewElement(xp, "total")
           CALL xml_addCharacters(xp, obj%total, fmt='s16')
        CALL xml_EndElement(xp, "total")
     END IF
     IF (obj%total_vec_ispresent) THEN
        CALL xml_NewElement(xp, "total_vec")
           CALL xml_addCharacters(xp, obj%total_vec, fmt='s16')
        CALL xml_EndElement(xp, "total_vec")
     END IF
     CALL xml_NewElement(xp, 'absolute')
        CALL xml_addCharacters(xp, obj%absolute, fmt='s16')
     CALL xml_EndElement(xp, 'absolute')
     IF (obj%Scalar_Site_Magnetic_Moments_ispresent) THEN
        CALL qes_write_scalmags (xp, obj%Scalar_Site_Magnetic_Moments)
     END IF
     IF (obj%Site_Magnetizations_ispresent) THEN
        CALL qes_write_d3mags (xp, obj%Site_Magnetizations)
     END IF
     IF (obj%do_magnetization_ispresent) THEN
        CALL xml_NewElement(xp, "do_magnetization")
           CALL xml_addCharacters(xp, obj%do_magnetization)
        CALL xml_EndElement(xp, "do_magnetization")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_magnetization

   SUBROUTINE qes_write_total_energy(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(total_energy_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'etot')
        CALL xml_addCharacters(xp, obj%etot, fmt='s16')
     CALL xml_EndElement(xp, 'etot')
     IF (obj%eband_ispresent) THEN
        CALL xml_NewElement(xp, "eband")
           CALL xml_addCharacters(xp, obj%eband, fmt='s16')
        CALL xml_EndElement(xp, "eband")
     END IF
     IF (obj%ehart_ispresent) THEN
        CALL xml_NewElement(xp, "ehart")
           CALL xml_addCharacters(xp, obj%ehart, fmt='s16')
        CALL xml_EndElement(xp, "ehart")
     END IF
     IF (obj%vtxc_ispresent) THEN
        CALL xml_NewElement(xp, "vtxc")
           CALL xml_addCharacters(xp, obj%vtxc, fmt='s16')
        CALL xml_EndElement(xp, "vtxc")
     END IF
     IF (obj%etxc_ispresent) THEN
        CALL xml_NewElement(xp, "etxc")
           CALL xml_addCharacters(xp, obj%etxc, fmt='s16')
        CALL xml_EndElement(xp, "etxc")
     END IF
     IF (obj%ewald_ispresent) THEN
        CALL xml_NewElement(xp, "ewald")
           CALL xml_addCharacters(xp, obj%ewald, fmt='s16')
        CALL xml_EndElement(xp, "ewald")
     END IF
     IF (obj%demet_ispresent) THEN
        CALL xml_NewElement(xp, "demet")
           CALL xml_addCharacters(xp, obj%demet, fmt='s16')
        CALL xml_EndElement(xp, "demet")
     END IF
     IF (obj%efieldcorr_ispresent) THEN
        CALL xml_NewElement(xp, "efieldcorr")
           CALL xml_addCharacters(xp, obj%efieldcorr, fmt='s16')
        CALL xml_EndElement(xp, "efieldcorr")
     END IF
     IF (obj%potentiostat_contr_ispresent) THEN
        CALL xml_NewElement(xp, "potentiostat_contr")
           CALL xml_addCharacters(xp, obj%potentiostat_contr, fmt='s16')
        CALL xml_EndElement(xp, "potentiostat_contr")
     END IF
     IF (obj%gatefield_contr_ispresent) THEN
        CALL xml_NewElement(xp, "gatefield_contr")
           CALL xml_addCharacters(xp, obj%gatefield_contr, fmt='s16')
        CALL xml_EndElement(xp, "gatefield_contr")
     END IF
     IF (obj%vdW_term_ispresent) THEN
        CALL xml_NewElement(xp, "vdW_term")
           CALL xml_addCharacters(xp, obj%vdW_term, fmt='s16')
        CALL xml_EndElement(xp, "vdW_term")
     END IF
     IF (obj%esol_ispresent) THEN
        CALL xml_NewElement(xp, "esol")
           CALL xml_addCharacters(xp, obj%esol, fmt='s16')
        CALL xml_EndElement(xp, "esol")
     END IF
     IF (obj%levelshift_contr_ispresent) THEN
        CALL xml_NewElement(xp, "levelshift_contr")
           CALL xml_addCharacters(xp, obj%levelshift_contr, fmt='s16')
        CALL xml_EndElement(xp, "levelshift_contr")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_total_energy

   SUBROUTINE qes_write_band_structure(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(band_structure_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'lsda')
        CALL xml_addCharacters(xp, obj%lsda)
     CALL xml_EndElement(xp, 'lsda')
     CALL xml_NewElement(xp, 'noncolin')
        CALL xml_addCharacters(xp, obj%noncolin)
     CALL xml_EndElement(xp, 'noncolin')
     CALL xml_NewElement(xp, 'spinorbit')
        CALL xml_addCharacters(xp, obj%spinorbit)
     CALL xml_EndElement(xp, 'spinorbit')
     IF (obj%nbnd_ispresent) THEN
        CALL xml_NewElement(xp, "nbnd")
           CALL xml_addCharacters(xp, obj%nbnd)
        CALL xml_EndElement(xp, "nbnd")
     END IF
     IF (obj%nbnd_up_ispresent) THEN
        CALL xml_NewElement(xp, "nbnd_up")
           CALL xml_addCharacters(xp, obj%nbnd_up)
        CALL xml_EndElement(xp, "nbnd_up")
     END IF
     IF (obj%nbnd_dw_ispresent) THEN
        CALL xml_NewElement(xp, "nbnd_dw")
           CALL xml_addCharacters(xp, obj%nbnd_dw)
        CALL xml_EndElement(xp, "nbnd_dw")
     END IF
     CALL xml_NewElement(xp, 'nelec')
        CALL xml_addCharacters(xp, obj%nelec, fmt='s16')
     CALL xml_EndElement(xp, 'nelec')
     IF (obj%fermi_energy_ispresent) THEN
        CALL xml_NewElement(xp, "fermi_energy")
           CALL xml_addCharacters(xp, obj%fermi_energy, fmt='s16')
        CALL xml_EndElement(xp, "fermi_energy")
     END IF
     IF (obj%highestOccupiedLevel_ispresent) THEN
        CALL xml_NewElement(xp, "highestOccupiedLevel")
           CALL xml_addCharacters(xp, obj%highestOccupiedLevel, fmt='s16')
        CALL xml_EndElement(xp, "highestOccupiedLevel")
     END IF
     IF (obj%lowestUnoccupiedLevel_ispresent) THEN
        CALL xml_NewElement(xp, "lowestUnoccupiedLevel")
           CALL xml_addCharacters(xp, obj%lowestUnoccupiedLevel, fmt='s16')
        CALL xml_EndElement(xp, "lowestUnoccupiedLevel")
     END IF
     IF (obj%two_fermi_energies_ispresent) THEN
        CALL xml_NewElement(xp, "two_fermi_energies")
           CALL xml_addCharacters(xp, obj%two_fermi_energies, fmt='s16')
        CALL xml_EndElement(xp, "two_fermi_energies")
     END IF
     CALL qes_write_k_points_IBZ (xp, obj%starting_k_points)
     CALL xml_NewElement(xp, 'nks')
        CALL xml_addCharacters(xp, obj%nks)
     CALL xml_EndElement(xp, 'nks')
     CALL qes_write_occupations (xp, obj%occupations_kind)
     IF (obj%smearing_ispresent) THEN
        CALL qes_write_smearing (xp, obj%smearing)
     END IF
     DO i = 1, obj%ndim_ks_energies
        CALL qes_write_ks_energies(xp, obj%ks_energies(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_band_structure

   SUBROUTINE qes_write_ks_energies(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(ks_energies_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_k_point (xp, obj%k_point)
     CALL xml_NewElement(xp, 'npw')
        CALL xml_addCharacters(xp, obj%npw)
     CALL xml_EndElement(xp, 'npw')
     CALL qes_write_vector (xp, obj%eigenvalues)
     CALL qes_write_vector (xp, obj%occupations)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_ks_energies

   SUBROUTINE qes_write_closed(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(closed_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%DATE_ispresent) CALL xml_addAttribute(xp, 'DATE', TRIM(obj%DATE) )
     IF (obj%TIME_ispresent) CALL xml_addAttribute(xp, 'TIME', TRIM(obj%TIME) )
        CALL xml_AddCharacters(xp, TRIM(obj%closed))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_closed

   SUBROUTINE qes_write_cpstatus(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(cpstatus_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL qes_write_cpnumstep (xp, obj%STEP)
     CALL qes_write_scalarQuantity (xp, obj%TIME)
     CALL xml_NewElement(xp, 'TITLE')
        CALL xml_addCharacters(xp, TRIM(obj%TITLE))
     CALL xml_EndElement(xp, 'TITLE')
     CALL qes_write_scalarQuantity (xp, obj%KINETIC_ENERGY)
     CALL qes_write_scalarQuantity (xp, obj%HARTREE_ENERGY)
     CALL qes_write_scalarQuantity (xp, obj%EWALD_TERM)
     CALL qes_write_scalarQuantity (xp, obj%GAUSS_SELFINT)
     CALL qes_write_scalarQuantity (xp, obj%LPSP_ENERGY)
     CALL qes_write_scalarQuantity (xp, obj%NLPSP_ENERGY)
     CALL qes_write_scalarQuantity (xp, obj%EXC_ENERGY)
     CALL qes_write_scalarQuantity (xp, obj%AVERAGE_POT)
     CALL qes_write_scalarQuantity (xp, obj%ENTHALPY)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_cpstatus

   SUBROUTINE qes_write_cpnumstep(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(cpnumstep_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%ITERATION_ispresent) CALL xml_addAttribute(xp, 'ITERATION', obj%ITERATION )
        CALL xml_AddCharacters(xp, TRIM(obj%cpnumstep))
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_cpnumstep

   SUBROUTINE qes_write_cptimesteps(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(cptimesteps_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%nt_ispresent) CALL xml_addAttribute(xp, 'nt', obj%nt )
     CALL qes_write_cpstep (xp, obj%STEP0)
     CALL qes_write_cpstep (xp, obj%STEPM)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_cptimesteps

   SUBROUTINE qes_write_cpstep(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(cpstep_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%ACCUMULATORS_ispresent) THEN
        CALL xml_NewElement(xp, "ACCUMULATORS")
           CALL xml_addCharacters(xp, obj%ACCUMULATORS, fmt='s16')
        CALL xml_EndElement(xp, "ACCUMULATORS")
     END IF
     CALL qes_write_cp_ionPos (xp, obj%IONS_POSITIONS)
     CALL qes_write_cp_ionsNose (xp, obj%IONS_NOSE)
     IF (obj%ekincm_ispresent) THEN
        CALL xml_NewElement(xp, "ekincm")
           CALL xml_addCharacters(xp, obj%ekincm, fmt='s16')
        CALL xml_EndElement(xp, "ekincm")
     END IF
     CALL qes_write_cp_elecNose (xp, obj%ELECTRONS_NOSE)
     CALL qes_write_cp_cell (xp, obj%CELL_PARAMETERS)
     CALL qes_write_cp_cellNose (xp, obj%CELL_NOSE)
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_cpstep

   SUBROUTINE qes_write_cp_ionPos(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(cp_ionPos_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'stau')
        CALL xml_addCharacters(xp, obj%stau, fmt='s16')
     CALL xml_EndElement(xp, 'stau')
     CALL xml_NewElement(xp, 'svel')
        CALL xml_addCharacters(xp, obj%svel, fmt='s16')
     CALL xml_EndElement(xp, 'svel')
     IF (obj%taui_ispresent) THEN
        CALL xml_NewElement(xp, "taui")
           CALL xml_addCharacters(xp, obj%taui, fmt='s16')
        CALL xml_EndElement(xp, "taui")
     END IF
     IF (obj%cdmi_ispresent) THEN
        CALL xml_NewElement(xp, "cdmi")
           CALL xml_addCharacters(xp, obj%cdmi, fmt='s16')
        CALL xml_EndElement(xp, "cdmi")
     END IF
     IF (obj%force_ispresent) THEN
        CALL xml_NewElement(xp, "force")
           CALL xml_addCharacters(xp, obj%force, fmt='s16')
        CALL xml_EndElement(xp, "force")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_cp_ionPos

   SUBROUTINE qes_write_cp_ionsNose(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(cp_ionsNose_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'nhpcl')
        CALL xml_addCharacters(xp, obj%nhpcl)
     CALL xml_EndElement(xp, 'nhpcl')
     CALL xml_NewElement(xp, 'nhpdim')
        CALL xml_addCharacters(xp, obj%nhpdim)
     CALL xml_EndElement(xp, 'nhpdim')
     CALL xml_NewElement(xp, 'xnhp')
        CALL xml_addCharacters(xp, obj%xnhp, fmt='s16')
     CALL xml_EndElement(xp, 'xnhp')
     IF (obj%vnhp_ispresent) THEN
        CALL xml_NewElement(xp, "vnhp")
           CALL xml_addCharacters(xp, obj%vnhp, fmt='s16')
        CALL xml_EndElement(xp, "vnhp")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_cp_ionsNose

   SUBROUTINE qes_write_cp_elecNose(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(cp_elecNose_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'xnhe')
        CALL xml_addCharacters(xp, obj%xnhe, fmt='s16')
     CALL xml_EndElement(xp, 'xnhe')
     IF (obj%vnhe_ispresent) THEN
        CALL xml_NewElement(xp, "vnhe")
           CALL xml_addCharacters(xp, obj%vnhe, fmt='s16')
        CALL xml_EndElement(xp, "vnhe")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_cp_elecNose

   SUBROUTINE qes_write_cp_cell(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(cp_cell_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'ht')
        CALL xml_addCharacters(xp, obj%ht, fmt='s16')
     CALL xml_EndElement(xp, 'ht')
     IF (obj%htvel_ispresent) THEN
        CALL xml_NewElement(xp, "htvel")
           CALL xml_addCharacters(xp, obj%htvel, fmt='s16')
        CALL xml_EndElement(xp, "htvel")
     END IF
     IF (obj%gvel_ispresent) THEN
        CALL xml_NewElement(xp, "gvel")
           CALL xml_addCharacters(xp, obj%gvel, fmt='s16')
        CALL xml_EndElement(xp, "gvel")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_cp_cell

   SUBROUTINE qes_write_cp_cellNose(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(cp_cellNose_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'xnhh')
        CALL xml_addCharacters(xp, obj%xnhh, fmt='s16')
     CALL xml_EndElement(xp, 'xnhh')
     IF (obj%vnhh_ispresent) THEN
        CALL xml_NewElement(xp, "vnhh")
           CALL xml_addCharacters(xp, obj%vnhh, fmt='s16')
        CALL xml_EndElement(xp, "vnhh")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_cp_cellNose

   SUBROUTINE qes_write_scalmags(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(scalmags_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%nat_ispresent) CALL xml_addAttribute(xp, 'nat', obj%nat )
     DO i = 1, obj%ndim_SiteMagnetization
        CALL qes_write_SiteMoment(xp, obj%SiteMagnetization(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_scalmags

   SUBROUTINE qes_write_d3mags(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(d3mags_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%nat_ispresent) CALL xml_addAttribute(xp, 'nat', obj%nat )
     DO i = 1, obj%ndim_SiteMagnetization
        CALL qes_write_SitMag(xp, obj%SiteMagnetization(i) )
     END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_d3mags

   SUBROUTINE qes_write_integerMatrix(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(integerMatrix_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_addAttribute(xp, 'rank', obj%rank )
     CALL xml_addAttribute(xp, 'dims', obj%dims )
     IF (obj%order_ispresent) CALL xml_addAttribute(xp, 'order', TRIM(obj%order) )
        CALL xml_addNewLine(xp)
        DO i = 1, obj%dims(2)
           CALL xml_AddCharacters(xp, obj%integerMatrix((i-1)*obj%dims(1)+1: i*obj%dims(1)) )
           CALL xml_addNewLine(xp)
        END DO
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_integerMatrix

   SUBROUTINE qes_write_scalarQuantity(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(scalarQuantity_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%Units_ispresent) CALL xml_addAttribute(xp, 'Units', TRIM(obj%Units) )
        CALL xml_AddCharacters(xp, obj%scalarQuantity, fmt='s16')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_scalarQuantity

   SUBROUTINE qes_write_rism3d(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(rism3d_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'nmol')
        CALL xml_addCharacters(xp, obj%nmol)
     CALL xml_EndElement(xp, 'nmol')
     IF (obj%molec_dir_ispresent) THEN
        CALL xml_NewElement(xp, "molec_dir")
           CALL xml_addCharacters(xp, TRIM(obj%molec_dir))
        CALL xml_EndElement(xp, "molec_dir")
     END IF
     DO i = 1, obj%ndim_solvent
        CALL qes_write_solvent(xp, obj%solvent(i) )
     END DO
     CALL xml_NewElement(xp, 'ecutsolv')
        CALL xml_addCharacters(xp, obj%ecutsolv, fmt='s16')
     CALL xml_EndElement(xp, 'ecutsolv')
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_rism3d

   SUBROUTINE qes_write_rismlaue(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(rismlaue_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     IF (obj%both_hands_ispresent) THEN
        CALL xml_NewElement(xp, "both_hands")
           CALL xml_addCharacters(xp, obj%both_hands)
        CALL xml_EndElement(xp, "both_hands")
     END IF
     IF (obj%nfit_ispresent) THEN
        CALL xml_NewElement(xp, "nfit")
           CALL xml_addCharacters(xp, obj%nfit)
        CALL xml_EndElement(xp, "nfit")
     END IF
     IF (obj%pot_ref_ispresent) THEN
        CALL xml_NewElement(xp, "pot_ref")
           CALL xml_addCharacters(xp, obj%pot_ref)
        CALL xml_EndElement(xp, "pot_ref")
     END IF
     IF (obj%charge_ispresent) THEN
        CALL xml_NewElement(xp, "charge")
           CALL xml_addCharacters(xp, obj%charge, fmt='s16')
        CALL xml_EndElement(xp, "charge")
     END IF
     IF (obj%right_start_ispresent) THEN
        CALL xml_NewElement(xp, "right_start")
           CALL xml_addCharacters(xp, obj%right_start, fmt='s16')
        CALL xml_EndElement(xp, "right_start")
     END IF
     IF (obj%right_expand_ispresent) THEN
        CALL xml_NewElement(xp, "right_expand")
           CALL xml_addCharacters(xp, obj%right_expand, fmt='s16')
        CALL xml_EndElement(xp, "right_expand")
     END IF
     IF (obj%right_buffer_ispresent) THEN
        CALL xml_NewElement(xp, "right_buffer")
           CALL xml_addCharacters(xp, obj%right_buffer, fmt='s16')
        CALL xml_EndElement(xp, "right_buffer")
     END IF
     IF (obj%right_buffer_u_ispresent) THEN
        CALL xml_NewElement(xp, "right_buffer_u")
           CALL xml_addCharacters(xp, obj%right_buffer_u, fmt='s16')
        CALL xml_EndElement(xp, "right_buffer_u")
     END IF
     IF (obj%right_buffer_v_ispresent) THEN
        CALL xml_NewElement(xp, "right_buffer_v")
           CALL xml_addCharacters(xp, obj%right_buffer_v, fmt='s16')
        CALL xml_EndElement(xp, "right_buffer_v")
     END IF
     IF (obj%left_start_ispresent) THEN
        CALL xml_NewElement(xp, "left_start")
           CALL xml_addCharacters(xp, obj%left_start, fmt='s16')
        CALL xml_EndElement(xp, "left_start")
     END IF
     IF (obj%left_expand_ispresent) THEN
        CALL xml_NewElement(xp, "left_expand")
           CALL xml_addCharacters(xp, obj%left_expand, fmt='s16')
        CALL xml_EndElement(xp, "left_expand")
     END IF
     IF (obj%left_buffer_ispresent) THEN
        CALL xml_NewElement(xp, "left_buffer")
           CALL xml_addCharacters(xp, obj%left_buffer, fmt='s16')
        CALL xml_EndElement(xp, "left_buffer")
     END IF
     IF (obj%left_buffer_u_ispresent) THEN
        CALL xml_NewElement(xp, "left_buffer_u")
           CALL xml_addCharacters(xp, obj%left_buffer_u, fmt='s16')
        CALL xml_EndElement(xp, "left_buffer_u")
     END IF
     IF (obj%left_buffer_v_ispresent) THEN
        CALL xml_NewElement(xp, "left_buffer_v")
           CALL xml_addCharacters(xp, obj%left_buffer_v, fmt='s16')
        CALL xml_EndElement(xp, "left_buffer_v")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_rismlaue

   SUBROUTINE qes_write_two_chem(xp, obj)
     !-----------------------------------------------------------------
     IMPLICIT NONE
     TYPE (xmlf_t),INTENT(INOUT)                      :: xp
     TYPE(two_chem_type),INTENT(IN)    :: obj
     ! 
     INTEGER                                          :: i 
     ! 
     IF ( .NOT. obj%lwrite ) RETURN 
     ! 
     CALL xml_NewElement(xp, TRIM(obj%tagname))
     CALL xml_NewElement(xp, 'twochem')
        CALL xml_addCharacters(xp, obj%twochem)
     CALL xml_EndElement(xp, 'twochem')
     CALL xml_NewElement(xp, 'nbnd_cond')
        CALL xml_addCharacters(xp, obj%nbnd_cond)
     CALL xml_EndElement(xp, 'nbnd_cond')
     CALL xml_NewElement(xp, 'degauss_cond')
        CALL xml_addCharacters(xp, obj%degauss_cond, fmt='s16')
     CALL xml_EndElement(xp, 'degauss_cond')
     CALL xml_NewElement(xp, 'nelec_cond')
        CALL xml_addCharacters(xp, obj%nelec_cond, fmt='s16')
     CALL xml_EndElement(xp, 'nelec_cond')
     IF (obj%ef_cond_ispresent) THEN
        CALL xml_NewElement(xp, "ef_cond")
           CALL xml_addCharacters(xp, obj%ef_cond, fmt='s16')
        CALL xml_EndElement(xp, "ef_cond")
     END IF
     CALL xml_EndElement(xp, TRIM(obj%tagname))
   END SUBROUTINE qes_write_two_chem

  !
END MODULE qes_write_module