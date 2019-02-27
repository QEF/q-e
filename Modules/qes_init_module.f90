!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE qes_init_module
  !
  ! Auto-generated code: don't edit or at least don't commit changes
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes-1.0
  !
  USE kinds, only: DP
  USE qes_types_module
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: qes_init
  !
  INTERFACE qes_init
    !
    MODULE PROCEDURE qes_init_espresso
    MODULE PROCEDURE qes_init_general_info
    MODULE PROCEDURE qes_init_parallel_info
    MODULE PROCEDURE qes_init_input
    MODULE PROCEDURE qes_init_step
    MODULE PROCEDURE qes_init_output
    MODULE PROCEDURE qes_init_timing
    MODULE PROCEDURE qes_init_clock
    MODULE PROCEDURE qes_init_control_variables
    MODULE PROCEDURE qes_init_xml_format
    MODULE PROCEDURE qes_init_creator
    MODULE PROCEDURE qes_init_created
    MODULE PROCEDURE qes_init_atomic_species
    MODULE PROCEDURE qes_init_species
    MODULE PROCEDURE qes_init_atomic_structure
    MODULE PROCEDURE qes_init_atomic_positions
    MODULE PROCEDURE qes_init_atom
    MODULE PROCEDURE qes_init_wyckoff_positions
    MODULE PROCEDURE qes_init_cell
    MODULE PROCEDURE qes_init_dft
    MODULE PROCEDURE qes_init_hybrid
    MODULE PROCEDURE qes_init_qpoint_grid
    MODULE PROCEDURE qes_init_dftU
    MODULE PROCEDURE qes_init_HubbardCommon
    MODULE PROCEDURE qes_init_HubbardJ
    MODULE PROCEDURE qes_init_starting_ns
    MODULE PROCEDURE qes_init_Hubbard_ns
    MODULE PROCEDURE qes_init_vdW
    MODULE PROCEDURE qes_init_spin
    MODULE PROCEDURE qes_init_bands
    MODULE PROCEDURE qes_init_smearing
    MODULE PROCEDURE qes_init_occupations
    MODULE PROCEDURE qes_init_basis
    MODULE PROCEDURE qes_init_basis_set
    MODULE PROCEDURE qes_init_basisSetItem
    MODULE PROCEDURE qes_init_reciprocal_lattice
    MODULE PROCEDURE qes_init_electron_control
    MODULE PROCEDURE qes_init_k_points_IBZ
    MODULE PROCEDURE qes_init_monkhorst_pack
    MODULE PROCEDURE qes_init_k_point
    MODULE PROCEDURE qes_init_ion_control
    MODULE PROCEDURE qes_init_bfgs
    MODULE PROCEDURE qes_init_md
    MODULE PROCEDURE qes_init_cell_control
    MODULE PROCEDURE qes_init_symmetry_flags
    MODULE PROCEDURE qes_init_boundary_conditions
    MODULE PROCEDURE qes_init_esm
    MODULE PROCEDURE qes_init_ekin_functional
    MODULE PROCEDURE qes_init_spin_constraints
    MODULE PROCEDURE qes_init_electric_field
    MODULE PROCEDURE qes_init_gate_settings
    MODULE PROCEDURE qes_init_atomic_constraints
    MODULE PROCEDURE qes_init_atomic_constraint
    MODULE PROCEDURE qes_init_inputOccupations
    MODULE PROCEDURE qes_init_outputElectricField
    MODULE PROCEDURE qes_init_BerryPhaseOutput
    MODULE PROCEDURE qes_init_dipoleOutput
    MODULE PROCEDURE qes_init_finiteFieldOut
    MODULE PROCEDURE qes_init_polarization
    MODULE PROCEDURE qes_init_ionicPolarization
    MODULE PROCEDURE qes_init_electronicPolarization
    MODULE PROCEDURE qes_init_phase
    MODULE PROCEDURE qes_init_gateInfo
    MODULE PROCEDURE qes_init_convergence_info
    MODULE PROCEDURE qes_init_scf_conv
    MODULE PROCEDURE qes_init_opt_conv
    MODULE PROCEDURE qes_init_algorithmic_info
    MODULE PROCEDURE qes_init_symmetries
    MODULE PROCEDURE qes_init_symmetry
    MODULE PROCEDURE qes_init_equivalent_atoms
    MODULE PROCEDURE qes_init_info
    MODULE PROCEDURE qes_init_outputPBC
    MODULE PROCEDURE qes_init_magnetization
    MODULE PROCEDURE qes_init_total_energy
    MODULE PROCEDURE qes_init_band_structure
    MODULE PROCEDURE qes_init_ks_energies
    MODULE PROCEDURE qes_init_closed
    MODULE PROCEDURE qes_init_vector
    MODULE PROCEDURE qes_init_integerVector
    MODULE PROCEDURE qes_init_matrix_1, qes_init_matrix_2, qes_init_matrix_3
    MODULE PROCEDURE qes_init_integerMatrix_1, qes_init_integerMatrix_2, qes_init_integerMatrix_3
    MODULE PROCEDURE qes_init_scalarQuantity
    !
  END INTERFACE qes_init
  !
  CONTAINS
  !
  !
  SUBROUTINE qes_init_espresso(obj, tagname, input, Units, general_info, parallel_info, step,&
                              output, status, cputime, timing_info, closed)
    !
    IMPLICIT NONE
    !
    TYPE(espresso_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: Units
    TYPE(general_info_type),OPTIONAL,INTENT(IN) :: general_info
    TYPE(parallel_info_type),OPTIONAL,INTENT(IN) :: parallel_info
    TYPE(input_type),INTENT(IN) :: input
    TYPE(step_type),OPTIONAL,DIMENSION(:),INTENT(IN) :: step
    TYPE(output_type),OPTIONAL,INTENT(IN) :: output
    INTEGER,OPTIONAL,INTENT(IN) :: status
    INTEGER,OPTIONAL,INTENT(IN) :: cputime
    TYPE(timing_type),OPTIONAL,INTENT(IN) :: timing_info
    TYPE(closed_type),OPTIONAL,INTENT(IN) :: closed
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    IF (PRESENT(Units)) THEN
      obj%Units_ispresent = .TRUE.
      obj%Units = Units
    ELSE 
      obj%Units_ispresent = .FALSE.
    END IF
    !
    IF ( PRESENT(general_info)) THEN 
      obj%general_info_ispresent = .TRUE. 
      obj%general_info = general_info
    ELSE 
      obj%general_info_ispresent = .FALSE.
    END IF
    IF ( PRESENT(parallel_info)) THEN 
      obj%parallel_info_ispresent = .TRUE. 
      obj%parallel_info = parallel_info
    ELSE 
      obj%parallel_info_ispresent = .FALSE.
    END IF
    obj%input = input
    IF ( PRESENT(step)) THEN 
      obj%step_ispresent = .TRUE.
      ALLOCATE(obj%step(SIZE(step)))
      obj%ndim_step = SIZE(step) 
      obj%step = step
    ELSE 
      obj%step_ispresent = .FALSE.
    END IF
    IF ( PRESENT(output)) THEN 
      obj%output_ispresent = .TRUE. 
      obj%output = output
    ELSE 
      obj%output_ispresent = .FALSE.
    END IF
    IF ( PRESENT(status)) THEN 
      obj%status_ispresent = .TRUE. 
      obj%status = status
    ELSE 
      obj%status_ispresent = .FALSE.
    END IF
    IF ( PRESENT(cputime)) THEN 
      obj%cputime_ispresent = .TRUE. 
      obj%cputime = cputime
    ELSE 
      obj%cputime_ispresent = .FALSE.
    END IF
    IF ( PRESENT(timing_info)) THEN 
      obj%timing_info_ispresent = .TRUE. 
      obj%timing_info = timing_info
    ELSE 
      obj%timing_info_ispresent = .FALSE.
    END IF
    IF ( PRESENT(closed)) THEN 
      obj%closed_ispresent = .TRUE. 
      obj%closed = closed
    ELSE 
      obj%closed_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_espresso 
  !
  !
  SUBROUTINE qes_init_general_info(obj, tagname, xml_format, creator, created, job)
    !
    IMPLICIT NONE
    !
    TYPE(general_info_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(xml_format_type),INTENT(IN) :: xml_format
    TYPE(creator_type),INTENT(IN) :: creator
    TYPE(created_type),INTENT(IN) :: created
    CHARACTER(LEN=*),INTENT(IN) :: job
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%xml_format = xml_format
    obj%creator = creator
    obj%created = created
    obj%job = job
    !
  END SUBROUTINE qes_init_general_info 
  !
  !
  SUBROUTINE qes_init_parallel_info(obj, tagname, nprocs, nthreads, ntasks, nbgrp, npool, ndiag)
    !
    IMPLICIT NONE
    !
    TYPE(parallel_info_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,INTENT(IN) :: nprocs
    INTEGER,INTENT(IN) :: nthreads
    INTEGER,INTENT(IN) :: ntasks
    INTEGER,INTENT(IN) :: nbgrp
    INTEGER,INTENT(IN) :: npool
    INTEGER,INTENT(IN) :: ndiag
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%nprocs = nprocs
    obj%nthreads = nthreads
    obj%ntasks = ntasks
    obj%nbgrp = nbgrp
    obj%npool = npool
    obj%ndiag = ndiag
    !
  END SUBROUTINE qes_init_parallel_info 
  !
  !
  SUBROUTINE qes_init_input(obj, tagname, control_variables, atomic_species, atomic_structure,&
                           dft, spin, bands, basis, electron_control, k_points_IBZ, ion_control,&
                           cell_control, symmetry_flags, boundary_conditions, ekin_functional,&
                           external_atomic_forces, free_positions, starting_atomic_velocities,&
                           electric_field, atomic_constraints, spin_constraints)
    !
    IMPLICIT NONE
    !
    TYPE(input_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(control_variables_type),INTENT(IN) :: control_variables
    TYPE(atomic_species_type),INTENT(IN) :: atomic_species
    TYPE(atomic_structure_type),INTENT(IN) :: atomic_structure
    TYPE(dft_type),INTENT(IN) :: dft
    TYPE(spin_type),INTENT(IN) :: spin
    TYPE(bands_type),INTENT(IN) :: bands
    TYPE(basis_type),INTENT(IN) :: basis
    TYPE(electron_control_type),INTENT(IN) :: electron_control
    TYPE(k_points_IBZ_type),INTENT(IN) :: k_points_IBZ
    TYPE(ion_control_type),INTENT(IN) :: ion_control
    TYPE(cell_control_type),INTENT(IN) :: cell_control
    TYPE(symmetry_flags_type),OPTIONAL,INTENT(IN) :: symmetry_flags
    TYPE(boundary_conditions_type),OPTIONAL,INTENT(IN) :: boundary_conditions
    TYPE(ekin_functional_type),OPTIONAL,INTENT(IN) :: ekin_functional
    TYPE(matrix_type),OPTIONAL,INTENT(IN) :: external_atomic_forces
    TYPE(integerMatrix_type),OPTIONAL,INTENT(IN) :: free_positions
    TYPE(matrix_type),OPTIONAL,INTENT(IN) :: starting_atomic_velocities
    TYPE(electric_field_type),OPTIONAL,INTENT(IN) :: electric_field
    TYPE(atomic_constraints_type),OPTIONAL,INTENT(IN) :: atomic_constraints
    TYPE(spin_constraints_type),OPTIONAL,INTENT(IN) :: spin_constraints
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%control_variables = control_variables
    obj%atomic_species = atomic_species
    obj%atomic_structure = atomic_structure
    obj%dft = dft
    obj%spin = spin
    obj%bands = bands
    obj%basis = basis
    obj%electron_control = electron_control
    obj%k_points_IBZ = k_points_IBZ
    obj%ion_control = ion_control
    obj%cell_control = cell_control
    IF ( PRESENT(symmetry_flags)) THEN 
      obj%symmetry_flags_ispresent = .TRUE. 
      obj%symmetry_flags = symmetry_flags
    ELSE 
      obj%symmetry_flags_ispresent = .FALSE.
    END IF
    IF ( PRESENT(boundary_conditions)) THEN 
      obj%boundary_conditions_ispresent = .TRUE. 
      obj%boundary_conditions = boundary_conditions
    ELSE 
      obj%boundary_conditions_ispresent = .FALSE.
    END IF
    IF ( PRESENT(ekin_functional)) THEN 
      obj%ekin_functional_ispresent = .TRUE. 
      obj%ekin_functional = ekin_functional
    ELSE 
      obj%ekin_functional_ispresent = .FALSE.
    END IF
    IF ( PRESENT(external_atomic_forces)) THEN 
      obj%external_atomic_forces_ispresent = .TRUE. 
      obj%external_atomic_forces = external_atomic_forces
    ELSE 
      obj%external_atomic_forces_ispresent = .FALSE.
    END IF
    IF ( PRESENT(free_positions)) THEN 
      obj%free_positions_ispresent = .TRUE. 
      obj%free_positions = free_positions
    ELSE 
      obj%free_positions_ispresent = .FALSE.
    END IF
    IF ( PRESENT(starting_atomic_velocities)) THEN 
      obj%starting_atomic_velocities_ispresent = .TRUE. 
      obj%starting_atomic_velocities = starting_atomic_velocities
    ELSE 
      obj%starting_atomic_velocities_ispresent = .FALSE.
    END IF
    IF ( PRESENT(electric_field)) THEN 
      obj%electric_field_ispresent = .TRUE. 
      obj%electric_field = electric_field
    ELSE 
      obj%electric_field_ispresent = .FALSE.
    END IF
    IF ( PRESENT(atomic_constraints)) THEN 
      obj%atomic_constraints_ispresent = .TRUE. 
      obj%atomic_constraints = atomic_constraints
    ELSE 
      obj%atomic_constraints_ispresent = .FALSE.
    END IF
    IF ( PRESENT(spin_constraints)) THEN 
      obj%spin_constraints_ispresent = .TRUE. 
      obj%spin_constraints = spin_constraints
    ELSE 
      obj%spin_constraints_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_input 
  !
  !
  SUBROUTINE qes_init_step(obj, tagname, n_step, scf_conv, atomic_structure, total_energy, forces,&
                          stress, FCP_force, FCP_tot_charge)
    !
    IMPLICIT NONE
    !
    TYPE(step_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER, INTENT(IN) :: n_step
    TYPE(scf_conv_type),INTENT(IN) :: scf_conv
    TYPE(atomic_structure_type),INTENT(IN) :: atomic_structure
    TYPE(total_energy_type),INTENT(IN) :: total_energy
    TYPE(matrix_type),INTENT(IN) :: forces
    TYPE(matrix_type),OPTIONAL,INTENT(IN) :: stress
    REAL(DP),OPTIONAL,INTENT(IN) :: FCP_force
    REAL(DP),OPTIONAL,INTENT(IN) :: FCP_tot_charge
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%n_step = n_step
    !
    obj%scf_conv = scf_conv
    obj%atomic_structure = atomic_structure
    obj%total_energy = total_energy
    obj%forces = forces
    IF ( PRESENT(stress)) THEN 
      obj%stress_ispresent = .TRUE. 
      obj%stress = stress
    ELSE 
      obj%stress_ispresent = .FALSE.
    END IF
    IF ( PRESENT(FCP_force)) THEN 
      obj%FCP_force_ispresent = .TRUE. 
      obj%FCP_force = FCP_force
    ELSE 
      obj%FCP_force_ispresent = .FALSE.
    END IF
    IF ( PRESENT(FCP_tot_charge)) THEN 
      obj%FCP_tot_charge_ispresent = .TRUE. 
      obj%FCP_tot_charge = FCP_tot_charge
    ELSE 
      obj%FCP_tot_charge_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_step 
  !
  !
  SUBROUTINE qes_init_output(obj, tagname, algorithmic_info, atomic_species, atomic_structure,&
                            basis_set, dft, magnetization, total_energy, band_structure, convergence_info,&
                            symmetries, boundary_conditions, forces, stress, electric_field,&
                            FCP_force, FCP_tot_charge)
    !
    IMPLICIT NONE
    !
    TYPE(output_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(convergence_info_type),OPTIONAL,INTENT(IN) :: convergence_info
    TYPE(algorithmic_info_type),INTENT(IN) :: algorithmic_info
    TYPE(atomic_species_type),INTENT(IN) :: atomic_species
    TYPE(atomic_structure_type),INTENT(IN) :: atomic_structure
    TYPE(symmetries_type),OPTIONAL,INTENT(IN) :: symmetries
    TYPE(basis_set_type),INTENT(IN) :: basis_set
    TYPE(dft_type),INTENT(IN) :: dft
    TYPE(outputPBC_type),OPTIONAL,INTENT(IN) :: boundary_conditions
    TYPE(magnetization_type),INTENT(IN) :: magnetization
    TYPE(total_energy_type),INTENT(IN) :: total_energy
    TYPE(band_structure_type),INTENT(IN) :: band_structure
    TYPE(matrix_type),OPTIONAL,INTENT(IN) :: forces
    TYPE(matrix_type),OPTIONAL,INTENT(IN) :: stress
    TYPE(outputElectricField_type),OPTIONAL,INTENT(IN) :: electric_field
    REAL(DP),OPTIONAL,INTENT(IN) :: FCP_force
    REAL(DP),OPTIONAL,INTENT(IN) :: FCP_tot_charge
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(convergence_info)) THEN 
      obj%convergence_info_ispresent = .TRUE. 
      obj%convergence_info = convergence_info
    ELSE 
      obj%convergence_info_ispresent = .FALSE.
    END IF
    obj%algorithmic_info = algorithmic_info
    obj%atomic_species = atomic_species
    obj%atomic_structure = atomic_structure
    IF ( PRESENT(symmetries)) THEN 
      obj%symmetries_ispresent = .TRUE. 
      obj%symmetries = symmetries
    ELSE 
      obj%symmetries_ispresent = .FALSE.
    END IF
    obj%basis_set = basis_set
    obj%dft = dft
    IF ( PRESENT(boundary_conditions)) THEN 
      obj%boundary_conditions_ispresent = .TRUE. 
      obj%boundary_conditions = boundary_conditions
    ELSE 
      obj%boundary_conditions_ispresent = .FALSE.
    END IF
    obj%magnetization = magnetization
    obj%total_energy = total_energy
    obj%band_structure = band_structure
    IF ( PRESENT(forces)) THEN 
      obj%forces_ispresent = .TRUE. 
      obj%forces = forces
    ELSE 
      obj%forces_ispresent = .FALSE.
    END IF
    IF ( PRESENT(stress)) THEN 
      obj%stress_ispresent = .TRUE. 
      obj%stress = stress
    ELSE 
      obj%stress_ispresent = .FALSE.
    END IF
    IF ( PRESENT(electric_field)) THEN 
      obj%electric_field_ispresent = .TRUE. 
      obj%electric_field = electric_field
    ELSE 
      obj%electric_field_ispresent = .FALSE.
    END IF
    IF ( PRESENT(FCP_force)) THEN 
      obj%FCP_force_ispresent = .TRUE. 
      obj%FCP_force = FCP_force
    ELSE 
      obj%FCP_force_ispresent = .FALSE.
    END IF
    IF ( PRESENT(FCP_tot_charge)) THEN 
      obj%FCP_tot_charge_ispresent = .TRUE. 
      obj%FCP_tot_charge = FCP_tot_charge
    ELSE 
      obj%FCP_tot_charge_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_output 
  !
  !
  SUBROUTINE qes_init_timing(obj, tagname, total, partial)
    !
    IMPLICIT NONE
    !
    TYPE(timing_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(clock_type),INTENT(IN) :: total
    TYPE(clock_type),OPTIONAL,DIMENSION(:),INTENT(IN) :: partial
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%total = total
    IF ( PRESENT(partial)) THEN 
      obj%partial_ispresent = .TRUE.
      ALLOCATE(obj%partial(SIZE(partial)))
      obj%ndim_partial = SIZE(partial) 
      obj%partial = partial
    ELSE 
      obj%partial_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_timing 
  !
  !
  SUBROUTINE qes_init_clock(obj, tagname, label, cpu, wall, calls)
    !
    IMPLICIT NONE
    !
    TYPE(clock_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, OPTIONAL, INTENT(IN) :: calls
    REAL(DP),INTENT(IN) :: cpu
    REAL(DP),INTENT(IN) :: wall
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%label = label
    IF (PRESENT(calls)) THEN
      obj%calls_ispresent = .TRUE.
      obj%calls = calls
    ELSE 
      obj%calls_ispresent = .FALSE.
    END IF
    !
    obj%cpu = cpu
    obj%wall = wall
    !
  END SUBROUTINE qes_init_clock 
  !
  !
  SUBROUTINE qes_init_control_variables(obj, tagname, title, calculation, restart_mode, prefix,&
                                       pseudo_dir, outdir, stress, forces, wf_collect, disk_io,&
                                       max_seconds, etot_conv_thr, forc_conv_thr, press_conv_thr,&
                                       verbosity, print_every, nstep)
    !
    IMPLICIT NONE
    !
    TYPE(control_variables_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),INTENT(IN) :: title
    CHARACTER(LEN=*),INTENT(IN) :: calculation
    CHARACTER(LEN=*),INTENT(IN) :: restart_mode
    CHARACTER(LEN=*),INTENT(IN) :: prefix
    CHARACTER(LEN=*),INTENT(IN) :: pseudo_dir
    CHARACTER(LEN=*),INTENT(IN) :: outdir
    LOGICAL,INTENT(IN) :: stress
    LOGICAL,INTENT(IN) :: forces
    LOGICAL,INTENT(IN) :: wf_collect
    CHARACTER(LEN=*),INTENT(IN) :: disk_io
    INTEGER,INTENT(IN) :: max_seconds
    INTEGER,OPTIONAL,INTENT(IN) :: nstep
    REAL(DP),INTENT(IN) :: etot_conv_thr
    REAL(DP),INTENT(IN) :: forc_conv_thr
    REAL(DP),INTENT(IN) :: press_conv_thr
    CHARACTER(LEN=*),INTENT(IN) :: verbosity
    INTEGER,INTENT(IN) :: print_every
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%title = title
    obj%calculation = calculation
    obj%restart_mode = restart_mode
    obj%prefix = prefix
    obj%pseudo_dir = pseudo_dir
    obj%outdir = outdir
    obj%stress = stress
    obj%forces = forces
    obj%wf_collect = wf_collect
    obj%disk_io = disk_io
    obj%max_seconds = max_seconds
    IF ( PRESENT(nstep)) THEN 
      obj%nstep_ispresent = .TRUE. 
      obj%nstep = nstep
    ELSE 
      obj%nstep_ispresent = .FALSE.
    END IF
    obj%etot_conv_thr = etot_conv_thr
    obj%forc_conv_thr = forc_conv_thr
    obj%press_conv_thr = press_conv_thr
    obj%verbosity = verbosity
    obj%print_every = print_every
    !
  END SUBROUTINE qes_init_control_variables 
  !
  !
  SUBROUTINE qes_init_xml_format(obj, tagname, NAME, VERSION, xml_format)
    !
    IMPLICIT NONE
    !
    TYPE(xml_format_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    CHARACTER(LEN=*), INTENT(IN) :: VERSION
    CHARACTER(LEN=*), INTENT(IN) :: xml_format
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%NAME = NAME
    obj%VERSION = VERSION
    !
    obj%xml_format = xml_format
    !
  END SUBROUTINE qes_init_xml_format 
  !
  !
  SUBROUTINE qes_init_creator(obj, tagname, NAME, VERSION, creator)
    !
    IMPLICIT NONE
    !
    TYPE(creator_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    CHARACTER(LEN=*), INTENT(IN) :: VERSION
    CHARACTER(LEN=*), INTENT(IN) :: creator
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%NAME = NAME
    obj%VERSION = VERSION
    !
    obj%creator = creator
    !
  END SUBROUTINE qes_init_creator 
  !
  !
  SUBROUTINE qes_init_created(obj, tagname, DATE, TIME, created)
    !
    IMPLICIT NONE
    !
    TYPE(created_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: DATE
    CHARACTER(LEN=*), INTENT(IN) :: TIME
    CHARACTER(LEN=*), INTENT(IN) :: created
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%DATE = DATE
    obj%TIME = TIME
    !
    obj%created = created
    !
  END SUBROUTINE qes_init_created 
  !
  !
  SUBROUTINE qes_init_atomic_species(obj, tagname, ntyp, species, pseudo_dir)
    !
    IMPLICIT NONE
    !
    TYPE(atomic_species_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER, INTENT(IN) :: ntyp
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: pseudo_dir
    TYPE(species_type),DIMENSION(:),INTENT(IN) :: species
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%ntyp = ntyp
    IF (PRESENT(pseudo_dir)) THEN
      obj%pseudo_dir_ispresent = .TRUE.
      obj%pseudo_dir = pseudo_dir
    ELSE 
      obj%pseudo_dir_ispresent = .FALSE.
    END IF
    !
    ALLOCATE( obj%species(SIZE(species))) 
    obj%ndim_species = SIZE(species)
    obj%species = species
    !
  END SUBROUTINE qes_init_atomic_species 
  !
  !
  SUBROUTINE qes_init_species(obj, tagname, name, pseudo_file, mass, starting_magnetization, spin_teta, spin_phi)
    !
    IMPLICIT NONE
    !
    TYPE(species_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: name
    REAL(DP),OPTIONAL,INTENT(IN) :: mass
    CHARACTER(LEN=*),INTENT(IN) :: pseudo_file
    REAL(DP),OPTIONAL,INTENT(IN) :: starting_magnetization
    REAL(DP),OPTIONAL,INTENT(IN) :: spin_teta
    REAL(DP),OPTIONAL,INTENT(IN) :: spin_phi
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%name = name
    !
    IF ( PRESENT(mass)) THEN 
      obj%mass_ispresent = .TRUE. 
      obj%mass = mass
    ELSE 
      obj%mass_ispresent = .FALSE.
    END IF
    obj%pseudo_file = pseudo_file
    IF ( PRESENT(starting_magnetization)) THEN 
      obj%starting_magnetization_ispresent = .TRUE. 
      obj%starting_magnetization = starting_magnetization
    ELSE 
      obj%starting_magnetization_ispresent = .FALSE.
    END IF
    IF ( PRESENT(spin_teta)) THEN 
      obj%spin_teta_ispresent = .TRUE. 
      obj%spin_teta = spin_teta
    ELSE 
      obj%spin_teta_ispresent = .FALSE.
    END IF
    IF ( PRESENT(spin_phi)) THEN 
      obj%spin_phi_ispresent = .TRUE. 
      obj%spin_phi = spin_phi
    ELSE 
      obj%spin_phi_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_species 
  !
  !
  SUBROUTINE qes_init_atomic_structure(obj, tagname, nat, cell, alat, bravais_index, atomic_positions,&
                                      wyckoff_positions, crystal_positions)
    !
    IMPLICIT NONE
    !
    TYPE(atomic_structure_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER, INTENT(IN) :: nat
    REAL(DP), OPTIONAL, INTENT(IN) :: alat
    INTEGER, OPTIONAL, INTENT(IN) :: bravais_index
    TYPE(atomic_positions_type),OPTIONAL,INTENT(IN) :: atomic_positions
    TYPE(wyckoff_positions_type),OPTIONAL,INTENT(IN) :: wyckoff_positions
    TYPE(atomic_positions_type),OPTIONAL,INTENT(IN) :: crystal_positions
    TYPE(cell_type),INTENT(IN) :: cell
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%nat = nat
    IF (PRESENT(alat)) THEN
      obj%alat_ispresent = .TRUE.
      obj%alat = alat
    ELSE 
      obj%alat_ispresent = .FALSE.
    END IF
    IF (PRESENT(bravais_index)) THEN
      obj%bravais_index_ispresent = .TRUE.
      obj%bravais_index = bravais_index
    ELSE 
      obj%bravais_index_ispresent = .FALSE.
    END IF
    !
    IF ( PRESENT(atomic_positions)) THEN 
      obj%atomic_positions_ispresent = .TRUE. 
      obj%atomic_positions = atomic_positions
    ELSE 
      obj%atomic_positions_ispresent = .FALSE.
    END IF
    IF ( PRESENT(wyckoff_positions)) THEN 
      obj%wyckoff_positions_ispresent = .TRUE. 
      obj%wyckoff_positions = wyckoff_positions
    ELSE 
      obj%wyckoff_positions_ispresent = .FALSE.
    END IF
    IF ( PRESENT(crystal_positions)) THEN 
      obj%crystal_positions_ispresent = .TRUE. 
      obj%crystal_positions = crystal_positions
    ELSE 
      obj%crystal_positions_ispresent = .FALSE.
    END IF
    obj%cell = cell
    !
  END SUBROUTINE qes_init_atomic_structure 
  !
  !
  SUBROUTINE qes_init_atomic_positions(obj, tagname, atom)
    !
    IMPLICIT NONE
    !
    TYPE(atomic_positions_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(atom_type),DIMENSION(:),INTENT(IN) :: atom
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    ALLOCATE( obj%atom(SIZE(atom))) 
    obj%ndim_atom = SIZE(atom)
    obj%atom = atom
    !
  END SUBROUTINE qes_init_atomic_positions 
  !
  !
  SUBROUTINE qes_init_atom(obj, tagname, name, atom, position, index)
    !
    IMPLICIT NONE
    !
    TYPE(atom_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: name
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: position
    INTEGER, OPTIONAL, INTENT(IN) :: index
    REAL(DP), DIMENSION(3), INTENT(IN) :: atom
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%name = name
    IF (PRESENT(position)) THEN
      obj%position_ispresent = .TRUE.
      obj%position = position
    ELSE 
      obj%position_ispresent = .FALSE.
    END IF
    IF (PRESENT(index)) THEN
      obj%index_ispresent = .TRUE.
      obj%index = index
    ELSE 
      obj%index_ispresent = .FALSE.
    END IF
    !
    obj%atom = atom
    !
  END SUBROUTINE qes_init_atom 
  !
  !
  SUBROUTINE qes_init_wyckoff_positions(obj, tagname, space_group, atom, more_options)
    !
    IMPLICIT NONE
    !
    TYPE(wyckoff_positions_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER, INTENT(IN) :: space_group
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: more_options
    TYPE(atom_type),DIMENSION(:),INTENT(IN) :: atom
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%space_group = space_group
    IF (PRESENT(more_options)) THEN
      obj%more_options_ispresent = .TRUE.
      obj%more_options = more_options
    ELSE 
      obj%more_options_ispresent = .FALSE.
    END IF
    !
    ALLOCATE( obj%atom(SIZE(atom))) 
    obj%ndim_atom = SIZE(atom)
    obj%atom = atom
    !
  END SUBROUTINE qes_init_wyckoff_positions 
  !
  !
  SUBROUTINE qes_init_cell(obj, tagname, a1, a2, a3)
    !
    IMPLICIT NONE
    !
    TYPE(cell_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP), DIMENSION(3),INTENT(IN) :: a1
    REAL(DP), DIMENSION(3),INTENT(IN) :: a2
    REAL(DP), DIMENSION(3),INTENT(IN) :: a3
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%a1 = a1
    obj%a2 = a2
    obj%a3 = a3
    !
  END SUBROUTINE qes_init_cell 
  !
  !
  SUBROUTINE qes_init_dft(obj, tagname, functional, hybrid, dftU, vdW)
    !
    IMPLICIT NONE
    !
    TYPE(dft_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),INTENT(IN) :: functional
    TYPE(hybrid_type),OPTIONAL,INTENT(IN) :: hybrid
    TYPE(dftU_type),OPTIONAL,INTENT(IN) :: dftU
    TYPE(vdW_type),OPTIONAL,INTENT(IN) :: vdW
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%functional = functional
    IF ( PRESENT(hybrid)) THEN 
      obj%hybrid_ispresent = .TRUE. 
      obj%hybrid = hybrid
    ELSE 
      obj%hybrid_ispresent = .FALSE.
    END IF
    IF ( PRESENT(dftU)) THEN 
      obj%dftU_ispresent = .TRUE. 
      obj%dftU = dftU
    ELSE 
      obj%dftU_ispresent = .FALSE.
    END IF
    IF ( PRESENT(vdW)) THEN 
      obj%vdW_ispresent = .TRUE. 
      obj%vdW = vdW
    ELSE 
      obj%vdW_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_dft 
  !
  !
  SUBROUTINE qes_init_hybrid(obj, tagname, qpoint_grid, ecutfock, exx_fraction, screening_parameter,&
                            exxdiv_treatment, x_gamma_extrapolation, ecutvcut, localization_threshold &
                            )
    !
    IMPLICIT NONE
    !
    TYPE(hybrid_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(qpoint_grid_type),OPTIONAL,INTENT(IN) :: qpoint_grid
    REAL(DP),OPTIONAL,INTENT(IN) :: ecutfock
    REAL(DP),OPTIONAL,INTENT(IN) :: exx_fraction
    REAL(DP),OPTIONAL,INTENT(IN) :: screening_parameter
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: exxdiv_treatment
    LOGICAL,OPTIONAL,INTENT(IN) :: x_gamma_extrapolation
    REAL(DP),OPTIONAL,INTENT(IN) :: ecutvcut
    REAL(DP),OPTIONAL,INTENT(IN) :: localization_threshold
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(qpoint_grid)) THEN 
      obj%qpoint_grid_ispresent = .TRUE. 
      obj%qpoint_grid = qpoint_grid
    ELSE 
      obj%qpoint_grid_ispresent = .FALSE.
    END IF
    IF ( PRESENT(ecutfock)) THEN 
      obj%ecutfock_ispresent = .TRUE. 
      obj%ecutfock = ecutfock
    ELSE 
      obj%ecutfock_ispresent = .FALSE.
    END IF
    IF ( PRESENT(exx_fraction)) THEN 
      obj%exx_fraction_ispresent = .TRUE. 
      obj%exx_fraction = exx_fraction
    ELSE 
      obj%exx_fraction_ispresent = .FALSE.
    END IF
    IF ( PRESENT(screening_parameter)) THEN 
      obj%screening_parameter_ispresent = .TRUE. 
      obj%screening_parameter = screening_parameter
    ELSE 
      obj%screening_parameter_ispresent = .FALSE.
    END IF
    IF ( PRESENT(exxdiv_treatment)) THEN 
      obj%exxdiv_treatment_ispresent = .TRUE. 
      obj%exxdiv_treatment = exxdiv_treatment
    ELSE 
      obj%exxdiv_treatment_ispresent = .FALSE.
    END IF
    IF ( PRESENT(x_gamma_extrapolation)) THEN 
      obj%x_gamma_extrapolation_ispresent = .TRUE. 
      obj%x_gamma_extrapolation = x_gamma_extrapolation
    ELSE 
      obj%x_gamma_extrapolation_ispresent = .FALSE.
    END IF
    IF ( PRESENT(ecutvcut)) THEN 
      obj%ecutvcut_ispresent = .TRUE. 
      obj%ecutvcut = ecutvcut
    ELSE 
      obj%ecutvcut_ispresent = .FALSE.
    END IF
    IF ( PRESENT(localization_threshold)) THEN 
      obj%localization_threshold_ispresent = .TRUE. 
      obj%localization_threshold = localization_threshold
    ELSE 
      obj%localization_threshold_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_hybrid 
  !
  !
  SUBROUTINE qes_init_qpoint_grid(obj, tagname, nqx1, nqx2, nqx3, qpoint_grid)
    !
    IMPLICIT NONE
    !
    TYPE(qpoint_grid_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER, INTENT(IN) :: nqx1
    INTEGER, INTENT(IN) :: nqx2
    INTEGER, INTENT(IN) :: nqx3
    CHARACTER(LEN=*), INTENT(IN) :: qpoint_grid
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%nqx1 = nqx1
    obj%nqx2 = nqx2
    obj%nqx3 = nqx3
    !
    obj%qpoint_grid = qpoint_grid
    !
  END SUBROUTINE qes_init_qpoint_grid 
  !
  !
  SUBROUTINE qes_init_dftU(obj, tagname, lda_plus_u_kind, Hubbard_U, Hubbard_J0, Hubbard_alpha,&
                          Hubbard_beta, Hubbard_J, starting_ns, Hubbard_ns, U_projection_type)
    !
    IMPLICIT NONE
    !
    TYPE(dftU_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,OPTIONAL,INTENT(IN) :: lda_plus_u_kind
    TYPE(HubbardCommon_type),OPTIONAL,DIMENSION(:),INTENT(IN) :: Hubbard_U
    TYPE(HubbardCommon_type),OPTIONAL,DIMENSION(:),INTENT(IN) :: Hubbard_J0
    TYPE(HubbardCommon_type),OPTIONAL,DIMENSION(:),INTENT(IN) :: Hubbard_alpha
    TYPE(HubbardCommon_type),OPTIONAL,DIMENSION(:),INTENT(IN) :: Hubbard_beta
    TYPE(HubbardJ_type),OPTIONAL,DIMENSION(:),INTENT(IN) :: Hubbard_J
    TYPE(starting_ns_type),OPTIONAL,DIMENSION(:),INTENT(IN) :: starting_ns
    TYPE(Hubbard_ns_type),OPTIONAL,DIMENSION(:),INTENT(IN) :: Hubbard_ns
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: U_projection_type
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(lda_plus_u_kind)) THEN 
      obj%lda_plus_u_kind_ispresent = .TRUE. 
      obj%lda_plus_u_kind = lda_plus_u_kind
    ELSE 
      obj%lda_plus_u_kind_ispresent = .FALSE.
    END IF
    IF ( PRESENT(Hubbard_U)) THEN 
      obj%Hubbard_U_ispresent = .TRUE.
      ALLOCATE(obj%Hubbard_U(SIZE(Hubbard_U)))
      obj%ndim_Hubbard_U = SIZE(Hubbard_U) 
      obj%Hubbard_U = Hubbard_U
    ELSE 
      obj%Hubbard_U_ispresent = .FALSE.
    END IF
    IF ( PRESENT(Hubbard_J0)) THEN 
      obj%Hubbard_J0_ispresent = .TRUE.
      ALLOCATE(obj%Hubbard_J0(SIZE(Hubbard_J0)))
      obj%ndim_Hubbard_J0 = SIZE(Hubbard_J0) 
      obj%Hubbard_J0 = Hubbard_J0
    ELSE 
      obj%Hubbard_J0_ispresent = .FALSE.
    END IF
    IF ( PRESENT(Hubbard_alpha)) THEN 
      obj%Hubbard_alpha_ispresent = .TRUE.
      ALLOCATE(obj%Hubbard_alpha(SIZE(Hubbard_alpha)))
      obj%ndim_Hubbard_alpha = SIZE(Hubbard_alpha) 
      obj%Hubbard_alpha = Hubbard_alpha
    ELSE 
      obj%Hubbard_alpha_ispresent = .FALSE.
    END IF
    IF ( PRESENT(Hubbard_beta)) THEN 
      obj%Hubbard_beta_ispresent = .TRUE.
      ALLOCATE(obj%Hubbard_beta(SIZE(Hubbard_beta)))
      obj%ndim_Hubbard_beta = SIZE(Hubbard_beta) 
      obj%Hubbard_beta = Hubbard_beta
    ELSE 
      obj%Hubbard_beta_ispresent = .FALSE.
    END IF
    IF ( PRESENT(Hubbard_J)) THEN 
      obj%Hubbard_J_ispresent = .TRUE.
      ALLOCATE(obj%Hubbard_J(SIZE(Hubbard_J)))
      obj%ndim_Hubbard_J = SIZE(Hubbard_J) 
      obj%Hubbard_J = Hubbard_J
    ELSE 
      obj%Hubbard_J_ispresent = .FALSE.
    END IF
    IF ( PRESENT(starting_ns)) THEN 
      obj%starting_ns_ispresent = .TRUE.
      ALLOCATE(obj%starting_ns(SIZE(starting_ns)))
      obj%ndim_starting_ns = SIZE(starting_ns) 
      obj%starting_ns = starting_ns
    ELSE 
      obj%starting_ns_ispresent = .FALSE.
    END IF
    IF ( PRESENT(Hubbard_ns)) THEN 
      obj%Hubbard_ns_ispresent = .TRUE.
      ALLOCATE(obj%Hubbard_ns(SIZE(Hubbard_ns)))
      obj%ndim_Hubbard_ns = SIZE(Hubbard_ns) 
      obj%Hubbard_ns = Hubbard_ns
    ELSE 
      obj%Hubbard_ns_ispresent = .FALSE.
    END IF
    IF ( PRESENT(U_projection_type)) THEN 
      obj%U_projection_type_ispresent = .TRUE. 
      obj%U_projection_type = U_projection_type
    ELSE 
      obj%U_projection_type_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_dftU 
  !
  !
  SUBROUTINE qes_init_HubbardCommon(obj, tagname, specie, HubbardCommon, label)
    !
    IMPLICIT NONE
    !
    TYPE(HubbardCommon_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: specie
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: label
    REAL(DP), INTENT(IN) :: HubbardCommon
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%specie = specie
    IF (PRESENT(label)) THEN
      obj%label_ispresent = .TRUE.
      obj%label = label
    ELSE 
      obj%label_ispresent = .FALSE.
    END IF
    !
    obj%HubbardCommon = HubbardCommon
    !
  END SUBROUTINE qes_init_HubbardCommon 
  !
  !
  SUBROUTINE qes_init_HubbardJ(obj, tagname, specie, label, HubbardJ)
    !
    IMPLICIT NONE
    !
    TYPE(HubbardJ_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: specie
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(DP), DIMENSION(3), INTENT(IN) :: HubbardJ
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%specie = specie
    obj%label = label
    !
    obj%HubbardJ = HubbardJ
    !
  END SUBROUTINE qes_init_HubbardJ 
  !
  !
  SUBROUTINE qes_init_starting_ns(obj, tagname, specie, label, spin, starting_ns)
    !
    IMPLICIT NONE
    !
    TYPE(starting_ns_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP), DIMENSION(:), INTENT(IN) :: starting_ns
    CHARACTER(LEN=*), INTENT(IN) :: specie
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(IN) :: spin
    !
    obj%tagname = TRIM(tagname)
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%specie = specie
    obj%label = label
    obj%spin = spin
    obj%size = size(starting_ns)
    ALLOCATE(obj%starting_ns(obj%size))
    obj%starting_ns = starting_ns
    !
  END SUBROUTINE qes_init_starting_ns 
  !
  !
  SUBROUTINE qes_init_Hubbard_ns(obj, tagname, specie, label, spin, index, order, Hubbard_ns)
    !
    IMPLICIT NONE
    !
    TYPE(Hubbard_ns_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: specie
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(IN) :: spin
    INTEGER, INTENT(IN) :: index
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: Hubbard_ns
    CHARACTER(LEN=*),INTENT(IN) :: order
    INTEGER :: length, i
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%specie = specie
    obj%label = label
    obj%spin = spin
    obj%index = index
    !
     
     
    
    obj%order = order 
    
    length = 1
    obj%rank = SIZE(shape(Hubbard_ns))
    ALLOCATE ( obj%dims(obj%rank))
    obj%dims = shape(Hubbard_ns)
    DO i = 1, obj%rank
      length = length * obj%dims(i)
    END DO
    ALLOCATE(obj%Hubbard_ns(length))
    obj%Hubbard_ns(1:length) = reshape(Hubbard_ns, [length])
    !
  END SUBROUTINE qes_init_Hubbard_ns
  !
  !
  SUBROUTINE qes_init_vdW(obj, tagname, vdw_corr, dftd3_version, dftd3_threebody, non_local_term,&
                         functional, total_energy_term, london_s6, ts_vdw_econv_thr, ts_vdw_isolated,&
                         london_rcut, xdm_a1, xdm_a2, london_c6)
    !
    IMPLICIT NONE
    !
    TYPE(vdW_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: vdw_corr
    INTEGER,OPTIONAL,INTENT(IN) :: dftd3_version
    LOGICAL,OPTIONAL,INTENT(IN) :: dftd3_threebody
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: non_local_term
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: functional
    REAL(DP),OPTIONAL,INTENT(IN) :: total_energy_term
    REAL(DP),OPTIONAL,INTENT(IN) :: london_s6
    REAL(DP),OPTIONAL,INTENT(IN) :: ts_vdw_econv_thr
    LOGICAL,OPTIONAL,INTENT(IN) :: ts_vdw_isolated
    REAL(DP),OPTIONAL,INTENT(IN) :: london_rcut
    REAL(DP),OPTIONAL,INTENT(IN) :: xdm_a1
    REAL(DP),OPTIONAL,INTENT(IN) :: xdm_a2
    TYPE(HubbardCommon_type),OPTIONAL,DIMENSION(:),INTENT(IN) :: london_c6
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(vdw_corr)) THEN 
      obj%vdw_corr_ispresent = .TRUE. 
      obj%vdw_corr = vdw_corr
    ELSE 
      obj%vdw_corr_ispresent = .FALSE.
    END IF
    IF ( PRESENT(dftd3_version)) THEN 
      obj%dftd3_version_ispresent = .TRUE. 
      obj%dftd3_version = dftd3_version
    ELSE 
      obj%dftd3_version_ispresent = .FALSE.
    END IF
    IF ( PRESENT(dftd3_threebody)) THEN 
      obj%dftd3_threebody_ispresent = .TRUE. 
      obj%dftd3_threebody = dftd3_threebody
    ELSE 
      obj%dftd3_threebody_ispresent = .FALSE.
    END IF
    IF ( PRESENT(non_local_term)) THEN 
      obj%non_local_term_ispresent = .TRUE. 
      obj%non_local_term = non_local_term
    ELSE 
      obj%non_local_term_ispresent = .FALSE.
    END IF
    IF ( PRESENT(functional)) THEN 
      obj%functional_ispresent = .TRUE. 
      obj%functional = functional
    ELSE 
      obj%functional_ispresent = .FALSE.
    END IF
    IF ( PRESENT(total_energy_term)) THEN 
      obj%total_energy_term_ispresent = .TRUE. 
      obj%total_energy_term = total_energy_term
    ELSE 
      obj%total_energy_term_ispresent = .FALSE.
    END IF
    IF ( PRESENT(london_s6)) THEN 
      obj%london_s6_ispresent = .TRUE. 
      obj%london_s6 = london_s6
    ELSE 
      obj%london_s6_ispresent = .FALSE.
    END IF
    IF ( PRESENT(ts_vdw_econv_thr)) THEN 
      obj%ts_vdw_econv_thr_ispresent = .TRUE. 
      obj%ts_vdw_econv_thr = ts_vdw_econv_thr
    ELSE 
      obj%ts_vdw_econv_thr_ispresent = .FALSE.
    END IF
    IF ( PRESENT(ts_vdw_isolated)) THEN 
      obj%ts_vdw_isolated_ispresent = .TRUE. 
      obj%ts_vdw_isolated = ts_vdw_isolated
    ELSE 
      obj%ts_vdw_isolated_ispresent = .FALSE.
    END IF
    IF ( PRESENT(london_rcut)) THEN 
      obj%london_rcut_ispresent = .TRUE. 
      obj%london_rcut = london_rcut
    ELSE 
      obj%london_rcut_ispresent = .FALSE.
    END IF
    IF ( PRESENT(xdm_a1)) THEN 
      obj%xdm_a1_ispresent = .TRUE. 
      obj%xdm_a1 = xdm_a1
    ELSE 
      obj%xdm_a1_ispresent = .FALSE.
    END IF
    IF ( PRESENT(xdm_a2)) THEN 
      obj%xdm_a2_ispresent = .TRUE. 
      obj%xdm_a2 = xdm_a2
    ELSE 
      obj%xdm_a2_ispresent = .FALSE.
    END IF
    IF ( PRESENT(london_c6)) THEN 
      obj%london_c6_ispresent = .TRUE.
      ALLOCATE(obj%london_c6(SIZE(london_c6)))
      obj%ndim_london_c6 = SIZE(london_c6) 
      obj%london_c6 = london_c6
    ELSE 
      obj%london_c6_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_vdW 
  !
  !
  SUBROUTINE qes_init_spin(obj, tagname, lsda, noncolin, spinorbit)
    !
    IMPLICIT NONE
    !
    TYPE(spin_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    LOGICAL,INTENT(IN) :: lsda
    LOGICAL,INTENT(IN) :: noncolin
    LOGICAL,INTENT(IN) :: spinorbit
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%lsda = lsda
    obj%noncolin = noncolin
    obj%spinorbit = spinorbit
    !
  END SUBROUTINE qes_init_spin 
  !
  !
  SUBROUTINE qes_init_bands(obj, tagname, occupations, nbnd, smearing, tot_charge, tot_magnetization, inputOccupations)
    !
    IMPLICIT NONE
    !
    TYPE(bands_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,OPTIONAL,INTENT(IN) :: nbnd
    TYPE(smearing_type),OPTIONAL,INTENT(IN) :: smearing
    REAL(DP),OPTIONAL,INTENT(IN) :: tot_charge
    REAL(DP),OPTIONAL,INTENT(IN) :: tot_magnetization
    TYPE(occupations_type),INTENT(IN) :: occupations
    TYPE(inputOccupations_type),OPTIONAL,DIMENSION(:),INTENT(IN) :: inputOccupations
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(nbnd)) THEN 
      obj%nbnd_ispresent = .TRUE. 
      obj%nbnd = nbnd
    ELSE 
      obj%nbnd_ispresent = .FALSE.
    END IF
    IF ( PRESENT(smearing)) THEN 
      obj%smearing_ispresent = .TRUE. 
      obj%smearing = smearing
    ELSE 
      obj%smearing_ispresent = .FALSE.
    END IF
    IF ( PRESENT(tot_charge)) THEN 
      obj%tot_charge_ispresent = .TRUE. 
      obj%tot_charge = tot_charge
    ELSE 
      obj%tot_charge_ispresent = .FALSE.
    END IF
    IF ( PRESENT(tot_magnetization)) THEN 
      obj%tot_magnetization_ispresent = .TRUE. 
      obj%tot_magnetization = tot_magnetization
    ELSE 
      obj%tot_magnetization_ispresent = .FALSE.
    END IF
    obj%occupations = occupations
    IF ( PRESENT(inputOccupations)) THEN 
      obj%inputOccupations_ispresent = .TRUE.
      ALLOCATE(obj%inputOccupations(SIZE(inputOccupations)))
      obj%ndim_inputOccupations = SIZE(inputOccupations) 
      obj%inputOccupations = inputOccupations
    ELSE 
      obj%inputOccupations_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_bands 
  !
  !
  SUBROUTINE qes_init_smearing(obj, tagname, degauss, smearing)
    !
    IMPLICIT NONE
    !
    TYPE(smearing_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP), INTENT(IN) :: degauss
    CHARACTER(LEN=*), INTENT(IN) :: smearing
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%degauss = degauss
    !
    obj%smearing = smearing
    !
  END SUBROUTINE qes_init_smearing 
  !
  !
  SUBROUTINE qes_init_occupations(obj, tagname, occupations, spin)
    !
    IMPLICIT NONE
    !
    TYPE(occupations_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER, OPTIONAL, INTENT(IN) :: spin
    CHARACTER(LEN=*), INTENT(IN) :: occupations
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    IF (PRESENT(spin)) THEN
      obj%spin_ispresent = .TRUE.
      obj%spin = spin
    ELSE 
      obj%spin_ispresent = .FALSE.
    END IF
    !
    obj%occupations = occupations
    !
  END SUBROUTINE qes_init_occupations 
  !
  !
  SUBROUTINE qes_init_basis(obj, tagname, ecutwfc, gamma_only, ecutrho, fft_grid, fft_smooth, fft_box)
    !
    IMPLICIT NONE
    !
    TYPE(basis_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    LOGICAL,OPTIONAL,INTENT(IN) :: gamma_only
    REAL(DP),INTENT(IN) :: ecutwfc
    REAL(DP),OPTIONAL,INTENT(IN) :: ecutrho
    TYPE(basisSetItem_type),OPTIONAL,INTENT(IN) :: fft_grid
    TYPE(basisSetItem_type),OPTIONAL,INTENT(IN) :: fft_smooth
    TYPE(basisSetItem_type),OPTIONAL,INTENT(IN) :: fft_box
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(gamma_only)) THEN 
      obj%gamma_only_ispresent = .TRUE. 
      obj%gamma_only = gamma_only
    ELSE 
      obj%gamma_only_ispresent = .FALSE.
    END IF
    obj%ecutwfc = ecutwfc
    IF ( PRESENT(ecutrho)) THEN 
      obj%ecutrho_ispresent = .TRUE. 
      obj%ecutrho = ecutrho
    ELSE 
      obj%ecutrho_ispresent = .FALSE.
    END IF
    IF ( PRESENT(fft_grid)) THEN 
      obj%fft_grid_ispresent = .TRUE. 
      obj%fft_grid = fft_grid
    ELSE 
      obj%fft_grid_ispresent = .FALSE.
    END IF
    IF ( PRESENT(fft_smooth)) THEN 
      obj%fft_smooth_ispresent = .TRUE. 
      obj%fft_smooth = fft_smooth
    ELSE 
      obj%fft_smooth_ispresent = .FALSE.
    END IF
    IF ( PRESENT(fft_box)) THEN 
      obj%fft_box_ispresent = .TRUE. 
      obj%fft_box = fft_box
    ELSE 
      obj%fft_box_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_basis 
  !
  !
  SUBROUTINE qes_init_basis_set(obj, tagname, ecutwfc, fft_grid, ngm, npwx, reciprocal_lattice,&
                               gamma_only, ecutrho, fft_smooth, fft_box, ngms)
    !
    IMPLICIT NONE
    !
    TYPE(basis_set_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    LOGICAL,OPTIONAL,INTENT(IN) :: gamma_only
    REAL(DP),INTENT(IN) :: ecutwfc
    REAL(DP),OPTIONAL,INTENT(IN) :: ecutrho
    TYPE(basisSetItem_type),INTENT(IN) :: fft_grid
    TYPE(basisSetItem_type),OPTIONAL,INTENT(IN) :: fft_smooth
    TYPE(basisSetItem_type),OPTIONAL,INTENT(IN) :: fft_box
    INTEGER,INTENT(IN) :: ngm
    INTEGER,OPTIONAL,INTENT(IN) :: ngms
    INTEGER,INTENT(IN) :: npwx
    TYPE(reciprocal_lattice_type),INTENT(IN) :: reciprocal_lattice
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(gamma_only)) THEN 
      obj%gamma_only_ispresent = .TRUE. 
      obj%gamma_only = gamma_only
    ELSE 
      obj%gamma_only_ispresent = .FALSE.
    END IF
    obj%ecutwfc = ecutwfc
    IF ( PRESENT(ecutrho)) THEN 
      obj%ecutrho_ispresent = .TRUE. 
      obj%ecutrho = ecutrho
    ELSE 
      obj%ecutrho_ispresent = .FALSE.
    END IF
    obj%fft_grid = fft_grid
    IF ( PRESENT(fft_smooth)) THEN 
      obj%fft_smooth_ispresent = .TRUE. 
      obj%fft_smooth = fft_smooth
    ELSE 
      obj%fft_smooth_ispresent = .FALSE.
    END IF
    IF ( PRESENT(fft_box)) THEN 
      obj%fft_box_ispresent = .TRUE. 
      obj%fft_box = fft_box
    ELSE 
      obj%fft_box_ispresent = .FALSE.
    END IF
    obj%ngm = ngm
    IF ( PRESENT(ngms)) THEN 
      obj%ngms_ispresent = .TRUE. 
      obj%ngms = ngms
    ELSE 
      obj%ngms_ispresent = .FALSE.
    END IF
    obj%npwx = npwx
    obj%reciprocal_lattice = reciprocal_lattice
    !
  END SUBROUTINE qes_init_basis_set 
  !
  !
  SUBROUTINE qes_init_basisSetItem(obj, tagname, nr1, nr2, nr3, basisSetItem)
    !
    IMPLICIT NONE
    !
    TYPE(basisSetItem_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER, INTENT(IN) :: nr1
    INTEGER, INTENT(IN) :: nr2
    INTEGER, INTENT(IN) :: nr3
    CHARACTER(LEN=*), INTENT(IN) :: basisSetItem
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%nr1 = nr1
    obj%nr2 = nr2
    obj%nr3 = nr3
    !
    obj%basisSetItem = basisSetItem
    !
  END SUBROUTINE qes_init_basisSetItem 
  !
  !
  SUBROUTINE qes_init_reciprocal_lattice(obj, tagname, b1, b2, b3)
    !
    IMPLICIT NONE
    !
    TYPE(reciprocal_lattice_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP), DIMENSION(3),INTENT(IN) :: b1
    REAL(DP), DIMENSION(3),INTENT(IN) :: b2
    REAL(DP), DIMENSION(3),INTENT(IN) :: b3
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%b1 = b1
    obj%b2 = b2
    obj%b3 = b3
    !
  END SUBROUTINE qes_init_reciprocal_lattice 
  !
  !
  SUBROUTINE qes_init_electron_control(obj, tagname, diagonalization, mixing_mode, mixing_beta,&
                                      conv_thr, mixing_ndim, max_nstep, tq_smoothing, tbeta_smoothing,&
                                      diago_thr_init, diago_full_acc, real_space_q, real_space_beta,&
                                      diago_cg_maxiter, diago_ppcg_maxiter, diago_david_ndim)
    !
    IMPLICIT NONE
    !
    TYPE(electron_control_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),INTENT(IN) :: diagonalization
    CHARACTER(LEN=*),INTENT(IN) :: mixing_mode
    REAL(DP),INTENT(IN) :: mixing_beta
    REAL(DP),INTENT(IN) :: conv_thr
    INTEGER,INTENT(IN) :: mixing_ndim
    INTEGER,INTENT(IN) :: max_nstep
    LOGICAL,OPTIONAL,INTENT(IN) :: real_space_q
    LOGICAL,OPTIONAL,INTENT(IN) :: real_space_beta
    LOGICAL,INTENT(IN) :: tq_smoothing
    LOGICAL,INTENT(IN) :: tbeta_smoothing
    REAL(DP),INTENT(IN) :: diago_thr_init
    LOGICAL,INTENT(IN) :: diago_full_acc
    INTEGER,OPTIONAL,INTENT(IN) :: diago_cg_maxiter
    INTEGER,OPTIONAL,INTENT(IN) :: diago_ppcg_maxiter
    INTEGER,OPTIONAL,INTENT(IN) :: diago_david_ndim
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%diagonalization = diagonalization
    obj%mixing_mode = mixing_mode
    obj%mixing_beta = mixing_beta
    obj%conv_thr = conv_thr
    obj%mixing_ndim = mixing_ndim
    obj%max_nstep = max_nstep
    IF ( PRESENT(real_space_q)) THEN 
      obj%real_space_q_ispresent = .TRUE. 
      obj%real_space_q = real_space_q
    ELSE 
      obj%real_space_q_ispresent = .FALSE.
    END IF
    IF ( PRESENT(real_space_beta)) THEN 
      obj%real_space_beta_ispresent = .TRUE. 
      obj%real_space_beta = real_space_beta
    ELSE 
      obj%real_space_beta_ispresent = .FALSE.
    END IF
    obj%tq_smoothing = tq_smoothing
    obj%tbeta_smoothing = tbeta_smoothing
    obj%diago_thr_init = diago_thr_init
    obj%diago_full_acc = diago_full_acc
    IF ( PRESENT(diago_cg_maxiter)) THEN 
      obj%diago_cg_maxiter_ispresent = .TRUE. 
      obj%diago_cg_maxiter = diago_cg_maxiter
    ELSE 
      obj%diago_cg_maxiter_ispresent = .FALSE.
    END IF
    IF ( PRESENT(diago_ppcg_maxiter)) THEN 
      obj%diago_ppcg_maxiter_ispresent = .TRUE. 
      obj%diago_ppcg_maxiter = diago_ppcg_maxiter
    ELSE 
      obj%diago_ppcg_maxiter_ispresent = .FALSE.
    END IF
    IF ( PRESENT(diago_david_ndim)) THEN 
      obj%diago_david_ndim_ispresent = .TRUE. 
      obj%diago_david_ndim = diago_david_ndim
    ELSE 
      obj%diago_david_ndim_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_electron_control 
  !
  !
  SUBROUTINE qes_init_k_points_IBZ(obj, tagname, monkhorst_pack, nk, k_point)
    !
    IMPLICIT NONE
    !
    TYPE(k_points_IBZ_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(monkhorst_pack_type),OPTIONAL,INTENT(IN) :: monkhorst_pack
    INTEGER,OPTIONAL,INTENT(IN) :: nk
    TYPE(k_point_type),OPTIONAL,DIMENSION(:),INTENT(IN) :: k_point
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(monkhorst_pack)) THEN 
      obj%monkhorst_pack_ispresent = .TRUE. 
      obj%monkhorst_pack = monkhorst_pack
    ELSE 
      obj%monkhorst_pack_ispresent = .FALSE.
    END IF
    IF ( PRESENT(nk)) THEN 
      obj%nk_ispresent = .TRUE. 
      obj%nk = nk
    ELSE 
      obj%nk_ispresent = .FALSE.
    END IF
    IF ( PRESENT(k_point)) THEN 
      obj%k_point_ispresent = .TRUE.
      ALLOCATE(obj%k_point(SIZE(k_point)))
      obj%ndim_k_point = SIZE(k_point) 
      obj%k_point = k_point
    ELSE 
      obj%k_point_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_k_points_IBZ 
  !
  !
  SUBROUTINE qes_init_monkhorst_pack(obj, tagname, nk1, nk2, nk3, k1, k2, k3, monkhorst_pack)
    !
    IMPLICIT NONE
    !
    TYPE(monkhorst_pack_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER, INTENT(IN) :: nk1
    INTEGER, INTENT(IN) :: nk2
    INTEGER, INTENT(IN) :: nk3
    INTEGER, INTENT(IN) :: k1
    INTEGER, INTENT(IN) :: k2
    INTEGER, INTENT(IN) :: k3
    CHARACTER(LEN=*), INTENT(IN) :: monkhorst_pack
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%nk1 = nk1
    obj%nk2 = nk2
    obj%nk3 = nk3
    obj%k1 = k1
    obj%k2 = k2
    obj%k3 = k3
    !
    obj%monkhorst_pack = monkhorst_pack
    !
  END SUBROUTINE qes_init_monkhorst_pack 
  !
  !
  SUBROUTINE qes_init_k_point(obj, tagname, k_point, weight, label)
    !
    IMPLICIT NONE
    !
    TYPE(k_point_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP), OPTIONAL, INTENT(IN) :: weight
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: label
    REAL(DP), DIMENSION(3), INTENT(IN) :: k_point
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    IF (PRESENT(weight)) THEN
      obj%weight_ispresent = .TRUE.
      obj%weight = weight
    ELSE 
      obj%weight_ispresent = .FALSE.
    END IF
    IF (PRESENT(label)) THEN
      obj%label_ispresent = .TRUE.
      obj%label = label
    ELSE 
      obj%label_ispresent = .FALSE.
    END IF
    !
    obj%k_point = k_point
    !
  END SUBROUTINE qes_init_k_point 
  !
  !
  SUBROUTINE qes_init_ion_control(obj, tagname, ion_dynamics, upscale, remove_rigid_rot, refold_pos,&
                                 bfgs, md)
    !
    IMPLICIT NONE
    !
    TYPE(ion_control_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),INTENT(IN) :: ion_dynamics
    REAL(DP),OPTIONAL,INTENT(IN) :: upscale
    LOGICAL,OPTIONAL,INTENT(IN) :: remove_rigid_rot
    LOGICAL,OPTIONAL,INTENT(IN) :: refold_pos
    TYPE(bfgs_type),OPTIONAL,INTENT(IN) :: bfgs
    TYPE(md_type),OPTIONAL,INTENT(IN) :: md
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%ion_dynamics = ion_dynamics
    IF ( PRESENT(upscale)) THEN 
      obj%upscale_ispresent = .TRUE. 
      obj%upscale = upscale
    ELSE 
      obj%upscale_ispresent = .FALSE.
    END IF
    IF ( PRESENT(remove_rigid_rot)) THEN 
      obj%remove_rigid_rot_ispresent = .TRUE. 
      obj%remove_rigid_rot = remove_rigid_rot
    ELSE 
      obj%remove_rigid_rot_ispresent = .FALSE.
    END IF
    IF ( PRESENT(refold_pos)) THEN 
      obj%refold_pos_ispresent = .TRUE. 
      obj%refold_pos = refold_pos
    ELSE 
      obj%refold_pos_ispresent = .FALSE.
    END IF
    IF ( PRESENT(bfgs)) THEN 
      obj%bfgs_ispresent = .TRUE. 
      obj%bfgs = bfgs
    ELSE 
      obj%bfgs_ispresent = .FALSE.
    END IF
    IF ( PRESENT(md)) THEN 
      obj%md_ispresent = .TRUE. 
      obj%md = md
    ELSE 
      obj%md_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_ion_control 
  !
  !
  SUBROUTINE qes_init_bfgs(obj, tagname, ndim, trust_radius_min, trust_radius_max, trust_radius_init,&
                          w1, w2)
    !
    IMPLICIT NONE
    !
    TYPE(bfgs_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,INTENT(IN) :: ndim
    REAL(DP),INTENT(IN) :: trust_radius_min
    REAL(DP),INTENT(IN) :: trust_radius_max
    REAL(DP),INTENT(IN) :: trust_radius_init
    REAL(DP),INTENT(IN) :: w1
    REAL(DP),INTENT(IN) :: w2
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%ndim = ndim
    obj%trust_radius_min = trust_radius_min
    obj%trust_radius_max = trust_radius_max
    obj%trust_radius_init = trust_radius_init
    obj%w1 = w1
    obj%w2 = w2
    !
  END SUBROUTINE qes_init_bfgs 
  !
  !
  SUBROUTINE qes_init_md(obj, tagname, pot_extrapolation, wfc_extrapolation, ion_temperature,&
                        timestep, tempw, tolp, deltaT, nraise)
    !
    IMPLICIT NONE
    !
    TYPE(md_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),INTENT(IN) :: pot_extrapolation
    CHARACTER(LEN=*),INTENT(IN) :: wfc_extrapolation
    CHARACTER(LEN=*),INTENT(IN) :: ion_temperature
    REAL(DP),INTENT(IN) :: timestep
    REAL(DP),INTENT(IN) :: tempw
    REAL(DP),INTENT(IN) :: tolp
    REAL(DP),INTENT(IN) :: deltaT
    INTEGER,INTENT(IN) :: nraise
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%pot_extrapolation = pot_extrapolation
    obj%wfc_extrapolation = wfc_extrapolation
    obj%ion_temperature = ion_temperature
    obj%timestep = timestep
    obj%tempw = tempw
    obj%tolp = tolp
    obj%deltaT = deltaT
    obj%nraise = nraise
    !
  END SUBROUTINE qes_init_md 
  !
  !
  SUBROUTINE qes_init_cell_control(obj, tagname, cell_dynamics, pressure, wmass, cell_factor,&
                                  fix_volume, fix_area, isotropic, free_cell)
    !
    IMPLICIT NONE
    !
    TYPE(cell_control_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),INTENT(IN) :: cell_dynamics
    REAL(DP),INTENT(IN) :: pressure
    REAL(DP),OPTIONAL,INTENT(IN) :: wmass
    REAL(DP),OPTIONAL,INTENT(IN) :: cell_factor
    LOGICAL,OPTIONAL,INTENT(IN) :: fix_volume
    LOGICAL,OPTIONAL,INTENT(IN) :: fix_area
    LOGICAL,OPTIONAL,INTENT(IN) :: isotropic
    TYPE(integerMatrix_type),OPTIONAL,INTENT(IN) :: free_cell
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%cell_dynamics = cell_dynamics
    obj%pressure = pressure
    IF ( PRESENT(wmass)) THEN 
      obj%wmass_ispresent = .TRUE. 
      obj%wmass = wmass
    ELSE 
      obj%wmass_ispresent = .FALSE.
    END IF
    IF ( PRESENT(cell_factor)) THEN 
      obj%cell_factor_ispresent = .TRUE. 
      obj%cell_factor = cell_factor
    ELSE 
      obj%cell_factor_ispresent = .FALSE.
    END IF
    IF ( PRESENT(fix_volume)) THEN 
      obj%fix_volume_ispresent = .TRUE. 
      obj%fix_volume = fix_volume
    ELSE 
      obj%fix_volume_ispresent = .FALSE.
    END IF
    IF ( PRESENT(fix_area)) THEN 
      obj%fix_area_ispresent = .TRUE. 
      obj%fix_area = fix_area
    ELSE 
      obj%fix_area_ispresent = .FALSE.
    END IF
    IF ( PRESENT(isotropic)) THEN 
      obj%isotropic_ispresent = .TRUE. 
      obj%isotropic = isotropic
    ELSE 
      obj%isotropic_ispresent = .FALSE.
    END IF
    IF ( PRESENT(free_cell)) THEN 
      obj%free_cell_ispresent = .TRUE. 
      obj%free_cell = free_cell
    ELSE 
      obj%free_cell_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_cell_control 
  !
  !
  SUBROUTINE qes_init_symmetry_flags(obj, tagname, nosym, nosym_evc, noinv, no_t_rev, force_symmorphic, use_all_frac)
    !
    IMPLICIT NONE
    !
    TYPE(symmetry_flags_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    LOGICAL,INTENT(IN) :: nosym
    LOGICAL,INTENT(IN) :: nosym_evc
    LOGICAL,INTENT(IN) :: noinv
    LOGICAL,INTENT(IN) :: no_t_rev
    LOGICAL,INTENT(IN) :: force_symmorphic
    LOGICAL,INTENT(IN) :: use_all_frac
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%nosym = nosym
    obj%nosym_evc = nosym_evc
    obj%noinv = noinv
    obj%no_t_rev = no_t_rev
    obj%force_symmorphic = force_symmorphic
    obj%use_all_frac = use_all_frac
    !
  END SUBROUTINE qes_init_symmetry_flags 
  !
  !
  SUBROUTINE qes_init_boundary_conditions(obj, tagname, assume_isolated, esm, fcp_opt, fcp_mu)
    !
    IMPLICIT NONE
    !
    TYPE(boundary_conditions_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),INTENT(IN) :: assume_isolated
    TYPE(esm_type),OPTIONAL,INTENT(IN) :: esm
    LOGICAL,OPTIONAL,INTENT(IN) :: fcp_opt
    REAL(DP),OPTIONAL,INTENT(IN) :: fcp_mu
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%assume_isolated = assume_isolated
    IF ( PRESENT(esm)) THEN 
      obj%esm_ispresent = .TRUE. 
      obj%esm = esm
    ELSE 
      obj%esm_ispresent = .FALSE.
    END IF
    IF ( PRESENT(fcp_opt)) THEN 
      obj%fcp_opt_ispresent = .TRUE. 
      obj%fcp_opt = fcp_opt
    ELSE 
      obj%fcp_opt_ispresent = .FALSE.
    END IF
    IF ( PRESENT(fcp_mu)) THEN 
      obj%fcp_mu_ispresent = .TRUE. 
      obj%fcp_mu = fcp_mu
    ELSE 
      obj%fcp_mu_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_boundary_conditions 
  !
  !
  SUBROUTINE qes_init_esm(obj, tagname, bc, nfit, w, efield)
    !
    IMPLICIT NONE
    !
    TYPE(esm_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),INTENT(IN) :: bc
    INTEGER,INTENT(IN) :: nfit
    REAL(DP),INTENT(IN) :: w
    REAL(DP),INTENT(IN) :: efield
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%bc = bc
    obj%nfit = nfit
    obj%w = w
    obj%efield = efield
    !
  END SUBROUTINE qes_init_esm 
  !
  !
  SUBROUTINE qes_init_ekin_functional(obj, tagname, ecfixed, qcutz, q2sigma)
    !
    IMPLICIT NONE
    !
    TYPE(ekin_functional_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP),INTENT(IN) :: ecfixed
    REAL(DP),INTENT(IN) :: qcutz
    REAL(DP),INTENT(IN) :: q2sigma
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%ecfixed = ecfixed
    obj%qcutz = qcutz
    obj%q2sigma = q2sigma
    !
  END SUBROUTINE qes_init_ekin_functional 
  !
  !
  SUBROUTINE qes_init_spin_constraints(obj, tagname, spin_constraints, lagrange_multiplier, target_magnetization)
    !
    IMPLICIT NONE
    !
    TYPE(spin_constraints_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),INTENT(IN) :: spin_constraints
    REAL(DP),INTENT(IN) :: lagrange_multiplier
    REAL(DP), DIMENSION(3),OPTIONAL,INTENT(IN) :: target_magnetization
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%spin_constraints = spin_constraints
    obj%lagrange_multiplier = lagrange_multiplier
    IF ( PRESENT(target_magnetization)) THEN 
      obj%target_magnetization_ispresent = .TRUE. 
      obj%target_magnetization = target_magnetization
    ELSE 
      obj%target_magnetization_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_spin_constraints 
  !
  !
  SUBROUTINE qes_init_electric_field(obj, tagname, electric_potential, dipole_correction, gate_settings,&
                                    electric_field_direction, potential_max_position, potential_decrease_width,&
                                    electric_field_amplitude, electric_field_vector, nk_per_string, n_berry_cycles &
                                    )
    !
    IMPLICIT NONE
    !
    TYPE(electric_field_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),INTENT(IN) :: electric_potential
    LOGICAL,OPTIONAL,INTENT(IN) :: dipole_correction
    TYPE(gate_settings_type),OPTIONAL,INTENT(IN) :: gate_settings
    INTEGER,OPTIONAL,INTENT(IN) :: electric_field_direction
    REAL(DP),OPTIONAL,INTENT(IN) :: potential_max_position
    REAL(DP),OPTIONAL,INTENT(IN) :: potential_decrease_width
    REAL(DP),OPTIONAL,INTENT(IN) :: electric_field_amplitude
    REAL(DP), DIMENSION(3),OPTIONAL,INTENT(IN) :: electric_field_vector
    INTEGER,OPTIONAL,INTENT(IN) :: nk_per_string
    INTEGER,OPTIONAL,INTENT(IN) :: n_berry_cycles
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%electric_potential = electric_potential
    IF ( PRESENT(dipole_correction)) THEN 
      obj%dipole_correction_ispresent = .TRUE. 
      obj%dipole_correction = dipole_correction
    ELSE 
      obj%dipole_correction_ispresent = .FALSE.
    END IF
    IF ( PRESENT(gate_settings)) THEN 
      obj%gate_settings_ispresent = .TRUE. 
      obj%gate_settings = gate_settings
    ELSE 
      obj%gate_settings_ispresent = .FALSE.
    END IF
    IF ( PRESENT(electric_field_direction)) THEN 
      obj%electric_field_direction_ispresent = .TRUE. 
      obj%electric_field_direction = electric_field_direction
    ELSE 
      obj%electric_field_direction_ispresent = .FALSE.
    END IF
    IF ( PRESENT(potential_max_position)) THEN 
      obj%potential_max_position_ispresent = .TRUE. 
      obj%potential_max_position = potential_max_position
    ELSE 
      obj%potential_max_position_ispresent = .FALSE.
    END IF
    IF ( PRESENT(potential_decrease_width)) THEN 
      obj%potential_decrease_width_ispresent = .TRUE. 
      obj%potential_decrease_width = potential_decrease_width
    ELSE 
      obj%potential_decrease_width_ispresent = .FALSE.
    END IF
    IF ( PRESENT(electric_field_amplitude)) THEN 
      obj%electric_field_amplitude_ispresent = .TRUE. 
      obj%electric_field_amplitude = electric_field_amplitude
    ELSE 
      obj%electric_field_amplitude_ispresent = .FALSE.
    END IF
    IF ( PRESENT(electric_field_vector)) THEN 
      obj%electric_field_vector_ispresent = .TRUE. 
      obj%electric_field_vector = electric_field_vector
    ELSE 
      obj%electric_field_vector_ispresent = .FALSE.
    END IF
    IF ( PRESENT(nk_per_string)) THEN 
      obj%nk_per_string_ispresent = .TRUE. 
      obj%nk_per_string = nk_per_string
    ELSE 
      obj%nk_per_string_ispresent = .FALSE.
    END IF
    IF ( PRESENT(n_berry_cycles)) THEN 
      obj%n_berry_cycles_ispresent = .TRUE. 
      obj%n_berry_cycles = n_berry_cycles
    ELSE 
      obj%n_berry_cycles_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_electric_field 
  !
  !
  SUBROUTINE qes_init_gate_settings(obj, tagname, use_gate, zgate, relaxz, block, block_1, block_2, block_height)
    !
    IMPLICIT NONE
    !
    TYPE(gate_settings_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    LOGICAL,INTENT(IN) :: use_gate
    REAL(DP),OPTIONAL,INTENT(IN) :: zgate
    LOGICAL,OPTIONAL,INTENT(IN) :: relaxz
    LOGICAL,OPTIONAL,INTENT(IN) :: block
    REAL(DP),OPTIONAL,INTENT(IN) :: block_1
    REAL(DP),OPTIONAL,INTENT(IN) :: block_2
    REAL(DP),OPTIONAL,INTENT(IN) :: block_height
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%use_gate = use_gate
    IF ( PRESENT(zgate)) THEN 
      obj%zgate_ispresent = .TRUE. 
      obj%zgate = zgate
    ELSE 
      obj%zgate_ispresent = .FALSE.
    END IF
    IF ( PRESENT(relaxz)) THEN 
      obj%relaxz_ispresent = .TRUE. 
      obj%relaxz = relaxz
    ELSE 
      obj%relaxz_ispresent = .FALSE.
    END IF
    IF ( PRESENT(block)) THEN 
      obj%block_ispresent = .TRUE. 
      obj%block = block
    ELSE 
      obj%block_ispresent = .FALSE.
    END IF
    IF ( PRESENT(block_1)) THEN 
      obj%block_1_ispresent = .TRUE. 
      obj%block_1 = block_1
    ELSE 
      obj%block_1_ispresent = .FALSE.
    END IF
    IF ( PRESENT(block_2)) THEN 
      obj%block_2_ispresent = .TRUE. 
      obj%block_2 = block_2
    ELSE 
      obj%block_2_ispresent = .FALSE.
    END IF
    IF ( PRESENT(block_height)) THEN 
      obj%block_height_ispresent = .TRUE. 
      obj%block_height = block_height
    ELSE 
      obj%block_height_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_gate_settings 
  !
  !
  SUBROUTINE qes_init_atomic_constraints(obj, tagname, num_of_constraints, tolerance, atomic_constraint)
    !
    IMPLICIT NONE
    !
    TYPE(atomic_constraints_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,INTENT(IN) :: num_of_constraints
    REAL(DP),INTENT(IN) :: tolerance
    TYPE(atomic_constraint_type),DIMENSION(:),INTENT(IN) :: atomic_constraint
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%num_of_constraints = num_of_constraints
    obj%tolerance = tolerance
    ALLOCATE( obj%atomic_constraint(SIZE(atomic_constraint))) 
    obj%ndim_atomic_constraint = SIZE(atomic_constraint)
    obj%atomic_constraint = atomic_constraint
    !
  END SUBROUTINE qes_init_atomic_constraints 
  !
  !
  SUBROUTINE qes_init_atomic_constraint(obj, tagname, constr_parms, constr_type, constr_target)
    !
    IMPLICIT NONE
    !
    TYPE(atomic_constraint_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP), DIMENSION(4),INTENT(IN) :: constr_parms
    CHARACTER(LEN=*),INTENT(IN) :: constr_type
    REAL(DP),INTENT(IN) :: constr_target
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%constr_parms = constr_parms
    obj%constr_type = constr_type
    obj%constr_target = constr_target
    !
  END SUBROUTINE qes_init_atomic_constraint 
  !
  !
  SUBROUTINE qes_init_inputOccupations(obj, tagname, ispin, spin_factor, inputOccupations)
    !
    IMPLICIT NONE
    !
    TYPE(inputOccupations_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP), DIMENSION(:), INTENT(IN) :: inputOccupations
    INTEGER, INTENT(IN) :: ispin
    REAL(DP), INTENT(IN) :: spin_factor
    !
    obj%tagname = TRIM(tagname)
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%ispin = ispin
    obj%spin_factor = spin_factor
    obj%size = size(inputOccupations)
    ALLOCATE(obj%inputOccupations(obj%size))
    obj%inputOccupations = inputOccupations
    !
  END SUBROUTINE qes_init_inputOccupations 
  !
  !
  SUBROUTINE qes_init_outputElectricField(obj, tagname, BerryPhase, finiteElectricFieldInfo, dipoleInfo, gateInfo)
    !
    IMPLICIT NONE
    !
    TYPE(outputElectricField_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(BerryPhaseOutput_type),OPTIONAL,INTENT(IN) :: BerryPhase
    TYPE(finiteFieldOut_type),OPTIONAL,INTENT(IN) :: finiteElectricFieldInfo
    TYPE(dipoleOutput_type),OPTIONAL,INTENT(IN) :: dipoleInfo
    TYPE(gateInfo_type),OPTIONAL,INTENT(IN) :: gateInfo
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    IF ( PRESENT(BerryPhase)) THEN 
      obj%BerryPhase_ispresent = .TRUE. 
      obj%BerryPhase = BerryPhase
    ELSE 
      obj%BerryPhase_ispresent = .FALSE.
    END IF
    IF ( PRESENT(finiteElectricFieldInfo)) THEN 
      obj%finiteElectricFieldInfo_ispresent = .TRUE. 
      obj%finiteElectricFieldInfo = finiteElectricFieldInfo
    ELSE 
      obj%finiteElectricFieldInfo_ispresent = .FALSE.
    END IF
    IF ( PRESENT(dipoleInfo)) THEN 
      obj%dipoleInfo_ispresent = .TRUE. 
      obj%dipoleInfo = dipoleInfo
    ELSE 
      obj%dipoleInfo_ispresent = .FALSE.
    END IF
    IF ( PRESENT(gateInfo)) THEN 
      obj%gateInfo_ispresent = .TRUE. 
      obj%gateInfo = gateInfo
    ELSE 
      obj%gateInfo_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_outputElectricField 
  !
  !
  SUBROUTINE qes_init_BerryPhaseOutput(obj, tagname, totalPolarization, totalPhase, ionicPolarization, electronicPolarization)
    !
    IMPLICIT NONE
    !
    TYPE(BerryPhaseOutput_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(polarization_type),INTENT(IN) :: totalPolarization
    TYPE(phase_type),INTENT(IN) :: totalPhase
    TYPE(ionicPolarization_type),DIMENSION(:),INTENT(IN) :: ionicPolarization
    TYPE(electronicPolarization_type),DIMENSION(:),INTENT(IN) :: electronicPolarization
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%totalPolarization = totalPolarization
    obj%totalPhase = totalPhase
    ALLOCATE( obj%ionicPolarization(SIZE(ionicPolarization))) 
    obj%ndim_ionicPolarization = SIZE(ionicPolarization)
    obj%ionicPolarization = ionicPolarization
    ALLOCATE( obj%electronicPolarization(SIZE(electronicPolarization))) 
    obj%ndim_electronicPolarization = SIZE(electronicPolarization)
    obj%electronicPolarization = electronicPolarization
    !
  END SUBROUTINE qes_init_BerryPhaseOutput 
  !
  !
  SUBROUTINE qes_init_dipoleOutput(obj, tagname, idir, dipole, ion_dipole, elec_dipole, dipoleField,&
                                  potentialAmp, totalLength)
    !
    IMPLICIT NONE
    !
    TYPE(dipoleOutput_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,INTENT(IN) :: idir
    TYPE(scalarQuantity_type),INTENT(IN) :: dipole
    TYPE(scalarQuantity_type),INTENT(IN) :: ion_dipole
    TYPE(scalarQuantity_type),INTENT(IN) :: elec_dipole
    TYPE(scalarQuantity_type),INTENT(IN) :: dipoleField
    TYPE(scalarQuantity_type),INTENT(IN) :: potentialAmp
    TYPE(scalarQuantity_type),INTENT(IN) :: totalLength
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%idir = idir
    obj%dipole = dipole
    obj%ion_dipole = ion_dipole
    obj%elec_dipole = elec_dipole
    obj%dipoleField = dipoleField
    obj%potentialAmp = potentialAmp
    obj%totalLength = totalLength
    !
  END SUBROUTINE qes_init_dipoleOutput 
  !
  !
  SUBROUTINE qes_init_finiteFieldOut(obj, tagname, electronicDipole, ionicDipole)
    !
    IMPLICIT NONE
    !
    TYPE(finiteFieldOut_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP), DIMENSION(3),INTENT(IN) :: electronicDipole
    REAL(DP), DIMENSION(3),INTENT(IN) :: ionicDipole
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%electronicDipole = electronicDipole
    obj%ionicDipole = ionicDipole
    !
  END SUBROUTINE qes_init_finiteFieldOut 
  !
  !
  SUBROUTINE qes_init_polarization(obj, tagname, polarization, modulus, direction)
    !
    IMPLICIT NONE
    !
    TYPE(polarization_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(scalarQuantity_type),INTENT(IN) :: polarization
    REAL(DP),INTENT(IN) :: modulus
    REAL(DP), DIMENSION(3),INTENT(IN) :: direction
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%polarization = polarization
    obj%modulus = modulus
    obj%direction = direction
    !
  END SUBROUTINE qes_init_polarization 
  !
  !
  SUBROUTINE qes_init_ionicPolarization(obj, tagname, ion, charge, phase)
    !
    IMPLICIT NONE
    !
    TYPE(ionicPolarization_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(atom_type),INTENT(IN) :: ion
    REAL(DP),INTENT(IN) :: charge
    TYPE(phase_type),INTENT(IN) :: phase
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%ion = ion
    obj%charge = charge
    obj%phase = phase
    !
  END SUBROUTINE qes_init_ionicPolarization 
  !
  !
  SUBROUTINE qes_init_electronicPolarization(obj, tagname, firstKeyPoint, phase, spin)
    !
    IMPLICIT NONE
    !
    TYPE(electronicPolarization_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(k_point_type),INTENT(IN) :: firstKeyPoint
    INTEGER,OPTIONAL,INTENT(IN) :: spin
    TYPE(phase_type),INTENT(IN) :: phase
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%firstKeyPoint = firstKeyPoint
    IF ( PRESENT(spin)) THEN 
      obj%spin_ispresent = .TRUE. 
      obj%spin = spin
    ELSE 
      obj%spin_ispresent = .FALSE.
    END IF
    obj%phase = phase
    !
  END SUBROUTINE qes_init_electronicPolarization 
  !
  !
  SUBROUTINE qes_init_phase(obj, tagname, phase, ionic, electronic, modulus)
    !
    IMPLICIT NONE
    !
    TYPE(phase_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP), OPTIONAL, INTENT(IN) :: ionic
    REAL(DP), OPTIONAL, INTENT(IN) :: electronic
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: modulus
    REAL(DP), INTENT(IN) :: phase
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    IF (PRESENT(ionic)) THEN
      obj%ionic_ispresent = .TRUE.
      obj%ionic = ionic
    ELSE 
      obj%ionic_ispresent = .FALSE.
    END IF
    IF (PRESENT(electronic)) THEN
      obj%electronic_ispresent = .TRUE.
      obj%electronic = electronic
    ELSE 
      obj%electronic_ispresent = .FALSE.
    END IF
    IF (PRESENT(modulus)) THEN
      obj%modulus_ispresent = .TRUE.
      obj%modulus = modulus
    ELSE 
      obj%modulus_ispresent = .FALSE.
    END IF
    !
    obj%phase = phase
    !
  END SUBROUTINE qes_init_phase 
  !
  !
  SUBROUTINE qes_init_gateInfo(obj, tagname, pot_prefactor, gate_zpos, gate_gate_term, gatefieldEnergy)
    !
    IMPLICIT NONE
    !
    TYPE(gateInfo_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP),INTENT(IN) :: pot_prefactor
    REAL(DP),INTENT(IN) :: gate_zpos
    REAL(DP),INTENT(IN) :: gate_gate_term
    REAL(DP),INTENT(IN) :: gatefieldEnergy
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%pot_prefactor = pot_prefactor
    obj%gate_zpos = gate_zpos
    obj%gate_gate_term = gate_gate_term
    obj%gatefieldEnergy = gatefieldEnergy
    !
  END SUBROUTINE qes_init_gateInfo 
  !
  !
  SUBROUTINE qes_init_convergence_info(obj, tagname, scf_conv, opt_conv)
    !
    IMPLICIT NONE
    !
    TYPE(convergence_info_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(scf_conv_type),INTENT(IN) :: scf_conv
    TYPE(opt_conv_type),OPTIONAL,INTENT(IN) :: opt_conv
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%scf_conv = scf_conv
    IF ( PRESENT(opt_conv)) THEN 
      obj%opt_conv_ispresent = .TRUE. 
      obj%opt_conv = opt_conv
    ELSE 
      obj%opt_conv_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_convergence_info 
  !
  !
  SUBROUTINE qes_init_scf_conv(obj, tagname, convergence_achieved, n_scf_steps, scf_error)
    !
    IMPLICIT NONE
    !
    TYPE(scf_conv_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    LOGICAL,INTENT(IN) :: convergence_achieved
    INTEGER,INTENT(IN) :: n_scf_steps
    REAL(DP),INTENT(IN) :: scf_error
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%convergence_achieved = convergence_achieved
    obj%n_scf_steps = n_scf_steps
    obj%scf_error = scf_error
    !
  END SUBROUTINE qes_init_scf_conv 
  !
  !
  SUBROUTINE qes_init_opt_conv(obj, tagname, convergence_achieved, n_opt_steps, grad_norm)
    !
    IMPLICIT NONE
    !
    TYPE(opt_conv_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    LOGICAL,INTENT(IN) :: convergence_achieved
    INTEGER,INTENT(IN) :: n_opt_steps
    REAL(DP),INTENT(IN) :: grad_norm
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%convergence_achieved = convergence_achieved
    obj%n_opt_steps = n_opt_steps
    obj%grad_norm = grad_norm
    !
  END SUBROUTINE qes_init_opt_conv 
  !
  !
  SUBROUTINE qes_init_algorithmic_info(obj, tagname, real_space_q, uspp, paw, real_space_beta)
    !
    IMPLICIT NONE
    !
    TYPE(algorithmic_info_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    LOGICAL,INTENT(IN) :: real_space_q
    LOGICAL,OPTIONAL,INTENT(IN) :: real_space_beta
    LOGICAL,INTENT(IN) :: uspp
    LOGICAL,INTENT(IN) :: paw
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%real_space_q = real_space_q
    IF ( PRESENT(real_space_beta)) THEN 
      obj%real_space_beta_ispresent = .TRUE. 
      obj%real_space_beta = real_space_beta
    ELSE 
      obj%real_space_beta_ispresent = .FALSE.
    END IF
    obj%uspp = uspp
    obj%paw = paw
    !
  END SUBROUTINE qes_init_algorithmic_info 
  !
  !
  SUBROUTINE qes_init_symmetries(obj, tagname, nsym, nrot, space_group, symmetry)
    !
    IMPLICIT NONE
    !
    TYPE(symmetries_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,INTENT(IN) :: nsym
    INTEGER,INTENT(IN) :: nrot
    INTEGER,INTENT(IN) :: space_group
    TYPE(symmetry_type),DIMENSION(:),INTENT(IN) :: symmetry
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%nsym = nsym
    obj%nrot = nrot
    obj%space_group = space_group
    ALLOCATE( obj%symmetry(SIZE(symmetry))) 
    obj%ndim_symmetry = SIZE(symmetry)
    obj%symmetry = symmetry
    !
  END SUBROUTINE qes_init_symmetries 
  !
  !
  SUBROUTINE qes_init_symmetry(obj, tagname, info, rotation, fractional_translation, equivalent_atoms)
    !
    IMPLICIT NONE
    !
    TYPE(symmetry_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(info_type),INTENT(IN) :: info
    TYPE(matrix_type),INTENT(IN) :: rotation
    REAL(DP), DIMENSION(3),OPTIONAL,INTENT(IN) :: fractional_translation
    TYPE(equivalent_atoms_type),OPTIONAL,INTENT(IN) :: equivalent_atoms
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%info = info
    obj%rotation = rotation
    IF ( PRESENT(fractional_translation)) THEN 
      obj%fractional_translation_ispresent = .TRUE. 
      obj%fractional_translation = fractional_translation
    ELSE 
      obj%fractional_translation_ispresent = .FALSE.
    END IF
    IF ( PRESENT(equivalent_atoms)) THEN 
      obj%equivalent_atoms_ispresent = .TRUE. 
      obj%equivalent_atoms = equivalent_atoms
    ELSE 
      obj%equivalent_atoms_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_symmetry 
  !
  !
  SUBROUTINE qes_init_equivalent_atoms(obj, tagname, nat, equivalent_atoms)
    !
    IMPLICIT NONE
    !
    TYPE(equivalent_atoms_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER, DIMENSION(:), INTENT(IN) :: equivalent_atoms
    INTEGER, INTENT(IN) :: nat
    !
    obj%tagname = TRIM(tagname)
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%nat = nat
    obj%size = size(equivalent_atoms)
    ALLOCATE(obj%equivalent_atoms(obj%size))
    obj%equivalent_atoms = equivalent_atoms
    !
  END SUBROUTINE qes_init_equivalent_atoms 
  !
  !
  SUBROUTINE qes_init_info(obj, tagname, info, name, class, time_reversal)
    !
    IMPLICIT NONE
    !
    TYPE(info_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: name
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: class
    LOGICAL, OPTIONAL, INTENT(IN) :: time_reversal
    CHARACTER(LEN=*), INTENT(IN) :: info
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    IF (PRESENT(name)) THEN
      obj%name_ispresent = .TRUE.
      obj%name = name
    ELSE 
      obj%name_ispresent = .FALSE.
    END IF
    IF (PRESENT(class)) THEN
      obj%class_ispresent = .TRUE.
      obj%class = class
    ELSE 
      obj%class_ispresent = .FALSE.
    END IF
    IF (PRESENT(time_reversal)) THEN
      obj%time_reversal_ispresent = .TRUE.
      obj%time_reversal = time_reversal
    ELSE 
      obj%time_reversal_ispresent = .FALSE.
    END IF
    !
    obj%info = info
    !
  END SUBROUTINE qes_init_info 
  !
  !
  SUBROUTINE qes_init_outputPBC(obj, tagname, assume_isolated)
    !
    IMPLICIT NONE
    !
    TYPE(outputPBC_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*),INTENT(IN) :: assume_isolated
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%assume_isolated = assume_isolated
    !
  END SUBROUTINE qes_init_outputPBC 
  !
  !
  SUBROUTINE qes_init_magnetization(obj, tagname, lsda, noncolin, spinorbit, total, absolute, do_magnetization)
    !
    IMPLICIT NONE
    !
    TYPE(magnetization_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    LOGICAL,INTENT(IN) :: lsda
    LOGICAL,INTENT(IN) :: noncolin
    LOGICAL,INTENT(IN) :: spinorbit
    REAL(DP),INTENT(IN) :: total
    REAL(DP),INTENT(IN) :: absolute
    LOGICAL,INTENT(IN) :: do_magnetization
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%lsda = lsda
    obj%noncolin = noncolin
    obj%spinorbit = spinorbit
    obj%total = total
    obj%absolute = absolute
    obj%do_magnetization = do_magnetization
    !
  END SUBROUTINE qes_init_magnetization 
  !
  !
  SUBROUTINE qes_init_total_energy(obj, tagname, etot, eband, ehart, vtxc, etxc, ewald, demet,&
                                  efieldcorr, potentiostat_contr, gatefield_contr, vdW_term)
    !
    IMPLICIT NONE
    !
    TYPE(total_energy_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP),INTENT(IN) :: etot
    REAL(DP),OPTIONAL,INTENT(IN) :: eband
    REAL(DP),OPTIONAL,INTENT(IN) :: ehart
    REAL(DP),OPTIONAL,INTENT(IN) :: vtxc
    REAL(DP),OPTIONAL,INTENT(IN) :: etxc
    REAL(DP),OPTIONAL,INTENT(IN) :: ewald
    REAL(DP),OPTIONAL,INTENT(IN) :: demet
    REAL(DP),OPTIONAL,INTENT(IN) :: efieldcorr
    REAL(DP),OPTIONAL,INTENT(IN) :: potentiostat_contr
    REAL(DP),OPTIONAL,INTENT(IN) :: gatefield_contr
    REAL(DP),OPTIONAL,INTENT(IN) :: vdW_term
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%etot = etot
    IF ( PRESENT(eband)) THEN 
      obj%eband_ispresent = .TRUE. 
      obj%eband = eband
    ELSE 
      obj%eband_ispresent = .FALSE.
    END IF
    IF ( PRESENT(ehart)) THEN 
      obj%ehart_ispresent = .TRUE. 
      obj%ehart = ehart
    ELSE 
      obj%ehart_ispresent = .FALSE.
    END IF
    IF ( PRESENT(vtxc)) THEN 
      obj%vtxc_ispresent = .TRUE. 
      obj%vtxc = vtxc
    ELSE 
      obj%vtxc_ispresent = .FALSE.
    END IF
    IF ( PRESENT(etxc)) THEN 
      obj%etxc_ispresent = .TRUE. 
      obj%etxc = etxc
    ELSE 
      obj%etxc_ispresent = .FALSE.
    END IF
    IF ( PRESENT(ewald)) THEN 
      obj%ewald_ispresent = .TRUE. 
      obj%ewald = ewald
    ELSE 
      obj%ewald_ispresent = .FALSE.
    END IF
    IF ( PRESENT(demet)) THEN 
      obj%demet_ispresent = .TRUE. 
      obj%demet = demet
    ELSE 
      obj%demet_ispresent = .FALSE.
    END IF
    IF ( PRESENT(efieldcorr)) THEN 
      obj%efieldcorr_ispresent = .TRUE. 
      obj%efieldcorr = efieldcorr
    ELSE 
      obj%efieldcorr_ispresent = .FALSE.
    END IF
    IF ( PRESENT(potentiostat_contr)) THEN 
      obj%potentiostat_contr_ispresent = .TRUE. 
      obj%potentiostat_contr = potentiostat_contr
    ELSE 
      obj%potentiostat_contr_ispresent = .FALSE.
    END IF
    IF ( PRESENT(gatefield_contr)) THEN 
      obj%gatefield_contr_ispresent = .TRUE. 
      obj%gatefield_contr = gatefield_contr
    ELSE 
      obj%gatefield_contr_ispresent = .FALSE.
    END IF
    IF ( PRESENT(vdW_term)) THEN 
      obj%vdW_term_ispresent = .TRUE. 
      obj%vdW_term = vdW_term
    ELSE 
      obj%vdW_term_ispresent = .FALSE.
    END IF
    !
  END SUBROUTINE qes_init_total_energy 
  !
  !
  SUBROUTINE qes_init_band_structure(obj, tagname, lsda, noncolin, spinorbit, nelec, wf_collected,&
                                    starting_k_points, nks, occupations_kind, ks_energies, nbnd,&
                                    nbnd_up, nbnd_dw, num_of_atomic_wfc, fermi_energy, highestOccupiedLevel,&
                                    lowestUnoccupiedLevel, two_fermi_energies, smearing)
    !
    IMPLICIT NONE
    !
    TYPE(band_structure_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    LOGICAL,INTENT(IN) :: lsda
    LOGICAL,INTENT(IN) :: noncolin
    LOGICAL,INTENT(IN) :: spinorbit
    INTEGER,OPTIONAL,INTENT(IN) :: nbnd
    INTEGER,OPTIONAL,INTENT(IN) :: nbnd_up
    INTEGER,OPTIONAL,INTENT(IN) :: nbnd_dw
    REAL(DP),INTENT(IN) :: nelec
    INTEGER,OPTIONAL,INTENT(IN) :: num_of_atomic_wfc
    LOGICAL,INTENT(IN) :: wf_collected
    REAL(DP),OPTIONAL,INTENT(IN) :: fermi_energy
    REAL(DP),OPTIONAL,INTENT(IN) :: highestOccupiedLevel
    REAL(DP),OPTIONAL,INTENT(IN) :: lowestUnoccupiedLevel
    REAL(DP), DIMENSION(2),OPTIONAL,INTENT(IN) :: two_fermi_energies
    TYPE(k_points_IBZ_type),INTENT(IN) :: starting_k_points
    INTEGER,INTENT(IN) :: nks
    TYPE(occupations_type),INTENT(IN) :: occupations_kind
    TYPE(smearing_type),OPTIONAL,INTENT(IN) :: smearing
    TYPE(ks_energies_type),DIMENSION(:),INTENT(IN) :: ks_energies
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%lsda = lsda
    obj%noncolin = noncolin
    obj%spinorbit = spinorbit
    IF ( PRESENT(nbnd)) THEN 
      obj%nbnd_ispresent = .TRUE. 
      obj%nbnd = nbnd
    ELSE 
      obj%nbnd_ispresent = .FALSE.
    END IF
    IF ( PRESENT(nbnd_up)) THEN 
      obj%nbnd_up_ispresent = .TRUE. 
      obj%nbnd_up = nbnd_up
    ELSE 
      obj%nbnd_up_ispresent = .FALSE.
    END IF
    IF ( PRESENT(nbnd_dw)) THEN 
      obj%nbnd_dw_ispresent = .TRUE. 
      obj%nbnd_dw = nbnd_dw
    ELSE 
      obj%nbnd_dw_ispresent = .FALSE.
    END IF
    obj%nelec = nelec
    IF ( PRESENT(num_of_atomic_wfc)) THEN 
      obj%num_of_atomic_wfc_ispresent = .TRUE. 
      obj%num_of_atomic_wfc = num_of_atomic_wfc
    ELSE 
      obj%num_of_atomic_wfc_ispresent = .FALSE.
    END IF
    obj%wf_collected = wf_collected
    IF ( PRESENT(fermi_energy)) THEN 
      obj%fermi_energy_ispresent = .TRUE. 
      obj%fermi_energy = fermi_energy
    ELSE 
      obj%fermi_energy_ispresent = .FALSE.
    END IF
    IF ( PRESENT(highestOccupiedLevel)) THEN 
      obj%highestOccupiedLevel_ispresent = .TRUE. 
      obj%highestOccupiedLevel = highestOccupiedLevel
    ELSE 
      obj%highestOccupiedLevel_ispresent = .FALSE.
    END IF
    IF ( PRESENT(lowestUnoccupiedLevel)) THEN 
      obj%lowestUnoccupiedLevel_ispresent = .TRUE. 
      obj%lowestUnoccupiedLevel = lowestUnoccupiedLevel
    ELSE 
      obj%lowestUnoccupiedLevel_ispresent = .FALSE.
    END IF
    IF ( PRESENT(two_fermi_energies)) THEN 
      obj%two_fermi_energies_ispresent = .TRUE. 
      obj%two_fermi_energies = two_fermi_energies
    ELSE 
      obj%two_fermi_energies_ispresent = .FALSE.
    END IF
    obj%starting_k_points = starting_k_points
    obj%nks = nks
    obj%occupations_kind = occupations_kind
    IF ( PRESENT(smearing)) THEN 
      obj%smearing_ispresent = .TRUE. 
      obj%smearing = smearing
    ELSE 
      obj%smearing_ispresent = .FALSE.
    END IF
    ALLOCATE( obj%ks_energies(SIZE(ks_energies))) 
    obj%ndim_ks_energies = SIZE(ks_energies)
    obj%ks_energies = ks_energies
    !
  END SUBROUTINE qes_init_band_structure 
  !
  !
  SUBROUTINE qes_init_ks_energies(obj, tagname, k_point, npw, eigenvalues, occupations)
    !
    IMPLICIT NONE
    !
    TYPE(ks_energies_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    TYPE(k_point_type),INTENT(IN) :: k_point
    INTEGER,INTENT(IN) :: npw
    TYPE(vector_type),INTENT(IN) :: eigenvalues
    TYPE(vector_type),INTENT(IN) :: occupations
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%k_point = k_point
    obj%npw = npw
    obj%eigenvalues = eigenvalues
    obj%occupations = occupations
    !
  END SUBROUTINE qes_init_ks_energies 
  !
  !
  SUBROUTINE qes_init_closed(obj, tagname, DATE, TIME, closed)
    !
    IMPLICIT NONE
    !
    TYPE(closed_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: DATE
    CHARACTER(LEN=*), INTENT(IN) :: TIME
    CHARACTER(LEN=*), INTENT(IN) :: closed
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%DATE = DATE
    obj%TIME = TIME
    !
    obj%closed = closed
    !
  END SUBROUTINE qes_init_closed 
  !
  !
  SUBROUTINE qes_init_vector(obj, tagname, vector)
    !
    IMPLICIT NONE
    !
    TYPE(vector_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    REAL(DP), DIMENSION(:), INTENT(IN) :: vector
    !
    obj%tagname = TRIM(tagname)
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%size = size(vector)
    ALLOCATE(obj%vector(obj%size))
    obj%vector = vector
    !
  END SUBROUTINE qes_init_vector 
  !
  !
  SUBROUTINE qes_init_integerVector(obj, tagname, integerVector)
    !
    IMPLICIT NONE
    !
    TYPE(integerVector_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER, DIMENSION(:), INTENT(IN) :: integerVector
    !
    obj%tagname = TRIM(tagname)
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    !
    obj%size = size(integerVector)
    ALLOCATE(obj%integerVector(obj%size))
    obj%integerVector = integerVector
    !
  END SUBROUTINE qes_init_integerVector 
  !

  !
  SUBROUTINE qes_init_matrix_1(obj, tagname, dims, mat, order)
    !
    IMPLICIT NONE
    !
    TYPE(matrix_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,DIMENSION(:),INTENT(IN) :: dims
    REAL(DP), INTENT(IN) :: mat(:)
    CHARACTER(LEN=*),OPTIONAL :: order
    INTEGER :: rank, length, i
    !
    obj%tagname = TRIM(tagname)
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    length = 1
    rank = SIZE(dims)
    DO i = 1, rank
      length = length * dims(i)
    END DO
    obj%rank = rank
    ALLOCATE(obj%matrix(length), obj%dims(rank) )
    obj%matrix(1:length) = mat(1:length)
    obj%dims = dims
    IF (PRESENT(order)) THEN
      obj%order = TRIM(order)
    ELSE
      obj%order = 'F'
    END IF
    !
  END SUBROUTINE qes_init_matrix_1
  !
  !
  SUBROUTINE qes_init_matrix_2(obj, tagname, dims, mat, order)
    !
    IMPLICIT NONE
    !
    TYPE(matrix_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,DIMENSION(:),INTENT(IN) :: dims
    REAL(DP), INTENT(IN) :: mat(:,:)
    CHARACTER(LEN=*),OPTIONAL :: order
    INTEGER :: rank, length, i
    !
    obj%tagname = TRIM(tagname)
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    length = 1
    rank = SIZE(dims)
    DO i = 1, rank
      length = length * dims(i)
    END DO
    obj%rank = rank
    ALLOCATE(obj%matrix(length), obj%dims(rank) )
    obj%matrix(1:length) = reshape(mat, [length])
    obj%dims = dims
    IF (PRESENT(order)) THEN
      obj%order = TRIM(order)
    ELSE
      obj%order = 'F'
    END IF
    !
  END SUBROUTINE qes_init_matrix_2
  !
  !
  SUBROUTINE qes_init_matrix_3(obj, tagname, dims, mat, order)
    !
    IMPLICIT NONE
    !
    TYPE(matrix_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,DIMENSION(:),INTENT(IN) :: dims
    REAL(DP), INTENT(IN) :: mat(:,:,:)
    CHARACTER(LEN=*),OPTIONAL :: order
    INTEGER :: rank, length, i
    !
    obj%tagname = TRIM(tagname)
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    length = 1
    rank = SIZE(dims)
    DO i = 1, rank
      length = length * dims(i)
    END DO
    obj%rank = rank
    ALLOCATE(obj%matrix(length), obj%dims(rank) )
    obj%matrix(1:length) = reshape(mat, [length])
    obj%dims = dims
    IF (PRESENT(order)) THEN
      obj%order = TRIM(order)
    ELSE
      obj%order = 'F'
    END IF
    !
  END SUBROUTINE qes_init_matrix_3
  !

  !
  SUBROUTINE qes_init_integerMatrix_1(obj, tagname, dims, mat, order)
    !
    IMPLICIT NONE
    !
    TYPE(integerMatrix_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,DIMENSION(:),INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: mat(:)
    CHARACTER(LEN=*),OPTIONAL :: order
    INTEGER :: rank, length, i
    !
    obj%tagname = TRIM(tagname)
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    length = 1
    rank = SIZE(dims)
    DO i = 1, rank
      length = length * dims(i)
    END DO
    obj%rank = rank
    ALLOCATE(obj%integerMatrix(length), obj%dims(rank) )
    obj%integerMatrix(1:length) = mat(1:length)
    obj%dims = dims
    IF (PRESENT(order)) THEN
      obj%order = TRIM(order)
    ELSE
      obj%order = 'F'
    END IF
    !
  END SUBROUTINE qes_init_integerMatrix_1
  !
  !
  SUBROUTINE qes_init_integerMatrix_2(obj, tagname, dims, mat, order)
    !
    IMPLICIT NONE
    !
    TYPE(integerMatrix_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,DIMENSION(:),INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: mat(:,:)
    CHARACTER(LEN=*),OPTIONAL :: order
    INTEGER :: rank, length, i
    !
    obj%tagname = TRIM(tagname)
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    length = 1
    rank = SIZE(dims)
    DO i = 1, rank
      length = length * dims(i)
    END DO
    obj%rank = rank
    ALLOCATE(obj%integerMatrix(length), obj%dims(rank) )
    obj%integerMatrix(1:length) = reshape(mat, [length])
    obj%dims = dims
    IF (PRESENT(order)) THEN
      obj%order = TRIM(order)
    ELSE
      obj%order = 'F'
    END IF
    !
  END SUBROUTINE qes_init_integerMatrix_2
  !
  !
  SUBROUTINE qes_init_integerMatrix_3(obj, tagname, dims, mat, order)
    !
    IMPLICIT NONE
    !
    TYPE(integerMatrix_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    INTEGER,DIMENSION(:),INTENT(IN) :: dims
    INTEGER, INTENT(IN) :: mat(:,:,:)
    CHARACTER(LEN=*),OPTIONAL :: order
    INTEGER :: rank, length, i
    !
    obj%tagname = TRIM(tagname)
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    length = 1
    rank = SIZE(dims)
    DO i = 1, rank
      length = length * dims(i)
    END DO
    obj%rank = rank
    ALLOCATE(obj%integerMatrix(length), obj%dims(rank) )
    obj%integerMatrix(1:length) = reshape(mat, [length])
    obj%dims = dims
    IF (PRESENT(order)) THEN
      obj%order = TRIM(order)
    ELSE
      obj%order = 'F'
    END IF
    !
  END SUBROUTINE qes_init_integerMatrix_3
  !
  !
  SUBROUTINE qes_init_scalarQuantity(obj, tagname, Units, scalarQuantity)
    !
    IMPLICIT NONE
    !
    TYPE(scalarQuantity_type), INTENT(OUT) :: obj
    CHARACTER(LEN=*), INTENT(IN) :: tagname
    CHARACTER(LEN=*), INTENT(IN) :: Units
    REAL(DP), INTENT(IN) :: scalarQuantity
    !
    obj%tagname = TRIM(tagname) 
    obj%lwrite = .TRUE.
    obj%lread = .TRUE.
    obj%Units = Units
    !
    obj%scalarQuantity = scalarQuantity
    !
  END SUBROUTINE qes_init_scalarQuantity 
  !
  !
END MODULE qes_init_module
