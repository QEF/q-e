!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE qes_read_module
  !
  ! Auto-generated code: don't edit or at least don't commit changes
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes-1.0 
  !
  USE FoX_dom
  USE qes_types_module
  !
  USE kinds, only: DP
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: qes_read
  !
  INTERFACE qes_read
    MODULE PROCEDURE qes_read_general_info
    MODULE PROCEDURE qes_read_parallel_info
    MODULE PROCEDURE qes_read_input
    MODULE PROCEDURE qes_read_step
    MODULE PROCEDURE qes_read_output
    MODULE PROCEDURE qes_read_control_variables
    MODULE PROCEDURE qes_read_xml_format
    MODULE PROCEDURE qes_read_creator
    MODULE PROCEDURE qes_read_created
    MODULE PROCEDURE qes_read_atomic_species
    MODULE PROCEDURE qes_read_species
    MODULE PROCEDURE qes_read_atomic_structure
    MODULE PROCEDURE qes_read_atomic_positions
    MODULE PROCEDURE qes_read_atom
    MODULE PROCEDURE qes_read_wyckoff_positions
    MODULE PROCEDURE qes_read_cell
    MODULE PROCEDURE qes_read_dft
    MODULE PROCEDURE qes_read_hybrid
    MODULE PROCEDURE qes_read_qpoint_grid
    MODULE PROCEDURE qes_read_dftU
    MODULE PROCEDURE qes_read_HubbardCommon
    MODULE PROCEDURE qes_read_HubbardJ
    MODULE PROCEDURE qes_read_starting_ns
    MODULE PROCEDURE qes_read_Hubbard_ns
    MODULE PROCEDURE qes_read_vdW
    MODULE PROCEDURE qes_read_spin
    MODULE PROCEDURE qes_read_bands
    MODULE PROCEDURE qes_read_smearing
    MODULE PROCEDURE qes_read_occupations
    MODULE PROCEDURE qes_read_basis
    MODULE PROCEDURE qes_read_basis_set
    MODULE PROCEDURE qes_read_basisSetItem
    MODULE PROCEDURE qes_read_reciprocal_lattice
    MODULE PROCEDURE qes_read_electron_control
    MODULE PROCEDURE qes_read_k_points_IBZ
    MODULE PROCEDURE qes_read_monkhorst_pack
    MODULE PROCEDURE qes_read_k_point
    MODULE PROCEDURE qes_read_ion_control
    MODULE PROCEDURE qes_read_bfgs
    MODULE PROCEDURE qes_read_md
    MODULE PROCEDURE qes_read_cell_control
    MODULE PROCEDURE qes_read_symmetry_flags
    MODULE PROCEDURE qes_read_boundary_conditions
    MODULE PROCEDURE qes_read_esm
    MODULE PROCEDURE qes_read_ekin_functional
    MODULE PROCEDURE qes_read_spin_constraints
    MODULE PROCEDURE qes_read_electric_field
    MODULE PROCEDURE qes_read_gate_settings
    MODULE PROCEDURE qes_read_atomic_constraints
    MODULE PROCEDURE qes_read_atomic_constraint
    MODULE PROCEDURE qes_read_inputOccupations
    MODULE PROCEDURE qes_read_outputElectricField
    MODULE PROCEDURE qes_read_BerryPhaseOutput
    MODULE PROCEDURE qes_read_dipoleOutput
    MODULE PROCEDURE qes_read_finiteFieldOut
    MODULE PROCEDURE qes_read_polarization
    MODULE PROCEDURE qes_read_ionicPolarization
    MODULE PROCEDURE qes_read_electronicPolarization
    MODULE PROCEDURE qes_read_phase
    MODULE PROCEDURE qes_read_gateInfo
    MODULE PROCEDURE qes_read_convergence_info
    MODULE PROCEDURE qes_read_scf_conv
    MODULE PROCEDURE qes_read_opt_conv
    MODULE PROCEDURE qes_read_algorithmic_info
    MODULE PROCEDURE qes_read_symmetries
    MODULE PROCEDURE qes_read_symmetry
    MODULE PROCEDURE qes_read_equivalent_atoms
    MODULE PROCEDURE qes_read_info
    MODULE PROCEDURE qes_read_outputPBC
    MODULE PROCEDURE qes_read_magnetization
    MODULE PROCEDURE qes_read_total_energy
    MODULE PROCEDURE qes_read_band_structure
    MODULE PROCEDURE qes_read_ks_energies
    MODULE PROCEDURE qes_read_closed
    MODULE PROCEDURE qes_read_vector
    MODULE PROCEDURE qes_read_integerVector
    MODULE PROCEDURE qes_read_matrix
    MODULE PROCEDURE qes_read_integerMatrix
    MODULE PROCEDURE qes_read_scalarQuantity
  END INTERFACE qes_read
  !
  CONTAINS
  !
  !
  SUBROUTINE qes_read_general_info(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(general_info_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "xml_format")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:general_infoType","xml_format: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:general_infoType","xml_format: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_xml_format(tmp_node, obj%xml_format, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "creator")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:general_infoType","creator: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:general_infoType","creator: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_creator(tmp_node, obj%creator, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "created")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:general_infoType","created: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:general_infoType","created: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_created(tmp_node, obj%created, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "job")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:general_infoType","job: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:general_infoType","job: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%job, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:general_infoType","error reading job")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:general_infoType","error reading job",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_general_info
  !
  !
  SUBROUTINE qes_read_parallel_info(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(parallel_info_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "nprocs")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:parallel_infoType","nprocs: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:parallel_infoType","nprocs: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nprocs, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:parallel_infoType","error reading nprocs")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:parallel_infoType","error reading nprocs",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nthreads")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:parallel_infoType","nthreads: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:parallel_infoType","nthreads: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nthreads, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:parallel_infoType","error reading nthreads")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:parallel_infoType","error reading nthreads",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ntasks")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:parallel_infoType","ntasks: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:parallel_infoType","ntasks: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ntasks, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:parallel_infoType","error reading ntasks")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:parallel_infoType","error reading ntasks",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nbgrp")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:parallel_infoType","nbgrp: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:parallel_infoType","nbgrp: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nbgrp, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:parallel_infoType","error reading nbgrp")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:parallel_infoType","error reading nbgrp",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "npool")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:parallel_infoType","npool: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:parallel_infoType","npool: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%npool, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:parallel_infoType","error reading npool")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:parallel_infoType","error reading npool",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ndiag")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:parallel_infoType","ndiag: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:parallel_infoType","ndiag: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ndiag, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:parallel_infoType","error reading ndiag")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:parallel_infoType","error reading ndiag",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_parallel_info
  !
  !
  SUBROUTINE qes_read_input(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(input_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "control_variables")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","control_variables: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","control_variables: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_control_variables(tmp_node, obj%control_variables, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "atomic_species")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","atomic_species: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","atomic_species: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_atomic_species(tmp_node, obj%atomic_species, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "atomic_structure")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","atomic_structure: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","atomic_structure: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_atomic_structure(tmp_node, obj%atomic_structure, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "dft")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","dft: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","dft: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_dft(tmp_node, obj%dft, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "spin")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","spin: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","spin: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_spin(tmp_node, obj%spin, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "bands")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","bands: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","bands: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_bands(tmp_node, obj%bands, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "basis")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","basis: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","basis: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_basis(tmp_node, obj%basis, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "electron_control")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","electron_control: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","electron_control: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_electron_control(tmp_node, obj%electron_control, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "k_points_IBZ")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","k_points_IBZ: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","k_points_IBZ: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_k_points_IBZ(tmp_node, obj%k_points_IBZ, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "ion_control")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","ion_control: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","ion_control: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_ion_control(tmp_node, obj%ion_control, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "cell_control")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","cell_control: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","cell_control: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_cell_control(tmp_node, obj%cell_control, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "symmetry_flags")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","symmetry_flags: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","symmetry_flags: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%symmetry_flags_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_symmetry_flags(tmp_node, obj%symmetry_flags, ierr )
    ELSE
       obj%symmetry_flags_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "boundary_conditions")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","boundary_conditions: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","boundary_conditions: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%boundary_conditions_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_boundary_conditions(tmp_node, obj%boundary_conditions, ierr )
    ELSE
       obj%boundary_conditions_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ekin_functional")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","ekin_functional: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","ekin_functional: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%ekin_functional_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_ekin_functional(tmp_node, obj%ekin_functional, ierr )
    ELSE
       obj%ekin_functional_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "external_atomic_forces")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","external_atomic_forces: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","external_atomic_forces: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%external_atomic_forces_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_matrix(tmp_node, obj%external_atomic_forces, ierr )
    ELSE
       obj%external_atomic_forces_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "free_positions")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","free_positions: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","free_positions: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%free_positions_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_integerMatrix(tmp_node, obj%free_positions, ierr )
    ELSE
       obj%free_positions_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "starting_atomic_velocities")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","starting_atomic_velocities: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","starting_atomic_velocities: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%starting_atomic_velocities_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_matrix(tmp_node, obj%starting_atomic_velocities, ierr )
    ELSE
       obj%starting_atomic_velocities_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "electric_field")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","electric_field: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","electric_field: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%electric_field_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_electric_field(tmp_node, obj%electric_field, ierr )
    ELSE
       obj%electric_field_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "atomic_constraints")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","atomic_constraints: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","atomic_constraints: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%atomic_constraints_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_atomic_constraints(tmp_node, obj%atomic_constraints, ierr )
    ELSE
       obj%atomic_constraints_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "spin_constraints")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:inputType","spin_constraints: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:inputType","spin_constraints: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%spin_constraints_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_spin_constraints(tmp_node, obj%spin_constraints, ierr )
    ELSE
       obj%spin_constraints_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_input
  !
  !
  SUBROUTINE qes_read_step(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(step_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "n_step")) THEN
      CALL extractDataAttribute(xml_node, "n_step", obj%n_step)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: stepType",&
                        "required attribute n_step not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: stepType",&
                      "required attribute n_step not found", 10 ) 
      END IF
    END IF
    !
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "scf_conv")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:stepType","scf_conv: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:stepType","scf_conv: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_scf_conv(tmp_node, obj%scf_conv, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "atomic_structure")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:stepType","atomic_structure: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:stepType","atomic_structure: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_atomic_structure(tmp_node, obj%atomic_structure, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "total_energy")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:stepType","total_energy: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:stepType","total_energy: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_total_energy(tmp_node, obj%total_energy, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "forces")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:stepType","forces: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:stepType","forces: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_matrix(tmp_node, obj%forces, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "stress")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:stepType","stress: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:stepType","stress: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%stress_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_matrix(tmp_node, obj%stress, ierr )
    ELSE
       obj%stress_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "FCP_force")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:stepType","FCP_force: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:stepType","FCP_force: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%FCP_force_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%FCP_force , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:stepType","error reading FCP_force")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:stepType","error reading FCP_force",10)
         END IF
      END IF
    ELSE
       obj%FCP_force_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "FCP_tot_charge")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:stepType","FCP_tot_charge: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:stepType","FCP_tot_charge: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%FCP_tot_charge_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%FCP_tot_charge , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:stepType","error reading FCP_tot_charge")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:stepType","error reading FCP_tot_charge",10)
         END IF
      END IF
    ELSE
       obj%FCP_tot_charge_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_step
  !
  !
  SUBROUTINE qes_read_output(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(output_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "convergence_info")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","convergence_info: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","convergence_info: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_convergence_info(tmp_node, obj%convergence_info, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "algorithmic_info")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","algorithmic_info: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","algorithmic_info: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_algorithmic_info(tmp_node, obj%algorithmic_info, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "atomic_species")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","atomic_species: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","atomic_species: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_atomic_species(tmp_node, obj%atomic_species, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "atomic_structure")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","atomic_structure: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","atomic_structure: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_atomic_structure(tmp_node, obj%atomic_structure, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "symmetries")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","symmetries: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","symmetries: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%symmetries_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_symmetries(tmp_node, obj%symmetries, ierr )
    ELSE
       obj%symmetries_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "basis_set")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","basis_set: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","basis_set: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_basis_set(tmp_node, obj%basis_set, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "dft")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","dft: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","dft: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_dft(tmp_node, obj%dft, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "boundary_conditions")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","boundary_conditions: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","boundary_conditions: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%boundary_conditions_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_outputPBC(tmp_node, obj%boundary_conditions, ierr )
    ELSE
       obj%boundary_conditions_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "magnetization")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","magnetization: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","magnetization: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_magnetization(tmp_node, obj%magnetization, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "total_energy")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","total_energy: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","total_energy: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_total_energy(tmp_node, obj%total_energy, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "band_structure")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","band_structure: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","band_structure: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_band_structure(tmp_node, obj%band_structure, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "forces")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","forces: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","forces: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%forces_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_matrix(tmp_node, obj%forces, ierr )
    ELSE
       obj%forces_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "stress")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","stress: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","stress: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%stress_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_matrix(tmp_node, obj%stress, ierr )
    ELSE
       obj%stress_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "electric_field")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","electric_field: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","electric_field: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%electric_field_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_outputElectricField(tmp_node, obj%electric_field, ierr )
    ELSE
       obj%electric_field_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "FCP_force")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","FCP_force: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","FCP_force: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%FCP_force_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%FCP_force , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:outputType","error reading FCP_force")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:outputType","error reading FCP_force",10)
         END IF
      END IF
    ELSE
       obj%FCP_force_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "FCP_tot_charge")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputType","FCP_tot_charge: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputType","FCP_tot_charge: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%FCP_tot_charge_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%FCP_tot_charge , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:outputType","error reading FCP_tot_charge")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:outputType","error reading FCP_tot_charge",10)
         END IF
      END IF
    ELSE
       obj%FCP_tot_charge_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_output
  !
  !
  SUBROUTINE qes_read_control_variables(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(control_variables_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "title")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","title: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","title: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%title, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading title")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading title",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "calculation")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","calculation: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","calculation: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%calculation, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading calculation")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading calculation",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "restart_mode")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","restart_mode: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","restart_mode: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%restart_mode, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading restart_mode")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading restart_mode",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "prefix")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","prefix: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","prefix: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%prefix, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading prefix")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading prefix",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "pseudo_dir")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","pseudo_dir: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","pseudo_dir: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%pseudo_dir, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading pseudo_dir")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading pseudo_dir",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "outdir")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","outdir: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","outdir: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%outdir, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading outdir")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading outdir",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "stress")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","stress: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","stress: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%stress, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading stress")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading stress",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "forces")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","forces: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","forces: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%forces, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading forces")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading forces",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "wf_collect")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","wf_collect: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","wf_collect: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%wf_collect, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading wf_collect")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading wf_collect",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "disk_io")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","disk_io: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","disk_io: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%disk_io, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading disk_io")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading disk_io",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "max_seconds")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","max_seconds: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","max_seconds: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%max_seconds, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading max_seconds")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading max_seconds",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nstep")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","nstep: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","nstep: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%nstep_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%nstep , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:control_variablesType","error reading nstep")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:control_variablesType","error reading nstep",10)
         END IF
      END IF
    ELSE
       obj%nstep_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "etot_conv_thr")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","etot_conv_thr: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","etot_conv_thr: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%etot_conv_thr, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading etot_conv_thr")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading etot_conv_thr",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "forc_conv_thr")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","forc_conv_thr: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","forc_conv_thr: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%forc_conv_thr, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading forc_conv_thr")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading forc_conv_thr",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "press_conv_thr")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","press_conv_thr: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","press_conv_thr: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%press_conv_thr, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading press_conv_thr")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading press_conv_thr",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "verbosity")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","verbosity: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","verbosity: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%verbosity, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading verbosity")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading verbosity",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "print_every")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:control_variablesType","print_every: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:control_variablesType","print_every: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%print_every, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:control_variablesType","error reading print_every")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:control_variablesType","error reading print_every",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_control_variables
  !
  !
  SUBROUTINE qes_read_xml_format(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(xml_format_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "NAME")) THEN
      CALL extractDataAttribute(xml_node, "NAME", obj%NAME)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: xml_formatType",&
                        "required attribute NAME not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: xml_formatType",&
                      "required attribute NAME not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "VERSION")) THEN
      CALL extractDataAttribute(xml_node, "VERSION", obj%VERSION)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: xml_formatType",&
                        "required attribute VERSION not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: xml_formatType",&
                      "required attribute VERSION not found", 10 ) 
      END IF
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%xml_format )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_xml_format
  !
  !
  SUBROUTINE qes_read_creator(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(creator_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "NAME")) THEN
      CALL extractDataAttribute(xml_node, "NAME", obj%NAME)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: creatorType",&
                        "required attribute NAME not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: creatorType",&
                      "required attribute NAME not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "VERSION")) THEN
      CALL extractDataAttribute(xml_node, "VERSION", obj%VERSION)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: creatorType",&
                        "required attribute VERSION not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: creatorType",&
                      "required attribute VERSION not found", 10 ) 
      END IF
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%creator )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_creator
  !
  !
  SUBROUTINE qes_read_created(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(created_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "DATE")) THEN
      CALL extractDataAttribute(xml_node, "DATE", obj%DATE)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: createdType",&
                        "required attribute DATE not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: createdType",&
                      "required attribute DATE not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "TIME")) THEN
      CALL extractDataAttribute(xml_node, "TIME", obj%TIME)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: createdType",&
                        "required attribute TIME not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: createdType",&
                      "required attribute TIME not found", 10 ) 
      END IF
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%created )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_created
  !
  !
  SUBROUTINE qes_read_atomic_species(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(atomic_species_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "ntyp")) THEN
      CALL extractDataAttribute(xml_node, "ntyp", obj%ntyp)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: atomic_speciesType",&
                        "required attribute ntyp not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: atomic_speciesType",&
                      "required attribute ntyp not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "pseudo_dir")) THEN
      CALL extractDataAttribute(xml_node, "pseudo_dir", obj%pseudo_dir)
      obj%pseudo_dir_ispresent = .TRUE.
    ELSE
      obj%pseudo_dir_ispresent = .FALSE.
    END IF
    !
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "species")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size < 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:atomic_speciesType","species: not enough elements")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:atomic_speciesType","species: not enough elements",10)
        END IF
    END IF
    !
    obj%ndim_species = tmp_node_list_size
    ALLOCATE(obj%species(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_species(tmp_node, obj%species(index), ierr )
    END DO
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_atomic_species
  !
  !
  SUBROUTINE qes_read_species(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(species_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "name")) THEN
      CALL extractDataAttribute(xml_node, "name", obj%name)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: speciesType",&
                        "required attribute name not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: speciesType",&
                      "required attribute name not found", 10 ) 
      END IF
    END IF
    !
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "mass")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:speciesType","mass: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:speciesType","mass: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%mass_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%mass , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:speciesType","error reading mass")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:speciesType","error reading mass",10)
         END IF
      END IF
    ELSE
       obj%mass_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "pseudo_file")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:speciesType","pseudo_file: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:speciesType","pseudo_file: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%pseudo_file, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:speciesType","error reading pseudo_file")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:speciesType","error reading pseudo_file",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "starting_magnetization")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:speciesType","starting_magnetization: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:speciesType","starting_magnetization: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%starting_magnetization_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%starting_magnetization , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:speciesType","error reading starting_magnetization")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:speciesType","error reading starting_magnetization",10)
         END IF
      END IF
    ELSE
       obj%starting_magnetization_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "spin_teta")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:speciesType","spin_teta: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:speciesType","spin_teta: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%spin_teta_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%spin_teta , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:speciesType","error reading spin_teta")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:speciesType","error reading spin_teta",10)
         END IF
      END IF
    ELSE
       obj%spin_teta_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "spin_phi")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:speciesType","spin_phi: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:speciesType","spin_phi: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%spin_phi_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%spin_phi , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:speciesType","error reading spin_phi")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:speciesType","error reading spin_phi",10)
         END IF
      END IF
    ELSE
       obj%spin_phi_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_species
  !
  !
  SUBROUTINE qes_read_atomic_structure(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(atomic_structure_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "nat")) THEN
      CALL extractDataAttribute(xml_node, "nat", obj%nat)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: atomic_structureType",&
                        "required attribute nat not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: atomic_structureType",&
                      "required attribute nat not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "alat")) THEN
      CALL extractDataAttribute(xml_node, "alat", obj%alat)
      obj%alat_ispresent = .TRUE.
    ELSE
      obj%alat_ispresent = .FALSE.
    END IF
    !
    IF (hasAttribute(xml_node, "bravais_index")) THEN
      CALL extractDataAttribute(xml_node, "bravais_index", obj%bravais_index)
      obj%bravais_index_ispresent = .TRUE.
    ELSE
      obj%bravais_index_ispresent = .FALSE.
    END IF
    !
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "atomic_positions")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:atomic_structureType","atomic_positions: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:atomic_structureType","atomic_positions: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%atomic_positions_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_atomic_positions(tmp_node, obj%atomic_positions, ierr )
    ELSE
       obj%atomic_positions_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "wyckoff_positions")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:atomic_structureType","wyckoff_positions: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:atomic_structureType","wyckoff_positions: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%wyckoff_positions_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_wyckoff_positions(tmp_node, obj%wyckoff_positions, ierr )
    ELSE
       obj%wyckoff_positions_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "crystal_positions")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:atomic_structureType","crystal_positions: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:atomic_structureType","crystal_positions: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%crystal_positions_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_atomic_positions(tmp_node, obj%crystal_positions, ierr )
    ELSE
       obj%crystal_positions_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "cell")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:atomic_structureType","cell: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:atomic_structureType","cell: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_cell(tmp_node, obj%cell, ierr )
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_atomic_structure
  !
  !
  SUBROUTINE qes_read_atomic_positions(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(atomic_positions_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "atom")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size < 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:atomic_positionsType","atom: not enough elements")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:atomic_positionsType","atom: not enough elements",10)
        END IF
    END IF
    !
    obj%ndim_atom = tmp_node_list_size
    ALLOCATE(obj%atom(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_atom(tmp_node, obj%atom(index), ierr )
    END DO
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_atomic_positions
  !
  !
  SUBROUTINE qes_read_atom(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(atom_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "name")) THEN
      CALL extractDataAttribute(xml_node, "name", obj%name)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: atomType",&
                        "required attribute name not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: atomType",&
                      "required attribute name not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "position")) THEN
      CALL extractDataAttribute(xml_node, "position", obj%position)
      obj%position_ispresent = .TRUE.
    ELSE
      obj%position_ispresent = .FALSE.
    END IF
    !
    IF (hasAttribute(xml_node, "index")) THEN
      CALL extractDataAttribute(xml_node, "index", obj%index)
      obj%index_ispresent = .TRUE.
    ELSE
      obj%index_ispresent = .FALSE.
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%atom )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_atom
  !
  !
  SUBROUTINE qes_read_wyckoff_positions(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(wyckoff_positions_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "space_group")) THEN
      CALL extractDataAttribute(xml_node, "space_group", obj%space_group)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: wyckoff_positionsType",&
                        "required attribute space_group not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: wyckoff_positionsType",&
                      "required attribute space_group not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "more_options")) THEN
      CALL extractDataAttribute(xml_node, "more_options", obj%more_options)
      obj%more_options_ispresent = .TRUE.
    ELSE
      obj%more_options_ispresent = .FALSE.
    END IF
    !
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "atom")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size < 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:wyckoff_positionsType","atom: not enough elements")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:wyckoff_positionsType","atom: not enough elements",10)
        END IF
    END IF
    !
    obj%ndim_atom = tmp_node_list_size
    ALLOCATE(obj%atom(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_atom(tmp_node, obj%atom(index), ierr )
    END DO
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_wyckoff_positions
  !
  !
  SUBROUTINE qes_read_cell(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(cell_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "a1")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:cellType","a1: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:cellType","a1: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%a1, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:cellType","error reading a1")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:cellType","error reading a1",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "a2")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:cellType","a2: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:cellType","a2: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%a2, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:cellType","error reading a2")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:cellType","error reading a2",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "a3")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:cellType","a3: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:cellType","a3: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%a3, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:cellType","error reading a3")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:cellType","error reading a3",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_cell
  !
  !
  SUBROUTINE qes_read_dft(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(dft_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "functional")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:dftType","functional: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:dftType","functional: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%functional, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:dftType","error reading functional")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:dftType","error reading functional",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "hybrid")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:dftType","hybrid: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:dftType","hybrid: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%hybrid_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_hybrid(tmp_node, obj%hybrid, ierr )
    ELSE
       obj%hybrid_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "dftU")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:dftType","dftU: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:dftType","dftU: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%dftU_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_dftU(tmp_node, obj%dftU, ierr )
    ELSE
       obj%dftU_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "vdW")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:dftType","vdW: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:dftType","vdW: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%vdW_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_vdW(tmp_node, obj%vdW, ierr )
    ELSE
       obj%vdW_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_dft
  !
  !
  SUBROUTINE qes_read_hybrid(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(hybrid_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "qpoint_grid")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:hybridType","qpoint_grid: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:hybridType","qpoint_grid: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_qpoint_grid(tmp_node, obj%qpoint_grid, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "ecutfock")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:hybridType","ecutfock: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:hybridType","ecutfock: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ecutfock, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:hybridType","error reading ecutfock")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:hybridType","error reading ecutfock",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "exx_fraction")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:hybridType","exx_fraction: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:hybridType","exx_fraction: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%exx_fraction, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:hybridType","error reading exx_fraction")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:hybridType","error reading exx_fraction",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "screening_parameter")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:hybridType","screening_parameter: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:hybridType","screening_parameter: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%screening_parameter, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:hybridType","error reading screening_parameter")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:hybridType","error reading screening_parameter",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "exxdiv_treatment")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:hybridType","exxdiv_treatment: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:hybridType","exxdiv_treatment: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%exxdiv_treatment, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:hybridType","error reading exxdiv_treatment")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:hybridType","error reading exxdiv_treatment",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "x_gamma_extrapolation")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:hybridType","x_gamma_extrapolation: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:hybridType","x_gamma_extrapolation: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%x_gamma_extrapolation, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:hybridType","error reading x_gamma_extrapolation")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:hybridType","error reading x_gamma_extrapolation",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ecutvcut")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:hybridType","ecutvcut: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:hybridType","ecutvcut: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ecutvcut, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:hybridType","error reading ecutvcut")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:hybridType","error reading ecutvcut",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_hybrid
  !
  !
  SUBROUTINE qes_read_qpoint_grid(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(qpoint_grid_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "nqx1")) THEN
      CALL extractDataAttribute(xml_node, "nqx1", obj%nqx1)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: qpoint_gridType",&
                        "required attribute nqx1 not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: qpoint_gridType",&
                      "required attribute nqx1 not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "nqx2")) THEN
      CALL extractDataAttribute(xml_node, "nqx2", obj%nqx2)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: qpoint_gridType",&
                        "required attribute nqx2 not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: qpoint_gridType",&
                      "required attribute nqx2 not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "nqx3")) THEN
      CALL extractDataAttribute(xml_node, "nqx3", obj%nqx3)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: qpoint_gridType",&
                        "required attribute nqx3 not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: qpoint_gridType",&
                      "required attribute nqx3 not found", 10 ) 
      END IF
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%qpoint_grid )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_qpoint_grid
  !
  !
  SUBROUTINE qes_read_dftU(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(dftU_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "lda_plus_u_kind")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:dftUType","lda_plus_u_kind: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:dftUType","lda_plus_u_kind: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%lda_plus_u_kind_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%lda_plus_u_kind , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:dftUType","error reading lda_plus_u_kind")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:dftUType","error reading lda_plus_u_kind",10)
         END IF
      END IF
    ELSE
       obj%lda_plus_u_kind_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "Hubbard_U")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    !
    IF (tmp_node_list_size>0) THEN
      obj%Hubbard_U_ispresent = .TRUE.
    ELSE
      obj%Hubbard_U_ispresent = .FALSE.
    END IF
    obj%ndim_Hubbard_U = tmp_node_list_size
    ALLOCATE(obj%Hubbard_U(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_HubbardCommon(tmp_node, obj%Hubbard_U(index), ierr )
    END DO
    !
    tmp_node_list => getElementsByTagname(xml_node, "Hubbard_J0")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    !
    IF (tmp_node_list_size>0) THEN
      obj%Hubbard_J0_ispresent = .TRUE.
    ELSE
      obj%Hubbard_J0_ispresent = .FALSE.
    END IF
    obj%ndim_Hubbard_J0 = tmp_node_list_size
    ALLOCATE(obj%Hubbard_J0(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_HubbardCommon(tmp_node, obj%Hubbard_J0(index), ierr )
    END DO
    !
    tmp_node_list => getElementsByTagname(xml_node, "Hubbard_alpha")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    !
    IF (tmp_node_list_size>0) THEN
      obj%Hubbard_alpha_ispresent = .TRUE.
    ELSE
      obj%Hubbard_alpha_ispresent = .FALSE.
    END IF
    obj%ndim_Hubbard_alpha = tmp_node_list_size
    ALLOCATE(obj%Hubbard_alpha(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_HubbardCommon(tmp_node, obj%Hubbard_alpha(index), ierr )
    END DO
    !
    tmp_node_list => getElementsByTagname(xml_node, "Hubbard_beta")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    !
    IF (tmp_node_list_size>0) THEN
      obj%Hubbard_beta_ispresent = .TRUE.
    ELSE
      obj%Hubbard_beta_ispresent = .FALSE.
    END IF
    obj%ndim_Hubbard_beta = tmp_node_list_size
    ALLOCATE(obj%Hubbard_beta(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_HubbardCommon(tmp_node, obj%Hubbard_beta(index), ierr )
    END DO
    !
    tmp_node_list => getElementsByTagname(xml_node, "Hubbard_J")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    !
    IF (tmp_node_list_size>0) THEN
      obj%Hubbard_J_ispresent = .TRUE.
    ELSE
      obj%Hubbard_J_ispresent = .FALSE.
    END IF
    obj%ndim_Hubbard_J = tmp_node_list_size
    ALLOCATE(obj%Hubbard_J(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_HubbardJ(tmp_node, obj%Hubbard_J(index), ierr )
    END DO
    !
    tmp_node_list => getElementsByTagname(xml_node, "starting_ns")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    !
    IF (tmp_node_list_size>0) THEN
      obj%starting_ns_ispresent = .TRUE.
    ELSE
      obj%starting_ns_ispresent = .FALSE.
    END IF
    obj%ndim_starting_ns = tmp_node_list_size
    ALLOCATE(obj%starting_ns(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_starting_ns(tmp_node, obj%starting_ns(index), ierr )
    END DO
    !
    tmp_node_list => getElementsByTagname(xml_node, "Hubbard_ns")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    !
    IF (tmp_node_list_size>0) THEN
      obj%Hubbard_ns_ispresent = .TRUE.
    ELSE
      obj%Hubbard_ns_ispresent = .FALSE.
    END IF
    obj%ndim_Hubbard_ns = tmp_node_list_size
    ALLOCATE(obj%Hubbard_ns(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_Hubbard_ns(tmp_node, obj%Hubbard_ns(index), ierr )
    END DO
    !
    tmp_node_list => getElementsByTagname(xml_node, "U_projection_type")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:dftUType","U_projection_type: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:dftUType","U_projection_type: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%U_projection_type_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%U_projection_type , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:dftUType","error reading U_projection_type")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:dftUType","error reading U_projection_type",10)
         END IF
      END IF
    ELSE
       obj%U_projection_type_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_dftU
  !
  !
  SUBROUTINE qes_read_HubbardCommon(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(HubbardCommon_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "specie")) THEN
      CALL extractDataAttribute(xml_node, "specie", obj%specie)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: HubbardCommonType",&
                        "required attribute specie not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: HubbardCommonType",&
                      "required attribute specie not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "label")) THEN
      CALL extractDataAttribute(xml_node, "label", obj%label)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: HubbardCommonType",&
                        "required attribute label not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: HubbardCommonType",&
                      "required attribute label not found", 10 ) 
      END IF
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%HubbardCommon )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_HubbardCommon
  !
  !
  SUBROUTINE qes_read_HubbardJ(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(HubbardJ_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "specie")) THEN
      CALL extractDataAttribute(xml_node, "specie", obj%specie)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: HubbardJType",&
                        "required attribute specie not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: HubbardJType",&
                      "required attribute specie not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "label")) THEN
      CALL extractDataAttribute(xml_node, "label", obj%label)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: HubbardJType",&
                        "required attribute label not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: HubbardJType",&
                      "required attribute label not found", 10 ) 
      END IF
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%HubbardJ )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_HubbardJ
  !
  !
  SUBROUTINE qes_read_starting_ns(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(starting_ns_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "specie")) THEN
      CALL extractDataAttribute(xml_node, "specie", obj%specie)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: starting_nsType",&
                        "required attribute specie not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: starting_nsType",&
                      "required attribute specie not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "label")) THEN
      CALL extractDataAttribute(xml_node, "label", obj%label)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: starting_nsType",&
                        "required attribute label not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: starting_nsType",&
                      "required attribute label not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "spin")) THEN
      CALL extractDataAttribute(xml_node, "spin", obj%spin)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: starting_nsType",&
                        "required attribute spin not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: starting_nsType",&
                      "required attribute spin not found", 10 ) 
      END IF
    END IF
    !
    !
    IF (hasAttribute(xml_node, "size"))  THEN
        CALL extractDataAttribute(xml_node, "size", obj%size)
    ELSE
        CALL errore ("qes_read: starting_nsType", &
                     "mandatory size attribute not found in "//TRIM(obj%tagname), 12)
    END IF
    !
    !
    !
    ALLOCATE (obj%starting_ns(obj%size))
    CALL extractDataContent(xml_node, obj%starting_ns )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_starting_ns
  !
  !
  SUBROUTINE qes_read_Hubbard_ns(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(Hubbard_ns_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    INTEGER :: i, length
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "specie")) THEN
      CALL extractDataAttribute(xml_node, "specie", obj%specie)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: Hubbard_nsType",&
                        "required attribute specie not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: Hubbard_nsType",&
                      "required attribute specie not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "label")) THEN
      CALL extractDataAttribute(xml_node, "label", obj%label)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: Hubbard_nsType",&
                        "required attribute label not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: Hubbard_nsType",&
                      "required attribute label not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "spin")) THEN
      CALL extractDataAttribute(xml_node, "spin", obj%spin)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: Hubbard_nsType",&
                        "required attribute spin not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: Hubbard_nsType",&
                      "required attribute spin not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "index")) THEN
      CALL extractDataAttribute(xml_node, "index", obj%index)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: Hubbard_nsType",&
                        "required attribute index not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: Hubbard_nsType",&
                      "required attribute index not found", 10 ) 
      END IF
    END IF
    !
    !
    IF (hasAttribute(xml_node, "rank"))  THEN
        CALL extractDataAttribute(xml_node, "rank", obj%rank)
    ELSE
        CALL errore ("qes_read: Hubbard_nsType",&
                      "required attribute rank not found, can't read further, stopping", 10 ) 
    END IF
    ALLOCATE (obj%dims(obj%rank))
    IF (hasAttribute(xml_node, "dims")) THEN
        CALL extractDataAttribute(xml_node, "dims", obj%dims)
    ELSE
        CALL errore ("qes_read: Hubbard_nsType",&
                      "required attribute dims not found, can't read further, stopping", 10 ) 
    END IF
    IF (hasAttribute(xml_node,"order")) THEN
        CALL extractDataAttribute(xml_node, "order", obj%order)
    ELSE
        obj%order = "F"
    END IF
    !
    !
    !
    length = 1
    DO i =1, obj%rank
       length = length * obj%dims(i)
    END DO
    ALLOCATE (obj%Hubbard_ns(length) )
    CALL extractDataContent(xml_node, obj%Hubbard_ns )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_Hubbard_ns
  !
  !
  SUBROUTINE qes_read_vdW(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(vdW_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "vdw_corr")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:vdWType","vdw_corr: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:vdWType","vdw_corr: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%vdw_corr, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:vdWType","error reading vdw_corr")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:vdWType","error reading vdw_corr",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "non_local_term")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:vdWType","non_local_term: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:vdWType","non_local_term: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%non_local_term_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%non_local_term , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:vdWType","error reading non_local_term")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:vdWType","error reading non_local_term",10)
         END IF
      END IF
    ELSE
       obj%non_local_term_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "london_s6")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:vdWType","london_s6: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:vdWType","london_s6: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%london_s6_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%london_s6 , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:vdWType","error reading london_s6")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:vdWType","error reading london_s6",10)
         END IF
      END IF
    ELSE
       obj%london_s6_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ts_vdw_econv_thr")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:vdWType","ts_vdw_econv_thr: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:vdWType","ts_vdw_econv_thr: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%ts_vdw_econv_thr_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%ts_vdw_econv_thr , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:vdWType","error reading ts_vdw_econv_thr")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:vdWType","error reading ts_vdw_econv_thr",10)
         END IF
      END IF
    ELSE
       obj%ts_vdw_econv_thr_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ts_vdw_isolated")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:vdWType","ts_vdw_isolated: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:vdWType","ts_vdw_isolated: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%ts_vdw_isolated_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%ts_vdw_isolated , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:vdWType","error reading ts_vdw_isolated")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:vdWType","error reading ts_vdw_isolated",10)
         END IF
      END IF
    ELSE
       obj%ts_vdw_isolated_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "london_rcut")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:vdWType","london_rcut: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:vdWType","london_rcut: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%london_rcut_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%london_rcut , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:vdWType","error reading london_rcut")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:vdWType","error reading london_rcut",10)
         END IF
      END IF
    ELSE
       obj%london_rcut_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "xdm_a1")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:vdWType","xdm_a1: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:vdWType","xdm_a1: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%xdm_a1_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%xdm_a1 , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:vdWType","error reading xdm_a1")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:vdWType","error reading xdm_a1",10)
         END IF
      END IF
    ELSE
       obj%xdm_a1_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "xdm_a2")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:vdWType","xdm_a2: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:vdWType","xdm_a2: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%xdm_a2_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%xdm_a2 , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:vdWType","error reading xdm_a2")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:vdWType","error reading xdm_a2",10)
         END IF
      END IF
    ELSE
       obj%xdm_a2_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "london_c6")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    !
    IF (tmp_node_list_size>0) THEN
      obj%london_c6_ispresent = .TRUE.
    ELSE
      obj%london_c6_ispresent = .FALSE.
    END IF
    obj%ndim_london_c6 = tmp_node_list_size
    ALLOCATE(obj%london_c6(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_HubbardCommon(tmp_node, obj%london_c6(index), ierr )
    END DO
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_vdW
  !
  !
  SUBROUTINE qes_read_spin(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(spin_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "lsda")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:spinType","lsda: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:spinType","lsda: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%lsda, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:spinType","error reading lsda")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:spinType","error reading lsda",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "noncolin")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:spinType","noncolin: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:spinType","noncolin: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%noncolin, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:spinType","error reading noncolin")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:spinType","error reading noncolin",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "spinorbit")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:spinType","spinorbit: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:spinType","spinorbit: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%spinorbit, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:spinType","error reading spinorbit")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:spinType","error reading spinorbit",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_spin
  !
  !
  SUBROUTINE qes_read_bands(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(bands_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "nbnd")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:bandsType","nbnd: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:bandsType","nbnd: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%nbnd_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%nbnd , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:bandsType","error reading nbnd")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:bandsType","error reading nbnd",10)
         END IF
      END IF
    ELSE
       obj%nbnd_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "smearing")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:bandsType","smearing: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:bandsType","smearing: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%smearing_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_smearing(tmp_node, obj%smearing, ierr )
    ELSE
       obj%smearing_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "tot_charge")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:bandsType","tot_charge: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:bandsType","tot_charge: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%tot_charge_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%tot_charge , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:bandsType","error reading tot_charge")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:bandsType","error reading tot_charge",10)
         END IF
      END IF
    ELSE
       obj%tot_charge_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "tot_magnetization")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:bandsType","tot_magnetization: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:bandsType","tot_magnetization: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%tot_magnetization_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%tot_magnetization , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:bandsType","error reading tot_magnetization")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:bandsType","error reading tot_magnetization",10)
         END IF
      END IF
    ELSE
       obj%tot_magnetization_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "occupations")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:bandsType","occupations: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:bandsType","occupations: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_occupations(tmp_node, obj%occupations, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "inputOccupations")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 2) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:bandsType","inputOccupations: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:bandsType","inputOccupations: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%inputOccupations_ispresent = .TRUE.
    ELSE
      obj%inputOccupations_ispresent = .FALSE.
    END IF
    obj%ndim_inputOccupations = tmp_node_list_size
    ALLOCATE(obj%inputOccupations(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_inputOccupations(tmp_node, obj%inputOccupations(index), ierr )
    END DO
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_bands
  !
  !
  SUBROUTINE qes_read_smearing(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(smearing_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "degauss")) THEN
      CALL extractDataAttribute(xml_node, "degauss", obj%degauss)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: smearingType",&
                        "required attribute degauss not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: smearingType",&
                      "required attribute degauss not found", 10 ) 
      END IF
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%smearing )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_smearing
  !
  !
  SUBROUTINE qes_read_occupations(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(occupations_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "spin")) THEN
      CALL extractDataAttribute(xml_node, "spin", obj%spin)
      obj%spin_ispresent = .TRUE.
    ELSE
      obj%spin_ispresent = .FALSE.
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%occupations )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_occupations
  !
  !
  SUBROUTINE qes_read_basis(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(basis_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "gamma_only")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basisType","gamma_only: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basisType","gamma_only: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%gamma_only_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%gamma_only , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:basisType","error reading gamma_only")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:basisType","error reading gamma_only",10)
         END IF
      END IF
    ELSE
       obj%gamma_only_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ecutwfc")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basisType","ecutwfc: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basisType","ecutwfc: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ecutwfc, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:basisType","error reading ecutwfc")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:basisType","error reading ecutwfc",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ecutrho")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basisType","ecutrho: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basisType","ecutrho: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%ecutrho_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%ecutrho , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:basisType","error reading ecutrho")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:basisType","error reading ecutrho",10)
         END IF
      END IF
    ELSE
       obj%ecutrho_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "fft_grid")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basisType","fft_grid: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basisType","fft_grid: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%fft_grid_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_basisSetItem(tmp_node, obj%fft_grid, ierr )
    ELSE
       obj%fft_grid_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "fft_smooth")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basisType","fft_smooth: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basisType","fft_smooth: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%fft_smooth_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_basisSetItem(tmp_node, obj%fft_smooth, ierr )
    ELSE
       obj%fft_smooth_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "fft_box")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basisType","fft_box: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basisType","fft_box: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%fft_box_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_basisSetItem(tmp_node, obj%fft_box, ierr )
    ELSE
       obj%fft_box_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_basis
  !
  !
  SUBROUTINE qes_read_basis_set(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(basis_set_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "gamma_only")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basis_setType","gamma_only: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basis_setType","gamma_only: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%gamma_only_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%gamma_only , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:basis_setType","error reading gamma_only")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:basis_setType","error reading gamma_only",10)
         END IF
      END IF
    ELSE
       obj%gamma_only_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ecutwfc")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basis_setType","ecutwfc: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basis_setType","ecutwfc: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ecutwfc, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:basis_setType","error reading ecutwfc")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:basis_setType","error reading ecutwfc",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ecutrho")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basis_setType","ecutrho: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basis_setType","ecutrho: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%ecutrho_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%ecutrho , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:basis_setType","error reading ecutrho")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:basis_setType","error reading ecutrho",10)
         END IF
      END IF
    ELSE
       obj%ecutrho_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "fft_grid")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basis_setType","fft_grid: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basis_setType","fft_grid: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_basisSetItem(tmp_node, obj%fft_grid, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "fft_smooth")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basis_setType","fft_smooth: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basis_setType","fft_smooth: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%fft_smooth_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_basisSetItem(tmp_node, obj%fft_smooth, ierr )
    ELSE
       obj%fft_smooth_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "fft_box")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basis_setType","fft_box: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basis_setType","fft_box: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%fft_box_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_basisSetItem(tmp_node, obj%fft_box, ierr )
    ELSE
       obj%fft_box_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ngm")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basis_setType","ngm: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basis_setType","ngm: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ngm, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:basis_setType","error reading ngm")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:basis_setType","error reading ngm",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ngms")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basis_setType","ngms: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basis_setType","ngms: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%ngms_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%ngms , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:basis_setType","error reading ngms")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:basis_setType","error reading ngms",10)
         END IF
      END IF
    ELSE
       obj%ngms_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "npwx")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basis_setType","npwx: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basis_setType","npwx: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%npwx, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:basis_setType","error reading npwx")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:basis_setType","error reading npwx",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "reciprocal_lattice")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:basis_setType","reciprocal_lattice: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:basis_setType","reciprocal_lattice: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_reciprocal_lattice(tmp_node, obj%reciprocal_lattice, ierr )
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_basis_set
  !
  !
  SUBROUTINE qes_read_basisSetItem(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(basisSetItem_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "nr1")) THEN
      CALL extractDataAttribute(xml_node, "nr1", obj%nr1)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: basisSetItemType",&
                        "required attribute nr1 not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: basisSetItemType",&
                      "required attribute nr1 not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "nr2")) THEN
      CALL extractDataAttribute(xml_node, "nr2", obj%nr2)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: basisSetItemType",&
                        "required attribute nr2 not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: basisSetItemType",&
                      "required attribute nr2 not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "nr3")) THEN
      CALL extractDataAttribute(xml_node, "nr3", obj%nr3)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: basisSetItemType",&
                        "required attribute nr3 not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: basisSetItemType",&
                      "required attribute nr3 not found", 10 ) 
      END IF
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%basisSetItem )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_basisSetItem
  !
  !
  SUBROUTINE qes_read_reciprocal_lattice(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(reciprocal_lattice_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "b1")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:reciprocal_latticeType","b1: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:reciprocal_latticeType","b1: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%b1, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:reciprocal_latticeType","error reading b1")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:reciprocal_latticeType","error reading b1",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "b2")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:reciprocal_latticeType","b2: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:reciprocal_latticeType","b2: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%b2, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:reciprocal_latticeType","error reading b2")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:reciprocal_latticeType","error reading b2",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "b3")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:reciprocal_latticeType","b3: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:reciprocal_latticeType","b3: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%b3, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:reciprocal_latticeType","error reading b3")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:reciprocal_latticeType","error reading b3",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_reciprocal_lattice
  !
  !
  SUBROUTINE qes_read_electron_control(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(electron_control_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "diagonalization")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","diagonalization: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","diagonalization: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%diagonalization, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:electron_controlType","error reading diagonalization")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:electron_controlType","error reading diagonalization",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "mixing_mode")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","mixing_mode: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","mixing_mode: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%mixing_mode, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:electron_controlType","error reading mixing_mode")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:electron_controlType","error reading mixing_mode",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "mixing_beta")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","mixing_beta: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","mixing_beta: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%mixing_beta, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:electron_controlType","error reading mixing_beta")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:electron_controlType","error reading mixing_beta",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "conv_thr")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","conv_thr: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","conv_thr: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%conv_thr, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:electron_controlType","error reading conv_thr")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:electron_controlType","error reading conv_thr",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "mixing_ndim")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","mixing_ndim: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","mixing_ndim: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%mixing_ndim, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:electron_controlType","error reading mixing_ndim")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:electron_controlType","error reading mixing_ndim",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "max_nstep")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","max_nstep: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","max_nstep: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%max_nstep, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:electron_controlType","error reading max_nstep")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:electron_controlType","error reading max_nstep",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "real_space_q")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","real_space_q: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","real_space_q: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%real_space_q, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:electron_controlType","error reading real_space_q")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:electron_controlType","error reading real_space_q",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "tq_smoothing")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","tq_smoothing: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","tq_smoothing: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%tq_smoothing, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:electron_controlType","error reading tq_smoothing")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:electron_controlType","error reading tq_smoothing",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "tbeta_smoothing")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","tbeta_smoothing: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","tbeta_smoothing: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%tbeta_smoothing, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:electron_controlType","error reading tbeta_smoothing")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:electron_controlType","error reading tbeta_smoothing",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "diago_thr_init")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","diago_thr_init: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","diago_thr_init: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%diago_thr_init, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:electron_controlType","error reading diago_thr_init")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:electron_controlType","error reading diago_thr_init",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "diago_full_acc")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","diago_full_acc: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","diago_full_acc: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%diago_full_acc, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:electron_controlType","error reading diago_full_acc")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:electron_controlType","error reading diago_full_acc",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "diago_cg_maxiter")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","diago_cg_maxiter: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","diago_cg_maxiter: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%diago_cg_maxiter_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%diago_cg_maxiter , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:electron_controlType","error reading diago_cg_maxiter")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:electron_controlType","error reading diago_cg_maxiter",10)
         END IF
      END IF
    ELSE
       obj%diago_cg_maxiter_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "diago_ppcg_maxiter")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","diago_ppcg_maxiter: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","diago_ppcg_maxiter: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%diago_ppcg_maxiter_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%diago_ppcg_maxiter , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:electron_controlType","error reading diago_ppcg_maxiter")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:electron_controlType","error reading diago_ppcg_maxiter",10)
         END IF
      END IF
    ELSE
       obj%diago_ppcg_maxiter_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "diago_david_ndim")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electron_controlType","diago_david_ndim: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electron_controlType","diago_david_ndim: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%diago_david_ndim_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%diago_david_ndim , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:electron_controlType","error reading diago_david_ndim")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:electron_controlType","error reading diago_david_ndim",10)
         END IF
      END IF
    ELSE
       obj%diago_david_ndim_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_electron_control
  !
  !
  SUBROUTINE qes_read_k_points_IBZ(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(k_points_IBZ_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "monkhorst_pack")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:k_points_IBZType","monkhorst_pack: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:k_points_IBZType","monkhorst_pack: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%monkhorst_pack_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_monkhorst_pack(tmp_node, obj%monkhorst_pack, ierr )
    ELSE
       obj%monkhorst_pack_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nk")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:k_points_IBZType","nk: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:k_points_IBZType","nk: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%nk_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%nk , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:k_points_IBZType","error reading nk")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:k_points_IBZType","error reading nk",10)
         END IF
      END IF
    ELSE
       obj%nk_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "k_point")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    !
    IF (tmp_node_list_size>0) THEN
      obj%k_point_ispresent = .TRUE.
    ELSE
      obj%k_point_ispresent = .FALSE.
    END IF
    obj%ndim_k_point = tmp_node_list_size
    ALLOCATE(obj%k_point(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_k_point(tmp_node, obj%k_point(index), ierr )
    END DO
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_k_points_IBZ
  !
  !
  SUBROUTINE qes_read_monkhorst_pack(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(monkhorst_pack_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "nk1")) THEN
      CALL extractDataAttribute(xml_node, "nk1", obj%nk1)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: monkhorst_packType",&
                        "required attribute nk1 not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: monkhorst_packType",&
                      "required attribute nk1 not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "nk2")) THEN
      CALL extractDataAttribute(xml_node, "nk2", obj%nk2)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: monkhorst_packType",&
                        "required attribute nk2 not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: monkhorst_packType",&
                      "required attribute nk2 not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "nk3")) THEN
      CALL extractDataAttribute(xml_node, "nk3", obj%nk3)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: monkhorst_packType",&
                        "required attribute nk3 not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: monkhorst_packType",&
                      "required attribute nk3 not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "k1")) THEN
      CALL extractDataAttribute(xml_node, "k1", obj%k1)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: monkhorst_packType",&
                        "required attribute k1 not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: monkhorst_packType",&
                      "required attribute k1 not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "k2")) THEN
      CALL extractDataAttribute(xml_node, "k2", obj%k2)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: monkhorst_packType",&
                        "required attribute k2 not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: monkhorst_packType",&
                      "required attribute k2 not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "k3")) THEN
      CALL extractDataAttribute(xml_node, "k3", obj%k3)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: monkhorst_packType",&
                        "required attribute k3 not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: monkhorst_packType",&
                      "required attribute k3 not found", 10 ) 
      END IF
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%monkhorst_pack )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_monkhorst_pack
  !
  !
  SUBROUTINE qes_read_k_point(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(k_point_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "weight")) THEN
      CALL extractDataAttribute(xml_node, "weight", obj%weight)
      obj%weight_ispresent = .TRUE.
    ELSE
      obj%weight_ispresent = .FALSE.
    END IF
    !
    IF (hasAttribute(xml_node, "label")) THEN
      CALL extractDataAttribute(xml_node, "label", obj%label)
      obj%label_ispresent = .TRUE.
    ELSE
      obj%label_ispresent = .FALSE.
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%k_point )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_k_point
  !
  !
  SUBROUTINE qes_read_ion_control(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(ion_control_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "ion_dynamics")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ion_controlType","ion_dynamics: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ion_controlType","ion_dynamics: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ion_dynamics, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:ion_controlType","error reading ion_dynamics")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:ion_controlType","error reading ion_dynamics",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "upscale")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ion_controlType","upscale: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ion_controlType","upscale: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%upscale_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%upscale , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:ion_controlType","error reading upscale")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:ion_controlType","error reading upscale",10)
         END IF
      END IF
    ELSE
       obj%upscale_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "remove_rigid_rot")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ion_controlType","remove_rigid_rot: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ion_controlType","remove_rigid_rot: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%remove_rigid_rot_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%remove_rigid_rot , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:ion_controlType","error reading remove_rigid_rot")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:ion_controlType","error reading remove_rigid_rot",10)
         END IF
      END IF
    ELSE
       obj%remove_rigid_rot_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "refold_pos")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ion_controlType","refold_pos: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ion_controlType","refold_pos: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%refold_pos_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%refold_pos , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:ion_controlType","error reading refold_pos")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:ion_controlType","error reading refold_pos",10)
         END IF
      END IF
    ELSE
       obj%refold_pos_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "bfgs")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ion_controlType","bfgs: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ion_controlType","bfgs: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%bfgs_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_bfgs(tmp_node, obj%bfgs, ierr )
    ELSE
       obj%bfgs_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "md")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ion_controlType","md: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ion_controlType","md: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%md_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_md(tmp_node, obj%md, ierr )
    ELSE
       obj%md_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_ion_control
  !
  !
  SUBROUTINE qes_read_bfgs(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(bfgs_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "ndim")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:bfgsType","ndim: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:bfgsType","ndim: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ndim, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:bfgsType","error reading ndim")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:bfgsType","error reading ndim",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "trust_radius_min")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:bfgsType","trust_radius_min: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:bfgsType","trust_radius_min: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%trust_radius_min, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:bfgsType","error reading trust_radius_min")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:bfgsType","error reading trust_radius_min",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "trust_radius_max")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:bfgsType","trust_radius_max: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:bfgsType","trust_radius_max: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%trust_radius_max, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:bfgsType","error reading trust_radius_max")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:bfgsType","error reading trust_radius_max",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "trust_radius_init")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:bfgsType","trust_radius_init: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:bfgsType","trust_radius_init: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%trust_radius_init, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:bfgsType","error reading trust_radius_init")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:bfgsType","error reading trust_radius_init",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "w1")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:bfgsType","w1: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:bfgsType","w1: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%w1, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:bfgsType","error reading w1")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:bfgsType","error reading w1",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "w2")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:bfgsType","w2: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:bfgsType","w2: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%w2, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:bfgsType","error reading w2")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:bfgsType","error reading w2",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_bfgs
  !
  !
  SUBROUTINE qes_read_md(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(md_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "pot_extrapolation")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:mdType","pot_extrapolation: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:mdType","pot_extrapolation: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%pot_extrapolation, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:mdType","error reading pot_extrapolation")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:mdType","error reading pot_extrapolation",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "wfc_extrapolation")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:mdType","wfc_extrapolation: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:mdType","wfc_extrapolation: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%wfc_extrapolation, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:mdType","error reading wfc_extrapolation")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:mdType","error reading wfc_extrapolation",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ion_temperature")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:mdType","ion_temperature: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:mdType","ion_temperature: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ion_temperature, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:mdType","error reading ion_temperature")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:mdType","error reading ion_temperature",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "timestep")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:mdType","timestep: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:mdType","timestep: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%timestep, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:mdType","error reading timestep")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:mdType","error reading timestep",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "tempw")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:mdType","tempw: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:mdType","tempw: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%tempw, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:mdType","error reading tempw")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:mdType","error reading tempw",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "tolp")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:mdType","tolp: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:mdType","tolp: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%tolp, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:mdType","error reading tolp")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:mdType","error reading tolp",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "deltaT")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:mdType","deltaT: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:mdType","deltaT: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%deltaT, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:mdType","error reading deltaT")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:mdType","error reading deltaT",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nraise")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:mdType","nraise: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:mdType","nraise: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nraise, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:mdType","error reading nraise")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:mdType","error reading nraise",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_md
  !
  !
  SUBROUTINE qes_read_cell_control(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(cell_control_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "cell_dynamics")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:cell_controlType","cell_dynamics: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:cell_controlType","cell_dynamics: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%cell_dynamics, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:cell_controlType","error reading cell_dynamics")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:cell_controlType","error reading cell_dynamics",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "pressure")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:cell_controlType","pressure: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:cell_controlType","pressure: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%pressure, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:cell_controlType","error reading pressure")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:cell_controlType","error reading pressure",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "wmass")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:cell_controlType","wmass: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:cell_controlType","wmass: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%wmass_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%wmass , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:cell_controlType","error reading wmass")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:cell_controlType","error reading wmass",10)
         END IF
      END IF
    ELSE
       obj%wmass_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "cell_factor")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:cell_controlType","cell_factor: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:cell_controlType","cell_factor: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%cell_factor_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%cell_factor , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:cell_controlType","error reading cell_factor")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:cell_controlType","error reading cell_factor",10)
         END IF
      END IF
    ELSE
       obj%cell_factor_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "fix_volume")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:cell_controlType","fix_volume: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:cell_controlType","fix_volume: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%fix_volume_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%fix_volume , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:cell_controlType","error reading fix_volume")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:cell_controlType","error reading fix_volume",10)
         END IF
      END IF
    ELSE
       obj%fix_volume_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "fix_area")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:cell_controlType","fix_area: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:cell_controlType","fix_area: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%fix_area_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%fix_area , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:cell_controlType","error reading fix_area")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:cell_controlType","error reading fix_area",10)
         END IF
      END IF
    ELSE
       obj%fix_area_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "isotropic")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:cell_controlType","isotropic: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:cell_controlType","isotropic: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%isotropic_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%isotropic , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:cell_controlType","error reading isotropic")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:cell_controlType","error reading isotropic",10)
         END IF
      END IF
    ELSE
       obj%isotropic_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "free_cell")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:cell_controlType","free_cell: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:cell_controlType","free_cell: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%free_cell_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_integerMatrix(tmp_node, obj%free_cell, ierr )
    ELSE
       obj%free_cell_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_cell_control
  !
  !
  SUBROUTINE qes_read_symmetry_flags(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(symmetry_flags_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "nosym")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetry_flagsType","nosym: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetry_flagsType","nosym: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nosym, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:symmetry_flagsType","error reading nosym")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:symmetry_flagsType","error reading nosym",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nosym_evc")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetry_flagsType","nosym_evc: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetry_flagsType","nosym_evc: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nosym_evc, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:symmetry_flagsType","error reading nosym_evc")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:symmetry_flagsType","error reading nosym_evc",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "noinv")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetry_flagsType","noinv: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetry_flagsType","noinv: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%noinv, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:symmetry_flagsType","error reading noinv")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:symmetry_flagsType","error reading noinv",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "no_t_rev")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetry_flagsType","no_t_rev: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetry_flagsType","no_t_rev: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%no_t_rev, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:symmetry_flagsType","error reading no_t_rev")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:symmetry_flagsType","error reading no_t_rev",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "force_symmorphic")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetry_flagsType","force_symmorphic: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetry_flagsType","force_symmorphic: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%force_symmorphic, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:symmetry_flagsType","error reading force_symmorphic")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:symmetry_flagsType","error reading force_symmorphic",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "use_all_frac")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetry_flagsType","use_all_frac: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetry_flagsType","use_all_frac: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%use_all_frac, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:symmetry_flagsType","error reading use_all_frac")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:symmetry_flagsType","error reading use_all_frac",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_symmetry_flags
  !
  !
  SUBROUTINE qes_read_boundary_conditions(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(boundary_conditions_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "assume_isolated")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:boundary_conditionsType","assume_isolated: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:boundary_conditionsType","assume_isolated: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%assume_isolated, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:boundary_conditionsType","error reading assume_isolated")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:boundary_conditionsType","error reading assume_isolated",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "esm")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:boundary_conditionsType","esm: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:boundary_conditionsType","esm: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%esm_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_esm(tmp_node, obj%esm, ierr )
    ELSE
       obj%esm_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "fcp_opt")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:boundary_conditionsType","fcp_opt: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:boundary_conditionsType","fcp_opt: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%fcp_opt_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%fcp_opt , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:boundary_conditionsType","error reading fcp_opt")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:boundary_conditionsType","error reading fcp_opt",10)
         END IF
      END IF
    ELSE
       obj%fcp_opt_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "fcp_mu")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:boundary_conditionsType","fcp_mu: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:boundary_conditionsType","fcp_mu: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%fcp_mu_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%fcp_mu , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:boundary_conditionsType","error reading fcp_mu")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:boundary_conditionsType","error reading fcp_mu",10)
         END IF
      END IF
    ELSE
       obj%fcp_mu_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_boundary_conditions
  !
  !
  SUBROUTINE qes_read_esm(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(esm_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "bc")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:esmType","bc: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:esmType","bc: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%bc, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:esmType","error reading bc")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:esmType","error reading bc",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nfit")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:esmType","nfit: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:esmType","nfit: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nfit, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:esmType","error reading nfit")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:esmType","error reading nfit",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "w")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:esmType","w: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:esmType","w: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%w, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:esmType","error reading w")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:esmType","error reading w",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "efield")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:esmType","efield: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:esmType","efield: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%efield, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:esmType","error reading efield")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:esmType","error reading efield",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_esm
  !
  !
  SUBROUTINE qes_read_ekin_functional(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(ekin_functional_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "ecfixed")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ekin_functionalType","ecfixed: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ekin_functionalType","ecfixed: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ecfixed, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:ekin_functionalType","error reading ecfixed")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:ekin_functionalType","error reading ecfixed",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "qcutz")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ekin_functionalType","qcutz: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ekin_functionalType","qcutz: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%qcutz, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:ekin_functionalType","error reading qcutz")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:ekin_functionalType","error reading qcutz",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "q2sigma")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ekin_functionalType","q2sigma: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ekin_functionalType","q2sigma: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%q2sigma, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:ekin_functionalType","error reading q2sigma")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:ekin_functionalType","error reading q2sigma",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_ekin_functional
  !
  !
  SUBROUTINE qes_read_spin_constraints(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(spin_constraints_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "spin_constraints")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:spin_constraintsType","spin_constraints: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:spin_constraintsType","spin_constraints: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%spin_constraints, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:spin_constraintsType","error reading spin_constraints")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:spin_constraintsType","error reading spin_constraints",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "lagrange_multiplier")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:spin_constraintsType","lagrange_multiplier: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:spin_constraintsType","lagrange_multiplier: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%lagrange_multiplier, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:spin_constraintsType","error reading lagrange_multiplier")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:spin_constraintsType","error reading lagrange_multiplier",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "target_magnetization")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:spin_constraintsType","target_magnetization: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:spin_constraintsType","target_magnetization: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%target_magnetization_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%target_magnetization , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:spin_constraintsType","error reading target_magnetization")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:spin_constraintsType","error reading target_magnetization",10)
         END IF
      END IF
    ELSE
       obj%target_magnetization_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_spin_constraints
  !
  !
  SUBROUTINE qes_read_electric_field(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(electric_field_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "electric_potential")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electric_fieldType","electric_potential: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electric_fieldType","electric_potential: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%electric_potential, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:electric_fieldType","error reading electric_potential")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:electric_fieldType","error reading electric_potential",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "dipole_correction")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electric_fieldType","dipole_correction: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electric_fieldType","dipole_correction: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%dipole_correction_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%dipole_correction , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:electric_fieldType","error reading dipole_correction")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:electric_fieldType","error reading dipole_correction",10)
         END IF
      END IF
    ELSE
       obj%dipole_correction_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "gate_settings")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electric_fieldType","gate_settings: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electric_fieldType","gate_settings: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%gate_settings_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_gate_settings(tmp_node, obj%gate_settings, ierr )
    ELSE
       obj%gate_settings_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "electric_field_direction")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electric_fieldType","electric_field_direction: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electric_fieldType","electric_field_direction: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%electric_field_direction_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%electric_field_direction , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:electric_fieldType","error reading electric_field_direction")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:electric_fieldType","error reading electric_field_direction",10)
         END IF
      END IF
    ELSE
       obj%electric_field_direction_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "potential_max_position")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electric_fieldType","potential_max_position: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electric_fieldType","potential_max_position: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%potential_max_position_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%potential_max_position , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:electric_fieldType","error reading potential_max_position")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:electric_fieldType","error reading potential_max_position",10)
         END IF
      END IF
    ELSE
       obj%potential_max_position_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "potential_decrease_width")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electric_fieldType","potential_decrease_width: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electric_fieldType","potential_decrease_width: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%potential_decrease_width_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%potential_decrease_width , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:electric_fieldType","error reading potential_decrease_width")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:electric_fieldType","error reading potential_decrease_width",10)
         END IF
      END IF
    ELSE
       obj%potential_decrease_width_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "electric_field_amplitude")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electric_fieldType","electric_field_amplitude: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electric_fieldType","electric_field_amplitude: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%electric_field_amplitude_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%electric_field_amplitude , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:electric_fieldType","error reading electric_field_amplitude")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:electric_fieldType","error reading electric_field_amplitude",10)
         END IF
      END IF
    ELSE
       obj%electric_field_amplitude_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "electric_field_vector")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electric_fieldType","electric_field_vector: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electric_fieldType","electric_field_vector: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%electric_field_vector_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%electric_field_vector , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:electric_fieldType","error reading electric_field_vector")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:electric_fieldType","error reading electric_field_vector",10)
         END IF
      END IF
    ELSE
       obj%electric_field_vector_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nk_per_string")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electric_fieldType","nk_per_string: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electric_fieldType","nk_per_string: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%nk_per_string_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%nk_per_string , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:electric_fieldType","error reading nk_per_string")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:electric_fieldType","error reading nk_per_string",10)
         END IF
      END IF
    ELSE
       obj%nk_per_string_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "n_berry_cycles")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electric_fieldType","n_berry_cycles: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electric_fieldType","n_berry_cycles: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%n_berry_cycles_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%n_berry_cycles , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:electric_fieldType","error reading n_berry_cycles")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:electric_fieldType","error reading n_berry_cycles",10)
         END IF
      END IF
    ELSE
       obj%n_berry_cycles_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_electric_field
  !
  !
  SUBROUTINE qes_read_gate_settings(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(gate_settings_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "use_gate")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gate_settingsType","use_gate: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gate_settingsType","use_gate: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%use_gate, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:gate_settingsType","error reading use_gate")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:gate_settingsType","error reading use_gate",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "zgate")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gate_settingsType","zgate: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gate_settingsType","zgate: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%zgate_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%zgate , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:gate_settingsType","error reading zgate")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:gate_settingsType","error reading zgate",10)
         END IF
      END IF
    ELSE
       obj%zgate_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "relaxz")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gate_settingsType","relaxz: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gate_settingsType","relaxz: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%relaxz_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%relaxz , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:gate_settingsType","error reading relaxz")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:gate_settingsType","error reading relaxz",10)
         END IF
      END IF
    ELSE
       obj%relaxz_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "block")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gate_settingsType","block: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gate_settingsType","block: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%block_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%block , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:gate_settingsType","error reading block")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:gate_settingsType","error reading block",10)
         END IF
      END IF
    ELSE
       obj%block_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "block_1")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gate_settingsType","block_1: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gate_settingsType","block_1: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%block_1_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%block_1 , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:gate_settingsType","error reading block_1")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:gate_settingsType","error reading block_1",10)
         END IF
      END IF
    ELSE
       obj%block_1_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "block_2")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gate_settingsType","block_2: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gate_settingsType","block_2: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%block_2_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%block_2 , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:gate_settingsType","error reading block_2")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:gate_settingsType","error reading block_2",10)
         END IF
      END IF
    ELSE
       obj%block_2_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "block_height")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gate_settingsType","block_height: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gate_settingsType","block_height: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%block_height_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%block_height , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:gate_settingsType","error reading block_height")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:gate_settingsType","error reading block_height",10)
         END IF
      END IF
    ELSE
       obj%block_height_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_gate_settings
  !
  !
  SUBROUTINE qes_read_atomic_constraints(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(atomic_constraints_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "num_of_constraints")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:atomic_constraintsType","num_of_constraints: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:atomic_constraintsType","num_of_constraints: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%num_of_constraints, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:atomic_constraintsType","error reading num_of_constraints")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:atomic_constraintsType","error reading num_of_constraints",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "tolerance")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:atomic_constraintsType","tolerance: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:atomic_constraintsType","tolerance: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%tolerance, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:atomic_constraintsType","error reading tolerance")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:atomic_constraintsType","error reading tolerance",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "atomic_constraint")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size < 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:atomic_constraintsType","atomic_constraint: not enough elements")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:atomic_constraintsType","atomic_constraint: not enough elements",10)
        END IF
    END IF
    !
    obj%ndim_atomic_constraint = tmp_node_list_size
    ALLOCATE(obj%atomic_constraint(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_atomic_constraint(tmp_node, obj%atomic_constraint(index), ierr )
    END DO
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_atomic_constraints
  !
  !
  SUBROUTINE qes_read_atomic_constraint(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(atomic_constraint_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "constr_parms")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:atomic_constraintType","constr_parms: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:atomic_constraintType","constr_parms: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%constr_parms, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:atomic_constraintType","error reading constr_parms")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:atomic_constraintType","error reading constr_parms",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "constr_type")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:atomic_constraintType","constr_type: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:atomic_constraintType","constr_type: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%constr_type, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:atomic_constraintType","error reading constr_type")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:atomic_constraintType","error reading constr_type",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "constr_target")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:atomic_constraintType","constr_target: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:atomic_constraintType","constr_target: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%constr_target, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:atomic_constraintType","error reading constr_target")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:atomic_constraintType","error reading constr_target",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_atomic_constraint
  !
  !
  SUBROUTINE qes_read_inputOccupations(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(inputOccupations_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "ispin")) THEN
      CALL extractDataAttribute(xml_node, "ispin", obj%ispin)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: inputOccupationsType",&
                        "required attribute ispin not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: inputOccupationsType",&
                      "required attribute ispin not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "spin_factor")) THEN
      CALL extractDataAttribute(xml_node, "spin_factor", obj%spin_factor)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: inputOccupationsType",&
                        "required attribute spin_factor not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: inputOccupationsType",&
                      "required attribute spin_factor not found", 10 ) 
      END IF
    END IF
    !
    !
    IF (hasAttribute(xml_node, "size"))  THEN
        CALL extractDataAttribute(xml_node, "size", obj%size)
    ELSE
        CALL errore ("qes_read: inputOccupationsType", &
                     "mandatory size attribute not found in "//TRIM(obj%tagname), 12)
    END IF
    !
    !
    !
    ALLOCATE (obj%inputOccupations(obj%size))
    CALL extractDataContent(xml_node, obj%inputOccupations )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_inputOccupations
  !
  !
  SUBROUTINE qes_read_outputElectricField(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(outputElectricField_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "BerryPhase")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputElectricFieldType","BerryPhase: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputElectricFieldType","BerryPhase: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%BerryPhase_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_BerryPhaseOutput(tmp_node, obj%BerryPhase, ierr )
    ELSE
       obj%BerryPhase_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "finiteElectricFieldInfo")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputElectricFieldType","finiteElectricFieldInfo: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputElectricFieldType","finiteElectricFieldInfo: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%finiteElectricFieldInfo_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_finiteFieldOut(tmp_node, obj%finiteElectricFieldInfo, ierr )
    ELSE
       obj%finiteElectricFieldInfo_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "dipoleInfo")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputElectricFieldType","dipoleInfo: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputElectricFieldType","dipoleInfo: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%dipoleInfo_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_dipoleOutput(tmp_node, obj%dipoleInfo, ierr )
    ELSE
       obj%dipoleInfo_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "gateInfo")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputElectricFieldType","gateInfo: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputElectricFieldType","gateInfo: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%gateInfo_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_gateInfo(tmp_node, obj%gateInfo, ierr )
    ELSE
       obj%gateInfo_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_outputElectricField
  !
  !
  SUBROUTINE qes_read_BerryPhaseOutput(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(BerryPhaseOutput_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "totalPolarization")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:BerryPhaseOutputType","totalPolarization: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:BerryPhaseOutputType","totalPolarization: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_polarization(tmp_node, obj%totalPolarization, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "totalPhase")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:BerryPhaseOutputType","totalPhase: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:BerryPhaseOutputType","totalPhase: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_phase(tmp_node, obj%totalPhase, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "ionicPolarization")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size < 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:BerryPhaseOutputType","ionicPolarization: not enough elements")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:BerryPhaseOutputType","ionicPolarization: not enough elements",10)
        END IF
    END IF
    !
    obj%ndim_ionicPolarization = tmp_node_list_size
    ALLOCATE(obj%ionicPolarization(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_ionicPolarization(tmp_node, obj%ionicPolarization(index), ierr )
    END DO
    !
    tmp_node_list => getElementsByTagname(xml_node, "electronicPolarization")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size < 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:BerryPhaseOutputType","electronicPolarization: not enough elements")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:BerryPhaseOutputType","electronicPolarization: not enough elements",10)
        END IF
    END IF
    !
    obj%ndim_electronicPolarization = tmp_node_list_size
    ALLOCATE(obj%electronicPolarization(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_electronicPolarization(tmp_node, obj%electronicPolarization(index), ierr )
    END DO
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_BerryPhaseOutput
  !
  !
  SUBROUTINE qes_read_dipoleOutput(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(dipoleOutput_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "idir")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:dipoleOutputType","idir: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:dipoleOutputType","idir: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%idir, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:dipoleOutputType","error reading idir")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:dipoleOutputType","error reading idir",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "dipole")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:dipoleOutputType","dipole: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:dipoleOutputType","dipole: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_scalarQuantity(tmp_node, obj%dipole, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "ion_dipole")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:dipoleOutputType","ion_dipole: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:dipoleOutputType","ion_dipole: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_scalarQuantity(tmp_node, obj%ion_dipole, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "elec_dipole")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:dipoleOutputType","elec_dipole: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:dipoleOutputType","elec_dipole: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_scalarQuantity(tmp_node, obj%elec_dipole, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "dipoleField")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:dipoleOutputType","dipoleField: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:dipoleOutputType","dipoleField: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_scalarQuantity(tmp_node, obj%dipoleField, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "potentialAmp")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:dipoleOutputType","potentialAmp: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:dipoleOutputType","potentialAmp: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_scalarQuantity(tmp_node, obj%potentialAmp, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "totalLength")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:dipoleOutputType","totalLength: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:dipoleOutputType","totalLength: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_scalarQuantity(tmp_node, obj%totalLength, ierr )
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_dipoleOutput
  !
  !
  SUBROUTINE qes_read_finiteFieldOut(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(finiteFieldOut_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "electronicDipole")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:finiteFieldOutType","electronicDipole: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:finiteFieldOutType","electronicDipole: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%electronicDipole, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:finiteFieldOutType","error reading electronicDipole")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:finiteFieldOutType","error reading electronicDipole",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ionicDipole")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:finiteFieldOutType","ionicDipole: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:finiteFieldOutType","ionicDipole: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%ionicDipole, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:finiteFieldOutType","error reading ionicDipole")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:finiteFieldOutType","error reading ionicDipole",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_finiteFieldOut
  !
  !
  SUBROUTINE qes_read_polarization(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(polarization_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "polarization")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:polarizationType","polarization: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:polarizationType","polarization: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_scalarQuantity(tmp_node, obj%polarization, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "modulus")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:polarizationType","modulus: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:polarizationType","modulus: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%modulus, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:polarizationType","error reading modulus")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:polarizationType","error reading modulus",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "direction")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:polarizationType","direction: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:polarizationType","direction: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%direction, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:polarizationType","error reading direction")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:polarizationType","error reading direction",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_polarization
  !
  !
  SUBROUTINE qes_read_ionicPolarization(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(ionicPolarization_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "ion")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ionicPolarizationType","ion: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ionicPolarizationType","ion: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_atom(tmp_node, obj%ion, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "charge")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ionicPolarizationType","charge: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ionicPolarizationType","charge: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%charge, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:ionicPolarizationType","error reading charge")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:ionicPolarizationType","error reading charge",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "phase")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ionicPolarizationType","phase: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ionicPolarizationType","phase: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_phase(tmp_node, obj%phase, ierr )
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_ionicPolarization
  !
  !
  SUBROUTINE qes_read_electronicPolarization(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(electronicPolarization_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "firstKeyPoint")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electronicPolarizationType","firstKeyPoint: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electronicPolarizationType","firstKeyPoint: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_k_point(tmp_node, obj%firstKeyPoint, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "spin")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electronicPolarizationType","spin: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electronicPolarizationType","spin: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%spin_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%spin , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:electronicPolarizationType","error reading spin")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:electronicPolarizationType","error reading spin",10)
         END IF
      END IF
    ELSE
       obj%spin_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "phase")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:electronicPolarizationType","phase: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:electronicPolarizationType","phase: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_phase(tmp_node, obj%phase, ierr )
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_electronicPolarization
  !
  !
  SUBROUTINE qes_read_phase(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(phase_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "ionic")) THEN
      CALL extractDataAttribute(xml_node, "ionic", obj%ionic)
      obj%ionic_ispresent = .TRUE.
    ELSE
      obj%ionic_ispresent = .FALSE.
    END IF
    !
    IF (hasAttribute(xml_node, "electronic")) THEN
      CALL extractDataAttribute(xml_node, "electronic", obj%electronic)
      obj%electronic_ispresent = .TRUE.
    ELSE
      obj%electronic_ispresent = .FALSE.
    END IF
    !
    IF (hasAttribute(xml_node, "modulus")) THEN
      CALL extractDataAttribute(xml_node, "modulus", obj%modulus)
      obj%modulus_ispresent = .TRUE.
    ELSE
      obj%modulus_ispresent = .FALSE.
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%phase )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_phase
  !
  !
  SUBROUTINE qes_read_gateInfo(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(gateInfo_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "pot_prefactor")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gateInfoType","pot_prefactor: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gateInfoType","pot_prefactor: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%pot_prefactor, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:gateInfoType","error reading pot_prefactor")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:gateInfoType","error reading pot_prefactor",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "gate_zpos")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gateInfoType","gate_zpos: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gateInfoType","gate_zpos: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%gate_zpos, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:gateInfoType","error reading gate_zpos")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:gateInfoType","error reading gate_zpos",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "gate_gate_term")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gateInfoType","gate_gate_term: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gateInfoType","gate_gate_term: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%gate_gate_term, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:gateInfoType","error reading gate_gate_term")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:gateInfoType","error reading gate_gate_term",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "gatefieldEnergy")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:gateInfoType","gatefieldEnergy: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:gateInfoType","gatefieldEnergy: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%gatefieldEnergy, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:gateInfoType","error reading gatefieldEnergy")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:gateInfoType","error reading gatefieldEnergy",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_gateInfo
  !
  !
  SUBROUTINE qes_read_convergence_info(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(convergence_info_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "scf_conv")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:convergence_infoType","scf_conv: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:convergence_infoType","scf_conv: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_scf_conv(tmp_node, obj%scf_conv, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "opt_conv")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:convergence_infoType","opt_conv: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:convergence_infoType","opt_conv: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%opt_conv_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_opt_conv(tmp_node, obj%opt_conv, ierr )
    ELSE
       obj%opt_conv_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_convergence_info
  !
  !
  SUBROUTINE qes_read_scf_conv(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(scf_conv_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "n_scf_steps")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:scf_convType","n_scf_steps: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:scf_convType","n_scf_steps: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%n_scf_steps, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:scf_convType","error reading n_scf_steps")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:scf_convType","error reading n_scf_steps",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "scf_error")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:scf_convType","scf_error: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:scf_convType","scf_error: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%scf_error, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:scf_convType","error reading scf_error")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:scf_convType","error reading scf_error",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_scf_conv
  !
  !
  SUBROUTINE qes_read_opt_conv(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(opt_conv_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "n_opt_steps")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:opt_convType","n_opt_steps: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:opt_convType","n_opt_steps: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%n_opt_steps, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:opt_convType","error reading n_opt_steps")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:opt_convType","error reading n_opt_steps",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "grad_norm")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:opt_convType","grad_norm: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:opt_convType","grad_norm: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%grad_norm, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:opt_convType","error reading grad_norm")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:opt_convType","error reading grad_norm",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_opt_conv
  !
  !
  SUBROUTINE qes_read_algorithmic_info(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(algorithmic_info_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "real_space_q")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:algorithmic_infoType","real_space_q: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:algorithmic_infoType","real_space_q: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%real_space_q, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:algorithmic_infoType","error reading real_space_q")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:algorithmic_infoType","error reading real_space_q",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "uspp")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:algorithmic_infoType","uspp: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:algorithmic_infoType","uspp: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%uspp, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:algorithmic_infoType","error reading uspp")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:algorithmic_infoType","error reading uspp",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "paw")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:algorithmic_infoType","paw: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:algorithmic_infoType","paw: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%paw, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:algorithmic_infoType","error reading paw")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:algorithmic_infoType","error reading paw",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_algorithmic_info
  !
  !
  SUBROUTINE qes_read_symmetries(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(symmetries_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "nsym")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetriesType","nsym: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetriesType","nsym: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nsym, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:symmetriesType","error reading nsym")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:symmetriesType","error reading nsym",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nrot")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetriesType","nrot: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetriesType","nrot: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nrot, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:symmetriesType","error reading nrot")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:symmetriesType","error reading nrot",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "space_group")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetriesType","space_group: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetriesType","space_group: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%space_group, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:symmetriesType","error reading space_group")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:symmetriesType","error reading space_group",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "symmetry")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size < 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetriesType","symmetry: not enough elements")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetriesType","symmetry: not enough elements",10)
        END IF
    END IF
    IF (tmp_node_list_size > 48) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetriesType","symmetry: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetriesType","symmetry: too many occurrences",10)
        END IF
    END IF
    !
    obj%ndim_symmetry = tmp_node_list_size
    ALLOCATE(obj%symmetry(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_symmetry(tmp_node, obj%symmetry(index), ierr )
    END DO
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_symmetries
  !
  !
  SUBROUTINE qes_read_symmetry(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(symmetry_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "info")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetryType","info: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetryType","info: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_info(tmp_node, obj%info, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "rotation")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetryType","rotation: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetryType","rotation: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_matrix(tmp_node, obj%rotation, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "fractional_translation")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetryType","fractional_translation: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetryType","fractional_translation: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%fractional_translation_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%fractional_translation , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:symmetryType","error reading fractional_translation")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:symmetryType","error reading fractional_translation",10)
         END IF
      END IF
    ELSE
       obj%fractional_translation_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "equivalent_atoms")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:symmetryType","equivalent_atoms: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:symmetryType","equivalent_atoms: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%equivalent_atoms_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_equivalent_atoms(tmp_node, obj%equivalent_atoms, ierr )
    ELSE
       obj%equivalent_atoms_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_symmetry
  !
  !
  SUBROUTINE qes_read_equivalent_atoms(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(equivalent_atoms_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "nat")) THEN
      CALL extractDataAttribute(xml_node, "nat", obj%nat)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: equivalent_atomsType",&
                        "required attribute nat not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: equivalent_atomsType",&
                      "required attribute nat not found", 10 ) 
      END IF
    END IF
    !
    !
    IF (hasAttribute(xml_node, "size"))  THEN
        CALL extractDataAttribute(xml_node, "size", obj%size)
    ELSE
        CALL errore ("qes_read: equivalent_atomsType", &
                     "mandatory size attribute not found in "//TRIM(obj%tagname), 12)
    END IF
    !
    !
    !
    ALLOCATE (obj%equivalent_atoms(obj%size))
    CALL extractDataContent(xml_node, obj%equivalent_atoms)
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_equivalent_atoms
  !
  !
  SUBROUTINE qes_read_info(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(info_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "name")) THEN
      CALL extractDataAttribute(xml_node, "name", obj%name)
      obj%name_ispresent = .TRUE.
    ELSE
      obj%name_ispresent = .FALSE.
    END IF
    !
    IF (hasAttribute(xml_node, "class")) THEN
      CALL extractDataAttribute(xml_node, "class", obj%class)
      obj%class_ispresent = .TRUE.
    ELSE
      obj%class_ispresent = .FALSE.
    END IF
    !
    IF (hasAttribute(xml_node, "time_reversal")) THEN
      CALL extractDataAttribute(xml_node, "time_reversal", obj%time_reversal)
      obj%time_reversal_ispresent = .TRUE.
    ELSE
      obj%time_reversal_ispresent = .FALSE.
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%info )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_info
  !
  !
  SUBROUTINE qes_read_outputPBC(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(outputPBC_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "assume_isolated")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:outputPBCType","assume_isolated: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:outputPBCType","assume_isolated: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%assume_isolated, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:outputPBCType","error reading assume_isolated")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:outputPBCType","error reading assume_isolated",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_outputPBC
  !
  !
  SUBROUTINE qes_read_magnetization(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(magnetization_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "lsda")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:magnetizationType","lsda: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:magnetizationType","lsda: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%lsda, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:magnetizationType","error reading lsda")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:magnetizationType","error reading lsda",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "noncolin")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:magnetizationType","noncolin: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:magnetizationType","noncolin: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%noncolin, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:magnetizationType","error reading noncolin")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:magnetizationType","error reading noncolin",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "spinorbit")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:magnetizationType","spinorbit: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:magnetizationType","spinorbit: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%spinorbit, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:magnetizationType","error reading spinorbit")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:magnetizationType","error reading spinorbit",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "total")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:magnetizationType","total: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:magnetizationType","total: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%total, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:magnetizationType","error reading total")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:magnetizationType","error reading total",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "absolute")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:magnetizationType","absolute: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:magnetizationType","absolute: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%absolute, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:magnetizationType","error reading absolute")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:magnetizationType","error reading absolute",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "do_magnetization")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:magnetizationType","do_magnetization: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:magnetizationType","do_magnetization: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%do_magnetization, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:magnetizationType","error reading do_magnetization")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:magnetizationType","error reading do_magnetization",10)
       END IF
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_magnetization
  !
  !
  SUBROUTINE qes_read_total_energy(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(total_energy_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "etot")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:total_energyType","etot: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:total_energyType","etot: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%etot, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:total_energyType","error reading etot")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:total_energyType","error reading etot",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "eband")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:total_energyType","eband: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:total_energyType","eband: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%eband_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%eband , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:total_energyType","error reading eband")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:total_energyType","error reading eband",10)
         END IF
      END IF
    ELSE
       obj%eband_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ehart")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:total_energyType","ehart: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:total_energyType","ehart: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%ehart_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%ehart , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:total_energyType","error reading ehart")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:total_energyType","error reading ehart",10)
         END IF
      END IF
    ELSE
       obj%ehart_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "vtxc")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:total_energyType","vtxc: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:total_energyType","vtxc: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%vtxc_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%vtxc , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:total_energyType","error reading vtxc")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:total_energyType","error reading vtxc",10)
         END IF
      END IF
    ELSE
       obj%vtxc_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "etxc")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:total_energyType","etxc: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:total_energyType","etxc: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%etxc_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%etxc , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:total_energyType","error reading etxc")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:total_energyType","error reading etxc",10)
         END IF
      END IF
    ELSE
       obj%etxc_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ewald")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:total_energyType","ewald: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:total_energyType","ewald: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%ewald_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%ewald , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:total_energyType","error reading ewald")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:total_energyType","error reading ewald",10)
         END IF
      END IF
    ELSE
       obj%ewald_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "demet")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:total_energyType","demet: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:total_energyType","demet: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%demet_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%demet , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:total_energyType","error reading demet")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:total_energyType","error reading demet",10)
         END IF
      END IF
    ELSE
       obj%demet_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "efieldcorr")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:total_energyType","efieldcorr: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:total_energyType","efieldcorr: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%efieldcorr_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%efieldcorr , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:total_energyType","error reading efieldcorr")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:total_energyType","error reading efieldcorr",10)
         END IF
      END IF
    ELSE
       obj%efieldcorr_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "potentiostat_contr")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:total_energyType","potentiostat_contr: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:total_energyType","potentiostat_contr: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%potentiostat_contr_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%potentiostat_contr , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:total_energyType","error reading potentiostat_contr")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:total_energyType","error reading potentiostat_contr",10)
         END IF
      END IF
    ELSE
       obj%potentiostat_contr_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "gatefield_contr")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:total_energyType","gatefield_contr: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:total_energyType","gatefield_contr: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%gatefield_contr_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%gatefield_contr , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:total_energyType","error reading gatefield_contr")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:total_energyType","error reading gatefield_contr",10)
         END IF
      END IF
    ELSE
       obj%gatefield_contr_ispresent = .FALSE.
    END IF
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_total_energy
  !
  !
  SUBROUTINE qes_read_band_structure(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(band_structure_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "lsda")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","lsda: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","lsda: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%lsda, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:band_structureType","error reading lsda")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:band_structureType","error reading lsda",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "noncolin")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","noncolin: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","noncolin: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%noncolin, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:band_structureType","error reading noncolin")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:band_structureType","error reading noncolin",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "spinorbit")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","spinorbit: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","spinorbit: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%spinorbit, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:band_structureType","error reading spinorbit")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:band_structureType","error reading spinorbit",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nbnd")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","nbnd: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","nbnd: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nbnd, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:band_structureType","error reading nbnd")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:band_structureType","error reading nbnd",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nbnd_up")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","nbnd_up: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","nbnd_up: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%nbnd_up_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%nbnd_up , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:band_structureType","error reading nbnd_up")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:band_structureType","error reading nbnd_up",10)
         END IF
      END IF
    ELSE
       obj%nbnd_up_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nbnd_dw")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","nbnd_dw: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","nbnd_dw: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%nbnd_dw_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%nbnd_dw , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:band_structureType","error reading nbnd_dw")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:band_structureType","error reading nbnd_dw",10)
         END IF
      END IF
    ELSE
       obj%nbnd_dw_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "nelec")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","nelec: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","nelec: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nelec, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:band_structureType","error reading nelec")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:band_structureType","error reading nelec",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "num_of_atomic_wfc")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","num_of_atomic_wfc: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","num_of_atomic_wfc: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%num_of_atomic_wfc_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%num_of_atomic_wfc , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:band_structureType","error reading num_of_atomic_wfc")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:band_structureType","error reading num_of_atomic_wfc",10)
         END IF
      END IF
    ELSE
       obj%num_of_atomic_wfc_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "wf_collected")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","wf_collected: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","wf_collected: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%wf_collected, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:band_structureType","error reading wf_collected")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:band_structureType","error reading wf_collected",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "fermi_energy")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","fermi_energy: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","fermi_energy: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%fermi_energy_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%fermi_energy , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:band_structureType","error reading fermi_energy")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:band_structureType","error reading fermi_energy",10)
         END IF
      END IF
    ELSE
       obj%fermi_energy_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "highestOccupiedLevel")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","highestOccupiedLevel: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","highestOccupiedLevel: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%highestOccupiedLevel_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%highestOccupiedLevel , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:band_structureType","error reading highestOccupiedLevel")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:band_structureType","error reading highestOccupiedLevel",10)
         END IF
      END IF
    ELSE
       obj%highestOccupiedLevel_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "two_fermi_energies")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","two_fermi_energies: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","two_fermi_energies: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%two_fermi_energies_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL extractDataContent(tmp_node, obj%two_fermi_energies , IOSTAT = iostat_)
      IF ( iostat_ /= 0 ) THEN
         IF ( PRESENT (ierr ) ) THEN 
            CALL infomsg("qes_read:band_structureType","error reading two_fermi_energies")
            ierr = ierr + 1
         ELSE 
            CALL errore ("qes_read:band_structureType","error reading two_fermi_energies",10)
         END IF
      END IF
    ELSE
       obj%two_fermi_energies_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "starting_k_points")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","starting_k_points: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","starting_k_points: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_k_points_IBZ(tmp_node, obj%starting_k_points, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "nks")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","nks: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","nks: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%nks, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:band_structureType","error reading nks")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:band_structureType","error reading nks",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "occupations_kind")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","occupations_kind: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","occupations_kind: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_occupations(tmp_node, obj%occupations_kind, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "smearing")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size > 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","smearing: too many occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","smearing: too many occurrences",10)
        END IF
    END IF
    !
    IF (tmp_node_list_size>0) THEN
      obj%smearing_ispresent = .TRUE.
      tmp_node => item(tmp_node_list, 0)
      CALL qes_read_smearing(tmp_node, obj%smearing, ierr )
    ELSE
       obj%smearing_ispresent = .FALSE.
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "ks_energies")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size < 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:band_structureType","ks_energies: not enough elements")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:band_structureType","ks_energies: not enough elements",10)
        END IF
    END IF
    !
    obj%ndim_ks_energies = tmp_node_list_size
    ALLOCATE(obj%ks_energies(tmp_node_list_size))
    DO index=1,tmp_node_list_size
        tmp_node => item( tmp_node_list, index-1 )
        CALL qes_read_ks_energies(tmp_node, obj%ks_energies(index), ierr )
    END DO
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_band_structure
  !
  !
  SUBROUTINE qes_read_ks_energies(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(ks_energies_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    !
    !
    tmp_node_list => getElementsByTagname(xml_node, "k_point")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ks_energiesType","k_point: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ks_energiesType","k_point: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_k_point(tmp_node, obj%k_point, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "npw")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ks_energiesType","npw: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ks_energiesType","npw: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL extractDataContent(tmp_node, obj%npw, IOSTAT = iostat_ )
    IF ( iostat_ /= 0 ) THEN
       IF ( PRESENT (ierr ) ) THEN 
          CALL infomsg("qes_read:ks_energiesType","error reading npw")
          ierr = ierr + 1
       ELSE 
          CALL errore ("qes_read:ks_energiesType","error reading npw",10)
       END IF
    END IF
    !
    tmp_node_list => getElementsByTagname(xml_node, "eigenvalues")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ks_energiesType","eigenvalues: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ks_energiesType","eigenvalues: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_vector(tmp_node, obj%eigenvalues, ierr )
    !
    tmp_node_list => getElementsByTagname(xml_node, "occupations")
    tmp_node_list_size = getLength(tmp_node_list)
    !
    IF (tmp_node_list_size /= 1) THEN
        IF (PRESENT(ierr) ) THEN 
           CALL infomsg("qes_read:ks_energiesType","occupations: wrong number of occurrences")
           ierr = ierr + 1 
        ELSE 
           CALL errore("qes_read:ks_energiesType","occupations: wrong number of occurrences",10)
        END IF
    END IF
    !
    tmp_node => item(tmp_node_list, 0)
    IF (ASSOCIATED(tmp_node))&
       CALL qes_read_vector(tmp_node, obj%occupations, ierr )
    !
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_ks_energies
  !
  !
  SUBROUTINE qes_read_closed(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(closed_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "DATE")) THEN
      CALL extractDataAttribute(xml_node, "DATE", obj%DATE)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: closedType",&
                        "required attribute DATE not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: closedType",&
                      "required attribute DATE not found", 10 ) 
      END IF
    END IF
    !
    IF (hasAttribute(xml_node, "TIME")) THEN
      CALL extractDataAttribute(xml_node, "TIME", obj%TIME)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: closedType",&
                        "required attribute TIME not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: closedType",&
                      "required attribute TIME not found", 10 ) 
      END IF
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%closed )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_closed
  !
  !
  SUBROUTINE qes_read_vector(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(vector_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "size"))  THEN
        CALL extractDataAttribute(xml_node, "size", obj%size)
    ELSE
        CALL errore ("qes_read: vectorType", &
                     "mandatory size attribute not found in "//TRIM(obj%tagname), 12)
    END IF
    !
    !
    !
    !
    ALLOCATE (obj%vector(obj%size))
    CALL extractDataContent(xml_node, obj%vector )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_vector
  !
  !
  SUBROUTINE qes_read_integerVector(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(integerVector_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "size"))  THEN
        CALL extractDataAttribute(xml_node, "size", obj%size)
    ELSE
        CALL errore ("qes_read: integerVectorType", &
                     "mandatory size attribute not found in "//TRIM(obj%tagname), 12)
    END IF
    !
    !
    !
    !
    ALLOCATE (obj%integerVector(obj%size))
    CALL extractDataContent(xml_node, obj%integerVector)
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_integerVector
  !
  !
  SUBROUTINE qes_read_matrix(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(matrix_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    INTEGER :: i, length
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "rank"))  THEN
        CALL extractDataAttribute(xml_node, "rank", obj%rank)
    ELSE
        CALL errore ("qes_read: matrixType",&
                      "required attribute rank not found, can't read further, stopping", 10 ) 
    END IF
    ALLOCATE (obj%dims(obj%rank))
    IF (hasAttribute(xml_node, "dims")) THEN
        CALL extractDataAttribute(xml_node, "dims", obj%dims)
    ELSE
        CALL errore ("qes_read: matrixType",&
                      "required attribute dims not found, can't read further, stopping", 10 ) 
    END IF
    IF (hasAttribute(xml_node,"order")) THEN
        CALL extractDataAttribute(xml_node, "order", obj%order)
    ELSE
        obj%order = "F"
    END IF
    !
    !
    !
    length = 1
    DO i =1, obj%rank
       length = length * obj%dims(i)
    END DO
    ALLOCATE (obj%matrix(length) )
    CALL extractDataContent(xml_node, obj%matrix )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_matrix
  !
  !
  SUBROUTINE qes_read_integerMatrix(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(integerMatrix_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    INTEGER :: i, length
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "rank"))  THEN
        CALL extractDataAttribute(xml_node, "rank", obj%rank)
    ELSE
        CALL errore ("qes_read: integerMatrixType",&
                      "required attribute rank not found, can't read further, stopping", 10 ) 
    END IF
    ALLOCATE (obj%dims(obj%rank))
    IF (hasAttribute(xml_node, "dims")) THEN
        CALL extractDataAttribute(xml_node, "dims", obj%dims)
    ELSE
        CALL errore ("qes_read: integerMatrixType",&
                      "required attribute dims not found, can't read further, stopping", 10 ) 
    END IF
    IF (hasAttribute(xml_node,"order")) THEN
        CALL extractDataAttribute(xml_node, "order", obj%order)
    ELSE
        obj%order = "F"
    END IF
    !
    !
    !
    length = 1
    DO i = 1, obj%rank
        length = length * obj%dims(i)
    END DO
    ALLOCATE(obj%integerMatrix(length))
    CALL extractDataContent(xml_node, obj%integerMatrix)
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_integerMatrix
  !
  !
  SUBROUTINE qes_read_scalarQuantity(xml_node, obj, ierr )
    !
    IMPLICIT NONE
    !
    TYPE(Node), INTENT(IN), POINTER                 :: xml_node
    TYPE(scalarQuantity_type), INTENT(OUT) :: obj
    INTEGER, OPTIONAL, INTENT(OUT)                  :: ierr 
    ! 
    TYPE(Node), POINTER :: tmp_node
    TYPE(NodeList), POINTER :: tmp_node_list
    INTEGER :: tmp_node_list_size, index, iostat_
    !
    obj%tagname = getTagName(xml_node)
    !
    IF (hasAttribute(xml_node, "Units")) THEN
      CALL extractDataAttribute(xml_node, "Units", obj%Units)
    ELSE
      IF ( PRESENT(ierr) ) THEN 
         CALL infomsg ( "qes_read: scalarQuantityType",&
                        "required attribute Units not found" )
         ierr = ierr + 1 
      ELSE 
         CALL errore ("qes_read: scalarQuantityType",&
                      "required attribute Units not found", 10 ) 
      END IF
    END IF
    !
    !
    !
    !
    !
    CALL extractDataContent(xml_node, obj%scalarQuantity )
    !
    obj%lwrite = .TRUE.
    !
  END SUBROUTINE qes_read_scalarQuantity
  !
  !
END MODULE qes_read_module
