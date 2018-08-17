!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE bcast_qes_types_module
  !
  ! Auto-generated code: don't edit this file
  !
  ! Quantum Espresso XSD namespace: http://www.quantum-espresso.org/ns/qes/qes-1.0 
  !
  USE kinds, ONLY: DP
  USE io_global, ONLY : ionode
  USE mp,    ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  PRIVATE 
  PUBLIC qes_bcast
  INTERFACE qes_bcast 
    MODULE PROCEDURE  bcast_xml_format_type
    MODULE PROCEDURE  bcast_creator_type
    MODULE PROCEDURE  bcast_created_type
    MODULE PROCEDURE  bcast_atom_type
    MODULE PROCEDURE  bcast_qpoint_grid_type
    MODULE PROCEDURE  bcast_HubbardCommon_type
    MODULE PROCEDURE  bcast_HubbardJ_type
    MODULE PROCEDURE  bcast_starting_ns_type
    MODULE PROCEDURE  bcast_hubbard_ns_type
    MODULE PROCEDURE  bcast_smearing_type
    MODULE PROCEDURE  bcast_occupations_type
    MODULE PROCEDURE  bcast_basissetitem_type
    MODULE PROCEDURE  bcast_monkhorst_pack_type
    MODULE PROCEDURE  bcast_k_point_type
    MODULE PROCEDURE  bcast_inputOccupations_type
    MODULE PROCEDURE  bcast_phase_type
    MODULE PROCEDURE  bcast_equivalent_atoms_type
    MODULE PROCEDURE  bcast_info_type
    MODULE PROCEDURE  bcast_closed_type
    MODULE PROCEDURE  bcast_vector_type
    MODULE PROCEDURE  bcast_integervector_type
    MODULE PROCEDURE  bcast_matrix_type
    MODULE PROCEDURE  bcast_integermatrix_type
    MODULE PROCEDURE  bcast_scalarquantity_type
    MODULE PROCEDURE  bcast_general_info_type
    MODULE PROCEDURE  bcast_parallel_info_type
    MODULE PROCEDURE  bcast_control_variables_type
    MODULE PROCEDURE  bcast_species_type
    MODULE PROCEDURE  bcast_atomic_positions_type
    MODULE PROCEDURE  bcast_wyckoff_positions_type
    MODULE PROCEDURE  bcast_cell_type
    MODULE PROCEDURE  bcast_hybrid_type
    MODULE PROCEDURE  bcast_dftu_type
    MODULE PROCEDURE  bcast_vdw_type
    MODULE PROCEDURE  bcast_spin_type
    MODULE PROCEDURE  bcast_bands_type
    MODULE PROCEDURE  bcast_basis_type
    MODULE PROCEDURE  bcast_reciprocal_lattice_type
    MODULE PROCEDURE  bcast_electron_control_type
    MODULE PROCEDURE  bcast_k_points_IBZ_type
    MODULE PROCEDURE  bcast_bfgs_type
    MODULE PROCEDURE  bcast_md_type
    MODULE PROCEDURE  bcast_cell_control_type
    MODULE PROCEDURE  bcast_symmetry_flags_type
    MODULE PROCEDURE  bcast_esm_type
    MODULE PROCEDURE  bcast_ekin_functional_type
    MODULE PROCEDURE  bcast_spin_constraints_type
    MODULE PROCEDURE  bcast_electric_field_type
    MODULE PROCEDURE  bcast_atomic_constraint_type
    MODULE PROCEDURE  bcast_dipoleOutput_type
    MODULE PROCEDURE  bcast_finiteFieldOut_type
    MODULE PROCEDURE  bcast_polarization_type
    MODULE PROCEDURE  bcast_ionicPolarization_type
    MODULE PROCEDURE  bcast_electronicPolarization_type
    MODULE PROCEDURE  bcast_scf_conv_type
    MODULE PROCEDURE  bcast_opt_conv_type
    MODULE PROCEDURE  bcast_algorithmic_info_type
    MODULE PROCEDURE  bcast_symmetry_type
    MODULE PROCEDURE  bcast_magnetization_type
    MODULE PROCEDURE  bcast_total_energy_type
    MODULE PROCEDURE  bcast_ks_energies_type
    MODULE PROCEDURE  bcast_atomic_species_type
    MODULE PROCEDURE  bcast_atomic_structure_type
    MODULE PROCEDURE  bcast_dft_type
    MODULE PROCEDURE  bcast_basis_set_type
    MODULE PROCEDURE  bcast_ion_control_type
    MODULE PROCEDURE  bcast_boundary_conditions_type
    MODULE PROCEDURE  bcast_atomic_constraints_type
    MODULE PROCEDURE  bcast_BerryPhaseOutput_type
    MODULE PROCEDURE  bcast_convergence_info_type
    MODULE PROCEDURE  bcast_symmetries_type
    MODULE PROCEDURE  bcast_band_structure_type
    MODULE PROCEDURE  bcast_input_type
    MODULE PROCEDURE  bcast_step_type
    MODULE PROCEDURE  bcast_outputElectricField_type
    MODULE PROCEDURE  bcast_output_type
END INTERFACE 

CONTAINS

SUBROUTINE bcast_xml_format_type(o, ionode_id, comm)

USE qes_types_module, ONLY : xml_format_type
IMPLICIT NONE
TYPE(xml_format_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%name, ionode_id, comm)
CALL mp_bcast(o%version, ionode_id, comm)
CALL mp_bcast(o%xml_format, ionode_id, comm)

RETURN
END SUBROUTINE bcast_xml_format_type

SUBROUTINE bcast_creator_type(o, ionode_id, comm)

USE qes_types_module, ONLY : creator_type
IMPLICIT NONE
TYPE(creator_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%name, ionode_id, comm)
CALL mp_bcast(o%version, ionode_id, comm)
CALL mp_bcast(o%creator, ionode_id, comm)

RETURN
END SUBROUTINE bcast_creator_type

SUBROUTINE bcast_created_type(o, ionode_id, comm)

USE qes_types_module, ONLY : created_type
IMPLICIT NONE
TYPE(created_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%date, ionode_id, comm)
CALL mp_bcast(o%time, ionode_id, comm)
CALL mp_bcast(o%created, ionode_id, comm)

RETURN
END SUBROUTINE bcast_created_type

SUBROUTINE bcast_atom_type(o, ionode_id, comm)

USE qes_types_module, ONLY : atom_type
IMPLICIT NONE
TYPE(atom_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%name, ionode_id, comm)
CALL mp_bcast(o%position_ispresent, ionode_id, comm)
IF (o%position_ispresent) CALL mp_bcast(o%position, ionode_id, comm)
CALL mp_bcast(o%index_ispresent, ionode_id, comm)
IF (o%index_ispresent) CALL mp_bcast(o%index, ionode_id, comm)
CALL mp_bcast(o%atom, ionode_id, comm)

RETURN
END SUBROUTINE bcast_atom_type

SUBROUTINE bcast_qpoint_grid_type(o, ionode_id, comm)

USE qes_types_module, ONLY : qpoint_grid_type
IMPLICIT NONE
TYPE(qpoint_grid_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%nqx1, ionode_id, comm)
CALL mp_bcast(o%nqx2, ionode_id, comm)
CALL mp_bcast(o%nqx3, ionode_id, comm)
CALL mp_bcast(o%qpoint_grid, ionode_id, comm)
    !
RETURN
END SUBROUTINE bcast_qpoint_grid_type

SUBROUTINE bcast_HubbardCommon_type(o, ionode_id, comm)

USE qes_types_module, ONLY : HubbardCommon_type
IMPLICIT NONE
TYPE(HubbardCommon_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%specie, ionode_id, comm)
CALL mp_bcast(o%label, ionode_id, comm)
CALL mp_bcast(o%HubbardCommon, ionode_id, comm)

RETURN
END SUBROUTINE bcast_HubbardCommon_type
    !
SUBROUTINE bcast_HubbardJ_type(o, ionode_id, comm)

USE qes_types_module, ONLY : HubbardJ_type
IMPLICIT NONE
TYPE(HubbardJ_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%specie, ionode_id, comm)
CALL mp_bcast(o%label, ionode_id, comm)
CALL mp_bcast(o%HubbardJ, ionode_id, comm)

RETURN
END SUBROUTINE bcast_HubbardJ_type

SUBROUTINE bcast_starting_ns_type(o, ionode_id, comm)

USE qes_types_module, ONLY : starting_ns_type
IMPLICIT NONE
TYPE(starting_ns_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%specie, ionode_id, comm)
CALL mp_bcast(o%label, ionode_id, comm)
CALL mp_bcast(o%spin, ionode_id, comm)
CALL mp_bcast(o%size, ionode_id, comm)
IF (.NOT.ionode) ALLOCATE(o%starting_ns(o%size))
CALL mp_bcast(o%starting_ns, ionode_id, comm)

RETURN
END SUBROUTINE bcast_starting_ns_type

SUBROUTINE bcast_hubbard_ns_type(o, ionode_id, comm)

USE qes_types_module, ONLY : hubbard_ns_type
IMPLICIT NONE
TYPE(hubbard_ns_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

INTEGER :: length, i

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%specie, ionode_id, comm)
CALL mp_bcast(o%label, ionode_id, comm)
CALL mp_bcast(o%spin, ionode_id, comm)
CALL mp_bcast(o%index, ionode_id, comm)
CALL mp_bcast(o%rank, ionode_id, comm)
IF (.NOT.ionode) ALLOCATE(o%dims(o%rank))
CALL mp_bcast(o%dims, ionode_id, comm)
CALL mp_bcast(o%order, ionode_id, comm)
IF (.NOT.ionode) THEN
   length = 1
   DO i =1, o%rank
      length = length * o%dims(i)
   END DO
   ALLOCATE (o%Hubbard_ns(length) )
ENDIF

CALL mp_bcast(o%Hubbard_ns, ionode_id, comm)
RETURN
END SUBROUTINE bcast_hubbard_ns_type

SUBROUTINE bcast_smearing_type(o, ionode_id, comm)

USE qes_types_module, ONLY : smearing_type
IMPLICIT NONE
TYPE(smearing_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%degauss, ionode_id, comm)
CALL mp_bcast(o%smearing, ionode_id, comm)

RETURN
END SUBROUTINE bcast_smearing_type

SUBROUTINE bcast_occupations_type(o, ionode_id, comm)

USE qes_types_module, ONLY : occupations_type
IMPLICIT NONE
TYPE(occupations_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%spin_ispresent, ionode_id, comm)
IF (o%spin_ispresent) &
   CALL mp_bcast(o%spin, ionode_id, comm)
CALL mp_bcast(o%occupations, ionode_id, comm)

RETURN
END SUBROUTINE bcast_occupations_type


SUBROUTINE bcast_basissetitem_type(o, ionode_id, comm)

USE qes_types_module, ONLY : basissetitem_type
IMPLICIT NONE
TYPE(basissetitem_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%nr1, ionode_id, comm)
CALL mp_bcast(o%nr2, ionode_id, comm)
CALL mp_bcast(o%nr3, ionode_id, comm)
CALL mp_bcast(o%basisSetItem, ionode_id, comm)

RETURN
END SUBROUTINE bcast_basissetitem_type

SUBROUTINE bcast_monkhorst_pack_type(o, ionode_id, comm)

USE qes_types_module, ONLY : monkhorst_pack_type
IMPLICIT NONE
TYPE(monkhorst_pack_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%nk1, ionode_id, comm)
CALL mp_bcast(o%nk2, ionode_id, comm)
CALL mp_bcast(o%nk3, ionode_id, comm)
CALL mp_bcast(o%k1, ionode_id, comm)
CALL mp_bcast(o%k2, ionode_id, comm)
CALL mp_bcast(o%k3, ionode_id, comm)
CALL mp_bcast(o%monkhorst_pack, ionode_id, comm)

RETURN
END SUBROUTINE bcast_monkhorst_pack_type
    !
SUBROUTINE bcast_k_point_type(o, ionode_id, comm)

USE qes_types_module, ONLY : k_point_type
IMPLICIT NONE
TYPE(k_point_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%weight_ispresent, ionode_id, comm)
IF (o%weight_ispresent) &
   CALL mp_bcast(o%weight, ionode_id, comm)
CALL mp_bcast(o%label_ispresent, ionode_id, comm)
IF (o%label_ispresent)&
   CALL mp_bcast(o%label, ionode_id, comm)
CALL mp_bcast(o%k_point, ionode_id, comm)

RETURN
END SUBROUTINE bcast_k_point_type

SUBROUTINE bcast_inputOccupations_type(o, ionode_id, comm)

USE qes_types_module, ONLY : inputOccupations_type
IMPLICIT NONE
TYPE(inputOccupations_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%ispin, ionode_id, comm)
CALL mp_bcast(o%spin_factor, ionode_id, comm)
CALL mp_bcast(o%size, ionode_id, comm)

IF (.NOT.ionode) ALLOCATE(o%inputOccupations(o%size))
CALL mp_bcast(o%inputOccupations, ionode_id, comm)

RETURN
END SUBROUTINE bcast_inputOccupations_type

SUBROUTINE bcast_phase_type(o, ionode_id, comm)

USE qes_types_module, ONLY : phase_type
IMPLICIT NONE
TYPE(phase_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%ionic_ispresent, ionode_id, comm)
IF (o%ionic_ispresent) &
   CALL mp_bcast(o%ionic, ionode_id, comm)
CALL mp_bcast(o%electronic_ispresent, ionode_id, comm)
IF (o%electronic_ispresent) &
   CALL mp_bcast(o%electronic, ionode_id, comm)
CALL mp_bcast(o%modulus_ispresent, ionode_id, comm)
IF (o%modulus_ispresent) &
   CALL mp_bcast(o%modulus, ionode_id, comm)

RETURN
END SUBROUTINE bcast_phase_type
    !
SUBROUTINE bcast_equivalent_atoms_type(o, ionode_id, comm)

USE qes_types_module, ONLY : equivalent_atoms_type
IMPLICIT NONE
TYPE(equivalent_atoms_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%nat, ionode_id, comm)
CALL mp_bcast(o%size, ionode_id, comm)
IF (.NOT.ionode) ALLOCATE(o%equivalent_atoms(o%size))
CALL mp_bcast(o%equivalent_atoms, ionode_id, comm)

RETURN
END SUBROUTINE bcast_equivalent_atoms_type

SUBROUTINE bcast_info_type(o, ionode_id, comm)

USE qes_types_module, ONLY : info_type
IMPLICIT NONE
TYPE(info_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%name_ispresent, ionode_id, comm)
IF (o%name_ispresent) &
   CALL mp_bcast(o%name, ionode_id, comm)
CALL mp_bcast(o%class_ispresent, ionode_id, comm)
IF (o%class_ispresent) &
   CALL mp_bcast(o%class, ionode_id, comm)
CALL mp_bcast(o%time_reversal_ispresent, ionode_id, comm)
IF (o%time_reversal_ispresent) &
   CALL mp_bcast(o%time_reversal, ionode_id, comm)

RETURN
END SUBROUTINE bcast_info_type

SUBROUTINE bcast_closed_type(o, ionode_id, comm)

USE qes_types_module, ONLY : closed_type
IMPLICIT NONE
TYPE(closed_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%date, ionode_id, comm)
CALL mp_bcast(o%time, ionode_id, comm)
CALL mp_bcast(o%closed, ionode_id, comm)

RETURN
END SUBROUTINE bcast_closed_type

SUBROUTINE bcast_vector_type(o, ionode_id, comm)

USE qes_types_module, ONLY : vector_type
IMPLICIT NONE
TYPE(vector_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%size, ionode_id, comm)
IF (.NOT.ionode) ALLOCATE(o%vector(o%size))
CALL mp_bcast(o%vector, ionode_id, comm)
RETURN
END SUBROUTINE bcast_vector_type

SUBROUTINE bcast_integervector_type(o, ionode_id, comm)

USE qes_types_module, ONLY : integervector_type
IMPLICIT NONE
TYPE(integervector_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%size, ionode_id, comm)
IF (.NOT.ionode) ALLOCATE(o%IntegerVector(o%size))
CALL mp_bcast(o%integerVector, ionode_id, comm)
RETURN
END SUBROUTINE bcast_integervector_type

SUBROUTINE bcast_matrix_type(o, ionode_id, comm)

USE qes_types_module, ONLY : matrix_type
IMPLICIT NONE
TYPE(matrix_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

INTEGER :: length, i

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%rank, ionode_id, comm)
IF (.NOT.ionode) ALLOCATE(o%dims(o%rank))
CALL mp_bcast(o%dims, ionode_id, comm)
CALL mp_bcast(o%order, ionode_id, comm)
IF (.NOT.ionode) THEN
   length = 1
   DO i =1, o%rank
      length = length * o%dims(i)
   END DO
   ALLOCATE(o%matrix(length))
ENDIF
CALL mp_bcast(o%matrix, ionode_id, comm)
RETURN
END SUBROUTINE bcast_matrix_type
    !
SUBROUTINE bcast_integermatrix_type(o, ionode_id, comm)

USE qes_types_module, ONLY : integermatrix_type
IMPLICIT NONE
TYPE(integermatrix_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

INTEGER :: length, i

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%rank, ionode_id, comm)
IF (.NOT.ionode) ALLOCATE(o%dims(o%rank))
CALL mp_bcast(o%dims, ionode_id, comm)
CALL mp_bcast(o%order, ionode_id, comm)
IF (.NOT.ionode) THEN
   length = 1
   DO i =1, o%rank
      length = length * o%dims(i)
   END DO
   ALLOCATE(o%integermatrix(length))
ENDIF

CALL mp_bcast(o%integerMatrix, ionode_id, comm)

RETURN
END SUBROUTINE bcast_integermatrix_type
    !
SUBROUTINE bcast_scalarquantity_type(o, ionode_id, comm)

USE qes_types_module, ONLY : scalarquantity_type
IMPLICIT NONE
TYPE(scalarquantity_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%units, ionode_id, comm)
CALL mp_bcast(o%scalarquantity, ionode_id, comm)

RETURN
END SUBROUTINE bcast_scalarquantity_type

SUBROUTINE bcast_general_info_type(o, ionode_id, comm)

USE qes_types_module, ONLY : general_info_type
IMPLICIT NONE
TYPE(general_info_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL bcast_xml_format_type(o%xml_format, ionode_id, comm)
CALL bcast_creator_type(o%creator, ionode_id, comm)
CALL bcast_created_type(o%created, ionode_id, comm)

RETURN
END SUBROUTINE bcast_general_info_type
    !

SUBROUTINE bcast_parallel_info_type(o, ionode_id, comm)

USE qes_types_module, ONLY : parallel_info_type
IMPLICIT NONE
TYPE(parallel_info_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%nprocs, ionode_id, comm)
CALL mp_bcast(o%nthreads, ionode_id, comm)
CALL mp_bcast(o%ntasks, ionode_id, comm)
CALL mp_bcast(o%nbgrp, ionode_id, comm)
CALL mp_bcast(o%npool, ionode_id, comm)
CALL mp_bcast(o%ndiag, ionode_id, comm)

RETURN
END SUBROUTINE bcast_parallel_info_type
    !
SUBROUTINE bcast_control_variables_type(o, ionode_id, comm)

USE qes_types_module, ONLY : control_variables_type
IMPLICIT NONE
TYPE(control_variables_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%title, ionode_id, comm)
CALL mp_bcast(o%calculation, ionode_id, comm)
CALL mp_bcast(o%restart_mode, ionode_id, comm)
CALL mp_bcast(o%prefix, ionode_id, comm)
CALL mp_bcast(o%pseudo_dir, ionode_id, comm)
CALL mp_bcast(o%outdir, ionode_id, comm)
CALL mp_bcast(o%stress, ionode_id, comm)
CALL mp_bcast(o%forces, ionode_id, comm)
CALL mp_bcast(o%wf_collect, ionode_id, comm)
CALL mp_bcast(o%disk_io, ionode_id, comm)
CALL mp_bcast(o%max_seconds, ionode_id, comm)
CALL mp_bcast(o%nstep_ispresent, ionode_id, comm)
IF (o%nstep_ispresent) &
   CALL mp_bcast(o%nstep, ionode_id, comm)
CALL mp_bcast(o%etot_conv_thr, ionode_id, comm)
CALL mp_bcast(o%forc_conv_thr, ionode_id, comm)
CALL mp_bcast(o%press_conv_thr, ionode_id, comm)
CALL mp_bcast(o%verbosity, ionode_id, comm)
CALL mp_bcast(o%print_every, ionode_id, comm)

RETURN
END SUBROUTINE bcast_control_variables_type
    !
SUBROUTINE bcast_species_type(o, ionode_id, comm)

USE qes_types_module, ONLY : species_type
IMPLICIT NONE
TYPE(species_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%name, ionode_id, comm)
CALL mp_bcast(o%mass_ispresent, ionode_id, comm)
IF (o%mass_ispresent) &
   CALL mp_bcast(o%mass, ionode_id, comm)
CALL mp_bcast(o%pseudo_file, ionode_id, comm)
CALL mp_bcast(o%starting_magnetization_ispresent, ionode_id, comm)
IF (o%starting_magnetization_ispresent) &
   CALL mp_bcast(o%starting_magnetization, ionode_id, comm)
CALL mp_bcast(o%spin_teta_ispresent, ionode_id, comm)
IF (o%spin_teta_ispresent) &
   CALL mp_bcast(o%spin_teta, ionode_id, comm)
CALL mp_bcast(o%spin_phi_ispresent, ionode_id, comm)
IF (o%spin_phi_ispresent) &
   CALL mp_bcast(o%spin_phi, ionode_id, comm)

RETURN
END SUBROUTINE bcast_species_type

SUBROUTINE bcast_atomic_positions_type(o, ionode_id, comm)

USE qes_types_module, ONLY : atomic_positions_type
IMPLICIT NONE
TYPE(atomic_positions_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

INTEGER :: na
CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%ndim_atom, ionode_id, comm)
IF (.NOT.ionode) ALLOCATE(o%atom(o%ndim_atom))
DO na=1,o%ndim_atom
   CALL bcast_atom_type(o%atom(na), ionode_id, comm)
ENDDO
RETURN
END SUBROUTINE bcast_atomic_positions_type

SUBROUTINE bcast_wyckoff_positions_type(o, ionode_id, comm)

USE qes_types_module, ONLY : wyckoff_positions_type
IMPLICIT NONE
TYPE(wyckoff_positions_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

INTEGER :: na
CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%space_group, ionode_id, comm)
CALL mp_bcast(o%more_options_ispresent, ionode_id, comm)
IF (o%more_options_ispresent) &
   CALL mp_bcast(o%more_options, ionode_id, comm)
CALL mp_bcast(o%ndim_atom, ionode_id, comm)
IF (.NOT.ionode) ALLOCATE(o%atom(o%ndim_atom))
DO na=1,o%ndim_atom
   CALL bcast_atom_type(o%atom(na), ionode_id, comm)
ENDDO

RETURN
END SUBROUTINE bcast_wyckoff_positions_type

SUBROUTINE bcast_cell_type(o, ionode_id, comm)

USE qes_types_module, ONLY : cell_type
IMPLICIT NONE
TYPE(cell_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%a1, ionode_id, comm)
CALL mp_bcast(o%a2, ionode_id, comm)
CALL mp_bcast(o%a3, ionode_id, comm)

RETURN
END SUBROUTINE bcast_cell_type

SUBROUTINE bcast_hybrid_type(o, ionode_id, comm)

USE qes_types_module, ONLY : hybrid_type
IMPLICIT NONE
TYPE(hybrid_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL bcast_qpoint_grid_type(o%qpoint_grid, ionode_id, comm)
CALL mp_bcast(o%ecutfock, ionode_id, comm)
CALL mp_bcast(o%exx_fraction, ionode_id, comm)
CALL mp_bcast(o%screening_parameter, ionode_id, comm)
CALL mp_bcast(o%exxdiv_treatment, ionode_id, comm)
CALL mp_bcast(o%x_gamma_extrapolation, ionode_id, comm)
CALL mp_bcast(o%ecutvcut, ionode_id, comm)

RETURN
END SUBROUTINE bcast_hybrid_type
    !
SUBROUTINE bcast_dftu_type(o, ionode_id, comm)

USE qes_types_module, ONLY : dftu_type
IMPLICIT NONE
TYPE(dftu_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

INTEGER :: i

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%lda_plus_u_kind_ispresent, ionode_id, comm)
IF (o%lda_plus_u_kind_ispresent) &
   CALL mp_bcast(o%lda_plus_u_kind, ionode_id, comm)

CALL mp_bcast(o%Hubbard_U_ispresent, ionode_id, comm)
IF (o%Hubbard_U_ispresent) THEN
   CALL mp_bcast(o%ndim_Hubbard_U, ionode_id, comm)
   IF (.NOT.ionode) ALLOCATE(o%Hubbard_U(o%ndim_Hubbard_U))
   DO i=1, o%ndim_Hubbard_U
      CALL bcast_HubbardCommon_type(o%Hubbard_U(i), ionode_id, comm)
   ENDDO
ENDIF

CALL mp_bcast(o%Hubbard_alpha_ispresent, ionode_id, comm)
IF (o%Hubbard_alpha_ispresent) THEN
   CALL mp_bcast(o%ndim_Hubbard_alpha, ionode_id, comm)
   IF (.NOT.ionode) ALLOCATE(o%Hubbard_alpha(o%ndim_Hubbard_alpha))
   DO i=1, o%ndim_Hubbard_alpha
      CALL bcast_HubbardCommon_type(o%Hubbard_alpha(i), ionode_id, comm)
   ENDDO
ENDIF

CALL mp_bcast(o%Hubbard_beta_ispresent, ionode_id, comm)
IF (o%Hubbard_beta_ispresent) THEN
   CALL mp_bcast(o%ndim_Hubbard_beta, ionode_id, comm)
   IF (.NOT.ionode) ALLOCATE(o%Hubbard_beta(o%ndim_Hubbard_beta))
   DO i=1, o%ndim_Hubbard_beta
      CALL bcast_HubbardCommon_type(o%Hubbard_beta(i), ionode_id, comm)
   ENDDO
ENDIF

CALL mp_bcast(o%Hubbard_j_ispresent, ionode_id, comm)
IF (o%Hubbard_j_ispresent) THEN
   CALL mp_bcast(o%ndim_Hubbard_j, ionode_id, comm)
   IF (.NOT.ionode) ALLOCATE(o%Hubbard_j(o%ndim_Hubbard_j))
   DO i=1, o%ndim_Hubbard_j
      CALL bcast_Hubbardj_type(o%Hubbard_j(i), ionode_id, comm)
   ENDDO
ENDIF

CALL mp_bcast(o%starting_ns_ispresent, ionode_id, comm)
IF (o%starting_ns_ispresent) THEN
   CALL mp_bcast(o%ndim_starting_ns, ionode_id, comm)
   IF (.NOT.ionode) ALLOCATE(o%starting_ns(o%ndim_starting_ns))
   DO i=1, o%ndim_starting_ns
      CALL bcast_starting_ns_type(o%starting_ns(i), ionode_id, comm)
   ENDDO
ENDIF

CALL mp_bcast(o%hubbard_ns_ispresent, ionode_id, comm)
IF (o%hubbard_ns_ispresent) THEN
   CALL mp_bcast(o%ndim_hubbard_ns, ionode_id, comm)
   IF (.NOT.ionode) ALLOCATE(o%hubbard_ns(o%ndim_hubbard_ns))
   DO i=1, o%ndim_hubbard_ns
      CALL bcast_hubbard_ns_type(o%hubbard_ns(i), ionode_id, comm)
   ENDDO
ENDIF

CALL mp_bcast(o%u_projection_type_ispresent, ionode_id, comm)
IF (o%u_projection_type_ispresent) &
   CALL mp_bcast(o%u_projection_type, ionode_id, comm)

RETURN
END SUBROUTINE bcast_dftu_type

SUBROUTINE bcast_vdw_type(o, ionode_id, comm)

USE qes_types_module, ONLY : vdw_type
IMPLICIT NONE
TYPE(vdw_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

INTEGER :: i

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%vdw_corr, ionode_id, comm)
CALL mp_bcast(o%non_local_term_ispresent, ionode_id, comm)
IF (o%non_local_term_ispresent) &
   CALL mp_bcast(o%non_local_term, ionode_id, comm)
CALL mp_bcast(o%london_s6_ispresent, ionode_id, comm)
IF (o%london_s6_ispresent) &
   CALL mp_bcast(o%london_s6, ionode_id, comm)
CALL mp_bcast(o%ts_vdw_econv_thr_ispresent, ionode_id, comm)
IF (o%ts_vdw_econv_thr_ispresent) &
   CALL mp_bcast(o%ts_vdw_econv_thr, ionode_id, comm)
CALL mp_bcast(o%ts_vdw_isolated_ispresent, ionode_id, comm)
IF (o%ts_vdw_isolated_ispresent) &
   CALL mp_bcast(o%ts_vdw_isolated, ionode_id, comm)
CALL mp_bcast(o%london_rcut_ispresent, ionode_id, comm)
IF (o%london_rcut_ispresent) &
   CALL mp_bcast(o%london_rcut, ionode_id, comm)
CALL mp_bcast(o%xdm_a1_ispresent, ionode_id, comm)
IF (o%xdm_a1_ispresent) &
   CALL mp_bcast(o%xdm_a1, ionode_id, comm)
CALL mp_bcast(o%xdm_a2_ispresent, ionode_id, comm)
IF (o%xdm_a2_ispresent) &
   CALL mp_bcast(o%xdm_a2, ionode_id, comm)

CALL mp_bcast(o%london_c6_ispresent , ionode_id, comm)
IF (o%london_c6_ispresent) THEN
   CALL mp_bcast(o%ndim_london_c6 , ionode_id, comm)
   IF (.NOT.ionode) ALLOCATE(o%london_c6(o%ndim_london_c6))
   DO i=1, o%ndim_london_c6
      CALL bcast_hubbardcommon_type(o%london_c6(i), ionode_id, comm)
   ENDDO
ENDIF

RETURN
END SUBROUTINE bcast_vdw_type

SUBROUTINE bcast_spin_type(o, ionode_id, comm)

USE qes_types_module, ONLY : spin_type
IMPLICIT NONE
TYPE(spin_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%lsda, ionode_id, comm)
CALL mp_bcast(o%noncolin, ionode_id, comm)
CALL mp_bcast(o%spinorbit, ionode_id, comm)

RETURN
END SUBROUTINE bcast_spin_type

SUBROUTINE bcast_bands_type(o, ionode_id, comm)

USE qes_types_module, ONLY : bands_type
IMPLICIT NONE
TYPE(bands_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

INTEGER :: i
CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%nbnd_ispresent, ionode_id, comm)
IF (o%nbnd_ispresent) &
   CALL mp_bcast(o%nbnd, ionode_id, comm)
CALL mp_bcast(o%smearing_ispresent, ionode_id, comm)
IF (o%smearing_ispresent) &
   CALL bcast_smearing_type(o%smearing, ionode_id, comm)
CALL mp_bcast(o%tot_charge_ispresent, ionode_id, comm)
IF (o%tot_charge_ispresent) &
   CALL mp_bcast(o%tot_charge, ionode_id, comm)
CALL mp_bcast(o%tot_magnetization_ispresent, ionode_id, comm)
IF (o%tot_magnetization_ispresent) &
   CALL mp_bcast(o%tot_magnetization, ionode_id, comm)
CALL bcast_occupations_type(o%occupations, ionode_id, comm)
CALL mp_bcast(o%inputoccupations_ispresent, ionode_id, comm)
IF (o%inputoccupations_ispresent) THEN
   CALL mp_bcast(o%ndim_inputOccupations, ionode_id, comm)
   IF (.NOT.ionode) ALLOCATE(o%inputoccupations(o%ndim_inputOccupations))
   DO i=1, o%ndim_inputOccupations
      CALL bcast_inputoccupations_type(o%inputoccupations(i), ionode_id, comm)
   ENDDO
ENDIF

RETURN
END SUBROUTINE bcast_bands_type
    !
SUBROUTINE bcast_basis_type(o, ionode_id, comm)

USE qes_types_module, ONLY : basis_type
IMPLICIT NONE
TYPE(basis_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%gamma_only_ispresent, ionode_id, comm)
IF (o%gamma_only_ispresent) &
   CALL mp_bcast(o%gamma_only, ionode_id, comm)
CALL mp_bcast(o%ecutwfc, ionode_id, comm)
CALL mp_bcast(o%ecutrho_ispresent, ionode_id, comm)
IF (o%ecutrho_ispresent) &
   CALL mp_bcast(o%ecutrho, ionode_id, comm)

CALL mp_bcast(o%fft_grid_ispresent, ionode_id, comm)
IF (o%fft_grid_ispresent) &
   CALL bcast_basissetitem_type(o%fft_grid, ionode_id, comm)
CALL mp_bcast(o%fft_smooth_ispresent, ionode_id, comm)
IF (o%fft_smooth_ispresent) &
   CALL bcast_basissetitem_type(o%fft_smooth, ionode_id, comm)
CALL mp_bcast(o%fft_box_ispresent, ionode_id, comm)
IF (o%fft_box_ispresent) &
   CALL bcast_basissetitem_type(o%fft_box, ionode_id, comm)

RETURN
END SUBROUTINE bcast_basis_type

SUBROUTINE bcast_reciprocal_lattice_type(o, ionode_id, comm)

USE qes_types_module, ONLY : reciprocal_lattice_type
IMPLICIT NONE
TYPE(reciprocal_lattice_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm
CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%b1, ionode_id, comm)
CALL mp_bcast(o%b2, ionode_id, comm)
CALL mp_bcast(o%b3, ionode_id, comm)

RETURN
END SUBROUTINE bcast_reciprocal_lattice_type

SUBROUTINE bcast_electron_control_type(o, ionode_id, comm)

USE qes_types_module, ONLY : electron_control_type
IMPLICIT NONE
TYPE(electron_control_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%diagonalization, ionode_id, comm)
CALL mp_bcast(o%mixing_mode, ionode_id, comm)
CALL mp_bcast(o%mixing_beta, ionode_id, comm)
CALL mp_bcast(o%conv_thr, ionode_id, comm)
CALL mp_bcast(o%mixing_ndim, ionode_id, comm)
CALL mp_bcast(o%max_nstep, ionode_id, comm)
CALL mp_bcast(o%real_space_q, ionode_id, comm)
CALL mp_bcast(o%tq_smoothing, ionode_id, comm)
CALL mp_bcast(o%tbeta_smoothing, ionode_id, comm)
CALL mp_bcast(o%diago_thr_init, ionode_id, comm)
CALL mp_bcast(o%diago_full_acc, ionode_id, comm)
CALL mp_bcast(o%diago_cg_maxiter_ispresent, ionode_id, comm)
CALL mp_bcast(o%diago_ppcg_maxiter_ispresent, ionode_id, comm)
IF (o%diago_cg_maxiter_ispresent) &
   CALL mp_bcast(o%diago_cg_maxiter, ionode_id, comm)
IF (o%diago_ppcg_maxiter_ispresent) &
   CALL mp_bcast(o%diago_ppcg_maxiter, ionode_id, comm)
CALL mp_bcast(o%diago_david_ndim_ispresent, ionode_id, comm)
IF (o%diago_david_ndim_ispresent) &
   CALL mp_bcast(o%diago_david_ndim, ionode_id, comm)

RETURN
END SUBROUTINE bcast_electron_control_type

SUBROUTINE bcast_k_points_IBZ_type(o, ionode_id, comm)

USE qes_types_module, ONLY : k_points_IBZ_type
IMPLICIT NONE
TYPE(k_points_IBZ_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

INTEGER :: i
CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%monkhorst_pack_ispresent, ionode_id, comm)
IF (o%monkhorst_pack_ispresent) &
   CALL bcast_monkhorst_pack_type(o%monkhorst_pack, ionode_id, comm)
CALL mp_bcast(o%nk_ispresent, ionode_id, comm)
IF (o%nk_ispresent) &
   CALL mp_bcast(o%nk, ionode_id, comm)
CALL mp_bcast(o%k_point_ispresent, ionode_id, comm)
IF (o%k_point_ispresent) THEN
   CALL mp_bcast(o%ndim_k_point, ionode_id, comm)
   IF (.NOT.ionode) ALLOCATE(o%k_point(o%ndim_k_point))
   DO i=1, o%ndim_k_point
      CALL bcast_k_point_type(o%k_point(i), ionode_id, comm)
   ENDDO
ENDIF

RETURN
END SUBROUTINE bcast_k_points_IBZ_type

SUBROUTINE bcast_bfgs_type(o, ionode_id, comm)

USE qes_types_module, ONLY : bfgs_type
IMPLICIT NONE
TYPE(bfgs_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%ndim, ionode_id, comm)
CALL mp_bcast(o%trust_radius_min, ionode_id, comm)
CALL mp_bcast(o%trust_radius_max, ionode_id, comm)
CALL mp_bcast(o%trust_radius_init, ionode_id, comm)

CALL mp_bcast(o%w1, ionode_id, comm)
CALL mp_bcast(o%w2, ionode_id, comm)

RETURN
END SUBROUTINE bcast_bfgs_type

SUBROUTINE bcast_md_type(o, ionode_id, comm)

USE qes_types_module, ONLY : md_type
IMPLICIT NONE
TYPE(md_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%pot_extrapolation, ionode_id, comm)
CALL mp_bcast(o%wfc_extrapolation, ionode_id, comm)
CALL mp_bcast(o%ion_temperature, ionode_id, comm)
CALL mp_bcast(o%timestep, ionode_id, comm)
CALL mp_bcast(o%tempw, ionode_id, comm)
CALL mp_bcast(o%tolp, ionode_id, comm)
CALL mp_bcast(o%deltat, ionode_id, comm)
CALL mp_bcast(o%nraise, ionode_id, comm)

RETURN
END SUBROUTINE bcast_md_type

SUBROUTINE bcast_cell_control_type(o, ionode_id, comm)

USE qes_types_module, ONLY : cell_control_type
IMPLICIT NONE
TYPE(cell_control_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%cell_dynamics, ionode_id, comm)
CALL mp_bcast(o%pressure, ionode_id, comm)
CALL mp_bcast(o%wmass_ispresent, ionode_id, comm)
IF (o%wmass_ispresent) &
   CALL mp_bcast(o%wmass, ionode_id, comm)
CALL mp_bcast(o%cell_factor_ispresent, ionode_id, comm)
IF (o%cell_factor_ispresent) &
   CALL mp_bcast(o%cell_factor, ionode_id, comm)
CALL mp_bcast(o%fix_volume_ispresent, ionode_id, comm)
IF (o%fix_volume_ispresent) &
   CALL mp_bcast(o%fix_volume, ionode_id, comm)
CALL mp_bcast(o%fix_area_ispresent, ionode_id, comm)
IF (o%fix_area_ispresent) &
   CALL mp_bcast(o%fix_area, ionode_id, comm)
CALL mp_bcast(o%isotropic_ispresent, ionode_id, comm)
IF (o%isotropic_ispresent) &
   CALL mp_bcast(o%isotropic, ionode_id, comm)
CALL mp_bcast(o%free_cell_ispresent, ionode_id, comm)
IF (o%free_cell_ispresent) &
   CALL bcast_integerMatrix_type(o%free_cell, ionode_id, comm)

RETURN
END SUBROUTINE bcast_cell_control_type

SUBROUTINE bcast_symmetry_flags_type(o, ionode_id, comm)

USE qes_types_module, ONLY : symmetry_flags_type
IMPLICIT NONE
TYPE(symmetry_flags_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%nosym, ionode_id, comm)
CALL mp_bcast(o%nosym_evc, ionode_id, comm)
CALL mp_bcast(o%noinv, ionode_id, comm)
CALL mp_bcast(o%no_t_rev, ionode_id, comm)
CALL mp_bcast(o%force_symmorphic, ionode_id, comm)
CALL mp_bcast(o%use_all_frac, ionode_id, comm)

RETURN
END SUBROUTINE bcast_symmetry_flags_type

SUBROUTINE bcast_esm_type(o, ionode_id, comm)

USE qes_types_module, ONLY : esm_type
IMPLICIT NONE
TYPE(esm_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%bc, ionode_id, comm)
CALL mp_bcast(o%nfit, ionode_id, comm)
CALL mp_bcast(o%w, ionode_id, comm)
CALL mp_bcast(o%efield, ionode_id, comm)

RETURN
END SUBROUTINE bcast_esm_type

SUBROUTINE bcast_ekin_functional_type(o, ionode_id, comm)

USE qes_types_module, ONLY : ekin_functional_type
IMPLICIT NONE
TYPE(ekin_functional_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%ecfixed, ionode_id, comm)
CALL mp_bcast(o%qcutz, ionode_id, comm)
CALL mp_bcast(o%q2sigma, ionode_id, comm)

RETURN
END SUBROUTINE bcast_ekin_functional_type

SUBROUTINE bcast_spin_constraints_type(o, ionode_id, comm)

USE qes_types_module, ONLY : spin_constraints_type
IMPLICIT NONE
TYPE(spin_constraints_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%spin_constraints, ionode_id, comm)
CALL mp_bcast(o%lagrange_multiplier, ionode_id, comm)
CALL mp_bcast(o%target_magnetization_ispresent, ionode_id, comm)
IF (o%target_magnetization_ispresent) &
   CALL mp_bcast(o%target_magnetization, ionode_id, comm)

RETURN
END SUBROUTINE bcast_spin_constraints_type

SUBROUTINE bcast_electric_field_type(o, ionode_id, comm)

USE qes_types_module, ONLY : electric_field_type
IMPLICIT NONE
TYPE(electric_field_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%electric_potential, ionode_id, comm)
CALL mp_bcast(o%dipole_correction_ispresent, ionode_id, comm)
IF (o%dipole_correction_ispresent) &
   CALL mp_bcast(o%dipole_correction, ionode_id, comm)
CALL mp_bcast(o%electric_field_direction_ispresent, ionode_id, comm)
IF (o%electric_field_direction_ispresent) &
   CALL mp_bcast(o%electric_field_direction, ionode_id, comm)

CALL mp_bcast(o%potential_max_position_ispresent, ionode_id, comm)
IF (o%potential_max_position_ispresent) &
   CALL mp_bcast(o%potential_max_position, ionode_id, comm)
CALL mp_bcast(o%potential_decrease_width_ispresent, ionode_id, comm)
IF (o%potential_decrease_width_ispresent) &
   CALL mp_bcast(o%potential_decrease_width, ionode_id, comm)
CALL mp_bcast(o%electric_field_amplitude_ispresent, ionode_id, comm)
IF (o%electric_field_amplitude_ispresent) &
   CALL mp_bcast(o%electric_field_amplitude, ionode_id, comm)
CALL mp_bcast(o%electric_field_vector_ispresent, ionode_id, comm)
IF (o%electric_field_vector_ispresent) &
   CALL mp_bcast(o%electric_field_vector, ionode_id, comm)
CALL mp_bcast(o%nk_per_string_ispresent, ionode_id, comm)
IF (o%nk_per_string_ispresent) &
   CALL mp_bcast(o%nk_per_string, ionode_id, comm)
CALL mp_bcast(o%n_berry_cycles_ispresent, ionode_id, comm)
IF (o%n_berry_cycles_ispresent) &
   CALL mp_bcast(o%n_berry_cycles, ionode_id, comm)

RETURN
END SUBROUTINE bcast_electric_field_type

SUBROUTINE bcast_atomic_constraint_type(o, ionode_id, comm)

USE qes_types_module, ONLY : atomic_constraint_type
IMPLICIT NONE
TYPE(atomic_constraint_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%constr_parms, ionode_id, comm)
CALL mp_bcast(o%constr_type, ionode_id, comm)
CALL mp_bcast(o%constr_target, ionode_id, comm)

RETURN
END SUBROUTINE bcast_atomic_constraint_type

SUBROUTINE bcast_dipoleOutput_type(o, ionode_id, comm)

USE qes_types_module, ONLY : dipoleOutput_type
IMPLICIT NONE
TYPE(dipoleOutput_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%idir, ionode_id, comm)
CALL bcast_scalarQuantity_type(o%dipole, ionode_id, comm)
CALL bcast_scalarQuantity_type(o%ion_dipole, ionode_id, comm)
CALL bcast_scalarQuantity_type(o%elec_dipole, ionode_id, comm)
CALL bcast_scalarQuantity_type(o%dipoleField, ionode_id, comm)
CALL bcast_scalarQuantity_type(o%potentialAmp, ionode_id, comm)
CALL bcast_scalarQuantity_type(o%totalLength, ionode_id, comm)

END SUBROUTINE bcast_dipoleOutput_type

SUBROUTINE bcast_finiteFieldOut_type(o, ionode_id, comm)

USE qes_types_module, ONLY : finiteFieldOut_type
IMPLICIT NONE
TYPE(finiteFieldOut_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%electronicDipole, ionode_id, comm)
CALL mp_bcast(o%ionicDipole, ionode_id, comm)

RETURN
END SUBROUTINE bcast_finiteFieldOut_type

SUBROUTINE bcast_polarization_type(o, ionode_id, comm)

USE qes_types_module, ONLY : polarization_type
IMPLICIT NONE
TYPE(Polarization_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL bcast_scalarQuantity_type(o%polarization, ionode_id, comm)
CALL mp_bcast(o%modulus, ionode_id, comm)
CALL mp_bcast(o%direction, ionode_id, comm)

RETURN
END SUBROUTINE bcast_polarization_type
    !
SUBROUTINE bcast_ionicPolarization_type(o, ionode_id, comm)

USE qes_types_module, ONLY : ionicPolarization_type
IMPLICIT NONE
TYPE(ionicPolarization_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm


CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL bcast_atom_type(o%ion, ionode_id, comm)
CALL mp_bcast(o%charge, ionode_id, comm)
CALL bcast_phase_type(o%phase, ionode_id, comm)
    !
RETURN
END SUBROUTINE bcast_ionicPolarization_type

SUBROUTINE bcast_electronicPolarization_type(o, ionode_id, comm)

USE qes_types_module, ONLY : electronicPolarization_type
IMPLICIT NONE
TYPE(electronicPolarization_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL bcast_k_point_type(o%firstKeyPoint, ionode_id, comm)
CALL mp_bcast(o%spin_ispresent, ionode_id, comm)
IF (o%spin_ispresent) &
   CALL mp_bcast(o%spin, ionode_id, comm)
CALL bcast_phase_type(o%phase, ionode_id, comm)
    !
RETURN
END SUBROUTINE bcast_electronicPolarization_type

SUBROUTINE bcast_scf_conv_type(o, ionode_id, comm)

USE qes_types_module, ONLY : scf_conv_type
IMPLICIT NONE
TYPE(scf_conv_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%n_scf_steps, ionode_id, comm)
CALL mp_bcast(o%scf_error, ionode_id, comm)
    !
RETURN
END SUBROUTINE bcast_scf_conv_type

SUBROUTINE bcast_opt_conv_type(o, ionode_id, comm)

USE qes_types_module, ONLY : opt_conv_type
IMPLICIT NONE
TYPE(opt_conv_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%n_opt_steps, ionode_id, comm)
CALL mp_bcast(o%grad_norm, ionode_id, comm)
    !
RETURN
END SUBROUTINE bcast_opt_conv_type

SUBROUTINE bcast_algorithmic_info_type(o, ionode_id, comm)

USE qes_types_module, ONLY : algorithmic_info_type
IMPLICIT NONE
TYPE(algorithmic_info_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%real_space_q, ionode_id, comm)
CALL mp_bcast(o%uspp, ionode_id, comm)
CALL mp_bcast(o%paw, ionode_id, comm)
    !
RETURN
END SUBROUTINE bcast_algorithmic_info_type

SUBROUTINE bcast_symmetry_type(o, ionode_id, comm)

USE qes_types_module, ONLY : symmetry_type
IMPLICIT NONE
TYPE(symmetry_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL bcast_info_type(o%info, ionode_id, comm)
CALL bcast_matrix_type(o%rotation, ionode_id, comm)
CALL mp_bcast(o%fractional_translation_ispresent, ionode_id, comm)
IF (o%fractional_translation_ispresent) &
   CALL mp_bcast(o%fractional_translation, ionode_id, comm)
CALL mp_bcast(o%equivalent_atoms_ispresent, ionode_id, comm)
IF (o%equivalent_atoms_ispresent) &
   CALL bcast_equivalent_atoms_type(o%equivalent_atoms, ionode_id, comm)
    !
RETURN
END SUBROUTINE bcast_symmetry_type

SUBROUTINE bcast_magnetization_type(o, ionode_id, comm)

USE qes_types_module, ONLY : magnetization_type
IMPLICIT NONE
TYPE(magnetization_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%lsda, ionode_id, comm)
CALL mp_bcast(o%noncolin, ionode_id, comm)
CALL mp_bcast(o%spinorbit, ionode_id, comm)
CALL mp_bcast(o%total, ionode_id, comm)
CALL mp_bcast(o%absolute, ionode_id, comm)
CALL mp_bcast(o%do_magnetization, ionode_id, comm)
    !
RETURN
END SUBROUTINE bcast_magnetization_type

SUBROUTINE bcast_total_energy_type(o, ionode_id, comm)

USE qes_types_module, ONLY : total_energy_type
IMPLICIT NONE
TYPE(total_energy_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%etot, ionode_id, comm)
CALL mp_bcast(o%eband_ispresent, ionode_id, comm)
IF (o%eband_ispresent) &
   CALL mp_bcast(o%eband, ionode_id, comm)
CALL mp_bcast(o%ehart_ispresent, ionode_id, comm)
IF (o%ehart_ispresent) &
   CALL mp_bcast(o%ehart, ionode_id, comm)
CALL mp_bcast(o%vtxc_ispresent, ionode_id, comm)
IF (o%vtxc_ispresent) &
   CALL mp_bcast(o%vtxc, ionode_id, comm)
CALL mp_bcast(o%etxc_ispresent, ionode_id, comm)
IF (o%etxc_ispresent) &
   CALL mp_bcast(o%etxc, ionode_id, comm)
CALL mp_bcast(o%ewald_ispresent, ionode_id, comm)
IF (o%ewald_ispresent) &
   CALL mp_bcast(o%ewald, ionode_id, comm)
CALL mp_bcast(o%demet_ispresent, ionode_id, comm)
IF (o%demet_ispresent) &
   CALL mp_bcast(o%demet, ionode_id, comm)
CALL mp_bcast(o%efieldcorr_ispresent, ionode_id, comm)
IF (o%efieldcorr_ispresent) &
   CALL mp_bcast(o%efieldcorr, ionode_id, comm)
CALL mp_bcast(o%potentiostat_contr_ispresent, ionode_id, comm)
IF (o%potentiostat_contr_ispresent) &
   CALL mp_bcast(o%potentiostat_contr, ionode_id, comm)

RETURN
END SUBROUTINE bcast_total_energy_type
    !
SUBROUTINE bcast_ks_energies_type(o, ionode_id, comm)

USE qes_types_module, ONLY : ks_energies_type
IMPLICIT NONE
TYPE(ks_energies_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL bcast_k_point_type(o%k_point, ionode_id, comm)
CALL mp_bcast(o%npw, ionode_id, comm)
CALL bcast_vector_type(o%eigenvalues, ionode_id, comm)
CALL bcast_vector_type(o%occupations, ionode_id, comm)

RETURN
END SUBROUTINE bcast_ks_energies_type

    !
SUBROUTINE bcast_atomic_species_type(o, ionode_id, comm)

USE qes_types_module, ONLY : atomic_species_type
IMPLICIT NONE
TYPE(atomic_species_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

INTEGER :: i

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%ntyp, ionode_id, comm)
CALL mp_bcast(o%pseudo_dir_ispresent, ionode_id, comm)
IF (o%pseudo_dir_ispresent) &
   CALL mp_bcast(o%pseudo_dir, ionode_id, comm)
CALL mp_bcast(o%ndim_species, ionode_id, comm)
IF (.NOT.ionode) ALLOCATE(o%species(o%ndim_species)) 
DO i=1, o%ndim_species
   CALL bcast_species_type(o%species(i), ionode_id, comm)
ENDDO

RETURN
END SUBROUTINE bcast_atomic_species_type

SUBROUTINE bcast_atomic_structure_type(o, ionode_id, comm)

USE qes_types_module, ONLY : atomic_structure_type
IMPLICIT NONE
TYPE(atomic_structure_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%nat, ionode_id, comm)
CALL mp_bcast(o%alat_ispresent, ionode_id, comm)
IF (o%alat_ispresent) &
   CALL mp_bcast(o%alat, ionode_id, comm)
CALL mp_bcast(o%bravais_index_ispresent, ionode_id, comm)
IF (o%bravais_index_ispresent) &
   CALL mp_bcast(o%bravais_index, ionode_id, comm)
CALL mp_bcast(o%atomic_positions_ispresent, ionode_id, comm)
IF (o%atomic_positions_ispresent) &
   CALL bcast_atomic_positions_type(o%atomic_positions, ionode_id, comm)
CALL mp_bcast(o%wyckoff_positions_ispresent, ionode_id, comm)
IF (o%wyckoff_positions_ispresent) &
   CALL bcast_wyckoff_positions_type(o%wyckoff_positions, ionode_id, comm)
CALL bcast_cell_type(o%cell, ionode_id, comm)

RETURN
END SUBROUTINE bcast_atomic_structure_type

SUBROUTINE bcast_dft_type(o, ionode_id, comm)

USE qes_types_module, ONLY : dft_type
IMPLICIT NONE
TYPE(dft_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%functional, ionode_id, comm)
CALL mp_bcast(o%hybrid_ispresent, ionode_id, comm)
IF (o%hybrid_ispresent) &
   CALL bcast_hybrid_type(o%hybrid, ionode_id, comm)
CALL mp_bcast(o%dftU_ispresent, ionode_id, comm)
IF (o%dftU_ispresent) &
   CALL bcast_dftU_type(o%dftU, ionode_id, comm)
CALL mp_bcast(o%vdW_ispresent, ionode_id, comm)
IF (o%vdW_ispresent) &
   CALL bcast_vdW_type(o%vdW, ionode_id, comm)

RETURN
END SUBROUTINE bcast_dft_type
    !

SUBROUTINE bcast_basis_set_type(o, ionode_id, comm)

USE qes_types_module, ONLY : basis_set_type
IMPLICIT NONE
TYPE(basis_set_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%gamma_only_ispresent, ionode_id, comm)
IF (o%gamma_only_ispresent) &
   CALL mp_bcast(o%gamma_only, ionode_id, comm)
CALL mp_bcast(o%ecutwfc, ionode_id, comm)
CALL mp_bcast(o%ecutrho_ispresent, ionode_id, comm)
IF (o%ecutrho_ispresent) &
   CALL mp_bcast(o%ecutrho, ionode_id, comm)
CALL bcast_basisSetItem_type(o%fft_grid, ionode_id, comm)
CALL mp_bcast(o%fft_smooth_ispresent, ionode_id, comm)
IF (o%fft_smooth_ispresent) &
   CALL bcast_basisSetItem_type(o%fft_smooth, ionode_id, comm)
CALL mp_bcast(o%fft_box_ispresent, ionode_id, comm)
IF (o%fft_box_ispresent) &
   CALL bcast_basisSetItem_type(o%fft_box, ionode_id, comm)
CALL mp_bcast(o%ngm, ionode_id, comm)
CALL mp_bcast(o%ngms_ispresent, ionode_id, comm)
IF (o%ngms_ispresent) &
   CALL mp_bcast(o%ngms, ionode_id, comm)
CALL mp_bcast(o%npwx, ionode_id, comm)
CALL bcast_reciprocal_lattice_type(o%reciprocal_lattice, ionode_id, comm)
    !
RETURN
END SUBROUTINE bcast_basis_set_type


SUBROUTINE bcast_ion_control_type(o, ionode_id, comm)

USE qes_types_module, ONLY : ion_control_type
IMPLICIT NONE
TYPE(ion_control_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%ion_dynamics, ionode_id, comm)
CALL mp_bcast(o%upscale_ispresent, ionode_id, comm)
IF (o%upscale_ispresent) &
   CALL mp_bcast(o%upscale, ionode_id, comm)
CALL mp_bcast(o%remove_rigid_rot_ispresent, ionode_id, comm)
IF (o%remove_rigid_rot_ispresent) &
   CALL mp_bcast(o%remove_rigid_rot, ionode_id, comm)
CALL mp_bcast(o%refold_pos_ispresent, ionode_id, comm)
IF (o%refold_pos_ispresent) &
   CALL mp_bcast(o%refold_pos, ionode_id, comm)
CALL mp_bcast(o%bfgs_ispresent, ionode_id, comm)
IF (o%bfgs_ispresent) &
   CALL bcast_bfgs_type(o%bfgs, ionode_id, comm)
CALL mp_bcast(o%md_ispresent, ionode_id, comm)
IF (o%md_ispresent) &
   CALL bcast_md_type(o%md, ionode_id, comm)

RETURN
END SUBROUTINE bcast_ion_control_type

SUBROUTINE bcast_boundary_conditions_type(o, ionode_id, comm)

USE qes_types_module, ONLY : boundary_conditions_type
IMPLICIT NONE
TYPE(boundary_conditions_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%assume_isolated, ionode_id, comm)
CALL mp_bcast(o%esm_ispresent, ionode_id, comm)
IF (o%esm_ispresent) &
   CALL bcast_esm_type(o%esm, ionode_id, comm)
CALL mp_bcast(o%fcp_opt_ispresent, ionode_id, comm)
IF (o%fcp_opt_ispresent) &
   CALL mp_bcast(o%fcp_opt, ionode_id, comm)
CALL mp_bcast(o%fcp_mu_ispresent, ionode_id, comm)
IF (o%fcp_mu_ispresent) &
   CALL mp_bcast(o%fcp_mu, ionode_id, comm)

RETURN
END SUBROUTINE bcast_boundary_conditions_type
    !
SUBROUTINE bcast_atomic_constraints_type(o, ionode_id, comm)

USE qes_types_module, ONLY : atomic_constraints_type
IMPLICIT NONE
TYPE(atomic_constraints_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm


INTEGER :: i
CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%num_of_constraints, ionode_id, comm)
CALL mp_bcast(o%tolerance, ionode_id, comm)
CALL mp_bcast(o%ndim_atomic_constraint, ionode_id, comm)
IF (.NOT.ionode) ALLOCATE(o%atomic_constraint(o%ndim_atomic_constraint))
DO i=1, o%ndim_atomic_constraint
   CALL bcast_atomic_constraint_type(o%atomic_constraint(i), ionode_id, comm)
ENDDO

RETURN
END SUBROUTINE bcast_atomic_constraints_type

    !
SUBROUTINE bcast_BerryPhaseOutput_type(o, ionode_id, comm)

USE qes_types_module, ONLY : BerryPhaseOutput_type
IMPLICIT NONE
TYPE(BerryPhaseOutput_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

INTEGER :: i

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL bcast_polarization_type(o%totalPolarization, ionode_id, comm)
CALL bcast_phase_type(o%totalPhase, ionode_id, comm)
CALL mp_bcast(o%ndim_ionicPolarization, ionode_id, comm)
IF (.NOT. ionode) ALLOCATE(o%ionicPolarization(o%ndim_ionicPolarization))
DO i=1, o%ndim_ionicPolarization
   CALL bcast_ionicPolarization_type(o%ionicPolarization(i), ionode_id, comm)
ENDDO
CALL mp_bcast(o%ndim_electronicPolarization, ionode_id, comm)
IF (.NOT. ionode) &
   ALLOCATE(o%electronicPolarization (o%ndim_electronicPolarization))
DO i=1, o%ndim_electronicPolarization
   CALL bcast_electronicPolarization_type(o%electronicPolarization(i), &
                                                            ionode_id, comm)
ENDDO

RETURN
END SUBROUTINE bcast_BerryPhaseOutput_type

SUBROUTINE bcast_convergence_info_type(o, ionode_id, comm)

USE qes_types_module, ONLY : convergence_info_type
IMPLICIT NONE
TYPE(convergence_info_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL bcast_scf_conv_type(o%scf_conv, ionode_id, comm)
CALL mp_bcast(o%opt_conv_ispresent, ionode_id, comm)
IF (o%opt_conv_ispresent) &
   CALL bcast_opt_conv_type(o%opt_conv, ionode_id, comm)

RETURN
END SUBROUTINE bcast_convergence_info_type

    !
SUBROUTINE bcast_symmetries_type(o, ionode_id, comm)

USE qes_types_module, ONLY : symmetries_type
IMPLICIT NONE
TYPE(symmetries_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

INTEGER :: i
CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%nsym, ionode_id, comm)
CALL mp_bcast(o%nrot, ionode_id, comm)
CALL mp_bcast(o%space_group, ionode_id, comm)
CALL mp_bcast(o%ndim_symmetry, ionode_id, comm)
IF (.NOT. ionode) ALLOCATE(o%symmetry(o%ndim_symmetry))
DO i=1,o%ndim_symmetry
   CALL bcast_symmetry_type(o%symmetry(i), ionode_id, comm)
ENDDO

RETURN
END SUBROUTINE bcast_symmetries_type

SUBROUTINE bcast_band_structure_type(o, ionode_id, comm)

USE qes_types_module, ONLY : band_structure_type
IMPLICIT NONE
TYPE(band_structure_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

INTEGER :: i

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%lsda, ionode_id, comm)
CALL mp_bcast(o%noncolin, ionode_id, comm)
CALL mp_bcast(o%spinorbit, ionode_id, comm)
CALL mp_bcast(o%nbnd, ionode_id, comm)
CALL mp_bcast(o%nbnd_up_ispresent, ionode_id, comm)
IF (o%nbnd_up_ispresent) &
   CALL mp_bcast(o%nbnd_up, ionode_id, comm)
CALL mp_bcast(o%nbnd_dw_ispresent, ionode_id, comm)
IF (o%nbnd_dw_ispresent) &
   CALL mp_bcast(o%nbnd_dw, ionode_id, comm)
CALL mp_bcast(o%nelec, ionode_id, comm)
CALL mp_bcast(o%num_of_atomic_wfc_ispresent, ionode_id, comm)
IF (o%num_of_atomic_wfc_ispresent) &
   CALL mp_bcast(o%num_of_atomic_wfc, ionode_id, comm)
CALL mp_bcast(o%wf_collected, ionode_id, comm)
CALL mp_bcast(o%fermi_energy_ispresent, ionode_id, comm)
IF (o%fermi_energy_ispresent) &
   CALL mp_bcast(o%fermi_energy, ionode_id, comm)
CALL mp_bcast(o%highestOccupiedLevel_ispresent, ionode_id, comm)
IF (o%highestOccupiedLevel_ispresent) &
   CALL mp_bcast(o%highestOccupiedLevel, ionode_id, comm)
CALL mp_bcast(o%two_fermi_energies_ispresent, ionode_id, comm)
IF (o%two_fermi_energies_ispresent) &
   CALL mp_bcast(o%two_fermi_energies, ionode_id, comm)
CALL bcast_k_points_IBZ_type(o%starting_k_points, ionode_id, comm)
CALL mp_bcast(o%nks, ionode_id, comm)
CALL bcast_occupations_type(o%occupations_kind, ionode_id, comm)
CALL mp_bcast(o%smearing_ispresent, ionode_id, comm)
IF (o%smearing_ispresent) &
   CALL bcast_smearing_type(o%smearing, ionode_id, comm)
CALL mp_bcast(o%ndim_ks_energies, ionode_id, comm)
IF (.NOT.ionode) ALLOCATE(o%ks_energies(o%ndim_ks_energies))
DO i=1, o%ndim_ks_energies
   CALL bcast_ks_energies_type(o%ks_energies(i), ionode_id, comm)
ENDDO

RETURN
END SUBROUTINE bcast_band_structure_type

SUBROUTINE bcast_input_type(o, ionode_id, comm)

USE qes_types_module, ONLY : input_type
IMPLICIT NONE
TYPE(input_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL bcast_control_variables_type(o%control_variables, ionode_id, comm)
CALL bcast_atomic_species_type(o%atomic_species, ionode_id, comm)
CALL bcast_atomic_structure_type(o%atomic_structure, ionode_id, comm)
CALL bcast_dft_type(o%dft, ionode_id, comm)
CALL bcast_spin_type(o%spin, ionode_id, comm)
CALL bcast_bands_type(o%bands, ionode_id, comm)
CALL bcast_basis_type(o%basis, ionode_id, comm)
CALL bcast_electron_control_type(o%electron_control, ionode_id, comm)
CALL bcast_k_points_IBZ_type(o%k_points_IBZ, ionode_id, comm)
CALL bcast_ion_control_type(o%ion_control, ionode_id, comm)
CALL bcast_cell_control_type(o%cell_control, ionode_id, comm)
CALL mp_bcast(o%symmetry_flags_ispresent, ionode_id, comm)
IF (o%symmetry_flags_ispresent) &
   CALL bcast_symmetry_flags_type(o%symmetry_flags, ionode_id, comm)
CALL mp_bcast(o%boundary_conditions_ispresent, ionode_id, comm)
IF (o%boundary_conditions_ispresent) &
   CALL bcast_boundary_conditions_type(o%boundary_conditions, ionode_id, comm)
CALL mp_bcast(o%ekin_functional_ispresent, ionode_id, comm)
IF (o%ekin_functional_ispresent) &
   CALL bcast_ekin_functional_type(o%ekin_functional, ionode_id, comm)
CALL mp_bcast(o%external_atomic_forces_ispresent, ionode_id, comm)
IF (o%external_atomic_forces_ispresent) &
   CALL bcast_matrix_type(o%external_atomic_forces, ionode_id, comm)
CALL mp_bcast(o%free_positions_ispresent, ionode_id, comm)
IF (o%free_positions_ispresent) &
   CALL bcast_integermatrix_type(o%free_positions, ionode_id, comm)
CALL mp_bcast(o%starting_atomic_velocities_ispresent, ionode_id, comm)
IF (o%starting_atomic_velocities_ispresent) &
   CALL bcast_matrix_type(o%starting_atomic_velocities, ionode_id, comm)
CALL mp_bcast(o%electric_field_ispresent, ionode_id, comm)
IF (o%electric_field_ispresent) &
   CALL bcast_electric_field_type(o%electric_field, ionode_id, comm)
CALL mp_bcast(o%atomic_constraints_ispresent, ionode_id, comm)
IF (o%atomic_constraints_ispresent) &
   CALL bcast_atomic_constraints_type(o%atomic_constraints, ionode_id, comm)
CALL mp_bcast(o%spin_constraints_ispresent, ionode_id, comm)
IF (o%spin_constraints_ispresent) &
   CALL bcast_spin_constraints_type(o%spin_constraints, ionode_id, comm)
    !
RETURN
END SUBROUTINE bcast_input_type

SUBROUTINE bcast_step_type(o, ionode_id, comm)

USE qes_types_module, ONLY : step_type
IMPLICIT NONE
TYPE(step_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm


CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)
CALL mp_bcast(o%n_step, ionode_id, comm)
CALL bcast_scf_conv_type(o%scf_conv, ionode_id, comm)
CALL bcast_atomic_structure_type(o%atomic_structure, ionode_id, comm)
CALL bcast_total_energy_type(o%total_energy, ionode_id, comm)
CALL bcast_matrix_type(o%forces, ionode_id, comm)
CALL mp_bcast(o%stress_ispresent, ionode_id, comm)
IF (o%stress_ispresent) &
   CALL bcast_matrix_type(o%stress, ionode_id, comm)
CALL mp_bcast(o%FCP_force_ispresent, ionode_id, comm)
IF (o%FCP_force_ispresent) &
   CALL mp_bcast(o%FCP_force, ionode_id, comm)
CALL mp_bcast(o%FCP_tot_charge_ispresent, ionode_id, comm)
IF (o%FCP_tot_charge_ispresent) &
   CALL mp_bcast(o%FCP_tot_charge, ionode_id, comm)

RETURN
END SUBROUTINE bcast_step_type

SUBROUTINE bcast_outputElectricField_type(o, ionode_id, comm)

USE qes_types_module, ONLY : outputElectricField_type
IMPLICIT NONE
TYPE(outputElectricField_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL mp_bcast(o%tagname, ionode_id, comm)
CALL mp_bcast(o%lwrite, ionode_id, comm)
CALL mp_bcast(o%lread, ionode_id, comm)

CALL mp_bcast(o%BerryPhase_ispresent, ionode_id, comm)
IF (o%BerryPhase_ispresent) &
   CALL bcast_BerryPhaseOutput_type(o%BerryPhase, ionode_id, comm)

CALL mp_bcast(o%finiteElectricFieldInfo_ispresent, ionode_id, comm)
IF (o%finiteElectricFieldInfo_ispresent) &
   CALL bcast_finiteFieldOut_type(o%finiteElectricFieldInfo, &
                                                          ionode_id, comm)
CALL mp_bcast(o%dipoleInfo_ispresent, ionode_id, comm)
IF (o%dipoleInfo_ispresent) &
   CALL bcast_dipoleOutput_type(o%dipoleInfo, ionode_id, comm)

RETURN
END SUBROUTINE bcast_outputElectricField_type

SUBROUTINE bcast_output_type(o, ionode_id, comm)

USE qes_types_module, ONLY : output_type
IMPLICIT NONE
TYPE(output_type), INTENT(INOUT) :: o
INTEGER, INTENT(IN) :: ionode_id, comm

CALL bcast_convergence_info_type(o%convergence_info, ionode_id, comm)
CALL bcast_algorithmic_info_type(o%algorithmic_info, ionode_id, comm)
CALL bcast_atomic_species_type(o%atomic_species, ionode_id, comm)
CALL bcast_atomic_structure_type(o%atomic_structure, ionode_id, comm)
CALL mp_bcast(o%symmetries_ispresent, ionode_id, comm)
IF (o%symmetries_ispresent) &
   CALL bcast_symmetries_type(o%symmetries, ionode_id, comm)
CALL bcast_basis_set_type(o%basis_set, ionode_id, comm)
CALL bcast_dft_type(o%dft, ionode_id, comm)
CALL bcast_magnetization_type(o%magnetization, ionode_id, comm)
CALL bcast_total_energy_type(o%total_energy, ionode_id, comm)
CALL bcast_band_structure_type(o%band_structure, ionode_id, comm)
CALL mp_bcast(o%forces_ispresent, ionode_id, comm)
IF (o%forces_ispresent) &
   CALL bcast_matrix_type(o%forces, ionode_id, comm)
CALL mp_bcast(o%stress_ispresent, ionode_id, comm)
IF (o%stress_ispresent) &
   CALL bcast_matrix_type(o%stress, ionode_id, comm)
CALL mp_bcast(o%electric_field_ispresent, ionode_id, comm)
IF (o%electric_field_ispresent) &
   CALL bcast_outputelectricfield_type(o%electric_field, ionode_id, comm)
CALL mp_bcast(o%fcp_force_ispresent, ionode_id, comm)
IF (o%fcp_force_ispresent) &
   CALL mp_bcast(o%fcp_force, ionode_id, comm)
CALL mp_bcast(o%fcp_tot_charge_ispresent, ionode_id, comm)
IF (o%fcp_tot_charge_ispresent) &
   CALL mp_bcast(o%fcp_tot_charge, ionode_id, comm)

RETURN
END SUBROUTINE bcast_output_type

  !
END MODULE bcast_qes_types_module
