!
! Copyright (C) 2003-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE qes_libs_module
! This module contains some basic subroutines initialize data_type used for
! reading and writing  XML files produced by Quantum ESPRESSO package.
!
! Written by Giovanni Borghi, A. Ferretti, ... (2015).
!

   USE qes_types_module
   USE iotk_module
   !
   INTEGER, PARAMETER      :: max_real_per_line=5
   CHARACTER(iotk_attlenx) :: attr
   CHARACTER(32)           :: fmtstr
   !
   PRIVATE :: attr, fmtstr
!
CONTAINS
!
SUBROUTINE qes_write_closed(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(closed_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'DATE', TRIM(obj%DATE))
   CALL iotk_write_attr(attr, 'TIME', TRIM(obj%TIME))

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr),new_line=.FALSE.)
      !
      WRITE(iun, '(A)',advance='no')  TRIM(obj%closed)
   CALL iotk_write_end(iun, TRIM(obj%tagname),indentation=.FALSE.)
   !
END SUBROUTINE qes_write_closed

SUBROUTINE qes_init_closed(obj, tagname, DATE, TIME, closed)
   IMPLICIT NONE

   TYPE(closed_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: DATE
   CHARACTER(len=*) :: TIME
   CHARACTER(len=*) :: closed

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%DATE = TRIM(DATE)


   obj%TIME = TRIM(TIME)

   obj%closed = closed

END SUBROUTINE qes_init_closed

SUBROUTINE qes_reset_closed(obj)
   IMPLICIT NONE
   TYPE(closed_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_closed


SUBROUTINE qes_write_status(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(status_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
         WRITE(iun, '(I12)') obj%status
         !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_status

SUBROUTINE qes_init_status(obj, tagname, status)
   IMPLICIT NONE

   TYPE(status_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: status

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%status = status

END SUBROUTINE qes_init_status

SUBROUTINE qes_reset_status(obj)
   IMPLICIT NONE
   TYPE(status_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_status


SUBROUTINE qes_write_scalarQuantity(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(scalarQuantity_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'Units', TRIM(obj%Units))

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      WRITE(iun, '(E24.16)') obj%scalarQuantity
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_scalarQuantity

SUBROUTINE qes_init_scalarQuantity(obj, tagname, Units, scalarQuantity)
   IMPLICIT NONE

   TYPE(scalarQuantity_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: Units
   REAL(DP) :: scalarQuantity

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%Units = TRIM(Units)

   obj%scalarQuantity = scalarQuantity

END SUBROUTINE qes_init_scalarQuantity

SUBROUTINE qes_reset_scalarQuantity(obj)
   IMPLICIT NONE
   TYPE(scalarQuantity_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_scalarQuantity


SUBROUTINE qes_write_finiteFieldOut(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(finiteFieldOut_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'electronicDipole')
         WRITE(iun, '(3(E24.16))') obj%electronicDipole
      CALL iotk_write_end(iun, 'electronicDipole')
      CALL iotk_write_begin(iun, 'ionicDipole')
         WRITE(iun, '(3(E24.16))') obj%ionicDipole
      CALL iotk_write_end(iun, 'ionicDipole')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_finiteFieldOut

SUBROUTINE qes_init_finiteFieldOut(obj, tagname, electronicDipole, ionicDipole)
   IMPLICIT NONE

   TYPE(finiteFieldOut_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   REAL(DP), DIMENSION(3) :: electronicDipole
   REAL(DP), DIMENSION(3) :: ionicDipole

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%electronicDipole = electronicDipole
   obj%ionicDipole = ionicDipole

END SUBROUTINE qes_init_finiteFieldOut

SUBROUTINE qes_reset_finiteFieldOut(obj)
   IMPLICIT NONE
   TYPE(finiteFieldOut_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_finiteFieldOut


SUBROUTINE qes_write_k_point(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(k_point_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   IF(obj%weight_ispresent) THEN
      CALL iotk_write_attr(attr, 'weight', obj%weight)
   END IF
   IF(obj%label_ispresent) THEN
      CALL iotk_write_attr(attr, 'label', TRIM(obj%label))
   END IF

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      WRITE(iun, '(3(E24.16))') obj%k_point
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_k_point

SUBROUTINE qes_init_k_point(obj, tagname, weight, weight_ispresent, label, label_ispresent, &
                              k_point)
   IMPLICIT NONE

   TYPE(k_point_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   LOGICAL  :: weight_ispresent
   REAL(DP), OPTIONAL :: weight
   LOGICAL  :: label_ispresent
   CHARACTER(len=*), OPTIONAL :: label
   REAL(DP), DIMENSION(3) :: k_point

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%weight_ispresent = weight_ispresent
   IF (obj%weight_ispresent) THEN
      obj%weight = weight
   ENDIF


   obj%label_ispresent = label_ispresent
   IF (obj%label_ispresent) THEN
      obj%label = TRIM(label)
   ENDIF

   obj%k_point = k_point

END SUBROUTINE qes_init_k_point

SUBROUTINE qes_reset_k_point(obj)
   IMPLICIT NONE
   TYPE(k_point_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_k_point


SUBROUTINE qes_write_atom(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(atom_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'name', TRIM(obj%name))
   IF(obj%position_ispresent) THEN
      CALL iotk_write_attr(attr, 'position', TRIM(obj%position))
   END IF
   IF(obj%index_ispresent) THEN
      CALL iotk_write_attr(attr, 'index', obj%index)
   END IF

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      WRITE(iun, '(3(E24.16))') obj%atom
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_atom

SUBROUTINE qes_init_atom(obj, tagname, name, position, position_ispresent, index, index_ispresent, &
                              atom)
   IMPLICIT NONE

   TYPE(atom_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: name
   LOGICAL  :: position_ispresent
   CHARACTER(len=*), OPTIONAL :: position
   LOGICAL  :: index_ispresent
   INTEGER , OPTIONAL :: index
   REAL(DP), DIMENSION(3) :: atom

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%name = TRIM(name)


   obj%position_ispresent = position_ispresent
   IF (obj%position_ispresent) THEN
      obj%position = TRIM(position)
   ENDIF


   obj%index_ispresent = index_ispresent
   IF (obj%index_ispresent) THEN
      obj%index = index
   ENDIF

   obj%atom = atom

END SUBROUTINE qes_init_atom

SUBROUTINE qes_reset_atom(obj)
   IMPLICIT NONE
   TYPE(atom_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_atom


SUBROUTINE qes_write_phase(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(phase_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   IF(obj%ionic_ispresent) THEN
      CALL iotk_write_attr(attr, 'ionic', obj%ionic)
   END IF
   IF(obj%electronic_ispresent) THEN
      CALL iotk_write_attr(attr, 'electronic', obj%electronic)
   END IF
   IF(obj%modulus_ispresent) THEN
      CALL iotk_write_attr(attr, 'modulus', TRIM(obj%modulus))
   END IF

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      WRITE(iun, '(E24.16)') obj%phase
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_phase

SUBROUTINE qes_init_phase(obj, tagname, ionic, ionic_ispresent, electronic, electronic_ispresent, &
                              modulus, modulus_ispresent, phase)
   IMPLICIT NONE

   TYPE(phase_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   LOGICAL  :: ionic_ispresent
   REAL(DP), OPTIONAL :: ionic
   LOGICAL  :: electronic_ispresent
   REAL(DP), OPTIONAL :: electronic
   LOGICAL  :: modulus_ispresent
   CHARACTER(len=*), OPTIONAL :: modulus
   REAL(DP) :: phase

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%ionic_ispresent = ionic_ispresent
   IF (obj%ionic_ispresent) THEN
      obj%ionic = ionic
   ENDIF


   obj%electronic_ispresent = electronic_ispresent
   IF (obj%electronic_ispresent) THEN
      obj%electronic = electronic
   ENDIF


   obj%modulus_ispresent = modulus_ispresent
   IF (obj%modulus_ispresent) THEN
      obj%modulus = TRIM(modulus)
   ENDIF

   obj%phase = phase

END SUBROUTINE qes_init_phase

SUBROUTINE qes_reset_phase(obj)
   IMPLICIT NONE
   TYPE(phase_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_phase


SUBROUTINE qes_write_dipoleOutput(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(dipoleOutput_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'idir')
         WRITE(iun, '(I12)') obj%idir
      CALL iotk_write_end(iun, 'idir')
      CALL qes_write_scalarQuantity(iun, obj%dipole)
      !
      CALL qes_write_scalarQuantity(iun, obj%ion_dipole)
      !
      CALL qes_write_scalarQuantity(iun, obj%elec_dipole)
      !
      CALL qes_write_scalarQuantity(iun, obj%dipoleField)
      !
      CALL qes_write_scalarQuantity(iun, obj%potentialAmp)
      !
      CALL qes_write_scalarQuantity(iun, obj%totalLength)
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_dipoleOutput

SUBROUTINE qes_init_dipoleOutput(obj, tagname, idir, dipole, ion_dipole, elec_dipole, &
                              dipoleField, potentialAmp, totalLength)
   IMPLICIT NONE

   TYPE(dipoleOutput_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: idir
   TYPE(scalarQuantity_type) :: dipole
   TYPE(scalarQuantity_type) :: ion_dipole
   TYPE(scalarQuantity_type) :: elec_dipole
   TYPE(scalarQuantity_type) :: dipoleField
   TYPE(scalarQuantity_type) :: potentialAmp
   TYPE(scalarQuantity_type) :: totalLength

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%idir = idir
   obj%dipole = dipole
   obj%ion_dipole = ion_dipole
   obj%elec_dipole = elec_dipole
   obj%dipoleField = dipoleField
   obj%potentialAmp = potentialAmp
   obj%totalLength = totalLength

END SUBROUTINE qes_init_dipoleOutput

SUBROUTINE qes_reset_dipoleOutput(obj)
   IMPLICIT NONE
   TYPE(dipoleOutput_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   CALL qes_reset_scalarQuantity(obj%dipole)
   CALL qes_reset_scalarQuantity(obj%ion_dipole)
   CALL qes_reset_scalarQuantity(obj%elec_dipole)
   CALL qes_reset_scalarQuantity(obj%dipoleField)
   CALL qes_reset_scalarQuantity(obj%potentialAmp)
   CALL qes_reset_scalarQuantity(obj%totalLength)

END SUBROUTINE qes_reset_dipoleOutput


SUBROUTINE qes_write_polarization(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(polarization_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL qes_write_scalarQuantity(iun, obj%polarization)
      !
      CALL iotk_write_begin(iun, 'modulus')
         WRITE(iun, '(E24.16)') obj%modulus
      CALL iotk_write_end(iun, 'modulus')
      CALL iotk_write_begin(iun, 'direction')
         WRITE(iun, '(3(E24.16))') obj%direction
      CALL iotk_write_end(iun, 'direction')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_polarization

SUBROUTINE qes_init_polarization(obj, tagname, polarization, modulus, direction)
   IMPLICIT NONE

   TYPE(polarization_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   TYPE(scalarQuantity_type) :: polarization
   REAL(DP) :: modulus
   REAL(DP), DIMENSION(3) :: direction

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%polarization = polarization
   obj%modulus = modulus
   obj%direction = direction

END SUBROUTINE qes_init_polarization

SUBROUTINE qes_reset_polarization(obj)
   IMPLICIT NONE
   TYPE(polarization_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   CALL qes_reset_scalarQuantity(obj%polarization)

END SUBROUTINE qes_reset_polarization


SUBROUTINE qes_write_ionicPolarization(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(ionicPolarization_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL qes_write_atom(iun, obj%ion)
      !
      CALL iotk_write_begin(iun, 'charge')
         WRITE(iun, '(E24.16)') obj%charge
      CALL iotk_write_end(iun, 'charge')
      CALL qes_write_phase(iun, obj%phase)
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_ionicPolarization

SUBROUTINE qes_init_ionicPolarization(obj, tagname, ion, charge, phase)
   IMPLICIT NONE

   TYPE(ionicPolarization_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   TYPE(atom_type) :: ion
   REAL(DP) :: charge
   TYPE(phase_type) :: phase

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%ion = ion
   obj%charge = charge
   obj%phase = phase

END SUBROUTINE qes_init_ionicPolarization

SUBROUTINE qes_reset_ionicPolarization(obj)
   IMPLICIT NONE
   TYPE(ionicPolarization_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   CALL qes_reset_atom(obj%ion)
   CALL qes_reset_phase(obj%phase)

END SUBROUTINE qes_reset_ionicPolarization


SUBROUTINE qes_write_electronicPolarization(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(electronicPolarization_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL qes_write_k_point(iun, obj%firstKeyPoint)
      !
      IF(obj%spin_ispresent) THEN
         CALL iotk_write_begin(iun, 'spin')
            WRITE(iun, '(I12)') obj%spin
         CALL iotk_write_end(iun, 'spin')
      ENDIF
      !
      CALL qes_write_phase(iun, obj%phase)
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_electronicPolarization

SUBROUTINE qes_init_electronicPolarization(obj, tagname, firstKeyPoint, spin_ispresent, &
                              spin, phase)
   IMPLICIT NONE

   TYPE(electronicPolarization_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   TYPE(k_point_type) :: firstKeyPoint
   LOGICAL  :: spin_ispresent
   INTEGER  :: spin
   TYPE(phase_type) :: phase

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%firstKeyPoint = firstKeyPoint
   obj%spin_ispresent = spin_ispresent
   IF(obj%spin_ispresent) THEN
      obj%spin = spin
   ENDIF
   obj%phase = phase

END SUBROUTINE qes_init_electronicPolarization

SUBROUTINE qes_reset_electronicPolarization(obj)
   IMPLICIT NONE
   TYPE(electronicPolarization_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   CALL qes_reset_k_point(obj%firstKeyPoint)
   IF(obj%spin_ispresent) THEN
      obj%spin_ispresent = .FALSE.
   ENDIF
   CALL qes_reset_phase(obj%phase)

END SUBROUTINE qes_reset_electronicPolarization


SUBROUTINE qes_write_BerryPhaseOutput(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(BerryPhaseOutput_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL qes_write_polarization(iun, obj%polarization)
      !
      CALL qes_write_phase(iun, obj%totalPhase)
      !
      DO i = 1, obj%ndim_ionicPolarization
         CALL qes_write_ionicPolarization(iun, obj%ionicPolarization(i))
         !
      END DO
      DO i = 1, obj%ndim_electronicPolarization
         CALL qes_write_electronicPolarization(iun, obj%electronicPolarization(i))
         !
      END DO
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_BerryPhaseOutput

SUBROUTINE qes_init_BerryPhaseOutput(obj, tagname, polarization, totalPhase, &
                              ndim_ionicPolarization, ionicPolarization, &
                              ndim_electronicPolarization, electronicPolarization)
   IMPLICIT NONE

   TYPE(BerryPhaseOutput_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   TYPE(polarization_type) :: polarization
   TYPE(phase_type) :: totalPhase
   INTEGER  :: ndim_ionicPolarization
   TYPE(ionicPolarization_type ), DIMENSION( ndim_ionicPolarization )  :: ionicPolarization
   INTEGER  :: ndim_electronicPolarization
   TYPE(electronicPolarization_type ), DIMENSION( ndim_electronicPolarization )  :: electronicPolarization

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%polarization = polarization
   obj%totalPhase = totalPhase
   ALLOCATE(obj%ionicPolarization(SIZE(ionicPolarization)))
   DO i = 1, SIZE(ionicPolarization)
      obj%ionicPolarization(i) = ionicPolarization(i)
   ENDDO
   obj%ndim_ionicPolarization = ndim_ionicPolarization
   ALLOCATE(obj%electronicPolarization(SIZE(electronicPolarization)))
   DO i = 1, SIZE(electronicPolarization)
      obj%electronicPolarization(i) = electronicPolarization(i)
   ENDDO
   obj%ndim_electronicPolarization = ndim_electronicPolarization

END SUBROUTINE qes_init_BerryPhaseOutput

SUBROUTINE qes_reset_BerryPhaseOutput(obj)
   IMPLICIT NONE
   TYPE(BerryPhaseOutput_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   CALL qes_reset_polarization(obj%polarization)
   CALL qes_reset_phase(obj%totalPhase)
   DO i = 1, SIZE(obj%ionicPolarization)
      CALL qes_reset_ionicPolarization(obj%ionicPolarization(i))
   ENDDO
   IF (ALLOCATED(obj%ionicPolarization)) DEALLOCATE(obj%ionicPolarization)
   DO i = 1, SIZE(obj%electronicPolarization)
      CALL qes_reset_electronicPolarization(obj%electronicPolarization(i))
   ENDDO
   IF (ALLOCATED(obj%electronicPolarization)) DEALLOCATE(obj%electronicPolarization)

END SUBROUTINE qes_reset_BerryPhaseOutput


SUBROUTINE qes_write_vector(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(vector_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun,'vec')
         WRITE(fmtstr,'(a)') '(5E24.16)'
         WRITE(iun, fmtstr) obj%vec
      CALL iotk_write_end(iun,'vec')
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_vector

SUBROUTINE qes_init_vector(obj, tagname, ndim_vec, vec)
   IMPLICIT NONE

   TYPE(vector_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: ndim_vec
   REAL(DP), DIMENSION(ndim_vec) :: vec

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   ALLOCATE(obj%vec(ndim_vec))
   obj%vec(:) = vec(:)
   obj%ndim_vec = ndim_vec

END SUBROUTINE qes_init_vector

SUBROUTINE qes_reset_vector(obj)
   IMPLICIT NONE
   TYPE(vector_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF (ALLOCATED(obj%vec))  DEALLOCATE(obj%vec)

END SUBROUTINE qes_reset_vector


SUBROUTINE qes_write_ks_energies(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(ks_energies_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL qes_write_k_point(iun, obj%k_point)
      !
      CALL iotk_write_begin(iun, 'npw')
         WRITE(iun, '(I12)') obj%npw
      CALL iotk_write_end(iun, 'npw')
      CALL iotk_write_begin(iun,'eigenvalues')
         WRITE(fmtstr,'(a)') '(5E24.16)'
         WRITE(iun, fmtstr) obj%eigenvalues
      CALL iotk_write_end(iun,'eigenvalues')
      !
      CALL iotk_write_begin(iun,'occupations')
         WRITE(fmtstr,'(a)') '(5E24.16)'
         WRITE(iun, fmtstr) obj%occupations
      CALL iotk_write_end(iun,'occupations')
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_ks_energies

SUBROUTINE qes_init_ks_energies(obj, tagname, k_point, npw, ndim_eigenvalues, eigenvalues, &
                              ndim_occupations, occupations)
   IMPLICIT NONE

   TYPE(ks_energies_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   TYPE(k_point_type) :: k_point
   INTEGER  :: npw
   INTEGER  :: ndim_eigenvalues
   REAL(DP), DIMENSION(ndim_eigenvalues) :: eigenvalues
   INTEGER  :: ndim_occupations
   REAL(DP), DIMENSION(ndim_occupations) :: occupations

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%k_point = k_point
   obj%npw = npw
   ALLOCATE(obj%eigenvalues(ndim_eigenvalues))
   obj%eigenvalues(:) = eigenvalues(:)
   obj%ndim_eigenvalues = ndim_eigenvalues
   ALLOCATE(obj%occupations(ndim_occupations))
   obj%occupations(:) = occupations(:)
   obj%ndim_occupations = ndim_occupations

END SUBROUTINE qes_init_ks_energies

SUBROUTINE qes_reset_ks_energies(obj)
   IMPLICIT NONE
   TYPE(ks_energies_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   CALL qes_reset_k_point(obj%k_point)
   IF (ALLOCATED(obj%eigenvalues))  DEALLOCATE(obj%eigenvalues)
   IF (ALLOCATED(obj%occupations))  DEALLOCATE(obj%occupations)

END SUBROUTINE qes_reset_ks_energies


SUBROUTINE qes_write_magnetization(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(magnetization_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'lsda',new_line=.FALSE.)
         IF (obj%lsda) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'lsda',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'noncolin',new_line=.FALSE.)
         IF (obj%noncolin) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'noncolin',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'spinorbit',new_line=.FALSE.)
         IF (obj%spinorbit) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'spinorbit',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'total')
         WRITE(iun, '(E24.16)') obj%total
      CALL iotk_write_end(iun, 'total')
      CALL iotk_write_begin(iun, 'absolute')
         WRITE(iun, '(E24.16)') obj%absolute
      CALL iotk_write_end(iun, 'absolute')
      CALL iotk_write_begin(iun, 'do_magnetization',new_line=.FALSE.)
         IF (obj%do_magnetization) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'do_magnetization',indentation=.FALSE.)
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_magnetization

SUBROUTINE qes_init_magnetization(obj, tagname, lsda, noncolin, spinorbit, total, absolute, &
                              do_magnetization)
   IMPLICIT NONE

   TYPE(magnetization_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   LOGICAL  :: lsda
   LOGICAL  :: noncolin
   LOGICAL  :: spinorbit
   REAL(DP) :: total
   REAL(DP) :: absolute
   LOGICAL  :: do_magnetization

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%lsda = lsda
   obj%noncolin = noncolin
   obj%spinorbit = spinorbit
   obj%total = total
   obj%absolute = absolute
   obj%do_magnetization = do_magnetization

END SUBROUTINE qes_init_magnetization

SUBROUTINE qes_reset_magnetization(obj)
   IMPLICIT NONE
   TYPE(magnetization_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_magnetization


SUBROUTINE qes_write_reciprocal_lattice(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(reciprocal_lattice_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'b1')
         WRITE(iun, '(3(E24.16))') obj%b1
      CALL iotk_write_end(iun, 'b1')
      CALL iotk_write_begin(iun, 'b2')
         WRITE(iun, '(3(E24.16))') obj%b2
      CALL iotk_write_end(iun, 'b2')
      CALL iotk_write_begin(iun, 'b3')
         WRITE(iun, '(3(E24.16))') obj%b3
      CALL iotk_write_end(iun, 'b3')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_reciprocal_lattice

SUBROUTINE qes_init_reciprocal_lattice(obj, tagname, b1, b2, b3)
   IMPLICIT NONE

   TYPE(reciprocal_lattice_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   REAL(DP), DIMENSION(3) :: b1
   REAL(DP), DIMENSION(3) :: b2
   REAL(DP), DIMENSION(3) :: b3

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%b1 = b1
   obj%b2 = b2
   obj%b3 = b3

END SUBROUTINE qes_init_reciprocal_lattice

SUBROUTINE qes_reset_reciprocal_lattice(obj)
   IMPLICIT NONE
   TYPE(reciprocal_lattice_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_reciprocal_lattice


SUBROUTINE qes_write_basisSetItem(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(basisSetItem_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'nr1', obj%nr1)
   CALL iotk_write_attr(attr, 'nr2', obj%nr2)
   CALL iotk_write_attr(attr, 'nr3', obj%nr3)

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr),new_line=.FALSE.)
      !
      WRITE(iun, '(A)',advance='no')  TRIM(obj%basisSetItem)
   CALL iotk_write_end(iun, TRIM(obj%tagname),indentation=.FALSE.)
   !
END SUBROUTINE qes_write_basisSetItem

SUBROUTINE qes_init_basisSetItem(obj, tagname, nr1, nr2, nr3, basisSetItem)
   IMPLICIT NONE

   TYPE(basisSetItem_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: nr1
   INTEGER  :: nr2
   INTEGER  :: nr3
   CHARACTER(len=*) :: basisSetItem

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%nr1 = nr1


   obj%nr2 = nr2


   obj%nr3 = nr3

   obj%basisSetItem = basisSetItem

END SUBROUTINE qes_init_basisSetItem

SUBROUTINE qes_reset_basisSetItem(obj)
   IMPLICIT NONE
   TYPE(basisSetItem_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_basisSetItem


SUBROUTINE qes_write_equivalent_atoms(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(equivalent_atoms_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'nat', obj%nat)

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
         WRITE(fmtstr,'(a)') '(12I6)'
         WRITE(iun, fmtstr) obj%index_list
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_equivalent_atoms

SUBROUTINE qes_init_equivalent_atoms(obj, tagname, nat, ndim_index_list, index_list)
   IMPLICIT NONE

   TYPE(equivalent_atoms_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: nat
   INTEGER  :: ndim_index_list
   INTEGER, DIMENSION(ndim_index_list) :: index_list

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%nat = nat

   ALLOCATE(obj%index_list(ndim_index_list))
   obj%index_list(:) = index_list(:)
   obj%ndim_index_list = ndim_index_list

END SUBROUTINE qes_init_equivalent_atoms

SUBROUTINE qes_reset_equivalent_atoms(obj)
   IMPLICIT NONE
   TYPE(equivalent_atoms_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF (ALLOCATED(obj%index_list))  DEALLOCATE(obj%index_list)

END SUBROUTINE qes_reset_equivalent_atoms


SUBROUTINE qes_write_info(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(info_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   IF(obj%name_ispresent) THEN
      CALL iotk_write_attr(attr, 'name', TRIM(obj%name))
   END IF
   IF(obj%class_ispresent) THEN
      CALL iotk_write_attr(attr, 'class', TRIM(obj%class))
   END IF
   IF(obj%time_reversal_ispresent) THEN
      IF (obj%time_reversal) THEN
         WRITE (fmtstr,'("true")')
      ELSE
         WRITE (fmtstr,'("false")')
END IF
      CALL iotk_write_attr(attr, 'time_reversal', TRIM(fmtstr) )
   END IF

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr),new_line=.FALSE.)
      !
      WRITE(iun, '(A)',advance='no')  TRIM(obj%info)
   CALL iotk_write_end(iun, TRIM(obj%tagname),indentation=.FALSE.)
   !
END SUBROUTINE qes_write_info

SUBROUTINE qes_init_info(obj, tagname, name, name_ispresent, class, class_ispresent, &
                              time_reversal, time_reversal_ispresent, info)
   IMPLICIT NONE

   TYPE(info_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   LOGICAL  :: name_ispresent
   CHARACTER(len=*), OPTIONAL :: name
   LOGICAL  :: class_ispresent
   CHARACTER(len=*), OPTIONAL :: class
   LOGICAL  :: time_reversal_ispresent
   LOGICAL , OPTIONAL :: time_reversal
   CHARACTER(len=*) :: info

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%name_ispresent = name_ispresent
   IF (obj%name_ispresent) THEN
      obj%name = TRIM(name)
   ENDIF


   obj%class_ispresent = class_ispresent
   IF (obj%class_ispresent) THEN
      obj%class = TRIM(class)
   ENDIF


   obj%time_reversal_ispresent = time_reversal_ispresent
   IF (obj%time_reversal_ispresent) THEN
      obj%time_reversal = time_reversal
   ENDIF

   obj%info = info

END SUBROUTINE qes_init_info

SUBROUTINE qes_reset_info(obj)
   IMPLICIT NONE
   TYPE(info_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_info


SUBROUTINE qes_write_matrix(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(matrix_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      WRITE(fmtstr,'(a)') '(5E24.16)'
      DO i = 1, SIZE(obj%mat,2)
         WRITE(iun, fmtstr) obj%mat(:,i)
      ENDDO
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_matrix

SUBROUTINE qes_init_matrix(obj, tagname, ndim1_mat, ndim2_mat, mat)
   IMPLICIT NONE

   TYPE(matrix_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: ndim1_mat
   INTEGER  :: ndim2_mat
   REAL(DP), DIMENSION(ndim1_mat,ndim2_mat) :: mat

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   ALLOCATE(obj%mat(ndim1_mat,ndim2_mat))
   obj%mat(:,:) = mat(:,:)
   obj%ndim1_mat = ndim1_mat
   obj%ndim2_mat = ndim2_mat

END SUBROUTINE qes_init_matrix

SUBROUTINE qes_reset_matrix(obj)
   IMPLICIT NONE
   TYPE(matrix_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF (ALLOCATED(obj%mat))  DEALLOCATE(obj%mat)

END SUBROUTINE qes_reset_matrix


SUBROUTINE qes_write_symmetry(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(symmetry_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL qes_write_info(iun, obj%info)
      !
      CALL qes_write_matrix(iun, obj%rotation)
      !
      IF(obj%fractional_translation_ispresent) THEN
         CALL iotk_write_begin(iun, 'fractional_translation')
            WRITE(iun, '(3(E24.16))') obj%fractional_translation
         CALL iotk_write_end(iun, 'fractional_translation')
      ENDIF
      !
      IF(obj%equivalent_atoms_ispresent) THEN
         CALL qes_write_equivalent_atoms(iun, obj%equivalent_atoms)
         !
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_symmetry

SUBROUTINE qes_init_symmetry(obj, tagname, info, rotation, &
                              fractional_translation_ispresent, fractional_translation, &
                              equivalent_atoms_ispresent, equivalent_atoms)
   IMPLICIT NONE

   TYPE(symmetry_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   TYPE(info_type) :: info
   TYPE(matrix_type) :: rotation
   LOGICAL  :: fractional_translation_ispresent
   REAL(DP), DIMENSION(3) :: fractional_translation
   LOGICAL  :: equivalent_atoms_ispresent
   TYPE(equivalent_atoms_type) :: equivalent_atoms

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%info = info
   obj%rotation = rotation
   obj%fractional_translation_ispresent = fractional_translation_ispresent
   IF(obj%fractional_translation_ispresent) THEN
      obj%fractional_translation = fractional_translation
   ENDIF
   obj%equivalent_atoms_ispresent = equivalent_atoms_ispresent
   IF(obj%equivalent_atoms_ispresent) THEN
      obj%equivalent_atoms = equivalent_atoms
   ENDIF

END SUBROUTINE qes_init_symmetry

SUBROUTINE qes_reset_symmetry(obj)
   IMPLICIT NONE
   TYPE(symmetry_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   CALL qes_reset_info(obj%info)
   CALL qes_reset_matrix(obj%rotation)
   IF(obj%fractional_translation_ispresent) THEN
      obj%fractional_translation_ispresent = .FALSE.
   ENDIF
   IF(obj%equivalent_atoms_ispresent) THEN
      CALL qes_reset_equivalent_atoms(obj%equivalent_atoms)
      obj%equivalent_atoms_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_symmetry


SUBROUTINE qes_write_algorithmic_info(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(algorithmic_info_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'real_space_q',new_line=.FALSE.)
         IF (obj%real_space_q) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'real_space_q',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'uspp',new_line=.FALSE.)
         IF (obj%uspp) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'uspp',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'paw',new_line=.FALSE.)
         IF (obj%paw) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'paw',indentation=.FALSE.)
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_algorithmic_info

SUBROUTINE qes_init_algorithmic_info(obj, tagname, real_space_q, uspp, paw)
   IMPLICIT NONE

   TYPE(algorithmic_info_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   LOGICAL  :: real_space_q
   LOGICAL  :: uspp
   LOGICAL  :: paw

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%real_space_q = real_space_q
   obj%uspp = uspp
   obj%paw = paw

END SUBROUTINE qes_init_algorithmic_info

SUBROUTINE qes_reset_algorithmic_info(obj)
   IMPLICIT NONE
   TYPE(algorithmic_info_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_algorithmic_info


SUBROUTINE qes_write_opt_conv(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(opt_conv_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'n_opt_steps')
         WRITE(iun, '(I12)') obj%n_opt_steps
      CALL iotk_write_end(iun, 'n_opt_steps')
      CALL iotk_write_begin(iun, 'grad_norm')
         WRITE(iun, '(E24.16)') obj%grad_norm
      CALL iotk_write_end(iun, 'grad_norm')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_opt_conv

SUBROUTINE qes_init_opt_conv(obj, tagname, n_opt_steps, grad_norm)
   IMPLICIT NONE

   TYPE(opt_conv_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: n_opt_steps
   REAL(DP) :: grad_norm

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%n_opt_steps = n_opt_steps
   obj%grad_norm = grad_norm

END SUBROUTINE qes_init_opt_conv

SUBROUTINE qes_reset_opt_conv(obj)
   IMPLICIT NONE
   TYPE(opt_conv_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_opt_conv


SUBROUTINE qes_write_scf_conv(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(scf_conv_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'n_scf_steps')
         WRITE(iun, '(I12)') obj%n_scf_steps
      CALL iotk_write_end(iun, 'n_scf_steps')
      CALL iotk_write_begin(iun, 'scf_error')
         WRITE(iun, '(E24.16)') obj%scf_error
      CALL iotk_write_end(iun, 'scf_error')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_scf_conv

SUBROUTINE qes_init_scf_conv(obj, tagname, n_scf_steps, scf_error)
   IMPLICIT NONE

   TYPE(scf_conv_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: n_scf_steps
   REAL(DP) :: scf_error

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%n_scf_steps = n_scf_steps
   obj%scf_error = scf_error

END SUBROUTINE qes_init_scf_conv

SUBROUTINE qes_reset_scf_conv(obj)
   IMPLICIT NONE
   TYPE(scf_conv_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_scf_conv


SUBROUTINE qes_write_species(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(species_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'name', TRIM(obj%name))

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      IF(obj%mass_ispresent) THEN
         CALL iotk_write_begin(iun, 'mass')
            WRITE(iun, '(E24.16)') obj%mass
         CALL iotk_write_end(iun, 'mass')
      ENDIF
      !
      CALL iotk_write_begin(iun, 'pseudo_file',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%pseudo_file)
      CALL iotk_write_end(iun, 'pseudo_file',indentation=.FALSE.)
      IF(obj%starting_magnetization_ispresent) THEN
         CALL iotk_write_begin(iun, 'starting_magnetization')
            WRITE(iun, '(E24.16)') obj%starting_magnetization
         CALL iotk_write_end(iun, 'starting_magnetization')
      ENDIF
      !
      IF(obj%spin_teta_ispresent) THEN
         CALL iotk_write_begin(iun, 'spin_teta')
            WRITE(iun, '(E24.16)') obj%spin_teta
         CALL iotk_write_end(iun, 'spin_teta')
      ENDIF
      !
      IF(obj%spin_phi_ispresent) THEN
         CALL iotk_write_begin(iun, 'spin_phi')
            WRITE(iun, '(E24.16)') obj%spin_phi
         CALL iotk_write_end(iun, 'spin_phi')
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_species

SUBROUTINE qes_init_species(obj, tagname, name, mass_ispresent, mass, pseudo_file, &
                              starting_magnetization_ispresent, starting_magnetization, &
                              spin_teta_ispresent, spin_teta, spin_phi_ispresent, spin_phi)
   IMPLICIT NONE

   TYPE(species_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: name
   LOGICAL  :: mass_ispresent
   REAL(DP) :: mass
   CHARACTER(len=*) :: pseudo_file
   LOGICAL  :: starting_magnetization_ispresent
   REAL(DP) :: starting_magnetization
   LOGICAL  :: spin_teta_ispresent
   REAL(DP) :: spin_teta
   LOGICAL  :: spin_phi_ispresent
   REAL(DP) :: spin_phi

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%name = TRIM(name)

   obj%mass_ispresent = mass_ispresent
   IF(obj%mass_ispresent) THEN
      obj%mass = mass
   ENDIF
   obj%pseudo_file = pseudo_file
   obj%starting_magnetization_ispresent = starting_magnetization_ispresent
   IF(obj%starting_magnetization_ispresent) THEN
      obj%starting_magnetization = starting_magnetization
   ENDIF
   obj%spin_teta_ispresent = spin_teta_ispresent
   IF(obj%spin_teta_ispresent) THEN
      obj%spin_teta = spin_teta
   ENDIF
   obj%spin_phi_ispresent = spin_phi_ispresent
   IF(obj%spin_phi_ispresent) THEN
      obj%spin_phi = spin_phi
   ENDIF

END SUBROUTINE qes_init_species

SUBROUTINE qes_reset_species(obj)
   IMPLICIT NONE
   TYPE(species_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%mass_ispresent) THEN
      obj%mass_ispresent = .FALSE.
   ENDIF
   IF(obj%starting_magnetization_ispresent) THEN
      obj%starting_magnetization_ispresent = .FALSE.
   ENDIF
   IF(obj%spin_teta_ispresent) THEN
      obj%spin_teta_ispresent = .FALSE.
   ENDIF
   IF(obj%spin_phi_ispresent) THEN
      obj%spin_phi_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_species


SUBROUTINE qes_write_total_energy(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(total_energy_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'etot')
         WRITE(iun, '(E24.16)') obj%etot
      CALL iotk_write_end(iun, 'etot')
      IF(obj%eband_ispresent) THEN
         CALL iotk_write_begin(iun, 'eband')
            WRITE(iun, '(E24.16)') obj%eband
         CALL iotk_write_end(iun, 'eband')
      ENDIF
      !
      IF(obj%ehart_ispresent) THEN
         CALL iotk_write_begin(iun, 'ehart')
            WRITE(iun, '(E24.16)') obj%ehart
         CALL iotk_write_end(iun, 'ehart')
      ENDIF
      !
      IF(obj%vtxc_ispresent) THEN
         CALL iotk_write_begin(iun, 'vtxc')
            WRITE(iun, '(E24.16)') obj%vtxc
         CALL iotk_write_end(iun, 'vtxc')
      ENDIF
      !
      IF(obj%etxc_ispresent) THEN
         CALL iotk_write_begin(iun, 'etxc')
            WRITE(iun, '(E24.16)') obj%etxc
         CALL iotk_write_end(iun, 'etxc')
      ENDIF
      !
      IF(obj%ewald_ispresent) THEN
         CALL iotk_write_begin(iun, 'ewald')
            WRITE(iun, '(E24.16)') obj%ewald
         CALL iotk_write_end(iun, 'ewald')
      ENDIF
      !
      IF(obj%demet_ispresent) THEN
         CALL iotk_write_begin(iun, 'demet')
            WRITE(iun, '(E24.16)') obj%demet
         CALL iotk_write_end(iun, 'demet')
      ENDIF
      !
      IF(obj%efieldcorr_ispresent) THEN
         CALL iotk_write_begin(iun, 'efieldcorr')
            WRITE(iun, '(E24.16)') obj%efieldcorr
         CALL iotk_write_end(iun, 'efieldcorr')
      ENDIF
      !
      IF(obj%potentiostat_contr_ispresent) THEN
         CALL iotk_write_begin(iun, 'potentiostat_contr')
            WRITE(iun, '(E24.16)') obj%potentiostat_contr
         CALL iotk_write_end(iun, 'potentiostat_contr')
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_total_energy

SUBROUTINE qes_init_total_energy(obj, tagname, etot, eband_ispresent, eband, &
                              ehart_ispresent, ehart, vtxc_ispresent, vtxc, etxc_ispresent, &
                              etxc, ewald_ispresent, ewald, demet_ispresent, demet, &
                              efieldcorr_ispresent, efieldcorr, &
                              potentiostat_contr_ispresent, potentiostat_contr)
   IMPLICIT NONE

   TYPE(total_energy_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
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

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%etot = etot
   obj%eband_ispresent = eband_ispresent
   IF(obj%eband_ispresent) THEN
      obj%eband = eband
   ENDIF
   obj%ehart_ispresent = ehart_ispresent
   IF(obj%ehart_ispresent) THEN
      obj%ehart = ehart
   ENDIF
   obj%vtxc_ispresent = vtxc_ispresent
   IF(obj%vtxc_ispresent) THEN
      obj%vtxc = vtxc
   ENDIF
   obj%etxc_ispresent = etxc_ispresent
   IF(obj%etxc_ispresent) THEN
      obj%etxc = etxc
   ENDIF
   obj%ewald_ispresent = ewald_ispresent
   IF(obj%ewald_ispresent) THEN
      obj%ewald = ewald
   ENDIF
   obj%demet_ispresent = demet_ispresent
   IF(obj%demet_ispresent) THEN
      obj%demet = demet
   ENDIF
   obj%efieldcorr_ispresent = efieldcorr_ispresent
   IF(obj%efieldcorr_ispresent) THEN
      obj%efieldcorr = efieldcorr
   ENDIF
   obj%potentiostat_contr_ispresent = potentiostat_contr_ispresent
   IF(obj%potentiostat_contr_ispresent) THEN
      obj%potentiostat_contr = potentiostat_contr
   ENDIF

END SUBROUTINE qes_init_total_energy

SUBROUTINE qes_reset_total_energy(obj)
   IMPLICIT NONE
   TYPE(total_energy_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%eband_ispresent) THEN
      obj%eband_ispresent = .FALSE.
   ENDIF
   IF(obj%ehart_ispresent) THEN
      obj%ehart_ispresent = .FALSE.
   ENDIF
   IF(obj%vtxc_ispresent) THEN
      obj%vtxc_ispresent = .FALSE.
   ENDIF
   IF(obj%etxc_ispresent) THEN
      obj%etxc_ispresent = .FALSE.
   ENDIF
   IF(obj%ewald_ispresent) THEN
      obj%ewald_ispresent = .FALSE.
   ENDIF
   IF(obj%demet_ispresent) THEN
      obj%demet_ispresent = .FALSE.
   ENDIF
   IF(obj%efieldcorr_ispresent) THEN
      obj%efieldcorr_ispresent = .FALSE.
   ENDIF
   IF(obj%potentiostat_contr_ispresent) THEN
      obj%potentiostat_contr_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_total_energy


SUBROUTINE qes_write_convergence_info(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(convergence_info_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL qes_write_scf_conv(iun, obj%scf_conv)
      !
      IF(obj%opt_conv_ispresent) THEN
         CALL qes_write_opt_conv(iun, obj%opt_conv)
         !
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_convergence_info

SUBROUTINE qes_init_convergence_info(obj, tagname, scf_conv, opt_conv_ispresent, opt_conv)
   IMPLICIT NONE

   TYPE(convergence_info_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   TYPE(scf_conv_type) :: scf_conv
   LOGICAL  :: opt_conv_ispresent
   TYPE(opt_conv_type) :: opt_conv

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%scf_conv = scf_conv
   obj%opt_conv_ispresent = opt_conv_ispresent
   IF(obj%opt_conv_ispresent) THEN
      obj%opt_conv = opt_conv
   ENDIF

END SUBROUTINE qes_init_convergence_info

SUBROUTINE qes_reset_convergence_info(obj)
   IMPLICIT NONE
   TYPE(convergence_info_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   CALL qes_reset_scf_conv(obj%scf_conv)
   IF(obj%opt_conv_ispresent) THEN
      CALL qes_reset_opt_conv(obj%opt_conv)
      obj%opt_conv_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_convergence_info


SUBROUTINE qes_write_outputElectricField(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(outputElectricField_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      IF(obj%BerryPhase_ispresent) THEN
         CALL qes_write_BerryPhaseOutput(iun, obj%BerryPhase)
         !
      ENDIF
      !
      IF(obj%finiteElectricFieldInfo_ispresent) THEN
         CALL qes_write_finiteFieldOut(iun, obj%finiteElectricFieldInfo)
         !
      ENDIF
      !
      IF(obj%dipoleInfo_ispresent) THEN
         CALL qes_write_dipoleOutput(iun, obj%dipoleInfo)
         !
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_outputElectricField

SUBROUTINE qes_init_outputElectricField(obj, tagname, BerryPhase_ispresent, BerryPhase, &
                              finiteElectricFieldInfo_ispresent, finiteElectricFieldInfo, &
                              dipoleInfo_ispresent, dipoleInfo)
   IMPLICIT NONE

   TYPE(outputElectricField_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   LOGICAL  :: BerryPhase_ispresent
   TYPE(BerryPhaseOutput_type) :: BerryPhase
   LOGICAL  :: finiteElectricFieldInfo_ispresent
   TYPE(finiteFieldOut_type) :: finiteElectricFieldInfo
   LOGICAL  :: dipoleInfo_ispresent
   TYPE(dipoleOutput_type) :: dipoleInfo

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%BerryPhase_ispresent = BerryPhase_ispresent
   IF(obj%BerryPhase_ispresent) THEN
      obj%BerryPhase = BerryPhase
   ENDIF
   obj%finiteElectricFieldInfo_ispresent = finiteElectricFieldInfo_ispresent
   IF(obj%finiteElectricFieldInfo_ispresent) THEN
      obj%finiteElectricFieldInfo = finiteElectricFieldInfo
   ENDIF
   obj%dipoleInfo_ispresent = dipoleInfo_ispresent
   IF(obj%dipoleInfo_ispresent) THEN
      obj%dipoleInfo = dipoleInfo
   ENDIF

END SUBROUTINE qes_init_outputElectricField

SUBROUTINE qes_reset_outputElectricField(obj)
   IMPLICIT NONE
   TYPE(outputElectricField_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%BerryPhase_ispresent) THEN
      CALL qes_reset_BerryPhaseOutput(obj%BerryPhase)
      obj%BerryPhase_ispresent = .FALSE.
   ENDIF
   IF(obj%finiteElectricFieldInfo_ispresent) THEN
      CALL qes_reset_finiteFieldOut(obj%finiteElectricFieldInfo)
      obj%finiteElectricFieldInfo_ispresent = .FALSE.
   ENDIF
   IF(obj%dipoleInfo_ispresent) THEN
      CALL qes_reset_dipoleOutput(obj%dipoleInfo)
      obj%dipoleInfo_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_outputElectricField


SUBROUTINE qes_write_spin_constraints(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(spin_constraints_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'spin_constraints',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%spin_constraints)
      CALL iotk_write_end(iun, 'spin_constraints',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'lagrange_multiplier')
         WRITE(iun, '(E24.16)') obj%lagrange_multiplier
      CALL iotk_write_end(iun, 'lagrange_multiplier')
      IF(obj%target_magnetization_ispresent) THEN
         CALL iotk_write_begin(iun, 'target_magnetization')
            WRITE(iun, '(3(E24.16))') obj%target_magnetization
         CALL iotk_write_end(iun, 'target_magnetization')
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_spin_constraints

SUBROUTINE qes_init_spin_constraints(obj, tagname, spin_constraints, lagrange_multiplier, &
                              target_magnetization_ispresent, target_magnetization)
   IMPLICIT NONE

   TYPE(spin_constraints_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: spin_constraints
   REAL(DP) :: lagrange_multiplier
   LOGICAL  :: target_magnetization_ispresent
   REAL(DP), DIMENSION(3) :: target_magnetization

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%spin_constraints = spin_constraints
   obj%lagrange_multiplier = lagrange_multiplier
   obj%target_magnetization_ispresent = target_magnetization_ispresent
   IF(obj%target_magnetization_ispresent) THEN
      obj%target_magnetization = target_magnetization
   ENDIF

END SUBROUTINE qes_init_spin_constraints

SUBROUTINE qes_reset_spin_constraints(obj)
   IMPLICIT NONE
   TYPE(spin_constraints_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%target_magnetization_ispresent) THEN
      obj%target_magnetization_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_spin_constraints


SUBROUTINE qes_write_constr_type(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(constr_type_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_constr_type

SUBROUTINE qes_init_constr_type(obj, tagname)
   IMPLICIT NONE

   TYPE(constr_type_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

END SUBROUTINE qes_init_constr_type

SUBROUTINE qes_reset_constr_type(obj)
   IMPLICIT NONE
   TYPE(constr_type_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_constr_type


SUBROUTINE qes_write_constr_parms_list(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(constr_parms_list_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_constr_parms_list

SUBROUTINE qes_init_constr_parms_list(obj, tagname)
   IMPLICIT NONE

   TYPE(constr_parms_list_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

END SUBROUTINE qes_init_constr_parms_list

SUBROUTINE qes_reset_constr_parms_list(obj)
   IMPLICIT NONE
   TYPE(constr_parms_list_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_constr_parms_list


SUBROUTINE qes_write_atomic_constraint(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(atomic_constraint_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'constr_parms')
         WRITE(iun, '(4(E24.16))') obj%constr_parms
      CALL iotk_write_end(iun, 'constr_parms')
      CALL iotk_write_begin(iun, 'constr_type',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%constr_type)
      CALL iotk_write_end(iun, 'constr_type',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'constr_target')
         WRITE(iun, '(E24.16)') obj%constr_target
      CALL iotk_write_end(iun, 'constr_target')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_atomic_constraint

SUBROUTINE qes_init_atomic_constraint(obj, tagname, constr_parms, constr_type, &
                              constr_target)
   IMPLICIT NONE

   TYPE(atomic_constraint_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   REAL(DP), DIMENSION(4) :: constr_parms
   CHARACTER(len=*) :: constr_type
   REAL(DP) :: constr_target

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%constr_parms = constr_parms
   obj%constr_type = constr_type
   obj%constr_target = constr_target

END SUBROUTINE qes_init_atomic_constraint

SUBROUTINE qes_reset_atomic_constraint(obj)
   IMPLICIT NONE
   TYPE(atomic_constraint_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_atomic_constraint


SUBROUTINE qes_write_atomic_constraints(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(atomic_constraints_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'num_of_constraints')
         WRITE(iun, '(I12)') obj%num_of_constraints
      CALL iotk_write_end(iun, 'num_of_constraints')
      CALL iotk_write_begin(iun, 'tolerance')
         WRITE(iun, '(E24.16)') obj%tolerance
      CALL iotk_write_end(iun, 'tolerance')
      DO i = 1, obj%ndim_atomic_constraint
         CALL qes_write_atomic_constraint(iun, obj%atomic_constraint(i))
         !
      END DO
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_atomic_constraints

SUBROUTINE qes_init_atomic_constraints(obj, tagname, num_of_constraints, tolerance, &
                              ndim_atomic_constraint, atomic_constraint)
   IMPLICIT NONE

   TYPE(atomic_constraints_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: num_of_constraints
   REAL(DP) :: tolerance
   INTEGER  :: ndim_atomic_constraint
   TYPE(atomic_constraint_type ), DIMENSION( ndim_atomic_constraint )  :: atomic_constraint

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%num_of_constraints = num_of_constraints
   obj%tolerance = tolerance
   ALLOCATE(obj%atomic_constraint(SIZE(atomic_constraint)))
   DO i = 1, SIZE(atomic_constraint)
      obj%atomic_constraint(i) = atomic_constraint(i)
   ENDDO
   obj%ndim_atomic_constraint = ndim_atomic_constraint

END SUBROUTINE qes_init_atomic_constraints

SUBROUTINE qes_reset_atomic_constraints(obj)
   IMPLICIT NONE
   TYPE(atomic_constraints_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   DO i = 1, SIZE(obj%atomic_constraint)
      CALL qes_reset_atomic_constraint(obj%atomic_constraint(i))
   ENDDO
   IF (ALLOCATED(obj%atomic_constraint)) DEALLOCATE(obj%atomic_constraint)

END SUBROUTINE qes_reset_atomic_constraints


SUBROUTINE qes_write_electric_potential(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(electric_potential_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_electric_potential

SUBROUTINE qes_init_electric_potential(obj, tagname)
   IMPLICIT NONE

   TYPE(electric_potential_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

END SUBROUTINE qes_init_electric_potential

SUBROUTINE qes_reset_electric_potential(obj)
   IMPLICIT NONE
   TYPE(electric_potential_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_electric_potential


SUBROUTINE qes_write_electric_field(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(electric_field_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'electric_potential',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%electric_potential)
      CALL iotk_write_end(iun, 'electric_potential',indentation=.FALSE.)
      IF(obj%dipole_correction_ispresent) THEN
         CALL iotk_write_begin(iun, 'dipole_correction',new_line=.FALSE.)
            IF (obj%dipole_correction) THEN
               WRITE(iun, '(A)',advance='no')  'true'
            ELSE
               WRITE(iun, '(A)',advance='no')  'false'
            ENDIF
         CALL iotk_write_end(iun, 'dipole_correction',indentation=.FALSE.)
      ENDIF
      !
      IF(obj%electric_field_direction_ispresent) THEN
         CALL iotk_write_begin(iun, 'electric_field_direction')
            WRITE(iun, '(I12)') obj%electric_field_direction
         CALL iotk_write_end(iun, 'electric_field_direction')
      ENDIF
      !
      IF(obj%potential_max_position_ispresent) THEN
         CALL iotk_write_begin(iun, 'potential_max_position')
            WRITE(iun, '(E24.16)') obj%potential_max_position
         CALL iotk_write_end(iun, 'potential_max_position')
      ENDIF
      !
      IF(obj%potential_decrease_width_ispresent) THEN
         CALL iotk_write_begin(iun, 'potential_decrease_width')
            WRITE(iun, '(E24.16)') obj%potential_decrease_width
         CALL iotk_write_end(iun, 'potential_decrease_width')
      ENDIF
      !
      IF(obj%electric_field_amplitude_ispresent) THEN
         CALL iotk_write_begin(iun, 'electric_field_amplitude')
            WRITE(iun, '(E24.16)') obj%electric_field_amplitude
         CALL iotk_write_end(iun, 'electric_field_amplitude')
      ENDIF
      !
      IF(obj%electric_field_vector_ispresent) THEN
         CALL iotk_write_begin(iun, 'electric_field_vector')
            WRITE(iun, '(3(E24.16))') obj%electric_field_vector
         CALL iotk_write_end(iun, 'electric_field_vector')
      ENDIF
      !
      IF(obj%nk_per_string_ispresent) THEN
         CALL iotk_write_begin(iun, 'nk_per_string')
            WRITE(iun, '(I12)') obj%nk_per_string
         CALL iotk_write_end(iun, 'nk_per_string')
      ENDIF
      !
      IF(obj%n_berry_cycles_ispresent) THEN
         CALL iotk_write_begin(iun, 'n_berry_cycles')
            WRITE(iun, '(I12)') obj%n_berry_cycles
         CALL iotk_write_end(iun, 'n_berry_cycles')
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_electric_field

SUBROUTINE qes_init_electric_field(obj, tagname, electric_potential, &
                              dipole_correction_ispresent, dipole_correction, &
                              electric_field_direction_ispresent, electric_field_direction, &
                              potential_max_position_ispresent, potential_max_position, &
                              potential_decrease_width_ispresent, potential_decrease_width, &
                              electric_field_amplitude_ispresent, electric_field_amplitude, &
                              electric_field_vector_ispresent, electric_field_vector, &
                              nk_per_string_ispresent, nk_per_string, &
                              n_berry_cycles_ispresent, n_berry_cycles)
   IMPLICIT NONE

   TYPE(electric_field_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: electric_potential
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

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%electric_potential = electric_potential
   obj%dipole_correction_ispresent = dipole_correction_ispresent
   IF(obj%dipole_correction_ispresent) THEN
      obj%dipole_correction = dipole_correction
   ENDIF
   obj%electric_field_direction_ispresent = electric_field_direction_ispresent
   IF(obj%electric_field_direction_ispresent) THEN
      obj%electric_field_direction = electric_field_direction
   ENDIF
   obj%potential_max_position_ispresent = potential_max_position_ispresent
   IF(obj%potential_max_position_ispresent) THEN
      obj%potential_max_position = potential_max_position
   ENDIF
   obj%potential_decrease_width_ispresent = potential_decrease_width_ispresent
   IF(obj%potential_decrease_width_ispresent) THEN
      obj%potential_decrease_width = potential_decrease_width
   ENDIF
   obj%electric_field_amplitude_ispresent = electric_field_amplitude_ispresent
   IF(obj%electric_field_amplitude_ispresent) THEN
      obj%electric_field_amplitude = electric_field_amplitude
   ENDIF
   obj%electric_field_vector_ispresent = electric_field_vector_ispresent
   IF(obj%electric_field_vector_ispresent) THEN
      obj%electric_field_vector = electric_field_vector
   ENDIF
   obj%nk_per_string_ispresent = nk_per_string_ispresent
   IF(obj%nk_per_string_ispresent) THEN
      obj%nk_per_string = nk_per_string
   ENDIF
   obj%n_berry_cycles_ispresent = n_berry_cycles_ispresent
   IF(obj%n_berry_cycles_ispresent) THEN
      obj%n_berry_cycles = n_berry_cycles
   ENDIF

END SUBROUTINE qes_init_electric_field

SUBROUTINE qes_reset_electric_field(obj)
   IMPLICIT NONE
   TYPE(electric_field_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%dipole_correction_ispresent) THEN
      obj%dipole_correction_ispresent = .FALSE.
   ENDIF
   IF(obj%electric_field_direction_ispresent) THEN
      obj%electric_field_direction_ispresent = .FALSE.
   ENDIF
   IF(obj%potential_max_position_ispresent) THEN
      obj%potential_max_position_ispresent = .FALSE.
   ENDIF
   IF(obj%potential_decrease_width_ispresent) THEN
      obj%potential_decrease_width_ispresent = .FALSE.
   ENDIF
   IF(obj%electric_field_amplitude_ispresent) THEN
      obj%electric_field_amplitude_ispresent = .FALSE.
   ENDIF
   IF(obj%electric_field_vector_ispresent) THEN
      obj%electric_field_vector_ispresent = .FALSE.
   ENDIF
   IF(obj%nk_per_string_ispresent) THEN
      obj%nk_per_string_ispresent = .FALSE.
   ENDIF
   IF(obj%n_berry_cycles_ispresent) THEN
      obj%n_berry_cycles_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_electric_field


SUBROUTINE qes_write_symmetries(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(symmetries_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'nsym',new_line=.FALSE.)
         WRITE(iun, '(I0)',advance='no') obj%nsym
      CALL iotk_write_end(iun, 'nsym',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'nrot',new_line=.FALSE.)
         WRITE(iun, '(I0)',advance='no') obj%nrot
      CALL iotk_write_end(iun, 'nrot',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'space_group',new_line=.FALSE.)
         WRITE(iun, '(I0)',advance='no') obj%space_group
      CALL iotk_write_end(iun, 'space_group',indentation=.FALSE.)
      DO i = 1, obj%ndim_symmetry
         CALL qes_write_symmetry(iun, obj%symmetry(i))
         !
      END DO
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_symmetries

SUBROUTINE qes_init_symmetries(obj, tagname, nsym, nrot, space_group, ndim_symmetry, &
                              symmetry)
   IMPLICIT NONE

   TYPE(symmetries_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER :: nsym
   INTEGER :: nrot
   INTEGER :: space_group
   INTEGER  :: ndim_symmetry
   TYPE(symmetry_type ), DIMENSION( ndim_symmetry )  :: symmetry

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%nsym = nsym
   obj%nrot = nrot
   obj%space_group = space_group
   ALLOCATE(obj%symmetry(SIZE(symmetry)))
   DO i = 1, SIZE(symmetry)
      obj%symmetry(i) = symmetry(i)
   ENDDO
   obj%ndim_symmetry = ndim_symmetry

END SUBROUTINE qes_init_symmetries

SUBROUTINE qes_reset_symmetries(obj)
   IMPLICIT NONE
   TYPE(symmetries_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   DO i = 1, SIZE(obj%symmetry)
      CALL qes_reset_symmetry(obj%symmetry(i))
   ENDDO
   IF (ALLOCATED(obj%symmetry)) DEALLOCATE(obj%symmetry)

END SUBROUTINE qes_reset_symmetries


SUBROUTINE qes_write_ekin_functional(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(ekin_functional_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'ecfixed')
         WRITE(iun, '(E24.16)') obj%ecfixed
      CALL iotk_write_end(iun, 'ecfixed')
      CALL iotk_write_begin(iun, 'qcutz')
         WRITE(iun, '(E24.16)') obj%qcutz
      CALL iotk_write_end(iun, 'qcutz')
      CALL iotk_write_begin(iun, 'q2sigma')
         WRITE(iun, '(E24.16)') obj%q2sigma
      CALL iotk_write_end(iun, 'q2sigma')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_ekin_functional

SUBROUTINE qes_init_ekin_functional(obj, tagname, ecfixed, qcutz, q2sigma)
   IMPLICIT NONE

   TYPE(ekin_functional_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   REAL(DP) :: ecfixed
   REAL(DP) :: qcutz
   REAL(DP) :: q2sigma

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%ecfixed = ecfixed
   obj%qcutz = qcutz
   obj%q2sigma = q2sigma

END SUBROUTINE qes_init_ekin_functional

SUBROUTINE qes_reset_ekin_functional(obj)
   IMPLICIT NONE
   TYPE(ekin_functional_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_ekin_functional


SUBROUTINE qes_write_esm(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(esm_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'bc',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%bc)
      CALL iotk_write_end(iun, 'bc',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'nfit')
         WRITE(iun, '(I12)') obj%nfit
      CALL iotk_write_end(iun, 'nfit')
      CALL iotk_write_begin(iun, 'w')
         WRITE(iun, '(E24.16)') obj%w
      CALL iotk_write_end(iun, 'w')
      CALL iotk_write_begin(iun, 'efield')
         WRITE(iun, '(E24.16)') obj%efield
      CALL iotk_write_end(iun, 'efield')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_esm

SUBROUTINE qes_init_esm(obj, tagname, bc, nfit, w, efield)
   IMPLICIT NONE

   TYPE(esm_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: bc
   INTEGER  :: nfit
   REAL(DP) :: w
   REAL(DP) :: efield

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%bc = bc
   obj%nfit = nfit
   obj%w = w
   obj%efield = efield

END SUBROUTINE qes_init_esm

SUBROUTINE qes_reset_esm(obj)
   IMPLICIT NONE
   TYPE(esm_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_esm


SUBROUTINE qes_write_boundary_conditions(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(boundary_conditions_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'assume_isolated',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%assume_isolated)
      CALL iotk_write_end(iun, 'assume_isolated',indentation=.FALSE.)
      IF(obj%esm_ispresent) THEN
         CALL qes_write_esm(iun, obj%esm)
         !
      ENDIF
      !
      IF(obj%fcp_opt_ispresent) THEN
         CALL iotk_write_begin(iun, 'fcp_opt',new_line=.FALSE.)
            IF (obj%fcp_opt) THEN
               WRITE(iun, '(A)',advance='no')  'true'
            ELSE
               WRITE(iun, '(A)',advance='no')  'false'
            ENDIF
         CALL iotk_write_end(iun, 'fcp_opt',indentation=.FALSE.)
      ENDIF
      !
      IF(obj%fcp_mu_ispresent) THEN
         CALL iotk_write_begin(iun, 'fcp_mu')
            WRITE(iun, '(E24.16)') obj%fcp_mu
         CALL iotk_write_end(iun, 'fcp_mu')
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_boundary_conditions

SUBROUTINE qes_init_boundary_conditions(obj, tagname, assume_isolated, esm_ispresent, esm, &
                              fcp_opt_ispresent, fcp_opt, fcp_mu_ispresent, fcp_mu)
   IMPLICIT NONE

   TYPE(boundary_conditions_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: assume_isolated
   LOGICAL  :: esm_ispresent
   TYPE(esm_type) :: esm
   LOGICAL  :: fcp_opt_ispresent
   LOGICAL  :: fcp_opt
   LOGICAL  :: fcp_mu_ispresent
   REAL(DP) :: fcp_mu

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%assume_isolated = assume_isolated
   obj%esm_ispresent = esm_ispresent
   IF(obj%esm_ispresent) THEN
      obj%esm = esm
   ENDIF
   obj%fcp_opt_ispresent = fcp_opt_ispresent
   IF(obj%fcp_opt_ispresent) THEN
      obj%fcp_opt = fcp_opt
   ENDIF
   obj%fcp_mu_ispresent = fcp_mu_ispresent
   IF(obj%fcp_mu_ispresent) THEN
      obj%fcp_mu = fcp_mu
   ENDIF

END SUBROUTINE qes_init_boundary_conditions

SUBROUTINE qes_reset_boundary_conditions(obj)
   IMPLICIT NONE
   TYPE(boundary_conditions_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%esm_ispresent) THEN
      CALL qes_reset_esm(obj%esm)
      obj%esm_ispresent = .FALSE.
   ENDIF
   IF(obj%fcp_opt_ispresent) THEN
      obj%fcp_opt_ispresent = .FALSE.
   ENDIF
   IF(obj%fcp_mu_ispresent) THEN
      obj%fcp_mu_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_boundary_conditions


SUBROUTINE qes_write_symmetry_flags(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(symmetry_flags_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'nosym',new_line=.FALSE.)
         IF (obj%nosym) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'nosym',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'nosym_evc',new_line=.FALSE.)
         IF (obj%nosym_evc) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'nosym_evc',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'noinv',new_line=.FALSE.)
         IF (obj%noinv) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'noinv',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'no_t_rev',new_line=.FALSE.)
         IF (obj%no_t_rev) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'no_t_rev',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'force_symmorphic',new_line=.FALSE.)
         IF (obj%force_symmorphic) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'force_symmorphic',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'use_all_frac',new_line=.FALSE.)
         IF (obj%use_all_frac) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'use_all_frac',indentation=.FALSE.)
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_symmetry_flags

SUBROUTINE qes_init_symmetry_flags(obj, tagname, nosym, nosym_evc, noinv, no_t_rev, &
                              force_symmorphic, use_all_frac)
   IMPLICIT NONE

   TYPE(symmetry_flags_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   LOGICAL  :: nosym
   LOGICAL  :: nosym_evc
   LOGICAL  :: noinv
   LOGICAL  :: no_t_rev
   LOGICAL  :: force_symmorphic
   LOGICAL  :: use_all_frac

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%nosym = nosym
   obj%nosym_evc = nosym_evc
   obj%noinv = noinv
   obj%no_t_rev = no_t_rev
   obj%force_symmorphic = force_symmorphic
   obj%use_all_frac = use_all_frac

END SUBROUTINE qes_init_symmetry_flags

SUBROUTINE qes_reset_symmetry_flags(obj)
   IMPLICIT NONE
   TYPE(symmetry_flags_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_symmetry_flags


SUBROUTINE qes_write_integerMatrix(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(integerMatrix_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      WRITE(fmtstr,'(a)') '(12I5)'
      DO i = 1, SIZE(obj%int_mat,2)
         WRITE(iun, fmtstr) obj%int_mat(:,i)
      ENDDO
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_integerMatrix

SUBROUTINE qes_init_integerMatrix(obj, tagname, ndim1_int_mat, ndim2_int_mat, int_mat)
   IMPLICIT NONE

   TYPE(integerMatrix_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: ndim1_int_mat
   INTEGER  :: ndim2_int_mat
   INTEGER, DIMENSION(ndim1_int_mat,ndim2_int_mat) :: int_mat

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   ALLOCATE(obj%int_mat(ndim1_int_mat,ndim2_int_mat))
   obj%int_mat(:,:) = int_mat(:,:)
   obj%ndim1_int_mat = ndim1_int_mat
   obj%ndim2_int_mat = ndim2_int_mat

END SUBROUTINE qes_init_integerMatrix

SUBROUTINE qes_reset_integerMatrix(obj)
   IMPLICIT NONE
   TYPE(integerMatrix_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF (ALLOCATED(obj%int_mat))  DEALLOCATE(obj%int_mat)

END SUBROUTINE qes_reset_integerMatrix


SUBROUTINE qes_write_cell_control(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(cell_control_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'cell_dynamics',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%cell_dynamics)
      CALL iotk_write_end(iun, 'cell_dynamics',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'pressure')
         WRITE(iun, '(E24.16)') obj%pressure
      CALL iotk_write_end(iun, 'pressure')
      IF(obj%wmass_ispresent) THEN
         CALL iotk_write_begin(iun, 'wmass')
            WRITE(iun, '(E24.16)') obj%wmass
         CALL iotk_write_end(iun, 'wmass')
      ENDIF
      !
      IF(obj%cell_factor_ispresent) THEN
         CALL iotk_write_begin(iun, 'cell_factor')
            WRITE(iun, '(E24.16)') obj%cell_factor
         CALL iotk_write_end(iun, 'cell_factor')
      ENDIF
      !
      IF(obj%fix_volume_ispresent) THEN
         CALL iotk_write_begin(iun, 'fix_volume',new_line=.FALSE.)
            IF (obj%fix_volume) THEN
               WRITE(iun, '(A)',advance='no')  'true'
            ELSE
               WRITE(iun, '(A)',advance='no')  'false'
            ENDIF
         CALL iotk_write_end(iun, 'fix_volume',indentation=.FALSE.)
      ENDIF
      !
      IF(obj%fix_area_ispresent) THEN
         CALL iotk_write_begin(iun, 'fix_area',new_line=.FALSE.)
            IF (obj%fix_area) THEN
               WRITE(iun, '(A)',advance='no')  'true'
            ELSE
               WRITE(iun, '(A)',advance='no')  'false'
            ENDIF
         CALL iotk_write_end(iun, 'fix_area',indentation=.FALSE.)
      ENDIF
      !
      IF(obj%isotropic_ispresent) THEN
         CALL iotk_write_begin(iun, 'isotropic',new_line=.FALSE.)
            IF (obj%isotropic) THEN
               WRITE(iun, '(A)',advance='no')  'true'
            ELSE
               WRITE(iun, '(A)',advance='no')  'false'
            ENDIF
         CALL iotk_write_end(iun, 'isotropic',indentation=.FALSE.)
      ENDIF
      !
      IF(obj%free_cell_ispresent) THEN
         CALL qes_write_integerMatrix(iun, obj%free_cell)
         !
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_cell_control

SUBROUTINE qes_init_cell_control(obj, tagname, cell_dynamics, pressure, wmass_ispresent, &
                              wmass, cell_factor_ispresent, cell_factor, &
                              fix_volume_ispresent, fix_volume, fix_area_ispresent, &
                              fix_area, isotropic_ispresent, isotropic, &
                              free_cell_ispresent, free_cell)
   IMPLICIT NONE

   TYPE(cell_control_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: cell_dynamics
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

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%cell_dynamics = cell_dynamics
   obj%pressure = pressure
   obj%wmass_ispresent = wmass_ispresent
   IF(obj%wmass_ispresent) THEN
      obj%wmass = wmass
   ENDIF
   obj%cell_factor_ispresent = cell_factor_ispresent
   IF(obj%cell_factor_ispresent) THEN
      obj%cell_factor = cell_factor
   ENDIF
   obj%fix_volume_ispresent = fix_volume_ispresent
   IF(obj%fix_volume_ispresent) THEN
      obj%fix_volume = fix_volume
   ENDIF
   obj%fix_area_ispresent = fix_area_ispresent
   IF(obj%fix_area_ispresent) THEN
      obj%fix_area = fix_area
   ENDIF
   obj%isotropic_ispresent = isotropic_ispresent
   IF(obj%isotropic_ispresent) THEN
      obj%isotropic = isotropic
   ENDIF
   obj%free_cell_ispresent = free_cell_ispresent
   IF(obj%free_cell_ispresent) THEN
      obj%free_cell = free_cell
   ENDIF

END SUBROUTINE qes_init_cell_control

SUBROUTINE qes_reset_cell_control(obj)
   IMPLICIT NONE
   TYPE(cell_control_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%wmass_ispresent) THEN
      obj%wmass_ispresent = .FALSE.
   ENDIF
   IF(obj%cell_factor_ispresent) THEN
      obj%cell_factor_ispresent = .FALSE.
   ENDIF
   IF(obj%fix_volume_ispresent) THEN
      obj%fix_volume_ispresent = .FALSE.
   ENDIF
   IF(obj%fix_area_ispresent) THEN
      obj%fix_area_ispresent = .FALSE.
   ENDIF
   IF(obj%isotropic_ispresent) THEN
      obj%isotropic_ispresent = .FALSE.
   ENDIF
   IF(obj%free_cell_ispresent) THEN
      CALL qes_reset_integerMatrix(obj%free_cell)
      obj%free_cell_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_cell_control


SUBROUTINE qes_write_md(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(md_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'pot_extrapolation',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%pot_extrapolation)
      CALL iotk_write_end(iun, 'pot_extrapolation',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'wfc_extrapolation',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%wfc_extrapolation)
      CALL iotk_write_end(iun, 'wfc_extrapolation',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'ion_temperature',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%ion_temperature)
      CALL iotk_write_end(iun, 'ion_temperature',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'timestep')
         WRITE(iun, '(E24.16)') obj%timestep
      CALL iotk_write_end(iun, 'timestep')
      CALL iotk_write_begin(iun, 'tempw')
         WRITE(iun, '(E24.16)') obj%tempw
      CALL iotk_write_end(iun, 'tempw')
      CALL iotk_write_begin(iun, 'tolp')
         WRITE(iun, '(E24.16)') obj%tolp
      CALL iotk_write_end(iun, 'tolp')
      CALL iotk_write_begin(iun, 'deltaT')
         WRITE(iun, '(E24.16)') obj%deltaT
      CALL iotk_write_end(iun, 'deltaT')
      CALL iotk_write_begin(iun, 'nraise')
         WRITE(iun, '(I12)') obj%nraise
      CALL iotk_write_end(iun, 'nraise')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_md

SUBROUTINE qes_init_md(obj, tagname, pot_extrapolation, wfc_extrapolation, ion_temperature, &
                              timestep, tempw, tolp, deltaT, nraise)
   IMPLICIT NONE

   TYPE(md_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: pot_extrapolation
   CHARACTER(len=*) :: wfc_extrapolation
   CHARACTER(len=*) :: ion_temperature
   REAL(DP) :: timestep
   REAL(DP) :: tempw
   REAL(DP) :: tolp
   REAL(DP) :: deltaT
   INTEGER  :: nraise

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%pot_extrapolation = pot_extrapolation
   obj%wfc_extrapolation = wfc_extrapolation
   obj%ion_temperature = ion_temperature
   obj%timestep = timestep
   obj%tempw = tempw
   obj%tolp = tolp
   obj%deltaT = deltaT
   obj%nraise = nraise

END SUBROUTINE qes_init_md

SUBROUTINE qes_reset_md(obj)
   IMPLICIT NONE
   TYPE(md_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_md


SUBROUTINE qes_write_bfgs(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(bfgs_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'ndim')
         WRITE(iun, '(I12)') obj%ndim
      CALL iotk_write_end(iun, 'ndim')
      CALL iotk_write_begin(iun, 'trust_radius_min')
         WRITE(iun, '(E24.16)') obj%trust_radius_min
      CALL iotk_write_end(iun, 'trust_radius_min')
      CALL iotk_write_begin(iun, 'trust_radius_max')
         WRITE(iun, '(E24.16)') obj%trust_radius_max
      CALL iotk_write_end(iun, 'trust_radius_max')
      CALL iotk_write_begin(iun, 'trust_radius_init')
         WRITE(iun, '(E24.16)') obj%trust_radius_init
      CALL iotk_write_end(iun, 'trust_radius_init')
      CALL iotk_write_begin(iun, 'w1')
         WRITE(iun, '(E24.16)') obj%w1
      CALL iotk_write_end(iun, 'w1')
      CALL iotk_write_begin(iun, 'w2')
         WRITE(iun, '(E24.16)') obj%w2
      CALL iotk_write_end(iun, 'w2')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_bfgs

SUBROUTINE qes_init_bfgs(obj, tagname, ndim, trust_radius_min, trust_radius_max, &
                              trust_radius_init, w1, w2)
   IMPLICIT NONE

   TYPE(bfgs_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: ndim
   REAL(DP) :: trust_radius_min
   REAL(DP) :: trust_radius_max
   REAL(DP) :: trust_radius_init
   REAL(DP) :: w1
   REAL(DP) :: w2

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%ndim = ndim
   obj%trust_radius_min = trust_radius_min
   obj%trust_radius_max = trust_radius_max
   obj%trust_radius_init = trust_radius_init
   obj%w1 = w1
   obj%w2 = w2

END SUBROUTINE qes_init_bfgs

SUBROUTINE qes_reset_bfgs(obj)
   IMPLICIT NONE
   TYPE(bfgs_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_bfgs


SUBROUTINE qes_write_ion_control(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(ion_control_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'ion_dynamics',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%ion_dynamics)
      CALL iotk_write_end(iun, 'ion_dynamics',indentation=.FALSE.)
      IF(obj%upscale_ispresent) THEN
         CALL iotk_write_begin(iun, 'upscale')
            WRITE(iun, '(E24.16)') obj%upscale
         CALL iotk_write_end(iun, 'upscale')
      ENDIF
      !
      IF(obj%remove_rigid_rot_ispresent) THEN
         CALL iotk_write_begin(iun, 'remove_rigid_rot',new_line=.FALSE.)
            IF (obj%remove_rigid_rot) THEN
               WRITE(iun, '(A)',advance='no')  'true'
            ELSE
               WRITE(iun, '(A)',advance='no')  'false'
            ENDIF
         CALL iotk_write_end(iun, 'remove_rigid_rot',indentation=.FALSE.)
      ENDIF
      !
      IF(obj%refold_pos_ispresent) THEN
         CALL iotk_write_begin(iun, 'refold_pos',new_line=.FALSE.)
            IF (obj%refold_pos) THEN
               WRITE(iun, '(A)',advance='no')  'true'
            ELSE
               WRITE(iun, '(A)',advance='no')  'false'
            ENDIF
         CALL iotk_write_end(iun, 'refold_pos',indentation=.FALSE.)
      ENDIF
      !
      IF(obj%bfgs_ispresent) THEN
         CALL qes_write_bfgs(iun, obj%bfgs)
         !
      ENDIF
      !
      IF(obj%md_ispresent) THEN
         CALL qes_write_md(iun, obj%md)
         !
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_ion_control

SUBROUTINE qes_init_ion_control(obj, tagname, ion_dynamics, upscale_ispresent, upscale, &
                              remove_rigid_rot_ispresent, remove_rigid_rot, &
                              refold_pos_ispresent, refold_pos, bfgs_ispresent, bfgs, &
                              md_ispresent, md)
   IMPLICIT NONE

   TYPE(ion_control_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: ion_dynamics
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

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%ion_dynamics = ion_dynamics
   obj%upscale_ispresent = upscale_ispresent
   IF(obj%upscale_ispresent) THEN
      obj%upscale = upscale
   ENDIF
   obj%remove_rigid_rot_ispresent = remove_rigid_rot_ispresent
   IF(obj%remove_rigid_rot_ispresent) THEN
      obj%remove_rigid_rot = remove_rigid_rot
   ENDIF
   obj%refold_pos_ispresent = refold_pos_ispresent
   IF(obj%refold_pos_ispresent) THEN
      obj%refold_pos = refold_pos
   ENDIF
   obj%bfgs_ispresent = bfgs_ispresent
   IF(obj%bfgs_ispresent) THEN
      obj%bfgs = bfgs
   ENDIF
   obj%md_ispresent = md_ispresent
   IF(obj%md_ispresent) THEN
      obj%md = md
   ENDIF

END SUBROUTINE qes_init_ion_control

SUBROUTINE qes_reset_ion_control(obj)
   IMPLICIT NONE
   TYPE(ion_control_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%upscale_ispresent) THEN
      obj%upscale_ispresent = .FALSE.
   ENDIF
   IF(obj%remove_rigid_rot_ispresent) THEN
      obj%remove_rigid_rot_ispresent = .FALSE.
   ENDIF
   IF(obj%refold_pos_ispresent) THEN
      obj%refold_pos_ispresent = .FALSE.
   ENDIF
   IF(obj%bfgs_ispresent) THEN
      CALL qes_reset_bfgs(obj%bfgs)
      obj%bfgs_ispresent = .FALSE.
   ENDIF
   IF(obj%md_ispresent) THEN
      CALL qes_reset_md(obj%md)
      obj%md_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_ion_control


SUBROUTINE qes_write_monkhorst_pack(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(monkhorst_pack_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'nk1', obj%nk1)
   CALL iotk_write_attr(attr, 'nk2', obj%nk2)
   CALL iotk_write_attr(attr, 'nk3', obj%nk3)
   CALL iotk_write_attr(attr, 'k1', obj%k1)
   CALL iotk_write_attr(attr, 'k2', obj%k2)
   CALL iotk_write_attr(attr, 'k3', obj%k3)

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr),new_line=.FALSE.)
      !
      WRITE(iun, '(A)',advance='no')  TRIM(obj%monkhorst_pack)
   CALL iotk_write_end(iun, TRIM(obj%tagname),indentation=.FALSE.)
   !
END SUBROUTINE qes_write_monkhorst_pack

SUBROUTINE qes_init_monkhorst_pack(obj, tagname, nk1, nk2, nk3, k1, k2, k3, monkhorst_pack)
   IMPLICIT NONE

   TYPE(monkhorst_pack_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: nk1
   INTEGER  :: nk2
   INTEGER  :: nk3
   INTEGER  :: k1
   INTEGER  :: k2
   INTEGER  :: k3
   CHARACTER(len=*) :: monkhorst_pack

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%nk1 = nk1


   obj%nk2 = nk2


   obj%nk3 = nk3


   obj%k1 = k1


   obj%k2 = k2


   obj%k3 = k3

   obj%monkhorst_pack = monkhorst_pack

END SUBROUTINE qes_init_monkhorst_pack

SUBROUTINE qes_reset_monkhorst_pack(obj)
   IMPLICIT NONE
   TYPE(monkhorst_pack_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_monkhorst_pack


SUBROUTINE qes_write_k_points_IBZ(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(k_points_IBZ_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      IF(obj%monkhorst_pack_ispresent) THEN
         CALL qes_write_monkhorst_pack(iun, obj%monkhorst_pack)
         !
      ENDIF
      !
      IF(obj%nk_ispresent) THEN
         CALL iotk_write_begin(iun, 'nk')
            WRITE(iun, '(I12)') obj%nk
         CALL iotk_write_end(iun, 'nk')
      ENDIF
      !
      IF(obj%k_point_ispresent) THEN
         DO i = 1, obj%ndim_k_point
            CALL qes_write_k_point(iun, obj%k_point(i))
            !
         END DO
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_k_points_IBZ

SUBROUTINE qes_init_k_points_IBZ(obj, tagname, monkhorst_pack_ispresent, monkhorst_pack, &
                              nk_ispresent, nk, k_point_ispresent, ndim_k_point, k_point)
   IMPLICIT NONE

   TYPE(k_points_IBZ_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   LOGICAL  :: monkhorst_pack_ispresent
   TYPE(monkhorst_pack_type) :: monkhorst_pack
   LOGICAL  :: nk_ispresent
   INTEGER  :: nk
   LOGICAL  :: k_point_ispresent
   INTEGER  :: ndim_k_point
   TYPE(k_point_type ), DIMENSION( ndim_k_point )  :: k_point

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%monkhorst_pack_ispresent = monkhorst_pack_ispresent
   IF(obj%monkhorst_pack_ispresent) THEN
      obj%monkhorst_pack = monkhorst_pack
   ENDIF
   obj%nk_ispresent = nk_ispresent
   IF(obj%nk_ispresent) THEN
      obj%nk = nk
   ENDIF
   obj%k_point_ispresent = k_point_ispresent
   IF(obj%k_point_ispresent) THEN
      ALLOCATE(obj%k_point(SIZE(k_point)))
      DO i = 1, SIZE(k_point)
         obj%k_point(i) = k_point(i)
      ENDDO
   obj%ndim_k_point = ndim_k_point
   ENDIF

END SUBROUTINE qes_init_k_points_IBZ

SUBROUTINE qes_reset_k_points_IBZ(obj)
   IMPLICIT NONE
   TYPE(k_points_IBZ_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%monkhorst_pack_ispresent) THEN
      CALL qes_reset_monkhorst_pack(obj%monkhorst_pack)
      obj%monkhorst_pack_ispresent = .FALSE.
   ENDIF
   IF(obj%nk_ispresent) THEN
      obj%nk_ispresent = .FALSE.
   ENDIF
   IF(obj%k_point_ispresent) THEN
      DO i = 1, SIZE(obj%k_point)
         CALL qes_reset_k_point(obj%k_point(i))
      ENDDO
      IF (ALLOCATED(obj%k_point)) DEALLOCATE(obj%k_point)
      obj%k_point_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_k_points_IBZ


SUBROUTINE qes_write_occupations(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(occupations_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   IF(obj%spin_ispresent) THEN
      CALL iotk_write_attr(attr, 'spin', obj%spin)
   END IF

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr),new_line=.FALSE.)
      !
      WRITE(iun, '(A)',advance='no')  TRIM(obj%occupations)
   CALL iotk_write_end(iun, TRIM(obj%tagname),indentation=.FALSE.)
   !
END SUBROUTINE qes_write_occupations

SUBROUTINE qes_init_occupations(obj, tagname, spin, spin_ispresent, occupations)
   IMPLICIT NONE

   TYPE(occupations_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   LOGICAL  :: spin_ispresent
   INTEGER , OPTIONAL :: spin
   CHARACTER(len=*) :: occupations

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%spin_ispresent = spin_ispresent
   IF (obj%spin_ispresent) THEN
      obj%spin = spin
   ENDIF

   obj%occupations = occupations

END SUBROUTINE qes_init_occupations

SUBROUTINE qes_reset_occupations(obj)
   IMPLICIT NONE
   TYPE(occupations_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_occupations


SUBROUTINE qes_write_mixingMode(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(mixingMode_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_mixingMode

SUBROUTINE qes_init_mixingMode(obj, tagname)
   IMPLICIT NONE

   TYPE(mixingMode_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

END SUBROUTINE qes_init_mixingMode

SUBROUTINE qes_reset_mixingMode(obj)
   IMPLICIT NONE
   TYPE(mixingMode_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_mixingMode


SUBROUTINE qes_write_diago(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(diago_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_diago

SUBROUTINE qes_init_diago(obj, tagname)
   IMPLICIT NONE

   TYPE(diago_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

END SUBROUTINE qes_init_diago

SUBROUTINE qes_reset_diago(obj)
   IMPLICIT NONE
   TYPE(diago_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_diago


SUBROUTINE qes_write_electron_control(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(electron_control_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'diagonalization',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%diagonalization)
      CALL iotk_write_end(iun, 'diagonalization',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'mixing_mode',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%mixing_mode)
      CALL iotk_write_end(iun, 'mixing_mode',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'mixing_beta')
         WRITE(iun, '(E24.16)') obj%mixing_beta
      CALL iotk_write_end(iun, 'mixing_beta')
      CALL iotk_write_begin(iun, 'conv_thr')
         WRITE(iun, '(E24.16)') obj%conv_thr
      CALL iotk_write_end(iun, 'conv_thr')
      CALL iotk_write_begin(iun, 'mixing_ndim')
         WRITE(iun, '(I12)') obj%mixing_ndim
      CALL iotk_write_end(iun, 'mixing_ndim')
      CALL iotk_write_begin(iun, 'max_nstep')
         WRITE(iun, '(I12)') obj%max_nstep
      CALL iotk_write_end(iun, 'max_nstep')
      CALL iotk_write_begin(iun, 'real_space_q',new_line=.FALSE.)
         IF (obj%real_space_q) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'real_space_q',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'tq_smoothing',new_line=.FALSE.)
         IF (obj%tq_smoothing) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'tq_smoothing',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'tbeta_smoothing',new_line=.FALSE.)
         IF (obj%tbeta_smoothing) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'tbeta_smoothing',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'diago_thr_init')
         WRITE(iun, '(E24.16)') obj%diago_thr_init
      CALL iotk_write_end(iun, 'diago_thr_init')
      CALL iotk_write_begin(iun, 'diago_full_acc',new_line=.FALSE.)
         IF (obj%diago_full_acc) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'diago_full_acc',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'diago_cg_maxiter')
         WRITE(iun, '(I12)') obj%diago_cg_maxiter
      CALL iotk_write_end(iun, 'diago_cg_maxiter')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_electron_control

SUBROUTINE qes_init_electron_control(obj, tagname, diagonalization, mixing_mode, &
                              mixing_beta, conv_thr, mixing_ndim, max_nstep, real_space_q, &
                              tq_smoothing, tbeta_smoothing, diago_thr_init, &
                              diago_full_acc, diago_cg_maxiter)
   IMPLICIT NONE

   TYPE(electron_control_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: diagonalization
   CHARACTER(len=*) :: mixing_mode
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

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%diagonalization = diagonalization
   obj%mixing_mode = mixing_mode
   obj%mixing_beta = mixing_beta
   obj%conv_thr = conv_thr
   obj%mixing_ndim = mixing_ndim
   obj%max_nstep = max_nstep
   obj%real_space_q = real_space_q
   obj%tq_smoothing = tq_smoothing
   obj%tbeta_smoothing = tbeta_smoothing
   obj%diago_thr_init = diago_thr_init
   obj%diago_full_acc = diago_full_acc
   obj%diago_cg_maxiter = diago_cg_maxiter

END SUBROUTINE qes_init_electron_control

SUBROUTINE qes_reset_electron_control(obj)
   IMPLICIT NONE
   TYPE(electron_control_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_electron_control


SUBROUTINE qes_write_basis_set(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(basis_set_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      IF(obj%gamma_only_ispresent) THEN
         CALL iotk_write_begin(iun, 'gamma_only',new_line=.FALSE.)
            IF (obj%gamma_only) THEN
               WRITE(iun, '(A)',advance='no')  'true'
            ELSE
               WRITE(iun, '(A)',advance='no')  'false'
            ENDIF
         CALL iotk_write_end(iun, 'gamma_only',indentation=.FALSE.)
      ENDIF
      !
      CALL iotk_write_begin(iun, 'ecutwfc')
         WRITE(iun, '(E24.16)') obj%ecutwfc
      CALL iotk_write_end(iun, 'ecutwfc')
      IF(obj%ecutrho_ispresent) THEN
         CALL iotk_write_begin(iun, 'ecutrho')
            WRITE(iun, '(E24.16)') obj%ecutrho
         CALL iotk_write_end(iun, 'ecutrho')
      ENDIF
      !
      CALL qes_write_basisSetItem(iun, obj%fft_grid)
      !
      IF(obj%fft_smooth_ispresent) THEN
         CALL qes_write_basisSetItem(iun, obj%fft_smooth)
         !
      ENDIF
      !
      IF(obj%fft_box_ispresent) THEN
         CALL qes_write_basisSetItem(iun, obj%fft_box)
         !
      ENDIF
      !
      CALL iotk_write_begin(iun, 'ngm')
         WRITE(iun, '(I12)') obj%ngm
      CALL iotk_write_end(iun, 'ngm')
      IF(obj%ngms_ispresent) THEN
         CALL iotk_write_begin(iun, 'ngms')
            WRITE(iun, '(I12)') obj%ngms
         CALL iotk_write_end(iun, 'ngms')
      ENDIF
      !
      CALL iotk_write_begin(iun, 'npwx')
         WRITE(iun, '(I12)') obj%npwx
      CALL iotk_write_end(iun, 'npwx')
      CALL qes_write_reciprocal_lattice(iun, obj%reciprocal_lattice)
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_basis_set

SUBROUTINE qes_init_basis_set(obj, tagname, gamma_only_ispresent, gamma_only, ecutwfc, &
                              ecutrho_ispresent, ecutrho, fft_grid, fft_smooth_ispresent, &
                              fft_smooth, fft_box_ispresent, fft_box, ngm, ngms_ispresent, &
                              ngms, npwx, reciprocal_lattice)
   IMPLICIT NONE

   TYPE(basis_set_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
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

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%gamma_only_ispresent = gamma_only_ispresent
   IF(obj%gamma_only_ispresent) THEN
      obj%gamma_only = gamma_only
   ENDIF
   obj%ecutwfc = ecutwfc
   obj%ecutrho_ispresent = ecutrho_ispresent
   IF(obj%ecutrho_ispresent) THEN
      obj%ecutrho = ecutrho
   ENDIF
   obj%fft_grid = fft_grid
   obj%fft_smooth_ispresent = fft_smooth_ispresent
   IF(obj%fft_smooth_ispresent) THEN
      obj%fft_smooth = fft_smooth
   ENDIF
   obj%fft_box_ispresent = fft_box_ispresent
   IF(obj%fft_box_ispresent) THEN
      obj%fft_box = fft_box
   ENDIF
   obj%ngm = ngm
   obj%ngms_ispresent = ngms_ispresent
   IF(obj%ngms_ispresent) THEN
      obj%ngms = ngms
   ENDIF
   obj%npwx = npwx
   obj%reciprocal_lattice = reciprocal_lattice

END SUBROUTINE qes_init_basis_set

SUBROUTINE qes_reset_basis_set(obj)
   IMPLICIT NONE
   TYPE(basis_set_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%gamma_only_ispresent) THEN
      obj%gamma_only_ispresent = .FALSE.
   ENDIF
   IF(obj%ecutrho_ispresent) THEN
      obj%ecutrho_ispresent = .FALSE.
   ENDIF
   CALL qes_reset_basisSetItem(obj%fft_grid)
   IF(obj%fft_smooth_ispresent) THEN
      CALL qes_reset_basisSetItem(obj%fft_smooth)
      obj%fft_smooth_ispresent = .FALSE.
   ENDIF
   IF(obj%fft_box_ispresent) THEN
      CALL qes_reset_basisSetItem(obj%fft_box)
      obj%fft_box_ispresent = .FALSE.
   ENDIF
   IF(obj%ngms_ispresent) THEN
      obj%ngms_ispresent = .FALSE.
   ENDIF
   CALL qes_reset_reciprocal_lattice(obj%reciprocal_lattice)

END SUBROUTINE qes_reset_basis_set


SUBROUTINE qes_write_basis(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(basis_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      IF(obj%gamma_only_ispresent) THEN
         CALL iotk_write_begin(iun, 'gamma_only',new_line=.FALSE.)
            IF (obj%gamma_only) THEN
               WRITE(iun, '(A)',advance='no')  'true'
            ELSE
               WRITE(iun, '(A)',advance='no')  'false'
            ENDIF
         CALL iotk_write_end(iun, 'gamma_only',indentation=.FALSE.)
      ENDIF
      !
      CALL iotk_write_begin(iun, 'ecutwfc')
         WRITE(iun, '(E24.16)') obj%ecutwfc
      CALL iotk_write_end(iun, 'ecutwfc')
      IF(obj%ecutrho_ispresent) THEN
         CALL iotk_write_begin(iun, 'ecutrho')
            WRITE(iun, '(E24.16)') obj%ecutrho
         CALL iotk_write_end(iun, 'ecutrho')
      ENDIF
      !
      IF(obj%fft_grid_ispresent) THEN
         CALL qes_write_basisSetItem(iun, obj%fft_grid)
         !
      ENDIF
      !
      IF(obj%fft_smooth_ispresent) THEN
         CALL qes_write_basisSetItem(iun, obj%fft_smooth)
         !
      ENDIF
      !
      IF(obj%fft_box_ispresent) THEN
         CALL qes_write_basisSetItem(iun, obj%fft_box)
         !
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_basis

SUBROUTINE qes_init_basis(obj, tagname, gamma_only_ispresent, gamma_only, ecutwfc, &
                              ecutrho_ispresent, ecutrho, fft_grid_ispresent, fft_grid, &
                              fft_smooth_ispresent, fft_smooth, fft_box_ispresent, fft_box)
   IMPLICIT NONE

   TYPE(basis_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
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

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%gamma_only_ispresent = gamma_only_ispresent
   IF(obj%gamma_only_ispresent) THEN
      obj%gamma_only = gamma_only
   ENDIF
   obj%ecutwfc = ecutwfc
   obj%ecutrho_ispresent = ecutrho_ispresent
   IF(obj%ecutrho_ispresent) THEN
      obj%ecutrho = ecutrho
   ENDIF
   obj%fft_grid_ispresent = fft_grid_ispresent
   IF(obj%fft_grid_ispresent) THEN
      obj%fft_grid = fft_grid
   ENDIF
   obj%fft_smooth_ispresent = fft_smooth_ispresent
   IF(obj%fft_smooth_ispresent) THEN
      obj%fft_smooth = fft_smooth
   ENDIF
   obj%fft_box_ispresent = fft_box_ispresent
   IF(obj%fft_box_ispresent) THEN
      obj%fft_box = fft_box
   ENDIF

END SUBROUTINE qes_init_basis

SUBROUTINE qes_reset_basis(obj)
   IMPLICIT NONE
   TYPE(basis_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%gamma_only_ispresent) THEN
      obj%gamma_only_ispresent = .FALSE.
   ENDIF
   IF(obj%ecutrho_ispresent) THEN
      obj%ecutrho_ispresent = .FALSE.
   ENDIF
   IF(obj%fft_grid_ispresent) THEN
      CALL qes_reset_basisSetItem(obj%fft_grid)
      obj%fft_grid_ispresent = .FALSE.
   ENDIF
   IF(obj%fft_smooth_ispresent) THEN
      CALL qes_reset_basisSetItem(obj%fft_smooth)
      obj%fft_smooth_ispresent = .FALSE.
   ENDIF
   IF(obj%fft_box_ispresent) THEN
      CALL qes_reset_basisSetItem(obj%fft_box)
      obj%fft_box_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_basis


SUBROUTINE qes_write_inputOccupations(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(inputOccupations_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'ispin', obj%ispin)
   CALL iotk_write_attr(attr, 'spin_factor', obj%spin_factor)

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
         WRITE(fmtstr,'(a)') '(5E24.16)'
         WRITE(iun, fmtstr) obj%vec
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_inputOccupations

SUBROUTINE qes_init_inputOccupations(obj, tagname, ispin, spin_factor, ndim_vec, vec)
   IMPLICIT NONE

   TYPE(inputOccupations_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: ispin
   REAL(DP) :: spin_factor
   INTEGER  :: ndim_vec
   REAL(DP), DIMENSION(ndim_vec) :: vec

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%ispin = ispin


   obj%spin_factor = spin_factor

   ALLOCATE(obj%vec(ndim_vec))
   obj%vec(:) = vec(:)
   obj%ndim_vec = ndim_vec

END SUBROUTINE qes_init_inputOccupations

SUBROUTINE qes_reset_inputOccupations(obj)
   IMPLICIT NONE
   TYPE(inputOccupations_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF (ALLOCATED(obj%vec))  DEALLOCATE(obj%vec)

END SUBROUTINE qes_reset_inputOccupations


SUBROUTINE qes_write_smearing(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(smearing_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'degauss', obj%degauss)

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr),new_line=.FALSE.)
      !
      WRITE(iun, '(A)',advance='no')  TRIM(obj%smearing)
   CALL iotk_write_end(iun, TRIM(obj%tagname),indentation=.FALSE.)
   !
END SUBROUTINE qes_write_smearing

SUBROUTINE qes_init_smearing(obj, tagname, degauss, smearing)
   IMPLICIT NONE

   TYPE(smearing_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   REAL(DP) :: degauss
   CHARACTER(len=*) :: smearing

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%degauss = degauss

   obj%smearing = smearing

END SUBROUTINE qes_init_smearing

SUBROUTINE qes_reset_smearing(obj)
   IMPLICIT NONE
   TYPE(smearing_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_smearing


SUBROUTINE qes_write_band_structure(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(band_structure_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'lsda',new_line=.FALSE.)
         IF (obj%lsda) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'lsda',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'noncolin',new_line=.FALSE.)
         IF (obj%noncolin) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'noncolin',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'spinorbit',new_line=.FALSE.)
         IF (obj%spinorbit) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'spinorbit',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'nbnd')
         WRITE(iun, '(I12)') obj%nbnd
      CALL iotk_write_end(iun, 'nbnd')
      IF(obj%nbnd_up_ispresent) THEN
         CALL iotk_write_begin(iun, 'nbnd_up')
            WRITE(iun, '(I12)') obj%nbnd_up
         CALL iotk_write_end(iun, 'nbnd_up')
      ENDIF
      !
      IF(obj%nbnd_dw_ispresent) THEN
         CALL iotk_write_begin(iun, 'nbnd_dw')
            WRITE(iun, '(I12)') obj%nbnd_dw
         CALL iotk_write_end(iun, 'nbnd_dw')
      ENDIF
      !
      CALL iotk_write_begin(iun, 'nelec')
         WRITE(iun, '(E24.16)') obj%nelec
      CALL iotk_write_end(iun, 'nelec')
      IF(obj%num_of_atomic_wfc_ispresent) THEN
         CALL iotk_write_begin(iun, 'num_of_atomic_wfc')
            WRITE(iun, '(I12)') obj%num_of_atomic_wfc
         CALL iotk_write_end(iun, 'num_of_atomic_wfc')
      ENDIF
      !
      CALL iotk_write_begin(iun, 'wf_collected',new_line=.FALSE.)
         IF (obj%wf_collected) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'wf_collected',indentation=.FALSE.)
      IF(obj%fermi_energy_ispresent) THEN
         CALL iotk_write_begin(iun, 'fermi_energy')
            WRITE(iun, '(E24.16)') obj%fermi_energy
         CALL iotk_write_end(iun, 'fermi_energy')
      ENDIF
      !
      IF(obj%highestOccupiedLevel_ispresent) THEN
         CALL iotk_write_begin(iun, 'highestOccupiedLevel')
            WRITE(iun, '(E24.16)') obj%highestOccupiedLevel
         CALL iotk_write_end(iun, 'highestOccupiedLevel')
      ENDIF
      !
      IF(obj%two_fermi_energies_ispresent) THEN
         CALL iotk_write_begin(iun,'two_fermi_energies')
            WRITE(fmtstr,'(a)') '(5E24.16)'
            WRITE(iun, fmtstr) obj%two_fermi_energies
         CALL iotk_write_end(iun,'two_fermi_energies')
         !
      ENDIF
      !
      CALL qes_write_k_points_IBZ(iun, obj%starting_k_points)
      !
      CALL iotk_write_begin(iun, 'nks')
         WRITE(iun, '(I12)') obj%nks
      CALL iotk_write_end(iun, 'nks')
      CALL qes_write_occupations(iun, obj%occupations_kind)
      !
      IF(obj%smearing_ispresent) THEN
         CALL qes_write_smearing(iun, obj%smearing)
         !
      ENDIF
      !
      DO i = 1, obj%ndim_ks_energies
         CALL qes_write_ks_energies(iun, obj%ks_energies(i))
         !
      END DO
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_band_structure

SUBROUTINE qes_init_band_structure(obj, tagname, lsda, noncolin, spinorbit, nbnd, &
                              nbnd_up_ispresent, nbnd_up, nbnd_dw_ispresent, nbnd_dw, &
                              nelec, num_of_atomic_wfc_ispresent, num_of_atomic_wfc, &
                              wf_collected, fermi_energy_ispresent, fermi_energy, &
                              highestOccupiedLevel_ispresent, highestOccupiedLevel, &
                              two_fermi_energies_ispresent, &
                              ndim_two_fermi_energies, two_fermi_energies, &
                              starting_k_points, nks, occupations_kind, smearing_ispresent, &
                              smearing, ndim_ks_energies, ks_energies)
   IMPLICIT NONE

   TYPE(band_structure_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
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
   LOGICAL  :: wf_collected
   LOGICAL  :: fermi_energy_ispresent
   REAL(DP) :: fermi_energy
   LOGICAL  :: highestOccupiedLevel_ispresent
   REAL(DP) :: highestOccupiedLevel
   LOGICAL  :: two_fermi_energies_ispresent
   INTEGER  :: ndim_two_fermi_energies
   REAL(DP), DIMENSION(ndim_two_fermi_energies) :: two_fermi_energies
   TYPE(k_points_IBZ_type) :: starting_k_points
   INTEGER  :: nks
   TYPE(occupations_type) :: occupations_kind
   LOGICAL  :: smearing_ispresent
   TYPE(smearing_type) :: smearing
   INTEGER  :: ndim_ks_energies
   TYPE(ks_energies_type ), DIMENSION( ndim_ks_energies )  :: ks_energies

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%lsda = lsda
   obj%noncolin = noncolin
   obj%spinorbit = spinorbit
   obj%nbnd = nbnd
   obj%nbnd_up_ispresent = nbnd_up_ispresent
   IF(obj%nbnd_up_ispresent) THEN
      obj%nbnd_up = nbnd_up
   ENDIF
   obj%nbnd_dw_ispresent = nbnd_dw_ispresent
   IF(obj%nbnd_dw_ispresent) THEN
      obj%nbnd_dw = nbnd_dw
   ENDIF
   obj%nelec = nelec
   obj%num_of_atomic_wfc_ispresent = num_of_atomic_wfc_ispresent
   IF(obj%num_of_atomic_wfc_ispresent) THEN
      obj%num_of_atomic_wfc = num_of_atomic_wfc
   ENDIF
   obj%wf_collected = wf_collected
   obj%fermi_energy_ispresent = fermi_energy_ispresent
   IF(obj%fermi_energy_ispresent) THEN
      obj%fermi_energy = fermi_energy
   ENDIF
   obj%highestOccupiedLevel_ispresent = highestOccupiedLevel_ispresent
   IF(obj%highestOccupiedLevel_ispresent) THEN
      obj%highestOccupiedLevel = highestOccupiedLevel
   ENDIF
   obj%two_fermi_energies_ispresent = two_fermi_energies_ispresent
   IF(obj%two_fermi_energies_ispresent) THEN
   ALLOCATE(obj%two_fermi_energies(ndim_two_fermi_energies))
   obj%two_fermi_energies(:) = two_fermi_energies(:)
   obj%ndim_two_fermi_energies = ndim_two_fermi_energies
   ENDIF
   obj%starting_k_points = starting_k_points
   obj%nks = nks
   obj%occupations_kind = occupations_kind
   obj%smearing_ispresent = smearing_ispresent
   IF(obj%smearing_ispresent) THEN
      obj%smearing = smearing
   ENDIF
   ALLOCATE(obj%ks_energies(SIZE(ks_energies)))
   DO i = 1, SIZE(ks_energies)
      obj%ks_energies(i) = ks_energies(i)
   ENDDO
   obj%ndim_ks_energies = ndim_ks_energies

END SUBROUTINE qes_init_band_structure

SUBROUTINE qes_reset_band_structure(obj)
   IMPLICIT NONE
   TYPE(band_structure_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%nbnd_up_ispresent) THEN
      obj%nbnd_up_ispresent = .FALSE.
   ENDIF
   IF(obj%nbnd_dw_ispresent) THEN
      obj%nbnd_dw_ispresent = .FALSE.
   ENDIF
   IF(obj%num_of_atomic_wfc_ispresent) THEN
      obj%num_of_atomic_wfc_ispresent = .FALSE.
   ENDIF
   IF(obj%fermi_energy_ispresent) THEN
      obj%fermi_energy_ispresent = .FALSE.
   ENDIF
   IF(obj%highestOccupiedLevel_ispresent) THEN
      obj%highestOccupiedLevel_ispresent = .FALSE.
   ENDIF
   IF(obj%two_fermi_energies_ispresent) THEN
   IF (ALLOCATED(obj%two_fermi_energies))  DEALLOCATE(obj%two_fermi_energies)
      obj%two_fermi_energies_ispresent = .FALSE.
   ENDIF
   CALL qes_reset_k_points_IBZ(obj%starting_k_points)
   CALL qes_reset_occupations(obj%occupations_kind)
   IF(obj%smearing_ispresent) THEN
      CALL qes_reset_smearing(obj%smearing)
      obj%smearing_ispresent = .FALSE.
   ENDIF
   DO i = 1, SIZE(obj%ks_energies)
      CALL qes_reset_ks_energies(obj%ks_energies(i))
   ENDDO
   IF (ALLOCATED(obj%ks_energies)) DEALLOCATE(obj%ks_energies)

END SUBROUTINE qes_reset_band_structure


SUBROUTINE qes_write_bands(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(bands_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      IF(obj%nbnd_ispresent) THEN
         CALL iotk_write_begin(iun, 'nbnd')
            WRITE(iun, '(I12)') obj%nbnd
         CALL iotk_write_end(iun, 'nbnd')
      ENDIF
      !
      IF(obj%smearing_ispresent) THEN
         CALL qes_write_smearing(iun, obj%smearing)
         !
      ENDIF
      !
      IF(obj%tot_charge_ispresent) THEN
         CALL iotk_write_begin(iun, 'tot_charge')
            WRITE(iun, '(E24.16)') obj%tot_charge
         CALL iotk_write_end(iun, 'tot_charge')
      ENDIF
      !
      IF(obj%tot_magnetization_ispresent) THEN
         CALL iotk_write_begin(iun, 'tot_magnetization')
            WRITE(iun, '(E24.16)') obj%tot_magnetization
         CALL iotk_write_end(iun, 'tot_magnetization')
      ENDIF
      !
      CALL qes_write_occupations(iun, obj%occupations)
      !
      IF(obj%inputOccupations_ispresent) THEN
         DO i = 1, obj%ndim_inputOccupations
            CALL qes_write_inputOccupations(iun, obj%inputOccupations(i))
            !
         END DO
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_bands

SUBROUTINE qes_init_bands(obj, tagname, nbnd_ispresent, nbnd, smearing_ispresent, smearing, &
                              tot_charge_ispresent, tot_charge, &
                              tot_magnetization_ispresent, tot_magnetization, occupations, &
                              inputOccupations_ispresent, ndim_inputOccupations, &
                              inputOccupations)
   IMPLICIT NONE

   TYPE(bands_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
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
   INTEGER  :: ndim_inputOccupations
   TYPE(inputOccupations_type ), DIMENSION( ndim_inputOccupations )  :: inputOccupations

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%nbnd_ispresent = nbnd_ispresent
   IF(obj%nbnd_ispresent) THEN
      obj%nbnd = nbnd
   ENDIF
   obj%smearing_ispresent = smearing_ispresent
   IF(obj%smearing_ispresent) THEN
      obj%smearing = smearing
   ENDIF
   obj%tot_charge_ispresent = tot_charge_ispresent
   IF(obj%tot_charge_ispresent) THEN
      obj%tot_charge = tot_charge
   ENDIF
   obj%tot_magnetization_ispresent = tot_magnetization_ispresent
   IF(obj%tot_magnetization_ispresent) THEN
      obj%tot_magnetization = tot_magnetization
   ENDIF
   obj%occupations = occupations
   obj%inputOccupations_ispresent = inputOccupations_ispresent
   IF(obj%inputOccupations_ispresent) THEN
      ALLOCATE(obj%inputOccupations(SIZE(inputOccupations)))
      DO i = 1, SIZE(inputOccupations)
         obj%inputOccupations(i) = inputOccupations(i)
      ENDDO
   obj%ndim_inputOccupations = ndim_inputOccupations
   ENDIF

END SUBROUTINE qes_init_bands

SUBROUTINE qes_reset_bands(obj)
   IMPLICIT NONE
   TYPE(bands_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%nbnd_ispresent) THEN
      obj%nbnd_ispresent = .FALSE.
   ENDIF
   IF(obj%smearing_ispresent) THEN
      CALL qes_reset_smearing(obj%smearing)
      obj%smearing_ispresent = .FALSE.
   ENDIF
   IF(obj%tot_charge_ispresent) THEN
      obj%tot_charge_ispresent = .FALSE.
   ENDIF
   IF(obj%tot_magnetization_ispresent) THEN
      obj%tot_magnetization_ispresent = .FALSE.
   ENDIF
   CALL qes_reset_occupations(obj%occupations)
   IF(obj%inputOccupations_ispresent) THEN
      DO i = 1, SIZE(obj%inputOccupations)
         CALL qes_reset_inputOccupations(obj%inputOccupations(i))
      ENDDO
      IF (ALLOCATED(obj%inputOccupations)) DEALLOCATE(obj%inputOccupations)
      obj%inputOccupations_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_bands


SUBROUTINE qes_write_spin(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(spin_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'lsda',new_line=.FALSE.)
         IF (obj%lsda) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'lsda',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'noncolin',new_line=.FALSE.)
         IF (obj%noncolin) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'noncolin',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'spinorbit',new_line=.FALSE.)
         IF (obj%spinorbit) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'spinorbit',indentation=.FALSE.)
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_spin

SUBROUTINE qes_init_spin(obj, tagname, lsda, noncolin, spinorbit)
   IMPLICIT NONE

   TYPE(spin_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   LOGICAL  :: lsda
   LOGICAL  :: noncolin
   LOGICAL  :: spinorbit

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%lsda = lsda
   obj%noncolin = noncolin
   obj%spinorbit = spinorbit

END SUBROUTINE qes_init_spin

SUBROUTINE qes_reset_spin(obj)
   IMPLICIT NONE
   TYPE(spin_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_spin


SUBROUTINE qes_write_HubbardCommon(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(HubbardCommon_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'specie', TRIM(obj%specie))
   CALL iotk_write_attr(attr, 'label', TRIM(obj%label))

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      WRITE(iun, '(E24.16)') obj%HubbardCommon
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_HubbardCommon

SUBROUTINE qes_init_HubbardCommon(obj, tagname, specie, label, HubbardCommon)
   IMPLICIT NONE

   TYPE(HubbardCommon_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: specie
   CHARACTER(len=*) :: label
   REAL(DP) :: HubbardCommon

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%specie = TRIM(specie)


   obj%label = TRIM(label)

   obj%HubbardCommon = HubbardCommon

END SUBROUTINE qes_init_HubbardCommon

SUBROUTINE qes_reset_HubbardCommon(obj)
   IMPLICIT NONE
   TYPE(HubbardCommon_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_HubbardCommon


SUBROUTINE qes_write_HubbardProj(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(HubbardProj_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_HubbardProj

SUBROUTINE qes_init_HubbardProj(obj, tagname)
   IMPLICIT NONE

   TYPE(HubbardProj_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

END SUBROUTINE qes_init_HubbardProj

SUBROUTINE qes_reset_HubbardProj(obj)
   IMPLICIT NONE
   TYPE(HubbardProj_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_HubbardProj


SUBROUTINE qes_write_Hubbard_ns(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(Hubbard_ns_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'specie', TRIM(obj%specie))
   CALL iotk_write_attr(attr, 'label', TRIM(obj%label))
   CALL iotk_write_attr(attr, 'spin', obj%spin)
   CALL iotk_write_attr(attr, 'index', obj%index)

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      WRITE(fmtstr,'(a)') '(5E24.16)'
      DO i = 1, SIZE(obj%mat,2)
         WRITE(iun, fmtstr) obj%mat(:,i)
      ENDDO
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_Hubbard_ns

SUBROUTINE qes_init_Hubbard_ns(obj, tagname, specie, label, spin, index, &
                              ndim1_mat, ndim2_mat, mat)
   IMPLICIT NONE

   TYPE(Hubbard_ns_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: specie
   CHARACTER(len=*) :: label
   INTEGER  :: spin
   INTEGER  :: index
   INTEGER  :: ndim1_mat
   INTEGER  :: ndim2_mat
   REAL(DP), DIMENSION(ndim1_mat,ndim2_mat) :: mat

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%specie = TRIM(specie)


   obj%label = TRIM(label)


   obj%spin = spin


   obj%index = index

   ALLOCATE(obj%mat(ndim1_mat,ndim2_mat))
   obj%mat(:,:) = mat(:,:)
   obj%ndim1_mat = ndim1_mat
   obj%ndim2_mat = ndim2_mat

END SUBROUTINE qes_init_Hubbard_ns

SUBROUTINE qes_reset_Hubbard_ns(obj)
   IMPLICIT NONE
   TYPE(Hubbard_ns_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF (ALLOCATED(obj%mat))  DEALLOCATE(obj%mat)

END SUBROUTINE qes_reset_Hubbard_ns


SUBROUTINE qes_write_starting_ns(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(starting_ns_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'specie', TRIM(obj%specie))
   CALL iotk_write_attr(attr, 'label', TRIM(obj%label))
   CALL iotk_write_attr(attr, 'spin', obj%spin)

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
         WRITE(fmtstr,'(a)') '(5E24.16)'
         WRITE(iun, fmtstr) obj%vec
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_starting_ns

SUBROUTINE qes_init_starting_ns(obj, tagname, specie, label, spin, ndim_vec, vec)
   IMPLICIT NONE

   TYPE(starting_ns_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: specie
   CHARACTER(len=*) :: label
   INTEGER  :: spin
   INTEGER  :: ndim_vec
   REAL(DP), DIMENSION(ndim_vec) :: vec

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%specie = TRIM(specie)


   obj%label = TRIM(label)


   obj%spin = spin

   ALLOCATE(obj%vec(ndim_vec))
   obj%vec(:) = vec(:)
   obj%ndim_vec = ndim_vec

END SUBROUTINE qes_init_starting_ns

SUBROUTINE qes_reset_starting_ns(obj)
   IMPLICIT NONE
   TYPE(starting_ns_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF (ALLOCATED(obj%vec))  DEALLOCATE(obj%vec)

END SUBROUTINE qes_reset_starting_ns


SUBROUTINE qes_write_HubbardJ(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(HubbardJ_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'specie', TRIM(obj%specie))
   CALL iotk_write_attr(attr, 'label', TRIM(obj%label))

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      WRITE(iun, '(3(E24.16))') obj%HubbardJ
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_HubbardJ

SUBROUTINE qes_init_HubbardJ(obj, tagname, specie, label, HubbardJ)
   IMPLICIT NONE

   TYPE(HubbardJ_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: specie
   CHARACTER(len=*) :: label
   REAL(DP), DIMENSION(3) :: HubbardJ

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%specie = TRIM(specie)


   obj%label = TRIM(label)

   obj%HubbardJ = HubbardJ

END SUBROUTINE qes_init_HubbardJ

SUBROUTINE qes_reset_HubbardJ(obj)
   IMPLICIT NONE
   TYPE(HubbardJ_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_HubbardJ


SUBROUTINE qes_write_vdW(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(vdW_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'vdw_corr',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%vdw_corr)
      CALL iotk_write_end(iun, 'vdw_corr',indentation=.FALSE.)
      IF(obj%non_local_term_ispresent) THEN
         CALL iotk_write_begin(iun, 'non_local_term',new_line=.FALSE.)
            WRITE(iun, '(A)',advance='no')  TRIM(obj%non_local_term)
         CALL iotk_write_end(iun, 'non_local_term',indentation=.FALSE.)
      ENDIF
      !
      IF(obj%london_s6_ispresent) THEN
         CALL iotk_write_begin(iun, 'london_s6')
            WRITE(iun, '(E24.16)') obj%london_s6
         CALL iotk_write_end(iun, 'london_s6')
      ENDIF
      !
      IF(obj%ts_vdw_econv_thr_ispresent) THEN
         CALL iotk_write_begin(iun, 'ts_vdw_econv_thr')
            WRITE(iun, '(E24.16)') obj%ts_vdw_econv_thr
         CALL iotk_write_end(iun, 'ts_vdw_econv_thr')
      ENDIF
      !
      IF(obj%ts_vdw_isolated_ispresent) THEN
         CALL iotk_write_begin(iun, 'ts_vdw_isolated',new_line=.FALSE.)
            IF (obj%ts_vdw_isolated) THEN
               WRITE(iun, '(A)',advance='no')  'true'
            ELSE
               WRITE(iun, '(A)',advance='no')  'false'
            ENDIF
         CALL iotk_write_end(iun, 'ts_vdw_isolated',indentation=.FALSE.)
      ENDIF
      !
      IF(obj%london_rcut_ispresent) THEN
         CALL iotk_write_begin(iun, 'london_rcut')
            WRITE(iun, '(E24.16)') obj%london_rcut
         CALL iotk_write_end(iun, 'london_rcut')
      ENDIF
      !
      IF(obj%xdm_a1_ispresent) THEN
         CALL iotk_write_begin(iun, 'xdm_a1')
            WRITE(iun, '(E24.16)') obj%xdm_a1
         CALL iotk_write_end(iun, 'xdm_a1')
      ENDIF
      !
      IF(obj%xdm_a2_ispresent) THEN
         CALL iotk_write_begin(iun, 'xdm_a2')
            WRITE(iun, '(E24.16)') obj%xdm_a2
         CALL iotk_write_end(iun, 'xdm_a2')
      ENDIF
      !
      IF(obj%london_c6_ispresent) THEN
         DO i = 1, obj%ndim_london_c6
            CALL qes_write_HubbardCommon(iun, obj%london_c6(i))
            !
         END DO
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_vdW

SUBROUTINE qes_init_vdW(obj, tagname, vdw_corr, non_local_term_ispresent, non_local_term, &
                              london_s6_ispresent, london_s6, ts_vdw_econv_thr_ispresent, &
                              ts_vdw_econv_thr, ts_vdw_isolated_ispresent, ts_vdw_isolated, &
                              london_rcut_ispresent, london_rcut, xdm_a1_ispresent, xdm_a1, &
                              xdm_a2_ispresent, xdm_a2, london_c6_ispresent, &
                              ndim_london_c6, london_c6)
   IMPLICIT NONE

   TYPE(vdW_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: vdw_corr
   LOGICAL  :: non_local_term_ispresent
   CHARACTER(len=*) :: non_local_term
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
   INTEGER  :: ndim_london_c6
   TYPE(HubbardCommon_type ), DIMENSION( ndim_london_c6 )  :: london_c6

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%vdw_corr = vdw_corr
   obj%non_local_term_ispresent = non_local_term_ispresent
   IF(obj%non_local_term_ispresent) THEN
      obj%non_local_term = non_local_term
   ENDIF
   obj%london_s6_ispresent = london_s6_ispresent
   IF(obj%london_s6_ispresent) THEN
      obj%london_s6 = london_s6
   ENDIF
   obj%ts_vdw_econv_thr_ispresent = ts_vdw_econv_thr_ispresent
   IF(obj%ts_vdw_econv_thr_ispresent) THEN
      obj%ts_vdw_econv_thr = ts_vdw_econv_thr
   ENDIF
   obj%ts_vdw_isolated_ispresent = ts_vdw_isolated_ispresent
   IF(obj%ts_vdw_isolated_ispresent) THEN
      obj%ts_vdw_isolated = ts_vdw_isolated
   ENDIF
   obj%london_rcut_ispresent = london_rcut_ispresent
   IF(obj%london_rcut_ispresent) THEN
      obj%london_rcut = london_rcut
   ENDIF
   obj%xdm_a1_ispresent = xdm_a1_ispresent
   IF(obj%xdm_a1_ispresent) THEN
      obj%xdm_a1 = xdm_a1
   ENDIF
   obj%xdm_a2_ispresent = xdm_a2_ispresent
   IF(obj%xdm_a2_ispresent) THEN
      obj%xdm_a2 = xdm_a2
   ENDIF
   obj%london_c6_ispresent = london_c6_ispresent
   IF(obj%london_c6_ispresent) THEN
      ALLOCATE(obj%london_c6(SIZE(london_c6)))
      DO i = 1, SIZE(london_c6)
         obj%london_c6(i) = london_c6(i)
      ENDDO
   obj%ndim_london_c6 = ndim_london_c6
   ENDIF

END SUBROUTINE qes_init_vdW

SUBROUTINE qes_reset_vdW(obj)
   IMPLICIT NONE
   TYPE(vdW_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%non_local_term_ispresent) THEN
      obj%non_local_term_ispresent = .FALSE.
   ENDIF
   IF(obj%london_s6_ispresent) THEN
      obj%london_s6_ispresent = .FALSE.
   ENDIF
   IF(obj%ts_vdw_econv_thr_ispresent) THEN
      obj%ts_vdw_econv_thr_ispresent = .FALSE.
   ENDIF
   IF(obj%ts_vdw_isolated_ispresent) THEN
      obj%ts_vdw_isolated_ispresent = .FALSE.
   ENDIF
   IF(obj%london_rcut_ispresent) THEN
      obj%london_rcut_ispresent = .FALSE.
   ENDIF
   IF(obj%xdm_a1_ispresent) THEN
      obj%xdm_a1_ispresent = .FALSE.
   ENDIF
   IF(obj%xdm_a2_ispresent) THEN
      obj%xdm_a2_ispresent = .FALSE.
   ENDIF
   IF(obj%london_c6_ispresent) THEN
      DO i = 1, SIZE(obj%london_c6)
         CALL qes_reset_HubbardCommon(obj%london_c6(i))
      ENDDO
      IF (ALLOCATED(obj%london_c6)) DEALLOCATE(obj%london_c6)
      obj%london_c6_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_vdW


SUBROUTINE qes_write_dftU(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(dftU_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      IF(obj%lda_plus_u_kind_ispresent) THEN
         CALL iotk_write_begin(iun, 'lda_plus_u_kind',new_line=.FALSE.)
            WRITE(iun, '(I0)',advance='no') obj%lda_plus_u_kind
         CALL iotk_write_end(iun, 'lda_plus_u_kind',indentation=.FALSE.)
      ENDIF
      !
      IF(obj%Hubbard_U_ispresent) THEN
         DO i = 1, obj%ndim_Hubbard_U
            CALL qes_write_HubbardCommon(iun, obj%Hubbard_U(i))
            !
         END DO
      ENDIF
      !
      IF(obj%Hubbard_J0_ispresent) THEN
         DO i = 1, obj%ndim_Hubbard_J0
            CALL qes_write_HubbardCommon(iun, obj%Hubbard_J0(i))
            !
         END DO
      ENDIF
      !
      IF(obj%Hubbard_alpha_ispresent) THEN
         DO i = 1, obj%ndim_Hubbard_alpha
            CALL qes_write_HubbardCommon(iun, obj%Hubbard_alpha(i))
            !
         END DO
      ENDIF
      !
      IF(obj%Hubbard_beta_ispresent) THEN
         DO i = 1, obj%ndim_Hubbard_beta
            CALL qes_write_HubbardCommon(iun, obj%Hubbard_beta(i))
            !
         END DO
      ENDIF
      !
      IF(obj%Hubbard_J_ispresent) THEN
         DO i = 1, obj%ndim_Hubbard_J
            CALL qes_write_HubbardJ(iun, obj%Hubbard_J(i))
            !
         END DO
      ENDIF
      !
      IF(obj%starting_ns_ispresent) THEN
         DO i = 1, obj%ndim_starting_ns
            CALL qes_write_starting_ns(iun, obj%starting_ns(i))
            !
         END DO
      ENDIF
      !
      IF(obj%Hubbard_ns_ispresent) THEN
         DO i = 1, obj%ndim_Hubbard_ns
            CALL qes_write_Hubbard_ns(iun, obj%Hubbard_ns(i))
            !
         END DO
      ENDIF
      !
      IF(obj%U_projection_type_ispresent) THEN
         CALL iotk_write_begin(iun, 'U_projection_type',new_line=.FALSE.)
            WRITE(iun, '(A)',advance='no')  TRIM(obj%U_projection_type)
         CALL iotk_write_end(iun, 'U_projection_type',indentation=.FALSE.)
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_dftU

SUBROUTINE qes_init_dftU(obj, tagname, lda_plus_u_kind_ispresent, lda_plus_u_kind, &
                              Hubbard_U_ispresent, ndim_Hubbard_U, Hubbard_U, &
                              Hubbard_J0_ispresent, ndim_Hubbard_J0, Hubbard_J0, &
                              Hubbard_alpha_ispresent, ndim_Hubbard_alpha, Hubbard_alpha, &
                              Hubbard_beta_ispresent, ndim_Hubbard_beta, Hubbard_beta, &
                              Hubbard_J_ispresent, ndim_Hubbard_J, Hubbard_J, &
                              starting_ns_ispresent, ndim_starting_ns, starting_ns, &
                              Hubbard_ns_ispresent, ndim_Hubbard_ns, Hubbard_ns, &
                              U_projection_type_ispresent, U_projection_type)
   IMPLICIT NONE

   TYPE(dftU_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   LOGICAL  :: lda_plus_u_kind_ispresent
   INTEGER :: lda_plus_u_kind
   LOGICAL  :: Hubbard_U_ispresent
   INTEGER  :: ndim_Hubbard_U
   TYPE(HubbardCommon_type ), DIMENSION( ndim_Hubbard_U )  :: Hubbard_U
   LOGICAL  :: Hubbard_J0_ispresent
   INTEGER  :: ndim_Hubbard_J0
   TYPE(HubbardCommon_type ), DIMENSION( ndim_Hubbard_J0 )  :: Hubbard_J0
   LOGICAL  :: Hubbard_alpha_ispresent
   INTEGER  :: ndim_Hubbard_alpha
   TYPE(HubbardCommon_type ), DIMENSION( ndim_Hubbard_alpha )  :: Hubbard_alpha
   LOGICAL  :: Hubbard_beta_ispresent
   INTEGER  :: ndim_Hubbard_beta
   TYPE(HubbardCommon_type ), DIMENSION( ndim_Hubbard_beta )  :: Hubbard_beta
   LOGICAL  :: Hubbard_J_ispresent
   INTEGER  :: ndim_Hubbard_J
   TYPE(HubbardJ_type ), DIMENSION( ndim_Hubbard_J )  :: Hubbard_J
   LOGICAL  :: starting_ns_ispresent
   INTEGER  :: ndim_starting_ns
   TYPE(starting_ns_type ), DIMENSION( ndim_starting_ns )  :: starting_ns
   LOGICAL  :: Hubbard_ns_ispresent
   INTEGER  :: ndim_Hubbard_ns
   TYPE(Hubbard_ns_type ), DIMENSION( ndim_Hubbard_ns )  :: Hubbard_ns
   LOGICAL  :: U_projection_type_ispresent
   CHARACTER(len=*) :: U_projection_type

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%lda_plus_u_kind_ispresent = lda_plus_u_kind_ispresent
   IF(obj%lda_plus_u_kind_ispresent) THEN
      obj%lda_plus_u_kind = lda_plus_u_kind
   ENDIF
   obj%Hubbard_U_ispresent = Hubbard_U_ispresent
   IF(obj%Hubbard_U_ispresent) THEN
      ALLOCATE(obj%Hubbard_U(SIZE(Hubbard_U)))
      DO i = 1, SIZE(Hubbard_U)
         obj%Hubbard_U(i) = Hubbard_U(i)
      ENDDO
   obj%ndim_Hubbard_U = ndim_Hubbard_U
   ENDIF
   obj%Hubbard_J0_ispresent = Hubbard_J0_ispresent
   IF(obj%Hubbard_J0_ispresent) THEN
      ALLOCATE(obj%Hubbard_J0(SIZE(Hubbard_J0)))
      DO i = 1, SIZE(Hubbard_J0)
         obj%Hubbard_J0(i) = Hubbard_J0(i)
      ENDDO
   obj%ndim_Hubbard_J0 = ndim_Hubbard_J0
   ENDIF
   obj%Hubbard_alpha_ispresent = Hubbard_alpha_ispresent
   IF(obj%Hubbard_alpha_ispresent) THEN
      ALLOCATE(obj%Hubbard_alpha(SIZE(Hubbard_alpha)))
      DO i = 1, SIZE(Hubbard_alpha)
         obj%Hubbard_alpha(i) = Hubbard_alpha(i)
      ENDDO
   obj%ndim_Hubbard_alpha = ndim_Hubbard_alpha
   ENDIF
   obj%Hubbard_beta_ispresent = Hubbard_beta_ispresent
   IF(obj%Hubbard_beta_ispresent) THEN
      ALLOCATE(obj%Hubbard_beta(SIZE(Hubbard_beta)))
      DO i = 1, SIZE(Hubbard_beta)
         obj%Hubbard_beta(i) = Hubbard_beta(i)
      ENDDO
   obj%ndim_Hubbard_beta = ndim_Hubbard_beta
   ENDIF
   obj%Hubbard_J_ispresent = Hubbard_J_ispresent
   IF(obj%Hubbard_J_ispresent) THEN
      ALLOCATE(obj%Hubbard_J(SIZE(Hubbard_J)))
      DO i = 1, SIZE(Hubbard_J)
         obj%Hubbard_J(i) = Hubbard_J(i)
      ENDDO
   obj%ndim_Hubbard_J = ndim_Hubbard_J
   ENDIF
   obj%starting_ns_ispresent = starting_ns_ispresent
   IF(obj%starting_ns_ispresent) THEN
      ALLOCATE(obj%starting_ns(SIZE(starting_ns)))
      DO i = 1, SIZE(starting_ns)
         obj%starting_ns(i) = starting_ns(i)
      ENDDO
   obj%ndim_starting_ns = ndim_starting_ns
   ENDIF
   obj%Hubbard_ns_ispresent = Hubbard_ns_ispresent
   IF(obj%Hubbard_ns_ispresent) THEN
      ALLOCATE(obj%Hubbard_ns(SIZE(Hubbard_ns)))
      DO i = 1, SIZE(Hubbard_ns)
         obj%Hubbard_ns(i) = Hubbard_ns(i)
      ENDDO
   obj%ndim_Hubbard_ns = ndim_Hubbard_ns
   ENDIF
   obj%U_projection_type_ispresent = U_projection_type_ispresent
   IF(obj%U_projection_type_ispresent) THEN
      obj%U_projection_type = U_projection_type
   ENDIF

END SUBROUTINE qes_init_dftU

SUBROUTINE qes_reset_dftU(obj)
   IMPLICIT NONE
   TYPE(dftU_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%lda_plus_u_kind_ispresent) THEN
      obj%lda_plus_u_kind_ispresent = .FALSE.
   ENDIF
   IF(obj%Hubbard_U_ispresent) THEN
      DO i = 1, SIZE(obj%Hubbard_U)
         CALL qes_reset_HubbardCommon(obj%Hubbard_U(i))
      ENDDO
      IF (ALLOCATED(obj%Hubbard_U)) DEALLOCATE(obj%Hubbard_U)
      obj%Hubbard_U_ispresent = .FALSE.
   ENDIF
   IF(obj%Hubbard_J0_ispresent) THEN
      DO i = 1, SIZE(obj%Hubbard_J0)
         CALL qes_reset_HubbardCommon(obj%Hubbard_J0(i))
      ENDDO
      IF (ALLOCATED(obj%Hubbard_J0)) DEALLOCATE(obj%Hubbard_J0)
      obj%Hubbard_J0_ispresent = .FALSE.
   ENDIF
   IF(obj%Hubbard_alpha_ispresent) THEN
      DO i = 1, SIZE(obj%Hubbard_alpha)
         CALL qes_reset_HubbardCommon(obj%Hubbard_alpha(i))
      ENDDO
      IF (ALLOCATED(obj%Hubbard_alpha)) DEALLOCATE(obj%Hubbard_alpha)
      obj%Hubbard_alpha_ispresent = .FALSE.
   ENDIF
   IF(obj%Hubbard_beta_ispresent) THEN
      DO i = 1, SIZE(obj%Hubbard_beta)
         CALL qes_reset_HubbardCommon(obj%Hubbard_beta(i))
      ENDDO
      IF (ALLOCATED(obj%Hubbard_beta)) DEALLOCATE(obj%Hubbard_beta)
      obj%Hubbard_beta_ispresent = .FALSE.
   ENDIF
   IF(obj%Hubbard_J_ispresent) THEN
      DO i = 1, SIZE(obj%Hubbard_J)
         CALL qes_reset_HubbardJ(obj%Hubbard_J(i))
      ENDDO
      IF (ALLOCATED(obj%Hubbard_J)) DEALLOCATE(obj%Hubbard_J)
      obj%Hubbard_J_ispresent = .FALSE.
   ENDIF
   IF(obj%starting_ns_ispresent) THEN
      DO i = 1, SIZE(obj%starting_ns)
         CALL qes_reset_starting_ns(obj%starting_ns(i))
      ENDDO
      IF (ALLOCATED(obj%starting_ns)) DEALLOCATE(obj%starting_ns)
      obj%starting_ns_ispresent = .FALSE.
   ENDIF
   IF(obj%Hubbard_ns_ispresent) THEN
      DO i = 1, SIZE(obj%Hubbard_ns)
         CALL qes_reset_Hubbard_ns(obj%Hubbard_ns(i))
      ENDDO
      IF (ALLOCATED(obj%Hubbard_ns)) DEALLOCATE(obj%Hubbard_ns)
      obj%Hubbard_ns_ispresent = .FALSE.
   ENDIF
   IF(obj%U_projection_type_ispresent) THEN
      obj%U_projection_type_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_dftU


SUBROUTINE qes_write_qpoint_grid(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(qpoint_grid_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'nqx1', obj%nqx1)
   CALL iotk_write_attr(attr, 'nqx2', obj%nqx2)
   CALL iotk_write_attr(attr, 'nqx3', obj%nqx3)

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr),new_line=.FALSE.)
      !
      WRITE(iun, '(A)',advance='no')  TRIM(obj%qpoint_grid)
   CALL iotk_write_end(iun, TRIM(obj%tagname),indentation=.FALSE.)
   !
END SUBROUTINE qes_write_qpoint_grid

SUBROUTINE qes_init_qpoint_grid(obj, tagname, nqx1, nqx2, nqx3, qpoint_grid)
   IMPLICIT NONE

   TYPE(qpoint_grid_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: nqx1
   INTEGER  :: nqx2
   INTEGER  :: nqx3
   CHARACTER(len=*) :: qpoint_grid

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%nqx1 = nqx1


   obj%nqx2 = nqx2


   obj%nqx3 = nqx3

   obj%qpoint_grid = qpoint_grid

END SUBROUTINE qes_init_qpoint_grid

SUBROUTINE qes_reset_qpoint_grid(obj)
   IMPLICIT NONE
   TYPE(qpoint_grid_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_qpoint_grid


SUBROUTINE qes_write_hybrid(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(hybrid_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL qes_write_qpoint_grid(iun, obj%qpoint_grid)
      !
      CALL iotk_write_begin(iun, 'ecutfock')
         WRITE(iun, '(E24.16)') obj%ecutfock
      CALL iotk_write_end(iun, 'ecutfock')
      CALL iotk_write_begin(iun, 'exx_fraction')
         WRITE(iun, '(E24.16)') obj%exx_fraction
      CALL iotk_write_end(iun, 'exx_fraction')
      CALL iotk_write_begin(iun, 'screening_parameter')
         WRITE(iun, '(E24.16)') obj%screening_parameter
      CALL iotk_write_end(iun, 'screening_parameter')
      CALL iotk_write_begin(iun, 'exxdiv_treatment',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%exxdiv_treatment)
      CALL iotk_write_end(iun, 'exxdiv_treatment',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'x_gamma_extrapolation',new_line=.FALSE.)
         IF (obj%x_gamma_extrapolation) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'x_gamma_extrapolation',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'ecutvcut')
         WRITE(iun, '(E24.16)') obj%ecutvcut
      CALL iotk_write_end(iun, 'ecutvcut')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_hybrid

SUBROUTINE qes_init_hybrid(obj, tagname, qpoint_grid, ecutfock, exx_fraction, &
                              screening_parameter, exxdiv_treatment, x_gamma_extrapolation, &
                              ecutvcut)
   IMPLICIT NONE

   TYPE(hybrid_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   TYPE(qpoint_grid_type) :: qpoint_grid
   REAL(DP) :: ecutfock
   REAL(DP) :: exx_fraction
   REAL(DP) :: screening_parameter
   CHARACTER(len=*) :: exxdiv_treatment
   LOGICAL  :: x_gamma_extrapolation
   REAL(DP) :: ecutvcut

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%qpoint_grid = qpoint_grid
   obj%ecutfock = ecutfock
   obj%exx_fraction = exx_fraction
   obj%screening_parameter = screening_parameter
   obj%exxdiv_treatment = exxdiv_treatment
   obj%x_gamma_extrapolation = x_gamma_extrapolation
   obj%ecutvcut = ecutvcut

END SUBROUTINE qes_init_hybrid

SUBROUTINE qes_reset_hybrid(obj)
   IMPLICIT NONE
   TYPE(hybrid_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   CALL qes_reset_qpoint_grid(obj%qpoint_grid)

END SUBROUTINE qes_reset_hybrid


SUBROUTINE qes_write_functional(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(functional_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_functional

SUBROUTINE qes_init_functional(obj, tagname)
   IMPLICIT NONE

   TYPE(functional_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

END SUBROUTINE qes_init_functional

SUBROUTINE qes_reset_functional(obj)
   IMPLICIT NONE
   TYPE(functional_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_functional


SUBROUTINE qes_write_dft(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(dft_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'functional',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%functional)
      CALL iotk_write_end(iun, 'functional',indentation=.FALSE.)
      IF(obj%hybrid_ispresent) THEN
         CALL qes_write_hybrid(iun, obj%hybrid)
         !
      ENDIF
      !
      IF(obj%dftU_ispresent) THEN
         CALL qes_write_dftU(iun, obj%dftU)
         !
      ENDIF
      !
      IF(obj%vdW_ispresent) THEN
         CALL qes_write_vdW(iun, obj%vdW)
         !
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_dft

SUBROUTINE qes_init_dft(obj, tagname, functional, hybrid_ispresent, hybrid, dftU_ispresent, &
                              dftU, vdW_ispresent, vdW)
   IMPLICIT NONE

   TYPE(dft_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: functional
   LOGICAL  :: hybrid_ispresent
   TYPE(hybrid_type) :: hybrid
   LOGICAL  :: dftU_ispresent
   TYPE(dftU_type) :: dftU
   LOGICAL  :: vdW_ispresent
   TYPE(vdW_type) :: vdW

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%functional = functional
   obj%hybrid_ispresent = hybrid_ispresent
   IF(obj%hybrid_ispresent) THEN
      obj%hybrid = hybrid
   ENDIF
   obj%dftU_ispresent = dftU_ispresent
   IF(obj%dftU_ispresent) THEN
      obj%dftU = dftU
   ENDIF
   obj%vdW_ispresent = vdW_ispresent
   IF(obj%vdW_ispresent) THEN
      obj%vdW = vdW
   ENDIF

END SUBROUTINE qes_init_dft

SUBROUTINE qes_reset_dft(obj)
   IMPLICIT NONE
   TYPE(dft_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%hybrid_ispresent) THEN
      CALL qes_reset_hybrid(obj%hybrid)
      obj%hybrid_ispresent = .FALSE.
   ENDIF
   IF(obj%dftU_ispresent) THEN
      CALL qes_reset_dftU(obj%dftU)
      obj%dftU_ispresent = .FALSE.
   ENDIF
   IF(obj%vdW_ispresent) THEN
      CALL qes_reset_vdW(obj%vdW)
      obj%vdW_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_dft


SUBROUTINE qes_write_d3vector(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(d3vector_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_d3vector

SUBROUTINE qes_init_d3vector(obj, tagname)
   IMPLICIT NONE

   TYPE(d3vector_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

END SUBROUTINE qes_init_d3vector

SUBROUTINE qes_reset_d3vector(obj)
   IMPLICIT NONE
   TYPE(d3vector_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_d3vector


SUBROUTINE qes_write_cell(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(cell_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'a1')
         WRITE(iun, '(3(E24.16))') obj%a1
      CALL iotk_write_end(iun, 'a1')
      CALL iotk_write_begin(iun, 'a2')
         WRITE(iun, '(3(E24.16))') obj%a2
      CALL iotk_write_end(iun, 'a2')
      CALL iotk_write_begin(iun, 'a3')
         WRITE(iun, '(3(E24.16))') obj%a3
      CALL iotk_write_end(iun, 'a3')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_cell

SUBROUTINE qes_init_cell(obj, tagname, a1, a2, a3)
   IMPLICIT NONE

   TYPE(cell_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   REAL(DP), DIMENSION(3) :: a1
   REAL(DP), DIMENSION(3) :: a2
   REAL(DP), DIMENSION(3) :: a3

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%a1 = a1
   obj%a2 = a2
   obj%a3 = a3

END SUBROUTINE qes_init_cell

SUBROUTINE qes_reset_cell(obj)
   IMPLICIT NONE
   TYPE(cell_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_cell


SUBROUTINE qes_write_wyckoff_positions(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(wyckoff_positions_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'space_group', obj%space_group)
   IF(obj%more_options_ispresent) THEN
      CALL iotk_write_attr(attr, 'more_options', TRIM(obj%more_options))
   END IF

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      DO i = 1, obj%ndim_atom
         CALL qes_write_atom(iun, obj%atom(i))
         !
      END DO
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_wyckoff_positions

SUBROUTINE qes_init_wyckoff_positions(obj, tagname, space_group, more_options, more_options_ispresent, &
                              ndim_atom, atom)
   IMPLICIT NONE

   TYPE(wyckoff_positions_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER :: space_group
   LOGICAL  :: more_options_ispresent
   CHARACTER(len=*), OPTIONAL :: more_options
   INTEGER  :: ndim_atom
   TYPE(atom_type ), DIMENSION( ndim_atom )  :: atom

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%space_group = space_group


   obj%more_options_ispresent = more_options_ispresent
   IF (obj%more_options_ispresent) THEN
      obj%more_options = TRIM(more_options)
   ENDIF

   ALLOCATE(obj%atom(SIZE(atom)))
   DO i = 1, SIZE(atom)
      obj%atom(i) = atom(i)
   ENDDO
   obj%ndim_atom = ndim_atom

END SUBROUTINE qes_init_wyckoff_positions

SUBROUTINE qes_reset_wyckoff_positions(obj)
   IMPLICIT NONE
   TYPE(wyckoff_positions_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   DO i = 1, SIZE(obj%atom)
      CALL qes_reset_atom(obj%atom(i))
   ENDDO
   IF (ALLOCATED(obj%atom)) DEALLOCATE(obj%atom)

END SUBROUTINE qes_reset_wyckoff_positions


SUBROUTINE qes_write_atomic_positions(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(atomic_positions_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      DO i = 1, obj%ndim_atom
         CALL qes_write_atom(iun, obj%atom(i))
         !
      END DO
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_atomic_positions

SUBROUTINE qes_init_atomic_positions(obj, tagname, ndim_atom, atom)
   IMPLICIT NONE

   TYPE(atomic_positions_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: ndim_atom
   TYPE(atom_type ), DIMENSION( ndim_atom )  :: atom

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   ALLOCATE(obj%atom(SIZE(atom)))
   DO i = 1, SIZE(atom)
      obj%atom(i) = atom(i)
   ENDDO
   obj%ndim_atom = ndim_atom

END SUBROUTINE qes_init_atomic_positions

SUBROUTINE qes_reset_atomic_positions(obj)
   IMPLICIT NONE
   TYPE(atomic_positions_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   DO i = 1, SIZE(obj%atom)
      CALL qes_reset_atom(obj%atom(i))
   ENDDO
   IF (ALLOCATED(obj%atom)) DEALLOCATE(obj%atom)

END SUBROUTINE qes_reset_atomic_positions


SUBROUTINE qes_write_atomic_structure(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(atomic_structure_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'nat', obj%nat)
   IF(obj%alat_ispresent) THEN
      CALL iotk_write_attr(attr, 'alat', obj%alat)
   END IF
   IF(obj%bravais_index_ispresent) THEN
      CALL iotk_write_attr(attr, 'bravais_index', obj%bravais_index)
   END IF

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      IF(obj%atomic_positions_ispresent) THEN
         CALL qes_write_atomic_positions(iun, obj%atomic_positions)
         !
      ENDIF
      !
      IF(obj%wyckoff_positions_ispresent) THEN
         CALL qes_write_wyckoff_positions(iun, obj%wyckoff_positions)
         !
      ENDIF
      !
      CALL qes_write_cell(iun, obj%cell)
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_atomic_structure

SUBROUTINE qes_init_atomic_structure(obj, tagname, nat, alat, alat_ispresent, &
                              bravais_index, bravais_index_ispresent, &
                              atomic_positions_ispresent, atomic_positions, &
                              wyckoff_positions_ispresent, wyckoff_positions, cell)
   IMPLICIT NONE

   TYPE(atomic_structure_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: nat
   LOGICAL  :: alat_ispresent
   REAL(DP), OPTIONAL :: alat
   LOGICAL  :: bravais_index_ispresent
   INTEGER , OPTIONAL :: bravais_index
   LOGICAL  :: atomic_positions_ispresent
   TYPE(atomic_positions_type) :: atomic_positions
   LOGICAL  :: wyckoff_positions_ispresent
   TYPE(wyckoff_positions_type) :: wyckoff_positions
   TYPE(cell_type) :: cell

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%nat = nat


   obj%alat_ispresent = alat_ispresent
   IF (obj%alat_ispresent) THEN
      obj%alat = alat
   ENDIF


   obj%bravais_index_ispresent = bravais_index_ispresent
   IF (obj%bravais_index_ispresent) THEN
      obj%bravais_index = bravais_index
   ENDIF

   obj%atomic_positions_ispresent = atomic_positions_ispresent
   IF(obj%atomic_positions_ispresent) THEN
      obj%atomic_positions = atomic_positions
   ENDIF
   obj%wyckoff_positions_ispresent = wyckoff_positions_ispresent
   IF(obj%wyckoff_positions_ispresent) THEN
      obj%wyckoff_positions = wyckoff_positions
   ENDIF
   obj%cell = cell

END SUBROUTINE qes_init_atomic_structure

SUBROUTINE qes_reset_atomic_structure(obj)
   IMPLICIT NONE
   TYPE(atomic_structure_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%atomic_positions_ispresent) THEN
      CALL qes_reset_atomic_positions(obj%atomic_positions)
      obj%atomic_positions_ispresent = .FALSE.
   ENDIF
   IF(obj%wyckoff_positions_ispresent) THEN
      CALL qes_reset_wyckoff_positions(obj%wyckoff_positions)
      obj%wyckoff_positions_ispresent = .FALSE.
   ENDIF
   CALL qes_reset_cell(obj%cell)

END SUBROUTINE qes_reset_atomic_structure


SUBROUTINE qes_write_step(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(step_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'n_step', obj%n_step)

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL qes_write_scf_conv(iun, obj%scf_conv)
      !
      CALL qes_write_atomic_structure(iun, obj%atomic_structure)
      !
      CALL qes_write_total_energy(iun, obj%total_energy)
      !
      CALL qes_write_matrix(iun, obj%forces)
      !
      IF(obj%stress_ispresent) THEN
         CALL qes_write_matrix(iun, obj%stress)
         !
      ENDIF
      !
      IF(obj%FCP_force_ispresent) THEN
         CALL iotk_write_begin(iun, 'FCP_force')
            WRITE(iun, '(E24.16)') obj%FCP_force
         CALL iotk_write_end(iun, 'FCP_force')
      ENDIF
      !
      IF(obj%FCP_tot_charge_ispresent) THEN
         CALL iotk_write_begin(iun, 'FCP_tot_charge')
            WRITE(iun, '(E24.16)') obj%FCP_tot_charge
         CALL iotk_write_end(iun, 'FCP_tot_charge')
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_step

SUBROUTINE qes_init_step(obj, tagname, n_step, scf_conv, atomic_structure, total_energy, &
                              forces, stress_ispresent, stress, FCP_force_ispresent, &
                              FCP_force, FCP_tot_charge_ispresent, FCP_tot_charge)
   IMPLICIT NONE

   TYPE(step_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
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

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%n_step = n_step

   obj%scf_conv = scf_conv
   obj%atomic_structure = atomic_structure
   obj%total_energy = total_energy
   obj%forces = forces
   obj%stress_ispresent = stress_ispresent
   IF(obj%stress_ispresent) THEN
      obj%stress = stress
   ENDIF
   obj%FCP_force_ispresent = FCP_force_ispresent
   IF(obj%FCP_force_ispresent) THEN
      obj%FCP_force = FCP_force
   ENDIF
   obj%FCP_tot_charge_ispresent = FCP_tot_charge_ispresent
   IF(obj%FCP_tot_charge_ispresent) THEN
      obj%FCP_tot_charge = FCP_tot_charge
   ENDIF

END SUBROUTINE qes_init_step

SUBROUTINE qes_reset_step(obj)
   IMPLICIT NONE
   TYPE(step_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   CALL qes_reset_scf_conv(obj%scf_conv)
   CALL qes_reset_atomic_structure(obj%atomic_structure)
   CALL qes_reset_total_energy(obj%total_energy)
   CALL qes_reset_matrix(obj%forces)
   IF(obj%stress_ispresent) THEN
      CALL qes_reset_matrix(obj%stress)
      obj%stress_ispresent = .FALSE.
   ENDIF
   IF(obj%FCP_force_ispresent) THEN
      obj%FCP_force_ispresent = .FALSE.
   ENDIF
   IF(obj%FCP_tot_charge_ispresent) THEN
      obj%FCP_tot_charge_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_step


SUBROUTINE qes_write_atomic_species(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(atomic_species_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'ntyp', obj%ntyp)

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      DO i = 1, obj%ndim_species
         CALL qes_write_species(iun, obj%species(i))
         !
      END DO
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_atomic_species

SUBROUTINE qes_init_atomic_species(obj, tagname, ntyp, ndim_species, species)
   IMPLICIT NONE

   TYPE(atomic_species_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: ntyp
   INTEGER  :: ndim_species
   TYPE(species_type ), DIMENSION( ndim_species )  :: species

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%ntyp = ntyp

   ALLOCATE(obj%species(SIZE(species)))
   DO i = 1, SIZE(species)
      obj%species(i) = species(i)
   ENDDO
   obj%ndim_species = ndim_species

END SUBROUTINE qes_init_atomic_species

SUBROUTINE qes_reset_atomic_species(obj)
   IMPLICIT NONE
   TYPE(atomic_species_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   DO i = 1, SIZE(obj%species)
      CALL qes_reset_species(obj%species(i))
   ENDDO
   IF (ALLOCATED(obj%species)) DEALLOCATE(obj%species)

END SUBROUTINE qes_reset_atomic_species


SUBROUTINE qes_write_output(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(output_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL qes_write_convergence_info(iun, obj%convergence_info)
      !
      CALL qes_write_algorithmic_info(iun, obj%algorithmic_info)
      !
      CALL qes_write_atomic_species(iun, obj%atomic_species)
      !
      CALL qes_write_atomic_structure(iun, obj%atomic_structure)
      !
      IF(obj%symmetries_ispresent) THEN
         CALL qes_write_symmetries(iun, obj%symmetries)
         !
      ENDIF
      !
      CALL qes_write_basis_set(iun, obj%basis_set)
      !
      CALL qes_write_dft(iun, obj%dft)
      !
      CALL qes_write_magnetization(iun, obj%magnetization)
      !
      CALL qes_write_total_energy(iun, obj%total_energy)
      !
      CALL qes_write_band_structure(iun, obj%band_structure)
      !
      IF(obj%forces_ispresent) THEN
         CALL qes_write_matrix(iun, obj%forces)
         !
      ENDIF
      !
      IF(obj%stress_ispresent) THEN
         CALL qes_write_matrix(iun, obj%stress)
         !
      ENDIF
      !
      IF(obj%electric_field_ispresent) THEN
         CALL qes_write_outputElectricField(iun, obj%electric_field)
         !
      ENDIF
      !
      IF(obj%FCP_force_ispresent) THEN
         CALL iotk_write_begin(iun, 'FCP_force')
            WRITE(iun, '(E24.16)') obj%FCP_force
         CALL iotk_write_end(iun, 'FCP_force')
      ENDIF
      !
      IF(obj%FCP_tot_charge_ispresent) THEN
         CALL iotk_write_begin(iun, 'FCP_tot_charge')
            WRITE(iun, '(E24.16)') obj%FCP_tot_charge
         CALL iotk_write_end(iun, 'FCP_tot_charge')
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_output

SUBROUTINE qes_init_output(obj, tagname, convergence_info, algorithmic_info, &
                              atomic_species, atomic_structure, symmetries_ispresent, &
                              symmetries, basis_set, dft, magnetization, total_energy, &
                              band_structure, forces_ispresent, forces, stress_ispresent, &
                              stress, electric_field_ispresent, electric_field, &
                              FCP_force_ispresent, FCP_force, FCP_tot_charge_ispresent, &
                              FCP_tot_charge)
   IMPLICIT NONE

   TYPE(output_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   TYPE(convergence_info_type) :: convergence_info
   TYPE(algorithmic_info_type) :: algorithmic_info
   TYPE(atomic_species_type) :: atomic_species
   TYPE(atomic_structure_type) :: atomic_structure
   LOGICAL  :: symmetries_ispresent
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

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%convergence_info = convergence_info
   obj%algorithmic_info = algorithmic_info
   obj%atomic_species = atomic_species
   obj%atomic_structure = atomic_structure
   obj%symmetries_ispresent = symmetries_ispresent
   IF(obj%symmetries_ispresent) THEN
      obj%symmetries = symmetries
   ENDIF
   obj%basis_set = basis_set
   obj%dft = dft
   obj%magnetization = magnetization
   obj%total_energy = total_energy
   obj%band_structure = band_structure
   obj%forces_ispresent = forces_ispresent
   IF(obj%forces_ispresent) THEN
      obj%forces = forces
   ENDIF
   obj%stress_ispresent = stress_ispresent
   IF(obj%stress_ispresent) THEN
      obj%stress = stress
   ENDIF
   obj%electric_field_ispresent = electric_field_ispresent
   IF(obj%electric_field_ispresent) THEN
      obj%electric_field = electric_field
   ENDIF
   obj%FCP_force_ispresent = FCP_force_ispresent
   IF(obj%FCP_force_ispresent) THEN
      obj%FCP_force = FCP_force
   ENDIF
   obj%FCP_tot_charge_ispresent = FCP_tot_charge_ispresent
   IF(obj%FCP_tot_charge_ispresent) THEN
      obj%FCP_tot_charge = FCP_tot_charge
   ENDIF

END SUBROUTINE qes_init_output

SUBROUTINE qes_reset_output(obj)
   IMPLICIT NONE
   TYPE(output_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   CALL qes_reset_convergence_info(obj%convergence_info)
   CALL qes_reset_algorithmic_info(obj%algorithmic_info)
   CALL qes_reset_atomic_species(obj%atomic_species)
   CALL qes_reset_atomic_structure(obj%atomic_structure)
   IF(obj%symmetries_ispresent) THEN
      CALL qes_reset_symmetries(obj%symmetries)
      obj%symmetries_ispresent = .FALSE.
   ENDIF
   CALL qes_reset_basis_set(obj%basis_set)

   CALL qes_reset_dft(obj%dft)

   CALL qes_reset_magnetization(obj%magnetization)
   CALL qes_reset_total_energy(obj%total_energy)
   CALL qes_reset_band_structure(obj%band_structure)

   IF(obj%forces_ispresent) THEN
      CALL qes_reset_matrix(obj%forces)
      obj%forces_ispresent = .FALSE.
   ENDIF
   IF(obj%stress_ispresent) THEN
      CALL qes_reset_matrix(obj%stress)
      obj%stress_ispresent = .FALSE.
   ENDIF
   IF(obj%electric_field_ispresent) THEN
      CALL qes_reset_outputElectricField(obj%electric_field)
      obj%electric_field_ispresent = .FALSE.
   ENDIF
   IF(obj%FCP_force_ispresent) THEN
      obj%FCP_force_ispresent = .FALSE.
   ENDIF
   IF(obj%FCP_tot_charge_ispresent) THEN
      obj%FCP_tot_charge_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_output


SUBROUTINE qes_write_lowhigh(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(lowhigh_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_lowhigh

SUBROUTINE qes_init_lowhigh(obj, tagname)
   IMPLICIT NONE

   TYPE(lowhigh_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

END SUBROUTINE qes_init_lowhigh

SUBROUTINE qes_reset_lowhigh(obj)
   IMPLICIT NONE
   TYPE(lowhigh_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_lowhigh


SUBROUTINE qes_write_controlRestartMode(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(controlRestartMode_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_controlRestartMode

SUBROUTINE qes_init_controlRestartMode(obj, tagname)
   IMPLICIT NONE

   TYPE(controlRestartMode_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

END SUBROUTINE qes_init_controlRestartMode

SUBROUTINE qes_reset_controlRestartMode(obj)
   IMPLICIT NONE
   TYPE(controlRestartMode_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_controlRestartMode


SUBROUTINE qes_write_calculation(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(calculation_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_calculation

SUBROUTINE qes_init_calculation(obj, tagname)
   IMPLICIT NONE

   TYPE(calculation_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

END SUBROUTINE qes_init_calculation

SUBROUTINE qes_reset_calculation(obj)
   IMPLICIT NONE
   TYPE(calculation_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_calculation


SUBROUTINE qes_write_control_variables(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(control_variables_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'title',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%title)
      CALL iotk_write_end(iun, 'title',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'calculation',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%calculation)
      CALL iotk_write_end(iun, 'calculation',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'restart_mode',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%restart_mode)
      CALL iotk_write_end(iun, 'restart_mode',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'prefix',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%prefix)
      CALL iotk_write_end(iun, 'prefix',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'pseudo_dir',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%pseudo_dir)
      CALL iotk_write_end(iun, 'pseudo_dir',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'outdir',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%outdir)
      CALL iotk_write_end(iun, 'outdir',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'stress',new_line=.FALSE.)
         IF (obj%stress) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'stress',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'forces',new_line=.FALSE.)
         IF (obj%forces) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'forces',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'wf_collect',new_line=.FALSE.)
         IF (obj%wf_collect) THEN
            WRITE(iun, '(A)',advance='no')  'true'
         ELSE
            WRITE(iun, '(A)',advance='no')  'false'
         ENDIF
      CALL iotk_write_end(iun, 'wf_collect',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'disk_io',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%disk_io)
      CALL iotk_write_end(iun, 'disk_io',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'max_seconds')
         WRITE(iun, '(I12)') obj%max_seconds
      CALL iotk_write_end(iun, 'max_seconds')
      IF(obj%nstep_ispresent) THEN
         CALL iotk_write_begin(iun, 'nstep')
            WRITE(iun, '(I12)') obj%nstep
         CALL iotk_write_end(iun, 'nstep')
      ENDIF
      !
      CALL iotk_write_begin(iun, 'etot_conv_thr')
         WRITE(iun, '(E24.16)') obj%etot_conv_thr
      CALL iotk_write_end(iun, 'etot_conv_thr')
      CALL iotk_write_begin(iun, 'forc_conv_thr')
         WRITE(iun, '(E24.16)') obj%forc_conv_thr
      CALL iotk_write_end(iun, 'forc_conv_thr')
      CALL iotk_write_begin(iun, 'press_conv_thr')
         WRITE(iun, '(E24.16)') obj%press_conv_thr
      CALL iotk_write_end(iun, 'press_conv_thr')
      CALL iotk_write_begin(iun, 'verbosity',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%verbosity)
      CALL iotk_write_end(iun, 'verbosity',indentation=.FALSE.)
      CALL iotk_write_begin(iun, 'print_every')
         WRITE(iun, '(I12)') obj%print_every
      CALL iotk_write_end(iun, 'print_every')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_control_variables

SUBROUTINE qes_init_control_variables(obj, tagname, title, calculation, restart_mode, &
                              prefix, pseudo_dir, outdir, stress, forces, wf_collect, &
                              disk_io, max_seconds, nstep_ispresent, nstep, etot_conv_thr, &
                              forc_conv_thr, press_conv_thr, verbosity, print_every)
   IMPLICIT NONE

   TYPE(control_variables_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: title
   CHARACTER(len=*) :: calculation
   CHARACTER(len=*) :: restart_mode
   CHARACTER(len=*) :: prefix
   CHARACTER(len=*) :: pseudo_dir
   CHARACTER(len=*) :: outdir
   LOGICAL  :: stress
   LOGICAL  :: forces
   LOGICAL  :: wf_collect
   CHARACTER(len=*) :: disk_io
   INTEGER  :: max_seconds
   LOGICAL  :: nstep_ispresent
   INTEGER  :: nstep
   REAL(DP) :: etot_conv_thr
   REAL(DP) :: forc_conv_thr
   REAL(DP) :: press_conv_thr
   CHARACTER(len=*) :: verbosity
   INTEGER  :: print_every

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
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
   obj%nstep_ispresent = nstep_ispresent
   IF(obj%nstep_ispresent) THEN
      obj%nstep = nstep
   ENDIF
   obj%etot_conv_thr = etot_conv_thr
   obj%forc_conv_thr = forc_conv_thr
   obj%press_conv_thr = press_conv_thr
   obj%verbosity = verbosity
   obj%print_every = print_every

END SUBROUTINE qes_init_control_variables

SUBROUTINE qes_reset_control_variables(obj)
   IMPLICIT NONE
   TYPE(control_variables_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   IF(obj%nstep_ispresent) THEN
      obj%nstep_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_control_variables


SUBROUTINE qes_write_input(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(input_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL qes_write_control_variables(iun, obj%control_variables)
      !
      CALL qes_write_atomic_species(iun, obj%atomic_species)
      !
      CALL qes_write_atomic_structure(iun, obj%atomic_structure)
      !
      CALL qes_write_dft(iun, obj%dft)
      !
      CALL qes_write_spin(iun, obj%spin)
      !
      CALL qes_write_bands(iun, obj%bands)
      !
      CALL qes_write_basis(iun, obj%basis)
      !
      CALL qes_write_electron_control(iun, obj%electron_control)
      !
      CALL qes_write_k_points_IBZ(iun, obj%k_points_IBZ)
      !
      CALL qes_write_ion_control(iun, obj%ion_control)
      !
      CALL qes_write_cell_control(iun, obj%cell_control)
      !
      IF(obj%symmetry_flags_ispresent) THEN
         CALL qes_write_symmetry_flags(iun, obj%symmetry_flags)
         !
      ENDIF
      !
      IF(obj%boundary_conditions_ispresent) THEN
         CALL qes_write_boundary_conditions(iun, obj%boundary_conditions)
         !
      ENDIF
      !
      IF(obj%ekin_functional_ispresent) THEN
         CALL qes_write_ekin_functional(iun, obj%ekin_functional)
         !
      ENDIF
      !
      IF(obj%external_atomic_forces_ispresent) THEN
         CALL qes_write_matrix(iun, obj%external_atomic_forces)
         !
      ENDIF
      !
      IF(obj%free_positions_ispresent) THEN
         CALL qes_write_integerMatrix(iun, obj%free_positions)
         !
      ENDIF
      !
      IF(obj%starting_atomic_velocities_ispresent) THEN
         CALL qes_write_matrix(iun, obj%starting_atomic_velocities)
         !
      ENDIF
      !
      IF(obj%electric_field_ispresent) THEN
         CALL qes_write_electric_field(iun, obj%electric_field)
         !
      ENDIF
      !
      IF(obj%atomic_constraints_ispresent) THEN
         CALL qes_write_atomic_constraints(iun, obj%atomic_constraints)
         !
      ENDIF
      !
      IF(obj%spin_constraints_ispresent) THEN
         CALL qes_write_spin_constraints(iun, obj%spin_constraints)
         !
      ENDIF
      !
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_input

SUBROUTINE qes_init_input(obj, tagname, control_variables, atomic_species, &
                              atomic_structure, dft, spin, bands, basis, electron_control, &
                              k_points_IBZ, ion_control, cell_control, &
                              symmetry_flags_ispresent, symmetry_flags, &
                              boundary_conditions_ispresent, boundary_conditions, &
                              ekin_functional_ispresent, ekin_functional, &
                              external_atomic_forces_ispresent, external_atomic_forces, &
                              free_positions_ispresent, free_positions, &
                              starting_atomic_velocities_ispresent, &
                              starting_atomic_velocities, electric_field_ispresent, &
                              electric_field, atomic_constraints_ispresent, &
                              atomic_constraints, spin_constraints_ispresent, &
                              spin_constraints)
   IMPLICIT NONE

   TYPE(input_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
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

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
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
   obj%symmetry_flags_ispresent = symmetry_flags_ispresent
   IF(obj%symmetry_flags_ispresent) THEN
      obj%symmetry_flags = symmetry_flags
   ENDIF
   obj%boundary_conditions_ispresent = boundary_conditions_ispresent
   IF(obj%boundary_conditions_ispresent) THEN
      obj%boundary_conditions = boundary_conditions
   ENDIF
   obj%ekin_functional_ispresent = ekin_functional_ispresent
   IF(obj%ekin_functional_ispresent) THEN
      obj%ekin_functional = ekin_functional
   ENDIF
   obj%external_atomic_forces_ispresent = external_atomic_forces_ispresent
   IF(obj%external_atomic_forces_ispresent) THEN
      obj%external_atomic_forces = external_atomic_forces
   ENDIF
   obj%free_positions_ispresent = free_positions_ispresent
   IF(obj%free_positions_ispresent) THEN
      obj%free_positions = free_positions
   ENDIF
   obj%starting_atomic_velocities_ispresent = starting_atomic_velocities_ispresent
   IF(obj%starting_atomic_velocities_ispresent) THEN
      obj%starting_atomic_velocities = starting_atomic_velocities
   ENDIF
   obj%electric_field_ispresent = electric_field_ispresent
   IF(obj%electric_field_ispresent) THEN
      obj%electric_field = electric_field
   ENDIF
   obj%atomic_constraints_ispresent = atomic_constraints_ispresent
   IF(obj%atomic_constraints_ispresent) THEN
      obj%atomic_constraints = atomic_constraints
   ENDIF
   obj%spin_constraints_ispresent = spin_constraints_ispresent
   IF(obj%spin_constraints_ispresent) THEN
      obj%spin_constraints = spin_constraints
   ENDIF

END SUBROUTINE qes_init_input

SUBROUTINE qes_reset_input(obj)
   IMPLICIT NONE
   TYPE(input_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

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
   IF(obj%symmetry_flags_ispresent) THEN
      CALL qes_reset_symmetry_flags(obj%symmetry_flags)
      obj%symmetry_flags_ispresent = .FALSE.
   ENDIF
   IF(obj%boundary_conditions_ispresent) THEN
      CALL qes_reset_boundary_conditions(obj%boundary_conditions)
      obj%boundary_conditions_ispresent = .FALSE.
   ENDIF
   IF(obj%ekin_functional_ispresent) THEN
      CALL qes_reset_ekin_functional(obj%ekin_functional)
      obj%ekin_functional_ispresent = .FALSE.
   ENDIF
   IF(obj%external_atomic_forces_ispresent) THEN
      CALL qes_reset_matrix(obj%external_atomic_forces)
      obj%external_atomic_forces_ispresent = .FALSE.
   ENDIF
   IF(obj%free_positions_ispresent) THEN
      CALL qes_reset_integerMatrix(obj%free_positions)
      obj%free_positions_ispresent = .FALSE.
   ENDIF
   IF(obj%starting_atomic_velocities_ispresent) THEN
      CALL qes_reset_matrix(obj%starting_atomic_velocities)
      obj%starting_atomic_velocities_ispresent = .FALSE.
   ENDIF
   IF(obj%electric_field_ispresent) THEN
      CALL qes_reset_electric_field(obj%electric_field)
      obj%electric_field_ispresent = .FALSE.
   ENDIF
   IF(obj%atomic_constraints_ispresent) THEN
      CALL qes_reset_atomic_constraints(obj%atomic_constraints)
      obj%atomic_constraints_ispresent = .FALSE.
   ENDIF
   IF(obj%spin_constraints_ispresent) THEN
      CALL qes_reset_spin_constraints(obj%spin_constraints)
      obj%spin_constraints_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_input


SUBROUTINE qes_write_parallel_info(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(parallel_info_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL iotk_write_begin(iun, 'nprocs')
         WRITE(iun, '(I12)') obj%nprocs
      CALL iotk_write_end(iun, 'nprocs')
      CALL iotk_write_begin(iun, 'nthreads')
         WRITE(iun, '(I12)') obj%nthreads
      CALL iotk_write_end(iun, 'nthreads')
      CALL iotk_write_begin(iun, 'ntasks')
         WRITE(iun, '(I12)') obj%ntasks
      CALL iotk_write_end(iun, 'ntasks')
      CALL iotk_write_begin(iun, 'nbgrp')
         WRITE(iun, '(I12)') obj%nbgrp
      CALL iotk_write_end(iun, 'nbgrp')
      CALL iotk_write_begin(iun, 'npool')
         WRITE(iun, '(I12)') obj%npool
      CALL iotk_write_end(iun, 'npool')
      CALL iotk_write_begin(iun, 'ndiag')
         WRITE(iun, '(I12)') obj%ndiag
      CALL iotk_write_end(iun, 'ndiag')
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_parallel_info

SUBROUTINE qes_init_parallel_info(obj, tagname, nprocs, nthreads, ntasks, nbgrp, npool, &
                              ndiag)
   IMPLICIT NONE

   TYPE(parallel_info_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: nprocs
   INTEGER  :: nthreads
   INTEGER  :: ntasks
   INTEGER  :: nbgrp
   INTEGER  :: npool
   INTEGER  :: ndiag

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%nprocs = nprocs
   obj%nthreads = nthreads
   obj%ntasks = ntasks
   obj%nbgrp = nbgrp
   obj%npool = npool
   obj%ndiag = ndiag

END SUBROUTINE qes_init_parallel_info

SUBROUTINE qes_reset_parallel_info(obj)
   IMPLICIT NONE
   TYPE(parallel_info_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_parallel_info


SUBROUTINE qes_write_created(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(created_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'DATE', TRIM(obj%DATE))
   CALL iotk_write_attr(attr, 'TIME', TRIM(obj%TIME))

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr),new_line=.FALSE.)
      !
      WRITE(iun, '(A)',advance='no')  TRIM(obj%created)
   CALL iotk_write_end(iun, TRIM(obj%tagname),indentation=.FALSE.)
   !
END SUBROUTINE qes_write_created

SUBROUTINE qes_init_created(obj, tagname, DATE, TIME, created)
   IMPLICIT NONE

   TYPE(created_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: DATE
   CHARACTER(len=*) :: TIME
   CHARACTER(len=*) :: created

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%DATE = TRIM(DATE)


   obj%TIME = TRIM(TIME)

   obj%created = created

END SUBROUTINE qes_init_created

SUBROUTINE qes_reset_created(obj)
   IMPLICIT NONE
   TYPE(created_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_created


SUBROUTINE qes_write_creator(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(creator_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'NAME', TRIM(obj%NAME))
   CALL iotk_write_attr(attr, 'VERSION', TRIM(obj%VERSION))

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr),new_line=.FALSE.)
      !
      WRITE(iun, '(A)',advance='no')  TRIM(obj%creator)
   CALL iotk_write_end(iun, TRIM(obj%tagname),indentation=.FALSE.)
   !
END SUBROUTINE qes_write_creator

SUBROUTINE qes_init_creator(obj, tagname, NAME, VERSION, creator)
   IMPLICIT NONE

   TYPE(creator_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: NAME
   CHARACTER(len=*) :: VERSION
   CHARACTER(len=*) :: creator

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%NAME = TRIM(NAME)


   obj%VERSION = TRIM(VERSION)

   obj%creator = creator

END SUBROUTINE qes_init_creator

SUBROUTINE qes_reset_creator(obj)
   IMPLICIT NONE
   TYPE(creator_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_creator


SUBROUTINE qes_write_xml_format(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(xml_format_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "
   CALL iotk_write_attr(attr, 'NAME', TRIM(obj%NAME))
   CALL iotk_write_attr(attr, 'VERSION', TRIM(obj%VERSION))

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr),new_line=.FALSE.)
      !
      WRITE(iun, '(A)',advance='no')  TRIM(obj%xml_format)
   CALL iotk_write_end(iun, TRIM(obj%tagname),indentation=.FALSE.)
   !
END SUBROUTINE qes_write_xml_format

SUBROUTINE qes_init_xml_format(obj, tagname, NAME, VERSION, xml_format)
   IMPLICIT NONE

   TYPE(xml_format_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: NAME
   CHARACTER(len=*) :: VERSION
   CHARACTER(len=*) :: xml_format

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%NAME = TRIM(NAME)


   obj%VERSION = TRIM(VERSION)

   obj%xml_format = xml_format

END SUBROUTINE qes_init_xml_format

SUBROUTINE qes_reset_xml_format(obj)
   IMPLICIT NONE
   TYPE(xml_format_type) :: obj
   INTEGER  :: i

   obj%tagname = ""


END SUBROUTINE qes_reset_xml_format


SUBROUTINE qes_write_general_info(iun, obj)
   IMPLICIT NONE

   INTEGER  :: iun
   TYPE(general_info_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   attr = " "

   CALL iotk_write_begin(iun, TRIM(obj%tagname), attr=TRIM(attr))
      !
      CALL qes_write_xml_format(iun, obj%xml_format)
      !
      CALL qes_write_creator(iun, obj%creator)
      !
      CALL qes_write_created(iun, obj%created)
      !
      CALL iotk_write_begin(iun, 'job',new_line=.FALSE.)
         WRITE(iun, '(A)',advance='no')  TRIM(obj%job)
      CALL iotk_write_end(iun, 'job',indentation=.FALSE.)
   CALL iotk_write_end(iun, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_general_info

SUBROUTINE qes_init_general_info(obj, tagname, xml_format, creator, created, job)
   IMPLICIT NONE

   TYPE(general_info_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   TYPE(xml_format_type) :: xml_format
   TYPE(creator_type) :: creator
   TYPE(created_type) :: created
   CHARACTER(len=*) :: job

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%xml_format = xml_format
   obj%creator = creator
   obj%created = created
   obj%job = job

END SUBROUTINE qes_init_general_info

SUBROUTINE qes_reset_general_info(obj)
   IMPLICIT NONE
   TYPE(general_info_type) :: obj
   INTEGER  :: i

   obj%tagname = ""

   CALL qes_reset_xml_format(obj%xml_format)
   CALL qes_reset_creator(obj%creator)
   CALL qes_reset_created(obj%created)

END SUBROUTINE qes_reset_general_info


END MODULE qes_libs_module
