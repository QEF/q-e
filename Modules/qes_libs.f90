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
   USE kinds,           ONLY: DP
   USE qes_types_module
   USE FoX_wxml,        ONLY: xmlf_t, xml_NewElement, xml_addAttribute, xml_addCharacters, xml_EndElement,&
                              xml_addnewline
   IMPLICIT NONE
   !
   INTEGER, PARAMETER      :: max_real_per_line=5
   CHARACTER(32)           :: fmtstr
   !
   PUBLIC
   PRIVATE :: fmtstr, DP
   PRIVATE :: xmlf_t, xml_NewElement, xml_addAttribute, xml_addCharacters, xml_EndElement, xml_addnewline
!
INTERFACE qes_init_matrix
   MODULE PROCEDURE  qes_init_matrix_1, qes_init_matrix_2, qes_init_matrix_3
END INTERFACE

INTERFACE qes_init_integerMatrix
  MODULE PROCEDURE  qes_init_integerMatrix_1, qes_init_integerMatrix_2, qes_init_integerMatrix_3
END INTERFACE
 
CONTAINS
!

SUBROUTINE qes_write_xml_format(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)          :: xp
   TYPE(xml_format_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
  
   CALL xml_NewElement (xp, TRIM(obj%tagname))
   CALL xml_addAttribute(xp, 'NAME', TRIM(obj%NAME))
   CALL xml_addAttribute(xp, 'VERSION', TRIM(obj%VERSION))

      !
      CALL xml_AddCharacters(xp, TRIM(obj%xml_format) )
      
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite = .FALSE.

END SUBROUTINE qes_reset_xml_format

SUBROUTINE qes_write_creator(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)       :: xp
   TYPE(creator_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   
   CALL xml_NewElement (xp, TRIM(obj%tagname) )
   CALL xml_addAttribute(xp, 'NAME', TRIM(obj%NAME))
   CALL xml_addAttribute(xp, 'VERSION', TRIM(obj%VERSION))

      !
      CALL xml_addCharacters(xp, TRIM(obj%creator) )
   CALL xml_endElement(xp, TRIM(obj%tagname) )
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
   obj%lwrite  = .FALSE.

END SUBROUTINE qes_reset_creator

SUBROUTINE qes_write_created(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)       :: xp
   TYPE(created_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
    
   CALL xml_NewElement(xp, TRIM(obj%tagname) )
 
   CALL xml_addAttribute( xp, 'DATE', TRIM(obj%DATE))
   CALL xml_addAttribute( xp, 'TIME', TRIM(obj%TIME))

      !
      CALL xml_addCharacters(xp, TRIM(obj%created) )
   CALL xml_EndElement(xp, TRIM(obj%tagname) )
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
   obj%lwrite = .FALSE.

END SUBROUTINE qes_reset_created

SUBROUTINE qes_write_atom(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(atom_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   
   CALL xml_NewElement(xp, TRIM(obj%tagname))
   CALL xml_addAttribute( xp, 'name', TRIM(obj%name))
   IF(obj%position_ispresent) THEN
      CALL xml_addAttribute( xp, 'position', TRIM(obj%position))
   END IF
   IF(obj%index_ispresent) THEN
      CALL xml_addAttribute( xp, 'index', obj%index)
   END IF

      !
      CALL xml_addCharacters(xp, obj%atom, fmt = 's16')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_atom

SUBROUTINE qes_init_atom(obj, tagname, name, position, position_ispresent, index, index_ispresent, &
                              atom)
   IMPLICIT NONE

   TYPE(atom_type) :: obj
   CHARACTER(len=*) :: tagname
   
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
   obj%lwrite  = .FALSE.

END SUBROUTINE qes_reset_atom


SUBROUTINE qes_write_qpoint_grid(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(qpoint_grid_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN

   CALL xml_NewElement  (xp, TRIM(obj%tagname) )
   CALL xml_addAttribute( xp, 'nqx1', obj%nqx1)
   CALL xml_addAttribute( xp, 'nqx2', obj%nqx2)
   CALL xml_addAttribute( xp, 'nqx3', obj%nqx3)

      !
      CALL xml_addCharacters(xp, TRIM(obj%qpoint_grid) )
   CALL xml_EndElement(xp, TRIM(obj%tagname) )
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
   obj%lwrite  = .FALSE.

END SUBROUTINE qes_reset_qpoint_grid

SUBROUTINE qes_write_HubbardCommon(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(HubbardCommon_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN

   CALL xml_NewElement(xp, TRIM(obj%tagname) )
   CALL xml_addAttribute( xp, 'specie', TRIM(obj%specie))
   CALL xml_addAttribute( xp, 'label', TRIM(obj%label))

      !
      CALL xml_addCharacters(xp, obj%HubbardCommon, fmt = 's16')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

END SUBROUTINE qes_reset_HubbardCommon

SUBROUTINE qes_write_HubbardJ(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(HubbardJ_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN


   CALL xml_NewElement(xp, TRIM(obj%tagname) )

   CALL xml_addAttribute( xp, 'specie', TRIM(obj%specie))
   CALL xml_addAttribute( xp, 'label', TRIM(obj%label))

      !
      CALL xml_addCharacters(xp, obj%HubbardJ, fmt = 's16')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

END SUBROUTINE qes_reset_HubbardJ

SUBROUTINE qes_write_starting_ns(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(starting_ns_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   
      
   CALL xml_NewElement(xp, TRIM(obj%tagname) )

   CALL xml_addAttribute( xp, 'specie', TRIM(obj%specie))
   CALL xml_addAttribute( xp, 'label', TRIM(obj%label))
   CALL xml_addAttribute( xp, 'spin', obj%spin)
   CALL xml_addAttribute( xp, 'size', obj%size)

      !
        CALL xml_addCharacters(xp, obj%starting_ns, fmt = 's16')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_starting_ns

SUBROUTINE qes_init_starting_ns(obj, tagname, specie, label, spin, vec)
   IMPLICIT NONE

   TYPE(starting_ns_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: specie
   CHARACTER(len=*) :: label
   INTEGER  :: spin
   REAL(DP), DIMENSION(:) :: vec

   INTEGER  :: ndim_vec
   ndim_vec = size(vec)
   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%specie = TRIM(specie)


   obj%label = TRIM(label)


   obj%spin = spin
   obj%size = ndim_vec

   ALLOCATE(obj%starting_ns(ndim_vec))
   obj%starting_ns(:) = vec(:)

END SUBROUTINE qes_init_starting_ns

SUBROUTINE qes_reset_starting_ns(obj)
   IMPLICIT NONE
   TYPE(starting_ns_type) :: obj
   INTEGER  :: i

   obj%tagname = ""
   obj%lwrite  = .FALSE.
   IF (ALLOCATED(obj%starting_ns))  DEALLOCATE(obj%starting_ns)

END SUBROUTINE qes_reset_starting_ns


SUBROUTINE qes_write_Hubbard_ns(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(Hubbard_ns_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   CALL xml_NewElement(xp, TRIM(obj%tagname) )
   CALL xml_addAttribute( xp, 'specie', TRIM(obj%specie))
   CALL xml_addAttribute( xp, 'label', TRIM(obj%label))
   CALL xml_addAttribute( xp, 'spin', obj%spin)
   CALL xml_addAttribute( xp, 'index', obj%index)
   CALL xml_addAttribute( xp, 'rank', obj%rank) 
   CALL xml_addAttribute( xp, 'dims', obj%dims)
   CALL xml_addAttribute( xp, 'order',TRIM(obj%order))

      !
   CALL xml_addNewLine (xp)
   DO i=1, obj%dims(2)
      CALL xml_addCharacters(xp, obj%Hubbard_ns((i-1)*obj%dims(1)+1:i*obj%dims(1) ),fmt ='s16')
      CALL xml_addNewLine(xp)
   END DO 
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_Hubbard_ns

SUBROUTINE qes_init_Hubbard_ns(obj, tagname, specie, label, spin, index, mat)
   IMPLICIT NONE

   TYPE(Hubbard_ns_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: specie
   CHARACTER(len=*) :: label
   INTEGER  :: spin
   INTEGER  :: index
   REAL(DP), DIMENSION(:,:) :: mat

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%specie = TRIM(specie)


   obj%label = TRIM(label)


   obj%spin = spin


   obj%index = index
   obj%rank = 2
     
   ALLOCATE(obj%Hubbard_ns(size(mat,1)*size(mat,2)), obj%dims(2) )
   obj%dims = shape(mat)
   obj%order = "F"
   DO i = 1, obj%dims(2)
     obj%Hubbard_ns((i-1)*obj%dims(1)+1:i*obj%dims(1) ) = mat(:,i)
   END DO

END SUBROUTINE qes_init_Hubbard_ns

SUBROUTINE qes_reset_Hubbard_ns(obj)
   IMPLICIT NONE
   TYPE(Hubbard_ns_type) :: obj
   INTEGER  :: i

   obj%tagname = ""
   obj%lwrite = .FALSE.
   IF (ALLOCATED(obj%Hubbard_ns))  DEALLOCATE(obj%Hubbard_ns)

END SUBROUTINE qes_reset_Hubbard_ns

SUBROUTINE qes_write_smearing(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(smearing_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   
   CALL xml_NewElement(xp, TRIM(obj%tagname))

   CALL xml_addAttribute( xp, 'degauss', obj%degauss)

      !
      CALL xml_addCharacters(xp, TRIM(obj%smearing) )
   CALL xml_EndElement(xp, TRIM(obj%tagname) )
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
   obj%lwrite = .FALSE.

END SUBROUTINE qes_reset_smearing

SUBROUTINE qes_write_occupations(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(occupations_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
    
   CALL xml_NewElement(xp, TRIM(obj%tagname))
   IF(obj%spin_ispresent) THEN
      CALL xml_addAttribute( xp, 'spin', obj%spin)
   END IF

      !
      CALL xml_addCharacters(xp, TRIM(obj%occupations) )
   CALL xml_EndElement(xp, TRIM(obj%tagname) )
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
   obj%lwrite = .FALSE.

END SUBROUTINE qes_reset_occupations


SUBROUTINE qes_write_basisSetItem(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(basisSetItem_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN


   CALL xml_NewElement(xp, TRIM(obj%tagname) )

   CALL xml_addAttribute( xp, 'nr1', obj%nr1)
   CALL xml_addAttribute( xp, 'nr2', obj%nr2)
   CALL xml_addAttribute( xp, 'nr3', obj%nr3)

      !
      IF (TRIM(obj%basisSetItem) /= "")  CALL xml_addCharacters(xp, TRIM(obj%basisSetItem) )
   CALL xml_EndElement(xp, TRIM(obj%tagname) )
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
   obj%lwrite = .FALSE.

END SUBROUTINE qes_reset_basisSetItem


SUBROUTINE qes_write_monkhorst_pack(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(monkhorst_pack_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   
   CALL xml_NewElement(xp, TRIM(obj%tagname) )

   CALL xml_addAttribute( xp, 'nk1', obj%nk1)
   CALL xml_addAttribute( xp, 'nk2', obj%nk2)
   CALL xml_addAttribute( xp, 'nk3', obj%nk3)
   CALL xml_addAttribute( xp, 'k1', obj%k1)
   CALL xml_addAttribute( xp, 'k2', obj%k2)
   CALL xml_addAttribute( xp, 'k3', obj%k3)

      !
      IF (obj%monkhorst_pack /="")   CALL xml_addCharacters(xp, TRIM(obj%monkhorst_pack) )
   CALL xml_EndElement(xp, TRIM(obj%tagname) )
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
   obj%lwrite = .FALSE.


END SUBROUTINE qes_reset_monkhorst_pack


SUBROUTINE qes_write_k_point(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(k_point_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
   IF(obj%weight_ispresent) THEN
      CALL xml_addAttribute( xp, 'weight', obj%weight)
   END IF
   IF(obj%label_ispresent) THEN
      CALL xml_addAttribute( xp, 'label', TRIM(obj%label))
   END IF

      !
      CALL xml_addCharacters(xp, obj%k_point, fmt = 's16')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

END SUBROUTINE qes_reset_k_point



SUBROUTINE qes_write_inputOccupations(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(inputOccupations_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN


   CALL xml_NewElement(xp, TRIM(obj%tagname) )

   CALL xml_addAttribute( xp, 'ispin', obj%ispin)
   CALL xml_addAttribute( xp, 'spin_factor', obj%spin_factor)

      !
         CALL xml_addCharacters(xp, obj%inputOccupations, fmt = 's16')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_inputOccupations

SUBROUTINE qes_init_inputOccupations(obj, tagname, ispin, spin_factor, vec)
   IMPLICIT NONE

   TYPE(inputOccupations_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: ispin
   REAL(DP) :: spin_factor
   INTEGER  :: ndim_vec
   REAL(DP), DIMENSION(:) :: vec

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%ispin = ispin


   obj%spin_factor = spin_factor

   ALLOCATE(obj%inputOccupations(size(vec)))
   obj%inputOccupations(:) = vec(:)

END SUBROUTINE qes_init_inputOccupations

SUBROUTINE qes_reset_inputOccupations(obj)
   IMPLICIT NONE
   TYPE(inputOccupations_type) :: obj
   INTEGER  :: i

   obj%tagname = ""
   obj%lwrite  = .FALSE.
   IF (ALLOCATED(obj%inputOccupations))  DEALLOCATE(obj%inputOccupations)

END SUBROUTINE qes_reset_inputOccupations


SUBROUTINE qes_write_phase(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(phase_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN

   CALL xml_NewElement(xp, TRIM(obj%tagname) )

   IF(obj%ionic_ispresent) THEN
      CALL xml_addAttribute( xp, 'ionic', obj%ionic)
   END IF
   IF(obj%electronic_ispresent) THEN
      CALL xml_addAttribute( xp, 'electronic', obj%electronic)
   END IF
   IF(obj%modulus_ispresent) THEN
      CALL xml_addAttribute( xp, 'modulus', TRIM(obj%modulus))
   END IF

      !
      CALL xml_addCharacters(xp, obj%phase, fmt = 's16')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.


END SUBROUTINE qes_reset_phase


SUBROUTINE qes_write_equivalent_atoms(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(equivalent_atoms_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN


   CALL xml_NewElement(xp, TRIM(obj%tagname) )

   CALL xml_addAttribute( xp, 'nat', obj%nat)
   ! 
   CALL xml_addAttribute( xp, 'size', SIZE(obj%equivalent_atoms))

      !
         CALL xml_addCharacters(xp, obj%equivalent_atoms )
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_equivalent_atoms

SUBROUTINE qes_init_equivalent_atoms(obj, tagname, nat, index_list)
   IMPLICIT NONE

   TYPE(equivalent_atoms_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   INTEGER  :: nat
   INTEGER, DIMENSION(:) :: index_list

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.

   obj%nat = nat

   ALLOCATE(obj%equivalent_atoms(size(index_list)))
   obj%equivalent_atoms(:) = index_list(:)

END SUBROUTINE qes_init_equivalent_atoms

SUBROUTINE qes_reset_equivalent_atoms(obj)
   IMPLICIT NONE
   TYPE(equivalent_atoms_type) :: obj
   INTEGER  :: i

   obj%tagname = ""
   obj%lwrite = .FALSE.
   IF (ALLOCATED(obj%equivalent_atoms))  DEALLOCATE(obj%equivalent_atoms)

END SUBROUTINE qes_reset_equivalent_atoms


SUBROUTINE qes_write_info(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(info_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN


   CALL xml_NewElement(xp, TRIM(obj%tagname) )

   IF(obj%name_ispresent) THEN
      CALL xml_addAttribute( xp, 'name', TRIM(obj%name))
   END IF
   IF(obj%class_ispresent) THEN
      CALL xml_addAttribute( xp, 'class', TRIM(obj%class))
   END IF
   IF(obj%time_reversal_ispresent) CALL xml_addAttribute(xp, 'time_reversal', obj%time_reversal)
      !
      CALL xml_addCharacters(xp, TRIM(obj%info) )
   CALL xml_EndElement(xp, TRIM(obj%tagname) )
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
   obj%lwrite  = .FALSE.

END SUBROUTINE qes_reset_info


SUBROUTINE qes_write_closed(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(closed_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   
     
   CALL xml_NewElement(xp, TRIM(obj%tagname) )

   CALL xml_addAttribute( xp, 'DATE', TRIM(obj%DATE))
   CALL xml_addAttribute( xp, 'TIME', TRIM(obj%TIME))

      !
      CALL xml_addCharacters(xp, TRIM(obj%closed) )
   CALL xml_EndElement(xp, TRIM(obj%tagname) )
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
   obj%lwrite  = .FALSE.


END SUBROUTINE qes_reset_closed


SUBROUTINE qes_write_vector( xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)       :: xp
   TYPE(vector_type)  :: obj
   ! 
   INTEGER            :: i 

   IF ( .NOT. obj%lwrite ) RETURN 

   CALL xml_NewElement( xp, TRIM(obj%tagname))
   CALL xml_addAttribute( xp, 'size', obj%size)
   CALL xml_addCharacters( xp, obj%vector, fmt = 's16')
   CALL xml_EndElement( xp, TRIM(obj%tagname)) 

END SUBROUTINE qes_write_vector

SUBROUTINE qes_init_vector( obj, tagname, vec)
   IMPLICIT NONE
   
   TYPE (vector_type)     :: obj
   CHARACTER(LEN=*)       :: tagname
   REAL(DP),DIMENSION(:)  :: vec
  
   obj%tagname = TRIM(tagname)
   obj%lwrite  = .TRUE.
   obj%lread   = .TRUE.
   obj%size =  size(vec)
   ALLOCATE ( obj%vector(obj%size)) 
   obj%vector = vec
END SUBROUTINE qes_init_vector

SUBROUTINE qes_reset_vector( obj )
   IMPLICIT NONE
   
   TYPE (vector_type)   :: obj

   obj%tagname = ""
   obj%lwrite  = .FALSE.
  
   IF (ALLOCATED(obj%vector)) DEALLOCATE (obj%vector)
END SUBROUTINE qes_reset_vector  


SUBROUTINE qes_write_integerVector( xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)       :: xp
   TYPE(integerVector_type)  :: obj
   ! 
   INTEGER            :: i 

   IF ( .NOT. obj%lwrite ) RETURN 

   CALL xml_NewElement( xp, TRIM(obj%tagname))
   CALL xml_addAttribute( xp, 'size', obj%size)
   CALL xml_addCharacters( xp, obj%integerVector ) 
   CALL xml_EndElement( xp, TRIM(obj%tagname)) 

END SUBROUTINE qes_write_integerVector

SUBROUTINE qes_init_integerVector( obj, tagname, vec)
   IMPLICIT NONE
   
   TYPE (integerVector_type)     :: obj
   CHARACTER(LEN=*)       :: tagname
   INTEGER, DIMENSION(:)  :: vec
  
   obj%tagname = TRIM(tagname)
   obj%lwrite  = .TRUE.
   obj%lread   = .TRUE.
   obj%size =  size(vec)
   ALLOCATE ( obj%integerVector(obj%size)) 
   obj%integerVector = vec
END SUBROUTINE qes_init_integerVector

SUBROUTINE qes_reset_integerVector( obj )
   IMPLICIT NONE

   TYPE (integerVector_type)   :: obj 

   obj%tagname = ""
   obj%lwrite  = .FALSE.
  
   IF (ALLOCATED(obj%integerVector)) DEALLOCATE (obj%integerVector)
END SUBROUTINE qes_reset_integerVector  


SUBROUTINE qes_write_matrix(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(matrix_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   

   CALL xml_NewElement(xp, TRIM(obj%tagname) )
   CALL xml_addAttribute (xp,'rank', obj%rank)
   CALL xml_addAttribute (xp,'dims', obj%dims)
   CALL xml_addAttribute( xp, 'order', TRIM(obj%order) )    
      IF (obj%rank == 2 ) THEN
          CALL xml_addNewLine(xp)
          DO i = 1, obj%dims(2)
             CALL xml_addCharacters(xp,obj%matrix((i-1)*obj%dims(1)+1 : i*obj%dims(1))  , fmt = 's16')
             CALL xml_addNewLine(xp)
          END DO
      ELSE
          !
          CALL xml_addCharacters(xp, obj%matrix, fmt = 's16')
          !
      END IF
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_matrix

SUBROUTINE qes_init_matrix_1(obj, tagname, dims, mat, order)
   IMPLICIT NONE

   TYPE(matrix_type),        INTENT(OUT)   :: obj
   CHARACTER(len=*),         INTENT(IN)    :: tagname
   INTEGER,DIMENSION(:),     INTENT(IN)    :: dims
   REAL(DP), DIMENSION(:),   INTENT(IN)    :: mat
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN)    :: order 
   INTEGER  :: i, rank, length

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   rank = size(dims) 
   length = 1
   DO i =1, rank
      length=length*dims(i)
   END DO
   ALLOCATE(obj%matrix(length), obj%dims(rank) )
   obj%matrix(1:length) = mat(1:length)
   obj%dims = dims
   obj%rank = rank
   IF (PRESENT (order) ) THEN 
      obj%order = TRIM(order)
   ELSE 
      obj%order = 'F'
   END IF
END SUBROUTINE qes_init_matrix_1

SUBROUTINE qes_init_matrix_2(obj, tagname, dims, mat)
   IMPLICIT NONE

   TYPE(matrix_type),        INTENT(OUT)   :: obj
   CHARACTER(len=*),         INTENT(IN)    :: tagname
   INTEGER,DIMENSION(:),     INTENT(IN)    :: dims
   REAL(DP), DIMENSION(:,:),   INTENT(IN)    :: mat
   INTEGER  :: i, rank, length

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   rank = size(dims) 
   length = 1
   DO i =1, rank
      length=length*dims(i)
   END DO
   ALLOCATE(obj%matrix(length), obj%dims(rank) )
   obj%matrix(1:length) = reshape(mat,[length]) 
   obj%dims = dims
   obj%rank = rank
   !
   obj%order = 'F'
END SUBROUTINE qes_init_matrix_2

SUBROUTINE qes_init_matrix_3(obj, tagname, dims, mat )
   IMPLICIT NONE

   TYPE(matrix_type),        INTENT(OUT)   :: obj
   CHARACTER(len=*),         INTENT(IN)    :: tagname
   INTEGER,DIMENSION(:),     INTENT(IN)    :: dims
   REAL(DP), DIMENSION(:,:,:),   INTENT(IN)    :: mat
   INTEGER  :: i, rank, length

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   rank = size(dims) 
   length = 1
   DO i =1, rank
      length=length*dims(i)
   END DO
   ALLOCATE(obj%matrix(length), obj%dims(rank) )
   obj%matrix(1:length) = reshape(mat,[length]) 
   obj%dims = dims
   obj%rank = rank
   !
   obj%order = 'F'
END SUBROUTINE qes_init_matrix_3



SUBROUTINE qes_reset_matrix(obj)
   IMPLICIT NONE
   TYPE(matrix_type) :: obj
   INTEGER  :: i

   obj%tagname = ""
   obj%lwrite  = .FALSE.
   IF (ALLOCATED(obj%matrix))  DEALLOCATE(obj%matrix)
   IF (ALLOCATED(obj%dims)) DEALLOCATE (obj%dims)
END SUBROUTINE qes_reset_matrix


SUBROUTINE qes_write_integerMatrix(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(integerMatrix_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   

   CALL xml_NewElement(xp, TRIM(obj%tagname) )
      !
      CALL xml_addAttribute ( xp, 'rank', obj%rank) 
      CALL xml_addAttribute ( xp, 'dims', obj%dims) 
      CALL xml_addAttribute ( xp, 'order', TRIM(obj%order) )
      IF (obj%rank == 2) THEN
         CALL xml_addNewLine (xp)
         DO i =1, obj%dims(2)
            CALL xml_addCharacters(xp, obj%integerMatrix( (i-1)*obj%dims(1)+1: i*obj%dims(1) ) )
            CALL xml_addNewLine(xp)
         END DO
      ELSE
         CALL xml_addCharacters(xp, obj%integerMatrix )
      END IF 
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_integerMatrix

SUBROUTINE qes_init_integerMatrix_1(obj, tagname, dims, int_mat, order )
   IMPLICIT NONE

   TYPE(integerMatrix_type),INTENT(OUT):: obj
   CHARACTER(len=*),        INTENT(IN) :: tagname
   INTEGER,DIMENSION(:),    INTENT(IN) :: dims
   INTEGER,DIMENSION(:),    INTENT(IN) :: int_mat
   CHARACTER(LEN=*),OPTIONAL,INTENT(IN):: order
   INTEGER  :: i,rank, length
  

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   rank = size(dims)
   length =1 
   DO i =1, rank
      length = length*dims(i)
   END DO
   ALLOCATE(obj%integerMatrix(length), obj%dims(rank) )
   obj%integerMatrix(1:length) = int_mat(1:length)
   obj%rank           = rank
   obj%dims           = dims
   IF (PRESENT(order) ) THEN 
      obj%order = TRIM(order) 
   ELSE 
      obj%order = 'F'
   END IF

END SUBROUTINE qes_init_integerMatrix_1

SUBROUTINE qes_init_integerMatrix_2(obj, tagname, dims, int_mat)
   IMPLICIT NONE

   TYPE(integerMatrix_type),INTENT(OUT):: obj
   CHARACTER(len=*),        INTENT(IN) :: tagname
   INTEGER,DIMENSION(:),    INTENT(IN) :: dims
   INTEGER,DIMENSION(:,:),    INTENT(IN) :: int_mat
   INTEGER  :: i,rank, length
  

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   rank = size(dims)
   length =1 
   DO i =1, rank
      length = length*dims(i)
   END DO
   ALLOCATE(obj%integerMatrix(length), obj%dims(rank) )
   obj%integerMatrix(1:length) = reshape(int_mat, [length])
   obj%rank           = rank
   obj%dims           = dims
   !
   obj%order = 'F'
   !
END SUBROUTINE qes_init_integerMatrix_2

SUBROUTINE qes_init_integerMatrix_3(obj, tagname, dims, int_mat )
   IMPLICIT NONE

   TYPE(integerMatrix_type),INTENT(OUT):: obj
   CHARACTER(len=*),        INTENT(IN) :: tagname
   INTEGER,DIMENSION(:),    INTENT(IN) :: dims
   INTEGER,DIMENSION(:,:,:),    INTENT(IN) :: int_mat
   INTEGER  :: i,rank, length
  

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   rank = size(dims)
   length =1 
   DO i =1, rank
      length = length*dims(i)
   END DO
   ALLOCATE(obj%integerMatrix(length), obj%dims(rank) )
   obj%integerMatrix(1:length) = reshape(int_mat, [length])
   obj%rank           = rank
   obj%dims           = dims
   !
   obj%order = 'F'
   !
END SUBROUTINE qes_init_integerMatrix_3



SUBROUTINE qes_reset_integerMatrix(obj)
   IMPLICIT NONE
   TYPE(integerMatrix_type) :: obj
   INTEGER  :: i

   obj%tagname = ""
   obj%lwrite = .FALSE.
   IF (ALLOCATED(obj%integerMatrix))  DEALLOCATE(obj%integerMatrix)
   IF (ALLOCATED(obj%dims))           DEALLOCATE(obj%dims)
END SUBROUTINE qes_reset_integerMatrix



SUBROUTINE qes_write_scalarQuantity(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(scalarQuantity_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
 
   CALL xml_NewElement(xp, TRIM(obj%tagname) )
      CALL xml_addAttribute( xp, 'Units', TRIM(obj%Units))

      !
      CALL xml_addCharacters(xp, obj%scalarQuantity, fmt ='s16')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite = .FALSE.


END SUBROUTINE qes_reset_scalarQuantity



SUBROUTINE qes_write_general_info(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(general_info_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   

   CALL xml_NewElement(xp, TRIM(obj%tagname) )
      !
      CALL qes_write_xml_format(xp, obj%xml_format)
      !
      CALL qes_write_creator(xp, obj%creator)
      !
      CALL qes_write_created(xp, obj%created)
      !
      CALL xml_NewElement(xp, 'job' )
         CALL xml_addCharacters(xp, TRIM(obj%job) )
      CALL xml_EndElement(xp, 'job' )
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%job = TRIM(job)

END SUBROUTINE qes_init_general_info

SUBROUTINE qes_reset_general_info(obj)
   IMPLICIT NONE
   TYPE(general_info_type) :: obj
   INTEGER  :: i

   obj%tagname = ""
   obj%lwrite = .FALSE.

   CALL qes_reset_xml_format(obj%xml_format)
   CALL qes_reset_creator(obj%creator)
   CALL qes_reset_created(obj%created)

END SUBROUTINE qes_reset_general_info


SUBROUTINE qes_write_parallel_info(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(parallel_info_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   

   CALL xml_NewElement(xp, TRIM(obj%tagname) )
      !
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
   obj%lwrite =.FALSE.

END SUBROUTINE qes_reset_parallel_info


SUBROUTINE qes_write_control_variables(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(control_variables_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN

   CALL xml_NewElement(xp, TRIM(obj%tagname) )
      !
      CALL xml_NewElement(xp, 'title' )
         CALL xml_addCharacters(xp, TRIM(obj%title) )
      CALL xml_EndElement(xp, 'title' )
      CALL xml_NewElement(xp, 'calculation' )
         CALL xml_addCharacters(xp, TRIM(obj%calculation) )
      CALL xml_EndElement(xp, 'calculation' )
      CALL xml_NewElement(xp, 'restart_mode' )
         CALL xml_addCharacters(xp, TRIM(obj%restart_mode))
      CALL xml_EndElement(xp, 'restart_mode')
      CALL xml_NewElement(xp, 'prefix')
         CALL xml_addCharacters(xp, TRIM(obj%prefix))
      CALL xml_EndElement(xp, 'prefix' )
      CALL xml_NewElement(xp, 'pseudo_dir' )
         CALL xml_addCharacters(xp, TRIM(obj%pseudo_dir))
      CALL xml_EndElement(xp, 'pseudo_dir' )
      CALL xml_NewElement(xp, 'outdir' )
         CALL xml_addCharacters(xp, TRIM(obj%outdir) )
      CALL xml_EndElement(xp, 'outdir' )
      CALL xml_NewElement(xp, 'stress' )
         CALL xml_addCharacters(xp, obj%stress)
      CALL xml_EndElement(xp, 'stress' )
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
         CALL xml_addCharacters(xp, obj%max_seconds )
      CALL xml_EndElement(xp, 'max_seconds')
      IF(obj%nstep_ispresent) THEN
         CALL xml_NewElement(xp, 'nstep')
            CALL xml_addCharacters(xp, obj%nstep)
         CALL xml_EndElement(xp, 'nstep')
      ENDIF
      !
      CALL xml_NewElement(xp, 'etot_conv_thr')
         CALL xml_addCharacters(xp,obj%etot_conv_thr )
      CALL xml_EndElement(xp, 'etot_conv_thr')
      CALL xml_NewElement(xp, 'forc_conv_thr')
         CALL xml_addCharacters(xp,obj%forc_conv_thr )
      CALL xml_EndElement(xp, 'forc_conv_thr')
      CALL xml_NewElement(xp, 'press_conv_thr')
         CALL xml_addCharacters(xp, obj%press_conv_thr )
      CALL xml_EndElement(xp, 'press_conv_thr')
      CALL xml_NewElement(xp, 'verbosity')
         CALL xml_addCharacters(xp, TRIM(obj%verbosity) )
      CALL xml_EndElement(xp, 'verbosity')
      CALL xml_NewElement(xp, 'print_every')
         CALL xml_addCharacters(xp, obj%print_every )
      CALL xml_EndElement(xp, 'print_every')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite = .FALSE.
   IF(obj%nstep_ispresent) THEN
      obj%nstep_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_control_variables


SUBROUTINE qes_write_species(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(species_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   
   CALL xml_NewElement(xp, TRIM(obj%tagname) )
   
      CALL xml_addAttribute( xp, 'name', TRIM(obj%name))

      !
      IF(obj%mass_ispresent) THEN
         CALL xml_NewElement(xp, 'mass')
            CALL xml_addCharacters(xp, obj%mass, fmt = 's16')
         CALL xml_EndElement(xp, 'mass')
      ENDIF
      !
      CALL xml_NewElement(xp, 'pseudo_file')
         CALL xml_addCharacters(xp,  TRIM(obj%pseudo_file) )
      CALL xml_EndElement(xp, 'pseudo_file')
      IF(obj%starting_magnetization_ispresent) THEN
         CALL xml_NewElement(xp, 'starting_magnetization')
            CALL xml_addCharacters(xp, obj%starting_magnetization, fmt ='s16')
         CALL xml_EndElement(xp, 'starting_magnetization')
      ENDIF
      !
      IF(obj%spin_teta_ispresent) THEN
         CALL xml_NewElement(xp, 'spin_teta')
            CALL xml_addCharacters(xp, obj%spin_teta, fmt ='s16')
         CALL xml_EndElement(xp, 'spin_teta')
      ENDIF
      !
      IF(obj%spin_phi_ispresent) THEN
         CALL xml_NewElement(xp, 'spin_phi')
            CALL xml_addCharacters(xp, obj%spin_phi, fmt = 's16')
         CALL xml_EndElement(xp, 'spin_phi')
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite = .FALSE.
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


SUBROUTINE qes_write_atomic_positions(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(atomic_positions_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   

   CALL xml_NewElement(xp, TRIM(obj%tagname) )
      !
      DO i = 1, obj%ndim_atom
         CALL qes_write_atom(xp, obj%atom(i))
         !
      END DO
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.
   IF (ALLOCATED (obj%atom)) THEN 
      DO i = 1, SIZE(obj%atom)
         CALL qes_reset_atom(obj%atom(i))
      ENDDO
      DEALLOCATE(obj%atom)
   END IF


END SUBROUTINE qes_reset_atomic_positions


SUBROUTINE qes_write_wyckoff_positions(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(wyckoff_positions_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   
   CALL xml_NewElement(xp, TRIM(obj%tagname) )

   CALL xml_addAttribute( xp, 'space_group', obj%space_group)
   IF(obj%more_options_ispresent) THEN
      CALL xml_addAttribute( xp, 'more_options', TRIM(obj%more_options))
   END IF

      !
      DO i = 1, obj%ndim_atom
         CALL qes_write_atom(xp, obj%atom(i))
         !
      END DO
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.
   IF (ALLOCATED(obj%atom)) THEN 
      DO i = 1, SIZE(obj%atom)
         CALL qes_reset_atom(obj%atom(i))
      ENDDO
      DEALLOCATE(obj%atom)
   END IF
END SUBROUTINE qes_reset_wyckoff_positions


SUBROUTINE qes_write_cell(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(cell_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   

   CALL xml_NewElement(xp, TRIM(obj%tagname) )
      !
      CALL xml_NewElement(xp, 'a1')
         CALL xml_addCharacters(xp, obj%a1, fmt = 's16')
      CALL xml_EndElement(xp, 'a1')
      CALL xml_NewElement(xp, 'a2')
         CALL xml_addCharacters(xp, obj%a2, fmt = 's16')
      CALL xml_EndElement(xp, 'a2')
      CALL xml_NewElement(xp, 'a3')
         CALL xml_addCharacters(xp, obj%a3, fmt = 's16')
      CALL xml_EndElement(xp, 'a3')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

END SUBROUTINE qes_reset_cell


SUBROUTINE qes_write_hybrid(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(hybrid_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   

   CALL xml_NewElement(xp, TRIM(obj%tagname) )
      !
      CALL qes_write_qpoint_grid(xp, obj%qpoint_grid)
      !
      CALL xml_NewElement(xp, 'ecutfock')
         CALL xml_addCharacters(xp,obj%ecutfock, fmt = 's16')
      CALL xml_EndElement(xp, 'ecutfock')
      CALL xml_NewElement(xp, 'exx_fraction')
         CALL xml_addCharacters(xp, obj%exx_fraction, fmt = 's16')
      CALL xml_EndElement(xp, 'exx_fraction')
      CALL xml_NewElement(xp, 'screening_parameter')
         CALL xml_addCharacters(xp, obj%screening_parameter, fmt = 's16')
      CALL xml_EndElement(xp, 'screening_parameter')
      CALL xml_NewElement(xp, 'exxdiv_treatment')
         CALL xml_addCharacters(xp, TRIM(obj%exxdiv_treatment) )
      CALL xml_EndElement(xp, 'exxdiv_treatment')
      CALL xml_NewElement(xp, 'x_gamma_extrapolation')
            CALL xml_addCharacters(xp, obj%x_gamma_extrapolation)
      CALL xml_EndElement(xp, 'x_gamma_extrapolation')
      CALL xml_NewElement(xp, 'ecutvcut')
         CALL xml_addCharacters(xp, obj%ecutvcut, fmt = 's16')
      CALL xml_EndElement(xp, 'ecutvcut')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  =.FALSE.
   CALL qes_reset_qpoint_grid(obj%qpoint_grid)

END SUBROUTINE qes_reset_hybrid


SUBROUTINE qes_write_dftU(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(dftU_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      IF(obj%lda_plus_u_kind_ispresent) THEN
         CALL xml_NewElement(xp, 'lda_plus_u_kind')
            CALL xml_addCharacters(xp, obj%lda_plus_u_kind)
         CALL xml_EndElement(xp, 'lda_plus_u_kind')
       ENDIF
      !
      IF(obj%Hubbard_U_ispresent) THEN
         DO i = 1, obj%ndim_Hubbard_U
            CALL qes_write_HubbardCommon(xp, obj%Hubbard_U(i))
            !
         END DO
      ENDIF
      !
      IF(obj%Hubbard_J0_ispresent) THEN
         DO i = 1, obj%ndim_Hubbard_J0
            CALL qes_write_HubbardCommon(xp, obj%Hubbard_J0(i))
            !
         END DO
      ENDIF
      !
      IF(obj%Hubbard_alpha_ispresent) THEN
         DO i = 1, obj%ndim_Hubbard_alpha
            CALL qes_write_HubbardCommon(xp, obj%Hubbard_alpha(i))
            !
         END DO
      ENDIF
      !
      IF(obj%Hubbard_beta_ispresent) THEN
         DO i = 1, obj%ndim_Hubbard_beta
            CALL qes_write_HubbardCommon(xp, obj%Hubbard_beta(i))
            !
         END DO
      ENDIF
      !
      IF(obj%Hubbard_J_ispresent) THEN
         DO i = 1, obj%ndim_Hubbard_J
            CALL qes_write_HubbardJ(xp, obj%Hubbard_J(i))
            !
         END DO
      ENDIF
      !
      IF(obj%starting_ns_ispresent) THEN
         DO i = 1, obj%ndim_starting_ns
            CALL qes_write_starting_ns(xp, obj%starting_ns(i))
            !
         END DO
      ENDIF
      !
      IF(obj%Hubbard_ns_ispresent) THEN
         DO i = 1, obj%ndim_Hubbard_ns
            CALL qes_write_Hubbard_ns(xp, obj%Hubbard_ns(i))
            !
         END DO
      ENDIF
      !
      IF(obj%U_projection_type_ispresent) THEN
         CALL xml_NewElement(xp, 'U_projection_type')
            CALL xml_addCharacters(xp, TRIM(obj%U_projection_type) )
         CALL xml_EndElement(xp, 'U_projection_type')
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite = .FALSE.
   IF(obj%lda_plus_u_kind_ispresent) THEN
      obj%lda_plus_u_kind_ispresent = .FALSE.
   ENDIF
   IF(obj%Hubbard_U_ispresent) THEN
      IF (ALLOCATED(obj%Hubbard_U)) THEN 
         DO i = 1, SIZE(obj%Hubbard_U)
            CALL qes_reset_HubbardCommon(obj%Hubbard_U(i))
         ENDDO
         DEALLOCATE(obj%Hubbard_U)
      END IF
      obj%Hubbard_U_ispresent = .FALSE.
   ENDIF
   IF(obj%Hubbard_J0_ispresent) THEN
      IF (ALLOCATED(obj%Hubbard_J0)) THEN 
         DO i = 1, SIZE(obj%Hubbard_J0)
            CALL qes_reset_HubbardCommon(obj%Hubbard_J0(i))
         ENDDO
         DEALLOCATE(obj%Hubbard_J0)
      END IF
      obj%Hubbard_J0_ispresent = .FALSE.
   ENDIF
   IF(obj%Hubbard_alpha_ispresent) THEN
      IF (ALLOCATED(obj%Hubbard_alpha)) THEN
         DO i = 1, SIZE(obj%Hubbard_alpha)
            CALL qes_reset_HubbardCommon(obj%Hubbard_alpha(i))
         ENDDO
         DEALLOCATE(obj%Hubbard_alpha)
      END IF
      obj%Hubbard_alpha_ispresent = .FALSE.
   ENDIF
   IF(obj%Hubbard_beta_ispresent) THEN
      IF (ALLOCATED(obj%Hubbard_beta)) THEN
         DO i = 1, SIZE(obj%Hubbard_beta)
            CALL qes_reset_HubbardCommon(obj%Hubbard_beta(i))
         ENDDO
         DEALLOCATE(obj%Hubbard_beta)
      END IF
      obj%Hubbard_beta_ispresent = .FALSE.
   ENDIF
   IF(obj%Hubbard_J_ispresent) THEN
      IF (ALLOCATED(obj%Hubbard_J)) THEN 
         DO i = 1, SIZE(obj%Hubbard_J)
            CALL qes_reset_HubbardJ(obj%Hubbard_J(i))
         ENDDO
         DEALLOCATE(obj%Hubbard_J)
      END IF
      obj%Hubbard_J_ispresent = .FALSE.
   ENDIF
   IF(obj%starting_ns_ispresent) THEN
      IF (ALLOCATED(obj%starting_ns)) THEN
         DO i = 1, SIZE(obj%starting_ns)
            CALL qes_reset_starting_ns(obj%starting_ns(i))
         ENDDO
         DEALLOCATE(obj%starting_ns)
      END IF
      obj%starting_ns_ispresent = .FALSE.
   ENDIF
   IF(obj%Hubbard_ns_ispresent) THEN
      IF (ALLOCATED(obj%Hubbard_ns)) THEN 
         DO i = 1, SIZE(obj%Hubbard_ns)
            CALL qes_reset_Hubbard_ns(obj%Hubbard_ns(i))
         ENDDO
         DEALLOCATE(obj%Hubbard_ns)
      END IF
      obj%Hubbard_ns_ispresent = .FALSE.
   ENDIF
   IF(obj%U_projection_type_ispresent) THEN
      obj%U_projection_type_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_dftU


SUBROUTINE qes_write_vdW(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(vdW_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'vdw_corr')
         CALL xml_addCharacters(xp,  TRIM(obj%vdw_corr) )
      CALL xml_EndElement(xp, 'vdw_corr')
      IF(obj%non_local_term_ispresent) THEN
         CALL xml_NewElement(xp, 'non_local_term')
            CALL xml_addCharacters(xp, TRIM(obj%non_local_term) )
         CALL xml_EndElement(xp, 'non_local_term')
      ENDIF
      !
      IF(obj%london_s6_ispresent) THEN
         CALL xml_NewElement(xp, 'london_s6')
            CALL xml_addCharacters(xp, obj%london_s6, fmt = 's16') 
         CALL xml_EndElement(xp, 'london_s6')
      ENDIF
      !
      IF(obj%ts_vdw_econv_thr_ispresent) THEN
         CALL xml_NewElement(xp, 'ts_vdw_econv_thr')
            CALL xml_addCharacters(xp, obj%ts_vdw_econv_thr, fmt ='s16')
         CALL xml_EndElement(xp, 'ts_vdw_econv_thr')
      ENDIF
      !
      IF(obj%ts_vdw_isolated_ispresent) THEN
         CALL xml_NewElement(xp, 'ts_vdw_isolated')
               CALL xml_addCharacters(xp, obj%ts_vdw_isolated )
         CALL xml_EndElement(xp, 'ts_vdw_isolated')
      ENDIF
      !
      IF(obj%london_rcut_ispresent) THEN
         CALL xml_NewElement(xp, 'london_rcut')
            CALL xml_addCharacters(xp, obj%london_rcut, fmt ='s16')
         CALL xml_EndElement(xp, 'london_rcut')
      ENDIF
      !
      IF(obj%xdm_a1_ispresent) THEN
         CALL xml_NewElement(xp, 'xdm_a1')
            CALL xml_addCharacters(xp, obj%xdm_a1, fmt = 's16')
         CALL xml_EndElement(xp, 'xdm_a1')
      ENDIF
      !
      IF(obj%xdm_a2_ispresent) THEN
         CALL xml_NewElement(xp, 'xdm_a2')
            CALL xml_addCharacters(xp, obj%xdm_a2, fmt = 's16')
         CALL xml_EndElement(xp, 'xdm_a2')
      ENDIF
      !
      IF(obj%london_c6_ispresent) THEN
         DO i = 1, obj%ndim_london_c6
            CALL qes_write_HubbardCommon(xp, obj%london_c6(i))
            !
         END DO
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

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
      IF (ALLOCATED(obj%london_c6)) THEN 
         DO i = 1, SIZE(obj%london_c6)
            CALL qes_reset_HubbardCommon(obj%london_c6(i))
         ENDDO
         DEALLOCATE(obj%london_c6)
      END IF
      obj%london_c6_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_vdW


SUBROUTINE qes_write_spin(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(spin_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'lsda')
            CALL xml_addCharacters(xp,  obj%lsda)      
      CALL xml_EndElement(xp, 'lsda')
      CALL xml_NewElement(xp, 'noncolin')
            CALL xml_addCharacters(xp, obj%noncolin) 
      CALL xml_EndElement(xp, 'noncolin')
      CALL xml_NewElement(xp, 'spinorbit' )
            CALL xml_addCharacters(xp,  obj%spinorbit)
      CALL xml_EndElement(xp, 'spinorbit')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite = .FALSE.

END SUBROUTINE qes_reset_spin

SUBROUTINE qes_write_bands(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(bands_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      IF(obj%nbnd_ispresent) THEN
         CALL xml_NewElement(xp, 'nbnd')
            CALL xml_addCharacters(xp, obj%nbnd )
         CALL xml_EndElement(xp, 'nbnd')
      ENDIF
      !
      IF(obj%smearing_ispresent) THEN
         CALL qes_write_smearing(xp, obj%smearing)
         !
      ENDIF
      !
      IF(obj%tot_charge_ispresent) THEN
         CALL xml_NewElement(xp, 'tot_charge')
            CALL xml_addCharacters(xp, obj%tot_charge, fmt = 's16' )
         CALL xml_EndElement(xp, 'tot_charge')
      ENDIF
      !
      IF(obj%tot_magnetization_ispresent) THEN
         CALL xml_NewElement(xp, 'tot_magnetization')
            CALL xml_addCharacters(xp, obj%tot_magnetization, fmt = 's16')
         CALL xml_EndElement(xp, 'tot_magnetization')
      ENDIF
      !
      CALL qes_write_occupations(xp, obj%occupations)
      !
      IF(obj%inputOccupations_ispresent) THEN
         DO i = 1, obj%ndim_inputOccupations
            CALL qes_write_inputOccupations(xp, obj%inputOccupations(i))
            !
         END DO
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

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
      IF (ALLOCATED(obj%inputOccupations)) THEN 
         DO i = 1, SIZE(obj%inputOccupations)
            CALL qes_reset_inputOccupations(obj%inputOccupations(i))
         ENDDO
         DEALLOCATE(obj%inputOccupations)
      END IF
      obj%inputOccupations_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_bands


SUBROUTINE qes_write_basis(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(basis_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      IF(obj%gamma_only_ispresent) THEN
         CALL xml_NewElement(xp, 'gamma_only' )
               CALL xml_addCharacters(xp, obj%gamma_only ) 
         CALL xml_EndElement(xp, 'gamma_only')
      ENDIF
      !
      CALL xml_NewElement(xp, 'ecutwfc')
         CALL xml_addCharacters(xp, obj%ecutwfc, fmt = 's16')
      CALL xml_EndElement(xp, 'ecutwfc')
      IF(obj%ecutrho_ispresent) THEN
         CALL xml_NewElement(xp, 'ecutrho')
            CALL xml_addCharacters(xp, obj%ecutrho, fmt = 's16')
         CALL xml_EndElement(xp, 'ecutrho')
      ENDIF
      !
      IF(obj%fft_grid_ispresent) THEN
         CALL qes_write_basisSetItem(xp, obj%fft_grid)
         !
      ENDIF
      !
      IF(obj%fft_smooth_ispresent) THEN
         CALL qes_write_basisSetItem(xp, obj%fft_smooth)
         !
      ENDIF
      !
      IF(obj%fft_box_ispresent) THEN
         CALL qes_write_basisSetItem(xp, obj%fft_box)
         !
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

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


SUBROUTINE qes_write_reciprocal_lattice(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(reciprocal_lattice_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'b1')
         CALL xml_addCharacters(xp, obj%b1, fmt = 's16')
      CALL xml_EndElement(xp, 'b1')
      CALL xml_NewElement(xp, 'b2')
         CALL xml_addCharacters(xp, obj%b2, fmt = 's16')
      CALL xml_EndElement(xp, 'b2')
      CALL xml_NewElement(xp, 'b3')
         CALL xml_addCharacters(xp, obj%b3, fmt = 's16')
      CALL xml_EndElement(xp, 'b3')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  =.FALSE.


END SUBROUTINE qes_reset_reciprocal_lattice


SUBROUTINE qes_write_electron_control(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(electron_control_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'diagonalization' )
         CALL xml_addCharacters(xp,  TRIM(obj%diagonalization))
      CALL xml_EndElement(xp, 'diagonalization')
      CALL xml_NewElement(xp, 'mixing_mode' )
         CALL xml_addCharacters(xp, TRIM(obj%mixing_mode))
      CALL xml_EndElement(xp, 'mixing_mode')
      CALL xml_NewElement(xp, 'mixing_beta')
         CALL xml_addCharacters(xp, obj%mixing_beta, fmt = 's16')
      CALL xml_EndElement(xp, 'mixing_beta')
      CALL xml_NewElement(xp, 'conv_thr')
         CALL xml_addCharacters(xp, obj%conv_thr, fmt = 's16')
      CALL xml_EndElement(xp, 'conv_thr')
      CALL xml_NewElement(xp, 'mixing_ndim')
         CALL xml_addCharacters(xp, obj%mixing_ndim)
      CALL xml_EndElement(xp, 'mixing_ndim')
      CALL xml_NewElement(xp, 'max_nstep')
         CALL xml_addCharacters(xp, obj%max_nstep)
      CALL xml_EndElement(xp, 'max_nstep')
      CALL xml_NewElement(xp, 'real_space_q' )
            CALL xml_addCharacters(xp, obj%real_space_q )    
      CALL xml_EndElement(xp, 'real_space_q')
      CALL xml_NewElement(xp, 'tq_smoothing' )
            CALL xml_addCharacters(xp, obj%tq_smoothing)
      CALL xml_EndElement(xp, 'tq_smoothing')
      CALL xml_NewElement(xp, 'tbeta_smoothing' )
            CALL xml_addCharacters(xp, obj%tbeta_smoothing) 
      CALL xml_EndElement(xp, 'tbeta_smoothing')
      CALL xml_NewElement(xp, 'diago_thr_init')
         CALL xml_addCharacters(xp, obj%diago_thr_init, fmt = 's16')
      CALL xml_EndElement(xp, 'diago_thr_init')
      CALL xml_NewElement(xp, 'diago_full_acc' )
            CALL xml_addCharacters(xp, obj%diago_full_acc )
      CALL xml_EndElement(xp, 'diago_full_acc')
      CALL xml_NewElement(xp, 'diago_cg_maxiter')
         CALL xml_addCharacters(xp, obj%diago_cg_maxiter)
      CALL xml_EndElement(xp, 'diago_cg_maxiter')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

END SUBROUTINE qes_reset_electron_control


SUBROUTINE qes_write_k_points_IBZ(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(k_points_IBZ_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      IF(obj%monkhorst_pack_ispresent) THEN
         CALL qes_write_monkhorst_pack(xp, obj%monkhorst_pack)
         !
      ENDIF
      !
      IF(obj%nk_ispresent) THEN
         CALL xml_NewElement(xp, 'nk')
            CALL xml_addCharacters(xp, obj%nk)
         CALL xml_EndElement(xp, 'nk')
      ENDIF
      !
      IF(obj%k_point_ispresent) THEN
         DO i = 1, obj%ndim_k_point
            CALL qes_write_k_point(xp, obj%k_point(i))
            !
         END DO
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite = .FALSE.

   IF(obj%monkhorst_pack_ispresent) THEN
      CALL qes_reset_monkhorst_pack(obj%monkhorst_pack)
      obj%monkhorst_pack_ispresent = .FALSE.
   ENDIF
   IF(obj%nk_ispresent) THEN
      obj%nk_ispresent = .FALSE.
   ENDIF
   IF(obj%k_point_ispresent) THEN
      IF (ALLOCATED(obj%k_point)) THEN  
         DO i = 1, SIZE(obj%k_point)
            CALL qes_reset_k_point(obj%k_point(i))
         ENDDO
         DEALLOCATE(obj%k_point)
      END IF
      obj%k_point_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_k_points_IBZ


SUBROUTINE qes_write_bfgs(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(bfgs_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'ndim')
         CALL xml_addCharacters(xp, obj%ndim)
      CALL xml_EndElement(xp, 'ndim')
      CALL xml_NewElement(xp, 'trust_radius_min')
         CALL xml_addCharacters(xp, obj%trust_radius_min, fmt ='s16')
      CALL xml_EndElement(xp, 'trust_radius_min')
      CALL xml_NewElement(xp, 'trust_radius_max')
         CALL xml_addCharacters(xp, obj%trust_radius_max, fmt = 's16')
      CALL xml_EndElement(xp, 'trust_radius_max')
      CALL xml_NewElement(xp, 'trust_radius_init')
         CALL xml_addCharacters(xp, obj%trust_radius_init, fmt = 's16')
      CALL xml_EndElement(xp, 'trust_radius_init')
      CALL xml_NewElement(xp, 'w1')
         CALL xml_addCharacters(xp, obj%w1, fmt = 's16')
      CALL xml_EndElement(xp, 'w1')
      CALL xml_NewElement(xp, 'w2')
         CALL xml_addCharacters(xp, obj%w2, fmt = 's16')
      CALL xml_EndElement(xp, 'w2')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE. 


END SUBROUTINE qes_reset_bfgs


SUBROUTINE qes_write_md(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(md_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'pot_extrapolation' )
         CALL xml_addCharacters(xp, TRIM(obj%pot_extrapolation) )
      CALL xml_EndElement(xp, 'pot_extrapolation')
      CALL xml_NewElement(xp, 'wfc_extrapolation' )
         CALL xml_addCharacters(xp, TRIM(obj%wfc_extrapolation) )
      CALL xml_EndElement(xp, 'wfc_extrapolation')
      CALL xml_NewElement(xp, 'ion_temperature' )
         CALL xml_addCharacters(xp, TRIM(obj%ion_temperature))
      CALL xml_EndElement(xp, 'ion_temperature')
      CALL xml_NewElement(xp, 'timestep')
         CALL xml_addCharacters(xp, obj%timestep, fmt = 's16')
      CALL xml_EndElement(xp, 'timestep')
      CALL xml_NewElement(xp, 'tempw')
         CALL xml_addCharacters(xp, obj%tempw, fmt = 's16')
      CALL xml_EndElement(xp, 'tempw')
      CALL xml_NewElement(xp, 'tolp')
         CALL xml_addCharacters(xp, obj%tolp, 's16')
      CALL xml_EndElement(xp, 'tolp')
      CALL xml_NewElement(xp, 'deltaT')
         CALL xml_addCharacters(xp, obj%deltaT, fmt = 's16')
      CALL xml_EndElement(xp, 'deltaT')
      CALL xml_NewElement(xp, 'nraise')
         CALL xml_addCharacters(xp, obj%nraise)
      CALL xml_EndElement(xp, 'nraise')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

END SUBROUTINE qes_reset_md


SUBROUTINE qes_write_cell_control(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(cell_control_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'cell_dynamics' )
         CALL xml_addCharacters(xp, TRIM(obj%cell_dynamics))
      CALL xml_EndElement(xp, 'cell_dynamics')
      CALL xml_NewElement(xp, 'pressure')
         CALL xml_addCharacters(xp, obj%pressure, fmt = 's16')
      CALL xml_EndElement(xp, 'pressure')
      IF(obj%wmass_ispresent) THEN
         CALL xml_NewElement(xp, 'wmass')
            CALL xml_addCharacters(xp, obj%wmass, fmt = 's16')
         CALL xml_EndElement(xp, 'wmass')
      ENDIF
      !
      IF(obj%cell_factor_ispresent) THEN
         CALL xml_NewElement(xp, 'cell_factor')
            CALL xml_addCharacters(xp, obj%cell_factor, fmt = 's16') 
         CALL xml_EndElement(xp, 'cell_factor')
      ENDIF
      !
      IF(obj%fix_volume_ispresent) THEN
         CALL xml_NewElement(xp, 'fix_volume' )
               CALL xml_addCharacters(xp, obj%fix_volume) 
         CALL xml_EndElement(xp, 'fix_volume')
      ENDIF
      !
      IF(obj%fix_area_ispresent) THEN
         CALL xml_NewElement(xp, 'fix_area' )
               CALL xml_addCharacters(xp, obj%fix_area) 
         CALL xml_EndElement(xp, 'fix_area')
      ENDIF
      !
      IF(obj%isotropic_ispresent) THEN
         CALL xml_NewElement(xp, 'isotropic' )
               CALL xml_addCharacters(xp, obj%isotropic) 
         CALL xml_EndElement(xp, 'isotropic')
      ENDIF
      !
      IF(obj%free_cell_ispresent) THEN
         CALL qes_write_integerMatrix(xp, obj%free_cell)
         !
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

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

SUBROUTINE qes_write_symmetry_flags(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(symmetry_flags_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'nosym' )
            CALL xml_addCharacters(xp, obj%nosym)
      CALL xml_EndElement(xp, 'nosym')
      CALL xml_NewElement(xp, 'nosym_evc' )
            CALL xml_addCharacters(xp, obj%nosym_evc)  
      CALL xml_EndElement(xp, 'nosym_evc')
      CALL xml_NewElement(xp, 'noinv' )
            CALL xml_addCharacters(xp, obj%noinv ) 
      CALL xml_EndElement(xp, 'noinv')
      CALL xml_NewElement(xp, 'no_t_rev' )
            CALL xml_addCharacters(xp, obj%no_t_rev ) 
      CALL xml_EndElement(xp, 'no_t_rev')
      CALL xml_NewElement(xp, 'force_symmorphic' )
            CALL xml_addCharacters(xp,obj%force_symmorphic)
      CALL xml_EndElement(xp, 'force_symmorphic')
      CALL xml_NewElement(xp, 'use_all_frac' )
            CALL xml_addCharacters(xp, obj%use_all_frac ) 
      CALL xml_EndElement(xp, 'use_all_frac')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

END SUBROUTINE qes_reset_symmetry_flags

SUBROUTINE qes_write_esm(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(esm_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'bc' )
         CALL xml_addCharacters(xp, TRIM(obj%bc) )
      CALL xml_EndElement(xp, 'bc')
      CALL xml_NewElement(xp, 'nfit')
         CALL xml_addCharacters(xp, obj%nfit)
      CALL xml_EndElement(xp, 'nfit')
      CALL xml_NewElement(xp, 'w')
         CALL xml_addCharacters(xp, obj%w, fmt = 's16') 
      CALL xml_EndElement(xp, 'w')
      CALL xml_NewElement(xp, 'efield')
         CALL xml_addCharacters(xp, obj%efield, fmt = 's16')
      CALL xml_EndElement(xp, 'efield')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.


END SUBROUTINE qes_reset_esm

SUBROUTINE qes_write_ekin_functional(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(ekin_functional_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'ecfixed')
         CALL xml_addCharacters(xp, obj%ecfixed, fmt = 's16')
      CALL xml_EndElement(xp, 'ecfixed')
      CALL xml_NewElement(xp, 'qcutz')
         CALL xml_addCharacters(xp, obj%qcutz, fmt = 's16')
      CALL xml_EndElement(xp, 'qcutz')
      CALL xml_NewElement(xp, 'q2sigma')
         CALL xml_addCharacters(xp, obj%q2sigma, fmt = 's16')
      CALL xml_EndElement(xp, 'q2sigma')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.


END SUBROUTINE qes_reset_ekin_functional

SUBROUTINE qes_write_spin_constraints(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(spin_constraints_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'spin_constraints' )
         CALL xml_addCharacters(xp, TRIM(obj%spin_constraints) )
      CALL xml_EndElement(xp, 'spin_constraints')
      CALL xml_NewElement(xp, 'lagrange_multiplier')
         CALL xml_addCharacters(xp, obj%lagrange_multiplier, fmt = 's16' ) 
      CALL xml_EndElement(xp, 'lagrange_multiplier')
      IF(obj%target_magnetization_ispresent) THEN
         CALL xml_NewElement(xp, 'target_magnetization')
            CALL xml_addCharacters(xp, obj%target_magnetization, fmt = 's16')
         CALL xml_EndElement(xp, 'target_magnetization')
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

   IF(obj%target_magnetization_ispresent) THEN
      obj%target_magnetization_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_spin_constraints

SUBROUTINE qes_write_gate_settings(xp, obj)
   IMPLICIT NONE
   TYPE(xmlf_t)      :: xp
   TYPE(gate_settings_type)  :: obj 
   !
   IF ( .NOT. obj%lwrite ) RETURN 
   CALL  xml_NewElement(xp,TRIM(obj%tagname)) 
      CALL  xml_NewElement(xp,"use_gate")
         CALL xml_addCharacters(xp, obj%use_gate)
      CALL xml_endElement(xp,"use_gate")
      IF ( obj%zgate_ispresent ) THEN 
         CALL xml_NewElement(xp, "zgate") 
            CALL xml_addCharacters(xp, obj%zgate)
         CALL xml_EndElement(xp,"zgate") 
      END IF
      IF ( obj%relaxz_ispresent) THEN 
         CALL xml_NewElement(xp,"relaxz")
            CALL xml_addCharacters(xp, obj%relaxz)
         CALL xml_EndElement(xp, "relaxz")
      ENDIF
      IF ( obj%block_ispresent) THEN 
         CALL xml_NewElement(xp,"block")
         CALL xml_addCharacters(xp, obj%block)
         CALL xml_EndElement(xp, "block")
      ENDIF
      IF ( obj%block_1_ispresent) THEN 
         CALL xml_NewElement(xp,"block_1")
         CALL xml_addCharacters(xp, obj%block_1)
         CALL xml_EndElement(xp, "block_1")
      ENDIF
      IF ( obj%block_2_ispresent) THEN 
         CALL xml_NewElement(xp,"block_2")
         CALL xml_addCharacters(xp, obj%block_2)
         CALL xml_EndElement(xp, "block_2")
      ENDIF
      IF ( obj%block_height_ispresent) THEN 
         CALL xml_NewElement(xp,"block_height")
         CALL xml_addCharacters(xp, obj%block_height)
         CALL xml_EndElement(xp, "block_height")
      ENDIF
   CALL xml_endElement(xp, TRIM(obj%tagname)) 
END SUBROUTINE qes_write_gate_settings

SUBROUTINE qes_init_gate_settings( obj, tagname, use_gate, zgate, relaxz, block, block_1, block_2, block_height )
   IMPLICIT NONE
   TYPE(gate_settings_type)    :: obj
   CHARACTER(LEN=*)            :: tagname
   LOGICAL                     :: use_gate
   REAL(DP),OPTIONAL           :: zgate
   LOGICAL,OPTIONAL            :: relaxz
   LOGICAL,OPTIONAL            :: block
   REAL(DP),OPTIONAL           :: block_1
   REAL(DP),OPTIONAL           :: block_2
   REAL(DP),OPTIONAL           :: block_height
   !
   obj%tagname = TRIM(tagname)
   obj%use_gate = use_gate
   obj%relaxz_ispresent = PRESENT(relaxz)
   IF( obj%relaxz_ispresent ) obj%relaxz = relaxz
   obj%zgate_ispresent = PRESENT(zgate)
   IF (obj%zgate_ispresent) obj%zgate = zgate
   obj%block_ispresent = PRESENT(block)
   IF( obj%block_ispresent) obj%block = block
   obj%block_1_ispresent = PRESENT(block_1)
   IF (obj%block_1_ispresent) obj%block_1 = block_1
   obj%block_2_ispresent = PRESENT(block_2)
   IF (obj%block_2_ispresent) obj%block_2 = block_2
   obj%block_height_ispresent = PRESENT(block_height)
   IF (obj%block_height_ispresent) obj%block_height = block_height
   ! 
   obj%lwrite = .TRUE.
   obj%lread =.TRUE.
END SUBROUTINE qes_init_gate_settings

SUBROUTINE qes_reset_gate_settings(obj)
   IMPLICIT NONE
   TYPE(gate_settings_type)  :: obj
   obj%tagname = ""
   obj%lwrite = .FALSE.
   obj%lread  = .FALSE.
END SUBROUTINE qes_reset_gate_settings


SUBROUTINE qes_write_electric_field(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(electric_field_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'electric_potential' )
         CALL xml_addCharacters(xp, TRIM(obj%electric_potential))
      CALL xml_EndElement(xp, 'electric_potential')
      IF(obj%dipole_correction_ispresent) THEN
         CALL xml_NewElement(xp, 'dipole_correction' )
               CALL xml_addCharacters(xp, obj%dipole_correction) 
         CALL xml_EndElement(xp, 'dipole_correction')
      ENDIF
      !
      IF(obj%electric_field_direction_ispresent) THEN
         CALL xml_NewElement(xp, 'electric_field_direction')
            CALL xml_addCharacters(xp, obj%electric_field_direction)
         CALL xml_EndElement(xp, 'electric_field_direction')
      ENDIF
      !
      IF(obj%potential_max_position_ispresent) THEN
         CALL xml_NewElement(xp, 'potential_max_position')
            CALL xml_addCharacters(xp, obj%potential_max_position, fmt = 's16')
         CALL xml_EndElement(xp, 'potential_max_position')
      ENDIF
      !
      IF(obj%potential_decrease_width_ispresent) THEN
         CALL xml_NewElement(xp, 'potential_decrease_width')
            CALL xml_addCharacters(xp, obj%potential_decrease_width, fmt = 's16')
         CALL xml_EndElement(xp, 'potential_decrease_width')
      ENDIF
      !
      IF(obj%electric_field_amplitude_ispresent) THEN
         CALL xml_NewElement(xp, 'electric_field_amplitude')
            CALL xml_addCharacters(xp, obj%electric_field_amplitude, fmt = 's16')
         CALL xml_EndElement(xp, 'electric_field_amplitude')
      ENDIF
      !
      IF(obj%electric_field_vector_ispresent) THEN
         CALL xml_NewElement(xp, 'electric_field_vector')
            CALL xml_addCharacters(xp, obj%electric_field_vector, fmt = 's16') 
         CALL xml_EndElement(xp, 'electric_field_vector')
      ENDIF
      !
      IF(obj%nk_per_string_ispresent) THEN
         CALL xml_NewElement(xp, 'nk_per_string')
            CALL xml_addCharacters(xp, obj%nk_per_string)
         CALL xml_EndElement(xp, 'nk_per_string')
      ENDIF
      !
      IF(obj%n_berry_cycles_ispresent) THEN
         CALL xml_NewElement(xp, 'n_berry_cycles')
            CALL xml_addCharacters(xp, obj%n_berry_cycles)
         CALL xml_EndElement(xp, 'n_berry_cycles')
      ENDIF
      IF ( obj%gate_settings_ispresent) CALL qes_write_gate_settings(xp, obj%gate_settings) 
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
                              n_berry_cycles_ispresent, n_berry_cycles, gate_settings)
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
   TYPE(gate_settings_type), OPTIONAL  :: gate_settings

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
   obj%gate_settings_ispresent = PRESENT( gate_settings) 
   IF ( obj%gate_settings_ispresent ) obj%gate_settings = gate_settings 
END SUBROUTINE qes_init_electric_field

SUBROUTINE qes_reset_electric_field(obj)
   IMPLICIT NONE
   TYPE(electric_field_type) :: obj
   INTEGER  :: i

   obj%tagname = ""
   obj%lwrite = .FALSE.

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
   IF (obj%gate_settings_ispresent) THEN
      obj%gate_settings_ispresent = .FALSE.
      CALL qes_reset_gate_settings(obj%gate_settings) 
   END IF
END SUBROUTINE qes_reset_electric_field


SUBROUTINE qes_write_atomic_constraint(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)                 :: xp
   TYPE(atomic_constraint_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
  
   CALL xml_NewElement(xp, TRIM(obj%tagname)) 
      !
      CALL xml_NewElement(xp, 'constr_parms') 
         CALL xml_addCharacters(xp, obj%constr_parms, fmt = 's16') 
      CALL xml_EndElement(xp, 'constr_parms') 
      CALL xml_NewElement(xp, 'constr_type') 
         CALL xml_addCharacters(xp, TRIM(obj%constr_type) ) 
      CALL xml_EndElement(xp, 'constr_type' )
      CALL xml_NewElement(xp, 'constr_target')
         CALL xml_addCharacters(xp, obj%constr_target, fmt = 's16') 
      CALL xml_EndElement(xp, 'constr_target')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.


END SUBROUTINE qes_reset_atomic_constraint



SUBROUTINE qes_write_dipoleOutput(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(dipoleOutput_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'idir')
         CALL xml_addCharacters(xp, obj%idir)
      CALL xml_EndElement(xp, 'idir')
      CALL qes_write_scalarQuantity(xp, obj%dipole)
      !
      CALL qes_write_scalarQuantity(xp, obj%ion_dipole)
      !
      CALL qes_write_scalarQuantity(xp, obj%elec_dipole)
      !
      CALL qes_write_scalarQuantity(xp, obj%dipoleField)
      !
      CALL qes_write_scalarQuantity(xp, obj%potentialAmp)
      !
      CALL qes_write_scalarQuantity(xp, obj%totalLength)
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

   CALL qes_reset_scalarQuantity(obj%dipole)
   CALL qes_reset_scalarQuantity(obj%ion_dipole)
   CALL qes_reset_scalarQuantity(obj%elec_dipole)
   CALL qes_reset_scalarQuantity(obj%dipoleField)
   CALL qes_reset_scalarQuantity(obj%potentialAmp)
   CALL qes_reset_scalarQuantity(obj%totalLength)

END SUBROUTINE qes_reset_dipoleOutput


SUBROUTINE qes_write_finiteFieldOut(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(finiteFieldOut_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'electronicDipole')
         CALL xml_addCharacters(xp, obj%electronicDipole, fmt = 's16')
      CALL xml_EndElement(xp, 'electronicDipole')
      CALL xml_NewElement(xp, 'ionicDipole')
         CALL xml_addCharacters(xp, obj%ionicDipole, fmt = 's16')
      CALL xml_EndElement(xp, 'ionicDipole')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE. 


END SUBROUTINE qes_reset_finiteFieldOut


SUBROUTINE qes_write_polarization(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(polarization_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL qes_write_scalarQuantity(xp, obj%polarization)
      !
      CALL xml_NewElement(xp, 'modulus')
         CALL xml_addCharacters(xp, obj%modulus, fmt = 's16')
      CALL xml_EndElement(xp, 'modulus')
      CALL xml_NewElement(xp, 'direction')
         CALL xml_addCharacters(xp, obj%direction, fmt = 's16')
      CALL xml_EndElement(xp, 'direction')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

   CALL qes_reset_scalarQuantity(obj%polarization)

END SUBROUTINE qes_reset_polarization


SUBROUTINE qes_write_ionicPolarization(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(ionicPolarization_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL qes_write_atom(xp, obj%ion)
      !
      CALL xml_NewElement(xp, 'charge')
         CALL xml_addCharacters(xp, obj%charge, fmt = 's16')
      CALL xml_EndElement(xp, 'charge')
      CALL qes_write_phase(xp, obj%phase)
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

   CALL qes_reset_atom(obj%ion)
   CALL qes_reset_phase(obj%phase)

END SUBROUTINE qes_reset_ionicPolarization


SUBROUTINE qes_write_electronicPolarization(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(electronicPolarization_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL qes_write_k_point(xp, obj%firstKeyPoint)
      !
      IF(obj%spin_ispresent) THEN
         CALL xml_NewElement(xp, 'spin')
            CALL xml_addCharacters(xp, obj%spin)
         CALL xml_EndElement(xp, 'spin')
      ENDIF
      !
      CALL qes_write_phase(xp, obj%phase)
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

   CALL qes_reset_k_point(obj%firstKeyPoint)
   IF(obj%spin_ispresent) THEN
      obj%spin_ispresent = .FALSE.
   ENDIF
   CALL qes_reset_phase(obj%phase)

END SUBROUTINE qes_reset_electronicPolarization


SUBROUTINE qes_write_scf_conv(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(scf_conv_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'n_scf_steps')
         CALL xml_addCharacters(xp, obj%n_scf_steps)
      CALL xml_EndElement(xp, 'n_scf_steps')
      CALL xml_NewElement(xp, 'scf_error')
         CALL xml_addCharacters(xp, obj%scf_error, fmt = 's16')
      CALL xml_EndElement(xp, 'scf_error')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.


END SUBROUTINE qes_reset_scf_conv


SUBROUTINE qes_write_opt_conv(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(opt_conv_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'n_opt_steps')
         CALL xml_addCharacters(xp, obj%n_opt_steps)
      CALL xml_EndElement(xp, 'n_opt_steps')
      CALL xml_NewElement(xp, 'grad_norm')
         CALL xml_addCharacters(xp, obj%grad_norm, fmt = 's16')
      CALL xml_EndElement(xp, 'grad_norm')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.


END SUBROUTINE qes_reset_opt_conv


SUBROUTINE qes_write_algorithmic_info(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(algorithmic_info_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'real_space_q' )
            CALL xml_addCharacters(xp, obj%real_space_q ) 
      CALL xml_EndElement(xp, 'real_space_q')
      CALL xml_NewElement(xp, 'uspp' )
            CALL xml_addCharacters(xp, obj%uspp)
      CALL xml_EndElement(xp, 'uspp')
      CALL xml_NewElement(xp, 'paw' )
            CALL xml_addCharacters(xp, obj%paw)
      CALL xml_EndElement(xp, 'paw')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

END SUBROUTINE qes_reset_algorithmic_info


SUBROUTINE qes_write_symmetry(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(symmetry_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL qes_write_info(xp, obj%info)
      !
      CALL qes_write_matrix(xp, obj%rotation)
      !
      IF(obj%fractional_translation_ispresent) THEN
         CALL xml_NewElement(xp, 'fractional_translation')
            CALL xml_addCharacters(xp, obj%fractional_translation, fmt = 's16')
         CALL xml_EndElement(xp, 'fractional_translation')
      ENDIF
      !
      IF(obj%equivalent_atoms_ispresent) THEN
         CALL qes_write_equivalent_atoms(xp, obj%equivalent_atoms)
         !
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

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


SUBROUTINE qes_write_magnetization(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(magnetization_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'lsda' )
            CALL xml_addCharacters(xp, obj%lsda) 
      CALL xml_EndElement(xp, 'lsda')
      CALL xml_NewElement(xp, 'noncolin' )
            CALL xml_addCharacters(xp, obj%noncolin)
      CALL xml_EndElement(xp, 'noncolin')
      CALL xml_NewElement(xp, 'spinorbit' )
            CALL xml_addCharacters(xp, obj%spinorbit)
      CALL xml_EndElement(xp, 'spinorbit')
      CALL xml_NewElement(xp, 'total')
         CALL xml_addCharacters(xp, obj%total, fmt = 's16')
      CALL xml_EndElement(xp, 'total')
      CALL xml_NewElement(xp, 'absolute')
         CALL xml_addCharacters(xp, obj%absolute, fmt = 's16')
      CALL xml_EndElement(xp, 'absolute')
      CALL xml_NewElement(xp, 'do_magnetization' )
            CALL xml_addCharacters(xp, obj%do_magnetization)
      CALL xml_EndElement(xp, 'do_magnetization')
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

END SUBROUTINE qes_reset_magnetization


SUBROUTINE qes_write_total_energy(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(total_energy_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'etot')
         CALL xml_addCharacters(xp, obj%etot, fmt = 's16')
      CALL xml_EndElement(xp, 'etot')
      IF(obj%eband_ispresent) THEN
         CALL xml_NewElement(xp, 'eband')
            CALL xml_addCharacters(xp, obj%eband, fmt ='s16')
         CALL xml_EndElement(xp, 'eband')
      ENDIF
      !
      IF(obj%ehart_ispresent) THEN
         CALL xml_NewElement(xp, 'ehart')
            CALL xml_addCharacters(xp, obj%ehart, fmt = 's16')
         CALL xml_EndElement(xp, 'ehart')
      ENDIF
      !
      IF(obj%vtxc_ispresent) THEN
         CALL xml_NewElement(xp, 'vtxc')
            CALL xml_addCharacters(xp, obj%vtxc, fmt = 's16')
         CALL xml_EndElement(xp, 'vtxc')
      ENDIF
      !
      IF(obj%etxc_ispresent) THEN
         CALL xml_NewElement(xp, 'etxc')
            CALL xml_addCharacters(xp, obj%etxc, fmt = 's16')
         CALL xml_EndElement(xp, 'etxc')
      ENDIF
      !
      IF(obj%ewald_ispresent) THEN
         CALL xml_NewElement(xp, 'ewald')
            CALL xml_addCharacters(xp, obj%ewald, fmt = 's16')
         CALL xml_EndElement(xp, 'ewald')
      ENDIF
      !
      IF(obj%demet_ispresent) THEN
         CALL xml_NewElement(xp, 'demet')
            CALL xml_addCharacters(xp, obj%demet, fmt = 's16')
         CALL xml_EndElement(xp, 'demet')
      ENDIF
      !
      IF(obj%efieldcorr_ispresent) THEN
         CALL xml_NewElement(xp, 'efieldcorr')
            CALL xml_addCharacters(xp, obj%efieldcorr, fmt = 's16')
         CALL xml_EndElement(xp, 'efieldcorr')
      ENDIF
      !
      IF(obj%potentiostat_contr_ispresent) THEN
         CALL xml_NewElement(xp, 'potentiostat_contr')
            CALL xml_addCharacters(xp, obj%potentiostat_contr, fmt = 's16')
         CALL xml_EndElement(xp, 'potentiostat_contr')
      ENDIF
      !
      IF (obj%gatefield_contr_ispresent) THEN 
         CALL xml_NewElement(xp, 'gatefield_contr') 
            CALL xml_addCharacters(xp, obj%gatefield_contr, fmt = 's16') 
         CALL xml_EndElement(xp, 'gatefield_contr') 
      END IF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_total_energy

SUBROUTINE qes_init_total_energy(obj, tagname, etot, eband, ehart, vtxc, etxc, ewald, demet, &
                                 efieldcorr, potentiostat_contr, gate_contribution)
   IMPLICIT NONE

   TYPE(total_energy_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   REAL(DP) :: etot
   REAL(DP),OPTIONAL :: eband
   REAL(DP),OPTIONAL :: ehart
   REAL(DP),OPTIONAL :: vtxc
   REAL(DP),OPTIONAL :: etxc
   REAL(DP),OPTIONAL :: ewald
   REAL(DP),OPTIONAL :: demet
   REAL(DP),OPTIONAL :: efieldcorr
   REAL(DP),OPTIONAL :: potentiostat_contr
   REAL(DP),OPTIONAL :: gate_contribution
   ! 
   !
   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%etot = etot
   obj%eband_ispresent = PRESENT(eband)
   IF(obj%eband_ispresent) THEN
      obj%eband = eband
   ENDIF
   obj%ehart_ispresent = PRESENT(ehart)
   IF(obj%ehart_ispresent) THEN
      obj%ehart = ehart
   ENDIF
   obj%vtxc_ispresent = PRESENT(vtxc)
   IF(obj%vtxc_ispresent) THEN
      obj%vtxc = vtxc
   ENDIF
   obj%etxc_ispresent = PRESENT(etxc)
   IF(obj%etxc_ispresent) THEN
      obj%etxc = etxc
   ENDIF
   obj%ewald_ispresent = PRESENT(ewald)
   IF(obj%ewald_ispresent) THEN
      obj%ewald = ewald
   ENDIF
   obj%demet_ispresent = PRESENT(demet)
   IF(obj%demet_ispresent) THEN
      obj%demet = demet
   ENDIF
   obj%efieldcorr_ispresent = PRESENT(efieldcorr)
   IF(obj%efieldcorr_ispresent) THEN
      obj%efieldcorr = efieldcorr
   ENDIF
   obj%potentiostat_contr_ispresent = PRESENT(potentiostat_contr) 
   IF(obj%potentiostat_contr_ispresent) THEN
      obj%potentiostat_contr = potentiostat_contr
   ENDIF
   obj%gatefield_contr_ispresent = PRESENT(gate_contribution)
   IF (obj%gatefield_contr_ispresent) obj%gatefield_contr=gate_contribution 
END SUBROUTINE qes_init_total_energy

SUBROUTINE qes_reset_total_energy(obj)
   IMPLICIT NONE
   TYPE(total_energy_type) :: obj
   INTEGER  :: i

   obj%tagname = ""
   obj%lwrite  = .FALSE.

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


SUBROUTINE qes_write_ks_energies(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(ks_energies_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL qes_write_k_point(xp, obj%k_point)
      !
      CALL xml_NewElement(xp, 'npw')
         CALL xml_addCharacters(xp, obj%npw)
      CALL xml_EndElement(xp, 'npw')
      !
      CALL qes_write_vector(xp, obj%eigenvalues)
      !
      CALL qes_write_vector(xp, obj%occupations)
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_ks_energies

SUBROUTINE qes_init_ks_energies(obj, tagname, k_point, npw, eigenvalues, occupations)
   IMPLICIT NONE

   TYPE(ks_energies_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   TYPE(k_point_type) :: k_point
   TYPE(vector_type)  :: eigenvalues
   TYPE(vector_type)  :: occupations
   INTEGER  :: npw

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%k_point = k_point
   obj%npw = npw
   obj%eigenvalues =  eigenvalues
   obj%occupations = occupations

END SUBROUTINE qes_init_ks_energies

SUBROUTINE qes_reset_ks_energies(obj)
   IMPLICIT NONE
   TYPE(ks_energies_type) :: obj
   INTEGER  :: i

   obj%tagname = ""
   obj%lwrite  = .FALSE.

   CALL qes_reset_k_point(obj%k_point)
   CALL qes_reset_vector(obj%eigenvalues)
   CALL qes_reset_vector(obj%occupations)

END SUBROUTINE qes_reset_ks_energies


SUBROUTINE qes_write_atomic_species(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(atomic_species_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
   CALL xml_addAttribute( xp, 'ntyp', obj%ntyp)
   IF (obj%pseudo_dir_ispresent) CALL xml_addAttribute( xp, 'pseudo_dir', TRIM ( obj%pseudo_dir) )

      !
      DO i = 1, obj%ndim_species
         CALL qes_write_species(xp, obj%species(i))
         !
      END DO
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  =.FALSE.
   IF (ALLOCATED(obj%species)) THEN
      DO i = 1, SIZE(obj%species)
         CALL qes_reset_species(obj%species(i))
      ENDDO
      DEALLOCATE(obj%species)
   END IF

END SUBROUTINE qes_reset_atomic_species


SUBROUTINE qes_write_atomic_structure(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(atomic_structure_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
   CALL xml_addAttribute( xp, 'nat', obj%nat)
   IF(obj%alat_ispresent) THEN
      CALL xml_addAttribute( xp, 'alat', obj%alat)
   END IF
   IF(obj%bravais_index_ispresent) THEN
      CALL xml_addAttribute( xp, 'bravais_index', obj%bravais_index)
   END IF

      !
      IF(obj%atomic_positions_ispresent) THEN
         CALL qes_write_atomic_positions(xp, obj%atomic_positions)
         !
      ENDIF
      !
      IF(obj%wyckoff_positions_ispresent) THEN
         CALL qes_write_wyckoff_positions(xp, obj%wyckoff_positions)
         !
      ENDIF
      !
      CALL qes_write_cell(xp, obj%cell)
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

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


SUBROUTINE qes_write_dft(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(dft_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'functional' )
         CALL xml_addCharacters(xp,  TRIM(obj%functional) )
      CALL xml_EndElement(xp, 'functional')
      IF(obj%hybrid_ispresent) THEN
         CALL qes_write_hybrid(xp, obj%hybrid)
         !
      ENDIF
      !
      IF(obj%dftU_ispresent) THEN
         CALL qes_write_dftU(xp, obj%dftU)
         !
      ENDIF
      !
      IF(obj%vdW_ispresent) THEN
         CALL qes_write_vdW(xp, obj%vdW)
         !
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

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


SUBROUTINE qes_write_basis_set(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(basis_set_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      IF(obj%gamma_only_ispresent) THEN
         CALL xml_NewElement(xp, 'gamma_only' )
               CALL xml_addCharacters(xp, obj%gamma_only) 
         CALL xml_EndElement(xp, 'gamma_only')
      ENDIF
      !
      CALL xml_NewElement(xp, 'ecutwfc')
         CALL xml_addCharacters(xp, obj%ecutwfc,    fmt = 's16')
      CALL xml_EndElement(xp, 'ecutwfc')
      IF(obj%ecutrho_ispresent) THEN
         CALL xml_NewElement(xp, 'ecutrho')
            CALL xml_addCharacters(xp, obj%ecutrho, fmt = 's16')
         CALL xml_EndElement(xp, 'ecutrho')
      ENDIF
      !
      CALL qes_write_basisSetItem(xp, obj%fft_grid)
      !
      IF(obj%fft_smooth_ispresent) THEN
         CALL qes_write_basisSetItem(xp, obj%fft_smooth)
         !
      ENDIF
      !
      IF(obj%fft_box_ispresent) THEN
         CALL qes_write_basisSetItem(xp, obj%fft_box)
         !
      ENDIF
      !
      CALL xml_NewElement(xp, 'ngm')
         CALL xml_addCharacters(xp, obj%ngm )
      CALL xml_EndElement(xp, 'ngm')
      IF(obj%ngms_ispresent) THEN
         CALL xml_NewElement(xp, 'ngms')
            CALL xml_addCharacters(xp, obj%ngms )
         CALL xml_EndElement(xp, 'ngms')
      ENDIF
      !
      CALL xml_NewElement(xp, 'npwx')
         CALL xml_addCharacters(xp, obj%npwx )
      CALL xml_EndElement(xp, 'npwx')
      CALL qes_write_reciprocal_lattice(xp, obj%reciprocal_lattice)
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

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


SUBROUTINE qes_write_ion_control(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(ion_control_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'ion_dynamics' )
         CALL xml_addCharacters(xp,  TRIM(obj%ion_dynamics))
      CALL xml_EndElement(xp, 'ion_dynamics')
      IF(obj%upscale_ispresent) THEN
         CALL xml_NewElement(xp, 'upscale')
            CALL xml_addCharacters(xp, obj%upscale, fmt = 's16')
         CALL xml_EndElement(xp, 'upscale')
      ENDIF
      !
      IF(obj%remove_rigid_rot_ispresent) THEN
         CALL xml_NewElement(xp, 'remove_rigid_rot' )
            CALL xml_addCharacters(xp,  obj%remove_rigid_rot)
         CALL xml_EndElement(xp, 'remove_rigid_rot')
      ENDIF
      !
      IF(obj%refold_pos_ispresent) THEN
         CALL xml_NewElement(xp, 'refold_pos' )
               CALL xml_addCharacters(xp, obj%refold_pos)  
         CALL xml_EndElement(xp, 'refold_pos')
      ENDIF
      !
      IF(obj%bfgs_ispresent) THEN
         CALL qes_write_bfgs(xp, obj%bfgs)
         !
      ENDIF
      !
      IF(obj%md_ispresent) THEN
         CALL qes_write_md(xp, obj%md)
         !
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.
   !
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


SUBROUTINE qes_write_boundary_conditions(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(boundary_conditions_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'assume_isolated' )
         CALL xml_addCharacters(xp,  TRIM(obj%assume_isolated))
      CALL xml_EndElement(xp, 'assume_isolated')
      IF(obj%esm_ispresent) THEN
         CALL qes_write_esm(xp, obj%esm)
         !
      ENDIF
      !
      IF(obj%fcp_opt_ispresent) THEN
         CALL xml_NewElement(xp, 'fcp_opt' )
               CALL xml_addCharacters(xp, obj%fcp_opt)
         CALL xml_EndElement(xp, 'fcp_opt')
      ENDIF
      !
      IF(obj%fcp_mu_ispresent) THEN
         CALL xml_NewElement(xp, 'fcp_mu')
            CALL xml_addCharacters(xp,obj%fcp_mu, fmt = 's16')
         CALL xml_EndElement(xp, 'fcp_mu')
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  =.FALSE.

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


SUBROUTINE qes_write_atomic_constraints(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(atomic_constraints_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'num_of_constraints')
         CALL xml_addCharacters(xp, obj%num_of_constraints)
      CALL xml_EndElement(xp, 'num_of_constraints')
      CALL xml_NewElement(xp, 'tolerance')
         CALL xml_addCharacters(xp, obj%tolerance, fmt = 's16')
      CALL xml_EndElement(xp, 'tolerance')
      DO i = 1, obj%ndim_atomic_constraint
         CALL qes_write_atomic_constraint(xp, obj%atomic_constraint(i))
         !
      END DO
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.
   IF (ALLOCATED(obj%atomic_constraint)) THEN
      DO i = 1, SIZE(obj%atomic_constraint)
         CALL qes_reset_atomic_constraint(obj%atomic_constraint(i))
      ENDDO
      DEALLOCATE(obj%atomic_constraint)
   END IF

END SUBROUTINE qes_reset_atomic_constraints


SUBROUTINE qes_write_BerryPhaseOutput(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(BerryPhaseOutput_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL qes_write_polarization(xp, obj%totalPolarization)
      !
      CALL qes_write_phase(xp, obj%totalPhase)
      !
      DO i = 1, obj%ndim_ionicPolarization
         CALL qes_write_ionicPolarization(xp, obj%ionicPolarization(i))
         !
      END DO
      DO i = 1, obj%ndim_electronicPolarization
         CALL qes_write_electronicPolarization(xp, obj%electronicPolarization(i))
         !
      END DO
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_BerryPhaseOutput

SUBROUTINE qes_init_BerryPhaseOutput(obj, tagname, totalPolarization, totalPhase, &
                              ndim_ionicPolarization, ionicPolarization, &
                              ndim_electronicPolarization, electronicPolarization)
   IMPLICIT NONE

   TYPE(BerryPhaseOutput_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   TYPE(polarization_type) :: totalPolarization
   TYPE(phase_type) :: totalPhase
   INTEGER  :: ndim_ionicPolarization
   TYPE(ionicPolarization_type ), DIMENSION( ndim_ionicPolarization )  :: ionicPolarization
   INTEGER  :: ndim_electronicPolarization
   TYPE(electronicPolarization_type ), DIMENSION( ndim_electronicPolarization )  :: electronicPolarization

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%totalPolarization = totalPolarization
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
   obj%lwrite  = .FALSE.

   CALL qes_reset_polarization(obj%totalPolarization)
   CALL qes_reset_phase(obj%totalPhase)
   IF (ALLOCATED(obj%ionicPolarization)) THEN 
      DO i = 1, SIZE(obj%ionicPolarization)
         CALL qes_reset_ionicPolarization(obj%ionicPolarization(i))
      ENDDO
      DEALLOCATE(obj%ionicPolarization)
   END IF
   
   IF (ALLOCATED(obj%electronicPolarization)) THEN
      DO i = 1, SIZE(obj%electronicPolarization)
         CALL qes_reset_electronicPolarization(obj%electronicPolarization(i))
      ENDDO
      DEALLOCATE(obj%electronicPolarization)
   END IF

END SUBROUTINE qes_reset_BerryPhaseOutput


SUBROUTINE qes_write_convergence_info(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(convergence_info_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL qes_write_scf_conv(xp, obj%scf_conv)
      !
      IF(obj%opt_conv_ispresent) THEN
         CALL qes_write_opt_conv(xp, obj%opt_conv)
         !
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

   CALL qes_reset_scf_conv(obj%scf_conv)
   IF(obj%opt_conv_ispresent) THEN
      CALL qes_reset_opt_conv(obj%opt_conv)
      obj%opt_conv_ispresent = .FALSE.
   ENDIF

END SUBROUTINE qes_reset_convergence_info


SUBROUTINE qes_write_symmetries(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(symmetries_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'nsym' )
         CALL xml_addCharacters(xp, obj%nsym)
      CALL xml_EndElement(xp, 'nsym')
      CALL xml_NewElement(xp, 'nrot' )
         CALL xml_addCharacters(xp, obj%nrot)
      CALL xml_EndElement(xp, 'nrot')
      CALL xml_NewElement(xp, 'space_group' )
         CALL xml_addCharacters(xp, obj%space_group)
      CALL xml_EndElement(xp, 'space_group')
      DO i = 1, obj%ndim_symmetry
         CALL qes_write_symmetry(xp, obj%symmetry(i))
         !
      END DO
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.
   IF (ALLOCATED(obj%symmetry)) THEN 
      DO i = 1, SIZE(obj%symmetry)
         CALL qes_reset_symmetry(obj%symmetry(i))
      ENDDO
      DEALLOCATE(obj%symmetry)
   END IF

END SUBROUTINE qes_reset_symmetries


SUBROUTINE qes_write_band_structure(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(band_structure_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'lsda' )
            CALL xml_addCharacters(xp, obj%lsda) 
      CALL xml_EndElement(xp, 'lsda')
      CALL xml_NewElement(xp, 'noncolin' )
            CALL xml_addCharacters(xp, obj%noncolin) 
      CALL xml_EndElement(xp, 'noncolin')
      CALL xml_NewElement(xp, 'spinorbit' )
            CALL xml_addCharacters(xp, obj%spinorbit) 
      CALL xml_EndElement(xp, 'spinorbit')
      CALL xml_NewElement(xp, 'nbnd')
         CALL xml_addCharacters(xp, obj%nbnd)
      CALL xml_EndElement(xp, 'nbnd')
      IF(obj%nbnd_up_ispresent) THEN
         CALL xml_NewElement(xp, 'nbnd_up')
            CALL xml_addCharacters(xp, obj%nbnd_up)
         CALL xml_EndElement(xp, 'nbnd_up')
      ENDIF
      !
      IF(obj%nbnd_dw_ispresent) THEN
         CALL xml_NewElement(xp, 'nbnd_dw')
            CALL xml_addCharacters(xp, obj%nbnd_dw)
         CALL xml_EndElement(xp, 'nbnd_dw')
      ENDIF
      !
      CALL xml_NewElement(xp, 'nelec')
         CALL xml_addCharacters(xp, obj%nelec, fmt = 's16')
      CALL xml_EndElement(xp, 'nelec')
      IF(obj%num_of_atomic_wfc_ispresent) THEN
         CALL xml_NewElement(xp, 'num_of_atomic_wfc')
            CALL xml_addCharacters(xp, obj%num_of_atomic_wfc)
         CALL xml_EndElement(xp, 'num_of_atomic_wfc')
      ENDIF
      !
      CALL xml_NewElement(xp, 'wf_collected' )
            CALL xml_addCharacters(xp, obj%wf_collected)
      CALL xml_EndElement(xp, 'wf_collected')
      IF(obj%fermi_energy_ispresent) THEN
         CALL xml_NewElement(xp, 'fermi_energy')
            CALL xml_addCharacters(xp, obj%fermi_energy, fmt = 's16')
         CALL xml_EndElement(xp, 'fermi_energy')
      ENDIF
      !
      IF(obj%highestOccupiedLevel_ispresent) THEN
         CALL xml_NewElement(xp, 'highestOccupiedLevel')
            CALL xml_addCharacters(xp, obj%highestOccupiedLevel, fmt = 's16')
         CALL xml_EndElement(xp, 'highestOccupiedLevel')
      ENDIF
      !
      IF(obj%two_fermi_energies_ispresent) THEN
         CALL xml_NewElement(xp,'two_fermi_energies')
            CALL xml_addCharacters(xp, obj%two_fermi_energies, fmt = 's16')
         CALL xml_EndElement(xp,'two_fermi_energies')
         !
      ENDIF
      !
      CALL qes_write_k_points_IBZ(xp, obj%starting_k_points)
      !
      CALL xml_NewElement(xp, 'nks')
         CALL xml_addCharacters(xp, obj%nks)
      CALL xml_EndElement(xp, 'nks')
      CALL qes_write_occupations(xp, obj%occupations_kind)
      !
      IF(obj%smearing_ispresent) THEN
         CALL qes_write_smearing(xp, obj%smearing)
         !
      ENDIF
      !
      DO i = 1, obj%ndim_ks_energies
         CALL qes_write_ks_energies(xp, obj%ks_energies(i))
         !
      END DO
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_band_structure

SUBROUTINE qes_init_band_structure(obj, tagname, lsda, noncolin, spinorbit, nbnd, &
                              nbnd_up_ispresent, nbnd_up, nbnd_dw_ispresent, nbnd_dw, &
                              nelec, num_of_atomic_wfc_ispresent, num_of_atomic_wfc, &
                              wf_collected, fermi_energy_ispresent, fermi_energy, &
                              highestOccupiedLevel_ispresent, highestOccupiedLevel, &
                              two_fermi_energies_ispresent, two_fermi_energies, &
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
   REAL(DP), DIMENSION(2) :: two_fermi_energies
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
      obj%two_fermi_energies(:) = two_fermi_energies(:)
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
      obj%two_fermi_energies_ispresent = .FALSE.
   ENDIF
   CALL qes_reset_k_points_IBZ(obj%starting_k_points)
   CALL qes_reset_occupations(obj%occupations_kind)
   IF(obj%smearing_ispresent) THEN
      CALL qes_reset_smearing(obj%smearing)
      obj%smearing_ispresent = .FALSE.
   ENDIF
   IF (ALLOCATED(obj%ks_energies)) THEN 
      DO i = 1, SIZE(obj%ks_energies)
         CALL qes_reset_ks_energies(obj%ks_energies(i))
      ENDDO
      DEALLOCATE(obj%ks_energies)
   END IF

END SUBROUTINE qes_reset_band_structure


SUBROUTINE qes_write_input(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(input_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL qes_write_control_variables(xp, obj%control_variables)
      !
      CALL qes_write_atomic_species(xp, obj%atomic_species)
      !
      CALL qes_write_atomic_structure(xp, obj%atomic_structure)
      !
      CALL qes_write_dft(xp, obj%dft)
      !
      CALL qes_write_spin(xp, obj%spin)
      !
      CALL qes_write_bands(xp, obj%bands)
      !
      CALL qes_write_basis(xp, obj%basis)
      !
      CALL qes_write_electron_control(xp, obj%electron_control)
      !
      CALL qes_write_k_points_IBZ(xp, obj%k_points_IBZ)
      !
      CALL qes_write_ion_control(xp, obj%ion_control)
      !
      CALL qes_write_cell_control(xp, obj%cell_control)
      !
      IF(obj%symmetry_flags_ispresent) THEN
         CALL qes_write_symmetry_flags(xp, obj%symmetry_flags)
         !
      ENDIF
      !
      IF(obj%boundary_conditions_ispresent) THEN
         CALL qes_write_boundary_conditions(xp, obj%boundary_conditions)
         !
      ENDIF
      !
      IF(obj%ekin_functional_ispresent) THEN
         CALL qes_write_ekin_functional(xp, obj%ekin_functional)
         !
      ENDIF
      !
      IF(obj%external_atomic_forces_ispresent) THEN
         CALL qes_write_matrix(xp, obj%external_atomic_forces)
         !
      ENDIF
      !
      IF(obj%free_positions_ispresent) THEN
         CALL qes_write_integerMatrix(xp, obj%free_positions)
         !
      ENDIF
      !
      IF(obj%starting_atomic_velocities_ispresent) THEN
         CALL qes_write_matrix(xp, obj%starting_atomic_velocities)
         !
      ENDIF
      !
      IF(obj%electric_field_ispresent) THEN
         CALL qes_write_electric_field(xp, obj%electric_field)
         !
      ENDIF
      !
      IF(obj%atomic_constraints_ispresent) THEN
         CALL qes_write_atomic_constraints(xp, obj%atomic_constraints)
         !
      ENDIF
      !
      IF(obj%spin_constraints_ispresent) THEN
         CALL qes_write_spin_constraints(xp, obj%spin_constraints)
         !
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite  = .FALSE.

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


SUBROUTINE qes_write_step(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(step_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !
   CALL xml_NewElement(xp, TRIM(obj%tagname))
   CALL xml_addAttribute( xp, 'n_step', obj%n_step)

      !
      CALL qes_write_scf_conv(xp, obj%scf_conv)
      !
      CALL qes_write_atomic_structure(xp, obj%atomic_structure)
      !
      CALL qes_write_total_energy(xp, obj%total_energy)
      !
      CALL qes_write_matrix(xp, obj%forces)
      !
      IF(obj%stress_ispresent) THEN
         CALL qes_write_matrix(xp, obj%stress)
         !
      ENDIF
      !
      IF(obj%FCP_force_ispresent) THEN
         CALL xml_NewElement(xp, 'FCP_force')
            CALL xml_addCharacters(xp, obj%FCP_force, fmt = 's16')
         CALL xml_EndElement(xp, 'FCP_force')
      ENDIF
      !
      IF(obj%FCP_tot_charge_ispresent) THEN
         CALL xml_NewElement(xp, 'FCP_tot_charge')
            CALL xml_addCharacters(xp, obj%FCP_tot_charge, fmt = 's16')
         CALL xml_EndElement(xp, 'FCP_tot_charge')
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
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
   obj%lwrite = .FALSE.

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

SUBROUTINE qes_write_gateInfo(xp, obj) 
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE (gateInfo_type) :: obj 
   ! 
   INTEGER :: i
   ! 
   IF ( .NOT. obj%lwrite ) RETURN
   ! 
   CALL xml_NewElement( xp, TRIM(obj%tagname) ) 
   !
   CALL xml_NewElement( xp, "pot_prefactor") 
   CALL xml_addCharacters( xp, obj%pot_prefactor, fmt='s16') 
   CALL xml_endElement(xp, "pot_prefactor") 
   ! 
   CALL xml_NewElement( xp, "gate_zpos") 
   CALL xml_addCharacters( xp, obj%gate_zpos, fmt = 's16') 
   CALL xml_EndElement(xp, "gate_zpos") 
   ! 
   CALL xml_NewElement(xp, "gate_gate_term" )
   CALL xml_AddCharacters(xp, obj%gate_gate_term, fmt = 's16')
   CALL xml_endElement(xp, "gate_gate_term") 
   ! 
   CALL xml_NewElement(xp, "gatefieldEnergy" )
   CALL xml_AddCharacters(xp, obj%gatefieldEnergy, fmt = 's16')
   CALL xml_endElement(xp, "gatefieldEnergy") 
   ! 
   CALL xml_EndElement(xp, TRIM(obj%tagname))
END SUBROUTINE qes_write_gateInfo 


SUBROUTINE qes_init_gateInfo(obj, tagname, pot_prefactor, gate_zpos, gate_gate_term, gatefieldEnergy)
   IMPLICIT NONE
   TYPE(gateInfo_type),INTENT(INOUT)    :: obj
   CHARACTER(LEN=*),INTENT(IN)          :: tagname
   REAL(DP),INTENT(IN)                  :: pot_prefactor
   REAL(DP),INTENT(IN)                  :: gate_zpos
   REAL(DP),INTENT(IN)                  :: gate_gate_term
   REAL(DP),INTENT(IN)                  :: gatefieldEnergy
   !
   obj%tagname = TRIM(tagname)
   obj%lwrite = .TRUE.
   obj%lread  = .TRUE.
   obj%pot_prefactor = pot_prefactor 
   obj%gate_zpos     = gate_zpos
   obj%gate_gate_term = gate_gate_term
   obj%gatefieldEnergy = gatefieldEnergy
END SUBROUTINE qes_init_gateInfo


SUBROUTINE qes_reset_gateInfo(obj) 
   IMPLICIT NONE
   TYPE(gateInfo_type),INTENT(OUT) :: obj
   obj%tagname= ""
   obj%lwrite = .FALSE.
   obj%lread  = .FALSE.
END SUBROUTINE qes_reset_gateInfo


SUBROUTINE qes_write_outputElectricField(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(outputElectricField_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      IF(obj%BerryPhase_ispresent) CALL qes_write_BerryPhaseOutput(xp, obj%BerryPhase)
      !
      IF(obj%finiteElectricFieldInfo_ispresent) &
         CALL qes_write_finiteFieldOut(xp, obj%finiteElectricFieldInfo)
      !
      IF(obj%dipoleInfo_ispresent) CALL qes_write_dipoleOutput(xp, obj%dipoleInfo)
         
      ! 
      IF (obj%gateInfo_ispresent) CALL qes_write_gateInfo(xp, obj%gateInfo) 
      
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_outputElectricField

SUBROUTINE qes_init_outputElectricField(obj, tagname, BerryPhase, finiteElectricFieldInfo,& 
                                        dipoleInfo, gateInfo)
   IMPLICIT NONE

   TYPE(outputElectricField_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   TYPE(BerryPhaseOutput_type),OPTIONAL   :: BerryPhase
   TYPE(finiteFieldOut_type),OPTIONAL     :: finiteElectricFieldInfo
   TYPE(dipoleOutput_type),OPTIONAL       :: dipoleInfo
   TYPE(gateInfo_type),OPTIONAL           :: gateInfo

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%BerryPhase_ispresent = PRESENT(BerryPhase)
   IF(obj%BerryPhase_ispresent) THEN
      obj%BerryPhase = BerryPhase
   ENDIF
   obj%finiteElectricFieldInfo_ispresent = PRESENT(finiteElectricFieldInfo)
   IF(obj%finiteElectricFieldInfo_ispresent) THEN
      obj%finiteElectricFieldInfo = finiteElectricFieldInfo
   ENDIF
   obj%dipoleInfo_ispresent = PRESENT(dipoleInfo)
   IF(obj%dipoleInfo_ispresent) THEN
      obj%dipoleInfo = dipoleInfo
   ENDIF
   obj%gateInfo_ispresent = PRESENT(gateInfo) 
   IF ( obj%gateInfo_ispresent) obj%gateInfo = gateInfo
END SUBROUTINE qes_init_outputElectricField

SUBROUTINE qes_reset_outputElectricField(obj)
   IMPLICIT NONE
   TYPE(outputElectricField_type) :: obj
   INTEGER  :: i

   obj%tagname = ""
   obj%lwrite  = .FALSE.

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

SUBROUTINE qes_write_outputPBC(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(outputPBC_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL xml_NewElement(xp, 'assume_isolated' )
         CALL xml_addCharacters(xp,  TRIM(obj%assume_isolated))
      CALL xml_EndElement(xp, 'assume_isolated')
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_outputPBC

SUBROUTINE qes_init_outputPBC(obj, tagname, assume_isolated)
   IMPLICIT NONE

   TYPE(outputPBC_type) :: obj
   CHARACTER(len=*) :: tagname
   INTEGER  :: i
   CHARACTER(len=*) :: assume_isolated

   obj%tagname = TRIM(tagname)
   obj%lwrite   = .TRUE.
   obj%lread    = .TRUE.
   obj%assume_isolated = assume_isolated
END SUBROUTINE qes_init_outputPBC

SUBROUTINE qes_reset_outputPBC(obj)
   IMPLICIT NONE
   TYPE(outputPBC_type) :: obj
   INTEGER  :: i

   obj%tagname = ""
   obj%lwrite  =.FALSE.

END SUBROUTINE qes_reset_outputPBC

SUBROUTINE qes_write_output(xp, obj)
   IMPLICIT NONE

   TYPE(xmlf_t)  :: xp
   TYPE(output_type) :: obj
   !
   INTEGER  :: i

   IF (.NOT. obj%lwrite) RETURN
   !

   CALL xml_NewElement(xp, TRIM(obj%tagname))
      !
      CALL qes_write_convergence_info(xp, obj%convergence_info)
      !
      CALL qes_write_algorithmic_info(xp, obj%algorithmic_info)
      !
      CALL qes_write_atomic_species(xp, obj%atomic_species)
      !
      CALL qes_write_atomic_structure(xp, obj%atomic_structure)
      !
      IF(obj%symmetries_ispresent) THEN
         CALL qes_write_symmetries(xp, obj%symmetries)
         !
      ENDIF
      !
      CALL qes_write_basis_set(xp, obj%basis_set)
      !
      CALL qes_write_dft(xp, obj%dft)
      !Tsoh
      IF(obj%boundary_conditions_ispresent) THEN
         CALL qes_write_outputPBC(xp,obj%boundary_conditions)
      ENDIF
      !
      CALL qes_write_magnetization(xp, obj%magnetization)
      !
      CALL qes_write_total_energy(xp, obj%total_energy)
      !
      CALL qes_write_band_structure(xp, obj%band_structure)
      !
      IF(obj%forces_ispresent) THEN
         CALL qes_write_matrix(xp, obj%forces)
         !
      ENDIF
      !
      IF(obj%stress_ispresent) THEN
         CALL qes_write_matrix(xp, obj%stress)
         !
      ENDIF
      !
      IF(obj%electric_field_ispresent) THEN
         CALL qes_write_outputElectricField(xp, obj%electric_field)
         !
      ENDIF
      !
      IF(obj%FCP_force_ispresent) THEN
         CALL xml_NewElement(xp, 'FCP_force')
            CALL xml_addCharacters(xp, obj%FCP_force, fmt = 's16')
         CALL xml_EndElement(xp, 'FCP_force')
      ENDIF
      !
      IF(obj%FCP_tot_charge_ispresent) THEN
         CALL xml_NewElement(xp, 'FCP_tot_charge')
            CALL xml_addCharacters(xp, obj%FCP_tot_charge, fmt = 's16')
         CALL xml_EndElement(xp, 'FCP_tot_charge')
      ENDIF
      !
   CALL xml_EndElement(xp, TRIM(obj%tagname))
   !
END SUBROUTINE qes_write_output

SUBROUTINE qes_init_output(obj, tagname, convergence_info, algorithmic_info, &
                              atomic_species, atomic_structure, symmetries_ispresent, &
                              symmetries, basis_set, dft, magnetization, total_energy, &
                              band_structure, forces_ispresent, forces, stress_ispresent, &
                              stress, electric_field_ispresent, electric_field, &
                              FCP_force_ispresent, FCP_force, FCP_tot_charge_ispresent, &
                              FCP_tot_charge,boundary_conditions_ispresent, boundary_conditions)
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
   LOGICAL  :: boundary_conditions_ispresent
   TYPE(outputPBC_type) :: boundary_conditions
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
   IF(obj%boundary_conditions_ispresent) THEN
      obj%boundary_conditions=boundary_conditions
   ENDIF
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
   obj%lwrite  = .FALSE.

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
   IF(obj%boundary_conditions_ispresent) THEN
      CALL qes_reset_outputPBC(obj%boundary_conditions)
      obj%boundary_conditions_ispresent=.FALSE.
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


END MODULE qes_libs_module
