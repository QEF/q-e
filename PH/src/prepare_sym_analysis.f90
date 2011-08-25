!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
  SUBROUTINE prepare_sym_analysis(nsym,sr,t_rev,magnetic_sym)

  USE kinds,    ONLY : DP
  USE rap_point_group,  ONLY : code_group, nclass, nelem, elem, which_irr,  &
                               char_mat, name_rap, gname, name_class, ir_ram
  USE rap_point_group_is, ONLY : code_group_is, gname_is

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nsym
  REAL(DP), INTENT(IN) :: sr(3,3,nsym)
  INTEGER, INTENT(IN) :: t_rev(nsym)
  LOGICAL, INTENT(IN) :: magnetic_sym

  INTEGER :: nsym_is, isym
  REAL(DP) :: sr_is(3,3,48)
!
!  Find the group name and sets its irreducible representation in the
!  rap_point_group module variables
!
  CALL find_group(nsym,sr,gname,code_group)
  CALL set_irr_rap(code_group,nclass,char_mat,name_rap,name_class,ir_ram)
  CALL divide_class(code_group,nsym,sr,nclass,nelem,elem,which_irr)
!
!  If some symmetry needs the time reversal check which group is formed
!  by the operations that do not need time reversal
!
  IF (magnetic_sym) THEN
     nsym_is=0
     DO isym=1,nsym
        IF (t_rev(isym)==0) THEN
           nsym_is=nsym_is+1
           sr_is(:,:,nsym_is) = sr(:,:,isym)
        ENDIF
     ENDDO
     CALL find_group(nsym_is,sr_is,gname_is,code_group_is)
  ENDIF

  RETURN
  END SUBROUTINE prepare_sym_analysis
