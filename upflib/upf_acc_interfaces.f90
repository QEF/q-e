!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE upf_acc_interfaces
  !
  INTERFACE sph_bes_acc
    SUBROUTINE sph_bes( msh, r, q, l, jl )
      !$acc routine vector
      IMPLICIT NONE
      INTEGER  :: msh, l
      REAL(8) :: r(msh), q, jl(msh)
    END SUBROUTINE
  END INTERFACE
  !
  INTERFACE
    SUBROUTINE simpson( mesh, func, rab, asum )
      !$acc routine vector
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: rab(mesh), func(mesh)
      REAL(8), INTENT(OUT):: asum
    END SUBROUTINE
  END INTERFACE
  !
END MODULE
