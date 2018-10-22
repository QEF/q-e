!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE hp_R_points
  !---------------------------------------------------------------------
  !
  ! This routine generates the R-points grid. Every R point
  ! corresponds to the position of primitive cell in a virtual
  ! supercell. R=0 is the origin and it corresponds to the 
  ! real primitive cell from which all virtual cells are generated.
  ! R appear in the phase factor in Eq. (42) in Ref. [1]
  ! [1] Phys. Rev. B 98, 085127 (2018)
  !
  USE cell_base,    ONLY : at
  USE ldaU_hp,      ONLY : nqsh, Rvect, nq1, nq2, nq3
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j, k, ipol, icell
  !
  ! Number of unit cells ( = number of q points)
  !
  ALLOCATE (Rvect(3,nqsh))
  !
  IF ( nqsh == 1 ) THEN
    !
    Rvect(:,1) = 0.0d0
    !
  ELSE
    !
    ! "at" are in units of alat
    !
    icell = 0
    !
    DO i = 1, nq1
      DO j = 1, nq2
        DO k = 1, nq3
           !
           icell = icell + 1
           !
           Rvect(:,icell) = DBLE(i-1) * at(:,1) + &
                            DBLE(j-1) * at(:,2) + &
                            DBLE(k-1) * at(:,3)
           !
        ENDDO
      ENDDO
    ENDDO
    !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE hp_R_points
