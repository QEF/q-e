!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE kcw_R_points
  !---------------------------------------------------------------------
  !
  !! This routine generates the R-points grid. Every R point
  !! corresponds to the position of primitive cell in a virtual
  !! supercell. R=0 is the origin and it corresponds to the 
  !! real primitive cell which all virtual cells are generated from.
  !
  USE cell_base,            ONLY : at
  USE control_kcw,          ONLY : Rvect, mp1, mp2 ,mp3, irvect
  USE klist,                ONLY : nkstot
  USE lsda_mod,             ONLY : nspin
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j, k, icell, num_R, mptot
  ! Number of unit cells ( = number of q points)
  !
  num_R = nkstot/nspin
  mptot=mp1*mp2*mp3
  IF (num_R .ne. mptot) &
     CALL errore('kcw_R_points', ' Mismatch between num of kpoints and MP grid from input', num_R)
  !
  ALLOCATE (Rvect(3,num_R))
  ALLOCATE (iRvect(3,num_R))
  !
  WRITE(stdout,'(/,5X, "INFO: total number of primitive cell", i5)') num_R
  !
  IF ( nkstot == 1 ) THEN
    !
    Rvect(:,1) = 0.0d0
    irvect(:,1) = (/0,0,0/)
    !
  ELSE
    !
    ! "at" are in units of alat
    !
    icell = 0
    !
    DO i = 1, mp1
      DO j = 1, mp2
        DO k = 1, mp3
           !
           icell = icell + 1
           !
           Rvect(:,icell) = DBLE(i-1) * at(:,1) + &
                            DBLE(j-1) * at(:,2) + &
                            DBLE(k-1) * at(:,3)
           irvect(:,icell) = (/i-1,j-1,k-1/)
           !
        ENDDO
      ENDDO
    ENDDO
    !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE kcw_R_points

