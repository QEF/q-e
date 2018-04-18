!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE g_psi_mod_gpum
  !
  ! ... These are the variables needed in g_psi
  !  
  USE kinds, only : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), DEVICE, ALLOCATABLE :: &
    h_diag_d (:,:),&   ! diagonal part of the Hamiltonian
    s_diag_d (:,:)     ! diagonal part of the overlap matrix
  
  LOGICAL :: h_diag_ood = .false.    ! used to flag out of date variables
  LOGICAL :: s_diag_ood = .false.
  LOGICAL :: h_diag_d_ood = .false.
  LOGICAL :: s_diag_d_ood = .false.
  !
  CONTAINS
  !
  SUBROUTINE using_h_diag(writing)
      USE g_psi_mod, ONLY : h_diag
      implicit none
      LOGICAL, INTENT(IN) :: writing
      !
      IF (h_diag_ood) THEN
          h_diag = h_diag_d
          h_diag_ood = .false.
      ENDIF
      IF (writing)    h_diag_d_ood = .true.
  END SUBROUTINE using_h_diag
  !
  SUBROUTINE using_h_diag_d(writing)
      USE g_psi_mod, ONLY : h_diag
      implicit none
      LOGICAL, INTENT(IN) :: writing
      !
      IF (.not. allocated(h_diag)) THEN
          IF (allocated(h_diag_d)) DEALLOCATE(h_diag_d)
          h_diag_d_ood = .false.
          RETURN
      END IF
      IF (h_diag_d_ood) THEN
          IF (.not. allocated(h_diag_d)) THEN
              ALLOCATE(h_diag_d, SOURCE=h_diag)
          ELSE
              h_diag_d = h_diag
          ENDIF
          h_diag_d_ood = .false.
      ENDIF
      IF (writing)    h_diag_ood = .true.
  END SUBROUTINE using_h_diag_d
  !
  SUBROUTINE using_s_diag(writing)
      USE g_psi_mod, ONLY : s_diag
      implicit none
      LOGICAL, INTENT(IN) :: writing
      !
      IF (s_diag_ood) THEN
          s_diag = s_diag_d
          s_diag_ood = .false.
      ENDIF
      IF (writing)    s_diag_d_ood = .true.
  END SUBROUTINE using_s_diag
  !
  SUBROUTINE using_s_diag_d(writing)
      USE g_psi_mod, ONLY : s_diag
      implicit none
      LOGICAL, INTENT(IN) :: writing
      !
      IF (.not. allocated(s_diag)) THEN
          IF (allocated(s_diag_d)) DEALLOCATE(s_diag_d)
          s_diag_d_ood = .false.
          RETURN
      END IF
      IF (s_diag_d_ood) THEN
          IF (.not. allocated(s_diag_d)) THEN
              ALLOCATE(s_diag_d, SOURCE=s_diag)
          ELSE
              s_diag_d = s_diag
          END IF
          s_diag_d_ood = .false.
      ENDIF
      IF (writing)    s_diag_ood = .true.
  END SUBROUTINE using_s_diag_d
  !
END MODULE g_psi_mod_gpum
