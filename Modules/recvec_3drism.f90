!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE gvec_3drism
  !--------------------------------------------------------------------------
  !
  ! ... the reciprocal lattice vectors for 3D-RISM
  !
  USE cell_base,     ONLY : at
  USE constants,     ONLY : eps8
  USE control_flags, ONLY : gamma_only
  USE fft_types,     ONLY : fft_type_descriptor
  USE gvect,         ONLY : gg, g, mill, gstart
  USE kinds,         ONLY : DP
  USE mp,            ONLY : mp_max, mp_sum
  USE recvec_subs,   ONLY : ggens
  USE parallel_include
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define variables
  TYPE gvec_type
    !
    REAL(DP) :: ecut
    REAL(DP) :: gcutm
    !
    INTEGER  :: ngm
    INTEGER  :: ngmx
    INTEGER  :: ngm_g
    !
    REAL(DP), POINTER :: gg(:)
    REAL(DP), POINTER :: g (:, :)
    INTEGER           :: gstart
    !
    INTEGER,  POINTER :: ig1(:)
    INTEGER,  POINTER :: ig2(:)
    INTEGER,  POINTER :: ig3(:)
    !
    INTEGER           :: ngl
    REAL(DP), POINTER :: gl(:)
    INTEGER,  POINTER :: igtongl(:)
    !
  END TYPE gvec_type
  !
  ! ... public components
  PUBLIC :: gvec_type
  PUBLIC :: gvec_init_3drism
  PUBLIC :: deallocate_gvec_3drism
  PUBLIC :: ggen_3drism
  PUBLIC :: gshells_3drism
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE gvec_init_3drism(gvec, ngm_, comm)
    !
    ! ... Set local and global dimensions, allocate arrays
    !
    IMPLICIT NONE
    !
    TYPE(gvec_type),  INTENT(INOUT) :: gvec
    INTEGER,          INTENT(IN)    :: ngm_
    INTEGER,          INTENT(IN)    :: comm  ! communicator of the group on which G-vecs are distributed
    !
    gvec%ngm = ngm_
    !
    ! ... calculate maximum over all processors
    !
    gvec%ngmx = ngm_
    CALL mp_max(gvec%ngmx, comm)
    !
    ! ... calculate sum over all processors
    !
    gvec%ngm_g = ngm_
    CALL mp_sum(gvec%ngm_g, comm)
    !
  END SUBROUTINE gvec_init_3drism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE deallocate_gvec_3drism(gvec)
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(gvec_type),  INTENT(INOUT) :: gvec
    !
    IF (ASSOCIATED(gvec%gg))      DEALLOCATE(gvec%gg)
    IF (ASSOCIATED(gvec%g))       DEALLOCATE(gvec%g)
    !
    IF (ASSOCIATED(gvec%ig1))     DEALLOCATE(gvec%ig1)
    IF (ASSOCIATED(gvec%ig2))     DEALLOCATE(gvec%ig2)
    IF (ASSOCIATED(gvec%ig3))     DEALLOCATE(gvec%ig3)
    !
    IF (ASSOCIATED(gvec%gl))      DEALLOCATE(gvec%gl)
    IF (ASSOCIATED(gvec%igtongl)) DEALLOCATE(gvec%igtongl)
    !
  END SUBROUTINE deallocate_gvec_3drism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ggen_3drism(gvec, dfft)
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    TYPE(gvec_type),           INTENT(INOUT) :: gvec
    TYPE(fft_type_descriptor), INTENT(INOUT) :: dfft
    !
    INTEGER :: ig
    INTEGER :: n1, n2, n3
    !
    ! ... set G-vectors
    !
    CALL ggens(dfft, gamma_only, at, g, gg, mill, gvec%gcutm, gvec%ngm, gvec%g, gvec%gg)
    !
    gvec%gstart = gstart
    !
    ! ... set ig1, ig2, ig3, which are 3d-indexes of FFT-box for the ig-th G-vector.
    !
    ALLOCATE(gvec%ig1(gvec%ngm))
    ALLOCATE(gvec%ig2(gvec%ngm))
    ALLOCATE(gvec%ig3(gvec%ngm))
    !
    DO ig = 1, gvec%ngm
      !
      n1 = NINT(SUM(gvec%g(:, ig) * at(:, 1))) + 1
      gvec%ig1(ig) = n1 - 1
      !
      n2 = NINT(SUM(gvec%g(:, ig) * at(:, 2))) + 1
      gvec%ig2(ig) = n2 - 1
      !
      n3 = NINT(SUM(gvec%g(:, ig) * at(:, 3))) + 1
      gvec%ig3(ig) = n3 - 1
      !
      IF (n1 < 1) n1 = n1 + dfft%nr1
      IF (n1 < 1 .OR. dfft%nr1 < n1) CALL errore('ggen_3drism', 'incorrect ig1', ig)
      !
      IF (n2 < 1) n2 = n2 + dfft%nr2
      IF (n2 < 1 .OR. dfft%nr2 < n2) CALL errore('ggen_3drism', 'incorrect ig2', ig)
      !
      IF (n3 < 1) n3 = n3 + dfft%nr3
      IF (n3 < 1 .OR. dfft%nr3 < n3) CALL errore('ggen_3drism', 'incorrect ig3', ig)
      !
    END DO
    !    
  END SUBROUTINE ggen_3drism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE gshells_3drism(gvec)
    !--------------------------------------------------------------------------
    !
    ! ... calculate number of G shells: ngl, and the index ig = igtongl(ig)
    ! ... that gives the shell index for ig-th G-vector.
    !
    IMPLICIT NONE
    !
    TYPE(gvec_type), INTENT(INOUT) :: gvec
    !
    INTEGER :: ig
    INTEGER :: igl
    !
    ! ... deallocate memory, if needed
    !
    IF(ASSOCIATED(gvec%gl))      DEALLOCATE(gvec%gl)
    IF(ASSOCIATED(gvec%igtongl)) DEALLOCATE(gvec%igtongl)
    !
    ! ... G-vectors are grouped in shells with the same norm
    !
    ALLOCATE(gvec%igtongl(gvec%ngm))
    !
    gvec%ngl        = 1
    gvec%igtongl(1) = 1
    !
    DO ig = 2, gvec%ngm
      !
      IF (gvec%gg(ig) > gvec%gg(ig - 1) + eps8) THEN
        gvec%ngl = gvec%ngl + 1
      END IF
      !
      gvec%igtongl(ig) = gvec%ngl
      !
    END DO
    !
    ALLOCATE(gvec%gl(gvec%ngl))
    !
    igl = 1
    gvec%gl(1) = gvec%gg(1)
    !
    DO ig = 2, gvec%ngm
      !
      IF (gvec%gg(ig) > gvec%gg(ig - 1) + eps8) THEN
        igl = igl + 1
        gvec%gl(igl) = gvec%gg(ig)
      END IF
      !
    END DO
    !
    IF (igl /= gvec%ngl) CALL errore ('gshells_3drism', 'igl <> ngl', gvec%ngl)
    !
  END SUBROUTINE gshells_3drism
  !
END MODULE gvec_3drism
