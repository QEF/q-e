!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE allocate_fft_3drism(dfft, gvec, ecutv, laue, mp_task)
  !---------------------------------------------------------------------------
  !
  ! ... initialize 3D-FFT for 3D-RISM
  !
  USE cell_base,     ONLY : tpiba2
  USE control_flags, ONLY : gamma_only
  USE fft_types,     ONLY : fft_type_descriptor
  USE gvec_3drism,   ONLY : gvec_type, ggen_3drism, gshells_3drism
  USE kinds,         ONLY : DP
  USE mp_rism,       ONLY : mp_rism_task
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor), INTENT(INOUT) :: dfft
  TYPE(gvec_type),           INTENT(INOUT) :: gvec
  REAL(DP),                  INTENT(IN)    :: ecutv
  LOGICAL,                   INTENT(IN)    :: laue
  TYPE(mp_rism_task),        INTENT(IN)    :: mp_task
  !
  gvec%ecut  = ecutv
  gvec%gcutm = ecutv / tpiba2
  !
  CALL data_structure_3drism(dfft, gvec, gamma_only, mp_task)
  !
  CALL ggen_3drism(gvec, dfft)
  !
  IF (.NOT. laue) THEN
    CALL gshells_3drism(gvec)
  END IF
  !
END SUBROUTINE allocate_fft_3drism
