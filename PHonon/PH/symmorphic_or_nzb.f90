!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
LOGICAL FUNCTION symmorphic_or_nzb()
!
! This function returns true if the small group of the current q is symmorphic
! or if the q point is not at zone border. Presently the routine that
! finds the mode symmetry works only when this function is .true..
!
  USE symm_base,    ONLY : ftau
  USE lr_symm_base, ONLY : gi, nsymq

  IMPLICIT NONE
  LOGICAL :: is_symmorphic, result_sym
  INTEGER :: isym

  is_symmorphic=.NOT.(ANY(ftau(:,1:nsymq) /= 0))
  IF (is_symmorphic) THEN
     symmorphic_or_nzb=.TRUE.
     RETURN
  ELSE
     result_sym=.TRUE.
     DO isym=1,nsymq
        result_sym=( result_sym .AND. (abs(gi(1,isym))<1.d-8).and.  &
                                      (abs(gi(2,isym))<1.d-8).and.  &
                                      (abs(gi(3,isym))<1.d-8) )
     END DO
     symmorphic_or_nzb=result_sym
  END IF
  RETURN
END FUNCTION symmorphic_or_nzb
