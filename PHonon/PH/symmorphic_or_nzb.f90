!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
LOGICAL FUNCTION symmorphic_or_nzb()
!
! This function returns .true. if the small group of the current q is symmorphic
! or if the q point is not at zone border or if all the phase factors
! e^{i G f} are equal to one. Presently the routine that
! finds the mode symmetry works only when this function is .true..
!
  USE kinds,        ONLY : DP
  USE cell_base,    ONLY : at
  USE fft_base,     ONLY : dfftp
  USE symm_base,    ONLY : ftau
  USE lr_symm_base, ONLY : gi, nsymq

  IMPLICIT NONE
  LOGICAL :: is_symmorphic, result_sym
  INTEGER :: isym, jsym
  REAL(DP) :: ft(3,nsymq)

  is_symmorphic=.NOT.(ANY(ftau(:,1:nsymq) /= 0))
  IF (is_symmorphic) THEN
     symmorphic_or_nzb=.TRUE.
     RETURN
  ELSE
     result_sym=.TRUE.
     DO isym = 1, nsymq
        ft(1,isym) = DBLE(ftau(1,isym)) / DBLE(dfftp%nr1)
        ft(2,isym) = DBLE(ftau(2,isym)) / DBLE(dfftp%nr2)
        ft(3,isym) = DBLE(ftau(3,isym)) / DBLE(dfftp%nr3)
     END DO
     CALL cryst_to_cart(nsymq, ft, at, 1)

     DO isym=1,nsymq
        DO jsym=1,nsymq
           result_sym=( result_sym.AND.(ABS( gi(1,isym)*ft(1,jsym) +  &
                                             gi(2,isym)*ft(2,jsym) +  &
                                             gi(3,isym)*ft(3,jsym) ) < 1.D-8) )
        END DO
     END DO
     symmorphic_or_nzb=result_sym
  END IF
  RETURN
END FUNCTION symmorphic_or_nzb
