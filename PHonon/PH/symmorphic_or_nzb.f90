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
  USE symm_base,    ONLY : ft
  USE lr_symm_base, ONLY : gi, nsymq

  IMPLICIT NONE
  LOGICAL :: is_symmorphic, result_sym
  INTEGER :: isym, jsym
  REAL(DP) :: ft_(3,nsymq)

  is_symmorphic=.NOT.(ANY( ABS(ft(:,1:nsymq)) > 1.0d-8 ) )
  IF (is_symmorphic) THEN
     symmorphic_or_nzb=.TRUE.
     RETURN
  ELSE
     result_sym=.TRUE.
     ft_(:,1:nsymq) = ft(:,1:nsymq)
     CALL cryst_to_cart(nsymq, ft_, at, 1)

     DO isym=1,nsymq
        DO jsym=1,nsymq
           result_sym=( result_sym.AND.(ABS( gi(1,isym)*ft_(1,jsym) +  &
                                             gi(2,isym)*ft_(2,jsym) +  &
                                             gi(3,isym)*ft_(3,jsym) ) < 1.D-8) )
        END DO
     END DO
     symmorphic_or_nzb=result_sym
  END IF
  RETURN
END FUNCTION symmorphic_or_nzb
