!
! Copyright (C) 2001-2004 Carlo Cavazzoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------

        INTEGER FUNCTION ig_local( ig, ig_l2g, sortedig_l2g, ng )
!
! This function computes the local index of the G vector whose
! global index is ig. If the G vector is not local to the current
! processor, then the function returns -1
!
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: ig
          INTEGER, INTENT(IN) :: ng
          INTEGER, INTENT(IN) :: ig_l2g( ng ), sortedig_l2g( ng )
          INTEGER :: lb, ub, i

          lb = 1   ! initialize search interval lower bound
          ub = ng  ! initialize search interval upper bound

          IF( ig < ig_l2g( sortedig_l2g(lb) ) .OR. ig > ig_l2g( sortedig_l2g(ub) ) )THEN
            ig_local = -1
            RETURN
          END IF

          BINARY_SEARCH: DO
            i = lb + (ub - lb)/2
            IF( ig >= ig_l2g( sortedig_l2g(i) ) )THEN
              lb = i
            ELSE IF( ig < ig_l2g( sortedig_l2g(i) ) )THEN
              ub = i
            ELSE
              lb = ub
            END IF
            IF( lb >= (ub-1) ) EXIT BINARY_SEARCH
          END DO BINARY_SEARCH

          IF( .NOT. ( (lb==ub) .OR. (lb==(ub-1)) ) )THEN
            CALL errore(' ig_local ',' algorithmic error ', 5)
          END IF
          IF( ig == ig_l2g( sortedig_l2g(lb) ) )THEN
            ig_local = sortedig_l2g(lb)
          ELSE IF( ig == ig_l2g( sortedig_l2g(ub) ) )THEN
            ig_local = sortedig_l2g(ub)
          ELSE
            ig_local = -1
          END IF

          RETURN
        END FUNCTION ig_local

