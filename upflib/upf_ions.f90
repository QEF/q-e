!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE upf_ions
  !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: n_atom_wfc
  
CONTAINS
  !
  !----------------------------------------------------------------------------
  FUNCTION n_atom_wfc( nat, ityp, noncolin )
  !----------------------------------------------------------------------------
  !
  ! ... Find number of starting atomic orbitals
  !
  USE uspp_param,   ONLY : upf
  !
  IMPLICIT NONE
  !
  INTEGER,            INTENT(IN) :: nat, ityp(nat)
  LOGICAL, OPTIONAL,  INTENT(IN) :: noncolin
  !
  INTEGER  :: n_atom_wfc
  !
  INTEGER  :: na, nt, n
  LOGICAL  :: non_col
  !
  !
  non_col = .FALSE.
  IF ( PRESENT (noncolin) ) non_col=noncolin
  n_atom_wfc = 0
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     DO n = 1, upf(nt)%nwfc
        !
        IF ( upf(nt)%oc(n) >= 0.D0 ) THEN
           !
           IF ( non_col ) THEN
              !
              IF ( upf(nt)%has_so ) THEN
                 !
                 n_atom_wfc = n_atom_wfc + 2 * upf(nt)%lchi(n)
                 !
                 IF ( ABS( upf(nt)%jchi(n)-upf(nt)%lchi(n) - 0.5D0 ) < 1.D-6 ) &
                    n_atom_wfc = n_atom_wfc + 2
                 !
              ELSE
                 !
                 n_atom_wfc = n_atom_wfc + 2 * ( 2 * upf(nt)%lchi(n) + 1 )
                 !
              END IF
              !
           ELSE
              !
              n_atom_wfc = n_atom_wfc + 2 * upf(nt)%lchi(n) + 1
              !
           END IF
        END IF
     END DO
  END DO
  !
  RETURN
  !
END FUNCTION n_atom_wfc


END MODULE upf_ions

