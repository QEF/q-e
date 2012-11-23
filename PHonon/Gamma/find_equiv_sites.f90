!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE find_equiv_sites (nat,nsym,irt,has_equivalent,     &
     n_diff_sites,n_equiv_atoms,equiv_atoms)
  !
  IMPLICIT NONE
  INTEGER :: nat, nsym, na, nb, ns, n_diff_sites, irt(48,nat),    &
       equiv_atoms(nat,nat), n_equiv_atoms(nat), has_equivalent(nat)
  !
  n_diff_sites = 0
  DO na = 1,nat
     has_equivalent(na) = 0
  ENDDO
  !
  DO na = 1,nat
     IF (has_equivalent(na)==0) THEN
        n_diff_sites = n_diff_sites + 1
        n_equiv_atoms (n_diff_sites) =  1
        equiv_atoms(n_diff_sites,1) = na
        !
        DO nb = na+1,nat
           DO ns = 1, nsym
              IF ( irt(ns,nb) == na) THEN
                 has_equivalent(nb) = 1
                 n_equiv_atoms (n_diff_sites) =  &
                      n_equiv_atoms (n_diff_sites) + 1
                 equiv_atoms(n_diff_sites, &
                      n_equiv_atoms(n_diff_sites)) = nb
                 GOTO 10
              ENDIF
           ENDDO
10         CONTINUE
        ENDDO
     ENDIF
  ENDDO
  !
  RETURN
END SUBROUTINE find_equiv_sites
