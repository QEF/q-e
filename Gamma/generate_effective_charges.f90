!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE generate_effective_charges (nat, nsym, s, invs, irt, at, bg, &
     n_diff_sites, equiv_atoms, has_equivalent, zstar)
  !-----------------------------------------------------------------------
  !
  ! generate all effective charges
  !
  USE kinds, ONLY : DP
  USE symme, ONLY : crys_to_cart, cart_to_crys
  IMPLICIT NONE
  INTEGER :: nat, nsym, n_diff_sites, irt(48,nat), equiv_atoms(nat,nat),&
       s(3,3,48), has_equivalent(nat), invs(48)
  INTEGER :: isym, na, ni, nj, sni, i, j, k, l
  real(DP) :: zstar(3,3,nat), at(3,3), bg(3,3)
  LOGICAL :: done(nat), no_equivalent_atoms
  !
  no_equivalent_atoms=.true.
  DO na = 1,nat
     no_equivalent_atoms = no_equivalent_atoms .and. has_equivalent(na)==0
  ENDDO
  IF (no_equivalent_atoms) RETURN
  ! transform to crystal axis
  DO na = 1,nat
     IF (has_equivalent(na)==0 ) THEN
        CALL cart_to_crys ( zstar(:,:,na) )
        done(na)=.true.
     ELSE
        zstar(:,:,na) = 0.d0
        done(na)=.false.
     ENDIF
  ENDDO
  !
  DO isym = 1,nsym
     DO na = 1,n_diff_sites
        ni = equiv_atoms(na,1)
        sni = irt(isym,ni)
        IF ( .not.done(sni) ) THEN
           DO i = 1,3
              DO j = 1,3
                 DO k = 1,3
                    DO l = 1,3
                       zstar(i,j,sni) =  zstar(i,j,sni) +  &
                            s(i,k,invs(isym))*s(j,l,invs(isym))*zstar(k,l,ni)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           done(sni)=.true.
        ENDIF
     ENDDO
  ENDDO
  ! back to cartesian axis
  DO na = 1,nat
     CALL crys_to_cart ( zstar(:,:,na) )
  ENDDO
  !
  RETURN
END SUBROUTINE generate_effective_charges
