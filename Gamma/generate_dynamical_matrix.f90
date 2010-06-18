!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE generate_dynamical_matrix   &
     (nat, nsym, s, invs, irt, at, bg, n_diff_sites, equiv_atoms, &
     has_equivalent, dyn)
  !-----------------------------------------------------------------------
  !
  !  generate the complete dynamical matrix from independent modes only
  !  Input: dyn = irreducible dyn.mat.  Output: dyn = complete dyn.mat.
  !
  USE kinds, ONLY : DP
  USE symme, ONLY : crys_to_cart, cart_to_crys
  IMPLICIT NONE
  INTEGER :: nat, nsym, n_diff_sites, irt(48,nat), invs(48), &
       equiv_atoms(nat,nat), s(3,3,48),  has_equivalent(nat)
  real(DP) :: dyn(3*nat,3*nat), at(3,3), bg(3,3)
  !
  INTEGER :: isym, na, nb, ni, nj, sni, snj, smu_i, smu_j,  &
       i, j, k, l, mu_k, mu_l
  real(DP), ALLOCATABLE :: irreducible_dyn(:,:)
  real(DP) :: work(3,3)
  LOGICAL :: no_equivalent_atoms
  INTEGER, ALLOCATABLE ::done(:,:)
  !
  no_equivalent_atoms=.true.
  DO na = 1,nat
     no_equivalent_atoms = no_equivalent_atoms .and. has_equivalent(na)==0
  ENDDO
  IF (no_equivalent_atoms) RETURN
  !
  ALLOCATE  ( irreducible_dyn( 3*nat, 3*nat))
  CALL dcopy(3*nat*3*nat,dyn,1,irreducible_dyn,1)
  !
  DO na = 1,nat
     IF (has_equivalent(na)==0 ) THEN
        DO nb = 1,nat
           DO i = 1,3
              DO j = 1,3
                 work(i,j) = irreducible_dyn(3*(na-1)+i,3*(nb-1)+j)
              ENDDO
           ENDDO
           !
           !  transform to crystal axis
           !
           CALL cart_to_crys ( work )
           DO i = 1,3
              DO j = 1,3
                 irreducible_dyn(3*(na-1)+i,3*(nb-1)+j) = work(i,j)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  !
  ALLOCATE  (done( 3*nat, 3*nat))
  DO smu_i = 1,3*nat
     DO smu_j = 1,3*nat
        dyn(smu_i,smu_j) = 0.0d0
        done(smu_i,smu_j)= 0
     ENDDO
  ENDDO
  !
  DO isym = 1,nsym
     DO na = 1,n_diff_sites
        ni = equiv_atoms(na,1)
        sni = irt(isym,ni)
        DO i = 1,3
           smu_i = 3*(sni-1)+i
           DO nj = 1,nat
              snj = irt(isym,nj)
              DO j = 1,3
                 smu_j = 3*(snj-1)+j
                 IF (done(smu_i,smu_j)==0) THEN
                    DO k = 1,3
                       mu_k = 3*(ni-1)+k
                       DO l = 1,3
                          mu_l = 3*(nj-1)+l
                          dyn(smu_i,smu_j) = dyn(smu_i,smu_j) + &
                               s(i,k,invs(isym)) * s(j,l,invs(isym)) * &
                               irreducible_dyn(mu_k,mu_l)
                          !  rotation matrices are S^-1
                       ENDDO
                    ENDDO
                    done(smu_i,smu_j)=1
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  DEALLOCATE(done)
  DEALLOCATE(irreducible_dyn)
  !
  DO na = 1,nat
     DO nb = 1,nat
        DO i = 1,3
           DO j = 1,3
              work(i,j) = dyn(3*(na-1)+i,3*(nb-1)+j)
           ENDDO
        ENDDO
        !  back to cartesian axes
        CALL crys_to_cart ( work )
        DO i = 1,3
           DO j = 1,3
              dyn(3*(na-1)+i,3*(nb-1)+j) = work(i,j)
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  RETURN
END SUBROUTINE generate_dynamical_matrix
