!
! Copyright (C) 2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
  SUBROUTINE rotate_pattern_add(nat, u, dyn, dynwrk)
  !--------------------------------------------------------------------
  !! This routine rotates the dynamical matrix dynwork written
  !! in cartesian basis to the basis of the patterns u and adds it to
  !! the dynamical matrix dyn that is supposed to be in the basis of the
  !! patterns.
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  COMPLEX(DP), INTENT(IN) :: u(3*nat, 3*nat)
  COMPLEX(DP), INTENT(INOUT) :: dyn(3*nat, 3*nat)
  COMPLEX(DP), INTENT(IN) :: dynwrk(3*nat, 3*nat)
  !
  COMPLEX(DP), ALLOCATABLE :: work(:, :)
  ALLOCATE(work(3*nat, 3*nat))
  !
  !INTEGER :: nu_i, nu_j, na_icart, na_jcart
  !
  !DO nu_i = 1, 3 * nat
  !   DO nu_j = 1, 3 * nat
  !      work = (0.0d0, 0.0d0)
  !      DO na_jcart = 1, 3 * nat
  !         DO na_icart = 1, 3 * nat
  !            work = work + CONJG (u (na_icart, nu_i) ) * &
  !                          dynwrk (na_icart, na_jcart) * &
  !                          u (na_jcart, nu_j)
  !         ENDDO
  !      ENDDO
  !      dyn (nu_i, nu_j) = dyn (nu_i, nu_j) + work
  !   ENDDO
  !ENDDO
  !
  CALL zgemm('N', 'N', 3*nat, 3*nat, 3*nat, (1.0d0, 0.d0), dynwrk, 3*nat, u,    3*nat, (0.d0, 0.d0), work, 3*nat)
  CALL zgemm('C', 'N', 3*nat, 3*nat, 3*nat, (1.0d0, 0.d0), u,      3*nat, work, 3*nat, (1.d0, 0.d0), dyn,  3*nat)
  !
  DEALLOCATE(work)
  RETURN 
  END SUBROUTINE rotate_pattern_add
!
!----------------------------------------------------------------------
  SUBROUTINE dyn_pattern_to_cart(nat, u, dyn, phi)
  !---------------------------------------------------------------------
  !! This routine transforms the dynamical matrix dyn, written in the basis
  !! of the pattern, in the dynamical matrix phi, in the cartesian basis.
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  COMPLEX(DP), INTENT(IN) :: u(3*nat, 3*nat)
  COMPLEX(DP), INTENT(IN) :: dyn(3*nat, 3*nat)
  COMPLEX(DP), INTENT(OUT) :: phi(3, 3, nat, nat)
  !
  COMPLEX(DP), ALLOCATABLE :: work1(:, :), work2(:, :)
  ALLOCATE(work1(3*nat, 3*nat), work2(3*nat, 3*nat))

  !INTEGER :: i, j, icart, jcart, na, nb, mu, nu
  !
  !DO i = 1, 3 * nat
  !   na = (i - 1) / 3 + 1
  !   icart = i - 3 * (na - 1)
  !   DO j = 1, 3 * nat
  !      nb = (j - 1) / 3 + 1
  !      jcart = j - 3 * (nb - 1)
  !      work = (0.d0, 0.d0)
  !      DO mu = 1, 3 * nat
  !         DO nu = 1, 3 * nat
  !            work = work + u (i, mu) * dyn (mu, nu) * CONJG(u (j, nu) )
  !         ENDDO
  !      ENDDO
  !      phi (icart, jcart, na, nb) = work
  !   ENDDO
  !ENDDO
  CALL zgemm('N', 'C', 3*nat, 3*nat, 3*nat, (1.0d0, 0.d0), dyn, 3*nat, u,    3*nat, (0.d0, 0.d0), work1, 3*nat)
  CALL zgemm('N', 'N', 3*nat, 3*nat, 3*nat, (1.0d0, 0.d0), u,   3*nat, work1,3*nat, (0.d0, 0.d0), work2, 3*nat)
  !
  CALL scompact_dyn(nat, work2, phi)
  DEALLOCATE(work1, work2)
  !
  RETURN 
  END SUBROUTINE dyn_pattern_to_cart
!
!-----------------------------------------------------------------------
  SUBROUTINE compact_dyn(nat, dyn, phi)
  !-----------------------------------------------------------------------
  !! This routine writes the dynamical matrix from a 3,3,nat,nat array
  !! to a 3*nat,3*nat array.
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  COMPLEX(DP), INTENT(IN) :: phi(3,3,nat,nat)
  COMPLEX(DP), INTENT(OUT) :: dyn(3*nat, 3*nat)

  INTEGER :: na, nb, icart, jcart, imode, jmode

  DO na = 1, nat
     DO icart = 1, 3
        imode = 3 * ( na - 1 ) + icart
        DO nb = 1, nat
           DO jcart = 1, 3
              jmode = 3 * ( nb - 1 ) + jcart
              dyn (imode, jmode) = phi (icart, jcart, na, nb)
           END DO 
        END DO
     END DO
  END DO
  RETURN
  END SUBROUTINE compact_dyn
!
!-----------------------------------------------------------------------
  SUBROUTINE scompact_dyn(nat, dyn, phi)
  !-----------------------------------------------------------------------
  !! This routine writes the dynamical matrix from a 3*nat,3*nat array
  !! to a 3,3,nat,nat array.
  !

  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  COMPLEX(DP), INTENT(OUT) :: phi(3,3,nat,nat)
  COMPLEX(DP), INTENT(IN) :: dyn(3*nat, 3*nat)

  INTEGER :: na, nb, icart, jcart, imode, jmode

  DO na = 1, nat
     DO icart = 1, 3
        imode = 3 * ( na - 1 ) + icart
        DO nb = 1, nat
           DO jcart = 1, 3
              jmode = 3 * ( nb - 1 ) + jcart
              phi (icart, jcart, na, nb) = dyn (imode, jmode) 
           END DO 
        END DO
     END DO
  END DO
  RETURN
  END SUBROUTINE scompact_dyn
