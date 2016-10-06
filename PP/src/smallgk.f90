!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE smallgk (xk, at, bg, s, ftau, t_rev, sname, nsym, sk, &
                    ftauk, gk, t_revk, invsk, snamek, nsymk)
!-----------------------------------------------------------------------
!
! This routine selects, among the symmetry matrices of the point group
! of a crystal, the symmetry operations which leave k unchanged.
!
!
USE kinds, ONLY : DP
IMPLICIT NONE
REAL(DP), PARAMETER :: accep=1.d-5
CHARACTER(len=45) :: snamek(48), sname(48)

REAL(DP) :: bg (3, 3), at (3, 3), xk (3)
! input: the reciprocal lattice vectors
! input: the direct lattice vectors
! input: the k point of the crystal

INTEGER :: s (3, 3, 48), ftau(3,48), t_rev(48), nsym, sk (3, 3, 48), &
           ftauk(3,48), t_revk(48), gk(3,48), invsk(48), nsymk
! input: the symmetry matrices
! input: fractional translation associated to each rotation
! input: possible time reversal associated to the rotation
! input: the inverse of each symmetry operation
! input: dimension of the point group
! output: the symmetry matrices of the small point group of k
! output: the fract. trans. associated to the operations of the small group of k
! output: the time reversal associated to the operations of the small group of k
! output: the G vector which connects k and the rotated k.
! output: the inverse of each operation or the small point group of k.

  REAL(DP) :: ak (3), rak (3), zero (3)
  ! k vector in crystal basis
  ! the rotated of the k vector
  ! the zero vector

  INTEGER :: isym, jsym, ipol, jpol
  ! counter on symmetry operations
  ! counter on symmetry operations
  ! counter on polarizations
  ! counter on polarizations
  INTEGER :: ss(3,3)

  LOGICAL :: eqvect, found
  ! logical function, check if two vectors are equal
  ! logical variable to check the inverse
  !
  !  Set to zero some variables and transform xq to the crystal basis
  !
  zero = 0.d0
  ak = xk
  CALL cryst_to_cart (1, ak, at, - 1)
  !
  !   test all symmetries to see if the operation S sends k in k+G ...
  !
  nsymk = 0
  DO isym = 1, nsym
     rak = 0.d0
     DO ipol = 1, 3
        DO jpol = 1, 3
           rak (ipol) = rak (ipol) + dble (s (ipol, jpol, isym) ) * &
                ak (jpol)
        ENDDO
     ENDDO
     IF ((t_rev(isym)==0 .and. eqvect(rak, ak, zero,accep)) .or. &
         (t_rev(isym)==1 .and. eqvect(rak, -ak, zero,accep)) ) THEN
        nsymk=nsymk+1
        sk(:,:,nsymk)=s(:,:,isym)
        ftauk(:,nsymk)=ftau(:,isym)
        snamek(nsymk)=sname(isym)
        t_revk(nsymk)=t_rev(isym)
        IF (t_rev(isym)==0) THEN
           gk(:,nsymk)=nint(rak(:)-ak(:))
        ELSEIF (t_rev(isym)==1) THEN
           gk(:,nsymk)=nint(rak(:)+ak(:))
        ELSE
           CALL errore('smallgk','wrong t_rev',1)
        ENDIF
     ENDIF
  ENDDO
!
!  Find the inverse of each element
!
  DO isym = 1, nsymk
     found = .FALSE.
     DO jsym = 1, nsymk
        !
        ss = MATMUL (sk(:,:,jsym),sk(:,:,isym))
        ! s(:,:,1) is the identity
        IF ( ALL ( sk(:,:,1) == ss(:,:) ) ) THEN
           invsk (isym) = jsym
           found = .TRUE.
        ENDIF
     ENDDO
     IF ( .NOT.found) CALL errore ('smallgk', ' Not a group', 1)
  ENDDO
  !
  RETURN
END SUBROUTINE smallgk

