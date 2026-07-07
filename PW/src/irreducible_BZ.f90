!
! Copyright (C) 2001-2011 Quantum ESPRESSO  group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE irreducible_BZ( nrot, s, nsym, minus_q, magnetic_sym, at, bg, &
                           npk, nks, xk, wk, t_rev )
   !-----------------------------------------------------------------------
   !! Given a set of k points, unfold it using symmetries of the Bravais lattice (nrot)
   !! and then reduce it using symmetries of the crystal (nsym < nrot).
   !! This subroutine can be used to start with the irreducible k points of the Bravais
   !! lattice and then expand it into the irreducible k points of the crystal.
   !!
   !! NOTE: The initial k points are all kept, even if there are redundant or symmetry-
   !! equivalent ones. This is to be consistent with previous versions of QE, and is needed
   !! e.g. in the case of calculating polarization where a string of k points along a given
   !! direction is needed.
   !
   USE kinds,   ONLY: DP
   !
   IMPLICIT NONE
   !
   INTEGER,  INTENT(IN) :: nrot
   !! order of the parent point group
   INTEGER,  INTENT(IN) :: nsym
   !! order of the subgroup
   INTEGER,  INTENT(IN) :: npk
   !! maximum number of special points
   INTEGER,  INTENT(IN) :: s(3,3,48)
   !! symmetry matrices, in crystal axis
   INTEGER,  INTENT(IN) :: t_rev(48)
   !! time reversal operation
   REAL(DP), INTENT(IN) :: at(3,3)
   !! basis vectors of the Bravais lattice
   REAL(DP), INTENT(IN) :: bg(3,3)
   !! basis vectors of the reciprocal lattice
   LOGICAL,  INTENT(IN) :: minus_q
   !! it is .TRUE. if symmetries q=-q+G are acceptable
   LOGICAL,  INTENT(IN) :: magnetic_sym
   !! magnetic_sym = noncolin .AND. domag
   INTEGER,  INTENT(INOUT) :: nks
   !! number of special points
   REAL(DP), INTENT(INOUT) :: xk(3,npk)
   !! special points
   REAL(DP), INTENT(INOUT) :: wk(npk)
   !! weights for special points
   !
   ! ... local variables
   !
   LOGICAL :: satm
   ! true if equivalent point found
   INTEGER :: table(48,48), invs(3,3,48)
   ! table: multiplication table of the group
   ! invs:  contains the inverse of each rotation
   INTEGER :: isym, jsym, nks0, jk, irot, jrot, ik
   ! nks0: used to save the initial number of k-points
   INTEGER :: count
   !! Number of k points equivalent to ik-th k point
   REAL(DP) :: xkg(3), xks(3), xkn(3), one, xk_new(3,npk), wk_new(npk), wk_sum
   ! coordinates of the k point in crystal axis
   ! coordinates of the rotated k point
   ! buffer which contains the weight of k points
   ! total weight of k-points
   REAL(DP) :: wk_for_ik
   !! used to sum the weights of k points equivalent to ik-th k point
   LOGICAL, ALLOCATABLE :: done(:)
   !! True if the k-point has been already processed
   INTEGER, ALLOCATABLE :: equivalent_with_ik(:)
   !! Index of k points equivalent to ik-th k point
   !
   ! ... We compute the multiplication table of the group.
   !
   CALL multable( nrot, s, table )
   !
   ! ... And we set the matrices of the inverse.
   !
   DO isym = 1, nrot
      DO jsym = 1, nrot
         IF (table(isym,jsym) == 1) invs(:,:,isym) = s(:,:,jsym)
      ENDDO
   ENDDO
   !
   wk_sum = SUM(wk(1:nks))
   nks0 = nks
   !
   ! First add all original k points to the new list. Do not reduce them by symmetry.
   !
   wk_new(1:nks) = wk(1:nks)
   xk_new(:, 1:nks) = xk(:, 1:nks)
   CALL cryst_to_cart(nks, xk_new, at, -1)  ! Convert from Cartesian to crystal axis
   !
   ! Now loop over original k points and symmetry of the parent group (nrot) to unfold
   ! the k points, and then reduce them using the symmetries of the subgroup (nsym)
   !
   DO jk = 1, nks0
      !
      ! ... The k point is first computed in crystal axis
      !
      ! xkg are the components of xk in the crystal base
      xkg(:) = at(1,:) * xk(1,jk) + &
               at(2,:) * xk(2,jk) + &
               at(3,:) * xk(3,jk)
      !
      ! Skip irot = 1 which is the identity (already all states are kept).
      !
      DO irot = 2, nrot
         !
         ! ... Then it is rotated with each symmetry of the global group.
         !
         xks(:) = invs(:,1,irot) * xkg(1) + &
                  invs(:,2,irot) * xkg(2) + &
                  invs(:,3,irot) * xkg(3)
         IF (magnetic_sym .AND. (t_rev(irot) == 1)) xks = -xks
         !
         !  ... Now check if there is an operation of the subgroup that
         !      makes xks equivalent to some other already found k point
         !
         DO jrot = 1, nsym
            xkn(:) = invs(:,1,jrot) * xks(1) + &
                     invs(:,2,jrot) * xks(2) + &
                     invs(:,3,jrot) * xks(3)
            IF (magnetic_sym .AND. (t_rev(jrot) == 1)) xkn = -xkn
            !
            DO ik = 1, nks
               satm = are_k_equivalent(xk_new(:,ik), xkn)
               IF (minus_q) satm = satm .OR. are_k_equivalent(xk_new(:,ik), -xkn)
               !
               IF (satm) THEN
                  wk_new(ik) = wk_new(ik) + wk(jk)
                  GOTO 100
               ENDIF
            ENDDO
         ENDDO
         nks = nks+1
         IF (nks > npk) CALL errore('irreducible_BZ', 'too many k points (step 2)', 1)
         xk_new(:,nks) = xks
         wk_new(nks) = wk(jk)
100      CONTINUE
      ENDDO ! irot
   ENDDO ! jk
   !
   ! Average weights of equivalent k points
   !
   ALLOCATE(equivalent_with_ik(nks))
   ALLOCATE(done(nks))
   done(:) = .FALSE.
   equivalent_with_ik(:) = -1
   !
   DO ik = 1, nks
      !
      IF (done(ik)) CYCLE
      !
      count = 1
      wk_for_ik = wk_new(ik)
      equivalent_with_ik(count) = ik
      done(ik) = .TRUE.
      !
      ! Find list of k points equivalent to ik
      !
      DO isym = 1, nsym
         !
         xkn(:) = invs(:,1,isym) * xk_new(1, ik) + &
                  invs(:,2,isym) * xk_new(2, ik) + &
                  invs(:,3,isym) * xk_new(3, ik)
         IF (magnetic_sym .AND. (t_rev(isym) == 1)) xkn = -xkn
         !
         DO jk = ik+1, nks
            IF (.NOT. done(jk)) THEN
               satm = are_k_equivalent(xk_new(:,jk), xkn)
               IF (minus_q) satm = satm .OR. are_k_equivalent(xk_new(:,jk), -xkn)
               !
               ! If the k points are equivalent, update the count and mark as done
               !
               IF (satm) THEN
                  !
                  IF (ABS(wk_new(jk)) < 2.d-8) THEN
                     ! If the weight is almost zero, skip this k point from averaging to                    ! avoid k points generated by add_additional_kpoints being averaged
                     ! with the original k points.
                     !
                     done(jk) = .TRUE.
                     !
                  ELSE
                     count = count + 1
                     wk_for_ik = wk_for_ik + wk_new(jk)
                     equivalent_with_ik(count) = jk
                     done(jk) = .TRUE.
                  ENDIF
                  !
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      !
      ! Set the average weight for all equivalent k points
      !
      IF (count > 1) THEN
         DO jk = 1, count
            wk_new(equivalent_with_ik(jk)) = wk_for_ik / REAL(count, DP)
         ENDDO
      ENDIF
      !
   ENDDO
   !
   DEALLOCATE(equivalent_with_ik)
   DEALLOCATE(done)
   !
   ! ... Divide the weights by nrot, because each k point is unfolded nrot times
   !
   wk_new = wk_new / nrot
   !
   ! Convert xk_new from crystal to cartesian axis, copy to xk and wk
   !
   wk(:) = 0.d0
   wk(1:nks) = wk_new(1:nks)
   xk(:, 1:nks) = xk_new(:, 1:nks)
   CALL cryst_to_cart(nks, xk, bg, +1)  ! Convert from crystal to Cartesian axis
   !
   ! Check the total weight did not change
   !
   IF ( ABS( SUM(wk(1:nks)) - wk_sum ) > 1.0d-10 ) THEN
      CALL errore('irreducible_BZ','weights changed',1)
   ENDIF
   !
   ! normalize weights to one
   !
   one = SUM( wk(1:nks) )
   IF ( one > 0.d0 ) wk(1:nks) = wk(1:nks) / one
   !
   RETURN
   !
CONTAINS
   !
   FUNCTION are_k_equivalent(k1, k2) RESULT(res)
      REAL(DP), INTENT(IN) :: k1(3), k2(3)
      LOGICAL :: res
      res = ABS(  k1(1) - k2(1) - NINT( k1(1) - k2(1) ) ) < 1.0d-5 .AND. &
            ABS(  k1(2) - k2(2) - NINT( k1(2) - k2(2) ) ) < 1.0d-5 .AND. &
            ABS(  k1(3) - k2(3) - NINT( k1(3) - k2(3) ) ) < 1.0d-5
   END FUNCTION are_k_equivalent
   !
END SUBROUTINE irreducible_BZ
