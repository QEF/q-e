!
! Copyright (C) 2001-2008 Quantum ESPRESSO  group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine irreducible_BZ (nrot, s, nsym, minus_q, at, bg, npk, nks, &
                           xk, wk, t_rev)
  !-----------------------------------------------------------------------
  !
  !     This routine finds the special points in the irreducible wedge of 
  !     the true point group (or small group of q) of the crystal, 
  !     starting from the points in the irreducible BZ wedge 
  !     of the point group of the Bravais lattice.
  !
  USE kinds, only : DP
  implicit none
  !
  integer,  intent(in) :: nrot, nsym, npk, s(3,3,48), t_rev(48)
  real(DP), intent(in) :: at (3,3), bg (3,3)
  logical,  intent(in) :: minus_q
  integer,  intent(inout) :: nks
  real(DP), intent(inout) :: xk (3, npk), wk (npk)
  !
  integer :: table (48, 48), invs (3, 3, 48), irg (48)
  ! table: multiplication table of the group
  ! invs : contains the inverse of each rotation
  ! irg  : gives the correspondence of symmetry operations forming a n-th coset
  integer :: isym, jsym
  logical :: sym(48)
  !
  !    We compute the multiplication table of the group
  !
  call multable (nrot, s, table)
  !
  !   And we set the matrices of the inverse
  !
  DO isym = 1, nrot
     DO jsym = 1, nrot
        IF (table (isym, jsym)==1) invs (:,:,isym) = s(:,:,jsym)
     ENDDO
  ENDDO
  !
  !    Find the coset in the point group of the Bravais lattice
  !
  sym(1:nsym) = .true.
  sym(nsym+1:)= .false.
  call coset (nrot, table, sym, nsym, irg)
  !
  !    here we set the k-points in the irreducible wedge of the point grou
  !    of the crystal
  !
  call irrek (at, bg, nrot, invs, nsym, irg, minus_q, npk, nks, xk, wk, t_rev)
  !
  return
  !
end subroutine irreducible_BZ 
!
!-----------------------------------------------------------------------
subroutine irrek (at, bg, nrot, invs, nsym, irg, minus_q, npk, &
                  nks, xk, wk, t_rev)
  !-----------------------------------------------------------------------
  !
  !  Given a set of special points in the Irreducible Wedge of some
  !  group, finds the equivalent special points in the IW of one of
  !  its subgroups.
  !
  USE kinds, only : DP
  implicit none
  !
  integer, intent(in) :: npk, nrot, nsym, invs (3, 3, 48), irg (nrot)
  ! maximum number of special points
  ! order of the parent point group
  ! order of the subgroup
  ! inverse of the elements of the symmetry group
  ! partition of the elements of the symmetry group into left cosets,
  ! as given by SUBROUTINE COSET
  integer, intent(inout) :: nks
  ! number of special points
  integer, intent(in) :: t_rev(48)
  real(DP), intent(in) :: at (3, 3), bg (3, 3)
  ! basis vectors of the Bravais and reciprocal lattice
  real(DP), intent(inout) :: xk (3, npk), wk (npk)
  ! special points and weights
  logical, intent(in) :: minus_q
  ! .true. if symmetries q = -q+G are acceptable
  !
  !    here the local variables
  !
  integer :: nks0, jk, kpol, irot, jrot, ncos, jc, ic, isym
  ! nks0: used to save the initial number of k-points
  ! ncos: total number of cosets
  real(DP) :: xkg (3), xks (3, 48), w (48), sw, one
  ! coordinates of the k point in crystal axis
  ! coordinates of the rotated k point
  ! weight of each coset
  ! buffer which contains the weight of k points
  ! total weight of k-points
  logical :: latm, satm
  ! true if a k-point is equivalent to a previous one
  ! true if equivalent point found

  nks0 = nks
  do jk = 1, nks0
     !
     !     The k point is first computed in crystal axis
     !
     do kpol = 1, 3
        ! xkg are the components ofx k in the crystal RL base
        xkg (kpol) = at (1, kpol) * xk (1, jk) + &
                     at (2, kpol) * xk (2, jk) + &
                     at (3, kpol) * xk (3, jk)
     enddo
     !
     !   Then it is rotated with each symmetry of the global group. Note that
     !   the irg vector is used to divide all the rotated vector in cosets
     !
     do irot = 1, nrot
        jrot = irg (irot)
        do kpol = 1, 3
           ! the rotated of xkg with respect to the group operations
           xks (kpol, irot) = invs (kpol, 1, jrot) * xkg (1) + &
                              invs (kpol, 2, jrot) * xkg (2) + &
                              invs (kpol, 3, jrot) * xkg (3)
        enddo
        IF (t_rev(jrot)==1)  xks (:, irot)=-xks(:, irot)
     enddo
     !
     !    For each coset one point is tested with all the preceding
     !
     ncos = nrot / nsym
     do ic = 1, ncos
        irot = (ic - 1) * nsym + 1
        latm = .false.
        !
        !  latm = .true. if the present k-vector is equivalent to some previous
        !
        do jc = 1, ic - 1
           do isym = 1, nsym
              !
              !   satm = .true. if the present symmetry operation makes 
              !   the ir and ik k-vectors equivalent ...
              !
              jrot = (jc - 1) * nsym + isym
              satm = abs (xks (1, irot) - xks (1, jrot) - &
                     nint (xks (1, irot) - xks (1, jrot) ) ) < 1.0d-5 .and. &
                     abs (xks (2, irot) - xks (2, jrot) - &
                     nint (xks (2, irot) - xks (2, jrot) ) ) < 1.0d-5 .and. &
                     abs (xks (3, irot) - xks (3, jrot) - &
                     nint (xks (3, irot) - xks (3, jrot) ) ) < 1.0d-5
              !
              !  .... or equivalent to minus each other when minus_q=.t.
              !
              if (minus_q) satm = satm .or. &
                   abs (xks (1, irot) + xks (1, jrot) - &
                   nint (xks (1, irot) + xks (1, jrot) ) ) < 1.0d-5 .and. &
                   abs (xks (2, irot) + xks (2, jrot) - &
                   nint (xks (2, irot) + xks (2, jrot) ) ) < 1.0d-5 .and. &
                   abs (xks (3, irot) + xks (3, jrot) - &
                   nint (xks (3, irot) + xks (3, jrot) ) ) < 1.0d-5
              latm = latm .or. satm
              if (satm .and. w (jc) /= 0.d0) then
                 w (jc) = w (jc) + 1.d0
                 goto 100
              endif
           enddo

        enddo
100     continue
        if (latm) then
           w (ic) = 0.d0
        else
           w (ic) = 1.d0
        endif
     enddo
     !
     !     here the k-point list is updated
     !
     sw = wk (jk) / SUM (w(1:ncos))
     wk (jk) = sw * w (1)
     do ic = 2, ncos
        irot = (ic - 1) * nsym + 1
        if (w (ic) /= 0.d0) then
           nks = nks + 1
           if (nks > npk) call errore ('irrek', 'too many k-points', nks)
           wk (nks) = sw * w (ic)
           do kpol = 1, 3
              xk (kpol, nks) = bg (kpol, 1) * xks (1, irot) + &
                               bg (kpol, 2) * xks (2, irot) + &
                               bg (kpol, 3) * xks (3, irot)
           enddo
        endif
     enddo

  enddo
  !
  ! normalize weights to one
  !
  one = SUM (wk(1:nks))
  if ( one > 0.d0 ) wk(1:nks) = wk(1:nks) / one
  !
  return
end subroutine irrek
