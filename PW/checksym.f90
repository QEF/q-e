!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine checksym (ir, nat, ityp, xau, rau, ft, sym, irt)
  !-----------------------------------------------------------------------
  !
  !   This routine receives as input all the atomic positions xau,
  !   and the rotated rau by the symmetry operation ir. It sets to true
  !   sym(ir) if for each atom na, it is possible to find an atom nb
  !   which is of the same type of na, and coincide with it after the
  !   symmetry operation. Fractional translations are allowed.
  !
  !   Revised layout 1 may 1995 by A. Dal Corso
  !
  USE kinds
  implicit none
  !
  !     first the dummy variables
  !
  integer :: nat, ityp (nat), irt (48, nat), ir
  ! input: the total number of atoms
  ! input: the type of each atom
  ! output: the rotated of each atom
  ! input: the rotation to be tested
  real(DP) :: xau (3, nat), rau (3, nat), ft (3)
  ! input: the initial vectors
  ! input: the rotated vectors
  ! input: the possible fractionary translat
  logical :: sym (48)
  ! output: if true this is a symmetry opera
  !
  !  few local variables
  !
  integer :: na, nb
  ! counter on atoms
  ! counter on atoms
  logical :: eqvect
  ! the testing function

  external eqvect
  do na = 1, nat
     do nb = 1, nat
        sym (ir) = ityp (na) .eq.ityp (nb) .and.eqvect (rau (1, na), &
             xau (1, nb), ft)
        if (sym (ir) ) then
           !
           ! the rotated atom does coincide with one of the like atoms
           ! keep track of which atom the rotated atom coincides with
           !
           irt (ir, na) = nb
           goto 10
        endif
     enddo
     !
     ! the rotated atom does not coincide with any of the like atoms
     ! s(ir) + ft is not a symmetry operation
     !
     return
10   continue
  enddo
  !
  ! s(ir) + ft is a symmetry operation
  !
  return
end subroutine checksym
