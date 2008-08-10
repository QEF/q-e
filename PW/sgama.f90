!
! Copyright (C) 2008 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE sgama2 ( nrot, nat, s, sname, t_rev, at, bg, tau, ityp,  &
                    nsym, nr1, nr2, nr3, irt, ftau, invsym,         &
                    magnetic_sym, m_loc )
  !-----------------------------------------------------------------------
  !
  !     This routine finds the point group of the crystal, by eliminating
  !     the symmetries of the Bravais lattice which are not allowed
  !     by the atomic positions (or by the magnetization if present)
  !
  USE kinds, only : DP
  implicit none
  !
  integer, intent(in) :: nrot, nat, ityp (nat), nr1, nr2, nr3  
  real(DP), intent(in) :: at (3,3), bg (3,3), tau (3,nat), m_loc(3,nat)
  logical, intent(in) :: magnetic_sym
  !
  character(len=45), intent(inout) :: sname (48)
  ! name of the rotation part of each symmetry operation
  integer, intent(inout) :: s(3,3,48)
  !
  integer, intent(out) :: nsym, irt (48, nat), ftau (3, 48)
  logical, intent(out) :: invsym
  !
  integer :: t_rev(48)
  ! for magnetic symmetries: if 1 there is time reversal operation
  logical :: sym (48)
  ! if true the corresponding operation is a symmetry operation
  integer, external :: copy_sym
  !
  !    Here we find the true symmetries of the crystal
  !
  CALL sgam_at (nrot, s, nat, tau, ityp, at, bg, nr1, nr2, nr3, &
                sym, irt, ftau)
  !
  !    Here we check for magnetic symmetries
  !
  IF ( magnetic_sym ) &
     CALL sgam_at_mag (nrot, s, nat, bg, irt, m_loc, sname, sym, t_rev)
  !
  !    Here we re-order all rotations in such a way that true sym.ops 
  !    are the first nsym; rotations that are not sym.ops. follow
  !
  nsym = copy_sym ( nrot, sym, s, sname, ftau, nat, irt, t_rev ) 
  !
  ! check if inversion (I) is a symmetry.
  ! If so, it should be the (nsym/2+1)-th operation of the group
  !
  invsym = ALL ( s(:,:,nsym/2+1) == -s(:,:,1) )
  !
  return
  !
END SUBROUTINE sgama2
!
!-----------------------------------------------------------------------
INTEGER FUNCTION copy_sym ( nrot, sym, s, sname, ftau, nat, irt, t_rev ) 
!-----------------------------------------------------------------------
  !
  implicit none
  integer,  intent(in) :: nrot, nat
  integer, intent(inout) :: s (3,3,48), irt (48, nat), ftau (3, 48), &
                            t_rev(48)
  character(len=45), intent(inout) :: sname (48)
  logical, intent(inout) :: sym(48)
  !
  integer :: stemp(3,3), irot, jrot
  !
  ! copy symm. operations in sequential order so that
  ! s(i,j,irot) , irot <= nsym          are the sym.ops. of the crystal
  !               nsym+1 < irot <= nrot are the sym.ops. of the lattice
  jrot = 0
  do irot = 1, nrot
     if (sym (irot) ) then
        jrot = jrot + 1
        stemp = s(:,:,jrot)
        s (:,:, jrot) = s (:,:, irot)
        s (:,:, irot) = stemp
        ftau (:, jrot) = ftau (:, irot)
        irt (jrot,1:nat) = irt (irot,1:nat)
        sname (jrot) = sname (irot)
        t_rev (jrot) = t_rev(irot)
     endif
  enddo
  sym (1:jrot) = .true.
  sym (jrot+1:nrot) = .false.
  !
  ! Sets to zero the first matrix that is not a symmetry of the crystal.
  ! This will be used by d3toten program (obsolete?)
  !
  if (nrot < 48) s(:,:, nrot+1) = 0
  !
  copy_sym = jrot
  return
  !
END FUNCTION copy_sym
!
!-----------------------------------------------------------------------
subroutine irreducible_BZ (nrot, s, nsym, at, bg, npk, nks, xk, wk, minus_q)
  !-----------------------------------------------------------------------
  !
  !     This routine finds the special points in the irreducible wedge of the
  !     true point group (or small group of q) of the crystal, starting
  !     from the points in the irreducible wedge of the point group
  !     of the Bravais lattice.
  !
  USE kinds, only : DP
  implicit none
  !
  integer,  intent(in) :: nrot, nsym, npk, s(3,3,48)
  real(DP), intent(in) :: at (3,3), bg (3,3)
  logical,  intent(in) :: minus_q
  integer,  intent(inout) :: nks
  real(DP), intent(inout) :: xk (3, npk), wk (npk)
  !
  integer :: table (48, 48), invs (3, 3, 48), irg (48)
  ! table: multiplication table of the group
  ! invs : contains the inverse of each rotation
  ! irg  : gives the correspondence of symmetry operations forming a n-th coset
  logical :: sym(48)
  !
  !    We compute the multiplication table of the group
  !
  call multable (nrot, s, table)
  !
  !   And we set the matrices of the inverse
  !
  call inverse_s (nrot, s, table, invs)
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
  call irrek (npk, nks, xk, wk, at, bg, nrot, invs, nsym, irg, minus_q)
  !
  return
  !
end subroutine irreducible_BZ 
