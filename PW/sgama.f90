!
! Copyright (C) 2008 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE sgama ( nrot, nat, s, sname, t_rev, at, bg, tau, ityp,  &
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
END SUBROUTINE sgama
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
  integer :: stemp(3,3), ftemp(3), irtemp(nat), ttemp, irot, jrot
  character(len=45) :: nametemp
  !
  ! copy symm. operations in sequential order so that
  ! s(i,j,irot) , irot <= nsym          are the sym.ops. of the crystal
  !               nsym+1 < irot <= nrot are the sym.ops. of the lattice
  jrot = 0
  do irot = 1, nrot
     if (sym (irot) ) then
        jrot = jrot + 1
        if ( irot > jrot ) then 
           stemp = s(:,:,jrot)
           s (:,:, jrot) = s (:,:, irot)
           s (:,:, irot) = stemp
           ftemp(:) = ftau(:,jrot)
           ftau (:, jrot) = ftau (:, irot)
           ftau (:, irot) = ftemp(:)
           irtemp (:) = irt (jrot,:)
           irt (jrot,:) = irt (irot,:)
           irt (irot,:) = irtemp (:)
           nametemp = sname (jrot)
           sname (jrot) = sname (irot)
           sname (irot) = nametemp
           ttemp = t_rev(jrot)
           t_rev(jrot) = t_rev(irot)
           t_rev(irot) = ttemp
        endif
     endif
  enddo
  sym (1:jrot) = .true.
  sym (jrot+1:nrot) = .false.
  !
  copy_sym = jrot
  return
  !
END FUNCTION copy_sym
