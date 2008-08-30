!
! Copyright (C) 2001-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine sgama_d3 (nrot, nsymq, nat, s, ityp, nr1, nr2, nr3, nsymg0, &
     irt, ftau, at, bg, tau)
  !-----------------------------------------------------------------------
  !
  ! It calculates and/or reorder:  nsymg0, s, irt, ftau
  !
  ! Matrices of the symmetry operations -s- are read from the data file.
  ! They are calculated by the "sgama" routine and are ordered in this way:
  !  a) the first nrot matrices are symmetries of the lattice
  !  b) the first nsymq matrices are symmetries for the small group of q
  !
  ! This routine finds which symmetries of the lattice are also symmetries
  ! of the crystal, calculates the order of the crystal group:   nsymg0
  ! and reorder the s matrices in this way:
  !  a) the first nsymg0 matrices are symmetries of the crystal
  !  b) the first nsymq matrices are symmetries for the small group of q
  !
#include "f_defs.h"
  USE kinds, only : DP
  implicit none

  integer, intent(in) :: nrot, nsymq, nat, ityp (nat), nr1, nr2, nr3
  ! nrot : number of symmetry operations of the bravais lattice
  ! nsymq: order of the small group of q
  ! nat  : number of atoms in the cell
  ! ityp : type of each atom
  ! nr1,2,3: dimension of the FFT mesh
  integer, intent(inout) :: s (3, 3, 48)
  ! s    : matrices of the symmetry operations
  integer, intent(out) :: nsymg0, irt (48, nat), ftau (3, 48)
  ! nsymg0: order of the crystal group
  ! irt   : for each atom, the rotated atom
  ! output: fractionary translation of each sym
  real (DP), intent(in) :: at (3, 3), bg (3, 3), tau (3, nat)
  ! at: direct lattice vectors
  ! bg: reciprocal lattice vectors
  ! tau: coordinates of atomic positions
  !
  integer :: irot, jrot, ipol, jpol, na
  ! counters: on rotations, on polarizations, on atoms
  logical :: sym (48)
  ! if true the symmetry is a true symmetry
  !
  !
  ! find the true symmetries of the crystal
  !
  call sgam_at (nrot, s, nat, tau, ityp, at, bg, nr1, nr2, nr3, sym, &
       irt, ftau)
  !
  ! copy symm. operation in sequential order so that:
  !               irot <= nsymq   are sym.ops. of the small group of q
  !    nsymq+1 <= irot <= nsymg0  are sym.ops. of the crystal
  !
  do irot = 1, nsymq
     if (.not.sym (irot) ) call errore ('sgama_d3', 'unexpected', 1)
  enddo
  jrot = nsymq
  do irot = nsymq + 1, nrot
     if (sym (irot) ) then
        jrot = jrot + 1
        do ipol = 1, 3
           do jpol = 1, 3
              s (ipol, jpol, jrot) = s (ipol, jpol, irot)
           enddo
           ftau (ipol, jrot) = ftau (ipol, irot)
        enddo
        do na = 1, nat
           irt (jrot, na) = irt (irot, na)
        enddo
     endif
  enddo

  nsymg0 = jrot
  return
end subroutine sgama_d3
