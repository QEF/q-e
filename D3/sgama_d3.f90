!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine sgama_d3 (nsymq, nat, s, ityp, nr1, nr2, nr3, nsymg0, &
 irt, ftau, at, bg, tau)
!-----------------------------------------------------------------------
!
! It calculates e/o reorder:  nsymg0, s, irt, ftau
!
! Matrices of the symmetry operations -s- are read from the iunpun file.
! They are calculated by the pw/lib/sgama.F routine and are
! ordered in the following way this way:
!  a) the first nrot matrices are symmetries of the lattice
!  b) the first nsymq matrices are symmetries for the small group of q
!
! This routine finds which symmetries of the lattice are also symmetries
! of the crystal,
! it calculates the order of the crystal group:   nsymg0
! and reorder the s matrices in this way:
!  a) the first nsymg0 matrices are symmetries of the crystal
!  b) the first nsymq matrices are symmetries for the small group of q
!
#include "f_defs.h"
USE kinds, only : DP
implicit none

integer :: nsymq, nat, s (3, 3, 48), ityp (nat), nr1, nr2, nr3, &
 nsymg0, irt (48, nat), ftau (3, 48)
                           ! input: order of the small group of q
                           ! input: number of atoms in the cell
                           ! in/out: matrices of the symmetry operations
                           ! input: type of each atom
                           ! input:
                           ! input: dimension of the FFT mesh
                           ! input:
                           ! output: order of the crystal group
                           ! output: for each atom gives the rotated ato
                           ! output: fractionary translation of each sym

real (DP) :: at (3, 3), bg (3, 3), tau (3, nat)
                           ! input: direct lattice vectors
                           ! input: reciprocal lattice vectors
                           ! input: coordinates of atomic positions
!
! local variables
!

integer :: nrot, irot, jrot, ipol, jpol, na
                           ! order of the lattice point group
                           ! counter on the rotations
                           ! counter on the rotations
                           ! counter on the polarizations
                           ! counter on the polarizations
                           ! counter on atoms

logical :: sym (48)
                           ! if true the symmetry is a true symmetry
!
! It calculates the order of the lattice group by finding the first
! singular matrice
!
nrot = 48
do irot = 1, 48
if ( (s (1, 1, irot) .eq.0) .and. (s (2, 1, irot) .eq.0) .and. (s &
 (3, 1, irot) .eq.0) ) then
   nrot = irot - 1
   goto 10
endif
enddo
   10 continue
!
! It finds the true symmetries of the crystal
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
