!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)
!-----------------------------------------------------------------------
!
!     This routine computes the vector rtau which contains for each
!     atom and each rotation the vector S\tau_a - \tau_b, where
!     b is the rotated a atom, given by the array irt. These rtau are
!     non zero only if fractional translations are present.
!
!     revised layout 2 may 1995 by Andrea Dal Corso
!
#include "machine.h"
USE kinds
implicit none
!
!     first the dummy variables
!
integer :: nsym, s (3, 3, 48), nat, irt (48, nat)
                                   ! input: symmetries of the point grou
                                   ! input: symmetry matrices
                                   ! input: number of atoms in the unit
                                   ! input: for each atom gives the rota

real(kind=DP) :: at (3, 3), bg (3, 3), tau (3, nat), rtau (3, 48, nat)
                                   ! input: direct lattice vectors
                                   ! input: reciprocal lattice vectors
                                   ! input: coordinates of the atoms
                                   ! output: the direct translations
logical :: sym (nsym)
                                   ! input: if true the symmetry exists
!
!    here the local variables
!
integer :: na, isym, nb, ipol
                                   ! counter on atoms
                                   ! counter on symmetry operations
                                   ! buffer for atom
                                   ! counter on polarization

real(kind=DP) , allocatable :: xau (:,:), rau (:,:)
real(kind=DP) :: ft (3)
                                   ! atomic coordinates in crystal axis
                                   ! rotated atomic coordinates
                                   ! fractionary translation

allocate (xau(3,nat))    
allocate (rau(3,nat))    
!
!   compute the atomic coordinates in crystal axis
!
do na = 1, nat
do ipol = 1, 3
                                                 ! crystal coordinates
xau (ipol, na) = bg (1, ipol) * tau (1, na) + bg (2, ipol) &
 * tau (2, na) + bg (3, ipol) * tau (3, na)
                                                 ! of the current atom
enddo

enddo
!
!    for each symmetry operation
!
do isym = 1, nsym
if (sym (isym) ) then
   do na = 1, nat
   nb = irt (isym, na)
   do ipol = 1, 3
   ft (ipol) = s (1, ipol, isym) * xau (1, na) + s (2, ipol, isym) &
    * xau (2, na) + s (3, ipol, isym) * xau (3, na) - xau (ipol, nb)
   enddo
   do ipol = 1, 3
   rtau (ipol, isym, na) = at (ipol, 1) * ft (1) + at (ipol, 2) &
    * ft (2) + at (ipol, 3) * ft (3)
   enddo
   enddo
endif
enddo
!
!    deallocate workspace
!
deallocate(rau)
deallocate(xau)
return
end subroutine sgam_ph

