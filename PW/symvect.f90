!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine symvect (nat, vect, nsym, s, irt)
   !-----------------------------------------------------------------------
   !
   ! This routine symmetrizes a vector (like forces which in the crystal
   ! axis is represented on the reciprocal lattice basis and which depends
   ! on an atomic position) in the crystal
   ! axis basis
   !
#include "f_defs.h"
   USE kinds
   implicit none
   !
   !    I/O variables first
   !
   integer :: nat,             & ! input: the number of atoms in the cell
              nsym,            & ! input: the number of symmetries
              irt (48, nat),   & ! input: for each atom gives the rotated
              s (3, 3, 48)       ! input: the rotation matrices
   real(kind=DP) :: vect(3,nat)  ! inp/out: the vector to rotate
   !
   !   the local variables
   !
   integer :: na,              & ! counter on atoms
              nar,             & ! the rotated of each atom
              isym,            & ! counter on symmetries
              ipol               ! counter on polarization

   real(kind=DP), allocatable :: work (:,:)

   external DCOPY, DSCAL

   if (nsym.eq.1) return

   allocate(work(3,nat))
   work(:,:) = 0.d0
   do na = 1, nat
      do isym = 1, nsym
         nar = irt (isym, na)
         do ipol = 1, 3
            work (ipol, na) = work (ipol, na) + &
                             s (ipol, 1, isym) * vect (1, nar) + &
                             s (ipol, 2, isym) * vect (2, nar) + &
                             s (ipol, 3, isym) * vect (3, nar)
         enddo
      enddo
   enddo
   call DSCAL (3 * nat, 1.d0 / nsym, work, 1)
   call DCOPY (3 * nat, work, 1, vect, 1)

   deallocate (work)
   return
end subroutine symvect

