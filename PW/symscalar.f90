!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine symscalar (nat, scalar, nsym, s, irt)
   !-----------------------------------------------------------------------
   !
   ! This routine symmetrizes an atom dependent scalar quantity
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
   real(kind=DP) :: scalar(nat)  ! inp/out: the scalar to rotate
   !
   !   the local variables
   !
   integer :: na,              & ! counter on atoms
              nar,             & ! the rotated of each atom
              isym,            & ! counter on symmetries
              ipol               ! counter on polarization

   real(kind=DP), allocatable :: work (:)

   external DCOPY, DSCAL

   if (nsym.eq.1) return

   allocate(work(nat))
   work(:) = 0.d0
   do na = 1, nat
      do isym = 1, nsym
         nar = irt (isym, na)
         work (na) = work (na) +  scalar(nar)
      enddo
   enddo
   call DSCAL (nat, 1.d0 / nsym, work, 1)
   call DCOPY (nat, work, 1, scalar, 1)

   deallocate (work)
   return
end subroutine symscalar

