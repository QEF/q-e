!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine qstar_d3 (d3dyn, at, bg, nat, nsym, s, invs, irt, rtau, &
     nq, sxq, isq, imq, iudyn, wrmode)
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  USE kinds, only : DP
  implicit none
  !
  ! input variables
  !

  integer :: nat, nsym, s (3, 3, 48), invs (48), irt (48, nat), &
       nq, isq (48), imq, iudyn
  ! number of atoms in the unit cell
  ! number of symmetry operations
  ! the symmetry operations
  ! index of the inverse operations
  ! index of the rotated atom
  ! degeneracy of the star of q
  ! symmetry op. giving the rotated q
  ! index of -q in the star (0 if nont present)
  ! unit number

  complex (kind = dp) :: d3dyn (3 * nat, 3 * nat, 3 * nat)
  ! the dynmatrix derivative

  real (kind = dp) :: at (3, 3), bg (3, 3), rtau (3, 48, nat), sxq (3, 48)
  ! direct lattice vectors
  ! reciprocal lattice vectors
  ! for eaxh atom and rotation gives the
  ! R vector involved
  ! list of q in the star
  logical :: wrmode (3 * nat )
  ! if .true. this mode is to be written
  !
  !  local variables
  !

  integer :: iq, nsq, isym, na, nb, nc, icar, jcar, kcar, i, j, k
  ! counter on q vectors
  ! number of sym.op. giving each q in the list
  ! index of a symm.op.
  ! counters on atoms
  ! counters on atoms
  ! counters on atoms
  ! cartesian coordinate counters
  ! cartesian coordinate counters
  ! cartesian coordinate counters
  ! generic counters

  complex (kind = dp), allocatable :: phi (:,:,:,:,:,:), phi2 (:,:,:,:,:,:)
  ! an auxiliary dyn. matrix.
  ! another one

  allocate  (phi (3,3,3,nat,nat,nat))    
  allocate  (phi2(3,3,3,nat,nat,nat))    
  !
  ! Sets number of symmetry operations giving each q in the list
  !
  nsq = nsym / nq
  if (nsq * nq.ne.nsym) call errore ('qstar_d3', 'wrong degeneracy', &
       1)
  !
  ! Writes dyn.mat d3dyn(3*nat,3*nat,3*nat)
  ! on the 6-index array phi(3,3,3,nat,nat,nat)
  !
  do i = 1, 3 * nat
     na = (i - 1) / 3 + 1
     icar = i - 3 * (na - 1)
     do j = 1, 3 * nat
        nb = (j - 1) / 3 + 1
        jcar = j - 3 * (nb - 1)
        do k = 1, 3 * nat
           nc = (k - 1) / 3 + 1
           kcar = k - 3 * (nc - 1)
           phi (icar, jcar, kcar, na, nb, nc) = d3dyn (i, j, k)
        enddo
     enddo
  enddo
  !
  ! Goes  to crystal coordinates
  !
  do na = 1, nat
     do nb = 1, nat
        do nc = 1, nat
           call trntnsc_3 (phi (1, 1, 1, na, nb, nc), at, bg, - 1)
        enddo
     enddo
  enddo
  !
  ! For each q of the star rotates phi with the appropriate sym.op. -> phi
  !
  do iq = 1, nq
     call setv (2 * 27 * nat * nat * nat, 0.d0, phi2, 1)
     do isym = 1, nsym
        if (isq (isym) .eq.iq) then
           call rotate_and_add_d3 (phi, phi2, nat, isym, s, invs, irt, &
                rtau, sxq (1, iq) )
        endif
     enddo
     call DSCAL (2 * 27 * nat * nat * nat, 1.d0 / nsq, phi2, 1)
     !
     ! Back to cartesian coordinates
     !
     do na = 1, nat
        do nb = 1, nat
           do nc = 1, nat
              call trntnsc_3 (phi2 (1, 1, 1, na, nb, nc), at, bg, + 1)
           enddo
        enddo
     enddo
     !
     ! Writes the dynamical matrix in cartesian coordinates on file
     !

     call write_d3dyn (sxq (1, iq), phi2, nat, iudyn, wrmode)
     if (imq.eq.0) then
        !
        ! if -q is not in the star recovers its matrix by time reversal
        !
        do na = 1, nat
           do nb = 1, nat
              do nc = 1, nat
                 do i = 1, 3
                    do j = 1, 3
                       do k = 1, 3
                          phi2 (i, j, k, na, nb, nc) = conjg (phi2 (i, j, k, na, nb, nc) &
                               )
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
        !
        ! and writes it (changing temporarily sign to q)
        !
        call DSCAL (3, - 1.d0, sxq (1, iq), 1)
        call write_d3dyn (sxq (1, iq), phi2, nat, iudyn, wrmode)

        call DSCAL (3, - 1.d0, sxq (1, iq), 1)

     endif

  enddo
  deallocate (phi)
  deallocate (phi2)
  return
end subroutine qstar_d3
