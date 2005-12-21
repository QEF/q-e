!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine dynmat0
  !-----------------------------------------------------------------------
  !
  !     This routine computes the part of the dynamical matrix which
  !     does not depend upon the change of the Bloch wavefunctions.
  !     It is a driver which calls the routines dynmat_## and d2ionq
  !     for computing respectively the electronic part and
  !     the ionic part
  !
  !
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat,ntyp => nsp, ityp, zv, tau
  USE cell_base, ONLY: alat, omega, at, bg
  USE gvect, ONLY: g, gg, ngm, gcutm
  USE symme, ONLY: irt, s
  USE control_flags, ONLY : modenum
  USE kinds,         ONLY : DP
  use phcom
  implicit none

  integer :: nu_i, nu_j, na_icart, nb_jcart
  ! counters

  complex(DP) :: wrk, dynwrk (3 * nat, 3 * nat)
  ! auxiliary space

  call start_clock ('dynmat0')
  call ZCOPY (9 * nat * nat, dyn00, 1, dyn, 1)
  !
  ! first electronic contribution arising from the term  <psi|d2v|psi>
  !
  call dynmat_us
  !
  !   Here the ionic contribution
  !
  call d2ionq (nat, ntyp, ityp, zv, tau, alat, omega, xq, at, bg, g, &
       gg, ngm, gcutm, nmodes, u, dyn)
  !
  !   Add non-linear core-correction (NLCC) contribution (if any)
  !
  call dynmatcc
  !
  !   Symmetrizes the dynamical matrix w.r.t. the small group of q and of
  !   mode. This is done here, because this part of the dynmical matrix is
  !   saved with recover and in the other runs the symmetry group might change
  !
  if (modenum .ne. 0) then

     call symdyn_munu (dyn, u, xq, s, invs, rtau, irt, irgq, at, bg, &
          nsymq, nat, irotmq, minus_q)
     !
     ! rotate again in the pattern basis
     !
     call ZCOPY (9 * nat * nat, dyn, 1, dynwrk, 1)
     do nu_i = 1, 3 * nat
        do nu_j = 1, 3 * nat
           wrk = (0.d0, 0.d0)
           do nb_jcart = 1, 3 * nat
              do na_icart = 1, 3 * nat
                 wrk = wrk + CONJG(u (na_icart, nu_i) ) * &
                             dynwrk (na_icart, nb_jcart) * &
                             u (nb_jcart, nu_j)
              enddo
           enddo
           dyn (nu_i, nu_j) = wrk
        enddo
     enddo
     call ZCOPY (9 * nat * nat, dyn, 1, dyn00, 1)
  endif
  !      call tra_write_matrix('dynmat0 dyn',dyn,u,nat)
  call stop_clock ('dynmat0')
  return
end subroutine dynmat0
