!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dynmatrix
  !-----------------------------------------------------------------------
  !
  ! This routine is a driver which computes the symmetrized dynamical
  ! matrix at q (and in the star of q) and diagonalizes it.
  ! It writes the result on a iudyn file and writes the eigenvalues on
  ! output.
  !
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, atm, amass
  USE io_global,     ONLY : stdout
  USE control_flags, ONLY : modenum, noinv
  USE cell_base,     ONLY : at, bg, celldm, ibrav, symm_type
  USE gvect,         ONLY : nr1, nr2, nr3
  USE printout_base, ONLY : title
  USE symme,         ONLY : s, irt, nsym
  use phcom
  USE ramanm,        ONLY: lraman, ramtns
  implicit none
  ! local variables

  integer :: nq, isq (48), imq, na, nt, imode0, jmode0, irr, jrr, &
       ipert, jpert, mu, nu, i, j
  ! nq :  degeneracy of the star of q
  ! isq: index of q in the star of a given sym.op.
  ! imq: index of -q in the star of q (0 if not present)

  real(DP) :: sxq (3, 48)
  ! list of vectors in the star of q
  !
  !     Puts all noncomputed elements to zero
  !
  call start_clock('dynmatrix')
  imode0 = 0
  do irr = 1, nirr
     jmode0 = 0
     do jrr = 1, nirr
        if (done_irr (irr) .eq.0.and.done_irr (jrr) .eq.0) then
           do ipert = 1, npert (irr)
              mu = imode0 + ipert
              do jpert = 1, npert (jrr)
                 nu = jmode0 + jpert
                 dyn (mu, nu) = CMPLX (0.d0, 0.d0)
              enddo
           enddo
        elseif (done_irr (irr) .eq.0.and.done_irr (jrr) .ne.0) then
           do ipert = 1, npert (irr)
              mu = imode0 + ipert
              do jpert = 1, npert (jrr)
                 nu = jmode0 + jpert
                 dyn (mu, nu) = CONJG(dyn (nu, mu) )
              enddo
           enddo
        endif
        jmode0 = jmode0 + npert (jrr)
     enddo
     imode0 = imode0 + npert (irr)
  enddo
  !
  !   Symmetrizes the dynamical matrix w.r.t. the small group of q
  !

  call symdyn_munu (dyn, u, xq, s, invs, rtau, irt, irgq, at, bg, &
       nsymq, nat, irotmq, minus_q)
  !
  !  if only one mode is computed write the dynamical matrix and stop
  !
  if (modenum .ne. 0) then
     WRITE( stdout, '(/,5x,"Dynamical matrix:")')
     do nu = 1, 3 * nat
        WRITE( stdout, '(5x,2i5,2f10.6)') modenum, nu, dyn (modenum, nu)
     enddo
     call stop_ph (.false.)

  endif
  !
  !   Generates the star of q
  !
  call star_q (xq, at, bg, ibrav, symm_type, nat, tau, ityp, nr1, &
       nr2, nr3, nsym, s, invs, irt, rtau, nq, sxq, isq, imq, noinv, &
       modenum)
  !
  ! write on file information on the system
  !
  write (iudyn, '(a)') title
  write (iudyn, '(a)') title_ph
  write (iudyn, '(i3,i5,i3,6f11.7)') ntyp, nat, ibrav, celldm
  if (ibrav==0) then
     write (iudyn,'(a)') symm_type
     write (iudyn,'(2x,3f15.9)') ((at(i,j),i=1,3),j=1,3)
  end if
  do nt = 1, ntyp
     write (iudyn, * ) nt, ' ''', atm (nt) , ' '' ', amass (nt)
  enddo
  do na = 1, nat
     write (iudyn, '(2i5,3f15.7)') na, ityp (na) , (tau (j, na) , j = 1, 3)
  enddo
  !
  !   Rotates and writes on iudyn the dynamical matrices of the star of q
  !
  call q2qstar_ph (dyn, at, bg, nat, nsym, s, invs, irt, rtau, nq, &
       sxq, isq, imq, iudyn)
  !
  !   Writes (if the case) results for quantities involving electric field
  !
  if (epsil) call write_epsilon_and_zeu (zstareu, epsilon, nat, iudyn)
  if (zue) call sym_and_write_zue
  if (lraman) call write_ramtns (iudyn, ramtns)
  !
  !   Diagonalizes the dynamical matrix at q
  !

  if (all_comp) call dyndia (xq, nmodes, nat, ntyp, ityp, amass, iudyn, dyn, w2)
  call stop_clock('dynmatrix')
  return
end subroutine dynmatrix
