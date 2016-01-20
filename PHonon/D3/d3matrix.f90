!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine d3matrix
  !-----------------------------------------------------------------------
  !
  ! This routine is driver which computes the symmetrized derivative
  ! of the dynamical matrix at q and in the star of q.
  ! The result is written on a iudyn file
  !
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp, tau, atm, amass
  USE run_info, ONLY : title
  USE kinds, only : DP
  use pwcom
  USE symm_base, ONLY : s, irt, invs
  USE control_flags, ONLY : modenum
  USE qpoint, ONLY : xq
  use phcom
  use d3com

  USE lr_symm_base, ONLY : nsymq, minus_q, irotmq, irgq, rtau

  implicit none

  integer :: nq, isq (48), imq, na, nt, j
  ! degeneracy of the star of q
  ! index of q in the star of a given sym.op.
  ! index of -q in the star of q (0 if not present)
  ! counter on atoms
  ! counter on atomic type
  ! generic counter

  real (DP) :: sxq (3, 48)
  ! list of vectors in the star of q
  !
  ! Symmetrizes the dynamical matrix w.r.t. the small group of q
  !
  call d3_symdyn (d3dyn, u, ug0, xq, s, invs, rtau, irt, irgq, at, &
       bg, nsymq, nat, irotmq, minus_q, npert_i, npert_f)
  !
  ! Generates the star of q
  !
  call star_q (xq, at, bg, nsymg0, s, invs, nq, sxq, isq, imq, .TRUE.)
  !
  ! Write on file information on the system
  !
  write (iudyn, '("Derivative of the force constants")')
  write (iudyn, '(a)') title
  write (iudyn, '(i3,i5,i3,6f11.7)') ntyp, nat, ibrav, celldm
  do nt = 1, ntyp
     write (iudyn, * ) nt, " '", atm (nt) , "' ", amass (nt)
  enddo
  do na = 1, nat
     write (iudyn, '(2i5,3f15.7)') na, ityp (na) , (tau (j, na) , j = &
          1, 3)
  enddo
  !
  ! Rotates and writes on iudyn the dyn.matrix derivative of the star of q
  !

  call qstar_d3 (d3dyn, at, bg, nat, nsymg0, s, invs, irt, rtau, nq, &
       sxq, isq, imq, iudyn, wrmode)
  return
end subroutine d3matrix
