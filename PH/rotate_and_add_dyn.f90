!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine rotate_and_add_dyn (phi, phi2, nat, isym, s, invs, irt, &
     rtau, sxq)
  !-----------------------------------------------------------------------
  !  Rotates a dynamical matrix (phi) in crystal coordinates according
  !  to the specified symmetry operation and add the rotated matrix
  !  to phi2.   phi is left unmodified.
  !
  USE kinds, only : DP
  USE constants, ONLY : tpi
  implicit none
  ! input variables

  integer :: nat, isym, s (3, 3, 48), invs (48), irt (48, nat)
  ! number of atoms in the unit cell
  ! index of the symm.op.
  ! the symmetry operations
  ! index of the inverse operations
  ! index of the rotated atom

  complex(DP) :: phi (3, 3, nat, nat), phi2 (3, 3, nat, nat)
  ! the input dyn.mat. in crystal coordinates
  ! the rotated dyn.mat. in crystal coordinates

  real(DP) :: rtau (3, 48, nat), sxq (3)
  ! for eaxh atom and rotation gives the R vector
  !involved
  ! the rotated q involved in this sym.op.
  !  local variables
  integer :: na, nb, sna, snb, ism1, i, j, k, l
  ! counters on atoms
  ! indices of rotated atoms
  ! index of the inverse symm.op.
  ! generic counters
  real(DP) :: arg
  ! argument of the phase
  complex(DP) :: phase, work


  ism1 = invs (isym)
  do na = 1, nat
     do nb = 1, nat
        sna = irt (isym, na)
        snb = irt (isym, nb)
        arg = (sxq (1) * (rtau (1, isym, na) - rtau (1, isym, nb) ) &
             + sxq (2) * (rtau (2, isym, na) - rtau (2, isym, nb) ) + sxq (3) &
             * (rtau (3, isym, na) - rtau (3, isym, nb) ) ) * tpi
        phase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
        do i = 1, 3
           do j = 1, 3
              work = CMPLX(0.d0, 0.d0,kind=DP)
              do k = 1, 3
                 do l = 1, 3
                    work = work + s (i, k, ism1) * s (j, l, ism1) * phi (k, l, na, nb) &
                         * phase
                 enddo
              enddo
              phi2 (i, j, sna, snb) = phi2 (i, j, sna, snb) + work
           enddo
        enddo
     enddo
  enddo
  !
  return
end subroutine rotate_and_add_dyn
