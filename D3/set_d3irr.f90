!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine set_d3irr
  !-----------------------------------------------------------------------
  !
  ! It computes a basis for all the irreducible representations of the
  ! group of the crystal, which are contained in the representation
  ! which has as basis the displacement vectors.
  ! This basis will be used for those quantities that depend on the
  ! q=0 perturbation.
  !
  ! Receives in input:   nsymg0, s, invs, irt, rtau
  ! Calculates:          ug0, tg0, npertg0, nirrg0, irgq
  !
  ! NB: It assumes that the phonon calculation for the q=0 case, has been
  ! performed with iswitch=-2. If this is not the case the following
  ! routine does not work.
  !
#include "machine.h"
  use pwcom
  use phcom
  use d3com

  implicit none
  integer :: w_nsymq, w_irotmq
  ! work array
  ! work array

  real (8) :: zero (3), w_gi (3, 48), w_gimq (3)
  ! a null vector
  ! work array
  complex (8) :: w_tmq (3, 3, 3 * nat)
  ! work array

  logical :: w_minus_q
  ! work array

  call setv (3, 0.d0, zero, 1)

  w_minus_q = .true.
  if (nsymg0.gt.1) then
     call io_pattern(fild0rho,nirrg0,npertg0,ug0,-1)
     call set_sym_irr (nat, at, bg, zero, s, invs, nsymg0, rtau, irt, &
          irgq, w_nsymq, w_minus_q, w_irotmq, tg0, w_tmq, ug0, npertg0, &
          nirrg0, w_gi, w_gimq, iverbosity)
  else
     call set_irr_nosym (nat, at, bg, zero, s, invs, nsymg0, rtau, &
          irt, irgq, w_nsymq, w_minus_q, w_irotmq, tg0, w_tmq, ug0, &
          npertg0, nirrg0, w_gi, w_gimq, iverbosity)
  endif

  return
end subroutine set_d3irr
