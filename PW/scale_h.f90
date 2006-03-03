!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine scale_h
  !-----------------------------------------------------------------------
  ! When variable cell calculation are performed this routine scales the
  ! quantities needed in the calculation of the hamiltonian using the
  ! new and old cell parameters.
  !
  USE io_global,  ONLY : stdout
  USE ions_base,  ONLY : ntyp => nsp
  USE cell_base,  ONLY : bg, omega
  USE cellmd,     ONLY : at_old, omega_old
  USE gvect,      ONLY : g, gg, ngm
  USE klist,      ONLY : xk, wk, nkstot
  USE us,         ONLY : nqxq, nqx, qrad, tab
  !
  implicit none
  !
  integer :: ig
  ! counter on G vectors

  integer :: ik, ipol
  !
  ! scale the k points
  !
  call cryst_to_cart (nkstot, xk, at_old, - 1)
  call cryst_to_cart (nkstot, xk, bg, + 1)
  WRITE( stdout, * ) ' NEW K-POINTS'
  do ik = 1, nkstot
     WRITE( stdout, '(3f12.7,f12.7)') (xk (ipol, ik) , ipol = 1, 3) , wk (ik)
  enddo
  !
  ! scale the g vectors (as well as gg and gl arrays)
  !
  call cryst_to_cart (ngm, g, at_old, - 1)
  call cryst_to_cart (ngm, g, bg, + 1)
  do ig = 1, ngm
     gg (ig) = g(1, ig) * g(1, ig) + g(2, ig) * g(2, ig) + g(3, ig) * g(3, ig)
  enddo
  !
  ! scale the non-local pseudopotential tables
  !
  tab(:,:,:) = tab(:,:,:) * sqrt (omega_old/omega)
  qrad(:,:,:,:) = qrad(:,:,:,:) * omega_old/omega
  !
  ! recalculate the local part of the pseudopotential
  !
  call init_vloc ( )
  !
  return
end subroutine scale_h

