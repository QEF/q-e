!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine hinit1
  !-----------------------------------------------------------------------
  !  Atomic configuration dependent hamiltonian initialization
  !
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,     ONLY : at, bg, omega, tpiba2
  USE cellmd,        ONLY : lmovecell 
  USE gvect,         ONLY : nr1, nr2, nr3, nrxx, ngm, g, &
                            eigts1, eigts2, eigts3
  USE gsmooth,       ONLY : doublegrid
  USE ldaU,          ONLY : lda_plus_u
  USE lsda_mod,      ONLY : nspin
  USE scf,           ONLY : vrs, vltot, vr
  USE vlocal,        ONLY : strf
  USE control_flags, ONLY : order
  !
  implicit none
  !
  !  update the potential
  !
  call update_pot
  !
  !  initialize structure factor array if it has not already been calculat
  !  update_pot ( this is done if order > 0 )
  !
  if (order.eq.0) then
     if (lmovecell) call scale_h
     call struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, &
          nr3, strf, eigts1, eigts2, eigts3)
     !
     !  calculate the core charge (if any) for the nonlinear core correction
     !
     call set_rhoc
  endif
  !
  ! calculate the total local potential
  !
  call setlocal
  !
  ! define the total local potential (external+scf)
  !
  call set_vrs (vrs, vltot, vr, nrxx, nspin, doublegrid)
  !
  ! update the D matrix
  !
  call newd
  !
  ! and recalculate the products of the S with the atomic wfcs used in LDA+U
  ! calculations
  !
  if (lda_plus_u) call orthoatwfc

  return
end subroutine hinit1

