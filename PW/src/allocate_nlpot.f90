!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE allocate_nlpot
  !-----------------------------------------------------------------------
  !
  ! This routine computes the dimension of the Hamiltonian matrix and
  ! allocates arrays containing the non-local part of the pseudopotential
  !
  ! It computes the following global quantities:
  !
  !     npwx          !  maximum number of plane waves
  !     nqx           !  number of points of the interpolation table
  !     nqxq          !  as above, for q-function interpolation table
  !
  !
  USE ions_base,        ONLY : nat, nsp, ityp
  USE cellmd,           ONLY : cell_factor
  USE gvect,            ONLY : ngm, gcutm, g
  USE klist,            ONLY : xk, wk, nks, qnorm
  USE lsda_mod,         ONLY : nspin
  USE ldaU,             ONLY : Hubbard_lmax
  USE scf,              ONLY : rho
  USE noncollin_module, ONLY : noncolin
  USE wvfct,            ONLY : npwx, g2kin
  USE gvecw,            ONLY : gcutw, ecutwfc
  USE us,               ONLY : qrad, tab, tab_d2y, tab_at, dq, nqx, &
                               nqxq, spline_ps
  USE uspp,             ONLY : indv, nhtol, nhtolm, ijtoh, qq, dvan, deeq, &
                               vkb, indv_ijkb0, okvan, nkb, nkbus, nhtoj, &
                               becsum, qq_so,dvan_so, deeq_nc
  USE uspp_param,       ONLY : upf, lmaxq, lmaxkb, nh, nhm, nbetam
  USE spin_orb,         ONLY : lspinorb, fcoef
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: n_plane_waves
  INTEGER :: nwfcm
  !
  !   calculate number of PWs for all kpoints
  !
  npwx = n_plane_waves (gcutw, nks, xk, g, ngm)
  !
  !   g2kin contains the kinetic energy \hbar^2(k+G)^2/2m
  !
  ALLOCATE (g2kin ( npwx ) )
  !
  ! Note: computation of the number of beta functions for
  ! each atomic type and the maximum number of beta functions
  ! and the number of beta functions of the solid has been
  ! moved to init_run.f90 : pre_init()
  !
  ALLOCATE (indv( nhm, nsp))
  ALLOCATE (nhtol(nhm, nsp))
  ALLOCATE (nhtolm(nhm, nsp))
  ALLOCATE (nhtoj(nhm, nsp))
  ALLOCATE (ijtoh(nhm, nhm, nsp))
  ALLOCATE (indv_ijkb0(nat))
  ALLOCATE (deeq( nhm, nhm, nat, nspin))
  IF (noncolin) THEN
     ALLOCATE (deeq_nc( nhm, nhm, nat, nspin))
  ENDIF
  ALLOCATE (qq(   nhm, nhm, nsp))
  IF (lspinorb) THEN
    ALLOCATE (qq_so(nhm, nhm, 4, nsp))
    ALLOCATE (dvan_so( nhm, nhm, nspin, nsp))
    ALLOCATE (fcoef(nhm,nhm,2,2,nsp))
  ELSE
    ALLOCATE (dvan( nhm, nhm, nsp))
  ENDIF
  ! GIPAW needs a slighly larger q-space interpolation for quantities calculated
  ! at k+q_gipaw, and I'm using the spline_ps=.true. flag to signal that
  IF (spline_ps .and. cell_factor <= 1.1d0) cell_factor = 1.1d0
  !
  ! This routine is called also by the phonon code, in which case it should
  ! allocate an array that includes q+G vectors up to |q+G|_max <= |Gmax|+|q|
  !
  nqxq = int( ( (sqrt(gcutm) + qnorm) / dq + 4) * cell_factor )
  lmaxq = 2*lmaxkb+1
  !
  IF (lmaxq > 0) ALLOCATE (qrad( nqxq, nbetam*(nbetam+1)/2, lmaxq, nsp))
  ALLOCATE (vkb( npwx,  nkb))
  ALLOCATE (becsum( nhm * (nhm + 1)/2, nat, nspin))
  !
  ! Calculate dimensions for array tab (including a possible factor
  ! coming from cell contraction during variable cell relaxation/MD)
  !
  nqx = int( (sqrt (ecutwfc) / dq + 4) * cell_factor )

  ALLOCATE (tab( nqx , nbetam , nsp))

  ! d2y is for the cubic splines
  IF (spline_ps) ALLOCATE (tab_d2y( nqx , nbetam , nsp))

  nwfcm = maxval ( upf(1:nsp)%nwfc )
  ALLOCATE (tab_at( nqx , nwfcm , nsp))

  RETURN
END SUBROUTINE allocate_nlpot

