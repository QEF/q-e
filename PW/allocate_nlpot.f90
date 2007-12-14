!
! Copyright (C) 2001-2007 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine allocate_nlpot
  !-----------------------------------------------------------------------
  !
  ! This routine computes the dimension of the Hamiltonian matrix and
  ! allocates arrays containing the non-local part of the pseudopotential
  !
  ! It computes the following global quantities:
  !
  !     ngk           !  number of plane waves (for each k point)
  !     npwx          !  maximum number of plane waves
  !     nkb           !  number of beta functions for the solid
  !     nqx           !  number of points of the interpolation table
  !     nh            !  number of beta functions for each atom type
  !     nhm           !  maximum number of different beta functions
  !     nbetam        !  maximum number of beta functions
  !
  !
  USE ions_base,        ONLY : nat, nsp, ityp
  USE cell_base,        ONLY : tpiba2
  USE cellmd,           ONLY : cell_factor
  USE gvect,            ONLY : ngm, gcutm, ecutwfc, g
  USE klist,            ONLY : xk, wk, ngk, nks, xqq
  USE lsda_mod,         ONLY : nspin
  USE ldaU,             ONLY : Hubbard_lmax
  USE scf,              ONLY : rho
  USE noncollin_module, ONLY : noncolin
  USE wvfct,            ONLY : npwx, npw, igk, g2kin
  USE us,               ONLY : qrad, tab, tab_d2y, tab_at, dq, nqx, &
                               nqxq, spline_ps
  USE uspp,             ONLY : indv, nhtol, nhtolm, ijtoh, qq, dvan, deeq, vkb, &
                               nkb, nkbus, nhtoj, becsum, qq_so, dvan_so, deeq_nc
  USE uspp_param,       ONLY : upf, lmaxq, lmaxkb, nh, nhm, nbetam
  USE spin_orb,         ONLY : lspinorb, fcoef
  USE paw_variables,    ONLY : okpaw
  !
  implicit none
  !
  !    a few local variables
  !
  integer :: nt, na, nb, nwfcm  
  ! counters on atom type, atoms, beta functions
  !
  !   calculate number of PWs for all kpoints
  !
  allocate (ngk( nks ))
  !
  call n_plane_waves (ecutwfc, tpiba2, nks, xk, g, ngm, npwx, ngk)
  !
  !   igk relates the index of PW k+G to index in the list of G vector
  !
  allocate (igk( npwx ), g2kin ( npwx ) )    
  !
  !     calculate the number of beta functions for each atomic type
  !
  lmaxkb = - 1
  do nt = 1, nsp
     nh (nt) = 0
     do nb = 1, upf(nt)%nbeta
        nh (nt) = nh (nt) + 2 * upf(nt)%lll(nb) + 1
        lmaxkb = max (lmaxkb, upf(nt)%lll(nb) )
     enddo
  enddo
  !
  ! calculate the maximum number of beta functions
  !
  nhm = MAXVAL (nh (1:nsp))
  nbetam = MAXVAL (upf(:)%nbeta)
  !
  ! calculate the number of beta functions of the solid
  !
  nkb = 0
  nkbus = 0
  do na = 1, nat
     nt = ityp(na)
     nkb = nkb + nh (nt)
     if (upf(nt)%tvanp) nkbus = nkbus + nh (nt)
  enddo
  !
  allocate (indv( nhm, nsp))    
  allocate (nhtol(nhm, nsp))    
  allocate (nhtolm(nhm, nsp))    
  allocate (nhtoj(nhm, nsp))    
  allocate (ijtoh(nhm, nhm, nsp))
  allocate (deeq( nhm, nhm, nat, nspin))    
  if (noncolin) then
     allocate (deeq_nc( nhm, nhm, nat, nspin))    
  endif
  allocate (qq(   nhm, nhm, nsp))    
  if (lspinorb) then
    allocate (qq_so(nhm, nhm, 4, nsp))    
    allocate (dvan_so( nhm, nhm, nspin, nsp))    
    allocate (fcoef(nhm,nhm,2,2,nsp))
  else
    allocate (dvan( nhm, nhm, nsp))    
  endif
  !
  nqxq = INT( ( (sqrt(gcutm) + sqrt(xqq(1)**2 + xqq(2)**2 + xqq(3)**2) ) &
          / dq + 4) * cell_factor )
  lmaxq = 2*lmaxkb+1
  !
  if (lmaxq > 0) allocate (qrad( nqxq, nbetam*(nbetam+1)/2, lmaxq, nsp))    
  if (nkb > 0) allocate (vkb( npwx,  nkb))    
  allocate (becsum( nhm * (nhm + 1)/2, nat, nspin))    
  ! In PAW becsum has to be treated self-consistently:
  if (okpaw) rho%bec => becsum
  !
  ! Calculate dimensions for array tab (including a possible factor
  ! coming from cell contraction during variable cell relaxation/MD)
  !
  nqx = INT( (sqrt (ecutwfc) / dq + 4) * cell_factor )

  allocate (tab( nqx , nbetam , nsp))

  ! d2y is for the cubic splines
  if (spline_ps) allocate (tab_d2y( nqx , nbetam , nsp))

  nwfcm = MAXVAL ( upf(1:nsp)%nwfc )
  allocate (tab_at( nqx , nwfcm , nsp))

  return
end subroutine allocate_nlpot

