!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine allocate_phq
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities needed for the linear
  ! response problem
  !
#include "f_defs.h"

  USE ions_base, ONLY : nat, ntyp => nsp
  use pwcom
  USE noncollin_module, ONLY : noncolin, npol
  USE wavefunctions_module,  ONLY: evc
  USE kinds, only : DP
  USE spin_orb, ONLY : lspinorb
  use phcom
  use el_phon
  USE becmod, ONLY: becp, becp_nc
  USE uspp_param, ONLY: nhm
  USE ramanm, ONLY: ramtns, lraman
  implicit none
  !
  !  allocate space for the quantities needed in the phonon program
  !
  if (lgamma) then
     !
     !  q=0  : evq and igkq are pointers to evc and igk
     !
     evq  => evc
     igkq => igk
  else
     !
     !  q!=0 : evq, igkq are allocated and calculated at point k+q
     !
     allocate (evq ( npwx*npol , nbnd))    
     allocate (igkq ( npwx))    
  endif
  !
  allocate (dvpsi ( npwx*npol , nbnd))    
  allocate ( dpsi ( npwx*npol , nbnd))    
  !
  allocate (vlocq ( ngm , ntyp))    
  allocate (dmuxc ( nrxx , nspin , nspin))    
  !
  allocate (eigqts ( nat))    
  allocate (rtau ( 3, 48, nat))    
  allocate (u ( 3 * nat, 3 * nat))    
  allocate (ubar ( 3 * nat))    
  allocate (dyn ( 3 * nat, 3 * nat))    
  allocate (dyn00 ( 3 * nat, 3 * nat))    
  allocate (w2 ( 3 * nat))    
  allocate (t (max_irr_dim, max_irr_dim, 48,3 * nat))    
  allocate (tmq (max_irr_dim, max_irr_dim, 3 * nat))    
  allocate (npert ( 3 * nat))    
  allocate (zstareu (3, 3,  nat))    
  allocate (zstareu0 (3, 3 * nat))    
  allocate (zstarue (3 , nat, 3))    
  allocate (zstarue0 (3 * nat, 3))    
  if (okvan) then
     allocate (int1 ( nhm, nhm, 3, nat, nspin))    
     allocate (int2 ( nhm , nhm , 3 , nat , nat))    
     allocate (int3 ( nhm , nhm , max_irr_dim , nat , nspin))    
     allocate (int4 ( nhm * (nhm + 1)/2,  3 , 3 , nat, nspin))    
     allocate (int5 ( nhm * (nhm + 1)/2 , 3 , 3 , nat , nat))    
     allocate (dpqq( nhm, nhm, 3, ntyp))    
     IF (noncolin) THEN
        ALLOCATE(int1_nc( nhm, nhm, 3, nat, nspin))    
        ALLOCATE(int3_nc( nhm, nhm, max_irr_dim , nat , nspin))    
        ALLOCATE(int4_nc( nhm, nhm, 3, 3, nat, nspin))    
        ALLOCATE(becsum_nc( nhm*(nhm+1)/2, nat, npol, npol))    
        ALLOCATE(alphasum_nc( nhm*(nhm+1)/2, 3, nat, npol, npol))    
        IF (lspinorb) THEN
           ALLOCATE(int2_so( nhm, nhm, 3, nat , nat, nspin))    
           ALLOCATE(int5_so( nhm, nhm, 3, 3, nat , nat, nspin))    
!           allocate(dpqq_so( nhm, nhm, nspin, 3, ntyp))    
        END IF
     END IF
     allocate (alphasum ( nhm * (nhm + 1)/2 , 3 , nat , nspin))    
     allocate (this_dvkb3_is_on_file(nksq))    
     this_dvkb3_is_on_file(:)=.false.
  endif
  allocate (this_pcxpsi_is_on_file(nksq,3))
  this_pcxpsi_is_on_file(:,:)=.false.
  IF (noncolin) THEN
     ALLOCATE(alphap_nc(nkb, npol, nbnd , 3 , nksq))    
     ALLOCATE(becp1_nc(nkb, npol, nbnd, nksq)) 
     ALLOCATE(becp_nc(nkb, npol, nbnd) )
  ELSE
     ALLOCATE( alphap ( nkb , nbnd , 3 , nksq) )    
     ALLOCATE( becp1 (nkb, nbnd, nksq) ) 
     ALLOCATE( becp(nkb, nbnd) )
  END IF
  if (elph) allocate (el_ph_mat( nbnd, nbnd, nks, 3*nat))    
  if (lraman) allocate ( ramtns (3, 3, 3, nat) )
  return
end subroutine allocate_phq
