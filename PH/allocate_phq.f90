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
#include "machine.h"

  USE ions_base, ONLY : nat, ntyp => nsp
  use pwcom
  USE wavefunctions_module,  ONLY: evc
  USE kinds, only : DP
  use phcom
  use el_phon
  USE becmod, ONLY: becp
  USE uspp_param, ONLY: nhm
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
     allocate (evq ( npwx , nbnd))    
     allocate (igkq ( npwx))    
  endif
  !
  allocate (dvpsi ( npwx , nbnd))    
  allocate ( dpsi ( npwx , nbnd))    
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
     allocate (int3 ( nhm , nhm , 3 , nat , nspin))    
     allocate (int4 ( nhm * (nhm + 1)/2,  3 , 3 , nat, nspin))    
     allocate (int5 ( nhm * (nhm + 1)/2 , 3 , 3 , nat , nat))    
     allocate (dpqq( nhm, nhm, 3, ntyp))    
     allocate (alphasum ( nhm * (nhm + 1)/2 , 3 , nat , nspin))    
     allocate (this_dvkb3_is_on_file(nksq))    
     this_dvkb3_is_on_file(:)=.false.
  endif
  allocate (this_pcxpsi_is_on_file(nksq,3))
  this_pcxpsi_is_on_file(:,:)=.false.
  allocate ( alphap ( nkb , nbnd , 3 , nksq))    
  allocate ( becp1 (nkb, nbnd, nksq), becp(nkb, nbnd) )
  if (elph) allocate (el_ph_mat( nbnd, nbnd, nks, 3*nat))    
  return
end subroutine allocate_phq
