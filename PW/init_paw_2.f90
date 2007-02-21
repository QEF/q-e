!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
subroutine init_paw_2 (npw_, igk_, q_, vkb_)
  !----------------------------------------------------------------------
  !
  !   Calculates paw_beta functions (paw projectors), with
  !   structure factor, for all atoms, in reciprocal space
  !
  USE kinds ,     ONLY : dp
  USE constants , ONLY : tpi
  USE wvfct ,     ONLY : npwx
  USE cell_base , ONLY : tpiba
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE gvect ,     ONLY : eigts1, eigts2, eigts3, g, ig1, ig2, ig3 
  USE us,         ONLY : dq
!  USE paw,        ONLY : paw_nkb, paw_lmaxkb, paw_nh, paw_nhtol, &
!                         paw_nhtom, paw_indv, paw_nbeta, paw_tab, paw_tab_d2y
  USE paw,        ONLY : paw_recon, paw_nkb, paw_lmaxkb
  USE us,         ONLY : nqxq, dq, spline_ps
  USE splinelib
  !
  implicit none
  !
  integer :: npw_, igk_ (npw_)
  ! input: number of PW's
  ! input: indices of q+G
  real(DP) :: q_(3)
  ! input: q vector
  complex(DP) :: vkb_ (npwx, paw_nkb)
  ! output: beta functions
  !
  !     Local variables
  !
  integer :: i0,i1,i2,i3, ig, l, lm, na, nt, nb, ih, jkb

  real(DP) :: px, ux, vx, wx, arg
  real(DP), allocatable :: gk (:,:), qg (:), vq (:), ylm (:,:), vkb1(:,:)

  complex(DP) :: phase, pref
  complex(DP), allocatable :: sk(:)

  real(DP), allocatable :: xdata(:)
  integer :: iq
  
  !
  
  if ( paw_lmaxkb < 0 ) return
  call start_clock ('init_paw_2')
  allocate (  sk( npw_))    
  allocate (  qg( npw_))    
  allocate (  vq( npw_))    
  allocate ( ylm( npw_, (paw_lmaxkb + 1) **2))    
  allocate (  gk( 3, npw_))    
  
  do ig = 1, npw_
     gk (1,ig) = q_(1) + g(1, igk_(ig) )
     gk (2,ig) = q_(2) + g(2, igk_(ig) )
     gk (3,ig) = q_(3) + g(3, igk_(ig) )
     qg (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  end do
  !
  call ylmr2 ( (paw_lmaxkb+1)**2, npw_, gk, qg, ylm )
  !
  ! set now qg=|q+G| in atomic units
  !
  do ig = 1, npw_
     qg(ig) = sqrt(qg(ig))*tpiba
  end do
  
  if (spline_ps) then
    allocate(xdata(nqxq))
    do iq = 1, nqxq
      xdata(iq) = (iq - 1) * dq
    end do
  end if

  jkb = 0
  do nt = 1, ntyp
     
     allocate ( vkb1(npw_,paw_recon(nt)%paw_nbeta) )
     
     ! calculate beta in G-space using an interpolation table
     do nb = 1, paw_recon(nt)%paw_nbeta
        do ig = 1, npw_
           if (spline_ps) then
              vq(ig) = splint ( xdata, paw_recon(nt)%paw_tab(:,nb), &
                   paw_recon(nt)%paw_tab_d2y(:,nb), qg(ig) )
           else
             px = qg (ig) / dq - int (qg (ig) / dq)
             ux = 1.d0 - px
             vx = 2.d0 - px
             wx = 3.d0 - px
             i0 = qg (ig) / dq + 1
             i1 = i0 + 1
             i2 = i0 + 2
             i3 = i0 + 3
             vq (ig) = paw_recon(nt)%paw_tab(i0,nb) * ux * vx * wx / 6.d0 + &
                       paw_recon(nt)%paw_tab(i1,nb) * px * vx * wx / 2.d0 - &
                       paw_recon(nt)%paw_tab(i2,nb) * px * ux * wx / 2.d0 + &
                       paw_recon(nt)%paw_tab(i3,nb) * px * ux * vx / 6.d0
           end if
        end do
        ! add spherical harmonic part
        do ih = 1, paw_recon(nt)%paw_nh
           if ( nb == paw_recon(nt)%paw_indv(ih) ) then
              l = paw_recon(nt)%paw_nhtol(ih)
              lm = l * l + paw_recon(nt)%paw_nhtom(ih)
              do ig = 1, npw_
                 vkb1(ig,ih) = ylm(ig,lm) * vq(ig)
              end do 
           end if
        end do
     end do
     !
     ! vkb1 contains all betas including angular part for type nt
     ! now add the structure factor and factor (-i)^l
     !
     do na = 1, nat
        ! ordering: first all betas for atoms of type 1
        !           then  all betas for atoms of type 2  and so on
        if ( ityp(na) == nt ) then
           arg = (q_(1) * tau (1, na) + &
                q_(2) * tau (2, na) + &
                q_(3) * tau (3, na) ) * tpi
           phase = CMPLX (cos (arg), - sin (arg) )
           do ig = 1, npw_
              sk (ig) = eigts1 (ig1(igk_(ig)), na) * &
                   eigts2 (ig2(igk_(ig)), na) * &
                        eigts3 (ig3(igk_(ig)), na)
           end do
           do ih = 1, paw_recon(nt)%paw_nh
              jkb = jkb + 1
              pref = (0.d0, -1.d0) ** paw_recon(nt)%paw_nhtol(ih) * phase
              do ig = 1, npw_
                 vkb_(ig, jkb) = vkb1 (ig,ih) * sk (ig) * pref
              end do
              
           end do
        end if
        
     end do
     
     deallocate (vkb1)
  end do
  
  deallocate (gk)
  deallocate (ylm)
  deallocate (vq)
  deallocate (qg)
  deallocate (sk)

  call stop_clock ('init_paw_2')
  return
end subroutine init_paw_2

