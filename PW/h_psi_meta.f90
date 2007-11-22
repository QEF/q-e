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
subroutine h_psi_meta (ldap, np, mp, psip, hpsi)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the Hubbard potential applied to the electronic
  ! of the current k-point, the result is added to hpsi
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba
  USE lsda_mod,  ONLY : nspin, current_spin
  USE wvfct,     ONLY : igk, current_k
  USE gsmooth,   ONLY : nls, nlsm, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE gvect,     ONLY : g
  USE scf,       ONLY : kedtau
  USE klist,     ONLY : xk
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  !
  implicit none
  !
  COMPLEX(DP), PARAMETER :: ci=(0.d0,1.d0)
  integer :: ldap, np, mp
  complex(DP) :: psip (ldap, mp), hpsi (ldap, mp)
  real (DP), allocatable :: kplusg (:)
!  complex (DP), allocatable :: psi(:)
  !
  integer :: im, j,ir
  !
  CALL start_clock( 'h_psi_meta' )
  allocate (kplusg(np))
  if (gamma_only) then
     !
     ! gamma algorithm
     !
     do im = 1, mp, 2
        do j =1,3
           psic(1:nrxxs) = ( 0.D0, 0.D0 )
           !
           kplusg (1:np) = (xk(j,current_k)+g(j,igk(1:np))) * tpiba
           if (im < mp ) then
              psic(nls(igk(1:np))) =  ci * kplusg(1:np) * &
                              ( psip (1:np,im) + ci * psip(1:np,im+1) )
              psic(nlsm(igk(1:np)))= -ci * kplusg(1:np) * &
                        CONJG ( psip (1:np,im) - ci * psip(1:np,im+1) )
           else
              psic(nls(igk(1:np))) =  ci * kplusg(1:np) *       psip(1:np,im) 
              psic(nlsm(igk(1:np)))= -ci * kplusg(1:np) * CONJG(psip(1:np,im))
           end if
           !
           CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2 )
           !
           psic(1:nrxxs) = kedtau(1:nrxxs,current_spin) * psic(1:nrxxs) 
           !
           CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2 )
           !
           if ( im < mp ) then
              hpsi(1:np,im)  = hpsi(1:np,im)   - ci * kplusg(1:np) * 0.5d0 * &
                       ( psic(nls(igk(1:np))) + CONJG(psic(nlsm(igk(1:np)))) )
              hpsi(1:np,im+1)= hpsi(1:np,im+1) - kplusg(1:np) * 0.5d0 * &
                       ( psic(nls(igk(1:np))) - CONJG(psic(nlsm(igk(1:np)))) )
           else
              hpsi(1:np,im) = hpsi(1:np,im) - ci * kplusg(1:np) * &
                                                      psic(nls(igk(1:np)))
           end if
        end do
     end do
  else
     !
     ! generic k algorithm
     !
     do im = 1, mp
        do j =1,3
           psic(1:nrxxs) = ( 0.D0, 0.D0 )
           !
           kplusg (1:np) = (xk(j,current_k)+g(j,igk(1:np))) * tpiba
           psic(nls(igk(1:np))) = CMPLX (0d0, kplusg(1:np)) * psip (1:np,im)
           !
           CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2 )
           !
           psic(1:nrxxs) = kedtau(1:nrxxs,current_spin) * psic(1:nrxxs) 
           !
           CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2 )
           !
           hpsi(1:np,im) = hpsi(1:np,im) - &
                           CMPLX (0d0, kplusg(1:np)) * psic(nls(igk(1:np)))
        end do
     end do
  end if
  deallocate (kplusg)
  CALL stop_clock( 'h_psi_meta' )


  return

end subroutine h_psi_meta
