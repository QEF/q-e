!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine h_psi_meta (ldap, np, mp, psip, hpsi)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the specific contribution from the meta-GGA
  ! potential to H*psi; the result is added to hpsi
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba
  USE lsda_mod,  ONLY : nspin, current_spin
  USE wvfct,     ONLY : current_k
  USE gvecs,     ONLY : nls, nlsm
  USE gvect,     ONLY : g
  USE scf,       ONLY : kedtau
  USE klist,     ONLY : xk, igk_k
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  !
  implicit none
  !
  COMPLEX(DP), PARAMETER :: ci=(0.d0,1.d0)
  integer :: ldap, np, mp
  complex(DP) :: psip (ldap, mp), hpsi (ldap, mp)
  real (DP), allocatable :: kplusg (:)
  !
  integer :: im, j, nrxxs
  !
  CALL start_clock( 'h_psi_meta' )
  nrxxs = dffts%nnr
  allocate (kplusg(np))
  if (gamma_only) then
     !
     ! gamma algorithm
     !
     do im = 1, mp, 2
        do j =1,3
           psic(1:nrxxs) = ( 0.D0, 0.D0 )
           !
           kplusg (1:np) = (xk(j,current_k)+g(j,1:np)) * tpiba
           if (im < mp ) then
              psic(nls (1:np)) =  ci * kplusg(1:np) * &
                              ( psip (1:np,im) + ci * psip(1:np,im+1) )
              psic(nlsm(1:np)) = -ci * kplusg(1:np) * &
                        CONJG ( psip (1:np,im) - ci * psip(1:np,im+1) )
           else
              psic(nls (1:np)) =  ci * kplusg(1:np) *       psip(1:np,im) 
              psic(nlsm(1:np)) = -ci * kplusg(1:np) * CONJG(psip(1:np,im))
           end if
           !
           CALL invfft ('Wave', psic, dffts)
           !
           psic(1:nrxxs) = kedtau(1:nrxxs,current_spin) * psic(1:nrxxs) 
           !
           CALL fwfft ('Wave', psic, dffts)
           !
           if ( im < mp ) then
              hpsi(1:np,im)  = hpsi(1:np,im)   - ci * kplusg(1:np) * 0.5d0 * &
                       ( psic(nls(1:np)) + CONJG(psic(nlsm(1:np))) )
              hpsi(1:np,im+1)= hpsi(1:np,im+1) - kplusg(1:np) * 0.5d0 * &
                       ( psic(nls(1:np)) - CONJG(psic(nlsm(1:np))) )
           else
              hpsi(1:np,im) = hpsi(1:np,im) - ci * kplusg(1:np) * &
                                              psic(nls(1:np))
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
           kplusg (1:np) = (xk(j,current_k)+g(j,igk_k(1:np,current_k)))*tpiba
           psic(nls(igk_k(1:np,current_k))) = CMPLX(0d0, kplusg(1:np),kind=DP)&
                                            * psip (1:np,im)
           !
           CALL invfft ('Wave', psic, dffts)
           !
           psic(1:nrxxs) = kedtau(1:nrxxs,current_spin) * psic(1:nrxxs) 
           !
           CALL fwfft ('Wave', psic, dffts)
           !
           hpsi(1:np,im) = hpsi(1:np,im) - CMPLX(0d0, kplusg(1:np),kind=DP) &
                                         * psic(nls(igk_k(1:np,current_k)))
        end do
     end do
  end if
  deallocate (kplusg)
  CALL stop_clock( 'h_psi_meta' )


  return

end subroutine h_psi_meta
