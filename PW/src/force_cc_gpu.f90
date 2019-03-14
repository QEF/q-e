!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine force_cc_gpu (forcecc)
  !----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE atom,                 ONLY : rgrid, msh
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : alat, omega, tpiba, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, g, gg, ngl, gl, igtongl, igtongl_d
  USE ener,                 ONLY : etxc, vtxc
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin
  USE wavefunctions, ONLY : psic
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  implicit none
  !
  !   first the dummy variable
  !
  real(DP) :: forcecc (3, nat)
  ! output: the local forces on atoms

  integer :: ipol, ig, ir, nt, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on FFT grid points
  ! counter on types of atoms
  ! counter on atoms


  real(DP), allocatable :: vxc (:,:), rhocg (:)
  ! exchange-correlation potential
  ! radial fourier transform of rho core
  real(DP)  ::  arg, fact

  !
  forcecc(:,:) = 0.d0
  if ( ANY ( upf(1:ntyp)%nlcc ) ) go to 15
  return
  !
15 continue
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  !
  ! recalculate the exchange-correlation potential
  !
  allocate ( vxc(dfftp%nnr,nspin) )
  !
  call v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxc)
  !
  psic=(0.0_DP,0.0_DP)
  if (nspin == 1 .or. nspin == 4) then
     do ir = 1, dfftp%nnr
        psic (ir) = vxc (ir, 1)
     enddo
  else
     do ir = 1, dfftp%nnr
        psic (ir) = 0.5d0 * (vxc (ir, 1) + vxc (ir, 2) )
     enddo
  endif
  deallocate (vxc)
  CALL fwfft ('Rho', psic, dfftp)
  !
  ! psic contains now Vxc(G)
  !
  allocate ( rhocg(ngl) )
  !
  ! core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
  ! g = 0 term gives no contribution
  !
  do nt = 1, ntyp
     if ( upf(nt)%nlcc ) then

        call drhoc (ngl, gl, omega, tpiba2, msh(nt), rgrid(nt)%r,&
             rgrid(nt)%rab, upf(nt)%rho_atc, rhocg)

        rhocg_d = rhocg
        do na = 1, nat
           if (nt.eq.ityp (na) ) then

              tau1 = tau(1, na)
              tau2 = tau(2, na)
              tau3 = tau(3, na)

              fcc1 = 0.d0
              fcc2 = 0.d0
              fcc3 = 0.d0

              !$cuf kernel do (1) <<<*, *>>>
              do ig = gstart, ngm
                 arg = (g_d (1, ig) * tau1 + g_d (2, ig) * tau2 &
                      + g_d (3, ig) * tau3 ) * tpi
                 prod = tpiba * omega * &
                      rhocg_d (igtongl_d (ig) ) * dble(CONJG(psic_d (nl_d (ig) ) ) * &
                      CMPLX( sin (arg), cos (arg), kind=DP))* fact

                 fcc1 = fcc1 + g_d (1, ig) * prod
                 fcc2 = fcc2 + g_d (2, ig) * prod
                 fcc3 = fcc3 + g_d (3, ig) * prod
              enddo

              forcecc(1, na) = forcecc(1, na) + fcc1
              forcecc(2, na) = forcecc(2, na) + fcc2
              forcecc(3, na) = forcecc(3, na) + fcc3
           endif
        enddo
     endif
  enddo

  !
  call mp_sum(  forcecc, intra_bgrp_comm )
  !
  deallocate (rhocg)
  !
  return
end subroutine force_cc_gpu
