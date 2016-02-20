!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
subroutine dynmatcc
  !--------------------------------------------------------------------
  !
  ! diagonal (q-independent) NLCC contribution to the dynamical matrix
  !
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba2, omega
  USE ions_base, ONLY : nat, ityp, tau
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE gvect,     ONLY : nl, ngm, g
  USE lsda_mod,  ONLY : nspin
  use scf,       ONLY : rho, rho_core, rhog_core
  USE modes,     ONLY : u
  USE qpoint,    ONLY : xq
  USE nlcc_ph,   ONLY : drc
  USE uspp,      ONLY : nlcc_any
  USE dynmat,    ONLY : dyn
  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum

  implicit none

  complex(DP) ::  dynwrk (3 * nat, 3 * nat),  wrk, exc
  complex(DP), allocatable :: vxc (:), work(:)
  ! exchange and correlation potential
  real(DP), allocatable :: v (:,:)
  real(DP) :: q0(3), arg, etxcd, vtxcd
  !
  integer :: i, j, ir, is, isup, isdw, ig, na, nta, na_i, na_j, &
       na_icart, nb_jcart, nu_i, nu_j
  !
  ! return immediately if nlcc is not present
  !
  if (.not.nlcc_any) return
  call start_clock ('dynmatcc')
  !
  ! allocate workspace
  !
  allocate (vxc( dfftp%nnr))
  allocate (v  ( dfftp%nnr , nspin))
  !
  call v_xc (rho, rho_core, rhog_core, etxcd, vtxcd, v)
  !
  if (nspin == 1 .OR. nspin==4) then
     is=1
     do ir = 1, dfftp%nnr
        vxc(ir) = v(ir,is)
     end do
  else
     isup=1
     isdw=2
     do ir = 1, dfftp%nnr
        vxc (ir) = (v(ir,isup) + v(ir,isdw))*0.5d0
     end do
  end if
  deallocate (v)
  !
  CALL fwfft ('Dense', vxc, dfftp)
  !
  ! vxc is the spin-averaged XC potential (in G-space)
  !
  q0 = 0.d0
  call set_drhoc (q0, drc)
  !
  ! set_drhoc produces drc=Drho_core(G)/DG , without struct.fact.
  !
  dynwrk (:,:) = (0.d0, 0.d0)
  allocate (work (dfftp%nnr))
  do na = 1, nat
     nta = ityp (na)
     work (:) = (0.d0, 0.d0)
     do ig = 1, ngm
        arg = tpi * (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) &
             + g (3, ig) * tau (3, na) )
        exc = CMPLX(cos (arg), - sin (arg) ,kind=DP) * tpiba2
        work (ig) = drc (ig, nta) * exc * CONJG(vxc (nl (ig) ) )
     enddo
     do i = 1, 3
        na_i = 3 * (na - 1) + i
        do j = 1, 3
           na_j = 3 * (na - 1) + j
           do ig = 1, ngm
              dynwrk (na_i, na_j) = dynwrk (na_i, na_j) - g(i, ig) * g(j, ig) &
                   * work (ig)
           enddo
        enddo
     enddo
  enddo
  call mp_sum (dynwrk,intra_bgrp_comm)
  !
  dynwrk = dynwrk * omega
  !
  ! calculate drc for later use in calculation of non diagonal term
  !
  call set_drhoc (xq, drc)
  !
  ! rotate in the pattern basis and add to dynmat
  !
  CALL rotate_pattern_add(nat, u, dyn, dynwrk)
  !
  deallocate (work)
  deallocate (vxc)
  call stop_clock ('dynmatcc')
  return
end subroutine dynmatcc
