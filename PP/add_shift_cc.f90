!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine add_shift_cc (shift_cc)
  !----------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds, ONLY : DP
  USE constants, ONLY: tpi
  USE atom, ONLY: rho_atc, numeric, mesh, r, rab, nlcc
  USE ions_base, ONLY: nat, ntyp => nsp, ityp, tau
  USE cell_base, ONLY: alat, omega, tpiba, tpiba2
  USE gvect, ONLY: ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
       nrxx, nl, g, gg, ngl, gl, igtongl
  USE ener, ONLY: etxc, vtxc
  USE lsda_mod, ONLY: nspin
  USE pseud, ONLY: a_nlcc, b_nlcc, alpha_nlcc
  USE scf, ONLY: rho, rhog, rho_core, rhog_core
  USE wvfct, ONLY: gamma_only
  USE wavefunctions_module,    ONLY : psic
  implicit none
  !
  !   first the dummy variable
  !

  real(DP) :: shift_cc (nat)
  ! output: the local forces on atoms

  integer :: ig, ir, nt, na
  ! counter on G vectors
  ! counter on FFT grid points
  ! counter on types of atoms
  ! counter on atoms


  real(DP), allocatable :: vxc (:,:), rhocg (:), shift_(:)
  ! exchange-correlation potential
  ! radial fourier trasform of rho core
  real(DP)  ::  arg, fact
  !
  do nt = 1, ntyp
     if (nlcc (nt) ) goto 15
  enddo
  return
  !
15 continue
   call infomsg ('add_shift_cc','BEWARE: shift with CC never tested !!!', -1 )
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  !
  ! recalculate the exchange-correlation potential
  !
  allocate ( vxc(nrxx,nspin), shift_(nat) )
  shift_(:) = 0.d0
  !
  call v_xc (rho, rhog, rho_core, rhog_core, etxc, vtxc, vxc)
  !
  if (nspin.eq.1) then
     do ir = 1, nrxx
        psic (ir) = vxc (ir, 1)
     enddo
  else
     do ir = 1, nrxx
        psic (ir) = 0.5d0 * (vxc (ir, 1) + vxc (ir, 2) )
     enddo
  endif
  deallocate (vxc)
  call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  !
  ! psic contains now Vxc(G)
  !
  allocate ( rhocg(ngl) )
  !
  ! core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
  ! g = 0 term gives no contribution
  !
  do nt = 1, ntyp
     if (nlcc (nt) ) then

        call drhoc (ngl, gl, omega, tpiba2, numeric (nt), a_nlcc (nt), &
             b_nlcc (nt), alpha_nlcc (nt), mesh (nt), r (1, nt), rab (1, nt), &
             rho_atc (1, nt), rhocg)
        do na = 1, nat
           if (nt.eq.ityp (na) ) then
              if (gstart.eq.2)  shift_(na) = omega * rhocg (igtongl (1) ) * &
                                                     CONJG(psic (nl (1) ) )
              do ig = gstart, ngm
                 arg = (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) &
                      + g (3, ig) * tau (3, na) ) * tpi
                 shift_ (na) = shift_( na) + omega * &
                         rhocg (igtongl (ig) ) * CONJG(psic (nl (ig) ) ) * &
                         CMPLX ( cos(arg), -sin(arg)) * fact
              enddo
           endif
        enddo
     endif
  enddo
#ifdef __PARA
  call reduce (3 * nat, shift_)
#endif
  shift_cc(:) = shift_cc(:) + shift_(:)
  deallocate (rhocg, shift_)
  !
  return
end subroutine add_shift_cc

