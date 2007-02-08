!
! Copyright (C) 2001 PWSCF group
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
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat, ityp, tau
  use pwcom
  USE kinds, only : DP
  use phcom
  implicit none

  complex(DP) ::  dynwrk (3 * nat, 3 * nat),  wrk, exc
  complex(DP), allocatable :: vxc (:), work(:)
  ! exchange and correlation potential
  real(DP), allocatable :: v (:,:)
  real(DP) :: q0(3), arg, etxcd, vtxcd
  !
  integer :: i, j, ir, is, isup, isdw, ig, na, nb, nta, ntb, na_i, na_j, &
       nb_j, na_icart, nb_jcart, nu_i, nu_j
  !
  ! return immediately if nlcc is not present
  !
  if (.not.nlcc_any) return
  call start_clock ('dynmatcc')
  !
  ! allocate workspace
  !
  allocate (vxc( nrxx))    
  allocate (v  ( nrxx , nspin))    
  !
  call v_xc (rho, rhog, rho_core, rhog_core, etxcd, vtxcd, v)
  !
  if (nspin == 1 .OR. nspin==4) then
     is=1
     do ir = 1, nrxx
        vxc(ir) = v(ir,is)
     end do
  else
     isup=1
     isdw=2
     do ir = 1, nrxx
        vxc (ir) = (v(ir,isup) + v(ir,isdw))*0.5d0
     end do
  end if
  deallocate (v)
  !
  call cft3 (vxc, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  !
  ! vxc is the spin-averaged XC potential (in G-space)
  !
  q0 = 0.d0
  call set_drhoc (q0)
  !
  ! set_drhoc produces drc=Drho_core(G)/DG , without struct.fact.
  !
  dynwrk (:,:) = (0.d0, 0.d0)
  allocate (work (nrxx))    
  do na = 1, nat
     nta = ityp (na)
     work (:) = (0.d0, 0.d0)
     do ig = 1, ngm
        arg = tpi * (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) &
             + g (3, ig) * tau (3, na) )
        exc = CMPLX (cos (arg), - sin (arg) ) * tpiba2
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
#ifdef __PARA
  call reduce (18 * nat * nat, dynwrk)
#endif
  !
  dynwrk = dynwrk * omega
  !
  ! calculate drc for later use in calculation of non diagonal term
  !
  call set_drhoc (xq)
  !
  ! rotate in the pattern basis and add to dynmat
  !
  do nu_i = 1, 3 * nat
     do nu_j = 1, 3 * nat
        wrk = (0.d0, 0.d0)
        do nb_jcart = 1, 3 * nat
           do na_icart = 1, 3 * nat
              wrk = wrk + CONJG(u (na_icart, nu_i) ) * dynwrk (na_icart, &
                   nb_jcart) * u (nb_jcart, nu_j)
           enddo
        enddo
        dyn (nu_i, nu_j) = dyn (nu_i, nu_j) + wrk
     enddo
  enddo
  !
  deallocate (work)
  deallocate (vxc)
  call stop_clock ('dynmatcc')
  return
end subroutine dynmatcc
