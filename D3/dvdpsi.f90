!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dvdpsi (nu_i, xq_, dvloc, vkb_, vkbq_, psi_, dvpsi_)  
!-----------------------------------------------------------------------
!
! Receives in input the variation of the local part of the KS-potential
! and calculates dV(xq_)_KS*psi_ in G_space, for all bands
!
#include "machine.h"
  use pwcom
  use phcom
  use d3com
  use allocate
!
  implicit none
  integer :: nu_i  
  ! input: the mode under consideration
  real (8) :: xq_ (3)  
  ! input: coordinates of the q point describing the perturbation
  complex (8) :: dvloc (nrxx), psi_ (npwx, nbnd), dvpsi_ (npwx, nbnd)
  ! input: local part of the KS potential
  ! input: wavefunction
  ! output: variation of the KS potential applied to psi_
  complex(8) :: vkb_(npwx,nkb), vkbq_(npwx,nkb)
  !
  ! Local variables
  !
  integer :: na, mu, ig, igg, ir, ibnd, nt, ikb, jkb  
  ! counter on atoms
  ! counter on modes
  ! counter on G vectors
  ! counter on G vectors
  ! counter on real space points
  ! counter on bands
  ! counter on atomic types
  ! counters on beta functions
  complex (8), pointer :: u_x (:,:), aux (:), ps (:,:), wrk2 (:)
  ! the transformation modes patterns
  ! work space
  complex (8) :: ZDOTC
  logical :: q_eq_zero
  !
  call mallocate (aux, nrxx)  
  call mallocate (ps, 2, nbnd)  
  call mallocate (wrk2, npwx)
  q_eq_zero = xq_ (1) .eq.0.d0.and.xq_ (2) .eq.0.d0.and.xq_ (3) .eq.0.d0
  if (q_eq_zero) then  
     u_x => ug0  
  else  
     u_x => u  
  endif
  !
  do ibnd = 1, nbnd  
     call setv (2 * nrxxs, 0.d0, aux, 1)  
     do ig = 1, npw  
        aux (nls (igk (ig) ) ) = psi_ (ig, ibnd)  
     enddo
     call cft3s (aux, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)  
     do ir = 1, nrxxs  
        aux (ir) = aux (ir) * dvloc (ir)  
     enddo
     call cft3s (aux, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)  
     do ig = 1, npwq  
        dvpsi_ (ig, ibnd) = aux (nls (igkq (ig) ) )  
     enddo
  enddo
  !
  !    Now the contribution of the non local part in the KB form
  !
  jkb=0
  do nt = 1, ntyp  
     do na = 1, nat  
        if (ityp (na).eq.nt) then  
           mu = 3 * (na - 1)  
           do ikb = 1, nh (nt)  
              jkb = jkb+1
              if (abs (u_x (mu + 1, nu_i) ) + abs (u_x (mu + 2, nu_i) ) + &
                  abs (u_x (mu + 3, nu_i) ) .gt.1.0d-12) then
           !
           ! first term: sum_l v_l beta_l(k+q+G) \sum_G' beta^*_l(k+G') (iG'*u) psi
           ! second term: sum_l E_l(-i(q+G)*u) beta_l(k+q+G)\sum_G'beta^*_l(k+G')ps
           !
                 do ig = 1, npw  
                    wrk2 (ig) = vkb_(ig,jkb) * &
                         conjg(DCMPLX(0.d0,1.d0) *tpiba * &
                                (g (1, igk (ig) ) * u_x (mu + 1, nu_i) + &
                                 g (2, igk (ig) ) * u_x (mu + 2, nu_i) + &
                                 g (3, igk (ig) ) * u_x (mu + 3, nu_i) ) )
                 enddo
                 do ibnd = 1, nbnd  
                    ps(1,ibnd) = dvan(ikb,ikb,nt) * &
                         ZDOTC(npw, wrk2, 1, psi_(1,ibnd), 1)
                    ps(2,ibnd) = dvan(ikb,ikb,nt) * &
                         ZDOTC(npw,vkb_(1,jkb),1,psi_(1,ibnd),1)
                 enddo
#ifdef PARA
                 call reduce (4 * nbnd, ps)  
#endif
                 do ig = 1, npwq  
                    wrk2 (ig) = vkbq_(ig,jkb) * DCMPLX(0.d0,-1.d0) * tpiba * &
                         ( (g (1, igkq (ig) ) + xq_ (1) ) * u_x (mu+1, nu_i) +&
                           (g (2, igkq (ig) ) + xq_ (2) ) * u_x (mu+2, nu_i) +&
                           (g (3, igkq (ig) ) + xq_ (3) ) * u_x (mu+3, nu_i) )
                 enddo

                 do ibnd = 1, nbnd  
                    call ZAXPY(npwq,ps(1,ibnd),vkbq_(1,jkb),1,dvpsi_(1,ibnd),1)
                    call ZAXPY(npwq,ps(2,ibnd),       wrk2, 1,dvpsi_(1,ibnd),1)
                 enddo
              endif
           enddo
        end if
     end do
  end do
  call mfree (wrk2)  
  call mfree (ps)  
  call mfree (aux)  
  return  
end subroutine dvdpsi
