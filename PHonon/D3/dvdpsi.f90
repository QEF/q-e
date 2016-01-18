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
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp
  USE cell_base,  ONLY : tpiba
  USE fft_base,   ONLY : dffts, dfftp
  USE fft_interfaces,  ONLY : fwfft, invfft
  USE gvect,      ONLY : g
  USE gvecs,    ONLY : nls
  USE wvfct,      ONLY : nbnd, npwx, npw, igk
  use qpoint,     ONLY : npwq, igkq
  use phcom
  use d3com
  USE uspp,       ONLY : nkb, dvan
  USE uspp_param, ONLY : nh
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum

!
  implicit none
  integer :: nu_i
  ! input: the mode under consideration
  real (DP) :: xq_ (3)
  ! input: coordinates of the q point describing the perturbation
  complex (DP) :: dvloc (dfftp%nnr), psi_ (npwx, nbnd), dvpsi_ (npwx, nbnd)
  ! input: local part of the KS potential
  ! input: wavefunction
  ! output: variation of the KS potential applied to psi_
  complex(DP) :: vkb_(npwx,nkb), vkbq_(npwx,nkb)
  !
  ! Local variables
  !
  integer :: na, mu, ig, igg, ir, ibnd, nt, ikb, jkb
  ! counters
  complex (DP), pointer :: u_x (:,:)
  ! the transformation modes patterns
  complex (DP), allocatable :: aux (:), ps (:,:), wrk2 (:)
  ! work space
  complex (DP) , external:: zdotc
  logical :: q_eq_zero
  !
  allocate  (aux( dfftp%nnr))
  allocate  (ps( 2, nbnd))
  allocate  (wrk2( npwx))
  q_eq_zero = xq_ (1) == 0.d0 .and. xq_ (2) == 0.d0 .and. xq_ (3) == 0.d0
  if (q_eq_zero) then
     u_x => ug0
  else
     u_x => u
  endif
  !
  do ibnd = 1, nbnd
     aux (:) = (0.d0, 0.d0)
     do ig = 1, npw
        aux (nls (igk (ig) ) ) = psi_ (ig, ibnd)
     enddo
     CALL invfft ('Wave', aux, dffts)
     do ir = 1, dffts%nnr
        aux (ir) = aux (ir) * dvloc (ir)
     enddo
     CALL fwfft ('Wave', aux, dffts)
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
                  abs (u_x (mu + 3, nu_i) ) > 1.0d-12) then
           !
           ! first term: sum_l v_l beta_l(k+q+G) \sum_G' beta^*_l(k+G') (iG'*u) psi
           ! second term: sum_l E_l(-i(q+G)*u) beta_l(k+q+G)\sum_G'beta^*_l(k+G')ps
           !
                 do ig = 1, npw
                    wrk2 (ig) = vkb_(ig,jkb) * &
                         CONJG(CMPLX(0.d0,1.d0,kind=DP) *tpiba * &
                                (g (1, igk (ig) ) * u_x (mu + 1, nu_i) + &
                                 g (2, igk (ig) ) * u_x (mu + 2, nu_i) + &
                                 g (3, igk (ig) ) * u_x (mu + 3, nu_i) ) )
                 enddo
                 do ibnd = 1, nbnd
                    ps(1,ibnd) = dvan(ikb,ikb,nt) * &
                         zdotc(npw, wrk2, 1, psi_(1,ibnd), 1)
                    ps(2,ibnd) = dvan(ikb,ikb,nt) * &
                         zdotc(npw,vkb_(1,jkb),1,psi_(1,ibnd),1)
                 enddo
                 !
                 ! when build is serial this call does nothing, we leave it there
                 !
                 call mp_sum ( ps, intra_pool_comm )

                 do ig = 1, npwq
                    wrk2 (ig) = vkbq_(ig,jkb) * CMPLX(0.d0,-1.d0,kind=DP) * tpiba * &
                         ( (g (1, igkq (ig) ) + xq_ (1) ) * u_x (mu+1, nu_i) +&
                           (g (2, igkq (ig) ) + xq_ (2) ) * u_x (mu+2, nu_i) +&
                           (g (3, igkq (ig) ) + xq_ (3) ) * u_x (mu+3, nu_i) )
                 enddo

                 do ibnd = 1, nbnd
                    call zaxpy(npwq,ps(1,ibnd),vkbq_(1,jkb),1,dvpsi_(1,ibnd),1)
                    call zaxpy(npwq,ps(2,ibnd),       wrk2, 1,dvpsi_(1,ibnd),1)
                 enddo
              endif
           enddo
        end if
     end do
  end do
  deallocate (wrk2)
  deallocate (ps)
  deallocate (aux)
  return
end subroutine dvdpsi
