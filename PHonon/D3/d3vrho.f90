!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine d3vrho()
  !-----------------------------------------------------------------------
  !
  !  This routine calculates the electronic term: <psi|V"'|psi>
  !  of the third order dynamical matrix.
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp, tau
  USE uspp,                 ONLY : dvan
  USE scf,                  ONLY : rho
  USE gvect,                ONLY : g, ngm, nl, igtongl
  USE wvfct,                ONLY : npw, npwx, nbnd, igk, wg
  USE vlocal,               ONLY : vloc
  USE klist,                ONLY : xk
  USE cell_base,            ONLY : omega, tpiba, tpiba2
  USE uspp_param,           ONLY : nh
  USE wavefunctions_module, ONLY : evc
  USE io_files,             ONLY : iunigk
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE phcom
  USE d3com

  USE qpoint,               ONLY : nksq, npwq, igkq
  USE control_lr,           ONLY : nbnd_occ, lgamma
  !
  implicit none
  integer :: icart, jcart, kcart, na_i, na_j, na_k, na, ng, ir, nt, &
       ik, ikk, ig, ibnd, ikb, jkb, ios, igg, ia
  ! counters

  real (DP) :: gtau, fac, wgg
  ! the product G*\tau_s
  ! auxiliary variable
  ! the true weight of a K point

  complex (DP) :: alpha (8), zdotc, work
  complex (DP), allocatable :: d3dynwrk (:,:,:), d3dynwrk2 (:,:,:), &
       rhog (:), work1 (:,:), work2 (:,:), work3 (:)

  allocate  (rhog( dfftp%nnr))
  allocate  (d3dynwrk( 3 * nat, 3 * nat, 3 * nat))
  allocate  (d3dynwrk2(3 * nat, 3 * nat, 3 * nat))
  allocate  (work1(  npwx, 3))
  allocate  (work2(  npwx, 3))
  allocate  (work3(  npwx))

  d3dynwrk (:,:,:) = (0.d0, 0.d0)
  do ir = 1, dfftp%nnr
     rhog (ir) = CMPLX(rho%of_r (ir, 1), 0.d0,kind=DP)
  enddo
  CALL fwfft ('Dense', rhog, dfftp)
  !
  !     Contribution deriving from the local part of the potential
  !
  do na_i = npert_i, npert_f
     na = (na_i - 1) / 3 + 1
     icart = na_i - 3 * (na - 1)
     do jcart = 1, 3
        na_j = 3 * (na - 1) + jcart
        do kcart = 1, 3
           na_k = 3 * (na - 1) + kcart
           do ng = 1, ngm
              gtau = tpi * (g (1, ng) * tau (1, na) + g (2, ng) * tau (2, na) &
                   + g (3, ng) * tau (3, na) )
              fac = vloc (igtongl (ng), ityp (na) ) * tpiba2 * tpiba * omega *&
                   (DBLE (rhog (nl (ng) ) ) * sin (gtau) + &
                   AIMAG (rhog (nl (ng) ) ) * cos (gtau) )
              d3dynwrk (na_i, na_j, na_k) = d3dynwrk (na_i, na_j, na_k) + &
                   fac * g (icart, ng) * g (jcart, ng) * g (kcart, ng)
           enddo
        enddo
     enddo
  enddo
#ifdef __MPI
  call mp_sum ( d3dynwrk, intra_pool_comm )
#endif
  !
  !     Non local Kleinman-Bylander potential contribution
  !
  rewind (unit = iunigk)

  do ik = 1, nksq
     read (iunigk, err = 100, iostat = ios) npw, igk
     if (lgamma) then
        ikk = ik
     else
        read (iunigk, err = 200, iostat = ios) npwq, igkq
        ikk = 2 * ik - 1
     endif
100  call errore ('d3vrho', 'reading igk', abs (ios) )
200  call errore ('d3vrho', 'reading igkq', abs (ios) )
     call davcio (evc, lrwfc, iuwfc, ikk, - 1)
     call init_us_2 (npw, igk, xk (1, ikk), vkb0)

     do kcart = 1, 3
        do icart = 1, 3
           do jcart = 1, 3
              do ibnd = 1, nbnd_occ (ikk)
                 wgg = wg (ibnd, ikk)
                 do ig = 1, npw
                    work3 (ig) = evc (ig, ibnd) * tpiba * g (icart, igk (ig) )&
                         * tpiba * g (jcart, igk (ig) ) * tpiba * g (kcart, igk (ig) )
                    work2 (ig, 1) = evc (ig, ibnd) * tpiba * g (icart, igk (ig) ) &
                         * tpiba * g (jcart, igk (ig) )
                    work2 (ig, 2) = evc (ig, ibnd) * tpiba * g (jcart, igk (ig) ) &
                         * tpiba * g (kcart, igk (ig) )
                    work2 (ig, 3) = evc (ig, ibnd) * tpiba * g (kcart, igk (ig) ) &
                         * tpiba * g (icart, igk (ig) )
                    work1 (ig, 1) = evc (ig, ibnd) * tpiba * g (kcart, igk (ig) )
                    work1 (ig, 2) = evc (ig, ibnd) * tpiba * g (icart, igk (ig) )
                    work1 (ig, 3) = evc (ig, ibnd) * tpiba * g (jcart, igk (ig) )
                 enddo
                 jkb=0
                 do nt = 1, ntyp
                    do na = 1, nat
                       if (ityp (na) == nt) then
                          na_k = 3 * (na - 1) + kcart
                          na_i = 3 * (na - 1) + icart
                          na_j = 3 * (na - 1) + jcart
                          do ikb = 1, nh (nt)
                             jkb=jkb+1
                             alpha (1) = zdotc (npw, work3, 1, vkb0(1,jkb), 1)
                             alpha (2) = zdotc (npw, vkb0(1,jkb), 1, evc (1, ibnd), 1)
                             alpha (3) = zdotc (npw,work1(1, 1),1,vkb0(1,jkb),1)
                             alpha (4) = zdotc (npw,vkb0(1,jkb),1,work2(1, 1),1)
                             alpha (5) = zdotc (npw,work1(1, 2),1,vkb0(1,jkb),1)
                             alpha (6) = zdotc (npw,vkb0(1,jkb),1,work2(1, 2),1)
                             alpha (7) = zdotc (npw,work1(1, 3),1,vkb0(1,jkb),1)
                             alpha (8) = zdotc (npw,vkb0(1,jkb),1,work2(1, 3),1)
#ifdef __MPI
                             call mp_sum ( alpha, intra_pool_comm )
#endif
                             d3dynwrk (na_k, na_i, na_j) = d3dynwrk (na_k, na_i, na_j) - &
                                  2.0d0 * dvan(ikb,ikb,nt) * wgg * &
                                  AIMAG(alpha(1)*alpha(2) + alpha(3)*alpha(4) +&
                                        alpha(5)*alpha(6) + alpha(7)*alpha(8))
                          enddo
                       endif
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
#ifdef __MPI
  call mp_sum( d3dynwrk, inter_pool_comm )
#endif
  !
  !   The dynamical matrix was computed in cartesian axis and now we put
  !   it on the basis of the modes
  !
  d3dynwrk2(:,:,:) = (0.d0, 0.d0)
  do na_k = npert_i, npert_f
     if (q0mode (na_k) ) then
        do na_i = 1, 3 * nat
           do na_j = 1, 3 * nat
              work = (0.d0, 0.d0)
              do kcart = 1, 3 * nat
                 do icart = 1, 3 * nat
                    do jcart = 1, 3 * nat
                       work = work + ug0 (kcart, na_k) * CONJG(u (icart, na_i) ) &
                            * d3dynwrk (kcart, icart, jcart) * u (jcart, na_j)
                    enddo
                 enddo
              enddo
              d3dynwrk2 (na_k, na_i, na_j) = work
           enddo
        enddo
     endif
  enddo
#ifdef __MPI
  call mp_sum( d3dynwrk2, inter_pool_comm )
#endif
  d3dyn (:,:,:) = d3dyn (:,:,:) +  d3dynwrk2 (:,:,:)
  d3dyn_aux1(:,:,:) = d3dynwrk2 (:,:,:)

  deallocate (work1)
  deallocate (work2)
  deallocate (work3)
  deallocate (d3dynwrk2)
  deallocate (d3dynwrk)
  deallocate (rhog)

  return
end subroutine d3vrho
