!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine d3dyn_cc
  !-----------------------------------------------------------------------
  !
  ! It calculates contribution due to non-linear-core-correction
  ! The variation of the density with respect to the perturbation must
  ! be corrected before calling this routine:
  ! while reading the variation of the density on unit iudrho and iud0rho
  ! it assumes it is the total density, i.e. sum of valence + core.
  !
  USE ions_base,  ONLY : nat, ityp, tau
  USE kinds, only : DP
  USE funct, only : xc
  USE fft_base,              ONLY : dfftp
  USE fft_interfaces,        ONLY : fwfft
  use pwcom
  use scf, only : rho, rho_core
  use qpoint, ONLY : xq
  use phcom
  use d3com
  USE mp_global,  ONLY : inter_pool_comm, intra_pool_comm
  USE mp,         ONLY : mp_sum

  implicit none
  integer :: na, nta, ig, ir, i_cart, j_cart, k_cart, na_i, na_j, &
       na_k, nu_i, nu_j, nu_k, na_icart, nb_jcart, nc_kcart

  real (DP) :: rhox, arhox, ex, ec, vx, vc, arg
  ! the total charge in each point
  ! the absolute value of the charge
  ! local exchange energy
  ! local correlation energy
  ! local exchange potential
  ! local correlation potential
  ! argument of the phase factor

  complex (DP) ::  exc, work, work0, work1, work2, work3
  complex (DP), allocatable :: drc_exp (:,:), aux (:), d3dyn0 (:,:,:), &
       d3dyn1 (:,:,:), d3dyn2 (:,:,:), d3dyn3 (:,:,:), d3dyn4 (:,:,:)

  if (.not.nlcc_any) return
  allocate  (aux ( dfftp%nnr))
  allocate  (drc_exp ( ngm, nat))
  allocate  (d3dyn0 ( 3 * nat, 3 * nat, 3 * nat))
  allocate  (d3dyn1 ( 3 * nat, 3 * nat, 3 * nat))
  allocate  (d3dyn2 ( 3 * nat, 3 * nat, 3 * nat))
  allocate  (d3dyn3 ( 3 * nat, 3 * nat, 3 * nat))
  allocate  (d3dyn4 ( 3 * nat, 3 * nat, 3 * nat))

  d3dyn0(:,:,:) = (0.d0, 0.d0)
  d3dyn1(:,:,:) = (0.d0, 0.d0)
  d3dyn2(:,:,:) = (0.d0, 0.d0)
  d3dyn3(:,:,:) = (0.d0, 0.d0)
  drc_exp(:,:) = (0.d0, 0.d0)

  do na = 1, nat
     nta = ityp (na)
     do ig = 1, ngm
        arg = - tpi * (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) &
             + g (3, ig) * tau (3, na) )
        exc = CMPLX(cos (arg), sin (arg) ,kind=DP)
        drc_exp (ig, na) = d0rc (ig, nta) * exc
     enddo
  enddo
  aux(:) = (0.d0, 0.d0)
  do ir = 1, dfftp%nnr
     rhox = rho%of_r (ir, 1) + rho_core (ir)
     arhox = abs (rhox)
     if (arhox > 1.0d-30) then
        call xc (arhox, ex, ec, vx, vc)
        aux (ir) = CMPLX(e2 * (vx + vc), 0.d0,kind=DP)
     endif
  enddo

  CALL fwfft ('Dense', aux, dfftp)
  do na_i = npert_i, npert_f
     na = (na_i - 1) / 3 + 1
     i_cart = na_i - 3 * (na - 1)
     do j_cart = 1, 3
        na_j = j_cart + 3 * (na - 1)
        do k_cart = 1, 3

           na_k = k_cart + 3 * (na - 1)
           work = (0.d0, 0.d0)
           do ig = 1, ngm
              work = work + (0.d0, 1.d0) * g (i_cart, ig) * g (j_cart, ig) &
                   * g (k_cart, ig) * CONJG(aux (nl (ig) ) ) * drc_exp (ig, na)

           enddo

           d3dyn0 (na_i, na_j, na_k) = work * omega * tpiba2 * tpiba
        enddo
     enddo
  enddo
#ifdef __MPI
  do nu_i = 1, 3 * nat
     call davcio_drho (aux, lrdrho, iud0rho, nu_i, - 1)
  enddo
  do nu_i = 1, npert_i - 1
     call davcio_drho (aux, lrdrho, iud0rho, nu_i, - 1)
  enddo
#endif
  do nu_i = npert_i, npert_f
     call davcio_drho (aux, lrdrho, iud0rho, nu_i, - 1)
     do ir = 1, dfftp%nnr
        aux (ir) = aux (ir) * dmuxc (ir, 1, 1)
     enddo

     CALL fwfft ('Dense', aux, dfftp)

     do na = 1, nat
        do i_cart = 1, 3
           na_i = i_cart + 3 * (na - 1)
           do j_cart = 1, 3
              na_j = j_cart + 3 * (na - 1)
              work = (0.d0, 0.d0)
              do ig = 1, ngm
                 work = work - CONJG(aux (nl (ig) ) ) * g (i_cart, ig) * g ( &
                      j_cart, ig) * drc_exp (ig, na)

              enddo
              d3dyn1 (nu_i, na_i, na_j) = work * tpiba2 * omega
           enddo
        enddo
     enddo

  enddo
#ifdef __MPI
  do nu_i = npert_f + 1, 3 * nat
     call davcio_drho (aux, lrdrho, iud0rho, nu_i, - 1)
  enddo
#endif
  drc_exp(:,:) = (0.d0, 0.d0)
  do na = 1, nat
     nta = ityp (na)
     do ig = 1, ngm
        arg = - tpi * ( (g (1, ig) + xq (1) ) * tau (1, na) + (g (2, ig) &
             + xq (2) ) * tau (2, na) + (g (3, ig) + xq (3) ) * tau (3, na) )
        exc = CMPLX(cos (arg), sin (arg) ,kind=DP)
        drc_exp (ig, na) = drc (ig, nta) * exc
     enddo
  enddo
#ifdef __MPI
  do nu_i = 1, 3 * nat
     call davcio_drho (aux, lrdrho, iudrho, nu_i, - 1)
  enddo
  do nu_i = 1, npert_i - 1
     call davcio_drho (aux, lrdrho, iudrho, nu_i, - 1)
  enddo
#endif
  do nu_i = npert_i, npert_f
     call davcio_drho (aux, lrdrho, iudrho, nu_i, - 1)
     do ir = 1, dfftp%nnr
        aux (ir) = aux (ir) * dmuxc (ir, 1, 1)
     enddo

     CALL fwfft ('Dense', aux, dfftp)

     do na = 1, nat
        do i_cart = 1, 3
           na_i = i_cart + 3 * (na - 1)
           do j_cart = 1, 3

              na_j = j_cart + 3 * (na - 1)
              work = (0.d0, 0.d0)
              do ig = 1, ngm
                 work = work - CONJG(aux (nl (ig) ) ) * drc_exp (ig, na) * &
                      (g (i_cart, ig) + xq (i_cart) ) * (g (j_cart, ig) + xq (j_cart) )

              enddo
              d3dyn2 (na_i, nu_i, na_j) = work * omega * tpiba2
              d3dyn3 (na_i, na_j, nu_i) = CONJG(work) * omega * tpiba2
           enddo
        enddo
     enddo
  enddo

#ifdef __MPI
  do nu_i = npert_f + 1, 3 * nat
     call davcio_drho (aux, lrdrho, iudrho, nu_i, - 1)
  enddo
  call mp_sum ( d3dyn0, intra_pool_comm )
  call mp_sum ( d3dyn1, intra_pool_comm )
  call mp_sum ( d3dyn2, intra_pool_comm )
  call mp_sum ( d3dyn3, intra_pool_comm )
  call mp_sum ( d3dyn0, inter_pool_comm )
  call mp_sum ( d3dyn1, inter_pool_comm )
  call mp_sum ( d3dyn2, inter_pool_comm )
  call mp_sum ( d3dyn3, inter_pool_comm )
#endif
  !
  !   The dynamical matrix was computed in cartesian axis and now we put
  !   it on the basis of the modes
  !
  d3dyn4(:,:,:) = (0.d0, 0.d0)
  do nu_k = npert_i, npert_f
     if (q0mode (nu_k) ) then
        do nu_i = 1, 3 * nat

           do nu_j = 1, 3 * nat
              work0 = (0.d0, 0.d0)
              do nc_kcart = 1, 3 * nat
                 do na_icart = 1, 3 * nat
                    do nb_jcart = 1, 3 * nat
                       work0 = work0 + ug0 (nc_kcart, nu_k) * &
                            CONJG(u (na_icart, nu_i) ) * &
                            d3dyn0 (nc_kcart, na_icart, nb_jcart) * &
                            u (nb_jcart, nu_j)
                    enddo
                 enddo
              enddo

              work1 = (0.d0, 0.d0)
              do na_icart = 1, 3 * nat
                 do nb_jcart = 1, 3 * nat
                    work1 = work1 + CONJG(u (na_icart, nu_i) ) * d3dyn1 (nu_k, &
                         na_icart, nb_jcart) * u (nb_jcart, nu_j)
                 enddo
              enddo

              work2 = (0.d0, 0.d0)
              do nc_kcart = 1, 3 * nat
                 do nb_jcart = 1, 3 * nat
                    work2 = work2 + ug0 (nc_kcart, nu_k) * d3dyn2 (nc_kcart, nu_i, &
                         nb_jcart) * u (nb_jcart, nu_j)
                 enddo
              enddo

              work3 = (0.d0, 0.d0)
              do nc_kcart = 1, 3 * nat
                 do na_icart = 1, 3 * nat
                    work3 = work3 + ug0 (nc_kcart, nu_k) * &
                         CONJG(u (na_icart, nu_i) ) * &
                         d3dyn3 (nc_kcart, na_icart, nu_j)
                 enddo
              enddo
              d3dyn4 (nu_k, nu_i, nu_j) = work0 + work1 + work2 + work3
           enddo
        enddo
     endif
  enddo
#ifdef __MPI
  call mp_sum( d3dyn4, inter_pool_comm )
#endif
  d3dyn (:,:,:) = d3dyn(:,:,:) + d3dyn4(:,:,:)
  d3dyn_aux8(:,:,:) = d3dyn4(:,:,:)

  deallocate (aux)
  deallocate (drc_exp)
  deallocate (d3dyn0)
  deallocate (d3dyn1)
  deallocate (d3dyn2)
  deallocate (d3dyn3)
  deallocate (d3dyn4)

  return
end subroutine d3dyn_cc
