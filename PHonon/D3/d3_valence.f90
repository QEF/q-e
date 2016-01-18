!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine d3_valence
  !-----------------------------------------------------------------------
  !
  USE ions_base,  ONLY : nat
  USE kinds, only : DP
  use pwcom
  use qpoint, ONLY: nksq
  use phcom
  use d3com
  USE mp_global,  ONLY : inter_pool_comm, intra_pool_comm
  USE mp,         ONLY : mp_sum

  implicit none
  integer :: ik, ikk, ikq, nu_i, nu_j, nu_k, ibnd, jbnd, kbnd, nrec

  real (DP) :: de1, de2, de3, wg1, wg2, wg3, wwg1, wwg2, d_dos, wrk, &
       wga (nbnd), wgq (nbnd), w0g (nbnd), w1g (nbnd)
   real (DP), external ::  wgauss, w0gauss, w_1gauss

  complex (DP) :: wrk1, aux (3 * nat)
  complex (DP), allocatable :: pdvp_i (:,:), pdvp_j (:,:),  dpsidvpsi (:,:),  &
       pdvp_k (:,:), aux1 (:,:,:), aux2 (:,:,:), aux3 (:,:,:), aux4 (:,:,:)

  if (degauss == 0.d0) return
  allocate  (pdvp_i( nbnd, nbnd))
  allocate  (pdvp_j( nbnd, nbnd))
  allocate  (pdvp_k( nbnd, nbnd))
  allocate  (aux1  ( 3 * nat, 3 * nat, 3 * nat))
  allocate  (aux2  ( 3 * nat, 3 * nat, 3 * nat))
  allocate  (aux3  ( 3 * nat, 3 * nat, 3 * nat))
  allocate  (aux4  ( 3 * nat, 3 * nat, 3 * nat))
  allocate  (dpsidvpsi( nbnd, nbnd))

  aux1(:,:,:) = (0.d0, 0.d0)
  aux2(:,:,:) = (0.d0, 0.d0)
  aux3(:,:,:) = (0.d0, 0.d0)
  aux4(:,:,:) = (0.d0, 0.d0)

  call read_ef

  do ik = 1, nksq
     if (lgamma) then
        ikk = ik
        ikq = ik
     else
        ikk = 2 * ik - 1
        ikq = 2 * ik

     endif
     do ibnd = 1, nbnd
        wga (ibnd) = wgauss ( (ef - et (ibnd, ikk) ) / degauss, ngauss)
        wgq (ibnd) = wgauss ( (ef - et (ibnd, ikq) ) / degauss, ngauss)
        w0g (ibnd) = w0gauss ( (ef - et (ibnd, ikk) ) / degauss, ngauss) &
             / degauss
        w1g (ibnd) = w_1gauss ( (ef - et (ibnd, ikk) ) / degauss, ngauss) &
             / (degauss**2)
     enddo

     do nu_i = 1, 3 * nat
        if (q0mode (nu_i) ) then
           nrec = nu_i + (ik - 1) * 3 * nat

           call davcio (pdvp_i, lrpdqvp, iupd0vp, nrec, - 1)
           do nu_j = 1, 3 * nat
              nrec = nu_j + (ik - 1) * 3 * nat

              call davcio (pdvp_j, lrpdqvp, iupdqvp, nrec, - 1)
              do nu_k = 1, 3 * nat
                 nrec = nu_k + (ik - 1) * 3 * nat

                 call davcio (pdvp_k, lrpdqvp, iupdqvp, nrec, - 1)
                 do ibnd = 1, nbnd
                    wg1 = wga (ibnd)

                    wwg1 = w0g (ibnd)
                    do jbnd = 1, nbnd
                       wg2 = wga (jbnd)

                       wwg2 = w0g (jbnd)
                       de1 = et (ibnd, ikk) - et (jbnd, ikk)
                       do kbnd = 1, nbnd
                          wg3 = wgq (kbnd)
                          de2 = et (jbnd, ikk) - et (kbnd, ikq)
                          de3 = et (kbnd, ikq) - et (ibnd, ikk)
                          if (abs (de1) < 2.0d-5 .and. abs (de2) < 2.0d-5 &
                               .and. abs (de3) < 2.0d-5) then
                             wrk = 0.5d0 * w1g (ibnd)
                          elseif (abs (de1) < 1.0d-5) then
                             wrk = ( (wg1 - wg3) / de2 + wwg1) / de3
                          elseif (abs (de2) < 1.0d-5) then
                             wrk = ( (wg2 - wg1) / de3 + wwg2) / de1
                          elseif (abs (de3) < 1.0d-5) then
                             wrk = ( (wg3 - wg2) / de1 + wwg1) / de2
                          else
                             wrk = - (wg1 * de2 + wg2 * de3 + wg3 * de1) / &
                                  (de1 * de2 * de3)
                          endif
                          aux1 (nu_i, nu_j, nu_k) = aux1 (nu_i, nu_j, nu_k) + &
                               2.d0 * wrk * wk (ikk) * pdvp_i (ibnd, jbnd) * &
                               CONJG(pdvp_j (kbnd, jbnd) ) * pdvp_k (kbnd, ibnd)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        endif
     enddo
  enddo

  do ik = 1, nksq
     if (lgamma) then
        ikk = ik
        ikq = ik
     else
        ikk = 2 * ik - 1
        ikq = 2 * ik
     endif

     do nu_j = 1, 3 * nat
        nrec = nu_j + (ik - 1) * 3 * nat

        call davcio (pdvp_j, lrpdqvp, iupdqvp, nrec, - 1)
        do nu_k = 1, 3 * nat
           nrec = nu_k + (ik - 1) * 3 * nat

           call davcio (pdvp_k, lrpdqvp, iupdqvp, nrec, - 1)
           nrec = nu_j + (nu_k - 1) * 3 * nat + (ik - 1) * 9 * nat * nat

           call davcio (dpsidvpsi, lrdpdvp, iudpdvp_1, nrec, - 1)
           do nu_i = 1, 3 * nat

              if (q0mode (nu_i) .or.lgamma) then
                 wrk1 = CMPLX(0.d0, 0.d0,kind=DP)
                 do ibnd = 1, nbnd
                    do jbnd = 1, nbnd
                       de1 = et (ibnd, ikk) - et (jbnd, ikq)
                       if (abs (de1) > 1.0d-5) then
                          wrk = (w0gauss ( (ef - et (ibnd, ikk) ) / degauss, ngauss) &
                               / degauss - w0gauss ( (ef - et (jbnd, ikq) ) / degauss, &
                               ngauss) / degauss) / de1
                       else
                          wrk = - w_1gauss ( (ef - et (ibnd, ikk) ) / degauss, ngauss) &
                               / (degauss**2)
                       endif
                       wrk1 = wrk1 + wk (ikk) * wrk * ef_sh (nu_i) * CONJG(pdvp_j ( &
                            jbnd, ibnd) ) * pdvp_k (jbnd, ibnd)
                    enddo
                 enddo

                 aux2 (nu_i, nu_j, nu_k) = aux2 (nu_i, nu_j, nu_k) + wrk1
                 if (lgamma) then
                    aux2 (nu_k, nu_i, nu_j) = aux2 (nu_k, nu_i, nu_j) + wrk1
                    aux2 (nu_j, nu_k, nu_i) = aux2 (nu_j, nu_k, nu_i) + wrk1
                 endif
                 wrk1 = CMPLX(0.d0, 0.d0,kind=DP)
                 do ibnd = 1, nbnd
                    wrk1 = wrk1 + wk (ikk) * ef_sh (nu_i) * dpsidvpsi (ibnd, ibnd) &
                         * w0gauss ( (ef - et (ibnd, ikk) ) / degauss, ngauss) / &
                         degauss
                 enddo
                 aux2 (nu_i, nu_j, nu_k) = aux2 (nu_i, nu_j, nu_k) + wrk1
                 aux2 (nu_i, nu_k, nu_j) = aux2 (nu_i, nu_k, nu_j) + &
                      CONJG(wrk1)
                 if (lgamma) then
                    aux2 (nu_k, nu_i, nu_j) = aux2 (nu_k, nu_i, nu_j) + wrk1
                    aux2 (nu_j, nu_i, nu_k) = aux2 (nu_j, nu_i, nu_k) + &
                         CONJG(wrk1)
                    aux2 (nu_j, nu_k, nu_i) = aux2 (nu_j, nu_k, nu_i) + wrk1
                    aux2 (nu_k, nu_j, nu_i) = aux2 (nu_k, nu_j, nu_i) + &
                         CONJG(wrk1)
                 endif
              endif
           enddo
        enddo
     enddo
  enddo
  if (lgamma) then
     do nu_i = 1, 3 * nat
        if (.not.q0mode (nu_i) ) then
           do nu_j = 1, 3 * nat
              do nu_k = 1, 3 * nat
                 aux2 (nu_i, nu_j, nu_k) = CMPLX(0.d0, 0.d0,kind=DP)
              enddo
           enddo
        endif
     enddo
  endif
  if (lgamma) then
     d_dos = 0.d0
     aux(:) = (0.d0, 0.d0)
     do ik = 1, nksq
        ikk = ik
        do ibnd = 1, nbnd
           d_dos = d_dos + wk (ikk) * w_1gauss ( (ef - et (ibnd, ikk) ) &
                / degauss, ngauss) / (degauss**2)
        enddo
        do nu_i = 1, 3 * nat
           nrec = nu_i + (ik - 1) * 3 * nat
           call davcio (pdvp_i, lrpdqvp, iupd0vp, nrec, - 1)
           do ibnd = 1, nbnd
              aux (nu_i) = aux (nu_i) + pdvp_i (ibnd, ibnd) * wk (ikk) &
                   * w_1gauss ( (ef - et (ibnd, ikk) ) / degauss, ngauss) / &
                   (degauss**2)
           enddo
        enddo
     enddo
     do nu_i = 1, 3 * nat
        if (q0mode (nu_i) ) then
           do nu_j = 1, 3 * nat
              do nu_k = 1, 3 * nat
                 aux3 (nu_i, nu_j, nu_k) = aux3 (nu_i, nu_j, nu_k) + &
                      ef_sh (nu_i) * ef_sh (nu_j) * aux (nu_k) + &
                      ef_sh (nu_j) * ef_sh (nu_k) * aux (nu_i) + &
                      ef_sh (nu_k) * ef_sh (nu_i) * aux (nu_j)
                 aux4 (nu_i, nu_j, nu_k) = aux4 (nu_i, nu_j, nu_k) - &
                      ef_sh (nu_i) * ef_sh (nu_j) * ef_sh (nu_k) * d_dos
              enddo
           enddo
        endif
     enddo
  endif
#ifdef __MPI
  call mp_sum( aux1, inter_pool_comm )
  call mp_sum( aux2, inter_pool_comm )
  if (lgamma) then
     call mp_sum( aux3, inter_pool_comm )
     call mp_sum( aux4, inter_pool_comm )
  endif
#endif
  d3dyn = d3dyn + aux1 + aux2 + aux3 + aux4
  d3dyn_aux7 = d3dyn_aux7 + aux1 + aux2 + aux3 + aux4

  deallocate (pdvp_i)
  deallocate (pdvp_j)
  deallocate (pdvp_k)
  deallocate (aux1)
  deallocate (aux2)
  deallocate (aux3)
  deallocate (aux4)
  deallocate (dpsidvpsi)

  return
end subroutine d3_valence
