!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine dvpsi_e2
  !-----------------------------------------------------------------------
  !
  ! This routine shold be called before the self-consistent cycle used to
  ! compute the second derivative of the wavefunctions with respect to
  ! electric-fields. It computes that part of the potential that remains
  ! constant during the cycle.
  !
  USE kinds,           ONLY : DP
  USE cell_base,       ONLY : omega
  USE klist,           ONLY : wk, ngk
  USE gvecs,           ONLY : doublegrid
  USE wvfct,           ONLY : npwx, nbnd
  USE wavefunctions_module, ONLY: evc
  USE buffers,         ONLY : get_buffer
  USE fft_base,        ONLY : dfftp, dffts
  USE scf,             ONLY : rho
  USE qpoint,          ONLY : nksq
  USE units_ph,        ONLY : lrdrho, iudrho, lrdwf, iudwf, lrwfc, iuwfc
  USE control_lr,      ONLY : nbnd_occ
  USE ramanm,          ONLY : lrba2, iuba2, lrchf, iuchf, a1j, a2j
  USE mp_pools,        ONLY : my_pool_id, inter_pool_comm
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp,        ONLY: mp_sum
  USE dv_of_drho_lr

  implicit none

  INTEGER :: npw, npwq
  integer :: ik, ipa, ipb, ir, ibnd, jbnd, nrec
  ! counter on k-points
  ! counter on polarizations
  ! counter on points of the real-space mesh
  ! counter on bands
  ! the record number
  real(DP), allocatable  :: raux6 (:,:), d2muxc (:)
  ! function on the real space smooth-mesh
  ! second derivative of the XC-potential
  real(DP) ::  d2mxc, rhotot
  ! external function
  ! total charge on a point
  complex(DP), allocatable :: depsi (:,:,:), auxg (:,:), auxs1 (:), &
               auxs2 (:), aux3s (:,:), aux3 (:,:), ps (:,:,:,:)
  ! d |psi> / dE  (E=electric field)
  ! chi-wavefunction
  ! function on the real space smooth-mesh
  ! function on the real space smooth-mesh
  ! function on the real space smooth-mesh
  ! function on the real space thick-mesh
  complex(DP), pointer :: aux6s (:,:), aux6 (:,:)
  ! function on the real space smooth-mesh
  ! function on the real space thick-mesh
  complex(DP) :: tmp, weight
  ! working space
  ! weight in k-point summation
  !
  call start_clock('dvpsi_e2')
  !
  ! First, calculates the second derivative of the charge-density
  ! -only the part that does not depend on the self-consistent cycle-
  !
  allocate (raux6  (dffts%nnr,6))
  allocate (depsi  (npwx,nbnd,3))
  allocate (aux3s  (dffts%nnr,3))
  allocate (ps     (nbnd,nbnd,3,3))

  raux6 (:,:) = 0.d0
  do ik = 1, nksq
     npw = ngk(ik)
     npwq = npw
     if (nksq.gt.1) call get_buffer (evc, lrwfc, iuwfc, ik)
     weight = 2.d0 * wk(ik) / omega

     do ipa = 1, 3
        nrec = (ipa - 1) * nksq + ik
        call get_buffer (depsi (1, 1, ipa), lrdwf, iudwf, nrec)
     enddo

     do ibnd = 1, nbnd_occ (ik)
        do ipa = 1, 3
           call cft_wave (ik, depsi (1, ibnd, ipa), aux3s (1, ipa), +1)
        enddo
        do ipa = 1, 6
           do ir = 1, dffts%nnr
              tmp = CONJG(aux3s (ir, a1j (ipa))) * aux3s (ir, a2j (ipa))
              raux6 (ir, ipa) = raux6 (ir, ipa) + weight *  DBLE (tmp)
           enddo
        enddo
     enddo

     do ipa = 1, 3
        do ipb = 1, 3
           CALL zgemm( 'C', 'N', nbnd_occ (ik), nbnd_occ (ik), npwq, &
                (1.d0,0.d0), depsi(1,1, ipa), npwx, depsi(1,1,ipb), npwx, &
                (0.d0,0.d0), ps(1,1,ipa,ipb), nbnd )
        enddo
     enddo
     call mp_sum ( ps, intra_bgrp_comm )

     do ibnd = 1, nbnd_occ (ik)
        call cft_wave (ik, evc (1, ibnd), aux3s (1,1), +1)
        do jbnd = 1, nbnd_occ (ik)
           call cft_wave (ik, evc (1, jbnd), aux3s (1,2), +1)
           do ipa = 1, 6
              do ir = 1, dffts%nnr
                 tmp =  aux3s (ir,1) *                           &
                        ps(ibnd, jbnd, a1j (ipa), a2j (ipa)) *   &
                        CONJG(aux3s (ir,2))
                 raux6 (ir, ipa) = raux6 (ir, ipa) - weight *  DBLE (tmp)
              enddo
           enddo
        enddo
     enddo

  enddo

  deallocate (depsi)
  deallocate (aux3s)
  deallocate (ps)

  !
  ! Multiplies the charge with the potential
  !
  if (doublegrid) then
     allocate (auxs1  (dffts%nnr))
     allocate (aux6   (dfftp%nnr,6))
  else
     allocate (aux6s  (dffts%nnr,6))
     aux6 => aux6s
  endif

  do ipa = 1, 6
     if (doublegrid) then
        do ir = 1, dffts%nnr
           auxs1 (ir) = CMPLX(raux6 (ir, ipa), 0.d0,kind=DP)
        enddo
        call cinterpolate (aux6 (1, ipa), auxs1, +1)
     else
        do ir = 1, dffts%nnr
           aux6 (ir, ipa) = CMPLX(raux6 (ir, ipa), 0.d0,kind=DP)
        enddo
     endif
     call dv_of_drho (aux6(:, ipa), .false.)
  enddo

  if (doublegrid) deallocate (auxs1)
  deallocate (raux6)

  !
  ! Calculates the term depending on the third derivative of the
  !                     Exchange-correlation energy
  !
  allocate (d2muxc (dfftp%nnr))
  allocate (aux3   (dfftp%nnr,3))
  do ipa = 1, 3
     call davcio_drho (aux3 (1, ipa), lrdrho, iudrho, ipa, -1)
  enddo

#if defined(__MPI)
  if (my_pool_id .ne. 0) goto 100
#endif
  d2muxc (:) = 0.d0
  do ir = 1, dfftp%nnr
!     rhotot = rho%of_r(ir,1) + rho_core(ir)
     rhotot = rho%of_r(ir,1)
     if ( rhotot.gt. 1.d-30 ) d2muxc(ir)= d2mxc( rhotot)
     if ( rhotot.lt.-1.d-30 ) d2muxc(ir)=-d2mxc(-rhotot)
  enddo

  do ipa = 1, 6
     do ir = 1, dfftp%nnr
        aux6 (ir, ipa) = aux6 (ir, ipa) + d2muxc (ir) * &
                   aux3 (ir, a1j (ipa)) * aux3 (ir, a2j (ipa))
     enddo
  enddo

 100  continue
  call mp_sum ( aux6, inter_pool_comm )
  call psyme2 (aux6)

  deallocate (d2muxc)
  deallocate (aux3)


  if (doublegrid) then
     allocate (aux6s  (dffts%nnr,6))
     do ipa = 1, 6
        call cinterpolate (aux6 (1, ipa), aux6s (1, ipa), -1)
     enddo
     deallocate (aux6)
  endif
  !
  ! Multiplies the obtained potential with the wavefunctions and
  ! writes the results on iuba2; a faster way of proceeding would
  ! be that of keeping the potential in memory and use it directly in
  ! solve_e2
  !
  allocate (auxg   (npwx,nbnd))
  allocate (auxs1  (dffts%nnr))
  allocate (auxs2  (dffts%nnr))

  do ik = 1, nksq
     npw = ngk(ik)
     npwq = npw
     if (nksq.gt.1) call get_buffer (evc, lrwfc, iuwfc, ik)
     do ipa = 1, 6
        nrec = (ipa - 1) * nksq + ik
        call davcio (auxg, lrchf, iuchf, nrec, -1)
        do ibnd = 1, nbnd_occ (ik)
           call cft_wave (ik, evc (1, ibnd), auxs1, +1)
           do ir = 1, dffts%nnr
              auxs2 (ir) = auxs1 (ir) * aux6s (ir, ipa)
           enddo
           call cft_wave (ik, auxg (1, ibnd), auxs2, -1)
        enddo
        nrec = (ipa - 1) * nksq + ik
        call davcio (auxg, lrba2, iuba2, nrec, +1)
     enddo
  enddo

  deallocate (auxg)
  deallocate (auxs1)
  deallocate (auxs2)

  deallocate (aux6s)
  call stop_clock('dvpsi_e2')

  return
end subroutine dvpsi_e2
