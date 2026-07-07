!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
!
! Define batch size for perturbations processing
! Can be overridden with -DPERTS_PER_BATCH=N at compile time
#ifndef PERTS_PER_BATCH
#define PERTS_PER_BATCH 8
#endif
!
subroutine drho
  !-----------------------------------------------------------------------
  !! Here we compute, for each mode the change of the charge density
  !! due to the displacement, at fixed wavefunctions. These terms
  !! are saved on disk. The orthogonality part is included in the
  !! computed change.
  !
  !
  !
  USE kinds,      ONLY : DP
  USE gvecs,      ONLY : doublegrid
  USE fft_base,   ONLY : dfftp, dffts
  USE lsda_mod,   ONLY : nspin
  USE cell_base,  ONLY : omega
  USE ions_base,  ONLY : nat
  USE buffers,    ONLY : save_buffer
  USE noncollin_module, ONLY : noncolin, npol, nspin_lsda, nspin_mag
  USE uspp_param, ONLY : upf, nhm
  USE uspp,       ONLY : okvan, nkb
  USE wvfct,      ONLY : nbnd
  USE paw_variables,    ONLY : okpaw
  USE control_ph, ONLY : all_done
  USE lrus,       ONLY : becp1
  USE qpoint,     ONLY : nksq
  USE control_lr, ONLY : lgamma, rec_code_read

  USE dynmat,     ONLY : dyn00
  USE modes,      ONLY : npertx, npert, nirr, u
  USE phus,       ONLY : becsumort, alphap
  USE units_ph,   ONLY : lrdrhous, iudrhous

  USE mp_pools,   ONLY : inter_pool_comm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  USE becmod,     ONLY : bec_type, allocate_bec_type, deallocate_bec_type
  USE fft_interfaces, ONLY : fft_interpolate

  implicit none

  integer :: mode, is, ir, irr, iper, npe, nrstot, nu_i, nu_j, ik, &
             ipol, irr_batch, irr_end, npe_total, npe_offset, mode_batch_start
  ! counter on modes
  ! counter on atoms and polarizations
  ! counter on atoms
  ! counter on spin
  ! counter on perturbations
  ! the number of points
  ! counter on modes
  ! counter on k-point
  ! counter on coordinates

  real(DP), allocatable :: wgg (:,:,:)
  ! the weight of each point


  complex(DP) :: wdyn (3 * nat, 3 * nat)
  type (bec_type), pointer :: becq(:), alpq(:,:)
  complex(DP), allocatable :: dvlocin (:), drhous (:,:,:),&
       drhoust (:,:,:), dbecsum(:,:,:,:), dbecsum_nc(:,:,:,:,:)
  ! auxiliary to store bec at k+q
  ! auxiliary to store alphap at
  ! the change of the local potential
  ! the change of the charge density
  ! the change of the charge density
  ! the derivative

!
!  The PAW case requires dbecsumort so we recalculate this starting part
!  This will be changed soon
!
  if (all_done) return
  if ((rec_code_read >=-20 .and..not.okpaw)) return

  dyn00(:,:) = (0.d0,0.d0)
  if (.not.okvan) return
  call start_clock ('drho')
  !
  !    first compute the terms needed for the change of the charge density
  !    due to the displacement of the augmentation charge
  !
  call compute_becsum_ph()
  !
  call compute_alphasum()
  !
  !    then compute the weights
  !
  call start_clock('compute_wgg')
  allocate (wgg (nbnd ,nbnd , nksq))
  if (lgamma) then
     becq => becp1
     alpq => alphap
  else
     allocate (becq ( nksq))
     allocate (alpq ( 3, nksq))
     do ik =1,nksq
        call allocate_bec_type (  nkb, nbnd, becq(ik))
        DO ipol=1,3
           CALL allocate_bec_type (  nkb, nbnd, alpq(ipol,ik))
        ENDDO
     end do
  endif
  call compute_weight (wgg)
  call stop_clock('compute_wgg')
  !
  !    becq and alpq are sufficient to compute the part of C^3 (See Eq. 37
  !    which does not contain the local potential
  !
  IF (.not.lgamma) call compute_becalp (becq, alpq)
  call start_clock('nldyntot')
  call compute_nldyn (dyn00, wgg, becq, alpq)
  call stop_clock('nldyntot')
  !
  !   now we compute the change of the charge density due to the change of
  !   the orthogonality constraint
  !
  allocate (drhous ( dffts%nnr, nspin_mag , 3 * nat))
  allocate (dbecsum( nhm * (nhm + 1) /2, nat, nspin_mag, 3 * nat))
  dbecsum=(0.d0,0.d0)
  call start_clock('drhous')
  IF (noncolin) THEN
     allocate (dbecsum_nc( nhm, nhm, nat, nspin, 3 * nat))
     dbecsum_nc=(0.d0,0.d0)
     call compute_drhous_nc (drhous, dbecsum_nc, wgg, becq, alpq)
  ELSE
     call compute_drhous (drhous, dbecsum, wgg, becq, alpq)
  ENDIF
  call stop_clock('drhous')

  if (.not.lgamma) then
     do ik=1,nksq
        call deallocate_bec_type(becq(ik))
        DO ipol=1,3
           call deallocate_bec_type(alpq(ipol,ik))
        ENDDO
     end do
     deallocate (becq)
     deallocate (alpq)
  endif
  deallocate (wgg)
  !
  !  The part of C^3 (Eq. 37) which contain the local potential can be
  !  evaluated with an integral of this change of potential and drhous
  !
  allocate (dvlocin(dffts%nnr))

  wdyn (:,:) = (0.d0, 0.d0)
  nrstot = dffts%nr1 * dffts%nr2 * dffts%nr3
  do nu_i = 1, 3 * nat
     call compute_dvloc (u(1, nu_i), .FALSE., dvlocin)
     do nu_j = 1, 3 * nat
        do is = 1, nspin_lsda
        ! FIXME: use zgemm instead of dot_product
           wdyn (nu_j, nu_i) = wdyn (nu_j, nu_i) + &
                dot_product (drhous(1:dffts%nnr,is,nu_j), dvlocin) * &
                omega / DBLE (nrstot)
        enddo
     enddo

  enddo
  !
  ! collect contributions from all pools (sum over k-points)
  !
  call mp_sum ( dyn00, inter_pool_comm )
  call mp_sum ( wdyn, inter_pool_comm )
  !
  ! collect contributions from nodes of a pool (sum over G & R space)
  !
  call mp_sum ( wdyn, intra_bgrp_comm )

  call zaxpy (3 * nat * 3 * nat, (1.d0, 0.d0), wdyn, 1, dyn00, 1)
  !
  !     force this term to be hermitean
  !
  do nu_i = 1, 3 * nat
     do nu_j = 1, nu_i
        dyn00(nu_i,nu_j) = 0.5d0*( dyn00(nu_i,nu_j) + CONJG(dyn00(nu_j,nu_i)))
        dyn00(nu_j,nu_i) = CONJG(dyn00(nu_i,nu_j))
     enddo
  enddo
  !      call tra_write_matrix('drho dyn00',dyn00,u,nat)
  !
  !    add the augmentation term to the charge density and save it
  !
  allocate (drhoust(dfftp%nnr, nspin_mag , PERTS_PER_BATCH*npertx))
  drhoust=(0.d0,0.d0)
  !
  !  The calculation of dbecsum is distributed across processors (see addusdbec)
  !  Sum over processors the contributions coming from each slice of bands
  !
  IF (noncolin) THEN
     call mp_sum ( dbecsum_nc, intra_bgrp_comm )
  ELSE
     call mp_sum ( dbecsum, intra_bgrp_comm )
  END IF

  IF (noncolin.and.okvan) CALL set_dbecsum_nc(dbecsum_nc, dbecsum, 3*nat)

  mode = 0
  if (okpaw) becsumort=(0.0_DP,0.0_DP)
  
  ! Process irr values in batches of PERTS_PER_BATCH
  do irr_batch = 1, nirr, PERTS_PER_BATCH
     irr_end = min(irr_batch + PERTS_PER_BATCH - 1, nirr)
     
     ! Calculate total perturbations in this batch
     npe_total = 0
     do irr = irr_batch, irr_end
        npe_total = npe_total + npert(irr)
     enddo
     
     ! Process all irr values in current batch
     npe_offset = 0
     do irr = irr_batch, irr_end
        npe = npert(irr)
        if (doublegrid) then
           do is = 1, nspin_mag
              do iper = 1, npe
                 call fft_interpolate (dffts, drhous(:,is,mode+iper), dfftp, &
                                     drhoust(:,is,npe_offset+iper))
              enddo
           enddo
        else
           call zcopy (dfftp%nnr*nspin_mag*npe, drhous(1,1,mode+1), 1, &
                      drhoust(1,1,npe_offset+1), 1)
        endif
        npe_offset = npe_offset + npe
        mode = mode + npe
     enddo
     
     ! Scale the density for all perturbations in batch
     call dscal (2*dfftp%nnr*nspin_mag*npe_total, 0.5d0, drhoust, 1)
     
     ! Process augmentation density for entire batch
     mode_batch_start = mode - npe_total
     call addusddens_pulay (drhoust, dbecsum(1,1,1,mode_batch_start+1), mode_batch_start, npe_total)
     
     ! Save buffers for all perturbations in batch
     do iper = 1, npe_total
        nu_i = mode_batch_start + iper
        call save_buffer (drhoust (1, 1, iper), lrdrhous, iudrhous, nu_i)
     enddo
  enddo
   !
   !  Collect the sum over k points in different pools.
   !
   IF (okpaw) call mp_sum ( becsumort, inter_pool_comm )

  deallocate (drhoust)
  deallocate (dvlocin)
  deallocate (dbecsum)
  if (noncolin) deallocate(dbecsum_nc)
  deallocate (drhous)

  call stop_clock ('drho')
  return
end subroutine drho
