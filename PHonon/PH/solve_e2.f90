!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine solve_e2
  !-----------------------------------------------------------------------
  !! Self consistent cycle to compute the second order derivatives
  !! of the wavefunctions with respect to electric fields.
  !
  USe kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout
  USE cell_base,             ONLY : tpiba2
  USE klist,                 ONLY : ltetra, lgauss, wk, xk, ngk, igk_k
  USE lsda_mod,              ONLY : lsda, nspin
  USE gvect,                 ONLY : g
  USE gvecs,                 ONLY : doublegrid
  USE fft_base,              ONLY : dfftp, dffts
  USE fft_interfaces,        ONLY : fft_interpolate
  USE wvfct,                 ONLY : npwx, nbnd, et
  USE noncollin_module,      ONLY : nspin_mag
  USE buffers,   ONLY: get_buffer
  USE ions_base, ONLY: nat
  USE uspp,      ONLY: okvan, nkb, vkb
  USE uspp_param,ONLY : nhm
  USE wavefunctions,  ONLY: evc
  USE units_lr,   ONLY : lrwfc, iuwfc
  USE ramanm,     ONLY : lrba2, iuba2, lrd2w, iud2w
  USE recover_mod, ONLY : read_rec, write_rec

  USE check_stop, ONLY: check_stop_now
  USE mp_pools,   ONLY : inter_pool_comm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum

  USE eqv,       ONLY : dpsi, dvpsi
  USE qpoint,    ONLY : nksq, ikks, ikqs
  USE control_lr, ONLY : nbnd_occ, lgamma, convt, nmix_ph, alpha_mix, tr2_ph, niter_ph, &
                         rec_code, rec_code_read, flmixdpot
  USE dv_of_drho_lr
  USE uspp_init,        ONLY : init_us_2
  USE dfpt_type,        ONLY : dfpt_data_type, allocate_dfpt_data, deallocate_dfpt_data, &
                               dfpt_dvscfp_to_dvscfs
  USE incdrhoscf_mod,   ONLY : incdrhoscf

  implicit none

  real(DP) ::  thresh, weight, avg_iter, dr2
  ! convergence threshold for the solution of the
  ! linear system
  ! used for summation over k points
  ! average number of iterations
  ! convergence limit
  complex(DP) , allocatable :: aux1 (:)
  ! auxiliary space
  COMPLEX(DP), ALLOCATABLE :: dvscftmp(:, :, :)
   !! change of scf potential (output before mixing)
  TYPE(dfpt_data_type) :: dfpt_data
  !! Data that describes linear response quantities

  logical :: exst
  ! used to open the recover file

  integer :: npw, npwq, ikk, ikq
  integer :: kter, iter0, ipol, ibnd, iter, ik, is, ig, iig, irr, ir, nrec, ios
  ! counter on iterations
  ! counter on perturbations
  ! counter on bands
  ! counter on iterations
  ! counter on k points
  ! counter on G vectors
  ! counter on g vectors
  ! counter on mesh points
  ! the record number
  ! integer variable for I/O control

  if (lsda) call errore ('solve_e2', ' LSDA not implemented', 1)
  if (okvan) call errore ('solve_e2', ' Ultrasoft PP not implemented', 1)

  call start_clock('solve_e2')
  CALL allocate_dfpt_data(dfpt_data, 6)
  ALLOCATE(dvscftmp(dfftp%nnr, nspin_mag, 6))
  allocate (aux1(dffts%nnr))
  convt=.FALSE.
  if (rec_code_read == -10) then
     ! restarting in Raman
     CALL read_rec(dr2, iter0, dfpt_data)
  else
     iter0 = 0
  end if
  if (convt) goto 155
  !
  if ( (lgauss .or. ltetra) .or..not.lgamma) &
        call errore ('solve_e2', 'called in the wrong case', 1)
  !
  !   The outside loop is over the iterations
  !

  do kter = 1, niter_ph

     iter = kter + iter0
     avg_iter = 0.d0

     dfpt_data%drhop = (0.d0, 0.d0)
     dfpt_data%drhos = (0.d0, 0.d0)
     dfpt_data%dbecsum = (0.d0, 0.d0)

     do ik = 1, nksq
        ! in this routine, ikk=ikq=ik always (q=0)
        ikk = ikks(ik)
        ikq = ikqs(ik)
        !
        ! reads unperturbed wavefunctions psi_k in G_space, for all bands
        !
        if (nksq.gt.1) then 
           call get_buffer(evc, lrwfc, iuwfc, ikk)
           !$acc update device(evc)
        endif
        npw = ngk(ikk)
        npwq= ngk(ikq)
        !
        ! compute beta functions and kinetic energy for k-point ikk
        ! needed by h_psi, called by pcgreen
        !
        CALL init_us_2 (npw, igk_k(1,ikk), xk (1, ikk), vkb)
        CALL g2_kin(ikk)
        !
        ! The counter on the polarizations runs only on the 6 inequivalent
        ! indexes --see the comment on raman.F--
        !
        do ipol = 1, 6
           nrec = (ipol - 1) * nksq + ik

           if (kter.eq.1) then
              dpsi (:,:) = (0.d0, 0.d0)
           else
              call davcio (dpsi, lrd2w, iud2w, nrec, -1)
           endif

           if (iter.eq.1) then
              dfpt_data%dvscfs (:,:,:) = (0.d0, 0.d0)
              call davcio (dvpsi, lrba2, iuba2, nrec, -1)
              thresh = 1.0d-2
           else
              call davcio (dvpsi, lrba2, iuba2, nrec, -1)
              do ibnd = 1, nbnd_occ (ik)
                 call cft_wave (ik, evc (1, ibnd), aux1, +1)
                 do ir = 1, dffts%nnr
                    aux1 (ir) = aux1 (ir) * dfpt_data%dvscfs(ir, 1, ipol)
                 enddo
                 call cft_wave (ik, dvpsi (1, ibnd), aux1, -1)
              enddo
              thresh = min (0.1d0 * sqrt(dr2), 1.0d-2)
           endif

           call pcgreen (avg_iter, thresh, ik, et (1, ik) )
           call davcio ( dpsi, lrd2w, iud2w, nrec, +1)
           !
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
           weight = wk (ik)
           call incdrhoscf (dfpt_data%drhos(:,1,ipol), weight, ik, &
                            dfpt_data%dbecsum (1,1,1,ipol), dpsi)
           enddo   ! on perturbations
        enddo      ! on k points
     call mp_sum ( dfpt_data%dbecsum, intra_bgrp_comm )
     if (doublegrid) then
        do is = 1, nspin
           do ipol = 1, 6
              call fft_interpolate (dffts, dfpt_data%drhos(:, is, ipol), dfftp, dfpt_data%drhop(:, is, ipol))
           enddo
        enddo
     else
        CALL zcopy(dffts%nnr*nspin_mag*6, dfpt_data%drhos, 1, dfpt_data%drhop, 1)
     endif
     !
     !   After the loop over the perturbations we have the change of the pote
     !   for all the modes, and we symmetrize this potential
     !
     call mp_sum ( dfpt_data%drhop, inter_pool_comm )
     do ipol = 1, 6
        CALL zcopy(dfftp%nnr*nspin_mag, dfpt_data%drhop(1, 1, ipol), 1, dvscftmp(1, 1, ipol), 1)
        call dv_of_drho(dvscftmp(1, 1, ipol))
     enddo

     call psyme2(dvscftmp)
     !
     ! Mixing with the old potential
     !
     call mix_potential (2 * 6 * dfftp%nnr* nspin, dvscftmp, dfpt_data%dvscfp,  &
                         alpha_mix (kter), dr2, 6 * tr2_ph, iter,  &
                         nmix_ph, flmixdpot, convt)
     !
     CALL dfpt_dvscfp_to_dvscfs(dfpt_data)

     write (6, "(//,5x,' iter # ',i3, &
          &      '   av.it.: ',f5.1)") iter, avg_iter / (6.d0 * nksq)
     dr2 = dr2 / 6
     write (6, "(5x,' thresh=',es10.3, ' alpha_mix = ',f6.3, &
          &      ' |ddv_scf|^2 = ',es10.3 )") thresh, alpha_mix (kter), dr2
     !
     FLUSH( stdout )
     !
     ! rec_code: state of the calculation
     ! rec_code=-10  to -19 Raman
     rec_code=-10
     CALL write_rec('solve_e2..', irr, dr2, iter, convt, dfpt_data)

     if ( check_stop_now() ) call stop_smoothly_ph (.false.)
     if ( convt ) goto 155

  enddo
155 continue
  CALL deallocate_dfpt_data(dfpt_data)
  deallocate (aux1 )

  call stop_clock('solve_e2')

  return
end subroutine solve_e2
