!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine solve_e
  !-----------------------------------------------------------------------
  !! This routine is a driver for the solution of the linear system which
  !! defines the change of the wavefunction due to an electric field.
  !! It performs the following tasks:  
  !! a) computes the bare potential term times \(|\psi\rangle \);  
  !! b) adds to it the screening term \(\Delta V_\text{SCF}|psi\rangle\).
  !!    If \(\text{lda_plus_u}=\text{TRUE}\) compute also the SCF part
  !!    of the response Hubbard potential;  
  !! c) applies \(P_c^+\) (orthogonalization to valence states);  
  !! d) calls \(\texttt{cgsolve_all}\) to solve the linear system;  
  !! e) computes \(\Delta \rho\), \(\Delta V_\text{SCF}|psi\rangle\) and
  !!    symmetrizes them;  
  !! f) if \(\text{lda_plus_u}=\text{TRUE}\) compute also the response
  !!    occupation matrices dnsscf.  
  !! Step b, c, d are done in \(\text{sternheimer_kernel}\).
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat
  USE io_global,             ONLY : stdout, ionode
  USE io_files,              ONLY : diropn
  USE mp,                    ONLY : mp_sum
  USE mp_pools,              ONLY : inter_pool_comm
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE klist,                 ONLY : ltetra, lgauss, xk, ngk, igk_k
  USE gvecs,                 ONLY : doublegrid
  USE fft_base,              ONLY : dfftp, dffts
  USE lsda_mod,              ONLY : nspin, lsda, current_spin, isk
  USE wvfct,                 ONLY : nbnd, npwx
  USE check_stop,            ONLY : check_stop_now
  USE buffers,               ONLY : get_buffer
  USE wavefunctions,         ONLY : evc
  USE uspp,                  ONLY : okvan, vkb
  USE uspp_param,            ONLY : nhm
  USE noncollin_module,      ONLY : noncolin, npol, nspin_mag, domag
  USE scf,                   ONLY : rho
  USE paw_variables,         ONLY : okpaw
  USE paw_onecenter,         ONLY : paw_dpotential
  USE paw_symmetry,          ONLY : paw_desymmetrize

  USE units_ph,              ONLY : lrdrho, iudrho, lrebar, iuebar
  USE units_lr,              ONLY : iuwfc, lrwfc
  USE output,                ONLY : fildrho
  USE control_ph,            ONLY : ext_recover, rec_code, &
                                    lnoloc, convt, tr2_ph, nmix_ph, &
                                    alpha_mix, lgamma_gamma, niter_ph, &
                                    flmixdpot, rec_code_read
  USE recover_mod,           ONLY : read_rec, write_rec
  USE lrus,                  ONLY : int3_paw
  USE qpoint,                ONLY : nksq, ikks
  USE control_lr,            ONLY : lgamma
  USE dv_of_drho_lr,         ONLY : dv_of_drho
  USE fft_interfaces,        ONLY : fft_interpolate
  USE ldaU,                  ONLY : lda_plus_u
  USE apply_dpot_mod,        ONLY : apply_dpot_allocate, apply_dpot_deallocate
  USE response_kernels,      ONLY : sternheimer_kernel
  USE uspp_init,             ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst
  !!
  LOGICAL :: all_conv
  !! True if sternheimer_kernel is converged at all k points and perturbations
  INTEGER :: ikk, npw, kter, iter0, ipol, iter, ik, is, ndim
  !! counters
  REAL(DP) :: thresh
  !! convergence threshold
  REAL(DP) :: averlt
  !! average number of iterations
  REAL(DP) :: dr2
  !! self-consistency error
  REAL(DP) :: tcpu, get_clock
  !! timing variables

  COMPLEX(DP), ALLOCATABLE, TARGET :: dvscfin (:,:,:)
  !! change of the scf potential (input)
  COMPLEX(DP), POINTER :: dvscfins (:,:,:)
  !! change of the scf potential (smooth)
  COMPLEX(DP), ALLOCATABLE :: dvscfout (:,:,:)
  !! change of the scf potential (output)
  COMPLEX(DP), ALLOCATABLE :: dbecsum(:,:,:,:)
  !! the becsum with dpsi
  COMPLEX(DP), ALLOCATABLE :: dbecsum_nc(:,:,:,:,:)
  !! the becsum with dpsi
  COMPLEX(DP), ALLOCATABLE :: mixin(:), mixout(:)
  !! auxiliary for paw mixing
  INTEGER :: nnr
  !
  call start_clock ('solve_e')
  !
  !  This routine is task group aware
  !
  allocate (dvscfin( dfftp%nnr, nspin_mag, 3))
  nnr = dfftp%nnr
  dvscfin=(0.0_DP,0.0_DP)
  if (doublegrid) then
     allocate (dvscfins(dffts%nnr, nspin_mag, 3))
     nnr = dffts%nnr
  else
     dvscfins => dvscfin
  endif
  !$acc enter data create(dvscfins(1:nnr, 1:nspin_mag, 1:3))
  allocate (dvscfout(dfftp%nnr, nspin_mag, 3))
  IF (okpaw) THEN
     ALLOCATE (mixin(dfftp%nnr*nspin_mag*3+(nhm*(nhm+1)*nat*nspin_mag*3)/2) )
     ALLOCATE (mixout(dfftp%nnr*nspin_mag*3+(nhm*(nhm+1)*nat*nspin_mag*3)/2) )
     mixin=(0.0_DP,0.0_DP)
  ENDIF
  allocate (dbecsum( nhm*(nhm+1)/2, nat, nspin_mag, 3))
  IF (noncolin) allocate (dbecsum_nc (nhm, nhm, nat, nspin, 3))
  CALL apply_dpot_allocate()

  if (rec_code_read == -20.AND.ext_recover) then
     ! restarting in Electric field calculation
     IF (okpaw) THEN
        CALL read_rec(dr2, iter0, 3, dvscfin, dvscfins, dvscfout, dbecsum)
        CALL setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
                    mixin, dvscfin, dbecsum, ndim, -1 )
     ELSE
        CALL read_rec(dr2, iter0, 3, dvscfin, dvscfins)
     ENDIF
  else if (rec_code_read > -20 .AND. rec_code_read <= -10) then
     ! restarting in Raman: proceed
     convt = .true.
  else
     convt = .false.
     iter0 = 0
  endif
  !
  IF ( ionode .AND. fildrho /= ' ') THEN
     INQUIRE (UNIT = iudrho, OPENED = exst)
     IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
     CALL diropn (iudrho, TRIM(fildrho)//'.E', lrdrho, exst)
  end if
  IF (rec_code_read > -20) convt=.TRUE.
  !
  if (convt) go to 155
  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !
  if ( (lgauss .or. ltetra) .or..not.lgamma) call errore ('solve_e', &
       'called in the wrong case', 1)
  !
  ! Compute P_c^+ x psi for all polarization and k points and store in buffer
  !
  DO ik = 1, nksq
     DO ipol = 1, 3
        ikk = ikks(ik)
        npw = ngk(ikk)
        IF (lsda) current_spin = isk(ikk)
        !
        ! reads unperturbed wavefunctions psi_k in G_space, for all bands
        !
        IF (nksq > 1) THEN
           CALL get_buffer(evc, lrwfc, iuwfc, ikk)
        ENDIF
        !
        CALL init_us_2(npw, igk_k(1, ikk), xk(1, ikk), vkb, .true.)
        !$acc update host(vkb)
        !
        ! computes P_c^+ x psi_kpoint, written to buffer iuebar.
        !
        CALL dvpsi_e(ik, ipol)
        !
     ENDDO ! ipol
  ENDDO ! ik
  !
  !   The outside loop is over the iterations
  !
  do kter = 1, niter_ph
     !
     FLUSH( stdout )
     iter = kter + iter0
     !
     dvscfout = (0.d0,0.d0)
     dbecsum = (0.d0,0.d0)
     IF (noncolin) dbecsum_nc = (0.d0,0.d0)
     !
     ! DFPT+U: at each iteration calculate dnsscf,
     ! i.e. the scf variation of the occupation matrix ns.
     !
     IF (lda_plus_u .AND. (iter /= 1)) CALL dnsq_scf(3, .false., 0, 1, .false.)
     !
     ! set threshold for the iterative solution of the linear system
     !
     IF (iter == 1) THEN
        thresh = 1.d-2
        IF (lnoloc) thresh = 1.d-5
     ELSE
        thresh = MIN(0.1d0 * SQRT(dr2), 1.0d-2)
     ENDIF
     !
     ! Compute dvscfout, the charge density response to the total potential
     !
     CALL sternheimer_kernel(iter==1, .FALSE., 3, lrebar, iuebar, thresh, dvscfins, &
                             all_conv, averlt, dvscfout, dbecsum, dbecsum_nc)
     !
     !  The calculation of dbecsum is distributed across processors
     !  (see addusdbec) - we sum over processors the contributions
     !  coming from each slice of bands
     !
     IF (noncolin) THEN
        call mp_sum ( dbecsum_nc, intra_bgrp_comm )
     ELSE
        call mp_sum ( dbecsum, intra_bgrp_comm )
     END IF

     if (doublegrid) then
        do is=1,nspin_mag
           do ipol=1,3
              call fft_interpolate (dffts, dvscfout(:,is,ipol), dfftp, dvscfout(:,is,ipol))
           enddo
        enddo
     endif
     !
     IF (noncolin.and.okvan) CALL set_dbecsum_nc(dbecsum_nc, dbecsum, 3)
     !
     call addusddense (dvscfout, dbecsum)
     !
     !   dvscfout contains the (unsymmetrized) linear charge response
     !   for the three polarizations - symmetrize it
     !
     call mp_sum ( dvscfout, inter_pool_comm )
     IF (okpaw) call mp_sum ( dbecsum, inter_pool_comm )
     if (.not.lgamma_gamma) then
        call psyme (dvscfout)
        IF ( noncolin.and.domag ) CALL psym_dmage(dvscfout)
     endif
     !
     !   save the symmetrized linear charge response to file
     !   calculate the corresponding linear potential response
     !
     do ipol=1,3
        if (fildrho.ne.' ') call davcio_drho(dvscfout(1,1,ipol),lrdrho, &
             iudrho,ipol,+1)
        IF (lnoloc) then
           dvscfout(:,:,ipol)=(0.d0,0.d0)
        ELSE
           call dv_of_drho (dvscfout (1, 1, ipol), .false.)
        ENDIF
     enddo
     !
     !   mix the new potential with the old
     !
     IF (okpaw) THEN
     !
     !  In this case we mix also dbecsum
     !
        call setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
                    mixout, dvscfout, dbecsum, ndim, -1 )
        call mix_potential (2*3*dfftp%nnr*nspin_mag+2*ndim, mixout, mixin, &
                         alpha_mix(kter), dr2, 3*tr2_ph/npol, iter, &
                         nmix_ph, flmixdpot, convt)
        call setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
                       mixin, dvscfin, dbecsum, ndim, 1 )
     ELSE
        call mix_potential (2*3*dfftp%nnr*nspin_mag, dvscfout, dvscfin, alpha_mix ( &
          kter), dr2, 3 * tr2_ph / npol, iter, nmix_ph, flmixdpot, convt)
     ENDIF
     if (doublegrid) then
        do is=1,nspin_mag
           do ipol = 1, 3
              call fft_interpolate (dfftp, dvscfin(:,is,ipol), dffts, dvscfins(:,is,ipol))
           enddo
        enddo
     endif

     IF (okpaw) THEN
        IF (noncolin) THEN
!           call PAW_dpotential(dbecsum_nc,becsum_nc,int3_paw,3)
        ELSE
!
!    The presence of c.c. in the formula gives a factor 2.0
!
           dbecsum=2.0_DP * dbecsum
           IF (.NOT. lgamma_gamma) CALL PAW_desymmetrize(dbecsum)
           call PAW_dpotential(dbecsum,rho%bec,int3_paw,3)
        ENDIF
     ENDIF

     call newdq(dvscfin,3)

     tcpu = get_clock ('PHONON')
     WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
          &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
     dr2 = dr2 / 3
     WRITE( stdout, "(5x,' thresh=',es10.3, ' alpha_mix = ',f6.3, &
          &      ' |ddv_scf|^2 = ',es10.3 )") thresh, alpha_mix (kter), dr2
     !
     FLUSH( stdout )
     !
     ! rec_code: state of the calculation
     ! rec_code=-20 Electric Field
     !
     rec_code=-20
     IF (okpaw) THEN
        CALL write_rec('solve_e...', 0, dr2, iter, convt, 3, dvscfin, &
                                                       dvscfout, dbecsum)
     ELSE
        CALL write_rec('solve_e...', 0, dr2, iter, convt, 3, dvscfin)
     ENDIF

     if (check_stop_now()) call stop_smoothly_ph (.false.)

     if (convt) goto 155

  enddo
155 continue
  !
  CALL apply_dpot_deallocate()
  deallocate (dbecsum)
  deallocate (dvscfout)
  IF (okpaw) THEN
     DEALLOCATE(mixin)
     DEALLOCATE(mixout)
  ENDIF
  !$acc exit data delete(dvscfins)
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)
  if (noncolin) deallocate(dbecsum_nc)

  call stop_clock ('solve_e')
  return
end subroutine solve_e
