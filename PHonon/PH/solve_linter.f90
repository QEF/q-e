!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE solve_linter (irr, imode0, npe, drhoscf)
  !-----------------------------------------------------------------------
  !! Driver routine for the solution of the linear system that
  !! defines the change of the wavefunction due to a lattice distorsion.
  !! It performs the following tasks:
  !! a) computes the bare potential term \(\Delta V | \psi \rangle\)
  !!    and an additional term in the case of US pseudopotentials.
  !!    If \(\text{lda_plus_u}=\text{TRUE}\) compute also the bare
  !!    potential term Delta \(V_\text{hub} | \psi \rangle\);
  !! b) adds to it the screening term \(\Delta V_\text{SCF} | \psi \rangle\).
  !!    If \(\text{lda_plus_u}=\text{TRUE}\) computes also the SCF part
  !!    of the response Hubbard potential;
  !! c) applies \(P_c^+\) (orthogonalization to valence states);
  !! d) calls \(\text{cgsolve_all}\) to solve the linear system;
  !! e) computes \(\Delta\rho\), \(\Delta V_\text{SCF}\) and symmetrizes
  !!    them;
  !! f) If \(\text{lda_plus_u}=\text{TRUE}\) compute also the response
  !!    occupation matrices \(\text{dnsscf}\);
  !! g) --Introduced in February 2020-- If \(\text{noncolin}=\text{TRUE}\)
  !!    and \(\text{domag}=\text{TRUE}\), the linear system is solved twice
  !!    (\(\text{nsolv}=2\), the case \(\text{isolv}=2\) needs the time-reversed
  !!    wave functions). For the theoretical background, please refer to:
  !!    Phys. Rev. B 100, 045115 (2019).
  !! Step b, c, d are done inside sternheimer_kernel.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : prefix, diropn
  USE check_stop,           ONLY : check_stop_now
  USE wavefunctions,        ONLY : evc
  USE cell_base,            ONLY : at
  USE klist,                ONLY : ltetra, lgauss, xk, ngk, igk_k
  USE gvecs,                ONLY : doublegrid
  USE fft_base,             ONLY : dfftp, dffts
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wvfct,                ONLY : nbnd, npwx
  USE scf,                  ONLY : rho, vrs
#if defined(__CUDA)
  USE scf_gpum,             ONLY : vrs_d
#endif
  USE uspp,                 ONLY : okvan, vkb, deeq_nc
  USE uspp_param,           ONLY : nhm
  USE noncollin_module,     ONLY : noncolin, domag, npol, nspin_mag
  USE paw_variables,        ONLY : okpaw
  USE paw_onecenter,        ONLY : paw_dpotential
  USE paw_symmetry,         ONLY : paw_dusymmetrize, paw_dumqsymmetrize
  USE buffers,              ONLY : save_buffer, get_buffer
  USE control_ph,           ONLY : rec_code, niter_ph, nmix_ph, tr2_ph, &
                                   lgamma_gamma, convt, &
                                   alpha_mix, rec_code_read, &
                                   where_rec, flmixdpot, ext_recover
  USE el_phon,              ONLY : elph
  USE uspp,                 ONLY : nlcc_any
  USE units_ph,             ONLY : iudrho, lrdrho, iubar, lrbar, &
                                   iudvscf, iuint3paw, lint3paw
  USE units_lr,             ONLY : iuwfc, lrwfc
  USE output,               ONLY : fildrho, fildvscf
  USE phus,                 ONLY : becsumort, alphap, int1_nc
  USE modes,                ONLY : npertx, u, t, tmq
  USE recover_mod,          ONLY : read_rec, write_rec
  ! used to write fildrho:
  USE dfile_autoname,       ONLY : dfile_name
  USE save_ph,              ONLY : tmp_dir_save
  ! used oly to write the restart file
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm, me_bgrp
  USE mp,                   ONLY : mp_sum
  USE efermi_shift,         ONLY : ef_shift, ef_shift_wfc, def
  USE lrus,                 ONLY : int3_paw, becp1, int3_nc
  USE lr_symm_base,         ONLY : irotmq, minus_q, nsymq, rtau
  USE eqv,                  ONLY : dvpsi
  USE qpoint,               ONLY : xq, nksq, ikks, ikqs
  USE qpoint_aux,           ONLY : ikmks, becpt, alphapt
  USE control_lr,           ONLY : lgamma
  USE dv_of_drho_lr,        ONLY : dv_of_drho
  USE fft_interfaces,       ONLY : fft_interpolate
  USE ldaU,                 ONLY : lda_plus_u
  USE nc_mag_aux,           ONLY : int1_nc_save, deeq_nc_save, int3_save
  USE apply_dpot_mod,       ONLY : apply_dpot_allocate, apply_dpot_deallocate
  USE response_kernels,     ONLY : sternheimer_kernel
  USE uspp_init,            ONLY : init_us_2
  USE sym_def_module,       ONLY : sym_def
  implicit none

  integer :: irr
  !! input: the irreducible representation
  integer :: npe
  !! input: the number of perturbation
  integer :: imode0
  !! input: the position of the modes
  complex(DP) :: drhoscf(dfftp%nnr,nspin_mag,npe)
  !! output: the change of the scf charge
  !
  ! ... local variables
  !
  real(DP) :: thresh, averlt, dr2
  ! thresh: convergence threshold
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  real(DP) :: dos_ef
  ! Misc variables for metals
  ! dos_ef: density of states at Ef

  complex(DP), allocatable, target :: dvscfin(:,:,:)
  ! change of the scf potential
  complex(DP), pointer :: dvscfins (:,:,:)
  ! change of the scf potential (smooth part only)
  complex(DP), allocatable :: drhoscfh (:,:,:), dvscfout (:,:,:)
  ! change of rho / scf potential (output)
  ! change of scf potential (output)
  complex(DP), allocatable :: ldos (:,:), ldoss (:,:), mixin(:), mixout(:), &
       dbecsum (:,:,:,:), dbecsum_nc(:,:,:,:,:,:), aux2(:,:), drhoc(:), &
       dbecsum_aux (:,:,:,:)
  ! Misc work space
  ! ldos : local density of states af Ef
  ! ldoss: as above, without augmentation charges
  ! dbecsum: the derivative of becsum
  ! drhoc: response core charge density
  REAL(DP), allocatable :: becsum1(:,:,:)

  LOGICAL :: all_conv
  !! True if sternheimer_kernel is converged at all k points and perturbations
  logical :: exst,       & ! used to open the recover file
             lmetq0,     & ! true if xq=(0,0,0) in a metal
             first_iter    ! true if first iteration where induced rho is not yet calculated

  integer :: kter,       & ! counter on iterations
             iter0,      & ! starting iteration
             ipert,      & ! counter on perturbations
             iter,       & ! counter on iterations
             ik, ikk,    & ! counter on k points
             ikq,        & ! counter on k+q points
             ndim,       &
             is,         & ! counter on spin polarizations
             nrec,       & ! the record number for dvpsi and dpsi
             mode,       & ! mode index
             isolv,      & ! counter on linear systems
             nsolv,      & ! number of linear systems
             ikmk          ! index of mk

  integer  :: npw, npwq
  integer  :: iq_dummy
  real(DP) :: tcpu, get_clock ! timing variables
  character(len=256) :: filename

  integer :: nnr
  !
  IF (rec_code_read > 20 ) RETURN

  call start_clock ('solve_linter')
!
!  This routine is task group aware
!
  nsolv=1
  IF (noncolin.AND.domag) nsolv=2

  allocate (dvscfin ( dfftp%nnr , nspin_mag , npe))
  nnr = dfftp%nnr
  dvscfin=(0.0_DP,0.0_DP)
  if (doublegrid) then
     allocate (dvscfins (dffts%nnr , nspin_mag , npe))
     nnr = dffts%nnr
  else
     dvscfins => dvscfin
  endif
  !$acc enter data create(dvscfins(1:nnr, 1:nspin_mag, 1:npe))
  allocate (drhoscfh ( dfftp%nnr, nspin_mag , npe))
  allocate (dvscfout ( dfftp%nnr, nspin_mag , npe))
  allocate (dbecsum ( (nhm * (nhm + 1))/2 , nat , nspin_mag , npe))
  IF (okpaw) THEN
     allocate (mixin(dfftp%nnr*nspin_mag*npe+(nhm*(nhm+1)*nat*nspin_mag*npe)/2) )
     allocate (mixout(dfftp%nnr*nspin_mag*npe+(nhm*(nhm+1)*nat*nspin_mag*npe)/2) )
     mixin=(0.0_DP,0.0_DP)
  ELSE
     ALLOCATE(mixin(1))
  ENDIF
  IF (noncolin) allocate (dbecsum_nc (nhm,nhm, nat , nspin , npe, nsolv))
  allocate (aux2(npwx*npol, nbnd))
  allocate (drhoc(dfftp%nnr))
  IF (noncolin.AND.domag.AND.okvan) THEN
     ALLOCATE (int3_save( nhm, nhm, nat, nspin_mag, npe, 2))
     ALLOCATE (dbecsum_aux ( (nhm * (nhm + 1))/2 , nat , nspin_mag , npe))
  ENDIF
  CALL apply_dpot_allocate()
  !
  if (rec_code_read == 10.AND.ext_recover) then
     ! restart from Phonon calculation
     IF (okpaw) THEN
        CALL read_rec(dr2, iter0, npe, dvscfin, dvscfins, drhoscfh, dbecsum)
        IF (convt) THEN
           CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,npe)
        ELSE
           CALL setmixout(npe*dfftp%nnr*nspin_mag,&
           (nhm*(nhm+1)*nat*nspin_mag*npe)/2,mixin,dvscfin,dbecsum,ndim,-1)
        ENDIF
     ELSE
        CALL read_rec(dr2, iter0, npe, dvscfin, dvscfins, drhoscfh)
     ENDIF
     rec_code=0
  else
    iter0 = 0
    convt =.FALSE.
    where_rec='no_recover'
  endif

  IF (ionode .AND. fildrho /= ' ') THEN
     INQUIRE (UNIT = iudrho, OPENED = exst)
     IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
     filename = dfile_name(xq, at, fildrho, TRIM(tmp_dir_save)//prefix, generate=.true., index_q=iq_dummy)
     CALL diropn (iudrho, filename, lrdrho, exst)
  END IF

  IF (convt) GOTO 155
  !
  ! if q=0 for a metal: allocate and compute local DOS at Ef
  !

  lmetq0 = (lgauss .OR. ltetra) .AND. lgamma
  if (lmetq0) then
     allocate ( ldos ( dfftp%nnr  , nspin_mag) )
     allocate ( ldoss( dffts%nnr , nspin_mag) )
     allocate (becsum1 ( (nhm * (nhm + 1))/2 , nat , nspin_mag))
     call localdos ( ldos , ldoss , becsum1, dos_ef )
     IF (.NOT.okpaw) deallocate(becsum1)
  endif
  !
  !
  ! In this case it has recovered after computing the contribution
  ! to the dynamical matrix. This is a new iteration that has to
  ! start from the beginning.
  !
  IF (iter0==-1000) iter0=0
  !
  ! Compute dV_bare * psi and write to buffer iubar
  !
  DO ik = 1, nksq
     !
     ikk  = ikks(ik)
     ikq  = ikqs(ik)
     npw  = ngk(ikk)
     npwq = ngk(ikq)
     !
     IF (lsda) current_spin = isk(ikk)
     !
     ! compute beta functions for k-point ikq
     !
     CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb)
     !
     DO isolv = 1, nsolv
        IF (isolv == 1) THEN
           ikmk = ikks(ik)
        ELSE
           ikmk = ikmks(ik)
        ENDIF
        !
        ! read unperturbed wavefunctions psi(k) and psi(k+q)
        !
        IF (nksq > 1 .OR. nsolv == 2) THEN
           CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
        ENDIF
        !
        DO ipert = 1, npe
           mode = imode0 + ipert
           nrec = (isolv-1) * npe * nksq + (ipert - 1) * nksq + ik
           !
           IF (isolv==1) THEN
              CALL dvqpsi_us(ik, u(1, mode), .FAlSE., becp1, alphap)
              !
              ! DFPT+U: At the first ph iteration the bare perturbed
              ! Hubbard potential dvbare_hub_q * psi_kpoint
              ! is calculated and added to dvpsi.
              !
              IF (lda_plus_u) CALL dvqhub_barepsi_us(ik, u(1,mode))
              !
           ELSE
              IF (okvan) THEN
                 deeq_nc(:,:,:,:) = deeq_nc_save(:,:,:,:,2)
                 !$acc update device(deeq_nc)
                 int1_nc(:,:,:,:,:) = int1_nc_save(:,:,:,:,:,2)
              ENDIF
              CALL dvqpsi_us(ik, u(1, mode), .FAlSE., becpt, alphapt)
              IF (okvan) THEN
                 deeq_nc(:,:,:,:) = deeq_nc_save(:,:,:,:,1)
                 !$acc update device(deeq_nc)
                 int1_nc(:,:,:,:,:) = int1_nc_save(:,:,:,:,:,1)
              ENDIF
           ENDIF
           !
           CALL save_buffer(dvpsi, lrbar, iubar, nrec)
           !
        ENDDO ! ipert
     ENDDO ! isolv
  ENDDO ! ik
  !
  !   The outside loop is over the iterations
  !
  do kter = 1, niter_ph
     !
     iter = kter + iter0
     !
     first_iter = .NOT. (where_rec == 'solve_lint' .OR. iter > 1)
     !
     drhoscf = (0.d0, 0.d0)
     dbecsum = (0.d0, 0.d0)
     IF (noncolin) dbecsum_nc = (0.d0, 0.d0)
     !
     ! DFPT+U: at each ph iteration calculate dnsscf,
     ! i.e. the scf variation of the occupation matrix ns.
     !
     IF (lda_plus_u .AND. (iter /= 1)) CALL dnsq_scf(npe, lmetq0, imode0, irr, .true.)
     !
     ! Start the loop on the two linear systems, one at B and one at -B
     !
     DO isolv = 1, nsolv
        !
        !  change the sign of the magnetic field if required
        !
        IF (isolv == 2) THEN
           IF (.NOT. first_iter) THEN
              dvscfins(:, 2:4, :) = -dvscfins(:, 2:4, :)
              IF (okvan) int3_nc(:,:,:,:,:) = int3_save(:,:,:,:,:,2)
           ENDIF
           vrs(:, 2:4) = -vrs(:, 2:4)
#if defined(__CUDA)
           vrs_d = vrs
#endif
           IF (okvan) THEN
                   deeq_nc(:,:,:,:) = deeq_nc_save(:,:,:,:,2)
                   !$acc update device(deeq_nc)
           ENDIF
        ENDIF
        !
        ! set threshold for iterative solution of the linear system
        !
        IF (first_iter) THEN
           thresh = 1.0d-2
        ELSE
           thresh = min (1.d-1 * sqrt (dr2), 1.d-2)
        ENDIF
        !
        ! Compute drhoscf, the charge density response to the total potential
        !
        CALL sternheimer_kernel(first_iter, isolv==2, npe, lrbar, iubar, &
            thresh, dvscfins, all_conv, averlt, drhoscf, dbecsum, &
            dbecsum_nc(:,:,:,:,:,isolv))
        !
        !  reset the original magnetic field if it was changed
        !
        IF (isolv == 2) THEN
           IF (.NOT. first_iter) THEN
              dvscfins(:, 2:4, :) = -dvscfins(:, 2:4, :)
              IF (okvan) int3_nc(:,:,:,:,:) = int3_save(:,:,:,:,:,1)
           ENDIF
           vrs(:, 2:4) = -vrs(:, 2:4)
#if defined(__CUDA)
           vrs_d = vrs
#endif
           IF (okvan) THEN
                   deeq_nc(:,:,:,:) = deeq_nc_save(:,:,:,:,1)
                   !$acc update device(deeq_nc)
           ENDIF
        ENDIF
        !
     END DO ! isolv
     !
     IF (nsolv==2) THEN
        drhoscf = drhoscf / 2.0_DP
        dbecsum = dbecsum / 2.0_DP
        dbecsum_nc = dbecsum_nc / 2.0_DP
     ENDIF
     !
     !  The calculation of dbecsum is distributed across processors (see addusdbec)
     !  Sum over processors the contributions coming from each slice of bands
     !
     IF (noncolin) THEN
        call mp_sum ( dbecsum_nc, intra_bgrp_comm )
     ELSE
        call mp_sum ( dbecsum, intra_bgrp_comm )
     ENDIF

     if (doublegrid) then
        do is = 1, nspin_mag
           do ipert = 1, npe
              call fft_interpolate (dffts, drhoscf(:,is,ipert), dfftp, drhoscfh(:,is,ipert))
           enddo
        enddo
     else
        call zcopy (npe*nspin_mag*dfftp%nnr, drhoscf, 1, drhoscfh, 1)
     endif
     !
     !  In the noncolinear, spin-orbit case rotate dbecsum
     !
     IF (noncolin.and.okvan) THEN
        CALL set_dbecsum_nc(dbecsum_nc, dbecsum, npe)
        IF (nsolv==2) THEN
           dbecsum_aux=(0.0_DP,0.0_DP)
           CALL set_dbecsum_nc(dbecsum_nc(1,1,1,1,1,2), dbecsum_aux, npe)
           dbecsum(:,:,1,:)=dbecsum(:,:,1,:)+dbecsum_aux(:,:,1,:)
           dbecsum(:,:,2:4,:)=dbecsum(:,:,2:4,:)-dbecsum_aux(:,:,2:4,:)
        ENDIF
     ENDIF
     !
     !    Now we compute for all perturbations the total charge and potential
     !
     call addusddens (drhoscfh, dbecsum, imode0, npe, 0)
     !
     !   Reduce the delta rho across pools
     !
     call mp_sum ( drhoscf, inter_pool_comm )
     call mp_sum ( drhoscfh, inter_pool_comm )
     IF (okpaw) call mp_sum ( dbecsum, inter_pool_comm )
     !
     IF (okpaw) THEN
        DO ipert=1,npe
           dbecsum(:,:,:,ipert)=2.0_DP *dbecsum(:,:,:,ipert) &
                               +becsumort(:,:,:,imode0+ipert)
        ENDDO
     ENDIF
     !
     ! q=0 in metallic case deserve special care (e_Fermi can shift)
     !
     IF (lmetq0) THEN
        IF (okpaw) THEN
           CALL ef_shift(npe, dos_ef, ldos, drhoscfh, dbecsum, becsum1, irr, sym_def)
        ELSE
           CALL ef_shift(npe, dos_ef, ldos, drhoscfh, irr=irr, sym_def=sym_def)
        ENDIF
     ENDIF
     !
     !   After the loop over the perturbations we have the linear change
     !   in the charge density for each mode of this representation.
     !   Here we symmetrize them ...
     !
     IF (.not.lgamma_gamma) THEN
        call psymdvscf (npe, irr, drhoscfh)
        IF ( noncolin.and.domag ) CALL psym_dmag( npe, irr, drhoscfh)
        IF (okpaw) THEN
           IF (minus_q) CALL PAW_dumqsymmetrize(dbecsum,npe,irr, npertx,irotmq,rtau,xq,tmq)
           CALL PAW_dusymmetrize(dbecsum,npe,irr,npertx,nsymq,rtau,xq,t)
        END IF
     ENDIF
     !
     !   ... save them on disk and
     !   compute the corresponding change in scf potential
     !
     do ipert = 1, npe
        if (fildrho.ne.' ') then
           call davcio_drho (drhoscfh(1,1,ipert), lrdrho, iudrho, imode0+ipert, +1)
!           close(iudrho)
        endif
        !
        call zcopy (dfftp%nnr*nspin_mag,drhoscfh(1,1,ipert),1,dvscfout(1,1,ipert),1)
        !
        ! Compute the response of the core charge density
        ! IT: Should the condition "imode0+ipert > 0" be removed?
        !
        if (imode0+ipert > 0) then
           call addcore (imode0+ipert, drhoc)
        else
           drhoc(:) = (0.0_DP,0.0_DP)
        endif
        !
        ! Compute the response HXC potential
        call dv_of_drho (dvscfout(1,1,ipert), .true., drhoc)
        !
     enddo
     !
     !   And we mix with the old potential
     !
     IF (okpaw) THEN
        !
        !  In this case we mix also dbecsum
        !
        call setmixout(npe*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*npe)/2, &
                    mixout, dvscfout, dbecsum, ndim, -1 )
        call mix_potential (2*npe*dfftp%nnr*nspin_mag+2*ndim, mixout, mixin, &
                         alpha_mix(kter), dr2, npe*tr2_ph/npol, iter, &
                         nmix_ph, flmixdpot, convt)
        call setmixout(npe*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*npe)/2, &
                       mixin, dvscfin, dbecsum, ndim, 1 )
     ELSE
        call mix_potential (2*npe*dfftp%nnr*nspin_mag, dvscfout, dvscfin, &
                         alpha_mix(kter), dr2, npe*tr2_ph/npol, iter, &
                         nmix_ph, flmixdpot, convt)
     ENDIF
     !
     IF (lmetq0 .AND. convt) CALL ef_shift_wfc(npe, ldoss, drhoscf)
     !
     ! check that convergent have been reached on ALL processors in this image
     CALL check_all_convt(convt)

     if (doublegrid) then
        do ipert = 1, npe
           do is = 1, nspin_mag
              call fft_interpolate (dfftp, dvscfin(:,is,ipert), dffts, dvscfins(:,is,ipert))
           enddo
        enddo
     endif
!
!   calculate here the change of the D1-~D1 coefficients due to the phonon
!   perturbation
!
     IF (okvan) THEN
        IF (okpaw) CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,npe)
        !
        !     with the new change of the potential we compute the integrals
        !     of the change of potential and Q
        !
        call newdq (dvscfin, npe)
        !
        !  In the noncollinear magnetic case computes the int3 coefficients with
        !  the opposite sign of the magnetic field. They are saved in int3_save,
        !  that must have been allocated by the calling routine
        !
        IF (noncolin.AND.domag) THEN
           int3_save(:,:,:,:,:,1)=int3_nc(:,:,:,:,:)
           IF (okpaw) rho%bec(:,:,2:4)=-rho%bec(:,:,2:4)
           DO ipert=1,npe
              dvscfin(:,2:4,ipert)=-dvscfin(:,2:4,ipert)
              IF (okpaw) dbecsum(:,:,2:4,ipert)=-dbecsum(:,:,2:4,ipert)
           ENDDO
           !
           !   if needed recompute the paw coeffients with the opposite sign of
           !   the magnetic field
           !
           IF (okpaw) CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,npe)
           !
           !   here compute the int3 integrals
           !
           CALL newdq (dvscfin, npe)
           int3_save(:,:,:,:,:,2)=int3_nc(:,:,:,:,:)
           !
           !  restore the correct sign of the magnetic field.
           !
           DO ipert=1,npe
              dvscfin(:,2:4,ipert)=-dvscfin(:,2:4,ipert)
              IF (okpaw) dbecsum(:,:,2:4,ipert)=-dbecsum(:,:,2:4,ipert)
           ENDDO
           IF (okpaw) rho%bec(:,:,2:4)=-rho%bec(:,:,2:4)
           !
           !  put into int3_nc the coefficient with +B
           !
           int3_nc(:,:,:,:,:)=int3_save(:,:,:,:,:,1)
        ENDIF
     END IF
     !
     tcpu = get_clock ('PHONON')

     WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
          &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
     dr2 = dr2 / npe
     WRITE( stdout, '(5x," thresh=",es10.3, " alpha_mix = ",f6.3, &
          &      " |ddv_scf|^2 = ",es10.3 )') thresh, alpha_mix (kter) , dr2
     !
     !    Here we save the information for recovering the run from this poin
     !
     FLUSH( stdout )
     !
     rec_code=10
     IF (okpaw) THEN
        CALL write_rec('solve_lint', irr, dr2, iter, convt, npe, &
                                               dvscfin, drhoscfh, dbecsum)
     ELSE
        CALL write_rec('solve_lint', irr, dr2, iter, convt, npe, &
                                               dvscfin, drhoscfh)
     ENDIF

     if (check_stop_now()) call stop_smoothly_ph (.false.)
     if (convt) goto 155
  enddo
155 iter0=0
  !
  !    A part of the dynamical matrix requires the integral of
  !    the self consistent change of the potential and the variation of
  !    the charge due to the displacement of the atoms.
  !    We compute it here.
  !
  if (convt) then
     call drhodvus (irr, imode0, dvscfin, npe)
     if (fildvscf.ne.' ') then
        do ipert = 1, npe
           if(lmetq0) then
                dvscfin(:,:,ipert) = dvscfin(:,:,ipert)-def(ipert)
                if (doublegrid) dvscfins(:,:,ipert) = dvscfins(:,:,ipert)-def(ipert)
           endif
           call davcio_drho ( dvscfin(1,1,ipert),  lrdrho, iudvscf, imode0 + ipert, +1 )
           IF (okpaw.AND.ionode) CALL davcio( int3_paw(:,:,:,:,ipert), lint3paw, &
                                              iuint3paw, imode0+ipert, + 1 )
        end do
        if (elph) call elphel (irr, npe, imode0, dvscfins)
     end if
  endif
  if (convt.and.nlcc_any) call addnlcc (imode0, drhoscfh, npe)
  !
  CALL apply_dpot_deallocate()
  if (allocated(ldoss)) deallocate (ldoss)
  if (allocated(ldos)) deallocate (ldos)
  deallocate (dbecsum)
  IF (okpaw) THEN
     if (allocated(becsum1)) deallocate (becsum1)
     deallocate (mixin)
     deallocate (mixout)
  ENDIF
  IF (noncolin) deallocate (dbecsum_nc)
  deallocate (dvscfout)
  deallocate (drhoscfh)
  !$acc exit data delete(dvscfins)
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)
  deallocate(aux2)
  deallocate(drhoc)
  IF (noncolin.AND.domag.AND.okvan) THEN
     DEALLOCATE (int3_save)
     DEALLOCATE (dbecsum_aux)
  ENDIF

  call stop_clock ('solve_linter')

  RETURN

END SUBROUTINE solve_linter

!------------------------------------------------------------------
SUBROUTINE check_all_convt( convt )
  !---------------------------------------------------------------
  !! Work out how many processes have converged.
  !
  USE mp,        ONLY : mp_sum
  USE mp_images, ONLY : nproc_image, me_image, intra_image_comm
  !
  IMPLICIT NONE
  !
  LOGICAL,INTENT(in) :: convt
  INTEGER            :: tot_conv
  !
  IF(nproc_image==1) RETURN
  !
  ! Work out how many processes have converged
  !
  tot_conv = 0
  IF(convt) tot_conv = 1
  CALL mp_sum(tot_conv, intra_image_comm)
  !
  IF ((tot_conv > 0) .and. (tot_conv < nproc_image)) THEN
    CALL errore('check_all_convt', 'Only some processors converged: '&
               &' either something is wrong with solve_linter, or a different'&
               &' parallelism scheme should be used.', 1)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE
