!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE solve_linter (irr, imode0, npe, drhoscf)
  !-----------------------------------------------------------------------
  !
  !    Driver routine for the solution of the linear system which
  !    defines the change of the wavefunction due to a lattice distorsion
  !    It performs the following tasks:
  !     a) computes the bare potential term Delta V | psi >
  !        and an additional term in the case of US pseudopotentials
  !     b) adds to it the screening term Delta V_{SCF} | psi >
  !     c) applies P_c^+ (orthogonalization to valence states)
  !     d) calls cgsolve_all to solve the linear system
  !     e) computes Delta rho, Delta V_{SCF} and symmetrizes them
  !

  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : prefix, diropn
  USE check_stop,           ONLY : check_stop_now
  USE wavefunctions_module, ONLY : evc
  USE constants,            ONLY : degspin
  USE cell_base,            ONLY : at, tpiba2
  USE klist,                ONLY : ltetra, lgauss, degauss, ngauss, &
                                   xk, wk, ngk, igk_k
  USE gvect,                ONLY : g
  USE gvecs,                ONLY : doublegrid
  USE fft_base,             ONLY : dfftp, dffts, dtgs
  USE fft_parallel,         ONLY : tg_cgather
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE spin_orb,             ONLY : domag
  USE wvfct,                ONLY : nbnd, npwx, g2kin,  et
  USE scf,                  ONLY : rho
  USE uspp,                 ONLY : okvan, vkb
  USE uspp_param,           ONLY : upf, nhm, nh
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
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
  USE units_ph,             ONLY : iudrho, lrdrho, iudwf, lrdwf, iubar, lrbar, &
                                   iuwfc, lrwfc, iudvscf, iuint3paw, lint3paw
  USE output,               ONLY : fildrho, fildvscf
  USE phus,                 ONLY : becsumort
  USE modes,                ONLY : npertx, npert, u, t, tmq

  USE recover_mod,          ONLY : read_rec, write_rec
  ! used to write fildrho:
  USE dfile_autoname,       ONLY : dfile_name
  USE save_ph,              ONLY : tmp_dir_save
  ! used oly to write the restart file
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm, ntask_groups, me_bgrp
  USE mp,                   ONLY : mp_sum
  USE efermi_shift,         ONLY : ef_shift, ef_shift_paw,  def

  USE lrus,         ONLY : int3_paw
  USE lr_symm_base, ONLY : irotmq, minus_q, nsymq, rtau
  USE eqv,          ONLY : dvpsi, dpsi, evq
  USE qpoint,       ONLY : xq, nksq, ikks, ikqs
  USE control_lr,   ONLY : alpha_pv, nbnd_occ, lgamma
  USE dv_of_drho_lr

  implicit none

  integer :: irr, npe, imode0
  ! input: the irreducible representation
  ! input: the number of perturbation
  ! input: the position of the modes

  complex(DP) :: drhoscf (dfftp%nnr, nspin_mag, npe)
  ! output: the change of the scf charge

  real(DP) , allocatable :: h_diag (:,:)
  ! h_diag: diagonal part of the Hamiltonian
  real(DP) :: thresh, anorm, averlt, dr2
  ! thresh: convergence threshold
  ! anorm : the norm of the error
  ! averlt: average number of iterations
  ! dr2   : self-consistency error
  real(DP) :: dos_ef, weight, aux_avg (2)
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
       dbecsum (:,:,:,:), dbecsum_nc(:,:,:,:,:), aux1 (:,:), tg_dv(:,:), &
       tg_psic(:,:), aux2(:,:), drhoc(:)
  ! Misc work space
  ! ldos : local density of states af Ef
  ! ldoss: as above, without augmentation charges
  ! dbecsum: the derivative of becsum
  ! drhoc: response core charge density
  REAL(DP), allocatable :: becsum1(:,:,:)

  logical :: conv_root,  & ! true if linear system is converged
             exst,       & ! used to open the recover file
             lmetq0        ! true if xq=(0,0,0) in a metal

  integer :: kter,       & ! counter on iterations
             iter0,      & ! starting iteration
             ipert,      & ! counter on perturbations
             ibnd,       & ! counter on bands
             iter,       & ! counter on iterations
             lter,       & ! counter on iterations of linear system
             ltaver,     & ! average counter
             lintercall, & ! average number of calls to cgsolve_all
             ik, ikk,    & ! counter on k points
             ikq,        & ! counter on k+q points
             ig,         & ! counter on G vectors
             ndim,       &
             is,         & ! counter on spin polarizations
             nrec,       & ! the record number for dvpsi and dpsi
             ios,        & ! integer variable for I/O control
             incr,       & ! used for tg
             v_siz,      & ! size of the potential
             ipol,       & ! counter on polarization
             mode          ! mode index

  integer  :: npw, npwq
  integer  :: iq_dummy
  real(DP) :: tcpu, get_clock ! timing variables
  character(len=256) :: filename

  external ch_psi_all, cg_psi
  !
  IF (rec_code_read > 20 ) RETURN

  call start_clock ('solve_linter')
!
!  This routine is task group aware
!
  allocate (dvscfin ( dfftp%nnr , nspin_mag , npe))
  if (doublegrid) then
     allocate (dvscfins (dffts%nnr , nspin_mag , npe))
  else
     dvscfins => dvscfin
  endif
  allocate (drhoscfh ( dfftp%nnr, nspin_mag , npe))
  allocate (dvscfout ( dfftp%nnr, nspin_mag , npe))
  allocate (dbecsum ( (nhm * (nhm + 1))/2 , nat , nspin_mag , npe))
  IF (okpaw) THEN
     allocate (mixin(dfftp%nnr*nspin_mag*npe+(nhm*(nhm+1)*nat*nspin_mag*npe)/2) )
     allocate (mixout(dfftp%nnr*nspin_mag*npe+(nhm*(nhm+1)*nat*nspin_mag*npe)/2) )
     mixin=(0.0_DP,0.0_DP)
  ENDIF
  IF (noncolin) allocate (dbecsum_nc (nhm,nhm, nat , nspin , npe))
  allocate (aux1 ( dffts%nnr, npol))
  allocate (h_diag ( npwx*npol, nbnd))
  allocate (aux2(npwx*npol, nbnd))
  allocate (drhoc(dfftp%nnr))
  incr=1
  IF ( dtgs%have_task_groups ) THEN
     !
     v_siz =  dtgs%tg_nnr * dtgs%nogrp
     ALLOCATE( tg_dv   ( v_siz, nspin_mag ) )
     ALLOCATE( tg_psic( v_siz, npol ) )
     incr = dtgs%nogrp
     !
  ENDIF
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
     call localdos_paw ( ldos , ldoss , becsum1, dos_ef )
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
  !   The outside loop is over the iterations
  !
  do kter = 1, niter_ph
     iter = kter + iter0

     ltaver = 0

     lintercall = 0
     drhoscf(:,:,:) = (0.d0, 0.d0)
     dbecsum(:,:,:,:) = (0.d0, 0.d0)
     IF (noncolin) dbecsum_nc = (0.d0, 0.d0)
     !
     do ik = 1, nksq
        !
        ikk = ikks(ik)
        ikq = ikqs(ik)
        npw = ngk(ikk)
        npwq= ngk(ikq)
        !
        if (lsda) current_spin = isk (ikk)
        !
        ! read unperturbed wavefunctions psi(k) and psi(k+q)
        !
        if (nksq.gt.1) then
           if (lgamma) then
              call get_buffer (evc, lrwfc, iuwfc, ikk)
           else
              call get_buffer (evc, lrwfc, iuwfc, ikk)
              call get_buffer (evq, lrwfc, iuwfc, ikq)
           endif

        endif
        !
        ! compute beta functions and kinetic energy for k-point ikq
        ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
        !
        CALL init_us_2 (npwq, igk_k(1,ikq), xk (1, ikq), vkb)
        CALL g2_kin (ikq) 
        !
        ! compute preconditioning matrix h_diag used by cgsolve_all
        !
        CALL h_prec (ik, evq, h_diag)
        !
        do ipert = 1, npe
           mode = imode0 + ipert
           nrec = (ipert - 1) * nksq + ik
           !
           !  and now adds the contribution of the self consistent term
           !
           if (where_rec =='solve_lint'.or.iter>1) then
              !
              ! After the first iteration dvbare_q*psi_kpoint is read from file
              !
              call get_buffer (dvpsi, lrbar, iubar, nrec)
              !
              ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
              ! dvscf_q from previous iteration (mix_potential)
              !
              call start_clock ('vpsifft')
              IF( dtgs%have_task_groups ) THEN
                 IF (noncolin) THEN
                    CALL tg_cgather( dffts, dtgs, dvscfins(:,1,ipert), &
                                                                tg_dv(:,1))
                    IF (domag) THEN
                       DO ipol=2,4
                          CALL tg_cgather( dffts, dtgs, dvscfins(:,ipol,ipert), &
                                                             tg_dv(:,ipol))
                       ENDDO
                    ENDIF
                 ELSE
                    CALL tg_cgather( dffts, dtgs, dvscfins(:,current_spin,ipert), &
                                                             tg_dv(:,1))
                 ENDIF
              ENDIF
              aux2=(0.0_DP,0.0_DP)
              do ibnd = 1, nbnd_occ (ikk), incr
                 IF( dtgs%have_task_groups ) THEN
                    call cft_wave_tg (ik, evc, tg_psic, 1, v_siz, ibnd, &
                                      nbnd_occ (ikk) )
                    call apply_dpot(v_siz, tg_psic, tg_dv, 1)
                    call cft_wave_tg (ik, aux2, tg_psic, -1, v_siz, ibnd, &
                                      nbnd_occ (ikk))
                 ELSE
                    call cft_wave (ik, evc (1, ibnd), aux1, +1)
                    call apply_dpot(dffts%nnr,aux1, dvscfins(1,1,ipert), current_spin)
                    call cft_wave (ik, aux2 (1, ibnd), aux1, -1)
                 ENDIF
              enddo
              dvpsi=dvpsi+aux2
              call stop_clock ('vpsifft')
              !
              !  In the case of US pseudopotentials there is an additional
              !  selfconsist term which comes from the dependence of D on
              !  V_{eff} on the bare change of the potential
              !
              call adddvscf (ipert, ik)
           else
              !
              !  At the first iteration dvbare_q*psi_kpoint is calculated
              !  and written to file
              !
              call dvqpsi_us (ik, u (1, mode),.false. )
              call save_buffer (dvpsi, lrbar, iubar, nrec)
           endif
           !
           ! Ortogonalize dvpsi to valence states: ps = <evq|dvpsi>
           ! Apply -P_c^+.
           !
           CALL orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, .false.)
           !
           if (where_rec=='solve_lint'.or.iter > 1) then
              !
              ! starting value for delta_psi is read from iudwf
              !
              call get_buffer( dpsi, lrdwf, iudwf, nrec)
              !
              ! threshold for iterative solution of the linear system
              !
              thresh = min (1.d-1 * sqrt (dr2), 1.d-2)
           else
              !
              !  At the first iteration dpsi and dvscfin are set to zero
              !
              dpsi(:,:) = (0.d0, 0.d0)
              dvscfin (:, :, ipert) = (0.d0, 0.d0)
              !
              ! starting threshold for iterative solution of the linear system
              !
              thresh = 1.0d-2
           endif

           !
           ! iterative solution of the linear system (H-eS)*dpsi=dvpsi,
           ! dvpsi=-P_c^+ (dvbare+dvscf)*psi , dvscf fixed.
           !
           conv_root = .true.

           call cgsolve_all (ch_psi_all, cg_psi, et(1,ikk), dvpsi, dpsi, &
                             h_diag, npwx, npwq, thresh, ik, lter, conv_root, &
                             anorm, nbnd_occ(ikk), npol )

           ltaver = ltaver + lter
           lintercall = lintercall + 1
           if (.not.conv_root) WRITE( stdout, '(5x,"kpoint",i4," ibnd",i4,  &
                &              " solve_linter: root not converged ",es10.3)') &
                &              ik , ibnd, anorm
           !
           ! writes delta_psi on iunit iudwf, k=kpoint,
           !
           !               if (nksq.gt.1 .or. npert(irr).gt.1)
           call save_buffer (dpsi, lrdwf, iudwf, nrec)
           !
           ! calculates dvscf, sum over k => dvscf_q_ipert
           !
           weight = wk (ikk)
           IF (noncolin) THEN
              call incdrhoscf_nc(drhoscf(1,1,ipert),weight,ik, &
                                 dbecsum_nc(1,1,1,1,ipert), dpsi)
           ELSE
              call incdrhoscf (drhoscf(1,current_spin,ipert), weight, ik, &
                            dbecsum(1,1,current_spin,ipert), dpsi)
           END IF
           ! on perturbations
        enddo
        ! on k-points
     enddo
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
              call cinterpolate (drhoscfh(1,is,ipert), drhoscf(1,is,ipert), 1)
           enddo
        enddo
     else
        call zcopy (npe*nspin_mag*dfftp%nnr, drhoscf, 1, drhoscfh, 1)
     endif
     !
     !  In the noncolinear, spin-orbit case rotate dbecsum
     !
     IF (noncolin.and.okvan) CALL set_dbecsum_nc(dbecsum_nc, dbecsum, npe)
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
     ! q=0 in metallic case deserve special care (e_Fermi can shift)
     !

     IF (okpaw) THEN
        IF (lmetq0) &
           call ef_shift_paw (drhoscfh, dbecsum, ldos, ldoss, becsum1, &
                                                  dos_ef, irr, npe, .false.)
        DO ipert=1,npe
           dbecsum(:,:,:,ipert)=2.0_DP *dbecsum(:,:,:,ipert) &
                               +becsumort(:,:,:,imode0+ipert)
        ENDDO
     ELSE
        IF (lmetq0) call ef_shift(drhoscfh,ldos,ldoss,dos_ef,irr,npe,.false.)
     ENDIF
     !
     !   After the loop over the perturbations we have the linear change
     !   in the charge density for each mode of this representation.
     !   Here we symmetrize them ...
     !
     IF (.not.lgamma_gamma) THEN
        call psymdvscf (npe, irr, drhoscfh)
        IF ( noncolin.and.domag ) &
           CALL psym_dmag( npe, irr, drhoscfh)
        IF (okpaw) THEN
           IF (minus_q) CALL PAW_dumqsymmetrize(dbecsum,npe,irr, &
                             npertx,irotmq,rtau,xq,tmq)
           CALL  &
              PAW_dusymmetrize(dbecsum,npe,irr,npertx,nsymq,rtau,xq,t)
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
        call mix_potential (2*npe*dfftp%nnr*nspin_mag+2*ndim, &
                         mixout, mixin, &
                         alpha_mix(kter), dr2, npe*tr2_ph/npol, iter, &
                         nmix_ph, flmixdpot, convt)
        call setmixout(npe*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*npe)/2, &
                       mixin, dvscfin, dbecsum, ndim, 1 )
        if (lmetq0.and.convt) &
           call ef_shift_paw (drhoscf, dbecsum, ldos, ldoss, becsum1, &
                                                  dos_ef, irr, npe, .true.)
     ELSE
        call mix_potential (2*npe*dfftp%nnr*nspin_mag, dvscfout, dvscfin, &
                         alpha_mix(kter), dr2, npe*tr2_ph/npol, iter, &
                         nmix_ph, flmixdpot, convt)
        if (lmetq0.and.convt) &
            call ef_shift (drhoscf, ldos, ldoss, dos_ef, irr, npe, .true.)
     ENDIF
     ! check that convergent have been reached on ALL processors in this image
     CALL check_all_convt(convt)

     if (doublegrid) then
        do ipert = 1, npe
           do is = 1, nspin_mag
              call cinterpolate (dvscfin(1,is,ipert), dvscfins(1,is,ipert), -1)
           enddo
        enddo
     endif
!
!   calculate here the change of the D1-~D1 coefficients due to the phonon
!   perturbation
!
     IF (okpaw) CALL PAW_dpotential(dbecsum,rho%bec,int3_paw,npe)
     !
     !     with the new change of the potential we compute the integrals
     !     of the change of potential and Q
     !
     call newdq (dvscfin, npe)
#if defined(__MPI)
     aux_avg (1) = DBLE (ltaver)
     aux_avg (2) = DBLE (lintercall)
     call mp_sum ( aux_avg, inter_pool_comm )
     averlt = aux_avg (1) / aux_avg (2)
#else
     averlt = DBLE (ltaver) / lintercall
#endif
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
           call davcio_drho ( dvscfin(1,1,ipert),  lrdrho, iudvscf, &
                         imode0 + ipert, +1 )
           IF (okpaw.AND.me_bgrp==0) CALL davcio( int3_paw(:,:,:,:,ipert), lint3paw, &
                                                  iuint3paw, imode0+ipert, + 1 )
        end do
        if (elph) call elphel (irr, npe, imode0, dvscfins)
     end if
  endif
  if (convt.and.nlcc_any) call addnlcc (imode0, drhoscfh, npe)
  if (allocated(ldoss)) deallocate (ldoss)
  if (allocated(ldos)) deallocate (ldos)
  deallocate (h_diag)
  deallocate (aux1)
  deallocate (dbecsum)
  IF (okpaw) THEN
     if (allocated(becsum1)) deallocate (becsum1)
     deallocate (mixin)
     deallocate (mixout)
  ENDIF
  IF (noncolin) deallocate (dbecsum_nc)
  deallocate (dvscfout)
  deallocate (drhoscfh)
  if (doublegrid) deallocate (dvscfins)
  deallocate (dvscfin)
  deallocate(aux2)
  deallocate(drhoc)
  IF ( dtgs%have_task_groups ) THEN
     DEALLOCATE( tg_dv )
     DEALLOCATE( tg_psic )
  ENDIF

  call stop_clock ('solve_linter')
END SUBROUTINE solve_linter


SUBROUTINE setmixout(in1, in2, mix, dvscfout, dbecsum, ndim, flag )
USE kinds, ONLY : DP
USE mp_bands, ONLY : intra_bgrp_comm
USE mp, ONLY : mp_sum
IMPLICIT NONE
INTEGER :: in1, in2, flag, ndim, startb, lastb
COMPLEX(DP) :: mix(in1+in2), dvscfout(in1), dbecsum(in2)

CALL divide (intra_bgrp_comm, in2, startb, lastb)
ndim=lastb-startb+1

IF (flag==-1) THEN
   mix(1:in1)=dvscfout(1:in1)
   mix(in1+1:in1+ndim)=dbecsum(startb:lastb)
ELSE
   dvscfout(1:in1)=mix(1:in1)
   dbecsum=(0.0_DP,0.0_DP)
   dbecsum(startb:lastb)=mix(in1+1:in1+ndim)
   CALL mp_sum(dbecsum, intra_bgrp_comm)
ENDIF
END SUBROUTINE setmixout

SUBROUTINE check_all_convt(convt)
  USE mp,        ONLY : mp_sum
  USE mp_images, ONLY : nproc_image, me_image, intra_image_comm
  IMPLICIT NONE
  LOGICAL,INTENT(in) :: convt
  INTEGER,ALLOCATABLE :: convt_check(:)
  !
  IF(nproc_image==1) RETURN
  !
  ALLOCATE(convt_check(nproc_image+1))
  !
  convt_check = 1
  IF(convt) convt_check(me_image+1) = 0
  !
  CALL mp_sum(convt_check, intra_image_comm)
  !CALL mp_sum(ios, inter_pool_comm)
  !CALL mp_sum(ios, intra_bgrp_comm)
  !
!  convt = ALL(convt_check==0)
  IF(ANY(convt_check==0).and..not.ALL(convt_check==0) ) THEN
    CALL errore('check_all_convt', 'Only some processors converged: '&
               &' something is wrong with solve_linter', 1)
  ENDIF
  !
  DEALLOCATE(convt_check)
  RETURN
  !
END SUBROUTINE
