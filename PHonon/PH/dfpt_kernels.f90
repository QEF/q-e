!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------------------------
MODULE dfpt_kernels
CONTAINS
!
SUBROUTINE dfpt_kernel(code, npert, iter0, lrdvpsi, iudvpsi, dr2, drhos, drhop, dvscfs, dvscfp, dbecsum, &
                       irr, imode0, option, drhoc)
   !------------------------------------------------------------------------------
   !! Driver routine for the solution of the linear system that computes the change
   !! of the electron density due to a generic perturbation to the Hamiltonian.
   !!
   !! This routine assumes the follwing input:
   !! a) the bare potential term multiplied to the wavefunction \(\Delta V | \psi \rangle\)
   !!    and an additional term in the case of US pseudopotentials.
   !!    For DFPT+U also the bare perturbation of the Hubbard potential should be included.
   !!    This data is stored in the buffer iudvpsi with length lrdvpsi.
   !!
   !! Then, it performs the following tasks:
   !! b) compute the wavefunction response \(\Delta \psi\) and the charge density response
   !!    \(\Delta \rho\text{SCF}\) by calling sternheimer_kernel.
   !!    Inside sternheimer_kernel, the following tasks are performed:
   !!       b-1) adds the screening term \(\Delta V_\text{SCF} | \psi \rangle\).
   !!       b-2) applies \(P_c^+\) (orthogonalization to valence states);
   !!       b-3) CALLs \(\text{cgsolve_all}\) to solve the linear system;
   !!       b-4) compute the
   !! c) For the USPP case, add additional contribution to \(\Delta\rho\text{SCF}\).
   !! d) Symmetrize the charge density response \(\Delta\rho\text{SCF}\).
   !! e) Compute \(\Delta V_\text{SCF}\) from \(\Delta\rho\text{SCF}\).
   !! f) If \(\text{lda_plus_u}=\text{TRUE}\) compute also the response
   !!    occupation matrices \(\text{dnsscf}\);
   !! g) Mix the potentials;
   !! h) If converged, return. Otherwise, return to step b.
   !!
   !! For systems with noncollinear magnetism (noncolin = true and domag = true), the linear
   !! system is solved twice (nsolv = 2).
   !! The first case \(\text{isolv}=1\) is the ordinary one, while the second case
   !! \(\text{isolv}=2\) needs the time-reversed wave functions. For the theoretical
   !! background, please refer to:
   !!    Phys. Rev. B 100, 045115 (2019).
   !!
   !! Input:
   !!    - code : Name of the code calling this subroutine
   !!    - npert : Number of perturbations
   !!    - iter0 : Number of iterations already performed (> 0 for restart, 0 otherwise)
   !!    - lrdvpsi : Record length for the buffer storing dV_bare * psi
   !!    - iudvpsi : Unit for the buffer storing dV_bare * psi
   !!    - option : Option that tells the type of perturbation. phonon / efield / ...
   !!
   !! Input/Output:
   !!    - drhos : change of the charge density (smooth part only, dffts)
   !!    - drhop : change of the charge density (smooth and hard parts, dfftp)
   !!    - dvscfs : change of the scf potential (smooth part only, dffts)
   !!    - dvscfp : change of the scf potential (smooth and hard parts, dfftp)
   !!    - dbecsum : change of the becsum (used if USPP or PAW)
   !! (When restarting, nonzero input is used. When starting from scratch, zero input is used.)
   !!
   !! Output (stored in modules, not in argument):
   !!    - dpsi : Wavefunction perturbation dpsi. Stored in buffer iudwf with length lrdwf
   !!    - def : Fermi energy shift (only for metals at q=0)
   !!
   !! variable "*s" : real-space quantity defined on the soft grid (dffts)
   !! variable "*h" : real-space quantity defined on the hard grid (dfftp)
   !! (If doublegrid == false, the two quantities are identical.)
   !!
   !! A short summary of the workflow:
   !!    dvscfs, dvscfp    -> (sternheimer_kernel)
   !! -> drhos,  drhop     -> (dv_of_drho)
   !! -> dvscftmp          -> (mix_potential)
   !! -> dvscfs, dvscfp    -> ...
   !!
   !! NOTE: The following features are not implemented.
   !!    - Custom convergence parameter
   !!    - Custom printing
   !!
   !!
   !! Currently the code needs branches for phonon / e-field (using variable option)
   !! in the following places
   !!    1. addusddens / addusddense : USPP contribution to drho
   !!    2. becsumort : alphasum contribution to dbecsum
   !!    3. PAW related stuff. symmetrization, factor of 2, ...
   !! The plan is to get rid of all these branches by designing a generic subroutine, or
   !! if that is not feasible, by using callback arguments.
   !!
   !------------------------------------------------------------------------------
   !
   USE kinds,                ONLY : DP
   USE mp_pools,             ONLY : inter_pool_comm
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   USE fft_interfaces,       ONLY : fft_interpolate
   USE fft_base,             ONLY : dfftp, dffts
   USE ions_base,            ONLY : nat
   USE io_global,            ONLY : stdout
   USE io_files,             ONLY : diropn
   USE check_stop,           ONLY : check_stop_now
   USE klist,                ONLY : ltetra, lgauss
   USE gvecs,                ONLY : doublegrid
   USE lsda_mod,             ONLY : nspin, lsda
   USE scf,                  ONLY : rho
   USE uspp,                 ONLY : okvan
   USE uspp_param,           ONLY : nhm
   USE uspp_init,            ONLY : init_us_2
   USE noncollin_module,     ONLY : noncolin, domag, npol, nspin_mag
   USE paw_variables,        ONLY : okpaw
   USE paw_onecenter,        ONLY : paw_dpotential
   USE paw_symmetry,         ONLY : paw_dusymmetrize, paw_dumqsymmetrize, paw_desymmetrize
   USE phus,                 ONLY : becsumort
   USE modes,                ONLY : npertx, t, tmq
   USE recover_mod,          ONLY : write_rec
   USE efermi_shift,         ONLY : ef_shift, ef_shift_wfc, def
   USE lrus,                 ONLY : int3_paw, int3_nc
   USE lr_symm_base,         ONLY : irotmq, minus_q, nsymq, rtau
   USE qpoint,               ONLY : xq
   USE control_ph,           ONLY : lnoloc
   USE control_lr,           ONLY : lgamma, niter_ph, nmix_ph, tr2_ph, alpha_mix, convt, &
                                    lgamma_gamma, flmixdpot, where_rec, lmultipole
   USE dv_of_drho_lr,        ONLY : dv_of_drho
   USE ldaU,                 ONLY : lda_plus_u
   USE lr_nc_mag,            ONLY : int3_nc_save
   USE apply_dpot_mod,       ONLY : apply_dpot_allocate, apply_dpot_deallocate
   USE response_kernels,     ONLY : sternheimer_kernel
   USE two_chem,             ONLY : twochem
   USE lr_two_chem,          ONLY : allocate_twochem, deallocate_twochem
   !
   IMPLICIT NONE
   !
   CHARACTER(LEN=*), INTENT(IN) :: code
   !! Name of the code calling this subroutine
   INTEGER, INTENT(IN) :: npert
   !! Number of perturbations
   INTEGER, INTENT(IN) :: iter0
   !! Number of iterations already performed (> 0 for restart, 0 otherwise)
   INTEGER, INTENT(IN) :: lrdvpsi
   !! record length for the buffer storing dV_bare * psi
   INTEGER, INTENT(IN) :: iudvpsi
   !! unit for the buffer storing dV_bare * psi
   REAL(DP), INTENT(INOUT) :: dr2
   ! self-consistency error. Input is used for restart.
   COMPLEX(DP), INTENT(INOUT) :: drhos(dffts%nnr, nspin_mag, npert)
   !! change of the charge density (smooth part only, but allocated with dfftp)
   COMPLEX(DP), INTENT(INOUT) :: drhop(dfftp%nnr, nspin_mag, npert)
   !! change of the charge density (smooth and hard parts, dfftp)
   COMPLEX(DP), POINTER, INTENT(INOUT) :: dvscfs(:, :, :)
   !! change of the scf potential (smooth part only, dffts)
   COMPLEX(DP), INTENT(INOUT) :: dvscfp(dfftp%nnr, nspin_mag, npert)
   !! change of the scf potential (smooth and hard parts, dfftp)
   COMPLEX(DP), INTENT(INOUT) :: dbecsum((nhm * (nhm + 1))/2, nat, nspin_mag , npert)
   !! change of becsum
   COMPLEX(DP), INTENT(IN), OPTIONAL :: drhoc(dfftp%nnr, npert)
   !! Change in the core charge due to the perturbation.
   !
   INTEGER, INTENT(IN) :: irr
   !! input: the irreducible representation (to be removed after moving from PH to LR_Modules)
   !! Set to 1 if not using phonon perturbation.
   INTEGER, INTENT(IN) :: imode0
   !! input: the position of the modes (to be removed after moving from PH to LR_Modules)
   !! Set to 0 if not using phonon perturbation.
   CHARACTER(LEN=*), INTENT(IN) :: option
   !! input: Option that tells the type of perturbation.
   !! Possible options: 'phonon', 'efield'
   !
   ! ... local variables
   !
   REAL(DP) :: thresh, averlt
   ! thresh: convergence threshold
   ! averlt: average number of iterations
   REAL(DP) :: dos_ef
   !! dos_ef: density of states at Ef
   COMPLEX(DP), ALLOCATABLE :: dvscftmp(:, :, :)
   !! change of scf potential (output before mixing)
   COMPLEX(DP), ALLOCATABLE :: ldos (:,:), ldoss (:,:), mixin(:), mixout(:), &
      dbecsum_nc(:,:,:,:,:,:), dbecsum_aux (:,:,:,:)
   ! Misc work space
   ! ldos : local density of states af Ef
   ! ldoss: as above, without augmentation charges
   ! dbecsum: the derivative of becsum
   REAL(DP), ALLOCATABLE :: becsum1(:,:,:)

   LOGICAL :: all_conv
   !! True if sternheimer_kernel is converged at all k points and perturbations
   logical :: lmetq0,     & ! true if xq=(0,0,0) in a metal
      first_iter    ! true if first iteration where induced rho is not yet calculated

   integer :: kter,       & ! counter on iterations
      ipert,      & ! counter on perturbations
      iter,       & ! counter on iterations
      ndim,       &
      ndim_pot,   &
      ndim_paw,   &
      is,         & ! counter on spin polarizations
      isolv,      & ! counter on linear systems
      nsolv         ! number of linear systems
   !
   REAL(DP) :: tcpu, get_clock
   !! timing variables
   !
   !  This routine is task group aware
   !
   CALL start_clock('dfpt_kernel')
   !
   IF (option /= 'phonon' .AND. option /= 'efield') THEN
      CALL errore('dfpt_kernel', 'Unknown option' // TRIM(option), 1)
   ENDIF
   !
   nsolv = 1
   IF (noncolin .AND. domag) nsolv=2
   !
   CALL apply_dpot_allocate()
   !
   ALLOCATE(dvscftmp(dfftp%nnr, nspin_mag, npert))
   IF (noncolin) ALLOCATE(dbecsum_nc(nhm, nhm, nat, nspin, npert, nsolv))
   IF (noncolin .AND. domag .AND. okvan) THEN
      ALLOCATE(int3_nc_save(nhm, nhm, nat, nspin_mag, npert, 2))
      ALLOCATE(dbecsum_aux((nhm * (nhm + 1))/2, nat, nspin_mag, npert))
   ENDIF
   !
   ! Allocate temporary variables for mixing
   !
   ndim_pot = npert * dfftp%nnr * nspin_mag
   !
   IF (okpaw) THEN
      ndim_paw = (nhm * (nhm+1) * nat * nspin_mag * npert) / 2
      ALLOCATE(mixin(ndim_pot + ndim_paw))
      ALLOCATE(mixout(ndim_pot + ndim_paw))
      !
      ! PAW: Set mixin from the intial dvscfp and dbecsum. Used as the input for mixing.
      CALL setmixout(ndim_pot, ndim_paw, mixin, dvscfp, dbecsum, ndim, -1)
      !
   ELSE
      ALLOCATE(mixin(1))
      ALLOCATE(mixout(1))
   ENDIF
   !
   ! if q=0 for a metal: allocate and compute local DOS at Ef
   !
   lmetq0 = (lgauss .OR. ltetra) .AND. lgamma
   IF (lmetq0) THEN
      ALLOCATE(ldos(dfftp%nnr, nspin_mag))
      ALLOCATE(ldoss(dffts%nnr, nspin_mag))
      ALLOCATE(becsum1((nhm * (nhm + 1))/2 , nat , nspin_mag))
      CALL localdos(ldos, ldoss, becsum1, dos_ef)
      IF (.NOT. okpaw) DEALLOCATE(becsum1)
   endif
   !
   IF (twochem) CALL allocate_twochem(npert, nsolv)
   !
   ! Loop over DFPT self-consistent iterations
   !
   do kter = 1, niter_ph
      !
      iter = kter + iter0
      !
      first_iter = (iter == 1)
      !
      drhos = (0.d0, 0.d0)
      dbecsum = (0.d0, 0.d0)
      IF (noncolin) dbecsum_nc = (0.d0, 0.d0)
      !
      ! DFPT+U: at each ph iteration calculate dnsscf,
      ! i.e. the scf variation of the occupation matrix ns.
      !
      IF (lda_plus_u .AND. (.NOT. first_iter)) CALL dnsq_scf(npert, lmetq0)
      !
      ! Start the loop on the two linear systems, one at B and one at -B
      !
      DO isolv = 1, nsolv
         !
         ! Set threshold for iterative solution of the linear system
         !
         IF (first_iter) THEN
            thresh = 1.0d-2
         ELSE
            thresh = min(1.d-1 * sqrt(dr2), 1.d-2)
         ENDIF
         !
         ! Compute drhos, the charge density response to the total potential
         !
         CALL sternheimer_kernel(first_iter, isolv==2, npert, lrdvpsi, iudvpsi, &
            thresh, dvscfs, all_conv, averlt, drhos, dbecsum, &
            dbecsum_nc(:,:,:,:,:,isolv))
         !
      ENDDO ! isolv
      !
      IF (nsolv == 2) THEN
         drhos = drhos / 2.0_DP
         dbecsum = dbecsum / 2.0_DP
         dbecsum_nc = dbecsum_nc / 2.0_DP
      ENDIF
      !
      !  The calculation of dbecsum is distributed across processors (see addusdbec)
      !  Sum over processors the contributions coming from each slice of bands
      !
      IF (noncolin) THEN
         CALL mp_sum(dbecsum_nc, intra_bgrp_comm)
      ELSE
         CALL mp_sum(dbecsum, intra_bgrp_comm)
      ENDIF
      !
      IF (doublegrid) THEN
         DO is = 1, nspin_mag
            DO ipert = 1, npert
               CALL fft_interpolate(dffts, drhos(:, is, ipert), dfftp, drhop(:, is, ipert))
            ENDDO
         ENDDO
      ELSE
         CALL zcopy(npert * nspin_mag * dfftp%nnr, drhos, 1, drhop, 1)
      ENDIF
      !
      !  In the noncolinear, spin-orbit case rotate dbecsum
      !
      IF (noncolin .AND. okvan) THEN
         CALL set_dbecsum_nc(dbecsum_nc, dbecsum, npert)
         IF (nsolv == 2) THEN
            dbecsum_aux = (0.0_DP, 0.0_DP)
            CALL set_dbecsum_nc(dbecsum_nc(1,1,1,1,1,2), dbecsum_aux, npert)
            dbecsum(:,:,1,:) = dbecsum(:,:,1,:) + dbecsum_aux(:,:,1,:)
            dbecsum(:,:,2:4,:) = dbecsum(:,:,2:4,:) - dbecsum_aux(:,:,2:4,:)
         ENDIF
      ENDIF
      !
      !    Now we compute for all perturbations the total charge and potential
      !
      IF (option == 'phonon') THEN
         CALL addusddens(drhop, dbecsum, imode0, npert, 0)
      ELSEIF (option == 'efield') THEN
         call addusddense(drhop, dbecsum)
      ELSE
         CALL errore('dfpt_kernel', 'Unknown option' // TRIM(option), 1)
      ENDIF
      !
      !   Reduce the delta rho across pools
      !
      CALL mp_sum(drhos, inter_pool_comm)
      CALL mp_sum(drhop, inter_pool_comm)
      IF (okpaw) CALL mp_sum(dbecsum, inter_pool_comm)
      !
      IF (okpaw) THEN
         IF (option == 'phonon') THEN
            ! For efield, factor 2 is multiplied after mixing
            DO ipert = 1, npert
               dbecsum(:,:,:,ipert) = 2.0_DP * dbecsum(:,:,:,ipert) &
                  + becsumort(:,:,:,imode0+ipert)
            ENDDO
         ENDIF
      ENDIF
      !
      ! Metallic case and q=0: add a correction to the response charge density
      ! due to the shift of the Fermi energy (see Eq.(75) in Rev. Mod. Phys. 73, 515 (2001)).
      ! This term is added to the response charge density (in order to obtain correct
      ! response HXC potential) and to the response occupation matrices.
      !
      IF (lmetq0) THEN
         IF (okpaw) THEN
            CALL ef_shift(npert, dos_ef, ldos, drhop, dbecsum, becsum1)
         ELSE
            CALL ef_shift(npert, dos_ef, ldos, drhop)
         ENDIF
      ENDIF
      !
      ! Repeat the above for the case of two chemical potentials
      !
      IF (twochem) THEN
         IF (okpaw) THEN
            CALL twochem_postproc_dfpt(npert, nsolv, imode0, lmetq0, &
                  convt, dos_ef, ldos, ldoss, drhop, dbecsum, becsum1)
         ELSE
            CALL twochem_postproc_dfpt(npert, nsolv, imode0, lmetq0, &
                  convt, dos_ef, ldos, ldoss, drhop, dbecsum)
         ENDIF
      ENDIF
      !
      !   After the loop over the perturbations we have the linear change
      !   in the charge density for each mode of this representation.
      !   Here we symmetrize them ...
      !
      IF (.NOT. lgamma_gamma) THEN
         CALL psymdvscf(drhop)
         IF (lmultipole) THEN !FM: density computed for the first representation only, needs to be symmetrized
            IF (doublegrid) then
               DO is = 1, nspin_mag
                  DO ipert = 1, npert
                     CALL fft_interpolate(dfftp, drhop(:, is, ipert), dffts, drhos(:, is, ipert))
                  ENDDO
               ENDDO
            ELSE
              CALL zcopy(npert * nspin_mag * dfftp%nnr, drhop, 1, drhos, 1)
            ENDIF
         ENDIF
         !
         IF (okpaw) THEN
            IF (option == 'phonon') THEN
               ! For efield, PAW symmetrization is done after mixing
               IF (minus_q) CALL PAW_dumqsymmetrize(dbecsum,npert,irr, npertx,irotmq,rtau,xq,tmq)
               CALL PAW_dusymmetrize(dbecsum,npert,irr,npertx,nsymq,rtau,xq,t)
            ENDIF
         ENDIF
      ENDIF
      !
      !   compute the corresponding change in scf potential : drhop -> dvscftmp
      !
      DO ipert = 1, npert
         !
         CALL zcopy(dfftp%nnr*nspin_mag, drhop(1,1,ipert), 1, dvscftmp(1,1,ipert), 1)
         !
         ! Compute the response HXC potential
         !
         IF (lnoloc) THEN
            ! No local field effect: set dvscf to 0
            dvscftmp(:, :, ipert) = (0.d0, 0.d0)
         ELSE
            IF (PRESENT(drhoc) .AND. .NOT. lmultipole) THEN
               CALL dv_of_drho(dvscftmp(1, 1, ipert), drhoc = drhoc(:, ipert))
            ELSE !FM: as the case of solve_e
               CALL dv_of_drho(dvscftmp(1, 1, ipert))
            ENDIF
         ENDIF
         !
      ENDDO
      !
      !   And we mix with the old potential: dvscftmp -> dvscfp
      !
      IF (okpaw) THEN
         !
         !  In this case we mix also dbecsum
         !
         CALL setmixout(ndim_pot, ndim_paw, mixout, dvscftmp, dbecsum, ndim, -1)
         CALL mix_potential(2 * (ndim_pot + ndim), mixout, mixin, alpha_mix(kter), dr2, &
                            npert * tr2_ph / npol, iter, nmix_ph, flmixdpot, convt)
         CALL setmixout(ndim_pot, ndim_paw, mixin, dvscfp, dbecsum, ndim, 1)
      ELSE
         CALL mix_potential(2 * ndim_pot, dvscftmp, dvscfp, alpha_mix(kter), dr2, &
                            npert * tr2_ph / npol, iter, nmix_ph, flmixdpot, convt)
      ENDIF
      !
      ! Check that convergence have been reached on ALL processors in this image.
      ! If convergence is reached on only some processors, call errore.
      !
      CALL check_all_convt(convt)
      !
      IF (doublegrid) THEN
         DO ipert = 1, npert
            DO is = 1, nspin_mag
               CALL fft_interpolate(dfftp, dvscfp(:, is, ipert), dffts, dvscfs(:, is, ipert))
            ENDDO
         ENDDO
      ENDIF
      !
      IF (option == 'efield') THEN
         IF (okpaw) THEN
            IF (noncolin) THEN
               ! call PAW_dpotential(dbecsum_nc, becsum_nc, int3_paw, npert)
            ELSE
               !
               ! The presence of c.c. in the formula gives a factor 2.0
               ! For phonons, this factor was multiplied before mixing
               !
               dbecsum = 2.0_DP * dbecsum
               IF (.NOT. lgamma_gamma) CALL PAW_desymmetrize(dbecsum)
            ENDIF
         ENDIF
      ENDIF
      !
      !   calculate here the change of the D1-~D1 coefficients due to the phonon
      !   perturbation
      !
      IF (okvan) THEN
         IF (okpaw) CALL PAW_dpotential(dbecsum, rho%bec, int3_paw, npert)
         !
         !     with the new change of the potential we compute the integrals
         !     of the change of potential and Q
         !
         CALL newdq(dvscfp, npert)
         !
         !  In the noncollinear magnetic case computes the int3 coefficients with
         !  the opposite sign of the magnetic field. They are saved in int3_nc_save,
         !  that must have been ALLOCATEd by the calling routine
         !
         IF (noncolin .AND. domag) THEN
            int3_nc_save(:,:,:,:,:,1) = int3_nc(:,:,:,:,:)
            IF (okpaw) rho%bec(:,:,2:4) = -rho%bec(:,:,2:4)
            DO ipert = 1, npert
               dvscfp(:, 2:4, ipert) = -dvscfp(:, 2:4, ipert)
               IF (okpaw) dbecsum(:, :, 2:4, ipert) = -dbecsum(:, :, 2:4, ipert)
            ENDDO
            !
            !   if needed recompute the paw coeffients with the opposite sign of
            !   the magnetic field
            !
            IF (okpaw) CALL PAW_dpotential(dbecsum, rho%bec, int3_paw, npert)
            !
            !   here compute the int3 integrals
            !
            CALL newdq(dvscfp, npert)
            int3_nc_save(:,:,:,:,:,2) = int3_nc(:,:,:,:,:)
            !
            !  restore the correct sign of the magnetic field.
            !
            DO ipert = 1, npert
               dvscfp(:, 2:4, ipert) = -dvscfp(:, 2:4, ipert)
               IF (okpaw) dbecsum(:, :, 2:4, ipert) = -dbecsum(:, :, 2:4, ipert)
            ENDDO
            IF (okpaw) rho%bec(:, :, 2:4) = -rho%bec(:, :, 2:4)
            !
            !  put into int3_nc the coefficient with +B
            !
            int3_nc(:,:,:,:,:)=int3_nc_save(:,:,:,:,:,1)
         ENDIF ! noncolin .AND. domag
      ENDIF ! okvan
      !
      tcpu = get_clock(code)
      !
      dr2 = dr2 / npert
      WRITE( stdout, '(/,5x," iter # ",i3," total cpu time :",f8.1, &
      &      " secs   av.it.: ",f5.1)') iter, tcpu, averlt
      WRITE( stdout, '(5x," thresh=",es10.3, " alpha_mix = ",f6.3, &
      &      " |ddv_scf|^2 = ",es10.3 )') thresh, alpha_mix(kter) , dr2
      FLUSH(stdout)
      !
      !    Here we save the information for recovering the run from this poin
      !
      IF (okpaw) THEN
         CALL write_rec(where_rec, irr, dr2, iter, convt, npert, dvscfp, drhop, dbecsum)
      ELSE
         CALL write_rec(where_rec, irr, dr2, iter, convt, npert, dvscfp, drhop)
      ENDIF
      !
      IF (check_stop_now()) CALL stop_smoothly_ph(.FALSE.)
      !
      IF (convt) EXIT
      !
   ENDDO ! iteration
   !
   ! Update the responses with the Fermi energy shift
   !
   IF (lmetq0) THEN
      !
      ! Update the response potential.
      !
      DO ipert = 1, npert
         dvscfp(:, 1, ipert) = dvscfp(:, 1, ipert) - def(ipert)
         IF (doublegrid) dvscfs(:, 1, ipert) = dvscfs(:, 1, ipert) - def(ipert)
         !
         IF (lsda) THEN
            dvscfp(:, 2, ipert) = dvscfp(:, 2, ipert) - def(ipert)
            IF (doublegrid) dvscfs(:, 2, ipert) = dvscfs(:, 2, ipert) - def(ipert)
         ENDIF
      ENDDO
      !
      ! Update the response wavefunction (dpsi, stored in buffer iudwf).
      ! Also update drhos. (drhop is updated by ef_shift.)
      !
      CALL ef_shift_wfc(npert, ldoss, drhos)
      !
   ENDIF ! lmetq0
   !
   CALL apply_dpot_deallocate()
   !
   DEALLOCATE(dvscftmp)
   IF (noncolin) DEALLOCATE(dbecsum_nc)
   IF (noncolin .AND. domag .AND. okvan) THEN
      DEALLOCATE(int3_nc_save)
      DEALLOCATE(dbecsum_aux)
   ENDIF
   DEALLOCATE(mixin)
   DEALLOCATE(mixout)
   !
   IF (lmetq0) THEN
      DEALLOCATE(ldoss)
      DEALLOCATE(ldos)
      IF (okpaw) DEALLOCATE(becsum1)
   ENDIF
   !
   IF (twochem) CALL deallocate_twochem()
   !
   CALL stop_clock ('dfpt_kernel')
   !
!------------------------------------------------------------------------------
END SUBROUTINE dfpt_kernel
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
SUBROUTINE check_all_convt( convt )
   !---------------------------------------------------------------------------
   !! Work out how many processes in the image have converged.
   !! If only some processors have converged, call errore.
   !---------------------------------------------------------------------------
   !
   USE mp,        ONLY : mp_sum
   USE mp_images, ONLY : nproc_image, intra_image_comm
   !
   IMPLICIT NONE
   !
   LOGICAL, INTENT(IN) :: convt
   INTEGER :: tot_conv
   !
   IF (nproc_image==1) RETURN
   !
   ! Work out how many processes have converged
   !
   tot_conv = 0
   IF (convt) tot_conv = 1
   CALL mp_sum(tot_conv, intra_image_comm)
   !
   IF ((tot_conv > 0) .AND. (tot_conv < nproc_image)) THEN
     CALL errore('check_all_convt', 'Only some processors converged: '&
                &' either something is wrong with dfpt_kernel, or a different'&
                &' parallelism scheme should be used.', 1)
   ENDIF
   !
!------------------------------------------------------------------------------
END SUBROUTINE check_all_convt
!------------------------------------------------------------------------------
!
END MODULE dfpt_kernels
!------------------------------------------------------------------------------
