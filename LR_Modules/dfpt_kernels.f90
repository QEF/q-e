!
! Copyright (C) 2025-2026 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------------------------
MODULE dfpt_kernels
IMPLICIT NONE
CONTAINS
SUBROUTINE dfpt_kernel(code, npert, iter0, lrdvpsi, iudvpsi, dr2, dfpt_data, &
                       irr, imode0, write_rec_callback, stop_callback, w_freq)
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
   !!
   !! Input/Output:
   !!    - drhos : change of the charge density (smooth part only, dffts)
   !!              (Note that this is different from drhop after fft_interpolate because
   !!               drhos does not include augmentation charges, while drhop does.)
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
   !! (If doublegrid == false and using NCPP, the two quantities are identical.)
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
   !------------------------------------------------------------------------------
   !
   USE kinds,                ONLY : DP
   USE mp_pools,             ONLY : inter_pool_comm
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_interfaces,       ONLY : fft_interpolate
   USE buffers,              ONLY : get_buffer
   USE ions_base,            ONLY : nat
   USE io_global,            ONLY : stdout
   USE io_files,             ONLY : diropn
   USE check_stop,           ONLY : check_stop_now
   USE klist,                ONLY : ltetra, lgauss
   USE gvecs,                ONLY : doublegrid
   USE lsda_mod,             ONLY : nspin, lsda
   USE scf,                  ONLY : rho
   USE ldaU,                 ONLY : lda_plus_u
   USE uspp,                 ONLY : okvan
   USE uspp_param,           ONLY : nhm
   USE noncollin_module,     ONLY : noncolin, domag, npol, nspin_mag
   USE paw_variables,        ONLY : okpaw
   USE paw_onecenter,        ONLY : paw_dpotential
   USE ener,                 ONLY : ef, ef_cond
   USE efermi_shift,         ONLY : ef_shift_new, ef_shift_wfc_new
   USE lrus,                 ONLY : int3_paw, int3_nc
   USE control_lr,           ONLY : lgamma, niter_ph, nmix_ph, tr2_ph, alpha_mix, convt, &
                                    lgamma_gamma, flmixdpot, where_rec, lnoloc
   USE dv_of_drho_lr,        ONLY : dv_of_drho
   USE lr_nc_mag,            ONLY : int3_nc_save
   USE apply_dpot_mod,       ONLY : apply_dpot_allocate, apply_dpot_deallocate
   USE dfpt_type,            ONLY : dfpt_data_type, dfpt_dvscfp_to_dvscfs, &
                                    dfpt_ldos_type, allocate_dfpt_ldos, deallocate_dfpt_ldos
   USE response_kernels,     ONLY : sternheimer_kernel, sternheimer_kernel_freq, &
                                    sternheimer_postprocess
   USE two_chem,             ONLY : twochem
   USE lr_restart,           ONLY : write_rec_interface, stop_callback_interface
   USE lr_two_chem,          ONLY : allocate_twochem, deallocate_twochem, &
                                    sternheimer_kernel_twochem, dbecsum_cond, &
                                    dbecsum_cond_nc, drhos_cond, ef_shift_wfc_twochem
   USE localdos_mod,         ONLY : localdos_new
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
   !! self-consistency error. Input is used for restart.
   TYPE(dfpt_data_type), INTENT(INOUT) :: dfpt_data
   !! Data that describes linear response quantities
   PROCEDURE(write_rec_interface), OPTIONAL :: write_rec_callback
   !! Callback subroutine for restart checkpointing
   PROCEDURE(stop_callback_interface), OPTIONAL :: stop_callback
   !! Callback subroutine for graceful stop
   COMPLEX(DP), INTENT(IN), OPTIONAL :: w_freq
   !! Frequency for finite-frequency DFPT.
   !
   INTEGER, INTENT(IN) :: irr
   !! input: the irreducible representation (to be removed after moving from PH to LR_Modules)
   !! Set to 1 if not using phonon perturbation.
   INTEGER, INTENT(IN) :: imode0
   !! input: the position of the modes (to be removed after moving from PH to LR_Modules)
   !! Set to 0 if not using phonon perturbation.
   !
   ! ... local variables
   !
   REAL(DP) :: thresh, averlt
   ! thresh: convergence threshold
   ! averlt: average number of iterations
   TYPE(dfpt_ldos_type) :: ldos_data
   !! Local density of states at Ef (encapsulates ldos, ldoss, becsum_dos, dos_ef)
   COMPLEX(DP), ALLOCATABLE :: dvscftmp(:, :, :)
   !! change of scf potential (output before mixing)
   COMPLEX(DP), ALLOCATABLE :: mixin(:), mixout(:)
   ! Misc work space
   COMPLEX(DP), ALLOCATABLE :: dbecsum_nc(:,:,:,:,:), dbecsum_nc_trev(:,:,:,:,:)
   ! dbecsum: the derivative of becsum
   !
   LOGICAL :: all_conv
   !! True if sternheimer_kernel is converged at all k points and perturbations
   logical :: lmetq0,     & ! true if xq=(0,0,0) in a metal
      first_iter    ! true if first iteration where induced rho is not yet calculated
   LOGICAL :: finite_freq
   !! True if solving a finite-frequency DFPT

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
   finite_freq = PRESENT(w_freq)
   !
   nsolv = 1
   IF ((noncolin .AND. domag) .OR. finite_freq) nsolv = 2
   !
   IF (twochem .AND. (noncolin .AND. domag)) CALL errore('dfpt_kernel', &
      'DFPT with twochem and noncollinear magnetism is not implemented', 1)
   IF (twochem .AND. finite_freq) CALL errore('dfpt_kernel', &
      'DFPT with twochem and finite frequency is not implemented', 1)
   IF ((noncolin .AND. domag) .AND. finite_freq) CALL errore('dfpt_kernel', &
      'DFPT with noncollinear magnetism and finite frequency is not implemented', 1)
   !
   CALL apply_dpot_allocate()
   !
   ALLOCATE(dvscftmp(dfftp%nnr, nspin_mag, npert))
   IF (noncolin) ALLOCATE(dbecsum_nc(nhm, nhm, nat, nspin, npert))
   IF (noncolin .AND. domag) ALLOCATE(dbecsum_nc_trev(nhm, nhm, nat, nspin, npert))
   IF (noncolin .AND. domag .AND. okvan) THEN
      ALLOCATE(int3_nc_save(nhm, nhm, nat, nspin_mag, npert, 2))
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
      CALL setmixout(ndim_pot, ndim_paw, mixin, dfpt_data%dvscfp, dfpt_data%dbecsum, ndim, -1)
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
      CALL allocate_dfpt_ldos(ldos_data)
      CALL localdos_new(ldos_data,ef)
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
      dfpt_data%drhos = (0.d0, 0.d0)
      dfpt_data%dbecsum = (0.d0, 0.d0)
      IF (noncolin) dbecsum_nc = (0.d0, 0.d0)
      IF (noncolin .AND. domag) dbecsum_nc_trev = (0.d0, 0.d0)
      !
      IF (lmetq0 .AND. twochem) THEN
         drhos_cond = (0.d0, 0.d0)
         dbecsum_cond = (0.d0, 0.d0)
         IF (noncolin) dbecsum_cond_nc = (0.d0, 0.d0)
      ENDIF
      !
      ! Set threshold for iterative solution of the linear system
      !
      IF (first_iter .OR. kter == 1) THEN
         ! If first iteration or first iteration after restart
         thresh = 1.0d-2
      ELSE
         thresh = min(1.d-1 * sqrt(dr2), thresh)
      ENDIF
      !
      ! DFPT+U: at each ph iteration calculate dnsscf,
      ! i.e. the scf variation of the occupation matrix ns.
      !
      IF (lda_plus_u .AND. (.NOT. first_iter)) CALL dnsq_scf(npert, lmetq0)
      !
      ! Start the loop on the two linear systems, one at B and one at -B
      !
      IF (finite_freq) THEN
         !
         ! Finite-frequency Sternheimer equation
         ! The isolv loop is taken inside sternheimer_kernel_freq
         !
         ! time_reversed (noncollinear magnetic) case with finite-frequency Sternheimer
         ! is not implemented. So the second argument is .FALSE.
         !
         CALL sternheimer_kernel_freq(first_iter, .FALSE., npert, lrdvpsi, iudvpsi, &
             thresh, w_freq, dfpt_data%dvscfs, all_conv, averlt, dfpt_data%drhos, &
             dfpt_data%dbecsum, dbecsum_nc)
      ELSE
         !
         ! Zero-frequency Sternheimer equation
         !
         DO isolv = 1, nsolv
            !
            ! Compute drhos, the charge density response to the total potential
            !
            IF (.NOT. (twochem .AND. lmetq0)) then
               !
               ! Sternheimer kernel for the normal case
               !
               IF (isolv == 1) THEN
                  CALL sternheimer_kernel(first_iter, isolv==2, npert, lrdvpsi, iudvpsi, &
                     thresh, dfpt_data%dvscfs, all_conv, averlt, dfpt_data%drhos, dfpt_data%dbecsum, &
                     dbecsum_nc)
               ELSE
                  CALL sternheimer_kernel(first_iter, isolv==2, npert, lrdvpsi, iudvpsi, &
                     thresh, dfpt_data%dvscfs, all_conv, averlt, dfpt_data%drhos, dfpt_data%dbecsum, &
                     dbecsum_nc_trev)
               ENDIF
               !
            ELSE
               !
               ! Sternheimer kernel for the twochem case
               !
               CALL sternheimer_kernel_twochem(first_iter, isolv==2, npert, lrdvpsi, iudvpsi, &
                  thresh, dfpt_data%dvscfs, all_conv, averlt, dfpt_data%drhos, dfpt_data%dbecsum, &
                  dbecsum_nc, drhos_cond, dbecsum_cond, dbecsum_cond_nc)
               !
               IF (isolv == 2) THEN
                  CALL errore("dfpt_kernel", "DFPT with twochem and noncollinear magnetism is not implemented", 1)
               ENDIF
               !
            ENDIF
            !
         ENDDO ! isolv
      ENDIF
      !
      IF (noncolin .AND. domag) THEN
         CALL sternheimer_postprocess(nsolv, dfpt_data, dbecsum_nc, dbecsum_nc_trev)
      ELSE
         CALL sternheimer_postprocess(nsolv, dfpt_data, dbecsum_nc)
      ENDIF
      !
      !    Now we compute for all perturbations the total charge and potential
      !
      ! If .NOT. okpaw, dbecsum is not used below. So, postprocess dbecsum only if okpaw.
      !
      IF (okpaw) THEN
         !
         ! The presence of c.c. in the formula gives a factor 2.0
         !
         dfpt_data%dbecsum = 2.0_DP * dfpt_data%dbecsum
         !
         ! Add Pulay correction to dbecsum
         IF (ALLOCATED(dfpt_data%dbecsum_pulay)) THEN
            dfpt_data%dbecsum = dfpt_data%dbecsum + dfpt_data%dbecsum_pulay
         ENDIF
      ENDIF
      !
      ! Metallic case and q=0: add a correction to the response charge density
      ! due to the shift of the Fermi energy (see Eq.(75) in Rev. Mod. Phys. 73, 515 (2001)).
      ! This term is added to the response charge density (in order to obtain correct
      ! response HXC potential) and to the response occupation matrices.
      !
      IF (lmetq0 .AND. (.NOT. twochem)) THEN
         CALL ef_shift_new(ldos_data, dfpt_data)
      ENDIF
      !
      ! Post processing for the case of two chemical potentials
      !
      IF (twochem) THEN
         CALL twochem_postproc_dfpt(npert, nsolv, imode0, lmetq0, ldos_data, dfpt_data)
      ENDIF
      !
      !   After the loop over the perturbations we have the linear change
      !   in the charge density for each mode of this representation.
      !   Here we symmetrize them ...
      !
      IF (.NOT. lgamma_gamma) THEN
         ! FIXME: Improve psymdvscf using idea from psymeq, which saves computation
         !        by each core computing only its own part.
         CALL psymdvscf(dfpt_data%drhos, dffts)
         CALL psymdvscf(dfpt_data%drhop, dfftp)
         IF (okpaw) CALL PAW_dsymmetrize(dfpt_data%dbecsum)
      ENDIF
      !
      !   compute the corresponding change in scf potential : drhop -> dvscftmp
      !
      IF (lnoloc) THEN
         ! No local field effect: set dvscf to 0
         dvscftmp(:, :, :) = (0.d0, 0.d0)
      ELSE
         ! Compute the response HXC potential
         DO ipert = 1, npert
            CALL zcopy(dfftp%nnr*nspin_mag, dfpt_data%drhop(1,1,ipert), 1, dvscftmp(1,1,ipert), 1)
            IF (ALLOCATED(dfpt_data%drhoc)) THEN
               ! Consider core charge perturbation if present
               CALL dv_of_drho(dvscftmp(1, 1, ipert), drhoc = dfpt_data%drhoc(:, ipert))
            ELSE
               CALL dv_of_drho(dvscftmp(1, 1, ipert))
            ENDIF
         ENDDO
      ENDIF
      !
      !   And we mix with the old potential: dvscftmp -> dvscfp
      !
      IF (okpaw) THEN
         !
         !  In this case we mix also dbecsum
         !
         CALL setmixout(ndim_pot, ndim_paw, mixout, dvscftmp, dfpt_data%dbecsum, ndim, -1)
         CALL mix_potential(2 * (ndim_pot + ndim), mixout, mixin, alpha_mix(kter), dr2, &
                            npert * tr2_ph / npol, iter, nmix_ph, flmixdpot, convt)
         CALL setmixout(ndim_pot, ndim_paw, mixin, dfpt_data%dvscfp, dfpt_data%dbecsum, ndim, 1)
      ELSE
         CALL mix_potential(2 * ndim_pot, dvscftmp, dfpt_data%dvscfp, alpha_mix(kter), dr2, &
                            npert * tr2_ph / npol, iter, nmix_ph, flmixdpot, convt)
      ENDIF
      !
      !   fft_interpolate or copy potential from hard to smooth grid: dvscfp -> dvscfs
      !
      CALL dfpt_dvscfp_to_dvscfs(dfpt_data)
      !
      ! Check that convergence have been reached on ALL processors in this image.
      ! If convergence is reached on only some processors, call errore.
      !
      CALL check_all_convt(convt)
      !
      !   calculate here the change of the D1-~D1 coefficients due to the phonon
      !   perturbation
      !
      IF (okvan) THEN
         IF (okpaw) CALL PAW_dpotential(dfpt_data%dbecsum, rho%bec, int3_paw, npert)
         !
         !     with the new change of the potential we compute the integrals
         !     of the change of potential and Q
         !
         CALL newdq(dfpt_data%dvscfp, npert)
         !
         !  In the noncollinear magnetic case computes the int3 coefficients with
         !  the opposite sign of the magnetic field. They are saved in int3_nc_save,
         !  that must have been ALLOCATEd by the calling routine
         !
         IF (noncolin .AND. domag) THEN
            int3_nc_save(:,:,:,:,:,1) = int3_nc(:,:,:,:,:)
            IF (okpaw) rho%bec(:,:,2:4) = -rho%bec(:,:,2:4)
            DO ipert = 1, npert
               dfpt_data%dvscfp(:, 2:4, ipert) = -dfpt_data%dvscfp(:, 2:4, ipert)
               IF (okpaw) dfpt_data%dbecsum(:, :, 2:4, ipert) = -dfpt_data%dbecsum(:, :, 2:4, ipert)
            ENDDO
            !
            !   if needed recompute the paw coeffients with the opposite sign of
            !   the magnetic field
            !
            IF (okpaw) CALL PAW_dpotential(dfpt_data%dbecsum, rho%bec, int3_paw, npert)
            !
            !   here compute the int3 integrals
            !
            CALL newdq(dfpt_data%dvscfp, npert)
            int3_nc_save(:,:,:,:,:,2) = int3_nc(:,:,:,:,:)
            !
            !  restore the correct sign of the magnetic field.
            !
            DO ipert = 1, npert
               dfpt_data%dvscfp(:, 2:4, ipert) = -dfpt_data%dvscfp(:, 2:4, ipert)
               IF (okpaw) dfpt_data%dbecsum(:, :, 2:4, ipert) = -dfpt_data%dbecsum(:, :, 2:4, ipert)
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
      !    Here we save the information for recovering the run from this point
      !
      IF (PRESENT(write_rec_callback)) THEN
         CALL write_rec_callback(where_rec, irr, dr2, iter, convt, dfpt_data)
      ENDIF
      !
      IF (check_stop_now() .AND. PRESENT(stop_callback)) CALL stop_callback(.FALSE.)
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
         dfpt_data%dvscfp(:, 1, ipert) = dfpt_data%dvscfp(:, 1, ipert) - dfpt_data%def(ipert)
         dfpt_data%dvscfs(:, 1, ipert) = dfpt_data%dvscfs(:, 1, ipert) - dfpt_data%def(ipert)
         !
         IF (lsda) THEN
            dfpt_data%dvscfp(:, 2, ipert) = dfpt_data%dvscfp(:, 2, ipert) - dfpt_data%def(ipert)
            dfpt_data%dvscfs(:, 2, ipert) = dfpt_data%dvscfs(:, 2, ipert) - dfpt_data%def(ipert)
         ENDIF
      ENDDO
      !
      ! Update the response wavefunction (dpsi, stored in buffer iudwf).
      ! Also update drhos. (drhop is updated by ef_shift.)
      !
      IF (.NOT. twochem) THEN
         CALL ef_shift_wfc_new(dfpt_data)
      ELSE
         CALL ef_shift_wfc_twochem(npert, ldos_data, dfpt_data%drhos)
      ENDIF
      !
   ENDIF ! lmetq0
   !
   CALL apply_dpot_deallocate()
   !
   DEALLOCATE(dvscftmp)
   IF (noncolin) DEALLOCATE(dbecsum_nc)
   IF (noncolin .AND. domag) DEALLOCATE(dbecsum_nc_trev)
   IF (noncolin .AND. domag .AND. okvan) THEN
      DEALLOCATE(int3_nc_save)
   ENDIF
   DEALLOCATE(mixin)
   DEALLOCATE(mixout)
   !
   IF (lmetq0) THEN
      CALL deallocate_dfpt_ldos(ldos_data)
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
