!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE  lr_two_chem
  USE kinds, ONLY : DP
  USE dfpt_type, ONLY : dfpt_ldos_type
  COMPLEX(DP),SAVE,PUBLIC :: def_val(3)
  COMPLEX(DP),SAVE,PUBLIC :: def_cond(3)
  ! the change of the Fermi energy for each pert. for valence and conduction manifold in the twochem case.
  COMPLEX(DP),ALLOCATABLE, SAVE,PUBLIC :: drhos_cond(:,:,:)
  !! output: the change of the scf charge
  COMPLEX(DP),ALLOCATABLE, SAVE,PUBLIC :: drhop_cond(:,:,:)
  COMPLEX(DP),ALLOCATABLE, SAVE,PUBLIC :: dbecsum_cond(:,:,:,:),dbecsum_cond_nc(:,:,:,:,:,:)
  TYPE(dfpt_ldos_type), SAVE, PUBLIC :: ldos_cond_data
  !! Local density of states of the conduction band states at the conduction band Fermi level
  !! Contains: dos_ef, ldos, ldoss, becsum_dos
  CONTAINS
!
!-----------------------------------------------------------------------
subroutine ef_shift_twochem (npert, ldos_data, drhop, drhop_cond, dbecsum, dbecsum_cond)
  !-----------------------------------------------------------------------
  !! This routine takes care of the effects of a shift of the two chemical potentials, due to the
  !! perturbation, that can take place in a metal at q=0, in the twochem case
  !! Optionally, update dbecsum using becsum_dos from ldos_data.
  !
  USE kinds,                ONLY : DP
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : gg
  USE noncollin_module,     ONLY : nspin_mag, nspin_lsda
  USE dfpt_type,            ONLY : dfpt_ldos_type
  !
  IMPLICIT NONE
  !
  ! input/output variables
  !
  INTEGER, INTENT(IN) :: npert
  !! the number of perturbation
  TYPE(dfpt_ldos_type), INTENT(IN) :: ldos_data
  !! Local density of states at Ef for valence states
  !! Contains: dos_ef, ldos, ldoss, becsum_dos
  !! Note: conduction DOS data comes from module variable ldos_cond_data
  COMPLEX(DP), INTENT(INOUT) :: drhop(dfftp%nnr, nspin_mag, npert)
  COMPLEX(DP), INTENT(INOUT) :: drhop_cond(dfftp%nnr, nspin_mag, npert)
  !! the change of the charge (with augmentation)
  COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum(:, :, :, :)
  !! input:  dbecsum = 2 <psi|beta> <beta|dpsi>
  !! output: dbecsum = 2 <psi|beta> <beta|dpsi> + def * ldos_data%becsum_dos
  COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum_cond(:, :, :, :)
  !! same as above, but restricted to conduction states
  !
  ! local variables
  !
  INTEGER :: is
  !! counter on spin polarizations
  INTEGER :: ipert
  !! counter on perturbations
  COMPLEX(DP) :: drhop_val(dfftp%nnr, nspin_mag,npert)
  COMPLEX(DP) :: delta_nv,delta_nc
  !! the change in electron number
  !! This may be complex since perturbation may be complex
  !
  REAL(DP), external :: w0gauss
  !! the smeared delta function
  !
  call start_clock ('ef_shift_twochem')
  !
  ! This routine is used only at q=Gamma where the dimension of irrep never exceeds 3
  IF (npert > 3) CALL errore("ef_shift_twochem", "npert exceeds 3", 1)
  !
  ! determines Fermi energy shift (such that each pertubation is neutral)
  !
  WRITE( stdout, * )
  drhop_val(:,:,:) = drhop(:,:,:)-drhop_cond(:,:,:)
  do ipert = 1, npert
     delta_nv = (0.d0, 0.d0)
     delta_nc = (0.d0, 0.d0)
     do is = 1, nspin_lsda
        CALL fwfft ('Rho', drhop_val(:,is,ipert), dfftp)
        if (gg(1) < 1.0d-8) delta_nv = delta_nv + omega*drhop_val(dfftp%nl(1),is,ipert)
        CALL invfft ('Rho', drhop_val(:,is,ipert), dfftp)
        !valence states
        CALL fwfft ('Rho', drhop_cond(:,is,ipert), dfftp)
        if (gg(1) < 1.0d-8) delta_nc = delta_nc + omega*drhop_cond(dfftp%nl(1),is,ipert)
        CALL invfft ('Rho', drhop_cond(:,is,ipert), dfftp)
     enddo
     call mp_sum ( delta_nv, intra_bgrp_comm )
     call mp_sum ( delta_nc, intra_bgrp_comm )
     IF ( ABS(ldos_data%dos_ef + ldos_cond_data%dos_ef) > 1.d-18 ) THEN
        def_val (ipert) = - delta_nv  / ldos_data%dos_ef
        def_cond(ipert) = - delta_nc / ldos_cond_data%dos_ef
     ELSE
        def_val (ipert) = 0.0_dp
        def_cond(ipert) = 0.0_dp
     ENDIF
  enddo
  !
  ! symmetrizes the Fermi energy shift
  !
  CALL sym_def(npert, def_val)
  WRITE( stdout, '(5x,"Pert. #",i3,": Fermi energy shift valence (Ry) =",2es15.4)')&
       (ipert, def_val (ipert) , ipert = 1, npert )
  CALL sym_def(npert, def_cond)
  WRITE( stdout, '(5x,"Pert. #",i3,": Fermi energy shift conduction (Ry) =",2es15.4)')&
       (ipert, def_cond (ipert) , ipert = 1, npert )
  !
  ! corrects the density response accordingly...
  !
  do ipert = 1, npert
     call zaxpy (dfftp%nnr*nspin_mag, def_val(ipert),  ldos_data%ldos,      1, drhop(1,1,ipert), 1)
     call zaxpy (dfftp%nnr*nspin_mag, def_cond(ipert), ldos_cond_data%ldos, 1, drhop(1,1,ipert), 1)
  enddo
  !
  ! In the PAW case there is also a metallic term
  !
  IF (PRESENT(dbecsum) .AND. ALLOCATED(ldos_data%becsum_dos)) THEN
     DO ipert = 1, npert
        dbecsum(:,:,:,ipert) = dbecsum(:,:,:,ipert) &
           + def_val(ipert)  * CMPLX(ldos_data%becsum_dos(:,:,:), 0.0_DP, KIND=DP) &
           + def_cond(ipert) * CMPLX(ldos_cond_data%becsum_dos(:,:,:), 0.0_DP, KIND=DP)
     ENDDO

  ENDIF
  !
  CALL stop_clock ('ef_shift_twochem')
  !
  end subroutine ef_shift_twochem
!-------------------------------------------------------------------------
!
!-------------------------------------------------------------------------
subroutine ef_shift_wfc_twochem(npert, ldos_data, drhos)
  !-----------------------------------------------------------------------
  !! This routine takes care of the effects of a shift of both chemical potentials, due to the
  !! perturbation, that can take place in a metal at q=0, on the wavefunctions.
  !
  USE kinds,                ONLY : DP
  USE mp,                   ONLY : mp_sum
  USE wavefunctions,        ONLY : evc
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE buffers,              ONLY : get_buffer, save_buffer
  USE wvfct,                ONLY : npwx, et, nbnd, nbnd_cond
  USE klist,                ONLY : degauss, ngauss, ngk, ltetra,degauss_cond
  USE ener,                 ONLY : ef,ef_cond
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag
  USE qpoint,               ONLY : nksq
  USE control_lr,           ONLY : nbnd_occ
  USE units_lr,             ONLY : iuwfc, lrwfc, lrdwf, iudwf
  USE eqv,                  ONLY : dpsi
  USE dfpt_tetra_mod,       ONLY : dfpt_tetra_delta
  USE dfpt_type,            ONLY : dfpt_ldos_type
  !
  IMPLICIT NONE
  !
  ! input/output variables
  !
  INTEGER, INTENT(IN) :: npert
  !! the number of perturbation
  TYPE(dfpt_ldos_type), INTENT(IN) :: ldos_data
  !! Local density of states at Ef for valence states
  !! Contains: dos_ef, ldos, ldoss, becsum_dos
  !! Note: conduction DOS data comes from module variable ldos_cond_data
  COMPLEX(DP), INTENT(INOUT) :: drhos(dffts%nnr, nspin_mag, npert)
  !! the change of the charge (with augmentation)
  !
  ! local variables
  !
  INTEGER :: npw, ibnd, ik, is, ipert, nrec, ikrec
  ! counter on occupied bands
  ! counter on k-point
  ! counter on spin polarizations
  ! counter on perturbations
  ! record number
  ! record position of wfc at k
  ! auxiliary for spin
  COMPLEX(DP) :: wfshift
  !! the shift coefficient for the wavefunction
  !! This may be complex since perturbation may be complex
  !
  REAL(DP), external :: w0gauss
  ! the smeared delta function
  !
  call start_clock ('ef_shift_wfc_twochem')
  !
  ! This routine is used only at q=Gamma where the dimension of irrep never exceeds 3
  IF (npert > 3) CALL errore("ef_shift_wfc_twochem", "npert exceeds 3", 1)
  !
  ! Update the perturbed wavefunctions according to the Fermi energy shift
  !
  do ik = 1, nksq
     npw = ngk (ik)
     !
     ! reads unperturbed wavefuctions psi_k in G_space, for all bands
     !
     ikrec = ik
     if (nksq > 1) call get_buffer (evc, lrwfc, iuwfc, ikrec)
     !
     ! reads delta_psi from iunit iudwf, k=kpoint
     !
     do ipert = 1, npert
        nrec = (ipert - 1) * nksq + ik
        IF (nksq > 1 .OR. npert > 1) CALL get_buffer(dpsi, lrdwf, iudwf, nrec)
        do ibnd = 1, nbnd_occ (ik)
                        if (ibnd.le.(nbnd-nbnd_cond)) then
                 wfshift = 0.5d0 * def_val(ipert) * &
                 w0gauss( (ef-et(ibnd,ik))/degauss, ngauss) / degauss
                 else
                 wfshift = 0.5d0 * def_cond(ipert) * &
                 w0gauss((ef_cond-et(ibnd,ik))/degauss_cond, ngauss) / degauss_cond
           end if
           !
           IF (noncolin) THEN
              call zaxpy (npwx*npol,wfshift,evc(1,ibnd),1,dpsi(1,ibnd),1)
           ELSE
              call zaxpy (npw, wfshift, evc(1,ibnd), 1, dpsi(1,ibnd), 1)
           ENDIF
        enddo
        !
        ! writes corrected delta_psi to iunit iudwf, k=kpoint,
        !
        IF (nksq > 1 .OR. npert > 1) CALL save_buffer(dpsi, lrdwf, iudwf, nrec)
     enddo
  enddo
  !
  do ipert = 1, npert
     do is = 1, nspin_mag
        call zaxpy (dffts%nnr, def_val(ipert),  ldos_data%ldoss(1,is),      1, drhos(1,is,ipert), 1)
        call zaxpy (dffts%nnr, def_cond(ipert), ldos_cond_data%ldoss(1,is), 1, drhos(1,is,ipert), 1)
     enddo
  enddo
  !
  CALL stop_clock ('ef_shift_wfc_twochem')
  !
  end subroutine ef_shift_wfc_twochem
  !
!
SUBROUTINE sternheimer_kernel_twochem(first_iter, time_reversed, npert, lrdvpsi, iudvpsi, &
         thresh, dvscfins, all_conv, avg_iter, drhoout, dbecsum, dbecsum_nc, &
         drhoout_cond,dbecsum_cond,dbecsum_cond_nc,exclude_hubbard)
   !----------------------------------------------------------------------------
   !This is a copy of the sternheimer kernel to be used when twochem and lmetq0.
   !In addition to the usual density response, it also calculates the conduction bands only
   !density response, needed to determine the conduction and valence manifold chemical potential
   !shift that can be present in the q=0 case
   !----------------------------------------------------------------------------
   USE kinds,                 ONLY : DP
   USE io_global,             ONLY : stdout
   USE mp,                    ONLY : mp_sum
   USE mp_pools,              ONLY : inter_pool_comm
   USE buffers,               ONLY : get_buffer, save_buffer
   USE fft_base,              ONLY : dffts
   USE ions_base,             ONLY : nat
   USE klist,                 ONLY : xk, wk, ngk, igk_k
   USE lsda_mod,              ONLY : lsda, nspin, current_spin, isk
   USE wvfct,                 ONLY : nbnd, npwx, et, nbnd_cond
   USE wavefunctions,         ONLY : evc
   USE noncollin_module,      ONLY : noncolin, domag, npol, nspin_mag
   USE uspp,                  ONLY : vkb
   USE uspp_param,            ONLY : nhm
   USE uspp_init,             ONLY : init_us_2
   USE ldaU,                  ONLY : lda_plus_u
   USE units_lr,              ONLY : iuwfc, lrwfc, lrdwf, iudwf
   USE control_lr,            ONLY : nbnd_occ, lgamma
   USE qpoint,                ONLY : nksq, ikks, ikqs
   USE qpoint_aux,            ONLY : ikmks, ikmkmqs, becpt
   USE eqv,                   ONLY : dpsi, dvpsi, evq
   USE apply_dpot_mod,        ONLY : apply_dpot_bands
   USE incdrhoscf_mod,        ONLY : incdrhoscf, incdrhoscf_nc
   !
   IMPLICIT NONE
   !
   LOGICAL, INTENT(IN) :: first_iter
   !! true if the first iteration.
   LOGICAL, INTENT(IN) :: time_reversed
   !! true if solving for time reversed wave functions
   LOGICAL, INTENT(IN), OPTIONAL :: exclude_hubbard
   !! true if ignoring the Hubbard response term
   INTEGER, INTENT(IN) :: npert
   !! number of perturbations
   INTEGER, INTENT(IN) :: lrdvpsi
   !! record length for the buffer storing dV_bare * psi
   INTEGER, INTENT(IN) :: iudvpsi
   !! unit for the buffer storing dV_bare * psi
   REAL(DP), INTENT(IN) :: thresh
   !! threshold for solving the linear equation
   LOGICAL, INTENT(OUT) :: all_conv
   !! True if converged at all k points and perturbations
   REAL(DP), INTENT(OUT) :: avg_iter
   !! average number of iterations for the linear equation solver
   COMPLEX(DP), INTENT(IN) :: dvscfins(dffts%nnr, nspin_mag, npert)
   !! dV_ind calculated in the previous iteration
   COMPLEX(DP), INTENT(INOUT) :: drhoout(dffts%nnr, nspin_mag, npert)
   !! induced charge density
   COMPLEX(DP), INTENT(INOUT) :: drhoout_cond(dffts%nnr, nspin_mag, npert)
   !! induced charge density
   COMPLEX(DP), INTENT(INOUT) :: dbecsum(nhm*(nhm+1)/2, nat, nspin_mag, npert)
   !! becsum with dpsi
   COMPLEX(DP), INTENT(INOUT) :: dbecsum_cond(nhm*(nhm+1)/2, nat, nspin_mag, npert)
   !! becsum with dpsi
   COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum_nc(nhm, nhm, nat, nspin, npert)
   !! becsum with dpsi. Used if noncolin is true.
   COMPLEX(DP), INTENT(INOUT), OPTIONAL :: dbecsum_cond_nc(nhm, nhm, nat, nspin, npert)
   !! becsum with dpsi. Used if noncolin is true.
   !
   LOGICAL :: conv_root
   !! true if linear system is converged
   LOGICAL :: exclude_hubbard_
   !! Local variable to set the default of exclude_hubbard to false
   INTEGER :: ikk, ikq, npw, npwq, ipert, num_iter, ik, nrec, ikmk, ikmkmq
   !! counters
   INTEGER :: tot_num_iter
   !! total number of iterations in cgsolve_all
   INTEGER :: tot_cg_calls
   !! total number of cgsolve_all calls
   REAL(DP) :: anorm
   !! the norm of the error of the conjugate gradient solution
   REAL(DP) :: rsign
   !! sign of the term in the magnetization
   REAL(DP), ALLOCATABLE :: h_diag(:, :)
   !! diagonal part of the Hamiltonian, used for preconditioning
   COMPLEX(DP) , ALLOCATABLE :: aux2(:, :)
   !! temporary storage used in apply_dpot_bands
   !
   EXTERNAL ch_psi_all, cg_psi
   !! functions passed to cgsolve_all
   !
   ! Initialization
   !
   CALL start_clock("sth_kernel")
   !
   exclude_hubbard_ = .FALSE.
   IF (PRESENT(exclude_hubbard)) exclude_hubbard_ = exclude_hubbard
   !
   ALLOCATE(h_diag(npwx*npol, nbnd))
   ALLOCATE(aux2(npwx*npol, nbnd))
   !
   !$acc enter data create(aux2(1:npwx*npol, 1:nbnd))
   !
   all_conv = .TRUE.
   tot_num_iter = 0
   tot_cg_calls = 0
   !
   DO ik = 1, nksq
      ikk  = ikks(ik)
      ikq  = ikqs(ik)
      npw  = ngk(ikk)
      npwq = ngk(ikq)
      !
      ! Set time-reversed k and k+q points
      !
      IF (time_reversed) THEN
         ikmk = ikmks(ik)
         ikmkmq = ikmkmqs(ik)
         rsign = -1.0_DP
      ELSE
         ikmk = ikk
         ikmkmq = ikq
         rsign = 1.0_DP
      ENDIF
      !
      IF (lsda) current_spin = isk(ikk)
      !
      ! reads unperturbed wavefunctions psi_k in G_space, for all bands
      ! if q=0, evq is a pointer to evc
      !
      IF (nksq > 1 .OR. (noncolin .AND. domag)) THEN
         IF (lgamma) THEN
            CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
         ELSE
            CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
            CALL get_buffer(evq, lrwfc, iuwfc, ikmkmq)
         ENDIF
      ENDIF
      !
      ! compute beta functions and kinetic energy for k-point ik
      ! needed by h_psi, called by ch_psi_all, called by cgsolve_all
      !
      CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb, .true.)
      !$acc update host(vkb)
      CALL g2_kin(ikq)
      !
      ! compute preconditioning matrix h_diag used by cgsolve_all
      !
      CALL h_prec(ik, evq, h_diag)
      !
      DO ipert = 1, npert
         !
         ! read P_c^+ x psi_kpoint into dvpsi.
         !
         nrec = (ipert - 1) * nksq + ik
         IF (time_reversed) nrec = nrec + npert * nksq
         !
         CALL get_buffer(dvpsi, lrdvpsi, iudvpsi, nrec)
         !
         IF (.NOT. first_iter) THEN
            !
            ! calculates dvscf_q*psi_k in G_space, for all bands, k=kpoint
            ! dvscf_q from previous iteration (mix_potential)
            !
            CALL apply_dpot_bands(ik, nbnd_occ(ikk), dvscfins(:, :, ipert), evc, aux2)
            dvpsi = dvpsi + aux2
            !
            !  In the case of US pseudopotentials there is an additional
            !  selfconsist term which comes from the dependence of D on
            !  V_{eff} on the bare change of the potential
            !
            CALL adddvscf(ipert, ik, time_reversed)
            !
            ! DFPT+U: add to dvpsi the scf part of the response
            ! Hubbard potential dV_hub
            !
            IF (lda_plus_u .AND. (.NOT. exclude_hubbard_)) CALL adddvhubscf(ipert, ik)
            !
         ENDIF
         !
         ! Orthogonalize dvpsi to valence states
         !
         CALL orthogonalize(dvpsi, evq, ikmk, ikmkmq, dpsi, npwq, .FALSE.)
         !
         ! Initial guess for dpsi
         !
         IF (first_iter) THEN
            !
            !  At the first iteration dpsi is set to zero
            !
            dpsi(:, :) = (0.d0,0.d0)
         ELSE
            !
            ! starting value for delta_psi is read from iudwf
            !
            CALL get_buffer(dpsi, lrdwf, iudwf, nrec)
         ENDIF
         !
         ! iterative solution of the linear system (H-e)*dpsi=dvpsi
         ! dvpsi=-P_c+ (dvbare+dvscf)*psi , dvscf fixed.
         !
         conv_root = .TRUE.
         !
         ! TODO: should nbnd_occ(ikk) be nbnd_occ(ikmk)?
         CALL cgsolve_all(ch_psi_all, cg_psi, et(1, ikmk), dvpsi, dpsi, h_diag, &
            npwx, npwq, thresh, ik, num_iter, conv_root, anorm, nbnd_occ(ikk), npol)
         !
         tot_num_iter = tot_num_iter + num_iter
         tot_cg_calls = tot_cg_calls + 1
         !
         IF (.NOT. conv_root) THEN
            all_conv = .FALSE.
            WRITE( stdout, "(5x, 'kpoint', i4, ' sternheimer_kernel: &
               &root not converged, thresh < ', es10.3)") ik, anorm
         ENDIF
         !
         ! writes delta_psi on iunit iudwf, k=kpoint,
         !
         CALL save_buffer(dpsi, lrdwf, iudwf, nrec)
         !
         ! calculates dvscf, sum over k => dvscf_q_ipert
         !
         IF (noncolin) THEN
            CALL incdrhoscf_nc(drhoout(1,1,ipert), wk(ikk), ik, &
                               dbecsum_nc(1,1,1,1,ipert), dpsi, rsign)
            CALL incdrhoscf_nc(drhoout_cond(1,1,ipert), wk(ikk), ik, &
                 dbecsum_cond_nc(1,1,1,1,ipert), dpsi, rsign, &
                 firstband=1+nbnd-nbnd_cond )
         ELSE
            CALL incdrhoscf(drhoout(1,current_spin,ipert), wk(ikk), &
                            ik, dbecsum(1,1,current_spin,ipert), dpsi)
            CALL incdrhoscf(drhoout_cond(1,current_spin,ipert), wk(ikk), &
                 ik, dbecsum_cond(1,1,current_spin,ipert), dpsi, &
                            firstband=1+nbnd-nbnd_cond )
         ENDIF
      ENDDO ! ipert
   ENDDO ! ik
   !
   CALL mp_sum(tot_num_iter, inter_pool_comm)
   CALL mp_sum(tot_cg_calls, inter_pool_comm)
   avg_iter = REAL(tot_num_iter, DP) / REAL(tot_cg_calls, DP)
   !
   !$acc exit data delete(aux2)
   !
   DEALLOCATE(aux2)
   DEALLOCATE(h_diag)
   !
   CALL stop_clock("sth_kernel")
   !
!----------------------------------------------------------------------------
END SUBROUTINE sternheimer_kernel_twochem
!
SUBROUTINE allocate_twochem(npe, nsolv)
   !
   ! allocate arrays for twochem calculation
   ! NOTE: also computes local dos
   !
   USE fft_base,             ONLY : dfftp, dffts
   USE ions_base,            ONLY : nat
   USE uspp_param,           ONLY : nhm
   USE lsda_mod,             ONLY : nspin
   USE paw_variables,        ONLY : okpaw
   USE noncollin_module,     ONLY : noncolin, nspin_mag
   USE dfpt_type,            ONLY : allocate_dfpt_ldos
   USE localdos_mod,         ONLY : localdos_new
   USE ener,                 ONLY : ef_cond
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: npe, nsolv
   !
   IF (noncolin) allocate (dbecsum_cond_nc (nhm,nhm, nat , nspin , npe, nsolv))
   allocate (drhos_cond ( dffts%nnr, nspin_mag , npe))
   allocate (drhop_cond ( dfftp%nnr, nspin_mag , npe))
   allocate (dbecsum_cond ( (nhm * (nhm + 1))/2 , nat , nspin_mag , npe))
   CALL allocate_dfpt_ldos(ldos_cond_data)
   !
   CALL localdos_new(ldos_cond_data, ef_cond)
   !
END SUBROUTINE allocate_twochem
!
SUBROUTINE deallocate_twochem
   USE noncollin_module,     ONLY : noncolin
   USE dfpt_type,            ONLY : deallocate_dfpt_ldos
   !
   IMPLICIT NONE
   !
   !deallocate for twochem calculation at gamma
   CALL deallocate_dfpt_ldos(ldos_cond_data)
   deallocate (dbecsum_cond)
   IF (noncolin) deallocate (dbecsum_cond_nc)
   deallocate (drhop_cond)
   deallocate (drhos_cond)
END SUBROUTINE deallocate_twochem
!
END MODULE lr_two_chem
