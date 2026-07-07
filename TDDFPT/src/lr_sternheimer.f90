!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE lr_sternheimer
  !
  !    This routine generalizes to finite complex frequencies and
  !    finite q vectors the routine solve_e of the Quantum ESPRESSO
  !    distribution.
  !
  !    This routine is a driver for the solution of the linear system which
  !    defines the change of the wavefunction due to an electric field
  !    of finite wavevector q and complex frequency omega.
  !    It performs the following tasks:
  !     a) computes the bare potential term  e^{iqr} | psi >
  !     b) adds to it the screening term Delta V_{SCF} | psi >
  !     c) applies P_c^+ (orthogonalization to valence states)
  !     d) calls cgsolve_all to solve the linear system at zero
  !        frequency or ccg_many_vectors
  !     e) computes Delta rho, Delta V_{SCF} and symmetrizes them
  !
CONTAINS

SUBROUTINE one_sternheimer_step(iu, flag)
    !
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : e2, fpi, rytoev
    USE ions_base,              ONLY : nat
    USE io_global,              ONLY : stdout
    USE io_files,               ONLY : diropn, nwordwfc
    USE cell_base,              ONLY : tpiba2
    USE fft_interfaces,         ONLY : fwfft
    USE klist,                  ONLY : lgauss, xk
    USE fft_base,               ONLY : dfftp
    USE lsda_mod,               ONLY : lsda, current_spin, isk
    USE klist,                  ONLY : ngk, igk_k
    USE buffers,                ONLY : get_buffer, save_buffer
    USE wavefunctions,          ONLY : evc
    USE uspp,                   ONLY : vkb
    USE uspp_param,             ONLY : nhm
    USE noncollin_module,       ONLY : nspin_mag
    USE gvect,                  ONLY : gg
    USE paw_variables,          ONLY : okpaw
    USE eqv,                    ONLY : dvpsi
    USE units_lr,               ONLY : lrwfc, iuwfc
    USE control_lr,             ONLY : alpha_pv, convt, nmix_ph, rec_code_read, niter_ph
    USE qpoint,                 ONLY : xq, nksq, ikks, ikqs
    USE mp_bands,               ONLY : intra_bgrp_comm
    USE mp,                     ONLY : mp_sum
    USE lr_variables,           ONLY : fru, fiu, iundvpsi, n_ipol, lr_verbosity, &
                                       chirr, chirz, chizr, chizz, epsm1, &
                                       current_w, itermax
    USE wavefunctions,          ONLY : psic
    USE dfpt_type,              ONLY : dfpt_data_type, allocate_dfpt_data, &
                                       deallocate_dfpt_data
    USE uspp_init,              ONLY : init_us_2
    USE lr_sym_mod,             ONLY : psymeq
    USE dfpt_kernels,           ONLY : dfpt_kernel
    USE lr_symm_base,           ONLY : nsymq, minus_q, lr_npert, upert
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: iu
    INTEGER, INTENT(IN) :: flag   ! if 1 compute the charge-charge and
                                  ! charge magnetization responses
                                  ! if 2 and lsda computes the magnetization
                                  ! magnetization response
    REAL(DP) ::  dr2
    ! dr2   : self-consistency error
    !
    COMPLEX(DP) , ALLOCATABLE ::   &
                   drhoscfout (:,:), & ! change of the scf charge (output)
                   mixin(:), mixout(:)  ! auxiliary for paw mixing
    !
    INTEGER :: iter0, ik, ikk, ikq, is, nrec, ndim, npw, npwq
    ! counters
    INTEGER :: npert
    !! number of perturbations
    INTEGER :: ipert
    !! Counter on perturbations
    INTEGER :: ndim_pot, ndim_paw
    !!
    REAL(DP) :: xqmod2, alpha_pv0
    !
    COMPLEX(DP) :: w  !frequency
    LOGICAL :: ldpsi1
    TYPE(dfpt_data_type) :: dfpt_data
    !! Data that describes linear response quantities
    !
    !
    CALL start_clock ('stern_step')
    !
    ! NOTE: In EELS, n_ipol = 1. In this code, n_ipol = 1 is assumed in a few places.
    !
    npert = n_ipol
    !
    nmix_ph = 20  ! TODO: Add this to input parameter, rename to nmix_dfpt ?
    niter_ph = itermax
    !
    ! Setup symmetry representation. We assume exp(i*q*r) is the perturbation, so that the
    ! symmetry representation is a trivial identity.
    ! upert_mq is not used since for finite-frequency DFPT, the time-reversal
    ! symmetry cannot be used.
    !
    lr_npert = npert
    ALLOCATE(upert(lr_npert, lr_npert, nsymq))
    upert(1, 1, :) = (1.d0, 0.d0)
    !
    w=CMPLX(fru(iu),fiu(iu))
    ldpsi1=ABS(w)>1.D-7
    alpha_pv0=alpha_pv
    alpha_pv=alpha_pv0 + REAL(w)
    !
    ALLOCATE (drhoscfout(dfftp%nnr, nspin_mag))
    !
    ndim_pot = dfftp%nnr * nspin_mag * 1
    IF (okpaw) THEN
       ndim_paw = (nhm * (nhm+1) * nat * nspin_mag * 1) / 2
       ALLOCATE(mixin(ndim_pot + ndim_paw))
       ALLOCATE(mixout(ndim_pot + ndim_paw))
       mixin = (0.0_DP, 0.0_DP)
    ELSE
       ALLOCATE(mixin(1))
       ALLOCATE(mixout(1))
    ENDIF
    !
    CALL allocate_dfpt_data(dfpt_data, 1)
    !
    !$acc enter data create(dfpt_data%dvscfs)
    dvpsi =(0.0d0, 0.0d0)

!    IF (rec_code_read == -20.AND.ext_recover) then
!       ! restarting in Electric field calculation
!       IF (okpaw) THEN
!          CALL read_rec(dr2, iter0, 1, dfpt_data%dvscfp, dfpt_data%dvscfs, dfpt_data%drhop, dfpt_data%dbecsum)
!          CALL setmixout(3*dfftp%nnr*nspin_mag,(nhm*(nhm+1)*nat*nspin_mag*3)/2, &
!                      mixin, dfpt_data%dvscfp, dfpt_data%dbecsum, ndim, -1 )
!       ELSE
!          CALL read_rec(dr2, iter0, 1, dfpt_data%dvscfp, dfpt_data%dvscfs)
!       ENDIF
!    ELSEIF (rec_code_read > -20 .AND. rec_code_read <= -10) then
!       ! restarting in Raman: proceed
!       convt = .true.
!    ELSE
       convt = .false.
       iter0 = 0
!    ENDIF
    !
!    IF ( ionode .AND. fildrho /= ' ') THEN
!       INQUIRE (UNIT = iudrho, OPENED = exst)
!       IF (exst) CLOSE (UNIT = iudrho, STATUS='keep')
!       CALL diropn (iudrho, TRIM(fildrho)//'.E', lrdrho, exst)
!    ENDIF
    IF (rec_code_read > -20) convt=.TRUE.
    !
    IF (convt) go to 155
    !
    IF ((lgauss.and..not.ldpsi1)) &
            CALL errore ('solve_eq', 'insert a finite frequency', 1)
    !
    IF (lr_verbosity > 5) THEN
       WRITE(stdout,'("<lr_sternheimer_one_step>")')
    ENDIF
    !
    IF (.NOT. ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
    !
    ! Calculate bare perturbation multiplied to the wavefunctions, save on buffer iundvpsi
    !
    DO ik = 1, nksq
       !
       ikk  = ikks(ik)
       ikq  = ikqs(ik)
       npw  = ngk(ikk)
       npwq = ngk(ikq)
       IF (lsda) current_spin = isk (ikk)
       !
       ! Read unperturbed wavefuctions evc (wfct at k)
       ! and evq (wfct at k+q)
       !
       IF (nksq > 1) THEN
          CALL get_buffer(evc, lrwfc, iuwfc, ikk)
       ENDIF
       !
       ! Calculate beta-functions vkb at k+q (Kleinman-Bylander projectors)
       ! The vkb's are needed for the non-local potential in h_psi,
       ! and for the ultrasoft term.
       !
       CALL init_us_2(npwq, igk_k(1, ikq), xk(1, ikq), vkb, .true.)
       !$acc update host(vkb)
       !
       ! do over polarization
       !
       DO ipert = 1, npert
          !
          nrec = (ipert - 1) * nksq + ik
          !  IF (isolv == 2) nrec = nrec + npert * nksq
          !
          CALL dveqpsi_us(ik)
          !
          !  with flag=2 the perturbation is a magnetic field along z
          !
          IF (lsda.AND.current_spin==2.AND.flag==2) dvpsi=-dvpsi
          !
          CALL save_buffer(dvpsi, nwordwfc, iundvpsi, nrec)
          !
       ENDDO ! ipert
    ENDDO ! ik
    !
    !    Solve DFPT self-consistent equation
    !
    CALL dfpt_kernel('turboEELS', npert, iter0, nwordwfc, iundvpsi, dr2, dfpt_data, 0, 0, w_freq = w)
    !
155 CONTINUE
    !
    drhoscfout(:,:) = dfpt_data%drhop(:,:,1)
    !
    !  compute here the susceptibility and the inverse of the dielectric
    !  constant
    !
    !  CALL compute_susceptibility(drhoscfout)
    !
    DO is=1,nspin_mag
       CALL fwfft ('Rho', drhoscfout(:,is), dfftp)
    ENDDO
    !
    IF (flag==1) THEN
       chirr(iu)=(0.0_DP,0.0_DP)
       chizr(iu)=(0.0_DP,0.0_DP)
       epsm1(iu)=(0.0_DP,0.0_DP)
    ELSE
       chirz(iu)=(0.0_DP,0.0_DP)
       chizz(iu)=(0.0_DP,0.0_DP)
    ENDIF
    !
    xqmod2=(xq(1)**2+xq(2)**2+xq(3)**2)*tpiba2
    !
    IF (ABS(gg(1))<1.d-8) THEN
       IF (flag==1) THEN
          chirr(iu) = drhoscfout(dfftp%nl(1),1)
          IF (lsda) chirr(iu) = chirr(iu) + drhoscfout(dfftp%nl(1),2)
          epsm1(iu) = CMPLX(1.0_DP,0.0_DP)+ chirr(iu)*fpi*e2/xqmod2
          IF (lsda) chizr(iu) = drhoscfout(dfftp%nl(1),1) - &
                                drhoscfout(dfftp%nl(1),2)
       ELSEIF (lsda) THEN
          chizz(iu)=drhoscfout(dfftp%nl(1),1)-drhoscfout(dfftp%nl(1),2)
          chirz(iu)=drhoscfout(dfftp%nl(1),1)+drhoscfout(dfftp%nl(1),2)
       ENDIF
    ENDIF
    !
    IF (flag==1) THEN
       CALL mp_sum(epsm1(iu),intra_bgrp_comm)
       CALL mp_sum(chirr(iu),intra_bgrp_comm)
       CALL mp_sum(chizr(iu),intra_bgrp_comm)
    ELSE
       CALL mp_sum(chizz(iu),intra_bgrp_comm)
       CALL mp_sum(chirz(iu),intra_bgrp_comm)
    ENDIF
    !
    IF (flag==1) THEN
       WRITE(stdout, '(/,6x,"Inverse dielectric constant at &
                          &frequency",f9.4," +",f9.4," i Ry")') fru(iu), fiu(iu)
       WRITE(stdout, '(46x,f9.4," +",f9.4," i eV")') current_w * rytoev
       WRITE(stdout,'(/,6x,"epsilon^-1(q,w) =",2f15.6)') epsm1(iu)
       !
       WRITE( stdout, '(/,6x,"Charge-charge susceptibility:")')
       !
       WRITE(stdout,'(/,6x,"chirr(q,w) =",2f15.6)') chirr(iu)
       IF (lsda) THEN
          WRITE(stdout,'(/,6x,"m_z-charge susceptibility:")')
          WRITE(stdout,'(/,6x,"chizr(q,w) =",2f15.6)') chizr(iu)
       ENDIF
       !
    ELSEIF (lsda) THEN
       WRITE( stdout, '(/,6x,"m_z - m_z susceptibility at &
                       &frequency",f9.4," +",f9.4," i Ry")') fru(iu), fiu(iu)
       WRITE( stdout, '(43x,f9.4," +",f9.4," i eV")') current_w * rytoev
       WRITE(stdout,'(/,6x,"chizz(q,w) =",2f15.6)') chizz(iu)
       WRITE(stdout,'(/,6x,"chirz(q,w) =",2f15.6)') chirz(iu)
    ENDIF
    !
    IF (okpaw) THEN
       DEALLOCATE(mixin)
       DEALLOCATE(mixout)
    ENDIF
    DEALLOCATE(upert)
    deallocate (drhoscfout)
    !$acc exit data delete(dfpt_data%dvscfs)
    CALL deallocate_dfpt_data(dfpt_data)
    !
    alpha_pv=alpha_pv0
    !
    CALL stop_clock ('stern_step')
    !
    RETURN
    !
END SUBROUTINE one_sternheimer_step

END MODULE lr_sternheimer
