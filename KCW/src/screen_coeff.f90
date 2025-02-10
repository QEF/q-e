! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO CMPLX(0.D0,0.D0,kind=DP)
#define ONE (0.D0,1.D0)
!#define DEBUG
!-----------------------------------------------------------------------
SUBROUTINE screen_coeff ()
  !---------------------------------------------------------------------
  !
  !! LR calculation of the screening coefficients. ADAPETD from do_phonon
  !
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE klist,                ONLY : nkstot
  USE mp,                   ONLY : mp_sum
  USE control_kcw,          ONLY : kcw_iverbosity, spin_component, num_wann, iorb_start, l_do_alpha, &
                                   iorb_end, alpha_final, nqstot, eps_inf, l_vcut, l_unique_manifold, &  
                                   group_alpha, tmp_dir_kcw, iurho_wann, tmp_dir_kcwq, x_q, tmp_dir_save, &
                                   i_orb, nrho, wq_ibz, fbz2ibz, irr_bz, setup_pw
  USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux, nspin_mag, npol
  USE buffers,              ONLY : get_buffer, save_buffer
  USE io_global,            ONLY : stdout, ionode
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE cell_base,            ONLY : at
  USE dv_of_drho_lr,        ONLY : dv_of_drho
  USE fft_interfaces,       ONLY : fwfft
  USE control_lr,           ONLY : lgamma
  USE lsda_mod,             ONLY : nspin
  USE gvecs,                ONLY : ngms
  USE io_files,             ONLY : tmp_dir, prefix
  USE solve_linter_koop_mod 
  !USE exx_base,             ONLY : exx_divergence
  USE coulomb,              ONLY : exxdiv, exxdiv_eps
  USE control_kcw,          ONLY : nsym_old
  USE symm_base,            ONLY : nsym
  USE cell_base,            ONLY : omega
  !
  IMPLICIT NONE
  ! 
  INTEGER :: iq, nqs, spin_ref, is, ip, iq_ibz
  !! Counter for the k/q points in the BZ, 
  !! nqs: total number of q points 
  !! spin_ref: spin of the Wannier function
  !! ip: spin/magnetizations index
  !! iq_ibz: qindex q point in the IBZ
  !
  INTEGER :: iwann, jwann, lrrho, iq_start, iun_res
  !! Band counter, leght of the rho record, starting iq (if restart), iunit partial results
  !
  COMPLEX(DP) :: rhowann(dffts%nnr, num_wann,nrho), rhor(dffts%nnr,nrho) 
  !! The periodic part of the wannier orbital density
  !
  COMPLEX (DP):: delta_vr(dffts%nnr,nspin_mag), delta_vr_(dffts%nnr,nspin_mag)
  !! The perturbing potential w and w/o q+G=0 contribution
  !
  COMPLEX(DP) :: vki_u(num_wann), sh(num_wann), vki_r(num_wann)
  ! ki unrelaxed and relaxed potential and self-hartree
  !
  COMPLEX(DP), ALLOCATABLE  :: rhog(:,:), delta_vg(:,:), vh_rhog(:), drhog_scf(:,:), drhor_scf(:,:), delta_vg_(:,:)
  ! wanier density, perturbing potential, hartree potential, density variation (g and r), perturbing pot with G0=0 
  !
  LOGICAL :: do_band 
  !
  COMPLEX(DP) :: phase(dffts%nnr), wann_c(dffts%nnr,num_wann,nrho), rho_c(dffts%nnr,num_wann,nrho)
  !! The phase associated to the hift k+q-> k'
  !
  COMPLEX(DP) :: int_rho, int_wann, pi_q_unrelax, pi_q_unrelax_, sh_q, pi_q_relax, pi_q_relax_rs
  COMPLEX(DP) :: struct_fact
  LOGICAL :: do_real_space = .false. 
  LOGICAL :: exst
  !
  COMPLEX(DP) :: drho_zero
  !
  REAL(DP) :: xq_(3), weight(nkstot)
  REAL(DP) :: alpha
  REAL(DP) :: div, div_eps
  !
  nqs = nqstot
  !
  IF (nqs == 1) do_real_space = .TRUE. 
  IF (do_real_space) THEN 
     ALLOCATE ( drhor_scf(dffts%nnr,nspin_mag) ) 
     drhor_scf = ZERO
  ENDIF
  ! If only 1 K point do also integral in REAL space (Extra check for the Super cell calculation)
  !
#ifdef DEBUG
  write(*,'(/,"DEBUG: The list of G  vectors")')
  do ig = 1, 10
     xq_ = g(:,ig)
     CALL cryst_to_cart(1, xq_, at, -1)
     write(*,'("i = ", i3, 3x, "G(i) = ", 3x, 3f8.4, "  [Cryst]",3x, "|G(i)| = ", f12.6 )') & 
                       ig, (xq_(iq), iq=1,3),  sqrt(sum (xq_(:)*xq_(:)))
  enddo
#endif
  !
  ! INITIALIZATIONs
  div = exxdiv
  div_eps = exxdiv_eps
  wann_c = ZERO
  vki_u = ZERO
  vki_r = ZERO
  sh = ZERO
  pi_q_relax = ZERO
  pi_q_relax_rs = ZERO
  pi_q_unrelax = ZERO
  iq_start = 1
  !
  drho_zero = ZERO
  !
  WRITE(stdout,'(/)')
  WRITE( stdout, '(5X,"INFO: LR CALCULATION ...")')
  !
  iun_res = 987
  OPEN (iun_res, file = TRIM(tmp_dir_kcw)//TRIM(prefix)//'.LR_res.txt')
  CALL restart_screen (num_wann, iq_start, vki_r, vki_u, sh, do_real_space)
  !
  nsym_old = nsym
  !
  !IF(irr_bz) CALL read_qlist_ibz()
  DO iq = iq_start, nqs
      !! For each q in the mesh 
    !
    CALL kcw_prepare_q ( do_band, setup_pw, iq ) 
    IF (kcw_iverbosity .gt. -1 ) WRITE(stdout,'(8x, "INFO: prepare_q DONE",/)') 
    !! Prepare the q point calculation
    !
    lrrho=num_wann*dffts%nnr*nrho
    tmp_dir = tmp_dir_save  ! the periodic part are written on the original outdir 
    CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
    tmp_dir = tmp_dir_kcwq       ! go back to the q-specific directory
    ! ... Retrive the rho_wann_q(r) from buffer in REAL space
    !
    IF (kcw_iverbosity .gt. -1 ) WRITE(stdout,'(8X, "INFO: rhowan_q(r) RETRIEVED"/)') 
    !
    ! The NSCF can be run only once for each qpoint if we are not using symmeties
    ! If using symmetry this nneds to be done inside the wannier loop as each wannier 
    ! as a different number of symmetries and hence k points (see below)
    IF (setup_pw .AND. .NOT. irr_bz) THEN 
       CALL kcw_run_nscf(do_band)
       IF (kcw_iverbosity .gt. -1) WRITE(stdout,'(/,8X, "INFO: NSCF calculation DONE",/)')
    ENDIF
    ! 
    IF (.NOT. irr_bz) THEN
       CALL kcw_initialize_ph ( ) 
       IF (kcw_iverbosity .gt. -1 ) WRITE(stdout,'(8X, "INFO: kcw_q initialization DONE",/)')
    ENDIF
    !
    ALLOCATE ( rhog (ngms,nrho) , delta_vg(ngms,nspin_mag), vh_rhog(ngms), drhog_scf (ngms, nspin_mag), delta_vg_(ngms,nspin_mag) )
    !
    IF ( lgamma .AND. .NOT. l_unique_manifold) CALL check_density (rhowann) 
    !
    DO iwann = iorb_start, iorb_end  ! for each wannier, that is actually the perturbation
      !
      IF ( .NOT. l_do_alpha (iwann)) CYCLE
      !
      WRITE( stdout, '(/, 8x,"Start LR calculation for the wannier #",i3, 3x, &
              "iq =", i3, 3x, "spin =", i3, /)') iwann, iq, spin_component
      !
      ! initialize all quantities that in general depend on the wannier function
      ! (because of symmetry)
      !
      IF(irr_bz) THEN
        IF( fbz2ibz(iq, iwann) .eq. -1 ) THEN
                WRITE(stdout, '(8X, "SYM : iq =", i5, 3x, "NOT DONE for WF =", i5, 3x, &
                        "as bz2ibz is -1")') iq, iwann
          CYCLE
        END IF
        CALL reset_symmetry_op(iwann) 
        iq_ibz = fbz2ibz(iq, iwann)
        weight(iq) = wq_ibz(iq_ibz, iwann)
      ELSE 
         iq_ibz = iq
         weight(iq) = 1.D0/nqs  !*nspin ! No SYMM  ! CHECK nspin??? 
      END IF
      !
      WRITE(stdout, '(8X, "SYM : iwann= ", i5, " iq = ", i5, " weight = ", f12.8)') iwann, iq, weight(iq)
       !
       ! consider only symmetries respected by wf iwann
       ! in arrays defined in symm_base
       !
       ! now the run_nscf happens inside the loop over wannier functions because 
       ! symmetries depend on wannier functions. Only if we use symmetries, otherwise
       ! The NSCF is done outside the Wannier loop as the k and q points are not going to change
       !
       IF (setup_pw .AND. irr_bz) THEN 
          CALL kcw_run_nscf(do_band)
          IF (kcw_iverbosity .gt. -1) WRITE(stdout,'(/,8X, "INFO: NSCF calculation DONE",/)')
       ENDIF
       ! ... IF needed run a nscf calculation for k+q
       ! ... NB: since the k+q is always a p we already have, this could be 
       !         avoided in principle (see compute_map and rho_of_q). 
       !         For the moment it's easier to  keep it. 
       !
       ! Skip LR calculation if this orbital match with a one already computed (see group_orbital)
       !
       drhog_scf (:,:) = ZERO
       rhog(:,:)         = ZERO
       delta_vg(:,:)   = ZERO
       delta_vg_(:,:)  = ZERO
       vh_rhog(:)      = ZERO
       rhor(:,:)         = ZERO
       !
       rhor(:,:) = rhowann(:,iwann,:)
       !! ... The periodic part of the orbital desity in real space
       !
       !
       !get arrays for the phase e^{iGr} when summing over small group of q
       !
       IF (irr_bz) CALL kcw_initialize_ph ( ) 
       CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_) 
       !! ... The periodic part of the perturbation
       !
       IF (nspin==2 .OR. nspin==1) THEN
        ! This can/should be an IF over noncolin
        pi_q_unrelax  = sum (CONJG(rhog (:,1)) * delta_vg(:,spin_component)) *weight(iq)*omega
        pi_q_unrelax_ = sum (CONJG(rhog (:,1)) * delta_vg_(:,spin_component))*weight(iq)*omega
       ELSEIF (nspin==4) THEN
        pi_q_unrelax = ZERO
        pi_q_unrelax_ = ZERO
        DO ip=1,nspin_mag
          pi_q_unrelax  = pi_q_unrelax + sum (CONJG(rhog (:,ip)) * delta_vg(:,ip)) *weight(iq)*omega
          pi_q_unrelax_ = pi_q_unrelax_+ sum (CONJG(rhog (:,ip)) * delta_vg_(:,ip))*weight(iq)*omega
        END DO
       ENDIF
       sh_q  = 0.5D0 * sum (CONJG(rhog (:,1)) * vh_rhog(:)) *weight(iq)*omega

       !
       CALL mp_sum (pi_q_unrelax,  intra_bgrp_comm)
       CALL mp_sum (pi_q_unrelax_, intra_bgrp_comm)
       CALL mp_sum (sh_q,          intra_bgrp_comm)
       !! Sum over different processes
       !
       vki_u(iwann) = vki_u(iwann) + pi_q_unrelax
       sh(iwann)    = sh(iwann)    + sh_q
       ! Accumulate over q points
       !
       !
       spin_ref = spin_component
       drhog_scf = CMPLX(0.D0,0.D0,kind=DP)
       !
       IF (do_real_space) THEN
         !
         drhor_scf = CMPLX(0.D0,0.D0,kind=DP)
         CALL solve_linter_koop ( spin_ref, iwann, delta_vr_, drhog_scf, delta_vg_, drhor_scf )
         ! Also the density response in Real space is passed back. Just for check 
         !
       ELSE
         !
         CALL solve_linter_koop ( spin_ref, iwann, delta_vr_, drhog_scf, delta_vg_ )
         !
       ENDIF
       !
       pi_q_relax    = ZERO
       pi_q_relax_rs = ZERO
       DO is =1, nspin_mag
         !
         pi_q_relax = pi_q_relax + sum (CONJG(drhog_scf (:,is)) * delta_vg(:,is))*weight(iq)*omega
         IF (do_real_space) pi_q_relax_rs = pi_q_relax_rs &  
               + sum (CONJG(drhor_scf (:,is)) * delta_vr(:,is))/( dffts%nr1*dffts%nr2*dffts%nr3 )*omega
         !
       ENDDO
       !
       CALL mp_sum (pi_q_relax, intra_bgrp_comm)
       !
       pi_q_relax = pi_q_relax + pi_q_unrelax_
       !
       IF (lgamma .and. l_vcut ) THEN 
          !
          !WRITE(stdout, 9013) iq, iwann, pi_q_relax, pi_q_unrelax_
          WRITE(stdout, '(/, 7X, "INFO: Result without q+G=0   ", 3x, "rPi*      ", 1F15.8, 3x,"uPi*      ", 1F15.8)') &
                         DBLE(pi_q_relax), DBLE(pi_q_unrelax_)
          WRITE(stdout, '(   7X, "INFO: The q+G=0 contribution", 3x, "rPi(q+G=0)", 1F15.8, 3x,"uPi(q+G=0)", 1F15.8)') &
                         -(div_eps)/omega/nqs, -div/omega/nqs
          pi_q_relax = pi_q_relax - (div_eps)/omega/nqs 
          !
       ENDIF
       !
       IF (do_real_space) THEN 
         !
         CALL mp_sum (pi_q_relax_rs, intra_bgrp_comm)
         pi_q_relax_rs = pi_q_relax_rs + pi_q_unrelax
         IF (lgamma .and. l_vcut ) pi_q_relax_rs = pi_q_relax_rs - (div_eps)/omega/nqs 
         !
         WRITE(stdout, 9010) iq, iwann, pi_q_relax, pi_q_relax_rs, pi_q_unrelax, sh_q
         IF( ionode) WRITE(iun_res,'(2I5,8F20.12)') iq, iwann, pi_q_relax, pi_q_relax_rs, pi_q_unrelax, sh_q
         !
       ELSE
         WRITE(stdout, 9011) iq, iwann, pi_q_relax, pi_q_unrelax, sh_q
         IF( ionode) WRITE(iun_res,'(2I5,6F20.12)') iq, iwann, pi_q_relax, pi_q_unrelax, sh_q
       ENDIF
       !
       vki_r(iwann) = vki_r(iwann) + pi_q_relax
       ! Sum over q points
       !
       !
       IF (irr_bz) CALL clean_pw_kcw()
       !
    ENDDO
    !
    IF (.NOT. irr_bz) CALL clean_pw_kcw( )
    DEALLOCATE ( rhog , delta_vg, vh_rhog, drhog_scf, delta_vg_ )
    !
    IF (ionode) THEN 
      INQUIRE(file=TRIM(tmp_dir_kcw)//TRIM(prefix)//'.alpha.status', exist=exst)
      IF (.NOT. exst) OPEN(986, file=TRIM(tmp_dir_kcw)//TRIM(prefix)//'.alpha.status')
      REWIND(986)
      WRITE(986,'(i5)') iq
    ENDIF
    !
  ENDDO ! qpoints
  !
  !
  WRITE(stdout,'(/)')
  WRITE( stdout, '(5X,"INFO: LR CALCULATION ... DONE")')
  WRITE(stdout,'(/)')
  !
  IF (ionode) CLOSE (986, STATUS='DELETE')
  IF (ionode) CLOSE (iun_res, STATUS='KEEP')
  !
  !
  IF (do_real_space) DEALLOCATE ( drhor_scf ) 
  !
  WRITE(stdout,'(/)') 
  !
  IF ( i_orb == -1 ) THEN 
    OPEN (876, file = TRIM(tmp_dir_kcw)//TRIM(prefix)//'.alpha.dat')
    WRITE(876,'(i5)') num_wann
  ENDIF
  !
  DO jwann = iorb_start, iorb_end
    !
    iwann = group_alpha(jwann)
    alpha = REAL(vki_r(iwann))/REAL(vki_u(iwann))
    !
    IF ( l_do_alpha(jwann) ) THEN 
      WRITE(stdout,'(/, 8x, "iwann  = ", i5, 3x, "relaxed = ", 1f15.8, 3x, "unrelaxed = " & 
                &  , 1f15.8, 3x, "alpha =", f12.8, 3x, "self Hartree =", f12.8 )') &
                &  jwann, REAL(vki_r(iwann)), REAL(vki_u(iwann)), alpha, REAL(sh(iwann))
    ELSE
      WRITE(stdout,'(/, 8x, "iwann* = ", i5, 3x, "relaxed = ", 1f15.8, 3x, "unrelaxed = " & 
                 &  , 1f15.8, 3x, "alpha =", f12.8, 3x, "self Hartree =", f12.8 )') &
                 &  jwann, REAL(vki_r(iwann)), REAL(vki_u(iwann)), alpha, REAL(sh(iwann))
    ENDIF
    !
    ! store the final value of alpha
    alpha_final(iwann) = alpha
    IF (i_orb == -1) WRITE(876,'(i5, 2(3x, F16.12))') iwann, alpha, REAL(sh(iwann))
    !
  ENDDO
  !
  WRITE(stdout, '(3/)')
  IF ( i_orb == -1 ) CLOSE (876)
  !
  CALL kcw_deallocate_symmetry_arrays()
  !
9010 FORMAT(/, 8x, "iq =", i4, 3x, "iwann =", i4, 3x, "rPi_q =", 2f15.8, 3x, & 
               "rPi_q_RS =", 2f15.8, 3x, "uPi_q =", 2f15.8, 3x, "Self Hartree =", 2f15.8)
9011 FORMAT(/, 8x, "iq =", i4, 3x, "iwann =", i4, 3x, "rPi_q =", 2f15.8, 3x, "uPi_q =", & 
               2f15.8, 3x, "SH_q =", 2f15.8)
9012 FORMAT(/, 8x, "iq =", i4, 3x, "iwann =", i4, 3x, "rPi_q =", 2f15.8, 3x, "uPi_q =", & 
               2f15.8)
9013 FORMAT(/, 8x, "iq =", i4, 3x, "iwann =", i4, 3x, "rPi_q* =", 2f15.8, 3x, "uPi_q* =", & 
               2f15.8)

END subroutine screen_coeff


!-----------------------------------------------------------------------
SUBROUTINE restart_screen (num_wann, iq_start, vki_r, vki_u, sh, do_real_space)
  !-----------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : intra_image_comm
  USE control_kcw,      ONLY : l_do_alpha, iorb_start, iorb_end, tmp_dir_kcw
  USE io_files,         ONLY : prefix
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: vki_u(num_wann), sh(num_wann), vki_r(num_wann)
  ! ki unrelaxed and relaxed potential and self-hartree
  COMPLEX(DP) :: pi_q_unrelax, sh_q, pi_q_relax, pi_q_relax_rs
  !
  INTEGER :: iq_start, num_wann, iwann, iwann_, iq, iq_, iun_res
  LOGICAL :: exst, do_real_space
  !
  INQUIRE(file=TRIM(tmp_dir_kcw)//TRIM(prefix)//'.alpha.status', exist=exst)
  IF( .NOT. exst) THEN
    RETURN
  ELSE 
    OPEN (986, file = TRIM(tmp_dir_kcw)//TRIM(prefix)//'.alpha.status')
  ENDIF 
  !
  iun_res = 987
  OPEN (iun_res, file = TRIM(tmp_dir_kcw)//TRIM(prefix)//'.LR_res.txt')
  IF (ionode) THEN 
    !
    READ(986, '(i5)') iq_start
    iq_start = iq_start+1
    !
    WRITE(stdout, '(5X, "restart FOUND. Results up to now:")')
    DO iq = 1, iq_start -1 
      WRITE (stdout, '(/)')
      DO iwann = iorb_start, iorb_end
        IF ( .NOT. l_do_alpha (iwann)) CYCLE
        IF (do_real_space) THEN 
          READ(iun_res, '(2I5,8F20.12)') iq_, iwann_ , pi_q_relax, pi_q_relax_rs, pi_q_unrelax, sh_q
          WRITE(stdout, 9011) iq, iwann, pi_q_relax, pi_q_relax_rs, pi_q_unrelax, sh_q
          !
          vki_u(iwann) = vki_u(iwann) + pi_q_unrelax
          vki_r(iwann) = vki_r(iwann) + pi_q_relax
          sh(iwann)    = sh(iwann)    + sh_q
          !
        ELSE
          READ(iun_res, '(2I5,6F20.12)') iq_, iwann_ , pi_q_relax, pi_q_unrelax, sh_q
          WRITE(stdout, 9011) iq, iwann, pi_q_relax, pi_q_unrelax, sh_q
          !
          vki_u(iwann) = vki_u(iwann) + pi_q_unrelax
          vki_r(iwann) = vki_r(iwann) + pi_q_relax
          sh(iwann)    = sh(iwann)    + sh_q
          !
        ENDIF
      ENDDO 
      !
    ENDDO
    ! 
  ENDIF
  ! 
  call mp_bcast ( iq_start,  ionode_id, intra_image_comm ) 
  call mp_bcast ( vki_u,     ionode_id, intra_image_comm ) 
  call mp_bcast ( vki_r,     ionode_id, intra_image_comm ) 
  call mp_bcast ( sh,        ionode_id, intra_image_comm ) 
  RETURN
  !
9010 FORMAT(/, 8x, "iq =", i4, 3x, "iwann =", i4, 3x, "rPi_q =", 2f15.8, 3x, & 
               "rPi_q_RS =", 2f15.8, 3x, "uPi_q =", 2f15.8, 3x, "Self Hartree =", 2f15.8)
9011 FORMAT(/, 8x, "iq =", i4, 3x, "iwann =", i4, 3x, "rPi_q =", 2f15.8, 3x, "uPi_q =", & 
               2f15.8, 3x, "SH_q =", 2f15.8)
  !
END SUBROUTINE restart_screen 

!-----------------------------------------------------------------------
SUBROUTINE check_density (rhowann)
  !-----------------------------------------------------------------------
  !
  ! This routine check that the density computed from rho_wann is consistent
  ! with the one read from PWSCF
  !
  USE kinds,                ONLY : DP
  USE control_kcw,          ONLY : num_wann_occ, kcw_iverbosity, num_wann, spin_component, nrho
  USE fft_base,             ONLY : dffts, dfftp
  USE mp,                   ONLY : mp_sum
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE cell_base,            ONLY : omega
  USE io_global,            ONLY : stdout
  USE scf,                  ONLY: rho
  USE lsda_mod,             ONLY : nspin
  !
  IMPLICIT NONE
  !
  INTEGER :: iwann,ii
  REAL(DP) :: int_rho, w1

  COMPLEX(DP) :: density(dffts%nnr)
  COMPLEX(DP), INTENT (IN) :: rhowann(dffts%nnr, num_wann,nrho)
  REAL(DP),    ALLOCATABLE :: rhoup(:), rhodw(:)
  REAL(DP),    ALLOCATABLE :: int_rhos(:), rhos(:,:), checks(:)
  COMPLEX(DP),    ALLOCATABLE  :: densities(:,:)
  !
  if (nspin==4) then
    ALLOCATE (rhos(dfftp%nnr, nrho))
    ALLOCATE (densities(dfftp%nnr, nrho))
    densities(:,:) = (0.D0,0.D0)
    ALLOCATE (int_rhos(nrho),checks(nrho))
  else
    ALLOCATE (rhoup(dfftp%nnr), rhodw(dfftp%nnr))
  endif
  !
  density = (0.D0,0.D0)
  DO iwann = 1, num_wann_occ ! sum over the bands (rhowann already summed over k) 
    if (nspin==4) then
        w1 = 1.D0/omega
        do ii=1,nspin 
            densities(:,ii) = densities(:,ii) + w1 * rhowann(:,iwann,ii)
        end do
    else
        w1 = 2.D0/nspin  / (omega)    ! true only if wg is the same for each k point (no symm)
        density(:) = density(:) + w1 * rhowann(:,iwann,1)
    endif 

  ENDDO
  !
  if (nspin == 1) then
     !
     rhoup(:) =  (rho%of_r (:, 1) )
     rhodw(:) =  (rho%of_r (:, 1) )
     !
  elseif (nspin == 2) then
     !
     rhoup(:) = ( rho%of_r(:, 1) + rho%of_r(:, 2) )*0.5d0
     rhodw(:) = ( rho%of_r(:, 1) - rho%of_r(:, 2) )*0.5d0
     !
  else 
     !
      rhos(:,1) =  (rho%of_r (:, 1) )
      rhos(:,2) =  (rho%of_r (:, 2) )
      rhos(:,3) =  (rho%of_r (:, 3) )
      rhos(:,4) =  (rho%of_r (:, 4) )
    
  endif
  !
  if (nspin==4) then
      do ii=1,nrho
        checks(ii) = SUM( REAL(densities(:,ii)))/(dffts%nr1*dffts%nr2*dffts%nr3)
        int_rhos(ii) = SUM( REAL(densities(:,ii))-rhos(:,ii))/(dffts%nr1*dffts%nr2*dffts%nr3)
        CALL mp_sum( int_rhos(ii), intra_bgrp_comm)
        CALL mp_sum( checks(ii), intra_bgrp_comm)
      enddo
  else
      IF (spin_component == 1 ) int_rho = SUM( REAL(density(:))-rhoup(:))/(dffts%nr1*dffts%nr2*dffts%nr3) 
      IF (spin_component == 2 ) int_rho = SUM( REAL(density(:))-rhodw(:))/(dffts%nr1*dffts%nr2*dffts%nr3) 
      CALL mp_sum( int_rho, intra_bgrp_comm )
  end if


  !
  if (nspin==4) then
    WRITE(stdout,'(8X, "DEBUG: check 1 = ",E18.6,/)') checks(1)
    WRITE(stdout,'(8X, "DEBUG: check 2 = ",E18.6,/)') checks(2)
    WRITE(stdout,'(8X, "DEBUG: check 3 = ",E18.6,/)') checks(3)
    WRITE(stdout,'(8X, "DEBUG: check 4 = ",E18.6,/)') checks(4)

    do ii=1,nrho
      IF (kcw_iverbosity > 1 ) WRITE(stdout,'(8X, "DEBUG: rho component (range 1-4)= ",I1, /)') ii
      IF (kcw_iverbosity > 1 ) WRITE(stdout,'(8X, "DEBUG: \int dr [rho - rho_PWSCF] = ",E18.6, /)') int_rhos(ii)
      IF ( ABS(int_rhos(ii)) .gt. 1e-4) THEN
        WRITE(stdout,'(8X, "DEBUG: rho component (range 1-4): ",I1,/)') ii
        WRITE(stdout,'(8X, "DEBUG: \int dr [rho - rho_PWSCF] = ",E18.6,/)') int_rhos(ii)
        CALL errore ('check_density','\int dr [rho - rho_PWSCF] > 1e-4; SOMETHING WRONG',1)
      ENDIF
    enddo
    DEALLOCATE (rhos)
  else
    IF (kcw_iverbosity > 1 ) WRITE(stdout,'(8X, "DEBUG: \int dr [rho - rho_PWSCF] = ",E18.6, /)') int_rho
    IF ( ABS(int_rho) .gt. 1e-8) THEN
       WRITE(stdout,'(8X, "DEBUG: \int dr [rho - rho_PWSCF] = ",E18.6,/)') int_rho
       CALL errore ('check_density','\int dr [rho - rho_PWSCF] > 1e-8; SOMETHING WRONG',1)
    ENDIF
    DEALLOCATE (rhoup, rhodw)
  endif
  !
END subroutine check_density

