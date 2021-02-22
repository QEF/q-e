! Copyright (C) 2001 PWSCF group
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
  USE control_kc_wann,      ONLY : kc_iverbosity, spin_component, num_wann, iorb_start, l_do_alpha, &
                                   iorb_end, alpha_final, nqstot, eps_inf, l_vcut, l_unique_manifold, group_alpha
  USE buffers,              ONLY : get_buffer, save_buffer
  USE io_global,            ONLY : stdout, ionode
  USE control_kc_wann,      ONLY : iurho_wann
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE cell_base,            ONLY : at
  USE dv_of_drho_lr,        ONLY : dv_of_drho
  USE fft_interfaces,       ONLY : fwfft
  USE control_lr,           ONLY : lgamma
  USE disp,                 ONLY : x_q
  USE lsda_mod,             ONLY : nspin
  USE grid_irr_iq,          ONLY : done_bands
  USE gvecs,                ONLY : ngms
  USE io_files,             ONLY : tmp_dir, prefix
  USE save_ph,              ONLY : tmp_dir_save
  USE control_ph,           ONLY : tmp_dir_phq
  USE solve_linter_koop_mod 
  USE exx_base,             ONLY : exx_divergence
  !
  !USE mp_world,             ONLY : mpime
  !
  USE cell_base,            ONLY : omega
  USE el_phon,              ONLY : elph_mat
  ! 
  !
  IMPLICIT NONE
  ! 
  INTEGER :: iq, nqs, spin_ref, is
  !! Counter for the k/q points in the BZ, total number of q points and number of pw for a given k (k+q) point
  !
  INTEGER :: iwann, jwann, lrrho, iq_start, iun_res
  !! Band counter, leght of the rho record, starting iq (if restart), iunit partial results
  !
  COMPLEX(DP) :: rhowann(dffts%nnr, num_wann), rhor(dffts%nnr), delta_vr(dffts%nnr,nspin), delta_vr_(dffts%nnr,nspin)
  !! The periodic part of the wannier orbital density
  !
  COMPLEX(DP) :: vki_u(num_wann), sh(num_wann), vki_r(num_wann)
  ! ki unrelaxed and relaxed potential and self-hartree
  !
  COMPLEX(DP), ALLOCATABLE  :: rhog(:), delta_vg(:,:), vh_rhog(:), drhog_scf(:,:), drhor_scf(:,:), delta_vg_(:,:)
  ! wanier density, perturbing potential, hartree potential, density variation (g and r), perturbing pot with G0=0 
  !
  LOGICAL :: do_band, do_iq, setup_pw, elph_mat_save
  !
  COMPLEX(DP) :: phase(dffts%nnr), wann_c(dffts%nnr,num_wann), rho_c(dffts%nnr,num_wann)
  !! The phase associated to the hift k+q-> k'
  !
  COMPLEX(DP) :: int_rho, int_wann, pi_q_unrelax, sh_q, pi_q_relax, pi_q_relax_rs
  COMPLEX(DP) :: struct_fact
  LOGICAL :: do_real_space = .false. 
  !
  COMPLEX(DP) :: drho_zero
  !
  REAL(DP) :: xq_(3), weight(nkstot)
  REAL(DP) :: alpha
  REAL(DP) :: div
  !
  nqs = nqstot
  !
  IF (nqs == 1) do_real_space = .TRUE. 
  IF (do_real_space) THEN 
     ALLOCATE ( drhor_scf(dffts%nnr,nspin) ) 
     drhor_scf = ZERO
  ENDIF
  ! If only 1 K point do also integral in REAL space (Extra check for the Super cell calculation)
  !
  ALLOCATE(done_bands(nqs))
  done_bands= .FALSE.
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
  div = exx_divergence()
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
  OPEN (iun_res, file = TRIM(tmp_dir)//TRIM(prefix)//'.LR_res.txt')
  CALL restart_screen (num_wann, iq_start, vki_r, vki_u, sh, do_real_space)
  !
  DO iq = iq_start, nqs
    !! For each q in the mesh 
    !
    CALL kc_prepare_q ( do_band, do_iq, setup_pw, iq )
    IF (kc_iverbosity .gt. -1 ) WRITE(stdout,'(8x, "INFO: prepare_q DONE",/)') 
    !! Prepare the q point calculation
    !
    lrrho=num_wann*dffts%nnr
    tmp_dir = tmp_dir_save  ! the periodic part are written on the original outdir 
    CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
    tmp_dir = tmp_dir_phq   ! go back to the q-specific directory
    ! ... Retrive the rho_wann_q(r) from buffer in REAL space
    !
    IF (kc_iverbosity .gt. -1 ) WRITE(stdout,'(8X, "INFO: rhowan_q(r) RETRIEVED"/)') 
    !
    ! FIXME: here elph_mat (global var) is used to skip the chek on equivalent k-points
    elph_mat_save = elph_mat
    elph_mat = .TRUE. 
    IF (setup_pw) CALL run_nscf(do_band, iq)
    elph_mat = elph_mat_save
    !
    IF (kc_iverbosity .gt. -1 .AND. setup_pw) WRITE(stdout,'(/,8X, "INFO: NSCF calculation DONE",/)')
    ! ... IF needed run a nscf calculation for k+q
    ! ... NB: since the k+q is always a p we already have, this could be 
    !         avoided in principle (see compute_map and rho_of_q). 
    !         For the moment it's easier to  keep it. 
    !
    CALL kc_initialize_ph ( ) 
    IF (kc_iverbosity .gt. -1 ) WRITE(stdout,'(8X, "INFO: kc_q initialization DONE",/)') 
    !
    ALLOCATE ( rhog (ngms) , delta_vg(ngms,nspin), vh_rhog(ngms), drhog_scf (ngms, nspin), delta_vg_(ngms,nspin) )
    !
    IF ( lgamma .AND. .NOT. l_unique_manifold) CALL check_density (rhowann) 
    !! ... For q==0 the sum over k and v should give the density. If not something wrong...
    !
    weight(iq) = 1.D0/nqs  !*nspin ! No SYMM  ! CHECK nspin???
    !
    WRITE(stdout, '("weight =", i5, f12.8)') iq, weight(iq)
    !
    DO iwann = iorb_start, iorb_end  ! for each wannier, that is actually the perturbation
       !
       IF ( .NOT. l_do_alpha (iwann)) CYCLE
       ! Skip LR calculation if this orbital match with a one already computed (see group_orbital)
       !
       drhog_scf (:,:) = ZERO
       rhog(:)         = ZERO
       delta_vg(:,:)   = ZERO
       delta_vg_(:,:)  = ZERO
       vh_rhog(:)      = ZERO
       rhor(:)         = ZERO
       !
       rhor(:) = rhowann(:,iwann)
       !! ... The periodic part of the orbital desity in real space
       !
       CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_) 
       !! ... The periodic part of the perturbation
       !
       pi_q_unrelax = sum (CONJG(rhog (:)) * delta_vg(:,spin_component))*weight(iq)*omega
       sh_q = 0.5D0 * sum (CONJG(rhog (:)) * vh_rhog(:)                )*weight(iq)*omega
       !
       CALL mp_sum (pi_q_unrelax, intra_bgrp_comm)
       CALL mp_sum (sh_q,         intra_bgrp_comm)
       !! Sum over different processes
       !
       vki_u(iwann) = vki_u(iwann) + pi_q_unrelax
       sh(iwann)    = sh(iwann)    + sh_q
       ! Accumulate over q points
       !
       WRITE( stdout, '(/, 8x,"Start linear response calculation for the wannier #",i3, "    spin =", i3)') iwann, spin_component
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
       DO is =1, nspin
         !
         pi_q_relax = pi_q_relax + sum (CONJG(drhog_scf (:,is)) * delta_vg(:,is))*weight(iq)*omega
         IF (do_real_space) pi_q_relax_rs = pi_q_relax_rs &  
               + sum (CONJG(drhor_scf (:,is)) * delta_vr(:,is))/( dffts%nr1*dffts%nr2*dffts%nr3 )*omega
         !
       ENDDO
       !
       CALL mp_sum (pi_q_relax, intra_bgrp_comm)
       !
       IF (lgamma .and. l_vcut ) THEN 
          IF ( ABS(eps_inf - 1.D0) .gt. 1e-6) THEN 
            WRITE(stdout, '(/, 5X, "INFO: Adding the q+G contribution", 3x, "bare", 1es15.5, 3x,"screen", 1es15.5)') &
                            -div/omega/nqs, -(1.D0/eps_inf -1.D0)*div/omega/nqs
            pi_q_relax = pi_q_relax - (1.D0/eps_inf -1.D0)*div/omega/nqs 
          ELSE
            ! This is in the case eps_inf is not provided in input but l_vcut=.true. .
            ! The LR is done with G0=0 while the unrelaxed with GB for G0=0. This remove the inconsistency
            pi_q_relax = pi_q_relax  + div/omega/nqs 
          ENDIF
       ENDIF
       !
       IF (do_real_space) THEN 
         !
         CALL mp_sum (pi_q_relax_rs, intra_bgrp_comm)
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
    ENDDO
    !
    DEALLOCATE ( rhog , delta_vg, vh_rhog, drhog_scf, delta_vg_ )
    !
    CALL clean_pw_kc(iq)
    !    
    IF (ionode) REWIND(986)
    IF (ionode) WRITE(986,'(i5)') iq
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
  DO jwann = iorb_start, iorb_end
    !
    iwann = group_alpha(jwann)
    alpha = REAL(vki_u(iwann)+vki_r(iwann))/REAL(vki_u(iwann))
    !
    IF (l_do_alpha(jwann) ) THEN 
      WRITE(stdout,'(/, 8x, "iwann  = ", i5, 3x, "relaxed = ", 1f15.8, 3x, "unrelaxed = " & 
                &  , 1f15.8, 3x, "alpha =", f12.8, 3x, "self Hartree =", f12.8 )') &
                &  jwann, REAL(vki_u(iwann)+vki_r(iwann)), REAL(vki_u(iwann)), alpha, REAL(sh(iwann))
    ELSE
      WRITE(stdout,'(/, 8x, "iwann* = ", i5, 3x, "relaxed = ", 1f15.8, 3x, "unrelaxed = " & 
                 &  , 1f15.8, 3x, "alpha =", f12.8, 3x, "self Hartree =", f12.8 )') &
                 &  jwann, REAL(vki_u(iwann)+vki_r(iwann)), REAL(vki_u(iwann)), alpha, REAL(sh(iwann))
    ENDIF
    !
    ! store the final value of alpha
    alpha_final(iwann) = alpha
    !
  ENDDO
  !
  WRITE(stdout, '(3/)')
  !
  !! This is TEMPORARY to check the rho_q(r) are correctly computed 
  rho_c = ZERO
  wann_c = ZERO
  DO iq = 1, nqs 
    !
    tmp_dir = tmp_dir_save  ! the periodic part are written on the original outdir 
    lrrho=num_wann*dffts%nnr
    CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
    tmp_dir = tmp_dir_phq   ! go back to the q-specific directory
    !
    xq_(:) = -x_q(:,iq)
    ! calculate_phase has a - sign inside
    phase=ZERO
    CALL calculate_phase(xq_, phase)
    CALL structure_factor(iq, struct_fact)
    CALL cryst_to_cart(1, xq_, at, -1)
    IF (kc_iverbosity .gt. 1) WRITE(stdout,'(8X, "INFO: iq = ", i5, 3x, "Structure Factor S(q) [Re, Im] = ", 2f12.8,/)') & 
                                                                iq, struct_fact
    !
     DO iwann=iorb_start, iorb_end
       !
       wann_c(:,iwann) = wann_c(:,iwann) + phase(:)*rhowann(:,iwann)*weight(iq)
       rho_c(:,iwann)  = rho_c(:,iwann)  + phase(:)*rhowann(:,iwann)*struct_fact*weight(iq)
       !
    ENDDO
    !
  ENDDO ! qpoints
  !
  IF (kc_iverbosity .gt. 1) THEN
  !  write(*,'(/,"DEBUG")')
    DO iwann= iorb_start, iorb_end
       int_rho = (0.D0,0.D0)
       int_rho = SUM (rho_c(:,iwann))/(dffts%nr1*dffts%nr2*dffts%nr3)
       int_wann = (0.D0,0.D0)
       int_wann = SUM (wann_c(:,iwann))/(dffts%nr1*dffts%nr2*dffts%nr3)
       CALL mp_sum( int_rho, intra_bgrp_comm )
       CALL mp_sum( int_wann, intra_bgrp_comm )
       WRITE(stdout,'(8X, "iwann= ", i3, 3x, "int rho_wann(r) [Re, Im] =", 2f12.6)') iwann, int_rho
       WRITE(stdout,'(8X, "iwann= ", i3, 3x, "int Im[rho_wann(r)]      =", 1f12.6)') iwann, AIMAG(int_wann)
    ENDDO
  ENDIF
  !


9010 FORMAT(/, 8x, "iq =", i4, 3x, "iwann =", i4, 3x, "rPi_q =", 2f15.8, 3x, & 
               "rPi_q_RS =", 2f15.8, 3x, "uPi_q =", 2f15.8, 3x, "Self Hartree =", 2f15.8)
9011 FORMAT(/, 8x, "iq =", i4, 3x, "iwann =", i4, 3x, "rPi_q =", 2f15.8, 3x, "uPi_q =", & 
               2f15.8, 3x, "SH_q =", 2f15.8)

END subroutine screen_coeff


!-----------------------------------------------------------------------
SUBROUTINE restart_screen (num_wann, iq_start, vki_r, vki_u, sh, do_real_space)
  !-----------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : ionode, ionode_id, stdout
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : intra_image_comm
  USE control_kc_wann,  ONLY : l_do_alpha, iorb_start, iorb_end
  USE io_files,         ONLY : tmp_dir, prefix
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
  INQUIRE(file="fort.986", exist=exst)
  IF( .NOT. exst) RETURN
  !
  iun_res = 987
  OPEN (iun_res, file = TRIM(tmp_dir)//TRIM(prefix)//'.LR_res.txt')
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
  ! This routine check that the density computed from rho_wann cis consisten
  ! with the one read from PWSCF
  !
  USE kinds,                ONLY : DP
  USE control_kc_wann,      ONLY : num_wann_occ, kc_iverbosity, num_wann, spin_component
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
  INTEGER :: iwann
  REAL(DP) :: int_rho, w1
  COMPLEX(DP) :: density(dffts%nnr)
  COMPLEX(DP), INTENT (IN) :: rhowann(dffts%nnr, num_wann)
  REAL(DP),    ALLOCATABLE :: rhoup(:), rhodw(:)
  !
  ALLOCATE (rhoup(dfftp%nnr), rhodw(dfftp%nnr))
  !
  density = (0.D0,0.D0)
  DO iwann = 1, num_wann_occ ! sum over the bands (rhowann already summed over k) 
     w1 = 2.D0/nspin  / (omega)    ! true only if wg is the same for each k point (no symm)
     density(:) = density(:) + w1 * rhowann(:,iwann)
  ENDDO
  !
  if (nspin == 1) then
     !
     rhoup(:) =  (rho%of_r (:, 1) )
     rhodw(:) =  (rho%of_r (:, 1) )
     !
  else
     !
     rhoup(:) = ( rho%of_r(:, 1) + rho%of_r(:, 2) )*0.5d0
     rhodw(:) = ( rho%of_r(:, 1) - rho%of_r(:, 2) )*0.5d0
     !
  endif
  !
  IF (spin_component == 1 ) int_rho = SUM( REAL(density(:))-rhoup(:))/(dffts%nr1*dffts%nr2*dffts%nr3) 
  IF (spin_component == 2 ) int_rho = SUM( REAL(density(:))-rhodw(:))/(dffts%nr1*dffts%nr2*dffts%nr3) 
  
  CALL mp_sum( int_rho, intra_bgrp_comm )
  !
  IF (kc_iverbosity > 1 ) WRITE(stdout,'(8X, "DEBUG: \int dr [rho - rho_PWSCF] = ",E18.6, /)') int_rho
  IF ( ABS(int_rho) .gt. 1e-8) THEN
     WRITE(stdout,'(8X, "DEBUG: \int dr [rho - rho_PWSCF] = ",E18.6,/)') int_rho
     CALL errore ('check_density','\int dr [rho - rho_PWSCF] > 1e-8; SOMETHING WRONG',1)
  ENDIF
  DEALLOCATE (rhoup, rhodw)
  !
END subroutine

