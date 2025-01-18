! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO CMPLX(0.D0,0.D0,kind=DP)
#define ONE CMPLX(0.D0,1.D0, kind=DP)
!#define DEBUG

!-----------------------------------------------------------------------
SUBROUTINE ham_koopmans_k (ik)
  !---------------------------------------------------------------------
  !
  !! This routine compute the KI correction to second order to the KS
  !! Hamiltonian at a given k point ik given from input
  !
  !
  USE kinds,                ONLY : DP
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE fft_wave,             ONLY : invfft_wave
  USE klist,                ONLY : nkstot, igk_k, ngk, xk
  USE mp,                   ONLY : mp_sum
  USE control_kcw,          ONLY : kcw_iverbosity, spin_component, num_wann, x_q, &
                                   alpha_final, num_wann_occ, evc0, iuwfc_wann, &
                                   map_ikq, shift_1bz, nqstot, Hamlt, l_alpha_corr, alpha_corr_done
  USE buffers,              ONLY : get_buffer, save_buffer
  USE io_global,            ONLY : stdout, ionode
  USE control_kcw,          ONLY : iurho_wann
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE dv_of_drho_lr,        ONLY : dv_of_drho
  USE control_lr,           ONLY : lgamma
  USE lsda_mod,             ONLY : nspin
  USE gvecs,                ONLY : ngms
  USE solve_linter_koop_mod 
  USE qpoint,               ONLY : xq
  USE wvfct,                ONLY : npwx 
  USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux, nspin_lsda, nspin_gga, nspin_mag, npol
  !
  !USE mp_world,             ONLY : mpime
  !
  USE cell_base,            ONLY : omega
  USE noncollin_module,     ONLY : npol
  !
  USE constants,            ONLY : rytoev
  ! 
  !
  IMPLICIT NONE
  ! 
  INTEGER, INTENT(IN) :: ik
  ! the k-point index in hte orignal mesh of k-points
  !
  INTEGER :: iq, nqs
  ! Counter for the k/q points in the BZ, total number of q points and number of pw for a given k (k+q) point
  !
  INTEGER :: iwann, jwann, lrrho, lrwfc, i
  ! Band counters, leght of the rho record
  !
  COMPLEX(DP) :: rhowann(dffts%nnr, num_wann), rhor(dffts%nnr), sh(num_wann)
  COMPLEX(DP) :: delta_vr(dffts%nnr,nspin_mag), delta_vr_(dffts%nnr,nspin_mag)
  ! The periodic part of the wannier orbital density in r space
  ! The self-hartree 
  ! The perturbig potential in real space
  ! The perturbig potential in real space (without the g=0 contribution) 
  !
  COMPLEX(DP), ALLOCATABLE  :: rhog(:), delta_vg(:,:), vh_rhog(:), delta_vg_(:,:)
  ! The periodic parte of the wannier density in G-space 
  ! The perturbig potential in G space
  ! The hartree potential due to the wannier in G-space 
  ! The perturbig potential in G space without the g=0 contribution 
  !
  COMPLEX(DP), ALLOCATABLE :: deltaH(:,:)
  ! The KC controbution to the Hamiltonian 
  !
  REAL(DP) :: weight(nqstot), alpha_(num_wann), alpha_fd
  ! weight of the q points, alpha from LR, alpha after the all-order correction. 
  !
  REAL(DP) second_der(num_wann) 
  !
  INTEGER :: ikq, npw_k, npw_kq
  ! Counter for the k/q points in the BZ, total number of q points and number of pw for a given k (k+q) point
  !
  COMPLEX(DP) ::  evc0_kq(npwx*npol, num_wann)
  ! Auxiliary vector to store the wfc at k+q
  !
  COMPLEX(DP) :: rho_r_nm(dffts%nnr), rho_g_nm(ngms), aux(dffts%nnr)
  REAL(DP) :: g_vect(3)
  ! G vector that shift the k+q inside the 1BZ
  !
  COMPLEX(DP) :: evc_k_g (npwx*npol), evc_k_r (dffts%nnr), phase(dffts%nnr)
  ! Auxiliary wfc in reciprocal and real space, the phase associated to the hift k+q-> k'
  !
  COMPLEX(DP) :: evc_kq_g (npwx*npol), evc_kq_r (dffts%nnr)
  ! Auxiliary wfc in reciprocal and real space
  !
  LOGICAL :: off_diag = .TRUE.
  ! compute Off-diagonal elements 
  !
  COMPLEX(DP) :: ham(num_wann,num_wann), eigvc(npwx,num_wann)
  ! the hamiltonian at a given k
  !
  REAL(DP) :: eigvl(num_wann), delta
  ! KC eigenvalues 
  ! auxiliary variable for alpha correction
  !
  LOGICAL :: corr_done=.FALSE.
  ! whether the correction to the current wannier was already done 
  !
  if (nspin_mag == 4) &
     CALL errore ('hamilt', ' ham_koopmans_k not implemented for non-collinear magnetic calculations ', 1)
  nqs = nkstot/nspin_mag !<--- non-collinear not implemented
  nqs = nqstot
  !
  ALLOCATE( deltaH(num_wann,num_wann) )
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
  deltaH = ZERO
  sh = ZERO
  !
  lrwfc = num_wann*npwx
  CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik )
  ! Retrive the ks function at k (in the Wannier Gauge)
  IF (kcw_iverbosity .gt. 0 ) WRITE(stdout,'(8X, "INFO: u_k(g) RETRIEVED"/)') 
  !
  CALL compute_map_ikq_single (ik)
  ! find tha map k+q --> k'+G and store the res 
  !
  WRITE(stdout,'(/)')
  WRITE( stdout, '(5X,"INFO: KC HAMILTONIAN CALCULATION ...")')
  !
  DO iq = 1, nqs
    !! Sum over the BZ 
    !
    xq(1:3)  = x_q(1:3,iq)
    !
    lgamma = ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 )
    !
    lrrho=num_wann*dffts%nnr
    CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
    ! Retrive the rho_wann_q(r) from buffer in REAL space
    IF (kcw_iverbosity .gt. 0 ) WRITE(stdout,'(8X, "INFO: rhowan_q(r) RETRIEVED"/)') 
    !
    ALLOCATE ( rhog (ngms) , delta_vg(ngms,nspin_mag), vh_rhog(ngms), delta_vg_(ngms,nspin_mag) )
    !
    IF ( lgamma ) CALL check_density (rhowann) 
    ! CHECK: For q==0 the sum over k and v should give the density. If not something wrong...
    !
    weight(iq) = 1.D0/nqs ! No SYMM 
    !
    !WRITE(stdout, '("weight =", i5, f12.8)') iq, weight(iq)
    !
    ikq = map_ikq(iq)
    g_vect(:) = shift_1bz(:,iq)
    ! The index ikq in the 1BZ corresponding at the k+q, and the vector G_bar defining the shift 
    ! xk(:,ikq)+G_bar = xk(:,k+q)
    ! see compute_map_ikq_single.f90
    !
    phase(:) = 0.D0
    CALL calculate_phase(g_vect, phase) 
    ! Calculate the phase associated to the k+q-> ikq map: exp[ -i(G_bar * r) ]
    !
    evc0_kq = ZERO
    lrwfc = num_wann * npwx 
    CALL get_buffer ( evc0_kq, lrwfc, iuwfc_wann, ikq )
    ! Retrive the ks function at k+q (in the Wannier Gauge): 
    !
    IF (kcw_iverbosity .gt. 0 ) WRITE(stdout,'(8X, "INFO: u_kq(g) RETRIEVED"/)') 
    !
    ! Occupied state first
    DO iwann = 1, num_wann_occ  
       !
       rhog(:)         = ZERO
       delta_vg(:,:)   = ZERO
       vh_rhog(:)      = ZERO
       rhor(:)         =ZERO
       !
       rhor(:) = rhowann(:,iwann)
       !! The periodic part of the orbital desity R=0, n=iwann in real space
       !
       CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ ) 
       ! The periodic part of the perturbation DeltaV_q(G)
       !
       sh(iwann) = sh(iwann) + 0.5D0 * sum (CONJG(rhog (:)) * vh_rhog(:)                )*weight(iq)*omega
       deltaH(iwann, iwann) = deltaH(iwann,iwann) -0.5D0 * sum (CONJG(rhog (:)) * delta_vg(:,spin_component))*weight(iq)*omega
       ! the diagonal term 
       !
    ENDDO
    !
    ! Empty states
    DO iwann = num_wann_occ+1, num_wann
       !
       rhog(:)         = ZERO
       delta_vg(:,:)   = ZERO
       vh_rhog(:)      = ZERO
       rhor(:)         = ZERO
       !
       rhor(:) = rhowann(:,iwann)
       ! The periodic part of the orbital desity R=0, n=iwann in real space
       !
       CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
       ! The periodic part of the perturbation DeltaV_q(G)
       !
       sh(iwann) = sh(iwann) + 0.5D0 * sum (CONJG(rhog (:)) * vh_rhog(:)                )*weight(iq)*omega
       deltaH(iwann, iwann) = deltaH(iwann,iwann) + 0.5D0 * sum (CONJG(rhog (:)) * delta_vg(:,spin_component))*weight(iq)*omega
       !! the diagonal term 
       !
       !
       npw_k = ngk(ik)
       evc_k_g(:) =  evc0(:,iwann)
       evc_k_r(:) = ZERO
       CALL invfft_wave (npwx, npw_k, igk_k (1,ik), evc_k_g , evc_k_r )
       !! The wfc R=0 n=iwann in R-space at k
       !
       ! Off-diagonal terms
       IF (off_diag) THEN
       DO jwann = iwann+1, num_wann 
          !
          rhog(:)         = ZERO
          delta_vg(:,:)   = ZERO
          vh_rhog(:)      = ZERO
          rhor(:)         = ZERO
          !
          rhor(:) = rhowann(:,jwann)
          ! The periodic part of the orbital desity R=0, n=iwann in real space
          !
          CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
          ! The periodic part of the perturbation DeltaV_q(G)
          !
          npw_kq = ngk(ikq)
          evc_kq_g = evc0_kq(:,jwann)
          evc_kq_r = ZERO
          CALL invfft_wave (npwx, npw_kq, igk_k (1,ikq), evc_kq_g , evc_kq_r )
          ! The wfc in R-space at k' <-- k+q where k' = (k+q)-G_bar
          ! evc_k+q(r) = sum_G exp[iG r] c_(k+q+G) = sum_G exp[iG r] c_k'+G_bar+G 
          !            = exp[-iG_bar r] sum_G' exp[iG'r] c_k'+G' = exp[-iG_bar r] *evc_k'(r)
          !
          rho_r_nm(:) = conjg(evc_k_r(:))*evc_kq_r(:)*phase(:)/nqs
          ! generalized density in r-space 
!          IF (jwann == iwann) THEN
!            WRITE(*,'("NICOLA evc i", 6F20.15)') evc_k_r(1:3)
!            WRITE(*,'("NICOLA evc j", 6F20.15)') evc_kq_r(1:3)
!            WRITE(*,'("NICOLA rho_ij", 6F20.15)') rho_r_nm(1:3)
!            WRITE(*,'("NICOLA rho_ii", 6F20.15)') rhor(1:3)
!          ENDIF
          !WRITE(*,'("NICOLA R", 2i5, 2F20.15)'), iwann, jwann, SUM( delta_vr(:,spin_component)*rho_r_nm(:) )/( dffts%nr1*dffts%nr2*dffts%nr3 )
          !
          aux(:) = rho_r_nm(:)/omega
          CALL fwfft ('Rho', aux, dffts)
          rho_g_nm(:) = aux(dffts%nl(:))
          ! generalized density in G-spage
!          IF (jwann == iwann) THEN
!            WRITE(*,'("NICOLA rho_ij", 6F20.15)') rho_g_nm(1:3)
!            WRITE(*,'("NICOLA rho_ii", 6F20.15)') rhog(1:3)
!            WRITE(*,'("NICOLA vi_g", 6F20.15)') delta_vg(1:3,spin_component)
!            WRITE(*,'("NICOLA i, q r ", i5, 10F16.8)') iwann, SUM (CONJG(rhog (:)) * delta_vg(:,spin_component))*weight(iq)*omega
!            WRITE(*,'("NICOLA i, q r ", i5, 10F16.8)') iwann, SUM (CONJG(rho_g_nm (:)) * delta_vg(:,spin_component))*weight(iq)*omega
!          ENDIF
          !
          deltaH(iwann, jwann) = deltaH(iwann,jwann) + SUM((rho_g_nm(:))*CONJG(delta_vg(:,spin_component)))*weight(iq)*omega
          !deltaH(jwann, iwann) = deltaH(jwann,iwann) + SUM(rho_g_nm(:)*CONJG(delta_vg(:,spin_component)))*weight(iq)*omega
          !WRITE(*,'("NICOLA G", 2i5, 2F20.15)') iwann, jwann, SUM (CONJG(rho_g_nm (:)) * delta_vg(:,spin_component))*weight(iq)*omega
          !WRITE(*,'("NICOLA G", 2i5, 2F20.15)') iwann, jwann, SUM (rho_g_nm (:) * CONJG(delta_vg(:,spin_component)))*weight(iq)*omega
          deltaH(jwann,iwann) = CONJG(deltaH(iwann,jwann))
          !
       ENDDO ! jwann
       ENDIF
       ! 
       !
    ENDDO ! iwann
    !
    DEALLOCATE ( rhog , delta_vg, vh_rhog, delta_vg_ )
    !
    !    
  ENDDO ! qpoints
  !
  CALL mp_sum (deltaH, intra_bgrp_comm)
  CALL mp_sum (sh, intra_bgrp_comm)
  ! Sum over different processes (G vectors) 
  !
  DO iwann = 1, num_wann_occ
    second_der(iwann) = -deltaH(iwann,iwann)
  ENDDO
  DO iwann = num_wann_occ+1, num_wann
    second_der(iwann) = deltaH(iwann,iwann)
  ENDDO
  ! Screening parameter and all-order correction
  ! NEED TO MOVE IT OUTSIDE THE k-point LOOP (DOES NOT DEPEND ON k) !
  !
  DO iwann = 1, num_wann 
    delta =0.D0 
    alpha_(iwann) = alpha_final(iwann) 
    corr_done = alpha_corr_done(iwann)
    !
    ! ... Only if l_alpha_corr =.true. and (for the moment) only for Occupied state ...
    !IF (l_alpha_corr) THEN 
    IF (l_alpha_corr .AND. iwann .le. num_wann_occ .AND. .NOT. corr_done) THEN 
    !IF (l_alpha_corr .AND. .NOT. corr_done) THEN 
      !
      ! ... Compute the difference between the parabolic extrapolation at N \pm 1 and the real 
      ! ... value of the energy in the frozen orbital approximation ...
      CALL alpha_corr (iwann, delta)
      !
      ! ... The new alpha that should be closer to the Finite-difference one ...
      ! ... Remember DeltaH is nothing but the second derivative wrt the orbital occupation ...
      alpha_fd = (alpha_final(iwann)*second_der(iwann) + delta)/ (second_der(iwann)+delta)
      !
      ! ... Since the ham in the frozen approximation is approximated to second
      ! ... order only, this is the alpha we want to use. Only the
      ! ... numerator matters.
      alpha_(iwann) = (alpha_final(iwann)*second_der(iwann) + delta)/second_der(iwann)
      !
      ! ... Write it just to compare with the FD one from CP ... 
      WRITE(stdout,'(5X, "INFO: iwann, LR-alpha, FD-alpha, alpha", i3, 3f12.8)') iwann, alpha_final(iwann),alpha_fd,  alpha_(iwann)
      !
      !WRITE(stdout,'("Nicola", i3, 6f12.8)') iwann, deltaH(iwann,iwann)
      !
      alpha_final(iwann) = alpha_(iwann) 
      WRITE( stdout, '(5X,"INFO: alpha RE-DEFINED ...", i5, f12.8)') iwann, alpha_final(iwann)
      alpha_corr_done (iwann) = .TRUE.
    ENDIF
    
    ! ... Apply the scaling parameter ...
    !deltaH(:,iwann) = alpha_ (iwann)* deltaH(:,iwann) 
    deltaH(:,iwann) = alpha_final (iwann)* deltaH(:,iwann) 
    !deltaH(:,iwann) = alpha_final(iwann) * deltaH(:,iwann) + delta
    !
  ENDDO
  !
  !
  ! The KS hamiltonian in the Wannier Gauge (just to check)
  ham(:,:) = Hamlt(ik,:,:) 
  CALL cdiagh( num_wann, ham, num_wann, eigvl, eigvc )
  !
  Hamlt(ik,:,:) = Hamlt(ik,:,:) + deltaH(:,:)
  !
  IF (ionode) THEN 
    WRITE(stdout, '("Self Hartree, Bare Pi, Relaxed Pi")')
    DO iwann = 1, num_wann
      WRITE(stdout, '(i6, 2f12.8, 3x, 2f12.8,3x, 2f12.8)') iwann, sh(iwann), ABS(deltaH(iwann,iwann))*2.D0/alpha_(iwann), &
                                                           ABS(deltaH(iwann,iwann))*2.D0
    ENDDO
  ENDIF 
  !
  IF (ionode) THEN 
    WRITE(stdout, '("KI Contribution to the Hamiltonian at k = ", i4)') ik
!    DO iwann = 1, num_wann
!      WRITE(stdout, '(200(2f8.4,2x))') (REAL(deltaH(iwann,jwann)),AIMAG(deltaH(iwann,jwann)), jwann=1,num_wann)
    DO iwann = num_wann_occ+1, num_wann_occ+4
      WRITE(stdout, '(200(2f20.15,2x))') (REAL(deltaH(iwann,jwann)),AIMAG(deltaH(iwann,jwann)), jwann=num_wann_occ+1,num_wann_occ+4)
    ENDDO
  ENDIF 
  !
  IF (ionode) THEN 
    WRITE(stdout, '(/, "KI Hamiltonian at k = ", i4)') ik
!    DO iwann = 1, num_wann
!      WRITE(stdout, '(200(2f8.4,2x))') (REAL(Hamlt(ik, iwann,jwann)),AIMAG(Hamlt(ik,iwann,jwann)), jwann=1,num_wann)
    DO iwann = num_wann_occ+1, num_wann_occ+4
      WRITE(stdout, '(200(2f8.4,2x))') (REAL(Hamlt(ik,iwann,jwann)),AIMAG(deltaH(iwann,jwann)), jwann=num_wann_occ+1,num_wann_occ+4)
    ENDDO
  ENDIF 
  !
  WRITE( stdout, 9020 ) ( xk(i,ik), i = 1, 3 )
  WRITE( stdout, '(10x, "KS  ",8F11.4)' ) (eigvl(iwann)*rytoev, iwann=1,num_wann)
  ham(:,:) = Hamlt(ik,:,:) 
  CALL cdiagh( num_wann, ham, num_wann, eigvl, eigvc )
  WRITE( stdout, '(10x, "KI  ",8F11.4)' ) (eigvl(iwann)*rytoev, iwann=1,num_wann)
  !
  !deltaHk(1:num_wann,1:num_wann,ik) = deltaH(1:num_wann,1:num_wan)
  !! Store the result in a global varible
  !
  WRITE(stdout,'(/)') 
  DEALLOCATE (deltaH)
  !
  !
9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
END subroutine 

