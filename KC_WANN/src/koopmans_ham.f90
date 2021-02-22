! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#define DEBUG
!-----------------------------------------------------------------------
SUBROUTINE koopmans_ham ()
  !---------------------------------------------------------------------
  !
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : DP
  USE klist,                 ONLY : nkstot, xk
  USE lsda_mod,              ONLY : nspin
  USE control_kc_wann,       ONLY : num_wann, Hamlt, nqstot, l_alpha_corr,&
                                    alpha_final, num_wann_occ, on_site_only
  USE constants,             ONLY : rytoev
  USE wvfct,                 ONLY : npwx
  !
  IMPLICIT NONE
  !
  ! The k point index 
  INTEGER :: ik
  !
  ! The scalar part (independent on k) <rho_0i|v_0i|rho_0i>delta_ij
  COMPLEX(DP) :: deltah_scal (num_wann, num_wann)
  !
  ! The Hamiltonian due to the real contribution to the potential \int f_Hxc(r,r') rho_0i(r')
  COMPLEX(DP) deltah_real (num_wann, num_wann)
  !
  ! the KI hamiltonian, the KI contribution, and the new eigenvectors at a given k-point
  COMPLEX(DP) :: ham(num_wann,num_wann), deltaH(num_wann,num_wann), eigvc(npwx,num_wann)
  !
  ! The new eigenalues 
  REAL(DP) :: eigvl(num_wann)
  !
  ! The correction to the diagonal term beyond second order
  REAL(DP) :: ddH(num_wann)
  !
  INTEGER :: i, iwann, jwann
  !
  !
  IF (on_site_only) WRITE(stdout, '(/,5X, "INFO: Skipping off-diag: only R=0 and i=j")') 
  !
  ! The scalar term 
  deltah_scal=CMPLX(0.D0,0.D0,kind=DP)
  CALL ham_scalar (deltah_scal)
  ! 
#ifdef DEBUG
  WRITE( stdout,'(/,5X," Scalar term Hamiltonian:")')
  DO iwann=1, num_wann
    WRITE(stdout,'(5X,10(2F10.6, 2x))') (deltah_scal(iwann,jwann), jwann=1,num_wann)
  ENDDO
  ! DEBUG
#endif
  !
  ! ... The correction beyond 2nd order 
  IF (l_alpha_corr) CALL beyond_2nd (deltah_scal, ddH)
  !
#ifdef DEBUG
  WRITE( stdout,'(/,5X," Scalar term Hamiltonian plus correction:")')
  ham = deltah_scal
  DO iwann=1, num_wann
    ham(iwann,iwann) = ham(iwann,iwann)-ddH(iwann)
    WRITE(stdout,'(5X,10(2F10.6, 2x))') (ham(iwann,jwann), jwann=1,num_wann)
  ENDDO
  ! DEBUG
#endif
  !
  DO ik = 1, nkstot/nspin
    !
    deltah_real = CMPLX(0.D0, 0.D0, kind = DP)
    !
    IF (.NOT. on_site_only) THEN 
       ! General routine for empty state hamiltonian at k
       CALL koopmans_ham_real_k ( ik, deltah_real )
    ELSE
      ! SKIP off-diagonal elements in REAL SPACE (i.e. R/=0 i/=j) 
      ! If only R=0 and i=j, then the "real" contribution for empty states
      ! is <wi| {\int f_Hxc wi^2} |wi> (which is simply -2.D0 times
      ! the scalar controbution -0.5*<wi^2|f_Hxc|wi^2>) 
      DO iwann = num_wann_occ+1, num_wann
        deltah_real(iwann,iwann) = -2.D0*deltah_scal(iwann,iwann)
      ENDDO
    ENDIF
    !
    !
#ifdef DEBUG
    WRITE( stdout,'(/,5X," Real term Hamiltonian ik =:", i5)') ik
    DO iwann=1, num_wann
      WRITE(stdout,'(5X,10(2F10.6, 2x))') (deltah_real(iwann,jwann), jwann=1,num_wann)
    ENDDO
#endif
    !
    ! The KS hamiltonian in the Wannier Gauge (just to check)
    ham(:,:) = Hamlt(ik,:,:) 
    CALL cdiagh( num_wann, ham, num_wann, eigvl, eigvc )
    !
    ! The KI contribution at k
    deltaH = deltah_scal + deltah_real
    !
    ! Apply the screening coefficient
!    DO iwann = 1, num_wann; deltaH(:,iwann) = alpha_final (iwann)* deltaH(:,iwann) ; ENDDO  
    DO iwann = 1, num_wann
      DO jwann = iwann , num_wann
        deltaH(iwann,jwann) = deltaH(iwann,jwann)*alpha_final(jwann)
        deltaH(jwann,iwann) = CONJG(deltaH(iwann,jwann))
      ENDDO
    ENDDO
  
    !
    ! Add to the KS Hamiltonian
    Hamlt(ik,:,:) = Hamlt(ik,:,:) + deltaH(:,:) 
    !
#ifdef DEBUG
    WRITE(stdout, '("KI Contribution to the Hamiltonian at k = ", i4)') ik
    DO iwann = 1, num_wann
      WRITE(stdout, '(200(2f8.4,2x))') (REAL(deltaH(iwann,jwann)),AIMAG(deltaH(iwann,jwann)), jwann=1,num_wann)
    ENDDO
    !
    WRITE(stdout, '(/, "KI Hamiltonian at k = ", i4)') ik
    DO iwann = 1, num_wann
      WRITE(stdout, '(200(2f8.4,2x))') (REAL(Hamlt(ik, iwann,jwann)),AIMAG(Hamlt(ik,iwann,jwann)), jwann=1,num_wann)
    ENDDO
#endif
    !
    WRITE( stdout, 9020 ) ( xk(i,ik), i = 1, 3 )
    WRITE( stdout, '(10x, "KS  ",8F9.4)' ) (eigvl(iwann)*rytoev, iwann=1,num_wann)
    !
    ham(:,:) = Hamlt(ik,:,:) 
    CALL cdiagh( num_wann, ham, num_wann, eigvl, eigvc )
    WRITE( stdout, '(10x, "KI  ",8F9.4)' ) (eigvl(iwann)*rytoev, iwann=1,num_wann)
    !
  ENDDO
  !
  ! For an isolated molecule the full Hamiltonian is also available
  ik=1
  IF (nkstot/nspin == 1) CALL full_ham ( ik )
  !
  !
  9020 FORMAT(/'          k =',3F7.4,'     band energies (ev):'/ )
  RETURN
  !
  CONTAINS 
  !
  !
  !----------------------------------------------------------------
  SUBROUTINE beyond_2nd (deltaH, ddH)
    !----------------------------------------------------------------
    !
    USE kinds,                 ONLY : DP
    USE control_kc_wann,       ONLY : num_wann, num_wann_occ, &
                                      alpha_final, alpha_final_full
    USE control_lr,            ONLY : lrpa
    !
    IMPLICIT NONE  
    !
    COMPLEX(DP), INTENT(INOUT) :: deltaH(num_wann,num_wann)
    REAL(DP), INTENT(OUT) :: ddH(num_wann)
    !
    REAL(DP) :: alpha_(num_wann), alpha_fd
    ! weight of the q points, alpha from LR, alpha after the all-order correction. 
    !
    REAL(DP) second_der(num_wann), delta 
    !
    INTEGER :: iwann
    !
    ddH = 0.D0
    !
    WRITE( stdout, '(/,5X, "INFO: Correction beyond 2nd order ...",/)')
    IF (lrpa) THEN 
      !
      WRITE(*, '(8X,"INFO: l_alpha_corr and lrpa are NOT consistent.At RPA")')
      WRITE(*, '(8X,"      level there is no contribution beyond 2nd order.")')
      WRITE(*, '(8X,"      Nothing to do here. RETURN")')
      !
      RETURN
      !
    ENDIF 
    !
    DO iwann = 1, num_wann
      second_der(iwann) = -REAL(deltaH(iwann,iwann))
    ENDDO
    !
    !DO iwann = 1, num_wann_occ
    DO iwann = 1, num_wann
      delta =0.D0 
      alpha_(iwann) = alpha_final(iwann) 
      !
      !
       ! ... Compute the difference between the parabolic extrapolation at N \pm 1 and the real 
       ! ... value of the energy in the frozen orbital approximation ...
       CALL alpha_corr (iwann, delta)
       ddH(iwann) = delta
       !deltaH(iwann,iwann) = deltaH(iwann,iwann)-ddH(iwann)
       !
       ! ... The new alpha that should be closer to the Finite-difference one ...
       ! ... Remember DeltaH is nothing but the second derivative wrt the orbital occupation ...
       alpha_fd = (alpha_final(iwann)*second_der(iwann) + delta)/ (second_der(iwann)+delta)
       IF(nkstot/nspin == 1) alpha_final_full(iwann) = alpha_fd
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
      !
    ENDDO
    !
  END SUBROUTINE beyond_2nd
  !
  !----------------------------------------------------------------
  SUBROUTINE ham_scalar (deltah_scal)
    !----------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout
    USE control_kc_wann,      ONLY : num_wann, nqstot, iurho_wann, &
                                     spin_component
    USE fft_base,             ONLY : dffts
    USE cell_base,            ONLY : omega
    USE gvecs,                ONLY : ngms
    USE mp_bands,             ONLY : intra_bgrp_comm
    USE mp,                   ONLY : mp_sum
    USE buffers,              ONLY : get_buffer
    !
    IMPLICIT NONE
    !
    ! The scalar contribution to the hamiltonian 
    COMPLEX(DP), INTENT(INOUT) :: deltah_scal (num_wann, num_wann)
    !
    ! Couters for the q point, wannier index. record length for the wannier density
    INTEGER :: iq, iwann, lrrho
    !
    ! The periodic part of the wannier orbital density
    COMPLEX(DP) :: rhowann(dffts%nnr, num_wann), rhor(dffts%nnr), delta_vr(dffts%nnr,nspin), delta_vr_(dffts%nnr,nspin)
    !
    ! The self Hartree
    COMPLEX(DP) :: sh(num_wann)
    !
    ! Auxiliary variables 
    COMPLEX(DP), ALLOCATABLE  :: rhog(:), delta_vg(:,:), vh_rhog(:), delta_vg_(:,:)
    !
    ! The weight of each q point
    REAL(DP) :: weight(nqstot)
    !
    WRITE( stdout, '(/,5X, "INFO: KC SCALAR TERM CALCULATION ... START")')
    !
    ALLOCATE ( rhog (ngms) , delta_vg(ngms,nspin), vh_rhog(ngms), delta_vg_(ngms,nspin) )
    !
    DO iq = 1, nqstot
      !
      lrrho=num_wann*dffts%nnr
      CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
      !! Retrive the rho_wann_q(r) from buffer in REAL space
      !
      weight(iq) = 1.D0/nqstot ! No SYMM 
      !
      DO iwann = 1, num_wann  ! for each band, that is actually the perturbation
         !
         rhog(:)         = CMPLX(0.D0,0.D0,kind=DP)
         delta_vg(:,:)   = CMPLX(0.D0,0.D0,kind=DP)
         vh_rhog(:)      = CMPLX(0.D0,0.D0,kind=DP)
         rhor(:)         = CMPLX(0.D0,0.D0,kind=DP)
         !
         rhor(:) = rhowann(:,iwann)
         !! The periodic part of the orbital desity in real space
         !
         CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
         !! The periodic part of the perturbation DeltaV_q(G)
         ! 
         sh(iwann) = sh(iwann) + 0.5D0 * sum (CONJG(rhog (:)) * vh_rhog(:)                )*weight(iq)*omega
         deltah_scal(iwann, iwann) = deltah_scal(iwann,iwann) - 0.5D0 * sum (CONJG(rhog (:)) * delta_vg(:,spin_component))*weight(iq)*omega
         !
      ENDDO
      ! 
    ENDDO ! qpoints
    WRITE( stdout, '(/,5X, "INFO: KC SCALAR TERM CALCULATION ... END")')
    !
    DEALLOCATE ( rhog , delta_vg, vh_rhog, delta_vg_ )
    !
    CALL mp_sum (deltah_scal, intra_bgrp_comm)
    CALL mp_sum (sh, intra_bgrp_comm)
   !
  END SUBROUTINE ham_scalar
  !
  !-----------------------------------------------------------------------
  SUBROUTINE koopmans_ham_real_k (ik, deltaH)
    !---------------------------------------------------------------------
    !
    ! This routine compute the KI real term correction to second order to the KS
    ! Hamiltonian at a given k point ik given from input
    !
    USE kinds,                ONLY : DP
    USE fft_base,             ONLY : dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE klist,                ONLY : igk_k, ngk
    USE mp,                   ONLY : mp_sum
    USE control_kc_wann,      ONLY : spin_component, num_wann, &
                                     num_wann_occ, evc0, iuwfc_wann, &
                                     map_ikq, shift_1bz
    USE buffers,              ONLY : get_buffer, save_buffer
    USE io_global,            ONLY : stdout
    USE control_kc_wann,      ONLY : iurho_wann
    USE mp_bands,             ONLY : intra_bgrp_comm
    USE dv_of_drho_lr,        ONLY : dv_of_drho
    USE control_lr,           ONLY : lgamma
    USE disp,                 ONLY : x_q
    USE lsda_mod,             ONLY : nspin
    USE gvecs,                ONLY : ngms
    USE solve_linter_koop_mod 
    USE qpoint,               ONLY : xq
    USE wvfct,                ONLY : npwx 
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
    ! The KI real term contribution to the Hamiltonian
    COMPLEX(DP), INTENT (INOUT) :: deltaH(num_wann, num_wann)
    !
    INTEGER :: iq, nqs
    ! Counter for the k/q points in the BZ, total number of q points and number of pw for a given k (k+q) point
    !
    INTEGER :: iwann, jwann, lrrho, lrwfc
    ! Band counters, leght of the rho record
    !
    COMPLEX(DP) :: rhowann(dffts%nnr, num_wann), rhor(dffts%nnr), delta_vr(dffts%nnr,nspin), sh(num_wann), delta_vr_(dffts%nnr,nspin)
    ! The periodic part of the wannier orbital density in r space
    ! The perturbig potential in real space
    ! The self-hartree 
    ! The perturbig potential in real space (without the g=0 contribution) 
    !
    COMPLEX(DP), ALLOCATABLE  :: rhog(:), delta_vg(:,:), vh_rhog(:), delta_vg_(:,:)
    ! The periodic parte of the wannier density in G-space 
    ! The perturbig potential in G space
    ! The hartree potential due to the wannier in G-space 
    ! The perturbig potential in G space without the g=0 contribution 
    !
    REAL(DP) :: weight(nqstot)
    ! weight of the q points, alpha from LR, alpha after the all-order correction. 
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
    ! compute Off-diagonal elements. NsC: not sure off_diag=.false. here makes sense: DO NOT CHANGE!!!!
    !
    WRITE( stdout, '(/,/,5X,"INFO: KC HAMILTONIAN CALCULATION ik= ", i4, " ...", /)') ik
    !
    nqs = nqstot
    !
    !
    deltaH = CMPLX(0.D0,0.D0,kind=DP)
    sh     = CMPLX(0.D0,0.D0,kind=DP)
    !
    lrwfc = num_wann*npwx
    CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik )
    ! Retrive the ks function at k (in the Wannier Gauge)
    ! IF (kc_iverbosity .gt. 0 ) WRITE(stdout,'(8X, "INFO: u_k(g) RETRIEVED"/)') 
    !
    CALL compute_map_ikq_single (ik)
    ! find tha map k+q --> k'+G and store the res 
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
      !IF (kc_iverbosity .gt. 0 ) WRITE(stdout,'(8X, "INFO: rhowan_q(r) RETRIEVED"/)') 
      !
      ALLOCATE ( rhog (ngms) , delta_vg(ngms,nspin), vh_rhog(ngms), delta_vg_(ngms,nspin) )
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
      evc0_kq = CMPLX(0.D0,0.D0,kind=DP)
      lrwfc = num_wann * npwx 
      CALL get_buffer ( evc0_kq, lrwfc, iuwfc_wann, ikq )
      ! Retrive the ks function at k+q (in the Wannier Gauge): 
      !
      !IF (kc_iverbosity .gt. 0 ) WRITE(stdout,'(8X, "INFO: u_kq(g) RETRIEVED")') 
      !
      DO iwann = num_wann_occ+1, num_wann
!      DO iwann = 1, num_wann
         !
         npw_k = ngk(ik)
         evc_k_g(:) =  evc0(:,iwann)
         evc_k_r(:) = CMPLX(0.D0,0.D0,kind=DP)
         CALL invfft_wave (npw_k, igk_k (1,ik), evc_k_g , evc_k_r )
         !! The wfc R=0 n=iwann in R-space at k
         !
         !DO jwann = iwann+1, num_wann 
         DO jwann = iwann, num_wann 
            !
            IF (.NOT. off_diag .AND. jwann /= iwann) CYCLE 
            !
            rhog(:)         = CMPLX(0.D0,0.D0,kind=DP)
            delta_vg(:,:)   = CMPLX(0.D0,0.D0,kind=DP)
            vh_rhog(:)      = CMPLX(0.D0,0.D0,kind=DP)
            rhor(:)         = CMPLX(0.D0,0.D0,kind=DP)
            !
            rhor(:) = rhowann(:,jwann)
            ! The periodic part of the orbital desity R=0, n=iwann in real space
            !
            CALL bare_pot ( rhor, rhog, vh_rhog, delta_vr, delta_vg, iq, delta_vr_, delta_vg_ )
            ! The periodic part of the perturbation DeltaV_q(G)
            !
            npw_kq = ngk(ikq)
            evc_kq_g = evc0_kq(:,jwann)
            evc_kq_r = CMPLX(0.D0,0.D0,kind=DP)
            CALL invfft_wave (npw_kq, igk_k (1,ikq), evc_kq_g , evc_kq_r )
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
         ! 
         !
      ENDDO ! iwann
      !
      DEALLOCATE ( rhog , delta_vg, vh_rhog, delta_vg_ )
      !
      !    
    ENDDO ! qpoints
    !WRITE( stdout, '(5X,"INFO: KC HAMILTONIAN CALCULATION ik= ", i4, " ... DONE")') ik
    !
    deltaH = nqstot*deltaH
    CALL mp_sum (deltaH, intra_bgrp_comm)
    CALL mp_sum (sh, intra_bgrp_comm)
    ! Sum over different processes (G vectors) 
    !
    RETURN 
    !
  END subroutine 
  ! 
END SUBROUTINE koopmans_ham
