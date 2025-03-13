!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.D0,0.D0)
#define ONE (0.D0,1.D0)
!#define DEBUG
!-----------------------------------------------------------------------
SUBROUTINE alpha_corr ( iwann, delta)
  !-----------------------------------------------------------------------
  !
  !! This routine compute a correction to the second order approximation
  !! It is based on the approximation that the difference between the full energy 
  !! at N /pm 1 and the 2nd order approximation is the same for the un-relaxed
  !! and relaxed picture: 
  !! E_u(N-1) = E(N) - e_nn + 0.5*[d^2E/df_n^2]_u + delta_u
  !! E_r(N-1) = E(N) - e_nn + 0.5*[d^2E/df_n^2]_r + delta_r
  !!
  !! Assumption delta_r= delta_u
  !! alpha = (E_r(N-1)-E(N) +e_nn )/(E_u(N-1)=E(N) + e_nn) =
  !!       = (0.5*[d^2E/df_n^2]_r + delta_r) / (0.5*[d^2E/df_n^2]_u + delta_u) ~
  !!       ~ (0.5*[d^2E/df_n^2]_r + delta_u) / (0.5*[d^2E/df_n^2]_u + delta_u) -->
  !! 
  !!       ~ (0.5*[d^2E/df_n^2]_r + delta_u) / (0.5*[d^2E/df_n^2]_u)
  !!
  !! NOTA BENE: the only relevant quantity is the numerator. If we still want to use 
  !!            the 2nd order approx to the Hamiltonian in kcw_ham.x we do not add delta_u
  !!            to the denominator so that the product of the bare kcw hamiltonian at 2nd
  !!            order times alpha give exactly the desired quantity, i.e. the numerator, i.e. 
  !!            the (approximated) DeltaSCF energy fully relaxed.
  !! 
  !! Only delta_u needs to be calculated.
  !! Since everithing is frozen and we look for terms beyond the 2nd order, only the xc controbutions remains. 
  !! delta_u = E_u(N-1) - E(N) + e_nn - 0.5*[d^2E/df_n^2]_u 
  !!         = E^{xc}_u(N-1) - E^{xc}(N) + < phi_0n |v_xc | phi_0n > - 0.5*< rho_0n | f_xc[rho] | rho_0n >
  !! 
  !! The subroutine computes these three quantity:
  !! E^{xc}(N)
  !! E^{xc}_u(N-1)
  !! < phi_0n |v_xc | phi_0n >
  !! < rho_0n | f_xc[rho] | rho_0n >
  !! where phi_0n(r) is the Wannier function and rho_0n its density. 
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE control_kcw,          ONLY : num_wann_occ, nqstot
  USE mp_pools,             ONLY : inter_pool_comm
  !USE funct,                ONLY : dft_is_gradient
  USE xc_lib,               ONLY : xclib_dft_is
  ! 
  IMPLICIT NONE 
  ! 
  INTEGER, INTENT (IN) :: iwann 
  REAL(DP), INTENT (OUT) :: delta 
  LOGICAL :: is_emp
  REAL(DP) :: en, en_pm1, eig, krnl
  ! 
  en     = 0.D0 
  en_pm1 = 0.D0 
  eig    = 0.D0
  krnl   = 0.d0
  ! 
  is_emp = .false.
  IF (iwann .gt. num_wann_occ) is_emp = .true.
  !
  IF ( xclib_dft_is('gradient') .AND. nqstot .gt. 1) & 
     CALL infomsg( 'alpha_corr', 'Plus/minus 1 contribution to the gradient correction DISREGARDED' )
  !
  CALL xc_energy_n ( iwann, en, eig, krnl ) 
  CALL xc_energy_npm1 ( iwann, en_pm1, is_emp )
  !
  delta = en_pm1 -en + eig - 0.5D0*krnl
  IF (is_emp) delta = en_pm1 - (en + eig + 0.5D0*krnl)
  !
  WRITE(stdout,'(5X,"INFO: iwann , e(N), de/df, d2e/df2, e(N-1), delta", i5, 5f20.12)') iwann, en, eig, krnl, en_pm1, delta 
  ! 
  RETURN
  ! 
  CONTAINS 
    !
    !------------------------------------------------------------------------
    SUBROUTINE xc_energy_n ( iwann, en, eig, krnl )
      !----------------------------------------------------------------------
      !
      USE scf,                   ONLY : rho, rho_core, rhog_core
      USE fft_base,              ONLY : dffts
      USE fft_interfaces,        ONLY : fwfft
      USE fft_wave,              ONLY : invfft_wave
      USE lsda_mod,              ONLY : nspin, isk, lsda
      USE klist,                 ONLY : nkstot, ngk, igk_k, nks
      USE control_kcw,           ONLY : evc0, iuwfc_wann, num_wann, spin_component, nqstot, iurho_wann,&
                                        nkstot_eff, nrho
      USE wvfct,                 ONLY : npwx
      USE cell_base,             ONLY : omega
      USE mp_bands,              ONLY : intra_bgrp_comm
      USE buffers,               ONLY : get_buffer
      USE mp,                    ONLY : mp_sum
      USE eqv,                   ONLY : dmuxc
      USE gvecs,                 ONLY : ngms
      USE noncollin_module,      ONLY : domag, noncolin, m_loc, angle1, angle2, &
                                        ux, nspin_lsda, nspin_gga, nspin_mag, npol
      USE dv_of_drho_lr,         ONLY : dv_of_drho_xc
      !
      !
      IMPLICIT NONE 
      !  
      INTEGER, INTENT (IN) :: iwann
      INTEGER :: ik, npw, lrrho, iq, ir, is, ip, iss
      REAL(DP), INTENT (OUT) :: en, eig, krnl 
      REAL(DP) ::  vtxc, etxc, vxc(dffts%nnr,nspin_mag), eig_k, krnl_q
      INTEGER :: lrwfc
      COMPLEX(DP) :: evc_g (npwx*npol), evc_r (dffts%nnr,npol)
      COMPLEX(DP) :: rhowann(dffts%nnr, num_wann, nrho), rhor(dffts%nnr,nrho), delta_vr(dffts%nnr,nspin_mag)
      COMPLEX(DP), ALLOCATABLE  :: rho_wann_g(:,:), delta_vg(:,:), aux(:)
      INTEGER :: ik_eff
      COMPLEX(DP), ALLOCATABLE :: rhor_(:,:)
      !
      ALLOCATE ( rho_wann_g (ngms,nrho) , delta_vg(ngms,nspin_mag), aux (dffts%nnr) )
      ! 
      ! The Exc(N) energy: just the xc energy in the PC times the SC dimension
      !
      etxc = 0.D0; vtxc = 0.D0; vxc(:,:) = 0.D0
      CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, vxc )
      en = etxc*nkstot_eff
      !
      ! The xc contribution to the expertation value of the DFT Ham on the Wann function
      ! eig = \sum_k <u_nk | vxc[\rho] u_nk > (A part from prefactors that needs to be checked 
      ! throughout the code .... )
      eig  = 0.D0
      krnl = 0.D0
      ! 
      DO ik = 1, nks
        !
        IF ( lsda .AND. isk(ik) /= spin_component) CYCLE
        npw = ngk(ik)
        lrwfc = num_wann * npwx * npol
        ik_eff = ik-(spin_component-1)*nkstot_eff
        CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik_eff )
        evc_g(:) =  evc0(:,iwann)
        !CALL get_buffer ( evc, nwordwfc, iuwfc, ik )
        !evc_g(:) =  evc(:,iwann)
        !
        evc_r(:,:) = ZERO
        CALL invfft_wave (npwx, npw, igk_k (1,ik), evc_g , evc_r )
        !! The wfc in R-space at k
        IF (nspin_mag==2) THEN
            eig_k = sum ( vxc(:,spin_component) * evc_r(:,1) * CONJG(evc_r(:,1) ) ) !check the (:,->1)
        ELSE
            eig_k = sum ( vxc(:,1) * (evc_r(:,1) * CONJG(evc_r(:,1) )+evc_r(:,2) * CONJG(evc_r(:,2)))) + & 
                    sum ( vxc(:,2) * (conjg(evc_r(:,1))*evc_r(:,2) + conjg(evc_r(:,2))*evc_r(:,1))) + & 
                    sum ( vxc(:,3) * (-CMPLX(0.D0,1.D0, kind=DP) * conjg(evc_r(:,1))*evc_r(:,2) & 
                                     + CMPLX(0.D0,1.D0, kind=DP) * conjg(evc_r(:,2))*evc_r(:,1))) + & 
                    sum ( vxc(:,4) * ( conjg(evc_r(:,1))*evc_r(:,1) - conjg(evc_r(:,2))*evc_r(:,2)))
        END IF
        eig_k = eig_k/( dffts%nr1*dffts%nr2*dffts%nr3 )
        CALL mp_sum (eig_k, intra_bgrp_comm) 
        eig = eig + eig_k
        !
      ENDDO
      CALL mp_sum (eig, inter_pool_comm )
      eig = eig/(nkstot_eff)
      !
      ! The kernel term (as in bare_pot.f90, but only xc contribution)
      !
      DO iq = 1, nqstot
        !
        lrrho=num_wann * dffts%nnr * nrho
        CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
        rhor(:,:) = rhowann(:,iwann,:)
        !
        DO ip=1,nrho
          aux(:) = rhor(:,ip)/omega
          CALL fwfft ('Rho', aux, dffts)  ! NsC: Dense or smooth grid?? I think smooth is the right one. 
          rho_wann_g(:,ip) = aux(dffts%nl(:))
        END DO
        !OLD IMPLEMENTATION FOR NSPIN==2
        !aux(:) = rhor(:)/omega
        !CALL fwfft ('Rho', aux, dffts)  ! NsC: Dense or smooth grid?? I think smooth is the right one. 
        !rho_wann_g(:) = aux(dffts%nl(:))
        !
        delta_vr = CMPLX(0.D0,0.D0, kind=DP)
        !
        ALLOCATE ( rhor_(dffts%nnr,nspin_mag) )
        rhor_ = CMPLX(0.D0,0.D0,kind=DP)
        IF (nspin_mag == 4) THEN
          rhor_=rhor/omega
        ELSE
          rhor_(:,spin_component) = rhor(:,1)/omega
        ENDIF
        CALL dv_of_drho_xc(delta_vr, rhor_)
        DEALLOCATE (rhor_)
        !
!        IF (nspin_mag==2) THEN
!          DO is = 1, nspin_mag
!            DO ir = 1, dffts%nnr
!              delta_vr(ir,is) = delta_vr(ir,is) + dmuxc(ir,is,spin_component) * rhor(ir,1)/omega
!            END DO
!          END DO
!        ELSE
!          DO is = 1, nspin_mag
!              DO ir = 1, dffts%nnr
!                  DO iss = 1, nspin_mag
!                delta_vr(ir,is) = delta_vr(ir,is) + dmuxc(ir,is,iss) * rhor(ir,iss)/omega
!              END DO
!            END DO
!          END DO
!        END IF
        ! In g-space
        DO is = 1, nspin_mag
          aux(:) = delta_vr(:,is)
          CALL fwfft ('Rho', aux, dffts)
          delta_vg(:,is) = aux(dffts%nl(:))
        ENDDO
        IF (nspin_mag==2) THEN
          krnl_q = sum (CONJG(rho_wann_g (:,1)) * delta_vg(:,spin_component))*omega
        ELSE
          krnl_q = 0.0_DP
          DO ip=1,nrho
            krnl_q = krnl_q + sum (CONJG(rho_wann_g (:,ip)) * delta_vg(:,ip))*omega
          END DO
        END IF
        CALL mp_sum (krnl_q, intra_bgrp_comm)
        krnl = krnl + krnl_q/nqstot
      ENDDO
      !
      DEALLOCATE ( rho_wann_g , delta_vg, aux )
      RETURN
      !
    END SUBROUTINE xc_energy_n
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE xc_energy_npm1 ( iwann, en_pm1, is_emp )
      !----------------------------------------------------------------------
      !
      USE scf,                   ONLY : rho, rho_core, rhog_core, create_scf_type, scf_type
      USE fft_base,              ONLY : dffts, dfftp
      USE fft_interfaces,        ONLY : fwfft
      USE lsda_mod,              ONLY : nspin
      USE klist,                 ONLY : nkstot
      USE control_kcw,           ONLY : num_wann, nqstot, iurho_wann, &
                                        rvect, x_q, nkstot_eff, nrho
      USE cell_base,             ONLY : omega
      USE buffers,               ONLY : get_buffer
      USE mp,                    ONLY : mp_sum
      USE constants,             ONLY : tpi
      !USE funct,                 ONLY : dft_is_gradient
      USE xc_lib,                ONLY : xclib_dft_is
      USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux, nspin_lsda, nspin_gga, nspin_mag, npol
      !
      IMPLICIT NONE 
      !
      INTEGER, INTENT (IN) :: iwann
      REAL(DP), INTENT (OUT) :: en_pm1
      LOGICAL, INTENT (IN) :: is_emp
      !
      INTEGER :: ir, iq, nspin_aux, segno
      TYPE (scf_type) :: rho_minus1
      COMPLEX (DP) :: rho_wann_ir(dffts%nnr,nrho), rho_wann(dffts%nnr,nrho)
      REAL(DP) :: xq(3), xq_(3)
      COMPLEX(DP) :: rhowann(dffts%nnr, num_wann,nrho)
      COMPLEX(DP) :: phase_sc, phase_pc(dffts%nnr) 
      COMPLEX(DP), ALLOCATABLE :: nq_r(:,:), aux(:)
      INTEGER :: lrrho
      LOGICAL :: lgamma
      REAL(DP) ::  vtxc, etxc, vxc(dffts%nnr,nspin_mag)
      !
      COMPLEX(DP) :: imag = (0.D0,1.D0)
      !
      ALLOCATE ( nq_r (dffts%nnr,nrho), aux(dfftp%nnr)  )
      segno=-1
      IF ( is_emp ) segno=+1
      !
      nspin_aux=nspin_mag
      IF (nspin_mag .ne. 4) THEN
        nspin_mag=2
      CALL create_scf_type (rho_minus1)
      ELSE
      CALL create_scf_type (rho_minus1)
      END IF
      nspin_mag=nspin_aux
      !
      rho_wann = CMPLX(0.D0, 0.D0, kind=DP)
      en_pm1 = 0.D0 
      !
      DO ir = 1, nkstot_eff
        !
        ! Initialize the spin-component for the N \pm 1 density
        IF (nspin_mag == 1 ) THEN
          rho_minus1%of_r(:,1) = rho%of_r(:,1)
          rho_minus1%of_r(:,2) = 0.D0
          rho_minus1%of_g(:,1) = rho%of_g(:,1)
          rho_minus1%of_g(:,2) = 0.D0
        ELSE IF (nspin_mag == 2) THEN
          rho_minus1%of_r(:,:) = rho%of_r(:,:)
          rho_minus1%of_g(:,:) = rho%of_g(:,:)
        ELSE
          rho_minus1%of_r(:,1) = rho%of_r(:,1)
          rho_minus1%of_r(:,2) = rho%of_r(:,2)
          rho_minus1%of_r(:,3) = rho%of_r(:,3)
          rho_minus1%of_r(:,4) = rho%of_r(:,4)
          rho_minus1%of_g(:,1) = rho%of_g(:,1)
          rho_minus1%of_g(:,2) = rho%of_g(:,2)
          rho_minus1%of_g(:,3) = rho%of_g(:,3)
          rho_minus1%of_g(:,4) = rho%of_g(:,4)                              
        ENDIF
#ifdef DEBUG
        WRITE(*,'("Int rho Im, Re", 1f15.8, 4i8, f15.6)') SUM(ABS(rho%of_r(:,1))), dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nnr, omega
        WRITE(*,'("Int rho Im, Re", 1f15.8, 4i8, f15.6)') SUM(ABS(rho%of_r(:,1)))*omega/dfftp%nnr
        WRITE(*,'("Int rho Im, Re", 1f15.8, 4i8, f15.6)') SUM(ABS(rho_minus1%of_r(:,1))), dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nnr, omega
        WRITE(*,'("Int rho Im, Re", 1f15.8, 4i8, f15.6)') SUM(ABS(rho_minus1%of_r(:,1)))*omega/dfftp%nnr
        WRITE(*,'("Int rho Im, Re", 1f15.8, 4i8, f15.6)') SUM(ABS(rho_minus1%of_r(:,2))), dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nnr, omega
        WRITE(*,'("Int rho Im, Re", 1f15.8, 4i8, f15.6)') SUM(ABS(rho_minus1%of_r(:,2)))*omega/dfftp%nnr
        WRITE(*,'("Re n_r",5f15.8)') REAL(rho%of_r(1:5,1))
        WRITE(*,'("Re n-1_r UP",5f15.8)') REAL(rho_minus1%of_r(1:5,1))
        WRITE(*,'("Re n-1_r DW",5f15.8)') REAL(rho_minus1%of_r(1:5,2))
        WRITE(*,'("Re n-1_g UP",5f15.8)') REAL(rho_minus1%of_g(1:5,1))
        WRITE(*,'("Re n-1_g DW",5f15.8)') REAL(rho_minus1%of_g(1:5,2))
        WRITE(*,'("Im n-1_g UP",5f15.8)') AIMAG(rho_minus1%of_g(1:5,1))
        WRITE(*,'("Im n-1_g DW",5f15.8)') AIMAG(rho_minus1%of_g(1:5,2))
        !
        WRITE(*,'("rvect",i3, 3f12.5)')  ir, rvect(:,ir) 
#endif
        rho_wann_ir = CMPLX(0.D0, 0.D0, kind=DP)
        DO iq = 1, nqstot
           !
           xq(1:3)  = x_q(1:3,iq)
           lgamma = ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 )
#ifdef DEBUG
           WRITE(*,'("q",i3, 3f12.5, 3x, l)') iq,  xq(:), lgamma
#endif
           !
           !CALL cryst_to_cart( 1, xq, at, -1 )
           !
           ! the peridoc part of the wannier density for this q
           lrrho=num_wann*dffts%nnr*nrho
           CALL get_buffer (rhowann, lrrho, iurho_wann, iq)
           nq_r(:,:) = rhowann(:,iwann,:)/omega
           ! plus or minus depending whether the wannier is empty or full
           nq_r(:,:) = segno * nq_r(:,:)
           !
#ifdef DEBUG
           WRITE(*,'("Re nq_r",5f15.8)') REAL(nq_r(1:5))
           WRITE(*,'("Im nq_r",5f15.8)') AIMAG(nq_r(1:5))
#endif
           ! exp(-iq*R)
           phase_sc = EXP( - imag * tpi * DOT_PRODUCT( xq, rvect(:,ir) ) ) 
           ! exp (iq*r), there is a minus sign inside calculate_phase
           xq_=-x_q(:,iq)
           phase_pc=CMPLX(0.D0,0.D0, kind=DP)
           CALL calculate_phase (xq_, phase_pc)
           !
#ifdef DEBUG
           WRITE(*,'("PHASES", 2f12.5, 3x, 10f12.5)') phase_sc, phase_pc(1:5)
#endif
           ! ... accumulate over q points ...
           rho_wann_ir(:,1) = rho_wann_ir(:,1) + phase_sc*phase_pc(:)*nq_r(:,1)/nqstot
           if (nspin_mag==4) then
            rho_wann_ir(:,2) = rho_wann_ir(:,2) + phase_sc*phase_pc(:)*nq_r(:,2)/nqstot
            rho_wann_ir(:,3) = rho_wann_ir(:,3) + phase_sc*phase_pc(:)*nq_r(:,3)/nqstot
            rho_wann_ir(:,4) = rho_wann_ir(:,4) + phase_sc*phase_pc(:)*nq_r(:,4)/nqstot
           end if
        !
        ENDDO ! q-points LOOP
        !
        ! ... accumulate over R. To compute the integral of the Wannier over the SC ...
        rho_wann (:,:) = rho_wann (:,:) + rho_wann_ir(:,:)  ! DEBUG
        !
#ifdef DEBUG
        WRITE(*,'("Re rhowann(R)", 5f15.8)') REAL(rho_wann_ir(1:5))
        WRITE(*,'("Im rhowann(R)", 5f15.8)') AIMAG(rho_wann_ir(1:5))
        WRITE(*,'("Int rhowann(R) Im, Re", 2f15.8, 4i8, f15.6)') SUM(ABS(AIMAG(rho_wann_ir(:)))), SUM(REAL(rho_wann_ir(:))), dffts%nr1, dffts%nr2, dffts%nr3, dffts%nnr, omega
        WRITE(*,'("Int rhowann(R) Im, Re", 2f15.8, 4i8, f15.6)') SUM(ABS(AIMAG(rho_wann_ir(:))))*omega/dffts%nnr, SUM(REAL(rho_wann_ir(:)))*omega/dffts%nnr
#endif
        !
        ! ... Add the wannier density to the GS density... 
        IF (nspin_mag==4) THEN
        rho_minus1%of_r(:,1) = REAL(rho_wann_ir(:,1)) + rho_minus1%of_r(:,1)
        rho_minus1%of_r(:,2) = REAL(rho_wann_ir(:,2)) + rho_minus1%of_r(:,2)
        rho_minus1%of_r(:,3) = REAL(rho_wann_ir(:,3)) + rho_minus1%of_r(:,3)
        rho_minus1%of_r(:,4) = REAL(rho_wann_ir(:,4)) + rho_minus1%of_r(:,4)
        ELSE
        rho_minus1%of_r(:,1) = REAL(rho_wann_ir(:,1)) + rho_minus1%of_r(:,1)
        rho_minus1%of_r(:,2) = REAL(rho_wann_ir(:,1)) + rho_minus1%of_r(:,2) !@Nicola Colonna please check
        END IF
        !
        ! NOTA BENE: The following is correct only in a 1 k-point calculation (Supercell) 
        ! This density is periodic in the SC and this FFT does not make sense if nq/=1 
        ! The rho in  G-space is needed to compute the Gradient correction (if present) 
        !IF ( .NOT. dft_is_gradient ( ) .OR. nqstot == 1 ) THEN 
        ! NB: THE CODE HERE BELOW HAS NOT EXTENDED TO NON-COLLINEAR MODE
        IF ( .NOT. xclib_dft_is('gradient') .OR. nqstot == 1 ) THEN
          !
          aux(:) = rho_minus1%of_r(:,1)
          CALL fwfft ('Rho', aux, dfftp)
          rho_minus1%of_g(:,1) = aux(dfftp%nl(:))
          !
          aux(:) = rho_minus1%of_r(:,2)
          CALL fwfft ('Rho', aux, dfftp)
          rho_minus1%of_g(:,2) = aux(dfftp%nl(:))
        ENDIF
        !
#ifdef DEBUG
        WRITE(*,'("Re n-1_r UP",5f15.8)') REAL(rho_minus1%of_r(1:5,1))
        WRITE(*,'("Re n-1_r DW",5f15.8)') REAL(rho_minus1%of_r(1:5,2))
        WRITE(*,'("Re n-1_g UP",5f15.8)') REAL(rho_minus1%of_g(1:5,1))
        WRITE(*,'("Re n-1_g DW",5f15.8)') REAL(rho_minus1%of_g(1:5,2))
        WRITE(*,'("Im n-1_g UP",5f15.8)') AIMAG(rho_minus1%of_g(1:5,1))
        WRITE(*,'("Im n-1_g DW",5f15.8)') AIMAG(rho_minus1%of_g(1:5,2))
#endif
        !
        IF (nspin_mag .ne. 4) THEN
        nspin_aux=nspin_mag; nspin_mag=2
        etxc = 0.D0; vtxc = 0.D0; vxc(:,:) = 0.D0
        CALL v_xc( rho_minus1, rho_core, rhog_core, etxc, vtxc, vxc )
        nspin_mag=nspin_aux
        ELSE
        etxc = 0.D0; vtxc = 0.D0; vxc(:,:) = 0.D0
        CALL v_xc( rho_minus1, rho_core, rhog_core, etxc, vtxc, vxc )
        END IF
        !
        en_pm1 = en_pm1 + etxc
        !
      ENDDO
#ifdef DEBUG
      WRITE(*,'("Int rhowann Im, Re", 2f15.8, 4i8, f15.6)') SUM(ABS(AIMAG(rho_wann(:)))), SUM(REAL(rho_wann(:))), dffts%nr1, dffts%nr2, dffts%nr3, dffts%nnr, omega
      WRITE(*,'("Int rhowann Im, Re", 2f15.8, 4i8, f15.6)') SUM(ABS(AIMAG(rho_wann(:)))), SUM(REAL(rho_wann(:)))*omega/dffts%nnr
#endif
      !
      DEALLOCATE ( nq_r, aux  )
      RETURN
      !
    END SUBROUTINE xc_energy_npm1
    !
END SUBROUTINE alpha_corr
