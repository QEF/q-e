!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#include "f_defs.h"
#define DEBUG
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!-----------------------------------------------------------------------
SUBROUTINE full_ham (ik)
  !-----------------------------------------------------------------------
  !
  ! This routine compute the KI correction to the spectrum
  !
  USE wavefunctions,         ONLY : psic
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout
  USE wvfct,                 ONLY : npw, nbnd, npwx
  USE wvfct,                 ONLY : et
  USE fft_base,              ONLY : dffts, dfftp
  USE lsda_mod,              ONLY : lsda, current_spin,nspin
  USE klist,                 ONLY : nks, ngk, init_igk, igk_k, nkstot
  USE gvect,                 ONLY : ngm
  USE buffers,               ONLY : get_buffer
  USE fft_interfaces,        ONLY : fwfft, invfft
  USE control_kc_wann,       ONLY : kc_at_ks, homo_only, alpha_final_full, hamlt, num_wann_occ, iuwfc_wann, &
                                    kc_iverbosity, qp_symm, evc0, kipz_corr, num_wann, spin_component
  USE control_lr,            ONLY : lrpa
  USE mp,                    ONLY : mp_sum
  USE mp_bands,              ONLY : intra_bgrp_comm
  USE cell_base,             ONLY : omega
  USE scf,                   ONLY : rho, rho_core, rhog_core, create_scf_type, scf_type
  USE constants,             ONLY : rytoev
  !
  USE exx,                   ONLY : vexx, exxinit, exxenergy2
  USE funct,                 ONLY : init_dft_exxrpa, dft_is_hybrid, exx_is_active, start_exx, stop_exx
  USE input_parameters,      ONLY : exxdiv_treatment
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  ! the k-point index in hte orignal mesh of k-points
  !
  INTEGER :: ibnd, ir, k, dim_ham, nspin_aux, n_orb, lrwfc
  INTEGER i_start, i_end, ik_eff
  !
  !
  REAL(DP) :: n_r(dfftp%nnr), num1, num2, sh, aux_r(dfftp%nnr) 
  ! ... orbital density in rela space
  !
  COMPLEX(DP) n_g(ngm), n_g_aux(ngm,nspin), aux_g(ngm)
  ! ... orbital density in G space
  !
  REAL(DP) :: ehart, v(dfftp%nnr,nspin), vxc_minus1(dfftp%nnr,2), vxc(dfftp%nnr,nspin), charge, w1
  !  
  TYPE (scf_type) :: rho_minus1
  !
  REAL(DP) ::  vtxc, etxc, vtxc_minus1, etxc_minus1, etmp, delta_eig(nbnd)
  COMPLEX(DP) :: etmp2
  REAL(DP) :: etmp1
  !
  COMPLEX(DP) , ALLOCATABLE :: psic_1(:) , eigvc_ki(:,:)
  COMPLEX(DP) , ALLOCATABLE :: ham (:,:), ham_up(:,:), ham_dw(:,:), vpsi(:), vpsi_r(:), ham_aux(:,:), v_ki(:,:)
  REAL(DP), ALLOCATABLE :: eigvl_ki(:), et_aux(:,:)
  !
  LOGICAL :: off_diag = .TRUE.
  !
  !
  ALLOCATE (psic_1( dfftp%nnr), vpsi_r(dffts%nnr), vpsi(npwx), v_ki(npwx,nbnd))
  ALLOCATE (et_aux(nbnd,nks))
  !
  !
  WRITE(stdout,'(/,2x,"#########",3x, "STARTING KI corrections calcualtion", 3x, "##########",/ )')
  !
  ! ... Loop over k_point: actually it's a loop over the spin ik=1 ==> spin_up; ik=2 ==> spin_dw ...
  !
  w1 = 1.D0 / omega
  !
  !alpha_final_full(:,:)=1.0  ! Just for debug
  !
  nspin_aux=nspin
  nspin=2
  CALL create_scf_type (rho_minus1)
  nspin=nspin_aux
  !
  lrwfc = num_wann*npwx
  CALL get_buffer ( evc0, lrwfc, iuwfc_wann, ik )
  ! Retrive the ks function at k 
  IF (kc_iverbosity .gt. 0 ) WRITE(stdout,'(8X, "INFO: u_k(g) RETRIEVED"/)')
  !
  CALL compute_map_ikq_single (ik)
  ! find tha map k+q --> k'+G and store the res 
  !
  dim_ham = num_wann
  ALLOCATE ( ham (dim_ham,dim_ham) )
  ALLOCATE ( ham_up (dim_ham,dim_ham) )
  ALLOCATE ( ham_dw (dim_ham,dim_ham) )
  ham (:,:) = (0.D0,0.D0)
  ham_up (:,:) = (0.D0,0.D0)
  ham_dw (:,:) = (0.D0,0.D0)
  !
  ! ... KS Hamiltonian ....
  ik_eff = ik + (spin_component -1)*nkstot/nspin
  CALL ks_hamiltonian (evc0, ik_eff, dim_ham)
  !
  v_ki(:,:) = (0.D0,0.D0)
  !GOTO 101
  !
  n_orb = num_wann
  orb_loop: DO ibnd = 1, n_orb
     !
     IF (kc_at_ks .AND. homo_only .AND. (ibnd .ne. num_wann) ) CYCLE orb_loop ! only homo orbilal (for fast debug) 
     ! 
     IF ( lsda ) current_spin = spin_component
     ! 
     ! ... Compute the orbital density ...
     !

     npw = ngk(ik)
     psic(:) = ( 0.D0, 0.D0 )
     psic(dffts%nl(igk_k(1:npw,ik))) = evc0(1:npw,ibnd)
     CALL invfft ('Wave', psic, dffts)
     !
     ! Store the result needed below
     psic_1(:) = psic(:) 
     !
     ! ... orbital density in real space ...
     num1 = 0.D0
     num2 = 0.D0
     n_r(:) = 0.0
     DO ir = 1, dffts%nnr
       n_r(ir) = n_r(ir) + ( DBLE( psic(ir) )**2 + AIMAG( psic(ir) )**2 )*w1
#ifdef DEBUG
       num1 = num1 + DBLE( DBLE( psic(ir) )**2 + AIMAG( psic(ir) )**2  )*w1
       num2 = num2 + rho%of_r(ir,1)
#endif
     ENDDO 
     !
#ifdef DEBUG
     CALL mp_sum (num1, intra_bgrp_comm) 
     CALL mp_sum (num2, intra_bgrp_comm) 
     WRITE(stdout,'(2x, "orbital charge", 2F18.12)') num1/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )*omega
     WRITE(stdout,'(2x, "spin-up charge", 2F18.12)') num2/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )*omega
#endif
     !
     ! ... orbital density in reciprocal space ...
     psic(:) = (0.D0, 0.D0)
     psic(:) =  CMPLX(n_r(:),0.D0,kind=dp)
     CALL fwfft ('Rho', psic, dfftp)
     n_g(:) = psic(dfftp%nl(:))
     !
     ! ... Compute Int[e^2/|r-r'|*n_i(r')] ...
     !
     sh = 0.D0
     v(:,:)=0.D0
     n_g_aux(:,:) = (0.D0, 0.D0)
     n_g_aux(:,1) = n_g(:)
     CALL v_h( n_g_aux, ehart, charge, v )
     sh = ehart
     !
!     WRITE(stdout,'("v_hatree", 2i5, 3F15.8)') ibnd, current_spin ,REAL(v(1:3,1))
     WRITE(stdout,'("self_hatree, EXX", 2i5, 1F15.8)') ibnd, current_spin, -sh
     !
#ifdef DEBUG
     WRITE(stdout,'(2x, "orbital=", i3, 2x, "Self-Hartree", F15.10, 3x, "Ry",/)') ibnd, sh
     WRITE(stdout,'(2x,"orbital charge from v_h",F15.12,/)') charge
#endif
     !
     ! .. Add the xc contribution ...
     !
     etxc = 0.D0; vtxc = 0.D0; vxc(:,:) = 0.D0
     CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, vxc )
     !
     etmp=0.D0
     etmp = sum ( vxc(1:dfftp%nnr,current_spin) * n_r(1:dfftp%nnr) )
     etmp = etmp/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )*omega
     CALL mp_sum (etmp, intra_bgrp_comm) 
     !
     IF (nspin == 1 ) THEN
        rho_minus1%of_r(:,1) = rho%of_r(:,1)
        rho_minus1%of_r(:,2) = 0.D0
        rho_minus1%of_g(:,1) = rho%of_g(:,1)
        rho_minus1%of_g(:,2) = 0.D0
     ELSE
        rho_minus1%of_r(:,:) = rho%of_r(:,:)
        rho_minus1%of_g(:,:) = rho%of_g(:,:)
     ENDIF
     !
     etmp1 = 0.D0
     delta_eig(ibnd) = 0.D0
     !
     IF ( ibnd .LE. num_wann_occ ) THEN
        !
        rho_minus1%of_r(:,1) = rho_minus1%of_r(:,1) - n_r(:)  ! denisty rho-rho_i in real space
        rho_minus1%of_r(:,2) = rho_minus1%of_r(:,2) - n_r(:)  ! denisty rho-rho_i in real space
        rho_minus1%of_g(:,1) = rho_minus1%of_g(:,1) - n_g(:)  ! denisty rho-rho_i in reciprocal scape
        rho_minus1%of_g(:,2) = rho_minus1%of_g(:,2) - n_g(:)  ! denisty rho-rho_i in reciprocal scape
        ! 
        etxc_minus1 = 0.D0; vtxc_minus1 = 0.D0; vxc_minus1(:,:) = 0.D0
        nspin_aux=nspin; nspin=2
        CALL v_xc( rho_minus1, rho_core, rhog_core, etxc_minus1, vtxc_minus1, vxc_minus1 )
        nspin=nspin_aux
        ! 
        !
        delta_eig(ibnd) = (-sh+(etxc-etxc_minus1-etmp)) 
        WRITE(stdout, '("NICOLA KI corr const term, sh[n_i], Exc[n], Exc[n-n_i], int{v_xc[n] n_i} ", 4F14.8)') sh, etxc, etxc_minus1, etmp 
        IF (lrpa) delta_eig(ibnd) = (-sh)  !! hartree only for debug
        !
        vpsi_r(:) = (0.D0, 0.D0)
        DO ir = 1, dffts%nnr
           vpsi_r (ir) = CMPLX( ( delta_eig(ibnd)),0.D0) * psic_1(ir)
        ENDDO
        !
        CALL fwfft ('Wave', vpsi_r, dffts)
        v_ki(:,ibnd) = vpsi_r(dffts%nl(igk_k(:,ik))) 
        !
!        WRITE(stdout,'("evc_occ", i5, 6F15.8)') ibnd, evc0(1:3,ibnd)
!        WRITE(stdout,'("vpsi", i5, 6F15.8)') ibnd, v_ki(1:3,ibnd)
        !
     ELSE
        !
        rho_minus1%of_r(:,1) = rho_minus1%of_r(:,1) + n_r(:)  ! denisty rho+rho_i in real space
        rho_minus1%of_r(:,2) = rho_minus1%of_r(:,2) + n_r(:)  ! denisty rho+rho_i in real space
        rho_minus1%of_g(:,1) = rho_minus1%of_g(:,1) + n_g(:)  ! denisty rho+rho_i in reciprocal scape
        rho_minus1%of_g(:,2) = rho_minus1%of_g(:,2) + n_g(:)  ! denisty rho+rho_i in reciprocal scape
        ! 
        etxc_minus1 = 0.D0; vtxc_minus1 = 0.D0; vxc_minus1(:,:) = 0.D0
        nspin_aux=nspin; nspin=2
        CALL v_xc( rho_minus1, rho_core, rhog_core, etxc_minus1, vtxc_minus1, vxc_minus1 )
        nspin=nspin_aux
        ! 
        etmp1 = sum ( vxc_minus1(1:dfftp%nnr,current_spin) * n_r(1:dfftp%nnr) )
        etmp1= etmp1/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )*omega
        CALL mp_sum (etmp1, intra_bgrp_comm)
        !
        ! Scalar-term correction for Diagonal elements only
        delta_eig(ibnd) = -sh + (etxc_minus1 - etxc - etmp1)
        WRITE(stdout, '("NICOLA KI corr const term, sh[n_i], Exc[n], Exc[n+n_i], int{v_xc[n] n_i} ", 4F14.8)') sh, etxc, etxc_minus1, etmp 
        IF (lrpa) delta_eig(ibnd) = (-sh)  !! hartree only for debug
!        delta_eig(ibnd) = -sh
!        WRITE(stdout,'("const_term_empty_hartree", i5, 3F15.8)') ibnd, -sh
!        WRITE(stdout,'("const_term_empty_xc", i5, 3F15.8)') ibnd, etxc_minus1, etxc, etmp1
!        WRITE(stdout,'("Vhartree", i5, 6F15.8)') ibnd, v(1:3,current_spin) 
        ! 
        ! Non-scalar component of the KI functional for both diagonal and off-diagonal elements
        ! This is the KI gradient (the non-scalar part)
        vpsi_r(:) = (0.D0, 0.D0)
        etmp2 = (0.D0, 0.D0)
        DO ir = 1, dffts%nnr
           vpsi_r (ir) = CMPLX( ( v(ir,current_spin) + vxc_minus1(ir,current_spin) - vxc(ir,current_spin) + delta_eig(ibnd)),0.D0) * psic_1(ir)
           IF(lrpa) vpsi_r (ir) = CMPLX( ( v(ir,current_spin) + delta_eig(ibnd)),0.D0) * psic_1(ir)
           etmp2 = etmp2 + CONJG(psic_1(ir))*vpsi_r(ir)
        ENDDO 
!!        WRITE(*,'("NICOLA i, vi_r", i5, 10F16.8)') ibnd, v(1:5,spin_component)
        etmp2=etmp2/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
        CALL mp_sum (etmp2, intra_bgrp_comm)
        delta_eig(ibnd) = etmp2 
        !
        ! 1) GO to G-space and store the ki gradient 
        CALL fwfft ('Wave', vpsi_r, dffts)
        v_ki(:,ibnd) = vpsi_r(dffts%nl(igk_k(:,ik)))
        !
!        WRITE(stdout,'("evc_empty", i5, 6F15.8)') ibnd, evc0(1:3,ibnd)
!        WRITE(stdout,'("vpsi", i5, 6F15.8)') ibnd, v_ki(1:3,ibnd)
        !
     ENDIF
     !
     IF (kipz_corr) THEN !! Add the PZ part of the KIPZ potential (this is in the spirit of Perturbation theory) 
        ! 
        WRITE(stdout,'(5x, "ADDING KIPZ CORRECTION ....")') 
        !
        ! Use rho_minus as a workspace 
        rho_minus1%of_r(:,:) = 0.D0
        rho_minus1%of_g(:,:) = (0.D0, 0.D0) 
        rho_minus1%of_r(:,1) = n_r(:)  ! orbital denisty rho_i in real space
        rho_minus1%of_r(:,2) = n_r(:)  ! orbital denisty rho_i in real space
        rho_minus1%of_g(:,1) = n_g(:)  ! orbital denisty rho_i in reciprocal scape
        rho_minus1%of_g(:,2) = n_g(:)  ! orbital denisty rho_i in reciprocal scape
        !
        etxc_minus1 = 0.D0; vtxc_minus1 = 0.D0; vxc_minus1(:,:) = 0.D0
        nspin_aux=nspin; nspin=2
        aux_r =0.D0 ; aux_g = (0.D0, 0.D0)
        CALL v_xc( rho_minus1, aux_r, aux_g, etxc_minus1, vtxc_minus1, vxc_minus1 )
        nspin=nspin_aux
        !
        ! The constant term 
        etmp1 =0.D0
        etmp1 = sum ( vxc_minus1(1:dfftp%nnr,current_spin) * n_r(1:dfftp%nnr) )
        etmp1= etmp1/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )*omega
        CALL mp_sum (etmp1, intra_bgrp_comm)
        WRITE(stdout , '("NICOLA PZ corr const term, sh[n_i], Exc[n_i], int{v_xc[n_i] n_i}, int{v_xc[n_i] n_i}", 4F15.8)'), sh, etxc_minus1, etmp1, vtxc_minus1
        etmp1 = + sh - etxc_minus1 + etmp1
        !
        vpsi_r(:) = (0.D0, 0.D0)
        etmp2 = (0.D0, 0.D0)
        DO ir = 1, dffts%nnr
           vpsi_r (ir) = CMPLX( ( etmp1 - v(ir,current_spin) - vxc_minus1(ir,current_spin) ),0.D0) * psic_1(ir)
           etmp2 = etmp2 + CONJG(psic_1(ir))*vpsi_r(ir)
!           vpsi_r (ir) = CMPLX( ( v(ir,current_spin) + delta_eig(ibnd)),0.D0) * psic_1(ir)
        ENDDO
        ! The diagonal term. It's the shift of the KS eigenvalue in the canonical representation
        etmp2=etmp2/( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
        CALL mp_sum (etmp2, intra_bgrp_comm)
        delta_eig(ibnd) = delta_eig(ibnd) + REAL(etmp2)
        WRITE(stdout,*) "NICOLA Delta KIPZ",  REAL(etmp2)
        !
        ! 1) GO to G-space and store the kipz gradient 
        CALL fwfft ('Wave', vpsi_r, dffts)
        v_ki(:,ibnd) = v_ki(:,ibnd) + vpsi_r(dffts%nl(igk_k(:,ik)))
        !
        WRITE(stdout,'(5x, ".... DONE ")') 
        !
     ENDIF
     !
     IF (alpha_final_full(ibnd) .gt. 1.02 ) THEN 
        WRITE(stdout,'("WARNING: alpha for orbital", i5, i3, "  bigger than 1.02.", F15.8, "Set it to 1.00",/)') ibnd, ik, alpha_final_full(ibnd)
         alpha_final_full(ibnd) = 1.D0
     ENDIF 
     !
     IF (alpha_final_full(ibnd) .lt. 0.00 ) THEN 
        WRITE(stdout,'("WARNING: alpha for orbital", i5, i3, "  smaller than 0.00.", F15.8, "Set it to 1.00",/)') ibnd, ik, alpha_final_full(ibnd)
         alpha_final_full(ibnd) = 1.D0
     ENDIF 
     !
     v_ki(:,ibnd) = v_ki(:,ibnd) * alpha_final_full(ibnd)
     !
     WRITE(stdout,'(3x, "orbital", i3, 3x, "spin", i3, 5x, "uKI_diag", F15.8 ," Ry", 3x, "rKI_diag", F15.8, " Ry", 3x, "alpha=", F15.8, 3x )') &
                           ibnd, current_spin, delta_eig(ibnd), delta_eig(ibnd)*alpha_final_full(ibnd), alpha_final_full(ibnd)
     !
  ENDDO orb_loop
  !
101  CONTINUE
  !
  ! ##### Build up the KI Hamiltonian 
  !
  ham_up(:,:)= (0.D0,0D0)
  ham_dw(:,:)= (0.D0,0D0)
  DO ibnd = 1, n_orb
     ! 
     DO k = ibnd, n_orb
        !
        etmp2 = 0.D0
        DO ir = 1, npw
!           etmp2 = etmp2 + CONJG(evc0(ir,ibnd)) * (v_ki(ir,k)+hpsi(ir,k))
           etmp2 = etmp2 + CONJG(evc0(ir,ibnd)) * (v_ki(ir,k))
        ENDDO
        CALL mp_sum (etmp2, intra_bgrp_comm)
!        WRITE(*,'("Nic: KI    matrix ele", 2i5, 3x, 2F20.15)') ibnd, k , etmp2
!        WRITE(*,'("Nic: KS    matrix ele", 2i5, 3x, 2F20.15)') ibnd, k , Hamlt(ik,ibnd,k)
!        WRITE(*,'("Nic: KS+KI matrix ele", 2i5, 3x, 2F20.15)') ibnd, k , Hamlt(ik,ibnd,k)+etmp2
        !
        ham_up(ibnd,k) = etmp2
        ! kill occ-empty matrix elements
        IF ( (ibnd .le. num_wann_occ)  .AND. (k .gt. num_wann_occ) ) ham_up (ibnd, k) = (0.D0,0D0)  
        ham_up(k,ibnd) = CONJG(ham_up(ibnd,k))
        !
     ENDDO
     !
     IF (qp_symm) THEN 
        !
        ! ... The lower part of the Hamiltonian ...
        ! ... In apost processing fashion the KI hamiltonian is not hermitian ...
        ! ... Here we build an Hermitian operator in the spirit of scQP GW ...
        DO k = 1, ibnd
           !
        etmp2 = 0.D0
           DO ir = 1, npw
!              etmp2 = etmp2 + CONJG(evc0(ir,ibnd)) * (v_ki(ir,k)+hpsi(ir,k))
              etmp2 = etmp2 + CONJG(evc0(ir,ibnd)) * (v_ki(ir,k))
           ENDDO
           CALL mp_sum (etmp2, intra_bgrp_comm)
           !
           ham_dw(ibnd,k) = etmp2
           IF (ibnd .gt. num_wann_occ .AND. k .le. num_wann_occ ) ham_dw (ibnd, k) = (0.D0,0D0)  ! kill occ-empty matrix elements
           ham_dw(k,ibnd) = CONJG(ham_dw(ibnd,k))
           !
        ENDDO
        !
     ENDIF
     !
     ! Shift of the KS eigenvalues (only diagonal KI)
!     delta_eig(ibnd) = DBLE( ham_up(ibnd,ibnd) )-et(ibnd,ik)
     delta_eig(ibnd) = DBLE( ham_up(ibnd,ibnd) )
!     WRITE(stdout,'("NICOLA", 3F15.8,/)')  DBLE( ham_up(ibnd,ibnd) ), DBLE(ham_dw(ibnd,ibnd))
     !
  ENDDO
  !
!!!!!! DEBUG
  IF ( .NOT. off_diag )THEN 
    ham_up = (0.D0,0.D0)
    ham_up = (0.D0,0.D0)
    DO ibnd = 1, n_orb; ham_up(ibnd,ibnd)= delta_eig(ibnd); enddo
  ENDIF
!!!!!!!!!!
  !
  WRITE(stdout,'(/,3x, "###  KI HAMILTONIAN EMPTY ###")')
  DO k = num_wann_occ+1, num_wann_occ+4
     WRITE(stdout,'(16F20.15)') ham_up(num_wann_occ+1:num_wann_occ+4,k)
  ENDDO
  WRITE(stdout,*) 
  !
  OPEN(987,FILE='hamiltonian_emp.dat',FORM='formatted',status='UNKNOWN')
  DO k = num_wann_occ+1, num_wann
     WRITE(987,'(2E18.10)') ham_up(num_wann_occ+1:num_wann, k)
  ENDDO
  CLOSE(987)
  !
  ham=(0.D0,0.D0)
  !

  !
  ! The KS hamiltonian in the Wannier Gauge (just to check)
  ham(:,:) = Hamlt(ik,:,:)+ham_up(:,:)
  IF (qp_symm)  ham(:,:) = Hamlt(ik,:,:)+0.5D0*(ham_up(:,:)+ham_dw(:,:))
!  WRITE(stdout,'("NICOLA", 2F15.8)')  ham(2,3), ham(3,2) 
  !
  ! Store the res in the global variable
  Hamlt(ik,:,:) = ham(:,:)
  !
  WRITE(stdout,'(/,/,2x,"#########",3x, "END KI corrections calcualtion", 3x, "##########",/ )')
  ! 
#ifdef DEBUG
  WRITE(stdout,'(/,3x, "###  KI HAMILTONIAN ###")')
  i_end = MIN(8,n_orb)
  DO k = 1, i_end
     WRITE(stdout,'(16F12.7)') ham(1:i_end,k)
  ENDDO
  WRITE(stdout,*) 
  !
  WRITE(stdout,'(/,3x, "###  KI HAMILTONIAN EMPTY ###")')
  DO k = num_wann_occ+1, num_wann_occ+4
     WRITE(stdout,'(16F12.7)') ham(num_wann_occ+1:num_wann_occ+4,k)
  ENDDO
  WRITE(stdout,*) 
#endif
  !
  ! Perturbative approach: ONLY diagonal correction
  et_aux(:,ik) = 0.D0
  DO ibnd = 1, n_orb
     et_aux(ibnd,ik) = DBLE(ham(ibnd,ibnd))
  ENDDO
  !
  ! Here diagonalize ham for empty states with 
  ! different dimension of the hilbert space
  !
  WRITE(stdout,'(/,3x, "###  Empty states spectrum ###")')
!  DO k = 1, dim_ham
  DO k = 1, dim_ham-num_wann_occ
     !
     i_start = num_wann_occ+1
     i_end = num_wann_occ+k
     !
     ALLOCATE (ham_aux(k,k), eigvl_ki(k), eigvc_ki(k,k))
     !ham_aux(1:k,1:k) = ham(1:k,1:k)
     ham_aux(1:k,1:k) = ham(i_start:i_end,i_start:i_end)
     !
     CALL cdiagh( k, ham_aux, k, eigvl_ki, eigvc_ki )
     !
     IF (k.le.10) THEN 
        WRITE(stdout,'(I3, 10F10.4)') k, eigvl_ki(1:k)*rytoev  ! First 10 eigenvalues
     ELSE
        WRITE(stdout,'(I3, 10F10.4)') k, eigvl_ki(1:10)*rytoev  ! First 10 eigenvalues
     ENDIF
     !
     DEALLOCATE (ham_aux)
     DEALLOCATE (eigvl_ki, eigvc_ki)
     !
  ENDDO
  ! 
  ALLOCATE (ham_aux(n_orb,n_orb), eigvl_ki(n_orb), eigvc_ki(n_orb,n_orb))
  ham_aux(1:n_orb,1:n_orb) = ham(1:n_orb,1:n_orb)
  CALL cdiagh( n_orb, ham_aux, n_orb, eigvl_ki, eigvc_ki )
  DO ibnd = 1, n_orb
     et(ibnd,ik) = eigvl_ki(ibnd)
  ENDDO
  !
  DEALLOCATE (ham) 
  DEALLOCATE (ham_up) 
  DEALLOCATE (ham_dw) 
  IF (ALLOCATED(eigvl_ki)) DEALLOCATE (eigvl_ki)
  IF (ALLOCATED(eigvc_ki)) DEALLOCATE (eigvc_ki)
  IF (ALLOCATED(ham_aux)) DEALLOCATE (ham_aux)
  !
  !
  write(*,*) n_orb
  WRITE( stdout, '("KI Full ",8F9.4)' ) (et(ibnd,ik)*rytoev, ibnd=1,n_orb) 
  IF (kc_at_ks ) WRITE( stdout, '("KI Pert ",8F9.4)' ) (et_aux(ibnd,ik)*rytoev, ibnd=1,n_orb)
  !
  DEALLOCATE (psic_1, vpsi_r, vpsi, v_ki) 
  DEALLOCATE (et_aux)
  !
END subroutine full_ham
