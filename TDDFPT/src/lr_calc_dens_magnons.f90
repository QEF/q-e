!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE lr_calc_dens_magnons (drhoscf, dpsi, L_dag)
  !-------------------------------------------------------------------------
  !
  ! Calculates response charge-density from linear-response
  ! orbitals and ground-state orbitals (noncollinear case).
  ! See Eq.(36) in B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007).
  !
  ! It does:
  ! 1) Calculates the response charge-density
  ! 2) Calculates the contribution to the charge-density due to US PP's
  ! 3) Sums up the normal and ultrasoft terms
  ! 4) Symmetrizes the total response charge-density
  !
  ! Written by Iurii Timrov (2013)
  !
  USE kinds,                 ONLY : DP
  USE ions_base,             ONLY : nat
  USE gvect,                 ONLY : ngm, g
  USE fft_base,              ONLY : dffts, dfftp
  USE klist,                 ONLY : xk, wk, igk_k, ngk
  USE wvfct,                 ONLY : nbnd, npwx
  USE gvecw,                 ONLY : gcutw
  USE qpoint,                ONLY : nksq, ikks, ikqs
  USE wavefunctions,         ONLY : evc
  USE noncollin_module,      ONLY : npol, nspin_mag
  USE uspp_param,            ONLY : nhm
  USE uspp,                  ONLY : okvan, vkb
  USE lsda_mod,              ONLY : nspin
  USE mp_global,             ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp,                    ONLY : mp_sum
  USE io_files,              ONLY : iunwfc, nwordwfc
  USE buffers,               ONLY : get_buffer
  USE fft_interfaces,        ONLY : fft_interpolate
  USE lr_variables,          ONLY : iunTwfc
  USE cell_base,             ONLY : omega
  USE control_lr,            ONLY : nbnd_occ, nbnd_occx
  USE fft_interfaces,        ONLY : invfft
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(in) :: dpsi(npwx*npol,nbnd_occx,nksq,2)
  ! input: the perturbed wavefunctions
  LOGICAL, INTENT(in)     :: L_dag
  ! input: in the case of L^+ computes rho(X,-Y)
  COMPLEX(DP), INTENT(out) :: drhoscf(dfftp%nnr,nspin_mag)
  ! input/output: the accumulated change of the charge density
  !
  COMPLEX(DP), ALLOCATABLE :: drhoscfh(:,:), dbecsum(:,:,:), dbecsum_nc(:,:,:,:)
  ! the accumulated change to the charge density on the smooth mesh
  ! the ultrasoft term
  INTEGER ::  is,  & 
              ik,  &
              ikk, & ! index of the point k
              ikq, & ! index of the point k+q
              npwq   ! number of the plane-waves at point k+q
  REAL(DP) :: weight ! weight of the k point
  !
  ! -k part
  COMPLEX(DP)  :: Tevc(npol*npwx,nbnd)      ! T-rev op. applied to u_{-k}
  !
  INTEGER :: imk
  !
  INTEGER :: ibnd, ig, npw, ir, itmp, v_siz
  REAL(DP) :: wgt
  COMPLEX(DP), ALLOCATABLE :: psi (:,:)
  COMPLEX(DP), ALLOCATABLE :: dpsic (:,:)
  !
  ! For device buffer
#if defined(__CUDA)
  INTEGER, POINTER, DEVICE :: nl_d(:)
  !
  nl_d  => dffts%nl_d
#else
  INTEGER, ALLOCATABLE :: nl_d(:)
  !
  ALLOCATE( nl_d(dffts%ngm) )
  nl_d  = dffts%nl
#endif
  !
  ALLOCATE (psi  (dffts%nnr, npol))
  ALLOCATE (dpsic(dffts%nnr, npol))
  !
  CALL start_clock_gpu ('lr_calc_dens')
  !
  ALLOCATE (drhoscfh(dffts%nnr,nspin_mag))
  !
  v_siz = dffts%nnr
  !
  !$acc data copyin(dpsi(1:npwx*npol,1:nbnd_occx,1:nksq,1:2), Tevc(1:npol*npwx,1:nbnd)) copyout(drhoscfh(1:v_siz,1:nspin_mag)) create(dpsic(1:v_siz,1:npol), psi(1:v_siz,1:npol)) present(igk_k) deviceptr(nl_d)
  !
  !$acc kernels
  drhoscfh(:,:) = (0.0d0, 0.0d0)
  !$acc end kernels
  !
  ! evc  contains the unperturbed wavefunctions of this k point
  !
  !-----------------------------
  ! RESONANT PART
  !-----------------------------
  !
  ! 
  DO ik = 1, nksq
     !
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npw  = ngk(ikk)
     npwq = ngk(ikq)
     !
     ! Read the unperturbed wavefuctions evc(k)
     !
     !IF (nksq > 1) 
     CALL get_buffer (evc, nwordwfc, iunwfc, ikk)
     !
     !$acc update device(evc)
     ! Calculation of the response charge density
     !
     wgt = wk(ikk)/omega
     !
     DO ibnd = 1, nbnd_occ(ikk)
        !
        ! FFT to R-space of the unperturbed wfct's evc
        !
        !$acc kernels
        psi(:,:) = (0.d0, 0.d0)
        !$acc end kernels
        !
        !$acc parallel loop 
        DO ig = 1, npw   ! should have +k behaviour. Check
           itmp = nl_d ( igk_k(ig,ikk) )
           psi(itmp,1) = evc(ig,ibnd)
           psi(itmp,2) = evc(ig+npwx,ibnd)
        ENDDO
        !
        !$acc host_data use_device(psi, psi)
        CALL invfft ('Wave', psi(:,1), dffts)
        CALL invfft ('Wave', psi(:,2), dffts)
        !$acc end host_data
        !
        ! FFT to R-space of the perturbed wfct's dpsi
        !
        !$acc kernels
        dpsic(:,:) = (0.d0, 0.d0)
        !$acc end kernels
        !$acc parallel loop
        DO ig = 1, npwq
           itmp = nl_d (igk_k(ig,ikq))
           dpsic(itmp,1) = dpsi(ig,ibnd,ik,1)
           dpsic(itmp,2) = dpsi(ig+npwx,ibnd,ik,1)
        ENDDO
        !
        !$acc host_data use_device(dpsic)
        CALL invfft ('Wave', dpsic(:,1), dffts)
        CALL invfft ('Wave', dpsic(:,2), dffts)
        !$acc end host_data
        !
        ! Calculation of the response charge density
        !
        !$acc parallel loop
        DO ir = 1, v_siz !dffts%nnr
           ! n'^{res}
           drhoscfh(ir,1) = drhoscfh(ir,1) + wgt * (CONJG(psi(ir,1))*dpsic(ir,1) + &
                                                 &  CONJG(psi(ir,2))*dpsic(ir,2) )
           ! mx'^{res}
           drhoscfh(ir,2) = drhoscfh(ir,2) + wgt * (CONJG(psi(ir,1))*dpsic(ir,2) + &
                                                &   CONJG(psi(ir,2))*dpsic(ir,1) )
           ! my'^{res}
           drhoscfh(ir,3) = drhoscfh(ir,3) + wgt * (CONJG(psi(ir,1))*dpsic(ir,2) - &
                                                &   CONJG(psi(ir,2))*dpsic(ir,1) ) * (0.d0,-1.d0)
           ! mz'^{res}
           drhoscfh(ir,4) = drhoscfh(ir,4) + wgt * (CONJG(psi(ir,1))*dpsic(ir,1)- &
                                                &   CONJG(psi(ir,2))*dpsic(ir,2) )
        ENDDO
        !
        !
     ENDDO ! loop on bands
     !
  ENDDO 
  !
  !-----------------------------
  ! ANTI-RESONANT PART
  !-----------------------------
  !
  DO ik = 1, nksq
     !
     ikk = ikks(ik)
     ikq = ikqs(ik)
     npwq = ngk(ikq)
     npw  = ngk(ikk)
     !
     if ( mod(ik,2) == 0) then
        imk  = ikk - 3      ! position of -k
     else
        imk  = ikk + 3      ! position of -k
     endif
     ! Read the unperturbed wavefuctions Tevc(-k)
     !
     Tevc(:,:) = (0.0d0, 0.0d0)
     !
     CALL get_buffer (Tevc, nwordwfc, iunTwfc, 2*ik-1)
     !
     !$acc update device(Tevc)
     !
     wgt = wk(imk)/omega
     !
     DO ibnd = 1, nbnd_occ(imk)
        !
        ! FFT to R-space of the unperturbed wfct's evc
        !
        !$acc kernels
        psi(:,:) = (0.d0, 0.d0)
        !$acc end kernels
        !$acc parallel loop
        DO ig = 1, npw   ! should have +k behaviour. Check
           itmp = nl_d ( igk_k(ig,ikk) )
           psi(itmp,1) = Tevc(ig,ibnd)
           psi(itmp,2) = Tevc(ig+npwx,ibnd)
        ENDDO
        !$acc host_data use_device(psi)
        CALL invfft ('Wave', psi(:,1), dffts)
        CALL invfft ('Wave', psi(:,2), dffts)
        !$acc end host_data
        !
        ! FFT to R-space of the perturbed wfct's dpsi
        !
        !$acc kernels
        dpsic(:,:) = (0.d0, 0.d0)
        !$acc end kernels
        !$acc parallel loop
        DO ig = 1, npwq 
           itmp = nl_d ( igk_k(ig,ikq) )
           dpsic(itmp,1) = dpsi(ig,ibnd,ik,2)
           dpsic(itmp,2) = dpsi(ig+npwx,ibnd,ik,2)
        ENDDO
        !$acc host_data use_device(dpsic)
        CALL invfft ('Wave', dpsic(:,1), dffts)
        CALL invfft ('Wave', dpsic(:,2), dffts)
        !$acc end host_data
        !
        ! In the L^+ case we compute rho(X,-Y)
        !
        !$acc kernels
        IF (L_dag) dpsic(:,:) = -dpsic(:,:)
        !$acc end kernels
        !
        ! Calculation of the response charge density
        !
        !$acc parallel loop
        DO ir = 1, v_siz ! dffts%nnr
           ! n'^{anti-res}
           drhoscfh(ir,1) = drhoscfh(ir,1) + wgt * (CONJG(psi(ir,1))*dpsic(ir,1) + &
                                                 &  CONJG(psi(ir,2))*dpsic(ir,2) )
           ! mx'^{anti-res}
           drhoscfh(ir,2) = drhoscfh(ir,2) - wgt * (CONJG(psi(ir,1))*dpsic(ir,2) + &
                                                &   CONJG(psi(ir,2))*dpsic(ir,1) )
           ! my'^{anti-res}
           drhoscfh(ir,3) = drhoscfh(ir,3) - wgt * (CONJG(psi(ir,1))*dpsic(ir,2) - &
                                                &   CONJG(psi(ir,2))*dpsic(ir,1) ) * (0.d0,-1.d0)
           ! mz'^{anti-res}
           drhoscfh(ir,4) = drhoscfh(ir,4) - wgt * (CONJG(psi(ir,1))*dpsic(ir,1)- &
                                                &   CONJG(psi(ir,2))*dpsic(ir,2) )
        ENDDO
        !
     ENDDO ! loop on bands
     !
  ENDDO ! loop on k points
  !
  ! Interpolate the charge density from the smooth mesh
  ! to a thicker mesh (if doublegrid=.true.)
  ! drhoscfh -> drhoscf
  !
  !$acc end data
  !
  DO is = 1, nspin_mag
     CALL fft_interpolate(dffts, drhoscfh(:,is), dfftp, drhoscf(:,is) )
  ENDDO
  !
  ! Reduce the charge density response across pools.
  !
#if defined(__MPI)
  CALL mp_sum (drhoscf, inter_pool_comm)
#endif
  !
  DEALLOCATE (drhoscfh)
  !
  DEALLOCATE (psi, dpsic )
  !
  CALL stop_clock_gpu ('lr_calc_dens')
  !
  RETURN
  !
END SUBROUTINE lr_calc_dens_magnons
