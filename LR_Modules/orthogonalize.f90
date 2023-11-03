!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE orthogonalize(dvpsi, evq, ikk, ikq, dpsi, npwq, dpsi_computed)
!------------------------------------------------------------------------
  !
  ! This routine orthogonalizes dvpsi to the valence states: ps = <evq|dvpsi>
  ! It should be quite general. It works for metals and insulators, with
  ! NC as well as with US PP, both SR or FR.
  ! Note that on output it changes sign. So it applies -P_c^+.
  ! P_c^+ = 1 - |dpsi><evq| = 1 - S|evq><evq|
  !
  ! NB: IN/OUT is dvpsi ; 
  ! dpsi is used as work space: dpsi = S*evq
  !
  ! If dpsi_computed=.true.  then dpsi=S*evq has been already computed;
  ! If dpsi_computed=.false. then dpsi=S*evq must be computed here.
  !
  ! The variable dpsi_computed has been introduced in order to make 
  ! this subroutine more general, because evq are the ground-state wfcts,
  ! which do not change during the linear-response calculation and 
  ! hence it can be useful to compute S*evq just once and use it 
  ! throughout the whole calculation (e.g., as in TDDFPT for k=0). 
  !
  USE kinds,            ONLY : DP
  USE klist,            ONLY : lgauss, degauss, ngauss, ltetra, wk
  USE noncollin_module, ONLY : noncolin, npol
  USE wvfct,            ONLY : npwx, nbnd, wg, et
  USE ener,             ONLY : ef
  USE becmod,           ONLY : bec_type, becp, calbec
  USE becmod_subs_gpum, ONLY : using_becp_auto
  USE uspp,             ONLY : vkb, okvan
  USE mp_bands,         ONLY : use_bgrp_in_hpsi, inter_bgrp_comm, intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  USE xc_lib,           ONLY : exx_is_active
  USE control_flags,    ONLY : gamma_only
  USE gvect,            ONLY : gstart
  USE control_lr,       ONLY : alpha_pv, nbnd_occ
  USE dfpt_tetra_mod,   ONLY : dfpt_tetra_beta
#if defined(__CUDA)
  USE becmod_subs_gpum, ONLY : calbec_gpu, using_becp_d_auto
  USE becmod_gpum,      ONLY : becp_d
  USE cublas
#endif
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ikk, ikq   ! the index of the k and k+q points
  INTEGER, INTENT(IN) :: npwq       ! the number of plane waves for q
  COMPLEX(DP), INTENT(IN)    :: evq(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: dvpsi(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: dpsi(npwx*npol,nbnd)
  LOGICAL, INTENT(IN) :: dpsi_computed
  !
  COMPLEX(DP), ALLOCATABLE :: ps(:,:)
  REAL(DP), ALLOCATABLE    :: ps_r(:,:)
  INTEGER :: ibnd, jbnd, nbnd_eff, n_start, n_end
  REAL(DP) :: wg1, w0g, wgp, wwg(nbnd), deltae, theta
  REAL(DP), EXTERNAL :: w0gauss, wgauss
  ! functions computing the delta and theta function
  !
  CALL start_clock ('ortho')
  !
  IF (gamma_only) THEN
     ALLOCATE(ps_r(nbnd,nbnd))
  ENDIF
  !
  ALLOCATE(ps(nbnd,nbnd))
  !
  !$acc data copyin(evq) copy(dvpsi,dpsi) create(ps(1:nbnd, 1:nbnd), ps_r(1:nbnd, 1:nbnd))
  IF (gamma_only) THEN
     !$acc kernels
     ps_r(:,:) = 0.0d0
     !$acc end kernels
  ENDIF        
  !$acc kernels
  ps(:,:) = (0.0d0, 0.0d0)
  !$acc end kernels
  !
  IF (ltetra .OR. lgauss) THEN
     !
     !  metallic case
     !
     IF (gamma_only) CALL errore ('orthogonalize', 'smearing or tetrahedra &
         & with gamma-point algorithm?',1)
     !
     !$acc host_data use_device(evq, dvpsi, ps)
     IF (noncolin) THEN
        CALL zgemm( 'C', 'N', nbnd, nbnd_occ (ikk), npwx*npol, (1.d0,0.d0), &
             evq, npwx*npol, dvpsi, npwx*npol, (0.d0,0.d0), ps, nbnd )
     ELSE
        CALL zgemm( 'C', 'N', nbnd, nbnd_occ (ikk), npwq, (1.d0,0.d0), &
             evq, npwx, dvpsi, npwx, (0.d0,0.d0), ps, nbnd )
     END IF
     !$acc end host_data
     !
     DO ibnd = 1, nbnd_occ (ikk)
        !
        IF ( lgauss ) THEN
           !
           wg1 = wgauss ((ef-et(ibnd,ikk)) / degauss, ngauss)
           w0g = w0gauss((ef-et(ibnd,ikk)) / degauss, ngauss) / degauss
           DO jbnd = 1, nbnd
              wgp = wgauss ( (ef - et (jbnd, ikq) ) / degauss, ngauss)
              deltae = et (jbnd, ikq) - et (ibnd, ikk)
              theta = wgauss (deltae / degauss, 0)
              wwg(jbnd) = wg1 * (1.d0 - theta) + wgp * theta
              IF (jbnd <= nbnd_occ (ikq) ) THEN
                 IF (abs (deltae) > 1.0d-5) THEN
                    wwg(jbnd) = wwg(jbnd) + alpha_pv * theta * (wgp - wg1) / deltae
                 ELSE
                    !
                    !  if the two energies are too close takes the limit
                    !  of the 0/0 ratio
                    !
                    wwg(jbnd) = wwg(jbnd) - alpha_pv * theta * w0g
                 ENDIF
              ENDIF
              !
           ENDDO
           !
           !$acc parallel loop copyin(wwg)
           DO jbnd = 1, nbnd
              ps(jbnd,ibnd) = wwg(jbnd) * ps(jbnd,ibnd)
           END DO
           !$acc end parallel loop
           !
        ELSE
           !
           wg1 = wg(ibnd,ikk) / wk(ikk)
           !
           DO jbnd = 1, nbnd
              !
              wwg(jbnd) = dfpt_tetra_beta(jbnd,ibnd,ikk)
              !
              !
           ENDDO
           !
           !$acc kernels copyin(wwg)
           DO jbnd = 1, nbnd
              ps(jbnd,ibnd) = wwg(jbnd) * ps(jbnd,ibnd)
           END DO
           !$acc end kernels
           !
        ENDIF
        !
        IF (noncolin) THEN
           !$acc kernels 
           dvpsi(1:npwx*npol, ibnd) = wg1 * dvpsi(1:npwx*npol, ibnd)
           !$acc end kernels
        ELSE
           !$acc kernels
           dvpsi(1:npwq, ibnd) = wg1 * dvpsi(1:npwq, ibnd)
           !$acc end kernels
        END IF
        !
     END DO
     !
     nbnd_eff=nbnd
     !   
  ELSE
     !
     !  insulators
     !
     !$acc host_data use_device(evq, dvpsi, ps)
     IF (noncolin) THEN
        CALL zgemm( 'C', 'N',nbnd_occ(ikq), nbnd_occ(ikk), npwx*npol, &
             (1.d0,0.d0), evq, npwx*npol, dvpsi, npwx*npol, &
             (0.d0,0.d0), ps, nbnd )
     ELSEIF (gamma_only) THEN
        !$acc host_data use_device(ps_r)
        CALL dgemm( 'C', 'N', nbnd_occ(ikq), nbnd_occ (ikk), 2*npwq, &
             2.0_DP, evq, 2*npwx, dvpsi, 2*npwx, &
             0.0_DP, ps_r, nbnd )
        IF (gstart == 2 ) THEN
           CALL mydger( nbnd_occ(ikq), nbnd_occ (ikk), -1.0_DP, evq, &
                & 2*npwq, dvpsi, 2*npwx, ps_r, nbnd )
        ENDIF
        !$acc end host_data
     ELSE
        CALL zgemm( 'C', 'N', nbnd_occ(ikq), nbnd_occ (ikk), npwq, &
             (1.d0,0.d0), evq, npwx, dvpsi, npwx, &
             (0.d0,0.d0), ps, nbnd )
     END IF
     !$acc end host_data
     !
     nbnd_eff=nbnd_occ(ikk)
     !
  END IF
  !
  IF (gamma_only) THEN
     !$acc host_data use_device(ps_r)
     CALL mp_sum(ps_r(:,:),intra_bgrp_comm)
     !$acc end host_data
  ELSE
     !$acc host_data use_device(ps)
     CALL mp_sum(ps(:,1:nbnd_eff),intra_bgrp_comm)
     !$acc end host_data
  ENDIF
  !
  ! dpsi is used as work space to store S*evq
  !
  IF (.NOT.dpsi_computed) THEN
     !
     IF (okvan) then
        if (use_bgrp_in_hpsi .AND. .NOT. exx_is_active() .AND. nbnd_eff > 1) then
           call divide(inter_bgrp_comm,nbnd_eff, n_start, n_end)
#if defined(__CUDA)
           if (n_end >= n_start) then
              CALL using_becp_d_auto(2)
              !$acc host_data use_device(vkb, evq)
              CALL calbec_gpu (npwq, vkb, evq(:,n_start:n_end), becp_d, n_end- n_start + 1)    !
              !$acc end host_data
           endif
#else
           if ( n_end >= n_start) then
              CALL using_becp_auto(2)
              CALL calbec ( npwq, vkb, evq(:,n_start:n_end), becp, n_end - n_start + 1 )
           endif
#endif
        else
#if defined(__CUDA)
           CALL using_becp_d_auto(2)
           !$acc host_data use_device(vkb, evq)
           CALL calbec_gpu(npwq, vkb, evq, becp_d, nbnd_eff)
           !$acc end host_data
#else
           CALL using_becp_auto(2)
           CALL calbec ( npwq, vkb, evq, becp, nbnd_eff )
#endif
        end if
     end if
     !
#if defined(__CUDA)
     !$acc host_data use_device(evq, dpsi)
     CALL s_psi_gpu (npwx, npwq, nbnd_eff, evq, dpsi)
     !$acc end host_data
#else
     CALL s_psi (npwx, npwq, nbnd_eff, evq, dpsi)
#endif
     !
  ENDIF
  !
  ! |dvspi> =  -(|dvpsi> - S|evq><evq|dvpsi>)
  !
  IF (gamma_only) THEN
     !$acc kernels 
     ps = CMPLX (ps_r,0.0_DP, KIND=DP)
     !$acc end kernels
  ENDIF        
  !
  !$acc host_data use_device(dpsi, ps, dvpsi)
  IF (lgauss .OR. ltetra ) THEN
     !
     !  metallic case
     !
     IF (noncolin) THEN
        CALL zgemm( 'N', 'N', npwx*npol, nbnd_occ(ikk), nbnd, &
             (1.d0,0.d0), dpsi, npwx*npol, ps, nbnd, (-1.0d0,0.d0), &
             dvpsi, npwx*npol )
     ELSE
        CALL zgemm( 'N', 'N', npwq, nbnd_occ(ikk), nbnd, &
             (1.d0,0.d0), dpsi, npwx, ps, nbnd, (-1.0d0,0.d0), &
             dvpsi, npwx )
     END IF
     !
  ELSE
     !
     !  Insulators: note that nbnd_occ(ikk)=nbnd_occ(ikq) in an insulator
     !
     IF (noncolin) THEN
        CALL zgemm( 'N', 'N', npwx*npol, nbnd_occ(ikk), nbnd_occ(ikk), &
             (1.d0,0.d0),dpsi,npwx*npol,ps,nbnd,(-1.0d0,0.d0), &
             dvpsi, npwx*npol )
     ELSEIF (gamma_only) THEN             
        CALL ZGEMM( 'N', 'N', npwq, nbnd_occ(ikk), nbnd_occ(ikk), &
             (1.d0,0.d0), dpsi, npwx, ps, nbnd, (-1.0d0,0.d0), &
             dvpsi, npwx )
     ELSE
        CALL zgemm( 'N', 'N', npwq, nbnd_occ(ikk), nbnd_occ(ikk), &
             (1.d0,0.d0), dpsi, npwx, ps, nbnd, (-1.0d0,0.d0), &
             dvpsi, npwx )
     END IF
     !
  ENDIF
  !$acc end host_data
  !
  !$acc end data
  !
  DEALLOCATE(ps)
  !
  IF (gamma_only) THEN
     DEALLOCATE(ps_r)
  ENDIF
  !
  CALL stop_clock ('ortho')
  !
  RETURN
  !
END SUBROUTINE orthogonalize
