!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE force_us_gpu( forcenl )
  !----------------------------------------------------------------------------
  !! The nonlocal potential contribution to forces.
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE cell_base,            ONLY : at, bg, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE klist,                ONLY : nks, xk, ngk, igk_k, igk_k_d
  USE gvect_gpum,           ONLY : g_d
  USE uspp,                 ONLY : nkb, qq_at, deeq, qq_so, deeq_nc, indv_ijkb0
  USE uspp_gpum,            ONLY : vkb_d, using_vkb, using_vkb_d
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wvfct,                ONLY : nbnd, npwx, wg, et
  USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
  USE symme,                ONLY : symvector
  USE wavefunctions,        ONLY : evc
  USE wavefunctions_gpum,   ONLY : evc_d, using_evc, using_evc_d
  USE noncollin_module,     ONLY : npol, noncolin
  USE spin_orb,             ONLY : lspinorb
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE becmod,               ONLY : calbec, becp, bec_type, allocate_bec_type, &
                                   deallocate_bec_type
  USE becmod_gpum,          ONLY : becp_d, bec_type_d
  USE becmod_subs_gpum,     ONLY : using_becp_d_auto, calbec_gpu
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_get_comm_null

  USE wavefunctions_gpum,   ONLY : using_evc
  USE wvfct_gpum,           ONLY : using_et
  USE uspp_gpum,            ONLY : using_vkb, using_indv_ijkb0, using_qq_at, &
                                   using_deeq
  USE becmod_subs_gpum,     ONLY : using_becp_auto, allocate_bec_type_gpu, &
                                    synchronize_bec_type_gpu
  USE device_fbuff_m,             ONLY : dev_buf
  USE control_flags,        ONLY : use_gpu
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: forcenl(3,nat)
  !! the nonlocal contribution
  !
  ! ... local variables
  !
  COMPLEX(DP), POINTER     :: vkb1_d(:,:)   ! contains g*|beta>
#if defined(__CUDA)
  attributes(DEVICE) :: vkb1_d
#endif
  COMPLEX(DP), ALLOCATABLE :: deff_nc(:,:,:,:)
  REAL(DP), ALLOCATABLE :: deff(:,:,:)
  TYPE(bec_type)           :: dbecp                 ! contains <dbeta|psi>
  TYPE(bec_type_d), TARGET :: dbecp_d               ! contains <dbeta|psi>
  INTEGER    :: npw, ik, ipol, ig, jkb
  INTEGER    :: ierr
  !
  forcenl(:,:) = 0.D0
  !
  call start_clock_gpu('fus_allbec') 
  CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm )
  CALL using_becp_auto(2)
  CALL allocate_bec_type ( nkb, nbnd, dbecp, intra_bgrp_comm )
  CALL allocate_bec_type_gpu ( nkb, nbnd, dbecp_d, intra_bgrp_comm )
  !
  CALL dev_buf%lock_buffer( vkb1_d, (/ npwx, nkb /), ierr )
  IF (noncolin) THEN
     ALLOCATE( deff_nc(nhm,nhm,nat,nspin) )
  ELSEIF (.NOT. gamma_only ) THEN
     ALLOCATE( deff(nhm,nhm,nat) )
  ENDIF
  !
  ! ... the forces are a sum over the K points and over the bands
  !
  CALL using_evc_d(0)
  call stop_clock_gpu('fus_allbec') 
  !
  DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     npw = ngk (ik)

     IF ( nks > 1 ) THEN
        CALL get_buffer( evc, nwordwfc, iunwfc, ik )
        CALL using_evc(1)
        IF ( nkb > 0 ) CALL using_vkb_d(1)
        IF ( nkb > 0 ) &
             CALL init_us_2_gpu( npw, igk_k_d(1,ik), xk(1,ik), vkb_d )
     ENDIF
     !
     CALL using_evc_d(0); CALL using_vkb_d(0); CALL using_becp_d_auto(2)
     CALL calbec_gpu ( npw, vkb_d, evc_d, becp_d )
     !
     CALL using_evc_d(0)
     DO ipol = 1, 3
!$cuf kernel do(2) <<<*,*>>>
        DO jkb = 1, nkb
           DO ig = 1, npw
              vkb1_d(ig,jkb) = vkb_d(ig,jkb) * (0.D0,-1.D0) * g_d(ipol,igk_k_d(ig,ik))
           ENDDO
        ENDDO
        !
        CALL calbec_gpu ( npw, vkb1_d, evc_d, dbecp_d )
        CALL synchronize_bec_type_gpu(dbecp_d, dbecp, 'h')
        !
        IF ( gamma_only ) THEN
           !
           CALL force_us_gamma_gpu( forcenl )
           !
        ELSE
           !
           CALL force_us_k_gpu( forcenl )
           !
        ENDIF
     ENDDO
  ENDDO
  !
  ! ... if sums over bands are parallelized over the band group
  !
  CALL using_becp_auto(0)
  IF ( becp%comm /= mp_get_comm_null() ) CALL mp_sum( forcenl, becp%comm )
  !
  IF (noncolin) THEN
     DEALLOCATE( deff_nc )
  ELSEIF ( .NOT. GAMMA_ONLY) THEN
     DEALLOCATE( deff )
  ENDIF
  CALL dev_buf%release_buffer( vkb1_d, ierr )
  CALL deallocate_bec_type ( dbecp ) 
  CALL deallocate_bec_type ( becp )
  CALL using_becp_auto(2)
  CALL using_becp_d_auto(2)
  !
  ! ... collect contributions across pools from all k-points
  !
  CALL mp_sum( forcenl, inter_pool_comm )
  !
  ! ... The total D matrix depends on the ionic position via the
  ! ... augmentation part \int V_eff Q dr, the term deriving from the 
  ! ... derivative of Q is added in the routine addusforce
  !
  IF (use_gpu) CALL addusforce_gpu(forcenl)  
  IF (.NOT. use_gpu) CALL addusforce( forcenl )
  !
  ! ... Since our summation over k points was only on the irreducible 
  ! ... BZ we have to symmetrize the forces.
  !
  CALL symvector ( nat, forcenl )
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE force_us_gamma_gpu( forcenl )
       !-----------------------------------------------------------------------
       !! Nonlocal contributiuon. Calculation at gamma.
       !
#if defined(__CUDA)
       USE cublas
#endif
       USE uspp_gpum,            ONLY : qq_at_d, using_qq_at_d, deeq_d, using_deeq_d
       USE wvfct_gpum,           ONLY : wg_d, using_wg_d, et_d, using_et_d
       IMPLICIT NONE
       !
       REAL(DP) :: forcenl(3,nat)
       !! the nonlocal contribution
       !
       ! ... local variables
       !
       REAL(DP), POINTER :: aux_d(:,:)
#if defined(__CUDA)
       attributes(DEVICE) :: aux_d
#endif
       INTEGER ::  nt, na, ibnd, ibnd_loc, ih, jh, ijkb0 ! counters
       !
       ! CUDA Fortran workarounds
       INTEGER :: nh_nt, becp_d_ibnd_begin, becp_d_nbnd_loc
       REAL(DP), POINTER :: dbecp_d_r_d(:,:), becp_d_r_d(:,:)
       REAL(DP) :: forcenl_ipol
#if defined(__CUDA)
       attributes(DEVICE) :: dbecp_d_r_d, becp_d_r_d
#endif
       !
       ! ... Important notice about parallelization over the band group of processors:
       ! ... 1) internally, "calbec" parallelises on plane waves over the band group
       ! ... 2) the results of "calbec" are distributed across processors of the band
       ! ...    group: the band index of becp, dbecp is distributed
       ! ... 3) the band group is subsequently used to parallelize over bands
       !
       !
       CALL using_et_d(0)
       CALL using_wg_d(0)
       CALL using_indv_ijkb0(0)
       CALL using_deeq_d(0)
       CALL using_deeq(0)
       CALL using_qq_at_d(0)
       !!!!! CHECK becp (set above)
       becp_d_ibnd_begin = becp_d%ibnd_begin
       becp_d_nbnd_loc = becp_d%nbnd_loc

       dbecp_d_r_d => dbecp_d%r_d
       becp_d_r_d  => becp_d%r_d

       DO nt = 1, ntyp
          IF ( nh(nt) == 0 ) CYCLE
          CALL dev_buf%lock_buffer(aux_d, (/nh(nt),becp_d%nbnd_loc/), ierr ) !ALLOCATE ( aux(nh(nt),becp%nbnd_loc) )
          nh_nt = nh(nt)
          
          DO na = 1, nat
             IF ( ityp(na) == nt ) THEN
                ijkb0 = indv_ijkb0(na)
                ! this is \sum_j q_{ij} <beta_j|psi>
                CALL DGEMM ('N','N', nh(nt), becp_d%nbnd_loc, nh(nt), &
                     1.0_dp, qq_at_d(1,1,na), nhm, becp_d%r_d(ijkb0+1,1),&
                     nkb, 0.0_dp, aux_d, nh(nt) )
                ! multiply by -\epsilon_n
!$cuf kernel do(2)
                DO ih = 1, nh_nt
                   DO ibnd_loc = 1, becp_d_nbnd_loc
                      ibnd = ibnd_loc + becp_d_ibnd_begin - 1
                      aux_d(ih,ibnd_loc) = - et_d(ibnd,ik) * aux_d(ih,ibnd_loc)
                   END DO
                END DO

                ! add  \sum_j d_{ij} <beta_j|psi>
                CALL DGEMM ('N','N', nh(nt), becp_d%nbnd_loc, nh(nt), &
                     1.0_dp, deeq_d(1,1,na,current_spin), nhm, &
                     becp_d%r_d(ijkb0+1,1), nkb, 1.0_dp, aux_d, nh(nt) )

                ! Auxiliary variable to perform the reduction with cuf kernels
                forcenl_ipol = 0.0_dp
!$cuf kernel do(2)
                DO ih = 1, nh_nt
                   DO ibnd_loc = 1, becp_d_nbnd_loc
                      ibnd = ibnd_loc + becp_d_ibnd_begin - 1
                      forcenl_ipol = forcenl_ipol - &
                           2.0_dp * tpiba * aux_d(ih,ibnd_loc) * &
                           dbecp_d_r_d(ijkb0+ih,ibnd_loc) * wg_d(ibnd,ik)
                   ENDDO
                ENDDO
                forcenl(ipol,na) = forcenl(ipol,na) + forcenl_ipol
                !
             ENDIF
          ENDDO
          CALL dev_buf%release_buffer(aux_d, ierr)
       ENDDO
       !
     END SUBROUTINE force_us_gamma_gpu
     !     
     !-----------------------------------------------------------------------
     SUBROUTINE force_us_k_gpu( forcenl )
       !-----------------------------------------------------------------------
       !! Nonlocal contributiuon. Calculation for k-points.
       !
       IMPLICIT NONE
       !
       REAL(DP) :: forcenl(3,nat)
       !! the nonlocal contribution
       !
       ! ... local variables
       !
       REAL(DP) :: fac
       INTEGER  :: ibnd, ih, jh, na, nt, ikb, jkb, ijkb0, is, js, ijs !counters
       !
       CALL using_et(0)
       CALL using_indv_ijkb0(0)
       CALL using_becp_auto(0);
       DO ibnd = 1, nbnd
          IF (noncolin) THEN
             CALL compute_deff_nc( deff_nc, et(ibnd,ik) )
          ELSE
             CALL compute_deff( deff, et(ibnd,ik) )
          ENDIF
          !
          fac = wg(ibnd,ik)*tpiba
          !
          DO nt = 1, ntyp
             DO na = 1, nat
                ijkb0 = indv_ijkb0(na)
                IF ( ityp(na) == nt ) THEN
                   DO ih = 1, nh(nt)
                      ikb = ijkb0 + ih
                      IF (noncolin) THEN
                         ijs=0
                         DO is = 1, npol
                            DO js = 1, npol
                               ijs=ijs+1
                               forcenl(ipol,na) = forcenl(ipol,na)- &
                                    deff_nc(ih,ih,na,ijs)*fac*(     &
                                    CONJG(dbecp%nc(ikb,is,ibnd))*   &
                                    becp%nc(ikb,js,ibnd)+           &
                                    CONJG(becp%nc(ikb,is,ibnd))*    &
                                    dbecp%nc(ikb,js,ibnd) )
                            ENDDO
                         ENDDO
                      ELSE
                         forcenl(ipol,na) = forcenl(ipol,na) -   &
                              2.D0 * fac * deff(ih,ih,na)*       &
                              DBLE( CONJG( dbecp%k(ikb,ibnd) ) * &
                              becp%k(ikb,ibnd) )
                      ENDIF
                   ENDDO
                   !
                   IF ( upf(nt)%tvanp .OR. upf(nt)%is_multiproj ) THEN
                      DO ih = 1, nh(nt)
                         ikb = ijkb0 + ih
                         !
                         ! ... in US case there is a contribution for jh<>ih. 
                         ! ... We use here the symmetry in the interchange 
                         ! ... of ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(nt)
                            jkb = ijkb0 + jh
                            IF (noncolin) THEN
                               ijs=0
                               DO is = 1, npol
                                  DO js = 1, npol
                                     ijs = ijs + 1
                                     forcenl(ipol,na) = forcenl(ipol,na)- &
                                          deff_nc(ih,jh,na,ijs)*fac*(     &
                                          CONJG(dbecp%nc(ikb,is,ibnd))*   &
                                          becp%nc(jkb,js,ibnd)+           &
                                          CONJG(becp%nc(ikb,is,ibnd))*    &
                                          dbecp%nc(jkb,js,ibnd))-         &
                                          deff_nc(jh,ih,na,ijs)*fac*(     &
                                          CONJG(dbecp%nc(jkb,is,ibnd))*   &
                                          becp%nc(ikb,js,ibnd)+           &
                                          CONJG(becp%nc(jkb,is,ibnd))*    &
                                          dbecp%nc(ikb,js,ibnd) )
                                  ENDDO
                               ENDDO
                            ELSE
                               forcenl(ipol,na) = forcenl(ipol,na) -     &
                                    2.D0 * fac * deff(ih,jh,na) *        &
                                    DBLE( CONJG( dbecp%k(ikb,ibnd) ) *   &
                                    becp%k(jkb,ibnd) + dbecp%k(jkb,ibnd) &
                                    * CONJG( becp%k(ikb,ibnd) ) )
                            ENDIF
                         ENDDO !jh
                      ENDDO !ih
                   ENDIF ! tvanp
                   !
                ENDIF ! ityp(na) == nt
             ENDDO ! nat
          ENDDO ! ntyp
       ENDDO ! nbnd
       !
       !
     END SUBROUTINE force_us_k_gpu
     !
     !
END SUBROUTINE force_us_gpu
