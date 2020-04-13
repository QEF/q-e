!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stres_us_gpu( ik, gk_d, sigmanlc )
  !----------------------------------------------------------------------------
  !! nonlocal (separable pseudopotential) contribution to the stress
  !! NOTICE: sum of partial results over procs is performed in calling routine
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE constants,            ONLY : eps8
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wvfct,                ONLY : npwx, nbnd, wg, et
  USE control_flags,        ONLY : gamma_only
  USE uspp_param,           ONLY : upf, lmaxkb, nh, nhm
  USE uspp,                 ONLY : nkb, vkb, deeq, deeq_nc
  USE wavefunctions,        ONLY : evc
  USE spin_orb,             ONLY : lspinorb
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : noncolin, npol
  USE mp_pools,             ONLY : me_pool, root_pool
  USE mp_bands,             ONLY : intra_bgrp_comm, me_bgrp, root_bgrp
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   bec_type, becp, calbec
  USE mp,                   ONLY : mp_sum, mp_get_comm_null, mp_circular_shift_left 
  USE wavefunctions_gpum, ONLY : using_evc, using_evc_d, evc_d
  USE wvfct_gpum,                ONLY : using_et
  USE uspp_gpum,                 ONLY : using_vkb, using_deeq
  USE becmod_subs_gpum,          ONLY : using_becp_auto
  !
  !^^^^^^^^^^^^gpu^^^^^^^^^
  !
  !USE becmod,    ONLY :  allocate_bec_type,deallocate_bec_type
  USE becmod_gpum,          ONLY : bec_type_d !, becp_d
  
  USE gbuffers,             ONLY : dev_buf
  USE device_util_m,        ONLY : dev_memcpy
  !^^^^^^^^^^^^^^^
  !
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)    :: ik
  !! k-point index
  REAL(DP), INTENT(IN)    :: gk_d(npwx,3)
  !! wave function components for fixed k-point
  REAL(DP), INTENT(INOUT) :: sigmanlc(3,3)
  !! stress tensor, non-local contribution
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: qm1_d(:)
  REAL(DP) :: q
  INTEGER  :: npw , iu, ierr1
  !
  !INTEGER :: ierr
  INTEGER :: na1, np1, nh_np1, ijkb01, itt
  LOGICAL :: ismulti_np
  INTEGER, ALLOCATABLE :: shift(:)
  INTEGER, ALLOCATABLE :: na_d(:), ih_d(:), ikb_d(:), nhnp_d(:), &
                          ityp_d(:), nh_d(:), ishift_d(:), shift_d(:)
  LOGICAL, ALLOCATABLE :: is_multinp_d(:)
  
  
  !
#if defined(__CUDA)
  attributes(DEVICE) :: gk_d, qm1_d, is_multinp_d, na_d, ishift_d, &
                        ih_d, ikb_d, nhnp_d, shift_d, ityp_d, nh_d
#endif 
  !
  CALL using_evc_d(0)
  CALL using_evc(0)
  !
  IF ( nkb == 0 ) RETURN
  !
  IF ( lsda ) current_spin = isk(ik)
  npw = ngk(ik)
  IF ( nks > 1 ) THEN
    CALL using_vkb(1)
    CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
  ENDIF
  !
  CALL allocate_bec_type( nkb, nbnd, becp, intra_bgrp_comm ) 
  
  CALL using_vkb(0) ; CALL using_becp_auto(2)
  
  CALL calbec( npw, vkb, evc, becp )
  !CALL calbec_gpu ( npw, vkb_d, evc_d, becp_d )

  !
  ALLOCATE( qm1_d( npwx ) )
  !
  !$cuf kernel do (1) <<<*,*>>>
  DO iu = 1, npw
     q = SQRT( gk_d(iu,1)**2 + gk_d(iu,2)**2 + gk_d(iu,3)**2 )
     IF ( q > eps8 ) THEN
        qm1_d(iu) = 1._DP / q
     ELSE
        qm1_d(iu) = 0._DP
     ENDIF
  ENDDO
  !
  !
  ! some atomic index arrays
  ALLOCATE( is_multinp_d(nat*nhm) )
  ALLOCATE( na_d(nat*nhm),  ih_d(nat*nhm)   )
  ALLOCATE( ikb_d(nat*nhm), nhnp_d(nat*nhm) )
  ALLOCATE( ishift_d(nat*nhm), shift_d(nat) )
  ALLOCATE( ityp_d(nat) )
  ALLOCATE( nh_d(ntyp) )
  !
  ityp_d = ityp
  nh_d = nh
  !
  
  ALLOCATE(shift(nat))
  
  ijkb01 = 0 
  DO iu = 1, ntyp
    DO na1 = 1, nat
      IF (ityp(na1) == iu ) THEN
         shift(na1) = ijkb01
         ijkb01 = ijkb01 + nh(iu)
      ENDIF
    ENDDO
  ENDDO
  shift_d = shift
  
  
  ijkb01 = 0
  itt=0
  DO na1 = 1, nat
     np1 = ityp(na1)
     ijkb01 = shift(na1) !ijkb01+nh(np1)
     nh_np1 = nh(np1)
     ismulti_np = upf(np1)%tvanp .OR. upf(np1)%is_multiproj
     !$cuf kernel do (1) <<<*,*>>>
     DO iu = itt+1, itt+nh_np1
       ishift_d(iu) = ijkb01
       nhnp_d(iu)  = nh_np1
       na_d(iu)    = na1
       ih_d(iu)    = iu-itt
       ikb_d(iu)   = ijkb01 + iu - itt
       is_multinp_d(iu) = ismulti_np
     ENDDO 
     itt = itt + nh_np1
  ENDDO
  !
  !
  IF ( gamma_only ) THEN
     !
     CALL stres_us_gamma()
     !
  ELSE
     !
     CALL stres_us_k()
     !
  END IF
  !
  DEALLOCATE( qm1_d )
  
  DEALLOCATE( is_multinp_d  )
  DEALLOCATE( na_d, ih_d    )
  DEALLOCATE( ikb_d, nhnp_d ) 
  DEALLOCATE( ishift_d, shift_d )
  DEALLOCATE( ityp_d )
  DEALLOCATE( nh_d )
  DEALLOCATE( shift )
  
  !
  CALL deallocate_bec_type( becp ) 
  CALL using_becp_auto(2)
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE stres_us_gamma()
       !-----------------------------------------------------------------------
       !! nonlocal contribution to the stress - gamma version
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER                  :: na, np, nt, ibnd, ipol, jpol, l, i, &
                                   ikb, jkb, ih, jh, ibnd_loc,ijkb0,nh_np, &
                                   nproc, nbnd_loc, nbnd_begin, icyc
       REAL(DP)                 :: fac, xyz(3,3), evps, ddot
       !
       !
       !REAL(DP), POINTER        :: deff_d(:,:,:), becpr_d(:,:)  
       COMPLEX(DP), ALLOCATABLE, DEVICE :: deff_d(:,:,:)
       REAL(DP), ALLOCATABLE, DEVICE :: becpr_d(:,:)
       COMPLEX(DP), ALLOCATABLE, DEVICE :: ps_d(:)
       !
       REAL(DP) :: dot11, dot21, dot31, dot22, dot32, dot33
       REAL(DP) :: qm1i, gk1, gk2, gk3
       COMPLEX(DP) :: worksum, cv, wsum1, wsum2, wsum3
       COMPLEX(DP), ALLOCATABLE, DEVICE :: dvkb_d(:,:,:)
       COMPLEX(DP) :: evci
       REAL(DP) :: wg_nk
       !INTEGER :: ierr
       
       integer :: nhmx
       complex(dp) :: becy, defy
       !^^PROVV^^
       !
       COMPLEX(DP), ALLOCATABLE :: dvkb(:,:,:)
       ! dvkb contains the derivatives of the kb potential
       COMPLEX(DP)              :: ps
       ! xyz are the three unit vectors in the x,y,z directions
       DATA xyz / 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 /
       !
       !
       IF( becp%comm /= mp_get_comm_null() ) THEN
          nproc      = becp%nproc
          nbnd_loc   = becp%nbnd_loc
          nbnd_begin = becp%ibnd_begin
          IF( ( nbnd_begin + nbnd_loc - 1 ) > nbnd ) nbnd_loc = nbnd - nbnd_begin + 1
       ELSE
          nproc      = 1
          nbnd_loc   = nbnd
          nbnd_begin = 1
       ENDIF

       ALLOCATE( deff_d(nhm,nhm,nat) )
       !
       ! ... diagonal contribution - if the result from "calbec" are not 
       ! ... distributed, must be calculated on a single processor
       !
       ! ... for the moment when using_gpu is true becp is always fully present in all processors
       !
       CALL using_et(0) ! compute_deff : intent(in)
       !
       ALLOCATE( ps_d(nkb) )
       !
       ALLOCATE( becpr_d(nkb,nbnd_loc) )
       becpr_d = becp%r
       !
       evps = 0.D0
       !
       compute_evps: IF ( .NOT. (nproc == 1 .AND. me_pool /= root_pool) ) THEN
         !
         DO ibnd_loc = 1, nbnd_loc
            ibnd = ibnd_loc + becp%ibnd_begin - 1 
            CALL compute_deff_gpu( deff_d, et(ibnd,ik) )
            wg_nk = wg(ibnd,ik)
            !
            !$cuf kernel do (1) <<<*,*>>>
            DO i = 1, itt
              IF (.NOT. is_multinp_d(i)) THEN
                 !
                 !ih = ih_d(i)
                 !na = na_d(i)
                 !ikb = ikb_d(i)
                 !
                 evps = evps + wg_nk * DBLE(deff_d(ih_d(i),ih_d(i),na_d(i))) * &
                                       ABS(becpr_d(ikb_d(i),ibnd_loc))**2
                 !
              ELSE
                 !
                 evps = evps + wg_nk * DBLE(deff_d(ih_d(i),ih_d(i),na_d(i)))*ABS(becpr_d(ikb_d(i),ibnd_loc))**2 + &
                 becpr_d(ikb_d(i),ibnd_loc)* wg_nk * 2._DP * &
                 SUM( DBLE(deff_d(ih_d(i),ih_d(i)+1:nhnp_d(i),na_d(i))) *  &
                 becpr_d(ishift_d(i)+ih_d(i)+1:ishift_d(i)+nhnp_d(i),ibnd_loc) )
                 !
              ENDIF
            ENDDO
            !
         ENDDO
         !
       ENDIF compute_evps
       !
       !
       ! ... non diagonal contribution - derivative of the bessel function
       !------------------------------------
       ALLOCATE( dvkb( npwx, nkb,4 ) )
       ALLOCATE( dvkb_d(npwx,nkb,4) )
       !
       CALL gen_us_dj( ik, dvkb(:,:,4) )
       IF ( lmaxkb > 0 ) THEN 
         DO ipol = 1, 3
           CALL gen_us_dy( ik, xyz(1,ipol), dvkb(:,:,ipol))
         ENDDO
       ENDIF
       dvkb_d = dvkb
       !
       !
       DO icyc = 0, nproc -1
          !
          DO ibnd_loc = 1, nbnd_loc
             !  
             ibnd = ibnd_loc + becp%ibnd_begin - 1 
             CALL compute_deff_gpu( deff_d, et(ibnd,ik) )
             !
             !$cuf kernel do (1) <<<*,*>>>
             DO i = 1, itt
               !
               IF (.NOT. is_multinp_d(i)) THEN
                  !
                  defy = deff_d(ih_d(i),ih_d(i),na_d(i))
                  becy = CMPLX(becpr_d(ikb_d(i),ibnd_loc))
                  !
                  ps_d(ikb_d(i)) =  becy * defy
                  !
               ELSE
                  !
                  ps_d(ikb_d(i)) = CMPLX(SUM( becpr_d(ishift_d(i)+1:ishift_d(i)+nhnp_d(i),ibnd_loc) * DBLE(deff_d(ih_d(i),1:nhnp_d(i),na_d(i)))))
                  !
               ENDIF
               !
             ENDDO
             !
             !
             dot11 = 0._DP ; dot21 = 0._DP ; dot31 = 0._DP
             dot22 = 0._DP ; dot32 = 0._DP ; dot33 = 0._DP
             !
             !$cuf kernel do(2) <<<*,*>>> 
             DO na =1, nat 
                DO i = 1, npw
                   worksum = (0._DP, 0._DP) 
                   np = ityp_d(na) 
                   ijkb0 = shift_d(na)
                   nh_np = nh_d(np)
                   DO ih = 1, nh_np
                      ikb = ijkb0 + ih  
                      worksum = worksum + ps_d(ikb) * dvkb_d(i,ikb,4) 
                   ENDDO
                   evci = evc_d(i,ibnd) 
                   gk1  = gk_d(i,1) 
                   gk2  = gk_d(i,2) 
                   gk3  = gk_d(i,3) 
                   qm1i = qm1_d(i) 
                   !  
                   cv = evci * CMPLX(gk1 * gk1  * qm1i)
                   dot11 = dot11 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv) 
                   !  
                   cv = evci * CMPLX(gk2 * gk1 * qm1i )
                   dot21 = dot21 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                   ! 
                   cv = evci * CMPLX(gk3 * gk1 * qm1i)
                   dot31 = dot31 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                   ! 
                   cv = evci * CMPLX(gk2 * gk2 * qm1i)
                   dot22 = dot22 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                   ! 
                   cv = evci * CMPLX(gk3 * gk2 * qm1i)
                   dot32 = dot32 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                   !
                   cv = evci * CMPLX(gk3 * gk3 * qm1i)
                   dot33 = dot33 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                END DO
             END DO
             ! ... a factor 2 accounts for the other half of the G-vector sphere
             sigmanlc(:,1) = sigmanlc(:,1) - 4._DP * wg(ibnd,ik) * [dot11, dot21, dot31] 
             sigmanlc(:,2) = sigmanlc(:,2) - 4._DP * wg(ibnd,ik) * [0._DP, dot22, dot32]
             sigmanlc(:,3) = sigmanlc(:,3) - 4._DP * wg(ibnd,ik) * [0._DP, 0._DP, dot33]                  

             !
             ! ... non diagonal contribution - derivative of the spherical harmonics
             ! ... (no contribution from l=0)
             !
             IF ( lmaxkb == 0 ) CYCLE 
             !
             !------------------------------------
             !
             dot11 = 0._DP ; dot21 = 0._DP ; dot31 = 0._DP
             dot22 = 0._DP ; dot32 = 0._DP ; dot33 = 0._DP
             !
             !$cuf kernel do(2) <<<*,*>>> 
             DO ikb = 1, nkb 
                DO i = 1, npw  
                   !
                   wsum1 = ps_d(ikb)*dvkb_d(i,ikb,1) 
                   wsum2 = ps_d(ikb)*dvkb_d(i,ikb,2) 
                   wsum3 = ps_d(ikb)*dvkb_d(i,ikb,3)      
                   !
                   evci = evc_d(i,ibnd) 
                   gk1 = gk_d(i,1)
                   gk2 = gk_d(i,2)
                   gk3 = gk_d(i,3) 
                   !
                   cv = evci * CMPLX(gk1 )
                   dot11 = dot11 + DBLE( wsum1)* DBLE(cv) + DIMAG(wsum1)*DIMAG(cv)
                   dot21 = dot21 + DBLE( wsum2)* DBLE(cv) + DIMAG(wsum2)*DIMAG(cv)
                   dot31 = dot31 + DBLE( wsum3)* DBLE(cv) + DIMAG(wsum3)*DIMAG(cv)
                   !
                   cv = evci * CMPLX( gk2)
                   dot22 = dot22 + DBLE( wsum2)* DBLE(cv) + DIMAG(wsum2)*DIMAG(cv) 
                   dot32 = dot32 + DBLE( wsum3)* DBLE(cv) + DIMAG(wsum3)*DIMAG(cv)
                   ! 
                   cv =  evci * CMPLX( gk3 )
                   dot33 = dot33 + DBLE( wsum3)* DBLE(cv) + DIMAG(wsum3)*DIMAG(cv)
                ENDDO
             ENDDO 
             !
             ! ... a factor 2 accounts for the other half of the G-vector sphere
             !
             sigmanlc(:,1) = sigmanlc(:,1) -4._DP * wg(ibnd, ik) * [dot11, dot21, dot31]
             sigmanlc(:,2) = sigmanlc(:,2) -4._DP * wg(ibnd, ik) * [0._DP, dot22, dot32]
             sigmanlc(:,3) = sigmanlc(:,3) -4._DP * wg(ibnd, ik) * [0._DP, 0._DP, dot33]
             IF ( nproc > 1 ) THEN
                 CALL errore ('stres_us_gamma_gpu line 303', &
                       'unexpected error nproc be 1 with GPU acceleration', 100) 
                !CALL mp_circular_shift_left(becp%r, icyc, becp%comm)
                !CALL mp_circular_shift_left(becp%ibnd_begin, icyc, becp%comm)
                !CALL mp_circular_shift_left(nbnd_loc, icyc, becp%comm)
             ENDIF
          ENDDO
       ENDDO
       !
10     CONTINUE
       !
       DO l = 1, 3
          sigmanlc(l,l) = sigmanlc(l,l) - evps
       ENDDO
       !
       !
       DEALLOCATE( dvkb )
       DEALLOCATE( deff_d  )
       DEALLOCATE( ps_d    )
       DEALLOCATE( becpr_d )
       DEALLOCATE( dvkb_d  )
       !
       !
       RETURN
       !
     END SUBROUTINE stres_us_gamma     
     !
     !
     !----------------------------------------------------------------------
     SUBROUTINE stres_us_k()
       !----------------------------------------------------------------------  
       !! nonlocal contribution to the stress - k-points version       
       !
#if defined(_OPENMP) && defined(__PGI)
       USE omp_lib
#endif
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: na, np, ibnd, ipol, jpol, l, i, ipw, &
                  ikb, jkb, ih, jh, is, js, ijs
       REAL(DP) :: fac, xyz(3,3), evps, ddot
       COMPLEX(DP), ALLOCATABLE :: dvkb(:,:,:)
       !
       
       
       !^^^^^^^^^^^^gpu
       COMPLEX(DP) :: qm1i
       COMPLEX(DP) :: gk1, gk2, gk3
       !
       REAL(DP) :: dot11, dot21, dot31, dot22, dot32, dot33
       !
       COMPLEX(DP) :: cv, cv1, cv2, worksum, worksum1, worksum2, evci, evc1i, evc2i, & 
                      ps1, ps2, ps1d1, ps1d2, ps1d3, ps2d1, ps2d2, ps2d3, psd1, psd2, psd3
       
       INTEGER :: nt, npol2
       INTEGER :: ierr
       
       !TYPE(bec_type_d), TARGET :: becp_d 
       
       INTEGER :: nh_np
       REAL(DP), ALLOCATABLE, DEVICE :: deeq_d(:,:,:,:)
       COMPLEX(DP), ALLOCATABLE, DEVICE :: dvkb_d(:,:,:)
       COMPLEX(DP), ALLOCATABLE, DEVICE :: deff_d(:,:,:), deff_nc_d(:,:,:,:)
       
       COMPLEX(DP), POINTER :: ps_d(:)
       COMPLEX(DP), ALLOCATABLE :: ps_nc_d(:,:) !metti pointer
       
       !COMPLEX(DP), POINTER ::  becpnc_d(:,:,:)
       COMPLEX(DP), ALLOCATABLE, DEVICE :: becpk_d(:,:)
       COMPLEX(DP), ALLOCATABLE, DEVICE :: becpnc_d(:,:)
       INTEGER :: nhmx !, npwt
       complex(dp) :: becy, deqy
       !^^PROVV^^
       
#if defined(__CUDA) 
  ATTRIBUTES(DEVICE) :: ps_d, ps_nc_d
#endif
       
       !^^^^^^^^^^^^^
       
       COMPLEX(DP), ALLOCATABLE :: deff_nc(:,:,:,:)
       ! dvkb contains the derivatives of the kb potential
       COMPLEX(DP) :: ps, ps_nc(2)
       ! xyz are the three unit vectors in the x,y,z directions
       DATA xyz / 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 /
       !
       ! WORKAROUND STARTS ==================================================
       !
       ! There seems to be a bug with the OpenMP code generated by PGI 18.5
       ! for this subroutine.
       !
#if defined(_OPENMP) && defined(__PGI)
       INTEGER :: num_threads
       num_threads=omp_get_max_threads()
       CALL omp_set_num_threads(1)
#endif
       ! WORKAROUND ENDS ====================================================
       !
       !
       IF (noncolin) ALLOCATE( deff_nc(nhm,nhm,nat,nspin) )
       !
       !
       evps = 0.D0
       ! ... diagonal contribution
       !
       !
       ALLOCATE( dvkb(npwx,nkb,4) )
       !^^^^
       !CALL gen_us_dj( ik, dvkb )
       !
       !^^^^^^^^^^^^^^^^^^^^^^
       
       ALLOCATE( dvkb_d(npwx,nkb,4) )
       CALL gen_us_dj( ik, dvkb(:,:,4) )
       IF ( lmaxkb > 0 ) THEN 
         DO ipol = 1, 3
           CALL gen_us_dy( ik, xyz(1,ipol), dvkb(:,:,ipol))
         ENDDO
       ENDIF
       dvkb_d = dvkb
       !
       !
       IF (noncolin) THEN 
          !CALL dev_buf%lock_buffer( ps_nc_d, [nkb, npol], ierr )
          ALLOCATE( ps_nc_d(nkb,npol) )
          
          !becpnc_d => becp_d%nc_d
          ALLOCATE( becpnc_d(nkb,npol) )
          !becpnc_d = becp%nc
          !CALL dev_buf%lock_buffer( becpnc_d, [nkb, npol], ierr )
          !CALL dev_memcpy( becpnc_d, becp%nc )
          !
          ALLOCATE( deff_nc_d(nhm,nhm,nat,nspin) )
       ELSE 
          CALL dev_buf%lock_buffer ( ps_d, nkb, ierr)
          !becpk_d => becp_d%k_d
          ALLOCATE( becpk_d(nkb,nbnd) )
          becpk_d = becp%k 
          
          ALLOCATE( deff_d(nhm,nhm,nat) )
       END IF
       !
       
       !
       CALL using_et(0) ! this is redundant
       CALL using_deeq(0)
       
       !
       !^^^^^^^^^^^^^^^^^^^^^^^^ 
       nhmx = SIZE(deeq(:,1,1,1))
       ALLOCATE( deeq_d(nhmx,nhmx,nat,nspin) )
       deeq_d = deeq
       !
       !
       CALL start_clock('knl_ciclo2')
       !
       CALL using_et(0) ! compute_deff : intent(in)
       
       IF ( me_bgrp /= root_bgrp ) GO TO 100
       !
       ! ... the contribution is calculated only on one processor because
       ! ... partial results are later summed over all processors
       !
       DO ibnd = 1, nbnd
          fac = wg(ibnd,ik)
          IF (ABS(fac) < 1.d-9) CYCLE
          IF (noncolin) THEN
             CALL compute_deff_nc(deff_nc,et(ibnd,ik))
             deff_nc_d = deff_nc
             becpnc_d = becp%nc(:,:,ibnd)
             
          ELSE
             CALL compute_deff_gpu( deff_d, et(ibnd,ik) )
          ENDIF
          !
          !
          IF (noncolin) THEN
            !
            !$cuf kernel do (1) <<<*,*>>>
            DO i = 1, itt
              IF (.NOT. is_multinp_d(i)) THEN
                 !
                 ijs = 0
                 DO is = 1, npol
                   DO js = 1, npol
                      ijs = ijs + 1
                      evps=evps+fac*DBLE(deff_nc_d(ih_d(i),ih_d(i),na_d(i),ijs)* &
                                 CONJG(becpnc_d(ikb_d(i),is))* &
                                      becpnc_d(ikb_d(i),js))
                   ENDDO
                 ENDDO
                 !
              ELSE
                 !
                 ijs=0
                 DO is=1,npol
                   DO js=1,npol
                      ijs=ijs+1
                      evps=evps+fac*DBLE(deff_nc_d(ih_d(i),ih_d(i),na_d(i),ijs)*   &
                                 CONJG(becpnc_d(ikb_d(i),is))* &
                                      becpnc_d(ikb_d(i),js))
                   ENDDO
                 ENDDO
                 !
                 !
                 DO jh = ( ih_d(i) + 1 ), nhnp_d(i)
                    jkb = ishift_d(i) + jh
                    ijs=0
                    DO is=1,npol
                      DO js=1,npol
                         ijs=ijs+1
                         evps = evps + &
                                +2._DP*fac &
                                *DBLE(deff_nc_d(ih_d(i),jh,na_d(i),ijs) *    &
                                (CONJG( becpnc_d(ikb_d(i),is) ) * &
                                        becpnc_d(jkb,js)) )
                      ENDDO
                    ENDDO
                 ENDDO
                 !
              ENDIF
            ENDDO
            !

            !
          ELSE
            !
            !$cuf kernel do (1) <<<*,*>>>
            DO i = 1, itt
              IF (.NOT. is_multinp_d(i)) THEN
                 !
                 evps = evps+fac*deff_d(ih_d(i),ih_d(i),na_d(i))*ABS(becpk_d(ikb_d(i),ibnd) )**2
                 !
              ELSE
                 !
                 evps = evps +fac*DBLE(deff_d(ih_d(i),ih_d(i),na_d(i))*ABS(becpk_d(ikb_d(i),ibnd) )**2) + &
                 DBLE(SUM( deff_d(ih_d(i),ih_d(i)+1:nhnp_d(i),na_d(i)) * fac * 2._DP * &
                               DBLE( CONJG( becpk_d(ikb_d(i),ibnd) ) * &
                               becpk_d(ishift_d(i)+ih_d(i)+1:ishift_d(i)+nhnp_d(i),ibnd) ) ))
                 !
              ENDIF
            ENDDO
            !
          ENDIF
          !
       END DO
       
       DO l = 1, 3
          sigmanlc(l,l) = sigmanlc(l,l) - evps
       END DO
       !
       
       
       CALL stop_clock('knl_ciclo2')
       
       
100    CONTINUE
       !
       ! ... non diagonal contribution - derivative of the bessel function
       !
       !
       CALL start_clock('knl_ciclo2')
       !
       
       !call using_becp_auto(0)
       !
       DO ibnd = 1, nbnd
         !
         
         IF ( noncolin ) THEN
            !
            !^^^^^^^^^^no gpu^^ 
            CALL compute_deff_nc( deff_nc, et(ibnd,ik) )
            deff_nc_d = deff_nc
            !^^
            !CALL compute_deff_nc_gpu(deff_nc_d,et_d(ibnd,ik))
            !^^^^^^^^^^^
            !
            becpnc_d = becp%nc(:,:,ibnd)
            
            
            !
            !$cuf kernel do (1) <<<*,*>>>
            DO i = 1, itt
              IF (.NOT. is_multinp_d(i)) THEN
                 !
                 ijs=0
                 ps_nc_d(ikb_d(i),:) = (0._DP,0._DP)
                 DO is=1,npol
                   DO js=1,npol
                      ijs=ijs+1

                      ps_nc_d(ikb_d(i),is) = ps_nc_d(ikb_d(i),is) + becpnc_d(ikb_d(i),js)* &
                                      deff_nc_d(ih_d(i),ih_d(i),na_d(i),ijs)
                   ENDDO
                 ENDDO
                 !
              ELSE
                 !
                 ps_nc_d(ikb_d(i),:) = (0._DP,0._DP) 
                 DO jh = 1, nhnp_d(i)
                    jkb = ishift_d(i) + jh
                    ijs=0
                    DO is=1,npol
                      DO js=1,npol
                        ijs=ijs+1

                        ps_nc_d(ikb_d(i),is)= ps_nc_d(ikb_d(i),is) + becpnc_d(jkb,js)* &
                                        deff_nc_d(ih_d(i),jh,na_d(i),ijs)
                      ENDDO
                    ENDDO
                 ENDDO   
                 !
              ENDIF
            ENDDO
            
            !
         ELSE
            !
            !^^^^^^
            CALL compute_deff_gpu( deff_d, et(ibnd,ik) )
            !
            !$cuf kernel do (1) <<<*,*>>>
            DO i = 1, itt
               !
               IF (.NOT. is_multinp_d(i)) THEN
                   !
                   deqy = CMPLX(deeq_d(ih_d(i),ih_d(i),na_d(i),current_spin))
                   becy = becpk_d(ikb_d(i),ibnd)
                   !
                   ps_d(ikb_d(i)) =  becy * deqy
                   !
               ELSE 
                   !
                   ps_d(ikb_d(i)) = SUM( becpk_d(ishift_d(i)+1:ishift_d(i)+nhnp_d(i),ibnd) * deff_d(ih_d(i),1:nhnp_d(i),na_d(i)) )
                   !
               ENDIF
               !
            ENDDO
            !
         ENDIF
         !
         dot11 = 0._DP ; dot21 = 0._DP ; dot31 = 0._DP
         dot22 = 0._DP ; dot32 = 0._DP ; dot33 = 0._DP
         !
         !
         IF (noncolin) THEN
            !$cuf kernel do(2) <<<*,*>>>    
            DO ikb =1, nkb
               DO i = 1, npw    
                  evc1i = evc_d(i, ibnd)
                  evc2i = evc_d(i+npwx,ibnd)
                  qm1i = CMPLX(qm1_d(i))
                  gk1 = CMPLX(gk_d(i,1))
                  gk2 = CMPLX(gk_d(i,2))
                  gk3 = CMPLX(gk_d(i,3))
                  worksum1 = ps_nc_d(ikb,1) * dvkb_d(i,ikb,4)     
                  worksum2 = ps_nc_d(ikb,2) * dvkb_d(i,ikb,4)   
                  !   
                  cv1 = evc1i * gk1 * gk1 * qm1i
                  cv2 = evc2i * gk1 * gk1 * qm1i
                  dot11 = dot11 + DBLE(worksum1)*DBLE(cv1)   + &
                                  DIMAG(worksum1)*DIMAG(cv1) + &
                                  DBLE(worksum2)*DBLE(cv2)   + &
                                  DIMAG(worksum2)*DIMAG(cv2)
                  !   
                  cv1 = evc1i * gk2 * gk1 * qm1i
                  cv2 = evc2i * gk2 * gk1 * qm1i
                  dot21 = dot21 + DBLE(worksum1)*DBLE(cv1)   + &
                                  DIMAG(worksum1)*DIMAG(cv1) + &
                                  DBLE(worksum2)*DBLE(cv2)   + &
                                  DIMAG(worksum2)*DIMAG(cv2)   
                  !   
                  cv1 = evc1i * gk3 * gk1 * qm1i
                  cv2 = evc2i * gk3 * gk1 * qm1i
                  dot31 = dot31 + DBLE(worksum1)*DBLE(cv1)   + &
                                  DIMAG(worksum1)*DIMAG(cv1) + &
                                  DBLE(worksum2)*DBLE(cv2)   + &
                                  DIMAG(worksum2)*DIMAG(cv2)   
                  !   
                  cv1 = evc1i * gk2 * gk2 * qm1i
                  cv2 = evc2i * gk2 * gk2 * qm1i
                  dot22 = dot22 + DBLE(worksum1)*DBLE(cv1)   + &
                                  DIMAG(worksum1)*DIMAG(cv1) + &
                                  DBLE(worksum2)*DBLE(cv2)   + &
                                  DIMAG(worksum2)*DIMAG(cv2)   
                  !   
                  cv1 = evc1i * gk3 * gk2 * qm1i
                  cv2 = evc2i * gk3 * gk2 * qm1i
                  dot32 = dot32 + DBLE(worksum1)*DBLE(cv1)   + &
                                  DIMAG(worksum1)*DIMAG(cv1) + &
                                  DBLE(worksum2)*DBLE(cv2)   + &
                                  DIMAG(worksum2)*DIMAG(cv2)
                  !   
                  cv1  = evc1i * gk3 * gk3 * qm1i
                  cv2  = evc2i * gk3 * gk3 * qm1i
                  dot33 = dot33 + DBLE(worksum1)*DBLE(cv1)   + &
                                  DIMAG(worksum1)*DIMAG(cv1) + &
                                  DBLE(worksum2)*DBLE(cv2)   + &
                                  DIMAG(worksum2)*DIMAG(cv2)
               END DO
            END DO
            !
         ELSE   
            !
            !
            !$cuf kernel do(2) <<<*,*>>>   
            DO ikb = 1, nkb
               DO i = 1, npw
                  !
                  worksum = ps_d(ikb) *dvkb_d(i,ikb,4)   
                  !    
                  evci = evc_d(i,ibnd)
                  qm1i = CMPLX(qm1_d(i))
                  gk1 = CMPLX(gk_d(i,1))
                  gk2 = CMPLX(gk_d(i,2))
                  gk3 = CMPLX(gk_d(i,3))
                  !
                  cv = evci * gk1 * gk1 * qm1i
                  dot11 = dot11 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                  !
                  cv = evci * gk2 * gk1 * qm1i
                  dot21 = dot21 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                  !
                  cv = evci * gk3 * gk1 * qm1i
                  dot31 = dot31 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                  !
                  cv = evci * gk2 * gk2 * qm1i
                  dot22 = dot22 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                  !
                  cv = evci * gk3 * gk2 * qm1i
                  dot32 = dot32 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                  !
                  cv = evci * gk3 * gk3 * qm1i
                  dot33 = dot33 + DBLE(worksum)*DBLE(cv) + DIMAG(worksum)*DIMAG(cv)
                  !   
               ENDDO   
            ENDDO   
            !   
            !IF ( me_bgrp /= root_bgrp ) stop     
            !   
         ENDIF
         !
         sigmanlc(:,1) = sigmanlc(:,1) - 2._DP * wg(ibnd,ik) * [dot11, dot21, dot31]
         sigmanlc(:,2) = sigmanlc(:,2) - 2._DP * wg(ibnd,ik) * [0._DP, dot22, dot32]
         sigmanlc(:,3) = sigmanlc(:,3) - 2._DP * wg(ibnd,ik) * [0._DP, 0._DP, dot33]
         !      
         ! ... non diagonal contribution - derivative of the spherical harmonics
         ! ... (no contribution from l=0)
         !      
         IF ( lmaxkb == 0 ) CYCLE       
         !      
         dot11 = 0._DP ;  dot21 = 0._DP      
         dot31 = 0._DP ;  dot22 = 0._DP      
         dot32 = 0._DP ;  dot33 = 0._DP      
         !
         IF (noncolin) THEN      
            !      
            !$cuf kernel do(2) <<<*,*>>>      
            DO ikb =1, nkb
               DO i = 1, npw
                  !       
                  gk1 = CMPLX(gk_d(i,1))
                  gk2 = CMPLX(gk_d(i,2))
                  gk3 = CMPLX(gk_d(i,3))
                  !
                  ps1 = ps_nc_d(ikb,1)
                  ps2 = ps_nc_d(ikb,2)
                  !
                  ps1d1 = ps1 * dvkb_d(i,ikb,1)
                  ps1d2 = ps1 * dvkb_d(i,ikb,2)       
                  ps1d3 = ps1 * dvkb_d(i,ikb,3)       
                  !
                  ps2d1 = ps2 * dvkb_d(i,ikb,1)
                  ps2d2 = ps2 * dvkb_d(i,ikb,2)
                  ps2d3 = ps2 * dvkb_d(i,ikb,3)
                  !
                  evc1i = evc_d(i,ibnd)       
                  evc2i = evc_d(i+npwx,ibnd)       
                  !      
                  cv1 = evc1i * gk1
                  cv2 = evc2i * gk1
                  dot11 = dot11 + DBLE(ps1d1)*DBLE(cv1)   + &
                                  DIMAG(ps1d1)*DIMAG(cv1) + &
                                  DBLE(ps2d1)*DBLE(cv2)   + &
                                  DIMAG(ps2d1)*DIMAG(cv2)
                  !
                  dot21 = dot21 + DBLE(ps1d2)*DBLE(cv1)   + &
                                  DIMAG(ps1d2)*DIMAG(cv1) + &
                                  DBLE(ps2d2)*DBLE(cv2)   + &
                                  DIMAG(ps2d2)*DIMAG(cv2)
                  !
                  dot31 = dot31 + DBLE(ps1d3)*DBLE(cv1)   + &
                                  DIMAG(ps1d3)*DIMAG(cv1) + &
                                  DBLE(ps2d3)*DBLE(cv2)   + &
                                  DIMAG(ps2d3)*DIMAG(cv2)       
                  !
                  cv1 = evc1i * gk2
                  cv2 = evc2i * gk2
                  dot22 = dot22 + DBLE(ps1d2)*DBLE(cv1)   + &
                                  DIMAG(ps1d2)*DIMAG(cv1) + &
                                  DBLE(ps2d2)*DBLE(cv2)   + &
                                  DIMAG(ps2d2)*DIMAG(cv2)
                  !
                  dot32 = dot32 + DBLE(ps1d3)*DBLE(cv1)   + &
                                  DIMAG(ps1d3)*DIMAG(cv1) + &
                                  DBLE(ps2d3)*DBLE(cv2)   + &
                                  DIMAG(ps2d3)*DIMAG(cv2)
                  !
                  cv1 = evc1i * gk3
                  cv2 = evc2i * gk3
                  dot33 = dot33 + DBLE(ps1d3)*DBLE(cv1)   + &
                                  DIMAG(ps1d3)*DIMAG(cv1) + &
                                  DBLE(ps2d3)*DBLE(cv2)   + &
                                  DIMAG(ps2d3)*DIMAG(cv2)      
                  !
               END DO
            END DO
            !
         ELSE
            !
            !$cuf kernel do(2) <<<*,*>>>
            DO ikb = 1, nkb
               DO i = 1, npw
                 ps   = ps_d(ikb)
                 psd1 = ps*dvkb_d(i,ikb,1)
                 psd2 = ps*dvkb_d(i,ikb,2)       
                 psd3 = ps*dvkb_d(i,ikb,3)
                 evci = evc_d(i,ibnd)
                 gk1  = CMPLX(gk_d(i,1))
                 gk2  = CMPLX(gk_d(i,2))
                 gk3  = CMPLX(gk_d(i,3))
                 !
                 cv = evci * gk1
                 dot11 = dot11 + DBLE(psd1)*DBLE(cv) + DIMAG(psd1)*DIMAG(cv)
                 dot21 = dot21 + DBLE(psd2)*DBLE(cv) + DIMAG(psd2)*DIMAG(cv)
                 dot31 = dot31 + DBLE(psd3)*DBLE(cv) + DIMAG(psd3)*DIMAG(cv)
                 !
                 cv = evci * gk2
                 dot22 = dot22 + DBLE(psd2)*DBLE(cv) + DIMAG(psd2)*DIMAG(cv)
                 dot32 = dot32 + DBLE(psd3)*DBLE(cv) + DIMAG(psd3)*DIMAG(cv)
                 !
                 cv = evci * gk3
                 dot33 = dot33 + DBLE(psd3)*DBLE(cv) + DIMAG(psd3)*DIMAG(cv)
               ENDDO
            ENDDO
            !   
         ENDIF      
         !       
         sigmanlc(:,1) = sigmanlc(:,1) -2._DP * wg(ibnd,ik) * [dot11, dot21, dot31]
         sigmanlc(:,2) = sigmanlc(:,2) -2._DP * wg(ibnd,ik) * [0._DP, dot22, dot32]
         sigmanlc(:,3) = sigmanlc(:,3) -2._DP * wg(ibnd,ik) * [0._DP, 0._DP, dot33]
         !
       ENDDO 
       !
10     CONTINUE
       !
       DEALLOCATE( dvkb )
       DEALLOCATE( dvkb_d )
       DEALLOCATE( deeq_d )
       IF (noncolin) THEN 
          DEALLOCATE( ps_nc_d )
          DEALLOCATE( becpnc_d )
          DEALLOCATE( deff_nc )
          DEALLOCATE( deff_nc_d )
       ELSE 
          CALL dev_buf%release_buffer( ps_d, ierr )
          DEALLOCATE( becpk_d )
          DEALLOCATE( deff_d )
       ENDIF
       !
       ! WORKAROUND STARTS ==================================================
       !
       ! ... and now restore the previous value.
       !
#if defined(_OPENMP) && defined(__PGI)
       CALL omp_set_num_threads(num_threads)
#endif
       ! WORKAROUND ENDS ====================================================
       !
       RETURN
       !
     END SUBROUTINE stres_us_k
     !
END SUBROUTINE stres_us_gpu
