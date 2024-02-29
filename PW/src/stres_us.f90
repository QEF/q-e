!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stres_us( ik, gk, sigmanlc )
  !----------------------------------------------------------------------------
  !! Nonlocal (separable pseudopotential) contribution to the stress
  !! NOTICE: sum of partial results over procs is performed in calling routine
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE constants,            ONLY : eps8
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wvfct,                ONLY : npwx, nbnd, wg, et
  USE control_flags,        ONLY : gamma_only, offload_type
  USE uspp_param,           ONLY : upf, lmaxkb, nh, nhm
  USE uspp,                 ONLY : nkb, vkb, deeq
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : noncolin, npol, lspinorb
  USE mp_pools,             ONLY : me_pool, root_pool
  USE mp_bands,             ONLY : intra_bgrp_comm, me_bgrp, root_bgrp
  USE becmod,               ONLY : allocate_bec_type_acc, deallocate_bec_type_acc, &
                                   bec_type, becp, calbec
  USE mp,                   ONLY : mp_sum, mp_get_comm_null, &
                                   mp_circular_shift_left 
  USE wavefunctions,        ONLY : evc
  USE uspp_init,            ONLY : init_us_2, gen_us_dj, gen_us_dy
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: ik
  !! k-point index
  REAL(DP), INTENT(IN) :: gk(npwx,3)
  !! wave function components for fixed k-point
  REAL(DP), INTENT(INOUT) :: sigmanlc(3,3)
  !! stress tensor, non-local contribution
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: qm1(:)
  COMPLEX(DP), ALLOCATABLE :: evcv(:)
  !
  REAL(DP) :: q
  INTEGER  :: npw , iu, np
  !
  INTEGER :: na1, np1, nh_np1, ijkb01, itot, levc
  LOGICAL :: ismulti_np
  INTEGER, ALLOCATABLE :: shift(:), na_list(:), nh_list(:), ih_list(:), &
                          ishift_list(:)
  LOGICAL, ALLOCATABLE :: is_multinp(:)
  !
  IF ( nkb == 0 ) RETURN
  !
  IF ( lsda ) current_spin = isk(ik)
  npw = ngk(ik)
  !
  IF ( nks > 1 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .TRUE. )
  !
  CALL allocate_bec_type_acc( nkb, nbnd, becp, intra_bgrp_comm )
  !$acc data present( vkb, evc, becp )
  CALL calbec( offload_type, npw, vkb, evc, becp )
  !$acc end data
  !
  ALLOCATE( qm1(npwx) )
  !$acc data create(qm1) present(gk)
  !
  !$acc parallel loop
  DO iu = 1, npw
     q = SQRT( gk(iu,1)**2 + gk(iu,2)**2 + gk(iu,3)**2 )
     IF ( q > eps8 ) THEN
        qm1(iu) = 1._DP / q
     ELSE
        qm1(iu) = 0._DP
     ENDIF
  ENDDO
  !
  ! ... define index arrays (type, atom, etc.) for kernel loops
  !
  ALLOCATE( is_multinp(nat*nhm) )
  ALLOCATE( na_list(nat*nhm), nh_list(nat*nhm), ih_list(nat*nhm), &
            ishift_list(nat*nhm) )
  !$acc data create(is_multinp,na_list,nh_list,ih_list,ishift_list) &
  !$acc&     copyin(ityp,nh)
  !
  ALLOCATE( shift(nat) )
  ijkb01 = 0
  DO iu = 1, ntyp
    DO na1 = 1, nat
      IF (ityp(na1) == iu ) THEN
         shift(na1) = ijkb01
         ijkb01 = ijkb01 + nh(iu)
      ENDIF
    ENDDO
  ENDDO
  !
  !$acc data copyin(shift)
  !
  ijkb01 = 0
  itot=0
  DO np = 1, ntyp
    DO na1 = 1, nat
      np1 = ityp(na1)
      IF (np /= np1) CYCLE
      ijkb01 = shift(na1)
      nh_np1 = nh(np1)
      ismulti_np = upf(np1)%tvanp .OR. upf(np1)%is_multiproj
      IF ( .NOT. ismulti_np .AND. noncolin .AND. lspinorb ) &
                                   CALL errore('stres_us','wrong case',1)
      !$acc parallel loop
      DO iu = itot+1, itot+nh_np1
        na_list(iu) = na1
        ih_list(iu) = iu-itot
        nh_list(iu) = nh_np1
        ishift_list(iu) = ijkb01
        is_multinp(iu) = ismulti_np
      ENDDO
      itot = itot + nh_np1
    ENDDO
  ENDDO
  !
  levc = SIZE(evc(:,1))
  ALLOCATE( evcv(1:levc) )
  !$acc data create(evcv)
  !
  IF ( gamma_only ) THEN
     !
     CALL stres_us_gamma()
     !
  ELSE
     !
     CALL stres_us_k()
     !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( evcv )
  !
  !$acc end data
  !$acc end data
  !$acc end data
  !
  DEALLOCATE( qm1 )
  DEALLOCATE( is_multinp )
  DEALLOCATE( na_list, nh_list, ih_list, ishift_list )
  DEALLOCATE( shift )
  !
  CALL deallocate_bec_type_acc( becp ) 
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE stres_us_gamma()
       !-----------------------------------------------------------------------
       !! nonlocal contribution to the stress - gamma version.
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER  :: na, np, nt, ibnd, ipol, jpol, l, i, ikb,  &
                   jkb, ih, jh, ibnd_loc,ijkb0,nh_np, nproc, &
                   nbnd_loc, nbnd_begin, icyc, ishift
       REAL(DP) :: dot11, dot21, dot31, dot22, dot32, dot33,  &
                   qm1i, gk1, gk2, gk3, wg_nk, fac, evps, aux,&
                   Re_worksum, Im_worksum
       COMPLEX(DP) :: worksum, cv, wsum1, wsum2, wsum3, evci
       !
       REAL(DP), ALLOCATABLE :: deff(:,:,:)
       COMPLEX(DP), ALLOCATABLE :: ps(:), dvkb(:,:,:)
       !
       REAL(DP) :: xyz(3,3)
       !
       ! xyz are the three unit vectors in the x,y,z directions
       DATA xyz / 1._DP, 0._DP, 0._DP, 0._DP, 1._DP, 0._DP, 0._DP, 0._DP, &
                  1._DP /
       REAL(DP), ALLOCATABLE :: becpr(:,:)
       !$acc declare device_resident(becpr)
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
       !
       ALLOCATE( deff(nhm,nhm,nat) )
       ALLOCATE( ps(nkb) )
       !$acc data create(deff,ps)
       !
       ! ... diagonal contribution - if the result from "calbec" are not 
       ! ... distributed, must be calculated on a single processor
       !
       ! ... for the moment when using_gpu is true becp is always fully present
       !     in all processors
       !
       ALLOCATE( becpr(nkb,nbnd_loc) )
       !$acc kernels present(becp%r)
       becpr = becp%r
       !$acc end kernels
       !
       evps = 0._DP
       !
       compute_evps: IF ( .NOT. (nproc==1 .AND. me_pool/=root_pool) ) THEN
         !
         DO ibnd_loc = 1, nbnd_loc
            ibnd = ibnd_loc+becp%ibnd_begin-1
            CALL compute_deff( deff, et(ibnd,ik) )
            wg_nk = wg(ibnd,ik)
            !
            !$acc parallel loop reduction(+:evps)
            DO i = 1, itot
              ih = ih_list(i)           ;   na = na_list(i)
              ishift = ishift_list(i)   ;   ikb = ishift + ih
              !
              IF (.NOT. is_multinp(i)) THEN
                 aux = wg_nk * deff(ih,ih,na) * ABS(becpr(ikb,ibnd_loc))**2
              ELSE
                 nh_np = nh_list(i)
                 !
                 aux = wg_nk * deff(ih,ih,na) * ABS(becpr(ikb,ibnd_loc))**2 &
                             +  becpr(ikb,ibnd_loc)* wg_nk * 2._DP  &
                                * SUM( deff(ih,ih+1:nh_np,na)       &
                                * becpr(ishift+ih+1:ishift+nh_np,ibnd_loc))
              ENDIF
              !
              evps = evps + aux
              !
            ENDDO
            !
         ENDDO
         !
       ENDIF compute_evps
       !
       ! ... non diagonal contribution - derivative of the Bessel function
       !
       ALLOCATE( dvkb(npwx,nkb,4) )
       !$acc data create(dvkb)
       !
       CALL gen_us_dj( ik, dvkb(:,:,4) )
       IF ( lmaxkb > 0 ) THEN
         DO ipol = 1, 3
           CALL gen_us_dy( ik, xyz(1,ipol), dvkb(:,:,ipol) )
         ENDDO
       ENDIF
       !
       DO icyc = 0, nproc-1
          !
          DO ibnd_loc = 1, nbnd_loc
             !
             ibnd = ibnd_loc + becp%ibnd_begin - 1
             CALL compute_deff( deff, et(ibnd,ik) )
             !
             !$acc kernels
             evcv(:) = evc(:,ibnd)
             !$acc end kernels
             !
             !$acc parallel loop
             DO i = 1, itot
               ih = ih_list(i)         ; na = na_list(i)
               ishift = ishift_list(i) ; ikb = ishift + ih
               !
               IF (.NOT. is_multinp(i)) THEN
                  ps(ikb) =  CMPLX(deff(ih,ih,na) * becpr(ikb,ibnd_loc), KIND=DP)
               ELSE
                  nh_np = nh_list(i)
                  !
                  ps(ikb) = CMPLX( SUM( becpr(ishift+1:ishift+nh_np,ibnd_loc) &
                                    * deff(ih,1:nh_np,na) ), KIND=DP )
               ENDIF
             ENDDO
             !
             dot11 = 0._DP ; dot21 = 0._DP ; dot31 = 0._DP
             dot22 = 0._DP ; dot32 = 0._DP ; dot33 = 0._DP
             !
             !$acc parallel loop collapse(2) reduction(+:dot11,dot21,dot31,&
             !$acc&                                      dot22,dot32,dot33)
             DO na =1, nat
                DO i = 1, npw
                   np = ityp(na)
                   ijkb0 = shift(na)
                   nh_np = nh(np)
                   worksum = (0._DP,0._DP)
                   DO ih = 1, nh_np
                      ikb = ijkb0 + ih
                      worksum = worksum + ps(ikb) * dvkb(i,ikb,4)
                   ENDDO
                   Re_worksum = DBLE(worksum) ;  Im_worksum = DIMAG(worksum)
                   evci = evcv(i)
                   gk1  = gk(i,1) ;  gk2 = gk(i,2) ;  gk3 = gk(i,3)
                   qm1i = qm1(i)
                   !
                   cv = evci * CMPLX(gk1 * gk1  * qm1i, KIND=DP)
                   dot11 = dot11 + Re_worksum*DBLE(cv) + Im_worksum*DIMAG(cv)
                   !
                   cv = evci * CMPLX(gk2 * gk1 * qm1i, KIND=DP)
                   dot21 = dot21 + Re_worksum*DBLE(cv) + Im_worksum*DIMAG(cv)
                   !
                   cv = evci * CMPLX(gk3 * gk1 * qm1i, KIND=DP)
                   dot31 = dot31 + Re_worksum*DBLE(cv) + Im_worksum*DIMAG(cv)
                   !
                   cv = evci * CMPLX(gk2 * gk2 * qm1i, KIND=DP)
                   dot22 = dot22 + Re_worksum*DBLE(cv) + Im_worksum*DIMAG(cv)
                   !
                   cv = evci * CMPLX(gk3 * gk2 * qm1i, KIND=DP)
                   dot32 = dot32 + Re_worksum*DBLE(cv) + Im_worksum*DIMAG(cv)
                   !
                   cv = evci * CMPLX(gk3 * gk3 * qm1i, KIND=DP)
                   dot33 = dot33 + Re_worksum*DBLE(cv) + Im_worksum*DIMAG(cv)
                ENDDO
             ENDDO
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
             dot11 = 0._DP ; dot21 = 0._DP ; dot31 = 0._DP
             dot22 = 0._DP ; dot32 = 0._DP ; dot33 = 0._DP
             !
             !$acc parallel loop collapse(2) reduction(+:dot11,dot21,dot31,dot22,dot32,dot33)
             DO ikb = 1, nkb
                DO i = 1, npw
                   wsum1 = ps(ikb)*dvkb(i,ikb,1)
                   wsum2 = ps(ikb)*dvkb(i,ikb,2)
                   wsum3 = ps(ikb)*dvkb(i,ikb,3)
                   !
                   evci = evcv(i)
                   gk1 = gk(i,1)
                   gk2 = gk(i,2)
                   gk3 = gk(i,3)
                   !
                   cv = evci * CMPLX(gk1, KIND=DP)
                   dot11 = dot11 + DBLE(wsum1)* DBLE(cv) + DIMAG(wsum1)*DIMAG(cv)
                   dot21 = dot21 + DBLE(wsum2)* DBLE(cv) + DIMAG(wsum2)*DIMAG(cv)
                   dot31 = dot31 + DBLE(wsum3)* DBLE(cv) + DIMAG(wsum3)*DIMAG(cv)
                   !
                   cv = evci * CMPLX(gk2, KIND=DP)
                   dot22 = dot22 + DBLE(wsum2)* DBLE(cv) + DIMAG(wsum2)*DIMAG(cv) 
                   dot32 = dot32 + DBLE(wsum3)* DBLE(cv) + DIMAG(wsum3)*DIMAG(cv)
                   ! 
                   cv =  evci * CMPLX(gk3, KIND=DP)
                   dot33 = dot33 + DBLE(wsum3)* DBLE(cv) + DIMAG(wsum3)*DIMAG(cv)
                ENDDO
             ENDDO 
             !
             ! ... a factor 2 accounts for the other half of the G-vector sphere
             !
             sigmanlc(:,1) = sigmanlc(:,1) -4._DP * wg(ibnd,ik) * [dot11, dot21, dot31]
             sigmanlc(:,2) = sigmanlc(:,2) -4._DP * wg(ibnd,ik) * [0._DP, dot22, dot32]
             sigmanlc(:,3) = sigmanlc(:,3) -4._DP * wg(ibnd,ik) * [0._DP, 0._DP, dot33]
          ENDDO
          !
          IF ( nproc > 1 ) THEN
#if defined(__CUDA) || defined(_OPENACC)
             CALL errore( 'stres_us_gamma', &
                          'unexpected error nproc be 1 with GPU acceleration', 100 )
#else
             CALL mp_circular_shift_left( becp%r, icyc, becp%comm )
             becpr = becp%r
             CALL mp_circular_shift_left( becp%ibnd_begin, icyc, becp%comm )
             CALL mp_circular_shift_left( nbnd_loc, icyc, becp%comm )
#endif
          ENDIF
          !
       ENDDO
       !
10     CONTINUE
       !
       DO l = 1, 3
          sigmanlc(l,l) = sigmanlc(l,l) - evps
       ENDDO
       !
       !$acc end data
       !$acc end data
       DEALLOCATE( deff, ps )
       DEALLOCATE( dvkb )
       DEALLOCATE( becpr )
       !
       RETURN
       !
     END SUBROUTINE stres_us_gamma
     !
     !
     !----------------------------------------------------------------------
     SUBROUTINE stres_us_k()
       !----------------------------------------------------------------------
       !! nonlocal contribution to the stress - k-points version.
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER  :: na, np, ibnd, ipol, jpol, l, i, nt, &
                   ikb, jkb, ih, jh, is, js, ijs, ishift, nh_np
       REAL(DP) :: fac, evps, dot11, dot21, dot31, dot22, dot32, dot33, aux, &
                   Re_worksum, Re_worksum1, Re_worksum2, Im_worksum,         &
                   Im_worksum1, Im_worksum2
       COMPLEX(DP) :: qm1i, gk1, gk2, gk3, pss
       COMPLEX(DP) :: cv, cv1, cv2, worksum, worksum1, worksum2, evci, evc1i, &
                      evc2i, ps1, ps2, ps1d1, ps1d2, ps1d3, ps2d1, ps2d2,     &
                      ps2d3, psd1, psd2, psd3
       !
       REAL(DP), ALLOCATABLE :: deff(:,:,:)
       COMPLEX(DP), ALLOCATABLE :: deff_nc(:,:,:,:)
       COMPLEX(DP), ALLOCATABLE :: ps(:), ps_nc(:,:), dvkb(:,:,:)
       !
       REAL(DP) :: xyz(3,3)
       ! xyz are the three unit vectors in the x,y,z directions
       DATA xyz / 1._DP, 0._DP, 0._DP, 0._DP, 1._DP, 0._DP, 0._DP, &
                  0._DP, 1._DP /
       !
       COMPLEX(DP), ALLOCATABLE :: becpnc(:,:,:), becpk(:,:)
       !$acc declare device_resident(becpnc, becpk)
       !
       evps = 0._DP
       ! ... diagonal contribution
       !
       ALLOCATE( dvkb(npwx,nkb,4) )
       !$acc data create( dvkb )
       !
       CALL gen_us_dj( ik, dvkb(:,:,4) )
       IF ( lmaxkb > 0 ) THEN
         DO ipol = 1, 3
           CALL gen_us_dy( ik, xyz(1,ipol), dvkb(:,:,ipol) )
         ENDDO
       ENDIF
       !
       IF (noncolin) THEN
          ALLOCATE( ps_nc(nkb,npol) )
          ALLOCATE( deff_nc(nhm,nhm,nat,nspin) )
          ALLOCATE( becpnc(nkb,npol,nbnd) )
          !$acc kernels present(becp%nc)
          becpnc = becp%nc
          !$acc end kernels
       ELSE
          ALLOCATE( ps(nkb) )
          ALLOCATE( deff(nhm,nhm,nat) )
          ALLOCATE( becpk(nkb,nbnd) )
          !$acc kernels present(becp%k)
          becpk = becp%k
          !$acc end kernels
       ENDIF
       !$acc data create( ps, ps_nc, deff, deff_nc )
       !
       ! ... the contribution is calculated only on one processor because
       ! ... partial results are later summed over all processors
       !
       IF ( me_bgrp /= root_bgrp ) GO TO 100
       !
       DO ibnd = 1, nbnd
          fac = wg(ibnd,ik)
          IF (ABS(fac) < 1.d-9) CYCLE
          IF (noncolin) THEN
             CALL compute_deff_nc( deff_nc, et(ibnd,ik) )
          ELSE
             CALL compute_deff( deff, et(ibnd,ik) )
          ENDIF
          !
          IF (noncolin) THEN
            !
#if defined(_OPENACC)
            !$acc parallel loop reduction(+:evps)
#else
            !$omp parallel do reduction(+:evps), private(ih,na,ishift,ikb,aux,&
            !$omp&            ijs,is,js,nh_np,jkb)
#endif
            DO i = 1, itot
              !
              ih = ih_list(i)         ; na = na_list(i)
              ishift = ishift_list(i) ; ikb = ishift + ih
              aux = 0.d0
              !
              IF (.NOT. is_multinp(i)) THEN
                 ijs = 0
                 !$acc loop seq collapse(2) reduction(+:aux)
                 DO is = 1, npol
                   DO js = 1, npol
                      ijs = ijs + 1
                      aux = aux + fac * DBLE(deff_nc(ih,ih,na,ijs) * &
                                             CONJG(becpnc(ikb,is,ibnd)) * &
                                             becpnc(ikb,js,ibnd))
                   ENDDO
                 ENDDO
              ELSE
                 nh_np = nh_list(i)
                 ijs = 0
                 !$acc loop seq collapse(2) reduction(+:aux)
                 DO is = 1, npol
                   DO js = 1, npol
                      ijs = ijs + 1
                      aux = aux + fac * DBLE(deff_nc(ih,ih,na,ijs) * &
                                             CONJG(becpnc(ikb,is,ibnd)) * &
                                             becpnc(ikb,js,ibnd))
                   ENDDO
                 ENDDO
                 !$acc loop seq
                 DO jh = ih+1, nh_np
                    jkb = ishift + jh
                    ijs = 0
                    !$acc loop seq collapse(2) reduction(+:aux)
                    DO is = 1, npol
                      DO js = 1, npol
                         ijs = ijs + 1
                         aux = aux + 2._DP*fac * &
                                       DBLE(deff_nc(ih,jh,na,ijs) * &
                                            (CONJG(becpnc(ikb,is,ibnd)) * &
                                             becpnc(jkb,js,ibnd)))
                      ENDDO
                    ENDDO
                 ENDDO
              ENDIF
              !
              evps = evps + aux
              !
            ENDDO
#if !defined(_OPENACC)
            !$omp end parallel do
#endif
            !
          ELSE
            !
            aux = 0.d0
#if defined(_OPENACC)
            !$acc parallel loop reduction(+:evps)
#else
            !$omp parallel do reduction(+:evps) private(ih,na,ishift,ikb,&
            !$omp&            aux,nh_np)
#endif
            DO i = 1, itot
              !
              ih = ih_list(i)         ; na = na_list(i)
              ishift = ishift_list(i) ; ikb = ishift + ih
              !
              IF (.NOT. is_multinp(i)) THEN
                 aux = fac * deff(ih,ih,na) * &
                               ABS(becpk(ikb,ibnd) )**2
              ELSE
                 nh_np = nh_list(i)
                 aux = fac * deff(ih,ih,na) * ABS(becpk(ikb,ibnd) )**2 + &
                             SUM( deff(ih,ih+1:nh_np,na) * &
                                  fac * 2._DP*DBLE( CONJG(becpk(ikb,ibnd)) &
                                  * becpk(ishift+ih+1:ishift+nh_np,ibnd) ) )
              ENDIF
              !
              evps = evps + aux
              !
            ENDDO
#if !defined(_OPENACC)
            !$omp end parallel do
#endif
            !
          ENDIF
          !
       ENDDO
       !
       DO l = 1, 3
          sigmanlc(l,l) = sigmanlc(l,l) - evps
       ENDDO
       !
100    CONTINUE
       !
       ! ... non diagonal contribution - derivative of the bessel function
       !
       DO ibnd = 1, nbnd
         !
         !$acc kernels
         evcv(:) = evc(:,ibnd)
         !$acc end kernels
         !
         IF ( noncolin ) THEN
            !
            CALL compute_deff_nc( deff_nc, et(ibnd,ik) )
            !
#if defined(_OPENACC)
            !$acc parallel loop
#else
            !$omp parallel do private(ih,na,ishift,ikb,is,ijs,nh_np)
#endif
            DO i = 1, itot
              !
              ih = ih_list(i)         ; na = na_list(i)
              ishift = ishift_list(i) ; ikb = ishift + ih
              !
              IF (.NOT. is_multinp(i)) THEN
                 !
                 !$acc loop seq
                 DO is = 1, npol
                   ijs = (is-1)*npol
                   ps_nc(ikb,is) = SUM( becpnc(ikb,1:npol,ibnd) * &
                                        deff_nc(ih,ih,na,ijs+1:ijs+npol) )
                 ENDDO
                 !
              ELSE
                 !
                 nh_np = nh_list(i)
                 !
                 !$acc loop seq
                 DO is = 1, npol
                   ijs = (is-1)*npol
                   ps_nc(ikb,is) = SUM( becpnc(ishift+1:ishift+nh_np,1:npol,ibnd) * &
                                        deff_nc(ih,1:nh_np,na,ijs+1:ijs+npol) )
                 ENDDO
                 !
              ENDIF
            ENDDO
#if !defined(_OPENACC)
            !$omp end parallel do
#endif
            !
         ELSE
            !
            CALL compute_deff( deff, et(ibnd,ik) )
            !
#if defined(_OPENACC)
            !$acc parallel loop
#else
            !$omp parallel do private(ih,na,ishift,ikb,nh_np)
#endif
            DO i = 1, itot
               !
               ih = ih_list(i)         ; na = na_list(i)
               ishift = ishift_list(i) ; ikb = ishift + ih
               !
               IF (.NOT. is_multinp(i)) THEN
                  ps(ikb) = CMPLX(deeq(ih,ih,na,current_spin), KIND=DP) * &
                                               becpk(ikb,ibnd)
               ELSE 
                  nh_np = nh_list(i)
                  !
                  ps(ikb) = SUM( becpk(ishift+1:ishift+nh_np,ibnd) * &
                                               deff(ih,1:nh_np,na) )
               ENDIF
               !
            ENDDO
#if !defined(_OPENACC)
            !$omp end parallel do
#endif
            !
         ENDIF
         !
         dot11=0._DP ; dot21=0._DP ; dot31=0._DP
         dot22=0._DP ; dot32=0._DP ; dot33=0._DP
         !
         IF (noncolin) THEN
#if defined(_OPENACC)
            !$acc parallel loop collapse(2) reduction(+:dot11,dot21,dot31,&
            !$acc&                                      dot22,dot32,dot33)
#else
            !$omp parallel do collapse(2) reduction(+:dot11,dot21,dot31,dot22,&
            !$omp&    dot32,dot33) shared(evcv,qm1,gk,ps_nc,dvkb)
#endif
            DO ikb = 1, nkb
               DO i = 1, npw
                  evc1i = evcv(i)
                  evc2i = evcv(i+npwx)
                  qm1i = CMPLX(qm1(i), KIND=DP)
                  gk1 = CMPLX(gk(i,1), KIND=DP)
                  gk2 = CMPLX(gk(i,2), KIND=DP)
                  gk3 = CMPLX(gk(i,3), KIND=DP)
                  worksum1 = ps_nc(ikb,1) * dvkb(i,ikb,4)
                  worksum2 = ps_nc(ikb,2) * dvkb(i,ikb,4)
                  Re_worksum1 = DBLE(worksum1) ;  Im_worksum1 = DIMAG(worksum1)
                  Re_worksum2 = DBLE(worksum2) ;  Im_worksum2 = DIMAG(worksum2)
                  !
                  cv1 = evc1i * gk1 * gk1 * qm1i
                  cv2 = evc2i * gk1 * gk1 * qm1i
                  dot11 = dot11 + Re_worksum1*DBLE(cv1) + Im_worksum1*DIMAG(cv1) + &
                                  Re_worksum2*DBLE(cv2) + Im_worksum2*DIMAG(cv2)
                  !
                  cv1 = evc1i * gk2 * gk1 * qm1i
                  cv2 = evc2i * gk2 * gk1 * qm1i
                  dot21 = dot21 + Re_worksum1*DBLE(cv1) + Im_worksum1*DIMAG(cv1) + &
                                  Re_worksum2*DBLE(cv2) + Im_worksum2*DIMAG(cv2)
                  !
                  cv1 = evc1i * gk3 * gk1 * qm1i
                  cv2 = evc2i * gk3 * gk1 * qm1i
                  dot31 = dot31 + Re_worksum1*DBLE(cv1) + Im_worksum1*DIMAG(cv1) + &
                                  Re_worksum2*DBLE(cv2) + Im_worksum2*DIMAG(cv2)
                  !
                  cv1 = evc1i * gk2 * gk2 * qm1i
                  cv2 = evc2i * gk2 * gk2 * qm1i
                  dot22 = dot22 + Re_worksum1*DBLE(cv1) + Im_worksum1*DIMAG(cv1) + &
                                  Re_worksum2*DBLE(cv2) + Im_worksum2*DIMAG(cv2)
                  !
                  cv1 = evc1i * gk3 * gk2 * qm1i
                  cv2 = evc2i * gk3 * gk2 * qm1i
                  dot32 = dot32 + Re_worksum1*DBLE(cv1) + Im_worksum1*DIMAG(cv1) + &
                                  Re_worksum2*DBLE(cv2) + Im_worksum2*DIMAG(cv2)
                  !
                  cv1  = evc1i * gk3 * gk3 * qm1i
                  cv2  = evc2i * gk3 * gk3 * qm1i
                  dot33 = dot33 + Re_worksum1*DBLE(cv1) + Im_worksum1*DIMAG(cv1) + &
                                  Re_worksum2*DBLE(cv2) + Im_worksum2*DIMAG(cv2)
               ENDDO
            ENDDO
#if !defined(_OPENACC)
            !$omp end parallel do
#endif
            !
         ELSE
            !
#if defined(_OPENACC)
            !$acc parallel loop collapse(2) reduction(+:dot11,dot21,dot31,&
            !$acc&                                      dot22,dot32,dot33)
#else
            !$omp parallel do collapse(2) reduction(+:dot11,dot21,dot31,dot22,&
            !$omp&    dot32,dot33) shared(evcv,qm1,gk,ps,dvkb)
#endif
            DO ikb = 1, nkb
               DO i = 1, npw
                  !
                  worksum = ps(ikb) *dvkb(i,ikb,4)
                  Re_worksum = DBLE(worksum) ;  Im_worksum = DIMAG(worksum)
                  !
                  evci = evcv(i)
                  qm1i = CMPLX(qm1(i), KIND=DP)
                  gk1 = CMPLX(gk(i,1), KIND=DP)
                  gk2 = CMPLX(gk(i,2), KIND=DP)
                  gk3 = CMPLX(gk(i,3), KIND=DP)
                  !
                  cv = evci * gk1 * gk1 * qm1i
                  dot11 = dot11 + Re_worksum*DBLE(cv) + Im_worksum*DIMAG(cv)
                  !
                  cv = evci * gk2 * gk1 * qm1i
                  dot21 = dot21 + Re_worksum*DBLE(cv) + Im_worksum*DIMAG(cv)
                  !
                  cv = evci * gk3 * gk1 * qm1i
                  dot31 = dot31 + Re_worksum*DBLE(cv) + Im_worksum*DIMAG(cv)
                  !
                  cv = evci * gk2 * gk2 * qm1i
                  dot22 = dot22 + Re_worksum*DBLE(cv) + Im_worksum*DIMAG(cv)
                  !
                  cv = evci * gk3 * gk2 * qm1i
                  dot32 = dot32 + Re_worksum*DBLE(cv) + Im_worksum*DIMAG(cv)
                  !
                  cv = evci * gk3 * gk3 * qm1i
                  dot33 = dot33 + Re_worksum*DBLE(cv) + Im_worksum*DIMAG(cv)
                  !
               ENDDO
            ENDDO
#if !defined(_OPENACC)
            !$omp end parallel do
#endif
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
#if defined(_OPENACC)
            !$acc parallel loop collapse(2) reduction(+:dot11,dot21,dot31,&
            !$acc&                                      dot22,dot32,dot33)
#else
            !$omp parallel do collapse(2) reduction(+:dot11,dot21,dot31,dot22,&
            !$omp&    dot32,dot33) shared(evcv,gk,ps_nc,dvkb)
#endif
            DO ikb =1, nkb
               DO i = 1, npw
                  !
                  gk1 = CMPLX(gk(i,1), KIND=DP)
                  gk2 = CMPLX(gk(i,2), KIND=DP)
                  gk3 = CMPLX(gk(i,3), KIND=DP)
                  !
                  ps1 = ps_nc(ikb,1)
                  ps2 = ps_nc(ikb,2)
                  !
                  ps1d1 = ps1 * dvkb(i,ikb,1)
                  ps1d2 = ps1 * dvkb(i,ikb,2)
                  ps1d3 = ps1 * dvkb(i,ikb,3)
                  !
                  ps2d1 = ps2 * dvkb(i,ikb,1)
                  ps2d2 = ps2 * dvkb(i,ikb,2)
                  ps2d3 = ps2 * dvkb(i,ikb,3)
                  !
                  evc1i = evcv(i)
                  evc2i = evcv(i+npwx)
                  !
                  cv1 = evc1i * gk1
                  cv2 = evc2i * gk1
                  dot11 = dot11 + DBLE(ps1d1)*DBLE(cv1) + DIMAG(ps1d1)*DIMAG(cv1) + &
                                  DBLE(ps2d1)*DBLE(cv2) + DIMAG(ps2d1)*DIMAG(cv2)
                  !
                  dot21 = dot21 + DBLE(ps1d2)*DBLE(cv1) + DIMAG(ps1d2)*DIMAG(cv1) + &
                                  DBLE(ps2d2)*DBLE(cv2) + DIMAG(ps2d2)*DIMAG(cv2)
                  !
                  dot31 = dot31 + DBLE(ps1d3)*DBLE(cv1) + DIMAG(ps1d3)*DIMAG(cv1) + &
                                  DBLE(ps2d3)*DBLE(cv2) + DIMAG(ps2d3)*DIMAG(cv2)
                  !
                  cv1 = evc1i * gk2
                  cv2 = evc2i * gk2
                  dot22 = dot22 + DBLE(ps1d2)*DBLE(cv1) + DIMAG(ps1d2)*DIMAG(cv1) + &
                                  DBLE(ps2d2)*DBLE(cv2) + DIMAG(ps2d2)*DIMAG(cv2)
                  !
                  dot32 = dot32 + DBLE(ps1d3)*DBLE(cv1) + DIMAG(ps1d3)*DIMAG(cv1) + &
                                  DBLE(ps2d3)*DBLE(cv2) + DIMAG(ps2d3)*DIMAG(cv2)
                  !
                  cv1 = evc1i * gk3
                  cv2 = evc2i * gk3
                  dot33 = dot33 + DBLE(ps1d3)*DBLE(cv1) + DIMAG(ps1d3)*DIMAG(cv1) + &
                                  DBLE(ps2d3)*DBLE(cv2) + DIMAG(ps2d3)*DIMAG(cv2)
                  !
               ENDDO
            ENDDO
#if !defined(_OPENACC)
            !$omp end parallel do
#endif
            !
         ELSE
            !
#if defined(_OPENACC)
            !$acc parallel loop collapse(2) reduction(+:dot11,dot21,dot31,&
            !$acc&                                      dot22,dot32,dot33)
#else
            !$omp parallel do collapse(2) reduction(+:dot11,dot21,dot31,dot22,&
            !$omp&    dot32,dot33) shared(evcv,gk,ps,dvkb)
#endif
            DO ikb = 1, nkb
               DO i = 1, npw
                 pss  = ps(ikb)
                 psd1 = pss*dvkb(i,ikb,1)
                 psd2 = pss*dvkb(i,ikb,2)
                 psd3 = pss*dvkb(i,ikb,3)
                 evci = evcv(i)
                 gk1  = CMPLX(gk(i,1), KIND=DP)
                 gk2  = CMPLX(gk(i,2), KIND=DP)
                 gk3  = CMPLX(gk(i,3), KIND=DP)
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
#if !defined(_OPENACC)
            !$omp end parallel do
#endif
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
       !$acc end data
       IF (noncolin) THEN
          DEALLOCATE( ps_nc )
          DEALLOCATE( deff_nc )
       ELSE
          DEALLOCATE( ps )
          DEALLOCATE( deff )
       ENDIF
       !
       IF ( noncolin ) THEN
         DEALLOCATE( becpnc )
       ELSE
         DEALLOCATE( becpk )
       ENDIF
       !$acc end data
       DEALLOCATE( dvkb )
       !
       RETURN
       !
     END SUBROUTINE stres_us_k
     !
     !
END SUBROUTINE stres_us
