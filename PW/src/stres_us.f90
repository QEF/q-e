!
! Copyright (C) 2001-2024 Quantum ESPRESSO group
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
  !! IMPORTANT NOTICE: sum of partial results over processors (inter_pool_comm
  !! and intra_bgrp_group) is done in the calling routine "stres_knl"
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE constants,            ONLY : eps8
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wvfct,                ONLY : npwx, nbnd, wg, et
  USE control_flags,        ONLY : gamma_only, offload_type
  USE uspp_param,           ONLY : upf, lmaxkb, nh, nhm
  USE uspp,                 ONLY : nkb, vkb, deeq, qq_at, ofsbeta, deeq_nc, qq_so
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : noncolin, npol, lspinorb
  USE mp_pools,             ONLY : me_pool, root_pool
  USE mp_bands,             ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp, root_bgrp
  USE becmod,               ONLY : allocate_bec_type_acc, deallocate_bec_type_acc, &
                                   bec_type, becp, calbec
  USE mp,                   ONLY : mp_sum
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
  REAL(DP) :: q
  INTEGER  :: npw, i
  !
  !
  IF ( nkb == 0 ) RETURN
  !
  IF ( lsda ) current_spin = isk(ik)
  npw = ngk(ik)
  !
  IF ( nks > 1 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .TRUE. )
  !
  CALL allocate_bec_type_acc( nkb, nbnd, becp )
  !$acc data present( vkb, evc, becp )
  CALL calbec( offload_type, npw, vkb, evc, becp )
  !$acc end data
  !
  ALLOCATE( qm1(npwx) )
  !$acc data create(qm1) present(gk)
  !
  !$acc parallel loop
  DO i = 1, npw
     q = SQRT( gk(i,1)**2 + gk(i,2)**2 + gk(i,3)**2 )
     IF ( q > eps8 ) THEN
        qm1(i) = 1._DP / q
     ELSE
        qm1(i) = 0._DP
     ENDIF
  ENDDO
  !
  !$acc data present(et) copyin( ityp, wg, nh, ofsbeta ) 
  IF ( gamma_only ) THEN
     !
     CALL stres_us_gamma()
     !
  ELSE IF ( noncolin ) THEN
     ! 
     CALL stres_us_nc()
     !
  ELSE
     !
     CALL stres_us_k()
    !
  ENDIF
  !$acc end data
  !
  !$acc end data
  !
  DEALLOCATE( qm1 )
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
                   jkb, ih, jh, ijkb0, na_s, na_e, mykey
       REAL(DP) :: sigmaij
       COMPLEX(DP), ALLOCATABLE :: dvkb(:,:)
       TYPE(bec_type) :: becd
       !
       !! xyz are the three unit vectors in the x,y,z directions
       REAL(DP) :: xyz(3,3) 
       DATA xyz / 1._DP, 0._DP, 0._DP, &
                  0._DP, 1._DP, 0._DP, &
                  0._DP, 0._DP, 1._DP /
       !
       sigmaij = 0._DP
       !
       ! ... Calls to calbec are parallelized over the bgrp group
       ! ... The rest of the calculation is parallelized by subdividing 
       ! ... the atoms over the bgrp group
       !
       CALL block_distribute( nat, me_bgrp, nproc_bgrp, na_s, na_e, mykey )
       !
       compute_diag: IF ( mykey == 0 ) THEN
          !$acc parallel loop collapse(2) present(deeq, qq_at, becp%r) &
          !$acc reduction(+:sigmaij)
          DO na = na_s, na_e
             DO ibnd = 1, nbnd
                np = ityp(na)
                ijkb0 = ofsbeta(na)
                !$acc loop seq collapse(2)
                DO ih = 1, nh(np)
                   DO jh = 1, nh(np)
                      !! NOTE: in deff(ih,jh)=deeq(ih,jh)-et(i)*qq_at(ih,jh)
                      !! a nondiagonal contribution (ih,jh) is present only
                      !! for US-PP or multiprojector PP:
                      !!   IF ( upf(np)%tvanp) .or. upf(np)%is_multiproj ) 
                      !! but it may not be worth to make two distinct cases
                      !! For norm-conserving PP, qq_at=0 but again, it may
                      !! not be worth to make two distinct cases
                      !! The same applies to the two similar loops below
                      !
                      ikb = ijkb0 + ih
                      jkb = ijkb0 + jh
                      sigmaij = sigmaij + ( deeq(ih,jh,na,current_spin) - &
                              et(ibnd,ik)*qq_at(ih,jh,na)) * wg(ibnd,ik) * &
                              becp%r(ikb,ibnd) * becp%r(jkb,ibnd)
                   END DO
                END DO
             END DO
          END DO
          !
          DO l = 1, 3
             sigmanlc(l,l) = sigmanlc(l,l) - sigmaij
          ENDDO
       END IF compute_diag
       !
       ! ... non diagonal contribution - derivative of the Bessel function
       !
       CALL allocate_bec_type_acc( nkb, nbnd, becd )
       !
       ALLOCATE( dvkb(npwx,nkb) )
       !$acc data create(dvkb)
       !$acc data present( dvkb, becd )
       CALL gen_us_dj( ik, dvkb )
       !
       DO ipol = 1,3
          DO jpol = ipol, 3
             sigmaij = 0.0_dp 
             !$acc parallel loop collapse(2) present(vkb, gk, qm1)
             DO ikb = 1, nkb
                DO i = 1, npw
                   ! ... vkb is used here as work space
                   vkb(i,ikb) = dvkb(i,ikb) * gk(i,ipol) * gk(i,jpol) * qm1(i)
                END DO
             END DO
             ! ... becd like becp with derivatives of beta functions in vkb
             CALL calbec( offload_type, npw, vkb, evc, becd )
             !
             compute_djl: IF ( mykey == 0 ) THEN
                !$acc parallel loop collapse(2) present(deeq, qq_at, becp%r, becd%r) &
                !$acc reduction(+:sigmaij)
                DO na = na_s, na_e
                   DO ibnd = 1, nbnd
                      np = ityp(na)
                      ijkb0 = ofsbeta(na)
                      !$acc loop seq collapse(2)
                      DO ih = 1, nh(np)
                         DO jh = 1, nh(np)
                            ! 
                            ikb = ijkb0 + ih
                            jkb = ijkb0 + jh
                            sigmaij = sigmaij + ( deeq(ih,jh,na,current_spin) -&
                              et(ibnd,ik)*qq_at(ih,jh,na) ) * &
                              wg(ibnd,ik) * becp%r(ikb,ibnd) * &
                              becd%r(jkb,ibnd)
                         END DO
                      END DO
                   END DO
                END DO
                sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - 2.0_dp*sigmaij
             END IF compute_djl
          END DO
       END DO
       !
       ! ... non diagonal contribution - derivative of the spherical harmonics
       IF ( lmaxkb <= 0 ) GO TO 10
       !
       DO ipol = 1,3
          CALL gen_us_dy( ik, xyz(1,ipol), dvkb )
          DO jpol = ipol, 3
             sigmaij = 0.0_dp 
             !$acc parallel loop collapse(2) present(vkb, gk, qm1)
             DO ikb = 1, nkb
                DO i = 1, npw
                   ! ... vkb is used here as work space
                   vkb(i,ikb) = dvkb(i,ikb) * gk(i,jpol)
                END DO
             END DO
             ! ... becd like becp with derivatives of beta functions in dvkb
             CALL calbec( offload_type, npw, vkb, evc, becd )
             !
             compute_dylm: IF ( mykey == 0 ) THEN
                !$acc parallel loop collapse(2) present(deeq, qq_at, becp%r, becd%r) &
                !$acc reduction(+:sigmaij)
                DO na = na_s, na_e
                   DO ibnd = 1, nbnd
                      np = ityp(na)
                      ijkb0 = ofsbeta(na)
                      !$acc loop seq collapse(2)
                      DO ih = 1, nh(np)
                         DO jh = 1, nh(np)
                            ikb = ijkb0 + ih
                            jkb = ijkb0 + jh
                            sigmaij = sigmaij + (deeq(ih,jh,na,current_spin) -&
                                 et(ibnd,ik)*qq_at(ih,jh,na) ) * &
                                 wg(ibnd,ik) * becp%r(ikb,ibnd) * &
                                 becd%r(jkb,ibnd)
                         END DO
                      END DO
                   END DO
                END DO
                sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - 2.0_dp*sigmaij
             END IF compute_dylm
          END DO
       END DO
10     CONTINUE
       !$acc end data
       !$acc end data
       DEALLOCATE( dvkb )
       CALL deallocate_bec_type_acc( becd ) 
       DO ipol = 1,3
          DO jpol = 1,ipol-1
             sigmanlc(ipol,jpol) = sigmanlc(jpol,ipol)
          END DO
       END DO
       !
     END SUBROUTINE stres_us_gamma
     !
     !-----------------------------------------------------------------------
     SUBROUTINE stres_us_k()
       !-----------------------------------------------------------------------
       !! nonlocal contribution to the stress - k-point version
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER  :: na, np, nt, ibnd, ipol, jpol, l, i, ikb,  &
                   jkb, ih, jh, ijkb0, na_s, na_e, mykey
       REAL(DP) :: sigmaij
       COMPLEX(DP), ALLOCATABLE :: dvkb(:,:)
       TYPE(bec_type) :: becd
       !
       ! xyz are the three unit vectors in the x,y,z directions
       REAL(DP) :: xyz(3,3) 
       DATA xyz / 1._DP, 0._DP, 0._DP, &
                  0._DP, 1._DP, 0._DP, &
                  0._DP, 0._DP, 1._DP /
       !
       sigmaij = 0._DP
       !
       ! ... Calls to calbec are parallelized over the bgrp group
       ! ... The rest of the calculation is parallelized by subdividing 
       ! ... the atoms over the bgrp group
       !
       CALL block_distribute( nat, me_bgrp, nproc_bgrp, na_s, na_e, mykey )
       !
       compute_diag: IF ( mykey == 0 ) THEN
          !$acc parallel loop collapse(2) present(deeq, qq_at, becp%k) &
          !$acc reduction(+:sigmaij)
          DO na = na_s, na_e
             DO ibnd = 1, nbnd
                np = ityp(na)
                ijkb0 = ofsbeta(na)
                !$acc loop seq collapse(2)
                DO ih = 1, nh(np)
                   DO jh = 1, nh(np)
                      ! 
                      ! see note in the analogous loop for gamma-only case
                      !
                      ikb = ijkb0 + ih
                      jkb = ijkb0 + jh
                      sigmaij = sigmaij + ( deeq(ih,jh,na,current_spin) - &
                              et(ibnd,ik)*qq_at(ih,jh,na)) * wg(ibnd,ik) * &
                              DBLE( CONJG(becp%k(ikb,ibnd)) * becp%k(jkb,ibnd))
                   END DO
                END DO
             END DO
          END DO
          !
          DO l = 1, 3
             sigmanlc(l,l) = sigmanlc(l,l) - sigmaij
          ENDDO
       END IF compute_diag
       !
       ! ... non diagonal contribution - derivative of the Bessel function
       !
       CALL allocate_bec_type_acc( nkb, nbnd, becd )
       !
       ALLOCATE( dvkb(npwx,nkb) )
       !$acc data create(dvkb)
       !$acc data present( dvkb, becd )
       CALL gen_us_dj( ik, dvkb )
       !
       DO ipol = 1,3
          DO jpol = ipol, 3
             sigmaij = 0.0_dp 
             !$acc parallel loop collapse(2) present(vkb, gk, qm1)
             DO ikb = 1, nkb
                DO i = 1, npw
                   ! ... vkb is used here as work space
                   vkb(i,ikb) = dvkb(i,ikb) * gk(i,ipol) * gk(i,jpol) * qm1(i)
                END DO
             END DO
             ! ... becd like becp with derivatives of beta functions in vkb
             CALL calbec( offload_type, npw, vkb, evc, becd )
             !
             compute_djl: IF ( mykey == 0 ) THEN
                !$acc parallel loop collapse(2) present(deeq, qq_at, becp%k, becd%k) &
                !$acc reduction(+:sigmaij)
                DO na = na_s, na_e
                   DO ibnd = 1, nbnd
                      np = ityp(na)
                      ijkb0 = ofsbeta(na)
                      !$acc loop seq collapse(2)
                      DO ih = 1, nh(np)
                         DO jh = 1, nh(np)
                            ! 
                            ikb = ijkb0 + ih
                            jkb = ijkb0 + jh
                            sigmaij = sigmaij + ( deeq(ih,jh,na,current_spin) -&
                              et(ibnd,ik)*qq_at(ih,jh,na) ) * &
                              wg(ibnd,ik) * DBLE ( CONJG(becp%k(ikb,ibnd)) * &
                              becd%k(jkb,ibnd) )
                         END DO
                      END DO
                   END DO
                END DO
                sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - 2.0_dp*sigmaij
             END IF compute_djl
          END DO
       END DO
       !
       ! ... non diagonal contribution - derivative of the spherical harmonics
       IF ( lmaxkb <= 0 ) GO TO 10
       !
       DO ipol = 1,3
          CALL gen_us_dy( ik, xyz(1,ipol), dvkb )
          DO jpol = ipol, 3
             sigmaij = 0.0_dp 
             !$acc parallel loop collapse(2) present(vkb, gk, qm1)
             DO ikb = 1, nkb
                DO i = 1, npw
                   ! ... vkb is used here as work space
                   vkb(i,ikb) = dvkb(i,ikb) * gk(i,jpol)
                END DO
             END DO
             ! ... becd like becp with derivatives of beta functions in dvkb
             CALL calbec( offload_type, npw, vkb, evc, becd )
             !
             compute_dylm: IF ( mykey == 0 ) THEN
                !$acc parallel loop collapse(2) present(deeq, qq_at, becp%k, becd%k) &
                !$acc reduction(+:sigmaij)
                DO na = na_s, na_e
                   DO ibnd = 1, nbnd
                      np = ityp(na)
                      ijkb0 = ofsbeta(na)
                      !$acc loop seq collapse(2)
                      DO ih = 1, nh(np)
                         DO jh = 1, nh(np)
                            ikb = ijkb0 + ih
                            jkb = ijkb0 + jh
                            sigmaij = sigmaij + (deeq(ih,jh,na,current_spin) -&
                                 et(ibnd,ik)*qq_at(ih,jh,na) ) * &
                                 wg(ibnd,ik) * DBLE ( CONJG(becp%k(ikb,ibnd)) *&
                                 becd%k(jkb,ibnd) )
                         END DO
                      END DO
                   END DO
                END DO
                sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - 2.0_dp*sigmaij
             END IF compute_dylm
          END DO
       END DO
10     CONTINUE
       !$acc end data
       !$acc end data
       DEALLOCATE( dvkb )
       CALL deallocate_bec_type_acc( becd ) 
       DO ipol = 1,3
          DO jpol = 1,ipol-1
             sigmanlc(ipol,jpol) = sigmanlc(jpol,ipol)
          END DO
       END DO
       !
     END SUBROUTINE stres_us_k
     !
     !-----------------------------------------------------------------------
     SUBROUTINE stres_us_nc()
       !-----------------------------------------------------------------------
       !! nonlocal contribution to the stress - noncollinear version
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER  :: na, np, nt, ibnd, ipol, jpol, l, i, ikb,  &
                   is, js, ijs, jkb, ih, jh, ijkb0, na_s, na_e, mykey
       REAL(DP) :: sigmaij
       REAL(DP) :: ps(4)
       ! NOTE: This variable ps, even if not strictly necessary, has been defined as a workaround
       !       for an issue in nvhpc-24.3 compiler (it fails to perform acc reductons seemingly when
       !       a >3 number of sequential nested loops are present in a main parallel loop)
       !
       COMPLEX(DP), ALLOCATABLE :: dvkb(:,:)
       COMPLEX(DP) :: deff_nc
       TYPE(bec_type) :: becd
       !
       ! xyz are the three unit vectors in the x,y,z directions
       REAL(DP) :: xyz(3,3) 
       DATA xyz / 1._DP, 0._DP, 0._DP, &
                  0._DP, 1._DP, 0._DP, &
                  0._DP, 0._DP, 1._DP /
       !
       sigmaij = 0._DP
       !
       ! ... Calls to calbec are parallelized over the bgrp group
       ! ... The rest of the calculation is parallelized by subdividing 
       ! ... the atoms over the bgrp group
       !
       CALL block_distribute( nat, me_bgrp, nproc_bgrp, na_s, na_e, mykey )
       !
       compute_diag: IF ( mykey == 0 ) THEN
          !$acc parallel loop collapse(2) present(deeq_nc, qq_at, qq_so, becp%nc) &
          !$acc reduction(+:sigmaij) private(ps) 
          DO na = na_s, na_e
             DO ibnd = 1, nbnd
                np = ityp(na)
                ijkb0 = ofsbeta(na)
                !$acc loop seq collapse(2)
                DO ih = 1, nh(np)
                   DO jh = 1, nh(np)
                      !
                      ps = 0._dp
                      !$acc loop seq collapse(2)
                      DO is = 1, npol
                         DO js = 1, npol
                            ijs = (is-1)*npol+js
                            ! 
                            ! see note in the analogous loop for gamma-only case
                            ! note that unlike deff, deff_nc is complex
                            !
                            ikb = ijkb0 + ih
                            jkb = ijkb0 + jh
                            deff_nc = deeq_nc(ih,jh,na,ijs)
                            IF (lspinorb) THEN
                               deff_nc = deff_nc - et(ibnd,ik)*qq_so(ih,jh,ijs,np)
                            ELSE IF (is == js) THEN
                               deff_nc = deff_nc - et(ibnd,ik)*qq_at(ih,jh,na)
                            END IF
                            ps(ijs) = wg(ibnd,ik) * DBLE(deff_nc * &
                                 CONJG(becp%nc(ikb,is,ibnd)) * &
                                 becp%nc(jkb,js,ibnd))
                         END DO
                      END DO
                      sigmaij = sigmaij + ps(1) + ps(2) + ps(3) + ps(4) 
                   END DO
                END DO
             END DO
          END DO
          !
          DO l = 1, 3
             sigmanlc(l,l) = sigmanlc(l,l) - sigmaij
          ENDDO
       END IF compute_diag
       !
       ! ... non diagonal contribution - derivative of the Bessel function
       !
       CALL allocate_bec_type_acc( nkb, nbnd, becd )
       !
       ALLOCATE( dvkb(npwx,nkb) )
       !$acc data create(dvkb)
       !$acc data present( dvkb, becd )
       CALL gen_us_dj( ik, dvkb )
       !
       DO ipol = 1,3
          DO jpol = ipol, 3
             sigmaij = 0.0_dp 
             !$acc parallel loop collapse(2) present(vkb, gk, qm1)
             DO ikb = 1, nkb
                DO i = 1, npw
                   ! ... vkb is used here as work space
                   vkb(i,ikb) = dvkb(i,ikb) * gk(i,ipol) * gk(i,jpol) * qm1(i)
                END DO
             END DO
             ! ... becd like becp with derivatives of beta functions in vkb
             CALL calbec( offload_type, npw, vkb, evc, becd )
             !
             compute_djl: IF ( mykey == 0 ) THEN
                !$acc parallel loop collapse(2) present(deeq, qq_at, becp%k, becd%k) &
                !$acc reduction(+:sigmaij) private(ps)  
                DO na = na_s, na_e
                   DO ibnd = 1, nbnd
                      np = ityp(na)
                      ijkb0 = ofsbeta(na)
                      !$acc loop seq collapse(2)
                      DO ih = 1, nh(np)
                         DO jh = 1, nh(np)
                            !
                            ps = 0._dp
                            !$acc loop seq collapse(2)
                            DO is = 1, npol
                               DO js = 1, npol
                                  ijs = (is-1)*npol+js
                                  ! 
                                  ikb = ijkb0 + ih
                                  jkb = ijkb0 + jh
                                  deff_nc = deeq_nc(ih,jh,na,ijs)
                                  IF (lspinorb) THEN
                                     deff_nc = deff_nc - et(ibnd,ik)*qq_so(ih,jh,ijs,np)
                                  ELSE IF (is == js) THEN
                                     deff_nc = deff_nc - et(ibnd,ik)*qq_at(ih,jh,na)
                                  END IF
                                  ps(ijs) = wg(ibnd,ik) * DBLE(deff_nc * &
                                       CONJG(becp%nc(ikb,is,ibnd)) * becd%nc(jkb,js,ibnd))
                               END DO
                            END DO
                            sigmaij = sigmaij + ps(1) + ps(2) + ps(3) + ps(4) 
                         END DO
                      END DO
                   END DO
                END DO
                sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - 2.0_dp*sigmaij
             END IF compute_djl
          END DO
       END DO
       !
       ! ... non diagonal contribution - derivative of the spherical harmonics
       IF ( lmaxkb <= 0 ) GO TO 10
       !
       DO ipol = 1,3
          CALL gen_us_dy( ik, xyz(1,ipol), dvkb )
          DO jpol = ipol, 3
             sigmaij = 0.0_dp 
             !$acc parallel loop collapse(2) present(vkb, gk, qm1)
             DO ikb = 1, nkb
                DO i = 1, npw
                   ! ... vkb is used here as work space
                   vkb(i,ikb) = dvkb(i,ikb) * gk(i,jpol)
                END DO
             END DO
             ! ... becd like becp with derivatives of beta functions in dvkb
             CALL calbec( offload_type, npw, vkb, evc, becd )
             !
             compute_dylm: IF ( mykey == 0 ) THEN
                !$acc parallel loop collapse(2) present(deeq, qq_at, becp%k, becd%k) &
                !$acc reduction(+:sigmaij) private(ps) 
                DO na = na_s, na_e
                   DO ibnd = 1, nbnd
                      np = ityp(na)
                      ijkb0 = ofsbeta(na)
                      !$acc loop seq collapse(2)
                      DO ih = 1, nh(np)
                         DO jh = 1, nh(np)
                            !
                            ps = 0._dp
                            !$acc loop seq collapse(2)
                            DO is = 1, npol
                               DO js = 1, npol
                                  ijs = (is-1)*npol+js
                                  ikb = ijkb0 + ih
                                  jkb = ijkb0 + jh
                                  deff_nc = deeq_nc(ih,jh,na,ijs)
                                  IF (lspinorb) THEN
                                     deff_nc = deff_nc - et(ibnd,ik)*qq_so(ih,jh,ijs,np)
                                  ELSE IF (is == js) THEN
                                     deff_nc = deff_nc - et(ibnd,ik)*qq_at(ih,jh,na)
                                  END IF
                                  ps(ijs) = wg(ibnd,ik) * DBLE(deff_nc * &
                                       CONJG(becp%nc(ikb,is,ibnd)) * becd%nc(jkb,js,ibnd))
                               END DO
                            END DO
                            sigmaij = sigmaij + ps(1) + ps(2) + ps(3) + ps(4)
                         END DO
                      END DO
                   END DO
                END DO
                sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - 2.0_dp*sigmaij
             END IF compute_dylm
          END DO
       END DO
10     CONTINUE
       !$acc end data
       !$acc end data
       DEALLOCATE( dvkb )
       CALL deallocate_bec_type_acc( becd ) 
       DO ipol = 1,3
          DO jpol = 1,ipol-1
             sigmanlc(ipol,jpol) = sigmanlc(jpol,ipol)
          END DO
       END DO
       !
     END SUBROUTINE stres_us_nc
     !
END SUBROUTINE stres_us
