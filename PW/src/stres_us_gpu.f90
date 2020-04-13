!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stres_us_gpu( ik, gk, sigmanlc )
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
  USE wavefunctions_gpum, ONLY : using_evc
  USE wvfct_gpum,                ONLY : using_et
  USE uspp_gpum,                 ONLY : using_vkb, using_deeq
  USE becmod_subs_gpum,          ONLY : using_becp_auto
  !
  !^^^^^^^^^^^^gpu^^^^^^^^^
  !
  USE becmod,    ONLY :  allocate_bec_type,deallocate_bec_type
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
  REAL(DP), INTENT(IN)    :: gk(npwx,3)
  !! wave function components for fixed k-point
  REAL(DP), INTENT(INOUT) :: sigmanlc(3,3)
  !! stress tensor, non-local contribution
  !
  REAL(DP), ALLOCATABLE   :: qm1(:)
  REAL(DP)                :: q
  INTEGER                 :: npw, i
  !
  CALL using_evc(0)
  !
  !
  IF ( nkb == 0 ) RETURN
  !
  IF ( lsda ) current_spin = isk(ik)
  npw = ngk(ik)
  IF ( nks > 1 ) CALL using_vkb(1)
  IF ( nks > 1 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
  !
  CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm ) 
  
  CALL using_vkb(0); CALL using_becp_auto(2)
  CALL calbec( npw, vkb, evc, becp )
  !
  ALLOCATE( qm1( npwx ) )
  DO i = 1, npw
     q = SQRT( gk(i, 1)**2 + gk(i, 2)**2 + gk(i, 3)**2 )
     IF ( q > eps8 ) THEN
        qm1(i) = 1.D0 / q
     ELSE
        qm1(i) = 0.D0
     END IF
  END DO
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
  DEALLOCATE( qm1 )
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
       INTEGER                  :: na, np, ibnd, ipol, jpol, l, i,    &
                                   ikb, jkb, ih, jh, ijkb0, ibnd_loc, &
                                   nproc, nbnd_loc, nbnd_begin, icyc
       REAL(DP)                 :: fac, xyz(3,3), evps, ddot
       REAL(DP), ALLOCATABLE    :: deff(:,:,:)
       COMPLEX(DP), ALLOCATABLE :: work1(:), work2(:), dvkb(:,:)
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
       END IF

       ALLOCATE( work1( npwx ), work2( npwx ) ) 
       ALLOCATE( deff(nhm,nhm,nat) )
       !
       ! ... diagonal contribution - if the result from "calbec" are not 
       ! ... distributed, must be calculated on a single processor
       !
       evps = 0.D0
       IF ( nproc == 1 .AND. me_pool /= root_pool ) GO TO 100
       !
       
       CALL using_et(0) ! compute_deff : intent(in)
       DO ibnd_loc = 1, nbnd_loc
          ibnd = ibnd_loc + becp%ibnd_begin - 1 
          CALL compute_deff ( deff, et(ibnd,ik) )
          fac = wg(ibnd,ik)
          ijkb0 = 0
          DO np = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == np ) THEN
                   DO ih = 1, nh(np)
                      ikb = ijkb0 + ih
                      evps = evps + fac * deff(ih,ih,na) * &
                                    ABS( becp%r(ikb,ibnd_loc) )**2
                      !
                      IF ( upf(np)%tvanp .OR. upf(np)%is_multiproj ) THEN
                         !
                         ! ... only in the US case there is a contribution 
                         ! ... for jh<>ih
                         ! ... we use here the symmetry in the interchange of 
                         ! ... ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(np)
                            jkb = ijkb0 + jh
                            evps = evps + deff(ih,jh,na) * fac * 2.D0 * &
                                 becp%r(ikb,ibnd_loc) * becp%r(jkb,ibnd_loc)
                         END DO
                      END IF
                   END DO
                   ijkb0 = ijkb0 + nh(np)
                END IF
             END DO
          END DO
       END DO
       !
100    CONTINUE
       !
       ! ... non diagonal contribution - derivative of the bessel function
       !------------------------------------
       ALLOCATE( dvkb( npwx, nkb ) )
       !
       CALL gen_us_dj( ik, dvkb )
       !
       CALL using_et(0) ! compute_deff : intent(in)
       CALL using_evc(0)
       DO icyc = 0, nproc -1
          !
          DO ibnd_loc = 1, nbnd_loc
             !  
             ibnd = ibnd_loc + becp%ibnd_begin - 1 
             CALL compute_deff ( deff, et(ibnd,ik) )
             work2(:) = (0.D0,0.D0)
             ijkb0 = 0
             DO np = 1, ntyp
                DO na = 1, nat
                   IF ( ityp(na) == np ) THEN
                      DO ih = 1, nh(np)
                         ikb = ijkb0 + ih
                         IF ( .NOT. ( upf(np)%tvanp .OR. upf(np)%is_multiproj ) ) THEN
                            ps = becp%r(ikb,ibnd_loc) * deff(ih,ih,na)
                         ELSE
                            !
                            ! ... in the US case there is a contribution 
                            ! ... also for jh<>ih
                            !
                            ps = (0.D0,0.D0)
                            DO jh = 1, nh(np)
                               jkb = ijkb0 + jh
                               ps = ps + becp%r(jkb,ibnd_loc) * deff(ih,jh,na)
                            END DO
                         END IF
                         CALL zaxpy( npw, ps, dvkb(1,ikb), 1, work2, 1 )
                      END DO
                      ijkb0 = ijkb0 + nh(np)
                   END IF
                END DO
             END DO
             !
             ! ... a factor 2 accounts for the other half of the G-vector sphere
             !
             DO ipol = 1, 3
                DO jpol = 1, ipol
                   DO i = 1, npw
                      work1(i) = evc(i,ibnd) * gk(i, ipol) * gk(i, jpol) * qm1(i)
                   END DO
                   sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - &
                        4.D0 * wg(ibnd,ik) * &
                        ddot( 2 * npw, work1, 1, work2, 1 )
                END DO
             END DO
          END DO
          IF ( nproc > 1 ) THEN
             CALL mp_circular_shift_left(becp%r, icyc, becp%comm)
             CALL mp_circular_shift_left(becp%ibnd_begin, icyc, becp%comm)
             CALL mp_circular_shift_left(nbnd_loc, icyc, becp%comm)
          END IF
       END DO
       !
       ! ... non diagonal contribution - derivative of the spherical harmonics
       ! ... (no contribution from l=0)
       !
       IF ( lmaxkb == 0 ) GO TO 10
       !
       !------------------------------------
       CALL using_evc(0); CALL using_et(0) ! compute_deff : intent(in) (this is redundant)
       DO ipol = 1, 3
          !
          CALL gen_us_dy( ik, xyz(1,ipol), dvkb )
          !
          DO icyc = 0, nproc -1
             !
             DO ibnd_loc = 1, nbnd_loc
                ibnd = ibnd_loc + becp%ibnd_begin - 1 
                CALL compute_deff ( deff, et(ibnd,ik) )
                work2(:) = (0.D0,0.D0)
                ijkb0 = 0
                DO np = 1, ntyp
                   DO na = 1, nat
                      IF ( ityp(na) == np ) THEN
                         DO ih = 1, nh(np)
                            ikb = ijkb0 + ih
                            IF ( .NOT. ( upf(np)%tvanp .OR. upf(np)%is_multiproj ) ) THEN
                               ps = becp%r(ikb,ibnd_loc) * deff(ih,ih,na)
                            ELSE 
                               !
                               ! ... in the US case there is a contribution 
                               ! ... also for jh<>ih
                               !
                               ps = (0.D0,0.D0)
                               DO jh = 1, nh(np)
                                  jkb = ijkb0 + jh
                                  ps = ps + becp%r(jkb,ibnd_loc)*deff(ih,jh,na)
                               END DO
                            END IF
                            CALL zaxpy( npw, ps, dvkb(1,ikb), 1, work2, 1 )
                         END DO
                         ijkb0 = ijkb0 + nh(np)
                      END IF
                   END DO
                END DO
                !
                ! ... a factor 2 accounts for the other half of the G-vector sphere
                !
                DO jpol = 1, ipol
                   DO i = 1, npw
                      work1(i) = evc(i,ibnd) * gk(i, jpol)
                   END DO
                   sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - &
                        4.D0 * wg(ibnd,ik) * &
                        ddot( 2 * npw, work1, 1, work2, 1 )
                END DO
             END DO

             IF ( nproc > 1 ) THEN
                CALL mp_circular_shift_left(becp%r, icyc, becp%comm)
                CALL mp_circular_shift_left(becp%ibnd_begin, icyc, becp%comm)
                CALL mp_circular_shift_left(nbnd_loc, icyc, becp%comm)
             END IF

          ENDDO
       END DO

10     CONTINUE
       !
       DO l = 1, 3
          sigmanlc(l,l) = sigmanlc(l,l) - evps
       END DO
       !
       DEALLOCATE( dvkb )
       DEALLOCATE( deff, work2, work1 )
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
                  ikb, jkb, ih, jh, ijkb0, is, js, ijs
       REAL(DP) :: fac, xyz (3, 3), evps, ddot
       COMPLEX(DP), ALLOCATABLE :: work1(:), work2(:), dvkb(:,:,:)
       !
       
       
       !^^^^^^^^^^^^gpu
       REAL(DP) :: wg_nk, qm1i
       COMPLEX(DP) :: gk1, gk2, gk3
       !
       REAL(DP) :: dot11, dot21, dot31, dot22, dot32, dot33
       !
       COMPLEX(DP) :: cv, cv1, cv2, worksum, worksum1, worksum2, evci, evc1i, evc2i, & 
                      ps1, ps2, ps1d1, ps1d2, ps1d3, ps2d1, ps2d2, ps2d3, psd1, psd2, psd3
       
       INTEGER :: nt, npol2
       INTEGER :: ierr
       LOGICAL :: ismulti_np
       
       TYPE(bec_type_d), TARGET :: becp_d 
       
       INTEGER :: nh_np
       REAL(DP), POINTER :: qm1_d(:)
       REAL(DP), allocatable, device :: gk_d(:,:)
       REAL(DP), ALLOCATABLE, DEVICE :: deeq_d(:,:,:,:)
       COMPLEX(DP), ALLOCATABLE, DEVICE :: evc_d(:,:)
       COMPLEX(DP), ALLOCATABLE, DEVICE :: dvkb_d(:,:,:)
       COMPLEX(DP), ALLOCATABLE, DEVICE :: deff_d(:,:,:), deff_nc_d(:,:,:,:)
       
       COMPLEX(DP), POINTER :: ps_d(:), ps_nc_d(:,:)
       !COMPLEX(DP), POINTER ::  becpnc_d(:,:,:)
       COMPLEX(DP), ALLOCATABLE, DEVICE :: becpk_d(:,:)
       COMPLEX(DP), ALLOCATABLE, DEVICE :: becpnc_d(:,:,:)
       INTEGER, ALLOCATABLE :: shift(:)
       INTEGER, POINTER :: shift_d(:), ityp_d(:), nh_d(:)
       LOGICAL,ALLOCATABLE :: tvanp_d(:) 
       
       LOGICAL, ALLOCATABLE, DEVICE :: is_multinp_d(:)
       INTEGER, ALLOCATABLE, DEVICE :: na_d(:),ih_d(:),ikb_d(:),nhnp_d(:)
       integer :: bbb, itt
       complex(dp) :: becy, deqy
       !^^PROVV^^
       
#if defined(__CUDA) 
  ATTRIBUTES(DEVICE) :: shift_d, tvanp_d, ps_d, ps_nc_d,&
                        qm1_d, ityp_d, nh_d
#endif
       
       !^^^^^^^^^^^^^
       
       
       COMPLEX(DP), ALLOCATABLE :: work2_nc(:,:)
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
       if (noncolin) then
          ALLOCATE( work2_nc(npwx,npol) )
          ALLOCATE( deff_nc(nhm,nhm,nat,nspin) )
       endif
       !
       ALLOCATE( work1(npwx), work2(npwx) )
       !
       evps = 0.D0
       ! ... diagonal contribution
       !
       !
       !^^^^^^^^^^^^^^^
       ALLOCATE (shift(nat)) 
       CALL dev_buf%lock_buffer(shift_d, nat, ierr)
       ijkb0 = 0 
       DO nt = 1, ntyp
         DO na = 1, nat
           IF (ityp(na) == nt ) THEN
             shift (na) = ijkb0
             ijkb0 = ijkb0 + nh(nt)
           ENDIF
         ENDDO
       ENDDO
       shift_d(1:nat) = shift(1:nat)
       !^^^^^^^^^^^
       ALLOCATE(tvanp_d(ntyp))
       DO nt = 1, ntyp 
         tvanp_d(nt) = upf(nt)%tvanp .OR. upf(nt)%is_multiproj
       ENDDO
       
       ALLOCATE( gk_d(npwx,3) )
       gk_d = gk
       
       CALL dev_buf%lock_buffer( qm1_d, npwx, ierr )
       !$cuf kernel do(1) <<<*,*>>> 
       DO i = 1, npw
          q = SQRT( gk_d(i,1)**2 + gk_d(i,2)**2 + gk_d(i,3)**2 )
          IF ( q > eps8 ) THEN
            qm1_d(i) = 1._DP / q
          ELSE
            qm1_d(i) = 0._DP
          ENDIF
       ENDDO

       !^^^^^^^^^^^^^
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
       
       CALL dev_buf%lock_buffer( ityp_d, nat, ierr )
       ityp_d(1:nat) = ityp(1:nat)
       !
       IF (noncolin) THEN 
          CALL dev_buf%lock_buffer( ps_nc_d, [nkb, npol], ierr )
          !becpnc_d => becp_d%nc_d
          ALLOCATE( becpnc_d(nkb,npol,nbnd) )
          becpnc_d = becp%nc
          !
          ALLOCATE( deff_nc_d(nhm,nhm,nat,nspin) )
       ELSE 
          CALL dev_buf%lock_buffer ( ps_d, nkb, ierr)
          !becpk_d => becp_d%k_d
          ALLOCATE( becpk_d(nkb,nbnd) )
          becpk_d = becp%k      !correggi---devi usare becp_d.....
          
          ALLOCATE( deff_d(nhm,nhm,nat) )
       END IF
       !
       CALL dev_buf%lock_buffer( nh_d, ntyp, ierr )
       nh_d(1:ntyp) = nh(1:ntyp)
       !^^^^^^^^^^^^^^^^^^^^^^
       
       !
       CALL using_evc(0); CALL using_et(0) ! this is redundant
       CALL using_deeq(0)
       
       !
       !^^^^^^^^^^^^^^^^^^^^^^^^    1) metti gli using_...
       ALLOCATE( deeq_d(nhm,nhm,nat,nspin) )
       ALLOCATE( evc_d(npw,nbnd) )
       evc_d = evc
       deeq_d = deeq
       
       
       
       IF (noncolin) THEN
         !
       ELSE
         !
         bbb = nat*nhm
         ALLOCATE( is_multinp_d(bbb) )
         ALLOCATE( na_d(bbb), ih_d(bbb)    )
         ALLOCATE( ikb_d(bbb), nhnp_d(bbb) )
         !
         itt=0
         DO na = 1, nat
           np = ityp(na)
           ijkb0 = shift(na)
           nh_np = nh(np)
           ismulti_np = upf(np)%tvanp .OR. upf(np)%is_multiproj
           !$cuf kernel do (1) <<<*,*>>>
           DO ih = itt+1, itt+nh_np
             nhnp_d(ih) = nh_np
             na_d(ih)  = na
             ih_d(ih)  = ih-itt
             ikb_d(ih) = ijkb0 + ih - itt
             is_multinp_d(ih) = ismulti_np
           ENDDO 
           itt = itt + nh_np
         ENDDO
         !
       ENDIF 
       
       
       !
       CALL start_clock('knl_ciclo2')
       !
       CALL using_et(0) ! compute_deff : intent(in)
       
       
       IF ( me_bgrp /= root_bgrp ) GO TO 100
       !
       ! ... the contribution is calculated only on one processor because
       ! ... partial results are later summed over all processors
       
       
       DO ibnd = 1, nbnd
          fac = wg(ibnd,ik)
          IF (ABS(fac) < 1.d-9) CYCLE
          IF (noncolin) THEN
             CALL compute_deff_nc(deff_nc,et(ibnd,ik))
          ELSE
             CALL compute_deff_gpu( deff_d, et(ibnd,ik) )
          ENDIF
          ijkb0 = 0
          
          IF (noncolin) THEN
            !
            DO np = 1, ntyp
              DO na = 1, nat
                IF ( ityp(na) == np ) THEN
                   DO ih = 1, nh(np)
                      ikb = ijkb0 + ih
                      ! 
                      ijs=0
                      DO is=1,npol
                        DO js=1,npol
                           ijs=ijs+1
                           evps=evps+fac*deff_nc(ih,ih,na,ijs)*   &
                                     CONJG(becp%nc(ikb,is,ibnd))* &
                                           becp%nc(ikb,js,ibnd)
                        END DO
                      END DO
                      !
                      IF ( upf(np)%tvanp .OR. upf(np)%is_multiproj ) THEN
                         !
                         DO jh = ( ih + 1 ), nh(np)
                            jkb = ijkb0 + jh
                            ijs=0
                            DO is=1,npol
                              DO js=1,npol
                                 ijs=ijs+1
                                 evps = evps+2.d0*fac&
                                       *DBLE(deff_nc(ih,jh,na,ijs)*     &
                                       (CONJG( becp%nc(ikb,is,ibnd) ) * &
                                               becp%nc(jkb,js,ibnd))  )
                              END DO
                            END DO
                         END DO
                      END IF
                   END DO
                   ijkb0 = ijkb0 + nh(np)
                END IF
              END DO
            END DO
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
                 DBLE(SUM( deff_d(ih_d(i),ih_d(i)+1:nhnp_d(i),na_d(i)) * fac * 2.D0 * &
                               DBLE( CONJG( becpk_d(ikb_d(i),ibnd) ) * &
                               becpk_d(shift_d(na_d(i))+ih_d(i)+1:shift_d(na_d(i))+nhnp_d(i),ibnd) ) ))
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
            !$cuf kernel do (1) <<<*,*>>>
            DO na = 1, nat
               np = ityp_d(na)
               ijkb0 = shift_d(na)
               nh_np = nh_d(np)
               ismulti_np = tvanp_d(np)
               IF (.NOT. ismulti_np) THEN
                  !
                  DO ih = 1, nh_np
                     ikb = ijkb0 + ih
                     DO ijs = 1, npol2 
                        ps_nc_d(ikb,:) = (0.D0, 0.D0)
                        js = MOD(ijs-1, npol)
                        is = (ijs-1)/npol + 1
                        ps_nc_d(ikb,is) = ps_nc_d(ikb,is) + becpnc_d(ikb,js,ibnd)* &
                                          deff_nc_d(ih,ih,na,ijs)
                     ENDDO
                  ENDDO
                  !
               ELSE 
                  !
                  DO ih =1, nh_np
                     ikb = ijkb0 + ih
                     ps_nc_d(ikb,:) = (0._DP, 0._DP) 
                     DO jh = 1, nh_np
                        jkb = ijkb0 + jh
                        DO ijs = 1, npol2 
                           js = MOD(ijs-1,npol) 
                           is = (ijs-1)/npol + 1
                           ps_nc_d(ikb,is)= ps_nc_d(ikb,is) + becpnc_d(jkb,js,ibnd)* &
                                            deff_nc_d(ih,jh,na,ijs)
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
               IF (.NOT. is_multinp_d(i)) THEN !....questo caso fuori ciclo
                   !
                   deqy = CMPLX(deeq_d(ih_d(i),ih_d(i),na_d(i),current_spin))
                   becy = becpk_d(ikb_d(i),ibnd)
                   !
                   ps_d(ikb_d(i)) =  becy * deqy
                   !
               ELSE 
                   !
                   ps_d(ikb_d(i)) = SUM( becpk_d(shift_d(na_d(i))+1:shift_d(na_d(i))+nhnp_d(i),ibnd) * deff_d(ih_d(i),1:nhnp_d(i),na_d(i)) )
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
                  gk1 = CMPLX(gk_d(i,1))
                  gk2 = CMPLX(gk_d(i,2))
                  gk3 = CMPLX(gk_d(i,3))
                  worksum1 = ps_nc_d(ikb,1) * dvkb_d(i,ikb,4)     
                  worksum2 = ps_nc_d(ikb,2) * dvkb_d(i,ikb,4)   
                  !   
                  cv1 = evc1i * gk1 * gk1
                  cv2 = evc2i * gk1 * gk1
                  dot11 = dot11 + DREAL(worksum1)*DREAL(cv1) + DIMAG (worksum1)*DIMAG(cv1)
                  dot11 = dot11 + DREAL(worksum2)*DREAL(cv2) + DIMAG (worksum2)*DIMAG(cv2)
                  !   
                  cv1 = evc1i * gk2 * gk1
                  cv2 = evc2i * gk2 * gk1
                  dot21 = dot21 + DREAL(worksum1)*DREAL(cv1) + DIMAG (worksum1)*DIMAG(cv1)   
                  dot21 = dot21 + DREAL(worksum2)*DREAL(cv2) + DIMAG (worksum2)*DIMAG(cv2)   
                  !   
                  cv1 = evc1i * gk3 * gk1
                  cv2 = evc2i * gk3 * gk1
                  dot31 = dot31 + DREAL(worksum1)*DREAL(cv1) + DIMAG (worksum1)*DIMAG(cv1)   
                  dot31 = dot31 + DREAL(worksum2)*DREAL(cv2) + DIMAG (worksum2)*DIMAG(cv2)   
                  !   
                  cv1 = evc1i * gk2 * gk2
                  cv2 = evc2i * gk2 * gk2
                  dot22 = dot22 + DREAL(worksum1)*DREAL(cv1) + DIMAG (worksum1)*DIMAG(cv1)   
                  dot22 = dot22 + DREAL(worksum2)*DREAL(cv2) + DIMAG (worksum2)*DIMAG(cv2)   
                  !   
                  cv1 = evc1i * gk3 * gk2
                  cv2 = evc2i * gk3 * gk2
                  dot32 = dot32 + DREAL(worksum1)*DREAL(cv1) + DIMAG (worksum1)*DIMAG(cv1)
                  dot32 = dot32 + DREAL(worksum2)*DREAL(cv2) + DIMAG (worksum2)*DIMAG(cv2)
                  !   
                  cv1  = evc1i * gk3 * gk3
                  cv2  = evc2i * gk3 * gk3
                  dot33 = dot33 + DREAL(worksum1)*DREAL(cv1) + DIMAG (worksum1)*DIMAG(cv1)
                  dot33 = dot33 + DREAL(worksum2)*DREAL(cv2) + DIMAG (worksum2)*DIMAG(cv2)
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
                  qm1i = qm1_d(i)
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
               
            !stop  !------- 
               
            !IF ( me_bgrp /= root_bgrp ) stop     
            !   
         ENDIF
         !
         sigmanlc(:,1) = sigmanlc(:,1) - 2._DP * wg(ibnd,ik) * [dot11, dot21, dot31]
         sigmanlc(:,2) = sigmanlc(:,2) - 2._DP * wg(ibnd,ik) * [0._DP, dot22, dot32]
         sigmanlc(:,3) = sigmanlc(:,3) - 2._DP * wg(ibnd,ik) * [0._DP, 0._DP, dot33]
         !
         !do i = 1, 3
         !  print *, 'stres-USUSU:', sigmanlc(i,1), sigmanlc(i,2), sigmanlc(i,3)
         !enddo
         !!
         !stop
         !===========^=========mmcpy
         
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
                  dot11 = dot11 + DREAL(ps1d1)*DREAL(cv1) + DIMAG(ps1d1)*DIMAG(cv1)      
                  dot21 = dot21 + DREAL(ps1d2)*DREAL(cv1) + DIMAG(ps1d2)*DIMAG(cv1)      
                  dot31 = dot31 + DREAL(ps1d3)*DREAL(cv1) + DIMAG(ps1d3)*DIMAG(cv1)
                  !
                  dot11 = dot11 + DREAL(ps2d1)*DREAL(cv2) + DIMAG(ps2d1)*DIMAG(cv2)
                  dot21 = dot21 + DREAL(ps2d2)*DREAL(cv2) + DIMAG(ps2d2)*DIMAG(cv2)
                  dot31 = dot31 + DREAL(ps2d3)*DREAL(cv2) + DIMAG(ps2d3)*DIMAG(cv2)       
                  !      
                  cv1 = evc1i * gk2
                  cv2 = evc2i * gk2
                  dot22 = dot22 + DREAL(ps1d2)*DREAL(cv1) + DIMAG(ps1d2)*DIMAG(cv1)
                  dot32 = dot32 + DREAL(ps1d3)*DREAL(cv1) + DIMAG(ps1d3)*DIMAG(cv1)
                  !       
                  dot22 = dot22 + DREAL(ps2d2)*DREAL(cv2) + DIMAG(ps2d2)*DIMAG(cv2)
                  dot32 = dot32 + DREAL(ps2d3)*DREAL(cv2) + DIMAG(ps2d3)*DIMAG(cv2)
                  !       
                  cv1 = evc1i * gk3
                  cv2 = evc2i * gk3
                  dot33 = dot33 + DREAL(ps1d3)*DREAL(cv1) + DIMAG(ps1d3)*DIMAG(cv1)      
                  dot33 = dot33 + DREAL(ps2d3)*DREAL(cv2) + DIMAG(ps2d3)*DIMAG(cv2)      
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
       
       
         
  do i = 1, 3
    print *, 'stres-USUSU:', sigmanlc(i,1), sigmanlc(i,2), sigmanlc(i,3)
  enddo
  stop
       
       
       !stop !------------------------------
       
       
       
       IF (noncolin) THEN
           CALL dev_buf%release_buffer( ps_nc_d, ierr ) 
           DEALLOCATE( deff_nc_d )
       ELSE
           CALL dev_buf%release_buffer( ps_d, ierr) 
           !DEALLOCATE( deff_d )
       ENDIF
       !
       IF (noncolin) THEN
           DEALLOCATE( work2_nc )
           DEALLOCATE( deff_nc )
       ELSE
           DEALLOCATE( work2 )
       ENDIF
       DEALLOCATE( dvkb )
       DEALLOCATE( work1 )
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
