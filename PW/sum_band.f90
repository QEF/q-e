!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE sum_band()
  !----------------------------------------------------------------------------
  !
  ! ... calculates the symmetrized charge density and sum of occupied
  ! ... eigenvalues.
  ! ... this version works also for metals (gaussian spreading technique)  
  !
  USE kinds,                ONLY : DP
  USE ener,                 ONLY : eband
  USE control_flags,        ONLY : diago_full_acc, gamma_only, tqr
  USE cell_base,            ONLY : at, bg, omega, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   ngm, g, nl, nlm
  USE gsmooth,              ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                                   nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk
  USE ldaU,                 ONLY : lda_plus_U
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho
  USE symme,                ONLY : nsym, s, ftau
  USE io_files,             ONLY : iunwfc, nwordwfc, iunigk
  USE buffers,              ONLY : get_buffer
  USE uspp,                 ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv, okvan
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions_module, ONLY : evc, psic, psic_nc
  USE noncollin_module,     ONLY : noncolin, npol
  USE spin_orb,             ONLY : lspinorb, domag, so, fcoef
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et, btype
  USE mp_global,            ONLY : intra_image_comm, me_image, &
                                   root_image, npool, my_pool_id, inter_pool_comm
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE funct,                ONLY : dft_is_meta
  USE paw_onecenter,        ONLY : PAW_symmetrize
  USE paw_variables,        ONLY : okpaw
  USE becmod,               ONLY : allocate_bec, deallocate_bec
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: ikb, jkb, ijkb0, ih, jh, ijh, na, np
    ! counters on beta functions, atoms, pseudopotentials  
  INTEGER :: ir, is, ig, ibnd, ik, j
    ! counter on 3D r points
    ! counter on spin polarizations
    ! counter on g vectors
    ! counter on bands
    ! counter on k points  
  real (DP), allocatable :: kplusg (:)
  !
  !
  CALL start_clock( 'sum_band' )
  !
  becsum(:,:,:) = 0.D0
  rho%of_r(:,:)      = 0.D0
  rho%of_g(:,:)      = 0.D0
  if ( dft_is_meta() ) then
     rho%kin_r(:,:)      = 0.D0
     rho%kin_g(:,:)      = 0.D0
     allocate (kplusg(npwx))
  end if
  eband         = 0.D0  

  !
  ! ... calculates weights of Kohn-Sham orbitals used in calculation of rho
  !
  CALL weights ( )
  !
  IF ( diago_full_acc ) THEN
     !
     ! ... for diagonalization purposes all the bands are considered occupied
     !
     btype(:,:) = 1
     !
  ELSE
     !
     ! ... for diagonalization purposes a band is considered empty when its 
     ! ... occupation is less than 1.0 %
     !
     btype(:,:) = 1
     !
     FORALL( ik = 1:nks, wk(ik) > 0.D0 ) 
        !
        WHERE( wg(:,ik) / wk(ik) < 0.01D0 ) btype(:,ik) = 0
        !
     END FORALL
     !
  END IF
  !
  ! ... Needed for LDA+U
  !
  IF ( lda_plus_u ) CALL new_ns(rho%ns)  
  !
  IF ( okvan ) CALL allocate_bec (nkb,nbnd)
  !     
  ! ... specific routines are called to sum for each k point the contribution
  ! ... of the wavefunctions to the charge
  !
  IF ( gamma_only ) THEN
     !
     CALL sum_band_gamma()
     !
  ELSE
     !
     CALL sum_band_k()
     !
  END IF    
  !
  IF ( okvan ) CALL deallocate_bec ( )
  !
  ! ... If a double grid is used, interpolate onto the fine grid
  !
  IF ( doublegrid ) THEN
     !
     DO is = 1, nspin
        !
        CALL interpolate( rho%of_r(1,is), rho%of_r(1,is), 1 )
        if (dft_is_meta()) CALL interpolate(rho%kin_r(1,is),rho%kin_r(1,is),1)
        !
     END DO
     !
  END IF
  !
  ! ... Here we add the Ultrasoft contribution to the charge
  !
  IF ( okvan ) CALL addusdens()
  !
  IF ( noncolin .AND. .NOT. domag ) rho%of_r(:,2:4)=0.D0
  !
  CALL mp_sum( eband, inter_pool_comm )
  !
  ! ... symmetrization of the charge density (and local magnetization)
  !
#if defined (__PARA)
  !
  ! ... reduce charge density across pools
  !
  CALL mp_sum( rho%of_r, inter_pool_comm )
  if (dft_is_meta() ) CALL mp_sum( rho%kin_r, inter_pool_comm )
  !
  IF ( noncolin ) THEN
     !
     CALL psymrho( rho%of_r(1,1), nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau )
     !
     IF ( domag ) &
        CALL psymrho_mag( rho%of_r(1,2), nrx1, nrx2, nrx3, &
                          nr1, nr2, nr3, nsym, s, ftau, bg, at )

     !
  ELSE
     !
     DO is = 1, nspin
        !
        CALL psymrho( rho%of_r(1,is), nrx1, nrx2, nrx3, &
                      nr1, nr2, nr3, nsym, s, ftau )
        if (dft_is_meta() ) CALL psymrho( rho%kin_r(1,is), nrx1, nrx2, nrx3, &
                                          nr1, nr2, nr3, nsym, s, ftau )
        !
     END DO
     !
  END IF
  !
#else
  !
  IF ( noncolin ) THEN
     !
     CALL symrho( rho%of_r(1,1), nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau )
     !
     IF ( domag ) &
        CALL symrho_mag( rho%of_r(1,2), nrx1, nrx2, nrx3, &
                         nr1, nr2, nr3, nsym, s, ftau, bg, at )
     !
  ELSE
     !
     DO is = 1, nspin
        !
        CALL symrho( rho%of_r(1,is), nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau )
        if (dft_is_meta() ) CALL symrho( rho%kin_r(1,is), nrx1, nrx2, nrx3, &
                                         nr1, nr2, nr3, nsym, s, ftau )
        !
     END DO
     !
  END IF
  !
#endif
  ! ... Needed for PAW: becsum has to be symmetrized so that they reflect a real integral
  ! in k-space, not only on the irreducible zone. For USPP there is no need to do this as
  ! becsums are only used to compute the density, which is symmetrized later.
  !
  IF ( okpaw ) CALL PAW_symmetrize(becsum)
  !
  if (dft_is_meta() ) deallocate (kplusg)
  !
  ! ... synchronize rho%of_g to the calculated rho%of_r
  !
  DO is = 1, nspin
     !
     ! use psic as work array
     psic(:) = rho%of_r(:,is)
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
     rho%of_g(:,is) = psic(nl(:))
     !
     IF ( okvan .AND. tqr ) THEN
        ! ... in case the augmentation charges are computed in real space
        ! ... we apply an FFT filter to the density in real space to
        ! ... remove features that are not compatible with the FFT grid.
        !
        psic(:) = ( 0.D0, 0.D0 )
        psic(nl(:)) = rho%of_g(:,is)
        IF ( gamma_only ) psic(nlm(:)) = CONJG( rho%of_g(:,is) )
        CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
        rho%of_r(:,is) = psic(:)
        !
     END IF
     !
  END DO
  !
  ! ... the same for rho%kin_r and rho%kin_g
  !
  IF ( dft_is_meta()) THEN
     DO is = 1, nspin
        !
        ! use psic as work array
        psic(:) = rho%kin_r(:,is)
        CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
        rho%kin_g(:,is) = psic(nl(:))
        !
        IF ( okvan .AND. tqr ) THEN
           ! ... in case the augmentation charges are computed in real space
           ! ... we apply an FFT filter to the density in real space to
           ! ... remove features that are not compatible with the FFT grid.
           !
           psic(:) = ( 0.D0, 0.D0 )
           psic(nl(:)) = rho%kin_g(:,is)
           IF ( gamma_only ) psic(nlm(:)) = CONJG( rho%kin_g(:,is) )
           CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
           rho%kin_r(:,is) = psic(:)
           !
        END IF
        !
     END DO
  END IF
  !
  CALL stop_clock( 'sum_band' )      
  !
  RETURN
  !
  CONTAINS
     !
     ! ... internal procedures
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_gamma()
       !-----------------------------------------------------------------------
       !
       ! ... gamma version
       !
       USE becmod,        ONLY : rbecp, calbec
       USE fft_parallel,  ONLY : tg_cft3s
       USE fft_base,      ONLY : dffts
       USE control_flags, ONLY : use_task_groups
       USE mp_global,     ONLY : nogrp, nolist, ogrp_comm, me_pool
       USE mp,            ONLY : mp_sum
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1, w2
         ! weights
       INTEGER  :: idx, ioff, incr, v_siz, j, ip
       COMPLEX(DP), ALLOCATABLE :: tg_psi(:)
       REAL(DP),    ALLOCATABLE :: tg_rho(:)
       LOGICAL  :: use_tg
       !
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       IF ( nks > 1 ) REWIND( iunigk )
       !
       use_tg = ( use_task_groups ) .AND. ( nbnd >= nogrp )
       !
       incr = 2
       !
       IF( use_tg ) THEN
          !
          IF( dft_is_meta() ) &
             CALL errore( ' sum_band ', ' task groups with meta dft, not yet implemented ', 1 )
          !
          v_siz = dffts%nnrx * nogrp
          !
          ALLOCATE( tg_psi( v_siz ) )
          ALLOCATE( tg_rho( v_siz ) )
          !
          tg_rho = 0.0_DP
          !
          incr  = 2 * nogrp
          !
       END IF
       !
       k_loop: DO ik = 1, nks
          !
          IF ( lsda ) current_spin = isk(ik)
          npw = ngk(ik)
          !
          IF ( nks > 1 ) THEN
             !
             READ( iunigk ) igk
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
             !
          END IF
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk, xk(1,ik), vkb )
          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
          DO ibnd = 1, nbnd
             !
             ! ... the sum of eband and demet is the integral for  
             ! ... e < ef of e n(e) which reduces for degauss=0 to the sum of 
             ! ... the eigenvalues.
             !
             eband = eband + et(ibnd,ik) * wg(ibnd,ik)
             !
          END DO
          !
          DO ibnd = 1, nbnd, incr
             !
             IF( use_tg ) THEN
                !
                tg_psi(:) = ( 0.D0, 0.D0 )
                ioff   = 0
                !
                DO idx = 1, 2*nogrp, 2
                   !
                   ! ... 2*nogrp ffts at the same time
                   !
                   IF( idx + ibnd - 1 < nbnd ) THEN
                      DO j = 1, npw
                         tg_psi(nls (igk(j))+ioff) =        evc(j,idx+ibnd-1) + (0.0d0,1.d0) * evc(j,idx+ibnd)
                         tg_psi(nlsm(igk(j))+ioff) = CONJG( evc(j,idx+ibnd-1) - (0.0d0,1.d0) * evc(j,idx+ibnd) )
                      END DO
                   ELSE IF( idx + ibnd - 1 == nbnd ) THEN
                      DO j = 1, npw
                         tg_psi(nls (igk(j))+ioff) =        evc(j,idx+ibnd-1)
                         tg_psi(nlsm(igk(j))+ioff) = CONJG( evc(j,idx+ibnd-1) )
                      END DO
                   END IF

                   ioff = ioff + dffts%nnrx

                END DO
                !
                CALL tg_cft3s ( tg_psi, dffts, 2, use_tg )
                !
                ! Now the first proc of the group holds the first two bands
                ! of the 2*nogrp bands that we are processing at the same time,
                ! the second proc. holds the third and fourth band
                ! and so on
                !
                ! Compute the proper factor for each band
                !
                DO idx = 1, nogrp
                   IF( nolist( idx ) == me_pool ) EXIT
                END DO
                !
                ! Remember two bands are packed in a single array :
                ! proc 0 has bands ibnd   and ibnd+1
                ! proc 1 has bands ibnd+2 and ibnd+3
                ! ....
                !
                idx = 2 * idx - 1
                !
                IF( idx + ibnd - 1 < nbnd ) THEN
                   w1 = wg( idx + ibnd - 1, ik) / omega
                   w2 = wg( idx + ibnd    , ik) / omega
                ELSE IF( idx + ibnd - 1 == nbnd ) THEN
                   w1 = wg( idx + ibnd - 1, ik) / omega
                   w2 = w1
                ELSE 
                   w1 = 0.0d0
                   w2 = w1
                END IF
                !
                DO ir = 1, dffts%tg_npp( me_pool + 1 ) * dffts%nr1x * dffts%nr2x
                   !
                   tg_rho(ir) = tg_rho(ir) + &
                                        w1 *  DBLE( tg_psi(ir) )**2 + &
                                        w2 * AIMAG( tg_psi(ir) )**2
                   !
                END DO
                !
             ELSE
                !
                psic(:) = ( 0.D0, 0.D0 )
                !
                IF ( ibnd < nbnd ) THEN
                   !
                   ! ... two ffts at the same time
                   !
                   psic(nls(igk(1:npw)))  = evc(1:npw,ibnd) + &
                                               ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1)
                   psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ibnd) - &
                                               ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                   !
                ELSE
                   !
                   psic(nls(igk(1:npw)))  = evc(1:npw,ibnd)
                   psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ibnd) )
                   !
                END IF
                !
                CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
                !
                w1 = wg(ibnd,ik) / omega
                !
                ! ... increment the charge density ...
                !
                IF ( ibnd < nbnd ) THEN
                   !
                   ! ... two ffts at the same time
                   !
                   w2 = wg(ibnd+1,ik) / omega
                   !
                ELSE
                   !
                   w2 = w1
                   !
                END IF
                !
                DO ir = 1, nrxxs
                   !
                   rho%of_r(ir,current_spin) = rho%of_r(ir,current_spin) + &
                                                 w1 *  DBLE( psic(ir) )**2 + &
                                                 w2 * AIMAG( psic(ir) )**2
                   !
                END DO
                !
             END IF
             !
             IF (dft_is_meta()) THEN
                DO j=1,3
                   psic(:) = ( 0.D0, 0.D0 )
                   !
                   kplusg (1:npw) = (xk(j,ik)+g(j,igk(1:npw))) * tpiba

                   IF ( ibnd < nbnd ) THEN
                      ! ... two ffts at the same time
                      psic(nls(igk(1:npw))) = CMPLX (0d0, kplusg(1:npw)) * &
                                            ( evc(1:npw,ibnd) + &
                                            ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                      psic(nlsm(igk(1:npw))) = CMPLX (0d0, -kplusg(1:npw)) * &
                                       CONJG( evc(1:npw,ibnd) - &
                                            ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                   ELSE
                      psic(nls(igk(1:npw))) = CMPLX (0d0, kplusg(1:npw)) * &
                                              evc(1:npw,ibnd) 
                      psic(nlsm(igk(1:npw))) = CMPLX (0d0, -kplusg(1:npw)) * &
                                       CONJG( evc(1:npw,ibnd) )
                   END IF
                   !
                   CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
                   !
                   ! ... increment the kinetic energy density ...
                   !
                   DO ir = 1, nrxxs
                      rho%kin_r(ir,current_spin) = &
                                           rho%kin_r(ir,current_spin) + &
                                           w1 *  DBLE( psic(ir) )**2 + &
                                           w2 * AIMAG( psic(ir) )**2
                   END DO
                   !
                END DO
             END IF
             !
             !
          END DO
          !
          IF( use_tg ) THEN
             !
             ! reduce the group charge
             !
             CALL mp_sum( tg_rho, gid = ogrp_comm )
             !
             ioff = 0
             DO idx = 1, nogrp
                IF( me_pool == nolist( idx ) ) EXIT
                ioff = ioff + dffts%nr1x * dffts%nr2x * dffts%npp( nolist( idx ) + 1 )
             END DO
             !
             ! copy the charge back to the processor location
             !
             DO ir = 1, nrxxs
                rho%of_r(ir,current_spin) = rho%of_r(ir,current_spin) + tg_rho(ir+ioff)
             END DO
     
          END IF
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( .NOT. okvan ) CYCLE k_loop
          !
          CALL calbec( npw, vkb, evc, rbecp )
          !
          CALL start_clock( 'sum_band:becsum' )
          !
          DO ibnd = 1, nbnd
             !
             w1 = wg(ibnd,ik)
             ijkb0 = 0
             !
             DO np = 1, ntyp
                !
                IF ( upf(np)%tvanp ) THEN
                   !
                   DO na = 1, nat
                      !
                      IF ( ityp(na) == np ) THEN
                         !
                         ijh = 1
                         !
                         DO ih = 1, nh(np)
                            !
                            ikb = ijkb0 + ih
                            !
                            becsum(ijh,na,current_spin) = &
                                            becsum(ijh,na,current_spin) + &
                                            w1 *rbecp(ikb,ibnd) *rbecp(ikb,ibnd)
                            !
                            ijh = ijh + 1
                            !
                            DO jh = ( ih + 1 ), nh(np)
                               !
                               jkb = ijkb0 + jh
                               !
                               becsum(ijh,na,current_spin) = &
                                     becsum(ijh,na,current_spin) + &
                                     w1 * 2.D0 *rbecp(ikb,ibnd) *rbecp(jkb,ibnd)
                               !
                               ijh = ijh + 1
                               !
                            END DO
                            !
                         END DO
                         !
                         ijkb0 = ijkb0 + nh(np)
                         !
                      END IF
                      !
                   END DO
                   !
                ELSE
                   !
                   DO na = 1, nat
                      !
                      IF ( ityp(na) == np ) ijkb0 = ijkb0 + nh(np)
                      !
                   END DO
                   !
                END IF
                !
             END DO
             !
          END DO
          !
          CALL stop_clock( 'sum_band:becsum' )
          !
       END DO k_loop
       !
       IF( use_tg ) THEN
          DEALLOCATE( tg_psi )
          DEALLOCATE( tg_rho )
       END IF
       !
       RETURN
       !
     END SUBROUTINE sum_band_gamma
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_k()
       !-----------------------------------------------------------------------
       !
       ! ... k-points version
       !
       USE becmod, ONLY : becp, becp_nc, calbec
       USE fft_parallel,  ONLY : tg_cft3s
       USE fft_base,      ONLY : dffts
       USE control_flags, ONLY : use_task_groups
       USE mp_global,     ONLY : nogrp, nolist, ogrp_comm, me_pool
       USE mp,            ONLY : mp_sum
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1
       ! weights
       COMPLEX(DP), ALLOCATABLE :: becsum_nc(:,:,:,:)
       !
       INTEGER :: ipol, kh, kkb, is1, is2, js
       !
       INTEGER  :: idx, ioff, incr, v_siz, j, ip
       COMPLEX(DP), ALLOCATABLE :: tg_psi(:)
       REAL(DP),    ALLOCATABLE :: tg_rho(:)
       LOGICAL  :: use_tg
       !

       IF (okvan .AND.  noncolin) THEN
          ALLOCATE(becsum_nc(nhm*(nhm+1)/2,nat,npol,npol))
          becsum_nc=(0.d0, 0.d0)
       ENDIF
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       IF ( nks > 1 ) REWIND( iunigk )
       !
       use_tg = ( use_task_groups ) .AND. ( nbnd >= nogrp )
       use_tg = use_tg .AND. ( .NOT. noncolin )
       use_tg = use_tg .AND. ( .NOT. dft_is_meta() )
       !
       incr = 1
       !
       IF( use_tg ) THEN
          !
          v_siz = dffts%nnrx * nogrp
          !
          ALLOCATE( tg_psi( v_siz ) )
          ALLOCATE( tg_rho( v_siz ) )
          !
          tg_rho = 0.0_DP
          !
          incr  = nogrp
          !
       END IF
       !
       k_loop: DO ik = 1, nks
          !
          IF ( lsda ) current_spin = isk(ik)
          npw = ngk (ik)
          !
          IF ( nks > 1 ) THEN
             !
             READ( iunigk ) igk
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
             !
          END IF
          !
          IF ( nkb > 0 ) &
             CALL init_us_2( npw, igk, xk(1,ik), vkb )
          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
          DO ibnd = 1, nbnd, incr
             !
             IF( use_tg ) THEN   
                DO idx = 1, nogrp
                   IF( idx + ibnd - 1 <= nbnd ) eband = eband + et( idx + ibnd - 1, ik ) * wg( idx + ibnd - 1, ik )
                END DO
             ELSE
                eband = eband + et( ibnd, ik ) * wg( ibnd, ik )
             END IF
             !
             ! ... the sum of eband and demet is the integral for e < ef of 
             ! ... e n(e) which reduces for degauss=0 to the sum of the 
             ! ... eigenvalues 
             w1 = wg(ibnd,ik) / omega
             !
             IF (noncolin) THEN
                psic_nc = (0.D0,0.D0)
                DO ig = 1, npw
                   psic_nc(nls(igk(ig)),1)=evc(ig     ,ibnd)
                   psic_nc(nls(igk(ig)),2)=evc(ig+npwx,ibnd)
                END DO
                call cft3s (psic_nc(1,1), nr1s, nr2s, nr3s, &
                                          nrx1s,nrx2s,nrx3s, 2)
                call cft3s (psic_nc(1,2), nr1s, nr2s, nr3s, &
                                          nrx1s,nrx2s,nrx3s, 2)
                !
                ! increment the charge density ...
                !
                DO ipol=1,npol
                   DO ir = 1, nrxxs
                      rho%of_r (ir, 1) = rho%of_r (ir, 1) + &
                      w1*( DBLE(psic_nc(ir,ipol))**2+AIMAG(psic_nc(ir,ipol))**2)
                   END DO
                END DO
                !
                ! In this case, calculate also the three
                ! components of the magnetization (stored in rho%of_r(ir,2-4))
                !
                IF (domag) THEN
                   DO ir = 1,nrxxs
                      rho%of_r(ir,2) = rho%of_r(ir,2) + w1*2.D0* &
                         (DBLE(psic_nc(ir,1))* DBLE(psic_nc(ir,2)) + &
                         AIMAG(psic_nc(ir,1))*AIMAG(psic_nc(ir,2)))

                      rho%of_r(ir,3) = rho%of_r(ir,3) + w1*2.D0* &
                         (DBLE(psic_nc(ir,1))*AIMAG(psic_nc(ir,2)) - &
                          DBLE(psic_nc(ir,2))*AIMAG(psic_nc(ir,1)))

                      rho%of_r(ir,4) = rho%of_r(ir,4) + w1* &
                         (DBLE(psic_nc(ir,1))**2+AIMAG(psic_nc(ir,1))**2 &
                         -DBLE(psic_nc(ir,2))**2-AIMAG(psic_nc(ir,2))**2)
                   END DO
                ELSE
                   rho%of_r(:,2:4)=0.0_DP
                END IF
                !
             ELSE
                !
                IF( use_tg ) THEN
                   !
                   tg_psi(:) = ( 0.D0, 0.D0 )
                   ioff   = 0
                   !
                   DO idx = 1, nogrp
                      !
                      ! ... nogrp ffts at the same time
                      !
                      IF( idx + ibnd - 1 <= nbnd ) THEN
                         DO j = 1, npw
                            tg_psi( nls( igk( j ) ) + ioff ) = evc( j, idx+ibnd-1 )
                         END DO
                      END IF

                      ioff = ioff + dffts%nnrx

                   END DO
                   !
                   CALL tg_cft3s ( tg_psi, dffts, 2, use_tg )
                   !
                   ! Now the first proc of the group holds the first band
                   ! of the nogrp bands that we are processing at the same time,
                   ! the second proc. holds the second and so on
                   !
                   ! Compute the proper factor for each band
                   !
                   DO idx = 1, nogrp
                      IF( nolist( idx ) == me_pool ) EXIT
                   END DO
                   !
                   ! Remember
                   ! proc 0 has bands ibnd
                   ! proc 1 has bands ibnd+1
                   ! ....
                   !
                   IF( idx + ibnd - 1 <= nbnd ) THEN
                      w1 = wg( idx + ibnd - 1, ik) / omega
                   ELSE 
                      w1 = 0.0d0
                   END IF
                   !
                   DO ir = 1, dffts%tg_npp( me_pool + 1 ) * dffts%nr1x * dffts%nr2x
                      !
                      tg_rho(ir) = tg_rho(ir) + w1 *  ( DBLE( tg_psi(ir) )**2 + AIMAG( tg_psi(ir) )**2 )
                      !
                   END DO
                   !
                ELSE
                   !
                   psic(:) = ( 0.D0, 0.D0 )
                   !
                   psic(nls(igk(1:npw))) = evc(1:npw,ibnd)
                   !
                   CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
                   !
                   ! ... increment the charge density ...
                   !
                   DO ir = 1, nrxxs
                      !
                      rho%of_r(ir,current_spin) = rho%of_r(ir,current_spin) + &
                                               w1 * ( DBLE( psic(ir) )**2 + &
                                                     AIMAG( psic(ir) )**2 )
                      !
                   END DO

                END IF
                !
                IF (dft_is_meta()) THEN
                   DO j=1,3
                      psic(:) = ( 0.D0, 0.D0 )
                      !
                      kplusg (1:npw) = (xk(j,ik)+g(j,igk(1:npw))) * tpiba
                      psic(nls(igk(1:npw))) = CMPLX (0d0, kplusg(1:npw)) * &
                                              evc(1:npw,ibnd)
                      !
                      CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
                      !
                      ! ... increment the kinetic energy density ...
                      !
                      DO ir = 1, nrxxs
                         rho%kin_r(ir,current_spin) = &
                                                rho%kin_r(ir,current_spin) + &
                                                w1 * ( DBLE( psic(ir) )**2 + &
                                                      AIMAG( psic(ir) )**2 )
                      END DO
                      !
                   END DO
                END IF
                !
             END IF
             !
          END DO
          !
          IF( use_tg ) THEN
             !
             ! reduce the group charge
             !
             CALL mp_sum( tg_rho, gid = ogrp_comm )
             !
             ioff = 0
             DO idx = 1, nogrp
                IF( me_pool == nolist( idx ) ) EXIT
                ioff = ioff + dffts%nr1x * dffts%nr2x * dffts%npp( nolist( idx ) + 1 )
             END DO
             !
             ! copy the charge back to the proper processor location
             !
             DO ir = 1, nrxxs
                rho%of_r(ir,current_spin) = rho%of_r(ir,current_spin) + tg_rho(ir+ioff)
             END DO
             !
          END IF 
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( .NOT. okvan ) CYCLE k_loop
          !
          IF (noncolin) THEN
             CALL calbec( npw, vkb, evc, becp_nc )
          ELSE
             CALL calbec( npw, vkb, evc, becp )
          ENDIF
          !
          CALL start_clock( 'sum_band:becsum' )
          !
          DO ibnd = 1, nbnd
             !
             w1 = wg(ibnd,ik)
             ijkb0 = 0
             !
             DO np = 1, ntyp
                !
                IF ( upf(np)%tvanp ) THEN
                   !
                   DO na = 1, nat
                      !
                      IF (ityp(na)==np) THEN
                         !
                         ijh = 1
                         !
                         DO ih = 1, nh(np)
                            !
                            ikb = ijkb0 + ih
                            !
                            IF (noncolin) THEN
                               !
                               DO is=1,npol
                                  !
                                  DO js=1,npol
                                     becsum_nc(ijh,na,is,js) =         &
                                         becsum_nc(ijh,na,is,js)+w1 *  &
                                          CONJG(becp_nc(ikb,is,ibnd)) * &
                                                becp_nc(ikb,js,ibnd)
                                  END DO
                                  !
                               END DO
                               !
                            ELSE
                               !
                               becsum(ijh,na,current_spin) = &
                                        becsum(ijh,na,current_spin) + &
                                        w1 * DBLE( CONJG( becp(ikb,ibnd) ) * &
                                                          becp(ikb,ibnd) )
                               !
                            END IF                       
                            !
                            ijh = ijh + 1
                            !
                            DO jh = ( ih + 1 ), nh(np)
                               !
                               jkb = ijkb0 + jh
                               !
                               IF (noncolin) THEN
                                  !
                                  DO is=1,npol
                                     !
                                     DO js=1,npol
                                        becsum_nc(ijh,na,is,js) =         &
                                           becsum_nc(ijh,na,is,js) + w1 * &
                                           CONJG(becp_nc(ikb,is,ibnd)) *  &
                                                 becp_nc(jkb,js,ibnd)
                                     END DO
                                     !
                                  END DO
                                  !
                               ELSE
                                  !
                                  becsum(ijh,na,current_spin) = &
                                     becsum(ijh,na,current_spin) + w1 * 2.D0 * &
                                     DBLE( CONJG( becp(ikb,ibnd) ) * &
                                                  becp(jkb,ibnd) )
                               ENDIF
                               !            
                               ijh = ijh + 1
                               !
                            END DO
                            !
                         END DO
                         !
                         ijkb0 = ijkb0 + nh(np)
                         !
                      END IF
                      !
                   END DO
                   !
                ELSE
                   !
                   DO na = 1, nat
                      !
                      IF ( ityp(na) == np ) ijkb0 = ijkb0 + nh(np)
                      !
                   END DO
                   !
                END IF
                !
             END DO
             !
          END DO
          !
          CALL stop_clock( 'sum_band:becsum' )
          !
       END DO k_loop

       IF( use_tg ) THEN
          DEALLOCATE( tg_psi )
          DEALLOCATE( tg_rho )
       END IF

       IF (noncolin.and.okvan) THEN
          DO np = 1, ntyp
             IF ( upf(np)%tvanp ) THEN
                DO na = 1, nat
                   IF (ityp(na)==np) THEN
                      IF (so(np)) THEN
                         CALL transform_becsum_so(becsum_nc,becsum,na)
                      ELSE
                         CALL transform_becsum_nc(becsum_nc,becsum,na)
                      END IF
                   END IF
                END DO
             END IF
          END DO
       END IF
       !
       IF ( ALLOCATED (becsum_nc) ) DEALLOCATE( becsum_nc )
       !
       RETURN
       !
     END SUBROUTINE sum_band_k
     !
END SUBROUTINE sum_band
