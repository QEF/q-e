!
! Copyright (C) 2001-2005 PWSCF group
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
  USE wvfct,                ONLY : gamma_only
  USE control_flags,        ONLY : diago_full_acc
  USE cell_base,            ONLY : at, bg, omega, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   ngm, g, nl, nlm
  USE gsmooth,              ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                                   nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk
  USE ldaU,                 ONLY : lda_plus_U
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE realus,               ONLY : tqr
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
                                   root_image, npool, my_pool_id
  USE mp,                   ONLY : mp_bcast
  USE funct,                ONLY : dft_is_meta
  USE rad_paw_routines,     ONLY : PAW_symmetrize
  USE grid_paw_variables,   ONLY : okpaw
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
  CALL poolreduce( 1, eband )
  !
  ! ... symmetrization of the charge density (and local magnetization)
  !
#if defined (__PARA)
  !
  ! ... reduce charge density across pools
  !
  CALL poolreduce( nspin * nrxx, rho%of_r )
  if (dft_is_meta() ) CALL poolreduce( nspin * nrxx, rho%kin_r )
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
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1, w2
         ! weights
       REAL(DP), ALLOCATABLE :: becp(:,:)
         ! contains <beta|psi>
       !
       !
       ALLOCATE( becp( nkb, nbnd ) )
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       IF ( nks > 1 ) REWIND( iunigk )
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
          DO ibnd = 1, nbnd, 2
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
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( .NOT. okvan ) CYCLE k_loop
          !
          IF ( nkb > 0 ) &
             CALL ccalbec( nkb, npwx, npw, nbnd, becp, vkb, evc )
          !
          CALL start_clock( 'becsum' )
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
                                            w1 * becp(ikb,ibnd) * becp(ikb,ibnd)
                            !
                            ijh = ijh + 1
                            !
                            DO jh = ( ih + 1 ), nh(np)
                               !
                               jkb = ijkb0 + jh
                               !
                               becsum(ijh,na,current_spin) = &
                                     becsum(ijh,na,current_spin) + &
                                     w1 * 2.D0 * becp(ikb,ibnd) * becp(jkb,ibnd)
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
          CALL stop_clock( 'becsum' )
          !
       END DO k_loop
       !
       DEALLOCATE( becp )
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
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1
       ! weights
       COMPLEX(DP), ALLOCATABLE :: becp(:,:), becp_nc(:,:,:)
       ! contains <beta|psi>
       !
       COMPLEX(DP), ALLOCATABLE :: becsum_nc(:,:,:,:)
       !
       INTEGER :: ipol, kh, kkb, is1, is2, js
       !

       IF (okvan) THEN
          IF (noncolin) THEN
             ALLOCATE(becsum_nc(nhm*(nhm+1)/2,nat,npol,npol))
             becsum_nc=(0.d0, 0.d0)
             ALLOCATE( becp_nc( nkb, npol, nbnd ) )
          ELSE
             ALLOCATE( becp( nkb, nbnd ) )
          END IF
       ENDIF
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       IF ( nks > 1 ) REWIND( iunigk )
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
          DO ibnd = 1, nbnd
             !
             eband = eband + et(ibnd,ik) * wg(ibnd,ik)
             !
             ! ... the sum of eband and demet is the integral for e < ef of 
             ! ... e n(e) which reduces for degauss=0 to the sum of the 
             ! ... eigenvalues 
             w1 = wg(ibnd,ik) / omega
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
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( .NOT. okvan ) CYCLE k_loop
          !
          IF (noncolin) THEN
             IF ( nkb > 0 ) &
                CALL ccalbec_nc( nkb, npwx, npw, npol, nbnd, &
                                                 becp_nc, vkb, evc )
          ELSE
             IF ( nkb > 0 ) &
                CALL ccalbec( nkb, npwx, npw, nbnd, becp, vkb, evc )
          ENDIF
          !
          CALL start_clock( 'becsum' )
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
          CALL stop_clock( 'becsum' )
          !
       END DO k_loop

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
       IF (okvan) THEN
          IF (noncolin) THEN
             DEALLOCATE( becsum_nc )
             DEALLOCATE( becp_nc )
          ELSE
             DEALLOCATE( becp )
          ENDIF
       END IF
       !
       RETURN
       !
     END SUBROUTINE sum_band_k
     !
END SUBROUTINE sum_band
