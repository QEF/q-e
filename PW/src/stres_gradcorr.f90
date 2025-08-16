!
! Copyright (C) 2001-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stres_gradcorr( rho, rho_core, rhog_core, nspin, domag, &
                           dfft, g, alat, omega, sigmaxc )
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY: DP
  USE xc_lib,           ONLY: xclib_dft_is, xclib_get_id, xc_gcx, xc_metagcx
  USE scf,              ONLY: scf_type
  USE mp_bands,         ONLY: intra_bgrp_comm
  USE mp,               ONLY: mp_sum
  USE fft_types,        ONLY: fft_type_descriptor
  USE fft_rho,          ONLY: rho_r2g
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor), INTENT(IN):: dfft
  INTEGER,  INTENT(IN) :: nspin
  LOGICAL,  INTENT(IN) :: domag
  TYPE(scf_type), INTENT(IN) :: rho
  REAL(DP), INTENT(IN) :: rho_core(dfft%nnr)
  REAL(DP), INTENT(IN) :: g(3,dfft%ngm)
  REAL(DP), INTENT(IN) :: alat, omega
  COMPLEX(DP), INTENT(IN) :: rhog_core(dfft%ngm)
  REAL(DP), INTENT(INOUT) :: sigmaxc(3,3)
  !
  ! ... local variables
  !
  INTEGER :: k, l, m, ipol, ir, ig, is, nspin0, np
  INTEGER :: nr1, nr2, nr3, nrxx, ngm
  REAL(DP), ALLOCATABLE :: grho(:,:,:), grho2(:,:), rhoaux(:,:), &
                           segni(:), kedtaue2(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogaux(:,:)
  !
  REAL(DP), ALLOCATABLE :: sx(:), sc(:)
  REAL(DP), ALLOCATABLE :: v1x(:,:), v2x(:,:), v3x(:,:)
  REAL(DP), ALLOCATABLE :: v1c(:,:), v2c(:,:), v3c(:,:), v2c_ud(:), v2cm(:,:,:)
  !
  REAL(DP), PARAMETER :: epsr = 1.0d-6, epsg = 1.0d-10, e2 = 2.d0
  REAL(DP) :: v2xc, v2xc_uu, v2xc_dd
  REAL(DP) :: sigma_gc11, sigma_gc31, sigma_gc21, &
              sigma_gc32, sigma_gc22, sigma_gc33
  REAL(DP) :: sigma_gradcorr(3,3)
  !
  IF ( .NOT. xclib_dft_is('gradient') .AND. .NOT. xclib_dft_is('meta') ) RETURN
  !
  IF ( xclib_dft_is('meta') .AND. nspin>1 )  CALL errore( 'stres_gradcorr', &
       'Meta-GGA stress does not work with spin polarization', 1 )
  !
  !$acc data present_or_copyin( rho, rho%of_r, rho%of_g, rho_core, rhog_core, g )
  !
  np = 1
  IF ( nspin==2 .AND. xclib_dft_is('meta') ) np = 3
  !
  nspin0 = nspin
  IF (nspin==4) nspin0 = 1
  IF (nspin==4.AND.domag) nspin0 = 2
  !
  sigma_gc11=0.0_DP ; sigma_gc21=0.0_DP ; sigma_gc22=0.0_DP
  sigma_gc31=0.0_DP ; sigma_gc32=0.0_DP ; sigma_gc33=0.0_DP
  !
  nr1 = dfft%nr1
  nr2 = dfft%nr2
  nr3 = dfft%nr3
  nrxx= dfft%nnr
  ngm = dfft%ngm
  !
  ALLOCATE( grho(3,nrxx,nspin0) )
  ALLOCATE( rhoaux(nrxx,nspin0) )
  ALLOCATE( rhogaux(ngm,nspin0) )
  IF (xclib_dft_is('meta')) ALLOCATE( kedtaue2(dfft%nnr,nspin) )
  !$acc data create( grho, rhoaux )
  !$acc data create( rhogaux )
  !
  ! calculate the gradient of rho+rhocore in real space
  ! For convenience rhoaux is in (up,down) format
  !
  IF ( nspin0 == 1 ) THEN
     !
     !$acc parallel loop
     DO k = 1, nrxx
       rhoaux(k,1) = rho%of_r(k,1) + rho_core(k)
     ENDDO
     !$acc parallel loop
     DO k = 1, ngm
       rhogaux(k,1) = rho%of_g(k,1) + rhog_core(k)
     ENDDO
     !
  ELSEIF ( nspin0 == 2 ) THEN
     !
     IF ( nspin == 4 .AND. domag ) THEN
        !
        ALLOCATE( segni( nrxx ) )
        !
        !$acc data copyout( segni )
        CALL compute_rho( rho%of_r, rhoaux, segni, nrxx )
        !$acc parallel loop
        DO k = 1, nrxx
          rhoaux(k,1) = rhoaux(k,1) + rho_core(k) / 2.0_DP
          rhoaux(k,2) = rhoaux(k,2) + rho_core(k) / 2.0_DP
        ENDDO
        CALL rho_r2g( dfft, rhoaux(:,1:nspin0), rhogaux(:,1:nspin0) )
        !$acc end data
     ELSE
        !$acc parallel loop
        DO k = 1, nrxx
          rhoaux(k,1)  = ( rho%of_r(k,1) + rho%of_r(k,2) + rho_core(k) ) / 2.0_DP
          rhoaux(k,2)  = ( rho%of_r(k,1) - rho%of_r(k,2) + rho_core(k) ) / 2.0_DP
        ENDDO
        !$acc parallel loop
        DO k = 1, ngm
          rhogaux(k,1)  = ( rho%of_g(k,1) + rho%of_g(k,2) + rhog_core(k) ) / 2.0_DP
          rhogaux(k,2)  = ( rho%of_g(k,1) - rho%of_g(k,2) + rhog_core(k) ) / 2.0_DP
        ENDDO
     ENDIF
  ENDIF
  !
  DO is = 1, nspin0
    CALL fft_gradient_g2r( dfft, rhogaux(:,is), g, grho(:,:,is) )
  ENDDO
  !
  !$acc end data
  DEALLOCATE( rhogaux )
  !
  ALLOCATE( grho2(nrxx,nspin0)  )
  ALLOCATE( v1x(nrxx,nspin0), v2x(nrxx,nspin0) )
  ALLOCATE( v1c(nrxx,nspin0), v2c(nrxx,nspin0) )
  ALLOCATE( sx(nrxx), sc(nrxx) )
  !
  IF ( xclib_dft_is('meta') ) &
        ALLOCATE( v2cm(np,nrxx,nspin0), v3x(nrxx,nspin0), v3c(nrxx,nspin0) )
  !$acc data create( grho2, sx, sc, v1x, v2x, v1c, v2c )
  !
  IF (nspin0==1) THEN
     !
     !    Spin-unpolarized case
     !
     ! sigma_gradcor_{alpha,beta} ==
     !     omega^-1 \int (grad_alpha rho) ( D(rho*Exc)/D(grad_alpha rho) ) d3
     !
     ! routine computing v1x_v and v2x_v is different for GGA and meta-GGA
     !
     !$acc parallel loop
     DO k = 1, nrxx
       grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
     ENDDO
     !
     IF ( xclib_dft_is('meta') .AND. xclib_get_id('MGGA','EXCH') /= 4 ) THEN
        !$acc data present_or_copyin(rho%kin_r) create( kedtaue2, v2cm, v3x, v3c )
        !$acc parallel loop
        DO k = 1, nrxx
          kedtaue2(k,1) = rho%kin_r(k,1) / e2
        ENDDO
        CALL xc_metagcx( nrxx, 1, np, rhoaux, grho, kedtaue2, sx, sc, &
                         v1x, v2x, v3x, v1c, v2cm, v3c, gpu_args_=.TRUE. )
        !$acc parallel loop
        DO k = 1, nrxx
          v2c(k,1) = v2cm(1,k,1)
        ENDDO
        !$acc end data
     ELSE
        CALL xc_gcx( nrxx, nspin0, rhoaux, grho, sx, sc, v1x, v2x, v1c, v2c, gpu_args_=.TRUE. )
     ENDIF
     !
     !$acc parallel loop reduction(+:sigma_gc11,sigma_gc21,sigma_gc22, &
     !$acc&                          sigma_gc31,sigma_gc32,sigma_gc33)
     DO k = 1, nrxx
       v2xc = e2 * (v2x(k,1) + v2c(k,1))
       sigma_gc11 = sigma_gc11 + grho(1,k,1)*grho(1,k,1) * v2xc
       sigma_gc21 = sigma_gc21 + grho(2,k,1)*grho(1,k,1) * v2xc
       sigma_gc22 = sigma_gc22 + grho(2,k,1)*grho(2,k,1) * v2xc
       sigma_gc31 = sigma_gc31 + grho(3,k,1)*grho(1,k,1) * v2xc
       sigma_gc32 = sigma_gc32 + grho(3,k,1)*grho(2,k,1) * v2xc
       sigma_gc33 = sigma_gc33 + grho(3,k,1)*grho(3,k,1) * v2xc
     ENDDO
     !
  ELSEIF (nspin0 == 2) THEN
     !
     !    Spin-polarized case
     !
     !$acc parallel loop
     DO k = 1, nrxx
       grho2(k,1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
       grho2(k,2) = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
     ENDDO
     !
     IF ( xclib_dft_is('meta') ) THEN
        !
        !$acc data present_or_copyin(rho%kin_r) create( kedtaue2, v2cm, v3x, v3c )
        !$acc parallel loop
        DO k = 1, nrxx
          kedtaue2(k,1:nspin0) = rho%kin_r(k,1:nspin0) / e2
        ENDDO
        CALL xc_metagcx( nrxx, nspin0, np, rhoaux, grho, kedtaue2, sx, sc, &
                         v1x, v2x, v3x, v1c, v2cm, v3c, gpu_args_=.TRUE. )
        !$acc parallel loop
        DO k = 1, nrxx
          v2c(k,:) = v2cm(1,k,:)
        ENDDO
        !$acc end data
        ! FIXME : what are we supposed to do now?
        !
     ELSE
        !
        ALLOCATE( v2c_ud(nrxx) )
        !$acc data create( v2c_ud )
        !
        CALL xc_gcx( nrxx, nspin0, rhoaux, grho, sx, sc, v1x, v2x, v1c, v2c, v2c_ud, gpu_args_=.TRUE. )
        !
        !$acc parallel loop reduction(+:sigma_gc11,sigma_gc21,sigma_gc22, &
        !$acc&                          sigma_gc31,sigma_gc32,sigma_gc33)
        DO k = 1, nrxx
          !
          v2xc_uu = e2 * (v2x(k,1)+v2c(k,1))
          v2xc_dd = e2 * (v2x(k,2)+v2c(k,2))
          !
          sigma_gc11 = sigma_gc11 + grho(1,k,1)*grho(1,k,1) * v2xc_uu + &
                                    grho(1,k,2)*grho(1,k,2) * v2xc_dd + &
                                   (grho(1,k,1)*grho(1,k,2) + &
                                    grho(1,k,2)*grho(1,k,1)) * v2c_ud(k) * e2
          sigma_gc21 = sigma_gc21 + grho(2,k,1)*grho(1,k,1) * v2xc_uu + &
                                    grho(2,k,2)*grho(1,k,2) * v2xc_dd + &
                                   (grho(2,k,1)*grho(1,k,2) + &
                                    grho(1,k,2)*grho(2,k,1)) * v2c_ud(k) * e2
          sigma_gc22 = sigma_gc22 + grho(2,k,1)*grho(2,k,1) * v2xc_uu + &
                                    grho(2,k,2)*grho(2,k,2) * v2xc_dd + &
                                   (grho(2,k,1)*grho(2,k,2) + &
                                    grho(2,k,2)*grho(2,k,1)) * v2c_ud(k) * e2
          sigma_gc31 = sigma_gc31 + grho(3,k,1)*grho(1,k,1) * v2xc_uu + &
                                    grho(3,k,2)*grho(1,k,2) * v2xc_dd + &
                                   (grho(3,k,1)*grho(1,k,2) + &
                                    grho(1,k,2)*grho(3,k,1)) * v2c_ud(k) * e2
          sigma_gc32 = sigma_gc32 + grho(3,k,1)*grho(2,k,1) * v2xc_uu + &
                                    grho(3,k,2)*grho(2,k,2) * v2xc_dd + &
                                   (grho(3,k,1)*grho(2,k,2) + &
                                    grho(2,k,2)*grho(3,k,1)) * v2c_ud(k) * e2
          sigma_gc33 = sigma_gc33 + grho(3,k,1)*grho(3,k,1) * v2xc_uu + &
                                    grho(3,k,2)*grho(3,k,2) * v2xc_dd + &
                                   (grho(3,k,1)*grho(3,k,2) + &
                                    grho(3,k,2)*grho(3,k,1)) * v2c_ud(k) * e2
        ENDDO
        !
        !$acc end data
        DEALLOCATE( v2c_ud )
        !
     ENDIF
     !
  ENDIF
  !
  sigma_gradcorr(1,1) = sigma_gc11
  sigma_gradcorr(2,1) = sigma_gc21
  sigma_gradcorr(2,2) = sigma_gc22
  sigma_gradcorr(3,1) = sigma_gc31
  sigma_gradcorr(3,2) = sigma_gc32
  sigma_gradcorr(3,3) = sigma_gc33 
  !
  !$acc end data
  !$acc end data
  DEALLOCATE( sc, sx )
  DEALLOCATE( v1c, v2c )
  DEALLOCATE( v1x, v2x )
  DEALLOCATE( grho, grho2  )
  DEALLOCATE( rhoaux )
  IF (xclib_dft_is('meta')) DEALLOCATE( kedtaue2, v2cm, v3x, v3c )
  !
  DO l = 1, 3
     DO m = 1, l - 1
        sigma_gradcorr(m,l) = sigma_gradcorr(l,m)
     ENDDO
  ENDDO
  !
  CALL mp_sum( sigma_gradcorr, intra_bgrp_comm )
  !
  sigmaxc(:,:) = sigmaxc(:,:) + sigma_gradcorr(:,:) / (nr1 * nr2 * nr3)
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE stres_gradcorr

