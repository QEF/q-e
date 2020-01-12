!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stres_gradcorr( rho, rhog, rho_core, rhog_core, kedtau, nspin, &
                           dfft, g, alat, omega, sigmaxc )
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY: DP
  USE funct,            ONLY: dft_is_gradient, dft_is_meta, get_meta
  USE xc_gga,           ONLY: xc_gcx
  USE xc_mgga,          ONLY: xc_metagcx
  USE spin_orb,         ONLY: domag
  USE mp_bands,         ONLY: intra_bgrp_comm
  USE mp,               ONLY: mp_sum
  USE fft_types,        ONLY: fft_type_descriptor
  USE fft_rho,          ONLY: rho_r2g
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor), INTENT(IN):: dfft
  INTEGER,  INTENT(IN) :: nspin
  REAL(DP), INTENT(IN) :: rho(dfft%nnr,nspin), rho_core(dfft%nnr)
  REAL(DP), INTENT(IN) :: g(3,dfft%ngm), alat, omega
  REAL(DP), INTENT(INOUT) :: kedtau(dfft%nnr, nspin) ! FIXME: should be INTENT(IN)
  COMPLEX(DP), INTENT(IN) :: rhog(dfft%ngm, nspin)
  COMPLEX(DP), INTENT(IN) :: rhog_core(dfft%ngm)
  REAL(DP), INTENT(INOUT) :: sigmaxc(3,3)
  !
  ! ... local variables
  !
  INTEGER :: k, l, m, ipol, ir, ig, is, nspin0, np
  INTEGER :: nr1, nr2, nr3, nrxx, ngm
  REAL(DP), ALLOCATABLE :: grho(:,:,:), grho2(:,:), rhoaux(:,:), segni(:)
  COMPLEX(DP), ALLOCATABLE :: rhogaux(:,:)
  !
  REAL(DP), ALLOCATABLE :: sx(:), sc(:)
  REAL(DP), ALLOCATABLE :: v1x(:,:), v2x(:,:), v3x(:,:), rhos(:)
  REAL(DP), ALLOCATABLE :: v1c(:,:), v2c(:,:,:), v3c(:,:), v2c_ud(:)
  !
  REAL(DP), PARAMETER :: epsr = 1.0d-6, epsg = 1.0d-10, e2 = 2.d0
  REAL(DP) :: sigma_gradcorr(3, 3)
  !
  !
  IF ( .NOT. dft_is_gradient() .AND. .NOT. dft_is_meta() ) RETURN
  !
  IF ( dft_is_meta() .and. nspin>1 )  CALL errore('stres_gradcorr', &
       'Meta-GGA stress does not work with spin polarization',1)
  !
  np = 1
  IF ( nspin==2 .AND. dft_is_meta() ) np=3
  !
  nspin0 = nspin
  IF (nspin==4) nspin0 = 1
  IF (nspin==4.and.domag) nspin0 = 2
  !
  sigma_gradcorr(:,:) = 0.0_DP
  !
  nr1 = dfft%nr1
  nr2 = dfft%nr2
  nr3 = dfft%nr3
  nrxx= dfft%nnr
  ngm = dfft%ngm
  !
  ALLOCATE( grho(3,nrxx,nspin0) )
  ALLOCATE (rhoaux(nrxx,nspin0))
  ALLOCATE( rhogaux(ngm,nspin0) )
  !
  ! calculate the gradient of rho+rhocore in real space
  ! For convenience rhoaux is in (up,down) format
  !
  IF ( nspin0 == 1 ) THEN
     !
     rhoaux(:,1)  = rho(:,1)  + rho_core(:)
     rhogaux(:,1) = rhog(:,1) + rhog_core(:)
     !
  ELSE IF ( nspin0 == 2 ) THEN
     !
     IF ( nspin == 4 .AND. domag ) THEN
        ALLOCATE( segni( nrxx ) )
        CALL compute_rho( rho, rhoaux, segni, nrxx )
        DEALLOCATE( segni )
        rhoaux(:,1) = rhoaux(:,1) + rho_core(:) / 2.0_DP
        rhoaux(:,2) = rhoaux(:,2) + rho_core(:) / 2.0_DP
        CALL rho_r2g ( dfft, rhoaux(:,1:nspin0), rhogaux(:,1:nspin0) )
     ELSE
        rhoaux(:,1)  = ( rho(:,1) + rho(:,2) + rho_core(:) ) / 2.0_DP
        rhoaux(:,2)  = ( rho(:,1) - rho(:,2) + rho_core(:) ) / 2.0_DP
        rhogaux(:,1)  = ( rhog(:,1) + rhog(:,2) + rhog_core(:) ) / 2.0_DP
        rhogaux(:,2)  = ( rhog(:,1) - rhog(:,2) + rhog_core(:) ) / 2.0_DP
     END IF
  ENDIF
  !
  DO is = 1, nspin0
     CALL fft_gradient_g2r( dfft, rhogaux(1,is), g, grho(1,1,is) )
  ENDDO
  DEALLOCATE (rhogaux)
  !
  ALLOCATE( grho2(nrxx,nspin0)  )
  ALLOCATE( v1x(nrxx,nspin0), v2x(nrxx,nspin0), v3x(nrxx,nspin0) )
  ALLOCATE( v1c(nrxx,nspin0), v2c(np,nrxx,nspin0), v3c(nrxx,nspin0) )
  ALLOCATE( sx(nrxx), sc(nrxx) )
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
     grho2(:,1) = grho(1,:,1)**2 + grho(2,:,1)**2 + grho(3,:,1)**2
     !
     IF ( dft_is_meta() .AND. get_meta() /= 4 ) THEN
        kedtau(:,1) = kedtau(:,1) / e2
        CALL xc_metagcx( nrxx, 1, np, rhoaux, grho, kedtau, sx, sc, &
                         v1x, v2x, v3x, v1c, v2c, v3c )
        kedtau(:,1) = kedtau(:,1) * e2
     ELSE
        CALL xc_gcx( nrxx, nspin0, rhoaux, grho, sx, sc, v1x, v2x, v1c, v2c(1,:,:) )
     ENDIF
     !
     DO l = 1, 3
       DO m = 1, l
         sigma_gradcorr(l,m) = sigma_gradcorr(l,m) + SUM( grho(l,:,1)*grho(m,:,1)*   &
                                                          e2 * (v2x(:,1) + v2c(1,:,1)) )
       ENDDO
     ENDDO
     !
  ELSEIF (nspin0 == 2) THEN
     !
     !    Spin-polarized case
     !
     grho2(:,:) = grho(1,:,:)**2 + grho(2,:,:)**2 + grho(3,:,:)**2
     !
     IF ( dft_is_meta() ) THEN
        !
        kedtau(:,1:nspin0) = kedtau(:,1:nspin0) / e2
        CALL xc_metagcx( nrxx, nspin0, np, rhoaux, grho, kedtau, sx, sc, &
                         v1x, v2x, v3x, v1c, v2c, v3c )
        kedtau(:,1:nspin0) = kedtau(:,1:nspin0) * e2
        ! FIXME : what are we supposed to do now?
        !
     ELSE
        !
        ALLOCATE( v2c_ud(nrxx) )
        !
        CALL xc_gcx( nrxx, nspin0, rhoaux, grho, sx, sc, v1x, v2x, v1c, v2c(1,:,:), v2c_ud )
        !
        DO l = 1, 3
           DO m = 1, l
              !
              ! ... exchange
              sigma_gradcorr(l,m) = &
                  SUM( grho(l,:,1) * grho(m,:,1) * e2 * v2x(:,1) + &
                       grho(l,:,2) * grho(m,:,2) * e2 * v2x(:,2) )
              !
              ! ... correlation
              sigma_gradcorr(l,m) = sigma_gradcorr(l,m)     +     &
                  SUM( grho(l,:,1) * grho(m,:,1) * v2c(1,:,1) +     &
                       grho(l,:,2) * grho(m,:,2) * v2c(1,:,2) +     &
                      (grho(l,:,1) * grho(m,:,2) +                &
                       grho(l,:,2) * grho(m,:,1)) * v2c_ud(:) ) * e2
           ENDDO
        ENDDO
        !
        DEALLOCATE( v2c_ud )
        !
     ENDIF
     !
  ENDIF
  !
  DEALLOCATE( sc, sx )
  DEALLOCATE( v1c, v2c, v3c )
  DEALLOCATE( v1x, v2x, v3x )
  DEALLOCATE( grho2  )
  DEALLOCATE( rhoaux )
  DEALLOCATE( grho   )
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
  RETURN
  !
END SUBROUTINE stres_gradcorr

