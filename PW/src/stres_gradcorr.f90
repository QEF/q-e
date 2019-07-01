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
  USE funct,            ONLY: dft_is_gradient, dft_is_meta, get_igcc,  &
                              tau_xc, tau_xc_spin, get_meta
  USE xc_gga,           ONLY: xc_gcx !gcxc, gcx_spin, gcc_spin, gcc_spin_more
  USE mp_bands,         ONLY: intra_bgrp_comm
  USE mp,               ONLY: mp_sum
  USE fft_types,        ONLY: fft_type_descriptor
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
  INTEGER :: k, l, m, ipol, is, nspin0
  INTEGER :: nr1, nr2, nr3, nrxx, ngm
  REAL(DP), ALLOCATABLE :: grho(:,:,:), rhoaux(:,:)
  REAL(DP), ALLOCATABLE :: grho2(:,:), grho_ud(:),grhor2(:)
  COMPLEX(DP), ALLOCATABLE :: rhogaux(:,:)
  !^^^
  REAL(DP), ALLOCATABLE :: zeta(:), rh(:)
  REAL(DP) :: sx(dfft%nnr), sc(dfft%nnr)
  REAL(DP), ALLOCATABLE :: v1x(:,:), v2x(:,:), rhos(:)
  REAL(DP), ALLOCATABLE :: v1c(:,:), v2c(:,:), v2c_ud(:)
  INTEGER,  ALLOCATABLE :: null_v(:)
  REAL(DP) :: vnull
  !
  REAL(DP), PARAMETER :: epsr = 1.0d-6, epsg = 1.0d-10, e2 = 2.d0
  REAL(DP) :: sigma_gradcorr(3, 3)
  LOGICAL :: igcc_is_lyp
  !
  ! ... dummy variables for meta-gga
  REAL(DP) :: v3x, v3c, v3xup, v3xdw, v3cup, v3cdw
  !
  !
  IF ( .NOT. dft_is_gradient() .AND. .NOT. dft_is_meta() ) RETURN
  !
  nspin0 = nspin
  !if (nspin==4) nspin0 = 1
  !if (nspin==4.and.domag) nspin0 = 2
  IF ( nspin==4 ) CALL errore('stres_gradcorr', &
                    'noncollinear stress + GGA not implemented',1)
  IF ( dft_is_meta() .and. (nspin>1) )  CALL errore('stres_gradcorr', &
                    'Meta-GGA stress not yet implemented with spin polarization',1)

  igcc_is_lyp = (get_igcc() == 3)

  sigma_gradcorr(:,:) = 0.0_DP

  nr1 = dfft%nr1
  nr2 = dfft%nr2
  nr3 = dfft%nr3
  nrxx= dfft%nnr
  ngm = dfft%ngm
  !
  ALLOCATE( grho(3,nrxx,nspin0) )
  ALLOCATE( grho2(nrxx,nspin0)  )
  ALLOCATE( rhogaux(ngm,nspin0) )
  !
  ALLOCATE( v1x(nrxx,nspin0), v2x(nrxx,nspin0) )
  ALLOCATE( v1c(nrxx,nspin0), v2c(nrxx,nspin0) )
  !
  ! calculate the gradient of rho+rhocore in real space
  ! in LSDA case rho is temporarily converted in (up,down) format
  !
  IF ( nspin == 1 ) THEN
     rhogaux(:,1)  = rhog(:,1) + rhog_core(:)
  ELSE
     rhogaux(:,1)  = ( rhog(:,1) + rhog(:,2) + rhog_core(:) ) / 2.0_DP
     rhogaux(:,2)  = ( rhog(:,1) - rhog(:,2) + rhog_core(:) ) / 2.0_DP
  ENDIF
  !
  DO is = 1, nspin
     CALL fft_gradient_g2r( dfft, rhogaux(1,is), g, grho(1,1,is) )
  ENDDO
  DEALLOCATE (rhogaux)
  !
  !
  ALLOCATE (rhoaux(nrxx,nspin))
  IF ( nspin == 1 ) THEN
     rhoaux(:,1)  = rho(:,1) + rho_core(:)
  ELSE
     rhoaux(:,1)  = ( rho(:,1) + rho(:,2) + rho_core(:) ) / 2.0_DP
     rhoaux(:,2)  = ( rho(:,1) - rho(:,2) + rho_core(:) ) / 2.0_DP
  ENDIF
  !
  !
  IF (nspin==1) THEN
     !
     !    This is the LDA case
     !
     ! sigma_gradcor_{alpha,beta} ==
     !     omega^-1 \int (grad_alpha rho) ( D(rho*Exc)/D(grad_alpha rho) ) d3
     !
     ALLOCATE( null_v(nrxx) )     
     !
     ! routine computing v1x_v and v2x_v is different for GGA and meta-GGA
     ! FIXME : inefficient implementation
     !
     null_v = 1
     grho2(:,1) = grho(1,:,1)**2 + grho(2,:,1)**2 + grho(3,:,1)**2
     !
     WHERE( .NOT. (ABS(rhoaux(:,1))>epsr .AND. grho2(:,1)>epsg) ) null_v(:)  = 0
     !
     IF ( .NOT. (dft_is_meta() .AND. get_meta() /= 4) ) &
       CALL xc_gcx( nrxx, nspin, rhoaux, grho, sx, sc, v1x, v2x, v1c, v2c )
     !
     DO k = 1, nrxx
        IF ( dft_is_meta() .AND. get_meta() /= 4 .AND. null_v(k) /= 0 ) THEN
           IF ( ABS(rhoaux(k,1))>epsr .AND. grho2(k,1)>epsg ) THEN
              !
              kedtau(k,1) = kedtau(k,1) / e2
              CALL tau_xc( rhoaux(k,1), grho2(k,1), kedtau(k,1), sx(k), sc(k), &
                                 v1x(k,1), v2x(k,1), v3x, v1c(k,1), v2c(k,1), v3c )
              kedtau(k,1) = kedtau(k,1) * e2
              !
           ENDIF
        ENDIF
        !
        DO l = 1, 3
           DO m = 1, l
              sigma_gradcorr(l,m) = sigma_gradcorr(l,m) + grho(l,k,1)*grho(m,k,1)* &
                                    e2 * (v2x(k,1) + v2c(k,1)) !* null_v(k)
           ENDDO
        ENDDO
        !
     ENDDO
     !
     DEALLOCATE( null_v )
     !
     !
  ELSEIF (nspin == 2) THEN
     !
     !    This is the LSDA case
     !
     grho2(:,:) = grho(1,:,:)**2 + grho(2,:,:)**2 + grho(3,:,:)**2
     !
     IF ( dft_is_meta() ) THEN
        !
        ! Not working with spin polarization ... FIXME
        ! call tau_xc_spin (rho(k,1), rho(k,2), grho2(1), grho2(2), &
        !            kedtau(k,1), kedtau(k,2), sx, sc, &
        !            v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw, &
        !            v1cup, v1cdw, v2cup, v2cdw, v3cup, v3cdw )
        !
     ELSE
        !
        ALLOCATE( v2c_ud(dfft%nnr) )
        !
        CALL xc_gcx( nrxx, nspin, rhoaux, grho, sx, sc, v1x, v2x, v1c, v2c, v2c_ud )
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
                  SUM( grho(l,:,1) * grho(m,:,1) * v2c(:,1) +     &
                       grho(l,:,2) * grho(m,:,2) * v2c(:,2) +     &
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
  DEALLOCATE( rhoaux )
  DEALLOCATE( grho   )
  DEALLOCATE( grho2  )
  !
  DEALLOCATE( v1x, v2x )
  DEALLOCATE( v1c, v2c )
  !
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

