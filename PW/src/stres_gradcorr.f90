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
  USE kinds,            ONLY : DP
  USE funct,            ONLY : gcxc, gcx_spin, gcc_spin, gcc_spin_more, &
                               dft_is_gradient, dft_is_meta, get_igcc, &
                               tau_xc, tau_xc_spin, get_meta
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  USE fft_types,        ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor), INTENT(in):: dfft
  INTEGER, INTENT(in) :: nspin
  real(DP), INTENT(in):: rho (dfft%nnr, nspin), rho_core (dfft%nnr), g(3,dfft%ngm), alat, omega
  real(DP), INTENT(inout):: kedtau(dfft%nnr, nspin) ! FIXME: should be intent(in)
  COMPLEX(DP), INTENT(in) :: rhog(dfft%ngm, nspin)
  COMPLEX(DP), INTENT(in) :: rhog_core(dfft%ngm)
  real(dp), INTENT(inout) :: sigmaxc(3, 3)
  !
  INTEGER :: k, l, m, ipol, is
  INTEGER :: nr1, nr2, nr3, nrxx, ngm
  real(DP) , ALLOCATABLE :: grho (:,:,:), rhoaux(:,:)
  COMPLEX(dp), ALLOCATABLE :: rhogaux(:,:)
  real(DP), PARAMETER :: epsr = 1.0d-6, epsg = 1.0d-10, e2 = 2.d0
  real(DP) :: grh2, grho2(2), sx, sc, v1x, v2x, v1c, v2c, &
       v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw, v2cup, v2cdw, v2cud, &
       zeta, rh, rup, rdw, grhoup, grhodw, grhoud, grup, grdw, &
       sigma_gradcorr(3, 3)
  LOGICAL :: igcc_is_lyp
  !dummy variables for meta-gga
  real(DP) :: v3x,v3c,v3xup,v3xdw,v3cup,v3cdw

  IF ( .not. dft_is_gradient() .and. .not. dft_is_meta() ) RETURN
  !
  IF ( nspin==4 ) CALL errore('stres_gradcorr', &
                    'noncollinear stress + GGA not implemented',1)
  IF ( dft_is_meta() .and. (nspin>1) )  CALL errore('stres_gradcorr', &
                    'Meta-GGA stress not yet implemented with spin polarization',1)

  igcc_is_lyp = (get_igcc() == 3)

  sigma_gradcorr(:,:) = 0.d0

  nr1 = dfft%nr1
  nr2 = dfft%nr2
  nr3 = dfft%nr3
  nrxx= dfft%nnr
  ngm = dfft%ngm
  !
  ALLOCATE (grho( 3, nrxx, nspin))
  ALLOCATE (rhogaux(ngm,nspin))
  !
  !    calculate the gradient of rho+rhocore in real space
  !    in LSDA case rho is temporarily converted in (up,down) format
  !
  IF ( nspin == 1 ) THEN
     rhogaux(:,1)  = rhog(:,1) + rhog_core(:)
  ELSE
     rhogaux(:,1)  = ( rhog(:,1) + rhog(:,2) + rhog_core(:) ) / 2.0_dp
     rhogaux(:,2)  = ( rhog(:,1) - rhog(:,2) + rhog_core(:) ) / 2.0_dp
  ENDIF
  DO is = 1, nspin
     CALL fft_gradient_g2r( dfft, rhogaux(1,is), g, grho(1,1,is) )
  ENDDO
  DEALLOCATE (rhogaux)
  !
  ALLOCATE (rhoaux(nrxx,nspin))
  IF ( nspin == 1 ) THEN
     rhoaux(:,1)  = rho(:,1) + rho_core(:)
  ELSE
     rhoaux(:,1)  = ( rho(:,1) + rho(:,2) + rho_core(:) ) / 2.0_dp
     rhoaux(:,2)  = ( rho(:,1) - rho(:,2) + rho_core(:) ) / 2.0_dp
  ENDIF

  IF (nspin==1) THEN
     !
     !    This is the LDA case
     !
     ! sigma_gradcor_{alpha,beta} ==
     !     omega^-1 \int (grad_alpha rho) ( D(rho*Exc)/D(grad_alpha rho) ) d3
     !
     DO k = 1, nrxx
        !
        grho2 (1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
        !
        IF (abs (rhoaux(k, 1) ) >epsr.and.grho2 (1) >epsg) THEN
           !
           ! routine computing v1x and v2x is different for GGA and meta-GGA
           ! FIXME : inefficient implementation
           !
           IF ( dft_is_meta() .and. get_meta() /= 4 ) THEN
              !
              kedtau(k,1) = kedtau(k,1) / e2
              CALL tau_xc (rhoaux(k,1), grho2(1),kedtau(k,1), sx, sc, v1x, v2x,v3x,v1c,v2c,v3c)
              kedtau(k,1) = kedtau(k,1) * e2
              !
           ELSE
              !
              CALL gcxc (rhoaux(k, 1), grho2(1), sx, sc, v1x, v2x, v1c, v2c)
              !
           ENDIF
           !
           DO l = 1, 3
              !
              DO m = 1, l
                 !
                 sigma_gradcorr (l, m) = sigma_gradcorr (l, m) + &
                                grho(l,k,1) * grho(m,k,1) * e2 * (v2x + v2c)
                 !
              ENDDO
              !
           ENDDO
           !
        ENDIF
        !
     ENDDO
     !
  ELSEIF (nspin == 2) THEN
     !
     !    This is the LSDA case
     !
     DO k = 1, nrxx
        grho2 (1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
        grho2 (2) = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
        !
        IF ( (abs (rhoaux(k, 1) ) >epsr.and.grho2 (1) >epsg) .and. &
             (abs (rhoaux(k, 2) ) >epsr.and.grho2 (2) >epsg) ) THEN
           !
           ! Spin polarization for metagga
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
              !
              CALL gcx_spin (rhoaux(k, 1), rhoaux(k, 2), grho2 (1), grho2 (2), &
                   sx, v1xup, v1xdw, v2xup, v2xdw)
              !
              rh = rhoaux(k, 1) + rhoaux(k, 2)
              !
              IF (rh>epsr) THEN
                 IF ( igcc_is_lyp ) THEN
                    rup = rhoaux(k, 1)
                    rdw = rhoaux(k, 2)
                    grhoup = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
                    grhodw = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
                    grhoud = grho(1,k,1) * grho(1,k,2) + &
                             grho(2,k,1) * grho(2,k,2) + &
                             grho(3,k,1) * grho(3,k,2)
                    CALL gcc_spin_more(rup, rdw, grhoup, grhodw, grhoud, sc, &
                                    v1cup, v1cdw, v2cup, v2cdw, v2cud)
                 ELSE
                    zeta = (rhoaux(k, 1) - rhoaux(k, 2) ) / rh

                    grh2 = (grho (1, k, 1) + grho (1, k, 2) ) **2 + &
                           (grho (2, k, 1) + grho (2, k, 2) ) **2 + &
                           (grho (3, k, 1) + grho (3, k, 2) ) **2
                    CALL gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
                    v2cup = v2c
                    v2cdw = v2c
                    v2cud = v2c
                 ENDIF
              ELSE
                 sc = 0.d0
                 v1cup = 0.d0
                 v1cdw = 0.d0
                 v2c = 0.d0
                 v2cup = 0.d0
                 v2cdw = 0.d0
                 v2cud = 0.d0
              ENDIF
              !
           ENDIF
           !
           DO l = 1, 3
              DO m = 1, l
                 !    exchange
                 sigma_gradcorr (l, m) = sigma_gradcorr (l, m) + &
                      grho (l, k, 1) * grho (m, k, 1) * e2 * v2xup + &
                      grho (l, k, 2) * grho (m, k, 2) * e2 * v2xdw
                 !    correlation
                 sigma_gradcorr (l, m) = sigma_gradcorr (l, m) + &
                     ( grho (l, k, 1) * grho (m, k, 1) * v2cup + &
                       grho (l, k, 2) * grho (m, k, 2) * v2cdw + &
                      (grho (l, k, 1) * grho (m, k, 2) +         &
                       grho (l, k, 2) * grho (m, k, 1) ) * v2cud ) * e2
              ENDDO
              !
           ENDDO
           !
        ENDIF
        !
     ENDDO
     !
  ENDIF
  DEALLOCATE (rhoaux)
  DEALLOCATE(grho)
  !
  DO l = 1, 3
     DO m = 1, l - 1
        sigma_gradcorr (m, l) = sigma_gradcorr (l, m)
     ENDDO

  ENDDO
  CALL mp_sum(  sigma_gradcorr, intra_bgrp_comm )

  sigmaxc(:,:) = sigmaxc(:,:) + sigma_gradcorr(:,:) / &
                                (nr1 * nr2 * nr3)
  !
  RETURN

END SUBROUTINE stres_gradcorr

