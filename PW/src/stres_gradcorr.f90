!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
subroutine stres_gradcorr( rho, rhog, rho_core, rhog_core, kedtau, nspin, &
                           dfft, g, alat, omega, sigmaxc )
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : noncolin
  USE spin_orb,         ONLY : domag
  use funct,            ONLY : gcxc, gcx_spin, gcc_spin, gcc_spin_more, &
                               dft_is_gradient, dft_is_meta, get_igcc, &
                               tau_xc, tau_xc_spin
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  USE fft_types,        ONLY : fft_type_descriptor
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE wavefunctions_module, ONLY : psic
  !
  IMPLICIT NONE
  !
  TYPE(fft_type_descriptor), INTENT(IN):: dfft
  integer, intent(in) :: nspin
  real(DP), intent(in):: rho (dfft%nnr, nspin)
  real(DP), intent(inout):: kedtau(dfft%nnr, nspin)
  ! FIXME: should be intent(in)
  real(dp), intent(in) :: rho_core (dfft%nnr), g(3,dfft%ngm), alat, omega
  complex(DP), intent(in) :: rhog(dfft%ngm, nspin)
  ! FIXME: should be intent(in)
  complex(DP), intent(in) :: rhog_core(dfft%ngm)
  real(dp), intent(inout) :: sigmaxc (3, 3)
  !
  integer :: k, l, m, ipol, is, nspin0
  integer :: nr1, nr2, nr3, nrxx, ngm
  real(DP) , allocatable :: grho (:,:,:), segni(:), rhoout(:,:)
  complex(DP) , allocatable :: rhogsum(:,:)
  real(DP), parameter :: epsr = 1.0d-6, epsg = 1.0d-10, e2 = 2.d0
  real(DP) :: grh2, grho2 (2), sx, sc, v1x, v2x, v1c, v2c, fac, &
       v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw, v2cup, v2cdw, v2cud, &
       zeta, rh, rup, rdw, grhoup, grhodw, grhoud, grup, grdw, &
       sigma_gradcorr (3, 3)
  logical :: igcc_is_lyp
  !dummy variables for meta-gga
  real(DP) :: v3x,v3c,v3xup,v3xdw,v3cup,v3cdw

  if ( .not. dft_is_gradient() .and. .not. dft_is_meta() ) return
  !
  !if (noncolin) call errore('stres_gradcorr', &
  !                 'noncollinear stress + GGA not implemented',1)
  if ( dft_is_meta() .and. (nspin>1) )  call errore('stres_gradcorr', &
                    'Meta-GGA stress not yet implemented with spin polarization',1)



  igcc_is_lyp = (get_igcc() == 3)

  sigma_gradcorr(:,:) = 0.d0

  nr1 = dfft%nr1
  nr2 = dfft%nr2
  nr3 = dfft%nr3
  nrxx= dfft%nnr
  ngm = dfft%ngm
  !
  if(nspin==4) RETURN
  allocate (grho( 3, nrxx, nspin))    
  nspin0=nspin
  if (nspin==4) nspin0=1
  if (nspin==4.and.domag) nspin0=2
  fac = 1.d0 / DBLE (nspin0)
  !
  !    calculate the gradient of rho+rhocore in real space
  ALLOCATE( rhoout( dfftp%nnr, nspin0) )
  IF (nspin==4.AND.domag) THEN
!     ALLOCATE( vgg( dfftp%nnr, nspin0 ) )
!     ALLOCATE( vsave( dfftp%nnr, nspin ) )
     ALLOCATE( segni( dfftp%nnr ) )
!     vsave=v
!     v=0.d0
  ENDIF
  !
  ALLOCATE( rhogsum( ngm, nspin0 ) )
  !
  ! ... calculate the gradient of rho + rho_core in real space
  !
  IF ( nspin == 4 .AND. domag ) THEN
     !
     CALL compute_rho(rho,rhoout,segni,dfftp%nnr)
     !
     ! ... bring starting rhoout to G-space
     !
     DO is = 1, nspin0
        !
        psic(:) = rhoout(:,is)
        !
        CALL fwfft ('Rho', psic, dfftp)
        !
        rhogsum(:,is) = psic(dfftp%nl(:))
        !
     END DO
  ELSE
     !
     rhoout(:,1:nspin0)  = rho(:,1:nspin0)
     rhogsum(:,1:nspin0) = rhog(:,1:nspin0)
     !
  ENDIF

  !
  DO is = 1, nspin0
     !
     rhoout(:,is)  = fac * rho_core(:)  + rhoout(:,is)
     rhogsum(:,is) = fac * rhog_core(:) + rhogsum(:,is)
     !
     CALL fft_gradient_g2r( dfft, rhogsum(1,is), g, grho(1,1,is) )
     !
  END DO
  !
  DEALLOCATE(rhogsum)
  !
  if (nspin.eq.1) then
     !
     !    This is the LDA case
     !
     ! sigma_gradcor_{alpha,beta} ==
     !     omega^-1 \int (grad_alpha rho) ( D(rho*Exc)/D(grad_alpha rho) ) d3
     !
     do k = 1, nrxx
        !
        grho2 (1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
        !
        if (abs (rhoout (k, 1) ) .gt.epsr.and.grho2 (1) .gt.epsg) then
           !
           ! routine computing v1x and v2x is different for GGA and meta-GGA
           ! FIXME : inefficient implementation
           !
           if ( dft_is_meta() ) then
              !
              kedtau(k,1) = kedtau(k,1) / e2
              call tau_xc (rhoout(k,1), grho2(1),kedtau(k,1), sx, sc, v1x, v2x,v3x,v1c,v2c,v3c)
              kedtau(k,1) = kedtau(k,1) * e2
              !
           else
              !
              call gcxc (rhoout (k, 1), grho2(1), sx, sc, v1x, v2x, v1c, v2c)
              !
           endif
           !
           do l = 1, 3
              !
              do m = 1, l
                 !
                 sigma_gradcorr (l, m) = sigma_gradcorr (l, m) + &
                                grho(l,k,1) * grho(m,k,1) * e2 * (v2x + v2c)
                 !
              enddo
              !
           enddo
           !
        endif
        !
     enddo
     !
  else
     !
     !    This is the LSDA case
     !
     do k = 1, nrxx
        grho2 (1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
        grho2 (2) = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
        !
        if ( (abs (rhoout (k, 1) ) .gt.epsr.and.grho2 (1) .gt.epsg) .and. &
             (abs (rhoout (k, 2) ) .gt.epsr.and.grho2 (2) .gt.epsg) ) then
           !
           ! Spin polarization for metagga
           !
           if ( dft_is_meta() ) then
           !
             ! Not working with spin polarization ... FIXME
             ! call tau_xc_spin (rho(k,1), rho(k,2), grho2(1), grho2(2), &
             !            kedtau(k,1), kedtau(k,2), sx, sc, &
             !            v1xup, v1xdw, v2xup, v2xdw, v3xup, v3xdw, &
             !            v1cup, v1cdw, v2cup, v2cdw, v3cup, v3cdw ) 

           !
           else
              ! 
              !
              call gcx_spin (rhoout (k, 1), rhoout (k, 2), grho2 (1), grho2 (2), &
                   sx, v1xup, v1xdw, v2xup, v2xdw)
              !
              rh = rhoout (k, 1) + rhoout (k, 2)
              !
              if (rh.gt.epsr) then
                 if ( igcc_is_lyp ) then
                    rup = rhoout (k, 1)
                    rdw = rhoout (k, 2)
                    grhoup = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
                    grhodw = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
                    grhoud = grho(1,k,1) * grho(1,k,2) + &
                             grho(2,k,1) * grho(2,k,2) + &
                             grho(3,k,1) * grho(3,k,2)
                    call gcc_spin_more(rup, rdw, grhoup, grhodw, grhoud, sc, &
                                    v1cup, v1cdw, v2cup, v2cdw, v2cud)
                 else
                    zeta = (rhoout (k, 1) - rhoout (k, 2) ) / rh
              
                    grh2 = (grho (1, k, 1) + grho (1, k, 2) ) **2 + &
                           (grho (2, k, 1) + grho (2, k, 2) ) **2 + &
                           (grho (3, k, 1) + grho (3, k, 2) ) **2
                    call gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
                    v2cup = v2c
                    v2cdw = v2c
                    v2cud = v2c
                 end if
              else
                 sc = 0.d0
                 v1cup = 0.d0
                 v1cdw = 0.d0
                 v2c = 0.d0
                 v2cup = 0.d0
                 v2cdw = 0.d0
                 v2cud = 0.d0
              endif
              !
           endif
           !
           do l = 1, 3
              do m = 1, l
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
              enddo
              !
           enddo
           !
        endif
        !
     enddo
     !
  endif
  !
  do l = 1, 3
     do m = 1, l - 1
        sigma_gradcorr (m, l) = sigma_gradcorr (l, m)
     enddo

  enddo
  call mp_sum(  sigma_gradcorr, intra_bgrp_comm )

  sigmaxc(:,:) = sigmaxc(:,:) + sigma_gradcorr(:,:) / &
                                (nr1 * nr2 * nr3)
  
!  DO is = 1, nspin0
!     !
!     rho(:,is)  = rho(:,is)  - fac * rho_core(:)
!     rhog(:,is) = rhog(:,is) - fac * rhog_core(:)
!     !
!  END DO
  !
  deallocate(grho, rhoout)
  IF(allocated(segni)) DEALLOCATE(segni)
  return

end subroutine stres_gradcorr

