!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE gradcorr( rho, rhog, rho_core, rhog_core, etxc, vtxc, v )
  !----------------------------------------------------------------------------
  !! Calls the xc GGA drivers and calculates total energy and potential.
  !
  USE constants,            ONLY : e2
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : ngm, g
  USE lsda_mod,             ONLY : nspin
  USE cell_base,            ONLY : omega
  USE funct,                ONLY : igcc_is_lyp, dft_is_gradient, get_igcc
  USE xc_gga,               ONLY : xc_gcx, gcx_spin, gcc_spin
  USE spin_orb,             ONLY : domag
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE fft_rho,              ONLY: rho_r2g
  !
  IMPLICIT NONE
  !
  REAL(DP),    INTENT(IN)    :: rho(dfftp%nnr,nspin), rho_core(dfftp%nnr)
  COMPLEX(DP), INTENT(IN)    :: rhog(ngm,nspin), rhog_core(ngm)
  REAL(DP),    INTENT(INOUT) :: v(dfftp%nnr,nspin)
  REAL(DP),    INTENT(INOUT) :: vtxc, etxc
  !
  INTEGER :: k, ipol, is, nspin0, ir, jpol
  !
  REAL(DP), ALLOCATABLE :: grho(:,:,:), h(:,:,:), dh(:)
  REAL(DP), ALLOCATABLE :: rhoaux(:,:), segni(:), vgg(:,:), vsave(:,:)
  REAL(DP), ALLOCATABLE :: gmag(:,:,:)

  COMPLEX(DP), ALLOCATABLE :: rhogaux(:,:)
  !
  REAL(DP), ALLOCATABLE :: grho2(:,:), grho_ud(:)
  REAL(DP), ALLOCATABLE :: rh(:), zeta(:)  
  REAL(DP), ALLOCATABLE :: v1x(:,:), v2x(:,:)
  REAL(DP), ALLOCATABLE :: v1c(:,:), v2c(:,:), v2c_ud(:)
  REAL(DP) :: vnull
  REAL(DP) :: sx(dfftp%nnr), sc(dfftp%nnr)
  !
  REAL(DP) :: sgn(2), etxcgc, vtxcgc, segno, fac, amag 
  !
  REAL(DP) :: grup, grdw, seg, gr
  !
  REAL(DP), PARAMETER :: epsr = 1.D-6, epsg = 1.D-10
  !
  !
  IF ( .NOT. dft_is_gradient() ) RETURN
  !
  etxcgc = 0.0_DP
  vtxcgc = 0.0_DP
  !
  nspin0 = nspin
  IF ( nspin==4 ) nspin0 = 1
  IF ( nspin==4 .AND. domag ) nspin0 = 2
  fac = 1.0_DP / DBLE( nspin0 )
  sgn(1) = 1._DP ;  sgn(2) = -1._DP
  !
  ALLOCATE( h(3,dfftp%nnr,nspin0)    )
  ALLOCATE( grho(3,dfftp%nnr,nspin0) )
  ALLOCATE( grho2(dfftp%nnr,nspin0)  )
  ALLOCATE( rhoaux(dfftp%nnr,nspin0) )
  !
  ALLOCATE( v1x(dfftp%nnr,nspin0), v2x(dfftp%nnr,nspin0) )
  ALLOCATE( v1c(dfftp%nnr,nspin0), v2c(dfftp%nnr,nspin0) )
  !
  IF ( nspin==4 .AND. domag ) THEN
     ALLOCATE( vgg(dfftp%nnr,nspin0)  )
     ALLOCATE( vsave(dfftp%nnr,nspin) )
     ALLOCATE( segni(dfftp%nnr) )
     vsave=v
     v=0._DP
  ENDIF
  !
  ALLOCATE( rhogaux( ngm, nspin0 ) )
  !
  ! ... calculate the gradient of rho + rho_core in real space
  !
  IF ( nspin == 4 .AND. domag ) THEN
     !
     CALL compute_rho( rho, rhoaux, segni, dfftp%nnr ) 
     !
     ! ... bring starting rhoaux to G-space
     !
     CALL rho_r2g ( dfftp, rhoaux(:,1:nspin0), rhogaux(:,1:nspin0) )
     !
  ELSE
     !
     ! ... for convenience rhoaux and rhogaux are in (up,down) format, when LSDA
     !
     DO is = 1, nspin0
       rhoaux(:,is)  = (  rho(:,1) + sgn(is) *  rho(:,nspin0) ) * 0.5_DP
       rhogaux(:,is) = ( rhog(:,1) + sgn(is) * rhog(:,nspin0) ) * 0.5_DP
     ENDDO
     !
  ENDIF
  !
  DO is = 1, nspin0
     !
     rhoaux(:,is)  = fac *  rho_core(:) +  rhoaux(:,is)
     rhogaux(:,is) = fac * rhog_core(:) + rhogaux(:,is)
     !
     CALL fft_gradient_g2r( dfftp, rhogaux(1,is), g, grho(1,1,is) )
     !
  ENDDO
  !
  DEALLOCATE( rhogaux )
  !
  !
  IF ( nspin0 == 1 ) THEN
     !
     ! ... This is the spin-unpolarised case
     !
     CALL xc_gcx( dfftp%nnr, nspin0, rhoaux, grho, sx, sc, v1x, v2x, v1c, v2c )
     !
     DO k = 1, dfftp%nnr
        !
        ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
        v(k,1) = v(k,1) + e2 * ( v1x(k,1) + v1c(k,1) ) ! * vnull
        !
        ! ... h contains:  D(rho*Exc) / D(|grad rho|) * (grad rho) / |grad rho|
        h(:,k,1) = e2 * ( v2x(k,1) + v2c(k,1) ) * grho(:,k,1)
        !
        vtxcgc = vtxcgc + e2 * ( v1x(k,1) + v1c(k,1) ) * &
                               (rhoaux(k,1) - rho_core(k) )
        !
        etxcgc = etxcgc + e2 * ( sx(k) + sc(k) )
        !
     ENDDO
     !
  ELSE
     !
     ! ... spin-polarised case
     !
     ALLOCATE( v2c_ud(dfftp%nnr) )
     !   
     !
     CALL xc_gcx( dfftp%nnr, nspin0, rhoaux, grho, sx, sc, v1x, v2x, v1c, v2c, v2c_ud )
     !
     !
     ! ... h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
     !
     DO k = 1, dfftp%nnr
        v(k,1:nspin0) = v(k,1:nspin0) + e2*( v1x(k,1:nspin0) + v1c(k,1:nspin0)) 
        DO ipol = 1, 3
           !
           grup = grho(ipol,k,1)
           grdw = grho(ipol,k,2)
           h(ipol,k,1) = e2*( (v2x(k,1) + v2c(k,1))*grup + v2c_ud(k)*grdw )
           h(ipol,k,2) = e2*( (v2x(k,2) + v2c(k,2))*grdw + v2c_ud(k)*grup )
           !
        ENDDO
        !
        vtxcgc = vtxcgc + &
                e2 * ( v1x(k,1) + v1c(k,1) ) * ( rhoaux(k,1) - rho_core(k) * fac )
        vtxcgc = vtxcgc + &
                e2 * ( v1x(k,2) + v1c(k,2) ) * ( rhoaux(k,2) - rho_core(k) * fac )
        etxcgc = etxcgc + e2 * ( sx(k) + sc(k) )
        !
     ENDDO
     !
     DEALLOCATE( v2c_ud )
     !
  ENDIF
  !
  !
  DO is = 1, nspin0
     !
     rhoaux(:,is) = rhoaux(:,is) - fac * rho_core(:)
     !
  END DO
  !
  DEALLOCATE( grho )
  DEALLOCATE( v1x, v2x )
  DEALLOCATE( v1c, v2c )
  !
  ALLOCATE( dh(dfftp%nnr) )    
  !
  ! ... second term of the gradient correction :
  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  DO is = 1, nspin0
     !
     CALL fft_graddot( dfftp, h(1,1,is), g, dh )
     !
     v(:,is) = v(:,is) - dh(:)
     !
     vtxcgc = vtxcgc - SUM( dh(:) * rhoaux(:,is) )
     !
  END DO
  !
  vtxc = vtxc + omega * vtxcgc / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  etxc = etxc + omega * etxcgc / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  !
  IF (nspin==4 .AND. domag) THEN
     DO is = 1, nspin0
        vgg(:,is) = v(:,is)
     ENDDO
     !
     v = vsave
     DO k = 1, dfftp%nnr
        v(k,1) = v(k,1) + 0.5d0*(vgg(k,1)+vgg(k,2))
        amag = SQRT(rho(k,2)**2+rho(k,3)**2+rho(k,4)**2)
        IF (amag > 1.d-12) THEN
           v(k,2) = v(k,2) + segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,2)/amag
           v(k,3) = v(k,3) + segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,3)/amag
           v(k,4) = v(k,4) + segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,4)/amag
        ENDIF
     ENDDO
  ENDIF
  !
  DEALLOCATE( dh )
  DEALLOCATE( h )
  DEALLOCATE( grho2 )
  DEALLOCATE( rhoaux )
  IF (nspin==4 .AND. domag) THEN
     DEALLOCATE( vgg )
     DEALLOCATE( vsave )
     DEALLOCATE( segni )
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE gradcorr
