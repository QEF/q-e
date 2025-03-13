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
  USE xc_lib,               ONLY : igcc_is_lyp, xclib_dft_is, xc_gcx
  USE noncollin_module,     ONLY : domag
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE fft_rho,              ONLY : rho_r2g
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
  COMPLEX(DP), ALLOCATABLE :: rhogaux(:,:)
  !
  REAL(DP), ALLOCATABLE :: v1x(:,:), v2x(:,:)
  REAL(DP), ALLOCATABLE :: v1c(:,:), v2c(:,:), v2c_ud(:)
  REAL(DP), ALLOCATABLE :: sx(:), sc(:)
  !
  REAL(DP) :: sgn_is, etxcgc, vtxcgc, fac, amag
  REAL(DP) :: grup, grdw
  REAL(DP), PARAMETER :: epsr = 1.D-6, epsg = 1.D-10
  !
  IF ( .NOT. xclib_dft_is('gradient') ) RETURN
  !
  !$acc data present( rho, rho_core, rhog, rhog_core, v )
  !
  etxcgc = 0.0_DP
  vtxcgc = 0.0_DP
  !
  nspin0 = nspin
  IF ( nspin==4 ) nspin0 = 1
  IF ( nspin==4 .AND. domag ) nspin0 = 2
  fac = 1.0_DP / DBLE( nspin0 )
  !
  ALLOCATE( h(3,dfftp%nnr,nspin0)    )
  ALLOCATE( grho(3,dfftp%nnr,nspin0) )
  ALLOCATE( rhoaux(dfftp%nnr,nspin0) )
  ALLOCATE( v1x(dfftp%nnr,nspin0), v2x(dfftp%nnr,nspin0) )
  ALLOCATE( v1c(dfftp%nnr,nspin0), v2c(dfftp%nnr,nspin0) )
  ALLOCATE( sx(dfftp%nnr), sc(dfftp%nnr) )
  !$acc data create( rhoaux, grho, sx, sc, v1x, v2x, v1c, v2c, h )
  !
  ALLOCATE( rhogaux(ngm,nspin0) )
  !$acc data create( rhogaux )
  !
  ! ... calculate the gradient of rho + rho_core in real space
  !
  IF ( nspin == 4 .AND. domag ) THEN
    ALLOCATE( vsave(dfftp%nnr,nspin) )
    ALLOCATE( segni(dfftp%nnr) )
    !
    !$acc data copyout( vsave )
    !$acc parallel loop collapse(2)
    DO is = 1, nspin
      DO ir = 1, dfftp%nnr
        vsave(ir,is) = v(ir,is)
        v(ir,is) = 0._DP
      ENDDO
    ENDDO
    !$acc end data
    !
    ! ... bring starting rhoaux to G-space
    CALL compute_rho( rho, rhoaux, segni, dfftp%nnr )
    CALL rho_r2g( dfftp, rhoaux(:,1:nspin0), rhogaux(:,1:nspin0) )
    !
  ELSE
    ! ... for convenience rhoaux and rhogaux are in (up,down) format, when LSDA
    !
    !$acc parallel loop collapse(2)
    DO is = 1, nspin0
      DO ir = 1, dfftp%nnr
        sgn_is = DBLE(3-2*is)
        rhoaux(ir,is) = ( rho(ir,1) + sgn_is * rho(ir,nspin0) ) * 0.5_DP
      ENDDO
    ENDDO
    !$acc parallel loop collapse(2)
    DO is = 1, nspin0
      DO ir = 1, ngm
        sgn_is = DBLE(3-2*is)
        rhogaux(ir,is) = ( rhog(ir,1) + CMPLX(sgn_is,KIND=DP) * rhog(ir,nspin0) ) &
                         * (0.5_DP,0._DP)
      ENDDO
    ENDDO
    !
  ENDIF
  !
  !$acc parallel loop collapse(2)
  DO is = 1, nspin0
    DO ir = 1, dfftp%nnr
      rhoaux(ir,is) = fac * rho_core(ir) + rhoaux(ir,is)
    ENDDO
  ENDDO
  !$acc parallel loop collapse(2)
  DO is = 1, nspin0
    DO ir = 1, ngm
      rhogaux(ir,is) = CMPLX(fac,kind=DP) * rhog_core(ir) + rhogaux(ir,is)
    ENDDO
  ENDDO
  !
  DO is = 1, nspin0
    CALL fft_gradient_g2r( dfftp, rhogaux(:,is), g, grho(:,:,is) )
  ENDDO
  !
  !$acc end data
  DEALLOCATE( rhogaux )
  !
  IF ( nspin0 == 1 ) THEN
     !
     ! ... This is the spin-unpolarised case
     !
     CALL xc_gcx( dfftp%nnr, nspin0, rhoaux, grho, sx, sc, v1x, v2x, v1c, v2c, &
                  gpu_args_=.TRUE. )
     !
     !$acc parallel loop reduction(+:etxcgc) reduction(+:vtxcgc)
     DO k = 1, dfftp%nnr
        ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
        v(k,1) = v(k,1) + e2 * ( v1x(k,1) + v1c(k,1) )
        ! ... h contains:  D(rho*Exc) / D(|grad rho|) * (grad rho) / |grad rho|
        DO ipol = 1, 3
          h(ipol,k,1) = e2 * ( v2x(k,1) + v2c(k,1) ) * grho(ipol,k,1)
        ENDDO
        !
        vtxcgc = vtxcgc + e2 * ( v1x(k,1) + v1c(k,1) ) * &
                               ( rhoaux(k,1) - rho_core(k) )
        etxcgc = etxcgc + e2 * ( sx(k) + sc(k) )
     ENDDO
     !
  ELSE
     !
     ! ... spin-polarised case
     !
     ALLOCATE( v2c_ud(dfftp%nnr) )
     !$acc data create( v2c_ud )
     !
     CALL xc_gcx( dfftp%nnr, nspin0, rhoaux, grho, sx, sc, v1x, v2x, v1c, v2c, &
                  v2c_ud, gpu_args_=.TRUE. )
     !
     ! ... h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
     !
     !$acc parallel loop reduction(+:etxcgc) reduction(+:vtxcgc)
     DO k = 1, dfftp%nnr
        !
        DO is = 1, nspin0
           v(k,is) = v(k,is) + e2*(v1x(k,is) + v1c(k,is))
        ENDDO
        !
        DO ipol = 1, 3
           grup = grho(ipol,k,1)
           grdw = grho(ipol,k,2)
           h(ipol,k,1) = e2*( (v2x(k,1) + v2c(k,1))*grup + v2c_ud(k)*grdw )
           h(ipol,k,2) = e2*( (v2x(k,2) + v2c(k,2))*grdw + v2c_ud(k)*grup )
        ENDDO
        !
        vtxcgc = vtxcgc + e2 * ( &
                 ( v1x(k,1) + v1c(k,1) ) * ( rhoaux(k,1) - rho_core(k) * fac ) + &
                 ( v1x(k,2) + v1c(k,2) ) * ( rhoaux(k,2) - rho_core(k) * fac ) )
        etxcgc = etxcgc + e2 * ( sx(k) + sc(k) )
        !
     ENDDO
     !
     !$acc end data
     DEALLOCATE( v2c_ud )
     !
  ENDIF
  !
  !$acc parallel loop collapse(2)
  DO is = 1, nspin0
    DO k = 1, dfftp%nnr
      rhoaux(k,is) = rhoaux(k,is) - fac * rho_core(k)
    ENDDO
  ENDDO
  !
  ALLOCATE( dh(dfftp%nnr) )
  !$acc data create( dh )
  !
  ! ... second term of the gradient correction :
  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  DO is = 1, nspin0
     CALL fft_graddot( dfftp, h(1,1,is), g, dh )
     !$acc parallel loop reduction(-:vtxcgc)
     DO k = 1, dfftp%nnr
       v(k,is) = v(k,is) - dh(k)
       vtxcgc = vtxcgc - dh(k) * rhoaux(k,is)
     ENDDO
  ENDDO
  !
  !$acc end data
  !
  vtxc = vtxc + omega * vtxcgc / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  etxc = etxc + omega * etxcgc / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  !
  IF (nspin==4 .AND. domag) THEN
     ALLOCATE( vgg(dfftp%nnr,nspin0)  )
     !$acc data copyin( segni, vsave ) create( vgg )
     !
     !$acc parallel loop collapse(2)
     DO is = 1, nspin
       DO ir = 1, dfftp%nnr
         IF (is<=nspin0) vgg(ir,is) = v(ir,is)
         v(ir,is) = vsave(ir,is)
       ENDDO
     ENDDO
     !
     !$acc parallel loop
     DO k = 1, dfftp%nnr
        v(k,1) = v(k,1) + 0.5d0*(vgg(k,1)+vgg(k,2))
        amag = SQRT(rho(k,2)**2+rho(k,3)**2+rho(k,4)**2)
        IF (amag > 1.d-12) THEN
           v(k,2) = v(k,2) + segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,2)/amag
           v(k,3) = v(k,3) + segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,3)/amag
           v(k,4) = v(k,4) + segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,4)/amag
        ENDIF
     ENDDO
     !
     !$acc end data
     DEALLOCATE( segni )
     DEALLOCATE( vgg, vsave )
  ENDIF
  !
  !$acc end data
  !
  DEALLOCATE( sc, sx )
  DEALLOCATE( rhoaux, grho )
  DEALLOCATE( v1x, v2x )
  DEALLOCATE( v1c, v2c )
  DEALLOCATE( dh, h )
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE gradcorr
