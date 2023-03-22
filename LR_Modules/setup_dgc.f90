!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE setup_dgc
  !-----------------------------------------------------------------------
  !! This subroutine computes \(\text{dvxc}\), the derivative of the XC 
  !! potential for the gradient correction (GGA).
  !
  !! GGA+LSDA is allowed. ADC (September 1999);  
  !! GGA+LSDA+NLCC is allowed. ADC (November 1999);  
  !! GGA+noncollinear+NLCC is allowed. ADC (June 2007).
  !
  USE constants,            ONLY : e2
  USE fft_base,             ONLY : dfftp
  USE fft_rho,              ONLY : rho_r2g
  USE gvect,                ONLY : ngm, g
  USE scf,                  ONLY : rho, rho_core, rhog_core, rhoz_or_updw
  USE noncollin_module,     ONLY : noncolin, domag, ux, nspin_gga, nspin_mag
  USE kinds,                ONLY : DP
  USE xc_lib,               ONLY : xclib_dft_is, xc_gcx, dgcxc
  USE uspp,                 ONLY : nlcc_any
  USE gc_lr,                ONLY : grho, gmag, dvxc_rr, dvxc_sr, &
                                   dvxc_ss, dvxc_s, vsgga, segni
  !
  IMPLICIT NONE
  !
  INTEGER :: k, is, ipol, jpol, ir, dfftp_nnr
  !
  REAL(DP), ALLOCATABLE :: grh(:,:,:)
  REAL(DP) :: fac, sgn(2)
  !
  REAL(DP), ALLOCATABLE :: v1x(:,:), v2x(:,:), v1c(:,:), v2c(:,:)
  REAL(DP), ALLOCATABLE :: v2c_ud(:)
  REAL(DP), ALLOCATABLE :: sc(:), sx(:)
  !
  REAL(DP), ALLOCATABLE :: rhoout(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogout(:,:)
  !
  REAL(DP), PARAMETER :: epsr=1.0d-6, epsg=1.0d-10
  !
  !
  IF ( .NOT. xclib_dft_is('gradient') ) RETURN
  !
  CALL start_clock( 'setup_dgc' )
  !
  dfftp_nnr = dfftp%nnr !to avoid unnecessary copies in acc loop
  !
  !$acc data copyin( rho )
  !$acc data copyin( rho_core, rhog_core, rho%of_r, rho%of_g )
  !
  IF (noncolin .AND. domag) THEN
     ALLOCATE( segni(dfftp%nnr) )
     ALLOCATE( vsgga(dfftp%nnr) )
     ALLOCATE( gmag(3,dfftp%nnr,nspin_mag) )
     gmag = 0.0_dp
  ENDIF
  !
  IF (.NOT.ALLOCATED(dvxc_rr)) ALLOCATE( dvxc_rr(dfftp%nnr,nspin_gga,nspin_gga) )
  IF (.NOT.ALLOCATED(dvxc_sr)) ALLOCATE( dvxc_sr(dfftp%nnr,nspin_gga,nspin_gga) )
  IF (.NOT.ALLOCATED(dvxc_ss)) ALLOCATE( dvxc_ss(dfftp%nnr,nspin_gga,nspin_gga) )
  IF (.NOT.ALLOCATED(dvxc_s) ) ALLOCATE( dvxc_s(dfftp%nnr,nspin_gga,nspin_gga)  )
  IF (.NOT.ALLOCATED(grho)   ) ALLOCATE( grho(3,dfftp%nnr,nspin_gga) )
  IF (.NOT.ALLOCATED(rhoout) ) ALLOCATE( rhoout(dfftp%nnr,nspin_gga) )
  !
  ALLOCATE( v1x(dfftp%nnr,nspin_gga), v2x(dfftp%nnr,nspin_gga) )
  ALLOCATE( v1c(dfftp%nnr,nspin_gga), v2c(dfftp%nnr,nspin_gga) )
  IF (nspin_gga == 2) ALLOCATE( v2c_ud(dfftp%nnr) )
  ALLOCATE( sx(dfftp%nnr), sc(dfftp%nnr) )
  ALLOCATE( grh(dfftp%nnr,3,nspin_gga) )
  !
  !$acc data create( rhoout, grh, sx, sc, v1x, v2x, v1c, v2c, v2c_ud )
  !$acc data copyout( dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, grho )
  !
  fac = 1.d0/DBLE(nspin_gga)
  !
  !$acc kernels
  dvxc_rr(:,:,:) = 0.d0
  dvxc_sr(:,:,:) = 0.d0
  dvxc_ss(:,:,:) = 0.d0
  dvxc_s(:,:,:)  = 0.d0
  grho(:,:,:) = 0.d0
  !$acc end kernels
  !
  IF (noncolin .AND. domag) THEN
     !
     ALLOCATE( rhogout(ngm,nspin_mag) )
     !$acc data create(rhogout)
     !
     CALL compute_rho( rho%of_r, rhoout, segni, dfftp%nnr )
     !
     DO is = 1, nspin_gga
        !
        IF (nlcc_any) THEN
          !$acc kernels
          rhoout(:,is) = rhoout(:,is) + fac*rho_core(:)
          !$acc end kernels
        ENDIF
        !
        CALL rho_r2g( dfftp, rhoout(:,is), rhogout(:,is:is) )
        CALL fft_gradient_g2r( dfftp, rhogout(1,is), g, grho(1,1,is) )
        !
     ENDDO
     !
     !$acc end data
     DEALLOCATE( rhogout )
     !
  ELSE
     ! ... for convenience, if LSDA, rhoout is kept in (up,down) format
     !
     sgn(1)=1.d0  ;   sgn(2)=-1.d0
     !
     !$acc parallel loop collapse(2) copyin(sgn) present(rho)
     DO is = 1, nspin_gga
       DO ir = 1, dfftp_nnr
         rhoout(ir,is) = ( rho%of_r(ir,1) + sgn(is)*rho%of_r(ir,nspin_gga) )*0.5d0
       ENDDO
     ENDDO
     !
     ! ... if LSDA rho%of_g is temporarily converted in (up,down) format
     !
     CALL rhoz_or_updw( rho, 'only_g', '->updw' )
     !
     IF (nlcc_any) THEN
        DO is = 1, nspin_gga
           !$acc kernels present(rho)
           rhoout(:,is) = fac * rho_core(:) + rhoout(:,is)
           rho%of_g(:,is) = fac * rhog_core(:) + rho%of_g(:,is)
           !$acc end kernels
        ENDDO
     ENDIF
     !
     DO is = 1, nspin_gga
        CALL fft_gradient_g2r( dfftp, rho%of_g(1,is), g, grho(1,1,is) )
     ENDDO
     !
  ENDIF
  !
  ! ... swap grho indices to match xc_gcx input (waiting for a better fix)
  !$acc parallel loop
  DO k = 1, dfftp_nnr
     grh(k,1:3,1) = grho(1:3,k,1)
     IF (nspin_gga==2) grh(k,1:3,2) = grho(1:3,k,2)
  ENDDO
  !
  !
  IF (nspin_gga == 1) THEN
     !
     CALL dgcxc( dfftp_nnr, 1, rhoout, grh, dvxc_rr, dvxc_sr, dvxc_ss, &
                 gpu_args_=.TRUE. )
     !
     !$acc parallel loop
     DO ir = 1, dfftp_nnr
       IF( rhoout(ir,1)<0.d0 ) rhoout(ir,1)=0.d0
     ENDDO
     !
     CALL xc_gcx( dfftp_nnr, nspin_gga, rhoout, grho, sx, sc, v1x, v2x, v1c, v2c, &
                  gpu_args_=.TRUE. )
     !
     !$acc kernels
     dvxc_s(:,1,1)  = e2 * (v2x(:,1) + v2c(:,1))
     !$acc end kernels
     !
  ELSE
     !
     CALL dgcxc( dfftp_nnr, nspin_gga, rhoout, grh, dvxc_rr, dvxc_sr, dvxc_ss, &
                 gpu_args_=.TRUE. )
     !
     CALL xc_gcx( dfftp_nnr, nspin_gga, rhoout, grho, sx, sc, v1x, v2x, v1c, v2c, &
                  v2c_ud, gpu_args_=.TRUE. )
     !
     !$acc parallel loop
     DO k = 1, dfftp_nnr
        IF ( rhoout(k,1)+rhoout(k,2) > epsr) THEN
           dvxc_s(k,1,1) = e2 * (v2x(k,1) + v2c(k,1))
           dvxc_s(k,1,2) = e2 * v2c(k,1)
           dvxc_s(k,2,1) = e2 * v2c(k,1)
           dvxc_s(k,2,2) = e2 * (v2x(k,2) + v2c(k,1))
        ENDIF
     ENDDO
     !
  ENDIF
  !
  IF (noncolin .AND. domag) THEN
     CALL compute_vsgga( rhoout, grho, vsgga )
  ELSE
     IF (nlcc_any) THEN
        DO is = 1, nspin_gga
           !$acc kernels present(rho)
           rho%of_g(:,is) = rho%of_g(:,is) - fac*rhog_core(:)
           !$acc end kernels
        ENDDO
     ENDIF
     !
     CALL rhoz_or_updw( rho, 'only_g', '->rhoz' )
     !
  ENDIF
  !
  !$acc end data
  !$acc end data
  !
  !$acc end data
  !$acc end data
  !
  DEALLOCATE( v1x, v2x, v1c, v2c )
  IF (nspin_gga == 2) DEALLOCATE( v2c_ud )
  DEALLOCATE( sx, sc )
  DEALLOCATE( grh )
  DEALLOCATE( rhoout )
  !
  CALL stop_clock( 'setup_dgc' )
  !
  RETURN
  !
END SUBROUTINE setup_dgc
