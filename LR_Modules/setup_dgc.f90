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
  !! Allocate and set up all variable needed in the gradient correction case.
  !
  !! GGA+LSDA is allowed. ADC (September 1999);  
  !! GGA+LSDA+NLCC is allowed. ADC (November 1999);  
  !! GGA+noncollinear+NLCC is allowed. ADC (June 2007).
  !
  USE constants,            ONLY : e2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, g
  USE spin_orb,             ONLY : domag
  USE scf,                  ONLY : rho, rho_core, rhog_core, rhoz_or_updw
  USE noncollin_module,     ONLY : noncolin, ux, nspin_gga, nspin_mag
  USE wavefunctions,        ONLY : psic
  USE kinds,                ONLY : DP
  USE funct,                ONLY : dft_is_gradient
  USE xc_gga,               ONLY : gcxc, gcx_spin, gcc_spin, libxc_switches_gga
  USE uspp,                 ONLY : nlcc_any
  USE gc_lr,                ONLY : grho, gmag, dvxc_rr, dvxc_sr, &
                                   dvxc_ss, dvxc_s, vsgga, segni
  !
  IMPLICIT NONE
  !
  INTEGER :: k, is, ipol, jpol, ir
  !
  REAL(DP), ALLOCATABLE :: rh(:), zeta(:)
  REAL(DP), ALLOCATABLE :: grh(:,:,:), grho2(:,:), grh2(:)
  REAL(DP) :: ne2, fac, sgn(2), null_v(dfftp%nnr)
  !
  REAL(DP), ALLOCATABLE :: v1x(:,:), v2x(:,:), v1c(:,:), v2c(:,:)
  REAL(DP), ALLOCATABLE :: vrrx(:,:), vsrx(:,:), vssx(:,:), vrrc(:,:)
  REAL(DP), ALLOCATABLE :: vsrc(:,:), vrzc(:,:), vssc(:), sc(:), sx(:)
  !
  REAL(DP), ALLOCATABLE :: rhoout(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogout(:,:)
  !
  REAL(DP), PARAMETER :: epsr=1.0d-6, epsg=1.0d-10
  REAL(DP), PARAMETER :: rho_trash=0.5d0, zeta_trash=0.2d0, &
                         grho2_trash=0.1d0
  !
  IF ( .NOT. dft_is_gradient() ) RETURN
  !
  CALL start_clock( 'setup_dgc' )
  ! 
  IF ( SUM(libxc_switches_gga(:)) /= 0 ) CALL errore( 'setup_dgc', 'libxc derivatives of &
                                                      &xc potentials for GGA not implemented yet', 1 )
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
  IF (.NOT.ALLOCATED(grho2)  ) ALLOCATE( grho2(dfftp%nnr,nspin_gga)  )
  IF (.NOT.ALLOCATED(rhoout) ) ALLOCATE( rhoout(dfftp%nnr,nspin_gga) )
  !
  ALLOCATE( rh(dfftp%nnr) )
  ALLOCATE( v1x(dfftp%nnr,nspin_gga), v2x(dfftp%nnr,nspin_gga) )
  ALLOCATE( v1c(dfftp%nnr,nspin_gga), v2c(dfftp%nnr,nspin_gga) )
  ALLOCATE( vrrx(dfftp%nnr,nspin_gga), vsrx(dfftp%nnr,nspin_gga) )
  ALLOCATE( vssx(dfftp%nnr,nspin_gga), vrrc(dfftp%nnr,nspin_gga) )
  ALLOCATE( vsrc(dfftp%nnr,nspin_gga), vssc(dfftp%nnr) )
  ALLOCATE( sx(dfftp%nnr), sc(dfftp%nnr) )
  IF (nspin_gga /= 1) THEN
     ALLOCATE( zeta(dfftp%nnr) )
     ALLOCATE( grh(dfftp%nnr,3,nspin_gga), grh2(dfftp%nnr) )
     ALLOCATE( vrzc(dfftp%nnr,nspin_gga) )
  ENDIF
  !
  dvxc_rr(:,:,:) = 0.d0
  dvxc_sr(:,:,:) = 0.d0
  dvxc_ss(:,:,:) = 0.d0
  dvxc_s(:,:,:) = 0.d0
  grho(:,:,:) = 0.d0
  !
  sgn(1)=1.d0  ;   sgn(2)=-1.d0
  fac = 1.d0/DBLE(nspin_gga)
  !
  IF (noncolin .AND. domag) THEN
     !
     ALLOCATE( rhogout(ngm,nspin_mag) )
     !
     CALL compute_rho( rho%of_r, rhoout, segni, dfftp%nnr )
     !
     DO is = 1, nspin_gga
        IF (nlcc_any) rhoout(:,is) = fac*rho_core(:) + rhoout(:,is)
        psic(:) = rhoout(:,is)
        CALL fwfft( 'Rho', psic, dfftp )
        rhogout(:,is) = psic(dfftp%nl(:))
        CALL fft_gradient_g2r( dfftp, rhogout(1,is), g, grho(1,1,is) )
     ENDDO
     !
     DEALLOCATE( rhogout )
     !
  ELSE
     ! ... for convenience, if LSDA, rhoout is kept in (up,down) format
     DO is = 1, nspin_gga
        rhoout(:,is) = ( rho%of_r(:,1) + sgn(is)*rho%of_r(:,nspin_gga) )*0.5d0
     ENDDO
     !
     ! ... if LSDA rho%of_g is temporarily converted in (up,down) format
     CALL rhoz_or_updw( rho, 'only_g', '->updw' )
     !
     IF (nlcc_any) THEN
        DO is = 1, nspin_gga
           rhoout(:,is) = fac * rho_core(:)  + rhoout(:,is)
           rho%of_g(:,is) = fac * rhog_core(:) + rho%of_g(:,is)
        ENDDO
     ENDIF
     !
     DO is = 1, nspin_gga
        CALL fft_gradient_g2r( dfftp, rho%of_g(1,is), g, grho(1,1,is) )
     ENDDO
     !
  ENDIF
  !
  grho2(:,1) = grho(1,:,1)**2 + grho(2,:,1)**2 + grho(3,:,1)**2
  !
  null_v = 1.0_DP
  !
  IF (nspin_gga == 1) THEN
     !
     rh(:) = rhoout(:,1)
     WHERE (ABS(rhoout(:,1))<=epsr .OR. grho2(:,1)<=epsg)
        rh = rho_trash
        grho2(:,1) = grho2_trash
        null_v = 0.0_DP
     END WHERE
     !
     CALL gcxc( dfftp%nnr, rh, grho2(:,1), sx, sc, v1x(:,1), &
                                v2x(:,1), v1c(:,1), v2c(:,1) )
     !
     CALL dgcxc( dfftp%nnr, rh, grho2(:,1), vrrx(:,1), vsrx(:,1), &
                       vssx(:,1), vrrc(:,1), vsrc(:,1), vssc )
     !
     dvxc_rr(:,1,1) = e2 * (vrrx(:,1) + vrrc(:,1)) * null_v
     dvxc_sr(:,1,1) = e2 * (vsrx(:,1) + vsrc(:,1)) * null_v
     dvxc_ss(:,1,1) = e2 * (vssx(:,1) + vssc(:)  ) * null_v
     dvxc_s(:,1,1)  = e2 * (v2x(:,1)  + v2c(:,1) ) * null_v
     !
  ELSE
     !
     grho2(:,2) = grho(1,:,2)**2 + grho(2,:,2)**2 + grho(3,:,2)**2
     !   
     CALL gcx_spin( dfftp%nnr, rhoout, grho2, sx, v1x, v2x )
     !
     ! ... swap grho indices to match dgcx_spin input (waiting for a better fix)
     DO k = 1, dfftp%nnr
        grh(k,1:3,1) = grho(1:3,k,1)
        grh(k,1:3,2) = grho(1:3,k,2)
     ENDDO
     !
     CALL dgcxc_spin( dfftp%nnr, rhoout, grh, vrrx, vsrx, vssx, vrrc, &
                      vsrc, vssc, vrzc )
     !
     rh = rhoout(:,1) + rhoout(:,2)
     grh2(:) = (grho(1,:,1)+grho(1,:,2))**2 + (grho(2,:,1) &   
               +grho(2,:,2))**2 + (grho(3,:,1)+grho(3,:,2))**2
     !
     WHERE (rh(:) > epsr)
        zeta = (rhoout(:,1)-rhoout(:,2)) / rh
     ELSEWHERE
        rh = rho_trash
        zeta = zeta_trash
        null_v = 0.0_DP
     END WHERE
     !
     CALL gcc_spin( dfftp%nnr, rh, zeta, grh2, sc, v1c, v2c(:,1) )
     !
     DO k = 1, dfftp%nnr
        ne2 = null_v(k) * e2
        !
        dvxc_rr(k,1,1) = ne2 * (vrrx(k,1) + vrrc(k,1) + vrzc(k,1) * (1.d0 - zeta(k)) / rh(k))
        dvxc_rr(k,1,2) = ne2 * (vrrc(k,1) - vrzc(k,1) * (1.d0 + zeta(k)) / rh(k))
        dvxc_rr(k,2,1) = ne2 * (vrrc(k,2) + vrzc(k,2) * (1.d0 - zeta(k)) / rh(k))
        dvxc_rr(k,2,2) = ne2 * (vrrx(k,2) + vrrc(k,2) - vrzc(k,2) * (1.d0 + zeta(k)) / rh(k))
        !
        dvxc_s(k,1,1) = ne2 * (v2x(k,1) + v2c(k,1))
        dvxc_s(k,1,2) = ne2 * v2c(k,1)
        dvxc_s(k,2,1) = ne2 * v2c(k,1)
        dvxc_s(k,2,2) = ne2 * (v2x(k,2) + v2c(k,1))
        !
        dvxc_sr(k,1,1) = e2 * (vsrx(k,1) + vsrc(k,1))   
        dvxc_sr(k,1,2) = e2 * vsrc(k,1)   
        dvxc_sr(k,2,1) = e2 * vsrc(k,2)   
        dvxc_sr(k,2,2) = e2 * (vsrx(k,2) + vsrc(k,2))   
        !
        dvxc_ss(k,1,1) = e2 * (vssx(k,1) + vssc(k))   
        dvxc_ss(k,1,2) = e2 * vssc(k)   
        dvxc_ss(k,2,1) = e2 * vssc(k)   
        dvxc_ss(k,2,2) = e2 * (vssx(k,2) + vssc(k))
     ENDDO
     !
  ENDIF
  !
  !
  IF (noncolin .AND. domag) THEN
     CALL compute_vsgga( rhoout, grho, vsgga )
  ELSE
     IF (nlcc_any) THEN
        DO is = 1, nspin_gga
           rho%of_g(:,is) = rho%of_g(:,is) - fac*rhog_core(:)
        ENDDO
     ENDIF
     !
     CALL rhoz_or_updw( rho, 'only_g', '->rhoz' )
     !
  ENDIF
  !
  DEALLOCATE( rh, grho2  )
  DEALLOCATE( v1x, v2x, v1c, v2c )
  DEALLOCATE( vrrx, vsrx, vssx, vrrc )
  DEALLOCATE( vsrc, vssc, sx, sc )
  IF (nspin_gga /= 1) THEN
     DEALLOCATE( zeta, grh, grh2 )
     DEALLOCATE( vrzc )
  ENDIF
  DEALLOCATE( rhoout )
  !
  CALL stop_clock( 'setup_dgc' )
  !
  RETURN
  !
END SUBROUTINE setup_dgc
