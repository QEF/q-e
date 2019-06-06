!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE cg_setupdgc
  !-----------------------------------------------------------------------
  !! Set up all arrays needed in the gradient correction case.  
  !! This version requires allocated arrays on input.
  !
  USE kinds,     ONLY: DP
  USE constants, ONLY: e2
  USE scf,       ONLY: rho, rho_core, rhog_core, rhoz_or_updw
  USE funct,     ONLY: dft_is_gradient, init_xc
  USE xc_gga,    ONLY: gcxc, gcx_spin, gcc_spin, libxc_switches_gga
  USE fft_base,  ONLY: dfftp
  USE gvect,     ONLY: ngm, g
  USE lsda_mod,  ONLY: nspin
  USE uspp,      ONLY: nlcc_any
  USE cgcom
  !
  IMPLICIT NONE
  !
  INTEGER k, is
  !
  REAL(DP), ALLOCATABLE :: rh(:), zeta(:)
  REAL(DP), ALLOCATABLE :: grh(:,:,:), grho2(:,:), grh2(:)
  REAL(DP) :: ne2, fac, sgn(2), null_v(dfftp%nnr)
  !
  REAL(DP), ALLOCATABLE :: v1x(:,:), v2x(:,:), v1c(:,:), v2c(:,:)
  REAL(DP), ALLOCATABLE :: vrrx(:,:), vsrx(:,:), vssx(:,:), vrrc(:,:)
  REAL(DP), ALLOCATABLE :: vsrc(:,:), vrzc(:,:), vssc(:), sc(:), sx(:)
  !
  REAL(DP), PARAMETER :: epsr=1.0d-6, epsg=1.0d-10
  REAL(DP), PARAMETER :: rho_trash=0.5d0, zeta_trash=0.2d0, &
                         grho2_trash=0.1d0
  !
  IF (.NOT. dft_is_gradient() ) RETURN
  !
  CALL start_clock( 'setup_dgc' )
  !
  dvxc_rr(:,:,:) = 0.d0
  dvxc_sr(:,:,:) = 0.d0
  dvxc_ss(:,:,:) = 0.d0
  dvxc_s(:,:,:) = 0.d0
  grho(:,:,:) = 0.d0
  !
  CALL init_xc( 'GGA' )
  !
  IF ( SUM(libxc_switches_gga(:)) /= 0 ) CALL errore( 'cg_setupdgc', 'libxc derivatives of &
                                                      &xc potentials for GGA not implemented yet', 1 )
  !
  ALLOCATE( rh(dfftp%nnr), grho2(dfftp%nnr,nspin) )
  ALLOCATE( v1x(dfftp%nnr,nspin), v2x(dfftp%nnr,nspin) )
  ALLOCATE( v1c(dfftp%nnr,nspin), v2c(dfftp%nnr,nspin) )
  ALLOCATE( vrrx(dfftp%nnr,nspin), vsrx(dfftp%nnr,nspin) )
  ALLOCATE( vssx(dfftp%nnr,nspin), vrrc(dfftp%nnr,nspin) )
  ALLOCATE( vsrc(dfftp%nnr,nspin), vssc(dfftp%nnr) )
  ALLOCATE( sx(dfftp%nnr), sc(dfftp%nnr) )
  IF (nspin /= 1) THEN
     ALLOCATE( zeta(dfftp%nnr) )
     ALLOCATE( grh(dfftp%nnr,3,nspin), grh2(dfftp%nnr) )
     ALLOCATE( vrzc(dfftp%nnr,nspin) )
  ENDIF
  !
  ! ... add rho_core to the charge density rho
  !
  IF (nlcc_any) THEN
     rho%of_r(:,1) = rho_core(:)  + rho%of_r(:,1)
     rho%of_g(:,1) = rhog_core(:) + rho%of_g(:,1)
  ENDIF
  !
  ! ... for LSDA, convert rho to (up,down) (gradient grho is (up,down))
  !
  IF (nspin == 2) CALL rhoz_or_updw( rho, 'r_and_g', '->updw' )
  !
  DO is = 1, nspin
     CALL fft_gradient_g2r( dfftp, rho%of_g(1,is), g, grho(1,1,is) )
  ENDDO
  !
  grho2(:,1) = grho(1,:,1)**2 + grho(2,:,1)**2 + grho(3,:,1)**2
  !
  IF (nspin == 1) THEN
     !
     rh(:) = rho%of_r(:,1)
     null_v(:) = 1.0_DP
     WHERE (ABS(rho%of_r(:,1))<=epsr .OR. grho2(:,1)<=epsg)
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
     CALL gcx_spin( dfftp%nnr, rho%of_r, grho2, sx, v1x, v2x )
     !
     ! ... swap grho indices to match dgcx_spin input (waiting for a better fix)
     DO k = 1, dfftp%nnr
        grh(k,1:3,1) = grho(1:3,k,1)
        grh(k,1:3,2) = grho(1:3,k,2)
     ENDDO
     !
     CALL dgcxc_spin( dfftp%nnr, rho%of_r, grh, vrrx, vsrx, vssx, vrrc, &
                      vsrc, vssc, vrzc )
     !
     rh = rho%of_r(:,1) + rho%of_r(:,2)
     grh2(:) = (grho(1,:,1)+grho(1,:,2))**2 + (grho(2,:,1) &   
               +grho(2,:,2))**2 + (grho(3,:,1)+grho(3,:,2))**2
     !
     null_v(:) = 1.0_DP
     WHERE (rh > epsr)
        zeta = (rho%of_r(:,1)-rho%of_r(:,2)) / rh
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
  ! ... restore rho to its input value
  IF (nspin == 2) CALL rhoz_or_updw( rho, 'r_and_g', '->rhoz' )
  !
  IF (nlcc_any) THEN
     rho%of_r(:,1) = rho%of_r(:,1) - rho_core(:)
     rho%of_g(:,1) = rho%of_g(:,1) - rhog_core(:)
  ENDIF
  !
  DEALLOCATE( rh, grho2 )
  DEALLOCATE( v1x, v2x, v1c, v2c )
  DEALLOCATE( vrrx, vsrx, vssx, vrrc )
  DEALLOCATE( vsrc, vssc, sx, sc )
  IF (nspin /= 1) THEN
     DEALLOCATE( zeta, grh, grh2 )
     DEALLOCATE( vrzc )
  ENDIF
  !
  CALL stop_clock( 'setup_dgc' )
  !
  RETURN
  !
END SUBROUTINE cg_setupdgc
