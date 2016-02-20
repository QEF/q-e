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
  ! Setup all arrays needed in the gradient correction case
  ! This version requires on input allocated array
  !
  USE kinds, ONLY: dp
  USE constants, ONLY: e2
  USE scf,   ONLY: rho, rho_core, rhog_core
  USE funct, ONLY: gcxc, gcx_spin, gcc_spin, dgcxc, dgcxc_spin, dft_is_gradient
  USE fft_base, ONLY: dfftp
  USE gvect,    ONLY: nl, ngm, g
  USE lsda_mod, ONLY: nspin
  USE uspp,     ONLY: nlcc_any
  USE cgcom
  !
  IMPLICIT NONE
  INTEGER k, is
  real(DP) &
       &       grho2(2), rh, zeta, grh2, epsr, epsg, fac,                 &
       &       sx,sc,v1x,v2x,v1c,v2c,vrrx,vsrx,vssx,                      &
       &       vrrc,vsrc,vssc,                                            &
       &       v1xup,v1xdw,v2xup,v2xdw,                                   &
       &       v1cup,v1cdw,                                               &
       &       vrrxup,vrrxdw,vrsxup,vrsxdw,vssxup,vssxdw,                 &
       &       vrrcup,vrrcdw,vrscup,vrscdw,                               &
       &       vrzcup,vrzcdw
  !
  PARAMETER (epsr=1.0d-6, epsg=1.0d-10)
  !
  IF (.not. dft_is_gradient() ) RETURN
  CALL start_clock('setup_dgc')
  !
  dvxc_rr(:,:,:) = 0.d0
  dvxc_sr(:,:,:) = 0.d0
  dvxc_ss(:,:,:) = 0.d0
  dvxc_s (:,:,:) = 0.d0
  grho (:,:,:) = 0.d0
  !
  !    add rho_core
  !
  fac=1.d0/dble(nspin)
  IF (nlcc_any) THEN
     DO is=1,nspin
        rho%of_r(:,is)  = fac * rho_core(:)  + rho%of_r(:,is)
        rho%of_g(:,is) = fac * rhog_core(:) + rho%of_g(:,is)
     ENDDO
  ENDIF
  DO is=1,nspin
     CALL gradrho (dfftp%nnr, rho%of_g(1,is), ngm, g, nl, grho(1,1,is))
  ENDDO
  !
  IF (nspin==1) THEN
     DO k = 1,dfftp%nnr
        grho2(1)=grho(1,k,1)**2+grho(2,k,1)**2+grho(3,k,1)**2
        IF (abs(rho%of_r(k,1))>epsr.and.grho2(1)>epsg) THEN
           CALL gcxc(rho%of_r(k,nspin),grho2(1),sx,sc,v1x,v2x,v1c,v2c)
           CALL dgcxc(rho%of_r(k,nspin),grho2(1),vrrx,vsrx,vssx,vrrc,vsrc,vssc)
           dvxc_rr(k,1,1) = e2 * ( vrrx + vrrc )
           dvxc_sr(k,1,1) = e2 * ( vsrx + vsrc )
           dvxc_ss(k,1,1) = e2 * ( vssx + vssc )
           dvxc_s (k,1,1) = e2 * ( v2x + v2c )
        ENDIF
     ENDDO
  ELSE
     DO k = 1,dfftp%nnr
        grho2(2)=grho(1,k,2)**2+grho(2,k,2)**2+grho(3,k,2)**2
        rh=rho%of_r(k,1)+rho%of_r(k,2)
        grh2= (grho(1,k,1)+grho(1,k,2))**2                          &
                        + (grho(2,k,1)+grho(2,k,2))**2              &
                        + (grho(3,k,1)+grho(3,k,2))**2
        !
        CALL gcx_spin(rho%of_r(k,1),rho%of_r(k,2),grho2(1),grho2(2),sx,       &
             v1xup,v1xdw,v2xup,v2xdw)
        !
        CALL dgcxc_spin(rho%of_r(k,1),rho%of_r(k,2),grho(1,k,1),grho(1,k,2),     &
             vrrxup,vrrxdw,vrsxup,vrsxdw,vssxup,vssxdw, &
             vrrcup,vrrcdw,vrscup,vrscdw,vssc,vrzcup,vrzcdw)
        !
        IF (rh>epsr) THEN
           zeta=(rho%of_r(k,1)-rho%of_r(k,2))/rh
           CALL gcc_spin(rh,zeta,grh2,sc,v1cup,v1cdw,v2c)
           !
           dvxc_rr(k,1,1)=e2*(vrrxup+vrrcup+vrzcup*(1.d0-zeta)/rh)
           dvxc_rr(k,1,2)=e2*(vrrcup-vrzcup*(1.d0+zeta)/rh)
           dvxc_rr(k,2,1)=e2*(vrrcdw+vrzcdw*(1.d0-zeta)/rh)
           dvxc_rr(k,2,2)=e2*(vrrxdw+vrrcdw-vrzcdw*(1.d0+zeta)/rh)
           !
           dvxc_s(k,1,1)=e2*(v2xup+v2c)
           dvxc_s(k,1,2)=e2*v2c
           dvxc_s(k,2,1)=e2*v2c
           dvxc_s(k,2,2)=e2*(v2xdw+v2c)
        ELSE
           dvxc_rr(k,1,1)=0.d0
           dvxc_rr(k,1,2)=0.d0
           dvxc_rr(k,2,1)=0.d0
           dvxc_rr(k,2,2)=0.d0
           !
           dvxc_s(k,1,1)=0.d0
           dvxc_s(k,1,2)=0.d0
           dvxc_s(k,2,1)=0.d0
           dvxc_s(k,2,2)=0.d0
        ENDIF
        dvxc_sr(k,1,1)=e2*(vrsxup+vrscup)
        dvxc_sr(k,1,2)=e2*vrscup
        dvxc_sr(k,2,1)=e2*vrscdw
        dvxc_sr(k,2,2)=e2*(vrsxdw+vrscdw)
        !
        dvxc_ss(k,1,1)=e2*(vssxup+vssc)
        dvxc_ss(k,1,2)=e2*vssc
        dvxc_ss(k,2,1)=e2*vssc
        dvxc_ss(k,2,2)=e2*(vssxdw+vssc)
     ENDDO
  ENDIF
  IF (nlcc_any) THEN
     DO is=1,nspin
        rho%of_r(:,is)  = rho%of_r(:,is)  - fac * rho_core(:)
        rho%of_g(:,is) = rho%of_g(:,is) - fac * rhog_core(:)
     ENDDO
  ENDIF
  CALL stop_clock('setup_dgc')
  !
  RETURN
END SUBROUTINE cg_setupdgc
