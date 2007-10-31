!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine cg_setupdgc
  !-----------------------------------------------------------------------
  ! Setup all arrays needed in the gradient correction case
  ! This version requires on input allocated array
  !
  USE kinds, only: DP
  use pwcom
  use scf, only : rho, rhog, rho_core, rhog_core
  use cgcom
  use funct, only: gcxc, gcx_spin, gcc_spin, dft_is_gradient
  !
  implicit none
  integer k, is
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
  parameter (epsr=1.0d-6, epsg=1.0d-10)
  !
  if (.not. dft_is_gradient() ) return
  call start_clock('setup_dgc')
  !
  dvxc_rr(:,:,:) = 0.d0
  dvxc_sr(:,:,:) = 0.d0
  dvxc_ss(:,:,:) = 0.d0
  dvxc_s (:,:,:) = 0.d0
  grho (:,:,:) = 0.d0
  !
  !    add rho_core
  !
  fac=1.d0/DBLE(nspin)
  if (nlcc_any) then
     do is=1,nspin
        rho(:,is)  = fac * rho_core(:)  + rho(:,is)
        rhog(:,is) = fac * rhog_core(:) + rhog(:,is)
     enddo
  endif
  do is=1,nspin
     call gradrho (nrx1,nrx2,nrx3,nr1,nr2,nr3,nrxx,rhog(1,is),   &
          ngm,g,nl,grho(1,1,is))
  enddo
  !
  if (nspin.eq.1) then
     do k = 1,nrxx
        grho2(1)=grho(1,k,1)**2+grho(2,k,1)**2+grho(3,k,1)**2
        if (abs(rho(k,1)).gt.epsr.and.grho2(1).gt.epsg) then
           call gcxc(rho(k,nspin),grho2(1),sx,sc,v1x,v2x,v1c,v2c)
           call dgcxc(rho(k,nspin),grho2,vrrx,vsrx,vssx,vrrc,vsrc,vssc)
           dvxc_rr(k,1,1) = e2 * ( vrrx + vrrc )
           dvxc_sr(k,1,1) = e2 * ( vsrx + vsrc )
           dvxc_ss(k,1,1) = e2 * ( vssx + vssc )
           dvxc_s (k,1,1) = e2 * ( v2x + v2c )
        endif
     end do
  else
     do k = 1,nrxx
        grho2(2)=grho(1,k,2)**2+grho(2,k,2)**2+grho(3,k,2)**2
        rh=rho(k,1)+rho(k,2)
        grh2= (grho(1,k,1)+grho(1,k,2))**2                          &
                        + (grho(2,k,1)+grho(2,k,2))**2              &
                        + (grho(3,k,1)+grho(3,k,2))**2
        !
        call gcx_spin(rho(k,1),rho(k,2),grho2(1),grho2(2),sx,       &
             v1xup,v1xdw,v2xup,v2xdw)
        !
        call dgcxc_spin(rho(k,1),rho(k,2),grho(1,k,1),grho(1,k,2),     &
             vrrxup,vrrxdw,vrsxup,vrsxdw,vssxup,vssxdw, &
             vrrcup,vrrcdw,vrscup,vrscdw,vssc,vrzcup,vrzcdw)
        !
        if (rh.gt.epsr) then
           zeta=(rho(k,1)-rho(k,2))/rh
           call gcc_spin(rh,zeta,grh2,sc,v1cup,v1cdw,v2c)
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
        else
           dvxc_rr(k,1,1)=0.d0
           dvxc_rr(k,1,2)=0.d0
           dvxc_rr(k,2,1)=0.d0
           dvxc_rr(k,2,2)=0.d0
           !
           dvxc_s(k,1,1)=0.d0
           dvxc_s(k,1,2)=0.d0
           dvxc_s(k,2,1)=0.d0
           dvxc_s(k,2,2)=0.d0
        endif
        dvxc_sr(k,1,1)=e2*(vrsxup+vrscup)
        dvxc_sr(k,1,2)=e2*vrscup
        dvxc_sr(k,2,1)=e2*vrscdw
        dvxc_sr(k,2,2)=e2*(vrsxdw+vrscdw)
        !
        dvxc_ss(k,1,1)=e2*(vssxup+vssc)
        dvxc_ss(k,1,2)=e2*vssc
        dvxc_ss(k,2,1)=e2*vssc
        dvxc_ss(k,2,2)=e2*(vssxdw+vssc)
     enddo
  endif
  if (nlcc_any) then
     do is=1,nspin
        rho(:,is)  = rho(:,is)  - fac * rho_core(:)
        rhog(:,is) = rhog(:,is) - fac * rhog_core(:)
     enddo
  endif
  call stop_clock('setup_dgc')
  !
  return
end subroutine cg_setupdgc
