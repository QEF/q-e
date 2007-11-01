
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define __OLD_NONCOLIN_GGA
!-----------------------------------------------------------------------
subroutine setup_dgc
  !-----------------------------------------------------------------------
  ! Allocate and setup all variable needed in the gradient correction case
  !
  !  GGA+LSDA is allowed. ADC (September 1999).
  !  GGA+LSDA+NLCC is allowed. ADC (November 1999).
  !  GGA+noncollinear+NLCC is allowed. ADC (June 2007).
  !
#include "f_defs.h"

  use pwcom
  use scf, only : rho, rhog, rho_core, rhog_core
  USE noncollin_module, ONLY : noncolin, ux
  USE wavefunctions_module, ONLY : psic
  USE kinds, only : DP
  use phcom
  use funct, only : dft_is_gradient, gcxc, gcx_spin, gcc_spin
  implicit none
  integer :: k, is, nspin0, ipol, jpol, ir
  real(DP) :: grho2 (2), rh, zeta, grh2, fac, sx, sc, &
       v1x, v2x, v1c, v2c, vrrx, vsrx, vssx, vrrc, vsrc, vssc, v1xup, &
       v1xdw, v2xup, v2xdw, v1cup, v1cdw, vrrxup, vrrxdw, vrsxup, vrsxdw, &
       vssxup, vssxdw, vrrcup, vrrcdw, vrscup, vrscdw, vrzcup, vrzcdw,   &
       amag, seg, seg0
  COMPLEX(DP), ALLOCATABLE :: rhogout(:,:)
  real(DP), allocatable :: rhoout(:,:)
  real (DP), parameter :: epsr = 1.0d-6, epsg = 1.0d-10

  if ( .not. dft_is_gradient() ) return

  nspin0=nspin
  IF (noncolin) THEN
     IF (domag) THEN
        allocate (segni (nrxx))    
        allocate (vsgga (nrxx))    
        allocate (gmag (3, nrxx, nspin))    
        nspin0=2
        gmag=0.0_dp
     ELSE
        nspin0=1
     ENDIF
  ENDIF

  allocate (dvxc_rr(  nrxx , nspin0 , nspin0))    
  allocate (dvxc_sr(  nrxx , nspin0 , nspin0))    
  allocate (dvxc_ss(  nrxx , nspin0 , nspin0))    
  allocate (dvxc_s (  nrxx , nspin0 , nspin0))    
  allocate (grho   (  3    , nrxx   , nspin0))    
  allocate (rhoout (  nrxx , nspin0))    

  dvxc_rr(:,:,:) = 0.d0
  dvxc_sr(:,:,:) = 0.d0
  dvxc_ss(:,:,:) = 0.d0
  dvxc_s (:,:,:) = 0.d0
  grho   (:,:,:) = 0.d0
  !
  !    add rho_core
  !
  fac = 1.d0 / DBLE (nspin0)
  IF (noncolin.and.domag) THEN
     allocate(rhogout(ngm,nspin))
#ifdef __OLD_NONCOLIN_GGA
     call compute_rho(rho%of_r,rhoout,segni,nrxx)
     DO is = 1, nspin0
        !
        if (nlcc_any) rhoout(:,is)  = fac * rho_core(:)  + rhoout(:,is)
       
        psic(:) = rhoout(:,is)
        !
        CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
        !
        rhogout(:,is) = psic(nl(:))
        !
        !
        CALL gradrho( nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, &
                   rhogout(1,is), ngm, g, nl, grho(1,1,is) )
        !
     END DO
#else
     call compute_rho_new(rho%of_r,rhoout,segni,nrxx,ux)
     do is=1,nspin
        rhogout(:,is) = rhog(:,is)
     enddo
     if (nlcc_any) then
        rhogout(:,1) = rhog_core(:) + rhog(:,1)
     endif
     do is = 1, nspin
        call gradrho (nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, rhogout(1, is), &
             ngm, g, nl, gmag (1, 1, is) )
     enddo
     DO is=1,nspin0
        IF (is==1) seg0=0.5_dp
        IF (is==2) seg0=-0.5_dp
        DO ipol=1,3
           grho(ipol,:,is) = 0.5_dp*gmag(ipol,:,1)
           DO ir=1,nrxx
              seg=seg0*segni(ir)
              rhoout(ir,is) = fac*rho_core(ir) + 0.5_dp*rho%of_r(ir,1)
              amag=sqrt(rho%of_r(ir,2)**2+rho%of_r(ir,3)**2+rho%of_r(ir,4)**2)
              IF (amag>1.d-12) THEN
                 rhoout(ir,is)=rhoout(ir,is)+seg*amag
                 DO jpol=2,4
                    grho(ipol,ir,is)=grho(ipol,ir,is)+ seg*rho%of_r(ir,jpol)* &
                                                 gmag(ipol,ir,jpol)/amag
                 END DO
              END IF
           END DO
        END DO
     END DO
!     write(6,*) 'setup dgc'
!     do k=2,2
!        write(6,'(3f20.5)') gmag(3,k,1), grho(3,k,1)
!     enddo
#endif
     DEALLOCATE(rhogout)
  ELSE
     do is = 1, nspin0
        rhoout(:,is)  =  rho%of_r(:,is)
     enddo
     if (nlcc_any) then
        do is = 1, nspin0
           rhoout(:,is)  = fac * rho_core(:)  + rho%of_r(:,is)
           rhog(:,is) = fac * rhog_core(:) + rhog(:,is)
        enddo
     endif
     do is = 1, nspin0
        call gradrho (nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, rhog (1, is), &
             ngm, g, nl, grho (1, 1, is) )
     enddo
  END IF


  do k = 1, nrxx
     grho2 (1) = grho (1, k, 1) **2 + grho (2, k, 1) **2 + grho (3, k, 1) **2
     if (nspin0 == 1) then
        if (abs (rhoout (k, 1) ) > epsr .and. grho2 (1) > epsg) then
           call gcxc (rhoout (k, 1), grho2(1), sx, sc, v1x, v2x, v1c, v2c)
           call dgcxc (rhoout (k, 1), grho2(1), vrrx, vsrx, vssx, vrrc, &
                vsrc, vssc)
           dvxc_rr (k, 1, 1) = e2 * (vrrx + vrrc)
           dvxc_sr (k, 1, 1) = e2 * (vsrx + vsrc)
           dvxc_ss (k, 1, 1) = e2 * (vssx + vssc)
           dvxc_s (k, 1, 1) = e2 * (v2x + v2c)
        endif
     else
        grho2 (2) = grho(1, k, 2) **2 + grho(2, k, 2) **2 + grho(3, k, 2) **2
        rh = rhoout (k, 1) + rhoout (k, 2)

        grh2 = (grho (1, k, 1) + grho (1, k, 2) ) **2 + (grho (2, k, 1) &
             + grho (2, k, 2) ) **2 + (grho (3, k, 1) + grho (3, k, 2) ) ** 2

        call gcx_spin (rhoout (k, 1), rhoout (k, 2), grho2 (1), grho2 (2), &
             sx, v1xup, v1xdw, v2xup, v2xdw)

        call dgcxc_spin (rhoout (k, 1), rhoout (k, 2), grho (1, k, 1), &
             grho (1, k, 2), vrrxup, vrrxdw, vrsxup, vrsxdw, vssxup, vssxdw, &
             vrrcup, vrrcdw, vrscup, vrscdw, vssc, vrzcup, vrzcdw)
        if (rh > epsr) then
           zeta = (rhoout (k, 1) - rhoout (k, 2) ) / rh
           call gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
           dvxc_rr (k, 1, 1) = e2 * (vrrxup + vrrcup + vrzcup * &
                (1.d0 - zeta) / rh)
           dvxc_rr (k, 1, 2) = e2 * (vrrcup - vrzcup * (1.d0 + zeta) / rh)
           dvxc_rr (k, 2, 1) = e2 * (vrrcdw + vrzcdw * (1.d0 - zeta) / rh)
           dvxc_rr (k, 2, 2) = e2 * (vrrxdw + vrrcdw - vrzcdw * &
                (1.d0 + zeta) / rh)
           dvxc_s (k, 1, 1) = e2 * (v2xup + v2c)
           dvxc_s (k, 1, 2) = e2 * v2c
           dvxc_s (k, 2, 1) = e2 * v2c
           dvxc_s (k, 2, 2) = e2 * (v2xdw + v2c)
        else
           dvxc_rr (k, 1, 1) = 0.d0
           dvxc_rr (k, 1, 2) = 0.d0
           dvxc_rr (k, 2, 1) = 0.d0
           dvxc_rr (k, 2, 2) = 0.d0
           dvxc_s (k, 1, 1) = 0.d0
           dvxc_s (k, 1, 2) = 0.d0
           dvxc_s (k, 2, 1) = 0.d0
           dvxc_s (k, 2, 2) = 0.d0
        endif
        dvxc_sr (k, 1, 1) = e2 * (vrsxup + vrscup)
        dvxc_sr (k, 1, 2) = e2 * vrscup
        dvxc_sr (k, 2, 1) = e2 * vrscdw
        dvxc_sr (k, 2, 2) = e2 * (vrsxdw + vrscdw)
        dvxc_ss (k, 1, 1) = e2 * (vssxup + vssc)
        dvxc_ss (k, 1, 2) = e2 * vssc
        dvxc_ss (k, 2, 1) = e2 * vssc
        dvxc_ss (k, 2, 2) = e2 * (vssxdw + vssc)
     endif
  enddo
  if (noncolin.and.domag) then
     call compute_vsgga(rhoout, grho, vsgga)
  else
     if (nlcc_any) then
        do is = 1, nspin0
           rhog(:,is) = rhog(:,is) - fac * rhog_core(:)
        enddo
     endif
  endif

!  write(6,*) 'setup dgc'
!  do k=2,2
!     write(6,'(3f20.5)') rhoout(k,1), grho(3,k,1)
!     write(6,'(3f20.5)') rhoout(k,2), grho(3,k,2)
!  enddo
!  write(6,*) 'exit setup dgc'

  DEALLOCATE(rhoout)

  return
end subroutine setup_dgc
