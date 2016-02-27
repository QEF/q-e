!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine setup_dgc
  !-----------------------------------------------------------------------
  !
  ! Allocate and setup all variable needed in the gradient correction case
  !
  ! GGA+LSDA is allowed. ADC (September 1999).
  ! GGA+LSDA+NLCC is allowed. ADC (November 1999).
  ! GGA+noncollinear+NLCC is allowed. ADC (June 2007).
  !
  USE constants,            ONLY : e2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, g, nl
  USE spin_orb,             ONLY : domag
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE noncollin_module,     ONLY : noncolin, ux, nspin_gga, nspin_mag
  USE wavefunctions_module, ONLY : psic
  USE kinds,                ONLY : DP
  USE funct,                ONLY : dft_is_gradient, gcxc, gcx_spin, &
                                   gcc_spin, dgcxc, dgcxc_spin
  USE uspp,                 ONLY : nlcc_any
  USE gc_lr,                ONLY : grho, gmag, dvxc_rr, dvxc_sr, &
                                   dvxc_ss, dvxc_s, vsgga, segni

  implicit none
  integer :: k, is, ipol, jpol, ir
  real(DP) :: grho2 (2), rh, zeta, grh2, fac, sx, sc, &
       v1x, v2x, v1c, v2c, vrrx, vsrx, vssx, vrrc, vsrc, vssc, v1xup, &
       v1xdw, v2xup, v2xdw, v1cup, v1cdw, vrrxup, vrrxdw, vrsxup, vrsxdw, &
       vssxup, vssxdw, vrrcup, vrrcdw, vrscup, vrscdw, vrzcup, vrzcdw,   &
       amag, seg, seg0
  COMPLEX(DP), ALLOCATABLE :: rhogout(:,:)
  real(DP), allocatable :: rhoout(:,:)
  real (DP), parameter :: epsr = 1.0d-6, epsg = 1.0d-10

  IF ( .NOT. dft_is_gradient() ) RETURN
  
  IF (noncolin.AND.domag) THEN
     allocate (segni (dfftp%nnr))
     allocate (vsgga (dfftp%nnr))
     allocate (gmag (3, dfftp%nnr, nspin_mag))
     gmag=0.0_dp
  ENDIF

  IF(.NOT.ALLOCATED(dvxc_rr)) ALLOCATE (dvxc_rr(dfftp%nnr, nspin_gga , nspin_gga))
  IF(.NOT.ALLOCATED(dvxc_sr)) ALLOCATE (dvxc_sr(dfftp%nnr, nspin_gga , nspin_gga))
  IF(.NOT.ALLOCATED(dvxc_ss)) ALLOCATE (dvxc_ss(dfftp%nnr, nspin_gga , nspin_gga))
  IF(.NOT.ALLOCATED(dvxc_s))  ALLOCATE (dvxc_s (dfftp%nnr, nspin_gga , nspin_gga))
  IF(.NOT.ALLOCATED(grho))    ALLOCATE (grho   (  3    , dfftp%nnr, nspin_gga))
  IF(.NOT.ALLOCATED(rhoout))  ALLOCATE (rhoout ( dfftp%nnr, nspin_gga))

  dvxc_rr(:,:,:) = 0.d0
  dvxc_sr(:,:,:) = 0.d0
  dvxc_ss(:,:,:) = 0.d0
  dvxc_s (:,:,:) = 0.d0
  grho   (:,:,:) = 0.d0
  !
  !    add rho_core
  !
  fac = 1.d0 / DBLE (nspin_gga)
  IF (noncolin.and.domag) THEN
     allocate(rhogout(ngm,nspin_mag))
     call compute_rho(rho%of_r,rhoout,segni,dfftp%nnr)
     DO is = 1, nspin_gga
        !
        if (nlcc_any) rhoout(:,is)  = fac * rho_core(:)  + rhoout(:,is)

        psic(:) = rhoout(:,is)
        !
        CALL fwfft ('Dense', psic, dfftp)
        !
        rhogout(:,is) = psic(nl(:))
        !
        !
        CALL gradrho(dfftp%nnr, rhogout(1,is), ngm, g, nl, grho(1,1,is) )
        !
     END DO
     DEALLOCATE(rhogout)
  ELSE
     do is = 1, nspin_gga
        rhoout(:,is)  =  rho%of_r(:,is)
     enddo
     if (nlcc_any) then
        do is = 1, nspin_gga
           rhoout(:,is)  = fac * rho_core(:)  + rho%of_r(:,is)
           rho%of_g(:,is) = fac * rhog_core(:) + rho%of_g(:,is)
        enddo
     endif
     do is = 1, nspin_gga
        call gradrho (dfftp%nnr, rho%of_g (1, is), ngm, g, nl, grho (1, 1, is) )
     enddo
  END IF


  do k = 1, dfftp%nnr
     grho2 (1) = grho (1, k, 1) **2 + grho (2, k, 1) **2 + grho (3, k, 1) **2
     if (nspin_gga == 1) then
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
        do is = 1, nspin_gga
           rho%of_g(:,is) = rho%of_g(:,is) - fac * rhog_core(:)
        enddo
     endif
  endif

  DEALLOCATE(rhoout)

  RETURN

end subroutine setup_dgc
