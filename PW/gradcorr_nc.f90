!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
subroutine gradcorr_nc (rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, &
     nrx3, nrxx, nl, ngm, g, alat, omega, e2, etxc, vtxc, v, nspin)
  !     ===================
  !--------------------------------------------------------------------
#include "f_defs.h"
  USE kinds
  use funct
  implicit none
  !
  integer :: nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, nl (ngm), &
       nspin

  real(kind=DP) :: rho(nrxx,nspin), &
     rho_core(nrxx),v(nrxx,nspin), &
     g(3,ngm),vtxc,etxc,e2,alat,omega,zeta,rh,grh2
  integer :: k, ipol, is, nspin0
  real(kind=DP), allocatable :: grho(:,:,:),h(:,:,:),dh(:),   &
                                segni(:),rhoout(:,:),vgg(:,:)
  real(kind=DP) :: grho2(2),sx,sc,v1x,v2x,v1c,v2c,v1xup,v1xdw,amag, &
       v2xup, v2xdw, v1cup, v1cdw , etxcgc, vtxcgc, segno, arho, fac
  real(kind=DP), parameter :: epsr = 1.0d-6, epsg = 1.0d-10

  if (igcx.eq.0.and.igcc.eq.0) return

  nspin0=nspin
  if (nspin0.eq.4) nspin0=2

  etxcgc = 0.d0
  vtxcgc = 0.d0

  allocate (h( 3, nrxx, nspin))    
  allocate (grho( 3, nrxx, nspin))    
  if (nspin == 4) allocate (segni(nrxx))
  allocate (rhoout(nrxx,nspin))
  allocate (vgg( nrxx, nspin))    

  vgg=0.d0
  if (nspin.eq.4) then
     call compute_rho(rho,rhoout,segni,nrxx)
  else
     call DCOPY(nspin*nrxx,rho,1,rhoout,1)
  endif 

  ! calculate the gradient of rho+rho_core in real space
  fac = 1.d0 / float (nspin0)
  do is = 1, nspin0
     call DAXPY (nrxx, fac, rho_core, 1, rhoout (1, is), 1)
     call gradient (nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, rhoout (1, is), &
          ngm, g, nl, alat, grho (1, 1, is) )
  enddo
  do k = 1, nrxx
     do is = 1, nspin0
        grho2 (is) = grho(1, k, is)**2 + grho(2, k, is)**2 + grho(3, k, is)**2
     enddo
     if (nspin.eq.1) then
        !
        !    This is the spin-unpolarised case
        !
        arho = abs (rhoout (k, 1) )
        segno = sign (1.d0, rhoout (k, 1) )
        if (arho.gt.epsr.and.grho2 (1) .gt.epsg) then

           call gcxc (arho, grho2, sx, sc, v1x, v2x, v1c, v2c)
           !
           ! first term of the gradient correction : D(rho*Exc)/D(rho)

           v (k, 1) = v (k, 1) + e2 * (v1x + v1c)
           ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
           do ipol = 1, 3
              h (ipol, k, 1) = e2 * (v2x + v2c) * grho (ipol, k, 1)
           enddo
           vtxcgc = vtxcgc + e2 * (v1x + v1c) * (rhoout(k, 1) - rho_core(k) )
           etxcgc = etxcgc + e2 * (sx + sc) * segno
        else
           do ipol = 1, 3
              h (ipol, k, 1) = 0.d0
           enddo
        endif
     else
        !
        !    spin-polarised case
        !
        call gcx_spin (rhoout (k, 1), rhoout (k, 2), grho2 (1), grho2 (2), &
             sx, v1xup, v1xdw, v2xup, v2xdw)
        rh = rhoout (k, 1) + rhoout (k, 2)
        if (rh.gt.epsr) then
           zeta = (rhoout (k, 1) - rhoout (k, 2) ) / rh
           if (nspin.eq.4) zeta=abs(zeta)*segni(k)
           grh2 = (grho (1, k, 1) + grho (1, k, 2) ) **2 + &
                  (grho (2, k, 1) + grho (2, k, 2) ) **2 + &
                  (grho (3, k, 1) + grho (3, k, 2) ) **2
           call gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
        else
           sc = 0.d0
           v1cup = 0.d0
           v1cdw = 0.d0
           v2c = 0.d0
        endif
        !
        ! first term of the gradient correction : D(rho*Exc)/D(rho)
        !
!!!!!!!!!
        vgg (k, 1) = vgg (k, 1) + e2 * (v1xup + v1cup)
        vgg (k, 2) = vgg (k, 2) + e2 * (v1xdw + v1cdw)
!!!!!!!!!
        !
        ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
        !
        do ipol = 1, 3
           h (ipol, k, 1) = e2 * ( (v2xup + v2c) * grho (ipol, k, 1) &
                + v2c * grho (ipol, k, 2) )
           h (ipol, k, 2) = e2 * ( (v2xdw + v2c) * grho (ipol, k, 2) &
                + v2c * grho (ipol, k, 1) )
        enddo
        vtxcgc = vtxcgc + e2 * (v1xup + v1cup) * (rhoout (k, 1) - &
             rho_core (k) * fac)
        vtxcgc = vtxcgc + e2 * (v1xdw + v1cdw) * (rhoout (k, 2) - &
             rho_core (k) * fac)
        etxcgc = etxcgc + e2 * (sx + sc)
     endif
  enddo
  do is = 1, nspin0
     call DAXPY (nrxx, - fac, rho_core, 1, rhoout (1, is), 1)
  enddo
  deallocate(grho)
  allocate (dh( nrxx))    
  !
  ! second term of the gradient correction :
  ! \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  do is = 1, nspin0
     call grad_dot (nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, h (1, 1, is), &
          ngm, g, nl, alat, dh)
     do k = 1, nrxx
!!!!!!!!
        vgg (k, is) = vgg (k, is) - dh (k)
!!!!!!!!
        vtxcgc = vtxcgc - dh (k) * rho (k, is)
     enddo
  enddo

  vtxc = vtxc + omega * vtxcgc / (nr1 * nr2 * nr3)
  etxc = etxc + omega * etxcgc / (nr1 * nr2 * nr3)

      if (nspin.eq.4) then
         do k=1,nrxx
            v(k,1)=v(k,1)+0.5d0*(vgg(k,1)+vgg(k,2))
            amag=sqrt(rho(k,2)**2+rho(k,3)**2+rho(k,4)**2)
            if (amag.gt.1.d-12) then
               v(k,2)=v(k,2)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,2)/amag
               v(k,3)=v(k,3)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,3)/amag
               v(k,4)=v(k,4)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,4)/amag
            endif
         enddo
      else
         call DAXPY(nspin*nrxx,1.d0,vgg,1,v,1)
      endif

  deallocate (dh)
  deallocate (vgg)    
  deallocate (h)
  if (nspin == 4) deallocate (segni)
  deallocate (rhoout)
  return

end subroutine gradcorr_nc
