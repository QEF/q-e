!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
subroutine dgradcorr (rho, grho, dvxc_rr, dvxc_sr, dvxc_ss, &
     dvxc_s, xq, drho, nrxx, nspin, nspin0, nl, ngm, g, alat, dvxc)
  !--------------------------------------------------------------------
  !
  !  Add gradient correction contribution to 
  !  the responce exchange-correlation potential dvxc.
  !  LSDA is allowed.         ADC (September 1999)
  !  Noncollinear is allowed. ADC (June 2007)
  !
  USE kinds,            ONLY : DP
  USE noncollin_module, ONLY : noncolin
  USE spin_orb,         ONLY : domag
  USE gc_lr,            ONLY : gmag, vsgga, segni

  implicit none
  integer :: nrxx, ngm, nl (ngm), &
       nspin, nspin0
  real(DP) :: rho (nrxx, nspin), grho (3, nrxx, nspin0), &
       dvxc_rr(nrxx, nspin0, nspin0), dvxc_sr (nrxx, nspin0, nspin0), &
       dvxc_ss (nrxx,nspin0, nspin0), dvxc_s (nrxx, nspin0, nspin0),&
       g (3, ngm), xq(3), alat
  complex(DP) :: drho (nrxx, nspin), dvxc (nrxx, nspin)

  real(DP), parameter :: epsr = 1.0d-6, epsg = 1.0d-10
  real(DP) :: grho2, seg, seg0, amag
  complex(DP) :: s1, fact, term
  complex(DP) :: a (2, 2, 2), b (2, 2, 2, 2), c (2, 2, 2), &
                      ps (2, 2), ps1 (3, 2, 2), ps2 (3, 2, 2, 2)
  complex(DP), allocatable  :: gdrho (:,:,:), h (:,:,:), dh (:)
  complex(DP), allocatable  :: gdmag (:,:,:), dvxcsave(:,:), vgg(:,:)
  complex(DP), allocatable  :: drhoout(:,:)
  real(DP), allocatable :: rhoout(:,:)
  integer :: k, ipol, jpol, is, js, ks, ls

  if (noncolin.and.domag) then
     allocate (gdmag(3, nrxx, nspin))
     allocate (dvxcsave(nrxx, nspin))
     allocate (vgg(nrxx, nspin0))
     dvxcsave=dvxc
     dvxc=(0.0_dp,0.0_dp)
  endif
  allocate (rhoout( nrxx, nspin0))
  allocate (drhoout( nrxx, nspin0))
  allocate (gdrho( 3, nrxx, nspin0))
  allocate (h( 3, nrxx, nspin0))
  allocate (dh( nrxx))

  h (:, :, :) = (0.d0, 0.d0)
  if (noncolin.and.domag) then
     do is = 1, nspin
        call qgradient (xq, nrxx, &
            drho (1, is), ngm, g, nl, alat, gdmag (1, 1, is) )
     enddo
     DO is=1,nspin0
        IF (is==1) seg0=0.5_dp
        IF (is==2) seg0=-0.5_dp
        rhoout(:,is) = 0.5_dp*rho(:,1)
        drhoout(:,is) = 0.5_dp*drho(:,1)
        DO ipol=1,3
           gdrho(ipol,:,is) = 0.5_dp*gdmag(ipol,:,1)
        ENDDO
        DO k=1,nrxx
           seg=seg0*segni(k)
           amag=sqrt(rho(k,2)**2+rho(k,3)**2+rho(k,4)**2)
           IF (amag>1.d-12) THEN
              rhoout(k,is) = rhoout(k,is)+seg*amag
              DO jpol=2,4
                 drhoout(k,is) = drhoout(k,is)+seg*rho(k,jpol)* &
                                                 drho(k,jpol)/amag
              END DO
              DO ipol=1,3
                 fact=(0.0_dp,0.0_dp)
                 DO jpol=2,4
                    fact=fact+rho(k,jpol)*drho(k,jpol)
                 END DO
                 DO jpol=2,4
                    gdrho(ipol,k,is) = gdrho(ipol,k,is)+ seg*( &
                        drho(k,jpol)*gmag(ipol,k,jpol)+ &
                        rho(k,jpol)*gdmag(ipol,k,jpol))/amag &
                        -seg*(rho(k,jpol)*gmag(ipol,k,jpol)*fact)/amag**3
                 END DO
              END DO
           END IF
        END DO
     END DO
  ELSE
     DO is = 1, nspin0
        CALL qgradient (xq, nrxx, &
            drho (1, is), ngm, g, nl, alat, gdrho (1, 1, is) )
        rhoout(:,is)=rho(:,is)
        drhoout(:,is)=drho(:,is)
     ENDDO
  ENDIF

  do k = 1, nrxx
     grho2 = grho(1, k, 1)**2 + grho(2, k, 1)**2 + grho(3, k, 1)**2
     if (nspin == 1) then
        !
        !    LDA case
        !
        if (abs (rho (k, 1) ) > epsr .and. grho2 > epsg) then
           s1 = grho (1, k, 1) * gdrho (1, k, 1) + &
                grho (2, k, 1) * gdrho (2, k, 1) + &
                grho (3, k, 1) * gdrho (3, k, 1)
           !
           ! linear variation of the first term
           !
           dvxc (k, 1) = dvxc (k, 1) + dvxc_rr (k, 1, 1) * drho (k, 1) &
                + dvxc_sr (k, 1, 1) * s1
           do ipol = 1, 3
              h (ipol, k, 1) = (dvxc_sr(k, 1, 1) * drho(k, 1) + &
                                dvxc_ss(k, 1, 1) * s1 )*grho(ipol, k, 1) + &
                                dvxc_s (k, 1, 1) * gdrho (ipol, k, 1)
           enddo
        else
           do ipol = 1, 3
              h (ipol, k, 1) = (0.d0, 0.d0)
           enddo
        endif
     else
        !
        !    LSDA case
        !
        ps (:,:) = (0.d0, 0.d0)
        do is = 1, nspin0
           do js = 1, nspin0
              do ipol = 1, 3
                 ps1(ipol, is, js) = drhoout (k, is) * grho (ipol, k, js)
                 ps(is, js) = ps(is, js) + grho(ipol,k,is)*gdrho(ipol,k,js)
              enddo
              do ks = 1, nspin0
                 if (is == js .and. js == ks) then
                    a (is, js, ks) = dvxc_sr (k, is, is)
                    c (is, js, ks) = dvxc_sr (k, is, is)
                 else
                    if (is == 1) then
                       a (is, js, ks) = dvxc_sr (k, 1, 2)
                    else
                       a (is, js, ks) = dvxc_sr (k, 2, 1)
                    endif
                    if (js == 1) then
                       c (is, js, ks) = dvxc_sr (k, 1, 2)
                    else
                       c (is, js, ks) = dvxc_sr (k, 2, 1)
                    endif
                 endif
                 do ipol = 1, 3
                    ps2 (ipol, is, js, ks) = ps (is, js) * grho (ipol, k, ks)
                 enddo
                 do ls = 1, nspin0
                    if (is == js .and. js == ks .and. ks == ls) then
                       b (is, js, ks, ls) = dvxc_ss (k, is, is)
                    else
                       if (is == 1) then
                          b (is, js, ks, ls) = dvxc_ss (k, 1, 2)
                       else
                          b (is, js, ks, ls) = dvxc_ss (k, 2, 1)
                       endif
                    endif
                 enddo
              enddo
           enddo
        enddo
        do is = 1, nspin0
           do js = 1, nspin0
              dvxc (k, is) = dvxc (k, is) + dvxc_rr (k,is,js)*drhoout(k, js)
              do ipol = 1, 3
                 h (ipol, k, is) = h (ipol, k, is) + &
                      dvxc_s (k, is, js) * gdrho(ipol, k, js)
              enddo
              do ks = 1, nspin0
                 dvxc (k, is) = dvxc (k, is) + a (is, js, ks) * ps (js, ks)
                 do ipol = 1, 3
                    h (ipol, k, is) = h (ipol, k, is) + &
                         c (is, js, ks) * ps1 (ipol, js, ks)
                 enddo
                 do ls = 1, nspin0
                    do ipol = 1, 3
                       h (ipol, k, is) = h (ipol, k, is) + &
                            b (is, js, ks, ls) * ps2 (ipol, js, ks, ls)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     endif
  enddo
  ! linear variation of the second term
  do is = 1, nspin0
     call qgrad_dot (xq, nrxx, h (1, 1, is), ngm, g, nl, alat, dh)
     do k = 1, nrxx
        dvxc (k, is) = dvxc (k, is) - dh (k)
     enddo
  enddo
  IF (noncolin.AND.domag) THEN
     DO is=1,nspin0
        vgg(:,is)=dvxc(:,is)
     ENDDO
     dvxc=dvxcsave
     DO k=1,nrxx
        dvxc(k,1)=dvxc(k,1)+0.5d0*(vgg(k,1)+vgg(k,2))
        amag=sqrt(rho(k,2)**2+rho(k,3)**2+rho(k,4)**2)
        IF (amag.GT.1.d-12) THEN
           DO is=2,4
              term=(0.0_dp,0.0_dp)
              DO jpol=2,4
                 term=term+rho(k,jpol)*drho(k,jpol)
              ENDDO
              term=term*rho(k,is)/amag**2
              dvxc(k,is)=dvxc(k,is)+0.5d0*segni(k)*((vgg(k,1)-vgg(k,2)) &
                                  *rho(k,is)+vsgga(k)*(drho(k,is)-term))/amag
           ENDDO
        ENDIF
     ENDDO
  ENDIF

  deallocate (dh)
  deallocate (h)
  deallocate (gdrho)
  deallocate (rhoout)
  deallocate (drhoout)
  if (noncolin.and.domag) then
     deallocate (gdmag)
     deallocate (dvxcsave)
     deallocate (vgg)
  endif

  return

end subroutine dgradcorr
!
!--------------------------------------------------------------------
subroutine qgradient (xq, nrxx, a, ngm, g, nl, alat, ga)
  !--------------------------------------------------------------------
  !
  ! Calculates ga = \grad a in R-space (a is also in R-space)
  ! Note: gamma_only is disregarded for phonon calculations
  !
  USE kinds,          ONLY : DP
  USE constants,      ONLY : tpi
  USE control_flags,  ONLY : gamma_only
  USE fft_base,       ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE gvect,          ONLY : nlm
  
  implicit none
  integer :: nrxx, ngm, nl (ngm)
  complex(DP) :: a (nrxx), ga (3, nrxx)
  real(DP) :: g (3, ngm), alat, xq (3)
  integer :: n, ipol
  real(DP) :: tpiba
  complex(DP), allocatable :: aux (:), gaux (:)

  allocate (gaux(  nrxx))
  allocate (aux (  nrxx))

  tpiba = tpi / alat
  ! bring a(r) to G-space, a(G) ...
  aux (:) = a(:)

  CALL fwfft ('Dense', aux, dfftp)
  ! multiply by i(q+G) to get (\grad_ipol a)(q+G) ...
  do ipol = 1, 3
     gaux (:) = (0.d0, 0.d0)
     do n = 1, ngm
        gaux(nl(n)) = CMPLX(0.d0, xq (ipol) + g (ipol, n),kind=DP) * aux (nl(n))
        if (gamma_only) gaux( nlm(n) ) = conjg( gaux( nl(n) ) )
     enddo
     ! bring back to R-space, (\grad_ipol a)(r) ...

     CALL invfft ('Dense', gaux, dfftp)
     ! ...and add the factor 2\pi/a  missing in the definition of q+G
     do n = 1, nrxx
        ga (ipol, n) = gaux (n) * tpiba
     enddo
  enddo
  deallocate (aux)
  deallocate (gaux)
  
  return

end subroutine qgradient
!--------------------------------------------------------------------
subroutine qgrad_dot (xq, nrxx, a, ngm, g, nl, alat, da)
  !--------------------------------------------------------------------
  !
  ! Calculates da = \sum_i \grad_i a_i in R-space
  ! Note: gamma_only is disregarded for phonon calculations
  !
  USE kinds,          ONLY : DP
  USE constants,      ONLY : tpi
  USE control_flags,  ONLY : gamma_only
  USE fft_base,       ONLY : dfftp
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvect,          ONLY : nlm
  
  implicit none
  integer ::  nrxx, ngm, nl (ngm)
  complex(DP) :: a (3, nrxx), da (nrxx)

  real(DP) :: xq (3), g (3, ngm), alat
  integer :: n, ipol
  real(DP) :: tpiba
  complex(DP), allocatable :: aux (:)

  allocate (aux (nrxx))
  tpiba = tpi / alat
  da(:) = (0.d0, 0.d0)
  do ipol = 1, 3
     ! copy a(ipol,r) to a complex array...
     do n = 1, nrxx
        aux (n) = a (ipol, n)
     enddo
     ! bring a(ipol,r) to G-space, a(G) ...
     CALL fwfft ('Dense', aux, dfftp)
     ! multiply by i(q+G) to get (\grad_ipol a)(q+G) ...
     do n = 1, ngm
        da (nl(n)) = da (nl(n)) + &
             CMPLX(0.d0, xq (ipol) + g (ipol, n),kind=DP) * aux(nl(n))
     enddo
  enddo
  if (gamma_only) then
     !
     do n = 1, ngm
        !
        da( nlm(n) ) = conjg( da( nl(n) ) )
        !
     end do
     !
  end if

  !  bring back to R-space, (\grad_ipol a)(r) ...
  CALL invfft ('Dense', da, dfftp)
  ! ...add the factor 2\pi/a  missing in the definition of q+G and sum
  da (:) = da (:) * tpiba
  deallocate (aux)

  return
 
end subroutine qgrad_dot
