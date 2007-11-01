!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine A_h(e,h,ah)
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds, only: DP
  USE cell_base,ONLY : alat, omega, tpiba2
  USE uspp,     ONLY : vkb, nkb
  USE lsda_mod, ONLY : current_spin, nspin
  USE wvfct, ONLY: nbnd, npwx, npw, g2kin, igk
  USE wavefunctions_module,  ONLY: evc, psic
  USE scf,      ONLY : vrs, rho 
  USE gvect,    ONLY : gstart, nl, nlm, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
       nrxx, ngm, g, gg
  USE constants,  ONLY: degspin, e2, fpi
  use becmod, only: rbecp
  use cgcom
  use funct, only: dft_is_gradient
  !
  implicit none
  integer :: j, jkb, ibnd, na,nt,ih
  real(DP) :: e(nbnd)
  complex(DP) :: h(npwx,nbnd), ah(npwx,nbnd)
  !
  complex(DP) :: fp, fm
  complex(DP), pointer :: dpsic(:), drhoc(:), dvxc(:)
  real(DP), pointer  :: dv(:), drho(:)
  !
  call start_clock('a_h')
  !
  drho  => auxr
  dpsic => aux2
  drhoc => aux3
  !
  drho(:) = 0.d0
  !
  ! [(k+G)^2 - e ]psi
  do ibnd = 1,nbnd
     ! set to zero the imaginary part of h at G=0
     ! needed for numerical stability
     if (gstart==2) h(1,ibnd) = CMPLX( DBLE(h(1,ibnd)),0.d0)
     do j = 1,npw
        ah(j,ibnd) = (g2kin(j)-e(ibnd)) * h(j,ibnd)
     end do
  end do
  !     V_Loc psi
  do ibnd = 1,nbnd, 2
     dpsic(:)= (0.d0, 0.d0)
     psic(:) = (0.d0, 0.d0)
     if (ibnd.lt.nbnd) then
        ! two ffts at the same time
        do j = 1,npw
           psic (nl (igk(j))) = evc(j,ibnd) + (0.0d0,1.d0)* evc(j,ibnd+1)
           dpsic(nl (igk(j))) =   h(j,ibnd) + (0.0d0,1.d0)*   h(j,ibnd+1)
           psic (nlm(igk(j)))= CONJG(evc(j,ibnd)-(0.0d0,1.d0)* evc(j,ibnd+1))
           dpsic(nlm(igk(j)))= CONJG(  h(j,ibnd)-(0.0d0,1.d0)*   h(j,ibnd+1))
        end do
     else
        do j = 1,npw
           psic (nl (igk(j))) = evc(j,ibnd)
           dpsic(nl (igk(j))) =   h(j,ibnd)
           psic (nlm(igk(j))) = CONJG( evc(j,ibnd))
           dpsic(nlm(igk(j))) = CONJG(   h(j,ibnd))
        end do
     end if
     call cft3s( psic,nr1,nr2,nr3,nrx1,nr2,nr3,2)
     call cft3s(dpsic,nr1,nr2,nr3,nrx1,nr2,nr3,2)
     do j = 1,nrxx
        drho(j) = drho(j) - 2.0d0*degspin/omega *   &
              DBLE(psic(j)*CONJG(dpsic(j)))
        dpsic(j) = dpsic(j) * vrs(j,current_spin)
     end do
     call cft3s(dpsic,nr1,nr2,nr3,nrx1,nr2,nr3,-2)
     if (ibnd.lt.nbnd) then
        ! two ffts at the same time
        do j = 1,npw
           fp = (dpsic (nl(igk(j))) + dpsic (nlm(igk(j))))*0.5d0
           fm = (dpsic (nl(igk(j))) - dpsic (nlm(igk(j))))*0.5d0
           ah(j,ibnd  ) = ah(j,ibnd)  +CMPLX( DBLE(fp), AIMAG(fm))
           ah(j,ibnd+1) = ah(j,ibnd+1)+CMPLX(AIMAG(fp),- DBLE(fm))
        end do
     else
        do j = 1,npw
           ah(j,ibnd) = ah(j,ibnd)  + dpsic (nl(igk(j)))
        end do
     end if
  end do
  !
  nullify(dpsic)
  ! V_NL psi
  call pw_gemm ('Y', nkb, nbnd, npw, vkb, npwx, h, npwx, rbecp, nkb)
  if (nkb.gt.0) call add_vuspsi (npwx, npw, nbnd, h, ah)
  !
  do j = 1,nrxx
     drhoc(j) = CMPLX(drho(j),0.d0)
  end do
  call cft3(drhoc,nr1,nr2,nr3,nrx1,nr2,nr3,-1)
  !
  ! drho is deltarho(r), drhoc is deltarho(g)
  !
  !  mu'(n(r)) psi(r) delta psi(r)
  !
  dvxc  => aux2
  do j = 1,nrxx
     dvxc(j) = drho(j)*dmuxc(j)
  end do
  !
  !  add gradient correction contribution (if any)
  !
  call start_clock('dgradcorr')
  if (dft_is_gradient() ) call dgradcor1  &
       (rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s,            &
        drho, drhoc, nr1,nr2,nr3, nrx1, nrx2, nrx3, nrxx, nspin, &
        nl, nlm, ngm, g, alat, omega, dvxc)
  call stop_clock('dgradcorr')
  nullify (drho)
  !
  !  1/|r-r'| * psi(r') delta psi(r')
  !
  ! gstart is the first nonzero G vector (needed for parallel execution)
  !
  if (gstart==2) drhoc(nl(1)) = 0.d0
  !
  do j = gstart,ngm
     drhoc(nl (j)) = e2*fpi*drhoc(nl(j))/ (tpiba2*gg(j))
     drhoc(nlm(j)) = CONJG(drhoc(nl (j)))
  end do
  call cft3(drhoc,nr1,nr2,nr3,nrx1,nr2,nr3,+1)
  !
  ! drhoc now contains deltaV_hartree
  !
  dv => auxr
  do j = 1,nrxx
     dv(j) = -  DBLE(dvxc(j)) - DBLE(drhoc(j))
  end do
  !
  call vloc_psi(npwx, npw, nbnd, evc, dv, ah)
  !
  ! set to zero the imaginary part of ah at G=0
  ! needed for numerical stability
  if (gstart.eq.2) then
     do ibnd = 1, nbnd
        ah(1,ibnd) = CMPLX( DBLE(ah(1,ibnd)),0.d0)
     end do
  end if
  !
  call stop_clock('a_h')
  !
  return
end subroutine A_h
