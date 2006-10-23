!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
subroutine stress
  !----------------------------------------------------------------------
  !
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : omega, alat, at, bg
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, zv
  USE constants,     ONLY : uakbar
  USE ener,          ONLY : etxc, vtxc
  USE force_mod,     ONLY : sigma
  USE gvect,         ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                            nrxx, nl, g, gg, gcutm
  USE ldaU,          ONLY : lda_plus_u
  USE lsda_mod,      ONLY : nspin
  USE scf,           ONLY : rho, rho_core
  USE control_flags, ONLY : iverbosity
  USE wvfct,         ONLY : gamma_only
  !
  implicit none
  !
  real(DP) :: sigmakin (3, 3), sigmaloc (3, 3), sigmahar (3, 3), &
       sigmaxc (3, 3), sigmaxcc (3, 3), sigmaewa (3, 3), sigmanlc (3, 3), &
       sigmabare (3, 3), sigmah (3, 3)
  integer :: l, m
  !
  call start_clock ('stress')
  WRITE( stdout, '(//5x,"entering subroutine stress ..."/)')

  !
  !   contribution from local  potential
  !
  call stres_loc (sigmaloc)
  !
  !  hartree contribution
  !
  call stres_har (sigmahar)
  !
  !  xc contribution (diagonal)
  !
  sigmaxc(:,:) = 0.d0
  do l = 1, 3
     sigmaxc (l, l) = - (etxc - vtxc) / omega
  enddo
  !
  !  xc contribution: add gradient corrections (non diagonal)
  !
  call stres_gradcorr (rho, rho_core, nspin, nr1, nr2, nr3, nrx1, &
       nrx2, nrx3, nrxx, nl, ngm, g, alat, omega, sigmaxc)
  !
  ! core correction contribution
  !
  call stres_cc (sigmaxcc)
  !
  !  ewald contribution
  !
  call stres_ewa (alat, nat, ntyp, ityp, zv, at, bg, tau, omega, g, &
       gg, ngm, gstart, gamma_only, gcutm, sigmaewa)
  !
  !  kinetic + nonlocal contribuition
  !
  call stres_knl (sigmanlc, sigmakin)
  !
  do l = 1, 3
     do m = 1, 3
        sigmabare (l, m) = sigmaloc (l, m) + sigmanlc (l, m)
     enddo
  enddo
  !
  !  Hubbard contribution
  !
  sigmah(:,:) = 0.d0
  if (lda_plus_u) call stres_hub(sigmah)
  !
  sigma(:,:) = sigmakin(:,:) + sigmaloc(:,:) + sigmahar(:,:) + &
               sigmaxc(:,:) + sigmaxcc(:,:) + sigmaewa(:,:) + &
               sigmanlc(:,:) + sigmah(:,:)
  !
  ! write results in Ryd/(a.u.)^3 and in kbar
  !

  WRITE( stdout, 9000) (sigma(1,1) + sigma(2,2) + sigma(3,3)) * uakbar / 3d0,  &
                  (sigma(l,1), sigma(l,2), sigma(l,3),                    &
            sigma(l,1)*uakbar, sigma(l,2)*uakbar, sigma(l,3)*uakbar, l=1,3)

  if (iverbosity.ge.1) WRITE( stdout, 9005) &
     (sigmakin(l,1)*uakbar,sigmakin(l,2)*uakbar,sigmakin(l,3)*uakbar, l=1,3),&
     (sigmaloc(l,1)*uakbar,sigmaloc(l,2)*uakbar,sigmaloc(l,3)*uakbar, l=1,3),&
     (sigmanlc(l,1)*uakbar,sigmanlc(l,2)*uakbar,sigmanlc(l,3)*uakbar, l=1,3),&
     (sigmahar(l,1)*uakbar,sigmahar(l,2)*uakbar,sigmahar(l,3)*uakbar, l=1,3),&
     (sigmaxc (l,1)*uakbar,sigmaxc (l,2)*uakbar,sigmaxc (l,3)*uakbar, l=1,3),&
     (sigmaxcc(l,1)*uakbar,sigmaxcc(l,2)*uakbar,sigmaxcc(l,3)*uakbar, l=1,3),&
     (sigmaewa(l,1)*uakbar,sigmaewa(l,2)*uakbar,sigmaewa(l,3)*uakbar, l=1,3),&
     (sigmah  (l,1)*uakbar,sigmah  (l,2)*uakbar,sigmah  (l,3)*uakbar, l=1,3)

  call stop_clock ('stress')

  return
9000 format (10x,'total   stress  (Ry/bohr**3) ',18x,'(kbar)', &
             &5x,'P=',f8.2/3 (3f13.8,4x,3f10.2/))
9005 format &
         &  (5x,'kinetic stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'local   stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'nonloc. stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'hartree stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'exc-cor stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'corecor stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'ewald   stress (kbar)',3f10.2/2(26x,3f10.2/)/ &
         &   5x,'hubbard stress (kbar)',3f10.2/2(26x,3f10.2/)/ )

  end subroutine stress

