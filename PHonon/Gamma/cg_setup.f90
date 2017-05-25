!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE cg_setup
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY: DP
  USE cell_base,  ONLY: bg
  USE ions_base,  ONLY: nat, ntyp => nsp, ityp, tau, amass
  USE scf,        ONLY: rho, rho_core, v, vltot, vrs, kedtau
  USE uspp,       ONLY: vkb, nlcc_any
  USE uspp_param, ONLY: upf
  USE mp_global,  ONLY: kunit
  USE wavefunctions_module,  ONLY: evc
  USE io_files,   ONLY: prefix, iunpun, iunres, diropn
  USE funct,      ONLY: dft_is_gradient, dmxc
  USE dfunct,     ONLY: newd
  USE fft_base,   ONLY: dfftp
  USE gvect,      ONLY: g, ngm, eigts1, eigts2, eigts3
  USE gvecs,      ONLY: doublegrid
  USE klist,      ONLY: xk, ngk, igk_k
  USE lsda_mod,   ONLY: nspin, current_spin
  USE vlocal,     ONLY: strf
  USE wvfct,      ONLY: nbnd, npwx
  USE gvecw,      ONLY: gcutw
  USE cgcom
  !
  IMPLICIT NONE
  !
  INTEGER :: i, l, nt, ik
  LOGICAL :: exst
  CHARACTER (len=256) :: filint
  REAL(DP) :: rhotot
  INTEGER       :: ndr, kunittmp, ierr
  REAL(DP) :: edum(1,1), wdum(1,1)
  !
  CALL start_clock('cg_setup')
  !
  !  sum self-consistent part (vr) and local part (vltot) of potential
  !
  CALL set_vrs(vrs,vltot,v%of_r,kedtau, v%kin_r, dfftp%nnr,nspin,doublegrid)
  !
  !  allocate memory for various arrays
  !
  ALLOCATE  (dmuxc( dfftp%nnr))
  ALLOCATE  (dvpsi( npwx, nbnd))
  ALLOCATE  ( dpsi( npwx, nbnd))
  ALLOCATE  ( auxr( dfftp%nnr))
  ALLOCATE  ( aux2( dfftp%nnr))
  ALLOCATE  ( aux3( dfftp%nnr))
  !
  !  allocate memory for gradient corrections (if needed)
  !
  IF ( dft_is_gradient() ) THEN
     ALLOCATE  ( dvxc_rr(dfftp%nnr,nspin,nspin))
     ALLOCATE  ( dvxc_sr(dfftp%nnr,nspin,nspin))
     ALLOCATE  ( dvxc_ss(dfftp%nnr,nspin,nspin))
     ALLOCATE  ( dvxc_s (dfftp%nnr,nspin,nspin))
     ALLOCATE  ( grho   (3, dfftp%nnr, nspin))
  ENDIF
  !
  !
  !  initialize structure factor array
  !
  CALL struc_fact ( nat, tau, ntyp, ityp, ngm, g, bg,               &
       &                  dfftp%nr1, dfftp%nr2, dfftp%nr3, strf, eigts1, eigts2, eigts3 )
  !
  !  compute drhocore/dtau for each atom type (if needed)
  !
  nlcc_any = any  ( upf(1:ntyp)%nlcc )
  !!! if (nlcc_any) call set_drhoc(xq, drc)
  !
  !  local potential
  !
  CALL init_vloc
  !
  CALL init_us_1
  !
  CALL newd
  !
  !  derivative of the xc potential
  !
  dmuxc(:) = 0.d0
  DO i = 1,dfftp%nnr
     rhotot = rho%of_r(i,current_spin)+rho_core(i)
     IF ( rhotot> 1.d-30 ) dmuxc(i)= dmxc( rhotot)
     IF ( rhotot<-1.d-30 ) dmuxc(i)=-dmxc(-rhotot)
  ENDDO
  !
  !  initialize data needed for gradient corrections
  !
  CALL cg_setupdgc
  !
  iunres=88
  !
  !   open the wavefunction file (already existing)
  !
  lrwfc=2*nbnd*npwx
  CALL diropn(iunpun, 'wfc',lrwfc,exst)
  IF(.not.exst) THEN
     CALL errore('main','file '//trim(prefix) // '.wfc not found',1)
  ENDIF
  !  read wave functions and calculate indices
  !
  ik=1
  CALL davcio(evc,lrwfc,iunpun,ik,-1)
  IF ( exst ) THEN
     CLOSE(unit=iunpun,status='keep')
  ELSE
     CLOSE(unit=iunpun,status='delete')
  ENDIF
  !
  !  Kleinman-Bylander PPs
  !
  CALL init_us_2 (ngk(ik), igk_k(1,ik), xk(1,ik), vkb)
  !
  CALL stop_clock('cg_setup')
  !
  RETURN
END SUBROUTINE cg_setup
