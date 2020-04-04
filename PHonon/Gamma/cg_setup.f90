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
  USE wavefunctions,  ONLY: evc
  USE io_files,   ONLY: prefix, iunpun, iunres, diropn
  USE funct,      ONLY: dft_is_gradient
  USE dfunct,     ONLY: newd
  USE fft_base,   ONLY: dfftp
  USE gvect,      ONLY: g, ngm, eigts1, eigts2, eigts3
  USE gvecs,      ONLY: doublegrid
  USE klist,      ONLY: xk, ngk, igk_k
  USE lsda_mod,   ONLY: nspin, current_spin
  USE vlocal,     ONLY: strf
  USE wvfct,      ONLY: nbnd, npwx
  USE gvecw,      ONLY: gcutw
  USE gc_lr, ONLY:  grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
  USE cgcom, ONLY: dmuxc, dvpsi, dpsi, auxr, aux2, aux3, lrwfc
  !
  IMPLICIT NONE
  !
  INTEGER :: i, l, nt, ik
  LOGICAL :: exst
  CHARACTER (len=256) :: filint
  INTEGER  :: ndr, ierr
  REAL(DP) :: edum(1,1), wdum(1,1)
  REAL(DP), DIMENSION(dfftp%nnr) :: rhotot
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
  !  compute drhocore/dtau for each atom type (if needed)
  !
  nlcc_any = any  ( upf(1:ntyp)%nlcc )
  !!! if (nlcc_any) call set_drhoc(xq, drc)
  !
  rhotot(:) = rho%of_r(:,1) + rho_core(:)
  !
  CALL dmxc_lda( dfftp%nnr, rhotot, dmuxc )
  !
  !  initialize data needed for gradient corrections
  !
  CALL setup_dgc( ) 
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
