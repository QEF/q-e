!
! Copyright (C) 2003-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE cg_setup
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY: DP
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, amass
  USE pwcom
  USE uspp_param, ONLY: upf
  USE mp_global,        ONLY : kunit
  USE wavefunctions_module,  ONLY: evc
  USE io_files, ONLY: prefix, iunpun, iunres
  USE cgcom
  USE funct, only : dft_is_gradient, dmxc
  !
  IMPLICIT NONE
  !
  INTEGER :: i, l, nt, kpoint
  LOGICAL :: exst
  CHARACTER (len=256) :: filint
  REAL(DP) :: rhotot
  INTEGER       :: ndr, kunittmp, ierr
  REAL(DP) :: edum(1,1), wdum(1,1)
  !
  CALL start_clock('cg_setup')
  !
  !  convert masses to atomic units
  !
  CALL DSCAL(ntyp,amconv,amass,1)
  !
  !  sum self-consistent part (vr) and local part (vltot) of potential
  !
  CALL set_vrs(vrs,vltot,vr,nrxx,nspin,doublegrid)
  !
  !  allocate memory for various arrays
  !
  ALLOCATE  (dmuxc( nrxx))    
  ALLOCATE  (dvpsi( npwx, nbnd))    
  ALLOCATE  ( dpsi( npwx, nbnd))    
  ALLOCATE  ( auxr( nrxx))    
  ALLOCATE  ( aux2( nrxx))    
  ALLOCATE  ( aux3( nrxx))    
  !
  !  allocate memory for gradient corrections (if needed)
  !
  IF ( dft_is_gradient() ) THEN
     ALLOCATE  ( dvxc_rr(nrxx,nspin,nspin))    
     ALLOCATE  ( dvxc_sr(nrxx,nspin,nspin))    
     ALLOCATE  ( dvxc_ss(nrxx,nspin,nspin))    
     ALLOCATE  ( dvxc_s (nrxx,nspin,nspin))    
     ALLOCATE  ( grho   (3, nrxx, nspin))    
  END IF
  !
  !
  !  initialize structure factor array
  !
  CALL struc_fact ( nat, tau, ntyp, ityp, ngm, g, bg,               &
       &                  nr1, nr2, nr3, strf, eigts1, eigts2, eigts3 )
  !
  !  compute drhocore/dtau for each atom type (if needed)
  !
  nlcc_any = ANY  ( upf(1:ntyp)%nlcc )
  !!! if (nlcc_any) call set_drhoc(xq)
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
  DO i = 1,nrxx
     rhotot = rho(i,current_spin)+rho_core(i)
     IF ( rhotot.GT. 1.d-30 ) dmuxc(i)= dmxc( rhotot)
     IF ( rhotot.LT.-1.d-30 ) dmuxc(i)=-dmxc(-rhotot)
  END DO
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
  IF(.NOT.exst) THEN 
     CALL errore('main','file '//trim(prefix) // '.wfc not found',1)
  END IF
  !  read wave functions and calculate indices
  !
  kpoint=1
  CALL davcio(evc,lrwfc,iunpun,kpoint,-1)
  IF ( exst ) THEN
     CLOSE(unit=iunpun,status='keep')
  ELSE 
     CLOSE(unit=iunpun,status='delete')
  END IF
  CALL gk_sort (xk(1,kpoint),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)
  !
  !  Kleinman-Bylander PPs
  !
  CALL init_us_2 (npw, igk, xk(1,kpoint), vkb)
  !
  CALL stop_clock('cg_setup')
  !
  RETURN
END SUBROUTINE cg_setup
