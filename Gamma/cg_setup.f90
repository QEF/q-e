!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine cg_setup
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  use parameters, only: DP
  use pwcom
  USE wavefunctions_module,  ONLY: evc
  use io_files, only: prefix, iunpun, iunres
  use cgcom
  use funct
  !
  implicit none
  !
  integer :: i, l, nt, kpoint
  logical :: exst
  character (len=256) :: filint
  real(kind=DP) :: rhotot, dmxc
  external dmxc
  !
  call start_clock('cg_setup')
  !
  !  convert masses to atomic units
  !
  call DSCAL(ntyp,amconv,amass,1)
  !
  !  sum self-consistent part (vr) and local part (vltot) of potential
  !
  call set_vrs(vrs,vltot,vr,nrxx,nspin,doublegrid)
  !
  !  allocate memory for various arrays
  !
  allocate  (dmuxc( nrxx))    
  allocate  (dvpsi( npwx, nbnd))    
  allocate  ( dpsi( npwx, nbnd))    
  allocate  ( auxr( nrxx))    
  allocate  ( aux2( nrxx))    
  allocate  ( aux3( nrxx))    
  !
  !  allocate memory for gradient corrections (if needed)
  !
  if (igcx.ne.0 .or. igcc.ne.0) then
     allocate  ( dvxc_rr(nrxx,nspin,nspin))    
     allocate  ( dvxc_sr(nrxx,nspin,nspin))    
     allocate  ( dvxc_ss(nrxx,nspin,nspin))    
     allocate  ( dvxc_s (nrxx,nspin,nspin))    
     allocate  ( grho   (3, nrxx, nspin))    
  end if
  !
  !
  !  initialize structure factor array
  !
  call struc_fact ( nat, tau, ntyp, ityp, ngm, g, bg,               &
       &                  nr1, nr2, nr3, strf, eigts1, eigts2, eigts3 )
  !
  !  compute drhocore/dtau for each atom type (if needed)
  !
  nlcc_any = .false.
  do nt=1,ntyp
     nlcc_any = nlcc_any .or. nlcc(nt)
  end do
!!!      if (nlcc_any) call set_drhoc(xq)
  !
  !  local potential
  !
  call init_vloc
  !
  call convert_to_num                                               &
       &   (ntyp,numeric,ndm,mesh,r,lmaxx,lmax,lloc,nnl,aps,alps,vnl)
  !
  call init_us_1
  !
  call newd
  !
  !  derivative of the xc potential
  !
  dmuxc(:) = 0.d0
  do i = 1,nrxx
     rhotot = rho(i,current_spin)+rho_core(i)
     if ( rhotot.gt. 1.d-30 ) dmuxc(i)= dmxc( rhotot)
     if ( rhotot.lt.-1.d-30 ) dmuxc(i)=-dmxc(-rhotot)
  end do
  !
  !  initialize data needed for gradient corrections
  !
  call cg_setupdgc
  !
  iunres=88
  !
  !   open the wavefunction file (already existing)
  !
  lrwfc=2*nbnd*npwx
  filint = trim(prefix) //'.wfc'
  call diropn(iunpun,filint,lrwfc,exst)
  if(.not.exst) call errore('main','file '//filint//' not found',1)
  !
  !  read wave functions and calculate indices
  !
  kpoint=1
  call davcio(evc,lrwfc,iunpun,kpoint,-1)
  close(unit=iunpun,status='keep')
  call gk_sort (xk(1,kpoint),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)
  !
  !  Kleinman-Bylander PPs
  !
  call init_us_2 (npw, igk, xk(1,kpoint), vkb)
  !
  call stop_clock('cg_setup')
  !
  return
end subroutine cg_setup
