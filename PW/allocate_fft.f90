!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine allocate_fft
  !-----------------------------------------------------------------------
  !     This routine computes the data structure associated to the FFT
  !     grid and allocate memory for all the arrays which depend upon
  !     these dimensions
  !
  USE io_global, ONLY : stdout
  USE gvect,     ONLY : nr1, nr2, nr3, nrxx, ngm, g, gg, nl, nlm, &
       ig1, ig2, ig3, eigts1, eigts2, eigts3, igtongl
  USE gsmooth,   ONLY : nr1s,nr2s,nr3s,nrxxs,ngms, nls, nlsm, doublegrid
  USE ions_base, ONLY : nat
  USE lsda_mod,  ONLY : nspin
  USE scf,       ONLY : rho, rhog, vr, vltot, vrs, rho_core, rhog_core, &
                        tauk, taukg, kedtau, kedtaur, allocate_scf_type
  USE vlocal,    ONLY : vnew
  USE wvfct,     ONLY : gamma_only
  USE noncollin_module, ONLY : pointlist, factlist, r_loc, &
      report, i_cons, noncolin, npol
  USE wavefunctions_module, ONLY : psic, psic_nc
  USE funct,     ONLY: dft_is_meta
  implicit none
  !
  !     determines the data structure for fft arrays
  !
  call data_structure( gamma_only )
  !
  if (nrxx.lt.ngm) then
     WRITE( stdout, '(/,4x," nr1=",i4," nr2= ", i4, " nr3=",i4, &
          &" nrxx = ",i8," ngm=",i8)') nr1, nr2, nr3, nrxx, ngm
     call errore ('allocate_fft', 'the nr"s are too small!', 1)

  endif
  if (nrxxs.lt.ngms) then
     WRITE( stdout, '(/,4x," nr1s=",i4," nr2s= ", i4, " nr3s=",i4, &
          &" nrxxs = ",i8," ngms=",i8)') nr1s, nr2s, nr3s, nrxxs, ngms
     call errore ('allocate_fft', 'the nrs"s are too small!', 1)

  endif
  if (ngm  <= 0) call errore ('allocate_fft', 'wrong ngm', 1)
  if (ngms <= 0) call errore ('allocate_fft', 'wrong ngms', 1)
  if (nrxx <= 0) call errore ('allocate_fft', 'wrong nrxx', 1)
  if (nrxxs<= 0) call errore ('allocate_fft', 'wrong nrxxs', 1)
  if (nspin<= 0) call errore ('allocate_fft', 'wrong nspin', 1)
  !
  !     Allocate memory for all kind of stuff.
  !
  allocate (g( 3, ngm))    
  allocate (gg( ngm))    
  allocate (nl(  ngm))    
  if (gamma_only) allocate (nlm(ngm))
  allocate (igtongl(  ngm))    
  allocate (ig1(  ngm))    
  allocate (ig2(  ngm))    
  allocate (ig3(  ngm))    

  call allocate_scf_type(rho)
  ALLOCATE( rhog( ngm, nspin ) )
  allocate (vr( nrxx,nspin))    
  allocate (vltot( nrxx))    
  allocate (vnew  ( nrxx, nspin))    
  allocate (rho_core( nrxx))
  if (dft_is_meta() ) then
     allocate ( tauk(nrxx,nspin) )
     allocate ( kedtau(nrxxs,nspin) )
     allocate ( kedtaur(nrxx,nspin) )
     ALLOCATE ( taukg( ngm, nspin ) )
  else
     allocate ( tauk(1,nspin) )
     allocate ( kedtau(1,nspin) )
     allocate ( kedtaur(1,nspin) )
     ALLOCATE ( taukg( 1, nspin ) )
  end if
  ALLOCATE( rhog_core( ngm ) )
  allocate (psic( nrxx))    
  allocate (vrs( nrxx, nspin))    
  if (doublegrid) then
     allocate (nls( ngms))    
     if (gamma_only) allocate (nlsm(ngm))
  else
     nls => nl
     if (gamma_only) nlsm=> nlm
  endif

  if (noncolin) allocate (psic_nc( nrxx, npol))    

  if ( ((report.ne.0).or.(i_cons.ne.0)) .and. (noncolin) .or. (i_cons.eq.1) ) then
!
! In order to print out local quantities, integrated around the atoms,
! we need the following variables
!
     allocate(pointlist(nrxx))
     allocate(factlist(nrxx))
     allocate(r_loc(nat))
     call make_pointlists
  endif

  return
end subroutine allocate_fft
