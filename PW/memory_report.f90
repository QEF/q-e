!
! Copyright (C) 2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE memory_report()
  !----------------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  USE wvfct,     ONLY : npwx, nbnd, nbndx
  USE basis,     ONLY : natomwfc
  USE gvect,     ONLY : ngl, nr1, nr2, nr3, nrxx, ngm
  USE uspp,      ONLY : nkb
  USE ldaU,      ONLY : lda_plus_u
  USE lsda_mod,  ONLY : nspin
  USE noncollin_module,     ONLY : npol
  USE wavefunctions_module, ONLY : evc
  USE control_flags, ONLY: isolve, nmix, gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: Mb=1024*1024, complex_size=16, real_size=8
  INTEGER :: size
  !
  !
  WRITE( stdout, '(/5x,"Largest allocated arrays",5x,"est. size (Mb)", &
                   &5x,"dimensions")')
  WRITE( stdout, '(8x,"Kohn-Sham Wavefunctions   ",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,")")') &
     DBLE(complex_size*nbnd*npol*npwx)/Mb, npwx*npol,nbnd
  IF ( lda_plus_u ) &
     WRITE( stdout, '(8x,"Atomic wavefunctions      ",f10.2," Mb", &
                    & 5x,"(",i7,",",i4,")")') &
     DBLE(complex_size*natomwfc*npol*npwx)/Mb, npwx*npol,natomwfc
  WRITE( stdout, '(8x,"NL pseudopotentials       ",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,")")') &
     DBLE(complex_size*nkb*npwx)/Mb, npwx, nkb
  IF ( nspin == 2 ) THEN
     WRITE( stdout, '(8x,"Each V/rho on FFT grid    ",f10.2," Mb", &
                    & 5x,"(",i7,",",i4,")")') &
                    DBLE(complex_size*nrxx*nspin)/Mb, nrxx, nspin
  ELSE
     WRITE( stdout, '(8x,"Each V/rho on FFT grid    ",f10.2," Mb", &
                    & 5x,"(",i7,")")') DBLE(complex_size*nrxx)/Mb, nrxx
  END IF
  WRITE( stdout, '(8x,"Each G-vector array       ",f10.2," Mb", &
                 & 5x,"(",i7,")")') DBLE(real_size*ngm)/Mb, ngm
  WRITE( stdout, '(8x,"G-vector shells           ",f10.2," Mb", &
                 & 5x,"(",i7,")")') DBLE(real_size*ngl)/Mb, ngl
  !
  WRITE( stdout, '(5x,"Largest temporary arrays",5x,"est. size (Mb)", &
                   &5x,"dimensions")')
  IF ( gamma_only) THEN
     size = real_size
  ELSE
     size = complex_size
  END IF
  !
  IF ( isolve == 0 ) THEN
     WRITE( stdout, '(8x,"Auxiliary wavefunctions   ",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,")")') &
     DBLE(size*nbndx*npol*npwx)/Mb, npwx*npol, nbndx
     WRITE( stdout, '(8x,"Each subspace H/S matrix  ",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,")")') &
     DBLE(size*nbndx*nbndx)/Mb, nbndx, nbndx
  ELSE
     WRITE( stdout, '(8x,"Each subspace H/S matrix  ",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,")")') &
     DBLE(size*nbnd*nbnd)/Mb, nbnd, nbnd
  ENDIF
  !
  IF ( npol > 1 ) THEN
     WRITE( stdout, '(8x,"Each <psi_i|beta_j> matrix",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,",",i4,")")') &
     DBLE(size*nkb*npol*nbnd)/Mb, nkb, npol, nbnd
  ELSE
     WRITE( stdout, '(8x,"Each <psi_i|beta_j> matrix",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,")")') &
     DBLE(size*nkb*nbnd)/Mb, nkb, nbnd
  END IF
  !
  WRITE( stdout, '(8x,"Arrays for rho mixing     ",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,")")') &
     DBLE(complex_size*nrxx*nmix)/Mb, nrxx, nmix
  !
  RETURN
  !
END subroutine memory_report
