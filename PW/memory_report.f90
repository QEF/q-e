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
  USE noncollin_module,     ONLY : npol
  USE wavefunctions_module, ONLY : evc
  USE control_flags, ONLY: isolve, nmix
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: Mb=1024*1024
  !
  !
  WRITE( stdout, '(/5x,"Largest allocated arrays,  est. size (Mb),  dimensions")')
  WRITE( stdout, '(8x,"Kohn-Sham Wavefunctions ",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,")")') &
     DBLE(16*nbnd*npol*npwx)/Mb, npwx*npol,nbnd
  IF ( lda_plus_u ) &
     WRITE( stdout, '(8x,"Atomic wavefunctions    ",f10.2," Mb", &
                    & 5x,"(",i7,",",i4,")")') &
     DBLE(16*natomwfc*npol*npwx)/Mb, npwx*npol,natomwfc
  WRITE( stdout, '(8x,"NL pseudopotentials     ",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,")")') &
     DBLE(16*nkb*npwx)/Mb, npwx, nkb
  WRITE( stdout, '(8x,"Each V/rho on FFT grid  ",f10.2," Mb", &
                 & 5x,"(",i7,")")') DBLE(16*nrxx)/Mb, nrxx
  WRITE( stdout, '(8x,"Each G-vector array     ",f10.2," Mb", &
                 & 5x,"(",i7,")")') DBLE(8*ngm)/Mb, ngm
  WRITE( stdout, '(8x,"G-vector shells         ",f10.2," Mb", &
                 & 5x,"(",i7,")")') DBLE(8*ngl)/Mb, ngl
  WRITE( stdout, '(5x,"Largest temporary arrays,  est. size (Mb),  dimensions")')
  IF ( isolve == 0 ) THEN
     WRITE( stdout, '(8x,"Auxiliary wavefunctions ",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,")")') &
     DBLE(16*nbndx*npol*npwx)/Mb, npwx*npol, nbndx
     WRITE( stdout, '(8x,"Each subspace H/S matrix",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,")")') &
     DBLE(16*nbndx*nbndx)/Mb, nbndx, nbndx
  ELSE
     WRITE( stdout, '(8x,"Each subspace H/S matrix",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,")")') &
     DBLE(16*nbnd*nbnd)/Mb, nbnd, nbnd
  ENDIF
  WRITE( stdout, '(8x,"Arrays for rho mixing   ",f10.2," Mb", &
                 & 5x,"(",i7,",",i4,")")') &
     DBLE(16*nrxx*nmix)/Mb, nrxx, nmix
  !
  RETURN
  !
END subroutine memory_report
