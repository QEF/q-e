!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE memory_report()
  !----------------------------------------------------------------------------
  !
  USE io_global, ONLY : stdout
  USE wvfct,     ONLY : npwx, nbnd, nbndx
  USE basis,     ONLY : natomwfc
  USE fft_base,  ONLY : dfftp
  USE gvect,     ONLY : ngl, ngm
  USE uspp,      ONLY : nkb
  USE ldaU,      ONLY : lda_plus_u, U_projection, nwfcU
  USE fixed_occ, ONLY : one_atom_occupations
  USE wannier_new,ONLY: use_wannier
  USE lsda_mod,  ONLY : nspin
  USE noncollin_module,     ONLY : npol
  USE control_flags, ONLY: isolve, nmix, gamma_only, lscf
  USE mp_diag,   ONLY : np_ortho
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: Mb=1024*1024, complex_size=16, real_size=8
  INTEGER :: g_size, nbnd_l
  !
  ! the conversions to double prevent integer overflow in very large run
  !
  WRITE( stdout, '(/5x,"Largest allocated arrays",5x,"est. size (Mb)", &
                   &5x,"dimensions")')
  WRITE( stdout, '(8x,"Kohn-Sham Wavefunctions   ",f10.2," Mb", &
                 & 5x,"(",i8,",",i5,")")') &
     complex_size*nbnd*npol*DBLE(npwx)/Mb, npwx*npol,nbnd
  IF ( one_atom_occupations .OR. use_wannier ) &
     WRITE( stdout, '(8x,"Atomic wavefunctions      ",f10.2," Mb", &
                    & 5x,"(",i8,",",i5,")")') &
   & complex_size*natomwfc*npol*DBLE(npwx)/Mb, npwx*npol,natomwfc
  IF ( lda_plus_u .AND. U_projection .NE. 'pseudo' ) &
     WRITE( stdout, '(8x,"Atomic Hubbard wavefuncts ",f10.2," Mb", &
                    & 5x,"(",i8,",",i5,")")') &
   & complex_size*nwfcU*npol*DBLE(npwx)/Mb, npwx*npol,nwfcU
  WRITE( stdout, '(8x,"NL pseudopotentials       ",f10.2," Mb", &
                 & 5x,"(",i8,",",i5,")")') &
     complex_size*nkb*DBLE(npwx)/Mb, npwx, nkb
  IF ( nspin == 2 ) THEN
     WRITE( stdout, '(8x,"Each V/rho on FFT grid    ",f10.2," Mb", &
                    & 5x,"(",i8,",",i4,")")') &
                    complex_size*nspin*DBLE(dfftp%nnr)/Mb, dfftp%nnr, nspin
  ELSE
     WRITE( stdout, '(8x,"Each V/rho on FFT grid    ",f10.2," Mb", &
                    & 5x,"(",i8,")")') DBLE(complex_size*dfftp%nnr)/Mb, dfftp%nnr
  END IF
  WRITE( stdout, '(8x,"Each G-vector array       ",f10.2," Mb", &
                 & 5x,"(",i8,")")') DBLE(real_size*ngm)/Mb, ngm
  WRITE( stdout, '(8x,"G-vector shells           ",f10.2," Mb", &
                 & 5x,"(",i8,")")') DBLE(real_size*ngl)/Mb, ngl
  !
  WRITE( stdout, '(5x,"Largest temporary arrays",5x,"est. size (Mb)", &
                   &5x,"dimensions")')
  IF ( gamma_only) THEN
     g_size = real_size
  ELSE
     g_size = complex_size
  END IF
  !
  IF ( isolve == 0 ) THEN
     WRITE( stdout, '(8x,"Auxiliary wavefunctions   ",f10.2," Mb", &
                 & 5x,"(",i8,",",i5,")")') &
     g_size*nbndx*npol*DBLE(npwx)/Mb, npwx*npol, nbndx
  ENDIF
  ! nbnd_l : estimated dimension of distributed matrices
  nbnd_l = nbndx/np_ortho(1)
  WRITE( stdout, '(8x,"Each subspace H/S matrix  ",f10.2," Mb", &
                 & 5x,"(",i8,",",i5,")")') &
     g_size*nbnd_l*DBLE(nbnd_l)/Mb, nbnd_l, nbnd_l
  !
  IF ( npol > 1 ) THEN
     WRITE( stdout, '(8x,"Each <psi_i|beta_j> matrix",f10.2," Mb", &
                 & 5x,"(",i8,",",i4,",",i5,")")') &
     g_size*nkb*DBLE(npol*nbnd)/Mb, nkb, npol, nbnd
  ELSE
     WRITE( stdout, '(8x,"Each <psi_i|beta_j> matrix",f10.2," Mb", &
                 & 5x,"(",i8,",",i5,")")') &
     g_size*nkb*DBLE(nbnd)/Mb, nkb, nbnd
  END IF
  !
  IF ( lscf) WRITE( stdout, &
     '(8x,"Arrays for rho mixing     ",f10.2," Mb", 5x,"(",i8,",",i5,")")') &
     complex_size*dfftp%nnr*DBLE(nmix)/Mb, dfftp%nnr, nmix
  !
  RETURN
  !
END subroutine memory_report
