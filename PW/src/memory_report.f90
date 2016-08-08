!
! Copyright (C) 2007-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE memory_report()
  !----------------------------------------------------------------------------
  !
  ! Very rough estimate of the dynamical memory allocated by the pw.x code
  ! Some dimensions are estimated since they are not known when this routine 
  ! is called (must be called after the first steps of initialization are done
  ! but before large arrays are actually allocated). 
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : dp 
  USE constants, ONLY : tpi, fpi
  USE wvfct,     ONLY : nbnd, nbndx
  USE basis,     ONLY : natomwfc
  USE cell_base, ONLY : omega
  USE fft_base,  ONLY : dfftp
  USE gvect,     ONLY : ngm
  USE gvecw,     ONLY : ecutwfc
  USE uspp,      ONLY : nkb
  USE ldaU,      ONLY : lda_plus_u, U_projection, nwfcU
  USE fixed_occ, ONLY : one_atom_occupations
  USE wannier_new,ONLY: use_wannier
  USE lsda_mod,  ONLY : nspin
  USE noncollin_module, ONLY : npol
  USE control_flags, ONLY: isolve, nmix, gamma_only, lscf
  USE mp_diag,   ONLY : np_ortho
  USE mp_bands,  ONLY : nproc_bgrp
  !
  IMPLICIT NONE
  !
  REAL(dp) :: omegabz
  INTEGER, PARAMETER :: Mb=1024*1024, complex_size=16, real_size=8
  INTEGER :: g_fact, nbnd_l, npwx_g, npwx_l
  !
  IF ( gamma_only) THEN
     g_fact = 2
  ELSE
     g_fact = 1
  END IF
  ! estimate npwx_l (max number of plane waves) on this processor 
  ! npwx_g is npwx_l summed over all processors
  !
  omegabz=tpi**3/omega
  npwx_g = fpi/3.0_dp * SQRT(ecutwfc)**3 / omegabz / g_fact
  npwx_l = npwx_g/nproc_bgrp
  ! 
  ! the conversions to double prevent integer overflow in very large runs
  !
  WRITE( stdout, '(/5x,"Largest allocated arrays",5x,"est. size (Mb)", &
                   &5x,"approx. dimensions")')
  WRITE( stdout, '(8x,"Kohn-Sham Wavefunctions   ",f10.2," Mb", &
                 & 5x,"(",i8,",",i5,")")') &
     complex_size*nbnd*npol*DBLE(npwx_l)/Mb, npwx_l*npol,nbnd
  IF ( one_atom_occupations .OR. use_wannier ) &
     WRITE( stdout, '(8x,"Atomic wavefunctions      ",f10.2," Mb", &
                    & 5x,"(",i8,",",i5,")")') &
   & complex_size*natomwfc*npol*DBLE(npwx_l)/Mb, npwx_l*npol,natomwfc
  IF ( lda_plus_u .AND. U_projection .NE. 'pseudo' ) &
     WRITE( stdout, '(8x,"Atomic Hubbard wavefuncts ",f10.2," Mb", &
                    & 5x,"(",i8,",",i5,")")') &
   & complex_size*nwfcU*npol*DBLE(npwx_l)/Mb, npwx_l*npol,nwfcU
  WRITE( stdout, '(8x,"NL pseudopotentials       ",f10.2," Mb", &
                 & 5x,"(",i8,",",i5,")")') &
     complex_size*nkb*DBLE(npwx_l)/Mb, npwx_l, nkb
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
  !
  WRITE( stdout, '(5x,"Largest temporary arrays",5x,"est. size (Mb)", &
                   &5x,"dimensions")')
  IF ( isolve == 0 ) THEN
     WRITE( stdout, '(8x,"Auxiliary wavefunctions   ",f10.2," Mb", &
                 & 5x,"(",i8,",",i5,")")') &
     complex_size*nbndx*npol*DBLE(npwx_l)/Mb, npwx_l*npol, nbndx
  ENDIF
  ! nbnd_l : estimated dimension of distributed matrices
  nbnd_l = nbndx/np_ortho(1)
  WRITE( stdout, '(8x,"Each subspace H/S matrix  ",f10.2," Mb", &
                 & 5x,"(",i8,",",i5,")")') &
     complex_size*nbnd_l*DBLE(nbnd_l)/g_fact/Mb, nbnd_l, nbnd_l
  !
  IF ( npol > 1 ) THEN
     WRITE( stdout, '(8x,"Each <psi_i|beta_j> matrix",f10.2," Mb", &
                 & 5x,"(",i8,",",i4,",",i5,")")') &
     complex_size*nkb*DBLE(npol*nbnd)/g_fact/Mb, nkb, npol, nbnd
  ELSE
     WRITE( stdout, '(8x,"Each <psi_i|beta_j> matrix",f10.2," Mb", &
                 & 5x,"(",i8,",",i5,")")') &
     complex_size*nkb*DBLE(nbnd)/g_fact/Mb, nkb, nbnd
  END IF
  !
  IF ( lscf) WRITE( stdout, &
     '(8x,"Arrays for rho mixing     ",f10.2," Mb", 5x,"(",i8,",",i5,")")') &
     complex_size*dfftp%nnr*DBLE(nmix)/Mb, dfftp%nnr, nmix
  !
  RETURN
  !
END subroutine memory_report
