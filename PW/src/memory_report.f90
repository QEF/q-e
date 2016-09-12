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
  ! Should be called after the first steps of initialization are done, but
  ! before large arrays are actually allocated. Does not cover cases like:
  ! real-space q/beta, vdW functionals, finite electric fields.
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : dp 
  USE constants, ONLY : tpi, fpi
  USE wvfct,     ONLY : nbnd, nbndx
  USE basis,     ONLY : natomwfc
  USE cell_base, ONLY : omega
  USE exx,       ONLY : ecutfock, nkqs
  USE fft_base,  ONLY : dffts, dfftp
  USE gvect,     ONLY : ngm, ngm_g
  USE gvecs,     ONLY : ngms, doublegrid
  USE gvecw,     ONLY : ecutwfc
  USE klist,     ONLY : nks, nkstot
  USE uspp,      ONLY : nkb, okvan
  USE funct,     ONLY : dft_is_meta, dft_is_hybrid
  USE ldaU,      ONLY : lda_plus_u, U_projection, nwfcU
  USE fixed_occ, ONLY : one_atom_occupations
  USE wannier_new,ONLY: use_wannier
  USE lsda_mod,  ONLY : nspin
  USE noncollin_module, ONLY : npol
  USE control_flags, ONLY: isolve, nmix, imix, gamma_only, lscf, io_level, &
       lxdm, smallmem
  USE ions_base, ONLY : ntyp=>nsp
  USE mp_diag,   ONLY : np_ortho
  USE mp_bands,  ONLY : nproc_bgrp, nbgrp
  USE mp_images, ONLY : nproc_image  
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: Mb=1024*1024
  INTEGER :: g_fact, mix_type_size, scf_type_size
  INTEGER :: nk, nbnd_l, npwx_g, npwx_l, ngxx_g, nexx_l
  !
  ! these quantities are real in order to prevent integer overflow
  !
  REAL(dp), PARAMETER :: complex_size=16_dp, real_size=8_dp, int_size=4_dp
  REAL(dp) :: ram, ram1, maxram, totram
  !
  IF ( gamma_only) THEN
     g_fact = 2  ! use half plane waves or G-vectors
  ELSE
     g_fact = 1  ! use  all plane waves or G-vectors
  END IF
  !
  ! npwx (max number of plane waves) is roughly estimated from the volume
  ! of the G-vector sphere, divided by the volume of the Brillouin Zone
  ! (the exact value can be computed only after G-vectors are computed)
  ! npwx_l is npwx on this processor (local)
  ! npwx_g is npwx summed over all processors (global)
  !
  npwx_g = NINT ( fpi/3.0_dp * SQRT(ecutwfc)**3 / (tpi**3/omega) / g_fact )
  npwx_l = npwx_g/nproc_bgrp
  !
  ! ram   = dynamically (permanently) allocated RAM, per process
  ! maxram= "high watermark": max ram needed during a run
  ! totram= max ram needed summed over all processors
  !
  ! Wavefunctions (including buffer space)
  IF ( io_level > 0 .OR. nks == 1) THEN
     nk = 1      ! store just one wavefunction array at the time
  ELSE
     nk = nks+1  ! store nks wavefunctions in memory (buffers)
  END IF
  ram =  complex_size * nbnd * npol * npwx_l * nk
  ! atomic wavefunctions 
  IF ( one_atom_occupations .OR. use_wannier ) &
     ram = ram + complex_size * natomwfc * npol * npwx_l 
  ! Hubbard wavefunctions
  IF ( lda_plus_u .AND. U_projection .NE. 'pseudo' ) &
     ram = ram + complex_size * nwfcU * npol * npwx_l 
  ! hybrid functionals
  IF ( dft_is_hybrid () ) THEN
     ! ngxx_g = estimated global number of G-vectors used in V_x\psi products
     ! nexx_l = estimated local size of the FFT grid used in V_x\psi products
     ngxx_g = NINT ( fpi/3.0_dp * SQRT(ecutfock)**3 / (tpi**3/omega) )
     nexx_l = 16*ngxx_g/nproc_bgrp
     ! nbnd_l : estimated number of bands per proc with band parallelization
     nbnd_l = NINT( DBLE(nbnd) / nbgrp )
     ! Stored wavefunctions in real space 
     ram = ram + complex_size/g_fact * nexx_l * npol * nbnd_l * nkqs
#if defined(__EXX_ACE)
     ! Projectors
     ram = ram + complex_size * npwx_l * npol * nbnd * nks
#endif
  END IF
  ! Nonlocal pseudopotentials V_NL (beta functions), reciprocal space
  ram =  ram + complex_size * nkb * npwx_l
  ! Charge density and potentials - see scf_type in scf_mod
  scf_type_size =  (complex_size * ngm + real_size * dfftp%nnr ) * nspin
  IF ( dft_is_meta() .or. lxdm ) scf_type_size =  2 * scf_type_size
  ! rho, v, vnew
  ram =  ram + 3 * scf_type_size
  ! vltot, vrs, rho_core, rhog_core, psic, strf, kedtau if needed
  ram =  ram + complex_size * ( dfftp%nnr + ngm *( 1 + ntyp ) ) + &
       real_size * dfftp%nnr*(2+nspin)
  IF ( dft_is_meta() ) ram = ram + real_size * dfftp%nnr*nspin
  ! arrays for rho mixing
  IF ( lscf ) THEN
     ! rhoin
     ram =  ram + scf_type_size
     ! see mix_type in scf_mod
     mix_type_size =  complex_size * ngm * nspin
     IF ( dft_is_meta() .or. lxdm ) mix_type_size =  2 * mix_type_size
     ! df, dv (if kept in memory)
     IF ( io_level < 2 ) &
          ram =  ram + mix_type_size * 2 * nmix
  END IF
  ! G-vectors: g, gg, mill, nl, nlm, ig_l2g, igtongl
  ram =  ram + real_size * ngm * 4 + int_size * ngm * 7
  ! double grid: nls, nlsm
  ram =  ram + int_size * ngms * 2
  !
  ! now scratch space that raises the "high watermark"
  !
  ! hpsi, spsi, matrices allocated in iterative diagonalization 
  ! nbnd_l : estimated dimension of distributed matrices
  !
  nbnd_l = nbndx/np_ortho(1)
  ram1 = complex_size/g_fact * ( 2*nbnd_l**2 + & ! hr, sr
                                 nkb*npol*nbnd ) ! <psi|beta>
  IF ( isolve == 0 ) THEN
     ram1 = ram1 + complex_size * nbndx * npol * npwx_l              ! hpsi
     IF ( okvan ) ram1 = ram1 + complex_size * nbndx * npol * npwx_l ! spsi
  END IF
  !
  ! arrays allocated in approx_screening2 during charge mixing 
  !
  IF ( lscf .AND. imix > 1 ) &
     ram1 = MAX( ram1, complex_size * ngm * 27 + real_size * dffts%nnr )

  maxram = ram + ram1
  !
  ! arrays used for global sorting in ggen:
  !    mill_g, mill_unsorted, igsrt, g2sort_g, total dimensions:
  !
  IF ( .NOT. smallmem ) maxram = MAX ( maxram, &
       int_size * 7 * ngm_g + real_size * ngm_g )

  totram = maxram * nproc_image
  
  WRITE( stdout, '(/5x,"Estimated max dynamical RAM per process > ", &
       & F10.2,"Mb")' ) maxram/Mb
  IF ( nproc_image > 1) WRITE( stdout, '(/5x, &
     & "Estimated total allocated dynamical RAM > ",F10.2,"Mb")' ) totram/Mb
  !
END subroutine memory_report
