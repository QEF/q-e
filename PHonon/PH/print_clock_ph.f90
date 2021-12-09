!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine print_clock_ph
  !-----------------------------------------------------------------------
  !! Print CPU and WALL time for each routine.
  !
  USE io_global,         ONLY : stdout
  USE uspp,              ONLY : okvan, nlcc_any
  USE control_ph,        ONLY : trans, zue, epsil
  USE ramanm,            ONLY : lraman, elop
  USE ldau,              ONLY : lda_plus_u
  USE el_phon,           ONLY : elph
  USE dvscf_interpolate, ONLY : ldvscf_interpolate
  USE ahc,        ONLY : elph_ahc
  USE el_phon,    ONLY : elph
  !
  implicit none
  !
  WRITE( stdout, * )
  call print_clock ('PHONON')
  WRITE( stdout,  * ) '    INITIALIZATION: '
  call print_clock ('phq_setup')
  call print_clock ('phq_init')
  WRITE( stdout, * )
  call print_clock ('phq_init')
  if (nlcc_any) call print_clock ('set_drhoc')
  call print_clock ('init_vloc')
  call print_clock ('init_us_1')
  call print_clock ('init_us_2')
  call print_clock ('newd')
  call print_clock ('dvanqq')
  call print_clock ('drho')
  if ((epsil.or.zue).and.okvan) call print_clock ('cmpt_qdipol')

  if(epsil) then
     WRITE( stdout, * )
     WRITE( stdout,  * ) '    DIELECTRIC CONSTANT AND EFFECTIVE CHARGES:'
     call print_clock ('solve_e')
     call print_clock ('dielec')
     call print_clock ('zstar_eu')
#ifdef TIMING_ZSTAR_US
     WRITE( stdout, * )
     call print_clock ('zstar_eu_us')
     call print_clock ('zstar_us_1')
     call print_clock ('zstar_us_2')
     call print_clock ('zstar_us_3')
     call print_clock ('zstar_us_4')
     call print_clock ('zstar_us_5')
#endif
#ifdef TIMING_ADD_DKMDS
     WRITE( stdout, * )
     call print_clock ('add_dkmds')
     call print_clock ('add_dkmds1')
     call print_clock ('add_dkmds2')
     call print_clock ('add_dkmds3')
     call print_clock ('add_dkmds4')
     call print_clock ('add_dkmds5')
     call print_clock ('add_dkmds6')
#endif
     if (lraman.OR.elop) then
        WRITE( stdout, * )
        WRITE( stdout,  * ) '    RAMAN COEFFICIENTS, THIRD-ORDER CHI:'
        call print_clock ('dhdrhopsi')
        if (elop) call print_clock ('el_opt')
     endif
     if (lraman) call print_clock ('dvpsi_e2')
     if (lraman) call print_clock ('solve_e2')
  endif
  if(trans) then
     WRITE( stdout, * )
     WRITE( stdout,  * ) '    DYNAMICAL MATRIX:'
     call print_clock ('dynmat0')
     call print_clock ('phqscf')
     call print_clock ('dynmatrix')
     WRITE( stdout, * )
     call print_clock ('phqscf')
     call print_clock ('solve_linter')
     call print_clock ('drhodv')
     if (zue) call print_clock('add_zstar_ue')
     if (zue) call print_clock('add_zstar_1')
     if (zue.and.okvan) call print_clock('add_zstar_us')
  endif
  WRITE( stdout, * )
  call print_clock ('dynmat0')
  call print_clock ('dynmat_us')
  call print_clock ('addusdynmat1')
  call print_clock ('d2ionq')
  if (nlcc_any) call print_clock ('dynmatcc')
  WRITE( stdout, * )
  call print_clock ('dynmat_us')
  call print_clock ('addusdynmat')
  WRITE( stdout, * )
  call print_clock ('phqscf')
  call print_clock ('solve_linter')
  WRITE( stdout, * )
  call print_clock ('solve_linter')
  call print_clock ('dvqpsi_us')
  call print_clock ('sth_kernel')
  call print_clock ('apply_dpot_b')
  call print_clock ('ortho')
  call print_clock ('cgsolve')
  call print_clock ('incdrhoscf')
  call print_clock ('addusddens')
  call print_clock ('vpsifft')
  call print_clock ('dv_of_drho')
  call print_clock ('mix_pot')
  call print_clock ('ef_shift')
  call print_clock ('ef_shift_wfc')
  call print_clock ('localdos')
#if defined(__MPI)
  call print_clock ('psymdvscf')
#else
  call print_clock ('symdvscf')
#endif
  call print_clock ('newdq')
  call print_clock ('adddvscf')


  call print_clock ('drhodvus')
  WRITE( stdout, * )
  call print_clock ('dvqpsi_us')

  call print_clock ('dvqpsi_us_on')
  WRITE( stdout, * )
  call print_clock ('cgsolve')

  call print_clock ('ch_psi')
  WRITE( stdout, * )
  call print_clock ('ch_psi')
  call print_clock ('first')
  call print_clock ('h_psi')

  call print_clock ('last')
  WRITE( stdout, * )
  call print_clock ('h_psi')
  call print_clock ('firstfft')
  call print_clock ('product')
  call print_clock ('secondfft')

  call print_clock ('add_vuspsi')
  WRITE( stdout, * )
  call print_clock ('incdrhoscf')

  call print_clock ('addusdbec')
  WRITE( stdout, * )
  call print_clock ('drhodvus')

  call print_clock ('addusddort')
  WRITE( stdout, * )
  
  IF (lda_plus_u) THEN
     call print_clock ('dnsq_bare')
     call print_clock ('dwfc')
     call print_clock ('swfc')
     call print_clock ('delta_sphi')
     call print_clock ('dnsq_orth')
     call print_clock ('d2nsq_bare_k')
     call print_clock ('dnsq_scf')
     call print_clock ('adddvhubscf')
     call print_clock ('dynmat_hub_bare')
     call print_clock ('dvqhub_barepsi_us')
     call print_clock ('dvqhub_barepsi_us2')
     call print_clock ('dynmat_hub_scf')
     call print_clock ('doubleprojqq')
     call print_clock ('doubleprojqq2')
     WRITE( stdout, * )
  ENDIF

  IF (elph) THEN
     WRITE( stdout, * ) '     Electron-phonon coupling'
     CALL print_clock('elphon')
     CALL print_clock('elphel')
     CALL print_clock('elphsum')
     WRITE( stdout, * )
  ENDIF

  IF (ldvscf_interpolate) THEN
     WRITE( stdout, * ) '     Fourier interpolation of dVscf'
     ! CALL print_clock('dvscf_setup')
     CALL print_clock('dvscf_r2q')
     CALL print_clock('dvscf_davcio')
     CALL print_clock('dvscf_scatgrid')
     ! CALL print_clock('dvscf_fourier')
     ! CALL print_clock('dvscf_mpsum')
     ! CALL print_clock('dvscf_shift')
     ! CALL print_clock('dvscf_iqr')
     ! CALL print_clock('dvscf_long')
     CALL print_clock('dvscf_bare')
     ! CALL print_clock('dvscf_cart2u')
     WRITE( stdout, * )
  ENDIF

  IF (elph_ahc) THEN
     WRITE( stdout, * ) '     El-ph coupling for electron self-energy'
     CALL print_clock('ahc_elph')
     CALL print_clock('ahc_upfan')
     CALL print_clock('ahc_dw')
     CALL print_clock('ahc_gkk')
     CALL print_clock('ahc_gauge')
     WRITE( stdout, * )
  ENDIF

  WRITE( stdout,  * ) '     General routines'
  call print_clock ('calbec')
  call print_clock ('fft')
  call print_clock ('ffts')
  call print_clock ('fftw')
  call print_clock ('cinterpolate')
  call print_clock ('davcio')
  call print_clock ('write_rec')
  WRITE( stdout, * )
  return
end subroutine print_clock_ph
