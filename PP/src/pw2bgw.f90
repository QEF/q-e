!
! Copyright (C) 2008-2012 Georgy Samsonidze
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Converts the output files produced by pw.x to the input files for BerkeleyGW.
! Last modified 2/27/2024 by Fangzhou Zhao:
! Full spinor with magnetization support for pw2bgw in QE7.2
! To use full spinor + magnetization, no symmetries are allowed
! disabled npool use in write_vkbg
!
! Usage: Copy pw2bgw_qe6.4_with_spinor.f90 to QE6.4/PP/src/pw2bgw.f90 and then make
! This version supports conversion for both non-spinor and spinor cases
! For full-spinor, only VXC is supported (vxc.dat is not yet)
!
!-------------------------------------------------------------------------------
!
! BerkeleyGW, Copyright (c) 2011, The Regents of the University of
! California, through Lawrence Berkeley National Laboratory (subject to
! receipt of any required approvals from the U.S. Dept. of Energy).
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
! (1) Redistributions of source code must retain the above copyright
! notice, this list of conditions and the following disclaimer.
!
! (2) Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! (3) Neither the name of the University of California, Lawrence
! Berkeley National Laboratory, U.S. Dept. of Energy nor the names of
! its contributors may be used to endorse or promote products derived
! from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! You are under no obligation whatsoever to provide any bug fixes,
! patches, or upgrades to the features, functionality or performance of
! the source code ("Enhancements") to anyone; however, if you choose to
! make your Enhancements available either publicly, or directly to
! Lawrence Berkeley National Laboratory, without imposing a separate
! written license agreement for such Enhancements, then you hereby grant
! the following license: a  non-exclusive, royalty-free perpetual
! license to install, use, modify, prepare derivative works, incorporate
! into other computer software, distribute, and sublicense such
! enhancements or derivative works thereof, in binary and source code
! form.
!
!-------------------------------------------------------------------------------
!
! pw2bgw subroutines:
!
! write_wfng  - generates complex wavefunctions in G-space (normalized to 1)
! real_wfng   - constructs real wavefunctions by applying the Gram-Schmidt
!               process (called from write_wfng)
! write_rhog  - generates real/complex charge density in G-space
!               (units of the number of electronic states per unit cell)
! calc_rhog   - computes charge density by summing over a subset of occupied
!               bands (called from write_rhog), destroys charge density
! write_vxcg  - generates real/complex exchange-correlation potential in
!               G-space (units of Rydberg) [only local part of Vxc]
! write_vxc0  - prints real/complex exchange-correlation potential at G=0
!               (units of eV) [only local part of Vxc]
! write_vxc_r - calculates matrix elements of exchange-correlation potential
!               in R-space (units of eV) [only local part of Vxc]
! write_vxc_g - calculates matrix elements of exchange-correlation potential
!               in G-space (units of eV) [only local part of Vxc]
! write_kih   - calculates matrix elements of KIH = Kinetic energy + Ionic potential
!               + Hartree (units of eV) [supports non-local Vxc] [supports metaGGA, 
!               hybrid functionals]  
! write_vhub_g - calculates matrix elements of Hubbard potential
! write_vxc_hybrid - calculates matrix elements of E_nk - KIH (which is the exact
!               exchange part + Vxc] (units of eV) [supports non-local Vxc]
!               [supports nspin = 2] !FZ: only used for hybrid functional calculations
!   [FZ: for hybrid functionals: Both write_kih and write_vxc_hybrid can be used,
!    write_kih is recommented. For non-hybrid functionals: write_kih can be used, should
!    not use write_vxc_hybrid. ]
! write_vscg  - generates real/complex self-consistent potential in G-space
!               (units of Rydberg) [only local part of Vsc]
! write_vkbg  - generates complex Kleinman-Bylander projectors in G-space
!               (units of Rydberg)
! check_inversion - checks whether real/complex version is appropriate
!               (called from everywhere)
!
! Quantum ESPRESSO stores the wavefunctions in is-ik-ib-ig order
! BerkeleyGW stores the wavefunctions in ik-ib-is-ig order
! the outer loop is over is(QE)/ik(BGW) and the inner loop is over ig
! ik = k-point index, is = spin index, ib = band index, ig = G-vector index
!
! write_wfng reverts the order of is and ik using smap and kmap arrays,
! distributes wavefunctions over processors by ig (either real case or
! spin-polarized case), calls real_wfng that applies the Gram-Schmidt
! process (real case), reverts the order of is and ib (spin-polarized
! case), and writes wavefunctions to disk
!
!-------------------------------------------------------------------------------

PROGRAM pw2bgw

  USE constants, ONLY : eps12
  USE control_flags, ONLY : gamma_only
  USE environment, ONLY : environment_start, environment_end
  USE io_files, ONLY : prefix, tmp_dir
  USE io_global, ONLY : ionode, ionode_id
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin, starting_magnetization !FZ: added starting_magnetization for spinors
  USE mp, ONLY : mp_bcast
  USE mp_world, ONLY : world_comm
  USE mp_global, ONLY : mp_startup
  USE paw_variables, ONLY : okpaw
  USE scf, ONLY : rho_core, rhog_core
  USE uspp, ONLY : okvan
  !USE funct, ONLY : dft_is_hybrid, dft_is_meta, dft_is_gradient !FZ: for hybrid functional, metaGGA, add dft_is_gradient for spinors
  USE xc_lib, ONLY : xclib_dft_is  !FZ: for qe6.7, for hybrid functional, metaGGA, add dft_is_gradient for spinors
  USE ldaU, ONLY : lda_plus_u
  USE scf,   ONLY : rho, rho_core, rhog_core, v, vrs!FZ: test spinors added vrs 
  USE fft_rho,  ONLY : rho_g2r  
  USE fft_base,  ONLY : dfftp    
  USE noncollin_module,   ONLY : noncolin, npol, m_loc, ux, angle1, angle2  !FZ: test for spinors
  USE ions_base,          ONLY : nat, tau, ntyp => nsp, ityp, zv  !FZ: test for spinors

  IMPLICIT NONE

  character(len=6) :: codename = 'PW2BGW'

  integer :: real_or_complex
  character ( len = 9 ) :: symm_type
  logical :: wfng_flag
  character ( len = 256 ) :: wfng_file
  logical :: wfng_kgrid
  integer :: wfng_nk1
  integer :: wfng_nk2
  integer :: wfng_nk3
  real (DP) :: wfng_dk1
  real (DP) :: wfng_dk2
  real (DP) :: wfng_dk3
  logical :: wfng_occupation
  integer :: wfng_nvmin
  integer :: wfng_nvmax
  logical :: rhog_flag
  character ( len = 256 ) :: rhog_file
  integer :: rhog_nvmin
  integer :: rhog_nvmax
  logical :: vxcg_flag
  character ( len = 256 ) :: vxcg_file
  logical :: vxc0_flag
  character ( len = 256 ) :: vxc0_file
  logical :: vxc_flag
  character ( len = 256 ) :: vxc_file
  character :: vxc_integral
  logical :: kih_flag                 !FZ: Kinetic energy + Ionic potential + Hartree 
  character ( len = 256 ) :: kih_file !FZ: 
  logical :: vxc_hybrid_flag                 !FZ: Only used for hybrid functional!
  character ( len = 256 ) :: vxc_hybrid_file !FZ: Band Energy - (Kinetic energy + Ionic potential + Hartree) 
  integer :: vxc_diag_nmin
  integer :: vxc_diag_nmax
  integer :: vxc_offdiag_nmin
  integer :: vxc_offdiag_nmax
  integer :: kih_diag_nmin        !FZ: for kih, setting kih_*_nm* = vxc_*_nm* 
  integer :: kih_diag_nmax        !FZ: for kih 
  integer :: kih_offdiag_nmin     !FZ: for kih
  integer :: kih_offdiag_nmax     !FZ: for kih
  !integer :: vxc_hybrid_diag_nmin        !FZ: Only for hybrid functional
  !integer :: vxc_hybrid_diag_nmax        !FZ: 
  !integer :: vxc_hybrid_offdiag_nmin     !FZ: 
  !integer :: vxc_hybrid_offdiag_nmax     !FZ: 
  logical :: vxc_zero_rho_core
  logical :: vscg_flag
  character ( len = 256 ) :: vscg_file
  logical :: vkbg_flag
  character ( len = 256 ) :: vkbg_file
  character ( len = 256 ) :: outdir
  logical :: vhub_flag  
  character ( len = 256 ) :: vhub_file
  integer :: vhub_diag_nmin, vhub_diag_nmax, vhub_offdiag_nmin, vhub_offdiag_nmax

  NAMELIST / input_pw2bgw / prefix, outdir, &
    real_or_complex, symm_type, wfng_flag, wfng_file, wfng_kgrid, &
    wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, wfng_dk3, &
    wfng_occupation, wfng_nvmin, wfng_nvmax, rhog_flag, rhog_file, &
    rhog_nvmin, rhog_nvmax, vxcg_flag, vxcg_file, vxc0_flag, vxc0_file, &
    vxc_flag, vxc_file, vxc_integral, vxc_diag_nmin, vxc_diag_nmax, &
    vxc_offdiag_nmin, vxc_offdiag_nmax, vxc_zero_rho_core, &
    vscg_flag, vscg_file, vkbg_flag, vkbg_file, kih_flag, & !FZ: for KIH
    kih_file, vxc_hybrid_flag, vxc_hybrid_file, & !FZ: for hybrid functional
    vhub_flag, vhub_file, vhub_diag_nmin, vhub_diag_nmax, vhub_offdiag_nmin, vhub_offdiag_nmax ! MW: for Hubbard potential
    !kih_diag_nmin, kih_diag_nmax, & !FZ: for KIH 
    !kih_offdiag_nmin, kih_offdiag_nmax, vxc_hybrid_diag_nmin, &  !FZ: for hybrid
    !vxc_hybrid_diag_nmax, vxc_hybrid_offdiag_nmin, vxc_hybrid_offdiag_nmax  !FZ: for hybrid functional
    !FZ: these variables are currently commented out because

  integer :: ii, ios
  character ( len = 256 ) :: output_file_name
  character ( len = 256 ) :: output_file_name2  !FZ:
  logical :: allow_spinors

  character (len=256), external :: trimcheck
  character (len=1), external :: lowercase
  REAL(DP) :: etxc, vtxc  !FZ: change062320
  INTEGER  :: na   !FZ: test for spinors

#if defined(__MPI)
  CALL mp_startup ( )
#endif

  CALL environment_start ( codename )

  prefix = 'prefix'
  CALL get_environment_variable ( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM ( outdir ) == ' ' ) outdir = './'
  real_or_complex = 2
  symm_type = 'cubic'
  wfng_flag = .FALSE.
  wfng_file = 'WFN'
  wfng_kgrid = .FALSE.
  wfng_nk1 = 0
  wfng_nk2 = 0
  wfng_nk3 = 0
  wfng_dk1 = 0.0D0
  wfng_dk2 = 0.0D0
  wfng_dk3 = 0.0D0
  wfng_occupation = .FALSE.
  wfng_nvmin = 0
  wfng_nvmax = 0
  rhog_flag = .FALSE.
  rhog_file = 'RHO'
  rhog_nvmin = 0
  rhog_nvmax = 0
  vxcg_flag = .FALSE.
  vxcg_file = 'VXC'
  vxc0_flag = .FALSE.
  vxc0_file = 'vxc0.dat'
  vxc_flag = .FALSE.
  vxc_file = 'vxc.dat'
  kih_flag = .FALSE.    !FZ: output kih energy
  kih_file = 'kih.dat'  !FZ:
  vxc_hybrid_flag = .FALSE.    !FZ: only used for hybrid functionals
  vxc_hybrid_file = 'vxc_hybrid.dat'  !FZ: 
  vxc_integral = 'g'
  vxc_diag_nmin = 0
  vxc_diag_nmax = 0
  vxc_offdiag_nmin = 0
  vxc_offdiag_nmax = 0
  vxc_zero_rho_core = .TRUE.
  vscg_flag = .FALSE.
  vscg_file = 'VSC'
  vkbg_flag = .FALSE.
  vkbg_file = 'VKB'
  vhub_flag = .false.  
  vhub_file = 'vhub.dat'
  vhub_diag_nmin = 0
  vhub_diag_nmax = 0
  vhub_offdiag_nmin = 0
  vhub_offdiag_nmax = 0

  IF ( ionode ) THEN
    CALL input_from_file ( )
    READ ( 5, input_pw2bgw, iostat = ios )
    IF ( ios /= 0 ) CALL errore ( codename, 'input_pw2bgw', abs ( ios ) )

    DO ii = 1, LEN_TRIM (symm_type)
      symm_type(ii:ii) = lowercase (symm_type(ii:ii))
    END DO
    DO ii = 1, LEN_TRIM (vxc_integral)
      vxc_integral(ii:ii) = lowercase (vxc_integral(ii:ii))
    END DO

    IF ( real_or_complex /= 1 .AND. real_or_complex /= 2 ) &
      CALL errore ( codename, 'real_or_complex', 1 )
    IF ( symm_type /= 'cubic' .AND. symm_type /= 'hexagonal' ) &
      CALL errore ( codename, 'symm_type', 1 )
    IF ( vxc_integral /= 'r' .AND. vxc_integral /= 'g' ) &
      CALL errore ( codename, 'vxc_integral', 1 )
    kih_diag_nmin = vxc_diag_nmin                 !FZ: For usage of KIH energy
    kih_diag_nmax = vxc_diag_nmax                 !FZ: temporarily setting all the
    kih_offdiag_nmin = vxc_offdiag_nmin           !FZ: number of diagonal and offdiag
    kih_offdiag_nmax = vxc_offdiag_nmax           !FZ: mtxel numbers equal to that
    !vxc_hybrid_diag_nmin = vxc_diag_nmin          !FZ: of vxc
    !vxc_hybrid_diag_nmax = vxc_diag_nmax          !FZ: 
    !vxc_hybrid_offdiag_nmin = vxc_offdiag_nmin    !FZ: 
    !vxc_hybrid_offdiag_nmax = vxc_offdiag_nmax    !FZ:
  ENDIF

  tmp_dir = trimcheck ( outdir )
  CALL mp_bcast ( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast ( prefix, ionode_id, world_comm )
  CALL mp_bcast ( real_or_complex, ionode_id, world_comm )
  CALL mp_bcast ( symm_type, ionode_id, world_comm )
  CALL mp_bcast ( wfng_flag, ionode_id, world_comm )
  CALL mp_bcast ( wfng_file, ionode_id, world_comm )
  CALL mp_bcast ( wfng_kgrid, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nk1, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nk2, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nk3, ionode_id, world_comm )
  CALL mp_bcast ( wfng_dk1, ionode_id, world_comm )
  CALL mp_bcast ( wfng_dk2, ionode_id, world_comm )
  CALL mp_bcast ( wfng_dk3, ionode_id, world_comm )
  CALL mp_bcast ( wfng_occupation, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nvmin, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nvmax, ionode_id, world_comm )
  CALL mp_bcast ( rhog_flag, ionode_id, world_comm )
  CALL mp_bcast ( rhog_file, ionode_id, world_comm )
  CALL mp_bcast ( rhog_nvmin, ionode_id, world_comm )
  CALL mp_bcast ( rhog_nvmax, ionode_id, world_comm )
  CALL mp_bcast ( vxcg_flag, ionode_id, world_comm )
  CALL mp_bcast ( vxcg_file, ionode_id, world_comm )
  CALL mp_bcast ( vxc0_flag, ionode_id, world_comm )
  CALL mp_bcast ( vxc0_file, ionode_id, world_comm )
  CALL mp_bcast ( vxc_flag, ionode_id, world_comm )
  CALL mp_bcast ( vxc_integral, ionode_id, world_comm )
  CALL mp_bcast ( vxc_file, ionode_id, world_comm )
  CALL mp_bcast ( kih_flag, ionode_id, world_comm )     !FZ: kih energy
  CALL mp_bcast ( kih_file, ionode_id, world_comm )     !FZ: 
  CALL mp_bcast ( vxc_hybrid_flag, ionode_id, world_comm )     !FZ: Only used for hybrid functionals
  CALL mp_bcast ( vxc_hybrid_file, ionode_id, world_comm )     !FZ: 
  CALL mp_bcast ( kih_diag_nmin, ionode_id, world_comm )     !FZ: 
  CALL mp_bcast ( kih_diag_nmax, ionode_id, world_comm )     !FZ: 
  CALL mp_bcast ( kih_offdiag_nmin, ionode_id, world_comm )     !FZ: 
  CALL mp_bcast ( kih_offdiag_nmax, ionode_id, world_comm )     !FZ: 
  !CALL mp_bcast ( vxc_hybrid_diag_nmin, ionode_id, world_comm )     !FZ: 
  !CALL mp_bcast ( vxc_hybrid_diag_nmax, ionode_id, world_comm )     !FZ: 
  !CALL mp_bcast ( vxc_hybrid_offdiag_nmin, ionode_id, world_comm )     !FZ: 
  !CALL mp_bcast ( vxc_hybrid_offdiag_nmax, ionode_id, world_comm )     !FZ: 
  CALL mp_bcast ( vxc_diag_nmin, ionode_id, world_comm )
  CALL mp_bcast ( vxc_diag_nmax, ionode_id, world_comm )
  CALL mp_bcast ( vxc_offdiag_nmin, ionode_id, world_comm )
  CALL mp_bcast ( vxc_offdiag_nmax, ionode_id, world_comm )
  CALL mp_bcast ( vxc_zero_rho_core, ionode_id, world_comm )
  CALL mp_bcast ( vscg_flag, ionode_id, world_comm )
  CALL mp_bcast ( vscg_file, ionode_id, world_comm )
  CALL mp_bcast ( vkbg_flag, ionode_id, world_comm )
  CALL mp_bcast ( vkbg_file, ionode_id, world_comm )
  CALL mp_bcast ( vhub_flag, ionode_id, world_comm )  
  CALL mp_bcast ( vhub_file, ionode_id, world_comm )
  CALL mp_bcast ( vhub_diag_nmin, ionode_id, world_comm )
  CALL mp_bcast ( vhub_diag_nmax, ionode_id, world_comm )
  CALL mp_bcast ( vhub_offdiag_nmin, ionode_id, world_comm )
  CALL mp_bcast ( vhub_offdiag_nmax, ionode_id, world_comm )
  CALL read_file ( )

  !CALL setup ()
  ALLOCATE( m_loc( 3, nat ) )
  IF ( noncolin ) THEN
     DO na = 1, nat
        !
        m_loc(1,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * COS( angle2(ityp(na)) )
        m_loc(2,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * SIN( angle2(ityp(na)) )
        m_loc(3,na) = starting_magnetization(ityp(na)) * &
                      COS( angle1(ityp(na)) )
     END DO
     !FZ:  initialize the quantization direction for gga
     ux=0.0_DP
     !if (dft_is_gradient()) call compute_ux(m_loc,ux,nat)
     if (xclib_dft_is('gradient')) call compute_ux(m_loc,ux,nat)  !FZ: for qe6.7
  ENDIF

  if (ionode) then
    if (MAX (MAXVAL (ABS (rho_core (:) ) ), MAXVAL (ABS (rhog_core (:) ) ) ) &
      .LT. eps12) then
      WRITE ( 6, '(/,5x,"NLCC is absent")' )
    else
      WRITE ( 6, '(/,5x,"NLCC is present")' )
    endif
  endif
  if (okvan) call errore ( 'pw2bgw', 'BGW cannot use USPP.', 3 )
  if (okpaw) call errore ( 'pw2bgw', 'BGW cannot use PAW.', 4 )
  if (gamma_only) call errore ( 'pw2bgw', 'BGW cannot use gamma-only run.', 5 )
  ! ZL: change to a if statement for internal version test
  !if (nspin == 4) call errore ( 'pw2bgw', 'BGW cannot use spinors.', 6 )
  allow_spinors = .true.
  if (nspin == 4) then
    if (.not. allow_spinors) then
      call errore ( 'pw2bgw', 'BGW cannot use spinors.', 6 )
    else
      WRITE ( 6, '(/,5x,"pw2bgw running with full spinors (non-collinear).")' ) 
    endif
  endif
  ! ZL: end modify

  ! FZ: changes for hybrid functional support
  !if (dft_is_hybrid()) then
  if (xclib_dft_is('hybrid')) then  !FZ: for qe6.7
    !FZ: commented 01/21
    !if (nspin == 4) then
    !  call errore ( 'pw2bgw', 'pw2bgw current do not support hybrid functional in spinor calculations.', 9 )
    !endif
    if ((vxc_flag .or. vxcg_flag) .or. vxc0_flag) then
      vxc_flag = .FALSE.
      vxcg_flag = .FALSE.
      vxc0_flag = .FALSE.
      WRITE ( 6, '(/,5x," setting vxc_flag, vxcg_flag, vxc0_flag = .false. automatically")' ) !FZ: if using hybrid functional, 
      call errore ( 'pw2bgw', ' hybrid functionals should and should only use kih.dat or vxc_hybrid.dat.', 10 ) !FZ: if using hybrid functional,
      !FZ: setting vxc_flag, vxcg_flag, vxc0_flag = .false. automatically, should not use the normal way to calculate vxc
      if (.not. (kih_flag .or. vxc_hybrid_flag)) then
        WRITE ( 6, '(/,5x," You are using hybrid functional, suggest setting kih_flag or vxc_hybrid_flag = .true.")' )
        call errore ( 'pw2bgw', ' hybrid functionals should and should only use kih.dat or vxc_hybrid.dat.', 10 )
      endif
    endif
  else
    if (vxc_hybrid_flag) then
      call errore ( 'pw2bgw', 'vxc_hybrid should only be used when using hybrid functionals.', 7 )
    endif
  endif
  !FZ: 1220
  !if (nspin == 4 .and. (((vxc_flag .or. vxc0_flag) .or. kih_flag) .or. vxc_hybrid_flag)) then
  !  call errore ( 'pw2bgw', 'pw2bgw cannot support spinors when calculating vxc.dat '// &
  !  'vxc0.dat, kih.dat and vxc_hybrid.dat. Please use VXC.', 8 )
  !endif
  !if (dft_is_meta() .and. nspin == 4) then
  if (xclib_dft_is('meta') .and. nspin == 4) then   !FZ: for qe6.7
    call errore ( 'pw2bgw', 'pw2bgw current do not support metaGGA in spinor calculations.', 10 )
  endif

  ! MW: open access to prefix.hub* (ready to write)
  CALL openfil()
  ! MW: hinit0() contains gk_sort(); open access to ./prefix.wfc* files for reading
  ! set nwordwfc = nbnd * npwx * npol here
  CALL hinit0()

  !IF (dft_is_meta()) then                           !FZ:  for metaGGA change062320
  IF (xclib_dft_is('meta')) then                           !FZ:  for qe6.7
     CALL rho_g2r ( dfftp, rho%kin_g, rho%kin_r )   !FZ:  for metaGGA
     CALL v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, v%of_r, v%kin_r )
  ENDIF
  ! FZ: end modify

  if (real_or_complex == 1 .AND. vxc_flag .AND. vxc_offdiag_nmax > 0) &
    call errore ( 'pw2bgw', 'Off-diagonal matrix elements of Vxc ' // &
    'with real wavefunctions are not implemented, compute them in ' // &
    'Sigma using VXC.', 7)

  CALL openfil_pp ( )

  ! MW: Atomic configuration dependent hamiltonian initialization
  IF ( lda_plus_u ) THEN
    IF (ionode) THEN
      write(*,*) "lda_plus_u = T, call orthoUwfc()"
    ENDIF
    CALL orthoUwfc()
  ENDIF
  
  if ( ionode ) WRITE ( 6, '("")' )

  IF ( wfng_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( wfng_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_wfng")' )
    CALL start_clock ( 'write_wfng' )
    CALL write_wfng ( output_file_name, real_or_complex, symm_type, &
      wfng_kgrid, wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, &
      wfng_dk3, wfng_occupation, wfng_nvmin, wfng_nvmax )
    CALL stop_clock ( 'write_wfng' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_wfng",/)' )
  ENDIF

  IF ( vxcg_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( vxcg_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_vxcg")' )
    CALL start_clock ( 'write_vxcg' )
    !IF (dft_is_meta()) THEN  !FZ: for metaGGA
    IF (xclib_dft_is('meta')) THEN  !FZ: for qe6.7
      CALL write_vxcg_meta ( output_file_name, real_or_complex, symm_type, &
        vxc_zero_rho_core )
    ELSE
      CALL write_vxcg ( output_file_name, real_or_complex, symm_type, &
        vxc_zero_rho_core )
    ENDIF   !FZ: end modify
    CALL stop_clock ( 'write_vxcg' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_vxcg",/)' )
  ENDIF

  IF ( vxc0_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( vxc0_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_vxc0")' )
    CALL start_clock ( 'write_vxc0' )
    !IF (dft_is_meta()) THEN  !FZ: for metaGGA
    IF (xclib_dft_is('meta')) THEN  !FZ: for qe6.7
      CALL write_vxc0_meta ( output_file_name, vxc_zero_rho_core )
    ELSE
      CALL write_vxc0 ( output_file_name, vxc_zero_rho_core )
    ENDIF  !FZ: end modify
    CALL stop_clock ( 'write_vxc0' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_vxc0",/)' )
  ENDIF

  IF ( vxc_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( vxc_file )
    IF ( vxc_integral .EQ. 'r' ) THEN
      IF ( ionode ) WRITE ( 6, '(5x,"call write_vxc_r")' )
      CALL start_clock ( 'write_vxc_r' )
      !IF (dft_is_meta()) THEN  !FZ: for metaGGA
      IF (xclib_dft_is('meta')) THEN  !FZ: for qe6.7
        CALL write_vxc_r_meta ( output_file_name, &
          vxc_diag_nmin, vxc_diag_nmax, &
          vxc_offdiag_nmin, vxc_offdiag_nmax, &
          vxc_zero_rho_core )
      ELSE
        CALL write_vxc_r ( output_file_name, &
          vxc_diag_nmin, vxc_diag_nmax, &
          vxc_offdiag_nmin, vxc_offdiag_nmax, &
          vxc_zero_rho_core )
      ENDIF  !FZ: end modify
      CALL stop_clock ( 'write_vxc_r' )
      IF ( ionode ) WRITE ( 6, '(5x,"done write_vxc_r",/)' )
    ENDIF
    IF ( vxc_integral .EQ. 'g' ) THEN
      IF ( ionode ) WRITE ( 6, '(5x,"call write_vxc_g")' )
      CALL start_clock ( 'write_vxc_g' )
      !IF (dft_is_meta()) THEN  !FZ: for metaGGA
      IF (xclib_dft_is('meta')) THEN  !FZ: for qe6.7
        CALL write_vxc_g_meta ( output_file_name, &
          vxc_diag_nmin, vxc_diag_nmax, &
          vxc_offdiag_nmin, vxc_offdiag_nmax, &
          vxc_zero_rho_core )
      ELSE
        CALL write_vxc_g ( output_file_name, &
          vxc_diag_nmin, vxc_diag_nmax, &
          vxc_offdiag_nmin, vxc_offdiag_nmax, &
          vxc_zero_rho_core )
      ENDIF  !FZ: end modify
      CALL stop_clock ( 'write_vxc_g' )
      IF ( ionode ) WRITE ( 6, '(5x,"done write_vxc_g",/)' )
    ENDIF
  ENDIF

  !FZ:   This block is for kih energy
  IF ( kih_flag .or. vxc_hybrid_flag) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( kih_file )
    output_file_name2 = TRIM ( tmp_dir ) // TRIM ( vxc_hybrid_file )
      IF ( ionode ) WRITE ( 6, '(5x,"call write_kih")' )
      IF ( ionode ) WRITE ( 6, '(5x,A)') &
           "If you use pw2bgw to output matrix elements of the KIH energy, please cite " 
           WRITE ( 6, '(5x,A)')  & 
           "F. Zhao and S. G. Louie, U. C. Berkeley PhD dissertation (2021)."    
      CALL start_clock ( 'write_kih' )
      CALL write_kih ( output_file_name, output_file_name2, &
        kih_diag_nmin, kih_diag_nmax, &
        kih_offdiag_nmin, kih_offdiag_nmax, &
        kih_flag, vxc_hybrid_flag)
      CALL stop_clock ( 'write_kih' )
      IF ( ionode ) WRITE ( 6, '(5x,"done write_kih",/)' )
  ENDIF

  IF ( vscg_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( vscg_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_vscg")' )
    CALL start_clock ( 'write_vscg' )
    CALL write_vscg ( output_file_name, real_or_complex, symm_type )
    CALL stop_clock ( 'write_vscg' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_vscg",/)' )
  ENDIF

  IF ( vkbg_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( vkbg_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_vkbg")' )
    CALL start_clock ( 'write_vkbg' )
    CALL write_vkbg ( output_file_name, symm_type, wfng_kgrid, wfng_nk1, &
      wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, wfng_dk3 )
    CALL stop_clock ( 'write_vkbg' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_vkbg",/)' )
  ENDIF

  ! since calc_rhog (called from write_rhog) destroys charge density,
  ! it must be called after v_xc (called from write_vxcg, write_vxc0,
  ! write_vxc_r, write_vxc_g)
  IF ( rhog_flag ) THEN
    output_file_name = TRIM ( tmp_dir ) // TRIM ( rhog_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_rhog")' )
    CALL start_clock ( 'write_rhog' )
    CALL write_rhog ( output_file_name, real_or_complex, symm_type, &
      rhog_nvmin, rhog_nvmax )
    CALL stop_clock ( 'write_rhog' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_rhog",/)' )
  ENDIF

  IF (lda_plus_u .and. vhub_flag) THEN
    output_file_name = TRIM ( outdir ) // '/' // TRIM ( vhub_file )
    IF ( ionode ) WRITE ( 6, '(5x,"call write_vhub")' )
    IF ( ionode ) WRITE ( 6, '(5x,A)') &
         "If you use pw2bgw to output matrix elements of the Hubbard potential, please cite Nat. Commun. 10, 2371 (2019)."
    CALL start_clock ( 'write_vhub' )
    CALL write_vhub_g ( output_file_name, vhub_diag_nmin, vhub_diag_nmax, vhub_offdiag_nmin, vhub_offdiag_nmax)
    CALL stop_clock ( 'write_vhub' )
    IF ( ionode ) WRITE ( 6, '(5x,"done write_vhub",/)' )
  ENDIF
  
  IF ( ionode ) WRITE ( 6, * )
  IF ( wfng_flag ) CALL print_clock ( 'write_wfng' )
  IF ( rhog_flag ) CALL print_clock ( 'write_rhog' )
  IF ( vxcg_flag ) CALL print_clock ( 'write_vxcg' )
  IF ( vxc0_flag ) CALL print_clock ( 'write_vxc0' )
  IF ( vxc_flag ) THEN
    IF ( vxc_integral .EQ. 'r' ) CALL print_clock ( 'write_vxc_r' )
    IF ( vxc_integral .EQ. 'g' ) CALL print_clock ( 'write_vxc_g' )
  ENDIF
  IF ( vscg_flag ) CALL print_clock ( 'write_vscg' )
  IF ( vkbg_flag ) CALL print_clock ( 'write_vkbg' )
  IF ( wfng_flag .AND. real_or_complex .EQ. 1 ) THEN
    IF ( ionode ) WRITE ( 6, '(/,5x,"Called by write_wfng:")' )
    CALL print_clock ( 'real_wfng' )
  ENDIF

  CALL environment_end ( codename )

  CALL stop_pp ( )

  ! this is needed because openfil is called above
  CALL close_files ( .false. )

  STOP

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE write_wfng ( output_file_name, real_or_complex, symm_type, &
  wfng_kgrid, wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, &
  wfng_dk3, wfng_occupation, wfng_nvmin, wfng_nvmax )

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  USE constants, ONLY : pi, tpi, eps6
  USE fft_base, ONLY : dfftp
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, g, mill, ecutrho
  USE io_files, ONLY : iunwfc, nwordwfc
  USE io_global, ONLY : ionode, ionode_id
  USE ions_base, ONLY : nat, atm, ityp, tau
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, wk, ngk, nks, nkstot, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum, mp_max, mp_get, mp_bcast, mp_barrier
  USE mp_pools, ONLY : me_pool, root_pool, npool, nproc_pool, &
    intra_pool_comm, inter_pool_comm
  USE mp_wave, ONLY : mergewf
  USE mp_world, ONLY : mpime, nproc, world_comm
  USE mp_bands, ONLY : intra_bgrp_comm, nbgrp
  USE start_k, ONLY : nk1, nk2, nk3, k1, k2, k3
  USE symm_base, ONLY : s, ft, nsym
  USE wavefunctions, ONLY : evc
  USE wvfct, ONLY : npwx, nbnd, et, wg
  USE gvecw, ONLY : ecutwfc
  USE matrix_inversion
#if defined(__MPI)
  USE parallel_include, ONLY : MPI_DOUBLE_COMPLEX
#endif

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  integer, intent (in) :: real_or_complex
  character ( len = 9 ), intent (in) :: symm_type
  logical, intent (in) :: wfng_kgrid
  integer, intent (in) :: wfng_nk1
  integer, intent (in) :: wfng_nk2
  integer, intent (in) :: wfng_nk3
  real (DP), intent (in) :: wfng_dk1
  real (DP), intent (in) :: wfng_dk2
  real (DP), intent (in) :: wfng_dk3
  logical, intent (in) :: wfng_occupation
  integer, intent (in) :: wfng_nvmin
  integer, intent (in) :: wfng_nvmax

  character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  logical :: proc_wf, bad_kgrid
  integer :: unit, i, j, k, cell_symmetry, nrecord
  integer :: id, ib, ik, iks, ike, is, isp, ig, ierr  ! ZL: added isp
  integer :: nd, ntran, nb, nk_l, nk_g, ns, nst, nsf, ng_l, ng_g, npw_g2  ! ZL: added nst, nsf, npw_g2
  integer :: npw, ngg, npw_g, npwx_g
  integer :: local_pw, ipsour, igwx, ngkdist_g, ngkdist_l
  real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: kmap ( : )
  integer, allocatable :: smap ( : )
  integer, allocatable :: ifmin ( : )
  integer, allocatable :: ifmax ( : )
  integer, allocatable :: itmp ( : )
  integer, allocatable :: ngk_g ( : )
  integer, allocatable :: ipmask ( : )
  integer, allocatable :: igwk ( : )
  integer, allocatable :: igwf_l2g ( : )
  integer, allocatable :: g_g ( :, : )
  integer, allocatable :: igk_l2g ( :, : )
  real (DP), allocatable :: et_g ( :, : )
  real (DP), allocatable :: wg_g ( :, : )
  real (DP), allocatable :: energy ( :, : )
  complex (DP), allocatable :: wfng ( : )
  ! ZL: full-spinor two-component wavefunctions, for non-collinear spin case
  complex (DP), allocatable :: wfng_nc ( :, : )
  complex (DP), allocatable :: wfng_buf ( :, : )
  complex (DP), allocatable :: wfng_dist ( :, :, : )

  INTEGER, EXTERNAL :: atomic_number, global_kpoint_index

  IF ( real_or_complex .EQ. 1 .OR. nspin .GT. 1 ) THEN
    proc_wf = .TRUE.
  ELSE
    proc_wf = .FALSE.
  ENDIF

  bad_kgrid = .FALSE.
  IF ( wfng_kgrid ) THEN
    IF ( wfng_nk1 .LE. 0 .OR. wfng_nk2 .LE. 0 .OR. wfng_nk3 .LE. 0 ) &
      bad_kgrid = .TRUE.
  ELSE
    IF ( nk1 .LE. 0 .OR. nk2 .LE. 0 .OR. nk3 .LE. 0 ) &
      bad_kgrid = .TRUE.
  ENDIF
  IF ( bad_kgrid .AND. ionode ) THEN
    WRITE ( 6, 101 )
  ENDIF

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  IF ( real_or_complex .EQ. 1 ) THEN
    WRITE ( stitle, '("WFN-Real",24X)' )
  ELSE
    WRITE ( stitle, '("WFN-Complex",21X)' )
  ENDIF

  unit = 4
  nrecord = 1
  nd = 3

  nb = nbnd
  nk_l = nks
  nk_g = nkstot
  ! ZL: comment out, instead, use set_spin()
  !ns = nspin
  call set_spin(ns, nst, nsf, nspin)
  ng_l = ngm
  ng_g = ngm_g

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  ALLOCATE ( kmap ( nk_g ) )
  ALLOCATE ( smap ( nk_g ) )

  DO i = 1, nk_g
    j = ( i - 1 ) / ns
    k = i - 1 - j * ns
    kmap ( i ) = j + k * ( nk_g / ns ) + 1
    smap ( i ) = k + 1
  ENDDO
  ierr = 0
  DO i = 1, nk_g
    ik = kmap ( i )
    is = smap ( i )
    IF ( ik .GE. iks .AND. ik .LE. ike .AND. is .NE. isk ( ik ) ) &
      ierr = ierr + 1
  ENDDO
  CALL mp_max ( ierr, world_comm )
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_wfng', 'smap', ierr )

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_wfng', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2 )
    t1 ( 1 ) = ft ( 1, i )
    t1 ( 2 ) = ft ( 2, i )
    t1 ( 3 ) = ft ( 3, i )
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true., translation )

  ALLOCATE ( et_g ( nb, nk_g ) )

  DO ik = 1, nk_l
    DO ib = 1, nb
      et_g ( ib, ik ) = et ( ib, ik )
    ENDDO
  ENDDO
#if defined(__MPI)
  CALL poolrecover ( et_g, nb, nk_g, nk_l )
  CALL mp_bcast ( et_g, ionode_id, world_comm )
#endif

  ALLOCATE ( wg_g ( nb, nk_g ) )
  ALLOCATE ( ifmin ( nk_g ) )
  ALLOCATE ( ifmax ( nk_g ) )

  IF ( wfng_occupation ) THEN

    DO ik = 1, nk_g
      DO ib = 1, nb
        IF ( ib .GE. wfng_nvmin .AND. ib .LE. wfng_nvmax ) THEN
          wg_g ( ib, ik ) = 1.0D0
        ELSE
          wg_g ( ib, ik ) = 0.0D0
        ENDIF
      ENDDO
    ENDDO
    DO ik = 1, nk_g
      ifmin ( ik ) = wfng_nvmin
    ENDDO
    DO ik = 1, nk_g
      ifmax ( ik ) = wfng_nvmax
    ENDDO

  ELSE

    DO ik = 1, nk_l
      DO ib = 1, nb
        IF ( wk(ik) == 0.D0 ) THEN
          wg_g(ib,ik) = wg(ib,ik)
        ELSE
          wg_g(ib,ik) = wg(ib,ik) / wk(ik)
        ENDIF
      ENDDO
    ENDDO
#if defined(__MPI)
    CALL poolrecover ( wg_g, nb, nk_g, nk_l )
#endif
    DO ik = 1, nk_g
      ifmin ( ik ) = 0
    ENDDO
    DO ik = 1, nk_g
      ifmax ( ik ) = 0
    ENDDO
    DO ik = 1, nk_g
      DO ib = 1, nb
        IF ( wg_g( ib, ik ) .GT. 0.5D0 ) THEN
          IF ( ifmin ( ik ) .EQ. 0 ) ifmin ( ik ) = ib
          ifmax ( ik ) = ib
        ENDIF
      ENDDO
    ENDDO

  ENDIF

  ALLOCATE ( g_g ( nd, ng_g ) )

  DO ig = 1, ng_g
    DO id = 1, nd
      g_g ( id, ig ) = 0
    ENDDO
  ENDDO
  DO ig = 1, ng_l
    g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO
  CALL mp_sum ( g_g, intra_bgrp_comm )

  ALLOCATE ( igk_l2g ( npwx, nk_l ) )

  DO ik = 1, nk_l
    npw = ngk ( ik )
    DO ig = 1, npw
      igk_l2g ( ig, ik ) = ig_l2g ( igk_k (ig, ik) )
    ENDDO
    DO ig = npw + 1, npwx
      igk_l2g ( ig, ik ) = 0
    ENDDO
  ENDDO

  ALLOCATE ( ngk_g ( nk_g ) )

  DO ik = 1, nk_g
    ngk_g ( ik ) = 0
  ENDDO
  DO ik = 1, nk_l
    ngk_g ( ik + iks - 1 ) = ngk ( ik )
  ENDDO
  CALL mp_sum( ngk_g, inter_pool_comm )
  CALL mp_sum( ngk_g, intra_pool_comm )
  ngk_g = ngk_g / nbgrp

  npw_g = MAXVAL ( igk_l2g ( :, : ) )
  CALL mp_max( npw_g, intra_pool_comm )
  CALL mp_max( npw_g, inter_pool_comm )

  npwx_g = MAXVAL ( ngk_g ( : ) )

  CALL cryst_to_cart ( nk_g / ns, xk, at, - 1 )

  IF ( ionode ) THEN
    OPEN ( unit = unit, file = TRIM ( output_file_name ), &
      form = 'unformatted', status = 'replace' )
    WRITE ( unit ) stitle, sdate, stime
    ! ZL: change ns to nsf
    !WRITE ( unit ) ns, ng_g, ntran, cell_symmetry, nat, ecutrho, &
    WRITE ( unit ) nsf, ng_g, ntran, cell_symmetry, nat, ecutrho, &
      nk_g / ns, nb, npwx_g, ecutwfc
    IF ( wfng_kgrid ) THEN
      WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3, wfng_nk1, wfng_nk2, wfng_nk3, &
        wfng_dk1, wfng_dk2, wfng_dk3
    ELSE
      WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3, nk1, nk2, nk3, &
        0.5D0 * dble ( k1 ), 0.5D0 * dble ( k2 ), 0.5D0 * dble ( k3 )
    ENDIF
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    WRITE ( unit ) ( ngk_g ( ik ), ik = 1, nk_g / ns )
    ! ZL: change ns to nst
    !WRITE ( unit ) ( wk ( ik ) * dble ( ns ) / 2.0D0, ik = 1, nk_g / ns )
    WRITE ( unit ) ( wk ( ik ) * dble ( nst ) / 2.0D0, ik = 1, nk_g / ns )
    WRITE ( unit ) ( ( xk ( id, ik ), id = 1, nd ), ik = 1, nk_g / ns )
    WRITE ( unit ) ( ifmin ( ik ), ik = 1, nk_g )
    WRITE ( unit ) ( ifmax ( ik ), ik = 1, nk_g )
    WRITE ( unit ) ( ( et_g ( ib, ik ), ib = 1, nb ), ik = 1, nk_g )
    WRITE ( unit ) ( ( wg_g ( ib, ik ), ib = 1, nb ), ik = 1, nk_g )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
  ENDIF

  CALL cryst_to_cart ( nk_g / ns, xk, bg, 1 )

  DEALLOCATE ( wg_g )
  DEALLOCATE ( ifmax )
  DEALLOCATE ( ifmin )

  ALLOCATE ( igwk ( npwx_g ) )

  IF ( proc_wf ) THEN
    IF ( MOD ( npwx_g, nproc ) .EQ. 0 ) THEN
      ngkdist_l = npwx_g / nproc
    ELSE
      ngkdist_l = npwx_g / nproc + 1
    ENDIF
    ngkdist_g = ngkdist_l * nproc
    IF ( real_or_complex .EQ. 1 ) &
      ALLOCATE ( energy ( nb, ns ) ) ! ZL note: for noncollinear, there is no spin up/down
    ! ZL: change ns to nst, because wavefunctions are now two-component spinors
    !ALLOCATE ( wfng_buf ( ngkdist_g, ns ) )
    !ALLOCATE ( wfng_dist ( ngkdist_l, nb, ns ) )
    ALLOCATE ( wfng_buf ( ngkdist_g, nst ) )
    ALLOCATE ( wfng_dist ( ngkdist_l, nb, nst ) )
  ENDIF

  DO i = 1, nk_g

    ik = kmap ( i )
    is = smap ( i )

    IF ( real_or_complex .EQ. 1 ) THEN
      DO ib = 1, nb
        energy ( ib, is ) = et_g ( ib, i )
      ENDDO
    ENDIF

    DO j = 1, npwx_g
      igwk ( j ) = 0
    ENDDO
    ALLOCATE ( itmp ( npw_g ) )
    DO j = 1, npw_g
      itmp ( j ) = 0
    ENDDO
    IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
      DO ig = 1, ngk ( ik - iks + 1 )
        itmp ( igk_l2g ( ig, ik - iks + 1 ) ) = igk_l2g ( ig, ik - iks + 1 )
      ENDDO
    ENDIF
    CALL mp_sum( itmp, intra_bgrp_comm )
    !XXX
    CALL mp_sum( itmp, inter_pool_comm )
    !XXX
    ngg = 0
    DO ig = 1, npw_g
      IF ( itmp ( ig ) .EQ. ig ) THEN
        ngg = ngg + 1
        igwk ( ngg ) = ig
      ENDIF
    ENDDO
    DEALLOCATE ( itmp )

    IF ( ionode ) THEN
      IF ( is .EQ. 1 ) THEN
        WRITE ( unit ) nrecord
        WRITE ( unit ) ngk_g ( ik )
        WRITE ( unit ) ( ( g_g ( id, igwk ( ig ) ), id = 1, nd ), &
          ig = 1, ngk_g ( ik ) )
      ENDIF
    ENDIF

    local_pw = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
      CALL davcio ( evc, 2*nwordwfc, iunwfc, ik - iks + 1, - 1 )
      local_pw = ngk ( ik - iks + 1 )
    ENDIF

    ALLOCATE ( igwf_l2g ( local_pw ) )

    DO ig = 1, local_pw
      igwf_l2g ( ig ) = 0
    ENDDO
    DO ig = 1, local_pw
      ngg = igk_l2g ( ig, ik - iks + 1 )
      DO j = 1, ngk_g ( ik )
        IF ( ngg .EQ. igwk ( j ) ) THEN
          igwf_l2g ( ig ) = j
          EXIT
        ENDIF
      ENDDO
    ENDDO

    ALLOCATE ( ipmask ( nproc ) )
    DO j = 1, nproc
      ipmask ( j ) = 0
    ENDDO
    ipsour = ionode_id
    IF ( npool .GT. 1 ) THEN
      IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
        IF ( me_pool .EQ. root_pool ) ipmask ( mpime + 1 ) = 1
      ENDIF
      CALL mp_sum ( ipmask, world_comm )
      DO j = 1, nproc
        IF ( ipmask ( j ) .EQ. 1 ) ipsour = j - 1
      ENDDO
    ENDIF
    DEALLOCATE ( ipmask )

    igwx = 0
    ierr = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) &
      igwx = MAXVAL ( igwf_l2g ( 1 : local_pw ) )
    CALL mp_max ( igwx, intra_pool_comm )
    IF ( ipsour .NE. ionode_id ) &
      CALL mp_get ( igwx, igwx, mpime, ionode_id, ipsour, 1, world_comm )
    ierr = 0
    IF ( ik .GE. iks .AND. ik .LE. ike .AND. igwx .NE. ngk_g ( ik ) ) &
      ierr = 1
    CALL mp_max ( ierr, world_comm )
    IF ( ierr .GT. 0 ) &
      CALL errore ( 'write_wfng', 'igwx ngk_g', ierr )

    ! ZL: add
    IF (nspin .NE. 4) THEN
      ALLOCATE ( wfng ( MAX ( 1, igwx ) ) )
    ELSE
      ALLOCATE ( wfng_nc ( MAX ( 1, igwx ), nst ) )
    ENDIF

    DO ib = 1, nb

      DO j = 1, igwx
        ! ZL: add
        IF (nspin .NE. 4) THEN
          wfng ( j ) = ( 0.0D0, 0.0D0 )
        ELSE
          DO isp = 1, nst
            wfng_nc ( j ,isp) = ( 0.0D0, 0.0D0 )
          ENDDO
        ENDIF
      ENDDO
      IF ( npool .GT. 1 ) THEN
        IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
          ! ZL: add for spinors
          IF ( nspin .NE. 4) THEN
            CALL mergewf ( evc ( :, ib ), wfng, local_pw, igwf_l2g, &
              me_pool, nproc_pool, root_pool, intra_pool_comm )
          ELSE
            DO isp = 1, nst
              CALL mergewf ( evc ( 1+npwx*(isp-1):npwx*isp, ib ), wfng_nc(:,isp), local_pw, igwf_l2g, &
                me_pool, nproc_pool, root_pool, intra_pool_comm )
            ENDDO
          ENDIF
        ENDIF
        IF ( ipsour .NE. ionode_id ) THEN
          ! ZL: add for spinors
          IF (nspin .NE. 4) THEN
            CALL mp_get ( wfng, wfng, mpime, ionode_id, ipsour, ib, &
              world_comm )
          ELSE
            DO isp=1, nst
              CALL mp_get ( wfng_nc(:,isp), wfng_nc(:,isp), mpime, ionode_id, ipsour, ib, &
                world_comm )
            ENDDO
          ENDIF
        ENDIF
      ELSE
        ! ZL: add for spinors
        IF ( nspin .NE. 4) THEN
          CALL mergewf ( evc ( :, ib ), wfng, local_pw, igwf_l2g, &
            mpime, nproc, ionode_id, world_comm )
        ELSE
          DO isp=1, nst
            CALL mergewf ( evc ( 1 + (isp-1)*npwx : isp*npwx, ib ), wfng_nc(:,isp), &
              local_pw, igwf_l2g, mpime, nproc, ionode_id, world_comm )
          ENDDO
        ENDIF
      ENDIF

      IF ( proc_wf ) THEN
        DO ig = 1, igwx
          ! ZL: add for spinors
          IF ( nspin .NE. 4) THEN
            wfng_buf ( ig, is ) = wfng ( ig )
          ELSE
            DO isp=1, nst
              wfng_buf ( ig, isp ) = wfng_nc ( ig ,isp)
            ENDDO
          ENDIF
        ENDDO
        DO ig = igwx + 1, ngkdist_g
          IF ( nspin .NE. 4) THEN
            wfng_buf ( ig, is ) = ( 0.0D0, 0.0D0 )
          ELSE
            DO isp=1, nst
              wfng_buf ( ig, isp ) = ( 0.0D0, 0.0D0 )
            ENDDO
          ENDIF
        ENDDO
#if defined(__MPI)
        CALL mp_barrier ( world_comm )
        ! ZL: add for spinors
        IF (nspin .NE. 4) THEN
          CALL MPI_Scatter ( wfng_buf ( :, is ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
          wfng_dist ( :, ib, is ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
          ionode_id, world_comm, ierr )
        ELSE
          DO isp = 1, nst
            CALL MPI_Scatter ( wfng_buf ( :, isp ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
            wfng_dist ( :, ib, isp ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
            ionode_id, world_comm, ierr )
          ENDDO
        ENDIF
        IF ( ierr .GT. 0 ) &
          CALL errore ( 'write_wfng', 'mpi_scatter', ierr )
#else
        DO ig = 1, ngkdist_g
          ! ZL: add for spinors
          IF (nspin .NE. 4) THEN
            wfng_dist ( ig, ib, is ) = wfng_buf ( ig, is )
          ELSE
            DO isp = 1, nst
              wfng_dist ( ig, ib, isp ) = wfng_buf ( ig, isp )
            ENDDO
          ENDIF
        ENDDO
#endif
      ELSE
        IF ( ionode ) THEN
          WRITE ( unit ) nrecord
          WRITE ( unit ) ngk_g ( ik )
          ! ZL: add for spinors
          IF (nspin.NE.4) THEN
            WRITE ( unit ) ( wfng ( ig ), ig = 1, igwx )
          ELSE
            WRITE ( unit ) (( wfng_nc ( ig, isp ), ig = 1, igwx), isp = 1, nst )
          ENDIF
        ENDIF
      ENDIF

    ENDDO

    ! ZL: add for spinors
    IF (nspin .NE. 4) THEN
      DEALLOCATE ( wfng )
    ELSE
      DEALLOCATE ( wfng_nc )
    ENDIF
    DEALLOCATE ( igwf_l2g )

    IF ( proc_wf .AND. is .EQ. ns ) THEN
      IF ( real_or_complex .EQ. 1 ) THEN
        CALL start_clock ( 'real_wfng' )
        CALL real_wfng ( ik, ngkdist_l, nb, ns, energy, wfng_dist )
        CALL stop_clock ( 'real_wfng' )
      ENDIF
      DO ib = 1, nb
        ! ZL: change ns to nst
        !DO is = 1, ns
        DO is = 1, nst
#if defined(__MPI)
          CALL mp_barrier ( world_comm )
          CALL MPI_Gather ( wfng_dist ( :, ib, is ), ngkdist_l, &
          MPI_DOUBLE_COMPLEX, wfng_buf ( :, is ), ngkdist_l, &
          MPI_DOUBLE_COMPLEX, ionode_id, world_comm, ierr )
          IF ( ierr .GT. 0 ) &
            CALL errore ( 'write_wfng', 'mpi_gather', ierr )
#else
          DO ig = 1, ngkdist_g
            wfng_buf ( ig, is ) = wfng_dist ( ig, ib, is )
          ENDDO
#endif
        ENDDO
        IF ( ionode ) THEN
          WRITE ( unit ) nrecord
          WRITE ( unit ) ngk_g ( ik )
          IF ( real_or_complex .EQ. 1 ) THEN
            WRITE ( unit ) ( ( dble ( wfng_buf ( ig, is ) ), &
              ! ZL: change ns to nst
              !ig = 1, igwx ), is = 1, ns )
              ig = 1, igwx ), is = 1, nst )
          ELSE
            WRITE ( unit ) ( ( wfng_buf ( ig, is ), &
              ! ZL: change ns to nst
              !ig = 1, igwx ), is = 1, ns )
              ig = 1, igwx ), is = 1, nst )
          ENDIF
        ENDIF
      ENDDO
    ENDIF

  ENDDO

  DEALLOCATE ( igwk )
  DEALLOCATE ( ngk_g )
  DEALLOCATE ( igk_l2g )
  DEALLOCATE ( et_g )

  IF ( proc_wf ) THEN
    IF ( real_or_complex .EQ. 1 ) &
    DEALLOCATE ( energy )
    DEALLOCATE ( wfng_buf )
    DEALLOCATE ( wfng_dist )
  ENDIF

  IF ( ionode ) THEN
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( g_g )
  DEALLOCATE ( smap )
  DEALLOCATE ( kmap )

  CALL mp_barrier ( world_comm )

  RETURN

101 FORMAT ( /, 5X, "WARNING: kgrid is set to zero in the wavefunction file.", &
             /, 14X, "The resulting file will only be usable as the fine grid in inteqp.", / )

END SUBROUTINE write_wfng

!-------------------------------------------------------------------------------

SUBROUTINE real_wfng ( ik, ngkdist_l, nb, ns, energy, wfng_dist )

  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode
  USE lsda_mod, ONLY : nspin  ! ZL: add for spinors information check
  USE mp, ONLY : mp_sum
  USE mp_world, ONLY : world_comm

  IMPLICIT NONE

  integer, intent (in) :: ik, ngkdist_l, nb, ns
  real (DP), intent (in) :: energy ( :, : ) ! ( nb, ns )
  complex (DP), intent (inout) :: wfng_dist ( :, :, : ) ! ( ngkdist_l, nb, ns )

  real (DP), PARAMETER :: eps2 = 1.0D-2
  real (DP), PARAMETER :: eps5 = 1.0D-5
  real (DP), PARAMETER :: eps6 = 1.0D-6

  character :: tmpstr*80
  integer :: i, j, k, is, ib, jb, ig, inum, deg, mdeg, inc
  integer :: dimension_span, reduced_span, ierr
  real (DP) :: x
  integer, allocatable :: imap ( :, : )
  integer, allocatable :: inums ( : )
  integer, allocatable :: inull ( : )
  integer, allocatable :: null_map ( :, : )
  real (DP), allocatable :: psi ( :, : )
  real (DP), allocatable :: phi ( :, : )
  real (DP), allocatable :: vec ( : )
  complex (DP), allocatable :: wfc ( : )

  IF (nspin.EQ.4)  CALL errore('real_wfng not allowed for full-spinor case','nspin',nspin)

  mdeg = 1
  DO is = 1, ns
    DO ib = 1, nb - 1
      deg = 1
      DO jb = ib + 1, nb
        IF ( abs ( energy ( ib, is ) - energy ( jb, is ) ) &
          .LT. eps5 * dble ( jb - ib + 1 ) ) deg = deg + 1
      ENDDO
      IF ( deg .GT. mdeg ) mdeg = deg
    ENDDO
  ENDDO
  mdeg = mdeg * 2

  ALLOCATE ( imap ( nb, ns ) )
  ALLOCATE ( inums ( ns ) )
  ALLOCATE ( inull ( nb ) )
  ALLOCATE ( null_map ( mdeg, nb ) )

  DO is = 1, ns
    inum = 1
    DO ib = 1, nb
      IF ( ib .EQ. nb ) THEN
        imap ( inum, is ) = ib
        inum = inum + 1
      ELSEIF ( abs ( energy ( ib, is ) - &
        energy ( ib + 1, is ) ) .GT. eps5 ) THEN
        imap ( inum, is ) = ib
        inum = inum + 1
      ENDIF
    ENDDO
    inum = inum - 1
    inums ( is ) = inum
  ENDDO

  ALLOCATE ( wfc ( ngkdist_l ) )
  ALLOCATE ( psi ( ngkdist_l, mdeg ) )
  ALLOCATE ( phi ( ngkdist_l, mdeg ) )
  ALLOCATE ( vec ( ngkdist_l ) )

  DO is = 1, ns
    inc = 1
    inum = inums ( is )
    DO i = 1, inum
      inull ( i ) = 1
      DO ib = inc, imap ( i, is )
        DO ig = 1, ngkdist_l
          wfc ( ig ) = wfng_dist ( ig, ib, is )
        ENDDO
        x = 0.0D0
        DO ig = 1, ngkdist_l
          x = x + dble ( wfc ( ig ) ) **2
        ENDDO
        CALL mp_sum ( x, world_comm )
        IF ( x .LT. eps2 ) null_map ( inull ( i ), i ) = 0
        IF ( x .GT. eps2 ) null_map ( inull ( i ), i ) = 1
        inull ( i ) = inull ( i ) + 1
        x = 0.0D0
        DO ig = 1, ngkdist_l
          x = x + aimag ( wfc ( ig ) ) **2
        ENDDO
        CALL mp_sum ( x, world_comm )
        IF ( x .LT. eps2 ) null_map ( inull ( i ), i ) = 0
        IF ( x .GT. eps2 ) null_map ( inull ( i ), i ) = 1
        inull ( i ) = inull ( i ) + 1
      ENDDO
      inull ( i ) = inull ( i ) - 1
      inc = imap ( i, is ) + 1
    ENDDO
    inc = 1
    ib = 1
    DO i = 1, inum
      k = 1
      DO j = 1, 2 * ( imap ( i, is ) - inc ) + 1, 2
        IF ( null_map ( j, i ) .EQ. 1 .OR. &
          null_map ( j + 1, i ) .EQ. 1 ) THEN
          DO ig = 1, ngkdist_l
            wfc ( ig ) = wfng_dist ( ig, ib, is )
          ENDDO
          IF ( null_map ( j, i ) .EQ. 1 ) THEN
            DO ig = 1, ngkdist_l
              phi ( ig, k ) = dble ( wfc ( ig ) )
            ENDDO
            k = k + 1
          ENDIF
          IF ( null_map ( j + 1, i ) .EQ. 1 ) THEN
            DO ig = 1, ngkdist_l
              phi ( ig, k ) = aimag ( wfc ( ig ) )
            ENDDO
            k = k + 1
          ENDIF
          ib = ib + 1
        ENDIF
      ENDDO
      dimension_span = k - 1
      IF ( dimension_span .EQ. 0 ) THEN
        ierr = 201
        WRITE ( tmpstr, 201 ) ik, is, inc
        CALL errore ( 'real_wfng', tmpstr, ierr )
      ENDIF
      DO j = 1, dimension_span
        x = 0.0D0
        DO ig = 1, ngkdist_l
          x = x + phi ( ig, j ) **2
        ENDDO
        CALL mp_sum ( x, world_comm )
        x = sqrt ( x )
        DO ig = 1, ngkdist_l
          phi ( ig, j ) = phi ( ig, j ) / x
        ENDDO
      ENDDO
!
! the Gram-Schmidt process begins
!
      reduced_span = 1
      DO ig = 1, ngkdist_l
        psi ( ig, 1 ) = phi ( ig, 1 )
      ENDDO
      DO j = 1, dimension_span - 1
        DO ig = 1, ngkdist_l
          vec ( ig ) = phi ( ig, j + 1 )
        ENDDO
        DO k = 1, reduced_span
          x = 0.0D0
          DO ig = 1, ngkdist_l
            x = x + phi ( ig, j + 1 ) * psi ( ig, k )
          ENDDO
          CALL mp_sum ( x, world_comm )
          DO ig = 1, ngkdist_l
            vec ( ig ) = vec ( ig ) - psi ( ig, k ) * x
          ENDDO
        ENDDO
        x = 0.0D0
        DO ig = 1, ngkdist_l
          x = x + vec ( ig ) **2
        ENDDO
        CALL mp_sum ( x, world_comm )
        x = sqrt ( x )
        IF ( x .GT. eps6 ) THEN
          reduced_span = reduced_span + 1
          DO ig = 1, ngkdist_l
            psi ( ig, reduced_span ) = vec ( ig ) / x
          ENDDO
        ENDIF
      ENDDO
!
! the Gram-Schmidt process ends
!
      IF ( reduced_span .LT. imap ( i, is ) - inc + 1 ) THEN
        ierr = 202
        WRITE ( tmpstr, 202 ) ik, is, inc
        CALL errore ( 'real_wfng', tmpstr, ierr )
      ENDIF
      DO ib = inc, imap ( i, is )
        DO ig = 1, ngkdist_l
          wfng_dist ( ig, ib, is ) = &
          CMPLX ( psi ( ig, ib - inc + 1 ), 0.0D0, KIND=dp )
        ENDDO
      ENDDO
      inc = imap ( i, is ) + 1
    ENDDO
  ENDDO

  DEALLOCATE ( vec )
  DEALLOCATE ( phi )
  DEALLOCATE ( psi )
  DEALLOCATE ( wfc )
  DEALLOCATE ( null_map )
  DEALLOCATE ( inull )
  DEALLOCATE ( inums )
  DEALLOCATE ( imap )

  RETURN

201 FORMAT("failed Gram-Schmidt dimension span for kpoint =",i6," spin =",i2," band =",i6)
202 FORMAT("failed Gram-Schmidt reduced span for kpoint =",i6," spin =",i2," band =",i6)

END SUBROUTINE real_wfng

!-------------------------------------------------------------------------------

SUBROUTINE write_rhog ( output_file_name, real_or_complex, symm_type, &
  rhog_nvmin, rhog_nvmax )

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  USE constants, ONLY : pi, tpi, eps6
  USE fft_base, ONLY : dfftp
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, mill, ecutrho
  USE io_global, ONLY : ionode
  USE ions_base, ONLY : nat, atm, ityp, tau 
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, rhoz_or_updw
  USE symm_base, ONLY : s, ft, nsym
  USE matrix_inversion

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  integer, intent (in) :: real_or_complex
  character ( len = 9 ), intent (in) :: symm_type
  integer, intent (in) :: rhog_nvmin
  integer, intent (in) :: rhog_nvmax

  character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  integer :: unit, id, is, ig, i, j, k, ierr
  integer :: nd, ns, ng_l, ng_g, nst, nsf  ! ZL: add nst, nsf
  integer :: ntran, cell_symmetry, nrecord
  real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: g_g ( :, : )
  complex (DP), allocatable :: rhog_g ( :, : )

  INTEGER, EXTERNAL :: atomic_number

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  IF ( real_or_complex .EQ. 1 ) THEN
    WRITE ( stitle, '("RHO-Real",24X)' )
  ELSE
    WRITE ( stitle, '("RHO-Complex",21X)' )
  ENDIF

  unit = 4
  nrecord = 1
  nd = 3

  ! ZL: use set_spin instead
  !ns = nspin
  call set_spin(ns, nst, nsf, nspin)
  ng_l = ngm
  ng_g = ngm_g

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_rhog', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2 )
    t1 ( 1 ) = ft ( 1, i )
    t1 ( 2 ) = ft ( 2, i )
    t1 ( 3 ) = ft ( 3, i )
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true., translation )

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  IF ( rhog_nvmin .NE. 0 .AND. rhog_nvmax .NE. 0 ) &
    CALL calc_rhog ( rhog_nvmin, rhog_nvmax )
    !
    IF ( nspin==2 ) CALL rhoz_or_updw(rho, 'only_g', '->updw')
    !
  ALLOCATE ( g_g ( nd, ng_g ) )
  ALLOCATE ( rhog_g ( ng_g, ns ) )

  DO ig = 1, ng_g
    DO id = 1, nd
      g_g ( id, ig ) = 0
    ENDDO
  ENDDO
  DO is = 1, ns
    DO ig = 1, ng_g
      rhog_g ( ig, is ) = ( 0.0D0, 0.0D0 )
    ENDDO
  ENDDO

  DO ig = 1, ng_l
    g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO
  DO is = 1, ns
    DO ig = 1, ng_l
      rhog_g ( ig_l2g ( ig ), is ) = rho%of_g ( ig, is )
    ENDDO
  ENDDO

  CALL mp_sum ( g_g, intra_bgrp_comm )
  CALL mp_sum ( rhog_g, intra_bgrp_comm )

  DO is = 1, ns
    DO ig = 1, ng_g
      rhog_g ( ig, is ) = rhog_g ( ig, is ) * CMPLX ( omega, 0.0D0, KIND=dp )
    ENDDO
  ENDDO

  IF ( ionode ) THEN
    OPEN ( unit = unit, file = TRIM ( output_file_name ), &
      form = 'unformatted', status = 'replace' )
    WRITE ( unit ) stitle, sdate, stime
    ! ZL: change ns to nsf
    !WRITE ( unit ) ns, ng_g, ntran, cell_symmetry, nat, ecutrho
    WRITE ( unit ) nsf, ng_g, ntran, cell_symmetry, nat, ecutrho
    WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    IF ( real_or_complex .EQ. 1 ) THEN
      WRITE ( unit ) ( ( dble ( rhog_g ( ig, is ) ), &
        ig = 1, ng_g ), is = 1, ns )
    ELSE
      WRITE ( unit ) ( ( rhog_g ( ig, is ), &
        ig = 1, ng_g ), is = 1, ns )
    ENDIF
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( rhog_g )
  DEALLOCATE ( g_g )

  RETURN

END SUBROUTINE write_rhog

!-------------------------------------------------------------------------------

SUBROUTINE calc_rhog (rhog_nvmin, rhog_nvmax)

! calc_rhog    Originally By Brad D. Malone    Last Modified (night before his thesis defense)
! Computes charge density by summing over a subset of occupied bands

  USE cell_base, ONLY : omega, tpiba2
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE gvect, ONLY : ngm, g
  USE io_files, ONLY : nwordwfc, iunwfc
  USE klist, ONLY : xk, nkstot, ngk, nks, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : nbgrp, inter_bgrp_comm
  USE noncollin_module, ONLY : nspin_mag
  USE scf, ONLY : rho
  USE symme, ONLY : sym_rho, sym_rho_init
  USE wavefunctions, ONLY : evc, psic, psic_nc  ! ZL: add psic_nc for spinors
  USE wvfct, ONLY : wg, npwx ! ZL: add npwx

  IMPLICIT NONE

  integer, intent (in) :: rhog_nvmin
  integer, intent (in) :: rhog_nvmax
  integer, external :: global_kpoint_index
  integer :: npw, ik, is, ib, ig, ir, iks, ike, isp, nst  ! ZL: add isp, nst

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  CALL weights ()

  rho%of_r (:, :) = 0.0D0

  ! take psi to R-space, compute rho in R-space
  DO ik = iks, ike
    is = isk (ik)
    npw = ngk ( ik - iks + 1 )
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    DO ib = rhog_nvmin, rhog_nvmax
      psic (:) = (0.0D0, 0.0D0)
      if (nspin == 4) psic_nc (:,:) = (0.0D0, 0.0D0)
      DO ig = 1, npw
        ! ZL: add for spinors
        IF (nspin == 4) THEN
          DO isp =1,nst
            psic_nc (dfftp%nl (igk_k (ig, ik-iks+1)), isp) = evc (ig+npwx*(isp-1), ib)
          ENDDO
        ELSE
          psic (dfftp%nl (igk_k (ig, ik-iks+1))) = evc (ig, ib)
        ENDIF
      ENDDO
      ! ZL: add for spinors
      IF (nspin == 4) THEN
        DO isp=1,nst
          CALL invfft ('Rho', psic_nc(:,isp), dfftp)
        ENDDO
      ELSE
        CALL invfft ('Rho', psic, dfftp)
      ENDIF
      DO ir = 1, dfftp%nnr
        ! ZL: add for spinors
        IF (nspin == 4) THEN
          DO isp=1,nst
            rho%of_r (ir, is) = rho%of_r (ir, is) + wg (ib, ik) / omega &
              * (dble (psic_nc (ir,isp)) **2 + aimag (psic_nc (ir,isp)) **2)
          ENDDO
        ELSE
          rho%of_r (ir, is) = rho%of_r (ir, is) + wg (ib, ik) / omega &
            * (dble (psic (ir)) **2 + aimag (psic (ir)) **2)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  CALL mp_sum (rho%of_r, inter_pool_comm)
  CALL mp_sum (rho%of_r, inter_bgrp_comm)
  rho%of_r = rho%of_r / nbgrp

  ! take rho to G-space
  ! ZL comment: in the above DO-loop, is = isk(ik) where isk(:) = 1 for non-collinear,
  ! see PW/src/setup.f90 . So for spinors, after the above loop, only rho%of_r(:,1) and the below
  ! rho%of_g(:, 1) are meaningful.  Then in the subroutine write_rhog, for spinors, actually
  ! only this component is being written.
  DO is = 1, nspin
    psic (:) = (0.0D0, 0.0D0)
    psic (:) = rho%of_r (:, is)
    CALL fwfft ('Rho', psic, dfftp)
    rho%of_g (:, is) = psic (dfftp%nl (:))
  ENDDO

  ! symmetrize rho (didn`t make a difference)
  CALL sym_rho_init (.False.)
  CALL sym_rho (nspin_mag, rho%of_g)

  RETURN

END SUBROUTINE calc_rhog

!-------------------------------------------------------------------------------

SUBROUTINE write_vxcg ( output_file_name, real_or_complex, symm_type, &
  vxc_zero_rho_core )

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  USE constants, ONLY : pi, tpi, eps6
  USE ener, ONLY : etxc, vtxc
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, mill, ecutrho
  USE io_global, ONLY : ionode
  USE ions_base, ONLY : nat, atm, ityp, tau 
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE symm_base, ONLY : s, ft, nsym
  USE wavefunctions, ONLY : psic
  USE matrix_inversion

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  integer, intent (in) :: real_or_complex
  character ( len = 9 ), intent (in) :: symm_type
  logical, intent (in) :: vxc_zero_rho_core

  character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  integer :: unit, id, is, ir, ig, i, j, k, ierr
  integer :: nd, ns, nr, ng_l, ng_g, nst, nsf  ! ZL: add nst, nsf
  integer :: ntran, cell_symmetry, nrecord
  real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: g_g ( :, : )
  real (DP), allocatable :: vxcr_g ( :, : )
  complex (DP), allocatable :: vxcg_g ( :, : )

  INTEGER, EXTERNAL :: atomic_number

  ! ZL: Disable VXC for spinor case
  !     this is because for SOC-mag case, the matrix elements of VXC by BGW is
  !     not supported yet
  if (nspin .eq. 4) then
    call errore('write_vxcg', 'full spinors do not support VXC, nspin =', nspin)
  endif

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  IF ( real_or_complex .EQ. 1 ) THEN
    WRITE ( stitle, '("VXC-Real",24X)' )
  ELSE
    WRITE ( stitle, '("VXC-Complex",21X)' )
  ENDIF

  unit = 4
  nrecord = 1
  nd = 3

  ! ZL: use set_spin for spinors
  !ns = nspin
  call set_spin(ns, nst, nsf, nspin)
  nr = dfftp%nnr
  ng_l = ngm
  ng_g = ngm_g

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_vxcg', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2 )
    t1 ( 1 ) = ft ( 1, i )
    t1 ( 2 ) = ft ( 2, i )
    t1 ( 3 ) = ft ( 3, i )
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true., translation )

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  ALLOCATE ( g_g ( nd, ng_g ) )
  ! ZL: change size of vxcr_g from ns to nsf
  !ALLOCATE ( vxcr_g ( nr, ns ) )
  ALLOCATE ( vxcr_g ( nr, nsf ) )
  ALLOCATE ( vxcg_g ( ng_g, ns ) )

  DO ig = 1, ng_g
    DO id = 1, nd
      g_g ( id, ig ) = 0
    ENDDO
  ENDDO
  DO is = 1, ns
    DO ig = 1, ng_g
      vxcg_g ( ig, is ) = ( 0.0D0, 0.0D0 )
    ENDDO
  ENDDO

  DO ig = 1, ng_l
    g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO
  vxcr_g ( :, : ) = 0.0D0
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF
  !
  CALL v_xc ( rho, rho_core, rhog_core, etxc, vtxc, vxcr_g )
  !
  DO is = 1, ns
    DO ir = 1, nr
      psic ( ir ) = CMPLX ( vxcr_g ( ir, is ), 0.0D0, KIND=dp )
    ENDDO
    CALL fwfft ( 'Rho', psic, dfftp )
    DO ig = 1, ng_l
      vxcg_g ( ig_l2g ( ig ), is ) = psic ( dfftp%nl ( ig ) )
    ENDDO
  ENDDO

  CALL mp_sum ( g_g, intra_bgrp_comm )
  CALL mp_sum ( vxcg_g, intra_bgrp_comm )

  IF ( ionode ) THEN
    OPEN ( unit = unit, file = TRIM ( output_file_name ), &
      form = 'unformatted', status = 'replace' )
    WRITE ( unit ) stitle, sdate, stime
    ! ZL: change ns to nsf for header
    !WRITE ( unit ) ns, ng_g, ntran, cell_symmetry, nat, ecutrho
    WRITE ( unit ) nsf, ng_g, ntran, cell_symmetry, nat, ecutrho
    WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    IF ( real_or_complex .EQ. 1 ) THEN
      WRITE ( unit ) ( ( dble ( vxcg_g ( ig, is ) ), &
        ig = 1, ng_g ), is = 1, ns )   !FZ:
    ELSE
      WRITE ( unit ) ( ( vxcg_g ( ig, is ), &
        ig = 1, ng_g ), is = 1, ns )   !FZ:
    ENDIF
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( vxcg_g )
  DEALLOCATE ( vxcr_g )
  DEALLOCATE ( g_g )

  RETURN

END SUBROUTINE write_vxcg

!-------------------------------------------------------------------------------
!FZ: added for metaGGA

SUBROUTINE write_vxcg_meta ( output_file_name, real_or_complex, symm_type, &
  vxc_zero_rho_core )

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  USE constants, ONLY : pi, tpi, eps6
  USE ener, ONLY : etxc, vtxc
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  !USE funct, ONLY : exx_is_active, dft_is_meta, get_meta     !FZ: Added for metaGGA
  USE xc_lib, ONLY : exx_is_active, xclib_dft_is     !FZ: for qe6.7
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, mill, ecutrho
  USE io_global, ONLY : ionode
  USE ions_base, ONLY : nat, atm, ityp, tau 
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE symm_base, ONLY : s, ft, nsym
  USE wavefunctions, ONLY : psic
  USE matrix_inversion

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  integer, intent (in) :: real_or_complex
  character ( len = 9 ), intent (in) :: symm_type
  logical, intent (in) :: vxc_zero_rho_core

  character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  integer :: unit, id, is, ir, ig, i, j, k, ierr
  integer :: nd, ns, nr, ng_l, ng_g, nst, nsf !FZ: added nst, nsf
  integer :: ntran, cell_symmetry, nrecord
  real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: g_g ( :, : )
  real (DP), allocatable :: vxcr_g ( :, : )  
  real (DP), allocatable :: kedtaur (:, :)   !FZ: for metaGGA 
  complex (DP), allocatable :: vxcg_g ( :, : )

  INTEGER, EXTERNAL :: atomic_number

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  IF ( real_or_complex .EQ. 1 ) THEN
    WRITE ( stitle, '("VXC-Real",24X)' )
  ELSE
    WRITE ( stitle, '("VXC-Complex",21X)' )
  ENDIF

  unit = 4
  nrecord = 1
  nd = 3

  ! ZL: disable full spinor
  if (nspin .eq. 4) then
    call errore('write_vxcg_meta', 'full spinors do not support VXC, nspin =', nspin)
  endif

  !FZ: use set_spin for spinors
  !ns = nspin
  call set_spin(ns, nst, nsf, nspin)
  nr = dfftp%nnr
  ng_l = ngm
  ng_g = ngm_g

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_vxcg', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2 )
    t1 ( 1 ) = ft ( 1, i ) 
    t1 ( 2 ) = ft ( 2, i ) 
    t1 ( 3 ) = ft ( 3, i ) 
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true., translation )

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  ALLOCATE ( g_g ( nd, ng_g ) )
  ! FZ: change size of vxcr_g from ns to nsf (nsf = 4 if nspin=4)
  !ALLOCATE ( vxcr_g ( nr, ns ) )
  ALLOCATE ( vxcr_g ( nr, nsf ) )
  ALLOCATE ( vxcg_g ( ng_g, ns ) )
  ALLOCATE (kedtaur (dfftp%nnr, nspin))   !FZ: for meta GGA 

  DO ig = 1, ng_g
    DO id = 1, nd
      g_g ( id, ig ) = 0
    ENDDO
  ENDDO
  DO is = 1, ns
    DO ig = 1, ng_g
      vxcg_g ( ig, is ) = ( 0.0D0, 0.0D0 )
    ENDDO
  ENDDO
  kedtaur (:, :) = 0.0D0    !FZ: for metaGGA

  DO ig = 1, ng_l
    g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO
  vxcr_g ( :, : ) = 0.0D0
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF
  !
  !IF (dft_is_meta() .and. (get_meta() /= 4)) then         !FZ: for metaGGA
  IF (xclib_dft_is('meta')) then         !FZ: for qe6.7
     CALL v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, vxcr_g, kedtaur )    !FZ: for metaGGA
  ELSE                           
     CALL v_xc ( rho, rho_core, rhog_core, etxc, vtxc, vxcr_g )
  ENDIF                         
  !
  DO is = 1, ns
    DO ir = 1, nr
      psic ( ir ) = CMPLX ( vxcr_g ( ir, is ), 0.0D0, KIND=dp )
    ENDDO
    CALL fwfft ( 'Rho', psic, dfftp )
    DO ig = 1, ng_l
      vxcg_g ( ig_l2g ( ig ), is ) = psic ( dfftp%nl ( ig ) )
    ENDDO
  ENDDO

  CALL mp_sum ( g_g, intra_bgrp_comm )
  CALL mp_sum ( vxcg_g, intra_bgrp_comm )

  IF ( ionode ) THEN
    OPEN ( unit = unit, file = TRIM ( output_file_name ), &
      form = 'unformatted', status = 'replace' )
    WRITE ( unit ) stitle, sdate, stime
    WRITE ( unit ) ns, ng_g, ntran, cell_symmetry, nat, ecutrho
    WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    IF ( real_or_complex .EQ. 1 ) THEN
      WRITE ( unit ) ( ( dble ( vxcg_g ( ig, is ) ), &   
        ig = 1, ng_g ), is = 1, ns )
    ELSE
      WRITE ( unit ) ( ( vxcg_g ( ig, is ), &
        ig = 1, ng_g ), is = 1, ns )
    ENDIF
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( vxcg_g )
  DEALLOCATE ( vxcr_g )
  DEALLOCATE ( g_g )
  DEALLOCATE (kedtaur)      !FZ: for metaGGA

  RETURN

END SUBROUTINE write_vxcg_meta
!FZ: end modify
!-------------------------------------------------------------------------------

SUBROUTINE write_vxc0 ( output_file_name, vxc_zero_rho_core )

  USE constants, ONLY : RYTOEV
  USE ener, ONLY : etxc, vtxc
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE gvect, ONLY : ngm, mill
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE wavefunctions, ONLY : psic

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  logical, intent (in) :: vxc_zero_rho_core

  integer :: unit
  integer :: is, ir, ig
  integer :: nd, ns, nr, ng_l
  real (DP), allocatable :: vxcr_g ( :, : )
  complex (DP), allocatable :: vxc0_g ( : )

  unit = 4
  nd = 3

  ns = nspin
  nr = dfftp%nnr
  ng_l = ngm

  ! ZL: disable full spinor
  if (nspin .eq. 4) then
    call errore('write_vxc0', 'full spinors do not support VXC, nspin =', nspin)
  endif

  !IF (nspin.EQ.4)  CALL errore('write_vxc0 not supported with full spinors','ns',ns)

  ALLOCATE ( vxcr_g ( nr, ns ) )
  ALLOCATE ( vxc0_g ( ns ) )

  DO is = 1, ns
    vxc0_g ( is ) = ( 0.0D0, 0.0D0 )
  ENDDO

  vxcr_g ( :, : ) = 0.0D0
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF
  !
  CALL v_xc ( rho, rho_core, rhog_core, etxc, vtxc, vxcr_g )
  !
  DO is = 1, ns
    DO ir = 1, nr
      psic ( ir ) = CMPLX ( vxcr_g ( ir, is ), 0.0D0, KIND=dp )
    ENDDO
    CALL fwfft ( 'Rho', psic, dfftp )
    DO ig = 1, ng_l
      IF ( mill ( 1, ig ) .EQ. 0 .AND. mill ( 2, ig ) .EQ. 0 .AND. &
        mill ( 3, ig ) .EQ. 0 ) vxc0_g ( is ) = psic ( dfftp%nl ( ig ) )
    ENDDO
  ENDDO

  CALL mp_sum ( vxc0_g, intra_bgrp_comm )

  DO is = 1, ns
    vxc0_g ( is ) = vxc0_g ( is ) * CMPLX ( RYTOEV, 0.0D0, KIND=dp )
  ENDDO

  IF ( ionode ) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    WRITE ( unit, 101 )
    DO is = 1, ns
      WRITE ( unit, 102 ) is, vxc0_g ( is )
    ENDDO
    WRITE ( unit, 103 )
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  DEALLOCATE ( vxcr_g )
  DEALLOCATE ( vxc0_g )

  RETURN

101 FORMAT ( /, 5X, "--------------------------------------------", &
             /, 5X, "spin    Re Vxc(G=0) (eV)    Im Vxc(G=0) (eV)", &
             /, 5X, "--------------------------------------------" )
102 FORMAT ( 5X, I1, 3X, 2F20.15 )
103 FORMAT ( 5X, "--------------------------------------------", / )

END SUBROUTINE write_vxc0

!-------------------------------------------------------------------------------
! FZ: added for metaGGA

SUBROUTINE write_vxc0_meta ( output_file_name, vxc_zero_rho_core )

  USE constants, ONLY : RYTOEV
  USE ener, ONLY : etxc, vtxc
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE gvect, ONLY : ngm, mill
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  !USE funct, ONLY : exx_is_active, dft_is_meta, get_meta     !FZ: Added for metaGGA
  USE xc_lib, ONLY : exx_is_active, xclib_dft_is     !FZ: Added for metaGGA
  USE mp, ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE wavefunctions, ONLY : psic

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  logical, intent (in) :: vxc_zero_rho_core

  integer :: unit
  integer :: is, ir, ig
  integer :: nd, ns, nr, ng_l
  real (DP), allocatable :: vxcr_g ( :, : )
  complex (DP), allocatable :: vxc0_g ( : )
  real (DP), allocatable :: kedtaur (:, :)   !FZ: for metaGGA   kedtau: local K energy density,  kedtaur (in realspace)

  unit = 4
  nd = 3

  ns = nspin
  nr = dfftp%nnr
  ng_l = ngm

  ! ZL: disable full spinor
  if (nspin .eq. 4) then
    call errore('write_vxc0_meta', 'full spinors not supported, nspin =', nspin)
  endif

  ALLOCATE ( vxcr_g ( nr, ns ) )
  ALLOCATE ( vxc0_g ( ns ) )
  ALLOCATE (kedtaur (dfftp%nnr, nspin))   !FZ: for meta GGA

  DO is = 1, ns
    vxc0_g ( is ) = ( 0.0D0, 0.0D0 )
  ENDDO
  kedtaur (:, :) = 0.0D0    !FZ: for metaGGA

  vxcr_g ( :, : ) = 0.0D0
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF
  !
  !IF (dft_is_meta() .and. (get_meta() /= 4)) then         !FZ: for metaGGA
  IF (xclib_dft_is('meta')) then         !FZ: for qe6.7
     CALL v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, vxcr_g, kedtaur )    !FZ: for metaGGA
  ELSE                          
     CALL v_xc ( rho, rho_core, rhog_core, etxc, vtxc, vxcr_g )  
  ENDIF                          
  !
  DO is = 1, ns
    DO ir = 1, nr
      psic ( ir ) = CMPLX ( vxcr_g ( ir, is ), 0.0D0, KIND=dp )
    ENDDO
    CALL fwfft ( 'Rho', psic, dfftp )
    DO ig = 1, ng_l
      IF ( mill ( 1, ig ) .EQ. 0 .AND. mill ( 2, ig ) .EQ. 0 .AND. &
        mill ( 3, ig ) .EQ. 0 ) vxc0_g ( is ) = psic ( dfftp%nl ( ig ) )
    ENDDO
  ENDDO

  CALL mp_sum ( vxc0_g, intra_bgrp_comm )

  DO is = 1, ns
    vxc0_g ( is ) = vxc0_g ( is ) * CMPLX ( RYTOEV, 0.0D0, KIND=dp )
  ENDDO

  IF ( ionode ) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    WRITE ( unit, 101 )
    DO is = 1, ns
      WRITE ( unit, 102 ) is, vxc0_g ( is )
    ENDDO
    WRITE ( unit, 103 )
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  DEALLOCATE ( vxcr_g )
  DEALLOCATE ( vxc0_g )
  DEALLOCATE (kedtaur)      !FZ: for metaGGA

  RETURN

101 FORMAT ( /, 5X, "--------------------------------------------", &
             /, 5X, "spin    Re Vxc(G=0) (eV)    Im Vxc(G=0) (eV)", &
             /, 5X, "--------------------------------------------" )
102 FORMAT ( 5X, I1, 3X, 2F20.15 )
103 FORMAT ( 5X, "--------------------------------------------", / )

END SUBROUTINE write_vxc0_meta
!FZ: end modify
!-------------------------------------------------------------------------------

SUBROUTINE write_vxc_r (output_file_name, diag_nmin, diag_nmax, &
  offdiag_nmin, offdiag_nmax, vxc_zero_rho_core)

  USE kinds, ONLY : DP
  USE constants, ONLY : rytoev, eps12  ! ZL: add eps12
  USE cell_base, ONLY : tpiba2, at, bg
  USE ener, ONLY : etxc, vtxc
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE gvect, ONLY : ngm, g
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global, ONLY : ionode
  USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE wavefunctions, ONLY : evc, psic, psic_nc !FZ: added psic_nc for spinors
  USE wvfct, ONLY : nbnd, npwx !FZ: added npwx
  USE symm_base, ONLY : nsym ! ZL: add
  USE noncollin_module, ONLY : magtot_nc ! ZL: add
  USE lsda_mod, ONLY : starting_magnetization ! ZL: add

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  integer, intent (inout) :: diag_nmin
  integer, intent (inout) :: diag_nmax
  integer, intent (inout) :: offdiag_nmin
  integer, intent (inout) :: offdiag_nmax
  logical, intent (in) :: vxc_zero_rho_core

  integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2
  integer, external :: global_kpoint_index
  real (DP) :: dummyr
  complex (DP) :: dummyc
  ! ZL: change mtxeld to complex to be compatible for SOC-mag case
  !real (DP), allocatable :: mtxeld (:, :)
  complex (DP), allocatable :: mtxeld (:, :)
  complex (DP), allocatable :: mtxelo (:, :, :)
  real (DP), allocatable :: vxcr (:, :)
  complex (DP), allocatable :: psic2 (:)
  complex (DP), allocatable :: psic_nc2 (:,:)  !FZ: added for spinor
  integer :: isp, ns, nst, nsf !FZ: added for spinor

  ! ZL: for vxcr_spinor_mag
  complex (DP), allocatable :: psic_nc_temp (:,:) ! ZL: added for vxcr_spinor_mag
  complex (DP), allocatable :: hpsi_nc_temp (:,:) ! ZL: added for vxcr_spinor_mag
  complex (DP) :: temp_psic_nc1, temp_psic_nc2  ! ZL: added for vxcr_spinor_mag
  complex (DP) :: dummy_aux

  ! ZL: full spinors don't support vxc.dat yet
  !if (nspin .eq. 4) then
  !  call errore('write_vxc_r', 'full spinors do not support vxc.dat, nspin =', nspin)
  !endif
  !FZ: deleted, full spinors only works when there is no magnetization
  ! ZL: -- updated Dec. 2023, now full spinors work with magnetization,
  !     but no symmetry is allowed in the SOC-mag case.

  ! ZL: soc-mag must not use symmetries
  if ( sum(abs(starting_magnetization)) .gt. eps12 ) then
    if (nsym .ne. 1) then
      call errore('write_vxc_r', 'Magnetization detected for full spinor. No symmetry is allowed, nsym =', nsym)
    endif
  endif

  if(diag_nmin > diag_nmax) then
    call errore ( 'write_vxc_r', 'diag_nmin > diag_nmax', diag_nmin )
  endif
  IF (diag_nmin .LT. 1) diag_nmin = 1
  IF (diag_nmax .GT. nbnd) then
    write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
    diag_nmax = nbnd
  ENDIF
  ndiag = MAX (diag_nmax - diag_nmin + 1, 0)

  if(offdiag_nmin > offdiag_nmax) then
    call errore ( 'write_vxc_r', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )
  endif
  IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
  IF (offdiag_nmax .GT. nbnd)  then
    write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
    offdiag_nmax = nbnd
  ENDIF
  noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

  IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

  unit = 4
  
  call set_spin(ns, nst, nsf, nspin) !FZ: added for spinor

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  IF (ndiag .GT. 0) THEN
    ALLOCATE (mtxeld (ndiag, nkstot))
    mtxeld (:, :) = (0.0D0, 0.0D0) ! ZL: change to complex
  ENDIF
  IF (noffdiag .GT. 0) THEN
    ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
    mtxelo (:, :, :) = (0.0D0, 0.0D0)
  ENDIF

  ALLOCATE (vxcr (dfftp%nnr, nspin))
  IF (noffdiag .GT. 0) ALLOCATE (psic2 (dfftp%nnr))
  IF (noffdiag .GT. 0) ALLOCATE (psic_nc2 (dfftp%nnr, 2))   !FZ: added for spinor

  vxcr (:, :) = 0.0D0
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF

  ! ZL: add
  ALLOCATE (psic_nc_temp (dfftp%nnr, 2))
  ALLOCATE (hpsi_nc_temp (dfftp%nnr, 2))

  !
  CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxcr)
  !
  DO ik = iks, ike
    npw = ngk ( ik - iks + 1 )
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    IF (ndiag .GT. 0) THEN
      DO ib = diag_nmin, diag_nmax
        psic (:) = (0.0D0, 0.0D0)

        if (nspin == 4) psic_nc (:,:) = (0.0D0, 0.0D0) ! FZ: added for spinors
        if (nspin == 4) psic_nc_temp(:, :) = (0.0D0, 0.0D0) ! ZL: add
        if (nspin == 4) hpsi_nc_temp(:, :) = (0.0D0, 0.0D0) ! ZL: add

        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              psic_nc (dfftp%nl (igk_k (ig, ik-iks+1)), isp) = evc(ig+npwx*(isp-1), ib)  !FZ: added for spinors
            ENDDO    !FZ: added for spinors
          ELSE
            psic (dfftp%nl (igk_k (ig,ik-iks+1))) = evc (ig, ib)
          ENDIF
        ENDDO
        IF (nspin == 4) THEN !FZ: added for spinors
          DO isp=1, nst  !FZ: added for spinors
            CALL invfft ('Rho', psic_nc(:,isp),dfftp)  !FZ: added for spinors
          ENDDO  !FZ: added for spinors
        ELSE  !FZ: added for spinors
          CALL invfft ('Rho', psic, dfftp)
        ENDIF  !FZ: added for spinors
        dummyr = 0.0D0
        DO ir = 1, dfftp%nnr
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp=1,nst  !FZ: added for spinors
              dummyr = dummyr + vxcr (ir, 1) &   !FZ: added for spinors
                * (dble (psic_nc (ir, isp)) **2 + aimag (psic_nc (ir,isp)) **2)   !FZ: added for spinors
            ENDDO  !FZ: added for spinors

            ! ZL: add for extra SOC-mag terms
            temp_psic_nc1 = psic_nc(ir, 1)
            temp_psic_nc2 = psic_nc(ir, 2)

            psic_nc_temp (ir, 1) =  temp_psic_nc1 * vxcr (ir, 4) + temp_psic_nc2 * ( vxcr (ir, 2) - (0.d0,1.d0) * (vxcr (ir, 3)) ) 
            psic_nc_temp (ir, 2) = -temp_psic_nc2 * vxcr (ir, 4) + temp_psic_nc1 * ( vxcr (ir, 2) + (0.d0,1.d0) * (vxcr (ir, 3)) ) 

          ELSE !FZ: added for spinors
            dummyr = dummyr + vxcr (ir, isk (ik)) &
              * (dble (psic (ir)) **2 + aimag (psic (ir)) **2)
          ENDIF !FZ: added for spinors
        ENDDO

        dummyr = dummyr * rytoev / dble (dfftp%nr1x * dfftp%nr2x * dfftp%nr3x)
        CALL mp_sum (dummyr, intra_bgrp_comm)

        ! ZL: we handle the vxc_spinor_mag part in G-space
        if (nspin == 4) then 
          do isp=1, nst
            CALL fwfft ('Rho', psic_nc_temp(:,isp),dfftp)
          enddo

          do ig = 1, npw
            do isp=1, nst
              hpsi_nc_temp (ig, isp) = psic_nc_temp (dfftp%nl (igk_k(ig,ik-iks+1)), isp)
            enddo
          enddo

          psic_nc (:,:) = (0.0D0, 0.0D0)
          do ig = 1, npw
            do isp = 1, nst
              psic_nc (ig, isp) = evc(ig+npwx*(isp-1), ib)
            enddo
          enddo

          dummy_aux = (0.0D0, 0.0D0)
          do ig = 1, npw
            do isp = 1, nst
              ! collecting terms for vxc_spinor_mag
              dummy_aux = dummy_aux + conjg (psic_nc(ig, isp)) * hpsi_nc_temp(ig, isp)
            enddo
          enddo

          dummy_aux = dummy_aux * CMPLX (rytoev, 0.0D0, KIND=dp)
          CALL mp_sum (dummy_aux, intra_bgrp_comm)

        endif
        ! ZL: done getting the matrix elements

        ! ZL: now constructing matrix elements
        if (nspin == 4) then
          mtxeld (ib - diag_nmin + 1, ik) = cmplx(dummyr, 0.0d0, kind=dp) + dummy_aux
        else
          mtxeld (ib - diag_nmin + 1, ik) = dummyr
        endif  ! ZL: done
      ENDDO
    ENDIF

    ! ZL: non-offdiag matrix elements
    IF (noffdiag .GT. 0) THEN
      DO ib = offdiag_nmin, offdiag_nmax
        psic (:) = (0.0D0, 0.0D0)

        if (nspin == 4) psic_nc (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        if (nspin == 4) psic_nc_temp (:,:) = (0.0D0, 0.0D0) ! ZL: add
        if (nspin == 4) hpsi_nc_temp (:,:) = (0.0D0, 0.0D0) ! ZL: add

        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              psic_nc (dfftp%nl (igk_k (ig, ik-iks+1)), isp) = evc(ig+npwx*(isp-1), ib)  !FZ: added for spinors
            ENDDO    !FZ: added for spinors
          ELSE   !FZ: added for spinors
            psic (dfftp%nl (igk_k (ig,ik-iks+1))) = evc (ig, ib)
          ENDIF  !FZ: added for spinors
        ENDDO
        IF (nspin == 4) THEN !FZ: added for spinors
          DO isp=1, nst  !FZ: added for spinors
            CALL invfft ('Rho', psic_nc(:,isp),dfftp)  !FZ: added for spinors
          ENDDO  !FZ: added for spinors
        ELSE  !FZ: added for spinors
          CALL invfft ('Rho', psic, dfftp)
        ENDIF  !FZ: added for spinors
        DO ib2 = offdiag_nmin, offdiag_nmax
          psic2 (:) = (0.0D0, 0.0D0)
          if (nspin == 4) psic_nc2 (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
          DO ig = 1, npw
            IF (nspin == 4) THEN !FZ: added for spinors
              DO isp = 1, nst !FZ: added for spinors
                psic_nc2 (dfftp%nl (igk_k (ig, ik-iks+1)), isp) = evc(ig+npwx*(isp-1), ib2)  !FZ: added for spinors
              ENDDO    !FZ: added for spinors
            ELSE   !FZ: added for spinors
              psic2 (dfftp%nl (igk_k (ig,ik-iks+1))) = evc (ig, ib2)
            ENDIF  !FZ: added for spinors
          ENDDO
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp=1, nst  !FZ: added for spinors
              CALL invfft ('Rho', psic_nc2(:,isp),dfftp)  !FZ: added for spinors
            ENDDO  !FZ: added for spinors
          ELSE  !FZ: added for spinors
            CALL invfft ('Rho', psic2, dfftp)
          ENDIF  !FZ: added for spinors
          dummyc = (0.0D0, 0.0D0)
          DO ir = 1, dfftp%nnr
            IF (nspin == 4) THEN !FZ: added for spinors
              DO isp=1,nst  !FZ: added for spinors
                dummyc = dummyc + CMPLX (vxcr (ir, 1), 0.0D0, KIND=dp) & !FZ: added for spinors
                  * conjg (psic_nc2 (ir, isp)) * psic_nc (ir, isp)   !FZ: added for spinors
              ENDDO  !FZ: added for spinors

              ! ZL: add for extra SOC-mag terms
              temp_psic_nc1 = psic_nc (ir, 1)
              temp_psic_nc2 = psic_nc (ir, 2)

              psic_nc_temp (ir, 1) =  temp_psic_nc1 * vxcr (ir, 4) + temp_psic_nc2 * ( vxcr (ir, 2) - (0.d0,1.d0) * vxcr (ir, 3) ) 
              psic_nc_temp (ir, 2) = -temp_psic_nc2 * vxcr (ir, 4) + temp_psic_nc1 * ( vxcr (ir, 2) + (0.d0,1.d0) * vxcr (ir, 3) ) 

            ELSE !FZ: added for spinors
              dummyc = dummyc + CMPLX (vxcr (ir, isk (ik)), 0.0D0, KIND=dp) &
                * conjg (psic2 (ir)) * psic (ir)
            ENDIF !FZ: added for spinors
          ENDDO

          dummyc = dummyc &
               * CMPLX (rytoev / dble (dfftp%nr1x * dfftp%nr2x * dfftp%nr3x), &
                        0.0D0, KIND=dp)
          CALL mp_sum (dummyc, intra_bgrp_comm)

          ! ZL: we handle the vxc_spinor_mag part in G-space
          if (nspin == 4) then
            do isp=1, nst
              CALL fwfft ('Rho', psic_nc_temp(:,isp),dfftp)
            enddo

            do ig = 1, npw
              do isp=1, nst
                hpsi_nc_temp (ig, isp) = psic_nc_temp (dfftp%nl (igk_k(ig,ik-iks+1)), isp)
              enddo
            enddo

            psic_nc2 (:,:) = (0.0D0, 0.0D0)
            do ig = 1, npw
              do isp = 1, nst
                psic_nc2 (ig, isp) = evc(ig+npwx*(isp-1), ib2)
              enddo
            enddo

            dummy_aux = (0.0D0, 0.0D0)
            do ig = 1, npw
              do isp=1, nst
                ! collecting terms for vxc_spinor_mag
                dummy_aux = dummy_aux + conjg (psic_nc2 (ig, isp)) * hpsi_nc_temp(ig, isp)
              enddo
            enddo

            dummy_aux = dummy_aux * CMPLX (rytoev, 0.0D0, KIND=dp)
            CALL mp_sum (dummy_aux, intra_bgrp_comm)

          endif ! ZL: done

          ! ZL: now constructing matrix elements
          if (nspin == 4) then
            mtxelo (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, ik) = dummyc + dummy_aux
          else
            mtxelo (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, ik) = dummyc
          endif  ! ZL: done
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  DEALLOCATE (vxcr)
  IF (noffdiag .GT. 0) DEALLOCATE (psic2)
  IF (noffdiag .GT. 0) DEALLOCATE (psic_nc2) !FZ: added for spinors

  ! ZL: add
  DEALLOCATE (psic_nc_temp)
  DEALLOCATE (hpsi_nc_temp)

  IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
  IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    IF (nspin == 4) THEN !FZ: added for spinors
      DO ik = 1, nkstot   !FZ: added for spinors
        WRITE (unit, 101) xk(:, ik), ndiag, &
          noffdiag **2
        IF (ndiag .GT. 0) THEN
          DO ib = diag_nmin, diag_nmax
            ! ZL: now mtxeld is complex
            !WRITE (unit, 102) 1, ib, mtxeld &
            !  (ib - diag_nmin + 1, ik), 0.0D0
            WRITE (unit, 102) 1, ib, mtxeld &
              (ib - diag_nmin + 1, ik)
          ENDDO
        ENDIF
        IF (noffdiag .GT. 0) THEN
          DO ib = offdiag_nmin, offdiag_nmax
            DO ib2 = offdiag_nmin, offdiag_nmax
              WRITE (unit, 103) 1, ib2, ib, mtxelo &
                (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, ik)
            ENDDO
          ENDDO
        ENDIF
      ENDDO !FZ: added for spinors
    ELSE !FZ: added for spinors
      DO ik = 1, nkstot / nspin
        WRITE (unit, 101) xk(:, ik), nspin * ndiag, &
          nspin * noffdiag **2
        IF (ndiag .GT. 0) THEN
          DO ib = diag_nmin, diag_nmax
            DO is = 1, nspin  !FZ: change order
              ! ZL: now mtxeld is complex
              !WRITE (unit, 102) is, ib, mtxeld &
              !  (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin), &
              !  0.0D0
              WRITE (unit, 102) is, ib, mtxeld &
                (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin)
            ENDDO
          ENDDO
        ENDIF
        IF (noffdiag .GT. 0) THEN
          DO ib = offdiag_nmin, offdiag_nmax
            DO ib2 = offdiag_nmin, offdiag_nmax
              DO is = 1, nspin  !FZ: change order
                WRITE (unit, 103) is, ib2, ib, mtxelo &
                  (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                  ik + (is - 1) * nkstot / nspin)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
  IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

  RETURN

  101 FORMAT (3F13.9, 2I8)
  102 FORMAT (2I8, 2F15.9)
  103 FORMAT (3I8, 2F15.9)

END SUBROUTINE write_vxc_r

!-------------------------------------------------------------------------------
!FZ: added for metaGGA
SUBROUTINE write_vxc_r_meta (output_file_name, diag_nmin, diag_nmax, &
  offdiag_nmin, offdiag_nmax, vxc_zero_rho_core)

  USE kinds, ONLY : DP
  USE constants, ONLY : rytoev
  USE cell_base, ONLY : tpiba2, at, bg
  USE ener, ONLY : etxc, vtxc
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : invfft
  USE gvect, ONLY : ngm, g
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global,        ONLY : ionode   !FZ:
  !USE funct, ONLY : exx_is_active, dft_is_meta, get_meta     !FZ: Added for metaGGA
  USE xc_lib, ONLY : exx_is_active, xclib_dft_is     !FZ: forqe6.7
  USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE wavefunctions, ONLY : evc, psic
  USE wvfct, ONLY : nbnd

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  integer, intent (inout) :: diag_nmin
  integer, intent (inout) :: diag_nmax
  integer, intent (inout) :: offdiag_nmin
  integer, intent (inout) :: offdiag_nmax
  logical, intent (in) :: vxc_zero_rho_core

  integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2
  integer, external :: global_kpoint_index
  real (DP) :: dummyr
  complex (DP) :: dummyc
  real (DP), allocatable :: mtxeld (:, :)
  complex (DP), allocatable :: mtxelo (:, :, :)
  real (DP), allocatable :: vxcr (:, :)
  real (DP), allocatable :: kedtaur (:, :)   !FZ: for metaGGA  
  complex (DP), allocatable :: psic2 (:)

  ! ZL: disable full spinor
  if (nspin .eq. 4) then
    call errore('write_vxc_r_meta', 'full spinors not supported, nspin =', nspin)
  endif

  if(diag_nmin > diag_nmax) then
    call errore ( 'write_vxc_r', 'diag_nmin > diag_nmax', diag_nmin )
  endif
  IF (diag_nmin .LT. 1) diag_nmin = 1
  IF (diag_nmax .GT. nbnd) then
    write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
    diag_nmax = nbnd
  ENDIF
  ndiag = MAX (diag_nmax - diag_nmin + 1, 0)

  if(offdiag_nmin > offdiag_nmax) then
    call errore ( 'write_vxc_r', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )
  endif
  IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
  IF (offdiag_nmax .GT. nbnd)  then
    write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
    offdiag_nmax = nbnd
  ENDIF
  noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

  IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

  unit = 4

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  IF (ndiag .GT. 0) THEN
    ALLOCATE (mtxeld (ndiag, nkstot))
    mtxeld (:, :) = 0.0D0
  ENDIF
  IF (noffdiag .GT. 0) THEN
    ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
    mtxelo (:, :, :) = (0.0D0, 0.0D0)
  ENDIF

  ALLOCATE (vxcr (dfftp%nnr, nspin))
  IF (noffdiag .GT. 0) ALLOCATE (psic2 (dfftp%nnr))

  ALLOCATE (kedtaur (dfftp%nnr, nspin))   !FZ: for meta GGA
  vxcr (:, :) = 0.0D0
  kedtaur (:, :) = 0.0D0    !FZ: for metaGGA
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF
  !
  !IF (dft_is_meta() .and. (get_meta() /= 4)) then         !FZ: for metaGGA
  IF (xclib_dft_is('meta')) then         !FZ: for qe6.7
     CALL v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, vxcr, kedtaur )    !FZ: for metaGGA
  ELSE                          
     CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxcr)
  ENDIF                        
  !
  DO ik = iks, ike
    npw = ngk ( ik - iks + 1 )
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    IF (ndiag .GT. 0) THEN
      DO ib = diag_nmin, diag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (dfftp%nl (igk_k (ig,ik-iks+1))) = evc (ig, ib)
        ENDDO
        CALL invfft ('Rho', psic, dfftp)
        dummyr = 0.0D0
        DO ir = 1, dfftp%nnr
          dummyr = dummyr + vxcr (ir, isk (ik)) &
            * (dble (psic (ir)) **2 + aimag (psic (ir)) **2)
        ENDDO
        dummyr = dummyr * rytoev / dble (dfftp%nr1x * dfftp%nr2x * dfftp%nr3x)
        CALL mp_sum (dummyr, intra_bgrp_comm)
        mtxeld (ib - diag_nmin + 1, ik) = dummyr
      ENDDO
    ENDIF
    IF (noffdiag .GT. 0) THEN
      DO ib = offdiag_nmin, offdiag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (dfftp%nl (igk_k (ig,ik-iks+1))) = evc (ig, ib)
        ENDDO
        CALL invfft ('Rho', psic, dfftp)
        DO ib2 = offdiag_nmin, offdiag_nmax
          psic2 (:) = (0.0D0, 0.0D0)
          DO ig = 1, npw
            psic2 (dfftp%nl (igk_k (ig,ik-iks+1))) = evc (ig, ib2)
          ENDDO
          CALL invfft ('Rho', psic2, dfftp)
          dummyc = (0.0D0, 0.0D0)
          DO ir = 1, dfftp%nnr
            dummyc = dummyc + CMPLX (vxcr (ir, isk (ik)), 0.0D0, KIND=dp) &
              * conjg (psic2 (ir)) * psic (ir)
          ENDDO
          dummyc = dummyc &
               * CMPLX (rytoev / dble (dfftp%nr1x * dfftp%nr2x * dfftp%nr3x), &
                        0.0D0, KIND=dp)
          CALL mp_sum (dummyc, intra_bgrp_comm)
          mtxelo (ib2 - offdiag_nmin + 1, ib - offdiag_nmin &
            + 1, ik) = dummyc
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  DEALLOCATE (vxcr)
  DEALLOCATE (kedtaur)      !FZ: for metaGGA
  IF (noffdiag .GT. 0) DEALLOCATE (psic2)

  IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
  IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    DO ik = 1, nkstot / nspin
      WRITE (unit, 101) xk(:, ik), nspin * ndiag, &
        nspin * noffdiag **2
      IF (ndiag .GT. 0) THEN
        DO ib = diag_nmin, diag_nmax
          DO is = 1, nspin  !FZ: change order
            WRITE (unit, 102) is, ib, mtxeld &
              (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin), &
              0.0D0
          ENDDO
        ENDDO
      ENDIF
      IF (noffdiag .GT. 0) THEN
        DO ib = offdiag_nmin, offdiag_nmax
          DO ib2 = offdiag_nmin, offdiag_nmax
            DO is = 1, nspin  !FZ: change order
              WRITE (unit, 103) is, ib2, ib, mtxelo &
                (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                ik + (is - 1) * nkstot / nspin)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
  IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

  RETURN

  101 FORMAT (3F13.9, 2I8)
  102 FORMAT (2I8, 2F15.9)
  103 FORMAT (3I8, 2F15.9)

END SUBROUTINE write_vxc_r_meta
!FZ: end modify
!-------------------------------------------------------------------------------

SUBROUTINE write_vxc_g (output_file_name, diag_nmin, diag_nmax, &
  offdiag_nmin, offdiag_nmax, vxc_zero_rho_core)

  USE constants, ONLY : rytoev, eps12
  USE cell_base, ONLY : tpiba2, at, bg
  USE ener, ONLY : etxc, vtxc
  USE exx, ONLY : vexx
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE xc_lib, ONLY: exx_is_active
  USE gvect, ONLY : ngm, g
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
  USE lsda_mod, ONLY : nspin, isk, starting_magnetization
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE wavefunctions, ONLY : evc, psic, psic_nc !FZ: added psic_nc for spinors
  USE wvfct, ONLY : npwx, nbnd
  USE symm_base, ONLY : nsym

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  integer, intent (inout) :: diag_nmin
  integer, intent (inout) :: diag_nmax
  integer, intent (inout) :: offdiag_nmin
  integer, intent (inout) :: offdiag_nmax
  logical, intent (in) :: vxc_zero_rho_core

  integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2, ikk
  integer, external :: global_kpoint_index
  complex (DP) :: dummy
  complex (DP), allocatable :: mtxeld (:, :)
  complex (DP), allocatable :: mtxelo (:, :, :)
  real (DP), allocatable :: vxcr (:, :)
  complex (DP), allocatable :: psic2 (:)
  complex (DP), allocatable :: psic_nc2 (:,:) !FZ: added for spinors
  complex (DP), allocatable :: hpsi (:)
  complex (DP), allocatable :: hpsi_nc (:,:)  !FZ: added for spinors
  integer :: isp, ns, nst, nsf !FZ: added for spinors
  complex (DP) :: temp_psic_nc1, temp_psic_nc2  ! ZL: added for vxcr_spinor_mag

  ! ZL: full spinors don't support vxc.dat yet
  !if (nspin .eq. 4) then
  !  call errore('write_vxc_g', 'full spinors do not support vxc.dat, nspin =', nspin)
  !endif
  !FZ: deleted, full spinors only works when there is no magnetization

  ! ZL: soc-mag must not use symmetries
  if ( sum(abs(starting_magnetization)) .gt. eps12 ) then
    if (nsym .ne. 1) then
      call errore('write_vxc_r', 'Magnetization detected for full spinor. No symmetry is allowed, nsym =', nsym)
    endif
  endif

  if(diag_nmin > diag_nmax) then
    call errore ( 'write_vxc_g', 'diag_nmin > diag_nmax', diag_nmin )
  endif
  IF (diag_nmin .LT. 1) diag_nmin = 1
  IF (diag_nmax .GT. nbnd) then
    write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
    diag_nmax = nbnd
  ENDIF
  ndiag = MAX (diag_nmax - diag_nmin + 1, 0)

  if(offdiag_nmin > offdiag_nmax) then
    call errore ( 'write_vxc_g', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )
  endif
  IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
  IF (offdiag_nmax .GT. nbnd)  then
    write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
    offdiag_nmax = nbnd
  ENDIF
  noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

  IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

  unit = 4

  call set_spin(ns, nst, nsf, nspin) !FZ: added for spinor

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  IF (ndiag .GT. 0) THEN
    ALLOCATE (mtxeld (ndiag, nkstot))
    mtxeld (:, :) = (0.0D0, 0.0D0)
  ENDIF
  IF (noffdiag .GT. 0) THEN
    ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
    mtxelo (:, :, :) = (0.0D0, 0.0D0)
  ENDIF

  ALLOCATE (vxcr (dfftp%nnr, nspin))
  IF (noffdiag .GT. 0) ALLOCATE (psic2 (dfftp%nnr))
  IF (noffdiag .GT. 0) ALLOCATE (psic_nc2 (dfftp%nnr, 2))   !FZ: added for spinor
  ALLOCATE (hpsi (dfftp%nnr))
  ALLOCATE (hpsi_nc (dfftp%nnr, 2))  !FZ: added for spinors

  vxcr (:, :) = 0.0D0
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF
  !
  CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxcr)
  !
  DO ik = iks, ike
    ikk = ik - iks + 1
    npw = ngk ( ik - iks + 1 )
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    IF (ndiag .GT. 0) THEN
      DO ib = diag_nmin, diag_nmax
        psic (:) = (0.0D0, 0.0D0)
        if (nspin == 4) psic_nc (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              psic_nc (dfftp%nl (igk_k (ig, ikk)), isp) = evc(ig+npwx*(isp-1), ib)  !FZ: added for spinors
            ENDDO    !FZ: added for spinors
          ELSE
            psic (dfftp%nl (igk_k(ig,ikk))) = evc (ig, ib)
          ENDIF
        ENDDO
        IF (nspin == 4) THEN !FZ: added for spinors
          DO isp=1, nst  !FZ: added for spinors
            CALL invfft ('Rho', psic_nc(:,isp),dfftp)  !FZ: added for spinors
          ENDDO  !FZ: added for spinors
        ELSE  !FZ: added for spinors
          CALL invfft ('Rho', psic, dfftp)
        ENDIF  !FZ: added for spinors
        DO ir = 1, dfftp%nnr
          IF (nspin == 4) THEN !FZ: added for spinors
            ! ZL: add for extra SOC-mag terms
            temp_psic_nc1 = psic_nc(ir, 1)
            temp_psic_nc2 = psic_nc(ir, 2)

            ! ZL: redo to include all terms
            !DO isp=1,nst  !FZ: added for spinors
            !  psic_nc (ir, isp) = psic_nc (ir, isp) * vxcr (ir, 1)  !FZ: addedfor spinors
            !ENDDO  !FZ: added for spinors

            ! ZL: include all terms
            psic_nc(ir, 1) = temp_psic_nc1 * vxcr (ir, 1) + temp_psic_nc1 * vxcr (ir, 4) + &
                             temp_psic_nc2 * ( vxcr (ir, 2) - (0.d0,1.d0) * (vxcr (ir, 3)) ) 
            psic_nc(ir, 2) = temp_psic_nc2 * vxcr (ir, 1) - temp_psic_nc2 * vxcr (ir, 4) + &
                             temp_psic_nc1 * ( vxcr (ir, 2) + (0.d0,1.d0) * (vxcr (ir, 3)) ) 

          ELSE  !FZ: added for spinors
            psic (ir) = psic (ir) * vxcr (ir, isk (ik))
          ENDIF  !FZ: added for spinors
        ENDDO
        IF (nspin == 4) THEN !FZ: added for spinors
          DO isp=1, nst  !FZ: added for spinors
            CALL fwfft ('Rho', psic_nc(:,isp),dfftp)  !FZ: added for spinors
          ENDDO  !FZ: added for spinors
        ELSE  !FZ: added for spinors
          CALL fwfft ('Rho', psic, dfftp)
        ENDIF  !FZ: added for spinors
        hpsi (:) = (0.0D0, 0.0D0)
        hpsi_nc (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp=1, nst  !FZ: added for spinors
              hpsi_nc (ig, isp) = psic_nc (dfftp%nl (igk_k(ig,ikk)), isp) !FZ:added for spinors
            ENDDO  !FZ: added for spinors
          ELSE  !FZ: added for spinors
            hpsi (ig) = psic (dfftp%nl (igk_k(ig,ikk)))
          ENDIF
        ENDDO
        psic (:) = (0.0D0, 0.0D0)
        if (nspin == 4) psic_nc (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              psic_nc (ig, isp) = evc(ig+npwx*(isp-1), ib)  !FZ: added for spinors
            ENDDO    !FZ: added for spinors
          ELSE   !FZ: added for spinors
            psic (ig) = evc (ig, ib)
          ENDIF  !FZ: added for spinors
        ENDDO
        !IF (exx_is_active ()) &
        !   CALL vexx (npwx, npw, 1, psic, hpsi)
        dummy = (0.0D0, 0.0D0)
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              dummy = dummy + conjg (psic_nc (ig, isp)) * hpsi_nc (ig, isp) !FZ: added for spinors
            ENDDO  !FZ: added for spinors
          ELSE   !FZ: added for spinors
            dummy = dummy + conjg (psic (ig)) * hpsi (ig)
          ENDIF  !FZ: added for spinors
        ENDDO
        dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
        CALL mp_sum (dummy, intra_bgrp_comm)
        mtxeld (ib - diag_nmin + 1, ik) = dummy
      ENDDO
    ENDIF
    IF (noffdiag .GT. 0) THEN
      DO ib = offdiag_nmin, offdiag_nmax
        psic (:) = (0.0D0, 0.0D0)
        if (nspin == 4) psic_nc (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              psic_nc (dfftp%nl (igk_k(ig,ikk)), isp) = evc (ig+npwx*(isp-1), ib)  !FZ: added for spinors
            ENDDO !FZ: added for spinors
          ELSE !FZ: added for spinors
            psic (dfftp%nl (igk_k(ig,ikk))) = evc (ig, ib)
          ENDIF !FZ: added for spinors
        ENDDO
        IF (nspin == 4) THEN !FZ: added for spinors
          DO isp=1, nst  !FZ: added for spinors
            CALL invfft ('Rho', psic_nc(:,isp),dfftp)  !FZ: added for spinors
          ENDDO  !FZ: added for spinors
        ELSE  !FZ: added for spinors
          CALL invfft ('Rho', psic, dfftp)
        ENDIF !FZ: added for spinors
        DO ir = 1, dfftp%nnr
          IF (nspin == 4) THEN !FZ: added for spinors
            ! ZL: add for extra SOC-mag terms
            temp_psic_nc1 = psic_nc(ir, 1)
            temp_psic_nc2 = psic_nc(ir, 2)

            ! ZL: redo to include all terms
            !DO isp=1,nst  !FZ: added for spinors
            !  psic_nc (ir, isp) = psic_nc (ir, isp) * vxcr (ir, 1)  !FZ: addedfor spinors
            !ENDDO  !FZ: added for spinors

            ! ZL: include all terms
            psic_nc(ir, 1) = temp_psic_nc1 * vxcr (ir, 1) + temp_psic_nc1 * vxcr (ir, 4) + &
                             temp_psic_nc2 * ( vxcr (ir, 2) - (0.d0,1.d0) * (vxcr (ir, 3)) ) 
            psic_nc(ir, 2) = temp_psic_nc2 * vxcr (ir, 1) - temp_psic_nc2 * vxcr (ir, 4) + &
                             temp_psic_nc1 * ( vxcr (ir, 2) + (0.d0,1.d0) * (vxcr (ir, 3)) ) 

          ELSE  !FZ: added for spinors
            psic (ir) = psic (ir) * vxcr (ir, isk (ik))
          ENDIF  !FZ: added for spinors
        ENDDO
        IF (nspin == 4) THEN !FZ: added for spinors
          DO isp=1, nst  !FZ: added for spinors
            CALL fwfft ('Rho', psic_nc(:,isp),dfftp)  !FZ: added for spinors
          ENDDO  !FZ: added for spinors
        ELSE  !FZ: added for spinors
          CALL fwfft ('Rho', psic, dfftp)
        ENDIF  !FZ: added for spinors
        hpsi (:) = (0.0D0, 0.0D0)
        hpsi_nc (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp=1, nst  !FZ: added for spinors
              hpsi_nc (ig, isp) = psic_nc (dfftp%nl (igk_k(ig,ikk)), isp) !FZ:added for spinors
            ENDDO  !FZ: added for spinors
          ELSE  !FZ: added for spinors
            hpsi (ig) = psic (dfftp%nl (igk_k (ig,ikk)))
          ENDIF
        ENDDO
        psic (:) = (0.0D0, 0.0D0)
        if (nspin == 4) psic_nc (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              psic_nc (ig, isp) = evc(ig+npwx*(isp-1), ib)  !FZ: added for spinors
            ENDDO    !FZ: added for spinors
          ELSE   !FZ: added for spinors
            psic (ig) = evc (ig, ib)
          ENDIF  !FZ: added for spinors
        ENDDO
        !IF (exx_is_active ()) &
        !   CALL vexx (npwx, npw, 1, psic, hpsi)
        DO ib2 = offdiag_nmin, offdiag_nmax
          psic2 (:) = (0.0D0, 0.0D0)
          if (nspin == 4) psic_nc2 (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
          DO ig = 1, npw
            IF (nspin == 4) THEN !FZ: added for spinors
              DO isp = 1, nst !FZ: added for spinors
                psic_nc2 (ig, isp) = evc(ig+npwx*(isp-1), ib2)  !FZ: added for spinors
              ENDDO    !FZ: added for spinors
            ELSE   !FZ: added for spinors
              psic2 (ig) = evc (ig, ib2)
            ENDIF  !FZ: added for spinors
          ENDDO
          dummy = (0.0D0, 0.0D0)
          DO ig = 1, npw
            IF (nspin == 4) THEN !FZ: added for spinors
              DO isp=1,nst  !FZ: added for spinors
                dummy = dummy + conjg (psic_nc2 (ig, isp)) * hpsi_nc (ig, isp)
              ENDDO  !FZ: added for spinors
            ELSE !FZ: added for spinors
              dummy = dummy + conjg (psic2 (ig)) * hpsi (ig)
            ENDIF !FZ: added for spinors
          ENDDO
          dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
          CALL mp_sum (dummy, intra_bgrp_comm)
          mtxelo (ib2 - offdiag_nmin + 1, ib - offdiag_nmin &
            + 1, ik) = dummy
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  DEALLOCATE (vxcr)
  IF (noffdiag .GT. 0) DEALLOCATE (psic2)
  IF (noffdiag .GT. 0) DEALLOCATE (psic_nc2)  !FZ: added for spinors
  DEALLOCATE (hpsi)
  DEALLOCATE (hpsi_nc)  !FZ: added for spinors

  IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
  IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    IF (nspin == 4) THEN !FZ: added for spinors
      DO ik = 1, nkstot   !FZ: added for spinors
        WRITE (unit, 101) xk(:, ik), ndiag, &
          noffdiag **2
        IF (ndiag .GT. 0) THEN
          DO ib = diag_nmin, diag_nmax
            WRITE (unit, 102) 1, ib, mtxeld &
              (ib - diag_nmin + 1, ik)
          ENDDO
        ENDIF
        IF (noffdiag .GT. 0) THEN
          DO ib = offdiag_nmin, offdiag_nmax
            DO ib2 = offdiag_nmin, offdiag_nmax
              WRITE (unit, 103) 1, ib2, ib, mtxelo &
                (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, ik)
            ENDDO
          ENDDO
        ENDIF
      ENDDO !FZ: added for spinors
    ELSE !FZ: added for spinors
      DO ik = 1, nkstot / nspin
        WRITE (unit, 101) xk(:, ik), nspin * ndiag, &
          nspin * noffdiag **2
        IF (ndiag .GT. 0) THEN
          DO ib = diag_nmin, diag_nmax
            DO is = 1, nspin  !FZ: change order
              WRITE (unit, 102) is, ib, mtxeld &
                (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin)
            ENDDO  !FZ change order
          ENDDO
        ENDIF
        IF (noffdiag .GT. 0) THEN
          DO ib = offdiag_nmin, offdiag_nmax
            DO ib2 = offdiag_nmin, offdiag_nmax
              DO is = 1, nspin
                WRITE (unit, 103) is, ib2, ib, mtxelo &
                  (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                  ik + (is - 1) * nkstot / nspin)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF !FZ: added for spinors
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
  IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

  RETURN

  101 FORMAT (3F13.9, 2I8)
  102 FORMAT (2I8, 2F15.9)
  103 FORMAT (3I8, 2F15.9)

END SUBROUTINE write_vxc_g

!-------------------------------------------------------------------------------
!FZ: added for metaGGA

SUBROUTINE write_vxc_g_meta (output_file_name, diag_nmin, diag_nmax, &
  offdiag_nmin, offdiag_nmax, vxc_zero_rho_core)

  USE constants, ONLY : rytoev
  USE cell_base, ONLY : tpiba2, at, bg
  USE ener, ONLY : etxc, vtxc
  USE exx, ONLY : vexx    !FZ 
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  !USE funct, ONLY : exx_is_active, dft_is_meta, get_meta, dft_is_hybrid     !FZ: Added for metaGGA, Added dft_is_hybrid for hybrid functional
  USE xc_lib, ONLY : exx_is_active, xclib_dft_is   !FZ: for qe6.7
  USE gvect, ONLY : ngm, g
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
  USE lsda_mod, ONLY : nspin, isk, current_spin, lsda !FZ: added current_spin, lsda
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm, intra_pool_comm  !FZ: added intra_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE mp_global,  ONLY : mp_startup   !FZ: added for Hybrid functional
  USE scf, ONLY : rho, rho_core, rhog_core
  USE wavefunctions, ONLY : evc, psic
  USE wvfct, ONLY : npwx, nbnd

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  integer, intent (inout) :: diag_nmin
  integer, intent (inout) :: diag_nmax
  integer, intent (inout) :: offdiag_nmin      
  integer, intent (inout) :: offdiag_nmax      
  logical, intent (in) :: vxc_zero_rho_core

  integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2, ikk
  integer, external :: global_kpoint_index
  complex (DP) :: dummy
  complex (DP), allocatable :: mtxeld (:, :)
  complex (DP), allocatable :: mtxelo (:, :, :)
  real (DP), allocatable :: vxcr (:, :)
  real (DP), allocatable :: kedtaur (:, :)   !FZ: for metaGGA  
  complex (DP), allocatable :: psic2 (:)
  complex (DP), allocatable :: hpsi (:)

  ! ZL: Disable VXC for spinor case
  !     this is because for SOC-mag case, the matrix elements of VXC by BGW is
  !     not supported yet
  if (nspin .eq. 4) then
    call errore('write_vxcg', 'full spinors do not support VXC, nspin =', nspin)
  endif

  if(diag_nmin > diag_nmax) then
    call errore ( 'write_vxc_g', 'diag_nmin > diag_nmax', diag_nmin )
  endif
  IF (diag_nmin .LT. 1) diag_nmin = 1
  IF (diag_nmax .GT. nbnd) then
    write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
    diag_nmax = nbnd
  ENDIF
  ndiag = MAX (diag_nmax - diag_nmin + 1, 0)

  if(offdiag_nmin > offdiag_nmax) then
    call errore ( 'write_vxc_g', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )
  endif
  IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
  IF (offdiag_nmax .GT. nbnd)  then
    write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
    offdiag_nmax = nbnd
  ENDIF
  noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

  IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

  unit = 4

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  IF (ndiag .GT. 0) THEN
    ALLOCATE (mtxeld (ndiag, nkstot))
    mtxeld (:, :) = (0.0D0, 0.0D0)
  ENDIF
  IF (noffdiag .GT. 0) THEN
    ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
    mtxelo (:, :, :) = (0.0D0, 0.0D0)
  ENDIF

  ALLOCATE (vxcr (dfftp%nnr, nspin))
  IF (noffdiag .GT. 0) ALLOCATE (psic2 (dfftp%nnr))
  ALLOCATE (hpsi (dfftp%nnr))

  ALLOCATE (kedtaur (dfftp%nnr, nspin))   !FZ: for meta GGA
  vxcr (:, :) = 0.0D0
  kedtaur (:, :) = 0.0D0    !FZ: for metaGGA
  IF ( vxc_zero_rho_core ) THEN
    rho_core ( : ) = 0.0D0
    rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  ENDIF
  !
  !IF (dft_is_meta() .and. (get_meta() /= 4)) then         !FZ: for metaGGA
  IF (xclib_dft_is('meta')) then         !FZ: for qe6.7
     CALL v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, vxcr, kedtaur )    !FZ: for metaGGA
  ELSE                         
     CALL v_xc(rho, rho_core, rhog_core, etxc, vtxc, vxcr)
  ENDIF                       
  !

  DO ik = iks, ike
    ikk = ik - iks + 1
    npw = ngk ( ik - iks + 1 )

    IF ( lsda ) current_spin = isk(ik)  !FZ: for LSDA specify spin

    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    IF (ndiag .GT. 0) THEN
      DO ib = diag_nmin, diag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (dfftp%nl (igk_k(ig,ikk))) = evc (ig, ib)     
        ENDDO
        CALL invfft ('Rho', psic, dfftp)
        DO ir = 1, dfftp%nnr
          psic (ir) = psic (ir) * vxcr (ir, isk (ik))
        ENDDO
        CALL fwfft ('Rho', psic, dfftp)
        hpsi (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          hpsi (ig) = psic (dfftp%nl (igk_k(ig,ikk)))
        ENDDO
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (ig) = evc (ig, ib)
        ENDDO
        !IF (exx_is_active ()) &    !FZ:  
        !   CALL vexx (npwx, npw, 1, psic, hpsi)  !FZ: commented!
        dummy = (0.0D0, 0.0D0)
        DO ig = 1, npw
          dummy = dummy + conjg (psic (ig)) * hpsi (ig)
        ENDDO
        dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
        CALL mp_sum (dummy, intra_bgrp_comm)
        mtxeld (ib - diag_nmin + 1, ik) = dummy
      ENDDO
    ENDIF
    IF (noffdiag .GT. 0) THEN
      DO ib = offdiag_nmin, offdiag_nmax
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (dfftp%nl (igk_k(ig,ikk))) = evc (ig, ib)
        ENDDO
        CALL invfft ('Rho', psic, dfftp)
        DO ir = 1, dfftp%nnr
          psic (ir) = psic (ir) * vxcr (ir, isk (ik))
        ENDDO
        CALL fwfft ('Rho', psic, dfftp)
        hpsi (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          hpsi (ig) = psic (dfftp%nl (igk_k (ig,ikk)))
        ENDDO
        psic (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          psic (ig) = evc (ig, ib)
        ENDDO
        !IF (exx_is_active ()) &  !FZ: 
        !   CALL vexx (npwx, npw, 1, psic, hpsi)   !FZ: commented 
        DO ib2 = offdiag_nmin, offdiag_nmax
          psic2 (:) = (0.0D0, 0.0D0)
          DO ig = 1, npw
            psic2 (ig) = evc (ig, ib2)
          ENDDO
          dummy = (0.0D0, 0.0D0)
          DO ig = 1, npw
            dummy = dummy + conjg (psic2 (ig)) * hpsi (ig)
          ENDDO
          dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
          CALL mp_sum (dummy, intra_bgrp_comm)
          mtxelo (ib2 - offdiag_nmin + 1, ib - offdiag_nmin &
            + 1, ik) = dummy
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  DEALLOCATE (vxcr)
  DEALLOCATE (kedtaur)      !FZ: for metaGGA
  IF (noffdiag .GT. 0) DEALLOCATE (psic2)
  DEALLOCATE (hpsi)

  IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
  IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'formatted', status = 'replace')
    DO ik = 1, nkstot / nspin
      WRITE (unit, 101) xk(:, ik), nspin * ndiag, &
        nspin * noffdiag **2
      IF (ndiag .GT. 0) THEN
        DO ib = diag_nmin, diag_nmax
          DO is = 1, nspin  !FZ: change order
            WRITE (unit, 102) is, ib, mtxeld &
              (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin)
          ENDDO
        ENDDO
      ENDIF
      IF (noffdiag .GT. 0) THEN
        DO ib = offdiag_nmin, offdiag_nmax
          DO ib2 = offdiag_nmin, offdiag_nmax
            DO is = 1, nspin  !FZ: change order
              WRITE (unit, 103) is, ib2, ib, mtxelo &
                (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                ik + (is - 1) * nkstot / nspin)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
  IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

  RETURN

  101 FORMAT (3F13.9, 2I8)
  102 FORMAT (2I8, 2F15.9)
  103 FORMAT (3I8, 2F15.9)

END SUBROUTINE write_vxc_g_meta
!FZ: end modify
!-------------------------------------------------------------------------------

!FZ: added for kih energy, hybrid functionals
SUBROUTINE write_kih (kih_file_name, vxc_hybrid_file_name, diag_nmin, diag_nmax, &     
  offdiag_nmin, offdiag_nmax, kih_flag, vxc_hybrid_flag)                     

  USE constants, ONLY : rytoev, eps12  ! ZL: add eps12
  USE cell_base, ONLY : tpiba2, at, bg
  USE ener, ONLY : etxc, vtxc, ehart !FZ:  added ehart
  USE exx, ONLY : vexx
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  !USE funct, ONLY : exx_is_active, dft_is_meta, get_meta     !FZ: Added for metaGGA commented in v_h
  USE xc_lib, ONLY : exx_is_active, xclib_dft_is     !FZ: for qe6.7
  USE gvect, ONLY : ngm, g
  USE gvecs, ONLY : doublegrid   !FZ: added for calculation of vrs
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
  !USE lsda_mod, ONLY : nspin, isk
  USE lsda_mod,    ONLY : current_spin, lsda, isk   !FZ: for LSDA
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : rho, vrs, vltot, v_of_0, v, kedtau, rho_core, rhog_core  !FZ: for vhartree 
  USE wavefunctions, ONLY : evc, psic, psic_nc !FZ: added psic_nc for spinors
  USE wvfct, ONLY : npwx, nbnd, g2kin, et  !FZ: added g2kin
  USE g_psi_mod,            ONLY : h_diag, s_diag 
  USE noncollin_module, ONLY: noncolin, npol  
  !USE uspp,                 ONLY : vkb, nkb, okvan, deeq, qq_at, qq_so, deeq_nc, indv_ijkb0 !FZ: 
  USE uspp,                 ONLY : vkb, nkb, okvan, deeq, qq_at, qq_so, deeq_nc, ofsbeta !FZ: 
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp 
  USE wvfct, ONLY: npwx, current_k !FZ: 
  USE uspp_param, ONLY: upf, nh    !FZ: 
  USE becmod,   ONLY : bec_type, becp, calbec, &  
                         allocate_bec_type, deallocate_bec_type 
  USE fft_base,      ONLY : dffts !FZ: test
  USE uspp_init,            ONLY : init_us_2
  USE symm_base, ONLY : nsym ! ZL: add
  USE noncollin_module, ONLY : magtot_nc ! ZL: add
  USE lsda_mod, ONLY: nspin, starting_magnetization

  IMPLICIT NONE

  character (len = 256), intent (in) :: kih_file_name
  character (len = 256), intent (in) :: vxc_hybrid_file_name
  integer, intent (inout) :: diag_nmin
  integer, intent (inout) :: diag_nmax
  integer, intent (inout) :: offdiag_nmin
  integer, intent (inout) :: offdiag_nmax
  logical :: kih_flag, vxc_hybrid_flag

  integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2, ikk
  integer, external :: global_kpoint_index
  complex (DP) :: dummy
  complex (DP), allocatable :: mtxeld (:, :)
  complex (DP), allocatable :: mtxelo (:, :, :)
  real (DP), allocatable :: vxcr (:, :)
  real (DP), allocatable :: vxcr_2 (:, :)    
  real (DP), allocatable :: vxcr_nlcc (:, :) 
  real (DP), allocatable :: v_har (:, :)      
  real (DP), allocatable :: kedtaur (:, :)   !FZ: for metaGGA  
  complex (DP), allocatable :: psic2 (:)
  complex (DP), allocatable :: psic_nc2 (:,:) !FZ: added for spinors
  complex (DP), allocatable :: psic_temp (:) 
  complex (DP), allocatable :: hpsi (:)
  complex (DP), allocatable :: hpsi_temp (:)  
  complex (DP), allocatable :: hpsi_nc (:,:)  !FZ: added for spinors
  complex (DP), allocatable :: psic_nc_temp (:,:)  !FZ: added for spinors
  complex (DP), allocatable :: hpsi_nc_temp (:,:)  !FZ: added for spinors
  integer :: isp, ns, nst, nsf !FZ: added for spinors
  complex (DP) :: temp_psic_nc1, temp_psic_nc2  !FZ: added for spinors
  REAL(DP)  ::  charge    
  character(len=20) :: Ctemp
  INTEGER :: ierr   
  COMPLEX(DP), allocatable :: hpsinl(:,:)  
  COMPLEX(DP), allocatable :: psinl(:,:)  

  ! ZL: soc-mag must not use symmetries
  if ( sum(abs(starting_magnetization)) .gt. eps12 ) then
    if (nsym .ne. 1) then
      call errore('write_vxc_r', 'Magnetization detected for full spinor. No symmetry is allowed, nsym =', nsym)
    endif
  endif

  if(diag_nmin > diag_nmax) then
    call errore ( 'write_IHK', 'diag_nmin > diag_nmax', diag_nmin )  
  endif
  IF (diag_nmin .LT. 1) diag_nmin = 1
  IF (diag_nmax .GT. nbnd) then
    write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
    diag_nmax = nbnd
  ENDIF
  ndiag = MAX (diag_nmax - diag_nmin + 1, 0)
  
  if(offdiag_nmin > offdiag_nmax) then
    call errore ( 'write_IHK', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )
  endif
  IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
  IF (offdiag_nmax .GT. nbnd)  then
    write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
    offdiag_nmax = nbnd
  ENDIF
  noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

  IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

  unit = 4

  call set_spin(ns, nst, nsf, nspin) !FZ: added for spinor

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  IF (ndiag .GT. 0) THEN
    ALLOCATE (mtxeld (ndiag, nkstot))
    mtxeld (:, :) = (0.0D0, 0.0D0)
  ENDIF
  IF (noffdiag .GT. 0) THEN
    ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
    mtxelo (:, :, :) = (0.0D0, 0.0D0)
  ENDIF

  ALLOCATE (vxcr (dfftp%nnr, nspin))
  ALLOCATE (vxcr_2 (dfftp%nnr, nspin))   
  ALLOCATE (vxcr_nlcc (dfftp%nnr, nspin)) 
  IF (noffdiag .GT. 0) ALLOCATE (psic2 (dfftp%nnr))
  IF (noffdiag .GT. 0) ALLOCATE (psic_nc2 (dfftp%nnr, 2))   !FZ: added for spinor
  ALLOCATE (v_har (dfftp%nnr, nspin))    
  ALLOCATE (hpsi (dfftp%nnr))
  ALLOCATE (hpsi_nc (dfftp%nnr, 2))  !FZ: added for spinors
  ALLOCATE (hpsi_temp (dfftp%nnr))
  ALLOCATE (psic_temp (dfftp%nnr))
  ALLOCATE (psic_nc_temp (dfftp%nnr, 2))  !FZ: added for spinors
  ALLOCATE (hpsi_nc_temp (dfftp%nnr, 2))  !FZ: added for spinors
  ALLOCATE (hpsinl (npwx*npol, nbnd))
  ALLOCATE (psinl (npwx*npol, nbnd))
  ALLOCATE (kedtaur (dfftp%nnr, nspin))   !FZ: for meta GGA
  CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm )   

  vxcr (:, :) = 0.0D0
  vxcr_2 (:, :) = 0.0D0   
  vxcr_nlcc (:, :) = 0.0D0 
  hpsinl(:,:) = (0.0D0, 0.0D0) 
  psinl(:,:) = (0.0D0, 0.0D0) 
  v_har (:, :) = 0.0D0       
  kedtaur (:, :) = 0.0D0    
  ehart  = 0.D0
  charge = 0.D0

  rho_core ( : ) = 0.0D0
  rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  !IF (dft_is_meta() .and. (get_meta() /= 4)) then         !FZ: for metaGGA
  IF (xclib_dft_is('meta')) then         !FZ: for qe6.7
     CALL v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, vxcr, kedtaur )    !FZ: for metaGGA
  ELSE                      
     CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxcr)
  ENDIF
  CALL set_rhoc()    
  !IF (dft_is_meta() .and. (get_meta() /= 4)) then         !FZ: for metaGGA
  IF (xclib_dft_is('meta')) then         !FZ: for qe6.7
     CALL v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, vxcr_2, kedtaur )    !FZ: for metaGGA
  ELSE                       
     CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxcr_2)
  ENDIF

  CALL v_h (rho%of_g(:,1), ehart, charge, v_har)   
  ALLOCATE( h_diag( npwx, npol ), STAT=ierr )    
  IF( ierr /= 0 ) &    
     CALL errore( ' diag_bands ', ' cannot allocate h_diag ', ABS(ierr) )   
  ALLOCATE( s_diag( npwx, npol ), STAT=ierr )   
  IF( ierr /= 0 ) &     
     CALL errore( ' diag_bands ', ' cannot allocate s_diag ', ABS(ierr) )    
  h_diag (:,:) = 0.0D0  
  s_diag (:,:) = 0.0D0  
  CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, &    
                         nspin, doublegrid )    
  
  DO ik = iks, ike
    ikk = ik - iks + 1
    npw = ngk ( ik - iks + 1 )

    IF ( lsda ) current_spin = isk(ik)  

    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    CALL threaded_memcpy(psinl, evc, nbnd*npol*npwx*2)
    vkb (:,:) = 0.0D0  
    IF ( nkb > 0 ) CALL init_us_2( ngk(ik), igk_k(1,ik), xk(1,ik), vkb ) 
    h_diag (:,:) = 0.0D0  
    s_diag (:,:) = 0.0D0  
    hpsinl (:,:) = 0.0D0  
    CALL usnldiag( npw, h_diag, s_diag )   
    g2kin(:) = 0.0D0  
    call g2_kin( ik )  
    CALL calbec ( npw, vkb, psinl, becp, nbnd )   
    CALL add_vuspsi( npwx, npw, nbnd, hpsinl )   
    current_k = ik     
    !if (dft_is_meta()) call h_psi_meta (npwx, npw, nbnd, psinl, hpsinl)  
    if (xclib_dft_is('meta')) call h_psi_meta (npwx, npw, nbnd, psinl, hpsinl) !FZ: for qe6.7 
    IF (ndiag .GT. 0) THEN
      DO ib = diag_nmin, diag_nmax
        psic (:) = (0.0D0, 0.0D0)
        hpsi (:) = (0.0D0, 0.0D0)  
        psic_temp (:) = (0.0D0, 0.0D0)
        hpsi_temp (:) = (0.0D0, 0.0D0)  
        if (nspin == 4) psic_nc (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        hpsi_nc (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        psic_nc_temp (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        hpsi_nc_temp (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              psic_nc (igk_k (ig, ikk), isp) = evc(ig+npwx*(isp-1), ib)  !FZ: added for spinors
            ENDDO    !FZ: added for spinors
          ELSE
            psic (igk_k(ig,ikk)) = evc (ig, ib)     
          ENDIF
        ENDDO
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              psic_nc_temp (dfftp%nl (igk_k (ig, ikk)), isp) = evc(ig+npwx*(isp-1), ib)  !FZ: added for spinors
            ENDDO    !FZ: added for spinors
          ELSE
            psic_temp (dfftp%nl (igk_k(ig,ikk))) = evc (ig, ib)  
          ENDIF  !FZ: added for spinors
        ENDDO
        IF (nspin == 4) THEN !FZ: added for spinors
          DO isp=1, nst  !FZ: added for spinors
            CALL invfft ('Rho', psic_nc_temp(:,isp),dfftp)  !FZ: added for spinors
          ENDDO  !FZ: added for spinors
        ELSE  !FZ: added for spinors
          CALL invfft ('Rho', psic_temp, dfftp)
        ENDIF  !FZ: added for spinors
        DO ir = 1, dfftp%nnr
          IF (nspin == 4) THEN !FZ: added for spinors
            temp_psic_nc1 = psic_nc_temp (ir, 1)
            temp_psic_nc2 = psic_nc_temp (ir, 2)

            ! ZL: vxc terms are not needed here
            !psic_nc_temp (ir, 1) = temp_psic_nc1 * (vltot (ir) + v_har (ir, 1) + vxcr_2 (ir, 1) - vxcr (ir, 1) + vxcr_2 (ir, 4))  &  !FZ: added for spinors
            !            + temp_psic_nc2 * ( (vxcr_2 (ir, 2) ) - (0.d0,1.d0) * (vxcr_2 (ir, 3)) ) 
            !psic_nc_temp (ir, 2) = temp_psic_nc2 * (vltot (ir) + v_har (ir, 1) + vxcr_2 (ir, 1) - vxcr (ir, 1) - vxcr_2 (ir, 4))  &  !FZ: added for spinors
            !            + temp_psic_nc1 * ( (vxcr_2 (ir, 2)) + (0.d0,1.d0) * (vxcr_2 (ir, 3)) ) 
            psic_nc_temp (ir, 1) = temp_psic_nc1 * (vltot (ir) + v_har (ir, 1))
            psic_nc_temp (ir, 2) = temp_psic_nc2 * (vltot (ir) + v_har (ir, 1))

          ELSE  !FZ: added for spinors
            ! ZL: no vxc
            !psic_temp (ir) = psic_temp (ir) * (vltot (ir) + v_har (ir, isk (ik)) + vxcr_2 (ir, isk (ik)) - vxcr (ir, isk (ik)))   
            psic_temp (ir) = psic_temp (ir) * (vltot (ir) + v_har (ir, isk (ik)))   
          ENDIF  !FZ: added for spinors
        ENDDO
        IF (nspin == 4) THEN !FZ: added for spinors
          DO isp=1, nst  !FZ: added for spinors
            CALL fwfft ('Rho', psic_nc_temp(:,isp),dfftp)  !FZ: added for spinors
          ENDDO  !FZ: added for spinors
        ELSE  !FZ: added for spinors
          CALL fwfft ('Rho', psic_temp, dfftp)
        ENDIF  !FZ: added for spinors
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp=1, nst  !FZ: added for spinors
              hpsi_nc_temp (ig, isp) = psic_nc_temp (dfftp%nl (igk_k(ig,ikk)), isp) !FZ:added for spinors
            ENDDO  !FZ: added for spinors
          ELSE  !FZ: added for spinors
            hpsi_temp (ig) = psic_temp (dfftp%nl (igk_k(ig,ikk)))
          ENDIF  !FZ: added for spinors
        ENDDO
        DO ig = 1, npw                                
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp=1, nst  !FZ: added for spinors
              hpsi_nc (ig, isp) = g2kin(ig) * psic_nc (igk_k(ig,ikk), isp) !FZ:added for spinors
            ENDDO  !FZ: added for spinors
          ELSE  !FZ: added for spinors
            hpsi (ig) = g2kin(ig) * psic (igk_k(ig,ikk)) 
          ENDIF  !FZ: added for spinors
        ENDDO                                           
        psic (:) = (0.0D0, 0.0D0)
        if (nspin == 4) psic_nc (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              psic_nc (ig, isp) = evc(ig+npwx*(isp-1), ib)  !FZ: added for spinors
            ENDDO    !FZ: added for spinors
          ELSE
            psic (ig) = evc (ig, ib)     
          ENDIF  !FZ: added for spinors
        ENDDO
        !IF (exx_is_active ()) &
        !   CALL vexx (npwx, npw, 1, psic, hpsi)
        dummy = (0.0D0, 0.0D0)
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              dummy = dummy + conjg (psic_nc(ig, isp)) * (hpsinl(ig+npwx*(isp-1), ib) + hpsi_nc(ig, isp) + hpsi_nc_temp(ig, isp))   !FZ: added for spinors
            ENDDO    !FZ: added for spinors
          ELSE
            dummy = dummy + conjg (psic (ig)) * (hpsinl (ig, ib) + hpsi(ig) + hpsi_temp (ig))    
          ENDIF  !FZ: added for spinors
        ENDDO
        dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
        CALL mp_sum (dummy, intra_bgrp_comm)
        mtxeld (ib - diag_nmin + 1, ik) = dummy
      ENDDO
    ENDIF
    IF (noffdiag .GT. 0) THEN
      DO ib = offdiag_nmin, offdiag_nmax
        psic (:) = (0.0D0, 0.0D0)
        hpsi (:) = (0.0D0, 0.0D0)  
        if (nspin == 4) psic_nc (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        hpsi_nc (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        if (nspin == 4) psic_nc_temp (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        hpsi_nc_temp (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              psic_nc (ig, isp) = evc(ig+npwx*(isp-1), ib)  !FZ: added for spinors
            ENDDO    !FZ: added for spinors
          ELSE
            psic (ig) = evc (ig, ib)  
          ENDIF
        ENDDO
        psic_temp (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              psic_nc_temp (dfftp%nl (igk_k (ig, ikk)), isp) = evc(ig+npwx*(isp-1), ib)  !FZ: added for spinors
            ENDDO    !FZ: added for spinors
          ELSE
            psic_temp (dfftp%nl (igk_k(ig,ikk))) = evc (ig, ib)
          ENDIF  !FZ: added for spinors
        ENDDO
        IF (nspin == 4) THEN !FZ: added for spinors
          DO isp=1, nst  !FZ: added for spinors
            CALL invfft ('Rho', psic_nc_temp(:,isp),dfftp)  !FZ: added for spinors
          ENDDO  !FZ: added for spinors
        ELSE  !FZ: added for spinors
          CALL invfft ('Rho', psic_temp, dfftp)
        ENDIF  !FZ: added for spinors
        DO ir = 1, dfftp%nnr
          IF (nspin == 4) THEN !FZ: added for spinors
            temp_psic_nc1 = psic_nc_temp (ir, 1)
            temp_psic_nc2 = psic_nc_temp (ir, 2)

            ! ZL: vxc terms are not needed here
            !psic_nc_temp (ir, 1) = temp_psic_nc1 * (vltot (ir) + v_har (ir, 1) + vxcr_2 (ir, 1) - vxcr (ir, 1) + vxcr_2 (ir, 4))  &  !FZ: added for spinors
            !            + temp_psic_nc2 * ( (vxcr_2 (ir, 2) ) - (0.d0,1.d0) * (vxcr_2 (ir, 3)) ) 
            !psic_nc_temp (ir, 2) = temp_psic_nc2 * (vltot (ir) + v_har (ir, 1) + vxcr_2 (ir, 1) - vxcr (ir, 1) - vxcr_2 (ir, 4))  &  !FZ: added for spinors
            !            + temp_psic_nc1 * ( (vxcr_2 (ir, 2)) + (0.d0,1.d0) * (vxcr_2 (ir, 3)) ) 
            psic_nc_temp (ir, 1) = temp_psic_nc1 * (vltot (ir) + v_har (ir, 1))
            psic_nc_temp (ir, 2) = temp_psic_nc2 * (vltot (ir) + v_har (ir, 1))

          ELSE  !FZ: added for spinors
            ! ZL: no vxc
            !psic_temp (ir) = psic_temp (ir) * (vltot (ir) + v_har (ir, isk (ik)) + vxcr_2 (ir, isk (ik)) - vxcr (ir, isk (ik)))
            psic_temp (ir) = psic_temp (ir) * (vltot (ir) + v_har (ir, isk (ik)))
          ENDIF  !FZ: added for spinors
        ENDDO
        IF (nspin == 4) THEN !FZ: added for spinors
          DO isp=1, nst  !FZ: added for spinors
            CALL fwfft ('Rho', psic_nc_temp(:,isp),dfftp)  !FZ: added for spinors
          ENDDO  !FZ: added for spinors
        ELSE  !FZ: added for spinors
          CALL fwfft ('Rho', psic_temp, dfftp)
        ENDIF  !FZ: added for spinors
        hpsi_temp (:) = (0.0D0, 0.0D0)
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp=1, nst  !FZ: added for spinors
              hpsi_nc_temp (ig, isp) = psic_nc_temp (dfftp%nl (igk_k(ig,ikk)), isp) !FZ:added for spinors
            ENDDO  !FZ: added for spinors
          ELSE  !FZ: added for spinors
            hpsi_temp (ig) = psic_temp (dfftp%nl (igk_k (ig,ikk)))
          ENDIF  !FZ: added for spinors
        ENDDO
        DO ig = 1, npw                        
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp=1, nst  !FZ: added for spinors
              hpsi_nc (ig, isp) = g2kin(ig) * psic_nc (ig, isp) !FZ:added for spinors
            ENDDO  !FZ: added for spinors
          ELSE  !FZ: added for spinors
            hpsi (ig) = g2kin(ig) * psic (ig)   
          ENDIF  !FZ: added for spinors
        ENDDO
        psic (:) = (0.0D0, 0.0D0)
        if (nspin == 4) psic_nc (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
        DO ig = 1, npw
          IF (nspin == 4) THEN !FZ: added for spinors
            DO isp = 1, nst !FZ: added for spinors
              psic_nc (ig, isp) = evc(ig+npwx*(isp-1), ib)  !FZ: added for spinors
            ENDDO    !FZ: added for spinors
          ELSE
            psic (ig) = evc (ig, ib)
          ENDIF  !FZ: added for spinors
        ENDDO
        !IF (exx_is_active ()) &
        !   CALL vexx (npwx, npw, 1, psic, hpsi)
        DO ib2 = offdiag_nmin, offdiag_nmax
          psic2 (:) = (0.0D0, 0.0D0)
          if (nspin == 4) psic_nc2 (:,:) = (0.0D0, 0.0D0)! FZ: added for spinors
          DO ig = 1, npw
            IF (nspin == 4) THEN !FZ: added for spinors
              DO isp = 1, nst !FZ: added for spinors
                psic_nc2 (ig, isp) = evc(ig+npwx*(isp-1), ib2)  !FZ: added for spinors
              ENDDO    !FZ: added for spinors
            ELSE   !FZ: added for spinors
              psic2 (ig) = evc (ig, ib2)
            ENDIF  !FZ: added for spinors
          ENDDO
          dummy = (0.0D0, 0.0D0)
          DO ig = 1, npw
            IF (nspin == 4) THEN !FZ: added for spinors
              DO isp=1,nst  !FZ: added for spinors
                dummy = dummy + conjg (psic_nc2 (ig, isp)) * (hpsinl (ig+npwx*(isp-1), ib) &
                        + hpsi_nc (ig, isp) + hpsi_nc_temp(ig, isp))
              ENDDO  !FZ: added for spinors
            ELSE !FZ: added for spinors
              dummy = dummy + conjg (psic2 (ig)) * (hpsinl (ig, ib) + hpsi(ig) + hpsi_temp (ig)) 
            ENDIF !FZ: added for spinors
          ENDDO
          dummy = dummy * CMPLX (rytoev, 0.0D0, KIND=dp)
          CALL mp_sum (dummy, intra_bgrp_comm)
          mtxelo (ib2 - offdiag_nmin + 1, ib - offdiag_nmin &
            + 1, ik) = dummy
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  DEALLOCATE (vxcr)
  DEALLOCATE (vxcr_2)
  DEALLOCATE (vxcr_nlcc)
  DEALLOCATE (v_har)          !FZ: for metaGGA 
  DEALLOCATE (kedtaur)      !FZ: for metaGGA
  DEALLOCATE (hpsi)
  DEALLOCATE (hpsinl) 
  DEALLOCATE (psinl)  
  DEALLOCATE (h_diag) 
  DEALLOCATE (s_diag) 
  IF (noffdiag .GT. 0) DEALLOCATE (psic2)
  CALL deallocate_bec_type ( becp )   
  DEALLOCATE (hpsi_temp)
  DEALLOCATE (psic_temp)

  IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
  IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (kih_flag) THEN
    IF (ionode) THEN
      OPEN (unit = unit, file = TRIM (kih_file_name), &
        form = 'formatted', status = 'replace')
      IF (nspin == 4) THEN !FZ: added for spinors
        DO ik = 1, nkstot !FZ: added for spinors
          WRITE (unit, 101) xk(:, ik), ndiag, &  
            noffdiag **2     
          IF (ndiag .GT. 0) THEN
            DO ib = diag_nmin, diag_nmax
              WRITE (unit, 102) 1, ib, mtxeld &
                (ib - diag_nmin + 1, ik)
            ENDDO
          ENDIF
          IF (noffdiag .GT. 0) THEN
            DO ib = offdiag_nmin, offdiag_nmax
              DO ib2 = offdiag_nmin, offdiag_nmax
                WRITE (unit, 103) 1, ib2, ib, mtxelo &
                  (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, ik)
              ENDDO
            ENDDO
          ENDIF
        ENDDO  !FZ: added for spinors
      ELSE !FZ: added for spinors
        DO ik = 1, nkstot / nspin
          WRITE (unit, 101) xk(:, ik), nspin * ndiag, &  
            nspin * noffdiag **2     
          IF (ndiag .GT. 0) THEN
            DO ib = diag_nmin, diag_nmax
              DO is = 1, nspin  !FZ: change order
                WRITE (unit, 102) is, ib, mtxeld &
                  (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin)
              ENDDO
            ENDDO
          ENDIF
          IF (noffdiag .GT. 0) THEN
            DO ib = offdiag_nmin, offdiag_nmax
              DO ib2 = offdiag_nmin, offdiag_nmax
                DO is = 1, nspin  !FZ: change order
                  WRITE (unit, 103) is, ib2, ib, mtxelo &
                    (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                    ik + (is - 1) * nkstot / nspin)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF !FZ: added for spinors
      CLOSE (unit = unit, status = 'keep')
    ENDIF
  ENDIF
  IF (vxc_hybrid_flag) THEN
    IF (ionode) THEN
      OPEN (unit = unit, file = TRIM (vxc_hybrid_file_name), &
        form = 'formatted', status = 'replace')
      IF (nspin == 4) THEN !FZ: added for spinors
        DO ik = 1, nkstot   !FZ: added for spinors
          WRITE (unit, 101) xk(:, ik), ndiag, &   
            noffdiag **2     
          IF (ndiag .GT. 0) THEN
            DO ib = diag_nmin, diag_nmax
              WRITE (unit, 102) 1, ib, et(ib, ik) * CMPLX (rytoev, 0.0D0, KIND=dp) - mtxeld &
                (ib - diag_nmin + 1, ik)
            ENDDO
          ENDIF
          IF (noffdiag .GT. 0) THEN
            DO ib = offdiag_nmin, offdiag_nmax
              DO ib2 = offdiag_nmin, offdiag_nmax
                IF (ib .EQ. ib2) THEN   !FZ: 
                  WRITE (unit, 103) 1, ib2, ib, et(ib, ik) * CMPLX (rytoev, 0.0D0, KIND=dp) - mtxelo &
                    (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, ik)
                ELSE        !FZ: 
                   WRITE (unit, 103) 1, ib2, ib, - mtxelo &
                    (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, ik)
                ENDIF       !FZ: 
              ENDDO
            ENDDO
          ENDIF
        ENDDO  !FZ: added for spinors
      ELSE !FZ: added for spinors
        DO ik = 1, nkstot / nspin
          WRITE (unit, 101) xk(:, ik), nspin * ndiag, &   
            nspin * noffdiag **2     
          IF (ndiag .GT. 0) THEN
            DO ib = diag_nmin, diag_nmax
              DO is = 1, nspin  !FZ: change order
                WRITE (unit, 102) is, ib, et(ib, ik + (is - 1) * nkstot / nspin) * CMPLX (rytoev, 0.0D0, KIND=dp) - mtxeld &
                  (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin)
              ENDDO
            ENDDO
          ENDIF
          IF (noffdiag .GT. 0) THEN
            DO ib = offdiag_nmin, offdiag_nmax
              DO ib2 = offdiag_nmin, offdiag_nmax
                DO is = 1, nspin  !FZ: change order
                  IF (ib .EQ. ib2) THEN   !FZ: 
                    WRITE (unit, 103) is, ib2, ib, et(ib, ik + (is - 1) * nkstot / nspin) * CMPLX (rytoev, 0.0D0, KIND=dp) & 
                      - mtxelo(ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                      ik + (is - 1) * nkstot / nspin)
                  ELSE        !FZ: 
                     WRITE (unit, 103) is, ib2, ib, - mtxelo &
                      (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                      ik + (is - 1) * nkstot / nspin)
                  ENDIF       !FZ: 
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF !FZ: added for spinors
      CLOSE (unit = unit, status = 'keep')
    ENDIF
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
  IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

  RETURN

  101 FORMAT (3F13.9, 2I8)   
  102 FORMAT (2I8, 2F15.9)
  103 FORMAT (3I8, 2F15.9)

END SUBROUTINE write_kih
!FZ: end modify
!-------------------------------------------------------------------------------

SUBROUTINE write_vscg ( output_file_name, real_or_complex, symm_type )

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  USE constants, ONLY : pi, tpi, eps6
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, mill, ecutrho
  USE io_global, ONLY : ionode
  USE ions_base, ONLY : nat, atm, ityp, tau 
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_sum
  USE mp_bands, ONLY : intra_bgrp_comm
  USE scf, ONLY : vltot, v
  USE symm_base, ONLY : s, ft, nsym
  USE wavefunctions, ONLY : psic
  USE matrix_inversion

  IMPLICIT NONE

  character ( len = 256 ), intent (in) :: output_file_name
  integer, intent (in) :: real_or_complex
  character ( len = 9 ), intent (in) :: symm_type

  character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  integer :: unit, id, is, ir, ig, i, j, k, ierr
  integer :: nd, ns, nr, ng_l, ng_g, nst, nsf  ! ZL: add nst, nsf
  integer :: ntran, cell_symmetry, nrecord
  real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: g_g ( :, : )
  real (DP), allocatable :: vscr_g ( :, : )
  complex (DP), allocatable :: vscg_g ( :, : )

  INTEGER, EXTERNAL :: atomic_number

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  ! this is supposed to be VSC-Real/Complex but BGW wfn_rho_vxc IO
  ! does not recognize VSC header so we are using VXC instead
  IF ( real_or_complex .EQ. 1 ) THEN
    WRITE ( stitle, '("VSC-Real",24X)' )
  ELSE
    WRITE ( stitle, '("VSC-Complex",21X)' )
  ENDIF

  unit = 4
  nrecord = 1
  nd = 3

  ! ZL: use set_spin for spinors
  !ns = nspin
  call set_spin(ns, nst, nsf, nspin)
  nr = dfftp%nnr
  ng_l = ngm
  ng_g = ngm_g

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_vscg', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2 )
    t1 ( 1 ) = ft ( 1, i )
    t1 ( 2 ) = ft ( 2, i )
    t1 ( 3 ) = ft ( 3, i )
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true., translation )

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  ALLOCATE ( g_g ( nd, ng_g ) )
  ALLOCATE ( vscr_g ( ng_g, ns ) )
  ALLOCATE ( vscg_g ( ng_g, ns ) )

  DO ig = 1, ng_g
    DO id = 1, nd
      g_g ( id, ig ) = 0
    ENDDO
  ENDDO
  DO is = 1, ns
    DO ig = 1, ng_g
      vscg_g ( ig, is ) = ( 0.0D0, 0.0D0 )
    ENDDO
  ENDDO

  DO ig = 1, ng_l
    g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO
  vscr_g ( :, : ) = 0.0D0
  DO is = 1, ns
    DO ir = 1, nr
      psic ( ir ) = CMPLX ( v%of_r ( ir, is ) + vltot ( ir ), 0.0D0, KIND=dp )
    ENDDO
    CALL fwfft ( 'Rho', psic, dfftp )
    DO ig = 1, ng_l
      vscg_g ( ig_l2g ( ig ), is ) = psic ( dfftp%nl ( ig ) )
    ENDDO
  ENDDO

  CALL mp_sum ( g_g, intra_bgrp_comm )
  CALL mp_sum ( vscg_g, intra_bgrp_comm )

  IF ( ionode ) THEN
    OPEN ( unit = unit, file = TRIM ( output_file_name ), &
      form = 'unformatted', status = 'replace' )
    WRITE ( unit ) stitle, sdate, stime
    ! ZL: change ns to nsf
    !WRITE ( unit ) ns, ng_g, ntran, cell_symmetry, nat, ecutrho
    WRITE ( unit ) nsf, ng_g, ntran, cell_symmetry, nat, ecutrho
    WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
    WRITE ( unit ) nrecord
    WRITE ( unit ) ng_g
    IF ( real_or_complex .EQ. 1 ) THEN
      WRITE ( unit ) ( ( dble ( vscg_g ( ig, is ) ), &
        ig = 1, ng_g ), is = 1, ns )
    ELSE
      WRITE ( unit ) ( ( vscg_g ( ig, is ), &
        ig = 1, ng_g ), is = 1, ns )
    ENDIF
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( vscg_g )
  DEALLOCATE ( vscr_g )
  DEALLOCATE ( g_g )

  RETURN

END SUBROUTINE write_vscg

!-------------------------------------------------------------------------------

SUBROUTINE write_vkbg (output_file_name, symm_type, wfng_kgrid, &
  wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, wfng_dk3)

  USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  USE constants, ONLY : pi, tpi, eps6
  USE fft_base, ONLY : dfftp
  USE gvect, ONLY : ngm, ngm_g, ig_l2g, g, mill, ecutrho
  USE io_global, ONLY : ionode, ionode_id
  USE ions_base, ONLY : nat, atm, ityp, tau, nsp
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, wk, ngk, nks, nkstot, igk_k
  USE lsda_mod, ONLY : nspin, isk
  USE mp, ONLY : mp_sum, mp_max, mp_get, mp_barrier
  USE mp_world, ONLY : mpime, nproc, world_comm
  USE mp_bands, ONLY : intra_bgrp_comm, nbgrp
  USE mp_pools, ONLY : me_pool, root_pool, npool, nproc_pool, &
    intra_pool_comm, inter_pool_comm
  USE mp_wave, ONLY : mergewf
  USE start_k, ONLY : nk1, nk2, nk3, k1, k2, k3
  USE symm_base, ONLY : s, ft, nsym
  USE uspp, ONLY : nkb, vkb, deeq, deeq_nc
  USE uspp_param, ONLY : nhm, nh
  USE wvfct, ONLY : npwx
  USE gvecw, ONLY : ecutwfc
  USE matrix_inversion
  USE uspp_init,            ONLY : init_us_2

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  character ( len = 9 ), intent (in) :: symm_type
  logical, intent (in) :: wfng_kgrid
  integer, intent (in) :: wfng_nk1
  integer, intent (in) :: wfng_nk2
  integer, intent (in) :: wfng_nk3
  real (DP), intent (in) :: wfng_dk1
  real (DP), intent (in) :: wfng_dk2
  real (DP), intent (in) :: wfng_dk3

  character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  integer :: i, j, k, ierr, ik, is, ig, ikb, iat, isp, ih, jh, &
    unit, iks, ike, npw, npw_g, npwx_g, ngg, ipsour, &
    igwx, local_pw, id, nd, ntran, cell_symmetry, nrecord
  real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  integer :: ns, nst, nsf ! ZL: add ns, nst, nsf
  real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  integer, allocatable :: kmap ( : )
  integer, allocatable :: smap ( : )
  integer, allocatable :: gvec ( :, : )
  integer, allocatable :: ngk_g ( : )
  integer, allocatable :: igk_l2g ( :, : )
  integer, allocatable :: itmp ( : )
  integer, allocatable :: igwk ( : )
  integer, allocatable :: igwf_l2g ( : )
  integer, allocatable :: ipmask ( : )
  complex (DP), allocatable :: vkb_g ( : )

  INTEGER, EXTERNAL :: atomic_number, global_kpoint_index

  IF ( nkb == 0 ) RETURN

  CALL date_and_tim ( cdate, ctime )
  WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  ! BGW wfn_rho_vxc IO does not recognize VKB header so this file
  ! is read directly by SAPO code in BerkeleyGW
  WRITE ( stitle, '("VKB-Complex",21X)' )

  unit = 4
  nrecord = 1
  nd = 3

  ! ZL: use set_spin for spinors
  call set_spin(ns, nst, nsf, nspin)

  ! FZ: npool is currently disabled since it may have some issues
  if ( npool .gt. 1 ) then
    call errore('write_vkbg', 'npool is not used for VKB calculations, please do not set npool when doing VKB calculation', npool)
  endif

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1 

  ierr = 0
  IF ( ibrav .EQ. 0 ) THEN
    IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
      cell_symmetry = 0
    ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
      cell_symmetry = 1
    ELSE
      ierr = 1
    ENDIF
  ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
    cell_symmetry = 0
  ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
    cell_symmetry = 1
  ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
    cell_symmetry = 0
  ELSE
    ierr = 1
  ENDIF
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_vkbg', 'cell_symmetry', ierr )

  ntran = nsym
  DO i = 1, ntran
    DO j = 1, nd
      DO k = 1, nd
        r1 ( k, j ) = dble ( s ( k, j, i ) )
      ENDDO
    ENDDO
    CALL invmat ( 3, r1, r2 )
    t1 ( 1 ) = ft ( 1, i )
    t1 ( 2 ) = ft ( 2, i )
    t1 ( 3 ) = ft ( 3, i )
    DO j = 1, nd
      t2 ( j ) = 0.0D0
      DO k = 1, nd
        t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
      ENDDO
      IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
      IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
        t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    ENDDO
    DO j = 1, nd
      translation ( j, i ) = t2 ( j ) * tpi
    ENDDO
  ENDDO

  alat2 = alat ** 2
  recvol = 8.0D0 * pi**3 / omega

  DO i = 1, nd
    DO j = 1, nd
      adot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        adot ( j, i ) = adot ( j, i ) + &
          at ( k, j ) * at ( k, i ) * alat2
      ENDDO
    ENDDO
  ENDDO

  DO i = 1, nd
    DO j = 1, nd
      bdot ( j, i ) = 0.0D0
    ENDDO
  ENDDO
  DO i = 1, nd
    DO j = 1, nd
      DO k = 1, nd
        bdot ( j, i ) = bdot ( j, i ) + &
          bg ( k, j ) * bg ( k, i ) * tpiba2
      ENDDO
    ENDDO
  ENDDO

  ALLOCATE ( kmap ( nkstot ) )
  ALLOCATE ( smap ( nkstot ) )

  DO i = 1, nkstot
    ! ZL: change nspin to ns
    !j = ( i - 1 ) / nspin
    !k = i - 1 - j * nspin
    !kmap ( i ) = j + k * ( nkstot / nspin ) + 1
    j = ( i - 1 ) / ns
    k = i - 1 - j * ns
    kmap ( i ) = j + k * ( nkstot / ns ) + 1
    smap ( i ) = k + 1
  ENDDO
  ierr = 0
  DO i = 1, nkstot
    ik = kmap ( i )
    is = smap ( i )
    IF ( ik .GE. iks .AND. ik .LE. ike .AND. is .NE. isk ( ik ) ) &
      ierr = ierr + 1
  ENDDO
  CALL mp_max ( ierr, world_comm )
  IF ( ierr .GT. 0 ) &
    CALL errore ( 'write_vkbg', 'smap', ierr )

  ALLOCATE ( gvec ( 3, ngm_g ) )
  gvec = 0
  DO ig = 1, ngm
    gvec ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
    gvec ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
    gvec ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  ENDDO
  CALL mp_sum ( gvec, intra_bgrp_comm )

  ALLOCATE ( ngk_g ( nkstot ) )
  ALLOCATE ( igk_l2g ( npwx, nks ) )
  ngk_g = 0
  igk_l2g = 0
  DO ik = 1, nks
    npw = ngk ( ik )
    DO ig = 1, npw
      igk_l2g ( ig, ik ) = ig_l2g ( igk_k ( ig, ik ) )
    ENDDO
  ENDDO
  DO ik = 1, nks
    ngk_g ( ik + iks - 1 ) = ngk ( ik )
  ENDDO
  CALL mp_sum ( ngk_g, inter_pool_comm )
  CALL mp_sum ( ngk_g, intra_pool_comm )
  ngk_g = ngk_g / nbgrp
  npw_g = MAXVAL ( igk_l2g ( :, : ) )
  CALL mp_max ( npw_g, intra_pool_comm )
  !XX
  CALL mp_max( npw_g, inter_pool_comm )
  !XX
  npwx_g = MAXVAL ( ngk_g ( : ) )

  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), &
      form = 'unformatted', status = 'replace')
    WRITE ( unit ) stitle, sdate, stime
    ! ZL: change nspin to nsf, and nkstot/nspin to nkstot/ns
    !WRITE ( unit ) nspin, ngm_g, ntran, cell_symmetry, nat, ecutrho, &
    !  nkstot / nspin, nsp, nkb, nhm, npwx_g, ecutwfc
    WRITE ( unit ) nsf, ngm_g, ntran, cell_symmetry, nat, ecutrho, &
      nkstot / ns, nsp, nkb, nhm, npwx_g, ecutwfc
    IF ( wfng_kgrid ) THEN
      WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3, wfng_nk1, wfng_nk2, wfng_nk3, &
        wfng_dk1, wfng_dk2, wfng_dk3
    ELSE
      WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3, nk1, nk2, nk3, &
        0.5D0 * dble ( k1 ), 0.5D0 * dble ( k2 ), 0.5D0 * dble ( k3 )
    ENDIF
    WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
      ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
    WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
    WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
    ! ZL: change
    !WRITE ( unit ) ( ngk_g ( ik ), ik = 1, nkstot / nspin )
    !WRITE ( unit ) ( wk ( ik ) * dble ( nspin ) / 2.0D0, ik = 1, nkstot / nspin )
    !WRITE ( unit ) ( ( xk ( id, ik ), id = 1, nd ), ik = 1, nkstot / nspin )
    WRITE ( unit ) ( ngk_g ( ik ), ik = 1, nkstot / ns )
    WRITE ( unit ) ( wk ( ik ) * dble ( nst ) / 2.0D0, ik = 1, nkstot / ns )
    WRITE ( unit ) ( ( xk ( id, ik ), id = 1, nd ), ik = 1, nkstot / ns )
    WRITE ( unit ) ( ityp ( iat ), iat = 1, nat )
    WRITE ( unit ) ( nh ( isp ), isp = 1, nsp )
    IF ( nsf == 4 ) THEN
      WRITE ( unit ) ( ( ( ( deeq_nc ( jh, ih, iat, is ), &
        jh = 1, nhm ), ih = 1, nhm ), iat = 1, nat ), is = 1, 4 )
    ELSE
      WRITE ( unit ) ( ( ( ( deeq ( jh, ih, iat, is ), &
        jh = 1, nhm ), ih = 1, nhm ), iat = 1, nat ), is = 1, ns )
    ENDIF
      !jh = 1, nhm ), ih = 1, nhm ), iat = 1, nat ), is = 1, nspin ) ! ZL: change in above line
    WRITE ( unit ) nrecord
    WRITE ( unit ) ngm_g
    WRITE ( unit ) ( ( gvec ( id, ig ), id = 1, nd ), ig = 1, ngm_g )
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  ALLOCATE ( igwk ( npwx_g ) )

  DO i = 1, nkstot

    ik = kmap ( i )
    is = smap ( i )

    igwk = 0

    ALLOCATE ( itmp ( npw_g ) )
    itmp = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
      DO ig = 1, ngk ( ik - iks + 1 )
        itmp ( igk_l2g ( ig, ik - iks + 1 ) ) = igk_l2g ( ig, ik - iks + 1 )
      ENDDO
    ENDIF
    CALL mp_sum ( itmp, intra_bgrp_comm )
    !XX
    CALL mp_sum ( itmp, inter_pool_comm )
    !XX
    ngg = 0
    DO ig = 1, npw_g
      IF ( itmp ( ig ) .EQ. ig ) THEN
        ngg = ngg + 1
        igwk ( ngg ) = ig
      ENDIF
    ENDDO
    DEALLOCATE ( itmp )

    IF ( ionode ) THEN
      IF ( is .EQ. 1 ) THEN
        WRITE ( unit ) nrecord
        WRITE ( unit ) ngk_g ( ik )
        WRITE ( unit ) ( ( gvec ( j, igwk ( ig ) ), j = 1, 3 ), &
          ig = 1, ngk_g ( ik ) )
      ENDIF
    ENDIF

    local_pw = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
      npw = ngk ( ik - iks + 1 )
      CALL init_us_2 ( npw, igk_k(1, ik-iks+1), xk ( 1, ik ), vkb )
      local_pw = npw
    ENDIF

    ALLOCATE ( igwf_l2g ( local_pw ) )
    igwf_l2g = 0
    DO ig = 1, local_pw
      ngg = igk_l2g ( ig, ik - iks + 1 )
      DO j = 1, ngk_g ( ik )
        IF ( ngg .EQ. igwk ( j ) ) THEN
          igwf_l2g ( ig ) = j
          EXIT
        ENDIF
      ENDDO
    ENDDO

    ALLOCATE ( ipmask ( nproc ) )
    ipmask = 0
    ipsour = ionode_id
    !FZ: npool is not used in VKB, since it may have induce some error
    !IF ( npool .GT. 1 ) THEN
    !  IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
    !    IF ( me_pool .EQ. root_pool ) ipmask ( mpime + 1 ) = 1
    !  ENDIF
    !  CALL mp_sum ( ipmask, world_comm )
    !  DO j = 1, nproc
    !    IF ( ipmask ( j ) .EQ. 1 ) ipsour = j - 1
    !  ENDDO
    !ENDIF
    DEALLOCATE ( ipmask )

    igwx = 0
    ierr = 0
    IF ( ik .GE. iks .AND. ik .LE. ike ) &
      igwx = MAXVAL ( igwf_l2g ( 1 : local_pw ) )
    CALL mp_max ( igwx, intra_pool_comm )
    IF ( ipsour .NE. ionode_id ) &
      CALL mp_get ( igwx, igwx, mpime, ionode_id, ipsour, 1, world_comm )
    ierr = 0
    IF ( ik .GE. iks .AND. ik .LE. ike .AND. igwx .NE. ngk_g ( ik ) ) &
      ierr = 1
    CALL mp_max ( ierr, world_comm )
    IF ( ierr .GT. 0 ) &
      CALL errore ( 'write_vkbg', 'igwx ngk_g', ierr )

    ALLOCATE ( vkb_g ( MAX ( 1, igwx ) ) )

    DO ikb = 1, nkb

      vkb_g = ( 0.0D0, 0.0D0 )
      !FZ: npool is not used in VKB, since it may have induce some error
      !IF ( npool .GT. 1 ) THEN
      !  IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
      !    CALL mergewf ( vkb ( :, ikb ), vkb_g, local_pw, igwf_l2g, &
      !      me_pool, nproc_pool, root_pool, intra_pool_comm )
      !  ENDIF
      !  IF ( ipsour .NE. ionode_id ) THEN
      !    CALL mp_get ( vkb_g, vkb_g, mpime, ionode_id, ipsour, ikb, &
      !      world_comm )
      !  ENDIF
      !ELSE
      CALL mergewf ( vkb ( :, ikb ), vkb_g, local_pw, igwf_l2g, &
        mpime, nproc, ionode_id, world_comm )
      !ENDIF

      IF ( ionode ) THEN
        WRITE ( unit ) nrecord
        WRITE ( unit ) ngk_g ( ik )
        WRITE ( unit ) ( vkb_g ( ig ), ig = 1, igwx )
      ENDIF

    ENDDO

    DEALLOCATE ( vkb_g )
    DEALLOCATE ( igwf_l2g )

  ENDDO

  IF ( ionode ) THEN
    CLOSE ( unit = unit, status = 'keep' )
  ENDIF

  DEALLOCATE ( igwk )
  DEALLOCATE ( igk_l2g )
  DEALLOCATE ( ngk_g )
  DEALLOCATE ( gvec )
  DEALLOCATE ( smap )
  DEALLOCATE ( kmap )

  RETURN

END SUBROUTINE write_vkbg

SUBROUTINE write_vhub_g (output_file_name, diag_nmin, diag_nmax, offdiag_nmin, offdiag_nmax)

  USE constants, ONLY : rytoev
  USE cell_base, ONLY : tpiba2, at, bg
  USE ener, ONLY : etxc, vtxc
  USE exx, ONLY : vexx
  USE fft_base, ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  !USE funct, ONLY : exx_is_active
  USE xc_lib, ONLY : exx_is_active !FZ: for qe6.7
  USE gvect, ONLY : ngm, g
  USE io_files, ONLY : nwordwfc, iunwfc, nwordwfcU, iunhub
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP
  USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
  USE lsda_mod, ONLY : nspin, isk, current_spin, lsda
  USE ldaU, ONLY : lda_plus_u, Hubbard_projectors, wfcU, offsetU, Hubbard_l, nwfcU, Hubbard_lmax, is_hubbard
  USE mp, ONLY : mp_sum
  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE mp_pools, ONLY : intra_pool_comm
  USE scf, ONLY : rho, rho_core, rhog_core
  USE wavefunctions, ONLY : evc, psic
  USE wvfct, ONLY : npwx, nbnd
  USE noncollin_module , ONLY : noncolin , npol
  USE buffers, ONLY : get_buffer
  USE scf, ONLY : v
  USE ions_base,            ONLY : nat, ityp
  USE becmod,               ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type

  IMPLICIT NONE

  character (len = 256), intent (in) :: output_file_name
  integer, intent (inout) :: diag_nmin
  integer, intent (inout) :: diag_nmax
  integer, intent (inout) :: offdiag_nmin
  integer, intent (inout) :: offdiag_nmax
  integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2, ikk
  integer, external :: global_kpoint_index
  complex (DP) :: dummy
  complex (DP), allocatable :: mtxeld (:, :)
  complex (DP), allocatable :: mtxelo (:, :, :)
  complex (DP) :: hpsi(npwx*npol,nbnd)
  complex (DP), allocatable :: hc(:,:)
  integer :: nspin_
  integer :: kdim, kdmx
  integer :: ldim, is1, ibnd, i, na, m1, nt
  character(LEN=20) :: ik_string, ib_string, is_string

  if (diag_nmin > diag_nmax) then
    if (ionode) then
      call errore ( 'write_vxc_g', 'diag_nmin > diag_nmax', diag_nmin )
    endif
  endif
  IF (diag_nmin .LT. 1) diag_nmin = 1
  IF (diag_nmax .GT. nbnd) then
    if (ionode) then
      write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
    endif
    diag_nmax = nbnd
  ENDIF

  ndiag = MAX (diag_nmax - diag_nmin + 1, 0)

  if (offdiag_nmin > offdiag_nmax) then
    if (ionode) then
      call errore ( 'write_vxc_g', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )
    endif
  endif
  IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
  IF (offdiag_nmax .GT. nbnd)  then
    if (ionode) then
      write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
    endif
    offdiag_nmax = nbnd
  ENDIF
  noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

  IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

  ALLOCATE( hc( nbnd, nbnd) )

  unit = 4

  iks = global_kpoint_index (nkstot, 1)
  ike = iks + nks - 1

  IF (ndiag .GT. 0) THEN
    ALLOCATE (mtxeld (ndiag, nkstot))
    mtxeld (:, :) = (0.0D0, 0.0D0)
  ENDIF
  IF (noffdiag .GT. 0) THEN
    ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
    mtxelo (:, :, :) = (0.0D0, 0.0D0)
  ENDIF

  ldim = 2 * Hubbard_lmax + 1

  DO ik = iks, ike

    hpsi(:,:) = (0.0D0, 0.0D0)
    evc(:,:) = (0.0D0, 0.0D0)
    hc(:,:) = (0.0D0,0.0D0)
    ikk = ik - iks + 1
    npw = ngk ( ik - iks + 1 )

    IF ( npol == 1 ) THEN
      kdim = npw
      kdmx = npwx
    ELSE
      kdim = npwx*npol
      kdmx = npwx*npol
    ENDIF
    IF ( lsda ) current_spin = isk(ikk)
    CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
    IF ( lda_plus_u .AND. (Hubbard_projectors .NE. 'pseudo') ) THEN
      CALL get_buffer ( wfcU, nwordwfcU, iunhub, ikk )
    ENDIF

    IF ( lda_plus_u .AND. Hubbard_projectors.NE."pseudo" ) THEN
      IF (noncolin) THEN
        CALL vhpsi_nc( npwx, npw, nbnd, evc, hpsi )
      ELSE
        CALL vhpsi( npwx, npw, nbnd, evc, hpsi )
      ENDIF
    ENDIF

    ! MW: Use MATMUL instead of ZGEMM
    ! CALL ZGEMM( 'C', 'N', nbnd, nbnd, kdim, ( 1.D0, 0.D0 ), evc, kdmx,  hpsi, kdmx, ( 0.D0, 0.D0 ), hc, nbnd )
    hc = MATMUL(CONJG(TRANSPOSE(evc)), hpsi)

    CALL mp_sum(  hc , intra_bgrp_comm )

    ! diagonal
    IF (ndiag .GT. 0) THEN
      DO ib = diag_nmin, diag_nmax
        ! mtxeld in units of eV
        mtxeld (ib - diag_nmin + 1, ik) = hc (ib, ib) * CMPLX (rytoev, 0.0D0)
      ENDDO
    ENDIF
    ! off-diagonal
    IF (noffdiag .gt. 0) THEN
      DO ib = offdiag_nmin, offdiag_nmax
        DO ib2 = offdiag_nmin, offdiag_nmax
          ! mtxeld in units of eV
          mtxelo (ib2-offdiag_nmin+1, ib-offdiag_nmin+1, ik) = hc (ib2, ib) * CMPLX (rytoev, 0.0D0)
        ENDDO
      ENDDO
    ENDIF
  ENDDO ! ik_global

  DEALLOCATE (hc)

  IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
  IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)
  ! xk(1:3) : cartesian ==> crytal
  CALL cryst_to_cart (nkstot, xk, at, -1)

  IF (ionode) THEN
    OPEN (unit = unit, file = TRIM (output_file_name), form = 'formatted', status = 'replace')

    IF (nspin .eq. 4) THEN
      nspin_ = 1
    ELSE
      nspin_ = nspin
    ENDIF

    DO ik = 1, nkstot / nspin_
      WRITE (unit, 101) xk(:, ik), nspin_ * ndiag, nspin_ * noffdiag **2
      IF (ndiag .GT. 0) THEN
        DO ib = diag_nmin, diag_nmax
          DO is = 1, nspin_   
            WRITE (unit, 102) is, ib, mtxeld(ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin_)
          ENDDO
        ENDDO
      ENDIF
      IF (noffdiag .GT. 0) THEN
        DO ib = offdiag_nmin, offdiag_nmax
          DO ib2 = offdiag_nmin, offdiag_nmax
            DO is = 1, nspin_
              WRITE (unit, 103) is, ib2, ib, mtxelo (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, ik + (is - 1) * nkstot / nspin_)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    CLOSE (unit = unit, status = 'keep')
  ENDIF

  CALL cryst_to_cart (nkstot, xk, bg, 1)

  IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
  IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

  RETURN

101 FORMAT (3F13.9, 2I8)
102 FORMAT (2I8, 2F15.9)
103 FORMAT (3I8, 2F15.9)

END SUBROUTINE write_vhub_g

!-------------------------------------------------------------------------------

subroutine check_inversion(real_or_complex, ntran, mtrx, nspin, warn, real_need_inv, tnp)

! check_inversion    Originally By David A. Strubbe    Last Modified 11/18/2013
! Check whether our choice of real/complex version is appropriate given the
! presence or absence of inversion symmetry.

  USE constants, ONLY : eps6
  USE io_global, ONLY : ionode
  USE kinds, ONLY : DP

  implicit none

  integer, intent(in) :: real_or_complex
  integer, intent(in) :: ntran
  integer, intent(in) :: mtrx(3, 3, 48) !< symmetry operations matrices
  integer, intent(in) :: nspin
  logical, intent(in) :: warn !< set to false to suppress warnings, for converters
  logical, intent(in) :: real_need_inv !< use for generating routines to block real without inversion
     !! this is not always true so that it is possible to run real without using symmetries
  real(DP), intent(in) :: tnp(3, 48) !< fractional translations.
     !! optional only to avoid changing external interface for library.

  integer :: invflag, isym, ii, jj, itest
  logical :: origin_inv

  invflag = 0
  origin_inv = .false.
  do isym = 1, ntran
    itest = 0
    do ii = 1, 3
      do jj = 1, 3
        if(ii .eq. jj) then
          itest = itest + (mtrx(ii, jj, isym) + 1)**2
        else
          itest = itest + mtrx(ii, jj, isym)**2
        endif
      enddo
    enddo
    if(itest .eq. 0) then
      invflag = invflag + 1
      !if(present(tnp)) then
        if(sum(abs(tnp(1:3, isym))) < eps6) origin_inv = .true.
      !else
      !  origin_inv = .true.
      !endif
    endif
  enddo
  if(invflag > 0 .and. .not. origin_inv) then
    write(0, '(a)') "WARNING: Inversion symmetry is present only with a fractional translation."
    write(0, '(a)') "Apply the translation so inversion is about the origin, to be able to use the real version."
  endif
  if(invflag .gt. 1) write(0, '(a)') "WARNING: More than one inversion symmetry operation is present."

!  if(invflag > 0 .and. .not. present(tnp)) then
!    write(0, '(a)') "WARNING: check_inversion did not receive fractional translations."
!    write(0, '(a)') "Cannot confirm that inversion symmetry is about the origin for use of real version."
!  endif

  if(real_or_complex .eq. 2) then
    if(origin_inv .and. warn .and. nspin == 1) then
      if(ionode) &
        write(0, '(a)') "WARNING: Inversion symmetry about the origin is present. The real version would be faster."
    endif
  else
    if(.not. origin_inv) then
      if(real_need_inv) then
        call errore("check_inversion", "The real version cannot be used without inversion symmetry about the origin.", -1)
      endif
      if(ionode) then
        write(0, '(a)') "WARNING: Inversion symmetry about the origin is absent in symmetries used to reduce k-grid."
        write(0, '(a)') "Be sure inversion about the origin is still a spatial symmetry, or you must use complex version instead."
      endif
    endif
    if(nspin > 1) then
      call errore("check_inversion", "Real version may only be used for spin-unpolarized calculations.", nspin)
    endif
  endif

  return

end subroutine check_inversion

!-------------------------------------------------------------------------------

! Set spin-related parameters ns, nst, nsf
!
subroutine set_spin(ns, nst, nsf, nspin)

  IMPLICIT NONE

  integer, intent(out) :: ns, nst, nsf
  integer, intent(in) :: nspin

  IF ( nspin == 4 ) THEN
    ! ZL: here spinors check has already been performed
    ns = 1
    nst = 2
    nsf = 4
  ELSE
    ns = nspin
    nst = nspin
    nsf = nspin
  ENDIF

  RETURN

end subroutine set_spin

!-------------------------------------------------------------------------------

END PROGRAM pw2bgw

