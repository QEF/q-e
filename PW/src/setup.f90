!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE setup()
  !----------------------------------------------------------------------------
  !! This routine is called once at the beginning of the calculation and:
  !
  !! 1) determines various parameters of the calculation:
  !
  !!  * zv:        charge of each atomic type;
  !!  * nelec:     total number of electrons (if not given in input);
  !!  * nbnd:      total number of bands (if not given in input);
  !!  * nbndx:     max number of bands used in iterative diagonalization;
  !!  * tpiba:     2 pi / a (a = lattice parameter);
  !!  * tpiba2:    square of tpiba;
  !!  * gcutm:     cut-off in g space for charge/potentials;
  !!  * gcutms:    cut-off in g space for smooth charge;
  !!  * ethr:      convergence threshold for iterative diagonalization;
  !
  !! 2) finds actual crystal symmetry:
  !
  !!  * s:         symmetry matrices in the direct lattice vectors basis;
  !!  * nsym:      number of crystal symmetry operations;
  !!  * nrot:      number of lattice symmetry operations;
  !!  * ft:        fractionary translations;
  !!  * irt:       for each atom gives the corresponding symmetric;
  !!  * invsym:    if true the system has inversion symmetry;
  !
  !! 3) generates k-points corresponding to the actual crystal symmetry;
  !
  !! 4) calculates various quantities used in magnetic, spin-orbit, PAW
  !!    electric-field, DFT+U(+V) calculations
  !!
  USE kinds,              ONLY : DP
  USE constants,          ONLY : eps8, e2, fpi, pi, degspin
  USE parameters,         ONLY : npk
  USE io_global,          ONLY : stdout, ionode, ionode_id
  USE io_files,           ONLY : xmlfile
  USE cell_base,          ONLY : at, bg, alat, tpiba, tpiba2, ibrav
  USE ions_base,          ONLY : nat, tau, ntyp => nsp, ityp, zv
  USE starting_scf,       ONLY : starting_pot
  USE fft_support,        ONLY : good_fft_order
  USE gvect,              ONLY : gcutm, ecutrho
  USE gvecw,              ONLY : gcutw, ecutwfc
  USE gvecs,              ONLY : doublegrid, gcutms, dual
  USE klist,              ONLY : xk, wk, nks, nelec, degauss, lgauss, &
                                 ltetra, lxkcry, nkstot, &
                                 nelup, neldw, two_fermi_energies, &
                                 tot_charge, tot_magnetization
  USE ener,               ONLY : ef, ef_up, ef_dw
  USE electrons_base,     ONLY : set_nelup_neldw
  USE start_k,            ONLY : nks_start, xk_start, wk_start, &
                                 nk1, nk2, nk3, k1, k2, k3
  USE ktetra,             ONLY : tetra_type, opt_tetra_init, tetra_init
  USE symm_base,          ONLY : s, t_rev, irt, nrot, nsym, invsym, nosym, &
                                 d1,d2,d3, time_reversal, sname, set_sym_bl, &
                                 find_sym, inverse_s, no_t_rev, fft_fact,  &
                                 allfrac
  USE wvfct,              ONLY : nbnd, nbndx
  USE control_flags,      ONLY : tr2, ethr, lscf, lbfgs, lmd, david, lecrpa,  &
                                 isolve, niter, noinv, ts_vdw, tstress, &
                                 lbands, gamma_only, restart, use_spinflip, symm_by_label 
  USE cellmd,             ONLY : calc
  USE upf_ions,           ONLY : n_atom_wfc
  USE uspp_param,         ONLY : upf
  USE uspp,               ONLY : okvan
  USE ldaU,               ONLY : lda_plus_u, init_hubbard, lda_plus_u_kind, orbital_resolved
  USE bp,                 ONLY : gdir, lberry, nppstr, lelfield, lorbm, nx_el,&
                                 nppstr_3d,l3dstring, efield
  USE fixed_occ,          ONLY : f_inp, tfixed_occ, one_atom_occupations
  USE mp_pools,           ONLY : kunit, npool
  USE mp_images,          ONLY : intra_image_comm
  USE mp,                 ONLY : mp_bcast
  USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk, &
                                 starting_magnetization
  USE noncollin_module,   ONLY : noncolin, domag, npol, i_cons, m_loc, &
                                 angle1, angle2, bfield, ux, nspin_lsda, &
                                 nspin_gga, nspin_mag, lspinorb, colin_mag
  USE qexsd_module,       ONLY : qexsd_readschema
  USE qexsd_copy,         ONLY : qexsd_copy_efermi
  USE qes_libs_module,    ONLY : qes_reset
  USE qes_types_module,   ONLY : output_type
  USE exx,                ONLY : ecutfock
  USE xc_lib,             ONLY : xclib_dft_is
  USE paw_variables,      ONLY : okpaw
  USE extfield,           ONLY : gate
  USE esm,                ONLY : esm_z_inv
  USE fcp_module,         ONLY : lfcp
  USE gcscf_module,       ONLY : lgcscf
  USE rism_module,        ONLY : lrism, rism_calc1d
  USE additional_kpoints, ONLY : add_additional_kpoints
  USE control_flags,      ONLY : sic
  USE sic_mod,            ONLY : init_sic, occ_f2fn, sic_energy
  USE random_numbers,     ONLY : set_random_seed
  USE dynamics_module,    ONLY : control_temp
#if defined (__OSCDFT)
  USE plugin_flags,       ONLY : use_oscdft
  USE oscdft_base,        ONLY : oscdft_ctx
#endif

  !
  IMPLICIT NONE
  !
  INTEGER  :: na, is, ierr, ibnd, ik, nrot_, nbnd_, nr3, nk_, natomwfc 
  LOGICAL  :: magnetic_sym, skip_equivalence=.FALSE.
  REAL(DP) :: iocc, ionic_charge, one
  !
  TYPE(output_type)  :: output_obj 
  !  
  ! ... okvan/okpaw = .TRUE. : at least one pseudopotential is US/PAW
  !
  okvan = ANY( upf(1:ntyp)%tvanp )
  okpaw = ANY( upf(1:ntyp)%tpawp )
  !
  ! ... check for features not implemented with US-PP or PAW
  !
  IF ( okvan .OR. okpaw ) THEN
     IF ( xclib_dft_is('meta') ) CALL errore( 'setup', &
                               'Meta-GGA not implemented with USPP/PAW', 1 )
     IF ( noncolin .AND. lberry)  CALL errore( 'setup', &
       'Noncolinear Berry Phase/electric not implemented with USPP/PAW', 1 )
     IF  (ts_vdw ) CALL errore ('setup',&
                   'Tkatchenko-Scheffler not implemented with USPP/PAW', 1 )
     IF ( lorbm ) CALL errore( 'setup', &
                  'Orbital Magnetization not implemented with USPP/PAW', 1 )
  END IF

  IF ( xclib_dft_is('hybrid') ) THEN
     IF ( lberry ) CALL errore( 'setup ', &
                         'hybrid XC not allowed in Berry-phase calculations',1 )
     IF ( lelfield ) CALL errore( 'setup ', &
                         'hybrid XC and electric fields untested',1 )
     IF ( allfrac ) CALL errore( 'setup ', &
                         'option use_all_frac incompatible with hybrid XC', 1 )
!     IF (.NOT. lscf) CALL errore( 'setup ', &
!                         'hybrid XC not allowed in non-scf calculations', 1 )
     IF ( ANY (upf(1:ntyp)%nlcc) ) CALL infomsg( 'setup ', 'BEWARE:' // &
               & ' nonlinear core correction is not consistent with hybrid XC')
     IF (okvan) THEN
        IF (ecutfock /= 4.0_dp*ecutwfc) THEN
           ecutfock = MIN(4.0_dp*ecutwfc,ecutrho)
           CALL infomsg ('setup', &
                    'Warning: ecutfock not valid for US/PAW, ignored')
        END IF
        IF ( lmd .OR. lbfgs ) CALL errore &
           ('setup','forces for hybrid functionals + US/PAW not implemented',1)
        IF ( noncolin ) CALL errore &
           ('setup','Noncolinear hybrid XC for USPP not implemented',1)
     END IF
     IF ( noncolin ) no_t_rev=.true.
  END IF
  !
  IF ( xclib_dft_is('meta') ) THEN
     IF ( noncolin )  CALL errore( 'setup', &
                               'Non-collinear Meta-GGA not implemented', 1 )
     IF ( ANY (upf(1:ntyp)%nlcc) ) CALL infomsg( 'setup ', 'BEWARE:' // &
               & ' nonlinear core correction is not consistent with meta-GGA')
  END IF
  !
  ! ... Compute the ionic charge for each atom type and the total ionic charge
  !
  zv(1:ntyp) = upf(1:ntyp)%zp
  !
#if defined (__PGI)
     ionic_charge = 0._DP
     DO na = 1, nat
        ionic_charge = ionic_charge + zv( ityp(na) )
     END DO
#else
     ionic_charge = SUM( zv(ityp(1:nat)) ) 
#endif
  !
  ! ... set the number of electrons 
  !
  nelec = ionic_charge - tot_charge
  !
  IF ( .NOT. lscf .OR. ( (lfcp .OR. lgcscf) .AND. restart ) ) THEN
     !
     ! ... in these cases, we need (or it is useful) to read the Fermi energy
     ! ... also, number of bands is needed for FCP/GC-SCF
     !
     IF (ionode) CALL qexsd_readschema ( xmlfile(), ierr, output_obj )
     CALL mp_bcast(ierr, ionode_id, intra_image_comm)
     IF ( ierr > 0 ) CALL errore ( 'setup', 'problem reading ef from file ' //&
          & TRIM(xmlfile()), ierr )
     IF (ionode) CALL qexsd_copy_efermi ( output_obj%band_structure, &
          nelec, ef, two_fermi_energies, ef_up, ef_dw, nbnd_ )
     ! convert to Ry a.u. 
     ef = ef*e2
     ef_up = ef_up*e2
     ef_dw = ef_dw*e2
     CALL mp_bcast(nelec, ionode_id, intra_image_comm)
     CALL mp_bcast(ef, ionode_id, intra_image_comm)
     CALL mp_bcast(two_fermi_energies, ionode_id, intra_image_comm)
     CALL mp_bcast(ef_up, ionode_id, intra_image_comm)
     CALL mp_bcast(ef_dw, ionode_id, intra_image_comm)
     CALL mp_bcast(nbnd_, ionode_id, intra_image_comm)
     CALL qes_reset  ( output_obj )
     !
  END IF
  !
  IF ( (lfcp .OR. lgcscf) .AND. restart ) THEN
     tot_charge = ionic_charge - nelec
     nbnd = nbnd_
  END IF
  !
  ! ... magnetism-related quantities
  !
  ! ... Set the domag variable to make a spin-orbit calculation with zero
  ! ... magnetization
  !
  IF ( noncolin  ) THEN
     domag = ANY ( ABS( starting_magnetization(1:ntyp) ) > 1.D-6 )
  ELSE
     domag = .false.
  END IF
  !
  !  Set the different spin indices
  !
  CALL set_spin_vars( lsda, noncolin, domag, &
         npol, nspin, nspin_lsda, nspin_mag, nspin_gga, current_spin )
  ! set colin_mag.
  IF (symm_by_label .AND. nspin == 2 .AND. (ANY ( ABS( starting_magnetization(1:ntyp) ) > 1.D-6)) ) THEN 
    IF (use_spinflip) THEN 
       colin_mag = 2
    ELSE 
       colin_mag = 1 
    END IF 
  ELSE IF (symm_by_label) THEN 
     colin_mag = 0
  END IF
  IF ( xclib_dft_is('hybrid') ) THEN
     IF ( colin_mag == 2 ) THEN
        CALL infomsg( 'setup', 'colin_mag=2 not implemented for hybrid' )
        colin_mag = 1
     ENDIF
  ENDIF
  IF (lda_plus_u .AND. lda_plus_u_kind == 2) THEN
     IF ( colin_mag == 2 ) THEN
        CALL infomsg( 'setup', 'colin_mag=2 not implemented for lda+U+V' )
        colin_mag = 1
     ENDIF
  ENDIF
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) THEN
     IF ( colin_mag == 2 ) THEN
        CALL infomsg( 'setup', 'colin_mag=2 not implemented for OSCDFT' )
        colin_mag = 1
     ENDIF
  ENDIF
#endif
  !
  ! time reversal operation is set up to 0 by default
  t_rev = 0
  !
  ALLOCATE( m_loc( 3, nat ) )
  IF ( noncolin ) THEN
     !
     ! gamma_only and noncollinear not allowed
     !
     if (gamma_only) call errore('setup', &
                                 'gamma_only and noncolin not allowed',1)
     !
     DO na = 1, nat
        !
        m_loc(1,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * COS( angle2(ityp(na)) )
        m_loc(2,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * SIN( angle2(ityp(na)) )
        m_loc(3,na) = starting_magnetization(ityp(na)) * &
                      COS( angle1(ityp(na)) )
     END DO
     !
     !  initialize the quantization direction for gga
     !
     ux=0.0_DP
     if (xclib_dft_is('gradient')) call compute_ux(m_loc,ux,nat)
     !
  ELSE
     !
     IF (lspinorb)  CALL errore( 'setup ',  &
         'spin orbit requires a non collinear calculation', 1 )
     !
     IF ( i_cons == 1) then
        do na=1,nat
           m_loc(1,na) = starting_magnetization(ityp(na))
        end do
     !  set initial magnetization for collinear case
     ELSE IF ( colin_mag >= 1 ) THEN
        DO na = 1, nat
            m_loc(1,na) = 0.0_dp
            m_loc(2,na) = 0.0_dp
            m_loc(3,na) = starting_magnetization(ityp(na))
        END DO
     ENDIF     

     IF ( i_cons /= 0 .AND. nspin==1 ) &
        CALL errore( 'setup', 'this i_cons requires a magnetic calculation ', 1 )
     IF ( i_cons /= 0 .AND. i_cons /= 1 ) &
          CALL errore( 'setup', 'this i_cons requires a non colinear run', 1 )
     !
  END IF
  !
  ! ... if this is not a spin-orbit calculation, all spin-orbit pseudopotentials
  ! ... are transformed into standard pseudopotentials
  !
  IF ( lspinorb ) THEN
     IF ( ALL ( .NOT. upf(:)%has_so ) ) CALL errore ('setup', &
         'Spin-orbit calculations require at least one spin-orbit pseudo',1)
     IF ( ANY ( .NOT. upf(:)%has_so ) ) CALL infomsg ('setup', &
         'Not all pseudopotentials have spin-orbit data')
  ELSE
     CALL average_pp ( ntyp )
  END IF
  !
  ! ... If the occupations are from input, check the consistency with the
  ! ... number of electrons
  !
  IF ( tfixed_occ ) THEN
     !
     iocc = 0
     !
     DO is = 1, nspin_lsda
        !
#if defined (__PGI)
        DO ibnd = 1, nbnd
           iocc = iocc + f_inp(ibnd,is)
        END DO
#else
        iocc = iocc + SUM( f_inp(1:nbnd,is) )
#endif
        !
        DO ibnd = 1, nbnd
           if (f_inp(ibnd,is) > 2.d0/nspin_lsda .or. f_inp(ibnd,is) < 0.d0) &
              call errore('setup','wrong fixed occupations',is)
        END DO
     END DO
     !
     IF ( ABS( iocc - nelec ) > 1D-5 ) &
        CALL errore( 'setup', 'strange occupations: '//&
                     'number of electrons from occupations is wrong.', 1 )
     !
  END IF
  !
  ! ... Check: if there is an odd number of electrons, the crystal is a metal
  !
  IF ( lscf .AND. ABS( NINT( nelec / 2.D0 ) - nelec / 2.D0 ) > eps8 &
            .AND. .NOT. lgauss .AND. .NOT. ltetra .AND. .NOT. tfixed_occ ) &
      CALL infomsg( 'setup', 'the system is metallic, specify occupations' )
  !
  ! ... Check: spin-polarized calculations require either broadening or
  !             fixed occupation
  !
  IF ( lscf .AND. lsda &
            .AND. .NOT. lgauss .AND. .NOT. ltetra &
            .AND. .NOT. tfixed_occ .AND. .NOT. two_fermi_energies ) &
      CALL errore( 'setup', 'spin-polarized system, specify occupations', 1 )
  !
  ! ... setting the number of up and down electrons
  !
  call set_nelup_neldw ( tot_magnetization, nelec, nelup, neldw )
  !
  ! ... Set the number of occupied bands if not given in input
  !
  IF ( nbnd == 0 ) THEN
     !
     IF (nat==0) CALL errore('setup','free electrons: nbnd required in input',1)
     !
     nbnd = MAX ( NINT( nelec / degspin ), NINT(nelup), NINT(neldw) )
     !
     IF ( lgauss .OR. ltetra ) THEN
        !
        ! ... metallic case: add 20% more bands, with a minimum of 4
        !
        nbnd = MAX( NINT( 1.2D0 * nelec / degspin ), &
                    NINT( 1.2D0 * nelup), NINT( 1.2d0 * neldw ), &
                    ( nbnd + 4 ) )
        !
     END IF
     !
     ! ... In the case of noncollinear magnetism, bands are NOT
     ! ... twofold degenerate :
     !
     IF ( noncolin ) nbnd = INT( degspin ) * nbnd
     !
  ELSE
     !
     IF ( nbnd < NINT( nelec / degspin ) .AND. lscf ) &
        CALL errore( 'setup', 'too few bands', 1 )
     !
     IF ( nbnd < NINT( nelup ) .AND. lscf ) &
        CALL errore( 'setup', 'too few spin up bands', 1 )
     IF ( nbnd < NINT( neldw ) .AND. lscf ) &
        CALL errore( 'setup', 'too few spin dw bands', 1 )
     !
     IF ( nbnd < NINT( nelec ) .AND. lscf .AND. noncolin ) &
        CALL errore( 'setup', 'too few bands noncolin case', 1 )
     !
  END IF
  !
  ! ... Here we  set the precision of the diagonalization for the first scf
  ! ... iteration of for the first ionic step
  ! ... for subsequent steps ethr is automatically updated in electrons
  !
  IF ( nat==0 ) THEN
     !
     ethr=1.0D-8
     !
  ELSE IF ( .NOT. lscf ) THEN
     !
     ! ... do not allow convergence threshold of scf and nscf to become too small 
     ! 
     IF ( ethr == 0.D0 ) ethr = MAX(1.D-13, 0.1D0 * MIN( 1.D-2, tr2 / nelec ))
     !
  ELSE
     !
     IF ( ethr == 0.D0 ) THEN
        !
        IF ( starting_pot == 'file' ) THEN
           !
           ! ... if you think that the starting potential is good
           ! ... do not spoil it with a lousy first diagonalization :
           ! ... set a strict ethr in the input file (diago_thr_init)
           !
           IF ( lgcscf .OR. orbital_resolved ) THEN
              !
              ethr = 1.D-8
              !
           ELSE
              !
              ethr = 1.D-5
              !
           END IF
           !
        ELSE
           !
           ! ... starting atomic potential is probably far from scf
           ! ... do not waste iterations in the first diagonalizations
           !
           IF ( lgcscf ) THEN
              !
              ethr = 1.0D-5
              !
           ELSE
              !
              ethr = 1.0D-2
              !
           END IF
           !
        END IF
        !
     END IF
     !
  END IF
  !
  IF ( .NOT. lscf ) niter = 1
  !
  !
  ! ... set the max number of bands used in iterative diagonalization
  !
  nbndx = nbnd
  IF ( isolve == 0  ) nbndx = david * nbnd 
  IF (isolve == 4 ) nbndx = 2 *nbnd 
  !
  ! ... Set the units in real and reciprocal space
  !
  tpiba  = 2.D0 * pi / alat
  tpiba2 = tpiba**2
  !
  ! ... Compute the cut-off of the G vectors
  !
  doublegrid = ( dual > 4.0_dp + eps8 )
  IF ( doublegrid .AND. ( .NOT.okvan .AND. .NOT.okpaw .AND. &
                          .NOT. ANY (upf(1:ntyp)%nlcc) .AND. &
                          .NOT. lrism ) ) &
       CALL infomsg ( 'setup', 'no reason to have ecutrho>4*ecutwfc' )
  IF ( ecutwfc > 10000.d0 .OR. ecutwfc < 1.d0 ) THEN
       WRITE(stdout,*) 'ECUTWFC = ', ecutwfc
       CALL errore ( 'setup', 'meaningless value for ecutwfc', 1)
  END IF
  gcutm = dual * ecutwfc / tpiba2
  gcutw = ecutwfc / tpiba2
  !
  IF ( doublegrid ) THEN
     !
     gcutms = 4.D0 * ecutwfc / tpiba2
     !
  ELSE
     !
     gcutms = gcutm
     !
  END IF
  !
  ! ... Test that atoms do not overlap
  !
  call check_atoms ( nat, tau, bg )
  !
  !  ... generate transformation matrices for the crystal point group
  !  ... First we generate all the symmetry matrices of the Bravais lattice
  !
  call set_sym_bl ( )
  !
  ! ... If lecrpa is true, nosym must be set to true also
  !
  IF ( lecrpa ) nosym = .TRUE.
  IF ( lecrpa ) skip_equivalence=.TRUE.
  !
  ! ... if nosym is true: do not reduce automatic k-point grids to the IBZ
  ! ... using the symmetries of the lattice (only k <-> -k symmetry is used)
  ! ... Does not change the number "nrot" of symmetries of the lattice
  !
  IF ( nosym ) THEN
     nrot_ = 1
  ELSE
     nrot_ = nrot
  END IF
  !
  ! ... time_reversal = use q=>-q symmetry for k-point generation
  !
  magnetic_sym = noncolin .AND. domag 
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  !
  ! ... Automatic generation of k-points (if required)
  !
  IF ( nks_start == 0 ) THEN
     !
     IF (lelfield .OR. lorbm) THEN
         !
        CALL kpoint_grid_efield (at,bg, npk, &
             k1,k2,k3, nk1,nk2,nk3, nkstot, xk, wk, nspin)
        nosym = .TRUE.
        nrot_ = 1
        !
     ELSE IF (lberry ) THEN
        !
        CALL kp_strings( nppstr, gdir, nrot_, s, bg, npk, &
                         k1, k2, k3, nk1, nk2, nk3, nkstot, xk, wk )
        nosym = .TRUE.
        nrot_ = 1
        !
     ELSE
        !
        CALL kpoint_grid ( nrot_,time_reversal, skip_equivalence, s, t_rev, bg,&
                           npk, k1,k2,k3, nk1,nk2,nk3, nkstot, xk, wk )
        !
     END IF
     !
  ELSE 
     nkstot = nks_start
     xk(:,1:nkstot) = xk_start(:,1:nks_start)
     wk(1:nkstot) = wk_start(1:nks_start)
     !
     IF( lelfield) THEN
        !
        IF(noncolin) THEN
           allocate(nx_el(nkstot,3))
        ELSE
           allocate(nx_el(nkstot*nspin,3))
        END IF

        IF ( gdir<1 .OR. gdir>3 ) CALL errore('setup','invalid gdir value'&
                                  &' (valid values: 1=x, 2=y, 3=z)',10) 
        DO ik=1,nkstot
           nx_el(ik,gdir)=ik
        END DO
        ! sanity check (when nkstot==1 we /could/ just set nppstr=1):
        IF(nppstr==0) CALL errore('setup', 'When lefield is true and kpoint are '&
                                 &'specified manually you MUST set nppstr',1)
        if(nspin==2) nx_el(nkstot+1:2*nkstot,:) = nx_el(1:nkstot,:) + nkstot
        nppstr_3d(gdir)=nppstr
        l3dstring=.false.
        nosym = .TRUE.
        nrot_ = 1
        !
     END IF
  END IF
  !
  CALL add_additional_kpoints(nkstot, xk, wk)
  !
  IF ( nat==0 ) THEN
     !
     nsym=nrot
     invsym=.true.
     CALL inverse_s ( ) 
     !
  ELSE
     !
     ! ... eliminate rotations that are not symmetry operations
     !
     CALL find_sym ( nat, tau, ityp, magnetic_sym, m_loc, &
                   & gate .OR. (.NOT. esm_z_inv(lrism)) )
     !
     ! ... do not force FFT grid to be commensurate with fractional translations
     !
     IF ( allfrac ) fft_fact(:) = 1 
     !
  END IF
  !
  ! ... nosym: do not use any point-group symmetry (s(:,:,1) is the identity)
  !
  IF ( nosym ) THEN
     nsym = 1
     invsym = .FALSE.
     fft_fact(:) = 1
  END IF
  !
  IF ( nsym > 1 .AND. ibrav == 0 ) CALL infomsg('setup', &
       'using ibrav=0 with symmetry is DISCOURAGED, use correct ibrav instead')
  !
  ! ... Input k-points are assumed to be  given in the IBZ of the Bravais
  ! ... lattice, with the full point symmetry of the lattice.
  ! ... If some symmetries of the lattice are missing in the crystal,
  ! ... "irreducible_BZ" computes the missing k-points.
  !
  IF ( .NOT. lbands ) THEN
     CALL irreducible_BZ (nrot_, s, nsym, time_reversal, &
                          magnetic_sym, at, bg, npk, nkstot, xk, wk, t_rev )
  ELSE
     one = SUM (wk(1:nkstot))
     IF ( one > 0.0_dp ) wk(1:nkstot) = wk(1:nkstot) / one
  END IF
  !
  ! ... if dynamics is done the system should have no symmetries
  ! ... (inversion symmetry alone is allowed)
  !
  IF ( lmd .AND. ( nsym == 2 .AND. .NOT. invsym .OR. nsym > 2 ) &
           .AND. .NOT. ( calc == 'mm' .OR. calc == 'nm' ) ) &
       CALL infomsg( 'setup', 'Dynamics, you should have no symmetries' )
  !
  IF ( ltetra ) THEN
     !
     ! ... Calculate quantities used in tetrahedra method
     !
     IF (nks_start /= 0) CALL errore( 'setup ', 'tetrahedra need automatic k-point grid',1)
     !
     IF (tetra_type == 0) then
        CALL tetra_init( nsym, s, time_reversal, t_rev, at, bg, npk, k1,k2,k3, &
             nk1, nk2, nk3, nkstot, xk )
     ELSE 
        CALL opt_tetra_init(nsym, s, time_reversal, t_rev, at, bg, npk, &
             k1, k2, k3, nk1, nk2, nk3, nkstot, xk, 1)
     END IF
     !
  END IF
  !
  IF ( lsda ) THEN
     !
     ! ... LSDA case: two different spin polarizations,
     ! ...            each with its own kpoints
     !
     if (nspin /= 2) call errore ('setup','nspin should be 2; check iosys',1)
     !
     CALL set_kup_and_kdw( xk, wk, isk, nkstot, npk )
     !
  ELSE IF ( noncolin ) THEN
     !
     ! ... noncolinear magnetism: potential and charge have dimension 4 (1+3)
     !
     if (nspin /= 4) call errore ('setup','nspin should be 4; check iosys',1)
     current_spin = 1
     isk(:) = 1
     !
  ELSE
     !
     ! ... LDA case: the two spin polarizations are identical
     !
     wk(1:nkstot)    = wk(1:nkstot) * degspin
     current_spin = 1
     isk(:) = 1
     !
     IF ( nspin /= 1 ) &
        CALL errore( 'setup', 'nspin should be 1; check iosys', 1 )
     !
  END IF
  !
  IF ( nkstot > npk ) CALL errore( 'setup', 'too many k points', nkstot )
  !
  IF (one_atom_occupations) THEN
     natomwfc = n_atom_wfc( nat, ityp, noncolin )
     DO ik=1,nkstot
        DO ibnd=natomwfc+1, nbnd
           IF (f_inp(ibnd,ik)> 0.0_DP) CALL errore('setup', &
               'no atomic wavefunction for some band',1)
        ENDDO
     ENDDO
  ENDIF
  !
  ! ... Set up Hubbard parameters for DFT+Hubbard
  !
  CALL init_hubbard ( upf(1:ntyp)%psd, nspin, noncolin )
  !
  ! ... initialize d1 and d2 to rotate the spherical harmonics
  !
  IF (lda_plus_u .or. okpaw .or. (okvan.and.xclib_dft_is('hybrid')) ) CALL d_matrix( d1, d2, d3 )
  !
  ! ... set parallelization strategy: need an estimate of FFT dimension along z
  !
  nr3 = int ( sqrt(gcutms)*sqrt (at(1, 3)**2 + at(2, 3)**2 + at(3, 3)**2) ) + 1
  nr3 = good_fft_order( 2*nr3, fft_fact(3) )
  !
  ! ... nk_ = number of k-points usable for k-point parallelization
  !           (set to 1 if k-point parallelization is not implemented)
  nk_ = nkstot
  IF ( lberry .OR. lelfield .OR. lorbm .OR. &
       ( xclib_dft_is('hybrid') .AND. tstress ) ) nk_ = 1
  !
  CALL setup_para ( nr3, nk_, nbnd )
  !
  ! ... distribute k-points across processors of a pool
  !
  kunit   = 1
  CALL divide_et_impera ( nkstot, xk, wk, isk, nks )
  !
  ! ... checks and initializations to be performed after parallelization setup
  !
  IF ( lberry .OR. lelfield .OR. lorbm ) THEN
     IF ( npool > 1 ) CALL errore( 'setup', &
          'Berry Phase/electric fields not implemented with pools', 1 )
  END IF
  IF ( gamma_only .AND. nkstot == 1 .AND. npool > 1 ) CALL errore( 'setup', &
          'Gamma-only calculations not allowed with pools', 1 )
  IF ( xclib_dft_is('hybrid') ) THEN
     IF ( nks == 0 ) CALL errore('setup','pools with no k-points' &
          & // ' not allowed for hybrid functionals',1)
     IF ( tstress .and. npool > 1 )  CALL errore('setup', &
         'stress for hybrid functionals not available with pools', 1)
     !!!IF ( tstress )  CALL errore('setup', &
     !!!    'stress for hybrid functionals not available', 1)
     !
     CALL setup_exx  ()
     !
  END IF
  !
  ! ... calculate solvent-solvent interaction (1D-RISM).
  !
  IF (lrism) CALL rism_calc1d()
  !
  ! ... SIC calculation
  !
  IF(sic) THEN
     CALL init_sic()
     IF (sic_energy) CALL occ_f2fn()
  END IF
  !
  ! ... next command prevents different MD runs to start
  ! ... with exactly the same "random" velocities
  !
  IF (lmd.AND.control_temp) CALL set_random_seed( )
  !
END SUBROUTINE setup
!
!----------------------------------------------------------------------------
SUBROUTINE setup_para ( nr3, nkstot, nbnd )
  !----------------------------------------------------------------------------
  !
  ! Initialize the various parallelization levels, trying to guess decent
  ! parameters for npool, ndiag, ntg, if not specified in the command line
  ! Must be called at the end of "setup" but before "setup_exx",
  ! only once per run
  !
  USE control_flags, ONLY : use_para_diag, use_gpu
  USE io_global, ONLY : stdout
  USE mp_bands,  ONLY : mp_start_bands, nbgrp, ntask_groups, nproc_bgrp,&
                        nyfft
  USE mp_pools,  ONLY : intra_pool_comm, mp_start_pools, npool
  USE mp_images, ONLY : intra_image_comm, nproc_image
  USE command_line_options, ONLY : npool_, ndiag_, nband_, ntg_, nyfft_, &
          nmany_, pencil_decomposition_
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nr3
  INTEGER, INTENT(in) :: nkstot
  INTEGER, INTENT(in) :: nbnd
  !
  LOGICAL, EXTERNAL  :: check_gpu_support
  LOGICAL, SAVE :: first = .TRUE.
  INTEGER :: maxtask, np
  !
  ! do not execute twice: unpredictable results may follow
  !
  IF ( .NOT.first ) RETURN
  first = .false.
  !
  ! GPUs (not sure it serves any purpose)
  !
  use_gpu = check_gpu_support( )
  !
  ! k-point parallelization first
  !
  IF ( npool_== 0 ) THEN
     npool_ = 1
     !
     if ( nproc_image > nr3/2 .and. nkstot > 1 ) then
        !
        ! if too many mpi processes for this fft dimension,
        ! use k-point parallelization if available
        !
pool:   do np = 2, nkstot
           ! npool should be a divisor of the number of processors
           if ( mod(nproc_image, np) /= 0 ) cycle
           npool_= np
           if ( nproc_image/np <= nr3/2 ) exit pool
        end do pool
     end if
  END IF
  CALL mp_start_pools ( npool_, intra_image_comm )
  !
  ! band parallelization (nband only; ntg and nyfft below are just copied)
  !
#if defined (__CUDA_OPTIMIZED)
  ! band parallelization is effectively disabled in this case
  CALL mp_start_bands ( 1 , ntg_, nyfft_, intra_pool_comm )
#else
  CALL mp_start_bands ( nband_, ntg_, nyfft_, intra_pool_comm )
#endif
  !
  ! Set "task groups", max 16, if too many processors for PW parallelization
  !
  IF ( ntask_groups == 0 ) THEN
     ntask_groups = 1
     if ( nproc_bgrp > nr3 ) THEN
        maxtask = min (nbnd, 16)
task:   do np = 2, maxtask
           ! ensure that ntask_group, that coincides with nyfft,
           ! is commensurate with the number of procs for PW parallelization
           if ( mod(nproc_bgrp,np) /= 0 ) cycle
           if ( nproc_bgrp/np < nr3/4 ) then 
              ntask_groups = np
              exit task 
           end if
        end do task
      end if
  END IF
  !
  ! Note that "task_groups" require to set pencil_decomposition to .true.
  ! Same if there are more processors than xy planes in the FFT
  !
  IF ( ntask_groups /= 1 .OR. nproc_bgrp > nr3 ) pencil_decomposition_ = .true. 
  !
  ! printout - same as in environment.f90
  !
  WRITE( stdout, * )
  IF ( npool > 1 ) WRITE( stdout, &
         '(5X,"K-points division:     npool     = ",I7)' ) npool
  IF ( nbgrp > 1 ) WRITE( stdout, &
         '(5X,"band groups division:  nbgrp     = ",I7)' ) nbgrp
  IF ( nproc_bgrp > 1 ) WRITE( stdout, &
         '(5X,"R & G space division:  proc/nbgrp/npool/nimage = ",I7)' ) nproc_bgrp
  IF ( nproc_bgrp > nr3 ) WRITE( stdout, &
         '(5X,"WARNING: too many processors for an effective parallelization!")' ) 
  IF ( nyfft > 1 ) WRITE( stdout, &
         '(5X,"wavefunctions fft division:  Y-proc x Z-proc = ",2I7)' ) &
         nyfft, nproc_bgrp / nyfft
  IF ( ntask_groups > 1 ) WRITE( stdout, &
         '(5X,"wavefunctions fft division:  task group distribution",/,34X,"#TG    x Z-proc = ",2I7)' ) &
         ntask_groups, nproc_bgrp / ntask_groups
  IF ( nmany_ > 1) WRITE( stdout, '(5X,"FFT bands division:     nmany     = ",I7)' ) nmany_
  !
  ! linear algebra - for GPUs, ndiag = 1; otherwise, a value ensuring that
  ! matrices nbnd*nbnd are distributed into blocks of size > 100x100
  !
  if ( ndiag_ == 0 .AND. use_gpu ) ndiag_ = 1
  if ( ndiag_ == 0 ) then
     do np = nint(nbnd/100.), 1, -1
         if ( np**2 <= nproc_bgrp ) exit
     end do
     ndiag_ = max(1,np**2)
  end if
  !
  CALL set_para_diag( nbnd, use_para_diag )
  !
END SUBROUTINE setup_para
!
!----------------------------------------------------------------------------
LOGICAL FUNCTION check_gpu_support( )
  ! Minimal case: returns true if compiled for GPUs
  IMPLICIT NONE
  !
#if defined(__CUDA)
  check_gpu_support = .TRUE.
#else
  check_gpu_support = .FALSE.
#endif
  !
END FUNCTION check_gpu_support
!
!----------------------------------------------------------------------------
SUBROUTINE setup_exx ( )
  !----------------------------------------------------------------------------
  !
  ! Must be called after setup_para, before init_run, only once
  !
  USE mp_exx,    ONLY : mp_start_exx
  USE mp_pools,  ONLY : intra_pool_comm
  USE exx_base,  ONLY : exx_grid_init, exx_mp_init, exx_div_check
  USE command_line_options, ONLY : nband_, ntg_
  !
  IMPLICIT NONE
  !
  CALL mp_start_exx ( nband_, ntg_, intra_pool_comm )
  CALL exx_grid_init()
  CALL exx_mp_init()
  CALL exx_div_check()
  !
END SUBROUTINE setup_exx
