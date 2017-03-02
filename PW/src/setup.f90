!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE setup()
  !----------------------------------------------------------------------------
  !
  ! ... This routine is called at the beginning of the calculation and
  ! ... 1) determines various parameters of the calculation:
  ! ...    zv        charge of each atomic type
  ! ...    nelec     total number of electrons (if not given in input)
  ! ...    nbnd      total number of bands (if not given in input)
  ! ...    nbndx     max number of bands used in iterative diagonalization
  ! ...    tpiba     2 pi / a (a = lattice parameter)
  ! ...    tpiba2    square of tpiba
  ! ...    gcutm     cut-off in g space for charge/potentials
  ! ...    gcutms    cut-off in g space for smooth charge
  ! ...    ethr      convergence threshold for iterative diagonalization
  ! ... 2) finds actual crystal symmetry:
  ! ...    s         symmetry matrices in the direct lattice vectors basis
  ! ...    nsym      number of crystal symmetry operations
  ! ...    nrot      number of lattice symmetry operations
  ! ...    ft        fractionary translations
  ! ...    irt       for each atom gives the corresponding symmetric
  ! ...    invsym    if true the system has inversion symmetry
  ! ... 3) generates k-points corresponding to the actual crystal symmetry
  ! ... 4) calculates various quantities used in magnetic, spin-orbit, PAW
  ! ...    electric-field, LDA+U calculations, and for parallelism
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : eps8, rytoev, fpi
  USE parameters,         ONLY : npk
  USE io_global,          ONLY : stdout
  USE io_files,           ONLY : tmp_dir, prefix, xmlpun, delete_if_present
  USE constants,          ONLY : pi, degspin
  USE cell_base,          ONLY : at, bg, alat, tpiba, tpiba2, ibrav, omega
  USE ions_base,          ONLY : nat, tau, ntyp => nsp, ityp, zv
  USE basis,              ONLY : starting_pot, natomwfc
  USE gvect,              ONLY : gcutm, ecutrho
  USE gvecw,              ONLY : gcutw, ecutwfc
  USE fft_base,           ONLY : dfftp
  USE fft_base,           ONLY : dffts
  USE fft_types,          ONLY : fft_type_init, fft_type_allocate
  USE fft_base,           ONLY : smap
  USE gvecs,              ONLY : doublegrid, gcutms, dual
  USE klist,              ONLY : xk, wk, nks, nelec, degauss, lgauss, &
                                 ltetra, lxkcry, nkstot, &
                                 nelup, neldw, two_fermi_energies, &
                                 tot_charge, tot_magnetization
  USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk, &
                                 starting_magnetization
  USE ener,               ONLY : ef
  USE electrons_base,     ONLY : set_nelup_neldw
  USE start_k,            ONLY : nks_start, xk_start, wk_start, &
                                 nk1, nk2, nk3, k1, k2, k3
  USE ktetra,             ONLY : tetra_type, opt_tetra_init, tetra_init
  USE symm_base,          ONLY : s, t_rev, irt, nrot, nsym, invsym, nosym, &
                                 d1,d2,d3, time_reversal, sname, set_sym_bl, &
                                 find_sym, inverse_s, no_t_rev &
                                 , allfrac, remove_sym
  USE wvfct,              ONLY : nbnd, nbndx
  USE control_flags,      ONLY : tr2, ethr, lscf, lmd, david, lecrpa,  &
                                 isolve, niter, noinv, ts_vdw, &
                                 lbands, use_para_diag, gamma_only, &
                                 restart
  USE cellmd,             ONLY : calc
  USE uspp_param,         ONLY : upf, n_atom_wfc
  USE uspp,               ONLY : okvan
  USE ldaU,               ONLY : lda_plus_u, init_lda_plus_u
  USE bp,                 ONLY : gdir, lberry, nppstr, lelfield, lorbm, nx_el,&
                                 nppstr_3d,l3dstring, efield
  USE fixed_occ,          ONLY : f_inp, tfixed_occ, one_atom_occupations
  USE funct,              ONLY : set_dft_from_name
  USE mp_pools,           ONLY : kunit
  USE mp_bands,           ONLY : intra_bgrp_comm
  USE spin_orb,           ONLY : lspinorb, domag
  USE noncollin_module,   ONLY : noncolin, npol, m_loc, i_cons, &
                                 angle1, angle2, bfield, ux, nspin_lsda, &
                                 nspin_gga, nspin_mag
#if defined(__XSD) 
  USE pw_restart_new,     ONLY : pw_readschema_file, init_vars_from_schema 
  USE qes_libs_module,    ONLY : qes_reset_output, qes_reset_parallel_info, qes_reset_general_info
  USE qes_types_module,   ONLY : output_type, parallel_info_type, general_info_type 
#else
  USE pw_restart,         ONLY : pw_readfile
#endif
  USE exx,                ONLY : ecutfock, exx_grid_init, exx_mp_init, exx_div_check
  USE funct,              ONLY : dft_is_meta, dft_is_hybrid, dft_is_gradient
  USE paw_variables,      ONLY : okpaw
  USE fcp_variables,      ONLY : lfcpopt, lfcpdyn
  USE extfield,           ONLY : monopole
  !
  IMPLICIT NONE
  !
  INTEGER  :: na, is, ierr, ibnd, ik
  LOGICAL  :: magnetic_sym, skip_equivalence=.FALSE.
  REAL(DP) :: iocc, ionic_charge, one
  !
  LOGICAL, EXTERNAL  :: check_para_diag
!
#if defined(__XSD)
  TYPE(output_type),ALLOCATABLE             :: output_obj 
  TYPE(parallel_info_type),ALLOCATABLE      :: parinfo_obj
  TYPE(general_info_type),ALLOCATABLE       :: geninfo_obj
#endif
!  
#if defined(__MPI)
  LOGICAL :: lpara = .true.
#else
  LOGICAL :: lpara = .false.
#endif

  !
  ! ... okvan/okpaw = .TRUE. : at least one pseudopotential is US/PAW
  !
  okvan = ANY( upf(1:ntyp)%tvanp )
  okpaw = ANY( upf(1:ntyp)%tpawp )
  !
  ! ... check for features not implemented with US-PP or PAW
  !
  IF ( okvan .OR. okpaw ) THEN
     IF ( dft_is_meta() ) CALL errore( 'setup', &
                               'Meta-GGA not implemented with USPP/PAW', 1 )
     IF ( noncolin .AND. lberry)  CALL errore( 'setup', &
       'Noncolinear Berry Phase/electric not implemented with USPP/PAW', 1 )
     IF  (ts_vdw ) CALL errore ('setup',&
                   'Tkatchenko-Scheffler not implemented with USPP/PAW', 1 )
     IF ( lorbm ) CALL errore( 'setup', &
                  'Orbital Magnetization not implemented with USPP/PAW', 1 )
  END IF

  IF ( dft_is_hybrid() ) THEN
     IF ( allfrac ) CALL errore( 'setup ', &
                         'option use_all_frac incompatible with hybrid XC', 1 )
     IF (.NOT. lscf) CALL errore( 'setup ', &
                         'hybrid XC not allowed in non-scf calculations', 1 )
     IF ( ANY (upf(1:ntyp)%nlcc) ) CALL infomsg( 'setup ', 'BEWARE:' // &
               & ' nonlinear core correction is not consistent with hybrid XC')
     !IF (okpaw) CALL errore('setup','PAW and hybrid XC not tested',1)
     IF (okvan) THEN
        IF (ecutfock /= 4*ecutwfc) CALL infomsg &
           ('setup','Warning: US/PAW use ecutfock=4*ecutwfc, ecutfock ignored')
        IF ( noncolin ) CALL errore &
           ('setup','Noncolinear hybrid XC for USPP not implemented',1)
     END IF
     IF ( noncolin ) no_t_rev=.true.
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
#if defined (__XSD)
  IF ( lbands .OR. ( (lfcpopt .OR. lfcpdyn ) .AND. restart )) THEN 
     ALLOCATE ( output_obj, parinfo_obj, geninfo_obj )
     CALL pw_readschema_file( ierr , output_obj, parinfo_obj, geninfo_obj )
  END IF
  !
  ! 
  IF (lfcpopt .AND. restart ) THEN  
     CALL init_vars_from_schema( 'fcpopt', ierr,  output_obj, parinfo_obj, geninfo_obj)
     tot_charge = ionic_charge - nelec
  END IF 
  IF (lfcpdyn .AND. restart ) THEN    
     CALL init_vars_from_schema( 'fcpdyn', ierr,  output_obj, parinfo_obj, geninfo_obj ) 
     tot_charge = ionic_charge - nelec 
  END IF
#else 
  IF ( lfcpopt .AND. restart ) THEN
     CALL pw_readfile( 'fcpopt', ierr )
     tot_charge = ionic_charge - nelec
  END IF
  !
  IF ( lfcpdyn .AND. restart ) THEN

     CALL pw_readfile( 'fcpdyn', ierr )
     tot_charge = ionic_charge - nelec
  END IF
#endif
  !
  ! ... magnetism-related quantities
  !
  ALLOCATE( m_loc( 3, nat ) )
  ! time reversal operation is set up to 0 by default
  t_rev = 0
  IF ( noncolin ) THEN
     !
     ! gamma_only and noncollinear not allowed
     !
     if (gamma_only) call errore('setup', &
                                 'gamma_only and noncolin not allowed',1)
     !
     ! ... wavefunctions are spinors with 2 components
     !
     npol = 2
     !
     ! ... Set the domag variable to make a spin-orbit calculation with zero
     ! ... magnetization
     !
     IF ( lspinorb ) THEN
        !
        domag = ANY ( ABS( starting_magnetization(1:ntyp) ) > 1.D-6 )
        !
     ELSE
        !
        domag = .TRUE.
        !
     END IF
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
     if (dft_is_gradient()) call compute_ux(m_loc,ux,nat)
     !
  ELSE
     !
     ! ... wavefunctions are scalars
     !
     IF (lspinorb)  CALL errore( 'setup ',  &
         'spin orbit requires a non collinear calculation', 1 )
     npol = 1
     !
     !
     IF ( i_cons == 1) then
        do na=1,nat
           m_loc(1,na) = starting_magnetization(ityp(na))
        end do
     end if
     IF ( i_cons /= 0 .AND. nspin ==1) &
        CALL errore( 'setup', 'this i_cons requires a magnetic calculation ', 1 )
     IF ( i_cons /= 0 .AND. i_cons /= 1 ) &
        CALL errore( 'setup', 'this i_cons requires a non colinear run', 1 )
  END IF
  !
  !  Set the different spin indices
  !
  nspin_mag  = nspin
  nspin_lsda = nspin
  nspin_gga  = nspin
  IF (nspin==4) THEN
     nspin_lsda=1
     IF (domag) THEN
        nspin_gga=2
     ELSE
        nspin_gga=1
        nspin_mag=1
     ENDIF
  ENDIF    
  !
  ! ... if this is not a spin-orbit calculation, all spin-orbit pseudopotentials
  ! ... are transformed into standard pseudopotentials
  !
  IF ( lspinorb .AND. ALL ( .NOT. upf(:)%has_so ) ) &
        CALL infomsg ('setup','At least one non s.o. pseudo')
  !
  IF ( .NOT. lspinorb ) CALL average_pp ( ntyp )
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
  ! ... setting nelup/neldw 
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
        CALL errore( 'setup', 'too few bands', 1 )
     !
  END IF
  !
  ! ... Here we  set the precision of the diagonalization for the first scf
  ! ... iteration of for the first ionic step
  ! ... for subsequent steps ethr is automatically updated in electrons
  !
  IF ( nat==0 ) THEN
     ethr=1.0D-8
  ELSE IF ( .NOT. lscf ) THEN
     !
     IF ( ethr == 0.D0 ) ethr = 0.1D0 * MIN( 1.D-2, tr2 / nelec )
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
           ethr = 1.D-5
           !
        ELSE
           !
           ! ... starting atomic potential is probably far from scf
           ! ... do not waste iterations in the first diagonalizations
           !
           ethr = 1.0D-2
           !
        END IF
        !
     END IF
     !
  END IF
  !
  IF ( .NOT. lscf ) niter = 1
  !
  ! ... set number of atomic wavefunctions
  !
  natomwfc = n_atom_wfc( nat, ityp, noncolin )
  !
  ! ... set the max number of bands used in iterative diagonalization
  !
  nbndx = nbnd
  IF ( isolve == 0 ) nbndx = david * nbnd
  !
#if defined(__MPI)
  use_para_diag = check_para_diag( nbnd )
#else
  use_para_diag = .FALSE.
#endif
  !
  ! ... Set the units in real and reciprocal space
  !
  tpiba  = 2.D0 * pi / alat
  tpiba2 = tpiba**2
  !
  ! ... Compute the cut-off of the G vectors
  !
  doublegrid = ( dual > 4.D0 )
  IF ( doublegrid .AND. (.NOT.okvan .AND. .not.okpaw) ) &
     CALL infomsg ( 'setup', 'no reason to have ecutrho>4*ecutwfc' )
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
  ! ... calculate dimensions of the FFT grid
  !
  ! ... if the smooth and dense grid must coincide, ensure that they do
  ! ... also if dense grid is set from input and smooth grid is not
  !
  IF ( ( dfftp%nr1 /= 0 .AND. dfftp%nr2 /= 0 .AND. dfftp%nr3 /= 0 ) .AND. &
       ( dffts%nr1 == 0 .AND. dffts%nr2 == 0 .AND. dffts%nr3 == 0 ) .AND. &
       .NOT. doublegrid ) THEN
     dffts%nr1 = dfftp%nr1
     dffts%nr2 = dfftp%nr2
     dffts%nr3 = dfftp%nr3
  END IF
  CALL fft_type_allocate ( dfftp, at, bg, gcutm, intra_bgrp_comm )
  CALL fft_type_allocate ( dffts, at, bg, gcutms, intra_bgrp_comm)
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
  ! ... If nosym is true do not use any point-group symmetry
  !
  IF ( nosym ) nrot = 1
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
        nrot  = 1
        nsym  = 1
        !
     ELSE IF (lberry ) THEN
        !
        CALL kp_strings( nppstr, gdir, nrot, s, bg, npk, &
                         k1, k2, k3, nk1, nk2, nk3, nkstot, xk, wk )
        nosym = .TRUE.
        nrot  = 1
        nsym  = 1
        !
     ELSE
        !
        CALL kpoint_grid ( nrot, time_reversal, skip_equivalence, s, t_rev, bg,&
                           npk, k1,k2,k3, nk1,nk2,nk3, nkstot, xk, wk)
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

        ! <AF>
        IF ( gdir<1 .OR. gdir>3 ) CALL errore('setup','invalid gdir value'&
                                  &' (valid values: 1=x, 2=y, 3=z)',10) 
        !
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
        nrot  = 1
        nsym  = 1
        !
     END IF
  END IF
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
     CALL find_sym ( nat, tau, ityp, magnetic_sym, m_loc, monopole )
     !
     IF ( .NOT. allfrac ) CALL remove_sym ( dfftp%nr1, dfftp%nr2, dfftp%nr3 )
     !
  END IF
  !
  ! ... Input k-points are assumed to be  given in the IBZ of the Bravais
  ! ... lattice, with the full point symmetry of the lattice.
  ! ... If some symmetries of the lattice are missing in the crystal,
  ! ... "irreducible_BZ" computes the missing k-points.
  !
  IF ( .NOT. lbands ) THEN
     CALL irreducible_BZ (nrot, s, nsym, time_reversal, &
                          magnetic_sym, at, bg, npk, nkstot, xk, wk, t_rev)
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
  IF ( lbands ) THEN
     !
     ! ... if calculating bands, we read the Fermi energy
     !
#if defined (__XSD)
     CALL init_vars_from_schema( 'ef',   ierr , output_obj, parinfo_obj, geninfo_obj)
#else
     CALL pw_readfile( 'reset', ierr )
     CALL pw_readfile( 'ef',   ierr )
#endif 
     CALL errore( 'setup ', 'problem reading ef from file ' // &
             & TRIM( tmp_dir ) // TRIM( prefix ) // '.save', ierr )

     !
  ELSE IF ( ltetra ) THEN
     !
     ! ... Calculate quantities used in tetrahedra method
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
#if defined(__XSD) 
  IF ( lbands .OR. ( (lfcpopt .OR. lfcpdyn ) .AND. restart ) ) THEN 
     CALL qes_reset_output ( output_obj ) 
     CALL qes_reset_parallel_info ( parinfo_obj ) 
     CALL qes_reset_general_info ( geninfo_obj ) 
     DEALLOCATE ( output_obj, parinfo_obj, geninfo_obj ) 
  END IF 
#endif
  !
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
  ! ... distribute k-points (and their weights and spin indices)
  !
  kunit = 1
  CALL divide_et_impera ( nkstot, xk, wk, isk, nks )
  !
  IF ( dft_is_hybrid() ) THEN
     CALL exx_grid_init()
     CALL exx_mp_init()
     CALL exx_div_check()
  ENDIF

  IF (one_atom_occupations) THEN
     DO ik=1,nkstot
        DO ibnd=natomwfc+1, nbnd
           IF (f_inp(ibnd,ik)> 0.0_DP) CALL errore('setup', &
               'no atomic wavefunction for some band',1)
        ENDDO
     ENDDO
  ENDIF
  !
  ! ... Set up Hubbard parameters for LDA+U calculation
  !
  CALL init_lda_plus_u ( upf(1:ntyp)%psd, noncolin )
  !
  ! ... initialize d1 and d2 to rotate the spherical harmonics
  !
  IF (lda_plus_u .or. okpaw .or. (okvan.and.dft_is_hybrid()) ) CALL d_matrix( d1, d2, d3 )
  !
  RETURN
  !
END SUBROUTINE setup
!
!----------------------------------------------------------------------------
LOGICAL FUNCTION check_para_diag( nbnd )
  !
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE mp_diag,          ONLY : np_ortho, ortho_parent_comm 
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp_pools,         ONLY : intra_pool_comm
  USE control_flags,    ONLY : gamma_only

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nbnd
  LOGICAL, SAVE :: first = .TRUE.
  LOGICAL, SAVE :: saved_value = .FALSE.

  IF( .NOT. first ) then
      check_para_diag = saved_value
      RETURN
  end if
  first = .FALSE.
  !
  IF( np_ortho(1) > nbnd ) &
     CALL errore ('check_para_diag', 'Too few bands for required ndiag',nbnd)
  !
  check_para_diag = ( np_ortho( 1 ) > 1 .AND. np_ortho( 2 ) > 1 )
  saved_value = check_para_diag
  !
  IF ( ionode ) THEN
     !
     WRITE( stdout, '(/,5X,"Subspace diagonalization in iterative solution ",&
                     &     "of the eigenvalue problem:")' ) 
     IF ( check_para_diag ) THEN
        IF (ortho_parent_comm .EQ. intra_pool_comm) THEN
           WRITE( stdout, '(5X,"one sub-group per k-point group (pool) will be used")' )
        ELSE IF (ortho_parent_comm .EQ. intra_bgrp_comm) THEN
           WRITE( stdout, '(5X,"one sub-group per band group will be used")' )
        ELSE
           CALL errore( 'setup','Unexpected sub-group communicator ', 1 )
        END IF
#if defined(__ELPA) || defined(__ELPA_2015) || defined(__ELPA_2016)
        WRITE( stdout, '(5X,"ELPA distributed-memory algorithm ", &
              & "(size of sub-group: ", I2, "*", I3, " procs)",/)') &
               np_ortho(1), np_ortho(2)
#elif defined(__SCALAPACK)
        WRITE( stdout, '(5X,"scalapack distributed-memory algorithm ", &
              & "(size of sub-group: ", I2, "*", I3, " procs)",/)') &
               np_ortho(1), np_ortho(2)
#else
        WRITE( stdout, '(5X,"custom distributed-memory algorithm ", &
              & "(size of sub-group: ", I2, "*", I3, " procs)",/)') &
               np_ortho(1), np_ortho(2)
#endif
     ELSE
        WRITE( stdout, '(5X,"a serial algorithm will be used",/)' )
     END IF
     !
  END IF
  !
  RETURN
END FUNCTION check_para_diag
