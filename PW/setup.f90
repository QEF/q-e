!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE setup()
  !----------------------------------------------------------------------------
  !
  ! ... This routine
  ! ... 1) determines various parameters of the calculation
  ! ... 2) finds actual crystal symmetry, determine lattice
  ! ... 3) generates k-points corresponding to the crystal symmetry
  !
  ! ... Calculated parameters:
  ! ...   zv        charge of each atomic type
  ! ...   nelec     total number of electrons
  ! ...   nbnd      total number of bands
  ! ...   nbndx     max number of bands used in iterative diagonalization
  ! ...   tpiba     2 pi / a (a = lattice parameter)
  ! ...   tpiba2    square of tpiba
  ! ...   gcutm     cut-off in g space
  ! ...   gcutms    cut-off in g space for smooth functions
  ! ...   ethr      convergence limit of iterative diagonalization
  ! ...   at        direct lattice vectors
  ! ...   omega     volume of the unit cell
  ! ...   bg        reciprocal lattice vectors
  ! ...   s         symmetry matrices in the direct lattice vectors basis
  ! ...   nsym      total number of symmetry operations
  ! ...   ftau      fractionary translations
  ! ...   irt       for each atom gives the corresponding symmetric
  ! ...   invsym    if true the system has inversion symmetry
  ! ... + non-collinear related quantities
  ! ... + spin-orbit related quantities
  ! ... + LDA+U-related quantities
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : eps8
  USE parameters,         ONLY : npk
  USE io_global,          ONLY : stdout
  USE io_files,           ONLY : tmp_dir, prefix, delete_if_present
  USE constants,          ONLY : pi, degspin
  USE cell_base,          ONLY : at, bg, alat, tpiba, tpiba2, ibrav, &
                                 symm_type, omega
  USE ions_base,          ONLY : nat, tau, ntyp => nsp, ityp, zv
  USE basis,              ONLY : startingpot, natomwfc
  USE gvect,              ONLY : gcutm, ecutwfc, dual, nr1, nr2, nr3
  USE gsmooth,            ONLY : doublegrid, gcutms
  USE klist,              ONLY : xk, wk, xqq, nks, nelec, degauss, lgauss, &
                                 lxkcry, nkstot, &
                                 nelup, neldw, two_fermi_energies, &
                                 tot_charge, tot_magnetization, multiplicity
  USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk, &
                                 starting_magnetization
  USE ener,               ONLY : ef
  USE electrons_base,     ONLY : set_nelup_neldw
  USE ktetra,             ONLY : nk1, nk2, nk3, k1, k2, k3, &
                                 tetra, ntetra, ltetra
  USE symme,              ONLY : s, t_rev, irt, ftau, nsym, invsym
  USE wvfct,              ONLY : nbnd, nbndx, gamma_only
  USE control_flags,      ONLY : tr2, ethr, lscf, lmd, lpath, lphonon, david,  &
                                 isolve, niter, noinv, nosym, modenum, lbands, &
                                 use_para_diago, use_distpara_diago
  USE control_flags,      ONLY : para_diago_dim   ! debug
  USE relax,              ONLY : starting_diag_threshold
  USE cellmd,             ONLY : calc
  USE uspp_param,         ONLY : zp, psd, nbeta, jjj, tvanp
  USE uspp,               ONLY : okvan
  USE ldaU,               ONLY : d1, d2, d3, lda_plus_u, Hubbard_U, &
                                 Hubbard_l, Hubbard_alpha, Hubbard_lmax
  USE bp,                 ONLY : gdir, lberry, nppstr
  USE fixed_occ,          ONLY : f_inp, tfixed_occ
  USE char,               ONLY : sname
  USE funct,              ONLY : set_dft_from_name
  USE mp_global,          ONLY : nimage, kunit
  USE spin_orb,           ONLY : lspinorb, domag, so
  USE noncollin_module,   ONLY : noncolin, npol, m_loc, i_cons, mcons, &
                                 angle1, angle2, bfield
  USE pw_restart,         ONLY : pw_readfile
  USE input_parameters,   ONLY : restart_mode
#if defined (EXX)
  USE exx,                ONLY : exx_grid_init
#endif
  USE funct,              ONLY : dft_is_meta, dft_is_hybrid
  !
  IMPLICIT NONE
  !
  INTEGER  :: na, nt, input_nks, nrot, irot, isym, tipo, is, nb, ierr, ibnd
  LOGICAL  :: minus_q, ltest
  REAL(DP) :: iocc, ionic_charge
  !
  INTEGER, EXTERNAL :: n_atom_wfc, set_Hubbard_l
  LOGICAL, EXTERNAL :: lchk_tauxk ! tests that atomic coordinates do not overlap
  !
  !
  IF (dft_is_meta()) THEN
     DO nt=1,ntyp
        IF ( tvanp(nt) ) &
             CALL errore( 'setup', 'US and Meta-GGA not yet implemented', 1 )
     END DO
  END IF

  ALLOCATE( m_loc( 3, nat ) )
  !
  IF ( nimage > 1 .AND. .NOT. lpath ) &
     CALL errore( 'setup', 'images parallelization not permitted', 1 )
  !
  ! ... Compute the ionic charge for each atom type
  !
  zv(1:ntyp) = zp(1:ntyp)
  !
#if defined (__PGI)
     ionic_charge = - tot_charge
     DO na = 1, nat
        ionic_charge = ionic_charge + zv( ityp(na) )
     END DO
#else
     ionic_charge = SUM( zv(ityp(1:nat)) ) - tot_charge
#endif
  !
  IF ( nelec == 0.D0 ) THEN
     !
     ! ... set the number of electrons equal to the total ionic charge
     !
     nelec = ionic_charge
     !
  END IF
  !
  ! ... If the occupations are from input, check the consistency with the
  ! ... number of electrons
  !
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
     ! ... transform angles to radiants
     !
!     DO nt = 1, ntyp
        !
!        angle1(nt) = pi * angle1(nt) / 180.D0
!        angle2(nt) = pi * angle2(nt) / 180.D0
!        !
!     END DO
     !
     ! ... Set the nomag variable to make a spin-orbit calculation with zero
     ! ... magnetization
     !
     IF ( lspinorb ) THEN
        !
        domag = .FALSE.
        !
        DO nt = 1, ntyp
           !
           domag = domag .OR. ( ABS( starting_magnetization(nt) ) > 1.D-6 )
           !
        END DO
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
     bfield=0.D0
     !
     IF ( i_cons == 2 ) THEN
        !
        ! ... angle theta between the magnetic moments and the z-axis is
        ! ... constrained. Transform theta to radiants
        !
        DO na = 1, ntyp
           !
           mcons(1,na) = pi * mcons(1,na) / 180.D0
           !
        END DO
        !
     ELSE IF ( i_cons == 4 ) THEN
        !
        bfield(:) = mcons(:,1)
        !
     END IF
     !
  ELSE
     !
     ! ... wavefunctions are scalars
     !
     npol = 1
     !
     IF ( i_cons == 5 ) THEN
        !
        nelup = ( nelec + mcons(3,1) ) * 0.5D0
        neldw = ( nelec - mcons(3,1) ) * 0.5D0
        !
     ENDIF
     !
     IF ( i_cons /= 0 .AND. i_cons /= 5 ) &
        CALL errore( 'setup', 'this i_cons requires a non colinear run', 1 )
     IF ( i_cons == 5 .AND. nspin /= 2 ) &
        CALL errore( 'setup', 'i_cons can be 5 only with nspin=2', 1 )
  END IF
  !
  IF ( tfixed_occ ) THEN
     !
     iocc = 0
     !
     IF ( noncolin ) THEN
        !
#if defined (__PGI)
        DO ibnd = 1, nbnd
           iocc = iocc + f_inp(ibnd,1)
        END DO
#else
        iocc = iocc + SUM( f_inp(1:nbnd,1) )
#endif
        DO ibnd = 1, nbnd
           if (f_inp(ibnd,1).gt.1.d0.or.f_inp(ibnd,1).lt.0.d0) call &
              errore('setup','wrong fixed occupations',1)
        END DO
        !
     ELSE
        !
        DO is = 1, nspin
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
              if (f_inp(ibnd,is) > 2.d0/nspin .or. f_inp(ibnd,is) < 0.d0) &
                   call errore('setup','wrong fixed occupations',1)
           END DO
        END DO
        !
     END IF
     !
     IF ( ABS( iocc - nelec ) > 1D-5 ) &
        CALL errore( 'setup', 'strange occupations', 1 )
     !
  END IF
  !
  ! ... For metals: check whether Gaussian broadening or Tetrahedron method
  ! ...             is used
  !
  lgauss = ( ( degauss /= 0.D0 ) .AND. ( .NOT. tfixed_occ ) )
  !
  ! ... Check: if there is an odd number of electrons, the crystal is a metal
  !
  IF ( lscf .AND. ABS( NINT( nelec / 2.D0 ) - nelec / 2.D0 ) > eps8 &
            .AND. .NOT. lgauss .AND. .NOT. ltetra .AND. .NOT. tfixed_occ ) &
      CALL infomsg( 'setup', 'the system is metallic, specify occupations' )
  !
  ! ... Check: spin-polarized calculations require tetrahedra or broadening
  !            or fixed occupation - the simple filling of levels is not
  !            implemented right now (it will yield an unpolarized system)
  !
  IF ( lscf .AND. lsda &
            .AND. .NOT. lgauss .AND. .NOT. ltetra &
            .AND. .NOT. tfixed_occ .AND. .NOT. two_fermi_energies ) &
      CALL errore( 'setup', 'spin-polarized system, specify occupations', 1 )
  !
  ! ... setting nelup/neldw if not set in the input
  !
  call set_nelup_neldw(nelec, nelup, neldw, tot_magnetization, &
       multiplicity)
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
  ltest = ( ethr == 0.D0 )
  !
  IF ( lphonon ) THEN
     !
     ! ... in the case of a phonon calculation ethr can not be specified
     ! ... in the input file
     !
     IF ( .NOT. ltest ) &
        WRITE( UNIT = stdout, &
             & FMT = '(5X,"diago_thr_init overwritten ", &
             &            "with conv_thr / nelec")' )
     !
     ethr = 0.1D0 * MIN( 1.D-2, tr2 / nelec )
     !
  ELSE IF ( .NOT. lscf ) THEN
     !
     IF ( ltest ) THEN
        !
        ethr = 0.1D0 * MIN( 1.D-2, tr2 / nelec )
        !
     END IF
     !
  ELSE
     !
     IF ( ltest ) THEN
        !
        IF ( startingpot == 'file' ) THEN
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
  IF (nat==0) THEN
     ethr=1.0D-8
!
!  In this case, print the Hatree-Fock energy of free electrons per cell
!  (not per electron).
!
     CALL set_dft_from_name('sla-noc-nogx-nogc')
  END IF
  !
  IF ( .NOT. lscf ) niter = 1
  !
  starting_diag_threshold = ethr
  !
  !   check if spin-orbit is possible
  !
  IF (lspinorb.and..not.noncolin)  CALL errore( 'setup ',  &
         'spin orbit requires a non collinear calculation', ierr )
  !
  ! ... if this is not a spin-orbit calculation, all spin-orbit pseudopotentials
  ! ... are transformed into standard pseudopotentials
  !
  DO nt = 1, ntyp
     !
     so(nt) = ( nbeta(nt) > 0 )
     !
     IF ( lspinorb .AND. ALL ( ABS( jjj(1:nbeta(nt),nt) ) < 1.D-7 ) ) &
        CALL infomsg ('setup','At least one non s.o. pseudo')
     !
     DO nb = 1, nbeta(nt)
        !
        so(nt) = so(nt) .AND. ( ABS( jjj(nb,nt) ) > 1.D-7 )
        !
     END DO
     !
  END DO
  !
  IF ( .NOT. lspinorb ) CALL average_pp ( ntyp )
  !
  ! ... set number of atomic wavefunctions
  !
  natomwfc = n_atom_wfc( nat, ityp )
  !
  ! ... set the max number of bands used in iterative diagonalization
  !
  nbndx = nbnd
  !
  IF ( isolve == 0 ) nbndx = david * nbnd
  !
  ! ... Set the units in real and reciprocal space
  !
  tpiba  = 2.D0 * pi / alat
  tpiba2 = tpiba**2
  !
  ! ... Compute the cut-off of the G vectors
  !
  gcutm = dual * ecutwfc / tpiba2
  !
  doublegrid = ( dual > 4.D0 )
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
  ! ... If  lxkcry = .TRUE. , the input k-point components in crystal
  ! ... axis are transformed in cartesian coordinates
  !
  IF ( lxkcry ) CALL cryst_to_cart( nkstot, xk, bg, 1 )
  !
  ! ... Test that atoms do not overlap
  !
  IF ( .NOT. ( lchk_tauxk( nat, tau, bg ) ) ) &
     CALL errore( 'setup', 'Wrong atomic coordinates ', 1 )
  !
  ! ... calculate dimensions of the FFT grid
  !
  CALL set_fft_dim()
  !
  !  ... generate transformation matrices for the crystal point group
  !  ... First we generate all the symmetry matrices of the Bravais lattice
  !
  IF ( ibrav == 4 .OR. ibrav == 5 ) THEN
     !
     ! ... here the hexagonal or trigonal bravais lattice
     !
     CALL hexsym( at, s, sname, nrot )
     !
     tipo = 2
     !
  ELSE IF ( ibrav >=1  .AND. ibrav <= 14 ) THEN
     !
     ! ... here for the cubic bravais lattice
     !
     CALL cubicsym( at, s, sname, nrot )
     !
     tipo = 1
     !
  ELSE IF ( ibrav == 0 ) THEN
     !
     IF ( symm_type == 'cubic' ) THEN
        !
        tipo = 1
        !
        CALL cubicsym( at, s, sname, nrot )
        !
     ELSE IF ( symm_type == 'hexagonal' ) THEN
        !
        tipo = 2
        !
        CALL hexsym( at, s, sname, nrot )
        !
     END IF
     !
  ELSE
     !
     CALL errore( 'setup', 'wrong ibrav', 1 )
     !
  END IF
  !
  ! ... if noinv is .TRUE. eliminate all symmetries which exchange z with -z
  !
  IF ( noinv ) THEN
     !
     irot = 0
     !
     DO isym = 1, nrot
        IF ( s(1,3,isym) == 0 .AND. s(3,1,isym) == 0 .AND. &
             s(2,3,isym) == 0 .AND. s(3,2,isym) == 0 .AND. &
             s(3,3,isym) == 1) THEN
           !
           irot = irot + 1
           !
           s(:,:,irot) = s(:,:,isym)
           !
           sname(irot) = sname(isym)
           !
        END IF
        !
     END DO
     !
     nrot = irot
     !
  END IF
  !
  ! ... If nosym is true do not use any point-group symmetry
  !
  IF ( nosym ) nrot = 1
  !
  ! ... Automatic generation of k-points (if required)
  !
  IF ( nkstot < 0 ) THEN
     !
     CALL setupkpoint( s, nrot, xk, wk, nkstot, npk, nk1, &
                       nk2, nk3, k1, k2, k3, at, bg, tipo )
     !
  ELSE IF ( nkstot == 0 ) THEN
     !
     IF ( lberry ) THEN
        !
        CALL kp_strings( nppstr, gdir, nrot, s, bg, npk, &
                         k1, k2, k3, nk1, nk2, nk3, nkstot, xk, wk )
        !
        nosym = .TRUE.
        nrot  = 1
        nsym  = 1
        !
     ELSE
        !
        CALL kpoint_grid( nrot, s, bg, npk, k1, k2, k3, &
                          nk1, nk2, nk3, nkstot, xk, wk )
        !
     END IF
     !
  END IF
  !
  ! ...  allocate space for irt
  !
  ALLOCATE( irt( 48, nat ) )
  !
  ! ... "sgama" eliminates rotations that are not symmetry operations
  ! ... Input k-points are assumed to be  given in the IBZ of the Bravais
  ! ... lattice, with the full point symmetry of the lattice.
  ! ... If some symmetries of the lattice are missing in the crystal,
  ! ... "sgama" computes the missing k-points.
  !
  input_nks = nkstot
  !
  IF (nat>0) THEN
     CALL sgama( nrot, nat, s, sname, t_rev, at, bg, tau, ityp, nsym, nr1,&
          nr2, nr3, irt, ftau, npk, nkstot, xk, wk, invsym, minus_q,  &
          xqq, modenum, noncolin, domag, m_loc )
  ELSE
     nsym=nrot
     invsym=.true.
  ENDIF
  !
  IF (nat>0) &
     CALL checkallsym( nsym, s, nat, tau, ityp, at, &
          bg, nr1, nr2, nr3, irt, ftau, alat, omega )
  !
  ! ... if dynamics is done the system should have no symmetries
  ! ... (inversion symmetry alone is allowed)
  !
  IF ( lmd .AND. ( nsym == 2 .AND. .NOT. invsym .OR. nsym > 2 ) &
           .AND. .NOT. ( calc == 'mm' .OR. calc == 'nm' ) ) &
       CALL infomsg( 'setup', 'Dynamics, you should have no symmetries' )
  !
  ntetra = 0
  !
  IF ( lbands ) THEN
     !
     ! ... if calculating bands, we leave k-points unchanged and read the
     ! Fermi energy
     !
     nkstot = input_nks
     CALL pw_readfile( 'reset', ierr )
     CALL pw_readfile( 'ef',   ierr )
     CALL errore( 'setup ', 'problem reading ef from file ' // &
             & TRIM( tmp_dir ) // TRIM( prefix ) // '.save', ierr )

     IF ( restart_mode == 'from_scratch' ) &
           CALL delete_if_present( TRIM( tmp_dir ) // TRIM( prefix ) &
                                     //'.save/data-file.xml' )
     !
  ELSE IF ( ltetra ) THEN
     !
     ! ... Calculate quantities used in tetrahedra method
     !
     ntetra = 6 * nk1 * nk2 * nk3
     !
     ALLOCATE( tetra( 4, ntetra ) )
     !
     CALL tetrahedra( nsym, s, minus_q, at, bg, npk, k1, k2, k3, &
          nk1, nk2, nk3, nkstot, xk, wk, ntetra, tetra )
     !
  END IF
  !
  ! ... phonon calculation: add k+q to the list of k
  !
  IF ( lphonon ) CALL set_kplusq( xk, wk, xqq, nkstot, npk )
  !
#if defined (EXX)
  IF ( dft_is_hybrid() ) CALL exx_grid_init()
#else
  IF ( dft_is_hybrid() ) CALL errore( 'setup ', &
                         'HYBRID XC not implemented in PWscf', 1 )
#endif

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
     !
  ELSE
     !
     ! ... LDA case: the two spin polarizations are identical
     !
     wk(1:nkstot)    = wk(1:nkstot) * degspin
     current_spin = 1
     !
     IF ( nspin /= 1 ) &
        CALL errore( 'setup', 'nspin should be 1; check iosys', 1 )
     !
  END IF
  !
  IF ( nkstot > npk ) CALL errore( 'setup', 'too many k points', nkstot )
  !
#ifdef __PARA
  !
  ! ... set the granularity for k-point distribution
  !
  IF ( ( ABS( xqq(1) ) < eps8 .AND. ABS( xqq(2) ) < eps8 .AND. &
         ABS( xqq(3) ) < eps8) .OR. ( .NOT. lphonon ) ) THEN
     !
     kunit = 1
     !
  ELSE
     !
     kunit = 2
     !
  ENDIF
  !
  ! ... distribute k-points (and their weights and spin indices)
  !
  CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
  !
  IF ( use_para_diago )      CALL check_para_diag_efficiency()
  !
  IF ( use_distpara_diago )  CALL check_distpara_availability()
  !
#else
  !
  nks = nkstot
  !
#endif
  !
  ! ... okvan = .TRUE. : at least one pseudopotential is US
  !
  okvan = ANY( tvanp(1:ntyp) )
  !
  ! ... Needed for LDA+U
  !
  ! ... initialize d1 and d2 to rotate the spherical harmonics
  !
  IF ( lda_plus_u ) THEN
     !
     Hubbard_lmax = -1
     !
     DO nt = 1, ntyp
        !
        IF ( Hubbard_U(nt) /= 0.D0 .OR. Hubbard_alpha(nt) /= 0.D0 ) THEN
           !
           Hubbard_l(nt) = set_Hubbard_l( psd(nt) )
           !
           Hubbard_lmax = MAX( Hubbard_lmax, Hubbard_l(nt) )
           !
           WRITE( UNIT = stdout, &
                  FMT = * ) ' HUBBARD L FOR TYPE ',psd(nt),' IS ', Hubbard_l(nt)
           !
        END IF
        !
     END DO
     !
     WRITE( UNIT = stdout, &
            FMT = * ) ' MAXIMUM HUBBARD L IS ', Hubbard_lmax
     !
     IF ( Hubbard_lmax == -1 ) &
        CALL errore( 'setup', &
                   & 'lda_plus_u calculation but Hubbard_l not set', 1 )
     !
     CALL d_matrix( d1, d2, d3 )
     !
  ELSE
     !
     Hubbard_lmax = 0
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE setup
!
!----------------------------------------------------------------------------
FUNCTION n_atom_wfc( nat, ityp )
  !----------------------------------------------------------------------------
  !
  ! ... Find number of starting atomic orbitals
  !
  USE atom,             ONLY : oc, chi, nchi, lchi, jchi
  USE noncollin_module, ONLY : noncolin
  USE spin_orb,         ONLY : so
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: nat, ityp(nat)
  !
  INTEGER  :: n_atom_wfc
  !
  INTEGER  :: na, nt, n
  !
  !
  n_atom_wfc = 0
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     DO n = 1, nchi(nt)
        !
        IF ( oc(n,nt) >= 0.D0 ) THEN
           !
           IF ( noncolin ) THEN
              !
              IF ( so(nt) ) THEN
                 !
                 n_atom_wfc = n_atom_wfc + 2 * lchi(n,nt)
                 !
                 IF ( ABS( jchi(n,nt) - lchi(n,nt) - 0.5D0 ) < 1.D-6 ) &
                    n_atom_wfc = n_atom_wfc + 2
                 !
              ELSE
                 !
                 n_atom_wfc = n_atom_wfc + 2 * ( 2 * lchi(n,nt) + 1 )
                 !
              END IF
              !
           ELSE
              !
              n_atom_wfc = n_atom_wfc + 2 * lchi(n,nt) + 1
              !
           END IF
        END IF
     END DO
  END DO
  !
  RETURN
  !
END FUNCTION n_atom_wfc
!
!----------------------------------------------------------------------------
SUBROUTINE check_para_diag_efficiency()
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE wvfct,            ONLY : nbnd, nbndx, gamma_only
  USE control_flags,    ONLY : use_para_diago, para_diago_dim, isolve, ortho_para
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE random_numbers,   ONLY : rranf
  USE mp_global,        ONLY : nproc_pool, me_pool, intra_pool_comm, init_ortho_group
  USE mp,               ONLY : mp_bcast, mp_barrier
  !
  IMPLICIT NONE
  !
  INTEGER                  :: dim1, dim_pool, i, j, m, m_min
  REAL(DP)                 :: time_para, time_serial, t_incr, t_incr_old
  LOGICAL                  :: lfirst
  REAL(DP),    ALLOCATABLE :: ar(:,:), br(:,:), vr(:,:)
  COMPLEX(DP), ALLOCATABLE :: ac(:,:), bc(:,:), vc(:,:)
  REAL(DP),    ALLOCATABLE :: e(:)
  !
  REAL(DP), EXTERNAL :: scnds
  !
  use_para_diago = .FALSE.
  !
  !  here we initialize the sub group of processors that will take part
  !  in the matrix diagonalization. 
  !  NOTE that the maximum number of processors could not be the optimal one,
  !  and ortho_para input keyword can be used to force a given number of proc
  !
  IF( ortho_para < 1 ) THEN
     ! 
     !  use all the available processors
     !
     CALL init_ortho_group( nproc_pool, intra_pool_comm )
     !
  ELSE
     !
     CALL init_ortho_group( ortho_para, intra_pool_comm )
     !
  END IF
  !
  IF ( isolve /= 0 .OR. nproc_pool == 1 ) RETURN
  !
  IF ( ionode ) &
     WRITE( stdout, '(/,5X,"looking for ", &
                     &     "the optimal diagonalization algorithm ...",/)' )
  !
  m_min = ( nbnd / nproc_pool )*nproc_pool
  !
  m = ( 100 / nproc_pool )*nproc_pool
  !
  IF ( m > nbndx .OR. nbndx < 200 )  THEN
     !
     IF ( ionode ) WRITE( stdout, '(5X,"a serial algorithm will be used",/)' )
     !
     RETURN
     !
  END IF
  !
  IF ( ionode ) &
     WRITE( stdout, '(5X,"dimension   time para (sec)   time serial (sec)")' )
  !
  lfirst = .TRUE.
  !
  dim_loop: DO dim1 = m_min, nbndx, m
     !
     dim_pool = dim1 / nproc_pool
     !
     ALLOCATE( e( dim1 ) )
     !
     IF ( gamma_only ) THEN
        !
        ALLOCATE( ar( dim1, dim1 ) )
        ALLOCATE( br( dim1, dim1 ) )
        ALLOCATE( vr( dim1, dim1 ) )
        !
        ar(:,:) = 0.D0
        br(:,:) = 0.D0
        !
        DO i = me_pool*dim_pool + 1, ( me_pool + 1 )*dim_pool
           !
           DO j = i, dim1
              !
              ar(i,j) = rranf() - 0.5D0
              ar(j,i) = ar(i,j)
              !
           END DO
           !
        END DO
        !
        CALL reduce( dim1*dim1, ar )
        !
        FORALL( i = 1:dim1 ) br(i,i) = 1.D0
        !
     ELSE
        !
        ALLOCATE( ac( dim1, dim1 ) )
        ALLOCATE( bc( dim1, dim1 ) )
        ALLOCATE( vc( dim1, dim1 ) )
        !
        ac(:,:) = ( 0.D0, 0.D0 )
        bc(:,:) = ( 0.D0, 0.D0 )
        !
        DO i = me_pool*dim_pool + 1, ( me_pool + 1 )*dim_pool
           !
           DO j = i, dim1
              !
              ac(i,j) = CMPLX( rranf() - 0.5D0, rranf() - 0.5D0 )
              ac(j,i) = ac(i,j)
              !
           END DO
           !
        END DO
        !
        CALL reduce( 2*dim1*dim1, ac )
        !
        FORALL( i = 1:dim1 ) bc(i,i) = ( 1.D0, 0.D0 )
        !
     END IF
     !
     CALL mp_barrier()
     !
     ! ... parallel diagonalizer
     !
     use_para_diago = .TRUE.
     para_diago_dim = 1
     !
     IF ( ionode ) time_para = scnds()
     !
     CALL do_test()
     !
     IF ( ionode ) &
        time_para = ( scnds() - time_para ) / DBLE( nbndx / dim1 )
     !
     ! ... serial diagonalizer
     !
     use_para_diago = .FALSE.
     !
     IF ( ionode ) time_serial = scnds()
     !
     CALL do_test()
     !
     IF ( ionode ) &
        time_serial = ( scnds() - time_serial ) / DBLE( nbndx / dim1 )
     !
     IF ( gamma_only ) THEN
        DEALLOCATE( ar, br, vr, e )
     ELSE
        DEALLOCATE( ac, bc, vc, e )
     END IF
     !
     IF ( ionode ) &
        WRITE( stdout, '(5X,I5,2(6X,F12.8))' ) dim1, time_para, time_serial
     !
     CALL mp_bcast( time_para,   ionode_id )
     CALL mp_bcast( time_serial, ionode_id )
     !
     t_incr = ( time_para - time_serial ) / time_serial
     !
     IF ( time_para < time_serial ) THEN
        !
        use_para_diago = .TRUE.
        para_diago_dim = dim1
        !
        EXIT dim_loop
        !
     ELSE IF ( .NOT.lfirst .AND. t_incr > t_incr_old ) THEN
        !
        ! ... the parallel diagonalizer is getting slower and slower
        !
        use_para_diago = .FALSE.
        !
        EXIT dim_loop
        !
     END IF
     !
     lfirst = .FALSE.
     !
     t_incr_old = t_incr
     !
  END DO dim_loop
  !
  ! use_para_diago = .TRUE.  !  debug
  ! para_diago_dim = 1       !  debug
  !
  IF ( ionode ) THEN
     !
     IF ( use_para_diago ) THEN
        WRITE( stdout, '(/,5X,"a parallel algorithm will be used for ", &
                        &     "matrices larger than ",I4,/)' ) para_diago_dim
     ELSE
        WRITE( stdout, '(/,5X,"a serial algorithm will be used",/)' )
     END IF
     !
  END IF
  !
  RETURN
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE do_test()
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      DO i = 1, nbndx / dim1
         !
         IF ( gamma_only ) THEN
            CALL rdiaghg( dim1, dim1, ar, br, dim1, e, vr )
         ELSE
            CALL cdiaghg( dim1, dim1, ac, bc, dim1, e, vc )
         END IF
      END DO
      !
      CALL mp_barrier()
      !
    END SUBROUTINE do_test
    !
END SUBROUTINE check_para_diag_efficiency


SUBROUTINE check_distpara_availability()
  !
  USE wvfct,            ONLY : nbnd, nbndx, gamma_only
  USE control_flags,    ONLY : use_distpara_diago, isolve, ortho_para, &
                               use_para_diago, para_diago_dim
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE mp_global,        ONLY : nproc_pool, init_ortho_group, np_ortho, intra_pool_comm

  IMPLICIT NONE

  use_distpara_diago = .TRUE.
  !
  !  here we initialize the sub group of processors that will take part
  !  in the matrix diagonalization. 
  !  NOTE that the maximum number of processors could not be the optimal one,
  !  and ortho_para input keyword can be used to force a given number of proc
  !
  IF( ortho_para < 1 ) THEN
     ! 
     !  use all the available processors
     !
     CALL init_ortho_group( nproc_pool, intra_pool_comm )
     !
  ELSE
     !
     CALL init_ortho_group( ortho_para, intra_pool_comm )
     !
  END IF

  IF ( isolve /= 0 .OR. nproc_pool == 1 ) RETURN

  IF ( ionode ) THEN
     !
     WRITE( stdout, '(/,5X,"Iterative solution of the eigenvalue problem")' ) 
  END IF

  !  we need at least 4 procs to use distributed algorithm

  IF( np_ortho( 1 ) == 1 .AND. np_ortho( 2 ) == 1 ) THEN
     use_distpara_diago = .FALSE.
     IF ( ionode ) WRITE( stdout, '(5X,"Too few procs, we need at least 4 procs")' )
  END IF

  !  we need to have at least 1 electronic band per block

  IF( np_ortho( 1 ) > nbnd ) THEN
     use_distpara_diago = .FALSE.
     IF ( ionode ) WRITE( stdout, '(5X,"Too few bands, we need at least as many bands as SQRT(nproc)")' )
  END IF

  IF( use_distpara_diago ) THEN
     ! 
     ! turn on parallel diagonalization, for all matrixes
     !
     use_para_diago = .TRUE.  
     para_diago_dim = 1       
     !
  END IF

  IF ( ionode ) THEN
     !
     IF ( use_distpara_diago ) THEN
        WRITE( stdout, '(/,5X,"a parallel distributed memory algorithm will be used",/)' ) 
     ELSE
        WRITE( stdout, '(/,5X,"a serial algorithm will be used",/)' )
     END IF
     !
  END IF

  RETURN
END SUBROUTINE check_distpara_availability
