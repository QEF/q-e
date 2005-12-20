!
! Copyright (C) 2001-2004 PWSCF group
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
  ! ...   msh       mesh point (atomic grid) for which R(msh) > Rcut = 10a.u.
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
  USE parameters,         ONLY : npsx, nchix, npk
  USE io_global,          ONLY : stdout
  USE constants,          ONLY : pi, degspin
  USE cell_base,          ONLY : at, bg, alat, tpiba, tpiba2, ibrav, symm_type
  USE ions_base,          ONLY : nat, tau, ntyp => nsp, ityp, zv
  USE basis,              ONLY : startingpot, natomwfc
  USE gvect,              ONLY : gcutm, ecutwfc, dual, nr1, nr2, nr3
  USE gsmooth,            ONLY : doublegrid, gcutms
  USE klist,              ONLY : xk, wk, xqq, nks, nelec, degauss, lgauss, &
                                 lxkcry, nkstot, b_length, lcart, &
                                 nelup, neldw, two_fermi_energies
  USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk, &
                                 starting_magnetization
  USE ktetra,             ONLY : nk1, nk2, nk3, k1, k2, k3, &
                                 tetra, ntetra, ltetra
  USE symme,              ONLY : s, t_rev, irt, ftau, nsym, invsym
  USE atom,               ONLY : r, oc, chi, nchi, lchi, jchi, mesh, msh
  USE pseud,              ONLY : zp
  USE wvfct,              ONLY : nbnd, nbndx, gamma_only
  USE control_flags,      ONLY : tr2, ethr, alpha0, beta0, lscf, lmd, lpath, &
                                 lphonon, david, isolve, niter, noinv, nosym, &
                                 modenum, lraman
  USE relax,              ONLY : starting_diag_threshold
  USE cellmd,             ONLY : calc
  USE uspp_param,         ONLY : psd, betar, nbeta, dion, jjj, lll, tvanp
  USE uspp,               ONLY : okvan
  USE ldaU,               ONLY : d1, d2, d3, lda_plus_u, Hubbard_U, Hubbard_l, &
                                 Hubbard_alpha, Hubbard_lmax
  USE bp,                 ONLY : gdir, lberry, nppstr
  USE fixed_occ,          ONLY : f_inp, tfixed_occ   
  USE char,               ONLY : sname
  USE mp_global,          ONLY : nimage, kunit
  USE spin_orb,           ONLY : lspinorb, domag
  USE noncollin_module,   ONLY : noncolin, npol, m_loc, i_cons, mcons, &
                                 angle1, angle2, bfield
#if defined (EXX)
  USE exx,                ONLY : exx_grid_init
  USE funct,              ONLY : dft_is_hybrid
#endif
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: &
      rcut = 10.D0,   &!  cut-off radius for radial integrations
      eps  = 1.0D-12   !  small number
  INTEGER :: & 
      na,             &!
      ir,             &!
      nt,             &!
      input_nks,      &!
      nrot,           &!
      iter,           &!
      ierr,           &!
      irot,           &!
      isym,           &!
      ipol,           &!
      jpol,           &!
      tipo,           &!
      is,             &!
      nb,             &!
      nbe,            &!
      ind, ind1,      &!
      l,              &!
      ibnd             !
  LOGICAL :: &
      so(npsx),       &!
      minus_q,        &!
      ltest            !
  REAL(DP) :: &
      vionl,       &   !
      iocc             !
  INTEGER, EXTERNAL :: &
      n_atom_wfc,     &!
      set_Hubbard_l
  LOGICAL, EXTERNAL :: &
      lchk_tauxk       ! tests that atomic coordinates do not overlap
  !
  !
  ALLOCATE( m_loc( 3, nat ) )
  !
  IF ( nimage > 1 .AND. .NOT. lpath ) &
     CALL errore( 'setup', 'images parallelization not permitted', 1 )
  !
  DO nt = 1, ntyp
     !
     DO ir = 1, mesh(nt)
        !
        IF ( r(ir,nt) > rcut ) THEN
           !
           msh(nt) = ir
           !
           GO TO 5
           !
        END IF
        !
     END DO
     !
     msh(nt) = mesh(nt)
     !
     ! ... force msh to be odd for simpson integration
     !
5    msh(nt) = 2 * ( ( msh(nt) + 1 ) / 2 ) - 1
     !
  END DO
  !
  ! ... Compute the ionic charge for each atom type
  !
  zv(1:ntyp) = zp(1:ntyp)
  !
  ! ... Set the number of electrons equal to the total ionic charge
  !
  IF ( nelec == 0.D0 ) THEN
     !
#if defined (__PGI)
     DO na = 1, nat
        nelec = nelec + zv( ityp(na) )
     END DO
#else
     nelec = SUM( zv(ityp(1:nat)) )
#endif
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
     DO nt = 1, ntyp
        !
        angle1(nt) = pi * angle1(nt) / 180.D0
        angle2(nt) = pi * angle2(nt) / 180.D0
        !
     END DO
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
      CALL infomsg( 'setup', 'the system is metallic, specify occupations', -1 )
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
  ! ... Set the number of occupied bands if not given in input
  !
  IF ( nbnd == 0 ) THEN
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
  IF ( lphonon .or. lraman ) THEN
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
  IF ( .NOT. lscf ) niter = 1
  !
  starting_diag_threshold = ethr
  !
  !   check if spin-orbit is possible
  !
  lspinorb = lspinorb .AND. noncolin
  !
  ! ... if this is not a spin-orbit calculation, all spin-orbit pseudopotentials
  ! ... are transformed into standard pseudopotentials
  !
  DO nt = 1, ntyp
     !
     so(nt) = .TRUE.
     !
     DO nb = 1, nbeta(nt)
        !
        so(nt) = so(nt) .AND. ( ABS( jjj(nb,nt) ) > 1.D-7 )
        !
     END DO
     !
  END DO
  !
  IF ( .NOT. lspinorb ) THEN
     !
     DO nt = 1, ntyp
        !
        IF ( so(nt) ) THEN
           !
           IF ( tvanp(nt) ) &
                CALL errore( 'setup', 'US j-average not yet implemented', 1 )
           !
           nbe = 0
           !
           DO nb = 1, nbeta(nt)
              !
              nbe = nbe + 1
              !
              IF ( lll(nb,nt) /= 0 .AND. &
                   ABS( jjj(nb,nt) - lll(nb,nt) - 0.5D0 ) < 1.D-7 ) nbe = nbe - 1
           END DO
           !
           nbeta(nt) = nbe
           !
           nbe = 0
           !
           DO nb = 1, nbeta(nt)
              !
              nbe = nbe + 1
              !
              l = lll(nbe,nt)
              !
              IF ( l /= 0 ) THEN
                 !
                 IF (ABS(jjj(nbe,nt)-lll(nbe,nt)+0.5d0).LT.1.d-7) THEN
                    IF (ABS(jjj(nbe+1,nt)-lll(nbe+1,nt)-0.5d0).GT.1.d-7) &
                       call errore('setup','wrong beta functions',1)
                    ind=nbe+1
                    ind1=nbe
                 ELSE
                    IF (ABS(jjj(nbe+1,nt)-lll(nbe+1,nt)+0.5d0).GT.1.d-7) &
                       call errore('setup','wrong beta functions',1)
                    ind=nbe
                    ind1=nbe+1
                 ENDIF
                 !
                 vionl = ( ( l + 1.D0 ) * dion(ind,ind,nt) + &
                      l * dion(ind1,ind1,nt) ) / ( 2.D0 * l + 1.D0 )
                 !
                 betar(1:mesh(nt),nb,nt) = 1.D0 / ( 2.D0 * l + 1.D0 ) * &
                      ( ( l + 1.D0 ) * SQRT( dion(ind,ind,nt) / vionl ) * &
                      betar(1:mesh(nt),ind,nt) + &
                      l * SQRT( dion(ind1,ind1,nt) / vionl ) * &
                      betar(1:mesh(nt),ind1,nt) )
                 !
                 dion(nb,nb,nt) = vionl
                 !
                 nbe = nbe + 1
                 !
              ELSE
                 !
                 betar(1:mesh(nt),nb,nt) = betar(1:mesh(nt),nbe,nt)
                 !
                 dion(nb,nb,nt) = dion(nbe,nbe,nt)
                 !
              END IF
              !
              lll(nb,nt)=lll(nbe,nt)
              !
           END DO
           !
           nbe = 0
           !
           DO nb = 1, nchi(nt)
              !
              nbe = nbe + 1
              !
              IF ( lchi(nb,nt) /= 0 .AND. &
                   ABS(jchi(nb,nt)-lchi(nb,nt)-0.5D0 ) < 1.D-7 ) nbe = nbe - 1
              !
           END DO
           !
           nchi(nt) = nbe
           ! 
           nbe = 0
           !
           do nb = 1, nchi(nt)
              !
              nbe = nbe + 1
              !
              l = lchi(nbe,nt)
              !
              IF ( l /= 0 ) THEN
                 !
                 IF (ABS(jchi(nbe,nt)-lchi(nbe,nt)+0.5d0).LT.1.d-7) THEN
                    IF (ABS(jchi(nbe+1,nt)-lchi(nbe+1,nt)-0.5d0).GT.1.d-7) &
                       call errore('setup','wrong chi functions',1)
                    ind=nbe+1
                    ind1=nbe
                 ELSE
                    IF (ABS(jchi(nbe+1,nt)-lchi(nbe+1,nt)+0.5d0).GT.1.d-7) &
                       call errore('setup','wrong chi functions',1)
                    ind=nbe
                    ind1=nbe+1
                 END IF
                 !
                 chi(1:mesh(nt),nb,nt)=((l+1.D0) * chi(1:mesh(nt),ind,nt)+ &
                      l * chi(1:mesh(nt),ind1,nt)) / ( 2.D0 * l + 1.D0 )

                 nbe = nbe + 1
                 !
              ELSE
                 !
                 chi(1:mesh(nt),nb,nt) = chi(1:mesh(nt),nbe,nt)
                 !
              END IF
              !
              lchi(nb,nt)= lchi(nbe,nt)
              !
           END DO
           !
        END IF
        !
     END DO
     !
  END IF
  !
  ! ... set number of atomic wavefunctions
  !
  natomwfc = n_atom_wfc( nat, npsx, ityp, nchix, nchi, oc, lchi, jchi )
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
  IF ( lxkcry ) CALL cryst_to_cart( nks, xk, bg, 1 )
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
  IF ( nks < 0 ) THEN
     !
     CALL setupkpoint( s, nrot, xk, wk, nks, npk, nk1, &
                       nk2, nk3, k1, k2, k3, at, bg, tipo )
     !
  ELSE IF ( nks == 0 ) THEN
     !
     IF ( lberry ) THEN
        !
        CALL kp_strings( nppstr, gdir, nrot, s, bg, npk, &
                         k1, k2, k3, nk1, nk2, nk3, nks, xk, wk )
        !
        nosym = .TRUE.
        nrot  = 1
        nsym  = 1
        !
     ELSE
        !
        CALL kpoint_grid( nrot, s, bg, npk, k1, k2, k3, &
                          nk1, nk2, nk3, nks, xk, wk )
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
  ! ... If some symmetries are missing in the crystal, "sgama" computes
  ! ... the missing k-points. If nosym is true (see above) we do not use
  ! ... any point-group symmetry and leave k-points unchanged.
  !
  input_nks = nks
  !
  CALL sgama( nrot, nat, s, sname, t_rev, at, bg, tau, ityp, nsym, nr1,&
              nr2, nr3, irt, ftau, npk, nks, xk, wk, invsym, minus_q,  &
              xqq, modenum, noncolin, m_loc )
  !
  CALL checkallsym( nsym, s, nat, tau, ityp, at, &
                    bg, nr1, nr2, nr3, irt, ftau )
  !
  ! ... if dynamics is done the system should have no symmetries
  ! ... (inversion symmetry alone is allowed)
  !
  IF ( lmd .AND. ( nsym == 2 .AND. .NOT. invsym .OR. nsym > 2 ) &
           .AND. .NOT. ( calc == 'mm' .OR. calc == 'nm' ) ) &
     CALL infomsg( 'setup', 'Dynamics, you should have no symmetries', -1 )
  !
  ! ... Calculate quantities used in tetrahedra method
  !
  IF ( ltetra ) THEN
     !
     ntetra = 6 * nk1 * nk2 * nk3
     !
     ALLOCATE( tetra( 4, ntetra ) )    
     !
     CALL tetrahedra( nsym, s, minus_q, at, bg, npk, k1, k2, k3, &
                      nk1, nk2, nk3, nks, xk, wk, ntetra, tetra )
     !
  ELSE
     !
     ntetra = 0
     !
  END IF
  !
  ! ... non scf calculation: do not change the number of k-points
  ! ... to account for reduced symmetry, unless you need to
  ! ... (as in phonon or raman or DOS calculations, or whenever the
  ! ...  Fermi energy has to be calculated)
  !
  ltest = ( nks /= input_nks ) .AND. & 
          ( .NOT. ( ltetra .OR. lgauss ) ) .AND. &
          ( .NOT. lscf ) .AND. ( .NOT. ( lphonon .OR. lraman ) ) 
  IF ( ltest ) THEN
     !
     WRITE( stdout, '(/,5X,"Only input k-points are used ", &
                         & "(inequivalent points not generated)",/)' )
     !
     nks = input_nks
     !
  END IF
  !
  ! ... phonon calculation: add k+q to the list of k
  !
  IF ( lphonon ) CALL set_kplusq( xk, wk, xqq, nks, npk ) 
  !
  ! ... raman calculation: add k+b to the list of k
  !
  IF ( lraman ) CALL set_kplusb(ibrav, xk, wk, b_length, nks, npk, lcart)
  !
#if defined (EXX)
  IF ( dft_is_hybrid() ) call exx_grid_init()
#endif
  !
  IF ( lsda ) THEN
     !
     ! ... LSDA case: two different spin polarizations, 
     ! ...            each with its own kpoints
     !
     if (nspin /= 2) call errore ('setup','nspin should be 2; check iosys',1)
     !
     CALL set_kup_and_kdw( xk, wk, isk, nks, npk )
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
     wk(1:nks)    = wk(1:nks) * degspin
     current_spin = 1
     if (nspin /= 1) call errore ('setup','nspin should be 1; check iosys',1)
     !
  END IF
  !
  IF ( nks > npk ) CALL errore( 'setup', 'too many k points', nks )
  !
#ifdef __PARA
  !
  ! ... set the granularity for k-point distribution
  !
  IF ( ( ABS( xqq(1) ) < eps .AND. ABS( xqq(2) ) < eps .AND. &
         ABS( xqq(3) ) < eps) .OR. ( .NOT. lphonon ) ) THEN
     !
     kunit = 1
     !
  ELSE
     !
     kunit = 2
     !
  ENDIF
  !
  IF ( lraman ) THEN
     !
     IF( lcart ) THEN
        !
        kunit = 7
        !
     ELSE
        !
        IF ( ibrav == 1 ) kunit =  7
        IF ( ibrav == 2 ) kunit =  9
        IF ( ibrav == 3 ) kunit = 13
        !
     END IF
     !
  END IF
  !
  ! ... distribute the k-points (and their weights and spin indices)
  !
  CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
  !
#else
  !
  ! ... set nkstot which is used to write results for all k-points
  !
  nkstot = nks
  !
#endif
  !
  ! ... okvan = .TRUE. : at least one pseudopotential is US
  !
  okvan = ANY( tvanp(:) )
  !
  ! ... initialize parameters for charge density extrapolation during dynamics
  !
  alpha0 = 1.D0
  beta0  = 0.D0
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
!
!----------------------------------------------------------------------------
FUNCTION n_atom_wfc( nat, npsx, ityp, nchix, nchi, oc, lchi, jchi )
  !----------------------------------------------------------------------------
  !
  ! ... Find max number of bands needed
  !
  USE kinds,            ONLY : DP
  use noncollin_module, ONLY : noncolin
  use spin_orb,         ONLY : lspinorb
  !
  IMPLICIT NONE
  !
  INTEGER       :: n_atom_wfc
  INTEGER       :: nat, npsx, ityp(nat), nchix, nchi(npsx), lchi(nchix,npsx)
  REAL(DP) :: oc(nchix,npsx), jchi(nchix,npsx)
  INTEGER       :: na, nt, n
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
              IF ( lspinorb ) THEN
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
           !
        END IF
        !
     END DO
     !
  END DO
  !
  RETURN
  !
END FUNCTION n_atom_wfc
