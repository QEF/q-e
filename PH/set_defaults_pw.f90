!
! Copyright (C) 2008 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE set_defaults_pw
  !-----------------------------------------------------------------------------
  !
  ! ...  this subroutine sets the default values for the variables
  ! ...  read from input by pw which are not saved into the xml file. 
  ! ...  It has to be called by programs that run "electrons" without 
  ! ...  reading the input data. It is possible to change the default
  ! ...  values initialized here. Variables in input_parameters are
  ! ...  initialized by that module  and do not need to be initialized here. 
  ! ...  Actually this routine should not be needed. All the variables
  ! ...  which are read from input should have a default value and should
  ! ...  be in a module that initializes them. The pw code should then
  ! ...  use those variables. Unfortunately the routine input now 
  ! ...  initializes several variables that are in pw modules and
  ! ...  have just a different name from the variable contained in 
  ! ...  input parameters, or are calculated starting from the variable
  ! ...  contained in input_parameters.  
  ! ...  Moreover many variables contained in control_flags are not 
  ! ...  initialized and need to be initialized here ... 
  !
  !
  USE kinds,         ONLY : DP
  USE bp,            ONLY : lberry,   &
                            lelfield
  !
  USE basis,         ONLY : startingwfc, &
                            startingpot
  !
  USE char,          ONLY : crystal
  !
  USE cellmd,        ONLY : calc, lmovecell
  !
  USE force_mod,     ONLY : lforce, lstres
  !
  USE gvect,         ONLY : ecfixed, qcutz, q2sigma
  !
  USE klist,         ONLY : lxkcry, tot_charge, &
                            tot_magnetization, &
                            multiplicity

  USE relax,         ONLY : starting_scf_threshold
  !
  USE control_flags, ONLY : isolve, max_cg_iter, tr2, imix, &
                            nmix, iverbosity, niter, pot_order, wfc_order, &
                            assume_isolated, &
                            diago_full_acc, &
                            mixing_beta, &
                            upscale, &
                            nstep, &
                            iprint, &
                            nosym, &
                            io_level, lscf, lbfgs, lmd, lpath, lneb,   &
                            lsmd, lphonon, ldamped, lbands, lmetadyn, llang, &
                            lconstrain, lcoarsegrained, restart, &
                            use_para_diag

  USE bfgs_module,   ONLY : bfgs_ndim, &
                            trust_radius_max, &
                            trust_radius_min, &
                            trust_radius_ini, &
                            w_1, &
                            w_2
  USE us, ONLY : spline_ps
  USE a2F, ONLY : la2F
  USE klist, ONLY : nkstot, wk, xk
  USE start_k,         ONLY : xk_start, wk_start, nks_start

  !
  IMPLICIT NONE
  !
  iprint = 100000
  lberry   = .FALSE.
  lelfield = .FALSE.
  lxkcry=.FALSE.
  tot_charge = 0.0_DP
  tot_magnetization = -1
  multiplicity = 0
  nosym = .FALSE.
  ecfixed = 0.0_DP
  qcutz   = 0.0_DP
  q2sigma = 0.01_DP
  !
  !  ... postprocessing of DOS & phonons & el-ph
  la2F = .FALSE.
  !
  ! ... non collinear program variables
  !
  assume_isolated = .FALSE.
  !
  spline_ps = .FALSE.
  !
  diago_full_acc = .FALSE.
  !
  upscale           = 10.0_DP
  mixing_beta       = 0.7
  !
  ! ... BFGS defaults
  !
  bfgs_ndim        = 1
  trust_radius_max = 0.8_DP   ! bohr
  trust_radius_min = 1.E-4_DP ! bohr
  trust_radius_ini = 0.5_DP   ! bohr
  w_1              = 0.01_DP
  w_2              = 0.50_DP
  !
  startingpot = 'file'
  startingwfc = 'atomic'
  !
  restart        = .FALSE.
  !
  io_level = 1
  !
  ! ... various initializations of control variables
  !
  lscf      = .FALSE.
  lmd       = .FALSE.
  lmetadyn  = .FALSE.
  lpath     = .FALSE.
  lneb      = .FALSE.
  lsmd      = .FALSE.
  lmovecell = .FALSE.
  lphonon   = .FALSE.
  lbands    = .FALSE.
  lbfgs     = .FALSE.
  ldamped   = .FALSE.
  lforce    = .FALSE.
  lstres    = .FALSE.
  !
  nstep = 1
  !
  isolve = 1
  max_cg_iter = 100
  use_para_diag = .FALSE.
  !
  niter = 1000
  !
  pot_order = 0
  wfc_order = 0
  !
  tr2=1.D-6
  starting_scf_threshold = tr2
  imix = 0
  nmix = 0
  !
  iverbosity = 0
  !
  crystal     = ' '
  calc      = ' '
  !
  !  Reset the k points
  !
  nkstot = nks_start
  !
  xk(:,1:nkstot) = xk_start(:,1:nkstot)
  wk(1:nkstot)   = wk_start(1:nkstot)
  !
  RETURN
  !
END SUBROUTINE set_defaults_pw
!
! Copyright (C) 2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE setup_ph()
  !----------------------------------------------------------------------------
  !
  ! ... This routine finds the symmetry group of the crystal that leaves the
  ! ... phonon q-vector (xqq) or the single atomic displacement (modenum)
  ! ... unchanged; determines the k- and k+q points in the irreducible BZ
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : eps8
  USE parameters,         ONLY : npk
  USE io_global,          ONLY : stdout
  USE constants,          ONLY : pi, degspin
  USE cell_base,          ONLY : at, bg, alat, tpiba, tpiba2, ibrav, &
                                 symm_type, omega
  USE ions_base,          ONLY : nat, tau, ntyp => nsp, ityp, zv
  USE basis,              ONLY : natomwfc
  USE gvect,              ONLY : nr1, nr2, nr3
  USE klist,              ONLY : xk, wk, xqq, nks, nelec, degauss, lgauss, &
                                 nkstot
  USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk
  USE symme,              ONLY : s, t_rev, irt, ftau, nsym, invsym, time_reversal
  USE wvfct,              ONLY : nbnd, nbndx
  USE control_flags,      ONLY : tr2, ethr, lphonon, isolve, david, &
                                 noinv, nosym, modenum, use_para_diag
  USE char,               ONLY : sname
  USE mp_global,          ONLY : kunit
  USE spin_orb,           ONLY : domag
  USE noncollin_module,   ONLY : noncolin, m_loc
  USE start_k,            ONLY : nks_start, xk_start, wk_start
  !
  IMPLICIT NONE
  !
  INTEGER  :: na, nt, nrot, irot, isym, tipo, is, nb, ierr, ik
  LOGICAL  :: minus_q, magnetic_sym
  !
  INTEGER, EXTERNAL :: n_atom_wfc
  !
  ! time reversal operation is set up to 0 by default
  t_rev = 0
  !
  IF ( lphonon ) THEN
     !
     ! ... in the case of a phonon calculation ethr can not be specified
     ! ... in the input file
     !
     IF ( ethr /= 0.D0 ) &
        WRITE( UNIT = stdout, &
             & FMT = '(5X,"diago_thr_init overwritten ", &
             &            "with conv_thr / nelec")' )
     !
     ethr = 0.1D0 * MIN( 1.D-2, tr2 / nelec )
     !
  END IF
  !
  ! ... set number of atomic wavefunctions
  !
  natomwfc = n_atom_wfc( nat, ityp )
  !
  ! ... set the max number of bands used in iterative diagonalization
  !
  nbndx = nbnd
  IF ( isolve == 0 ) nbndx = david * nbnd
  !
#ifdef __PARA
  IF ( use_para_diag )  CALL check_para_diag( nelec )
#else
  use_para_diag = .FALSE.
#endif
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
  ! ... If nosym is true do not use any point-group symmetry
  !
  IF ( nosym ) nrot = 1
  !
  ! ... time_reversal = use q=>-q symmetry for k-point generation
  !
  magnetic_sym = noncolin .AND. domag 
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  !
  ! ...  allocate space for irt
  !
  ALLOCATE( irt( 48, nat ) )
  !
  ! ... "sgama" eliminates rotations that are not symmetry operations
  !
  CALL sgama( nrot, nat, s, sname, t_rev, at, bg, tau, ityp, nsym, &
              nr1, nr2, nr3, irt, ftau, invsym, minus_q, xqq, &
              modenum, time_reversal, magnetic_sym, m_loc)
  !
  CALL checkallsym( nsym, s, nat, tau, ityp, at, &
          bg, nr1, nr2, nr3, irt, ftau, alat, omega )
  !
  ! ... Input k-points are assumed to be  given in the IBZ of the Bravais
  ! ... lattice, with the full point symmetry of the lattice.
  ! ... If some symmetries of the lattice are missing in the crystal,
  ! ... "irreducible_BZ" computes the missing k-points.
  !
  CALL irreducible_BZ (nrot, s, nsym, at, bg, npk, nkstot, xk, wk, minus_q)
  !
  ! ... phonon calculation: add k+q to the list of k
  !
  IF ( lphonon ) CALL set_kplusq( xk, wk, xqq, nkstot, npk )
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
#else
  !
  nks = nkstot
  !
#endif
  !
  RETURN
  !
END SUBROUTINE setup_ph
