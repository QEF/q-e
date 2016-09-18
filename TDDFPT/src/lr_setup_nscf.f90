!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE lr_setup_nscf ()
  !---------------------------------------------------------------------
  !
  ! This subroutine initializes variables for the non-scf calculations at k
  ! and k+q required by the linear response calculation at finite q.
  ! Determines k and k+q points in the irreducible BZ.
  ! Needed on input (read from data file):
  ! "nkstot" k-points in the irreducible BZ wrt lattice symmetry.
  ! Inspired by PH/set_defaults_pw.f90
  !
  ! Created by Iurii Timrov (2013)
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : eps8
  USE parameters,         ONLY : npk
  USE constants,          ONLY : pi, degspin
  USE cell_base,          ONLY : at, bg, alat, tpiba, tpiba2, ibrav, omega
  USE ions_base,          ONLY : nat, tau, ityp, zv
  USE force_mod,          ONLY : force
  USE basis,              ONLY : natomwfc
  USE klist,              ONLY : xk, wk, nks, nelec, degauss, lgauss, &
                                 nkstot, qnorm
  USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk
  USE symm_base,          ONLY : s, t_rev, nrot, nsym, time_reversal
  USE wvfct,              ONLY : nbnd, nbndx
  USE control_flags,      ONLY : ethr, isolve, david, use_para_diag, &
                                 & noinv, max_cg_iter
  USE mp_pools,           ONLY : kunit
  USE spin_orb,           ONLY : domag
  USE noncollin_module,   ONLY : noncolin
  USE start_k,            ONLY : nks_start, xk_start, wk_start, &
                                 nk1, nk2, nk3, k1, k2, k3
  USE uspp_param,         ONLY : n_atom_wfc
  USE lr_symm_base,       ONLY : nsymq, minus_q
  USE qpoint,             ONLY : xq
  ! 
  IMPLICIT NONE
  !
  LOGICAL :: magnetic_sym 
  !
  CALL start_clock( 'lr_setup_nscf' )
  ! 
  IF ( .NOT. ALLOCATED( force ) ) ALLOCATE( force( 3, nat ) )
  !
  ! ... threshold for diagonalization ethr - should be good for all cases
  !
  ethr = 1.0D-9 / nelec
  !
  ! ... variables for iterative diagonalization (Davidson is assumed)
  !
  isolve = 0
  david  = 4
  nbndx  = david*nbnd
  max_cg_iter = 20
  natomwfc = n_atom_wfc( nat, ityp, noncolin )
  !
#if defined(__MPI)
  IF ( use_para_diag )  CALL check_para_diag( nbnd )
#else
  use_para_diag = .FALSE.
#endif
  !
  ! Symmetry section
  !
  ! time_reversal = use q=>-q symmetry for k-point generation
  !
  magnetic_sym = noncolin .AND. domag
  !
  !time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  !
  ! Determine the small group of q : S q = q + G
  !
  CALL lr_smallgq (xq)
  !
  !  K points section   
  !
  ! Input k-points are assumed to be given in the IBZ of the Bravais
  ! lattice, with the full point symmetry of the lattice.
  !
  IF ( nks_start > 0 ) THEN
     !
     ! In this case keep the same points of the charge-density calculation
     !
     nkstot = nks_start
     xk(:,1:nkstot) = xk_start(:,1:nkstot)
     wk(1:nkstot)   = wk_start(1:nkstot)
     !
  ELSE
     !
     ! Generate a new set of k-points
     !
     CALL kpoint_grid ( nrot, time_reversal,.false., s, t_rev, &
                      bg, nk1*nk2*nk3, k1,k2,k3, nk1,nk2,nk3, nkstot, xk, wk)
     !
  ENDIF
  !
  ! If some symmetries of the lattice are missing in the crystal,
  ! "irreducible_BZ" computes the missing k-points.
  !
  CALL irreducible_BZ (nrot, s, nsymq, minus_q, magnetic_sym, &
                       at, bg, npk, nkstot, xk, wk, t_rev)
  !
  ! Add k+q to the list of k
  !
  CALL set_kplusq( xk, wk, xq, nkstot, npk )
  !
  IF ( lsda ) THEN
     !
     ! ... LSDA case: two different spin polarizations,
     ! ...            each with its own kpoints
     !
     IF (nspin /= 2) CALL errore ('lr_setup_nscf','nspin should be 2; check iosys',1)
     !
     CALL set_kup_and_kdw( xk, wk, isk, nkstot, npk )
     !
  ELSEIF ( noncolin ) THEN
     !
     ! ... noncolinear magnetism: potential and charge have dimension 4 (1+3)
     !
     IF (nspin /= 4) CALL errore ('lr_setup_nscf','nspin should be 4; check iosys',1)
     current_spin = 1
     !
  ELSE
     !
     ! ... LDA case: the two spin polarizations are identical
     !
     wk(1:nkstot) = wk(1:nkstot) * degspin
     current_spin = 1
     !
     IF ( nspin /= 1 ) &
        CALL errore( 'lr_setup_nscf', 'nspin should be 1; check iosys', 1 )
     !
  ENDIF
  !
  IF ( nkstot > npk ) CALL errore( 'lr_setup_nscf', 'too many k points', nkstot )
  !
  ! ...notice: qnorm is used by allocate_nlpot to determine
  ! the correct size of the interpolation table "qrad"
  !
  qnorm = sqrt(xq(1)**2 + xq(2)**2 + xq(3)**2)
  !
  ! ... set the granularity for k-point distribution
  !
  IF ( ABS( xq(1) ) < eps8 .AND. ABS( xq(2) ) < eps8 .AND. &
       ABS( xq(3) ) < eps8 ) THEN
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
  CALL divide_et_impera( nkstot, xk, wk, isk, nks )
  !
  CALL stop_clock( 'lr_setup_nscf' )
  !
  RETURN
  !
END SUBROUTINE lr_setup_nscf
