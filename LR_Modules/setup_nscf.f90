!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE setup_nscf ( newgrid, xq, elph_mat )
  !----------------------------------------------------------------------------
  !
  ! ... This routine initializes variables for the non-scf calculations at k
  ! ... and k+q required by the linear response calculation at finite q.
  ! ... In particular: finds the symmetry group of the crystal that leaves
  ! ... the phonon q-vector (xq) or the single atomic displacement (modenum)
  ! ... unchanged; determines the k- and k+q points in the irreducible BZ
  ! ... Needed on input (read from data file):
  ! ... "nsym" crystal symmetries s, ftau, t_rev, "nrot" lattice symetries "s"
  ! ... "nkstot" k-points in the irreducible BZ wrt lattice symmetry
  ! ... Produced on output:
  ! ... symmetries ordered with the "nsymq" phonon symmetries first
  ! ... "nkstot" k- and k+q-points in the IBZ calculated for the phonon sym.)
  ! ... Misc. data needed for running the non-scf calculation
  !
  USE kinds,              ONLY : DP
  USE parameters,         ONLY : npk
  USE io_global,          ONLY : stdout
  USE constants,          ONLY : pi, degspin
  USE cell_base,          ONLY : at, bg
  USE ions_base,          ONLY : nat, tau, ityp, zv
  USE force_mod,          ONLY : force
  USE basis,              ONLY : natomwfc
  USE klist,              ONLY : xk, wk, nks, nelec, degauss, lgauss, &
                                 nkstot, qnorm
  USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk
  USE symm_base,          ONLY : s, t_rev, irt, ftau, nrot, nsym, &
                                 time_reversal
  USE wvfct,              ONLY : nbnd, nbndx
  USE control_flags,      ONLY : ethr, isolve, david, max_cg_iter, &
                                 noinv, use_para_diag
!!!!!!!!!  USE el_phon,            ONLY : elph_mat
  USE mp_pools,           ONLY : kunit
  USE spin_orb,           ONLY : domag
  USE noncollin_module,   ONLY : noncolin
  USE start_k,            ONLY : nks_start, xk_start, wk_start, &
                                 nk1, nk2, nk3, k1, k2, k3
  USE paw_variables,      ONLY : okpaw
  USE uspp_param,         ONLY : n_atom_wfc
 
  USE lr_symm_base, ONLY : nsymq, invsymq, minus_q
  USE control_lr,   ONLY : lgamma
  !
  IMPLICIT NONE
  !
  REAL (DP), INTENT(IN) :: xq(3)
  LOGICAL, INTENT (IN) :: newgrid
  LOGICAL, INTENT (IN) :: elph_mat  ! used to be passed through a module. 
  !
  REAL (DP), ALLOCATABLE :: rtau (:,:,:)
  LOGICAL  :: magnetic_sym, sym(48)
  LOGICAL  :: skip_equivalence
  !
  IF ( .NOT. ALLOCATED( force ) ) ALLOCATE( force( 3, nat ) )
  !
  ! ... threshold for diagonalization ethr - should be good for all cases
  !
  ethr= 1.0D-9 / nelec
  !
  ! ... variables for iterative diagonalization
  ! ... Davdson: isolve=0, david=4 ; CG: isolve=1, david=1
  isolve = 0
  david  = 4
  nbndx = david*nbnd
  max_cg_iter=20
  natomwfc = n_atom_wfc( nat, ityp, noncolin )
  !
#if defined(__MPI)
  IF ( use_para_diag )  CALL check_para_diag( nbnd )
#else
  use_para_diag = .FALSE.
#endif
  !
  ! ... Symmetry and k-point section
  !
  ! ... time_reversal = use q=>-q symmetry for k-point generation
  !
  magnetic_sym = noncolin .AND. domag
  !
  ! ... smallg_q flags in symmetry operations of the crystal
  ! ... that are not symmetry operations of the small group of q
  !
  CALL set_small_group_of_q(nsymq,invsymq,minus_q)
  !
  ! ... Input k-points are assumed to be  given in the IBZ of the Bravais
  ! ... lattice, with the full point symmetry of the lattice.
  !
  if( nks_start > 0 .AND. .NOT. newgrid ) then
     !
     !  In this case I keep the same points of the Charge density
     !  calculations
     !
     nkstot = nks_start
     xk(:,1:nkstot) = xk_start(:,1:nkstot)
     wk(1:nkstot)   = wk_start(1:nkstot)
  else
     !
     ! In this case I generate a new set of k-points
     !
     ! In the case of electron-phonon matrix element with wannier functions 
     ! (and possibly in other cases as well) the k-points should not be reduced
     !
     skip_equivalence = elph_mat
     CALL kpoint_grid ( nrot, time_reversal, skip_equivalence, s, t_rev, &
                      bg, nk1*nk2*nk3, k1,k2,k3, nk1,nk2,nk3, nkstot, xk, wk)
  endif

  !
  ! ... If some symmetries of the lattice are missing in the crystal,
  ! ... "irreducible_BZ" computes the missing k-points.
  !
  if(.not.elph_mat) &
  CALL irreducible_BZ (nrot, s, nsymq, minus_q, magnetic_sym, &
                       at, bg, npk, nkstot, xk, wk, t_rev)
  !
  ! ... add k+q to the list of k
  !
  CALL set_kplusq( xk, wk, xq, nkstot, npk )
  !
  IF ( lsda ) THEN
     !
     ! ... LSDA case: two different spin polarizations,
     ! ...            each with its own kpoints
     !
     if (nspin /= 2) call errore ('setup_nscf','nspin should be 2; check iosys',1)
     !
     CALL set_kup_and_kdw( xk, wk, isk, nkstot, npk )
     !
  ELSE IF ( noncolin ) THEN
     !
     ! ... noncolinear magnetism: potential and charge have dimension 4 (1+3)
     !
     if (nspin /= 4) call errore ('setup_nscf','nspin should be 4; check iosys',1)
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
        CALL errore( 'setup_nscf', 'nspin should be 1; check iosys', 1 )
     !
  END IF
  !
  IF ( nkstot > npk ) CALL errore( 'setup_nscf', 'too many k points', nkstot )
  !
  ! ...notice: qnorm is used by allocate_nlpot to determine
  ! the correct size of the interpolation table "qrad"
  !
  qnorm = sqrt(xq(1)**2 + xq(2)**2 + xq(3)**2)
  !
  ! ... set the granularity for k-point distribution
  !
  IF ( lgamma  ) THEN
     kunit = 1
  ELSE
     kunit = 2
  ENDIF
  !
  ! ... distribute k-points (and their weights and spin indices)
  !
  CALL divide_et_impera( nkstot, xk, wk, isk, nks )
  !
  RETURN
  !
END SUBROUTINE setup_nscf
