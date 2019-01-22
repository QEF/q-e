!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_setup_q()
  !-----------------------------------------------------------------------
  !
  !  This subroutine prepares several variables which are needed in the
  !  HP program at fixed q point:
  !  1) computes the total local potential (external+scf) on the smooth
  !     grid to be used in h_psi (and other places?)
  !  2) Set the nonlinear core correction variable 
  !  3) Allocate the variable of the local magnetization
  !  4) Compute the derivative of the XC potential (dmuxc)
  !  5) Setup gradient correction stuff
  !  6) Compute the inverse of each matrix of the crystal symmetry group
  !  7) Computes the number of occupied bands for each k point
  !  8) Compute alpha_pv
  !  9) Set various symmetry-related variables
  !       time_reversal true if there is time-reversal symmetry
  !       gi            the G associated to each symmetry operation
  !       gimq          the G of the q -> -q+G symmetry
  !       nsymq         the order of the small group of q
  !       irotmq        the index of the q->-q+G symmetry
  !       minus_q       true if there is a symmetry sending q -> -q+G
  !       rtau          rtau = S\tau_a - \tau_b
  ! 10) Setup the parameters alpha_mix needed for the 
  !     solution of the linear system  
  ! 11) Initialize d1, d2, d3 to rotate the spherical harmonics
  !
  !  IMPORTANT NOTE ABOUT SYMMETRIES:
  !  nrot  is the number of sym.ops. of the Bravais lattice
  !        read from data file, only used in set_default_pw
  !  nsym  is the number of sym.ops. of the crystal symmetry group
  !        read from data file, should never be changed
  !  nsymq is the number of sym.ops. of the small group of q
  !        it is calculated in set_defaults_pw for each q
  !  The matrices "s" of sym.ops are ordered as follows:
  !   first the nsymq sym.ops. of the small group of q
  !   (the ordering is done in subroutine copy_sym in set_defaults_pw),
  !   followed by the remaining nsym-nsymq sym.ops. of the crystal group,
  !   followed by the remaining nrot-nsym sym.ops. of the Bravais  group
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : tau, nat, ntyp => nsp, ityp
  USE cell_base,        ONLY : at, bg
  USE io_global,        ONLY : stdout
  USE lsda_mod,         ONLY : nspin
  USE scf,              ONLY : v, vrs, vltot, rho, kedtau
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm
  USE gvecs,            ONLY : doublegrid
  USE symm_base,        ONLY : nrot, nsym, s, ftau, irt, time_reversal, &
                               inverse_s, d1, d2, d3
  USE uspp_param,       ONLY : upf
  USE uspp,             ONLY : nlcc_any
  USE spin_orb,         ONLY : domag
  USE constants,        ONLY : degspin, pi, rytoev
  USE noncollin_module, ONLY : noncolin, m_loc, nspin_mag
  USE wvfct,            ONLY : nbnd, et
  USE control_flags,    ONLY : noinv
  USE eqv,              ONLY : dmuxc
  USE qpoint,           ONLY : xq
  USE control_lr,       ONLY : lgamma
  USE lr_symm_base,     ONLY : gi, gimq, irotmq, minus_q, invsymq, nsymq, rtau
  USE ldaU_hp,          ONLY : niter_max, search_sym, alpha_mix, skip_equivalence_q
  !
  IMPLICIT NONE
  INTEGER :: ir, isym, ik, it
  LOGICAL :: sym(48), magnetic_sym, is_symmorphic
  !
  CALL start_clock ('hp_setup_q')
  !
  ! 1) Compute the total local potential (external+scf) on the smooth grid
  !
  CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
  !
  ! 2) Set the nonlinear core correction variable
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  !
  ! 3) Allocate the variable of the local magnetization. 
  !    This is needed in find_sym
  !
  IF (.NOT.ALLOCATED(m_loc)) ALLOCATE( m_loc( 3, nat ) )
  ! 
  ! 4) Compute the derivative of the XC potential (dmuxc)
  !
  CALL setup_dmuxc()
  !
  ! 5) Setup gradient correction stuff
  !
  CALL setup_dgc()
  !
  ! 6) Compute the inverse of each matrix of the crystal symmetry group
  !
  CALL inverse_s()
  !
  ! 7) Computes the number of occupied bands for each k point
  !
  CALL setup_nbnd_occ()
  !
  ! 8) Compute alpha_pv
  !
  CALL setup_alpha_pv()
  !
  ! 9) Set various symmetry-related variables
  !
  magnetic_sym = noncolin .AND. domag
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  !
  ! The small group of q was already determined. At q\=0 it is calculated
  ! by set_nscf, at q=0 it coincides with the point group and we take nsymq=nsym
  !
  IF (lgamma) THEN
     !
     nsymq   = nsym
     !
     IF ( time_reversal ) THEN
        minus_q = .TRUE.
     ELSE
        minus_q = .FALSE.
     ENDIF
     !
  ENDIF
  !
  ! Calculate rtau (the Bravais lattice vector associated to a rotation) 
  ! with the new symmetry order
  !
  CALL sgam_lr (at, bg, nsym, s, irt, tau, rtau, nat)
  !
  ! Calculate the vectors G associated to the symmetry Sq = q + G
  ! If minus_q=.true. calculate also irotmq and the G associated to Sq=-q+G
  !
  CALL set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)
  !
  ! Check if there are fractional translations
  ! Note: Try to use PH/symmorphic_or_nzb ?
  !
  is_symmorphic = .NOT.(ANY(ftau(:,1:nsymq) /= 0))
  !
  IF (skip_equivalence_q) THEN
     search_sym = .FALSE.
  ELSE
     search_sym = .TRUE.
     IF (.NOT.is_symmorphic) THEN
        DO isym = 1, nsymq
           search_sym = ( search_sym.AND.(ABS(gi(1,isym))<1.d-8).and.  &
                                         (ABS(gi(2,isym))<1.d-8).and.  &
                                         (ABS(gi(3,isym))<1.d-8) )
        ENDDO
     ENDIF
  ENDIF
  !
  ! 10) Setup the parameters alpha_mix
  !
  DO it = 2, niter_max
     IF (alpha_mix(it).eq.0.d0) alpha_mix(it) = alpha_mix(it - 1)
  ENDDO
  !
  ! 11) Since the order of the S matrices is changed (for q\=0) 
  !     we need to re-initialize d1, d2, d3 to rotate the spherical harmonics
  !
  CALL d_matrix( d1, d2, d3 )
  !
  CALL stop_clock ('hp_setup_q')
  !
  RETURN
  !
END SUBROUTINE hp_setup_q
