!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine kcw_q_setup
  !-----------------------------------------------------------------------
  !
  !!  This routine is adapted from phq_setup. I removed all the stuff 
  !!  explicitely related to phonons. 
  ! 
  !
  !!  This subroutine prepares several variables which are needed for the
  !!  LR calculation of the screening coefficients. 
  !!  1) computes the total local potential (external+scf) on the smooth
  !!     grid to be used in h_psi and similia
  !!  2) computes the local magnetization (if necessary)
  !!  3) computes dmuxc (with GC if needed)
  !!  4) set the inverse of every matrix invs
  !!  5) for metals sets the occupied bands
  !!  6) computes alpha_pv
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
  !
  USE kinds,            ONLY : DP
  USE ions_base,        ONLY : nat, ityp
  USE lsda_mod,         ONLY : nspin, starting_magnetization
  USE scf,              ONLY : v, vrs, vltot, kedtau
  USE fft_base,         ONLY : dfftp
  USE gvecs,            ONLY : doublegrid
  USE symm_base,        ONLY : nrot, nsym, s, ft, irt, time_reversal, &
                               inverse_s, d1, d2, d3
  USE lr_symm_base,     ONLY : gi, gimq, irotmq, minus_q, invsymq, nsymq, rtau
  USE qpoint,           ONLY : xq
  USE control_lr,       ONLY : lgamma
  USE noncollin_module, ONLY : domag, noncolin, m_loc, angle1, angle2, ux
  !USE funct,            ONLY : dft_is_gradient
  USE xc_lib,           ONLY : xclib_dft_is
  USE control_kcw,      ONLY : niter, alpha_mix
  USE symm_base,        ONLY : time_reversal
  USE control_flags,    ONLY : noinv
  USE noncollin_module,  ONLY : domag, noncolin, m_loc, angle1, angle2, ux, nspin_lsda, nspin_gga, nspin_mag, npol
  !
  IMPLICIT NONE
  !
  INTEGER :: na, it 
  ! counters
  LOGICAL :: magnetic_sym
  !
  call start_clock ('kcw_q_setup')
  !
  ! 1) Computes the total local potential (external+scf) on the smooth grid
  !
  call set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin_mag, doublegrid)
  !
  ! 2) If necessary calculate the local magnetization. This information is
  !      needed in find_sym
  !
  IF (.not.ALLOCATED(m_loc)) ALLOCATE( m_loc( 3, nat ) )
  IF (noncolin.and.domag) THEN
     DO na = 1, nat
        !
        m_loc(1,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * COS( angle2(ityp(na)) )
        m_loc(2,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * SIN( angle2(ityp(na)) )
        m_loc(3,na) = starting_magnetization(ityp(na)) * &
                      COS( angle1(ityp(na)) )
     END DO
     ux=0.0_DP
     !if (dft_is_gradient()) call compute_ux(m_loc,ux,nat)
     IF ( xclib_dft_is('gradient') ) CALL compute_ux(m_loc,ux,nat)
  ENDIF
  !
  ! 3) Computes the derivative of the XC potential
  !
  call setup_dmuxc()
  !
  ! Setup all gradient correction stuff
  !
  call setup_dgc()
  !
  ! 4) Computes the inverse of each matrix of the crystal symmetry group
  !
  call inverse_s()
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  call setup_nbnd_occ()
  !
  ! 6) Computes alpha_pv
  !
  call setup_alpha_pv()
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
  ! Calculate the vectors G associated to the symmetry Sq = q + G
  ! If minus_q=.true. calculate also irotmq and the G associated to Sq=-q+G
  !
  CALL set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)

  !
  !  set the alpha_mix parameter
  !
  do it = 2, niter
     if (alpha_mix (it) .eq.0.d0) alpha_mix (it) = alpha_mix (it - 1)
     !if (abs(alpha_mix (it)) .lt. 1.0d-6) alpha_mix (it) = alpha_mix (it - 1)
  enddo
  !
  ! NsC: Not sure the next two lines are needed
  magnetic_sym = noncolin .AND. domag
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  !
  CALL stop_clock ('kcw_q_setup')
  !
  RETURN 
  !
END SUBROUTINE kcw_q_setup
