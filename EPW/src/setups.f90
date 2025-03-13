  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  ! Adapted from the code PH/phq_setups - Quantum-ESPRESSO group
  !-----------------------------------------------------------------------
  SUBROUTINE setups()
  !-----------------------------------------------------------------------
  !!
  !! EPW setups.
  !!
  !! RM - Nov 2014: Noncolinear case implemented
  !! RM - Nov 2018: Updated based on QE 6.3
  !!
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, ityp
  USE cell_base,     ONLY : at, bg
  USE lsda_mod,      ONLY : nspin, starting_magnetization
  USE scf,           ONLY : v, vrs, vltot, kedtau
  USE gvect,         ONLY : ngm
  USE symm_base,     ONLY : nsym, s, irt, t_rev, time_reversal, sr, &
                            inverse_s
  USE eqv,           ONLY : dmuxc
  USE uspp_param,    ONLY : upf
  USE ep_constants,  ONLY : zero, czero
  USE nlcc_ph,       ONLY : drc
  USE uspp,          ONLY : nlcc_any
  USE control_ph,    ONLY : search_sym, u_from_file
  USE modes,         ONLY : npertx, npert, nirr, nmodes, num_rap_mode, u, name_rap_mode
  USE lr_symm_base,  ONLY : gi, gimq, irotmq, minus_q, nsymq, invsymq, rtau
  USE qpoint,        ONLY : xq
  USE control_flags, ONLY : modenum, noinv
  USE xc_lib,        ONLY : xclib_dft_is
  USE fft_base,      ONLY : dfftp
  USE gvecs,         ONLY : doublegrid
  USE noncollin_module, ONLY : noncolin, domag, m_loc, angle1, angle2, ux, nspin_mag
  !
  IMPLICIT NONE
  !
  LOGICAL :: magnetic_sym
  !! the symmetry operations
  LOGICAL :: symmorphic_or_nzb
  !! Symmorphic operation
  INTEGER :: irr
  !! counter on irrepr
  INTEGER :: na
  !! counter on atoms
  INTEGER :: ierr
  !! Error status
  !
  CALL start_clock('setups')
  !
  ! 1) Computes the total local potential (external+scf) on the smooth grid
  !
  CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
  !
  ! Set non linear core correction stuff
  !
  nlcc_any = ANY(upf(1:ntyp)%nlcc)
  IF (nlcc_any) ALLOCATE(drc(ngm, ntyp))
  !
  !  2) If necessary calculate the local magnetization. This information is
  !      needed in sgama
  !
  IF (noncolin .AND. domag) THEN
    ALLOCATE(m_loc(3, nat), STAT = ierr)
    IF (ierr /= 0) CALL errore('setups', 'Error allocating m_loc', 1)
    DO na = 1, nat
      m_loc(1, na) = starting_magnetization(ityp(na)) * SIN(angle1(ityp(na))) * COS(angle2(ityp(na)))
      m_loc(2, na) = starting_magnetization(ityp(na)) * SIN(angle1(ityp(na))) * SIN(angle2(ityp(na)))
      m_loc(3, na) = starting_magnetization(ityp(na)) * COS(angle1(ityp(na)))
    ENDDO
    ux = zero
    IF (xclib_dft_is('gradient')) THEN
      CALL compute_ux(m_loc,ux,nat)
    ENDIF
    DEALLOCATE(m_loc, STAT = ierr)
    IF (ierr /= 0) CALL errore('setups', 'Error deallocating m_loc', 1)
  ENDIF
  !
  ! 3) Computes the derivative of the xc potential
  !
  ALLOCATE(dmuxc(dfftp%nnr, nspin_mag, nspin_mag), STAT = ierr)
  IF (ierr /= 0) CALL errore('setups', 'Error allocating dmuxc', 1)
  CALL setup_dmuxc()
  !
  ! 3.1) Setup all gradient correction stuff
  !
  CALL setup_dgc
  !
  ! 4) Computes the inverse of each matrix of the crystal symmetry group
  !
  CALL inverse_s()
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  CALL setup_nbnd_occ()
  !
  ! 6) Computes alpha_pv
  !
  CALL setup_alpha_pv()
  !
  ! 7) set all the variables needed to use the pattern representation
  !
  magnetic_sym = noncolin .AND. domag
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  !
  nmodes = 3 * nat
  !
  !   If the code arrives here and nsymq is still 0 the small group of q has
  !   not been calculated by set_nscf because this is a recover run.
  !   We recalculate here the small group of q.
  !
  IF (nsymq == 0) THEN
    CALL set_small_group_of_q(nsymq, invsymq, minus_q)
  ENDIF
  IF (.NOT. time_reversal) THEN
    minus_q = .FALSE.
  ENDIF
  !
  IF (modenum > 0) THEN
    search_sym = .FALSE.
    minus_q = .FALSE.
  ENDIF
  !
  ! Allocate and calculate rtau, the Bravais lattice vector associated to a rotation
  !
  ALLOCATE(rtau(3, 48, nat), STAT = ierr)
  IF (ierr /= 0) CALL errore('setups', 'Error allocating rtau', 1)
  ALLOCATE(npert(3 * nat), STAT = ierr)
  IF (ierr /= 0) CALL errore('setups', 'Error allocating npert', 1)
  CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
  !
  !    and calculate the vectors G associated to the symmetry Sq = q + G
  !    if minus_q is true calculate also irotmq and the G associated to Sq=-g+G
  !
  CALL set_giq(xq, s, nsymq, nsym, irotmq, minus_q, gi, gimq)
  !
  search_sym = search_sym .AND. symmorphic_or_nzb()
  !
  ALLOCATE(num_rap_mode(3 * nat), STAT = ierr)
  IF (ierr /= 0) CALL errore('setups', 'Error allocating num_rap_mode', 1)
  num_rap_mode = -1
  IF (search_sym) THEN
    CALL prepare_sym_analysis(nsymq, sr, t_rev, magnetic_sym)
  ENDIF
  !
  ALLOCATE(name_rap_mode(3 * nat), STAT = ierr)
  IF (ierr /= 0) CALL errore('setups', 'Error allocating name_rap_mode', 1)
  ALLOCATE(u(3 * nat, 3 * nat), STAT = ierr)
  IF (ierr /= 0) CALL errore('setups', 'Error allocating u', 1)
  u(:, :) = czero
  IF (.NOT. u_from_file) THEN
     ! SP: These calls set the u
     CALL find_irrep()
  ENDIF
  CALL find_irrep_sym()
  !
  DEALLOCATE(num_rap_mode, STAT = ierr)
  IF (ierr /= 0) CALL errore('setups', 'Error deallocating num_rap_mode', 1)
  DEALLOCATE(name_rap_mode, STAT = ierr)
  IF (ierr /= 0) CALL errore('setups', 'Error deallocating name_rap_mode', 1)
  !
  !  8) set max perturbation
  !
  npertx = 0
  DO irr = 1, nirr
    npertx = MAX(npertx, npert(irr))
  ENDDO
  !
  CALL stop_clock('setups')
  RETURN
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE setups
  !-----------------------------------------------------------------------
