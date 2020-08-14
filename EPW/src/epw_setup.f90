  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  ! Adapted from the code PH/phq_setup - Quantum-ESPRESSO group
  !-----------------------------------------------------------------------
  SUBROUTINE epw_setup()
  !-----------------------------------------------------------------------
  !!
  !! EPW setup.
  !!
  !! RM - Nov 2014: Noncolinear case implemented
  !! RM - Nov 2018: Updated based on QE 6.3
  !!
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, ityp
  USE cell_base,     ONLY : at, bg
  USE klist,         ONLY : nkstot
  USE lsda_mod,      ONLY : nspin, starting_magnetization
  USE scf,           ONLY : v, vrs, vltot, kedtau
  USE gvect,         ONLY : ngm
  USE symm_base,     ONLY : nsym, s, irt, t_rev, time_reversal, sr, &
                            inverse_s
  USE eqv,           ONLY : dmuxc
  USE uspp_param,    ONLY : upf
  USE spin_orb,      ONLY : domag
  USE constants_epw, ONLY : zero, eps5, czero, ryd2ev, kelvin2ev
  USE nlcc_ph,       ONLY : drc
  USE uspp,          ONLY : nlcc_any
  USE control_ph,    ONLY : search_sym, u_from_file
  USE modes,         ONLY : npertx, npert, nirr, nmodes, num_rap_mode, u, name_rap_mode
  USE lr_symm_base,  ONLY : gi, gimq, irotmq, minus_q, nsymq, invsymq, rtau
  USE qpoint,        ONLY : xq
  USE control_flags, ONLY : modenum, noinv
  USE funct,         ONLY : dft_is_gradient
  USE mp_global,     ONLY : world_comm
  USE mp,            ONLY : mp_bcast
  USE epwcom,        ONLY : scattering, nkc1, nkc2, nkc3
  USE klist_epw,     ONLY : xk_cryst
  USE fft_base,      ONLY : dfftp
  USE gvecs,         ONLY : doublegrid
  USE noncollin_module, ONLY : noncolin, m_loc, angle1, angle2, ux, nspin_mag
  ! ---------------------------------------------------------------------------------
  ! Added for polaron calculations. Originally by Danny Sio, modified by Chao Lian.
  ! Shell implementation for future use.
  USE epwcom,        ONLY: polaron_wf
  ! ---------------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  LOGICAL :: magnetic_sym
  !! the symmetry operations
  LOGICAL :: symmorphic_or_nzb
  !! Symmorphic operation
  INTEGER :: jk
  !! counter on k points
  INTEGER :: irr
  !! counter on irrepr
  INTEGER :: na
  !! counter on atoms
  INTEGER :: ierr
  !! Error status
  REAL(KIND = DP) :: xx_c, yy_c, zz_c
  !! k-points in crystal coords. in multiple of nkc1, nkc2, nkc3
  !
  CALL start_clock('epw_setup')
  !
  IF (.NOT. polaron_wf)  THEN
    DO jk = 1, nkstot
      xx_c = xk_cryst(1, jk) * nkc1
      yy_c = xk_cryst(2, jk) * nkc2
      zz_c = xk_cryst(3, jk) * nkc3
      !
      ! check that the k-mesh was defined in the positive region of 1st BZ
      !
      IF (xx_c < -eps5 .OR. yy_c < -eps5 .OR. zz_c < -eps5) &
        CALL errore('epw_setup', 'coarse k-mesh needs to be strictly positive in 1st BZ', 1)
      !
    ENDDO
  ENDIF ! not polaron_wf
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
    IF (ierr /= 0) CALL errore('epw_setup', 'Error allocating m_loc', 1)
    DO na = 1, nat
      m_loc(1, na) = starting_magnetization(ityp(na)) * SIN(angle1(ityp(na))) * COS(angle2(ityp(na)))
      m_loc(2, na) = starting_magnetization(ityp(na)) * SIN(angle1(ityp(na))) * SIN(angle2(ityp(na)))
      m_loc(3, na) = starting_magnetization(ityp(na)) * COS(angle1(ityp(na)))
    ENDDO
    ux = zero
    IF (dft_is_gradient()) THEN
      CALL compute_ux(m_loc,ux,nat)
    ENDIF
    DEALLOCATE(m_loc, STAT = ierr)
    IF (ierr /= 0) CALL errore('epw_setup', 'Error deallocating m_loc', 1)
  ENDIF
  !
  ! 3) Computes the derivative of the xc potential
  !
  ALLOCATE(dmuxc(dfftp%nnr, nspin_mag, nspin_mag), STAT = ierr)
  IF (ierr /= 0) CALL errore('epw_setup', 'Error allocating dmuxc', 1)
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
  IF (ierr /= 0) CALL errore('epw_setup', 'Error allocating rtau', 1)
  ALLOCATE(npert(3 * nat), STAT = ierr)
  IF (ierr /= 0) CALL errore('epw_setup', 'Error allocating npert', 1)
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
  IF (ierr /= 0) CALL errore('epw_setup', 'Error allocating num_rap_mode', 1)
  num_rap_mode = -1
  IF (search_sym) THEN
    CALL prepare_sym_analysis(nsymq, sr, t_rev, magnetic_sym)
  ENDIF
  !
  ALLOCATE(name_rap_mode(3 * nat), STAT = ierr)
  IF (ierr /= 0) CALL errore('epw_setup', 'Error allocating name_rap_mode', 1)
  ALLOCATE(u(3 * nat, 3 * nat), STAT = ierr)
  IF (ierr /= 0) CALL errore('epw_setup', 'Error allocating u', 1)
  u(:, :) = czero
  IF (.NOT. u_from_file) THEN
     ! SP: These calls set the u
     CALL find_irrep()
  ENDIF
  CALL find_irrep_sym()
  !
  DEALLOCATE(num_rap_mode, STAT = ierr)
  IF (ierr /= 0) CALL errore('epw_setup', 'Error deallocating num_rap_mode', 1)
  DEALLOCATE(name_rap_mode, STAT = ierr)
  IF (ierr /= 0) CALL errore('epw_setup', 'Error deallocating name_rap_mode', 1)
  !
  !  8) set max perturbation
  !
  npertx = 0
  DO irr = 1, nirr
    npertx = MAX(npertx, npert(irr))
  ENDDO
  !
  CALL stop_clock('epw_setup')
  RETURN
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE epw_setup
  !-----------------------------------------------------------------------
