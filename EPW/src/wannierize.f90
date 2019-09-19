  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------------------
  SUBROUTINE wann_run()
  !---------------------------------------------------------------------
  !!
  !!  This is the routine which controls the w90 run.  Primarily,        
  !!  we get the phases to remove degeneracies in the wfs, and 
  !!  call pw2wan90epw 
  !!  
  !
  USE kinds,          ONLY : DP
  USE io_global,      ONLY : stdout, ionode_id
  USE wvfct,          ONLY : nbnd
  USE epwcom,         ONLY : nkc1, nkc2, nkc3
  USE pwcom,          ONLY : nkstot 
  USE klist_epw,      ONLY : xk_cryst           
  USE wannierEPW,     ONLY : mp_grid, n_wannier, kpt_latt
  USE mp,             ONLY : mp_bcast
  USE mp_world,       ONLY : world_comm
  !
  IMPLICIT NONE
  !
  ! Local variables
  INTEGER :: num_kpts
  !! number of k-points in wannierization 
  INTEGER :: ierr
  !! Error status
  !
  CALL start_clock('WANNIER')
  ! 
  mp_grid(1) = nkc1
  mp_grid(2) = nkc2
  mp_grid(3) = nkc3
  num_kpts = mp_grid(1) * mp_grid(2) * mp_grid(3)
  !
  IF (num_kpts /= nkstot) CALL errore('wannierize', 'inconsistent nscf and elph k-grids', 1) 
  IF (nbnd < n_wannier)  CALL errore('wannierize', 'Must have as many or more bands than Wannier functions', 1) 
  !
  ALLOCATE(kpt_latt(3, num_kpts), STAT = ierr)
  IF (ierr /= 0) CALL errore('wann_run', 'Error allocating kpt_latt', 1)
  !
  WRITE(stdout, '(5x,a)') REPEAT("-", 67)
  WRITE(stdout, '(a, i2, a, i2, a, i2, a)') "     Wannierization on ", nkc1, " x ", nkc2, " x ", nkc3 , " electronic grid"
  WRITE(stdout, '(5x, a)') REPEAT("-",67)
  !
  kpt_latt = xk_cryst(:, 1:num_kpts)
  CALL mp_bcast(kpt_latt, ionode_id, world_comm)
  !
  ! write the short input file for the wannier90 code
  !
  CALL write_winfil()
  !
  ! run the wannier90 code to create MLWFs
  !
  CALL pw2wan90epw()
  !
  ! project the Wannier functions onto energy space
  !
  DEALLOCATE(kpt_latt, STAT = ierr)
  IF (ierr /= 0) CALL errore('wann_run', 'Error deallocating kpt_latt', 1)
  !
  WRITE(stdout, '(5x, a)') REPEAT("-", 67)
  CALL print_clock('WANNIER')
  WRITE(stdout, '(5x, a)') REPEAT("-", 67)
  !
  !------------------------------------------------------------
  END SUBROUTINE wann_run
  !------------------------------------------------------------
  !
  !------------------------------------------------------------
  SUBROUTINE write_winfil()
  !------------------------------------------------------------
  !!
  !!  This routine writes the prefix.win file which wannier90.x
  !!  needs to run.  Primarily it contains information about the 
  !!  windows used for the disentanglement, and the initial projections.
  !!  JN - 10/2008  projections now in elph.in file  
  !
  USE kinds,       ONLY : DP
  USE io_files,    ONLY : prefix
  USE io_var,      ONLY : iuwinfil
  USE io_global,   ONLY : meta_ionode
  USE pwcom,       ONLY : et, nbnd, nkstot, nks
  USE epwcom,      ONLY : nbndsub, nwanxx, proj, iprint, dis_win_min, &
                          dis_win_max, dis_froz_min, dis_froz_max, num_iter, &
                          bands_skipped, wdata, vme, auto_projections
  USE constants_epw, ONLY : ryd2ev
  USE poolgathering, ONLY : poolgather
  !
  IMPLICIT NONE
  !
  ! Local variables
  LOGICAL :: random
  !! Random 
  INTEGER :: i
  !! Band index
  REAL(KIND = DP) :: et_tmp(nbnd, nkstot)
  !! eigenvalues on full coarse k-mesh
  !
  CALL poolgather(nbnd, nkstot, nks, et(1:nbnd, 1:nks), et_tmp)
  et_tmp = et_tmp * ryd2ev
  !
  IF (meta_ionode) THEN
    !
    IF (nbndsub > nwanxx) CALL errore('write_winfil', 'Too many wannier bands', nbndsub)
    !
    OPEN(UNIT = iuwinfil, FILE = TRIM(prefix) // ".win", FORM = 'formatted')
    !    
    !  more input and options for interfacing with w90 can/will be added later
    IF (auto_projections) THEN
      WRITE(iuwinfil, '(a)') "auto_projections = .true."
    ELSE
      WRITE(iuwinfil, '(a)') "begin projections"
      !
      random = .TRUE.
      DO i = 1, nbndsub
        IF (proj(i) /= ' ') THEN
          WRITE(iuwinfil, *) TRIM(proj(i))
          random = .FALSE.
        ENDIF
      ENDDO
      !
      IF (random) WRITE(iuwinfil, '(a)') "random" 
      !
      WRITE(iuwinfil, '(a)') "end projections"
    ENDIF
    !
    IF (bands_skipped /= ' ') WRITE(iuwinfil, *) bands_skipped
    !
    WRITE(iuwinfil, '("num_wann = ", i3)') nbndsub
    WRITE(iuwinfil, '("iprint = ", i3)') iprint
    !
    ! SP: You can have more bands in nscf.in than in 
    !     nbndskip+nbndsub. In which case the dis_win_max can be larger than 
    !     nbndskip+nbndsub. This is crucial for disantanglement. 
    IF (dis_froz_min < MINVAL(et_tmp)) dis_froz_min = MINVAL(et_tmp)
    IF (dis_froz_max > MAXVAL(et_tmp)) dis_froz_max = MAXVAL(et_tmp)
    !
    WRITE(iuwinfil, '("dis_win_min = ", f18.12)')  dis_win_min
    WRITE(iuwinfil, '("dis_win_max = ", f18.12)')  dis_win_max
    WRITE(iuwinfil, '("dis_froz_min = ", f18.12)') dis_froz_min
    WRITE(iuwinfil, '("dis_froz_max = ", f18.12)') dis_froz_max
    WRITE(iuwinfil, '("num_iter = ", i7)')         num_iter
    IF (vme) WRITE(iuwinfil, '(a)') "write_bvec = .true."
    !
    ! write any extra parameters to the prefix.win file
    DO i = 1, nwanxx
      IF (wdata(i) /= ' ') WRITE(iuwinfil, *) wdata(i)
    ENDDO
    !
    CLOSE(iuwinfil)
    !
  ENDIF ! meta_ionode
  !
  !------------------------------------------------------------
  END SUBROUTINE write_winfil
  !------------------------------------------------------------
  ! 
  !------------------------------------------------------------
  SUBROUTINE proj_w90()
  !------------------------------------------------------------
  !!
  !! This routine computes the energy projections of
  !! the computed Wannier functions
  !! 07/2010  Needs work.  Right now this sub is nearly worthless  
  !
  USE kinds,       ONLY : DP
  USE io_files,    ONLY : prefix 
  USE io_var,      ONLY : iuprojfil
  USE mp_global,   ONLY : inter_pool_comm
  USE io_global,   ONLY : stdout, meta_ionode
  USE mp,          ONLY : mp_sum
  USE epwcom,      ONLY : dis_win_max, dis_win_min
  USE constants_epw, ONLY : ryd2ev, zero
  USE wannierEPW,  ONLY : n_wannier
  USE wvfct,       ONLY : nbnd, et
  USE klist,       ONLY : nks, nkstot
  USE elph2,       ONLY : xkq
  !
  IMPLICIT NONE
  !
  LOGICAL :: lwin(nbnd, nks)
  !! Bands within the energy window
  LOGICAL :: lwinq(nbnd, nks)
  !! k+q bands within the energy window
  LOGICAL :: exband(nbnd)
  !! Exclude bands
  INTEGER :: ik
  !! K-point index
  INTEGER :: ibnd
  !! 
  INTEGER :: ne 
  !!
  INTEGER :: ie
  !!
  INTEGER :: iwann
  !!
  INTEGER :: ierr
  !! Error status
  REAL(KIND = DP) :: dE
  !! 
  REAL(KIND = DP) :: sigma
  !! 
  REAL(KIND = DP) :: argv
  !! 
  REAL(KIND = DP) :: en
  !! 
  REAL(KIND = DP) :: xxq(3)
  !! Current q-point
  REAL(KIND = DP), ALLOCATABLE :: proj_wf(:, :)
  !! Projection
  COMPLEX(KIND = DP), ALLOCATABLE :: cu(:, :, :)
  !! k rotation matrix
  COMPLEX(KIND = DP), ALLOCATABLE :: cuq(:, :, :)
  !! k+q rotation matrix
  !
  WRITE(stdout, '(5x, "Computing energy projections")')
  ! dummy var
  xxq = zero
  !
  ! tmp value
  dE = 0.05d0
  sigma = 2.d0 * dE
  !
  ! maxvalue = dis_win_max + 1
  ! minvalue = dis_win_min - 1
  ne = INT((dis_win_max - dis_win_min + 1) / dE)
  IF (ne < 1) CALL errore('proj_wan', 'Problem with disentanglement window', 1)
  !
  ALLOCATE(proj_wf(n_wannier, ne + 1), STAT = ierr)
  IF (ierr /= 0) CALL errore('proj_w90', 'Error allocating proj_wf', 1)
  proj_wf = 0.d0
  !
  ALLOCATE(cu(nbnd, n_wannier, nks), STAT = ierr)
  IF (ierr /= 0) CALL errore('proj_w90', 'Error allocating cu', 1)
  ALLOCATE(cuq(nbnd, n_wannier, nks), STAT = ierr)
  IF (ierr /= 0) CALL errore('proj_w90', 'Error allocating cuq', 1)
  ALLOCATE(xkq(3, nks), STAT = ierr)
  IF (ierr /= 0) CALL errore('proj_w90', 'Error allocating xkq', 1)
  CALL loadumat(nbnd, n_wannier, nks, nkstot, xxq, cu, cuq, lwin, lwinq, exband) 
  DEALLOCATE(xkq, STAT = ierr)
  IF (ierr /= 0) CALL errore('proj_w90', 'Error deallocating xkq', 1)
  !
  DO iwann = 1, n_wannier
    DO ie = 1, ne
      en = DBLE(ie) / DBLE(ne) * (dis_win_max - dis_win_min + 1.0) + dis_win_min
      DO ik = 1, nks
        DO ibnd = 1, nbnd
          argv = (et(ibnd, ik) * ryd2ev - en)**2 / (2 * sigma **2)
          proj_wf(iwann, ie) = proj_wf(iwann, ie) + EXP(-argv) * REAL(cu(ibnd, iwann, ik) * CONJG(cu(ibnd, iwann, ik)))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  ! sum the contributions from all k-points
  CALL mp_sum(proj_wf, inter_pool_comm)
  !
  IF (meta_ionode) THEN
    !
    OPEN(UNIT = iuprojfil, FILE = TRIM(prefix) // ".projw90", FORM = 'formatted')
    !
    WRITE(iuprojfil, '(5x, "Wannier energy projections")')
    !
    DO ie = 1, ne
      en =  DBLE(ie) / DBLE(ne) * (dis_win_max - dis_win_min + 1) + dis_win_min
      WRITE(iuprojfil, '(f9.3, 25f8.4)') en, proj_wf(:, ie)
    ENDDO
    !
    CLOSE(iuprojfil)
  ENDIF
  DEALLOCATE(proj_wf, STAT = ierr)
  IF (ierr /= 0) CALL errore('proj_w90', 'Error deallocating proj_wf', 1)
  DEALLOCATE(cu, STAT = ierr)
  IF (ierr /= 0) CALL errore('proj_w90', 'Error deallocating cu', 1)
  DEALLOCATE(cuq, STAT = ierr)
  IF (ierr /= 0) CALL errore('proj_w90', 'Error deallocating cuq', 1)
  !
  !------------------------------------------------------------
  END SUBROUTINE proj_w90
  !------------------------------------------------------------
