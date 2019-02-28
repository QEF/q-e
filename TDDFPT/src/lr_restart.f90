!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE lr_restart(iter_restart,rflag)
  !---------------------------------------------------------------------
  !
  ! Restart the Lanczos recursion
  !
  ! Modified by Osman Baris Malcioglu (2009)
  ! Modified by Iurii Timrov (2013)
  !
  USE kinds,                ONLY : DP 
  USE io_global,            ONLY : stdout, ionode_id
  USE control_flags,        ONLY : gamma_only
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE io_files,             ONLY : tmp_dir, prefix, diropn, wfc_dir
  USE lr_variables,         ONLY : itermax, evc1, evc1_old, &
                                   restart, nwordrestart, iunrestart,project,nbnd_total,F, &
                                   bgz_suffix, beta_store, gamma_store, zeta_store, norm0, &
                                   lr_verbosity, charge_response, LR_polarization, n_ipol, &
                                   eels, sum_rule
  USE charg_resp,           ONLY : resonance_condition, rho_1_tot,rho_1_tot_im
  USE wvfct,                ONLY : nbnd, npwx
  USE becmod,               ONLY : bec_type, becp, calbec
  USE uspp,                 ONLY : vkb 
  USE io_global,            ONLY : ionode
  USE mp,                   ONLY : mp_bcast
  USE mp_world,             ONLY : world_comm
  USE fft_base,             ONLY : dfftp
  USE noncollin_module,     ONLY : nspin_mag, npol
  USE qpoint,               ONLY : nksq

  IMPLICIT NONE
  !
  INTEGER, INTENT(OUT) :: iter_restart
  LOGICAL, INTENT(OUT) :: rflag
  !
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  ! local variables
  !
  INTEGER :: i,ibnd,ibnd_occ,ibnd_virt,temp
  INTEGER :: ik, ig, ip
  LOGICAL :: exst
  CHARACTER(len=256) :: tempfile, filename, tmp_dir_saved
  INTEGER :: pol_index
  !
  IF (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_restart>")')
  ENDIF
  !
  pol_index = 1
  !
  IF ( n_ipol /= 1 ) pol_index = LR_polarization
  !
  IF (.not.restart) RETURN
  !
  rflag = .false.
  !
  ! Optical case: recompute the kintic-energy g2kin and 
  ! beta functions vkb (needed only in the US case).
  ! Note, this is done only in the gamma_only case,
  ! because in the k-points version all is recomputed
  ! on-the-fly for every k point.
  !
  IF (gamma_only) THEN
     CALL g2_kin(1)
     CALL init_us_2(ngk(1),igk_k(:,1),xk(:,1),vkb)
  ENDIF
  !
  ! Reading Lanczos coefficients
  !
  IF (eels) THEN
    filename = trim(prefix) // trim(bgz_suffix) // trim("dat")
  ELSE
    filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(LR_polarization))
  ENDIF
  tempfile = trim(tmp_dir) // trim(filename)
  !
  INQUIRE (file = tempfile, exist = exst)
  !
  IF (.not.exst) THEN
     !
     WRITE( stdout,*) "WARNING: " // trim(filename) // " does not exist"
     rflag = .true.
     RETURN
     !
  ENDIF
  !
  ! Ionode only reads
  ! Note: ionode file I/O is done in tmp_dir
  !
#if defined(__MPI)
  IF (ionode) THEN
#endif
  !
  ! Read and broadcast beta, gamma, and zeta.
  !
  OPEN (158, file = tempfile, form = 'formatted', status = 'old')
  !
  READ(158,*,end=301,err=303) iter_restart
  !
  IF ( iter_restart >= itermax ) iter_restart = itermax
  !
  READ(158,*,end=301,err=303) norm0(pol_index)
  !
  ! X. Ge: New structure for the pseudo-Hermitian algorithm.
  !
  beta_store(pol_index,1) = norm0(pol_index)
  !
  DO i=1,7
     READ(158,*,end=301,err=303)
  ENDDO
  ! 
  DO i=1,iter_restart-1
     READ(158,*,end=301,err=303) beta_store(pol_index,i+1)
     READ(158,*,end=301,err=303) gamma_store(pol_index,i+1)
     READ(158,*,end=301,err=303) zeta_store (pol_index,:,i)
  ENDDO
  !
  READ(158,*,end=301,err=303) beta_store(pol_index,iter_restart)
  READ(158,*,end=301,err=303) gamma_store(pol_index,iter_restart)
  READ(158,*,end=301,err=303) zeta_store (pol_index,:,iter_restart)
  !
  CLOSE(158)
  !
#if defined(__MPI)
  ENDIF
  CALL mp_bcast (iter_restart, ionode_id, world_comm)
  CALL mp_bcast (norm0(pol_index), ionode_id, world_comm)
  CALL mp_bcast (beta_store(pol_index,:), ionode_id, world_comm)
  CALL mp_bcast (gamma_store(pol_index,:), ionode_id, world_comm)
  CALL mp_bcast (zeta_store(pol_index,:,:), ionode_id, world_comm)
#endif
  !
  ! Optical case: read projection
  !
  IF (project .and. .not.eels) THEN
#if defined(__MPI)
  IF (ionode) THEN
#endif
    !
    filename = trim(prefix) // ".projection." // trim(int_to_char(LR_polarization))
    tempfile = trim(tmp_dir) // trim(filename)
    !
    OPEN (158, file = tempfile, form = 'formatted', status = 'unknown')
    !
    READ(158,*,end=301,err=303) temp
    !
    IF (temp /= iter_restart) CALL errore ('lr_restart', 'Iteration mismatch reading projections', 1 )
    !
    READ(158,*,end=301,err=303) temp   !number of filled bands
    !
    IF (temp /= nbnd) CALL errore ('lr_restart', 'NBND mismatch reading projections', 1 )
    !
    READ(158,*,end=301,err=303) temp !total number of bands
    !
    IF (temp /= nbnd_total) CALL errore ('lr_restart', 'Total number of bands mismatch reading projections', 1 )
    !
    DO ibnd_occ=1,nbnd
       DO ibnd_virt=1,(nbnd_total-nbnd)
        READ(158,*,end=301,err=303) F(ibnd_occ,ibnd_virt,pol_index)
       ENDDO
    ENDDO
    !
    CLOSE(158)
#if defined(__MPI)
  ENDIF
  CALL mp_bcast (F, ionode_id, world_comm)
#endif
  ENDIF
  !
  iter_restart = iter_restart + 1
  !
  ! Parallel reading
  ! Note: Restart files are always in outdir
  ! Reading Lanczos vectors
  !
  nwordrestart = 2 * nbnd * npwx * npol * nksq
  !
  CALL diropn ( iunrestart, 'restart_lanczos.'//trim(int_to_char(LR_polarization)), nwordrestart, exst)
  !
  CALL davcio(evc1(:,:,:,1),nwordrestart,iunrestart,1,-1)
  CALL davcio(evc1(:,:,:,2),nwordrestart,iunrestart,2,-1)
  CALL davcio(evc1_old(:,:,:,1),nwordrestart,iunrestart,3,-1)
  CALL davcio(evc1_old(:,:,:,2),nwordrestart,iunrestart,4,-1)
  !
  CLOSE( unit = iunrestart)
  !
  ! Optical case: read the response charge density
  !
  IF (charge_response == 1 .and. .not.eels) THEN
     IF (resonance_condition) THEN
         CALL diropn ( iunrestart, 'restart_lanczos-rho_tot.'//trim(int_to_char(LR_polarization)), 2*dfftp%nnr*nspin_mag, exst)
        CALL davcio(rho_1_tot_im(:,:),2*dfftp%nnr*nspin_mag,iunrestart,1,-1)
        CLOSE( unit = iunrestart)
    ELSE
        CALL diropn ( iunrestart, 'restart_lanczos-rho_tot.'//trim(int_to_char(LR_polarization)), 2*dfftp%nnr*nspin_mag, exst)
        CALL davcio(rho_1_tot(:,:),2*dfftp%nnr*nspin_mag,iunrestart,1,-1)
        CLOSE( unit = iunrestart)
    ENDIF
 ENDIF
  !
  ! End of all file I/O for restart.
  !
  RETURN
  !
  301 CALL errore ('restart', 'A File is corrupted, file ended unexpectedly', 1 )
  303 CALL errore ('restart', 'A File is corrupted, error in reading data', 1)
  !
END SUBROUTINE lr_restart
