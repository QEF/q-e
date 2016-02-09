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
  ! ... restart the Lanczos recursion
  !---------------------------------------------------------------------
  !
  ! Modified by Osman Baris Malcioglu (2009)
  ! Modified by Iurii Timrov (2013)
  !
  USE kinds,                ONLY : DP 
  USE io_global,            ONLY : stdout, ionode_id
  USE control_flags,        ONLY : gamma_only
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE cell_base,            ONLY : tpiba2
  USE gvect,                ONLY : g
  USE io_files,             ONLY : tmp_dir, prefix, diropn, wfc_dir
  USE lr_variables,         ONLY : itermax,evc1, evc1_new, sevc1_new,evc1_old, &
                                   restart, nwordrestart, iunrestart,project,nbnd_total,F, &
                                   bgz_suffix, beta_store, gamma_store, zeta_store, norm0, &
                                   lr_verbosity, charge_response, LR_polarization, n_ipol, eels, sum_rule
  USE charg_resp,           ONLY : resonance_condition, rho_1_tot,rho_1_tot_im
  USE wvfct,                ONLY : npw, igk, nbnd, g2kin, npwx
  USE becmod,               ONLY : bec_type, becp, calbec
  USE uspp,                 ONLY : vkb, nkb, okvan
  USE io_global,            ONLY : ionode
  USE mp,                   ONLY : mp_bcast
  USE mp_world,             ONLY : world_comm
  USE realus,               ONLY : real_space, invfft_orbital_gamma, initialisation_level, &
                                   fwfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                   v_loc_psir, s_psir_gamma, real_space_debug
  USE fft_base,             ONLY : dfftp
  USE noncollin_module,     ONLY : nspin_mag

  IMPLICIT NONE
  !
  CHARACTER(len=6), EXTERNAL :: int_to_char
  INTEGER, INTENT(out) :: iter_restart
  LOGICAL, INTENT(out) :: rflag
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
  ! Optical case: restarting kintic-energy and ultrasoft
  !
  IF (gamma_only) THEN
     !
     DO ig=1,npwx
        !
        g2kin(ig)=((xk(1,1)+g(1,igk_k(ig,1)))**2 &
                 + (xk(2,1)+g(2,igk_k(ig,1)))**2 &
                 + (xk(3,1)+g(3,igk_k(ig,1)))**2)*tpiba2
        !
     ENDDO
     !
     CALL init_us_2(npw,igk,xk(1,1),vkb)
     !
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
#ifdef __MPI
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
#ifdef __MPI
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
#ifdef __MPI
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
#ifdef __MPI
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
  !nwordrestart = 2 * nbnd * npwx * nks
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
  IF (charge_response == 1 .and. .not. (eels )) THEN
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
  ! Reinitializing sevc1_new vector
  !
  !IF (eels) THEN
  !   CALL lr_recalc_sevc1_new_eels()
  !ELSE
  !   CALL lr_recalc_sevc1_new_optical()
  !ENDIF
  !
  RETURN
  !
  301 CALL errore ('restart', 'A File is corrupted, file ended unexpectedly', 1 )
  303 CALL errore ('restart', 'A File is corrupted, error in reading data', 1)
  !
CONTAINS
  !
SUBROUTINE lr_recalc_sevc1_new_optical()
  !
  ! Optical case: 
  ! This subroutine recalculates sevc1_new(:,:,ik,1) and sevc1_new(:,:,ik,2)
  !
  USE realus,         ONLY : real_space, invfft_orbital_gamma, initialisation_level, &
                             fwfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                             v_loc_psir, s_psir_gamma, real_space_debug
  !
  if (gamma_only) then
     !
     ! Calculation of sevc1_new(:,:,1,1)
     !
     if ( nkb > 0 .and. okvan ) then
        !
        if (real_space_debug>6) then
           !
           do ibnd=1,nbnd,2
              !
              call invfft_orbital_gamma(evc1_new(:,:,1,1),ibnd,nbnd)
              call calbec_rs_gamma(ibnd,nbnd,becp%r)
              call s_psir_gamma(ibnd,nbnd)
              call fwfft_orbital_gamma(sevc1_new(:,:,1,1),ibnd,nbnd)
              !
           enddo
           !
        else
           call calbec(ngk(1),vkb,evc1_new(:,:,1,1),becp)
           call s_psi(npwx,ngk(1),nbnd,evc1_new(:,:,1,1),sevc1_new(:,:,1,1))
        endif
        !
     else
        !
        ! nkb = 0, not real space
        !
        call s_psi(npwx,ngk(1),nbnd,evc1_new(:,:,1,1),sevc1_new(:,:,1,1))
        !
     endif
     !
     ! Calculation of sevc1_new(:,:,1,2)
     !
     if ( nkb > 0 .and. okvan ) then
        ! 
        if (real_space_debug>6) then
           !
           do ibnd=1,nbnd,2
              !
              call invfft_orbital_gamma(evc1_new(:,:,1,2),ibnd,nbnd)
              call calbec_rs_gamma(ibnd,nbnd,becp%r)
              call s_psir_gamma(ibnd,nbnd)
              call fwfft_orbital_gamma(sevc1_new(:,:,1,2),ibnd,nbnd)
              !  
           enddo
           !
        else
           call calbec(ngk(1),vkb,evc1_new(:,:,1,2),becp%r)
           call s_psi(npwx,ngk(1),nbnd,evc1_new(:,:,1,2),sevc1_new(:,:,1,2))
        endif
        !
     else
        !
        ! nkb = 0, not real space
        !
        call s_psi(npwx,ngk(1),nbnd,evc1_new(:,:,1,2),sevc1_new(:,:,1,2))
        !
     endif
     !
  else
     !
     ! K points case: calculation of sevc1_new(:,:,ik,1) and sevc1_new(:,:,ik,2)
     !
     do ik=1,nks
        !
        if ( nkb > 0 .and. okvan ) then
           !
           call init_us_2(ngk(ik),igk_k(1,ik),xk(1,ik),vkb)
           call calbec(ngk(ik), vkb, evc1_new(:,:,ik,1), becp)
           !
        endif
        !
        call s_psi(npwx,ngk(ik),nbnd,evc1_new(:,:,ik,1),sevc1_new(:,:,ik,1))
        !
        if (nkb > 0 .and. okvan) call calbec(ngk(ik), vkb, evc1_new(:,:,ik,2),becp)
        !
        call s_psi(npwx,ngk(ik),nbnd,evc1_new(:,:,ik,2),sevc1_new(:,:,ik,2))
        !
     enddo
     !
  end if
  !
  return
  !
END SUBROUTINE lr_recalc_sevc1_new_optical

SUBROUTINE lr_recalc_sevc1_new_eels()
    !
    ! EELS: This subroutine recalculates sevc1_new(:,:,ik,1) and
    ! sevc1_new(:,:,ik,2)
    !
    use becmod,              only : becp, calbec
    use lr_variables,        only : lr_periodic
    use qpoint,              only : nksq, npwq, igkq, ikks, ikqs
    use gvect,               only : ngm, g
    use wvfct,               only : g2kin
    use gvecw,               only : gcutw
    use control_lr,          only : nbnd_occ

    implicit none
    integer :: ikk, ikq
    !
    do ik = 1, nksq
       !
       if (lr_periodic) then
          ikk = ik
          ikq = ik
       else
          ikk = ikks(ik)
          ikq = ikqs(ik)
       endif
       !
       ! Determination of npwq, igkq; g2kin is used here as a workspace.
       !
       CALL gk_sort( xk(1,ikq), ngm, g, gcutw, npwq, igkq, g2kin)
       !
       if ( okvan .and. nkb > 0 ) then
          !
          ! Calculate beta-functions vkb at point k+q
          !
          call init_us_2(npwq, igkq, xk(1,ikq), vkb)
          !
          ! Calculate the product of beta-functions vkb with
          ! the response orbitals evc1_new : becp%k = <vkb|evc1_new>
          !
          call calbec(npwq, vkb, evc1_new(:,:,ik,1), becp, nbnd_occ(ikk))
          !
       endif
       !
       call s_psi(npwx, npwq, nbnd_occ(ikk), evc1_new(:,:,ik,1), sevc1_new(:,:,ik,1))
       !
       if ( okvan .and. nkb > 0 ) then
          !
          ! Calculate the product of beta-functions vkb with
          ! the response orbitals evc1_new : becp%k = <vkb|evc1_new>
          !
          call calbec(npwq, vkb, evc1_new(:,:,ik,2), becp, nbnd_occ(ikk))
          !
       endif
       !
       call s_psi(npwx, npwq, nbnd_occ(ikk), evc1_new(:,:,ik,2), sevc1_new(:,:,ik,2))
       !
    enddo
    !
    return
    !
END SUBROUTINE lr_recalc_sevc1_new_eels
  !
END SUBROUTINE lr_restart
!-----------------------------------------------------------------------
