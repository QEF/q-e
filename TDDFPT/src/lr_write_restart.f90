!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
SUBROUTINE lr_write_restart()
  !---------------------------------------------------------------------
  ! 
  ! This subroutine reads in and stores vectors necessary to
  ! restart the Lanczos recursion.
  !
  ! Modified by Osman Baris Malcioglu (2009)
  ! Modified by Xiaochuan Ge (May. 2013) to adapt pseudo-hermitian
  !
  USE kinds,                ONLY : DP
  USE io_files,             ONLY : tmp_dir, prefix, diropn
  USE lr_variables,         ONLY : beta_store, gamma_store, zeta_store, norm0, &
                                   LR_polarization, LR_iteration, n_ipol,F,project,&
                                   evc1,evc1_new,evc1_old,iunrestart, nwordrestart, &
                                   nbnd_total, charge_response,lr_verbosity,&
                                   bgz_suffix, eels, q1, q2, q3, sum_rule, tmp_dir_lr
  USE charg_resp,           ONLY : resonance_condition, rho_1_tot, rho_1_tot_im
  USE wvfct,                ONLY : nbnd, npwx
  USE fft_base,             ONLY : dfftp
  USE io_global,            ONLY : ionode, stdout
  USE klist,                ONLY : nks, nelec
  USE noncollin_module,     ONLY : nspin_mag, noncolin, npol
  use lsda_mod,             ONLY : nspin
  USE cell_base,            ONLY : alat, omega
  USE qpoint,               ONLY : nksq
  !
  IMPLICIT NONE
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  ! local variables
  !
  INTEGER :: i, j, pol_index,ibnd_occ,ibnd_virt
  CHARACTER(len=256) :: tempfile, filename
  LOGICAL :: exst
  real(kind=dp) :: degspin
  !
  IF (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_write_restart>")')
  ENDIF
  !
  ! Note: ionode only operations are carried out in tmp_dir not wfc_dir
  !
  ! If there is only one polarization dir, storage is one rank less.
  !
  pol_index = 1
  ! 
  IF ( n_ipol /= 1 ) pol_index = LR_polarization
  !
  IF (eels) tmp_dir = tmp_dir_lr
  !
#if defined(__MPI)
  IF (ionode) THEN
#endif
  !
  ! Writing beta, gamma and zeta coefficients.
  !
  IF (eels) THEN
     filename = trim(prefix) // trim(bgz_suffix) // trim("dat")
  ELSE
     filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(LR_polarization))
  ENDIF
  tempfile = trim(tmp_dir) // trim(filename)
  !
  OPEN (158, file = tempfile, form = 'formatted', status = 'unknown')
  WRITE(158,*) LR_iteration
  !
  norm0(pol_index) = beta_store(pol_index,1)
  WRITE(158,*) norm0(pol_index)
  !
  IF (nspin==2) THEN
        degspin = 1.0d0
  ELSE
        degspin = 2.0d0
  ENDIF
  IF (noncolin) degspin = 1.0d0
  !
  ! Write the degenaracy wrt spin
  !
  WRITE(158,*) degspin
  !
  ! ------ Needed for EELS ----------
  !
  ! Write the lattice parameter
  !
  WRITE(158,*) alat
  !
  ! Write the unit-cell volume
  !
  WRITE(158,*) omega
  !
  ! Write the number of valence (and semicore electrons) in the unit cell
  !
  WRITE(158,*) nelec
  !
  ! Write the components of the transferred momentum
  !
  WRITE(158,*) q1
  WRITE(158,*) q2
  WRITE(158,*) q3
  !
  !-----------------------------------
  !
  DO i=1,LR_iteration-1
     !
     WRITE(158,*) beta_store(pol_index,i+1)
     WRITE(158,*) gamma_store(pol_index,i+1)
     !
     ! This is absolutely necessary for cross platform compatibilty
     !
     DO j=1,n_ipol
      WRITE(158,*) zeta_store (pol_index,j,i)
     ENDDO
     !
  ENDDO
  !
  ! X. Ge: Compatable with the old version. The beta & gamma will not be used in 
  ! the spectrum calculation.
  !
  WRITE(158,*) beta_store(pol_index,LR_iteration)             
  WRITE(158,*) gamma_store(pol_index,LR_iteration)             
  DO j=1,n_ipol                                               
    WRITE(158,*) zeta_store (pol_index,j,LR_iteration)        
  ENDDO
  !
  CLOSE(158)
  !
  ! Optical case: writing F
  !
  IF (project .AND. .NOT.eels) THEN
     !
     filename = trim(prefix) // ".projection." // trim(int_to_char(LR_polarization))
     tempfile = trim(tmp_dir) // trim(filename)
     !
     OPEN (158, file = tempfile, form = 'formatted', status = 'unknown')
     WRITE(158,*) LR_iteration
     WRITE(158,*) nbnd        ! number of filled bands
     WRITE(158,*) nbnd_total  !total number of bands
     !
     DO ibnd_occ=1,nbnd
        DO ibnd_virt=1,(nbnd_total-nbnd)
           WRITE(158,*) F(ibnd_occ,ibnd_virt,pol_index)
        ENDDO
     ENDDO
     !
     CLOSE(158)
     !
  ENDIF
  !
#if defined(__MPI)
  ENDIF
#endif
    !
    ! Parallel writing operations
    !
    ! Note: Restart files are writen in outdir.
    ! If you do not want them to be written,
    ! just disable restart saving completely.
    !
    ! Writing wavefuncion files for restart
    !
    nwordrestart = 2 * nbnd * npwx * npol * nksq
    !
    CALL diropn ( iunrestart, 'restart_lanczos.'//trim(int_to_char(LR_polarization)), nwordrestart, exst)
    !
    CALL davcio(evc1(:,:,:,1),nwordrestart,iunrestart,1,1)
    CALL davcio(evc1(:,:,:,2),nwordrestart,iunrestart,2,1)
    CALL davcio(evc1_old(:,:,:,1),nwordrestart,iunrestart,3,1)
    CALL davcio(evc1_old(:,:,:,2),nwordrestart,iunrestart,4,1)
    !
    CLOSE( unit = iunrestart)
    !
    ! Optical case: Writing charge response density for restart
    !
    IF (charge_response == 1 .AND. .NOT.eels) THEN
       !
       IF (resonance_condition) THEN
          CALL diropn ( iunrestart, 'restart_lanczos-rho_tot.'//trim(int_to_char(LR_polarization)), 2*dfftp%nnr*nspin_mag, exst)
          CALL davcio(rho_1_tot_im(:,:),2*dfftp%nnr*nspin_mag,iunrestart,1,1)
          CLOSE( unit = iunrestart)
       ELSE
          CALL diropn ( iunrestart, 'restart_lanczos-rho_tot.'//trim(int_to_char(LR_polarization)), 2*dfftp%nnr*nspin_mag, exst)
          CALL davcio(rho_1_tot(:,:),2*dfftp%nnr*nspin_mag,iunrestart,1,1)
          CLOSE( unit = iunrestart)
       ENDIF
       !
    ENDIF
    IF (sum_rule == -2 .AND. .NOT.eels) THEN
       !
       CALL diropn ( iunrestart, 'restart_lanczos-sum-2.'//trim(int_to_char(LR_polarization)), 2*dfftp%nnr*nspin_mag, exst)
       CALL davcio(rho_1_tot_im(:,:),2*dfftp%nnr*nspin_mag,iunrestart,1,1)
       CLOSE( unit = iunrestart)
       !
    ENDIF

    !
    RETURN
    !
END SUBROUTINE lr_write_restart
!-----------------------------------------------------------------------
