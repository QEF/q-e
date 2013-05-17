!
! Copyright (C) 2004-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------

SUBROUTINE lr_write_restart()
  !---------------------------------------------------------------------
  ! ... reads in and stores the vectors necessary to
  ! ... restart the Lanczos recursion
  !---------------------------------------------------------------------
  !
  ! Modified by Osman Baris Malcioglu (2009)
  ! Modified by Xiaochuan Ge (May. 2013) to adapt pseudo-hermitian
  !
  USE io_files,             ONLY : tmp_dir, prefix, diropn
  USE lr_variables,         ONLY : beta_store, gamma_store, zeta_store, norm0, &
                                   LR_polarization, LR_iteration, n_ipol,F,project,&
                                   evc1,evc1_new,evc1_old,iunrestart, nwordrestart, &
                                   nbnd_total, charge_response,lr_verbosity,&
                                   bgz_suffix
  USE charg_resp,           ONLY : resonance_condition, rho_1_tot, rho_1_tot_im
  USE wvfct,                ONLY : nbnd, npwx, npw
  USE fft_base,             ONLY : dfftp
  USE io_global,            ONLY : ionode
  USE klist,                ONLY : nks
  USE noncollin_module,     ONLY : nspin_mag
  USE io_global,      ONLY : stdout
  !
  IMPLICIT NONE
  !
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  !integer, intent(in) :: pol, iter
  !
  ! local variables
  !
  INTEGER :: i, j, pol_index,ibnd_occ,ibnd_virt
  CHARACTER(len=256) :: tempfile, filename
  LOGICAL :: exst
  !
  IF (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_write_restart>")')
  ENDIF

  !
  !ionode only operations:
  ! Note: ionode only operations are carried out in tmp_dir not wfc_dir
  !
  pol_index=1 !if there is only one polarization dir, storage is one rank less
  IF ( n_ipol /= 1 ) pol_index=LR_polarization

#ifdef __MPI
  IF (ionode) THEN
#endif
  !
  !Writing beta gamma and zeta
  !
  !
  filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(LR_polarization))
  tempfile = trim(tmp_dir) // trim(filename)
  !
  !
  OPEN (158, file = tempfile, form = 'formatted', status = 'unknown')
  !
  WRITE(158,*) LR_iteration
  !
  norm0(pol_index)=beta_store(pol_index,1)
  WRITE(158,*) norm0(pol_index)
  !
  DO i=1,LR_iteration-1
     !
     WRITE(158,*) beta_store(pol_index,i+1)
     WRITE(158,*) gamma_store(pol_index,i+1)
     !This is absolutely necessary for cross platform compatibilty
     DO j=1,n_ipol
      WRITE(158,*) zeta_store (pol_index,j,i)
     ENDDO
     !
  ENDDO
  ! XC: compatable with the old version. The beta&gamma will not be used in 
  ! spectrum calculation
  WRITE(158,*) beta_store(pol_index,LR_iteration)             
  WRITE(158,*) gamma_store(pol_index,LR_iteration)             
  DO j=1,n_ipol                                               
    WRITE(158,*) zeta_store (pol_index,j,LR_iteration)        
  ENDDO
  !
  CLOSE(158)
  !
  !Writing F
  !
  IF (project) THEN
    filename = trim(prefix) // ".projection." // trim(int_to_char(LR_polarization))
    tempfile = trim(tmp_dir) // trim(filename)
    !
    !
    OPEN (158, file = tempfile, form = 'formatted', status = 'unknown')
    !
    WRITE(158,*) LR_iteration
    !
    WRITE(158,*) nbnd   !number of filled bands
    !
    WRITE(158,*) nbnd_total !total number of bands
    !
    DO ibnd_occ=1,nbnd
       DO ibnd_virt=1,(nbnd_total-nbnd)
        WRITE(158,*) F(ibnd_occ,ibnd_virt,pol_index)
       ENDDO
    ENDDO
    !
    CLOSE(158)
  ENDIF
  !
#ifdef __MPI
  ENDIF
#endif
    !
    ! Parallel writing operations
    !
    ! Note: Restart files are writen in outdir, if you do not want them to be
    ! written just disable restart saving completely
    !
    !
    ! Writing wavefuncion files for restart
    !
       !
       nwordrestart = 2 * nbnd * npwx * nks
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
    ! Writing charge response density for restart
    !
       IF (charge_response == 1 ) THEN
        IF (resonance_condition) THEN
         CALL diropn ( iunrestart, 'restart_lanczos-rho_tot.'//trim(int_to_char(LR_polarization)), 2*dfftp%nnr*nspin_mag, exst)
         CALL davcio(rho_1_tot_im(:,:),2*dfftp%nnr*nspin_mag,iunrestart,1,1)
         CLOSE( unit = iunrestart)
        ELSE
         CALL diropn ( iunrestart, 'restart_lanczos-rho_tot.'//trim(int_to_char(LR_polarization)), 2*dfftp%nnr*nspin_mag, exst)
         CALL davcio(rho_1_tot(:,:),2*dfftp%nnr*nspin_mag,iunrestart,1,1)
         CLOSE( unit = iunrestart)
        ENDIF
       ENDIF
       !
  RETURN
END SUBROUTINE lr_write_restart
!-----------------------------------------------------------------------
