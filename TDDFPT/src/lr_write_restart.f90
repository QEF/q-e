!-----------------------------------------------------------------------
subroutine lr_write_restart()
  !---------------------------------------------------------------------
  ! ... reads in and stores the vectors necessary to 
  ! ... restart the Lanczos recursion
  !---------------------------------------------------------------------
  !
  ! Modified by Osman Baris Malcioglu (2009)
#include "f_defs.h"
  !
  use io_files,             only : tmp_dir, prefix, diropn
  use lr_variables,         only : beta_store, gamma_store, zeta_store, norm0, &
                                   LR_polarization, LR_iteration, n_ipol,F,project,&
                                   evc1,evc1_new,iunrestart, nwordrestart, rho_1_tot, rho_1_tot_im, &
                                   nbnd_total, charge_response,lr_verbosity,&
                                   bgz_suffix
  use charg_resp,           only : resonance_condition
  use wvfct,                only : nbnd, npwx, npw
  use grid_dimensions,      only : nrxx
  USE io_global,            ONLY : ionode
  use klist,                only : nks
  USE noncollin_module,     ONLY : nspin_mag
  USE io_global,      ONLY : stdout
  !
  implicit none
  !
  character(len=6), external :: int_to_char
  !
  !integer, intent(in) :: pol, iter
  !
  ! local variables
  !
  integer :: i, j, pol_index,ibnd_occ,ibnd_virt
  character(len=256) :: tempfile, filename
  logical :: exst
  !
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_write_restart>")')
  endif
  
  !
  !ionode only operations:
  ! Note: ionode only operations are carried out in tmp_dir not wfc_dir
  !
  pol_index=1 !if there is only one polarization dir, storage is one rank less
  if ( n_ipol /= 1 ) pol_index=LR_polarization
  
#ifdef __PARA
  if (ionode) then
#endif
  !
  !Writing beta gamma and zeta
  !
  !
  filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(LR_polarization))
  tempfile = trim(tmp_dir) // trim(filename)
  !
  !
  open (158, file = tempfile, form = 'formatted', status = 'unknown')
  !
  write(158,*) LR_iteration
  !
  write(158,*) norm0(pol_index)
  !
  do i=1,LR_iteration
     !
     write(158,*) beta_store(pol_index,i)
     write(158,*) gamma_store(pol_index,i)
     !This is absolutely necessary for cross platform compatibilty
     do j=1,n_ipol                                    
      write(158,*) zeta_store (pol_index,j,i)
     end do
     !
  end do
  !
  close(158)
  !
  !Writing F
  !
  if (project) then
    filename = trim(prefix) // ".projection." // trim(int_to_char(LR_polarization))
    tempfile = trim(tmp_dir) // trim(filename)
    !
    !
    open (158, file = tempfile, form = 'formatted', status = 'unknown')
    !
    write(158,*) LR_iteration
    !
    write(158,*) nbnd   !number of filled bands
    !
    write(158,*) nbnd_total !total number of bands
    !
    do ibnd_occ=1,nbnd
       do ibnd_virt=1,(nbnd_total-nbnd)
        write(158,*) F(ibnd_occ,ibnd_virt,pol_index)
       enddo
    enddo
    !
    close(158)
  endif
  !
#ifdef __PARA
  end if
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
       call diropn ( iunrestart, 'restart_lanczos.'//trim(int_to_char(LR_polarization)), nwordrestart, exst)
       !
       call davcio(evc1(:,:,:,1),nwordrestart,iunrestart,1,1)
       call davcio(evc1(:,:,:,2),nwordrestart,iunrestart,2,1)
       call davcio(evc1_new(:,:,:,1),nwordrestart,iunrestart,3,1)
       call davcio(evc1_new(:,:,:,2),nwordrestart,iunrestart,4,1)
       !
       close( unit = iunrestart)
    !
    ! Writing charge response density for restart
    !
       if (charge_response == 1 ) then 
        if (resonance_condition) then
         call diropn ( iunrestart, 'restart_lanczos-rho_tot.'//trim(int_to_char(LR_polarization)), 2*nrxx*nspin_mag, exst)
         call davcio(rho_1_tot_im(:,:),2*nrxx*nspin_mag,iunrestart,1,1)
         close( unit = iunrestart)
        else
         call diropn ( iunrestart, 'restart_lanczos-rho_tot.'//trim(int_to_char(LR_polarization)), 2*nrxx*nspin_mag, exst)
         call davcio(rho_1_tot(:,:),2*nrxx*nspin_mag,iunrestart,1,1)
         close( unit = iunrestart)
        endif
       endif
       !
  return
end subroutine lr_write_restart
!-----------------------------------------------------------------------
