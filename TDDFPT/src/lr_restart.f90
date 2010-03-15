!-----------------------------------------------------------------------
subroutine lr_restart(iter_restart,rflag)
  !---------------------------------------------------------------------
  ! ... restart the Lanczos recursion
  !---------------------------------------------------------------------
  !
  ! OBM :
  ! 050608 Modified for calbec interface in v4.0
  !        gamma_only correction
#include "f_defs.h"
  !
  use io_global,            only : stdout, ionode_id
  use control_flags,                only : gamma_only
  use klist,                only : nks, xk 
  use cell_base,            only : tpiba2
  use gvect,                only : g
  use io_files,             only : tmp_dir, prefix
  use lr_variables,         only : itermax,evc1, evc1_new, sevc1_new, rho_1_tot ,&
                                   restart, nwordrestart, iunrestart,project,nbnd_total,F
  use wvfct,                only : npw, igk, nbnd, g2kin, npwx
  use lr_variables,         only : beta_store, gamma_store, zeta_store, norm0!,real_space
  use becmod,               only : bec_type, becp, calbec
  use uspp,                 only : vkb, nkb, okvan
  USE io_global,            ONLY : ionode
  use mp,                   only : mp_bcast
  !use real_beta,            only : ccalbecr_gamma,s_psir,fft_orbital_gamma,bfft_orbital_gamma
  USE realus,               ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                                    bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                    v_loc_psir, s_psir_gamma,igk_k,npw_k, &
                                    real_space_debug 
  use gvect,                only : nrxx
  USE lr_variables,         ONLY : lr_verbosity, charge_response, LR_polarization, n_ipol
  USE noncollin_module,     ONLY : nspin_mag


  !
  implicit none
  !
  character(len=6), external :: int_to_char
  !
  !integer, intent(in) :: pol
  integer, intent(out) :: iter_restart
  logical, intent(out) :: rflag
  !
  ! local variables
  !
  integer :: i,ibnd,ibnd_occ,ibnd_virt,temp
  integer :: ik, ig, ip
  logical :: exst
  character(len=256) :: tempfile, filename
  integer :: pol_index
  !
  !
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_restart>")')
  endif

  pol_index=1
  if ( n_ipol /= 1 ) pol_index=LR_polarization

  if (.not.restart) return
  !
  rflag = .false.
  !
  ! Restarting kintic-energy and ultrasoft
  !
  if (gamma_only) then
     !
     do ig=1,npwx
        !
        g2kin(ig)=((xk(1,1)+g(1,igk_k(ig,1)))**2 &
                 +(xk(2,1)+g(2,igk_k(ig,1)))**2 &
                 +(xk(3,1)+g(3,igk_k(ig,1)))**2)*tpiba2
        !
     enddo
     !
     call init_us_2(npw,igk,xk(1,1),vkb)
     !
  end if
  !
  ! Reading Lanczos coefficients
  !
  !
  filename = trim(prefix) // ".beta_gamma_z." // trim(int_to_char(LR_polarization))
  tempfile = trim(tmp_dir) // trim(filename)
  !
  inquire (file = tempfile, exist = exst)
  !
  if (.not.exst) then
     !
     WRITE( stdout,*) "WARNING: " // trim(filename) // " does not exist"
     rflag = .true.
     return
     !
  end if
  !
  !
  !Ionode only reads
  !
#ifdef __PARA
  if (ionode) then
#endif
  !
  ! Read and broadcast beta gamma zeta 
  !
  open (158, file = tempfile, form = 'formatted', status = 'old')
  !
  read(158,*,end=301,err=303) iter_restart
  !
  if ( iter_restart .ge. itermax ) iter_restart = itermax
  !
  read(158,*,end=301,err=303) norm0(pol_index)
  !
  do i=1,iter_restart
     !
     read(158,*,end=301,err=303) beta_store(pol_index,i)
     read(158,*,end=301,err=303) gamma_store(pol_index,i)
     read(158,*,end=301,err=303) zeta_store (pol_index,:,i)
     !
  end do
  !
  close(158)
#ifdef __PARA
  call mp_bcast (iter_restart, ionode_id)
  call mp_bcast (norm0(pol_index), ionode_id)
  call mp_bcast (beta_store(pol_index,:), ionode_id)
  call mp_bcast (gamma_store(pol_index,:), ionode_id)
  call mp_bcast (zeta_store(pol_index,:,:), ionode_id)
#endif
  !
  !
  ! Read projection
  !
  if (project) then
    filename = trim(prefix) // ".projection." // trim(int_to_char(LR_polarization))
    tempfile = trim(tmp_dir) // trim(filename)
    !
    !
    open (158, file = tempfile, form = 'formatted', status = 'unknown')
    !
    read(158,*,end=301,err=303) temp
    !
    if (temp /= iter_restart) call errore ('lr_restart', 'Iteration mismatch reading projections', 1 )
    !
    read(158,*,end=301,err=303) temp   !number of filled bands
    !
    if (temp /= nbnd) call errore ('lr_restart', 'NBND mismatch reading projections', 1 )
    !
    read(158,*,end=301,err=303) temp !total number of bands
    !
    if (temp /= nbnd_total) call errore ('lr_restart', 'Total number of bands mismatch reading projections', 1 )
    !
    do ibnd_occ=1,nbnd
       do ibnd_virt=1,(nbnd_total-nbnd)
        read(158,*,end=301,err=303) F(ibnd_occ,ibnd_virt,pol_index)
       enddo
    enddo
    !
    close(158)
  endif
#ifdef __PARA
  call mp_bcast (F, ionode_id)
#endif


#ifdef __PARA
  end if
#endif
  !
  iter_restart = iter_restart + 1
  !
  ! Reading Lanczos vectors
  !
  nwordrestart = 2 * nbnd * npwx * nks
  !
  call diropn ( iunrestart, 'restart_lanczos.'//trim(int_to_char(LR_polarization)), nwordrestart, exst)
  !
  call davcio(evc1(:,:,:,1),nwordrestart,iunrestart,1,-1)
  call davcio(evc1(:,:,:,2),nwordrestart,iunrestart,2,-1)
  call davcio(evc1_new(:,:,:,1),nwordrestart,iunrestart,3,-1)
  call davcio(evc1_new(:,:,:,2),nwordrestart,iunrestart,4,-1)
  !
  close( unit = iunrestart)
  if (charge_response == 2 ) then 
         call diropn ( iunrestart, 'restart_lanczos-rho_tot.'//trim(int_to_char(LR_polarization)), 2*nrxx, exst)
         call davcio(rho_1_tot(:,:),2*nrxx*nspin_mag,iunrestart,1,-1)
         close( unit = iunrestart)
       endif


  !
  ! Reinitializing sevc1_new vector
  !
  if (gamma_only) then
     !
     if ( nkb > 0 .and. okvan ) then 
      if (real_space_debug>6) then
       do ibnd=1,nbnd,2
        call fft_orbital_gamma(evc1_new(:,:,1,1),ibnd,nbnd)
        call calbec_rs_gamma(ibnd,nbnd,becp%r)
        call s_psir_gamma(ibnd,nbnd)
        call bfft_orbital_gamma(sevc1_new(:,:,1,1),ibnd,nbnd)
       enddo
      else
       call calbec(npw_k(1),vkb,evc1_new(:,:,1,1),becp)
       !call pw_gemm('Y',nkb,nbnd,npw_k(1),vkb,npwx,evc1_new(1,1,1,1),npwx,rbecp,nkb)
       call s_psi(npwx,npw_k(1),nbnd,evc1_new(:,:,1,1),sevc1_new(:,:,1,1))
      endif
     else
      !nkb = 0 not real space
       !
       call s_psi(npwx,npw_k(1),nbnd,evc1_new(:,:,1,1),sevc1_new(:,:,1,1))
       !
      !
     endif
     !
     if ( nkb > 0 .and. okvan ) then 
      if (real_space_debug>6) then  
        do ibnd=1,nbnd,2
        call fft_orbital_gamma(evc1_new(:,:,1,2),ibnd,nbnd)
        call calbec_rs_gamma(ibnd,nbnd,becp%r)
        call s_psir_gamma(ibnd,nbnd)
        call bfft_orbital_gamma(sevc1_new(:,:,1,2),ibnd,nbnd)
       enddo
     else
       call calbec(npw_k(1),vkb,evc1_new(:,:,1,2),becp%r)
      !call pw_gemm('Y',nkb,nbnd,npw_k(1),vkb,npwx,evc1_new(1,1,1,2),npwx,rbecp,nkb)
       call s_psi(npwx,npw_k(1),nbnd,evc1_new(:,:,1,2),sevc1_new(:,:,1,2))
      endif
     endif
     !call s_psi(npwx,npw_k(1),nbnd,evc1_new(:,:,1,2),sevc1_new(:,:,1,2))
     !
  else
     !
     do ik=1,nks
        !
        if ( nkb > 0 .and. okvan ) then
           call init_us_2(npw_k(ik),igk_k(1,ik),xk(1,ik),vkb)
           !call ccalbec(nkb,npwx,npw_k(ik),nbnd,becp,vkb,evc1_new(1,1,ik,1))
           call calbec(npw_k(ik), vkb, evc1_new(:,:,ik,1), becp)
        endif
        call s_psi(npwx,npw_k(ik),nbnd,evc1_new(:,:,ik,1),sevc1_new(:,:,ik,1))
        !
        if (nkb > 0) call calbec(npw_k(ik), vkb, evc1_new(:,:,ik,2), becp)
        !call ccalbec(nkb,npwx,npw_k(ik),nbnd,becp,vkb,evc1_new(1,1,ik,2))
        call s_psi(npwx,npw_k(ik),nbnd,evc1_new(:,:,ik,2),sevc1_new(:,:,ik,2))
        !
     enddo
     !
  end if
  !
  !
  return
  301 call errore ('restart', 'A File is corrupted, file ended unexpectedly', 1 ) 
  303 call errore ('restart', 'A File is corrupted, error in reading data', 1)
end subroutine lr_restart
!-----------------------------------------------------------------------
