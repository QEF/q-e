!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE lr_solve_e
  !-----------------------------------------------------------------------
  !
  ! bwalker:   This routine is a driver for the solution of the linear
  ! bwalker:   system which defines the change of the wavefunction
  ! bwalker:   due to an electric field.
  ! bwalker:   Calculates the initial starting vectors for use in the
  ! bwalker:   block Lanczos.
  ! bwalker:   We have to solve to find the action of the electric field
  ! bwalker:   operator on the initial state.
  ! bwalker:   Inspired by PHONON subroutine "solve_e".
  !
  !-----------------------------------------------------------------------
  !
  ! Modified by Osman Baris Malcioglu (2009)
#include "f_defs.h"
  USE kinds,                ONLY : dp
  USE gvect,                ONLY : gstart
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : diropn, tmp_dir, wfc_dir
  USE klist,                ONLY : nks, xk, degauss
  USE lr_variables,         ONLY : nwordd0psi, iund0psi,LR_polarization, test_case_no
  USE lr_variables,         ONLY : n_ipol, evc0, d0psi, evc1, lr_verbosity
  USE realus,               ONLY : igk_k,npw_k
  USE lsda_mod,             ONLY : lsda, isk, current_spin
  USE uspp,                 ONLY : vkb
  USE wvfct,                ONLY : igk, nbnd, npwx, npw, et
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_max,mp_min
  USE realus,               ONLY : real_space, real_space_debug!, dvpsir_e
  USE control_ph,           ONLY : alpha_pv
  !
  IMPLICIT NONE
  !
  ! counter on bands
  ! counter on k points
  ! counter on spins
  ! counter on polarizations
  INTEGER :: ibnd, ik, is, ip
  !
  !OBM!! this has been moved to lr_init_nfo
  !! variables for calculating lr_alpha_pv
  !real(kind=dp) :: emin, emax
  !
  CHARACTER(len=6), EXTERNAL :: int_to_char
  LOGICAL :: exst
  real (kind=dp) :: anorm
  CHARACTER(len=256) :: tmp_dir_saved
  !
  IF (lr_verbosity > 5) WRITE(stdout,'("<lr_solve_e>")')
  !if ( lsda ) call errore ( 'lr_solve_e' , ' LSDA not implemented' , 1)
  !
  CALL start_clock ('lr_solve_e')

  !OBM!!! This has been moved to lr_init_nfo
!  !!
!  !!   Calculate spread in eigenvalues, corresponds to alpha_pv in PHONON
!  !!
!  !if (degauss /= 0) then
!     !
!     call errore(' lr_solve_e ','degauss not equal to 0 ',1)
!     !
!  else
!     !
!     emin=minval(et(:,:))
!     !
!#ifdef __MPI
!     !   Find the minimum across pools
!     !call poolextreme(emin,-1)
!     call mp_min(emin, inter_pool_comm)
!#endif
!     !
!     emax=maxval(et(:,:))
!     !
!#ifdef __MPI
!     !   Find the maximum across pools
!     !call poolextreme(emax,+1)
!     call mp_max(emax, inter_pool_comm)
!#endif
  !   !
  !   lr_alpha_pv=2.0d0*(emax-emin)
  !   !   Avoid zero value for alpha_pv
  !   lr_alpha_pv = max(lr_alpha_pv,1.0d-2)
  !   !
  !endif
  !!
  IF( lr_verbosity > 1 ) &
       WRITE(stdout,'(5X,"lr_solve_e: alpha_pv=",1X,e12.5)') alpha_pv
  !
  !
  !if ( real_space_debug > 8 .and. gamma_only) then
  !print *, "Experimental, non-vkb electric field operator"
  !     evc(:,:)=evc0(:,:,1)
  !     if ( n_ipol==3 ) then
  !        !
  !         do ip=1,3
  !           !
  !           call dvpsir_e(ik,ip,d0psi(:,:,1,ip),lr_alpha_pv)
  !           !
  !        end do
  !        !
  !     else if ( n_ipol==1 ) then
  !       !
  !       call dvpsir_e(ik,ipol,d0psi(:,:,1,1),lr_alpha_pv)
  !       !
  !     end if
  !else
  !  print *, "Vkb electric field operator"
    DO ik=1,nks
        !
        IF ( lsda ) current_spin = isk(ik)
       !
      evc(:,:)=evc0(:,:,ik)
      !
       npw=npw_k(ik)
       igk(:)=igk_k(:,ik)
       !
       CALL init_us_2(npw,igk,xk(1,ik),vkb)
       !
       !   Computes/reads P_c^+ x psi_kpoint into d0psi array
       !
       IF ( n_ipol==3 ) THEN
          !
           DO ip=1,3
             !
             CALL lr_dvpsi_e(ik,ip,d0psi(:,:,ik,ip))
             !
          ENDDO
          !
       ELSEIF ( n_ipol==1 ) THEN
         !
         CALL lr_dvpsi_e(ik,LR_polarization,d0psi(:,:,ik,1))
         !
       ENDIF
       !
       !print *, "lr_solve_e, after lr_dvpsi_e"
       !CALL lr_normalise( d0psi(:,:,1,1), anorm)
    ENDDO
  !endif
  !
  IF (gstart == 2 .and. gamma_only) d0psi(1,:,:,:) = cmplx(dble(d0psi(1,:,:,:)),0.0d0,dp)
!OBM!!! debug
IF (test_case_no == 2) THEN
          PRINT *,"dumping d0psi"
          OPEN(UNIT=47,FILE="d0psi.dump",STATUS='NEW',ACCESS = 'SEQUENTIAL')
          WRITE(unit=47,FMT=*) "Kpoint --- band --- plane wave --- value for pol1 --- value for pol2 --- value for pol3"
          DO ik=1,nks
           DO ibnd=1,nbnd
            DO ip=1, npw
             WRITE(unit=47,FMT='(I3," ",2(I7," "), 3("(",E14.5," ",E14.5,"i)"))') ik,  &
             ibnd, ip, d0psi(ip,ibnd,ik,1), d0psi(ip,ibnd,ik,3), d0psi(ip,ibnd,ik,3)
            ENDDO
           ENDDO
          ENDDO
          CLOSE(47)
          PRINT *, "dump complete"
ENDIF
!OBM!!! end of debug

       !print *, "lr_solve_e before dump"
       !CALL lr_normalise( d0psi(:,:,1,1), anorm)
  !
  ! Writing d0psi for restart
  !
  nwordd0psi = 2 * nbnd * npwx * nks
  !
  !   Reading of files:
  !   This is a parallel read, done in wfc_dir
  tmp_dir_saved = tmp_dir
  IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir

  DO ip = 1, n_ipol
     !
     IF (n_ipol==1) CALL diropn ( iund0psi, 'd0psi.'//trim(int_to_char(LR_polarization)), nwordd0psi, exst)
     IF (n_ipol==3) CALL diropn ( iund0psi, 'd0psi.'//trim(int_to_char(ip)), nwordd0psi, exst)
     !
     CALL davcio(d0psi(1,1,1,ip),nwordd0psi,iund0psi,1,1)
     !
     CLOSE( UNIT = iund0psi)
     !
  ENDDO
  ! End of file i/o
  tmp_dir = tmp_dir_saved
  !
  ! end writing
  !
  CALL stop_clock ('lr_solve_e')
  !
  WRITE(stdout,'(5X,"lr_wfcinit_spectrum: finished lr_solve_e")')
  !
  RETURN
  !
END SUBROUTINE lr_solve_e
!-------------------------------------------------------------------------
