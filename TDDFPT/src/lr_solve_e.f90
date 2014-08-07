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
  !
  USE kinds,                ONLY : dp
  USE gvect,                ONLY : gstart
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : diropn, tmp_dir, wfc_dir
  USE klist,                ONLY : nks, xk, degauss
  USE lr_variables,         ONLY : nwordd0psi, iund0psi,LR_polarization, test_case_no
  USE lr_variables,         ONLY : n_ipol, evc0, d0psi, evc1,lr_verbosity,d0psi_rs
  USE realus,               ONLY : igk_k,npw_k
  USE lsda_mod,             ONLY : lsda, isk, current_spin
  USE uspp,                 ONLY : vkb
  USE wvfct,                ONLY : igk, nbnd, npwx, npw, et, current_k
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp,                   ONLY : mp_max,mp_min,mp_barrier
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
  !
  CHARACTER(len=6), EXTERNAL :: int_to_char
  LOGICAL :: exst
  real (kind=dp) :: anorm
  CHARACTER(len=256) :: tmp_dir_saved
  !
  IF (lr_verbosity > 5) WRITE(stdout,'("<lr_solve_e>")')
  !
  CALL start_clock ('lr_solve_e')

  IF( lr_verbosity > 1 ) &
       WRITE(stdout,'(5X,"lr_solve_e: alpha_pv=",1X,e12.5)') alpha_pv

  DO ik=1,nks
       current_k=ik
       IF ( lsda ) current_spin = isk(ik)
       evc(:,:)=evc0(:,:,ik)

       npw=npw_k(ik)
       igk(:)=igk_k(:,ik)

       CALL init_us_2(npw,igk,xk(1,ik),vkb)

       !   Computes/reads P_c^+ x psi_kpoint into d0psi array
       IF ( n_ipol==3 ) THEN
         DO ip=1,3
           CALL lr_dvpsi_e(ik,ip,d0psi(:,:,ik,ip))
         ENDDO
       ELSEIF ( n_ipol==1 ) THEN
         CALL lr_dvpsi_e(ik,LR_polarization,d0psi(:,:,ik,1))
       ENDIF
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

  ! Writing d0psi for restart
  nwordd0psi = 2 * nbnd * npwx * nks

  !   Reading of files:
  !   This is a parallel read, done in wfc_dir
  tmp_dir_saved = tmp_dir
  IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir

  DO ip = 1, n_ipol
     IF (n_ipol==1) CALL diropn ( iund0psi, 'd0psi.'//trim(int_to_char(LR_polarization)), nwordd0psi, exst)
     IF (n_ipol==3) CALL diropn ( iund0psi, 'd0psi.'//trim(int_to_char(ip)), nwordd0psi, exst)
     CALL davcio(d0psi(1,1,1,ip),nwordd0psi,iund0psi,1,1)
     CLOSE( UNIT = iund0psi)
  ENDDO
  ! End of file i/o
  tmp_dir = tmp_dir_saved
  
  if(d0psi_rs .eqv. .true.) call compute_d0psi_rs(n_ipol)

  CALL stop_clock ('lr_solve_e')
  WRITE(stdout,'(5X,"lr_wfcinit_spectrum: finished lr_solve_e")')
  RETURN

contains
!--------------------------------------------------------------------
SUBROUTINE compute_d0psi_rs(  n_ipol )
!--------------------------------------------------------------------
!
! ... original code from routine compute_dipole
! ... modified to calculate the d0psi in the real space
! ... by Xiaochuan Ge, Oct, 2013
!
  USE kinds,            ONLY : DP
  USE cell_base,        ONLY : at, bg, alat, omega
  USE fft_base,         ONLY : dfftp,dffts
  USE mp_global,        ONLY : me_bgrp, intra_bgrp_comm
  use mp,               only : mp_barrier
  use io_global,    only : stdout
  USE wvfct,            ONLY : nbnd,npwx
  use klist,            only : nks
  use lr_variables,     only : evc0,sevc0,d0psi,lshift_d0psi
  use wavefunctions_module, only : psic
  use uspp,           only : okvan
  use realus,          only :  bfft_orbital_gamma,fft_orbital_gamma
  !
  IMPLICIT NONE
  !
  ! ... Define variables
  complex(dp) :: wfck(npwx,1)

  ! ... Local variables
  REAL(DP),allocatable :: r(:,:)
  complex(dp),allocatable :: psic_temp(:)
  INTEGER  :: i, j, k, ip, ir, ir_end, index0,ib,n_ipol
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3

  if(.not. allocated(psic)) allocate(psic(dfftp%nnr))
  allocate(psic_temp(dfftp%nnr))
  allocate(r(dfftp%nnr,n_ipol))

  ! ... Initialization
  write(stdout,'(5x,"Calculating d0psi in the real space."//)')
  if(okvan) then
      write(stdout,'(10x,"At this moment d0psi_rs is not available for USPP !!!",//)')
#ifdef __MPI
      call mp_barrier(intra_bgrp_comm)
 !     call mp_stop(100)
#endif
      stop
  endif

  ! Calculat r
  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )

#if defined (__MPI)
  index0 = dfftp%nr1x*dfftp%nr2x*SUM(dfftp%npp(1:me_bgrp))
  ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npp(me_bgrp+1))
#else
  index0 = 0
  ir_end = dfftp%nnr
#endif
  !
  DO ir = 1, ir_end 
     !
     ! ... three dimensional indexes
     i = index0 + ir - 1
     k = i / (dfftp%nr1x*dfftp%nr2x)
     i = i - (dfftp%nr1x*dfftp%nr2x)*k
     j = i / dfftp%nr1x
     i = i - dfftp%nr1x*j

     DO ip = 1, n_ipol
        r(ir,ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                DBLE( j )*inv_nr2*at(ip,2) + &
                DBLE( k )*inv_nr3*at(ip,3)
     END DO
  enddo

  if(lshift_d0psi) call shift_d0psi(r,n_ipol)

   do ib = 1, nbnd 
     wfck(:,1)=evc0(:,ib,1)
     call fft_orbital_gamma(wfck(:,:),1,1)
     psic_temp(:)=psic(:)
     DO ip = 1, n_ipol
       ! Apply dipole operator
       do ir = 1, dfftp%nnr
         psic(ir)=r(ir,ip)*alat*psic_temp(ir)
       enddo
       ! Convert to g space
       call bfft_orbital_gamma(wfck(:,:),1,1) 
       d0psi(:,ib,1,ip)=wfck(:,1)
     enddo
   enddo
  
   ! orthogonalized batch orbitals to occupied minifold
   do ip = 1, n_ipol
     CALL lr_ortho(d0psi(:,:,:,ip), evc0(:,:,1), 1, 1, sevc0(:,:,1),.true.)
   enddo

  deallocate(r)
  deallocate(psic_temp)
  RETURN

  END SUBROUTINE compute_d0psi_rs
!----------------------------------------------------------------------------


!--------------------------------------------------------------------
  SUBROUTINE shift_d0psi( r, n_ipol )
!--------------------------------------------------------------------
  ! Shift the molecule to the center of the cell for a more proper
  ! calculation of d0psi in r-space
  USE fft_base,         ONLY : dfftp
  use kinds,            only : dp
  use ions_base,        only : nat, tau
  USE io_global,            ONLY : stdout
  use cell_base,        only : alat, at

  REAL(dp),intent(inout) :: r(dfftp%nnr,n_ipol)
  integer,intent(in) :: n_ipol
  real(dp) ::  mmin(3), mmax(3), center(3), origin(3), check_cell
  integer :: ip, iatm, ir, ip2

  WRITE(stdout,'(/,5X,"Wavefunctions are shifted to the center of cell for the calculation of d0psi. ")')

  check_cell = 0.d0
  do ip = 1, 3
    do ip2 = 1,3
      if(.not. ip .eq. ip2) check_cell = check_cell + at(ip,ip2)**2
    enddo
  enddo
  if(check_cell .gt. 1.d-5) call errore('lr_read_wfc'," I am not sure that this type of super cell is supported now. ",1)

  mmin(:) = 2000.d0
  mmax(:)= -2000.d0

  do ip = 1, n_ipol
    do iatm = 1, nat
      mmin(ip) = min(mmin(ip), tau(ip,iatm))
      mmax(ip) = max(mmax(ip), tau(ip,iatm))
    enddo
  enddo

  center(:)= 0.5d0*(mmin(:)+mmax(:))
  do ip = 1, n_ipol
    origin(ip)= center(ip)-0.5d0*at(ip,ip)
  enddo

  do ir = 1, dfftp%nnr
    r(ir,:)= r(ir,:) - origin(:)
    do ip = 1, n_ipol
      if(r(ir,ip) .lt. 0) r(ir,ip)=r(ir,ip)+at(ip,ip)
      if(r(ir,ip) .gt. at(ip,ip)) r(ir,ip)=r(ir,ip)-at(ip,ip)
    enddo
  enddo

  return
  END SUBROUTINE shift_d0psi
!----------------------------------------------------------------------------

END SUBROUTINE lr_solve_e
