!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE lr_solve_e
  !-----------------------------------------------------------------------
  !
  !  This subroutine is a driver for the solution of the linear
  !  system which defines the change of the wavefunction due to
  !  an external perturbation.
  !  Calculates the initial starting vectors for use in the Lanczos
  !  algorithm.
  !  Inspired by PHonon subroutine "solve_e".
  !
  !  Written by Brent Walker (2004)
  !  Modified by Osman Baris Malcioglu (2009)
  !  Modified by Iurii Timrov (2013) 
  !
  USE kinds,                ONLY : dp
  USE gvect,                ONLY : gstart
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : diropn, tmp_dir, wfc_dir
  USE klist,                ONLY : nks, xk, ngk, igk_k, degauss
  USE lr_variables,         ONLY : nwordd0psi, iund0psi,LR_polarization, test_case_no, &
                                   & n_ipol, evc0, d0psi, d0psi2, evc1, lr_verbosity, &
                                   & d0psi_rs, eels
  USE lsda_mod,             ONLY : lsda, isk, current_spin
  USE uspp,                 ONLY : vkb
  USE wvfct,                ONLY : nbnd, npwx, et, current_k
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp,                   ONLY : mp_max, mp_min, mp_barrier
  USE realus,               ONLY : real_space, real_space_debug 
  USE control_lr,           ONLY : alpha_pv
  USE qpoint,               ONLY : nksq
  !
  IMPLICIT NONE
  INTEGER :: ibnd, ik, is, ip
  ! counter on bands
  ! counter on k points
  ! counter on spins
  ! counter on polarizations
  CHARACTER(len=6), EXTERNAL :: int_to_char
  LOGICAL :: exst
  real (kind=dp) :: anorm
  CHARACTER(len=256) :: tmp_dir_saved
  !
  IF (lr_verbosity > 5) WRITE(stdout,'("<lr_solve_e>")')
  !
  CALL start_clock ('lr_solve_e')
  !
  IF( lr_verbosity > 1 ) &
       WRITE(stdout,'(5X,"lr_solve_e: alpha_pv=",1X,e12.5)') alpha_pv
  !
  IF (eels) THEN
     !
     ! EELS
     !
     DO ik = 1, nksq
        CALL lr_dvpsi_eels(ik, d0psi(:,:,ik,1), d0psi2(:,:,ik,1))
     ENDDO
     !
  ELSE
     !
     ! Optical case
     !
     DO ik = 1, nks
        !
        current_k = ik
        !
        IF ( lsda ) current_spin = isk(ik)
        !
        evc(:,:) = evc0(:,:,ik)
        !
        ! Ultrasoft case: calculate beta-functions vkb.
        !
        CALL init_us_2(ngk(ik), igk_k(:,ik), xk(:,ik), vkb)
        !
        !   Computes/reads P_c^+ x psi_kpoint into d0psi array
        !
        IF ( n_ipol==3 ) THEN
           DO ip=1,3
              CALL lr_dvpsi_e(ik,ip,d0psi(:,:,ik,ip))
           ENDDO
        ELSEIF ( n_ipol==1 ) THEN
           CALL lr_dvpsi_e(ik,LR_polarization,d0psi(:,:,ik,1))
        ENDIF
        !
     ENDDO
     !
  ENDIF
  !
  IF (gstart == 2 .and. gamma_only) d0psi(1,:,:,:) = cmplx(dble(d0psi(1,:,:,:)),0.0d0,dp)
  !
  ! X. Ge: compute d0psi in real-space, will overwrite d0psi calculated before.
  !
  IF (.NOT. eels) THEN
     IF (d0psi_rs) CALL compute_d0psi_rs(n_ipol)
  ENDIF
  !
  ! Writing of d0psi to the file. 
  ! This is a parallel writing, done in wfc_dir 
  !
  tmp_dir_saved = tmp_dir
  !
  IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir
  !
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
  !
  ! EELS: Writing of d0psi2 to the file.
  !
  IF (eels) THEN
     !
     CALL diropn ( iund0psi, 'd0psi2.'//trim(int_to_char(LR_polarization)), nwordd0psi, exst)
     CALL davcio(d0psi2(1,1,1,1),nwordd0psi,iund0psi,1,1)
     CLOSE( UNIT = iund0psi)
     !
  ENDIF
  !
  tmp_dir = tmp_dir_saved
  !
  CALL stop_clock ('lr_solve_e')
  !
  RETURN
  !
CONTAINS

SUBROUTINE compute_d0psi_rs( n_ipol )
  !
  ! ... original code from routine compute_dipole
  ! ... modified to calculate the d0psi in the real space
  ! ... by Xiaochuan Ge, Oct, 2013
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at, bg, alat, omega
  USE fft_base,             ONLY : dfftp,dffts
  USE mp_global,            ONLY : me_bgrp, intra_bgrp_comm
  USE mp,                   ONLY : mp_barrier
  USE io_global,            ONLY : stdout
  USE wvfct,                ONLY : nbnd,npwx
  USE klist,                ONLY : nks
  USE lr_variables,         ONLY : evc0,sevc0,d0psi,lshift_d0psi
  USE wavefunctions_module, ONLY : psic
  USE uspp,                 ONLY : okvan
  USE realus,               ONLY : fwfft_orbital_gamma,invfft_orbital_gamma
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
  !
  IF (.not. ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
  !
  ALLOCATE(psic_temp(dfftp%nnr))
  ALLOCATE(r(dfftp%nnr,n_ipol))
  !
  ! ... Initialization
  WRITE(stdout,'(/,5x,"Calculating d0psi in the real space.")')
  !
  IF (okvan) THEN
     !
     WRITE(stdout,'(10x,"At this moment d0psi_rs is not available for USPP !!!",//)')
#ifdef __MPI
      call mp_barrier(intra_bgrp_comm)
!     call mp_stop(100)
#endif
      stop
  endif
  !
  IF (nks > 1) CALL errore( 'compute_d0psi_rs', 'nks>1 is not supported', 1 )
  !
  ! Calculate r
  !
  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
  !
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
     !
     i = index0 + ir - 1
     k = i / (dfftp%nr1x*dfftp%nr2x)
     i = i - (dfftp%nr1x*dfftp%nr2x)*k
     j = i / dfftp%nr1x
     i = i - dfftp%nr1x*j
     !
     DO ip = 1, n_ipol
        r(ir,ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                   DBLE( j )*inv_nr2*at(ip,2) + &
                   DBLE( k )*inv_nr3*at(ip,3)
     ENDDO
     !
  ENDDO
  !
  If (lshift_d0psi) call shift_d0psi(r,n_ipol)
  !
  DO ib = 1, nbnd 
     !
     wfck(:,1) = evc0(:,ib,1)
     !
     CALL invfft_orbital_gamma(wfck(:,:),1,1)
     !
     psic_temp(:) = psic(:)
     !
     DO ip = 1, n_ipol
       !
       ! Apply the dipole operator
       !
       DO ir = 1, dfftp%nnr
          psic(ir) = r(ir,ip) * alat * psic_temp(ir)
       ENDDO
       !
       ! Convert to G-space
       !
       CALL fwfft_orbital_gamma(wfck(:,:),1,1) 
       !
       d0psi(:,ib,1,ip) = wfck(:,1)
       !
     ENDDO
     !
  ENDDO
  ! 
  ! Orthogonalized batch orbitals to occupied minifold
  ! 
  DO ip = 1, n_ipol
     CALL orthogonalize(d0psi(:,:,1,ip), evc0(:,:,1), 1, 1, sevc0(:,:,1), ngk(1), .true.)
     d0psi(:,:,1,ip) = -d0psi(:,:,1,ip)
  ENDDO
  !
  DEALLOCATE(r)
  DEALLOCATE(psic_temp)
  !
  RETURN
  !
END SUBROUTINE compute_d0psi_rs

SUBROUTINE shift_d0psi( r, n_ipol )
  !
  ! Shift a position operator r to the center of the molecule 
  ! for a proper calculation of d0psi in R-space.
  !
  USE fft_base,         ONLY : dfftp
  use kinds,            only : dp
  use ions_base,        only : nat, tau
  USE io_global,        ONLY : stdout
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

  RETURN

END SUBROUTINE shift_d0psi

END SUBROUTINE lr_solve_e
