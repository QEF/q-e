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

  IF ( .NOT. d0psi_rs ) THEN
     !
     ! G-space: compute commutator [H,x]
     !
     ! Note: If you are performing a hybrid calculation ( e.g. PBE0, B3LYP),
     ! then this option is not full because the commutator [V_EXX,r] 
     ! is not implemented. 
     !
     WRITE(stdout,'(/,5X,"Calculation of the dipole in reciprocal space")')
     !
     DO ik = 1, nks
        !
        current_k = ik
        !
        IF ( lsda ) current_spin = isk(ik)
        ! 
        evc(:,:) = evc0(:,:,ik)
        !
        npw = npw_k(ik)
        igk(:) = igk_k(:,ik)
        !
        CALL init_us_2(npw,igk,xk(1,ik),vkb)
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
     IF (gstart == 2 .and. gamma_only) d0psi(1,:,:,:) = cmplx(dble(d0psi(1,:,:,:)),0.0d0,dp)
     !
  ELSE
     !
     ! R-space
     !
     WRITE(stdout,'(/,5X,"Calculation of the dipole in real space")')
     !
     CALL compute_d0psi_rs(n_ipol)
     !
  ENDIF
  !
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
     IF (n_ipol==1) CALL diropn ( iund0psi, 'd0psi.'//trim(int_to_char(LR_polarization)), nwordd0psi, exst)
     IF (n_ipol==3) CALL diropn ( iund0psi, 'd0psi.'//trim(int_to_char(ip)), nwordd0psi, exst)
     CALL davcio(d0psi(1,1,1,ip),nwordd0psi,iund0psi,1,1)
     CLOSE( UNIT = iund0psi)
  ENDDO
  ! End of file i/o
  tmp_dir = tmp_dir_saved

  CALL stop_clock ('lr_solve_e')
  !WRITE(stdout,'(5X,"lr_wfcinit_spectrum: finished lr_solve_e")')
  RETURN

contains
!--------------------------------------------------------------------
SUBROUTINE compute_d0psi_rs(  n_ipol )
!--------------------------------------------------------------------
!
! ... original code from routine compute_dipole
! ... modified to calculate the d0psi in the real space
! ... by Xiaochuan Ge, Oct, 2013
! ... modified by I. Timrov, June 2014
!
! ... Warning! If you are using this routine, the molecule
! ... must be in the center of the supercell (not in the corner)!
! ... TODO: improve this routine using the logic of compute_dipole
!
  USE kinds,            ONLY : DP
  USE cell_base,        ONLY : at, bg, alat, omega, tpiba
  USE fft_base,         ONLY : dfftp,dffts
  USE mp_global,        ONLY : me_bgrp, intra_bgrp_comm
  use mp,               only : mp_barrier
  use io_global,        only : stdout
  USE wvfct,            ONLY : nbnd,npwx,current_k,igk,npw, g2kin
  use klist,            only : xk,nks
  use lr_variables,     only : evc0,sevc0,d0psi, lshift_cell
  use wavefunctions_module, only : psic, evc
  use uspp,             only : okvan
  USE lsda_mod,         ONLY : lsda, isk
  USE uspp,             ONLY : vkb
  USE control_flags,    ONLY : gamma_only
  USE becmod,           ONLY : becp, calbec
  USE gvect,            ONLY : g
  use realus,           only :  bfft_orbital_gamma,fft_orbital_gamma, &
                              & fft_orbital_k, bfft_orbital_k, &
                              & igk_k,npw_k
  !
  IMPLICIT NONE
  ! ... Local variables
  REAL(DP),allocatable :: r(:,:), gk(:,:)
  INTEGER  :: ik, i, j, k, ip, ir, ig, ir_end, index0,ibnd,n_ipol
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3

  if (.not.allocated(psic)) allocate(psic(dfftp%nnr))
  ! 
  allocate(r(dfftp%nnr,n_ipol))
  !
  if (okvan) then
      write(stdout,'(10x,"At this moment d0psi_rs is not available for USPP !!!",//)')
#ifdef __MPI
      call mp_barrier(intra_bgrp_comm)
#endif
      stop
  endif
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
  DO ik = 1, nks
     !
     current_k = ik
     !
     IF ( lsda ) current_spin = isk(ik)
     !
     evc(:,:) = evc0(:,:,ik)
     !
     npw = npw_k(ik)
     igk(:) = igk_k(:,ik)
     !
     CALL init_us_2(npw,igk,xk(1,ik),vkb)
     !
     ! Compute the kinetic energy which is needed for the whole code.
     !
     allocate (gk(3,npwx))
     do ig = 1, npw
        gk(1:3,ig) = (xk(1:3,ik) + g(1:3,igk(ig))) * tpiba
        g2kin(ig) = SUM(gk(1:3,ig)**2 )
     enddo
     deallocate(gk)
     !
     if(lshift_cell) call shift_cell(r,n_ipol)
     IF (gamma_only) THEN
        !
        DO ibnd = 1,nbnd,2 
           !
           ! FFT to R-space: evc0 -> psic
           !
           CALL fft_orbital_gamma(evc0(:,:,1), ibnd, nbnd)
           !
           DO ip = 1, n_ipol
              !
              ! Apply the dipole operator
              !
              DO ir = 1, dfftp%nnr
                 psic(ir) = r(ir,ip)*alat*psic(ir)
              ENDDO
              !
              ! FFT to G-space: psic -> d0psi
              !
              CALL bfft_orbital_gamma(d0psi(:,:,1,ip), ibnd, nbnd)
              !
           ENDDO
           !
        ENDDO
        !
     ELSE
        !
        ! k-points algorithm
        !
        DO ibnd = 1,nbnd
           !
           ! FFT to R-space: evc0 -> psic
           !
           CALL fft_orbital_k(evc0(:,:,ik),ibnd,nbnd,ik)
           !
           DO ip = 1, n_ipol
              !
              ! Apply dipole operator
              !
              DO ir = 1, dfftp%nnr
                 psic(ir) = r(ir,ip)*alat*psic(ir)
              ENDDO
              !
              ! FFT to G-space: psic -> d0psi
              !
              CALL bfft_orbital_k(d0psi(:,:,ik,ip),ibnd,nbnd,ik)
              !
           ENDDO
           !
        ENDDO
        !
     ENDIF
     !
   ENDDO ! k points
   !
   ! Orthogonalized batch orbitals to occupied minifold
   !
   do ik=1,nks
    do ip = 1, n_ipol
      CALL lr_ortho(d0psi(:,:,ik,ip), evc0(:,:,ik), ik, ik, sevc0(:,:,ik),.true.)
    enddo
   enddo
   !
   deallocate(r)
   !
   RETURN
   !
 END SUBROUTINE compute_d0psi_rs
!----------------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE shift_cell( r, n_ipol )
!--------------------------------------------------------------------
  ! Shift the molecule to the center of the cell for a more proper
  ! calculation of d0psi in r-space
  USE fft_base,         ONLY : dfftp
  use kinds,            only : dp
  use ions_base,        only : nat, tau
  USE io_global,            ONLY : stdout

  REAL(dp),intent(inout) :: r(dfftp%nnr,n_ipol)
  integer,intent(in) :: n_ipol
  real(dp) ::  mmin(3), mmax(3), center(3), origin(3)
  integer :: ip, iatm, ir

  WRITE(stdout,'(/,5X,"Wavefunctions are shifted to the center of cell for the calculation of d0psi. ")')

  mmin(:) = 2.d0
  mmax(:)= -1.d0

  do ip = 1, n_ipol
    do iatm = 1, nat
      mmin(ip) = min(mmin(ip), tau(ip,iatm))
      mmax(ip) = max(mmax(ip), tau(ip,iatm))
    enddo
  enddo

  center(:)= 0.5d0*(mmin(:)+mmax(:))
  origin(:)= center(:)-0.5d0
  
  do ir = 1, dfftp%nnr
    r(ir,:)= r(ir,:) - origin(:)
    do ip = 1, n_ipol
      if(r(ir,ip) .lt. 0) r(ir,ip)=r(ir,ip)+1.d0
      if(r(ir,ip) .gt. 1) r(ir,ip)=r(ir,ip)-1.d0
    enddo
  enddo

  return
  END SUBROUTINE shift_cell
!----------------------------------------------------------------------------


END SUBROUTINE lr_solve_e
