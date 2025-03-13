!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
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
                                   & d0psi_rs, eels, lr_exx,  magnons, &
                                   & V0psi, ipol, O_psi, n_op
  USE lsda_mod,             ONLY : lsda, isk, current_spin,nspin
  USE uspp,                 ONLY : vkb, okvan
  USE wvfct,                ONLY : nbnd, npwx, et, current_k
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : evc
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp,                   ONLY : mp_max, mp_min, mp_barrier
  USE control_lr,           ONLY : alpha_pv, nbnd_occx
  USE qpoint,               ONLY : nksq
  USE noncollin_module,     ONLY : npol,noncolin
  USE uspp_param,           ONLY : nhm
  USE ions_base,            ONLY : nat
  USE lrus,                 ONLY : intq, intq_nc
  USE uspp_init,            ONLY : init_us_2

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
  INTEGER :: pol_index
  !
  CALL start_clock ('lr_solve_e')
  !
  IF (eels) THEN
     !
     ! EELS case
     !
     IF (okvan) THEN
        ALLOCATE (intq (nhm, nhm, nat))
        IF (noncolin) ALLOCATE(intq_nc( nhm, nhm, nat, nspin))
        CALL lr_compute_intq()
     ENDIF
     !
     DO ik = 1, nksq
        CALL lr_dvpsi_eels(ik, d0psi(:,:,ik,1), d0psi2(:,:,ik,1))
     ENDDO
     !
     IF (okvan) THEN
        DEALLOCATE (intq)
        IF (noncolin) DEALLOCATE(intq_nc)
     ENDIF
  ELSE IF (magnons) THEN
     !
     ! MAGNONS case
     !
     WRITE(stdout,'(5X,"magnon calculation, n_ipol =",1X,i3,1x,"n_op =", 1X,i3)') n_ipol, n_op
     !
     V0psi = (0.0d0,0.0d0)
     O_psi = (0.0d0,0.0d0)
     !
     DO ik = 1, nksq
        !
        DO ip = 1, n_ipol
           !
           IF ( n_ipol == 1 ) THEN
              pol_index = ipol 
           ELSE 
              pol_index = ip
           ENDIF
           !
           CALL lr_dvpsi_magnons(ik, pol_index, V0psi(:,:,ik,:,ip))
        ENDDO
        !
        DO ip = 1, n_op
           CALL lr_Opsi_magnons(ik, ip, O_psi(:,:,ik,:,ip))
        ENDDO
        !
     ENDDO
     !
  ELSE
     !
     ! Optical case
     !
     IF (.NOT.d0psi_rs) THEN
        !
        ! Compute dipole in the G space using the commutator [H,r]
        !
        ! TDDFPT+hybrids: Current implementation of the commutator [H,r] 
        ! (commutator_Hx_psi) does not contain the contribution [V_EXX,r], 
        ! hence the relative intensities of the peaks in the spectrum
        ! will be wrong. One should use d0psi_rs=.true. in this case
        ! (see the documentation and the explanation in CPC 185, 2080 (2014)).
        ! 
        IF (lr_exx) WRITE( stdout, '(/5X,"WARNING!!!",              &
             & " TDDFPT + hybrids will give wrong intensities.",    &
             & /5x,"Set d0psi_rs=.true. in order to have correct intensities.")' )
        !
        DO ik = 1, nks
           !
           current_k = ik
           !
           IF ( lsda ) current_spin = isk(ik)
           !
           evc(:,:) = evc0(:,:,ik)
           !
           ! US case: calculate beta-functions vkb.
           !
           CALL init_us_2(ngk(ik), igk_k(:,ik), xk(:,ik), vkb, .true.)
           !$acc update host(vkb)
           !
           ! Compute d0psi = P_c^+ r psi_k 
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
      ELSE
        !
        ! Compute dipole in the R space. This option can be used 
        ! only for finite systems (e.g. molecules).
        !
        CALL compute_d0psi_rs(n_ipol)
        !
        ! In the gamma_only case we compute the kinetic energy
        ! and the beta functions vkb here just once, and then 
        ! they are used throughout the calculation. In the 
        ! k-points version these quantities are always recomputed.
        !
        IF (gamma_only) THEN
           CALL g2_kin(1)
           CALL init_us_2(ngk(1), igk_k(:,1), xk(:,1), vkb, .true.)
           !$acc update host(vkb)
        ENDIF
        !
      ENDIF 
      !
  ENDIF
  !
  IF (gstart==2 .and. gamma_only) d0psi(1,:,:,:) = &
                         & CMPLX(DBLE(d0psi(1,:,:,:)),0.0d0,dp)
  !
  ! Writing of d0psi to the file. 
  ! This is a parallel writing, done in wfc_dir 
  !
  nwordd0psi = 2 * nbnd * npwx * npol * nksq
  !
  tmp_dir_saved = tmp_dir
  !
  IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir
  !
  IF (.not. magnons) THEN
     !
     DO ip = 1, n_ipol
        !
        IF (n_ipol==1) CALL diropn ( iund0psi, 'd0psi.'// &
                    & trim(int_to_char(LR_polarization)), nwordd0psi, exst)
        IF (n_ipol==3) CALL diropn ( iund0psi, 'd0psi.'// &
                    & trim(int_to_char(ip)), nwordd0psi, exst)
        !
        CALL davcio(d0psi(1,1,1,ip),nwordd0psi,iund0psi,1,1)
        !
        CLOSE( UNIT = iund0psi)
        !
     ENDDO
     !
  ELSE
     !
     ! MAGNONS: Writing of V0psi and O_psi to the files
     !
     nwordd0psi = 4 * nbnd_occx * npwx * npol * nksq
     !
     DO ip = 1, n_ipol
        !
        IF (n_ipol==1) CALL diropn ( iund0psi, 'V0psi.'//trim(int_to_char(ipol)), nwordd0psi, exst)
        IF (n_ipol==3) CALL diropn ( iund0psi, 'V0psi.'//trim(int_to_char(ip)), nwordd0psi, exst)
        !
        CALL davcio(V0psi(:,:,:,:,ip),nwordd0psi,iund0psi,1,1)
        !
        CLOSE( UNIT = iund0psi)
        !
     ENDDO
     !
     DO ip = 1, n_op
        !
        CALL diropn ( iund0psi, 'O_psi.'//trim(int_to_char(ip)), nwordd0psi, exst)
        !
        CALL davcio(O_psi(:,:,:,:,ip),nwordd0psi,iund0psi,1,1)
        !
        CLOSE( UNIT = iund0psi)
        !
     ENDDO
     !
  ENDIF
  !
  ! EELS: Writing of d0psi2 to the file.
  !
  IF (eels) THEN
     !
     CALL diropn ( iund0psi, 'd0psi2.'// &
          & trim(int_to_char(LR_polarization)), nwordd0psi, exst)
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
  ! ... Original code from routine compute_dipole
  ! ... modified to calculate the d0psi in the real space
  ! ... by Xiaochuan Ge, Oct, 2013
  ! ... Modified by I. Timrov, June 2016
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at, bg, alat, omega
  USE fft_base,             ONLY : dfftp,dffts
  USE fft_types,            ONLY : fft_index_to_3d
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_barrier
  USE io_global,            ONLY : stdout
  USE wvfct,                ONLY : nbnd,npwx
  USE klist,                ONLY : nks
  USE lr_variables,         ONLY : evc0,sevc0,d0psi,lshift_d0psi
  USE wavefunctions, ONLY : psic
  USE uspp,                 ONLY : okvan
  USE realus,               ONLY : fwfft_orbital_gamma, invfft_orbital_gamma, &
                                   fwfft_orbital_k, invfft_orbital_k
  !
  IMPLICIT NONE
  !
  ! ... Local variables
  !
  REAL(DP),    ALLOCATABLE :: r(:,:)
  COMPLEX(DP), ALLOCATABLE :: psic_temp(:), evc_temp(:,:)
  INTEGER  :: i, j, k, ip, ir, ir_end, ibnb, n_ipol
  REAL(DP) :: inv_nr1, inv_nr2, inv_nr3
  LOGICAL  :: offrange
  !
  WRITE(stdout,'(/,5X,"Calculation of the dipole in real space")')
  !
  IF (.NOT. ALLOCATED(psic)) ALLOCATE(psic(dfftp%nnr))
  ALLOCATE(evc_temp(npwx,nbnd))
  ALLOCATE(psic_temp(dfftp%nnr))
  ALLOCATE(r(dfftp%nnr,n_ipol))
  evc_temp(:,:) = (0.0d0, 0.0d0)
  psic_temp(:)  = (0.0d0, 0.0d0)
  r(:,:) = 0.0d0 
  !
  IF (okvan) THEN
     !
     WRITE(stdout,'(10x,"Real space dipole + USPP is not supported",/)')
#if defined(__MPI)
     CALL mp_barrier(intra_bgrp_comm)
#endif
     STOP
     !
  ENDIF
  !
  ! Calculate r 
  !
  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
  !
#if defined (__MPI)
  ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%my_nr2p*dfftp%my_nr3p)
#else
  ir_end = dfftp%nnr
#endif
  !
  DO ir = 1, ir_end 
     !
     ! ... three dimensional indexes
     !
     CALL fft_index_to_3d (ir, dfftp, i,j,k, offrange)
     IF ( offrange ) CYCLE
     !
     DO ip = 1, n_ipol
        r(ir,ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                   DBLE( j )*inv_nr2*at(ip,2) + &
                   DBLE( k )*inv_nr3*at(ip,3)
     ENDDO
     !
  ENDDO
  !
  IF (lshift_d0psi) CALL shift_d0psi(r,n_ipol)
  !
  ! Calculate the product r * psi(r)
  !
  IF (gamma_only) THEN
     !
     ! gamma_only version
     ! Note: two bands are treated at the same time.
     !
     DO ibnd = 1, nbnd, 2
        !
        ! FFT to R-space: evc_temp -> psic
        !
        CALL invfft_orbital_gamma(evc0(:,:,1), ibnd, nbnd)
        !
        psic_temp(:) = psic(:)
        !
        DO ip = 1, n_ipol
           !
           ! Multiple in R space the dipole and
           ! the ground state wavefunction 
           !
           DO ir = 1, dfftp%nnr
              psic(ir) = r(ir,ip) * alat * psic_temp(ir)
           ENDDO
           !
           ! FFT to G space: psic -> evc_temp
           !
           CALL fwfft_orbital_gamma(evc_temp, ibnd, nbnd) 
           !
           d0psi(:,ibnd,  1,ip) = evc_temp(:,ibnd)
           d0psi(:,ibnd+1,1,ip) = evc_temp(:,ibnd+1)
           !
        ENDDO
        !
     ENDDO
     !
  ELSE
     !
     ! k-points version
     ! Note: one band is treated at a time
     !
     DO ik = 1, nks
        !
        DO ibnd = 1, nbnd
          !
          ! FFT to R-space: evc0 -> psic
          !
          CALL invfft_orbital_k(evc0(:,:,ik), ibnd, nbnd, ik)
          !
          psic_temp(:) = psic(:)
          !
          DO ip = 1, n_ipol
             !
             ! Multiple in R space the dipole and
             ! the ground state wavefunction 
             !
             DO ir = 1, dfftp%nnr
                psic(ir) = r(ir,ip) * alat * psic_temp(ir)
             ENDDO
             !
             ! FFT to G space: psic -> evc_temp
             !
             CALL fwfft_orbital_k(evc_temp, ibnd, nbnd, ik)
             !
             d0psi(:,ibnd,ik,ip) = evc_temp(:,ibnd)
             ! 
          ENDDO
          !
        ENDDO 
        !
     ENDDO 
     !
  ENDIF
  ! 
  ! Orthogonalization: apply P_c to d0psi
  ! Output: d0psi = P_c r * psi(r)
  ! 
  DO ip = 1, n_ipol
     DO ik = 1, nks
        !
        CALL orthogonalize(d0psi(:,:,ik,ip), evc0(:,:,ik), &
                   & ik, ik, sevc0(:,:,ik), ngk(ik), .true.)
        d0psi(:,:,ik,ip) = -d0psi(:,:,ik,ip)
        !
     ENDDO
  ENDDO
  !
  DEALLOCATE(evc_temp)
  DEALLOCATE(psic_temp)
  DEALLOCATE(r)
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
  !
  IMPLICIT NONE
  REAL(dp),intent(inout) :: r(dfftp%nnr,n_ipol)
  integer,intent(in) :: n_ipol
  real(dp) ::  mmin(3), mmax(3), center(3), origin(3), check_cell
  integer :: ip, iatm, ir, ip2
  !
  WRITE(stdout,'(/,5X,"Dipole is shifted to the center of cell", &
                     & " for the calculation of d0psi. ")')
  !
  check_cell = 0.d0
  do ip = 1, 3
    do ip2 = 1,3
      if(.not. ip .eq. ip2) check_cell = check_cell + at(ip,ip2)**2
    enddo
  enddo
  !
  ! XG: I am not sure that this type of super cell is supported now.
  !
  if (check_cell .gt. 1.d-5) call errore('shift_d0psi', &
        & "This type of the supercell is not supported",1)
  !
  mmin(:) = 2000.d0
  mmax(:)= -2000.d0
  center(:) = 0.0d0
  origin(:) = 0.0d0
  !
  do ip = 1, n_ipol
    do iatm = 1, nat
      mmin(ip) = min(mmin(ip), tau(ip,iatm))
      mmax(ip) = max(mmax(ip), tau(ip,iatm))
    enddo
  enddo
  !
  do ip = 1, n_ipol
    center(ip) = 0.5d0*(mmin(ip)+mmax(ip))
    origin(ip) = center(ip)-0.5d0*at(ip,ip)
  enddo
  !
  do ir = 1, dfftp%nnr
    do ip = 1, n_ipol
      r(ir,ip)= r(ir,ip) - origin(ip)
      if (r(ir,ip).lt.0)         r(ir,ip) = r(ir,ip) + at(ip,ip)
      if (r(ir,ip).gt.at(ip,ip)) r(ir,ip) = r(ir,ip) - at(ip,ip)
    enddo
  enddo
  !
  RETURN
  !
END SUBROUTINE shift_d0psi

END SUBROUTINE lr_solve_e
