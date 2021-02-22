!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE kc_prepare_q(do_band, do_iq, setup_pw, iq)
  !-----------------------------------------------------------------------
  !
  !!  This routine is ADAPTED from PH/prepare_q.f90
  !!  This routine prepares a few variables that are needed to control
  !!  the LR run after the q point has been decided, but before
  !!  doing the band calculation. In particular it sets:
  !!  xq : the q point for the phonon calculation
  !!  lgamma : if this is a gamma point calculation
  !!  current_iq : the current q point
  !!  do_iq : if .true. q point has to be calculated
  !!  setup_pw : if .true. the pw_setup has to be run
  !!  do_band : if .true. the bands need to be calculated before phonon
  !!  nbnd: REDEFINE the minimum number of bands for the nscf 
  !!        calculation (only occupied manifold)
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE qpoint,               ONLY : xq
  USE control_ph,           ONLY : current_iq, tmp_dir_ph, tmp_dir_phq, u_from_file!, newgrid
  USE control_lr,           ONLY : lgamma
  USE io_files,             ONLY : prefix
  USE output,               ONLY : fildyn
  USE ph_restart,           ONLY : ph_writefile
  USE io_files,             ONLY : prefix, create_directory
  USE mp,                   ONLY : mp_bcast
  USE disp,                 ONLY : x_q
  USE cell_base,            ONLY : at
  USE wvfct,                ONLY : nbnd
  USE klist,                ONLY : nelup, neldw, nelec, lgauss, ltetra
  USE start_k,              ONLY : reset_grid
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  LOGICAL, INTENT(OUT) :: do_band, do_iq, setup_pw
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  INTEGER :: ierr, nbnd_old, degspin
  REAL(DP) :: xq_(3)
  !
  do_iq=.TRUE.
  !
  current_iq = iq
  !
  tmp_dir_phq=tmp_dir_ph
  !
  !
  ! ... set the name for the output file
  !
  fildyn = TRIM( prefix ) // TRIM( int_to_char( iq ) )
  !
  ! ... set the q point
  !
  xq(1:3)  = x_q(1:3,iq)
  !
  !  Check if it is lgamma
  !
  lgamma = ( xq(1) == 0.D0 .AND. xq(2) == 0.D0 .AND. xq(3) == 0.D0 )
  !
  ! ... each q /= gamma is saved on a different directory
  !
  IF (.NOT.lgamma) &
     tmp_dir_phq= TRIM (tmp_dir_ph) // TRIM(prefix) // '_q' &
                   & // TRIM(int_to_char(iq))//'/'
  !
  !  Save the current status of the run: all the flags, the list of q,
  !  and the current q, the fact that we are before the bands
  !
  !CALL ph_writefile('init',0,0,ierr)
  !
  ! In the case of q != 0, we make first a non selfconsistent run
  ! But first we reset the grid of irrebucible Kpoints (NOT WORKING YET!!)
  ! The following line should be there to enable only non-equivalent kpoint (if nrot .gt. 1), BUT STILL NOT WORKING
  !
  !newgrid = reset_grid ( nk1, nk2, nk3, k1, k2, k3 )
  !IF (nrot .gt. 1) newgrid = reset_grid ( 2, 2, 2, 0, 0, 0  )
  !modenum = 0
  setup_pw = (.NOT.lgamma)
  !
  ! Only the occupied bands neeed to be computed. 
  ! In general we might have run PWSCF with a huge number of Empty states 
  ! needed for the Wannierization. Here we reset the number of bands to be 
  ! computed at k+q. 
  ! 
  IF( setup_pw) THEN 
     !
     nbnd_old = nbnd 
     degspin = 2
     !
     nbnd = MAX ( NINT( nelec / degspin ), NINT(nelup), NINT(neldw) ) + 2
     !
     IF ( lgauss .OR. ltetra ) THEN
        !
        ! ... metallic case: add 20% more bands, with a minimum of 4
        !
        nbnd = MAX( NINT( 1.2D0 * nelec / degspin ), &
                    NINT( 1.2D0 * nelup), NINT( 1.2d0 * neldw ), &
                    ( nbnd + 4 ) )
        !
    ENDIF
    !
    IF (nbnd /= nbnd_old) WRITE(stdout,'(/,8X, "INFO: nbnd REDIFINED", i5, " --> ", i5)') nbnd_old, nbnd
    !
  ENDIF
  !
  do_band=.TRUE.
  !
  u_from_file=.FALSE. ! NsC 
  !
  xq_ = x_q(:,iq)
  CALL cryst_to_cart(1, xq_, at, -1)
  WRITE(stdout,'(/,/,5X, 78("="))')
  WRITE( stdout, '(5X,"Calculation of q = ",3F12.7, "  [Cart ]")') x_q(:,iq)
  WRITE( stdout, '(5X,"Calculation of q = ",3F12.7, "  [Cryst]")') xq_(:)
  WRITE(stdout,'(5X, 78("="),/)')
  !
  RETURN
  !
END SUBROUTINE kc_prepare_q
