!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE kcw_prepare_q(do_band, setup_pw, iq)
  !-----------------------------------------------------------------------
  !
  !!  This routine is ADAPTED from PH/prepare_q.f90
  !!  This routine prepares a few variables that are needed to control
  !!  the LR run after the q point has been decided, but before
  !!  doing the band calculation. In particular it sets:
  !!  xq : the q point for the LR calculation
  !!  lgamma : if this is a gamma point calculation
  !!  setup_pw : if .true. the pw_setup has to be run
  !!  do_band : if .true. the bands need to be calculated before phonon
  !!  nbnd: REDEFINE the minimum number of bands for the nscf 
  !!        calculation (only occupied manifold)
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE qpoint,               ONLY : xq
  USE control_kcw,          ONLY : tmp_dir_kcw, tmp_dir_kcwq, x_q 
  USE control_lr,           ONLY : lgamma
  USE io_files,             ONLY : create_directory
  USE mp,                   ONLY : mp_bcast
  USE cell_base,            ONLY : at
  USE wvfct,                ONLY : nbnd
  USE klist,                ONLY : nelup, neldw, nelec, lgauss, ltetra
  USE start_k,              ONLY : reset_grid
  USE noncollin_module,     ONLY : domag, noncolin
  USE lsda_mod,             ONLY : lsda
  USE control_kcw,          ONLY : irr_bz
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  LOGICAL, INTENT(OUT) :: do_band, setup_pw
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  INTEGER :: nbnd_old, degspin
  REAL(DP) :: xq_(3)
  !
  tmp_dir_kcwq=tmp_dir_kcw
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
  IF (.NOT.lgamma .OR. (domag .AND. noncolin) .OR. irr_bz) &
     tmp_dir_kcwq= TRIM (tmp_dir_kcw) // 'q' &
                   & // TRIM(int_to_char(iq))//'/'
  !
  ! In the case of q != 0, we make first a non selfconsistent run
  ! But first we reset the grid of irrebucible Kpoints (NOT WORKING YET!!)
  ! The following line should be there to enable only non-equivalent kpoint (if nrot .gt. 1), BUT STILL NOT WORKING
  !
  !newgrid = reset_grid ( nk1, nk2, nk3, k1, k2, k3 )
  !IF (nrot .gt. 1) newgrid = reset_grid ( 2, 2, 2, 0, 0, 0  )
  !modenum = 0
  setup_pw = (.NOT.lgamma)
  IF(irr_bz) setup_pw = .TRUE.
  !
  IF (noncolin.AND.domag) setup_pw=.true. !! NsC need to check this is needed (see comment in PH/prepare_q.f90)
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
     !IF (noncolin .OR. .NOT. lsda) degspin = 1
     IF (noncolin) degspin = 1
     !
     nbnd = MAX ( NINT( nelec / degspin ), NINT(nelup), NINT(neldw) ) + 3
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
  xq_ = x_q(:,iq)
  CALL cryst_to_cart(1, xq_, at, -1)
  WRITE(stdout,'(/,/,5X, 78("="))')
  WRITE( stdout, '(5X,"Calculation of q = ",3F12.7, "  [Cart ]")') x_q(:,iq)
  WRITE( stdout, '(5X,"Calculation of q = ",3F12.7, "  [Cryst]")') xq_(:)
  WRITE(stdout,'(5X, 78("="),/)')
  !
  RETURN
  !
END SUBROUTINE kcw_prepare_q
