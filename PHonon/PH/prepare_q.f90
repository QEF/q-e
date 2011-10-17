!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE prepare_q(auxdyn, do_band, do_iq, setup_pw, iq)
  !-----------------------------------------------------------------------
  !
  !  This routine prepares a few variables that are needed to control
  !  the phonon run after the q point has been decided, but before
  !  doing the band calculation. In particular if ldisp=true it sets:
  !  xq : the q point for the phonon calculation
  !  fildyn : the name of the dynamical matrix
  !  lgamma : if this is a gamma point calculation
  !  epsil and zue : if epsil and zue need to be calculated
  !  In all cases it sets:
  !  current_iq : the current q point
  !  do_iq : if .true. q point has to be calculated
  !  setup_pw : if .true. the pw_setup has to be run
  !  do_band : if .true. the bands need to be calculated before phonon
  !
  USE control_flags,   ONLY : modenum
  USE io_global,       ONLY : stdout, ionode
  USE klist,           ONLY : lgauss
  USE qpoint,          ONLY : xq
  USE disp,            ONLY : x_q, done_iq, rep_iq, done_rep_iq, comp_iq
  USE control_ph,      ONLY : ldisp, lgamma, epsil, trans, zue, zeu, &
                              start_irr, last_irr, current_iq, newgrid, &
                              done_bands, tmp_dir_ph, tmp_dir_phq, lqdir
  USE io_files,        ONLY : prefix
  USE ramanm,          ONLY : lraman, elop
  USE freq_ph,         ONLY : fpol
  USE output,          ONLY : fildyn, fildvscf
  USE ph_restart,      ONLY : ph_writefile
  USE el_phon,         ONLY : elph_mat, wan_index_dyn, auxdvscf
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  LOGICAL, INTENT(OUT) :: do_band, do_iq, setup_pw
  CHARACTER (LEN=256), INTENT(IN) :: auxdyn
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  INTEGER :: irr
  !
  do_iq=.TRUE.
  !
  ! Case 1) This q point is not calculated because not requested in this run
  !
  IF ( comp_iq(iq)==0 ) THEN
     do_iq=.FALSE.
     RETURN
  ENDIF
  !
  !  Case 2) This q point is not calculated because it has too few
  !          representation and the starting representation is larger
  !          than the number of available representations
  !
  IF (start_irr>rep_iq(iq)) THEN
     WRITE(6,'(5x,"Exiting... start_irr,",i4,&
            & " > number of representations,",i4   )') &
            start_irr, rep_iq(iq)
     do_iq=.FALSE.
     RETURN
  ENDIF
  !
  current_iq = iq
  !
  tmp_dir_phq=tmp_dir_ph
  !
  IF ( ldisp ) THEN
     !
     ! ... set the name for the output file
     !
     if(elph_mat) then
        fildyn = TRIM( auxdyn ) // TRIM( int_to_char( wan_index_dyn(iq) ) )
        fildvscf = TRIM( auxdvscf ) // TRIM( int_to_char( iq ) ) // '_'
     else
        fildyn = TRIM( auxdyn ) // TRIM( int_to_char( iq ) )
     endif
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
     IF (.NOT.lgamma.AND.lqdir) &
         tmp_dir_phq= TRIM (tmp_dir_ph) // TRIM(prefix) // '_q' &
                       & // TRIM(int_to_char(iq))//'/'
     !
     IF ( lgamma ) THEN
        !
        IF ( .NOT. lgauss ) THEN
           !
           ! ... in the case of an insulator at q=0 one has to calculate
           ! ... the dielectric constant and the Born eff. charges
           ! ... the other flags depend on input
           !
           epsil = .TRUE.
           zeu = .TRUE.
           zue = .TRUE.
           !
        ELSE
           !
           ! For a metal no electric field perturbation is available
           !
           epsil = .FALSE.
           zeu   = .FALSE.
           zue   = .FALSE.
           elop  = .FALSE.
           lraman = .FALSE.
           fpol =.FALSE.
           !
        END IF
        !
     ELSE
        !
        ! ... for q /= 0 no calculation of the dielectric tensor,
        ! ...            Born eff. charges, electro-optic, raman or
        ! ...            frequency dependent tensor
        !
        epsil = .FALSE.
        zue   = .FALSE.
        zeu   = .FALSE.
        elop  = .FALSE.
        lraman = .FALSE.
        fpol =.FALSE.
        !
        !
     END IF
  ENDIF
  !
  !  Save the current status of the run: all the flags, the list of q,
  !  and the current q, the fact that we are before the bands
  !
  CALL ph_writefile('init',0)
  !
  ! ... In the case:
  !     of q = 0 and one of nk1, nk2 or nk3 = 0 (newgrid=.false.)
  !              we do not make a non selfconsistent run
  !     of q = 0 and nk1*nk2*nk3 \=0 (newgrid = .true.) 
  !              we do make first a nscf run
  !     of q \= 0   we do make first a nscf run
  !
  setup_pw = (.NOT.lgamma.OR.modenum /= 0 .OR. newgrid).AND..NOT. done_bands
  do_band=.FALSE.
  DO irr=start_irr, MIN(ABS(last_irr),rep_iq(iq))
     IF (done_rep_iq(irr,iq) /= 1) THEN
        do_band=.TRUE.
        EXIT
     ENDIF
  ENDDO
!
!  If this q has been already calculated we only diagonalize the dynamical 
!  matrix
!

  IF ( done_iq(iq) == 1 ) do_band=.FALSE.

  WRITE( stdout, '(/,5X,"Calculation of q = ",3F12.7)') xq

  RETURN
  !
END SUBROUTINE prepare_q
