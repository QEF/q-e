!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE nose_ions
!=----------------------------------------------------------------------------=!

        USE kinds
        USE control_flags, ONLY: tnosep
        USE ions_nose, ONLY: &
          xnos0 => xnhp0, &
          xnosm => xnhpm, &
          xnos2m => xnhpm2, &
          xnosp => xnhpp, &
          tempw, & 
          fnosep, &
          gkbt

        IMPLICIT NONE
        SAVE

        PRIVATE

        REAL(dbl) :: ekincm
        REAL(dbl) :: enosep
        REAL(dbl) :: wnosep
        REAL(dbl) :: qnosep

        PUBLIC :: enosep, get_nose_ions, set_nose_ions
        PUBLIC :: nose_ions_setup, nosep_velocity
        PUBLIC :: nosepinit, update_nose_ions, movenosep

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!

      subroutine get_nose_ions(XNOS0_out,XNOSM_out,XNOS2M_out,XNOSP_out)
        REAL(dbl), intent(out) :: XNOS0_out,XNOSM_out,XNOS2M_out,XNOSP_out
        XNOS0_out = XNOS0
        XNOSM_out = XNOSM
        XNOS2M_out = XNOS2M
        XNOSP_out  = XNOSP
        return
      end subroutine get_nose_ions

      subroutine set_nose_ions(XNOS0_inp,XNOSM_inp,XNOS2M_inp,XNOSP_inp)
        REAL(dbl), intent(in) :: XNOS0_inp,XNOSM_inp,XNOS2M_inp,XNOSP_inp
        XNOS0 = XNOS0_inp
        XNOSM = XNOSM_inp
        XNOS2M = XNOS2M_inp
        XNOSP  = XNOSP_inp
        return
      end subroutine set_nose_ions

      subroutine nose_ions_setup( fnosep_inp, tempw_inp )
        REAL(dbl),  intent(in) :: fnosep_inp
        REAL(dbl),  intent(in) :: tempw_inp
        fnosep = fnosep_inp
        tempw  = tempw_inp
        return
      end subroutine nose_ions_setup

      subroutine update_nose_ions
        XNOS2M = XNOSM
        XNOSM  = XNOS0
        XNOS0  = XNOSP
        return
      end subroutine update_nose_ions

      REAL(dbl) FUNCTION nosep_velocity()
        USE time_step, ONLY: delt
        nosep_velocity = ( 3.D0 * xnos0 - 4.D0 * xnosm + xnos2m ) / ( 2.0d0 * delt )
        RETURN
      END FUNCTION nosep_velocity

      SUBROUTINE movenosep( ekinp )
        USE time_step, only: delt
        IMPLICIT NONE
        REAL(dbl), INTENT(IN) :: ekinp
        xnosp  = 2.D0 * xnos0 - xnosm + 2.D0 * delt * delt / qnosep * ( ekinp -  0.5D0 * gkbt )
        enosep = 0.5D0 * qnosep * ( ( xnosp - xnosm ) / ( 2.0d0 * delt ) )**2 + gkbt * xnos0
        RETURN
      END SUBROUTINE movenosep

!=----------------------------------------------------------------------------=!

      SUBROUTINE NOSEPINIT(GLIB)

      use constants, only: FACTEM,TERAHERTZ,PI
      use time_step, only: delt
      USE io_global, ONLY: stdout

      IMPLICIT NONE      

      REAL(dbl) GLIB
      INTEGER IS,NSVAR

      ENOSEP = 0.D0
      XNOS2M = 0.D0
      XNOSM  = 0.D0
      XNOS0  = 0.D0
      XNOSP  = 0.D0

      IF( TNOSEP ) THEN
        IF( FNOSEP <= 0.D0) CALL errore('IOSYS','FNOSEP.LE.0',0)
        GKBT   = GLIB*TEMPW/FACTEM
        WNOSEP = FNOSEP * ( 2.D0 * PI ) * TERAHERTZ
        QNOSEP = 2.D0*GKBT/WNOSEP**2
        XNOS2M = 0.D0
        XNOSM  = 0.D0
        XNOS0  = 0.D0
        NSVAR  = (2.D0*PI)/(WNOSEP*DELT)

        WRITE( stdout,100) 
        WRITE( stdout,110) QNOSEP,TEMPW
        WRITE( stdout,120) GLIB
        WRITE( stdout,130) NSVAR
      END IF

 100  FORMAT(//' * TEMPERATURE CONTROL OF IONS WITH NOSE THERMOSTAT'/)
 110  FORMAT(3X,'NOSE MASS:',F12.4,' TEMPERATURE (K):',F12.4)
 120  FORMAT(3X,'IONIC DEGREES OF FREEDOM: ',F5.0)
 130  FORMAT(3X,'== ',I5,' TIME STEPS PER NOSE OSCILL'//)

      RETURN
      END SUBROUTINE NOSEPINIT

!=----------------------------------------------------------------------------=!
   END MODULE nose_ions
!=----------------------------------------------------------------------------=!
