!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      module nose_electrons

        USE kinds
        USE electrons_nose, ONLY: ekinw => ekincw, fnosee, &
          xenos0 => xnhe0, &
          xenosm => xnhem, &
          xenosp => xnhep, &
          xenos2m => xnhem2

        IMPLICIT NONE
        SAVE

        PRIVATE

        REAL(dbl) :: gkbt
        REAL(dbl) :: wnosee 
        REAL(dbl) :: qnosee
        REAL(dbl) :: enosee

        PUBLIC :: get_nose_electrons, enosee
        PUBLIC :: nose_electrons_setup, set_nose_electrons
        PUBLIC :: noseeinit, update_nose_electrons
        PUBLIC :: nosee_velocity, movenosee

      contains

      subroutine get_nose_electrons(XENOS0_out,XENOSM_out,XENOS2M_out, &
                 XENOSP_out)
        REAL(dbl), intent(out) :: XENOS0_out,XENOSM_out,XENOS2M_out,XENOSP_out
        XENOS0_out = XENOS0
        XENOSM_out = XENOSM
        XENOS2M_out = XENOS2M
        XENOSP_out  = XENOSP
        return
      end subroutine get_nose_electrons

      subroutine set_nose_electrons(XENOS0_inp,XENOSM_inp,XENOS2M_inp, &
                 XENOSP_inp)
        REAL(dbl), intent(in) :: XENOS0_inp,XENOSM_inp,XENOS2M_inp,XENOSP_inp
        XENOS0 = XENOS0_inp
        XENOSM = XENOSM_inp
        XENOS2M = XENOS2M_inp
        XENOSP  = XENOSP_inp
        return
      end subroutine set_nose_electrons


      subroutine nose_electrons_setup(fnosee_inp,ekinw_inp)
        REAL(dbl), intent(in) :: fnosee_inp,ekinw_inp
        fnosee = fnosee_inp
        ekinw  = ekinw_inp
        return
      end subroutine nose_electrons_setup

      subroutine update_nose_electrons
        XENOS2M = XENOSM
        XENOSM  = XENOS0
        XENOS0  = XENOSP
        return
      end subroutine update_nose_electrons

      REAL(dbl) function nosee_velocity()
        use time_step, only: delt
        nosee_velocity = (3.D0*XENOS0-4.D0*XENOSM+XENOS2M)/(2.0d0*delt)
        return
      end function nosee_velocity

      SUBROUTINE MOVENOSEE(EKINC)
        use time_step, only: delt
        REAL(dbl) EKINC
        XENOSP = 2.D0*XENOS0-XENOSM+2.D0*delt*delt/QNOSEE*(EKINC-0.5d0*GKBT)
        ENOSEE = 0.5D0*QNOSEE*((XENOSP-XENOSM)/(2.0d0*delt))**2+GKBT*XENOS0
        return
      end SUBROUTINE MOVENOSEE


      SUBROUTINE NOSEEINIT(ANNEE)

      use constants, only: FACTEM,TERAHERTZ,PI
      use time_step, only: delt
      USE control_flags, ONLY: tnosee
      USE io_global, ONLY: stdout

      IMPLICIT NONE 
      REAL(dbl) ANNEE

      INTEGER NSVAR

      ENOSEE=0.D0
      XENOS0 = 0.D0
      XENOSP = 0.D0
      XENOSM = 0.0d0 
      XENOS2M= 0.D0

      IF( TNOSEE ) THEN
        IF(FNOSEE <= 0.D0) CALL errore('IOSYS','FNOSEE.LE.0',0)
        WNOSEE=FNOSEE*(2.D0*PI)*TERAHERTZ
        QNOSEE=4.D0*EKINW/WNOSEE**2
        GKBT  = 2.D0*EKINW
        XENOS0 = 0.D0
        XENOSM = -ANNEE/DELT
        XENOS2M= 2.D0*XENOSM
        NSVAR=(2.D0*PI)/(WNOSEE*DELT)
        WRITE( stdout,140)
        WRITE( stdout,150) QNOSEE
        WRITE( stdout,160) EKINW
        WRITE( stdout,170) NSVAR
      END IF

 140  FORMAT(//' * TEMPERATURE CONTROL OF',' ELECTRONS WITH NOSE THERMOSTAT')
 150  FORMAT(3X,'NOSE MASS:',F14.7)
 160  FORMAT(3X,'AVERAGED KINETIC ENERGY: ',F12.7)
 170  FORMAT(3X,'== ',I5,' TIME STEPS PER NOSE OSCILL'//)
      
      RETURN 
      END SUBROUTINE noseeinit

      end module nose_electrons
