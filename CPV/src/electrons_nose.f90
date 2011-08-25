!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
  MODULE electrons_nose
!------------------------------------------------------------------------------!

      USE kinds, ONLY: DP
!
      IMPLICIT NONE
      SAVE

      REAL(DP) :: fnosee   = 0.0_DP   !  frequency of the thermostat ( in THz )
      REAL(DP) :: qne      = 0.0_DP   !  mass of teh termostat
      REAL(DP) :: ekincw   = 0.0_DP   !  kinetic energy to be kept constant

      REAL(DP) :: xnhe0   = 0.0_DP  
      REAL(DP) :: xnhep   = 0.0_DP
      REAL(DP) :: xnhem   = 0.0_DP
      REAL(DP) :: vnhe    = 0.0_DP
      
!
!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!

  subroutine electrons_nose_init( ekincw_ , fnosee_ )
     USE constants, ONLY: pi, au_terahertz
     REAL(DP), INTENT(IN) :: ekincw_, fnosee_
     ! set thermostat parameter for electrons
     qne     = 0.0_DP
     ekincw  = ekincw_
     fnosee  = fnosee_
     xnhe0   = 0.0_DP
     xnhep   = 0.0_DP
     xnhem   = 0.0_DP
     vnhe    = 0.0_DP
     if( fnosee > 0.0_DP ) qne = 4.0_DP * ekincw / ( fnosee * ( 2.0_DP * pi ) * au_terahertz )**2
    return
  end subroutine electrons_nose_init


  function electrons_nose_nrg( xnhe0, vnhe, qne, ekincw )
    !  compute energy term for nose thermostat
    implicit none
    real(dp) :: electrons_nose_nrg
    real(dp), intent(in) :: xnhe0, vnhe, qne, ekincw
      !
      electrons_nose_nrg = 0.5_DP * qne * vnhe * vnhe + 2.0_DP * ekincw * xnhe0
      !
    return
  end function electrons_nose_nrg

  subroutine electrons_nose_shiftvar( xnhep, xnhe0, xnhem )
    !  shift values of nose variables to start a new step
    implicit none
    real(dp), intent(out) :: xnhem
    real(dp), intent(inout) :: xnhe0
    real(dp), intent(in) :: xnhep
      !
      xnhem = xnhe0
      xnhe0 = xnhep
      !
    return
  end subroutine electrons_nose_shiftvar

  subroutine electrons_nosevel( vnhe, xnhe0, xnhem, delt )
    implicit none
    real(dp), intent(inout) :: vnhe
    real(dp), intent(in) :: xnhe0, xnhem, delt 
    vnhe=2.0_DP*(xnhe0-xnhem)/delt-vnhe
    return
  end subroutine electrons_nosevel

  subroutine electrons_noseupd( xnhep, xnhe0, xnhem, delt, qne, ekinc, ekincw, vnhe )
    implicit none
    real(dp), intent(out) :: xnhep, vnhe
    real(dp), intent(in) :: xnhe0, xnhem, delt, qne, ekinc, ekincw
    xnhep = 2.0_DP * xnhe0 - xnhem + 2.0_DP * ( delt**2 / qne ) * ( ekinc - ekincw )
    vnhe  = ( xnhep - xnhem ) / ( 2.0_DP * delt )
    return
  end subroutine electrons_noseupd


  SUBROUTINE electrons_nose_info( delt)

      use constants,     only: au_terahertz, pi
      USE io_global,     ONLY: stdout
      USE control_flags, ONLY: tnosee

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: delt

      INTEGER   :: nsvar
      REAL(DP) :: wnosee

      IF( tnosee ) THEN
        !
        IF( fnosee <= 0.0_DP) &
          CALL errore(' electrons_nose_info ', ' fnosee less than zero ', 1)
        IF( delt <= 0.0_DP) &
          CALL errore(' electrons_nose_info ', ' delt less than zero ', 1)

        wnosee = fnosee * ( 2.0_DP * pi ) * au_terahertz
        nsvar  = ( 2.0_DP * pi ) / ( wnosee * delt )

        WRITE( stdout,563) ekincw, nsvar, fnosee, qne
      END IF

 563  format( //, &
     &       3X,'electrons dynamics with nose` temperature control:', /, &
     &       3X,'Kinetic energy required   = ', f10.5, ' (a.u.) ', /, &
     &       3X,'time steps per nose osc.  = ', i5, /, &
     &       3X,'nose` frequency           = ', f10.3, ' (THz) ', /, &
     &       3X,'nose` mass(es)            = ', 20(1X,f10.3),//)


    RETURN
  END SUBROUTINE electrons_nose_info

!------------------------------------------------------------------------------!
  END MODULE electrons_nose
!------------------------------------------------------------------------------!
