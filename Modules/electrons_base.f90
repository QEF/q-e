!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
  MODULE electrons_base
!------------------------------------------------------------------------------!

      USE kinds, ONLY: dbl
!
      IMPLICIT NONE
      SAVE

      INTEGER :: nbnd       = 0    !  total number of electronic states
      INTEGER :: nbndx      = 0    !  array dimension nbndx >= nbnd
      INTEGER :: nel(2)     = 0    !  number of electrons (up, down)
      INTEGER :: nelt       = 0    !  total number of electrons ( up + down )
      INTEGER :: nupdwn(2)  = 0    !  number of states with spin up (1) and down (2)
      INTEGER :: iupdwn(2)  = 0    !  first state with spin (1) and down (2)
      INTEGER :: nspin      = 0    !  nspin = number of spins (1=no spin, 2=LSDA)
      INTEGER :: nudx       = 0    !  max (nupdw(1),nupdw(2))

      REAL(dbl), ALLOCATABLE :: f(:)   ! occupation numbers ( at gamma )
      REAL(dbl) :: qbac = 0.0d0        ! background neutralizing charge
      INTEGER, ALLOCATABLE :: fspin(:) ! spin of each state
!
!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!

    SUBROUTINE electrons_base_init()
      IF (nspin.EQ.1) THEN 
        nelt=nel(1)
        nudx=nupdwn(1)
      ELSE
        nelt=nel(1)+nel(2)
        nudx=MAX(nupdwn(1),nupdwn(2))
      END IF
      RETURN
    END SUBROUTINE

    SUBROUTINE deallocate_elct()
      IF( ALLOCATED( f ) ) DEALLOCATE( f )
      IF( ALLOCATED( fspin ) ) DEALLOCATE( fspin )
      RETURN
    END SUBROUTINE


!------------------------------------------------------------------------------!
  END MODULE electrons_base
!------------------------------------------------------------------------------!



!------------------------------------------------------------------------------!
  MODULE electrons_nose
!------------------------------------------------------------------------------!

      USE kinds, ONLY: dbl
!
      IMPLICIT NONE
      SAVE

      REAL(dbl) :: fnosee   = 0.0d0   !  frequency of the thermostat ( in THz )
      REAL(dbl) :: qne      = 0.0d0   !  mass of teh termostat
      REAL(dbl) :: ekincw   = 0.0d0   !  kinetic energy to be kept constant

      REAL(dbl) :: xnhe0   = 0.0d0   
      REAL(dbl) :: xnhep   = 0.0d0   
      REAL(dbl) :: xnhem   = 0.0d0   
      REAL(dbl) :: xnhem2  = 0.0d0   
      REAL(dbl) :: vnhe    = 0.0d0   
!
!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!

  subroutine electrons_nose_init( ekincw_ , fnosee_ )
     USE constants, ONLY: factem, pi, terahertz
     REAL(dbl), INTENT(IN) :: ekincw_, fnosee_
     ! set thermostat parameter for electrons
     qne    = 0.0d0
     ekincw  = ekincw_
     fnosee = fnosee_
     if( fnosee > 0.0d0 ) qne = 4.d0*ekincw/(fnosee*(2.d0*pi)*terahertz)**2
    return
  end subroutine


  function electrons_nose_nrg( xnhe0, vnhe, qne, ekincw )
    !  compute energy term for nose thermostat
    implicit none
    real(kind=8) :: electrons_nose_nrg
    real(kind=8), intent(in) :: xnhe0, vnhe, qne, ekincw
      !
      electrons_nose_nrg = 0.5d0 * qne * vnhe * vnhe + 2.0d0 * ekincw * xnhe0
      !
    return
  end function

  subroutine electrons_nose_shiftvar( xnhep, xnhe0, xnhem )
    !  shift values of nose variables to start a new step
    implicit none
    real(kind=8), intent(out) :: xnhem, xnhe0
    real(kind=8), intent(in) :: xnhep
      !
      xnhem = xnhe0
      xnhe0 = xnhep
      !
    return
  end subroutine



!------------------------------------------------------------------------------!
  END MODULE electrons_nose
!------------------------------------------------------------------------------!
