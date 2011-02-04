!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
  MODULE cell_nose
!------------------------------------------------------------------------------!

      USE kinds, ONLY : DP
!
      IMPLICIT NONE
      SAVE

      REAL(DP) :: xnhh0(3,3) = 0.0_DP
      REAL(DP) :: xnhhm(3,3) = 0.0_DP
      REAL(DP) :: xnhhp(3,3) = 0.0_DP
      REAL(DP) :: vnhh(3,3)  = 0.0_DP
      REAL(DP) :: temph      = 0.0_DP  !  Thermostat temperature (from input)
      REAL(DP) :: fnoseh     = 0.0_DP  !  Thermostat frequency (from input)
      REAL(DP) :: qnh        = 0.0_DP  !  Thermostat mass (computed)

CONTAINS

  subroutine cell_nose_init( temph_init, fnoseh_init )
     USE constants, ONLY: pi, au_terahertz, k_boltzmann_au
     REAL(DP), INTENT(IN) :: temph_init, fnoseh_init
     ! set thermostat parameter for cell
     qnh    = 0.0_DP
     temph  = temph_init
     fnoseh = fnoseh_init
     if( fnoseh > 0.0_DP ) qnh = 2.0_DP * ( 3 * 3 ) * temph * k_boltzmann_au / &
          (fnoseh*(2.0_DP*pi)*au_terahertz)**2
    return
  end subroutine cell_nose_init

  subroutine cell_nosezero( vnhh, xnhh0, xnhhm )
    real(DP), intent(out) :: vnhh(3,3), xnhh0(3,3), xnhhm(3,3)
    xnhh0=0.0_DP
    xnhhm=0.0_DP
    vnhh =0.0_DP
    return
  end subroutine cell_nosezero

  subroutine cell_nosevel( vnhh, xnhh0, xnhhm, delt )
    implicit none
    REAL(DP), intent(inout) :: vnhh(3,3)
    REAL(DP), intent(in) :: xnhh0(3,3), xnhhm(3,3), delt
    vnhh(:,:)=2.0_DP*(xnhh0(:,:)-xnhhm(:,:))/delt-vnhh(:,:)
    return
  end subroutine cell_nosevel

  subroutine cell_noseupd( xnhhp, xnhh0, xnhhm, delt, qnh, temphh, temph, vnhh )
    use constants, only: k_boltzmann_au
    implicit none
    REAL(DP), intent(out) :: xnhhp(3,3), vnhh(3,3)
    REAL(DP), intent(in) :: xnhh0(3,3), xnhhm(3,3), delt, qnh, temphh(3,3), temph
    integer :: i, j
    do j=1,3
      do i=1,3
        xnhhp(i,j) = 2.0_DP*xnhh0(i,j)-xnhhm(i,j) + &
             (delt**2/qnh)* k_boltzmann_au * (temphh(i,j)-temph)
        vnhh(i,j) =(xnhhp(i,j)-xnhhm(i,j))/( 2.0_DP * delt )
      end do
    end do
    return
  end subroutine cell_noseupd

  
  REAL(DP) function cell_nose_nrg( qnh, xnhh0, vnhh, temph, iforceh )
    use constants, only: k_boltzmann_au
    implicit none
    REAL(DP) :: qnh, vnhh( 3, 3 ), temph, xnhh0( 3, 3 )
    integer :: iforceh( 3, 3 )
    integer :: i, j
    REAL(DP) :: enij
    cell_nose_nrg = 0.0_DP
    do i=1,3
      do j=1,3
        enij = 0.5_DP*qnh*vnhh(i,j)*vnhh(i,j)+temph*k_boltzmann_au*xnhh0(i,j)
        cell_nose_nrg = cell_nose_nrg + iforceh( i, j ) * enij
      enddo
    enddo
    return
  end function cell_nose_nrg

  subroutine cell_nose_shiftvar( xnhhp, xnhh0, xnhhm )
    !  shift values of nose variables to start a new step
    implicit none
    REAL(DP), intent(out) :: xnhhm(3,3)
    REAL(DP), intent(inout) :: xnhh0(3,3)
    REAL(DP), intent(in) :: xnhhp(3,3)
      xnhhm = xnhh0
      xnhh0 = xnhhp
    return
  end subroutine cell_nose_shiftvar


  SUBROUTINE cell_nose_info ( delt )

      use constants,     only: au_terahertz, pi
      USE io_global,     ONLY: stdout
      USE control_flags, ONLY: tnoseh

      IMPLICIT NONE

      REAL(DP), INTENT (IN) :: delt

      INTEGER   :: nsvar
      REAL(DP) :: wnoseh

      IF( tnoseh ) THEN
        !
        IF( fnoseh <= 0.0_DP) &
          CALL errore(' cell_nose_info ', ' fnoseh less than zero ', 1)
        IF( delt <= 0.0_DP) &
          CALL errore(' cell_nose_info ', ' delt less than zero ', 1)

        wnoseh = fnoseh * ( 2.0_DP * pi ) * au_terahertz
        nsvar  = ( 2.0_DP * pi ) / ( wnoseh * delt )

        WRITE( stdout,563) temph, nsvar, fnoseh, qnh
      END IF

 563  format( //, &
     &       3X,'cell dynamics with nose` temperature control:', /, &
     &       3X,'Kinetic energy required   = ', f10.5, ' (Kelvin) ', /, &
     &       3X,'time steps per nose osc.  = ', i5, /, &
     &       3X,'nose` frequency           = ', f10.3, ' (THz) ', /, &
     &       3X,'nose` mass(es)            = ', 20(1X,f10.3),//)

    RETURN
  END SUBROUTINE cell_nose_info


!
!------------------------------------------------------------------------------!
   END MODULE cell_nose
!------------------------------------------------------------------------------!

