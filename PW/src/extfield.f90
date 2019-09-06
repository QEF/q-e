!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE extfield
  !! The quantities needed in calculations with external field.
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  LOGICAL :: tefield
  !! if .TRUE. a finite electric field is added to the
  !! local potential
  LOGICAL :: dipfield
  !! if .TRUE. the dipole field is subtracted
  ! TB
  LOGICAL :: relaxz
  !! relax in z direction
  LOGICAL :: block
  !! add potential barrier
  LOGICAL :: gate
  !! if .TRUE. and system is charged, charge is represented
  !! with charged plate (gate)
  !
  INTEGER :: edir
  !! direction of the field
  REAL(DP) :: emaxpos
  !! position of the maximum of the field (0<emaxpos<1)
  REAL(DP) :: eopreg
  !! amplitude of the inverse region (0<eopreg<1)
  REAL(DP) :: eamp
  !! field amplitude (in a.u.) (1 a.u. = 51.44 10^11 V/m)
  REAL(DP) :: etotefield
  !! energy correction due to the field
  REAL(DP) :: el_dipole
  !! electronic_dipole, used when dipole correction is on
  REAL(DP) :: ion_dipole
  !! ionic_dipole, used when dipole correction is on
  REAL(DP) :: tot_dipole
  !! total dipole, used when dipole correction is on
  ! TB
  REAL(DP) :: zgate
  !! position of charged plate
  REAL(DP) :: block_1
  !! ...blocking potential
  REAL(DP) :: block_2
  !! ...blocking potential
  REAL(DP) :: block_height
  !! ...blocking potential
  REAL(DP) :: etotgatefield
  !! energy correction due to the gate
  REAL(DP), ALLOCATABLE :: forcefield(:,:)
  !
  REAL(DP), ALLOCATABLE :: forcegate(:,:)
  !! TB gate forces
  !
CONTAINS
!
!------------------------------------------------------------------------------!
!
      FUNCTION saw( emaxpos, eopreg, x ) RESULT ( sawout )
        !
        IMPLICIT NONE
        !
        REAL(DP) :: emaxpos, eopreg, x
        REAL(DP) :: y, sawout, z
        !
        z = x - emaxpos 
        y = z - FLOOR(z)
        !
        IF (y <= eopreg) THEN
            !
            sawout = (0.5_DP - y/eopreg) * (1._DP-eopreg)
            !
        ELSE
            !
            ! I would use:   sawout = y - 0.5_DP * ( 1.0_DP + eopreg )
            !
            sawout = (-0.5_DP + (y-eopreg)/(1._DP-eopreg)) * (1._DP-eopreg)
            !
        ENDIF
        !
      END FUNCTION saw

!TB - start
!------------------------------------------------------------------------------!
!
      FUNCTION mopopla( zgate, x, kflag ) RESULT ( mopoplaout )
        !------------------------------------------------------
        !! Add a potential of a charged plate (kflag = .true.)
        !! or the compensating background charge (kflag = .false.).
        !! Cite PRB 89, 245406 (2014).
        !  I split those in order to plot both independently
        !
        IMPLICIT NONE
        !
        REAL(DP) :: zgate,x
        REAL(DP) :: mopoplaout, z
        LOGICAL  :: kflag
        !
        DO ! is x within the cell?
          IF ( x>1.0 ) x=x-1.0
          IF ( x<0.0 ) x=x+1.0
          IF ( x<=1.0.AND.x>=0.0 ) EXIT
        ENDDO
        !
        z = (x - zgate)
        !
        !Charged-plate
        ! if z < 0, we are below the plane
        !    z > 0, above
        !    z < -0.5, the potential is again the same as for z > 0
        !              in order to make it periodic
        !    z > 0.5, the same as for z < 0
        !
        IF ( z<=-0.5 ) z=z+1
        IF ( z>=0.5  ) z=z-1
        IF ( z<=0    ) THEN
           IF (kflag) THEN
              mopoplaout = ( 1.0_DP*z )
           ELSE
              mopoplaout = ( z**2 )
           ENDIF
        ELSE
           IF ( kflag ) THEN
              mopoplaout = ( -1.0_DP*z )
           ELSE
              mopoplaout = ( z**2 )
           ENDIF
        ENDIF
        !
      END FUNCTION mopopla
      !
END MODULE extfield
