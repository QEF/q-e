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
  !
  ! ... The quantities needed in calculations with external field
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  LOGICAL :: &
       tefield,      &! if .TRUE. a finite electric field is added to the
                      ! local potential
       dipfield,     &! if .TRUE. the dipole field is subtracted
      ! TB
       relaxz,       &! relax in z direction
       block,        &! add potential barrier
       monopole       ! if .TRUE. and system is charged, charge is represented
                      ! with monopole plane (gate)
  INTEGER :: &
       edir           ! direction of the field
  REAL(DP) :: &
      emaxpos,       &! position of the maximum of the field (0<emaxpos<1)
      eopreg,        &! amplitude of the inverse region (0<eopreg<1)
      eamp,          &! field amplitude (in a.u.) (1 a.u. = 51.44 10^11 V/m)
      etotefield,    &! energy correction due to the field
      el_dipole,     &! electronic_dipole used when dipole correction is on
      ion_dipole,    &! ionic_dipole      used when dipole correction is on
      tot_dipole,    &! total dipole      used when dipole correction is on
      ! TB
      zmon,          &! position of monopole plane
      block_1,       &! blocking potential
      block_2,       &
      block_height,  &
      etotmonofield   ! energy correction due to the monopole
  REAL(DP), ALLOCATABLE :: &
      forcefield(:,:), &
      forcemono(:,:) ! TB monopole forces
  !
CONTAINS
!
!------------------------------------------------------------------------------!
!
      FUNCTION saw(emaxpos,eopreg,x) RESULT (sawout)
        IMPLICIT NONE
        REAL(DP) :: emaxpos,eopreg,x
        REAL(DP) :: y, sawout, z
        
        z = x - emaxpos 
        y = z - floor(z)
        
        if (y.le.eopreg) then
        
            sawout = (0.5_DP - y/eopreg) * (1._DP-eopreg)
        
        else 
!
! I would use:   sawout = y - 0.5_DP * ( 1.0_DP + eopreg )
!
            sawout = (-0.5_DP + (y-eopreg)/(1._DP-eopreg)) * (1._DP-eopreg)
        
        end if
        
      END FUNCTION saw

!TB - start
!------------------------------------------------------------------------------!
!mopopla - add a potential of a monopole plane (kflag = .true.)
!          or the compensating background charge (kflag = .false.)
!          I split those in order to plot both independently
! cite PRB 89, 245406 (2014)
!
      FUNCTION mopopla(zmon,x,kflag) RESULT (mopoplaout)
        IMPLICIT NONE
        REAL(DP) :: zmon,x
        REAL(DP) :: mopoplaout, z
        LOGICAL  :: kflag

        DO ! is x within the cell?
          IF (x>1.0) x=x-1.0
          IF (x<0.0) x=x+1.0
          IF (x<=1.0.and.x>=0.0) EXIT
        ENDDO

        z = (x - zmon)

        !Monopole-plane
        ! if z < 0, we are below the plane
        !    z > 0, above
        !    z < -0.5, the potential is again the same as for z > 0
        !              in order to make it periodic
        !    z > 0.5, the same as for z < 0
        !
        IF (z<=-0.5) z=z+1
        IF (z>=0.5) z=z-1
        IF (z.LE.0) THEN
           IF (kflag) THEN
              mopoplaout = ( 1.0_DP*z )
           ELSE
              mopoplaout = ( z**2 )
           ENDIF
        ELSE
           IF (kflag) THEN
              mopoplaout = ( -1.0_DP*z )
           ELSE
              mopoplaout = ( z**2 )
           ENDIF
        ENDIF

      END FUNCTION mopopla

END MODULE extfield
