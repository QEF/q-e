!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE dynamics()
  !----------------------------------------------------------------------------
  !
  ! ... This routine performs one step of molecular dynamics evolution
  ! ... using the Verlet algorithm. 
  !
  ! ... Parameters:
  ! ... mass         mass of the atoms
  ! ... dt           time step
  ! ... temperature  starting temperature
  ! ...              The starting velocities of atoms are set accordingly
  ! ...              to the starting temperature, in random directions.
  ! ...              The initial velocity distribution is therefore a constant
  !
  ! ... delta_t, nraise are used to change the temperature as follows:
  !
  ! ... delta_t = 1 :                  every 'nraise' step the actual 
  ! ...                                temperature is rescaled to the initial 
  ! ...                                value
  ! ... delta_t /= 1 and delta_t > 0 : at each step the actual temperature is
  ! ...                                multiplied by delta_t; this is done
  ! ...                                rescaling all the velocities
  ! ... delta_t < 0 :                  every 'nraise' step the temperature
  ! ...                                reduced by -delta_t
  !
  ! ... Dario Alfe' 1997  and  Carlo Sbraccia 2004-2005
  !
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  USE constants,     ONLY : amconv, eps8, convert_E_to_temp
  USE ions_base,     ONLY : nat, nsp, ityp, tau, if_pos, atm
  USE cell_base,     ONLY : alat
  USE dynam,         ONLY : amass, temperature, dt, delta_t, nraise
  USE ener,          ONLY : etot
  USE force_mod,     ONLY : force
  USE relax,         ONLY : epse, epsf
  USE control_flags, ONLY : istep, nstep, conv_ions, &
                            lconstrain, ldamped, lfixatom, lrescale_t
  USE io_files,      ONLY : prefix
  !
  USE constraints_module
  USE basic_algebra_routines
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: tau_old(:,:), tau_new(:,:), vel(:,:), acc(:,:)
  REAL(DP), ALLOCATABLE :: mass(:)
  REAL(DP)              :: ekin, etotold
  REAL(DP)              :: total_mass, temp_new, elapsed_time, norm_of_dtau
  REAL(DP)              :: ml(3), mlt
  INTEGER               :: i, na
  LOGICAL               :: exst
  !
  !
  ALLOCATE( mass( nat ) )
  !
  ALLOCATE( tau_old( 3, nat ) )
  ALLOCATE( tau_new( 3, nat ) )
  ALLOCATE( vel(     3, nat ) )
  ALLOCATE( acc(     3, nat ) )
  !
  tau_old = tau
  tau_new = 0.D0
  vel     = 0.D0
  acc     = 0.D0
  !
  IF ( istep == 1 ) THEN
     !
     IF ( ldamped ) THEN
        !
        WRITE( UNIT = stdout, &
               FMT = '(/,5X,"Damped Dynamics Calculation")' )
        !
     ELSE
        !
        WRITE( UNIT = stdout, &
               FMT = '(/,5X,"Molecular Dynamics Calculation")' )
        !
     END IF
     !
  END IF
  !
  ! ... one Ryd a.u. of time is 4.84*10^-17 seconds, i.e. 0.0484  femtoseconds
  !
  CALL seqopn( 4, 'md', 'FORMATTED', exst )
  !
  IF ( .NOT. exst ) THEN
     !
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
     WRITE( UNIT = stdout, &
            FMT = '(/,5X,"Starting temperature",T27," = ",F8.2," K")' ) &
         temperature
     !
     DO na = 1, nsp
        !
        WRITE( UNIT = stdout, &
               FMT = '(5X,"mass ",A2,T27," = ",F8.2)' ) atm(na), amass(na)
        !
     END DO
     !
     WRITE( UNIT = stdout, &
            FMT = '(5X,"Time step",T27," = ",F8.2," a.u.,",F8.4, &
                     & " femto-seconds")' ) dt, ( dt * 0.0484D0 )
     !
     ! ...  masses in rydberg atomic units
     !
     total_mass = 0.D0
     !
     DO na = 1, nat
        !
        mass(na) = amass( ityp(na) ) * amconv
        !
        total_mass = total_mass + mass(na)
        !
     END DO
     !
     IF ( lrescale_t ) THEN
        !
        ! ... initial thermalization. N.B. tau is in units of alat
        !
        CALL start_therm()
        !
        ! ... vel is used already multiplied by dt (it is used as an auxiliary 
        ! ... variable in the Verlet scheme)
        !
        vel = dt * vel
        !
     ELSE
        !
        tau_old = tau
        vel     = 0.D0
        !
     END IF
     !
     elapsed_time = 0.D0
     !
     temp_new = temperature
     !
     istep = 0
     !
  ELSE
     !
     ! ... the file is read
     !
     READ( UNIT = 4, FMT = * ) &
        etotold, temp_new, mass, total_mass, elapsed_time, istep, tau_old
     !
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     !
     ! ... the velocity is computed here (used as an auxiliary variable)
     !
     vel = ( tau - tau_old ) * DBLE( if_pos )
     ! 
  END IF
  !
  ! ... elapsed_time is in picoseconds
  !
  elapsed_time = elapsed_time + dt * 0.0000484D0
  !
  istep = istep + 1
  !
  IF ( lrescale_t ) THEN
     !
     IF ( MOD( istep, nraise ) == 0 ) THEN
        !
        IF ( delta_t == 1.D0 ) THEN
           !
           CALL thermalize( temp_new, temperature )
           !
        ELSE IF ( delta_t < 0 ) THEN
           !
           WRITE( UNIT = stdout, &
                  FMT = '(/,5X,"Thermalization: delta_t = ",F6.3, &
                          & ", T = ",F6.1)' )  - delta_t, ( temp_new - delta_t )
           !
           CALL thermalize( temp_new, ( temp_new - delta_t ) )
           !
        END IF
        !
     ELSE IF ( delta_t /= 1.D0 .AND. delta_t >= 0 ) THEN
        !
        WRITE( stdout, '(/,5X,"Thermalization: delta_t = ",F6.3, &
                            & ", T = ",F6.1)' ) delta_t, temp_new * delta_t
        !
        CALL thermalize( temp_new, temp_new * delta_t )
        !
     END IF
     !
  END IF
  !
  WRITE( UNIT = stdout, &
         FMT = '(/,5X,"Entering Dynamics:",T28,"iteration",T37," = ",I5,/, &
                & T28,"time",T37," =  ",F8.5," pico-seconds")' ) &
      istep, elapsed_time
  !
  ! ... here starts the molecular dynamics :
  !
  ! ... calculate accelerations in a.u. units / alat
  !
  FORALL( na = 1:nat ) acc(:,na) = force(:,na) / mass(na) / alat
  !
  ! ... Damped dynamics ( based on the quick-min algorithm )
  !
  IF ( ldamped ) CALL project_velocity()
  !
  ! ... the old positions are updated to reflect the new velocities
  ! ... NB: vel is actually dt*( tau - tau_old )
  !
  tau_old = tau - vel
  !
  ! ... atoms are moved accordingly to the classical equation of motion.
  ! ... Verlet integration scheme.
  !
  tau_new = 2.D0 * tau - tau_old + dt**2 * acc
  !
  IF ( lconstrain ) THEN
     !
     ! ... check if the new positions satisfy the constrain equation
     !
     CALL check_constraint( nat, tau_new, tau, &
                            force, if_pos, ityp, alat, dt, amconv )
     !
     WRITE( stdout, '(/5X,"Corrected atomic positions:")')
     !
     CALL output_tau( .FALSE. )
     !
     WRITE( stdout, '(/,5X,"Constrained forces (Ry/au):",/)')
     !
     DO na = 1, nat
        !
        WRITE( stdout, '(5X,"atom ",I3," type ",I2,3X,"force = ",3F14.8)' ) &
            na, ityp(na), force(:,na)
        !
     END DO
     !
  END IF
  !
  ! ... check if convergence for structural minimization is achieved
  !
  IF ( ldamped ) THEN
     !
     conv_ions = ( etotold - etot ) < epse
     conv_ions = conv_ions .AND. ( MAXVAL( ABS( force ) ) < epsf )
     !
     IF ( conv_ions ) THEN
        !
        WRITE( UNIT = stdout, &
               FMT = '(/,5X,"Damped Dynamics: convergence achieved in ",I3, &
                          & " steps")' ) istep
        WRITE( UNIT = stdout, &
               FMT = '(/,5X,"End of damped dynamics calculation")' )
        WRITE( UNIT = stdout, &
               FMT = '(/,5X,"Final energy = ",F18.10," ryd"/)' ) etot                 
        !
        CALL output_tau( .TRUE. )
        !
        RETURN
        !
     END IF     
     !
  ELSE
     !
     IF ( istep == nstep ) THEN
        !
        conv_ions = .TRUE.
        !
        WRITE( UNIT = stdout, &
               FMT = '(/,5X,"End of molecular dynamics calculation")' )
        !
        RETURN
        !
     END IF
     !   
  END IF
  !
  ! ... the linear momentum and the kinetic energy are computed here
  !
  IF ( istep == 1 ) THEN
     !
     vel = ( tau_new - tau ) / dt
     !
  ELSE
     !
     vel = ( tau_new - tau_old ) / ( 2.D0 * dt )
     !
  END IF
  !
  ml   = 0.D0
  ekin = 0.D0  
  !
  DO na = 1, nat 
     ! 
     ml(:) = ml(:) + vel(:,na) * mass(na)
     ekin  = ekin + &
             0.5D0 * mass(na) * ( vel(1,na)**2 + vel(2,na)**2 + vel(3,na)**2 )
     !
  END DO  
  !
  ekin = ekin * alat**2
  !
  ! ... find the new temperature
  !
  temp_new = 2.D0 / 3.D0 * ekin / nat * convert_E_to_temp
  !
  ! ... save on file all the needed quantities
  !
  CALL seqopn( 4, 'md', 'FORMATTED',  exst )
  !
  WRITE( UNIT = 4, FMT = * ) &
      etot, temp_new, mass, total_mass, elapsed_time, istep, tau
  !
  CLOSE( UNIT = 4, STATUS = 'KEEP' )
  !
  ! ... here the tau are shifted
  !
  tau(:,:) = tau_new(:,:)
  !
  ! ... infos are written on the standard output
  !
  CALL output_tau( .FALSE. )
  !
  IF ( istep == 1 ) THEN
     !
     WRITE( stdout, '(5X,"kinetic energy (Ekin) = ",F14.8," Ry",/, &
                    & 5X,"temperature           = ",F14.8," K", /, &
                    & 5X,"Ekin + Etot (const)   = ",F14.8," Ry")' ) &
         temperature * 3.D0 / 2.D0 * nat / convert_E_to_temp, &
         temperature, &
         temperature * 3.D0 / 2.D0 * nat / convert_E_to_temp + etot
     !
  ELSE
     !        
     WRITE( stdout, '(5X,"kinetic energy (Ekin) = ",F14.8," Ry",/, &
                    & 5X,"temperature           = ",F14.8," K", /, &
                    & 5X,"Ekin + Etot (const)   = ",F14.8," Ry")' ) &
         ekin, temp_new, ( ekin  + etot )
      !
  END IF    
  !
  ! ... total linear momentum must be zero if all atoms move
  !
  mlt = norm( ml(:) )
  !
  IF ( ( mlt > eps8 ) .AND. &
       .NOT. ( ldamped .OR. lconstrain .OR. lfixatom ) ) &
     CALL infomsg ( 'dynamics', 'Total linear momentum <> 0', -1 )
  !
  WRITE( stdout, '(/,5X,"Linear momentum :",3(2X,F14.10))' ) ml
  !
  DEALLOCATE( mass )
  DEALLOCATE( tau_old )
  DEALLOCATE( tau_new )
  DEALLOCATE( vel )         
  DEALLOCATE( acc )
  !
  RETURN
  !
  CONTAINS
     !
     ! ... internal procedure
     !
     !-----------------------------------------------------------------------
     SUBROUTINE project_velocity()
       !-----------------------------------------------------------------------
       !
       USE constants, ONLY : eps32
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: norm_acc, acc_versor(3,nat)
       !
       ! ... external functions
       !
       REAL(DP), EXTERNAL :: DNRM2, DDOT
       !
       !
       norm_acc = DNRM2( 3*nat, acc, 1 )
       !
       IF ( norm_acc > eps32 ) THEN
          !
          acc_versor = acc / norm_acc
          !
          vel = acc_versor * &
                MAX( 0.D0, DDOT( 3*nat, vel, 1, acc_versor, 1 ) )
          !
       ELSE
          !
          vel = 0.D0
          !
       END IF       
       !
     END SUBROUTINE project_velocity 
     !
     !-----------------------------------------------------------------------
     SUBROUTINE start_therm()
       !-----------------------------------------------------------------------
       !
       ! ... Starting thermalization of the system
       !
       USE symme, ONLY : invsym, nsym, irt
       !
       IMPLICIT NONE
       !
       INTEGER  :: na, nb
       REAL(DP) :: total_mass, aux, velox, ek, &
                   ml(3), dir_x, dir_y, dir_z, module
       !  
       REAL(DP), EXTERNAL :: rndm
       !
       !    
       aux = temperature / convert_E_to_temp
       !
       ! ... velocity in random direction, with modulus accordingly to mass 
       ! ... and temperature: 3/2KT = 1/2mv^2
       !
       DO na = 1, nat
          !
          ! ... N.B. velox is in a.u. units /alat
          !
          velox = SQRT( 3.D0 * aux / mass(na) ) / alat
          !
          dir_x = rndm() - 0.5D0
          dir_y = rndm() - 0.5D0
          dir_z = rndm() - 0.5D0
          !
          module = 1.D0 / SQRT( dir_x**2 + dir_y**2 + dir_z**2 )
          !
          vel(1,na) = velox * dir_x * module
          vel(2,na) = velox * dir_y * module
          vel(3,na) = velox * dir_z * module
          !
       END DO
       !
       ! ... the velocity of fixed ions must be zero
       !
       vel = vel * DBLE( if_pos )
       !
       ! ... if there is inversion symmetry, equivalent atoms have 
       ! ... opposite velocities
       !
       ml(:) = 0.D0
       !
       IF ( invsym ) THEN
          !
          DO na = 1, nat
             !
             nb = irt( ( nsym / 2 + 1 ), na )
             !
             IF ( nb > na ) vel(:,nb) = - vel(:,na)
             !
             ! ... the atom on the inversion center is kept fixed
             !
             IF ( na == nb ) vel(:,na) = 0.D0
             !
          END DO
          !
       ELSE
          !
          ! ... put total linear momentum equal zero if all atoms move
          !
          IF ( .NOT. lfixatom ) THEN
             !
             total_mass = 0.D0
             !
             DO na = 1, nat
                !
                total_mass = total_mass + mass(na)
                !
                ml(:) = ml(:) + mass(na) * vel(:,na)
                !
             END DO
             !
             ml(:) = ml(:) / total_mass
             !
          END IF
          !
       END IF
       !
       ek = 0.D0
       !
       DO na = 1, nat
          !
          vel(:,na) = vel(:,na) - ml(:)
          !
          ek = ek + 0.5D0 * mass(na) * ( ( vel(1,na) - ml(1) )**2 + &
                                         ( vel(2,na) - ml(2) )**2 + &
                                         ( vel(3,na) - ml(3) )**2 )
          !
       END DO   
       !
       ! ... after the velocity of the center of mass has been subtracted the
       ! ... temperature is usually changed. Set again the temperature to the
       ! ... right value.
       !
       temp_new = 2.D0 * ek / ( 3.D0 * nat ) * alat**2 * convert_E_to_temp
       !
       CALL thermalize( temp_new, temperature )
       !
       RETURN
       !
     END SUBROUTINE start_therm 
     !
     !-----------------------------------------------------------------------
     SUBROUTINE thermalize( system_temp, required_temp )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: system_temp, required_temp
       !
       REAL(DP) :: aux
       !
       !
       ! ... rescale the velocities by a factor 3 / 2KT / Ek
       !
       IF ( system_temp > 0.D0 .AND. required_temp > 0.D0 ) THEN
          !
          aux = SQRT( required_temp / system_temp )
          !
       ELSE
          !
          aux = 0.D0
          !
       END IF
       !
       vel = vel * aux
       !
       RETURN
       !
     END SUBROUTINE thermalize         
     !  
END SUBROUTINE dynamics
