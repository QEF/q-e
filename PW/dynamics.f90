!
! Copyright (C) 2001-2004 PWSCF group
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
  ! ... This routine performs one step of molecular dynamics evolution using
  ! ... the Velocity Verlet algorithm. 
  !
  ! ... Parameters:
  ! ... mass         mass of the atoms
  ! ... dt           time step
  ! ... temperature  starting temperature
  ! ...              The starting velocities of atoms are set accordingly
  ! ...              to the starting temperature, in random directions.
  ! ...              The initial velocity distribution is therefore a constant
  !
  ! ... delta_T, nraise are used to change the temperature as follows:
  !
  ! ... delta_T = 1 :                  nothing is done.
  ! ... delta_T /= 1 and delta_T > 0 : at each step the actual temperature is
  ! ...                                multiplied by delta_T; this is done
  ! ...                                rescaling all the velocities
  ! ... delta_T < 0 :                  every 'nraise' step the temperature
  ! ...                                reduced by -delta_T
  !
  ! ... Dario Alfe 1997  and  Carlo Sbraccia 2004
  !
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  USE constants,     ONLY : amconv, eps8
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, atm
  USE cell_base,     ONLY : alat
  USE dynam,         ONLY : amass, temperature, dt, delta_t, nraise
  USE ener,          ONLY : etot
  USE force_mod,     ONLY : force
  USE klist,         ONLY : nelec
  USE ions_base,     ONLY : if_pos
  USE relax,         ONLY : epse, epsf
  USE control_flags, ONLY : alpha0, beta0, istep, nstep, conv_ions, &
                            lconstrain, ldamped, lfixatom
  USE io_files,      ONLY : prefix
  !
  USE basic_algebra_routines
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  REAL(KIND=DP), ALLOCATABLE :: vel(:,:), acc(:,:), accold(:,:)
    ! velocities, accelerations and accelerations of the previous step
  REAL(KIND=DP), ALLOCATABLE :: mass(:)    
    ! masses of atoms
  REAL(KIND=DP) :: ekin, etotold
    ! ionic kinetic energy 
    ! ionic potential energy (of the previous step)
  REAL(KIND=DP) :: total_mass, temp_new, elapsed_time, norm_of_dtau
  REAL(KIND=DP) :: ml(3), mlt
    ! total linear momentum and its modulus
  INTEGER :: i, na
    ! counters
  LOGICAL :: exst
  REAL(KIND=DP), PARAMETER :: convert_E_to_temp = 315642.28D0 * 0.5D0
  !
  !
  ALLOCATE( mass( nat ) )
  ALLOCATE( vel( 3, nat ) )         
  ALLOCATE( acc( 3, nat ) )
  ALLOCATE( accold( 3, nat ) )
  !
  vel    = 0.D0
  acc    = 0.D0
  accold = 0.D0
  !
  IF ( istep == 1 ) THEN
     IF ( ldamped ) THEN
        WRITE( UNIT = stdout, &
               FMT = '(/,5X,"Damped Dynamics Calculation")' )
     ELSE
        WRITE( UNIT = stdout, &
               FMT = '(/,5X,"Molecular Dynamics Calculation")' )
     END IF
  END IF
  !
  ! ... one Ryd a.u. of time is 4.84*10^-17 seconds, i.e. 0.0484  femtoseconds
  !
  CALL seqopn( 4, TRIM( prefix ) // '.md', 'FORMATTED', exst )
  !
  IF ( .NOT. exst ) THEN
     !
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
     WRITE( UNIT = stdout, &
            FMT = '(/,5X,"Starting temperature",T27," = ",F8.2," K")' ) &
         temperature
     !
     DO na = 1, ntyp
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
        mass(na)   = amass( ityp(na) ) * amconv
        total_mass = total_mass + mass(na)
        !
     END DO
     !
     ! ... initial thermalization. N.B. tau is in units of alat
     !
     CALL start_therm()
     !
     elapsed_time = 0.D0
     !
     temp_new = temperature
     !
     istep = 0
     !
  ELSE
     !
     READ( UNIT = 4, FMT = * ) &
        etotold, temp_new, mass, total_mass, &
        elapsed_time, istep, tau, vel, accold
     !
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     !
  END IF
  !
  elapsed_time = elapsed_time + dt * 0.0000484D0
  !
  istep = istep + 1
  !
  IF ( MOD( istep, nraise ) == 0 .AND. delta_T < 0 ) THEN
     !
     WRITE( stdout, '(/,5X,"Thermalization: delta_T = ",F6.3, &
                         & ", T = ",F6.1)' )  - delta_T, ( temp_new - delta_T )
     !
     CALL thermalize( temp_new, ( temp_new - delta_T ) )
     !
  ELSE IF ( delta_T /= 1.D0 .AND. delta_T >= 0 ) THEN
     !
     WRITE( stdout, '(/,5X,"Thermalization: delta_T = ",F6.3, &
                         & ", T = ",F6.1)' ) delta_T, temp_new * delta_T
     !
     CALL thermalize( temp_new, temp_new * delta_T )
     !
  END IF
  !
  ! ... check if convergence for structural minimization is achieved
  !
  IF ( ldamped ) THEN
     !
     conv_ions = ( etotold - etot ) < epse
     !
     DO i = 1, 3
        DO na = 1, nat
           conv_ions = conv_ions .AND. ( ABS( force(i,na) ) < epsf )
        END DO
     END DO
     !
     IF ( conv_ions ) THEN
        !
        WRITE( UNIT = stdout, &
               FMT = '(/,5X,"Damped Dynamics: convergence achieved in ",I3, &
                          & " steps")' ) istep
        WRITE( UNIT = stdout, &
               FMT = '(/,5X,"End of damped dynamics calculation")' )
        WRITE( UNIT = stdout, &
               FMT = '(/,5X,"Efinal = ",F15.8,/)' ) etot                 
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
     END IF
     !   
  END IF
  !
  WRITE( UNIT = stdout, &
         FMT = '(/,5X,"Entering Dynamics:",T28,"iteration",T37," = ",I5,/, &
                & T28,"time",T37," =  ",F8.5," pico-seconds")' ) &
      istep, elapsed_time
  !
  ! ... calculate accelerations in a.u. units / alat
  !
  FORALL( na = 1 : nat ) &
     acc(:,na) = force(:,na) / mass(na) / alat
  !
  ! ... atoms are moved accordingly to the classical equation of motion.
  ! ... Damped dynamics ( based on the quick-min algorithm ) is also
  ! ... done here.
  !
  vel = vel + 0.5D0 * dt * ( accold + acc )
  !
  IF ( ldamped ) CALL project_velocity()
  !
  ! ... constrains ( atoms kept fixed ) are reinforced 
  !
  vel = vel * DBLE( if_pos )
  !
  tau = tau + dt * vel + 0.5D0 * dt**2 * acc
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
  ! ... find the new temperature
  !
  temp_new = 2.D0 / 3.D0 * ekin * alat**2 / nat * convert_E_to_temp
  !
  ! ... save on file needed quantity
  !
  CALL seqopn( 4, TRIM( prefix ) // '.md', 'FORMATTED',  exst )
  !
  WRITE( UNIT = 4, FMT = * ) &
      etot, temp_new, mass, total_mass, &
      elapsed_time, istep, tau, vel, acc
  !
  CLOSE( UNIT = 4, STATUS = 'KEEP' )
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
         ekin*alat**2, &
         temp_new, &
         ekin*alat**2 + etot
      !
  END IF    
  !
  ! ... total linear momentum must be zero if all atoms move
  !
  mlt = norm( ml(:) )
  !
  IF ( ( mlt > eps8 ) .AND. &
       .NOT. ( ldamped .OR. lconstrain .OR. lfixatom ) ) &
     CALL errore( 'dynamics', 'Total linear momentum <> 0', - 1 )
  !
  WRITE( stdout, '(/,5X,"Linear momentum :",3(2X,F14.10))' ) ml
  !
  DEALLOCATE( mass )
  DEALLOCATE( vel )         
  DEALLOCATE( acc )
  DEALLOCATE( accold )
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
       REAL (KIND=DP) :: norm_acc, acc_versor(3,nat)
       !
       ! ... external functions
       !
       REAL (KIND=DP), EXTERNAL :: DNRM2, DDOT
       !
       !
       norm_acc = DNRM2( 3*nat, acc, 1 )
       !
       IF ( norm_acc > eps32 ) THEN
          !
          acc_versor = acc / norm_acc
          !
          vel = MAX( 0.D0, DDOT( 3*nat, vel, 1, acc_versor, 1 ) ) * acc_versor
          !
       ELSE
          !
          vel = 0.D0
          !
       END IF       
       !
     END SUBROUTINE project_velocity 
     !
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
       ! ... local variables
       !
       INTEGER                :: na, nb
       REAL(KIND=DP)          :: total_mass, aux, velox, ek, &
                                 ml(3), dir_x, dir_y, dir_z, module
         ! ek = kinetic energy
       !  
       REAL(KIND=DP),EXTERNAL :: rndm
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
     !
     !-----------------------------------------------------------------------
     SUBROUTINE thermalize( temp_old, temp_new )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       ! ... INPUT variables
       !
       REAL(KIND=DP) :: temp_new, temp_old
       !
       ! ... local variables
       !
       INTEGER       :: na
       REAL(KIND=DP) :: aux
       !
       !
       ! ... rescale the velocities by a factor 3 / 2KT / Ek
       !
       IF ( temp_new > 0.D0 .AND. temp_old > 0.D0 ) THEN
          !
          aux = SQRT( temp_new / temp_old )
          !
       ELSE
          !
          aux = 0.D0
          !
       END IF
       !
       vel = vel * aux * DBLE( if_pos )
       !
       RETURN
       !
     END SUBROUTINE thermalize         
     !  
END SUBROUTINE dynamics
