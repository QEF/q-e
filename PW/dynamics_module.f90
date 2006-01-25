!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#define __BFGS
!
!----------------------------------------------------------------------------
MODULE dynamics_module
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : amass
  USE io_global, ONLY : stdout
  USE io_files,  ONLY : prefix, tmp_dir
  USE constants, ONLY : tpi, fpi
  USE constants, ONLY : amconv, convert_E_to_temp, au_ps, bohr_radius_cm
  USE constants, ONLY : eps8, eps16
  !
  USE basic_algebra_routines
  !
  IMPLICIT NONE
  !
  SAVE
  !
  REAL(DP) :: &
       dt,            &! time step
       temperature,   &! starting temperature
       delta_t         ! rate of thermalization
  INTEGER :: &
       nraise          ! the frequency of temperature raising
  !
  REAL(DP), ALLOCATABLE :: tau_old(:,:), tau_new(:,:), tau_ref(:,:)
  REAL(DP), ALLOCATABLE :: vel(:,:), acc(:,:)
  REAL(DP), ALLOCATABLE :: mass(:)
  REAL(DP), ALLOCATABLE :: diff_coeff(:)
  !
  INTEGER, PARAMETER :: hist_len = 1000
  REAL(DP)           :: pair_distr(1:hist_len)
  !
  CONTAINS
    !
    ! ... public methods
    !
    !------------------------------------------------------------------------
    SUBROUTINE allocate_dyn_vars()
      !------------------------------------------------------------------------
      !
      USE ions_base, ONLY : nat
      !
      ALLOCATE( mass( nat ) )
      !
      ALLOCATE( tau_old( 3, nat ) )
      ALLOCATE( tau_new( 3, nat ) )
      ALLOCATE( tau_ref( 3, nat ) )
      !
      ALLOCATE( vel( 3, nat ) )
      ALLOCATE( acc( 3, nat ) )
      !
      ALLOCATE( diff_coeff( nat ) )
      !
      RETURN
      !
    END SUBROUTINE allocate_dyn_vars
    !
    !------------------------------------------------------------------------
    SUBROUTINE deallocate_dyn_vars()
      !------------------------------------------------------------------------
      !
      IF ( ALLOCATED( mass ) )        DEALLOCATE( mass )
      IF ( ALLOCATED( tau_old ) )     DEALLOCATE( tau_old )
      IF ( ALLOCATED( tau_new ) )     DEALLOCATE( tau_new )
      IF ( ALLOCATED( tau_ref ) )     DEALLOCATE( tau_ref )
      IF ( ALLOCATED( vel )  )        DEALLOCATE( vel )         
      IF ( ALLOCATED( acc )  )        DEALLOCATE( acc )
      IF ( ALLOCATED( diff_coeff ) )  DEALLOCATE( diff_coeff )
      !
      RETURN
      !
    END SUBROUTINE deallocate_dyn_vars
    !
    !------------------------------------------------------------------------
    SUBROUTINE dynamics()
      !------------------------------------------------------------------------
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
      ! ...              The initial velocity distribution is therefore a 
      ! ...              constant.
      !
      ! ... delta_t, nraise are used to change the temperature as follows:
      !
      ! ... delta_t = 1 :                   every 'nraise' step the actual 
      ! ...                                 temperature is rescaled to the
      ! ...                                 initial value.
      ! ... delta_t /= 1 and delta_t > 0 :  at each step the actual temperature
      ! ...                                 is multiplied by delta_t; this is
      ! ...                                 done rescaling all the velocities.
      ! ... delta_t < 0 :                   every 'nraise' step the temperature
      ! ...                                 reduced by -delta_t.
      !
      ! ... Dario Alfe' 1997  and  Carlo Sbraccia 2004-2005
      !
      USE ions_base,     ONLY : nat, nsp, ityp, tau, if_pos, atm
      USE cell_base,     ONLY : alat
      USE ener,          ONLY : etot
      USE force_mod,     ONLY : force
      USE relax,         ONLY : epse, epsf
      USE control_flags, ONLY : istep, nstep, conv_ions, &
                                lconstrain, ldamped, lfixatom, lrescale_t
      !
      USE constraints_module
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      REAL(DP) :: ekin, etotold
      REAL(DP) :: total_mass, temp_new, elapsed_time
      REAL(DP) :: ml(3), mlt
      INTEGER  :: i, na
      LOGICAL  :: file_exists
      !
      REAL(DP), EXTERNAL :: DNRM2
      !
      tau_old(:,:) = tau(:,:)
      tau_new(:,:) = 0.D0
      vel(:,:)     = 0.D0
      acc(:,:)     = 0.D0
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
      CALL seqopn( 4, 'md', 'FORMATTED', file_exists )
      !
      IF ( file_exists ) THEN
         !
         ! ... the file is read
         !
         READ( UNIT = 4, FMT = * ) &
            etotold, temp_new, mass(:), total_mass, &
            elapsed_time, istep, tau_old(:,:), tau_ref(:,:)
         !
         CLOSE( UNIT = 4, STATUS = 'KEEP' )
         !
      ELSE
         !
         CLOSE( UNIT = 4, STATUS = 'DELETE' )
         !
         ! ... atoms are refold in the central box
         !
         CALL refold_tau()
         !
         ! ... reference positions
         !
         tau_ref(:,:) = tau(:,:)
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
                         & " femto-seconds")' ) dt, ( dt * 2.D+3 * au_ps )
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
            temp_new = temperature
            !
            ! ... the old positions are updated to reflect the initial 
            ! ... velocities ( notice that vel is not the real velocity,
            ! ... but just a displacement vector )
            !
            tau_old(:,:) = tau(:,:) - vel(:,:)
            !
         ELSE
            !
            tau_old(:,:) = tau(:,:)
            !
         END IF
         !
         elapsed_time = 0.D0
         !
         istep = 0
         !
      END IF
      !
      ! ... elapsed_time is in picoseconds
      !
      elapsed_time = elapsed_time + dt * 2.D0 * au_ps
      !
      istep = istep + 1
      !
      IF ( lrescale_t ) THEN
         !
         vel(:,:) = tau(:,:) - tau_old(:,:)
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
                          & ", T = ",F6.1)' )  delta_t, ( temp_new + delta_t )
               !
               CALL thermalize( temp_new, ( temp_new + delta_t ) )
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
         ! ... the old positions are updated to reflect the new velocities
         ! ... ( notice that vel is not the real velocity, but just a 
         ! ... displacement vector )
         !
         tau_old(:,:) = tau(:,:) - vel(:,:)
         !
      END IF
      !
      WRITE( UNIT = stdout, &
             FMT = '(/,5X,"Entering Dynamics:",T28,"iteration",T37," = ",&
                    &I5,/,T28,"time",T37," = ",F8.4," pico-seconds")' ) &
          istep, elapsed_time
      !
      IF ( lconstrain ) THEN
         !
         ! ... we first remove the component of the force along the constraint
         ! ... gradient (this constitutes the initial guess for the lagrange
         ! ... multipliers)
         !
         CALL remove_constr_force( nat, tau, if_pos, ityp, alat, force )
         !
      END IF
      !
      IF ( ldamped ) THEN
         !
         ! ... Damped dynamics ( based on the quick-min algorithm )
         !
         vel(:,:) = tau(:,:) - tau_old(:,:)
         !
         acc(:,:) = force(:,:) / alat
         !
         CALL force_precond( acc, etotold )
         !
         acc(:,:) = acc(:,:) / amconv
         !
         CALL project_velocity()
         !
         ! ... the old positions are updated to reflect the new velocities
         !
         tau_old(:,:) = tau(:,:) - vel(:,:)
         !
      ELSE
         !
         ! ... calculate accelerations in a.u. units / alat
         !
         FORALL( na = 1:nat ) acc(:,na) = force(:,na) / mass(na) / alat
         !
      END IF
      !
      ! ... atoms are moved accordingly to the classical equation of motion.
      ! ... Verlet integration scheme.
      !
      tau_new(:,:) = 2.D0 * tau(:,:) - tau_old(:,:) + dt**2 * acc(:,:)
      !
      IF ( lconstrain ) THEN
         !
         ! ... check if the new positions satisfy the constrain equation
         !
         CALL check_constraint( nat, tau_new, tau, &
                                force, if_pos, ityp, alat, dt, amconv )
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
         WRITE( stdout, '(/5X,"Total force = ",F12.6)') DNRM2( 3*nat, force, 1 )
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
            CALL output_tau( .TRUE. )
            !
            CALL print_averages()
            !
            RETURN
            !
         END IF
         !   
      END IF
      !
      ! ... the linear momentum and the kinetic energy are computed here
      !
      IF ( istep > 1 .OR. lrescale_t ) &
         vel = ( tau_new - tau_old ) / ( 2.D0 * dt ) * DBLE( if_pos )
      !
      ml   = 0.D0
      ekin = 0.D0  
      !
      DO na = 1, nat 
         ! 
         ml(:) = ml(:) + vel(:,na) * mass(na)
         ekin  = ekin + 0.5D0 * mass(na) * &
                        ( vel(1,na)**2 + vel(2,na)**2 + vel(3,na)**2 )
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
      CALL seqopn( 4, 'md', 'FORMATTED',  file_exists )
      !
      WRITE( UNIT = 4, FMT = * ) &
          etot, temp_new, mass(:), total_mass, &
          elapsed_time, istep, tau(:,:), tau_ref(:,:)
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
      IF ( .NOT. ldamped ) THEN
         !
         WRITE( stdout, '(5X,"kinetic energy (Ekin) = ",F14.8," Ry",/,  &
                        & 5X,"temperature           = ",F14.8," K ",/,  &
                        & 5X,"Ekin + Etot (const)   = ",F14.8," Ry")' ) &
             ekin, temp_new, ( ekin  + etot )
         !
         ! ... total linear momentum must be zero if all atoms move
         !
         mlt = norm( ml(:) )
         !
         IF ( mlt > eps8 .AND. .NOT.( lconstrain .OR. lfixatom ) ) &
            CALL infomsg ( 'dynamics', 'Total linear momentum <> 0', -1 )
         !
         WRITE( stdout, '(/,5X,"Linear momentum :",3(2X,F14.10))' ) ml
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE dynamics
    !
    !-----------------------------------------------------------------------
    SUBROUTINE refold_tau()
      !-----------------------------------------------------------------------
      !
      USE ions_base,          ONLY : nat, tau
      USE cell_base,          ONLY : alat
      USE constraints_module, ONLY : pbc
      !
      IMPLICIT NONE
      !
      INTEGER :: ia
      !
      !
      DO ia = 1, nat
         !
         tau(:,ia) = pbc( tau(:,ia) * alat ) / alat
         !
      END DO
      !
    END SUBROUTINE refold_tau
    !
    !-----------------------------------------------------------------------
    SUBROUTINE compute_averages( istep )
      !-----------------------------------------------------------------------
      !
      USE ions_base,          ONLY : nat, tau
      USE cell_base,          ONLY : alat, at
      USE constraints_module, ONLY : pbc
      USE parser,             ONLY : delete_if_present
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: istep
      !
      INTEGER               :: i, j, index
      REAL(DP)              :: dx, dy, dz, distsq
      REAL(DP)              :: dtau(3)
      REAL(DP)              :: dmax
      REAL(DP), ALLOCATABLE :: rmsd(:)
      REAL(DP), PARAMETER   :: one_sixth = 1.D0 / 6.D0
      !
      ! ... RMSD and diffusion coefficient
      !
      ALLOCATE( rmsd( nat ) )
      !
      IF ( istep == 1 ) THEN
         !
         diff_coeff(:) = 0.D0
         !
         pair_distr(:) = 0.D0
         !
         CALL delete_if_present( TRIM( tmp_dir ) // &
                               & TRIM( prefix ) // ".rmsd.dat" )
         !
      END IF
      !
      DO i = 1, nat
         !
         dx = ( tau(1,i) - tau_ref(1,i) ) * alat
         dy = ( tau(2,i) - tau_ref(2,i) ) * alat
         dz = ( tau(3,i) - tau_ref(3,i) ) * alat
         !
         distsq = dx*dx + dy*dy + dz*dz
         !
         rmsd(i) = distsq
         !
      END DO
      !
      diff_coeff(:) = diff_coeff(:) + one_sixth * rmsd(:)
      !
      OPEN( UNIT = 777, POSITION = 'APPEND', &
            FILE = TRIM( tmp_dir ) // TRIM( prefix ) // ".rmsd.dat" )
      !
      WRITE( 777, '(2(2X,F16.8))' ) &
          (istep * dt * 2.D0 * au_ps), SQRT( SUM( rmsd(:) ) / DBLE( nat - 1 ) )
      !
      CLOSE( UNIT = 777, STATUS = 'KEEP' )
      !
      DEALLOCATE( rmsd )
      !
      ! ... pair distribution function g(r)
      !
      dmax = norm( MATMUL( at(:,:), (/ 0.5D0, 0.5D0, 0.5D0 /) ) ) * alat
      !
      DO i = 1, nat
         !
         DO j = 1, nat
            !
            IF ( i == j ) CYCLE
            !
            dtau(:) = pbc( ( tau(:,i) - tau(:,j) ) * alat )
            !
            index = ANINT( norm( dtau(:) ) / dmax * DBLE( hist_len ) )
            !
            pair_distr(index) = pair_distr(index) + 1.D0
            !
         END DO
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE compute_averages
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_averages()
      !-----------------------------------------------------------------------
      !
      USE control_flags,      ONLY : nstep
      USE cell_base,          ONLY : omega, at, alat
      USE ions_base,          ONLY : nat
      !
      IMPLICIT NONE
      !
      INTEGER  :: i, index
      REAL(DP) :: dist, dmax
      !
      ! ... diffusion coefficient
      !
      WRITE( UNIT = stdout, &
             FMT = '(/,5X,"diffusion coefficients :")' )
      !
      diff_coeff(:) = diff_coeff(:) / ( DBLE( nstep ) * dt )
      !
      diff_coeff(:) = diff_coeff(:) * bohr_radius_cm**2 / ( 2.D-12 * au_ps )
      !
      DO i = 1, nat
          !
          WRITE( UNIT = stdout, &
                 FMT = '(5X,"atom ",I5,"   D = ",F16.8," cm^2/s")' ) &
              i, diff_coeff(i)
          !
      END DO
      !
      WRITE( UNIT = stdout, FMT = '(/,5X,"< D > = ",F16.8," cm^2/s")' ) &
          SUM( diff_coeff(:) ) / DBLE( nat - 1 )
      !
      ! ... pair distribution function g(r)
      !
      dmax = norm( MATMUL( at(:,:), (/ 0.5D0, 0.5D0, 0.5D0 /) ) ) * alat
      !
      pair_distr(:) = pair_distr(:) * omega / DBLE( nat*nat ) / fpi
      !
      pair_distr(:) = pair_distr(:) / ( dmax / DBLE( hist_len ) )
      !
      pair_distr(:) = pair_distr(:) / DBLE( nstep )
      !
      OPEN( UNIT = 4, FILE = TRIM( tmp_dir ) // TRIM( prefix ) // ".pd.dat" )
      !
      DO index = 1, hist_len
         !
         dist = DBLE( index ) / DBLE( hist_len ) * dmax
         !
         WRITE( 4, '(2(2X,F16.8))' ) dist, pair_distr(index) / dist**2
         !
      END DO
      !
      CLOSE( UNIT = 4 )
      !
      RETURN
      !
    END SUBROUTINE print_averages
    !
    ! ... private methods
    !
    !-----------------------------------------------------------------------
    SUBROUTINE force_precond( force, etotold )
      !-----------------------------------------------------------------------
      !
      ! ... this routine computes an estimate of H^-1 by using the BFGS
      ! ... algorithm and the preconditioned gradient  pg = H^-1 * g
      ! ... ( it works in units of alat )
      !
      USE ener,      ONLY : etot
      USE cell_base, ONLY : alat
      USE ions_base, ONLY : nat, tau
      USE io_files,  ONLY : iunbfgs, iunbroy, tmp_dir
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(INOUT) :: force(:,:)
      REAL(DP), INTENT(IN)    :: etotold
      !
#if defined (__BFGS)
      !
      REAL(DP), ALLOCATABLE :: pos(:), pos_p(:)
      REAL(DP), ALLOCATABLE :: grad(:), grad_p(:), precond_grad(:)
      REAL(DP), ALLOCATABLE :: inv_hess(:,:)
      REAL(DP), ALLOCATABLE :: y(:), s(:)
      REAL(DP), ALLOCATABLE :: Hs(:), Hy(:), yH(:)
      REAL(DP)              :: sdoty, pg_norm
      INTEGER               :: dim
      CHARACTER(LEN=256)    :: bfgs_file
      LOGICAL               :: file_exists
      REAL(DP), PARAMETER   :: max_step = 0.8D0  ! in bohr
      !
      !
      dim = 3 * nat
      !
      ALLOCATE( pos( dim ), pos_p( dim ) )
      ALLOCATE( grad( dim ), grad_p( dim ), precond_grad( dim ) )
      ALLOCATE( y( dim ), s( dim ) )
      ALLOCATE( inv_hess( dim, dim ) )
      ALLOCATE( Hs( dim ), Hy( dim ), yH( dim ) )       
      !
      pos(:)  =   RESHAPE( tau,   (/ dim /) )
      grad(:) = - RESHAPE( force, (/ dim /) )
      !
      bfgs_file = TRIM( tmp_dir ) // TRIM( prefix ) // '.bfgs'
      !
      INQUIRE( FILE = TRIM( bfgs_file ) , EXIST = file_exists )
      !
      IF ( file_exists ) THEN
         !
         OPEN( UNIT = iunbfgs, &
               FILE = TRIM( bfgs_file ), STATUS = 'OLD', ACTION = 'READ' )
         !
         READ( iunbfgs, * ) pos_p
         READ( iunbfgs, * ) grad_p
         READ( iunbfgs, * ) inv_hess
         !
         CLOSE( UNIT = iunbfgs )
         !
         IF ( etot < etotold ) THEN
            !
            ! ... BFGS update
            !
            s(:) = pos(:)  - pos_p(:)
            y(:) = grad(:) - grad_p(:)
            !
            sdoty = ( s(:) .dot. y(:) )
            !
            IF ( sdoty > eps8 ) THEN
               !
               Hs(:) = ( inv_hess(:,:) .times. s(:) )
               Hy(:) = ( inv_hess(:,:) .times. y(:) )
               yH(:) = ( y(:) .times. inv_hess(:,:) )
               !
               inv_hess = inv_hess + 1.D0 / sdoty * &
                        ( ( 1.D0 + ( y .dot. Hy ) / sdoty ) * matrix( s, s ) - &
                          ( matrix( s, yH ) +  matrix( Hy, s ) ) )
               !
            END IF
            !
         END IF
         !
      ELSE
         !
         inv_hess(:,:) = identity( dim )
         !
      END IF
      !
      precond_grad(:) = ( inv_hess(:,:) .times. grad(:) )
      !
      IF ( ( precond_grad(:) .dot. grad(:) ) < 0.D0 ) THEN
         !
         WRITE( UNIT = stdout, &
                FMT = '(5X,/,"uphill step: resetting bfgs history",/)' )
         !
         precond_grad(:) = grad(:)
         !
         inv_hess(:,:) = identity( dim )
         !
      END IF
      !
      OPEN( UNIT = iunbfgs, &
            FILE = TRIM( bfgs_file ), STATUS = 'UNKNOWN', ACTION = 'WRITE' )
      !
      WRITE( iunbfgs, * ) pos(:)
      WRITE( iunbfgs, * ) grad(:)
      WRITE( iunbfgs, * ) inv_hess(:,:)
      !
      CLOSE( UNIT = iunbfgs )
      !
      ! ... the length of the step is always shorter than pg_norm
      !
      pg_norm = norm( precond_grad(:) )
      !
      precond_grad(:) = precond_grad(:) / pg_norm * &
                        MIN( max_step / alat, pg_norm )
      !
      force(:,:) = - RESHAPE( precond_grad(:), (/ 3, nat /) )
      !
      DEALLOCATE( pos, pos_p )
      DEALLOCATE( grad, grad_p, precond_grad )
      DEALLOCATE( inv_hess )
      DEALLOCATE( y, s )
      DEALLOCATE( Hs, Hy, yH )
      !
#endif
      !
      RETURN
      !
    END SUBROUTINE force_precond
    !
    !-----------------------------------------------------------------------
    SUBROUTINE project_velocity()
      !-----------------------------------------------------------------------
      !
      ! ... quick-min algorithm
      !
      USE control_flags, ONLY : istep
      USE ions_base,     ONLY : nat
      !
      IMPLICIT NONE
      !
      REAL(DP)              :: norm_acc, projection
      REAL(DP), ALLOCATABLE :: acc_versor(:,:)
      !
      REAL(DP), EXTERNAL :: DNRM2, DDOT
      !
      !
      IF ( istep == 1 ) RETURN
      !
      ALLOCATE( acc_versor( 3, nat ) )
      !
      norm_acc = DNRM2( 3*nat, acc(:,:), 1 )
      !
      acc_versor(:,:) = acc(:,:) / norm_acc
      !
      projection = DDOT( 3*nat, vel(:,:), 1, acc_versor(:,:), 1 )
      !
      WRITE( UNIT = stdout, FMT = '(/,5X,"<vel(dt)|acc(dt)> = ",F12.8)' ) &
          projection / DNRM2( 3*nat, vel, 1 )
      !
      vel(:,:) = acc_versor(:,:) * MAX( 0.D0, projection )
      !
      DEALLOCATE( acc_versor )
      !
      RETURN
      !
    END SUBROUTINE project_velocity 
    !
    !-----------------------------------------------------------------------
    SUBROUTINE start_therm()
      !-----------------------------------------------------------------------
      !
      ! ... Starting thermalization of the system
      !
      USE symme,          ONLY : invsym, nsym, irt
      USE control_flags,  ONLY : lfixatom
      USE cell_base,      ONLY : alat
      USE ions_base,      ONLY : nat, if_pos
      USE random_numbers, ONLY : gauss_dist
      !
      IMPLICIT NONE
      !
      INTEGER  :: na, nb
      REAL(DP) :: total_mass, KT, sigma, coeff, ek, ml(3), system_temp
      INTEGER  :: hist(-1000:1000), index
      !
      !
      KT = temperature / convert_E_to_temp
      !
      ! ... starting velocities have a Maxwell-Boltzmann distribution
      !
      DO na = 1, nat
         !
         coeff = ( mass(na) / ( tpi * KT) )**( 3.D0 / 2.D0 )
         !
         sigma = SQRT( KT / mass(na) )
         !
         ! ... N.B. velocities must in a.u. units of alat
         !
         vel(:,na) = coeff * gauss_dist( 0.D0, sigma, 3 ) / alat
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
      system_temp = 2.D0 * ek / ( 3.D0 * nat ) * alat**2 * convert_E_to_temp
      !
      CALL thermalize( system_temp, temperature )
      !
      ! ... vel is used already multiplied by the time step
      !
      vel(:,:) = vel(:,:) * dt
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
      vel(:,:) = vel(:,:) * aux
      !
      RETURN
      !
    END SUBROUTINE thermalize         
    !  
END MODULE dynamics_module
