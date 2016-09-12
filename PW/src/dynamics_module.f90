!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! TB
! included test if relaxz=.true. to allow for movement of the center of mass
! search for 'TB' 
!----------------------------------------------------------------------------
!
#undef __NPT
#if defined (__NPT)
#define RELAXTIME 2000.D0
#define TARGPRESS 2.39D0
#endif
!
!----------------------------------------------------------------------------
MODULE dynamics_module
   !----------------------------------------------------------------------------
   !
   USE kinds,     ONLY : DP
   USE ions_base, ONLY : amass
   USE io_global, ONLY : stdout
   USE io_files,  ONLY : prefix, tmp_dir, seqopn
   USE constants, ONLY : tpi, fpi
   USE constants, ONLY : amu_ry, ry_to_kelvin, au_ps, bohr_radius_cm, ry_kbar
   USE constants, ONLY : eps8
   USE control_flags, ONLY : tolp
   !
   USE basic_algebra_routines
   !
   IMPLICIT NONE
   !
   SAVE
   PRIVATE
   PUBLIC :: verlet, proj_verlet, terminate_verlet, &
             langevin_md, smart_MC, allocate_dyn_vars, deallocate_dyn_vars
   PUBLIC :: temperature, refold_pos, vel
   PUBLIC :: dt, delta_t, nraise, control_temp, thermostat
   !
   REAL(DP) :: &
         dt,          &! time step
         temperature, &! starting temperature
         virial,      &! virial (used for the pressure)
         delta_t       ! parameter used in thermalization
   INTEGER :: &
         nraise,      &! parameter used in thermalization
         ndof,        &! the number of degrees of freedom
         num_accept=0  ! Number of the accepted proposal in Smart_MC
   LOGICAL :: &
         vel_defined,  &! if true, vel is used rather than tau_old to do the next step
         control_temp, &! if true a thermostat is used to control the temperature
         refold_pos,   &! if true the positions are refolded into the supercell
         first_iter=.true. ! if this is the first ionic iteration
   CHARACTER(len=10) &
         thermostat    ! the thermostat used to control the temperature
   ! tau_smart and force_smart is used for smart Monte Carlo to store the atomic position of the
   ! previous step.
   REAL(DP), ALLOCATABLE :: tau_smart(:,:), force_smart(:,:)
   real(dp) :: etot_smart
   REAL(DP), ALLOCATABLE :: tau_old(:,:), tau_new(:,:), tau_ref(:,:)
   REAL(DP), ALLOCATABLE :: vel(:,:), acc(:,:), chi(:,:)
   REAL(DP), ALLOCATABLE :: mass(:)
   REAL(DP), ALLOCATABLE :: diff_coeff(:)
   REAL(DP), ALLOCATABLE :: radial_distr(:,:)
   !
   INTEGER, PARAMETER :: hist_len = 1000
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
      IF ( .not.allocated( mass ) ) ALLOCATE( mass( nat ) )
      !
      IF ( .not.allocated( tau_old ) ) ALLOCATE( tau_old( 3, nat ) )
      IF ( .not.allocated( tau_new ) ) ALLOCATE( tau_new( 3, nat ) )
      IF ( .not.allocated( tau_ref ) ) ALLOCATE( tau_ref( 3, nat ) )
      !
      IF ( .not.allocated( vel ) ) ALLOCATE( vel( 3, nat ) )
      IF ( .not.allocated( acc ) ) ALLOCATE( acc( 3, nat ) )
      IF ( .not.allocated( chi ) ) ALLOCATE( chi( 3, nat ) )
      !
      IF ( .not.allocated( diff_coeff ) ) ALLOCATE( diff_coeff( nat ) )
      !
      IF ( .not.allocated( radial_distr ) ) &
         ALLOCATE( radial_distr( hist_len , nat ) )
      !
   END SUBROUTINE allocate_dyn_vars
   !
   !------------------------------------------------------------------------
   SUBROUTINE deallocate_dyn_vars()
      !------------------------------------------------------------------------
      !
      IF ( allocated( mass ) )          DEALLOCATE( mass )
      IF ( allocated( tau_old ) )       DEALLOCATE( tau_old )
      IF ( allocated( tau_new ) )       DEALLOCATE( tau_new )
      IF ( allocated( tau_ref ) )       DEALLOCATE( tau_ref )
      IF ( allocated( vel )  )          DEALLOCATE( vel )
      IF ( allocated( acc )  )          DEALLOCATE( acc )
      IF ( allocated( chi )  )          DEALLOCATE( chi )
      IF ( allocated( diff_coeff ) )    DEALLOCATE( diff_coeff )
      IF ( allocated( radial_distr ) )  DEALLOCATE( radial_distr )
      !
   END SUBROUTINE deallocate_dyn_vars
    !
    !------------------------------------------------------------------------
   SUBROUTINE verlet()
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
      ! ... Dario Alfe' 1997  and  Carlo Sbraccia 2004-2006
      !
      USE ions_base,      ONLY : nat, nsp, ityp, tau, if_pos, atm
      USE cell_base,      ONLY : alat, omega
      USE ener,           ONLY : etot
      USE force_mod,      ONLY : force, lstres
      USE control_flags,  ONLY : istep, lconstrain, tv0rd
      !
      USE constraints_module, ONLY : nconstr, check_constraint
      USE constraints_module, ONLY : remove_constr_force, remove_constr_vec
      !
      IMPLICIT NONE
      !
      REAL(DP) :: ekin, etotold
      REAL(DP) :: total_mass, temp_new, temp_av, elapsed_time
      REAL(DP) :: delta(3), ml(3), mlt
      INTEGER  :: na
      ! istep counts all MD steps, including those of previous runs
#if defined (__NPT)
      REAL(DP) :: chi, press_new
#endif
      LOGICAL  :: file_exists, leof
      REAL(DP), EXTERNAL :: dnrm2
      REAL(DP) :: kstress(3,3)
      INTEGER :: i, j
      !
      ! ... the number of degrees of freedom
      !
      IF ( any( if_pos(:,:) == 0 ) ) THEN
         !
         ndof = 3*nat - count( if_pos(:,:) == 0 ) - nconstr
         !
      ELSE
         !
         ndof = 3*nat - 3 - nconstr
         !
      ENDIF
      !
      vel_defined  = .true.
      temp_av      = 0.D0
      !
      CALL seqopn( 4, 'md', 'FORMATTED', file_exists )
      !
      IF ( file_exists ) THEN
         !
         ! ... the file is read :  simulation is continuing
         !
         READ( UNIT = 4, FMT = * ) etotold, istep, tau_old(:,:), leof
         !
         IF ( leof ) THEN
            !
            ! ... the file was created by projected_verlet:  Ignore it
            !
            CALL md_init()
            !
         ELSE
            !
            vel_defined = .false.
            !
            READ( UNIT = 4, FMT = * ) &
               temp_new, temp_av, mass(:), total_mass, elapsed_time, &
               tau_ref(:,:)
            !
         ENDIF
         !
         CLOSE( UNIT = 4, STATUS = 'KEEP' )
         !
      ELSE
         !
         CLOSE( UNIT = 4, STATUS = 'DELETE' )
         !
         ! ... the file is absent :  simulation is starting from scratch
         !
         CALL md_init()
         !
      ENDIF
      !
      ! ... elapsed_time is in picoseconds
      !
      elapsed_time = elapsed_time + dt*2.D0*au_ps
      !
      istep = istep + 1
      !
      WRITE( UNIT = stdout, &
             FMT = '(/,5X,"Entering Dynamics:",T28,"iteration",T37," = ", &
                    &I5,/,T28,"time",T37," = ",F8.4," pico-seconds",/)' ) &
          istep, elapsed_time
      !
      IF ( control_temp ) CALL apply_thermostat()
      !
      ! ... we first remove the component of the force along the
      ! ... constraint gradient ( this constitutes the initial
      ! ... guess for the calculation of the lagrange multipliers )
      !
      IF ( lconstrain ) &
         CALL remove_constr_force( nat, tau, if_pos, ityp, alat, force )
      !
      ! ... calculate accelerations in a.u. units / alat
      !
      FORALL( na = 1:nat ) acc(:,na) = force(:,na) / mass(na) / alat
      !
      ! ... Verlet integration scheme
      !
      IF (vel_defined) THEN
         !
         ! ... remove the component of the velocity along the
         ! ... constraint gradient
         !
         IF ( lconstrain ) &
            CALL remove_constr_vec( nat, tau, if_pos, ityp, alat, vel )
         !
         tau_new(:,:) = tau(:,:) + vel(:,:) * dt + 0.5_DP * acc(:,:) * dt**2
         tau_old(:,:) = tau(:,:) - vel(:,:) * dt + 0.5_DP * acc(:,:) * dt**2
         !
      ELSE
         !
         tau_new(:,:) = 2.D0*tau(:,:) - tau_old(:,:) + acc(:,:) * dt**2
         !
      ENDIF
      !
      IF ( .not. any( if_pos(:,:) == 0 ) ) THEN
         !
         ! ... if no atom has been fixed  we compute the displacement of the
         ! ... center of mass and we subtract it from the displaced positions
         !
         delta(:) = 0.D0
         DO na = 1, nat
            delta(:) = delta(:) + mass(na)*( tau_new(:,na) - tau(:,na) )
         ENDDO
         delta(:) = delta(:) / total_mass
         FORALL( na = 1:nat ) tau_new(:,na) = tau_new(:,na) - delta(:)
         !
         IF (vel_defined) THEN
            delta(:) = 0.D0
            DO na = 1, nat
               delta(:) = delta(:) + mass(na)*( tau_old(:,na) - tau(:,na) )
            ENDDO
            delta(:) = delta(:) / total_mass
            FORALL( na = 1:nat ) tau_old(:,na) = tau_old(:,na) - delta(:)
         ENDIF
         !
      ENDIF
      !
      IF ( lconstrain ) THEN
         !
         ! ... check if the new positions satisfy the constrain equation
         !
         CALL check_constraint( nat, tau_new, tau, &
                                force, if_pos, ityp, alat, dt, amu_ry )
         !
#if ! defined (__REDUCE_OUTPUT)
         !
         WRITE( stdout, '(/,5X,"Constrained forces (Ry/au):",/)')
         !
         DO na = 1, nat
            !
            WRITE( stdout, &
                   '(5X,"atom ",I3," type ",I2,3X,"force = ",3F14.8)' ) &
                na, ityp(na), force(:,na)
            !
         ENDDO
         !
         WRITE( stdout, '(/5X,"Total force = ",F12.6)') dnrm2( 3*nat, force, 1 )
         !
#endif

         IF (vel_defined) THEN
            CALL check_constraint( nat, tau_old, tau, &
                                   force, if_pos, ityp, alat, dt, amu_ry )
         ENDIF
         !
      ENDIF
      !
      ! ... the linear momentum and the kinetic energy are computed here
      !
      vel = ( tau_new - tau_old ) / ( 2.D0*dt ) * dble( if_pos )
      !
      ml   = 0.D0
      ekin = 0.D0
      kstress = 0.d0
      !
      DO na = 1, nat
         !
         ml(:) = ml(:) + vel(:,na) * mass(na)
         ekin  = ekin + 0.5D0 * mass(na) * &
                        ( vel(1,na)**2 + vel(2,na)**2 + vel(3,na)**2 )
         DO i = 1, 3
             DO j = 1, 3
                 kstress(i,j) = kstress(i,j) + mass(na)*vel(i,na)*vel(j,na)
             ENDDO
         ENDDO
         !
      ENDDO
      !
      ekin = ekin*alat**2
      kstress = kstress * alat**2 / omega
      !
      ! ... find the new temperature and update the average
      !
      temp_new = 2.D0 / dble( ndof ) * ekin * ry_to_kelvin
      !
      temp_av = temp_av + temp_new
      !
#if defined (__NPT)
      !
      ! ... find the new pressure (in Kbar)
      !
      press_new = ry_kbar*( nat*temp_new/ry_to_kelvin + virial ) / omega
      !
      chi = 1.D0 - dt / RELAXTIME*( TARGPRESS - press_new )
      !
      omega = chi * omega
      alat  = chi**(1.D0/3.D0) * alat
      !
      WRITE( stdout, '(/,5X,"NEW ALAT = ",F8.5,2X,"Bohr"  )' ) alat
      WRITE( stdout, '(  5X,"PRESSURE = ",F8.5,2X,"Kbar",/)' ) press_new
      !
#endif
      !
      ! ... save all the needed quantities on file
      !
      CALL seqopn( 4, 'md', 'FORMATTED',  file_exists )
      !
      leof = .false.
      WRITE( UNIT = 4, FMT = * ) etot, istep, tau(:,:), leof
      !
      WRITE( UNIT = 4, FMT = * ) &
          temp_new, temp_av, mass(:), total_mass, elapsed_time, tau_ref(:,:)
      !
      CLOSE( UNIT = 4, STATUS = 'KEEP' )
      !
      ! ... here the tau are shifted
      !
      tau(:,:) = tau_new(:,:)
      !
#if ! defined (__REDUCE_OUTPUT)
      !
      CALL output_tau( .false., .false. )
      !
#endif
      !
      ! ... infos are written on the standard output
      !
      WRITE( stdout, '(5X,"kinetic energy (Ekin) = ",F14.8," Ry",/,  &
                     & 5X,"temperature           = ",F14.8," K ",/,  &
                     & 5X,"Ekin + Etot (const)   = ",F14.8," Ry")' ) &
          ekin, temp_new, ( ekin  + etot )
      IF (lstres) WRITE ( stdout, &
      '(5X,"Ions kinetic stress = ",F10.2," (kbar)",/3(27X,3F10.2/)/)') &
              ((kstress(1,1)+kstress(2,2)+kstress(3,3))/3.d0*ry_kbar), &
              (kstress(i,1)*ry_kbar,kstress(i,2)*ry_kbar,kstress(i,3)*ry_kbar, i=1,3)
      !
      IF ( .not.( lconstrain .or. any( if_pos(:,:) == 0 ) ) ) THEN
         !
         ! ... total linear momentum must be zero if all atoms move
         !
         mlt = norm( ml(:) )
         !
         IF ( mlt > eps8 ) &
            CALL infomsg( 'dynamics', 'Total linear momentum <> 0' )
         !
         WRITE( stdout, '(/,5X,"Linear momentum :",3(2X,F14.10))' ) ml(:)
         !
      ENDIF
      !
      ! ... compute the average quantities
      !
      CALL compute_averages( istep )
      !
   CONTAINS
      !
      !--------------------------------------------------------------------
      SUBROUTINE md_init()
         !--------------------------------------------------------------------
         !
         IMPLICIT NONE
         !
         istep = 0
         !
         WRITE( UNIT = stdout, &
               FMT = '(/,5X,"Molecular Dynamics Calculation")' )
         !
         ! ... atoms are refold in the central box if required
         !
         IF ( refold_pos ) CALL refold_tau()
         !
         ! ... reference positions
         !
         tau_ref(:,:) = tau(:,:)
         !
         IF ( control_temp ) THEN
            !
            WRITE( stdout, &
                  '(/,5X,"Starting temperature",T27," = ",F8.2," K")' ) &
               temperature
            !
            SELECT CASE( trim( thermostat ) )
               !
            CASE( 'andersen', 'Andersen' )
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"temperature is controlled by Andersen ", &
                              &   "thermostat",/,5x,"Collision frequency =",&
                              &    f7.4,"/timestep")' ) 1.0_dp/nraise
               !
            CASE( 'berendsen', 'Berendsen' )
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"temperature is controlled by soft ", &
                            &     "(Berendsen) velocity rescaling",/,5x,&
                            &     "Characteristic time =",i3,"*timestep")') &
                               nraise
               !
            CASE( 'initial', 'Initial' )
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"temperature is set once at start"/)' )
               !
            CASE DEFAULT
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"temperature is controlled by ", &
                              &     "velocity rescaling (",A,")"/)' )&
                              trim( thermostat )
               !
            END SELECT
            !
         ENDIF
         !
         DO na = 1, nsp
            !
            WRITE( UNIT = stdout, &
                  FMT = '(5X,"mass ",A2,T27," = ",F8.2)' ) atm(na), amass(na)
            !
         ENDDO
         !
         WRITE( UNIT = stdout, &
               FMT = '(5X,"Time step",T27," = ",F8.2," a.u.,",F8.4, &
                        & " femto-seconds")' ) dt, dt*2.D+3*au_ps
         !
         ! ... masses in rydberg atomic units
         !
         total_mass = 0.D0
         !
         DO na = 1, nat
            !
            mass(na) = amass( ityp(na) ) * amu_ry
            !
            total_mass = total_mass + mass(na)
            !
         ENDDO
         !
         IF ( tv0rd ) THEN ! initial velocities available from input file
            !
            vel(:,:) = vel(:,:) / alat
            !
         ELSEIF ( control_temp ) THEN
            !
            ! ... initial thermalization. N.B. tau is in units of alat
            !
            CALL start_therm()
            vel_defined = .true.
            !
            temp_new = temperature
            !
            temp_av = 0.D0
            !
         ELSE
            !
            vel(:,:) = 0.0_DP
            vel_defined = .true.
            !
         ENDIF
         !
         elapsed_time = 0.D0
         !
      END SUBROUTINE md_init
      !
      !--------------------------------------------------------------------
      SUBROUTINE apply_thermostat()
         !--------------------------------------------------------------------
         !
         USE random_numbers, ONLY : randy, gauss_dist
         !
         IMPLICIT NONE
         !
         INTEGER :: nat_moved
         REAL(DP) :: sigma, kt
         !
         IF(.not.vel_defined)THEN
            vel(:,:) = (tau(:,:) - tau_old(:,:)) / dt
         ENDIF
         !
         SELECT CASE( trim( thermostat ) )
         CASE( 'rescaling' )
            IF ( abs (temp_new-temperature) > tolp ) THEN
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"Velocity rescaling: T (",F6.1,"K) ", &
                                 & "out of range, reset to " ,F6.1)' ) &
                           temp_new, temperature
               CALL thermalize( 0, temp_new, temperature )
               !
            ENDIF
         CASE( 'rescale-v', 'rescale-V', 'rescale_v', 'rescale_V' )
            IF ( mod( istep, nraise ) == 0 ) THEN
               !
               temp_av = temp_av / dble( nraise )
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"Velocity rescaling: average T on ",i3, &
                                 &" steps (",F6.1,"K) reset to ",F6.1)' )  &
                           nraise, temp_av, temperature
               !
               CALL thermalize( 0, temp_new, temperature )
               !
               temp_av = 0.D0
               !
            ENDIF
         CASE( 'rescale-T', 'rescale-t', 'rescale_T', 'rescale_t' )
            IF ( delta_t > 0 ) THEN
               !
               temperature = temp_new*delta_t
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"Thermalization: T (",F6.1,"K) rescaled ",&
                                 & "by a factor ",F6.3)' ) temp_new, delta_t
               !
               CALL thermalize( 0, temp_new, temperature )
               !
            ENDIF
         CASE( 'reduce-T', 'reduce-t', 'reduce_T', 'reduce_t' )
            IF ( mod( istep, nraise ) == 0 .and. delta_t < 0 ) THEN
               !
               temperature = temp_new + delta_t
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"Thermalization: T (",F6.1,"K) reduced ",&
                                 & "by ",F6.3)' ) temp_new, -delta_t
               !
               CALL thermalize( 0, temp_new, temperature )
               !
            ENDIF
            !
         CASE( 'berendsen', 'Berendsen' )
            !
            WRITE( UNIT = stdout, &
                FMT = '(/,5X,"Soft (Berendsen) velocity rescaling")' )
            !
            CALL thermalize( nraise, temp_new, temperature )
            !
         CASE( 'andersen', 'Andersen' )
            !
            kt = temperature / ry_to_kelvin
            nat_moved = 0
            !
            DO na = 1, nat
               !
               IF ( randy() < 1.D0 / dble( nraise ) ) THEN
                  !
                  nat_moved = nat_moved + 1
                  sigma = sqrt( kt / mass(na) )
                  !
                  ! ... N.B. velocities must in a.u. units of alat and are zero
                  ! ...      for fixed ions
                  !
                  vel(:,na) = dble( if_pos(:,na) ) * &
                              gauss_dist( 0.D0, sigma, 3 ) / alat
                  !
               ENDIF
               !
            ENDDO
            !
            IF ( nat_moved > 0) WRITE( UNIT = stdout, &
               FMT = '(/,5X,"Andersen thermostat: ",I4," collisions")' ) &
                     nat_moved
            !
         CASE( 'initial', 'Initial' )
            !
            CONTINUE
            !
         END SELECT
         !
         ! ... the old positions are updated to reflect the new velocities
         !
         IF(.not.vel_defined)THEN
            tau_old(:,:) = tau(:,:) - vel(:,:) * dt
         ENDIF
         !
      END SUBROUTINE apply_thermostat
      !
      !-----------------------------------------------------------------------
      SUBROUTINE start_therm()
         !-----------------------------------------------------------------------
         !
         ! ... Starting thermalization of the system
         !
         USE symm_base,      ONLY : invsym, nsym, irt
         USE cell_base,      ONLY : alat
         USE ions_base,      ONLY : nat, if_pos
         USE random_numbers, ONLY : gauss_dist, set_random_seed
         !
         IMPLICIT NONE
         !
         INTEGER  :: na, nb
         REAL(DP) :: total_mass, kt, sigma, ek, ml(3), system_temp
         !
         ! ... next command prevents different MD runs to start
         ! ... with exactly the same "random" velocities
         !
         CALL set_random_seed ( )
         kt = temperature / ry_to_kelvin
         !
         ! ... starting velocities have a Maxwell-Boltzmann distribution
         !
         DO na = 1, nat
            !
            sigma = sqrt( kt / mass(na) )
            !
            ! ... N.B. velocities must in a.u. units of alat
            !
            vel(:,na) = gauss_dist( 0.D0, sigma, 3 ) / alat
            !
         ENDDO
         !
         ! ... the velocity of fixed ions must be zero
         !
         vel = vel * dble( if_pos )
         !
         IF ( lconstrain ) THEN
            !
            ! ... remove the component of the velocity along the
            ! ... constraint gradient
            !
            CALL remove_constr_vec( nat, tau, if_pos, ityp, alat, vel )
            !
         ENDIF
         !
         IF ( invsym ) THEN
            !
            ! ... if there is inversion symmetry, equivalent atoms have
            ! ... opposite velocities
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
            ENDDO
            !
         ELSE
            !
            ! ... put total linear momentum equal zero if all atoms
            ! ... are free to move
            !
            ml(:) = 0.D0
            !
            IF ( .not. any( if_pos(:,:) == 0 ) ) THEN
               !
               total_mass = sum ( mass(1:nat) )
               DO na = 1, nat
                  ml(:) = ml(:) + mass(na)*vel(:,na)
               ENDDO
               ml(:) = ml(:) / total_mass
               !
            ENDIF
            !
         ENDIF
         !
         ek = 0.D0
         !
         DO na = 1, nat
            !
            vel(:,na) = vel(:,na) - ml(:)
            !
            ek = ek + 0.5D0 * mass(na) * &
                     ( ( vel(1,na) )**2 + ( vel(2,na) )**2 + ( vel(3,na) )**2 )
            !
         ENDDO
         !
         ! ... after the velocity of the center of mass has been subtracted the
         ! ... temperature is usually changed. Set again the temperature to the
         ! ... right value.
         !
         system_temp = 2.D0 / dble( ndof ) * ek * alat**2 * ry_to_kelvin
         !
         CALL thermalize( 0, system_temp, temperature )
         !
      END SUBROUTINE start_therm
      !
   END SUBROUTINE verlet
   !
   !------------------------------------------------------------------------
   SUBROUTINE terminate_verlet
     !------------------------------------------------------------------------
     !
     USE io_global, ONLY : stdout
     !
     WRITE( UNIT = stdout, &
          FMT = '(/,5X,"The maximum number of steps has been reached.")' )
     WRITE( UNIT = stdout, &
          FMT = '(/,5X,"End of molecular dynamics calculation")' )
     !
     CALL print_averages()
     !
   END SUBROUTINE terminate_verlet
   !
   !------------------------------------------------------------------------
   SUBROUTINE proj_verlet( conv_ions )
      !------------------------------------------------------------------------
      !
      ! ... This routine performs one step of structural relaxation using
      ! ... the preconditioned-projected-Verlet algorithm.
      !
      USE ions_base,     ONLY : nat, ityp, tau, if_pos
      USE cell_base,     ONLY : alat
      USE ener,          ONLY : etot
      USE force_mod,     ONLY : force
      USE relax,         ONLY : epse, epsf
      USE control_flags, ONLY : istep, lconstrain
      !
      USE constraints_module, ONLY : remove_constr_force, check_constraint
      !
      ! TB
      USE extfield,      ONLY : relaxz
      !
      IMPLICIT NONE
      LOGICAL, INTENT(OUT) :: conv_ions
      !
      REAL(DP), ALLOCATABLE :: step(:,:)
      REAL(DP)              :: norm_step, etotold, delta(3)
      INTEGER               :: na
      LOGICAL               :: file_exists,leof
      !
      REAL(DP), PARAMETER :: step_max = 0.6D0  ! bohr
      !
      REAL(DP), EXTERNAL :: dnrm2
      !
      !
      ALLOCATE( step( 3, nat ) )
      !
      tau_old(:,:) = tau(:,:)
      tau_new(:,:) = 0.D0
      vel(:,:)     = 0.D0
      acc(:,:)     = 0.D0
      conv_ions = .FALSE.
      !
      CALL seqopn( 4, 'md', 'FORMATTED', file_exists )
      !
      IF ( file_exists ) THEN
         !
         ! ... the file is read
         !
         READ( UNIT = 4, FMT = * ) etotold, istep, tau_old(:,:)
         !
         CLOSE( UNIT = 4, STATUS = 'KEEP' )
         !
      ELSE
         !
         CLOSE( UNIT = 4, STATUS = 'DELETE' )
         !
         ! ... atoms are refold in the central box
         !
         IF ( refold_pos ) CALL refold_tau()
         !
         tau_old(:,:) = tau(:,:)
         etotold = etot
         istep = 0
         WRITE( UNIT = stdout, &
                FMT = '(/,5X,"Damped Dynamics Calculation")' )
         !
      ENDIF
      !
      IF ( lconstrain ) THEN
         !
         ! ... we first remove the component of the force along the
         ! ... constraint gradient (this constitutes the initial guess
         ! ... for the calculation of the lagrange multipliers)
         !
         CALL remove_constr_force( nat, tau, if_pos, ityp, alat, force )
         !
#if ! defined (__REDUCE_OUTPUT)
         !
         WRITE( stdout, '(/,5X,"Constrained forces (Ry/au):",/)')
         !
         DO na = 1, nat
            !
            WRITE( stdout, &
                  '(5X,"atom ",I3," type ",I2,3X,"force = ",3F14.8)' ) &
               na, ityp(na), force(:,na)
            !
         ENDDO
         !
         WRITE( stdout, &
               '(/5X,"Total force = ",F12.6)') dnrm2( 3*nat, force, 1 )
         !
#endif
         !
      ENDIF
      !
      ! ... check if convergence for structural minimization is achieved
      !
      conv_ions = ( etotold - etot ) < epse
      conv_ions = conv_ions .and. ( maxval( abs( force ) ) < epsf )
      !
      IF ( conv_ions ) THEN
         !
         WRITE( UNIT = stdout, &
                FMT = '(/,5X,"Damped Dynamics: convergence achieved in " &
                       & ,I3," steps")' ) istep
         WRITE( UNIT = stdout, &
                FMT = '(/,5X,"End of damped dynamics calculation")' )
         WRITE( UNIT = stdout, &
                FMT = '(/,5X,"Final energy = ",F18.10," Ry"/)' ) etot
         !
         CALL output_tau( .true., .true. )
         !
         RETURN
         !
      ENDIF
      !
      istep = istep + 1
      WRITE( stdout, '(/,5X,"Entering Dynamics:",&
                      & T28,"iteration",T37," = ",I5)' ) istep
      !
      ! ... Damped dynamics ( based on the projected-Verlet algorithm )
      !
      vel(:,:) = tau(:,:) - tau_old(:,:)
      !
      CALL force_precond( istep, force, etotold )
      !
      acc(:,:) = force(:,:) / alat / amu_ry
      !
      CALL project_velocity()
      !
      step(:,:) = vel(:,:) + dt**2 * acc(:,:)
      !
      norm_step = dnrm2( 3*nat, step, 1 )
      !
      step(:,:) = step(:,:) / norm_step
      !
      tau_new(:,:) = tau(:,:) + step(:,:)*min( norm_step, step_max / alat )
      !
      ! TB
      !IF ( .not. any( if_pos(:,:) == 0 ) ) THEN
      IF ( .not. any( if_pos(:,:) == 0 ) .AND. (relaxz) ) THEN
         WRITE( stdout, '("relaxz = .true. => displacement of the center of mass is not subtracted")')
      ENDIF
      IF ( (.not. any( if_pos(:,:) == 0 )) .AND. (.not. relaxz) ) THEN
         !
         ! ... if no atom has been fixed  we compute the displacement of the
         ! ... center of mass and we subtract it from the displaced positions
         !
         delta(:) = 0.D0
         !
         DO na = 1, nat
            !
            delta(:) = delta(:) + ( tau_new(:,na) - tau(:,na) )
            !
         ENDDO
         !
         delta(:) = delta(:) / dble( nat )
         !
         FORALL( na = 1:nat ) tau_new(:,na) = tau_new(:,na) - delta(:)
         !
      ENDIF
      !
      IF ( lconstrain ) THEN
         !
         ! ... check if the new positions satisfy the constrain equation
         !
         CALL check_constraint( nat, tau_new, tau, &
                                force, if_pos, ityp, alat, dt, amu_ry )
         !
      ENDIF
      !
      ! ... save on file all the needed quantities
      !
      CALL seqopn( 4, 'md', 'FORMATTED',  file_exists )
      !
      leof = .true.
      WRITE( UNIT = 4, FMT = * ) etot, istep, tau(:,:), leof
      !
      CLOSE( UNIT = 4, STATUS = 'KEEP' )
      !
      ! ... here the tau are shifted
      !
      tau(:,:) = tau_new(:,:)
      !
#if ! defined (__REDUCE_OUTPUT)
      !
      CALL output_tau( .false., .false. )
      !
#endif
      !
      DEALLOCATE( step )
      !
   END SUBROUTINE proj_verlet
   !------------------------------------------------------------------------
   SUBROUTINE langevin_md()
      !------------------------------------------------------------------------
      !
      ! ...
      !
      USE ions_base,      ONLY : nat, ityp, tau, if_pos
      USE cell_base,      ONLY : alat
      USE ener,           ONLY : etot
      USE force_mod,      ONLY : force
      USE control_flags,  ONLY : istep, lconstrain
      USE random_numbers, ONLY : gauss_dist
      !
      USE constraints_module, ONLY : nconstr
      USE constraints_module, ONLY : remove_constr_force, check_constraint
      !
      IMPLICIT NONE
      !
      REAL(DP) :: sigma, kt
      REAL(DP) :: delta(3)
      INTEGER  :: na
      LOGICAL  :: file_exists
      !
      REAL(DP), EXTERNAL :: dnrm2
      !
      CALL seqopn( 4, 'md', 'FORMATTED', file_exists )
      !
      IF ( file_exists ) THEN
         !
         ! ... the file is read :  simulation is continuing
         !
         READ( UNIT = 4, FMT = * ) istep
         !
         CLOSE( UNIT = 4, STATUS = 'KEEP' )
         !
      ELSE
         !
         CLOSE( UNIT = 4, STATUS = 'DELETE' )
         !
         ! ... the file is absent :  simulation is starting from scratch
         !
         istep = 0
         !
         WRITE( UNIT = stdout, &
               FMT = '(/,5X,"Over-damped Langevin Dynamics Calculation")' )
         !
         ! ... atoms are refold in the central box if required
         !
         IF ( refold_pos ) CALL refold_tau()
         !
         WRITE( UNIT = stdout, &
                FMT = '(5X,"Integration step",T27," = ",F8.2," a.u.,")' ) dt
         !
      ENDIF
      !
      istep = istep + 1
      WRITE( UNIT = stdout, &
             FMT = '(/,5X,"Entering Dynamics:",T28, &
                    &     "iteration",T37," = ",I5,/)' ) istep
      !
      IF ( lconstrain ) THEN
         !
         ! ... we first remove the component of the force along the
         ! ... constraint gradient ( this constitutes the initial
         ! ... guess for the calculation of the lagrange multipliers )
         !
         CALL remove_constr_force( nat, tau, if_pos, ityp, alat, force )
         !
      ENDIF
      !
      ! ... compute the stochastic term
      !
      kt = temperature / ry_to_kelvin
      !
      sigma = sqrt( 2.D0*dt*kt )
      !
      delta(:) = 0.D0
      !
      DO na = 1, nat
         !
         chi(:,na) = gauss_dist( 0.D0, sigma, 3 )*dble( if_pos(:,na) )
         !
         delta(:) = delta(:) + chi(:,na)
         !
      ENDDO
      !
      FORALL( na = 1:nat ) chi(:,na) = chi(:,na) - delta(:) / dble( nat )
      !
      PRINT *, "|F|   = ", dt*dnrm2( 3*nat, force, 1 )
      PRINT *, "|CHI| = ", dnrm2( 3*nat, chi, 1 )
      !
      ! ... over-damped Langevin dynamics
      !
      tau_new(:,:) = tau(:,:) + ( dt*force(:,:) + chi(:,:) ) / alat
      !
      IF ( .not. any( if_pos(:,:) == 0 ) ) THEN
         !
         ! ... here we compute the displacement of the center of mass and we
         ! ... subtract it from the displaced positions
         !
         delta(:) = 0.D0
         !
         DO na = 1, nat
            !
            delta(:) = delta(:) + ( tau_new(:,na) - tau(:,na) )
            !
         ENDDO
         !
         FORALL( na = 1:nat ) tau_new(:,na) = tau_new(:,na) - delta(:)
         !
      ENDIF
      !
      IF ( lconstrain ) THEN
         !
         ! ... check if the new positions satisfy the constrain equation
         !
         CALL check_constraint( nat, tau_new, tau, &
                                force, if_pos, ityp, alat, dt, amu_ry )
         !
#if ! defined (__REDUCE_OUTPUT)
         !
         WRITE( stdout, '(/,5X,"Constrained forces (Ry/au):",/)')
         !
         DO na = 1, nat
            !
            WRITE( stdout, &
                   '(5X,"atom ",I3," type ",I2,3X,"force = ",3F14.8)' ) &
                na, ityp(na), force(:,na)
            !
         ENDDO
         !
         WRITE( stdout, '(/5X,"Total force = ",F12.6)') dnrm2( 3*nat, force, 1 )
         !
#endif
         !
      ENDIF
      !
      ! ... save all the needed quantities on file
      !
      CALL seqopn( 4, 'md', 'FORMATTED',  file_exists )
      !
      WRITE( UNIT = 4, FMT = * ) istep
      !
      CLOSE( UNIT = 4, STATUS = 'KEEP' )
      !
      ! ... here the tau are shifted
      !
      tau(:,:) = tau_new(:,:)
      !
#if ! defined (__REDUCE_OUTPUT)
      !
      CALL output_tau( .false., .false. )
      !
#endif
      !
   END SUBROUTINE langevin_md
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
      ENDDO
      !
   END SUBROUTINE refold_tau
   !
   !-----------------------------------------------------------------------
   SUBROUTINE compute_averages( istep )
      !-----------------------------------------------------------------------
      !
      USE ions_base,          ONLY : nat, tau, fixatom
      USE cell_base,          ONLY : alat, at
      USE constraints_module, ONLY : pbc
      USE io_files,           ONLY : delete_if_present
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: istep
      !
      INTEGER               :: i, j, idx
      REAL(DP)              :: dx, dy, dz
      REAL(DP)              :: dtau(3)
      REAL(DP)              :: inv_dmax
      REAL(DP), ALLOCATABLE :: msd(:)
      REAL(DP), PARAMETER   :: max_dist(3) = (/ 0.5D0, 0.5D0, 0.5D0 /)
      !
      ! ... MSD and diffusion coefficient
      !
      ALLOCATE( msd( nat ) )
      !
      IF ( istep == 1 ) THEN
         !
         radial_distr(:,:) = 0.D0
         !
         CALL delete_if_present( trim( tmp_dir ) // &
                               & trim( prefix ) // ".msd.dat" )
         !
      ENDIF
      !
      DO i = 1, nat
         !
         dx = ( tau(1,i) - tau_ref(1,i) ) * alat
         dy = ( tau(2,i) - tau_ref(2,i) ) * alat
         dz = ( tau(3,i) - tau_ref(3,i) ) * alat
         !
         msd(i) = dx*dx + dy*dy + dz*dz
         !
      ENDDO
      !
      diff_coeff(:) = msd(:) / ( 6.D0*dble( istep )*dt )
      !
      ! ... conversion from Rydberg atomic units to cm^2/sec
      !
      diff_coeff(:) = diff_coeff(:) * bohr_radius_cm**2 / ( 2.D-12*au_ps )
      !
      OPEN( UNIT = 4, POSITION = 'APPEND', &
            FILE = trim( tmp_dir ) // trim( prefix ) // ".msd.dat" )
      !
      WRITE( 4, '(2(2X,F16.8))' ) &
          ( istep*dt*2.D0*au_ps ), sum( msd(:) ) / dble( nat-fixatom )
      !
      CLOSE( UNIT = 4, STATUS = 'KEEP' )
      !
      DEALLOCATE( msd )
      !
      ! ... radial distribution function g(r)
      !
      inv_dmax = 1.D0 / ( norm( matmul( at(:,:), max_dist(:) ) ) * alat )
      !
      DO i = 1, nat
         !
         DO j = 1, nat
            !
            IF ( i == j ) CYCLE
            !
            dtau(:) = pbc( ( tau(:,i) - tau(:,j) ) * alat )
            !
            idx = anint( norm( dtau(:) ) * inv_dmax * dble( hist_len ) )
            !
            IF( idx > 0 .and. idx <= size( radial_distr, 1 ) ) &
               radial_distr(idx,i) = radial_distr(idx,i) + 1.D0
            !
         ENDDO
         !
      ENDDO
      !
   END SUBROUTINE compute_averages
   !
   !-----------------------------------------------------------------------
   SUBROUTINE print_averages()
      !-----------------------------------------------------------------------
      !
      USE control_flags, ONLY : nstep
      USE cell_base,     ONLY : omega, at, alat
      USE ions_base,     ONLY : nat, fixatom
      !
      IMPLICIT NONE
      !
      INTEGER             :: i, idx
      REAL(DP)            :: dist, dmax
      REAL(DP), PARAMETER :: max_dist(3) = (/ 0.5D0, 0.5D0, 0.5D0 /)
      !
      ! ... diffusion coefficient
      !
      WRITE( UNIT = stdout, &
             FMT = '(/,5X,"diffusion coefficients :")' )
      !
      DO i = 1, nat
          !
          WRITE( UNIT = stdout, &
                 FMT = '(5X,"atom ",I5,"   D = ",F16.8," cm^2/s")' ) &
              i, diff_coeff(i)
          !
      ENDDO
      !
      WRITE( UNIT = stdout, FMT = '(/,5X,"< D > = ",F16.8," cm^2/s")' ) &
          sum( diff_coeff(:) ) / dble( nat-fixatom )
      !
      ! ... radial distribution function g(r)
      !
      dmax = norm( matmul( at(:,:), max_dist(:) ) ) * alat
      !
      radial_distr(:,:) = radial_distr(:,:) * omega / dble( nat ) / fpi
      !
      radial_distr(:,:) = radial_distr(:,:) / ( dmax / dble( hist_len ) )
      !
      radial_distr(:,:) = radial_distr(:,:) / dble( nstep )
      !
      OPEN( UNIT = 4, FILE = trim( tmp_dir ) // trim( prefix ) // ".rdf.dat" )
      !
      DO idx = 1, hist_len
         !
         dist = dble( idx ) / dble( hist_len ) * dmax
         !
         IF ( dist > dmax / sqrt( 3.0d0 ) ) CYCLE
         !
         radial_distr(idx,:) = radial_distr(idx,:) / dist**2
         !
         WRITE( 4, '(2(2X,F16.8))' ) &
             dist, sum( radial_distr(idx,:) ) / dble( nat )
         !
      ENDDO
      !
      CLOSE( UNIT = 4 )
      !
   END SUBROUTINE print_averages
   !
   !-----------------------------------------------------------------------
   SUBROUTINE force_precond( istep, force, etotold )
      !-----------------------------------------------------------------------
      !
      ! ... this routine computes an estimate of H^-1 by using the BFGS
      ! ... algorithm and the preconditioned gradient  pg = H^-1 * g
      ! ... ( it works in atomic units )
      !
      USE ener,      ONLY : etot
      USE cell_base, ONLY : alat
      USE ions_base, ONLY : nat, tau
      USE io_files,  ONLY : iunbfgs, tmp_dir
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(in)    :: istep
      REAL(DP), INTENT(inout) :: force(:,:)
      REAL(DP), INTENT(in)    :: etotold
      !
      REAL(DP), ALLOCATABLE :: pos(:), pos_p(:)
      REAL(DP), ALLOCATABLE :: grad(:), grad_p(:), precond_grad(:)
      REAL(DP), ALLOCATABLE :: inv_hess(:,:)
      REAL(DP), ALLOCATABLE :: y(:), s(:)
      REAL(DP), ALLOCATABLE :: Hy(:), yH(:)
      REAL(DP)              :: sdoty, pg_norm
      INTEGER               :: dim
      CHARACTER(len=256)    :: bfgs_file
      LOGICAL               :: file_exists
      !
      INTEGER,  PARAMETER   :: nrefresh    = 25
      REAL(DP), PARAMETER   :: max_pg_norm = 0.8D0
      !
      !
      dim = 3 * nat
      !
      ALLOCATE( pos( dim ), pos_p( dim ) )
      ALLOCATE( grad( dim ), grad_p( dim ), precond_grad( dim ) )
      ALLOCATE( y( dim ), s( dim ) )
      ALLOCATE( inv_hess( dim, dim ) )
      ALLOCATE( Hy( dim ), yH( dim ) )
      !
      pos(:)  =   reshape( tau,   (/ dim /) ) * alat
      grad(:) = - reshape( force, (/ dim /) )
      !
      bfgs_file = trim( tmp_dir ) // trim( prefix ) // '.bfgs'
      !
      INQUIRE( FILE = trim( bfgs_file ) , EXIST = file_exists )
      !
      IF ( file_exists ) THEN
         !
         OPEN( UNIT = iunbfgs, &
               FILE = trim( bfgs_file ), STATUS = 'OLD', ACTION = 'READ' )
         !
         READ( iunbfgs, * ) pos_p
         READ( iunbfgs, * ) grad_p
         READ( iunbfgs, * ) inv_hess
         !
         CLOSE( UNIT = iunbfgs )
         !
         ! ... the approximate inverse hessian is reset to one every nrefresh
         ! ... iterations: this is one to clean-up the memory of the starting
         ! ... configuration
         !
         IF ( mod( istep, nrefresh ) == 0 ) inv_hess(:,:) = identity( dim )
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
               Hy(:) = ( inv_hess(:,:) .times. y(:) )
               yH(:) = ( y(:) .times. inv_hess(:,:) )
               !
               inv_hess = inv_hess + 1.D0 / sdoty * &
                        ( ( 1.D0 + ( y .dot. Hy ) / sdoty ) * matrix( s, s ) - &
                          ( matrix( s, yH ) +  matrix( Hy, s ) ) )
               !
            ENDIF
            !
         ENDIF
         !
      ELSE
         !
         inv_hess(:,:) = identity( dim )
         !
      ENDIF
      !
      precond_grad(:) = ( inv_hess(:,:) .times. grad(:) )
      !
      IF ( ( precond_grad(:) .dot. grad(:) ) < 0.D0 ) THEN
         !
         WRITE( UNIT = stdout, &
                FMT = '(/,5X,"uphill step: resetting bfgs history",/)' )
         !
         precond_grad(:) = grad(:)
         !
         inv_hess(:,:) = identity( dim )
         !
      ENDIF
      !
      OPEN( UNIT = iunbfgs, &
            FILE = trim( bfgs_file ), STATUS = 'UNKNOWN', ACTION = 'WRITE' )
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
      precond_grad(:) = precond_grad(:) / pg_norm
      precond_grad(:) = precond_grad(:) * min( pg_norm, max_pg_norm )
      !
      force(:,:) = - reshape( precond_grad(:), (/ 3, nat /) )
      !
      DEALLOCATE( pos, pos_p )
      DEALLOCATE( grad, grad_p, precond_grad )
      DEALLOCATE( inv_hess )
      DEALLOCATE( y, s )
      DEALLOCATE( Hy, yH )
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
      REAL(DP), EXTERNAL :: dnrm2, ddot
      !
      !
      IF ( istep == 1 ) RETURN
      !
      ALLOCATE( acc_versor( 3, nat ) )
      !
      norm_acc = dnrm2( 3*nat, acc(:,:), 1 )
      !
      acc_versor(:,:) = acc(:,:) / norm_acc
      !
      projection = ddot( 3*nat, vel(:,:), 1, acc_versor(:,:), 1 )
      !
      WRITE( UNIT = stdout, FMT = '(/,5X,"<vel(dt)|acc(dt)> = ",F12.8)' ) &
          projection / dnrm2( 3*nat, vel, 1 )
      !
      vel(:,:) = acc_versor(:,:) * max( 0.D0, projection )
      !
      DEALLOCATE( acc_versor )
      !
   END SUBROUTINE project_velocity
   !
   !-----------------------------------------------------------------------
   SUBROUTINE thermalize( nraise, system_temp, required_temp )
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(in) :: system_temp, required_temp
      INTEGER, INTENT(in) :: nraise
      !
      REAL(DP) :: aux
      !
      IF ( nraise > 0 ) THEN
         !
         ! ... Berendsen rescaling (Eq. 7.59 of Allen & Tildesley)
         ! ... the "rise time" is tau=nraise*dt so dt/tau=1/nraise
         ! ... Equivalent to traditional rescaling if nraise=1
         !
         IF ( system_temp > 0.D0 .and. required_temp > 0.D0 ) THEN
            !
            aux = sqrt( 1.d0 + (required_temp / system_temp - 1.d0) * &
                                (1.D0/dble (nraise) ) )
            !
         ELSE
            !
            aux = 0.d0
            !
         ENDIF
         !
      ELSE
         !
         ! ... rescale the velocities by a factor 3 / 2KT / Ek
         !
         IF ( system_temp > 0.D0 .and. required_temp > 0.D0 ) THEN
            !
            aux = sqrt( required_temp / system_temp )
            !
         ELSE
            !
            aux = 0.d0
            !
         ENDIF
         !
      ENDIF
      !
      vel(:,:) = vel(:,:) * aux
      !
   END SUBROUTINE thermalize
   !

   !-----------------------------------------------------------------------
   SUBROUTINE smart_MC()
      !-----------------------------------------------------------------------
      ! Routine to apply smart_MC
      ! Implemented by Xiaochuan Ge, Jul., 2013
      !
      ! At this moment works only with langevin dynamics !!
      ! For the formula see R.J.Rossky, JCP, 69, 4628(1978)

     USE ions_base,      ONLY : nat, ityp, tau, if_pos,atm
     USE cell_base,      ONLY : alat
     USE ener,           ONLY : etot
     USE force_mod,      ONLY : force
     USE control_flags,  ONLY : istep, lconstrain
     USE constraints_module, ONLY : remove_constr_force, check_constraint
     USE random_numbers, ONLY : randy
     USE io_files,      ONLY : prefix
     USE io_global,      ONLY : ionode
     USE constants, ONLY : bohr_radius_angs

     IMPLICIT NONE

     LOGICAL :: accept
     real(dp) :: kt,sigma2,&
                 T_ij,T_ji,boltzman_ji,& ! boltzman_ji=exp[-(etot_new-etot_old)/kt]
                 temp,p_smc                  ! *_smart means *_old, the quantity of the
                                        ! previous step

     INTEGER :: ia, ip

      IF ( lconstrain ) THEN
         ! ... we first remove the component of the force along the
         ! ... constraint gradient ( this constitutes the initial
         ! ... guess for the calculation of the lagrange multipliers )
         CALL remove_constr_force( nat, tau, if_pos, ityp, alat, force )
      ENDIF

     IF(first_iter) THEN ! For the first iteration
       ALLOCATE(tau_smart(3,nat))
       ALLOCATE(force_smart(3,nat))
       tau_smart=tau
       etot_smart=etot
       force_smart=force
       first_iter=.false.
       RETURN
     ENDIF

     kt = temperature / ry_to_kelvin
     sigma2 =  2.D0*dt*kt

     T_ij=0.0d0
     T_ji=0.0d0
     DO ia=1,nat
       DO ip = 1, 3
         T_ij=T_ij+((tau(ip,ia)-tau_smart(ip,ia))*alat-dt*force_smart(ip,ia))**2
         T_ji=T_ji+((tau_smart(ip,ia)-tau(ip,ia))*alat-dt*force(ip,ia))**2
       ENDDO
     ENDDO
     T_ij=exp(-T_ij/(2*sigma2))
     T_ji=exp(-T_ji/(2*sigma2))

     boltzman_ji=exp(-(etot-etot_smart)/kt)

     p_smc=T_ji*boltzman_ji/T_ij

     WRITE(stdout, '(5x,"The old energy is:",3x,F17.8," Ry")') etot_smart
     WRITE(stdout, '(5x,"The new energy is:",3x,F17.8," Ry")') etot
     WRITE(stdout, '(5x,"The possibility to accept this step is:",3x,F10.7/)') p_smc
     WRITE(stdout, '(5x,"Nervously waiting for the fate ..."/)')

     ! Decide if accept the new config
     temp = randy()
     WRITE(stdout, '(5x,"The fate says:",5x,F10.7)') temp
     IF(temp <= p_smc) THEN
       WRITE(stdout, '(5x,"The new config is accepted")')
       num_accept=num_accept+1
       tau_smart=tau
       etot_smart=etot
       force_smart=force
     ELSE
       WRITE(stdout, '(5x,"The new config is not accepted")')
       tau=tau_smart
       etot=etot_smart
       force=force_smart
     ENDIF

     WRITE (stdout, '(5x,"The current acceptance is :",3x,F10.6)') dble(num_accept)/istep

     ! Print the trajectory
#if defined(__MPI)
     IF(ionode) THEN
#endif
     OPEN(117,file="trajectory-"//trim(prefix)//".xyz",status="unknown",position='APPEND')
     WRITE(117,'(I5)') nat
     WRITE(117,'("# Step: ",I5,5x,"Total energy: ",F17.8,5x,"Ry")') istep-1, etot
     DO ia = 1, nat
       WRITE( 117, '(A3,3X,3F14.9)') atm(ityp(ia)),tau(:,ia)*alat*bohr_radius_angs
     ENDDO
     CLOSE(117)
#if defined(__MPI)
     ENDIF
#endif

     RETURN
   END SUBROUTINE smart_MC

END MODULE dynamics_module
