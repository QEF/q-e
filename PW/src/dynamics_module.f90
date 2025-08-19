!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! TB
! included test if relaxz=.TRUE. to allow for movement of the center of mass
! search for 'TB' 
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
MODULE dynamics_module
   !----------------------------------------------------------------------------
   !! It includes all the quantities and procedures related to the molecular
   !! dynamics.
   !
   USE kinds,           ONLY : DP
   USE ions_base,       ONLY : amass
   USE io_global,       ONLY : stdout
   USE io_files,        ONLY : prefix, tmp_dir, seqopn
   USE constants,       ONLY : tpi, fpi
   USE constants,       ONLY : amu_ry, ry_to_kelvin, au_ps, bohr_radius_cm, &
                               ry_kbar, HARTREE_SI, RYDBERG_SI  
   USE constants,       ONLY : eps8
   USE control_flags,   ONLY : tolp
   !
   USE basic_algebra_routines
   !
   IMPLICIT NONE
   !
   SAVE
   PRIVATE
   PUBLIC :: verlet, verlet_read_tau_from_conf, proj_verlet, terminate_verlet, &
             fire, langevin_md, smart_MC, allocate_dyn_vars, deallocate_dyn_vars,& 
             get_ndof
   PUBLIC :: velocity_verlet 
   PUBLIC :: temperature, refold_pos, vel
   PUBLIC :: dt, delta_t, nraise, control_temp, thermostat, elapsed_time
   ! FIRE parameters
   PUBLIC :: fire_nmin, fire_f_inc, fire_f_dec, fire_alpha_init, fire_falpha, fire_dtmax
   !
   !
   REAL(DP) :: dt
   !! time step
   REAL(DP) :: temperature
   !! starting temperature
   REAL(DP) :: virial
   !! virial (used for the pressure)
   REAL(DP) :: delta_t
   !! parameter used in thermalization
   INTEGER :: nraise
   !! parameter used in thermalization
   INTEGER :: ndof
   !! the number of degrees of freedom
   INTEGER :: num_accept=0
   !! Number of the accepted proposal in Smart_MC
   LOGICAL :: vel_defined
   !! if true, vel is used rather than tau_old to do the next step
   LOGICAL :: control_temp
   !! if true a thermostat is used to control the temperature
   LOGICAL :: refold_pos
   !! if true the positions are refolded into the supercell
   LOGICAL :: first_iter=.TRUE.
   !! true if this is the first ionic iteration
   CHARACTER(len=10) :: thermostat
   !! the thermostat used to control the temperature
   INTEGER ::  fire_nmin
   !! FIRE: minimum number of steps for time step increase 
   REAL(DP) :: fire_f_inc
   !! FIRE: factor for time step increase  
   REAL(DP) :: fire_f_dec
   !! FIRE: factor for time step decrease
   REAL(DP) :: fire_alpha_init
   !! FIRE: initial value of mixing factor
   REAL(DP) :: fire_falpha
   !! FIRE: modify the mixing factor
   REAL(DP) :: fire_dtmax
   !! FIRE: factor for max time step 
   REAL(DP), ALLOCATABLE :: tau_smart(:,:)
   !! used for smart Monte Carlo to store the atomic position of the
   !! previous step.
   REAL(DP), ALLOCATABLE :: force_smart(:,:)
   !! used for smart Monte Carlo to store the force of the
   !! previous step.
   REAL(DP) :: etot_smart
   !! used to keep the energy of the previous iteration
   REAL(DP), ALLOCATABLE :: tau_old(:,:)
   !! the atomic positions at the previous iteration
   REAL(DP), ALLOCATABLE :: tau_new(:,:)
   !! the atomic positions at the new iteration
   REAL(DP), ALLOCATABLE :: tau_ref(:,:)
   !! reference atomic positions
   REAL(DP), ALLOCATABLE :: vel(:,:)
   !! velocities
   REAL(DP), ALLOCATABLE :: acc(:,:)
   !! accelerations
   REAL(DP), ALLOCATABLE :: chi(:,:)
   !! chi
   REAL(DP), ALLOCATABLE :: mass(:)
   !! atomic masses
   REAL(DP), ALLOCATABLE :: diff_coeff(:)
   !! diffusion coefficients
   REAL(DP), ALLOCATABLE :: radial_distr(:,:)
   !! radial distribution
   !
   REAL(DP)  :: elapsed_time
   !! elapsed time in ps (picoseconds)
   REAL(DP), PARAMETER, PUBLIC  :: RyDt_to_HaDt = RYDBERG_SI/HARTREE_SI, HaddT_to_RyddT = HARTREE_SI/RYDBERG_SI,& 
                           Ha_to_Ry = HARTREE_SI/RYDBERG_SI
   !! 1/2  and 2 factors used to convert dt from Ry to Ha,  and  velocities from Ha to Ry 
   INTEGER, PARAMETER :: hist_len = 1000
   !
   ! Restart type
   INTEGER, PARAMETER :: restart_verlet = 1, restart_proj_verlet = 2, &
                         restart_fire = 3, restart_langevin = 4,      & 
                         restart_velverlet = 5
   !
CONTAINS
   !
   ! ... public methods
   !
   !------------------------------------------------------------------------
   SUBROUTINE allocate_dyn_vars()
      !------------------------------------------------------------------------
      !! Allocates dynamics variables
      !
      USE ions_base, ONLY : nat
      !
      IF ( .NOT.ALLOCATED( mass ) ) ALLOCATE( mass( nat ) )
      !
      IF ( .NOT.ALLOCATED( tau_old ) ) ALLOCATE( tau_old( 3, nat ) )
      IF ( .NOT.ALLOCATED( tau_new ) ) ALLOCATE( tau_new( 3, nat ) )
      IF ( .NOT.ALLOCATED( tau_ref ) ) ALLOCATE( tau_ref( 3, nat ) )
      !
      IF ( .NOT.ALLOCATED( vel ) ) ALLOCATE( vel( 3, nat ) )
      IF ( .NOT.ALLOCATED( acc ) ) ALLOCATE( acc( 3, nat ) )
      IF ( .NOT.ALLOCATED( chi ) ) ALLOCATE( chi( 3, nat ) )
      !
      IF ( .NOT.ALLOCATED( diff_coeff ) ) ALLOCATE( diff_coeff( nat ) )
      !
      IF ( .NOT.ALLOCATED( radial_distr ) ) &
         ALLOCATE( radial_distr( hist_len , nat ) )
      !
   END SUBROUTINE allocate_dyn_vars
   !
   !
   !------------------------------------------------------------------------
   SUBROUTINE deallocate_dyn_vars()
      !------------------------------------------------------------------------
      !! Deallocates dynamics variables.
      !
      IF ( ALLOCATED( mass ) )          DEALLOCATE( mass )
      IF ( ALLOCATED( tau_old ) )       DEALLOCATE( tau_old )
      IF ( ALLOCATED( tau_new ) )       DEALLOCATE( tau_new )
      IF ( ALLOCATED( tau_ref ) )       DEALLOCATE( tau_ref )
      IF ( ALLOCATED( vel )  )          DEALLOCATE( vel )
      IF ( ALLOCATED( acc )  )          DEALLOCATE( acc )
      IF ( ALLOCATED( chi )  )          DEALLOCATE( chi )
      IF ( ALLOCATED( diff_coeff ) )    DEALLOCATE( diff_coeff )
      IF ( ALLOCATED( radial_distr ) )  DEALLOCATE( radial_distr )
      !
   END SUBROUTINE deallocate_dyn_vars
    !
    !
    !------------------------------------------------------------------------
   SUBROUTINE verlet()
      !------------------------------------------------------------------------
      !! This routine performs one step of molecular dynamics evolution
      !! using the Verlet algorithm. Parameters:
      !
      !! * mass: mass of the atoms;
      !! * dt: time step;
      !! * temperature: starting temperature.
      !
      !! The starting velocities of atoms are set accordingly to the
      !! starting temperature, in random directions. 
      !! The initial velocity distribution is therefore a constant.
      !
      !! Must be run on a single processor, results broadcast to all other procs
      !! Unpredictable behavior may otherwise result due to I/O operations
      !!
      !! Original code: Dario Alfe' 1997  and  Carlo Sbraccia 2004-2006.
      !
      USE ions_base,          ONLY : nat, nsp, ityp, tau, if_pos, atm
      USE ions_nose,          ONLY : vnhp, atm2nhp, ions_nose_energy 
      USE cell_base,          ONLY : alat, omega
      USE ener,               ONLY : etot
      USE force_mod,          ONLY : force
      USE control_flags,      ONLY : istep, lconstrain, tv0rd, tstress, tnosep
      ! istep counts all MD steps, including those of previous runs
      USE constraints_module, ONLY : check_constraint, remove_constr_force, &
         remove_constr_vec, check_wall_constraint
      USE input_parameters,   ONLY : nextffield
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      REAL(DP) :: ekin, etotold
      REAL(DP) :: total_mass, temp_new, temp_av
      REAL(DP) :: delta(3), ml(3), mlt
      INTEGER  :: na
! FIXME: is it useful to keep trace of this possibility?
#undef __NPT
#if defined (__NPT)
#define RELAXTIME 2000.D0
#define TARGPRESS 2.39D0
      REAL(DP) :: chi, press_new
#endif
      LOGICAL  :: is_restart
      REAL(DP), EXTERNAL :: dnrm2
      REAL(DP) :: kstress(3,3), tau_tmp(3, nat)
      INTEGER :: i, j, restart_id
      !
      ndof = get_ndof()
      !
      vel_defined  = .TRUE.
      temp_av      = 0.D0
      !
      CALL seqopn( 4, 'md', 'FORMATTED', is_restart )
      !
      IF ( is_restart ) THEN
         !
         ! ... the file is read :  simulation is continuing
         !
         READ( UNIT = 4, FMT = * ) restart_id
         !
         IF ( restart_id .EQ. restart_verlet ) THEN
            ! Restarting...
            vel_defined = .FALSE.
            !
            ! tau_tmp is read here but not used. It is used for restart in
            ! verlet_read_tau_from_conf
            !
            READ( UNIT = 4, FMT = * ) istep, etotold, tau_tmp(:,:), &
               tau_old(:,:), temp_new, temp_av, mass(:), total_mass, &
               elapsed_time, tau_ref(:,:)
            !
            CLOSE( UNIT = 4, STATUS = 'KEEP' )
            !
         ELSE
            is_restart = .FALSE.
         ENDIF
      ENDIF
      !
      IF (.NOT.is_restart) THEN
         !
         CLOSE( UNIT = 4, STATUS = 'DELETE' )
         !
         ! ... the file is absent :  simulation is starting from scratch
         !
         CALL md_init(total_mass, ekin, temp_new, temp_av, tv0rd)
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
      IF ( control_temp ) CALL apply_thermostat(temp_new, temp_av)
      !
      ! ... Update forces if potential wall constraint was requested
      !
      IF (lconstrain) &
         CALL check_wall_constraint( nat, tau, if_pos, ityp, alat, force )
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
      IF (tnosep) THEN 
         DO na = 1, nat
           acc(:, na) = acc(:,na) - HaddT_to_RyddT * vnhp(atm2nhp(na)) * vel(:,na)
         END DO
      END IF

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
      IF ( .NOT. ANY( if_pos(:,:) == 0 ) .AND. nextffield == 0 ) THEN
         !
         ! ... if no atom has been fixed  we compute the displacement of the
         ! ... center of mass and we subtract it from the displaced positions
         !
         ! ... bypassed if external ionic force fields are activated
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
      ! ... kinetic energy and new temperature are computed here
      !
      vel = ( tau_new - tau_old ) / ( 2.D0*dt ) * DBLE( if_pos )
      !
      CALL compute_ekin( ekin, temp_new )
      !
      ! ... update average temperature
      !
      temp_av = temp_av + temp_new
      !
      ! ... linear momentum and kinetic stress
      !
      ml   = 0.D0
      kstress = 0.d0
      !
      DO na = 1, nat
         !
         ml(:) = ml(:) + vel(:,na) * mass(na)
         DO i = 1, 3
             DO j = 1, 3
                 kstress(i,j) = kstress(i,j) + mass(na)*vel(i,na)*vel(j,na)
             ENDDO
         ENDDO
         !
      ENDDO
      !
      kstress = kstress * alat**2 / omega
      !
      ! FIXME: still useful?
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
      CALL seqopn( 4, 'md', 'FORMATTED',  is_restart )
      !
      WRITE( UNIT = 4, FMT = * ) restart_verlet
      WRITE( UNIT = 4, FMT = * ) istep, etot, tau_new(:,:), tau(:,:), &
         temp_new, temp_av, mass(:), total_mass, elapsed_time, tau_ref(:,:)
      !
      CLOSE( UNIT = 4, STATUS = 'KEEP' )
      !
      CALL dump_trajectory_frame( elapsed_time, temperature )
      !
      ! ... here the tau are shifted
      !
      tau(:,:) = tau_new(:,:)
      !
#if ! defined (__REDUCE_OUTPUT)
      !
      CALL output_tau( .FALSE., .FALSE. )
      !
#endif
      !
      ! ... infos are written on the standard output
      !
      WRITE( stdout, '(5X,"kinetic energy (Ekin) = ",F20.8," Ry",/,  &
                     & 5X,"temperature           = ",F20.8," K ",/,  &
                     & 5X,"Ekin + Etot (const)   = ",F20.8," Ry")' ) &
             ekin, temp_new, ( ekin  + etot )
      IF (tnosep) THEN  
        WRITE (stdout, '(5X,"Ions Nose Energy   = ",    F20.8," Ry",/,  & 
                     &   5X,"Ekin + Etot + Ions Nose =",F20.8," Ry")'), &
                     & Ha_to_Ry *ions_nose_energy, (ekin + etot + Ha_to_Ry * ions_nose_energy)
      END IF  
      !
      IF (tstress) WRITE ( stdout, &
      '(5X,"Ions kinetic stress = ",F15.2," (kbar)",/3(27X,3F15.2/)/)') &
              ((kstress(1,1)+kstress(2,2)+kstress(3,3))/3.d0*ry_kbar), &
              (kstress(i,1)*ry_kbar,kstress(i,2)*ry_kbar,kstress(i,3)*ry_kbar, i=1,3)
      !
      IF ( .NOT.( lconstrain .or. ANY( if_pos(:,:) == 0 ) ) ) THEN
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
      !      !
   END SUBROUTINE verlet
   !
    !------------------------------------------------------------------------
   SUBROUTINE velocity_verlet()
      !------------------------------------------------------------------------
      !! This routine performs two half  steps of molecular dynamics evolution
      !! using the velocity Verlet algorithm. Parameters:
      !
      !! * mass: mass of the atoms;
      !! * dt: time step;
      !! * temperature: starting temperature.
      !! 
      !! reads from previous step tau(t) and vel(t-dt/2)
      !! computes:
      !!   (1)     v(t) = v(t-dt/2) + forces(t)/mass * dt/2 
      !!   (2)     tau(t+dt) = tau(t) + v(t) * dt + forces(t)/mass * dt^2/2
      !!   (3)     v(t+dt/2) = v(t) + forces(t)/mass * dt/2 
      !! saves for next step v(t+dt/2), tau(t + dt) 
      !! :
      !! The starting velocities of atoms are set accordingly to the
      !! starting temperature, in random directions. 
      !! The initial velocity distribution is therefore a constant.
      !
      !! Must be run on a single processor, results broadcast to all other procs
      !! Unpredictable behavior may otherwise result due to I/O operations
      !! 
      !!
      !! Author: Pietro Delugas 2025; based on the Verlet subroutine above: Dario Alfe' 1997  and  Carlo Sbraccia 2004-2006.
      !
      USE ions_base,          ONLY : nat, nsp, ityp, tau, if_pos, atm
      USE ions_nose,          ONLY : vnhp, atm2nhp, ions_nose_energy, ions_nosevel, ions_noseupd, &
                                     ions_nose_shiftvar, xnhp0, xnhpp, xnhpm, qnp, gkbt2nhp, kbt, & 
                                     nhpcl, nhpdim, ekin2nhp, gkbt2nhp, nhpbeg, nhpend, ions_nose_nrg    
      USE cell_base,          ONLY : alat, omega
      USE ener,               ONLY : etot
      USE force_mod,          ONLY : force
      USE control_flags,      ONLY : istep, lconstrain, tv0rd, tstress, tnosep
      ! istep counts all MD steps, including those of previous runs
      USE constraints_module, ONLY : check_constraint, remove_constr_force, &
         remove_constr_vec, check_wall_constraint
      USE input_parameters,   ONLY : nextffield
      !
      IMPLICIT NONE
      !
      ! ... local variables
      !
      REAL(DP) :: ekin, etotold
      REAL(DP) :: total_mass, temp_new, temp_av
      REAL(DP) :: delta(3), ml(3), mlt
      INTEGER  :: na
      LOGICAL  :: is_restart
      REAL(DP), EXTERNAL :: dnrm2
      REAL(DP) :: kstress(3,3), tau_tmp(3, nat), vel_tmp(3, nat)
      INTEGER :: i, j, restart_id
      REAL(DP) :: nose_energy
      !
      ndof = get_ndof()
      !
      vel_defined  = .TRUE.
      temp_av      = 0.D0
      !
      !       
      CALL seqopn( 4, 'md', 'FORMATTED', is_restart )
      !
      IF ( is_restart ) THEN
         !
         ! ... the file is read :  simulation is continuing
         !
         READ( UNIT = 4, FMT = * ) restart_id
         !
         IF ( restart_id .EQ. restart_velverlet ) THEN
            ! Restarting...
                        !
            ! tau_tmp is read here but not used. It is used for restart in
            ! verlet_read_tau_from_conf
            !
            READ( UNIT = 4, FMT = * ) istep, etotold, tau_tmp(:,:), &
               vel(:,:), temp_new, temp_av, mass(:), total_mass, &
               elapsed_time, tau_ref(:,:)
            !
            CLOSE( UNIT = 4, STATUS = 'KEEP' )
            ! for velocity verlet velocity is always defined 
            ! here is vel(t-dt/2) which is the same as (tau(t)-tau(t-dt)/dt in standard verlet
            vel_defined = .TRUE.
            !
         ELSE
            is_restart = .FALSE.
         ENDIF
      ENDIF
      !
      IF (.NOT.is_restart) THEN
         !
         CLOSE( UNIT = 4, STATUS = 'DELETE' )
         !
         ! ... the file is absent :  simulation is starting from scratch
         !
         CALL md_init(total_mass, ekin, temp_new, temp_av, tv0rd)
         !
      ENDIF
      ! this is the nose friction at time t computed with kinetic energy at time (t-dt/2) 
      IF (tnosep) THEN 
         call ions_nosevel(vnhp, xnhp0, xnhpm, RyDt_to_HaDt * 0.5 * dt, nhpcl, nhpdim)
         nose_energy = Ha_to_Ry * ions_nose_nrg(xnhp0, vnhp, qnp, gkbt2nhp, kbt, nhpcl, nhpdim) 
      END IF
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
      ! thermostat applied to v(t-dt/2) here    
      IF ( control_temp ) CALL apply_thermostat(temp_new, temp_av)
      !
      ! ... Update forces if potential wall constraint was requested
      !
      IF (lconstrain) &
         CALL check_wall_constraint( nat, tau, if_pos, ityp, alat, force )
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
      ! ... Velocity Verlet integration scheme
      !
      ! we first compute v(t) 
      !
      FORALL( na = 1:nat) vel(:,na) = vel(:,na) + 0.5_dp * dt * acc(:,na)
      IF (tnosep) THEN 
         FORALL ( na = 1:nat) vel(:,na) = vel(:,na)/(1 + 0.5 * dt * HaddT_to_RyddT * vnhp(atm2nhp(na)))
      END IF 
      !
      ! ... remove the component of the velocity along the
      ! ... constraint gradient
      !
      IF ( lconstrain ) CALL remove_constr_vec( nat, tau, if_pos, ityp, alat, vel )
      !
      IF (tnosep) THEN 
         DO na = 1, nat
           acc(:, na) = acc(:,na) - HaddT_to_RyddT * vnhp(atm2nhp(na)) * vel(:,na)
         END DO
      END IF
      !
      ! then we compute tau(t+dt)
      !
      tau_new(:,:) = tau(:,:) + vel(:,:) * 0.5 * dt + acc(:,:) * dt**2
      !   

      !
      IF ( .NOT. ANY( if_pos(:,:) == 0 ) .AND. nextffield == 0 ) THEN
         !
         ! ... if no atom has been fixed  we compute the displacement of the
         ! ... center of mass and we subtract it from the displaced positions
         !
         ! ... bypassed if external ionic force fields are activated
         ! 
         delta(:) = 0.D0
         DO na = 1, nat
            delta(:) = delta(:) + mass(na)*( tau_new(:,na) - tau(:,na) )
         ENDDO
         delta(:) = delta(:) / total_mass
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
      ! ... we compute v(t+dt/2)
      !
      FORALL (na=1:nat) vel_tmp(:,na) = vel(:,na) + 0.5_dp * dt * acc(:,na) * dble(if_pos(:,na))
      !
      ! ... kinetic energy and new temperature are computed here  with v(t)
      !
      CALL compute_ekin( ekin, temp_new )
      !   ... and we perform the first half-step for the Nose-Hoover chain 
      IF (tnosep) THEN 
         CALL ions_noseupd(xnhpp, xnhp0, xnhpm, 0.5_dp * RyDt_to_HaDt * dt, qnp, ekin2nhp, gkbt2nhp, &
                           vnhp, kbt, nhpcl, nhpdim, nhpbeg, nhpend)
         CALL ions_nose_shiftvar(xnhpp, xnhp0, xnhpm)
      END IF
      ! 
      ! ... update average temperature
      !
      temp_av = temp_av + temp_new
      !
      ! ... linear momentum and kinetic stress
      !
      ml   = 0.D0
      kstress = 0.d0
      !
      DO na = 1, nat
         !
         ml(:) = ml(:) + vel(:,na) * mass(na)
         DO i = 1, 3
             DO j = 1, 3
                 kstress(i,j) = kstress(i,j) + mass(na)*vel(i,na)*vel(j,na)
             ENDDO
         ENDDO
         !
      ENDDO
      !
      kstress = kstress * alat**2 / omega
      !
      !
      ! ... save all the needed quantities on file
      !
      CALL seqopn( 4, 'md', 'FORMATTED',  is_restart )
      !
      WRITE( UNIT = 4, FMT = * ) restart_velverlet
      WRITE( UNIT = 4, FMT = * ) istep, etot, tau_new(:,:), vel_tmp(:,:), &
         temp_new, temp_av, mass(:), total_mass, elapsed_time, tau_ref(:,:)
      !
      CLOSE( UNIT = 4, STATUS = 'KEEP' )
      !
      CALL dump_trajectory_frame( elapsed_time, temperature )
      !
      ! ... here the tau are shifted
      !
      tau(:,:) = tau_new(:,:)
      vel(:,:) = vel_tmp(:,:)

      !
#if ! defined (__REDUCE_OUTPUT)
      !
      CALL output_tau( .FALSE., .FALSE. )
      !
#endif
      !
      ! ... infos are written on the standard output
      !
      WRITE( stdout, '(5X,"kinetic energy (Ekin) = ",F20.8," Ry",/,  &
                     & 5X,"temperature           = ",F20.8," K ",/,  &
                     & 5X,"Ekin + Etot (const)   = ",F20.8," Ry")' ) &
             ekin, temp_new, ( ekin  + etot )
      IF (tnosep) THEN  
        WRITE (stdout, '(5X,"Ions Nose Energy   = ",    F20.8," Ry",/,  & 
                     &   5X,"Ekin + Etot + Ions Nose =",F20.8," Ry")'), &
                     & Ha_to_Ry * nose_energy, (ekin + etot + Ha_to_Ry * nose_energy)
      END IF  
      !
      IF (tstress) WRITE ( stdout, &
      '(5X,"Ions kinetic stress = ",F15.2," (kbar)",/3(27X,3F15.2/)/)') &
              ((kstress(1,1)+kstress(2,2)+kstress(3,3))/3.d0*ry_kbar), &
              (kstress(i,1)*ry_kbar,kstress(i,2)*ry_kbar,kstress(i,3)*ry_kbar, i=1,3)
      !
      IF ( .NOT.( lconstrain .or. ANY( if_pos(:,:) == 0 ) ) ) THEN
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
      ! after printout of quantities if needed we perform the second half-step of Nose-Hoover chains with vel = v(t+dt/2)
      IF (tnosep) THEN 
         CALL compute_ekin(ekin, temp_new)
         CALL ions_noseupd(xnhpp, xnhp0, xnhpm, 0.5_dp * RyDt_to_HaDt * dt, qnp, ekin2nhp, gkbt2nhp,&
                           vnhp, kbt, nhpcl, nhpdim, nhpbeg, nhpend) 
         CALL ions_nose_shiftvar(xnhpp, xnhp0, xnhpm) 
      END IF 
      !
      !      !
   END SUBROUTINE velocity_verlet

   !--------------------------------------------------------------------
   SUBROUTINE md_init(total_mass, ekin, temp_new, temp_av, tv0rd)
      !--------------------------------------------------------------------
      !
      USE ions_base, ONLY: tau, nsp, ityp, nat, atm, amass
      USE control_flags, ONLY: istep
      USE cell_base, ONLY: alat
      USE ions_nose, ONLY: ions_nose_info
      IMPLICIT NONE
      REAL(DP), INTENT(OUT) :: total_mass, ekin, temp_new, temp_av
      LOGICAL, INTENT(IN)   :: tv0rd
      INTEGER :: na
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
         SELECT CASE( TRIM( thermostat ) )
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
         CASE( 'svr', 'Svr', 'SVR' )
            !
            WRITE( UNIT = stdout, &
                  FMT = '(/,5X,"temperature is controlled by ", &
                         &     "stochastic velocity rescaling",/,5x,&
                         &     "Characteristic time   =",i3,"*timestep")') &
                            nraise
         CASE( 'initial', 'Initial' )
            !
            WRITE( UNIT = stdout, &
                  FMT = '(/,5X,"temperature is set once at start"/)' )
            !
         CASE ('nose') 
            WRITE( UNIT = stdout, &
                  FMT = '(/,5X,"temperature is controlled by ", &
                         &     "Nose Hoover thermostat",/,5x)') 
            CALL ions_nose_info(dt)
         CASE DEFAULT
            !
            WRITE( UNIT = stdout, &
                  FMT = '(/,5X,"temperature is controlled by ", &
                           &     "velocity rescaling (",A,")"/)' )&
                           TRIM( thermostat )
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
      IF ( tv0rd ) THEN
         !
         ! ... initial velocities available from input file
         !
         vel(:,:) = vel(:,:) / alat
         vel_defined = .true.
         !
         CALL compute_ekin ( ekin, temp_new )
         !
      ELSEIF ( control_temp ) THEN
         !
         ! ... initial thermalization. N.B. tau is in units of alat
         !
         CALL start_therm()
         vel_defined = .TRUE.
         !
         temp_new = temperature
         !
         temp_av = 0.D0
         !
      ELSE
         !
         vel(:,:) = 0.0_DP
         vel_defined = .TRUE.
         !
      ENDIF
      !
      elapsed_time = 0.D0
      !
   END SUBROUTINE md_init

   !-----------------------------------------------------------------------
   SUBROUTINE start_therm()
      !-----------------------------------------------------------------------
      !! Starting thermalization of the system.
      !
      USE control_flags,  ONLY: lconstrain
      USE symm_base,      ONLY : invsym, nsym, irt
      USE cell_base,      ONLY : alat
      USE ions_base,      ONLY : nat, if_pos, ityp, tau
      USE constraints_module, ONLY: remove_constr_vec
      USE random_numbers, ONLY : gauss_dist
      !
      IMPLICIT NONE
      !
      INTEGER  :: na_, nb
      REAL(DP) :: total_mass_, kt, sigma, ek, ml_(3), system_temp
      !
      kt = temperature / ry_to_kelvin
      !
      ! ... starting velocities have a Maxwell-Boltzmann distribution
      !
      DO na_ = 1, nat
         !
         sigma = SQRT( kt / mass(na_) )
         !
         ! ... N.B. velocities must in a.u. units of alat
         !
         vel(:,na_) = gauss_dist( 0.D0, sigma, 3 ) / alat
         !
      ENDDO
      !
      ! ... the velocity of fixed ions must be zero
      !
      vel = vel * DBLE( if_pos )
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
         DO na_ = 1, nat
            !
            nb = irt( ( nsym / 2 + 1 ), na_ )
            !
            IF ( nb > na_ ) vel(:,nb) = - vel(:,na_)
            !
            ! ... the atom on the inversion center is kept fixed
            !
            IF ( na_ == nb ) vel(:,na_) = 0.D0
            !
         ENDDO
         !
      ELSE
         !
         ! ... put total linear momentum equal zero if all atoms
         ! ... are free to move
         !
         ml_(:) = 0.D0
         !
         IF ( .NOT. ANY( if_pos(:,:) == 0 ) ) THEN
            !
            total_mass_ = SUM( mass(1:nat) )
            DO na_ = 1, nat
               ml_(:) = ml_(:) + mass(na_)*vel(:,na_)
            ENDDO
            ml_(:) = ml_(:) / total_mass_
            !
         ENDIF
         !
      ENDIF
      !
      ek = 0.D0
      !
      DO na_ = 1, nat
         !
         vel(:,na_) = vel(:,na_) - ml_(:)
         !
         ek = ek + 0.5D0 * mass(na_) * &
                  ( ( vel(1,na_) )**2 + ( vel(2,na_) )**2 + ( vel(3,na_) )**2 )
         !
      ENDDO
      !
      ! ... after the velocity of the center of mass has been subtracted the
      ! ... temperature is usually changed. Set again the temperature to the
      ! ... right value.
      !
      system_temp = 2.D0 / DBLE( ndof ) * ek * alat**2 * ry_to_kelvin
      !
      CALL thermalize( 0, system_temp, temperature )
      !
   END SUBROUTINE start_therm
      !--------------------------------------------------------------------
   SUBROUTINE apply_thermostat(temp_new, temp_av)
      !--------------------------------------------------------------------
      !
      USE random_numbers,    ONLY : randy, gauss_dist
      USE ions_base,         ONLY:  tau, if_pos, nat
      USE control_flags,     ONLY:  istep
      USE cell_base,         ONLY:  alat
      !
      IMPLICIT NONE
      REAL(DP), INTENT(IN)  :: temp_new
      REAL(DP), INTENT(OUT) :: temp_av
      INTEGER :: nat_moved
      REAL(DP) :: sigma, kt
      INTEGER  :: na
      !
      IF (.NOT. vel_defined) THEN
         vel(:,:) = ( tau(:,:) - tau_old(:,:) ) / dt * dble( if_pos(:,:) )
      ENDIF
      !
      SELECT CASE( TRIM( thermostat ) )
      CASE( 'rescaling' )
         !
         IF ( ABS(temp_new-temperature) > tolp ) THEN
            WRITE( UNIT = stdout, &
                  FMT = '(/,5X,"Velocity rescaling: T (",F6.1,"K) ", &
                              & "out of range, reset to " ,F6.1)' ) &
                        temp_new, temperature
            CALL thermalize( 0, temp_new, temperature )
         ENDIF
         !
      CASE( 'rescale-v', 'rescale-V', 'rescale_v', 'rescale_V' )
         !
         IF ( MOD( istep, nraise ) == 0 ) THEN
            !
            temp_av = temp_av / DBLE( nraise )
            !
            WRITE( UNIT = stdout, &
                  FMT = '(/,5X,"Velocity rescaling: average T on ",i3, &
                              &" steps (",F6.1,"K) reset to ",F6.1)' )  &
                        nraise, temp_av, temperature
            !
            CALL thermalize( 0, temp_new, temperature )
            temp_av = 0.D0
         ENDIF
         !
      CASE( 'rescale-T', 'rescale-t', 'rescale_T', 'rescale_t' )
         ! Clearly it makes sense to check for positive delta_t
         ! If a negative delta_t is given, I suggest to have a message
         ! printed, that delta_t is ignored (TODO)
         IF ( delta_t > 0 ) THEN
            !
            temperature = temperature*delta_t
            !
            WRITE( UNIT = stdout, &
                  FMT = '(/,5X,"Thermalization: T (",F6.1,"K) rescaled ",&
                              & "by a factor ",F6.3)' ) temp_new, delta_t
            !
            CALL thermalize( 0, temp_new, temperature )
            !
         ENDIF
      CASE( 'reduce-T', 'reduce-t', 'reduce_T', 'reduce_t' )
         IF ( mod( istep, nraise ) == 0 ) THEN
            !
            ! First printing message, than reduce target temperature:
            !
            IF ( delta_t > 0 ) THEN
              WRITE( UNIT = stdout, &
                  FMT = '(/,5X,"Thermalization: T (",F6.1,"K) augmented ",&
                              & "by ",F6.3)' ) temperature, delta_t

            ELSE
              WRITE( UNIT = stdout, &
                  FMT = '(/,5X,"Thermalization: T (",F6.1,"K) reduced ",&
                              & "by ",F6.3)' ) temperature, -delta_t
            ENDIF
            ! I check whether the temperature is negative, so that I avoid
            ! nonsensical behavior:
            IF (temperature < 0.0D0 ) CALL errore('apply_thermostat','Negative target temperature',1)
            !
            temperature = temperature + delta_t
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
      CASE( 'svr', 'Svr', 'SVR' )
         !
         WRITE( UNIT = stdout, &
             FMT = '(/,5X,"Canonical sampling velocity rescaling")' )
         !
         CALL thermalize_resamp_vscaling( nraise, temp_new, temperature )
         !
      CASE( 'andersen', 'Andersen' )
         !
         kt = temperature / ry_to_kelvin
         nat_moved = 0
         !
         DO na = 1, nat
            !
            IF ( randy() < 1.D0 / DBLE( nraise ) ) THEN
               !
               nat_moved = nat_moved + 1
               sigma = SQRT( kt / mass(na) )
               !
               ! ... N.B. velocities must in a.u. units of alat and are zero
               ! ...      for fixed ions
               !
               vel(:,na) = DBLE( if_pos(:,na) ) * &
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
      IF (.NOT. vel_defined) THEN
         tau_old(:,:) = tau(:,:) - vel(:,:) * dt
      ENDIF
      !
   END SUBROUTINE apply_thermostat
      !


   !------------------------------------------------------------------------
   SUBROUTINE compute_ekin ( ekin, temp_new )
     !------------------------------------------------------------------------
     USE cell_base,      ONLY : alat
     USE ions_base,      ONLY : nat
     USE ions_nose,      ONLY : ekin2nhp, atm2nhp
     USE control_flags,  ONLY : tnosep
     IMPLICIT NONE  
     REAL (dp), INTENT (out) :: ekin, temp_new
     INTEGER :: na
     REAL(DP), parameter :: Ry_to_Ha = 1._dp/Ha_to_Ry 
     REAL(dp) :: ekin_at 
     !
     ekin = 0.0_dp
     ekin_at = 0.0_dp 
     if (tnosep) ekin2nhp = 0.0_dp 
     DO na = 1, nat
        ekin_at  =  0.5_dp * mass(na) * &
             ( vel(1,na)**2 + vel(2,na)**2 + vel(3,na)**2 )
        ekin = ekin + ekin_at
        IF  (tnosep)  ekin2nhp(atm2nhp(na)) = ekin2nhp(atm2nhp(na)) + Ry_to_Ha * ekin_at*alat**2 
     END DO
     ekin = ekin*alat**2
     temp_new = 2.D0 / DBLE( ndof ) * ekin * ry_to_kelvin
     !
   END SUBROUTINE compute_ekin
   !
   !------------------------------------------------------------------------
   SUBROUTINE terminate_verlet
     !------------------------------------------------------------------------
     !! Terminate Verlet molecular dynamics calculation.
     !
     USE io_global, ONLY : stdout
     USE control_flags, ONLY : istep
     !
     WRITE( UNIT = stdout, &
          FMT = '(/,5X,"The maximum number of steps has been reached.")' )
     WRITE( UNIT = stdout, &
          FMT = '(/,5X,"End of molecular dynamics calculation")' )
     !
     IF (istep > 0) CALL print_averages()
     !
   END SUBROUTINE terminate_verlet
   !
   !
   !------------------------------------------------------------------------
   SUBROUTINE proj_verlet( conv_ions )
      !------------------------------------------------------------------------
      !! This routine performs one step of structural relaxation using
      !! the preconditioned-projected-Verlet algorithm. 
      !
      USE ions_base,          ONLY : nat, ityp, tau, if_pos
      USE cell_base,          ONLY : alat
      USE ener,               ONLY : etot
      USE force_mod,          ONLY : force
      USE relax,              ONLY : epse, epsf
      USE control_flags,      ONLY : istep, lconstrain
      !
      USE constraints_module, ONLY : remove_constr_force, check_constraint
      USE input_parameters,   ONLY : nextffield
      ! TB
      USE extfield,           ONLY : relaxz
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(OUT) :: conv_ions
      !
      REAL(DP), ALLOCATABLE :: step(:,:)
      REAL(DP)              :: norm_step, etotold, delta(3)
      INTEGER               :: na, restart_id
      LOGICAL               :: is_restart
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
      CALL seqopn( 4, 'md', 'FORMATTED', is_restart )
      !
      IF ( is_restart ) THEN
         !
         ! ... the file is read
         !
         READ( UNIT = 4, FMT = * ) restart_id
         !
         IF ( restart_id .EQ. restart_proj_verlet ) THEN
            !
            READ( UNIT = 4, FMT = * ) istep, etotold, tau_old(:,:)
            CLOSE( UNIT = 4, STATUS = 'KEEP' )
            !
         ELSE
            is_restart = .FALSE.
         END IF
      END iF
      !
      IF (.NOT.is_restart) THEN
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
      conv_ions = conv_ions .and. ( MAXVAL( ABS( force ) ) < epsf )
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
         CALL output_tau( .TRUE., .TRUE. )
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
      tau_new(:,:) = tau(:,:) + step(:,:)*MIN( norm_step, step_max / alat )
      !
      ! TB
      !IF ( .NOT. ANY( if_pos(:,:) == 0 ) ) THEN
      IF ( .NOT. ANY( if_pos(:,:) == 0 ) .AND. (relaxz) ) THEN
         WRITE( stdout, '("relaxz = .TRUE. => displacement of the center of mass is not subtracted")')
      ENDIF
      IF ( (.NOT. ANY( if_pos(:,:) == 0 )) .AND. (.NOT. relaxz) .AND. nextffield==0 ) THEN
         !
         ! ... if no atom has been fixed  we compute the displacement of the
         ! ... center of mass and we subtract it from the displaced positions
         !
         ! ... also bypassed if external ionic force fields are activated
         ! 
         delta(:) = 0.D0
         !
         DO na = 1, nat
            !
            delta(:) = delta(:) + ( tau_new(:,na) - tau(:,na) )
            !
         ENDDO
         !
         delta(:) = delta(:) / DBLE( nat )
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
      CALL seqopn( 4, 'md', 'FORMATTED',  is_restart )
      !
      WRITE( UNIT = 4, FMT = * ) restart_proj_verlet
      WRITE( UNIT = 4, FMT = * ) istep, etot, tau(:,:)
      !
      CLOSE( UNIT = 4, STATUS = 'KEEP' )
      !
      ! ... here the tau are shifted
      !
      tau(:,:) = tau_new(:,:)
      !
#if ! defined (__REDUCE_OUTPUT)
      !
      CALL output_tau( .FALSE., .FALSE. )
      !
#endif
      !
      DEALLOCATE( step )
      !
   END SUBROUTINE proj_verlet

   SUBROUTINE fire( conv_ions )

     !------------------------------------------------------------------------
     !! This routine performs one step of structural relaxation using
     !! the FIRE (Fast Inertial Relaxation Engine) algorithm using the
     !! semi-implicit Euler integration scheme with an energy monitor;
     !! 
     !! References: (1) Bitzek et al., Phys. Rev. Lett.,  97, 170201, (2006),
     !!                  doi: 10.1103/PhysRevLett.97.170201
     !!             (2) Shuang et al., Comp. Mater. Sci., 156, 135-141 (2019),
     !!                  doi: 10.1016/j.commatsci.2018.09.049
     !!             (3) (FIRE 2.0) Guénolé et al., Comp. Mater. Sci., 175, 109584, (2020),
     !!                  doi: 10.1016/j.commatsci.2020.109584
     !! 
     !!  There are seven global parameters that can be modified by the user:
     !!
     !!  dt .... initial time step of calculation
     !!  fire_nmin ... number of steps with P > 0 before dt is increased 
     !!  fire_f_inc ... factor for the increase of the time step
     !!  fire_f_dec ... factor for the decrease of the time step
     !!  fire_alpha_init ... initial value of the velocity mixing factor
     !!  fire_falpha ... modifies the velocity mixing factor
     !!  fire_dtmax ... maximum time step; calculated as dtmax = fire_dtmax*dt 
     !!
     !! Defaults are (taken from ref (2)):
     !!   dt = 20.0  (the default time step in PW ==  20.0 a.u. or 0.9674 fs )
     !!   fire_f_inc = 1.1 
     !!   fire_f_dec = 0.5
     !!   fire_f_alpha_init = 0.2
     !!   fire_falpha = 0.99
     !!----------------------------------------------------------------------- 
     USE ions_base,          ONLY : nat, ityp, tau, if_pos
     USE cell_base,          ONLY : alat
     USE ener,               ONLY : etot
     USE force_mod,          ONLY : force
     USE relax,              ONLY : epse, epsf
     USE control_flags,      only : istep, lconstrain
     !
     USE constraints_module, ONLY : remove_constr_force, remove_constr_vec, check_constraint
     ! TB
     USE extfield,      ONLY : relaxz
     !
     IMPLICIT NONE
     !
     LOGICAL, INTENT(OUT) :: conv_ions
     !
     REAL(DP), ALLOCATABLE :: step(:,:)
     REAL(DP)              :: norm_step, etotold, delta(3)
     INTEGER               :: na, i, restart_id
     LOGICAL               :: is_restart
     !
     REAL(DP), PARAMETER :: step_max = 0.6D0  ! bohr
     ! 
     REAL(DP), EXTERNAL :: dnrm2,ddot
     !
     ! FIRE local variables 
     !
     REAL(DP) :: P, alpha, fmax, dt_max, dt_curr
     INTEGER :: nsteppos
     ! 
     ! FIRE  parameters; read from input ...   
     !
     INTEGER ::  Nmin  ! minimum number of steps for time step increase 
     REAL(DP) :: f_inc ! factor for time step increase  
     REAL(DP) :: f_dec ! factor for time step decrease
     REAL(DP) :: alpha_init ! initial value of mixing factor
     REAL(DP) :: falpha ! modify the mixing factor 
     !
     ALLOCATE( step( 3, nat ) )
     !
     ! set up local variables for global input parameters (better readability) ... 
     !
     Nmin = fire_nmin
     f_inc = fire_f_inc
     f_dec = fire_f_dec
     alpha_init = fire_alpha_init
     falpha = fire_falpha
     !
     ! initialize alpha
     ! 
     tau_new(:,:) = 0.D0
     alpha = alpha_init
     ! maximum time step 
     dt_curr = dt
     dt_max = fire_dtmax*dt
     ! acc_old and vel_old are here to keep track of acceleration/velocity in the previous time step
     nsteppos = 0
     conv_ions = .FALSE.
     !
     CALL seqopn( 4, 'fire', 'FORMATTED', is_restart )
     !
     IF ( is_restart ) THEN
        !
        ! ... the file is read ...
        !
        READ( UNIT = 4, FMT = * ) restart_id
        !
        IF (restart_id .EQ. restart_fire) THEN
           !
           READ( UNIT = 4, FMT = * ) etotold, nsteppos, dt_curr, alpha
           CLOSE( UNIT = 4, STATUS = 'KEEP' )
           !
        ELSE
           !
           is_restart = .FALSE.
           !
        END IF
        !
     END IF
     !
     IF (.NOT.is_restart) THEN
        !
        CLOSE( UNIT = 4, STATUS = 'DELETE' )
        !
        ! ... atoms are refolded in the central box
        !
        IF ( refold_pos ) CALL refold_tau()
        !
        vel(:,:) = 0.D0
        acc(:,:) = 0.D0
        etotold = etot
        istep = 0
        WRITE( UNIT = stdout, &
             FMT = '(/,5X,"Minimization using the FIRE algorithm")' )
        !
        ! write out the input parameters 
        WRITE (UNIT = stdout, &
             FMT = '(/,5X,"FIRE input parameters:")')
        WRITE (UNIT = stdout, &
             FMT = '(/,5X," fire_nmin = ",I2," fire_f_inc = ", F4.2, & 
             " fire_f_dec = ",F4.2," fire_alpha = ",F4.2, & 
             " fire_falpha = ",F4.2, " dtmax = ",F5.1 )') &
             fire_nmin, fire_f_inc, fire_f_dec, &
             fire_alpha_init, fire_falpha, dt_max 
        !
     ENDIF
     !
     IF ( lconstrain ) THEN
        !
        ! ... we first remove the component of the force along the
        ! ... constraint gradient (this constitutes the initial guess
        ! ... for the calculation of the lagrange multipliers)
        !
        write (*,*) "Called remove constr" 
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
     conv_ions = conv_ions .and. ( MAXVAL( ABS( force ) ) < epsf )
     !
     IF ( conv_ions ) THEN
        !
        WRITE( UNIT = stdout, &
             FMT = '(/,5X,"FIRE: convergence achieved in " &
             & ,I3," steps")' ) istep
        WRITE( UNIT = stdout, &
             FMT = '(/,5X,"End of FIRE minimization")' )
        WRITE( UNIT = stdout, &
             FMT = '(/,5X,"Final energy = ",F18.10," Ry"/)' ) etot
        !
        CALL output_tau( .TRUE., .TRUE. )
        !
        RETURN
        !
     ENDIF
     !
     istep = istep + 1
     WRITE( stdout, '(/,5X,"Entering FIRE :",&
          & T28,"iteration",T37," = ",I5)' ) istep

     !
     ! calculate acceleration
     !
     acc(:,:) = force(:,:) / alat / amu_ry
     ! 
     ! calculate the projection of the velocity on the force 
     P = ddot(3*nat,force, 1, vel, 1)
     ! 
     step(:,:) = 0.0_DP
     ! 
     IF ( P < 0.0_DP  )  THEN 
        ! FIRE 2.0 algorithm: if P < 0 go back by half a step 
        ! for details see reference (2), doi: 10.1016/j.commatsci.2020.109584
        step(:,:) = step(:,:) - 0.5_DP*dt_curr*vel(:,:)
     ENDIF
     !
     ! ... manipulate the time step ... 
     !
     ! NOTEs:
     ! in original FIRE the condition is P > 0,
     ! however to prevent the time step decrease in the first step where v=0
     ! (P=0 and etot=etotold) the equality was changed to P >= 0
     ! The energy difference criterion is also added to prevent
     ! the minimization from going uphill  
     ! 

     IF ( P >= 0.0_DP  .AND. (etot - etotold) <= 0.D0 ) THEN
        ! 
        ! 
        nsteppos = nsteppos + 1
        ! increase time step and modify mixing factor only after Nmin steps in positive direction
        IF ( nsteppos > Nmin ) THEN 
           dt_curr = MIN(dt_curr*f_inc, dt_max )
           alpha = alpha*falpha
        END IF
     ELSE   
        !
        ! set velocity to 0; return alpha to the initial value; reduce time step
        !
        vel(:,:) = 0.D0
        alpha = alpha_init
        nsteppos = 0
        dt_curr = dt_curr*f_dec        
     END IF
     ! report current parameters 
     WRITE (stdout, '(/,5X, "FIRE Parameters: P = ", F10.8, ", dt = ", F5.2, & 
          &  ", alpha = ", F5.3, ", nsteppos = ", I3, " at step ", I3, /)' ) &
          & P, dt_curr, alpha, nsteppos, istep
     !
     ! calculate v(t+dt) = v(t) + a(t)*dt  
     !
     vel(:,:) = vel(:,:) + dt_curr*acc(:,:)
     ! 
     ! velocity mixing 
     !
     vel(:,:) = (1.D0 - alpha)*vel(:,:) + alpha*force(:,:)*dnrm2(3*nat,vel,1)/dnrm2(3*nat,force,1)
     !
     !
     ! ... the velocity of fixed ions must be zero
     !
     vel = vel * DBLE( if_pos )
     ! 
     IF ( lconstrain )  THEN
        !
        ! apply constraints to the velocity as well  
        !
        CALL remove_constr_vec( nat, tau, if_pos, ityp, alat, vel )
     ENDIF
     ! 
     ! 
     ! calculate the displacement x(t+dt) = x(t) + v(t+dt)*dt 
     ! 
     step(:,:) = step(:,:) +  vel(:,:)*dt_curr
     !
     norm_step = dnrm2( 3*nat, step, 1 )
     !
     step(:,:) = step(:,:) / norm_step
     !
     ! keep the step within a threshold (taken from damped dynamics)  
     !
     tau_new(:,:) = tau(:,:) + step(:,:)*MIN(norm_step, step_max/alat)
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
     ! ... save the needed quantities to a file 
     !
     CALL seqopn( 4, 'fire', 'FORMATTED',  is_restart )
     !
     WRITE( UNIT = 4, FMT = * ) restart_fire
     WRITE( UNIT = 4, FMT = * ) etot, nsteppos, dt_curr, alpha      !
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     !
     ! ... displace tau, and output new positions   
     !
     tau(:,:) = tau_new(:,:)
     !
     !
#if ! defined (__REDUCE_OUTPUT)
     !
     CALL output_tau( .FALSE., .FALSE. )
     !
#endif
     !
     DEALLOCATE( step )
     !
   END SUBROUTINE fire 
   !
   !
   !------------------------------------------------------------------------
   SUBROUTINE langevin_md()
      !------------------------------------------------------------------------
      !! Langevin molecular dynamics.
      !
      USE ions_base,      ONLY : nat, ityp, tau, if_pos
      USE cell_base,      ONLY : alat
      USE ener,           ONLY : etot
      USE force_mod,      ONLY : force
      USE control_flags,  ONLY : istep, lconstrain
      USE random_numbers, ONLY : gauss_dist
      !
      USE constraints_module, ONLY : remove_constr_force, check_constraint
      USE input_parameters,   ONLY : nextffield
      !
      IMPLICIT NONE
      !
      REAL(DP) :: sigma, kt
      REAL(DP) :: delta(3)
      INTEGER  :: na, restart_id
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
         READ( UNIT = 4, FMT = * ) restart_id
         IF ( restart_id == restart_langevin ) THEN
            READ( UNIT = 4, FMT = * ) istep
            CLOSE( UNIT = 4, STATUS = 'KEEP' )
         ELSE
            file_exists = .FALSE.
            CLOSE( UNIT = 4, STATUS = 'DELETE' )
         END IF
         !
      END IF
      !
      IF ( .NOT.file_exists ) THEN
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
      sigma = SQRT( 2.D0*dt*kt )
      !
      delta(:) = 0.D0
      !
      DO na = 1, nat
         !
         chi(:,na) = gauss_dist( 0.D0, sigma, 3 )*DBLE( if_pos(:,na) )
         !
         delta(:) = delta(:) + chi(:,na)
         !
      ENDDO
      !
      FORALL( na = 1:nat ) chi(:,na) = chi(:,na) - delta(:) / DBLE( nat )
      !
      PRINT *, "|F|   = ", dt*dnrm2( 3*nat, force, 1 )
      PRINT *, "|CHI| = ", dnrm2( 3*nat, chi, 1 )
      !
      ! ... over-damped Langevin dynamics
      !
      tau_new(:,:) = tau(:,:) + ( dt*force(:,:) + chi(:,:) ) / alat
      !
      IF ( .NOT. ANY( if_pos(:,:) == 0 ) .AND. nextffield==0) THEN
         !
         ! ... here we compute the displacement of the center of mass and we
         ! ... subtract it from the displaced positions
         !
         ! ... also bypassed if external ionic force fields are activated
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
      WRITE( UNIT = 4, FMT = * ) restart_langevin
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
      CALL output_tau( .FALSE., .FALSE. )
      !
#endif
      !
   END SUBROUTINE langevin_md
   !
   !
   !-----------------------------------------------------------------------
   SUBROUTINE refold_tau()
      !-----------------------------------------------------------------------
      !! Refold atomic positions.
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
   !
   !-----------------------------------------------------------------------
   SUBROUTINE compute_averages( istep )
      !-----------------------------------------------------------------------
      !! Molecular dynamics - compute averages.
      !
      USE ions_base,          ONLY : nat, tau, fixatom
      USE cell_base,          ONLY : alat, at
      USE constraints_module, ONLY : pbc
      USE io_files,           ONLY : delete_if_present
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: istep
      !! md step
      !
      ! ... local variables
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
         CALL delete_if_present( TRIM( tmp_dir ) // &
                               & TRIM( prefix ) // ".msd.dat" )
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
      diff_coeff(:) = msd(:) / ( 6.D0*DBLE( istep )*dt )
      !
      ! ... conversion from Rydberg atomic units to cm^2/sec
      !
      diff_coeff(:) = diff_coeff(:) * bohr_radius_cm**2 / ( 2.D-12*au_ps )
      !
      OPEN( UNIT = 4, POSITION = 'APPEND', &
            FILE = TRIM( tmp_dir ) // TRIM( prefix ) // ".msd.dat" )
      !
      WRITE( 4, '(2(2X,F16.8))' ) &
          ( istep*dt*2.D0*au_ps ), SUM( msd(:) ) / DBLE( nat-fixatom )
      !
      CLOSE( UNIT = 4, STATUS = 'KEEP' )
      !
      DEALLOCATE( msd )
      !
      ! ... radial distribution function g(r)
      !
      inv_dmax = 1.D0 / ( norm( MATMUL( at(:,:), max_dist(:) ) ) * alat )
      !
      DO i = 1, nat
         !
         DO j = 1, nat
            !
            IF ( i == j ) CYCLE
            !
            dtau(:) = pbc( ( tau(:,i) - tau(:,j) ) * alat )
            !
            idx = ANINT( norm( dtau(:) ) * inv_dmax * DBLE( hist_len ) )
            !
            IF( idx > 0 .and. idx <= SIZE( radial_distr, 1 ) ) &
               radial_distr(idx,i) = radial_distr(idx,i) + 1.D0
            !
         ENDDO
         !
      ENDDO
      !
   END SUBROUTINE compute_averages
   !
   !
   !-----------------------------------------------------------------------
   SUBROUTINE print_averages()
      !-----------------------------------------------------------------------
      !! Molecular dynamics - print averages.
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
          sum( diff_coeff(:) ) / DBLE( nat-fixatom )
      !
      ! ... radial distribution function g(r)
      !
      dmax = norm( matmul( at(:,:), max_dist(:) ) ) * alat
      !
      radial_distr(:,:) = radial_distr(:,:) * omega / DBLE( nat ) / fpi
      !
      radial_distr(:,:) = radial_distr(:,:) / ( dmax / DBLE( hist_len ) )
      !
      radial_distr(:,:) = radial_distr(:,:) / DBLE( nstep )
      !
      OPEN( UNIT = 4, FILE = TRIM( tmp_dir ) // TRIM( prefix ) // ".rdf.dat" )
      !
      DO idx = 1, hist_len
         !
         dist = DBLE( idx ) / DBLE( hist_len ) * dmax
         !
         IF ( dist > dmax / SQRT( 3.0d0 ) ) CYCLE
         !
         radial_distr(idx,:) = radial_distr(idx,:) / dist**2
         !
         WRITE( 4, '(2(2X,F16.8))' ) &
             dist, SUM( radial_distr(idx,:) ) / DBLE( nat )
         !
      ENDDO
      !
      CLOSE( UNIT = 4 )
      !
   END SUBROUTINE print_averages
   !
   !
   !-----------------------------------------------------------------------
   SUBROUTINE force_precond( istep, force, etotold )
      !-----------------------------------------------------------------------
      !! This routine computes an estimate of \(H^{-1}\) by using the BFGS
      !! algorithm and the preconditioned gradient \(\text{pg} = H^{-1} g\).
      !! It works in atomic units.
      !
      USE ener,        ONLY : etot
      USE cell_base,   ONLY : alat
      USE ions_base,   ONLY : nat, tau
      !
      IMPLICIT NONE
      !
      INTEGER,  INTENT(IN)    :: istep
      REAL(DP), INTENT(INOUT) :: force(:,:)
      REAL(DP), INTENT(IN)    :: etotold
      !
      REAL(DP), ALLOCATABLE :: pos(:), pos_p(:)
      REAL(DP), ALLOCATABLE :: grad(:), grad_p(:), precond_grad(:)
      REAL(DP), ALLOCATABLE :: inv_hess(:,:)
      REAL(DP), ALLOCATABLE :: y(:), s(:)
      REAL(DP), ALLOCATABLE :: Hy(:), yH(:)
      REAL(DP)              :: sdoty, pg_norm
      INTEGER               :: dim
      INTEGER               :: iunbfgs
      CHARACTER(LEN=256)    :: bfgs_file
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
      pos(:)  =   RESHAPE( tau,   (/ dim /) ) * alat
      grad(:) = - RESHAPE( force, (/ dim /) )
      !
      bfgs_file = TRIM( tmp_dir ) // TRIM( prefix ) // '.bfgs'
      !
      INQUIRE( FILE = TRIM( bfgs_file ) , EXIST = file_exists )
      !
      IF ( file_exists ) THEN
         !
         OPEN( NEWUNIT = iunbfgs, &
               FILE = TRIM( bfgs_file ), STATUS = 'OLD', ACTION = 'READ' )
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
         IF ( MOD( istep, nrefresh ) == 0 ) inv_hess(:,:) = identity( dim )
         !
         IF ( etot < etotold ) THEN
            !
            ! ... BFGS update
            !
            s(:) = pos(:)  - pos_p(:)
            y(:) = grad(:) - grad_p(:)
            !
            sdoty = ( s(:) .DOT. y(:) )
            !
            IF ( sdoty > eps8 ) THEN
               !
               Hy(:) = ( inv_hess(:,:) .TIMES. y(:) )
               yH(:) = ( y(:) .TIMES. inv_hess(:,:) )
               !
               inv_hess = inv_hess + 1.D0 / sdoty * &
                        ( ( 1.D0 + ( y .DOT. Hy ) / sdoty ) * matrix( s, s ) - &
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
      precond_grad(:) = ( inv_hess(:,:) .TIMES. grad(:) )
      !
      IF ( ( precond_grad(:) .DOT. grad(:) ) < 0.D0 ) THEN
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
      OPEN( NEWUNIT = iunbfgs, &
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
      precond_grad(:) = precond_grad(:) / pg_norm
      precond_grad(:) = precond_grad(:) * MIN( pg_norm, max_pg_norm )
      !
      force(:,:) = - RESHAPE( precond_grad(:), (/ 3, nat /) )
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
      !! Quick-min algorithm.
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
      vel(:,:) = acc_versor(:,:) * MAX( 0.D0, projection )
      !
      DEALLOCATE( acc_versor )
      !
   END SUBROUTINE project_velocity
   !
   !-----------------------------------------------------------------------
   SUBROUTINE thermalize( nraise_, system_temp, required_temp )
      !-----------------------------------------------------------------------
      !! * Berendsen rescaling (Eq. 7.59 of Allen & Tildesley);
      !! * rescale the velocities by a factor 3 / 2KT / Ek.
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: system_temp, required_temp
      INTEGER, INTENT(IN) :: nraise_
      !
      REAL(DP) :: aux
      !
      IF ( nraise_ > 0 ) THEN
         !
         ! ... Berendsen rescaling (Eq. 7.59 of Allen & Tildesley)
         ! ... the "rise time" is tau=nraise*dt so dt/tau=1/nraise
         ! ... Equivalent to traditional rescaling if nraise=1
         !
         IF ( system_temp > 0.D0 .AND. required_temp > 0.D0 ) THEN
            aux = SQRT( 1.d0 + (required_temp / system_temp - 1.d0) * &
                                (1.D0/DBLE(nraise) ) )
         ELSE
            aux = 0.d0
         ENDIF
         !
      ELSE
         !
         ! ... rescale the velocities by a factor 3 / 2KT / Ek
         !
         IF ( system_temp > 0.D0 .AND. required_temp > 0.D0 ) THEN
            aux = SQRT( required_temp / system_temp )
         ELSE
            aux = 0.d0
         ENDIF
         !
      ENDIF
      !
      vel(:,:) = vel(:,:) * aux
      !
   END SUBROUTINE thermalize
   !
   !
   !-----------------------------------------------------------------------
   SUBROUTINE smart_MC()
     !-----------------------------------------------------------------------
     !! Routine to apply smart_MC.  
     !! Implemented by Xiaochuan Ge, Jul., 2013
     !
     !! At this moment works only with Langevin dynamics.  
     !! For the formula see R.J.Rossky, JCP, 69, 4628(1978).
     ! 
     USE ions_base,           ONLY : nat, ityp, tau, if_pos,atm
     USE cell_base,           ONLY : alat
     USE ener,                ONLY : etot
     USE force_mod,           ONLY : force
     USE control_flags,       ONLY : istep, lconstrain
     USE constraints_module,  ONLY : remove_constr_force, check_constraint
     USE random_numbers,      ONLY : randy
     USE io_files,            ONLY : prefix
     USE constants,           ONLY : bohr_radius_angs
     !
     IMPLICIT NONE
     !
     LOGICAL :: accept
     REAL(DP) :: kt,sigma2,             &
                 T_ij,T_ji,boltzman_ji, &   ! boltzman_ji=exp[-(etot_new-etot_old)/kt]
                 temp,p_smc                 ! *_smart means *_old, the quantity of the
                                            ! previous step
     !
     INTEGER :: ia, ip
     !
     IF ( lconstrain ) THEN
        ! ... we first remove the component of the force along the
        ! ... constraint gradient ( this constitutes the initial
        ! ... guess for the calculation of the lagrange multipliers )
        CALL remove_constr_force( nat, tau, if_pos, ityp, alat, force )
     ENDIF
     !
     IF (first_iter) THEN ! For the first iteration
        ALLOCATE( tau_smart(3,nat) )
        ALLOCATE( force_smart(3,nat) )
        tau_smart = tau
        etot_smart = etot
        force_smart = force
        first_iter = .FALSE.
        RETURN
     ENDIF
     !
     kt = temperature / ry_to_kelvin
     sigma2 =  2.D0*dt*kt
     !
     T_ij=0.0d0
     T_ji=0.0d0
     DO ia = 1, nat
        DO ip = 1, 3
           T_ij = T_ij + ( (tau(ip,ia)-tau_smart(ip,ia))*alat-dt*force_smart(ip,ia) )**2
           T_ji = T_ji + ( (tau_smart(ip,ia)-tau(ip,ia))*alat-dt*force(ip,ia) )**2
        ENDDO
     ENDDO
     T_ij = EXP(-T_ij/(2*sigma2))
     T_ji = EXP(-T_ji/(2*sigma2))
     !
     boltzman_ji = EXP(-(etot-etot_smart)/kt)
     !
     p_smc = T_ji*boltzman_ji/T_ij
     !
     WRITE(stdout, '(5x,"The old energy is:",3x,F17.8," Ry")') etot_smart
     WRITE(stdout, '(5x,"The new energy is:",3x,F17.8," Ry")') etot
     WRITE(stdout, '(5x,"The possibility to accept this step is:",3x,F10.7/)') p_smc
     WRITE(stdout, '(5x,"Nervously waiting for the fate ..."/)')
     !
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
     !
     WRITE (stdout, '(5x,"The current acceptance is :",3x,F10.6)') DBLE(num_accept)/istep
     !
     ! Print the trajectory
     !
     OPEN(117,file="trajectory-"//TRIM(prefix)//".xyz", STATUS="unknown", POSITION='APPEND')
     WRITE(117,'(I5)') nat
     WRITE(117,'("# Step: ",I5,5x,"Total energy: ",F17.8,5x,"Ry")') istep-1, etot
     DO ia = 1, nat
        WRITE( 117, '(A3,3X,3F14.9)') atm(ityp(ia)),tau(:,ia)*alat*bohr_radius_angs
     ENDDO
     CLOSE(117)
     !
     RETURN
     !
   END SUBROUTINE smart_MC
   !
   !-----------------------------------------------------------------------
   SUBROUTINE thermalize_resamp_vscaling( nraise_, system_temp, required_temp )
      !-----------------------------------------------------------------------
      !! Sample velocities using stochastic velocity rescaling, based on:
      !! Bussi, Donadio, Parrinello, J. Chem. Phys. 126, 014101 (2007),
      !! doi: 10.1063/1.2408420
      !
      !! Implemented (2019) by Leonid Kahle and Ngoc Linh Nguyen,
      !! Theory and Simulations of Materials Laboratory, EPFL.
      !
      USE ions_base,          ONLY : nat, if_pos
      USE cell_base,          ONLY : alat
      USE random_numbers,     ONLY : gauss_dist, sum_of_gaussians2
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(in) :: system_temp, required_temp
      INTEGER,  INTENT(in) :: nraise_
      !
      INTEGER  :: i
      REAL(DP) :: factor, rr
      REAL(DP) :: aux, aux2
      real(DP), external :: gasdev, sumnoises
      INTEGER  :: na
      !
      ndof = get_ndof()
      !
      IF ( nraise_ > 0 ) THEN
         !
         ! ... the "rise time" is tau=nraise*dt so dt/tau=1/nraise
         ! ... Equivalent to traditional rescaling if nraise=1
         !
         factor = exp(-1.0/nraise_)
      ELSE
         !
         factor = 0.0
         !
      ENDIF
      !
      IF ( system_temp > 0.D0 .and. required_temp > 0.D0 ) THEN
         !
         ! Applying Eq. (A7) from J. Chem. Phys. 126, 014101 (2007)
         !
         rr = gauss_dist(0.0D0, 1.0D0)
         aux2 = factor + (1.0D0-factor)*( sum_of_gaussians2(ndof - 1) +rr**2) &
                * required_temp/(ndof * system_temp) &
                + 2*rr*sqrt((factor*(1.0D0-factor)*required_temp)/(ndof * system_temp))
         !
         aux  = sqrt(aux2)

      ELSE
         !
         aux = 0.d0
         !
      ENDIF
      !
      ! Global rescaling applied to velocities
      vel(:,:) = vel(:,:) * aux
      !
   END SUBROUTINE thermalize_resamp_vscaling
   !
   !-----------------------------------------------------------------------
   SUBROUTINE dump_trajectory_frame( time, temp )
      !-----------------------------------------------------------------------
      !! Dump trajectory frame into a file in tmp_dir with name:
      !! prefix.istep.mdtrj. Don't append, create a new file for each step,
      !! let the caller worry about merging. Safer for limited wall times,
      !! when the job can be killed anytime.
      !
      USE cell_base,   ONLY : alat, at
      USE constants,   ONLY : bohr_radius_angs
      USE ener,        ONLY : etot
      USE io_files,    ONLY : prefix, tmp_dir
      USE ions_base,   ONLY : nat, tau
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(in) :: time, temp ! in ps, K
      INTEGER              :: iunit, i, k
      !
      OPEN(NEWUNIT = iunit, FILE = TRIM( tmp_dir ) // TRIM( prefix ) // ".mdtrj", &
         STATUS="unknown", POSITION="APPEND")
      !
      ! Time (ps), temp (K), total energy (Ry), unit cell (9 values, Ang),
      ! atom coordinates (3 * nat values, Ang)
      !
      WRITE(iunit, *) time, temp, etot, &
         ( ( at(i,k) * alat * bohr_radius_angs, i = 1, 3), k = 1, 3 ), &
         ( ( tau(i,k) * alat * bohr_radius_angs, i = 1, 3), k = 1, nat )
      WRITE(iunit, *) ! new line
      !
      CLOSE(iunit)
      !
   END SUBROUTINE dump_trajectory_frame
   !
   !--------------------------------------------------------------------------
   SUBROUTINE verlet_read_tau_from_conf( )
      !-----------------------------------------------------------------------
      !! Try to set tau from restart prefix.md file.
      !
      USE io_files,  ONLY : prefix, seqopn
      USE ions_base, ONLY : nat, tau
      USE io_global, ONLY : ionode, ionode_id
      USE mp,        ONLY : mp_bcast
      USE mp_images, ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      LOGICAL  :: is_restart = .FALSE.
      INTEGER  :: restart_id = 0
      REAL(DP) :: etotold = 0.0_DP
      REAL(DP) :: tau_tmp(3, nat)
      INTEGER  :: istep
      !
      ! Try to read restart file on the ionode only
      IF ( ionode ) THEN
         !
         CALL seqopn( 4, 'md', 'FORMATTED', is_restart )
         !
         IF ( is_restart ) THEN
            !
            READ( UNIT = 4, FMT = * ) restart_id
            !
            IF ( restart_id .EQ. restart_verlet ) THEN
               !
               READ( UNIT = 4, FMT = * ) istep, etotold, tau_tmp(:,:)
               !
               IF ( SUM ( (tau_tmp(:,1:nat)-tau(:,1:nat))**2 ) > eps8 ) THEN
                  !
                  tau(:, 1:nat) = tau_tmp(:, 1:nat)
                  !
                  WRITE( stdout, '(/5X,"Atomic positions read from:", &
                     /,5X,A)') TRIM(prefix) // ".md"
               END IF
               !
            END IF
            !
            CLOSE( UNIT = 4 )
            !
         ELSE
            !
            CLOSE( UNIT = 4, STATUS = 'DELETE' )
            !
         END IF
      !
      END IF
      !
      CALL mp_bcast( tau, ionode_id, intra_image_comm )
      !
   END SUBROUTINE verlet_read_tau_from_conf
   !
   !-----------------------------------------------------------------------
   FUNCTION get_ndof()
      !-----------------------------------------------------------------------
      !! Get the number of degrees of freedom. Use number of constraints
      !! requested from the constraints_module.
      !
      USE ions_base,          ONLY : nat, if_pos
      USE constraints_module, ONLY : nconstr_ndof
      !
      IMPLICIT NONE
      !
      REAL(DP) :: get_ndof
      !
      ! ... the number of degrees of freedom
      !
      IF ( ANY( if_pos(:,:) == 0 ) ) THEN
         !
         get_ndof = 3*nat - count( if_pos(:,:) == 0 ) - nconstr_ndof
         !
      ELSE
         !
         get_ndof = 3*nat - 3 - nconstr_ndof
         !
      ENDIF
      !
   END FUNCTION get_ndof
   !
END MODULE dynamics_module
