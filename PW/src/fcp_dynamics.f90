!
! Copyright (C) 2001-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE fcp_dynamics
  !--------------------------------------------------------------------------
  !
  ! ... Original version by Minoru Otani (AIST) and Nicephore Bonnet (AIST).
  ! ...
  ! ... This module controls the Fictitious Charge Particle (FCP) for constant-mu
  ! ... method developed by N. Bonnet, T. Morishita, O. Sugino, and M. Otani
  ! ... (see PRL 109, 266101 [2012]).
  ! ...
  ! ... Constant-mu scheme with the boundary condition 'bc2' and 'bc3' enables
  ! ... description of the system connected to a potentiostat which preserves
  ! ... the Fermi energy of the system as the target Fermi energy (mu).
  ! ...
  ! ... Newton and BFGS algorithms are implemented by S. Nishihara (2016-2017)
  ! ...
  ! ... . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ! ...   This module performes dynamics of FCP.
  ! ... . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  !
  USE constants,      ONLY : RYTOEV, ry_to_kelvin
  USE control_flags,  ONLY : iverbosity
  USE ener,           ONLY : ef
  USE io_files,       ONLY : seqopn
  USE io_global,      ONLY : stdout
  USE ions_base,      ONLY : nat, ityp, zv
  USE kinds,          ONLY : DP
  USE klist,          ONLY : nelec, tot_charge
  USE random_numbers, ONLY : randy, gauss_dist, set_random_seed
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define parameters
  INTEGER, PARAMETER :: IDYN_NULL        = 0
  INTEGER, PARAMETER :: IDYN_VERLET      = 1
  INTEGER, PARAMETER :: IDYN_PROJ_VERLET = 2
  !
  CHARACTER(LEN=8), PARAMETER :: EXT_FCP = 'fcp'
  !
  ! ... define variables
  INTEGER           :: idyn          ! type of dynamics
  INTEGER           :: iter          ! number of iteration
  REAL(DP)          :: epsf          ! convergence threshold of projected-Verlet (in Ry)
  REAL(DP)          :: step_max      ! maximum step of projected-Verlet (in number of charge)
  REAL(DP)          :: nelec_old     ! old number of electrons
  REAL(DP)          :: vel           ! velocity of FCP
  REAL(DP)          :: vel_init      ! initial velocity of FCP
  REAL(DP)          :: acc           ! acceleration of FCP
  REAL(DP)          :: mass          ! mass of FCP
  LOGICAL           :: vel_verlet    ! if TRUE, Velocity-Verlet algorithm is used
  LOGICAL           :: vel_definit   ! if TRUE, initial velocity is defined
  LOGICAL           :: vel_defined   ! if TRUE, vel is used rather than nelec_old to do the next step
  LOGICAL           :: control_temp  ! if TRUE, a thermostat is used to control the temperature
  CHARACTER(LEN=10) :: thermostat    ! the thermostat used to control the temperature
  REAL(DP)          :: temperature   ! starting temperature
  REAL(DP)          :: tolp          ! tolerance for temperature variation
  REAL(DP)          :: delta_t       ! parameter used in thermalization
  INTEGER           :: nraise        ! parameter used in thermalization
  !
  ! ... public components
  PUBLIC :: fcpdyn_init
  PUBLIC :: fcpdyn_final
  PUBLIC :: fcpdyn_prm_mass
  PUBLIC :: fcpdyn_prm_velocity
  PUBLIC :: fcpdyn_prm_temp
  PUBLIC :: fcpdyn_set_verlet
  PUBLIC :: fcpdyn_set_velocity_verlet
  PUBLIC :: fcpdyn_set_proj_verlet
  PUBLIC :: fcpdyn_update
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcpdyn_init()
    !----------------------------------------------------------------------------
    !
    ! ... initialize this module
    !
    IMPLICIT NONE
    !
    idyn         = IDYN_NULL
    iter         = 0
    epsf         = 0.0_DP
    step_max     = 0.0_DP
    nelec_old    = 0.0_DP
    vel          = 0.0_DP
    vel_init     = 0.0_DP
    acc          = 0.0_DP
    mass         = 10000.0_DP
    vel_verlet   = .FALSE.
    vel_definit  = .FALSE.
    vel_defined  = .FALSE.
    control_temp = .FALSE.
    thermostat   = ''
    temperature  = 300.0_DP
    tolp         = 100.0_DP
    delta_t      = 1.0_DP
    nraise       = 1
    !
  END SUBROUTINE fcpdyn_init
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcpdyn_final()
    !----------------------------------------------------------------------------
    !
    ! ... finalize this module
    !
    IMPLICIT NONE
    !
    IF (idyn == IDYN_VERLET .AND. iter > 0) THEN
       !
       CALL fcp_dyn_printavg(iter)
       !
    END IF
    !
  END SUBROUTINE fcpdyn_final
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcpdyn_prm_mass(mass_)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: mass_
    !
    IF (mass_ > 0.0_DP) THEN
       mass = mass_
    END IF
    !
  END SUBROUTINE fcpdyn_prm_mass
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcpdyn_prm_velocity(vel_)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: vel_
    !
    vel_definit = .TRUE.
    vel_init    = vel_
    !
  END SUBROUTINE fcpdyn_prm_velocity
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcpdyn_prm_temp(thermostat_, tempw, tolp_, delta_t_, nraise_)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: thermostat_
    REAL(DP),         INTENT(IN) :: tempw
    REAL(DP),         INTENT(IN) :: tolp_
    REAL(DP),         INTENT(IN) :: delta_t_
    INTEGER,          INTENT(IN) :: nraise_
    !
    SELECT CASE (TRIM(thermostat_))
       !
    CASE ('not_controlled', 'not-controlled', 'not controlled')
       !
       control_temp = .FALSE.
       !
    CASE ('rescaling')
       !
       control_temp = .TRUE.
       thermostat = TRIM(thermostat_)
       !
       IF (tempw > 0.0_DP) THEN
          temperature = tempw
       END IF
       !
       IF (tolp_ >= 0.0_DP) THEN
          tolp = tolp_
       END IF
       !
    CASE ('rescale-v', 'rescale-V', 'rescale_v', 'rescale_V')
       !
       control_temp = .TRUE.
       thermostat = TRIM(thermostat_)
       !
       IF (tempw > 0.0_DP) THEN
          temperature = tempw
       END IF
       !
       IF (nraise_ > 0) THEN
          nraise = nraise_
       END IF
       !
    CASE ('rescale-T', 'rescale-t', 'rescale_T', 'rescale_t')
       !
       control_temp = .TRUE.
       thermostat = TRIM(thermostat_)
       !
       IF (tempw > 0.0_DP) THEN
          temperature = tempw
       END IF
       !
       IF (delta_t_ > 0.0_DP) THEN
          delta_t = delta_t_
       END IF
       !
    CASE ('reduce-T', 'reduce-t', 'reduce_T', 'reduce_t')
       !
       control_temp = .TRUE.
       thermostat = TRIM(thermostat_)
       !
       IF (tempw > 0.0_DP) THEN
          temperature = tempw
       END IF
       !
       IF (nraise_ > 0) THEN
          nraise = nraise_
       END IF
       !
       IF (delta_t_ < 0.0_DP) THEN
          delta_t = delta_t_
       END IF
       !
    CASE ('berendsen', 'Berendsen')
       !
       control_temp = .TRUE.
       thermostat = TRIM(thermostat_)
       !
       IF (tempw > 0.0_DP) THEN
          temperature = tempw
       END IF
       !
       IF (nraise_ > 0) THEN
          nraise = nraise_
       END IF
       !
    CASE ('andersen', 'Andersen')
       !
       control_temp = .TRUE.
       thermostat = TRIM(thermostat_)
       !
       IF (tempw > 0.0_DP) THEN
          temperature = tempw
       END IF
       !
       IF (nraise_ > 0) THEN
          nraise = nraise_
       END IF
       !
    CASE ('initial', 'Initial')
       !
       control_temp = .TRUE.
       thermostat = TRIM(thermostat_)
       !
       IF (tempw > 0.0_DP) THEN
          temperature = tempw
       END IF
       !
    CASE DEFAULT
       !
       CALL errore('iosys', 'unknown fcp_temperature ' // TRIM(thermostat), 1)
       !
    END SELECT
    !
  END SUBROUTINE fcpdyn_prm_temp
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcpdyn_set_verlet()
    !----------------------------------------------------------------------------
    !
    ! ... set the type of dynamics to Verlet
    !
    IMPLICIT NONE
    !
    idyn = IDYN_VERLET
    !
    vel_verlet = .FALSE.
    !
  END SUBROUTINE fcpdyn_set_verlet
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcpdyn_set_velocity_verlet()
    !----------------------------------------------------------------------------
    !
    ! ... set the type of dynamics to Velocity-Verlet
    !
    IMPLICIT NONE
    !
    idyn = IDYN_VERLET
    !
    vel_verlet = .TRUE.
    !
  END SUBROUTINE fcpdyn_set_velocity_verlet
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcpdyn_set_proj_verlet(eps, smax)
    !----------------------------------------------------------------------------
    !
    ! ... set the type of dynamics to projected-Verlet
    ! ...
    ! ... Variables:
    ! ...   eps:  convergende threshold (in Ry/e)
    ! ...   smax: maximum step (in e)
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: eps
    REAL(DP), INTENT(IN) :: smax
    !
    IF (eps < 0.0_DP) THEN
       !
       CALL errore('fcpdyn_set_proj_verlet', 'eps is negative', 1)
       !
    END IF
    !
    IF (smax <= 0.0_DP) THEN
       !
       CALL errore('fcpdyn_set_proj_verlet', 'smax is not positive', 1)
       !
    END IF
    !
    idyn = IDYN_PROJ_VERLET
    !
    epsf     = eps
    step_max = smax
    !
  END SUBROUTINE fcpdyn_set_proj_verlet
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcpdyn_update(mu, dt, conv)
    !----------------------------------------------------------------------------
    !
    ! ... update number of electrons
    ! ...
    ! ... Variables:
    ! ...   mu:   target Fermi energy (in Ry)
    ! ...   dt:   time step (in a.u.)
    ! ...   conv: converged, or not ?
    !
    IMPLICIT NONE
    !
    REAL(DP),          INTENT(IN)    :: mu
    REAL(DP),          INTENT(IN)    :: dt
    LOGICAL, OPTIONAL, INTENT(INOUT) :: conv
    !
    LOGICAL  :: conv_
    REAL(DP) :: force
    REAL(DP) :: tot_charge_
    REAL(DP) :: ekin, temp
    !
    conv_ = .FALSE.
    IF (PRESENT(conv)) THEN
      conv_ = conv
    END IF
    !
    ! ... set variables
    !
    force = mu - ef
    !
    tot_charge_ = tot_charge
    !
    ! ... update nelec
    !
    IF (idyn == IDYN_VERLET) THEN
       !
       CALL verlet(mu, dt)
       !
    ELSE IF (idyn == IDYN_PROJ_VERLET) THEN
       !
       CALL proj_verlet(mu, dt, conv_)
       !
    ELSE
       !
       CALL errore('fcpdyn_update', 'idyn is incorrect', 1)
       !
    END IF
    !
    ! ... update tot_charge
    !
    tot_charge = SUM(zv(ityp(1:nat))) - nelec
    !
    ! ... eval temperature
    !
    ekin = 0.5_DP * mass * vel * vel
    !
    temp = 2.0_DP * ekin * ry_to_kelvin
    !
    ! ... write information
    !
    IF (.NOT. conv_) THEN
       WRITE(stdout, '(/,5X,"FCP: iteration #",I5)') iter
       WRITE(stdout, '(  5X,"FCP: Total Charge = ",F12.6,"  -> ",F12.6)') tot_charge_, tot_charge
    ELSE
       WRITE(stdout, '(/,5X,"FCP: Total Charge = ",F12.6)') tot_charge_
    END IF
    !
    WRITE(stdout, '(5X,"FCP: Velocity     = ",1PE12.2," a.u.")') vel
    WRITE(stdout, '(5X,"FCP: Acceleration = ",1PE12.2," a.u.")') acc
    WRITE(stdout, '(5X,"FCP: Temperature  = ",F12.3," K")') temp
    WRITE(stdout, '(5X,"FCP: Fermi Energy = ",F12.6," Ry (",F12.6," eV)")') ef,    ef    * RYTOEV
    WRITE(stdout, '(5X,"FCP: Target Level = ",F12.6," Ry (",F12.6," eV)")') mu,    mu    * RYTOEV
    WRITE(stdout, '(5X,"FCP: Force on FCP = ",F12.6," Ry (",F12.6," eV)")') force, force * RYTOEV
    !
    IF (idyn == IDYN_PROJ_VERLET) THEN
       WRITE(stdout, '(5X,"FCP: Force Thr.   = ",F12.6," Ry (",F12.6," eV)")') epsf, epsf * RYTOEV
    END IF
    !
    WRITE(stdout, '(/)')
    !
    IF (PRESENT(conv)) THEN
      conv = conv_
    END IF
    !
  END SUBROUTINE fcpdyn_update
  !
  !----------------------------------------------------------------------------
  SUBROUTINE verlet(mu, dt)
    !----------------------------------------------------------------------------
    !
    ! ... perform one step of molecular dynamics evolution, using the Verlet algorithm.
    ! ... if vel_verlet is .TRUE., Velocity-Verlet algorithm is used.
    ! ... [NOTE: original code is from PW/src/dynamics_module.f90]
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: mu
    REAL(DP), INTENT(IN) :: dt
    !
    REAL(DP) :: ekin
    REAL(DP) :: force
    REAL(DP) :: temp_av
    REAL(DP) :: temp_new
    REAL(DP) :: nelec_new
    REAL(DP) :: vel_mod
    INTEGER  :: iun
    LOGICAL  :: leof
    LOGICAL  :: file_exists
    !
    INTEGER, EXTERNAL :: find_free_unit
    !
    vel_defined = .TRUE.
    temp_av     = 0.0_DP
    temp_new    = 0.0_DP
    vel_mod     = 0.0_DP
    !
    ! ... check and read the file (*.fcp)
    !
    iun = find_free_unit()
    CALL seqopn(iun, TRIM(EXT_FCP), 'FORMATTED', file_exists)
    !
    IF (file_exists) THEN
       !
       ! ... the file is read :  simulation is continuing
       !
       READ(UNIT = iun, FMT = *) iter, nelec_old, leof
       !
       IF (leof) THEN
          !
          ! ... the file was created by projected_verlet:  Ignore it
          !
          CALL md_init(temp_new, temp_av)
          !
       ELSE
          !
          vel_defined = .FALSE.
          !
          READ(UNIT = iun, FMT = *) vel_mod, temp_new, temp_av, mass
          !
       ENDIF
       !
       CLOSE(UNIT = iun, STATUS = 'KEEP')
       !
    ELSE
       !
       CLOSE(UNIT = iun, STATUS = 'DELETE')
       !
       ! ... the file is absent :  simulation is starting from scratch
       !
       CALL md_init(temp_new, temp_av)
       !
    END IF
    !
    ! ... update counter
    !
    iter = iter + 1
    !
    ! ... calculate force and acceleration
    !
    force = mu - ef
    !
    acc   = force / mass
    !
    IF (vel_verlet) THEN
       !
       ! >>> Velocity-Verlet integration scheme
       !
       IF (.NOT. vel_defined) THEN
          !
          ! ... update velocity
          !
          vel = vel_mod + 0.5_DP * acc * dt
          !
          ! ... calculate kinetic energy
          !
          ekin = 0.5_DP * mass * vel * vel
          !
          ! ... find the new temperature and update the average
          !
          temp_new = 2.0_DP * ekin * ry_to_kelvin
          !
          temp_av  = temp_av + temp_new
          !
       END IF
       !
       ! ... control the temperature
       !
       IF (control_temp) THEN
          !
          CALL apply_thermostat(dt, temp_new, temp_av, .TRUE.)
          !
       END IF
       !
       ! ... calculate new number of electrons
       !
       nelec_new = nelec + vel * dt + 0.5_DP * acc * dt * dt
       !
    ELSE
       !
       ! >>> Verlet integration scheme
       !
       ! ... control the temperature
       !
       IF (control_temp) THEN
          !
          CALL apply_thermostat(dt, temp_new, temp_av, vel_defined)
          !
       END IF
       !
       ! ... calculate new number of electrons
       !
       IF (vel_defined) THEN
          !
          nelec_new = nelec + vel * dt + 0.5_DP * acc * dt * dt
          nelec_old = nelec - vel * dt + 0.5_DP * acc * dt * dt
          !
       ELSE
          !
          nelec_new = 2.0_DP * nelec - nelec_old + acc * dt * dt
          !
       ENDIF
       !
       ! ... update velocity
       !
       vel = (nelec_new - nelec_old) / (2.0_DP * dt)
       !
       ! ... calculate kinetic energy
       !
       ekin = 0.5_DP * mass * vel * vel
       !
       ! ... find the new temperature and update the average
       !
       temp_new = 2.0_DP * ekin * ry_to_kelvin
       !
       temp_av  = temp_av + temp_new
       !
    END IF
    !
    ! ... modified velocity for Velocity-Verlet
    !
    vel_mod = vel + 0.5_DP * acc * dt
    !
    ! ... save all the needed quantities on file
    !
    CALL seqopn(iun, TRIM(EXT_FCP), 'FORMATTED', file_exists)
    !
    leof = .FALSE.
    WRITE(UNIT = iun, FMT = *) iter, nelec, leof
    !
    WRITE(UNIT = iun, FMT = *) vel_mod, temp_new, temp_av, mass
    !
    CLOSE(UNIT = iun, STATUS = 'KEEP')
    !
    ! ... evaluate some properties
    !
    CALL fcp_dyn_calcavg(iter, nelec, vel, acc, force, mass, temp_new)
    !
    ! ... here the nelec are shifted
    !
    nelec = nelec_new
    !
  END SUBROUTINE verlet
  !
  !----------------------------------------------------------------------------
  SUBROUTINE md_init(temp_new, temp_av)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: temp_new
    REAL(DP), INTENT(OUT) :: temp_av
    !
    REAL(DP) :: ekin
    !
    ! ... print initial message
    !
    WRITE(stdout, '(/,5X,"FCP Dynamics Calculation")')
    !
    IF (control_temp) THEN
       !
       WRITE(stdout, '(/,5X,"Starting temperature",T27," = ",F8.2," K")' ) temperature
       !
       SELECT CASE (TRIM(thermostat))
          !
       CASE ('andersen', 'Andersen')
          !
          WRITE(stdout, '(/,5X,"temperature is controlled by Andersen ", &
                           &   "thermostat",/,5X,"Collision frequency =",&
                           &    F7.4,"/timestep")' ) 1.0_DP / DBLE(nraise)
          !
       CASE ('berendsen', 'Berendsen')
          !
          WRITE(stdout, '(/,5X,"temperature is controlled by soft ", &
                           &   "(Berendsen) velocity rescaling",/,5X,&
                           &   "Characteristic time =",I3,"*timestep")') nraise
          !
       CASE ('initial', 'Initial')
          !
          WRITE(stdout, '(/,5X,"temperature is set once at start")')
          !
       CASE DEFAULT
          !
          WRITE(stdout, '(/,5X,"temperature is controlled by ", &
                           &   "velocity rescaling (",A,")")' ) TRIM(thermostat)
          !
       END SELECT
       !
    END IF
    !
    IF (vel_verlet) THEN
       WRITE(stdout, '(/,5X,"FCP: Velocity-Verlet Algorithm is used.")')
    ELSE
       WRITE(stdout, '(/,5X,"FCP: Verlet Algorithm is used.")')
    END IF
    !
    WRITE(stdout, '(5X,"FCP: Mass of FCP  = ",1PE12.2," a.u.")') mass
    !
    ! ... initialize counter
    !
    iter = 0
    !
    ! ... initialize velocity
    !
    IF (vel_definit) THEN
       !
       ! ... initial velocities available from input file
       !
       vel = vel_init
       vel_defined = .TRUE.
       !
    ELSE IF (control_temp) THEN
       !
       ! ... initial thermalization.
       !
       CALL start_therm()
       vel_defined = .TRUE.
       !
    ELSE
       !
       vel = 0.0_DP
       vel_defined = .TRUE.
       !
    END IF
    !
    ! ... calculate initial temperature
    !
    ekin = 0.5_DP * mass * vel * vel
    !
    temp_new = 2.0_DP * ekin * ry_to_kelvin
    !
    temp_av  = temp_new
    !
  END SUBROUTINE md_init
  !
  !----------------------------------------------------------------------------
  SUBROUTINE apply_thermostat(dt, temp_new, temp_av, has_vel)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)    :: dt
    REAL(DP), INTENT(IN)    :: temp_new
    REAL(DP), INTENT(INOUT) :: temp_av
    LOGICAL,  INTENT(IN)    :: has_vel
    !
    REAL(DP) :: kt
    REAL(DP) :: sigma
    LOGICAL  :: lmoved
    !
    IF (.NOT. has_vel) THEN
      vel = (nelec - nelec_old) / dt
    END IF
    !
    SELECT CASE (TRIM(thermostat))
       !
    CASE ('rescaling')
       !
       IF (ABS(temp_new - temperature) > tolp) THEN
          !
          WRITE(stdout, '(/,5X,"FCP: Velocity rescaling: T (",F6.1,"K) ", &
                           &   "out of range, reset to " ,F6.1)' ) temp_new, temperature
          !
          CALL thermalize(0, temp_new, temperature)
          !
       END IF
       !
    CASE ('rescale-v', 'rescale-V', 'rescale_v', 'rescale_V')
       !
       IF (MOD(iter, nraise) == 0) THEN
          !
          temp_av = temp_av / DBLE(nraise)
          !
          WRITE(stdout, '(/,5X,"FCP: Velocity rescaling: average T on ",I3, &
                           &   " steps (",F6.1,"K) reset to ",F6.1)') &
                           &   nraise, temp_av, temperature
          !
          CALL thermalize(0, temp_new, temperature)
          !
          temp_av = 0.0_DP
          !
       END IF
       !
    CASE ('rescale-T', 'rescale-t', 'rescale_T', 'rescale_t')
       !
       IF (delta_t > 0.0_DP) THEN
          !
          temperature = temp_new * delta_t
          !
          WRITE(stdout, '(/,5X,"FCP: Thermalization: T (",F6.1,"K) rescaled ",&
                           &   "by a factor ",F6.3)' ) temp_new, delta_t
          !
          CALL thermalize(0, temp_new, temperature)
          !
       END IF
       !
    CASE ('reduce-T', 'reduce-t', 'reduce_T', 'reduce_t')
       !
       IF (MOD(iter, nraise) == 0 .AND. delta_t < 0.0_DP) THEN
          !
          temperature = temp_new + delta_t
          !
          WRITE(stdout, '(/,5X,"FCP: Thermalization: T (",F6.1,"K) reduced ",&
                           &   "by ",F6.3)' ) temp_new, -delta_t
          !
          CALL thermalize(0, temp_new, temperature)
          !
       END IF
       !
    CASE ('berendsen', 'Berendsen')
       !
       WRITE(stdout, '(/,5X,"FCP: Soft (Berendsen) velocity rescaling")')
       !
       CALL thermalize(nraise, temp_new, temperature)
       !
    CASE ('andersen', 'Andersen')
       !
       kt = temperature / ry_to_kelvin
       lmoved = .FALSE.
       !
       IF (randy() < 1.0_DP / DBLE(nraise)) THEN
          !
          lmoved = .TRUE.
          sigma  = SQRT(kt / mass)
          vel    = gauss_dist(0.0_DP, sigma)
          !
       END IF
       !
       IF (lmoved) THEN
          WRITE(stdout, '(/,5X,"FCP Andersen thermostat: a collision has done")')
       END IF
       !
    CASE ('initial', 'Initial')
       !
       ! NOP
       !
    END SELECT
    !
    ! ... the old number of electron is updated to reflect the new velocities,
    ! ... only in the first
    !
    IF (.NOT. has_vel) THEN
       nelec_old = nelec - vel * dt
    END IF
    !
  END SUBROUTINE apply_thermostat
  !
  !----------------------------------------------------------------------------
  SUBROUTINE start_therm()
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP) :: kt
    REAL(DP) :: sigma
    !
    ! ... next command prevents different MD runs to start
    ! ... with exactly the same "random" velocities
    !
    CALL set_random_seed()
    !
    ! ... starting velocity
    !
    kt = temperature / ry_to_kelvin
    sigma = SQRT(kt / mass)
    !
    IF (randy() <= 0.5_DP) THEN
       vel = +sigma
    ELSE
       vel = -sigma
    END IF
    !
  END SUBROUTINE start_therm
  !
  !----------------------------------------------------------------------------
  SUBROUTINE thermalize(nraise, system_temp, required_temp)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(in) :: nraise
    REAL(DP), INTENT(in) :: system_temp
    REAL(DP), INTENT(in) :: required_temp
    !
    REAL(DP) :: aux
    !
    IF (nraise > 0) THEN
       !
       ! ... Berendsen rescaling (Eq. 7.59 of Allen & Tildesley)
       ! ... the "rise time" is tau=nraise*dt so dt/tau=1/nraise
       ! ... Equivalent to traditional rescaling if nraise=1
       !
       IF (system_temp > 0.0_DP .AND. required_temp > 0.0_DP) THEN
          !
          aux = SQRT(1.0_DP + (required_temp / system_temp - 1.0_DP) * (1.0_DP / DBLE(nraise)))
          !
       ELSE
          !
          aux = 0.0_DP
          !
       END IF
       !
    ELSE
       !
       ! ... rescale the velocities by a factor 1 / 2KT / Ek
       !
       IF (system_temp > 0.0_DP .AND. required_temp > 0.0_DP) THEN
          !
          aux = SQRT(required_temp / system_temp)
          !
       ELSE
          !
          aux = 0.0_DP
          !
       END IF
       !
    END IF
    !
    vel = vel * aux
    !
  END SUBROUTINE thermalize
  !
  !----------------------------------------------------------------------------
  SUBROUTINE proj_verlet(mu, dt, conv)
    !----------------------------------------------------------------------------
    !
    ! ... perform one step of FCP's relaxation using
    ! ... the projected-Verlet algorithm.
    ! ... [NOTE: original code is from PW/src/dynamics_module.f90]
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)    :: mu
    REAL(DP), INTENT(IN)    :: dt
    LOGICAL,  INTENT(INOUT) :: conv
    !
    REAL(DP) :: force
    REAL(DP) :: nelec_new
    REAL(DP) :: step
    REAL(DP) :: norm_step
    REAL(DP) :: sign_step
    REAL(DP) :: ionic_charge
    INTEGER  :: iun
    LOGICAL  :: leof
    LOGICAL  :: file_exists
    !
    INTEGER, EXTERNAL :: find_free_unit
    !
    nelec_old = nelec
    !
    ! ... check and read the file (*.fcp)
    !
    iun = find_free_unit()
    CALL seqopn(iun, TRIM(EXT_FCP), 'FORMATTED', file_exists)
    !
    IF (file_exists) THEN
       !
       ! ... the file is read
       !
       READ(UNIT = iun, FMT = *) iter, nelec_old
       !
       CLOSE(UNIT = iun, STATUS = 'KEEP')
       !
    ELSE
       !
       CLOSE(UNIT = iun, STATUS = 'DELETE')
       !
       iter = 0
       !
       WRITE(stdout, '(/,5X,"Damped FCP Dynamics Calculation")')
       WRITE(stdout, '(/,5X,"FCP: Mass of FCP  = ",1PE12.2," a.u.")') mass
       !
    END IF
    !
    ! ... update counter
    !
    iter = iter + 1
    !
    ! ... calculate force and acceleration
    !
    force = mu - ef
    !
    acc   = force / mass
    !
    ! ... check if convergence for FCP minimization is achieved
    !
    conv = conv .AND. (ABS(force) < epsf)
    !
    IF (conv) THEN
       !
       WRITE(stdout, '(/,5X,"Damped FCP Dynamics: convergence achieved in " &
                      & ,I3," steps")' ) (iter - 1)
       WRITE(stdout, '(/,5X,"End of damped FCP dynamics calculation")' )
       !
       RETURN
       !
    END IF
    !
    ! ... Damped dynamics ( based on the projected-Verlet algorithm )
    !
    vel = nelec - nelec_old
    !
    CALL project_velocity()
    !
    step = vel + dt * dt * acc
    !
    norm_step = ABS(step)
    !
    IF (norm_step > 0.0_DP) THEN
       sign_step = step / norm_step
    ELSE
       sign_step = 0.0_DP
    END IF
    !
    nelec_new = nelec + sign_step * MIN(norm_step, step_max)
    !
    ! ... save on file all the needed quantities
    !
    CALL seqopn(iun, TRIM(EXT_FCP), 'FORMATTED', file_exists)
    !
    leof = .TRUE.
    WRITE(UNIT = iun, FMT = *) iter, nelec, leof
    !
    CLOSE(UNIT = iun, STATUS = 'KEEP')
    !
    ! ... here the nelec are shifted
    !
    IF (iverbosity > 0) THEN
       !
       ionic_charge = SUM(zv(ityp(1:nat)))
       !
       WRITE(stdout,'(5X,"FCP: Original charge = ",F12.6)') ionic_charge - nelec
       WRITE(stdout,'(5X,"FCP: Expected charge = ",F12.6)') ionic_charge - (nelec + step)
       WRITE(stdout,'(5X,"FCP: Next charge     = ",F12.6)') ionic_charge - nelec_new
       !
    END IF
    !
    nelec = nelec_new
    !
  END SUBROUTINE proj_verlet
  !
  !----------------------------------------------------------------------------
  SUBROUTINE project_velocity()
    !----------------------------------------------------------------------------
    !
    ! ... quick-min algorithm
    !
    IMPLICIT NONE
    !
    REAL(DP) :: norm_acc
    REAL(DP) :: sign_acc
    REAL(DP) :: projection
    !
    IF (iter <= 1) RETURN
    !
    norm_acc = ABS(acc)
    !
    IF (norm_acc > 0.0_DP) THEN
       sign_acc = acc / norm_acc
    ELSE
       sign_acc = 0.0_DP
    END IF
    !
    projection = vel * sign_acc
    !
    IF (projection < 0.0_DP) THEN
       !
       WRITE(stdout, '(/,5X,"FCP: velocity and acceleration are opposite to each other")')
       WRITE(stdout, '(  5X,"FCP: -> velocity is set to zero")')
       !
       vel = 0.0_DP
       !
    END IF
    !
  END SUBROUTINE project_velocity
  !
END MODULE fcp_dynamics
