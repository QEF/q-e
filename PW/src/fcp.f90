!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE fcp
   !----------------------------------------------------------------------------
   !! Original version by Minoru Otani (AIST) and Nicephore Bonnet (AIST).
   !
   !! This module controls the Fictitious Charge Particle (FCP) for constant-mu
   !! method developed by N. Bonnet, T. Morishita, O. Sugino, and M. Otani
   !! (see PRL 109, 266101 [2012]).
   !
   !! Constant-mu scheme with the boundary condition 'bc2' and 'bc3' enables
   !! description of the system connected to a potentiostat which preserves
   !! the Fermi energy of the system as the target Fermi energy (mu).
   !
   !! MDIIS algorithm is implemented for relax-calculation by
   !! Satomichi Nishihara (2016).
   !
   USE kinds,             ONLY: DP
   USE io_global,         ONLY: stdout
   USE io_files,          ONLY: seqopn
   USE constants,         ONLY: amu_ry, ry_to_kelvin, au_ps, rytoev, e2, fpi
   USE control_flags,     ONLY: tolp
   USE dynamics_module,   ONLY: dt, delta_t, nraise, &
                                control_temp, thermostat
   USE fcp_variables,     ONLY: fcp_mu, fcp_mass, fcp_temperature, fcp_relax,   &
                                fcp_relax_step, fcp_relax_crit, fcp_mdiis_size, &
                                fcp_mdiis_step
   USE mdiis,             ONLY: mdiis_type, allocate_mdiis, deallocate_mdiis,   &
                                update_by_mdiis
   !
   IMPLICIT NONE
   !
   SAVE
   !
   PRIVATE 
   ! 
   LOGICAL :: vel_defined
   !! if .TRUE., vel is used rather than tau_old to do the next step
   REAL(DP) :: tau
   !! present position
   REAL(DP) :: tau_old
   !! previous step position
   REAL(DP) :: tau_new
   !! next step position
   REAL(DP) :: vel
   !! velocity
   REAL(DP) :: acc
   !! acceleration
   !
   ! ... data of MDIIS
   LOGICAL          :: init_mdiis = .FALSE.
   !! 
   TYPE(mdiis_type) :: mdiist
   !! 
   !
   PUBLIC :: fcp_line_minimisation, fcp_verlet, fcp_summary
   PUBLIC :: fcp_mdiis_update, fcp_mdiis_end
   !
 CONTAINS
   !
   ! ... public methods
   !
   !------------------------------------------------------------------------------
   SUBROUTINE fcp_verlet()
      !-----------------------------------------------------------------------------
      !! This routine performs one step of molecular dynamics evolution by
      !! using the Verlet algorithm.
      !
      !! Parameters (in \(\textrm{dynamics_module}\) and \(\textrm{fcp_variables}\)):
      !
      !! * mass: mass of the atoms;
      !! * dt: time step;
      !! * temperature: starting temperature. The starting velocities of atoms
      !!   are set accordingly to the starting temperature, in random directions.
      !!   Therefore the initial velocity distribution is a constant.
      !
      !! Dario Alfe' 1997  and  Carlo Sbraccia 2004-2006
      !
      USE cell_base,      ONLY: alat
      USE control_flags,  ONLY: istep
      USE ener,           ONLY: ef
      USE klist,          ONLY: nelec, tot_charge
      USE ions_base,      ONLY: nat, ityp, zv
      !
      IMPLICIT NONE
      !
      REAL(DP) :: ekin, tau
      REAL(DP) :: temp_new, temp_av, elapsed_time, force
      ! istep counts all MD steps, including those of previous runs
      LOGICAL  :: file_exists
      !
      tau = nelec
      !
      vel_defined  = .TRUE.
      temp_av = 0.0_DP
      !
      CALL seqopn( 4, 'fcp_md', 'FORMATTED', file_exists )
      !
      IF ( file_exists ) THEN
         !
         ! ... the file is read:  simulation is continuing
         !
         READ( UNIT = 4, FMT = * ) istep, tau_old
         !
         vel_defined = .FALSE.
         !
         READ( UNIT = 4, FMT = * ) temp_new, temp_av, fcp_mass, elapsed_time
         !
         CLOSE( UNIT = 4, STATUS = 'KEEP' )
         !
      ELSE
         !
         CLOSE( UNIT = 4, STATUS = 'DELETE' )
         !
         ! ... the file is absent:  simulation is starting from scratch
         !
         CALL md_init()
         !
      ENDIF
      !
      ! ... elapsed_time is in picoseconds
      !
      elapsed_time = elapsed_time + dt*2.0_DP*au_ps
      istep = istep + 1
      !
      IF ( control_temp ) CALL apply_thermostat()
      !
      ! ... we first remove the component of the force along the
      ! ... constraint gradient ( this constitutes the initial
      ! ... guess for the calculation of the lagrange multipliers )
      !
      force = fcp_mu - ef
      acc = force / fcp_mass / alat
      !
      ! ... Verlet integration scheme
      !
      IF (vel_defined) THEN
         !
         tau_new = tau + vel * dt + 0.5_DP * acc * dt**2
         tau_old = tau - vel * dt + 0.5_DP * acc * dt**2
         !
      ELSE
         !
         tau_new = 2.0_DP*tau - tau_old + acc * dt**2
         !
      ENDIF
      !
      ! ... the linear momentum and the kinetic energy are computed here
      !
      vel = ( tau_new - tau_old ) / ( 2.0_DP*dt )
      !
      ekin = 0.5_DP * fcp_mass * vel**2 
      !
      ekin = ekin*alat**2
      !
      ! ... find the new temperature and update the average
      !
      temp_new = 2.0_DP * ekin * ry_to_kelvin
      !
      temp_av = temp_av + temp_new
      !
      !
      ! ... save all the needed quantities on file
      !
      CALL seqopn( 4, 'fcp_md', 'FORMATTED',  file_exists )
      !
      WRITE( UNIT = 4, FMT = * ) istep, tau
      !
      WRITE( UNIT = 4, FMT = * ) temp_new, temp_av, fcp_mass, elapsed_time
      !
      CLOSE( UNIT = 4, STATUS = 'KEEP' )
      !
      ! ... here the tau are shifted
      !
      tau = tau_new
      nelec = tau
      tot_charge = SUM( zv(ityp(1:nat)) ) - nelec
      !
      WRITE( stdout, '(/,5X,"FCP : Fermi Energy = ",F12.6," eV")') ef     * rytoev
      WRITE( stdout, '(  5X,"FCP : Target Mu    = ",F12.6," eV")') fcp_mu * rytoev
      WRITE( stdout, '(  5X,"FCP : tot_charge   = ",F12.6      )') tot_charge
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
         ! ... reference positions
         !
         IF ( control_temp ) THEN
            !
            WRITE( stdout, &
                  '(/,5X,"Starting fcp temperature",T27," = ",F8.2," K")' ) &
               fcp_temperature
            !
            SELECT CASE( TRIM( thermostat ) )
            CASE( 'andersen', 'Andersen' )
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"fcp temperature is controlled by Andersen ", &
                              &   "thermostat",/,5x,"Collision frequency =",&
                              &    f7.4,"/timestep")' ) 1.0_dp/nraise
               !
            CASE( 'berendsen', 'Berendsen' )
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"fcp temperature is controlled by soft ", &
                            &     "(Berendsen) velocity rescaling",/,5x,&
                            &     "Characteristic time =",i3,"*timestep")') &
                               nraise
               !
            CASE( 'initial', 'Initial' )
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"fcp temperature is set once at start"/)' )
               !
            CASE DEFAULT
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"fcp temperature is controlled by ", &
                              &     "velocity rescaling (",A,")"/)' )&
                              TRIM( thermostat )
               !
            END SELECT
            !
         ENDIF
         !
         WRITE( UNIT = stdout, FMT = '(5X,"fcp_mass = ",F8.2)' ) fcp_mass
         !
         WRITE( UNIT = stdout, &
                FMT = '(5X,"Time step",T27," = ",F8.2," a.u.,",F8.4, &
                        & " femto-seconds")' ) dt, dt*2.D+3*au_ps
         !
         ! ... masses in rydberg atomic units
         !
         IF ( control_temp ) THEN
            !
            ! ... initial thermalization. N.B. tau is in units of alat
            !
            CALL start_therm()
            vel_defined = .TRUE.
            !
            temp_new = fcp_temperature
            !
            temp_av = 0.0_DP
            !
         ELSE
            !
            vel = 0.0_DP
            vel_defined = .TRUE.
            !
         ENDIF
         !
         elapsed_time = 0.0_DP
         !
      END SUBROUTINE md_init
      !
      !
      !--------------------------------------------------------------------
      SUBROUTINE apply_thermostat()
         !--------------------------------------------------------------------
         !
         USE random_numbers,   ONLY: randy, gauss_dist
         !
         IMPLICIT NONE
         !
         INTEGER :: nat_moved
         REAL(DP) :: sigma, kt
         !
         IF (.NOT. vel_defined) THEN
            vel = (tau - tau_old) / dt
         ENDIF
         !
         SELECT CASE( TRIM( thermostat ) )
         CASE( 'rescaling' )
            !
            IF ( ABS(temp_new-fcp_temperature) > tolp ) THEN
               !
               WRITE( UNIT = stdout, &
                      FMT = '(/,5X,"FCP Velocity rescaling: T (",F6.1,"K) ", &
                                 & "out of range, reset to " ,F6.1)' ) &
                                   temp_new, fcp_temperature
               CALL thermalize( 0, temp_new, fcp_temperature )
               !
            ENDIF
            !
         CASE( 'rescale-v', 'rescale-V', 'rescale_v', 'rescale_V' )
            !
            IF ( MOD(istep, nraise) == 0 ) THEN
               !
               temp_av = temp_av / DBLE( nraise )
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"FCP Velocity rescaling: average T on ",i3, &
                                 &" steps (",F6.1,"K) reset to ",F6.1)' )  &
                           nraise, temp_av, fcp_temperature
               CALL thermalize( 0, temp_new, fcp_temperature )
               !
               temp_av = 0.0_DP
               !
            ENDIF
            !
         CASE( 'rescale-T', 'rescale-t', 'rescale_T', 'rescale_t' )
            ! Clearly it makes sense to check for positive delta_t
            ! If a negative delta_t is given, I suggest to have a message
            ! printed, that delta_t is ignored (TODO)
            IF ( delta_t > 0 ) THEN
               !
               fcp_temperature = fcp_temperature*delta_t
               !
               WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"Thermalization: T (",F6.1,"K) rescaled ",&
                                 & "by a factor ",F6.3)' ) temp_new, delta_t
               !
               CALL thermalize( 0, temp_new, fcp_temperature )
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
                                 & "by ",F6.3)' ) fcp_temperature, delta_t

               ELSE
                 WRITE( UNIT = stdout, &
                     FMT = '(/,5X,"Thermalization: T (",F6.1,"K) reduced ",&
                                 & "by ",F6.3)' ) fcp_temperature, -delta_t
               ENDIF
               ! I check whether the temperature is negative, so that I avoid
               ! nonsensical behavior:
               IF (fcp_temperature < 0.0D0 ) CALL errore('apply_thermostat','Negative target temperature',1)
               !
               fcp_temperature = fcp_temperature + delta_t
               !
               CALL thermalize( 0, temp_new, fcp_temperature )
               !
            ENDIF
            !

         CASE( 'berendsen', 'Berendsen' )
            !
            WRITE( UNIT = stdout, &
                FMT = '(/,5X,"FCP Soft (Berendsen) velocity rescaling")' )
            !
            CALL thermalize( nraise, temp_new, fcp_temperature )
            !
         CASE( 'andersen', 'Andersen' )
            !
            kt = fcp_temperature / ry_to_kelvin
            nat_moved = 0
            !
            IF ( randy() < 1.D0 / DBLE( nraise ) ) THEN
               !
               nat_moved = nat_moved + 1
               sigma = SQRT( kt / fcp_mass )
               !
               ! ... N.B. velocities must in a.u. units of alat and are zero
               ! ...      for fixed ions
               !
               vel = gauss_dist( 0.D0, sigma ) / alat
               !
            ENDIF
            !
            IF ( nat_moved > 0) WRITE( UNIT = stdout, &
               FMT = '(/,5X,"FCP Andersen thermostat: ",I4," collisions")' ) &
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
            tau_old = tau - vel * dt
         ENDIF
         !
      END SUBROUTINE apply_thermostat
      !
      !
      !-----------------------------------------------------------------------
      SUBROUTINE start_therm()
         !-----------------------------------------------------------------------
         !! Starting thermalization of the system.
         !
         USE symm_base,      ONLY : invsym, nsym, irt
         USE cell_base,      ONLY : alat
         USE random_numbers, ONLY : gauss_dist, set_random_seed
         !
         IMPLICIT NONE
         !
         REAL(DP) :: kt, sigma, ek, system_temp
         !
         ! ... next command prevents different MD runs to start
         ! ... with exactly the same "random" velocities
         !
         CALL set_random_seed()
         !
         kt = fcp_temperature / ry_to_kelvin
         !
         ! ... starting velocities have a Maxwell-Boltzmann distribution
         !
         sigma = SQRT( kt / fcp_mass )
         !
         ! ... N.B. velocities must in a.u. units of alat
         !
         vel = gauss_dist( 0.D0, sigma ) / alat
         !
         ! ... the velocity of fixed ions must be zero
         !
         ek = 0.5D0 * fcp_mass * vel**2
         !
         ! ... after the velocity of the center of mass has been subtracted the
         ! ... temperature is usually changed. Set again the temperature to the
         ! ... right value.
         !
         system_temp = 2.D0 * ek * alat**2 * ry_to_kelvin
         !
         CALL thermalize( 0, system_temp, fcp_temperature )
         !
      END SUBROUTINE start_therm
      !
   END SUBROUTINE fcp_verlet
   !
   !-----------------------------------------------------------------------
   SUBROUTINE thermalize( nraise, system_temp, required_temp )
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: system_temp
      !! Temperature of the system
      REAL(DP), INTENT(IN) :: required_temp
      !! Required temperature
      INTEGER, INTENT(IN) :: nraise
      !! number of steps
      !
      ! ... local variables
      !
      REAL(DP) :: aux
      !
      IF ( nraise > 0 ) THEN
         !
         ! ... Berendsen rescaling (Eq. 7.59 of Allen & Tildesley)
         ! ... the "rise time" is tau=nraise*dt so dt/tau=1/nraise
         ! ... Equivalent to traditional rescaling if nraise=1
         !
         IF ( system_temp > 0.D0 .AND. required_temp > 0.D0 ) THEN
            !
            aux = SQRT( 1.d0 + (required_temp / system_temp - 1.d0) * &
                                ( 1.D0/DBLE(nraise) ) )
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
         IF ( system_temp > 0.D0 .AND. required_temp > 0.D0 ) THEN
            !
            aux = SQRT( required_temp / system_temp )
            !
         ELSE
            !
            aux = 0.d0
            !
         ENDIF
         !
      ENDIF
      !
      vel = vel * aux
      !
   END SUBROUTINE thermalize
   !
   !------------------------------------------------------------------------
   SUBROUTINE fcp_line_minimisation( conv_fcp )
      !------------------------------------------------------------------------
      !! This routine performs one step of relaxation by using the
      !! steepest descent algorithm.
      !
      USE control_flags,  ONLY : iverbosity
      USE ener,           ONLY : ef
      USE klist,          ONLY : nelec, tot_charge
      USE ions_base,      ONLY : nat, ityp, zv
      USE cell_base,      ONLY : at, alat
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(OUT) :: conv_fcp
      !!
      !
      ! ... local variables
      !
      REAL(DP) :: force, n_tmp, max_tot_charge, capacitance, &
                  ionic_charge
      !
      LOGICAL , SAVE :: firstcall = .TRUE.
      REAL(DP), SAVE :: force0 = 0.0_DP
      REAL(DP), SAVE :: nelec0 = 0.0_DP
      !
      force = fcp_mu - ef
      !
      ! ... assumption: capacitance with vacuum gives the upper bound of 
      ! ... tot_charge difference.
      !
      capacitance = (at(1,1) * at(2,2) - at(2,1) * at(1,2)) * alat**2 &
                    / (alat * at(3,3) / 2._DP ) / fpi
      !
      max_tot_charge = ABS( capacitance * force / e2 )
      !
      IF ( firstcall .OR. ABS( force0 - force ) < 1.0D-20 ) THEN
         firstcall = .FALSE.
         nelec0 = nelec
         force0 = force
         nelec = nelec + fcp_relax_step * force
      ELSE
         n_tmp = nelec
         nelec = ( nelec * force0 - nelec0 * force ) / ( force0 - force )
         nelec0 = n_tmp
         force0 = force
      ENDIF
      !
      ionic_charge = SUM( zv(ityp(1:nat)) )
      !
      IF ( iverbosity > 1 ) THEN
         WRITE( stdout,'(/,5X,"Upper bound for tot_charge:",F12.6)') &
                 max_tot_charge
         WRITE( stdout,'(5X,"Original:",F12.6," Expected:",F12.6)') &
                 ionic_charge - nelec0, ionic_charge - nelec
      ENDIF
      !
      IF ( nelec-nelec0 < -max_tot_charge ) nelec = nelec0 - max_tot_charge
      IF ( nelec-nelec0 >  max_tot_charge ) nelec = nelec0 + max_tot_charge
      tot_charge = ionic_charge - nelec
      IF ( iverbosity > 1 ) THEN
         WRITE( stdout,'(5X,"Next tot_charge:",F12.6)') tot_charge
      ENDIF
      !
      conv_fcp = .FALSE.
      !
      IF( ABS( force ) < fcp_relax_crit ) THEN
         conv_fcp = .TRUE.
         nelec = nelec0
         tot_charge = ionic_charge - nelec
      ENDIF
      !
      WRITE( stdout, FMT = 9001 ) force
      !
9001  FORMAT(/,5X,'FCP Optimisation: Force acting on FCP =',F12.6,' Ry',/)
      !
   END SUBROUTINE fcp_line_minimisation
   !
   !
   !------------------------------------------------------------------------
   SUBROUTINE fcp_mdiis_update( conv_fcp )
      !------------------------------------------------------------------------
      !! This routine performs one step of relaxation by
      !! using the MDIIS algorithm.
      !
      USE control_flags,  ONLY : iverbosity
      USE ener,           ONLY : ef
      USE klist,          ONLY : nelec, tot_charge
      USE ions_base,      ONLY : nat, ityp, zv
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(OUT) :: conv_fcp
      !! .TRUE. if converged
      !
      ! ... local variables
      !
      REAL(DP) :: ionic_charge
      REAL(DP) :: nelec0
      REAL(DP) :: nelec1(1)
      REAL(DP) :: force
      REAL(DP) :: force1(1)
      !
      IF ( .NOT. init_mdiis ) THEN
         init_mdiis = .TRUE.
         CALL allocate_mdiis( mdiist, fcp_mdiis_size, 1, fcp_mdiis_step, 1 )
      ENDIF
      !
      nelec0    = nelec
      nelec1(1) = nelec
      force     = fcp_mu - ef
      force1(1) = force
      !
      CALL update_by_mdiis( mdiist, nelec1, force1 )
      !
      nelec = nelec1(1)
      !
      ionic_charge = SUM( zv(ityp(1:nat)) )
      !
      IF ( iverbosity > 1 ) THEN
         write( stdout,'(5X,"Original:",F12.6," Expected:",F12.6)') &
                 ionic_charge - nelec0, ionic_charge - nelec
      ENDIF
      !
      tot_charge = ionic_charge - nelec
      IF ( iverbosity > 1 ) THEN
         write( stdout,'(5X,"Next tot_charge:",F12.6)') tot_charge
      ENDIF
      !
      conv_fcp = .FALSE.
      !
      IF( ABS( force ) < fcp_relax_crit ) THEN
         !
         conv_fcp = .TRUE.
         nelec = nelec0
         tot_charge = ionic_charge - nelec
         !
      ENDIF
      !
      WRITE( stdout, FMT = 9002 ) force
      !
9002  FORMAT(/,5X,'FCP Optimisation: Force acting on FCP =',F12.6,' Ry',/)
      !
   END SUBROUTINE fcp_mdiis_update
   !
   !------------------------------------------------------------------------
   SUBROUTINE fcp_mdiis_end()
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      IF ( init_mdiis ) THEN
         CALL deallocate_mdiis( mdiist )
      ENDIF
      !
   END SUBROUTINE fcp_mdiis_end
   !
   !
   !------------------------------------------------------------------------
   SUBROUTINE fcp_summary()
      !------------------------------------------------------------------------
      !
      USE io_global,        ONLY : stdout, ionode
      USE constants,        ONLY : rytoev, BOHR_RADIUS_ANGS
      USE klist,            ONLY : tot_charge
      USE fcp_variables,    ONLY : lfcpopt, lfcpdyn
      !
      IMPLICIT NONE
      !
      IF ( .NOT. ionode ) RETURN
      !
      IF ( lfcpopt ) THEN
         WRITE( UNIT = stdout, FMT  = '(5x,"-->FCP optimiser activated<--")' )
         WRITE( UNIT = stdout, FMT = 9056 ) fcp_mu*rytoev, fcp_mu
         WRITE( UNIT = stdout, FMT = 9057 ) tot_charge
         IF ( TRIM(fcp_relax) == 'lm' ) THEN
            WRITE( UNIT = stdout, FMT = 9058 )
            WRITE( UNIT = stdout, FMT = 9059 ) fcp_relax_step
         ELSEIF ( TRIM(fcp_relax) == 'mdiis' ) THEN
            WRITE( UNIT = stdout, FMT = 9060 )
            WRITE( UNIT = stdout, FMT = 9061 ) fcp_mdiis_size
            WRITE( UNIT = stdout, FMT = 9062 ) fcp_mdiis_step
         ENDIF
         WRITE( UNIT = stdout, FMT = 9063 ) fcp_relax_crit*rytoev, &
                                            fcp_relax_crit
      ELSEIF ( lfcpdyn ) THEN
         WRITE( UNIT = stdout, FMT  = '(5x,"-->FCP optimiser activated<--")' )
         WRITE( UNIT = stdout, FMT = 9064 ) fcp_mu*rytoev, fcp_mu
         WRITE( UNIT = stdout, FMT = 9065 ) tot_charge
      ENDIF
9056  FORMAT( '     Target Fermi energy              = ', F9.4,' eV' &
             /'                                      = ', F9.4,' Ry')
9057  FORMAT( '     Initial tot_charge               = ', F9.6)
9058  FORMAT( '     FCP relax algorism               =  Line-Minimisation')
9059  FORMAT( '     FCP relax step                   = ', F9.2)
9060  FORMAT( '     FCP relax algorism               =  MDIIS')
9061  FORMAT( '     FCP MDIIS-size                   = ', I6)
9062  FORMAT( '     FCP MDIIS-step                   = ', F9.2)
9063  FORMAT( '     FCP force convergence threshold  = ', 1PE9.1,' V' &
             /'                                      = ', 1PE9.1,' Ry')
9064  FORMAT( '     FCP temperature                  = ', F9.6)
9065  FORMAT( '     FCP mass                         = ', F9.6)
      !
   END SUBROUTINE fcp_summary
   !
   !
END MODULE fcp
