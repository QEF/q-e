!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE metadyn_base
  !----------------------------------------------------------------------------
  !
  ! ... this module contains the core methods used to implement meta-dynamics
  !
  ! ... meta-dynamics is implemented following these two references:
  !
  ! ... 1) A. Laio and M. Parrinello; PNAS 99, 12562 (2002);
  ! ... 2) C. Micheletti, A. Laio, and M Parrinello; PRL 92, 17061 (2004).
  !
  ! ... code written by Carlo Sbraccia (2005)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE metadyn_init( progname, tau )
      !------------------------------------------------------------------------
      !
      USE kinds,              ONLY : DP
      USE input_parameters,   ONLY : restart_mode
      USE constraints_module, ONLY : target
      USE control_flags,      ONLY : nstep, ndr
      USE constants,          ONLY : bohr_radius_angs
      USE cell_base,          ONLY : at, alat
      USE metadyn_vars,       ONLY : ncolvar, g_amplitude, fe_step, &
                                     max_metadyn_iter, metadyn_fmt, &
                                     first_metadyn_iter, gaussian_pos
      USE metadyn_io,         ONLY : read_metadyn_restart
      USE io_files,           ONLY : tmp_dir, scradir, prefix, iunaxsf, &
                                     iunmeta, delete_if_present
      USE io_global,          ONLY : stdout, ionode
      USE mp,                 ONLY : mp_bcast
      USE xml_io_base,        ONLY : restart_dir
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)    :: progname
      REAL(DP),         INTENT(INOUT) :: tau(:,:)
      !
      CHARACTER(LEN=256) :: dirname
      CHARACTER(LEN=4)   :: c_ncolvar
      CHARACTER(LEN=16)  :: fe_step_fmt
      !
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      !
      c_ncolvar  = int_to_char( ncolvar )
      !
      metadyn_fmt = '(I5,' // TRIM( c_ncolvar ) // '(2X,F10.5),2X,F14.8,' // &
                  & TRIM( c_ncolvar ) // '(2X,F10.5),' // &
                  & TRIM( c_ncolvar ) // '(2X,F10.7))'
      !
      IF ( nstep < 1 ) CALL errore( 'metadyn_init', 'nstep < 1', 1 )
      !
      max_metadyn_iter = nstep
      !
      IF ( restart_mode == 'from_scratch' ) THEN
         !
         IF ( ionode ) THEN
            !
            OPEN( UNIT = iunaxsf, &
                  FILE = TRIM( prefix ) // ".axsf", STATUS = 'UNKNOWN' )
            !
            WRITE( UNIT = iunaxsf, &
                   FMT = '(" ANIMSTEPS ",I5)' ) max_metadyn_iter
            !
            WRITE( UNIT = iunaxsf, FMT = '(" CRYSTAL ")' )
            WRITE( UNIT = iunaxsf, FMT = '(" PRIMVEC ")' )
            WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
                at(1,1) * alat * bohr_radius_angs, &
                at(2,1) * alat * bohr_radius_angs, &
                at(3,1) * alat * bohr_radius_angs
            WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
                at(1,2) * alat * bohr_radius_angs, &
                at(2,2) * alat * bohr_radius_angs, &
                at(3,2) * alat * bohr_radius_angs
            WRITE( UNIT = iunaxsf, FMT = '(3F14.10)' ) &
                at(1,3) * alat * bohr_radius_angs, &
                at(2,3) * alat * bohr_radius_angs, &
                at(3,3) * alat * bohr_radius_angs
            !
         END IF
         !
         CALL delete_if_present( TRIM( prefix ) // '.metadyn' )
         !
         IF ( ionode ) THEN
            !
            OPEN( UNIT = iunmeta, &
                  FILE = TRIM( prefix ) // '.metadyn', STATUS = 'NEW' )
            !
            WRITE( iunmeta, '(2(2X,I5))' ) ncolvar, max_metadyn_iter
            WRITE( iunmeta, '(2(2X,F12.8))' ) g_amplitude
            !
            fe_step_fmt = '(' // TRIM( c_ncolvar ) // '(2X,F12.8))'
            !
            WRITE( iunmeta, fe_step_fmt ) fe_step(:)
            !
         END IF
         !
         first_metadyn_iter = 0
         !
      ELSE
         !
         ! ... restarting from file
         !
         IF ( progname == 'PW' ) THEN
            !
            dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
            !
         ELSE IF ( progname == 'CP' ) THEN
            !
            dirname = restart_dir( scradir, ndr )
            !
         ELSE
            !
            CALL errore( 'metadyn_init', &
                         'wrong calling program: ' // TRIM( progname ), 1 )
            !
         END IF
         !
         CALL read_metadyn_restart( dirname, tau, alat )
         !
         IF ( ionode ) THEN
            !
            OPEN( UNIT = iunaxsf, FILE = TRIM( prefix ) // ".axsf", &
                  STATUS = 'UNKNOWN', ACTION = 'WRITE', POSITION = 'APPEND' )
            OPEN( UNIT = iunmeta, FILE = TRIM( prefix ) // '.metadyn', &
                  STATUS = 'UNKNOWN', ACTION = 'WRITE', POSITION = 'APPEND' )
            !
         END IF
         !
      END IF
      !
      IF ( first_metadyn_iter == max_metadyn_iter ) THEN
         !
         WRITE( stdout, '(/,5X,"Simulation already completed",/)' )
         !
         CLOSE( UNIT = iunmeta, STATUS = 'KEEP' )
         !
         CALL stop_run( .FALSE. )
         !
      END IF
      !
      gaussian_pos(:) = target(1:ncolvar)
      !
      RETURN
      !
    END SUBROUTINE metadyn_init
    !
    !------------------------------------------------------------------------
    SUBROUTINE add_gaussians( iter )
      !------------------------------------------------------------------------
      !
      USE metadyn_vars, ONLY : ncolvar, metadyn_history, fe_grad, fe_step, &
                               dfe_acc, g_amplitude
      USE basic_algebra_routines
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: iter
      !
      INTEGER               :: i
      REAL(DP), ALLOCATABLE :: delta(:)
      !
      ! ... history dependent term
      !
      IF ( iter == 1 ) RETURN
      !
      ALLOCATE( delta( ncolvar ) )
      !
      dfe_acc = 0.D0
      !
      DO i = 1, iter - 1
         !
         delta = metadyn_history(:,iter) - metadyn_history(:,i)
         !
         dfe_acc(:) = dfe_acc(:) + delta(:) / fe_step(:)**2 * &
                      EXP( - SUM( delta(:)**2 / ( 2.D0 * fe_step(:)**2 ) ) )
         !
      END DO
      !
      fe_grad(:) = fe_grad(:) - g_amplitude * dfe_acc(:)
      !
      DEALLOCATE( delta )
      !
      RETURN
      !
    END SUBROUTINE add_gaussians
    !
    !------------------------------------------------------------------------
    SUBROUTINE add_domain_potential()
      !------------------------------------------------------------------------
      !
      ! ... a repulsive potential is added to confine the collective variables
      ! ... within the appropriate domain :
      !
      ! ... V(s) = A * ( sigma / s )^12
      !
      ! ... where A is the amplitude of the gaussians used for meta-dynamics
      !
      USE constraints_module, ONLY : target, constr_type, dmax
      USE metadyn_vars,       ONLY : ncolvar, fe_grad, g_amplitude
      !
      IMPLICIT NONE
      !
      INTEGER  :: i
      REAL(DP) :: A, inv_s
      !
      REAL(DP), PARAMETER :: coord_sigma = 0.050D0
      REAL(DP), PARAMETER :: dist_sigma  = 0.050D0
      REAL(DP), PARAMETER :: stfac_sigma = 0.005D0
      !
      !
      A = 12.D0 * g_amplitude
      !
      DO i = 1, ncolvar
         !
         SELECT CASE( constr_type(i) )
         CASE( 1, 2 )
            !
            ! ... coordination must always be larger than a minimum threshold
            !
            inv_s = 1.D0 / target(i)
            !
            fe_grad(i) = fe_grad(i) - A * inv_s * ( coord_sigma * inv_s )**11
            !
         CASE( 3 )
            !
            ! ... a distance can never be larger than dmax ( check file
            ! ... constraints_module.f90 for its definition )
            !
            inv_s = 1.D0 / ( dmax - target(i) )
            !
            fe_grad(i) = fe_grad(i) - A * inv_s * ( dist_sigma * inv_s )**11
            !
         CASE( 6 )
            !
            ! ... the square modulus of the structure factor is never negative
            ! ... or larger than one
            !
            inv_s = 1.D0 / target(i)
            !
            fe_grad(i) = fe_grad(i) - A * inv_s * ( stfac_sigma * inv_s )**11
            !
            inv_s = 1.D0 / ( 1.D0 - target(i) )
            !
            fe_grad(i) = fe_grad(i) - A * inv_s * ( stfac_sigma * inv_s )**11
            !
         END SELECT
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE add_domain_potential
    !
    !------------------------------------------------------------------------
    SUBROUTINE evolve_collective_vars( norm_fe_grad )
      !------------------------------------------------------------------------
      !
      ! ... the collective variables are evolved taking care of the
      ! ... additional constraints imposed by the domain definition
      !
      USE constants,          ONLY : eps32
      USE constraints_module, ONLY : target
      USE metadyn_vars,       ONLY : ncolvar, fe_grad, fe_step, new_target, &
                                     to_target, shake_nstep, gaussian_pos, &
                                     g_amplitude
      USE random_numbers,     ONLY : rndm
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: norm_fe_grad
      !
      INTEGER  :: i
      REAL(DP) :: step
      !
      !
      IF ( norm_fe_grad < eps32 ) &
         CALL errore( 'evolve_collective_vars', 'norm( fe_grad ) = 0', 1 )
      !
      IF ( g_amplitude > 0.D0 ) fe_grad(:) = fe_grad(:) / norm_fe_grad
      !
      DO i = 1, ncolvar
         !
         gaussian_pos(i) = target(i) - fe_step(i) * fe_grad(i)
         !
         step = ( 1.D0 + 0.5D0*rndm() ) * fe_step(i)
         !
         new_target(i) = target(i) - step * fe_grad(i)
         !
      END DO
      !
      CALL impose_domain_constraints()
      !
      to_target(:) = ( new_target(:) - target(1:ncolvar) ) / DBLE( shake_nstep )
      !
      RETURN
      !
    END SUBROUTINE evolve_collective_vars
    !
    !------------------------------------------------------------------------
    SUBROUTINE impose_domain_constraints()
      !------------------------------------------------------------------------
      !
      USE constraints_module, ONLY : constr_type, dmax
      USE metadyn_vars,       ONLY : ncolvar, new_target
      !
      IMPLICIT NONE
      !
      INTEGER :: i
      !
      !
      DO i = 1, ncolvar
         !
         SELECT CASE( constr_type(i) )
         CASE( 1, 2 )
            !
            ! ... coordination must always be larger than zero
            !
            new_target(i) = ABS( new_target(i) )
            !
         CASE( 3 )
            !
            ! ... a distance can never be larger than dmax ( check file
            ! ... constraints_module.f90 for its definition )
            !
            IF ( new_target(i) > dmax ) &
               new_target(i) = 2.D0 * dmax - new_target(i)
            !
         CASE( 6 )
            !
            ! ... the square modulus of the structure factor is never 
            ! ... negative or larger than one
            !
            new_target(i) = ABS( new_target(i) )
            !
            IF ( new_target(i) > 1.D0 ) new_target(i) = 2.D0 - new_target(i)
            !
         CASE( 7 )
            !
            ! ... the spherical average of the structure factor must be within
            ! ... -1 and 1
            !
            IF ( new_target(i) > +1.D0 ) new_target(i) = +2.D0 - new_target(i)
            IF ( new_target(i) < -1.D0 ) new_target(i) = -2.D0 - new_target(i)
            !
         END SELECT
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE impose_domain_constraints
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_target()
      !------------------------------------------------------------------------
      !
      USE metadyn_vars,       ONLY : ncolvar, to_target, to_new_target
      USE constraints_module, ONLY : target
      !
      !
      IF ( to_new_target ) &
         target(1:ncolvar) = target(1:ncolvar) + to_target(:)
      !
      RETURN
      !
    END SUBROUTINE set_target
    !
END MODULE metadyn_base
