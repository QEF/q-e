!
! Copyright (C) 2002-2005 Quantum ESPRESSO group
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
    SUBROUTINE add_domain_potential()
      !------------------------------------------------------------------------
      !
      ! ... a repulsive potential is added to confine the collective variables
      ! ... within the appropriate domain (used to avoid singularities):
      !
      ! ... V(s) = a*( sigma / s )^12
      !
      ! ... where a is the amplitude of the gaussians used for meta-dynamics
      !
      USE constraints_module, ONLY : constr_target, constr_type, dmax
      USE metadyn_vars,       ONLY : ncolvar, fe_grad, g_amplitude
      !
      IMPLICIT NONE
      !
      INTEGER  :: i
      REAL(DP) :: a, inv_s
      !
      REAL(DP), PARAMETER :: coord_sigma = 0.050_DP
      REAL(DP), PARAMETER :: stfac_sigma = 0.005_DP
      !
      !
      a = 12.0_DP*g_amplitude
      !
      DO i = 1, ncolvar
         !
         SELECT CASE( constr_type(i) )
         CASE( 1, 2 )
            !
            ! ... coordination must always be larger than a minimum threshold
            !
            inv_s = 1.0_DP / constr_target(i)
            !
            fe_grad(i) = fe_grad(i) - a*inv_s*( coord_sigma*inv_s )**11
            !
         CASE( 6 )
            !
            ! ... the square modulus of the structure factor is never negative
            ! ... or larger than one
            !
            inv_s = 1.0_DP / constr_target(i)
            !
            fe_grad(i) = fe_grad(i) - a*inv_s*( stfac_sigma*inv_s )**11
            !
            inv_s = 1.0_DP / ( 1.0_DP - constr_target(i) )
            !
            fe_grad(i) = fe_grad(i) - a*inv_s*( stfac_sigma*inv_s )**11
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
               new_target(i) = 2.0_DP*dmax - new_target(i)
            !
         CASE( 4, 5 )
            !
            ! ... the cosine of the angle (planar or torsional) must be
            ! ... within -1 and 1
            !
            IF ( new_target(i) > +1.0_DP ) new_target(i) = +2.0_DP - new_target(i)
            IF ( new_target(i) < -1.0_DP ) new_target(i) = -2.0_DP - new_target(i)
            !
         CASE( 6 )
            !
            ! ... the square modulus of the structure factor is never 
            ! ... negative or larger than one
            !
            new_target(i) = ABS( new_target(i) )
            !
            IF ( new_target(i) > 1.0_DP ) new_target(i) = 2.0_DP - new_target(i)
            !
         CASE( 7 )
            !
            ! ... the spherical average of the structure factor must be within
            ! ... -1 and 1
            !
            IF ( new_target(i) > +1.0_DP ) new_target(i) = +2.0_DP - new_target(i)
            IF ( new_target(i) < -1.0_DP ) new_target(i) = -2.0_DP - new_target(i)
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
      USE constraints_module, ONLY : constr_target
      !
      !
      IF ( to_new_target ) &
         constr_target(1:ncolvar) = constr_target(1:ncolvar) + to_target(:)
      !
      RETURN
      !
    END SUBROUTINE set_target
    !
    !------------------------------------------------------------------------
    SUBROUTINE mean_force( step, etot, energy_units )
      !------------------------------------------------------------------------
      !
      USE io_global,          ONLY : stdout
      USE metadyn_vars,       ONLY : dfe_acc, etot_av, ncolvar, eq_nstep
      USE constraints_module, ONLY : lagrange
      !
      INTEGER,  INTENT(IN) :: step
      REAL(DP), INTENT(IN) :: etot, energy_units
      CHARACTER(LEN=80)    :: meanfor_fmt
      !
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      !
      IF ( step <= eq_nstep ) RETURN
      !
      etot_av = etot_av + etot
      !
      dfe_acc(:) = dfe_acc(:) - lagrange(1:ncolvar)
      !
      meanfor_fmt = '(/,5X,"MEAN-FORCE ESTIMATE ",' // &
                  & TRIM( int_to_char( ncolvar ) ) // '(X,F10.6),/)'
      !
      WRITE( stdout, meanfor_fmt ) &
          dfe_acc(:) / DBLE( step - eq_nstep ) / energy_units
      !
      RETURN
      !
    END SUBROUTINE mean_force
    !
END MODULE metadyn_base
