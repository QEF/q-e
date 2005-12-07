!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#define __REMOVE_CONSTRAINT_FORCE
!#define __DEBUG_CONSTRAINTS
!#define __USE_PBC
!
!----------------------------------------------------------------------------
MODULE constraints_module
  !----------------------------------------------------------------------------
  ! 
  ! ... variables and methods for constraint Molecular Dynamics and
  ! ... constrained ionic relaxations (the SHAKE algorithm based on 
  ! ... lagrange multipliers) are defined here.
  !
  ! ... written by Carlo Sbraccia ( 24/02/2004 )
  !
  ! ... references :
  !
  ! ... 1) M. P. Allen and D. J. Tildesley, Computer Simulations of Liquids,
  ! ...    Clarendon Press - Oxford (1986)
  !
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps32, tpi
  USE io_global, ONLY : stdout
  !
  USE basic_algebra_routines, ONLY : norm
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  ! ... public methods
  !
  PUBLIC :: init_constraint,     &
            check_constraint,    &
            remove_constr_force, &
            deallocate_constraint
  !
  ! ... public variables (assigned in the CONSTRAINTS input card)
  !
  PUBLIC :: nconstr,     &
            constr_tol,  &
            constr_type, &
            constr,      &
            lagrange,    &
            target,      &
            dmax
  !
  ! ... global variables
  !
  INTEGER               :: nconstr 
  REAL(DP)              :: constr_tol
  INTEGER,  ALLOCATABLE :: constr_type(:)
  REAL(DP), ALLOCATABLE :: constr(:,:)
  REAL(DP), ALLOCATABLE :: target(:)
  REAL(DP), ALLOCATABLE :: lagrange(:)
  REAL(DP)              :: dmax
  !
  !
  CONTAINS
     !
     ! ... public methods
     !
     !-----------------------------------------------------------------------
     SUBROUTINE init_constraint( nat, tau, tau_units, ityp )
       !-----------------------------------------------------------------------
       !
       USE input_parameters, ONLY : nconstr_inp, constr_tol_inp, &
                                    constr_type_inp, constr_inp, &
                                    constr_target, constr_target_set
       USE parser,           ONLY : int_to_char
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN) :: nat
       REAL(DP), INTENT(IN) :: tau(3,nat)
       REAL(DP), INTENT(IN) :: tau_units
       INTEGER,  INTENT(IN) :: ityp(nat)
       !
       INTEGER  :: ia, ia1, ia2, ia3, n_type_coord1
       REAL(DP) :: r21(3), r23(3)
       REAL(DP) :: k, r_c
       INTEGER  :: type_coord1, type_coord2
       REAL(DP) :: dtau(3), norm_dtau
       !
       !
       nconstr    = nconstr_inp
       constr_tol = constr_tol_inp
       !
       ALLOCATE( lagrange(    nconstr ) )
       ALLOCATE( target(      nconstr ) )
       ALLOCATE( constr_type( nconstr ) )
       !
       ALLOCATE( constr( 4, nconstr ) )
       !
       constr_type(:) = constr_type_inp(1:nconstr)
       constr(:,:)    = constr_inp(:,1:nconstr)
       !
       ! ... set the largest possible distance among two atoms within
       ! ... the supercell
       !
       IF ( ANY( constr_type(:) == 3 ) ) CALL compute_dmax()
       !
       ! ... target value of the constrain ( in bohr )
       !
       DO ia = 1, nconstr
          !
          SELECT CASE ( constr_type(ia) )
          CASE( 1 )
             !
             ! ... constraint on global coordination-number, i.e. the average 
             ! ... number of atoms of type B surrounding the atoms of type A
             !
             IF ( constr_target_set(ia) ) THEN
                !
                target(ia) = constr_target(ia)
                !
                CYCLE
                !
             END IF
             !
             type_coord1 = ANINT( constr(1,ia) )
             type_coord2 = ANINT( constr(2,ia) )
             !
             r_c = constr(3,ia)
             k   = constr(4,ia)
             !
             target(ia) = 0.D0
             !
             n_type_coord1 = 0
             !
             DO ia1 = 1, nat
                !
                IF ( ityp(ia1) /= type_coord1 ) CYCLE
                !
                DO ia2 = 1, nat
                   !
                   IF ( ia2 == ia1 ) CYCLE
                   !
                   IF ( ityp(ia2) /= type_coord2 ) CYCLE
                   !
                   dtau = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
                   !
                   norm_dtau = norm( dtau )
                   !
                   target(ia) = target(ia) + &
                                1.D0 / ( EXP( k * ( norm_dtau - r_c ) ) + 1.D0 )
                   !
                END DO
                !
                n_type_coord1 = n_type_coord1 + 1
                !
             END DO
             !
             target(ia) = target(ia) / DBLE( n_type_coord1 )
             !
          CASE( 2 )
             !
             ! ... constraint on local coordination-number, i.e. the average 
             ! ... number of atoms of type A surrounding a specific atom
             !
             IF ( constr_target_set(ia) ) THEN
                !
                target(ia) = constr_target(ia)
                !
                CYCLE
                !
             END IF
             !
             ia1         = ANINT( constr(1,ia) )
             type_coord1 = ANINT( constr(2,ia) )
             !
             r_c = constr(3,ia)
             k   = constr(4,ia)
             !
             target(ia) = 0.D0
             !
             DO ia2 = 1, nat
                !
                IF ( ia2 == ia1 ) CYCLE
                !
                IF ( ityp(ia2) /= type_coord1 ) CYCLE
                !
                dtau = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
                !
                norm_dtau = norm( dtau )
                !
                target(ia) = target(ia) + &
                             1.D0 / ( EXP( k * ( norm_dtau - r_c ) ) + 1.D0 )
                !
             END DO
             !
          CASE( 3 )
             !
             ! ... constraint on distance
             !
             IF ( constr_target_set(ia) ) THEN
                !
                target(ia) = constr_target(ia)
                !
             ELSE
                !
                ia1 = ANINT( constr(1,ia) )
                ia2 = ANINT( constr(2,ia) )
                !
                dtau = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
                !
                target(ia) = norm( dtau )
                !
             END IF
             !
             IF ( target(ia) > dmax ) &
                CALL errore( 'init_constraint', 'the target for constraint ' //&
                           & TRIM( int_to_char( ia ) ) // ' is larger than ' //&
                           & 'the largest possible value', 1 )
             !
          CASE( 4 )
             !
             ! ... constraint on planar angle
             !
             IF ( constr_target_set(ia) ) THEN
                !
                ! ... in the input target for the angle (in degrees) is 
                ! ... converted to the cosine of the angle
                !
                target(ia) = COS( constr_target(ia) * tpi / 360.D0 )
                !
                CYCLE
                !
             END IF
             !
             ia1 = ANINT( constr(1,ia) )
             ia2 = ANINT( constr(2,ia) )
             ia3 = ANINT( constr(3,ia) )
             !
             r21 = pbc( ( tau(:,ia2) - tau(:,ia1) ) * tau_units )
             r23 = pbc( ( tau(:,ia2) - tau(:,ia3) ) * tau_units )
             !
             r21 = r21 / norm( r21 )
             r23 = r23 / norm( r23 )
             !
             target(ia) = DOT_PRODUCT( r21, r23 )
             !
          CASE DEFAULT
             !
             CALL errore( 'init_constraint', &
                          'constrain type not implemented', 1 )
             !
          END SELECT
          !
       END DO
       !
       RETURN
       !
     END SUBROUTINE init_constraint
     !
     !-----------------------------------------------------------------------
     SUBROUTINE constraint_grad( index, nat, tau, &
                                 if_pos, ityp, tau_units, g, dg )
       !-----------------------------------------------------------------------
       !
       ! ... this routine defines the constraint equation:
       !
       ! ...  g(tau,dist) = 0
       !
       ! ... where tau are the atomic positions ( in tau_units ) and dist is
       ! ... the distance of two atoms ( in this case atom 1 and atom 2 ) which 
       ! ... is, in this case, a one dimensional constrain. 
       ! ... dg is in output the value of the gradient of g and dg2 is its 
       ! ... square modulus.
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN)  :: index
       INTEGER,  INTENT(IN)  :: nat
       REAL(DP), INTENT(IN)  :: tau(:,:)
       INTEGER,  INTENT(IN)  :: if_pos(:,:)
       INTEGER,  INTENT(IN)  :: ityp(:)
       REAL(DP), INTENT(IN)  :: tau_units
       REAL(DP), INTENT(OUT) :: dg(:,:)
       REAL(DP), INTENT(OUT) :: g
       !
       ! ... local variables
       !
       INTEGER  :: ia, ia1, ia2, ia3, n_type_coord1
       REAL(DP) :: r21(3), r23(3)
       REAL(DP) :: norm_r21, norm_r23, cos123, sin123
       REAL(DP) :: k, r_c
       INTEGER  :: type_coord1, type_coord2
       REAL(DP) :: dtau(3), norm_dtau, expo
       !
       ! ... external function
       !
       REAL(DP), EXTERNAL :: DDOT
       !
       !
       dg(:,:) = 0.D0
       !
       SELECT CASE ( constr_type(index) )
       CASE( 1 )
          !
          ! ... constraint on global coordination
          !
          type_coord1 = ANINT( constr(1,index) )
          type_coord2 = ANINT( constr(2,index) )
          !
          r_c = constr(3,index)
          k   = constr(4,index)
          !
          g = 0.D0
          !
          n_type_coord1 = 0
          !
          DO ia1 = 1, nat
             !
             IF ( ityp(ia1) /= type_coord1 ) CYCLE
             !
             DO ia2 = 1, nat
                !
                IF ( ia2 == ia1 ) CYCLE
                !
                IF ( ityp(ia2) /= type_coord2 ) CYCLE
                !
                dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
                !
                norm_dtau = norm( dtau(:) )
                !
                dtau(:) = dtau(:) / norm_dtau
                !
                expo = EXP( k * ( norm_dtau - r_c ) )
                !
                g = g + 1.D0 / ( expo + 1.D0 )
                !
                dtau(:) = dtau(:) * k * expo / ( expo + 1.D0 )**2
                !
                dg(:,ia2) = dg(:,ia2) + dtau(:)
                dg(:,ia1) = dg(:,ia1) - dtau(:)
                !
             END DO
             !
             n_type_coord1 = n_type_coord1 + 1
             !
          END DO
          !
          g  = g  / DBLE( n_type_coord1 )
          dg = dg / DBLE( n_type_coord1 )
          !
          g = ( g - target(index) )
          !
       CASE( 2 )
          !
          ! ... constraint on local coordination
          !
          ia          = ANINT( constr(1,index) )
          type_coord1 = ANINT( constr(2,index) )
          !
          r_c = constr(3,index)
          k   = constr(4,index)
          !
          g = 0.D0
          !
          DO ia1 = 1, nat
             !
             IF ( ia1 == ia ) CYCLE
             !
             IF ( ityp(ia1) /= type_coord1 ) CYCLE
             !
             dtau(:) = pbc( ( tau(:,ia) - tau(:,ia1) ) * tau_units )
             !
             norm_dtau = norm( dtau(:) )
             !
             dtau(:) = dtau(:) / norm_dtau
             !
             expo = EXP( k * ( norm_dtau - r_c ) )
             !
             g = g + 1.D0 / ( expo + 1.D0 )
             !
             dtau(:) = dtau(:) * k * expo / ( expo + 1.D0 )**2
             !
             dg(:,ia1) = dg(:,ia1) + dtau(:)
             dg(:,ia)  = dg(:,ia)  - dtau(:)
             !
          END DO
          !
          g = ( g - target(index) )
          !
       CASE( 3 )
          !
          ! ... constraint on distances
          !
          ia1 = ANINT( constr(1,index) )
          ia2 = ANINT( constr(2,index) )
          !
          dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
          !
          norm_dtau = norm( dtau(:) )
          !
          g = ( norm_dtau - target(index) )
          !
          dg(:,ia1) = dtau(:) / norm_dtau
          !
          dg(:,ia2) = - dg(:,ia1)
          !
       CASE( 4 )
          !
          ! ... constraint on planar angles
          !
          ia1 = ANINT( constr(1,index) )
          ia2 = ANINT( constr(2,index) )
          ia3 = ANINT( constr(3,index) )
          !
          r21 = pbc( ( tau(:,ia2) - tau(:,ia1) ) * tau_units )
          r23 = pbc( ( tau(:,ia2) - tau(:,ia3) ) * tau_units )
          !
          norm_r21 = norm( r21 )
          norm_r23 = norm( r23 )
          !
          r21 = r21 / norm_r21
          r23 = r23 / norm_r23
          !
          cos123 = DOT_PRODUCT( r21, r23 )
          sin123 = SQRT( 1.D0 - cos123**2 )
          !
          g = ( cos123 - target(index) )
          !
          dg(:,ia1) = ( cos123 * r21 - r23 ) / ( sin123 * norm_r21 )
          dg(:,ia3) = ( cos123 * r23 - r21 ) / ( sin123 * norm_r23 )
          dg(:,ia2) = - dg(:,ia1) - dg(:,ia3)
          !
       CASE DEFAULT
          !
          CALL errore( 'constraint_grad', &
                       'constrain type not implemented', 1 )
          !
       END SELECT
       !
       dg = dg * DBLE( if_pos )
       !
       RETURN
       !
     END SUBROUTINE constraint_grad
     !
     !-----------------------------------------------------------------------
     SUBROUTINE check_constraint( nat, taup, tau0, &
                                  force, if_pos, ityp, tau_units, dt, massconv )
       !-----------------------------------------------------------------------
       !
       ! ... update tau so that the constraint equation g=0 is satisfied,
       ! ... use the recursion formula:
       !
       ! ...                       g(taup)
       ! ... taup = taup - ----------------------- * dg(tau0)
       ! ...               M^-1<dg(taup)|dg(tau0)>
       !
       ! ... in normal cases the constraint equation should be always 
       ! ... satisfied at the very first iteration.
       !
       USE ions_base, ONLY : amass
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN)    :: nat
       REAL(DP), INTENT(INOUT) :: taup(3,nat)
       REAL(DP), INTENT(IN)    :: tau0(3,nat)
       INTEGER,  INTENT(IN)    :: if_pos(3,nat)
       REAL(DP), INTENT(INOUT) :: force(3,nat)
       INTEGER,  INTENT(IN)    :: ityp(nat)
       REAL(DP), INTENT(IN)    :: tau_units
       REAL(DP), INTENT(IN)    :: dt
       REAL(DP), INTENT(IN)    :: massconv
       !
       INTEGER               :: na, i, index
       REAL(DP), ALLOCATABLE :: gp(:), dgp(:,:), dg0(:,:,:)
       REAL(DP)              :: g0
       REAL(DP)              :: lambda, fac, invdtsq
       LOGICAL               :: ltest(nconstr), global_test
       INTEGER, PARAMETER    :: maxiter = 100
       !
       REAL(DP), EXTERNAL :: DDOT
       !
       !
       ALLOCATE( dgp( 3, nat ) )
       ALLOCATE( dg0( 3, nat, nconstr ) )
       ALLOCATE( gp( nconstr ) )
       !
       invdtsq  = 1.D0 / dt**2
       !
       DO index = 1, nconstr
          !
          CALL constraint_grad( index, nat, tau0, &
                                if_pos, ityp, tau_units, g0, dg0(:,:,index) )
          !
       END DO
       !
       outer_loop: DO i = 1, maxiter
          !
          inner_loop: DO index = 1, nconstr
             !
             ltest(index) = .FALSE.
             !
             CALL constraint_grad( index, nat, taup, &
                                   if_pos, ityp, tau_units, gp(index), dgp )
             !
             ! ... check if gp = 0
             !
#if defined (__DEBUG_CONSTRAINTS)
             WRITE( stdout, * ) i, index, ABS( gp(index) )
#endif
             !             
             IF ( ABS( gp(index) ) < constr_tol ) THEN
                !
                ltest(index) = .TRUE.
                !
                CYCLE inner_loop
                !
             END IF
             !
             ! ... if  gp <> 0  find new taup check again 
             ! ... ( gp is in bohr and taup in tau_units )
             !
             DO na = 1, nat
                !
                dgp(:,na) = dgp(:,na) / ( amass( ityp(na) ) * massconv )
                !
             END DO
             !
             lambda = gp(index) / DDOT( 3 * nat, dgp, 1, dg0(:,:,index), 1 )
             !
             DO na = 1, nat
                !
                fac = amass( ityp(na) ) * massconv * tau_units
                !
                taup(:,na) = taup(:,na) - lambda * dg0(:,na,index) / fac
                !
             END DO
             !
             force(:,:) = force(:,:) - lambda * dg0(:,:,index) * invdtsq
             !
             lagrange(index) = lagrange(index) + lambda * invdtsq
             !
          END DO inner_loop
          !
          global_test = ALL( ltest(:) )
          !
          ! ... all constraints are satisfied
          !
          IF ( global_test ) EXIT outer_loop
          !
       END DO outer_loop
       !
       IF ( .NOT. global_test ) THEN
          !
          ! ... error messages
          !
          WRITE( stdout, '(/,5X,"Number of step(s): ",I3)') MIN( i, maxiter )
          WRITE( stdout, '(/,5X,"target convergence: ")' )
          !
          DO i = 1, nconstr
             !
             WRITE( stdout, '(5X,"target = ",I3,2X,L1,2(2X,F16.10))' ) &
                 i, ltest(i), ABS( gp(i) ), target(i)
             !
          END DO
          !
          CALL errore( 'check_constrain', &
                       'on some constraint g = 0 is not satisfied', 1 )
          !
       END IF
       !
       DEALLOCATE( dgp )
       DEALLOCATE( dg0 )
       DEALLOCATE( gp )
       !
       RETURN
       !
     END SUBROUTINE check_constraint
     !
     !-----------------------------------------------------------------------
     SUBROUTINE remove_constr_force( nat, tau, &
                                     if_pos, ityp, tau_units, force )
       !-----------------------------------------------------------------------
       !
       ! ... the component of the force that is orthogonal to the
       ! ... ipersurface defined by the constraint equations is removed
       ! ... and the corresponding value of the lagrange multiplier computed
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN)    :: nat
       REAL(DP), INTENT(IN)    :: tau(:,:)
       INTEGER,  INTENT(IN)    :: if_pos(:,:)
       INTEGER,  INTENT(IN)    :: ityp(:)
       REAL(DP), INTENT(IN)    :: tau_units
       REAL(DP), INTENT(INOUT) :: force(:,:)
       !
       INTEGER               :: i, j
       REAL(DP)              :: g, norm_before, norm_after
       REAL(DP), ALLOCATABLE :: dg(:,:,:)
       REAL(DP), ALLOCATABLE :: dg_matrix(:,:)
       INTEGER,  ALLOCATABLE :: iwork(:)
       !
       REAL(DP), EXTERNAL :: DDOT, DNRM2
       !
       !
       lagrange(:) = 0.D0
       !
#if defined (__REMOVE_CONSTRAINT_FORCE)
       !
       norm_before = DNRM2( 3 * nat, force, 1 )
       !
       ALLOCATE( dg( 3, nat, nconstr ) )
       ALLOCATE( dg_matrix( nconstr, nconstr ) )
       ALLOCATE( iwork( nconstr ) )
       !
       DO i = 1, nconstr
          !
          CALL constraint_grad( i, nat, tau, if_pos, &
                                ityp, tau_units, g, dg(:,:,i) )
          !
       END DO
       !
       DO i = 1, nconstr
          !
          dg_matrix(i,i) = DDOT( 3 * nat, dg(:,:,i), 1, dg(:,:,i), 1 )
          !
          lagrange(i) = DDOT( 3 * nat, force, 1, dg(:,:,i), 1 )
          !
          DO j = i + 1, nconstr
             !
             dg_matrix(i,j) = DDOT( 3 * nat, dg(:,:,i), 1, dg(:,:,j), 1 )
             dg_matrix(j,i) = dg_matrix(i,j)
             !
          END DO
          !
       END DO
       !
       CALL DGESV( nconstr, 1, dg_matrix, nconstr, iwork, lagrange, nconstr, i )
       !
       IF ( i /= 0 ) &
          CALL errore( 'remove_constr_force', &
                       'error in the solution of the linear system', 1 )
       !
       DO i = 1, nconstr
          !
          force(:,:) = force(:,:) - lagrange(i) * dg(:,:,i)
          !
       END DO
       !
#if defined (__DEBUG_CONSTRAINTS)
       !
       WRITE( stdout, '(/,5X,"Intermediate forces (Ry/au):",/)')
       !
       DO i = 1, nat
          !
          WRITE( stdout, '(5X,"atom ",I3," type ",I2,3X,"force = ",3F14.8)' ) &
              i, ityp(i), force(:,i)
          !
       END DO
       !
#endif
       !
       norm_after = DNRM2( 3 * nat, force, 1 )
       !
       IF ( norm_before < norm_after ) THEN
          !
          WRITE( stdout, '(/,5X,"norm before = ",F16.10)' ) norm_before
          WRITE( stdout, '(  5X,"norm after  = ",F16.10)' ) norm_after
          !
          CALL errore( 'remove_constr_force', &
                       'norm(F) before < norm(F) after', 1 )
          !
       END IF
       !
       DEALLOCATE( dg )
       DEALLOCATE( dg_matrix )
       DEALLOCATE( iwork )
       !
#endif
       !
     END SUBROUTINE remove_constr_force
     !
     !-----------------------------------------------------------------------
     SUBROUTINE deallocate_constraint()
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       !
       IF ( ALLOCATED( lagrange ) )    DEALLOCATE( lagrange )
       IF ( ALLOCATED( constr ) )      DEALLOCATE( constr )
       IF ( ALLOCATED( constr_type ) ) DEALLOCATE( constr_type )
       IF ( ALLOCATED( target ) )      DEALLOCATE( target )       
       !
       RETURN
       !
     END SUBROUTINE deallocate_constraint
     !
     !-----------------------------------------------------------------------
     FUNCTION pbc( vect )
       !-----------------------------------------------------------------------
       !
       ! ... periodic boundary conditions ( vect is assumed to be given
       ! ... in cartesian units )
       !
       USE cell_base, ONLY : at, bg, alat
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: vect(3)
       REAL(DP)             :: pbc(3)
       !
       !
#if defined (__USE_PBC)
       !
       pbc(:) = MATMUL( vect(:), bg(:,:) ) / alat
       !
       pbc(:) = pbc(:) - ANINT( pbc(:) )
       !
       pbc(:) = MATMUL( at(:,:), pbc(:) ) * alat
       !
       !
#else
       !
       pbc(:) = vect(:)
       !
#endif
       RETURN
       !
     END FUNCTION pbc
     !
     !-----------------------------------------------------------------------
     SUBROUTINE compute_dmax()
       !-----------------------------------------------------------------------
       !
       USE cell_base, ONLY : at, alat
       !
       IMPLICIT NONE
       !
       dmax = norm( MATMUL( at(:,:), (/ 0.5D0, 0.5D0, 0.5D0 /) ) ) * alat
       !
       RETURN
       !
     END SUBROUTINE compute_dmax
     !
END MODULE constraints_module
