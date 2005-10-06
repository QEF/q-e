!
! Copyright (C) 2004-2005 PWSCF-CP90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!#define __DEBUG_CONSTRAINTS
#define __REMOVE_CONSTRAINT_FORCE
!
!----------------------------------------------------------------------------
MODULE constraints_module
  !----------------------------------------------------------------------------
  ! 
  ! ... variables and methods for constraint Molecular Dynamics and
  ! ... constraint ionic relaxations (the SHAKE algorithm based on 
  ! ... lagrange multipliers) are defined here.
  !
  ! ... Written by Carlo Sbraccia ( 24/02/2004 )
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
  PUBLIC :: init_constraint,         &
            check_constraint,        &
            remove_constraint_force, &
            deallocate_constraint
  !
  ! ... public variables (assigned in the CONSTRAINTS input card)
  !
  PUBLIC :: nconstr,     &
            constr_tol,  &
            constr_type, &
            constr,      &
            lagrange,    &
            target
  !
  ! ... global variables
  !
  INTEGER               :: nconstr 
  REAL(DP)              :: constr_tol
  INTEGER,  ALLOCATABLE :: constr_type(:)
  REAL(DP), ALLOCATABLE :: constr(:,:)
  REAL(DP), ALLOCATABLE :: target(:)
  REAL(DP), ALLOCATABLE :: lagrange(:)
  !
  !
  CONTAINS
     !
     ! ... public methods
     !
     !-----------------------------------------------------------------------
     SUBROUTINE init_constraint( nat, tau, alat, ityp )
       !-----------------------------------------------------------------------
       !
       USE input_parameters, ONLY : nconstr_inp, constr_tol_inp, &
                                    constr_type_inp, constr_inp, &
                                    constr_target, constr_target_set
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN) :: nat
       REAL(DP), INTENT(IN) :: tau(3,nat)
       REAL(DP), INTENT(IN) :: alat
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
       ! ... target value of the constrain ( in bohr )
       !
       DO ia = 1, nconstr
          !
          SELECT CASE ( constr_type(ia) )
          CASE( 1 )
             !
             ! ... constraint on global coordination number
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
                   dtau = ( tau(:,ia1) - tau(:,ia2) ) * alat
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
             ! ... constraint on local coordination number
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
                dtau = ( tau(:,ia1) - tau(:,ia2) ) * alat
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
                CYCLE
                !
             END IF
             !
             ia1 = ANINT( constr(1,ia) )
             ia2 = ANINT( constr(2,ia) )
             !
             target(ia) = norm( tau(:,ia1) - tau(:,ia2) ) * alat
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
             r21 = ( tau(:,ia2) - tau(:,ia1) ) * alat
             r23 = ( tau(:,ia2) - tau(:,ia3) ) * alat
             !
             r21 = r21 / norm( r21 )
             r23 = r23 / norm( r23 )
             !
             target(ia) = DOT_PRODUCT( r21, r23 )
             !
          CASE DEFAULT
             !
             CALL errore( 'init_constraint', &
                          'constrain type not implemented ', 1 )
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
     SUBROUTINE constraint_grad( index, nat, tau, if_pos, ityp, alat, g, dg )
       !-----------------------------------------------------------------------
       !
       ! ... this routine defines the constraint equation:
       !
       ! ...  g(tau,dist) = 0
       !
       ! ... where tau are the atomic positions ( in alat units ) and dist is
       ! ... the distance of two atoms ( in this case atom 1 and atom 2 ) which 
       ! ... is, in this case, a one dimensional constrain. 
       ! ... dg is in output the value of the gradient of g and dg2 is its 
       ! ... square modulus.
       !
       ! ... Dario Alfe 1997  and  Carlo Sbraccia 2004
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN)  :: index
       INTEGER,  INTENT(IN)  :: nat
       REAL(DP), INTENT(IN)  :: tau(:,:)
       INTEGER,  INTENT(IN)  :: if_pos(:,:)
       INTEGER,  INTENT(IN)  :: ityp(:)
       REAL(DP), INTENT(IN)  :: alat
       REAL(DP), INTENT(OUT) :: dg(:,:)
       REAL(DP), INTENT(OUT) :: g
       !
       ! ... local variables
       !
       REAL(DP) :: x1, x2, y1, y2, z1, z2
       REAL(DP) :: dist0
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
                dtau(:) = ( tau(:,ia1) - tau(:,ia2) ) * alat
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
             dtau(:) = ( tau(:,ia) - tau(:,ia1) ) * alat
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
          x1 = tau(1,ia1) * alat
          y1 = tau(2,ia1) * alat
          z1 = tau(3,ia1) * alat
          x2 = tau(1,ia2) * alat
          y2 = tau(2,ia2) * alat
          z2 = tau(3,ia2) * alat
          !
          ! ... the actual distance between the two atoms ( in bohr )
          !
          dist0 = SQRT( ( x1 - x2 )**2 + ( y1 - y2 )**2 + ( z1 - z2 )**2 )
          !
          g = ( dist0 - target(index) )
          !
          dg(1,ia1) = ( x1 - x2 ) / dist0
          dg(2,ia1) = ( y1 - y2 ) / dist0
          dg(3,ia1) = ( z1 - z2 ) / dist0
          !
          dg(1,ia2) = - dg(1,ia1)
          dg(2,ia2) = - dg(2,ia1)
          dg(3,ia2) = - dg(3,ia1)
          !
       CASE( 4 )
          !
          ! ... constraint on planar angles
          !
          ia1 = ANINT( constr(1,index) )
          ia2 = ANINT( constr(2,index) )
          ia3 = ANINT( constr(3,index) )
          !
          r21 = ( tau(:,ia2) - tau(:,ia1) ) * alat
          r23 = ( tau(:,ia2) - tau(:,ia3) ) * alat
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
          CALL errore( 'dist_constrain', &
                       'constrain type not implemented ', 1 )
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
                                  force, if_pos, ityp, alat, dt, massconv )
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
       REAL(DP), INTENT(IN)    :: alat
       REAL(DP), INTENT(IN)    :: dt
       REAL(DP), INTENT(IN)    :: massconv
       !
       INTEGER               :: na, i, index
       REAL(DP), ALLOCATABLE :: dgp(:,:), dg0(:,:)
       REAL(DP)              :: gp, g0
       REAL(DP)              :: lambda, fac, invdtsq
       LOGICAL               :: ltest(nconstr), global_test
       INTEGER, PARAMETER    :: maxiter = 100
       !
       REAL(DP), EXTERNAL :: DDOT
       !
       !
       ALLOCATE( dgp( 3, nat ) )
       ALLOCATE( dg0( 3, nat ) )
       !
       invdtsq  = 1.D0 / dt**2
       !
       outer_loop: DO i = 1, maxiter
          !
          inner_loop: DO index = 1, nconstr
             !
             ltest(index) = .FALSE.
             !
             CALL constraint_grad( index, nat, taup, &
                                   if_pos, ityp, alat, gp, dgp )
             !
             ! ... check if gp = 0
             !
#if defined (__DEBUG_CONSTRAINTS)
             WRITE( stdout, * ) i, index, ABS( gp )
#endif
             !             
             IF ( ABS( gp ) < constr_tol ) THEN
                !
                ltest(index) = .TRUE.
                !
                CYCLE inner_loop
                !
             END IF
             !
             ! ... if  gp <> 0  find new  taup check again 
             ! ... ( gp is in bohr and taup in alat units )
             !
             DO na = 1, nat
                !
                dgp(:,na) = dgp(:,na) / ( amass( ityp(na) ) * massconv )
                !
             END DO
             !
             CALL constraint_grad( index, nat, tau0, &
                                   if_pos, ityp, alat, g0, dg0 )
             !
             lambda = gp / DDOT( 3 * nat, dgp, 1, dg0, 1 )
             !
             DO na = 1, nat
                !
                fac = amass( ityp(na) ) * massconv * alat
                !
                taup(:,na) = taup(:,na) - lambda * dg0(:,na) / fac
                !
                force(:,na) = force(:,na) - lambda * dg0(:,na) * invdtsq
                !
             END DO
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
          CALL infomsg( 'check_constrain', 'g = 0 is not satisfied', -1 )
          !
          WRITE( stdout, '(5X,"Number of step(s): ",I3)') MIN( i, maxiter )
          !
       END IF
       !
       DEALLOCATE( dgp )
       DEALLOCATE( dg0 )
       !
       RETURN
       !
     END SUBROUTINE check_constraint
     !
     !-----------------------------------------------------------------------
     SUBROUTINE remove_constraint_force( nat, tau, if_pos, ityp, alat, force )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN)    :: nat
       REAL(DP), INTENT(IN)    :: tau(:,:)
       INTEGER,  INTENT(IN)    :: if_pos(:,:)
       INTEGER,  INTENT(IN)    :: ityp(:)
       REAL(DP), INTENT(IN)    :: alat
       REAL(DP), INTENT(INOUT) :: force(:,:)
       !
       INTEGER               :: index
       REAL(DP)              :: g
       REAL(DP), ALLOCATABLE :: dg(:,:)
       !
       REAL(DP), EXTERNAL :: DDOT, DNRM2
       !
       !
       lagrange(:) = 0.D0
       !
#if defined (__REMOVE_CONSTRAINT_FORCE)
       !
       ALLOCATE( dg( 3, nat ) )
       !
       DO index = 1, nconstr
          !
          CALL constraint_grad( index, nat, tau, if_pos, ityp, alat, g, dg )
          !
          dg = dg / DNRM2( 3 * nat, dg, 1 )
          !
          lagrange(index) = DDOT( 3 * nat, force, 1, dg, 1 )
          !
          force = force - lagrange(index) * dg
          !
       END DO
       !
       DEALLOCATE( dg )
       !
#endif
       !
     END SUBROUTINE remove_constraint_force     
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
END MODULE constraints_module
