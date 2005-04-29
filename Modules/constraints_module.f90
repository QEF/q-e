!
! Copyright (C) 2004-2005 PWSCF-CP90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
  USE constants, ONLY : eps16
  USE io_global, ONLY : stdout
  !
  USE basic_algebra_routines
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  ! ... public methods
  !
  PUBLIC :: init_constraint, &
            check_constrain, &
            remove_constraint_force, &
            deallocate_constraint
  !
  ! ... public variables (assigned in the CONSTRAINTS input card)
  !
  PUBLIC :: nconstr,     &
            constr_tol,  &
            constr_type, &
            constr,      &
            target
  !
  ! ... global variables
  !
  INTEGER                     :: nconstr 
  REAL (KIND=DP)              :: constr_tol
  INTEGER,        ALLOCATABLE :: constr_type(:)
  INTEGER,        ALLOCATABLE :: constr(:,:)
  REAL (KIND=DP), ALLOCATABLE :: target(:)
  !
  !
  CONTAINS
     !
     ! ... public methods
     !
     !-----------------------------------------------------------------------
     SUBROUTINE init_constraint( nat, tau, if_pos )
       !-----------------------------------------------------------------------
       !
       USE input_parameters, ONLY : nconstr_inp, constr_tol_inp, &
                                    constr_type_inp, constr_inp, &
                                    constr_target, constr_target_set
       !
       IMPLICIT NONE
       !
       INTEGER,        INTENT(IN) :: nat
       REAL (KIND=DP), INTENT(IN) :: tau(3,nat)
       INTEGER,        INTENT(IN) :: if_pos(3,nat)
       !
       INTEGER       :: ia, ia1, ia2, ia3
       REAL(KIND=DP) :: r12(3), r23(3)
       LOGICAL       :: ltest
       !
       nconstr    = nconstr_inp
       constr_tol = constr_tol_inp
       !
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
          IF ( constr_target_set(ia) ) THEN
             !
             target(ia) = constr_target(ia)
             !
             CYCLE
             !
          END IF
          !
          SELECT CASE ( constr_type(ia) )
          CASE( 1 )
             !
             target(ia) = norm( tau(:,constr(1,ia)) - tau(:,constr(2,ia)) )
             !
             ltest = .FALSE.
             ltest = ltest .OR. ANY( if_pos(:,constr(1,ia)) == 0 )
             ltest = ltest .OR. ANY( if_pos(:,constr(2,ia)) == 0 )
             !
          CASE( 2 )
             !
             ia1 = constr(1,ia)
             ia2 = constr(2,ia)
             ia3 = constr(3,ia)
             !
             r12 = tau(:,ia2) - tau(:,ia1)
             r23 = tau(:,ia2) - tau(:,ia3)
             !
             r12 = r12 / norm( r12 )
             r23 = r23 / norm( r23 )
             !
             target(ia) = DOT_PRODUCT( r12, r23 )
             !
             ltest = .FALSE.
             ltest = ltest .OR. ANY( if_pos(:,constr(1,ia)) == 0 )
             ltest = ltest .OR. ANY( if_pos(:,constr(2,ia)) == 0 )
             ltest = ltest .OR. ANY( if_pos(:,constr(3,ia)) == 0 )
             !
          CASE DEFAULT
             !
             CALL errore( 'init_constraint', &
                          'constrain type not implemented ', 1 )
             !
          END SELECT
          !
          IF ( ltest ) &
             CALL errore( 'init_constraint', &
                        & 'constraints cannot be set on fixed atoms', 1 )
          !
       END DO
       !
       RETURN
       !
     END SUBROUTINE init_constraint
     !
     !-----------------------------------------------------------------------
     SUBROUTINE dist_constrain( index, nat, tau, alat, g, dg, dg2 )
       !-----------------------------------------------------------------------
       ! 
       ! ... this routine defines the constrain equation:
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
       INTEGER,        INTENT(IN) :: index
       INTEGER,        INTENT(IN) :: nat
       REAL (KIND=DP), INTENT(IN) :: tau(:,:)
       REAL (KIND=DP), INTENT(IN) :: alat
       REAL (KIND=DP), INTENT(OUT):: dg(:,:)
       REAL (KIND=DP), INTENT(OUT):: dg2, g
         ! constrain terms ( in bohr )
       !
       ! ... local variables
       !
       REAL(KIND=DP) :: x1, x2, y1, y2, z1, z2
       REAL(KIND=DP) :: dist0
       INTEGER       :: ia1, ia2, ia3
       REAL(KIND=DP) :: r12(3), r23(3)
       REAL(KIND=DP) :: norm_r12, norm_r23, cos123, sin123
       !
       ! ... external function
       !
       REAL(KIND=DP), EXTERNAL :: DDOT
       !
       !
       dg(:,:) = 0.D0
       !
       SELECT CASE ( constr_type(index) )
       CASE( 1 )
          !
          ! ... constratint on distances
          !
          ia1 = constr(1,index)
          ia2 = constr(2,index)
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
          IF ( dist0 > eps16 ) THEN
             !           
             dg(1,ia1) = ( x1 - x2 ) / dist0
             dg(1,ia2) = ( x2 - x1 ) / dist0
             dg(2,ia1) = ( y1 - y2 ) / dist0
             dg(2,ia2) = ( y2 - y1 ) / dist0
             dg(3,ia1) = ( z1 - z2 ) / dist0
             dg(3,ia2) = ( z2 - z1 ) / dist0
             !
          END IF
          !
          dg2 = DDOT( 3 * nat, dg, 1, dg, 1 )
          !
       CASE( 2 )
          !
          ia1 = constr(1,index)
          ia2 = constr(2,index)
          ia3 = constr(3,index)
          !
          r12 = tau(:,ia2) - tau(:,ia1)
          r23 = tau(:,ia2) - tau(:,ia3)
          !
          norm_r12 = norm( r12 )
          norm_r23 = norm( r23 )
          !
          r12 = r12 / norm_r12
          r23 = r23 / norm_r23
          !
          cos123 = DOT_PRODUCT( r12, r23 )
          sin123 = SQRT( 1.D0 - cos123**2 )
          !
          g = ( cos123 - target(index) )
          !
          dg(:,ia1) = ( cos123 * r12 - r23 ) / ( sin123 * norm_r12 )
          dg(:,ia3) = ( cos123 * r23 - r12 ) / ( sin123 * norm_r23 )
          dg(:,ia2) = - dg(:,ia1) - dg(:,ia3)
          !
       CASE DEFAULT
          !
          CALL errore( 'dist_constrain', &
                       'constrain type not implemented ', 1 )
          !
       END SELECT
       !
       dg2 = DDOT( 3 * nat, dg, 1, dg, 1 )
       !
       RETURN
       !
     END SUBROUTINE dist_constrain
     !
     !-----------------------------------------------------------------------
     SUBROUTINE check_constrain( nat, tau, alat )
       !-----------------------------------------------------------------------
       !
       ! ... update tau so that the constraint equation g=0 is satisfied,
       ! ... use the recursion formula:
       !
       ! ...                  g(tau)
       ! ... tau' = tau -  ------------ * dg(tau)
       ! ...                |dg(tau)|^2
       !
       ! ... in normal cases the constraint equation should be always 
       ! ... satisfied the very first iteration.
       !
       IMPLICIT NONE
       !
       INTEGER,        INTENT(IN)    :: nat
       REAL (KIND=DP), INTENT(INOUT) :: tau(:,:)
       REAL (KIND=DP), INTENT(IN)    :: alat
       !
       INTEGER                    :: na, i, index
       REAL(KIND=DP), ALLOCATABLE :: dg(:,:)
       REAL(KIND=DP)              :: dg2, g
       LOGICAL                    :: ltest(nconstr), global_test
       INTEGER, PARAMETER         :: maxiter = 100
       !
       !
       ALLOCATE( dg( 3, nat ) )
       !
       outer_loop: DO i = 1, maxiter
          !
          inner_loop: DO index = 1, nconstr
             !
             ltest(index) = .FALSE.             
             !
             CALL dist_constrain( index, nat, tau, alat, g, dg, dg2 )
             !
             ! ... check if g = 0
             !
#if defined (__DEBUG_CONSTRAINTS)
             WRITE( stdout, * ) i, index, ABS( g )
#endif
             !             
             IF ( ABS( g ) < constr_tol ) THEN
                !
                ltest(index) = .TRUE.
                !
                CYCLE inner_loop
                !
             END IF   
             !
             ! ... if  g <> 0  find new  tau = tau - g * dg / dg2  and 
             ! ... check again ( g is in bohr and tau in alat units )
             !
             tau(:,:) = tau(:,:) - g * dg(:,:) / dg2 / alat
             !
          END DO inner_loop
          !
          global_test = .TRUE.
          !
          DO index = 1, nconstr
            !
            global_test = global_test .AND. ltest(index)
            !
          END DO
          !
          ! ... all constraints are satisfied
          !
          IF ( global_test ) EXIT outer_loop
          !
       END DO outer_loop
       !
       IF ( .NOT. global_test ) THEN
          !
          CALL errore( 'check_constrain', 'g = 0 is not satisfied g = ', - 1 )
          !
          WRITE( stdout, '(5X,"Number of step(s): ",I3)') MIN( i, maxiter )
          !
       END IF
       !
       DEALLOCATE( dg )
       !
       RETURN
       !
     END SUBROUTINE check_constrain         
     !
     !-----------------------------------------------------------------------
     SUBROUTINE new_force( nat, dg, dg2, force )
       !-----------------------------------------------------------------------
       !
       ! ... find the lagrange multiplier lambda for the problem with one 
       ! ... constrain
       !
       ! ...          force * dg
       ! ... lambda = ---------- ,
       ! ...            |dg|^2
       !
       ! ... and redefine the forces:
       !
       ! ... force = force - lambda * dg
       !
       ! ... where dg is the gradient of the constraint function
       !
       IMPLICIT NONE
       !
       INTEGER,        INTENT(IN)    :: nat
       REAL (KIND=DP), INTENT(IN)    :: dg(:,:), dg2
       REAL (KIND=DP), INTENT(INOUT) :: force(:,:)
       !
       INTEGER        :: na, i, ipol
       REAL (KIND=DP) :: lambda, sum(3)
       !
       ! ... external function
       !
       REAL(KIND=DP), EXTERNAL :: DDOT  
       !
       !
       lambda = 0.D0
       !
       IF ( dg2 > eps16 ) THEN
          !
          lambda = DDOT( 3 * nat, force, 1, dg, 1 ) / dg2
          !
          force = force - lambda * dg
          !
          IF ( DDOT( 3 * nat, force, 1, dg, 1 )**2 > eps16 ) THEN
             !
             CALL errore( 'new_force', &
                        & 'force is not orthogonal to constrain', 1 )
             WRITE( stdout, * ) DDOT( 3 * nat, force, 1, dg, 1 )**2
             !
          END IF
          !
          sum(:) = 0.D0
          !
          DO na = 1, nat
             !
             sum(:) = sum(:) + force(:,na)
             !
          END DO
          !
          ! ... impose total force = 0
          !
          DO na = 1, nat
             !
             force(:,na) = force(:,na) - sum(:) / nat
             !
          END DO
          !
       END IF
       !
       RETURN
       !
     END SUBROUTINE new_force
     !
     !-----------------------------------------------------------------------
     SUBROUTINE remove_constraint_force( nat, tau, alat, force )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       INTEGER,        INTENT(IN)    :: nat
       REAL (KIND=DP), INTENT(IN)    :: tau(:,:)
       REAL (KIND=DP), INTENT(IN)    :: alat
       REAL (KIND=DP), INTENT(INOUT) :: force(:,:)
       !
       INTEGER                     :: index, na
       REAL (KIND=DP)              :: gv
       REAL (KIND=DP), ALLOCATABLE :: dgv(:,:)
       REAL (KIND=DP)              :: dgv2
         ! gv = 0 defines the constrain
         ! the gradient of gv
         ! its square modulus       
       !
       !
       ALLOCATE( dgv( 3, nat ) )
       !
       ! ... find the constrained forces
       !
       DO index = 1, nconstr
          !
          CALL dist_constrain( index, nat, tau, alat, gv, dgv, dgv2 )
          !
          CALL new_force( nat, dgv, dgv2, force )
          !
       END DO
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
       IF ( ALLOCATED( constr ) )      DEALLOCATE( constr )
       IF ( ALLOCATED( constr_type ) ) DEALLOCATE( constr_type )
       IF ( ALLOCATED( target ) )      DEALLOCATE( target )       
       !
       RETURN
       !
     END SUBROUTINE deallocate_constraint
     !
END MODULE constraints_module
