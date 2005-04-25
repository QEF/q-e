!
! Copyright (C) 2004 PWSCF-CP90 groups
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
  ! ... constraint ionic relaxations (based on a penalty function or lagrange
  ! ... multipliers) are defined here.
  !
  ! ... Written by Carlo Sbraccia ( 24/02/2004 )
  !
  ! ... references :
  !
  ! ... 1) M. P. Allen and D. J. Tildesley, Computer Simulations of Liquids,
  ! ...    Clarendon Press - Oxford (1986)
  !
  !
  USE kinds, ONLY: DP  
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  ! ... public methods
  !
  PUBLIC :: dist_constrain,  &
            check_constrain, &
            new_force,       &
            remove_constraint_force
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
     SUBROUTINE dist_constrain( index, tau, g, dg, dg2 )
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
       USE constants, ONLY : eps32
       USE cell_base, ONLY : alat
       USE ions_base, ONLY : nat
       !
       IMPLICIT NONE
       !
       INTEGER,        INTENT(IN) :: index
       REAL (KIND=DP), INTENT(IN) :: tau(3,nat)
         ! the index of the required constrain
       REAL (KIND=DP), INTENT(OUT):: dg(3,nat), dg2, g
         ! constrain terms ( in bohr )
       !
       ! ... local variables
       !
       REAL(KIND=DP) :: x1, x2, y1, y2, z1, z2
       REAL(KIND=DP) :: dist0
       INTEGER       :: ia1, ia2
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
          IF ( dist0 > eps32 ) THEN
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
       CASE DEFAULT
          !
          CALL errore( 'dist_constrain', &
                       'constrain type not implemented ', 1 )
          !
       END SELECT
       !
       RETURN
       !
     END SUBROUTINE dist_constrain
     !
     !-----------------------------------------------------------------------
     SUBROUTINE check_constrain( tau )
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
       USE io_global,  ONLY : stdout
       USE constants,  ONLY : eps16
       USE cell_base,  ONLY : alat
       USE ions_base,  ONLY : nat, ityp, atm
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), INTENT(INOUT) :: tau(3,nat)
       !
       INTEGER                    :: na, i, index
       REAL(KIND=DP), ALLOCATABLE :: dg(:,:)
       REAL(KIND=DP)              :: dg2, g
       LOGICAL                    :: ltest(nconstr), global_test
       INTEGER, PARAMETER         :: maxiter = 250
       !
       !
       ALLOCATE( dg(3,nat) )
       !
       outer_loop: DO i = 1, maxiter
          !
          inner_loop: DO index = 1, nconstr
             !
             ltest(index) = .FALSE.             
             !
             CALL dist_constrain( index, tau, g, dg, dg2 )
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
       IF ( .NOT. global_test ) &
          CALL errore( 'check_constrain', 'g = 0 is not satisfied g = ', - 1 )
       !
       WRITE( stdout, '(5X,"Number of step(s): ",I3)') i - 1
       !
       ! ... if the atomic positions have been corrected write them on output
       !
       IF ( i > 1 ) THEN
          !
          WRITE( stdout, '(/5X,"Corrected atomic positions:")')
          !
          CALL output_tau( .FALSE. )
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
     SUBROUTINE new_force( dg, dg2, force )
       !-----------------------------------------------------------------------
       !
       ! ... find the lagrange multiplier lambda for the problem with one 
       ! ... constrain
       !
       ! ...            force * dg
       ! ... lambda = - ---------- ,
       ! ...              |dg|^2
       !
       ! ... and redefine the forces:
       !
       ! ... force = force + lambda * dg
       !
       ! ... where dg is the gradient of the constraint function
       !
       USE io_global, ONLY : stdout
       USE constants, ONLY : eps16, eps32
       USE ions_base, ONLY : nat
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), INTENT(IN)  :: dg(3,nat), dg2
       REAL (KIND=DP), INTENT(OUT) :: force(3,nat)
       !
       INTEGER        :: na, i, ipol
       REAL (KIND=DP) :: lambda, sum
       !
       ! ... external function
       !
       REAL(KIND=DP), EXTERNAL :: DDOT  
       !
       !
       lambda = 0.D0
       !
       IF ( dg2 > eps32 ) THEN
          !
          lambda = - DDOT( 3 * nat, force, 1, dg, 1 ) / dg2
          !
          force = lambda * dg + force
          !
          IF ( DDOT( 3 * nat, force, 1, dg, 1 )**2 > eps32 ) THEN
             !
             CALL errore( 'new_force', &
                        & 'force is not orthogonal to constrain', - 1 )
             WRITE( stdout, * ) DDOT( 3 * nat, force, 1, dg, 1 )**2
             !
          END IF
          !
          DO ipol = 1, 3
             !
             sum = 0.D0
             !
             DO na = 1, nat
                !
                sum = sum + force(ipol,na)
                !
             END DO
             !
             ! ... impose total force = 0
             !
             DO na = 1, nat
                !
                force(ipol,na) = force(ipol,na) - sum / nat
                !
             END DO
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
     SUBROUTINE remove_constraint_force( tau, force )
       !-----------------------------------------------------------------------
       !
       USE io_global, ONLY : stdout
       USE ions_base, ONLY : ityp, nat
       !
       IMPLICIT NONE
       !
       REAL (KIND=DP), INTENT(IN)  :: tau(3,nat)
       REAL (KIND=DP), INTENT(OUT) :: force(3,nat)
       !
       INTEGER        :: index, na
       REAL (KIND=DP) :: gv
       REAL (KIND=DP) :: dgv(3,nat)
       REAL (KIND=DP) :: dgv2
         ! gv = 0 defines the constrain
         ! the gradient of gv
         ! its square modulus       
       !
       ! ... find the constrained forces
       !
       DO index = 1, nconstr
          !
          CALL dist_constrain( index, tau, gv, dgv, dgv2 )
          !
          CALL new_force( dgv, dgv2, force )
          !
       END DO
       !
       WRITE( stdout, '(/,5X,"Constrained forces (Ry/au):",/)')
       !
       DO na = 1, nat
          !
          WRITE( UNIT = stdout, FMT = 9000 ) na, ityp(na), force(:,na)
          !
       END DO
       !
9000 FORMAT(5X,'atom ',I3,' type ',I2,'   force = ',3F14.8) 
       !
     END SUBROUTINE remove_constraint_force
     !
END MODULE constraints_module
