!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE constraints_module
  !----------------------------------------------------------------------------
  ! 
  ! ... variables and methods for constraint Molecular Dynamics and
  ! ... constraint ionic relaxations (based on a penalty function) are 
  ! ... defined here.
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
            compute_penalty
  !
  ! ... public variables (assigned in the CONSTRAINTS input card)
  !
  PUBLIC :: nconstr,    &
            constr_tol, &
            constr,     &
            target
  !
  ! ... global variables
  !
  INTEGER                     :: nconstr 
  REAL (KIND=DP)              :: constr_tol
  INTEGER,        ALLOCATABLE :: constr(:,:)
  REAL (KIND=DP), ALLOCATABLE :: target(:)
  !
  !
  CONTAINS
     !
     ! ... public methods
     !
     !-----------------------------------------------------------------------
     SUBROUTINE dist_constrain( index, g, dg, dg2 )
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
       USE ions_base, ONLY : nat, tau
       !
       IMPLICIT NONE
       !
       INTEGER,       INTENT(IN) :: index
         ! the index of the required constrain
       REAL(KIND=DP), INTENT(OUT):: dg(3,nat), dg2, g
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
       RETURN
       !
     END SUBROUTINE dist_constrain
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE check_constrain( )
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
       USE ions_base,  ONLY : nat, ityp, tau, atm
       !
       IMPLICIT NONE
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
             CALL dist_constrain( index, g, dg, dg2 )
             !
             ! ... check if g = 0
             !
#if defined (__DEBUG_CONSTRAINS)
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
          IF ( global_test ) EXIT outer_loop
          !
       END DO outer_loop
       !
       IF ( .NOT. global_test ) &
          CALL errore( 'new_dtau', 'g = 0 is not satisfied g = ', - 1 )
       !
       WRITE( stdout, '(5X,"Number of step(s): ",I3)') i - 1
       !
       ! ... if the atomic positions have been corrected write them on output
       !
       IF ( i > 1 ) THEN
          !
          WRITE( stdout, '(/5X,"Corrected atomic positions:",/)')
          !
          DO na = 1, nat
             !
             WRITE( stdout,'(A3,3X,3F14.9)') atm(ityp(na)), tau(:,na)
             !
          END DO
          !
       END IF
       !
       DEALLOCATE( dg )
       !
       RETURN
       !
     END SUBROUTINE check_constrain         
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE new_force( dg, dg2 )
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
       USE cell_base, ONLY : at, bg
       USE force_mod, ONLY : force
       USE symme,     ONLY : s, nsym, irt
       !
       IMPLICIT NONE
       !
       INTEGER       :: na, i, ipol
       REAL(KIND=DP) :: dg(3,nat), lambda, dg2, sum
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
          ! ... resymmetrize (should not be needed, but...)
          !
          IF ( nsym > 1 ) THEN
             !
             DO na = 1, nat
                !
                CALL trnvect( force(1,na), at, bg, - 1 )
                !
             END DO
             !
             CALL symvect( nat, force, nsym, s, irt )
             !
             DO na = 1, nat
                !
                CALL trnvect( force(1,na), at, bg, 1 )
                !
             END DO
             !
          END IF
          !
       END IF
       !
       RETURN
       !
     END SUBROUTINE new_force 
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE compute_penalty( g, dg, dg2 )
       !-----------------------------------------------------------------------
       ! 
       ! ... this routine defines the penalty equation:
       !
       ! ...  g(tau,dist) = 0
       !
       ! ... where tau are the atomic positions ( in alat units ) and dist is
       ! ... the distance of two atoms ( in this case atom 1 and atom 2 ) which 
       ! ... is, in this case, a one dimensional constrain. 
       ! ... dg is in output the value of the gradient of g and dg2 is its 
       ! ... square modulus.
       !
       USE constants, ONLY : eps32
       USE cell_base, ONLY : alat
       USE ions_base, ONLY : nat, tau
       !
       IMPLICIT NONE
       !
       REAL(KIND=DP), INTENT(OUT):: dg(3,nat), dg2, g
         ! constrain terms ( in bohr )
       !
       ! ... local variables
       !
       REAL(KIND=DP) :: x1, x2, y1, y2, z1, z2
       REAL(KIND=DP) :: dist0
       INTEGER       :: ia1, ia2, index
       !
       ! ... external function
       !
       REAL(KIND=DP), EXTERNAL :: DDOT
       !
       !
       dg(:,:) = 0.D0
       g       = 0.D0
       !
       DO index = 1, nconstr
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
          g = g + ( dist0 - target(index) )
          !
          IF ( dist0 > eps32 ) THEN
             !           
             dg(1,ia1) = dg(1,ia1) + ( x1 - x2 ) / dist0
             dg(1,ia2) = dg(1,ia2) + ( x2 - x1 ) / dist0
             dg(2,ia1) = dg(2,ia1) + ( y1 - y2 ) / dist0
             dg(2,ia2) = dg(2,ia2) + ( y2 - y1 ) / dist0
             dg(3,ia1) = dg(3,ia1) + ( z1 - z2 ) / dist0
             dg(3,ia2) = dg(3,ia2) + ( z2 - z1 ) / dist0
             !
          END IF   
          !
       END DO
       !
       dg2 = DDOT( 3 * nat, dg, 1, dg, 1 )
       !
       RETURN
       !
     END SUBROUTINE compute_penalty
     !
END MODULE constraints_module
