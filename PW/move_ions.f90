!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE move_ions()
  !----------------------------------------------------------------------------
  !
  !     This routine moves the ions according to the requested scheme:
  !
  !     iswitch = 1      bfgs minimizations or conjugate gradient
  !     iswitch = 2      constrained bfgs minimization:
  !                      the user must supply the routine 'constrain' which
  !                      defines the constraint equation and the gradient
  !                      the constraint function gv(tau), dgv(i,tau) such
  !                      that:
  !
  !                                gv({tau}) - target = 0,
  !
  !                      and
  !
  !                                             D gv( {tau} )
  !                                dgv(i,na) = ---------------.
  !                                             D tau(i,na)
  !
  !     iswitch = 3      molecular dynamics, (verlet of vcsmd)
  !     iswitch = 4      molecular dynamics with one constraint,
  !                      the same conventions as iswitch = 2
  !
  USE io_global,   ONLY : stdout
  USE io_files,    ONLY : tmp_dir
  USE bfgs_module, ONLY : new_bfgs => bfgs, lin_bfgs
  USE bfgs_module, ONLY : lbfgs_ndim
  USE parameters,  ONLY : DP
  USE brilz,       ONLY : alat, at, bg
  USE basis,       ONLY : nat, ityp, tau, atm
  USE gvect,       ONLY : nr1, nr2, nr3
  USE klist,       ONLY : nelec
  USE symme,       ONLY : s, ftau, nsym, irt
  USE ener,        ONLY : etot
  USE force_mod,   ONLY : force
  USE varie,       ONLY : tr2, upscale, iswitch, lbfgs, lnewbfgs, &
                          conv_ions
  USE relax,       ONLY : epse, epsf, starting_scf_threshold
  USE cellmd,      ONLY : lmovecell, calc
#if defined (__PARA)
  USE para,        ONLY : me, mypool
  USE io_global,   ONLY : ionode_id
  USE mp,          ONLY : mp_bcast
#endif
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP)              :: dummy, gv
  REAL(KIND=DP), ALLOCATABLE :: dgv (:,:)
  REAL(KIND=DP)              :: dgv2, theta0
  ! auxiliar variable :
  ! gv=0 defines the constrain
  ! the gradient of gv
  ! its square modulus
  ! the value of the one-dimensional const
  REAL(KIND=DP)              :: energy_error, gradient_error
  LOGICAL                    :: step_accepted
  REAL(KIND=DP), ALLOCATABLE :: pos(:), gradient(:)
  !
  !
  conv_ions = .FALSE.
  !
  IF ( iswitch == 2 .OR. iswitch == 4 ) THEN
     !
     ALLOCATE( dgv(3,nat) )
     !
     ! ... gv is the function which define the constrain, now first of all we
     ! ... find the constrain theta0 such that gv=0, and find the gradient of
     ! ... gv, dgv
     !
     dummy = 0.D0
     !
     CALL constrain( theta0, gv, dgv, dgv2, dummy, nat, tau, alat )
     !
     ! ... find the constrained forces
     !
     CALL new_force( dgv, dgv2 )
     !
     DEALLOCATE( dgv )
     !
  END IF
  !
  ! ... do the minimization / dynamics step
  !
  IF ( lmovecell .AND. ( iswitch == 2 .OR. iswitch == 4 ) ) &
     CALL errore( 'move_ions', &
                & 'variable cell and constrain not implemented', 1 )
  !
  IF ( lnewbfgs ) THEN
     !
     ! ... the new bfgs procedure is used
     !
#if defined (__PARA)
     !
     ! ... only one node does the calculation in the parallel case
     !
     IF ( me == 1 .AND. mypool == 1 ) THEN 
#endif     
        !
        ALLOCATE( pos( 3 * nat ) )
        ALLOCATE( gradient( 3 * nat ) )
        !
        pos      =   RESHAPE( SOURCE = tau, SHAPE = (/ 3 * nat /) ) * alat
        gradient = - RESHAPE( SOURCE = force, SHAPE = (/ 3 * nat /) )
        !
        IF ( lbfgs_ndim == 1 ) THEN
           !
           CALL new_bfgs( pos, etot, gradient, tmp_dir, stdout, epse, epsf, &
                      energy_error, gradient_error, step_accepted, conv_ions )
           !
        ELSE
           ! 
           CALL lin_bfgs( pos, etot, gradient, tmp_dir, stdout, epse, &
                          epsf, energy_error, gradient_error, step_accepted, &
                          conv_ions )
           !
        END IF
        !
        IF ( .NOT. conv_ions ) THEN
           !
           ! ... if a new bfgs step is done, new thresholds are computed
           !
           IF ( step_accepted ) THEN
              !
              tr2  = starting_scf_threshold * &
                     MIN( 1.D0, ( energy_error / ( epse * upscale ) ), &
                                ( gradient_error / ( epsf * upscale ) ) )
              tr2  = MAX( ( starting_scf_threshold / upscale ), tr2 ) 
              !       
           END IF       
           !
           WRITE( stdout, '(5X,"new conv_thr",T30,"= ",F18.10,/)' ) tr2
           !
        END IF
        !   
        tau   =   RESHAPE( SOURCE = pos, SHAPE = (/ 3 , nat /) ) / alat
        force = - RESHAPE( SOURCE = gradient, SHAPE = (/ 3 , nat /) )
        !
        CALL output_tau()
        !
        DEALLOCATE( pos )
        DEALLOCATE( gradient ) 
#if defined (__PARA)
        !
     END IF
     !  
     ! ... broadcast calculated quantities to all nodes
     !
     CALL mp_bcast( tau, ionode_id )
     CALL mp_bcast( force, ionode_id )
     CALL mp_bcast( tr2, ionode_id )
     CALL mp_bcast( conv_ions, ionode_id )
#endif 
     !
  ELSE
    !
    IF ( iswitch == 1 .OR. iswitch == 2 ) CALL bfgs()
    !   
  END IF   
  !
  IF ( iswitch == 3 .OR.iswitch == 4 ) THEN
     !
     IF ( calc == ' ' ) CALL dynamics()  ! verlet dynamics
     IF ( calc /= ' ' ) CALL vcsmd()     ! variable cell shape md
     !
  END IF
  !
  IF ( iswitch > 4 .OR. iswitch <= 0 ) THEN
     !
     CALL errore( 'move_ions', 'iswitch value not implemented or wrong', 1 )
     !
  END IF
  !
  ! ... check if the new positions satisfy the constrain equation, in
  ! ... the CP case this is done inside the routine "cp"
  !
  IF ( iswitch == 2 .OR. iswitch == 4 ) &
     CALL check_constrain( alat, tau, atm, ityp, theta0, nat )
  !
  ! ... before leaving check that the new positions still transform
  ! ... according to the symmetry of the system.
  !
  CALL checkallsym( nsym, s, nat, tau, ityp, at, bg, nr1, nr2, nr3, irt, ftau )
  !
  RETURN
  !
END SUBROUTINE move_ions
!
!
!----------------------------------------------------------------------------
SUBROUTINE new_force( dg, dg2 )
  !----------------------------------------------------------------------------
  !
  !     find the lagrange multiplier lambda for the problem with one const
  !
  !                force*dg
  !     lambda = - --------,
  !                 |dg|^2
  !
  !     and redefine the forces:
  !
  !     force = force + lambda*dg
  !
  !     where dg is the gradient of the constraint function
  !
  USE io_global,  ONLY : stdout
  USE pwcom
  !
  IMPLICIT NONE
  !
  INTEGER       :: na, i, ipol
  REAL(KIND=DP) :: dg (3, nat), lambda, dg2, sum
  REAL(KIND=DP) :: DDOT
  EXTERNAL         DDOT
  !
  !
  lambda = 0.D0
  !
  IF ( dg2 /= 0.D0 ) THEN
     !
     lambda = - DDOT( 3 * nat, force, 1, dg, 1 ) / dg2
     !
     CALL DAXPY( 3 * nat, lambda, dg, 1, force, 1)
     !
     IF ( DDOT( 3 * nat, force, 1, dg, 1 )**2 > 1.D-30 ) THEN
        !
        CALL errore( 'new_force', 'force is not orthogonal to constrain', - 1 )
        PRINT *, DDOT( 3 * nat, force, 1, dg, 1 )**2
        !
     END IF
     !
     DO ipol = 1, 3
        sum = 0.D0
        DO na = 1, nat
           sum = sum + force(ipol,na)
        END DO
        !
        ! ... impose total force = 0
        !
        DO na = 1, nat
           force(ipol,na) = force(ipol,na) - sum / nat
        END DO
     END DO
     !
     ! ... resymmetrize (should not be needed, but...)
     !
     IF ( nsym > 1 ) THEN
        !
        DO na = 1, nat
           CALL trnvect( force(1,na), at, bg, - 1 )
        END DO
        !
        CALL symvect( nat, force, nsym, s, irt )
        !
        DO na = 1, nat
           CALL trnvect( force(1,na), at, bg, 1 )
        END DO
        !
     END IF
     !
     WRITE( stdout, '(/5x,"Constrained forces")')
     !
     DO na = 1, nat
        WRITE( stdout, '(3F14.8)') ( force(i,na) , i = 1, 3 )
     END DO
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE new_force
!
!
!---------------------------------------------------------------------
SUBROUTINE check_constrain( alat, tau, atm, ityp, theta0, nat )
  !---------------------------------------------------------------------
  !
  !     update tau so that the constraint equation g=0 is satisfied,
  !     use the recursion formula:
  !
  !                      g(tau)
  !     tau' = tau -  ------------ * dg(tau)
  !                    |dg(tau)|^2
  !
  !     in normal cases the constraint equation should be always satisfied
  !     the very first iteration.
  !
  USE io_global,  ONLY : stdout
  USE parameters
  !
  IMPLICIT NONE
  !
  INTEGER                    :: ityp(:), nat, na, i, maxiter
  CHARACTER(LEN=3)           :: atm(:)
  REAL(KIND=DP)              :: tau (3,nat)
  REAL(KIND=DP), ALLOCATABLE :: dg (:,:)
  REAL(KIND=DP)              :: alat, dg2, g, theta0, dummy, eps
  !
  PARAMETER ( eps = 1.D-15, maxiter = 250 )
  !
  !
  ALLOCATE( dg(3,nat) )
  !
  CALL constrain( dummy, g, dg, dg2, theta0, nat, tau, alat )
  !
  WRITE( stdout, '(5X,"G = ",1PE9.2," iteration # ",I3)' ) g, 0
  !
  DO i = 1, maxiter
     !
     ! ... check if g=0
     !
     IF ( ABS( g ) < eps ) GO TO 14
     !
     ! ... if g<>0 find new tau = tau - g*dg/dg2 and check again
     !
     CALL DAXPY( 3 * nat, - g / dg2, dg, 1, tau, 1 )
     !
     CALL constrain( dummy, g, dg, dg2, theta0, nat, tau, alat )
     !
     WRITE( stdout, '(5X,"G = ",1PE9.2," iteration # ",I3)' ) g, i
     !
  END DO
  !
  CALL errore( 'new_dtau', 'g=0 is not satisfied g=', - 1 )
  !
14 CONTINUE
  !
  !     WRITE( stdout,'(5x,"G = ",1pe9.2)')g
  WRITE( stdout, '(5X,"Number of step(s): ",I3)') i - 1
  !
  ! ... if the atomic positions have been corrected write them on output
  !
  IF ( i > 1 ) THEN
     !
     WRITE( stdout, '(/5X,"Corrected atomic positions:",/)')
     DO na = 1, nat
        WRITE( stdout,'(A3,3X,3F14.9)') atm(ityp(na)), ( tau(i,na), i = 1, 3 )
     END DO
     !
  END IF
  !
  DEALLOCATE( dg )
  !
  RETURN
  !
END SUBROUTINE check_constrain
