!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE move_ions()
  !----------------------------------------------------------------------------
  !
  ! ... This routine moves the ions according to the requested scheme:
  !
  ! ... iswitch = 1      bfgs minimizations
  ! ... iswitch = 2      constrained bfgs minimization:
  ! ...                  the user must supply the routine 'constrain' which
  ! ...                  defines the constraint equation and the gradient
  ! ...                  the constraint function gv(tau), dgv(i,tau) such
  ! ...                  that:
  !
  ! ...                            gv({tau}) - target = 0,
  !
  ! ...                  and
  !
  ! ...                                         D gv( {tau} )
  ! ...                            dgv(i,na) = ---------------.
  ! ...                                         D tau(i,na)
  !
  ! ... iswitch = 3      molecular dynamics, ( verlet of vcsmd )
  ! ... iswitch = 4      molecular dynamics with one constraint,
  ! ...                  the same conventions as iswitch = 2
  !
  ! ... coefficients for potential and wavefunctions extrapolation are
  ! ... also computed here
  !
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : tmp_dir, prefix, iunupdate
  USE bfgs_module,   ONLY : lbfgs_ndim, new_bfgs => bfgs, lin_bfgs
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : alat, at, bg
  USE basis,         ONLY : nat, ityp, tau, atm
  USE gvect,         ONLY : nr1, nr2, nr3
  USE klist,         ONLY : nelec
  USE symme,         ONLY : s, ftau, nsym, irt
  USE ener,          ONLY : etot
  USE force_mod,     ONLY : force
  USE control_flags, ONLY : upscale, lbfgs, loldbfgs, lconstrain, &
                            lmd, conv_ions, history, alpha0, beta0, tr2
  USE relax,         ONLY : epse, epsf, starting_scf_threshold
  USE cellmd,        ONLY : lmovecell, calc
  USE mp_global,     ONLY : intra_image_comm
  USE io_global,     ONLY : ionode_id, ionode
  USE mp,            ONLY : mp_bcast
  !
  ! ... external procedures
  !
  USE constraints_module,     ONLY : dist_constrain, check_constrain, &
                                     new_force, compute_penalty
  USE basic_algebra_routines, ONLY : norm
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  REAL(KIND=DP), ALLOCATABLE :: tauold(:,:,:)
    ! previous positions of atoms  
  REAL(KIND=DP), SAVE        :: lambda = 0.5D0    
  INTEGER                    :: na  
  REAL(KIND=DP)              :: energy_error, gradient_error
  LOGICAL                    :: step_accepted, exst
  REAL(KIND=DP), ALLOCATABLE :: pos(:), gradient(:)
  !
  !
  ! ... only one node does the calculation in the parallel case
  !
  IF ( ionode ) THEN 
     !
     conv_ions = .FALSE.
     !
     ALLOCATE( tauold( 3, nat, 3 ) )   
     !
     ! ... constrains are imposed here
     !  
     IF ( lconstrain ) CALL impose_constrains()
     !
     ! ... the file containing old positions is opened 
     ! ... ( needed for extrapolation )
     !
     CALL seqopn( iunupdate, TRIM( prefix ) // '.update', 'FORMATTED', exst ) 
     !
     IF ( exst ) THEN
        !
        READ( UNIT = iunupdate, FMT = * ) history
        READ( UNIT = iunupdate, FMT = * ) tauold
        !
     ELSE
        !
        history = 0
        tauold  = 0.D0
        !
     END IF
     !
     CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
     !
     ! ... save the previous two steps ( a total of three steps is saved )
     !
     tauold(:,:,3) = tauold(:,:,2)
     tauold(:,:,2) = tauold(:,:,1)
     tauold(:,:,1) = tau(:,:)
     !
     ! ... do the minimization / dynamics step
     !
     IF ( lmovecell .AND. lconstrain ) &
        CALL errore( 'move_ions', &
                   & 'variable cell and constrain not implemented', 1 )
     !
     ! ... BFGS algorithm is used to minimize ionic configuration
     !
     IF ( lbfgs ) THEN
        !
        ! ... the new bfgs procedure is used
        !  
        ALLOCATE( pos( 3 * nat ) )
        ALLOCATE( gradient( 3 * nat ) )
        !
        pos      =   RESHAPE( SOURCE = tau,   SHAPE = (/ 3 * nat /) ) * alat
        gradient = - RESHAPE( SOURCE = force, SHAPE = (/ 3 * nat /) )
        !
        IF ( lbfgs_ndim == 1 ) THEN
           !
           CALL new_bfgs( pos, etot, gradient, tmp_dir, stdout, epse,        &
                          epsf, energy_error, gradient_error, step_accepted, &
                          conv_ions )
           !
        ELSE
           !
           CALL lin_bfgs( pos, etot, gradient, tmp_dir, stdout, epse,        &
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
        tau   =   RESHAPE( SOURCE = pos, SHAPE = (/ 3, nat /) ) / alat
        force = - RESHAPE( SOURCE = gradient, SHAPE = (/ 3, nat /) )
        !
        CALL output_tau( conv_ions )
        !
        DEALLOCATE( pos )
        DEALLOCATE( gradient ) 
        !
     ELSE IF ( loldbfgs ) THEN
        !
        ! ... the old bfgs scheme is used
        !
        CALL bfgs()
        !   
     END IF
     !
     ! ... molecular dynamics schemes are used
     !
     IF ( lmd ) THEN
        !
        IF ( calc == ' ' ) CALL dynamics()  ! verlet dynamics
        IF ( calc /= ' ' ) CALL vcsmd()     ! variable cell shape md
        !
     END IF
     !
     ! ... check if the new positions satisfy the constrain equation
     !
     IF ( lconstrain ) CALL check_constrain()
     !
     ! ... before leaving check that the new positions still transform
     ! ... according to the symmetry of the system.
     !
     CALL checkallsym( nsym, s, nat, tau, ityp, &
                       at, bg, nr1, nr2, nr3, irt, ftau )
     !
     ! ... find the best coefficients for the extrapolation of the potential
     !
     CALL find_alpha_and_beta( nat, tau, tauold, alpha0, beta0 )
     !
     history = MIN( 3, ( history + 1 ) )
     !
     CALL seqopn( iunupdate, TRIM( prefix ) // '.update', 'FORMATTED', exst ) 
     !
     WRITE( UNIT = iunupdate, FMT = * ) history
     WRITE( UNIT = iunupdate, FMT = * ) tauold
     !
     CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
     !  
     DEALLOCATE( tauold )
     !
  END IF
  !  
  ! ... broadcast calculated quantities to all nodes
  !
  CALL mp_bcast( conv_ions, ionode_id, intra_image_comm )
  CALL mp_bcast( tau,       ionode_id, intra_image_comm )
  CALL mp_bcast( force,     ionode_id, intra_image_comm )
  CALL mp_bcast( tr2,       ionode_id, intra_image_comm )
  CALL mp_bcast( conv_ions, ionode_id, intra_image_comm )
  CALL mp_bcast( alpha0,    ionode_id, intra_image_comm )
  CALL mp_bcast( beta0,     ionode_id, intra_image_comm )
  CALL mp_bcast( history,   ionode_id, intra_image_comm )
  ! 
  RETURN
  !
  CONTAINS
     !
     ! ... internal procedures   
     !  
     !-----------------------------------------------------------------------
     SUBROUTINE impose_constrains()
       !-----------------------------------------------------------------------
       !
       USE constraints_module, ONLY : nconstr
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER       :: index, na
       REAL(KIND=DP) :: gv
       REAL(KIND=DP) :: dgv(3,nat)
       REAL(KIND=DP) :: dgv2
         ! gv = 0 defines the constrain
         ! the gradient of gv
         ! its square modulus       
       !
       !
       IF ( lbfgs ) THEN
          !
          ! ... BFGS case: a penalty function is used
          !
          CALL compute_penalty( gv, dgv, dgv2 )
          !
          etot = etot + lambda * gv**2
          !
          force(:,:) = force(:,:) - 2.D0 * lambda * gv * dgv(:,:)
          !
       ELSE IF ( lmd ) THEN
          !
          ! ... molecular dynamics case: lagrange multipliers are used
          !
          ! ... find the constrained forces
          !
          DO index = 1, nconstr
             !
             CALL dist_constrain( index, gv, dgv, dgv2 )
             !
             CALL new_force( dgv, dgv2 )
             !
          END DO
          !
          WRITE( stdout, '(/5x,"Constrained forces")')
          !
          DO na = 1, nat
             !
             WRITE( stdout, '(3F14.8)') force(:,na)
             !
          END DO
          !   
       END IF       
       !
     END SUBROUTINE impose_constrains
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE compute_lambda()
       !-----------------------------------------------------------------------
       !
       USE constraints_module, ONLY : constr_tol
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       LOGICAL       :: ltest       
       REAL(KIND=DP) :: gv
       REAL(KIND=DP) :: dgv(3,nat)
       REAL(KIND=DP) :: dgv2
         ! gv = 0 defines the constrain
         ! the gradient of gv
         ! its square modulus             
       !
       !
       CALL compute_penalty( gv, dgv, dgv2 )
       !
       IF ( step_accepted ) THEN
          !
          lambda_loop: DO
             !
             IF ( ABS( gv ) > constr_tol ) lambda = lambda * 1.1D0        
             !
             ltest = .TRUE. 
             !
             DO na = 1, nat
                !
                IF ( 2.D0 * lambda * gv * norm( dgv(:,na) ) > 0.05D0 ) &
                   ltest = .FALSE.
                !
             END DO
             !
             IF ( ltest ) EXIT lambda_loop 
             !
             lambda = lambda * 0.5D0
             !
          END DO lambda_loop
          !
       END IF
       !
       WRITE( stdout, '("LAMBDA  = ",F14.10)' ) lambda
       WRITE( stdout, '("GV      = ",F14.10)' ) gv 
       WRITE( stdout, '("PENALTY = ",F14.10)' ) lambda * gv**2       
       !
       RETURN
       !
     END SUBROUTINE compute_lambda
     !
END SUBROUTINE move_ions     
!
! ... this routine is used also by compute_scf (NEB)
!
!----------------------------------------------------------------------------
SUBROUTINE find_alpha_and_beta( nat, tau, tauold, alpha0, beta0 )
  !----------------------------------------------------------------------------
  !
  ! ... This routine finds the best coefficients alpha0 and beta0 so that
  !
  ! ...    | tau(t+dt) - tau' | is minimum, where
  !
  ! ...    tau' = tau(t) + alpha0 * ( tau(t) - tau(t-dt) )
  ! ...                  + beta0 * ( tau(t-dt) -tau(t-2*dt) )
  !
  USE constants, ONLY : eps8
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER       :: nat, na, ipol
  REAL(KIND=DP) :: chi, alpha0, beta0, tau(3,nat), tauold(3,nat,3)
  REAL(KIND=DP) :: a11, a12, a21, a22, b1, b2, c, det
  !
  ! ... solution of the linear system
  !
  a11 = 0.D0
  a12 = 0.D0
  a21 = 0.D0
  a22 = 0.D0
  b1  = 0.D0
  b2  = 0.D0
  c   = 0.D0
  !
  DO na = 1, nat
     !
     DO ipol = 1, 3
        !
        a11 = a11 + ( tauold(ipol,na,1) - tauold(ipol,na,2) )**2
        !
        a12 = a12 + ( tauold(ipol,na,1) - tauold(ipol,na,2) ) * &
                    ( tauold(ipol,na,2) - tauold(ipol,na,3) )
        !
        a22 = a22 + ( tauold(ipol,na,2) - tauold(ipol,na,3) )**2
        !
        b1 = b1 - ( tauold(ipol,na,1) - tau(ipol,na) ) * &
                  ( tauold(ipol,na,1) - tauold(ipol,na,2) )
        !
        b2 = b2 - ( tauold(ipol,na,1) - tau(ipol,na) ) * &
                  ( tauold(ipol,na,2) - tauold(ipol,na,3) )
        ! 
        c = c + ( tauold(ipol,na,1) - tau(ipol,na) )**2
        !
     END DO
     !
  END DO
  !
  a21 = a12
  !
  det = a11 * a22 - a12 * a21
  !
  IF ( det < 0.D0 ) THEN
     !
     alpha0 = 0.D0
     beta0  = 0.D0
     !
     CALL errore( 'find_alpha_and_beta', ' det < 0', -1 )
     !
  END IF   
  !
  ! ... case det > 0:  a well defined minimum exists
  !
  IF ( det > eps8 ) THEN
     !
     alpha0 = ( b1 * a22 - b2 * a12 ) / det
     beta0  = ( a11 * b2 - a21 * b1 ) / det
     !
  ELSE
     !
     ! ... case det = 0 : the two increments are linearly dependent, 
     ! ...                chose solution with alpha = 0 and beta = 0 
     ! ...                ( discard oldest configuration )
     !
     alpha0 = 0.D0
     beta0  = 0.D0
     !
     IF ( a11 > 0.D0 ) alpha0 = b1 / a11
     !
  END IF
  !
  chi = 0.D0
  !
  DO na = 1, nat
     !
     DO ipol = 1, 3
        !
        chi = chi + ( ( 1.D0  + alpha0 ) * tauold(ipol,na,1) + &
                      ( beta0 - alpha0 ) * tauold(ipol,na,2) - & 
                      beta0 * tauold(ipol, na, 3) - tau(ipol,na) )**2
        !
     END DO
     !
  END DO
  !
#if defined (__DEBUG_EXTR)  
  PRINT *, ""
  PRINT *, "chi   = ", chi,    "  det  = ", det
  PRINT *, "alpha = ", alpha0, "  beta = ", beta0
  PRINT *, ""
  PRINT *, "PREDICTED POSITIONS:"
  PRINT *, tauold(1,:,1) + alpha0 * ( tauold(1,:,1) - tauold(1,:,2) ) + &
                            beta0 * ( tauold(1,:,2) - tauold(1,:,3) )
  PRINT *, ""
#endif  
  !
  RETURN
  !
END SUBROUTINE find_alpha_and_beta
