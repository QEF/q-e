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
  ! ... l(old)bfgs       bfgs minimizations
  ! ... lmd              molecular dynamics ( verlet of vcsmd )
  ! ... lmd+lconstrain   molecular dynamics with one constraint,
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
  ! ... coefficients for potential and wavefunctions extrapolation are
  ! ... also computed here
  !
  USE constants,     ONLY : eps8
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : tmp_dir, prefix, iunupdate
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : alat, at, bg, omega
  USE cellmd,        ONLY : omega_old, at_old, lmovecell, calc
  USE ions_base,     ONLY : nat, ityp, tau, atm
  USE gvect,         ONLY : nr1, nr2, nr3
  USE symme,         ONLY : s, ftau, nsym, irt
  USE ener,          ONLY : etot
  USE force_mod,     ONLY : force
  USE bfgs_module,   ONLY : lbfgs_ndim
  USE control_flags, ONLY : upscale, lbfgs, loldbfgs, lconstrain, &
                            lmd, conv_ions, history, alpha0, beta0, tr2, istep
  USE relax,         ONLY : epse, epsf, starting_scf_threshold
  USE lsda_mod,      ONLY : lsda, absmag
  USE mp_global,     ONLY : intra_image_comm
  USE io_global,     ONLY : ionode_id, ionode
  USE mp,            ONLY : mp_bcast
  !
  ! ... external procedures
  !
  USE bfgs_module,            ONLY : new_bfgs => bfgs, lin_bfgs, terminate_bfgs
  USE constraints_module,     ONLY : dist_constrain, check_constrain, &
                                     new_force
  USE basic_algebra_routines, ONLY : norm
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  LOGICAL, SAVE              :: lcheck_mag = .TRUE.
    ! .TRUE. if the check of zero absolute magnetization is required
  REAL(KIND=DP), ALLOCATABLE :: tauold(:,:,:)
    ! previous positions of atoms
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
     ! ... constraints are imposed here
     !  
     IF ( lconstrain ) CALL impose_constraints()
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
        WRITE( UNIT = iunupdate, FMT = * ) history
        WRITE( UNIT = iunupdate, FMT = * ) tauold
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
        ALLOCATE( pos(      3 * nat ) )
        ALLOCATE( gradient( 3 * nat ) )
        !
        pos      =   RESHAPE( SOURCE = tau,   SHAPE = (/ 3 * nat /) ) * alat
        gradient = - RESHAPE( SOURCE = force, SHAPE = (/ 3 * nat /) )
        !
        IF ( lbfgs_ndim == 1 ) THEN
           !
           ! ... standard BFGS 
           !
           CALL new_bfgs( pos, etot, gradient, tmp_dir, stdout, epse,        &
                          epsf, energy_error, gradient_error, step_accepted, &
                          conv_ions )
           !
        ELSE
           !
           ! ... linear scaling BFGS
           !
           CALL lin_bfgs( pos, etot, gradient, tmp_dir, stdout, epse,        &
                          epsf, energy_error, gradient_error, step_accepted, &
                          conv_ions )
           !
        END IF
        !
        IF ( conv_ions ) THEN
           !
           IF ( ( lsda .AND. ( absmag < eps8 ) .AND. lcheck_mag ) ) THEN
              !
              ! ... lsda relaxation :  a final configuration with zero 
              ! ...                    absolute magnetization has been found
              !
              ! ... here we check if it is really the minimum energy structure
              ! ... by performing a new scf iteration without any "electronic"
              ! ... history
              !
              WRITE( UNIT = stdout, FMT = 9010 )
              WRITE( UNIT = stdout, FMT = 9020 )
              !
              CALL hinit0()
              CALL potinit()
              CALL newd()
              CALL wfcinit()
              !
              ! ... this check is performed only once
              !
              lcheck_mag = .FALSE.
              !
              ! ... conv_ions is set to .FALSE. to perform a final scf cycle
              !
              conv_ions = .FALSE.
              ! 
           ELSE
              !
              CALL terminate_bfgs( etot, stdout, tmp_dir )
              !
           END IF   
           !
        ELSE
           !
           ! ... if a new bfgs step is done, new threshold is computed
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
           ! ... the logical flag lcheck_mag is set again to .TRUE. (needed if 
           ! ... a new configuration with zero zero absolute magnetization is 
           ! ... identified in the following steps of the relaxation)
           !
           lcheck_mag = .TRUE.
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
     ! ... history is updated (a new ionic step has been done)
     !
     history = MIN( 3, ( history + 1 ) )
     !
     ! ... find the best coefficients for the extrapolation of the potential
     !
     CALL find_alpha_and_beta( nat, tau, tauold, alpha0, beta0 )
     !
     ! ... old positions are written on file
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
  CALL mp_bcast( istep,     ionode_id, intra_image_comm )
  CALL mp_bcast( tau,       ionode_id, intra_image_comm )
  CALL mp_bcast( force,     ionode_id, intra_image_comm )
  CALL mp_bcast( tr2,       ionode_id, intra_image_comm )
  CALL mp_bcast( conv_ions, ionode_id, intra_image_comm )
  CALL mp_bcast( alpha0,    ionode_id, intra_image_comm )
  CALL mp_bcast( beta0,     ionode_id, intra_image_comm )
  CALL mp_bcast( history,   ionode_id, intra_image_comm )
  !
  IF (lmovecell) THEN
     !
     CALL mp_bcast( at,        ionode_id, intra_image_comm )
     CALL mp_bcast( at_old,    ionode_id, intra_image_comm )
     CALL mp_bcast( omega,     ionode_id, intra_image_comm )
     CALL mp_bcast( omega_old, ionode_id, intra_image_comm )
     CALL mp_bcast( bg,        ionode_id, intra_image_comm )
     !
  END IF
  ! 
  RETURN
  !
9010 FORMAT( /5X,'lsda relaxation :  a final configuration with zero', &
           & /5X,'                   absolute magnetization has been found' )
9020 FORMAT( /5X,'the program is checking if it is really ', &
           &     'the minimum energy structure',             &
           & /5X,'by performing a new scf iteration',        & 
           &     'without any "electronic" history' )               
  !
  CONTAINS
     !
     ! ... internal procedures   
     !  
     !-----------------------------------------------------------------------
     SUBROUTINE impose_constraints()
       !-----------------------------------------------------------------------
       !
       USE constraints_module, ONLY : nconstr
       USE ions_base,          ONLY : ityp
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
       ! ... molecular dynamics: lagrange multipliers are used
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
     END SUBROUTINE impose_constraints
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
  USE constants,     ONLY : eps16
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE control_flags, ONLY : order, history
  !
  IMPLICIT NONE
  !
  INTEGER       :: nat, na, ipol
  REAL(KIND=DP) :: chi, alpha0, beta0, tau(3,nat), tauold(3,nat,3)
  REAL(KIND=DP) :: a11, a12, a21, a22, b1, b2, c, det
  !
  !
  IF ( MIN( history, order ) < 2 ) THEN
     !
     RETURN
     !
  ELSE IF ( MIN( history, order ) == 2 ) THEN  
     !
     alpha0 = 1.D0
     beta0  = 0.D0
     !
     RETURN
     !
  END IF
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
  IF ( det < - eps16 ) THEN
     !
     alpha0 = 0.D0
     beta0  = 0.D0
     !
     WRITE( UNIT = stdout, &
            FMT = '(5X,"WARNING: in find_alpha_and_beta  det = ",F10.6)' ) det
     !
  END IF   
  !
  ! ... case det > 0:  a well defined minimum exists
  !
  IF ( det > eps16 ) THEN
     !
     alpha0 = ( b1 * a22 - b2 * a12 ) / det
     beta0  = ( a11 * b2 - a21 * b1 ) / det
     !
  ELSE
     !
     ! ... case det = 0 : the two increments are linearly dependent, 
     ! ...                chose solution with alpha = b1 / a11 and beta = 0 
     ! ...                ( discard oldest configuration )
     !
     alpha0 = 0.D0
     beta0  = 0.D0
     !
     IF ( a11 /= 0.D0 ) alpha0 = b1 / a11
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE find_alpha_and_beta
