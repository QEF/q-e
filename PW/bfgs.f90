!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!-----------------------------------------------------------------------
SUBROUTINE bfgs()
  !-----------------------------------------------------------------------
  !
  ! ... ionic relaxation through broyden-fletcher-goldfarb-shanno 
  ! ... minimization
  ! ... this version saves data at each iteration
  !
  USE kinds,         ONLY : DP
  USE cell_base,     ONLY : alat
  USE ions_base,     ONLY : nat, tau
  USE force_mod,     ONLY : force
  USE control_flags, ONLY : conv_ions, upscale, imix, tr2, ethr, istep
  USE relax,         ONLY : restart_bfgs, epse, epsf, starting_diag_threshold, &
                            starting_scf_threshold, dtau_ref
  USE ener,          ONLY : etot
  USE klist,         ONLY : nelec
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : prefix
  !
  IMPLICIT NONE
  !
  INTEGER :: iunit          ! unit for file containing bfgs info
  INTEGER :: nax3,        & ! 3 times the total number of atoms (nat)
       na, i          ! counters
  REAL(KIND=DP), ALLOCATABLE :: &
       hessm1 (:,:),& ! current estimate of hessian^-1
       dtau (:,:),  & ! direction (versor) for line minimization
       oldforce (:,:) ! gradient along minimization direction
  REAL(KIND=DP) :: &
       detot,       & ! forces at previous line minimum (LM)
       deold,       & ! as above at previous LM
       eold,        & ! energy at previous LM
       xold,        & ! position along dtau at previous LM
       x,           & ! present position along dtau
       xnew           ! next position along dtau

  LOGICAL :: exst,       & ! test variable on existence of bfgs file
       minimum_ok    ! true if linmin found a good line minimum
  !
  REAL(KIND=DP) :: DDOT
  !
  !
  ! ... the constrain on fixed coordinates is implemented setting to zero
  ! ... the forces in forces.f90  ( C.S. 15/10/2003 )
  !
  nax3 = 3 * nat
  conv_ions = .FALSE.
  !
  ALLOCATE( hessm1(nax3,nax3), dtau(3,nat), oldforce(3,nat) )
  !
  iunit = 4
  !
  CALL seqopn( iunit, TRIM( prefix )//'.bfgs', 'UNFORMATTED', exst )
  !
  ! ... exst flags whether restarting from preceding iterations
  ! ... do not restart from existing data unless explicitely required
  !
  exst = ( exst .AND. restart_bfgs )
  !
  IF ( .NOT. exst ) THEN
     !
     ! ... file not found: starting iteration
     !
     CLOSE( UNIT = iunit, STATUS = 'DELETE' )
     minimum_ok = .FALSE.
     !
     CALL estimate( hessm1, nax3, nat, nax3 )
     !
     WRITE( stdout, '(/5X,"EPSE = ",E9.2,"    EPSF = ",E9.2, &
          & "    UPSCALE = ",F6.2)') epse, epsf, upscale
     !
  ELSE
     !
     ! ... file found: restart from preceding iterations
     !
     READ( iunit ) &
          minimum_ok, xnew, starting_scf_threshold, starting_diag_threshold
     READ( iunit ) dtau
     READ( iunit ) hessm1
     READ( iunit ) xold, eold, deold, oldforce
     !
     CLOSE( UNIT = iunit, STATUS = 'KEEP' )
     !
  END IF
  !
20 CONTINUE
  !
  IF ( exst .AND. .NOT. minimum_ok ) THEN
     !
     ! ... Search for a new line minimum
     !
     x = xnew
     detot = - DDOT( nax3, force, 1, dtau, 1 )
     !
     ! ... line minimization with 3rd order interpolation formula
     !
     CALL linmin( xold, eold, deold, x, etot, detot, xnew, minimum_ok )
     !
     ! ... xnew close to x: line minimization already achieved
     !
     IF ( ABS( ( xnew - x ) / x ) < 0.05D0 .AND. minimum_ok ) GOTO 20
     !
     ! ... new positions ( hopefully close to the line minimum ) :
     !
     CALL DAXPY( nax3, ( xnew - x ) / alat, dtau, 1, tau, 1 )
     !
     IF ( .NOT. minimum_ok ) THEN
        !
        ! ... line minimum was not found take another step and reset the 
        ! ... starting point of the line minimization to the present position
        !
        xold  = x
        eold  = etot
        deold = detot
        !
        CALL DCOPY( nax3, force, 1, oldforce, 1 )
        !
     END IF
     !
  ELSE
     !
     IF ( exst ) THEN
        !
        ! ... We (hopefully) are at the line minimum: convergence check
        !
        conv_ions = ( ( eold - etot ) < epse )
        !
        DO i = 1, 3
           DO na = 1, nat
              conv_ions = ( conv_ions .AND. ( ABS( force(i,na) ) < epsf ) )
           END DO
        END DO
        !
        ! ... update the inverse hessian
        ! ... set dtau to the true displacements from previous to present LM
        !
        CALL DSCAL( nax3, ( xnew - xold ), dtau, 1 )
        CALL updathes( nax3, nax3, oldforce, force, hessm1, dtau )
        !
     END IF
     !
     ! ... find new minimization direction dtau
     !
     x          = 0.D0
     minimum_ok = .FALSE.
     dtau(:,:)  = 0.D0
     !
     CALL DGEMV( 'N', nax3, nax3, 1.D0, hessm1, nax3, force, 1, &
          0.D0, dtau, 1 )
     xnew = SQRT( DDOT( nax3, dtau, 1, dtau, 1 ) )
     CALL DSCAL( nax3, 1.D0 / xnew, dtau, 1 )
     !
     ! ... and the gradient along the minimization direction dtau
     !
     detot = - DDOT( nax3, force, 1, dtau, 1 )
     !
     IF ( detot > 0.D0 ) THEN
        !
        WRITE( stdout, '("uphill direction! de/dx =",E10.4)') detot
        WRITE( stdout, '("try steepest descent direction instead!")')
        !
        CALL DCOPY( nax3, force, 1, dtau, 1 )
        xnew = SQRT( DDOT( nax3, dtau, 1, dtau, 1 ) )
        CALL DSCAL( nax3, 1.D0 / xnew, dtau, 1 )
        detot = - DDOT( nax3, force, 1, dtau, 1 )
        !
     END IF
     !
     ! ... update atomic positions. NB: tau in units of alat!
     !
     CALL DAXPY( nax3, ( xnew - x ) / alat, dtau, 1, tau, 1 )
     !
     ! ... save values of variables at line minimum for later use
     !
     xold  = x
     eold  = etot
     deold = detot
     !
     CALL DCOPY( nax3, force, 1, oldforce, 1 )
     !
  END IF
  !
  ! ... set appropriate (?) thresholds for self-consistency
  !
  ! ... ethr is now updated in electrons with a different procedure
  ! ... this value of ethr is overwritten apart when the old style
  ! ... update is used (see OLDSTYLE precompiler variable in electrons)
  !     
  IF ( imix < 0 ) THEN
     tr2 = ( starting_scf_threshold * MAX( 1.D0 / upscale, MIN( 1.D0, &
          ABS( xnew - x ) / dtau_ref ) ) )**2
     ethr = starting_diag_threshold * MAX( 1.D0 / upscale, MIN( 1.D0, &
          ABS( xnew - x ) / dtau_ref ) )**2
  ELSE
     ethr = tr2 / nelec
     tr2 = starting_scf_threshold * &
          MAX( 1.D0/upscale, ABS( xnew - x ) / dtau_ref )
  END IF
  !
  ! ... report
  !
  IF ( conv_ions ) THEN
     WRITE( stdout, '(/5X,"BFGS: convergence achieved, Efinal=",F15.8)') etot
     WRITE( stdout, '(/72("-")//5X,"Final estimate of positions")')
  ELSE
     WRITE( stdout, '(/72("-")//5X, &
          &"Search of equilibrium positions: iteration # ",I4, &
          &", scf threshold ",1PE8.2/)') istep, tr2
  END IF
  !
  CALL output_tau (conv_ions)
  !
  ! ... save all quantities needed at the following iterations
  !
  IF ( .NOT. conv_ions ) THEN
     !
     CALL seqopn( iunit, TRIM( prefix )//'.bfgs', 'UNFORMATTED', exst )
     WRITE( iunit ) &
          minimum_ok, xnew, starting_scf_threshold, starting_diag_threshold
     WRITE( iunit ) dtau
     WRITE( iunit ) hessm1
     WRITE( iunit ) xold, eold, deold, oldforce
     CLOSE( UNIT = iunit, STATUS = 'KEEP' )
     !
     ! ... at next iteration read from file
     !
     restart_bfgs = .TRUE.
     !
  END IF
  !
  DEALLOCATE( hessm1, dtau, oldforce )
  !
  !
  RETURN
  !
END SUBROUTINE bfgs

