!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!----------------------------------------------------------------------------
SUBROUTINE update_pot()
  !----------------------------------------------------------------------------
  !
  ! ... update potential, use the integer variable order to decide the way
  !
  !     order = 0       copy the old potential (nothing is done)
  !
  !     order = 1       subtract old atomic charge density and sum the new
  !                     if dynamics is done the routine extrapolates also
  !                     the difference between the the scf charge and the
  !                     atomic one,
  !
  !     order = 2       extrapolate the wavefunctions:
  !                       |psi(t+dt)> = 2*|psi(t)> - |psi(t-dt)>
  !
  !     order = 3       extrapolate the wavefunctions with the second-order
  !                     formula:
  !                       |psi(t+dt)> = |psi(t) +
  !                                   + alpha0*(|psi(t)> -    |psi(t-dt)>
  !                                   + beta0* (|psi(t-dt)> - |psi(t-2*dt)>
  !
  !                     where alpha0 and beta0 are calculated in "dynamics" so
  !                     that |tau'-tau(t+dt)| is minimum; tau' and tau(t+dt)
  !                     are respectively the atomic positions at time t+dt
  !                     and  the extrapolated one:
  !                       tau(t+dt) = tau(t) +
  !                                    + alpha0*( tau(t) - tau(t-dt) )
  !                                    + beta0*( tau(t-dt) -tau(t-2*dt) )
  !
  !
  USE control_flags, ONLY : lbfgs, lneb, order
  !
  IMPLICIT NONE
  !
  !
  CALL start_clock( 'update_pot' )
  !
  IF ( order == 0 ) RETURN
  !
  IF ( order > 2 .AND. ( lbfgs .OR. lneb ) ) THEN
     !
     order = 2
     CALL errore( 'update_pot', 'order > 2 not allowed in bfgs', -1 )
     !
  END IF
  !
  CALL extrapolate_charge()
  !
  IF ( order >= 2 ) CALL extrapolate_wfcs()
  !
  CALL stop_clock( 'update_pot' )
  !
  RETURN
  !
END SUBROUTINE update_pot
!
!
!----------------------------------------------------------------------------
SUBROUTINE extrapolate_charge()
  !----------------------------------------------------------------------------
  !
  USE io_global,   ONLY :  stdout
  USE kinds,  ONLY :  DP
  USE brilz,       ONLY :  omega, bg, alat
  USE basis,       ONLY :  nat, tau, ntyp, ityp
  USE gvect,       ONLY :  nrxx, ngm, g, gg, gstart,  nr1, nr2, nr3, nl, &
                           eigts1, eigts2, eigts3, nrx1, nrx2, nrx3
  USE lsda_mod,    ONLY :  lsda, nspin
  USE scf,         ONLY :  rho, rho_core, vr
  USE control_flags,       ONLY :  lbfgs, istep, alpha0, beta0, imix
  USE ener,        ONLY :  ehart, etxc, vtxc
  USE cellmd,      ONLY :  lmovecell, omega_old
  USE vlocal,      ONLY :  strf
  USE io_files,    ONLY :  prefix
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP), ALLOCATABLE :: work(:), work1(:)
    ! work is the difference between charge density and atomic charge 
    !   at time t
    ! work1 is the same thing at time t-dt
  REAL(KIND=DP) :: charge
  !
  !
  IF ( istep == 0 ) RETURN 
  !
  ALLOCATE( work(nrxx) )
  work(:) = 0.D0
  !
  ! ... if order = 1 update the potential subtracting to the charge density
  ! ... the "old" atomic charge and summing the new one
  !
  WRITE( stdout,'(/5X,"NEW-OLD atomic charge density approx. for the potential")' )
  !
  ! ... in the lsda case the magnetization will follow rigidly the density
  ! ... keeping fixed the value of zeta = mag / rho_tot. 
  ! ... zeta is set here and put in rho(* ??? while rho(*,1) will contain the 
  ! ... total valence charge
  !
  IF ( lsda ) CALL rho2zeta( rho, rho_core, nrxx, nspin, 1 )
  !
  ! ... subtract the old atomic charge density
  !
  CALL atomic_rho( work, 1 )
  !
  CALL DAXPY( nrxx, -1.D0, work, 1, rho, 1 )
  !
  IF ( lmovecell ) CALL DSCAL( nrxx, omega_old, rho, 1 )
  !
  ! ... if dynamics extrapolate the difference between the atomic charge a
  ! ... the self-consistent one
  !
  IF ( .NOT. lbfgs ) THEN
     ! 
     IF ( istep == 1 ) THEN
        !
        CALL io_pot( + 1, TRIM( prefix )//'.oldrho', rho, 1 )
        !
     ELSE
        !
        ALLOCATE( work1(nrxx) )
        !
        work1(:) = 0.D0
        !
        CALL io_pot( - 1, TRIM( prefix )//'.oldrho', work, 1 )
        CALL io_pot( + 1, TRIM( prefix )//'.oldrho', rho, 1 )
        !
        IF ( istep == 2 ) &
           CALL io_pot( + 1, TRIM( prefix )//'.oldrho2', work, 1 )
        !
        CALL io_pot( - 1, TRIM( prefix )//'.oldrho2', work1, 1 )
        CALL io_pot( + 1, TRIM( prefix )//'.oldrho2', work, 1 )
        !
        ! ... alpha0 and beta0 have been calculated in dynamics 
        ! ... or in vcsmd subs.
        !
        rho(:,1) = rho(:,1) + alpha0 * ( rho(:,1) - work(:) ) + &
                               beta0 * ( work(:) - work1(:) )
        !
        DEALLOCATE( work1 )
        !
     END IF
     !
  END IF
  !
  IF ( lmovecell ) CALL DSCAL( nrxx, 1.0D0 / omega, rho, 1 )
  !
  ! ... calculate structure factors for the new positions
  !
  IF ( lmovecell ) CALL scale_h()
  !
  CALL struc_fact( nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, nr3, &
                   strf, eigts1, eigts2, eigts3 )
  !
  ! ... add atomic charges in the new positions
  !
  CALL atomic_rho( work, 1 )
  CALL DAXPY( nrxx, 1.D0, work, 1, rho, 1 )
  CALL set_rhoc()
  !
  ! ... reset up and down charge densities in the LSDA case
  !
  IF ( lsda ) CALL rho2zeta( rho, rho_core, nrxx, nspin, -1 )
  !
  CALL v_of_rho( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                 nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
                 ehart, etxc, vtxc, charge, vr )
  !
  ! ... write potential (and rho) on file
  !
  IF ( imix >= 0 ) CALL io_pot( + 1, TRIM( prefix )//'.rho', rho, nspin )
  !
  CALL io_pot( + 1, TRIM( prefix )//'.pot', vr, nspin )
  !
  DEALLOCATE( work )
  !
  RETURN
  !
END SUBROUTINE extrapolate_charge
!
!
!-----------------------------------------------------------------------
SUBROUTINE extrapolate_wfcs()
  !-----------------------------------------------------------------------
  !
  ! ... This routine extrapolate the wfc's after a "parallel alignment"
  ! ... of the basis of the t-dt and t time steps, according to the Mead
  ! ... recipe, see Rev. Mod. Phys., vol 64, pag. 51 (1992), eqs. 3.20-3.2
  !
#define ONE (1.D0,0.D0)
#define ZERO (0.D0,0.D0)  
  !
  USE io_global,            ONLY : stdout
  USE kinds,           ONLY : DP
  USE klist,                ONLY : nks
  USE control_flags,                ONLY : isolve, istep, order, alpha0, beta0
  USE basis,                ONLY : startingwfc
  USE wvfct,                ONLY : nbnd, npw, npwx, igk
  USE io_files,             ONLY : nwordwfc, iunigk, iunwfc, iunoldwfc, &
                                   iunoldwfc2
  USE wavefunctions_module, ONLY : evc
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: j, i, ik, zero_ew
    ! do-loop variables
    ! counter on k-points
    ! number of zero eigenvalues of the sp_m * s_m matrix
  COMPLEX(KIND=DP), ALLOCATABLE :: u_m(:,:), s_m(:,:), sp_m(:,:), temp(:,:)
    ! the unitary matrix (eq. 3.21)
    ! the overlap matrix s (eq. 3.24)
    ! its dagger
    ! workspace
  COMPLEX(KIND=DP), ALLOCATABLE :: evcold(:,:)
    ! wavefunctions at previous iteration
  REAL(KIND=DP), ALLOCATABLE :: ew(:)
    ! the eigenvalues of sp_m*s_m
  LOGICAL :: first
    ! Used for initialization
  DATA first / .TRUE. /
  SAVE first
  !
  !
  ! ... istep = 0 when extrapolate_wfcs() is called in neb 
  !
  IF ( istep == 0 ) RETURN   
  !
  IF ( first ) THEN
     !
     first = .FALSE.
     !
     IF ( isolve == 1 .AND. startingwfc == 'atomic' ) THEN
        !
        DEALLOCATE( evc )
        ALLOCATE( evc(npwx,nbnd) )
        !
     END IF
     !
  END IF
  !
  ALLOCATE( evcold(npwx,nbnd) )
  !
  IF ( istep == 1 ) THEN
     !
     IF ( nks > 1 ) REWIND( iunigk )
     !
     DO ik = 1, nks
        !
        IF ( nks > 1 ) READ( iunigk ) npw, igk
        !
        CALL davcio( evc, nwordwfc, iunwfc, ik, - 1 )
        !
        CALL ZCOPY( npwx * nbnd, evc, 1, evcold, 1 )
        !
        CALL davcio( evcold, nwordwfc, iunoldwfc, ik, 1 )
        !
     END DO
     !
  ELSE
     !
     IF ( order == 2 ) THEN
        WRITE( stdout, '(5X,"Extrapolating wave-functions (first order) ...")' )
     ELSE
        WRITE( stdout, '(5X,"Extrapolating wave-functions (second order) ...")' )
     END IF
     !
     ALLOCATE( u_m(nbnd,nbnd), s_m(nbnd,nbnd), sp_m(nbnd,nbnd), &
               temp(nbnd,nbnd), ew(nbnd) )
     !
     IF ( nks > 1 ) REWIND( iunigk )
     !
     zero_ew = 0
     !
     DO ik = 1, nks
        !
        IF ( nks > 1 ) READ( iunigk ) npw, igk
        !
        CALL davcio( evcold, nwordwfc, iunoldwfc, ik, - 1 )
        CALL davcio( evc, nwordwfc, iunwfc, ik, - 1 )
        !
        IF ( istep == 2 .AND. order > 2 ) &
           CALL davcio( evcold, nwordwfc, iunoldwfc2, ik, 1 )
        !
        ! ... construct s_m = <evcold|evc>
        !
        CALL ZGEMM( 'C', 'N', nbnd, nbnd, npw, ONE, evcold, npwx, evc, &
                    npwx, ZERO, s_m, nbnd )
#ifdef __PARA
        CALL reduce( 2 * nbnd * nbnd, s_m )
#endif
        !
        ! ... temp = sp_m * s_m
        !
        CALL ZGEMM( 'C', 'N', nbnd, nbnd, nbnd, ONE, s_m, nbnd, s_m, &
                    nbnd, ZERO, temp, nbnd )
        !
        ! ... diagonalize temp, use u_m as workspace to accomodate the 
        ! ... eigenvect matrix which diagonalizes temp, sp_m is its hermitean 
        ! ... conjugate
        !
        CALL cdiagh( nbnd, temp, nbnd, ew, u_m )
        !
        ! ... check on eigenvalues
        !
        DO i = 1, nbnd
          !
          IF ( ew(i) < 0.1D0 ) zero_ew = zero_ew + 1
          !
        END DO
        !
        DO i = 1, nbnd
          ! 
          sp_m(:,i) = CONJG( u_m (i,:) ) / SQRT( ABS( ew(:) ) )
          !
        END DO
        !
        ! ... temp = u_m * sp_m
        !
        CALL ZGEMM( 'N', 'N', nbnd, nbnd, nbnd, ONE, u_m, nbnd, sp_m, &
                    nbnd, ZERO, temp, nbnd )
        !
        ! ... temp = [ sp_m * s_m ]^(-1/2)
        !
        CALL ZGEMM( 'N', 'C', nbnd, nbnd, nbnd, ONE, temp, nbnd, s_m, &
                    nbnd, ZERO, u_m, nbnd )
        !
        ! ... and u_m is the unitary matrix [ sp_m * s_m ]^(-1/2)*sp_m (eq.3.29)
        ! ... now use evcold as workspace to calculate
        !
        ! ... evcold_i = sum_j evc_j*u_m_ji
        !
        CALL ZGEMM( 'N', 'N', npw, nbnd, nbnd, ONE, evc, npwx, u_m, &
                    nbnd, ZERO, evcold, npwx )
        !
        ! ... and copy evcold in evc
        !
        CALL ZCOPY( npwx * nbnd, evcold, 1, evc, 1 )
        !
        ! ... save on file evc
        !
        CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
        !
        ! ... re-read from file the right evcold
        !
        CALL davcio( evcold, nwordwfc, iunoldwfc, ik, - 1 )
        !
        ! ... extrapolate the wfc's, if order=3 use the second order 
        ! ... extrapolati formula, alpha0 and beta0 are calculated in 
        ! ... "dynamics"
        !
        IF ( order > 2 ) THEN
           !
           evc = ( 1 + alpha0 ) * evc + ( beta0 - alpha0 ) * evcold
           !
           CALL davcio( evcold, nwordwfc, iunoldwfc2, ik, - 1 )
           !
           evc = evc - beta0 * evcold
           !
        ELSE
           !
           evc = 2 * evc - evcold
           !
        END IF
        !
        ! ... move the files: "old" -> "old1" and "now" -> "old"
        !
        IF ( order > 2 ) THEN
           CALL davcio( evcold, nwordwfc, iunoldwfc, ik, - 1 )
           CALL davcio( evcold, nwordwfc, iunoldwfc2, ik, 1 )
        END IF
        !
        CALL davcio( evcold, nwordwfc, iunwfc, ik, - 1 )
        CALL davcio( evcold, nwordwfc, iunoldwfc, ik, 1 )
        !
        ! ... save evc on file iunwfc
        !
        CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
        !
     END DO
     !
     IF ( zero_ew > 0 ) &
        WRITE( stdout, '(/,5X,"from extrapolate_wfcs : WARNING",/,     &
                        &  5X,"the matrix <psi(t-dt)|psi(t)> has ",I2, &
                        &     " zero eigenvalues",/,                   &
                        &  5X,"change the number of bands ( nbnd )")' ) zero_ew     
     !
     DEALLOCATE( u_m, s_m, sp_m, temp, ew )
     !
  END IF
  !
  DEALLOCATE( evcold )
  !
  RETURN
  !
END SUBROUTINE extrapolate_wfcs
