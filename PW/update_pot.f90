!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#define ONE  (1.D0,0.D0)
#define ZERO (0.D0,0.D0)
!
!----------------------------------------------------------------------------
SUBROUTINE update_pot()
  !----------------------------------------------------------------------------
  !
  ! ... update the potential extrapolating the charge density and extrapolates
  ! ... the wave-functions
  !
  ! ... charge density extrapolation :
  !
  ! ... pot_order = 0   copy the old potential (nothing is done)
  !
  ! ... pot_order = 1   subtract old atomic charge density and sum the new
  ! ...                 if dynamics is done the routine extrapolates also
  ! ...                 the difference between the the scf charge and the
  ! ...                 atomic one,
  !
  ! ... pot_order = 2   first order extrapolation :
  !
  ! ...                   rho(t+dt) = 2*rho(t) - rho(t-dt)
  !
  ! ... pot_order = 3   second order extrapolation :
  !
  ! ...                   rho(t+dt) = rho(t) +
  ! ...                               + alpha0*( rho(t) - rho(t-dt) )
  ! ...                               + beta0* ( rho(t-dt) - rho(t-2*dt) )
  !
  !
  ! ... wave-functions extrapolation :
  !
  ! ... wfc_order = 0   nothing is done
  !  
  ! ... wfc_order = 2   first order extrapolation :
  !
  ! ...                   |psi(t+dt)> = 2*|psi(t)> - |psi(t-dt)>
  !
  ! ... wfc_order = 3   second order extrapolation :
  !
  ! ...                   |psi(t+dt)> = |psi(t)> +
  ! ...                               + alpha0*( |psi(t)> - |psi(t-dt)> )
  ! ...                               + beta0* ( |psi(t-dt)> - |psi(t-2*dt)> )
  !
  !
  ! ...  alpha0 and beta0 are calculated in "find_alpha_and_beta()" so that 
  ! ...  |tau'-tau(t+dt)| is minimum; 
  ! ...  tau' and tau(t+dt) are respectively the atomic positions at time 
  ! ...  t+dt and the extrapolated one:
  !
  ! ...  tau(t+dt) = tau(t) + alpha0*( tau(t) - tau(t-dt) )
  ! ...                     + beta0*( tau(t-dt) -tau(t-2*dt) )
  !
  !
  USE control_flags, ONLY : pot_order, wfc_order, history
  USE io_files,      ONLY : prefix, tmp_dir, nd_nmbr
  USE io_global,     ONLY : ionode, ionode_id
  USE mp,            ONLY : mp_bcast
  USE mp_global,     ONLY : intra_image_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: rho_extr, wfc_extr
  LOGICAL :: exists
  !
  !
  CALL start_clock( 'update_pot' )
  !
  IF ( pot_order == 0 .AND. wfc_order == 0 ) THEN
     !
     CALL stop_clock( 'update_pot' )
     !
     RETURN
     !
  END IF 
  !
  ! ... determines the maximum effective order of the extrapolation on the 
  ! ... basis of the files that are really available
  !
  ! ... for the charge density
  !
  IF ( ionode ) THEN
     !
     rho_extr = MIN( 1, history, pot_order )
     !
     INQUIRE( FILE = TRIM( tmp_dir ) // &
            & TRIM( prefix ) // '.oldrho', EXIST = exists )
     !
     IF ( exists ) THEN
        !
        rho_extr = MIN( 2, history, pot_order )
        !
        INQUIRE( FILE = TRIM( tmp_dir ) // &
               & TRIM( prefix ) // '.old2rho', EXIST = exists )     
        !
        IF ( exists ) rho_extr = MIN( 3, history, pot_order )
        !
     END IF
     !
  END IF
  !
  CALL mp_bcast( rho_extr, ionode_id, intra_image_comm )
  !
  IF ( pot_order > 0 ) CALL extrapolate_charge( rho_extr )
  !
  ! ... for the wavefunctions
  !
  IF ( ionode ) THEN
     !
     wfc_extr = MIN( 1, history, wfc_order ) 
     !
     INQUIRE( FILE = TRIM( tmp_dir ) // &
            & TRIM( prefix ) // '.oldwfc' // nd_nmbr, EXIST = exists )
     !
     IF ( exists ) THEN
        !
        wfc_extr = MIN( 2, history, wfc_order  )
        !
        INQUIRE( FILE = TRIM( tmp_dir ) // &
               & TRIM( prefix ) // '.old2wfc' // nd_nmbr , EXIST = exists )     
        !
        IF ( exists ) wfc_extr = MIN( 3, history, wfc_order )
        !
     END IF
     !
  END IF
  !
  CALL mp_bcast( wfc_extr, ionode_id, intra_image_comm )
  !
  IF ( wfc_order > 0 ) CALL extrapolate_wfcs( wfc_extr )
  !
  CALL stop_clock( 'update_pot' )
  !
  RETURN
  !
END SUBROUTINE update_pot
!
!----------------------------------------------------------------------------
SUBROUTINE extrapolate_charge( rho_extr )
  !----------------------------------------------------------------------------
  !
  USE constants,        ONLY : eps32
  USE io_global,        ONLY : stdout
  USE kinds,            ONLY : DP
  USE cell_base,        ONLY : omega, bg, alat
  USE ions_base,        ONLY : nat, tau, nsp, ityp
  USE gvect,            ONLY : nrxx, ngm, g, gg, gstart,  nr1, nr2, nr3, nl, &
                               eigts1, eigts2, eigts3, nrx1, nrx2, nrx3
  USE lsda_mod,         ONLY : lsda, nspin
  USE scf,              ONLY : rho, rho_core, vr
  USE control_flags,    ONLY : alpha0, beta0
  USE ener,             ONLY : ehart, etxc, vtxc
  USE extfield,         ONLY : etotefield
  USE cellmd,           ONLY : lmovecell, omega_old
  USE vlocal,           ONLY : strf
  USE noncollin_module, ONLY : noncolin
  USE noncollin_module, ONLY : factlist, pointlist, pointnum, mcons,&
                               i_cons, lambda, vtcon, report
  USE io_files,         ONLY : prefix
  USE klist,            ONLY : nelec
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: rho_extr
  !
  REAL(DP), ALLOCATABLE :: work(:), work1(:)
    ! work is the difference between charge density and atomic charge 
    !   at time t
    ! work1 is the same thing at time t-dt
  REAL(DP) :: charge
  !
  INTEGER :: ir, is
  !
  !
  IF ( rho_extr == 0 ) THEN
     !
     ! ... calculate structure factors for the new positions
     !
     IF ( lmovecell ) CALL scale_h()
     !
     CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, &
                      nr1, nr2, nr3, strf, eigts1, eigts2, eigts3 )
     !
     RETURN
     !
  END IF
  !
  ALLOCATE( work(nrxx) )
  !
  work(:) = 0.D0
  !
  ! ... in the lsda case the magnetization will follow rigidly the density
  ! ... keeping fixed the value of zeta = mag / rho_tot. 
  ! ... zeta is set here and put in rho(:,2) while rho(:,1) will contain the 
  ! ... total valence charge
  !
  IF ( lsda ) CALL rho2zeta( rho, rho_core, nrxx, nspin, 1 )
  !
  IF ( noncolin ) THEN
     !
     DO is = 2, nspin
        !
        WHERE( rho(:,1) > eps32 )
           !
           rho(:,is) = rho(:,is) / rho(:,1)
           !
        ELSEWHERE
           !
           rho(:,is) = 0.D0
           !
        END WHERE
        !
     END DO
     !
  END IF
  !
  ! ... subtract the old atomic charge density
  !
  CALL atomic_rho( work, 1 )
  !
  rho(:,1) = rho(:,1) - work(:)
  !
  IF ( lmovecell ) rho(:,1) = rho(:,1) * omega_old
  !
  ! ... extrapolate the difference between the atomic charge and
  ! ... the self-consistent one
  !
  IF ( rho_extr == 1 ) THEN
     !
     ! ... if rho_extr = 1  update the potential subtracting to the charge 
     ! ...                  density the "old" atomic charge and summing the 
     ! ...                  new one
     !
     WRITE( stdout, &
            '(/5X,"NEW-OLD atomic charge density approx. for the potential")' )
     !
     CALL io_pot( + 1, 'oldrho', rho, 1 )
     !
  ELSE IF ( rho_extr == 2 ) THEN
     !
     WRITE( UNIT = stdout, &
            FMT = '(/5X,"first order charge density extrapolation")' )
     !
     ! ...   oldrho  ->  work
     !
     CALL io_pot( - 1, 'oldrho',  work, 1 )
     !
     ! ...   rho   ->  oldrho          
     ! ...   work  ->  oldrho2     
     !
     CALL io_pot( + 1, 'oldrho',  rho,  1 )
     CALL io_pot( + 1, 'old2rho', work, 1 )
     !
     ! ... extrapolation
     !
     rho(:,1) = 2.D0 * rho(:,1) - work(:)
     !
  ELSE IF ( rho_extr == 3 ) THEN  
     !
     WRITE( UNIT = stdout, &
            FMT = '(/5X,"second order charge density extrapolation")' )
     !
     ALLOCATE( work1(nrxx) )
     !
     work1(:) = 0.D0
     !
     ! ...   oldrho2  ->  work1
     ! ...   oldrho   ->  work
     !
     CALL io_pot( - 1, 'old2rho', work1, 1 )
     CALL io_pot( - 1, 'oldrho',  work,  1 )
     !
     ! ...   rho   ->  oldrho     
     ! ...   work  ->  oldrho2     
     !
     CALL io_pot( + 1, 'oldrho',  rho,  1 )
     CALL io_pot( + 1, 'old2rho', work, 1 )
     !
     ! ... alpha0 and beta0 have been calculated in move_ions
     !
     rho(:,1) = rho(:,1) + alpha0 * ( rho(:,1) - work(:) ) + &
                            beta0 * ( work(:) - work1(:) )
     !
     DEALLOCATE( work1 )
     !
  END IF
  !
  IF ( lmovecell ) rho(:,1) = rho(:,1) / omega
  !
  ! ... calculate structure factors for the new positions
  !
  IF ( lmovecell ) CALL scale_h()
  !
  CALL struc_fact( nat, tau, nsp, ityp, ngm, g, bg, &
                   nr1, nr2, nr3, strf, eigts1, eigts2, eigts3 )
  !
  ! ... add atomic charges in the new positions
  !
  CALL atomic_rho( work, 1 )
  !
  rho(:,1) = rho(:,1) + work(:)
  !
  CALL set_rhoc()
  !
  ! ... reset up and down charge densities in the LSDA case
  !
  IF ( lsda ) CALL rho2zeta( rho, rho_core, nrxx, nspin, -1 )
  !
  IF ( noncolin ) THEN
     !
     DO is = 2, nspin
        !
        WHERE( rho(:,1) > eps32 )
           !
           rho(:,is) = rho(:,is) * rho(:,1)
           !
        ELSEWHERE
           !
           rho(:,is) = 0.D0
           !
        END WHERE
        !
     END DO
     !
  END IF
  !
  CALL v_of_rho( rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                 nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
                 ehart, etxc, vtxc, etotefield, charge, vr )
  !
  IF ( ABS( charge - nelec ) / charge > 1.D-7 ) THEN
     !
     WRITE( stdout, &
            '(/,5X,"extrapolated charge ",F10.5,", renormalised to ",F10.5)') &
         charge, nelec
     !
     rho = rho / charge * nelec
     !
  END IF
  !
  DEALLOCATE( work )
  !
  RETURN
  !
END SUBROUTINE extrapolate_charge
!
!-----------------------------------------------------------------------
SUBROUTINE extrapolate_wfcs( wfc_extr )
  !-----------------------------------------------------------------------
  !
  ! ... This routine extrapolate the wfc's after a "parallel alignment"
  ! ... of the basis of the t-dt and t time steps, according to a recipe
  ! ... by Mead, Rev. Mod. Phys., vol 64, pag. 51 (1992), eqs. 3.20-3.29
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE klist,                ONLY : nks
  USE control_flags,        ONLY : isolve, alpha0, beta0, wfc_order
  USE wvfct,                ONLY : nbnd, npw, npwx, igk
  USE io_files,             ONLY : nwordwfc, iunigk, iunwfc, iunoldwfc, &
                                   iunoldwfc2, prefix
  USE noncollin_module,     ONLY : noncolin
  USE wavefunctions_module, ONLY : evc
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: wfc_extr
  !
  INTEGER :: j, i, ik, zero_ew, lwork, info
    ! do-loop variables
    ! counter on k-points
    ! number of zero 'eigenvalues' of the s_m matrix
    ! used by singular value decomposition (ZGESVD)
    ! flag returned by ZGESVD
  COMPLEX(DP), ALLOCATABLE :: s_m(:,:), sp_m(:,:), u_m(:,:), w_m(:,:), work(:)
    ! the overlap matrix s (eq. 3.24)
    ! its dagger
    ! left unitary matrix in the SVD of sp_m
    ! right unitary matrix in the SVD of sp_m
    ! workspace for ZGESVD
  COMPLEX(DP), ALLOCATABLE :: evcold(:,:)
    ! wavefunctions at previous iteration
  REAL(DP), ALLOCATABLE :: ew(:), rwork(:)
    ! the eigenvalues of s_m
    ! workspace for ZGESVD
  LOGICAL :: exst
  !
  !
  IF ( wfc_extr == 1 ) THEN
     !
     CALL diropn( iunoldwfc, 'oldwfc', nwordwfc, exst )
     !
     DO ik = 1, nks
        !
        ! ... "now"  -> "old"
        !
        CALL davcio( evc, nwordwfc, iunwfc,    ik, - 1 )
        CALL davcio( evc, nwordwfc, iunoldwfc, ik, + 1 )
        !
     END DO
     !
     CLOSE( UNIT = iunoldwfc, STATUS = 'KEEP' )
     !
  ELSE IF ( wfc_extr == 2 ) THEN
     !
     CALL diropn( iunoldwfc, 'oldwfc', nwordwfc, exst )
     !
     IF ( wfc_order > 2 ) &
        CALL diropn( iunoldwfc2, 'old2wfc', nwordwfc, exst )
     !
     ALLOCATE( evcold(npwx,nbnd) )
     !
     WRITE( UNIT = stdout, &
            FMT = '(5X,"first order wave-functions extrapolation")' )
     !
     lwork = 5 * nbnd
     !
     ALLOCATE( s_m(nbnd,nbnd), sp_m(nbnd,nbnd), u_m(nbnd,nbnd), &
               w_m(nbnd,nbnd), work(lwork), ew(nbnd), rwork(lwork) )
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
        CALL davcio( evc,    nwordwfc, iunwfc,    ik, - 1 )
        !
        ! ... construct s_m = <evcold|evc>
        !
        CALL ZGEMM( 'C', 'N', nbnd, nbnd, npw, ONE, evcold, npwx, evc, &
                    npwx, ZERO, s_m, nbnd )
        !            
        CALL reduce( 2 * nbnd * nbnd, s_m )
        !
        ! ... construct sp_m
        !
        DO i = 1, nbnd
          ! 
          sp_m(:,i) = CONJG( s_m (i,:) )
          !
        END DO
        !
        ! ... the unitary matrix [sp_m*s_m]^(-1/2)*sp_m (eq. 3.29)
        ! ... by means the singular value decomposition (SVD) of
        ! ... sp_m = u_m * diag(ew) * w_m 
        ! ... becomes u_m * w_m
        !
        CALL ZGESVD( 'A', 'A', nbnd, nbnd, sp_m, nbnd, ew, u_m, &
                     nbnd, w_m, nbnd, work, lwork, rwork, info )
        !
        ! ... check on eigenvalues
        !
        zero_ew = COUNT( ew(:) < 0.1D0 )
        !
        ! ... use sp_m to store u_m * w_m
        !
        CALL ZGEMM( 'N', 'N', nbnd, nbnd, nbnd, ONE, &
                    u_m, nbnd, w_m, nbnd, ZERO, sp_m, nbnd )
        !
        ! ... now use evcold as workspace to calculate "aligned" wavefcts:
        !
        ! ... evcold_i = sum_j evc_j*sp_m_ji (eq.3.21)
        !
        CALL ZGEMM( 'N', 'N', npw, nbnd, nbnd, ONE, &
                    evc, npwx, sp_m, nbnd, ZERO, evcold, npwx )
        !
        ! ... save on file the aligned wavefcts
        !
        CALL davcio( evcold, nwordwfc, iunwfc, ik, + 1 )
        !
        ! ... re-read from file the wavefcts at (t-dt)
        !
        CALL davcio( evc, nwordwfc, iunoldwfc, ik, - 1 )
        !
        ! ... extrapolate the wfc's (note that evcold contains wavefcts 
        ! ... at (t) and evc contains wavefcts at (t-dt) )
        !
        evc = 2.D0 * evcold - evc
        !
        ! ... move the files: "old" -> "old1" and "now" -> "old"
        !
        IF ( wfc_order > 2 ) THEN
           !
           CALL davcio( evcold, nwordwfc, iunoldwfc,  ik, - 1 )
           CALL davcio( evcold, nwordwfc, iunoldwfc2, ik, + 1 )
           !
        END IF
        !
        CALL davcio( evcold, nwordwfc, iunwfc,    ik, - 1 )
        CALL davcio( evcold, nwordwfc, iunoldwfc, ik, + 1 )
        !
        ! ... save evc on file iunwfc
        !
        CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
        !
     END DO
     !
     IF ( zero_ew > 0 ) &
        WRITE( stdout, '(/,5X,"Message from extrapolate_wfcs: ",/,     &
                        &  5X,"the matrix <psi(t-dt)|psi(t)> has ",I2, &
                        &     " zero eigenvalues")' ) zero_ew     
     !
     DEALLOCATE( s_m, sp_m, u_m, w_m, work, ew, rwork )
     !
     DEALLOCATE( evcold )
     !
     CLOSE( UNIT = iunoldwfc, STATUS = 'KEEP' )
     IF ( wfc_order > 2 ) &
        CLOSE( UNIT = iunoldwfc2, STATUS = 'KEEP' )     
     !
  ELSE
     !
     ! ... case :  wfc_extr = 3
     !
     CALL diropn( iunoldwfc,  'oldwfc',  nwordwfc, exst )
     CALL diropn( iunoldwfc2, 'old2wfc', nwordwfc, exst )
     !
     ALLOCATE( evcold(npwx,nbnd) )
     !
     WRITE( UNIT = stdout, &
            FMT = '(5X,"second order wave-functions extrapolation")' )
     !
     lwork = 5 * nbnd
     !
     ALLOCATE( s_m(nbnd,nbnd), sp_m(nbnd,nbnd), u_m(nbnd,nbnd), &
               w_m(nbnd,nbnd), work(lwork), ew(nbnd), rwork(lwork) )
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
        CALL davcio( evc,    nwordwfc, iunwfc,    ik, - 1 )
        !
        ! ... construct s_m = <evcold|evc>
        !
        CALL ZGEMM( 'C', 'N', nbnd, nbnd, npw, ONE, evcold, npwx, evc, &
                    npwx, ZERO, s_m, nbnd )
        !            
        CALL reduce( 2 * nbnd * nbnd, s_m )
        !
        ! ... construct sp_m
        !
        DO i = 1, nbnd
          ! 
          sp_m(:,i) = CONJG( s_m (i,:) )
          !
        END DO
        !
        ! ... the unitary matrix [sp_m*s_m]^(-1/2)*sp_m (eq. 3.29)
        ! ... by means the singular value decomposition (SVD) of
        ! ... sp_m = u_m * diag(ew) * w_m 
        ! ... becomes u_m * w_m
        !
        CALL ZGESVD( 'A', 'A', nbnd, nbnd, sp_m, nbnd, ew, u_m, &
                     nbnd, w_m, nbnd, work, lwork, rwork, info )
        !
        ! ... check on eigenvalues
        !
        zero_ew = COUNT( ew(:) < 0.1D0 )
        !
        ! ... use sp_m to store u_m * w_m
        !
        CALL ZGEMM( 'N', 'N', nbnd, nbnd, nbnd, ONE, &
                    u_m, nbnd, w_m, nbnd, ZERO, sp_m, nbnd )
        !
        ! ... now use evcold as workspace to calculate "aligned" wavefcts:
        !
        ! ... evcold_i = sum_j evc_j*sp_m_ji (eq.3.21)
        !
        CALL ZGEMM( 'N', 'N', npw, nbnd, nbnd, ONE, &
                    evc, npwx, sp_m, nbnd, ZERO, evcold, npwx )
        !
        ! ... save on file the aligned wavefcts
        !
        CALL davcio( evcold, nwordwfc, iunwfc, ik, + 1 )
        !
        ! ... re-read from file the wavefcts at (t-dt)
        !
        CALL davcio( evc, nwordwfc, iunoldwfc, ik, - 1 )
        !
        ! ... extrapolate the wfc's,
        ! ... if wfc_extr == 3 use the second order extrapolation formula
        ! ... alpha0 and beta0 are calculated in "move_ions"
        !
        evc = ( 1 + alpha0 ) * evcold + ( beta0 - alpha0 ) * evc
        !
        CALL davcio( evcold, nwordwfc, iunoldwfc2, ik, - 1 )
        !
        evc = evc - beta0 * evcold 
        !
        ! ... move the files: "old" -> "old1" and "now" -> "old"
        !
        CALL davcio( evcold, nwordwfc, iunoldwfc,  ik, - 1 )
        CALL davcio( evcold, nwordwfc, iunoldwfc2, ik, + 1 )
        CALL davcio( evcold, nwordwfc, iunwfc,     ik, - 1 )
        CALL davcio( evcold, nwordwfc, iunoldwfc,  ik, + 1 )
        !
        ! ... save evc on file iunwfc
        !
        CALL davcio( evc, nwordwfc, iunwfc, ik, 1 )
        !
     END DO
     !
     IF ( zero_ew > 0 ) &
        WRITE( stdout, '(/,5X,"Message from extrapolate_wfcs: ",/,     &
                        &  5X,"the matrix <psi(t-dt)|psi(t)> has ",I2, &
                        &     " zero eigenvalues")' ) zero_ew     
     !
     DEALLOCATE( s_m, sp_m, u_m, w_m, work, ew, rwork )
     !
     DEALLOCATE( evcold )
     !
     CLOSE( UNIT = iunoldwfc,  STATUS = 'KEEP' )
     CLOSE( UNIT = iunoldwfc2, STATUS = 'KEEP' )
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE extrapolate_wfcs
