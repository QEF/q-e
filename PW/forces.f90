!
! Copyright (C) 2001-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE forces()
  !----------------------------------------------------------------------------
  !
  ! ... This routine is a driver routine which computes the forces
  ! ... acting on the atoms. The complete expression of the forces
  ! ... contains four parts which are computed by different routines:
  !
  ! ...  a)  force_lc,    local contribution to the forces
  ! ...  b)  force_cc,    contribution due to NLCC
  ! ...  c)  force_ew,    contribution due to the electrostatic ewald term
  ! ...  d)  force_us,    contribution due to the non-local potential
  ! ...  e)  force_corr,  correction term for incomplete self-consistency
  ! ...  f)  force_hub,   contribution due to the Hubbard term
  !
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE cell_base,     ONLY : at, bg, alat, omega  
  USE ions_base,     ONLY : nat, nsp, ityp, tau, zv, amass
  USE gvect,         ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                            nrxx, ngl, nl, igtongl, g, gg, gcutm
  USE lsda_mod,      ONLY : nspin
  USE symme,         ONLY : s, nsym, irt
  USE vlocal,        ONLY : strf, vloc
  USE wvfct,         ONLY : gamma_only
  USE force_mod,     ONLY : force, lforce
  USE scf,           ONLY : rho
  USE ions_base,     ONLY : if_pos
  USE ldaU,          ONLY : lda_plus_u
  USE extfield,      ONLY : tefield, forcefield
  USE control_flags, ONLY : remove_rigid_rot
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: forcenl(:,:), &
                           forcelc(:,:), &
                           forcecc(:,:), &
                           forceion(:,:), &
                           forcescc(:,:), &
                           forceh(:,:)
    ! nonlocal, local, core-correction, ewald, scf correction terms, and hubbard
  REAL(DP) :: sum, sumscf
  INTEGER  :: ipol, na
    ! counter on polarization
    ! counter on atoms
  !
  !
  CALL start_clock( 'forces' )
  !
  ALLOCATE( forcenl( 3, nat ), forcelc( 3, nat ), forcecc( 3, nat ), &
            forceh( 3, nat ), forceion( 3, nat ), forcescc( 3, nat ) )
  !    
  forcescc(:,:) = 0.D0
  forceh(:,:)   = 0.D0
  !
  WRITE( stdout, '(/,5x,"Forces acting on atoms (Ry/au):", / )')
  !
  ! ... The nonlocal contribution is computed here
  !
  CALL force_us( forcenl )
  !
  ! ... The local contribution
  !
  CALL force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, &
                 nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, g, rho, nl, &
                 nspin, gstart, gamma_only, vloc, forcelc )
  !
  ! ... The NLCC contribution
  !
  CALL force_cc( forcecc )
  !
  ! ... The Hubbard contribution
  !
  IF ( lda_plus_u ) CALL force_hub( forceh )
  !
  ! ... The ionic contribution is computed here
  !
  CALL force_ew( alat, nat, nsp, ityp, zv, at, bg, tau, omega, g, &
                 gg, ngm, gstart, gamma_only, gcutm, strf, forceion )
  !
  ! ... The SCF contribution
  !
  CALL force_corr( forcescc )
  !
  ! ... here we sum all the contributions and compute the total force acting
  ! ... on the crstal
  !
  DO ipol = 1, 3
     !
     sum = 0.D0
     !
     DO na = 1, nat
        !
        force(ipol,na) = forcenl(ipol,na)  + &
                         forceion(ipol,na) + &
                         forcelc(ipol,na)  + &
                         forcecc(ipol,na)  + &
                         forceh(ipol,na)   + &
                         forcescc(ipol,na)
        !
        IF ( tefield ) force(ipol,na) = force(ipol,na) + forcefield(ipol,na)
        !
        sum = sum + force(ipol,na)
        !
     END DO
     !
     ! ... impose total force = 0
     !
     DO na = 1, nat
        !
        force(ipol,na) = force(ipol,na) - sum / DBLE( nat )
        !
     END DO
     !
  END DO
  !
  ! ... resymmetrize (should not be needed, but ...)
  !
  IF ( nsym >= 1 ) THEN
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
  IF ( remove_rigid_rot ) &
     CALL remove_tot_torque( nat, tau, amass(ityp(:)), force  )
  !
  ! ... write on output the forces
  !
  DO na = 1, nat
     !
     WRITE( stdout, 9035) na, ityp(na), force(:,na)
     !
  END DO
  !
  ! ... forces on fixed coordinates are set to zero ( C.S. 15/10/2003 )
  !
  force(:,:)    = force(:,:)    * DBLE( if_pos )
  forcescc(:,:) = forcescc(:,:) * DBLE( if_pos )
  !
#if defined (DEBUG)
  !
  WRITE( stdout, '(5x,"The non-local contrib.  to forces")')
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), ( forcenl(ipol,na), ipol = 1, 3 )
  END DO
  WRITE( stdout, '(5x,"The ionic contribution  to forces")')
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), ( forceion(ipol,na), ipol = 1, 3 )
  END DO
  WRITE( stdout, '(5x,"The local contribution  to forces")')
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), ( forcelc(ipol,na), ipol = 1, 3 )
  END DO
  WRITE( stdout, '(5x,"The Hubbard contrib.    to forces")')
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), ( forceh(ipol,na), ipol = 1, 3 )
  END DO
  WRITE( stdout, '(5x,"The SCF correction term to forces")')
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), ( forcescc(ipol,na), ipol = 1, 3 )
  END DO
  !
#endif
  !
  sum    = 0.D0
  sumscf = 0.D0
  !
  DO na = 1, nat
     !
     sum = sum + &
           force(1,na)**2 + force(2,na)**2 + force(3,na)**2
     !
     sumscf = sumscf + &
              forcescc(1,na)**2 + forcescc(2,na)**2+ forcescc(3,na)**2
     !
  END DO
  !
  WRITE( stdout, '(/5x,"Total force = ",F12.6,5X, &
              &  "Total SCF correction = ",F12.6)') SQRT( sum ), SQRT( sumscf )
  !
  DEALLOCATE( forcenl, forcelc, forcecc, forceh, forceion, forcescc )
  !
  lforce = .TRUE.
  !
  CALL stop_clock( 'forces' )
  !
  RETURN
  !
9035 FORMAT(5X,'atom ',I3,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE forces
