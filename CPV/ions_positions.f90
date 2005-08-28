!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
MODULE ions_positions
!------------------------------------------------------------------------------!
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : natx
  !
  IMPLICIT NONE
  !
  ! ... Atomic positions arrays used in the cp codes during the dynamic
  !
  REAL(DP) :: tau0(3,natx), taum(3,natx),  taup(3,natx)
  REAL(DP) :: taus(3,natx), tausm(3,natx), tausp(3,natx)
  REAL(DP) :: vels(3,natx), velsm(3,natx), velsp(3,natx)
  !
  CONTAINS 
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ions_hmove( taus, tausm, iforce, pmass, fion, ainv, delt, na, nsp )
    !--------------------------------------------------------------------------
    !
    REAL(DP), INTENT(IN)  :: tausm(:,:), pmass(:), fion(:,:)
    INTEGER,        INTENT(IN)  :: iforce(:,:)
    REAL(DP), INTENT(IN)  :: ainv(3,3), delt
    REAL(DP), INTENT(OUT) :: taus(:,:) 
    INTEGER,        INTENT(IN)  :: na(:), nsp
    INTEGER                     :: is, ia, i, isa
    REAL(DP)              :: dt2by2, fac, fions(3)
    !
    !
    dt2by2 = 0.5D0 * delt * delt
    !
    isa = 0
    !
    DO is = 1, nsp
       !
       fac = dt2by2  / pmass(is)
       !
       DO ia = 1, na(is)
          !
          isa = isa + 1
          !
          DO i = 1, 3
             !
             fions(i) = fion(1,isa) * ainv(i,1) + &
                        fion(2,isa) * ainv(i,2) + &
                        fion(3,isa) * ainv(i,3)
             !
          END DO
          !
          taus(:,isa) = tausm(:,isa) + iforce(:,isa) * fac * fions(:)
          !
       END DO
       
    END DO
    !
    RETURN
    !
  END SUBROUTINE ions_hmove
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ions_move( tausp, taus, tausm, iforce, pmass, fion, ainv, &
                        delt, na, nsp, fricp, hgamma, vels, tsdp, tnosep, &
                        fionm, vnhp, velsp, velsm, nhpcl, nhpdim, atm2nhp )
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)    :: taus(:,:), tausm(:,:), pmass(:), fion(:,:)
    INTEGER,        INTENT(IN)    :: iforce(:,:)
    REAL(DP), INTENT(IN)    :: ainv(3,3), delt
    REAL(DP), INTENT(OUT)   :: tausp(:,:)
    INTEGER,        INTENT(IN)    :: na(:), nsp, nhpcl, nhpdim, atm2nhp(:)
    REAL(DP), INTENT(IN)    :: fricp, hgamma(3,3), vels(:,:)
    LOGICAL,        INTENT(IN)    :: tsdp, tnosep
    REAL(DP), INTENT(INOUT) :: fionm(:,:)
    REAL(DP), INTENT(IN)    :: vnhp(nhpcl,nhpdim)
    REAL(DP), INTENT(OUT)   :: velsp(:,:)
    REAL(DP), INTENT(IN)    :: velsm(:,:)
    INTEGER                       :: is, ia, i, isa
    REAL(DP)                :: dt2by2, fac, dt2, twodel
    REAL(DP)                :: verl1, verl2, verl3
    !
    !
    dt2by2 = 0.5D0 * delt * delt
    dt2    = delt * delt
    twodel = 2.D0 * delt
    !
    verl1 = 2.D0 / ( 1.D0 + fricp )
    verl2 = 1.D0 - verl1
    verl3 = dt2 / ( 1.D0 + fricp )
    !
    IF ( tsdp ) THEN
       !
       isa = 0
       !
       DO is = 1, nsp
          !
          DO ia = 1, na(is)
             !
             isa = isa + 1
             !
             tausp(:,isa) = taus(:,isa) - pmass(is) * &
                                             ( hgamma(i,1) * vels(1,isa) + &
                                               hgamma(i,2) * vels(2,isa) + &
                                               hgamma(i,3) * vels(3,isa) ) + &
                            iforce(:,isa) * dt2by2 / pmass(is) * &
                                             ( fion(1,isa) * ainv(i,1) + &
                                               fion(2,isa) * ainv(i,2) + &
                                               fion(3,isa) * ainv(i,3) )
             !
          END DO
          !
       END DO
       !
    ELSE IF ( tnosep ) THEN
       !
       isa = 0
       !
       DO is = 1, nsp
          !
          DO ia = 1, na(is)
             !
             isa = isa + 1
             !
             DO i = 1, 3
                !
                fionm(i,isa) = ainv(i,1) * fion(1,isa) + &
                               ainv(i,2) * fion(2,isa) + &
                               ainv(i,3) * fion(3,isa) - &
                               vnhp(1,atm2nhp(isa)) * vels(i,isa) * pmass(is) - &
                               pmass(is) * ( hgamma(i,1) * vels(1,isa) + &
                                             hgamma(i,2) * vels(2,isa) + &
                                             hgamma(i,3) * vels(3,isa) )
                !
             END DO
             !
             tausp(:,isa) = 2.D0 * taus(:,isa) - tausm(:,isa) + &
                            dt2 * iforce(:,isa) * fionm(:,isa) / pmass(is)
             !
             velsp(:,isa) = velsm(:,isa) + twodel * fionm(:,isa) / pmass(is)
             !
          END DO
          !
       END DO
       !
    ELSE
       !
       isa = 0
       !
       DO is = 1, nsp
          !
          DO ia = 1, na(is)
             !
             isa = isa + 1
             !
             DO i = 1, 3
                !
                tausp(i,isa) = verl1 * taus(i,isa) + verl2 * tausm(i,isa) + &
                               verl3 / pmass(is) * iforce(i,isa) * &
                                               ( ainv(i,1) * fion(1,isa) + &
                                                 ainv(i,2) * fion(2,isa) + &
                                                 ainv(i,3) * fion(3,isa) ) - &
                               verl3 * iforce(i,isa) * &
                                               ( hgamma(i,1) * vels(1,isa) + &
                                                 hgamma(i,2) * vels(2,isa) + &
                                                 hgamma(i,3) * vels(3,isa) )
                !
                velsp(i,isa) = velsm(i,isa) - 4.D0 * fricp * vels(i,isa) + &
                               twodel / pmass(is) * iforce(i,isa) * &
                                               ( ainv(i,1) * fion(1,isa) + &
                                                 ainv(i,2) * fion(2,isa) + &
                                                 ainv(i,3) * fion(3,isa) ) - &
                               twodel * iforce(i,isa) * &
                                               ( hgamma(i,1) * vels(1,isa) + &
                                                 hgamma(i,2) * vels(2,isa) + &
                                                 hgamma(i,3) * vels(3,isa) )
                !
             END DO
             !
          END DO
          !
       END DO
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE ions_move
  !
!------------------------------------------------------------------------------!
END MODULE ions_positions
!------------------------------------------------------------------------------!
