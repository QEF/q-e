!
! Copyright (C) 2002-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
MODULE ions_positions
!------------------------------------------------------------------------------!
  !! Module containing atomic positions arrays used in the CP codes during
  !! the dynamics.
  !
  USE kinds,             ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), TARGET, ALLOCATABLE :: tau0(:,:), taum(:,:),  taup(:,:)
  REAL(DP), TARGET, ALLOCATABLE :: taus(:,:), tausm(:,:), tausp(:,:)
  REAL(DP), TARGET, ALLOCATABLE :: vels(:,:), velsm(:,:), velsp(:,:)
  REAL(DP), TARGET, ALLOCATABLE :: fion(:,:), fionm(:,:), fionp(:,:)
  !INTEGER,  TARGET, ALLOCATABLE :: ityp(:), mobil(:,:)
  !
  CONTAINS 
  !
  !  ... meaning of some variables appearing in the following subs.
  !
  !  nsp   number of atomic species
  !  nax   maximum number of atoms per specie
  !  nat   total number of atoms
  !  na(:) number of atoms per specie
  !  pmass(:)   mass (converted to a.u.) of ions
  !
  !
  SUBROUTINE allocate_ions_positions( nsp, nat )
     !! Allocate ionic dynamics arrays
     !
     INTEGER, INTENT(IN) :: nsp
     !! number of atomic species
     INTEGER, INTENT(IN) :: nat
     !! total number of atoms
     !
     IF( ALLOCATED( tau0  ) ) DEALLOCATE( tau0  )
     IF( ALLOCATED( taum  ) ) DEALLOCATE( taum  ) 
     IF( ALLOCATED( taup  ) ) DEALLOCATE( taup  )
     IF( ALLOCATED( taus  ) ) DEALLOCATE( taus  )
     IF( ALLOCATED( tausm ) ) DEALLOCATE( tausm )
     IF( ALLOCATED( tausp ) ) DEALLOCATE( tausp )
     IF( ALLOCATED( vels  ) ) DEALLOCATE( vels  )
     IF( ALLOCATED( velsm ) ) DEALLOCATE( velsm )
     IF( ALLOCATED( velsp ) ) DEALLOCATE( velsp )
     IF( ALLOCATED( fion  ) ) DEALLOCATE( fion  )
     IF( ALLOCATED( fionm ) ) DEALLOCATE( fionm )
     IF( ALLOCATED( fionp ) ) DEALLOCATE( fionp )
     !IF( ALLOCATED( ityp  ) ) DEALLOCATE( ityp  )
     !IF( ALLOCATED( mobil ) ) DEALLOCATE( mobil )
     !
     ALLOCATE( tau0( 3, nat ) )
     ALLOCATE( taum( 3, nat ) )
     ALLOCATE( taup( 3, nat ) )
     ALLOCATE( taus( 3, nat ) )
     ALLOCATE( tausm( 3, nat ) )
     ALLOCATE( tausp( 3, nat ) )
     ALLOCATE( vels( 3, nat )  )
     ALLOCATE( velsm( 3, nat ) )
     ALLOCATE( velsp( 3, nat ) )
     ALLOCATE( fion( 3, nat )  )
     ALLOCATE( fionm( 3, nat ) )
     ALLOCATE( fionp( 3, nat ) )
     !ALLOCATE( ityp( nat ) )
     !ALLOCATE( mobil( 3, nat ) )
     !
     RETURN
  END SUBROUTINE allocate_ions_positions

  !--------------------------------------------------------------------------

  SUBROUTINE deallocate_ions_positions( )
     !! Deallocate ionic dynamics arrays.
     !
     IF( ALLOCATED( tau0  ) ) DEALLOCATE( tau0  )
     IF( ALLOCATED( taum  ) ) DEALLOCATE( taum  )
     IF( ALLOCATED( taup  ) ) DEALLOCATE( taup  )
     IF( ALLOCATED( taus  ) ) DEALLOCATE( taus  )
     IF( ALLOCATED( tausm ) ) DEALLOCATE( tausm )
     IF( ALLOCATED( tausp ) ) DEALLOCATE( tausp )
     IF( ALLOCATED( vels  ) ) DEALLOCATE( vels  )
     IF( ALLOCATED( velsm ) ) DEALLOCATE( velsm )
     IF( ALLOCATED( velsp ) ) DEALLOCATE( velsp )
     IF( ALLOCATED( fion  ) ) DEALLOCATE( fion  )
     IF( ALLOCATED( fionm ) ) DEALLOCATE( fionm )
     IF( ALLOCATED( fionp ) ) DEALLOCATE( fionp )
     !IF( ALLOCATED( ityp  ) ) DEALLOCATE( ityp  )
     !IF( ALLOCATED( mobil ) ) DEALLOCATE( mobil )
     RETURN
  END SUBROUTINE deallocate_ions_positions

  !--------------------------------------------------------------------------
  SUBROUTINE ions_hmove( taus, tausm, iforce, pmass, fion, ainv, delt, ityp, nat )
    !--------------------------------------------------------------------------
    !
    REAL(DP), INTENT(IN)  :: tausm(:,:), pmass(:), fion(:,:)
    INTEGER,  INTENT(IN)  :: iforce(:,:)
    REAL(DP), INTENT(IN)  :: ainv(3,3), delt
    REAL(DP), INTENT(OUT) :: taus(:,:) 
    INTEGER,  INTENT(IN)  :: ityp(:), nat
    INTEGER               :: ia, i
    REAL(DP)              :: dt2by2, fac, fions(3)
    !
    !
    dt2by2 = 0.5D0 * delt * delt
    !
    DO ia = 1, nat
       !
       fac = dt2by2  / pmass(ityp(ia))
       !
       DO i = 1, 3
          !
          fions(i) = fion(1,ia) * ainv(i,1) + &
                     fion(2,ia) * ainv(i,2) + &
                     fion(3,ia) * ainv(i,3)
          !
       END DO
       !
       taus(:,ia) = tausm(:,ia) + iforce(:,ia) * fac * fions(:)
       !
    END DO
    !
    RETURN
    !
  END SUBROUTINE ions_hmove
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ions_move( tausp, taus, tausm, iforce, pmass, fion, ainv, &
                        delt, ityp, nat, fricp, hgamma, vels, tsdp, tnosep, &
                        fionm, vnhp, velsp, velsm, nhpcl, nhpdim, atm2nhp )
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)    :: taus(:,:), tausm(:,:), pmass(:), fion(:,:)
    INTEGER,  INTENT(IN)    :: iforce(:,:)
    REAL(DP), INTENT(IN)    :: ainv(3,3), delt
    REAL(DP), INTENT(OUT)   :: tausp(:,:)
    INTEGER,  INTENT(IN)    :: ityp(:), nat, nhpcl, nhpdim, atm2nhp(:)
    REAL(DP), INTENT(IN)    :: fricp, hgamma(3,3), vels(:,:)
    LOGICAL,  INTENT(IN)    :: tsdp, tnosep
    REAL(DP), INTENT(INOUT) :: fionm(:,:)
    REAL(DP), INTENT(IN)    :: vnhp(nhpcl,nhpdim)
    REAL(DP), INTENT(OUT)   :: velsp(:,:)
    REAL(DP), INTENT(IN)    :: velsm(:,:)
    INTEGER                 :: is, ia, i
    REAL(DP)                :: dt2by2, fac, dt2, twodel
    REAL(DP)                :: verl1, verl2, verl3
    REAL(DP)                :: ftmp(3)
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
       DO ia = 1, nat
          is = ityp(ia)
          !
          DO i = 1, 3
             !
             tausp(i,ia) = taus(i,ia) - pmass(is) * &
                                          ( hgamma(i,1) * vels(1,ia) + &
                                            hgamma(i,2) * vels(2,ia) + &
                                            hgamma(i,3) * vels(3,ia) ) + &
                            iforce(i,ia) * dt2 / pmass(is) * &
                                          ( fion(1,ia) * ainv(i,1) + &
                                            fion(2,ia) * ainv(i,2) + &
                                            fion(3,ia) * ainv(i,3) )
             !
          END DO
          !
       END DO
       !
    ELSE IF ( tnosep ) THEN
       !
       DO ia = 1, nat
          is = ityp(ia)
          !
          DO i = 1, 3
             !
             fionm(i,ia) = ainv(i,1) * fion(1,ia) + &
                            ainv(i,2) * fion(2,ia) + &
                            ainv(i,3) * fion(3,ia) - &
                            vnhp(1,atm2nhp(ia)) * vels(i,ia) * pmass(is) - &
                            pmass(is) * ( hgamma(i,1) * vels(1,ia) + &
                                          hgamma(i,2) * vels(2,ia) + &
                                          hgamma(i,3) * vels(3,ia) )
             !
          END DO
          !
          tausp(:,ia) = 2.D0 * taus(:,ia) - tausm(:,ia) + &
                         dt2 * iforce(:,ia) * fionm(:,ia) / pmass(is)
          !
          velsp(:,ia) = velsm(:,ia) + twodel * fionm(:,ia) / pmass(is)
          !
       END DO
       !
    ELSE
       !
       DO ia = 1, nat
          is = ityp(ia)
          !
          DO i = 1, 3
             !
             tausp(i,ia) = verl1 * taus(i,ia) + verl2 * tausm(i,ia) + &
                            verl3 / pmass(is) * iforce(i,ia) * &
                                            ( ainv(i,1) * fion(1,ia) + &
                                              ainv(i,2) * fion(2,ia) + &
                                              ainv(i,3) * fion(3,ia) ) - &
                            verl3 * iforce(i,ia) * &
                                            ( hgamma(i,1) * vels(1,ia) + &
                                              hgamma(i,2) * vels(2,ia) + &
                                              hgamma(i,3) * vels(3,ia) )
             !
             velsp(i,ia) = velsm(i,ia) - 4.D0 * fricp * vels(i,ia) + &
                            twodel / pmass(is) * iforce(i,ia) * &
                                            ( ainv(i,1) * fion(1,ia) + &
                                              ainv(i,2) * fion(2,ia) + &
                                              ainv(i,3) * fion(3,ia) ) - &
                            twodel * iforce(i,ia) * &
                                            ( hgamma(i,1) * vels(1,ia) + &
                                              hgamma(i,2) * vels(2,ia) + &
                                              hgamma(i,3) * vels(3,ia) )
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
  !
  SUBROUTINE set_velocities( tausm, taus0, vels, iforce, nat, delt)
     USE kinds, ONLY : DP
     IMPLICIT NONE
     INTEGER,  INTENT(IN) :: nat
     REAL(DP)             :: tausm( 3, nat ), taus0( 3, nat )
     REAL(DP), INTENT(IN) :: delt
     REAL(DP), INTENT(IN) :: vels( 3, nat )
     INTEGER,  INTENT(IN) :: iforce( 3, nat )
     INTEGER :: i, ia
     DO ia = 1, nat
       tausm( :, ia ) = taus0( :, ia )
       DO i = 1, 3
         IF( iforce( i, ia ) > 0 ) THEN
           taus0( i, ia ) = taus0( i, ia ) + vels( i, ia ) * delt
         END IF
       ENDDO
     END DO
     RETURN
  END SUBROUTINE set_velocities
  !
!------------------------------------------------------------------------------!
END MODULE ions_positions
!------------------------------------------------------------------------------!
