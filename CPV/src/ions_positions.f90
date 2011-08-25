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
  USE kinds,             ONLY : DP
  USE atoms_type_module, ONLY : atoms_type, atoms_type_init
  !
  IMPLICIT NONE
  !
  ! ... Atomic positions arrays used in the cp codes during the dynamic
  !
  REAL(DP), TARGET, ALLOCATABLE :: tau0(:,:), taum(:,:),  taup(:,:)
  REAL(DP), TARGET, ALLOCATABLE :: taus(:,:), tausm(:,:), tausp(:,:)
  REAL(DP), TARGET, ALLOCATABLE :: vels(:,:), velsm(:,:), velsp(:,:)
  REAL(DP), TARGET, ALLOCATABLE :: fion(:,:), fionm(:,:), fionp(:,:)
  INTEGER,  TARGET, ALLOCATABLE :: ityp(:), mobil(:,:)
  !
  TYPE (atoms_type) :: atoms0, atomsp, atomsm
  !
  CONTAINS 
  !
  !  ... meaning of some variables appearing in the folloving subs.
  !
  !  nsp   number of atomic species
  !  nax   maximum number of atoms per specie
  !  nat   total number of atoms
  !  na(:) number of atoms per specie
  !  pmass(:)   mass (converted to a.u.) of ions
  !
  !
  !
  SUBROUTINE allocate_ions_positions( nsp, nat )
     INTEGER, INTENT(IN) :: nsp, nat
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
     IF( ALLOCATED( ityp  ) ) DEALLOCATE( ityp  )
     IF( ALLOCATED( mobil ) ) DEALLOCATE( mobil )
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
     ALLOCATE( ityp( nat ) )
     ALLOCATE( mobil( 3, nat ) )
     !
     NULLIFY( atoms0 % taur   )
     NULLIFY( atoms0 % taus   )
     NULLIFY( atoms0 % vels   )
     NULLIFY( atoms0 % for    )
     NULLIFY( atoms0 % mobile )
     NULLIFY( atoms0 % ityp   )
     NULLIFY( atomsm % taur   )
     NULLIFY( atomsm % taus   )
     NULLIFY( atomsm % vels   )
     NULLIFY( atomsm % for    )
     NULLIFY( atomsm % mobile )
     NULLIFY( atomsm % ityp   )
     NULLIFY( atomsp % taur   )
     NULLIFY( atomsp % taus   )
     NULLIFY( atomsp % vels   )
     NULLIFY( atomsp % for    )
     NULLIFY( atomsp % mobile )
     NULLIFY( atomsp % ityp   )
     !
     atoms0 % taur   => tau0
     atoms0 % taus   => taus
     atoms0 % vels   => vels
     atoms0 % for    => fion
     atoms0 % mobile => mobil
     atoms0 % ityp   => ityp
     atomsm % taur   => taum
     atomsm % taus   => tausm
     atomsm % vels   => velsm
     atomsm % for    => fionm
     atomsm % mobile => mobil
     atomsm % ityp   => ityp
     atomsp % taur   => taup
     atomsp % taus   => tausp
     atomsp % vels   => velsp
     atomsp % for    => fionp
     atomsp % mobile => mobil
     atomsp % ityp   => ityp
     !
     RETURN
  END SUBROUTINE allocate_ions_positions

  !--------------------------------------------------------------------------

  SUBROUTINE deallocate_ions_positions( )
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
     IF( ALLOCATED( ityp  ) ) DEALLOCATE( ityp  )
     IF( ALLOCATED( mobil ) ) DEALLOCATE( mobil )
     NULLIFY( atoms0 % taur   )
     NULLIFY( atoms0 % taus   )
     NULLIFY( atoms0 % vels   )
     NULLIFY( atoms0 % for    )
     NULLIFY( atoms0 % mobile )
     NULLIFY( atoms0 % ityp   )
     NULLIFY( atomsm % taur   )
     NULLIFY( atomsm % taus   )
     NULLIFY( atomsm % vels   )
     NULLIFY( atomsm % for    )
     NULLIFY( atomsm % mobile )
     NULLIFY( atomsm % ityp   )
     NULLIFY( atomsp % taur   )
     NULLIFY( atomsp % taus   )
     NULLIFY( atomsp % vels   )
     NULLIFY( atomsp % for    )
     NULLIFY( atomsp % mobile )
     NULLIFY( atomsp % ityp   )
     RETURN
  END SUBROUTINE deallocate_ions_positions



  !--------------------------------------------------------------------------
  SUBROUTINE ions_hmove( taus, tausm, iforce, pmass, fion, ainv, delt, na, nsp )
    !--------------------------------------------------------------------------
    !
    REAL(DP), INTENT(IN)  :: tausm(:,:), pmass(:), fion(:,:)
    INTEGER,  INTENT(IN)  :: iforce(:,:)
    REAL(DP), INTENT(IN)  :: ainv(3,3), delt
    REAL(DP), INTENT(OUT) :: taus(:,:) 
    INTEGER,  INTENT(IN)  :: na(:), nsp
    INTEGER               :: is, ia, i, isa
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
    INTEGER,  INTENT(IN)    :: iforce(:,:)
    REAL(DP), INTENT(IN)    :: ainv(3,3), delt
    REAL(DP), INTENT(OUT)   :: tausp(:,:)
    INTEGER,  INTENT(IN)    :: na(:), nsp, nhpcl, nhpdim, atm2nhp(:)
    REAL(DP), INTENT(IN)    :: fricp, hgamma(3,3), vels(:,:)
    LOGICAL,  INTENT(IN)    :: tsdp, tnosep
    REAL(DP), INTENT(INOUT) :: fionm(:,:)
    REAL(DP), INTENT(IN)    :: vnhp(nhpcl,nhpdim)
    REAL(DP), INTENT(OUT)   :: velsp(:,:)
    REAL(DP), INTENT(IN)    :: velsm(:,:)
    INTEGER                 :: is, ia, i, isa
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
                tausp(i,isa) = taus(i,isa) - pmass(is) * &
                                             ( hgamma(i,1) * vels(1,isa) + &
                                               hgamma(i,2) * vels(2,isa) + &
                                               hgamma(i,3) * vels(3,isa) ) + &
                               iforce(i,isa) * dt2 / pmass(is) * &
                                             ( fion(1,isa) * ainv(i,1) + &
                                               fion(2,isa) * ainv(i,2) + &
                                               fion(3,isa) * ainv(i,3) )
                !
             END DO
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
  !
  !
  !

       SUBROUTINE atoms_init(atoms_m, atoms_0, atoms_p, stau, ind_srt, if_pos, atml, h, nat, nsp, na, pmass )

         !   Allocate and fill the three atoms structure using scaled position an cell

         USE printout_base, ONLY : printout_pos
         USE io_global,     ONLY : ionode, stdout

         IMPLICIT NONE

         TYPE (atoms_type)    :: atoms_0, atoms_p, atoms_m
         REAL(DP), INTENT(IN) :: h( 3, 3 )
         REAL(DP), INTENT(IN) :: stau(:,:)
         CHARACTER(LEN=3), INTENT(IN) :: atml(:)
         INTEGER, INTENT(IN) :: ind_srt( : )
         INTEGER, INTENT(IN) :: if_pos( :, : )
         INTEGER, INTENT(IN) :: nat, nsp
         INTEGER, INTENT(IN) :: na( : )
         REAL(DP), INTENT(IN) :: pmass( : )
          
         CHARACTER(LEN=3), ALLOCATABLE :: labels( : )
         LOGICAL, ALLOCATABLE :: ismb(:,:)
         INTEGER              :: ia, is, isa
         LOGICAL              :: nofx

         ALLOCATE( ismb( 3, nat ) )

         ismb = .TRUE.
         nofx = .TRUE.
         DO isa = 1, nat
           ia = ind_srt( isa )
           ismb( 1, isa ) = ( if_pos( 1, ia ) /= 0 )
           ismb( 2, isa ) = ( if_pos( 2, ia ) /= 0 )
           ismb( 3, isa ) = ( if_pos( 3, ia ) /= 0 )
           nofx = nofx .AND. ismb( 1, isa )
           nofx = nofx .AND. ismb( 2, isa )
           nofx = nofx .AND. ismb( 3, isa )
         END DO

         CALL atoms_type_init(atoms_m, stau, ismb, atml, pmass, na, nsp, h)
         CALL atoms_type_init(atoms_0, stau, ismb, atml, pmass, na, nsp, h)
         CALL atoms_type_init(atoms_p, stau, ismb, atml, pmass, na, nsp, h)

         IF( ionode ) THEN
            !
            ALLOCATE( labels(  nat ) )
            !
            isa = 0
            DO is = 1, nsp
               DO ia = 1, na( is )
                  isa = isa + 1
                  labels( isa ) = atml( is )
               END DO
            END DO

            WRITE( stdout, * )

            CALL printout_pos( stdout, stau, nat, label = labels, &
                 head = 'Scaled positions from standard input' )

            IF( .NOT. nofx ) THEN
               WRITE( stdout, 10 )
 10            FORMAT( /, &
                      3X, 'Position components with 0 are kept fixed', /, &
                      3X, '  ia  x  y  z ' )
               DO isa = 1, nat
                  ia = ind_srt( isa )
                  WRITE( stdout, 20 ) isa, if_pos( 1, ia ), if_pos( 2, ia ), if_pos( 3, ia )
               END DO
 20            FORMAT( 3X, I4, I3, I3, I3 )
            END IF

            DEALLOCATE( labels )

         END IF

         DEALLOCATE( ismb )

         RETURN
       END SUBROUTINE atoms_init

!  --------------------------------------------------------------------------   

      SUBROUTINE ions_shiftval(atoms_m, atoms_0, atoms_p)

        !  Update ionic positions and velocities in atoms structures

        IMPLICIT NONE
        TYPE(atoms_type) :: atoms_m, atoms_0, atoms_p
        INTEGER :: is, ia, ub

          ub = atoms_m%nat
          atoms_m%taus(1:3,1:ub) = atoms_0%taus(1:3,1:ub)
          atoms_m%vels(1:3,1:ub) = atoms_0%vels(1:3,1:ub)
          atoms_m%for(1:3,1:ub) = atoms_0%for(1:3,1:ub)
          atoms_0%taus(1:3,1:ub) = atoms_p%taus(1:3,1:ub)
          atoms_0%vels(1:3,1:ub) = atoms_p%vels(1:3,1:ub)
          atoms_0%for(1:3,1:ub) = atoms_p%for(1:3,1:ub)

        RETURN
      END SUBROUTINE ions_shiftval




      REAL(DP) FUNCTION max_ion_forces( atoms )

        IMPLICIT NONE 
        TYPE (atoms_type) :: atoms
        INTEGER :: ia
        REAL(DP) :: fmax
        fmax = 0.0d0
        DO ia = 1, atoms%nat
          IF( atoms%mobile(1, ia) > 0 ) fmax = MAX( fmax, ABS( atoms%for(1, ia) ) )
          IF( atoms%mobile(2, ia) > 0 ) fmax = MAX( fmax, ABS( atoms%for(2, ia) ) )
          IF( atoms%mobile(3, ia) > 0 ) fmax = MAX( fmax, ABS( atoms%for(3, ia) ) )
        END DO
        max_ion_forces = fmax
        RETURN
      END FUNCTION max_ion_forces

!
!

      SUBROUTINE resort_position( pos, fion, atoms, isrt, ht )

         
        ! This subroutine copys positions and forces into
        ! array "pos" and "for" using the same atoms sequence
        ! as in the input file

        USE cell_base, ONLY: s_to_r
        USE cell_base, ONLY: boxdimensions

        IMPLICIT NONE

        REAL(DP), INTENT(OUT) :: pos(:,:), fion(:,:)
        TYPE (atoms_type), INTENT(IN) :: atoms
        TYPE (boxdimensions), INTENT(IN) :: ht
        INTEGER, INTENT(IN) :: isrt( : )
        INTEGER :: ia, is, isa, ipos

        isa = 0
        DO is = 1, atoms%nsp
          DO ia = 1, atoms%na(is)
            isa  = isa + 1
            ipos = isrt( isa )
            CALL s_to_r( atoms%taus( : , isa ), pos( :, ipos ), ht )
            fion( :, ipos ) = atoms%for( : , isa )
          END DO
        END DO

        RETURN
      END SUBROUTINE resort_position
     


  !
!------------------------------------------------------------------------------!
END MODULE ions_positions
!------------------------------------------------------------------------------!
