!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
  MODULE ions_base
!------------------------------------------------------------------------------!

      USE kinds, ONLY: dbl
      USE parameters, ONLY: nsx, natx, ntypx
!
      IMPLICIT NONE
      SAVE

      !     nsp       = number of species
      !     na(is)    = number of atoms of species is
      !     nax       = max number of atoms of a given species
      !     nat       = total number of atoms of all species

      INTEGER :: nsp     = 0
      INTEGER :: na(nsx) = 0    
      INTEGER :: nax     = 0
      INTEGER :: nat     = 0

      !     zv(is)    = (pseudo-)atomic charge
      !     pmass(is) = mass (converted to a.u.) of ions
      !     rcmax(is) = Ewald radius (for ion-ion interactions)

      REAL(dbl) :: zv(nsx)    = 0.0d0
      REAL(dbl) :: pmass(nsx) = 0.0d0
      REAL(dbl) :: amass(nsx) = 0.0d0
      REAL(dbl) :: rcmax(nsx) = 0.0d0

      !     ityp( i ) = the type of i-th atom in stdin
      !     atm( j )  = name of the type of the j-th atomic specie
      !     tau( 1:3, i ) = position of the i-th atom

      INTEGER,   ALLOCATABLE :: ityp(:)
      REAL(dbl), ALLOCATABLE :: tau(:,:)      !  initial position
      REAL(dbl), ALLOCATABLE :: vel(:,:)      !  initial velocities
      REAL(dbl), ALLOCATABLE :: tau_srt(:,:)  !  tau sorted by specie
      REAL(dbl), ALLOCATABLE :: vel_srt(:,:)  !  vel sorted by specie
      INTEGER,   ALLOCATABLE :: ind_srt( : )  !  index of tau sorted by specie
      CHARACTER(LEN=3  ) :: atm( ntypx ) 
      CHARACTER(LEN=80 ) :: tau_units

      ! if if_pos( x, i ) = 0 then 
      !    x coordinate of the i-th atom will be kept fixed

      INTEGER, ALLOCATABLE :: if_pos(:,:)  
      INTEGER :: fixatom  !!! to be removed

      INTEGER :: ind_localisation(natx) = 0   ! true if we want to know the localization arount the atom
      INTEGER :: nat_localisation = 0 
      LOGICAL :: print_localisation = .FALSE. ! Calculates hartree energy around specified atoms
      INTEGER :: self_interaction = 0 
      REAL(dbl) :: si_epsilon = 0.0d0
      REAL(dbl) :: rad_localisation = 0.0d0
      REAL(dbl), ALLOCATABLE :: pos_localisation(:,:)

      LOGICAL :: tions_base_init = .FALSE.
      LOGICAL, PRIVATE :: tdebug = .FALSE.
      
      INTERFACE ions_vel
         MODULE PROCEDURE ions_vel3, ions_vel2
      END INTERFACE
!

!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!

    SUBROUTINE packtau( taup, tau, na, nsp )
      IMPLICIT NONE
      REAL(dbl), INTENT(OUT) :: taup( :, : )
      REAL(dbl), INTENT(IN) :: tau( :, :, : )
      INTEGER, INTENT(IN) :: na( : ), nsp
      INTEGER :: is, ia, isa
      isa = 0
      DO is = 1, nsp
        DO ia = 1, na( is )
          isa = isa + 1
          taup( :, isa ) = tau( :, ia, is )
        END DO
      END DO
      RETURN
    END SUBROUTINE

    SUBROUTINE unpacktau( tau, taup, na, nsp )
      IMPLICIT NONE
      REAL(dbl), INTENT(IN) :: taup( :, : )
      REAL(dbl), INTENT(OUT) :: tau( :, :, : )
      INTEGER, INTENT(IN) :: na( : ), nsp
      INTEGER :: is, ia, isa
      isa = 0
      DO is = 1, nsp
        DO ia = 1, na( is )
          isa = isa + 1
          tau( :, ia, is ) = taup( :, isa )
        END DO
      END DO
      RETURN
    END SUBROUTINE


    SUBROUTINE sort_tau( tausrt, isrt, tau, isp, nat, nsp )
      IMPLICIT NONE
      REAL(dbl), INTENT(OUT) :: tausrt( :, : )
      INTEGER, INTENT(OUT) :: isrt( : )
      REAL(dbl), INTENT(IN) :: tau( :, : )
      INTEGER, INTENT(IN) :: nat, nsp, isp( : )
      INTEGER :: ina( nsp ), na( nsp )
      INTEGER :: is, ia

      ! ... count the atoms for each specie
      na  = 0
      DO ia = 1, nat
        is  =  isp( ia )
        IF( is < 1 .OR. is > nsp ) &
          CALL errore(' sorttau ', ' wrong species index for positions ', ia )
        na( is ) = na( is ) + 1
      END DO

      ! ... compute the index of the first atom in each specie
      ina( 1 ) = 0
      DO is = 2, nsp
        ina( is ) = ina( is - 1 ) + na( is - 1 )
      END DO

      ! ... sort the position according to atomic specie
      na  = 0
      DO ia = 1, nat
        is  =  isp( ia )
        na( is ) = na( is ) + 1
        tausrt( :, na(is) + ina(is) ) = tau(:, ia )
        isrt  (    na(is) + ina(is) ) = ia
      END DO
      RETURN
    END SUBROUTINE


    SUBROUTINE unsort_tau( tau, tausrt, isrt, nat )
      IMPLICIT NONE
      REAL(dbl), INTENT(IN) :: tausrt( :, : )
      INTEGER, INTENT(IN) :: isrt( : )
      REAL(dbl), INTENT(OUT) :: tau( :, : )
      INTEGER, INTENT(IN) :: nat
      INTEGER :: isa, ia
      DO isa = 1, nat
        ia  =  isrt( isa )
        tau( :, ia ) = tausrt( :, isa )
      END DO
      RETURN
    END SUBROUTINE





    SUBROUTINE ions_base_init( nsp_ , nat_ , na_ , ityp_ , tau_ , vel_, amass_ , &
        atm_ , if_pos_ , tau_units_ , id_loc_ , sic_ , sic_epsilon_, sic_rloc_ )

      USE constants, ONLY: scmass
      USE io_base, ONLY: stdout

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nsp_ , nat_ , na_ (:) , ityp_ (:)
      REAL(dbl), INTENT(IN) :: tau_(:,:)
      REAL(dbl), INTENT(IN) :: vel_(:,:)
      REAL(dbl), INTENT(IN) :: amass_(:)
      CHARACTER(LEN=*), INTENT(IN) :: atm_ (:)
      CHARACTER(LEN=*), INTENT(IN) :: tau_units_
      INTEGER, INTENT(IN) :: if_pos_ (:,:)
      INTEGER, OPTIONAL, INTENT(IN) :: id_loc_ (:)
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: sic_
      REAL(dbl), OPTIONAL, INTENT(IN) :: sic_epsilon_
      REAL(dbl), OPTIONAL, INTENT(IN) :: sic_rloc_
      INTEGER :: i, ia

      nsp = nsp_
      nat = nat_

      if(nat < 1) &
        call errore(' ions_base_init ', ' NAX OUT OF RANGE ',1)
      if(nsp < 1) &
        call errore(' ions_base_init ',' NSP OUT OF RANGE ',1)
      if(nsp > SIZE( na ) ) &
        call errore(' ions_base_init ',' NSP too large, increase NSX parameter ',1)

      na( 1:nsp ) = na_ ( 1:nsp )
      nax = MAXVAL( na( 1:nsp ) )

      atm( 1:nsp ) = atm_ ( 1:nsp )
      tau_units    = TRIM( tau_units_ )

      if ( nat /= SUM( na( 1:nsp ) ) ) &
        call errore(' ions_base_init ',' inconsistent NAT and NA ',1)

      ALLOCATE( ityp( nat ) )
      ALLOCATE( tau( 3, nat ) )
      ALLOCATE( vel( 3, nat ) )
      ALLOCATE( tau_srt( 3, nat ) )
      ALLOCATE( vel_srt( 3, nat ) )
      ALLOCATE( ind_srt( nat ) )
      ALLOCATE( if_pos( 3, nat ) )

      ityp( 1:nat )      = ityp_ ( 1:nat )
      tau( : , 1:nat )   = tau_ ( : , 1:nat )
      vel( : , 1:nat )   = vel_ ( : , 1:nat )
      if_pos( :, 1:nat ) = if_pos_ ( : , 1:nat )

! ...     tau_srt : atomic species are ordered according to
! ...     the ATOMIC_SPECIES input card. Within each specie atoms are ordered
! ...     according to the ATOMIC_POSITIONS input card.
! ...     ind_srt : can be used to restore the origina position

      CALL sort_tau( tau_srt, ind_srt, tau, ityp, nat, nsp )

      DO ia = 1, nat
        vel_srt( :, ia ) = vel( :, ind_srt( ia ) )
      END DO

      IF( tdebug ) THEN
        WRITE( stdout, * ) 'ions_base_init: unsorted position and velocities'
        DO ia = 1, nat
          WRITE( stdout, fmt="(A3,3D12.4,3X,3D12.4)") &
            atm( ityp( ia ) ), tau(1:3, ia), vel(1:3,ia)
        END DO
        WRITE( stdout, * ) 'ions_base_init: sorted position and velocities'
        DO ia = 1, nat
          WRITE( stdout, fmt="(A3,3D12.4,3X,3D12.4)") &
            atm( ityp( ind_srt( ia ) ) ), tau_srt(1:3, ia), vel_srt(1:3,ia)
        END DO
      END IF

      !
      ! ... The constrain on fixed coordinates is implemented using the array
      ! ... if_pos whose value is 0 when the coordinate is to be kept fixed, 1
      ! ... otherwise. fixatom is maintained for compatibility. ( C.S. 15/10/2003 )
      !
      if_pos = 1
      if_pos(:,:) = if_pos_ (:,1:nat)

      IF( PRESENT( sic_ ) ) THEN
        select case ( TRIM( sic_ ) )
        case ( 'sic_pz' ) 
          self_interaction = 1
        case ( 'sic_mac' )
          self_interaction = 2
        case ( 'only_sich' )
          self_interaction = 3
        case ( 'only_sicxc_pz' )
          self_interaction = -1
        case ( 'only_sicxc_mac' )
          self_interaction = -2
        case default
          self_interaction = 0
        end select
      END IF
      IF( PRESENT( sic_epsilon_ ) ) THEN
        si_epsilon       = sic_epsilon_
      END IF
      IF( PRESENT( sic_rloc_ ) ) THEN
        rad_localisation = sic_rloc_
      END IF
      IF( PRESENT( id_loc_ ) ) THEN
        ind_localisation(1:nat) = id_loc_ ( 1:nat )
        nat_localisation = COUNT( ind_localisation > 0 ) 
        ALLOCATE( pos_localisation( 4, nat_localisation ) )
      !counting the atoms around which i want to calculate the charge localization
      ELSE
        ind_localisation(1:nat) = 0
        nat_localisation = 0
      END IF
      !
      IF( nat_localisation > 0 ) print_localisation = .TRUE.
      !
      ! ... TEMP: calculate fixatom (to be removed)
      !
      fixatom = 0
      fix1: DO ia = nat, 1, -1
        IF ( if_pos(1,ia) /= 0 .OR. &
             if_pos(2,ia) /= 0 .OR. &
             if_pos(3,ia) /= 0 ) EXIT fix1
        fixatom = fixatom + 1
      END DO fix1

      amass( 1:nsp ) = amass_ ( 1:nsp )
      IF( ANY( amass( 1:nsp ) <= 0.0d0 ) ) &
        CALL errore( ' ions_base_init ', ' invalid  mass ', 1 ) 
      pmass( 1:nsp ) = amass_ ( 1:nsp ) * scmass

      tions_base_init = .TRUE.

      RETURN
    END SUBROUTINE
    

    SUBROUTINE deallocate_ions_base()
      IMPLICIT NONE
      IF( ALLOCATED( ityp ) ) DEALLOCATE( ityp )
      IF( ALLOCATED( tau ) ) DEALLOCATE( tau )
      IF( ALLOCATED( vel ) ) DEALLOCATE( vel )
      IF( ALLOCATED( tau_srt ) ) DEALLOCATE( tau_srt )
      IF( ALLOCATED( vel_srt ) ) DEALLOCATE( vel_srt )
      IF( ALLOCATED( ind_srt ) ) DEALLOCATE( ind_srt )
      IF( ALLOCATED( if_pos ) ) DEALLOCATE( if_pos )
      IF( ALLOCATED( pos_localisation ) ) DEALLOCATE( pos_localisation )
      tions_base_init = .FALSE.
      RETURN
    END SUBROUTINE


    SUBROUTINE ions_vel3( vel, taup, taum, na, nsp, dt )
      IMPLICIT NONE
      REAL(dbl) :: vel(:,:,:), taup(:,:,:), taum(:,:,:)
      INTEGER :: na(:), nsp
      REAL(dbl) :: dt
      INTEGER :: ia, is, i
      REAL(dbl) :: fac
      fac  = 1.0d0 / ( dt * 2.0d0 )
      DO is = 1, nsp
        DO ia = 1, na(is)
          DO i = 1, 3
            vel(i,ia,is) = ( taup(i,ia,is) - taum(i,ia,is) ) * fac
          END DO
        END DO
      END DO
      RETURN
    END SUBROUTINE


    SUBROUTINE ions_vel2( vel, taup, taum, nat, dt )
      IMPLICIT NONE
      REAL(dbl) :: vel(:,:), taup(:,:), taum(:,:)
      INTEGER :: nat
      REAL(dbl) :: dt
      INTEGER :: ia, i
      REAL(dbl) :: fac
      fac  = 1.0d0 / ( dt * 2.0d0 )
      DO ia = 1, nat
        DO i = 1, 3
          vel(i,ia) = ( taup(i,ia) - taum(i,ia) ) * fac
        END DO
      END DO
      RETURN
    END SUBROUTINE
 

!------------------------------------------------------------------------------!
  END MODULE ions_base
!------------------------------------------------------------------------------!
