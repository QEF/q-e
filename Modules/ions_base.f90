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
      INTEGER :: isort_pos(natx,nsx) = 0

      !     zv(is)    = (pseudo-)atomic charge
      !     pmass(is) = mass (converted to a.u.) of ions
      !     rcmax(is) = Ewald radius (for ion-ion interactions)

      REAL(dbl) :: zv(nsx)    = 0.0d0
      REAL(dbl) :: pmass(nsx) = 0.0d0
      REAL(dbl) :: amass(nsx) = 0.0d0
      REAL(dbl) :: rcmax(nsx) = 0.0d0

      INTEGER :: ipp(nsx) = 0

      !     ityp( i ) = the type of i-th atom in stdin
      !     atm( j )  = name of the type of the j-th atomic specie
      !     tau( 1:3, i ) = position of the i-th atom

      INTEGER, ALLOCATABLE :: ityp(:)
      REAL(dbl), ALLOCATABLE :: tau(:,:)
      CHARACTER(LEN=3 ) :: atm(ntypx) 

      ! if if_pos( x, i ) = 0 then 
      !    x coordinate of the i-th atom will be kept fixed

      INTEGER, ALLOCATABLE :: if_pos(:,:)  
      INTEGER :: fixatom  !!! to be removed


      LOGICAL :: tions_base_init = .FALSE.
      
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



    SUBROUTINE ions_base_init( nsp_ , nat_ , na_ , ityp_ , tau_ , amass_ , &
        atm_ , if_pos_ )
      USE constants, ONLY: scmass
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nsp_ , nat_ , na_ (:) , ityp_ (:)
      REAL(dbl), INTENT(IN) :: tau_(:,:)
      REAL(dbl), INTENT(IN) :: amass_(:)
      CHARACTER(LEN=*), INTENT(IN) :: atm_ (:)
      INTEGER, INTENT(IN) :: if_pos_ (:,:)
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

      if ( nat /= SUM( na( 1:nsp ) ) ) &
        call errore(' ions_base_init ',' inconsistent NAT and NA ',1)

      ALLOCATE( ityp( nat ) )
      ALLOCATE( tau( 3, nat ) )
      ALLOCATE( if_pos( 3, nat ) )
      ityp( 1:nat ) = ityp_ ( 1:nat )
      tau( : , 1:nat ) = tau_ ( : , 1:nat )
      if_pos( :, 1:nat ) = if_pos_ ( : , 1:nat )

      !
      ! ... The constrain on fixed coordinates is implemented using the array
      ! ... if_pos whose value is 0 when the coordinate is to be kept fixed, 1
      ! ... otherwise. fixatom is maintained for compatibility. ( C.S. 15/10/2003 )
      !
      if_pos = 1
      if_pos(:,:) = if_pos_ (:,1:nat)

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
      IF( ALLOCATED( if_pos ) ) DEALLOCATE( if_pos )
      RETURN
    END SUBROUTINE

!------------------------------------------------------------------------------!
  END MODULE ions_base
!------------------------------------------------------------------------------!
