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
      REAL(dbl), ALLOCATABLE :: tau(:,:)      !  initial positions read from stdin (in bohr)
      REAL(dbl), ALLOCATABLE :: vel(:,:)      !  initial velocities read from stdin (in bohr)
      REAL(dbl), ALLOCATABLE :: tau_srt(:,:)  !  tau sorted by specie in bohr
      REAL(dbl), ALLOCATABLE :: vel_srt(:,:)  !  vel sorted by specie in bohr
      INTEGER,   ALLOCATABLE :: ind_srt( : )  !  index of tau sorted by specie
      CHARACTER(LEN=3  ) :: atm( ntypx ) 
      CHARACTER(LEN=80 ) :: tau_units


      INTEGER, ALLOCATABLE :: if_pos(:,:)     ! if if_pos( x, i ) = 0 then  x coordinate of 
                                              ! the i-th atom will be kept fixed
      INTEGER, ALLOCATABLE :: iforce(:,:)     ! if_pos sorted by specie 
      INTEGER :: fixatom   = -1               !!! to be removed
      INTEGER :: ndofp     = -1               ! ionic degree of freedom

      REAL(dbl) :: fricp   ! friction parameter for damped dynamics
      REAL(dbl) :: greasp  ! friction parameter for damped dynamics

      LOGICAL :: tions_base_init = .FALSE.
      LOGICAL, PRIVATE :: tdebug = .FALSE.

      
      INTERFACE ions_vel
         MODULE PROCEDURE ions_vel3, ions_vel2
      END INTERFACE

      INTERFACE ions_cofmass
         MODULE PROCEDURE cofmass1, cofmass2
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

!------------------------------------------------------------------------------!

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

!------------------------------------------------------------------------------!

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

!------------------------------------------------------------------------------!

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


!------------------------------------------------------------------------------!



    SUBROUTINE ions_base_init( nsp_ , nat_ , na_ , ityp_ , tau_ , vel_, amass_ , &
        atm_ , if_pos_ , tau_units_ , alat_ , a1_ , a2_ , a3_ , rcmax_ )

      USE constants, ONLY: scmass
      USE io_global, ONLY: stdout

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nsp_ , nat_ , na_ (:) , ityp_ (:)
      REAL(dbl), INTENT(IN) :: tau_(:,:)
      REAL(dbl), INTENT(IN) :: vel_(:,:)
      REAL(dbl), INTENT(IN) :: amass_(:)
      CHARACTER(LEN=*), INTENT(IN) :: atm_ (:)
      CHARACTER(LEN=*), INTENT(IN) :: tau_units_
      INTEGER, INTENT(IN) :: if_pos_ (:,:)
      REAL(dbl), INTENT(IN) :: alat_ , a1_(3) , a2_(3) , a3_(3)
      REAL(dbl), INTENT(IN) :: rcmax_ (:)
      INTEGER :: i, ia, is

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
      ALLOCATE( iforce( 3, nat ) )

      ityp( 1:nat )      = ityp_ ( 1:nat )
      vel( : , 1:nat )   = vel_ ( : , 1:nat )
      if_pos( :, 1:nat ) = if_pos_ ( : , 1:nat )

      ! ...  radii, masses 

      DO is = 1, nsp_
         rcmax ( is ) = rcmax_ ( is )
         IF( rcmax( is ) <= 0.d0 ) THEN
            CALL errore(' ions_base_init ',' invalid rcmax ', is)
         END IF
      END DO

      SELECT CASE ( tau_units )
         !
         !  convert input atomic positions to internally used format:
         !  tau in atomic units
         !
         CASE ('alat')
            !
            !  input atomic positions are divided by a0
            !
            tau( :, 1:nat ) = tau_ ( :, 1:nat ) * alat_
            vel( :, 1:nat ) = vel_ ( :, 1:nat ) * alat_
            !
         CASE ('bohr')
            !
            !  input atomic positions are in a.u.: do nothing
            !
            tau( : , 1:nat )   = tau_ ( : , 1:nat )
            vel( : , 1:nat )   = vel_ ( : , 1:nat )
            !
         CASE ('crystal')
            !
            !  input atomic positions are in crystal axis ("scaled"):
            !
            do ia = 1, nat
              do i = 1, 3
                tau ( i, ia ) = a1_ (i) * tau_( 1, ia ) &
                              + a2_ (i) * tau_( 2, ia ) &
                              + a3_ (i) * tau_( 3, ia )
               end do
            end do
            !
            do ia = 1, nat
              do i = 1, 3
                vel ( i, ia ) = a1_ (i) * vel_( 1, ia ) &
                              + a2_ (i) * vel_( 2, ia ) &
                              + a3_ (i) * vel_( 3, ia )
               end do
            end do
            !
         CASE ('angstrom')
            !
            !  atomic positions in A
            !
            tau( :, 1:nat ) = tau_ ( :, 1:nat ) / 0.529177
            vel( :, 1:nat ) = vel_ ( :, 1:nat ) / 0.529177
            !
         CASE DEFAULT
            !
            CALL errore(' ions_base_init ',' tau_units='//TRIM(tau_units)// &
              ' not implemented ', 1 )
            !
      END SELECT


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

      iforce = 0
      !
      DO ia = 1, nat
        iforce ( :, ia ) = if_pos ( :, ind_srt( ia ) )
      END DO
      !
      ndofp = COUNT( iforce /= 0 )

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
    
!------------------------------------------------------------------------------!

    SUBROUTINE deallocate_ions_base()
      IMPLICIT NONE
      IF( ALLOCATED( ityp ) ) DEALLOCATE( ityp )
      IF( ALLOCATED( tau ) ) DEALLOCATE( tau )
      IF( ALLOCATED( vel ) ) DEALLOCATE( vel )
      IF( ALLOCATED( tau_srt ) ) DEALLOCATE( tau_srt )
      IF( ALLOCATED( vel_srt ) ) DEALLOCATE( vel_srt )
      IF( ALLOCATED( ind_srt ) ) DEALLOCATE( ind_srt )
      IF( ALLOCATED( if_pos ) ) DEALLOCATE( if_pos )
      IF( ALLOCATED( iforce ) ) DEALLOCATE( iforce )
      tions_base_init = .FALSE.
      RETURN
    END SUBROUTINE

!------------------------------------------------------------------------------!

    SUBROUTINE ions_vel3( vel, taup, taum, na, nsp, dt )
      IMPLICIT NONE
      REAL(dbl) :: vel(:,:), taup(:,:), taum(:,:)
      INTEGER :: na(:), nsp
      REAL(dbl) :: dt
      INTEGER :: ia, is, i, isa
      REAL(dbl) :: fac
      fac  = 1.0d0 / ( dt * 2.0d0 )
      isa = 0
      DO is = 1, nsp
        DO ia = 1, na(is)
          isa = isa + 1
          DO i = 1, 3
            vel(i,isa) = ( taup(i,isa) - taum(i,isa) ) * fac
          END DO
        END DO
      END DO
      RETURN
    END SUBROUTINE

!------------------------------------------------------------------------------!

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

    SUBROUTINE cofmass1( tau, pmass, na, nsp, cdm )
      IMPLICIT NONE
      REAL(dbl), INTENT(IN) :: tau(:,:,:), pmass(:)
      REAL(dbl), INTENT(OUT) :: cdm(3)
      INTEGER, INTENT(IN) :: na(:), nsp

      REAL(dbl) :: tmas
      INTEGER :: is, i, ia
!
      tmas=0.0
      do is=1,nsp
         tmas=tmas+na(is)*pmass(is)
      end do
!
      do i=1,3
         cdm(i)=0.0
         do is=1,nsp
            do ia=1,na(is)
               cdm(i)=cdm(i)+tau(i,ia,is)*pmass(is)
            end do
         end do
         cdm(i)=cdm(i)/tmas
      end do
!
      RETURN
    END SUBROUTINE

!------------------------------------------------------------------------------!

    SUBROUTINE cofmass2( tau, pmass, na, nsp, cdm )
      IMPLICIT NONE
      REAL(dbl), INTENT(IN) :: tau(:,:), pmass(:)
      REAL(dbl), INTENT(OUT) :: cdm(3)
      INTEGER, INTENT(IN) :: na(:), nsp

      REAL(dbl) :: tmas
      INTEGER :: is, i, ia, isa
!
      tmas=0.0
      do is=1,nsp
         tmas=tmas+na(is)*pmass(is)
      end do
!
      do i=1,3
         cdm(i)=0.0
         isa = 0
         do is=1,nsp
            do ia=1,na(is)
               isa = isa + 1
               cdm(i)=cdm(i)+tau(i,isa)*pmass(is)
            end do
         end do
         cdm(i)=cdm(i)/tmas
      end do
!
      RETURN
    END SUBROUTINE

!------------------------------------------------------------------------------!


!  BEGIN manual -------------------------------------------------------------

      SUBROUTINE randpos(tau, na, nsp, tranp, amprp, hinv, ifor )

!  Randomize ionic position
!  --------------------------------------------------------------------------
!  END manual ---------------------------------------------------------------

         USE cell_base, ONLY: r_to_s
         USE io_global, ONLY: stdout

         IMPLICIT NONE
         REAL(dbl) :: hinv(3,3)
         REAL(dbl) :: tau(:,:)
         INTEGER, INTENT(IN) :: ifor(:,:), na(:), nsp
         LOGICAL, INTENT(IN) :: tranp(:)
         REAL(dbl), INTENT(IN) :: amprp(:)
         REAL(dbl) :: oldp(3), rand_disp(3), rdisp(3)
         INTEGER :: k, is, isa, isa_s, isa_e, isat

         WRITE( stdout, 600 )

         isat = 0
         DO is = 1, nsp
           isa_s = isat + 1
           isa_e = isat + na(is) - 1
           IF( tranp(is) ) THEN
             WRITE( stdout,610) is, na(is)
             WRITE( stdout,615)
             DO isa = isa_s, isa_e
               oldp = tau(:,isa)
               CALL RANDOM_NUMBER( rand_disp )
               rand_disp = amprp(is) * ( rand_disp - 0.5d0 )
               rdisp     = rand_disp
               CALL r_to_s( rdisp(:), rand_disp(:), hinv )
               DO k = 1, 3
                 tau(k,isa) = tau(k,isa) + rand_disp(k) * ifor(k,isa)
               END DO
               WRITE( stdout,620) (oldp(k),k=1,3), (tau(k,isa),k=1,3)
             END DO
           END IF
           isat = isat + na(is)
         END DO

 600     FORMAT(//,3X,'Randomization of SCALED ionic coordinates')
 610     FORMAT(   3X,'Species ',I3,' atoms = ',I4)
 615     FORMAT(   3X,'     Old Positions               New Positions')
 620     FORMAT(   3X,3F10.6,2X,3F10.6)
       RETURN
       END SUBROUTINE randpos

!------------------------------------------------------------------------------!

  SUBROUTINE ions_kinene( ekinp, vels, na, nsp, h, pmass )
    IMPLICIT NONE
    real( kind=8 ), intent(out) :: ekinp     !  ionic kinetic energy
    real( kind=8 ), intent(in) :: vels(:,:)  !  scaled ionic velocities
    real( kind=8 ), intent(in) :: pmass(:)   !  ionic masses
    real( kind=8 ), intent(in) :: h(:,:)     !  simulation cell
    integer, intent(in) :: na(:), nsp
    integer :: i, j, is, ia, ii, isa
    ekinp = 0.0d0
    isa = 0
    do is=1,nsp
      do ia=1,na(is)
        isa = isa + 1
        do j=1,3
          do i=1,3
            do ii=1,3
              ekinp=ekinp+pmass(is)* h(j,i)*vels(i,isa)* h(j,ii)*vels(ii,isa)
            end do
          end do
        end do
      end do
    end do
    ekinp=0.5d0*ekinp
    return
  END SUBROUTINE

!------------------------------------------------------------------------------!

  subroutine ions_temp( tempp, temps, ekinpr, vels, na, nsp, h, pmass, ndega )
    use constants, only: factem
    implicit none
    real( kind=8 ), intent(out) :: ekinpr, tempp
    real( kind=8 ), intent(out) :: temps(:)
    real( kind=8 ), intent(in) :: vels(:,:)
    real( kind=8 ), intent(in) :: pmass(:)
    real( kind=8 ), intent(in) :: h(:,:)
    integer, intent(in) :: na(:), nsp, ndega
    integer :: nat, i, j, is, ia, ii, isa
    real( kind=8 ) :: cdmvel(3), eks
    call ions_cofmass( vels, pmass, na, nsp, cdmvel )
    nat = SUM( na(1:nsp) )
    ekinpr = 0.0d0
    temps( 1:nsp ) = 0.0d0
    do i=1,3
      do j=1,3
        do ii=1,3
          isa = 0
          do is=1,nsp
            eks = 0.0d0
            do ia=1,na(is)
              isa = isa + 1
              eks=eks+pmass(is)*h(j,i)*(vels(i,isa)-cdmvel(i))*h(j,ii)*(vels(ii,isa)-cdmvel(ii))
            end do
            ekinpr    = ekinpr    + eks
            temps(is) = temps(is) + eks
          end do
        end do
      end do
    end do
    do is = 1, nsp
      temps( is ) = temps( is ) * 0.5d0
      temps( is ) = temps( is ) * factem / ( 1.5d0 * na(is) )
    end do
    ekinpr = 0.5 * ekinpr
    tempp  = ekinpr * factem * 2.0d0 / DBLE( ndega )
    return
  end subroutine

!------------------------------------------------------------------------------!

  subroutine ions_thermal_stress( stress, pmass, omega, h, vels, nsp, na )
    real(kind=8), intent(inout) :: stress(3,3)
    real(kind=8), intent(in)  :: pmass(:), omega, h(3,3), vels(:,:)
    integer, intent(in) :: nsp, na(:)
    integer :: i, j, is, ia, isa
    isa    = 0
    do is = 1, nsp
      do ia = 1, na(is)
        isa = isa + 1
        do i = 1, 3
          do j = 1, 3
            stress(i,j) = stress(i,j) + pmass(is) / omega *           &
     &        (  (h(i,1)*vels(1,isa)+h(i,2)*vels(2,isa)+h(i,3)*vels(3,isa)) *  &
                 (h(j,1)*vels(1,isa)+h(j,2)*vels(2,isa)+h(j,3)*vels(3,isa))  )
          enddo
        enddo
      enddo
    enddo
    return
  end subroutine

!------------------------------------------------------------------------------!

  subroutine ions_vrescal( tcap, tempw, tempp, taup, tau0, taum, na, nsp, fion, iforce, &
                           pmass, delt )
    use constants, only: pi, factem
    implicit none
    logical, intent(in) :: tcap
    real(kind=8), intent(inout) :: taup(:,:)
    real(kind=8), intent(in) :: tau0(:,:), taum(:,:), fion(:,:)
    real(kind=8), intent(in) :: delt, pmass(:), tempw, tempp
    integer, intent(in) :: na(:), nsp
    integer, intent(in) :: iforce(:,:)

    real(kind=8) :: alfap, qr(3), alfar, gausp
    real(kind=8) :: dt2by2, ftmp
    real(kind=8) :: randy
    integer :: i, ia, is, nat, isa

    dt2by2 = .5d0 * delt * delt
    gausp = delt * sqrt( tempw / factem )
    nat = SUM( na( 1:nsp ) )

    if(.not.tcap) then
      alfap=.5d0*sqrt(tempw/tempp)
      isa = 0
      do is=1,nsp
        do ia=1,na(is)
          isa = isa + 1
          do i=1,3
            taup(i,isa) = tau0(i,isa) +                 &
     &                      alfap*(taup(i,isa)-taum(i,isa)) +      &
     &                      dt2by2/pmass(is)*fion(i,isa)*iforce(i,isa)
          end do
        end do
      end do
    else
      do i=1,3
        qr(i)=0.d0
        isa = 0
        do is=1,nsp
          do ia=1,na(is)
            isa = isa + 1
            alfar=gausp/sqrt(pmass(is))*cos(2.d0*pi*randy())*sqrt(-2.d0*log(randy()))
            taup(i,isa)=alfar
            qr(i)=qr(i)+alfar
          end do
        end do
        qr(i)=qr(i)/nat
      end do
      isa = 0
      do is=1,nsp
        do ia=1,na(is)
          isa = isa + 1
          do i=1,3
            alfar=taup(i,isa)-qr(i)
            taup(i,isa)=tau0(i,isa)+iforce(i,isa)*     &
     &                    (alfar+dt2by2/pmass(is)*fion(i,isa))
          end do
        end do
      end do
    end if
    return
  end subroutine

!------------------------------------------------------------------------------!

  subroutine ions_shiftvar( varp, var0, varm )
    implicit none
    real(kind=8), intent(in) :: varp(:,:)
    real(kind=8), intent(out) :: varm(:,:), var0(:,:)
      varm = var0
      var0 = varp
    return
  end subroutine


!------------------------------------------------------------------------------!
  END MODULE ions_base
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
  MODULE ions_positions
!------------------------------------------------------------------------------!

      USE kinds, ONLY: dbl
      USE parameters, ONLY: natx
!
      IMPLICIT NONE

! ... Atomic positions arrays used in the cp codes during the dynamic
      REAL(dbl) :: tau0(3,natx), taum(3,natx),  taup(3,natx)
      REAL(dbl) :: taus(3,natx), tausm(3,natx), tausp(3,natx)
      REAL(dbl) :: vels(3,natx), velsm(3,natx), velsp(3,natx)

!------------------------------------------------------------------------------!
  CONTAINS 
!------------------------------------------------------------------------------!

  subroutine ions_hmove( taus, tausm, iforce, pmass, fion, ainv, delt, na, nsp )
    real(kind=8), intent(in) :: tausm(:,:), pmass(:), fion(:,:)
    integer, intent(in) :: iforce(:,:)
    real(kind=8), intent(in) :: ainv(3,3), delt
    real(kind=8), intent(out) :: taus(:,:) 
    integer, intent(in) :: na(:), nsp
    integer :: is, ia, i, isa
    real(kind=8) :: dt2by2, fac, fions(3)

    dt2by2 = .5d0 * delt * delt

    isa = 0
    do is=1,nsp
      fac = dt2by2/pmass(is)
      do ia=1,na(is)
        isa = isa + 1
        do i=1,3
          fions( i ) = fion(1,isa)*ainv(i,1) + fion(2,isa)*ainv(i,2) + fion(3,isa)*ainv(i,3)
        end do
        do i=1,3
          taus(i,isa) = tausm(i,isa) + iforce(i,isa) * fac * fions( i )
        end do
      end do
    end do
    return
  end subroutine


  subroutine ions_move( tausp, taus, tausm, iforce, pmass, fion, ainv, delt, na, nsp, &
                        fricp, hgamma, vels, tsdp, tnosep, fionm, vnhp, velsp, velsm )
    implicit none
    real(kind=8), intent(in) :: taus(:,:), tausm(:,:), pmass(:), fion(:,:)
    integer, intent(in) :: iforce(:,:)
    real(kind=8), intent(in) :: ainv(3,3), delt
    real(kind=8), intent(out) :: tausp(:,:)
    integer, intent(in) :: na(:), nsp
    real(kind=8), intent(in) :: fricp, hgamma(3,3), vels(:,:)
    logical, intent(in) :: tsdp, tnosep
    real(kind=8), intent(inout) :: fionm(:,:)
    real(kind=8), intent(in) :: vnhp(:)
    real(kind=8), intent(out) :: velsp(:,:)
    real(kind=8), intent(in) :: velsm(:,:)
    integer :: is, ia, i, isa
    real(kind=8) :: dt2by2, fac, fions(3), dt2, twodel
    real(kind=8) :: verl1, verl2, verl3

    dt2by2 = .5d0 * delt * delt
    dt2    = delt * delt
    twodel = 2.0d0 * delt

         verl1=2./(1.+fricp)
         verl2=1.-verl1
         verl3=dt2/(1.+fricp)
!
         if(tsdp) then
            isa = 0
            do is=1,nsp
               do ia=1,na(is)
                  isa = isa + 1
                  do i=1,3
                     tausp(i,isa) = taus(i,isa) +                   &
     &                    iforce(i,isa)*dt2by2/pmass(is)*             &
     &        (ainv(i,1)*fion(1,isa)+ainv(i,2)*fion(2,isa)+         &
     &         ainv(i,3)*fion(3,isa) ) -                              &
     &                    pmass(is)*(hgamma(i,1)*vels(1,isa)+         &
     &         hgamma(i,2)*vels(2,isa)+hgamma(i,3)*vels(3,isa))
                  end do
               end do
            end do
         else if (tnosep) then
            isa = 0
            do is=1,nsp
               do ia=1,na(is)
                  isa = isa + 1
                  do i=1,3
                     fionm(i,isa) = (ainv(i,1)*fion(1,isa)          &
     &                                +ainv(i,2)*fion(2,isa)          &
     &                                +ainv(i,3)*fion(3,isa))         &
     &                              - vnhp(1)*vels(i,isa)*pmass(is)      &
     &                    - pmass(is)*(hgamma(i,1)*vels(1,isa)        &
     &                                +hgamma(i,2)*vels(2,isa)        &
     &                                +hgamma(i,3)*vels(3,isa))
                     tausp(i,isa)=-tausm(i,isa)+2.*taus(i,isa)+   &
     &                   iforce(i,isa)*dt2*fionm(i,isa)/pmass(is)
                     velsp(i,isa) = velsm(i,isa) +                  &
     &                    twodel*fionm(i,isa)/pmass(is)
                  end do
               end do
            end do
         else
            isa = 0
            do is=1,nsp
               do ia=1,na(is)
                  isa = isa + 1
                  do i=1,3
                     tausp(i,isa) = verl1*taus(i,isa)               &
     &                    + verl2*tausm(i,isa)                        &
     &        + verl3/pmass(is)*iforce(i,isa) * (ainv(i,1)*fion(1,isa)&
     &        + ainv(i,2)*fion(2,isa) + ainv(i,3)*fion(3,isa))      &
     &        - verl3*iforce(i,isa) * (hgamma(i,1)*vels(1,isa)      &
     &        + hgamma(i,2)*vels(2,isa) + hgamma(i,3)*vels(3,isa))
                     velsp(i,isa)=velsm(i,isa)                      &
     &        - 4.*fricp*vels(i,isa)                                  &
     &        + twodel/pmass(is)*iforce(i,isa)*(ainv(i,1)*fion(1,isa) &
     &        + ainv(i,2)*fion(2,isa) + ainv(i,3)*fion(3,isa))      &
     &        - twodel*iforce(i,isa) * (hgamma(i,1)*vels(1,isa)     &
     &        + hgamma(i,2)*vels(2,isa) + hgamma(i,3)*vels(3,isa))
                  end do
               end do
            end do
         endif
    return
  end subroutine



!------------------------------------------------------------------------------!
  END MODULE ions_positions
!------------------------------------------------------------------------------!


!------------------------------------------------------------------------------!
  MODULE ions_nose
!------------------------------------------------------------------------------!

      USE kinds, ONLY: dbl
      USE parameters, ONLY: nhclm
!
      IMPLICIT NONE

      INTEGER   :: nhpcl, ndega
      REAL(dbl) :: vnhp( nhclm ) = 0.0d0
      REAL(dbl) :: xnhp0( nhclm ) = 0.0d0
      REAL(dbl) :: xnhpm( nhclm ) = 0.0d0
      REAL(dbl) :: xnhpp( nhclm ) = 0.0d0
      REAL(dbl) :: xnhpm2( nhclm ) = 0.0d0
      REAL(dbl) :: qnp( nhclm ) = 0.0d0
      REAL(dbl) :: gkbt = 0.0d0
      REAL(dbl) :: kbt = 0.0d0
      REAL(dbl) :: tempw = 0.0d0
      REAL(dbl) :: fnosep( nhclm ) = 0.0d0

!------------------------------------------------------------------------------!
  CONTAINS 
!------------------------------------------------------------------------------!

  subroutine ions_nose_init( tempw_ , fnosep_ , nhpcl_ , ndega_ , nat )
    use constants,      only: factem, pi, terahertz
    use control_flags,  only: tnosep
    use ions_base,      only: ndofp, tions_base_init
    implicit none
    real(kind=8), intent(in)  :: tempw_ , fnosep_(:) 
    integer, intent(in) :: nhpcl_ , ndega_
    integer, intent(in) :: nat
    integer :: nsvar, i

    IF( .NOT. tions_base_init ) &
      CALL errore(' ions_nose_init ', ' you should call ions_base_init first ', 1 )

    vnhp  = 0.0d0
    xnhp0 = 0.0d0
    xnhpm = 0.0d0 
    xnhpm2 = 0.0d0 
    xnhpp = 0.0d0

    tempw     = tempw_
    fnosep    = 0.0d0

    qnp   = 0.0d0
    nhpcl = MAX( nhpcl_ , 1 )

    IF( nhpcl > nhclm ) &
      CALL errore(' ions_nose_init ', ' nhpcl out of range ', nhpcl )

    !  Setup Nose-Hoover chain masses
    !
    if ( ndega_ > 0 ) then
       ndega = ndega_
    else if ( ndega_ < 0 ) then
       ndega = ndofp - ( - ndega_ )
    else
       ndega = ndofp
    endif

    !  there is no need to initialize any Nose variables except for nhpcl
    !  and ndega if the thermostat is not used
    !

    IF( tnosep ) THEN

      IF( nhpcl > SIZE( fnosep_ ) ) &
        CALL errore(' ions_nose_init ', ' fnosep size too small ', nhpcl )

      gkbt = DBLE( ndega ) * tempw / factem
      kbt  = tempw / factem

      fnosep(1) = fnosep_ (1)
      if( fnosep(1) > 0.0d0 ) then
        qnp(1) = 2.d0 * gkbt / ( fnosep(1) * ( 2.d0 * pi ) * terahertz )**2
      end if

      if ( nhpcl > 1 ) then
        do i = 2, nhpcl
          fnosep(i) = fnosep_ (i)
          if( fnosep(i) > 0.0d0 ) then
            qnp(i) = 2.d0 * tempw / factem / ( fnosep(i) * ( 2.d0 * pi ) * terahertz )**2
          else
            qnp(i) = qnp(1) / dble(ndega)
          endif
        enddo
      endif

    END IF


    !    WRITE( stdout,100)
    !    WRITE( stdout,110) QNOSEP,TEMPW
    !    WRITE( stdout,120) GLIB
    !    WRITE( stdout,130) NSVAR
 100  FORMAT(//' * Temperature control of ions with nose thermostat'/)
 110  FORMAT(3X,'nose mass:',F12.4,' temperature (K):',F12.4)
 120  FORMAT(3X,'ionic degrees of freedom:        ',F5.0)
 130  FORMAT(3X,'time steps per nose oscillation: ',I5,//)

    return
  end subroutine


  SUBROUTINE ions_nose_info()

      use constants,     only: factem, terahertz, pi
      use time_step,     only: delt
      USE io_global,     ONLY: stdout
      USE control_flags, ONLY: tnosep

      IMPLICIT NONE

      INTEGER   :: nsvar, i
      REAL(dbl) :: wnosep

      IF( tnosep ) THEN
        !
        IF( fnosep(1) <= 0.D0) &
          CALL errore(' ions_nose_info ', ' fnosep less than zero ', 1)
        IF( delt <= 0.D0) &
          CALL errore(' ions_nose_info ', ' delt less than zero ', 1)

        wnosep = fnosep(1) * ( 2.d0 * pi ) * terahertz
        nsvar  = ( 2.d0 * pi ) / ( wnosep * delt )

        WRITE( stdout,563) tempw, nhpcl, ndega, (qnp(i),i=1,nhpcl), nsvar
      END IF

 563  format( //, &
            & 3X,'ion dynamics with nose` temperature control:', /, &
            & 3X,'temperature required      = ', f10.5, ' (kelvin) ', /, &
            & 3X,'NH chain length           = ', i3, /, &
            & 3X,'active degrees of freedom = ', i3, /, &
            & 3X,'time steps per nose osc.  = ', i5, /, &
            & 3X,'nose` mass(es)            = ', 20(1X,f10.3),//) 


    RETURN
  END SUBROUTINE ions_nose_info

  

  subroutine ions_nosevel( vnhp, xnhp0, xnhpm, delt, nhpcl )
    implicit none
    real(kind=8), intent(inout) :: vnhp(:)
    real(kind=8), intent(in) :: xnhp0(:), xnhpm(:), delt
    integer, intent(in) :: nhpcl
    integer :: i
    do i=1,nhpcl
       vnhp(i)=2.*(xnhp0(i)-xnhpm(i))/delt-vnhp(i)
    enddo
    return
  end subroutine


 subroutine ions_noseupd( xnhpp, xnhp0, xnhpm, delt, qnp, ekinpr, gkbt, vnhp, kbt, nhpcl )
    implicit none
    real(kind=8), intent(out) :: xnhpp(:), vnhp(:)
    real(kind=8), intent(in) :: xnhp0(:), xnhpm(:), delt, qnp(:), ekinpr, gkbt, kbt
    integer, intent(in) :: nhpcl
    integer :: i
    real(kind=8) :: dt2, zetfrc

    zetfrc = 2.0d0*ekinpr-gkbt
    dt2 = delt**2
    If (nhpcl.gt.1) then
       do i=1,(nhpcl-1)
          xnhpp(i)=(4.d0*xnhp0(i)-(2.d0-delt*vnhp(i+1))*xnhpm(i)+2.0d0*dt2*zetfrc/qnp(i))&
               &   /(2.d0+delt*vnhp(i+1))
          vnhp(i) =(xnhpp(i)-xnhpm(i))/( 2.0d0 * delt )
          zetfrc = (qnp(i)*vnhp(i)**2-kbt)
       enddo
    endif
    ! Last variable
    xnhpp(nhpcl)=2.d0*xnhp0(nhpcl)-xnhpm(nhpcl)+( delt**2 / qnp(nhpcl) )*zetfrc
    i = nhpcl
    vnhp(i) =(xnhpp(i)-xnhpm(i))/( 2.0d0 * delt )
    ! Update velocities
!    do i=1,nhpcl
!       vnhp(i) =(xnhpp(i)-xnhpm(i))/( 2.0d0 * delt )
!    end do
    ! These are the original expressions from cpr.f90
    !       xnhpp(1)=2.*xnhp0(1)-xnhpm(1)+2.*( delt**2 / qnp(1) )*(ekinpr-gkbt/2.)
    !       vnhp(1) =(xnhpp(1)-xnhpm(1))/( 2.0d0 * delt )
    return
  end subroutine ions_noseupd

  

  real(kind=8) function ions_nose_nrg( xnhp0, vnhp, qnp, gkbt, kbt, nhpcl )
    implicit none
    integer :: nhpcl
    real(kind=8) :: gkbt,qnp(:),vnhp(:),xnhp0(:),kbt
    integer :: i
    real(kind=8) :: stmp
    !
    stmp = 0.5d0 * qnp(1) * vnhp(1) * vnhp(1) +     gkbt * xnhp0(1)
    if (nhpcl > 1) then
       do i=2,nhpcl
          stmp=stmp+0.5*qnp(i)*vnhp(i)*vnhp(i) + kbt*xnhp0(i)
       enddo
    endif
    ions_nose_nrg = stmp
    return
  end function ions_nose_nrg



  subroutine ions_nose_shiftvar( xnhpp, xnhp0, xnhpm, xnhpm2 )
    !  shift values of nose variables to start a new step
    implicit none
    real(kind=8), optional, intent(out) :: xnhpm2(:)
    real(kind=8), intent(inout) :: xnhpm(:), xnhp0(:)
    real(kind=8), intent(in)  :: xnhpp(:)
      !
      IF( PRESENT( xnhpm2 ) ) xnhpm2 = xnhpm
      xnhpm = xnhp0
      xnhp0 = xnhpp
      !
    return
  end subroutine
  

!------------------------------------------------------------------------------!
  END MODULE ions_nose
!------------------------------------------------------------------------------!
