!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
  MODULE ions_base
!------------------------------------------------------------------------------!

      USE kinds,      ONLY : DP
      USE parameters, ONLY : ntypx
      !
      IMPLICIT NONE
      SAVE

      !     nsp       = number of species
      !     na(is)    = number of atoms of species is
      !     nax       = max number of atoms of a given species
      !     nat       = total number of atoms of all species

      INTEGER :: nsp     = 0
      INTEGER :: na(ntypx) = 0
      INTEGER :: nax     = 0
      INTEGER :: nat     = 0

      !     zv(is)    = (pseudo-)atomic charge
      !     amass(is) = mass of ions, in atomic mass units
      !     rcmax(is) = Ewald radius (for ion-ion interactions)

      REAL(DP) :: zv(ntypx)    = 0.0_DP
      REAL(DP) :: amass(ntypx) = 0.0_DP
      REAL(DP) :: rcmax(ntypx) = 0.0_DP

      !     ityp( i ) = the type of i-th atom in stdin
      !     atm( j )  = name of the type of the j-th atomic specie
      !     tau( 1:3, i ) = position of the i-th atom

      INTEGER,  ALLOCATABLE :: ityp(:)
      REAL(DP), ALLOCATABLE :: tau(:,:)     !  initial positions read from stdin (in bohr)
      REAL(DP), ALLOCATABLE :: vel(:,:)     !  initial velocities read from stdin (in bohr)
      CHARACTER(LEN=3)      :: atm( ntypx )
      CHARACTER(LEN=80)     :: tau_format   ! format of input atomic positions:
                                            ! 'alat','crystal','bohr','angstrom'

      ! if_pos( x, i ) = 0 : x coordinate of i-th atom will be kept fixed
      INTEGER, ALLOCATABLE :: if_pos(:,:)  ! allowed values: 0 or 1 only
      INTEGER, ALLOCATABLE :: iforce(:,:)  ! if_pos sorted by specie 
      INTEGER :: fixatom   = 0            ! number of frozen atoms
      INTEGER :: ndofp     =-1            ! ionic degree of freedom
      INTEGER :: ndfrz     = 0            ! frozen degrees of freedom

      REAL(DP) :: fricp   ! friction parameter for damped dynamics
      REAL(DP) :: greasp  ! friction parameter for damped dynamics

      ! ... taui = real ionic positions in the center of mass reference
      ! ... system at istep = 0
      ! ... this array is used to compute mean square displacements,
      ! ... it is initialized when NBEG = -1, NBEG = 0 and TAURDR = .TRUE.
      ! ... first index: x,y,z, second index: atom sorted by specie with respect input
      ! ... this array is saved in the restart file

      REAL(DP), ALLOCATABLE :: taui(:,:)

      ! ... cdmi = center of mass reference system (related to the taui)
      ! ... this vector is computed when NBEG = -1, NBEG = 0 and TAURDR = .TRUE.
      ! ... this array is saved in the restart file

      REAL(DP) :: cdmi(3), cdm(3)

      ! ... cdms = center of mass computed for scaled positions (taus)

      REAL(DP) :: cdms(3)
      !
      REAL(DP), ALLOCATABLE :: extfor(:,:)     !  external forces on atoms

      LOGICAL :: tions_base_init = .FALSE.
      LOGICAL, PRIVATE :: tdebug = .FALSE.

      
!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!

    !-------------------------------------------------------------------------
    SUBROUTINE ions_base_init( nsp_, nat_, na_, ityp_, tau_, vel_, amass_,&
                               atm_, if_pos_, tau_format_, alat_, at_,    & 
                               rcmax_ , extfor_ )
      !-------------------------------------------------------------------------
      !
      USE constants, ONLY: bohr_radius_angs
      USE io_global, ONLY: stdout
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN) :: nsp_, nat_, na_(:), ityp_(:)
      REAL(DP),         INTENT(IN) :: tau_(:,:)
      REAL(DP),         INTENT(IN) :: vel_(:,:)
      REAL(DP),         INTENT(IN) :: amass_(:)
      CHARACTER(LEN=*), INTENT(IN) :: atm_(:)
      CHARACTER(LEN=*), INTENT(IN) :: tau_format_
      INTEGER,          INTENT(IN) :: if_pos_(:,:)
      REAL(DP),         INTENT(IN) :: alat_, at_(3,3)
      REAL(DP),         INTENT(IN) :: rcmax_(:)
      REAL(DP),         INTENT(IN) :: extfor_(:,:)
      !
      INTEGER :: i, ia, is
      !
      !
      nsp = nsp_
      nat = nat_
      !
      IF ( nat < 1 ) &
         CALL errore( 'ions_base_init ', 'nax out of range', 1 )
      IF ( nsp < 1 ) &
         CALL errore( 'ions_base_init ', 'nsp out of range', 1 )
      IF ( nsp > SIZE( na ) ) &
         CALL errore( 'ions_base_init ', &
                    & 'nsp too large, increase ntypx parameter ', 1 )
      !
      na(1:nsp) = na_(1:nsp)
      nax       = MAXVAL( na(1:nsp) )
      !
      atm(1:nsp) = atm_(1:nsp)
      tau_format = TRIM( tau_format_ )
      !
      IF ( nat /= SUM( na(1:nsp) ) ) &
         CALL errore( 'ions_base_init ','inconsistent nat and na ', 1 )
      !
      CALL deallocate_ions_base()
      !
      ALLOCATE( ityp( nat ) )
      ALLOCATE( tau( 3, nat ) )
      ALLOCATE( vel( 3, nat ) )
      ALLOCATE( if_pos( 3, nat ) )
      ALLOCATE( iforce( 3, nat ) )
      ALLOCATE( taui( 3, nat ) )
      ALLOCATE( extfor( 3, nat ) )
      !
      ityp(1:nat)     = ityp_(1:nat)
      vel(:,1:nat)    = vel_(:,1:nat)
      if_pos(:,1:nat) = if_pos_(:,1:nat)
      !
      ! ... radii, masses 
      !
      DO is = 1, nsp_
         !
         rcmax(is) = rcmax_(is)
         !
         IF( rcmax(is) <= 0.0_DP ) &
            CALL errore( 'ions_base_init ', 'invalid rcmax', is )
         !
      END DO
      !
      SELECT CASE ( TRIM( tau_format ) )
         !
         ! ... convert input atomic positions to internally used format:
         ! ... tau in atomic units
         !
         CASE( 'alat' )
            !
            ! ... input atomic positions are divided by a0
            !
            tau(:,1:nat) = tau_(:,1:nat) * alat_
            vel(:,1:nat) = vel_(:,1:nat) * alat_
            !
         CASE( 'bohr' )
            !
            ! ... input atomic positions are in a.u.: do nothing
            !
            tau(:,1:nat) = tau_(:,1:nat)
            vel(:,1:nat) = vel_(:,1:nat)
            !
         CASE( 'crystal' )
            !
            ! ... input atomic positions are in crystal axis ("scaled")
            !
            DO ia = 1, nat
               !
               DO i = 1, 3
                  !
                  tau(i,ia) = at_(i,1)*alat_ * tau_(1,ia) + &
                              at_(i,2)*alat_ * tau_(2,ia) + &
                              at_(i,3)*alat_ * tau_(3,ia)
                  !
                  vel(i,ia) = at_(i,1)*alat_ * vel_(1,ia) + &
                              at_(i,2)*alat_ * vel_(2,ia) + &
                              at_(i,3)*alat_ * vel_(3,ia)
               
               END DO
               !
            END DO
            !
         CASE( 'angstrom' )
            !
            ! ... atomic positions in A
            !
            tau(:,1:nat) = tau_(:,1:nat) / bohr_radius_angs
            vel(:,1:nat) = vel_(:,1:nat) / bohr_radius_angs
            !
         CASE DEFAULT
            !
            CALL errore( 'ions_base_init',' tau_format = ' // &
                       & TRIM( tau_format ) // ' not implemented ', 1 )
            !
      END SELECT
      !
      extfor( :, 1:nat ) = extfor_( :, 1:nat )
      !
      IF( tdebug ) THEN
        WRITE( stdout, * ) 'ions_base_init: unsorted position and velocities'
        DO ia = 1, nat
          WRITE( stdout, fmt="(A3,3D12.4,3X,3D12.4)") &
            atm( ityp( ia ) ), tau(1:3, ia), vel(1:3,ia)
        END DO
        WRITE( stdout, * ) 'ions_base_init: sorted position and velocities'
        DO ia = 1, nat
          WRITE( stdout, fmt="(A3,3D12.4,3X,3D12.4)") &
            atm( ityp( ia ) ), tau(1:3, ia), vel(1:3,ia)
        END DO
      END IF
      !
      ! ... The constrain on fixed coordinates is implemented using the array
      ! ... if_pos whose value is 0 when the coordinate is to be kept fixed, 1
      ! ... otherwise. 
      !
      if_pos = 1
      if_pos(:,:) = if_pos_(:,1:nat)
      !
      iforce = 0
      iforce(:,:) = if_pos(:,:)
      !
      fixatom=COUNT( if_pos(1,:)==0 .AND. if_pos(2,:)==0 .AND. if_pos(3,:)==0 )
      ndofp = COUNT( iforce == 1 )
      ndfrz = 3*nat - ndofp
      !
      amass(1:nsp) = amass_(1:nsp)
      !
      IF ( ANY( amass(1:nsp) <= 0.0_DP ) ) &
         CALL errore( 'ions_base_init ', 'invalid  mass', 1 ) 
      !
      CALL ions_cofmass( tau, amass, nat, ityp, cdmi )
      !
      DO ia = 1, nat
         !
         taui(1:3,ia) = tau(1:3,ia) - cdmi(1:3)
         !
      END DO
      !
      tions_base_init = .TRUE.
      !
      RETURN
      !
    END SUBROUTINE ions_base_init
    !
    !-------------------------------------------------------------------------
    SUBROUTINE deallocate_ions_base()
      !-------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      IF ( ALLOCATED( ityp ) )    DEALLOCATE( ityp )
      IF ( ALLOCATED( tau ) )     DEALLOCATE( tau )
      IF ( ALLOCATED( vel ) )     DEALLOCATE( vel )
      IF ( ALLOCATED( if_pos ) )  DEALLOCATE( if_pos )
      IF ( ALLOCATED( iforce ) )  DEALLOCATE( iforce )
      IF ( ALLOCATED( taui ) )    DEALLOCATE( taui )
      IF ( ALLOCATED( extfor ) )  DEALLOCATE( extfor )
      !
      tions_base_init = .FALSE.
      !
      RETURN
      !
    END SUBROUTINE deallocate_ions_base
    !
    !-------------------------------------------------------------------------
    SUBROUTINE ions_vel( vel, taup, taum, dt )
      !-------------------------------------------------------------------------
      USE constants, ONLY : eps8
      IMPLICIT NONE
      REAL(DP) :: vel(:,:), taup(:,:), taum(:,:)
      REAL(DP) :: dt
      REAL(DP) :: fac
      IF( dt < eps8 ) &
         CALL errore( ' ions_vel ', ' dt <= 0 ', 1 )
      fac  = 1.0_DP / ( dt * 2.0_DP )
      vel(:,:) = ( taup(:,:) - taum(:,:) ) * fac
      RETURN
    END SUBROUTINE ions_vel

!------------------------------------------------------------------------------!

    SUBROUTINE ions_cofmass( tau, pmass, nat, ityp, cdm )
      USE constants, ONLY : eps8
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: tau(:,:), pmass(:)
      REAL(DP), INTENT(OUT) :: cdm(3)
      INTEGER, INTENT(IN) :: nat, ityp(:)

      REAL(DP) :: tmas
      INTEGER :: is, i, ia, isa
!
      tmas=0.0_DP
      cdm=0.0_DP
      do ia=1,nat
         do i=1,3
            cdm(i)=cdm(i)+tau(i,ia)*pmass(ityp(ia))
         end do
         tmas=tmas+pmass(ityp(ia))
      end do

      if( tmas < eps8 ) &
         call errore(' ions_cofmass ', ' total mass <= 0 ', 1 )
!
      cdm=cdm/tmas
!
      RETURN
    END SUBROUTINE ions_cofmass

!------------------------------------------------------------------------------!

      SUBROUTINE randpos(tau, nat, ityp, tranp, amprp, hinv, ifor )
        
         USE cell_base, ONLY: r_to_s
         USE io_global, ONLY: stdout
         USE random_numbers, ONLY: randy

         IMPLICIT NONE
         REAL(DP) :: hinv(3,3)
         REAL(DP) :: tau(:,:)
         INTEGER, INTENT(IN) :: ifor(:,:), nat, ityp(:)
         LOGICAL, INTENT(IN) :: tranp(:)
         REAL(DP), INTENT(IN) :: amprp(:)
         REAL(DP) :: oldp(3), rand_disp(3), rdisp(3)
         INTEGER :: k, ia

         WRITE( stdout, 600 )
         WRITE( stdout,615)
         DO ia = 1, nat
           IF( tranp(ityp(ia)) ) THEN
             oldp = tau(:,ia)
             rand_disp(1) = randy ()
             rand_disp(2) = randy ()
             rand_disp(3) = randy ()
             rand_disp = amprp(ityp(ia)) * ( rand_disp - 0.5_DP )
             rdisp     = rand_disp
             CALL r_to_s( rdisp(:), rand_disp(:), hinv )
             DO k = 1, 3
                tau(k,ia) = tau(k,ia) + rand_disp(k) * ifor(k,ia)
             END DO
             WRITE( stdout,620) (oldp(k),k=1,3), (tau(k,ia),k=1,3)
           END IF
         END DO

 600     FORMAT(//,3X,'Randomization of SCALED ionic coordinates')
 615     FORMAT(   3X,'     Old Positions               New Positions')
 620     FORMAT(   3X,3F10.6,2X,3F10.6)
       RETURN
       END SUBROUTINE randpos

!------------------------------------------------------------------------------!

  SUBROUTINE ions_kinene( ekinp, vels, nat, ityp, h, pmass )
    IMPLICIT NONE
    REAL(DP), intent(out) :: ekinp     !  ionic kinetic energy
    REAL(DP), intent(in) :: vels(:,:)  !  scaled ionic velocities
    REAL(DP), intent(in) :: pmass(:)   !  ionic masses
    REAL(DP), intent(in) :: h(:,:)     !  simulation cell
    integer, intent(in) :: nat, ityp(:)
    integer :: i, j, ia, ii
    ekinp = 0.0_DP
    do ia=1,nat
      do j=1,3
        do i=1,3
          do ii=1,3
            ekinp=ekinp+pmass(ityp(ia))* h(j,i)*vels(i,ia)* h(j,ii)*vels(ii,ia)
          end do
        end do
      end do
    end do
    ekinp=0.5_DP*ekinp
    return
  END SUBROUTINE ions_kinene

!------------------------------------------------------------------------------!

  subroutine ions_temp( tempp, temps, ekinpr, vels, nsp, na, nat, ityp, h, pmass, ndega, nhpdim, atm2nhp, ekin2nhp )
    !
    use constants, only: k_boltzmann_au
    !
    implicit none
    !
    REAL(DP), intent(out) :: ekinpr, tempp
    REAL(DP), intent(out) :: temps(:)
    REAL(DP), intent(out) :: ekin2nhp(:)
    REAL(DP), intent(in)  :: vels(:,:)
    REAL(DP), intent(in)  :: pmass(:)
    REAL(DP), intent(in)  :: h(:,:)
    integer,        intent(in)  :: nsp, na(:), nat, ityp(:), ndega, nhpdim, atm2nhp(:)
    !
    integer        :: i, j, ia, ii
    REAL(DP) :: cdmvel(3), eks1
    !
    call ions_cofmass( vels, pmass, nat, ityp, cdmvel )
    !
    ekinpr             = 0.0_DP
    temps( 1:nsp )     = 0.0_DP
    ekin2nhp(1:nhpdim) = 0.0_DP
    !
    do i=1,3
      do j=1,3
        do ii=1,3
          do ia=1,nat
            eks1 = pmass(ityp(ia))*h(j,i)*(vels(i,ia)-cdmvel(i))*h(j,ii)*(vels(ii,ia)-cdmvel(ii))
            ekin2nhp(atm2nhp(ia)) = ekin2nhp(atm2nhp(ia)) + eks1
            ekinpr    = ekinpr    + eks1
            temps(ityp(ia)) = temps(ityp(ia)) + eks1
          end do
        end do
      end do
    end do
    !
    ekin2nhp(1:nhpdim) = ekin2nhp(1:nhpdim) * 0.5_DP
    !
    !
    temps( 1:nsp ) = 0.5_DP * temps( 1:nsp ) / k_boltzmann_au / ( 1.5_DP * na(1:nsp) )
    !
    ekinpr = 0.5_DP * ekinpr
    !
    IF( ndega < 1 ) THEN 
       tempp = 0.0_DP
    ELSE
       tempp  = ekinpr / k_boltzmann_au * 2.0_DP / DBLE( ndega )
    END IF
    !
    return
  end subroutine ions_temp

!------------------------------------------------------------------------------!

  subroutine ions_thermal_stress( stress, nstress, pmass, omega, h, vels, nat, ityp )
    USE constants, ONLY : eps8
    REAL(DP), intent(inout) :: stress(3,3),nstress(3,3)
    REAL(DP), intent(in)  :: pmass(:), omega, h(3,3), vels(:,:)
    integer, intent(in) :: nat, ityp(:)
    integer :: i, j, ia

    nstress = 0.0_DP ! BS
    if( omega < eps8 ) call errore(' ions_thermal_stress ', ' omega <= 0 ', 1 )
    do ia = 1, nat
      do i = 1, 3
        do j = 1, 3
            stress(i,j) = stress(i,j) + pmass(ityp(ia)) / omega *           &
     &        (  (h(i,1)*vels(1,ia)+h(i,2)*vels(2,ia)+h(i,3)*vels(3,ia)) *  &
                 (h(j,1)*vels(1,ia)+h(j,2)*vels(2,ia)+h(j,3)*vels(3,ia))  )

     ! BS : to print Pressure of Nuclei at each step
            nstress(i,j) = nstress(i,j) + pmass(ityp(ia)) / omega *           &
     &        (  (h(i,1)*vels(1,ia)+h(i,2)*vels(2,ia)+h(i,3)*vels(3,ia)) *  &
                 (h(j,1)*vels(1,ia)+h(j,2)*vels(2,ia)+h(j,3)*vels(3,ia))  )
     ! BS

        enddo
      enddo
    enddo

    return
  end subroutine ions_thermal_stress

!------------------------------------------------------------------------------!

  subroutine randvel( tempw, tau0, taum, nat, ityp, iforce, amass, delt )
    use constants, only: pi, k_boltzmann_au, amu_au
    USE random_numbers, ONLY : randy
    implicit none
    REAL(DP), intent(out) :: taum(:,:)
    REAL(DP), intent(in) :: tau0(:,:)
    REAL(DP), intent(in) :: delt, amass(:), tempw
    integer, intent(in) :: nat, ityp(:)
    integer, intent(in) :: iforce(:,:)

    REAL(DP) :: alfap, qr(3), alfar, gausp
    REAL(DP) :: dt2by2
    integer :: i, ia

    dt2by2 = 0.5_DP * delt * delt
    gausp = delt * sqrt( tempw * k_boltzmann_au )

      ! randomize velocities, calculate center of mass vel.
      do i=1,3
        qr(i)=0.0_DP
        do ia=1,nat
          alfar=gausp/sqrt(amass(ityp(ia))*amu_au)*cos(2.0_DP*pi*randy())*sqrt(-2.0_DP*log(randy()))
          taum(i,ia)=alfar
          qr(i)=qr(i)+alfar
        end do
        qr(i)=qr(i)/nat
      end do
      !subtract center of mass velocity
      do ia=1,nat
        do i=1,3
            alfar=taum(i,ia)-qr(i)
            taum(i,ia)=tau0(i,ia)-iforce(i,ia)*alfar
        end do
      end do
    return
  end subroutine randvel

!------------------------------------------------------------------------------!

  subroutine ions_vrescal( tcap, tempw, tempp, taup, tau0, taum, nat, ityp, fion, iforce, pmass, delt )
    use constants, only: pi, k_boltzmann_au, eps8
    USE random_numbers, ONLY : randy
    implicit none
    logical, intent(in) :: tcap
    REAL(DP), intent(inout) :: taup(:,:)
    REAL(DP), intent(in) :: tau0(:,:), taum(:,:), fion(:,:)
    REAL(DP), intent(in) :: delt, pmass(:), tempw, tempp
    integer, intent(in) :: nat, ityp(:)
    integer, intent(in) :: iforce(:,:)

    REAL(DP) :: alfap, qr(3), alfar, gausp
    REAL(DP) :: dt2by2
    integer :: i, ia

    dt2by2 = 0.5_DP * delt * delt
    gausp = delt * sqrt( tempw * k_boltzmann_au )

    if(.not.tcap) then
      if( tempp < eps8 ) call errore(' ions_vrescal ', ' tempp <= 0 ', 1 )
      alfap = 0.5_DP * sqrt(tempw/tempp)
      do ia=1,nat
        do i=1,3
            taup(i,ia) = tau0(i,ia) +                 &
     &                      alfap*(taup(i,ia)-taum(i,ia)) +      &
     &                      dt2by2/pmass(ityp(ia))*fion(i,ia)*iforce(i,ia)
        end do
      end do
    else
      ! randomize velocities, calculate center of mass vel.
      do i=1,3
        qr(i)=0.0_DP
        do ia=1,nat
           alfar=gausp/sqrt(pmass(ityp(ia)))*cos(2.0_DP*pi*randy())*sqrt(-2.0_DP*log(randy()))
           taup(i,ia)=alfar
           qr(i)=qr(i)+alfar
        end do
        qr(i)=qr(i)/nat
      end do
      !subtract center of mass velocity
      do ia=1,nat
          do i=1,3
            alfar=taup(i,ia)-qr(i)
            taup(i,ia)=tau0(i,ia)+iforce(i,ia)*(alfar+dt2by2/pmass(ityp(ia))*fion(i,ia))
          end do
      end do
    end if
    return
  end subroutine ions_vrescal

!------------------------------------------------------------------------------!

  subroutine ions_shiftvar( varp, var0, varm )
    implicit none
    REAL(DP), intent(in) :: varp(:,:)
    REAL(DP), intent(out) :: varm(:,:), var0(:,:)
      varm = var0
      var0 = varp
    return
  end subroutine ions_shiftvar

!------------------------------------------------------------------------------!

  SUBROUTINE ions_reference_positions( tau )

     !  Calculate the real position of atoms relative to the center of mass (cdm)
     !  and store them in taui
     !    cdmi: initial position of the center of mass (cdm) in cartesian coor.  

     IMPLICIT NONE

     REAL(DP) :: tau( :, : )

     INTEGER  :: ia

     CALL ions_cofmass( tau, amass, nat, ityp, cdmi )
     DO ia = 1, nat
        taui(:,ia) = tau(:,ia) - cdmi(:)
     END DO

     RETURN
  END SUBROUTINE ions_reference_positions


!------------------------------------------------------------------------------!


   SUBROUTINE ions_displacement( dis, tau, nsp, nat, ityp )

      !  Calculate the sum of the quadratic displacements of the atoms in the ref.
      !    of cdm respect to the initial positions.
      !    taui: initial positions in real units in the ref. of cdm

      !  ----------------------------------------------
      !  att!     tau_ref: starting position in center-of-mass ref. in real units
      !  ----------------------------------------------

      IMPLICIT NONE

      REAL (DP), INTENT(OUT) :: dis(:)
      REAL (DP), INTENT(IN)  :: tau(:,:)
      INTEGER, INTENT(IN)  :: nsp, nat, ityp(:)

      REAL(DP) :: rdist(3), cdm(3)
      INTEGER  :: ia

      ! ...   Compute the current value of cdm "Center of Mass"
      !
      CALL ions_cofmass(tau, amass, nat, ityp, cdm )
      !
      dis  = 0.0_DP
      DO ia = 1, nat
         rdist = tau(:,ia) - cdm
         dis(ityp(ia)) = dis(ityp(ia)) + SUM( ( rdist(:) - taui(:,ia) )**2 )
      END DO
      dis(1:nsp) = dis(1:nsp) / DBLE(na(1:nsp))

      RETURN
   END SUBROUTINE ions_displacement

  !--------------------------------------------------------------------------
  SUBROUTINE ions_cofmsub( tausp, iforce, nat, cdm, cdm0 )
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: tausp(:,:)
    INTEGER,  INTENT(IN)    :: iforce(:,:)
    INTEGER,  INTENT(IN)    :: nat
    REAL(DP), INTENT(IN)    :: cdm(:), cdm0(:)
    !
    INTEGER :: i, ia
    !
    DO ia = 1, nat
       !
       DO i = 1, 3
          !
          tausp(i,ia) = tausp(i,ia) + DBLE( iforce(i,ia) ) * ( cdm0(i) - cdm(i) )
          !
       END DO
       !
    END DO
    !
    RETURN
    !
  END SUBROUTINE ions_cofmsub


  REAL(DP) FUNCTION compute_eextfor( tau0 )
     IMPLICIT NONE
     REAL(DP), OPTIONAL, INTENT(IN) :: tau0(:,:)
     INTEGER :: i
     REAL(DP) :: e
     compute_eextfor = 0.0d0
     e = 0.0d0
     IF( PRESENT( tau0 ) ) THEN
        DO i = 1, SIZE( extfor,2 )
          e = e + extfor( 3, i ) *  tau0( 3, i ) &
                + extfor( 2, i ) *  tau0( 2, i ) &
                + extfor( 1, i ) *  tau0( 1, i )
        END DO
     ELSE
        DO i = 1, SIZE( extfor,2 )
          e = e + extfor( 3, i ) *  tau( 3, i ) &
                + extfor( 2, i ) *  tau( 2, i ) &
                + extfor( 1, i ) *  tau( 1, i )
        END DO
     END IF
     compute_eextfor = - e
     RETURN
  END FUNCTION compute_eextfor
  


!------------------------------------------------------------------------------!
  END MODULE ions_base
!------------------------------------------------------------------------------!
