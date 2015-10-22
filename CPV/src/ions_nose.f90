!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
  MODULE ions_nose
!------------------------------------------------------------------------------!

      USE kinds, ONLY: DP
!
      IMPLICIT NONE
! Some comments are in order on how Nose-Hoover chains work here (K.N. Kudin)
! the present code allows one to use "massive" Nose-Hoover chains:
! TOBIAS DJ, MARTYNA GJ, KLEIN ML
! JOURNAL OF PHYSICAL CHEMISTRY 97 (49): 12959-12966 DEC 9 1993
!
! one chain for the whole system is specified by nhptyp=0 (or nothing)
! currently input options allow one chain per atomic type (nhptyp=1),
! one chain per atom (nhptyp=2), and fancy stuff with nhptyp=3 (& nhgrp)
!
! nhpdim is the total number of the resulting NH chains
! nhpend is 1 if there is a chain above all chains, otherwise it is 0
! nhpbeg is usually 0, however, if using groups (nhptyp = 3), it might
! be desirable to have atoms with uncontrolled temperature, so then
! nhpbeg becomes 1, and all the "uncontrolled" atoms are assigned to the
! 1st thermostat that is always zero in velocity and so it does not
! affect ionic motion
! array atm2nhp(1:nat) gives the chain number from the atom list (which
! is sorted by type)
! anum2nhp is the number of degrees of freedom per chain (now just 3*nat_i)
! ekin2nhp is the kinetic energy of the present chain
! gkbt2nhp are the NH chain parameters
! qnp are the chain masses, qnp_ is a temporary array for now
! see subroutine ions_nose_allocate on what are the dimensions of these
! variables
!
      INTEGER   :: nhpcl=1, ndega, nhpdim=1, nhptyp=0, nhpbeg=0, nhpend=0
      INTEGER, ALLOCATABLE   :: atm2nhp(:)
      INTEGER, ALLOCATABLE   :: anum2nhp(:)
      REAL(DP), ALLOCATABLE :: vnhp(:), xnhp0(:), xnhpm(:), xnhpp(:), &
      ekin2nhp(:), gkbt2nhp(:), scal2nhp(:), qnp(:), qnp_(:), fnosep(:)

      REAL(DP) :: gkbt = 0.0_DP
      REAL(DP) :: kbt = 0.0_DP
      REAL(DP) :: tempw = 0.0_DP

!------------------------------------------------------------------------------!
  CONTAINS 
!------------------------------------------------------------------------------!

  subroutine ions_nose_init( tempw_ , fnosep_ , nhpcl_ , nhptyp_ , ndega_ , nhgrp_ , fnhscl_)
    use constants,      only: k_boltzmann_au, pi, au_terahertz
    use control_flags,  only: tnosep
    use ions_base,      only: ndofp, tions_base_init, nsp, nat, na, amass, ityp
    real(DP), intent(in)  :: tempw_ , fnosep_(:), fnhscl_(:) 
    integer, intent(in) :: nhpcl_ , nhptyp_ , ndega_ , nhgrp_(:)
    integer :: i, j, iat, is, ia
    REAL(DP) :: amass_mean

    IF( .NOT. tions_base_init ) &
      CALL errore(' ions_nose_init ', ' you should call ions_base_init first ', 1 )
    !
    tempw     = tempw_
    !
    IF( ALLOCATED( atm2nhp ) ) DEALLOCATE( atm2nhp )
    ALLOCATE( atm2nhp( nat ) )
    !
    atm2nhp(1:nat) = 1
    !
    if (tnosep) then
       nhpcl = MAX( nhpcl_ , 1 )
       if (abs(nhptyp_).eq.1) then
          nhptyp = 1
          if (nhptyp_.gt.0) nhpend = 1
          nhpdim = nsp
          iat = 0
          do is=1,nsp
             do ia=1,na(is)
                iat = iat+1
                atm2nhp(iat) = is
             enddo
          enddo
       elseif (abs(nhptyp_).eq.2) then
          nhptyp = 2
          if (nhptyp_.gt.0) nhpend = 1
          nhpdim = nat
          do i=1,nat
             atm2nhp(i) = i
          enddo
       elseif (abs(nhptyp_).eq.3) then
          nhptyp = 3
          if (nhptyp_.gt.0) nhpend = 1
          call set_atmnhp(nhgrp_,atm2nhp,nhpdim,nhpbeg)
       endif
       ! Add one more chain on top if needed
       nhpdim = nhpdim + nhpend

    endif
    !
    CALL ions_nose_allocate()
    !
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

      ! count the number of atoms per thermostat and set the value
      anum2nhp = 0
      iat = 0
      ! Here we shall check if the scaling factors are provided
      If (maxval(fnhscl_(1:nsp)).lt.0.0d0) then
         scal2nhp = DBLE(ndega)/DBLE(3*nat)
      else
         scal2nhp = -1.0_DP
      endif
      !
      do is=1,nsp
         do ia=1,na(is)
            iat = iat+1
            anum2nhp(atm2nhp(iat)) = anum2nhp(atm2nhp(iat)) + 3
            if (scal2nhp(atm2nhp(iat)).lt.0.0_DP) &
                 scal2nhp(atm2nhp(iat)) = fnhscl_(is)
         enddo
      enddo
      if (nhpend.eq.1) anum2nhp(nhpdim) = nhpdim - 1 - nhpbeg
      ! set gkbt2nhp for each thermostat
      do is=1,nhpdim
         gkbt2nhp(is) = DBLE(anum2nhp(is)) * tempw * k_boltzmann_au
      enddo
      ! scale the target energy by some factor convering 3*nat to ndega
      if (nhpdim.gt.1) then
         do is=1,(nhpdim-nhpend)
            if (scal2nhp(is).lt.0.0_DP) scal2nhp(is) = 1.0_DP
            gkbt2nhp(is) = gkbt2nhp(is)*scal2nhp(is)
         enddo
      endif
      !
      gkbt = DBLE( ndega ) * tempw * k_boltzmann_au
      if (nhpdim.eq.1) gkbt2nhp(1) = gkbt
      kbt  = tempw * k_boltzmann_au

      fnosep(1) = fnosep_ (1)
      if( fnosep(1) > 0.0_DP ) then
        qnp_(1) = 2.0_DP * gkbt / ( fnosep(1) * ( 2.0_DP * pi ) * au_terahertz )**2
      end if

      if ( nhpcl > 1 ) then
        do i = 2, nhpcl
          fnosep(i) = fnosep_ (i)
          if( fnosep(i) > 0.0_DP ) then
            qnp_(i) = 2.0_DP * tempw  * k_boltzmann_au / &
                 ( fnosep(i) * ( 2.0_DP * pi ) * au_terahertz )**2
          else
            qnp_(i) = qnp_(1) / DBLE(ndega)
          endif
        enddo
      endif
      !
      ! set the NH masses for all the chains
      !
      IF((nhptyp.EQ.1).OR.(nhptyp.EQ.2)) THEN 
        !===============================================================
        ! for nhptyp = 1 or 2, the default NH masses are multiplied by 
        ! the atomic masses so that the relative rates of thermostat 
        ! fluctuations are inversely proportional to the masses of the
        ! of the particles to which they are coupled (JOURNAL OF PHYSICAL
        ! CHEMISTRY 97 (49): 12959-12966 DEC 9 1993). Helps in faster
        ! equipartitioning of thermal energy. When nhpend = 1, the mass
        ! of the common thermostat is multiplied by mean atomic mass of 
        ! the system, amass _mean, = total atomic mass / number of atoms  
        !
        DO j=1,(nhpdim-nhpend)
          ! 
          IF(nhptyp.EQ.2) qnp((j-1)*nhpcl+1)=amass(ityp(j))*qnp_(1)*gkbt2nhp(j)/gkbt  
          IF(nhptyp.EQ.1) qnp((j-1)*nhpcl+1)=amass(j)*qnp_(1)*gkbt2nhp(j)/gkbt  
          !
          IF(nhpcl.GT.1) THEN
            !
            DO i=2,nhpcl
              !
              IF(nhptyp.EQ.2) qnp((j-1)*nhpcl+i)=amass(ityp(j))*qnp_(i)
              IF(nhptyp.EQ.1) qnp((j-1)*nhpcl+i)=amass(j)*qnp_(i)
              !
            END DO
            !
          END IF
          ! 
        END DO !j
        !
        IF(nhpend.EQ.1) THEN
          !
          amass_mean=0.0_DP
          !
          DO is=1,nat
            ! 
            amass_mean=amass_mean+amass(ityp(is))
            !
          END DO
          !
          amass_mean=amass_mean/DBLE(nat)
          !
          qnp((nhpdim-1)*nhpcl+1)=amass_mean*qnp_(1)*gkbt2nhp(nhpdim)/gkbt
          !
          IF(nhpcl.GT.1) THEN
            !
            DO i=2,nhpcl
              !
              qnp((nhpdim-1)*nhpcl+i)=amass_mean*qnp_(i)
              !
            END DO
            !
          END IF
          !
        END IF
        !
        !===============================================================
        !
      ELSE
        !
        !===============================================================
        ! Default code : for nhptyp = 1 or 3 
        DO j=1,nhpdim
           qnp((j-1)*nhpcl+1) = qnp_(1)*gkbt2nhp(j)/gkbt
           IF (nhpcl > 1) THEN
              DO i=2,nhpcl
                 qnp((j-1)*nhpcl+i) = qnp_(i)
              END DO
           END IF
        END DO
        ! Default code end 
        !
      END IF !nhptyp
      !
      !
   END IF !tnosep


    !    WRITE( stdout,100)
    !    WRITE( stdout,110) QNOSEP,TEMPW
    !    WRITE( stdout,120) GLIB
    !    WRITE( stdout,130) NSVAR
! 100  FORMAT(//' * Temperature control of ions with nose thermostat'/)
! 110  FORMAT(3X,'nose mass:',F12.4,' temperature (K):',F12.4)
! 120  FORMAT(3X,'ionic degrees of freedom:        ',F5.0)
! 130  FORMAT(3X,'time steps per nose oscillation: ',I5,//)

    return
  end subroutine ions_nose_init

  subroutine set_atmnhp(nhgrp,atm2nhp,nhpdim,nhpbeg)
    !
    use ions_base,      only: nsp, nat, na
    IMPLICIT NONE
    integer, intent(in) :: nhgrp(:)
    integer, intent(out) :: nhpdim, nhpbeg, atm2nhp(:)
    !
    integer :: i,iat,is,ia,igrpmax,ith
    INTEGER, ALLOCATABLE   :: indtmp(:)
    !
    ! find maximum group
    igrpmax = max(maxval(nhgrp(1:nsp)),1)
    ! find out which groups are assigned (assuming gaps)
    allocate(indtmp(igrpmax))
    indtmp=0
    do is=1,nsp
       if (nhgrp(is).gt.0) indtmp(nhgrp(is)) = 1 
    enddo
    ! assign thermostat index to requested groups
    ith = 0
    ! make the 1st thermostat idle if there are negative groups
    if (minval(nhgrp(1:nsp)).lt.0) ith = 1
    nhpbeg = ith
    !
    do i=1,igrpmax
       if (indtmp(i).gt.0) then
          ith = ith + 1
          indtmp(i) = ith
       endif
    enddo
    ! assign thermostats to atoms depending on what is requested
    iat = 0
    do is=1,nsp
       do ia=1,na(is)
          iat = iat+1
          if (nhgrp(is).gt.0) then
             atm2nhp(iat) = indtmp(nhgrp(is))
          elseif (nhgrp(is).eq.0) then
             ith = ith + 1
             atm2nhp(iat) = ith
          else
             atm2nhp(iat) = 1
          endif
       enddo
    enddo
    nhpdim = ith
    deallocate(indtmp)
    return
    !
  end subroutine set_atmnhp

  SUBROUTINE ions_nose_allocate()
    !
    IMPLICIT NONE
    !
    IF ( .NOT. ALLOCATED( vnhp ) )     ALLOCATE( vnhp(  nhpcl*nhpdim ) )
    IF ( .NOT. ALLOCATED( xnhp0 ) )    ALLOCATE( xnhp0( nhpcl*nhpdim ) )
    IF ( .NOT. ALLOCATED( xnhpm ) )    ALLOCATE( xnhpm( nhpcl*nhpdim ) )
    IF ( .NOT. ALLOCATED( xnhpp ) )    ALLOCATE( xnhpp( nhpcl*nhpdim ) )
    IF ( .NOT. ALLOCATED( ekin2nhp ) ) ALLOCATE( ekin2nhp( nhpdim ) )
    IF ( .NOT. ALLOCATED( gkbt2nhp ) ) ALLOCATE( gkbt2nhp( nhpdim ) )
    IF ( .NOT. ALLOCATED( scal2nhp ) ) ALLOCATE( scal2nhp( nhpdim ) )
    IF ( .NOT. ALLOCATED( anum2nhp ) ) ALLOCATE( anum2nhp( nhpdim ) )
    IF ( .NOT. ALLOCATED( qnp ) )      ALLOCATE( qnp( nhpcl*nhpdim ) )
    IF ( .NOT. ALLOCATED( qnp_ ) )     ALLOCATE( qnp_( nhpcl ) )
    IF ( .NOT. ALLOCATED( fnosep ) )   ALLOCATE( fnosep( nhpcl ) )
    !
    vnhp  = 0.0_DP
    xnhp0 = 0.0_DP
    xnhpm = 0.0_DP
    xnhpp = 0.0_DP
    qnp   = 0.0_DP
    qnp_  = 0.0_DP
    !
    RETURN
    !
  END SUBROUTINE ions_nose_allocate

  SUBROUTINE ions_nose_deallocate()
    !
    IMPLICIT NONE
    !
    IF ( ALLOCATED( vnhp ) )     DEALLOCATE( vnhp )
    IF ( ALLOCATED( xnhp0 ) )    DEALLOCATE( xnhp0 )
    IF ( ALLOCATED( xnhpm ) )    DEALLOCATE( xnhpm )
    IF ( ALLOCATED( xnhpp ) )    DEALLOCATE( xnhpp )
    IF ( ALLOCATED( ekin2nhp ) ) DEALLOCATE( ekin2nhp )
    IF ( ALLOCATED( gkbt2nhp ) ) DEALLOCATE( gkbt2nhp )
    IF ( ALLOCATED( scal2nhp ) ) DEALLOCATE( scal2nhp )
    IF ( ALLOCATED( anum2nhp ) ) DEALLOCATE( anum2nhp )
    IF ( ALLOCATED( qnp ) )      DEALLOCATE( qnp )
    IF ( ALLOCATED( qnp_ ) )     DEALLOCATE( qnp_ )
    IF ( ALLOCATED( fnosep ) )   DEALLOCATE( fnosep )
    !
    IF( ALLOCATED( atm2nhp ) ) DEALLOCATE( atm2nhp )
    !
    RETURN
    !
  END SUBROUTINE ions_nose_deallocate

  SUBROUTINE ions_nose_info(delt)

      use constants,     only: au_terahertz, pi
      USE io_global,     ONLY: stdout
      USE control_flags, ONLY: tnosep
      use ions_base,  only: nat

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: delt

      INTEGER   :: nsvar, i, j
      REAL(DP) :: wnosep

      IF( tnosep ) THEN
        !
        IF( fnosep(1) <= 0.0_DP) &
          CALL errore(' ions_nose_info ', ' fnosep less than zero ', 1)
        IF( delt <= 0.0_DP) &
          CALL errore(' ions_nose_info ', ' delt less than zero ', 1)

        wnosep = fnosep(1) * ( 2.0_DP * pi ) * au_terahertz
        nsvar  = ( 2.0_DP * pi ) / ( wnosep * delt )

        WRITE( stdout,563) tempw, nhpcl, ndega, nsvar
        WRITE( stdout,564) (fnosep(i),i=1,nhpcl)
        WRITE( stdout,565) nhptyp, (nhpdim-nhpend), nhpend , nhpbeg, &
             (anum2nhp(j),j=1,nhpdim)
        IF(nhptyp.EQ.1.OR.nhptyp.EQ.2)WRITE( stdout,'(//, &
        & "*** default NH masses are multiplied by atomic masses ***")')
        do j=1,nhpdim
           WRITE( stdout,566) j,(qnp((j-1)*nhpcl+i),i=1,nhpcl)
        enddo
        WRITE( stdout,567)
        do j=1,nat,20
           WRITE( stdout,568) atm2nhp(j:min(nat,j+19))
        enddo
     END IF

 563  format( //, &
     &       3X,'ion dynamics with nose` temperature control:', /, &
     &       3X,'temperature required      = ', f10.5, ' (kelvin) ', /, &
     &       3X,'NH chain length           = ', i3, /, &
     &       3X,'active degrees of freedom = ', i6, /, &
     &       3X,'time steps per nose osc.  = ', i6 )
 564  format( //, &
     &       3X,'nose` frequency(es)       = ', 20(1X,f10.3) ) 
! 565  format( //, &
!     &       3X,'nose` mass(es)            = ', 20(1X,f10.3), // ) 
 565  FORMAT( //, &
     &       3X,'the requested type of NH chains is ',I5, /, &
     &       3X,'total number of thermostats used ',I5,1X,I1,1X,I1, /, &
     &       3X,'ionic degrees of freedom for each chain ',20(1X,I3)) 
 566  format( //, &
     &       3X,'nose` mass(es) for chain ',i4,' = ', 20(1X,f10.3)) 
567  format( //, &
     &       3X,'atom i (in sorted order) is assigned to this thermostat :')
568  format(20(1X,I3))
    RETURN
  END SUBROUTINE ions_nose_info

  

  subroutine ions_nosevel( vnhp, xnhp0, xnhpm, delt, nhpcl, nhpdim )
    implicit none
    integer, intent(in) :: nhpcl, nhpdim
    real(DP), intent(inout) :: vnhp(nhpcl,nhpdim)
    real(DP), intent(in) :: xnhp0(nhpcl,nhpdim), xnhpm(nhpcl,nhpdim), delt
    integer :: i,j
    do j=1,nhpdim
       do i=1,nhpcl
          vnhp(i,j)=2.0_DP * (xnhp0(i,j)-xnhpm(i,j)) / delt-vnhp(i,j)
       end do
    end do
        !
        !  this is equivalent to:
        !  velocity = ( 3.0_DP * xnos0(1) - 4.0_DP * xnosm(1) + xnos2m(1) ) / ( 2.0_DP * delt )
        !  but we do not need variables at time t-2dt ( xnos2m )
        !
    return
  end subroutine ions_nosevel


 subroutine ions_noseupd( xnhpp, xnhp0, xnhpm, delt, qnp, ekin2nhp, gkbt2nhp, vnhp, kbt, nhpcl, nhpdim, nhpbeg, nhpend )
    implicit none
    integer, intent(in) :: nhpcl, nhpdim, nhpbeg, nhpend
    real(DP), intent(out) :: xnhpp(nhpcl,nhpdim)
    real(DP), intent(in) :: xnhp0(nhpcl,nhpdim), xnhpm(nhpcl,nhpdim), delt, qnp(nhpcl,nhpdim), gkbt2nhp(:), kbt
    real(DP), intent(inout) :: ekin2nhp(:), vnhp(nhpcl,nhpdim)
    integer :: i, j
    real(DP) :: dt2, zetfrc, vp1dlt, ekinend, vp1dend


    ekinend = 0.0_DP
    vp1dend = 0.0_DP
    if ( nhpend == 1 ) vp1dend = 0.5_DP * delt * vnhp(1,nhpdim)
    dt2 = delt**2
    if (nhpbeg.gt.0) then
       xnhpp(:,1:nhpbeg) = 0.0_DP
       vnhp (:,1:nhpbeg) = 0.0_DP
    endif
    do j=(1+nhpbeg),nhpdim
    zetfrc = dt2 * ( 2.0_DP * ekin2nhp(j) - gkbt2nhp(j) )
    if ( nhpcl > 1 ) then
       do i=1,(nhpcl-1)
          vp1dlt = 0.5_DP * delt * vnhp(i+1,j)
          xnhpp(i,j)=(2.0_DP * xnhp0(i,j)-(1.0_DP-vp1dlt)*xnhpm(i,j)+zetfrc/qnp(i,j))&
               &   /(1.0_DP+vp1dlt)
!           xnhpp(i,j)=(4.d0*xnhp0(i,j)-(2.d0-delt*vnhp(i+1,j))*xnhpm(i,j)+2.0d0*dt2*zetfrc/qnp(i,j))&
!                &   /(2.d0+delt*vnhp(i+1,j))
          vnhp(i,j) =(xnhpp(i,j)-xnhpm(i,j))/( 2.0_DP * delt )
          zetfrc = dt2*(qnp(i,j)*vnhp(i,j)**2-kbt)
       end do
    end if
    ! Last variable
    i = nhpcl
    if ( nhpend == 0 ) then
       xnhpp(i,j)=2.0_DP * xnhp0(i,j)-xnhpm(i,j) + zetfrc / qnp(i,j)
       vnhp(i,j) =(xnhpp(i,j)-xnhpm(i,j))/( 2.0_DP * delt )
    elseif (nhpend == 1) then
       xnhpp(i,j)=(2.0_DP*xnhp0(i,j)-(1.0_DP-vp1dend)*xnhpm(i,j)+zetfrc/qnp(i,j))&
            &   /(1.0_DP+vp1dend)       
       vnhp(i,j) =(xnhpp(i,j)-xnhpm(i,j))/( 2.0_DP * delt )
       ekinend = ekinend + (qnp(i,j)*vnhp(i,j)**2)
       if (j.eq.(nhpdim-nhpend)) then
          ekin2nhp(nhpdim) = 0.5_DP*ekinend       
          vp1dend = 0.0_DP
       endif
    endif
    enddo
    ! Update velocities
!    do i=1,nhpcl
!       vnhp(i) =(xnhpp(i)-xnhpm(i))/( 2.0d0 * delt )
!    end do
    ! These are the original expressions from cpr.f90
    !       xnhpp(1)=2.*xnhp0(1)-xnhpm(1)+2.*( delt**2 / qnp(1) )*(ekinpr-gkbt/2.)
    !       vnhp(1) =(xnhpp(1)-xnhpm(1))/( 2.0d0 * delt )
    return
  end subroutine ions_noseupd

  

  real(DP) function ions_nose_nrg( xnhp0, vnhp, qnp, gkbt2nhp, kbt, nhpcl, nhpdim )
    implicit none
    integer :: nhpcl, nhpdim
    real(DP) :: gkbt2nhp(:), qnp(nhpcl,nhpdim),vnhp(nhpcl,nhpdim),xnhp0(nhpcl,nhpdim),kbt
    integer :: i,j
    real(DP) :: stmp
    !
    stmp = 0.0_DP
    do j=1,nhpdim
    stmp = stmp + 0.5_DP * qnp(1,j) * vnhp(1,j) * vnhp(1,j) + gkbt2nhp(j) * xnhp0(1,j)
    if (nhpcl > 1) then
       do i=2,nhpcl
          stmp = stmp + 0.5_DP * qnp(i,j) * vnhp(i,j) * vnhp(i,j) + kbt * xnhp0(i,j)
       enddo
    endif
    enddo
    ions_nose_nrg = stmp
    return
  end function ions_nose_nrg



  subroutine ions_nose_shiftvar( xnhpp, xnhp0, xnhpm )
    !  shift values of nose variables to start a new step
    implicit none
    real(DP), intent(inout) :: xnhpm(:), xnhp0(:)
    real(DP), intent(in)  :: xnhpp(:)
      !
      xnhpm = xnhp0
      xnhp0 = xnhpp
      !
    return
  end subroutine ions_nose_shiftvar
  
!------------------------------------------------------------------------------!
  END MODULE ions_nose
!------------------------------------------------------------------------------!
