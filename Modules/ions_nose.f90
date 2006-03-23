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
      USE parameters, ONLY: nhclm
!
      IMPLICIT NONE
! Some comments are in order on how Nose-Hoover chains work here (K.N. Kudin)
! the present code allows one to use "massive" Nose-Hoover chains:
! TOBIAS DJ, MARTYNA GJ, KLEIN ML
! JOURNAL OF PHYSICAL CHEMISTRY 97 (49): 12959-12966 DEC 9 1993
!
! currently "easy" input options allow one chain per atomic type (nhptyp=1)
! and one chain per atom (nhptyp=2), but other options could be added too
! one chain for the whole system is specified by nhptyp=0 (or nothing)
!
! nhpdim is the total number of the resulting NH chains
! nhpend is 1 if there is a chain above all chains, otherwise it is 0
! array atm2nhp(1:nat) gives the chain number from the atom list (which
! is sorted by type)
! anum2nhp is the number of degrees of freedom per chain (now just 3*nat_i)
! ekin2nhp is the kinetic energy of the present chain
! gkbt2nhp are the NH chain parameters
! qnp are the chain masses, qnp_ is a temporary array for now
! see subroutine ions_nose_allocate on what are the dimensions of these
! variables
! nhclm is now mostly not used, needs to be cleaned up at some point
!
      INTEGER   :: nhpcl=1, ndega, nhpdim=1, nhptyp=0, nhpend=0
      INTEGER, ALLOCATABLE   :: atm2nhp(:)
      INTEGER, ALLOCATABLE   :: anum2nhp(:)
      REAL(DP), ALLOCATABLE :: vnhp(:), xnhp0(:), xnhpm(:), xnhpp(:), &
      ekin2nhp(:), gkbt2nhp(:), qnp(:), qnp_(:)

      REAL(DP) :: gkbt = 0.0d0
      REAL(DP) :: kbt = 0.0d0
      REAL(DP) :: tempw = 0.0d0
      REAL(DP) :: fnosep( nhclm ) = 0.0d0

!------------------------------------------------------------------------------!
  CONTAINS 
!------------------------------------------------------------------------------!

  subroutine ions_nose_init( tempw_ , fnosep_ , nhpcl_ , nhptyp_ , ndega_  )
    use constants,      only: factem, pi, terahertz
    use control_flags,  only: tnosep
    use ions_base,      only: ndofp, tions_base_init, nsp, nat, na
    real(DP), intent(in)  :: tempw_ , fnosep_(:) 
    integer, intent(in) :: nhpcl_ , nhptyp_ , ndega_
    integer :: nsvar, i, j, iat, is, ia

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
       endif
       ! Add one more chain on top if needed
       nhpdim = nhpdim + nhpend

       IF( nhpcl > nhclm ) &
            CALL errore(' ions_nose_init ', ' nhpcl out of range ', nhpcl )
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
      do is=1,nsp
         do ia=1,na(is)
            iat = iat+1
            anum2nhp(atm2nhp(iat)) = anum2nhp(atm2nhp(iat)) + 3
         enddo
      enddo
      if (nhpend.eq.1) anum2nhp(nhpdim) = nhpdim - 1
      ! set gkbt2nhp for each thermostat
      do is=1,nhpdim
         gkbt2nhp(is) = DBLE(anum2nhp(is)) * tempw / factem
      enddo
      ! scale the target energy by some factor convering 3*nat to ndega
      if (nhpdim.gt.1) then
         do is=1,(nhpdim-nhpend)
            gkbt2nhp(is) = gkbt2nhp(is)*DBLE(ndega)/DBLE(3*nat)
         enddo
      endif
      !
      gkbt = DBLE( ndega ) * tempw / factem
      if (nhpdim.eq.1) gkbt2nhp(1) = gkbt
      kbt  = tempw / factem

      fnosep(1) = fnosep_ (1)
      if( fnosep(1) > 0.0d0 ) then
        qnp_(1) = 2.d0 * gkbt / ( fnosep(1) * ( 2.d0 * pi ) * terahertz )**2
      end if

      if ( nhpcl > 1 ) then
        do i = 2, nhpcl
          fnosep(i) = fnosep_ (i)
          if( fnosep(i) > 0.0d0 ) then
            qnp_(i) = 2.d0 * tempw / factem / ( fnosep(i) * ( 2.d0 * pi ) * terahertz )**2
          else
            qnp_(i) = qnp_(1) / DBLE(ndega)
          endif
        enddo
      endif
      ! set the NH masses for all the chains
      do j=1,nhpdim
         qnp((j-1)*nhpcl+1) = qnp_(1)*gkbt2nhp(j)/gkbt
         If (nhpcl > 1) then
            do i=2,nhpcl
               qnp((j-1)*nhpcl+i) = qnp_(i)
            enddo
         endif
      enddo
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
  end subroutine ions_nose_init

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
    IF ( .NOT. ALLOCATED( anum2nhp ) ) ALLOCATE( anum2nhp( nhpdim ) )
    IF ( .NOT. ALLOCATED( qnp ) )      ALLOCATE( qnp( nhpcl*nhpdim ) )
    IF ( .NOT. ALLOCATED( qnp_ ) )     ALLOCATE( qnp_( nhpcl ) )
    !
    vnhp  = 0.D0
    xnhp0 = 0.D0
    xnhpm = 0.D0
    xnhpp = 0.D0
    qnp   = 0.D0
    qnp_  = 0.D0
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
    IF ( ALLOCATED( anum2nhp ) ) DEALLOCATE( anum2nhp )
    IF ( ALLOCATED( qnp ) )      DEALLOCATE( qnp )
    IF ( ALLOCATED( qnp_ ) )     DEALLOCATE( qnp_ )
    !
    IF( ALLOCATED( atm2nhp ) ) DEALLOCATE( atm2nhp )
    !
    RETURN
    !
  END SUBROUTINE ions_nose_deallocate

  SUBROUTINE ions_nose_info()

      use constants,     only: factem, terahertz, pi
      use time_step,     only: delt
      USE io_global,     ONLY: stdout
      USE control_flags, ONLY: tnosep

      IMPLICIT NONE

      INTEGER   :: nsvar, i, j
      REAL(DP) :: wnosep

      IF( tnosep ) THEN
        !
        IF( fnosep(1) <= 0.D0) &
          CALL errore(' ions_nose_info ', ' fnosep less than zero ', 1)
        IF( delt <= 0.D0) &
          CALL errore(' ions_nose_info ', ' delt less than zero ', 1)

        wnosep = fnosep(1) * ( 2.d0 * pi ) * terahertz
        nsvar  = ( 2.d0 * pi ) / ( wnosep * delt )

        WRITE( stdout,563) tempw, nhpcl, ndega, nsvar
        WRITE( stdout,564) (fnosep(i),i=1,nhpcl)
        WRITE( stdout,565) nhptyp, (nhpdim-nhpend), nhpend , &
             (anum2nhp(j),j=1,nhpdim)
        do j=1,nhpdim
           WRITE( stdout,566) j,(qnp((j-1)*nhpcl+i),i=1,nhpcl)
        enddo
     END IF

 563  format( //, &
     &       3X,'ion dynamics with nose` temperature control:', /, &
     &       3X,'temperature required      = ', f10.5, ' (kelvin) ', /, &
     &       3X,'NH chain length           = ', i3, /, &
     &       3X,'active degrees of freedom = ', i3, /, &
     &       3X,'time steps per nose osc.  = ', i5 )
 564  format( //, &
     &       3X,'nose` frequency(es)       = ', 20(1X,f10.3) ) 
! 565  format( //, &
!     &       3X,'nose` mass(es)            = ', 20(1X,f10.3), // ) 
 565  FORMAT( //, &
     &       3X,'the requested type of NH chains is ',I5, /, &
     &       3X,'total number of thermostats used ',I5,1X,I1, /, &
     &       3X,'ionic degrees of freedom for each chain ',20(1X,I3)) 
 566  format( //, &
     &       3X,'nose` mass(es) for chain ',i4,' = ', 20(1X,f10.3)) 
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
          vnhp(i,j)=2.*(xnhp0(i,j)-xnhpm(i,j))/delt-vnhp(i,j)
       enddo
    enddo
        !
        !  this is equivalent to:
        !  velocity = ( 3.D0 * xnos0(1) - 4.D0 * xnosm(1) + xnos2m(1) ) / ( 2.0d0 * delt )
        !  but we do not need variables at time t-2dt ( xnos2m )
        !
    return
  end subroutine ions_nosevel


 subroutine ions_noseupd( xnhpp, xnhp0, xnhpm, delt, qnp, ekin2nhp, gkbt2nhp, vnhp, kbt, nhpcl, nhpdim, nhpend )
    implicit none
    integer, intent(in) :: nhpcl, nhpdim, nhpend
    real(DP), intent(out) :: xnhpp(nhpcl,nhpdim)
    real(DP), intent(in) :: xnhp0(nhpcl,nhpdim), xnhpm(nhpcl,nhpdim), delt, qnp(nhpcl,nhpdim), gkbt2nhp(:), kbt
    real(DP), intent(inout) :: ekin2nhp(:), vnhp(nhpcl,nhpdim)
    integer :: i, j
    real(DP) :: dt2, zetfrc, vp1dlt, ekinend, vp1dend


    ekinend = 0.0d0
    vp1dend = 0.0d0
    if (nhpend.eq.1) vp1dend = 0.5d0*delt*vnhp(1,nhpdim)
    dt2 = delt**2
    do j=1,nhpdim
    zetfrc = dt2*(2.0d0*ekin2nhp(j)-gkbt2nhp(j))
    If (nhpcl.gt.1) then
       do i=1,(nhpcl-1)
          vp1dlt = 0.5d0*delt*vnhp(i+1,j)
          xnhpp(i,j)=(2.d0*xnhp0(i,j)-(1.d0-vp1dlt)*xnhpm(i,j)+zetfrc/qnp(i,j))&
               &   /(1.d0+vp1dlt)
!           xnhpp(i,j)=(4.d0*xnhp0(i,j)-(2.d0-delt*vnhp(i+1,j))*xnhpm(i,j)+2.0d0*dt2*zetfrc/qnp(i,j))&
!                &   /(2.d0+delt*vnhp(i+1,j))
          vnhp(i,j) =(xnhpp(i,j)-xnhpm(i,j))/( 2.0d0 * delt )
          zetfrc = dt2*(qnp(i,j)*vnhp(i,j)**2-kbt)
       enddo
    endif
    ! Last variable
    i = nhpcl
    if (nhpend.eq.0) then
       xnhpp(i,j)=2.d0*xnhp0(i,j)-xnhpm(i,j)+ zetfrc / qnp(i,j)
       vnhp(i,j) =(xnhpp(i,j)-xnhpm(i,j))/( 2.0d0 * delt )
    elseif (nhpend.eq.1) then
       xnhpp(i,j)=(2.d0*xnhp0(i,j)-(1.d0-vp1dend)*xnhpm(i,j)+zetfrc/qnp(i,j))&
            &   /(1.d0+vp1dend)       
       vnhp(i,j) =(xnhpp(i,j)-xnhpm(i,j))/( 2.0d0 * delt )
       ekinend = ekinend + (qnp(i,j)*vnhp(i,j)**2)
       if (j.eq.(nhpdim-nhpend)) then
          ekin2nhp(nhpdim) = 0.5d0*ekinend       
          vp1dend = 0.0d0
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
    stmp = 0.0d0
    do j=1,nhpdim
    stmp = stmp + 0.5d0 * qnp(1,j) * vnhp(1,j) * vnhp(1,j) + gkbt2nhp(j) * xnhp0(1,j)
    if (nhpcl > 1) then
       do i=2,nhpcl
          stmp=stmp+0.5*qnp(i,j)*vnhp(i,j)*vnhp(i,j) + kbt*xnhp0(i,j)
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
