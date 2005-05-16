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

      USE kinds, ONLY: dbl
      USE parameters, ONLY: nhclm
!
      IMPLICIT NONE

      INTEGER   :: nhpcl, ndega
      REAL(dbl) :: vnhp( nhclm ) = 0.0d0
      REAL(dbl) :: xnhp0( nhclm ) = 0.0d0
      REAL(dbl) :: xnhpm( nhclm ) = 0.0d0
      REAL(dbl) :: xnhpp( nhclm ) = 0.0d0
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
    real(KIND=dbl), intent(in)  :: tempw_ , fnosep_(:) 
    integer, intent(in) :: nhpcl_ , ndega_
    integer, intent(in) :: nat
    integer :: nsvar, i

    IF( .NOT. tions_base_init ) &
      CALL errore(' ions_nose_init ', ' you should call ions_base_init first ', 1 )

    vnhp  = 0.0d0
    xnhp0 = 0.0d0
    xnhpm = 0.0d0 
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

        WRITE( stdout,563) tempw, nhpcl, ndega, nsvar
        WRITE( stdout,564) (fnosep(i),i=1,nhpcl)
        WRITE( stdout,565) (qnp(i),i=1,nhpcl)
      END IF

 563  format( //, &
            & 3X,'ion dynamics with nose` temperature control:', /, &
            & 3X,'temperature required      = ', f10.5, ' (kelvin) ', /, &
            & 3X,'NH chain length           = ', i3, /, &
            & 3X,'active degrees of freedom = ', i3, /, &
            & 3X,'time steps per nose osc.  = ', i5 )
 564  format( //, &
            & 3X,'nose` frequency(es)       = ', 20(1X,f10.3) ) 
 565  format( //, &
            & 3X,'nose` mass(es)            = ', 20(1X,f10.3), // ) 


    RETURN
  END SUBROUTINE ions_nose_info

  

  subroutine ions_nosevel( vnhp, xnhp0, xnhpm, delt, nhpcl )
    implicit none
    real(KIND=dbl), intent(inout) :: vnhp(:)
    real(KIND=dbl), intent(in) :: xnhp0(:), xnhpm(:), delt
    integer, intent(in) :: nhpcl
    integer :: i
    do i=1,nhpcl
       vnhp(i)=2.*(xnhp0(i)-xnhpm(i))/delt-vnhp(i)
    enddo
        !
        !  this is equivalent to:
        !  velocity = ( 3.D0 * xnos0(1) - 4.D0 * xnosm(1) + xnos2m(1) ) / ( 2.0d0 * delt )
        !  but we do not need variables at time t-2dt ( xnos2m )
        !
    return
  end subroutine


 subroutine ions_noseupd( xnhpp, xnhp0, xnhpm, delt, qnp, ekinpr, gkbt, vnhp, kbt, nhpcl )
    implicit none
    real(KIND=dbl), intent(out) :: xnhpp(:), vnhp(:)
    real(KIND=dbl), intent(in) :: xnhp0(:), xnhpm(:), delt, qnp(:), ekinpr, gkbt, kbt
    integer, intent(in) :: nhpcl
    integer :: i
    real(KIND=dbl) :: dt2, zetfrc

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

  

  real(KIND=dbl) function ions_nose_nrg( xnhp0, vnhp, qnp, gkbt, kbt, nhpcl )
    implicit none
    integer :: nhpcl
    real(KIND=dbl) :: gkbt,qnp(:),vnhp(:),xnhp0(:),kbt
    integer :: i
    real(KIND=dbl) :: stmp
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



  subroutine ions_nose_shiftvar( xnhpp, xnhp0, xnhpm )
    !  shift values of nose variables to start a new step
    implicit none
    real(KIND=dbl), intent(inout) :: xnhpm(:), xnhp0(:)
    real(KIND=dbl), intent(in)  :: xnhpp(:)
      !
      xnhpm = xnhp0
      xnhp0 = xnhpp
      !
    return
  end subroutine
  
!------------------------------------------------------------------------------!
  END MODULE ions_nose
!------------------------------------------------------------------------------!
