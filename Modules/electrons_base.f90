!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------!
  MODULE electrons_base
!------------------------------------------------------------------------------!

      USE kinds, ONLY: DP
!
      IMPLICIT NONE
      SAVE

      INTEGER :: nbnd       = 0    !  number electronic bands, each band contains
                                   !  two spin states
      INTEGER :: nbndx      = 0    !  array dimension nbndx >= nbnd
      INTEGER :: nspin      = 0    !  nspin = number of spins (1=no spin, 2=LSDA)
      INTEGER :: nel(2)     = 0    !  number of electrons (up, down)
      INTEGER :: nelt       = 0    !  total number of electrons ( up + down )
      INTEGER :: nupdwn(2)  = 0    !  number of states with spin up (1) and down (2)
      INTEGER :: iupdwn(2)  = 0    !  first state with spin (1) and down (2)
      INTEGER :: nudx       = 0    !  max (nupdw(1),nupdw(2))
      INTEGER :: nbsp       = 0    !  total number of electronic states 
                                   !  (nupdwn(1)+nupdwn(2))
      INTEGER :: nbspx      = 0    !  array dimension nbspx >= nbsp

      LOGICAL :: telectrons_base_initval = .FALSE.

      REAL(DP), ALLOCATABLE :: f(:)   ! occupation numbers ( at gamma )
      REAL(DP) :: qbac = 0.0d0        ! background neutralizing charge
      INTEGER, ALLOCATABLE :: fspin(:) ! spin of each state
!
!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!


    SUBROUTINE electrons_base_initval( zv_ , na_ , nsp_ , nelec_ , nelup_ , neldw_ , nbnd_ , &
               nspin_ , occupations_ , f_inp, tot_charge_, multiplicity_ )

      USE constants, ONLY: eps8
      USE io_global, ONLY: stdout
      USE control_flags, ONLY: iprsta

      REAL(DP), INTENT(IN) :: zv_ (:), tot_charge_
      INTEGER, INTENT(IN) :: na_ (:) , nsp_, multiplicity_
      REAL(DP), INTENT(IN) :: nelec_ , nelup_ , neldw_
      INTEGER, INTENT(IN) :: nbnd_ , nspin_
      CHARACTER(LEN=*), INTENT(IN) :: occupations_
      REAL(DP), INTENT(IN) :: f_inp(:,:)
      REAL(DP) :: nelec, nelup, neldw, ocp, fsum
      INTEGER   :: iss, i, in

      nspin = nspin_

      IF( nelec_ /= 0 ) THEN
        nelec = nelec_
      ELSE
        nelec = 0.0d0
        DO i = 1, nsp_
          nelec = nelec + na_ ( i ) * zv_ ( i )
        END DO 
        nelec = nelec - tot_charge_
      END IF

      IF( nelup_ > 0.0d0 .AND. neldw_ > 0.0d0 ) THEN
        nelup = nelup_
        neldw = neldw_
      ELSE IF( nelup_ > 0.0d0 .AND. neldw_ == 0.0d0 ) THEN
        nelup = nelup_
        neldw = nelec - nelup_
      ELSE IF( nelup_ == 0.0d0 .AND. neldw_ > 0.0d0 ) THEN
        neldw = neldw_
        nelup = nelec - neldw_
      ELSE
         IF ( multiplicity_ == 0 ) THEN
            ! default when multiplicity is unspecified
            nelup = INT( nelec + 1 ) / 2
            neldw = nelec - nelup
         else
            ! non singlet multiplicity requires nspin=2
            if ( (multiplicity_ > 1) .and. (nspin==1) ) &
                 CALL errore(' electrons_base_initval ', &
                 ' spin multiplicity incosistent with nspin=1 ', 2 )
            ! odd multiplicity requires an even number of electrons
            ! even multiplicity requires an odd number of electrons
            if ( ((MOD(multiplicity_,2) == 0) .and. (MOD(NINT(nelec),2)==0)) .or.   &
                 ((MOD(multiplicity_,2) == 1) .and. (MOD(NINT(nelec),2)==1))      ) &
                 CALL errore(' electrons_base_initval ',                           &
                 ' spin multiplicity incosistent with total number of electrons ', 2 )

            nelup = ( INT(nelec) + ( multiplicity_-1 ) ) / 2
            neldw = ( INT(nelec) - ( multiplicity_-1 ) ) / 2
         end IF
      END IF

      IF( ABS( nelec - ( nelup + neldw ) ) > eps8 ) THEN
         CALL errore(' electrons_base_initval ',' inconsistent n. of electrons ', 2 )
      END IF
      !
      !   Compute the number of bands
      !
      IF( nbnd_ /= 0 ) THEN
        nbnd  = nbnd_                          ! nbnd is given from input
      ELSE
        nbnd  = NINT( MAX( nelup, neldw ) )    ! take the maximum between up and down states
      END IF


      IF( nelec < 1 ) THEN
         CALL errore(' electrons_base_initval ',' nelec less than 1 ', 1 )
      END IF
      !
      IF( ABS( NINT( nelec ) - nelec ) > eps8 ) THEN
         CALL errore(' electrons_base_initval ',' nelec must be integer', 2 )
      END IF
      !
      IF( nbnd < 1 ) &
        CALL errore(' electrons_base_initval ',' nbnd out of range ', 1 )
      !

      IF ( nspin /= 1 .AND. nspin /= 2 ) THEN
        WRITE( stdout, * ) 'nspin = ', nspin
        CALL errore( ' electrons_base_initval ', ' nspin out of range ', 1 )
      END IF

      IF( MOD( nbnd, 2 ) == 0 ) THEN
         nbspx = nbnd * nspin
      ELSE
         nbspx = ( nbnd + 1 ) * nspin
      END IF

      ALLOCATE( f     ( nbspx ) )
      ALLOCATE( fspin ( nbspx ) )
      f     = 0.0d0
      fspin = 0

      iupdwn ( 1 ) = 1
      nel          = 0

      SELECT CASE ( TRIM(occupations_) )
      CASE ('bogus')
         !
         ! empty-states calculation: occupancies have a (bogus) finite value
         !
         ! bogus to ensure \sum_i f_i = Nelec  (nelec is integer)
         !
         f ( : )    = nelec / nbspx
         nel (1)    = nint( nelec )
         nupdwn (1) = nbspx
         if ( nspin == 2 ) then
            !
            ! bogus to ensure Nelec = Nup + Ndw
            !
            nel (1) = ( nint(nelec) + 1 ) / 2
            nel (2) =   nint(nelec)       / 2
            nupdwn (1)=nbnd
            nupdwn (2)=nbnd
            iupdwn (2)=nbnd+1
         end if
         !
      CASE ('from_input')
         !
         ! occupancies have been read from input
         !
         ! count electrons
         !
         IF( nspin == 1 ) THEN
            nelec = SUM( f_inp( :, 1 ) )
            nelup = nelec / 2.0d0
            neldw = nelec / 2.0d0
         ELSE
            nelup = SUM ( f_inp ( :, 1 ) )
            neldw = SUM ( f_inp ( :, 2 ) )
            nelec = nelup + neldw 
         END IF
         !
         ! count bands
         !
         nupdwn (1) = 0
         nupdwn (2) = 0
         DO i = 1, nbnd
           IF( nspin == 1 ) THEN
             IF( f_inp( i, 1 ) > 0.0d0 ) nupdwn (1) = nupdwn (1) + 1
           ELSE
             IF( f_inp( i, 1 ) > 0.0d0 ) nupdwn (1) = nupdwn (1) + 1
             IF( f_inp( i, 2 ) > 0.0d0 ) nupdwn (2) = nupdwn (2) + 1
           END IF
         END DO
         !
         if( nspin == 1 ) then
           nel (1)    = nint( nelec )
           iupdwn (1) = 1
         else
           nel (1)    = nint(nelup)
           nel (2)    = nint(neldw)
           iupdwn (1) = 1
           iupdwn (2) = nupdwn (1) + 1
         end if
         !
         DO iss = 1, nspin
           DO in = iupdwn ( iss ), iupdwn ( iss ) - 1 + nupdwn ( iss )
             f( in ) = f_inp( in - iupdwn ( iss ) + 1, iss )
           END DO
         END DO
         !
      CASE ('fixed')

         if( nspin == 1 ) then
            nel(1)    = nint(nelec)
            nupdwn(1) = nbnd
            iupdwn(1) = 1
         else
            IF ( nelup + neldw /= nelec  ) THEN
               CALL errore(' electrons_base_initval ',' wrong # of up and down spin', 1 )
            END IF
            nel(1)    = nint(nelup)
            nel(2)    = nint(neldw)
            nupdwn(1) = nint(nelup)
            nupdwn(2) = nint(neldw)
            iupdwn(1) = 1
            iupdwn(2) = nupdwn(1) + 1
         end if

         if( (nspin == 1) .and. MOD( nint(nelec), 2 ) /= 0 ) &
              CALL errore(' electrons_base_initval ', &
              ' must use nspin=2 for odd number of electrons', 1 )
         
         ! ocp = 2 for spinless systems, ocp = 1 for spin-polarized systems
         ocp = 2.d0 / nspin
         !
         ! default filling: attribute ocp electrons to each states
         !                  until the good number of electrons is reached
         do iss = 1, nspin
            fsum = 0.0d0
            do in = iupdwn ( iss ), iupdwn ( iss ) - 1 + nupdwn ( iss )
               if ( fsum + ocp < nel ( iss ) + 0.0001 ) then
                  f (in) = ocp
               else
                  f (in) = max( nel ( iss ) - fsum, 0.d0 )
               end if
                fsum = fsum + f(in)
            end do
         end do
         !
      CASE ('ensemble','ensemble-dft','edft')

          if ( nspin == 1 ) then
            !
            f ( : )    = nelec / nbnd
            nel (1)    = nint(nelec)
            nupdwn (1) = nbnd
            !
          else
            !
            if (nelup.ne.0) then
              if ((nelup+neldw).ne.nelec) then
                 CALL errore(' electrons_base_initval ',' nelup+neldw .ne. nelec', 1 )
              end if
              nel (1) = nelup
              nel (2) = neldw
            else
              nel (1) = ( nint(nelec) + 1 ) / 2
              nel (2) =   nint(nelec)       / 2
            end if
            !
            nupdwn (1) = nbnd
            nupdwn (2) = nbnd
            iupdwn (2) = nbnd+1
            !
            do iss = 1, nspin
             do i = iupdwn ( iss ), iupdwn ( iss ) - 1 + nupdwn ( iss )
                f (i) =  nel (iss) / DBLE (nupdwn (iss))
             end do
            end do
            !
          end if

      CASE DEFAULT
         CALL errore(' electrons_base_initval ',' occupation method not implemented', 1 )
      END SELECT


      do iss = 1, nspin
         do in = iupdwn(iss), iupdwn(iss) - 1 + nupdwn(iss)
            fspin(in) = iss
         end do
      end do

      nbndx = MAXVAL( nupdwn )
      nbsp  = nupdwn (1) + nupdwn (2)

      IF ( nspin == 1 ) THEN 
        nelt = nel(1)
        nudx = nupdwn(1)
      ELSE
        nelt = nel(1) + nel(2)
        nudx = MAX( nupdwn(1), nupdwn(2) )
      END IF

      IF( nbnd < MAX( nupdwn(1), nupdwn(2) ) ) &
        CALL errore(' electrons_base_initval ',' inconsistent nbnd and nupdwn(1) or nupdwn(2) ', 1 )

      IF( ( 2 * nbnd ) < nelt ) &
        CALL errore(' electrons_base_initval ',' too few states ',  1  )

      IF( nbsp < INT( nelec * nspin / 2.0d0 ) ) &
        CALL errore(' electrons_base_initval ',' too many electrons ', 1 )


      telectrons_base_initval = .TRUE.

      RETURN

    END SUBROUTINE electrons_base_initval



    SUBROUTINE deallocate_elct()
      IF( ALLOCATED( f ) ) DEALLOCATE( f )
      IF( ALLOCATED( fspin ) ) DEALLOCATE( fspin )
      telectrons_base_initval = .FALSE.
      RETURN
    END SUBROUTINE deallocate_elct


!------------------------------------------------------------------------------!
  END MODULE electrons_base
!------------------------------------------------------------------------------!



!------------------------------------------------------------------------------!
  MODULE electrons_nose
!------------------------------------------------------------------------------!

      USE kinds, ONLY: DP
!
      IMPLICIT NONE
      SAVE

      REAL(DP) :: fnosee   = 0.0d0   !  frequency of the thermostat ( in THz )
      REAL(DP) :: qne      = 0.0d0   !  mass of teh termostat
      REAL(DP) :: ekincw   = 0.0d0   !  kinetic energy to be kept constant

      REAL(DP) :: xnhe0   = 0.0d0   
      REAL(DP) :: xnhep   = 0.0d0   
      REAL(DP) :: xnhem   = 0.0d0   
      REAL(DP) :: vnhe    = 0.0d0
      
!
!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!

  subroutine electrons_nose_init( ekincw_ , fnosee_ )
     USE constants, ONLY: factem, pi, terahertz
     REAL(DP), INTENT(IN) :: ekincw_, fnosee_
     ! set thermostat parameter for electrons
     qne     = 0.0d0
     ekincw  = ekincw_
     fnosee  = fnosee_
     xnhe0   = 0.0d0
     xnhep   = 0.0d0
     xnhem   = 0.0d0
     vnhe    = 0.0d0
     if( fnosee > 0.0d0 ) qne = 4.d0 * ekincw / ( fnosee * ( 2.d0 * pi ) * terahertz )**2
    return
  end subroutine electrons_nose_init


  function electrons_nose_nrg( xnhe0, vnhe, qne, ekincw )
    !  compute energy term for nose thermostat
    implicit none
    real(8) :: electrons_nose_nrg
    real(8), intent(in) :: xnhe0, vnhe, qne, ekincw
      !
      electrons_nose_nrg = 0.5d0 * qne * vnhe * vnhe + 2.0d0 * ekincw * xnhe0
      !
    return
  end function electrons_nose_nrg

  subroutine electrons_nose_shiftvar( xnhep, xnhe0, xnhem )
    !  shift values of nose variables to start a new step
    implicit none
    real(8), intent(out) :: xnhem
    real(8), intent(inout) :: xnhe0
    real(8), intent(in) :: xnhep
      !
      xnhem = xnhe0
      xnhe0 = xnhep
      !
    return
  end subroutine electrons_nose_shiftvar

  subroutine electrons_nosevel( vnhe, xnhe0, xnhem, delt )
    implicit none
    real(8), intent(inout) :: vnhe
    real(8), intent(in) :: xnhe0, xnhem, delt 
    vnhe=2.*(xnhe0-xnhem)/delt-vnhe
    return
  end subroutine electrons_nosevel

  subroutine electrons_noseupd( xnhep, xnhe0, xnhem, delt, qne, ekinc, ekincw, vnhe )
    implicit none
    real(8), intent(out) :: xnhep, vnhe
    real(8), intent(in) :: xnhe0, xnhem, delt, qne, ekinc, ekincw
    xnhep = 2.0d0 * xnhe0 - xnhem + 2.0d0 * ( delt**2 / qne ) * ( ekinc - ekincw )
    vnhe  = ( xnhep - xnhem ) / ( 2.0d0 * delt )
    return
  end subroutine electrons_noseupd


  SUBROUTINE electrons_nose_info()

      use constants,     only: factem, terahertz, pi
      use time_step,     only: delt
      USE io_global,     ONLY: stdout
      USE control_flags, ONLY: tnosee

      IMPLICIT NONE

      INTEGER   :: nsvar, i
      REAL(DP) :: wnosee

      IF( tnosee ) THEN
        !
        IF( fnosee <= 0.D0) &
          CALL errore(' electrons_nose_info ', ' fnosee less than zero ', 1)
        IF( delt <= 0.D0) &
          CALL errore(' electrons_nose_info ', ' delt less than zero ', 1)

        wnosee = fnosee * ( 2.d0 * pi ) * terahertz
        nsvar  = ( 2.d0 * pi ) / ( wnosee * delt )

        WRITE( stdout,563) ekincw, nsvar, fnosee, qne
      END IF

 563  format( //, &
     &       3X,'electrons dynamics with nose` temperature control:', /, &
     &       3X,'Kinetic energy required   = ', f10.5, ' (a.u.) ', /, &
     &       3X,'time steps per nose osc.  = ', i5, /, &
     &       3X,'nose` frequency           = ', f10.3, ' (THz) ', /, &
     &       3X,'nose` mass(es)            = ', 20(1X,f10.3),//)


    RETURN
  END SUBROUTINE electrons_nose_info



!------------------------------------------------------------------------------!
  END MODULE electrons_nose
!------------------------------------------------------------------------------!
