! cp_autopilot.f90
!********************************************************************************
! cp_autopilot.f90                             Copyright (c) 2005 Targacept, Inc.
!********************************************************************************
!   The Autopilot Feature suite is a user level enhancement that enables the 
! following features:  
!      automatic restart of a job; 
!      preconfiguration of job parameters; 
!      on-the-fly changes to job parameters; 
!      and pausing of a running job.  
!
! For more information, see AUTOPILOT in document directory.
!
! This program is free software; you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful, but WITHOUT ANY 
! WARRANTY; without even the implied warranty of MERCHANTABILITY FOR A PARTICULAR 
! PURPOSE.  See the GNU General Public License at www.gnu.or/copyleft/gpl.txt for 
! more details.
! 
! THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  
! EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES 
! PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESS OR IMPLIED, 
! INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND THE 
! PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, 
! YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
!
! IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING, 
! WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR REDISTRIBUTE 
! THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY 
! GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR 
! INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA 
! BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A 
! FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH HOLDER 
! OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
!
! You should have received a copy of the GNU General Public License along with 
! this program; if not, write to the 
! Free Software Foundation, Inc., 
! 51 Franklin Street, 
! Fifth Floor, 
! Boston, MA  02110-1301, USA.
! 
! Targacept's address is 
! 200 East First Street, Suite 300
! Winston-Salem, North Carolina USA 27101-4165 
! Attn: Molecular Design. 
! Email: atp@targacept.com
!
! This work was supported by the Advanced Technology Program of the 
! National Institute of Standards and Technology (NIST), Award No. 70NANB3H3065 
!
!********************************************************************************
 

MODULE cp_autopilot
  !---------------------------------------------------------------------------
  !
  ! This module handles the Autopilot Feature Suite
  ! Written by Lee Atkinson, with help from the ATP team at Targacept, Inc 
  ! Created June 2005
  ! Modified by Yonas Abrahm Sept 2006
  !
  !   The address for Targacept, Inc. is:
  !     200 East First Street, Suite
  !     300, Winston-Salem, North Carolina 27101; 
  !     Attn: Molecular Design.
  !
  ! See README.AUTOPILOT in the Doc directory for more information.
  !---------------------------------------------------------------------------

  USE kinds
  USE parser, ONLY :  read_line
  USE autopilot, ONLY : current_nfi, pilot_p, pilot_unit, pause_p,auto_error, &
        &  parse_mailbox, rule_isave, rule_iprint, rule_dt, rule_emass,       &
        &  rule_electron_dynamics, rule_electron_damping, rule_ion_dynamics, &
        &  rule_ion_damping,  rule_ion_temperature, rule_tempw
  USE autopilot, ONLY : event_index, event_step, event_isave, event_iprint, &
        &  event_dt, event_emass, event_electron_dynamics, event_electron_damping, &
        &  event_ion_dynamics, event_ion_damping, event_ion_temperature, event_tempw 

  IMPLICIT NONE
  SAVE


  PRIVATE
  PUBLIC ::  pilot, employ_rules

CONTAINS


  !-----------------------------------------------------------------------
  ! EMPLOY_RULES
  !-----------------------------------------------------------------------
  SUBROUTINE employ_rules()
    USE input_parameters, ONLY :   dt, & 
         &  electron_dynamics, electron_damping, &
         & ion_dynamics, ion_damping, &
         & ion_temperature, fnosep, nhpcl, nhptyp, nhgrp, fnhscl, ndega, nat
    use ions_nose, ONLY: tempw
    USE control_flags, only: tsde, tsdp, tfor, tcp, tnosep, isave,iprint,&
                             tconvthrs, tolp,ldamped, &
                             ekin_conv_thr, forc_conv_thr, etot_conv_thr
    USE control_flags, only: tsteepdesc_    => tsteepdesc, &
                             tdamp_         => tdamp,      &
                             tdampions_ => tdampions
    use wave_base, only: frice
    use ions_base, only: fricp
    USE ions_nose, ONLY: ions_nose_init
    USE io_global, ONLY: ionode, ionode_id
    USE time_step,          ONLY : set_time_step
    USE cp_electronic_mass, ONLY: emass
    IMPLICIT NONE


    ! This is notification to stdout
    ! It helps the user to identify 
    ! when rules are employed
    write(*,*)
    write(*,*) '========================================'
    write(*,*) 'EMPLOY RULES:'
    write(*,*) '  CURRENT_NFI=', current_nfi
    write(*,*) '  event_index=', event_index
    write(*,*) '  event_step==', event_step(event_index)
    write(*,*) '========================================'
    write(*,*)
    call flush_unit(6)


    !----------------------------------------
    !     &CONTROL
    !----------------------------------------

    ! ISAVE
    if (event_isave(event_index)) then
       isave            = rule_isave(event_index)
       write(*,*) 'RULE EVENT: isave', isave
    endif

    ! IPRINT
    if (event_iprint(event_index)) then
       iprint           = rule_iprint(event_index)
       write(*,*) 'RULE EVENT: iprint', iprint
    endif

    if (event_dt(event_index)) then
       dt               = rule_dt(event_index)
       CALL set_time_step( dt )
       write(*,*) 'RULE EVENT: dt', dt
    endif

    !----------------------------------------
    !     &SYSTEM
    !----------------------------------------

    !----------------------------------------
    !     &ELECTRONS
    !----------------------------------------

    ! EMASS    
    if (event_emass(event_index)) then
       emass            = rule_emass(event_index)
       write(*,*) 'RULE EVENT: emass', emass
    endif

    ! ELECTRON_DYNAMICS   

    if (event_electron_dynamics(event_index)) then
       electron_dynamics= rule_electron_dynamics(event_index)
      tdamp_          = .FALSE.
      tsteepdesc_     = .FALSE.
      frice = 0.d0
       select case ( electron_dynamics ) 
       case ('SD')
          tsde  = .true.
       case ('VERLET')
          tsde  = .false.
       case ('DAMP')
          tsde  = .false.
          tdamp_  = .TRUE.
          frice = electron_damping
       case ('NONE')
          tsde  = .false.
       case default
          call auto_error(' autopilot ',' unknown electron_dynamics '//trim(electron_dynamics) )
       end select

       write(*,*) 'RULE EVENT: electron_dynamics', electron_dynamics

    endif


    ! ELECTRON_DAMPING   
    if (event_electron_damping(event_index)) then
       ! meaningful only if " electron_dynamics = 'damp' "
       electron_damping = rule_electron_damping(event_index)
       frice = electron_damping
       write(*,*) 'RULE EVENT: electron_damping', electron_damping
    endif

    !----------------------------------------
    !     &IONS
    !----------------------------------------


    ! ION_DYNAMICS   
    ! ion_dynamics = 'default' | 'sd' | 'cg' | 'damp' | 'md' | 'none' | 'diis' 
    if (event_ion_dynamics(event_index)) then
       ion_dynamics= rule_ion_dynamics(event_index)
      tdampions_       = .FALSE.
      tconvthrs%active = .FALSE.
      tconvthrs%nstep  = 1
      tconvthrs%ekin   = 0.0d0
      tconvthrs%derho  = 0.0d0
      tconvthrs%force  = 0.0d0

       select case ( ion_dynamics ) 
       case ('SD')
          tsdp = .true.
          tfor = .true.
          fricp= 0.d0
          tconvthrs%ekin   = ekin_conv_thr
          tconvthrs%derho  = etot_conv_thr
          tconvthrs%force  = forc_conv_thr
          tconvthrs%active = .TRUE.
          tconvthrs%nstep  = 1
       case ('VERLET')
          tsdp = .false.
          tfor = .true.
          fricp= 0.d0
       case ('DAMP')
          ldamped    = .TRUE.
          tsdp = .false.
          tfor = .true.
          tdampions_ = .TRUE.
          fricp= ion_damping
          tconvthrs%ekin   = ekin_conv_thr
          tconvthrs%derho  = etot_conv_thr
          tconvthrs%force  = forc_conv_thr
          tconvthrs%active = .TRUE.
          tconvthrs%nstep  = 1
       case ('NONE')
          tsdp = .false.
          tfor = .false.
          fricp= 0.d0
       case default
          call auto_error(' iosys ',' unknown ion_dynamics '//trim(ion_dynamics) )
       end select

       write(*,*) 'RULE EVENT: ion_dynamics', ion_dynamics
    endif


    ! ION_DAMPING   
    if (event_ion_damping(event_index)) then
       ! meaningful only if " ion_dynamics = 'damp' "
       ion_damping = rule_ion_damping(event_index)
       write(*,*) 'RULE EVENT: ion_damping', ion_damping
    endif


    ! ION_TEMPERATURE 
    if (event_ion_temperature(event_index)) then
       ion_temperature = rule_ion_temperature(event_index)
      tcp      = .FALSE.
      tnosep   = .FALSE.
      tolp     = tolp
       select case ( ion_temperature ) 
          !         temperature control of ions via nose' thermostat
          !         tempw (real(DP))  frequency (in which units?)
          !         fnosep (real(DP))  temperature (in which units?)
       case ('NOSE')
          tnosep = .true.
          tcp = .false.
       case ('NOT_CONTROLLED')
          tnosep = .false.
          tcp = .false.
       case ('RESCALING' )
          tnosep = .false.
          tcp = .true.
       case default
          call auto_error(' iosys ',' unknown ion_temperature '//trim(ion_temperature) )
       end select

       write(*,*) 'RULE EVENT: ion_temperature', ion_temperature

    endif

    ! TEMPW
    if (event_tempw(event_index)) then
       tempw  = rule_tempw(event_index)
       ! The follwiong is a required side effect
       ! when resetting tempw
       CALL ions_nose_init( tempw, fnosep, nhpcl, nhptyp, ndega, nhgrp, fnhscl)
       write(*,*) 'RULE EVENT: tempw', tempw
    endif

    !----------------------------------------
    !     &CELL
    !----------------------------------------

    !----------------------------------------
    !     &PHONON
    !----------------------------------------


  END SUBROUTINE employ_rules



  !-----------------------------------------------------------------------
  ! PILOT
  !
  ! Here is the main pilot routine called in CPR, at the top
  ! of the basic dynamics loop just after nose hoover update 
  !-----------------------------------------------------------------------
  subroutine pilot (nfi)
    USE parser, ONLY: parse_unit
    USE io_global, ONLY: ionode, ionode_id
    USE mp,        ONLY : mp_bcast, mp_barrier
    IMPLICIT NONE
    INTEGER :: nfi
    LOGICAL :: file_p
    CHARACTER (LEN=256) :: mbfile = "pilot.mb"

    ! Dynamics Loop Started
    pilot_p   = .TRUE.

    ! This is so we can usurp the exiting parser
    ! that defaults to stdin (unit=5)
    ! We have to do it this way if we are to 
    ! call (reuse) the card_autopilot that is called
    ! by read_cards
    parse_unit = pilot_unit

    ! Our own local for nfi
    current_nfi = nfi

    ! Great for Debugging
    !IF( ionode ) THEN
    !write(*,*)
    !write(*,*) '========================================'
    !write(*,*) 'Autopilot (Dynamic Rules) Implementation'
    !write(*,*) '  CURRENT_NFI=', current_nfi
    !write(*,*) '  event_index=', event_index
    !write(*,*) '  event_step==', event_step(event_index)
    !write(*,*) '========================================'
    !write(*,*)
    !call flush_unit(6)
    !END IF


    ! This allows one pass. Calling parse_mailbox will either:
    ! 1) call init_auto_pilot, which will always set this modules global PAUSE_P variable to FALSE
    ! 2) detect a pause indicator, setting PAUSE_P to TRUE until a new mailbox overrides.
    pause_loop: do

       file_p = .FALSE.
       IF ( ionode ) INQUIRE( FILE = TRIM( mbfile ), EXIST = file_p )
       call mp_bcast(file_p, ionode_id)     

       IF ( file_p ) THEN

          WRITE(*,*) 
          WRITE(*,*) 'Pilot: Mailbox Found!'
          WRITE(*,*) '       CURRENT_NFI=', current_nfi
          call flush_unit(6)

          ! Open the mailbox
          IF ( ionode ) OPEN( UNIT = pilot_unit, FILE = TRIM( mbfile ) )

          ! Will reset PAUSE_P to false unless there is a PAUSE cmd
          ! The following call is MPI safe! It only generates side effects
          CALL parse_mailbox()
          !call mp_barrier()
          WRITE(*,*) 'return from parse_mailbox' 

          ! Perhaps instead of deleting move the file as an input log     
          IF( ionode ) CLOSE( UNIT = pilot_unit, STATUS = 'DELETE' )

       END IF

       IF( .NOT. pause_p ) THEN
          EXIT pause_loop
       ELSE
          write(*,*) 'SLEEPING .... send another pilot.mb'
          call sleep (5)
       END if

    end do pause_loop

    ! Autopilot (Dynamic Rules) Implementation
    ! When nfi has passed (is greater than
    ! the next event, then employ rules
    ! Mailbox may have issued several rules
    ! Attempt to catch up! 
    do while (current_nfi >= event_step(event_index) )         

       write(*,*) 'in while: event_index ', event_index 
       call employ_rules()
       call mp_barrier()

       ! update event_index to current
       event_index = event_index + 1
       write(*,*) 'in while after: event_index ', event_index 

    enddo

  end subroutine pilot

END MODULE cp_autopilot

