! autopilot.f90
!********************************************************************************
! autopilot.f90 				Copyright (c) 2005 Targacept, Inc.
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
 

MODULE autopilot
  !---------------------------------------------------------------------------
  !
  ! This module handles the Autopilot Feature Suite
  ! Written by Lee Atkinson, with help from the ATP team at Targacept, Inc 
  ! Created June 2005
  !
  !   The address for Targacept, Inc. is:
  !     200 East First Street, Suite
  !     300, Winston-Salem, North Carolina 27101; 
  !     Attn: Molecular Design.
  !
  ! See autopilot.html in the Doc directory for more information.
  !---------------------------------------------------------------------------

  USE kinds
  USE parser, ONLY :  read_line

  IMPLICIT NONE
  SAVE

  INTEGER, parameter :: MAX_INT = huge(1)  
  INTEGER, parameter :: max_event_step = 32  !right now there can be upto 32 Autopilot Events
  INTEGER, parameter :: n_auto_vars = 10     !right now there are only 10 Autopilot Variables

  INTEGER   :: n_events 
  INTEGER   :: event_index = 0
  INTEGER   :: max_rules = 320 !(max_event_step * n_auto_vars)
  INTEGER   :: n_rules 
  INTEGER   :: event_step(max_event_step)
  INTEGER   :: current_nfi
  LOGICAL   :: pilot_p   = .FALSE.    ! pilot property
  LOGICAL   :: restart_p = .FALSE.    ! restart property
  LOGICAL   :: pause_p   = .FALSE.    ! pause   property
  INTEGER   :: pilot_unit = 42   ! perhaps move this to io_files
  CHARACTER(LEN=256) :: pilot_type

  ! AUTOPILOT VARIABLES
  ! These are the variable tables which change the actual variable
  ! dynamically during the course of a simulation. There are many 
  ! parameters which govern a simulation, yet only these are allowed 
  ! to be changed using the event rule mechanism.
  ! Each of these tables are ytped according to their variable 
  ! and begin with event_

  !     &CONTROL
  INTEGER   :: rule_isave(max_event_step)
  INTEGER   :: rule_iprint(max_event_step)
  REAL(DP) :: rule_dt(max_event_step)
  !     &SYSTEM

  !     &ELECTRONS
  REAL(DP)         :: rule_emass(max_event_step)
  CHARACTER(LEN=80) :: rule_electron_dynamics(max_event_step)
  REAL(DP)         :: rule_electron_damping(max_event_step)


  !     &IONS
  CHARACTER(LEN=80) :: rule_ion_dynamics(max_event_step)
  REAL(DP)         :: rule_ion_damping(max_event_step)
  CHARACTER(LEN=80) :: rule_ion_temperature(max_event_step)

  REAL(DP) :: rule_tempw(max_event_step)
  !     &CELL

  !     &PHONON


  ! Each rule also needs to be correlated (flagged) against the event time table
  ! This is used to flag the a given variable is to be changed on an
  ! event. Initially all set to false, a rule against an event makes it true
  ! Each of these flags are LOGICALs and begin with event_
  !     &CONTROL
  LOGICAL :: event_isave(max_event_step)           
  LOGICAL :: event_iprint(max_event_step)          
  LOGICAL :: event_dt(max_event_step)              
  !     &SYSTEM

  !     &ELECTRONS
  LOGICAL :: event_emass(max_event_step)   
  LOGICAL :: event_electron_dynamics(max_event_step)   
  LOGICAL :: event_electron_damping(max_event_step)

  !     &IONS
  LOGICAL :: event_ion_dynamics(max_event_step)   
  LOGICAL :: event_ion_damping(max_event_step)
  LOGICAL :: event_ion_temperature(max_event_step)   

  LOGICAL :: event_tempw(max_event_step)           
  !     &CELL

  !     &PHONON


  PRIVATE
  PUBLIC :: auto_check, init_autopilot, card_autopilot, pilot, employ_rules, add_rule, & 
       & assign_rule, restart_p, max_event_step, event_index, event_step, rule_isave, &
       & rule_iprint, rule_dt, rule_emass, rule_electron_dynamics, rule_electron_damping, &
       & rule_ion_dynamics, rule_ion_damping, rule_ion_temperature, rule_tempw, &
       & event_isave, event_iprint, event_dt, event_emass, &
       & event_electron_dynamics, event_electron_damping, event_ion_dynamics, &
       & event_ion_damping, event_ion_temperature, event_tempw

CONTAINS

  !----------------------------------------------------------------------------
  SUBROUTINE auto_error( calling_routine, message)
    !----------------------------------------------------------------------------
    ! This routine calls errore based upon the pilot property flag.
    ! If the flag is TRUE, the simulation will not stop,
    !    but the pause property flag is set to TRUE,
    !    causing pilot to force a state of sleep, 
    !    until the user can mail proper commands.
    ! Otherwise, its assumed that dynamics have not started
    ! and this call is an indirect result of read_cards.
    ! Thus the simulation will stop.
    ! Either way errore will always issues a warning message.

    USE io_global, ONLY : ionode_id
    USE mp,        ONLY : mp_bcast

    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: calling_routine, message
    ! the name of the calling calling_routinee
    ! the output message
    INTEGER :: ierr
    ! the error flag

    IF (pilot_p) THEN
       ! if ierr < 0 errore writes the message but does not stop
       ierr = -1
       pause_p = .TRUE.
       !call mp_bcast(pause_p, ionode_id)
    ELSE
       ! if ierr > 0 it stops
       ierr = 1
    ENDIF

    CALL errore( calling_routine, message, ierr )

  END SUBROUTINE auto_error


  !-----------------------------------------------------------------------
  ! AUTO (restart) MODE
  !
  ! Checks if restart files are present, just like what readfile_cp does, 
  ! which calls cp_readfile which create a path to restart.xml. 
  ! This could be checked with inquire, as in check_restartfile. 
  ! If restart_mode=auto, and restart.xml is present, 
  ! then restart_mode="restart", else its "from_scratch".
  ! Set other associated vars, appropriately.
  !
  ! Put this in setcontrol_flags on the select statement
  !-----------------------------------------------------------------------
  LOGICAL FUNCTION auto_check (ndr, scradir)
    USE io_global, ONLY: ionode, ionode_id
    USE mp,        ONLY : mp_bcast
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ndr    !  I/O unit number
    CHARACTER(LEN=*), INTENT(IN) :: scradir
    CHARACTER(LEN=256) :: dirname, filename
    CHARACTER (LEN=6), EXTERNAL :: int_to_char
    LOGICAL :: restart_p = .FALSE.
    INTEGER :: strlen
    ! right now cp_readfile is called with scradir = ' '
    ! so, in keeping with the current design,
    ! the responsibility of setting 
    ! ndr and scradir is the calling program


    IF (ionode) THEN
       dirname = 'RESTART' // int_to_char( ndr )
       IF ( LEN( scradir ) > 1 ) THEN
          strlen  = index(scradir,' ') - 1
          dirname = scradir(1:strlen) // '/' // dirname
       END IF

       filename = TRIM( dirname ) // '/' // 'restart.xml'
       INQUIRE( FILE = TRIM( filename ), EXIST = restart_p )

       auto_check = restart_p
    END IF
    CALL mp_bcast(auto_check, ionode_id)

    return

  END FUNCTION auto_check


  !-----------------------------------------------------------------------
  ! INITIALIZE AUTOPILOT
  !
  ! Must be done, even if not in use.
  ! Add this as a call in read_cards.f90 sub: card_default_values
  !-----------------------------------------------------------------------
  SUBROUTINE init_autopilot
    IMPLICIT NONE
    integer :: event

    pause_p = .FALSE.

    ! Initialize all events to an iteration that should never occur
    DO event=1,max_event_step
       event_step(event) = MAX_INT
    ENDDO

    n_events = 0 
    n_rules  = 0
    event_index = 1

    ! Initialize here
    !     &CONTROL
    event_isave(:)            = .false.
    event_iprint(:)           = .false.
    event_dt(:)               = .false.
    !     &SYSTEM
    !     &ELECTRONS
    event_emass(:)             = .false.
    event_electron_dynamics(:) = .false.
    event_electron_damping(:)  = .false.

    !     &IONS
    event_ion_dynamics(:)      = .false.
    event_ion_damping(:)       = .false.
    event_ion_temperature(:)   = .false.
    event_tempw(:)             = .false.
    !     &CELL
    !     &PHONON

    rule_isave(:)             = 0
    rule_iprint(:)            = 0
    rule_dt(:)                = 0.0
    rule_emass(:)             = 0.0
    rule_electron_dynamics(:) = 'NONE'
    rule_electron_damping(:)  = 0.0
    rule_ion_dynamics(:)      = 'NONE'
    rule_ion_damping(:)       = 0.0
    rule_ion_temperature(:)   = 'NOT_CONTROLLED'
    rule_tempw(:)             = 0.01

  END SUBROUTINE init_autopilot



  !-----------------------------------------------------------------------
  ! subroutine CARD_AUTOPILOT
  !
  ! called in READ_CARDS and in PARSE_MAILBOX
  !-----------------------------------------------------------------------  
  SUBROUTINE card_autopilot( input_line )
    USE io_global, ONLY: ionode, ionode_id
    USE mp,        ONLY : mp_bcast
    IMPLICIT NONE
    INTEGER :: i, j, linelen
    CHARACTER(LEN=256) :: input_line
    LOGICAL            :: process_this_line = .FALSE.
    LOGICAL            :: endrules = .FALSE.
    LOGICAL            :: tend = .FALSE.
    LOGICAL, SAVE      :: tread = .FALSE.
    LOGICAL, EXTERNAL  :: matches
    CHARACTER(LEN=1), EXTERNAL :: capital

    ! This is so we do not read a autopilot card twice from the input file
    IF (( .NOT. pilot_p ) .and. tread ) THEN
       CALL errore( ' card_autopilot  ', ' two occurrence ', 2 )
    END IF

    ! If this routined has been called from parse_mailbox
    ! the pilot_type should be set
    IF ( pilot_p ) THEN
       ! IF its MANUAL then process this line first! 
       ! other skip this line and get the next
       IF (TRIM(pilot_type) .eq. 'MANUAL') THEN
          process_this_line = .TRUE.
       ELSE IF (TRIM(pilot_type) .eq. 'PILOT') THEN
          process_this_line = .FALSE.
       ELSE IF (TRIM(pilot_type) .eq. 'AUTO') THEN
          process_this_line = .FALSE.
       ELSE
          WRITE(*,*) 'AUTOPILOT: UNRECOGNIZED PILOT TYPE!', TRIM(pilot_type), '===='
          GO TO 10
       END IF
    ELSE
       ! this routine is called from read_cards
       pilot_type = 'AUTO'
       process_this_line = .FALSE.
    END IF

    j=0
    ! must use a local (j) since n_rules may not get updated correctly
    DO WHILE ((.NOT. endrules) .and. j<=max_rules)
       j=j+1

       IF (j > max_rules) THEN
          CALL auto_error( ' AutoPilot ','Maximum Number of Dynamic Rules May Have Been Execced!')
          go to 10
       END IF


       WRITE(*,*) 'card_autopilot 1: input_line ', input_line

       ! Assume that pilot_p is an indicator and that
       ! this call to card_autopilot is from parse_mailbox
       ! and not from read_cards
       IF(pilot_p) THEN
          IF ( .NOT. process_this_line ) THEN
             ! get the next one
             CALL read_line( input_line, end_of_file = tend)
          END IF
       ELSE
          ! from read_cards
          CALL read_line( input_line, end_of_file = tend)
       END IF

       linelen = LEN_TRIM( input_line )

       DO i = 1, linelen
          input_line( i : i ) = capital( input_line( i : i ) )
       END DO

       WRITE(*,*) 'card_autopilot 2: input_line ', input_line

       ! If ENDRULES isnt found, add_rule will fail
       ! and we run out of line anyway
       IF ( tend .or. matches( 'ENDRULES', input_line ) ) GO TO 10

       ! Assume this is a rule
       write(*,*) 'about to add_rule: input_line ', input_line
       CALL ADD_RULE(input_line)
       ! now, do not process this line anymore
       ! make sure we get the next one
       process_this_line = .FALSE.

    END DO

    WRITE(*,*) 'AUTOPILOT SET'

10  CONTINUE

  END SUBROUTINE card_autopilot





  !-----------------------------------------------------------------------
  ! ADD RULE
  !-----------------------------------------------------------------------
  SUBROUTINE add_rule( input_line )
    USE io_global, ONLY: ionode, ionode_id
    !USE mp,        ONLY : mp_bcast
    IMPLICIT NONE
    integer :: i, j, linelen
    integer :: eq1_pos, eq2_pos, plus_pos, colon_pos
    CHARACTER(LEN=256) :: input_line
    CHARACTER(LEN=32)  :: var_label
    CHARACTER(LEN=32)  :: value_str
    INTEGER            :: on_step, now_step, plus_step
    integer            :: ios
    integer            :: event

    LOGICAL, SAVE      :: tread = .FALSE.

    LOGICAL, EXTERNAL  :: matches
    CHARACTER(LEN=1), EXTERNAL :: capital


    ! this is a temporary local variable
    event = n_events

    ! important for parsing
    i=0
    j=0
    eq1_pos   = 0
    eq2_pos   = 0
    plus_pos  = 0
    colon_pos = 0

    linelen=LEN_TRIM( input_line )

    write(*,*) 'ADD_RULE: pilot_type', pilot_type


    ! Attempt to get PLUS SYMBOL
    i = 1
    do while( (plus_pos .eq. 0) .and. (i <= linelen) )
       i = i + 1
       if(input_line( i : i ) .eq. '+') then
          plus_pos = i
       endif
    end do

    ! Attempt to get a colon
    i = 1
    do while( (colon_pos .eq. 0) .and. (i <= linelen) )
       i = i + 1
       if(input_line( i : i ) .eq. ':') then
          colon_pos = i
       endif
    end do

    ! Attempt to get first equals
    i = 1
    do while( (eq1_pos .eq. 0) .and. (i <= linelen) )
       i = i + 1
       if(input_line( i : i ) .eq. '=') then
          eq1_pos = i
       endif
    end do


    ! Attempt to get second equals
    if ((eq1_pos .ne. 0) .and. (eq1_pos < colon_pos)) then
       i = colon_pos + 1
       do while( (eq2_pos .eq. 0) .and. (i <= linelen) )
          i = i + 1
          if(input_line( i : i ) .eq. '=') then
             eq2_pos = i
          endif
       end do
    endif

    ! Complain if there is a bad parse
    if (colon_pos .eq. 0) then
       call auto_error( ' AutoPilot ','Missing colon separator')
       go to 20
    endif

    if (eq1_pos .eq. 0) then
       call auto_error( ' AutoPilot ','Missing equals sign')
       go to 20
    endif

    if ((plus_pos > 0) .and. (eq1_pos < colon_pos)) then
       call auto_error( ' AutoPilot ','equals and plus found prior to colon')
       go to 20
    endif


    !================================================================================
    ! Detect events
    IF ( (pilot_type .eq. 'MANUAL') .or. (pilot_type .eq. 'PILOT') ) THEN
       !-------------------------------------------
       !Do the equivalent of the following:
       !READ(input_line, *) now_label, plus_label1, plus_step, colon_label, var_label, eq_label2, value_str
       !Format:
       ! [NOW [+ plus_step] :] var_label = value_str
       !-------------------------------------------

       write(*,*) 'ADD_RULE: MANUAL STEERING'

       ! if there is a NOW, get it and try to get plus and plus_step

       IF ( matches( 'NOW', input_line ) ) THEN
          ! Attempt to get PLUS_STEP
          plus_step = 0
          ! if all is non-trivial, read a PLUS_STEP
          if ((plus_pos > 0) .and. (colon_pos > plus_pos)) then
             ! assume a number lies between
             read(input_line(plus_pos+1:colon_pos-1),*,iostat=ios) plus_step
             if(ios .ne. 0) then
                CALL auto_error( ' AutoPilot ','Value Type Mismatch on NOW line!')
                go to 20
             end if
          end if
          ! set NOW event
          now_step = current_nfi + plus_step
       ELSE
          ! set NOW event
          now_step = current_nfi
       END IF


       !================================================================================
       ! set event
       !
       ! Heres where it get interesting
       ! We may have a new event , or not! :)

       IF ( ((event-1) .gt. 0) .and. ( now_step .lt. event_step(event-1)) ) THEN
          write(*,*) ' AutoPilot: current input_line', input_line 
          CALL auto_error( ' AutoPilot ','Dynamic Rule Event Out of Order!')
          go to 20
       ENDIF

       IF ( (event .eq. 0) .or. ( now_step .gt. event_step(event)) ) THEN
          ! new event
          event = event + 1

          IF (event > max_event_step) THEN
             write(*,*) ' AutoPilot: current input_line', input_line 
             CALL auto_error( ' AutoPilot ','Maximum Number of Dynamic Rule Event Exceeded!')
             go to 20
          ENDIF

          event_step(event) = now_step
          n_events = event       
       ENDIF


    ELSE IF ( matches( 'ON_STEP', input_line ) ) THEN
       ! Assuming pilot_type is AUTO
       ! if it isnt and ON_STEP these rules wont take anyway

       !-------------------------------------------
       !Do the equivalent of the following:
       !READ(input_line, *) on_step_label, eq_label1, on_step, colon_label, var_label, eq_label2, value_str
       !Format:
       ! ON_STEP = on_step : var_label = value_str
       !-------------------------------------------

       write(*,*) 'ADD_RULE: POWER STEERING'

       ! Attempt to get ON_STEP
       on_step = MAX_INT
       ! if all is non-trivial, read a PLUS_STEP
       if ((eq1_pos > 0) .and. (colon_pos > eq1_pos)) then
          ! assume a number lies between
          read(input_line(eq1_pos+1:colon_pos-1),*,iostat=ios) on_step
          if(ios .ne. 0) then
             CALL auto_error( ' AutoPilot ','Value Type Mismatch on ON_STEP line!')
             go to 20
          end if
       end if
       


       !================================================================================
       ! set event
       !
       ! Heres where it get interesting
       ! We may have a new event , or not! :)       


       IF ( ((event-1) .gt. 0) .and. ( on_step .lt. event_step(event-1)) ) THEN
          write(*,*) ' AutoPilot: current input_line', input_line 
          CALL auto_error( ' AutoPilot ','Dynamic Rule Event Out of Order!')
          go to 20
       ENDIF


       IF ( (event .eq. 0) .or. (on_step .gt. event_step(event)) ) THEN
          ! new event
          event = event + 1
          IF (event > max_event_step) THEN
             write(*,*) ' AutoPilot: current input_line', input_line
             CALL auto_error( ' AutoPilot ','Maximum Number of Dynamic Rule Event Exceeded!')
             go to 20
          ENDIF
          event_step(event) = on_step
          n_events = event       
       ENDIF

    END IF ! Event Detection Complete


    !-------------------------------------
    ! Now look for a label and a value
    !-------------------------------------

    if (eq2_pos .eq. 0) then
       var_label = input_line(colon_pos+1: eq1_pos-1)    
       read( input_line(eq1_pos+1:linelen), *, iostat=ios) value_str  
       if(ios .ne. 0) then
          CALL auto_error( ' AutoPilot ','Value Type Mismatch on NOW_STEP line!')
          go to 20
       end if
    else
       var_label = input_line(colon_pos+1: eq2_pos-1)    
       read( input_line(eq2_pos+1:linelen), *, iostat=ios) value_str  
       if(ios .ne. 0) then
          CALL auto_error( ' AutoPilot ','Value Type Mismatch on ON_STEP line!')
          go to 20
       end if
    endif

    ! The Assignment must lie outside the new event scope since
    ! there can exists more than one rule per event

    IF ( (n_rules+1) .gt. max_rules) THEN
       write(*,*) ' AutoPilot: current n_rules', n_rules
       CALL auto_error( ' AutoPilot ', ' invalid number of rules ')
       go to 20
    END IF

    call assign_rule(event, var_label, value_str)    

    write(*,*) 'n_rules=', n_rules

    CALL flush_unit(6)

20  CONTINUE

  END SUBROUTINE add_rule


  !-----------------------------------------------------------------------
  ! ASSIGN_RULE
  !-----------------------------------------------------------------------
  SUBROUTINE assign_rule(event, var, value)
    USE input_parameters, ONLY : isave, iprint, dt, tempw
    USE io_global, ONLY: ionode, ionode_id
    IMPLICIT NONE
    INTEGER :: i, event, varlen
    CHARACTER(LEN=32) :: var
    CHARACTER(LEN=32) :: value
    INTEGER   :: int_value
    REAL      :: real_value
    REAL(DP) :: realDP_value
    LOGICAL   :: logical_value
    LOGICAL   :: assigned
    LOGICAL, EXTERNAL  :: matches
    CHARACTER(LEN=1), EXTERNAL :: capital


    var = TRIM(var)
    varlen = LEN_TRIM(var)

    DO i = 1, varlen
       var( i : i ) = capital( var( i : i ) )
    END DO


    write(*,*) 'ASSIGNING RULE: event var value', event, var, value
    assigned = .TRUE.

    IF ( matches( "ISAVE", var ) ) THEN
       read(value, *) int_value
       rule_isave(event)  = int_value
       event_isave(event) = .true.
    ELSEIF ( matches( "IPRINT", var ) ) THEN
       read(value, *) int_value
       rule_iprint(event)  = int_value
       event_iprint(event) = .true.
    ELSEIF ( matches( "DT", var ) ) THEN
       read(value, *) real_value
       rule_dt(event)  = real_value
       event_dt(event) = .true.
       write(*,*) 'RULE_DT', rule_dt(event), 'EVENT', event
    ELSEIF ( matches( "EMASS", var ) ) THEN
       read(value, *) realDP_value
       rule_emass(event)  = realDP_value
       event_emass(event) = .true.
    ELSEIF ( matches( "ELECTRON_DYNAMICS", var ) ) THEN
       read(value, *) value
       if ((value .ne. 'SD') .and. (value .ne. 'VERLET') .and. (value .ne. 'DAMP') .and. (value .ne. 'NONE')) then
          call auto_error(' autopilot ',' unknown electron_dynamics '//trim(value) )
          assigned = .FALSE.
          go to 40
       endif
       rule_electron_dynamics(event)  = value
       event_electron_dynamics(event) = .true.
    ELSEIF ( matches( "ELECTRON_DAMPING", var ) ) THEN
       read(value, *) realDP_value
       rule_electron_damping(event)  = realDP_value
       event_electron_damping(event) = .true.
    ELSEIF ( matches( "ION_DYNAMICS", var ) ) THEN
       read(value, *) value
       if ((value .ne. 'SD') .and. (value .ne. 'VERLET') .and. (value .ne. 'DAMP') .and. (value .ne. 'NONE')) then
          call auto_error(' autopilot ',' unknown ion_dynamics '//trim(value) )
          assigned = .FALSE.
          go to 40
       endif
       rule_ion_dynamics(event)  = value
       event_ion_dynamics(event) = .true.
    ELSEIF ( matches( "ION_DAMPING", var ) ) THEN
       read(value, *) realDP_value
       rule_ion_damping(event)  = realDP_value
       event_ion_damping(event) = .true.
    ELSEIF ( matches( "ION_TEMPERATURE", var ) ) THEN
       read(value, *) value
       if ((value .ne. 'NOSE') .and. (value .ne. 'NOT_CONTROLLED') .and. (value .ne. 'RESCALING')) then
          call auto_error(' autopilot ',' unknown ion_temperature '//trim(value) )
          assigned = .FALSE.
          go to 40
       endif
       rule_ion_temperature(event)  = value
       event_ion_temperature(event) = .true.
    ELSEIF ( matches( "TEMPW", var ) ) THEN
       read(value, *) realDP_value
       rule_tempw(event)  = realDP_value
       event_tempw(event) = .true.
    ELSE
       CALL auto_error( 'autopilot', ' ASSIGN_RULE: FAILED  '//trim(var)//' '//trim(value) )
    END IF

40  if (assigned) then
       n_rules   = n_rules + 1
       write(*,*) 'Autopilot: Rule Assigned ', n_rules
    else
       write(*,*) 'Autopilot: Rule Assignment Failure '
       CALL auto_error( 'autopilot', ' ASSIGN_RULE: FAILED  '//trim(var)//' '//trim(value) )
    endif

  END SUBROUTINE assign_rule


  !-----------------------------------------------------------------------
  ! EMPLOY_RULES
  !-----------------------------------------------------------------------
  SUBROUTINE employ_rules()
    USE input_parameters, ONLY : isave, iprint, dt, & 
         & emass, electron_dynamics, electron_damping, &
         & ion_dynamics, ion_damping, &
         & ion_temperature, fnosep, nhpcl, ndega, nat
    use ions_nose, ONLY: tempw
    USE control_flags, only: tsde, tsdp, tfor, tcp, tnosep
    use wave_base, only: frice
    use ions_base, only: fricp
    USE ions_nose, ONLY: ions_nose_init
    USE io_global, ONLY: ionode, ionode_id
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
       select case ( electron_dynamics ) 
       case ('SD')
          tsde  = .true.
          frice = 0.d0
       case ('VERLET')
          tsde  = .false.
          frice = 0.d0
       case ('DAMP')
          tsde  = .false.
          frice = electron_damping
       case ('NONE')
          tsde  = .false.
          frice = 0.d0
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
       select case ( ion_dynamics ) 
       case ('SD')
          tsdp = .true.
          tfor = .true.
          fricp= 0.d0
       case ('VERLET')
          tsdp = .false.
          tfor = .true.
          fricp= 0.d0
       case ('DAMP')
          tsdp = .false.
          tfor = .true.
          fricp= ion_damping
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
       CALL ions_nose_init( tempw, fnosep, nhpcl, ndega, nat )
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
  ! PARSE_MAILBOX
  !
  ! Read the mailbox with a mailbox parser
  ! if it starts with ON_STEP, then apply to event table etc
  ! if not the try to establish that its a variable to set right now
  !-----------------------------------------------------------------------
  SUBROUTINE parse_mailbox ()
    !use ifport, only: sleep
    USE io_global, ONLY: ionode, ionode_id
    USE mp,        ONLY : mp_bcast, mp_barrier
    IMPLICIT NONE
    INTEGER :: i,j
    CHARACTER(LEN=256) :: input_line
    CHARACTER(LEN=80)  :: card
    LOGICAL            :: tend

    CHARACTER(LEN=1), EXTERNAL :: capital
    LOGICAL, EXTERNAL  :: matches


    ! we can use this parser routine, since parse_unit=pilot_unit
    CALL read_line( input_line, end_of_file=tend )
    IF (tend) GO TO 50

    DO i = 1, LEN_TRIM( input_line )
       input_line( i : i ) = capital( input_line( i : i ) )
    END DO

    ! This conditional implements the PAUSE feature calling init_auto_pilot, 
    ! will reset this modules global PAUSE_P variable to FALSE
    IF ( matches( "PAUSE", input_line ) .or. &
         matches( "SLEEP", input_line ) .or. &
         matches( "HOVER", input_line ) .or. &
         matches( "WAIT",  input_line ) .or. &
         matches( "HOLD",  input_line ) ) THEN

       write(*,*) 'SLEEPING'
       write(*,*) 'INPUT_LINE=', input_line
       pause_p = .TRUE.

       ! Now just quit this subroutine

    ELSE
       ! Also, We didnt see a PAUSE cmd!
       pause_p = .FALSE.

       ! now lets detect the mode for card_autopilot
       ! even though this line will be passed to it the first time

       IF ( matches( "AUTOPILOT", TRIM(input_line) ) ) THEN
          WRITE(*,*) 'NEW AUTOPILOT COURSE DETECTED' 
          pilot_type ='AUTO'
       ELSE IF (matches( "PILOT", TRIM(input_line) ) ) THEN
          WRITE(*,*) 'RELATIVE PILOT COURSE CORRECTION DETECTED'
          pilot_type ='PILOT'
       ELSE IF (matches( "NOW", TRIM(input_line) ) ) THEN
          WRITE(*,*) 'MANUAL PILOTING DETECTED'
          pilot_type ='MANUAL'
       ELSE
          ! Well lets just pause since this guys is throwing trash
          WRITE(*,*) 'MAILBOX CONTENTS NOT UNDERSTOOD: pausing'
          pause_p = .TRUE.
       ENDIF

    END IF

    IF (pause_p) GO TO 50


    ! ok if one adds a rule during steering`
    ! event table must be cleared (from steer point) forward
    !
    ! Every nodes gets this (and the call to card_autopilot
    ! which calls add_rule, which calls assign_rule, etc
    ! In this way we sync the event table
    ! Then we shouldn't have to sync employ_rules variable
    ! changes, or their subroutine side effects (like ions_nose_init)  

    CALL init_autopilot()


    WRITE(*,*) 'parse_mailbox: about to call card_autopilot: pilot_type', pilot_type
    write(*,*) 'input_line=', input_line    
    CALL card_autopilot( input_line )


    ! this is needed just befor we end this subroutine
50  CONTINUE
    !call mp_barrier()    

    WRITE(*,*) 'end of parse' 

  end subroutine parse_mailbox


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

END MODULE autopilot

