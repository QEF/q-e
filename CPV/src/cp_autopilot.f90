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

!********************************************************************************
! 04/2018
! Implemented switch from CG to Verlet and from Verlet to CG.
! The code automagically set the correct wfc velocity after the last CG step,
! and inizialize again the mass of the electron mu(k) after the CG.
! Note that CG modifies the array of the electron mass.
!
! Riccardo Bertossa, Federico Grasselli
!   (SISSA - via Bonomea, 265 - 34136 Trieste ITALY)
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
  USE autopilot, ONLY : current_nfi, pilot_p, pilot_unit, pause_p,auto_error, &
        &  parse_mailbox, rule_isave, rule_iprint, rule_dt, rule_emass,       &
        &  rule_electron_dynamics, rule_electron_damping, rule_ion_dynamics, &
        &  rule_ion_damping,  rule_ion_temperature, rule_tempw, &
        &  rule_electron_orthogonalization, rule_tprint
  USE autopilot, ONLY : event_index, event_step, event_isave, event_iprint, &
        &  event_dt, event_emass, event_electron_dynamics, event_electron_damping, &
        &  event_ion_dynamics, event_ion_damping, event_ion_temperature, event_tempw, &
        & event_electron_orthogonalization, &
        & event_tprint

  IMPLICIT NONE
  SAVE


  PRIVATE
  PUBLIC ::  pilot, employ_rules
  LOGICAL, PRIVATE :: had_tcg_true = .false.
  LOGICAL, PRIVATE :: had_tens_true = .false.
CONTAINS

  


  !-----------------------------------------------------------------------
  ! EMPLOY_RULES
  !-----------------------------------------------------------------------
  SUBROUTINE employ_rules()
    USE input_parameters, ONLY :   dt, & 
         &  electron_dynamics, electron_damping, &
         & ion_dynamics, ion_damping, &
         & ion_temperature, fnosep, nhpcl, nhptyp, nhgrp, fnhscl, ndega, nat, &
         & orthogonalization
    use ions_nose, ONLY: tempw
    USE control_flags, only: tsde, tsdp, tfor, tcp, tnosep, isave,iprint,&
                             tconvthrs, tolp, &
                             ekin_conv_thr, forc_conv_thr, etot_conv_thr,&
                             tortho, tfirst, tlast, tprint
    use wave_base, only: frice
    use ions_base, only: fricp
    USE ions_nose, ONLY: ions_nose_init
    USE io_global, ONLY: ionode, ionode_id
    USE time_step,                ONLY : set_time_step
    USE cp_electronic_mass,       ONLY: emass, emass_cutoff, emass_precond
    USE cg_module,                ONLY : tcg,allocate_cg,cg_info, &
                                         nfi_firstcg,c0old
    USE wavefunctions_module,     ONLY : cm_bgrp
    USE ensemble_dft,             ONLY : tens,allocate_ensemble_dft
    USE uspp,                     ONLY : nkb, nkbus
    USE electrons_base,           ONLY : nspin, nbsp, nbspx, nudx
    USE gvecw,                    ONLY : ngw
    USE fft_base,                 ONLY : dffts
    USE cp_main_variables,        ONLY : descla
    USE ions_base,                ONLY : nat_ions_base => nat 
    USE cell_base,                ONLY : tpiba2
    USE cp_main_variables,        ONLY : ema0bg
    USE gvecw,                    ONLY : g2kin, ngw 

    IMPLICIT NONE

    !----------------------------------------
    !     &CONTROL
    !----------------------------------------

    ! ISAVE
    if (event_isave(event_index)) then
       isave            = rule_isave(event_index)
       IF ( ionode ) write(*,'(4X,A,15X,I10)') 'Rule event: isave', isave
    endif

    ! IPRINT
    if (event_iprint(event_index)) then
       iprint           = rule_iprint(event_index)
       IF ( ionode ) write(*,'(4X,A,13X,I10)') 'Rule event: iprint', iprint
    endif

    ! TPRINT
    if (event_tprint(event_index)) then
       tprint           = rule_tprint(event_index)
       IF ( ionode ) write(*,*) 'NUOVO Rule event: tprint', tprint
    endif

    ! DT
    if (event_dt(event_index)) then
       dt               = rule_dt(event_index)
       CALL set_time_step( dt )  ! This will make the old velocities wrong!
       IF ( ionode ) write(*,'(4X,A,18X,F10.4)') 'Rule event: dt', dt
    endif

    !----------------------------------------
    !     &SYSTEM
    !----------------------------------------

    !----------------------------------------
    !     &ELECTRONS
    !----------------------------------------
    
    
    if (event_electron_orthogonalization(event_index)) then
       orthogonalization=rule_electron_orthogonalization(event_index)
       select case (orthogonalization)
       case ('ORTHO')
           tortho=.true.
           IF ( ionode ) write(*,*) 'Wow, setting tortho=.true. !' 
       case ('GRAM-SCHMIDT')
           tortho=.false.
           IF ( ionode ) write(*,*) 'Wow, setting tortho=.false. ! (Ma cossa xe sta monada?)' 
       case default 
           call auto_error(' autopilot ',' unknown orthogonalization'//trim(orthogonalization) )
       end select

    endif

    ! EMASS    
    if (event_emass(event_index)) then
       emass            = rule_emass(event_index)
       IF ( ionode ) write(*,'(4X,A,15X,F10.4)') 'Rule event: emass', emass
    endif

    ! ELECTRON_DYNAMICS   
    ! electron_dynamics = 'sd' | 'verlet' | 'damp' | 'none'
    if (event_electron_dynamics(event_index)) then
       electron_dynamics= rule_electron_dynamics(event_index)
       frice = 0.d0
       ! the cg algorithm uses cm for its internal purposes. The true cm is
       ! c0old, so now I have to copy it into cm
       if ( tcg ) then
           cm_bgrp=c0old
           ! tfirst=.true.   !check if this is needed (I don't think so)
           ! the conjugate gradient method modifies the electron mass mu(k),
           ! by calling the routine
           ! emass_precond_tpa( ema0bg, tpiba2, emass_cutoff )
           ! so I call here the standard routine to set mu(k)
           CALL emass_precond( ema0bg, g2kin, ngw, tpiba2, emass_cutoff)

       endif
       select case ( electron_dynamics ) 
       case ('SD')
          tsde  = .true.
          tcg=.false.
       case ('VERLET')
          tsde  = .false.
          tcg=.false.
       case ('DAMP')
          tsde  = .false.
          frice = electron_damping
          tcg=.false.
       case ('NONE')
          tsde  = .false.
          tcg=.false.
       case ('CG')
          IF ( ionode ) write(*,*) 'Wow, setting tcg=.true. at step '&
                 ,current_nfi,' ! (La ghe domandi a mia molie)'
          if (.not. tcg) then
              nfi_firstcg=current_nfi  !the first step is different!
          end if
          if (.not. had_tcg_true) then
              had_tcg_true=.true.
              if (.not. had_tens_true) then
                  CALL allocate_ensemble_dft( nkb, nbsp, ngw, nudx, nspin, nbspx, &
                                 dffts%nnr, nat_ions_base, descla )
              endif
              CALL allocate_cg( ngw, nbspx,nkbus )              
          endif
          CALL cg_info()
          tcg = .true.
       case default
          call auto_error(' autopilot ',' unknown electron_dynamics '//trim(electron_dynamics) )
       end select

       IF ( ionode ) write(*,'(4X,A,2X,A10)') 'Rule event: electron_dynamics', electron_dynamics

    endif


    ! ELECTRON_DAMPING   
    if (event_electron_damping(event_index)) then
       ! meaningful only if " electron_dynamics = 'damp' "
       electron_damping = rule_electron_damping(event_index)
       frice = electron_damping
       IF ( ionode ) write(*,'(4X,A,4X,F10.4)') 'Rule event: electron_damping', electron_damping
    endif

    !----------------------------------------
    !     &IONS
    !----------------------------------------


    ! ION_DYNAMICS   
    ! ion_dynamics = 'sd' | 'verlet' | 'damp' | 'none'
    if (event_ion_dynamics(event_index)) then
      ion_dynamics= rule_ion_dynamics(event_index)
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
          tsdp = .false.
          tfor = .true.
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

    endif


    ! ION_DAMPING   
    if (event_ion_damping(event_index)) then
       ! meaningful only if " ion_dynamics = 'damp' "
       ion_damping = rule_ion_damping(event_index)
       IF ( ionode ) write(*,'(4X,A,9X,F10.4)') 'Rule event: ion_damping', ion_damping
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
          if ( .not. event_tempw(event_index)) then
              IF ( ionode ) write(*,*) 'WARNING: missing tempw event (if undefined can make the code bananas) '
          end if
          if ( .not. tfor) then
              IF ( ionode ) write(*,*) 'WARNING: not doing Verlet on ions? (am I correct?)'
          end if
       case ('NOT_CONTROLLED')
          tnosep = .false.
          tcp = .false.
       case ('RESCALING' )
          tnosep = .false.
          tcp = .true.
       case default
          call auto_error(' iosys ',' unknown ion_temperature '//trim(ion_temperature) )
       end select

       IF ( ionode ) write(*,'(4X,A,5X,A)') 'Rule event: ion_temperature', ion_temperature

    endif

    ! TEMPW
    if (event_tempw(event_index)) then
       tempw  = rule_tempw(event_index)
       ! The follwiong is a required side effect
       ! when resetting tempw
       CALL ions_nose_init( tempw, fnosep, nhpcl, nhptyp, ndega, nhgrp, fnhscl)
       IF ( ionode ) write(*,'(4X,A,15X,F10.4)') 'Rule event: tempw', tempw
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
    USE mp_world,  ONLY : world_comm
#if defined (__NAG)
    USE f90_unix_proc
#endif
    USE ensemble_dft,             ONLY : tens
    USE cg_module, ONLY : tcg
    USE control_flags , only : tprint

    IMPLICIT NONE
    INTEGER :: nfi
    LOGICAL :: file_p
    CHARACTER (LEN=256) :: mbfile = "pilot.mb"

    ! The tcg control variable, that enables conjugate gradient electron
    ! dynamics, allocates some arrays when used in the input file. So I check
    ! if this variable is defined, and set had_tcg_true to true.

    if (tcg) had_tcg_true = .true.
    if (tens) had_tens_true = .true.

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


    ! This allows one pass. Calling parse_mailbox will either:
    ! 1) call init_auto_pilot, which will always set this modules global PAUSE_P variable to FALSE
    ! 2) detect a pause indicator, setting PAUSE_P to TRUE until a new mailbox overrides.
    pause_loop: do

       file_p = .FALSE.
       IF ( ionode ) INQUIRE( FILE = TRIM( mbfile ), EXIST = file_p )
       call mp_bcast(file_p, ionode_id,world_comm)

       IF ( file_p ) THEN

          IF ( ionode ) THEN
            WRITE(*,*)
            WRITE(*,*) '****************************************************'
            WRITE(*,*) '  Autopilot: Mailbox found at nfi=', current_nfi
          END IF
          FLUSH(6)

          ! Open the mailbox
          IF ( ionode ) OPEN( UNIT = pilot_unit, FILE = TRIM( mbfile ) )

          ! Will reset PAUSE_P to false unless there is a PAUSE cmd
          ! The following call is MPI safe! It only generates side effects
          CALL parse_mailbox()
          call mp_barrier( world_comm )
          
          IF ( ionode ) THEN
            WRITE(*,*) '  Autopilot: Done reading mailbox' 
            WRITE(*,*) '****************************************************'
            WRITE(*,*)
          END IF

          ! Perhaps instead of deleting move the file as an input log     
          IF( ionode ) CLOSE( UNIT = pilot_unit, STATUS = 'DELETE' )

       END IF

       IF( .NOT. pause_p ) THEN
          EXIT pause_loop
       ELSE
          IF( ionode ) write(*,*) 'SLEEPING .... send another pilot.mb'
          call sleep (5)
       END if

    end do pause_loop




    ! Autopilot (Dynamic Rules) Implementation
    ! When nfi has passed (is greater than
    ! the next event, then employ rules
    ! Mailbox may have issued several rules
    ! Attempt to catch up! 
    do while (current_nfi >= event_step(event_index) )         

       IF ( ionode ) THEN
          WRITE(*,*)
          WRITE(*,*) '****************************************************'
          WRITE(*,*) '  Autopilot employ rules: '
       END IF
       call employ_rules()
       IF ( ionode ) THEN
          WRITE(*,*) '****************************************************'
          WRITE(*,*)
       END IF
       call mp_barrier( world_comm )

       ! update event_index to current
       event_index = event_index + 1

    enddo
    
    ! During the last step of electron conjugate gradient,
    ! set tprint to .true. if next step has electron_dynamics = Verlet
    ! This is needed to calculate wavefunctions at time t and (t - dt)
    if (need_tprint_true() ) then
       tprint = .true.
       if (ionode) then
            WRITE(*,*) '=================================================='
            WRITE(*,*) ' Setting tprint=.true. for this step (last of CG)'
            WRITE(*,*) '=================================================='
       endif
    endif

  end subroutine pilot

  function  need_tprint_true()
  ! Check if next step is a Verlet. In such case returns .true.
  ! Used whenever tprint == .true. is needed, e.g. to let CG
  ! calculate wavefunctions at time t and (t - dt) via projections
  ! onto the occupied manifold.
  ! This function has to be called after 'call employ_rules()'!
      USE cg_module, ONLY : tcg
      LOGICAL :: need_tprint_true
      INTEGER :: event_idx

      need_tprint_true = .FALSE.
      event_idx = event_index
      do while ( event_idx .le. size(event_step) .and. (event_step(event_idx)==current_nfi+1) ) 

          if ( tcg  .and. event_electron_dynamics(event_idx) .and. &
            &  (rule_electron_dynamics(event_idx)=='VERLET') ) then
            
               need_tprint_true = .TRUE.

          end if 
	  event_idx = event_idx + 1
      enddo
      RETURN
  end function need_tprint_true

END MODULE cp_autopilot

