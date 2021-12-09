    Copyright (c) Targacept, Inc.

    Targacept, Inc., 
    200 East First Street, 
    Suite 300, 
    Winston-Salem, NC, USA 27101 
    atp@targacept.com

This file describes the Autopilot Feature Suite as introduced and used by 
Targacept, Inc.   This documentation accompanies free software; The software
is subject to the terms of the GNU General Public License as published by the 
Free Software Foundation; either version 2 of the License, or (at your option) 
any later version. See the GNU General Public License at 
www.gnu.or/copyleft/gpl.txt for more details.

This documentation, like the software it accompanies, is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
warranty of MERCHANTABILITY FOR A PARTICULAR PURPOSE.  

[[_TOC_]]

--------------------------------------------------------------------------------
AUTOPILOT DOCUMENTATION
--------------------------------------------------------------------------------

The Autopilot Feature Suite is a user level enhancement for directing 
Car-Parrinello simulations based on CP.X packaged in ESPRESSO. 

The following features are incorporated: 

 - Auto Restart Mode 
 - Autopilot Course Configuration (Dynamic Rules) 
 - Autopilot Course Correction (Steering) 


--------------------------------------------------------------------------------
Auto Restart Mode 
--------------------------------------------------------------------------------

Auto Restart Mode is an extension of restart_mode declared in the CONTROL section 
of the input file. When restart mode is set to "auto", control determines if the 
current run is "from_scratch" or a valid "restart" based on the presence of a 
restart file associated with unit NDR. When NDR, the unit number for input, and 
NDW, the unit number for output, are the same, a simulation that is system 
terminated can be restarted without significant loss, providing that ISAVE, the 
parameter that indicates the frequency at which intermediate data are saved, is 
not large. 

Auto Restart Mode implements an effective "upto" mode and is also designed for 
use on remote machines where simulations may frequently be terminated and 
restarted.  Auto Restart Mode is especially useful in connection with Autopilot's 
Dynamic Rules capability. When they are used together, only one segment of a 
simulation is necessary, thereby reducing run_script volume and errors, and 
placing more control with the user. 

    restart_mode   CHARACTER ( default = 'restart' )
           from_scratch = from scratch.  NEB only: the starting path is 
                             obtained with a linear interpolation between 
                             the images specified in the ATOMIC_POSITIONS 
                             card.  Note that in the linear interpolation,
                             periodic boundary conditions ARE NOT USED.
           restart  = continue a previous simulation and perform  
                             "nstep" new steps.
           reset_counters  = continue a previous simulation, perform  
                             "nstep" new steps, resetting the counter 
                             and averages.
           auto = automatically detect "from_scratch" or "restart"; 
                             continue any previous simulation, and stop 
                             when the counter value is equal to "nstep".




--------------------------------------------------------------------------------
Autopilot Course Configuration (Dynamic Rules) 
--------------------------------------------------------------------------------

Autopilot Course Configuration (Dynamic Rules) is a method that allows select 
input parameters (Autopilot variables) to change during the course of a 
simulation. This method allows the user to create a more concise set of 
instructions that are easier to read and maintain and enables a more continuous 
execution on remote resources. 

Typically and historically, a user issues a run_script that creates a sequence of
input files, each with fixed parameter values. This run_script then calls cp.x 
against each input file in the sequence, such that, after the first, each 
execution continues with the next input file as well as restart information from
the previous execution. 

The Autopilot Course Configuration effectively consolidates multiple input files 
into one, allowing the user to specify at what time step a parameter should change
along with its new value. Thus a run_script becomes much shorter, and the user can 
easily see the projected path of the simulation. 

The Autopilot Course Configuration feature is implemented by adding a new card type 
to the "CARDS" section of the input file. The Autopilot card must be placed after 
the "NAMELIST" section but otherwise may appear before or after any other card. 
A favorable place is as the first card. 

Sytnax is as follows: 

    CARDS ...

    AUTOPILOT

      optional card :  read dynamic rules to set parameters on an absolute
                       timestep (iteration) from either standard input or mailbox 
                       (pilot.mb)
      Syntax:

    AUTOPILOT
      ON_STEP = ith_event_STEP : varname = value  
      ON_STEP = jth_event_STEP : varname = value  
      ON_STEP = jth_event_STEP : varname = value  

    ...
      ON_STEP = nth_event_STEP : varname = value  
      ENDRULES

Description:

 - `ON_STEP` 	   LABEL, must be in numerical timestep order, otherwise rule 
    is ignored

 - `ith_event_STEP`   INTEGER, iteration (NFI) when rule is to be employed

 - `varname`	   Autopilot variable, currently limited to one of the 
    following: isave,iprint,dt,emass, electron_dynamics, 
    electron_damping, ion_dynamics, ion_damping, 
    ion_temperature, tempw.

 - `value`            Must be valid value of variable type 
    (for example: isave, iprint must have a value of type 
    INTEGER, while dt must have a value of type REAL)

 - `ENDRULES`         Required only for input (STDIN) if other cards follow.

The event specification (`ON_STEP`) should precede the variable assignment. The 
colon separator between the event assignment and the variable assignment is 
required, as are the equal signs. No semi-colon or comma should appear after the
variable assignment. There can be multiple rules per event but only one variable 
assignment per rule (and only one rule per line). Within one event, there should 
be only one assignment per variable.  If multiple assignments are made for the 
same variable for the same event, only the last assignment will be accepted. 
Rules for which event specifications are not in numerical order will be ignored. 
If syntax errors are present in the AUTOPILOT card during start-up, the 
execution will stop. 

Example Syntax: 

      AUTOPILOT
        ON_STEP = 200 : tempw = 500.0
        ON_STEP = 200 : dt = 3.0
        ON_STEP = 250 : ISAVE = 50
      ENDRULES

Currently there is a maximum of 32 supported events and 10 supported Autopilot 
variables. Events that are out of timestep order are ignored. A user may establish 
up to 10 rules (one for each Autopilot variable) per event. Currently implemented 
Autopilot variables are: isave, iprint, dt, emass, electron_dynamics, 
electron_damping, ion_dynamics, ion_damping, ion_temperature, and tempw. 
If desired, users may implement other Autopilot variables. See Appendix below for 
an explanation of "Adding an Autopilot Variable". 

IMPORTANT: Variables should have values in accordance with their TYPE, or a 
runtime error may occur.





--------------------------------------------------------------------------------
Autopilot Course Correction (Steering) 
--------------------------------------------------------------------------------

Autopilot Course Correction (Steering) provides a run-time method of changing 
Autopilot variables on the fly, after the simulation is underway. Autopilot 
Course Correction (Steering) can be applied through any of the following 
sub-features: New Course (power steering), Manual Steering, and Pause. 

Steering utilizes a new mailbox file:  pilot.mb. This file can be created via the
user's favorite text editor and can be "mailed" by placing the file in the 
"results" directory. The user can also very quickly implement a single course 
correction command with UNIX redirect to the pilot.mb file. 

When a pilot.mb mailbox file is detected, the current event table is cleared to 
prepare for the new course.  The mailbox file is then parsed, and Autopilot 
processes the command(s) before deleting the mailbox file. If Autopilot cannot 
parse a command, it issues a warning and goes into PAUSE mode (see below). 

The Steering subfeatures, including pilot.mb syntax are described here: 

 - New Course or 'power steering' is implemented with the same syntax as the 
   INPUT file card for Autopilot. Remember that ON_STEP represents an absolute 
   iteration (NFI) step. 

   For example: 

         AUTOPILOT                               -required
         ON_STEP=400 : ISAVE = 50                -events must be ordered by step
         ON_STEP=400 : DT = 5.0                  -use valid variable types (or die)
         ON_STEP = 600:IONS_TEMPERATURE='damped' -indention optional     
         ON_STEP = 600: TEMPW=350.0              -white spaces are ignored
         ENDRULES                                -optional

   In this example, when NFI reaches 400, the value of ISAVE will be reset to 50 
   and the value of DT to 5.0.  Then, when NFI reaches 600, IONS_TEMPERATURE and 
   TEMPW will be reset to the indicated values.

 - Manual Steering is implemented with a similar syntax except that the card type 
   is PILOT instead of AUTOPILOT and the user specifies a timestep relative to 
   the time the mailbox is read, rather than an absolute timestep. The relative 
   timestep allows the user to set a rule for a near future event without having 
   to judge the current absolute NFI value. The user may also pre-write multiple 
   mailboxes using relative event steps without regard to absolute iteration (NFI) 
   values.

   For example, assume mailbox contents are:

         NOW:ISAVE=50
         NOW+100:TEMPW=600.0

   Assume further that the mailbox is saved to the "results" directory and then 
   read when the NFI is 380.  Manual Steering will reset the value of ISAVE on 
   the next event that is modulo 50, and an ISAVE event will occur twice 
   (at 400 and again at 450) before TEMPW is reset to 600.0 on step 480. Compare 
   this with the syntax that specifies an absolute timestep:

         ON_STEP=400:ISAVE=50
         ON_STEP=500;TEMPW=600.0


   In this example, if the NFI is less than 400 when the mailbox is read, ISAVE 
   becomes 50 on step 400 and TEMPW becomes 600.0 on step 500, and ISAVE is 
   performed twice before TEMPW is reset, just as in the previous example that 
   uses relative indexing.  

   However, if the user misjudges the momentary NFI, and it is 530 when the 
   mailbox is read, then both rules are implemented immediately and simultaneously.
   Furthermore, the ISAVE rule takes effect after the NFI specified. Neither of 
   these effects may have been intended by the user.

   Following is an example of a Manual Steering mailbox to change temperature from 
   a relative iteration (NFI) step: 

   Example syntax for a Manual Steering mailbox is as follows: 

         PILOT                                -optional for single line
           NOW : ISAVE = 50                   -events must be ordered
           NOW : DT = 5.0                     -use valid variable types (or die)
           NOW+50 :IONS_TEMPERATURE='damped'  -offsets from NOW are supported  
           NOW + 150: TEMPW=350.0             -white spaces are ignored
         ENDRULES                             -optional

  Example format for a quick mailbox change using a single rule is as follows:

                                      	  -defaults to PILOT
     NOW + 250: TEMPW=450.0               -single line with NOW 

c) Pause is a steering sub-feature that allows the user to suspend the simulation
until the user can decide on a future course. Pause is very helpful when the user 
knows that a change should be imposed but needs time to establish rules and 
create an appropriate mailbox. Steering then resumes as AUTOPILOT or PILOT upon 
receiving another pilot.mb mailbox. The syntax is a single line with one of the 
following: 

      PAUSE
      SLEEP
      HOLD
      HOVER
      WAIT

All of the above perform the same PAUSE mechanism.  The user can issue the command
quickly through UNIX redirect:

  >echo "PAUSE" > results/pilot.mb

Any mailbox not correctly identified with a `AUTOPILOT`, `PILOT`, `NOW`, or a `PAUSE` 
command, will result in a warning to standard output (STDOUT), and the simulation 
will pause. 





--------------------------------------------------------------------------------
TESTING 
--------------------------------------------------------------------------------

The entire Autopilot Feature Suite issues directives to slave nodes under 
MPI, with the fewest broadcasted parameters. All features have been tested 
under Intel 8.1 with MKL 7.0.1 libraries on a Linux-32 single processor and under 
PGI 5.2 with MPI on Linux-64 with 1, 2 and 4 processors. 





--------------------------------------------------------------------------------
ADDING AN AUTOPILOT VARIABLE
--------------------------------------------------------------------------------

See `autopilot.f90` for examples.
- Select the input parameter from the list in file INPUT_CP 
- Identify parameter dependencies, initializations, assignments, etc 
- Edit autopilot.f90 to add the following, 
    where VARNAME is the name of the new Autopilot variable:
        VARTYPE :: rule_VARNAME(max_event_step) at module scope 
        LOGICAL :: event_VARNAME(max_event_step) at module scope
* Remember to add to the PUBLIC block as well
        event_VARNAME(:) = .false. to init_autopilot subroutine 
        rule_VARNAME(:) = VARDEFAULT to init_autopilot subroutine 
* Import VARNAME with USE to employ_rules subroutine
* In employ_rules, add conditional clause on event_VARNAME to assign VARNAME: 
         ! VARNAME
         if (event_VARNAME(event_index)) then
           VARNAME  = rule_VARNAME(event_index)
           CALL init_other_VARNAME_dependent_variables( VARNAME)
           write(*,*) 'RULE EVENT: VARNAME', VARNAME
         endif
* Import VARNAME with USE to assign_rule subroutine
* In assign_rule, add condition clause matching the VARNAME create rule as so: 
         ELSEIF ( matches( "VARNAME", var ) ) THEN
                     read(value, *) VARTYPE_value
                     rule_VARNAME(event)  = VARTYPE_value
                     event_VARNAME(event) = .true.
* TEST  

WARNING: Some Autopilot variables may create "side-effects".  For example, the 
inclusion of a rule for TEMPW rules invokes a side-effect call to ions_nose_init.  
The user is cautioned to be aware of possible side-effects when adding other 
Autopilot variables. 

