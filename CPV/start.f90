!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-2000
!  Last modified: Fri Feb 11 14:05:28 MET 2000
!  ----------------------------------------------
!  BEGIN manual

    PROGRAM start

!  this is the main routine
!  it calls preliminary subroutines that read input data and do some
!  initialization, then runs (one or more times) the main Car-Parrinello
!  or DIIS-minimization subroutines
!
!  see the input module for information about the input files
!  see subroutine cpmain for information about input/output units
!  ----------------------------------------------
!  this version features:
!  Parrinello-Rahman dynamics
!  generic k-points calculation
!  Nose' thermostat for ions and electrons
!  velocity rescaling for ions
!  Kleinman-Bylander fully non-local pseudopotentials
!  support for local and s, p and d nonlocality
!  generalized gradient corrections
!  core corrections
!  calculus of polarizability
!  DIIS minimization for electrons
!  ions dynamics with DIIS electronic minimization at each step 
!  --------------------------------------------
!  END manual

! ... declare modules
      USE kinds
      USE environment, ONLY: environment_start, environment_end
      USE input_fpmd, ONLY : read_input_file
      USE mp, ONLY: mp_start, mp_end, mp_env
      USE mp_global, ONLY: mp_global_start
      USE io_global, ONLY: io_global_start, io_global_getionode
      USE control_flags, ONLY: lneb, program_name
      USE version

      IMPLICIT NONE      
      
! ... declare variables
      INTEGER :: mpime, nproc, gid, root, ionode_id
      LOGICAL :: ionode
      INTEGER :: iloop

! ... end of declarations
!  ----------------------------------------------

! ... initialize MPI (parallel processing handling)

      root = 0
      CALL mp_start()
      CALL mp_env( nproc, mpime, gid )
      CALL mp_global_start( root, mpime, gid, nproc )

! ... mpime = processor number, starting from 0
! ... nproc = number of processors
! ... gid   = group index
! ... root  = index of the root processor

      program_name = 'FPMD'

! ... initialize input output

      CALL io_global_start( mpime, root )
      CALL io_global_getionode( ionode, ionode_id )

! ... get environment variables (through the C library function getenv)

      cp_version = TRIM (version_number) // " - " // TRIM (version_date)
      CALL environment_start( ionode, mpime, nproc, cp_version )

! ... read in the input file

      CALL read_input_file( lneb )

      IF( lneb ) THEN
        CALL neb_loop( 0 )
      ELSE
        CALL fpmd_loop( 0 )
      END IF

! ... Close the environment

      CALL environment_end( ionode, mpime )

! ... terminate MPI

      CALL mp_end()

      STOP 'FPMD'
    END PROGRAM start
