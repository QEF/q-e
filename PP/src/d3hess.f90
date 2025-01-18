!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
program d3hess
  USE io_global,        ONLY: ionode, ionode_id
  USE io_files,         ONLY: prefix, tmp_dir
  USE kinds,            ONLY: DP
  USE mp,               ONLY: mp_bcast
  USE mp_global,        ONLY: mp_startup
  USE mp_world,         ONLY: world_comm
  USE environment,      ONLY: environment_start, environment_end
  USE d3hess_mod,       ONLY: q_gamma, debug, step, d3hess_sub
  !
  IMPLICIT NONE
  INTEGER :: ios
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(len=256) :: filhess, outdir
  LOGICAL :: needwf = .FALSE.
  !
  NAMELIST /input/ prefix, outdir, step, q_gamma, filhess, debug
  !
  ! Initialize environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'd3hess' )
  !
  ios = 0
  !
  IF ( ionode ) THEN
     !
     !   set default values for variables in namelist
     !
     CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
     IF ( trim( outdir ) == ' ' ) outdir = './'
     prefix ='pwscf'
     filhess=' ' 
     q_gamma=.false. ! whether to use a much cheaper algorithm when q=0,0,0
     debug=.false.   ! whether to check consistency between hessian, forces and energies
     step=2.d-5      ! step for numerical differentiation
     !
     CALL input_from_file ( )
     READ (5,input,IOSTAT=ios)
     !
     tmp_dir = trimcheck (outdir)
     IF ( filhess == ' ' ) filhess = trim(prefix)//'.hess'
     filhess = TRIM(tmp_dir)//TRIM(filhess)
     !
  ENDIF
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast(ios, ionode_id, world_comm)
  CALL errore('d3hess', 'reading input namelist', ABS(ios))

  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( step, ionode_id, world_comm )
  CALL mp_bcast( q_gamma, ionode_id, world_comm )
  CALL mp_bcast( debug, ionode_id, world_comm )
  !
  CALL read_file_new ( needwf )
  !
  ! DFT-D3 functional dependent parameters have been set in read_file_new
  !
  CALL d3hess_sub(filhess)
  !
  CALL environment_end ( 'd3hess' )
  !
  CALL stop_pp
  !
end program
