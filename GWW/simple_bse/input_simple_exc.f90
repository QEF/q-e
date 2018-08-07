MODULE input_simple_exc
!this module provides input file routines 

  USE kinds, ONLY: DP

   TYPE input_options
      CHARACTER(len=256) :: prefix = 'prefix'!prefix to designate the files same as in PW        
      CHARACTER(len=256) :: outdir = './'!outdir to designate the files same as in PW    
      INTEGER :: task !0 find eigenvectors/value , 1 calculate spectrum
      INTEGER :: diago!for task==0, diagonalization scheme: 0=steepest descent, 1=conjugate gradient
      INTEGER :: nvec!for task=0 number of eigen vectors to be found
      REAL(kind=DP) :: thr_evc!threshold for the eigenvectors/value
      INTEGER :: max_nstep!maximum number of steps during minimization
      REAL(kind=DP) :: l_step!length of trial steps for minimization
      INTEGER :: h_level!level of excitonic hamiltonian: 0 RPA (diagonal), 1 TD-H (+exchange), 2 TD-HF(+direct-bare), 3-BSE (+direct-correlation)
      INTEGER :: spin_state!only for non spin polarised systems: 0:triplet, 1:singlet, 2:full (spin-orbit)
      INTEGER :: spectrum_points!number of points in the spectrum
      REAL(KIND=DP) :: omega_min!max omega in the spectrum
      REAL(KIND=DP) :: omega_max!min omega in the spectrum
      REAL(KIND=DP) :: eta!infinitesimal
      INTEGER :: lanczos_step!!number of steps for haydoch recursive method
      REAL(kind=DP) :: scissor!Gap scissor (eV)
   END TYPE input_options

   CONTAINS

     SUBROUTINE  read_input_simple_exc( simple_in )
       USE io_global,            ONLY : stdout, ionode, ionode_id
       USE mp,                   ONLY : mp_bcast
       USE mp_world,             ONLY : world_comm
       USE io_files,             ONLY :  tmp_dir,prefix
       

       implicit none

       CHARACTER(LEN=256), EXTERNAL :: trimcheck
       TYPE(input_options) :: simple_in  !in output the input parameters

       NAMELIST/inputsimple/simple_in

	CHARACTER(LEN=256) :: outdir

       if(ionode) then
          read(*,NML=inputsimple)
          outdir = trimcheck(simple_in%outdir)
          tmp_dir = outdir
          prefix = trim(simple_in%prefix)
       endif

       CALL mp_bcast( outdir,ionode_id, world_comm )
       CALL mp_bcast( tmp_dir,ionode_id, world_comm )
       CALL mp_bcast( prefix,ionode_id, world_comm )
       
       CALL mp_bcast( simple_in%task,ionode_id, world_comm )
       CALL mp_bcast( simple_in%diago,ionode_id, world_comm )
       CALL mp_bcast( simple_in%nvec,ionode_id, world_comm )
       CALL mp_bcast( simple_in%thr_evc,ionode_id, world_comm )
       CALL mp_bcast( simple_in%max_nstep,ionode_id, world_comm )
       CALL mp_bcast( simple_in%l_step,ionode_id, world_comm )
       CALL mp_bcast( simple_in%h_level,ionode_id, world_comm )
       CALL mp_bcast( simple_in%spin_state,ionode_id, world_comm )
       CALL mp_bcast( simple_in%spectrum_points,ionode_id, world_comm )
       CALL mp_bcast( simple_in%omega_max,ionode_id, world_comm )
       CALL mp_bcast( simple_in%omega_min,ionode_id, world_comm )
       CALL mp_bcast( simple_in%eta,ionode_id, world_comm )
       CALL mp_bcast( simple_in%lanczos_step,ionode_id, world_comm )
       CALL mp_bcast( simple_in%scissor,ionode_id, world_comm )      
 
       return

     END SUBROUTINE read_input_simple_exc

END MODULE input_simple_exc
