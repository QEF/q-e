!
! Copyright (C) 2004-2009 Dario Alfe' and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM pw2casino
  !-----------------------------------------------------------------------

  ! This subroutine writes the file "prefix".pwfn.data containing the
  ! plane wave coefficients and other stuff needed by the QMC code CASINO.
  ! May be useful to anybody desiring to extract data from Quantum ESPRESSO
  ! runs, since the output data is quite easy to understand.
  ! If you want to save the Fermi energy and state occupancies as well,
  ! look at tags KN (courtesy of Karoly Nemeth, Argonne)
  ! Not guaranteed to work in parallel execution ! If you want to read data
  ! written by a parallel run, ensure that the data file was saved in a
  ! portable format (option "wf_collect=.true." for PWscf), run pw2casino
  ! serially. Alternatively: run in the same number of processors and pools
  ! of the previous pw.x calculation.
  ! Usage:
  ! * run first a scf calculation with pw.x
  ! * run pw2casino.x with the following input:
  !      &inputpp prefix='...', outdir ='...' /
  !   where prefix and outdir are the same as those used in the scf calculation
  !   (you may use environment variable ESPRESSO_TMPDIR instead of outdir)
  ! * move all your files named prefix.pwfn.data? to pwfn.data?,
  !   merge the pwfn.data? files using the CASINO utility MERGE_PWFN.
  ! * convert to blips running the BLIP utility.
  ! You do not necessarily have to use casino PP's, but you can if you want;
  ! there is a conversion utility in the upftools directory of the espresso
  ! distribution.

  USE io_files,   ONLY : prefix, outdir, tmp_dir, trimcheck
  USE io_global,  ONLY : ionode, ionode_id
  USE mp,         ONLY : mp_bcast
  USE mp_global,  ONLY : mp_startup, npool, nimage, nogrp, npgrp
  USE environment,ONLY : environment_start
  !
  IMPLICIT NONE
  INTEGER :: ios
  LOGICAL :: casino_gather

  NAMELIST / inputpp / prefix, outdir, casino_gather
  !
  ! initialise environment
  !
#ifdef __PARA
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'PW2CASINO' )

  IF ( npool > 1 .or. nimage > 1) THEN
     CALL errore('pw2casino', 'pool or image parallelization not (yet) implemented')
  ENDIF

  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  casino_gather = .false.
  ios = 0
  IF ( ionode )  THEN
     !
     READ (5, inputpp, iostat=ios)
     tmp_dir = trimcheck (outdir)
     !
  ENDIF
  !
  CALL mp_bcast( ios, ionode_id )
  IF ( ios/=0 ) CALL errore('pw2casino', 'reading inputpp namelist', abs(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast(tmp_dir, ionode_id )
  CALL mp_bcast(casino_gather, ionode_id )
  !
  CALL read_file
  CALL openfil_pp
  !
  CALL write_casino_pwfn(casino_gather)
  !
  CALL stop_pp
  STOP

END PROGRAM pw2casino
