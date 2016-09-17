!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
PROGRAM pp
  !-----------------------------------------------------------------------
  !
  !    Program for data analysis and plotting. The two basic steps are:
  !    1) read the output file produced by pw.x, extract and calculate
  !       the desired quantity (rho, V, ...)
  !    2) write the desired quantity to file in a suitable format for
  !       various types of plotting and various plotting programs
  !    The two steps can be performed independently. Intermediate data
  !    can be saved to file in step 1 and read from file in step 2.
  !
  !    DESCRIPTION of the INPUT : see file Doc/INPUT_PP.*
  !
  USE io_global,  ONLY : ionode
  USE mp_global,  ONLY : mp_startup
  USE environment,ONLY : environment_start, environment_end

  !
  IMPLICIT NONE
  !
  CHARACTER(len=256) :: filplot
  INTEGER :: plot_num
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'POST-PROC' )
  !
  IF ( ionode )  CALL input_from_file ( )
  !
  CALL extract (filplot, plot_num)
  !
  CALL chdens (filplot, plot_num)
  !
  CALL environment_end ( 'POST-PROC' )
  !
  CALL stop_pp()
  !
END PROGRAM pp
!
!-----------------------------------------------------------------------
SUBROUTINE extract (filplot,plot_num)
  !-----------------------------------------------------------------------
  !
  !    This subroutine reads the data for the output file produced by pw.x
  !    extracts and calculates the desired quantity (rho, V, ...)
  !    writes it to a file for further processing or plotting
  !
  !    DESCRIPTION of the INPUT: see file Doc/INPUT_PP
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : bg
  USE ener,      ONLY : ef
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp, tau
  USE gvect
  USE fft_base,  ONLY : dfftp
  USE klist,     ONLY : two_fermi_energies
  USE vlocal,    ONLY : strf
  USE io_files,  ONLY : tmp_dir, prefix
  USE io_global, ONLY : ionode, ionode_id
  USE mp_global,     ONLY : nproc_pool, nproc_file, nproc_pool_file
  USE control_flags, ONLY : twfcollect
  USE noncollin_module, ONLY : i_cons
  USE paw_variables, ONLY : okpaw
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  USE constants, ONLY : rytoev

  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER(len=256), INTENT(out) :: filplot
  INTEGER, INTENT(out) :: plot_num

  INTEGER :: kpoint, kband, spin_component, ios
  LOGICAL :: lsign, needwf

  REAL(DP) :: emin, emax, sample_bias, z, dz, epsilon
  ! directory for temporary files
  CHARACTER(len=256) :: outdir

  NAMELIST / inputpp / outdir, prefix, plot_num, sample_bias, &
       spin_component, z, dz, emin, emax, kpoint, kband, &
       filplot, lsign, epsilon

  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filplot = 'tmp.pp'
  plot_num = -1
  spin_component = 0
  sample_bias = 0.01d0
  z = 1.d0
  dz = 0.05d0
  lsign=.false.
  emin = -999.0d0
  emax = +999.0d0
  epsilon=1.d0
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     !     reading the namelist inputpp
     !
     READ (5, inputpp, iostat = ios)
     !
     tmp_dir = trimcheck ( outdir )
     !
  ENDIF
  !
  CALL mp_bcast (ios, ionode_id, world_comm)
  !
  IF ( ios /= 0) CALL errore ('postproc', 'reading inputpp namelist', abs(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( plot_num, ionode_id, world_comm )
  CALL mp_bcast( sample_bias, ionode_id, world_comm )
  CALL mp_bcast( spin_component, ionode_id, world_comm )
  CALL mp_bcast( z, ionode_id, world_comm )
  CALL mp_bcast( dz, ionode_id, world_comm )
  CALL mp_bcast( emin, ionode_id, world_comm )
  CALL mp_bcast( emax, ionode_id, world_comm )
  CALL mp_bcast( kband, ionode_id, world_comm )
  CALL mp_bcast( kpoint, ionode_id, world_comm )
  CALL mp_bcast( filplot, ionode_id, world_comm )
  CALL mp_bcast( lsign, ionode_id, world_comm )
  CALL mp_bcast( epsilon, ionode_id, world_comm )
  !
  ! no task specified: do nothing and return
  !
  IF (plot_num == -1) RETURN
  !
  IF (plot_num < 0 .or. plot_num > 21) CALL errore ('postproc', &
          'Wrong plot_num', abs (plot_num) )

  IF (plot_num == 7 .or. plot_num == 13 .or. plot_num==18) THEN
     IF  (spin_component < 0 .or. spin_component > 3) CALL errore &
          ('postproc', 'wrong spin_component', 1)
  ELSEIF (plot_num == 10) THEN
     IF  (spin_component < 0 .or. spin_component > 2) CALL errore &
          ('postproc', 'wrong spin_component', 2)
  ELSE
     IF (spin_component < 0 ) CALL errore &
         ('postproc', 'wrong spin_component', 3)
  ENDIF
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  needwf=(plot_num==3).or.(plot_num==4).or.(plot_num==5).or.(plot_num==7).or. &
         (plot_num==8).or.(plot_num==10)
  IF ( needwf ) THEN
     CALL read_file ( )
     IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect)  &
        CALL errore('postproc', &
        'pw.x run with a different number of procs/pools. Use wf_collect=.true.',1)
     CALL openfil_pp ( )
  ELSE
     CALL read_xml_file ( )
  END IF
  !
  IF ( ( two_fermi_energies .or. i_cons /= 0) .and. &
       ( plot_num==3 .or. plot_num==4 .or. plot_num==5 ) ) &
     CALL errore('postproc',&
     'Post-processing with constrained magnetization is not available yet',1)
  !
  ! The following line sets emax to its default value if not set
  ! It is done here because Ef must be read from file
  !
  IF (emax == +999.0d0) emax = ef
  IF (plot_num == 10) THEN
     emin = emin / rytoev
     emax = emax / rytoev
  ENDIF
  !
  !   Now do whatever you want
  !
  CALL punch_plot (filplot, plot_num, sample_bias, z, dz, &
       emin, emax, kpoint, kband, spin_component, lsign, epsilon)
  !
END SUBROUTINE extract
