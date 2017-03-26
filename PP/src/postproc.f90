!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE pp_module
CONTAINS
!-----------------------------------------------------------------------
SUBROUTINE extract (plot_files,plot_num)
  !-----------------------------------------------------------------------
  !
  !    Reads data produced by pw.x, computes the desired quantity (rho, V, ...)
  !    and writes it to a file (or multiple files) for further processing or
  !    plotting
  !
  !    On return, plot_files contains a list of all written files.
  !
  !    DESCRIPTION of the INPUT: see file Doc/INPUT_PP
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : bg
  USE ener,      ONLY : ef
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp, tau
  USE gvect
  USE fft_base,  ONLY : dfftp
  USE klist,     ONLY : two_fermi_energies, degauss
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
  USE parameters, ONLY : npk
  USE io_global, ONLY : stdout

  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER(len=256), DIMENSION(:), ALLOCATABLE, INTENT(out) :: plot_files
  INTEGER, INTENT(out) :: plot_num

  CHARACTER (len=2), DIMENSION(0:3) :: spin_desc = &
       (/ '  ', '_X', '_Y', '_Z' /)

  INTEGER :: kpoint(2), kband(2), spin_component(3), ios
  LOGICAL :: lsign, needwf

  REAL(DP) :: emin, emax, sample_bias, z, dz
  
  REAL(DP) :: degauss_ldos, delta_e
  CHARACTER(len=256) :: filplot
  INTEGER :: plot_nkpt, plot_nbnd, plot_nspin, nplots
  INTEGER :: iplot, ikpt, ibnd, ispin

  ! directory for temporary files
  CHARACTER(len=256) :: outdir

  NAMELIST / inputpp / outdir, prefix, plot_num, sample_bias, &
      spin_component, z, dz, emin, emax, delta_e, degauss_ldos, kpoint, kband, &
      filplot, lsign
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filplot = 'tmp.pp'
  plot_num = -1
  kpoint(2) = 0
  kband(2) = 0
  spin_component = 0
  sample_bias = 0.01d0
  z = 1.d0
  dz = 0.05d0
  lsign=.false.
  emin = -999.0d0
  emax = +999.0d0
  delta_e=0.1d0
  degauss_ldos=-999.0d0
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
  CALL mp_bcast( degauss_ldos, ionode_id, world_comm )
  CALL mp_bcast( delta_e, ionode_id, world_comm )
  CALL mp_bcast( kband, ionode_id, world_comm )
  CALL mp_bcast( kpoint, ionode_id, world_comm )
  CALL mp_bcast( filplot, ionode_id, world_comm )
  CALL mp_bcast( lsign, ionode_id, world_comm )
  !
  ! no task specified: do nothing and return
  !
  IF (plot_num == -1) THEN
     ALLOCATE( plot_files(0) )
     RETURN
  ENDIF
  !
  IF (plot_num < 0 .or. plot_num > 22) CALL errore ('postproc', &
          'Wrong plot_num', abs (plot_num) )

  IF (plot_num == 7 .or. plot_num == 13 .or. plot_num==18) THEN
     IF  (spin_component(1) < 0 .or. spin_component(1) > 3) CALL errore &
          ('postproc', 'wrong spin_component', 1)
  ELSEIF (plot_num == 10) THEN
     IF  (spin_component(1) < 0 .or. spin_component(1) > 2) CALL errore &
          ('postproc', 'wrong spin_component', 2)
  ELSE
     IF (spin_component(1) < 0 ) CALL errore &
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
  ! Set default values for emin, emax, degauss_ldos
  ! Done here because ef, degauss must be read from file
  IF (emin > emax) CALL errore('postproc','emin > emax',0)
  IF (plot_num == 10) THEN
      IF (emax == +999.0d0) emax = ef * rytoev
  ELSEIF (plot_num == 3) THEN
      IF (emin == -999.0d0) emin = ef * rytoev
      IF (emax == +999.0d0) emax = ef * rytoev
      IF (degauss_ldos == -999.0d0) THEN
          WRITE(stdout, &
              '(/5x,"degauss_ldos not set, defaults to degauss = ",f6.4, " eV")') &
             degauss * rytoev
          degauss_ldos = degauss * rytoev
      ENDIF
  ENDIF
  ! transforming all back to Ry units
  emin = emin / rytoev
  emax = emax / rytoev
  delta_e = delta_e / rytoev
  degauss_ldos = degauss_ldos / rytoev

  ! Number of output files depends on input
  nplots = 1
  IF (plot_num == 3) THEN
     nplots=(emax-emin)/delta_e + 1
  ELSEIF (plot_num == 7) THEN
      IF (kpoint(2) == 0)  kpoint(2) = kpoint(1)
      plot_nkpt = kpoint(2) - kpoint(1) + 1
      IF (kband(2) == 0)  kband(2) = kband(1)
      plot_nbnd = kband(2) - kband(1) + 1
      IF (spin_component(2) == 0)  spin_component(2) = spin_component(1)
      plot_nspin = spin_component(2) - spin_component(1) + 1

      nplots = plot_nbnd * plot_nkpt * plot_nspin
  ENDIF
  ALLOCATE( plot_files(nplots) )
  plot_files(1) = filplot

  ! 
  ! First handle plot_nums with multiple calls to punch_plot
  !
  IF (nplots > 1 .AND. plot_num == 3) THEN
  ! Local density of states on energy grid of spacing delta_e within [emin, emax]
    DO iplot=1,nplots
      WRITE(plot_files(iplot),'(A, I0.3)') TRIM(filplot), iplot
      CALL punch_plot (TRIM(plot_files(iplot)), plot_num, sample_bias, z, dz, &
        emin, degauss_ldos, kpoint, kband, spin_component, lsign)
      emin=emin+delta_e
    ENDDO
  ELSEIF (nplots > 1 .AND. plot_num == 7) THEN
  ! Plot multiple KS orbitals in one go
    iplot = 1
    DO ikpt=kpoint(1), kpoint(2)
      DO ibnd=kband(1), kband(2)
        DO ispin=spin_component(1), spin_component(2)
          WRITE(plot_files(iplot),"(A,A,I0.3,A,I0.3,A)") &
            TRIM(filplot), "_K", ikpt, "_B", ibnd, TRIM(spin_desc(ispin))
          CALL punch_plot (TRIM(plot_files(iplot)), plot_num, sample_bias, z, dz, &
            emin, emax, ikpt, ibnd, ispin, lsign)
          iplot = iplot + 1
        ENDDO
      ENDDO
    ENDDO

  ELSE
  ! Single call to punch_plot
    IF (plot_num == 3) THEN
       CALL punch_plot (filplot, plot_num, sample_bias, z, dz, &
           emin, degauss_ldos, kpoint, kband, spin_component, lsign)
     ELSE
       CALL punch_plot (filplot, plot_num, sample_bias, z, dz, &
          emin, emax, kpoint, kband, spin_component, lsign)
     ENDIF

  ENDIF
  !
END SUBROUTINE extract

END MODULE pp_module
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
  USE chdens_module, ONLY : chdens
  USE pp_module, ONLY : extract

  !
  IMPLICIT NONE
  !
  !CHARACTER(len=256) :: filplot
  CHARACTER(len=256), DIMENSION(:), ALLOCATABLE :: plot_files
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
  CALL extract (plot_files, plot_num)
  !
  CALL chdens (plot_files, plot_num)
  !
  CALL environment_end ( 'POST-PROC' )
  !
  CALL stop_pp()
  !
END PROGRAM pp
