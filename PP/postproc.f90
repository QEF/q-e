!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
PROGRAM postproc
  !-----------------------------------------------------------------------
  !
  !    This program reads the output file produced by pw.x
  !    extracts and calculates the desired quantity (rho, V, ...)
  !    writes it to a file for further processing or plotting
  !
  !    DESCRIPTION of the INPUT: see file pwdocs/INPUT_PP
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : bg
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp, tau
  USE gvect
  USE ener,      ONLY : ef
  USE vlocal,    ONLY : strf
  USE io_files,  ONLY : tmp_dir, nd_nmbr, prefix
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_bcast

  IMPLICIT NONE
  CHARACTER(len=256) :: filplot

  INTEGER :: plot_num, kpoint, kband, spin_component, ios
  LOGICAL :: stm_wfc_matching, lsign

  REAL(kind=DP) :: emin, emax, sample_bias, z, dz
  ! directory for temporary files
  CHARACTER(len=256) :: outdir

  NAMELIST / inputpp / outdir, prefix, plot_num, stm_wfc_matching, &
       sample_bias, spin_component, z, dz, emin, emax, kpoint, kband,&
       filplot, lsign

  CHARACTER (LEN=256) :: input_file
  INTEGER             :: nargs, iiarg, ierr, ILEN
  INTEGER, EXTERNAL   :: iargc

  !
  CALL start_postproc (nd_nmbr)
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  outdir = './'
  filplot = 'pp.out' 
  plot_num = 0
  spin_component = 0
  sample_bias = 0.01d0
  z = 1.d0
  dz = 0.05d0
  stm_wfc_matching = .TRUE.
  lsign=.FALSE.
  emin = - 999.0d0
  emax = ef*13.6058d0

  IF ( ionode )  THEN
     !
     ! ... Input from file ?
     !
     nargs = iargc()
     !
     DO iiarg = 1, ( nargs - 1 )
        !
        CALL getarg( iiarg, input_file )
        IF ( TRIM( input_file ) == '-input' .OR. &
             TRIM( input_file ) == '-inp'   .OR. &
             TRIM( input_file ) == '-in' ) THEN
           !
           CALL getarg( ( iiarg + 1 ) , input_file )
           OPEN ( UNIT = 5, FILE = input_file, FORM = 'FORMATTED', &
                STATUS = 'OLD', IOSTAT = ierr )
           CALL errore( 'iosys', 'input file ' // TRIM( input_file ) // &
                & ' not found' , ierr )
           !
        END IF
        !
     END DO
     !
     !     reading the namelist inputpp
     !
     READ (5, inputpp, err = 200, iostat = ios)
200  CALL errore ('postproc', 'reading inputpp namelist', ABS (ios) )
     tmp_dir = TRIM(outdir)
     !
  END IF
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( plot_num, ionode_id )
  CALL mp_bcast( stm_wfc_matching, ionode_id )
  CALL mp_bcast( sample_bias, ionode_id )
  CALL mp_bcast( spin_component, ionode_id )
  CALL mp_bcast( z, ionode_id )
  CALL mp_bcast( dz, ionode_id )
  CALL mp_bcast( emin, ionode_id )
  CALL mp_bcast( emax, ionode_id )
  CALL mp_bcast( kband, ionode_id )
  CALL mp_bcast( kpoint, ionode_id )
  CALL mp_bcast( filplot, ionode_id )
  CALL mp_bcast( lsign, ionode_id )
  !
  !     Check of namelist variables
  !
  IF (plot_num < 0 .OR. plot_num > 13) CALL errore ('postproc', &
          'Wrong plot_num', ABS (plot_num) )

  IF ( (plot_num == 0 .OR. plot_num == 1) .AND.  &
       (spin_component < 0 .OR. spin_component > 2) ) CALL errore &
         ('postproc', 'wrong value of spin_component', 1)

  IF ( (plot_num == 13) .AND.   &
       (spin_component < 0 .OR. spin_component > 3) ) CALL errore &
          ('postproc', 'wrong spin_component', 1)


  IF (plot_num == 10) THEN
     emin = emin / 13.6058d0
     emax = emax / 13.6058d0
  END IF

  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file
  CALL openfil_pp
  CALL struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, nr3, &
       strf, eigts1, eigts2, eigts3)
  CALL init_us_1
  !
  !   Now do whatever you want
  !
  CALL punch_plot (filplot, plot_num, sample_bias, z, dz, &
       stm_wfc_matching, emin, emax, kpoint, kband, spin_component, lsign)
  !
  CALL stop_pp
  STOP
END PROGRAM postproc
