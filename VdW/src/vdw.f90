!
! Copyright (C) 2001-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
PROGRAM vdw
  !-----------------------------------------------------------------------
  !
  !    Program for calculating dynamic polarizability of a generic finite system
  !    (e.g., atoms, molecules,...) using TF-vW approximation.
  !    The two basic steps are:
  !    1) read the output file produced by pw.
  !    2) calculate the effective potential
  !    3) solve the modified Sternheimer equation to determine density response,
  !       hence polarizability
  !
  USE kinds,           ONLY : DP
  USE environment,     ONLY : environment_start
  USE mp_global,       ONLY : mp_startup
  USE io_global,       ONLY : ionode, stdout
  USE pwcom
  USE scf,             ONLY : rho
  USE uspp,            ONLY : nkb, vkb
  USE phcom
  USE qpoint,          ONLY : nksq
  USE control_vdw,     ONLY : thresh_veff
  USE check_stop,      ONLY : check_stop_init
  !
  IMPLICIT NONE
  !
  REAL (kind=DP) :: charge, vstart
  INTEGER :: i
  !
  ! initialise environment
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'VdW' )
  !
  IF ( ionode )  CALL input_from_file ( )
  !
  CALL vdw_init ( )
  !
  CALL check_stop_init()
  !
  !!! workaround to prevent array mismatch: vkb is allocated in allocate_nlpot,
  !!! called by read_file, called by vdw_init, with nkb set to the "correct"
  !!! value, but we need nkb=0 and vkb array consistently allocated
  nkb  = 0
  DEALLOCATE (vkb)
  ALLOCATE (vkb(npwx,nkb))
  !!!
  nbnd = 1
  nbndx = 4 * nbnd
  nksq   = 1
  !
  ! Allocate needed variables
  !
  CALL allocate_vdw()
  !
  ! Start polarizability calculation
  !
  WRITE( stdout,'(/,5x,"Effective Potential Calculation",/)')
  !
  ! Calculate the effective potential
  !
  CALL eff_pot (rho%of_r, nspin, alat, omega, &
                charge, vstart, thresh_veff)
  !
  WRITE( stdout,'(/,5x,"End of Effective Potential Calculation",/)')
  !
  ! Electric field calculation
  !
  lgamma = .true.
  CALL allocate_phq()
  ALLOCATE(ikks(nksq), ikqs(nksq))
  ikks(1)=1
  ikqs(1)=1
!  call newd()     ! don't remember why this routine is here. But it seem not to do
                   ! any thing.
!  call openfilq()
  !
  WRITE( stdout,'(/,5x,"Frequency Dependent Polarizability Calculation",/)')
  !
  !
  i = nfs
  freq_loop : DO WHILE ( i >= 1 )
     !
     CALL solve_e_vdw ( fiu(i) )
     IF ( convt ) CALL polariz_vdw ( fiu(i) )
     i = i - 1
     !
  ENDDO freq_loop
  !
  WRITE( stdout,'(/,5x,"End of Frequency Dependent Polarizability Calculation",/)')
  !
  CALL deallocate_phq()
  !
  !
  WRITE( stdout, '(/,5x,"End of vdW calculation")')
  !
  !
  CALL clean_pw( .true. )
  !
  CALL stop_vdw()
  !
  CALL stop_clock( 'VdW' )
  !
END PROGRAM vdw
!
!-----------------------------------------------------------------------
SUBROUTINE vdw_init ( )
  !-----------------------------------------------------------------------
  !
  !    This subroutine reads the data for the output file produced by pw.x
  !
  !    DESCRIPTION of the INPUT: see file Docs/INPUT_VdW
  !
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : bg
  USE ener,      ONLY : ef
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp, tau
  USE gvect
  USE fft_base,  ONLY : dfftp
  USE vlocal,    ONLY : strf
  USE io_files,  ONLY : tmp_dir, prefix
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE parser,    ONLY : read_line
  USE freq_ph
  USE control_vdw

  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: plot_num, kpoint, kband, spin_component, ios, flen
  LOGICAL :: stm_wfc_matching, lsign

  REAL(DP) :: emin, emax, sample_bias, z, dz, epsilon
  ! directory for temporary files
  CHARACTER(len=256) :: outdir
  CHARACTER(len=256)         :: input_line
  CHARACTER(len=80)          :: card
  CHARACTER(len=1), EXTERNAL :: capital
  LOGICAL                    :: tend
  LOGICAL                    :: end_of_file
  INTEGER                    :: i
  NAMELIST / inputvdw / outdir, prefix, tr2_vdw, thresh_veff, nmix_vdw, &
                        al_mix_vdw, niter_vdw
  ! prefix       : the prefix of files produced by pwscf
  ! outdir       : directory where input, output, temporary files reside
  ! tr2_vdw      : convergence threshold
  ! thresh_veff  : thresh_hold for iterative optimization of Veff
  ! nmix_vdw     : number of previous iterations used in mixing
  ! al_mix_vdw   : the mixing parameter
  ! niter_vdw    : maximum number of iterations
  !
  !   set default values for variables in namelist
  !
  CALL start_clock( 'vdw_init' )
  !
  prefix      = 'pwscf'
  outdir      = './'
  tr2_vdw     = 1.d-12
  thresh_veff = 4.d-3
  nmix_vdw    = 4
  al_mix_vdw  = 0.1d0
  niter_vdw   = 50
  !
  IF ( ionode )  THEN
     !
     !     reading the namelist inputpp
     !
     READ (5, inputvdw, err = 200, iostat = ios)
200  CALL errore ('vdw', 'reading inputvdw namelist', abs (ios) )
     tmp_dir = trimcheck ( outdir )
     !
     !   reading frequencies
     !
     READ (5, *, err = 10, iostat = ios) card
     READ (5, *, err = 10, iostat = ios) nfs
     ALLOCATE(fiu(nfs))
     DO i = 1, nfs
        READ (5, *, err = 10, iostat = ios) fiu(i)
     ENDDO
10   CALL errore ('extract', 'error or eof while reading frequencies', abs(ios) )
  ENDIF
  !
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( nfs, ionode_id )
  CALL mp_bcast( fiu, ionode_id )
  CALL mp_bcast( thresh_veff, ionode_id )
  CALL mp_bcast( tr2_vdw, ionode_id )
  CALL mp_bcast( nmix_vdw, ionode_id )
  CALL mp_bcast( al_mix_vdw, ionode_id )
  CALL mp_bcast( niter_vdw, ionode_id )
  !
  ! no task specified: do nothing and return
  !
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file ( )
  CALL openfil_pp ( )
  CALL struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
       strf, eigts1, eigts2, eigts3)
  CALL init_us_1 ( )
  !
  CALL stop_clock( 'vdw_init' )
  !
END SUBROUTINE vdw_init
