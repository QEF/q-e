!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE kc_pp_readin()
  !----------------------------------------------------------------------------
  !
  !    This routine reads the control variables for the program odd_apha.x
  !    from standard input (unit 5).
  !    A second routine readfile reads the variables saved on a file
  !    by the self-consistent program.
  !
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : ionode_id, stdout
  USE mp,            ONLY : mp_bcast
  USE klist,         ONLY : nkstot
  USE lsda_mod,      ONLY : nspin
  USE run_info,      ONLY : title
  USE io_files,      ONLY : tmp_dir, prefix, check_tempdir
  USE noncollin_module, ONLY : noncolin
  USE read_cards_module, ONLY : read_cards
  USE io_global,     ONLY : ionode
  USE mp_global,     ONLY :  intra_image_comm 
  USE paw_variables, ONLY : okpaw
  USE uspp,          ONLY : okvan
  USE control_kc_wann
  USE control_flags, ONLY : iverbosity
  USE symm_base,     ONLY : nrot
  USE mp_pools,      ONLY : npool
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ios 
    ! integer variable for I/O control
    ! counter on polarizations
    ! counter on iterations
    ! counter on atoms
    ! counter on types
  CHARACTER (LEN=256) :: outdir
  !
  CHARACTER(LEN=1), EXTERNAL :: capital
  INTEGER, EXTERNAL  :: atomic_number
  REAL(DP), EXTERNAL :: atom_weight
  LOGICAL, EXTERNAL  :: imatches
  LOGICAL, EXTERNAL  :: has_xml
  !
  ! kc_iverbosity   : verbosity control
  ! modenum      : single mode calculation
  ! fildyn       : output file for the dynamical matrix
  ! fildvscf     : output file containing deltavsc
  ! fildrho      : output file containing deltarho
  ! start_q      : in q list does the q points from start_q to last_q
  ! last_q       : 
  ! start_irr    : does the irred. representation from start_irr to last_irr
  ! last_irr     : 
  ! nogg         : if .true. lgamma_gamma tricks are not used
  !
  NAMELIST / KC_PP /    outdir, prefix, mp1, mp2, mp3, num_wann, seedname, use_ws_distance, &
                        num_wann_occ, num_wann_emp
  !
  ! prefix       : the prefix of files produced by pwscf
  ! outdir       : directory where input, output, temporary files reside
  ! seedname     : seedname for the Wannier calculation
  ! num_wann     : number of occupied wannier
  ! use_ws_distance : as in W90, if .true. the Wannier centers are considered in the interpolation
  ! 
  IF (ionode) THEN
  !
  ! ... Input from file ?
  !
     CALL input_from_file ( )
  !
  ! ... Read the first line of the input file
  !
     READ( 5, '(A)', IOSTAT = ios ) title
  !
  ENDIF
  !
  CALL mp_bcast(ios, ionode_id, intra_image_comm )
  CALL errore( 'KC_PP', 'reading title ', ABS( ios ) )
  ! 
  call mp_bcast ( title, ionode_id, intra_image_comm )
  !
  ! Rewind the input if the title is actually the beginning of inputph namelist
  IF( imatches("&kc_pp", title)) THEN
    WRITE(*, '(6x,a)') "Title line not specified: using 'default'."
    title='default'
    REWIND(5, iostat=ios)
    CALL errore('KC_readin', 'Title line missing from input.', abs(ios))
  ENDIF
  !
  ! ... set default values for variables in namelist
  !
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  !
  prefix              = 'kc_wann'
  seedname            = 'wann'
  num_wann_occ        = 0
  num_wann_emp        = 0 
  num_wann            = 0 
  mp1                 = -1
  mp2                 = -1
  mp3                 = -1
  use_ws_distance     = .TRUE.
  ! 
  ! ...  reading the namelist inputki
  !
  IF (ionode) READ( 5, KC_PP, IOSTAT = ios )
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore( 'kc_pp_readin', 'reading KC_PP namelist', ABS( ios ) )
  !
  IF (ionode) tmp_dir = trimcheck (outdir)
  !
  ! ... broadcasting all input variables to other nodes
  !
  CALL input_pp_summary ()
  CALL bcast_kc_pp_input ( ) 
  !
  ! .. READ the card with the K-point PATH
  !
  CALL read_cards( 'PW' )
  CALL convert_kpts_names( )
  !
  !
  IF (kc_iverbosity .gt. 1) iverbosity = 1
  
  IF (npool .gt. 1) &
      CALL errore ('kc_readin','Pool Parallelization not implemented. Re-run without pools.',1)
  !
  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !
  !   amass will also be read from file:
  !   save its content in auxiliary variables
  !
  !
  ! read data produced by pwscf
  !
  WRITE( stdout, '(5X,"INFO: Reading pwscf data")')
  CALL read_file ( )
  !
  CALL read_hr ( )
  !
  IF (okpaw.or.okvan) CALL errore('kc_readin',&
     'The kc_wann code with US or PAW is not available yet',1)

  IF (noncolin) CALL errore('kc_readin',&
   'The kc_wann code with non colliner spin is not available yet',1)

  IF (nrot .gt. 1 .AND. nkstot/nspin .gt. 1) &
    CALL errore('kc_readin','kc_wann requires nosym to be true ',nrot)
  !
  RETURN
  !
END SUBROUTINE kc_pp_readin
