!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-------------------------------------------------------------------
PROGRAM compute_self_hartree
  !-----------------------------------------------------------------
  !
  !!  This simple code read the PWSCF and Wannier output and 
  !!  compute the self-Hartree of the orbitals
  !!   
  !!  Code written by Nicola Colonna.
  !
  USE mp_global,             ONLY : mp_startup,mp_global_end
  USE environment,           ONLY : environment_start, environment_end
  USE check_stop,            ONLY : check_stop_init
  USE kinds,                 ONLY : DP
  USE klist,                 ONLY : nkstot
  USE lsda_mod,              ONLY : nspin
  USE io_global,             ONLY : ionode_id, stdout, ionode
  USE io_files,              ONLY : tmp_dir, prefix, check_tempdir
  USE control_kcw
  USE run_info,              ONLY : title
  USE mp_global,             ONLY : intra_image_comm 
  USE mp,                    ONLY : mp_bcast
  USE control_lr,            ONLY : lrpa
  USE input_parameters,      ONLY : assume_isolated
  !
  !
  IMPLICIT NONE
  COMPLEX(DP) :: sh_i
  INTEGER :: i
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER :: ios, nkstot_
  CHARACTER (LEN=256) :: outdir
  LOGICAL, EXTERNAL  :: imatches
  ! 
  NAMELIST / KCW_PP /    outdir, prefix, mp1, mp2, mp3, num_wann, seedname, kcw_iverbosity, &
                        l_vcut, assume_isolated, io_sp, io_real_space
  !
  ! prefix       : the prefix of files produced by pwscf
  ! outdir       : directory where input, output, temporary files reside
  ! num_wann     : number of occupied wannier
  ! the interpolation
  ! 
  CHARACTER(LEN=18) :: code='KC_PP_Self-Hartree'
  !
  CALL mp_startup ( )
  CALL environment_start ( code )
  !
  IF (ionode) THEN
     CALL input_from_file ( )
     READ( 5, '(A)', IOSTAT = ios ) title
  ENDIF
  !
  CALL mp_bcast(ios, ionode_id, intra_image_comm )
  CALL errore( 'KC_PP', 'reading title ', ABS( ios ) )
  ! 
  call mp_bcast ( title, ionode_id, intra_image_comm )
  !
  ! Rewind the input if the title is actually the beginning of inputph namelist
  IF( imatches("&kcw_pp", title)) THEN
    WRITE(stdout, '(6x,a)') "Title line not specified: using 'default'."
    title='default'
    REWIND(5, iostat=ios)
    CALL errore('conmpute_self_hartree', 'Title line missing from input.', abs(ios))
  ENDIF
  !
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  !
  prefix              = 'kcw_wann'
  seedname            = 'wann'
  num_wann            = 0 
  mp1                 = -1
  mp2                 = -1
  mp3                 = -1
  kcw_iverbosity       = 0
  l_vcut              = .false.
  assume_isolated     = "none" 
  io_sp               = .FALSE.
  io_real_space       = .FALSE.
  ! 
  ! ...  reading the namelist inputki
  !
  IF (ionode) READ( 5, KCW_PP, IOSTAT = ios )
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore( 'compute_self_hartree', 'reading KC_PP namelist', ABS( ios ) )
  !
  CALL input_pp_summary ()
  !
  CALL mp_bcast(outdir, ionode_id, intra_image_comm)
  CALL mp_bcast(prefix, ionode_id, intra_image_comm)
  CALL mp_bcast(seedname, ionode_id, intra_image_comm)
  CALL mp_bcast(num_wann, ionode_id, intra_image_comm)
  CALL mp_bcast(mp1, ionode_id, intra_image_comm)
  CALL mp_bcast(mp2, ionode_id, intra_image_comm)
  CALL mp_bcast(mp3, ionode_id, intra_image_comm)
  CALL mp_bcast(kcw_iverbosity, ionode_id, intra_image_comm)
  CALL mp_bcast(l_vcut, ionode_id, intra_image_comm)
  CALL mp_bcast(assume_isolated, ionode_id, intra_image_comm)
  CALL mp_bcast(io_sp, ionode_id, intra_image_comm)
  CALL mp_bcast(io_real_space, ionode_id, intra_image_comm)
  !
  !
  !
  tmp_dir = trimcheck (outdir)
  !
  tmp_dir_kcw= TRIM (tmp_dir) // 'kcw' //'/'
  ! 
  WRITE( stdout, '(5X,"INFO: Reading pwscf data")')
  CALL read_file ( )
  !
  IF (nspin == 4) THEN
    nkstot_ = nkstot
    nrho = 4
  ELSE
    nkstot_ = nkstot/nspin
    nrho = 1
  ENDIF
  !
  IF ( mp1*mp2*mp3 /= nkstot_ ) &
     CALL errore('compute_self_hartree', ' WRONG number of k points from input, check mp1, mp2, mp3', 1)
  !
  CALL sh_setup () 
  !
  lrpa = .true.
  WRITE( stdout,'(/,5X,"INFO: WANNIER orbital SH ",/)')
  DO i = 1, num_wann
    ! ... Compute the Self_hartree for each Wannier 
    !
    sh_i = CMPLX(0.D0, 0.D0, kind= DP)
    CALL self_hartree ( i, sh_i)
    WRITE(stdout,'(5X, "orb, Self hartree ", 1i5, 3x, 1F10.6)') i, REAL(sh_i)
    !
  ENDDO
  !
  ! Clean and Close 
  CALL mp_global_end()
  CALL environment_end( code )
  !
END PROGRAM compute_self_hartree
