!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
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
  USE control_kc_wann
  USE run_info,              ONLY : title
  USE mp_global,             ONLY : intra_image_comm 
  USE mp,                    ONLY : mp_bcast
  USE control_ph,            ONLY : tmp_dir_ph
  USE mp_images,             ONLY : my_image_id
  USE control_lr,            ONLY : lrpa
  USE input_parameters,      ONLY : assume_isolated
  !
  !
  IMPLICIT NONE
  COMPLEX(DP) :: sh_i
  INTEGER :: i
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER :: ios
  CHARACTER (LEN=256) :: outdir
  LOGICAL, EXTERNAL  :: imatches
  CHARACTER(LEN=6)   :: int_to_char
  !
  CHARACTER(LEN=18) :: code='KC_PP_Self-Hartree'
  !
  CALL mp_startup ( )
  CALL environment_start ( code )
  ! 
  NAMELIST / KC_PP /    outdir, prefix, mp1, mp2, mp3, num_wann, seedname, kc_iverbosity, &
                        l_vcut, assume_isolated
  !
  ! prefix       : the prefix of files produced by pwscf
  ! outdir       : directory where input, output, temporary files reside
  ! num_wann     : number of occupied wannier
  ! the interpolation
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
  IF( imatches("&kc_pp", title)) THEN
    WRITE(*, '(6x,a)') "Title line not specified: using 'default'."
    title='default'
    REWIND(5, iostat=ios)
    CALL errore('KC_readin', 'Title line missing from input.', abs(ios))
  ENDIF
  !
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  !
  prefix              = 'kc_wann'
  seedname            = 'wann'
  num_wann            = 0 
  mp1                 = -1
  mp2                 = -1
  mp3                 = -1
  kc_iverbosity       = 0
  l_vcut              = .false.
  assume_isolated     = "none" 
  ! 
  ! ...  reading the namelist inputki
  !
  IF (ionode) READ( 5, KC_PP, IOSTAT = ios )
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  CALL errore( 'kc_pp_readin', 'reading KC_PP namelist', ABS( ios ) )
  !
  CALL mp_bcast(outdir, ionode_id, intra_image_comm)
  CALL mp_bcast(prefix, ionode_id, intra_image_comm)
  CALL mp_bcast(seedname, ionode_id, intra_image_comm)
  CALL mp_bcast(num_wann, ionode_id, intra_image_comm)
  CALL mp_bcast(mp1, ionode_id, intra_image_comm)
  CALL mp_bcast(mp2, ionode_id, intra_image_comm)
  CALL mp_bcast(mp3, ionode_id, intra_image_comm)
  CALL mp_bcast(kc_iverbosity, ionode_id, intra_image_comm)
  CALL mp_bcast(l_vcut, ionode_id, intra_image_comm)
  CALL mp_bcast(assume_isolated, ionode_id, intra_image_comm)
  !
  !
  !
  tmp_dir = trimcheck (outdir)
  !
  tmp_dir_ph= TRIM (tmp_dir) // '_ph' // TRIM(int_to_char(my_image_id)) //'/'
  ! 
  WRITE( stdout, '(5X,"INFO: Reading pwscf data")')
  CALL read_file ( )
  !
  IF (calculation /= 'wann2kc' .AND. (mp1*mp2*mp3 /= nkstot/nspin) ) &
     CALL errore('kc_readin', ' WRONG number of k points from input, check mp1, mp2, mp3', 1)
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
