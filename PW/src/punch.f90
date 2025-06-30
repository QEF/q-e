!
! Copyright (C) 2001-2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE punch( what )
  !----------------------------------------------------------------------------
  !! This routine is called at the end of the run to save on a file
  !! the information needed for further processing (phonon etc.).
  !
  !! * what = 'all' : write xml data file and, for normal or high I/O options,
  !!                  tha charge density and wavefunctions.
  !!                  Used for writing final data.
  !! * what = 'config' : write xml data file and, for low I/O options,
  !!                     the charge density; also, for nks=1, wavefunctions
  !!                     in plain binary format (see comments below for reason) 
  !!                     Used for writing intermediate or incomplete results.
  !! * what = 'config-only' : write xml data file only
  !!                          Used for saving restart information.
  !! * what = 'config-init' : write xml data file only excluding final results
  !!                          Used for dry run, at the start of calculation.
  !
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : iunwfc, nwordwfc, diropn, tmp_dir, prefix, &
                                   restart_dir, xmlfile, create_directory
  USE control_flags,        ONLY : io_level, lscf, lxdm
  USE klist,                ONLY : nks
  USE io_files,             ONLY : psfile, pseudo_dir, molfile
  USE clib_wrappers,        ONLY : f_copy
  USE noncollin_module,     ONLY : lforcet
  USE scf,                  ONLY : rho
  USE lsda_mod,             ONLY : nspin
  USE ions_base,            ONLY : nsp
  USE pw_restart_new,       ONLY : pw_write_schema, write_collected_wfc
  USE qexsd_module,         ONLY : qexsd_reset_steps
  USE io_rho_xml,           ONLY : write_scf
  USE xc_lib,               ONLY : xclib_dft_is
  USE a2F,                  ONLY : la2F, a2Fsave
  USE wavefunctions,        ONLY : evc
  USE xdm_module,           ONLY : write_xdmdat
  USE rism3d_facade,        ONLY : lrism3d, rism3d_write_to_restart
  USE solvmol,              ONLY : nsolV
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: what
  !! what is what? see main comment
  !
  ! ... local variables
  !
  LOGICAL :: exst, only_init, wf_collect
  CHARACTER(LEN=320) :: cp_source, cp_dest
  INTEGER            :: cp_status, nt, isolV
  !
  !
  IF ( io_level < -2 ) THEN
     !! disk_io='none': no I/O at all, do not write anything
     CALL infomsg('punch',"Nothing written, as per user request (disk_io='none')")
     RETURN
  ELSE
     WRITE( stdout, '(/,5X,"Writing ",A," to output data dir ",A," :")' ) &
          TRIM ( what ), TRIM ( restart_dir ( ) )
  END IF
  !
  ! ... create the main restart directory (if needed)
  !
  CALL create_directory( restart_dir( ) )
  !
  ! ... wf_collect keeps track whether wfcs are written in portable format
  !
  wf_collect = ( TRIM(what) == 'all' )
  only_init  = ( TRIM(what) == 'config-init' )
  CALL pw_write_schema( only_init, wf_collect )
  WRITE( stdout, '(5X,"XML data file")', ADVANCE ='NO' )
  !
  ! ... nothing left to do for 'minimal' I/O option: return
  !
  IF ( io_level < -1 ) RETURN
  !
  IF (TRIM(what) == 'all' .OR. TRIM(what) == 'config' ) THEN
     !
     ! ... charge density - also writes rho%ns if lda+U and rho%bec if PAW
     ! ... do not overwrite the scf charge density with a non-scf one
     ! ... (except in the 'force theorem' calculation of MAE where the
     ! ...  charge density differs from the one read from disk)
     !
     IF ( lscf .OR. lforcet ) THEN
        CALL write_scf( rho, nspin )
        WRITE( stdout, '(", charge density")', ADVANCE ='NO' )
        IF ( xclib_dft_is('meta') ) &
           WRITE( stdout, '(" with meta-GGA kinetic term")', ADVANCE ='NO' )
     END IF 
     !
     ! ... correlation functions of 3D-RISM.
     ! ... do not overwrite them, in case of non-scf
     !
     IF ( lrism3d .AND. lscf ) THEN
        CALL rism3d_write_to_restart()
        WRITE( stdout, '(", RISM restart")', ADVANCE ='NO' )
     END IF
     !
  END IF
  !
  ! ... nothing left to do for 'nowf' I/O option: return
  !
  IF ( io_level < 0 ) RETURN
  !
  IF (TRIM(what) == 'all') THEN 
     !
     ! ... copy xml file one level up (FIXME: why?),
     ! ... copy pseudopotential files into the .save directory
     ! ... copy molecular files into the .save directory (for 3D-RISM)
     !
     IF (ionode) THEN
        !
        cp_source = xmlfile ( )
        cp_dest   = TRIM(tmp_dir)//TRIM(prefix)//'.xml'
        cp_status = f_copy(cp_source, cp_dest)
        !
        DO nt = 1, nsp
           cp_source = TRIM(pseudo_dir)//psfile(nt)
           cp_dest   = TRIM(restart_dir ( ) ) //psfile(nt)
           IF ( TRIM(cp_source) /= TRIM(cp_dest) ) &
                cp_status = f_copy(cp_source, cp_dest)
        ENDDO
        WRITE( stdout, '(", pseudopotentials")', ADVANCE ='NO' )
        !
        IF ( lrism3d ) THEN
           !
           DO isolV = 1, nsolV
              cp_source = TRIM(pseudo_dir)//molfile(isolV)
              cp_dest   = TRIM(restart_dir ( ) ) //molfile(isolV)
              IF ( TRIM(cp_source) /= TRIM(cp_dest) ) &
                   cp_status = f_copy(cp_source, cp_dest)
           ENDDO
           WRITE( stdout, '(", RISM molfiles")', ADVANCE ='NO' )
           !
        ENDIF
        !
        ! write XDM dispersion data (coefficients and vdw radii) to xdm.dat
        IF (lxdm) THEN
           CALL write_xdmdat()
           WRITE( stdout, '(", XDM data")', ADVANCE ='NO' )
        ENDIF
     ENDIF
     !
     ! ... wavefunctions in "collected" format - also G- and k+G-vectors
     !
     CALL write_collected_wfc( )
     WRITE( stdout, '(", collected wavefunctions")', ADVANCE ='NO' )
     !
     ! ... if allocated, deallocate variables containing info on ionic steps 
     ! 
     CALL qexsd_reset_steps()
     !
  ELSEIF ( TRIM(what) == 'config' .AND.  nks == 1 ) THEN
     !
     ! ... here we are stopping an incomplete calculations - wavefunctions are 
     ! ... stored in buffers and saved when buffers are closed. For 1 k-point 
     ! ... however there is no buffer: wavefunctions must be saved to file here
     ! ... (if disk_io='low' no file is open, must be opened and closed here)
     !
     IF (io_level < 1) CALL diropn( iunwfc, 'wfc', 2*nwordwfc, exst )
     CALL davcio ( evc, 2*nwordwfc, iunwfc, nks, 1 )
     IF (io_level < 1) CLOSE ( UNIT=iunwfc, STATUS='keep' )
     WRITE( stdout, '(" distributed wavefunctions")', ADVANCE ='NO' )
     !
  ENDIF
  WRITE( stdout, *)
  !
  ! ... FIXME: for electron-phonon calculations - data should be read from xml file!
  ! 
  IF ( la2F ) CALL a2Fsave()
  !
  RETURN
  !
END SUBROUTINE punch
