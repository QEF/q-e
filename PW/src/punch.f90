!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE punch( what )
  !----------------------------------------------------------------------------
  !
  ! ... This routine is called at the end of the run to save to a file
  ! ... the information needed for further processing (phonon etc.)
  ! ... what = 'all'          write xml data file, charge density, wavefunctions
  ! ...                       (for final data)
  ! ... what = 'config'       write xml data file and charge density; also,
  !                           for nks=1, wavefunctions in plain binary format
  !                           (see why in comments below)
  ! ...                       (for intermediate or incomplete results)
  ! ... what = 'init-config'  write xml data file only excluding final results
  ! ...                       (for dry run, can be called at early stages)
  !
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : iunpun, iunwfc, nwordwfc, diropn, &
       tmp_dir, prefix, postfix, create_directory
  USE control_flags,        ONLY : io_level, lscf
  USE klist,                ONLY : nks
  USE io_files,             ONLY : xmlpun_schema, psfile, pseudo_dir
  USE wrappers,             ONLY : f_copy
  USE spin_orb,             ONLY : lforcet
  USE scf,                  ONLY : rho
  USE lsda_mod,             ONLY : nspin
  USE ions_base,            ONLY : nsp
  USE funct,                ONLY : get_inlc
  USE kernel_table,         ONLY : vdw_table_name, kernel_file_name
  USE pw_restart_new,       ONLY : pw_write_schema, pw_write_binaries
  USE qexsd_module,         ONLY : qexsd_reset_steps
  USE io_rho_xml,           ONLY : write_scf
  USE a2F,                  ONLY : la2F, a2Fsave
  USE wavefunctions, ONLY : evc
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: what
  !
  LOGICAL :: exst, wf_collect
  CHARACTER(LEN=320) :: cp_source, cp_dest
  INTEGER            :: cp_status, nt, inlc
  !
  !
  IF (io_level < 0 ) RETURN
  !
  WRITE( UNIT = stdout, FMT = '(/,5X,"Writing output data file ",A)' ) &
      TRIM( prefix ) // postfix
  iunpun = 4
  !
  ! ...New-style I/O with xml schema and (optionally) hdf5 binaries
  !
  ! ... create the main restart directory (if needed)
  !
  CALL create_directory( TRIM( tmp_dir ) // TRIM( prefix ) // postfix )
  !
  ! ... wf_collect keeps track whether wfcs are written in portable format
  !
  wf_collect = ( TRIM(what) == 'all' )
  CALL pw_write_schema( what, wf_collect )
  !
  ! ... charge density - also writes rho%ns if lda+U and rho%bec if PAW
  ! ... do not overwrite the scf charge density with a non-scf one
  ! ... (except in the 'force theorem' calculation of MAE where the
  ! ...  charge density differs from the one read from disk)
  !
  IF (TRIM(what) == 'all' .OR. TRIM(what) == 'config' ) THEN
     IF ( lscf .OR. lforcet ) CALL write_scf( rho, nspin )
  END IF
  !
  IF (TRIM(what) == 'all') THEN 
     !
     ! ... make a copy of xml file one level up (FIXME: why?)
     !
     IF (ionode) THEN
        cp_source = TRIM(tmp_dir)//TRIM(prefix)//postfix//xmlpun_schema
        cp_dest   = TRIM(tmp_dir)//TRIM(prefix)//'.xml'
        cp_status = f_copy(cp_source, cp_dest)
     END IF
     !
     ! ... wavefunctions in "collected" format - also G- and k+G-vectors
     !
     CALL pw_write_binaries( )
     !
     ! ... copy pseudopotential files into the .save directory
     !
     DO nt = 1, nsp
        cp_source = TRIM(pseudo_dir)//psfile(nt)
        cp_dest   = TRIM(tmp_dir)//TRIM(prefix)//postfix//psfile(nt)
        IF ( TRIM(cp_source) /= TRIM(cp_dest) ) &
             cp_status = f_copy(cp_source, cp_dest)
     END DO
     !
     ! ... copy kernal table for vdW functionals if needed
     !
     inlc = get_inlc()
     IF ( inlc > 0 ) THEN 
        cp_source = TRIM(kernel_file_name)
        cp_dest = TRIM(tmp_dir)//TRIM(prefix)//postfix//TRIM(vdw_table_name)
        IF ( TRIM(cp_source) /= TRIM(cp_dest) ) & 
           cp_status = f_copy(cp_source, cp_dest)
     END IF  
     !
     ! ... if allocated, deallocate variables containing info on ionic steps 
     ! 
     CALL qexsd_reset_steps()
     !
  ELSE IF ( TRIM(what) == 'config' .AND.  nks == 1 ) THEN
     !
     ! ... here we are stopping an incomplete calculations - wavefunctions are 
     ! ... stored in buffers and saved when buffers are closed. For 1 k-point 
     ! ... however there is no buffer: wavefunctions must be saved to file here
     !
     IF (io_level < 1) CALL diropn( iunwfc, 'wfc', 2*nwordwfc, exst )
     CALL davcio ( evc, 2*nwordwfc, iunwfc, nks, 1 )
     IF (io_level < 1) CLOSE ( UNIT=iunwfc, STATUS='keep' )
     CALL infomsg('punch','wavefunctions written to file')
     !
  END IF
  !
  ! ... FIXME: for electron-phonon calculations - data should be read from xml file!
  ! 
  IF ( la2F ) CALL a2Fsave()
  !
  RETURN
  !
END SUBROUTINE punch
