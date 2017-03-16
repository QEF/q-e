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
  !
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : iunpun, iunwfc, nwordwfc, diropn, &
       tmp_dir, prefix
  USE control_flags,        ONLY : io_level, twfcollect, io_level, lscf
  USE klist,                ONLY : nks
  USE io_files,             ONLY : xmlpun_schema, psfile, pseudo_dir
  USE wrappers,             ONLY : f_copy
  USE xml_io_base,          ONLY : create_directory
  USE io_rho_xml,           ONLY : write_rho
  USE spin_orb,             ONLY : lforcet
  USE scf,                  ONLY : rho
  USE lsda_mod,             ONLY : nspin
  USE ions_base,            ONLY : nsp
  USE funct,                ONLY : get_inlc
  USE kernel_table,         ONLY : vdw_table_name, kernel_file_name
#if defined (__OLDXML) 
  USE pw_restart,           ONLY : pw_writefile
#else
  USE pw_restart_new,       ONLY : pw_write_schema, pw_write_binaries
#endif
  USE a2F,                  ONLY : la2F, a2Fsave
  USE wavefunctions_module, ONLY : evc
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*) :: what
  LOGICAL :: exst
  CHARACTER(LEN=320) :: cp_source, cp_dest
  INTEGER            :: cp_status, nt, inlc
  !
  !
  IF (io_level < 0 ) RETURN
  !
  WRITE( UNIT = stdout, FMT = '(/,5X,"Writing output data file ",A)' ) &
      TRIM( prefix ) // '.save'
  !
  ! ... if wavefunctions are stored in "distributed" format,
  ! ... save here wavefunctions to file if never saved before
  !
  IF ( .NOT. twfcollect .AND. nks == 1 ) THEN
     IF (io_level < 1) CALL diropn( iunwfc, 'wfc', 2*nwordwfc, exst )
     CALL davcio ( evc, 2*nwordwfc, iunwfc, nks, 1 )
     IF (io_level < 1) CLOSE ( UNIT=iunwfc, STATUS='keep' )
  END IF
  iunpun = 4
  !
#if defined(__OLDXML)
  !
  CALL pw_writefile( TRIM( what ) )
  !
#else
  !
  ! ...New-style I/O with xml schema and (optionally) hdf5 binaries
  !
  ! ... create the main restart directory (if needed)
  !
  CALL create_directory( TRIM( tmp_dir ) // TRIM( prefix ) // '.save' )
  !
  CALL pw_write_schema( )
  !
  ! ... charge density - also writes rho%ns if lda+U and rho%bec if PAW
  ! ... do not overwrite the scf charge density with a non-scf one
  ! ... (except in the 'force theorem' calculation of MAE where the
  ! ...  charge density differs from the one read from disk)
  !
  IF ( lscf .OR. lforcet ) CALL write_rho( rho, nspin )
  !
  IF (TRIM(what) == 'all') THEN 
     !
     ! ... make a copy of xml file one level up (FIXME: why?)
     !
     IF (ionode) THEN
        cp_source = TRIM(tmp_dir)//TRIM(prefix)//'.save/'//xmlpun_schema
        cp_dest   = TRIM(tmp_dir)//TRIM(prefix)//'.xml'
        cp_status = f_copy(cp_source, cp_dest)
     END IF
     !
     ! ... wavefunctions in "collected" format - also G- and k+G-vectors
     !
     IF ( twfcollect ) CALL pw_write_binaries( )
     !
     ! ... copy pseudopotential files into the .save directory
     !
     DO nt = 1, nsp
        cp_source = TRIM(pseudo_dir)//psfile(nt)
        cp_dest   = TRIM(tmp_dir)//TRIM(prefix)//'.save/'//psfile(nt)
        IF ( TRIM(cp_source) /= TRIM(cp_dest) ) &
             cp_status = f_copy(cp_source, cp_dest)
     END DO
     inlc = get_inlc()
     IF ( inlc > 0 ) THEN 
        cp_source = TRIM(kernel_file_name)
        cp_dest = TRIM(tmp_dir)//TRIM(prefix)//'.save/'//TRIM(vdw_table_name)
        IF ( TRIM(cp_source) /= TRIM(cp_dest) ) & 
           cp_status = f_copy(cp_source, cp_dest)
     END IF  
      !
  END IF
  !
#endif
  !
  IF ( la2F ) CALL a2Fsave()
  !
  RETURN
  !
END SUBROUTINE punch
