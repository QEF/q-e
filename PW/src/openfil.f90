!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE openfil()
  !----------------------------------------------------------------------------
  !
  ! ... This routine opens some files needed to the self consistent run,
  ! ... sets various file names, units, record lengths
  ! ... All units are set in Modules/io_files.f90
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE basis,            ONLY : natomwfc, starting_wfc
  USE wvfct,            ONLY : nbnd, npwx
  USE fixed_occ,        ONLY : one_atom_occupations
  USE klist,            ONLY : nks
  USE ldaU,             ONLY : lda_plus_U, U_projection
  USE io_files,         ONLY : prefix, iunpun, iunat, iunsat, iunwfc, iunigk, &
                               nwordwfc, nwordatwfc, iunefield, diropn, &
                               tmp_dir, wfc_dir, iunefieldm, iunefieldp, seqopn
  USE pw_restart,       ONLY : pw_readfile
  USE noncollin_module, ONLY : npol
  USE bp,               ONLY : lelfield
  USE buffers,          ONLY : open_buffer, init_buffer
  USE control_flags,    ONLY : io_level, twfcollect
  USE wannier_new,            ONLY : use_wannier
  !
  IMPLICIT NONE
  !
  LOGICAL            :: exst
  INTEGER            :: ierr
  CHARACTER(LEN=256) :: tmp_dir_save
  !
  ! ... tmp_dir may be replaced by wfc_dir  for large files
  !
  tmp_dir_save = tmp_dir
  !
  IF ( wfc_dir /= 'undefined' ) THEN
     !
     WRITE( stdout, '(5X,"writing wfc files to a dedicated directory")' )
     !
     tmp_dir = wfc_dir
     !
  END IF
  !
  ! ... nwordwfc is the record length (IN COMPLEX WORDS)
  ! ... for the direct-access file containing wavefunctions
  !
  nwordwfc = nbnd*npwx*npol
  !
  ! ... iunwfc=10: read/write wfc from/to file
  ! ... iunwfc=-1: copy wfc to/from RAM 
  !
  IF ( io_level > 0 ) THEN
     iunwfc = 10
  ELSE
     iunwfc = -1
  END IF
  CALL open_buffer( iunwfc, 'wfc', nwordwfc, nks, exst )
  !
  IF ( TRIM(starting_wfc) == 'file' .AND. .NOT. exst)  THEN
     !
     ierr = 1 
     IF ( twfcollect ) THEN
        !
        ! ... wavefunctions are read from the "save" file and rewritten
        ! ... (directly in pw_readfile) using the internal format
        !
        CALL pw_readfile( 'wave', ierr )
        !
     ELSE
        !
        ! ... wavefunctions are read into memory
        !
        CALL init_buffer ( iunwfc, exst, ierr )
        !
     END IF

     IF ( ierr > 0 ) THEN
        !
        WRITE( stdout, '(5X,"Cannot read wfc : file not found")' )
        !
        starting_wfc = 'atomic'
        !
     END IF
     !
  END IF
  !
  ! ... Needed for LDA+U
  !
  ! ... iunat  contains the (orthogonalized) atomic wfcs 
  ! ... iunsat contains the (orthogonalized) atomic wfcs * S
  ! ... iunocc contains the atomic occupations computed in new_ns
  ! ... it is opened and closed for each reading-writing operation  
  !
  nwordatwfc = 2*npwx*natomwfc*npol
  !
  IF ( ( lda_plus_u .AND. (U_projection.NE.'pseudo') ) .OR. &
        use_wannier .OR. one_atom_occupations ) THEN
     CALL diropn( iunat,  'atwfc',  nwordatwfc, exst )
     CALL diropn( iunsat, 'satwfc', nwordatwfc, exst )
  END IF
  !
  ! ... iunigk contains the number of PW and the indices igk
  ! ... Note that unit 15 is reserved for error messages 
  !
  CALL seqopn( iunigk, 'igk', 'UNFORMATTED', exst )
  !
  ! ... open units for electric field calculations
  !
  IF ( lelfield ) THEN
      CALL diropn( iunefield , 'ewfc' , 2*nwordwfc, exst )
      CALL diropn( iunefieldm, 'ewfcm', 2*nwordwfc, exst )
      CALL diropn( iunefieldp, 'ewfcp', 2*nwordwfc, exst )
  END IF
  !
  tmp_dir = tmp_dir_save
  !
  RETURN
  !
END SUBROUTINE openfil
