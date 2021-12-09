!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE openfil()
  !----------------------------------------------------------------------------
  !! This routine opens some files needed to the self consistent run,
  !! sets various file names, units, record lengths. 
  !! All units are set in Modules/io_files.f90
  !
  USE kinds,            ONLY : DP
  USE buffers,          ONLY : open_buffer
  USE control_flags,    ONLY : io_level
  USE basis,            ONLY : natomwfc
  USE wvfct,            ONLY : nbnd, npwx
  USE fixed_occ,        ONLY : one_atom_occupations
  USE ldaU,             ONLY : lda_plus_U, U_projection, nwfcU
  USE io_files,         ONLY : prefix, iunpun, iunsat, &
                               iunhub, nwordwfcU, nwordwfc, nwordatwfc, &
                               iunefield, iunefieldm, iunefieldp, seqopn
  USE noncollin_module, ONLY : npol
  USE bp,               ONLY : lelfield
  USE wannier_new,      ONLY : use_wannier
#if defined(__HDF5) && defined(__MPI) 
  USE hdf5_qe,          ONLY : initialize_hdf5
#endif 
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst
  !
  ! ... Files needed for DFT+U(+V)
  ! ... iunsat contains the (orthogonalized) atomic wfcs * S
  ! ... iunhub  as above, only wfcs * S with a U correction
  !
  ! ... nwordwfc is the record length (IN COMPLEX WORDS)
  ! ... for the direct-access file containing wavefunctions
  ! ... nwordatwfc/nwordwfcU as above for atomic/U-manifold wavefunctions
  !
  nwordwfc  = nbnd*npwx*npol
  nwordatwfc= npwx*natomwfc*npol
  nwordwfcU = npwx*nwfcU*npol
  !
  IF ( lda_plus_u .AND. (U_projection.NE.'pseudo') ) THEN
     CALL open_buffer( iunhub,  'hub',  nwordwfcU, io_level, exst )
  ENDIF
  IF ( use_wannier .OR. one_atom_occupations ) &
     CALL open_buffer( iunsat, 'satwfc', nwordatwfc, io_level, exst )
  !
  ! ... open units for electric field calculations
  !
  IF ( lelfield ) THEN
      CALL open_buffer( iunefield , 'ewfc' , nwordwfc, io_level, exst )
      CALL open_buffer( iunefieldm, 'ewfcm', nwordwfc, io_level, exst )
      CALL open_buffer( iunefieldp, 'ewfcp', nwordwfc, io_level, exst )
  ENDIF
  !
#if defined(__HDF5) && defined(__MPI) 
  ! calls h5open_f mandatory in any application using hdf5
  CALL initialize_hdf5()
#endif 
  !
END SUBROUTINE openfil
