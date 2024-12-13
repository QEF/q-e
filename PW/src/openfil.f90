!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
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
  USE ldaU,             ONLY : lda_plus_U, Hubbard_projectors, nwfcU
  USE io_files,         ONLY : nwordwfc, iunhub, iunhub_noS, nwordwfcU, &
                               iunefield, iunefieldm, iunefieldp
  USE noncollin_module, ONLY : npol
  USE bp,               ONLY : lelfield
#if defined(__HDF5) && defined(__MPI) 
  USE hdf5_qe,          ONLY : initialize_hdf5
#endif 
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst
  !
  ! ... Files needed for DFT+U(+V)
  ! ... iunhub contains S times the (orthogonalized) atomic wfcs
  ! ...        (only atomic wfcs with a U correction, in variable wfcU)
  ! ... iunhub_nos as above, without S
  !
  ! ... nwordwfc  is the record length (IN COMPLEX WORDS)
  ! ...           for the direct-access file containing wavefunctions
  ! ... nwordwfcU as above for atomic/U-manifold wavefunctions
  !
  nwordwfc  = nbnd*npwx*npol
  nwordwfcU = npwx*nwfcU*npol
  !
  IF ( lda_plus_u .AND. (Hubbard_projectors.NE.'pseudo') ) THEN
     CALL open_buffer( iunhub,  'hub',  nwordwfcU, io_level, exst )
     IF (io_level>=1) CALL open_buffer( iunhub_noS,  'hubnoS',  nwordwfcU, io_level, exst )
  ENDIF
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
