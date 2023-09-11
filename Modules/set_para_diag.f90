!
! Copyright (C) 2021 Quantum ESPRESSO Fondation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE set_para_diag( nbnd, use_para_diag )
  !-----------------------------------------------------------------------------
  !! Sets up the communicator used for parallel diagonalization in LAXlib.
  !! Merges previous code executed at startup with function "check_para_diag".
  !! To be called after the initialization of variables is completed and
  !! the dimension of matrices to be diagonalized is known
  !
  USE io_global,            ONLY : stdout, ionode
  USE mp_bands,             ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp_pools,             ONLY : intra_pool_comm
  USE mp_world,             ONLY : world_comm
  USE mp_exx,               ONLY : negrp
  USE command_line_options, ONLY : ndiag_

  IMPLICIT NONE

  INCLUDE 'laxlib.fh'

  INTEGER, INTENT(IN) :: nbnd
  !! dimension of matrices to be diagonalized (number of bands)
  LOGICAL, INTENT(INOUT) :: use_para_diag
  !! true if parallel linear algebra is to be used
  !
  LOGICAL, SAVE :: first = .TRUE.
  LOGICAL :: do_diag_in_band_group = .TRUE.
  INTEGER :: np_ortho(2), ortho_parent_comm

  IF( .NOT. first ) RETURN
  first = .FALSE.
  !
  IF( negrp > 1 .OR. do_diag_in_band_group ) THEN
     ! one diag group per bgrp with strict hierarchy: POOL > BAND > DIAG
     ! if using exx groups from mp_exx,  always use this diag method
     CALL laxlib_start ( ndiag_, intra_bgrp_comm, .TRUE. )
  ELSE
     ! one diag group per pool ( individual k-point level )
     ! with band group and diag group both being children of POOL comm
     CALL laxlib_start ( ndiag_, intra_pool_comm, .FALSE. )
  END IF
  CALL set_mpi_comm_4_solvers( intra_pool_comm, intra_bgrp_comm, &
       inter_bgrp_comm )
  !
#if defined(__MPI)
  CALL laxlib_getval( np_ortho = np_ortho, ortho_parent_comm = ortho_parent_comm )
  !
  IF( np_ortho(1) > nbnd ) &
     CALL errore ('set_para_diag', 'Too few bands for required ndiag',nbnd)
  !
  use_para_diag = ( np_ortho( 1 ) > 1 .AND. np_ortho( 2 ) > 1 )
  !
  IF ( ionode ) THEN
     !
     WRITE( stdout, '(5X,"Subspace diagonalization in iterative solution ",&
                     &   "of the eigenvalue problem:")' )
     IF ( use_para_diag ) THEN
        IF (ortho_parent_comm == intra_pool_comm) THEN
           WRITE( stdout, '(5X,"one sub-group per k-point group (pool) will be used")' )
        ELSE IF (ortho_parent_comm == intra_bgrp_comm) THEN
           WRITE( stdout, '(5X,"one sub-group per band group will be used")' )
        ELSE
           CALL errore( 'setup','Unexpected sub-group communicator ', 1 )
        END IF
#if defined(__ELPA) || defined(__ELPA_2015) || defined(__ELPA_2016)
        WRITE( stdout, '(5X,"ELPA distributed-memory algorithm ", &
              & "(size of sub-group: ", I2, "*", I3, " procs)",/)') &
               np_ortho(1), np_ortho(2)
#elif defined(__SCALAPACK)
        WRITE( stdout, '(5X,"scalapack distributed-memory algorithm ", &
              & "(size of sub-group: ", I2, "*", I3, " procs)",/)') &
               np_ortho(1), np_ortho(2)
#else
        WRITE( stdout, '(5X,"custom distributed-memory algorithm ", &
              & "(size of sub-group: ", I2, "*", I3, " procs)",/)') &
               np_ortho(1), np_ortho(2)
#endif
     ELSE
        WRITE( stdout, '(5X,"a serial algorithm will be used",/)' )
     END IF
     !
  END IF
  !
#else
  !
  use_para_diag = .FALSE.
  !
#endif
  !
  RETURN
END SUBROUTINE set_para_diag
