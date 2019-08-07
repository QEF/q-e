!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE set_mpi_comm_4_solvers(parent_comm, intra_bgrp_comm_, inter_bgrp_comm_ )
  !----------------------------------------------------------------------------
  !
  USE mp_bands_util
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: parent_comm, intra_bgrp_comm_, inter_bgrp_comm_
  ! local variables
  INTEGER :: parent_nproc, parent_mype, ortho_parent_comm_, ierr
  !
  !write(*,*) ' enter set_mpi_comm_4_davidson'
  intra_bgrp_comm   = intra_bgrp_comm_
  inter_bgrp_comm   = inter_bgrp_comm_
  !
#if defined (__MPI)
    !
    CALL mpi_comm_size(parent_comm,parent_nproc,ierr)
    IF (ierr/=0) CALL errore( ' set_mpi_comm_4_solvers ', ' problem getting MPI size ', 1 )
    CALL mpi_comm_rank(parent_comm,parent_mype,ierr)
    IF (ierr/=0) CALL errore( ' set_mpi_comm_4_solvers ', ' problem getting MPI rank ', 1 )
    !
    ! ... Set number of processors per band group
    !
    CALL mpi_comm_size(intra_bgrp_comm,nproc_bgrp,ierr)
    IF (ierr/=0) CALL errore( ' set_mpi_comm_4_solvers ', ' problem getting MPI size ', 1 )
    !
    nbgrp = parent_nproc / nproc_bgrp

    IF ( nbgrp < 1 .OR. nbgrp > parent_nproc ) &
       CALL errore( 'set_mpi_comm_4_solvers','invalid number of band groups, out of range', 1 )
    IF ( MOD( parent_nproc, nbgrp ) /= 0 ) &
       CALL errore( 'set_mpi_comm_4_solvers','n. of band groups  must be divisor of parent_nproc', 1 )
    !
    ! set logical flag so that band parallelization in H\psi is allowed
    ! (can be disabled before calling H\psi if not desired)
    !
    use_bgrp_in_hpsi = ( nbgrp > 1 )
    !
    ! ... set index of band group for this processor   ( 0 : nbgrp - 1 )
    !
    my_bgrp_id = parent_mype / nproc_bgrp
    !
    ! ... set index of processor within the image ( 0 : nproc_image - 1 )
    !
    me_bgrp    = MOD( parent_mype, nproc_bgrp )
    !
    CALL mpi_barrier( parent_comm, ierr )
    IF (ierr/=0) &
       CALL errore( 'set_mpi_comm_4_solvers','n. of band groups  must be divisor of parent_nproc', 1 )
    !
#else
    parent_nproc = 1
    parent_mype  = 0
    nproc_bgrp   = 1
    nbgrp        = 1 
    use_bgrp_in_hpsi = .false.
    my_bgrp_id   = 0
    me_bgrp      = 0
#endif
    !write(*,*) ' exit set_mpi_comm_4_davidson'
    RETURN
  !
END SUBROUTINE set_mpi_comm_4_solvers
