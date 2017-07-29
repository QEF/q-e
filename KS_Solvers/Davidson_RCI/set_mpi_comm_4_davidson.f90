!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE set_mpi_comm_4_davidson_rci(parent_comm, intra_bgrp_comm_, inter_bgrp_comm_ )
  !----------------------------------------------------------------------------
  !
  USE david_param,      ONLY : DP
  USE mp_bands_davidson
  USE mp,               ONLY : mp_size, mp_rank
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: parent_comm, intra_bgrp_comm_, inter_bgrp_comm_
  ! local variables
  INTEGER :: parent_nproc, parent_mype, ortho_parent_comm_
  !
  intra_bgrp_comm   = intra_bgrp_comm_
  inter_bgrp_comm   = inter_bgrp_comm_
  !
#if defined (__MPI)
    !
    parent_nproc = mp_size( parent_comm )
    parent_mype  = mp_rank( parent_comm )
    !
    ! ... Set number of processors per band group
    !
    nproc_bgrp = mp_size( intra_bgrp_comm )
    !
    nbgrp = parent_nproc / nproc_bgrp

    IF ( nbgrp < 1 .OR. nbgrp > parent_nproc ) CALL errore( 'mp_start_bands',&
                          'invalid number of band groups, out of range', 1 )
    IF ( MOD( parent_nproc, nbgrp ) /= 0 ) CALL errore( 'mp_start_bands', &
        'n. of band groups  must be divisor of parent_nproc', 1 )
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
    CALL mp_barrier( parent_comm )

    !
#endif
    RETURN
  !
END SUBROUTINE set_mpi_comm_4_davidson_rci
!----------------------------------------------------------------------------
SUBROUTINE unset_mpi_comm_4_davidson_rci()
  !----------------------------------------------------------------------------
  !
  use mp_diag
  IMPLICIT NONE
#if defined (__MPI)
  CALL clean_ortho_group ( )
#endif
    RETURN
  !
END SUBROUTINE unset_mpi_comm_4_davidson_rci
