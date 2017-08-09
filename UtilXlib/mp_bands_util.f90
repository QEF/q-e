!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_bands_util
  !----------------------------------------------------------------------------
  !
  USE mp, ONLY : mp_barrier, mp_bcast, mp_size, mp_rank, mp_comm_split
  USE parallel_include
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... Band groups (processors within a pool of bands)
  ! ... Subdivision of pool group, used for parallelization over bands
  !
  INTEGER :: nbgrp       = 1  ! number of band groups
  INTEGER :: nproc_bgrp  = 1  ! number of processors within a band group
  INTEGER :: me_bgrp     = 0  ! index of the processor within a band group
  INTEGER :: root_bgrp   = 0  ! index of the root processor within a band group
  INTEGER :: my_bgrp_id  = 0  ! index of my band group
  INTEGER :: root_bgrp_id     = 0  ! index of root band group
  INTEGER :: inter_bgrp_comm  = 0  ! inter band group communicator
  INTEGER :: intra_bgrp_comm  = 0  ! intra band group communicator  
  ! Next variable is .T. if band parallelization is performed inside H\psi 
  ! and S\psi, .F. otherwise (band parallelization can be performed outside
  ! H\psi and S\psi, though)  
  LOGICAL :: use_bgrp_in_hpsi = .FALSE.
  !
  ! ... "task" groups (for band parallelization of FFT)
  !
  INTEGER :: ntask_groups = 1  ! number of proc. in an orbital "task group"
  !
#if defined (__MPI)
  ! ...  variable gstart2 must be properly set to be used with the real algorithm
  !
  INTEGER :: gstart = -1  ! variable not set yet 
#else
  INTEGER :: gstart =  2  ! appropriate value for serial execution
#endif

CONTAINS
  !
  SUBROUTINE set_bgrp_indices(nbnd, ib_start, ib_end)
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbnd
    INTEGER, INTENT(OUT) :: ib_start, ib_end

    INTEGER :: rest, nbnd_per_bgrp

    rest = mod ( nbnd, nbgrp )
    nbnd_per_bgrp = int( nbnd / nbgrp ) 

    IF (rest > my_bgrp_id) THEN 
       ib_start =  my_bgrp_id    * (nbnd_per_bgrp+1) + 1
       ib_end   = (my_bgrp_id+1) * (nbnd_per_bgrp+1) 
    ELSE
       ib_start =  my_bgrp_id    * nbnd_per_bgrp + rest + 1
       ib_end   = (my_bgrp_id+1) * nbnd_per_bgrp + rest 
    ENDIF

  END SUBROUTINE set_bgrp_indices


END MODULE mp_bands_util
