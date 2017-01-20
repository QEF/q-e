!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_bands
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
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE mp_start_bands( nband_, ntg_, parent_comm )
    !---------------------------------------------------------------------------
    !
    ! ... Divide processors (of the "parent_comm" group) into nband_ pools
    ! ... Requires: nband_, read from command line
    ! ...           parent_comm, typically processors of a k-point pool
    ! ...           (intra_pool_comm)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nband_, parent_comm
    INTEGER, INTENT(IN), OPTIONAL :: ntg_
    !
    INTEGER :: parent_nproc = 1, parent_mype = 0
    !
#if defined (__MPI)
    !
    parent_nproc = mp_size( parent_comm )
    parent_mype  = mp_rank( parent_comm )
    !
    ! ... nband_ must have been previously read from command line argument
    ! ... by a call to routine get_command_line
    !
    nbgrp = nband_
    !
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
    ! ... Set number of processors per band group
    !
    nproc_bgrp = parent_nproc / nbgrp
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
    ! ... the intra_bgrp_comm communicator is created
    !
    CALL mp_comm_split( parent_comm, my_bgrp_id, parent_mype, intra_bgrp_comm )
    !
    CALL mp_barrier( parent_comm )
    !
    ! ... the inter_bgrp_comm communicator is created                     
    !     
    CALL mp_comm_split( parent_comm, me_bgrp, parent_mype, inter_bgrp_comm )  
    !
    IF ( PRESENT(ntg_) ) THEN
       ntask_groups = ntg_
    END IF
    !
#endif
    RETURN
    !
  END SUBROUTINE mp_start_bands
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

  INTEGER FUNCTION bgrp_start(nbnd)
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbnd

    INTEGER :: rest, nbnd_per_bgrp

    rest = mod ( nbnd, nbgrp )
    nbnd_per_bgrp = int( nbnd / nbgrp ) 

    IF (rest > my_bgrp_id) THEN 
       bgrp_start =  my_bgrp_id    * (nbnd_per_bgrp+1) + 1
    ELSE
       bgrp_start =  my_bgrp_id    * nbnd_per_bgrp + rest + 1
    ENDIF

  END FUNCTION bgrp_start

  INTEGER FUNCTION bgrp_end(nbnd)
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbnd

    INTEGER :: rest, nbnd_per_bgrp

    rest = mod ( nbnd, nbgrp )
    nbnd_per_bgrp = int( nbnd / nbgrp ) 

    IF (rest > my_bgrp_id) THEN 
       bgrp_end   = (my_bgrp_id+1) * (nbnd_per_bgrp+1) 
    ELSE
       bgrp_end   = (my_bgrp_id+1) * nbnd_per_bgrp + rest 
    ENDIF

  END FUNCTION bgrp_end

END MODULE mp_bands
!
!     
MODULE mp_bands_TDDFPT
!
! NB: These two variables used to be in mp_bands and are loaded from mp_global in TDDFPT 
!     I think they would better stay in a TDDFPT specific module but leave them here not to
!     be too invasive on a code I don't know well. SdG
!     
  INTEGER :: ibnd_start = 0              ! starting band index used in bgrp parallelization
  INTEGER :: ibnd_end = 0                ! ending band index used in bgrp parallelization
!     
END MODULE mp_bands_TDDFPT
!     
