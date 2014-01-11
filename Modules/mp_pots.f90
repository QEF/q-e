!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE mp_pots
  !----------------------------------------------------------------------------
  !
  USE mp, ONLY : mp_barrier, mp_size, mp_rank, mp_comm_split
  USE parallel_include
  !
  IMPLICIT NONE 
  SAVE
  !
  ! ... Pot groups (processors within a cooking-pot)
  ! ... Used only in a specialized calculation under development 
  !
  INTEGER :: npot       = 1  ! number of pots
  INTEGER :: nproc_pot  = 1  ! number of processors within a pot
  INTEGER :: me_pot     = 0  ! index of the processor within a pot
  INTEGER :: root_pot   = 0  ! index of the root processor within a pot
  INTEGER :: my_pot_id  = 0  ! index of my pot
  INTEGER :: inter_pot_comm  = 0  ! inter pot communicator
  INTEGER :: intra_pot_comm  = 0  ! intra pot communicator
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE mp_start_pots ( npot_, parent_comm )
    !---------------------------------------------------------------------------
    !
    ! ... Divide processors (of the "parent_comm" group) into "pots"
    ! ... Requires: npot_, read from command line
    ! ...           parent_comm, typically processors of an "image"
    ! ...           (intra_image_comm)
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: npot_, parent_comm
    !
    INTEGER :: parent_nproc = 1, parent_mype  = 0
    !
#if defined (__MPI)
    !
    parent_nproc = mp_size( parent_comm )
    parent_mype  = mp_rank( parent_comm )
    !
    ! ... npot_ must have been previously read from command line argument
    ! ... by a call to routine get_command_line
    !
    npot = npot_
    !
    IF ( npot < 1 .OR. npot > parent_nproc ) CALL errore( 'mp_start_pots',&
                          'invalid number of pot groups, out of range', 1 )

    IF ( MOD( parent_nproc, npot ) /= 0 ) CALL errore( 'mp_start_pots', &
           'invalid number of pots, parent_nproc /= nproc_pot * npot', 1 )  
    !
    ! ... number of cpus per pot (they are created inside each parent group)
    !
    nproc_pot = parent_nproc / npot
    !
    !
    ! ... my_pot_id  =  pot index for this processor    ( 0 : npot - 1 )
    ! ... me_pot     =  processor index within the pot  ( 0 : nproc_pot - 1 )
    !
    my_pot_id = parent_mype / nproc_pot    
    me_pot    = MOD( parent_mype, nproc_pot )
    !
    CALL mp_barrier( parent_comm )
    !
    ! ... the intra_pot_comm communicator is created
    !
    CALL mp_comm_split( parent_comm, my_pot_id, parent_mype, intra_pot_comm )
    !
    CALL mp_barrier( parent_comm )
    !
    ! ... the inter_pot_comm communicator is created
    !
    CALL mp_comm_split( parent_comm, me_pot, parent_mype, inter_pot_comm )
    !
#endif
    !
    RETURN
  END SUBROUTINE mp_start_pots
  !
END MODULE mp_pots
