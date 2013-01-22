!
! Copyright (C) 2002-2004 FPMD & PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE io_global
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  PUBLIC :: io_global_start, io_global_getionode
  PUBLIC :: stdout, ionode, ionode_id, meta_ionode, meta_ionode_id
  PUBLIC :: xmlinputunit, xmloutputunit, xmltmpunit
  !
  INTEGER :: stdout = 6            ! unit connected to standard output
  ! For parallel execution: I/O within an image
  INTEGER :: ionode_id = 0         ! index of the i/o node for this image
  LOGICAL :: ionode = .TRUE.       ! true if this processor is a i/o node
                                   ! for this image 
  ! For parallel execution: global I/O node (for NEB, PHonon, etc)
  INTEGER :: meta_ionode_id = 0    ! index of the global i/o node
  LOGICAL :: meta_ionode = .TRUE.  ! true if this processor is global i/o node
  LOGICAL :: first = .TRUE.
  !
  INTEGER :: xmlinputunit          ! unit connected to the xml input
  INTEGER :: xmloutputunit = 51    ! unit connected to the xml output
  INTEGER :: xmltmpunit    = 52    ! unit connected to the temp xml output
  !    
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE io_global_start( mpime, ionode_set )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: mpime, ionode_set
       !
       !
       IF ( mpime == ionode_set ) THEN
          ionode      = .TRUE.
       ELSE
          ionode      = .FALSE.
       END IF
       ionode_id      = ionode_set
       !
       first = .FALSE.
       !
       RETURN
       !
     END SUBROUTINE io_global_start
     !
     !-----------------------------------------------------------------------
     SUBROUTINE io_global_getionode( ionode_out, ionode_id_out )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       LOGICAL, INTENT(OUT) :: ionode_out
       INTEGER, INTENT(OUT) :: ionode_id_out
       !
       !
       IF ( first ) &
          CALL errore( ' io_global_getionode ', ' ionode not yet defined ', 1 )
       !
       ionode_out    = ionode
       ionode_id_out = ionode_id
       !
       RETURN
       !
     END SUBROUTINE io_global_getionode
     !  
     !
END MODULE io_global
