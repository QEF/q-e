!
! Copyright (C) 2002-2013 Quantum ESPRESSO group
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
  PUBLIC :: stdin, stdout, qestdin
  PUBLIC :: ionode, ionode_id, meta_ionode, meta_ionode_id
  !
  INTEGER, PARAMETER :: stdin  = 5    ! unit connected to standard input
  INTEGER :: qestdin= 9    ! unit connected to input file (xml or text)
  INTEGER :: stdout = 6    ! unit connected to standard output
  !
  ! For parallel execution: I/O within an image
  ! These are set at startup by calling mp_world_start
  !
  INTEGER :: ionode_id = 0         ! index of the i/o node for this image
  LOGICAL :: ionode = .TRUE.       ! true if this processor is a i/o node
                                   ! for this image 
  ! For parallel execution: global I/O node (for NEB, PHonon, etc)
  ! These are set at startup by calling mp_image_start
  !
  INTEGER :: meta_ionode_id = 0    ! index of the global i/o node
  LOGICAL :: meta_ionode = .TRUE.  ! true if this processor is global i/o node
  !
  INTEGER :: xmloutputunit = 51    ! unit connected to the xml output
  !    
END MODULE io_global
