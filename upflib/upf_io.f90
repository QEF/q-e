!
! Copyright (C) 2002-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE upf_io
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  PUBLIC :: stdin, stdout
  !
  INTEGER, PARAMETER :: stdin  = 5    ! unit connected to standard input
  INTEGER            :: stdout = 6    ! unit connected to standard output
  !
  !
END MODULE upf_io

