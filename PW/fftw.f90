!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module fftw_mod
  integer :: FFTW_FORWARD, FFTW_BACKWARD
  parameter (FFTW_FORWARD = - 1, FFTW_BACKWARD = 1)
  integer :: FFTW_ESTIMATE, FFTW_MEASURE
  parameter (FFTW_ESTIMATE = 0, FFTW_MEASURE = 1)
  integer :: FFTW_OUT_OF_PLACE, FFTW_IN_PLACE, FFTW_USE_WISDOM
  parameter (FFTW_OUT_OF_PLACE = 0)
  parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)
end module fftw_mod
