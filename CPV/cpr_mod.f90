!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

module stre
  implicit none 
  save
  real(kind=8) stress(3,3)
end module stre

module dqrad_mod
  implicit none 
  save
  real(kind=8),allocatable:: dqrad(:,:,:,:,:,:,:)
end module dqrad_mod

module betax
  implicit none 
  save
  integer, parameter:: mmx=5001
  real(kind=8) :: refg
  real(kind=8),allocatable:: betagx(:,:,:), dbetagx(:,:,:), &
                       qradx(:,:,:,:,:), dqradx(:,:,:,:,:)
end module betax


