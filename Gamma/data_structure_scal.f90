!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine data_structure_scal
  !-----------------------------------------------------------------------
  ! this routine sets the data structure for the fft arrays.
  ! This version computes also the smooth and hard mesh
  !
#include "machine.h"
  use pwcom
  use mp, only: mp_sum
  use mp_global, only: intra_pool_comm
  implicit none
  integer :: n1, n2, n3, i1, i2, i3
  ! counters on G space
  !
  integer :: good_fft_dimension
  ! a function with obvious meaning

  real(kind=DP) :: amod
  ! modulus of G vectors
  !
  nrx1 = good_fft_dimension (nr1)
  nrx1s = good_fft_dimension (nr1s)
  !
  !     nrx2 and nrx3 are there just for compatibility
  !
  nrx2 = nr2
  nrx3 = nr3

  nrxx = nrx1 * nrx2 * nrx3
  nrx2s = nr2s
  nrx3s = nr3s
  nrxxs = nrx1s * nrx2s * nrx3s
  !
  !     compute the number of g necessary to the calculation
  !
  n1 = nr1 + 1
  n2 = nr2 + 1
  n3 = nr3 + 1

  ngm = 0
  ngms = 0
  !
  !     exclude space with x<0
  !
  do i1 = 0, n1
     do i2 = - n2, n2
        !
        !     exclude plane with x=0, y<0
        !
        if(i1.eq.0.and.i2.lt.0) go to 10
        do i3 = - n3, n3
           !
           !     exclude line with x=0, y=0, z<0
           !
           if(i1.eq.0.and.i2.eq.0.and.i3.lt.0) go to 20
           amod = (i1 * bg (1, 1) + i2 * bg (1, 2) + i3 * bg (1, 3) ) **2 + &
                  (i1 * bg (2, 1) + i2 * bg (2, 2) + i3 * bg (2, 3) ) **2 + &
                  (i1 * bg (3, 1) + i2 * bg (3, 2) + i3 * bg (3, 3) ) **2
           if (amod.le.gcutm) ngm = ngm + 1
           if (amod.le.gcutms) ngms = ngms + 1
20         continue
        enddo
10      continue
     enddo
  enddo
  !
  !     compute the global number of g, i.e. the sum over all processors
  !     whithin a pool
  !
  ngm_l  = ngm
  ngms_l = ngms
  ngm_g  = ngm
  ngms_g = ngms
  call mp_sum( ngm_g , intra_pool_comm )
  call mp_sum( ngms_g, intra_pool_comm )
  !
  return
end subroutine data_structure_scal

