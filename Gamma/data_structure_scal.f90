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
  use mp_global, only: intra_pool_comm, nproc_pool
  use fft_scalar, only: good_fft_dimension
  use fft_types

  implicit none
  integer :: n1, n2, n3, i1, i2, i3
  ! counters on G space
  !
  ! a function with obvious meaning

  real(kind=DP) :: amod
  real(kind=DP) :: gkcut
  ! modulus of G vectors

  integer, allocatable :: stw(:,:)
  integer :: ub(3), lb(3), i, j, k, ii, jj
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

  CALL fft_dlay_allocate( dfftp, nproc_pool, MAX(nrx1, nrx3),  nrx2  )
  CALL fft_dlay_allocate( dffts, nproc_pool, MAX(nrx1s, nrx3s) , nrx2s )

  !
  !     compute the number of g necessary to the calculation
  !
  n1 = nr1 + 1
  n2 = nr2 + 1
  n3 = nr3 + 1

  ub = (/  n1,  n2,  n3 /)
  lb = (/ -n1, -n2, -n3 /)

  allocate( stw( lb(2) : ub(2), lb(3) : ub(3) ) )
  stw = 0

  gkcut = ecutwfc / tpiba2

  ngm  = 0
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
           if ( amod <= gkcut ) then
             stw(  i2,  i3 ) = 1
           end if
20         continue
        enddo
10      continue
     enddo
  enddo

  call fft_dlay_scalar( dfftp, ub, lb, nr1, nr2, nr3, nrx1, nrx2, nrx3, stw )
  call fft_dlay_scalar( dffts, ub, lb, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, stw )

  deallocate( stw )

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

