!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine rgen (dtau_in, rmax, mxr, at, bg, r, r2, nrm)
  !-----------------------------------------------------------------------
  !
  !   generates neighbours shells (in units of alat) with length
  !   less than rmax,and returns them in order of increasing length.
  !      r=i*a1+j*a2+k*a3-dtau,
  !   where a1,a2,a3 are the vectors defining the lattice
#include "f_defs.h"
  USE kinds
  implicit none
  !
  !    first the dummy variables
  !
  integer :: nrm, mxr
  ! output: the number of vectors in the spher
  ! input: the maximum number of vectors
  real(DP) :: r (3, mxr), r2 (mxr), at (3, 3), bg (3, 3), dtau_in (3), &
              rmax
  ! output: coordinates of vectors R+tau_s-tau
  ! output: square modulus of vectors R+tau_s-
  ! input: direct lattice vectors
  ! input: reciprocal lattice vectors
  ! input: the vector tau_s-tau_s'
  ! input: the radius of the sphere in real sp
  !
  !    and here the local variables
  !
  integer, allocatable :: irr (:)
  integer ::  nm1, nm2, nm3, i, j, k, ipol, ir, indsw, &
       iswap
  ! index on R vectors for order
  !
  !  maximum values for trial vectors
  !
  !
  !  counters on trial vectors
  !
  ! counter on polarizations
  ! counter on R vectors
  ! index of swapping
  ! used for swapping

  real(DP) :: ds(3), dtau(3)
  real(DP) :: t (3), tt, swap, DNRM2
  ! buffer contains the actual r
  ! buffer cotains the modulus of actual r
  ! used for swapping
  ! function to find the norm of a vector
  external DNRM2

  nrm = 0
  if (rmax.eq.0.d0) return

  ! convert dtau_in to the equivalent vector closest to the origin.
  ds(:) = MATMUL( dtau_in(:), bg(:,:) )
  ds(:) = ds(:) - anint(ds(:))
  dtau(:) = MATMUL( at(:,:), ds(:) )

  allocate (irr( mxr))    
  nm1 = int (DNRM2 (3, bg (1, 1), 1) * rmax) + 2
  nm2 = int (DNRM2 (3, bg (1, 2), 1) * rmax) + 2
  nm3 = int (DNRM2 (3, bg (1, 3), 1) * rmax) + 2
  !
  do i = - nm1, nm1
     do j = - nm2, nm2
        do k = - nm3, nm3
           tt = 0.d0
           do ipol = 1, 3
              t (ipol) = i * at (ipol, 1) + j * at (ipol, 2) + k * at (ipol, 3) &
                   - dtau (ipol)
              tt = tt + t (ipol) * t (ipol)
           enddo
           if (tt.le.rmax**2.and.abs (tt) .gt.1.d-10) then
              nrm = nrm + 1
              if (nrm.gt.mxr) call errore ('rgen', 'too many r-vectors', nrm)
              do ipol = 1, 3
                 r (ipol, nrm) = t (ipol)
              enddo
              r2 (nrm) = tt
           endif
        enddo
     enddo
  enddo
  !
  !   reorder the vectors in order of increasing magnitude
  !
  !   initialize the index inside sorting routine
  !
  irr (1) = 0
  if (nrm.gt.1) call hpsort (nrm, r2, irr)
  do ir = 1, nrm - 1
20   indsw = irr (ir)
     if (indsw.ne.ir) then
        do ipol = 1, 3
           swap = r (ipol, indsw)
           r (ipol, indsw) = r (ipol, irr (indsw) )
           r (ipol, irr (indsw) ) = swap
        enddo
        iswap = irr (ir)
        irr (ir) = irr (indsw)
        irr (indsw) = iswap
        goto 20
     endif

  enddo

  deallocate(irr)
  return
end subroutine rgen

