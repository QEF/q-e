!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

subroutine work_function (wf)
#include "machine.h"
  !
  ! Print out the workfunction, calculated as the difference between the
  ! potential energy and the fermi energy.
  ! Written for supercells with the main axis along z.
  !
  use pwcom
#ifdef PARA
  use para
#endif
  implicit none
#ifdef PARA
  include 'mpif.h'
#endif

  integer :: ierr

  real(kind=DP) :: wmean1, wmean2, meancharge, wx1, wx2, wxm, vx, vc, ex, &
       ec, rhox, rs, vcca, wf

  integer :: n1, n2, ni, nmean

  logical :: exst
  real(kind=DP), allocatable :: raux1 (:), vaux1 (:), aux (:)
  ! auxiliary vectors for charge and potential


  allocate (raux1( nrx1 * nrx2 * nrx3))    
  allocate (vaux1( nrx1 * nrx2 * nrx3))    
  allocate (aux  ( nrxx))    

  if (.not.lscf) call sum_band
  call DCOPY (nrxx, rho, 1, aux, 1)

  call DAXPY (nrxx, 1.d0, rho_core, 1, aux, 1)
#ifdef PARA
  call gather (aux, raux1)
#else
  call DCOPY (nrxx, aux, 1, raux1, 1)
#endif
  call DCOPY (nrxx, vltot, 1, aux, 1)

  call DAXPY (nrxx, 1.d0, vr, 1, aux, 1)
#ifdef PARA
  call gather (aux, vaux1)
#else
  call DCOPY (nrxx, aux, 1, vaux1, 1)
#endif
#ifdef PARA
  if (me.eq.1.and.mypool.eq.1) then
#endif
     call seqopn (17, 'workf', 'formatted', exst)
     call seqopn (19, 'charge', 'formatted', exst)
     nmean = (nr3 + 1) / 2
     do nmean = 1, nr3
        wmean1 = 0.d0
        wmean2 = 0.d0
        meancharge = 0.d0
        wx1 = 0.d0
        wx2 = 0.d0
        wxm = 0.d0
        do n2 = 1, nr2
           do n1 = 1, nr1
              ni = n1 + (n2 - 1) * nrx1 + (nmean - 1) * nrx1 * nrx2
              meancharge = meancharge+raux1 (ni)
              wxm = wxm + raux1 (ni) **2
              wmean1 = wmean1 + vaux1 (ni)
              wx1 = wx1 + vaux1 (ni) **2
              rhox = abs (raux1 (ni) )
              vx = 0.d0
              if (rhox.gt.1.0e-30) then
                 call xc (rhox, ex, ec, vx, vc)
                 vx = e2 * (vx + vc)
              endif
              wmean2 = wmean2 + vaux1 (ni) - vx
              wx2 = wx2 + (vaux1 (ni) - vx) **2
           enddo
        enddo
        wmean1 = wmean1 / dfloat (nr1 * nr2)
        wmean2 = wmean2 / dfloat (nr1 * nr2)
        meancharge = meancharge / dfloat (nr1 * nr2)
        wx1 = dsqrt (wx1 / dfloat (nr1 * nr2) - wmean1 * wmean1)
        wx2 = dsqrt (wx2 / dfloat (nr1 * nr2) - wmean2 * wmean2)
        wxm = dsqrt (wxm / dfloat (nr1 * nr2) - meancharge**2)
        if (nmean.eq. (nr3 + 1) / 2) wf = wmean2 - ef
        write (17, * ) nmean, (wmean1 - ef) * rytoev, wx1 * rytoev, &
             (wmean2 - ef) * rytoev, wx2 * rytoev
        write (19, * ) nmean, meancharge, wxm
        if (nmean.eq. (nr3 + 1) / 2) then
           write (6, 9130) rytoev * (wmean1 - ef), wx1 * rytoev, &
                rytoev * (wmean2 - ef), wx2 * rytoev
        endif
     enddo
#ifdef PARA
  endif

  call mpi_bcast (wf, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
#endif
  write (6, '(/5x,"Work function written on file workf")')

    write (6, '(5x,"Planar mean charge written on file charge")')

9130 format (/'     workfunction     = ',f10.4,' +- ',f6.4,' eV', &
         &        /'     without exchcorr = ',f10.4,' +- ',f6.4,' eV')
    close (17)

    close (19)
    deallocate(raux1)
    deallocate(vaux1)
    deallocate(aux)
    return

  end subroutine work_function
