!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine d3_exc  
  !-----------------------------------------------------------------------
  !
  !    Calculates the contribution to the derivative of the dynamical
  !    matrix due to the third derivative of the exchange and correlation
  !    energy
#include "machine.h"
  use pwcom
  use phcom
  use d3com
  use allocate
#ifdef PARA
  use para
#endif
  implicit none
  integer :: errcode, ir, ipert, jpert, kpert, npert1, npert2  
  real (8) :: d2mxc, rhotot, xq0 (3)
  real (8), pointer :: d2muxc (:)
  complex (8) :: aux
  complex (8), pointer :: work1 (:), work2 (:), work3 (:), d3dyn1 (:,:,:)

  call mallocate(d2muxc, nrxx)  
  call mallocate(work1 , nrxx)  
  call mallocate(work2 , nrxx)  
  call mallocate(work3 , nrxx)  
  call mallocate(d3dyn1, 3*nat, 3*nat, 3*nat)  
#ifdef PARA
  if (mypool.ne.1) goto 100  
#endif
  !
  ! Calculates third derivative of Exc
  !
  call setv (nrxx, 0.d0, d2muxc, 1)  
  do ir = 1, nrxx  
     rhotot = rho (ir, 1) + rho_core (ir)  
     if (rhotot.gt.1.d-30) d2muxc (ir) = d2mxc (rhotot)  
     if (rhotot.lt. - 1.d-30) d2muxc (ir) = - d2mxc ( - rhotot)  
  enddo
  !
  ! Calculates the contribution to d3dyn
  !
  call setv (2 * 27 * nat * nat * nat, 0.d0, d3dyn1, 1)  
  do ipert = 1, 3 * nat  
     if (q0mode (ipert) ) then  
        call davcio_drho (work1, lrdrho, iud0rho, ipert, - 1)  
        do jpert = 1, 3 * nat  
           call davcio_drho (work2, lrdrho, iudrho, jpert, - 1)  
           do kpert = 1, 3 * nat  
              call davcio_drho (work3, lrdrho, iudrho, kpert, - 1)  
              aux = DCMPLX (0.d0, 0.d0)  
              do ir = 1, nrxx  
                 aux = aux + d2muxc (ir) * work1 (ir) * conjg (work2 (ir) ) &
                      * work3 (ir)
              enddo
#ifdef PARA
              call reduce (2, aux)  
#endif
              d3dyn1 (ipert, jpert, kpert) = omega * aux / (nr1 * nr2 * nr3)  
           enddo
        enddo
     endif

  enddo
#ifdef PARA
100 continue  
  call poolbcast (2 * 27 * nat * nat * nat, d3dyn1)  
#endif

  call DAXPY (2 * 27 * nat * nat * nat, 1.d0, d3dyn1, 1, d3dyn, 1)  
  call ZCOPY (27 * nat * nat * nat, d3dyn1, 1, d3dyn_aux9, 1)

  call mfree (d2muxc)  
  call mfree (work1)  
  call mfree (work2)  
  call mfree (work3)  
  call mfree (d3dyn1) 
 
  return  
end subroutine d3_exc
