!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
!-----------------------------------------------------------------------
subroutine d3_exc
  !-----------------------------------------------------------------------
  !
  !    Calculates the contribution to the derivative of the dynamical
  !    matrix due to the third derivative of the exchange and correlation
  !    energy
  !
  USE ions_base,  ONLY : nat
  USE kinds, ONLY : DP
  use pwcom
  use phcom
  use d3com
  use para
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm, my_image_id
  USE mp,        ONLY : mp_bcast  

  implicit none
  
  integer :: errcode, ir, ipert, jpert, kpert, npert1, npert2
  real (kind = dp) :: d2mxc, rhotot, xq0 (3)
  real (kind = dp), allocatable :: d2muxc (:)
  complex (kind = dp) :: aux
  complex (kind = dp), allocatable :: work1 (:), work2 (:), work3 (:), d3dyn1 (:,:,:)

  allocate (d2muxc( nrxx))    
  allocate (work1 ( nrxx))    
  allocate (work2 ( nrxx))    
  allocate (work3 ( nrxx))    
  allocate (d3dyn1( 3*nat, 3*nat, 3*nat))    
#ifdef __PARA
  if (mypool.ne.1) goto 100
#endif
  !
  ! Calculates third derivative of Exc
  !
  d2muxc(:) = 0.d0
  do ir = 1, nrxx
     rhotot = rho (ir, 1) + rho_core (ir)
     if (rhotot > 1.d-30) d2muxc (ir) = d2mxc (rhotot)
     if (rhotot < - 1.d-30) d2muxc (ir) = - d2mxc ( - rhotot)
  enddo
  !
  ! Calculates the contribution to d3dyn
  !
  d3dyn1 (:,:,:) = (0.d0, 0.d0)
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
#ifdef __PARA
              call reduce (2, aux)
#endif
              d3dyn1 (ipert, jpert, kpert) = omega * aux / (nr1 * nr2 * nr3)
           enddo
        enddo
     endif

  enddo
#ifdef __PARA
100 continue  
  IF ( npool /= 1 ) &
     CALL mp_bcast( d3dyn1, ionode_id, inter_pool_comm )
#endif

  d3dyn = d3dyn  + d3dyn1
  d3dyn_aux9 = d3dyn1

  deallocate (d2muxc)
  deallocate (work1)
  deallocate (work2)
  deallocate (work3)
  deallocate (d3dyn1)

  return
end subroutine d3_exc
