!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine drho_drc (iudrho_x, u_x, xq_x, drc_x, scale)
  !-----------------------------------------------------------------------
  !  Reads the variation of the charge saved on a file and changes
  !  it according to the variation of the core_charge
  !  It is used by drho_cc. Have a look there for more explanation
  !
#include "machine.h"
  use pwcom
  use phcom
  use d3com
#ifdef PARA
  use para
#endif
  implicit none
#ifdef PARA
  include 'mpif.h'
#endif

  integer :: iudrho_x
  !input: the unit containing the charge variation
  real (8) :: xq_x (3), scale
  !input: q point
  !input: drhocore will be added to the valence charge scaled by this factor
  complex (8) :: u_x (3 * nat, 3 * nat), drc_x (ngm, ntyp)
  !input: the transformation modes patterns
  !input: contain the rhoc (without structu

  integer :: ipert, na, mu, nt, ig, errcode
  real (8) :: gtau
  complex (8) :: guexp
  complex (8), allocatable :: drhoc (:), drhov (:), uact (:)


  allocate  (drhoc( nrxx))    
  allocate  (drhov( nrxx))    
  allocate  (uact( 3 * nat))    
#ifdef PARA
!  if (mypool.ne.1) goto 100
#endif

  do ipert = 1, 3 * nat
     call setv (2 * nrxx, 0.d0, drhoc, 1)

     call ZCOPY (3 * nat, u_x (1, ipert), 1, uact, 1)
     do na = 1, nat
        mu = 3 * (na - 1)
        if (abs (uact (mu + 1) ) + abs (uact (mu + 2) ) + abs (uact (mu + &
             3) ) .gt.1.0d-12) then
           nt = ityp (na)
           if (nlcc (nt) ) then
              do ig = 1, ngm
                 gtau = tpi * ( (g (1, ig) + xq_x (1) ) * tau (1, na) &
                      + (g (2, ig) + xq_x (2) ) * tau (2, na) + (g (3, ig) &
                      + xq_x (3) ) * tau (3, na) )
                 guexp = tpiba * ( (g (1, ig) + xq_x (1) ) * uact (mu + 1) &
                      + (g (2, ig) + xq_x (2) ) * uact (mu + 2) + (g (3, ig) &
                      + xq_x (3) ) * uact (mu + 3) ) * DCMPLX (0.d0, - 1.d0) &
                      * DCMPLX (cos (gtau), - sin (gtau) )
                 drhoc (nl (ig) ) = drhoc (nl (ig) ) + drc_x (ig, nt) &
                      * guexp
              enddo
           endif
        endif

     enddo
     call cft3 (drhoc, nr1, nr2, nr3, nrx1, nrx2, nrx3, + 1)
     call davcio_drho2 (drhov, lrdrho, iudrho_x, ipert, - 1)
     call DAXPY (2 * nrxx, scale, drhoc, 1, drhov, 1)

     call davcio_drho2 (drhov, lrdrho, iudrho_x, ipert, + 1)
  enddo
#ifdef PARA
100 continue
  call MPI_barrier (MPI_COMM_WORLD, errcode)

  call error ('drho_drc', 'at barrier', errcode)
#endif
  deallocate (drhoc)
  deallocate (drhov)
  deallocate (uact)
  return
end subroutine drho_drc
