!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvqpsi_us_only (ik, mode, uact)
  !----------------------------------------------------------------------
  !
  ! This routine calculates dV_bare/dtau * psi for one perturbation
  ! with a given q. The displacements are described by a vector uact.
  ! The result is stored in dvpsi. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  !
#include "machine.h"

  use pwcom
  USE kinds, only : DP
  use phcom
  implicit none
  !
  !   The dummy variables
  !

  integer :: ik, mode
  ! input: the k point
  ! input: the actual perturbation
  complex(kind=DP) :: uact (3 * nat)
  ! input: the pattern of displacements
  !
  !   And the local variables
  !

  integer :: na, nb, mu, nu, ikk, ikq, ig, igg, nt, ibnd, ijkb0, &
       ikb, jkb, ih, jh, ipol
  ! counter on atoms
  ! counter on modes
  ! the point k
  ! the point k+q
  ! counter on G vectors
  ! auxiliary counter on G vectors
  ! counter on atomic types
  ! counter on bands
  ! auxiliary variable for counting
  ! counter on becp functions
  ! counter on becp functions
  ! counter on n index
  ! counter on m index
  ! counter on polarizations

  real(kind=DP), parameter :: eps = 1.d-12

  complex(kind=DP), allocatable :: ps1 (:,:), ps2 (:,:,:), aux (:)
  ! work space

  logical :: ok

  call start_clock ('dvqpsi_us_on')
  allocate (ps1 ( nkb , nbnd))    
  allocate (ps2 ( nkb , nbnd , 3))    
  allocate (aux ( npwx))    
  if (lgamma) then
     ikk = ik
     ikq = ik
  else
     ikk = 2 * ik - 1
     ikq = ikk + 1
  endif
  if (lsda) current_spin = isk (ikk)
  !
  !   we first compute the coefficients of the vectors
  !
  ps1(:,:)   = (0.d0, 0.d0)
  ps2(:,:,:) = (0.d0, 0.d0)
  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
           mu = 3 * (na - 1)
           do ih = 1, nh (nt)
              ikb = ijkb0 + ih
              do jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 do ipol = 1, 3
                    do ibnd = 1, nbnd
                       if ( abs (uact (mu + 1) ) + &
                            abs (uact (mu + 2) ) + &
                            abs (uact (mu + 3) ) > eps) then
                          ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                               (deeq (ih, jh, na, current_spin) - &
                                et (ibnd, ikk) * qq (ih, jh, nt) ) * &
                                alphap(jkb, ibnd, ipol, ik) * uact (mu + ipol)
                          ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) +&
                               (deeq (ih,jh, na, current_spin) - &
                                et (ibnd, ikk) * qq (ih, jh, nt) ) * &
                                (0.d0, -1.d0) * becp1 (jkb, ibnd, ik) * &
                                uact (mu + ipol) * tpiba
                          if (okvan) then
                             ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                                  (int1 (ih, jh, ipol,na, current_spin) * &
                                  becp1 (jkb, ibnd, ik) ) * uact (mu +ipol)
                          endif
                       endif
                       if (okvan) then
                          do nb = 1, nat
                             nu = 3 * (nb - 1)
                             ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                                  (int2 (ih, jh, ipol, nb, na) * &
                                   becp1 (jkb, ibnd, ik) ) * uact (nu + ipol)
                          enddo
                       endif
                    enddo
                 enddo
              enddo
           enddo
           ijkb0 = ijkb0 + nh (nt)
        endif
     enddo
  enddo
  !
  !      This term is proportional to beta(k+q+G)
  !
  if (nkb.gt.0) call ZGEMM ('N', 'N', npwq, nbnd, nkb, &
       (1.d0, 0.d0) , vkb, npwx, ps1, nkb, (1.d0, 0.d0) , dvpsi, npwx)
  !
  !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
  !
  do ikb = 1, nkb
     do ipol = 1, 3
        ok = .false.
        do ibnd = 1, nbnd
           ok = ok.or. (abs (ps2 (ikb, ibnd, ipol) ) .gt.eps)
        enddo
        if (ok) then
           do ig = 1, npwq
              igg = igkq (ig)
              aux (ig) =  vkb(ig, ikb) * (xk(ipol, ikq) + g(ipol, igg) )
           enddo
           do ibnd = 1, nbnd
              call ZAXPY (npwq, ps2(ikb,ibnd,ipol), aux, 1, dvpsi(1,ibnd), 1)
           enddo
        endif
     enddo

  enddo
  deallocate (aux)
  deallocate (ps2)
  deallocate (ps1)

  call stop_clock ('dvqpsi_us_on')
  return
end subroutine dvqpsi_us_only
