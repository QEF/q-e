!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine incdrhous (drhoscf, weight, ik, dbecsum, evcr, wgg, becq, &
     alpq, mode)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due
  !     to the displacement of the augmentation charge. Only the
  !     smooth part is computed here.
  !
#include "f_defs.h"
  USE ions_base, ONLY : ntyp => nsp, nat, ityp
  use pwcom
  USE kinds, only : DP
  USE uspp_param, ONLY: nhm, nh
  use phcom
  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum

  implicit none

  integer :: ik, mode
  ! input: the k point
  ! input: the mode which is computed

  real(DP) :: weight, wgg (nbnd, nbnd, nksq)
  ! input: the weight of the k point
  ! input: the weights

  complex(DP) :: evcr (nrxxs, nbnd), drhoscf (nrxxs), &
       dbecsum(nhm * (nhm + 1) / 2, nat), becq (nkb, nbnd, nksq), &
       alpq (nkb, nbnd, 3, nksq)
  ! input: the wavefunctions at k in real
  ! output: the change of the charge densi
  ! inp/out: the accumulated dbec
  ! input: the becp with psi_{k+q}
  ! input: the alphap with psi_{k+
  !
  !   here the local variable
  !
  real(DP) :: wgt
  ! the effective weight of the k point

  complex(DP), allocatable :: ps1 (:,:), dpsir (:)
  ! auxiliary space
  ! the change of wavefunctions in real sp

  integer :: ibnd, jbnd, nt, na, mu, ih, jh, ikb, jkb, ijkb0, &
       startb, lastb, ipol, ikk, ir, ig
  ! counters

  call start_clock ('incdrhous')
  allocate (dpsir( nrxxs))    
  allocate (ps1  ( nbnd , nbnd))    

  call divide (nbnd, startb, lastb)
  ps1 (:,:) = (0.d0, 0.d0)
  if (lgamma) then
     ikk = ik
  else
     ikk = 2 * ik - 1
  endif
  !
  !   Here we prepare the two terms
  !
  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) == nt) then
           mu = 3 * (na - 1)
           if (abs(u(mu+1,mode)) + abs(u(mu+2,mode)) &
                                 + abs(u(mu+3,mode)) > 1.0d-12) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ibnd = 1, nbnd
                       do jbnd = startb, lastb
                          do ipol = 1, 3
                             mu = 3 * (na - 1) + ipol
                             ps1(ibnd,jbnd) = ps1(ibnd,jbnd) - qq(ih,jh,nt) * &
                      ( alphap(ikb,ibnd,ipol,ik) * CONJG(becq(jkb,jbnd,ik)) + &
                        becp1(ikb,ibnd,ik) * CONJG(alpq(jkb,jbnd,ipol,ik)) ) * &
                        wgg (ibnd, jbnd, ik) * u (mu, mode)
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
           endif
           ijkb0 = ijkb0 + nh (nt)
        endif
     enddo
  enddo
#ifdef __PARA
  call mp_sum (ps1,intra_pool_comm)
#endif
  dpsi (:,:) = (0.d0, 0.d0)
  wgt = 2.d0 * weight / omega
  do ibnd = 1, nbnd_occ (ikk)
     do jbnd = 1, nbnd
        call ZAXPY (npwq, ps1(ibnd,jbnd), evq(1,jbnd), 1, dpsi(1,ibnd), 1)
     enddo
     dpsir(:) = (0.d0, 0.d0)
     do ig = 1, npwq
        dpsir(nls(igkq(ig))) = dpsi (ig, ibnd)
     enddo
     call cft3s (dpsir, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
     do ir = 1, nrxxs
        drhoscf(ir) = drhoscf(ir) + wgt * dpsir(ir) * CONJG(evcr(ir,ibnd))
     enddo
  enddo

  call addusdbec (ik, wgt, dpsi, dbecsum)
  deallocate (ps1)
  deallocate (dpsir)

  call stop_clock ('incdrhous')
  return
end subroutine incdrhous
