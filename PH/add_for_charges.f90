!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine add_for_charges (ik, uact)
!----------===============-----------------------------------------------
  !
  ! This subroutine calculates dS/du P_c [x, H-eS] |psi>
  !
#include "f_defs.h"

  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  use pwcom
  USE kinds, only : DP
  USE uspp_param, only: nh
  use phcom
  implicit none
  !
  !   The dummy variables
  !

  integer :: ik, mode
  ! input: the k point
  ! input: the actual perturbation
  complex(DP) :: uact (3 * nat)
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

  real(DP), parameter :: eps = 1.d-12

  complex(DP), allocatable :: ps1 (:,:), ps2 (:,:,:), aux (:)
  ! the scalar product
  ! the scalar product
  ! a mesh space for psi
  complex(DP), allocatable :: bedp(:,:), aux1(:,:), alphapp(:,:,:)

  logical :: ok
  ! used to save time

  allocate (ps1 ( nkb , nbnd))    
  allocate (ps2 ( nkb , nbnd , 3))    
  allocate (aux ( npwx))
  allocate (aux1( npwx, nbnd))
  allocate (bedp( nkb, nbnd) )
  allocate (alphapp (nkb,nbnd,3))
  if (lgamma) then
     ikk = ik
     ikq = ik
  else
     call infomsg ('add_for_charges', 'called for lgamma .eq. false', -1)
  endif
  if (lsda) current_spin = isk (ikk)
  !
  !   we first compute the coefficients of the vectors
  !
  ps1   = (0.d0, 0.d0)
  ps2   = (0.d0, 0.d0)
  aux1  = (0.d0, 0.d0)
  alphapp = (0.d0,0.d0)
  bedp = (0.d0,0.d0)

  !
  ! first we calculate the products of the beta functions with dvpsi 
  !
  call ccalbec (nkb, npwx, npw, nbnd, bedp, vkb, dpsi)
  do ipol = 1, 3
     do ibnd = 1, nbnd
        do ig = 1, npw
           aux1 (ig, ibnd) = dpsi(ig,ibnd) *           &
                tpiba * (0.d0,1.d0) *                  & 
                ( xk(ipol,ikk) + g(ipol,igk(ig)) )
        enddo
     enddo
     call ccalbec (nkb, npwx, npw, nbnd, alphapp(1,1,ipol), vkb, aux1)
  enddo


  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
           mu = 3 * (na - 1)
           if ( abs (uact (mu + 1) ) + &
                abs (uact (mu + 2) ) + &
                abs (uact (mu + 3) ) > eps) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    do ipol = 1, 3
                       do ibnd = 1, nbnd
                          ps1 (ikb, ibnd) = ps1 (ikb, ibnd) +     &
                               qq (ih, jh, nt) *                  &
                               alphapp(jkb, ibnd, ipol) *         &
                               uact (mu + ipol)
                          ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) + &
                               qq (ih, jh, nt) *                          &
                                (0.d0, -1.d0) *                           &
                                bedp (jkb, ibnd) *                        &
                                uact (mu + ipol) * tpiba
                       enddo
                    enddo
                 enddo
              enddo
           endif
           ijkb0 = ijkb0 + nh (nt)
        endif
     enddo
  enddo
  !
  !      This term is proportional to beta(k+q+G)
  !
  dvpsi = matmul(vkb,ps1) + dvpsi
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
              dvpsi(:,ibnd) = ps2(ikb,ibnd,ipol) * aux + dvpsi(:,ibnd)
           enddo
        endif
     enddo
  enddo
!
!    Now dvpsi contains dS/du x |psi>
!

  deallocate (aux)
  deallocate (aux1)
  deallocate (ps2)
  deallocate (ps1)
  deallocate (bedp)
  deallocate (alphapp)

  return
end subroutine add_for_charges

