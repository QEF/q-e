!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine psidspsi (ik, uact, pdsp)
!----------========----------------------------------------------------
  !
  ! This routine calculates <psi_v'|ds/du|psi_v>
  ! at q=0. The displacements are described by a vector uact.
  ! The result is stored in pdsp. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  !
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE pwcom
  USE kinds, ONLY : DP
  USE wavefunctions_module,    ONLY : evc
  USE uspp_param,    ONLY : nh
  USE phcom
  implicit none
  !
  !   The dummy variables
  !

  integer, intent(in) :: ik
  ! input: the k point
  complex(DP) :: uact (3 * nat), pdsp(nbnd,nbnd)
  ! input: the pattern of displacements
  ! output: <psi|ds/du|psi>
  !
  !   And the local variables
  !

  integer :: na, nb, mu, nu, ikk, ikq, ig, igg, nt, ibnd, jbnd, ijkb0, &
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

  complex(DP), ALLOCATABLE :: ps1 (:,:), ps2 (:,:,:), aux (:), dspsi(:,:)
  ! the scalar product
  ! the scalar product
  ! a mesh space for psi
  ! the matrix dspsi

  logical :: ok
  ! used to save time

  allocate (ps1 ( nkb , nbnd ))    
  allocate (ps2 ( nkb , 3, nbnd))
  allocate (dspsi (npwx,nbnd))
  allocate (aux ( npwx))    

  if (lgamma) then
     ikk = ik
     ikq = ik
  else
     call infomsg ('psidspsi', 'called for lgamma .eq. false', -1)
  endif
  if (lsda) current_spin = isk (ikk)

  ps1(:,:)   = (0.d0, 0.d0)
  ps2(:,:,:) = (0.d0, 0.d0)
  pdsp(:,:)   = (0.d0, 0.d0)
  dspsi = (0.d0,0.d0)
  !
  !
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
                          ps1 (ikb, ibnd) = ps1 (ikb, ibnd) +    &
                               qq (ih, jh, nt) *                 &
                               alphap(jkb, ibnd, ipol, ik) *     &
                               uact (mu + ipol)
                          ps2 (ikb, ipol, ibnd) = ps2 (ikb, ipol, ibnd) + &
                               qq (ih, jh, nt) *                          &
                               (0.d0, -1.d0) *                            &
                               becp1 (jkb, ibnd, ik) *                    &
                               uact (mu + ipol) * tpiba
                       enddo
                    enddo
                 enddo
              enddo
           endif
           ijkb0= ijkb0 + nh (nt)
        endif
     enddo
  enddo
  !
  !      This term is proportional to beta(k+q+G)
  !
  dspsi = matmul(vkb,ps1)+ dspsi
  !
  !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
  !
  do ikb = 1, nkb
     do ipol = 1, 3
        ok = .false.
        do ibnd = 1, nbnd
           ok = ok.or. (abs (ps2 (ikb, ipol, ibnd) ) .gt.eps)
        enddo
        if (ok) then
           do ig = 1, npw
              igg = igk (ig)
              aux (ig) =  vkb(ig, ikb) *    &
                   (xk(ipol, ik) + g(ipol, igg) )
           enddo
           do ibnd = 1, nbnd
              dspsi(1:npw,ibnd) = ps2(ikb,ipol,ibnd) * aux(1:npw) &
                   + dspsi(1:npw,ibnd)
           enddo
        endif
     enddo

  enddo
  do ibnd = 1, nbnd
     do jbnd=1, nbnd
        pdsp(ibnd,jbnd) = dot_product(evc(1:npw,ibnd),dspsi(1:npw,jbnd))
     enddo
  enddo

  if (allocated(aux)) deallocate (aux)
  if (allocated(ps2)) deallocate (ps2)
  if (allocated(ps1)) deallocate (ps1)
  if (allocated(dspsi)) deallocate (dspsi)

  return
end subroutine psidspsi
