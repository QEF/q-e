!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine drhodv (nu_i0, nper, drhoscf)
  !-----------------------------------------------------------------------
  !
  !    This subroutine computes the electronic term
  !    <psi|dv|dpsi> of the dynamical matrix
  !
#include "machine.h"
  !
  USE ions_base, ONLY : nat
  use pwcom
  USE kinds, only : DP
  USE io_files, ONLY: iunigk
  use phcom
  implicit none

  integer :: nper, nu_i0
  ! input: number of perturbations of this represent
  ! input: the initial position of the mode

  complex(kind=DP) :: drhoscf (nrxx, nspin, npertx)
  ! the change of density due to perturbations

  integer :: mu, ik, ikq, ig, nu_i, nu_j, na_jcart, ibnd, nrec, &
       ipol, ikk
  ! counters
  ! ikk: record position for wfc at k

  complex(kind=DP) :: fact, ps, dynwrk (3 * nat, 3 * nat), &
       wdyn (3 * nat, 3 * nat), ZDOTC
  complex(kind=DP), allocatable ::  aux (:,:), dbecq (:,:,:), &
       dalpq (:,:,:,:)
  ! work space
  !
  !   Initialize the auxiliary matrix wdyn
  !
  call start_clock ('drhodv')
  allocate (dbecq ( nkb , nbnd, nper))    
  allocate (dalpq ( nkb , nbnd ,3 ,nper))    
  allocate (aux   ( npwx , nbnd))    
  dynwrk(:,:) = (0.d0, 0.d0)
  wdyn  (:,:) = (0.d0, 0.d0)
  !
  !   We need a sum over all k points ...
  !
  if (nksq > 1) rewind (unit = iunigk)
  do ik = 1, nksq
     if (nksq > 1) read (iunigk) npw, igk
     if (lgamma) then
        ikk = ik
        ikq = ik
        npwq = npw
     else
        ikk = 2 * ik - 1
        ikq = ikk + 1
        if (nksq > 1) read (iunigk) npwq, igkq
     endif
     if (lsda) current_spin = isk (ikk)
     call init_us_2 (npwq, igkq, xk (1, ikq), vkb)

     do mu = 1, nper
        nrec = (mu - 1) * nksq + ik
        if (nksq > 1 .or. nper > 1) call davcio(dpsi, lrdwf, iudwf, nrec,-1)
        call ccalbec (nkb, npwx, npwq, nbnd, dbecq (1, 1, mu), vkb, dpsi)
        do ipol = 1, 3
           do ibnd = 1, nbnd
              do ig = 1, npwq
                 aux (ig, ibnd) = dpsi (ig, ibnd) * &
                      (xk (ipol, ikq) + g (ipol, igkq (ig) ) )
              enddo
           enddo
           call ccalbec (nkb, npwx, npwq, nbnd, dalpq(1,1,ipol,mu), vkb, aux)
        enddo
     enddo
     fact = DCMPLX (0.d0, tpiba)
     dalpq = dalpq * fact
     call drhodvnl (ik, ikk, nper, nu_i0, dynwrk, dbecq, dalpq)
  enddo
  !
  !   put in the basis of the modes
  !
  do nu_i = 1, 3 * nat
     do nu_j = 1, 3 * nat
        ps = (0.0d0, 0.0d0)
        do na_jcart = 1, 3 * nat
           ps = ps + dynwrk (nu_i, na_jcart) * u (na_jcart, nu_j)
        enddo
        wdyn (nu_i, nu_j) = wdyn (nu_i, nu_j) + ps
     enddo

  enddo
#ifdef __PARA
  !
  ! collect contributions from all pools (sum over k-points)
  !
  call poolreduce (18 * nat * nat, wdyn)
#endif
  !
  ! add the contribution of the local part of the perturbation
  !
  call drhodvloc (nu_i0, nper, drhoscf, wdyn)
  !
  ! add to the rest of the dynamical matrix
  !
  !      WRITE( stdout,*) 'drhodv dyn, wdyn'
  !      call tra_write_matrix('drhodv dyn',dyn,u,nat)
  !      call tra_write_matrix('drhodv wdyn',wdyn,u,nat)

  dyn (:,:) = dyn (:,:) + wdyn (:,:) 

  deallocate (aux)
  deallocate (dalpq)
  deallocate (dbecq)

  call stop_clock ('drhodv')
  return
end subroutine drhodv
