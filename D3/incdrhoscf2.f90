!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine incdrhoscf2 (drhoscf, weight, ik, dbecsum, mode, flag)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the charge density due to the
  !     perturbation. It is called at the end of the computation of the
  !     change of the wavefunction for a given k point.
  !
  !
#include "machine.h"
  use pwcom
  use phcom

  implicit none
  integer :: ik
  ! input: the k point

  real (8) :: weight
  ! input: the weight of the k point
  complex (8) :: drhoscf (nrxxs), dbecsum (nhm * (nhm + 1) / 2, nat)
  ! output: the change of the charge densit
  ! inp/out: the accumulated dbec
  integer :: mode, flag
  ! flag =1 if dpsi is used (in solve_linte
  ! flag!=1 if dpsi is not used (in addusdd
  !
  !   here the local variable
  !

  real (8) :: wgt
  ! the effective weight of the k point

  complex (8), allocatable :: psi (:), dpsic (:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real sp

  integer :: ibnd, jbnd, ikk, ir, ig
  ! counter on bands
  ! counter on bands
  ! the record ik
  ! counter on mesh points
  ! counter on G vectors

  call start_clock ('incdrhoscf')
  allocate  (dpsic( nrxxs))    
  allocate  (psi  ( nrxxs))    
  wgt = 2.d0 * weight / omega
  if (lgamma) then
     ikk = ik
  else
     ikk = 2 * ik - 1
  endif
  !
  ! dpsi contains the   perturbed wavefunctions of this k point
  ! evc  contains the unperturbed wavefunctions of this k point
  !
  !      do ibnd = 1,nbnd_occ(ikk)
  do ibnd = 1, nbnd
     call setv (2 * nrxxs, 0.d0, psi, 1)
     do ig = 1, npw
        psi (nls (igk (ig) ) ) = evc (ig, ibnd)
     enddo
     call cft3s (psi, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
     call setv (2 * nrxxs, 0.d0, dpsic, 1)
     !
     !    here we add the term in the valence due to the change of the
     !    constraint. dvpsi is used as work space, dpsi is unchanged
     !
     if (flag.eq.1) then
        call ZCOPY (npwx, dpsi (1, ibnd), 1, dvpsi (1, ibnd), 1)
     else
        call setv (2 * npwx, 0.d0, dvpsi (1, ibnd), 1)
     endif
     !         call ZGEMM('N','N', npwq, nbnd, nbnd, (1.d0,0.d0),
     !     +              evq, npwx, prodval(1,1,mode),nbnd,
     !     +             (1.d0,0.d0),dvpsi,npwx)
     if (okvan) then
        call error ('incdrhoscf2', 'US not allowed', 1)
        !            do jbnd=1,nbnd
        !               call ZAXPY(npwq,prodval(jbnd,ibnd,mode),
        !     +           evq(1,jbnd),1,dvpsi(1,ibnd),1)
        !            enddo
     endif
     do ig = 1, npwq
        dpsic (nls (igkq (ig) ) ) = dvpsi (ig, ibnd)
     enddo

     call cft3s (dpsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
     do ir = 1, nrxxs
        drhoscf (ir) = drhoscf (ir) + wgt * conjg (psi (ir) ) * dpsic (ir)
        !            if (ir.lt.20) write (6,*)   drhoscf(ir)
     enddo

  enddo
  call addusdbec (ik, wgt, dvpsi, dbecsum)
  !      write(6,*) '*********************'
  !      do ig=1,20
  !         write(6,*) dbecsum(ig,1)
  !      enddo
  !      call stoallocate  (ph(.true.))    
  deallocate (psi)
  deallocate (dpsic)

  call stop_clock ('incdrhoscf')
  return
end subroutine incdrhoscf2
