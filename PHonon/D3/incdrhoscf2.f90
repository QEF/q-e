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
  USE ions_base,  ONLY : nat
  USE kinds, only : DP
  USE fft_base,   ONLY : dffts
  USE fft_interfaces, ONLY : invfft
  use pwcom
  USE wavefunctions_module,  ONLY: evc
  USE uspp, ONLY: okvan
  USE uspp_param, ONLY: nhm
  USE qpoint,     ONLY: npwq, igkq
  use phcom

  implicit none
  integer :: ik
  ! input: the k point

  real (DP) :: weight
  ! input: the weight of the k point
  complex (DP) :: drhoscf (dffts%nnr), dbecsum (nhm * (nhm + 1) / 2, nat)
  ! output: the change of the charge densit
  ! inp/out: the accumulated dbec
  integer :: mode, flag
  ! flag =1 if dpsi is used (in solve_linte
  ! flag!=1 if dpsi is not used (in addusdd
  !
  !   here the local variable
  !

  real (DP) :: wgt
  ! the effective weight of the k point

  complex (DP), allocatable :: psi (:), dpsic (:)
  ! the wavefunctions in real space
  ! the change of wavefunctions in real space

  integer :: ibnd, jbnd, ikk, ir, ig
  ! counters

  call start_clock ('incdrhoscf')
  allocate  (dpsic( dffts%nnr))
  allocate  (psi  ( dffts%nnr))
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
     psi (:) = (0.d0, 0.d0)
     do ig = 1, npw
        psi (nls (igk (ig) ) ) = evc (ig, ibnd)
     enddo
     CALL invfft ('Wave', psi, dffts)
     dpsic(:) =(0.d0, 0.d0)
     !
     !    here we add the term in the valence due to the change of the
     !    constraint. dvpsi is used as work space, dpsi is unchanged
     !
     if (flag == 1) then
        dvpsi (:, ibnd) = dpsi (:, ibnd)
     else
        dvpsi (:, ibnd) = (0.d0, 0.d0)
     endif
     !         call zgemm('N','N', npwq, nbnd, nbnd, (1.d0,0.d0),
     !     +              evq, npwx, prodval(1,1,mode),nbnd,
     !     +             (1.d0,0.d0),dvpsi,npwx)
     if (okvan) then
        call errore ('incdrhoscf2', 'US not allowed', 1)
        !            do jbnd=1,nbnd
        !               call zaxpy(npwq,prodval(jbnd,ibnd,mode),
        !     +           evq(1,jbnd),1,dvpsi(1,ibnd),1)
        !            enddo
     endif
     do ig = 1, npwq
        dpsic (nls (igkq (ig) ) ) = dvpsi (ig, ibnd)
     enddo

     CALL invfft ('Wave', dpsic, dffts)
     do ir = 1, dffts%nnr
        drhoscf (ir) = drhoscf (ir) + wgt * CONJG(psi (ir) ) * dpsic (ir)
        !            if (ir.lt.20) WRITE( stdout,*)   drhoscf(ir)
     enddo

  enddo
  call addusdbec (ik, wgt, dvpsi, dbecsum)
  !      WRITE( stdout,*) '*********************'
  !      do ig=1,20
  !         WRITE( stdout,*) dbecsum(ig,1)
  !      enddo
  !      call stoallocate  (ph(.true.))
  deallocate (psi)
  deallocate (dpsic)

  call stop_clock ('incdrhoscf')
  return
end subroutine incdrhoscf2
