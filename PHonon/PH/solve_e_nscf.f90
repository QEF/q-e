!
! Copyright (C) 2001-208 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine solve_e_nscf( avg_iter, thresh, ik, ipol, dvscfs, auxr )
  !-----------------------------------------------------------------------
  !
  !   Solve the linear system which defines the change of the wavefunctions
  !   due to the electric field for a given k_point in a non self-consistent
  !   way. The self-consistent variation of the potential has been computed
  !   previously and is in dvscfs.
  !
  use kinds,                 ONLY : DP
  USE cell_base,             ONLY : tpiba2
  USE klist,                 ONLY : xk, ngk, igk_k
  USE fft_base,              ONLY : dffts
  USE fft_interfaces,        ONLY : fwfft, invfft
  USE buffers,               ONLY : get_buffer
  USE gvect,                 ONLY : g
  USE gvecs,                 ONLY : nls
  USE wvfct,                 ONLY : et
  USE wavefunctions_module,  ONLY : evc
  USE eqv,                   ONLY : dpsi, dvpsi
  USE units_ph,              ONLY : this_pcxpsi_is_on_file, lrdwf, iudwf
  USE qpoint,                ONLY : nksq
  USE control_lr,            ONLY : nbnd_occ
  implicit none

  !
  !  Input variables
  !
  integer :: ik, ipol
  ! input: k-point under consideration
  ! input: polarization of the electric field

  real(DP) :: thresh, avg_iter
  ! input: convergence threshold
  ! in/out: # of diagonalization iterations

  complex(DP) :: dvscfs (dffts%nnr, 3), auxr(dffts%nnr)
  ! input: potential on the smooth grid
  ! auxiliary work space

  !
  !  Local variables
  !
  integer :: npw, npwq, ibnd, ir, ig, nrec
  ! counter on bands
  ! counter on mesh points
  ! counter on G-points
  ! the record number

  !
  ! Calculates [H,x]*psi_kpoint
  !
  dpsi (:,:) = (0.d0, 0.d0)
  this_pcxpsi_is_on_file(:,:)=.false.
  call dvpsi_e (ik, ipol)
  npw = ngk(ik)
  npwq= npw     ! note: q=0
  !
  CALL g2_kin (ik) ! needed by pcgreen
  !
  ! Calculates dvscf*psi_k in G_space
  !
  do ibnd = 1, nbnd_occ (ik)
     auxr (:) = (0.d0, 0.d0)
     do ig = 1, npw
        auxr (nls (igk_k (ig,ik))) = evc (ig, ibnd)
     end do
     CALL invfft ('Wave', auxr, dffts)
     do ir = 1, dffts%nnr
        auxr (ir) = auxr(ir) * dvscfs(ir, ipol)
     end do
     CALL fwfft ('Wave', auxr, dffts)
     do ig = 1, npwq
        ! note: q=0
        dvpsi (ig, ibnd) = dvpsi(ig, ibnd) + auxr(nls (igk_k (ig,ik)))
     enddo
  enddo
  !
  ! starting value for  delta_psi is read from iudwf
  !
  nrec = (ipol - 1) * nksq + ik
  call get_buffer (dpsi, lrdwf, iudwf, nrec)
  call pcgreen (avg_iter, thresh, ik, et (1, ik))
!
!  The pcxpsi on file could be at k+dk and cannot be used by the following
!  codes that require pcxpsi at k.
!
  this_pcxpsi_is_on_file(ik,ipol)=.false.

  return
end subroutine solve_e_nscf
