!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine ch_psi_all2 (n, h, ah, e, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !
  USE kinds, only : DP
  use pwcom
  USE uspp,      ONLY: vkb
  use becmod
  use phcom
  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum

  use control_lr, ONLY : alpha_pv, lgamma

  implicit none
  integer :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point

  real (DP) :: e (m)
  ! input: the eigenvalue

  complex (DP) :: h (npwx, m), ah (npwx, m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  !   local variables
  !

  integer :: ibnd, ikq, ig
  ! counter on bands
  ! the point k+q
  ! counter on G vetors

  complex (DP), allocatable :: ps (:,:), hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h


  call start_clock ('ch_psi')
  allocate (ps( nbnd, m))
  allocate (hpsi( npwx, m))
  allocate (spsi( npwx, m))

  hpsi = (0.d0, 0.d0)
  spsi = (0.d0, 0.d0)
  !
  !   compute the product of the hamiltonian with the h vector
  !
  call h_psiq (npwx, n, m, h, hpsi, spsi)

  call start_clock ('last')
  !
  !   then we compute the operator H-epsilon S
  !
  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd) = hpsi (ig, ibnd) - e (ibnd) * spsi (ig, ibnd)
     enddo
  enddo
  !
  !   Here we compute the projector in the valence band
  !
  hpsi = (0.d0, 0.d0)
  if (lgamma) then
     ikq = ik
  else
     ikq = 2 * ik
  endif
  ps = (0.d0, 0.d0)

  call zgemm ('C', 'N', nbnd, m, n, (1.d0, 0.d0) , evq, npwx, spsi, &
       npwx, (0.d0, 0.d0) , ps, nbnd)
  ps = ps * alpha_pv
#ifdef __MPI
  call mp_sum(  ps, intra_pool_comm )
#endif

  call zgemm ('N', 'N', n, m, nbnd, (1.d0, 0.d0) , evq, npwx, ps, &
       nbnd, (1.d0, 0.d0) , hpsi, npwx)
  spsi = hpsi
  !
  !    And apply S again
  !
  call calbec (n, vkb, hpsi, becp, m)
  call s_psi (npwx, n, m, hpsi, spsi)
  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
     enddo
  enddo

  deallocate (spsi)
  deallocate (hpsi)
  deallocate (ps)
  call stop_clock ('last')
  call stop_clock ('ch_psi')
  return
end subroutine ch_psi_all2
