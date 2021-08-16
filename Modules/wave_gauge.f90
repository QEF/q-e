!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

MODULE wave_gauge
   USE kinds, ONLY: DP

   contains

   subroutine compute_dot_evc_parallel_gauge(t_minus, t_zero, t_plus, dot_evc, nbnd, npw, npwx, gstart)
      !! compute the derivative of the wavefunction multiplied by dt in the parallel transport gauge
      !! t_minus, t_zero, t_plus are the wavefunctions. t_minus and t_plus are separated by a time dt
      !! t_zero is computed in the middle. The gauge of the result is the same of t_zero.
      !! This works also if you put t_zero = t_minus or t_zero = t_plus (be careful to the gauge!)
      !! on output you have dot_evc*dt (the derivative multiplied by dt) with
      !! the gauge of t_zero. Gauges of t_minus and t_plus should not enter the
      !! result and can be arbitrary
      use mp, only: mp_sum
      USE mp_pools, ONLY: intra_pool_comm

      integer, intent(in) :: nbnd, npw, npwx, gstart
      COMPLEX(DP), intent(inout) ::  dot_evc(:, :)
      complex(dp), intent(in) :: t_minus(:,:), t_zero(:,:), t_plus(:,:)

      real(DP) :: sa(nbnd, nbnd), ssa(nbnd, nbnd), sb(nbnd, nbnd), ssb(nbnd, nbnd)

      sa = 0.d0
      sb = 0.d0

! computed sb_ij = < t_plus_i, t_zero_j >, sa_ij = < t_minus_i, t_zero_j > remove contribution at G=0
! sa is a multiple of the identity if 2 points are used (t_plus==t_zero)
      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, t_plus, 2*npwx, t_zero, 2*npwx, 0.d0, sa, nbnd)
      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, t_minus, 2*npwx, t_zero, 2*npwx, 0.d0, sb, nbnd)
      if (gstart == 2) then
         do ibnd = 1, nbnd
            do jbnd = 1, nbnd
               sa(ibnd, jbnd) = sa(ibnd, jbnd) - dble(conjg(t_plus(1, ibnd))*t_zero(1, jbnd))
               sb(ibnd, jbnd) = sb(ibnd, jbnd) - dble(conjg(t_minus(1, ibnd))*t_zero(1, jbnd))
            end do
         end do
      end if
      call mp_sum(sa, intra_pool_comm)
      call mp_sum(sb, intra_pool_comm)
! computes ssa = sa ^ 2 and ssb = sb ^ 2. The gauge cancels out here

      call dgemm('T', 'N', nbnd, nbnd, nbnd, 1.d0, sa, nbnd, sa, nbnd, 0.d0, ssa, nbnd)
      call dgemm('T', 'N', nbnd, nbnd, nbnd, 1.d0, sb, nbnd, sb, nbnd, 0.d0, ssb, nbnd)

      dot_evc = 0.d0
      do ibnd = 1, nbnd
         do jbnd = 1, nbnd
            do ig = 1, npw

            dot_evc(ig, ibnd) = dot_evc(ig, ibnd) + t_plus(ig, jbnd)*sa(jbnd, ibnd)  &
                                                  - t_zero(ig, jbnd)*ssa(ibnd, jbnd) &
                                                  - t_minus(ig, jbnd)*sb(jbnd, ibnd) &
                                                  + t_zero(ig, jbnd)*ssb(ibnd, jbnd)

            end do
         end do
      end do



   end subroutine


end module
