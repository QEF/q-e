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
! this module contains routines that are used to compute the numerical
! derivative of the wavefunctions fixing a particular gauge

   USE kinds, ONLY: DP
   implicit none


   contains

   subroutine compute_dot_evc_parallel_gauge(t_minus, t_zero, t_plus, dot_evc, nbnd, npw, npwx, gstart)
      ! compute the numerical derivative using the wavefunctions computed at t-dt/2, t, t+dt/2 
      ! in the parallel transport gauge

      integer, intent(in) :: nbnd, npw, npwx, gstart
      COMPLEX(DP), intent(inout) ::  dot_evc(:, :)
      complex(dp), intent(in) :: t_minus(:,:), t_zero(:,:), t_plus(:,:)
      !$acc declare present(t_minus, t_zero, t_plus, dot_evc)

      !$acc parallel
      dot_evc=(0.0_dp, 0.0_dp)
      !$acc end parallel
      call project_parallel_gauge(t_minus, t_zero, t_plus, dot_evc, nbnd, npw, npwx, gstart, 1.0_dp,&
                                  .true., .true.)
      
   end subroutine

   subroutine project_parallel_gauge_2(t_minus, t_zero, t_zero_proj , nbnd, npw, npwx, gstart)
      ! project t_zero over the manifold of t_minus in the parallel transport gauge.
      ! Do not project over the not occupied manifold

      integer, intent(in) :: nbnd, npw, npwx, gstart
      COMPLEX(DP), intent(inout) ::   t_zero_proj(:,:)
      complex(dp), intent(in) :: t_minus(:, :),  t_zero(:,:)
      complex(dp), allocatable :: dummy(:,:)
      !$acc declare device_resident(dummy)
      !$acc declare present(t_minus, t_zero, t_zero_proj)

      !$acc parallel
      t_zero_proj=0.0_dp
      !$acc end parallel
      call project_parallel_gauge(t_zero, t_minus, dummy, t_zero_proj,  nbnd, npw, npwx,&
                                  gstart, -1.0_dp, .false., .false. )

   end subroutine 

   subroutine project_parallel_gauge(t_minus, t_zero, t_plus, dot_evc, nbnd, npw, npwx, gstart, &
                                     factor, use_t_plus, project_conduction)
      !! implemented only in the plain dft norm conserving pseudopotential case 
      !! let P be the projector over the occupied manifold space:
      !! P = \sum_v |c_v><c_v|
      !! then (no surprise)
      !! |c> = P |c>
      !! doing a time derivative:
      !! \dot |c> = \dot P |c> + P \dot |c>
      !! the parallel transport gauge fixes P \dot |c> to zero (the derivative
      !! has no components over the occupied space). \dot P can be computed numerically
      !! using wavefunctions that possibly have different gauges.
      !! Then, to make sure that there are no components over the occupied space,
      !! a projector over the conduction band is added if project_conduction is .true.:
      !!
      !!  (1-P) \dot P |c> 
      !! 
      !! t_minus, t_zero, t_plus are the wavefunctions. t_minus and t_plus are separated by a time dt
      !! t_zero is computed in the middle. The gauge of the result is the same of t_zero.
      !! This works also if you put t_zero = t_minus or t_zero = t_plus (be careful to the gauge!).
      !! This routine adds \dot |c>*dt*factor (the derivative multiplied by dt multiplied by the specified factor)
      !! to the output array dot_evc with the gauge of t_zero. Gauges of t_minus and t_plus should not enter the
      !! result and can be arbitrary. If use_t_plus is false, t_plus = t_zero is assumed and some cpu time
      !! is saved. If project_conduction is true the (1-P) projector is computed, and this involves an additional
      !! matrix-matrix product with dimensions nbnd x nbnd
      use mp, only: mp_sum
      USE mp_pools, ONLY: intra_pool_comm

      integer, intent(in) :: nbnd, npw, npwx, gstart
      COMPLEX(DP), intent(inout) ::  dot_evc(:, :)
      complex(dp), intent(in) :: t_minus(:,:), t_zero(:,:), t_plus(:,:)
      real(dp), intent(in) :: factor
      logical, intent(in) :: use_t_plus, project_conduction
      !$acc declare present(t_minus, t_zero, t_plus, dot_evc)

      real(DP), allocatable :: sa(:,:), ssa(:,:), sb(:,:), ssb(:,:)
      !$acc declare device_resident(ssa, ssb)
      complex(dp) :: tmp,tmp2
      integer :: ibnd, jbnd, ig


      allocate(sb(nbnd, nbnd))
! \dot P |c> = 
! computed sa_ij = < t_plus_i, t_zero_j >, sb_ij = < t_minus_i, t_zero_j >
! remove contribution at G=0 due to c(-g)*=c(g), and only half vector is stored
! (gamma trick)
! sa is the identity if 2 points are used (t_plus==t_zero)
      if (use_t_plus) &
         allocate(sa(nbnd, nbnd))
      !$acc data create(sb, sa)
      if (use_t_plus) then
#if defined (__CUDA)  
         !$acc host_data use_device(t_plus,t_zero, sa)
         call mydgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, t_plus, 2*npwx, t_zero, 2*npwx, 0.d0, sa, nbnd)
         !$acc end host_data
#else
         call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, t_plus, 2*npwx, t_zero, 2*npwx, 0.d0, sa, nbnd) 
#endif 
      end if
#if defined (__CUDA) 
      !$acc host_data use_device(t_minus,t_zero, sb)
      call mydgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, t_minus, 2*npwx, t_zero, 2*npwx, 0.d0, sb, nbnd)
      !$acc end host_data
#else 
      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, t_minus, 2*npwx, t_zero, 2*npwx, 0.d0, sb, nbnd)
#endif 
      if (gstart == 2) then
         !$acc parallel loop collapse(2) present(sa,sb,t_plus,t_minus,t_zero)
         do ibnd = 1, nbnd
            do jbnd = 1, nbnd
               if (use_t_plus) &
                  sa(ibnd, jbnd) = sa(ibnd, jbnd) - dble(conjg(t_plus(1, ibnd))*t_zero(1, jbnd))
               sb(ibnd, jbnd) = sb(ibnd, jbnd) - dble(conjg(t_minus(1, ibnd))*t_zero(1, jbnd))
            end do
         end do
      end if
      !$acc host_data use_device(sb)
      call mp_sum(sb, intra_pool_comm)
      !$acc end host_data
      if (use_t_plus) then
         !$acc host_data use_device(sa)
         call mp_sum(sa, intra_pool_comm)
         !$acc end host_data
      end if

      ! compute scalar products that appear due to the ( 1 - P ) projector over
      ! the empty bands, that are 
      ! ssa = sa ^ 2 and ssb = sb ^ 2. The gauge cancels out here
      if (project_conduction) then
         allocate(ssb(nbnd, nbnd))
         if (use_t_plus) then
             allocate(ssa(nbnd, nbnd))
#if defined (__CUDA) 
             !$acc host_data use_device(sa,ssa)
             call mydgemm('T', 'N', nbnd, nbnd, nbnd, 1.d0, sa, nbnd, sa, nbnd, 0.d0, ssa, nbnd)
             !$acc end host_data
#else 
             call dgemm('T', 'N', nbnd, nbnd, nbnd, 1.d0, sa, nbnd, sa, nbnd, 0.d0, ssa, nbnd)
#endif 
         end if
#if defined (__CUDA) 
         !$acc host_data use_device(sb,ssb)
         call mydgemm('T', 'N', nbnd, nbnd, nbnd, 1.d0, sb, nbnd, sb, nbnd, 0.d0, ssb, nbnd)
         !$acc end host_data
#else 
         call dgemm('T', 'N', nbnd, nbnd, nbnd, 1.d0, sb, nbnd, sb, nbnd, 0.d0, ssb, nbnd)
#endif
      end if
      ! compute final projection
      !$acc parallel loop present(t_minus,sb,t_plus,sa,t_zero,ssa,ssb,dot_evc) copyin(factor) private(tmp2,tmp)
      do ibnd = 1, nbnd
         do jbnd = 1, nbnd
            do ig = 1, npw
                tmp=0.0_dp
                tmp = tmp - t_minus(ig, jbnd)*sb(jbnd, ibnd)
                if (use_t_plus) then
                   tmp = tmp + t_plus(ig, jbnd)*sa(jbnd, ibnd)
                endif
                if (project_conduction) then
                   if (use_t_plus)&
                       tmp = tmp - t_zero(ig, jbnd)*ssa(ibnd, jbnd)
                   tmp = tmp  + t_zero(ig, jbnd)*ssb(ibnd, jbnd)
                endif

                dot_evc(ig, ibnd) = dot_evc(ig, ibnd) + factor*tmp
            end do
         end do
      end do
      !$acc end data
      if (project_conduction) then
         deallocate(ssb)
         if (use_t_plus)&
             deallocate(ssa)
      end if
      if (use_t_plus)&
          deallocate(sa)
      deallocate(sb)

   end subroutine


end module
