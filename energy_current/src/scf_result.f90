module scf_result_mod
   USE kinds, ONLY: DP

type scf_result
      complex(kind=dp), allocatable :: evc(:,:)
      real(kind=dp), allocatable :: tau(:,:), vel(:,:)
   end type
   type multiple_scf_result
      type(scf_result) :: t_minus, t_zero, t_plus
   end type
 contains

   subroutine scf_result_allocate(t, npwx, nbnd, natoms)
       implicit none
       type(scf_result), intent(inout) :: t
       integer, intent(in) :: npwx, nbnd, natoms
       allocate(t%evc(npwx, nbnd))
       allocate(t%tau(3,natoms))
       allocate(t%vel(3,natoms))
   end subroutine

   subroutine scf_result_deallocate(t)
       implicit none
       type(scf_result), intent(inout) :: t
       deallocate(t%evc)
       deallocate(t%vel)
       deallocate(t%tau)
   end subroutine

   subroutine scf_result_set_from_global_variables(t)
       use wavefunctions, only: evc
       use ions_base, only : tau
       use dynamics_module, only: vel
       implicit none
       type(scf_result), intent(inout) :: t

       t%evc=evc
       t%tau=tau
       t%vel=vel

   end subroutine

   subroutine scf_result_set_tau_vel_global_variables(t)

       use ions_base, only : tau
       use dynamics_module, only: vel
       implicit none
       type(scf_result), intent(in) :: t

       tau=t%tau
       vel=t%vel

   end subroutine
   subroutine scf_result_set_global_variables(t)

       use wavefunctions, only: evc
       use ions_base, only : tau
       use dynamics_module, only: vel
       implicit none
       type(scf_result), intent(in) :: t

       evc=t%evc
       tau=t%tau
       vel=t%vel

   end subroutine

   subroutine multiple_scf_result_allocate(t,allocate_zero)
      use wvfct, only: nbnd, npwx
      use ions_base, only: nat
      implicit none
      type(multiple_scf_result), intent(inout) :: t
      logical, intent(in) :: allocate_zero

      call multiple_scf_result_allocate_(t, npwx, nbnd, nat, allocate_zero)

   end subroutine

   subroutine multiple_scf_result_allocate_(t, npwx, nbnd, natoms, allocate_zero)
      implicit none
      type(multiple_scf_result), intent(inout) :: t
      integer, intent(in) :: npwx, nbnd, natoms
      logical, intent(in) :: allocate_zero

      if (allocate_zero) call scf_result_allocate(t%t_zero, npwx, nbnd, natoms)
      call scf_result_allocate(t%t_minus, npwx, nbnd, natoms)
      call scf_result_allocate(t%t_plus, npwx, nbnd, natoms)

   end subroutine

   subroutine multiple_scf_result_deallocate(t)
      implicit none
      type(multiple_scf_result), intent(inout) :: t
      if (allocated(t%t_zero%evc )) call scf_result_deallocate(t%t_zero)
      call scf_result_deallocate(t%t_minus)
      call scf_result_deallocate(t%t_plus)

   end subroutine

end module

