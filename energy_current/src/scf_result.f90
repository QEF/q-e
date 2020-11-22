module scf_result_mod
   USE kinds, ONLY: DP
   use becmod, only : bec_type, allocate_bec_type, deallocate_bec_type,&
                      beccopy

   type scf_result
      complex(kind=dp), allocatable :: evc(:, :)
      REAL(DP), ALLOCATABLE :: vrs(:,:)
      type(bec_type) :: becp
      complex(dp),allocatable :: vkb(:,:) 
      real(kind=dp), allocatable :: tau(:, :), vel(:, :)
   end type
   type multiple_scf_result
      type(scf_result) :: t_minus, t_zero, t_plus
   end type
contains

   subroutine scf_result_allocate(t, npwx, nbnd, natoms, nnr, nspin, nkb)
      implicit none
      type(scf_result), intent(inout) :: t
      integer, intent(in) :: npwx, nbnd, natoms, nnr, nspin, nkb
      allocate (t%evc(npwx, nbnd))
      allocate (t%vrs(nnr,nspin))
      allocate (t%tau(3, natoms))
      allocate (t%vel(3, natoms))
      allocate (t%vkb(npwx, nkb))
      call allocate_bec_type(nkb, nbnd, t%becp)
   end subroutine

   subroutine scf_result_deallocate(t)
      implicit none
      type(scf_result), intent(inout) :: t
      deallocate (t%evc)
      deallocate (t%vrs)
      call deallocate_bec_type(t%becp)
      deallocate (t%vel)
      deallocate (t%tau)
      deallocate (t%vkb)
   end subroutine

   subroutine scf_result_set_from_global_variables(t)
      use wavefunctions, only: evc
      use ions_base, only: tau
      use dynamics_module, only: vel
      use scf, only : vrs
      use uspp, only: vkb
      use becmod, only : becp
      implicit none
      type(scf_result), intent(inout) :: t

      t%evc = evc
      t%vrs = vrs
      t%vkb = vkb
      if (allocated(becp%r)) &
         t%becp%r = becp%r
      if (allocated(becp%k)) &
         t%becp%k = becp%k
      if (allocated(becp%nc)) &
         t%becp%nc = becp%nc
      t%becp%comm = becp%comm
      t%becp%nbnd = becp%nbnd
      t%becp%nproc = becp%nproc
      t%becp%mype = becp%mype
      t%becp%nbnd_loc = becp%nbnd_loc
      t%becp%ibnd_begin = becp%ibnd_begin
      t%tau = tau
      t%vel = vel

   end subroutine

   subroutine scf_result_set_tau_vel_global_variables(t)

      use ions_base, only: tau
      use dynamics_module, only: vel
      implicit none
      type(scf_result), intent(in) :: t

      tau = t%tau
      vel = t%vel

   end subroutine
   subroutine scf_result_set_global_variables(t)

      use wavefunctions, only: evc
      use ions_base, only: tau
      use dynamics_module, only: vel
      use scf, only : vrs
      use uspp, only : vkb
      use becmod, only : becp
      implicit none
      type(scf_result), intent(in) :: t

      evc = t%evc
      vrs = t%vrs
      vkb = t%vkb
      if (allocated(becp%r)) &
         becp%r = t%becp%r
      if (allocated(becp%k)) &
         becp%k = t%becp%k
      if (allocated(becp%nc)) &
         becp%nc = t%becp%nc
      becp%comm =       t%becp%comm
      becp%nbnd =       t%becp%nbnd
      becp%nproc=       t%becp%nproc
      becp%mype =       t%becp%mype
      becp%nbnd_loc =   t%becp%nbnd_loc
      becp%ibnd_begin = t%becp%ibnd_begin
      tau = t%tau
      vel = t%vel

   end subroutine

   subroutine multiple_scf_result_allocate(t, allocate_zero)
      use wvfct, only: nbnd, npwx
      use ions_base, only: nat
      use fft_base, only: dffts
      use uspp, ONLY: nkb
      USE lsda_mod, ONLY : nspin
      implicit none
      type(multiple_scf_result), intent(inout) :: t
      logical, intent(in) :: allocate_zero

      call multiple_scf_result_allocate_(t, npwx, nbnd, nat, dffts%nnr, nspin, nkb, allocate_zero)

   end subroutine

   subroutine multiple_scf_result_allocate_(t, npwx, nbnd, natoms, nnr, nspin, nkb, allocate_zero)
      implicit none
      type(multiple_scf_result), intent(inout) :: t
      integer, intent(in) :: npwx, nbnd, natoms, nnr, nspin, nkb
      logical, intent(in) :: allocate_zero

      if (allocate_zero) call scf_result_allocate(t%t_zero, npwx, nbnd, natoms, nnr, nspin, nkb)
      call scf_result_allocate(t%t_minus, npwx, nbnd, natoms, nnr, nspin, nkb)
      call scf_result_allocate(t%t_plus, npwx, nbnd, natoms, nnr, nspin, nkb)

   end subroutine

   subroutine multiple_scf_result_deallocate(t)
      implicit none
      type(multiple_scf_result), intent(inout) :: t
      if (allocated(t%t_zero%evc)) call scf_result_deallocate(t%t_zero)
      call scf_result_deallocate(t%t_minus)
      call scf_result_deallocate(t%t_plus)

   end subroutine

end module

