MODULE hartree_mod

   USE kinds, ONLY: DP

   SAVE

   type scf_result
      complex(kind=dp), allocatable :: evc(:,:)
      real(kind=dp), allocatable :: tau(:,:), vel(:,:)
   end type
   type multiple_scf_result
      type(scf_result) :: t_minus, t_zero, t_plus
   end type

   CHARACTER(len=256) :: file_output, trajdir = '' !, init_linear
   real(kind=DP) ::J_kohn(3), J_kohn_a(3), J_kohn_b(3), J_hartree(3), J_xc(3), J_electron(3)
   real(kind=DP), allocatable :: v_cm(:, :)

   real(kind=DP) ::delta_t, ethr_small_step, ethr_big_step
   !complex(kind=DP), allocatable :: evc_uno(:, :), evc_due(:, :) ! TODO: maybe save one allocation of npwx*nb
   type(multiple_scf_result) :: scf_all

   complex(kind=DP), allocatable :: dvpsi_save(:,:,:) ! to save the solution of the system between iterations
   logical :: save_dvpsi = .true. ! if true dvpsi_save is allocated and used

   integer :: first_step, last_step, step_mul, step_rem, n_repeat_every_step
   logical :: restart ! if true try to read last calculated step from output and set first_step
   logical :: subtract_cm_vel ! if true do velocity renormalization
   logical :: ec_test ! activates tests for debugging purposes
   logical :: re_init_wfc_1=.false., re_init_wfc_2=.false. ! initialize again evc before scf step number 1 or 2



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
       type(scf_result), intent(inout) :: t

       tau=t%tau
       vel=t%vel

   end subroutine
   subroutine scf_result_set_global_variables(t)

       use wavefunctions, only: evc
       use ions_base, only : tau
       use dynamics_module, only: vel
       implicit none
       type(scf_result), intent(inout) :: t

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

   subroutine compute_hartree_vpoint_vmean(v_mean,v_point, g, gstart, omega, ngm, e2, fpi, tpiba2, &
                   delta_t, charge_g,charge_g_due, charge_g_tre)
      use kinds, only: DP
      implicit none
      integer, intent(in) :: gstart,ngm
      real(kind=dp), intent(in) :: g(:,:), omega, e2, fpi, tpiba2,delta_t
      complex(kind=dp), intent(in) ::  charge_g(:), charge_g_due(:)
      complex(kind=dp), intent(in), optional ::  charge_g_tre(:)
      complex(kind=dp), intent(inout) :: v_mean(:), v_point(:)
      real(kind=DP), allocatable :: fac(:)
      real(kind=dp) :: qq_fact
      integer :: igm
   
   
      allocate (fac(ngm))
   
   ! calculation of fac
   ! fac(G) = e2*fpi/(tpiba2*G^2*omega)
   ! v(G) = charge(G)*fac
   
      if (gstart == 2) fac(1) = 0.d0
      do igm = gstart, ngm
         qq_fact = g(1, igm)**2.d0 + g(2, igm)**2.d0 + g(3, igm)**2.d0
         fac(igm) = (e2*fpi/(tpiba2*qq_fact))
      end do
      fac(:) = fac(:)/omega
   
   
   !-------STEP4----- numerical derivatives of Hartree potentials
   
   !We compute v_point and v_mean
      do igm = 1, ngm
         v_point(igm) = (charge_g(igm) - charge_g_due(igm))*fac(igm)/delta_t ! v(t+dt)-v(t-dt)
      end do
      if (present(charge_g_tre)) then ! set v_mean as v(t)
         do igm = 1, ngm
            v_mean(igm) = charge_g_tre(igm)*fac(igm)
         end do
          
      else ! set v_mean as (v(t)+v(t+dt))/2.0
         do igm = 1, ngm
            v_mean(igm) = (charge_g(igm) + charge_g_due(igm))*fac(igm)/2.d0
         end do
      endif
   
      deallocate (fac)
   
   end subroutine

END MODULE hartree_mod

