MODULE hartree_mod

   USE kinds, ONLY: DP

   SAVE

   CHARACTER(len=256) :: file_output, trajdir = '' !, init_linear
   real(kind=DP) ::J_kohn(3), J_kohn_a(3), J_kohn_b(3), J_hartree(3), J_xc(3), J_electron(3)
   real(kind=DP), allocatable :: v_cm(:, :)

   real(kind=DP) ::delta_t, ethr_small_step, ethr_big_step
   complex(kind=DP), allocatable :: evc_uno(:, :), evc_due(:, :) ! TODO: maybe save one allocation of npwx*nb
   !evc_uno will be a copy of last evc
   !evc_due will be a copy of evc of the first electrons() call
   complex(kind=DP), allocatable :: dvpsi_save(:,:,:) ! to save the solution of the system between iterations
   logical :: save_dvpsi = .true. ! if true dvpsi_save is allocated and used

   integer :: first_step, last_step, step_mul, step_rem, n_repeat_every_step
   logical :: restart ! if true try to read last calculated step from output and set first_step
   logical :: subtract_cm_vel ! if true do velocity renormalization
   logical :: ec_test ! activates tests for debugging purposes
   logical :: re_init_wfc_1=.false., re_init_wfc_2=.false. ! initialize again evc before scf step number 1 or 2

   contains


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

