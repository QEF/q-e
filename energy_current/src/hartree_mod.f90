MODULE hartree_mod

   USE kinds, ONLY: DP
   use averages, only : online_average
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
   type(online_average) :: ave_cur
   real(kind=DP), allocatable :: v_cm(:, :)

   real(kind=DP) ::delta_t, ethr_small_step, ethr_big_step
   type(multiple_scf_result) :: scf_all

   complex(kind=DP), allocatable :: dvpsi_save(:,:,:) ! to save the solution of the system between iterations
   logical :: save_dvpsi = .true. ! if true dvpsi_save is allocated and used

   integer :: first_step, last_step, step_mul, step_rem, n_repeat_every_step
   logical :: restart ! if true try to read last calculated step from output and set first_step
   logical :: subtract_cm_vel ! if true do velocity renormalization
   logical :: ec_test ! activates tests for debugging purposes
   logical :: re_init_wfc_1=.false., re_init_wfc_2=.false. ! initialize again evc before scf step number 1 or 2
   logical :: re_init_wfc_3=.false. ! initialize again evc before scf step number 1 or 2
   logical :: three_point_derivative=.true. ! compute hartree derivative with 3 points


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

subroutine init_hartree()
   use qpoint, only: ikqs
   use wvfct, only: nbnd, npwx
   implicit none

   if (.not. allocated(ikqs)) allocate (ikqs(1))
   if (save_dvpsi) then
       if (.not. allocated(dvpsi_save)) then
               allocate (dvpsi_save(npwx, nbnd,3))
               dvpsi_save = (0.d0, 0.d0)
       end if
   end if


   ikqs(1) = 1

end subroutine



subroutine current_hartree_xc(nbnd, npw, npwx, dffts, psic, g, ngm, gstart, &
                tpiba, omega, tpiba2)
   use kinds, only: DP
   !use wvfct, only: nbnd, npw, npwx
   !use wavefunctions, only: psic, evc
   !use gvect, only: g, ngm, gstart
   !USE cell_base, ONLY: tpiba, omega, tpiba2, alat, at
   !use ions_base, only: tau
   !use uspp, ONLY: vkb, nkb
   !use klist, only: xk, igk_k
   !use wvfct, ONLY: g2kin, et
   use fft_base, only: fft_type_descriptor !,dffts
   use mp, only: mp_sum
   use io_global, only: stdout, ionode
   USE constants, ONLY: e2, fpi, pi
   USE fft_interfaces, ONLY: fwfft, invfft
   use becmod
   USE mp_pools, ONLY: intra_pool_comm
   USE funct, ONLY : get_igcx, get_igcc
   use compute_charge_mod, only : compute_charge

   implicit none

   INTEGER, intent(in) :: nbnd,  npwx
   INTEGER, intent(inout) :: npw
   TYPE ( fft_type_descriptor ),intent(inout) :: dffts
   COMPLEX(DP), intent(inout) :: psic(:)
   INTEGER, intent(in) :: ngm, gstart
   REAL(DP), intent(in) :: tpiba, omega, tpiba2
   REAL(DP), intent(inout) ::  g(:,:)

   complex(kind=DP), allocatable ::  evp(:, :), tmp(:, :)
   integer :: igm
   logical :: do_xc_curr
   real(kind=DP), allocatable ::charge(:)
   real(kind=DP), allocatable ::excharge_r(:), exgradcharge_r(:, :), exdotcharge_r(:)
   real(kind=DP) :: update(1:3), update_a(1:3), update_b(1:3)

!
   complex(kind=DP), allocatable ::charge_g(:),  v_point(:), v_mean(:), charge_g_due(:), &
           charge_g_tre(:)
!
   complex(kind=DP), allocatable ::exgradcharge_g(:, :)
!
   integer ::icoord, ir
!
!   logical  :: l_test
 !  real(DP) :: s(nbnd, nbnd), ss(nbnd, nbnd), &
 !  &
 !              sa(nbnd, nbnd), ssa(nbnd, nbnd), sb(nbnd, nbnd), ssb(nbnd,nbnd)


   integer, external :: find_free_unit

   if (ionode) write (stdout, *) 'BEGIN: HARTREE & KOHN'

   if ( (get_igcx() /= 3) .or. (get_igcc() /=4 ) ) then
           do_xc_curr=.false.
           if ( (get_igcx() /= 0) .or. (get_igcc() /=0 ) ) then
                call errore('ENERGY CURRENT', 'XC NOT PBE OR LDA. ABORT.', 1)
           end if
   else
           do_xc_curr=.true.
   end if

   npw = npwx ! only gamma
   allocate (tmp(npwx, nbnd))
   allocate (evp(npwx, nbnd))
   allocate (charge_g(ngm))
   allocate (charge_g_due(ngm))
   if (three_point_derivative) &
       allocate (charge_g_tre(ngm))
   allocate (charge(dffts%nnr))
   allocate (v_point(ngm))
   allocate (v_mean(ngm))
!
   allocate (excharge_r(dffts%nnr))
   allocate (exgradcharge_r(3, dffts%nnr))
   allocate (exgradcharge_g(3, ngm))
   allocate (exdotcharge_r(dffts%nnr))
   call start_clock('hartree_current')
   call init_hartree()

!-------------------------------HARTREE---------------------------------------------------------
! Calculation of Hartree, exchange and Kohn-Sham currents require that there is access in memory
! simultaneously to the wavefunctions evaluated at two (or three) different time steps (t-Dt), t and t+Dt
! , where Dt is the time step used for the numerical derivative (if 3 points are used with a factor 1/2).
! wavefunctions are stored in the scf_all type. scf_all%t_minus,t_zero,t_plus have are types that store
! wfc and position and velocities of atoms that were used to compute the wfcs.
! if only 2 points are used, t_minus == t_zero

!-------STEP1: reading and allocation (npwx = number of plane waves (npwx>npw), nbnd = number of bands (n_electrons/2 for insulators))

!-------STEP2.1: calculation of charge_g, the charge in reciprocal space at time t+Dt.

call compute_charge(psic, scf_all%t_plus%evc, npw, nbnd, ngm, dffts, charge, charge_g)

!!!!!!!!!!!!------------ Saving some quantities for the evaluation of the Exchange-correlation current 1/2  -------------!!!!!!!!!!!!!!!!!!

   if (do_xc_curr) then
! copying of charge in real space
     excharge_r(1:dffts%nnr) = charge(1:dffts%nnr)
     do icoord = 1, 3
!
       do igm = 1, ngm
          exgradcharge_g(icoord, igm) = charge_g(igm)*(0.d0, 1.d0)*g(icoord, igm)*tpiba
       end do
!
       psic = 0.d0
       psic(dffts%nl(1:ngm)) = exgradcharge_g(icoord, 1:ngm)
       psic(dffts%nlm(1:ngm)) = CONJG(exgradcharge_g(icoord, 1:ngm))
!
       call invfft('Rho', psic, dffts)
!
       exgradcharge_r(icoord, 1:dffts%nnr) = dble(psic(1:dffts%nnr))
     end do
   end if
!
!!!!!!!!!!!!------------- end saving quantities for XC current----------------- !!!!!!!!!!!!!!!!!!!!


!-------STEP2.2-------inizializing chargeg_due, charge at tempo t-Dt.
call compute_charge(psic, scf_all%t_minus%evc, npw, nbnd, ngm, dffts, charge, charge_g_due)

!!!!!!!!!!!!------------Saving quantities for XC current 2/2 -------------!!!!!!!!!!!!!!!!!!
!
   if (do_xc_curr) then
       exdotcharge_r(1:dffts%nnr) = (excharge_r(1:dffts%nnr) - charge(1:dffts%nnr))/delta_t
   end if
       !
!!!!!!!!!!!!------------- End saving quantities for XC current----------------- !!!!!!!!!!!!!!!!!!!!


!-------STEP3----- computation of Hartree potentials from the charges just computed.
!-------STEP4-----  and numerical derivatives of Hartree potentials
!the difference with three_point_derivative is that I use the potential at t_zero as v_mean
!in place of the average between t_plus and t_minus
if (three_point_derivative) then
    call compute_charge(psic, scf_all%t_zero%evc, npw, nbnd, ngm, dffts, charge, charge_g_tre)
    call compute_hartree_vpoint_vmean(v_mean, v_point, g, gstart, omega, ngm, e2, fpi, tpiba2, delta_t,&
            charge_g, charge_g_due, charge_g_tre)
else
    call compute_hartree_vpoint_vmean(v_mean, v_point, g, gstart, omega, ngm, e2, fpi, tpiba2, delta_t, charge_g, charge_g_due)
endif

!-------STEP 5----- Application of final formula.

   J_hartree = 0.d0
   do igm = gstart, ngm
      J_hartree(1:3) = J_hartree(1:3) - omega/(4.d0*pi)*dimag(CONJG(v_point(igm))*v_mean(igm))*g(1:3, igm)*tpiba
   end do
   call mp_sum(J_hartree, intra_pool_comm)

   call stop_clock('hartree_current')
   call print_clock('hartree_current')
   if (ionode) print *, 'HARTREE CURRENT CALCULATED'

!-----------------EXCHANGE-CORRELATION-----------------------------------------------

   call start_clock('xc_current')

   J_xc(1:3) = 0.d0

   if (do_xc_curr) then
       exdotcharge_r(1:dffts%nnr) = exdotcharge_r(1:dffts%nnr)/omega
       excharge_r(1:dffts%nnr) = excharge_r(1:dffts%nnr)/omega
       exgradcharge_r(1:3, 1:dffts%nnr) = exgradcharge_r(1:3, 1:dffts%nnr)/omega
!
!
       do ir = 1, dffts%nnr
          if (excharge_r(ir) > 1.0E-10_DP) then
             call pbex_current(abs(excharge_r(ir)), exgradcharge_r(1:3, ir), 1, update_a(1:3))
             call pbec_current(abs(excharge_r(ir)), exgradcharge_r(1:3, ir), 1, update_b(1:3))
             update(1:3) = e2*(update_a(1:3) + update_b(1:3))
             update(1:3) = -update(1:3)*exdotcharge_r(ir)*excharge_r(ir)
             J_xc(1:3) = J_xc(1:3) + update(1:3)
          end if
       end do
!
       call mp_sum(J_xc, intra_pool_comm)
!
!Volume element
!
       J_xc(1:3) = J_xc(1:3)*omega/(dffts%nr1*dffts%nr2*dffts%nr3)
!
   end if

   call stop_clock('xc_current')
   call print_clock('xc_current')
   if (ionode) print *, 'X-C CURRENT CALCULATED'
   deallocate (charge)
   deallocate (charge_g)
   deallocate (charge_g_due)
   if (allocated(charge_g_tre)) deallocate (charge_g_tre)
   deallocate (v_point)
   deallocate (v_mean)
   call deallocate_bec_type(becp)
   deallocate (excharge_r)
   deallocate (exgradcharge_r)
   deallocate (exgradcharge_g)
   deallocate (exdotcharge_r)


   end subroutine


   subroutine current_kohn_sham(nbnd, npw, npwx, dffts, evc, g, ngm, gstart, &
                tpiba, tpiba2,  at, vkb, nkb, xk, igk_k, g2kin, et)
   use kinds, only: DP
   !use wvfct, only: nbnd, npw, npwx
   !use wavefunctions, only: psic, evc
   !use gvect, only: g, ngm, gstart
   !USE cell_base, ONLY: tpiba, omega, tpiba2, alat, at
   !use ions_base, only: tau
   !use uspp, ONLY: vkb, nkb
   !use klist, only: xk, igk_k
   !use wvfct, ONLY: g2kin, et
   use fft_base, only: fft_type_descriptor !,dffts
   use mp, only: mp_sum
   use io_global, only: ionode
   USE constants, ONLY: e2, fpi, pi
   USE fft_interfaces, ONLY: fwfft, invfft
   use becmod
   USE eqv, ONLY: dpsi, dvpsi !TODO: those variables looks like they can be local variables
   USE mp_pools, ONLY: intra_pool_comm
   USE funct, ONLY : get_igcx, get_igcc
   use compute_charge_mod, only : compute_charge

   implicit none

   INTEGER, intent(in) :: nbnd,  npwx
   INTEGER, intent(inout) :: npw, igk_k(:,:)
   TYPE ( fft_type_descriptor ),intent(inout) :: dffts
   COMPLEX(DP), intent(inout) ::  evc(:,:)
   INTEGER, intent(in) :: ngm, gstart, nkb
   REAL(DP), intent(in) :: tpiba, tpiba2
   REAL(DP), intent(inout) ::  g(:,:), g2kin(:), at(:,:), &
                           xk(:,:), et(:,:)
   COMPLEX(DP), intent(inout) :: vkb(:,:)

   !character(LEN=20) :: dft_name
   complex(kind=DP), allocatable ::  evp(:, :), tmp(:, :)
   integer ::  ibnd
   real(kind=DP), allocatable ::charge(:)
   real(kind=DP), allocatable ::excharge_r(:), exgradcharge_r(:, :), exdotcharge_r(:)


   real(kind=DP) :: amodulus
!
   complex(kind=DP), allocatable ::charge_g(:),  v_point(:), v_mean(:), charge_g_due(:), &
           charge_g_tre(:)
!
   complex(kind=DP), allocatable ::exgradcharge_g(:, :)
!


!
!   logical  :: l_test
   real(DP) :: emme(nbnd, nbnd), kcurrent(3),  &
               kcurrent_a(3), kcurrent_b(3), ecurrent(3), &
               sa(nbnd, nbnd), ssa(nbnd, nbnd), sb(nbnd, nbnd), ssb(nbnd,nbnd)

   integer  :: jbnd, ig, ipol
   integer, external :: find_free_unit


   allocate (tmp(npwx, nbnd))
   allocate (evp(npwx, nbnd))
   allocate (charge_g(ngm))
   allocate (charge_g_due(ngm))
   if (three_point_derivative) &
       allocate (charge_g_tre(ngm))
   allocate (charge(dffts%nnr))
   allocate (v_point(ngm))
   allocate (v_mean(ngm))
!
   allocate (excharge_r(dffts%nnr))
   allocate (exgradcharge_r(3, dffts%nnr))
   allocate (exgradcharge_g(3, ngm))
   allocate (exdotcharge_r(dffts%nnr))


   !---------------------------------KOHN------------------------------------------------
   call start_clock('kohn_current')
   call scf_result_set_global_variables(scf_all%t_zero) ! this set evc, tau and vel from saved values
   allocate (dpsi(npwx, nbnd))
   allocate (dvpsi(npwx, nbnd))
! For preconditioning:
   do ig = 1, npw
      igk_k(ig, 1) = ig
      g2kin(ig) = ((xk(1, 1) + g(1, igk_k(ig, 1)))**2 + &
                   (xk(2, 1) + g(2, igk_k(ig, 1)))**2 + &
                   (xk(3, 1) + g(3, igk_k(ig, 1)))**2)*tpiba2
   end do
! init potentials needed to evaluate  H|psi>
   call init_us_1()
   call init_us_2(npw, igk_k(1, 1), xk(1, 1), vkb)
   call allocate_bec_type(nkb, nbnd, becp)
   call calbec(npw, vkb, evc, becp)

!   l_test = .true.
!   if (.not. l_test) then
!
!      s = 0.d0
!! computed s = < evc_due, evc_uno >, remove contribution at G=0
!      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, scf_all%t_plus%evc, 2*npwx, scf_all%t_minus%evc, 2*npwx, 0.d0, s, nbnd)
!      if (gstart == 2) then
!         do ibnd = 1, nbnd
!            do jbnd = 1, nbnd
!               s(ibnd, jbnd) = s(ibnd, jbnd) - dble(conjg(scf_all%t_plus%evc(1, ibnd))*scf_all%t_minus%evc(1, jbnd))
!            end do
!         end do
!      end if
!      call mp_sum(s, intra_pool_comm)
!! computes phi_v^c_punto using <s,s>
!      call dgemm('T', 'N', nbnd, nbnd, nbnd, 1.d0, s, nbnd, s, nbnd, 0.d0, ss, nbnd)
!      evp = 0.d0
!      do ibnd = 1, nbnd
!         do jbnd = 1, nbnd
!            do ig = 1, npw
!               evp(ig, ibnd) = evp(ig, ibnd) - scf_all%t_plus%evc(ig, jbnd)*s(jbnd, ibnd)
!               evp(ig, ibnd) = evp(ig, ibnd) + scf_all%t_minus%evc(ig, jbnd)*ss(ibnd, jbnd)
!            end do
!         end do
!      end do
!      evp(:, :) = evp(:, :)/delta_t
!
!   else

      sa = 0.d0
      sb = 0.d0

! computed sb = < evc_due, evc_uno >, sa = <evc_uno, evc_uno> remove contribution at G=0
! sb is the old s
! sb is a multiple of the identity if 2 points are used (t_minus==t_zero)

      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, scf_all%t_plus%evc, 2*npwx, scf_all%t_zero%evc, 2*npwx, 0.d0, sa, nbnd)
      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, scf_all%t_minus%evc, 2*npwx, scf_all%t_zero%evc, 2*npwx, 0.d0, sb, nbnd)
      if (gstart == 2) then
         do ibnd = 1, nbnd
            do jbnd = 1, nbnd
               sa(ibnd, jbnd) = sa(ibnd, jbnd) - dble(conjg(scf_all%t_plus%evc(1, ibnd))*scf_all%t_zero%evc(1, jbnd))
               sb(ibnd, jbnd) = sb(ibnd, jbnd) - dble(conjg(scf_all%t_minus%evc(1, ibnd))*scf_all%t_zero%evc(1, jbnd))
            end do
         end do
      end if
      call mp_sum(sa, intra_pool_comm)
      call mp_sum(sb, intra_pool_comm)
! computes phi_v^c_punto using <s,s>

      call dgemm('T', 'N', nbnd, nbnd, nbnd, 1.d0, sa, nbnd, sa, nbnd, 0.d0, ssa, nbnd)
      call dgemm('T', 'N', nbnd, nbnd, nbnd, 1.d0, sb, nbnd, sb, nbnd, 0.d0, ssb, nbnd)

      evp = 0.d0
      do ibnd = 1, nbnd
         do jbnd = 1, nbnd
            do ig = 1, npw

               evp(ig, ibnd) = evp(ig, ibnd) + scf_all%t_plus%evc(ig, jbnd)*sa(jbnd, ibnd)
               evp(ig, ibnd) = evp(ig, ibnd) - scf_all%t_zero%evc(ig, jbnd)*ssa(ibnd, jbnd)

               evp(ig, ibnd) = evp(ig, ibnd) - scf_all%t_minus%evc(ig, jbnd)*sb(jbnd, ibnd)
               evp(ig, ibnd) = evp(ig, ibnd) + scf_all%t_zero%evc(ig, jbnd)*ssb(ibnd, jbnd)

            end do
         end do
      end do
      evp(:, :) = evp(:, :)/delta_t

!   end if

   tmp(:, :) = (0.d0, 0.d0)
   !does H|evp>  i.e. H|phi_v^c_punto> and saves it into tmp:
   call h_psi(npwx, npw, nbnd, evp, tmp)
   kcurrent = 0.d0
   kcurrent_a = 0.d0
   kcurrent_b = 0.d0
   ecurrent = 0.d0
   polariz: do ipol = 1, 3
      call start_clock('project')
      ! computes projection
      call project(ipol)
      call stop_clock('project')
      call print_clock('project')

      call start_clock('kohn-first')
      emme = 0.d0

      ! scalar product
      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, dvpsi, 2*npwx, tmp, 2*npwx, 0.d0, emme, nbnd)
      if (gstart == 2) then
         do ibnd = 1, nbnd
            do jbnd = 1, nbnd
               emme(ibnd, jbnd) = emme(ibnd, jbnd) - dble(conjg(dvpsi(1, ibnd))*tmp(1, jbnd))
            end do
         end do
      end if
      call mp_sum(emme, intra_pool_comm)

      do ibnd = 1, nbnd
         kcurrent(ipol) = kcurrent(ipol) + emme(ibnd, ibnd)
         kcurrent_a(ipol) = kcurrent_a(ipol) + emme(ibnd, ibnd)
      end do
      call stop_clock('kohn-first')
      call print_clock('kohn-first')

      call start_clock('kohn-second')
      emme = 0.d0
      call dgemm('T', 'N', nbnd, nbnd, 2*npw, 2.d0, evp, 2*npwx, dvpsi, 2*npwx, 0.d0, emme, nbnd)
      if (gstart == 2) then
         do ibnd = 1, nbnd
            do jbnd = 1, nbnd
               emme(ibnd, jbnd) = emme(ibnd, jbnd) - dble(conjg(evp(1, ibnd))*dvpsi(1, jbnd))
            end do
         end do
      end if
      call mp_sum(emme, intra_pool_comm)

      ! et containes the eigvenvalues
      do ibnd = 1, nbnd
         kcurrent(ipol) = kcurrent(ipol) + et(ibnd, 1)*emme(ibnd, ibnd)
         kcurrent_b(ipol) = kcurrent_b(ipol) + et(ibnd, 1)*emme(ibnd, ibnd)
         ecurrent(ipol) = ecurrent(ipol) + emme(ibnd, ibnd)
      end do

      call stop_clock('kohn-second')
      call print_clock('kohn-second')

      if (ionode) write (*, "('  KOHN POLARIZATION ',I3,' COMPLETED')") ipol

   end do polariz
   J_kohn = 0.d0
   J_kohn_a = 0.d0
   J_kohn_b = 0.d0
   J_electron = 0.d0
   do ipol = 1, 3
!     at(:, ipol) / amodulus is the versor along direction ipol
!      amodulus = sqrt(at(1, ipol)**2 + at(2, ipol)**2 + at(3, ipol)**2)
      amodulus = at(1, ipol)**2 + at(2, ipol)**2 + at(3, ipol)**2
      J_kohn(:) = J_kohn(:) + 2.d0*at(:, ipol)*real(kcurrent(ipol))/amodulus
      J_kohn_a(:) = J_kohn_a(:) + 2.d0*at(:, ipol)*real(kcurrent_a(ipol))/amodulus
      J_kohn_b(:) = J_kohn_b(:) + 2.d0*at(:, ipol)*real(kcurrent_b(ipol))/amodulus
      J_electron(:) = J_electron(:) + 2.d0*2.d0*at(:, ipol)*real(ecurrent(ipol))/amodulus
!    J_kohn(ipol)=J_kohn(ipol)+2.d0*real(kcurrent(ipol))
   end do
   call stop_clock('kohn_current')
   call print_clock('kohn_current')
   !npw=npwold
   if (ionode) print *, 'KOHN CURRENT CALCULATED'
!---------------------------------------------------------------------------
   deallocate (charge)
   deallocate (charge_g)
   deallocate (charge_g_due)
   if (allocated(charge_g_tre)) deallocate (charge_g_tre)
   deallocate (v_point)
   deallocate (v_mean)
   call deallocate_bec_type(becp)
   deallocate (dpsi)
   deallocate (dvpsi)
   deallocate (excharge_r)
   deallocate (exgradcharge_r)
   deallocate (exgradcharge_g)
   deallocate (exdotcharge_r)

end subroutine



END MODULE hartree_mod

