module hartree_xc_mod
   USE kinds, ONLY: DP

contains

   subroutine compute_hartree_vpoint_vmean(v_mean, v_point, g, gstart, omega, ngm, e2, fpi, tpiba2, &
                                           delta_t, charge_g_plus, charge_minus, charge_zero)
      implicit none
      integer, intent(in) :: gstart, ngm
      real(kind=dp), intent(in) :: g(:, :), omega, e2, fpi, tpiba2, delta_t
      complex(kind=dp), intent(in) ::  charge_g_plus(:), charge_minus(:)
      complex(kind=dp), intent(in), optional ::  charge_zero(:)
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

      ! numerical derivatives of Hartree potentials

      !We compute v_point and v_mean
      do igm = 1, ngm
         v_point(igm) = (charge_g_plus(igm) - charge_minus(igm))*fac(igm)/delta_t ! v(t+dt)-v(t-dt)
      end do
      if (present(charge_zero)) then ! set v_mean as v(t)
         do igm = 1, ngm
            v_mean(igm) = charge_zero(igm)*fac(igm)
         end do
      else
         do igm = 1, ngm
            !v_mean(igm) = (charge_g_plus(igm) + charge_minus(igm))*fac(igm)/2.d0
            v_mean(igm) = charge_minus(igm)*fac(igm)/2.d0
         end do
      endif

      deallocate (fac)

   end subroutine

   subroutine current_hartree_xc(three_point_derivative, delta_t, scf_all, &
                                 j_hartree, j_xc, nbnd, npw, npwx, dffts, psic, g, ngm, gstart, &
                                 tpiba, omega, tpiba2)
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
!   use becmod
      USE mp_pools, ONLY: intra_pool_comm
      USE funct, ONLY: get_igcx, get_igcc
      use compute_charge_mod, only: compute_charge
      use scf_result_mod, only: scf_result, multiple_scf_result

      implicit none

      INTEGER, intent(in) :: nbnd, npwx
      INTEGER, intent(inout) :: npw
      TYPE(fft_type_descriptor), intent(inout) :: dffts
      COMPLEX(DP), intent(inout) :: psic(:)
      INTEGER, intent(in) :: ngm, gstart
      REAL(DP), intent(in) :: tpiba, omega, tpiba2, delta_t
      REAL(DP), intent(inout) ::  g(:, :)
      real(dp), intent(out) :: j_hartree(3), j_xc(3)
      type(multiple_scf_result), intent(in) :: scf_all
      logical, intent(in) :: three_point_derivative

      integer :: igm
      logical :: do_xc_curr
      real(kind=DP), allocatable ::charge_minus_or_zero(:)
      real(kind=DP), allocatable ::charge_plus(:), exgradcharge_r(:, :), exdotcharge_r(:)
      real(kind=DP) :: update(1:3), update_a(1:3), update_b(1:3)

!
      complex(kind=DP), allocatable ::charge_g_plus(:), v_point(:), v_mean(:), charge_g_minus(:), &
                                       charge_g_zero(:)
!
      complex(kind=DP), allocatable ::exgradcharge_g(:, :)
!
      integer ::icoord, ir
!

      if (ionode) write (stdout, *) 'BEGIN: HARTREE & KOHN'

      if ((get_igcx() /= 3) .or. (get_igcc() /= 4)) then
         do_xc_curr = .false.
         if ((get_igcx() /= 0) .or. (get_igcc() /= 0)) then
            call errore('ENERGY CURRENT', 'XC NOT PBE OR LDA. ABORT.', 1)
         end if
      else
         do_xc_curr = .true.
      end if

      npw = npwx ! only gamma
      allocate (charge_g_plus(ngm))
      allocate (charge_g_minus(ngm))
      if (three_point_derivative) &
         allocate (charge_g_zero(ngm))
      allocate (charge_minus_or_zero(dffts%nnr))
      allocate (v_point(ngm))
      allocate (v_mean(ngm))
!
      allocate (charge_plus(dffts%nnr))
      allocate (exgradcharge_r(3, dffts%nnr))
      allocate (exgradcharge_g(3, ngm))
      allocate (exdotcharge_r(dffts%nnr))
      call start_clock('hartree_current')

!-------------------------------HARTREE---------------------------------------------------------
! Calculation of Hartree, exchange and Kohn-Sham currents require that there is access in memory
! simultaneously to the wavefunctions evaluated at two (or three) different time steps (t-Dt), t and t+Dt
! , where Dt is the time step used for the numerical derivative (if 3 points are used with a factor 1/2).
! wavefunctions are stored in the scf_all type. scf_all%t_minus,t_zero,t_plus have are types that store
! wfc and position and velocities of atoms that were used to compute the wfcs.
! if only 2 points are used, t_minus == t_zero


      call compute_charge(psic, scf_all%t_plus%evc, npw, nbnd, ngm, dffts, charge_plus, charge_g_plus)
      call compute_charge(psic, scf_all%t_minus%evc, npw, nbnd, ngm, dffts, charge_minus_or_zero, charge_g_minus)

      if (do_xc_curr) then
         exdotcharge_r(1:dffts%nnr) = (charge_plus(1:dffts%nnr) - charge_minus_or_zero(1:dffts%nnr))/delta_t
      end if
      !

!-------STEP3----- computation of Hartree potentials from the charges just computed.
!-------STEP4-----  and numerical derivatives of Hartree potentials
!the difference with three_point_derivative is that I use the potential at t_zero as v_mean
!in place of the average between t_minus
      if (three_point_derivative) then
         !use charge_minus_or_zero!
         call compute_charge(psic, scf_all%t_zero%evc, npw, nbnd, ngm, dffts, charge_minus_or_zero, charge_g_zero)
         call compute_hartree_vpoint_vmean(v_mean, v_point, g, gstart, omega, ngm, e2, fpi, tpiba2, delta_t, &
                                           charge_g_plus, charge_g_minus, charge_g_zero)
      else
         call compute_hartree_vpoint_vmean(v_mean, v_point, g, gstart, omega, ngm, e2, fpi, tpiba2,&
                 delta_t, charge_g_plus, charge_g_minus)
      endif
      !now charge_minus_or_zero is charge at t-Dt if two point derivative is used; otherwise is charge at time t!

      if (do_xc_curr) then
         do icoord = 1, 3
            if (three_point_derivative) then
            do igm = 1, ngm
               exgradcharge_g(icoord, igm) = charge_g_zero(igm)*(0.d0, 1.d0)*g(icoord, igm)*tpiba
            end do
    else
            do igm = 1, ngm
               exgradcharge_g(icoord, igm) = charge_g_minus(igm)*(0.d0, 1.d0)*g(icoord, igm)*tpiba
            end do

    end if

            psic = 0.d0
            psic(dffts%nl(1:ngm)) = exgradcharge_g(icoord, 1:ngm)
            psic(dffts%nlm(1:ngm)) = CONJG(exgradcharge_g(icoord, 1:ngm))
            call invfft('Rho', psic, dffts)
            exgradcharge_r(icoord, 1:dffts%nnr) = dble(psic(1:dffts%nnr))
         end do
      end if
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
         !note that if three_point_derivative == true, charge_minus_or_zero is charge at time t
         charge_minus_or_zero(1:dffts%nnr) = charge_minus_or_zero(1:dffts%nnr)/omega
         exgradcharge_r(1:3, 1:dffts%nnr) = exgradcharge_r(1:3, 1:dffts%nnr)/omega
!
!
         do ir = 1, dffts%nnr
            if (charge_minus_or_zero(ir) > 1.0E-10_DP) then
               call pbex_current(abs(charge_minus_or_zero(ir)), exgradcharge_r(1:3, ir), 1, update_a(1:3))
               call pbec_current(abs(charge_minus_or_zero(ir)), exgradcharge_r(1:3, ir), 1, update_b(1:3))
               update(1:3) = e2*(update_a(1:3) + update_b(1:3))
               update(1:3) = -update(1:3)*exdotcharge_r(ir)*charge_minus_or_zero(ir)
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
      deallocate (charge_minus_or_zero)
      deallocate (charge_g_plus)
      deallocate (charge_g_minus)
      if (allocated(charge_g_zero)) deallocate (charge_g_zero)
      deallocate (v_point)
      deallocate (v_mean)
      deallocate (charge_plus)
      deallocate (exgradcharge_r)
      deallocate (exgradcharge_g)
      deallocate (exdotcharge_r)

   end subroutine

end module
