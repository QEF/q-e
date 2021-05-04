!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
module hartree_xc_mod
   USE kinds, ONLY: DP

contains

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
      USE xc_lib, ONLY: xclib_get_id
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
      real(kind=DP), allocatable ::charge_plus_or_zero(:), exgradcharge_r(:, :), exdotcharge_r(:)
      real(kind=DP) :: update(1:3), update_a(1:3), update_b(1:3)
      real(kind=DP) :: fac, qq_fact

!
      complex(kind=DP), allocatable ::charge_g_plus_or_zero(:), v_point(:), v_mean(:), charge_g_minus_or_zero(:)
!
      complex(kind=DP), allocatable ::exgradcharge_g(:, :)
!
      integer ::icoord, ir
!

      if (ionode) write (stdout, *) 'BEGIN: HARTREE & KOHN'

      if ((xclib_get_id('GGA', 'EXCH') /= 3) .or. (xclib_get_id('GGA', 'CORR') /= 4)) then
         do_xc_curr = .false.
         if ((xclib_get_id('GGA', 'EXCH') /= 0) .or. (xclib_get_id('GGA', 'CORR') /= 0)) then
            call errore('ENERGY CURRENT', 'XC NOT PBE OR LDA. ABORT.', 1)
         end if
      else
         do_xc_curr = .true.
      end if

      npw = npwx ! only gamma
      allocate (charge_g_plus_or_zero(ngm))
      allocate (charge_g_minus_or_zero(ngm))
      allocate (charge_minus_or_zero(dffts%nnr))
      allocate (v_point(ngm))
      allocate (v_mean(ngm))
!
      allocate (charge_plus_or_zero(dffts%nnr))
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

      call compute_charge(psic, scf_all%t_plus%evc, npw, nbnd, ngm, dffts, charge_plus_or_zero, charge_g_plus_or_zero)
      call compute_charge(psic, scf_all%t_minus%evc, npw, nbnd, ngm, dffts, charge_minus_or_zero, charge_g_minus_or_zero)

      if (do_xc_curr) then
         exdotcharge_r(1:dffts%nnr) = (charge_plus_or_zero(1:dffts%nnr) - charge_minus_or_zero(1:dffts%nnr))/delta_t
      end if
      !hartree derivative
      if (gstart == 2) then
         v_point(1) = 0.d0
      end if
      do igm = gstart, ngm
         qq_fact = g(1, igm)**2.d0 + g(2, igm)**2.d0 + g(3, igm)**2.d0
         fac = e2*fpi/(tpiba2*qq_fact*omega)
         v_point(igm) = (charge_g_plus_or_zero(igm) - charge_g_minus_or_zero(igm))*fac/delta_t ! v(t+dt)-v(t-dt)
      end do

      !we computed all the numerical derivatives: charge_plus_or_zero/minus are not needed anymore, we have to compute stuff
      !in t=0.
      !if 2pt derivative is used, t corresponds to scf_all%t_plus,
      !if 3pt derivative is used, t corresponds to scf_all%t_zero. I put the result in charge_g_plus_or_zero and charge_plus_or_zero
      !in any case, after the following two lines, charge_g_plus_or_zero and charge_plus_or_zero contains data related to t
      if (three_point_derivative) &
         call compute_charge(psic, scf_all%t_zero%evc, npw, nbnd, ngm, dffts, charge_plus_or_zero, charge_g_plus_or_zero)
      !hartree potential
      !charge_g_plus_or_zero is charge(r) at time t
      !charge_plus_or_zero is charge(G) at time t
      if (gstart == 2) then
         v_mean(1) = 0.d0
      end if
      do igm = gstart, ngm
         qq_fact = g(1, igm)**2.d0 + g(2, igm)**2.d0 + g(3, igm)**2.d0
         fac = e2*fpi/(tpiba2*qq_fact*omega)
         v_mean(igm) = charge_g_plus_or_zero(igm)*fac
      end do

      if (do_xc_curr) then
         do icoord = 1, 3
            do igm = 1, ngm
               exgradcharge_g(icoord, igm) = charge_g_plus_or_zero(igm)*(0.d0, 1.d0)*g(icoord, igm)*tpiba
            end do

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
         charge_plus_or_zero(1:dffts%nnr) = charge_plus_or_zero(1:dffts%nnr)/omega
         exgradcharge_r(1:3, 1:dffts%nnr) = exgradcharge_r(1:3, 1:dffts%nnr)/omega
!
!
         do ir = 1, dffts%nnr
            if (charge_plus_or_zero(ir) > 1.0E-10_DP) then
               call pbex_current(abs(charge_plus_or_zero(ir)), exgradcharge_r(1:3, ir), 1, update_a(1:3))
               call pbec_current(abs(charge_plus_or_zero(ir)), exgradcharge_r(1:3, ir), 1, update_b(1:3))
               update(1:3) = e2*(update_a(1:3) + update_b(1:3))
               update(1:3) = -update(1:3)*exdotcharge_r(ir)*charge_plus_or_zero(ir)
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
      deallocate (charge_g_plus_or_zero)
      deallocate (charge_g_minus_or_zero)
      deallocate (v_point)
      deallocate (v_mean)
      deallocate (charge_plus_or_zero)
      deallocate (exgradcharge_r)
      deallocate (exgradcharge_g)
      deallocate (exdotcharge_r)

   end subroutine

end module
