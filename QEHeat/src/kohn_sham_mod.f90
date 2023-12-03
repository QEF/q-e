!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
MODULE kohn_sham_mod

   USE kinds, ONLY: DP
   use scf_result_mod, only: scf_result, multiple_scf_result, scf_result_set_global_variables, &
                             scf_result_set_tau_vel_global_variables

contains

   subroutine init_kohn_sham()
      use qpoint, only: ikqs
      use wvfct, only: nbnd, npwx
      implicit none

      if (.not. allocated(ikqs)) allocate (ikqs(1))
      ikqs(1) = 1

   end subroutine



   subroutine current_kohn_sham(J, J_a, J_b, J_el, dt, scf_all, &
                                dvpsi_save, save_dvpsi, &
                                nbnd, npw, npwx, dffts, evc, g, ngm, gstart, &
                                tpiba2, at, vkb, nkb, xk, igk_k, g2kin, et, hpsi_test, &
                                omega, gg, intra_bgrp_comm, &
                                nat, ityp)
      use kinds, only: DP
      use fft_base, only: fft_type_descriptor !,dffts
      use mp, only: mp_sum
      use io_global, only: ionode
      USE constants, ONLY: e2, fpi, pi
      USE fft_interfaces, ONLY: fwfft, invfft
      use becmod
      USE eqv, ONLY: dpsi, dvpsi !TODO: those variables looks like they can be local variables
      USE mp_pools, ONLY: intra_pool_comm
      use compute_charge_mod, only: compute_charge
      use project_mod, only: project
      use extrapolation, only: update_pot
      use wave_gauge, only: compute_dot_evc_parallel_gauge

      use test_h_psi, only: test
      USE uspp_init,            ONLY : init_us_2

      implicit none

      type(multiple_scf_result), intent(in) :: scf_all
      logical, intent(in) :: save_dvpsi, hpsi_test
      complex(dp), intent(inout) :: dvpsi_save(:, :, :)

      INTEGER, intent(in) :: nbnd, npwx, intra_bgrp_comm, &
                             nat, ityp(:)
      INTEGER, intent(inout) :: npw, igk_k(:, :)
      TYPE(fft_type_descriptor), intent(inout) :: dffts
      COMPLEX(DP), intent(inout) ::  evc(:, :)
      INTEGER, intent(in) :: ngm, gstart, nkb
      REAL(DP), intent(inout) :: g2kin(:)
      REAL(DP), intent(out) :: J(3), J_a(3), J_b(3), J_el(3)
      REAL(DP), intent(in) ::  tpiba2, g(:, :), at(:, :), &
                              xk(:, :), et(:, :), dt, gg(:), &
                              omega
      COMPLEX(DP), intent(inout) :: vkb(:, :)

      !character(LEN=20) :: dft_name
      complex(kind=DP), allocatable ::  evp(:, :), tmp(:, :)
      integer ::  ibnd
      real(kind=DP) :: amodulus
      real(DP) :: emme(nbnd, nbnd), kcurrent(3), &
                  kcurrent_a(3), kcurrent_b(3), ecurrent(3), &
                  sa(nbnd, nbnd), ssa(nbnd, nbnd), sb(nbnd, nbnd), ssb(nbnd, nbnd)

      integer  :: jbnd, ig, ipol
      allocate (tmp(npwx, nbnd))
      allocate (evp(npwx, nbnd))

      !---------------------------------KOHN------------------------------------------------
      call start_clock('kohn_current')
      allocate (dpsi(npwx, nbnd))
      allocate (dvpsi(npwx, nbnd))
! init potentials needed to evaluate  H|psi>
      call update_pot()
      call hinit1()
      call init_us_2(npw, igk_k(:, 1), xk(:, 1), vkb)
      call sum_band()
      call allocate_bec_type(nkb, nbnd, becp)
      call calbec(npw, vkb, evc, becp)
      if (hpsi_test) &
         call test
! For preconditioning:
      do ig = 1, npw
         igk_k(ig, 1) = ig
         g2kin(ig) = ((xk(1, 1) + g(1, igk_k(ig, 1)))**2 + &
                      (xk(2, 1) + g(2, igk_k(ig, 1)))**2 + &
                      (xk(3, 1) + g(3, igk_k(ig, 1)))**2)*tpiba2
      end do


      call compute_dot_evc_parallel_gauge(scf_all%t_minus%evc, scf_all%t_zero%evc, scf_all%t_plus%evc,&
                                    evp, nbnd, npw, npwx, gstart)

      evp(:, :) = evp(:, :)/dt

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
         call project(ipol, dvpsi_save, save_dvpsi)
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
      J = 0.d0
      J_a = 0.d0
      J_b = 0.d0
      J_el = 0.d0
      do ipol = 1, 3
!     at(:, ipol) / amodulus is the versor along direction ipol
!      amodulus = sqrt(at(1, ipol)**2 + at(2, ipol)**2 + at(3, ipol)**2)
         amodulus = at(1, ipol)**2 + at(2, ipol)**2 + at(3, ipol)**2
         J(:) = J(:) + 2.d0*at(:, ipol)*real(kcurrent(ipol))/amodulus
         J_a(:) = J_a(:) + 2.d0*at(:, ipol)*real(kcurrent_a(ipol))/amodulus
         J_b(:) = J_b(:) + 2.d0*at(:, ipol)*real(kcurrent_b(ipol))/amodulus
         J_el(:) = J_el(:) + 2.d0*2.d0*at(:, ipol)*real(ecurrent(ipol))/amodulus
!    J(ipol)=J(ipol)+2.d0*real(kcurrent(ipol))
      end do
      call stop_clock('kohn_current')
      call print_clock('kohn_current')
      !npw=npwold
      if (ionode) print *, 'KOHN CURRENT CALCULATED'
!---------------------------------------------------------------------------
      call deallocate_bec_type(becp)
      deallocate (dpsi)
      deallocate (dvpsi)
      deallocate (tmp)
      deallocate (evp)

   end subroutine

END MODULE

