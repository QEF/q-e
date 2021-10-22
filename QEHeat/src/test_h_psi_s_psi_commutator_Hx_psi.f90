!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
module test_h_psi
   ! this module is a small benchmark that is activated only id hpsi_test
   ! is provided in the input namelist (not documented in the user input manual)
   ! it detects if somwthing is messing up stuff between calls of h_psi, s_psi or
   ! commutator_Hx_psi, checking that they give the same result when called
   ! in different parts of the code
   use kinds, only: dp
   USE wvfct, ONLY: npw, npwx, nbnd
   use klist, only: xk, igk_k
   USE becmod, ONLY: bec_type, becp, calbec, &
                     allocate_bec_type, deallocate_bec_type
   USE uspp, ONLY: nkb, vkb
   use cell_base, only: omega, at
   use mp_bands, only: intra_bgrp_comm
   use gvect, only: ngm, gg, g
   USE ions_base, ONLY: nat, nsp, ityp, tau
   USE uspp_init,            ONLY : init_us_2

   implicit none

   complex(dp), allocatable :: hpsi(:, :), spsi(:, :), chxpsi(:, :), & ! reference results
                               input(:, :) ! input vector

contains

   subroutine init_test(input_vector)

      implicit none
      complex(dp), intent(in) :: input_vector(:, :)
      type(bec_type) ::becp2
      complex(dp), allocatable :: tmp(:, :) ! temporary working vector

      if (.not. allocated(hpsi)) then
         allocate (hpsi(npwx, nbnd))
         allocate (spsi(npwx, nbnd))
         allocate (chxpsi(npwx, nbnd))
         allocate (input(npwx, nbnd))
      end if
      allocate (tmp(npwx, nbnd))

      input = input_vector
      hpsi = (0.d0, 0.d0)
      spsi = (0.d0, 0.d0)
      chxpsi = (0.d0, 0.d0)
      npw = npwx
      call init_us_1(nat, ityp, omega, ngm, g, gg, intra_bgrp_comm)
      call init_us_2(npw, igk_k(:, 1), xk(1, 1), vkb)

      call allocate_bec_type(nkb, nbnd, becp)
      call calbec(npw, vkb, input, becp)
      call h_psi(npwx, npw, nbnd, input, hpsi)
      call s_psi(npwx, npw, nbnd, input, spsi)

      call allocate_bec_type(nkb, nbnd, becp2)
      call commutator_Hx_psi(1, nbnd, at(:,1), becp, becp2, chxpsi)
      call deallocate_bec_type(becp2)

      call deallocate_bec_type(becp)
      deallocate (tmp)
   end subroutine

   subroutine test

      implicit none
      type(bec_type) ::becp2
      complex(dp), allocatable :: tmp(:, :), res(:, :) ! temporary working vector

      allocate (res(npwx, nbnd))
      allocate (tmp(npwx, nbnd))

      npw = npwx

      res = (0.d0, 0.d0)
      call h_psi(npwx, npw, nbnd, input, res)
      if (.not. test_equal(res, hpsi)) &
         write (*, *) 'ERROR: test failed for h_psi'
      res = (0.d0, 0.d0)
      call s_psi(npwx, npw, nbnd, input, res)
      if (.not. test_equal(res, spsi)) &
         write (*, *) 'ERROR: test failed for s_psi'

      call allocate_bec_type(nkb, nbnd, becp2)
      res = (0.d0, 0.d0)
      call commutator_Hx_psi(1, nbnd, at(:,1), becp, becp2, res)
      if (.not. test_equal(res, chxpsi)) &
         write (*, *) 'ERROR: test failed for commutator_Hx_psi'
      call deallocate_bec_type(becp2)

      deallocate (tmp)
      deallocate (res)

   end subroutine

   function test_equal(a, b) result(res)
      implicit none
      logical :: res
      complex(dp), intent(in) :: a(:, :), b(:, :)

      complex(dp) :: m, ma, mb

      m = cmplx(maxval(real(a - b)), maxval(aimag(a - b)), dp)
      ma = cmplx(maxval(real(a)), maxval(aimag(a)), dp)
      mb = cmplx(maxval(real(b)), maxval(aimag(b)), dp)

      write (*, *) m, m/ma, m/mb, ma, mb
      res = .true.
      if (abs(m) > 1e-10) res = .false.
   end function

end module
