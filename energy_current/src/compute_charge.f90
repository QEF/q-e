!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
module compute_charge_mod
   implicit none

contains
   subroutine compute_charge(psic, evc, npw, nbnd, ngm, dffts, charge, charge_g)
      use kinds, only: DP
      USE fft_types, ONLY: fft_type_descriptor
      USE fft_interfaces, ONLY: fwfft, invfft
      use io_global, only: ionode
      USE mp_pools, ONLY: intra_pool_comm
      use mp, only: mp_sum
      implicit none

      complex(kind=dp), intent(inout) :: psic(:)
      complex(kind=dp), intent(inout) :: charge_g(:)
      complex(kind=dp), intent(in) :: evc(:, :)
      real(kind=dp), intent(inout) :: charge(:)
      integer, intent(in) :: npw, nbnd, ngm
      type(fft_type_descriptor), intent(in) ::dffts
      real(kind=dp) :: q_tot
      integer :: iv, i

! charge_g is obtained from a  FFT of |evc(r)|^2, where evc(r)=IFFT(evc)
! Optimization: 2 bands done with a single IFFT

      charge = 0.d0
      do iv = 1, nbnd, 2
         psic = 0.d0
         if (iv == nbnd) then
            psic(dffts%nl(1:npw)) = evc(1:npw, iv)
            psic(dffts%nlm(1:npw)) = CONJG(evc(1:npw, iv))
         else
            psic(dffts%nl(1:npw)) = evc(1:npw, iv) + (0.D0, 1.D0)*evc(1:npw, iv + 1)
            psic(dffts%nlm(1:npw)) = CONJG(evc(1:npw, iv) - (0.D0, 1.D0)*evc(1:npw, iv + 1))
         end if
         call invfft('Wave', psic, dffts)
         charge(1:dffts%nnr) = charge(1:dffts%nnr) + dble(psic(1:dffts%nnr))**2.0
         if (iv /= nbnd) then
            charge(1:dffts%nnr) = charge(1:dffts%nnr) + dimag(psic(1:dffts%nnr))**2.0 !is dimag standard fortran?
         end if
      end do
      q_tot = 0.
      do i = 1, dffts%nnr
         q_tot = q_tot + charge(i)
      end do
      q_tot = q_tot/(dffts%nr1*dffts%nr2*dffts%nr3)
      call mp_sum(q_tot, intra_pool_comm)
      IF (ionode) THEN
         print *, 'check_charge', q_tot
      ENDIF

!We need to multiply by 2 for spin degeneracy
      charge(1:dffts%nnr) = charge(1:dffts%nnr)*2.d0
!

!computation of charge in reciprocal space (FFT of psic)
      psic = 0.d0
      psic(1:dffts%nnr) = dcmplx(charge(1:dffts%nnr), 0.d0)
      call fwfft('Rho', psic, dffts)
      charge_g(1:ngm) = psic(dffts%nl(1:ngm))

   end subroutine

end module
