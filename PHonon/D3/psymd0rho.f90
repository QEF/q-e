!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine psymd0rho (nper, irr, dvtosym)
  !-----------------------------------------------------------------------
  !  p-symmetrize the charge density.
  !

#ifdef __MPI
  USE kinds, ONLY : DP
  USE ions_base, ONLY : nat
  USE symm_base, ONLY : s, ftau
  USE pwcom
  USE phcom
  USE d3com
  USE mp_global, ONLY : me_pool
  USE fft_base,  ONLY : dfftp
  USE scatter_mod,  ONLY : cgather_sym

  USE lr_symm_base, ONLY : irgq

  IMPLICIT NONE

  integer :: nper, irr
  ! the number of perturbations
  ! the representation under consideration

  complex (DP) :: dvtosym (dfftp%nnr, nper)
  ! the potential to symmetrize

  ! local variables

  integer :: i, iper, npp0
  complex (DP),pointer  :: ddvtosym (:,:)
  ! the potential to symmetrize

!  if (nsymq.eq.1.and. (.not.minus_q) ) return

  call start_clock ('psymd0rho')

  allocate ( ddvtosym(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, nper))
  npp0 = 0
  do i = 1, me_pool
     npp0 = npp0 + dfftp%npp (i)
  enddo

  npp0 = npp0 * dfftp%nnp + 1
  do iper = 1, nper
     call cgather_sym (dfftp, dvtosym (:, iper), ddvtosym (:, iper) )
  enddo

  call symd0rho (npertx, nper, irr, ddvtosym, s, ftau, nsymg0, irgq, tg0, &
       nat, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x)
  do iper = 1, nper
     call zcopy (dfftp%npp (me_pool+1) * dfftp%nnp, ddvtosym (npp0, iper), 1, dvtosym &
          (1, iper), 1)
  enddo
  deallocate(ddvtosym)

  call stop_clock ('psymd0rho')
#endif
  return
end subroutine psymd0rho
