!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine drhod2v
  !-----------------------------------------------------------------------
  ! It calls the routines which calculate the term containing the first
  ! variation of the charge and the secon variation of the potential with
  ! respect to the perturbation.
  ! d0rhod2v: contains the terms depending on the first variation of the c
  ! with respect to a perturbaation at q=0
  ! dqrhod2v: contains the terms depending on the first variation of the c
  ! with respect to a perturbaation at a generic q
  ! The variation of the charge can be read from a file or calculated dire
  ! --this last option is to be used for testing pourposes--
  !
  USE ions_base,  ONLY : nat
  USE kinds, only : DP
  USE fft_base, ONLY : dfftp
  use pwcom
  use phcom
  use d3com
  !
  implicit none
  integer :: irr, irr1, imode0, ipert, ir
  real (DP) :: xq0 (3)
  complex (DP), allocatable :: drhoscf (:)
  ! the change of density due to perturbations

  allocate  (drhoscf( dfftp%nnr))

  call read_ef
  if (.not.allmodes) then
     do ipert = 1, 3 * nat
        call davcio_drho (drhoscf, lrdrho, iudrho, ipert, - 1)
        call dqrhod2v (ipert, drhoscf)
     enddo

  endif
  do ipert = 1, 3 * nat
     if (q0mode (ipert) ) then
        call davcio_drho (drhoscf, lrdrho, iud0rho, ipert, - 1)
        call d0rhod2v (ipert, drhoscf)
     endif

  enddo

  deallocate (drhoscf)
  return

end subroutine drhod2v
