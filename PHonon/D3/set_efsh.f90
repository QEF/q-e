!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine set_efsh (drhoscf, imode0, irr, npe)
  !-----------------------------------------------------------------------
  !  This routine calculates the FermiEnergy shift
  !   and stores it in the variable ef_sh
  !
  USE kinds, only : DP
  USE io_global,  ONLY : stdout
  USE fft_base,   ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  use pwcom
  use qpoint,     ONLY :  nksq
  use phcom
  use d3com
  USE mp_global,  ONLY : inter_pool_comm, intra_pool_comm
  USE mp,         ONLY : mp_sum
  implicit none
  integer :: npe, imode0, irr
  ! input: the number of perturbation
  ! input: the position of the current mode
  ! input: index of the current irr. rep.
  complex (DP) :: drhoscf (dfftp%nnr, npe)
  ! input: variation of the charge density

  integer :: ipert, ik, ikk, ibnd
  ! counters
  complex (DP) :: delta_n, def (npertx)
  ! the change in electron number
  ! the change of the Fermi energy for each perturbation
  real (DP) :: weight, wdelta
  ! kpoint weight
  ! delta function weight
  real (DP), save :: dos_ef
  ! density of states at Ef
  real (DP), external :: w0gauss
  logical, save :: first = .true.
  ! Used for initialization
  !
  ! first call: calculates density of states at Ef
  !
  if (first) then
     first = .false.
     dos_ef = 0.d0
     do ik = 1, nksq
        if (lgamma) then
           ikk = ik
        else
           ikk = 2 * ik - 1
        endif
        weight = wk (ikk)
        do ibnd = 1, nbnd
           wdelta = w0gauss ( (ef - et (ibnd, ikk) ) / degauss, ngauss) &
                / degauss
           dos_ef = dos_ef + weight * wdelta
        enddo
     enddo
#ifdef __MPI
     call mp_sum( dos_ef, inter_pool_comm )
#endif
  endif
  !
  ! determines Fermi energy shift (such that each pertubation is neutral)
  !
  WRITE( stdout, * )
  do ipert = 1, npe
     CALL fwfft ('Dense', drhoscf (:, ipert), dfftp)
#ifdef __MPI
     delta_n = (0.d0, 0.d0)
     if (gg (1) < 1.0d-8) delta_n = omega * drhoscf (nl (1), ipert)
     call mp_sum ( delta_n, intra_pool_comm )
#else
     delta_n = omega * drhoscf (nl (1), ipert)
#endif
     def (ipert) = - delta_n / dos_ef
  enddo
  !
  ! symmetrizes the Fermi energy shift
  !
  call sym_def1 (def, irr)
  do ipert = 1, npe
     ef_sh (imode0 + ipert) =  DBLE (def (ipert) )
  enddo

  WRITE( stdout, '(5x,"Pert. #",i3,": Fermi energy shift (Ry) =", &
       &            2f10.4)')  (ipert, def (ipert) , ipert = 1, npe)
  return
end subroutine set_efsh
