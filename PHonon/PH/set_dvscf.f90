!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
      subroutine set_dvscf (dvscfs)
  !-----------------------------------------------------------------------
  !! Read the variation of the charge and calculates the local part of
  !! the scf potential.
  !

  use kinds, only : DP
  USE gvecs, ONLY : doublegrid
  USE fft_base, ONLY : dfftp, dffts
  USE fft_interfaces, ONLY : fft_interpolate
  USE lsda_mod,ONLY : nspin
  USE units_ph, ONLY : iudrho, lrdrho
  USE output,   ONLY : fildrho
  USE dv_of_drho_lr
  implicit none

  complex(DP) :: dvscfs (dffts%nnr,3)
  complex(DP) , allocatable :: derho (:,:)
  integer :: ipl
  !  counter on the polarizations

  allocate (derho ( dfftp%nnr, nspin))

  if ( fildrho.eq.' ') call errore ('set_dvscf','where is fildrho?',1)
  !
  do ipl = 1, 3
     !
     ! read from file the variation of the charge
     !
     call davcio_drho (derho (1, 1), lrdrho, iudrho, ipl, -1)
     !
     ! Calculates the local part of the scf potential
     !
     call dv_of_drho (derho (1, 1), .false.)
     !
     if (doublegrid) then
        call fft_interpolate(dfftp, derho (:, 1), dffts, dvscfs (:, ipl) )
     else
        call zcopy (dfftp%nnr, derho (1, 1), 1, dvscfs (1, ipl), 1)
     endif
  end do

  deallocate (derho)

  return
end subroutine set_dvscf
