!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine compute_dvloc (mode, dvlocin)
  !----------------------------------------------------------------------
  !! This routine calculates \(dV_\text{bare}/\text{dtau}\cdot\text{psi}\)
  !! for one perturbation with a given q. The displacements are described 
  !! by a vector u.  
  !! The result is stored in \(\text{dvpsi}\). The routine is called for
  !! each k point and for each pattern u. It computes simultaneously all
  !! the bands.
  !
  !
  USE kinds,     ONLY : DP
  USE fft_base,  ONLY : dffts
  USE fft_interfaces, ONLY: invfft
  USE gvect,     ONLY : eigts1, eigts2, eigts3, mill, g
  USE gvecs,   ONLY : ngms
  USE cell_base, ONLY : tpiba
  USE ions_base, ONLY : nat, ityp
  USE modes,     ONLY : u
  USE qpoint,    ONLY : xq, eigqts
  USE eqv,       ONLY : evq, vlocq
  USE Coul_cut_2D,    ONLY: do_cutoff_2D 
  USE Coul_cut_2D_ph, ONLY: cutoff_localq 
  implicit none
  !

  integer :: mode
  !! input: the actual perturbation
  complex(DP) :: dvlocin(dffts%nnr)
  !! output: the change of the local potential
  !
  ! ... local variables
  !
  integer :: na, nt, mu, ig
  ! counters
  complex(DP) :: gtau, gu, fact, u1, u2, u3, gu0
  ! auxiliary variables

  call start_clock ('com_dvloc')
  dvlocin (:) = (0.d0, 0.d0)
  do na = 1, nat
     fact = tpiba * (0.d0, -1.d0) * eigqts (na)
     mu = 3 * (na - 1)
     if ( abs (u (mu + 1, mode) ) + abs (u (mu + 2, mode) ) + &
          abs (u (mu + 3, mode) ) > 1.0d-12) then
        nt = ityp (na)
        u1 = u (mu + 1, mode)
        u2 = u (mu + 2, mode)
        u3 = u (mu + 3, mode)
        gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
        do ig = 1, ngms
           gtau = eigts1 (mill(1,ig), na) * eigts2 (mill(2,ig), na) * &
                  eigts3 (mill(3,ig), na)
           gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
           dvlocin (dffts%nl (ig) ) = dvlocin (dffts%nl (ig) ) + vlocq (ig, nt) &
                * gu * fact * gtau
        enddo
        IF (do_cutoff_2D) then
           call cutoff_localq( dvlocin, fact, u1, u2, u3, gu0, nt, na) ! 2D: re-add cutoff LR part
        ENDIF
     endif
  enddo
  !
  ! Now we compute dV_loc/dtau in real space
  !

  CALL invfft ('Rho', dvlocin, dffts)

  call stop_clock ('com_dvloc')
  return
end subroutine compute_dvloc
