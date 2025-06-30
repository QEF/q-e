!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------
SUBROUTINE addnlcc_zstar_eu_us( drhop )
!----------===================-------------------

  USE kinds, ONLY : DP
  USE mp_pools, ONLY : my_pool_id
  USE cell_base, ONLY : omega
  USE fft_base, ONLY : dfftp
  USE noncollin_module, ONLY : nspin_lsda,  nspin_mag
  USE efield_mod, ONLY : zstareu0
  USE uspp,   ONLY : nlcc_any
  USE modes,  ONLY : npert, nirr, u
  USE dv_of_drho_lr,    ONLY : dv_of_drho_xc

  IMPLICIT NONE

  COMPLEX(DP) :: drhop (dfftp%nnr,nspin_mag,3)


  INTEGER :: nrtot, ipert, is, irr, mode
  INTEGER :: imode0, npe, ipol

  REAL(DP) :: fac

  COMPLEX(DP), DIMENSION(dfftp%nnr) :: drhoc
  COMPLEX(DP), DIMENSION(dfftp%nnr,nspin_mag) :: dvaux

  IF (.NOT.nlcc_any) RETURN

  IF ( my_pool_id /= 0 ) RETURN

  DO ipol = 1, 3
     imode0 = 0
     DO irr = 1, nirr
        npe = npert(irr)
        !
        !  compute the exchange and correlation potential for this mode
        !
        nrtot = dfftp%nr1 * dfftp%nr2 * dfftp%nr3
        fac = 1.d0 / DBLE (nspin_lsda)
        DO ipert = 1, npe
           mode = imode0 + ipert

           dvaux = (0.0_dp,0.0_dp)
           CALL addcore(u(1, mode), drhoc)
           !
           CALL dv_of_drho_xc(dvaux, drho = drhop(1, 1, ipol))
           !
           DO is = 1, nspin_lsda
              zstareu0(ipol,mode) = zstareu0(ipol,mode) -                  &
                   omega * fac / REAL(nrtot, DP) *         &
                   DOT_PRODUCT(dvaux(1:dfftp%nnr,is),drhoc(1:dfftp%nnr))
           END DO
        END DO
        imode0 = imode0 + npe
     END DO
  END DO

  RETURN
END SUBROUTINE addnlcc_zstar_eu_us
