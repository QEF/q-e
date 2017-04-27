!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------
SUBROUTINE addnlcc_zstar_eu_us( drhoscf )
!----------===================-------------------

  USE kinds, ONLY : DP
  USE funct, only : dft_is_gradient, dft_is_nonlocc
  USE scf, only : rho, rho_core
  USE cell_base, ONLY : omega, alat
  USE gvect, ONLY : ngm, nl, g
  USE fft_base, ONLY : dfftp
  USE noncollin_module, ONLY : nspin_lsda, nspin_gga, nspin_mag
  USE efield_mod, ONLY : zstareu0
  USE uspp,   ONLY : nlcc_any
  USE modes,  ONLY : npert, nirr

  USE mp_pools, ONLY : my_pool_id

  USE qpoint, ONLY : xq
  USE eqv,     ONLY : dmuxc
  USE gc_lr,   ONLY: grho, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s

  IMPLICIT NONE

  COMPLEX(DP) :: drhoscf (dfftp%nnr,nspin_mag,3)


  INTEGER :: nrtot, ipert, jpert, is, is1, irr, ir, mode, mode1
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
           CALL addcore (mode, drhoc)

           DO is = 1, nspin_lsda
              rho%of_r(:,is) = rho%of_r(:,is) + fac * rho_core
           END DO

           DO is = 1, nspin_mag
              DO is1 = 1, nspin_mag
                 DO ir = 1, dfftp%nnr
                    dvaux (ir, is) = dvaux (ir, is) +     &
                         dmuxc (ir, is, is1) *            &
                         drhoscf (ir, is1, ipol)
                 ENDDO
              ENDDO
           END DO
           !
           ! add gradient correction to xc, NB: if nlcc is true we need to add here
           ! its contribution. grho contains already the core charge
           !

           IF ( dft_is_gradient() ) &
                CALL dgradcorr (rho%of_r, grho, &
                    dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, drhoscf (1,1,ipol),&
                    dfftp%nnr, nspin_mag, nspin_gga, nl, ngm, g, alat, dvaux)
           if (dft_is_nonlocc()) &
                call dnonloccorr(rho%of_r, drhoscf (1, 1, ipol), xq, dvaux)

           DO is = 1, nspin_lsda
              rho%of_r(:,is) = rho%of_r(:,is) - fac * rho_core
           END DO

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
