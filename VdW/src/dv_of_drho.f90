!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dv_of_drho_vdw (mode, dvscf, flag)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the self consistent potential
  !     due to the perturbation.
  !
  USE funct, ONLY : dft_is_gradient, dmxc
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE pwcom
  USE scf, ONLY : rho, rho_core
  USE kinds, ONLY : DP
  USE phcom
  USE eff_v,  ONLY : rho_veff
  USE uspp,   ONLY : nlcc_any
  IMPLICIT NONE

  INTEGER :: mode
  ! input: the mode to do

  COMPLEX(kind=DP) :: dvscf (dfftp%nnr, nspin)
  ! input: the change of the charge,
  ! output: change of the potential

  LOGICAL :: flag
  ! input: if true add core charge

  INTEGER :: ir, is, is1, ig
  ! counter on r vectors
  ! counter on spin polarizations
  ! counter on g vectors

  real(kind=DP) :: qg2, fac, ttd
  ! the modulus of (q+G)^2
  ! the structure factor
  ! constant 2/3

  COMPLEX(kind=DP), ALLOCATABLE :: dvaux (:,:), drhoc (:),&
                                   dv_tfvw (:,:)
  ! auxiliary variable for potential
  !  the change of the core charge
  ! auxiliary variable for derivation of charge density
  !
  real(kind=dp) :: rhotot
  !
  CALL start_clock ('dv_of_drho')
  ALLOCATE (dvaux( dfftp%nnr,  nspin), dv_tfvw( dfftp%nnr, nspin))
  ALLOCATE (drhoc( dfftp%nnr))
  !
  dv_tfvw = dvscf
  !
  ! the exchange-correlation contribution is computed in real space
  !
  dvaux (:,:) = (0.d0, 0.d0)

  fac = 1.d0 / dble (nspin)
  IF (nlcc_any.and.flag) THEN
     CALL addcore (mode, drhoc)
     DO is = 1, nspin
        rho%of_r(:, is) = rho%of_r(:, is) + fac * rho_core (:)
        dvscf(:, is) = dvscf(:, is) + fac * drhoc (:)
     ENDDO
  ENDIF
!  allocate ( dmuxc(dfftp%nnr, nspin, nspin) )
  dmuxc(:,:,:) = 0.d0
  DO ir = 1, dfftp%nnr
     rhotot = rho%of_r (ir, nspin) + rho_core (ir)
     IF (rhotot>1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)
     IF (rhotot< - 1.d-30) dmuxc (ir, 1, 1) = - dmxc ( - rhotot)
  ENDDO
  DO is = 1, nspin
     DO is1 = 1, nspin
        DO ir = 1, dfftp%nnr
           dvaux(ir,is) = dvaux(ir,is) + dmuxc(ir,is,is1) * dvscf(ir,is1)
        ENDDO
     ENDDO
  ENDDO
!  deallocate ( dmuxc )
  !
  ! add gradient correction to xc, NB: if nlcc is true we need to add here
  ! its contribution. grho contains already the core charge
  !
!  if (igcx /= 0 .or. igcc /= 0) call dgradcorr &
!       (rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
!       dvscf, nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx, nspin, nl, ngm, g, &
!       alat, omega, dvaux)
  IF (nlcc_any.and.flag) THEN
     DO is = 1, nspin
        rho%of_r(:, is) = rho%of_r(:, is) - fac * rho_core (:)
        dvscf(:, is) = dvscf(:, is) - fac * drhoc (:)
     ENDDO
  ENDIF
  !
  ! copy the total (up+down) delta rho in dvscf(*,1) and go to G-space
  !
  IF (nspin == 2) THEN
     dvscf(:,1) = dvscf(:,1) + dvscf(:,2)
  ENDIF
  !
  CALL fwfft ('Dense', dvscf(:,1), dfftp)
  !
  ! hartree contribution is computed in reciprocal space
  !
  DO is = 1, nspin
     CALL fwfft ('Dense', dvaux (:, is), dfftp)
     DO ig = 1, ngm
        qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
        IF (qg2 > 1.d-8) THEN
           dvaux(nl(ig),is) = dvaux(nl(ig),is) + &
                              e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
!           dvaux(nl(ig),is) = e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
        ENDIF
     ENDDO
     !
     !  and transformed back to real space
     !
     CALL invfft ('Dense', dvaux (:, is), dfftp)
  ENDDO
  !
  ! TFvW contribution is computed in real space
  !
  ttd = 2.d0/3.d0
  is = 1
  dv_tfvw(:, is) = ttd * (0.125d0/ttd*fpi**2)**ttd * dv_tfvw (:, is) &
                       * ( abs(rho_veff(:,is))** (-ttd/2.d0) )
!  dv_tfvw(:, is) = ttd * (0.125d0/ttd*fpi**2)**ttd * dv_tfvw (:, is) &
!                       * ( abs(rhof_r(:,is)+rho_core(:))** (-ttd/2.d0) )
  !
  ! at the end the three contributes are added
  !
  dvscf (:,:) = dvaux (:,:) + dv_tfvw (:,:)
  !
  DEALLOCATE (drhoc)
  DEALLOCATE (dvaux)

  CALL stop_clock ('dv_of_drho')
  RETURN
END SUBROUTINE dv_of_drho_vdw
