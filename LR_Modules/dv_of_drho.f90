!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE dv_of_drho_lr

PUBLIC :: dv_of_drho
PUBLIC :: dv_of_drho_xc

CONTAINS

!-----------------------------------------------------------------------
subroutine dv_of_drho (dvscf, drhoc)
  !-----------------------------------------------------------------------
  !
  !  This routine computes the change of the self consistent potential
  !  (Hartree and XC) due to the perturbation.
  !  Note: gamma_only is disregarded for PHonon calculations, 
  !  TDDFPT purposes only.
  !
  USE kinds,             ONLY : DP
  USE constants,         ONLY : e2, fpi
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : fwfft, invfft
  USE gvect,             ONLY : ngm, g, gstart
  USE cell_base,         ONLY : tpiba2, omega
  USE noncollin_module,  ONLY : nspin_lsda, nspin_mag
  USE funct,             ONLY : dft_is_nonlocc
  USE xc_lib,            ONLY : xclib_dft_is
  USE control_flags,     ONLY : gamma_only
  USE martyna_tuckerman, ONLY : wg_corr_h, do_comp_mt
  USE Coul_cut_2D,       ONLY : do_cutoff_2D  
  USE Coul_cut_2D_ph,    ONLY : cutoff_dv_of_drho 
  USE qpoint,            ONLY : xq
  USE control_lr,        ONLY : lrpa, lnolr

  IMPLICIT NONE
  COMPLEX(DP), INTENT(INOUT) :: dvscf(dfftp%nnr, nspin_mag)
  ! input:  response charge density
  ! output: response Hartree-and-XC potential
  COMPLEX(DP), INTENT(IN), OPTIONAL :: drhoc(dfftp%nnr)
  ! input: response core charge density 
  ! (needed only for PHonon when add_nlcc=.true.)
  
  LOGICAL :: add_nlcc
  ! if true add core charge density
  INTEGER :: is, ig
  ! counter on r vectors
  ! counter on spin polarizations
  ! counter on g vectors
  REAL(DP) :: qg2, eh_corr, g2
  ! qg2: the modulus of (q+G)^2
  ! eh_corr: the correction to response Hartree energy due 
  ! to Martyna-Tuckerman correction (calculated, but not used).
  ! g2: modulus of G^2
  COMPLEX(DP), ALLOCATABLE :: dvaux(:,:), dvhart(:,:), & 
                              dvaux_mt(:), rgtot(:)
  ! dvaux: response XC potential 
  ! dvhart: response Hartree potential
  ! dvaux_mt: auxiliary array for Martyna-Tuckerman correction
  ! rgtot: total response density  

  CALL start_clock ('dv_of_drho')
  !
  allocate (dvaux( dfftp%nnr, nspin_mag))
  dvaux (:,:) = (0.d0, 0.d0)
  !
  add_nlcc = PRESENT(drhoc)
  !
  ! 1) The exchange-correlation contribution is computed in real space
  !
  IF (.NOT. lrpa) THEN
     IF (add_nlcc) THEN
        CALL dv_of_drho_xc(dvaux, drho = dvscf, drhoc = drhoc)
     ELSE
        CALL dv_of_drho_xc(dvaux, drho = dvscf)
     ENDIF
  ENDIF
  !
  ! 2) Hartree contribution is computed in reciprocal space
  !
  ! Copy the total (up+down) delta rho in dvscf(*,1) and go to G-space
  !
  if (nspin_mag == 2) then
     dvscf(:,1) = dvscf(:,1) + dvscf(:,2)
  end if
  !
  CALL fwfft ('Rho', dvscf(:,1), dfftp)
  !
  IF (do_comp_mt) THEN
      !
      ! Response Hartree potential with the Martyna-Tuckerman correction
      !
      allocate(dvhart(dfftp%nnr,nspin_mag))
      dvhart(:,:) = (0.d0,0.d0)
      !
      do is = 1, nspin_lsda
        do ig = gstart, ngm
          qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
          dvhart(dfftp%nl(ig),is) = e2 * fpi * dvscf(dfftp%nl(ig),1) / (tpiba2 * qg2)
        enddo
      enddo 
      !
      ! Add Martyna-Tuckerman correction to response Hartree potential
      !
      allocate( dvaux_mt( ngm ), rgtot(ngm) )
      !
      ! Total response density
      !
      do ig = 1, ngm
         rgtot(ig) = dvscf(dfftp%nl(ig),1)
      enddo
      !
      CALL wg_corr_h (omega, ngm, rgtot, dvaux_mt, eh_corr)
      !
      do is = 1, nspin_lsda
        !
        do ig = 1, ngm
           dvhart(dfftp%nl(ig),is)  = dvhart(dfftp%nl(ig),is)  + dvaux_mt(ig)
        enddo
        if (gamma_only) then
           do ig = 1, ngm
              dvhart(dfftp%nlm(ig),is) = conjg(dvhart(dfftp%nl(ig),is))
           enddo
        endif
        !
        ! Transform response Hartree potential to real space
        !
        CALL invfft ('Rho', dvhart (:,is), dfftp)
        !
      enddo
      !
      ! At the end the two contributions (XC+Hartree) are added
      ! 
      dvscf = dvaux + dvhart
      !
      deallocate( dvaux_mt, rgtot ) 
      deallocate(dvhart)
      !
  ELSE
   !
   ! Response Hartree potential (without Martyna-Tuckerman correction)
   !
   If (gamma_only) then
      !
      ! Gamma_only case
      !
      allocate(dvhart(dfftp%nnr,nspin_mag))
      dvhart(:,:) = (0.d0,0.d0)
      !
      do is = 1, nspin_lsda
        do ig = 1, ngm
           qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
           if (qg2 > 1.d-8) then
              dvhart(dfftp%nl(ig),is) = e2 * fpi * dvscf(dfftp%nl(ig),1) / (tpiba2 * qg2)
              dvhart(dfftp%nlm(ig),is)=conjg(dvhart(dfftp%nl(ig),is))
           endif
        enddo
        !
        ! Transformed back to real space
        !
        CALL invfft ('Rho', dvhart (:, is), dfftp)
        !
      enddo
      !
      ! At the end the two contributes are added
      !
      dvscf = dvaux + dvhart
      !
      deallocate(dvhart)
      !
   else
      !
      ! General k points implementation
      !
      do is = 1, nspin_lsda
         CALL fwfft ('Rho', dvaux (:, is), dfftp)
         IF (do_cutoff_2D) THEN 
            call cutoff_dv_of_drho(dvaux, is, dvscf)
         ELSE
            do ig = 1, ngm
               qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
               g2  = g(1,ig)**2 + g(2,ig)**2 + g(3,ig)**2
               IF (lnolr) THEN
                 !
                 IF (g2 > 1d-8 .AND. qg2 > 1.d-8) THEN
                    dvaux(dfftp%nl(ig),is) = dvaux(dfftp%nl(ig),is) + &
                                   & e2 * fpi * dvscf(dfftp%nl(ig),1) / (tpiba2 * qg2)
                 ENDIF
                 ! 
               ELSE
                 if (qg2 > 1.d-8) then
                    dvaux(dfftp%nl(ig),is) = dvaux(dfftp%nl(ig),is) + &
                                   & e2 * fpi * dvscf(dfftp%nl(ig),1) / (tpiba2 * qg2)
                 endif
               ENDIF
            enddo
         ENDIF
         !
         ! Transformed back to real space
         !
         CALL invfft ('Rho', dvaux (:, is), dfftp)
         !
      enddo
      !
      ! At the end the two contributes are added
      !
      dvscf (:,:) = dvaux (:,:)
      !
   endif
   !
  ENDIF
  !
  deallocate (dvaux)
  !
  CALL stop_clock ('dv_of_drho')
  !
  RETURN
  !
end subroutine dv_of_drho

!-----------------------------------------------------------------------
subroutine dv_of_drho_xc (dv, drho, drhoc)
  !-----------------------------------------------------------------------
  !
  !  This routine computes the change of the self consistent potential
  !  (Hartree and XC) due to the perturbation.
  !  Note: gamma_only is disregarded for PHonon calculations,
  !  TDDFPT purposes only.
  !
  USE kinds,             ONLY : DP
  USE constants,         ONLY : e2, fpi
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : fwfft, invfft
  USE gvect,             ONLY : g
  USE noncollin_module,  ONLY : nspin_lsda, nspin_mag, nspin_gga
  USE funct,             ONLY : dft_is_nonlocc
  USE xc_lib,            ONLY : xclib_dft_is
  USE scf,               ONLY : rho, rho_core
  USE uspp,              ONLY : nlcc_any
  USE Coul_cut_2D_ph,    ONLY : cutoff_dv_of_drho
  USE qpoint,            ONLY : xq
  USE gc_lr,             ONLY : grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
  USE eqv,               ONLY : dmuxc

  IMPLICIT NONE
  COMPLEX(DP), INTENT(INOUT) :: dv(dfftp%nnr, nspin_mag)
  ! output: response XC potential
  COMPLEX(DP), INTENT(IN), OPTIONAL :: drho(dfftp%nnr, nspin_mag)
  ! input:  response charge density
  COMPLEX(DP), INTENT(IN), OPTIONAL :: drhoc(dfftp%nnr)
  ! input: response core charge density
  ! (needed only for PHonon when add_nlcc=.true.)

  INTEGER :: ir, is, is1
  ! counter on r vectors
  ! counter on spin polarizations
  REAL(DP) :: fac
  ! fac: 1 / nspin_lsda
  COMPLEX(DP), ALLOCATABLE :: drhotot(:, :)
  ! Total charge density. Sum of electronic (drho) and nlcc (drhoc) terms.
  !
  ALLOCATE(drhotot(dfftp%nnr, nspin_mag))
  !
  IF (PRESENT(drho)) THEN
     drhotot = drho
  ELSE
     drhotot = (0.d0, 0.d0)
  ENDIF
  !
  IF (PRESENT(drhoc)) THEN
     ! Add core charge
     fac = 1.d0 / DBLE (nspin_lsda)
     DO is = 1, nspin_lsda
        drhotot(:, is) = drhotot(:, is) + fac * drhoc(:)
     ENDDO
  ENDIF
  !
  do is = 1, nspin_mag
     do is1 = 1, nspin_mag
        do ir = 1, dfftp%nnr
           dv(ir,is) = dv(ir,is) + dmuxc(ir,is,is1) * drhotot(ir,is1)
        enddo
     enddo
  enddo
  !
  ! Add gradient correction to the response XC potential.
  ! NB: If nlcc=.true. we need to add here its contribution.
  ! grho contains already the core charge
  !
  if (nlcc_any) rho%of_r(:, 1) = rho%of_r(:, 1) + rho_core (:)
  !
  if ( xclib_dft_is('gradient') ) call dgradcorr(dfftp, rho%of_r, grho, dvxc_rr, &
                                dvxc_sr, dvxc_ss, dvxc_s, xq, drhotot, &
                                nspin_mag, nspin_gga, g, dv)
  !
  if ( dft_is_nonlocc() )  call dnonloccorr(rho%of_r, drhotot, xq, dv)
  !
  if (nlcc_any) rho%of_r(:, 1) = rho%of_r(:, 1) - rho_core (:)
  !
end subroutine dv_of_drho_xc

END MODULE dv_of_drho_lr
