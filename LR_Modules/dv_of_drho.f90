!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE dv_of_drho_lr

CONTAINS

!-----------------------------------------------------------------------
subroutine dv_of_drho (dvscf, add_nlcc, drhoc)
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
  USE gvect,             ONLY : nl, ngm, g,nlm, gstart
  USE cell_base,         ONLY : alat, tpiba2, omega
  USE noncollin_module,  ONLY : nspin_lsda, nspin_mag, nspin_gga
  USE funct,             ONLY : dft_is_gradient, dft_is_nonlocc
  USE scf,               ONLY : rho, rho_core
  USE uspp,              ONLY : nlcc_any
  USE control_flags,     ONLY : gamma_only
  USE martyna_tuckerman, ONLY : wg_corr_h, do_comp_mt
  USE qpoint,            ONLY : xq
  USE gc_lr,             ONLY : grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
  USE control_lr,        ONLY : lrpa
  USE eqv,               ONLY : dmuxc

  IMPLICIT NONE
  COMPLEX(DP), INTENT(INOUT) :: dvscf(dfftp%nnr, nspin_mag)
  ! input:  response charge density
  ! output: response Hartree-and-XC potential
  LOGICAL, INTENT(IN) :: add_nlcc
  ! input: if true add core charge density
  COMPLEX(DP), INTENT(IN), OPTIONAL :: drhoc(dfftp%nnr)
  ! input: response core charge density 
  ! (needed only for PHonon when add_nlcc=.true.)
  
  INTEGER :: ir, is, is1, ig
  ! counter on r vectors
  ! counter on spin polarizations
  ! counter on g vectors
  REAL(DP) :: qg2, fac, eh_corr
  ! qg2: the modulus of (q+G)^2
  ! fac: the structure factor
  ! eh_corr: the correction to response Hartree energy due 
  ! to Martyna-Tuckerman correction (calculated, but not used).
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
  if (add_nlcc .and. .not.present(drhoc)) &
     & CALL errore( 'dv_of_drho', 'drhoc is not present in the input of the routine', 1 )   
  !
  ! 1) The exchange-correlation contribution is computed in real space
  !
  if (lrpa) goto 111
  !
  fac = 1.d0 / DBLE (nspin_lsda)
  !
  if (nlcc_any.and.add_nlcc) then
     do is = 1, nspin_lsda
        rho%of_r(:, is) = rho%of_r(:, is) + fac * rho_core (:)
        dvscf(:, is) = dvscf(:, is) + fac * drhoc (:)
     enddo
  endif
  !
  do is = 1, nspin_mag
     do is1 = 1, nspin_mag
        do ir = 1, dfftp%nnr
           dvaux(ir,is) = dvaux(ir,is) + dmuxc(ir,is,is1) * dvscf(ir,is1)
        enddo
     enddo
  enddo
  !
  ! Add gradient correction to the response XC potential.
  ! NB: If nlcc=.true. we need to add here its contribution. 
  ! grho contains already the core charge
  !
  if ( dft_is_gradient() ) call dgradcorr &
       (rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
       dvscf, dfftp%nnr, nspin_mag, nspin_gga, nl, ngm, g, alat, dvaux)
  !
  if (dft_is_nonlocc()) then
     call dnonloccorr(rho%of_r, dvscf, xq, dvaux)
  endif

  if (nlcc_any.and.add_nlcc) then
     do is = 1, nspin_lsda
        rho%of_r(:, is) = rho%of_r(:, is) - fac * rho_core (:)
        dvscf(:, is) = dvscf(:, is) - fac * drhoc (:)
     enddo
  endif
  !
111 continue
  !
  ! Copy the total (up+down) delta rho in dvscf(*,1) and go to G-space
  !
  if (nspin_mag == 2) then
     dvscf(:,1) = dvscf(:,1) + dvscf(:,2)
  end if
  !
  CALL fwfft ('Dense', dvscf(:,1), dfftp)
  !
  ! 2) Hartree contribution is computed in reciprocal space
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
          dvhart(nl(ig),is) = e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
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
         rgtot(ig) = dvscf(nl(ig),1)
      enddo
      !
      CALL wg_corr_h (omega, ngm, rgtot, dvaux_mt, eh_corr)
      !
      do is = 1, nspin_lsda
        !
        do ig = 1, ngm
           dvhart(nl(ig),is)  = dvhart(nl(ig),is)  + dvaux_mt(ig)
        enddo
        if (gamma_only) then
           do ig = 1, ngm
              dvhart(nlm(ig),is) = conjg(dvhart(nl(ig),is))
           enddo
        endif
        !
        ! Transform response Hartree potential to real space
        !
        CALL invfft ('Dense', dvhart (:,is), dfftp)
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
              dvhart(nl(ig),is) = e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
              dvhart(nlm(ig),is)=conjg(dvhart(nl(ig),is))
           endif
        enddo
        !
        ! Transformed back to real space
        !
        CALL invfft ('Dense', dvhart (:, is), dfftp)
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
         CALL fwfft ('Dense', dvaux (:, is), dfftp)
         do ig = 1, ngm
            qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
            if (qg2 > 1.d-8) then
               dvaux(nl(ig),is) = dvaux(nl(ig),is) + &
                              & e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
            endif
         enddo
         !
         ! Transformed back to real space
         !
         CALL invfft ('Dense', dvaux (:, is), dfftp)
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

END MODULE dv_of_drho_lr
