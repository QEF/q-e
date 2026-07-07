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
  INTEGER :: tpnnr

  CALL start_clock ('dv_of_drho')
  !
  tpnnr = dfftp%nnr
  !
  allocate (dvaux( dfftp%nnr, nspin_mag))
  !
!$acc data present_or_copy(dvscf(1:tpnnr,1:nspin_mag)) create(dvaux(1:tpnnr,1:nspin_mag)) 
  !
!$acc kernels  
  dvaux (:,:) = (0.d0, 0.d0)
!$acc end kernels
  !
  add_nlcc = PRESENT(drhoc)
  !
  ! 1) The exchange-correlation contribution is computed in real space
  !
  IF (.NOT. lrpa) THEN
     IF (add_nlcc) THEN
!$acc data present_or_copyin(drhoc(1:tpnnr))
        CALL dv_of_drho_xc(dvaux, drho = dvscf, drhoc = drhoc)
!$acc end data        
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
!$acc kernels          
     dvscf(:,1) = dvscf(:,1) + dvscf(:,2)
!$acc end kernels
  end if
  !
  !$acc host_data use_device(dvscf)
  CALL fwfft ('Rho', dvscf(:,1), dfftp)
  !$acc end host_data
  !
  IF (do_comp_mt) THEN
     !
     ! Response Hartree potential with the Martyna-Tuckerman correction
     !
     allocate(dvhart(dfftp%nnr,nspin_mag))
     allocate( dvaux_mt( ngm ), rgtot(ngm) )
     !
!$acc enter data create(dvhart, dvaux_mt, rgtot)     
     !
!$acc kernels     
     dvhart(:,:) = (0.d0,0.d0)
!$acc end kernels
     !
!$acc parallel loop collapse(2)     
     do is = 1, nspin_lsda
        do ig = gstart, ngm
           qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
           dvhart(dfftp%nl(ig),is) = e2 * fpi * dvscf(dfftp%nl(ig),1) / (tpiba2 * qg2)
        enddo
     enddo 
     !
     ! Add Martyna-Tuckerman correction to response Hartree potential
     !
     ! Total response density
     !
!$acc parallel loop     
     do ig = 1, ngm
        rgtot(ig) = dvscf(dfftp%nl(ig),1)
     enddo
     !
     CALL wg_corr_h (omega, ngm, rgtot, dvaux_mt, eh_corr)
     !
     do is = 1, nspin_lsda
        !
!$acc parallel loop       
        do ig = 1, ngm
           dvhart(dfftp%nl(ig),is)  = dvhart(dfftp%nl(ig),is)  + dvaux_mt(ig)
        enddo
        if (gamma_only) then
!$acc parallel loop               
           do ig = 1, ngm
              dvhart(dfftp%nlm(ig),is) = conjg(dvhart(dfftp%nl(ig),is))
           enddo
        endif
        !
        ! Transform response Hartree potential to real space
        !
        !$acc host_data use_device(dvhart)
        CALL invfft ('Rho', dvhart (:,is), dfftp)
        !$acc end host_data
        !
     enddo
     !
     ! At the end the two contributions (XC+Hartree) are added
     !
!$acc kernels 
     dvscf(:,:) = dvaux(:,:) + dvhart(:,:)
!$acc end kernels     
     !
!$acc exit data delete(dvhart,dvaux_mt,rgtot)
     !
     deallocate( dvaux_mt, rgtot ) 
     deallocate( dvhart )
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
!$acc enter data create(dvhart)
        !
!$acc kernels        
        dvhart(:,:) = (0.d0,0.d0)
!$acc end kernels        
        !
        do is = 1, nspin_lsda
           !$acc parallel loop
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
           !$acc host_data use_device(dvhart)
           CALL invfft ('Rho', dvhart (:, is), dfftp)
           !$acc end host_data
           !
        enddo
        !
        ! At the end the two contributes are added
        !
!$acc kernels        
        dvscf(:,:) = dvaux(:,:) + dvhart(:,:)
!$acc end kernels        
        !
!$acc exit data delete(dvhart)        
        deallocate(dvhart)
        !
     else
        !
        ! General k points implementation
        !
        do is = 1, nspin_lsda
           !$acc host_data use_device(dvaux)
           CALL fwfft ('Rho', dvaux (:, is), dfftp)
           !$acc end host_data
           IF (do_cutoff_2D) THEN 
              call cutoff_dv_of_drho(dvaux, is, dvscf)
           ELSE
!$acc parallel loop                  
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
           !$acc host_data use_device(dvaux)
           CALL invfft ('Rho', dvaux (:, is), dfftp)
           !$acc end Host_data
           !
        enddo
        !
        ! At the end the two contributes are added
        !
!$acc kernels        
        dvscf (:,:) = dvaux (:,:)
!$acc end kernels        
        !
     endif
     !
  ENDIF
  !
!$acc end data  
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
  INTEGER :: tpnnr
  !
  tpnnr = dfftp%nnr
  ALLOCATE(drhotot(dfftp%nnr, nspin_mag))
  !
!!!$acc data present_or_copy(dv(1:tpnnr, 1:nspin_mag)) create(drhotot(1:tpnnr, 1:nspin_mag))
!!$acc data copy(dv(1:tpnnr, 1:nspin_mag)) copyout(drhotot(1:tpnnr, 1:nspin_mag)) copyin(dmuxc(1:tpnnr,1:nspin_mag, 1:nspin_mag))
!$acc data present_or_copy(dv) create(drhotot) present_or_copyin(dmuxc)
  !
  IF (PRESENT(drho)) THEN
!$acc data present_or_copyin(drho(1:tpnnr, 1:nspin_mag))          
!$acc kernels          
     drhotot(:,:) = drho(:,:)
!$acc end kernels     
!$acc end data
  ELSE
!$acc kernels          
     drhotot(:,:) = (0.d0, 0.d0)
!$acc end kernels     
  ENDIF
  !
  IF (PRESENT(drhoc)) THEN
!$acc data present_or_copyin(drhoc(1:tpnnr))          
     ! Add core charge
     fac = 1.d0 / DBLE (nspin_lsda)
!$acc parallel loop collapse(2)     
     DO is = 1, nspin_lsda
        do ir = 1, tpnnr
           drhotot(ir, is) = drhotot(ir, is) + fac * drhoc(ir)
        enddo   
     ENDDO
!$acc end data     
  ENDIF
  !
  do is1 = 1, nspin_mag
!$acc parallel loop collapse(2)
     do is = 1, nspin_mag
        do ir = 1, tpnnr
           dv(ir,is) = dv(ir,is) + dmuxc(ir,is,is1) * drhotot(ir,is1)
        enddo
     enddo
  enddo
!!!$acc update self(dv)
!!!$acc end data
  
  
  !
  ! Add gradient correction to the response XC potential.
  ! NB: If nlcc=.true. we need to add here its contribution.
  ! grho contains already the core charge
  !

  if ( xclib_dft_is('gradient') ) THEN
     !$acc update self(dv,drhotot)     
     if (nlcc_any) rho%of_r(:, 1) = rho%of_r(:, 1) + rho_core (:)
     !
     call dgradcorr(dfftp, rho%of_r, grho, dvxc_rr, &
                    dvxc_sr, dvxc_ss, dvxc_s, xq, drhotot, &
                    nspin_mag, nspin_gga, g, dv) 
     !
     if ( dft_is_nonlocc() )  call dnonloccorr(rho%of_r, drhotot, xq, dv)
     !     
     if (nlcc_any) rho%of_r(:, 1) = rho%of_r(:, 1) - rho_core (:)
     !
     !$acc update device(dv)
  elseif ( dft_is_nonlocc() ) THEN
     !$acc update self(dv,drhotot)
     if (nlcc_any) rho%of_r(:, 1) = rho%of_r(:, 1) + rho_core (:)
     !
     call dnonloccorr(rho%of_r, drhotot, xq, dv)
     !
     if (nlcc_any) rho%of_r(:, 1) = rho%of_r(:, 1) - rho_core (:)
     !
     !$acc update device(dv)
  endif        
  !
!$acc end data
  !
!  if (nlcc_any) rho%of_r(:, 1) = rho%of_r(:, 1) + rho_core (:)
!  !
!  if ( xclib_dft_is('gradient') ) call dgradcorr(dfftp, rho%of_r, grho, dvxc_rr, &
!                                dvxc_sr, dvxc_ss, dvxc_s, xq, drhotot, &
!                                nspin_mag, nspin_gga, g, dv)
!  !
!  if ( dft_is_nonlocc() )  call dnonloccorr(rho%of_r, drhotot, xq, dv)
!  !
!  if (nlcc_any) rho%of_r(:, 1) = rho%of_r(:, 1) - rho_core (:)
  !
  DEALLOCATE(drhotot)
  !
end subroutine dv_of_drho_xc

END MODULE dv_of_drho_lr
