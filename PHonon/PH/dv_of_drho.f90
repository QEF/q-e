!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dv_of_drho (mode, dvscf, add_nlcc)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the change of the self consistent potential
  !     due to the perturbation.
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : e2, fpi
  USE fft_base,  ONLY: dfftp
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvect,     ONLY : nl, ngm, g,nlm, gstart
  USE cell_base, ONLY : alat, tpiba2, omega
  USE noncollin_module, ONLY : nspin_lsda, nspin_mag, nspin_gga
  USE funct,     ONLY : dft_is_gradient
  USE scf,       ONLY : rho, rho_core
  USE eqv,       ONLY : dmuxc
  USE nlcc_ph,   ONLY : nlcc_any
  USE qpoint,    ONLY : xq
  USE gc_ph,     ONLY : grho, dvxc_rr,  dvxc_sr,  dvxc_ss, dvxc_s
  USE control_ph, ONLY : lrpa
  USE control_flags, only : gamma_only, tddfpt
  USE martyna_tuckerman, ONLY : wg_corr_h, do_comp_mt
  !OBM: gamma_only is disregarded for phonon calculations, TDDFPT purposes only

  implicit none

  integer :: mode
  ! input: the mode to do

  complex(DP), intent(inout):: dvscf (dfftp%nnr, nspin_mag)
  ! input: the change of the charge,
  ! output: change of the potential

  logical :: add_nlcc
  ! input: if true add core charge

  integer :: ir, is, is1, ig
  ! counter on r vectors
  ! counter on spin polarizations
  ! counter on g vectors

  real(DP) :: qg2, fac
  ! the modulus of (q+G)^2
  ! the structure factor

  complex(DP), allocatable :: dvaux (:,:), drhoc (:)
  !  the change of the core charge
  complex(DP), allocatable :: dvhart (:,:) 
  complex(DP), allocatable :: dvaux_mt(:), rgtot(:)
  ! auxiliary array for Martyna-Tuckerman correction in TDDFPT
  ! total response density  
  real(DP) :: eh_corr
  ! Correction to response Hartree energy due to Martyna-Tuckerman correction 
  ! (only TDDFT). Not used.

  call start_clock ('dv_of_drho')
  allocate (dvaux( dfftp%nnr,  nspin_mag))
  dvaux (:,:) = (0.d0, 0.d0)
  if (add_nlcc) allocate (drhoc( dfftp%nnr))
  !
  ! the exchange-correlation contribution is computed in real space
  !
  if (lrpa) goto 111
  fac = 1.d0 / DBLE (nspin_lsda)
  if (nlcc_any.and.add_nlcc) then
     if (mode > 0) call addcore (mode, drhoc)
     do is = 1, nspin_lsda
        rho%of_r(:, is) = rho%of_r(:, is) + fac * rho_core (:)
        dvscf(:, is) = dvscf(:, is) + fac * drhoc (:)
     enddo
  endif
  do is = 1, nspin_mag
     do is1 = 1, nspin_mag
        do ir = 1, dfftp%nnr
           dvaux(ir,is) = dvaux(ir,is) + dmuxc(ir,is,is1) * dvscf(ir,is1)
        enddo
     enddo
  enddo


  !
  ! add gradient correction to xc, NB: if nlcc is true we need to add here
  ! its contribution. grho contains already the core charge
  !
  if ( dft_is_gradient() ) call dgradcorr &
       (rho%of_r, grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s, xq, &
       dvscf, dfftp%nnr, nspin_mag, nspin_gga, nl, ngm, g, alat, dvaux)
  if (nlcc_any.and.add_nlcc) then
     do is = 1, nspin_lsda
        rho%of_r(:, is) = rho%of_r(:, is) - fac * rho_core (:)
        dvscf(:, is) = dvscf(:, is) - fac * drhoc (:)
     enddo
  endif


111 continue
  !
  ! copy the total (up+down) delta rho in dvscf(*,1) and go to G-space
  !
  if (nspin_mag == 2) then
     dvscf(:,1) = dvscf(:,1) + dvscf(:,2)
  end if
  !
  CALL fwfft ('Dense', dvscf(:,1), dfftp)
  !
  ! Hartree contribution is computed in reciprocal space
  !
  IF (tddfpt .and. do_comp_mt) THEN
      !
      ! TDDFPT plus Martyna-Tuckerman correction
      ! (gamma_only and general k-points are supported)
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
      ! At the end the two contributions (Hartree+XC) are added
      ! 
      dvscf = dvaux + dvhart
      !
      deallocate( dvaux_mt, rgtot ) 
      deallocate(dvhart)
      !
  ELSE
   !
   ! PHonon, and TDDFPT (without Martyna-Tuckerman correction)
   !
   if (gamma_only) then
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
       !  and transformed back to real space
       !
       CALL invfft ('Dense', dvhart (:, is), dfftp)
      enddo
      !
      ! at the end the two contributes are added
      dvscf  = dvaux  + dvhart
      !OBM : Again not totally convinced about this trimming.
      !dvscf (:,:) = cmplx(DBLE(dvscf(:,:)),0.0d0,dp)
      deallocate(dvhart)
   else
    do is = 1, nspin_lsda
       CALL fwfft ('Dense', dvaux (:, is), dfftp)
       do ig = 1, ngm
          qg2 = (g(1,ig)+xq(1))**2 + (g(2,ig)+xq(2))**2 + (g(3,ig)+xq(3))**2
          if (qg2 > 1.d-8) then
             dvaux(nl(ig),is) = dvaux(nl(ig),is) + &
                                e2 * fpi * dvscf(nl(ig),1) / (tpiba2 * qg2)
          endif
       enddo
       !
       !  and transformed back to real space
       !
       CALL invfft ('Dense', dvaux (:, is), dfftp)
    enddo
    !
    ! at the end the two contributes are added
    dvscf (:,:) = dvaux (:,:)
   endif
   !
  ENDIF
  !
  if (add_nlcc) deallocate (drhoc)
  deallocate (dvaux)
  call stop_clock ('dv_of_drho')
  return
end subroutine dv_of_drho
