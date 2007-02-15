! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE compute_sigma_bare(chi_bare, sigma_bare)
  !-----------------------------------------------------------------------
  !
  ! ... Compute the bare contribution to the chemical shift at the
  ! ... position of the nuclei, given the induced field
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk
  USE wvfct,                ONLY : nbnd, npwx, npw, igk  
  USE gvect,                ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, &
                                   nrx3, nrxx, nl, nlm, g, gg, ecutwfc, gcutm
  USE ions_base,            ONLY : nat, tau, atm, ityp
  USE io_global,            ONLY : stdout
  USE symme,                ONLY : s, nsym, irt
  USE pwcom
  USE gipaw_module,         ONLY : use_nmr_macroscopic_shape, &
                                   nmr_macroscopic_shape, b_ind, tens_fmt

  ! Arguments
  REAL(DP), INTENT(IN) :: chi_bare(3,3)
  real(dp), intent(out) :: sigma_bare(3,3,nat)
  
  ! Local
  integer :: na, ig
  real(dp) :: arg, tr_sigma
  complex(dp) :: tmp_sigma(3,3)
  
  write(stdout,'(5X,''NMR chemical bare shifts in ppm:'')')
  write(stdout,*)
  
  do na = 1, nat
    tmp_sigma(:,:) = 0.0_dp
    
    do ig = gstart, ngm              
      arg = (g(1,ig)*tau(1,na) + g(2,ig)*tau(2,na) + g(3,ig)*tau(3,na)) * tpi
      tmp_sigma(:,:) = tmp_sigma(:,:) &
           + b_ind(ig,:,:) * cmplx(cos(arg),sin(arg))
    enddo
    
    if ( use_nmr_macroscopic_shape ) then
       ! this is the G = 0 term
       if (gstart == 2) then
          tmp_sigma(:,:) = tmp_sigma(:,:) &
               - (4.0_dp*pi) * nmr_macroscopic_shape(:,:) * chi_bare(:,:)
       end if
    end if
    
    sigma_bare(:,:,na) = real(tmp_sigma(:,:))
  enddo
#ifdef __PARA
  call reduce(9*na, sigma_bare)
#endif
  
#if 0
  ! symmetrize tensors ??
  do na = 1, nat
    call trntns (sigma_bare(1,1,na), at, bg, -1)
  enddo
  call symz(sigma_bare, nsym, s, nat, irt)
  do na = 1, nat
    call trntns (sigma_bare(1,1,na), at, bg, 1)
  enddo
#endif
  
  do na = 1, nat
    tr_sigma = (sigma_bare(1,1,na)+sigma_bare(2,2,na)+sigma_bare(3,3,na))/3.0_dp
    write(stdout, &
     '(5X,''Atom'',I3,2X,A3,'' pos: ('',3(F10.6),'')  sigma: '', F14.4)') &
     na, atm(ityp(na)), tau(:,na), tr_sigma*1d6
    write(stdout, tens_fmt) sigma_bare(:,:,na) * 1d6
  enddo
  
end subroutine compute_sigma_bare

!-----------------------------------------------------------------------
SUBROUTINE compute_sigma_diamagnetic( sigma_diamagnetic )
  !-----------------------------------------------------------------------
  !
  ! ... Compute the diamagnetic contribution to the chemical shift at the
  ! ... position of the nuclei, given the induced field
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk
  USE wvfct,                ONLY : nbnd, npwx, npw, igk  
  USE gvect,                ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, &
                                   nrx3, nrxx, nl, nlm, g, gg, ecutwfc, gcutm
  USE ions_base,            ONLY : nat, tau, atm, ityp
  USE io_global,       ONLY : stdout
  USE symme,     ONLY : s, nsym, irt
  USE pwcom
  USE gipaw_module

  ! Arguments
  real(dp), intent(inout) :: sigma_diamagnetic(3,3,nat)
  
  ! Local
  integer :: na
  real(dp) :: tr_sigma
  

  write(stdout,'(5X,''NMR chemical diamagnetic shifts in ppm:'')')
  write(stdout,*)
  
  ! symmetrize tensors
  do na = 1, nat
    call trntns (sigma_diamagnetic(1,1,na), at, bg, -1)
  enddo
  call symz(sigma_diamagnetic, nsym, s, nat, irt)
  do na = 1, nat
    call trntns (sigma_diamagnetic(1,1,na), at, bg, 1)
  enddo
  
  do na = 1, nat
    tr_sigma = (sigma_diamagnetic(1,1,na)+sigma_diamagnetic(2,2,na) &
         +sigma_diamagnetic(3,3,na))/3.0_dp
    write(stdout, &
      '(5X,''Atom'',I3,2X,A3,'' pos: ('',3(F10.6),'')  sigma: '',F14.4)') &
       na, atm(ityp(na)), tau(:,na), tr_sigma*1d6
    write(stdout, tens_fmt) sigma_diamagnetic(:,:,na) * 1d6
  enddo

end subroutine compute_sigma_diamagnetic

!-----------------------------------------------------------------------
SUBROUTINE compute_sigma_paramagnetic( sigma_paramagnetic )
  !-----------------------------------------------------------------------
  !
  ! ... Compute the paramagnetic contribution to the chemical shift at the
  ! ... position of the nuclei, given the induced field
  USE kinds,                ONLY : DP
  USE klist,                ONLY : xk
  USE wvfct,                ONLY : nbnd, npwx, npw, igk  
  USE gvect,                ONLY : ngm, gstart, nr1, nr2, nr3, nrx1, nrx2, &
                                   nrx3, nrxx, nl, nlm, g, gg, ecutwfc, gcutm
  USE ions_base,            ONLY : nat, tau, atm, ityp
  USE io_global,       ONLY : stdout
  USE symme,     ONLY : s, nsym, irt
  USE pwcom
  USE gipaw_module

  ! Arguments
  real(dp), intent(inout) :: sigma_paramagnetic(3,3,nat)
  
  ! Local
  integer :: na
  real(dp) :: tr_sigma

  write(stdout,'(5X,''NMR chemical paramagnetic shifts in ppm:'')')
  write(stdout,*)
  
  ! symmetrize tensors
  do na = 1, nat
    call trntns (sigma_paramagnetic(1,1,na), at, bg, -1)
  enddo
  call symz(sigma_paramagnetic, nsym, s, nat, irt)
  do na = 1, nat
    call trntns (sigma_paramagnetic(1,1,na), at, bg, 1)
  enddo
  
  do na = 1, nat
    tr_sigma = (sigma_paramagnetic(1,1,na)+sigma_paramagnetic(2,2,na) &
         +sigma_paramagnetic(3,3,na))/3.0_dp
    write(stdout, &
      '(5X,''Atom'',I3,2X,A3,'' pos: ('',3(F10.6),'')  sigma: '',F14.4)') &
       na, atm(ityp(na)), tau(:,na), tr_sigma*1d6
    write(stdout, tens_fmt) sigma_paramagnetic(:,:,na) * 1d6
  enddo

end subroutine compute_sigma_paramagnetic
