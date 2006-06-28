! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE compute_sigma_bare(chi_bare)
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
  USE io_global,       ONLY : stdout
  USE symme,     ONLY : s, nsym, irt
  USE pwcom
  USE nmr_module

  REAL(DP), INTENT(IN) :: chi_bare(3,3)
  REAL(DP) :: macroscopic_shape(3,3)
  real(dp) :: sigma_bare(3,3,nat)
  complex(dp) :: tmp_sigma(3,3)
  real(dp) :: arg, tr_sigma
  integer :: na, ig

  macroscopic_shape(:,:) = 2.d0/3.d0
  ! like in paratec:
  macroscopic_shape(:,:) = 0.d0
  do na = 1, 3 
    macroscopic_shape(na,na) = 2.d0/3.d0
  enddo

  write(stdout,'(5X,''NMR chemical shifts in ppm:'')')
  write(stdout,*)
  
  do na = 1, nat
    tmp_sigma(:,:) = 0.d0

    do ig = gstart, ngm              
      arg = (g(1,ig)*tau(1,na) + g(2,ig)*tau(2,na) + g(3,ig)*tau(3,na)) * tpi
      tmp_sigma(:,:) = tmp_sigma(:,:) + b_ind(ig,:,:) * cmplx(cos(arg),sin(arg))
    enddo

#if 0
    ! this is the G = 0 term
    if (gstart == 2) &
      tmp_sigma(:,:) = tmp_sigma(:,:) + &
        (4.d0*pi) * macroscopic_shape(:,:) * chi_bare(:,:)
#endif

    sigma_bare(:,:,na) = real(tmp_sigma(:,:))
  enddo

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
    tr_sigma = (sigma_bare(1,1,na)+sigma_bare(2,2,na)+sigma_bare(3,3,na))/3.d0
    write(stdout,'(5X,''Atom'',I3,2X,A3,'' pos: ('',3(F10.6),&
          '')  sigma: '',F14.4)') na, atm(ityp(na)), tau(:,na), tr_sigma*1d6

    write(stdout, '(3(5X,3(F14.4,2X)/))') sigma_bare(:,:,na) * 1d6
    !!!call sym_cart_tensor(sigma_bare(:,:,na))
    !!!write(stdout, '(3(5X,3(F12.6,2X)/))') sigma_bare(:,:,na) * 1d6
  enddo

end subroutine compute_sigma_bare



