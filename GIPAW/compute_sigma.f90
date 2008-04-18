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
  USE pwcom,                ONLY : pi, tpi
  USE gipaw_module,         ONLY : use_nmr_macroscopic_shape, &
                                   nmr_macroscopic_shape, b_ind, tens_fmt
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum

  IMPLICIT NONE

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
  call mp_sum( sigma_bare, intra_pool_comm )
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
    write(stdout,'(5X,"Atom",I3,2X,A3," pos: (",3(F10.6),")  sigma: ",F14.4)') &
       na, atm(ityp(na)), tau(:,na), tr_sigma*1e6_dp
    write(stdout, tens_fmt) sigma_bare(:,:,na) * 1e6_dp
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
  USE io_global,            ONLY : stdout
  USE symme,                ONLY : s, nsym, irt
  USE pwcom,                ONLY : at, bg
  USE gipaw_module,         ONLY : tens_fmt
  IMPLICIT NONE

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
    write(stdout,'(5X,"Atom",I3,2X,A3," pos: (",3(F10.6),")  sigma: ",F14.4)') &
        na, atm(ityp(na)), tau(:,na), tr_sigma*1e6_dp
    write(stdout, tens_fmt) sigma_diamagnetic(:,:,na) * 1e6_dp
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
  USE io_global,            ONLY : stdout
  USE symme,                ONLY : s, nsym, irt
  USE pwcom,                ONLY : at, bg
  USE gipaw_module,         ONLY : tens_fmt
  IMPLICIT NONE

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
    write(stdout,'(5X,"Atom",I3,2X,A3," pos: (",3(F10.6),")  sigma: ",F14.4)') &
       na, atm(ityp(na)), tau(:,na), tr_sigma*1e6_dp
    write(stdout, tens_fmt) sigma_paramagnetic(:,:,na) * 1e6_dp
  enddo

end subroutine compute_sigma_paramagnetic



!-----------------------------------------------------------------------
SUBROUTINE print_sigma_total(sigma_bare, sigma_paramagnetic, sigma_diamagnetic)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, tau, atm, ityp
  USE io_global,            ONLY : stdout
  USE gipaw_module,         ONLY : tens_fmt, iverbosity, nmr_shift_core
  IMPLICIT NONE

  ! Arguments
  REAL(dp), INTENT(IN) :: sigma_bare(3,3,nat)
  REAL(dp), INTENT(IN) :: sigma_paramagnetic(3,3,nat)
  REAL(dp), INTENT(IN) :: sigma_diamagnetic(3,3,nat)
  
  ! Local
  INTEGER :: na, i, ordering(3)
  REAL(dp) :: tmp(3,3), tr_sigma, nmr_eig(3), nmr_vect(3,3)
  REAL(dp) :: v(3), eta, tmp2(3,3)
  
  IF ( iverbosity > 5 ) THEN
     write(stdout,*)
     write(stdout,'(5X,''Summed isotropic NMR chemical shifts in ppm:'')')
     write(stdout,*)
     do na = 1, nat
        tmp(:,:) = sigma_bare(:,:,na) &
             + sigma_paramagnetic(:,:,na) &
             + sigma_diamagnetic(:,:,na)
        tr_sigma = (tmp(1,1) + tmp(2,2) + tmp(3,3))/3.0_dp
        write(stdout,'(5X,"Atom",I3,2X,A3," pos: (",3F10.6,")  sigma: ",F14.4)')&
             na, atm(ityp(na)), tau(:,na), tr_sigma*1e6_dp
        write(stdout, tens_fmt) tmp(:,:) * 1e6_dp
        write(stdout,*)
        write(stdout,*)
     end do
  END IF
  
  write(stdout,*)
  write(stdout,'(5X,''Total NMR chemical shifts in ppm:'')')
  write(stdout,*)
  
  do na = 1, nat
     tmp(:,:) = sigma_bare(:,:,na) &
          + sigma_paramagnetic(:,:,na) &
          + sigma_diamagnetic(:,:,na)
     tmp(1,1) = tmp(1,1) + nmr_shift_core(ityp(na))
     tmp(2,2) = tmp(2,2) + nmr_shift_core(ityp(na))
     tmp(3,3) = tmp(3,3) + nmr_shift_core(ityp(na))
     
     tr_sigma = ( tmp(1,1) + tmp(2,2) + tmp(3,3) ) / 3.0_dp
     write(stdout,'(5X,"Atom",I3,2X,A3," pos: (",3(F10.6),")  sigma: ",F14.4)')&
          na, atm(ityp(na)), tau(:,na), tr_sigma * 1e6_dp
     write(stdout, tens_fmt) tmp(:,:) * 1e6_dp
     
!     if ( abs ( v(1) ) < 1e-5 ) then
!        eta = 0.0_dp
!     else
!        eta = ( v(2) - v(3) ) / v(1)
!     end if
!     write ( stdout, '(/,5X,"eigenvalues: ",3(F14.4), / )' ), v(1:3) * 1e6
!     write ( stdout, '(/,5X,"eigenvalues: ",3(F14.4),/,5X, "anisotropy: ",F12.6, / )' ), v(1:3) * 1e6, eta
     
     tmp2 = 0.5d0 * ( tmp + TRANSPOSE ( tmp ) )
     
     !
     ! diagonalise the tensor to extract the eigenvalues and anisotropy
     !
     call rdiagh ( 3, tmp2, 3, nmr_eig, nmr_vect )
     
     v(2) = nmr_eig(2)
     ordering(2) = 2
     if ( abs(nmr_eig(1)) > abs(nmr_eig(3)) ) then
        v(1) = nmr_eig(1)
        v(3) = nmr_eig(3)
        ordering(1) = 1
        ordering(3) = 3
     else
        v(1) = nmr_eig(3)
        v(3) = nmr_eig(1)
        ordering(1) = 3
        ordering(3) = 1
     end if
     
     write ( stdout, '( 5X, "Symmetric tensor" )' )
     write(stdout, tens_fmt) tmp2 * 1e6_dp
     
     DO i = 1, 3
        write ( stdout, '(5X,"eigenvalue:  ", F14.4 )' ) & 
             v(i) * 1e6_dp
        write ( stdout, '(5X,"eigenvector: ",3(F14.4)/ )' ) &
             nmr_vect(:,ordering(i))
     END DO
     
     write ( stdout, '( 5X, "Anti-symmetric tensor" )' )
     write(stdout, tens_fmt) &
          0.5d0 * ( tmp(:,:) - TRANSPOSE ( tmp(:,:) ) ) * 1e6_dp
     
     write ( stdout, '( / )' )
     
  end do
  
  write(stdout,*)
  write(stdout,*)
  
end subroutine print_sigma_total
