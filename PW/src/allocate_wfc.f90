!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE allocate_wfc()
  !----------------------------------------------------------------------------
  !! Dynamical allocation of arrays: wavefunctions.  
  !! Requires dimensions: \(\text{npwx}\), \(\text{nbnd}\), \(\text{npol}\)
  !
#if defined (__CUDA)
  use, intrinsic :: iso_c_binding
  use cudafor
#endif
  USE io_global,           ONLY : stdout
  USE wvfct,               ONLY : npwx, nbnd
  USE noncollin_module,    ONLY : npol
  USE wavefunctions,       ONLY : evc
  USE control_flags,       ONLY : use_gpu
  !
  IMPLICIT NONE
  INTEGER :: istat
  !
  !
  ALLOCATE( evc(npwx*npol,nbnd) )
!civn: PIN evc memory here
#if defined(__CUDA)
  IF(use_gpu) istat = cudaHostRegister(C_LOC(evc(1,1)), sizeof(evc), cudaHostRegisterMapped)
  !$acc enter data create(evc)
#endif
  !
  RETURN
  !
END SUBROUTINE allocate_wfc
!
!----------------------------------------------------------------------------
SUBROUTINE check_wfc( label, what )
  !----------------------------------------------------------------------------
  !! Check consistency between host and device copies of evc 
  !! by computing the scalar product with matcalc_k
  !
  USE kinds,               ONLY :  DP
  USE wavefunctions,       ONLY : evc
  USE wvfct,               ONLY : npwx, nbnd
  USE noncollin_module,    ONLY : npol
  !
  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: label
  !! it specifies the meaning of the output
  CHARACTER(len=2), INTENT(IN) :: what 
  !! it specifies the meaning of the output
  COMPLEX(DP), ALLOCATABLE :: evc_copy(:,:)
  !
  ALLOCATE( evc_copy(npwx*npol,nbnd) )
  !
  IF(what == 'HH') THEN
    evc_copy(:,:) = evc(:,:)
  ELSE IF(what == 'DH') THEN
    !$acc kernels copyout(evc_copy)
    evc_copy(:,:) = evc(:,:)
    !$acc end kernels
  ELSE
    Call errore('check_wfc', 'wrong what input value', 1)
  END IF
  !
  ! compute the scalar product between host and device values of evc
  ! and check orthonormality
  Call wrapmatcalc( label , npwx*npol, nbnd, nbnd, evc_copy, evc)  
  !
  DEALLOCATE( evc_copy )
  !
  RETURN       
  !
END SUBROUTINE check_wfc
!
!----------------------------------------------------------------------------
SUBROUTINE allocate_wfc_k()
  !----------------------------------------------------------------------------
  !! Dynamical allocation of k-point-dependent arrays: wavefunctions, betas
  !! kinetic energy, k+G indices. Computes max no. of plane waves \(\text{npwx}\)
  !! and k+G indices \(\text{igk_k}\) (needs G-vectors and cutoff \(\text{gcutw}\)).  
  !! Requires dimensions \(\text{nbnd}\), \(\text{npol}\), \(\text{natomwfc}\),
  !! \(\text{nwfcU}\).  
  !! Requires that k-points are set up and distributed (if parallelized).
  !
  USE wvfct,            ONLY : npwx, g2kin
  USE uspp,             ONLY : vkb, nkb
  USE gvecw,            ONLY : gcutw
  USE gvect,            ONLY : ngm, g
  USE klist,            ONLY : xk, nks, init_igk
  !
  IMPLICIT NONE
  !
  INTEGER, EXTERNAL :: n_plane_waves
  !
  !   calculate number of PWs for all kpoints
  !
  npwx = n_plane_waves( gcutw, nks, xk, g, ngm )
  !
  !   compute indices j=igk(i) such that (k+G)_i = k+G_j, for all k
  !   compute number of plane waves ngk(ik) as well
  !
  CALL init_igk( npwx, ngm, g, gcutw )
  !
  CALL allocate_wfc()
  !
  !   beta functions
  !
  ALLOCATE( vkb(npwx,nkb) )
  !$acc enter data create(vkb) 
  !
  !   g2kin contains the kinetic energy \hbar^2(k+G)^2/2m
  !
  ALLOCATE( g2kin(npwx) )
  !$acc enter data create(g2kin) 
  !
  !
  RETURN
  !
END SUBROUTINE allocate_wfc_k
