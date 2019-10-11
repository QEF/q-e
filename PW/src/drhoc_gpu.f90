!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE drhoc_gpu( ngl, gl_d, omega, tpiba2, mesh, r_d, rab_d, rhoc_d, rhocg_d )
  !-----------------------------------------------------------------------
  !
  USE kinds
  USE constants, ONLY : pi, fpi
  !
  USE simpsn_gpum,  ONLY: simpsn_gpu_dev
  USE sph_bes_gpum, ONLY : sph_bes_gpu
  USE gbuffers,     ONLY : dev_buf
  implicit none
  !
  !    first the dummy variables
  !
  integer :: ngl, mesh
  ! input: the number of g shell
  ! input: the number of radial mesh points

  real(DP) :: gl_d (ngl), r_d (mesh), rab_d (mesh), rhoc_d (mesh), omega, &
       tpiba2, rhocg_d (ngl)
#if defined(__CUDA)
  attributes(DEVICE) :: gl_d, r_d, rab_d, rhoc_d, rhocg_d
#endif
  ! input: the number of G shells
  ! input: the radial mesh
  ! input: the derivative of the radial mesh
  ! input: the radial core charge
  ! input: the volume of the unit cell
  ! input: 2 times pi / alat
  ! output: the fourier transform of the core charge
  !
  !     here the local variables
  !
  real(DP) :: gx, rhocg1
  ! the modulus of g for a given shell
  ! the fourier transform
  real(DP), contiguous, pointer :: aux_d (:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: aux_d
#endif
  ! auxiliary memory for integration

  integer :: ir, igl, ierr
  ! counter on radial mesh points
  ! counter on g shells
  ! lower limit for loop on ngl
  !
  CALL dev_buf%lock_buffer(aux_d, (/ mesh, ngl /), ierr) ! allocate (aux_d( mesh, ngl))
  !
  ! G=0 and G <> 0 term done together (in the end is just an avoided bessel call and multiply by 1.0)
  !
!$cuf kernel do(1)
  DO igl = 1, ngl
     gx = SQRT(gl_d(igl) * tpiba2)
     CALL sph_bes_gpu( mesh, r_d, gx, 0, aux_d(:,igl) )
     DO ir = 1, mesh
        aux_d(ir,igl) = r_d(ir)**2 * rhoc_d(ir) * aux_d(ir,igl)
     ENDDO
     !
     CALL simpsn_gpu_dev( mesh, aux_d(:,igl), rab_d, rhocg_d(igl) )
     rhocg_d(igl) =    fpi * rhocg_d(igl) / omega
  ENDDO
  !
  CALL dev_buf%release_buffer(aux_d, ierr)
  !
  return
end subroutine drhoc_gpu
