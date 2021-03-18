!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE simpson_gpu(mesh, func_d, rab_d, asum)
  !-----------------------------------------------------------------------
  !
  !     simpson's rule integration. On input:
  !       mesh = the number of grid points (should be odd)
  !       func(i)= function to be integrated
  !       rab(i) = r(i) * dr(i)/di * di
  !     For the logarithmic grid not including r=0 :
  !       r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
  !     For the logarithmic grid including r=0 :
  !       r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
  !     Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
  !     where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, INTENT(in) :: mesh
  real(DP), INTENT(in) :: rab_d (mesh), func_d (mesh)
  real(DP), INTENT(out):: asum
#if defined(__CUDA)
  attributes(DEVICE) :: rab_d, func_d
#endif
  !
  real(DP) :: f1, f2, f3, r12
  INTEGER :: i
  !
  asum = 0.0d0
  r12 = 1.0d0 / 3.0d0

  !$cuf kernel do (1) <<<*, *>>>
  do i = 2, mesh-1, 2
    asum = asum + func_d(i-1)*rab_d(i-1) + 4.0d0*func_d(i)*rab_d(i) + func_d(i+1)*rab_d(i+1)
  end do

  asum = asum*r12
  !
  ! if mesh is not odd, use open formula instead:
  ! ... 2/3*f(n-5) + 4/3*f(n-4) + 13/12*f(n-3) + 0*f(n-2) + 27/12*f(n-1)
  !!! Under testing
  !
  !IF ( MOD(mesh,2) == 0 ) THEN
  !   print *, 'mesh even: correction:', f1*5.d0/4.d0-4.d0*f2+23.d0*f3/4.d0, &
  !                                      func(mesh)*rab(mesh), asum
  !   asum = asum + f1*5.d0/4.d0 - 4.d0*f2 + 23.d0*f3/4.d0
  !END IF

  RETURN
END SUBROUTINE simpson_gpu

module simpsn_gpum
contains
#if defined(__CUDA)
attributes(DEVICE) &
#endif
SUBROUTINE simpsn_gpu_dev(mesh, func_d, rab_d, asum)
  !-----------------------------------------------------------------------
  !
  !     simpson's rule integration. On input:
  !       mesh = the number of grid points (should be odd)
  !       func(i)= function to be integrated
  !       rab(i) = r(i) * dr(i)/di * di
  !     For the logarithmic grid not including r=0 :
  !       r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
  !     For the logarithmic grid including r=0 :
  !       r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
  !     Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
  !     where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, VALUE :: mesh
  real(DP) :: rab_d (mesh), func_d (mesh)
  real(DP) :: asum
  !
  real(DP) :: f1, f2, f3, r12
  INTEGER :: i
  !
  asum = 0.0d0
  r12 = 1.0d0 / 3.0d0

  do i = 2, mesh-1, 2
    asum = asum + func_d(i-1)*rab_d(i-1) + 4.0d0*func_d(i)*rab_d(i) + func_d(i+1)*rab_d(i+1)
  end do

  asum = asum*r12
  !
  ! if mesh is not odd, use open formula instead:
  ! ... 2/3*f(n-5) + 4/3*f(n-4) + 13/12*f(n-3) + 0*f(n-2) + 27/12*f(n-1)
  !!! Under testing
  !
  !IF ( MOD(mesh,2) == 0 ) THEN
  !   print *, 'mesh even: correction:', f1*5.d0/4.d0-4.d0*f2+23.d0*f3/4.d0, &
  !                                      func(mesh)*rab(mesh), asum
  !   asum = asum + f1*5.d0/4.d0 - 4.d0*f2 + 23.d0*f3/4.d0
  !END IF

  RETURN
END SUBROUTINE simpsn_gpu_dev
end module
