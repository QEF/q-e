!
! Copyright (C) 2018 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE plot_sphere_bspline(nr, lebedev, m1, x0, rhor, alat, ounit, laue)
  !--------------------------------------------------------------------------
  !
  ! ... Use B-spline interpolation instead of Fourier,
  ! ... and calculate spherical average with Lebedev Quadrature.
  !
  USE constants, ONLY : fpi
  USE io_global, ONLY : stdout, ionode
  USE kinds,     ONLY : DP
  USE fft_base,  ONLY : dfftp
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: nr
  INTEGER,  INTENT(IN) :: lebedev
  REAL(DP), INTENT(IN) :: m1
  REAL(DP), INTENT(IN) :: x0(3)
  REAL(DP), INTENT(IN) :: rhor(dfftp%nr1x, dfftp%nr2x, dfftp%nr3x)
  REAL(DP), INTENT(IN) :: alat
  INTEGER,  INTENT(IN) :: ounit
  LOGICAL,  INTENT(IN) :: laue
  !
  INTEGER  :: ir
  INTEGER  :: il
  INTEGER  :: itot
  INTEGER  :: ilebedev
  INTEGER  :: nlebedev
  INTEGER  :: lebedev_
  REAL(DP) :: weight
  REAL(DP) :: r, dr
  REAL(DP) :: e(3)
  REAL(DP) :: deltar
  REAL(DP) :: rhoint
  !
  REAL(DP), ALLOCATABLE :: lg(:,:) ! Lebedev grid
  REAL(DP), ALLOCATABLE :: rg(:,:) ! (x,y,z) on total grid
  REAL(DP), ALLOCATABLE :: vg(:)   ! value on total grid
  REAL(DP), ALLOCATABLE :: carica(:)
  !
  include "lebedev.f90"
  !
  ! ... decide order of Lebedev grid
  lebedev_ = -1
  nlebedev = SIZE(LEBEDEV_ORDER)
  !
  DO ilebedev = 1, nlebedev
    IF (lebedev <= LEBEDEV_ORDER(ilebedev)) THEN
      lebedev_ = LEBEDEV_ORDER(ilebedev)
      EXIT
    END IF
  END DO
  !
  IF (lebedev_ < 1) THEN
    lebedev_ = LEBEDEV_ORDER(nlebedev)
  END IF
  !
  WRITE(stdout, '(5X, "#Grid of Lebedev: ", I5)') lebedev_
  !
  ! ... allocate data
  ALLOCATE(lg(4, lebedev_))
  ALLOCATE(rg(3, nr * lebedev_))
  ALLOCATE(vg(   nr * lebedev_))
  ALLOCATE(carica(nr))
  !
  ! ... set Lebedev grid
  CALL set_lebedev_grid(lebedev_)
  !
  ! ... grid in cartesian coordinates, in units of alat
  deltar = m1 / DBLE(nr - 1)
  !
  DO ir = 1, nr
    r = DBLE(ir - 1) * deltar
    !
    DO il = 1, lebedev_
      e(1:3) = lg(1:3, il)
      !
      itot = il + (ir - 1) * lebedev_
      !
      rg(1, itot) = x0(1) + r * e(1)
      rg(2, itot) = x0(2) + r * e(2)
      rg(3, itot) = x0(3) + r * e(3)
    END DO
  END DO
  !
  ! ... interpolate
  CALL bspline_interpolation(nr * lebedev_, rg, rhor, vg, laue)
  !
  ! ... calculate spherical average
  DO ir = 1, nr
    carica(ir) = 0.0_DP
    !
    DO il = 1, lebedev_
      weight = fpi * lg(4, il)
      !
      itot = il + (ir - 1) * lebedev_
      !
      carica(ir) = carica(ir) + vg(itot) * weight
    END DO
  END DO
  !
  ! ... we print the charge on output
  !
  IF (ionode) THEN
    !
    dr = deltar * alat
    rhoint = 0.0_DP
    !
    DO ir = 1, nr
      !
      ! ... simple trapezoidal rule: rhoint=int carica(i) r^2(i) dr
      !
      r  = DBLE(ir - 1) * deltar * alat
      rhoint = rhoint + carica(ir) * r * r * dr
      !
      WRITE(ounit, '(3F20.10)') r / alat, carica(ir), rhoint
    END DO
    !
  END IF
  !
  ! ... deallocate data
  DEALLOCATE(lg)
  DEALLOCATE(rg)
  DEALLOCATE(vg)
  DEALLOCATE(carica)
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE set_lebedev_grid(ileb)
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ileb
    !
    SELECT CASE (ileb)
    CASE (6)
      lg = LEBEDEV6
    CASE (14)
      lg = LEBEDEV14
    CASE (26)
      lg = LEBEDEV26
    CASE (38)
      lg = LEBEDEV38
    CASE (50)
      lg = LEBEDEV50
    CASE (74)
      lg = LEBEDEV74
    CASE (86)
      lg = LEBEDEV86
    CASE (110)
      lg = LEBEDEV110
    CASE (146)
      lg = LEBEDEV146
    CASE (170)
      lg = LEBEDEV170
    CASE (194)
      lg = LEBEDEV194
    CASE (230)
      lg = LEBEDEV230
    CASE (266)
      lg = LEBEDEV266
    CASE (302)
      lg = LEBEDEV302
    CASE (350)
      lg = LEBEDEV350
    CASE (434)
      lg = LEBEDEV434
    CASE (590)
      lg = LEBEDEV590
    CASE (770)
      lg = LEBEDEV770
    CASE (974)
      lg = LEBEDEV974
    CASE (1202)
      lg = LEBEDEV1202
    CASE DEFAULT
      CALL errore('plot_sphere_bspline', 'incorrect Lebedev grid', 1)
    END SELECT
    !
  END SUBROUTINE set_lebedev_grid
  !
END SUBROUTINE plot_sphere_bspline
