!
! Copyright (C) 2001-2020 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine scale_sym_ops (nsym, s, ft, nr1, nr2, nr3, s_scaled, ftau)
  !----------------------------------------------------------------------
  !
  !    Generate rotation matrices and fractional translations in scaled
  !    crystallographic axis from symmetries, check their consistency.
  !    For use in real-space symmetrization routine "rotate_grud_point"
  !
  use kinds, only : dp
  implicit none
  !
  integer, intent(in) :: nsym
  !! number of symmetries
  integer, intent(in) :: s (3,3, nsym)
  !! rotation matrices in crystallographic axis
  real(dp), intent(in) :: ft(3, nsym)
  !! fractional translations in crystallographic axis
  integer, intent(in) :: nr1
  !! FFT dimension along axis 1
  integer, intent(in) :: nr2
  !! FFT dimension along axis 2
  integer, intent(in) :: nr3
  !! FFT dimension along axis 3
  integer, intent(out) :: s_scaled (3,3, nsym)
  !! rotation matrices in scaled crystallographic axis
  integer, intent(out) :: ftau (3, nsym)
  !! rotation matrices in scaled crystallographic axis
  !
  REAL(dp), PARAMETER :: eps2 = 1.0d-5
  REAL(dp) :: ft_(3)
  integer :: isym, i, j, bad
  !
  bad = 0
  DO isym = 1, nsym
     !
     ! check if rotation sends the FFT grid into itself
     !
     IF ( MOD( s(2,1,isym)*nr1, nr2) /= 0 .OR. &
          MOD( s(3,1,isym)*nr1, nr3) /= 0 .OR. &
          MOD( s(1,2,isym)*nr2, nr1) /= 0 .OR. &
          MOD( s(3,2,isym)*nr2, nr3) /= 0 .OR. &
          MOD( s(1,3,isym)*nr3, nr1) /= 0 .OR. &
          MOD( s(2,3,isym)*nr3, nr2) /= 0 ) THEN
        call infomsg('scale_sym_ops','found rotation not compatible with FFT grid')
        bad = bad + 1
     ENDIF
     !
     s_scaled (1,1,isym) = s (1,1,isym)
     s_scaled (2,1,isym) = s (2,1,isym) * nr1 / nr2
     s_scaled (3,1,isym) = s (3,1,isym) * nr1 / nr3
     s_scaled (1,2,isym) = s (1,2,isym) * nr2 / nr1
     s_scaled (2,2,isym) = s (2,2,isym)
     s_scaled (3,2,isym) = s (3,2,isym) * nr2 / nr3
     s_scaled (1,3,isym) = s (1,3,isym) * nr3 / nr1
     s_scaled (2,3,isym) = s (2,3,isym) * nr3 / nr2
     s_scaled (3,3,isym) = s (3,3,isym)
     !
     ft_(1) = ft(1,isym) * nr1
     ft_(2) = ft(2,isym) * nr2
     ft_(3) = ft(3,isym) * nr3
     !
     ! check if the fractional translations are commensurate
     ! with the FFT grid
     IF ( ANY ( ABS(ft_(:)-NINT(ft_(:))) > eps2 ) ) THEN
        call infomsg('scale_sym_ops','found fractional translation not compatible with FFT grid')
        bad = bad + 1
     ENDIF
     !
     ftau(:,isym) = NINT(ft_(:))
     !
  ENDDO
  IF (bad > 0) CALL errore('scale_sym_ops','incompatible FFT grid',bad)
  !
end subroutine scale_sym_ops

!----------------------------------------------------------------------
subroutine rotate_grid_point (s_scaled, ftau, i, j, k, nr1, nr2, nr3, &
     ri, rj, rk)
  !----------------------------------------------------------------------
  !
  !    This routine finds the symmetry-rotated point (ri, rj, rk)
  !    from a point (i, j, k) of the FFT real-space grid.
  !    The symmetry operation can include a fractional translation
  !
  implicit none
  !
  integer, intent(in) :: s_scaled (3, 3)
  !! rotation matrices in scaled crystallographic axis
  integer, intent(in) :: ftau (3)
  !! fractional translations in scaled crystallographic axis
  integer, intent(in) :: i
  !! point to rotate along axis 1:  i=1, ... , nr1
  integer, intent(in) :: j
  !! as above for axis 2
  integer, intent(in) :: k
  !! as above for axis 3
  integer, intent(in) :: nr1
  !! FFT dimension along axis 1
  integer, intent(in) :: nr2
  !! FFT dimension along axis 2
  integer, intent(in) :: nr3
  !! FFT dimension along axis 3
  integer, intent(out) :: ri
  !! ritated point along axis 1:  ri=1, ... , nr1
  integer, intent(out) :: rj
  !! as above for axis 2
  integer, intent(out) :: rk
  !! as above for axis 3
  !
  ri = s_scaled (1, 1) * (i-1) + &
       s_scaled (2, 1) * (j-1) + &
       s_scaled (3, 1) * (k-1) - ftau (1)
  ri = mod (ri, nr1) + 1
  if (ri < 1) ri = ri + nr1
  rj = s_scaled (1, 2) * (i-1) + &
       s_scaled (2, 2) * (j-1) + &
       s_scaled (3, 2) * (k-1) - ftau (2)
  rj = mod (rj, nr2) + 1
  if (rj < 1) rj = rj + nr2
  rk = s_scaled (1, 3) * (i-1) + &
       s_scaled (2, 3) * (j-1) + &
       s_scaled (3, 3) * (k-1) - ftau (3)
  rk = mod (rk, nr3) + 1
  if (rk < 1) rk = rk + nr3
  !
end subroutine rotate_grid_point
