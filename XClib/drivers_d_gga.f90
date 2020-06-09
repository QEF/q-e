!---------------------------------------------------------------------------
SUBROUTINE dgcxc_unpol_l( length, r_in, s2_in, vrrx, vsrx, vssx, vrrc, vsrc, vssc )
  !-------------------------------------------------------------------------
  !! This routine computes the derivative of the exchange and correlation
  !! potentials.
  !
  USE kind_l,         ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  REAL(DP), INTENT(IN), DIMENSION(length) :: r_in, s2_in
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vrrx, vsrx, vssx
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vrrc, vsrc, vssc
  !
  ! ... local variables
  !
  INTEGER :: i1, i2, i3, i4, f1, f2, f3, f4
  REAL(DP), DIMENSION(length) :: dr, s, ds
  REAL(DP), DIMENSION(4*length) :: raux, s2aux
  REAL(DP), ALLOCATABLE :: v1x(:), v2x(:), v1c(:), v2c(:)
  REAL(DP), ALLOCATABLE :: sx(:), sc(:)
  REAL(DP), PARAMETER :: small = 1.E-30_DP
  !
  ALLOCATE( v1x(4*length), v2x(4*length), sx(4*length) )
  ALLOCATE( v1c(4*length), v2c(4*length), sc(4*length) )
  !
  i1 = 1     ;   f1 = length     !4 blocks:  [ rho+dr ,    grho2    ]
  i2 = f1+1  ;   f2 = 2*length   !           [ rho-dr ,    grho2    ]
  i3 = f2+1  ;   f3 = 3*length   !           [ rho    , (grho+ds)^2 ]
  i4 = f3+1  ;   f4 = 4*length   !           [ rho    , (grho-ds)^2 ]
  !
  s  = SQRT(s2_in)
  dr = MIN(1.d-4, 1.d-2*r_in)
  ds = MIN(1.d-4, 1.d-2*s)
  !
  raux(i1:f1) = r_in+dr  ;   s2aux(i1:f1) = s2_in
  raux(i2:f2) = r_in-dr  ;   s2aux(i2:f2) = s2_in
  raux(i3:f3) = r_in     ;   s2aux(i3:f3) = (s+ds)**2
  raux(i4:f4) = r_in     ;   s2aux(i4:f4) = (s-ds)**2
  !
  CALL gcxc_l( length*4, raux, s2aux, sx, sc, v1x, v2x, v1c, v2c )
  !
  ! ... to avoid NaN in the next operations
  WHERE( r_in<=small .OR. s2_in<=small )
    dr = 1._DP ; ds = 1._DP ; s = 1._DP
  END WHERE
  !
  vrrx = 0.5_DP * (v1x(i1:f1) - v1x(i2:f2)) / dr
  vrrc = 0.5_DP * (v1c(i1:f1) - v1c(i2:f2)) / dr
  !
  vsrx = 0.25_DP * ((v2x(i1:f1) - v2x(i2:f2)) / dr + &
                    (v1x(i3:f3) - v1x(i4:f4)) / ds / s)
  vsrc = 0.25_DP * ((v2c(i1:f1) - v2c(i2:f2)) / dr + &
                    (v1c(i3:f3) - v1c(i4:f4)) / ds / s)
  !
  vssx = 0.5_DP * (v2x(i3:f3) - v2x(i4:f4)) / ds / s
  vssc = 0.5_DP * (v2c(i3:f3) - v2c(i4:f4)) / ds / s
  !
  DEALLOCATE( v1x, v2x, sx )
  DEALLOCATE( v1c, v2c, sc )
  !
  RETURN
  !
END SUBROUTINE dgcxc_unpol_l
!
