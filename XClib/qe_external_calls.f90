!
!------------------------------------------------------------------------
SUBROUTINE lsd_glyp_ext( rho_in_up, rho_in_dw, grho_up, grho_dw, grho_ud, sc, v1c_up, v1c_dw, v2c_up, v2c_dw, v2c_ud )                     !<GPU:DEVICE>
  !----------------------------------------------------------------------
  !! Lee, Yang, Parr: gradient correction part.
  !
  USE kind_l, ONLY: DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rho_in_up, rho_in_dw
  !! the total charge density
  REAL(DP), INTENT(IN) :: grho_up, grho_dw
  !! the gradient of the charge squared
  REAL(DP), INTENT(IN) :: grho_ud
  !! gradient off-diagonal term up-down
  REAL(DP), INTENT(OUT) :: sc
  !! correlation energy
  REAL(DP), INTENT(OUT) :: v1c_up, v1c_dw
  !! derivative of correlation wr. rho
  REAL(DP), INTENT(OUT) :: v2c_up, v2c_dw
  !! derivatives of correlation wr. grho
  REAL(DP), INTENT(OUT) :: v2c_ud
  !! derivative of correlation wr. grho, off-diag. term
  !
  ! ... local variables
  !
  REAL(DP) :: ra, rb, rho, grhoaa, grhoab, grhobb
  REAL(DP) :: rm3, dr, or, dor, der, dder
  REAL(DP) :: dlaa, dlab, dlbb, dlaaa, dlaab, dlaba, dlabb, dlbba, dlbbb
  REAL(DP), PARAMETER :: a=0.04918_DP, b=0.132_DP, c=0.2533_DP, d=0.349_DP
  !
  ra = rho_in_up
  rb = rho_in_dw
  rho = ra + rb
  rm3 = rho**(-1._DP/3._DP)
  !
  dr = ( 1._DP + d*rm3 )
  or = EXP(-c*rm3) / dr * rm3**11._DP
  dor = -1._DP/3._DP * rm3**4 * or * (11._DP/rm3-c-d/dr)
  !
  der  = c*rm3 + d*rm3/dr
  dder = 1._DP/3._DP * (d*d*rm3**5/dr/dr - der/rho)
  !
  dlaa = -a*b*or*( ra*rb/9._DP*(1._DP-3*der-(der-11._DP)*ra/rho)-rb*rb )
  dlab = -a*b*or*( ra*rb/9._DP*(47._DP-7._DP*der)-4._DP/3._DP*rho*rho  )
  dlbb = -a*b*or*( ra*rb/9._DP*(1._DP-3*der-(der-11._DP)*rb/rho)-ra*ra )
  !
  dlaaa = dor/or*dlaa-a*b*or*(rb/9._DP*(1._DP-3*der-(der-11._DP)*ra/rho)-     &
          ra*rb/9._DP*((3._DP+ra/rho)*dder+(der-11._DP)*rb/rho/rho))
  dlaab = dor/or*dlaa-a*b*or*(ra/9._DP*(1._DP-3._DP*der-(der-11._DP)*ra/rho)- &
          ra*rb/9._DP*((3._DP+ra/rho)*dder-(der-11._DP)*ra/rho/rho)-2._DP*rb)
  dlaba = dor/or*dlab-a*b*or*(rb/9._DP*(47._DP-7._DP*der)-7._DP/9._DP*ra*rb*dder- &
          8._DP/3._DP*rho)
  dlabb = dor/or*dlab-a*b*or*(ra/9._DP*(47._DP-7._DP*der)-7._DP/9._DP*ra*rb*dder- &
          8._DP/3._DP*rho)
  dlbba = dor/or*dlbb-a*b*or*(rb/9._DP*(1._DP-3._DP*der-(der-11._DP)*rb/rho)- &
          ra*rb/9._DP*((3._DP+rb/rho)*dder-(der-11._DP)*rb/rho/rho)-2._DP*ra)
  dlbbb = dor/or*dlbb-a*b*or*(ra/9._DP*(1._DP-3*der-(der-11._DP)*rb/rho)-     &
          ra*rb/9._DP*((3._DP+rb/rho)*dder+(der-11._DP)*ra/rho/rho))
  !
  grhoaa = grho_up
  grhobb = grho_dw
  grhoab = grho_ud
  !
  sc     = dlaa *grhoaa + dlab *grhoab + dlbb *grhobb
  v1c_up = dlaaa*grhoaa + dlaba*grhoab + dlbba*grhobb
  v1c_dw = dlaab*grhoaa + dlabb*grhoab + dlbbb*grhobb
  v2c_up = 2._DP*dlaa
  v2c_dw = 2._DP*dlbb
  v2c_ud = dlab
  !
  !
  RETURN
  !
END SUBROUTINE lsd_glyp_ext

!-------------------------------------------------------------------------
SUBROUTINE tpsscxc_ext( rho, grho, tau, sx, sc, v1x, v2x, v3x, v1c, v2c, v3c )                    !<GPU:DEVICE>
  !-----------------------------------------------------------------------
  !! TPSS metaGGA corrections for exchange and correlation - Hartree a.u.
  !
  !! Definition:  E_x = \int E_x(rho,grho) dr
  !
  USE kind_l,      ONLY : DP
  USE metagga,     ONLY : metax, metac
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rho
  !! the charge density
  REAL(DP), INTENT(IN) :: grho
  !! grho = |\nabla rho|^2
  REAL(DP), INTENT(IN) :: tau
  !! kinetic energy density
  REAL(DP), INTENT(OUT) :: sx
  !! sx = E_x(rho,grho)
  REAL(DP), INTENT(OUT) :: sc
  !! sc = E_c(rho,grho)
  REAL(DP), INTENT(OUT) :: v1x
  !! v1x = D(E_x)/D(rho)
  REAL(DP), INTENT(OUT) :: v2x
  !! v2x = D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  REAL(DP), INTENT(OUT) :: v3x
  !! v3x = D(E_x)/D(tau)
  REAL(DP), INTENT(OUT) :: v1c
  !! v1c = D(E_c)/D(rho)
  REAL(DP), INTENT(OUT) :: v2c
  !! v2c = D(E_c)/D( D rho/D r_alpha ) / |\nabla rho|
  REAL(DP), INTENT(OUT) :: v3c
  !! v3c = D(E_c)/D(tau)
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: small = 1.E-10_DP
  !
  IF (rho <= small) THEN
     sx  = 0.0_DP
     v1x = 0.0_DP
     v2x = 0.0_DP
     sc  = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
     v3x = 0.0_DP
     v3c = 0.0_DP
     RETURN
  ENDIF
  !
  ! exchange
  CALL metax( rho, grho, tau, sx, v1x, v2x, v3x )                   !<GPU:metax=>metax_d>
  ! correlation
  CALL metac( rho, grho, tau, sc, v1c, v2c, v3c )                   !<GPU:metac=>metac_d>
  !
  RETURN
  !
END SUBROUTINE tpsscxc_ext
