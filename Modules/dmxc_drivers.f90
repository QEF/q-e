!
! Copyright (C) 2004-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE dmxc( length, rho, dmuxc )
  !---------------------------------------------------------------------
  !! Computes the derivative of the xc potential with respect to the 
  !! local density.
  !
  USE xc_lda_lsda,  ONLY: xc, iexch_l, icorr_l
  USE funct,        ONLY: init_lda_xc
  USE kinds,        ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the input/output arrays
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho
  !! the charge density ( positive )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: dmuxc
  !! the derivative of the xc potential
  !
  ! ... local variables
  !
  REAL(DP), DIMENSION(length) :: ex, vx
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: rs, dr
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: ec, vc
  !
  REAL(DP), EXTERNAL :: dpz
  INTEGER :: iflg, ir
  !
  REAL(DP), PARAMETER :: small = 1.E-30_DP, e2 = 2.0_DP, &
       pi34 = 0.75_DP / 3.141592653589793_DP, third = 1.0_DP /3.0_DP
  !
  CALL init_lda_xc()
  !
  dmuxc = 0.0_DP
  !
  !    first case: analytical derivatives available
  !
  IF (iexch_l == 1 .AND. icorr_l == 1) THEN
     !
     ALLOCATE( rs(length) )
     !
     rs = 1.0_DP
     WHERE ( rho > small ) rs = (pi34 / rho)**third
     !
     !..exchange
     CALL slater( length, rs, ex, vx )
     dmuxc = vx / (3.0_DP * rho)
     !
     !..correlation
     DO ir = 1, length
        iflg = 2
        if (rs(ir) < 1.0_DP) iflg = 1
        dmuxc(ir) = dmuxc(ir) + dpz(rs(ir), iflg)
     ENDDO
     !
     WHERE ( rho < small ) dmuxc = 0._DP
     !
     DEALLOCATE( rs )
     !
  ELSE
     !
     !     second case: numerical derivatives
     !
     ALLOCATE( ec(length), vc(length) )
     ALLOCATE( dr(length) )
     !
     DO ir = 1, length
        !if (rho(ir) > small) then  !... cut already inside xc()
           dr(ir) = MIN( 1.E-6_DP, 1.E-4_DP * rho(ir) )
        !else
        !   dr(ir) = -0.5_dp 
        !endif
     ENDDO
     !
     CALL xc( length, rho+dr, ex, ec, vx, vc )
     !
     dmuxc = vx + vc
     !
     CALL xc( length, rho-dr, ex, ec, vx, vc )
     !
     dmuxc = (dmuxc - vx - vc) / (2.0_DP * dr)
     !
     DEALLOCATE( ec, vc )
     DEALLOCATE( dr )
     !
     ! where (rho<small) dmuxc=0._dp
     !
  ENDIF
  !
  ! bring to rydberg units
  !
  dmuxc = e2 * dmuxc
  !
  RETURN
  !
END SUBROUTINE dmxc
!
!
!-----------------------------------------------------------------------
SUBROUTINE dmxc_spin( length, rho, dmuxc )
!-----------------------------------------------------------------------
  !! Computes the derivative of the xc potential with respect to the 
  !! local density in the spin-polarized case.
  !
  USE xc_lda_lsda,    ONLY: xc_spin, iexch_l, icorr_l
  USE funct,          ONLY: init_lda_xc
  USE kinds,          ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the input/output arrays
  REAL(DP), INTENT(IN), DIMENSION(length,2) :: rho
  !! spin-up and spin-down charge density
  REAL(DP), INTENT(OUT), DIMENSION(length,2,2) :: dmuxc
  !! u-u, u-d, d-u, d-d derivatives of the XC functional
  !
  ! ... local variables
  !
  REAL(DP), DIMENSION(length) :: zeta, rhotot
  REAL(DP), DIMENSION(length) :: aux1, null_v
  !
  REAL(DP), ALLOCATABLE, DIMENSION(:)   :: rs, zeta_eff, ds
  REAL(DP), ALLOCATABLE, DIMENSION(:)   :: ecu, ecp, aux2
  REAL(DP), ALLOCATABLE, DIMENSION(:)   :: vcu, vcp
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: vx, vc
  !
  REAL(DP) :: fz, fz1, fz2, dmcu, dmcp, aa, bb, cc
  !
  REAL(DP), EXTERNAL :: dpz, dpz_polarized
  !
  INTEGER :: ir, is, iflg
  !
  REAL(DP), PARAMETER :: small = 1.E-30_DP, e2 = 2.0_DP, &
       pi34 = 0.75_DP / 3.141592653589793_DP, third = 1.0_DP/3.0_DP, &
       p43 = 4.0_DP / 3.0_DP, p49 = 4.0_DP / 9.0_DP, m23 = -2.0_DP / 3.0_DP
  !
  dmuxc  = 0.0_DP
  null_v = 1.0_DP
  zeta   = 0.5_DP
  rhotot(:) = rho(:,1) + rho(:,2)
  !
  CALL init_lda_xc()  !^^^
  !
  WHERE (rhotot <= small) !however the cut is already inside xc_spin()
     rhotot = 0.5_DP
     null_v = 0.0_DP
  ELSEWHERE
     zeta(:) = (rho(:,1) - rho(:,2)) / rhotot(:)
  END WHERE
  !
  WHERE (ABS(zeta) > 1.0_DP)
     zeta   = 0.5_DP
     null_v = 0.0_DP
  END WHERE
  !
  IF (iexch_l == 1 .AND. icorr_l == 1) THEN
     !
     !    first case: analytical derivative available
     !
     ALLOCATE( rs(length) )
     ALLOCATE( vcu(length), vcp(length), vx(length,1) )
     ALLOCATE( ecu(length), ecp(length) )
     !
     ! ... exchange
     !
     rs(:) = ( pi34 / (2.0_DP * rho(:,1)) )**third
     CALL slater( length, rs, aux1, vx(:,1) )
     !
     dmuxc(:,1,1) = vx(:,1) / (3.0_DP * rho(:,1))
     !
     rs(:) = ( pi34 / (2.0_DP * rho(:,2)) )**third
     CALL slater( length, rs, aux1, vx(:,1) )
     !
     dmuxc(:,2,2) = vx(:,1) / (3.0_DP * rho(:,2))
     !
     ! ... correlation
     !
     rs(:) = (pi34 / rhotot(:))**third
     !
     CALL pz( length, rs, 1, ecu, vcu )
     CALL pz_polarized( length, rs, ecp, vcp )
     !
     DO ir = 1, length
        fz  = ( (1.0_DP + zeta(ir))**p43 + (1.0_DP - zeta(ir))**p43 - 2.0_DP ) &
              / (2.0_DP**p43 - 2.0_DP)
        fz1 = p43 * ( (1.0_DP + zeta(ir))**third - (1.0_DP - zeta(ir))**third) &
              / (2.0_DP**p43 - 2.0_DP)
        fz2 = p49 * ( (1.0_DP + zeta(ir))**m23   + (1.0_DP - zeta(ir))**m23  ) &
              / (2.0_DP**p43 - 2.0_DP)
        !
        iflg = 2
        IF (rs(ir) < 1.0_DP) iflg = 1
        !
        dmcu = dpz( rs(ir), iflg )
        dmcp = dpz_polarized( rs(ir), iflg )
        !
        aa = dmcu + fz * (dmcp - dmcu)
        bb = 2.0_DP * fz1 * (vcp(ir) - vcu(ir) - (ecp(ir) - ecu(ir)) ) / rhotot(ir)
        cc = fz2 * (ecp(ir) - ecu(ir)) / rhotot(ir)
        !
        dmuxc(ir,1,1) = dmuxc(ir,1,1) + aa + (1.0_DP - zeta(ir)) * bb + &
                                                       (1.0_DP - zeta(ir))**2 * cc
        dmuxc(ir,2,1) = dmuxc(ir,2,1) + aa + (-zeta(ir)) * bb +         &
                                                       (zeta(ir)**2 - 1.0_DP) * cc
        dmuxc(ir,1,2) = dmuxc(ir,2,1)
        dmuxc(ir,2,2) = dmuxc(ir,2,2) + aa - (1.0_DP + zeta(ir)) * bb + &
                                                       (1.0_DP + zeta(ir))**2 * cc
     ENDDO
     !
     DEALLOCATE( rs )
     DEALLOCATE( vcu, vcp, vx )
     DEALLOCATE( ecu, ecp )
     !
  ELSE
     !
     ALLOCATE( vx(length,2), vc(length,2) )
     ALLOCATE( aux2(length) )
     ALLOCATE( zeta_eff(length), ds(length) )
     !
     ds(:) = MIN( 1.E-6_DP, 1.E-4_DP * rhotot(:) ) !here ds is drho
     !
     CALL xc_spin( length, rhotot+ds, zeta, aux1, aux2, vx, vc )
     !
     dmuxc(:,1,1) = vx(:,1) + vc(:,1)
     dmuxc(:,2,2) = vx(:,2) + vc(:,2)
     !
     CALL xc_spin( length, rhotot-ds, zeta, aux1, aux2, vx, vc )
     !
     dmuxc(:,1,1) = (dmuxc(:,1,1) - vx(:,1) - vc(:,1)) / (2.0_DP * ds(:))
     dmuxc(:,2,1) = dmuxc(:,1,1)
     dmuxc(:,2,2) = (dmuxc(:,2,2) - vx(:,2) - vc(:,2)) / (2.0_DP * ds(:))
     dmuxc(:,1,2) = dmuxc(:,2,2)
     !
     ! ds() = min (1.d-6, 1.d-4 * abs(zeta(:)) )
     ds(:) = 1.E-6_DP  ! now ds is dzeta
     !
     ! If zeta is too close to +-1, the derivative is computed at a slightly
     ! smaller zeta
     !
     zeta_eff(:) = SIGN( MIN( ABS(zeta(:)), (1.0_DP-2.0_DP*ds(:)) ), zeta(:) )
     !
     CALL xc_spin( length, rhotot, zeta_eff+ds, aux1, aux2, vx, vc )
     !
     aux1 = 1.0_DP / rhotot / (2.0_DP*ds)
     DO is = 1, 2
        vx(:,is) = (vx(:,is) + vc(:,is)) * aux1(:)  ! vx as workspace here
     ENDDO
     dmuxc(:,1,1) = dmuxc(:,1,1) + vx(:,1) * (1.0_DP-zeta(:))
     dmuxc(:,1,2) = dmuxc(:,1,2) + vx(:,2) * (1.0_DP-zeta(:))
     dmuxc(:,2,1) = dmuxc(:,2,1) - vx(:,1) * (1.0_DP+zeta(:))
     dmuxc(:,2,2) = dmuxc(:,2,2) - vx(:,2) * (1.0_DP+zeta(:))
     !
     CALL xc_spin( length, rhotot, zeta_eff-ds, aux1, aux2, vx, vc )
     !
     aux1 = 1.0_DP / rhotot / (2.0_DP*ds)
     DO is=1,2
        vx(:,is) = (vx(:,is) + vc(:,is)) * aux1(:)
     ENDDO
     dmuxc(:,1,1) = dmuxc(:,1,1) - vx(:,1) * (1.0_DP-zeta(:))
     dmuxc(:,1,2) = dmuxc(:,1,2) - vx(:,2) * (1.0_DP-zeta(:))
     dmuxc(:,2,1) = dmuxc(:,2,1) + vx(:,1) * (1.0_DP+zeta(:))
     dmuxc(:,2,2) = dmuxc(:,2,2) + vx(:,2) * (1.0_DP+zeta(:))
     !
     DEALLOCATE( vx, vc )
     DEALLOCATE( aux2 )
     DEALLOCATE( zeta_eff, ds )
     !
  ENDIF
  !
  ! bring to rydberg units
  !
  dmuxc(:,1,1) = e2 * dmuxc(:,1,1) * null_v  !up-up
  dmuxc(:,2,1) = e2 * dmuxc(:,2,1) * null_v  !down-up
  dmuxc(:,1,2) = e2 * dmuxc(:,1,2) * null_v  !up-down
  dmuxc(:,2,2) = e2 * dmuxc(:,2,2) * null_v  !down-down
  !
  RETURN
  !
END SUBROUTINE dmxc_spin

!-----------------------------------------------------------------------
SUBROUTINE dmxc_nc( length, rho_in, m, dmuxc )
!-----------------------------------------------------------------------
  !! Computes the derivative of the xc potential with respect to the 
  !! local density and magnetization in the non-collinear case.
  !
  USE xc_lda_lsda,  ONLY: xc_spin
  USE funct,        ONLY: init_lda_xc
  USE kinds,        ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the input/output arrays
  REAL(DP), INTENT(IN), DIMENSION(length) :: rho_in
  !! total charge density
  REAL(DP), INTENT(IN), DIMENSION(length,3) :: m
  !! magnetization vector
  REAL(DP), INTENT(OUT), DIMENSION(length,4,4) :: dmuxc
  !! derivative of XC functional
  !
  ! ... local variables
  !
  REAL(DP), DIMENSION(length)   :: rho, amag, zeta
  REAL(DP), DIMENSION(length)   :: zeta_eff, ds, null_v, null_m
  REAL(DP), DIMENSION(length)   :: vs, aux1, aux2
  REAL(DP), DIMENSION(length,2) :: vx, vxm, vxp
  REAL(DP), DIMENSION(length,2) :: vc, vcm, vcp
  REAL(DP), DIMENSION(length)   :: dvxc_rho, dbx_rho, dby_rho, dbz_rho
  !
  REAL(DP) :: dvxc_mx, dvxc_my, dvxc_mz, &
              dbx_mx, dbx_my, dbx_mz,    &
              dby_mx, dby_my, dby_mz,    &
              dbz_mx, dbz_my, dbz_mz
  REAL(DP) :: rnull
  !
  INTEGER :: i
  !
  REAL(DP), PARAMETER :: small = 1.E-30_DP, e2 = 2.0_DP
  !
  dmuxc = 0.0_DP
  !
  rho    = rho_in
  zeta   = 0.5_DP
  amag   = 0.025_DP
  null_v = 1.0_DP
  null_m = 1.0_DP
  !
  CALL init_lda_xc()
  !
  WHERE (rho_in <= small)
     rho = 0.5_DP
     null_v = 0.0_DP
  ELSEWHERE
     amag = SQRT( m(:,1)**2 + m(:,2)**2 + m(:,3)**2 )
     zeta = amag / rho
  END WHERE
  !
  WHERE (ABS(zeta) > 1.0_DP)
     zeta   = 0.5_DP
     null_v = 0.0_DP
  END WHERE
  !
  CALL xc_spin( length, rho, zeta, aux1, aux2, vx, vc )
  !
  vs = 0.5_DP*( vx(:,1)+vc(:,1)-vx(:,2)-vc(:,2) )
  !
  ! Here ds is drho
  ds = MIN( 1.E-6_DP, 1.E-4_DP * rho )
  !
  CALL xc_spin( length, rho-ds, zeta, aux1, aux2, vxm, vcm )
  !
  CALL xc_spin( length, rho+ds, zeta, aux1, aux2, vxp, vcp )
  !
  dvxc_rho = ((vxp(:,1) + vcp(:,1) - vxm(:,1) - vcm(:,1)) + &
              (vxp(:,2) + vcp(:,2) - vxm(:,2) - vcm(:,2))) / (4.0_DP*ds)
  !
  WHERE (amag < 1.E-10_DP)
     rho  = 0.5_DP
     zeta = 0.5_DP
     amag = 0.025_DP
     null_m = 0.0_DP
  END WHERE
  !
  !
  aux2(:) =  vxp(:,1) + vcp(:,1) - vxm(:,1) - vcm(:,1) - &
           ( vxp(:,2) + vcp(:,2) - vxm(:,2) - vcm(:,2) )
  !
  dbx_rho(:) = aux2 * m(:,1) / (4.0_DP*ds*amag)
  dby_rho(:) = aux2 * m(:,2) / (4.0_DP*ds*amag)
  dbz_rho(:) = aux2 * m(:,3) / (4.0_DP*ds*amag)
  !
  ! Now ds is dzeta
  ! ds = min (1.d-6, 1.d-4 * abs (zeta) )
  ds = 1.0E-6_DP
  !
  ! If zeta is too close to +-1, the derivative is computed at a slightly
  ! smaller zeta
  !
  DO i = 1, length
     zeta_eff(i) = SIGN( MIN( ABS( zeta(i) ), ( 1.0_DP - 2.0_DP*ds(i) ) ) , zeta(i) )
  ENDDO
  !
  CALL xc_spin( length, rho, zeta_eff-ds, aux1, aux2, vxm, vcm )
  !
  CALL xc_spin( length, rho, zeta_eff+ds, aux1, aux2, vxp, vcp )
  !
  !  The variables are rho and m, so zeta depends on rho
  !
  aux1(:) =  vxp(:,1) + vcp(:,1) - vxm(:,1) - vcm(:,1) + &
             vxp(:,2) + vcp(:,2) - vxm(:,2) - vcm(:,2)
  aux2(:) =  vxp(:,1) + vcp(:,1) - vxm(:,1) - vcm(:,1) - &
           ( vxp(:,2) + vcp(:,2) - vxm(:,2) - vcm(:,2) )
  !
  dvxc_rho(:) = dvxc_rho - aux1 * zeta/rho / (4.0_DP*ds) * null_m
  dbx_rho(:)  = dbx_rho  - aux2 * m(:,1) * zeta/rho / (4.0_DP*ds*amag)
  dby_rho(:)  = dby_rho  - aux2 * m(:,2) * zeta/rho / (4.0_DP*ds*amag)
  dbz_rho(:)  = dbz_rho  - aux2 * m(:,3) * zeta/rho / (4.0_DP*ds*amag)
  !
  dmuxc(:,1,1) = dvxc_rho * null_v
  dmuxc(:,2,1) = dbx_rho  * null_v
  dmuxc(:,3,1) = dby_rho  * null_v
  dmuxc(:,4,1) = dbz_rho  * null_v
  !
  ! Here the derivatives with respect to m
  !
  DO i = 1, length
     !
     rnull=null_v(i)
     !
     dvxc_mx = aux1(i) * m(i,1) / rho(i) / (4.0_DP*ds(i)*amag(i))
     dvxc_my = aux1(i) * m(i,2) / rho(i) / (4.0_DP*ds(i)*amag(i))
     dvxc_mz = aux1(i) * m(i,3) / rho(i) / (4.0_DP*ds(i)*amag(i))
     !
     dbx_mx  = (aux2(i) * m(i,1) * m(i,1) * amag(i)/rho(i) / (4.0_DP*ds(i)) + &
                vs(i) * (m(i,2)**2+m(i,3)**2)) / amag(i)**3
     dbx_my  = (aux2(i) * m(i,1) * m(i,2) * amag(i)/rho(i) / (4.0_DP*ds(i)) - &
                vs(i) * m(i,1) * m(i,2) ) / amag(i)**3
     dbx_mz  = (aux2(i) * m(i,1) * m(i,3) * amag(i)/rho(i) / (4.0_DP*ds(i)) - &
                vs(i) * m(i,1) * m(i,3) ) / amag(i)**3
     !
     dby_mx  = dbx_my
     dby_my  = (aux2(i) * m(i,2) * m(i,2) * amag(i)/rho(i) / (4.0_DP*ds(i)) + &
                vs(i) * (m(i,1)**2 + m(i,3)**2)) / amag(i)**3
     dby_mz  = (aux2(i) * m(i,2) * m(i,3) * amag(i)/rho(i) / (4.0_DP*ds(i)) - &
                vs(i) * m(i,2) * m(i,3)) / amag(i)**3
     !
     dbz_mx  = dbx_mz
     dbz_my  = dby_mz
     dbz_mz  = (aux2(i) * m(i,3) * m(i,3) * amag(i)/rho(i) / (4.0_DP*ds(i)) + &
                vs(i)*(m(i,1)**2 + m(i,2)**2)) / amag(i)**3
     !
     dmuxc(i,1,2) = dvxc_mx * rnull
     dmuxc(i,1,3) = dvxc_my * rnull 
     dmuxc(i,1,4) = dvxc_mz * rnull
     !
     dmuxc(i,2,2) = dbx_mx * rnull
     dmuxc(i,2,3) = dbx_my * rnull
     dmuxc(i,2,4) = dbx_mz * rnull
     !
     dmuxc(i,3,2) = dby_mx * rnull
     dmuxc(i,3,3) = dby_my * rnull
     dmuxc(i,3,4) = dby_mz * rnull
     !
     dmuxc(i,4,2) = dbz_mx * rnull
     dmuxc(i,4,3) = dbz_my * rnull
     dmuxc(i,4,4) = dbz_mz * rnull
     !
  ENDDO
  !
  ! bring to rydberg units
  !
  dmuxc = e2 * dmuxc
  !
  RETURN
  !
END SUBROUTINE dmxc_nc
