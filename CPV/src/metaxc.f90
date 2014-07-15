
!
! Copyright (C) 2005-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE tpssmeta(nnr, nspin,grho,rho,kedtau,etxc)
  !     ===================
  !--------------------------------------------------------------------
  use kinds, only: dp
  use funct, only: tau_xc, tau_xc_spin
  IMPLICIT NONE
  !
  ! input
  integer nspin , nnr
  real(dp)  grho(nnr,3,nspin), rho(nnr,nspin),kedtau(nnr,nspin)
  ! output: excrho: exc * rho ;  E_xc = \int excrho(r) d_r
  ! output: rhor:   contains the exchange-correlation potential
  real(dp)  etxc
  REAL(dp) :: zeta, rh, grh2
  INTEGER :: k, ipol, is
  REAL(dp) :: grho2 (2), sx, sc, v1x, v2x, v3x,v1c, v2c, v3c, &
       v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw ,v2cup(3),v2cdw(3), &
       v3xup, v3xdw,grhoup(3),grhodw(3),v3cup, v3cdw, segno, arho, atau
  REAL(dp), PARAMETER :: epsr = 1.0d-6, epsg = 1.0d-10
  etxc = 0.d0
  ! calculate the gradient of rho+rho_core in real space
  DO k = 1, nnr
     DO is = 1, nspin
        grho2 (is) = grho(k,1, is)**2 + grho(k,2,is)**2 + grho(k,3, is)**2
     ENDDO
     IF (nspin == 1) THEN
        !
        !    This is the spin-unpolarised case
        !
        arho = ABS (rho (k, 1) )
        segno = SIGN (1.d0, rho (k, 1) )
        atau = kedtau(k,1)
        IF (arho.GT.epsr.AND.grho2 (1) .GT.epsg.AND.ABS(atau).GT.epsr) THEN
           CALL tau_xc (arho, grho2(1), atau, sx, sc, &
                v1x, v2x, v3x, v1c, v2c, v3c)
           rho (k, 1) =  (v1x + v1c )
           kedtau(k,1)=  (v3x + v3c) *0.5d0
           ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
           DO ipol = 1, 3  
              grho(k,ipol,1) =  (v2x + v2c)*grho (k,ipol,1) 
           ENDDO
           etxc = etxc +  (sx + sc) * segno 
        ELSE  
           DO ipol = 1, 3  
              grho (k, ipol, 1) = 0.d0  
           ENDDO
           kedtau(k,1)=0.d0
        ENDIF
     ELSE
        !
        !    spin-polarised case
        !
        !CALL tpsscx_spin(rho (k, 1), rho (k, 2), grho2 (1), grho2 (2), &
        !     kedtau(k,1),kedtau(k,2),sx, &
        !     v1xup,v1xdw,v2xup,v2xdw,v3xup,v3xdw)
        rh = rho (k, 1) + rho (k, 2)
        IF (rh.GT.epsr) THEN
           !zeta = (rho (k, 1) - rho (k, 2) ) / rh
           DO ipol=1,3
              grhoup(ipol)=grho(k,ipol,1)
              grhodw(ipol)=grho(k,ipol,2)
           END DO
           ! atau=kedtau(k,1)+kedtau(k,2)
           call tau_xc_spin (rho(k,1), rho(k,2), grhoup, grhodw, &
                kedtau(k,1), kedtau(k,2), sx, sc, v1xup, v1xdw, v2xup, v2xdw, &
                v3xup, v3xdw, v1cup, v1cdw, v2cup, v2cdw,&
                v3cup, v3cdw)
           !CALL tpsscc_spin(rh,zeta,grhoup,grhodw, &
           !     atau,sc,v1cup,v1cdw,v2cup,v2cdw,v3c)
        ELSE
           sx = 0.d0  
           sc = 0.d0  
           v1xup = 0.d0  
           v1xdw = 0.d0  
           v2xup=0.d0
           v2xdw=0.d0
           v3xup=0.d0
           v3xdw=0.d0
           v1cup = 0.d0  
           v1cdw = 0.d0  
           v2cup=0.d0
           v2cdw=0.d0
           v3cup=0.d0
           v3cdw=0.d0
           !
        ENDIF
        !
        ! first term of the gradient correction : D(rho*Exc)/D(rho)
        !
        rho(k, 1) =  (v1xup + v1cup)
        rho(k, 2) =  (v1xdw + v1cdw) 
        !
        ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
        !
        DO ipol = 1, 3  
           grho(k,ipol,1) = (v2xup*grho(k,ipol,1) + v2cup(ipol))
           grho(k,ipol,2) = (v2xdw*grho(k,ipol,2) + v2cdw(ipol)) 
        ENDDO
        kedtau(k,1)=  (v3xup + v3cup) *0.5d0
        kedtau(k,2)=  (v3xdw + v3cdw) *0.5d0
        etxc = etxc +  (sx + sc)
     ENDIF
  ENDDO
  RETURN
END SUBROUTINE tpssmeta

!-----------------------------------------------------------------------
