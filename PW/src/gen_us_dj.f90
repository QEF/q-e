!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE gen_us_dj( ik, dvkb )
  !----------------------------------------------------------------------
  !! Calculates the beta function pseudopotentials with
  !! the derivative of the Bessel functions.
  !
  USE kinds,      ONLY: DP
  USE constants,  ONLY: tpi
  USE ions_base,  ONLY: nat, ntyp => nsp, ityp, tau
  USE cell_base,  ONLY: tpiba, omega
  USE klist,      ONLY: xk, ngk, igk_k
  USE gvect,      ONLY: mill, eigts1, eigts2, eigts3, g
  USE wvfct,      ONLY: npwx
  USE uspp,       ONLY: nkb, indv, nhtol, nhtolm
  USE us,         ONLY: nqx, tab, tab_d2y, dq, spline_ps
  USE m_gth,      ONLY: mk_dffnl_gth
  USE splinelib
  USE uspp_param, ONLY: upf, lmaxkb, nbetam, nh
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik
  !! k-point index
  COMPLEX(DP), INTENT(OUT) :: dvkb(npwx, nkb)
  !! the beta function pseudopotential
  !
  ! ... local variables
  !
  INTEGER :: npw, ikb, nb, ih, ig, i0, i1, i2, i3, nt
  ! counter on beta functions
  ! counter on beta functions
  ! counter on beta functions
  ! counter on G vectors
  ! index of the first nonzero point in the r
  ! counter on atomic type
  !
  REAL(DP) :: arg, px, ux, vx, wx
  ! argument of the atomic phase factor
  !
  COMPLEX(DP) :: phase, pref
  ! atomic phase factor
  ! prefactor
  !
  INTEGER :: na, l, iig, lm, iq
  REAL(DP), ALLOCATABLE :: djl(:,:,:), ylm(:,:), q(:), gk(:,:)
  REAL(DP) :: qt
  COMPLEX(DP), ALLOCATABLE :: sk(:)
  REAL(DP), ALLOCATABLE :: xdata(:)
  !
  IF (nkb == 0) RETURN
  !
  CALL start_clock( 'stres_us31' )
  !
  npw = ngk(ik)
  ALLOCATE( djl(npw,nbetam,ntyp)   )    
  ALLOCATE( ylm(npw,(lmaxkb+1)**2) )    
  ALLOCATE( gk(3,npw) )    
  ALLOCATE( q(npw)    )    
  !
  DO ig = 1, npw
     iig = igk_k(ig,ik)
     gk(1, ig) = xk(1, ik) + g(1, iig)
     gk(2, ig) = xk(2, ik) + g(2, iig)
     gk(3, ig) = xk(3, ik) + g(3, iig)
     q(ig) = gk(1, ig)**2 + gk(2, ig)**2 + gk(3, ig)**2
  ENDDO
  !
  CALL stop_clock( 'stres_us31' )
  CALL start_clock( 'stres_us32' )
  CALL ylmr2( (lmaxkb+1)**2, npw, gk, q, ylm )
  CALL stop_clock( 'stres_us32' )
  CALL start_clock( 'stres_us33' )
  !
  IF (spline_ps) THEN
    ALLOCATE( xdata(nqx) )
    DO iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    ENDDO
  ENDIF
  !
  DO nt = 1, ntyp
     DO nb = 1, upf(nt)%nbeta
        !
        IF ( upf(nt)%is_gth ) THEN
           CALL mk_dffnl_gth( nt, nb, npw, omega, tpiba, q, djl(1,nb,nt) )
           CYCLE
        ENDIF
        DO ig = 1, npw
           qt = SQRT(q (ig)) * tpiba
           IF (spline_ps) THEN
             djl(ig,nb,nt) = splint_deriv(xdata, tab(:,nb,nt), & 
                                                 tab_d2y(:,nb,nt), qt)
           ELSE
             px = qt / dq - INT(qt/dq)
             ux = 1.d0 - px
             vx = 2.d0 - px
             wx = 3.d0 - px
             i0 = qt / dq + 1
             i1 = i0 + 1
             i2 = i0 + 2
             i3 = i0 + 3
             djl(ig,nb,nt) = ( tab(i0, nb, nt) * (-vx*wx-ux*wx-ux*vx)/6.d0 + &
                               tab(i1, nb, nt) * (+vx*wx-px*wx-px*vx)/2.d0 - &
                               tab(i2, nb, nt) * (+ux*wx-px*wx-px*ux)/2.d0 + &
                               tab(i3, nb, nt) * (+ux*vx-px*vx-px*ux)/6.d0 )/dq
           ENDIF
        ENDDO
        !
     ENDDO
  ENDDO
  !
  CALL stop_clock( 'stres_us33' )
  CALL start_clock( 'stres_us34' )
  !
  DEALLOCATE ( q  )
  DEALLOCATE ( gk )
  !
  ALLOCATE ( sk(npw) )
  ikb = 0
  DO nt = 1, ntyp
     DO na = 1, nat
        !
        IF (ityp(na) == nt) THEN
           arg = ( xk(1, ik) * tau(1,na) + &
                   xk(2, ik) * tau(2,na) + &
                   xk(3, ik) * tau(3,na) ) * tpi
           phase = CMPLX( COS(arg), -SIN(arg) ,KIND=DP )
           DO ig = 1, npw
              iig = igk_k(ig,ik)
              sk (ig) = eigts1(mill (1,iig), na) * &
                        eigts2(mill (2,iig), na) * &
                        eigts3(mill (3,iig), na) * phase
           ENDDO
           DO ih = 1, nh(nt)
              nb = indv(ih, nt)
              l = nhtol(ih, nt)
              lm= nhtolm(ih, nt)
              ikb = ikb + 1
              pref = (0.d0, -1.d0) **l
              !
              DO ig = 1, npw
                 dvkb(ig, ikb) = djl(ig, nb, nt) * sk(ig) * ylm(ig, lm) &
                      * pref
              ENDDO
           ENDDO
        ENDIF
        !
     ENDDO
  ENDDO
  !
  CALL stop_clock('stres_us34')
  !
  IF (ikb /= nkb) CALL errore('gen_us_dj', 'unexpected error', 1)
  DEALLOCATE( sk  )
  DEALLOCATE( ylm )
  DEALLOCATE( djl )
  IF (spline_ps) DEALLOCATE( xdata )
  !
  RETURN
  !
END SUBROUTINE gen_us_dj

