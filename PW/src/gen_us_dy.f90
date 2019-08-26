!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE gen_us_dy( ik, u, dvkb )
  !----------------------------------------------------------------------
  !! Calculates the kleinman-bylander pseudopotentials with the
  !! derivative of the spherical harmonics projected on vector u
  !
  USE kinds,       ONLY: DP
  USE io_global,   ONLY: stdout
  USE constants,   ONLY: tpi
  USE ions_base,   ONLY: nat, ntyp => nsp, ityp, tau
  USE cell_base,   ONLY: tpiba
  USE klist,       ONLY: xk, ngk, igk_k
  USE gvect,       ONLY: mill, eigts1, eigts2, eigts3, g
  USE wvfct,       ONLY: npwx
  USE uspp,        ONLY: nkb, indv, nhtol, nhtolm
  USE us,          ONLY: nqx, tab, tab_d2y, dq, spline_ps
  USE splinelib
  USE uspp_param,  ONLY: upf, lmaxkb, nbetam, nh
  !
  IMPLICIT NONE
  !
  INTEGER  :: ik
  !! input: k-point index
  REAL(DP) :: u(3)
  !! input: projection vector
  COMPLEX(DP) :: dvkb(npwx,nkb)
  !! output: kleinman-bylander pseudopotential
  !
  ! ... local variables
  !
  INTEGER :: na, nt, nb, ih, l, lm, ikb, iig, ipol, i0, i1, i2, &
             i3, ig, npw
  REAL(DP), ALLOCATABLE :: gk(:,:), q(:)
  REAL(DP) :: px, ux, vx, wx, arg
  !
  REAL(DP), ALLOCATABLE :: vkb0(:,:,:), dylm(:,:), dylm_u(:,:)
  ! dylm = d Y_lm/dr_i in cartesian axes
  ! dylm_u as above projected on u
  !
  COMPLEX(DP), ALLOCATABLE :: sk(:)
  COMPLEX(DP) :: phase, pref
  !
  INTEGER :: iq
  REAL(DP), ALLOCATABLE :: xdata(:)
  !
  dvkb(:,:) = (0.d0, 0.d0)
  IF (lmaxkb <= 0) RETURN
  !
  npw = ngk(ik)
  ALLOCATE( vkb0(npw,nbetam,ntyp), dylm_u(npw,(lmaxkb+1)**2), gk(3,npw) )
  ALLOCATE( q(npw) )
  !
  DO ig = 1, npw
     iig = igk_k(ig,ik)
     gk(1, ig) = xk(1, ik) + g(1, iig)
     gk(2, ig) = xk(2, ik) + g(2, iig)
     gk(3, ig) = xk(3, ik) + g(3, iig)
     q(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  ENDDO
  !
  ALLOCATE( dylm(npw,(lmaxkb+1)**2) )
  dylm_u(:,:) = 0.d0
  DO ipol = 1, 3
     CALL dylmr2( (lmaxkb+1)**2, npw, gk, q, dylm, ipol )
     CALL daxpy( npw * (lmaxkb + 1) **2, u(ipol), dylm, 1, dylm_u, 1 )
  ENDDO
  DEALLOCATE( dylm )
  !
  DO ig = 1, npw
     q(ig) = SQRT(q(ig)) * tpiba
  ENDDO
  !
  IF ( spline_ps ) THEN
    ALLOCATE( xdata(nqx) )
    DO iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    ENDDO
  ENDIF
  !
  DO nt = 1, ntyp
     ! calculate beta in G-space using an interpolation table
     DO nb = 1, upf(nt)%nbeta
        DO ig = 1, npw
           IF ( spline_ps ) THEN
             vkb0(ig,nb,nt) = splint( xdata, tab(:,nb,nt), &
                                      tab_d2y(:,nb,nt), q(ig) )
           ELSE
             px = q(ig)/dq - INT(q(ig)/dq)
             ux = 1.d0 - px
             vx = 2.d0 - px
             wx = 3.d0 - px
             i0 = q(ig)/dq + 1
             i1 = i0 + 1
             i2 = i0 + 2
             i3 = i0 + 3
             vkb0(ig, nb, nt) = tab(i0, nb, nt) * ux * vx * wx / 6.d0 + &
                                tab(i1, nb, nt) * px * vx * wx / 2.d0 - &
                                tab(i2, nb, nt) * px * ux * wx / 2.d0 + &
                                tab(i3, nb, nt) * px * ux * vx / 6.d0
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  DEALLOCATE( q )
  !
  ALLOCATE( sk(npw) )
  !
  ikb = 0
  DO nt = 1, ntyp
     DO na = 1, nat
        IF ( ityp(na) == nt ) THEN
           arg = (xk(1, ik) * tau(1, na) + xk(2, ik) * tau(2, na) &
                + xk(3, ik) * tau(3, na) ) * tpi
           phase = CMPLX( COS(arg), -SIN(arg), KIND=DP )
           DO ig = 1, npw
              iig = igk_k(ig,ik)
              sk(ig) = eigts1(mill (1,iig), na) * &
                       eigts2(mill (2,iig), na) * &
                       eigts3(mill (3,iig), na) * phase
           ENDDO
           !
           DO ih = 1, nh(nt)
              nb = indv(ih, nt)
              l = nhtol(ih, nt)
              lm = nhtolm(ih, nt)
              ikb = ikb + 1
              pref = (0.d0, -1.d0)**l
              !
              DO ig = 1, npw
                 dvkb(ig, ikb) = vkb0(ig, nb, nt) * sk(ig) * dylm_u(ig, lm) &
                      * pref / tpiba
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
  IF (ikb /= nkb) THEN
     WRITE( stdout, * ) ikb, nkb
     CALL errore( 'gen_us_dy', 'unexpected error', 1 )
  ENDIF
  !
  DEALLOCATE( sk )
  DEALLOCATE( vkb0, dylm_u, gk )
  IF (spline_ps) DEALLOCATE( xdata )
  !
  RETURN
  !
END SUBROUTINE gen_us_dy
!
