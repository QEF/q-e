!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE gen_us_dy_gpu( ik, u, dvkb_d )
  !----------------------------------------------------------------------
  !! Calculates the kleinman-bylander pseudopotentials with the
  !! derivative of the spherical harmonics projected on vector u
  !
  USE kinds,       ONLY: DP
  USE io_global,   ONLY: stdout
  USE constants,   ONLY: tpi
  USE ions_base,   ONLY: nat, ntyp => nsp, ityp, tau
  USE cell_base,   ONLY: tpiba
  USE klist,       ONLY: xk, ngk, igk_k_d
  USE wvfct,       ONLY: npwx
  USE uspp,        ONLY: nkb, indv, nhtol, nhtolm
  USE us,          ONLY: nqx, tab, tab_d2y, dq, spline_ps
  USE splinelib
  USE uspp_param,  ONLY: upf, lmaxkb, nbetam, nh, nhm
  !
  USE us_gpum,     ONLY: using_tab, using_tab_d2y, &
                         using_tab_d, tab_d
  USE gvect_gpum,  ONLY: mill_d, eigts1_d, eigts2_d, eigts3_d, g_d
  USE device_fbuff_m,    ONLY: dev_buf
  !
  IMPLICIT NONE
  !
  INTEGER  :: ik
  !! input: k-point index
  REAL(DP) :: u(3)
  !! input: projection vector
  COMPLEX(DP) :: dvkb_d(npwx,nkb)
  !! output: kleinman-bylander pseudopotential
  !
  ! ... local variables
  !
  INTEGER :: na, nt, nb, ih, l, lm, ikb, iig, ipol, i0, i1, i2, &
             i3, ig, npw, nbm, iq, mil1, mil2, mil3, ikb_t,     &
             nht, ina, lmx2
  INTEGER :: nas(nat), ierr(4)
  !
  INTEGER, ALLOCATABLE :: ityp_d(:), ih_d(:), na_d(:),         &
                          nas_d(:), indv_d(:,:), nhtol_d(:,:), &
                          nhtolm_d(:,:)
  !
  REAL(DP), ALLOCATABLE :: q(:), vkb0(:,:,:), dylm(:,:)
  REAL(DP), ALLOCATABLE :: xdata(:), tau_d(:,:), q_d(:)
  !
  REAL(DP), POINTER :: gk_d(:,:)
  REAL(DP), POINTER :: vkb0_d(:,:,:), dylm_u_d(:,:), dylm_d(:,:,:)
  ! dylm = d Y_lm/dr_i in cartesian axes
  ! dylm_u as above projected on u
  COMPLEX(DP), ALLOCATABLE :: phase_d(:), sk_d(:,:)
  !
  REAL(DP) :: px, ux, vx, wx, arg, u_ipol1, u_ipol2, u_ipol3, xk1, xk2, xk3
  COMPLEX(DP) :: pref
  !
#if defined(__CUDA)
  attributes(DEVICE) :: dvkb_d, gk_d, q_d, sk_d, vkb0_d, &
                        dylm_u_d, dylm_d, indv_d, nhtol_d, nhtolm_d, &
                        ityp_d, phase_d, ih_d, na_d, tau_d, nas_d
#endif
  !
  dvkb_d = (0._DP,0._DP)
  !
  IF (lmaxkb <= 0) RETURN
  !
  CALL using_tab(0)
  CALL using_tab_d(0)
  !
  IF (spline_ps) CALL using_tab_d2y(0)
  !
  npw = ngk(ik)
  lmx2 = (lmaxkb+1)**2
  !
  CALL dev_buf%lock_buffer( dylm_u_d, (/ npw,lmx2 /), ierr(1) )
  CALL dev_buf%lock_buffer( vkb0_d, (/ npw,nbetam,ntyp /), ierr(2) )
  CALL dev_buf%lock_buffer( gk_d, (/ 3,npw /), ierr(3) )
  IF (ANY(ierr /= 0)) CALL errore( 'gen_us_dy_gpu', 'cannot allocate buffers', -1 )
  ALLOCATE( q_d(npw) )
  !
  xk1 = xk(1,ik)
  xk2 = xk(2,ik)
  xk3 = xk(3,ik)
  !
  !$cuf kernel do <<<*,*>>>
  DO ig = 1, npw
     iig = igk_k_d(ig,ik)
     gk_d(1,ig) = xk1 + g_d(1,iig)
     gk_d(2,ig) = xk2 + g_d(2,iig)
     gk_d(3,ig) = xk3 + g_d(3,iig)
     q_d(ig) = gk_d(1,ig)**2 +  gk_d(2,ig)**2 + gk_d(3,ig)**2
  ENDDO
  !
  CALL dev_buf%lock_buffer( dylm_d, (/npw,lmx2,3/), ierr(4) )
  DO ipol = 1, 3
     CALL dylmr2_gpu( lmx2, npw, gk_d, q_d, dylm_d(:,:,ipol), ipol )
  ENDDO   
  !
  u_ipol1 = u(1) ; u_ipol2 = u(2) ; u_ipol3 = u(3)
  !
  !$cuf kernel do (2) <<<*,*>>>
  DO lm = 1, lmx2
    DO ig = 1, npw
      dylm_u_d(ig,lm) = u_ipol1*dylm_d(ig,lm,1) + &
                        u_ipol2*dylm_d(ig,lm,2) + &
                        u_ipol3*dylm_d(ig,lm,3)
    ENDDO
  ENDDO
  CALL dev_buf%release_buffer( dylm_d, ierr(4) )
  !
  !
  !$cuf kernel do (1) <<<*,*>>>
  DO ig = 1, npw
     q_d(ig) = SQRT(q_d(ig)) * tpiba
  ENDDO
  !
  !
  IF ( spline_ps ) THEN
    !
    ALLOCATE( q(npw), xdata(nqx), vkb0(npw,nbetam,ntyp) )
    q = q_d
    DO iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    ENDDO
    !
    DO nt = 1, ntyp
      ! calculate beta in G-space using an interpolation table
      DO nb = 1, upf(nt)%nbeta
        DO ig = 1, npw
           vkb0(ig,nb,nt) = splint( xdata, tab(:,nb,nt), &
                                    tab_d2y(:,nb,nt), q(ig) )
        ENDDO
      ENDDO
    ENDDO
    vkb0_d = vkb0
    !
    DEALLOCATE( q, xdata, vkb0 )
    !
  ELSE
    !
    DO nt = 1, ntyp
      nbm = upf(nt)%nbeta
      !$cuf kernel do (2) <<<*,*>>>
      DO nb = 1, nbm
        DO ig = 1, npw
           px = q_d(ig)/dq - DBLE(INT(q_d(ig)/dq))
           ux = 1._DP - px
           vx = 2._DP - px
           wx = 3._DP - px
           i0 = INT(q_d(ig)/dq) + 1
           i1 = i0 + 1
           i2 = i0 + 2
           i3 = i0 + 3
           vkb0_d(ig,nb,nt) = tab_d(i0,nb,nt) * ux * vx * wx / 6._DP + &
                              tab_d(i1,nb,nt) * px * vx * wx / 2._DP - &
                              tab_d(i2,nb,nt) * px * ux * wx / 2._DP + &
                              tab_d(i3,nb,nt) * px * ux * vx / 6._DP
        ENDDO
      ENDDO
    ENDDO
    !
  ENDIF
  !
  DEALLOCATE( q_d )
  !
  ALLOCATE( ityp_d(nat), nas_d(nat) )
  ityp_d = ityp
  ALLOCATE( tau_d(3,nat) )
  tau_d  = tau
  ALLOCATE( phase_d(nat) )
  !
  ina = 0
  DO nt = 1, ntyp
  DO na = 1, nat
    IF ( ityp(na) == nt ) THEN
      ina = ina + 1
      nas(ina) = na
    ENDIF
  ENDDO
  ENDDO
  nas_d = nas
  !
  !$cuf kernel do (1) <<<*,*>>>
  DO ina = 1, nat
     na = nas_d(ina)
     arg = (xk1 * tau_d(1,na) + xk2 * tau_d(2,na) &
          + xk3 * tau_d(3,na) ) * tpi
     phase_d(na) = CMPLX( COS(arg), -SIN(arg), KIND=DP )
  ENDDO
  !
  DEALLOCATE( tau_d )
  !
  ALLOCATE( sk_d(npw,nat) )
  !
  !$cuf kernel do (2) <<<*,*>>>
  DO ina = 1, nat
    DO ig = 1, npw
      !
      na = nas_d(ina)
      iig = igk_k_d(ig,ik)
      mil1 = mill_d(1,iig)
      mil2 = mill_d(2,iig)
      mil3 = mill_d(3,iig)
      sk_d(ig,na) = eigts1_d(mil1,na) * &
                    eigts2_d(mil2,na) * &
                    eigts3_d(mil3,na) * phase_d(na)
    ENDDO
  ENDDO
  !
  DEALLOCATE( phase_d )
  !
  ALLOCATE( ih_d(nat*nhm), na_d(nat*nhm) )
  !
  ikb_t = 0
  DO ina = 1, nat
    na = nas(ina)
    nht = nh(ityp(na))
    !$cuf kernel do (1) <<<*,*>>>
    DO ih = 1, nht
       ih_d(ikb_t+ih) = ih
       na_d(ikb_t+ih) = na
    ENDDO
    ikb_t = ikb_t + nht
  ENDDO
  !
  ALLOCATE( indv_d(nhm,ntyp), nhtol_d(nhm,ntyp), nhtolm_d(nhm,ntyp) )
  !
  indv_d   = indv
  nhtol_d  = nhtol
  nhtolm_d = nhtolm
  !
  !$cuf kernel do (2) <<<*,*>>>
  DO ikb = 1, ikb_t
    DO ig = 1, npw
      ih = ih_d(ikb)
      na = na_d(ikb)
      nt = ityp_d(na)
      nb = indv_d(ih,nt)
      l  = nhtol_d(ih,nt)
      lm = nhtolm_d(ih,nt)
      pref = (0._DP, -1._DP)**l
      !
      dvkb_d(ig,ikb) = CMPLX(vkb0_d(ig,nb,nt)) * sk_d(ig,na) * &
                       CMPLX(dylm_u_d(ig,lm))  * pref / CMPLX(tpiba) 
    ENDDO
  ENDDO  
  !
  DEALLOCATE( sk_d )
  !
  IF (ikb_t /= nkb) THEN
     WRITE( stdout, * ) ikb_t, nkb
     CALL errore( 'gen_us_dy', 'unexpected error', 1 )
  ENDIF
  !
  CALL dev_buf%release_buffer( dylm_u_d, ierr(1) )
  CALL dev_buf%release_buffer( vkb0_d, ierr(2) )
  CALL dev_buf%release_buffer( gk_d, ierr(3) )
  !
  DEALLOCATE( ih_d, na_d, nas_d )
  DEALLOCATE( indv_d, nhtol_d, nhtolm_d )
  !
  !
  RETURN
  !
END SUBROUTINE gen_us_dy_gpu
!
