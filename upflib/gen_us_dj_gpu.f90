!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE gen_us_dj_gpu_ &
     ( npw, npwx, igk_d, xk, nat, tau, ityp, ntyp, tpiba, &
       omega, nr1, nr2, nr3, eigts1_d, eigts2_d, eigts3_d, mill_d, g_d, dvkb_d )
  !----------------------------------------------------------------------
  !! Calculates the beta function pseudopotentials with
  !! the derivative of the Bessel functions - GFPU version
  !
  ! AF: more gpu-resident variables can be used, avoiding local GPU-alloc
  !     and host2dev transfers
  !
  USE upf_kinds,   ONLY: dp
  USE upf_const,   ONLY: tpi
  USE uspp,        ONLY: nkb, indv_d, nhtol_d, nhtolm_d
  USE uspp_data,   ONLY: nqx, tab, tab_d, dq
  USE uspp_param,  ONLY: upf, lmaxkb, nbetam, nh, nhm
  USE device_fbuff_m,   ONLY: dev_buf
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw
  !! number ok plane waves 
  INTEGER, INTENT(IN) :: npwx
  !! max number ok plane waves across k-points
  INTEGER, INTENT(IN) :: igk_d(npw)
  !! indices of plane waves k+G
  REAL(dp), INTENT(IN) :: xk(3)
  !! k-point
  INTEGER, INTENT(IN) :: nat
  !! number of atoms
  INTEGER, INTENT(IN) :: ityp(nat)
  !! index of type per atom
  INTEGER, INTENT(IN) :: ntyp
  !! number of atomic types
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! atomic positions (cc alat units)
  REAL(DP), INTENT(IN) :: tpiba
  !! rec.lattice units 2pi/a
  REAL(DP), INTENT(IN) :: omega
  !! cell volume
  INTEGER, INTENT(IN) :: nr1,nr2,nr3
  !! fft dims (dense grid)
  COMPLEX(DP), INTENT(IN) :: eigts1_d(-nr1:nr1,nat)
  !! structure factor 1
  COMPLEX(DP), INTENT(IN) :: eigts2_d(-nr2:nr2,nat)
  !! structure factor 2
  COMPLEX(DP), INTENT(IN) :: eigts3_d(-nr3:nr3,nat)
  !! structure factor 3
  INTEGER, INTENT(IN) :: mill_d(3,*)
  !! miller index map
  REAL(DP), INTENT(IN) :: g_d(3,*)
  !! g vectors (2pi/a units)
  COMPLEX(DP), INTENT(OUT) :: dvkb_d(npwx, nkb)
  !! the beta function pseudopotential
  !
  ! ... local variables
  !
  INTEGER  :: na, nt, nb, ih, l, lm, ikb, iig, i0, i1, i2, &
              i3, ig, nbm, iq, mil1, mil2, mil3,      &
              ikb_t, nht, ina, nas(nat), ierr(3)
  REAL(DP) :: px, ux, vx, wx, arg, u_ipol, xk1, xk2, xk3, qt
  COMPLEX(DP) :: pref
  INTEGER,  ALLOCATABLE :: ityp_d(:), ih_d(:), na_d(:), nas_d(:)
  !
  REAL(DP), POINTER :: gk_d(:,:), djl_d(:,:,:),  ylm_d(:,:)
  REAL(DP), ALLOCATABLE :: q_d(:), tau_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: phase_d(:), sk_d(:,:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: igk_d, mill_d, eigts1_d, eigts2_d, eigts3_d, g_d
  attributes(DEVICE) :: gk_d, q_d, sk_d, djl_d, ylm_d, &
                        ityp_d, phase_d, ih_d, na_d, tau_d, nas_d
  attributes(DEVICE) :: dvkb_d
  !
  IF (nkb == 0) RETURN
  !
  CALL dev_buf%lock_buffer( ylm_d, (/ npw,(lmaxkb+1)**2 /), ierr(1) )
  CALL dev_buf%lock_buffer( djl_d, (/ npw,nbetam,ntyp /), ierr(2) )
  CALL dev_buf%lock_buffer( gk_d,  (/ 3,npw /), ierr(3) )
  ALLOCATE( q_d(npw) )
  !
  xk1 = xk(1)
  xk2 = xk(2)
  xk3 = xk(3)
  !
  !$cuf kernel do (1) <<<*,*>>>
  DO ig = 1, npw
     iig = igk_d(ig)
     gk_d(1,ig) = xk1 + g_d(1,iig)
     gk_d(2,ig) = xk2 + g_d(2,iig)
     gk_d(3,ig) = xk3 + g_d(3,iig)
     q_d(ig) = gk_d(1,ig)**2 +  gk_d(2,ig)**2 + gk_d(3,ig)**2
  ENDDO
  !
  CALL ylmr2_gpu( (lmaxkb+1)**2, npw, gk_d, q_d, ylm_d )
  !
  DO nt = 1, ntyp
      nbm = upf(nt)%nbeta
      !$cuf kernel do (2) <<<*,*>>>
      DO nb = 1, nbm
         DO ig = 1, npw
            qt = SQRT(q_d(ig)) * tpiba
            px = qt/dq - DBLE(INT(qt/dq))
            ux = 1._DP - px
            vx = 2._DP - px
            wx = 3._DP - px
            i0 = INT(qt/dq) + 1
            i1 = i0 + 1
            i2 = i0 + 2
            i3 = i0 + 3
            djl_d(ig,nb,nt) = (tab_d(i0,nb,nt) * (-vx*wx-ux*wx-ux*vx)/6._DP + &
                               tab_d(i1,nb,nt) * (+vx*wx-px*wx-px*vx)/2._DP - &
                               tab_d(i2,nb,nt) * (+ux*wx-px*wx-px*ux)/2._DP + &
                               tab_d(i3,nb,nt) * (+ux*vx-px*vx-px*ux)/6._DP)/dq
         ENDDO
      ENDDO
  ENDDO
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
      iig = igk_d(ig)
      mil1 = mill_d(1,iig)
      mil2 = mill_d(2,iig)
      mil3 = mill_d(3,iig)
      sk_d(ig,na) = eigts1_d(mil1,na) * &
                    eigts2_d(mil2,na) * &
                    eigts3_d(mil3,na) * phase_d(na)
    ENDDO
  ENDDO
  !
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
  !$cuf kernel do (2) <<<*,*>>>
  DO ikb = 1, ikb_t
    DO ig = 1, npw
      ih = ih_d(ikb)
      na = na_d(ikb)
      nt = ityp_d(na)
      nb = indv_d(ih,nt)
      l  = nhtol_d(ih,nt)
      lm = nhtolm_d(ih,nt)
      pref = (0._DP,-1._DP)**l
      !
      dvkb_d(ig,ikb) = CMPLX(djl_d(ig,nb,nt)) * sk_d(ig,na) * &
                       CMPLX(ylm_d(ig,lm))  * pref
    ENDDO
  ENDDO  
  !
  DEALLOCATE( sk_d )
  !
  IF (ikb_t /= nkb) CALL upf_error( 'gen_us_dj', 'unexpected error', 1 )
  !
  CALL dev_buf%release_buffer( ylm_d, ierr(1) )
  CALL dev_buf%release_buffer( djl_d, ierr(2) )
  CALL dev_buf%release_buffer( gk_d, ierr(3) )
  !
  DEALLOCATE( ih_d, na_d, nas_d )
#endif
  !
  RETURN
  !
END SUBROUTINE gen_us_dj_gpu_
