!
! Copyright (C) 2001-2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE gen_us_dj_base( npw, npwx, igk, xk, nat, tau, ityp, ntyp, tpiba, &
                           omega, nr1, nr2, nr3, eigts1, eigts2, eigts3,    &
                           mill, g, dvkb )
  !----------------------------------------------------------------------
  !! Calculates the beta function pseudopotentials with
  !! the derivative of the Bessel functions.
  !
  USE upf_kinds,  ONLY: dp
  USE upf_const,  ONLY: tpi
  USE uspp,       ONLY: nkb, indv, nhtol, nhtolm
  USE uspp_param, ONLY: lmaxkb, nbetam, nh, nhm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw
  !! number ok plane waves
  INTEGER, INTENT(IN) :: npwx
  !! max number ok plane waves across k-points
  INTEGER, INTENT(IN) :: igk(npw)
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
  INTEGER, INTENT(IN) :: nr1, nr2, nr3
  !! fft dims (dense grid)
  COMPLEX(DP), INTENT(IN) :: eigts1(-nr1:nr1,nat)
  !! structure factor 1
  COMPLEX(DP), INTENT(IN) :: eigts2(-nr2:nr2,nat)
  !! structure factor 2
  COMPLEX(DP), INTENT(IN) :: eigts3(-nr3:nr3,nat)
  !! structure factor 3
  INTEGER, INTENT(IN) :: mill(3,*)
  !! miller index map
  REAL(DP), INTENT(IN) :: g(3,*)
  !! g vectors (2pi/a units)
  COMPLEX(DP), INTENT(OUT) :: dvkb(npwx,nkb)
  !! the beta function pseudopotential
  !
  ! ... local variables
  !
  INTEGER :: ikb, nb, ih, ig, nt
  ! counter on beta functions
  ! counter on beta functions
  ! counter on beta functions
  ! counter on G vectors
  ! counter on atomic type
  !
  INTEGER :: ina, na, l, iig, lm, ikb_t, nht
  REAL(DP) :: arg
  ! argument of the atomic phase factor
  COMPLEX(DP) :: pref
  ! prefactor
  !
  INTEGER,     ALLOCATABLE :: nas(:), ihv(:), nav(:)
  REAL(DP),    ALLOCATABLE :: djl(:,:,:), ylm(:,:), q(:), gk(:,:)
  COMPLEX(DP), ALLOCATABLE :: sk(:,:), phase(:)
  ! phase: atomic phase factor
  !
  IF (nkb == 0) RETURN
  !
  !$acc data present_or_copyin(igk,eigts1,eigts2,eigts3,mill,g) &
  !$acc      present_or_copyout(dvkb) copyin(ityp,tau,xk)
  !
  !CALL start_clock( 'stres_us31' )
  !
  ALLOCATE( djl(npw,nbetam,ntyp), ylm(npw,(lmaxkb+1)**2) )
  ALLOCATE( gk(3,npw), q(npw) )
  !$acc data create( djl, ylm )
  !$acc data create( gk, q )
  !
  !$acc parallel loop
  DO ig = 1, npw
     iig = igk(ig)
     gk(1,ig) = xk(1) + g(1,iig)
     gk(2,ig) = xk(2) + g(2,iig)
     gk(3,ig) = xk(3) + g(3,iig)
     q(ig) = gk(1,ig)**2 + gk(2,ig)**2 + gk(3,ig)**2
  ENDDO
  !
#if defined(__CUDA)
  !$acc host_data use_device(gk,q,ylm)
  CALL ylmr2_gpu( (lmaxkb+1)**2, npw, gk, q, ylm )
  !$acc end host_data
#else
  !$acc update self(gk,q)
  CALL ylmr2( (lmaxkb+1)**2, npw, gk, q, ylm )
  !$acc update device(ylm)
#endif
  !
  !$acc parallel loop
  DO ig = 1, npw
     q(ig) = SQRT(q(ig)) * tpiba
  ENDDO
  DO nt = 1, ntyp
     CALL interp_dbeta( nt, npw, q, djl(:,:,nt) )
  ENDDO
  !
  !CALL stop_clock( 'stres_us33' )
  !CALL start_clock( 'stres_us34' )
  !
  !$acc end data
  DEALLOCATE( q, gk )
  !
  !
  ALLOCATE( phase(nat), sk(npw,nat) )
  ALLOCATE( nas(nat), ihv(nat*nhm), nav(nat*nhm) )
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
  !
  !$acc data create( phase, sk, ihv, nav ) copyin( nas )
  !
  !$acc parallel loop
  DO ina = 1, nat
     na = nas(ina)
     arg = (xk(1) * tau(1,na) + xk(2) * tau(2,na) &
          + xk(3) * tau(3,na) ) * tpi
     phase(na) = CMPLX( COS(arg), -SIN(arg), KIND=DP )
  ENDDO
  !
  !$acc parallel loop collapse(2)
  DO ina = 1, nat
    DO ig = 1, npw
      !
      na = nas(ina)
      iig = igk(ig)
      sk(ig,na) = eigts1(mill(1,iig),na) * &
                  eigts2(mill(2,iig),na) * &
                  eigts3(mill(3,iig),na) * phase(na)
    ENDDO
  ENDDO
  !
  ikb_t = 0
  DO ina = 1, nat
    na = nas(ina)
    nht = nh(ityp(na))
    !$acc kernels
    DO ih = 1, nht
       ihv(ikb_t+ih) = ih
       nav(ikb_t+ih) = na
    ENDDO
    !$acc end kernels
    ikb_t = ikb_t + nht
  ENDDO
  !
  !
  !$acc parallel loop collapse(2) copyin(indv,nhtol,nhtolm)
  DO ikb = 1, ikb_t
    DO ig = 1, npw
      ih = ihv(ikb)
      na = nav(ikb)
      nt = ityp(na)
      nb = indv(ih,nt)
      l  = nhtol(ih,nt)
      lm = nhtolm(ih,nt)
      pref = (0._DP,-1._DP)**l
      !
      dvkb(ig,ikb) = CMPLX(djl(ig,nb,nt),KIND=DP) * sk(ig,na) * &
                     CMPLX(ylm(ig,lm),KIND=DP) * pref
    ENDDO
  ENDDO
  !
  !ALL stop_clock('stres_us34')
  !
  !$acc end data
  !$acc end data
  DEALLOCATE( phase, sk )
  DEALLOCATE( nas, ihv, nav )
  !
  IF (ikb_t /= nkb) CALL upf_error('gen_us_dj', 'unexpected error', 1)
  !
  !$acc end data
  DEALLOCATE( djl, ylm )
  !
  RETURN
  !
END SUBROUTINE gen_us_dj_base
!
!----------------------------------------------------------------------
SUBROUTINE interp_dbeta( nt, npw, qg, vq )
  !----------------------------------------------------------------------
  !
  USE upf_kinds,  ONLY : dp
  USE uspp_param, ONLY : upf, nbetam
  USE uspp_data,  ONLY : nqx, dq, tab_beta
  !
  implicit none
  integer, intent(in) :: nt, npw
  real(dp), intent(in ) :: qg(npw)
  real(dp), intent(out) :: vq(npw,nbetam)
  !
  integer :: i0, i1, i2, i3, nbnt, nb, ig
  real(dp):: qgr, px, ux, vx, wx
  !
  nbnt = upf(nt)%nbeta
  !$acc data present (tab_beta) present_or_copyin (qg) present_or_copyout (vq)
  !$acc parallel loop collapse(2)
  DO nb = 1, nbnt
     DO ig = 1, npw
        qgr = qg(ig)
        px = qgr / dq - INT(qgr/dq)
        ux = 1.0_dp - px
        vx = 2.0_dp - px
        wx = 3.0_dp - px
        i0 = qgr / dq + 1
        i1 = i0 + 1
        i2 = i0 + 2
        i3 = i0 + 3
        IF ( i3 <= nqx ) THEN
            vq(ig,nb) = ( tab_beta(i0,nb,nt) * (-vx*wx-ux*wx-ux*vx)/6.0_dp + &
                          tab_beta(i1,nb,nt) * (+vx*wx-px*wx-px*vx)/2.0_dp - &
                          tab_beta(i2,nb,nt) * (+ux*wx-px*wx-px*ux)/2.0_dp + &
                          tab_beta(i3,nb,nt) * (+ux*vx-px*vx-px*ux)/6.0_dp )/dq
        ELSE
            vq(ig,nb) = 0.0_dp 
        END IF
     ENDDO
  END DO
  !$acc end data
  !----------------------------------------------------------------------
END SUBROUTINE interp_dbeta
