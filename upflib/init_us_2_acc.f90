!
! Copyright (C) 2021-2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE init_us_2_acc( npw, npwx, igk, q, nat, tau, ityp, &
     tpiba, omega, nr1, nr2, nr3, eigts1, eigts2, eigts3, mill, g, &
     vkb )
  !----------------------------------------------------------------------
  !! Calculates beta functions (Kleinman-Bylander projectors), with
  !! structure factor, for all atoms, in reciprocal space.
  !
  USE upf_kinds,   ONLY: dp
  USE upf_const,   ONLY: tpi
  USE uspp,        ONLY: nkb, indv, nhtol, nhtolm
  USE uspp_data,   ONLY: nqx, tab_beta, dq
  USE uspp_param,  ONLY: upf, lmaxkb, nbetam, nh, nhm, nsp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw
  !! number of plane waves
  INTEGER, INTENT(IN) :: npwx
  !! max number of plane waves and leading dimension of vkb
  INTEGER, INTENT(IN) :: igk(npw)
  !! indices of G in the list of q+G vectors
  REAL(dp), INTENT(IN) :: q(3)
  !! q vector (2pi/a units)
  INTEGER, INTENT(IN) :: nat
  !! number of atoms
  INTEGER, INTENT(IN) :: ityp(nat)
  !! index of type per atom
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! atomic positions (cc alat units)
  REAL(DP), INTENT(IN) :: tpiba
  !! rec.lattice units 2pi/a
  REAL(DP), INTENT(IN) :: omega
  !! cell volume
  INTEGER, INTENT(IN) :: nr1,nr2,nr3
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
  COMPLEX(DP), INTENT(OUT) :: vkb(npwx, nkb)
  !! the beta functions (npw <= npwx)
  !
  ! ... local variables
  !
  INTEGER :: na, nt, nb, ih, l, lm, ikb, iig, ipol, i0, i1, i2, &
             i3, ig, nbm, iq, mil1, mil2, mil3, ikb_t,     &
             nht, ina, lmx2
  !
  INTEGER, ALLOCATABLE :: nas(:), ihv(:), nav(:)
  !
  REAL(DP), ALLOCATABLE :: ylm(:,:)
  REAL(DP), ALLOCATABLE :: qg(:), gk(:,:), vkb0(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: phase(:), sk(:,:)
  !
  REAL(DP) :: px, ux, vx, wx, arg, q1, q2, q3
  COMPLEX(DP) :: pref
  !
  !$acc kernels present_or_copyout(vkb)
  vkb = (0._DP,0._DP)
  !$acc end kernels
  !
  IF (lmaxkb < 0) RETURN
  !
  !$acc data present( igk, eigts1, eigts2, eigts3, mill, g, vkb )
  !
  lmx2 = (lmaxkb+1)**2
  !
  ALLOCATE( gk(3,npw) )
  ALLOCATE( ylm(npw,lmx2) )
  ALLOCATE( vkb0(npw,nbetam,nsp) )
  ALLOCATE( qg(npw) )
  !$acc data create( ylm, vkb0 )
  !$acc data create( qg, gk )
  !
  q1 = q(1)
  q2 = q(2)
  q3 = q(3)
  !
  !$acc parallel loop
  DO ig = 1, npw
     iig = igk(ig)
     gk(1,ig) = q1 + g(1,iig)
     gk(2,ig) = q2 + g(2,iig)
     gk(3,ig) = q3 + g(3,iig)
     qg(ig) = gk(1,ig)**2 + gk(2,ig)**2 + gk(3,ig)**2
  ENDDO
  !
#if defined(__CUDA)
  !$acc host_data use_device( gk, qg, ylm )
  CALL ylmr2_gpu( lmx2, npw, gk, qg, ylm )
  !$acc end host_data
#else
  !$acc update self( gk, qg )
  CALL ylmr2( lmx2, npw, gk, qg, ylm )
  !$acc update device( ylm )
#endif
  !
  !$acc kernels
  qg(:) = SQRT(qg(:)) * tpiba
  !$acc end kernels
  !
  !$acc data present ( tab_beta )
  DO nt = 1, nsp
     nbm = upf(nt)%nbeta
     !$acc parallel loop collapse(2)
     DO nb = 1, nbm
        DO ig = 1, npw
           px = qg(ig)/dq - DBLE(INT(qg(ig)/dq))
           ux = 1._DP - px
           vx = 2._DP - px
           wx = 3._DP - px
           i0 = INT(qg(ig)/dq) + 1
           i1 = i0 + 1
           i2 = i0 + 2
           i3 = i0 + 3
           vkb0(ig,nb,nt) = tab_beta(i0,nb,nt) * ux * vx * wx / 6._DP + &
                            tab_beta(i1,nb,nt) * px * vx * wx / 2._DP - &
                            tab_beta(i2,nb,nt) * px * ux * wx / 2._DP + &
                            tab_beta(i3,nb,nt) * px * ux * vx / 6._DP
       ENDDO
    ENDDO
  ENDDO
  !$acc end data
  !
  !$acc end data
  DEALLOCATE( gk, qg )
  !
  ALLOCATE( nas(nat), phase(nat) )
  !
  ina = 0
  DO nt = 1, nsp
     DO na = 1, nat
        IF ( ityp(na) == nt ) THEN
           ina = ina + 1
           nas(ina) = na
        ENDIF
     ENDDO
  ENDDO
  !
  ALLOCATE( sk(npw,nat) )
  !$acc data create( sk )
  !
  !$acc data create( phase ) copyin( nas )
  !
  !$acc parallel loop copyin( tau )
  DO ina = 1, nat
     na = nas(ina)
     arg = ( q1 * tau(1,na) &
           + q2 * tau(2,na) &
           + q3 * tau(3,na) ) * tpi
     phase(na) = CMPLX( COS(arg), -SIN(arg), KIND=DP )
  ENDDO
  !
  !$acc parallel loop collapse(2)
  DO ina = 1, nat
    DO ig = 1, npw
      !
      na = nas(ina)
      iig = igk(ig)
      mil1 = mill(1,iig)
      mil2 = mill(2,iig)
      mil3 = mill(3,iig)
      sk(ig,na) = eigts1(mil1,na) * &
                  eigts2(mil2,na) * &
                  eigts3(mil3,na) * phase(na)
    ENDDO
  ENDDO
  !
  !$acc end data
  !
  ALLOCATE( ihv(nat*nhm), nav(nat*nhm) )
  !$acc data create( ihv, nav )
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
  IF (ikb_t /= nkb) CALL upf_error( 'init_us_2', 'unexpected error', 1 )
  !
  !$acc parallel loop collapse(2) copyin(ityp,indv,nhtol,nhtolm)
  DO ikb = 1, nkb
    DO ig = 1, npw
      ih = ihv(ikb)
      na = nav(ikb)
      nt = ityp(na)
      nb = indv(ih,nt)
      l  = nhtol(ih,nt)
      lm = nhtolm(ih,nt)
      pref = (0._DP,-1._DP)**l
      !
      vkb(ig,ikb) = CMPLX(vkb0(ig,nb,nt),KIND=DP) * sk(ig,na) * &
                    CMPLX(ylm(ig,lm),KIND=DP)  * pref
    ENDDO
  ENDDO
  !
  !$acc end data
  !$acc end data
  DEALLOCATE( ihv, nav )
  DEALLOCATE( phase, nas )
  DEALLOCATE( sk )
  !
  !$acc end data
  DEALLOCATE( vkb0 )
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE init_us_2_acc
