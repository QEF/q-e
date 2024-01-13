!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE gen_us_dy_base( npw, npwx, igk, xk, nat, tau, ityp, ntyp, tpiba, &
                           omega, nr1, nr2, nr3, eigts1, eigts2, eigts3,    &
                           mill, g, u, dvkb )
  !----------------------------------------------------------------------
  !! Calculates the Kleinman-Bylander pseudopotentials with the
  !! derivative of the spherical harmonics projected on vector u.
  !
  ! AF: more extensive use of GPU-resident vars possible
  !
  USE upf_kinds,   ONLY: dp
  USE upf_const,   ONLY: tpi
  USE uspp,        ONLY: nkb, indv, nhtol, nhtolm
  USE uspp_param,  ONLY: upf, lmaxkb, nbetam, nh, nhm
  USE beta_mod,    ONLY: interp_beta
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
  REAL(DP), INTENT(IN) :: u(3)
  !! projection vector
  COMPLEX(DP), INTENT(OUT) :: dvkb(npwx, nkb)
  !! the beta function pseudopotential
  !
  ! ... local variables
  !
  INTEGER :: na, nt, nb, ih, l, lm, ikb, iig, ipol, &
             ig, nbm, iq, mil1, mil2, mil3, ikb_t,     &
             nht, ina, lmx2
  !
  INTEGER, ALLOCATABLE :: nas(:), ihv(:), nav(:)
  !
  REAL(DP), ALLOCATABLE :: dylm(:,:,:), dylm_u(:,:)
  REAL(DP), ALLOCATABLE :: q(:), gk(:,:), vkb0(:,:,:)
  ! dylm = d Y_lm/dr_i in cartesian axes
  ! dylm_u as above projected on u
  COMPLEX(DP), ALLOCATABLE :: phase(:), sk(:,:)
  !
  REAL(DP) :: arg, u_ipol1, u_ipol2, u_ipol3, xk1, xk2, xk3
  COMPLEX(DP) :: pref
  !
  !$acc kernels present_or_copyout(dvkb)
  dvkb = (0._DP,0._DP)
  !$acc end kernels
  !
  IF (lmaxkb <= 0) RETURN
  !
  !$acc data present_or_copyin(igk,eigts1,eigts2,eigts3,mill,g) present(dvkb)
  !
  lmx2 = (lmaxkb+1)**2
  !
  ALLOCATE( gk(3,npw) )
  ALLOCATE( dylm_u(npw,lmx2) )
  ALLOCATE( vkb0(npw,nbetam,ntyp) )
  ALLOCATE( q(npw) )
  !$acc data create( dylm_u, vkb0 )
  !$acc data create( q, gk )
  !
  xk1 = xk(1)
  xk2 = xk(2)
  xk3 = xk(3)
  !
  !$acc parallel loop
  DO ig = 1, npw
     iig = igk(ig)
     gk(1,ig) = xk1 + g(1,iig)
     gk(2,ig) = xk2 + g(2,iig)
     gk(3,ig) = xk3 + g(3,iig)
     q(ig) = gk(1,ig)**2 + gk(2,ig)**2 + gk(3,ig)**2
  ENDDO
  !
  ALLOCATE( dylm(npw,(lmaxkb+1)**2,3) )
  !$acc data create( dylm )
  !
  DO ipol = 1, 3
     CALL dylmr2( lmx2, npw, gk, q, dylm(:,:,ipol), ipol )
  ENDDO
  !
  u_ipol1 = u(1) ; u_ipol2 = u(2) ; u_ipol3 = u(3)
  !
  !$acc parallel loop collapse(2)
  DO lm = 1, lmx2
    DO ig = 1, npw
      dylm_u(ig,lm) = u_ipol1*dylm(ig,lm,1) + &
                      u_ipol2*dylm(ig,lm,2) + &
                      u_ipol3*dylm(ig,lm,3)
    ENDDO
  ENDDO
  !$acc end data
  DEALLOCATE( dylm )
  !
  !$acc kernels
  q(:) = SQRT(q(:)) * tpiba
  !$acc end kernels
  !
  DO nt = 1, ntyp
     CALL interp_beta ( nt, npw, q, vkb0(:,:,nt))
  ENDDO
  !
  !$acc end data
  DEALLOCATE( gk, q )
  !
  ALLOCATE( nas(nat), phase(nat) )
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

  ALLOCATE( sk(npw,nat) )
  !$acc data create( sk )
  !
  !$acc data create( phase ) copyin( nas )
  !
  !$acc parallel loop copyin( tau )
  DO ina = 1, nat
     na = nas(ina)
     arg = ( xk1 * tau(1,na) + xk2 * tau(2,na) &
           + xk3 * tau(3,na) ) * tpi
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
  !$acc parallel loop collapse(2) copyin(ityp,indv,nhtol,nhtolm)
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
      dvkb(ig,ikb) = CMPLX(vkb0(ig,nb,nt),KIND=DP) * sk(ig,na) * &
                     CMPLX(dylm_u(ig,lm),KIND=DP)  * pref / CMPLX(tpiba,KIND=DP)
    ENDDO
  ENDDO
  !
  !$acc end data
  !$acc end data
  DEALLOCATE( ihv, nav )
  DEALLOCATE( phase, nas )
  DEALLOCATE( sk )
  !
  IF (ikb_t /= nkb) CALL upf_error( 'gen_us_dy', 'unexpected error', 1 )
  !
  !$acc end data
  DEALLOCATE( dylm_u, vkb0 )
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE gen_us_dy_base
