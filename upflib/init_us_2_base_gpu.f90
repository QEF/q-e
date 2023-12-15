!
! Copyright (C) 2021-2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE init_us_2_base_gpu( npw_, npwx, igk_, q_, nat, tau, ityp, &
     tpiba, omega, nr1, nr2, nr3, eigts1, eigts2, eigts3, mill, g, &
     vkb_ )
  !----------------------------------------------------------------------
  !! Calculates beta functions (Kleinman-Bylander projectors), with
  !! structure factor, for all atoms, in reciprocal space.
  !
  USE upf_kinds,    ONLY : DP
  USE upf_const,    ONLY : tpi
  USE uspp,         ONLY : nkb, nhtol, nhtolm, indv
  USE uspp_param,   ONLY : upf, lmaxkb, nhm, nh, nsp
  !
  implicit none
  !
  INTEGER,  INTENT(IN) :: npw_
  !! number of plane wave
  INTEGER,  INTENT(IN) :: npwx
  !! leading dim of vkb_
  INTEGER,  INTENT(IN) :: igk_(npw_)
  !! indices of G in the list of q+G vectors
  REAL(dp), INTENT(IN) :: q_(3)
  !! q vector (2pi/a units)
  INTEGER, INTENT(IN) :: nat
  !! number of atoms
  INTEGER, INTENT(IN) :: ityp(nat)
  !! index of type per atom
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! atomic positions (cc alat units)
  REAL(DP), INTENT(IN) :: tpiba, omega
  !! reclat units and cell volume
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
  COMPLEX(dp), INTENT(OUT) :: vkb_(npwx, nkb)
  !! beta functions (npw_ <= npwx)
  !
  !     Local variables
  !
  integer :: ig, lm, na, nt, nb, ih, jkb
  integer :: iv_d
  real(DP) :: arg, q1, q2, q3

  complex(DP) :: phase, pref
  real(DP), allocatable :: gk (:,:), qg (:), ylm(:,:), vq(:), vkb1(:,:)
  complex(DP), allocatable:: sk(:)
  !
  CALL start_clock( 'init_us_2:gpu' )
  !
  if (lmaxkb<0) return
  
  allocate (vkb1( npw_,nhm))
  allocate (  sk( npw_))
  allocate (  qg( npw_))
  allocate (  vq( npw_))
  allocate ( ylm( npw_, (lmaxkb + 1) **2))
  allocate (  gk( 3, npw_))
  !
  q1 = q_(1)
  q2 = q_(2)
  q3 = q_(3)

  !$acc data create(qg, gk, ylm, vq, vkb1, sk) present(g, igk_, eigts1, eigts2, eigts3, mill, vkb_) 
  !$acc parallel loop
  do ig = 1, npw_
     iv_d = igk_(ig)
     gk (1,ig) = q1 + g(1, iv_d )
     gk (2,ig) = q2 + g(2, iv_d )
     gk (3,ig) = q3 + g(3, iv_d )
     qg (ig) = gk(1, ig)*gk(1, ig) + &
               gk(2, ig)*gk(2, ig) + &
               gk(3, ig)*gk(3, ig)
  enddo
  !
  !$acc host_data use_device (gk, qg, ylm)
  call ylmr2_gpu ((lmaxkb+1)**2, npw_, gk, qg, ylm)
  !$acc end host_data 
  !
  ! set now qg=|q+G| in atomic units
  !
  !$acc parallel loop
  do ig = 1, npw_
     qg(ig) = sqrt(qg(ig))*tpiba
  enddo

  ! |beta_lm(q)> = (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
  jkb = 0
  do nt = 1, nsp     
     do nb = 1, upf(nt)%nbeta
        CALL interp_beta ( nt, nb, npw_, qg, vq )
        ! add spherical harmonic part  (Y_lm(q)*f_l(q)) 
        do ih = 1, nh (nt)
           if (nb.eq.indv (ih, nt) ) then
              lm =nhtolm (ih, nt)
              !$acc parallel loop
              do ig = 1, npw_
                 vkb1 (ig,ih) = ylm (ig, lm) * vq (ig)
              enddo
           endif
        enddo
     enddo
     !
     ! vkb1 contains all betas including angular part for type nt
     ! now add the structure factor and factor (-i)^l
     !
     do na = 1, nat
        ! ordering: first all betas for atoms of type 1
        !           then  all betas for atoms of type 2  and so on
        if (ityp (na) .eq.nt) then
           arg = (q1 * tau (1, na) + &
                  q2 * tau (2, na) + &
                  q3 * tau (3, na) ) * tpi
           phase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
           !
           !$acc parallel loop
           do ig = 1, npw_
              iv_d = igk_(ig)
              sk (ig) = eigts1 (mill(1,iv_d), na) * &
                        eigts2 (mill(2,iv_d), na) * &
                        eigts3 (mill(3,iv_d), na)
           enddo
           !
           do ih = 1, nh (nt)
              jkb = jkb + 1
              !l = nhtol (ih, nt)
              pref = (0.d0, -1.d0) **nhtol (ih, nt) * phase
              !$acc parallel loop
              do ig = 1, npw_
                 vkb_(ig, jkb) = vkb1 (ig,ih) * sk (ig) * pref
              enddo
              !$acc parallel loop
              do ig = npw_+1, npwx
                 vkb_(ig, jkb) = (0.0_dp, 0.0_dp)
              enddo
           enddo
        endif
     enddo
  enddo
  !$acc end data

  deallocate(gk)
  deallocate(ylm)
  deallocate(vq)
  deallocate(qg)
  deallocate(sk)
  deallocate(vkb1)
  !
  CALL stop_clock( 'init_us_2:gpu' )
  !
  return
end subroutine init_us_2_base_gpu
!
!-----------------------------------------------------------------------
SUBROUTINE interp_beta ( nt, nb, npw, qg, vq )
  !-----------------------------------------------------------------------
  !
  ! computes vq: radial fourier transform of beta functions
  !
  USE upf_kinds,  ONLY : dp
  USE uspp_data,  ONLY : dq, tab_beta
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: nt
  INTEGER, INTENT(IN)  :: npw
  INTEGER, INTENT(IN)  :: nb
  REAL(dp), INTENT(IN) :: qg(npw)
  REAL(dp), INTENT(OUT):: vq(npw)
  !
  INTEGER :: ig
  INTEGER :: i0, i1, i2, i3
  REAL(dp):: qgr, px, ux, vx, wx
  !
  !$acc data present(tab_beta, qg, vq)
  !$acc parallel loop
  DO ig = 1, npw
     qgr = qg(ig)
     px = qgr / dq - DBLE(INT(qgr/dq))
     ux = 1.d0 - px
     vx = 2.d0 - px
     wx = 3.d0 - px
     i0 = INT(qgr/dq) + 1
     i1 = i0 + 1
     i2 = i0 + 2
     i3 = i0 + 3
     vq(ig) = &
          tab_beta(i0,nb,nt) * ux * vx * wx / 6.d0 + &
          tab_beta(i1,nb,nt) * px * vx * wx / 2.d0 - &
          tab_beta(i2,nb,nt) * px * ux * wx / 2.d0 + &
          tab_beta(i3,nb,nt) * px * ux * vx / 6.d0
  END DO
  !$acc end data
END SUBROUTINE interp_beta
