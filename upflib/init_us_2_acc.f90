!
! Copyright (C) 2021-2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE init_us_2_acc( npw_, npwx, igk_, q_, nat, tau, ityp, &
     tpiba, omega, nr1, nr2, nr3, eigts1, eigts2, eigts3, mill, g, &
     vkb_ )
  !----------------------------------------------------------------------
  !! Calculates beta functions (Kleinman-Bylander projectors), with
  !! structure factor, for all atoms, in reciprocal space. NOTE:
  !! - variables ending with _ may vary during a run
  !! - see acc data present declaration below for ACC-required definitions
  !
  USE upf_kinds,    ONLY : DP
  USE upf_const,    ONLY : tpi
  USE uspp,         ONLY : nkb, nhtol, nhtolm, indv
  USE uspp_param,   ONLY : lmaxkb, nbetam, nhm, nh, nsp
  USE beta_mod,     ONLY : interp_beta
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
  !! k or k+q vector (2pi/a units)
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
  integer :: ig, lm, na, nt, nb, ih, ikb, jkb, nhnt
  integer :: iv_d, l
  real(DP) :: arg, q1, q2, q3

  complex(dp) :: pref
  !complex(dp) :: pref(0:4) = [ ( 1.0_dp,  0.0_dp), &
  !        &                    ( 0.0_dp, -1.0_dp), &
  !        &                    (-1.0_dp,  0.0_dp), &
  !        &                    ( 0.0_dp,  1.0_dp), &
  !        &                    ( 1.0_dp,  0.0_dp)] ! (-i)^l 
  real(DP), allocatable :: gk (:,:), qg (:), ylm(:,:), vq(:,:), vkb1(:,:)
  complex(DP), allocatable:: sk(:)
  !
  if (lmaxkb < 0) return
  !
  allocate (vkb1( npw_,nhm))
  allocate (  sk( npw_))
  allocate (  qg( npw_))
  allocate (  vq( npw_, nbetam ))
  allocate ( ylm( npw_, (lmaxkb + 1) **2))
  allocate (  gk( 3, npw_))
  !
  q1 = q_(1)
  q2 = q_(2)
  q3 = q_(3)
  !
  !$acc data create(qg, gk, ylm, vq, vkb1, sk) &
  !$acc      present(g, igk_, eigts1, eigts2, eigts3, mill, vkb_) &
  !$acc      copyin(nhtol, nhtolm, indv)
  !$acc kernels
  vkb_(:,:) = (0.0_dp, 0.0_dp)
  !$acc end kernels
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
  call ylmr2 ((lmaxkb+1)**2, npw_, gk, qg, ylm)
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
     !
     CALL interp_beta ( nt, npw_, qg, vq )
     ! add spherical harmonic part  (Y_lm(q)*f_l(q)) 
     nhnt = nh(nt)
     !$acc parallel loop collapse(2)
     do ih = 1, nhnt
        do ig = 1, npw_
           nb = indv (ih, nt)
           lm = nhtolm (ih, nt)
           vkb1 (ig,ih) = ylm (ig, lm) * vq (ig,nb)
        enddo
     enddo
     !
     ! vkb1 contains all betas including angular part for type nt
     ! now add the structure factor and factor (-i)^l
     !
     do na = 1, nat
        ! ordering: first all betas for atoms of type 1
        !           then  all betas for atoms of type 2  and so on
        if (ityp (na) == nt) then
           arg = (q1 * tau (1, na) + &
                  q2 * tau (2, na) + &
                  q3 * tau (3, na) ) * tpi
           !
           !$acc parallel loop
           do ig = 1, npw_
              iv_d = igk_(ig)
              sk (ig) = eigts1 (mill(1,iv_d), na) * &
                        eigts2 (mill(2,iv_d), na) * &
                        eigts3 (mill(3,iv_d), na) * &
                        CMPLX(cos (arg), -sin (arg) ,kind=DP)
           enddo
           ! 
           !$acc parallel loop collapse(2)
           do ih = 1, nhnt
              do ig = 1, npw_
                 !l = nhtol (ih, nt)
                 pref = (0.d0, -1.d0) **nhtol (ih, nt)
                 vkb_(ig, jkb+ih) = vkb1 (ig,ih) * sk (ig) * pref
              enddo
           enddo
           jkb=jkb+nhnt
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
  return
  !----------------------------------------------------------------------
end subroutine init_us_2_acc
!----------------------------------------------------------------------
