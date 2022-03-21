!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE init_us_2_base_gpu( npw_, npwx, igk__d, q_, nat, tau, ityp, &
     tpiba, omega, nr1, nr2, nr3, eigts1_d, eigts2_d, eigts3_d, mill_d, g_d, &
     vkb__d )
  !----------------------------------------------------------------------
  !! Calculates beta functions (Kleinman-Bylander projectors), with
  !! structure factor, for all atoms, in reciprocal space.
  !
  USE upf_kinds,    ONLY : DP
  USE upf_const,    ONLY : tpi
  USE uspp_data,    ONLY : nqx, dq, tab_d
  USE uspp,         ONLY : nkb, nhtol, nhtolm, indv
  USE uspp_param,   ONLY : upf, lmaxkb, nhm, nh, nsp
  USE device_fbuff_m,   ONLY : dev_buf
  !
  implicit none
  !
  INTEGER,  INTENT(IN) :: npw_
  INTEGER,  INTENT(IN) :: npwx
  !! leading dim of vkb_
  INTEGER,  INTENT(IN) :: igk__d(npw_)
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
  COMPLEX(dp), INTENT(OUT) :: vkb__d(npwx, nkb)
  !! beta functions (npw_ <= npwx)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: igk__d, vkb__d
  attributes(DEVICE) :: eigts1_d, eigts2_d, eigts3_d
  attributes(DEVICE) :: mill_d, g_d
#endif
  !
  !     Local variables
  !
  integer :: i0,i1,i2,i3, ig, lm, na, nt, nb, ih, jkb
  integer :: istat(6)
  integer :: iv_d
  real(DP) :: px, ux, vx, wx, arg, q1, q2, q3
  real(DP), pointer :: gk_d (:,:), qg_d (:), vq_d(:), ylm_d(:,:), vkb1_d(:,:)
  real(DP) :: rv_d

  complex(DP) :: phase, pref
  complex(DP), pointer :: sk_d(:)

  integer :: iq
#if defined(__CUDA)
  attributes(DEVICE) :: gk_d, qg_d, vq_d, ylm_d, vkb1_d, sk_d
  !
  CALL start_clock( 'init_us_2:gpu' )
  !
  if (lmaxkb<0) return
  
  ! JR Eventually replace with smarter allocation/deallocation of GPU temp arrays
  ! PB use buffer class here
  !allocate (vkb1_d( npw_,nhm))
  !allocate (  sk_d( npw_))
  !allocate (  qg_d( npw_))
  !allocate (  vq_d( npw_))
  !allocate ( ylm_d( npw_, (lmaxkb + 1) **2))
  !allocate (  gk_d( 3, npw_))
  CALL dev_buf%lock_buffer(vkb1_d, (/ npw_, nhm /), istat(1) )
  CALL dev_buf%lock_buffer(  sk_d, npw_, istat(2) )
  CALL dev_buf%lock_buffer(  qg_d, npw_, istat(3) )
  CALL dev_buf%lock_buffer(  vq_d, npw_, istat(4) )
  CALL dev_buf%lock_buffer( ylm_d, (/ npw_, (lmaxkb + 1) **2 /), istat(5) )
  CALL dev_buf%lock_buffer(  gk_d, (/ 3, npw_ /), istat(6) )
  IF (ANY(istat /= 0)) CALL upf_error( 'init_us_2_gpu', 'cannot allocate buffers', -1 )
  !
  q1 = q_(1)
  q2 = q_(2)
  q3 = q_(3)

  !$cuf kernel do(1) <<<*,*>>>
  do ig = 1, npw_
     iv_d = igk__d(ig)
     gk_d (1,ig) = q1 + g_d(1, iv_d )
     gk_d (2,ig) = q2 + g_d(2, iv_d )
     gk_d (3,ig) = q3 + g_d(3, iv_d )
     qg_d (ig) = gk_d(1, ig)*gk_d(1, ig) + &
                 gk_d(2, ig)*gk_d(2, ig) + &
                 gk_d(3, ig)*gk_d(3, ig)
  enddo
  !
  call ylmr2_gpu ((lmaxkb+1)**2, npw_, gk_d, qg_d, ylm_d)
  !
  ! set now qg=|q+G| in atomic units
  !
  !$cuf kernel do(1) <<<*,*>>>
  do ig = 1, npw_
     qg_d(ig) = sqrt(qg_d(ig))*tpiba
  enddo

  ! |beta_lm(q)> = (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
  jkb = 0
  do nt = 1, nsp
     do nb = 1, upf(nt)%nbeta
        !$cuf kernel do(1) <<<*,*>>>
        do ig = 1, npw_
           rv_d = qg_d(ig)
           px = rv_d / dq - int (rv_d / dq)
           ux = 1.d0 - px
           vx = 2.d0 - px
           wx = 3.d0 - px
           i0 = INT( rv_d / dq ) + 1
           i1 = i0 + 1
           i2 = i0 + 2
           i3 = i0 + 3
           vq_d (ig) = ux * vx * (wx * tab_d(i0, nb, nt) + px * tab_d(i3, nb, nt)) / 6.d0 + &
                       px * wx * (vx * tab_d(i1, nb, nt) - ux * tab_d(i2, nb, nt)) * 0.5d0
                          
        enddo

        ! add spherical harmonic part  (Y_lm(q)*f_l(q)) 
        do ih = 1, nh (nt)
           if (nb.eq.indv (ih, nt) ) then
              !l = nhtol (ih, nt)
              lm =nhtolm (ih, nt)

              !$cuf kernel do(1) <<<*,*>>>
              do ig = 1, npw_
                 vkb1_d (ig,ih) = ylm_d (ig, lm) * vq_d (ig)
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
           arg = (q_(1) * tau (1, na) + &
                  q_(2) * tau (2, na) + &
                  q_(3) * tau (3, na) ) * tpi
           phase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
           !
           !$cuf kernel do(1) <<<*,*>>>
           do ig = 1, npw_
              sk_d (ig) = eigts1_d (mill_d(1,igk__d(ig)), na) * &
                          eigts2_d (mill_d(2,igk__d(ig)), na) * &
                          eigts3_d (mill_d(3,igk__d(ig)), na)
           enddo
           !
           do ih = 1, nh (nt)
              jkb = jkb + 1
              pref = (0.d0, -1.d0) **nhtol (ih, nt) * phase
              !$cuf kernel do(1) <<<*,*>>>
              do ig = 1, npw_
                 vkb__d(ig, jkb) = vkb1_d (ig,ih) * sk_d (ig) * pref
              enddo
              !$cuf kernel do(1) <<<*,*>>>
              do ig = npw_+1, npwx
                 vkb__d(ig, jkb) = (0.0_dp, 0.0_dp)
              enddo
           enddo
        endif
     enddo
  enddo

  !deallocate(gk_d)
  !deallocate(ylm_d)
  !deallocate(vq_d)
  !deallocate(qg_d)
  !deallocate(sk_d)
  !deallocate(vkb1_d)
  CALL dev_buf%release_buffer(vkb1_d, istat(1) )
  CALL dev_buf%release_buffer(  sk_d, istat(2) )
  CALL dev_buf%release_buffer(  qg_d, istat(3) )
  CALL dev_buf%release_buffer(  vq_d, istat(4) )
  CALL dev_buf%release_buffer( ylm_d, istat(5) )
  CALL dev_buf%release_buffer(  gk_d, istat(6) )
  !
  CALL stop_clock( 'init_us_2:gpu' )
#endif
  !
  return
end subroutine init_us_2_base_gpu
