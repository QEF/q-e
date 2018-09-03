!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine init_us_2_gpu (npw_, igk__d, q_, vkb__d)
  !----------------------------------------------------------------------
  !
  !   Calculates beta functions (Kleinman-Bylander projectors), with
  !   structure factor, for all atoms, in reciprocal space. On input:
  !      npw_       : number of PWs 
  !      igk_(npw_) : indices of G in the list of q+G vectors
  !      q_(3)      : q vector (2pi/a units)
  !  On output:
  !      vkb_(npwx,nkb) : beta functions (npw_ <= npwx)
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,  ONLY : tpiba
  USE constants,  ONLY : tpi
  USE gvect_gpum, ONLY : eigts1_d, eigts2_d, eigts3_d, mill_d, g_d
  USE wvfct,      ONLY : npwx
  USE us,         ONLY : nqx, dq, spline_ps
  USE us_gpum,    ONLY : tab_d, tab_d2y_d
  USE m_gth,      ONLY : mk_ffnl_gth
  USE splinelib_gpum, ONLY : splint_eq_gpu
  USE uspp,       ONLY : nkb, nhtol, nhtolm, indv
  USE uspp_param, ONLY : upf, lmaxkb, nhm, nh
!  USE ylmr2_gpum, ONLY : ylmr2_gpu
  !
  USE us_gpum,    ONLY : using_tab_d, using_tab_d2y_d
  USE gbuffers,   ONLY : dev_buf
#if defined(__CUDA)
  USE cudafor
#endif
  !
  implicit none
  !
  INTEGER, INTENT (IN) :: npw_
  INTEGER, INTENT (IN) :: igk__d (npw_)
  REAL(dp), INTENT(IN) :: q_(3)
  COMPLEX(dp), INTENT(OUT) :: vkb__d (npwx, nkb)
#if defined(__CUDA)
  attributes(DEVICE) :: igk__d, vkb__d
#endif
  !
  !     Local variables
  !
  integer :: i0,i1,i2,i3, ig, lm, na, nt, nb, ih, jkb
  integer :: istat
  integer :: iv_d
  real(DP) :: px, ux, vx, wx, arg, q1, q2, q3
  real(DP), allocatable :: gk_d (:,:), qg_d (:), vq_d(:), ylm_d(:,:), vkb1_d(:,:)
  real(DP), allocatable :: qg_h (:), vq_h(:)
  real(DP) :: rv_d

  complex(DP) :: phase, pref
  complex(DP), allocatable :: sk_d(:)

  logical :: is_gth
  integer :: iq
#if defined(__CUDA)
  attributes(DEVICE) :: gk_d, qg_d, vq_d, ylm_d, vkb1_d, sk_d
  attributes(PINNED) :: qg_h, vq_h
#endif
  !
  !
  if (lmaxkb.lt.0) return
  call start_clock ('init_us_2')
  
  call using_tab_d(0)
  if (spline_ps) call using_tab_d2y_d(0)

  ! JR Eventually replace with smarter allocation/deallocation of GPU temp arrays
  ! PB use buffer class here
  allocate (vkb1_d( npw_,nhm))
  allocate (  sk_d( npw_))
  allocate (  qg_d( npw_))
  allocate (  vq_d( npw_))
  allocate ( ylm_d( npw_, (lmaxkb + 1) **2))
  allocate (  gk_d( 3, npw_))
  !CALL dev_buf%lock_buffer(vkb1_d, (/ npw_,nhm/), istat )
  !CALL dev_buf%lock_buffer(  sk_d, npw_, istat )
  !CALL dev_buf%lock_buffer(  qg_d, npw_, istat )
  !CALL dev_buf%lock_buffer(  vq_d, npw_, istat )
  !CALL dev_buf%lock_buffer( ylm_d, (/ npw_, (lmaxkb + 1) **2 /), istat )
  !CALL dev_buf%lock_buffer(  gk_d, (/ 3, npw_ /), istat )

  is_gth = .false.
  do nt = 1, ntyp
     is_gth = upf(nt)%is_gth
     if (is_gth) then
        allocate (  qg_h( npw_))    
        allocate (  vq_h( npw_)) 
        is_gth = .true.
        exit
     end if
  end do
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

  ! JR Don't need this when using splint_eq_gpu
  !if (spline_ps) then
  !  allocate(xdata(nqx))
  !  do iq = 1, nqx
  !    xdata(iq) = (iq - 1) * dq
  !  enddo
  !endif

  ! |beta_lm(q)> = (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
  jkb = 0
  do nt = 1, ntyp
     do nb = 1, upf(nt)%nbeta
        if ( upf(nt)%is_gth ) then
           qg_h = qg_d
           call mk_ffnl_gth( nt, nb, npw_, qg_h, vq_h )
           vq_d = vq_h
        else if (spline_ps) then
           call splint_eq_gpu(dq, tab_d(:,nb,nt), tab_d2y_d(:,nb,nt), qg_d, vq_d)
        else
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
                          
              !vq_d (ig) = tab_d (i0, nb, nt) * ux * vx * wx / 6.d0 + &
              !            tab_d (i1, nb, nt) * px * vx * wx / 2.d0 - &
              !            tab_d (i2, nb, nt) * px * ux * wx / 2.d0 + &
              !            tab_d (i3, nb, nt) * px * ux * vx / 6.d0
           enddo
        endif

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

  deallocate(gk_d)
  deallocate(ylm_d)
  deallocate(vq_d)
  deallocate(qg_d)
  deallocate(sk_d)
  deallocate(vkb1_d)
  !CALL dev_buf%release_buffer(vkb1_d, istat )
  !CALL dev_buf%release_buffer(  sk_d, istat )
  !CALL dev_buf%release_buffer(  qg_d, istat )
  !CALL dev_buf%release_buffer(  vq_d, istat )
  !CALL dev_buf%release_buffer( ylm_d, istat )
  !CALL dev_buf%release_buffer(  gk_d, istat )
  IF (is_gth) THEN
     deallocate ( qg_h, vq_h )
  END IF

  call stop_clock ('init_us_2')
  return
end subroutine init_us_2_gpu
