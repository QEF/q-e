!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE uspp_init
  !
  PRIVATE
  PUBLIC :: init_us_2, gen_us_dj, gen_us_dy
  PUBLIC :: gen_us_dj_gpu, gen_us_dy_gpu
  !
CONTAINS
  !----------------------------------------------------------------------
  SUBROUTINE init_us_2( npw_, igk_, q_, vkb_, run_on_gpu_)
    !----------------------------------------------------------------------
    !! wrapper to call init_us_2_base
    !! Calculates beta functions (Kleinman-Bylander projectors), with
    !! structure factor, for all atoms, in reciprocal space.
    !
    USE kinds,        ONLY : DP
    USE ions_base,    ONLY : nat, ityp, tau
    USE cell_base,    ONLY : tpiba, omega
    USE gvect,        ONLY : eigts1, eigts2, eigts3, mill, g
    USE wvfct,        ONLY : npwx
    USE uspp,         ONLY : nkb
    USE fft_base ,    ONLY : dfftp
    USE control_flags, ONLY : use_gpu
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: npw_
    !! number of PWs 
    INTEGER, INTENT(IN) :: igk_(npw_)
    !! indices of G in the list of q+G vectors
    REAL(DP), INTENT(IN) :: q_(3)
    !! q vector (2pi/a units)
    COMPLEX(DP), INTENT(OUT) :: vkb_(npwx,nkb)
    !! beta functions (npw_ <= npwx)
    LOGICAL, OPTIONAL, INTENT(IN) :: run_on_gpu_
    !! whether you wish to run on gpu in case use_gpu is true 
    !!
    LOGICAL :: run_on_gpu 
    !
    CALL start_clock( 'init_us_2' )
    !
    run_on_gpu = .false.
    if(present(run_on_gpu_)) run_on_gpu = run_on_gpu_
    !
    if(use_gpu.and.run_on_gpu) then   
      !
      !$acc data present(igk_(1:npw_), mill(:,:), g(:,:), vkb_(1:npwx,1:nkb), eigts1(:,:), eigts2(:,:), eigts3(:,:))
      !$acc host_data use_device(eigts1, eigts2, eigts3, mill, g, igk_, vkb_)
      CALL init_us_2_base_gpu(npw_, npwx, igk_, q_, nat, tau, ityp, tpiba, omega,&
        dfftp%nr1, dfftp%nr2, dfftp%nr3, eigts1, eigts2, eigts3, mill, g,&
        vkb_ )
      !$acc end host_data
      !$acc end data
      !
    else
      CALL init_us_2_base(npw_, npwx, igk_, q_, nat, tau, ityp, tpiba, omega, &
              dfftp%nr1, dfftp%nr2, dfftp%nr3, eigts1, eigts2, eigts3, mill, g,&
              vkb_ )
    end if 
    !
    CALL stop_clock( 'init_us_2' )
    !
  END SUBROUTINE init_us_2
  !
  !----------------------------------------------------------------------
  SUBROUTINE gen_us_dj ( ik, dvkb )
    !----------------------------------------------------------------------
    !! wrapper to call gen_us_dj: same as init_us_2, but with the derivative
    !! of the Bessel functions dj_l/dq instead of j_l(qr) in the integral
    !
    USE kinds,        ONLY : DP
    USE ions_base,    ONLY : nat, ntyp=>nsp, ityp, tau
    USE cell_base,    ONLY : tpiba, omega
    USE gvect,        ONLY : eigts1, eigts2, eigts3, mill, g
    USE wvfct,        ONLY : npwx
    USE uspp,         ONLY : nkb
    USE fft_base ,    ONLY : dfftp
    USE klist,        ONLY : xk, ngk, igk_k
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ik
    !! k-point index
    COMPLEX(DP), INTENT(OUT) :: dvkb(npwx,nkb)
    !! beta functions computed with dj_l/dq
    !
    ! CALL start_clock( 'gen_us_dj' )
    !
    CALL gen_us_dj_base (ngk(ik), npwx, igk_k(1,ik), xk(1,ik), nat, tau, &
            ityp, ntyp, tpiba, omega, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
            eigts1, eigts2, eigts3, mill, g, dvkb )
    !
    ! CALL stop_clock( 'gen_us_dj' )
    !
  END SUBROUTINE gen_us_dj
  !
  !----------------------------------------------------------------------
  SUBROUTINE gen_us_dy ( ik, u, dvkb )
    !----------------------------------------------------------------------
    !! wrapper to call gen_us_dj: same as init_us_2, but with the derivative
    !! of the spherical harmonics dY_lm/dq instead of Y_lm(q) in the integral
    !
    USE kinds,        ONLY : DP
    USE ions_base,    ONLY : nat, ntyp=>nsp, ityp, tau
    USE cell_base,    ONLY : tpiba, omega
    USE gvect,        ONLY : eigts1, eigts2, eigts3, mill, g
    USE wvfct,        ONLY : npwx
    USE uspp,         ONLY : nkb
    USE fft_base ,    ONLY : dfftp
    USE klist,        ONLY : xk, ngk, igk_k
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ik
    !! k-point index
    REAL(dp), INTENT(IN) :: u(3)
    !! k-point index
    COMPLEX(DP), INTENT(OUT) :: dvkb(npwx,nkb)
    !! beta functions computed with dY_lm/dq
    !
    ! CALL start_clock( 'gen_us_dy' )
    !
    CALL gen_us_dy_base (ngk(ik), npwx, igk_k(1,ik), xk(1,ik), nat, tau, &
            ityp, ntyp, tpiba, omega, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
            eigts1, eigts2, eigts3, mill, g, u, dvkb )
    !
    ! CALL stop_clock( 'gen_us_d' )
    !
  END SUBROUTINE gen_us_dy
  !
  !----------------------------------------------------------------------
  SUBROUTINE gen_us_dj_gpu ( ik, dvkb )
    !----------------------------------------------------------------------
    !! wrapper to call gen_us_dj: same as init_us_2, but with the derivative
    !! of the Bessel functions dj_l/dq instead of j_l(qr) in the integral
    !! (GPU version, FIXME: to be merged)
    !
    USE kinds,        ONLY : dp
    USE ions_base,    ONLY : nat, ntyp=>nsp, ityp, tau
    USE cell_base,    ONLY : tpiba, omega
    USE gvect,        ONLY : eigts1_d, eigts2_d, eigts3_d, mill_d, g_d
    USE wvfct,        ONLY : npwx
    USE uspp,         ONLY : nkb
    USE fft_base ,    ONLY : dfftp
    USE klist,        ONLY : xk, ngk, igk_k_d
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ik
    !! k-point index
    COMPLEX(DP), INTENT(OUT) :: dvkb(npwx,nkb)
    !! beta functions computed with dj_l/dq
#if defined(__CUDA)
  attributes(DEVICE) :: dvkb
#endif
    !
    ! CALL start_clock( 'gen_us_dj' )
    !
    CALL gen_us_dj_gpu_ (ngk(ik), npwx, igk_k_d(1,ik), xk(1,ik), nat, tau, &
            ityp, ntyp, tpiba, omega, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
            eigts1_d, eigts2_d, eigts3_d, mill_d, g_d, dvkb )
    !
    ! CALL stop_clock( 'gen_us_dj' )
    !
  END SUBROUTINE gen_us_dj_gpu
  !
  !----------------------------------------------------------------------
  SUBROUTINE gen_us_dy_gpu ( ik, u, dvkb )
    !----------------------------------------------------------------------
    !! wrapper to call gen_us_dj: same as init_us_2, but with the derivative
    !! of the spherical harmonics dY_lm/dq instead of Y_lm(q) in the integral
    !! (GPU version, FIXME: to be merged)
    !
    USE kinds,        ONLY : DP
    USE ions_base,    ONLY : nat, ntyp=>nsp, ityp, tau
    USE cell_base,    ONLY : tpiba, omega
    USE gvect,        ONLY : eigts1_d, eigts2_d, eigts3_d, mill_d, g_d
    USE wvfct,        ONLY : npwx
    USE uspp,         ONLY : nkb
    USE fft_base ,    ONLY : dfftp
    USE klist,        ONLY : xk, ngk, igk_k_d
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ik
    !! k-point index
    REAL(dp), INTENT(IN) :: u(3)
    !! k-point index
    COMPLEX(DP), INTENT(OUT) :: dvkb(npwx,nkb)
    !! beta functions computed with dY_lm/dq
#if defined(__CUDA)
  attributes(DEVICE) :: dvkb
#endif
    !
    ! CALL start_clock( 'gen_us_dy' )
    !
    CALL gen_us_dy_gpu_ (ngk(ik), npwx, igk_k_d(1,ik), xk(1,ik), nat, tau, &
            ityp, ntyp, tpiba, omega, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
            eigts1_d, eigts2_d, eigts3_d, mill_d, g_d, u, dvkb )
    !
    ! CALL stop_clock( 'gen_us_dy' )
    !
  END SUBROUTINE gen_us_dy_gpu
  !
END MODULE
