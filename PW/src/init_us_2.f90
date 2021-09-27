!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE uspp_init
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
      !$acc data present(igk_(npw_), mill(:,:), g(:,:), vkb_(npwx,nkb), eigts1(:,:), eigts2(:,:), eigts3(:,:))
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
  !
END MODULE
