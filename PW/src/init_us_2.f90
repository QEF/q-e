!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE init_us_2( npw_, igk_, q_, vkb_ )
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
  !
  CALL start_clock( 'init_us_2' )
  !
  CALL init_us_2_base(npw_, npwx, igk_, q_, nat, tau, ityp, tpiba, omega, &
          dfftp%nr1, dfftp%nr2, dfftp%nr3, eigts1, eigts2, eigts3, mill, g,&
          vkb_ )
  !
  CALL stop_clock( 'init_us_2' )
  !
END SUBROUTINE init_us_2

!----------------------------------------------------------------------
SUBROUTINE init_us_2_gpu( npw_, igk__d, q_, vkb__d )
  !----------------------------------------------------------------------
  !
  USE kinds,        ONLY : DP
  USE ions_base,    ONLY : nat, ityp, tau
  USE cell_base,    ONLY : tpiba, omega
  USE gvect,        ONLY : eigts1_d, eigts2_d, eigts3_d, mill_d, g_d
  USE wvfct,        ONLY : npwx
  USE uspp,         ONLY : nkb 
  USE fft_base ,    ONLY : dfftp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw_
  INTEGER, INTENT(IN) :: igk__d(npw_)
  REAL(DP), INTENT(IN) :: q_(3)
  COMPLEX(DP), INTENT(OUT) :: vkb__d(npwx,nkb)
#ifdef __CUDA
  attributes(DEVICE) :: igk__d, vkb__d
#endif
  !
  CALL start_clock( 'init_us_2_gpu' )
  !
  CALL init_us_2_base_gpu(npw_, npwx, igk__d, q_, nat, tau, ityp, tpiba, omega,&
    dfftp%nr1, dfftp%nr2, dfftp%nr3, eigts1_d, eigts2_d, eigts3_d, mill_d, g_d,&
    vkb__d )
  !
  CALL stop_clock( 'init_us_2_gpu' )
  !
END SUBROUTINE init_us_2_gpu

