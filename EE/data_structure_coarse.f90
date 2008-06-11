!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by I. Dabo and N. Marzari (MIT)
!
! contributions by E. Lamas and S. de Gironcoli (SISSA/DEMOCRITOS)
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine data_structure_coarse( lgamma, nr1, nr2, nr3, ecutwfc )
  !-----------------------------------------------------------------------
  ! this routine sets the data structure for the fft arrays.
  ! In the parallel case distributes columns to processes, too
  ! This version computes also the smooth and hard mesh
  !
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : bg,  tpiba
  USE gcoarse,    ONLY : dfftc, nr1c, nr2c, nr3c, nrx1c, nrx2c, nrx3c, nrxxc, &
                         ngmc, ngmc_l, ngmc_g, gcutmc
  USE mp,         ONLY : mp_sum
  USE mp_global,  ONLY : intra_pool_comm, nproc_pool, me_pool, my_image_id
  USE stick_base
  USE fft_scalar, ONLY : good_fft_dimension
  USE fft_types,  ONLY : fft_dlay_allocate, fft_dlay_set, fft_dlay_scalar
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: lgamma
  INTEGER, INTENT(IN) :: nr1, nr2, nr3
  REAL(DP), INTENT(IN) :: ecutwfc
  !
  integer :: n1, n2, n3, i1, i2, i3
  ! counters on G space
  !

  real(DP) :: amod
  ! modulus of G vectors

  integer, allocatable :: st(:,:), stw(:,:), sts(:,:) 
  ! sticks maps

  integer :: ub(3), lb(3)  
  ! upper and lower bounds for maps

  ! cut-off for the wavefunctions

  integer :: np, nps1, nq, nqs, max1, min1, max2, min2, kpoint, m1, &
       m2, i, mc, nct_, ic, ics
  !
  ! ... Sets the dimensions of the coarse grdi
  !
  nrx1c = good_fft_dimension( nr1c )
  nrx2c = nr2c
  nrx3c = nr3c
  nrxxc = nrx1c * nrx2c * nrx3c
  !
  ! ...
  !
  CALL fft_dlay_allocate( dfftc, nproc_pool, MAX( nrx1c, nrx3c ),  nrx2c )
  ! 
  ! ... Computes the number of g necessary to the calculation
  !
  n1 = nr1 + 1
  n2 = nr2 + 1
  n3 = nr3 + 1
  !
  ngmc = 0
  !
  ub =  (/  n1,  n2,  n3 /)
  lb =  (/ -n1, -n2, -n3 /)
  !
  ALLOCATE( stw ( lb( 2 ):ub( 2 ), lb( 3 ):ub( 3 ) ) )
  stw = 0
  !
  DO i1 = - n1, n1
    DO i2 = - n2, n2
      DO i3 = - n3, n3
        amod = (i1 * bg (1, 1) + i2 * bg (1, 2) + i3 * bg (1, 3) ) **2 + &
               (i1 * bg (2, 1) + i2 * bg (2, 2) + i3 * bg (2, 3) ) **2 + &
               (i1 * bg (3, 1) + i2 * bg (3, 2) + i3 * bg (3, 3) ) **2
        IF( amod <= gcutmc )  ngmc  = ngmc  + 1
        IF( amod <= gcutmc )  stw( i2, i3 ) = 1
      END DO
    END DO
  END DO
  !
  CALL fft_dlay_scalar( dfftc, ub, lb, nr1c, nr2c, nr3c, nrx1c, nrx2c, nrx3c, stw )
  !
  DEALLOCATE( stw )
  !
  ! ... 
  !
  ngmc_l = ngmc
  ngmc_g = ngmc
  !
  RETURN
  !

end subroutine data_structure_coarse

