  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the SUBROUTINE poolscatter in PW/para - Quantum-ESPRESSO group
  !
  !----------------------------------------------------------------------
  MODULE poolgathering
  !----------------------------------------------------------------------
  !! 
  !! This module contains the routines related to k-point or q-point grid 
  !! generation as well as selection of k/q points. 
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !--------------------------------------------------------------------
    SUBROUTINE poolgather(nsize, nkstot, nks, f_in, f_out)
    !--------------------------------------------------------------------
    !!
    !! Gather the kpoints and the electronic eigenvalues
    !! across the pools 
    !! doesn't work with the double grid (k and k+q)
    !!
    USE kinds,     only : DP
    USE mp_global, ONLY : my_pool_id, inter_pool_comm, kunit,npool, my_pool_id
    USE mp,        ONLY : mp_barrier, mp_bcast,mp_sum
    !
    IMPLICIT NONE  
    !
    INTEGER, INTENT(in) :: nsize
    !! first dimension of vectors f_in and f_out
    INTEGER, INTENT(in) :: nks
    !! number of k-points per pool
    INTEGER, INTENT(in) :: nkstot
    !! total number of k-points
    REAL(KIND = DP), INTENT(in) :: f_in(nsize, nks) 
    !! input ( only for k-points of mypool )
    REAL(KIND = DP), INTENT(out) :: f_out(nsize, nkstot)
    !! output  ( contains values for all k-point )
    !
    ! Local variables
#if defined(__MPI)
    INTEGER :: rest
    !! the rest of the INTEGER division nkstot / npo
    INTEGER :: nbase
    ! the position in the original list
    !
    rest = nkstot / kunit - (nkstot / kunit / npool) * npool 
    !
    nbase = nks * my_pool_id
    !
    IF ((my_pool_id + 1) > rest) nbase = nbase + rest * kunit
    f_out = 0.d0
    f_out(:, (nbase + 1):(nbase + nks)) = f_in(:, 1:nks)
    !
    ! Reduce across the pools
    CALL mp_sum(f_out, inter_pool_comm) 
    !
#else
    f_out(:, :) = f_in(:, :)
    !
#endif
    !
    !--------------------------------------------------------------------
    END SUBROUTINE poolgather
    !--------------------------------------------------------------------
    !
    !--------------------------------------------------------------------
    SUBROUTINE poolgather2(nsize, nkstot, nks, f_in, f_out)
    !--------------------------------------------------------------------
    !!
    !! Gather the kpoints and the electronic eigenvalues
    !! across the pools 
    !! works with the double grid (k and k+q)
    !! define rest and nbase as in loadkmesh_para subroutine
    !!
    !--------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE mp_global, ONLY : my_pool_id,    &
                          inter_pool_comm, npool, my_pool_id
    USE mp,        ONLY : mp_barrier, mp_bcast,mp_sum
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nsize
    !! first dimension of vectors f_in and f_out
    INTEGER, INTENT(in) :: nks
    !! number of k-points per pool
    INTEGER, INTENT(in) :: nkstot
    !! total number of k-points  
    REAL(KIND = DP), INTENT(in) :: f_in(nsize, nks) 
    !! input ( only for k-points of mypool )
    REAL(KIND = DP), INTENT(out) :: f_out(nsize, nkstot)
    !! output  ( contains values for all k-point )
    !
    ! Local variables
#if defined(__MPI)
    INTEGER :: rest
    !! the rest of the INTEGER division nkstot / npo
    INTEGER :: nbase
    !! the position in the original list
    INTEGER :: nkst
    !! Nb de kpt
    !
    nkst = 2 * (nkstot / 2 / npool)
    rest = (nkstot - nkst * npool) / 2
    IF (my_pool_id < rest) THEN
      nkst = nkst + 2
      nbase = my_pool_id * nkst
    ELSE
      nbase = rest * (nkst + 2) + (my_pool_id - rest) * nkst
    ENDIF
    f_out = 0.d0
    f_out(:, (nbase + 1):(nbase + nks)) = f_in(:, 1:nks)
    !
    ! Reduce across the pools
    !
    CALL mp_sum(f_out, inter_pool_comm)
    !
#else
    f_out(:, :) = f_in(:, :)
    !
#endif
    !
    !--------------------------------------------------------------------
    END SUBROUTINE poolgather2
    !--------------------------------------------------------------------
    !
    !--------------------------------------------------------------------
    SUBROUTINE poolgatherc4(nsize1, nsize2, nsize3, nkstot, nks, f_in, f_out)
    !--------------------------------------------------------------------
    !!
    !! gather the kpoints and the electronic eigenvalues
    !! across the pools 
    !! works with the double grid (k and k+q)
    !! define rest and nbase as in loadkmesh_para subroutine
    !!
    !--------------------------------------------------------------------
    USE kinds,     ONLY : DP
    USE mp_global, ONLY : my_pool_id, inter_pool_comm, npool
    USE mp,        ONLY : mp_barrier, mp_bcast,mp_sum
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nsize1
    !! first dimension of vectors f_in and f_out
    INTEGER, INTENT(in) :: nsize2
    !! second dimension of vectors f_in and f_out
    INTEGER, INTENT(in) :: nsize3
    !! third dimension of vectors f_in and f_out
    INTEGER, INTENT(in) :: nks
    !! number of k-points per pool
    INTEGER, INTENT(in) :: nkstot
    !! total number of k-points  
    COMPLEX(KIND = DP), INTENT(in) :: f_in(nsize1, nsize2, nsize3, nks)
    ! input ( only for k-points of mypool )
    COMPLEX(KIND = DP), INTENT(out)  :: f_out(nsize1, nsize2, nsize3, nkstot)
    ! output  ( contains values for all k-point )
    !
    ! Local variables
#if defined(__MPI)
    INTEGER :: rest
    !! the rest of the INTEGER division nkstot / npo
    INTEGER :: nbase
    !! the position in the original list
    INTEGER :: nkst
    !! Nb de kpt
    ! 
    nkst = 2 * (nkstot / 2 / npool)
    rest = (nkstot - nkst * npool) / 2
    IF (my_pool_id < rest) THEN
      nkst = nkst + 2
      nbase = my_pool_id * nkst
    ELSE
      nbase = rest * (nkst + 2) + (my_pool_id - rest) * nkst
    ENDIF
    f_out = 0.d0
    f_out(:, :, :, (nbase + 1):(nbase + nks)) = f_in(:, :, :, 1:nks)
    !
    ! Reduce across the pools
    CALL mp_sum(f_out, inter_pool_comm)
    !
#else
    f_out(:, :, :, :) = f_in(:, :, :, :)
    !
#endif
    !
    !--------------------------------------------------------------------
    END SUBROUTINE poolgatherc4
    !--------------------------------------------------------------------
    !
    !--------------------------------------------------------------------
    SUBROUTINE poolgather_int1(nkstot, nks, f_in, f_out)
    !--------------------------------------------------------------------
    !!
    !! Gather the kpoints and the electronic eigenvalues
    !! across the pools 
    !! works with the double grid (k and k+q)
    !! define rest and nbase as in loadkmesh_para subroutine
    !!
    !--------------------------------------------------------------------
    USE mp_global, ONLY : my_pool_id, inter_pool_comm, kunit,npool, my_pool_id
    USE mp,        ONLY : mp_barrier, mp_bcast,mp_sum
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nks
    !! number of k-points per pool
    INTEGER, INTENT(in) :: nkstot
    !! total number of k-points
    INTEGER, INTENT(in) :: f_in(nks) 
    !! input ( only for k-points of mypool )
    INTEGER, INTENT(out) :: f_out(nkstot)
    !! output  ( contains values for all k-point )
    !
    ! Local variables 
#if defined(__MPI)
    INTEGER :: rest
    !! the rest of the INTEGER division nkstot / npo
    INTEGER :: nbase
    !! the position in the original list
    !
    rest = nkstot / kunit - (nkstot / kunit / npool) * npool
    !
    nbase = nks * my_pool_id
    !
    IF ((my_pool_id + 1) > rest ) nbase = nbase + rest * kunit
    f_out = 0
    f_out((nbase + 1):(nbase + nks)) = f_in(1:nks)
    !
    ! Reduce across the pools
    !
    CALL mp_sum(f_out, inter_pool_comm)
    !
#else
    f_out(:) = f_in(:)
    !
#endif
    !
    !--------------------------------------------------------------------
    END SUBROUTINE poolgather_int1
    !--------------------------------------------------------------------
    !                                                            
    !--------------------------------------------------------------------
    SUBROUTINE poolgather_int(nsize, nkstot, nks, f_in, f_out)
    !--------------------------------------------------------------------
    !!
    !! Gather the kpoints and the electronic eigenvalues
    !! across the pools 
    !! works with the double grid (k and k+q)
    !! define rest and nbase as in loadkmesh_para subroutine
    !!
    !--------------------------------------------------------------------
    USE mp_global, ONLY : my_pool_id, inter_pool_comm, kunit,npool, my_pool_id
    USE mp,        ONLY : mp_barrier, mp_bcast,mp_sum
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nsize
    !! first dimension of vectors f_in and f_out
    INTEGER, INTENT(in) :: nks
    !! number of k-points per pool
    INTEGER, INTENT(in) :: nkstot
    !! total number of k-points
    INTEGER, INTENT(in) :: f_in(nsize, nks)
    !! input ( only for k-points of mypool )
    INTEGER, INTENT(out) :: f_out(nsize, nkstot)
    !! output ( contains values for all k-point )
    !
#if defined(__MPI)
    INTEGER :: rest
    !! the rest of the INTEGER division nkstot / npo
    INTEGER :: nbase
    !! the position in the original list
    !
    rest = nkstot / kunit - (nkstot / kunit / npool) * npool
    !
    nbase = nks * my_pool_id
    !
    IF ((my_pool_id + 1) > rest) nbase = nbase + rest * kunit
    f_out = 0
    f_out(:, (nbase + 1):(nbase + nks)) = f_in(:, 1:nks)
    !
    ! Reduce across the pools
    !
    CALL mp_sum(f_out, inter_pool_comm)
    !
#else
    f_out(:, :) = f_in(:, :)
    !
#endif
    !
    !--------------------------------------------------------------------
    END SUBROUTINE poolgather_int
    !--------------------------------------------------------------------
  !-----------------------------------------------------------------------
  END MODULE poolgathering
  !-----------------------------------------------------------------------


