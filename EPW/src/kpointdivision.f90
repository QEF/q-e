  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !---------------------------------------------------------------------
  subroutine kpointdivision ( ik0 )
  !---------------------------------------------------------------------
  !!
  !! This is just to find the first kpoint block in the pool
  !!
  !---------------------------------------------------------------------
  !
  USE mp_global,   ONLY : my_pool_id,npool
  use pwcom,       ONLY : nkstot
  ! 
  implicit none
  integer :: ik0
  integer :: nkl, nkr, iks
  !
#if defined(__MPI)
  !
  !   number of kpoint blocks, kpoints per pool and reminder
  !
  nkl   = 1 * ( nkstot / npool )
  nkr   = ( nkstot - nkl * npool ) / 1
  !
  !  the reminder goes to the first nkr pools  
  !   (0...nkr-1)
  !
  IF ( my_pool_id < nkr ) nkl = nkl + 1 !kunit
  !
  !  the index of the first k point in this pool
  !
  iks = nkl * my_pool_id + 1
  IF ( my_pool_id >= nkr ) iks = iks + nkr * 1 !kunit
  !
  !  the index of the first k point block in this pool - 1
  !  (I will need the index of ik, not ikk)
  !
  ik0 = ( iks - 1 )
  !
#else
  ik0 = 0
#endif
  !
  end subroutine kpointdivision

  !-----------------------------------------------------------------------
  SUBROUTINE fkbounds( nktot, lower_bnd, upper_bnd )
  !-----------------------------------------------------------------------
  !!
  !!   Subroutine finds the lower and upper bounds a k-grid in parallel
  !!
  !-----------------------------------------------------------------------
  !
  USE mp_global,    ONLY: my_pool_id, npool
  !
  IMPLICIT NONE
  !
  INTEGER :: nktot, lower_bnd, upper_bnd 
  ! input: nktot k-points splited over pools
  ! output : lower_bnd, upper_bnd lower and upper bound for each pool 
  !
#if defined(__MPI)
  !
  INTEGER :: nkl, nkr
  !
  ! find the bounds of k-dependent arrays in the parallel case
  ! number of kpoint blocks, kpoints per pool and reminder
  !
  nkl = nktot / npool
  nkr = nktot - nkl * npool
  !
  ! the reminder goes to the first nkr pools (0...nkr-1)
  !
  IF (my_pool_id < nkr ) nkl = nkl + 1
  !
  ! the index of the first k point in this pool
  !
  lower_bnd = my_pool_id * nkl + 1
  IF ( my_pool_id >= nkr ) lower_bnd = my_pool_id * nkl + 1 + nkr
  !
  ! the index of the last k point in this pool
  !
  upper_bnd = lower_bnd + nkl - 1
  !
#else  
  !     
  ! In serial the definitions are much easier 
  !     
  lower_bnd = 1
  upper_bnd = nktot
  !     
#endif 
  !
  RETURN
  !
  END SUBROUTINE fkbounds
  !
  !-----------------------------------------------------------------------
