  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the subroutine poolscatter in PW/para - Quantum-ESPRESSO group
  !
  !--------------------------------------------------------------------
  subroutine poolgather (nsize, nkstot, nks, f_in, f_out)
  !--------------------------------------------------------------------
  !
  !  gather the kpoints and the electronic eigenvalues
  !  across the pools 
  !  doesn't work with the double grid (k and k+q)
  !
  !--------------------------------------------------------------------
  USE kinds,     only : DP
#ifdef __PARA
  USE mp_global, ONLY : my_pool_id, nproc_pool,    &
                        inter_pool_comm, me_pool, root_pool, kunit,npool, my_pool_id
  USE mp,        ONLY : mp_barrier, mp_bcast,mp_sum
  USE mp_world,  ONLY : mpime
#endif
  implicit none  
  !
  INTEGER :: nsize, nks, nkstot
    ! first dimension of vectors f_in and f_out
    ! number of k-points per pool
    ! total number of k-points
  REAL (KIND=DP) :: f_in(nsize,nks), f_out(nsize,nkstot)
    ! input ( only for k-points of mypool )
    ! output  ( contains values for all k-point )
  !
#ifdef __PARA
  INTEGER :: rest, nbase
  ! the rest of the integer division nkstot / npo
  ! the position in the original list
  !
  rest = nkstot / kunit - ( nkstot / kunit / npool ) * npool 
  !
  nbase = nks * my_pool_id
  !
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest * kunit
  f_out = 0.d0
  f_out(:,(nbase+1):(nbase+nks)) = f_in(:,1:nks)
  !
  ! ... reduce across the pools
  !
  CALL mp_sum(f_out,inter_pool_comm) 
  !
#else
  f_out(:,:) = f_in(:,:)
  !
#endif
  !
  end subroutine poolgather
  !
  !--------------------------------------------------------------------
  subroutine poolgather2 (nsize, nkstot, nks, f_in, f_out)
  !--------------------------------------------------------------------
  !
  !  gather the kpoints and the electronic eigenvalues
  !  across the pools 
  !  works with the double grid (k and k+q)
  !  define rest and nbase as in loadkmesh_para subroutine
  !
  !--------------------------------------------------------------------
  USE kinds,     ONLY : DP
#ifdef __PARA
  USE mp_global, ONLY : my_pool_id, nproc_pool,    &
                        inter_pool_comm, me_pool, root_pool, kunit,npool, my_pool_id
  USE mp,        ONLY : mp_barrier, mp_bcast,mp_sum
  USE mp_world,  ONLY : mpime
#endif
  implicit none
  !
  INTEGER :: nsize, nks, nkstot
    ! first dimension of vectors f_in and f_out
    ! number of k-points per pool
    ! total number of k-points
  REAL (KIND=DP) :: f_in(nsize,nks), f_out(nsize,nkstot)
    ! input ( only for k-points of mypool )
    ! output  ( contains values for all k-point )
  !
#ifdef __PARA
  INTEGER :: rest, nbase, nkst
    ! the rest of the integer division nkstot / npo
    ! the position in the original list
  !
  nkst = 2 * ( nkstot / 2 / npool )
  rest = ( nkstot - nkst * npool ) / 2
  IF (my_pool_id < rest ) THEN
     nkst = nkst + 2
     nbase = my_pool_id*nkst
  ELSE
     nbase = rest*(nkst+2)+(my_pool_id-rest)*nkst
  ENDIF
  f_out = 0.d0
  f_out(:,(nbase+1):(nbase+nks)) = f_in(:,1:nks)
  !
  ! ... reduce across the pools
  !
  CALL mp_sum(f_out,inter_pool_comm)
  !
#else
  f_out(:,:) = f_in(:,:)
  !
#endif
  !
  end subroutine poolgather2
  !
                                                            
 
