  !                                                                            
  ! Copyright (C) 2010-2018 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! 
  !-----------------------------------------------------------------------
  SUBROUTINE load_rebal
  !-----------------------------------------------------------------------
  !
  ! Routine used to rebalance the load on k-points.  
  ! At the moment this routines is only called in the case of IBTE
  ! using k-point symmetry and an homogeneous grid. 
  ! The k-point that are within the fshick are equally splitted among cores. 
  ! 
  !-----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE elph2,         ONLY : etf, nkf, nkqtotf, xkf, wkf, etf, map_rebal, map_rebal_inv
  USE epwcom,        ONLY : fsthick, nbndsub, mp_mesh_k
  USE pwcom,         ONLY : ef
  USE mp_global,     ONLY : my_pool_id, npool
  USE constants_epw, ONLY : zero
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : inter_pool_comm, mp_bcast
  USE transportcom,  ONLY : lower_bnd
  !
  implicit none

  INTEGER :: pool_index(npool) 
  !! Index of the current pool
  INTEGER :: ipool
  !! Pool loop
  INTEGER :: ik, ikk, ikq, ikpt, ikpt2, maxin, maxout, rest, tot
  !! Total k-point index
  INTEGER :: kpt_in(nkqtotf)
  !! K-points that are within the fshick windows
  INTEGER :: kpt_out(nkqtotf)
  !! K-points that are outside of the fshick windows
  INTEGER :: map_rebal_tmp(nkqtotf/2)
  !! Temporary map between the initial ordering of k-point and the rebalanced one
  INTEGER :: map_rebal_inv_tmp(nkqtotf/2)
  !! Temporary inverse map between the initial ordering of k-point and the rebalanced one
  !
  REAL(kind=DP) :: xkf_all(3,nkqtotf)
  !! Collect k-point coordinate (and k+q) from all pools in parallel case
  REAL(kind=DP) :: wkf_all(nkqtotf)
  !! Collect k-point coordinate (and k+q) from all pools in parallel case
  REAL(kind=DP) :: etf_all(nbndsub, nkqtotf)
  !! Collect all the eigenenergies
  !
  ! Gather all the k-point coordinate and all the weights from all the pools
  ! Gather also all the eigeneneriges
  xkf_all(:,:) = zero
  wkf_all(:) = zero
  etf_all(:,:) = zero
  ! 
#ifdef __MPI
  CALL poolgather2 ( 1, nkqtotf, 2*nkf, wkf, wkf_all )
  CALL poolgather2 ( 3, nkqtotf, 2*nkf, xkf, xkf_all )
  CALL poolgather2 ( nbndsub, nkqtotf, 2*nkf, etf, etf_all )
#else
  xkf_all = xkf
  wkf_all = wkf
  etf_all = etf
#endif 
  ! 
  ALLOCATE(map_rebal(nkqtotf/2))
  ALLOCATE(map_rebal_inv(nkqtotf/2))
  ! 
  kpt_in(:) = 0 
  kpt_out(:) = 0 
  ! 
  !The sorting is done by master only
  IF ( my_pool_id == 0 ) THEN
    ! 
    ikpt = 0
    ikpt2 = 0
    ! 
    DO ik = 1, nkqtotf/2
      ikk = 2 * ik - 1
      ikq = ikk + 1
      IF ( minval ( abs(etf_all(:, ikk) - ef) ) < fsthick ) THEN
        ikpt = ikpt + 1  
        kpt_in( ikpt) = ikk
      ELSE 
        ikpt2 = ikpt2 + 1
        kpt_out( ikpt2 ) = ikk
      ENDIF ! fsthick
      ! 
    ENDDO ! ik
    ! In case of IBZ create a map of old to new points
    IF (mp_mesh_k) THEN
      map_rebal(:) = 0
      map_rebal_inv(:) = 0
      DO ik = 1, ikpt
        map_rebal( (kpt_in(ik) + 1) / 2 ) = ik
        map_rebal_inv(ik) = (kpt_in(ik) + 1) / 2 
      ENDDO
      DO ik = 1, ikpt2
        map_rebal( ( kpt_out(ik) + 1 ) /2 ) =  ikpt + ik
        map_rebal_inv( ikpt + ik ) = ( kpt_out(ik) + 1 ) /2  
      ENDDO
    ENDIF
  ENDIF ! my_pool_id
  ! 
  CALL mp_bcast (kpt_in, ionode_id, inter_pool_comm)
  CALL mp_bcast (kpt_out, ionode_id, inter_pool_comm)
  CALL mp_bcast (ikpt, ionode_id, inter_pool_comm)
  CALL mp_bcast (ikpt2, ionode_id, inter_pool_comm)
  CALL mp_bcast (map_rebal, ionode_id, inter_pool_comm)
  CALL mp_bcast (map_rebal_inv, ionode_id, inter_pool_comm)
  ! 
  ! map_rebal contains an array with all the k-point in the window and then
  ! all the k-points out of the window. 
  ! We then split those k-points such that the first core has the first k-point, 
  ! the second core has the second k-point etc 
  ! 
  tot = ( nkqtotf / (2 * npool) )
  rest = ( (nkqtotf/2) - tot * npool )
  ! 
  DO ipool=1, npool
    DO ik=1,  tot
      map_rebal_inv_tmp(ik + (ipool - 1) * tot) = map_rebal_inv( npool * ik - (npool-ipool) )
    ENDDO
  ENDDO
  ! Do the rest
  DO ik=1, rest
    map_rebal_inv_tmp(ik + npool * tot) = map_rebal_inv( npool * tot + ik )
  ENDDO
  map_rebal_inv(:) = map_rebal_inv_tmp(:) 
  ! 
  ! Now recontruct map_rebal so that it is the inverse mapping as map_rebal_inv
  DO ik=1,nkqtotf/2
    map_rebal( map_rebal_inv(ik) ) = ik
  ENDDO
  ! 
  ! We then assign this new order to the local k-points and weights on each cpu
  DO ik=1,nkf
    ikk = 2 * ik - 1
    xkf(:,ikk)   = xkf_all(:,2 * map_rebal_inv( ik + lower_bnd - 1  )-1 )
    xkf(:,ikk+1) = xkf_all(:,2 * map_rebal_inv( ik + lower_bnd - 1  )  )
    wkf(ikk)   = wkf_all( 2 * map_rebal_inv( ik + lower_bnd - 1 )-1)
    wkf(ikk+1) = wkf_all( 2 * map_rebal_inv( ik + lower_bnd - 1 ) )
  ENDDO
  !
  END SUBROUTINE load_rebal

