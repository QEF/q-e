  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------------
  subroutine loadumat ( nbnd, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq )
  !----------------------------------------------------------------------------
  !!
  !!   wannier interpolation of e-p vertex:
  !!   load rotation matrix on coarse mesh and distribute
  !!
  !!   07/2010 JN changed the way this was done.  Previously
  !!   a pool scatter call that didn't work on hbar was performed
  !!   and some crazy packing scheme was needed to use poolscatter
  !!   subsequently it is just a bcast followed by an appropriate assignment
  !! 
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE klist_epw,     ONLY : kmap   
  USE epwcom,        ONLY : filukk
  USE constants_epw, ONLY : czero
  USE io_epw,        ONLY : iunukk
  USE io_global, ONLY : ionode_id, meta_ionode
  USE mp_global, ONLY : inter_pool_comm
  USE mp,        ONLY : mp_sum, mp_barrier, mp_bcast
  USE mp_world,  ONLY : mpime
  implicit none
  ! 
  LOGICAL, INTENT (out) :: lwin( nbnd, nks )
  !! Band windows at k
  LOGICAL, INTENT (out) :: lwinq( nbnd, nks )
  !! Band windows at k+q
  !
  INTEGER, INTENT (in) :: nbnd
  !! Number of bands
  INTEGER, INTENT (in) :: nbndsub
  !! number of bands in the optimal subspace
  INTEGER, INTENT (in) :: nks
  !! number of kpoints 
  INTEGER, INTENT (in) :: nkstot
  !! total number of kpoints across pools
  ! 
  REAL(kind=DP), INTENT (in) :: xxq(3)
  !! the qpoint for folding of U
  !
  COMPLEX(kind=DP), INTENT (out) :: cu ( nbnd, nbndsub, nks)
  !! U matrix
  COMPLEX(kind=DP), INTENT (out) :: cuq ( nbnd, nbndsub, nks)
  !! U(k+q) matrix
  ! 
  ! work variables 
  !
  integer :: ik, ibnd, jbnd, ios, ik_start, ik_stop
  complex(kind=DP) :: cu_big ( nbnd, nbndsub, nkstot), cuq_big ( nbnd, nbndsub, nkstot)
  logical :: lwin_big( nbnd, nkstot ), lwinq_big( nbnd, nkstot )
  !
  cu_big = czero
  cuq_big = czero
  IF (meta_ionode) then
    !
    ! first proc read rotation matrix (coarse mesh) from file
    !
    open ( unit = iunukk, file = filukk, status = 'old', form = 'formatted',iostat=ios)
    IF (ios /=0) call errore ('loadumat', 'error opening ukk file',iunukk)

    !
    DO ik = 1, nkstot
      DO ibnd = 1, nbnd
        DO jbnd = 1, nbndsub
           READ(iunukk, *) cu_big (ibnd, jbnd, ik)
        ENDDO
      ENDDO
    ENDDO
    DO ik = 1,nkstot
       DO ibnd = 1,nbnd
          READ (iunukk,*) lwin_big(ibnd,ik)
       ENDDO
    ENDDO
    !
    close ( iunukk )
    !
    !  generate U(k+q) through the map 
    !
    !
    ! here we create the map k+q --> k
    ! (first proc has a copy of all kpoints)
    !
    !  generates kmap(ik) for this xxq
    !
    CALL createkmap2 ( xxq )
    !
    !  and we generate the matrix for the q-displaced mesh
    !
    DO ik = 1, nkstot
       cuq_big (:, :, ik) = cu_big (:, :, kmap(ik) )
       lwinq_big (:, ik) = lwin_big (:, kmap(ik) )
    ENDDO
    !
  ENDIF
  CALL mp_sum (cu_big,  inter_pool_comm)
  CALL mp_sum (cuq_big, inter_pool_comm)
  CALL mp_bcast (lwin_big, ionode_id, inter_pool_comm)
  CALL mp_bcast (lwinq_big, ionode_id, inter_pool_comm)
  !
  CALL ckbounds(ik_start, ik_stop)
  IF ( (ik_stop-ik_start+1) .ne. nks) call errore('loadumat',"Improper parallel ukk load",1)
  cu = cu_big (:, :, ik_start:ik_stop)
  cuq = cuq_big (:, :, ik_start:ik_stop)
  lwin = lwin_big (:, ik_start:ik_stop)
  lwinq = lwin_big (:, ik_start:ik_stop)
  !
 end subroutine loadumat

