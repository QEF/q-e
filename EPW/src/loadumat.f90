  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------------
  SUBROUTINE loadumat(nbndep, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq, &
                      exband, w_centers)
  !----------------------------------------------------------------------------
  !!
  !! Wannier interpolation of e-p vertex:
  !! load rotation matrix on coarse mesh and distribute
  !!
  !! 07/2010 JN changed the way this was done.  Previously
  !! a pool scatter call that didn't work on hbar was performed
  !! and some crazy packing scheme was needed to use poolscatter
  !! subsequently it is just a bcast followed by an appropriate assignment
  !!
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE klist_epw,     ONLY : kmap
  USE epwcom,        ONLY : filukk
  USE constants_epw, ONLY : czero, zero
  USE io_var,        ONLY : iunukk
  USE io_global,     ONLY : meta_ionode_id, meta_ionode
  USE mp,            ONLY : mp_sum, mp_barrier, mp_bcast
  USE division,      ONLY : fkbounds
  USE elph2,         ONLY : xkq
  USE kfold,         ONLY : createkmap2
  USE pwcom,         ONLY : nbnd
  USE mp_world,      ONLY : world_comm
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(out) :: lwin(nbndep, nks)
  !! Band windows at k
  LOGICAL, INTENT(out) :: lwinq(nbndep, nks)
  !! Band windows at k+q
  LOGICAL, INTENT(out) :: exband(nbnd)
  !! Band excluded
  INTEGER, INTENT(in) :: nbndep
  !! Number of bands
  INTEGER, INTENT(in) :: nbndsub
  !! number of bands in the optimal subspace
  INTEGER, INTENT(in) :: nks
  !! number of kpoints
  INTEGER, INTENT(in) :: nkstot
  !! total number of kpoints across pools
  REAL(KIND = DP), INTENT(in) :: xxq(3)
  !! the qpoint for folding of U
  REAL(KIND = DP), INTENT(inout) :: w_centers(3, nbndsub)
  !! Wannier centers
  COMPLEX(KIND = DP), INTENT(out) :: cu(nbndep, nbndsub, nks)
  !! U(k) matrix for k-points in the pool
  COMPLEX(KIND = DP), INTENT(out) :: cuq(nbndep, nbndsub, nks)
  !! U(k+q) matrix for k+q-points in the pool
  !
  ! Local variables
  LOGICAL :: lwin_big(nbndep, nkstot)
  !! .TRUE. if the band ibnd lies within the outer window at k-point ik
  LOGICAL :: lwinq_big(nbndep, nkstot)
  !! .TRUE. if the band ibnd lies within the outer window at k+qpoint ikq
  INTEGER :: ik
  !! Counter of k-point index
  INTEGER :: iw
  !! Counter on Wannier centers
  INTEGER :: ibnd
  !! Counter on band index
  INTEGER :: jbnd
  !! Counter on wannierized bands
  INTEGER :: ios
  !! INTEGER variable for I/O control
  INTEGER :: ik_start
  !! Index of first k-point in the pool
  INTEGER :: ik_stop
  !! Index of last k-point in the pool
  COMPLEX(KIND = DP) :: cu_big(nbndep, nbndsub, nkstot)
  !! U(k) matrix for all k-points
  COMPLEX(KIND = DP) :: cuq_big(nbndep, nbndsub, nkstot)
  !! U(k+q) matrix for all k+q-points
  !
  cu_big = czero
  cuq_big = czero
  IF (meta_ionode) THEN
    !
    ! First proc read rotation matrix (coarse mesh) from file
    !
    OPEN(iunukk, FILE = filukk, STATUS = 'old', FORM = 'formatted', IOSTAT = ios)
    IF (ios /=0) CALL errore('loadumat', 'error opening ukk file', iunukk)
    !
    ! dummy operation for skipping unnecessary data (ibndstart and ibndend) here
    !
    READ(iunukk, *) ibnd, jbnd
    !
    DO ik = 1, nkstot
      DO ibnd = 1, nbndep
        DO jbnd = 1, nbndsub
          READ(iunukk, *) cu_big(ibnd, jbnd, ik)
        ENDDO
      ENDDO
    ENDDO
    DO ik = 1, nkstot
      DO ibnd = 1, nbndep
        READ(iunukk, *) lwin_big(ibnd, ik)
      ENDDO
    ENDDO
    DO ibnd = 1, nbnd
      READ(iunukk, *) exband(ibnd)
    ENDDO
    ! Read the Wannier centers
    DO iw = 1, nbndsub
      READ(iunukk, *) w_centers(:, iw)
    ENDDO
    !
    CLOSE(iunukk)
    !
    ! Generate U(k+q) through the map
    ! here we create the map k+q --> k
    ! (first proc has a copy of all kpoints)
    ! Generates kmap(ik) for this xxq
    !
    xkq(:, :) = zero
    CALL createkmap2(xxq)
    !
    ! and we generate the matrix for the q-displaced mesh
    !
    DO ik = 1, nkstot
      cuq_big(:, :, ik) = cu_big(:, :, kmap(ik))
      lwinq_big(:, ik) = lwin_big(:, kmap(ik))
    ENDDO
  ENDIF ! meta_ionode
  !
  CALL mp_bcast(cu_big, meta_ionode_id, world_comm)
  CALL mp_bcast(cuq_big, meta_ionode_id, world_comm)
  CALL mp_bcast(lwin_big, meta_ionode_id, world_comm)
  CALL mp_bcast(lwinq_big, meta_ionode_id, world_comm)
  CALL mp_bcast(exband, meta_ionode_id, world_comm)
  CALL mp_bcast(w_centers, meta_ionode_id, world_comm)
  !
  CALL fkbounds(nkstot, ik_start, ik_stop)
  !
  IF ((ik_stop - ik_start + 1) /= nks) CALL errore('loadumat', "Improper parallel ukk load", 1)
  !
  cu = cu_big(:, :, ik_start:ik_stop)
  cuq = cuq_big(:, :, ik_start:ik_stop)
  lwin = lwin_big(:, ik_start:ik_stop)
  lwinq = lwinq_big(:, ik_start:ik_stop)
  !
  RETURN
  !-----------------------------------------------------------------------
  END SUBROUTINE loadumat
  !-----------------------------------------------------------------------
