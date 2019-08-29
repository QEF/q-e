  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !-----------------------------------------------------------------------
  SUBROUTINE kpoint_grid_epw(nrot, time_reversal, skip_equivalence, s, t_rev, &
                           bg, nk1, nk2, nk3, BZtoIBZ, s_BZtoIBZ)
  !-----------------------------------------------------------------------
  !!
  !!  Automatic generation of a uniform grid of k-points with symmetry. 
  !!  Routine copied from PW/src/kpoint_grid.f90.
  !!  We had to duplicate because the BZtoIBZ array was deallocated and is needed in
  !!  EPW 
  !!
  USE kinds,            ONLY : DP
  USE division,         ONLY : fkbounds
  USE mp,               ONLY : mp_barrier, mp_sum, mp_bcast
  USE mp_world,         ONLY : mpime
  USE mp_global,        ONLY : world_comm, my_pool_id, npool, inter_pool_comm
  USE io_global,        ONLY : ionode_id, stdout
  USE kinds_epw,        ONLY : SIK2
  USE constants_epw,    ONLY : eps6
  #if defined(__MPI)
  USE parallel_include, ONLY : MPI_INTEGER
  #endif
  
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nrot
  !! Number of Bravais symmetry
  INTEGER, INTENT(in) :: nk1, nk2, nk3
  !! k-grid 
  INTEGER, INTENT(in) :: t_rev(48)
  !! Time-reversal sym
  INTEGER, INTENT(in) :: s(3, 3, 48)
  !! Symmetry matrice. 
  INTEGER(SIK2), INTENT(inout) :: s_BZtoIBZ(nk1 * nk2 * nk3)
  !! Symeetry matrix that links an point to its IBZ friend.
  INTEGER, INTENT(inout) :: BZtoIBZ(nk1 * nk2 * nk3)
  !! Number of rotation
  LOGICAL, INTENT(in) :: time_reversal
  !! True if time reversal
  LOGICAL, INTENT(in) :: skip_equivalence
  !! True if equivalent point
  REAL(KIND = DP), INTENT(in) :: bg(3, 3)
  !! Reciprocal space vectors
  !
  ! Local variables
  LOGICAL :: in_the_list
  !! Is the current point in the list
  INTEGER(SIK2) :: s_save(nk1 * nk2 * nk3)
  !! Temporary symmetry matrix
  INTEGER :: nkr
  !! Total number of points
  INTEGER :: i, j, k
  !! Index on grid size
  INTEGER :: ns
  !! Index on symmetry operations
  INTEGER :: n
  !! Global k-point index
  INTEGER :: nk
  !! Equivalent point
  INTEGER :: equiv(nk1 * nk2 * nk3)
  !! Equivalent k-points
  INTEGER :: ik
  !! K-point index 
  INTEGER :: lower_bnd
  !! K-point paralelization (lower-bound index)
  INTEGER :: upper_bnd
  !! K-point paralelization (upper-bound index) 
  INTEGER :: cumul_nks
  !! Sum of points
  INTEGER :: BZtoIBZ_tmp(nk1 * nk2 * nk3)
  !! Temporrary BZtoIBZ map
  INTEGER :: ierr
  !! Error 
  INTEGER, ALLOCATABLE :: nkspar(:)
  !! Number of irr points (IBZ)
  REAL(KIND = DP) :: xkr(3)
  !! Current point
  REAL(KIND = DP) :: xx, yy, zz
  !! Current point coordinate
  REAL(KIND = DP), ALLOCATABLE :: xkg(:, :)
  !! Current point
  REAL(KIND = DP), ALLOCATABLE :: wkk(:)
  !! Weight of the k-point
  !
  ALLOCATE(nkspar(npool))
  nkr = nk1 * nk2 * nk3
  ALLOCATE(xkg(3, nkr))
  ALLOCATE(wkk(nkr))
  equiv(:) = 0
  s_save(:) = 0
  !
  DO i = 1, nk1
     DO j = 1, nk2
        DO k = 1, nk3
           !  this is nothing but consecutive ordering
           n = (k - 1) + ( j- 1 ) * nk3 + (i - 1) * nk2 * nk3 + 1
           !  xkg are the components of the complete grid in crystal axis
           xkg(1, n) = DBLE(i - 1) / nk1 
           xkg(2, n) = DBLE(j - 1) / nk2 
           xkg(3, n) = DBLE(k - 1) / nk3 
        ENDDO
     ENDDO
  ENDDO
  !  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
  !  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)
  DO nk = 1, nkr
    equiv(nk) = nk
  ENDDO
  !
  IF (skip_equivalence) THEN
    CALL infomsg('kpoint_grid', 'ATTENTION: skip check of k-points equivalence')
    wkk = 1.d0
  ELSE
    DO nk = 1, nkr
      !  check if this k-point has already been found equivalent to another
      IF (equiv(nk) == nk) THEN
        wkk(nk) = 1.0d0
        !  check if there are equivalent k-point to this in the list
        !  (excepted those previously found to be equivalent to another)
        !  check both k and -k
        DO ns = 1, nrot
          DO i = 1, 3
            xkr(i) = s(i, 1, ns) * xkg(1, nk) &
                   + s(i, 2, ns) * xkg(2, nk) &
                   + s(i, 3, ns) * xkg(3, nk)
            xkr(i) = xkr(i) - NINT(xkr(i))
          ENDDO
          IF(t_rev(ns) == 1) xkr = -xkr
          xx = xkr(1) * nk1 
          yy = xkr(2) * nk2 
          zz = xkr(3) * nk3 
          in_the_list = ABS(xx - NINT(xx)) <= eps6 .AND. &
                        ABS(yy - NINT(yy)) <= eps6 .AND. &
                        ABS(zz - NINT(zz)) <= eps6
          IF (in_the_list) THEN
            i = MOD(NINT(xkr(1) * nk1 + 2 * nk1), nk1) + 1
            j = MOD(NINT(xkr(2) * nk2 + 2 * nk2), nk2) + 1
            k = MOD(NINT(xkr(3) * nk3 + 2 * nk3), nk3) + 1
            n = (k - 1) + (j - 1) * nk3 + (i - 1) * nk2 * nk3 + 1
            IF (n > nk .AND. equiv(n) == n) THEN
              equiv(n) = nk
              wkk(nk) = wkk(nk) + 1.0d0
              s_save(n) = ns
            ELSE
              IF (equiv(n) /= nk .OR. n < nk ) CALL errore('kpoint_grid', &
                 'something wrong in the checking algorithm', 1)
            ENDIF
          ENDIF
  !        IF (time_reversal) THEN
  !           xx =-xkr(1)*nk1 
  !           yy =-xkr(2)*nk2 
  !           zz =-xkr(3)*nk3 
  !           in_the_list=ABS(xx-NINT(xx))<=eps.AND.ABS(yy-NINT(yy))<=eps &
  !                                              .AND. ABS(zz-NINT(zz))<=eps
  !           IF (in_the_list) THEN
  !              i = mod ( nint (-xkr(1)*nk1  + 2*nk1), nk1 ) + 1
  !              j = mod ( nint (-xkr(2)*nk2  + 2*nk2), nk2 ) + 1
  !              k = mod ( nint (-xkr(3)*nk3  + 2*nk3), nk3 ) + 1
  !              n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
  !              IF (n>nk .AND. equiv(n)==n) THEN
  !                 equiv(n) = nk
  !                 wkk(nk)=wkk(nk)+1.0d0
  !                 s_save(:,:,n) = -s(:,:,ns)
  !              ELSE
  !                 IF (equiv(n)/=nk.OR.n<nk) CALL errore('kpoint_grid', &
  !                 'something wrong in the checking algorithm',2)
  !              ENDIF
  !           ENDIF
  !        ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  ! 
  !  count irreducible points and order them
  nkspar(:) = 0
  DO nk = 1, nkr
    BZtoIBZ(nk) = equiv(nk)
  ENDDO
  !
  CALL fkbounds(nkr, lower_bnd, upper_bnd)
  DO nk = lower_bnd, upper_bnd
    IF (equiv(nk) == nk) THEN
      nkspar(my_pool_id + 1) = nkspar(my_pool_id + 1) + 1
      IF (nkspar(my_pool_id + 1) > nkr) CALL errore('kpoint_grid', 'Too many k-points', 1)
      BZtoIBZ(nk) = nkspar(my_pool_id + 1)
      ! Change all the one above
      DO ik = nk, nkr
        IF (equiv(ik) == nk) THEN
          BZtoIBZ(ik) = nkspar(my_pool_id + 1)
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  !Now recompose the vector with the right order 
  CALL mp_sum(nkspar, world_comm)
  cumul_nks = 0
  IF (my_pool_id > 0) THEN
    DO i = 1, my_pool_id
      cumul_nks = cumul_nks + nkspar(i)
    ENDDO
  ENDIF
  DO ik = 1, nkr
    IF((BZtoIBZ(ik) > nkspar(my_pool_id + 1)) .OR. (ik < lower_bnd)) THEN
      BZtoIBZ (ik) = 0
    ELSE
      BZtoIBZ(ik) = BZtoIBZ(ik) + cumul_nks
    ENDIF
  ENDDO
  BZtoIBZ_tmp(:) = 0
  DO i = 1, npool
    IF (my_pool_id + 1 == i) THEN
      DO ik = 1, nkr
        IF (BZtoIBZ_tmp(ik) == 0) THEN
          BZtoIBZ_tmp(ik) = BZtoIBZ(ik)
        ENDIF
      ENDDO
    ENDIF
    CALL mp_bcast(BZtoIBZ_tmp, i - 1, inter_pool_comm)
  ENDDO
  !
  BZtoIBZ = BZtoIBZ_tmp
  !
  ! Now do the symmetry mapping. 
  DO nk = 1, nkr
    ! If its an irreducible point 
    IF (equiv(nk) == nk ) THEN
      ! Then you have the identity matrix
      s_BZtoIBZ(nk) = 1
    ELSE
      s_BZtoIBZ(nk) = s_save(nk)  
    ENDIF
  ENDDO
  ! 
  DEALLOCATE(xkg)
  DEALLOCATE(wkk)
  DEALLOCATE(nkspar)
  RETURN
  END SUBROUTINE kpoint_grid_epw
  
