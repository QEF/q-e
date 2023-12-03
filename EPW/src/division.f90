  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE division
  !----------------------------------------------------------------------
  !!
  !! This module contains various core splitting routines
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !--------------------------------------------------------------------
    SUBROUTINE para_bounds(lower, upper, total)
    !--------------------------------------------------------------------
    !!
    !! Subroutine finds the lower and upper
    !! bounds if we split some quantity over pools
    !!
    !---------------------------------
    !
    USE mp_global,   ONLY : my_pool_id, npool
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: lower
    !! Lower bound
    INTEGER, INTENT(out) :: upper
    !! Upper bound
    INTEGER, INTENT(in)  :: total
    !! Total quantity
    !
#if defined(__MPI)
    !
    ! Local variables
    INTEGER :: rest
    !! Remaining
    INTEGER :: nrst
    !!
    !
    IF (total <= npool) THEN
      IF (my_pool_id < total) THEN
        lower = my_pool_id + 1
        upper = lower
      ELSE
        lower = 1
        upper  = 0
      ENDIF
    ELSE
      nrst = total / npool
      rest = total - nrst * npool
      IF (my_pool_id < rest) THEN
        nrst = nrst + 1
        lower = my_pool_id * nrst + 1
        upper = lower + nrst - 1
      ELSE
        lower = rest * (nrst + 1) + (my_pool_id - rest) * nrst + 1
        upper = lower + nrst - 1
      ENDIF
    ENDIF
#else
    lower = 1
    upper = total
#endif
    !
    ! --------------------------------------------------------------------
    END SUBROUTINE para_bounds
    !---------------------------------------------------------------------
    !
    !---------------------------------------------------------------------
    SUBROUTINE kpointdivision(ik0)
    !---------------------------------------------------------------------
    !!
    !! This is just to find the first kpoint block in the pool
    !!
    !---------------------------------------------------------------------
    !
    USE mp_global,   ONLY : my_pool_id,npool
    USE pwcom,       ONLY : nkstot
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(out) :: ik0
    !! Return the first kpoint in the pool
    !
    ! Local variables
    INTEGER :: nkl
    !! Number of kpoints block
    INTEGER :: nkr
    !! Reminder
    INTEGER :: iks
    !! Index of the first kpoint in the pool
    !
#if defined(__MPI)
    !
    ! Number of kpoint blocks, kpoints per pool and reminder
    !
    nkl = 1 * (nkstot / npool)
    nkr = (nkstot - nkl * npool) / 1
    !
    ! The reminder goes to the first nkr pools (0...nkr-1)
    !
    IF (my_pool_id < nkr) nkl = nkl + 1 !kunit
    !
    ! The index of the first k point in this pool
    !
    iks = nkl * my_pool_id + 1
    IF (my_pool_id >= nkr) iks = iks + nkr * 1 !kunit
    !
    !  the index of the first k point block in this pool - 1
    !  (I will need the index of ik, not ikk)
    !
    ik0 = (iks - 1)
#else
    ik0 = 0
#endif
    !-----------------------------------------------------------------------
    END SUBROUTINE kpointdivision
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE fkbounds(nktot, lower_bnd, upper_bnd)
    !-----------------------------------------------------------------------
    !!
    !!   Subroutine finds the lower and upper bounds a k-grid in parallel
    !!
    !! @ Note:
    !!    If you have 19 kpts and 2 pool, this routine will return
    !!    lower_bnd= 1 and upper_bnd=10 for the first pool
    !!    lower_bnd=11 and upper_bnd=19 for the second pool
    !-----------------------------------------------------------------------
    !
    USE mp_global,    ONLY : my_pool_id, npool
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nktot
    !! nktot k-points splited over pools
    INTEGER, INTENT(out) :: lower_bnd
    !! Lower kpt bounds for that image pool
    INTEGER, INTENT(out) :: upper_bnd
    !! Upper kpt for that image pool
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
    IF (my_pool_id >= nkr ) lower_bnd = my_pool_id * nkl + 1 + nkr
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
    !-----------------------------------------------------------------------
    END SUBROUTINE fkbounds
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE fkbounds2(nktot, lower_bnd, upper_bnd)
    !-----------------------------------------------------------------------
    !!
    !!   Subroutine finds the lower and upper bounds a k-grid in parallel
    !!
    !! @ Note:
    !!    If you have 19 kpts and 2 pool, this routine will return
    !!    lower_bnd= 1 and upper_bnd=10 for the first pool
    !!    lower_bnd=11 and upper_bnd=19 for the second pool
    !-----------------------------------------------------------------------
    !
    USE mp_global,        ONLY : my_pool_id, npool
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_OFFSET_KIND
#endif
    !
    IMPLICIT NONE
    !
#if defined(__MPI)
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(in) :: nktot
    !! nktot k-points splited over pools
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(out) :: lower_bnd
    !! Lower kpt bounds for that image pool
    INTEGER(KIND = MPI_OFFSET_KIND), INTENT(out) :: upper_bnd
    !! Upper kpt for that image pool
#else
    INTEGER(KIND = 8), INTENT(in)  :: nktot
    !! nktot k-points splited over pools
    INTEGER(KIND = 8), INTENT(out) :: lower_bnd
    !! Lower kpt bounds for that image pool
    INTEGER(KIND = 8), INTENT(out) :: upper_bnd
    !! Upper kpt for that image pool
#endif
    !
#if defined(__MPI)
    !
    ! Local variables
    INTEGER(KIND = MPI_OFFSET_KIND) :: nkl
    !! Number of k-point block
    INTEGER(KIND = MPI_OFFSET_KIND) :: nkr
    !! Reminder
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
    IF (my_pool_id >= nkr ) lower_bnd = my_pool_id * nkl + 1 + nkr
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
    !-----------------------------------------------------------------------
    END SUBROUTINE fkbounds2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE fkbounds_bnd(nbnd, lower_bnd, upper_bnd)
    !-----------------------------------------------------------------------
    !!
    !!   Subroutine finds the lower and upper bounds in band parallelization
    !!
    !! @ Note:
    !!    If you have 20 bands and 2 images, this routine will return
    !!    lower_bnd= 1 and upper_bnd=10 for the first pool
    !!    lower_bnd= 11 and upper_bnd=19 for the second pool
    !-----------------------------------------------------------------------
    !
    USE mp_images,  ONLY : nimage, my_image_id
    ! number of images, number of processor within an image, index of the proc within an image
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: nbnd
    !! Total number of band to be splitted among images
    INTEGER, INTENT(out) :: lower_bnd
    !! Lower band bounds for that image pool
    INTEGER, INTENT(out) :: upper_bnd
    !! Upper band bound for that image pool
    !
#if defined(__MPI)
    !
    ! Local variables
    INTEGER :: nkl
    !! ADD
    INTEGER :: nkr
    !!
    !
    ! find the bounds of k-dependent arrays in the parallel case
    ! number of kpoint blocks, kpoints per pool and reminder
    !
    nkl = nbnd / nimage
    nkr = nbnd - nkl * nimage
    !
    ! the reminder goes to the first nkr pools (0...nkr-1)
    !
    IF (my_image_id < nkr ) nkl = nkl + 1
    !
    ! the index of the first k point in this pool
    !
    lower_bnd = my_image_id * nkl + 1
    IF (my_image_id >= nkr ) lower_bnd = my_image_id * nkl + 1 + nkr
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
    upper_bnd = nbnd
    !
#endif
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE fkbounds_bnd
    !-----------------------------------------------------------------------
    !
  !-----------------------------------------------------------------------
  END MODULE division
  !-----------------------------------------------------------------------

