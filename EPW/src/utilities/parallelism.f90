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
  MODULE parallelism
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
    !--------------------------------------------------------------------
    SUBROUTINE image_division(totq, selecq, image_array, size_image_arr, size_image)
    !--------------------------------------------------------------------
    !!
    !! Added by S. Tiwari
    !! Subroutine finds the q-points for each image
    !!
    !---------------------------------
    !
    USE mp_global,   ONLY : my_pool_id, inter_image_comm
    USE mp_images,   ONLY : my_image_id, nimage
    USE mp,          ONLY : mp_barrier
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: totq
    !! Total number of q-points inside the window
    INTEGER, INTENT(in) :: selecq(totq)
    !! q-points inside the window
    INTEGER, INTENT(inout) :: size_image_arr(nimage)
    !! Number of q-points in each image
    INTEGER, INTENT(inout), ALLOCATABLE :: image_array(:)
    !! q-points in each image
    INTEGER, INTENT(out) :: size_image
    !! Total number of q-points inside the image
    !
    !! Size if q-points in each images 
    INTEGER :: image
    !! Counter over images
    INTEGER :: iq
    !! Counter over q-points
    INTEGER :: ierr
    !! Index for error
    !
    size_image = CEILING(REAL(totq) / nimage) 
    ALLOCATE(image_array(size_image), STAT = ierr)
    IF (ierr /= 0) CALL errore("image_division", "image_array not built", 1)
    image_array(:) = -1
    image = 1   
    size_image_arr(:) = 0 
    DO iq = 1, totq
      IF (MOD(iq, nimage) == my_image_id) THEN
        image_array(image) = selecq(iq)
        size_image_arr(my_image_id+1) = image
        image = image + 1
      ENDIF
    ENDDO    
    size_image = image-1
    ! 
    !---------------------------------------------------------------------
    END SUBROUTINE image_division
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
  END MODULE parallelism
  !-----------------------------------------------------------------------
