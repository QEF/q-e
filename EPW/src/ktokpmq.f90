  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------
  subroutine ktokpmq ( xk, xq, sign, ipool, nkq, nkq_abs)
  !--------------------------------------------------------
  !!
  !!   For a given k point in cart coord, find the index 
  !!   of the corresponding (k + sign*q) point
  !!
  !!   In the parallel case, determine also the pool number
  !!   nkq is the in-pool index, nkq_abs is the absolute
  !!   index
  !!
  !--------------------------------------------------------
  !
  USE kinds,          only : DP
  use pwcom,          only : nkstot
  USE cell_base,      ONLY : at
  USE start_k,        ONLY : nk1, nk2, nk3
  use epwcom,         only : xk_cryst
  USE mp_global,      only : nproc_pool, npool
  USE mp_images,      ONLY : nproc_image
  USE mp,             only : mp_barrier, mp_bcast
  implicit none
  !
  INTEGER, INTENT (in) :: sign
  !! +1 for searching k+q, -1 for k-q
  INTEGER, INTENT (out) :: nkq
  !! in the parallel case, the pool hosting the k+-q point    
  INTEGER, INTENT (out) :: nkq_abs
  !! the index of k+sign*q
  REAL(kind=DP), INTENT (in) :: xk(3)
  !! coordinates of k points and q points
  REAL(kind=DP), INTENT (in) :: xq(3)
  !! Coordinates of k+q point

  !
  ! work variables
  !
  INTEGER :: kunit
  !! the absolute index of k+sign*q (in the full k grid)
  real(kind=DP) :: xxk (3), xxq (3)
  integer ::  n,  ik, ipool
  real(kind=DP) :: xx, yy, zz, xx_c, yy_c, zz_c, eps
  logical :: in_the_list, found
  !
  integer :: iks, nkl, nkr, jpool
  !
  kunit =1 
  ! loosy tolerance, no problem since we use integer comparisons
  eps = 1.d-5
  IF (abs(sign).ne.1) call errore('ktokpmq','sign must be +1 or -1',1)
  !
  ! bring k and q in crystal coordinates
  !
  xxk = xk
  xxq = xq
  !
  CALL cryst_to_cart (1, xxk, at, -1)
  CALL cryst_to_cart (1, xxq, at, -1)
  !
  !  check that k is actually on a uniform mesh centered at gamma
  !
  xx = xxk(1)*nk1
  yy = xxk(2)*nk2
  zz = xxk(3)*nk3
  in_the_list = abs(xx-nint(xx)).le.eps .and. &
                abs(yy-nint(yy)).le.eps .and. &
                abs(zz-nint(zz)).le.eps
  IF (.not.in_the_list) call errore('ktokpmq','is this a uniform k-mesh?',1)
  !
  IF ( xx .lt. -eps .or. yy .lt. -eps .or. zz .lt. -eps ) &
     call errore('ktokpmq','coarse k-mesh needs to be strictly positive in 1st BZ',1)
  !
  !  now add the phonon wavevector and check that k+q falls again on the k grid
  !
  xxk = xxk + dble(sign) * xxq
  !
  xx = xxk(1)*nk1
  yy = xxk(2)*nk2
  zz = xxk(3)*nk3
  in_the_list = abs(xx-nint(xx)).le.eps .and. &
                abs(yy-nint(yy)).le.eps .and. &
                abs(zz-nint(zz)).le.eps
  IF (.not.in_the_list) call errore('ktokpmq','k+q does not fall on k-grid',1)
  !
  !  find the index of this k+q in the k-grid
  !
  !  make sure xx, yy and zz are in 1st BZ
  !
  CALL backtoBZ( xx, yy, zz, nk1, nk2, nk3 )
  !
  n = 0
  found = .false.
  DO ik = 1, nkstot
     xx_c = xk_cryst(1,ik)*nk1
     yy_c = xk_cryst(2,ik)*nk2
     zz_c = xk_cryst(3,ik)*nk3
     !
     ! check that the k-mesh was defined in the positive region of 1st BZ
     !
     IF ( xx_c .lt. -eps .or. yy_c .lt. -eps .or. zz_c .lt. -eps ) &
        call errore('ktokpmq','coarse k-mesh needs to be strictly positive in 1st BZ',1)
     !
     found = nint(xx_c) .eq. nint(xx) .and. &
             nint(yy_c) .eq. nint(yy) .and. &
             nint(zz_c) .eq. nint(zz)
     IF (found) THEN  
        n = ik
        EXIT
     ENDIF
  ENDDO
  !
  !  26/06/2012 RM
  !  since coarse k- and q- meshes are commensurate, one can easily find n
  !
!  n = nint(xx) * nk2 * nk3 + nint(yy) * nk3 + nint(zz) + 1
  !
  IF (n .eq. 0) call errore('ktokpmq','problem indexing k+q',1)
  !
  !  Now n represents the index of k+sign*q in the original k grid.
  !  In the parallel case we have to find the corresponding pool 
  !  and index in the pool
  !
#if defined(__MPI)
  !
  npool = nproc_image/nproc_pool
  !
  DO jpool = 0, npool-1
    !
    nkl   = kunit * ( nkstot / npool )
    nkr   = ( nkstot - nkl * npool ) / kunit
    !
    !  the reminder goes to the first nkr pools (0...nkr-1)
    !
    IF ( jpool < nkr ) nkl = nkl + kunit
    !
    !  the index of the first k point in this pool
    !
    iks = nkl * jpool + 1
    IF ( jpool >= nkr ) iks = iks + nkr * kunit
    !
    IF (n.ge.iks) then
      ipool = jpool+1
      nkq = n - iks + 1
    ENDIF
    !
  ENDDO
  !
#else
  !
  nkq = n
  !
#endif
  !
  nkq_abs = n
  !
  end subroutine ktokpmq

!---------------------------------
subroutine ckbounds(lower, upper)
!---------------------------------
!!
!!   Subroutine finds the lower and upper
!!   bounds of the coarse k-grid in parallel
!!
!---------------------------------
  !
  use pwcom,         ONLY : nkstot
  use mp_global,     ONLY : my_pool_id, npool
  !
  implicit none
  integer, intent(out):: lower, upper
  !
#if defined(__MPI)
  !
  integer :: nkst, nkstott, rest
  !
  nkstott = nkstot
  !
  nkst = ( nkstott/ npool )
  rest = ( nkstott - nkst * npool ) 
 IF (my_pool_id < rest ) then
    nkst=nkst+1
    lower = my_pool_id*nkst + 1
    upper = lower + nkst - 1
 ELSE
    lower = rest*(nkst+1)+(my_pool_id-rest)*nkst + 1
    upper = lower + nkst - 1
 ENDIF
#else
  lower = 1
  upper = nkstot
#endif

end subroutine

!---------------------------------
subroutine para_bounds(lower, upper, total)
!---------------------------------
!!
!!   Subroutine finds the lower and upper
!!   bounds if we split some quantity over pools
!!
!---------------------------------
  !
  use mp_global,   ONLY : my_pool_id, npool 
  !
  implicit none
  integer, intent(out):: lower, upper
  integer, intent(in):: total
  !
#if defined(__MPI)
  !  
  integer :: rest, nrst
  !
  IF (total .le. npool) THEN
     IF (my_pool_id .lt. total) THEN
        lower = my_pool_id + 1
        upper = lower
     ELSE
        lower = 1
        upper  = 0
     ENDIF
  ELSE
     nrst = total / npool
     rest = total - nrst * npool
     IF (my_pool_id < rest ) then
        nrst=nrst+1
        lower = my_pool_id*nrst + 1
        upper = lower + nrst - 1
     ELSE
        lower = rest*(nrst+1)+(my_pool_id-rest)*nrst + 1
        upper = lower + nrst - 1
     ENDIF
  ENDIF
#else
  lower = 1
  upper = total
#endif
end subroutine

!---------------------------------
subroutine backtoBZ ( xx, yy, zz, n1, n2, n3 )
!---------------------------------
  !!
  !!  Brings xx, yy, and zz  into first BZ 
  !!
  !---------------------------------
  !
  USE kinds,  ONLY : DP
  !
  implicit none
  integer, intent(in) :: n1, n2, n3
  real(DP), intent(inout) :: xx, yy, zz
  integer :: ib
  !
  ! more translations are needed to go back to the first BZ when the unit cell
  ! is far from cubic
  !
  DO ib = -2,0
     IF (nint(xx) .lt. ib*n1) xx = xx + (-ib+1)*n1
     IF (nint(yy) .lt. ib*n2) yy = yy + (-ib+1)*n2
     IF (nint(zz) .lt. ib*n3) zz = zz + (-ib+1)*n3
  ENDDO
  DO ib = 2,1,-1
     IF (nint(xx) .ge. ib*n1) xx = xx - ib*n1
     IF (nint(yy) .ge. ib*n2) yy = yy - ib*n2
     IF (nint(zz) .ge. ib*n3) zz = zz - ib*n3
  ENDDO
  !
end subroutine

