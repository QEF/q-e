!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#ifdef T3D

#ifdef AIX
#define CGERV2D zgerv2d
#define CGESD2D zgesd2d
#define CGEBR2D zgebr2d
#define CGEBS2D zgebs2d
#endif

!-----------------------------------------------------------------------

subroutine gridsetup_local (nproc, nprow, npcol)  
  !-----------------------------------------------------------------------
  !
  ! This subroutine factorizes the number of processors (NPROC)
  ! into NPROW and NPCOL,  that are the sizes of the 2D processors mesh.
  !
  ! Written by Carlo Cavazzoni
  !
  implicit none  

  integer :: nproc, nprow, npcol  
  ! input: number of processors
  ! output: number of rows
  ! output: number of column

  integer :: sqrtnp, i  
  ! the maximum size to check
  ! counter
  sqrtnp = int (sqrt (dble (nproc) ) + 1)  
  do i = 1, sqrtnp  
     if (mod (nproc, i) .eq.0) nprow = i  
  enddo
  npcol = nproc / nprow  
  return  

end subroutine gridsetup_local
subroutine blockset_priv (nb, nbuser, n, nprow, npcol)  
  !
  !     This subroutine try to choose an optimal block size
  !     for the distributed matrix.
  !
  !     Written by Carlo Cavazzoni
  !
  implicit none  
  integer :: nb, n, nprow, npcol, nbuser  
  nb = min (n / nprow, n / npcol)  
  if (nbuser.gt.0) then  
     nb = min (nb, nbuser)  
  endif
  nb = min (nb, n)  
  if (n.lt.10) nb = 1  
  return  
end subroutine blockset_priv
!
!----------------------------------------------------------------------
!
subroutine eigen (n, aout, ldaout, a, desca, work)  
  !
  implicit none  
  !     This routine accumulates the eigenvectors on the root
  !     processor and then broadcast them to the other processors
  !
  !     .. Scalar Arguments ..
  integer :: ia, icprnt, irprnt, ja, m, n, nout, ldaout  
  !     ..
  !     .. Array Arguments ..
  integer :: desca ( * )  
  complex (kind=8) :: aout (ldaout, * ), a ( * ), work ( * )  
  !     ..
  !     .. Parameters ..
  integer :: block_cyclic_2d, csrc_, ctxt_, dlen_, dtype_, lld_, &
       mb_, m_, nb_, n_, rsrc_
  parameter (block_cyclic_2d = 1, dlen_ = 9, dtype_ = 1, ctxt_ = 2, &
       m_ = 3, n_ = 4, mb_ = 5, nb_ = 6, rsrc_ = 7, csrc_ = 8, lld_ = 9)
  !     ..
  !     .. Local Scalars ..
  integer :: h, i, iacol, iarow, ib, ictxt, icurcol, icurrow, ii, &
       iia, in, j, jb, jj, jja, jn, k, lda, mycol, myrow, npcol, nprow
  !     ..
  !     .. External Subroutines ..
  external blacs_barrier, blacs_gridinfo, infog2l, CGERV2D, CGESD2D  
  !     ..
  !     .. External Functions ..
  integer :: iceil  
  external iceil  
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic min  
  !     ..
  !     .. Executable Statements ..
  !
  !     Get grid parameters
  !
  icprnt = 0  
  irprnt = 0  
  ia = 1  
  ja = 1  

  m = n  

  call setv (2 * ldaout * n, 0.d0, aout, 1)  
  ictxt = desca (ctxt_)  
  call blacs_gridinfo (ictxt, nprow, npcol, myrow, mycol)  
  !
  call infog2l (ia, ja, desca, nprow, npcol, myrow, mycol, iia, jja, &
       iarow, iacol)
  icurrow = iarow  
  icurcol = iacol  
  ii = iia  
  jj = jja  
  lda = desca (lld_)  
  !
  !     Handle the first block of column separately
  !
  jn = min (iceil (ja, desca (nb_) ) * desca (nb_), ja + n - 1)  
  jb = jn - ja + 1  
  do h = 0, jb - 1  
     in = min (iceil (ia, desca (mb_) ) * desca (mb_), ia + m - 1)  
     ib = in - ia + 1  
     if (icurrow.eq.irprnt.and.icurcol.eq.icprnt) then  
        if (myrow.eq.irprnt.and.mycol.eq.icprnt) then  
           do k = 0, ib - 1  
              aout (ia + k, ja + h) = a (ii + k + (jj + h - 1) * lda)  
           enddo
        endif
     else  
        if (myrow.eq.icurrow.and.mycol.eq.icurcol) then  
           call CGESD2D (ictxt, ib, 1, a (ii + (jj + h - 1) * lda), &
                lda, irprnt, icprnt)
        elseif (myrow.eq.irprnt.and.mycol.eq.icprnt) then  
           call CGERV2D (ictxt, ib, 1, work, desca (mb_), icurrow, &
                icurcol)
           do k = 1, ib  
              aout (ia + k - 1, ja + h) = work (k)  
           enddo
        endif
     endif
     if (myrow.eq.icurrow) ii = ii + ib  
     icurrow = mod (icurrow + 1, nprow)  
     call blacs_barrier (ictxt, 'All')  
     !
     !        Loop over remaining block of rows
     !
     do i = in + 1, ia + m - 1, desca (mb_)  
        ib = min (desca (mb_), ia + m - i)  
        if (icurrow.eq.irprnt.and.icurcol.eq.icprnt) then  
           if (myrow.eq.irprnt.and.mycol.eq.icprnt) then  
              do k = 0, ib - 1  
                 aout (i + k, ja + h) = a (ii + k + (jj + h - 1) * lda)  
              enddo
           endif
        else  
           if (myrow.eq.icurrow.and.mycol.eq.icurcol) then  
              call CGESD2D (ictxt, ib, 1, a (ii + (jj + h - 1) * lda), &
                   lda, irprnt, icprnt)
           elseif (myrow.eq.irprnt.and.mycol.eq.icprnt) then  
              call CGERV2D (ictxt, ib, 1, work, desca (mb_), icurrow, &
                   icurcol)
              do k = 1, ib  
                 aout (i + k - 1, ja + h) = work (k)  
              enddo
           endif
        endif
        if (myrow.eq.icurrow) ii = ii + ib  
        icurrow = mod (icurrow + 1, nprow)  
        call blacs_barrier (ictxt, 'All')  
     enddo
     !
     ii = iia  
     icurrow = iarow  

  enddo
  !
  if (mycol.eq.icurcol) jj = jj + jb  
  icurcol = mod (icurcol + 1, npcol)  
  call blacs_barrier (ictxt, 'All')  
  !
  !     Loop over remaining column blocks
  !
  do j = jn + 1, ja + n - 1, desca (nb_)  
     jb = min (desca (nb_), ja + n - j)  
     do h = 0, jb - 1  
        in = min (iceil (ia, desca (mb_) ) * desca (mb_), ia + m - 1)  
        ib = in - ia + 1  
        if (icurrow.eq.irprnt.and.icurcol.eq.icprnt) then  
           if (myrow.eq.irprnt.and.mycol.eq.icprnt) then  
              do k = 0, ib - 1  
                 aout (ia + k, j + h) = a (ii + k + (jj + h - 1) * lda)  
              enddo
           endif
        else  
           if (myrow.eq.icurrow.and.mycol.eq.icurcol) then  
              call CGESD2D (ictxt, ib, 1, a (ii + (jj + h - 1) * lda), &
                   lda, irprnt, icprnt)
           elseif (myrow.eq.irprnt.and.mycol.eq.icprnt) then  
              call CGERV2D (ictxt, ib, 1, work, desca (mb_), icurrow, &
                   icurcol)
              do k = 1, ib  
                 aout (ia + k - 1, j + h) = work (k)  
              enddo
           endif
        endif
        if (myrow.eq.icurrow) ii = ii + ib  
        icurrow = mod (icurrow + 1, nprow)  
        call blacs_barrier (ictxt, 'All')  
        !
        !           Loop over remaining block of rows
        !
        do i = in + 1, ia + m - 1, desca (mb_)  
           ib = min (desca (mb_), ia + m - i)  
           if (icurrow.eq.irprnt.and.icurcol.eq.icprnt) then  
              if (myrow.eq.irprnt.and.mycol.eq.icprnt) then  
                 do k = 0, ib - 1  
                    aout (i + k, j + h) = a (ii + k + (jj + h - 1) * lda)  
                 enddo
              endif
           else  
              if (myrow.eq.icurrow.and.mycol.eq.icurcol) then  
                 call CGESD2D (ictxt, ib, 1, a (ii + (jj + h - 1) * lda), &
                      lda, irprnt, icprnt)
              elseif (myrow.eq.irprnt.and.mycol.eq.icprnt) then  
                 call CGERV2D (ictxt, ib, 1, work, desca (mb_), icurrow, &
                      icurcol)
                 do k = 1, ib  
                    aout (i + k - 1, j + h) = work (k)  
                 enddo
              endif
           endif
           if (myrow.eq.icurrow) ii = ii + ib  
           icurrow = mod (icurrow + 1, nprow)  
           call blacs_barrier (ictxt, 'All')  
        enddo
        !
        ii = iia  
        icurrow = iarow  
     enddo
     !
     if (mycol.eq.icurcol) jj = jj + jb  
     icurcol = mod (icurcol + 1, npcol)  
     call blacs_barrier (ictxt, 'All')  
     !
  enddo
  !
  if ( (myrow.eq.0) .and. (mycol.eq.0) ) then  
     call CGEBS2D (ictxt, 'All', 'h', n, n, aout, ldaout)  
  else  
     call CGEBR2D (ictxt, 'All', 'h', n, n, aout, ldaout, 0, 0)  
  endif
  !
  call blacs_barrier (ictxt, 'All')  
  !
  return  
  !
end subroutine eigen
#else
subroutine scaladummy  
  return  
end subroutine scaladummy
#endif
