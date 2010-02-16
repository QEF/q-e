!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine hpsort_eps (n, ra, ind, eps)
  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm,
  ! and considering two elements being equal if their values differ
  ! for less than "eps".
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  use kinds
  implicit none  
  !-input/output variables
  integer, intent(in) :: n  
  integer, intent(inout) :: ind (*)  
  real(DP), intent(inout) :: ra (*)
  real(DP), intent(in) :: eps
  !-local variables
  integer :: i, ir, j, l, iind  
  real(DP) :: rra  
  ! initialize index array
  if (ind (1) .eq.0) then  
     do i = 1, n  
        ind (i) = i  
     enddo
  endif
  ! nothing to order
  if (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  

  ir = n  

  sorting: do 
  
    ! still in hiring phase
    if ( l .gt. 1 ) then  
       l    = l - 1  
       rra  = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    else  
       ! clear a space at the end of the array
       rra  = ra (ir)  
       !
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       !
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       if ( ir .eq. 1 ) then  
          ! the least competent worker at all !
          ra (1)  = rra  
          !
          ind (1) = iind  
          exit sorting  
       endif
    endif
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    do while ( j .le. ir )  
       if ( j .lt. ir ) then  
          ! compare to better underling
          if ( abs(ra(j)-ra(j+1)).ge.eps ) then  
             if (ra(j).lt.ra(j+1)) j = j + 1
          else
             ! this means ra(j) == ra(j+1) within tolerance
             if (ind (j) .lt.ind (j + 1) ) j = j + 1
          endif
       endif
       ! demote rra
       if ( abs(rra - ra(j)).ge.eps ) then  
          if (rra.lt.ra(j)) then
             ra (i) = ra (j)  
             ind (i) = ind (j)  
             i = j  
             j = j + j  
          else
             ! set j to terminate do-while loop
             j = ir + 1  
          end if
       else
          !this means rra == ra(j) within tolerance
          ! demote rra
          if (iind.lt.ind (j) ) then
             ra (i) = ra (j)
             ind (i) = ind (j)
             i = j
             j = j + j
          else
             ! set j to terminate do-while loop
             j = ir + 1
          endif
       end if
    enddo
    ra (i) = rra  
    ind (i) = iind  

  end do sorting    
  !
end subroutine hpsort_eps

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine hpsort (n, ra, ind)  
  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm.
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  use kinds
  implicit none  
  !-input/output variables
  integer :: n  
  integer :: ind (*)  
  real(DP) :: ra (*)  
  !-local variables
  integer :: i, ir, j, l, iind  
  real(DP) :: rra  
  ! initialize index array
  if (ind (1) .eq.0) then  
     do i = 1, n  
        ind (i) = i  
     enddo
  endif
  ! nothing to order
  if (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  
  ir = n  
10 continue  
  ! still in hiring phase
  if (l.gt.1) then  
     l = l - 1  
     rra = ra (l)  
     iind = ind (l)  
     ! in retirement-promotion phase.
  else  
     ! clear a space at the end of the array
     rra = ra (ir)  
     !
     iind = ind (ir)  
     ! retire the top of the heap into it
     ra (ir) = ra (1)  
     !
     ind (ir) = ind (1)  
     ! decrease the size of the corporation
     ir = ir - 1  
     ! done with the last promotion
     if (ir.eq.1) then  
        ! the least competent worker at all !
        ra (1) = rra  
        !
        ind (1) = iind  
        return  
     endif
  endif
  ! wheter in hiring or promotion phase, we
  i = l  
  ! set up to place rra in its proper level
  j = l + l  
  !
  do while (j.le.ir)  
     if (j.lt.ir) then  
        ! compare to better underling
        if (ra (j) .lt.ra (j + 1) ) then  
           j = j + 1  
        elseif (ra (j) .eq.ra (j + 1) ) then  
           if (ind (j) .lt.ind (j + 1) ) j = j + 1  
        endif
     endif
     ! demote rra
     if (rra.lt.ra (j) ) then  
        ra (i) = ra (j)  
        ind (i) = ind (j)  
        i = j  
        j = j + j  
     elseif (rra.eq.ra (j) ) then  
        ! demote rra
        if (iind.lt.ind (j) ) then  
           ra (i) = ra (j)  
           ind (i) = ind (j)  
           i = j  
           j = j + j  
        else  
           ! set j to terminate do-while loop
           j = ir + 1  
        endif
        ! this is the right place for rra
     else  
        ! set j to terminate do-while loop
        j = ir + 1  
     endif
  enddo
  ra (i) = rra  
  ind (i) = iind  
  goto 10  
  !
end subroutine hpsort


!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine ihpsort (n, ia, ind)  
  !---------------------------------------------------------------------
  ! sort an integer array ia(1:n) into ascending order using heapsort algorithm.
  ! n is input, ia is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ia).
  ! in case of equal values in the data array (ia) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  implicit none  
  !-input/output variables
  integer :: n  
  integer :: ind (*)  
  integer :: ia (*)  
  !-local variables
  integer :: i, ir, j, l, iind  
  integer :: iia  
  ! initialize index array
  if (ind (1) .eq.0) then  
     do i = 1, n  
        ind (i) = i  
     enddo
  endif
  ! nothing to order
  if (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  
  ir = n  
10 continue  
  ! still in hiring phase
  if (l.gt.1) then  
     l = l - 1  
     iia = ia (l)  
     iind = ind (l)  
     ! in retirement-promotion phase.
  else  
     ! clear a space at the end of the array
     iia = ia (ir)  
     !
     iind = ind (ir)  
     ! retire the top of the heap into it
     ia (ir) = ia (1)  
     !
     ind (ir) = ind (1)  
     ! decrease the size of the corporation
     ir = ir - 1  
     ! done with the last promotion
     if (ir.eq.1) then  
        ! the least competent worker at all !
        ia (1) = iia  
        !
        ind (1) = iind  
        return  
     endif
  endif
  ! wheter in hiring or promotion phase, we
  i = l  
  ! set up to place iia in its proper level
  j = l + l  
  !
  do while (j.le.ir)  
     if (j.lt.ir) then  
        ! compare to better underling
        if (ia (j) .lt.ia (j + 1) ) then  
           j = j + 1  
        elseif (ia (j) .eq.ia (j + 1) ) then  
           if (ind (j) .lt.ind (j + 1) ) j = j + 1  
        endif
     endif
     ! demote iia
     if (iia.lt.ia (j) ) then  
        ia (i) = ia (j)  
        ind (i) = ind (j)  
        i = j  
        j = j + j  
     elseif (iia.eq.ia (j) ) then  
        ! demote iia
        if (iind.lt.ind (j) ) then  
           ia (i) = ia (j)  
           ind (i) = ind (j)  
           i = j  
           j = j + j  
        else  
           ! set j to terminate do-while loop
           j = ir + 1  
        endif
        ! this is the right place for iia
     else  
        ! set j to terminate do-while loop
        j = ir + 1  
     endif
  enddo
  ia (i) = iia  
  ind (i) = iind  
  goto 10  
  !
end subroutine ihpsort
!-------------------------------------------------------------------------
      subroutine kb07ad_cp90(count,n,idx)
!-------------------------------------------------------------------------
!     
!             kb07ad      handles double precision variables
!  standard fortran 66 (a verified pfort subroutine)
!  the work-space 'mark' of length 50 permits up to 2**(50/2) numbers
!  to be sorted.
      implicit none
      integer :: n, idx(*)
      real(8) :: count(*)
      real(8) :: av, x
      integer :: k1, ifk, lngth, ip, k, it, ifka, intest, iy
      integer :: i, m, la, is, idf, mloop, is1, j, mark(50)
!  set index array to original order .
      do i=1,n
         idx(i)=i
      end do
!  check that a trivial case has not been entered .
      if(n.eq.1) go to 10
      if(n.gt.1) go to 30
      write(6,20)
   20 format(///20x,'***kb07ad***no numbers to be sorted ** return to', &
     & ' calling program' ) 
      goto 10
!  'm' is the length of segment which is short enough to enter
!  the final sorting routine. it may be easily changed.
   30 m=12
!  set up initial values.
      la=2
      is=1
      idf=n
      do 190 mloop=1,n
!  if segment is short enough sort with final sorting routine .
      ifka=idf-is
      if((ifka+1).gt.m)goto 70
!********* final sorting ***
!  ( a simple bubble sort )
      is1=is+1
      do 60 j=is1,idf
      i=j
   40 if(count(i-1).lt.count(i))goto 60
      if(count(i-1).gt.count(i))goto 50
      if(idx(i-1).lt.idx(i))goto 60
   50 av=count(i-1)
      count(i-1)=count(i)
      count(i)=av
      it=idx(i-1)
      idx(i-1)=idx(i)
      idx(i)=it
      i=i-1
      if(i.gt.is)goto  40
   60 continue
      la=la-2
      goto 170
!             *******  quicksort  ********
!  select the number in the central position in the segment as
!  the test number.replace it with the number from the segment's
!  highest address.
   70 iy=(is+idf)/2
      x=count(iy)
      intest=idx(iy)
      count(iy)=count(idf)
      idx(iy)=idx(idf)
!  the markers 'i' and 'ifk' are used for the beginning and end
!  of the section not so far tested against the present value
!  of x .
      k=1
      ifk=idf
!  we alternate between the outer loop that increases i and the
!  inner loop that reduces ifk, moving numbers and indices as
!  necessary, until they meet .
      do 110 i=is,idf
      if(x.gt.count(i))goto 110
      if(x.lt.count(i))goto 80
      if(intest.gt.idx(i))goto 110
   80 if(i.ge.ifk)goto 120
      count(ifk)=count(i)
      idx(ifk)=idx(i)
      k1=k
      do 100 k=k1,ifka
      ifk=idf-k
      if(count(ifk).gt.x)goto 100
      if(count(ifk).lt.x)goto 90
      if(intest.le.idx(ifk))goto 100
   90 if(i.ge.ifk)goto 130
      count(i)=count(ifk)
      idx(i)=idx(ifk)
      go to 110
  100 continue
      goto 120
  110 continue
!  return the test number to the position marked by the marker
!  which did not move last. it divides the initial segment into
!  2 parts. any element in the first part is less than or equal
!  to any element in the second part, and they may now be sorted
!  independently .
  120 count(ifk)=x
      idx(ifk)=intest
      ip=ifk
      goto 140
  130 count(i)=x
      idx(i)=intest
      ip=i
!  store the longer subdivision in workspace.
  140 if((ip-is).gt.(idf-ip))goto 150
      mark(la)=idf
      mark(la-1)=ip+1
      idf=ip-1
      goto 160
  150 mark(la)=ip-1
      mark(la-1)=is
      is=ip+1
!  find the length of the shorter subdivision.
  160 lngth=idf-is
      if(lngth.le.0)goto 180
!  if it contains more than one element supply it with workspace .
      la=la+2
      goto 190
  170 if(la.le.0)goto 10
!  obtain the address of the shortest segment awaiting quicksort
  180 idf=mark(la)
      is=mark(la-1)
  190 continue
   10 return
      end subroutine kb07ad_cp90

