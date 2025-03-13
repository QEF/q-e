!
! Copyright (C) 2001-2022 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
subroutine hpsort_eps (n, ra, ind, eps)
  !---------------------------------------------------------------------
  !! Sort an array ra(1:n) into ascending order using heapsort algorithm,
  !! and considering two elements being equal if their values differ
  !! for less than "eps". IMPORTANT NOTICE (PG February 2022):
  !! Assume you have in input a,b,c with c < b < a and a-b < eps, b-c < eps,
  !! but a-c > eps. The resulting output order should be c,a,b, but may turn
  !! out to be a,b,c instead. I think this is a bug. I don't know how to fix
  !! it, I am not sure it does any harm, but re-ordering k+G with this same
  !! routine may yield a different ordering for k+G and G vectors even if k=0.
  !! This is a bug that has been around for years. The current work-around 
  !! (in routine gk_sort) is to avoid recomputing indices for k+G if k=0.
  !!
  !! \(\text{n}\) is input, \(\text{ra}\) is replaced on output by its 
  !! sorted rearrangement.  
  !! Create an index table (ind) by making an exchange in the index array
  !! whenever an exchange is made on the sorted data array (\(\text{ra}\)).  
  !! In case of equal values in the data array (\(\text{ra}\)) the values
  !! in the index array (ind) are used to order the entries.  
  !! If on input ind(1) = 0 then indices are initialized in the routine,
  !! if on input ind(1) != 0 then indices are assumed to have been
  !! initialized before entering the routine and these indices are carried
  !! around during the sorting process.
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  use kinds, only : DP
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
!---------------------------------------------------------------------
subroutine hpsort (n, ra, ind)  
  !---------------------------------------------------------------------
  !! Sort an array ra(1:n) into ascending order using heapsort algorithm.
  !! As hpsort_eps, without the "eps" trick (or equivalently, eps=0)
  !
  use kinds, only : DP
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
!---------------------------------------------------------------------
subroutine ihpsort (n, ia, ind)  
  !---------------------------------------------------------------------
  !! As "hpsort", for integer array ia(1:n)
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
