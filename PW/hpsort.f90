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
  use parameters
  implicit none
  !-input/output variables
  integer :: n
  integer :: ind (n)
  real(kind=DP) :: ra (n)
  !-local variables
  integer :: i, ir, j, l, iind
  real(kind=DP) :: rra
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
