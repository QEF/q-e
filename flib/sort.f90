!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      logical function cpgt(a,b)
      USE kinds
      USE constants, only: eps8
      implicit none
      REAL(DP) :: a, b, r
      r = abs(a-b)
      if( r .lt.  eps8 ) then
        cpgt = .false.
      else
        cpgt = ( a .gt. b )
      end if
      end function cpgt

      logical function cplt(a,b)
      USE kinds
      USE constants, only: eps8
      implicit none
      REAL(DP) :: a, b, r
      r = abs(a-b)
      if( r .lt.  eps8 ) then
        cplt = .false.
      else
        cplt = ( a .lt. b )
      end if
      end function cplt
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
          if ( hslt( ra (j),  ra (j + 1) ) ) then  
             j = j + 1  
          else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
             ! this means ra(j) == ra(j+1) within tolerance
             if (ind (j) .lt.ind (j + 1) ) j = j + 1
          endif
       endif
       ! demote rra
       if ( hslt( rra, ra (j) ) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
       else if ( .not. hslt ( ra(j) , rra ) ) then
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
          ! this is the right place for rra
       else
          ! set j to terminate do-while loop
          j = ir + 1  
       endif
    enddo
    ra (i) = rra  
    ind (i) = iind  

  end do sorting    

contains 

  !  internal function 
  !  compare two real number and return the result

  logical function hslt( a, b )
    REAL(DP) :: a, b
    if( abs(a-b) <  eps ) then
      hslt = .false.
    else
      hslt = ( a < b )
    end if
  end function hslt

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

!     ==================================================================
      SUBROUTINE gqsort(COUNT,N,IDX)
!     ==--------------------------------------------------------------==
!     == Sorting routine for the reciprocal space vectors (g)         ==
!     == Warning, this is not an exact SORT!! This routine has been   ==
!     == designed to give always the same order for the G vectors of  ==
!     == a given shell, independently of the processor                ==
!     == THE WORK-SPACE 'MARK' OF LENGTH 50 PERMITS UP TO 2**(50/2)   ==
!     ==--------------------------------------------------------------==
      USE kinds
      
      INTEGER :: N, MARK, I, M, LA, IS, IF, MLOOP, IFKA, IS1, J, INT, &
                 IY, INTEST, K, IFK, K1, IP, LNGTH
      
      logical :: cpgt,cplt
      REAL(DP) :: COUNT(*),AV,X
      
      INTEGER :: IDX(*)
      DIMENSION :: MARK(50)
!     ==--------------------------------------------------------------==
!     ==  SET INDEX ARRAY TO ORIGINAL ORDER .                         ==
!     ==--------------------------------------------------------------==
      DO I=1,N
        IDX(I)=I
      ENDDO
!     ==--------------------------------------------------------------==
!     == CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED.              ==
!     ==--------------------------------------------------------------==
      IF(N.EQ.1)GOTO 200
      IF(N.GE.1)GOTO 30
      GOTO 200
!     ==--------------------------------------------------------------==
!     == 'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER  ==
!     == THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.         ==
!     ==--------------------------------------------------------------==
   30 M=12
!     ==--------------------------------------------------------------==
!     == SET UP INITIAL VALUES.                                       ==
!     ==--------------------------------------------------------------==
      LA=2
      IS=1
      IF=N
      DO 190 MLOOP=1,N
!     ==--------------------------------------------------------------==
!     ==  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE. ==
!     ==--------------------------------------------------------------==
        IFKA=IF-IS
        IF((IFKA+1).GT.M)GOTO 70
!     ==--------------------------------------------------------------==
!     == FINAL SORTING  ( A SIMPLE BUBBLE SORT )                      ==
!     ==--------------------------------------------------------------==
        IS1=IS+1
        DO 60 J=IS1,IF
          I=J
   40     IF(cplt(COUNT(I-1),COUNT(I)) )GOTO 60
          IF(cpgt(COUNT(I-1),COUNT(I)) )GOTO 50
          IF(IDX(I-1).LT.IDX(I))GOTO 60
   50     AV=COUNT(I-1)
          COUNT(I-1)=COUNT(I)
          COUNT(I)=AV
          INT=IDX(I-1)
          IDX(I-1)=IDX(I)
          IDX(I)=INT
          I=I-1
          IF(I.GT.IS)GOTO  40
   60   CONTINUE
        LA=LA-2
        GOTO 170
!     ==--------------------------------------------------------------==
!     ==                *******  QUICKSORT  ********                  ==
!     == SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS  ==
!     == THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENT'S==
!     == HIGHEST ADDRESS.                                             ==
!     ==--------------------------------------------------------------==
   70   IY=(IS+IF)/2
        X=COUNT(IY)
        INTEST=IDX(IY)
        COUNT(IY)=COUNT(IF)
        IDX(IY)=IDX(IF)
!     ==--------------------------------------------------------------==
!     == THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END ==
!     == OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE   ==
!     == OF X .                                                       ==
!     ==--------------------------------------------------------------==
        K=1
        IFK=IF
!     ==--------------------------------------------------------------==
!     == WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE ==
!     == INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AND INDICES AS   ==
!     == NECESSARY, UNTIL THEY MEET .                                 ==
!     ==--------------------------------------------------------------==
        DO 110 I=IS,IF
          IF(cpgt(X,COUNT(I)))GOTO 110
          IF(cplt(X,COUNT(I)))GOTO 80
          IF(INTEST.GT.IDX(I))GOTO 110
   80     IF(I.GE.IFK)GOTO 120
          COUNT(IFK)=COUNT(I)
          IDX(IFK)=IDX(I)
          K1=K
          DO 100 K=K1,IFKA
            IFK=IF-K
            IF(cpgt(COUNT(IFK),X))GOTO 100
            IF(cplt(COUNT(IFK),X))GOTO 90
            IF(INTEST.LE.IDX(IFK))GOTO 100
   90       IF(I.GE.IFK)GOTO 130
            COUNT(I)=COUNT(IFK)
            IDX(I)=IDX(IFK)
            GO TO 110
  100     CONTINUE
          GOTO 120
  110   CONTINUE
!     ==--------------------------------------------------------------==
!     == RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER  ==
!     == WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO ==
!     == 2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN OR EQUAL ==
!     == TO ANY ELEMENT IN THE SECOND PART, AND THEY MAY NOW BE SORTED==
!     == INDEPENDENTLY .                                              ==
!     ==--------------------------------------------------------------==
  120   COUNT(IFK)=X
        IDX(IFK)=INTEST
        IP=IFK
        GOTO 140
  130   COUNT(I)=X
        IDX(I)=INTEST
        IP=I
!     ==--------------------------------------------------------------==
!     ==  STORE THE LONGER SUBDIVISION IN WORKSPACE.                  ==
!     ==--------------------------------------------------------------==
  140   IF((IP-IS).GT.(IF-IP))GOTO 150
        MARK(LA)=IF
        MARK(LA-1)=IP+1
        IF=IP-1
        GOTO 160
  150   MARK(LA)=IP-1
        MARK(LA-1)=IS
        IS=IP+1
!     ==--------------------------------------------------------------==
!     == FIND THE LENGTH OF THE SHORTER SUBDIVISION.                  ==
!     ==--------------------------------------------------------------==
  160   LNGTH=IF-IS
        IF(LNGTH.LE.0)GOTO 180
!     ==--------------------------------------------------------------==
!     == IF IT CONTAINS MORE THAN ONE ELEMENT SUPPLY IT WITH WORKSPACE==
!     ==--------------------------------------------------------------==
        LA=LA+2
        GOTO 190
  170   IF(LA.LE.0)GOTO 200
!     ==--------------------------------------------------------------==
!     == OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT==
!     ==--------------------------------------------------------------==
  180   IF=MARK(LA)
        IS=MARK(LA-1)
  190 CONTINUE
!     ==--------------------------------------------------------------==
  200 RETURN
      END SUBROUTINE gqsort
!     ==================================================================


!     ==================================================================
      SUBROUTINE iqsort(COUNT,N,IDX)
!     ==--------------------------------------------------------------==
!     == same as rqsort but for array of integers                     ==
!     ==--------------------------------------------------------------==
      USE kinds
      
      INTEGER :: N, I, M, LA, IS, IF, MLOOP, IFKA, IS1, J, INT, &
                 IY, INTEST, K, IFK, K1, IP, LNGTH
            
      
      INTEGER :: COUNT(*),AV,X
      INTEGER :: IDX(*),MARK(50)
!     ==--------------------------------------------------------------==
!     ==  SET INDEX ARRAY TO ORIGINAL ORDER .                         ==
!     ==--------------------------------------------------------------==
      DO I=1,N
        IDX(I)=I
      ENDDO
!     ==--------------------------------------------------------------==
!     == CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED.              ==
!     ==--------------------------------------------------------------==
      IF(N.EQ.1)GOTO 200
      IF(N.GE.1)GOTO 30
      GOTO 200
!     ==--------------------------------------------------------------==
!     == 'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER  ==
!     == THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.         ==
!     ==--------------------------------------------------------------==
   30 M=12
!     ==--------------------------------------------------------------==
!     == SET UP INITIAL VALUES.                                       ==
!     ==--------------------------------------------------------------==
      LA=2
      IS=1
      IF=N
      DO 190 MLOOP=1,N
!     ==--------------------------------------------------------------==
!     ==  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE. ==
!     ==--------------------------------------------------------------==
        IFKA=IF-IS
        IF((IFKA+1).GT.M)GOTO 70
!     ==--------------------------------------------------------------==
!     == FINAL SORTING  ( A SIMPLE BUBBLE SORT )                      ==
!     ==--------------------------------------------------------------==
        IS1=IS+1
        DO 60 J=IS1,IF
          I=J
   40     IF((COUNT(I-1).LT.COUNT(I)) )GOTO 60
          IF((COUNT(I-1).GT.COUNT(I)) )GOTO 50
          IF(IDX(I-1).LT.IDX(I))GOTO 60
   50     AV=COUNT(I-1)
          COUNT(I-1)=COUNT(I)
          COUNT(I)=AV
          INT=IDX(I-1)
          IDX(I-1)=IDX(I)
          IDX(I)=INT
          I=I-1
          IF(I.GT.IS)GOTO  40
   60   CONTINUE
        LA=LA-2
        GOTO 170
!     ==--------------------------------------------------------------==
!     ==                *******  QUICKSORT  ********                  ==
!     == SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS  ==
!     == THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENT'S==
!     == HIGHEST ADDRESS.                                             ==
!     ==--------------------------------------------------------------==
   70   IY=(IS+IF)/2
        X=COUNT(IY)
        INTEST=IDX(IY)
        COUNT(IY)=COUNT(IF)
        IDX(IY)=IDX(IF)
!     ==--------------------------------------------------------------==
!     == THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END ==
!     == OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE   ==
!     == OF X .                                                       ==
!     ==--------------------------------------------------------------==
        K=1
        IFK=IF
!     ==--------------------------------------------------------------==
!     == WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE ==
!     == INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AND INDICES AS   ==
!     == NECESSARY, UNTIL THEY MEET .                                 ==
!     ==--------------------------------------------------------------==
        DO 110 I=IS,IF
          IF((X.GT.COUNT(I)))GOTO 110
          IF((X.LT.COUNT(I)))GOTO 80
          IF(INTEST.GT.IDX(I))GOTO 110
   80     IF(I.GE.IFK)GOTO 120
          COUNT(IFK)=COUNT(I)
          IDX(IFK)=IDX(I)
          K1=K
          DO 100 K=K1,IFKA
            IFK=IF-K
            IF((COUNT(IFK).GT.X))GOTO 100
            IF((COUNT(IFK).LT.X))GOTO 90
            IF(INTEST.LE.IDX(IFK))GOTO 100
   90       IF(I.GE.IFK)GOTO 130
            COUNT(I)=COUNT(IFK)
            IDX(I)=IDX(IFK)
            GO TO 110
  100     CONTINUE
          GOTO 120
  110   CONTINUE
!     ==--------------------------------------------------------------==
!     == RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER  ==
!     == WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO ==
!     == 2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN OR EQUAL ==
!     == TO ANY ELEMENT IN THE SECOND PART, AND THEY MAY NOW BE SORTED==
!     == INDEPENDENTLY .                                              ==
!     ==--------------------------------------------------------------==
  120   COUNT(IFK)=X
        IDX(IFK)=INTEST
        IP=IFK
        GOTO 140
  130   COUNT(I)=X
        IDX(I)=INTEST
        IP=I
!     ==--------------------------------------------------------------==
!     ==  STORE THE LONGER SUBDIVISION IN WORKSPACE.                  ==
!     ==--------------------------------------------------------------==
  140   IF((IP-IS).GT.(IF-IP))GOTO 150
        MARK(LA)=IF
        MARK(LA-1)=IP+1
        IF=IP-1
        GOTO 160
  150   MARK(LA)=IP-1
        MARK(LA-1)=IS
        IS=IP+1
!     ==--------------------------------------------------------------==
!     == FIND THE LENGTH OF THE SHORTER SUBDIVISION.                  ==
!     ==--------------------------------------------------------------==
  160   LNGTH=IF-IS
        IF(LNGTH.LE.0)GOTO 180
!     ==--------------------------------------------------------------==
!     == IF IT CONTAINS MORE THAN ONE ELEMENT SUPPLY IT WITH WORKSPACE==
!     ==--------------------------------------------------------------==
        LA=LA+2
        GOTO 190
  170   IF(LA.LE.0)GOTO 200
!     ==--------------------------------------------------------------==
!     == OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT==
!     ==--------------------------------------------------------------==
  180   IF=MARK(LA)
        IS=MARK(LA-1)
  190 CONTINUE
!     ==--------------------------------------------------------------==
  200 RETURN
      END SUBROUTINE iqsort
!     ==================================================================



!     ==================================================================
      SUBROUTINE rqsort(COUNT,N,IDX)
!     ==--------------------------------------------------------------==
!     == Sorting routine for the double precison arrayis              ==
!     == THE WORK-SPACE 'MARK' OF LENGTH 50 PERMITS UP TO 2**(50/2)   ==
!     ==--------------------------------------------------------------==
      USE kinds
      
      INTEGER :: N, I, M, LA, IS, IF, MLOOP, IFKA, IS1, J, INT, &
                 IY, INTEST, K, IFK, K1, IP, LNGTH
      
      REAL(DP) :: COUNT(*),AV,X
      INTEGER :: IDX(*), MARK(50)
!     ==--------------------------------------------------------------==
!     ==  SET INDEX ARRAY TO ORIGINAL ORDER .                         ==
!     ==--------------------------------------------------------------==
      DO I=1,N
        IDX(I)=I
      ENDDO
!     ==--------------------------------------------------------------==
!     == CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED.              ==
!     ==--------------------------------------------------------------==
      IF(N.EQ.1)GOTO 200
      IF(N.GE.1)GOTO 30
      GOTO 200
!     ==--------------------------------------------------------------==
!     == 'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER  ==
!     == THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.         ==
!     ==--------------------------------------------------------------==
   30 M=12
!     ==--------------------------------------------------------------==
!     == SET UP INITIAL VALUES.                                       ==
!     ==--------------------------------------------------------------==
      LA=2
      IS=1
      IF=N
      DO 190 MLOOP=1,N
!     ==--------------------------------------------------------------==
!     ==  IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE. ==
!     ==--------------------------------------------------------------==
        IFKA=IF-IS
        IF((IFKA+1).GT.M)GOTO 70
!     ==--------------------------------------------------------------==
!     == FINAL SORTING  ( A SIMPLE BUBBLE SORT )                      ==
!     ==--------------------------------------------------------------==
        IS1=IS+1
        DO 60 J=IS1,IF
          I=J
   40     IF( (COUNT(I-1) .LT. COUNT(I)) )GOTO 60
          IF( (COUNT(I-1) .GT. COUNT(I)) )GOTO 50
          IF(IDX(I-1).LT.IDX(I))GOTO 60
   50     AV=COUNT(I-1)
          COUNT(I-1)=COUNT(I)
          COUNT(I)=AV
          INT=IDX(I-1)
          IDX(I-1)=IDX(I)
          IDX(I)=INT
          I=I-1
          IF(I.GT.IS)GOTO  40
   60   CONTINUE
        LA=LA-2
        GOTO 170
!     ==--------------------------------------------------------------==
!     ==                *******  QUICKSORT  ********                  ==
!     == SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS  ==
!     == THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENT'S==
!     == HIGHEST ADDRESS.                                             ==
!     ==--------------------------------------------------------------==
   70   IY=(IS+IF)/2
        X=COUNT(IY)
        INTEST=IDX(IY)
        COUNT(IY)=COUNT(IF)
        IDX(IY)=IDX(IF)
!     ==--------------------------------------------------------------==
!     == THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END ==
!     == OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE   ==
!     == OF X .                                                       ==
!     ==--------------------------------------------------------------==
        K=1
        IFK=IF
!     ==--------------------------------------------------------------==
!     == WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE ==
!     == INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AND INDICES AS   ==
!     == NECESSARY, UNTIL THEY MEET .                                 ==
!     ==--------------------------------------------------------------==
        DO 110 I=IS,IF
          IF((X .GT. COUNT(I)))GOTO 110
          IF((X .LT. COUNT(I)))GOTO 80
          IF(INTEST.GT.IDX(I))GOTO 110
   80     IF(I.GE.IFK)GOTO 120
          COUNT(IFK)=COUNT(I)
          IDX(IFK)=IDX(I)
          K1=K
          DO 100 K=K1,IFKA
            IFK=IF-K
            IF((COUNT(IFK) .GT. X))GOTO 100
            IF((COUNT(IFK) .LT. X))GOTO 90
            IF(INTEST.LE.IDX(IFK))GOTO 100
   90       IF(I.GE.IFK)GOTO 130
            COUNT(I)=COUNT(IFK)
            IDX(I)=IDX(IFK)
            GO TO 110
  100     CONTINUE
          GOTO 120
  110   CONTINUE
!     ==--------------------------------------------------------------==
!     == RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER  ==
!     == WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO ==
!     == 2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN OR EQUAL ==
!     == TO ANY ELEMENT IN THE SECOND PART, AND THEY MAY NOW BE SORTED==
!     == INDEPENDENTLY .                                              ==
!     ==--------------------------------------------------------------==
  120   COUNT(IFK)=X
        IDX(IFK)=INTEST
        IP=IFK
        GOTO 140
  130   COUNT(I)=X
        IDX(I)=INTEST
        IP=I
!     ==--------------------------------------------------------------==
!     ==  STORE THE LONGER SUBDIVISION IN WORKSPACE.                  ==
!     ==--------------------------------------------------------------==
  140   IF((IP-IS).GT.(IF-IP))GOTO 150
        MARK(LA)=IF
        MARK(LA-1)=IP+1
        IF=IP-1
        GOTO 160
  150   MARK(LA)=IP-1
        MARK(LA-1)=IS
        IS=IP+1
!     ==--------------------------------------------------------------==
!     == FIND THE LENGTH OF THE SHORTER SUBDIVISION.                  ==
!     ==--------------------------------------------------------------==
  160   LNGTH=IF-IS
        IF(LNGTH.LE.0)GOTO 180
!     ==--------------------------------------------------------------==
!     == IF IT CONTAINS MORE THAN ONE ELEMENT SUPPLY IT WITH WORKSPACE==
!     ==--------------------------------------------------------------==
        LA=LA+2
        GOTO 190
  170   IF(LA.LE.0)GOTO 200
!     ==--------------------------------------------------------------==
!     == OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT==
!     ==--------------------------------------------------------------==
  180   IF=MARK(LA)
        IS=MARK(LA-1)
  190 CONTINUE
!     ==--------------------------------------------------------------==
  200 RETURN
      END SUBROUTINE rqsort
!     ==================================================================

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
      integer :: k1, ifk, lngth, ip, k, int, ifka, intest, iy
      integer :: i, m, la, is, if, mloop, ifca, is1, j, mark(50)
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
      if=n
      do 190 mloop=1,n
!  if segment is short enough sort with final sorting routine .
      ifka=if-is
      if((ifka+1).gt.m)goto 70
!********* final sorting ***
!  ( a simple bubble sort )
      is1=is+1
      do 60 j=is1,if
      i=j
   40 if(count(i-1).lt.count(i))goto 60
      if(count(i-1).gt.count(i))goto 50
      if(idx(i-1).lt.idx(i))goto 60
   50 av=count(i-1)
      count(i-1)=count(i)
      count(i)=av
      int=idx(i-1)
      idx(i-1)=idx(i)
      idx(i)=int
      i=i-1
      if(i.gt.is)goto  40
   60 continue
      la=la-2
      goto 170
!             *******  quicksort  ********
!  select the number in the central position in the segment as
!  the test number.replace it with the number from the segment's
!  highest address.
   70 iy=(is+if)/2
      x=count(iy)
      intest=idx(iy)
      count(iy)=count(if)
      idx(iy)=idx(if)
!  the markers 'i' and 'ifk' are used for the beginning and end
!  of the section not so far tested against the present value
!  of x .
      k=1
      ifk=if
!  we alternate between the outer loop that increases i and the
!  inner loop that reduces ifk, moving numbers and indices as
!  necessary, until they meet .
      do 110 i=is,if
      if(x.gt.count(i))goto 110
      if(x.lt.count(i))goto 80
      if(intest.gt.idx(i))goto 110
   80 if(i.ge.ifk)goto 120
      count(ifk)=count(i)
      idx(ifk)=idx(i)
      k1=k
      do 100 k=k1,ifka
      ifk=if-k
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
  140 if((ip-is).gt.(if-ip))goto 150
      mark(la)=if
      mark(la-1)=ip+1
      if=ip-1
      goto 160
  150 mark(la)=ip-1
      mark(la-1)=is
      is=ip+1
!  find the length of the shorter subdivision.
  160 lngth=if-is
      if(lngth.le.0)goto 180
!  if it contains more than one element supply it with workspace .
      la=la+2
      goto 190
  170 if(la.le.0)goto 10
!  obtain the address of the shortest segment awaiting quicksort
  180 if=mark(la)
      is=mark(la-1)
  190 continue
   10 return
      end subroutine kb07ad_cp90

