!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
      MODULE interfaces

        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE

        INTERFACE
          SUBROUTINE simpsn(inte,n,cl,s)
            USE kinds
            INTEGER, INTENT(IN) ::  N
            REAL(DP), INTENT(IN)  ::  INTE(N), cl
            REAL(DP), INTENT(OUT) ::  S
          END SUBROUTINE simpsn
        END INTERFACE 

!        INTERFACE
!          SUBROUTINE splinedx(xmin,xmax,y,n,yp1,ypn,y2)
!            USE kinds
!            INTEGER, INTENT(IN)  :: n
!            REAL(DP),  INTENT(IN)  :: yp1,ypn,xmin,xmax,y(n)
!            REAL(DP),  INTENT(OUT) :: y2(n) 
!          END SUBROUTINE splinedx
!        END INTERFACE 
!
!        INTERFACE
!          SUBROUTINE splintdx(xmin,xmax,ya,y2a,n,x,y)
!            USE kinds
!            INTEGER, INTENT(IN)  :: n
!            REAL(DP),  INTENT(IN)  :: x,xmin,xmax,ya(n),y2a(n)
!            REAL(DP),  INTENT(OUT) :: y 
!          END SUBROUTINE splintdx
!        END INTERFACE 

        INTERFACE
          FUNCTION miller2gsqr(i,j,k,b1,b2,b3)
            USE kinds
            REAL(DP) :: miller2gsqr
            INTEGER, INTENT(IN) :: i,j,k
            REAL(DP),  INTENT(IN) :: b1(3), b2(3), b3(3) 
          END FUNCTION miller2gsqr
        END INTERFACE 

        INTERFACE
          SUBROUTINE miller2nxh(i,j,k,n1h,n2h,n3h,nr1,nr2,nr3)  
            INTEGER, INTENT(IN)  :: i,j,k
            INTEGER, INTENT(OUT) :: n1h,n2h,n3h
            INTEGER, INTENT(IN)  :: nr1,nr2,nr3
          END SUBROUTINE miller2nxh
        END INTERFACE 

        INTERFACE
          SUBROUTINE miller2inxh(i,j,k,in1h,in2h,in3h,nr1,nr2,nr3)
            INTEGER, INTENT(IN)  :: i,j,k
            INTEGER, INTENT(OUT) :: in1h,in2h,in3h
            INTEGER, INTENT(IN)  :: nr1,nr2,nr3
          END SUBROUTINE miller2inxh
        END INTERFACE 

        INTERFACE
          SUBROUTINE miller2indx(i,j,k,ind1,ind2,ind3,nr1,nr2,nr3)
            INTEGER, INTENT(IN)  :: i,j,k
            INTEGER, INTENT(OUT) :: ind1,ind2,ind3
            INTEGER, INTENT(IN)  :: nr1,nr2,nr3
          END SUBROUTINE miller2indx
        END INTERFACE 

        INTERFACE
          SUBROUTINE inxh2miller(in1h,in2h,in3h,i,j,k,nr1,nr2,nr3)
            INTEGER, INTENT(OUT) :: i,j,k
            INTEGER, INTENT(IN)  :: in1h,in2h,in3h
            INTEGER, INTENT(IN)  :: nr1,nr2,nr3
          END SUBROUTINE inxh2miller
        END INTERFACE 

        INTERFACE
          SUBROUTINE inxh2nxh(in1h,in2h,in3h,n1h,n2h,n3h,nr1,nr2,nr3)
            INTEGER, INTENT(OUT) :: n1h,n2h,n3h
            INTEGER, INTENT(IN)  :: in1h,in2h,in3h
            INTEGER, INTENT(IN)  :: nr1,nr2,nr3  
          END SUBROUTINE inxh2nxh
        END INTERFACE 

        INTERFACE
          SUBROUTINE inxh2indx(in1h,in2h,in3h,ind1,ind2,ind3,nr1,nr2,nr3)
            INTEGER, INTENT(OUT) :: ind1,ind2,ind3
            INTEGER, INTENT(IN)  :: in1h,in2h,in3h
            INTEGER, INTENT(IN)  :: nr1,nr2,nr3 
          END SUBROUTINE inxh2indx
        END INTERFACE 

!        INTERFACE
!        END INTERFACE 

!        INTERFACE
!        END INTERFACE 


      END MODULE interfaces
