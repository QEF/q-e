!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this contains routines to fit a function on the positive part
!of the imaginary axes with a multipole expansion


  SUBROUTINE  fit_polynomial(n,m,z,s,a,b,delta,thres,maxiter,maxcycle)
!fits with the function f(z)=(Prod_i(z-a_i)/Prod_j(z-b_j)
!the values z_j,s_j

   USE kinds,         ONLY : DP
   USE io_global,     ONLY : stdout

   implicit none


   INTEGER, INTENT(in) :: n!numer of sampled values
   INTEGER, INTENT(in) :: m!number of parameters a
   COMPLEX(kind=DP), INTENT(in) :: z(n)!where 
   COMPLEX(kind=DP), INTENT(in) :: s(n)!values s(z_j) to be fitted
   COMPLEX(kind=DP), INTENT(inout) :: a(m)
   COMPLEX(kind=DP), INTENT(inout) :: b(m)
   REAL(kind=DP), INTENT(in) :: delta!parameter for steepest descend
   REAL(kind=DP), INTENT(in) :: thres!threshold for convergence
   INTEGER, INTENT(in) :: maxiter!maximum number of iterations
   INTEGER, INTENT(in) :: maxcycle!maximum number of cycles


   REAL(kind=DP)  :: chi0, chi1
   INTEGER :: i,j,k,it,ic
   COMPLEX(kind=DP) :: cc, grad
   COMPLEX(kind=DP) ::  new_a(m), new_b(m), old_a(m),old_b(m)
   COMPLEX(kind=DP) :: num,den
   REAL(kind=DP) :: dd

   dd = delta
  
!calculates initial chi

   write(*,*) 'a', a
   write(*,*) 'b', b
   write(*,*) 'z,s' , z(1),s(1),func(z(1))
   write(*,*) 'z,s' , z(n),s(n),func(z(n))

   do ic=1,maxcycle

!now fit a's
     chi0=0.d0
     do i=1,n
       cc = func(z(i))-s(i)
       chi0=chi0+cc*conjg(cc)
     enddo
     do it=1,maxiter

!updates a(:)

     do j=1,m
       grad=(0.d0,0.d0)
       do i=1,n
         den=(1.d0,0.d0)
         do k=1,m
           den=den*(conjg(z(i))-conjg(b(k)))
         enddo
         num=(1.d0,0.d0)
         do k=1,m
            if(k/=j) num=num*(conjg(z(i))-conjg(a(k)))
         enddo
         grad=grad+(func(z(i))-s(i))*(-1.d0,0.d0)*num/den
         if(it==1) write(*,*) 'Grad a', grad,num,den!ATTENZIONE
       enddo
       new_a(j)=a(j)-grad*dd
       if(it==1) write(*,*) 'Grad a', grad,num,den!ATTENZIONE
     enddo

!calculates new chi

     a(:)=new_a(:)
     
     chi1=0.d0
     do i=1,n
       cc = func(z(i))-s(i)
       chi1=chi1+cc*conjg(cc)
     enddo
   
     if(chi1 > chi0) then
        write(stdout,*) 'Routine fit_multipole: chi1 > chi0 '
        !return
        dd=dd*0.1d0
     endif

     if((chi0-chi1) <= thres) return
     chi0=chi1
     write(*,*) 'chi0',chi0!ATTENZIONE
    enddo
!now fit b's
     chi0=0.d0
     do i=1,n
       cc = (1.d0/func(z(i)))-(1.d0/s(i))
       chi0=chi0+cc*conjg(cc)
     enddo
     do it=1,maxiter

!updates b(:)

     do j=1,m
       grad=(0.d0,0.d0)
       do i=1,n
         den=(1.d0,0.d0)
         do k=1,m
           den=den*(conjg(z(i))-conjg(a(k)))
         enddo
         num=(1.d0,0.d0)
         do k=1,m
            if(k/=j) num=num*(conjg(z(i))-conjg(b(k)))
         enddo
         grad=grad+((1.d0/func(z(i)))-(1.d0/s(i)))*(-1.d0,0.d0)*num/den
       enddo
       new_b(j)=b(j)-grad*dd
       if(it==1) write(*,*) 'Grad b', grad!ATTENZIONE
     enddo

!calculates new chi
     old_b(:)=b(:)
     b(:)=new_b(:)

     chi1=0.d0
     do i=1,n
       cc = (1.d0/func(z(i)))-(1.d0/s(i))
       chi1=chi1+cc*conjg(cc)
     enddo

     if(chi1 > chi0) then
        write(stdout,*) 'Routine fit_multipole: chi1 > chi0 '
        !return
        b(:)=old_b(:)
        dd=dd*0.1d0
     else if((chi0-chi1) <= thres) then
       return
     else 
       chi0=chi1
     endif
     write(*,*) 'chi0',chi0!ATTENZIONE
   enddo
   enddo


   write(stdout,*) 'Routine fit_multipole: maxcycle reached '  
   return

   CONTAINS

  FUNCTION func(zz)

  COMPLEX(kind=DP) :: func
  COMPLEX(kind=DP) :: n,d
  COMPLEX(kind=DP) :: zz
  INTEGER :: ii,kk

  n=(1.d0,0.d0)
  do kk=1,m
    n=n*(zz-a(kk))
  enddo
  d=(1.d0,0.d0)
  do kk=1,m
    d=d*(zz-b(kk))
  enddo
  func=n/d

  return
  END FUNCTION func

  END SUBROUTINE fit_polynomial
  
   
