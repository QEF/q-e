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

  MODULE global_minpack
!this module conatins global variables(sigh!) for using old FORTRAN77 
! minpack routine
  USE kinds, ONLY : DP

  IMPLICIT NONE

  SAVE
   
  INTEGER, PARAMETER :: maxm=400!max number of samples
  INTEGER, PARAMETER :: maxpole=30
  INTEGER :: n_poles
  COMPLEX(kind=DP) :: c_target(maxm)
  REAL(kind=DP) ::  freq(maxm)


  END MODULE




  SUBROUTINE  fit_multipole(n,m,z,s,a_0,a,b,delta,thres,maxiter)
!fits with the function f(z)=a_0+\sum_{i=1,m} a_i/(z-b_i)
!the values z_j,s_j

   USE kinds,         ONLY : DP
   USE io_global,     ONLY : stdout

   implicit none


   INTEGER, INTENT(in) :: n!numer of sampled values
   INTEGER, INTENT(in) :: m!number of parameters a
   COMPLEX(kind=DP), INTENT(in) :: z(n)!where 
   COMPLEX(kind=DP), INTENT(in) :: s(n)!values s(z_j) to be fitted
   COMPLEX(kind=DP), INTENT(inout) :: a_0
   COMPLEX(kind=DP), INTENT(inout) :: a(m)
   COMPLEX(kind=DP), INTENT(inout) :: b(m)
   REAL(kind=DP), INTENT(in) :: delta!parameter for steepest descend
   REAL(kind=DP), INTENT(in) :: thres!threshold for convergence
   INTEGER, INTENT(in) :: maxiter!maximum number of iterations


   REAL(kind=DP)  :: chi0, chi1
   INTEGER :: i,j,it
   COMPLEX(kind=DP) :: cc, grad
   COMPLEX(kind=DP) :: new_a_0, new_a(m), new_b(m), old_b(m), old_a(m), old_a_0
   REAL(kind=DP) :: rr,rc
   REAL(kind=DP) :: dd,ddb
   LOGICAL :: random

   INTEGER :: ip, im

   ddb= 0.01d0 
   dd = delta
   random=.true. 
   ip=1
   im=1
   
!calculates initial chi

   chi0=0.d0
   do i=1,n
     cc = func(z(i))-s(i)
     !cc = (func(z(i))-s(i))/s(i)
     chi0=chi0+cc*conjg(cc)
   enddo
  write(stdout,*) 'a_0', a_0!ATTENZIONE
   write(stdout,*) 'a', a
   write(stdout,*) 'b', b
   write(stdout,*) 'z,s' , z(1),s(1),func(z(1))
   write(stdout,*) 'z,s' , z(n),s(n),func(z(n))

   do it=1,maxiter

!updates a_0
     grad=(0.d0,0.d0)
     do i=1,n
       grad=grad+(func(z(i))-s(i))
     enddo
     new_a_0=a_0-dd*grad
     if(it==1) write(stdout,*) 'Grad a_0', grad!ATTENZIONE
!updates a(:)

    do j=1,m
      grad=(0.d0,0.d0)
      do i=1,n
        grad=grad+(func(z(i))-s(i))/(conjg(z(i))-conjg(b(j)))
      enddo
      new_a(j)=a(j)-grad*dd
      if(it==1) write(stdout,*) 'Grad a', grad!ATTENZIONE
    enddo

!updates b(:)

    if(.not.random) then
      do j=1,m
        grad=(0.d0,0.d0)
        do i=1,n
          grad=grad+(func(z(i))-s(i))*conjg(a(j))/((conjg(z(i))-conjg(b(j)))*(conjg(z(i))-conjg(b(j))))
        enddo
        new_b(j)=b(j)-grad*dd
        if(it==1) write(stdout,*) 'Grad b', grad!ATTENZIONE
      enddo
    endif

!calculates new chi

    old_a_0=a_0
    a_0=new_a_0
    old_a(:)=a(:)
    a(:)=new_a(:)
    if(.not. random) b(:)=new_b(:)
     
    chi1=0.d0
    do i=1,n
      cc = func(z(i))-s(i)
      !cc = (func(z(i))-s(i))/s(i)
      chi1=chi1+cc*conjg(cc)
    enddo
   
    if(chi1 > chi0) then
       a_0=old_a_0
       a(:)=old_a(:)
       write(stdout,*) 'Routine fit_multipole: chi1 > chi0 '
       !return
       dd=dd*0.1d0
    endif

!    if((chi0-chi1) <= thres) then
!       write(stdout,*)  'Reached threshold', chi0
!       return
!    endif

    chi0=chi1
    if(random) then!minimize b in a random fashion
!      if(mod(it,50000) == 1) ddb=ddb*0.1d0
       if(mod(ip,10)==0) then
          ddb=ddb*0.1d0
          ip=1
          write(stdout,*) 'Random plus'
       endif
       if(mod(im,100)==0) then
          ddb=ddb*10.d0
          im=1
          write(stdout,*) 'Random minus'
       endif

      do j=1,m
          old_b(:)=b(:)
          call random_number(rr)
          call random_number(rc)
          b(j)=b(j)+ddb*cmplx(rr,rc)
          if(aimag(b(j)) >= 0.d0) b(j)=cmplx(real(b(j)),-aimag(b(j)))
          chi1=0.d0
          do i=1,n
            cc = func(z(i))-s(i)
            !cc = (func(z(i))-s(i))/s(i)
           chi1=chi1+cc*conjg(cc)
         enddo
         if(chi1<chi0) then
      !    write(*,*) 'PASSED'
          chi0=chi1
          ip=ip+1
          im=1
        else
          b(:)=old_b(:)
          ip=1
          im=im+1
        endif
       enddo
    endif
    !write(stdout,*) 'chi0',chi0!ATTENZIONE
 enddo
  write(stdout,*) 'Routine fit_multipole: maxiter reached ', chi0  
  return

  CONTAINS

  FUNCTION func(zz)

  COMPLEX(kind=DP) :: func

  COMPLEX(kind=DP) :: zz
  INTEGER :: ii

  func=a_0
  do ii=1,m
    func=func+a(ii)/(zz-b(ii))
  enddo

  return
END FUNCTION func

END SUBROUTINE fit_multipole
  
    SUBROUTINE  fit_multipole_verlet(n,m,z,s,a_0,a,b,thres,maxiter,ma_0,ma,mb,dt,frice)
!fits with the function f(z)=a_0+\sum_{i=1,m} a_i/(z-b_i)
!the values z_j,s_j

   USE kinds,         ONLY : DP
   USE io_global,     ONLY : stdout

   implicit none


   INTEGER, INTENT(in) :: n!numer of sampled values
   INTEGER, INTENT(in) :: m!number of parameters a
   COMPLEX(kind=DP), INTENT(in) :: z(n)!where
   COMPLEX(kind=DP), INTENT(in) :: s(n)!values s(z_j) to be fitted
   COMPLEX(kind=DP), INTENT(inout) :: a_0
   COMPLEX(kind=DP), INTENT(inout) :: a(m)
   COMPLEX(kind=DP), INTENT(inout) :: b(m)
   REAL(kind=DP), INTENT(in) :: thres!threshold for convergence
   INTEGER, INTENT(in) :: maxiter!maximum number of iterations
   REAL(kind=DP) :: ma_0!fictitious mass relative to a_0
   REAL(kind=DP) :: ma!fictitious mass relative to a's
   REAL(kind=DP) :: mb!fictitious mass relative to b's
   REAL(kind=DP) :: dt!
   REAL(kind=DP)  :: frice

   REAL(kind=DP)  :: chi0, chi1
   INTEGER :: i,j,it
   COMPLEX(kind=DP) :: cc
   COMPLEX(kind=DP) :: a_0p,a_0m,ap(m),am(m),bp(m),bm(m)
   COMPLEX(kind=DP) :: va_0,va(m),vb(m),fa_0,fa(m),fb(m),va_0m,vam(m),vbm(m),va_0p,vap(m),vbp(m)

!calculates initial chi

   chi0=0.d0
   do i=1,n
     cc = func(z(i))-s(i)
     chi0=chi0+cc*conjg(cc)
   enddo
  write(*,*) 'a_0', a_0!ATTENZIONE
   write(*,*) 'a', a
   write(*,*) 'b', b
   write(*,*) 'z,s' , z(1),s(1),func(z(1))
   write(*,*) 'z,s' , z(n),s(n),func(z(n))
   write(*,*) 'M a_0:',ma_0
   write(*,*) 'M a:',ma
   write(*,*) 'M b:',mb
   write(*,*) 'Dt:', dt
   write(*,*) 'frice:', frice
!first SD step


!updates a_0
     fa_0=(0.d0,0.d0)
     do i=1,n
       fa_0=fa_0-(func(z(i))-s(i))
     enddo
!updates a(:)

    do j=1,m
      fa(j)=(0.d0,0.d0)
      do i=1,n
        fa(j)=fa(j)-(func(z(i))-s(i))/(conjg(z(i))-conjg(b(j)))
      enddo
    enddo

!updates b(:)

    do j=1,m
      fb(j)=(0.d0,0.d0)
      do i=1,n
        fb(j)=fb(j)-(func(z(i))-s(i))*conjg(a(j))/((conjg(z(i))-conjg(b(j)))*(conjg(z(i))-conjg(b(j))))
      enddo
    enddo
    
!step
   
   a_0p=a_0+0.5d0*fa_0*dt*dt/ma_0
   ap(:)=a(:)+0.5d0*fa(:)*(dt)**2.d0/ma
   bp(:)=b(:)+0.5d0*fb(:)*(dt)**2.d0/mb
    write(*,*) 'fa_0', fa_0
    write(*,*) 'fa', fa
    write(*,*) 'fb', fb
   a_0m=a_0
   a_0=a_0p
   am(:)=a(:)
   a=ap(:)
   bm(:)=b(:)
   b=bp(:)
   write(*,*) 'a_0', a_0,ma_0
    write(*,*) 'a', a
    write(*,*) 'b', b
  
!set intial velocities
   va_0m=(0.d0,0.d0)
   vam(:)=(0.d0,0.d0)
   vbm(:)=(0.d0,0.d0)
   va_0=fa_0*dt/ma_0
   va(:)=fa(:)*dt/ma
   vb(:)=fb(:)*dt/mb 


   do it=1,maxiter

!calculates forces

!updates a_0
     fa_0=(0.d0,0.d0)
     do i=1,n
       fa_0=fa_0-(func(z(i))-s(i))
     enddo
      fa_0=(0.d0,0.d0)!ATTENZIONE
!updates a(:)

    do j=1,m
      fa(j)=(0.d0,0.d0)
      do i=1,n
        fa(j)=fa(j)-(func(z(i))-s(i))/(conjg(z(i))-conjg(b(j)))
      enddo
    enddo

!updates b(:)

    do j=1,m
      fb(j)=(0.d0,0.d0)
      do i=1,n
        fb(j)=fb(j)-(func(z(i))-s(i))*conjg(a(j))/((conjg(z(i))-conjg(b(j)))*(conjg(z(i))-conjg(b(j))))
      enddo
    enddo

!add of frice
!write(*,*) 'fa_0', fa_0
!    write(*,*) 'fa', fa
!    write(*,*) 'fb', fb

    fa_0=fa_0-frice*va_0
    fa(:)=fa(:)-frice*va(:)
    fb(:)=fb(:)-frice*vb(:)
!    write(*,*) 'fa_0', fa_0
!    write(*,*) 'fa', fa
!    write(*,*) 'fb', fb
    
!update positions

    a_0p=2.d0*a_0-a_0m+(fa_0/ma_0)*(dt)**2.d0
    ap(:)=2.d0*a(:)-am(:)+(fa(:)/ma)*(dt)**2.d0
    bp(:)=2.d0*b(:)-bm(:)+(fb(:)/mb)*(dt)**2.d0

    a_0m=a_0
    am(:)=a(:)
    bm(:)=b(:)
    a_0=a_0p
    a(:)=ap(:)
    b(:)=bp(:)

!calculates new chi


    chi1=0.d0
    do i=1,n
      cc = func(z(i))-s(i)
      chi1=chi1+cc*conjg(cc)
    enddo

!calculates velocities but at t-dt/2 ...
    va_0=(a_0-a_0m)/dt
    va(:)=(a(:)-am(:))/dt
    vb(:)=(b(:)-bm(:))/dt



!    if((chi0-chi1) <= thres) return
    chi0=chi1
    if(mod(it,100)==1)  write(*,*) chi0!ATTENZIONE
  enddo
  write(stdout,*) 'Routine fit_multipole: maxiter reached '
  return

  CONTAINS

  FUNCTION func(zz)

  COMPLEX(kind=DP) :: func

  COMPLEX(kind=DP) :: zz
  INTEGER :: ii

  func=a_0
  do ii=1,m
    func=func+a(ii)/(zz-b(ii))
  enddo

  return
  END FUNCTION func

  END SUBROUTINE fit_multipole_verlet


!the following routines are for MINPACK
  SUBROUTINE fcn_set(m, npoles, omegas,ctarget)
!added internal parmeters to set up function
!the parameters are in the order re(a_0),im(a_0),re(a_1)...re(a_n),im(a_1)..im(a_n),re(b_1)..re(b_n), im(b_1)..im(b_n)
  use kinds, ONLY : DP
  use io_global, ONLY : stdout
  use global_minpack


  implicit none

  INTEGER  :: m !number of variables
  INTEGER :: npoles
  REAL(kind=DP) :: omegas(m)
  COMPLEX(kind=DP) :: ctarget(m)

  n_poles=npoles
  freq(1:m)=omegas(1:m)
  c_target(1:m)=ctarget(1:m)
  return

 END SUBROUTINE


  SUBROUTINE fcn(m,n,x,fvec,iflag)
!added internal parmeters to set up function
!the parameters are in the order re(a_0),im(a_0),re(a_1)...re(a_n),im(a_1)..im(a_n),re(b_1)..re(b_n), im(b_1)..im(b_n)
  use kinds, ONLY : DP
  use io_global, ONLY : stdout
  use global_minpack
  

  implicit none

  INTEGER  :: m !number of variables
  INTEGER  :: n!total number of parameters
  REAL(kind=DP) :: x(n)!parameters
  REAL(kind=DP) :: fvec(m)!evaluated error function
  INTEGER :: iflag!not used

  COMPLEX(kind=DP) a_0, a(maxpole), b(maxpole)
  INTEGER i,j
  COMPLEX(kind=DP) :: func,zz


  if(m>maxm) then
    write(stdout,*) 'FCN: MAXN TOO SMALL'
    stop
  endif



!check number of parameters

  if(n /= 4*n_poles+2) then
    write(stdout,*) 'FCN: WRONG NUMBER OF PARAMETERS',n,n_poles
    stop
  endif

  if(n_poles>maxpole) then
   write(stdout,*) 'FCN: MAXPOLE TOO SMALL'
    stop
  endif
!set up parameters
  a_0=cmplx(x(1),x(2))
  do i=1,n_poles
    a(i)=cmplx(x(i*2+1),x(i*2+2))
  enddo
  do i=1,n_poles
    !b(i)=cmplx(x((i+n_poles)*2+1),-(x((i+n_poles)*2+2))**2.d0)
     b(i)=cmplx(x((i+n_poles)*2+1),x((i+n_poles)*2+2))!ATTENZIONE
  enddo  

!perform calculaation

 do i=1,m
   fvec(i)=0.d0
   func=a_0
   zz=cmplx(0.d0,freq(i))

   do j=1,n_poles
     func=func+a(j)/(zz-b(j))
   enddo
   func=func-c_target(i)
   !fvec(i)=sqrt(real(func*conjg(func)))
   fvec(i)=dble(func*conjg(func))

!   do j=1,n_poles
!     func=func+a(j)/(zz-b(j))
!   enddo
!   func=(func-c_target(i))/c_target(i)
!   fvec(i)=real(func*conjg(func))


 enddo
!  write(*,*) 'fcn', fvec(1)
  return
  END SUBROUTINE fcn

  SUBROUTINE  fit_multipole_verlet2(n,m,z,s,a_0,a,b,thres,n_max_iterations, chi, dt, frice)
!fits with the function f(z)=a_0+\sum_{i=1,m} a_i/(z-b_i)
!the values z_j,s_j

   USE kinds,         ONLY : DP
   USE io_global,     ONLY : stdout

   implicit none


   INTEGER, INTENT(in) :: n!numer of sampled values
   INTEGER, INTENT(in) :: m!number of parameters a
   COMPLEX(kind=DP), INTENT(in) :: z(n)!where
   COMPLEX(kind=DP), INTENT(in) :: s(n)!values s(z_j) to be fitted
   COMPLEX(kind=DP), INTENT(inout) :: a_0
   COMPLEX(kind=DP), INTENT(inout) :: a(m)
   COMPLEX(kind=DP), INTENT(inout) :: b(m)
   REAL(kind=DP), INTENT(in) :: thres!threshold for convergence
   INTEGER, INTENT(in) :: n_max_iterations!maximum number of search iterations
   REAL(kind=DP), INTENT(out) :: chi!the final chi
   REAL(kind=DP), INTENT(in) :: dt!time step
   REAL(kind=DP), INTENT(in) :: frice!frice


   REAL(kind=DP)  :: chi0, chi1, dtt
   INTEGER :: i,j,it
   COMPLEX(kind=DP) :: cc
   REAL(kind=DP), ALLOCATABLE :: omegas(:)
   REAL(kind=DP), ALLOCATABLE :: x0(:),x(:),v(:),f(:),ma(:),x1(:)
   INTEGER :: iflag
   INTEGER :: np

   np=2+4*m!number of parameters 
   allocate(omegas(n))
   allocate(x0(np),x(np),v(np),f(np),ma(np),x1(np))
   omegas(1:n)=aimag(z(1:n))

!calculates initial chi

   chi0=0.d0
   do i=1,n
     cc = func(z(i))-s(i)
     chi0=chi0+cc*conjg(cc)
   enddo

   write(stdout,*) 'Chi0 initial:', chi0

   !set in variables
   x(1)=real(a_0)
   x(2)=aimag(a_0)
   do i=1,m
     x(i*2+1)=real(a(i))
     x(i*2+2)=aimag(a(i))
     x((i+m)*2+1)=real(b(i))
     x((i+m)*2+2)=aimag(b(i))
  enddo

!set up fcn function paramters
  call fcn_set(n, m, omegas,s)

!intial values

  ma(:)=1.d0
  x0(:)=x(:)
  call fcn_point(n,np,x,chi1,f)
  write(stdout,*) 'VERLET2', chi1
  x(:)=x0(:)+f(:)*dt/ma(:)
  v(:)=0.d0
  chi0=chi1
  dtt=dt
  do it=1,n_max_iterations
     call fcn_point(n,np,x,chi1,f)
     if(chi1 >= chi0) dtt=dtt/10.d0
     chi0=chi1
     if(mod(it,1000)==1) write(stdout,*) 'VERLET2', it, chi1

     !x1(:)=2.d0*x(:)-x0(:)+f(:)*(dt**2.d0)/ma(:)-frice*v(:)*(dt**2.d0)
     
     x1(:)=x(:)+f(:)*dtt/ma(:)

     v(:)=(x1(:)-x0(:))/(2.d0*dt)

     x0(:)=x(:)
     x(:)=x1(:)
  enddo

  !set back parameters
  a_0=cmplx(x(1),x(2))
  do i=1,m
    a(i)=cmplx(x(i*2+1),x(i*2+2))
  enddo
  do i=1,m
    !b(i)=cmplx(x((i+m)*2+1),(-1.d0)*((x((i+m)*2+2))**2.d0))
     b(i)=cmplx(x((i+m)*2+1),x((i+m)*2+2))             
 enddo
  chi=chi1

  write(stdout,*) 'FINAL CHI', chi!ATTENZIONE





  deallocate(x,v,f,x0,ma,x1)

  return


  CONTAINS

  FUNCTION func(zz)

  COMPLEX(kind=DP) :: func

  COMPLEX(kind=DP) :: zz
  INTEGER :: ii

  func=a_0
  do ii=1,m
    func=func+a(ii)/(zz-b(ii))
  enddo

  return
  
END FUNCTION func
  

END SUBROUTINE fit_multipole_verlet2

  SUBROUTINE  fit_multipole_minpack(n,m,z,s,a_0,a,b,thres,n_max_iterations, chi)
!fits with the function f(z)=a_0+\sum_{i=1,m} a_i/(z-b_i)
!the values z_j,s_j

   USE kinds,         ONLY : DP
   USE io_global,     ONLY : stdout

   implicit none


   INTEGER, INTENT(in) :: n!numer of sampled values
   INTEGER, INTENT(in) :: m!number of parameters a
   COMPLEX(kind=DP), INTENT(in) :: z(n)!where
   COMPLEX(kind=DP), INTENT(in) :: s(n)!values s(z_j) to be fitted
   COMPLEX(kind=DP), INTENT(inout) :: a_0
   COMPLEX(kind=DP), INTENT(inout) :: a(m)
   COMPLEX(kind=DP), INTENT(inout) :: b(m)
   REAL(kind=DP), INTENT(in) :: thres!threshold for convergence
   INTEGER, INTENT(in) :: n_max_iterations!maximum number of search iterations
   REAL(kind=DP), INTENT(out) :: chi!the final chi

   REAL(kind=DP)  :: chi0, chi1
   INTEGER :: i,j,it
   COMPLEX(kind=DP) :: cc
   REAL(kind=DP), ALLOCATABLE :: omegas(:)
   REAL(kind=DP), ALLOCATABLE :: variables(:)
   INTEGER :: iflag
   INTEGER :: info
   INTEGER :: lwa
   INTEGER :: np
   INTEGER, ALLOCATABLE :: iwa(:)
   INTEGER, ALLOCATABLE :: ipvt(:)
   REAL(kind=DP), ALLOCATABLE :: wa(:)
   REAL(kind=DP), ALLOCATABLE :: fvec(:),fjac(:,:)
   EXTERNAL fcn,fcn2,fcnj
   INTEGER :: ldfjac

   ldfjac=n

   allocate(omegas(n))
   allocate(variables(2+4*m))
   allocate(iwa(2+4*m))
   lwa=n*(2+4*m)+5*(2+4*m)+n
   allocate(wa(lwa))
   allocate(fvec(n))

   np=2+4*m
   allocate(ipvt(np))
   allocate(fjac(n,np))
   write(stdout,*) 'Allocated'
   FLUSH(stdout)

   omegas(1:n)=aimag(z(1:n))

!calculates initial chi

   chi0=0.d0
   do i=1,n
     cc = func(z(i))-s(i)
     chi0=chi0+cc*conjg(cc)
   enddo
 ! write(*,*) 'a_0', a_0!ATTENZIONE
 !  write(*,*) 'a', a
 !  write(*,*) 'b', b
 !  write(*,*) 'z,s' , z(1),s(1),func(z(1))
 !  write(*,*) 'z,s' , z(n),s(n),func(z(n))
   write(stdout,*) 'Chi0 initial:', chi0
   FLUSH(stdout)
!set in variables
   variables(1)=real(a_0)
   variables(2)=aimag(a_0)
   do i=1,m
     variables(i*2+1)=real(a(i))
     variables(i*2+2)=aimag(a(i))
     variables((i+m)*2+1)=real(b(i))
     !variables((i+m)*2+2)=sqrt(abs(aimag(b(i))))
     variables((i+m)*2+2)=aimag(b(i))!ATTENZIONE 
   enddo   

!set up fcn function paramters
   call fcn_set(n, m, omegas,s)

!call minpack driver routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   info=1
   call fcn(n,np,variables,fvec,info)
   !write(*,*) 'FVEC', fvec(1:10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !call  lmdif1(fcn,n,np,n_max_iterations,variables,fvec,thres,info,iwa,wa,lwa)

   call lmder1(fcnj,n,np,variables,fvec,fjac,ldfjac,thres,info,ipvt,wa,lwa,n_max_iterations)



   write(stdout,*) ' INFO :', info,thres!ATTENZIONE
!set back parameters
  a_0=dcmplx(variables(1),variables(2))
  do i=1,m
    a(i)=dcmplx(variables(i*2+1),variables(i*2+2))
  enddo
  do i=1,m
    !b(i)=dcmplx(variables((i+m)*2+1),(-1.d0)*((variables((i+m)*2+2))**2.d0))
     b(i)=dcmplx(variables((i+m)*2+1),variables((i+m)*2+2))
  enddo

!recalculate chi and write result on stdout
!calculates initial chi

   chi0=0.d0
   do i=1,n
     cc = (func(z(i))-s(i))
     chi0=chi0+cc*conjg(cc)
   enddo
   !write(stdout,*) 'a_0', a_0!ATTENZIONE
   !write(stdout,*) 'a', a
   !write(stdout,*) 'b', b
   !write(stdout,*) 'z,s' , z(1),s(1),func(z(1))
   !write(stdout,*) 'z,s' , z(n),s(n),func(z(n))
   write(stdout,*) 'Minpack fit chi0 :', chi0
   chi=chi0


  deallocate(omegas)
  deallocate(variables)
  deallocate(iwa)
  deallocate(wa)
  deallocate(fvec)
  deallocate(ipvt)
  deallocate(fjac)

  return

  CONTAINS

  FUNCTION func(zz)

  COMPLEX(kind=DP) :: func

  COMPLEX(kind=DP) :: zz
  INTEGER :: ii

  func=a_0
  do ii=1,m
    func=func+a(ii)/(zz-b(ii))
  enddo

  return
  END FUNCTION func

  END SUBROUTINE fit_multipole_minpack
  SUBROUTINE fcn2(m,n,x,fvec,iflag)
!added internal parmeters to set up function
!the parameters are in the order re(a_0),im(a_0),re(a_1)...re(a_n),im(a_1)..im(a_n),re(b_1)..re(b_n), im(b_1)..im(b_n)
  use kinds, ONLY : DP
  use io_global, ONLY : stdout

  implicit none

  INTEGER  :: m !number of variables
  INTEGER  :: n!total number of parameters
  REAL(kind=DP) :: x(n)!parameters
  REAL(kind=DP) :: fvec(m)!evaluated error function
  INTEGER :: iflag!not used
  fvec(1:m)=0.d0

 END SUBROUTINE


  SUBROUTINE fcnj(m,n,x,fvec,fjac,ldfjac,iflag)
!this version calculates also the jacobian
!the parameters are in the order re(a_0),im(a_0),re(a_1)...re(a_n),im(a_1)..im(a_n),re(b_1)..re(b_n), im(b_1)..im(b_n)
  use kinds, ONLY : DP
  use io_global, ONLY : stdout
  use global_minpack


  implicit none

  INTEGER  :: m !number of variables
  INTEGER  :: n!total number of parameters
  REAL(kind=DP) :: x(n)!parameters
  REAL(kind=DP) :: fvec(m)!evaluated error function
  REAL(kind=DP) :: fjac(ldfjac,n)
  INTEGER :: ldfjac!leading dimension of fjac
  INTEGER :: iflag! =1 calculate fvec, =2 calculate fjac

  COMPLEX(kind=DP) a_0, a(maxpole), b(maxpole),g,h
  INTEGER i,j
  COMPLEX(kind=DP) :: func,zz


  if(m>maxm) then
    write(stdout,*) 'FCN: MAXN TOO SMALL'
    stop
  endif
!set up parameters
  a_0=dcmplx(x(1),x(2))
  do i=1,n_poles
    a(i)=dcmplx(x(i*2+1),x(i*2+2))
  enddo
  do i=1,n_poles
    !b(i)=dcmplx(x((i+n_poles)*2+1),-(x((i+n_poles)*2+2))**2.d0)
     b(i)=dcmplx(x((i+n_poles)*2+1),x((i+n_poles)*2+2))
  enddo

  if(iflag==1) then
!perform calculaation

     do i=1,m
        fvec(i)=0.d0
        func=a_0
        zz=dcmplx(0.d0,freq(i))
        
        do j=1,n_poles
           func=func+a(j)/(zz-b(j))
        enddo
        func=func-c_target(i)
        fvec(i)=dble(func*conjg(func))

     enddo

     else if(iflag==2) then
         do j=1,m
            fjac(j,:)=0.d0
!calculate g_j
            g=a_0
            zz=cmplx(0.d0,freq(j))
            do i=1,n_poles
               g=g+a(i)/(zz-b(i))
            enddo
            g=g-c_target(j)

!now term a_0
            fjac(j,1)=2.d0*real(g)
            fjac(j,2)=2.d0*aimag(g)
!now terms a_i
            do i=1, n_poles
               h=(1.d0,0.d0)/(zz-b(i))
               fjac(j,i*2+1)=2.d0*real(h*conjg(g))
               fjac(j,i*2+2)=-2.d0*aimag(h*conjg(g))
            enddo
!now terms b_i
            do i=1, n_poles
               h=a(i)/((zz-b(i))**2.d0)
               fjac(j,(i+n_poles)*2+1)=2.d0*real(h*conjg(g))
              !fjac(j,(i+n_poles)*2+2)=-2.d0*aimag(h*conjg(g))*(-2.d0*x((i+n_poles)*2+2)) 
              fjac(j,(i+n_poles)*2+2)=-2.d0*aimag(h*conjg(g))
            enddo

         enddo
     endif
  return
END SUBROUTINE fcnj


  SUBROUTINE fcn_point(m,n,x,value,nabla)
!this version calculates also the jacobian
!the parameters are in the order re(a_0),im(a_0),re(a_1)...re(a_n),im(a_1)..im(a_n),re(b_1)..re(b_n), im(b_1)..im(b_n)
  use kinds, ONLY : DP
  use io_global, ONLY : stdout
  use global_minpack


  implicit none

  INTEGER  :: m !number of variables
  INTEGER  :: n!total number of parameters
  REAL(kind=DP) :: x(n)!parameters
  REAL(kind=DP) :: value!total error
  REAL(kind=DP) :: nabla(n)!derivatives of total error respect to parameters

  COMPLEX(kind=DP) a_0, a(maxpole), b(maxpole),g,h
  INTEGER i,j
  COMPLEX(kind=DP) :: func,zz


  if(m>maxm) then
    write(stdout,*) 'FCN: MAXN TOO SMALL'
    stop
  endif
!set up parameters
  a_0=cmplx(x(1),x(2))
  do i=1,n_poles
    a(i)=cmplx(x(i*2+1),x(i*2+2))
  enddo
  do i=1,n_poles
  !   b(i)=cmplx(x((i+n_poles)*2+1),-(x((i+n_poles)*2+2))**2.d0)
     b(i)=cmplx(x((i+n_poles)*2+1),x((i+n_poles)*2+2))  
  enddo


!perform calculation of value

  value=0.d0
  do i=1,m
     func=a_0
     zz=cmplx(0.d0,freq(i))

     do j=1,n_poles
        func=func+a(j)/(zz-b(j))
     enddo
     func=func-c_target(i)
     value=value + dble(func*conjg(func))

  enddo

  nabla(:)=0.d0
  do j=1,m
!calculate g_j
     g=a_0
     zz=cmplx(0.d0,freq(j))
     do i=1,n_poles
        g=g+a(i)/(zz-b(i))
     enddo
     g=g-c_target(j)

!now term a_0
     nabla(1)=nabla(1)+2.d0*real(g)
     nabla(2)=nabla(2)+2.d0*aimag(g)
!now terms a_i
     do i=1, n_poles
        h=(1.d0,0.d0)/(zz-b(i))
        nabla(i*2+1)=nabla(i*2+1)+2.d0*real(h*conjg(g))
        nabla(i*2+2)=nabla(i*2+1)-2.d0*aimag(h*conjg(g))
     enddo
!now terms b_i
     do i=1, n_poles
        h=a(i)/((zz-b(i))**2.d0)
        nabla((i+n_poles)*2+1)=nabla((i+n_poles)*2+1)+2.d0*real(h*conjg(g))
        nabla((i+n_poles)*2+2)=nabla((i+n_poles)*2+2)-2.d0*aimag(h*conjg(g))
     enddo

  enddo
  nabla(:)=-nabla(:)

  return
END SUBROUTINE fcn_point
