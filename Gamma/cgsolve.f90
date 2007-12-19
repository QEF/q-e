!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine cgsolve (operator,npw,evc,npwx,nbnd,overlap,      &
     &              nbndx,orthonormal,precondition,diagonal, &
     &              startwith0,e,b,u,h,Ah,pu,niter,eps,iter,x)
  !-----------------------------------------------------------------------
  !
  !  conjugate-gradient solution of a system of constrained linear equations
  !  "operator" is the linear operator - diagonal preconditioning allowed
  !  x = solution, u = gradient, h = conjugate gradient, Ah = operator*h
  !
#include "f_defs.h"
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE becmod,    ONLY : calbec
  implicit none
  integer npw, npwx, nbnd, nbndx, niter, iter
  real(DP) :: diagonal(npw), e(nbnd), overlap(nbndx,nbnd)
  complex(DP) :: x(npwx,nbnd), b(npwx,nbnd), u(npwx,nbnd),          &
       h(npwx,nbnd),Ah(npwx,nbnd),evc(npwx,nbnd), pu(npwx,nbnd)
  logical :: orthonormal, precondition,startwith0
  !
  integer :: ibnd, jbnd, i, info
  real(DP) :: lagrange(nbnd,nbnd)
  real(DP) :: lambda, u_u, uu0, u_A_h, alfa, eps, uu(nbnd), DDOT
  external DDOT, operator
  !
  call start_clock('cgsolve')
  !
  ! starting gradient |u> = (A|x>-|b>)-lambda|psi> (lambda=<Ax-b|psi_i>)
  !
  if (.not.startwith0) then
     call operator(e,x,u)
  else
     u (:,:) = (0.d0, 0.d0)
     ! note that we assume x=0 on input
  end if
  !
  call DAXPY(2*npwx*nbnd,-1.d0,b,1,u,1)
  if (precondition) then
     call zvscal(npw,npwx,nbnd,diagonal,u,pu)
     call calbec ( npw, evc, pu, lagrange )
  else
     call calbec ( npw, evc,  u, lagrange )
   end if
  if (.not. orthonormal) &
       call DPOTRS('U',nbnd,nbnd,overlap,nbndx,lagrange,nbnd,info)
  if (info.ne.0) call errore('cgsolve','error in potrs',info)
  !
  call DGEMM ('N', 'N', 2*npw, nbnd, nbnd, -1.d0, evc, &
       2*npwx, lagrange, nbndx, 1.d0, u, 2*npwx)
  !
  ! starting conjugate gradient |h> = |u>
  if (precondition) then
     call zvscal(npw,npwx,nbnd,diagonal,u,h)
  else
     call ZCOPY(npwx,nbnd,u,1,h,1)
  end if
  ! uu = <u|h>
  call pw_dot('Y',npw,nbnd,u,npwx,h,npwx,uu)
  u_u = 0.0d0
  do ibnd=1,nbnd
     u_u = u_u + uu(ibnd)
  end do
  !
  !      print '("  iter # ",i3,"  u_u = ",e10.4)', 0, u_u
  !
  !   main iteration loop
  !
  do iter = 1, niter
     !
     ! calculate A|h>
     !
     call operator(e,h,Ah)
     !
     ! u_A_h = <u|A|h> (NB: must be equal to <h|A|h>)
     if (precondition) then
        call zvscal(npw,npwx,nbnd,diagonal,u,pu)
        ! uu = <u|PA|h>
        call pw_dot('Y',npw,nbnd,pu,npwx,Ah,npwx,uu)
     else
        ! uu = <u|A|h>
        call pw_dot('Y',npw,nbnd, u,npwx,Ah,npwx,uu)
     end if
     u_A_h = 0.0d0
     do ibnd=1,nbnd
        u_A_h = u_A_h + uu(ibnd)
     end do
     !
     lambda = - u_u / u_A_h
     ! update the gradient and the trial solution
     uu0 = u_u
     u_u = 0.0d0
     call DAXPY(2*npwx*nbnd,lambda, h,1,x,1)
     call DAXPY(2*npwx*nbnd,lambda,Ah,1,u,1)
     ! lagrange multipliers ensure orthogonality of the solution
     if (precondition) then
        call zvscal(npw,npwx,nbnd,diagonal,u,pu)
        call calbec ( npw, evc, pu, lagrange )
     else
        call calbec ( npw, evc,  u, lagrange )
     end if
     if (.not. orthonormal) &
          call DPOTRS('U',nbnd,nbnd,overlap,nbndx,lagrange,nbnd,info)
     if (info.ne.0) call errore('cgsolve','error in potrs',info)
     call DGEMM ('N', 'N', 2*npw, nbnd, nbnd,-1.d0, evc, &
          2*npwx, lagrange, nbndx, 1.d0, u, 2*npwx)
     if (precondition) then
        call zvscal(npw,npwx,nbnd,diagonal,u,pu)
        ! uu = <u|A|u>
        call pw_dot('Y',npw,nbnd, u,npwx,pu,npwx,uu)
     else
        ! uu = <u|u>
        call pw_dot('Y',npw,nbnd, u,npwx, u,npwx,uu)
     end if
     u_u = 0.0d0
     do ibnd=1,nbnd
        u_u = u_u + uu(ibnd)
     end do
     !         print '("  iter # ",i3,"  u_u = ",e10.4)', iter, u_u
     !
     if( u_u .le. eps) go to 10
     if (iter.eq.niter) then
        WRITE( stdout,'("   *** Conjugate Gradient minimization",   &
             &    " not converged after ",i3," iterations"/ &
             &    " residual norm |Ax-b|^2 : ",e10.4)') iter,u_u
        go to 10
     end if
     !   update the conjugate gradient
     alfa =  u_u / uu0
     do ibnd = 1,nbnd
        if (precondition) then
           do i=1,npw
              h(i,ibnd) = alfa*h(i,ibnd) + u(i,ibnd)*diagonal(i)
           end do
        else
           do i=1,npw
              h(i,ibnd) = alfa*h(i,ibnd) + u(i,ibnd)
           end do
        end if
     end do
  end do
  !
10 continue
  call stop_clock('cgsolve')
  !
  return
end subroutine cgsolve
