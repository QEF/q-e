!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE cgsolve (npw,evc,npwx,nbnd,overlap,      &
     &              nbndx,orthonormal,precondition,diagonal, &
     &              startwith0,e,b,u,h,Ah,pu,niter,eps,iter,x)
  !-----------------------------------------------------------------------
  !
  !  conjugate-gradient solution of a system of constrained linear equations
  !  "operator" is the linear operator - diagonal preconditioning allowed
  !  x = solution, u = gradient, h = conjugate gradient, Ah = operator*h
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE becmod,    ONLY : calbec
  USE io_global,   ONLY : ionode,stdout
  USE gvect,                 ONLY : gstart
  USE mp, ONLY : mp_sum, mp_barrier
  USE mp_world,             ONLY : world_comm

  IMPLICIT NONE
  INTEGER npw, npwx, nbnd, nbndx, niter, iter
  real(DP) :: diagonal(npw), e(nbnd), overlap(nbndx,nbnd)
  COMPLEX(DP) :: x(npwx,nbnd), b(npwx,nbnd), u(npwx,nbnd),          &
       h(npwx,nbnd),Ah(npwx,nbnd),evc(npwx,nbnd), pu(npwx,nbnd)
  LOGICAL :: orthonormal, precondition,startwith0
  !
  INTEGER :: ibnd, jbnd, i, info
  real(DP) :: lagrange(nbnd,nbnd)
  real(DP) :: lambda, u_u, uu0, u_A_h, alfa, eps, uu(nbnd), ddot
  logical ::  debug
  real(kind=DP), allocatable :: omat(:,:)
  integer j
  EXTERNAL ddot, H_h
  !
  CALL start_clock('cgsolve')
  debug=.false.
  !
  ! starting gradient |u> = (A|x>-|b>)-lambda|psi> (lambda=<Ax-b|psi_i>)
  !
  IF (.not.startwith0) THEN
     CALL H_h(e,x,u)
  ELSE
     u (:,:) = (0.d0, 0.d0)
     ! note that we assume x=0 on input
  ENDIF
  !

  allocate(omat(nbnd,nbnd))
  omat(1:nbnd,1:nbnd)=0.d0
  if(debug) then
     call dgemm('T','N',nbnd,nbnd,2*npw,2.d0,b,2*npwx,b,2*npwx,0.d0,omat,nbnd)
     if(gstart==2) then
        do i=1,nbnd
           do j=1,nbnd
              omat(i,j)=omat(i,j)-dble(b(1,i)*conjg(b(1,j)))
           enddo
        enddo
     endif
   
     call mp_sum(omat,world_comm)  
     do i=1,nbnd
        do j=1,nbnd
           if(ionode) write(*,*) 'nbnd1,nbnd2,omat4(1,2)', i,j,omat(i,j)
        enddo
     enddo
  endif


 
  CALL daxpy(2*npwx*nbnd,-1.d0,b,1,u,1)

  IF (precondition) THEN
     CALL zvscal(npw,npwx,nbnd,diagonal,u,pu)
     CALL calbec ( npw, evc, pu, lagrange )
  ELSE
     CALL calbec ( npw, evc,  u, lagrange )
   ENDIF
  IF (.not. orthonormal) &
       CALL DPOTRS('U',nbnd,nbnd,overlap,nbndx,lagrange,nbnd,info)
  IF (info/=0) CALL errore('cgsolve','error in potrs',info)
  !
  CALL dgemm ('N', 'N', 2*npw, nbnd, nbnd, -1.d0, evc, &
       2*npwx, lagrange, nbnd, 1.d0, u, 2*npwx)
  !
  ! starting conjugate gradient |h> = |u>
  IF (precondition) THEN
     CALL zvscal(npw,npwx,nbnd,diagonal,u,h)
  ELSE
     CALL zcopy(npwx,nbnd,u,1,h,1)
  ENDIF
  ! uu = <u|h>
  CALL pw_dot('Y',npw,nbnd,u,npwx,h,npwx,uu)
  u_u = 0.0d0
  DO ibnd=1,nbnd
     u_u = u_u + uu(ibnd)
  ENDDO
  !
  !      print '("  iter # ",i3,"  u_u = ",e10.4)', 0, u_u
  !
  !   main iteration loop
  !
  DO iter = 1, niter
  
     !
     ! calculate A|h>
     !
     CALL H_h(e,h(1,1),Ah(1,1))
     !
     ! u_A_h = <u|A|h> (NB: must be equal to <h|A|h>)
     IF (precondition) THEN
        CALL zvscal(npw,npwx,nbnd,diagonal,u,pu)
        ! uu = <u|PA|h>
        CALL pw_dot('Y',npw,nbnd,pu,npwx,Ah,npwx,uu)
     ELSE
        ! uu = <u|A|h>
        CALL pw_dot('Y',npw,nbnd, u,npwx,Ah,npwx,uu)
     ENDIF
     u_A_h = 0.0d0
     DO ibnd=1,nbnd
        u_A_h = u_A_h + uu(ibnd)
     ENDDO
     !
     lambda = - u_u / u_A_h
     ! update the gradient and the trial solution
     uu0 = u_u
     u_u = 0.0d0
     CALL daxpy(2*npwx*nbnd,lambda, h,1,x,1)
     CALL daxpy(2*npwx*nbnd,lambda,Ah,1,u,1)
     ! lagrange multipliers ensure orthogonality of the solution
     IF (precondition) THEN
        CALL zvscal(npw,npwx,nbnd,diagonal,u,pu)
        CALL calbec ( npw, evc, pu, lagrange )
     ELSE
        CALL calbec ( npw, evc,  u, lagrange )
     ENDIF
     IF (.not. orthonormal) &
          CALL DPOTRS('U',nbnd,nbnd,overlap,nbndx,lagrange,nbnd,info)
     IF (info/=0) CALL errore('cgsolve','error in potrs',info)
     CALL dgemm ('N', 'N', 2*npw, nbnd, nbnd,-1.d0, evc, &
          2*npwx, lagrange, nbndx, 1.d0, u, 2*npwx)
     IF (precondition) THEN
        CALL zvscal(npw,npwx,nbnd,diagonal,u,pu)
        ! uu = <u|A|u>
        CALL pw_dot('Y',npw,nbnd, u,npwx,pu,npwx,uu)
     ELSE
        ! uu = <u|u>
        CALL pw_dot('Y',npw,nbnd, u,npwx, u,npwx,uu)
     ENDIF
     u_u = 0.0d0
     DO ibnd=1,nbnd
        u_u = u_u + uu(ibnd)
     ENDDO
     !         print '("  iter # ",i3,"  u_u = ",e10.4)', iter, u_u
     !
     IF( u_u <= eps) GOTO 10
     IF (iter==niter) THEN
        WRITE( stdout,'("   *** Conjugate Gradient minimization",   &
             &    " not converged after ",i3," iterations"/ &
             &    " residual norm |Ax-b|^2 : ",e10.4)') iter,u_u
        GOTO 10
     ENDIF
     !   update the conjugate gradient
     alfa =  u_u / uu0
     DO ibnd = 1,nbnd
        IF (precondition) THEN
           DO i=1,npw
              h(i,ibnd) = alfa*h(i,ibnd) + u(i,ibnd)*diagonal(i)
           ENDDO
        ELSE
           DO i=1,npw
              h(i,ibnd) = alfa*h(i,ibnd) + u(i,ibnd)
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
10 CONTINUE
  CALL stop_clock('cgsolve')
  !
  deallocate(omat)
  RETURN
END SUBROUTINE cgsolve
