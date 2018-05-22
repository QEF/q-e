      subroutine lmdif1(fcn,m,n,n_max_iter,x,fvec,tol,info,iwa,wa,lwa)
      integer m,n,info,lwa
      integer iwa(n)
      double precision tol
      double precision x(n),fvec(m),wa(lwa)
      external fcn
!     **********
!
!     subroutine lmdif1
!
!     the purpose of lmdif1 is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of the
!     levenberg-marquardt algorithm. this is done by using the more
!     general least-squares solver lmdif. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!     the subroutine statement is
!
!       subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
!
!     where
!
!       fcn is the name of the user-supplied subroutine which
!         calculates the functions. fcn must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         subroutine fcn(m,n,x,fvec,iflag)
!         integer m,n,iflag
!         double precision x(n),fvec(m)
!         ----------
!         calculate the functions at x and
!         return this vector in fvec.
!         ----------
!         return
!         end
!
!         the value of iflag should not be changed by fcn unless
!         the user wants to terminate execution of lmdif1.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an array of length n. on input x must contain
!         an initial estimate of the solution vector. on output x
!         contains the final estimate of the solution vector.
!
!       fvec is an output array of length m which contains
!         the functions evaluated at the output x.
!
!       tol is a nonnegative input variable. termination occurs
!         when the algorithm estimates either that the relative
!         error in the sum of squares is at most tol or that
!         the relative error between x and the solution is at
!         most tol.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  algorithm estimates that the relative error
!                   in the sum of squares is at most tol.
!
!         info = 2  algorithm estimates that the relative error
!                   between x and the solution is at most tol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  fvec is orthogonal to the columns of the
!                   jacobian to machine precision.
!
!         info = 5  number of calls to fcn has reached or
!                   exceeded 200*(n+1).
!
!         info = 6  tol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  tol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!       iwa is an integer work array of length n.
!
!       wa is a work array of length lwa.
!
!       lwa is a positive integer input variable not less than
!         m*n+5*n+m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... lmdif
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer maxfev,mode,mp5n,nfev,nprint,n_max_iter,iflga
      double precision epsfcn,factor,ftol,gtol,xtol,zero
      info = 0
!
!     check the input parameters for errors.
!
      if (n .le. 0 .or. m .lt. n .or. tol .lt. zero                     &
     &    .or. lwa .lt. m*n + 5*n + m) go to 10
!
!     call lmdif.
!
      factor = 1.0d3
      zero = 0.0d0
      maxfev = n_max_iter*(n + 1)
      ftol = tol
      xtol = tol
      gtol = zero
      epsfcn = zero
      epsfcn = 1.d-9
      mode = 1
      nprint = 0
      mp5n = m + 5*n
! ATTENZIONE
      call fcn(m,n,x,fvec,iflga)
!      write(*,*) 'fvec',fvec(1:10)
!
!
      call lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,wa(1),     &
     &           mode,factor,nprint,info,nfev,wa(mp5n+1),m,iwa,         &
     &           wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
      if (info .eq. 8) info = 4
   10 continue
      return
!
!     last card of subroutine lmdif1.
!
      end
