      subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,     &
     &                 diag,mode,factor,nprint,info,nfev,fjac,ldfjac,   &
     &                 ipvt,qtf,wa1,wa2,wa3,wa4)
      integer m,n,maxfev,mode,nprint,info,nfev,ldfjac
      integer ipvt(n)
      double precision ftol,xtol,gtol,epsfcn,factor
      double precision x(n),fvec(m),diag(n),fjac(ldfjac,n),qtf(n),      &
     &                 wa1(n),wa2(n),wa3(n),wa4(m)
      external fcn
!     **********
!
!     subroutine lmdif
!
!     the purpose of lmdif is to minimize the sum of the squares of
!     m nonlinear functions in n variables by a modification of
!     the levenberg-marquardt algorithm. the user must provide a
!     subroutine which calculates the functions. the jacobian is
!     then calculated by a forward-difference approximation.
!
!     the subroutine statement is
!
!       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
!                        diag,mode,factor,nprint,info,nfev,fjac,
!                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
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
!         the user wants to terminate execution of lmdif.
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
!       ftol is a nonnegative input variable. termination
!         occurs when both the actual and predicted relative
!         reductions in the sum of squares are at most ftol.
!         therefore, ftol measures the relative error desired
!         in the sum of squares.
!
!       xtol is a nonnegative input variable. termination
!         occurs when the relative error between two consecutive
!         iterates is at most xtol. therefore, xtol measures the
!         relative error desired in the approximate solution.
!
!       gtol is a nonnegative input variable. termination
!         occurs when the cosine of the angle between fvec and
!         any column of the jacobian is at most gtol in absolute
!         value. therefore, gtol measures the orthogonality
!         desired between the function vector and the columns
!         of the jacobian.
!
!       maxfev is a positive integer input variable. termination
!         occurs when the number of calls to fcn is at least
!         maxfev by the end of an iteration.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       diag is an array of length n. if mode = 1 (see
!         below), diag is internally set. if mode = 2, diag
!         must contain positive entries that serve as
!         multiplicative scale factors for the variables.
!
!       mode is an integer input variable. if mode = 1, the
!         variables will be scaled internally. if mode = 2,
!         the scaling is specified by the input diag. other
!         values of mode are equivalent to mode = 1.
!
!       factor is a positive input variable used in determining the
!         initial step bound. this bound is set to the product of
!         factor and the euclidean norm of diag*x if nonzero, or else
!         to factor itself. in most cases factor should lie in the
!         interval (.1,100.). 100. is a generally recommended value.
!
!       nprint is an integer input variable that enables controlled
!         printing of iterates if it is positive. in this case,
!         fcn is called with iflag = 0 at the beginning of the first
!         iteration and every nprint iterations thereafter and
!         immediately prior to return, with x and fvec available
!         for printing. if nprint is not positive, no special calls
!         of fcn with iflag = 0 are made.
!
!       info is an integer output variable. if the user has
!         terminated execution, info is set to the (negative)
!         value of iflag. see description of fcn. otherwise,
!         info is set as follows.
!
!         info = 0  improper input parameters.
!
!         info = 1  both actual and predicted relative reductions
!                   in the sum of squares are at most ftol.
!
!         info = 2  relative error between two consecutive iterates
!                   is at most xtol.
!
!         info = 3  conditions for info = 1 and info = 2 both hold.
!
!         info = 4  the cosine of the angle between fvec and any
!                   column of the jacobian is at most gtol in
!                   absolute value.
!
!         info = 5  number of calls to fcn has reached or
!                   exceeded maxfev.
!
!         info = 6  ftol is too small. no further reduction in
!                   the sum of squares is possible.
!
!         info = 7  xtol is too small. no further improvement in
!                   the approximate solution x is possible.
!
!         info = 8  gtol is too small. fvec is orthogonal to the
!                   columns of the jacobian to machine precision.
!
!       nfev is an integer output variable set to the number of
!         calls to fcn.
!
!       fjac is an output m by n array. the upper n by n submatrix
!         of fjac contains an upper triangular matrix r with
!         diagonal elements of nonincreasing magnitude such that
!
!                t     t           t
!               p *(jac *jac)*p = r *r,
!
!         where p is a permutation matrix and jac is the final
!         calculated jacobian. column j of p is column ipvt(j)
!         (see below) of the identity matrix. the lower trapezoidal
!         part of fjac contains information generated during
!         the computation of r.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       ipvt is an integer output array of length n. ipvt
!         defines a permutation matrix p such that jac*p = q*r,
!         where jac is the final calculated jacobian, q is
!         orthogonal (not stored), and r is upper triangular
!         with diagonal elements of nonincreasing magnitude.
!         column j of p is column ipvt(j) of the identity matrix.
!
!       qtf is an output array of length n which contains
!         the first n elements of the vector (q transpose)*fvec.
!
!       wa1, wa2, and wa3 are work arrays of length n.
!
!       wa4 is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      integer i,iflag,iter,j,l
      double precision actred,delta,dirder,epsmch,fnorm,fnorm1,gnorm,   &
     &                 one,par,pnorm,prered,p1,p5,p25,p75,p0001,ratio,  &
     &                 sum,temp,temp1,temp2,xnorm,zero
      double precision dpmpar,enorm
      data one,p1,p5,p25,p75,p0001,zero                                 &
     &     /1.0d0,1.0d-1,5.0d-1,2.5d-1,7.5d-1,1.0d-4,0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
      info = 0
      iflag = 0
      nfev = 0
!
!     check the input parameters for errors.
!
      if (n .le. 0 .or. m .lt. n .or. ldfjac .lt. m                     &
     &    .or. ftol .lt. zero .or. xtol .lt. zero .or. gtol .lt. zero   &
     &    .or. maxfev .le. 0 .or. factor .le. zero) go to 300
      if (mode .ne. 2) go to 20
      do 10 j = 1, n
         if (diag(j) .le. zero) go to 300
   10    continue
   20 continue
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      call fcn(m,n,x,fvec,iflag)
      nfev = 1
      if (iflag .lt. 0) go to 300
      fnorm = enorm(m,fvec)
!
!     initialize levenberg-marquardt parameter and iteration counter.
!
      par = zero
      iter = 1
!
!     beginning of the outer loop.
!
   30 continue
!
!        calculate the jacobian matrix.
!
         iflag = 2
         call fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa4)
         nfev = nfev + n
         if (iflag .lt. 0) go to 300
!
!        if requested, call fcn to enable printing of iterates.
!
         if (nprint .le. 0) go to 40
         iflag = 0
!         if (mod(iter-1,nprint) .eq. 0) c
         call fcn(m,n,x,fvec,iflag)
!         write(*,*) 'ATTENZIONE l', fvec(1:10)
         if (iflag .lt. 0) go to 300
   40    continue
         call fcn(m,n,x,fvec,iflag)
!         write(*,*) 'ATTENZIONE l', fvec(1:10)
!
!
!        compute the qr factorization of the jacobian.
!
         call qrfac(m,n,fjac,ldfjac,.true.,ipvt,n,wa1,wa2,wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
         if (iter .ne. 1) go to 80
         if (mode .eq. 2) go to 60
         do 50 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) .eq. zero) diag(j) = one
   50       continue
   60    continue
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
         do 70 j = 1, n
            wa3(j) = diag(j)*x(j)
   70       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta .eq. zero) delta = factor
   80    continue
!
!        form (q transpose)*fvec and store the first n components in
!        qtf.
!
         do 90 i = 1, m
            wa4(i) = fvec(i)
   90       continue
         do 130 j = 1, n
            if (fjac(j,j) .eq. zero) go to 120
            sum = zero
            do 100 i = j, m
               sum = sum + fjac(i,j)*wa4(i)
  100          continue
            temp = -sum/fjac(j,j)
            do 110 i = j, m
               wa4(i) = wa4(i) + fjac(i,j)*temp
  110          continue
  120       continue
            fjac(j,j) = wa1(j)
            qtf(j) = wa4(j)
  130       continue
!
!        compute the norm of the scaled gradient.
!
         gnorm = zero
         if (fnorm .eq. zero) go to 170
         do 160 j = 1, n
            l = ipvt(j)
            if (wa2(l) .eq. zero) go to 150
            sum = zero
            do 140 i = 1, j
               sum = sum + fjac(i,j)*(qtf(i)/fnorm)
  140          continue
            gnorm = dmax1(gnorm,dabs(sum/wa2(l)))
  150       continue
  160       continue
  170    continue
!
!        test for convergence of the gradient norm.
!
!     ATTENZIONE
!         write(*,*) 'ATTENZIONE0',wa2(1:5)
!         write(*,*) 'ATTENZIONE',fnorm,gnorm,gtol,epsmch
!
         if (gnorm .le. gtol) info = 4
         if (info .ne. 0) go to 300
!
!        rescale if necessary.
!
         if (mode .eq. 2) go to 190
         do 180 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  180       continue
  190    continue
!
!        beginning of the inner loop.
!
  200    continue
!
!           determine the levenberg-marquardt parameter.
!
            call lmpar(n,fjac,ldfjac,ipvt,diag,qtf,delta,par,wa1,wa2,   &
     &                 wa3,wa4)
!
!           store the direction p and x + p. calculate the norm of p.
!
            do 210 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  210          continue
            pnorm = enorm(n,wa3)
!
!           on the first iteration, adjust the initial step bound.
!
            if (iter .eq. 1) delta = dmin1(delta,pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
            iflag = 1
            call fcn(m,n,wa2,wa4,iflag)
            nfev = nfev + 1
            if (iflag .lt. 0) go to 300
            fnorm1 = enorm(m,wa4)
!
!           compute the scaled actual reduction.
!
            actred = -one
            if (p1*fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
!
!           compute the scaled predicted reduction and
!           the scaled directional derivative.
!
            do 230 j = 1, n
               wa3(j) = zero
               l = ipvt(j)
               temp = wa1(l)
               do 220 i = 1, j
                  wa3(i) = wa3(i) + fjac(i,j)*temp
  220             continue
  230          continue
            temp1 = enorm(n,wa3)/fnorm
            temp2 = (dsqrt(par)*pnorm)/fnorm
            prered = temp1**2 + temp2**2/p5
            dirder = -(temp1**2 + temp2**2)
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
            ratio = zero
            if (prered .ne. zero) ratio = actred/prered
!
!           update the step bound.
!
            if (ratio .gt. p25) go to 240
               if (actred .ge. zero) temp = p5
               if (actred .lt. zero)                                    &
     &            temp = p5*dirder/(dirder + p5*actred)
               if (p1*fnorm1 .ge. fnorm .or. temp .lt. p1) temp = p1
               delta = temp*dmin1(delta,pnorm/p1)
               par = par/temp
               go to 260
  240       continue
               if (par .ne. zero .and. ratio .lt. p75) go to 250
               delta = pnorm/p5
               par = p5*par
  250          continue
  260       continue
!
!           test for successful iteration.
!
            if (ratio .lt. p0001) go to 290
!
!           successful iteration. update x, fvec, and their norms.
!
            do 270 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
  270          continue
            do 280 i = 1, m
               fvec(i) = wa4(i)
  280          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  290       continue
!
!           tests for convergence.
!
            if (dabs(actred) .le. ftol .and. prered .le. ftol           &
     &          .and. p5*ratio .le. one) info = 1
            if (delta .le. xtol*xnorm) info = 2
            if (dabs(actred) .le. ftol .and. prered .le. ftol           &
     &          .and. p5*ratio .le. one .and. info .eq. 2) info = 3
            if (info .ne. 0) go to 300
!
!           tests for termination and stringent tolerances.
!
            if (nfev .ge. maxfev) info = 5
            if (dabs(actred) .le. epsmch .and. prered .le. epsmch       &
     &          .and. p5*ratio .le. one) info = 6
            if (delta .le. epsmch*xnorm) info = 7
            if (gnorm .le. epsmch) info = 8
            if (info .ne. 0) go to 300
!
!           end of the inner loop. repeat if iteration unsuccessful.
!
            if (ratio .lt. p0001) go to 200
!
!        end of the outer loop.
!
         go to 30
  300 continue
!
!     termination, either normal or user imposed.
!
      if (iflag .lt. 0) info = iflag
      iflag = 0
      if (nprint .gt. 0) call fcn(m,n,x,fvec,iflag)
      return
!
!     last card of subroutine lmdif.
!
      end
