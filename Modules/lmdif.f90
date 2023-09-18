! Minpack Copyright Notice (1999) University of Chicago.  All rights reserved
! 
! Redistribution and use in source and binary forms, with or
! without modification, are permitted provided that the
! following conditions are met:
! 
! 1. Redistributions of source code must retain the above
! copyright notice, this list of conditions and the following
! disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above
! copyright notice, this list of conditions and the following
! disclaimer in the documentation and/or other materials
! provided with the distribution.
! 
! 3. The end-user documentation included with the
! redistribution, if any, must include the following
! acknowledgment:
! 
!    "This product includes software developed by the
!    University of Chicago, as Operator of Argonne National
!    Laboratory.
! 
! Alternately, this acknowledgment may appear in the software
! itself, if and wherever such third-party acknowledgments
! normally appear.
! 
! 4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
! WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
! UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
! THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
! OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
! OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
! USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
! THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
! DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
! UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
! BE CORRECTED.
! 
! 5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
! HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
! ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
! INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
! ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
! PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
! SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
! (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
! EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
! POSSIBILITY OF SUCH LOSS OR DAMAGES.
!
! -- End of copyright notice
!
! Contains lmdif1 method and dependencies, converted to fortran90
! syntax and module interface, September 2020
!
! Notes on usage: 
!  lmdif1 minimizes the sum of squares of m functions of n parameters.
!  For some reason, n cannot be larger than m. This means that if you want to
!  minimize a scalar value (i.e. chi-square), you'll have to pad it with zeros 
!  to form and array of size n.
!    
!
MODULE lmdif_module

   PRIVATE
   PUBLIC :: lmdif
   PUBLIC :: lmdif1 ! easy to use interface, scroll down for documentation
   PUBLIC :: lmdif0 ! easier to use interface

   CONTAINS

   SUBROUTINE lmdif0(fcn, m, n, x, fvec, tol, info)
      INTEGER m, n, info
      DOUBLEPRECISION tol
      DOUBLEPRECISION x (n), fvec (m)
      EXTERNAL fcn
      ! internal variables
      INTEGER ipvt(n), maxfev, mode
      DOUBLEPRECISION qtf(n), fjac(m,n), diag(n)
      DOUBLEPRECISION wa1 (n), wa2(n), wa3(n), wa4(m)
      INTEGER iwa (n), nprint, nfev
      DOUBLEPRECISION factor, epsdiff
      !
      IF(n>m)THEN
            PRINT*, "LMDIF expects n<=m"
            PRINT*, "Hint: give it f_fit-f_real, it will compute chi^2 internally"
            STOP 1
      ENDIF
      !
      diag= 1.d0        ! all the variables have the same importance
      mode = 1          ! set diag automatically
      factor = 1.d0     ! initial step factor
      epsdiff = 0d0     ! precision of fcn (used fo finite difference differentiation)
      nprint = 0
      maxfev = huge(1)  ! take as many iterations as needed

      CALL lmdif(fcn,m,n,x,fvec,tol,tol,0d0,maxfev,epsdiff, &
                 diag,mode,factor,0,info,nfev,fjac,  &
                 m,ipvt,qtf,wa1,wa2,wa3,wa4)

   END SUBROUTINE

      DOUBLEPRECISION function dpmpar (i)
      INTEGER i
!     **********
!
!     Function dpmpar
!
!     This function provides double precision machine parameters
!     when the appropriate set of data statements is activated (by
!     removing the c from column 1) and all other data statements are
!     rendered inactive. Most of the parameter values were obtained
!     from the corresponding Bell Laboratories Port Library function.
!
!     The function statement is
!
!       double precision function dpmpar(i)
!
!     where
!
!       i is an integer input variable set to 1, 2, or 3 which
!         selects the desired machine parameter. If the machine has
!         t base b digits and its smallest and largest exponents are
!         emin and emax, respectively, then these parameters are
!
!         dpmpar(1) = b**(1 - t), the machine precision,
!
!         dpmpar(2) = b**(emin - 1), the smallest magnitude,
!
!         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
!
!     Argonne National Laboratory. MINPACK Project. November 1996.
!     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'
!
!     **********
      INTEGER mcheps (4)
      INTEGER minmag (4)
      INTEGER maxmag (4)
      DOUBLEPRECISION dmach (3)
      EQUIVALENCE (dmach (1), mcheps (1) )
      EQUIVALENCE (dmach (2), minmag (1) )
      EQUIVALENCE (dmach (3), maxmag (1) )
!
!     Machine constants for the IBM 360/370 series,
!     the Amdahl 470/V6, the ICL 2900, the Itel AS/6,
!     the Xerox Sigma 5/7/9 and the Sel systems 85/86.
!
!     data mcheps(1),mcheps(2) / z34100000, z00000000 /
!     data minmag(1),minmag(2) / z00100000, z00000000 /
!     data maxmag(1),maxmag(2) / z7fffffff, zffffffff /
!
!     Machine constants for the Honeywell 600/6000 series.
!
!     data mcheps(1),mcheps(2) / o606400000000, o000000000000 /
!     data minmag(1),minmag(2) / o402400000000, o000000000000 /
!     data maxmag(1),maxmag(2) / o376777777777, o777777777777 /
!
!     Machine constants for the CDC 6000/7000 series.
!
!     data mcheps(1) / 15614000000000000000b /
!     data mcheps(2) / 15010000000000000000b /
!
!     data minmag(1) / 00604000000000000000b /
!     data minmag(2) / 00000000000000000000b /
!
!     data maxmag(1) / 37767777777777777777b /
!     data maxmag(2) / 37167777777777777777b /
!
!     Machine constants for the PDP-10 (KA processor).
!
!     data mcheps(1),mcheps(2) / "114400000000, "000000000000 /
!     data minmag(1),minmag(2) / "033400000000, "000000000000 /
!     data maxmag(1),maxmag(2) / "377777777777, "344777777777 /
!
!     Machine constants for the PDP-10 (KI processor).
!
!     data mcheps(1),mcheps(2) / "104400000000, "000000000000 /
!     data minmag(1),minmag(2) / "000400000000, "000000000000 /
!     data maxmag(1),maxmag(2) / "377777777777, "377777777777 /
!
!     Machine constants for the PDP-11.
!
!     data mcheps(1),mcheps(2) /   9472,      0 /
!     data mcheps(3),mcheps(4) /      0,      0 /
!
!     data minmag(1),minmag(2) /    128,      0 /
!     data minmag(3),minmag(4) /      0,      0 /
!
!     data maxmag(1),maxmag(2) /  32767,     -1 /
!     data maxmag(3),maxmag(4) /     -1,     -1 /
!
!     Machine constants for the Burroughs 6700/7700 systems.
!
!     data mcheps(1) / o1451000000000000 /
!     data mcheps(2) / o0000000000000000 /
!
!     data minmag(1) / o1771000000000000 /
!     data minmag(2) / o7770000000000000 /
!
!     data maxmag(1) / o0777777777777777 /
!     data maxmag(2) / o7777777777777777 /
!
!     Machine constants for the Burroughs 5700 system.
!
!     data mcheps(1) / o1451000000000000 /
!     data mcheps(2) / o0000000000000000 /
!
!     data minmag(1) / o1771000000000000 /
!     data minmag(2) / o0000000000000000 /
!
!     data maxmag(1) / o0777777777777777 /
!     data maxmag(2) / o0007777777777777 /
!
!     Machine constants for the Burroughs 1700 system.
!
!     data mcheps(1) / zcc6800000 /
!     data mcheps(2) / z000000000 /
!
!     data minmag(1) / zc00800000 /
!     data minmag(2) / z000000000 /
!
!     data maxmag(1) / zdffffffff /
!     data maxmag(2) / zfffffffff /
!
!     Machine constants for the Univac 1100 series.
!
!     data mcheps(1),mcheps(2) / o170640000000, o000000000000 /
!     data minmag(1),minmag(2) / o000040000000, o000000000000 /
!     data maxmag(1),maxmag(2) / o377777777777, o777777777777 /
!
!     Machine constants for the Data General Eclipse S/200.
!
!     Note - it may be appropriate to include the following card -
!     static dmach(3)
!
!     data minmag/20k,3*0/,maxmag/77777k,3*177777k/
!     data mcheps/32020k,3*0/
!
!     Machine constants for the Harris 220.
!
!     data mcheps(1),mcheps(2) / '20000000, '00000334 /
!     data minmag(1),minmag(2) / '20000000, '00000201 /
!     data maxmag(1),maxmag(2) / '37777777, '37777577 /
!
!     Machine constants for the Cray-1.
!
!     data mcheps(1) / 0376424000000000000000b /
!     data mcheps(2) / 0000000000000000000000b /
!
!     data minmag(1) / 0200034000000000000000b /
!     data minmag(2) / 0000000000000000000000b /
!
!     data maxmag(1) / 0577777777777777777777b /
!     data maxmag(2) / 0000007777777777777776b /
!
!     Machine constants for the Prime 400.
!
!     data mcheps(1),mcheps(2) / :10000000000, :00000000123 /
!     data minmag(1),minmag(2) / :10000000000, :00000100000 /
!     data maxmag(1),maxmag(2) / :17777777777, :37777677776 /
!
!     Machine constants for the VAX-11.
!
!     data mcheps(1),mcheps(2) /   9472,  0 /
!     data minmag(1),minmag(2) /    128,  0 /
!     data maxmag(1),maxmag(2) / -32769, -1 /
!
!     Machine constants for IEEE machines.
!
      DATA dmach (1) / 2.22044604926d-16 /
      DATA dmach (2) / 2.22507385852d-308 /
      DATA dmach (3) / 1.79769313485d+308 /
!
      dpmpar = dmach (i)
      RETURN
!
!     Last card of function dpmpar.
!
      END FUNCTION dpmpar
      DOUBLEPRECISION function enorm (n, x)
      INTEGER n
      DOUBLEPRECISION x (n)
!     **********
!
!     function enorm
!
!     given an n-vector x, this function calculates the
!     euclidean norm of x.
!
!     the euclidean norm is computed by accumulating the sum of
!     squares in three different sums. the sums of squares for the
!     small and large components are scaled so that no overflows
!     occur. non-destructive underflows are permitted. underflows
!     and overflows do not occur in the computation of the unscaled
!     sum of squares for the intermediate components.
!     the definitions of small, intermediate and large components
!     depend on two constants, rdwarf and rgiant. the main
!     restrictions on these constants are that rdwarf**2 not
!     underflow and rgiant**2 not overflow. the constants
!     given here are suitable for every known computer.
!
!     the function statement is
!
!       double precision function enorm(n,x)
!
!     where
!
!       n is a positive integer input variable.
!
!       x is an input array of length n.
!
!     subprograms called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i
      DOUBLEPRECISION agiant, floatn, one, rdwarf, rgiant, s1, s2, s3,  &
      xabs, x1max, x3max, zero
      DATA one, zero, rdwarf, rgiant / 1.0d0, 0.0d0, 3.834d-20,         &
      1.304d19 /
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = n
      agiant = rgiant / floatn
      DO 90 i = 1, n
         xabs = dabs (x (i) )
         IF (xabs.gt.rdwarf.and.xabs.lt.agiant) goto 70
         IF (xabs.le.rdwarf) goto 30
!
!              sum for large components.
!
         IF (xabs.le.x1max) goto 10
         s1 = one+s1 * (x1max / xabs) **2
         x1max = xabs
         GOTO 20
   10    CONTINUE
         s1 = s1 + (xabs / x1max) **2
   20    CONTINUE
         GOTO 60
   30    CONTINUE
!
!              sum for small components.
!
         IF (xabs.le.x3max) goto 40
         s3 = one+s3 * (x3max / xabs) **2
         x3max = xabs
         GOTO 50
   40    CONTINUE
         IF (xabs.ne.zero) s3 = s3 + (xabs / x3max) **2
   50    CONTINUE
   60    CONTINUE
         GOTO 80
   70    CONTINUE
!
!           sum for intermediate components.
!
         s2 = s2 + xabs**2
   80    CONTINUE
   90 END DO
!
!     calculation of norm.
!
      IF (s1.eq.zero) goto 100
      enorm = x1max * dsqrt (s1 + (s2 / x1max) / x1max)
      GOTO 130
  100 CONTINUE
      IF (s2.eq.zero) goto 110
      IF (s2.ge.x3max) enorm = dsqrt (s2 * (one+ (x3max / s2) * (x3max *&
      s3) ) )
      IF (s2.lt.x3max) enorm = dsqrt (x3max * ( (s2 / x3max) + (x3max * &
      s3) ) )
      GOTO 120
  110 CONTINUE
      enorm = x3max * dsqrt (s3)
  120 CONTINUE
  130 CONTINUE
      RETURN
!
!     last card of function enorm.
!
      END FUNCTION enorm
      SUBROUTINE fdjac2 (fcn, m, n, x, fvec, fjac, ldfjac, iflag,       &
      epsfcn, wa)
      INTEGER m, n, ldfjac, iflag
      DOUBLEPRECISION epsfcn
      DOUBLEPRECISION x (n), fvec (m), fjac (ldfjac, n), wa (m)
!     **********
!
!     subroutine fdjac2
!
!     this subroutine computes a forward-difference approximation
!     to the m by n jacobian matrix associated with a specified
!     problem of m functions in n variables.
!
!     the subroutine statement is
!
!       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
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
!         the user wants to terminate execution of fdjac2.
!         in this case set iflag to a negative integer.
!
!       m is a positive integer input variable set to the number
!         of functions.
!
!       n is a positive integer input variable set to the number
!         of variables. n must not exceed m.
!
!       x is an input array of length n.
!
!       fvec is an input array of length m which must contain the
!         functions evaluated at x.
!
!       fjac is an output m by n array which contains the
!         approximation to the jacobian matrix evaluated at x.
!
!       ldfjac is a positive integer input variable not less than m
!         which specifies the leading dimension of the array fjac.
!
!       iflag is an integer variable which can be used to terminate
!         the execution of fdjac2. see description of fcn.
!
!       epsfcn is an input variable used in determining a suitable
!         step length for the forward-difference approximation. this
!         approximation assumes that the relative errors in the
!         functions are of the order of epsfcn. if epsfcn is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       wa is a work array of length m.
!
!     subprograms called
!
!       user-supplied ...... fcn
!
!       minpack-supplied ... dpmpar
!
!       fortran-supplied ... dabs,dmax1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i, j
      DOUBLEPRECISION eps, epsmch, h, temp, zero
      !DOUBLEPRECISION dpmpar
      DATA zero / 0.0d0 /
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar (1)
!
      eps = dsqrt (dmax1 (epsfcn, epsmch) )
      DO 20 j = 1, n
         temp = x (j)
         h = eps * dabs (temp)
         IF (h.eq.zero) h = eps
         x (j) = temp + h
         CALL fcn (m, n, x, wa, iflag)
         IF (iflag.lt.0) goto 30
         x (j) = temp
         DO 10 i = 1, m
            fjac (i, j) = (wa (i) - fvec (i) ) / h
   10    END DO
   20 END DO
   30 CONTINUE
      RETURN
!
!     last card of subroutine fdjac2.
!
      END SUBROUTINE fdjac2
      SUBROUTINE lmdif (fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev,   &
      epsfcn, diag, mode, factor, nprint, info, nfev, fjac, ldfjac,     &
      ipvt, qtf, wa1, wa2, wa3, wa4)
      INTEGER m, n, maxfev, mode, nprint, info, nfev, ldfjac
      INTEGER ipvt (n)
      DOUBLEPRECISION ftol, xtol, gtol, epsfcn, factor
      DOUBLEPRECISION x (n), fvec (m), diag (n), fjac (ldfjac, n),      &
      qtf (n), wa1 (n), wa2 (n), wa3 (n), wa4 (m)
      EXTERNAL fcn
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
      INTEGER i, iflag, iter, j, l
      DOUBLEPRECISION actred, delta, dirder, epsmch, fnorm, fnorm1,     &
      gnorm, one, par, pnorm, prered, p1, p5, p25, p75, p0001, ratio,   &
      sum, temp, temp1, temp2, xnorm, zero
      !DOUBLEPRECISION dpmpar, enorm
      DATA one, p1, p5, p25, p75, p0001, zero / 1.0d0, 1.0d-1, 5.0d-1,  &
      2.5d-1, 7.5d-1, 1.0d-4, 0.0d0 /
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar (1)
!
      info = 0
      iflag = 0
      nfev = 0
!
!     check the input parameters for errors.
!
      IF (n.le.0.or.m.lt.n.or.ldfjac.lt.m.or.ftol.lt.zero.or.xtol.lt.zer&
     &o.or.gtol.lt.zero.or.maxfev.le.0.or.factor.le.zero) goto 300
      IF (mode.ne.2) goto 20
      DO 10 j = 1, n
         IF (diag (j) .le.zero) goto 300
   10 END DO
   20 CONTINUE
!
!     evaluate the function at the starting point
!     and calculate its norm.
!
      iflag = 1
      CALL fcn (m, n, x, fvec, iflag)
      nfev = 1
      IF (iflag.lt.0) goto 300
      fnorm = enorm (m, fvec)
!
!     initialize levenberg-marquardt parameter and iteration counter.
!
      par = zero
      iter = 1
!
!     beginning of the outer loop.
!
   30 CONTINUE
!
!        calculate the jacobian matrix.
!
      iflag = 2
      CALL fdjac2 (fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn, wa4)
      nfev = nfev + n
      IF (iflag.lt.0) goto 300
!
!        if requested, call fcn to enable printing of iterates.
!
      IF (nprint.le.0) goto 40
      iflag = 0
      IF (mod (iter - 1, nprint) .eq.0) call fcn (m, n, x, fvec, iflag)
      IF (iflag.lt.0) goto 300
   40 CONTINUE
!
!        compute the qr factorization of the jacobian.
!
      CALL qrfac (m, n, fjac, ldfjac, .true., ipvt, n, wa1, wa2, wa3)
!
!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.
!
      IF (iter.ne.1) goto 80
      IF (mode.eq.2) goto 60
      DO 50 j = 1, n
         diag (j) = wa2 (j)
         IF (wa2 (j) .eq.zero) diag (j) = one
   50 END DO
   60 CONTINUE
!
!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.
!
      DO 70 j = 1, n
         wa3 (j) = diag (j) * x (j)
   70 END DO
      xnorm = enorm (n, wa3)
      delta = factor * xnorm
      IF (delta.eq.zero) delta = factor
   80 CONTINUE
!
!        form (q transpose)*fvec and store the first n components in
!        qtf.
!
      DO 90 i = 1, m
         wa4 (i) = fvec (i)
   90 END DO
      DO 130 j = 1, n
         IF (fjac (j, j) .eq.zero) goto 120
         sum = zero
         DO 100 i = j, m
            sum = sum + fjac (i, j) * wa4 (i)
  100    END DO
         temp = - sum / fjac (j, j)
         DO 110 i = j, m
            wa4 (i) = wa4 (i) + fjac (i, j) * temp
  110    END DO
  120    CONTINUE
         fjac (j, j) = wa1 (j)
         qtf (j) = wa4 (j)
  130 END DO
!
!        compute the norm of the scaled gradient.
!
      gnorm = zero
      IF (fnorm.eq.zero) goto 170
      DO 160 j = 1, n
         l = ipvt (j)
         IF (wa2 (l) .eq.zero) goto 150
         sum = zero
         DO 140 i = 1, j
            sum = sum + fjac (i, j) * (qtf (i) / fnorm)
  140    END DO
         gnorm = dmax1 (gnorm, dabs (sum / wa2 (l) ) )
  150    CONTINUE
  160 END DO
  170 CONTINUE
!
!        test for convergence of the gradient norm.
!
      IF (gnorm.le.gtol) info = 4
      IF (info.ne.0) goto 300
!
!        rescale if necessary.
!
      IF (mode.eq.2) goto 190
      DO 180 j = 1, n
         diag (j) = dmax1 (diag (j), wa2 (j) )
  180 END DO
  190 CONTINUE
!
!        beginning of the inner loop.
!
  200 CONTINUE
!
!           determine the levenberg-marquardt parameter.
!
      CALL lmpar (n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1,    &
      wa2, wa3, wa4)
!
!           store the direction p and x + p. calculate the norm of p.
!
      DO 210 j = 1, n
         wa1 (j) = - wa1 (j)
         wa2 (j) = x (j) + wa1 (j)
         wa3 (j) = diag (j) * wa1 (j)
  210 END DO
      pnorm = enorm (n, wa3)
!
!           on the first iteration, adjust the initial step bound.
!
      IF (iter.eq.1) delta = dmin1 (delta, pnorm)
!
!           evaluate the function at x + p and calculate its norm.
!
      iflag = 1
      CALL fcn (m, n, wa2, wa4, iflag)
      nfev = nfev + 1
      IF (iflag.lt.0) goto 300
      fnorm1 = enorm (m, wa4)
!
!           compute the scaled actual reduction.
!
      actred = - one
      IF (p1 * fnorm1.lt.fnorm) actred = one- (fnorm1 / fnorm) **2
!
!           compute the scaled predicted reduction and
!           the scaled directional derivative.
!
      DO 230 j = 1, n
         wa3 (j) = zero
         l = ipvt (j)
         temp = wa1 (l)
         DO 220 i = 1, j
            wa3 (i) = wa3 (i) + fjac (i, j) * temp
  220    END DO
  230 END DO
      temp1 = enorm (n, wa3) / fnorm
      temp2 = (dsqrt (par) * pnorm) / fnorm
      prered = temp1**2 + temp2**2 / p5
      dirder = - (temp1**2 + temp2**2)
!
!           compute the ratio of the actual to the predicted
!           reduction.
!
      ratio = zero
      IF (prered.ne.zero) ratio = actred / prered
!
!           update the step bound.
!
      IF (ratio.gt.p25) goto 240
      IF (actred.ge.zero) temp = p5
      IF (actred.lt.zero) temp = p5 * dirder / (dirder + p5 * actred)
      IF (p1 * fnorm1.ge.fnorm.or.temp.lt.p1) temp = p1
      delta = temp * dmin1 (delta, pnorm / p1)
      par = par / temp
      GOTO 260
  240 CONTINUE
      IF (par.ne.zero.and.ratio.lt.p75) goto 250
      delta = pnorm / p5
      par = p5 * par
  250 CONTINUE
  260 CONTINUE
!
!           test for successful iteration.
!
      IF (ratio.lt.p0001) goto 290
!
!           successful iteration. update x, fvec, and their norms.
!
      DO 270 j = 1, n
         x (j) = wa2 (j)
         wa2 (j) = diag (j) * x (j)
  270 END DO
      DO 280 i = 1, m
         fvec (i) = wa4 (i)
  280 END DO
      xnorm = enorm (n, wa2)
      fnorm = fnorm1
      iter = iter + 1
  290 CONTINUE
!
!           tests for convergence.
!
      IF (dabs (actred) .le.ftol.and.prered.le.ftol.and.p5 *            &
      ratio.le.one) info = 1
      IF (delta.le.xtol * xnorm) info = 2
      IF (dabs (actred) .le.ftol.and.prered.le.ftol.and.p5 *            &
      ratio.le.one.and.info.eq.2) info = 3
      IF (info.ne.0) goto 300
!
!           tests for termination and stringent tolerances.
!
      IF (nfev.ge.maxfev) info = 5
      IF (dabs (actred) .le.epsmch.and.prered.le.epsmch.and.p5 *        &
      ratio.le.one) info = 6
      IF (delta.le.epsmch * xnorm) info = 7
      IF (gnorm.le.epsmch) info = 8
      IF (info.ne.0) goto 300
!
!           end of the inner loop. repeat if iteration unsuccessful.
!
      IF (ratio.lt.p0001) goto 200
!
!        end of the outer loop.
!
      GOTO 30
  300 CONTINUE
!
!     termination, either normal or user imposed.
!
      IF (iflag.lt.0) info = iflag
      iflag = 0
      IF (nprint.gt.0) call fcn (m, n, x, fvec, iflag)
      RETURN
!
!     last card of subroutine lmdif.
!
      END SUBROUTINE lmdif

! =====================================================================

      SUBROUTINE lmdif1 (fcn, m, n, x, fvec, tol, info, iwa, wa, lwa)
      INTEGER m, n, info, lwa
      INTEGER iwa (n)
      DOUBLEPRECISION tol
      DOUBLEPRECISION x (n), fvec (m), wa (lwa)
      EXTERNAL fcn
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
      INTEGER maxfev, mode, mp5n, nfev, nprint
      DOUBLEPRECISION epsfcn, factor, ftol, gtol, xtol, zero
      DATA factor, zero / 1.0d2, 0.0d0 /
      info = 0
!
!     check the input parameters for errors.
!
      IF (n.le.0.or.m.lt.n.or.tol.lt.zero.or.lwa.lt.m * n + 5 * n + m)  &
      goto 10
!
!     call lmdif.
!
      maxfev = 5000 * (n + 1)
      ftol = tol
      xtol = tol
      gtol = zero
      epsfcn = zero
      mode = 1
      nprint = 0
      mp5n = m + 5 * n
      CALL lmdif (fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn, &
      wa (1), mode, factor, nprint, info, nfev, wa (mp5n + 1), m, iwa,  &
      wa (n + 1), wa (2 * n + 1), wa (3 * n + 1), wa (4 * n + 1),       &
      wa (5 * n + 1) )
      IF (info.eq.8) info = 4
   10 CONTINUE
      RETURN
!
!     last card of subroutine lmdif1.
!
      END SUBROUTINE lmdif1
      SUBROUTINE lmpar (n, r, ldr, ipvt, diag, qtb, delta, par, x,      &
      sdiag, wa1, wa2)
      INTEGER n, ldr
      INTEGER ipvt (n)
      DOUBLEPRECISION delta, par
      DOUBLEPRECISION r (ldr, n), diag (n), qtb (n), x (n), sdiag (n),  &
      wa1 (n), wa2 (n)
!     **********
!
!     subroutine lmpar
!
!     given an m by n matrix a, an n by n nonsingular diagonal
!     matrix d, an m-vector b, and a positive number delta,
!     the problem is to determine a value for the parameter
!     par such that if x solves the system
!
!           a*x = b ,     sqrt(par)*d*x = 0 ,
!
!     in the least squares sense, and dxnorm is the euclidean
!     norm of d*x, then either par is zero and
!
!           (dxnorm-delta) .le. 0.1*delta ,
!
!     or par is positive and
!
!           abs(dxnorm-delta) .le. 0.1*delta .
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization, with column pivoting, of a. that is, if
!     a*p = q*r, where p is a permutation matrix, q has orthogonal
!     columns, and r is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then lmpar expects
!     the full upper triangle of r, the permutation matrix p,
!     and the first n components of (q transpose)*b. on output
!     lmpar also provides an upper triangular matrix s such that
!
!            t   t                   t
!           p *(a *a + par*d*d)*p = s *s .
!
!     s is employed within lmpar and may be of separate interest.
!
!     only a few iterations are generally needed for convergence
!     of the algorithm. if, however, the limit of 10 iterations
!     is reached, then the output par will contain the best
!     value obtained so far.
!
!     the subroutine statement is
!
!       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
!                        wa1,wa2)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. on input the full upper triangle
!         must contain the full upper triangle of the matrix r.
!         on output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix s.
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       ipvt is an integer input array of length n which defines the
!         permutation matrix p such that a*p = q*r. column j of p
!         is column ipvt(j) of the identity matrix.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.
!
!       delta is a positive input variable which specifies an upper
!         bound on the euclidean norm of d*x.
!
!       par is a nonnegative variable. on input par contains an
!         initial estimate of the levenberg-marquardt parameter.
!         on output par contains the final estimate.
!
!       x is an output array of length n which contains the least
!         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
!         for the output par.
!
!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix s.
!
!       wa1 and wa2 are work arrays of length n.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar,enorm,qrsolv
!
!       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i, iter, j, jm1, jp1, k, l, nsing
      DOUBLEPRECISION dxnorm, dwarf, fp, gnorm, parc, parl, paru, p1,   &
      p001, sum, temp, zero
      !DOUBLEPRECISION dpmpar, enorm
      DATA p1, p001, zero / 1.0d-1, 1.0d-3, 0.0d0 /
!
!     dwarf is the smallest positive magnitude.
!
      dwarf = dpmpar (2)
!
!     compute and store in x the gauss-newton direction. if the
!     jacobian is rank-deficient, obtain a least squares solution.
!
      nsing = n
      DO 10 j = 1, n
         wa1 (j) = qtb (j)
         IF (r (j, j) .eq.zero.and.nsing.eq.n) nsing = j - 1
         IF (nsing.lt.n) wa1 (j) = zero
   10 END DO
      IF (nsing.lt.1) goto 50
      DO 40 k = 1, nsing
         j = nsing - k + 1
         wa1 (j) = wa1 (j) / r (j, j)
         temp = wa1 (j)
         jm1 = j - 1
         IF (jm1.lt.1) goto 30
         DO 20 i = 1, jm1
            wa1 (i) = wa1 (i) - r (i, j) * temp
   20    END DO
   30    CONTINUE
   40 END DO
   50 CONTINUE
      DO 60 j = 1, n
         l = ipvt (j)
         x (l) = wa1 (j)
   60 END DO
!
!     initialize the iteration counter.
!     evaluate the function at the origin, and test
!     for acceptance of the gauss-newton direction.
!
      iter = 0
      DO 70 j = 1, n
         wa2 (j) = diag (j) * x (j)
   70 END DO
      dxnorm = enorm (n, wa2)
      fp = dxnorm - delta
      IF (fp.le.p1 * delta) goto 220
!
!     if the jacobian is not rank deficient, the newton
!     step provides a lower bound, parl, for the zero of
!     the function. otherwise set this bound to zero.
!
      parl = zero
      IF (nsing.lt.n) goto 120
      DO 80 j = 1, n
         l = ipvt (j)
         wa1 (j) = diag (l) * (wa2 (l) / dxnorm)
   80 END DO
      DO 110 j = 1, n
         sum = zero
         jm1 = j - 1
         IF (jm1.lt.1) goto 100
         DO 90 i = 1, jm1
            sum = sum + r (i, j) * wa1 (i)
   90    END DO
  100    CONTINUE
         wa1 (j) = (wa1 (j) - sum) / r (j, j)
  110 END DO
      temp = enorm (n, wa1)
      parl = ( (fp / delta) / temp) / temp
  120 CONTINUE
!
!     calculate an upper bound, paru, for the zero of the function.
!
      DO 140 j = 1, n
         sum = zero
         DO 130 i = 1, j
            sum = sum + r (i, j) * qtb (i)
  130    END DO
         l = ipvt (j)
         wa1 (j) = sum / diag (l)
  140 END DO
      gnorm = enorm (n, wa1)
      paru = gnorm / delta
      IF (paru.eq.zero) paru = dwarf / dmin1 (delta, p1)
!
!     if the input par lies outside of the interval (parl,paru),
!     set par to the closer endpoint.
!
      par = dmax1 (par, parl)
      par = dmin1 (par, paru)
      IF (par.eq.zero) par = gnorm / dxnorm
!
!     beginning of an iteration.
!
  150 CONTINUE
      iter = iter + 1
!
!        evaluate the function at the current value of par.
!
      IF (par.eq.zero) par = dmax1 (dwarf, p001 * paru)
      temp = dsqrt (par)
      DO 160 j = 1, n
         wa1 (j) = temp * diag (j)
  160 END DO
      CALL qrsolv (n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2)
      DO 170 j = 1, n
         wa2 (j) = diag (j) * x (j)
  170 END DO
      dxnorm = enorm (n, wa2)
      temp = fp
      fp = dxnorm - delta
!
!        if the function is small enough, accept the current value
!        of par. also test for the exceptional cases where parl
!        is zero or the number of iterations has reached 10.
!
      IF (dabs (fp) .le.p1 * delta.or.parl.eq.zero.and.fp.le.temp.and.te&
     &mp.lt.zero.or.iter.eq.10) goto 220
!
!        compute the newton correction.
!
      DO 180 j = 1, n
         l = ipvt (j)
         wa1 (j) = diag (l) * (wa2 (l) / dxnorm)
  180 END DO
      DO 210 j = 1, n
         wa1 (j) = wa1 (j) / sdiag (j)
         temp = wa1 (j)
         jp1 = j + 1
         IF (n.lt.jp1) goto 200
         DO 190 i = jp1, n
            wa1 (i) = wa1 (i) - r (i, j) * temp
  190    END DO
  200    CONTINUE
  210 END DO
      temp = enorm (n, wa1)
      parc = ( (fp / delta) / temp) / temp
!
!        depending on the sign of the function, update parl or paru.
!
      IF (fp.gt.zero) parl = dmax1 (parl, par)
      IF (fp.lt.zero) paru = dmin1 (paru, par)
!
!        compute an improved estimate for par.
!
      par = dmax1 (parl, par + parc)
!
!        end of an iteration.
!
      GOTO 150
  220 CONTINUE
!
!     termination.
!
      IF (iter.eq.0) par = zero
      RETURN
!
!     last card of subroutine lmpar.
!
      END SUBROUTINE lmpar
      SUBROUTINE qrfac (m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm,&
      wa)
      INTEGER m, n, lda, lipvt
      INTEGER ipvt (lipvt)
      LOGICAL pivot
      DOUBLEPRECISION a (lda, n), rdiag (n), acnorm (n), wa (n)
!     **********
!
!     subroutine qrfac
!
!     this subroutine uses householder transformations with column
!     pivoting (optional) to compute a qr factorization of the
!     m by n matrix a. that is, qrfac determines an orthogonal
!     matrix q, a permutation matrix p, and an upper trapezoidal
!     matrix r with diagonal elements of nonincreasing magnitude,
!     such that a*p = q*r. the householder transformation for
!     column k, k = 1,2,...,min(m,n), is of the form
!
!                           t
!           i - (1/u(k))*u*u
!
!     where u has zeros in the first k-1 positions. the form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding linpack subroutine.
!
!     the subroutine statement is
!
!       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
!
!     where
!
!       m is a positive integer input variable set to the number
!         of rows of a.
!
!       n is a positive integer input variable set to the number
!         of columns of a.
!
!       a is an m by n array. on input a contains the matrix for
!         which the qr factorization is to be computed. on output
!         the strict upper trapezoidal part of a contains the strict
!         upper trapezoidal part of r, and the lower trapezoidal
!         part of a contains a factored form of q (the non-trivial
!         elements of the u vectors described above).
!
!       lda is a positive integer input variable not less than m
!         which specifies the leading dimension of the array a.
!
!       pivot is a logical input variable. if pivot is set true,
!         then column pivoting is enforced. if pivot is set false,
!         then no column pivoting is done.
!
!       ipvt is an integer output array of length lipvt. ipvt
!         defines the permutation matrix p such that a*p = q*r.
!         column j of p is column ipvt(j) of the identity matrix.
!         if pivot is false, ipvt is not referenced.
!
!       lipvt is a positive integer input variable. if pivot is false,
!         then lipvt may be as small as 1. if pivot is true, then
!         lipvt must be at least n.
!
!       rdiag is an output array of length n which contains the
!         diagonal elements of r.
!
!       acnorm is an output array of length n which contains the
!         norms of the corresponding columns of the input matrix a.
!         if this information is not needed, then acnorm can coincide
!         with rdiag.
!
!       wa is a work array of length n. if pivot is false, then wa
!         can coincide with rdiag.
!
!     subprograms called
!
!       minpack-supplied ... dpmpar,enorm
!
!       fortran-supplied ... dmax1,dsqrt,min0
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i, j, jp1, k, kmax, minmn
      DOUBLEPRECISION ajnorm, epsmch, one, p05, sum, temp, zero
      !DOUBLEPRECISION dpmpar, enorm
      DATA one, p05, zero / 1.0d0, 5.0d-2, 0.0d0 /
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar (1)
!
!     compute the initial column norms and initialize several arrays.
!
      DO 10 j = 1, n
         acnorm (j) = enorm (m, a (1, j) )
         rdiag (j) = acnorm (j)
         wa (j) = rdiag (j)
         IF (pivot) ipvt (j) = j
   10 END DO
!
!     reduce a to r with householder transformations.
!
      minmn = min0 (m, n)
      DO 110 j = 1, minmn
         IF (.not.pivot) goto 40
!
!        bring the column of largest norm into the pivot position.
!
         kmax = j
         DO 20 k = j, n
            IF (rdiag (k) .gt.rdiag (kmax) ) kmax = k
   20    END DO
         IF (kmax.eq.j) goto 40
         DO 30 i = 1, m
            temp = a (i, j)
            a (i, j) = a (i, kmax)
            a (i, kmax) = temp
   30    END DO
         rdiag (kmax) = rdiag (j)
         wa (kmax) = wa (j)
         k = ipvt (j)
         ipvt (j) = ipvt (kmax)
         ipvt (kmax) = k
   40    CONTINUE
!
!        compute the householder transformation to reduce the
!        j-th column of a to a multiple of the j-th unit vector.
!
         ajnorm = enorm (m - j + 1, a (j, j) )
         IF (ajnorm.eq.zero) goto 100
         IF (a (j, j) .lt.zero) ajnorm = - ajnorm
         DO 50 i = j, m
            a (i, j) = a (i, j) / ajnorm
   50    END DO
         a (j, j) = a (j, j) + one
!
!        apply the transformation to the remaining columns
!        and update the norms.
!
         jp1 = j + 1
         IF (n.lt.jp1) goto 100
         DO 90 k = jp1, n
            sum = zero
            DO 60 i = j, m
               sum = sum + a (i, j) * a (i, k)
   60       END DO
            temp = sum / a (j, j)
            DO 70 i = j, m
               a (i, k) = a (i, k) - temp * a (i, j)
   70       END DO
            IF (.not.pivot.or.rdiag (k) .eq.zero) goto 80
            temp = a (j, k) / rdiag (k)
            rdiag (k) = rdiag (k) * dsqrt (dmax1 (zero, one-temp**2) )
            IF (p05 * (rdiag (k) / wa (k) ) **2.gt.epsmch) goto 80
            rdiag (k) = enorm (m - j, a (jp1, k) )
            wa (k) = rdiag (k)
   80       CONTINUE
   90    END DO
  100    CONTINUE
         rdiag (j) = - ajnorm
  110 END DO
      RETURN
!
!     last card of subroutine qrfac.
!
      END SUBROUTINE qrfac
      SUBROUTINE qrsolv (n, r, ldr, ipvt, diag, qtb, x, sdiag, wa)
      INTEGER n, ldr
      INTEGER ipvt (n)
      DOUBLEPRECISION r (ldr, n), diag (n), qtb (n), x (n), sdiag (n),  &
      wa (n)
!     **********
!
!     subroutine qrsolv
!
!     given an m by n matrix a, an n by n diagonal matrix d,
!     and an m-vector b, the problem is to determine an x which
!     solves the system
!
!           a*x = b ,     d*x = 0 ,
!
!     in the least squares sense.
!
!     this subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     qr factorization, with column pivoting, of a. that is, if
!     a*p = q*r, where p is a permutation matrix, q has orthogonal
!     columns, and r is an upper triangular matrix with diagonal
!     elements of nonincreasing magnitude, then qrsolv expects
!     the full upper triangle of r, the permutation matrix p,
!     and the first n components of (q transpose)*b. the system
!     a*x = b, d*x = 0, is then equivalent to
!
!                  t       t
!           r*z = q *b ,  p *d*p*z = 0 ,
!
!     where x = p*z. if this system does not have full rank,
!     then a least squares solution is obtained. on output qrsolv
!     also provides an upper triangular matrix s such that
!
!            t   t               t
!           p *(a *a + d*d)*p = s *s .
!
!     s is computed within qrsolv and may be of separate interest.
!
!     the subroutine statement is
!
!       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
!
!     where
!
!       n is a positive integer input variable set to the order of r.
!
!       r is an n by n array. on input the full upper triangle
!         must contain the full upper triangle of the matrix r.
!         on output the full upper triangle is unaltered, and the
!         strict lower triangle contains the strict upper triangle
!         (transposed) of the upper triangular matrix s.
!
!       ldr is a positive integer input variable not less than n
!         which specifies the leading dimension of the array r.
!
!       ipvt is an integer input array of length n which defines the
!         permutation matrix p such that a*p = q*r. column j of p
!         is column ipvt(j) of the identity matrix.
!
!       diag is an input array of length n which must contain the
!         diagonal elements of the matrix d.
!
!       qtb is an input array of length n which must contain the first
!         n elements of the vector (q transpose)*b.
!
!       x is an output array of length n which contains the least
!         squares solution of the system a*x = b, d*x = 0.
!
!       sdiag is an output array of length n which contains the
!         diagonal elements of the upper triangular matrix s.
!
!       wa is a work array of length n.
!
!     subprograms called
!
!       fortran-supplied ... dabs,dsqrt
!
!     argonne national laboratory. minpack project. march 1980.
!     burton s. garbow, kenneth e. hillstrom, jorge j. more
!
!     **********
      INTEGER i, j, jp1, k, kp1, l, nsing
      DOUBLEPRECISION cos, cotan, p5, p25, qtbpj, sin, sum, tan, temp,  &
      zero
      DATA p5, p25, zero / 5.0d-1, 2.5d-1, 0.0d0 /
!
!     copy r and (q transpose)*b to preserve input and initialize s.
!     in particular, save the diagonal elements of r in x.
!
      DO 20 j = 1, n
         DO 10 i = j, n
            r (i, j) = r (j, i)
   10    END DO
         x (j) = r (j, j)
         wa (j) = qtb (j)
   20 END DO
!
!     eliminate the diagonal matrix d using a givens rotation.
!
      DO 100 j = 1, n
!
!        prepare the row of d to be eliminated, locating the
!        diagonal element using p from the qr factorization.
!
         l = ipvt (j)
         IF (diag (l) .eq.zero) goto 90
         DO 30 k = j, n
            sdiag (k) = zero
   30    END DO
         sdiag (j) = diag (l)
!
!        the transformations to eliminate the row of d
!        modify only a single element of (q transpose)*b
!        beyond the first n, which is initially zero.
!
         qtbpj = zero
         DO 80 k = j, n
!
!           determine a givens rotation which eliminates the
!           appropriate element in the current row of d.
!
            IF (sdiag (k) .eq.zero) goto 70
            IF (dabs (r (k, k) ) .ge.dabs (sdiag (k) ) ) goto 40
            cotan = r (k, k) / sdiag (k)
            sin = p5 / dsqrt (p25 + p25 * cotan**2)
            cos = sin * cotan
            GOTO 50
   40       CONTINUE
            tan = sdiag (k) / r (k, k)
            cos = p5 / dsqrt (p25 + p25 * tan**2)
            sin = cos * tan
   50       CONTINUE
!
!           compute the modified diagonal element of r and
!           the modified element of ((q transpose)*b,0).
!
            r (k, k) = cos * r (k, k) + sin * sdiag (k)
            temp = cos * wa (k) + sin * qtbpj
            qtbpj = - sin * wa (k) + cos * qtbpj
            wa (k) = temp
!
!           accumulate the tranformation in the row of s.
!
            kp1 = k + 1
            IF (n.lt.kp1) goto 70
            DO 60 i = kp1, n
               temp = cos * r (i, k) + sin * sdiag (i)
               sdiag (i) = - sin * r (i, k) + cos * sdiag (i)
               r (i, k) = temp
   60       END DO
   70       CONTINUE
   80    END DO
   90    CONTINUE
!
!        store the diagonal element of s and restore
!        the corresponding diagonal element of r.
!
         sdiag (j) = r (j, j)
         r (j, j) = x (j)
  100 END DO
!
!     solve the triangular system for z. if the system is
!     singular, then obtain a least squares solution.
!
      nsing = n
      DO 110 j = 1, n
         IF (sdiag (j) .eq.zero.and.nsing.eq.n) nsing = j - 1
         IF (nsing.lt.n) wa (j) = zero
  110 END DO
      IF (nsing.lt.1) goto 150
      DO 140 k = 1, nsing
         j = nsing - k + 1
         sum = zero
         jp1 = j + 1
         IF (nsing.lt.jp1) goto 130
         DO 120 i = jp1, nsing
            sum = sum + r (i, j) * wa (i)
  120    END DO
  130    CONTINUE
         wa (j) = (wa (j) - sum) / sdiag (j)
  140 END DO
  150 CONTINUE
!
!     permute the components of z back to components of x.
!
      DO 160 j = 1, n
         l = ipvt (j)
         x (l) = wa (j)
  160 END DO
      RETURN
!
!     last card of subroutine qrsolv.
!
      END SUBROUTINE qrsolv


END MODULE
!
