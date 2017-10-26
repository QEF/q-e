      subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
      integer m,n,lda,lipvt
      integer ipvt(lipvt)
      logical pivot
      double precision a(lda,n),rdiag(n),acnorm(n),wa(n)
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
      integer i,j,jp1,k,kmax,minmn
      double precision ajnorm,epsmch,one,p05,sum,temp,zero
      double precision dpmpar,enorm
      data one,p05,zero /1.0d0,5.0d-2,0.0d0/
!
!     epsmch is the machine precision.
!
      epsmch = dpmpar(1)
!
!     compute the initial column norms and initialize several arrays.
!
      do 10 j = 1, n
         acnorm(j) = enorm(m,a(1,j))
         rdiag(j) = acnorm(j)
         wa(j) = rdiag(j)
         if (pivot) ipvt(j) = j
   10    continue
!
!     reduce a to r with householder transformations.
!
      minmn = min0(m,n)
      do 110 j = 1, minmn
         if (.not.pivot) go to 40
!
!        bring the column of largest norm into the pivot position.
!
         kmax = j
         do 20 k = j, n
            if (rdiag(k) .gt. rdiag(kmax)) kmax = k
   20       continue
         if (kmax .eq. j) go to 40
         do 30 i = 1, m
            temp = a(i,j)
            a(i,j) = a(i,kmax)
            a(i,kmax) = temp
   30       continue
         rdiag(kmax) = rdiag(j)
         wa(kmax) = wa(j)
         k = ipvt(j)
         ipvt(j) = ipvt(kmax)
         ipvt(kmax) = k
   40    continue
!
!        compute the householder transformation to reduce the
!        j-th column of a to a multiple of the j-th unit vector.
!
         ajnorm = enorm(m-j+1,a(j,j))
         if (ajnorm .eq. zero) go to 100
         if (a(j,j) .lt. zero) ajnorm = -ajnorm
         do 50 i = j, m
            a(i,j) = a(i,j)/ajnorm
   50       continue
         a(j,j) = a(j,j) + one
!
!        apply the transformation to the remaining columns
!        and update the norms.
!
         jp1 = j + 1
         if (n .lt. jp1) go to 100
         do 90 k = jp1, n
            sum = zero
            do 60 i = j, m
               sum = sum + a(i,j)*a(i,k)
   60          continue
            temp = sum/a(j,j)
            do 70 i = j, m
               a(i,k) = a(i,k) - temp*a(i,j)
   70          continue
            if (.not.pivot .or. rdiag(k) .eq. zero) go to 80
            temp = a(j,k)/rdiag(k)
            rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
            if (p05*(rdiag(k)/wa(k))**2 .gt. epsmch) go to 80
            rdiag(k) = enorm(m-j,a(jp1,k))
            wa(k) = rdiag(k)
   80       continue
   90       continue
  100    continue
         rdiag(j) = -ajnorm
  110    continue
      return
!
!     last card of subroutine qrfac.
!
      end
