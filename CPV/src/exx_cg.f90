SUBROUTINE hpotcg(np_in_sp_me, n, rho, pot, fullgrid, mvstep)
    !=======================================================================================
    ! Code Version 1.0 (Princeton University, September 2014)
    !=======================================================================================
    !-------------------------------------------------------------------------------
    !HPOTCG -- using Conjugate Gradient method to compute the Hartree Potential.
    ! Modified from the corresponding subroutine in PARSEC, see http://parsec.ices.utexas.edu/
    ! Lingzhu Kong
    !-------------------------------------------------------------------------------
    !
    USE kinds,                   ONLY  :  DP
    USE wannier_base,            ONLY  : poisson_eps
    USE parallel_include
    !
    IMPLICIT REAL(DP) (a-h,o-z)
    !
    INTEGER  np_in_sp_me,n
    LOGICAL fullgrid
    REAL(DP)  rho(n), pot(n)
    REAL(DP), ALLOCATABLE  :: wk(:), wk2(:)
    !
    INTEGER  i, iou, ipar(16), lwk, mvstep
    REAL*8   fpar(16), dnrm2
    EXTERNAL dnrm2
    !
    lwk = n*5
    ALLOCATE( wk(lwk) )
    ALLOCATE( wk2(np_in_sp_me) )
    !
    !     set up the parameter arrays
    !
    ipar(1) = 0
    ipar(2) = 0
    ipar(3) = 1
    ipar(4) = lwk
    ipar(5) = 5
    ipar(6) = 500
    fpar(1) = poisson_eps
    fpar(2) = poisson_eps
    fpar(11) = 0.0D0
    !
    iou = 6
    mvstep = 0
    10   CALL cg(n,rho,pot,ipar,fpar,wk)
    !
    IF (ipar(1).eq.1) THEN
      !
      CALL dcopy(n, wk(ipar(8)), 1, wk2, 1)
      CALL start_clock('lapmv')
      IF(fullgrid) THEN
        !
        CALL lapmvs(np_in_sp_me, n, wk2, wk(ipar(9)) )
        !
      END IF
      CALL stop_clock('lapmv')
      mvstep = mvstep + 1
      fpar(11) = fpar(11) + 74*n
      !write (iou, *) '# ipar(7) = ', ipar(7),  '   fpar(5) = ', fpar(5), '   cgstep = ', mvstep
      GOTO 10
      !
    ELSE IF (ipar(1).LE.0) THEN
      !  if (ipar(1).eq.0) then
      !     print *, 'Iterative sovler has satisfied convergence test.'
      IF (ipar(1).EQ.-1) THEN
        !
        PRINT *, 'Iterative potver has iterated too many times.'
        !
      ELSE IF (ipar(1).EQ.-2) THEN
        !
        PRINT *, 'Iterative potver was not given enough work space.'
        PRINT *, 'The work space should at least have ', ipar(4),  &
            &           ' elements.'
        !
      ELSE IF (ipar(1).EQ.-3) THEN
        !
        PRINT *, 'Iterative sovler is facing a break-down.'
        PRINT *, 'ipar(12) =', ipar(12)
        !   else
        !      print *, 'Iterative potver terminated. code =', ipar(1)
      ENDIF
      !
    END IF
    !     from_scratchwrite (iou, *) ipar(7), DBLE(fpar(6))
    !     write (iou, *) '# ',
    !    +     ipar(7), ' MATVECs   ', DBLE(fpar(11)), ' OPS'
    !      write (iou, *) '# return code = ', ipar(1),  '   cgstep = ', mvstep
    !      write (iou, *) 'fpar follows'
    !      write (iou, *) (fpar(i),i=1,7)
    !
    !     check the error
    !     CALL lapmvs(pot, wk)
    !     do i = 1, n
    !        wk(i) = wk(i) - rho(i)
    !     enddo
    !     write (iou, *) '# the residual norm ', DBLE(dnrm2(n,wk,1)),
    !    +     DBLE(fpar(5))
    !
    !
    DEALLOCATE( wk, wk2 )
    !
    RETURN
END SUBROUTINE hpotcg
!-----end-of-hpotcg

!-----------------------------------------------------------------------
!FUNCTION distdot(n,x,ix,y,iy)
!    !
!    IMPLICIT NONE
!    INTEGER n, ix, iy
!    REAL*8 distdot, x(*), y(*), ddot
!    EXTERNAL ddot
!    !
!    distdot = ddot(n,x,ix,y,iy)
!    !
!    RETURN
!    !
!END FUNCTION distdot
!-----end-of-distdot

FUNCTION distdot(n,x,ix,y,iy)
    !
    IMPLICIT NONE
    INTEGER n, ix, iy, i
    REAL*8 distdot, x(*), y(*), dtemp
    !
    distdot=0.0d0
    dtemp=0.0d0
    ! 
    !$omp parallel do reduction(+:dtemp)
    DO i = 1,n
      dtemp = dtemp + x(i)*y(i)
    END DO
    !$omp end parallel do
    !
    distdot = dtemp
    ! 
    RETURN
    !
END FUNCTION distdot

!-----------------------------------------------------------------------------
! The matrix-vector multiplication routine for HPOTCG(given p, return q)
! A big assumption is made here: the finite difference neighbors of the
! point in the sphere is still inside the total box 
!-----------------------------------------------------------------------
SUBROUTINE lapmvs(np_in_sp_me, n,p,q)
    !
    USE kinds,            ONLY  :  DP
    USE exx_module,       ONLY  :  odtothd_in_sp,  thdtood_in_sp
    USE exx_module,       ONLY  :  coeke, nord2
    USE constants,        ONLY  :  eps12
    !
    IMPLICIT NONE
    !
    INTEGER np_in_sp_me,n
    REAL(DP)  p(np_in_sp_me), q(n)
    INTEGER   i, ish, ii, jj, kk
    REAL(DP)  tmp, p1, p2, p3
    !
    ! set wave function outside the domain be zero
    !$omp parallel do 
    DO i = n+1, np_in_sp_me
      !
      p(i) = 0.d0 
      !
    ENDDO
    !$omp end parallel do 
    !
    !  --------------------------------------------------------    
    !  diagonal part 
    tmp = coeke(0,1,1) + coeke(0,2,2) + coeke(0,3,3)
    !$omp parallel do 
    DO i = 1, n
      !
      q(i) = tmp * p(i)
      !
    END DO     
    !$omp end parallel do 
    !
    !  kinetic energy part
    !
    !$omp parallel do private(ii,jj,kk,p1,p2,p3) 
    DO i = 1, n
      !
      ii = odtothd_in_sp(1,i)
      jj = odtothd_in_sp(2,i)        
      kk = odtothd_in_sp(3,i)
      !
      DO ish = 1, nord2
        !
        p1 = p( thdtood_in_sp( ii-ish, jj,     kk))    + &
            p( thdtood_in_sp( ii+ish, jj,     kk))
        !
        p2 = p( thdtood_in_sp( ii,     jj-ish, kk))    + &
            p( thdtood_in_sp( ii,     jj+ish, kk))
        !
        p3 = p( thdtood_in_sp( ii,     jj,     kk-ish))+ &
            p( thdtood_in_sp( ii,     jj,     kk+ish))
        !
        q(i) = q(i)+coeke(ish,1,1)*p1+coeke(ish,2,2)*p2+coeke(ish,3,3)*p3 ! stencil on axes
        !
      END DO
      !
    END DO  
    !$omp end parallel do 
    !
    ! cross derivatives
    !
    if (abs(coeke(1,1,2)).gt.eps12) then
      !$omp parallel do private(ii,jj,kk,p1,p2,p3)
      DO i = 1, n
        !
        ii = odtothd_in_sp(1,i)
        jj = odtothd_in_sp(2,i)
        kk = odtothd_in_sp(3,i)
        !
        DO ish = 1, nord2
          !
          p1 = p( thdtood_in_sp( ii+ish, jj+ish, kk    )) &
              -p( thdtood_in_sp( ii+ish, jj-ish, kk    )) &
              -p( thdtood_in_sp( ii-ish, jj+ish, kk    )) &
              +p( thdtood_in_sp( ii-ish, jj-ish, kk    ))
          !
          ! the X stencil for cross derivatives
          !
          q(i) = q(i)+coeke(ish,1,2)*p1
          !
        END DO
        !
      END DO
      !$omp end parallel do
    end if
    !
    if (abs(coeke(1,1,3)).gt.eps12) then
      !$omp parallel do private(ii,jj,kk,p1,p2,p3)
      DO i = 1, n
        !
        ii = odtothd_in_sp(1,i)
        jj = odtothd_in_sp(2,i)
        kk = odtothd_in_sp(3,i)
        !
        DO ish = 1, nord2
          !
          p2 = p( thdtood_in_sp( ii+ish, jj,     kk+ish)) &
              -p( thdtood_in_sp( ii+ish, jj,     kk-ish)) &
              -p( thdtood_in_sp( ii-ish, jj,     kk+ish)) &
              +p( thdtood_in_sp( ii-ish, jj,     kk-ish))
          !
          ! the X stencil for cross derivatives
          !
          q(i) = q(i)+coeke(ish,1,3)*p2
          !
        END DO
        !
      END DO
      !$omp end parallel do
    END IF
    !
    IF (abs(coeke(1,2,3)).GT.eps12) THEN
      !$omp parallel do private(ii,jj,kk,p1,p2,p3)
      DO i = 1, n
        !
        ii = odtothd_in_sp(1,i)
        jj = odtothd_in_sp(2,i)
        kk = odtothd_in_sp(3,i)
        !
        DO ish = 1, nord2
          !
          p3 = p( thdtood_in_sp( ii,     jj+ish, kk+ish)) &
              -p( thdtood_in_sp( ii,     jj+ish, kk-ish)) &
              -p( thdtood_in_sp( ii,     jj-ish, kk+ish)) &
              +p( thdtood_in_sp( ii,     jj-ish, kk-ish))
          !
          ! the X stencil for cross derivatives
          !
          q(i) = q(i)+coeke(ish,2,3)*p3
          !
        END DO
        !
      END DO
      !$omp end parallel do
    END IF
    !
    RETURN
    !
END SUBROUTINE lapmvs
!----------------------------------------------------------------------c

!----------------------------------------------------------------------c
!                          S P A R S K I T                             c
!----------------------------------------------------------------------c
!         Basic Iterative Solvers with Reverse Communication           c
!----------------------------------------------------------------------c
!     They all have the following calling sequence:
!      subroutine solver(n, rhs, sol, ipar, fpar, w)
!      integer n, ipar(16)
!      real*8 rhs(n), sol(n), fpar(16), w(*)
!     Where
!     (1) 'n' is the size of the linear system,
!     (2) 'rhs' is the right-hand side of the linear system,
!     (3) 'sol' is the solution to the linear system,
!     (4) 'ipar' is an integer parameter array for the reverse
!     communication protocol,
!     (5) 'fpar' is an floating-point parameter array storing
!     information to and from the iterative solvers.
!     (6) 'w' is the work space (size is specified in ipar)
!
!     They are preconditioned iterative solvers with reverse
!     communication. The preconditioners can be applied from either
!     from left or right or both (specified by ipar(2), see below).
!
!     Author: Kesheng John Wu (kewu@mail.cs.umn.edu) 1993
!
!     NOTES:
!
!     (1) Work space required by each of the iterative solver
!     routines is as follows:
!       CG      == 5 * n
!       CGNR    == 5 * n
!       BCG     == 7 * n
!       DBCG    == 11 * n
!       BCGSTAB == 8 * n
!       TFQMR   == 11 * n
!       FOM     == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
!       GMRES   == (n+3)*(m+2) + (m+1)*m/2 (m = ipar(5), default m=15)
!       FGMRES  == n*(2m+1) + (m+1)*m/2 + 3*m + 2 (m = ipar(5),
!                  default m=15)
!       DQGMRES == n + lb * (2*n+4) (lb=ipar(5)+1, default lb = 16)
!
!     (2) ALL iterative solvers require a user-supplied DOT-product
!     routine named DISTDOT. The prototype of DISTDOT is
!
!     real*8 function distdot(n,x,ix,y,iy)
!     integer n, ix, iy
!     real*8 x(1+(n-1)*ix), y(1+(n-1)*iy)
!
!     This interface of DISTDOT is exactly the same as that of
!     DDOT (or SDOT if real == real*8) from BLAS-1. It should have
!     same functionality as DDOT on a single processor machine. On a
!     parallel/distributed environment, each processor can perform
!     DDOT on the data it has, then perform a summation on all the
!     partial results.
!
!     (3) To use this set of routines under SPMD/MIMD program paradigm,
!     several things are to be noted: (a) 'n' should be the number of
!     vector elements of 'rhs' that is present on the local processor.
!     (b) if RHS(i) is on processor j, it is expected that SOL(i)
!     will be on the same processor, i.e. the vectors are distributed
!     to each processor in the same way. (c) the preconditioning and
!     stopping criteria specifications have to be the same on all
!     processor involved, ipar and fpar have to be the same on each
!     processor. (d) DISTDOT should be replaced by a distributed
!     dot-product function.
!
!     ..................................................................
!     Reverse Communication Protocols
!
!     When a reverse-communication routine returns, it could be either
!     that the routine has terminated or it simply requires the caller
!     to perform one matrix-vector multiplication. The possible matrices
!     that involve in the matrix-vector multiplications are:
!     A       (the matrix of the linear system),
!     A^T     (A transposed),
!     Ml^{-1} (inverse of the left preconditioner),
!     Ml^{-T} (inverse of the left preconditioner transposed),
!     Mr^{-1} (inverse of the right preconditioner),
!     Mr^{-T} (inverse of the right preconditioner transposed).
!     For all the matrix vector multiplication, v = A u. The input and
!     output vectors are supposed to be part of the work space 'w', and
!     the starting positions of them are stored in ipar(8:9), see below.
!
!     The array 'ipar' is used to store the information about the solver.
!     Here is the list of what each element represents:
!
!     ipar(1) -- status of the call/return.
!     A call to the solver with ipar(1) == 0 will initialize the
!     iterative solver. On return from the iterative solver, ipar(1)
!     carries the status flag which indicates the condition of the
!     return. The status information is divided into two categories,
!     (1) a positive value indicates the solver requires a matrix-vector
!     multiplication,
!     (2) a non-positive value indicates termination of the solver.
!     Here is the current definition:
!       1 == request a matvec with A,
!       2 == request a matvec with A^T,
!       3 == request a left preconditioner solve (Ml^{-1}),
!       4 == request a left preconditioner transposed solve (Ml^{-T}),
!       5 == request a right preconditioner solve (Mr^{-1}),
!       6 == request a right preconditioner transposed solve (Mr^{-T}),
!      10 == request the caller to perform stopping test,
!       0 == normal termination of the solver, satisfied the stopping
!            criteria,
!      -1 == termination because iteration number is greater than the
!            preset limit,
!      -2 == return due to insufficient work space,
!      -3 == return due to anticipated break-down / divide by zero,
!            in the case where Arnoldi procedure is used, additional
!            error code can be found in ipar(12), where ipar(12) is
!            the error code of orthogonalization procedure MGSRO:
!               -1: zero input vector
!               -2: input vector contains abnormal numbers
!               -3: input vector is a linear combination of others
!               -4: trianguler system in GMRES/FOM/etc. has rank 0 (zero)
!      -4 == the values of fpar(1) and fpar(2) are both <= 0, the valid
!            ranges are 0 <= fpar(1) < 1, 0 <= fpar(2), and they can
!            not be zero at the same time
!      -9 == while trying to detect a break-down, an abnormal number is
!            detected.
!     -10 == return due to some non-numerical reasons, e.g. invalid
!            floating-point numbers etc.
!
!     ipar(2) -- status of the preconditioning:
!       0 == no preconditioning
!       1 == left preconditioning only
!       2 == right preconditioning only
!       3 == both left and right preconditioning
!
!     ipar(3) -- stopping criteria (details of this will be
!     discussed later).
!
!     ipar(4) -- number of elements in the array 'w'. if this is less
!     than the desired size, it will be over-written with the minimum
!     requirement. In which case the status flag ipar(1) = -2.
!
!     ipar(5) -- size of the Krylov subspace (used by GMRES and its
!     variants), e.g. GMRES(ipar(5)), FGMRES(ipar(5)),
!     DQGMRES(ipar(5)).
!
!     ipar(6) -- maximum number of matrix-vector multiplies, if not a
!     positive number the iterative solver will run till convergence
!     test is satisfied.
!
!     ipar(7) -- current number of matrix-vector multiplies. It is
!     incremented after each matrix-vector multiplication. If there
!     is preconditioning, the counter is incremented after the
!     preconditioning associated with each matrix-vector multiplication.
!
!     ipar(8) -- pointer to the input vector to the requested matrix-
!     vector multiplication.
!
!     ipar(9) -- pointer to the output vector of the requested matrix-
!     vector multiplication.
!
!     To perform v = A * u, it is assumed that u is w(ipar(8):ipar(8)+n-1)
!     and v is stored as w(ipar(9):ipar(9)+n-1).
!
!     ipar(10) -- the return address (used to determine where to go to
!     inside the iterative solvers after the caller has performed the
!     requested services).
!
!     ipar(11) -- the result of the external convergence test
!     On final return from the iterative solvers, this value
!     will be reflected by ipar(1) = 0 (details discussed later)
!
!     ipar(12) -- error code of MGSRO, it is
!                  1 if the input vector to MGSRO is linear combination
!                    of others,
!                  0 if MGSRO was successful,
!                 -1 if the input vector to MGSRO is zero,
!                 -2 if the input vector contains invalid number.
!
!     ipar(13) -- number of initializations. During each initilization
!                 residual norm is computed directly from M_l(b - A x).
!
!     ipar(14) to ipar(16) are NOT defined, they are NOT USED by
!     any iterative solver at this time.
!
!     Information about the error and tolerance are stored in the array
!     FPAR. So are some internal variables that need to be saved from
!     one iteration to the next one. Since the internal variables are
!     not the same for each routine, we only define the common ones.
!
!     The first two are input parameters:
!     fpar(1) -- the relative tolerance,
!     fpar(2) -- the absolute tolerance (details discussed later),
!
!     When the iterative solver terminates,
!     fpar(3) -- initial residual/error norm,
!     fpar(4) -- target residual/error norm,
!     fpar(5) -- current residual norm (if available),
!     fpar(6) -- current residual/error norm,
!     fpar(7) -- convergence rate,
!
!     fpar(8:10) are used by some of the iterative solvers to save some
!     internal information.
!
!     fpar(11) -- number of floating-point operations. The iterative
!     solvers will add the number of FLOPS they used to this variable,
!     but they do NOT initialize it, nor add the number of FLOPS due to
!     matrix-vector multiplications (since matvec is outside of the
!     iterative solvers). To insure the correct FLOPS count, the
!     caller should set fpar(11) = 0 before invoking the iterative
!     solvers and account for the number of FLOPS from matrix-vector
!     multiplications and preconditioners.
!
!     fpar(12:16) are not used in current implementation.
!
!     Whether the content of fpar(3), fpar(4) and fpar(6) are residual
!     norms or error norms depends on ipar(3). If the requested
!     convergence test is based on the residual norm, they will be
!     residual norms. If the caller want to test convergence based the
!     error norms (estimated by the norm of the modifications applied
!     to the approximate solution), they will be error norms.
!     Convergence rate is defined by (Fortran 77 statement)
!     fpar(7) = log10(fpar(3) / fpar(6)) / (ipar(7)-ipar(13))
!     If fpar(7) = 0.5, it means that approximately every 2 (= 1/0.5)
!     steps the residual/error norm decrease by a factor of 10.
!
!     ..................................................................
!     Stopping criteria,
!
!     An iterative solver may be terminated due to (1) satisfying
!     convergence test; (2) exceeding iteration limit; (3) insufficient
!     work space; (4) break-down. Checking of the work space is
!     only done in the initialization stage, i.e. when it is called with
!     ipar(1) == 0. A complete convergence test is done after each
!     update of the solutions. Other conditions are monitored
!     continuously.
!
!     With regard to the number of iteration, when ipar(6) is positive,
!     the current iteration number will be checked against it. If
!     current iteration number is greater the ipar(6) than the solver
!     will return with status -1. If ipar(6) is not positive, the
!     iteration will continue until convergence test is satisfied.
!
!     Two things may be used in the convergence tests, one is the
!     residual 2-norm, the other one is 2-norm of the change in the
!     approximate solution. The residual and the change in approximate
!     solution are from the preconditioned system (if preconditioning
!     is applied). The DQGMRES and TFQMR use two estimates for the
!     residual norms. The estimates are not accurate, but they are
!     acceptable in most of the cases. Generally speaking, the error
!     of the TFQMR's estimate is less accurate.
!
!     The convergence test type is indicated by ipar(3). There are four
!     type convergence tests: (1) tests based on the residual norm;
!     (2) tests based on change in approximate solution; (3) caller
!     does not care, the solver choose one from above two on its own;
!     (4) caller will perform the test, the solver should simply continue.
!     Here is the complete definition:
!      -2 == || dx(i) || <= rtol * || rhs || + atol
!      -1 == || dx(i) || <= rtol * || dx(1) || + atol
!       0 == solver will choose test 1 (next)
!       1 == || residual || <= rtol * || initial residual || + atol
!       2 == || residual || <= rtol * || rhs || + atol
!     999 == caller will perform the test
!     where dx(i) denote the change in the solution at the ith update.
!     ||.|| denotes 2-norm. rtol = fpar(1) and atol = fpar(2).
!
!     If the caller is to perform the convergence test, the outcome
!     should be stored in ipar(11).
!     ipar(11) = 0 -- failed the convergence test, iterative solver
!     should continue
!     ipar(11) = 1 -- satisfied convergence test, iterative solver
!     should perform the clean up job and stop.
!
!     Upon return with ipar(1) = 10,
!     ipar(8)  points to the starting position of the change in
!              solution Sx, where the actual solution of the step is
!              x_j = x_0 + M_r^{-1} Sx.
!              Exception: ipar(8) < 0, Sx = 0. It is mostly used by
!              GMRES and variants to indicate (1) Sx was not necessary,
!              (2) intermediate result of Sx is not computed.
!     ipar(9)  points to the starting position of a work vector that
!              can be used by the caller.
!
!     NOTE: the caller should allow the iterative solver to perform
!     clean up job after the external convergence test is satisfied,
!     since some of the iterative solvers do not directly
!     update the 'sol' array. A typical clean-up stage includes
!     performing the final update of the approximate solution and
!     computing the convergence information (e.g. values of fpar(3:7)).
!
!     NOTE: fpar(4) and fpar(6) are not set by the accelerators (the
!     routines implemented here) if ipar(3) = 999.
!
!     ..................................................................
!     Usage:
!
!     To start solving a linear system, the user needs to specify
!     first 6 elements of the ipar, and first 2 elements of fpar.
!     The user may optionally set fpar(11) = 0 if one wants to count
!     the number of floating-point operations. (Note: the iterative
!     solvers will only add the floating-point operations inside
!     themselves, the caller will have to add the FLOPS from the
!     matrix-vector multiplication routines and the preconditioning
!     routines in order to account for all the arithmetic operations.)
!
!     Here is an example:
!     ipar(1) = 0       ! always 0 to start an iterative solver
!     ipar(2) = 2       ! right preconditioning
!     ipar(3) = 1       ! use convergence test scheme 1
!     ipar(4) = 10000   ! the 'w' has 10,000 elements
!     ipar(5) = 10      ! use *GMRES(10) (e.g. FGMRES(10))
!     ipar(6) = 100     ! use at most 100 matvec's
!     fpar(1) = 1.0D-6  ! relative tolerance 1.0D-6
!     fpar(2) = 1.0D-10 ! absolute tolerance 1.0D-10
!     fpar(11) = 0.0    ! clearing the FLOPS counter
!
!     After the above specifications, one can start to call an iterative
!     solver, say BCG. Here is a piece of pseudo-code showing how it can
!     be done,
!
! 10   call bcg(n,rhs,sol,ipar,fpar,w)
!      if (ipar(1).eq.1) then
!         call amux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
!         goto 10
!      else if (ipar(1).eq.2) then
!         call atmux(n,w(ipar(8)),w(ipar(9)),a,ja,ia)
!         goto 10
!      else if (ipar(1).eq.3) then
!         left preconditioner solver
!         goto 10
!      else if (ipar(1).eq.4) then
!         left preconditioner transposed solve
!         goto 10
!      else if (ipar(1).eq.5) then
!         right preconditioner solve
!         goto 10
!      else if (ipar(1).eq.6) then
!         right preconditioner transposed solve
!         goto 10
!      else if (ipar(1).eq.10) then
!         call my own stopping test routine
!         goto 10
!      else if (ipar(1).gt.0) then
!         ipar(1) is an unspecified code
!      else
!         the iterative solver terminated with code = ipar(1)
!      endif
!
!     This segment of pseudo-code assumes the matrix is in CSR format,
!     AMUX and ATMUX are two routines from the SPARSKIT MATVEC module.
!     They perform matrix-vector multiplications for CSR matrices,
!     where w(ipar(8)) is the first element of the input vectors to the
!     two routines, and w(ipar(9)) is the first element of the output
!     vectors from them. For simplicity, we did not show the name of
!     the routine that performs the preconditioning operations or the
!     convergence tests.
!-----------------------------------------------------------------------
subroutine cg(n, rhs, sol, ipar, fpar, w)
    implicit none
    integer n, ipar(16)
    real*8 rhs(n), sol(n), fpar(16), w(n,*)
    !-----------------------------------------------------------------------
    !     This is a implementation of the Conjugate Gradient (CG) method
    !     for solving linear system.
    !
    !     NOTE: This is not the PCG algorithm. It is a regular CG algorithm.
    !     To be consistent with the other solvers, the preconditioners are
    !     applied by performing Ml^{-1} A Mr^{-1} P in place of A P in the
    !     CG algorithm. The PCG uses its preconditioners very differently.
    !
    !     fpar(7) is used here internally to store <r, r>.
    !     w(:,1) -- residual vector
    !     w(:,2) -- P, the conjugate direction
    !     w(:,3) -- A P, matrix multiply the conjugate direction
    !     w(:,4) -- temporary storage for results of preconditioning
    !     w(:,5) -- change in the solution (sol) is stored here until
    !               termination of this solver
    !-----------------------------------------------------------------------
    !     external functions used
    real*8 distdot
    logical stopbis, brkdn
    external distdot, stopbis, brkdn, bisinit
    !
    !     local variables
    !
    integer i
    real*8 alpha
    logical lp,rp
    save
    !
    !     check the status of the call
    !
    if (ipar(1).le.0) ipar(10) = 0
    goto (10, 20, 40, 50, 60, 70, 80), ipar(10)
    !
    !     initialization
    !
    call bisinit(ipar,fpar,5*n,1,lp,rp,w)
    if (ipar(1).lt.0) return
    !
    !     request for matrix vector multiplication A*x in the initialization
    !
    ipar(1) = 1
    ipar(8) = n+1
    ipar(9) = ipar(8) + n
    ipar(10) = 1
    !$omp parallel do 
    do i = 1, n
      w(i,2) = sol(i)
    enddo
    !$omp end parallel do 
    return
    10   ipar(7) = ipar(7) + 1
    ipar(13) = 1
    !$omp parallel do 
    do i = 1, n
      w(i,2) = rhs(i) - w(i,3)
    enddo
    !$omp end parallel do 
    fpar(11) = fpar(11) + n
    !
    !     if left preconditioned
    !
    if (lp) then
      ipar(1) = 3
      ipar(9) = 1
      ipar(10) = 2
      return
    endif
    !
    20   if (lp) then
      !$omp parallel do 
      do i = 1, n
        w(i,2) = w(i,1)
      enddo
      !$omp end parallel do 
    else
      !$omp parallel do 
      do i = 1, n
        w(i,1) = w(i,2)
      enddo
      !$omp end parallel do 
    endif
    !
    fpar(7) = distdot(n,w,1,w,1)
    fpar(11) = fpar(11) + 2 * n
    fpar(3) = sqrt(fpar(7))
    fpar(5) = fpar(3)
    if (abs(ipar(3)).eq.2) then
      fpar(4) = fpar(1) * sqrt(distdot(n,rhs,1,rhs,1)) + fpar(2)
      fpar(11) = fpar(11) + 2 * n
    else if (ipar(3).ne.999) then
      fpar(4) = fpar(1) * fpar(3) + fpar(2)
    endif
    !
    !     before iteration can continue, we need to compute A * p, which
    !     includes the preconditioning operations
    !
    30   if (rp) then
      ipar(1) = 5
      ipar(8) = n + 1
      if (lp) then
        ipar(9) = ipar(8) + n
      else
        ipar(9) = 3*n + 1
      endif
      ipar(10) = 3
      return
    endif
    !
    40   ipar(1) = 1
    if (rp) then
      ipar(8) = ipar(9)
    else
      ipar(8) = n + 1
    endif
    if (lp) then
      ipar(9) = 3*n+1
    else
      ipar(9) = n+n+1
    endif
    ipar(10) = 4
    return
    !
    50   if (lp) then
      ipar(1) = 3
      ipar(8) = ipar(9)
      ipar(9) = n+n+1
      ipar(10) = 5
      return
    endif
    !
    !     continuing with the iterations
    !
    60   ipar(7) = ipar(7) + 1
    alpha = distdot(n,w(1,2),1,w(1,3),1)
    fpar(11) = fpar(11) + 2*n
    if (brkdn(alpha,ipar)) goto 900
    alpha = fpar(7) / alpha
    !$omp parallel do 
    do i = 1, n
      w(i,5) = w(i,5) + alpha * w(i,2)
      w(i,1) = w(i,1) - alpha * w(i,3)
    enddo
    !$omp end parallel do 
    fpar(11) = fpar(11) + 4*n
    !
    !     are we ready to terminate ?
    !
    if (ipar(3).eq.999) then
      ipar(1) = 10
      ipar(8) = 4*n + 1
      ipar(9) = 3*n + 1
      ipar(10) = 6
      return
    endif
    70   if (ipar(3).eq.999) then
      if (ipar(11).eq.1) goto 900
    else if (stopbis(n,ipar,1,fpar,w,w(1,2),alpha)) then
      goto 900
    endif
    !
    !     continue the iterations
    !
    alpha = fpar(5)*fpar(5) / fpar(7)
    fpar(7) = fpar(5)*fpar(5)
    !$omp parallel do 
    do i = 1, n
      w(i,2) = w(i,1) + alpha * w(i,2)
    enddo
    !$omp end parallel do 
    fpar(11) = fpar(11) + 2*n
    goto 30
    !
    !     clean up -- necessary to accommodate the right-preconditioning
    !
    900  if (rp) then
      if (ipar(1).lt.0) ipar(12) = ipar(1)
      ipar(1) = 5
      ipar(8) = 4*n + 1
      ipar(9) = ipar(8) - n
      ipar(10) = 7
      return
    endif
    80   if (rp) then
      call tidycg(n,ipar,fpar,sol,w(1,4))
    else
      call tidycg(n,ipar,fpar,sol,w(1,5))
    endif
    !
    return
end subroutine cg
!-----end-of-cg

!-----------------------------------------------------------------------
logical function stopbis(n,ipar,mvpi,fpar,r,delx,sx)
    implicit none
    integer n,mvpi,ipar(16)
    real*8 fpar(16), r(n), delx(n), sx, distdot
    external distdot
    !-----------------------------------------------------------------------
    !     function for determining the stopping criteria. return value of
    !     true if the stopbis criteria is satisfied.
    !-----------------------------------------------------------------------
    if (ipar(11) .eq. 1) then
      stopbis = .true.
    else
      stopbis = .false.
    endif
    if (ipar(6).gt.0 .and. ipar(7).ge.ipar(6)) then
      ipar(1) = -1
      stopbis = .true.
    endif
    if (stopbis) return
    !
    !     computes errors
    !
    fpar(5) = sqrt(distdot(n,r,1,r,1))
    fpar(11) = fpar(11) + 2 * n
    if (ipar(3).lt.0) then
      !
      !     compute the change in the solution vector
      !
      fpar(6) = sx * sqrt(distdot(n,delx,1,delx,1))
      fpar(11) = fpar(11) + 2 * n
      if (ipar(7).lt.mvpi+mvpi+1) then
        !
        !     if this is the end of the first iteration, set fpar(3:4)
        !
        fpar(3) = fpar(6)
        if (ipar(3).eq.-1) then
          fpar(4) = fpar(1) * fpar(3) + fpar(2)
        endif
        !
      endif
      !
    else
      !
      fpar(6) = fpar(5)
      !
    endif
    !
    !     .. the test is struct this way so that when the value in fpar(6)
    !       is not a valid number, STOPBIS is set to .true.
    !
    if (fpar(6).gt.fpar(4)) then
      stopbis = .false.
      ipar(11) = 0
    else
      stopbis = .true.
      ipar(11) = 1
    endif
    !
    return
end function stopbis
!-----end-of-stopbis

!-----------------------------------------------------------------------
subroutine tidycg(n,ipar,fpar,sol,delx)
    implicit none
    integer i,n,ipar(16)
    real*8 fpar(16),sol(n),delx(n)
    !-----------------------------------------------------------------------
    !     Some common operations required before terminating the CG routines
    !-----------------------------------------------------------------------
    real*8 zero
    parameter(zero=0.0D0)
    !
    if (ipar(12).ne.0) then
      ipar(1) = -3
    else if (ipar(1).gt.0) then
      if ((ipar(3).eq.999 .and. ipar(11).eq.1) .or. fpar(6).le.fpar(4)) then
        ipar(1) = 0
      else if (ipar(7).ge.ipar(6) .and. ipar(6).gt.0) then
        ipar(1) = -1
      else
        ipar(1) = -10
      endif
    endif
    if (fpar(3).gt.zero .and. fpar(6).gt.zero .and. ipar(7).gt.ipar(13)) then
      fpar(7) = log10(fpar(3) / fpar(6)) / dble(ipar(7)-ipar(13))
    else
      fpar(7) = zero
    endif
    !$omp parallel do 
    do i = 1, n
      sol(i) = sol(i) + delx(i)
    enddo
    !$omp end parallel do 
    return
end subroutine tidycg
!-----end-of-tidycg

!-----------------------------------------------------------------------
logical function brkdn(alpha, ipar)
    implicit none
    integer ipar(16)
    real*8 alpha, beta, zero, one
    parameter (zero=0.0D0, one=1.0D0)
    !-----------------------------------------------------------------------
    !     test whether alpha is zero or an abnormal number, if yes,
    !     this routine will return .true.
    !
    !     If alpha == 0, ipar(1) = -3,
    !     if alpha is an abnormal number, ipar(1) = -9.
    !-----------------------------------------------------------------------
    brkdn = .false.
    if (alpha.gt.zero) then
      beta = one / alpha
      if (.not. beta.gt.zero) then
        brkdn = .true.
        ipar(1) = -9
      endif
    else if (alpha.lt.zero) then
      beta = one / alpha
      if (.not. beta.lt.zero) then
        brkdn = .true.
        ipar(1) = -9
      endif
    else if (alpha.eq.zero) then
      brkdn = .true.
      ipar(1) = -3
    else
      brkdn = .true.
      ipar(1) = -9
    endif
    return
end function brkdn
!-----end-of-brkdn

!-----------------------------------------------------------------------
subroutine bisinit(ipar,fpar,wksize,dsc,lp,rp,wk)
    implicit none
    integer i,ipar(16),wksize,dsc
    logical lp,rp
    real*8  fpar(16),wk(*)
    !-----------------------------------------------------------------------
    !     some common initializations for the iterative solvers
    !-----------------------------------------------------------------------
    real*8 zero, one
    parameter(zero=0.0D0, one=1.0D0)
    !
    !     ipar(1) = -2 inidcate that there are not enough space in the work
    !     array
    !
    if (ipar(4).lt.wksize) then
      ipar(1) = -2
      ipar(4) = wksize
      return
    endif
    !
    if (ipar(2).gt.2) then
      lp = .true.
      rp = .true.
    else if (ipar(2).eq.2) then
      lp = .false.
      rp = .true.
    else if (ipar(2).eq.1) then
      lp = .true.
      rp = .false.
    else
      lp = .false.
      rp = .false.
    endif
    if (ipar(3).eq.0) ipar(3) = dsc
    !     .. clear the ipar elements used
    ipar(7) = 0
    ipar(8) = 0
    ipar(9) = 0
    ipar(10) = 0
    ipar(11) = 0
    ipar(12) = 0
    ipar(13) = 0
    !
    !     fpar(1) must be between (0, 1), fpar(2) must be positive,
    !     fpar(1) and fpar(2) can NOT both be zero
    !     Normally return ipar(1) = -4 to indicate any of above error
    !
    if (fpar(1).lt.zero .or. fpar(1).ge.one .or. fpar(2).lt.zero .or.  &
        &     (fpar(1).eq.zero .and. fpar(2).eq.zero)) then
    if (ipar(1).eq.0) then
      ipar(1) = -4
      return
    else
      fpar(1) = 1.0D-6
      fpar(2) = 1.0D-16
    endif
  endif
  !     .. clear the fpar elements
  do i = 3, 10
    fpar(i) = zero
  enddo
  if (fpar(11).lt.zero) fpar(11) = zero
  !     .. clear the used portion of the work array to zero
  do i = 1, wksize
    wk(i) = zero
  enddo
  !
  return
end subroutine bisinit
!-----end-of-bisinit
!-----------------------------------------------------------------------
