!
! Copyright (C) 2016 - present Andrea Dal Corso
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE linear_solvers

  USE io_global, ONLY : stdout
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE

  PUBLIC  ccg_many_vectors, cg_many_vectors, linsolvx, linsolvx_sym, &
          linsolvms, linsolvsvd, min_sqr_solve


CONTAINS

!----------------------------------------------------------------------
SUBROUTINE cg_many_vectors (apply_a, precond, scal_prod, d0psi, dpsi, &
                h_diag, ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the linear system:
  !
  !                  A  * dpsi = d0psi                      (1)
  !
  !     where A is an hermitean positive definite matrix, dpsi and d0psi are 
  !     complex vectors. It uses the conjugate gradient method.
  !     If you have an general complex operator you have to apply the 
  !     slower ccg_many_vector routine.
  !
  !     on input:
  !                 apply_a  EXTERNAL  name of a subroutine:
  !                          apply_a(ndmx,ndim,psi,psip,)
  !                          Calculates  H*psi products.
  !                          Vectors psi and psip should be dimensined
  !                          (ndmx,nvec). nvec=1 is used!
  !
  !                 cg_psi   EXTERNAL  name of a subroutine:
  !                          g_psi(ndmx,ndim,notcnv,psi,e)
  !                          which calculates (h-e)^-1 * psi, with
  !                          some approximation, e.g. (diag(h)-e)
  !
  !                 scal_prod EXTERNAL name of a subroutine
  !                          scal_prod(ndmx,ndim,psi1,psi2)
  !                          which calculate the scalar product
  !                          psi1 ^+ psi2
  !              
  !                 dpsi     contains an estimate of the solution
  !                          vector.
  !
  !                 d0psi    contains the right hand side vector
  !                          of the system.
  !
  !                 ndmx     integer row dimension of dpsi, ecc.
  !
  !                 ndim     integer actual row dimension of dpsi
  !
  !                 ethr     real     convergence threshold. solution
  !                          improvement is stopped when the error in
  !                          eq (1), defined as l.h.s. - r.h.s., becomes
  !                          less than ethr in norm.
  !
  !     on output:  dpsi     contains the refined estimate of the
  !                          solution vector.
  !
  !                 d0psi    is corrupted on exit
  !
  !   written 29/01/2016  by A. Dal Corso modifying the old routine
  !   contained in cgsolve_all of Quantum ESPRESSO.
  !
  USE kinds,          ONLY : DP
  USE io_global,       ONLY : stdout

  IMPLICIT NONE
  !
  !   first the I/O variables
  !
  INTEGER :: ndmx, & ! input: the maximum dimension of the vectors
             ndim, & ! input: the actual dimension of the vectors
             kter, & ! output: counter on iterations
             nbnd, & ! input: the number of bands
             ik      ! input: the k point

  REAL(DP) ::         &
             anorm,   & ! output: the norm of the error in the solution
             h_diag(ndmx,nbnd), & ! input: an estimate of ( H - \epsilon )
             ethr       ! input: the required precision

  COMPLEX(DP) ::                &
             dpsi (ndmx, nbnd), & ! output: the solution of the linear syst
             d0psi (ndmx, nbnd)   ! input: the known term

  LOGICAL :: conv_root ! output: if true the root is converged
  EXTERNAL apply_a     ! input: the routine computing A
  EXTERNAL precond     ! input: the routine applying the preconditioning
  EXTERNAL scal_prod   ! input: the routine computing the scalar product
  !
  !  here the local variables
  !
  INTEGER, PARAMETER :: maxter = 200
  ! the maximum number of iterations
  INTEGER :: iter, ibnd, lbnd
  ! counters on iteration, bands
  INTEGER , ALLOCATABLE :: conv (:)
  ! if 1 the root is converged

  COMPLEX(DP), ALLOCATABLE :: g (:,:), t (:,:), h (:,:), hold (:,:)
  !  the gradient of psi
  !  the preconditioned gradient
  !  the delta gradient
  !  the conjugate gradient
  !  work space
  REAL(DP) ::  dcgamma, dclambda
  !  the ratio between rho
  !  step length
  REAL(DP), ALLOCATABLE :: rho (:), rhoold (:), a(:), c(:)
  ! the residue
  ! auxiliary for h_diag
  !
  COMPLEX(DP) :: scal_prod
  !
  INTEGER, ALLOCATABLE :: indb(:) 
  !
  REAL(DP) :: kter_eff
  ! account the number of iterations with b
  !
  CALL start_clock ('cg_many_vectors')

  ALLOCATE ( g(ndmx,nbnd), t(ndmx,nbnd), h(ndmx,nbnd),hold(ndmx ,nbnd) )
  ALLOCATE (a(nbnd), c(nbnd))
  ALLOCATE (conv (nbnd))
  ALLOCATE (rho(nbnd),rhoold(nbnd))
  ALLOCATE (indb(nbnd))
  !
  !      WRITE( stdout,*) g,t,h,hold
  !
  kter_eff = 0.d0
  conv= 0
  g=(0.d0,0.d0)
  t=(0.d0,0.d0)
  h=(0.d0,0.d0)
  hold=(0.d0,0.d0)
  DO ibnd=1,nbnd
     indb(ibnd)=ibnd
  END DO
  DO iter = 1, maxter
     !
     !    compute the gradient. can reuse information from previous step
     !
     IF (iter == 1) THEN
        CALL apply_a (ndmx, ndim, dpsi, g, ik, nbnd, indb, 1)
        DO ibnd = 1, nbnd
           CALL zaxpy (ndmx, (-1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd), 1)
        ENDDO
     ENDIF
     !
     !    compute preconditioned residual vectors and convergence check
     !
     lbnd = 0
     DO ibnd = 1, nbnd
        IF (conv (ibnd) == 0) THEN
           lbnd = lbnd+1
           CALL zcopy (ndmx, g (1, ibnd), 1, h (1, ibnd), 1)
           CALL precond(ndmx, ndim, 1, h(1,ibnd), h_diag(1,ibnd) )
           rho(lbnd) = scal_prod (ndmx, ndim, h(1,ibnd), g(1,ibnd) )
        ENDIF
     ENDDO
     kter_eff = kter_eff + DBLE (lbnd) / DBLE (nbnd)
     DO ibnd = nbnd, 1, -1
        IF (conv(ibnd)==0) THEN
           rho(ibnd)=rho(lbnd)
           lbnd = lbnd -1
           anorm = SQRT (rho (ibnd) )
!           WRITE(stdout,*) iter, ibnd, anorm
           IF (anorm < ethr) conv (ibnd) = 1
        ENDIF
     ENDDO
!
     conv_root = .true.
     DO ibnd = 1, nbnd
        conv_root = conv_root.AND. (conv (ibnd) .eq.1)
     ENDDO
     IF (conv_root) GOTO 100
     !
     !        compute the step directions h and hs. 
     !
     lbnd = 0
     DO ibnd = 1, nbnd
        IF (conv (ibnd) == 0) THEN
!
!          change sign to h and hs
!
           CALL dscal (2 * ndmx, - 1.d0, h (1, ibnd), 1)
           IF (iter /= 1) THEN
              dcgamma = rho (ibnd) / rhoold (ibnd)
              CALL zaxpy (ndmx, dcgamma, hold (1, ibnd), 1, h (1, ibnd), 1)
           ENDIF

!
! here hold is used as auxiliary vector in order to efficiently compute t = A*h
! it is later set to the current (becoming old) value of h
!
           lbnd = lbnd+1
           CALL zcopy (ndmx, h (1, ibnd), 1, hold (1, lbnd), 1)
           indb (lbnd) = ibnd
        ENDIF
     ENDDO
     !
     !        compute t = A*h  
     !
     CALL apply_a (ndmx, ndim, hold, t, ik, lbnd, indb, 1)
     !
     !        compute the coefficients a and c for the line minimization
     !        compute step length lambda
     lbnd=0
     DO ibnd = 1, nbnd
        IF (conv (ibnd) == 0) THEN
           lbnd=lbnd+1
           a(lbnd) = scal_prod (ndmx, ndim, h(1,ibnd), g(1,ibnd))
           c(lbnd) = scal_prod (ndmx, ndim, h(1,ibnd), t(1,lbnd))
        END IF
     END DO
     lbnd=0
     DO ibnd = 1, nbnd
        IF (conv (ibnd) == 0) THEN
           lbnd=lbnd+1
           dclambda = CMPLX( - a(lbnd) / c(lbnd), 0.d0,kind=DP)
           !
           !    move to new position
           !
           CALL zaxpy (ndmx, dclambda, h(1,ibnd), 1, dpsi(1,ibnd), 1)
           !
           !    update to get the gradient
           !
           !g=g+lam
           CALL zaxpy (ndmx, dclambda, t(1,lbnd), 1, g(1,ibnd), 1)
           !
           !    save current (now old) h and rho for later use
           !
           CALL zcopy (ndmx, h(1,ibnd), 1, hold(1,ibnd), 1)
           rhoold (ibnd) = rho (ibnd)
        ENDIF
     ENDDO
  ENDDO
100 CONTINUE
  kter = kter_eff
  DEALLOCATE (indb)
  DEALLOCATE (rho, rhoold)
  DEALLOCATE (conv)
  DEALLOCATE (a,c)
  DEALLOCATE (g, t, h, hold)

  CALL stop_clock ('cg_many_vectors')
  RETURN
END SUBROUTINE cg_many_vectors

!----------------------------------------------------------------------
SUBROUTINE ccg_many_vectors (apply_a, precond, scal_prod, d0psi, dpsi, &
                h_diag, ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the linear system:
  !
  !                  A  * dpsi = d0psi                      (1)
  !
  !     where A is a complex general matrix, dpsi and d0psi are 
  !     complex vectors. It uses the biconjugate gradient method.
  !     Note that this routine is for general matrices and can deal with 
  !     many vectors d0psi.
  !     If you have an hermitean positive definite operator you can
  !     apply the much faster cg_many_vector routine.
  !     If you have only one vector d0psi you can apply the 
  !     cg_one_vector ccg_one_vector routines.
  !
  !     on input:
  !                 apply_a  EXTERNAL  name of a subroutine:
  !                          apply_a(ndmx,ndim,psi,psip,)
  !                          Calculates  H*psi products.
  !                          Vectors psi and psip should be dimensined
  !                          (ndmx,nvec). nvec=1 is used!
  !
  !                 cg_psi   EXTERNAL  name of a subroutine:
  !                          g_psi(ndmx,ndim,notcnv,psi,e)
  !                          which calculates (h-e)^-1 * psi, with
  !                          some approximation, e.g. (diag(h)-e)
  !
  !                 scal_prod EXTERNAL name of a subroutine
  !                          scal_prod(ndmx,ndim,psi1,psi2)
  !                          which calculate the scalar product
  !                          psi1 ^+ psi2
  !              
  !                 dpsi     contains an estimate of the solution
  !                          vector.
  !
  !                 d0psi    contains the right hand side vector
  !                          of the system.
  !
  !                 ndmx     integer row dimension of dpsi, ecc.
  !
  !                 ndim     integer actual row dimension of dpsi
  !
  !                 ethr     real     convergence threshold. solution
  !                          improvement is stopped when the error in
  !                          eq (1), defined as l.h.s. - r.h.s., becomes
  !                          less than ethr in norm.
  !
  !     on output:  dpsi     contains the refined estimate of the
  !                          solution vector.
  !
  !                 d0psi    is corrupted on exit
  !
  !   written 29/01/2016  by A. Dal Corso 
  !
  USE kinds,          ONLY : DP
  USE io_global,       ONLY : stdout

  IMPLICIT NONE
  !
  !   first the I/O variables
  !
  INTEGER :: ndmx, & ! input: the maximum dimension of the vectors
             ndim, & ! input: the actual dimension of the vectors
             kter, & ! output: counter on iterations
             nbnd, & ! input: the number of bands
             ik      ! input: the k point

  REAL(DP) :: &
             anorm,   & ! output: the norm of the error in the solution
             ethr       ! input: the required precision

  COMPLEX(DP) :: &
             h_diag(ndmx,nbnd), & ! input: an estimate of ( H - \epsilon- w )
                                  ! w can be complex
             dpsi (ndmx, nbnd), & ! output: the solution of the linear syst
             d0psi (ndmx, nbnd)   ! input: the known term

  LOGICAL :: conv_root ! output: if true the root is converged
  EXTERNAL apply_a     ! input: the routine computing A
  EXTERNAL precond     ! input: the routine applying the preconditioning
  EXTERNAL scal_prod   ! input: the routine computing the scalar product
  !
  !  here the local variables
  !
  INTEGER, PARAMETER :: maxter = 2000
  ! the maximum number of iterations
  INTEGER :: iter, ibnd, lbnd
  ! counters on iteration, bands
  INTEGER , ALLOCATABLE :: conv (:)
  ! if 1 the root is converged

  COMPLEX(DP), ALLOCATABLE :: g (:,:), t (:,:), h (:,:), hold (:,:)
  COMPLEX(DP), ALLOCATABLE :: gs (:,:), ts (:,:), hs (:,:), hsold (:,:)
  !  the gradient of psi
  !  the preconditioned gradient
  !  the delta gradient
  !  the conjugate gradient
  !  work space
  COMPLEX(DP) ::  dcgamma, dcgamma1, dclambda, dclambda1
  !  the ratio between rho
  !  step length
  COMPLEX(DP), ALLOCATABLE :: rho (:), rhoold (:), a(:), c(:)
  ! the residue
  ! auxiliary for h_diag
  !
  COMPLEX(DP) :: scal_prod
  !
  INTEGER, ALLOCATABLE :: indb(:) 
  !
  REAL(DP) :: kter_eff
  ! account the number of iterations with b
  !
  CALL start_clock ('ccg_many_vectors')

  ALLOCATE ( g(ndmx,nbnd), t(ndmx,nbnd), h(ndmx,nbnd),hold(ndmx ,nbnd) )
  ALLOCATE ( gs(ndmx,nbnd), ts(ndmx,nbnd), hs(ndmx,nbnd), hsold(ndmx,nbnd) )
  ALLOCATE (a(nbnd), c(nbnd))
  ALLOCATE (conv (nbnd))
  ALLOCATE (rho(nbnd),rhoold(nbnd))
  ALLOCATE (indb(nbnd))
  !      WRITE( stdout,*) g,t,h,hold

  kter_eff = 0.d0
  conv= 0
  g=(0.d0,0.d0)
  t=(0.d0,0.d0)
  h=(0.d0,0.d0)
  hold=(0.d0,0.d0)
  gs=(0.d0,0.d0)
  ts=(0.d0,0.d0)
  hs=(0.d0,0.d0)
  hsold=(0.d0,0.d0)
  DO ibnd=1,nbnd
     indb(ibnd)=ibnd
  END DO
  DO iter = 1, maxter
     !
     !    compute the gradient. can reuse information from previous step
     !
     IF (iter == 1) THEN
        CALL apply_a (ndmx, ndim, dpsi, g, ik, nbnd, indb, 1)
        DO ibnd = 1, nbnd
           CALL zaxpy (ndmx, (-1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd), 1)
        ENDDO
        gs(:,:) = CONJG(g(:,:))
     ENDIF
     !
     !    compute preconditioned residual vectors and convergence check
     !
     lbnd = 0
     DO ibnd = 1, nbnd
        IF (conv (ibnd) == 0) THEN
           lbnd = lbnd+1
           CALL zcopy (ndmx, g (1, ibnd), 1, h (1, ibnd), 1)
           CALL zcopy (ndmx, gs (1, ibnd), 1, hs (1, ibnd), 1)
           CALL precond(ndmx, ndim, 1, h(1,ibnd), h_diag(1,ibnd), 1 )
           CALL precond(ndmx, ndim, 1, hs(1,ibnd), h_diag(1,ibnd), -1 )
           rho(lbnd) = scal_prod (ndmx, ndim, hs(1,ibnd), g(1,ibnd))
        ENDIF
     ENDDO
     kter_eff = kter_eff + DBLE (lbnd) / DBLE (nbnd)
     DO ibnd = nbnd, 1, -1
        IF (conv(ibnd)==0) THEN
           rho(ibnd)=rho(lbnd)
           lbnd = lbnd -1
           anorm = SQRT (ABS (rho (ibnd)) )
!           IF (MOD(iter,20)==0) WRITE(stdout,*) iter, ibnd, anorm
           IF (anorm < ethr) conv (ibnd) = 1
        ENDIF
     ENDDO
!
     conv_root = .true.
     DO ibnd = 1, nbnd
        conv_root = conv_root.AND. (conv (ibnd) .eq.1)
     ENDDO
     IF (conv_root) GOTO 100
     !
     !        compute the step directions h and hs. 
     !
     lbnd = 0
     DO ibnd = 1, nbnd
        IF (conv (ibnd) == 0) THEN
!
!          change sign to h and hs
!
           CALL dscal (2 * ndmx, - 1.d0, h (1, ibnd), 1)
           CALL dscal (2 * ndmx, - 1.d0, hs (1, ibnd), 1)
           IF (iter /= 1) THEN
              dcgamma = rho (ibnd) / rhoold (ibnd)
              dcgamma1 = CONJG(dcgamma)
              CALL zaxpy (ndmx, dcgamma, hold (1, ibnd), 1, h (1, ibnd), 1)
              CALL zaxpy (ndmx, dcgamma1, hsold (1, ibnd), 1, hs (1, ibnd), 1)
           ENDIF

!
! here hold is used as auxiliary vector in order to efficiently compute t = A*h
! it is later set to the current (becoming old) value of h
!
           lbnd = lbnd+1
           CALL zcopy (ndmx, h (1, ibnd), 1, hold (1, lbnd), 1)
           CALL zcopy (ndmx, hs (1, ibnd), 1, hsold (1, lbnd), 1)
           indb (lbnd) = ibnd
        ENDIF
     ENDDO
     !
     !        compute t = A*h  ts= A^+ * h 
     !
     CALL apply_a (ndmx, ndim, hold, t, ik, lbnd, indb, 1)
     CALL apply_a (ndmx, ndim, hsold, ts, ik, lbnd, indb, -1)
     !
     !        compute the coefficients a and c for the line minimization
     !        compute step length lambda
     lbnd=0
     DO ibnd = 1, nbnd
        IF (conv (ibnd) == 0) THEN
           lbnd=lbnd+1
           a(lbnd) = scal_prod (ndmx, ndim, hs(1,ibnd), g(1,ibnd))
           c(lbnd) = scal_prod (ndmx, ndim, hs(1,ibnd), t(1,lbnd))
        END IF
     END DO
     lbnd=0
     DO ibnd = 1, nbnd
        IF (conv (ibnd) == 0) THEN
           lbnd=lbnd+1
           dclambda = - a(lbnd) / c(lbnd)
           dclambda1 = CONJG(dclambda)
           !
           !    move to new position
           !
           CALL zaxpy (ndmx, dclambda, h(1,ibnd), 1, dpsi(1,ibnd), 1)
           !
           !    update to get the gradient
           !
           !g=g+lam
           CALL zaxpy (ndmx, dclambda, t(1,lbnd), 1, g(1,ibnd), 1)
           CALL zaxpy (ndmx, dclambda1, ts(1,lbnd), 1, gs(1,ibnd), 1)
           !
           !    save current (now old) h and rho for later use
           !
           CALL zcopy (ndmx, h(1,ibnd), 1, hold(1,ibnd), 1)
           CALL zcopy (ndmx, hs(1,ibnd), 1, hsold(1,ibnd), 1)
           rhoold (ibnd) = rho (ibnd)
        ENDIF
     ENDDO
  ENDDO
100 CONTINUE
  kter = kter_eff
  DEALLOCATE (indb)
  DEALLOCATE (rho, rhoold)
  DEALLOCATE (conv)
  DEALLOCATE (a,c)
  DEALLOCATE (gs, ts, hs, hsold)
  DEALLOCATE (g, t, h, hold)

  CALL stop_clock ('ccg_many_vectors')
  RETURN
END SUBROUTINE ccg_many_vectors
!
!-------------------------------------------------------------------
SUBROUTINE linsolvms(hc,m,n,vc,alpha)
!-------------------------------------------------------------------
!
!    This routine is a driver for the correspondent lapack routines
!    which solve a linear system of equations with real coefficients. 
!    The system is assumed to be overdetermined, so the number of equation
!    is larger than the number of unknown. The program gives the solution
!    that minimize the || Ax - b ||^2. 
!    input the matrix is contained in hc, and the known part in vc, on 
!    output the solution is on alpha.
!
!
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER, INTENT(IN)  :: m, n        ! input: logical dimensions of hc

REAL(DP), INTENT(IN)  ::  hc(m,n),  &  ! input: the matrix to solve
                          vc(m)        ! input: the known part of the system

REAL(DP), INTENT(INOUT) ::  alpha(n)     ! output: the solution

REAL(DP) :: aux(m,1)                     ! auxiliary space
REAL(DP), ALLOCATABLE    :: work(:)
REAL(DP) :: rwork
INTEGER :: iwork
INTEGER :: info

aux(1:m,1)=vc(1:m)
iwork=-1
CALL dgels('N',m,n,1,hc,m,aux,m,rwork,iwork,info)
CALL errore('linsolvms','error finding optimal size',abs(info))
iwork=NINT(rwork)
ALLOCATE(work(iwork))
CALL dgels('N',m,n,1,hc,m,aux,m,work,iwork,info)
CALL errore('linsolvms','error in solving',abs(info))
alpha(1:n)=aux(1:n,1)
 
DEALLOCATE( work )

RETURN
END SUBROUTINE linsolvms

!-------------------------------------------------------------------
SUBROUTINE linsolvsvd(hc,m,n,vc,alpha)
!-------------------------------------------------------------------
!
!    This routine is a driver for the correspondent lapack routines
!    which solve a linear system of equations with real coefficients. 
!    The system is assumed to be overdetermined, so the number of equation
!    is larger than the number of unknown. The program gives the solution
!    that minimize the || Ax - b ||^2 using the singular value decomposition
!    method.
!    input the matrix is contained in hc, and the known part in vc, on 
!    output the solution is on alpha.
!
!
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER, INTENT(IN)  :: m, n        ! input: logical dimensions of hc

REAL(DP), INTENT(IN)  ::  hc(m,n),  &  ! input: the matrix to solve
                          vc(m)        ! input: the known part of the system

REAL(DP), INTENT(INOUT) ::  alpha(n)     ! output: the solution

REAL(DP) :: aux(m,1), w2(m)               ! auxiliary space
REAL(DP), ALLOCATABLE    :: work(:)
REAL(DP) :: rwork, rcond
INTEGER :: iwork, rank
INTEGER :: info

aux(1:m,1)=vc(1:m)
rcond=-1.0_DP
iwork=-1
CALL dgelss(m,n,1,hc,m,aux,m,w2,rcond,rank,rwork,iwork,info)
CALL errore('linsolvsvd','error finding optimal size',abs(info))
iwork=NINT(rwork)
ALLOCATE(work(iwork))
CALL dgelss(m,n,1,hc,m,aux,m,w2,rcond,rank,work,iwork,info)
CALL errore('linsolvsvd','error in solving',abs(info))
alpha(1:n)=aux(1:n,1)
WRITE(stdout,'(/5x,"In linsolvsvd m, n, and rank are: ",3i6)') m, n, rank
 
DEALLOCATE( work )

RETURN
END SUBROUTINE linsolvsvd

!
!-------------------------------------------------------------------
SUBROUTINE linsolvx(hc,n,vc,alpha)
!-------------------------------------------------------------------
!
!    This routine was initially in the Quantum ESPRESSO distribution.
!    
!    This routine is a driver for the correspondent lapack routines
!    which solve a linear system of equations with real coefficients. On 
!    input the matrix is contained in hc, and the known part in vc, on 
!    output the solution is on alpha.
!
!
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER, INTENT(IN)  :: n        ! input: logical dimension of hc

REAL(DP), INTENT(IN)  ::  hc(n,n),  &  ! input: the matrix to solve
                          vc(n)        ! input: the known part of the system

REAL(DP), INTENT(INOUT) ::  alpha(n)     ! output: the solution

INTEGER, ALLOCATABLE    :: iwork(:)
INTEGER :: info

ALLOCATE(iwork(n))

CALL dgetrf(n,n,hc,n,iwork,info)
CALL errore('linsolvx','error in factorization',abs(info))
alpha=vc
CALL dgetrs('N',n,1,hc,n,iwork,alpha,n,info)
CALL errore('linsolvx','error in solving',abs(info))
 
DEALLOCATE( iwork )

RETURN
END SUBROUTINE linsolvx
!
!-------------------------------------------------------------------
SUBROUTINE linsolvx_sym(hc,n,vc,alpha)
!-------------------------------------------------------------------
!
!    This routine is a driver for the corresponding lapack routine
!    which solves a linear system of equations with real coefficients
!    for a symmetric matrix hc. On input the matrix is contained in the
!    upper triangular part of hc, the right hand side in vc, on output 
!    the solution is on alpha.
!
!
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER, INTENT(IN)  :: n        ! input: logical dimension of hc

REAL(DP), INTENT(IN)  ::  hc(n,n),  &  ! input: the matrix to solve
                          vc(n)        ! input: the known part of the system

REAL(DP), INTENT(INOUT) ::  alpha(n)     ! output: the solution

REAL(DP) :: work(n)
INTEGER :: ipiv(n)
INTEGER :: info

alpha=vc
CALL dsysv('U',n,1,hc,n,ipiv,alpha,n,work,n,info)
CALL errore('linsolvx_sym','error in factorization',abs(info))
 
RETURN
END SUBROUTINE linsolvx_sym

SUBROUTINE min_sqr_solve(ndata, ncoeff, amat, f, coeff, lsolve)
!
!   This routine receives a linear system where in general the number of
!   rows ndata is larger than the number column ncoeff and finds the
!   solution of the system that minimize the norm ||AX-f||^2.
!   It can do it with three different algorithms according to the
!   value of lsolve.
!   A X = f      where A is a matrix of dimension ndata x ncoeff
!                Using 1 a matrix ncoeff x ncoeff is calculated as
!   A^T A X = A^T f   and the linear system provides X, the ncoeff 
!                coefficients of the polynomial.
!                Using 2 the overdetemined linear system AX=f is solved
!                using QR or LQ factorization.
!                Using 3 the overdetermined linear system AX=f is solved
!                using SVD decomposition.
!   If lsolve is not one of these values method 2 is used.
!
USE kinds, ONLY : DP
IMPLICIT NONE

INTEGER, INTENT(IN) :: ndata, ncoeff, lsolve
REAL(DP), INTENT(IN) :: amat(ndata,ncoeff), f(ndata) 
REAL(DP), INTENT(INOUT) :: coeff(ncoeff)

REAL(DP), ALLOCATABLE :: aa(:,:), b(:)
INTEGER :: idata, ivar, jvar, lsolve_

ALLOCATE(aa(ncoeff,ncoeff))
ALLOCATE(b(ncoeff))

aa=0.0_DP
b =0.0_DP
DO ivar=1,ncoeff
   DO jvar=1,ncoeff
      DO idata=1,ndata
         aa(ivar,jvar)= aa(ivar,jvar) + amat(idata,ivar) * amat(idata,jvar)
      END DO
   END DO
   DO idata=1,ndata
      b(ivar) = b(ivar) + amat(idata,ivar) * f(idata)
   END DO
END DO
!
!   solve the linear system and find the coefficients
!
coeff=0.0_DP
lsolve_=lsolve
IF (lsolve_<1.OR.lsolve_>3) lsolve_=2
IF (lsolve_==1) THEN
   WRITE(stdout,'(5x,"Finding the quartic polynomial using &
                                              &ncoeff x ncoeff matrix")')  
   CALL linsolvx(aa,ncoeff,b,coeff)
!   CALL linsolvx_sym(aa,ncoeff,b,coeff)

ELSEIF(lsolve_==2) THEN
   WRITE(stdout,'(5x,"Finding the quartic polynomial using &
                                                   &QR factorization")')  
   CALL linsolvms(amat,ndata,ncoeff,f,coeff)
ELSEIF(lsolve_==3) THEN
   WRITE(stdout,'(5x,"Finding the quartic polynomial using SVD")')  
   CALL linsolvsvd(amat,ndata,ncoeff,f,coeff)
ENDIF

DEALLOCATE(aa)
DEALLOCATE(b)

RETURN
END SUBROUTINE min_sqr_solve

END MODULE linear_solvers
