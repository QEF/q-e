!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine cgsolve_all (h_psi, cg_psi, e, d0psi, dpsi, h_diag, &
     ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the linear system:
  !
  !                 ( h - e + Q ) * dpsi = d0psi                      (1)
  !
  !                 where h is a complex hermitean matrix, e is a real sca
  !                 dpsi and d0psi are complex vectors
  !
  !     on input:
  !                 h_psi             EXTERNAL  name of a subroutine:
  !                                   h_psi(ndim,psi,psip)
  !                                   Calculates  H*psi products.
  !                                   Vectors psi and psip shoul be dimens
  !                                   (ndmx,nvec). nvec=1 is used!
  !
  !                 cg_psi             EXTERNAL  name of a subroutine:
  !                                   g_psi(ndmx,ndim,notcnv,psi,e)
  !                                   which calculates (h-e)^-1 * psi, wit
  !                                   some approximation, e.g. (diag(h)-e)
  !
  !                 e                 real     unperturbed eigenvalue.
  !
  !                 dpsi              contains an estimate of the solution
  !                                   vector.
  !
  !                 d0psi             contains the right hand side vector
  !                                   of the system.
  !
  !                 ndmx              integer row dimension of dpsi, ecc.
  !
  !                 ndim              integer actual row dimension of dpsi
  !
  !                 ethr              real     convergence threshold. solu
  !                                   improvement is stopped when the erro
  !                                   eq (1), defined as l.h.s. - r.h.s.,
  !                                   less than ethr in norm.
  !
  !     on output:  dpsi              contains the refined estimate of the
  !                                   solution vector.
  !
  !                 d0psi             is corrupted on exit
  !
  !
  !   Last revised 6 apr. 1997 by A. Dal Corso & F. Mauri
  !
  !
#include "machine.h"
  use parameters, only : DP
  implicit none
  !
  !   first the dummy variables
  !

  integer :: ndmx, ndim, kter, nbnd, ik
  ! input: the maximum dimension of the vectors
  ! input: the actual dimension of the vectors
  ! output: counter on iterations
  ! input: the number of bands
  ! input: the k point

  real(kind=DP) :: e (nbnd), anorm, h_diag (ndmx, nbnd), ethr
  ! input: the actual eigenvalue
  ! output: the norm of the error in the solution
  ! input: an estimate of ( H - \epsilon
  ! input: the required precision

  complex(kind=DP) :: dpsi (ndmx, nbnd), d0psi (ndmx, nbnd)
  ! output: the solution of the linear syst
  ! input: the known term

  logical :: conv_root
  ! output: if true the root is converged
  external h_psi, cg_psi
  ! input: the routine computing h_psi
  ! input: the routine computing cg_psi
  !
  !  three parameters
  !

  integer :: maxter
  ! the maximum number of iterations
  parameter (maxter = 200)
  !
  !  here the local variables
  !

  integer :: iter, ibnd, lbnd
  integer , allocatable :: conv (:)
  ! counter on iteration
  ! counter on bands
  ! if 1 the root is converged

  complex(kind=DP), allocatable :: g (:,:), gp (:,:), t (:,:), &
       h (:,:), hold (:,:), aux (:,:), aux1 (:,:)
  complex(kind=DP) ::  dcgamma, dclambda, ZDOTC
  !  the gradient of psi
  !  the preconditioned gradient
  !  the delta gradient
  !  the conjugate gradient
  !  the the old h
  !  the the old h
  !  the the old h
  !  the ratio between rho
  !  step lenght
  !  the scalar product


  real(kind=DP), allocatable :: rho (:), rhoold (:), auxr (:,:), eu (:)
  real(kind=DP) :: kter_eff, a, c
  ! the residue
  ! auxiliary for h_diag
  ! account the number of iterations with b
  ! coefficient of quadratic form
  !
  call start_clock ('cgsolve')
  allocate (g (  ndmx , nbnd))    
  allocate (gp ( ndmx , nbnd))    
  allocate (t ( ndmx , nbnd))    
  allocate (h (  ndmx , nbnd))    
  allocate (hold ( ndmx , nbnd))    
  allocate (aux1 ( ndmx , nbnd))    
  allocate (aux (  ndmx , nbnd))    
  allocate (conv ( nbnd))    
  allocate (rho (  nbnd))    
  allocate (rhoold (  nbnd))    
  allocate (auxr (  ndmx , nbnd))    
  allocate (eu (  nbnd))    
  !      write(6,*) g,gp,t,h,hold
  kter = 0

  kter_eff = 0.d0
  do ibnd = 1, nbnd
     conv (ibnd) = 0

  enddo
  do iter = 1, maxter
     !
     !        compute the gradient. can reuse information from previous step
     !
     !
     kter = kter + 1
     if (iter.eq.1) then
        call h_psi (ndim, dpsi, g, e, ik, nbnd)
        do ibnd = 1, nbnd
           call ZAXPY (ndim, ( - 1.d0, 0.d0), d0psi (1, ibnd), 1, g (1, &
                ibnd), 1)
        enddo
     endif
     !
     !        compute residual
     !
     lbnd = 0
     do ibnd = 1, nbnd
        if (conv (ibnd) .eq.0) then
           call ZCOPY (ndim, g (1, ibnd), 1, gp (1, ibnd), 1)
           lbnd = lbnd+1
           call ZCOPY (ndim, gp (1, ibnd), 1, aux (1, lbnd), 1)
           call DCOPY (ndmx, h_diag (1, ibnd), 1, auxr (1, lbnd), 1)
        endif
     enddo
     call cg_psi (ndmx, ndim, lbnd, aux, auxr)
     kter_eff = kter_eff + float (lbnd) / float (nbnd)
     lbnd = 0
     do ibnd = 1, nbnd
        if (conv (ibnd) .eq.0) then
           lbnd = lbnd+1
           call ZCOPY (ndim, aux (1, lbnd), 1, gp (1, ibnd), 1)
        endif

     enddo
     do ibnd = 1, nbnd
        if (conv (ibnd) .eq.0) then
           rho (ibnd) = ZDOTC (ndim, gp (1, ibnd), 1, g (1, ibnd), &
                1)
#ifdef __PARA

           call reduce (1, rho (ibnd) )
#endif
           anorm = sqrt (rho (ibnd) )
           !               write(6,'(2i5,e20.5)') iter,ibnd,anorm
#ifdef FLUSH
           !               call flush(6)
#endif
           if (anorm.lt.ethr) conv (ibnd) = 1
        endif
     enddo
     conv_root = .true.
     do ibnd = 1, nbnd
        conv_root = conv_root.and. (conv (ibnd) .eq.1)
     enddo
     if (conv_root) goto 100
     !
     !        compute the step direction h. Conjugate it to previous step
     !
     lbnd = 0
     do ibnd = 1, nbnd
        if (conv (ibnd) .eq.0) then
           call ZCOPY (ndim, gp (1, ibnd), 1, h (1, ibnd), 1)
           call DSCAL (2 * ndim, - 1.d0, h (1, ibnd), 1)
           if (iter.ne.1) then
              dcgamma = rho (ibnd) / rhoold (ibnd)
              call ZAXPY (ndim, dcgamma, hold (1, ibnd), 1, h (1, ibnd), &
                   1)

           endif

           rhoold (ibnd) = rho (ibnd)
           call ZCOPY (ndim, h (1, ibnd), 1, hold (1, ibnd), 1)
           lbnd = lbnd+1
           call ZCOPY (ndim, h (1, ibnd), 1, aux (1, lbnd), 1)
           eu (lbnd) = e (ibnd)
        endif
     enddo
     !
     !        compute t = A*h
     !
     call h_psi (ndim, aux, aux1, eu, ik, lbnd)
     lbnd = 0
     do ibnd = 1, nbnd
        if (conv (ibnd) .eq.0) then
           lbnd = lbnd+1
           call ZCOPY (ndim, aux1 (1, lbnd), 1, t (1, ibnd), 1)
        endif


     enddo
     !
     !        compute the coefficients a and c for the line minimization
     !        compute step length lambda
     do ibnd = 1, nbnd
        if (conv (ibnd) .eq.0) then
           a = ZDOTC (ndim, h (1, ibnd), 1, g (1, ibnd), 1)
           c = ZDOTC (ndim, h (1, ibnd), 1, t (1, ibnd), 1)
#ifdef __PARA
           call reduce (1, a)
           call reduce (1, c)
#endif
           dclambda = DCMPLX ( - a / c, 0.d0)
           !
           !           move to new position
           !
           call ZAXPY (ndim, dclambda, h (1, ibnd), 1, dpsi (1, ibnd), &
                1)
           !
           !           update to get the gradient
           !
           !g=g+lam
           call ZAXPY (ndim, dclambda, t (1, ibnd), 1, g (1, ibnd), &
                1)
        endif
     enddo

  enddo
100 continue
  kter = kter_eff
  deallocate (eu)
  deallocate (auxr)
  deallocate (rho)
  deallocate (rhoold)
  deallocate (conv)
  deallocate (aux)
  deallocate (aux1)
  deallocate (hold)
  deallocate (h)
  deallocate (t)
  deallocate (gp)
  deallocate (g)

  call stop_clock ('cgsolve')
  return
end subroutine cgsolve_all
