!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine gmressolve_all (h_psi, cg_psi, e, d0psi, dpsi, h_diag, &
     ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd, m)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the linear system by GMRES(m) method:
  !
  !                 ( h - e + Q ) * dpsi = d0psi                      (1)
  !
  !                 where h is a complex hermitean matrix, e is a complex sca
  !                 dpsi and d0psi are complex vectors
  !
  !     on input:
  !                 h_psi    EXTERNAL  name of a subroutine:
  !                          h_psi(ndim,psi,psip)
  !                          Calculates  H*psi products.
  !                          Vectors psi and psip should be dimensined
  !                          (ndmx,nvec). nvec=1 is used!
  !
  !                 cg_psi   EXTERNAL  name of a subroutine:
  !                          g_psi(ndmx,ndim,notcnv,psi,e)
  !                          which calculates (h-e)^-1 * psi, with
  !                          some approximation, e.g. (diag(h)-e)
  !
  !                 e        complex     unperturbed eigenvalue plus
  !                          imaginary frequency.
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
  !                 m        integer  # of basis vectors
  !
  !     on output:  dpsi     contains the refined estimate of the
  !                          solution vector.
  !
  !                 d0psi    is corrupted on exit
  !
  !   revised (extensively)       6 Apr 1997 by A. Dal Corso & F. Mauri
  !   revised (to reduce memory) 29 May 2004 by S. de Gironcoli
  !
  USE kinds, only : DP
  USE mp_bands, ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum

  implicit none
  !
  !   first the I/O variables
  !
  integer :: ndmx, & ! input: the maximum dimension of the vectors
             ndim, & ! input: the actual dimension of the vectors
             kter, & ! output: counter on iterations
             nbnd, & ! input: the number of bands
             ik,   & ! input: the k point
             m       ! # of basic vector

  real(kind=DP) :: &
             anorm,   & ! output: the norm of the error in the solution
             ethr       ! input: the required precision
  complex(kind=DP) ::  h_diag(ndmx,nbnd) ! input: an estimate of ( H - \epsilon - iu )

  complex(kind=DP) :: &
             e(nbnd), & ! input: the actual eigenvalue plus imaginary freq.
             dpsi (ndmx, nbnd), & ! output: the solution of the linear syst
             d0psi (ndmx, nbnd)   ! input: the known term

  logical :: conv_root ! output: if true the root is converged
  external h_psi       ! input: the routine computing h_psi
  external cg_psi      ! input: the routine computing cg_psi
  !
  !  here the local variables
  !
  integer, parameter :: maxter = 5000
  ! the maximum number of iterations
  integer :: iter, ibnd, i, j, bnd
  ! counters on iteration, bands
  ! control variables
  integer , allocatable :: conv (:)
  ! if 1 the root is converged

  complex(kind=DP), allocatable :: r (:,:), v(:,:,:), w (:,:)!, zz(:,:), p(:,:), pp(:,:)
  !  the gradient of psi
  !  the preconditioned gradient
  !  the delta gradient
  !  the conjugate gradient
  !  work space
  complex(kind=DP) ::  bk, ak
  !  the ratio between rho
  !  step length
  complex(kind=DP), external ::  zdotc
  !  the scalar product
  real(kind=DP) :: t
  complex(kind=DP):: c, s, ei
  real(kind=DP), allocatable :: bet (:)
  real(kind=DP), allocatable :: res (:)
  complex(kind=DP) :: hm (m+1,m),   & ! the Hessenberg matrix
                      e1(m+1)          ! unit vector
  complex(kind=DP) :: hm4para(1)      ! temp variable for hm in paralell calculation
!  real(kind=DP), allocatable :: rho (:), rhoold (:), eu (:), a(:), c(:)
  ! the residue
  ! auxiliary for h_diag
  real(kind=DP) :: kter_eff
  ! account the number of iterations with b
  ! coefficient of quadratic form
  !
  integer :: lbnd
  !
  !
  !
  call start_clock ('gmres_solve')
  !
  if (m .lt. 1) then
     write(*,*) '# of basis vectors is less than 1. Stop'
     stop
  else if (m .gt. 30) then
     write(*,*) '# of basis vectors is too large. Stop'
     stop
  endif
  !
  allocate ( r(ndmx,nbnd), v(ndmx,nbnd,m+1), w(ndmx,nbnd))
  allocate (conv ( nbnd))
  allocate (bet(nbnd), res(nbnd))
  !      WRITE( stdout,*) g,t,h,hold

  kter_eff = 0.d0
  do ibnd = 1, nbnd
     conv (ibnd) = 0
  enddo
  !
  do iter = 1, maxter
     !
!print*, 'iter=', iter
     do ibnd = 1, nbnd   ! loop over bands
        !
        if (conv(ibnd) .eq. 0) then
           !
           !    preliminary step to construct the basis set
           !
           !  r = H*dpsi
           call h_psi (ndim, dpsi(1,ibnd), r(1,ibnd), e(ibnd), ik, 1)
!print*,'dpsi',sum(dpsi),sum(d0psi)
           !
           ! r = H*dpsi - d0psi
           call zaxpy (ndim, (-1.d0,0.d0), d0psi(1,ibnd), 1, r(1,ibnd), 1)
!print*,'r1',sum(dpsi),sum(d0psi)
           ! change the size of r : r = d0psi - H*dpsi
           call dscal (2 * ndim, - 1.d0, r (1, ibnd), 1)
!print*,'r2',sum(dpsi),sum(d0psi)
           ! compute the preconditioned r  : r = M^-1*r
           call cg_psi(ndmx, ndim, 1, r(1,ibnd), h_diag(1,ibnd), 1 )
!print*,'r3',sum(dpsi),sum(d0psi)
           ! norm of pre. r : bet = |r|
           bet(ibnd) = zdotc (ndim, r(1,ibnd), 1, r(1,ibnd), 1)
#if defined(__MPI)
           call mp_sum ( bet(ibnd), intra_bgrp_comm  )
#endif
           bet(ibnd) = sqrt( bet(ibnd) )
           !
        endif
        !
     enddo
     !
     !   check the convergence
     !
     lbnd = 0
     do ibnd = 1, nbnd
        !
        if ( conv(ibnd) .eq. 0 ) then
           lbnd = lbnd + 1
!if (mod(iter,10) .eq. 0) print*, iter, bet(ibnd), ethr
           if (bet(ibnd) .lt. ethr) conv(ibnd) = 1
        endif
        !
     enddo
     kter_eff = kter_eff + DBLE (lbnd) / DBLE (nbnd)
     !
     conv_root = .true.
     do ibnd = 1, nbnd
        conv_root = conv_root .and. (conv (ibnd) .eq. 1)
     enddo
     if (conv_root) goto 100
     !
     !
     !
     do ibnd = 1, nbnd
       !
       if ( conv(ibnd) .eq. 0 ) then
        !
        hm (:,:) = (0.d0, 0.d0)
        ! normalize pre. r and keep in v(1)
        call dscal (2 * ndim, 1.d0/bet(ibnd), r (1, ibnd), 1)
        j = 1
        call zcopy (ndim, r (1, ibnd), 1, v (1, ibnd, j), 1)
!print*,'v',sum(r(1:ndim,ibnd))
        !
        !
        !   loop to construct basis set
        !
        !
        do j = 1, m
           ! w = A*v
           call h_psi (ndim, v(1,ibnd,j), w(1,ibnd), e, ik, 1)       ! NEED to be checked
!print*,'w1',sum(w(:,ibnd))
           !
           ! compute w = M^-1*A*v
           call cg_psi(ndmx, ndim, 1, w(1,ibnd), h_diag(1,ibnd), 1 )
!print*,'w2',sum(w(:,ibnd))
!print*,'h_diag',sum(h_diag)
           !
           do i = 1, j
              !
              ! compute hm(i,j)
!              hm(i,j) = zdotc (ndim, w(1,ibnd), 1, v(1,ibnd,i), 1)
              hm4para(1) = zdotc (ndim, w(1,ibnd), 1, v(1,ibnd,i), 1)
#if defined(__MPI)
              call mp_sum ( hm4para, intra_bgrp_comm )
#endif
              hm(i,j) = hm4para(1)
              ! w = w - hm_ij*v_i
              call zaxpy (ndim, -hm(i,j), v(1,ibnd,i), 1, w(1,ibnd), 1)
              !
           enddo
           !   compute hm(j+1,j)
!           hm(j+1,j) = zdotc (ndim, w(1,ibnd), 1, w(1,ibnd), 1)
           hm4para(1) = zdotc (ndim, w(1,ibnd), 1, w(1,ibnd), 1)
#if defined(__MPI)
           call mp_sum ( hm4para, intra_bgrp_comm )
#endif
           hm(j+1,j) = hm4para(1)
           !   compute v(j+1)
           call dscal (2 * ndim, 1.d0/real(hm(j+1,j)), w (1, ibnd), 1)
           call zcopy (ndim, w (1, ibnd), 1, v (1, ibnd, j+1), 1)
           !
        enddo
        !
        !   compute ym that minimize |beta*e_1 -hm*y|
        !
        !   initilize vector e1
        e1(1) = 1.d0 * bet(ibnd)
        e1(2:m+1) = 0.d0
        !
        !  transform hm to upper triangle matrix
        do i = 1, m
           !
           t = sqrt( abs(hm(i,i))**2 + abs(hm(i+1,i))**2 )
           c = hm(i,i) / t
           s = hm(i+1,i) / t
           !
           do j = i, m
              !
              ei = hm(i,j)
              hm(i,j) = hm(i,j) * c + hm(i+1,j) * s
              hm(i+1,j) = - s * ei + c * hm(i+1,j)
           enddo
           !
           ei = e1(i)
           e1(i) = e1(i)*c + e1(i+1)*s
           e1(i+1) = - ei*s + e1(i+1)*c
           !
        enddo
        !
        res(ibnd) = e1(m+1)
        !
        !  back subtitution to find ym (kept in e1)
        e1(m+1) = (0.d0, 0.d0)
        e1(m) = e1(m) / hm(m,m)
        !
        do i = m-1, 1, -1
           do j = m, i+1, -1
              e1(i) = e1(i) - e1(j)*hm(i,j)
           enddo
           e1(i) = e1(i) / hm(i,i)
        enddo
        !
        !   compute the new dpsi
        do i = 1, m
           do j = 1, ndmx
              dpsi(j, ibnd) = dpsi(j, ibnd) + e1(i)*v(j,ibnd,i)
           enddo
        enddo
        !
       end if
       !
     enddo   ! of loop over bands
     !
  enddo   ! loop over iteration
  !
100 continue
  kter = kter_eff
  !
  deallocate (bet, res)
  deallocate (conv)
  deallocate (r, v, w)
  !
  call stop_clock ('gmres_solve')
  !
  return
  !
end subroutine gmressolve_all
