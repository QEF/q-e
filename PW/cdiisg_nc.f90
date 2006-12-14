!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine cdiisg_nc (ndim, ndmx, nvec, nvecx, evc, e, ethr, &
      btype, notcnv, diis_iter, iter, npol)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the eigenvalue problem:
  !
  !     ( H - e S ) * evc = 0
  !
  !     where H is an hermitean operator, e is a real scalar,
  !     S is an overlap matrix, evc is a complex vector.
  !     The band-by-band RMM-DIIS method is used. 
#include "f_defs.h"
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  USE g_psi_mod
  USE uspp, ONLY: okvan

  implicit none
  ! on INPUT
  integer :: ndim, ndmx, nvec, nvecx, btype(nvec), iter, npol
  ! dimension of the matrix to be diagonalized
  ! leading dimension of matrix evc, as declared in the calling pgm unit
  ! integer number of searched low-lying roots
  ! maximum dimension of the reduced basis set
  !    (the basis set is refreshed when its dimension would exceed nvecx)
  ! band type (0=occupied, 1=empty)
  ! scf iteration
  ! number of coordinates of wfc
  real(DP) :: ethr
  ! energy threshold for convergence
  !   root improvement is stopped, when two consecutive estimates of the root
  !   differ by less than ethr.
  ! on OUTPUT
  complex(DP) :: evc (ndmx, npol, nvec)
  !  evc   contains the  refined estimates of the eigenvectors
  real(DP) :: e (nvec), diis_iter
  ! contains the estimated roots.
  ! average number of iterations performed per band
  integer :: notcnv
  ! number of unconverged roots

  ! LOCAL variables
  !
  integer, parameter :: maxter=20
  ! maximum number of iterations
  !
  integer :: kter, minter, nbase, ib, n, m 
  ! counter on iterations
  ! lower extreme for iteration loop
  ! dimension of the reduced basis
  ! counter on the reduced basis vectors
  ! do-loop counters
  complex(DP), allocatable :: rc (:,:),  hc (:,:), sc (:,:), &
       vc (:), vcn(:,:)
  ! <res_i|res_j> matrix
  ! H matrix on the reduced basis
  ! S matrix on the reduced basis
  ! the eigenvectors of the Hamiltonian
  ! workspace
  complex(DP), allocatable :: psi(:,:,:),hpsi(:,:,:),spsi(:,:,:),res(:,:,:)
  ! work space, contains psi
  ! the product of H and psi
  ! the product of S and psi
  ! residual vector
  complex(DP) :: hevc(ndmx, npol, nvec), sevc(ndmx, npol, nvec)
  ! the product of H and the best estimate of the eigenvectors evc
  ! the product of S and the best estimate of the eigenvectors evc
  real(DP), allocatable :: ew (:)
  ! eigenvalues of the reduced hamiltonian
  real(DP) :: ec, snorm, snorm0, lam, ew0
  ! dummy variable
  ! squared norm of current residual
  ! squared norm of initial residual calculated with the old evc
  ! initial eigenvalue from previous iteration
  ! variables for teter preconditioning (not used)
  logical :: verb, test_new_preconditioning_nc 
  ! controlling verbosity of printout
  complex(DP), external ::  ZDOTC
  external h_1psi, cdiagh

  call start_clock ('diis')

  test_new_preconditioning_nc = .true.
  verb = .false.

  !
  ! allocate the work arrays
  !
  allocate( psi (ndmx, npol, nvecx))
  allocate(hpsi (ndmx, npol, nvecx))
  allocate(spsi (ndmx, npol, nvecx))
  allocate( res (ndmx, npol, nvecx))
  allocate( rc (nvecx, nvecx))
  allocate( hc (nvecx, nvecx))
  allocate( sc (nvecx, nvecx))
  allocate( vc (nvecx))
  allocate( vcn(nvecx, nvecx))
  allocate( ew (nvecx))

  !
  ! rotate
  !
  call rotate_wfc_nc (ndmx, ndim, nvec, nvec, evc, npol, okvan, evc, e)

  notcnv = nvec
  lam = 1.d0
  minter = 1

  ! Loop over bands
  do ib = 1, nvec

10   continue
     !
     !     prepare the hamiltonian for the first iteration
     !
     nbase = 1     
     ew  (:) = 0.d0
     vc  (:) = (0.d0, 0.d0)
     vcn (:, :) = (0.d0, 0.d0)
     psi (:, :, :) = (0.d0, 0.d0)
     res (:, :, :) = (0.d0, 0.d0)
     hc  (:, :) = (0.d0, 0.d0)
     sc  (:, :) = (0.d0, 0.d0)
     rc  (:, :) = (0.d0, 0.d0)
     if (minter.eq.1) then
        ew0 = 1.d10
        snorm0 = 1.d-20
     endif
     ! |psi_1> is the best approximated eigenvector
     psi  = (0.d0, 0.d0)
     hpsi = (0.d0, 0.d0)
     psi(:, :, nbase) = evc(:, :, ib)
     vc(nbase) = (1.d0, 0.d0)
     !
     !     calculate hpsi=H|psi_1> and spsi=S|psi_1> 
     !
     call h_1psi (ndmx, ndim, psi(1,1,nbase), hpsi(1,1,nbase), &
                  spsi(1,1,nbase))
     !  
     !     calculate the first element of the reduced hamiltonian 
     !     and overlap matrices
     !     hc(1,1)=<psi_1|H|psi_1>    sc(1,1)=<psi_1|S|psi_1>
     !
     IF (npol == 1) THEN
        hc (1, 1) = ZDOTC (ndim, psi (1, 1, 1), 1, hpsi (1, 1, 1), 1)
        sc (1, 1) = ZDOTC (ndim, psi (1, 1, 1), 1, spsi (1, 1, 1), 1)
     ELSE
        hc (1, 1) = ZDOTC (ndmx*npol, psi (1, 1, 1), 1, hpsi (1, 1, 1), 1)
        sc (1, 1) = ZDOTC (ndmx*npol, psi (1, 1, 1), 1, spsi (1, 1, 1), 1)
     ENDIF
#ifdef __PARA
     call reduce (2, hc(1, 1))
     call reduce (2, sc(1, 1))
#endif
     !
     !   calculate the residual vector |R>=H|psi> - e S|psi>
     !
     IF (npol == 1) THEN
        call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0.d0), spsi, &
             ndmx, vc , nvecx, (0.d0, 0.d0), res (1, 1, nbase), ndmx)
        call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0d0), hpsi, &
             ndmx, vc , nvecx, CMPLX(- e(ib), 0.d0), res (1, 1, nbase), ndmx)
     ELSE
        call ZGEMM ('n', 'n', ndmx*npol, 1, nbase, (1.d0, 0.d0), spsi, &
             ndmx*npol, vc , nvecx, (0.d0, 0.d0), res (1, 1, nbase), ndmx*npol)
        call ZGEMM ('n', 'n', ndmx*npol, 1, nbase, (1.d0, 0d0), hpsi, &
                     ndmx*npol, vc , nvecx, CMPLX(- e(ib), 0.d0), & 
                     res (1, 1, nbase), ndmx*npol)
     ENDIF
     !
     !   calculate the first element of the <res_i|res_j> matrix
     IF (npol == 1) THEN
        rc (1, 1) = ZDOTC (ndim, res (1, 1, 1), 1, res (1, 1, 1), 1)
     ELSE
        rc (1, 1) = ZDOTC (ndmx*npol, res (1, 1, 1), 1, res (1, 1, 1), 1)
     ENDIF
#ifdef __PARA
     call reduce (2 , rc(1, 1))
#endif
     !
     !  iterate
     !
     do kter = minter, maxter
     !
     !   preconditionate the  residual vector |P_n>= K*|R_n>
     !   and add it to the basis => |psi_n+1>
     !
        IF (npol == 1) THEN
           call ZGEMM ('n', 'n', ndim , 1, nbase, (1.d0, 0.d0) , res, &
                ndmx, vc , nvecx, (0.d0, 0.d0) , psi (1, 1, nbase+1), ndmx)

        ELSE
           call ZGEMM ('n', 'n', ndmx*npol , 1, nbase, (1.d0, 0.d0) , res, &
                ndmx*npol, vc , nvecx, (0.d0, 0.d0) , &
                psi (1, 1, nbase+1), ndmx*npol)

        ENDIF
           call g_psi_nc( ndmx, ndim, 1, npol, psi (1, 1, nbase+1), e(ib) )

     ! 
     ! add a new vector to the basis: kresse method
     ! |psi_n+1> = |psi_n> + \lambda * K * |res_n>
!        call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0.d0), psi, &
!             ndmx, vc , nvecx, CMPLX (lam, 0.d0) , &
!             psi (1, 1, nbase+1), ndmx)

        nbase = nbase + 1
     !
     ! normalize new basis vector
        IF (npol == 1) THEN
           ec =  DBLE(ZDOTC (ndim, psi (1, 1, nbase), 1, psi (1, 1, nbase), 1))
        ELSE
           ec =  DBLE(ZDOTC(ndmx*npol,psi (1, 1,nbase),1,psi(1, 1, nbase),1))
        ENDIF
#ifdef __PARA
     call reduce (1 , ec)
#endif
        IF (npol == 1) THEN
           call ZSCAL (ndim, CMPLX (1/dsqrt(ec), 0.d0), psi (1,1,nbase), 1)
        ELSE
           call ZSCAL (ndmx*npol, CMPLX(1/dsqrt(ec),0.d0), psi (1,1,nbase), 1)
        ENDIF
     ! new eigenvector, normalize eigenvectors
        vc(nbase) = (1.d0, 0.d0)
        ec =  DBLE(ZDOTC (nbase, vc , 1, vc , 1))
        call ZSCAL (nvecx, CMPLX (1/dsqrt(ec), 0.d0), vc, 1)
     !
     !     calculate hpsi=H|psi> and spsi=S|psi> 
     !
        call h_1psi (ndmx, ndim, psi(1,1,nbase), hpsi(1,1,nbase), &
                     spsi(1,1,nbase))
     !
     !     orthogonalize
     !
           call cgramg1_nc (ndmx, nvecx, ndim, 1, nbase, psi, spsi, hpsi, npol)
     !
     !     calculate the new elements of the reduced hamiltonian
     !     and overlap matrices
     !     hc(i,j) =<psi_i|H|psi_j>  and sc(i,j)=<psi_i|S|psi_j>
     !
        IF (npol == 1) THEN
           call ZGEMM ('c', 'n', nbase, 1, ndim, (1.d0, 0.d0) , psi, &
                ndmx, hpsi (1, 1, nbase) , ndmx, (0.d0, 0.d0) , hc (1, nbase), &
                nvecx)
           call ZGEMM ('c', 'n', nbase, 1, ndim, (1.d0, 0.d0) , psi, &
                ndmx, spsi (1, 1, nbase) , ndmx, (0.d0, 0.d0) , sc (1, nbase), &
                nvecx)
        ELSE
           call ZGEMM ('c', 'n', nbase, 1, ndmx*npol, (1.d0, 0.d0) , psi, &
                ndmx*npol, hpsi (1, 1, nbase) , ndmx*npol, (0.d0, 0.d0) , &
                hc (1, nbase), nvecx)
           call ZGEMM ('c', 'n', nbase, 1, ndmx*npol, (1.d0, 0.d0) , psi, &
                ndmx*npol, spsi (1, 1, nbase) , ndmx*npol, &
                (0.d0, 0.d0) , sc (1, nbase), &
                nvecx)
        ENDIF
#ifdef __PARA
     call reduce (2 * nvecx, hc(1, nbase))
     call reduce (2 * nvecx, sc(1, nbase))
#endif
     !
     !   calculate the residual vector |R>=H|psi> - e S|psi>
     !
        IF (npol == 1) THEN
           call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0 , 0.d0), spsi, &
                ndmx, vc , nvecx, (0.d0, 0.d0), res (1, 1, nbase), ndmx)
           call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0d0), hpsi, &
                ndmx, vc , nvecx, CMPLX (-e(ib), 0.d0), res (1,1,nbase), ndmx)
        !
        !   calculate the new elements of the <res_i|res_j> matrix
           call ZGEMM ('c', 'n', nbase, 1, ndim, (1.d0, 0.d0) , &
                res , ndmx, res (1, 1, nbase) , ndmx, (0.d0, 0.d0) , &
                rc (1, nbase) , nvecx)
        ELSE
           call ZGEMM ('n', 'n', ndmx*npol, 1, nbase, (1.d0 , 0.d0), spsi, &
             ndmx*npol, vc , nvecx, (0.d0, 0.d0), res (1, 1, nbase), ndmx*npol)
           call ZGEMM ('n', 'n', ndmx*npol, 1, nbase, (1.d0, 0d0), hpsi, &
                ndmx*npol, vc , nvecx, CMPLX (-e(ib), 0.d0), &
                res (1,1,nbase), ndmx*npol)
        !
        !   calculate the new elements of the <res_i|res_j> matrix
           call ZGEMM ('c', 'n', nbase, 1, ndmx*npol, (1.d0, 0.d0) , &
                res , ndmx*npol, res (1,1,nbase) , ndmx*npol, (0.d0, 0.d0) , &
                rc (1, nbase) , nvecx)
        ENDIF
        ew(nbase) = rc (nbase, nbase)
#ifdef __PARA
     call reduce (2 * nvecx, rc(1, nbase))
#endif
     !
     !     rc, hc, and sc are hermitian
     !
     do n = 1, nbase
        !  the diagonal of rc and sc must be strictly real 
        rc (n, n) = CMPLX ( DBLE (rc (n, n) ), 0.d0)
        hc (n, n) = CMPLX ( DBLE (hc (n, n) ), 0.d0)
        sc (n, n) = CMPLX ( DBLE (sc (n, n) ), 0.d0)
        do m = n + 1, nbase
           rc (m, n) = CONJG (rc (n, m) )
           hc (m, n) = CONJG (hc (n, m) )
           sc (m, n) = CONJG (sc (n, m) )
        enddo
     enddo
     if (verb) then
        WRITE( stdout,*) 'overlap' 
        WRITE( stdout,*) ((m,n,sc(n,m), n=1,nbase), m=1,nbase)
        WRITE( stdout,*)  
        WRITE( stdout,*) 'rc' 
        WRITE( stdout,*) ((m,n,rc(n,m), n=1,nbase), m=1,nbase)
        WRITE( stdout,*)  
        WRITE( stdout,*) 'eigval'
     endif
     

     !
     !     diagonalize the reduced hamiltonian
     !
     call cdiaghg (nbase, 1, rc, sc, nvecx, ew, vcn )
     call ZCOPY (nvecx, vcn (1, 1), 1, vc, 1 )

     if (verb) then
        do n=1, nbase
           WRITE( stdout,*) n,ew(n)
        enddo
        WRITE( stdout,*)  
        do n=1,nbase
        enddo
        WRITE( stdout,*) 'eigvec' 
        do n=1, nbase
           WRITE( stdout,*) n, vc(n)
        enddo
     endif

     ! squared norm of current residual
     snorm = ew(1)
     !
     ! calculate new eigenvalues
     vcn(:,:)= (0.d0,0.d0)
     call ZGEMM ('n', 'n', nbase, 1, nbase, (1.d0 , 0.d0), hc, &
          nvecx, vc , nvecx, (0.d0, 0.d0), vcn (1, 1), nvecx)
     ec =  DBLE( ZDOTC (nvecx, vc, 1, vcn (1, 1), 1) )
     call ZGEMM ('n', 'n', nbase, 1, nbase, (1.d0 , 0.d0), sc, &
          nvecx, vc , nvecx, (0.d0, 0.d0), vcn (1, 1), nvecx)
     ec = ec /  DBLE( ZDOTC (nvecx, vc, 1, vcn (1, 1), 1) )

     if (verb) WRITE( stdout,*) 'NORM RES=',snorm,'DELTA EIG=',ec-e(ib)

     !
     ! Convergence?
     ! Non occupied levels are converged with a lower precision than 
     ! occupied ones.
     !
     if (btype(ib) .eq. 0) then 
        if ( (snorm.lt.snorm0*0.3d0 .or. abs(ec-e(ib)).lt.ethr) &
             .and. ec.le.ew0) goto 20
     else
        if ( (snorm.lt.snorm0*0.3d0 .or. abs(ec-e(ib)).lt.max(ethr*50.0d0,1.D-4)) &
             .and. ec.le.ew0) goto 20
     endif
        
     if (minter.eq.1) then
        snorm0 = snorm
        ew0 = ec
     endif
     e(ib) = ec
     
     !
     ! Size of reduced basis is exceeded: refresh
     !
     if (nbase.ge.nvecx) then
        IF (npol ==1) THEN
           call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0d0), psi, &
                ndmx, vc , nvecx, (0.d0, 0.d0), evc (1, 1, ib), ndmx)
        ELSE
           call ZGEMM ('n', 'n', ndmx*npol, 1, nbase, (1.d0, 0d0), psi, &
                ndmx*npol, vc , nvecx, (0.d0, 0.d0), evc (1, 1, ib), ndmx*npol)
        ENDIF

        if (verb) WRITE( stdout,*) 'rotate band ',ib
        minter = kter + 1
        goto 10
     endif

     enddo ! iterate
20   continue

     e(ib) = ec
     diis_iter = diis_iter + kter

     if (kter .gt. maxter) then
    !    WRITE( stdout, '("   WARNING: eigenvalue ",i5," not converged")') &
    !         ib
     else
        notcnv = notcnv - 1
        minter = 1
        if (verb) then
           WRITE( stdout,*) 'BAND ',ib, ' CONVERGED'
           WRITE( stdout,*)
        endif
     endif
     !
     ! calculate best approximated wavefunction and corresponding hpsi and spsi
     !
     IF (npol == 1) THEN
        call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0d0), psi, &
             ndmx, vc , nvecx, (0.d0, 0.d0), evc (1, 1, ib), ndmx)
        call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0d0), hpsi, &
             ndmx, vc , nvecx, (0.d0, 0.d0), hevc (1, 1, ib), ndmx)
        call ZGEMM ('n', 'n', ndim, 1, nbase, (1.d0, 0d0), spsi, &
             ndmx, vc , nvecx, (0.d0, 0.d0), sevc (1, 1, ib), ndmx)
     ELSE
        call ZGEMM ('n', 'n', ndmx*npol, 1, nbase, (1.d0, 0d0), psi, &
             ndmx*npol, vc , nvecx, (0.d0, 0.d0), evc (1, 1, ib), ndmx*npol)
        call ZGEMM ('n', 'n', ndmx*npol, 1, nbase, (1.d0, 0d0), hpsi, &
             ndmx*npol, vc , nvecx, (0.d0, 0.d0), hevc (1, 1, ib), ndmx*npol)
        call ZGEMM ('n', 'n', ndmx*npol, 1, nbase, (1.d0, 0d0), spsi, &
             ndmx*npol, vc , nvecx, (0.d0, 0.d0), sevc (1, 1, ib), ndmx*npol)
     ENDIF

  enddo ! loop over bands     

  diis_iter = diis_iter / nvec

  deallocate( psi )
  deallocate(hpsi )
  deallocate(spsi )
  deallocate( res )
  deallocate( rc )
  deallocate( hc )
  deallocate( sc )
  deallocate( vc )
  deallocate(vcn )
  deallocate( ew )

  !
  ! orthonormalize bands
  !
!  if (mod(iter,6).eq.0) &
       call cgramg1_nc (ndmx, nvec, ndim, 1, nvec, evc, sevc, hevc, npol)

  call stop_clock ('diis')
  return
end subroutine cdiisg_nc


