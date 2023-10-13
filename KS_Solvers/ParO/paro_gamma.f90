!
! Copyright (C) 2015-2016 Aihui Zhou's group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------------
!
! We propose some parallel orbital updating based plane wave basis methods
! for electronic structure calculations, which aims to the solution of the corresponding eigenvalue
! problems. Compared to the traditional plane wave methods, our methods have the feature of two level
! parallelization, which make them have great advantage in large-scale parallelization.
!
! The approach following Algorithm is the parallel orbital updating algorithm:
! 1. Choose initial $E_{\mathrm{cut}}^{(0)}$ and then obtain $V_{N_G^{0}}$, use the SCF method to solve
!    the Kohn-Sham equation in $V_{G_0}$ and get the initial $(\lambda_i^{0},u_i^{0}), i=1, \cdots, N$ 
!    and let $n=0$.
! 2. For $i=1,2,\ldots,N$, find $e_i^{n+1/2}\in V_{G_n}$ satisfying
!    $$a(\rho_{in}^{n}; e_i^{n+1/2}, v) = -[(a(\rho_{in}^{n}; u_i^{n}, v) - \lambda_i^{n} (u_i^{n}, v))]  $$
!    in parallel , where $\rho_{in}^{n}$ is the input charge density obtained by the orbits obtained in the 
!    $n$-th iteration or the former iterations.
! 3. Find $\{\lambda_i^{n+1},u_i^{n+1}\} \in \mathbf{R}\times \tilde{V}_N$   satisfying
!      $$a(\tilde{\rho}; u_i^{n+1}, v) = ( \lambda_i^{n+1}u_i^{n+1}, v) \quad  \forall v \in \tilde{V}_N$$
!      where $\tilde{V}_N = \mathrm{span}\{e_1^{n+1/2},\ldots,e_N^{n+1/2},u_1^{n},\ldots,u_N^{n}\}$, 
!      $\tilde{\rho}(x)$ is the input charge density obtained from the previous orbits.
! 4. Convergence check: if not converged, set $n=n+1$, go to step 2; else,  stop.
!
! You can see the detailed information through
!  X. Dai, X. Gong, A. Zhou, J. Zhu,
!   A parallel orbital-updating approach for electronic structure calculations, arXiv:1405.0260 (2014).
! X. Dai, Z. Liu, X. Zhang, A. Zhou,
!  A Parallel Orbital-updating Based Optimization Method for Electronic Structure Calculations, 
!   arXiv:1510.07230 (2015).
! Yan Pan, Xiaoying Dai, Xingao Gong, Stefano de Gironcoli, Gian-Marco Rignanese, and Aihui Zhou,
!  A Parallel Orbital-updating Based Plane Wave Basis Method. J. Comp. Phys. 348, 482-492 (2017).
!
! The file is written mainly by Stefano de Gironcoli and Yan Pan.
!
!-------------------------------------------------------------------------------
SUBROUTINE paro_gamma( h_psi_ptr, s_psi_ptr, hs_1psi_ptr, g_1psi_ptr, overlap, &
                 npwx, npw, nbnd, evc, eig, btype, ethr, notconv, nhpsi )
  !-------------------------------------------------------------------------------
  !paro_flag = 1: modified parallel orbital-updating method

  ! global variables
  USE util_param,          ONLY : DP, stdout
  USE mp_bands_util,       ONLY : my_bgrp_id, inter_bgrp_comm
  USE mp,                  ONLY : mp_sum

  IMPLICIT NONE
  !
  INCLUDE 'laxlib.fh'

  ! I/O variables
  LOGICAL, INTENT(IN)        :: overlap
  INTEGER, INTENT(IN)        :: npw, npwx, nbnd
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx,nbnd)
  REAL(DP), INTENT(IN)       :: ethr
  REAL(DP), INTENT(INOUT)    :: eig(nbnd)
  INTEGER, INTENT(IN)        :: btype(nbnd)
  INTEGER, INTENT(OUT)       :: notconv, nhpsi
!  INTEGER, INTENT(IN)        :: paro_flag
  
  ! local variables (used in the call to cegterg )
  !------------------------------------------------------------------------
  EXTERNAL h_psi_ptr, s_psi_ptr, hs_1psi_ptr, g_1psi_ptr
  ! subroutine h_psi_ptr  (npwx,npw,nvec,evc,hpsi)  computes H*evc  using band parallelization
  ! subroutine s_psi_ptr  (npwx,npw,nvec,evc,spsi)  computes S*evc  using band parallelization
  ! subroutine hs_1psi_ptr(npwx,npw,evc,hpsi,spsi)  computes H*evc and S*evc for a single band
  ! subroutine g_1psi_ptr (npwx,npw,psi,eig)       computes g*psi  for a single band

  !
  ! ... local variables
  !
  INTEGER :: iter, itry, paro_ntr, nconv, nextra, nactive, nbase, ndiag, nproc_ortho
  REAL(DP), ALLOCATABLE    :: ew(:)
  COMPLEX(DP), ALLOCATABLE :: psi2(:,:)
  LOGICAL, ALLOCATABLE     :: conv(:)

  INTEGER :: ibnd, ibnd_start, ibnd_end, lbnd
  !
  ! ... init local variables
  !
  CALL laxlib_getval( nproc_ortho = nproc_ortho )
  iter = 0
  nhpsi = 0
  paro_ntr = 20
  !
!  write (6,*) ' paro_flag = ', paro_flag
!  if (paro_flag /= 1) WRITE(stdout,*) 'wrong setting of paro_flag!! '

  ALLOCATE ( psi2(npwx,2*nbnd), ew(2*nbnd), conv(nbnd) )

  conv(:) =  .FALSE. ; nconv = COUNT ( conv(:) )
  psi2(:,1:nbnd) = evc(:,1:nbnd) ! copy input evc into work vector
  ew(1:nbnd)     = eig(1:nbnd)   ! copy input eigenvalues into work vector

  ParO_loop : &
  DO itry = 1,paro_ntr

     !write (6,*) ' paro_itry =', itry, ethr

     nactive = nbnd - (nconv+1)/2 ! number of correction vectors to be computed (<nbnd)
     notconv = nbnd - nconv       ! number of needed roots
     nextra  = nactive - notconv  ! number of extra vectors
     nbase   = nconv + nactive    ! number of orbitals the correction should be orthogonal to (<2*nbnd)
     ndiag   = nbase + nactive    ! dimension of the matrix to be diagonalized at this iteration (<2*nbnd)

     !write (*,*) itry, notconv, conv
     !write (6,*) ' nbnd, nconv, notconv, nextra, nactive, nbase, ndiag  =', nbnd, nconv, notconv, nextra, nactive, nbase, ndiag
     
     call s_psi_ptr  (npwx,npw,nbnd,psi2,evc) ! computes S*psi needed to ortogonalize to nbase
     lbnd = nbase
     DO ibnd = 1, nbnd ! pack unconverged roots
        IF (.NOT.conv(ibnd) ) THEN
           lbnd = lbnd+1
           psi2(:,lbnd) = psi2(:,ibnd)
           eig(lbnd-nbase) = ew(ibnd)
        END IF
     END DO
     DO ibnd = nbnd+1, nbase
        lbnd = lbnd + 1
        psi2(:,lbnd) = psi2(:,ibnd)
        eig(lbnd-nbase) = eig(lbnd-nbase-1)
     END DO
   
     !write (6,*) ' check nactive = ', lbnd-nbase, nactive
     if (lbnd .ne. nbase+nactive ) stop ' nactive check FAILED '

     CALL divide(inter_bgrp_comm,nactive,ibnd_start,ibnd_end)
     IF ( ibnd_start > 1  ) psi2(:, nbase+1:nbase+ibnd_start-1 ) = (0.0_dp,0.0_dp)
     DO ibnd=ibnd_start,ibnd_end
        !write (*,*) ' calling pcg for ibnd = ', ibnd, eig(ibnd)
        CALL pcg_gamma(hs_1psi_ptr, g_1psi_ptr, psi2, evc, npw, npwx, nbnd, psi2(:,nbase+ibnd), ethr, iter, eig(ibnd), nhpsi)
     END DO
     IF ( ibnd_end < nactive ) psi2(:, nbase+ibnd_end+1:nbase+nactive) = (0.0_dp,0.0_dp)
     CALL mp_sum(psi2(:,nbase+1:nbase+nactive),inter_bgrp_comm)

     eig(1:nbnd) = ew(1:nbnd)    ! reset first nbnd eigenvalues in their order
#if defined(__MPI)
     IF ( nproc_ortho == 1 ) THEN
#endif
        CALL rotate_wfc_gamma ( h_psi_ptr, s_psi_ptr, overlap, npwx, npw, ndiag, ndiag, psi2, psi2, ew )
#if defined(__MPI)
     ELSE
        CALL protate_wfc_gamma( h_psi_ptr, s_psi_ptr, overlap, npwx, npw, ndiag, ndiag, psi2, psi2, ew )
     ENDIF
#endif

     IF (my_bgrp_id==0) nhpsi = nhpsi + ndiag

     ! only the first nbnd eigenvalues are relevant for convergence
     conv(1:nbnd) = ABS(ew(1:nbnd)-eig(1:nbnd)).LT.ethr ; nconv = COUNT(conv(:)) ; notconv = nbnd - nconv
     IF ( nconv == nbnd ) EXIT ParO_loop

  END DO ParO_loop

  evc(:,1:nbnd) = psi2(:,1:nbnd)
  eig(1:nbnd)   = ew(1:nbnd)

  CALL mp_sum(nhpsi,inter_bgrp_comm)

  DEALLOCATE ( ew, conv, psi2 )

END SUBROUTINE paro_gamma
