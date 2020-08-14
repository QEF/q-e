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
SUBROUTINE paro_gamma_new( h_psi, s_psi, hs_psi, g_1psi, overlap, &
                 npwx, npw, nbnd, evc, eig, btype, ethr, notconv, nhpsi )
  !-------------------------------------------------------------------------------
  !paro_flag = 1: modified parallel orbital-updating method

  ! global variables
  USE util_param,          ONLY : DP, stdout
  USE mp_bands_util,       ONLY : inter_bgrp_comm, nbgrp, my_bgrp_id
  USE mp,                  ONLY : mp_sum, mp_allgather, mp_barrier, &
                                  mp_type_create_column_section, mp_type_free

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
  EXTERNAL h_psi, s_psi, hs_psi, g_1psi
  ! subroutine h_psi  (npwx,npw,nvec,evc,hpsi)  computes H*evc  using band parallelization
  ! subroutine s_psi  (npwx,npw,nvec,evc,spsi)  computes S*evc  using band parallelization
  ! subroutine hs_1psi(npwx,npw,evc,hpsi,spsi)  computes H*evc and S*evc for a single band
  ! subroutine g_1psi  (npwx,npw,psi,eig)       computes g*psi  for a single band

  !
  ! ... local variables
  !
  INTEGER :: itry, paro_ntr, nconv, nextra, nactive, nbase, ntrust, ndiag, nvecx, nproc_ortho
  REAL(DP), ALLOCATABLE    :: ew(:)
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
  LOGICAL, ALLOCATABLE     :: conv(:)

  REAL(DP), PARAMETER      :: extra_factor = 0.5 ! workspace is at most this factor larger than nbnd
  INTEGER, PARAMETER       :: min_extra = 4      ! but at least this lager

  INTEGER :: ibnd, ibnd_start, ibnd_end, how_many, lbnd, kbnd, last_unconverged, &
             recv_counts(nbgrp), displs(nbgrp), column_type
  !
  ! ... init local variables
  !
  CALL laxlib_getval( nproc_ortho = nproc_ortho )
  paro_ntr = 20

  nvecx = nbnd + max ( nint ( extra_factor * nbnd ), min_extra )
  !
  CALL start_clock( 'paro_gamma' ); !write (6,*) ' enter paro diag'

  CALL mp_type_create_column_section(evc(1,1), 0, npwx, npwx, column_type)

  ALLOCATE ( psi(npwx,nvecx), hpsi(npwx,nvecx), spsi(npwx,nvecx), ew(nvecx), conv(nbnd) )

  CALL start_clock( 'paro:init' ); 
  conv(:) =  .FALSE. ; nconv = COUNT ( conv(:) )
  psi(:,1:nbnd) = evc(:,1:nbnd) ! copy input evc into work vector
  call h_psi  (npwx,npw,nbnd,psi,hpsi) ! computes H*psi
  call s_psi  (npwx,npw,nbnd,psi,spsi) ! computes S*psi
  nhpsi = 0 ; IF (my_bgrp_id==0) nhpsi = nbnd
  CALL stop_clock( 'paro:init' ); 

#if defined(__MPI)
  IF ( nproc_ortho == 1 ) THEN
#endif
     CALL rotate_HSpsi_gamma (  npwx, npw, nbnd, nbnd, psi, hpsi, overlap, spsi, eig )
#if defined(__MPI)
  ELSE
     CALL protate_HSpsi_gamma(  npwx, npw, nbnd, nbnd, psi, hpsi, overlap, spsi, eig )
  ENDIF
#endif
  !write (6,'(10f10.4)') psi(1:5,1:3)

  !write (6,*) eig(1:nbnd)

  ParO_loop : &
  DO itry = 1,paro_ntr

     !write (6,*) ' paro_itry =', itry, ethr

     !----------------------------

     nactive = nbnd - (nconv+1)/2 ! number of correction vectors to be computed (<nbnd)
     notconv = nbnd - nconv       ! number of needed roots
     nextra  = nactive - notconv  ! number of extra vectors
     nbase   = nconv + nactive    ! number of orbitals the correction should be orthogonal to (<2*nbnd)
     ndiag   = nbase + nactive    ! dimension of the matrix to be diagonalized at this iteration (<2*nbnd)

     !----------------------------

     nactive = min ( (nvecx-nconv)/2, nvecx-nbnd) ! number of corrections there is space for
     notconv = nbnd - nconv                       ! number of needed roots
     nextra  = max ( nactive - notconv, 0 )       ! number of extra vectors, if any
     nbase   = max ( nconv + nactive , nbnd )     ! number of orbitals to be orthogonal to  (<nvecx)
     ntrust  = min ( nconv + nactive , nbnd )     ! number of orbitals that will be actually corrected
     ndiag   = nbase + nactive       ! dimension of the matrix to be diagonalized at this iteration (<nvecx)

     !write (6,*) itry, notconv, conv
     !write (6,*) ' nvecx, nbnd, nconv, notconv, nextra, nactive, nbase, ntrust, ndiag  =', nvecx, nbnd, nconv, notconv, nextra, nactive, nbase, ntrust, ndiag
     
     CALL divide_all(inter_bgrp_comm,nactive,ibnd_start,ibnd_end,recv_counts,displs)
     how_many = ibnd_end - ibnd_start + 1
     !write (6,*) nactive, ibnd_start, ibnd_end, recv_counts, displs

     CALL start_clock( 'paro:pack' ); 
     lbnd = 1; kbnd = 1
     DO ibnd = 1, ntrust ! pack unconverged roots in the available space
        IF (.NOT.conv(ibnd) ) THEN
           psi (:,nbase+kbnd)  = psi(:,ibnd)
           hpsi(:,nbase+kbnd) = hpsi(:,ibnd)
           spsi(:,nbase+kbnd) = spsi(:,ibnd)
           ew(kbnd) = eig(ibnd) ; last_unconverged = ibnd
           lbnd=lbnd+1 ; kbnd=kbnd+recv_counts(mod(lbnd-2,nbgrp)+1); if (kbnd>nactive) kbnd=kbnd+1-nactive
        END IF
     END DO
     DO ibnd = nbnd+1, nbase   ! add extra vectors if it is the case
        psi (:,nbase+kbnd)  = psi(:,ibnd)
        hpsi(:,nbase+kbnd) = hpsi(:,ibnd)
        spsi(:,nbase+kbnd) = spsi(:,ibnd)
        ew(kbnd) = eig(last_unconverged)
        lbnd=lbnd+1 ; kbnd=kbnd+recv_counts(mod(lbnd-2,nbgrp)+1); if (kbnd>nactive) kbnd=kbnd+1-nactive
     END DO
     psi (:,nbase+1:nbase+how_many) = psi (:,nbase+ibnd_start:nbase+ibnd_end)
     hpsi(:,nbase+1:nbase+how_many) = hpsi(:,nbase+ibnd_start:nbase+ibnd_end)
     spsi(:,nbase+1:nbase+how_many) = spsi(:,nbase+ibnd_start:nbase+ibnd_end)
     ew(1:how_many) = ew(ibnd_start:ibnd_end)
     CALL stop_clock( 'paro:pack' ); 
   
!     write (6,*) ' check nactive = ', lbnd, nactive
     if (lbnd .ne. nactive+1 ) stop ' nactive check FAILED '


     CALL bpcg_gamma(hs_psi, g_1psi, psi, spsi, npw, npwx, nbnd, how_many, &
                psi(:,nbase+1), hpsi(:,nbase+1), spsi(:,nbase+1), ethr, ew(1), nhpsi)

     CALL start_clock( 'paro:mp_bar' ); 
     CALL mp_barrier(inter_bgrp_comm)
     CALL stop_clock( 'paro:mp_bar' ); 
     CALL start_clock( 'paro:mp_sum' ); 
     psi (:,nbase+ibnd_start:nbase+ibnd_end) = psi (:,nbase+1:nbase+how_many) 
     CALL mp_allgather(psi (:,nbase+1:ndiag), column_type, recv_counts, displs, inter_bgrp_comm)
     hpsi(:,nbase+ibnd_start:nbase+ibnd_end) = hpsi(:,nbase+1:nbase+how_many) 
     CALL mp_allgather(hpsi(:,nbase+1:ndiag), column_type, recv_counts, displs, inter_bgrp_comm)
     spsi(:,nbase+ibnd_start:nbase+ibnd_end) = spsi(:,nbase+1:nbase+how_many) 
     CALL mp_allgather(spsi(:,nbase+1:ndiag), column_type, recv_counts, displs, inter_bgrp_comm)
     CALL stop_clock( 'paro:mp_sum' ); 

#if defined(__MPI)
     IF ( nproc_ortho == 1 ) THEN
#endif
        CALL rotate_HSpsi_gamma (  npwx, npw, ndiag, ndiag, psi, hpsi, overlap, spsi, ew )
#if defined(__MPI)
     ELSE
        CALL protate_HSpsi_gamma(  npwx, npw, ndiag, ndiag, psi, hpsi, overlap, spsi, ew )
     ENDIF
#endif

     !write (6,*) ' ew : ', ew(1:nbnd)
     ! only the first nbnd eigenvalues are relevant for convergence
     ! but only those that have actually been corrected should be trusted
     conv(1:nbnd) = .FALSE.
     conv(1:ntrust) = ABS(ew(1:ntrust)-eig(1:ntrust)).LT.ethr ; nconv = COUNT(conv(1:ntrust)) ; notconv = nbnd - nconv
     eig(1:nbnd)  = ew(1:nbnd)
     IF ( nconv == nbnd ) EXIT ParO_loop

  END DO ParO_loop

  evc(:,1:nbnd) = psi(:,1:nbnd)

  CALL mp_sum(nhpsi,inter_bgrp_comm)

  DEALLOCATE ( ew, conv, psi, hpsi, spsi )
  CALL mp_type_free( column_type )

  CALL stop_clock( 'paro_gamma' ); !write (6,*) ' exit paro diag'

END SUBROUTINE paro_gamma_new
