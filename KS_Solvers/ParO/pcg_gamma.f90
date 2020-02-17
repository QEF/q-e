! Copyright (C) 2015-2016 Aihui Zhou's group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
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
! The following file is for solving step 2 of the parallel orbital updating algorithm.
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE pcg_gamma( hs_1psi, g_1psi, psi0, spsi0, npw, npwx, nbnd, psi, ethr, iter, e, nhpsi )
  !----------------------------------------------------------------------------
  !
  !  ... solve the linear system
  !
  !      [ H - e S + lambda Pv ]|\tilde\psi> = [e S - H ] |psi>
  !      Pc [ H - e S ]|\tilde\psi> = Pc [ e S - H ] |psi>
  !
  ! the solution is sought until the residual norm is a fixed fraction of the RHS norm
  ! in this way the more accurate is the original problem the more accuratly the correction is computed
  !
  USE util_param,    ONLY : DP, stdout
  USE mp_bands_util, ONLY : intra_bgrp_comm, gstart
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! Following varibales are temporary

  COMPLEX(DP),INTENT(IN) :: psi0(npwx,nbnd)  ! psi0  needed to compute the Pv projection
  COMPLEX(DP),INTENT(IN) :: spsi0(npwx,nbnd) ! Spsi0  needed to compute the Pv projection
  INTEGER,  INTENT(IN)   :: npw, npwx, nbnd, iter ! input dimensions and iteration count
  REAL(DP), INTENT(IN)   :: ethr             ! threshold for convergence.
  REAL(DP), INTENT(IN)   :: e                ! current estimate of the target eigenvalue
  COMPLEX(DP),INTENT(INOUT) :: psi(npwx)     ! input: the current estimate of the eigenvector,
                                             ! output: the estimated correction vector
  INTEGER, INTENT(INOUT) :: nhpsi            ! (updated) number of Hpsi evaluations
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 5 ! maximum number of CG iterations 
  !
  COMPLEX(DP), ALLOCATABLE ::  hpsi(:), & ! the product of H and psi
                               spsi(:), & ! the product of S and psi
                               b(:),    & ! RHS for testing
                               r(:), p(:), sp(:), w(:),z(:) ! additional working vetors

  REAL(DP), ALLOCATABLE    ::  spsi0vec (:) ! the product of spsi0 and a trial vector

  REAL(DP) :: g0, g1, g2, beta, alpha, gamma, ethr_cg, ff, ff0
  INTEGER  :: npw2, npwx2, cg_iter, ibnd
  !
  REAL(DP), EXTERNAL :: DDOT
  COMPLEX(DP), EXTERNAL :: ZDOTC

  EXTERNAL  hs_1psi, g_1psi
  ! hs_1psi( npwx, npw, psi, hpsi, spsi )
  !
  CALL start_clock( 'pcg' )
  !write (6,*) ' enter pcg' , e
  !
  npw2  = 2*npw
  npwx2 = 2*npwx
  !
  ALLOCATE( hpsi( npwx ),  spsi( npwx ), r( npwx ), z( npwx ), b( npwx ) )
  ALLOCATE( spsi0vec(nbnd) )
  !
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  IF ( gstart == 2 ) psi(1) = CMPLX( DBLE( psi(1) ), 0.D0 ,kind=DP)
  !
  CALL start_clock( 'pcg:hs_1psi' )
  CALL hs_1psi( npwx, npw, psi, hpsi, spsi ) ! apply H and S to a single wavefunction (no bgrp parallelization inside)
  CALL stop_clock( 'pcg:hs_1psi' )
  !
  ! define CG algorithm RHS and initial solution
  r(:) = e * spsi(:) - hpsi(:)               ! initial gradient
  z(:) = r(:) ; call g_1psi(npwx,npw,z,e)    ! initial preconditioned gradient
!- project on conduction bands
  CALL start_clock( 'pcg:ortho' )
  CALL DGEMV( 'T', npw2, nbnd, 2.0_DP, spsi0, npwx2,       z, 1, 0.0_DP, spsi0vec, 1 )
  IF ( gstart == 2 ) spsi0vec(:) = spsi0vec(:) - CONJG(spsi0(1,:))*z(1)
  CALL mp_sum( spsi0vec, intra_bgrp_comm )
  CALL DGEMV( 'N', npw2, nbnd, -1.D0, psi0, npwx2, spsi0vec, 1, 1.0_DP,       z, 1 )
  CALL stop_clock( 'pcg:ortho' )
!-

  g0 = 2.d0 * DDOT( npw2, z ,1 ,r ,1) ; IF ( gstart == 2 ) g0 = g0 - CONJG (z(1)) * r(1)
  CALL mp_sum( g0, intra_bgrp_comm )         ! g0 = < initial z | initial r >
  ff = 0.d0 ; ff0 = ff
  !write (6,*) 0, g0, ff

  ALLOCATE( p( npwx ), sp( npwx ), w( npwx ) )

  ! ethr_cg = ethr  ! CG convergence threshold could be set from input but it is not ...
  ethr_cg = 1.0D-2  ! it makes more sense to fix the convergence of the CG solution to a 
                    ! fixed function of the RHS (see ethr_cg update later).
  ethr_cg = max ( 0.01*ethr, ethr_cg * g0 ) ! here we set the convergence of the correction

! save RHS for later
  b(:) = r(:)
  ! zero the trial solution: comment next line is you are looking for |\psi_new> = |\psi> + |\tilde \psi>
  psi(:) = ZERO
  ! initial search direction
  p(:) = z(:)

  iterate: &
  DO cg_iter = 1, maxter

     CALL start_clock( 'pcg:hs_1psi' )
     CALL hs_1psi( npwx, npw, p, w, sp ) ! apply H to a single wavefunction (no bgrp parallelization here!)
     CALL stop_clock( 'pcg:hs_1psi' )
     w = w - e* sp 

     gamma = 2.d0 * ddot( npw2, p ,1 ,w ,1) ; IF ( gstart == 2 ) gamma = gamma - CONJG(p(1)) * w(1)
     CALL mp_sum( gamma, intra_bgrp_comm )
     alpha = g0/gamma

     psi(:) = psi(:) + alpha * p(:)          ! updated solution
     r(:) = r(:) - alpha * w(:)              ! updated gradient
     g2 = 2.d0 * DDOT( npw2, z ,1 ,r ,1) ; IF ( gstart == 2 ) g2 = g2 - CONJG (z(1)) * r(1)
     CALL mp_sum( g2, intra_bgrp_comm )      ! g2 = < old z | new r >
     z(:) = r(:) ; call g_1psi(npwx,npw,z,e) ! updated preconditioned gradient
!- project on conduction bands
  CALL start_clock( 'pcg:ortho' )
  CALL DGEMV( 'T', npw2, nbnd, 2.0_DP, spsi0, npwx2,       z, 1, 0.0_DP, spsi0vec, 1 )
  IF ( gstart == 2 ) spsi0vec(:) = spsi0vec(:) - CONJG(spsi0(1,:))*z(1)
  CALL mp_sum( spsi0vec, intra_bgrp_comm )
  CALL DGEMV( 'N', npw2, nbnd, -1.D0, psi0, npwx2, spsi0vec, 1, 1.0_DP,       z, 1 )
  CALL stop_clock( 'pcg:ortho' )
!-
     g1 =  2.d0 * ddot( npw2, z, 1, r ,1) ; IF ( gstart == 2 ) g1 = g1 - CONJG(z(1)) * r(1)
     CALL mp_sum( g1, intra_bgrp_comm )      ! g1 = < new z | new r >
! evaluate the function
     ff = - ddot( npw2, psi, 1, r ,1) - ddot( npw2, psi, 1, b ,1) ; IF ( gstart == 2 ) ff = ff + 0.5d0* CONJG(psi(1)) * (r(1)+b(1))
     CALL mp_sum( ff, intra_bgrp_comm )
     !write (6,*) cg_iter, g1, ff,  gamma
     if ( ff > ff0 .AND. ff0 < 0.d0 ) psi(:) = psi(:) - alpha * p(:) ! fallback solution if last iteration failed to improve the function... exit and hope next time it'll be better
     IF ( ABS ( g1 ) < ethr_cg .OR. ( ff > ff0 ) ) EXIT iterate

     beta = (g1-g2)/g0                       ! Polak - Ribiere style update
     g0 = g1                                 ! < new z | new r >  ->  < old z | old r >
     p(:) = z(:) + beta * p(:)               ! updated search direction

     ff0 = ff                                ! updated minimum value reached by the function

  END DO iterate
  !write (6,*) ' exit  pcg loop'
! orthogonalize to psi0 ... 
! actually we are not doing that.. it would require both psi0 AND spsi0 to be computed and would
! remove an occupied orb contribution which is taken care of by the following rotate_wfc routine anyway

  !if ( cg_iter == maxter.and. ABS(g1) > ethr_cg) &
  !                           write (*,*) 'CG not converged maxter exceeded', cg_iter, g1, g0, ethr_cg
  !IF ( ABS ( g1 ) < ethr_cg) write (*,*) 'CG correction converged ', cg_iter, g1, ethr_cg
  !IF ( ABS ( g1 ) > g0     ) write (*,*) 'CG not converged ', cg_iter, g1, g0, ethr_cg

  nhpsi = nhpsi + cg_iter + 1
  !
  DEALLOCATE( spsi0vec )
  DEALLOCATE( r, p, sp, w, z )
  DEALLOCATE( hpsi,  spsi )
  !
  CALL stop_clock( 'pcg' )
  !
  RETURN
  !
END SUBROUTINE pcg_gamma
