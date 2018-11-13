!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )

FUNCTION KSDdot( n, A, incx, B, incy) result( res )
  !
  USE cg_param,       ONLY : DP
#if defined(__CUDA)
  USE cudafor
  USE cublas
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: incx,incy,n
  !
#if defined(__CUDA)
  REAL(DP), DEVICE, INTENT(IN) :: A(n), B(n)
#else
  REAL(DP),         INTENT(IN) :: A(n), B(n)
  REAL(DP),         EXTERNAL   :: ddot
#endif
  !
  REAL(DP) :: res
  !
#if defined(__CUDA)
  res = cublasDdot( n, A, incx, B, incy )
#else
  res = ddot( n, A, incx, B, incy )
#endif
  !
  RETURN
  !
END FUNCTION KSDdot

! define __VERBOSE to print a message after each eigenvalue is computed
!
!----------------------------------------------------------------------------
SUBROUTINE ccgdiagg_gpu( hs_1psi_gpu, s_1psi_gpu, precondition_d, &
     npwx, npw, nbnd, npol, psi_d, e_d, btype, &
     ethr, maxter, reorder, notconv, avg_iter )
  !----------------------------------------------------------------------------
  !
  ! ... "poor man" iterative diagonalization of a complex hermitian matrix
  ! ... through preconditioned conjugate gradient algorithm
  ! ... Band-by-band algorithm with minimal use of memory
  ! ... Calls hs_1psi and s_1psi to calculate H|psi> + S|psi> and S|psi>
  ! ... Works for generalized eigenvalue problem (US pseudopotentials) as well
  !
#if defined(__CUDA)
  USE cudafor
  USE cublas
#endif
  USE cg_param,       ONLY : DP
  USE mp_bands_util,  ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp,             ONLY : mp_sum
#if defined(__VERBOSE)
  USE cg_param,     ONLY : stdout
#endif
  !
  IMPLICIT NONE
  !
  ! ... Mathematical constants
  ! 
  REAL(DP), PARAMETER :: pi     = 3.14159265358979323846_DP
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npwx, npw, nbnd, npol, maxter
  INTEGER,     INTENT(IN)    :: btype(nbnd)
  REAL(DP),    INTENT(IN)    :: precondition_d(npwx*npol), ethr
  COMPLEX(DP), INTENT(INOUT) :: psi_d(npwx*npol,nbnd)
  REAL(DP),    INTENT(INOUT) :: e_d(nbnd)
  INTEGER,     INTENT(OUT)   :: notconv
  REAL(DP),    INTENT(OUT)   :: avg_iter
#if defined(__CUDA)
  attributes(DEVICE) :: precondition_d, psi_d, e_d
#endif
  !
  ! ... local variables
  !
  INTEGER                  :: i, j, k, m, m_start, m_end, iter, moved
  COMPLEX(DP), ALLOCATABLE :: hpsi_d(:), spsi_d(:), g_d(:), cg_d(:)
  COMPLEX(DP), ALLOCATABLE :: scg_d(:), ppsi_d(:), g0_d(:), lagrange_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: hpsi_d, spsi_d, g_d, cg_d, scg_d, ppsi_d, g0_d, lagrange_d
#endif
  COMPLEX(DP), ALLOCATABLE :: lagrange(:)
  REAL(DP)                 :: gamma, ddot_temp, es_1, es(2)
  REAL(DP), ALLOCATABLE    :: e(:)
  REAL(DP)                 :: a0, b0, gg0, gg, gg1, cg0, e0, psi_norm, sint, cost
  REAL(DP)                 :: theta, cos2t, sin2t
  LOGICAL                  :: reorder
  INTEGER                  :: kdim, kdmx, kdim2, ierr, istat
  REAL(DP)                 :: empty_ethr, ethr_m
  !
  ! ... external functions
  !
  REAL (DP), EXTERNAL :: ksDdot
  EXTERNAL  hs_1psi_gpu, s_1psi_gpu
  ! hs_1psi( npwx, npw, psi, hpsi, spsi )
  ! s_1psi( npwx, npw, psi, spsi )
  !
  CALL start_clock( 'ccgdiagg' )
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx * npol
     kdmx = npwx * npol
     !
  END IF
  !
  kdim2 = 2 * kdim
  !
  ALLOCATE(  hpsi_d(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate hpsi_d ', ABS(ierr) )
  ALLOCATE(  spsi_d(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate spsi_d ', ABS(ierr) )
  ALLOCATE(  g_d(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate g_d ', ABS(ierr) )
  ALLOCATE(  cg_d(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate cg_d ', ABS(ierr) )
  ALLOCATE(  scg_d(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate scg_d ', ABS(ierr) )
  ALLOCATE(  ppsi_d(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate ppsi_d ', ABS(ierr) )
  ALLOCATE(  g0_d(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate g0_d ', ABS(ierr) )
  !
  ALLOCATE(  lagrange_d(nbnd), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate lagrange_d ', ABS(ierr) )
  !
  ALLOCATE(  lagrange(nbnd), STAT=ierr  )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate lagrange ', ABS(ierr) )
  ALLOCATE(  e(nbnd) )
  !
  avg_iter = 0.D0
  notconv  = 0
  moved    = 0
  !
  ! ... every eigenfunction is calculated separately
  !
  DO m = 1, nbnd
     !
     IF ( btype(m) == 1 ) THEN
        !
        ethr_m = ethr
        !
     ELSE
        !
        ethr_m = empty_ethr
        !
     END IF
     !
     spsi_d     = ZERO
     scg_d      = ZERO
     hpsi_d     = ZERO
     g_d        = ZERO
     cg_d       = ZERO
     g0_d       = ZERO
     ppsi_d     = ZERO
     lagrange_d = ZERO
     !
     ! ... calculate S|psi>
     !
     CALL s_1psi_gpu( npwx, npw, psi_d(1,m), spsi_d )
     !
     ! ... orthogonalize starting eigenfunction to those already calculated
     !
     call divide(inter_bgrp_comm,m,m_start,m_end); !write(*,*) m,m_start,m_end
     lagrange = ZERO
     if(m_start.le.m_end) &
          CALL ZGEMV( 'C', kdim, m_end-m_start+1, ONE, psi_d(1,m_start), &
                      kdmx, spsi_d, 1, ZERO, lagrange_d(m_start), 1 )
     if(m_start.le.m_end) lagrange(m_start:m_end) = lagrange_d(m_start:m_end)

     CALL mp_sum( lagrange( 1:m ), inter_bgrp_comm )
     !
     CALL mp_sum( lagrange( 1:m ), intra_bgrp_comm )
     !
     psi_norm = DBLE( lagrange(m) )
     lagrange_d(1:m) = lagrange(1:m)
     !
     DO j = 1, m - 1
        !$cuf kernel do(1) <<<*,*>>>
        DO i = 1, kdmx
           !
           psi_d(i,m)  = psi_d(i,m) - lagrange_d(j) * psi_d(i,j)
           !
        END DO
        !
        psi_norm = psi_norm - ( DBLE( lagrange(j) )**2 + AIMAG( lagrange(j) )**2 )
        !
     END DO
     !
     psi_norm = SQRT( psi_norm )
     !
     !$cuf kernel do(1) <<<*,*>>>
     DO i = 1, kdmx
        psi_d(i,m) = psi_d(i,m) / psi_norm
     END DO
     !
     ! ... calculate starting gradient (|hpsi> = H|psi>) ...
     !
     CALL hs_1psi_gpu( npwx, npw, psi_d(1,m), hpsi_d, spsi_d )
     !
     ! ... and starting eigenvalue (e = <y|PHP|y> = <psi|H|psi>)
     !
     ! ... NB:  ddot(2*npw,a,1,b,1) = REAL( zdotc(npw,a,1,b,1) )
     !
     e(m) = ksDdot( kdim2, psi_d(1,m), 1, hpsi_d, 1 )
     !
     CALL mp_sum( e(m), intra_bgrp_comm )
     !
     ! ... start iteration for this band
     !
     iterate: DO iter = 1, maxter
        !
        ! ... calculate  P (PHP)|y>
        ! ... ( P = preconditioning matrix, assumed diagonal )
        !
        !$cuf kernel do(1) <<<*,*>>>
        DO i = 1, kdmx
           g_d(i)    = hpsi_d(i) / precondition_d(i)
           ppsi_d(i) = spsi_d(i) / precondition_d(i)
        END DO
        !
        ! ... ppsi is now S P(P^2)|y> = S P^2|psi>)
        !
        es(1) = ksDdot( kdim2, spsi_d(1), 1, g_d(1), 1 )
        es(2) = ksDdot( kdim2, spsi_d(1), 1, ppsi_d(1), 1 )
        !
        CALL mp_sum( es , intra_bgrp_comm )
        !
        es(1) = es(1) / es(2)
        es_1 = es(1)
        !
        !$cuf kernel do(1) <<<*,*>>>
        DO i = 1, kdmx
           g_d(i) = g_d(i) - es_1 * ppsi_d(i)
        END DO
        !
        ! ... e1 = <y| S P^2 PHP|y> / <y| S S P^2|y>  ensures that 
        ! ... <g| S P^2|y> = 0
        ! ... orthogonalize to lowest eigenfunctions (already calculated)
        !
        ! ... scg is used as workspace
        !
        CALL s_1psi_gpu( npwx, npw, g_d(1), scg_d(1) )
        !
        lagrange(1:m-1) = ZERO
        call divide(inter_bgrp_comm,m-1,m_start,m_end); !write(*,*) m-1,m_start,m_end
        if(m_start.le.m_end) &
        CALL ZGEMV( 'C', kdim, m_end-m_start+1, ONE, psi_d(1,m_start), &
                    kdmx, scg_d, 1, ZERO, lagrange_d(m_start), 1 )
        if(m_start.le.m_end) lagrange(m_start:m_end) = lagrange_d(m_start:m_end)
        CALL mp_sum( lagrange( 1:m-1 ), inter_bgrp_comm )
        !
        CALL mp_sum( lagrange( 1:m-1 ), intra_bgrp_comm )
        !
        lagrange_d(1:m) = lagrange(1:m)
        !
        DO j = 1, ( m - 1 )
           !
           !$cuf kernel do(1) <<<*,*>>>
           DO i = 1, kdmx
              g_d(i)   = g_d(i)   - lagrange_d(j) * psi_d(i,j)
              scg_d(i) = scg_d(i) - lagrange_d(j) * psi_d(i,j)
           END DO
           !
        END DO
        !
        IF ( iter /= 1 ) THEN
           !
           ! ... gg1 is <g(n+1)|S|g(n)> (used in Polak-Ribiere formula)
           !
           gg1 = ksDdot( kdim2, g_d(1), 1, g0_d(1), 1 )
           !
           CALL mp_sum( gg1, intra_bgrp_comm )
           !
        END IF
        !
        ! ... gg is <g(n+1)|S|g(n+1)>
        !
        !$cuf kernel do(1) <<<*,*>>>
        DO i = 1, kdmx
           g0_d(i) = scg_d(i) * precondition_d(i)
        END DO
        !
        gg = ksDdot( kdim2, g_d(1), 1, g0_d(1), 1 )
        !
        CALL mp_sum( gg, intra_bgrp_comm )
        !
        IF ( iter == 1 ) THEN
           !
           ! ... starting iteration, the conjugate gradient |cg> = |g>
           !
           gg0 = gg
           !
           !$cuf kernel do(1) <<<*,*>>>
           DO i = 1, kdmx
              cg_d(i) = g_d(i)
           END DO
           !
        ELSE
           !
           ! ... |cg(n+1)> = |g(n+1)> + gamma(n) * |cg(n)>
           !
           ! ... Polak-Ribiere formula :
           !
           gamma = ( gg - gg1 ) / gg0
           gg0   = gg
           !
           ! 
           !
           ! See comment below
           !!DO i = 1, kdmx
           !!   cg_d(i) = g_d(i) + cg_d(i) * gamma
           !!END DO
           !
           ! ... The following is needed because <y(n+1)| S P^2 |cg(n+1)> 
           ! ... is not 0. In fact :
           ! ... <y(n+1)| S P^2 |cg(n)> = sin(theta)*<cg(n)|S|cg(n)>
           !
           psi_norm = gamma * cg0 * sint
           !
           !$cuf kernel do(1) <<<*,*>>>
           DO i = 1, kdmx
              !          v== this breaks the logic, done here for performance
              cg_d(i) = (g_d(i) + cg_d(i) * gamma) - psi_norm * psi_d(i,m)
           END DO
           !
        END IF
        !
        ! ... |cg> contains now the conjugate gradient
        !
        ! ... |scg> is S|cg>
        !
        CALL hs_1psi_gpu( npwx, npw, cg_d(1), ppsi_d(1), scg_d(1) )
        !
        cg0 = ksDdot( kdim2, cg_d(1), 1, scg_d(1), 1 )
        !
        CALL mp_sum(  cg0 , intra_bgrp_comm )
        !
        cg0 = SQRT( cg0 )
        !
        ! ... |ppsi> contains now HP|cg>
        ! ... minimize <y(t)|PHP|y(t)> , where :
        ! ...                         |y(t)> = cos(t)|y> + sin(t)/cg0 |cg>
        ! ... Note that  <y|P^2S|y> = 1, <y|P^2S|cg> = 0 ,
        ! ...           <cg|P^2S|cg> = cg0^2
        ! ... so that the result is correctly normalized :
        ! ...                           <y(t)|P^2S|y(t)> = 1
        !
        a0 = 2.D0 *  ksDdot( kdim2, psi_d(1,m), 1, ppsi_d(1), 1 ) / cg0
        !
        CALL mp_sum(  a0 , intra_bgrp_comm )
        !
        b0 = ksDdot( kdim2, cg_d(1), 1, ppsi_d(1), 1 ) / cg0**2
        !
        CALL mp_sum(  b0 , intra_bgrp_comm )
        !
        e0 = e(m)
        !
        theta = 0.5D0 * ATAN( a0 / ( e0 - b0 ) )
        !
        cost = COS( theta )
        sint = SIN( theta )
        !
        cos2t = cost*cost - sint*sint
        sin2t = 2.D0*cost*sint
        !
        es(1) = 0.5D0 * (   ( e0 - b0 ) * cos2t + a0 * sin2t + e0 + b0 )
        es(2) = 0.5D0 * ( - ( e0 - b0 ) * cos2t - a0 * sin2t + e0 + b0 )
        !
        ! ... there are two possible solutions, choose the minimum
        !
        IF ( es(2) < es(1) ) THEN
           !
           theta = theta + 0.5D0 * pi
           !
           cost = COS( theta )
           sint = SIN( theta )
           !
        END IF
        !
        ! ... new estimate of the eigenvalue
        !
        e(m) = MIN( es(1), es(2) )
        ! ... upgrade |psi>
        !
        !$cuf kernel do(1) <<<*,*>>>
        DO i = 1, kdmx
           psi_d(i,m) = cost * psi_d(i,m) + sint / cg0 * cg_d(i)
        END DO
        !
        ! ... here one could test convergence on the energy
        !
        IF ( ABS( e(m) - e0 ) < ethr_m ) EXIT iterate
        !
        ! ... upgrade H|psi> and S|psi>
        !
        !$cuf kernel do(1) <<<*,*>>>
        DO i = 1, kdmx
           spsi_d(i) = cost * spsi_d(i) + sint / cg0 * scg_d(i)
           !
           hpsi_d(i) = cost * hpsi_d(i) + sint / cg0 * ppsi_d(i)
        END DO
        !
     END DO iterate
     !
#if defined(__VERBOSE)
     IF ( iter >= maxter ) THEN
        WRITE(stdout,'("e(",i4,") = ",f12.6," eV  (not converged after ",i3,&
             & " iterations)")') m, e(m)*13.6058, iter
     ELSE
        WRITE(stdout,'("e(",i4,") = ",f12.6," eV  (",i3," iterations)")') &
             m, e(m)*13.6058, iter
     END IF
     FLUSH (stdout)
#endif
     IF ( iter >= maxter ) notconv = notconv + 1
     !
     avg_iter = avg_iter + iter + 1
     ! ... reorder eigenvalues if they are not in the right order
     ! ... ( this CAN and WILL happen in not-so-special cases )
     !
     !
     IF ( m > 1 .AND. reorder ) THEN
        !
        IF ( e(m) - e(m-1) < - 2.D0 * ethr_m ) THEN
           ! ... if the last calculated eigenvalue is not the largest...
           !
           DO i = m - 2, 1, - 1
              !
              IF ( e(m) - e(i) > 2.D0 * ethr_m ) EXIT
              !
           END DO
           !
           i = i + 1
           !
           moved = moved + 1
           !
           ! ... last calculated eigenvalue should be in the 
           ! ... i-th position: reorder
           !
           e0 = e(m)
           !
           !$cuf kernel do(1) <<<*,*>>>
           DO k = 1, kdmx
              ppsi_d(k) = psi_d(k,m)
           END DO
           !
           DO j = m, i + 1, - 1
              !
              e(j) = e(j-1)
              !
              !$cuf kernel do(1) <<<*,*>>>
              DO k = 1, kdmx
                 psi_d(k,j) = psi_d(k,j-1)
              END DO
              !
           END DO
           !
           e(i) = e0
           !
           !$cuf kernel do(1) <<<*,*>>>
           DO k = 1, kdmx
              psi_d(k,i) = ppsi_d(k)
           END DO
           !
           ! ... this procedure should be good if only a few inversions occur,
           ! ... extremely inefficient if eigenvectors are often in bad order
           ! ... ( but this should not happen )
           !
        END IF
        !
     END IF
     !
  END DO
  !
  avg_iter = avg_iter / DBLE( nbnd )
  !
  ! STORING e in e_d since eigenvalues are always on the host
  e_d = e
  !
  DEALLOCATE( lagrange )
  DEALLOCATE( e )
  DEALLOCATE( lagrange_d )
  DEALLOCATE( ppsi_d )
  DEALLOCATE( g0_d )
  DEALLOCATE( cg_d )
  DEALLOCATE( g_d )
  DEALLOCATE( hpsi_d )
  DEALLOCATE( scg_d )
  DEALLOCATE( spsi_d )
  !
  CALL stop_clock( 'ccgdiagg' )
  !
  RETURN
  !
END SUBROUTINE ccgdiagg_gpu
