!
! Copyright (C) 2002-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! define __VERBOSE to print a message after each eigenvalue is computed
!----------------------------------------------------------------------------
SUBROUTINE rcgdiagg( npwx, npw, nbnd, psi, e, btype, precondition, &
                     ethr, maxter, reorder, notconv, avg_iter )
  !----------------------------------------------------------------------------
  !
  ! ... "poor man" iterative diagonalization of a complex hermitian matrix
  ! ... through preconditioned conjugate gradient algorithm
  ! ... Band-by-band algorithm with minimal use of memory
  ! ... Calls h_1psi and s_1psi to calculate H|psi> and S|psi>
  ! ... Works for generalized eigenvalue problem (US pseudopotentials) as well
  !
  USE constants, ONLY : pi
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : gstart
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
#if defined(__VERBOSE)
  USE io_global, only : stdout
#endif
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,      INTENT(IN)    :: npwx, npw, nbnd, maxter
  INTEGER,      INTENT(IN)    :: btype(nbnd)
  REAL (DP),    INTENT(IN)    :: precondition(npw), ethr
  COMPLEX (DP), INTENT(INOUT) :: psi(npwx,nbnd)
  REAL (DP),    INTENT(INOUT) :: e(nbnd)
  INTEGER,      INTENT(OUT)   :: notconv
  REAL (DP),    INTENT(OUT)   :: avg_iter
  !
  ! ... local variables
  !
  INTEGER                   :: i, j, m, iter, moved
  REAL (DP),    ALLOCATABLE :: lagrange(:)
  COMPLEX (DP), ALLOCATABLE :: hpsi(:), spsi(:), g(:), cg(:), &
                               scg(:), ppsi(:), g0(:)  
  REAL (DP)                 :: psi_norm, a0, b0, gg0, gamma, gg, gg1, &
                               cg0, e0, es(2)
  REAL (DP)                 :: theta, cost, sint, cos2t, sin2t
  LOGICAL                   :: reorder
  INTEGER                   :: npw2, npwx2
  REAL (DP)                 :: empty_ethr
  !
  ! ... external functions
  !
  REAL (DP), EXTERNAL :: ddot
  !
  !
  CALL start_clock( 'rcgdiagg' )
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  npw2 = 2 * npw
  npwx2 = 2 * npwx
  !
  ALLOCATE( spsi( npwx ) )
  ALLOCATE( scg(  npwx ) )
  ALLOCATE( hpsi( npwx ) )
  ALLOCATE( g(    npwx ) )
  ALLOCATE( cg(   npwx ) )
  ALLOCATE( g0(   npwx ) )
  ALLOCATE( ppsi( npwx ) )
  !    
  ALLOCATE( lagrange( nbnd ) )
  !
  avg_iter = 0.D0
  notconv  = 0
  moved    = 0
  !
  ! ... every eigenfunction is calculated separately
  !
  DO m = 1, nbnd
     !
     ! ... calculate S|psi>
     !
     CALL s_1psi( npwx, npw, psi(1,m), spsi )
     !
     ! ... orthogonalize starting eigenfunction to those already calculated
     !
     CALL DGEMV( 'T', npw2, m, 2.D0, psi, npwx2, spsi, 1, 0.D0, lagrange, 1 )
     !
     IF ( gstart == 2 ) lagrange(1:m) = lagrange(1:m) - psi(1,1:m) * spsi(1)
     !
     CALL mp_sum( lagrange( 1:m ), intra_bgrp_comm )
     !
     psi_norm = lagrange(m)
     !
     DO j = 1, m - 1
        !
        psi(:,m)  = psi(:,m) - lagrange(j) * psi(:,j)
        !
        psi_norm = psi_norm - lagrange(j)**2
        !
     END DO
     !
     psi_norm = SQRT( psi_norm )
     !
     psi(:,m) = psi(:,m) / psi_norm
     ! ... set Im[ psi(G=0) ] -  needed for numerical stability
     IF ( gstart == 2 ) psi(1,m) = CMPLX( DBLE(psi(1,m)), 0.D0 ,kind=DP)
     !
     ! ... calculate starting gradient (|hpsi> = H|psi>) ...
     !
     CALL h_1psi( npwx, npw, psi(1,m), hpsi, spsi )
     !
     ! ... and starting eigenvalue (e = <y|PHP|y> = <psi|H|psi>)
     !
     ! ... NB:  ddot(2*npw,a,1,b,1) = DBLE( zdotc(npw,a,1,b,1) )
     !
     e(m) = 2.D0 * ddot( npw2, psi(1,m), 1, hpsi, 1 )
     !
     IF ( gstart == 2 ) e(m) = e(m) - psi(1,m) * hpsi(1)
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
        g(1:npw)    = hpsi(1:npw) / precondition(:)
        ppsi(1:npw) = spsi(1:npw) / precondition(:)
        !
        ! ... ppsi is now S P(P^2)|y> = S P^2|psi>)
        !
        es(1) = 2.D0 * ddot( npw2, spsi(1), 1, g(1), 1 )
        es(2) = 2.D0 * ddot( npw2, spsi(1), 1, ppsi(1), 1 )
        !
        IF ( gstart == 2 ) THEN
           !
           es(1) = es(1) - spsi(1) * g(1)
           es(2) = es(2) - spsi(1) * ppsi(1)
           !
        END IF
        !
        CALL mp_sum(  es , intra_bgrp_comm )
        !
        es(1) = es(1) / es(2)
        !
        g(:) = g(:) - es(1) * ppsi(:)
        !
        ! ... e1 = <y| S P^2 PHP|y> / <y| S S P^2|y>  ensures that  
        ! ... <g| S P^2|y> = 0
        !
        ! ... orthogonalize to lowest eigenfunctions (already calculated)
        !
        ! ... scg is used as workspace
        !
        CALL s_1psi( npwx, npw, g(1), scg(1) )
        !
        CALL DGEMV( 'T', npw2, ( m - 1 ), 2.D0, &
                    psi, npwx2, scg, 1, 0.D0, lagrange, 1 )
        !
        IF ( gstart == 2 ) &
           lagrange(1:m-1) = lagrange(1:m-1) - psi(1,1:m-1) * scg(1)
        !
        CALL mp_sum( lagrange( 1 : m-1 ), intra_bgrp_comm )
        !
        DO j = 1, ( m - 1 )
           !
           g(:)   = g(:)   - lagrange(j) * psi(:,j)
           scg(:) = scg(:) - lagrange(j) * psi(:,j)
           !
        END DO
        !
        IF ( iter /= 1 ) THEN
           !
           ! ... gg1 is <g(n+1)|S|g(n)> (used in Polak-Ribiere formula)
           !
           gg1 = 2.D0 * ddot( npw2, g(1), 1, g0(1), 1 )
           !
           IF ( gstart == 2 ) gg1 = gg1 - g(1) * g0(1)
           !
           CALL mp_sum(  gg1 , intra_bgrp_comm )
           !
        END IF
        !
        ! ... gg is <g(n+1)|S|g(n+1)>
        !
        g0(:) = scg(:)
        !
        g0(1:npw) = g0(1:npw) * precondition(:)
        !
        gg = 2.D0 * ddot( npw2, g(1), 1, g0(1), 1 )
        !
        IF ( gstart == 2 ) gg = gg - g(1) * g0(1)
        !
        CALL mp_sum(  gg , intra_bgrp_comm )
        !
        IF ( iter == 1 ) THEN
           !
           ! ... starting iteration, the conjugate gradient |cg> = |g>
           !
           gg0 = gg
           !
           cg(:) = g(:)
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
           cg(:) = cg(:) * gamma
           cg(:) = g + cg(:)
           !
           ! ... The following is needed because <y(n+1)| S P^2 |cg(n+1)> 
           ! ... is not 0. In fact :
           ! ... <y(n+1)| S P^2 |cg(n)> = sin(theta)*<cg(n)|S|cg(n)>
           !
           psi_norm = gamma * cg0 * sint
           !
           cg(:) = cg(:) - psi_norm * psi(:,m)
           !
        END IF
        !
        ! ... |cg> contains now the conjugate gradient
        ! ... set Im[ cg(G=0) ] -  needed for numerical stability
        IF ( gstart == 2 ) cg(1) = CMPLX( DBLE(cg(1)), 0.D0 ,kind=DP)
        !
        ! ... |scg> is S|cg>
        !
        CALL h_1psi( npwx, npw, cg(1), ppsi(1), scg(1) )
        !
        cg0 = 2.D0 * ddot( npw2, cg(1), 1, scg(1), 1 )
        !
        IF ( gstart == 2 ) cg0 = cg0 - cg(1) * scg(1)
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
        a0 = 4.D0 * ddot( npw2, psi(1,m), 1, ppsi(1), 1 )
        !
        IF ( gstart == 2 ) a0 = a0 - 2.D0 * psi(1,m) * ppsi(1)
        !
        a0 = a0 / cg0
        !
        CALL mp_sum(  a0 , intra_bgrp_comm )
        !
        b0 = 2.D0 * ddot( npw2, cg(1), 1, ppsi(1), 1 )
        !
        IF ( gstart == 2 ) b0 = b0 - cg(1) * ppsi(1)
        !
        b0 = b0 / cg0**2
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
        !
        ! ... upgrade |psi>
        !
        psi(:,m) = cost * psi(:,m) + sint / cg0 * cg(:)
        !
        ! ... here one could test convergence on the energy
        !
        IF ( btype(m) == 1 ) THEN
           !
           IF ( ABS( e(m) - e0 ) < ethr ) EXIT iterate
           !
        ELSE
           !
           IF ( ABS( e(m) - e0 ) < empty_ethr ) EXIT iterate
           !
        END IF
        !
        ! ... upgrade H|psi> and S|psi>
        !
        spsi(:) = cost * spsi(:) + sint / cg0 * scg(:)
        !
        hpsi(:) = cost * hpsi(:) + sint / cg0 * ppsi(:)
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
     !
     ! ... reorder eigenvalues if they are not in the right order
     ! ... ( this CAN and WILL happen in not-so-special cases )
     !
     IF ( m > 1 .AND. reorder ) THEN
        !
        IF ( e(m) - e(m-1) < - 2.D0 * ethr ) THEN
           !
           ! ... if the last calculated eigenvalue is not the largest...
           !
           DO i = m - 2, 1, - 1
              !
              IF ( e(m) - e(i) > 2.D0 * ethr ) EXIT
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
           ppsi(:) = psi(:,m)
           !
           DO j = m, i + 1, - 1
              !
              e(j) = e(j-1)
              !
              psi(:,j) = psi(:,j-1)
              !
           END DO
           !
           e(i) = e0
           !
           psi(:,i) = ppsi(:)
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
  DEALLOCATE( lagrange )
  DEALLOCATE( ppsi )
  DEALLOCATE( g0 )
  DEALLOCATE( cg )
  DEALLOCATE( g )
  DEALLOCATE( hpsi )
  DEALLOCATE( scg )
  DEALLOCATE( spsi )
  !
  CALL stop_clock( 'rcgdiagg' )
  !
  RETURN
  !
END SUBROUTINE rcgdiagg
