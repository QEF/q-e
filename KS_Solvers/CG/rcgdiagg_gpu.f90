!
! Copyright (C) 2002-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


SUBROUTINE cgcudaDGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#if defined(__CUDA)
    use cudafor
    use cublas
#endif
    IMPLICIT NONE
    DOUBLE PRECISION :: ALPHA,BETA
    INTEGER :: INCX,INCY,LDA,M,N
    CHARACTER :: TRANS
    DOUBLE PRECISION :: A(LDA,*),X(*),Y(*)
#if defined(__CUDA)
    attributes(device) :: A, X, Y
#endif
    !
    call DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    !
END SUBROUTINE cgcudaDGEMV

! define __VERBOSE to print a message after each eigenvalue is computed
!----------------------------------------------------------------------------
SUBROUTINE rcgdiagg_gpu( hs_1psi_gpu, s_1psi_gpu, precondition_d, &
                     npwx, npw, nbnd, psi_d, e_d, btype, &
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
  USE mp_bands_util,  ONLY : intra_bgrp_comm, inter_bgrp_comm, gstart
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
  INTEGER,     INTENT(IN)    :: npwx, npw, nbnd, maxter
  INTEGER,     INTENT(IN)    :: btype(nbnd)
  REAL(DP),    INTENT(IN)    :: precondition_d(npw), ethr
  COMPLEX(DP), INTENT(INOUT) :: psi_d(npwx,nbnd)
  REAL(DP),    INTENT(INOUT) :: e_d(nbnd)
  INTEGER,     INTENT(OUT)   :: notconv
  REAL(DP),    INTENT(OUT)   :: avg_iter
#if defined(__CUDA)
  attributes(DEVICE) :: precondition_d, psi_d, e_d
#endif
  !
  ! ... local variables
  !
  INTEGER                  :: i, j, l, m, m_start, m_end, iter, moved
  REAL(DP),    ALLOCATABLE :: lagrange_d(:)
  REAL(DP),    ALLOCATABLE :: lagrange(:), e(:)
  COMPLEX(DP), ALLOCATABLE :: hpsi_d(:), spsi_d(:), g_d(:), cg_d(:), &
                              scg_d(:), ppsi_d(:), g0_d(:), psi_aux(:)
  COMPLEX(DP)              :: psi1, hpsi1, spsi1, ppsi1, scg1, cg1, g1, g01
  REAL(DP)                 :: psi_norm, a0, b0, gg0, gamma, gg, gg1, &
                              cg0, e0, es(2), aux
  REAL(DP)                 :: es1
  REAL(DP)                 :: theta, cost, sint, cos2t, sin2t
  LOGICAL                  :: reorder
  INTEGER                  :: npw2, npwx2
  REAL(DP)                 :: empty_ethr, ethr_m
#if defined(__CUDA)
  attributes(DEVICE) :: lagrange_d, hpsi_d, spsi_d, g_d, cg_d, scg_d, ppsi_d, g0_d
#endif
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: cgddot
  EXTERNAL  hs_1psi_gpu,    s_1psi_gpu 
  ! hs_1psi( npwx, npw, psi, hpsi, spsi )
  ! s_1psi( npwx, npw, psi, spsi )
  !
  !
  CALL start_clock( 'rcgdiagg' )
  !
  IF ( gstart == -1 ) CALL errore( 'regter', 'gstart variable not initialized', 1 )
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  npw2 = 2 * npw
  npwx2 = 2 * npwx
  !
  ALLOCATE( spsi_d( npwx ) )
  ALLOCATE( scg_d(  npwx ) )
  ALLOCATE( hpsi_d( npwx ) )
  ALLOCATE( g_d(    npwx ) )
  ALLOCATE( cg_d(   npwx ) )
  ALLOCATE( g0_d(   npwx ) )
  ALLOCATE( ppsi_d( npwx ) )
  !
  ALLOCATE( lagrange_d( nbnd ) )
  ALLOCATE( lagrange  ( nbnd ) )
  ALLOCATE( e         ( nbnd ) )
  ALLOCATE( psi_aux   ( nbnd ) )
  !
  ! Sync eigenvalues that will remain on the Host
  e(1:nbnd) = e_d(1:nbnd)
  !print *, 'init ', e(1:nbnd)
  !
  avg_iter = 0.D0
  notconv  = 0
  moved    = 0
  !
  ! ... every eigenfunction is calculated separately
  !
  DO m = 1, nbnd

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
     ! ... calculate S|psi>
     !
     CALL s_1psi_gpu( npwx, npw, psi_d(1,m), spsi_d )
     !
     ! ... orthogonalize starting eigenfunction to those already calculated
     !
     call divide(inter_bgrp_comm,m,m_start,m_end); !write(*,*) m,m_start,m_end
     lagrange = 0.d0
     if(m_start.le.m_end) &
     CALL cgcudaDGEMV( 'T', npw2, m_end-m_start+1, 2.D0, psi_d(1,m_start), npwx2, spsi_d, 1, 0.D0, lagrange_d(m_start), 1 )
     if(m_start.le.m_end) lagrange( m_start:m_end ) = lagrange_d( m_start:m_end )
     !print *, 'lagrange ', lagrange(1:m)

     CALL mp_sum( lagrange( 1:m ), inter_bgrp_comm )
     IF ( gstart == 2 ) THEN
        psi_aux(1:m) = psi_d(1,1:m)
        spsi1 = spsi_d(1)
        lagrange(1:m) = lagrange(1:m) - psi_aux(1:m) * spsi1
     END IF
     !
     CALL mp_sum( lagrange( 1:m ), intra_bgrp_comm )
     !
     psi_norm = lagrange(m)
     lagrange_d(1:m) = lagrange(1:m)
     !
     DO j = 1, m - 1
        !
        !$cuf kernel do(1) <<<*,*>>>
        DO i = 1, npwx
           psi_d(i,m)  = psi_d(i,m) - lagrange_d(j) * psi_d(i,j)
        END DO
        !
        !print *, 'psi_norm ', j, psi_norm
        psi_norm = psi_norm - lagrange(j)**2
        !
     END DO
     !
     psi_norm = SQRT( psi_norm )
     !print *, 'psi_norm 178', psi_norm
     !
     !$cuf kernel do(1) <<<*,*>>>
     DO i = 1, npwx
        psi_d(i,m) = psi_d(i,m) / psi_norm
        ! ... set Im[ psi(G=0) ] -  needed for numerical stability
        IF (i == 1) THEN
          IF ( gstart == 2 ) psi_d(1,m) = CMPLX( DBLE(psi_d(1,m)), 0.D0 ,kind=DP)
        END IF
     END DO
     !
     ! ... calculate starting gradient (|hpsi> = H|psi>) ...
     !
     CALL hs_1psi_gpu( npwx, npw, psi_d(1,m), hpsi_d, spsi_d )
     !
     ! ... and starting eigenvalue (e = <y|PHP|y> = <psi|H|psi>)
     !
     ! ... NB:  ddot(2*npw,a,1,b,1) = DBLE( zdotc(npw,a,1,b,1) )
     !
     e(m) = 2.D0 * cgddot( npw2, psi_d(1,m), 1, hpsi_d, 1 )
     !print *, 'e(m)', e(m)
     IF ( gstart == 2 ) THEN
        psi1  = psi_d(1,m)
        hpsi1 = hpsi_d(1)
        e(m) = e(m) - psi1 * hpsi1
     END IF
     !
     CALL mp_sum( e(m), intra_bgrp_comm )
     !print *, 'before iterate', psi1, hpsi1, spsi1, e(1:nbnd)
     !
     ! ... start iteration for this band
     !
     iterate: DO iter = 1, maxter
        !
        ! ... calculate  P (PHP)|y>
        ! ... ( P = preconditioning matrix, assumed diagonal )
        !
        !$cuf kernel do(1) <<<*,*>>>
        DO i = 1, npw
           g_d(i)  = hpsi_d(i) / precondition_d(i)
           ppsi_d(i) = spsi_d(i) / precondition_d(i)
        END DO
        !
        ! ... ppsi is now S P(P^2)|y> = S P^2|psi>)
        !
        es(1) = 2.D0 * cgddot( npw2, spsi_d(1), 1, g_d(1), 1 )
        es(2) = 2.D0 * cgddot( npw2, spsi_d(1), 1, ppsi_d(1), 1 )
        !
        IF ( gstart == 2 ) THEN
           !
           g1 = g_d(1); ppsi1 = ppsi_d(1); spsi1 = spsi_d(1)
           !
           es(1) = es(1) - spsi1 * g1
           es(2) = es(2) - spsi1 * ppsi1
           !
        END IF
        !
        CALL mp_sum( es , intra_bgrp_comm )
        !
        es(1) = es(1) / es(2)
        !
        es1 = es(1)
        !$cuf kernel do
        DO i=1, npwx
           g_d(i) = g_d(i) - es1 * ppsi_d(i)
        END DO
        !
        ! ... e1 = <y| S P^2 PHP|y> / <y| S S P^2|y>  ensures that  
        ! ... <g| S P^2|y> = 0
        !
        ! ... orthogonalize to lowest eigenfunctions (already calculated)
        !
        ! ... scg is used as workspace
        !
        CALL s_1psi_gpu( npwx, npw, g_d(1), scg_d(1) )
        !
        lagrange(1:m-1) = 0.d0
        call divide(inter_bgrp_comm,m-1,m_start,m_end); !write(*,*) m-1,m_start,m_end
        if(m_start.le.m_end) &
        CALL cgcudaDGEMV( 'T', npw2, m_end-m_start+1, 2.D0, psi_d(1,m_start), npw2, scg_d, 1, 0.D0, lagrange_d(m_start), 1 )
        if(m_start.le.m_end) lagrange( m_start:m_end ) = lagrange_d( m_start:m_end )
        CALL mp_sum( lagrange( 1:m-1 ), inter_bgrp_comm )
        IF ( gstart == 2 ) THEN
           psi_aux(1:m-1) = psi_d(1,1:m-1)
           scg1 = scg_d(1)
           lagrange(1:m-1) = lagrange(1:m-1) - psi_aux(1:m-1) * scg1
        END IF
        !
        CALL mp_sum( lagrange( 1:m-1 ), intra_bgrp_comm )
        !
        DO j = 1, ( m - 1 )
           !
           aux = lagrange(j)
           !$cuf kernel do(1)
           DO i = 1, npwx
              g_d(i)   = g_d(i)   - aux * psi_d(i,j)
              scg_d(i) = scg_d(i) - aux * psi_d(i,j)
           END DO
           !
        END DO
        !
        IF ( iter /= 1 ) THEN
           !
           ! ... gg1 is <g(n+1)|S|g(n)> (used in Polak-Ribiere formula)
           !
           gg1 = 2.D0 * cgddot( npw2, g_d(1), 1, g0_d(1), 1 )
           IF ( gstart == 2 ) THEN 
              g1 = g_d(1) ; g01 = g0_d(1)
              gg1 = gg1 - g1 * g01
           END IF
           !
           CALL mp_sum(  gg1 , intra_bgrp_comm )
           !
        END IF
        !
        ! ... gg is <g(n+1)|S|g(n+1)>
        !
        !$cuf kernel do
        do i=1, npwx
           g0_d(i) = scg_d(i)
        end do
        !
        !$cuf kernel do
        do i=1, npw
          g0_d(i) = g0_d(i) * precondition_d(i)
        end do
        !
        gg = 2.D0 * cgddot( npw2, g_d(1), 1, g0_d(1), 1 )
        IF ( gstart == 2 ) THEN 
           g1 = g_d(1) ; g01 = g0_d(1)
           gg = gg - g1*g01
        END IF
        !
        CALL mp_sum(  gg , intra_bgrp_comm )
        !
        IF ( iter == 1 ) THEN
           !
           ! ... starting iteration, the conjugate gradient |cg> = |g>
           !
           gg0 = gg
           !
           !$cuf kernel do
           DO i=1, npwx
              cg_d(i) = g_d(i)
              ! ... |cg> contains now the conjugate gradient
              ! ... set Im[ cg(G=0) ] -  needed for numerical stability
              IF ( gstart == 2 .and. i == 1 ) cg_d(1) = CMPLX( DBLE(cg_d(1)), 0.D0 ,kind=DP)
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
           !$cuf kernel do
           do i=1, npwx
              cg_d(i) = g_d(i) + cg_d(i) * gamma
           end do
           !
           ! ... The following is needed because <y(n+1)| S P^2 |cg(n+1)> 
           ! ... is not 0. In fact :
           ! ... <y(n+1)| S P^2 |cg(n)> = sin(theta)*<cg(n)|S|cg(n)>
           !
           psi_norm = gamma * cg0 * sint
           !
           !$cuf kernel do
           do i=1, npwx
              cg_d(i) = cg_d(i) - psi_norm * psi_d(i,m)
              ! ... |cg> contains now the conjugate gradient
              ! ... set Im[ cg(G=0) ] -  needed for numerical stability
              IF ( gstart == 2 .and. i == 1 ) cg_d(1) = CMPLX( DBLE(cg_d(1)), 0.D0 ,kind=DP)
           end do
           !
        END IF
        !
        ! ... |scg> is S|cg>
        !
        CALL hs_1psi_gpu( npwx, npw, cg_d(1), ppsi_d(1), scg_d(1) )
        !
        cg0 = 2.D0 * cgddot( npw2, cg_d(1), 1, scg_d(1), 1 )
        IF ( gstart == 2 ) THEN 
           cg1 = cg_d(1) ; scg1 = scg_d(1)
           cg0 = cg0 - cg1*scg1
        END IF
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
        a0 = 4.D0 * cgddot( npw2, psi_d(1,m), 1, ppsi_d(1), 1 )
        IF ( gstart == 2 ) THEN 
           psi1 = psi_d(1,m)
           ppsi1 = ppsi_d(1)
           a0 = a0 - 2.D0 * psi1 * ppsi1
        END IF
        !
        a0 = a0 / cg0
        !
        CALL mp_sum(  a0 , intra_bgrp_comm )
        !
        b0 = 2.D0 * cgddot( npw2, cg_d(1), 1, ppsi_d(1), 1 )
        IF ( gstart == 2 ) THEN 
           cg1 = cg_d(1)
           ppsi1 = ppsi_d(1)
           b0 = b0 - cg1 * ppsi1
        END IF
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
        !$cuf kernel do
        do i=1, npwx
           psi_d(i,m) = cost * psi_d(i,m) + sint / cg0 * cg_d(i)
        end do
        !
        ! ... here one could test convergence on the energy
        !
        IF ( ABS( e(m) - e0 ) < ethr_m ) EXIT iterate
        !
        ! ... upgrade H|psi> and S|psi>
        !
        !$cuf kernel do
        do i=1, npwx
           spsi_d(i) = cost * spsi_d(i) + sint / cg0 * scg_d(i)
           !
           hpsi_d(i) = cost * hpsi_d(i) + sint / cg0 * ppsi_d(i)
        end do
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
        IF ( e(m) - e(m-1) < - 2.D0 * ethr_m ) THEN
           !
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
           !$cuf kernel do
           do l=1, npwx
              ppsi_d(l) = psi_d(l,m)
           end do
           !
           DO j = m, i + 1, - 1
              !
              e(j) = e(j-1)
              !
              !$cuf kernel do
              do l=1, npwx
                 psi_d(l,j) = psi_d(l,j-1)
              end do
              !
           END DO
           !
           e(i) = e0
           !
           !$cuf kernel do
           do l=1, npwx
              psi_d(l,i) = ppsi_d(l)
           end do
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
  e_d(1:nbnd) = e(1:nbnd)
  !
  DEALLOCATE( lagrange )
  DEALLOCATE( lagrange_d )
  DEALLOCATE( e )
  DEALLOCATE( psi_aux )
  DEALLOCATE( ppsi_d )
  DEALLOCATE( g0_d )
  DEALLOCATE( cg_d )
  DEALLOCATE( g_d )
  DEALLOCATE( hpsi_d )
  DEALLOCATE( scg_d )
  DEALLOCATE( spsi_d )
  !
  CALL stop_clock( 'rcgdiagg' )
  !
  RETURN
  !
END SUBROUTINE rcgdiagg_gpu
