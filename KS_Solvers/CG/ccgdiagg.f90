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
! define __VERBOSE to print a message after each eigenvalue is computed
!
!civn: Dec, 19th 2024
!      merging rcgdiagg_gpu with rcgdiagg in a unique source file 
!      with OpenACC directives. 
!      go back to commit 751be151c67f2edea352f9a34107fe87706258cb 
!      if you need to see the old two separate files
!
!----------------------------------------------------------------------------
SUBROUTINE ccgdiagg( hs_1psi_ptr, s_1psi_ptr, precondition, &
     npwx, npw, nbnd, npol, psi, e, btype, &
     ethr, maxter, reorder, notconv, avg_iter )
  !----------------------------------------------------------------------------
  !
  ! ... "poor man" iterative diagonalization of a complex hermitian matrix
  ! ... through preconditioned conjugate gradient algorithm
  ! ... Band-by-band algorithm with minimal use of memory
  ! ... Calls hs_1psi and s_1psi to calculate H|psi> + S|psi> and S|psi>
  ! ... Works for generalized eigenvalue problem (US pseudopotentials) as well
  !
  USE util_param,     ONLY : DP
  USE mp_bands_util,  ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp,             ONLY : mp_sum
#if defined(__VERBOSE)
  USE util_param,     ONLY : stdout
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
  REAL(DP),    INTENT(IN)    :: precondition(npwx*npol), ethr
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nbnd)
  REAL(DP),    INTENT(INOUT) :: e(nbnd)
  INTEGER,     INTENT(OUT)   :: notconv
  REAL(DP),    INTENT(OUT)   :: avg_iter
  !
  ! ... local variables
  !
  INTEGER                  :: i, j, m, m_start, m_end, iter, moved
  COMPLEX(DP), ALLOCATABLE :: lagrange(:)
  COMPLEX(DP), ALLOCATABLE :: hpsi(:), spsi(:), g(:), cg(:), &
                              scg(:), ppsi(:), g0(:)  
  REAL(DP)                 :: psi_norm, a0, b0, gg0, gamma, gg, gg1, &
                              cg0, e0, es(2), es_1
  REAL(DP)                 :: theta, cost, sint, cos2t, sin2t
  LOGICAL                  :: reorder
  INTEGER                  :: kdim, kdmx, kdim2, ierr, istat
  REAL(DP)                 :: empty_ethr, ethr_m
  !
  ! ... external functions
  !
  REAL (DP), EXTERNAL :: MYDDOT
  EXTERNAL  hs_1psi_ptr, s_1psi_ptr
  ! hs_1psi_ptr( npwx, npw, psi, hpsi, spsi )
  ! s_1psi_ptr( npwx, npw, psi, spsi )
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
  ALLOCATE(  hpsi(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate hpsi ', ABS(ierr) )
  ALLOCATE(  spsi(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate spsi ', ABS(ierr) )
  ALLOCATE(  g(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate g ', ABS(ierr) )
  ALLOCATE(  cg(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate cg ', ABS(ierr) )
  ALLOCATE(  scg(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate scg ', ABS(ierr) )
  ALLOCATE(  ppsi(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate ppsi ', ABS(ierr) )
  ALLOCATE(  g0(kdmx), STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate g0 ', ABS(ierr) )
  ALLOCATE(  lagrange(nbnd), STAT=ierr  )
  IF( ierr /= 0 ) &
       CALL errore( ' ccgdiagg ',' cannot allocate lagrange ', ABS(ierr) )
  !$acc enter data create(hpsi, spsi, g, g0, cg, scg, ppsi, lagrange)
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
     !$acc kernels
     spsi     = ZERO
     scg      = ZERO
     hpsi     = ZERO
     g        = ZERO
     cg       = ZERO
     g0       = ZERO
     ppsi     = ZERO
     lagrange = ZERO
     !$acc end kernels
     !
     ! ... calculate S|psi>
     !
     CALL s_1psi_ptr( npwx, npw, psi(1,m), spsi )
     !
     ! ... orthogonalize starting eigenfunction to those already calculated
     !
     call divide(inter_bgrp_comm,m,m_start,m_end); !write(*,*) m,m_start,m_end
     !$acc host_data use_device(psi, spsi, lagrange)
     IF(m_start.le.m_end) THEN
       CALL MYZGEMV( 'C', kdim, m_end-m_start+1, ONE, psi(1,m_start), &
                   kdmx, spsi, 1, ZERO, lagrange(m_start), 1 )
     END IF
     CALL mp_sum( lagrange, 1, m, inter_bgrp_comm )
     CALL mp_sum( lagrange, 1, m, intra_bgrp_comm )
     !$acc end host_data
     !
     !$acc kernels copyin(m)
     psi_norm = DBLE( lagrange(m) )
     DO j = 1, m - 1
        psi_norm = psi_norm - ( DBLE( lagrange(j) )**2 + AIMAG( lagrange(j) )**2 )
     END DO
     psi_norm = SQRT( psi_norm )
     !$acc end kernels
     !
     !$acc parallel
     !$acc loop gang
     DO i = 1, kdmx
        !$acc loop seq 
        DO j = 1, m - 1
           psi(i,m)  = psi(i,m) - lagrange(j) * psi(i,j)
        END DO
     END DO
     !$acc end parallel 
     !
     !$acc kernels
     psi(:,m) = psi(:,m) / psi_norm
     !$acc end kernels
     !
     ! ... calculate starting gradient (|hpsi> = H|psi>) ...
     !
     CALL hs_1psi_ptr( npwx, npw, psi(1,m), hpsi, spsi )
     !
     ! ... and starting eigenvalue (e = <y|PHP|y> = <psi|H|psi>)
     !
     ! ... NB:  ddot(2*npw,a,1,b,1) = REAL( zdotc(npw,a,1,b,1) )
     !
     !$acc host_data use_device(psi, hpsi)
     e(m) = MYDDOT( kdim2, psi(1,m), 1, hpsi, 1 )
     !$acc end host_data
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
        !$acc kernels
        g(:)    = hpsi(:) / precondition(:)
        ppsi(:) = spsi(:) / precondition(:)
        !$acc end kernels
        !
        ! ... ppsi is now S P(P^2)|y> = S P^2|psi>)
        !
        !$acc host_data use_device(spsi, g, ppsi)
        es(1) = MYDDOT( kdim2, spsi(1), 1, g(1), 1 )
        es(2) = MYDDOT( kdim2, spsi(1), 1, ppsi(1), 1 )
        !$acc end host_data
        !
        CALL mp_sum( es , intra_bgrp_comm )
        !
        es(1) = es(1) / es(2)
        es_1 = es(1)
        !
        !$acc kernels
        g(:) = g(:) - es_1 * ppsi(:)
        !$acc end kernels
        !
        ! ... e1 = <y| S P^2 PHP|y> / <y| S S P^2|y>  ensures that 
        ! ... <g| S P^2|y> = 0
        ! ... orthogonalize to lowest eigenfunctions (already calculated)
        !
        ! ... scg is used as workspace
        !
        CALL s_1psi_ptr( npwx, npw, g(1), scg(1) )
        !
        !$acc kernels
        lagrange(1:m-1) = ZERO
        !$acc end kernels
        call divide(inter_bgrp_comm,m-1,m_start,m_end); !write(*,*) m-1,m_start,m_end
        !$acc host_data use_device(psi, scg, lagrange)
        IF(m_start.le.m_end) THEN
          CALL MYZGEMV( 'C', kdim, m_end-m_start+1, ONE, psi(1,m_start), &
                      kdmx, scg, 1, ZERO, lagrange(m_start), 1 )
        END IF
        CALL mp_sum( lagrange, 1, m-1, inter_bgrp_comm )
        CALL mp_sum( lagrange, 1, m-1, intra_bgrp_comm )
        !$acc end host_data
        !
        !$acc host_data use_device(psi, lagrange, g, scg)
        Call MYZGEMV('N', kdmx, (m-1), -(1.0d0, 0.0d0), psi, kdmx, lagrange, 1, (1.0d0, 0.0d0), g,   1)
        Call MYZGEMV('N', kdmx, (m-1), -(1.0d0, 0.0d0), psi, kdmx, lagrange, 1, (1.0d0, 0.0d0), scg, 1)
        !$acc end host_data
        !
        IF ( iter /= 1 ) THEN
           !
           ! ... gg1 is <g(n+1)|S|g(n)> (used in Polak-Ribiere formula)
           !
           !$acc host_data use_device(g, g0)
           gg1 = MYDDOT( kdim2, g(1), 1, g0(1), 1 )
           !$acc end host_data
           !
           CALL mp_sum( gg1, intra_bgrp_comm )
           !
        END IF
        !
        ! ... gg is <g(n+1)|S|g(n+1)>
        !
        !$acc kernels
        g0(:) = scg(:) * precondition(:)
        !$acc end kernels
        !
        !$acc host_data use_device(g, g0)
        gg = MYDDOT( kdim2, g(1), 1, g0(1), 1 )
        !$acc end host_data
        !
        CALL mp_sum( gg, intra_bgrp_comm )
        !
        IF ( iter == 1 ) THEN
           !
           ! ... starting iteration, the conjugate gradient |cg> = |g>
           !
           gg0 = gg
           !
           !$acc kernels 
           cg(:) = g(:)
           !$acc end kernels
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
           ! See comment below
           !!DO i = 1, kdmx
           !!   cg(i) = g(i) + cg(i) * gamma
           !!END DO
           !
           ! ... The following is needed because <y(n+1)| S P^2 |cg(n+1)> 
           ! ... is not 0. In fact :
           ! ... <y(n+1)| S P^2 |cg(n)> = sin(theta)*<cg(n)|S|cg(n)>
           !
           psi_norm = gamma * cg0 * sint
           !
           ! v== this breaks the logic, done here for performance
           !$acc kernels
           cg(:) = (g(:) + cg(:) * gamma) - psi_norm * psi(:,m)
           !$acc end kernels
           !
        END IF
        !
        ! ... |cg> contains now the conjugate gradient
        !
        ! ... |scg> is S|cg>
        !
        CALL hs_1psi_ptr( npwx, npw, cg(1), ppsi(1), scg(1) )
        !
        !$acc host_data use_device(scg, cg)
        cg0 = MYDDOT( kdim2, cg(1), 1, scg(1), 1 )
        !$acc end host_data
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
        !$acc host_data use_device(psi, ppsi)
        a0 = 2.D0 *  MYDDOT( kdim2, psi(1,m), 1, ppsi(1), 1 ) / cg0
        !$acc end host_data
        !
        CALL mp_sum(  a0 , intra_bgrp_comm )
        !
        !$acc host_data use_device(cg, ppsi)
        b0 = MYDDOT( kdim2, cg(1), 1, ppsi(1), 1 ) / cg0**2
        !$acc end host_data
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
        !$acc kernels
        psi(:,m) = cost * psi(:,m) + sint / cg0 * cg(:)
        !$acc end kernels
        !
        ! ... here one could test convergence on the energy
        !
        IF ( ABS( e(m) - e0 ) < ethr_m ) EXIT iterate
        !
        ! ... upgrade H|psi> and S|psi>
        !
        !$acc kernels
        spsi(:) = cost * spsi(:) + sint / cg0 * scg(:)
        !
        hpsi(:) = cost * hpsi(:) + sint / cg0 * ppsi(:)
        !$acc end kernels
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
           !$acc kernels
           ppsi(:) = psi(:,m)
           !$acc end kernels
           !
           DO j = m, i + 1, - 1
              !
              e(j) = e(j-1)
              !
              !$acc kernels 
              psi(:,j) = psi(:,j-1)
              !$acc end kernels
              !
           END DO
           !
           e(i) = e0
           !
           !$acc kernels 
           psi(:,i) = ppsi(:)
           !$acc end kernels
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
  !$acc exit data delete(hpsi, spsi, g, g0, cg, scg, ppsi, lagrange)
  DEALLOCATE( lagrange )
  DEALLOCATE( ppsi )
  DEALLOCATE( g0 )
  DEALLOCATE( cg )
  DEALLOCATE( g )
  DEALLOCATE( hpsi )
  DEALLOCATE( scg )
  DEALLOCATE( spsi )
  !
  CALL stop_clock( 'ccgdiagg' )
  !
  RETURN
  !
END SUBROUTINE ccgdiagg
