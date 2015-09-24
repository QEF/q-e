!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE diago_cg(ndim,omat,maxter,max_state,e,ovec,cutoff,ethr,found_state,l_para)
  !----------------------------------------------------------------------------
  !
  ! ... "poor man" iterative diagonalization of a real symmetric  matrix O
  ! ... through preconditioned conjugate gradient algorithm
  ! ... Band-by-band algorithm with minimal use of memory
  !
  USE constants, ONLY : pi
  USE kinds,     ONLY : DP
  USE io_global,        ONLY : stdout
  USE mp_world, ONLY : mpime,nproc,world_comm
  USE mp, ONLY : mp_sum
  USE random_numbers, ONLY : randy
  
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !

  INTEGER, INTENT(in) :: ndim!matrix dimension
  REAL(kind=DP), INTENT(in) :: omat(ndim,ndim)!matrix to be diagonalized
  INTEGER, INTENT(in) ::maxter!maximum number of iterations
  INTEGER, INTENT(in) :: max_state!maximum number of eigenvectors to be found
  REAL(kind=DP),INTENT(inout) :: e(ndim)!eigenvalues
  REAL(kind=DP), INTENT(inout) :: ovec(ndim,max_state)!eigenvector
  REAL(kind=DP),INTENT(in) :: cutoff!found eigenvalues larger than cutoff
  REAL (DP),    INTENT(IN)    ::  ethr!threshold for convergence
  INTEGER, INTENT(out) :: found_state!number of states found
  LOGICAL, INTENT(in) :: l_para!if true omat is distributed among processors

  !
  ! ... local variables
  !
  INTEGER                   :: i, j, m, iter, moved, iw, ig
  REAL (DP),    ALLOCATABLE :: lagrange(:)
  REAL (DP), ALLOCATABLE :: hpsi(:), spsi(:), g(:), cg(:), &
                               scg(:), ppsi(:), g0(:)  
  REAL (DP)                 :: psi_norm, a0, b0, gg0, gamma, gg, gg1, &
                               cg0, e0, es(2)
  REAL (DP)                 :: theta, cost, sint, cos2t, sin2t
  LOGICAL                   :: reorder=.true.
  LOGICAL                   :: l_all_ok, l_first_out
  INTEGER                   :: m_first_out, delta_first_out=10000
  INTEGER :: l_blk,nbegin,nend,nsize
  REAL(kind=DP)::avg_iter
  INTEGER :: notconv
  REAL(kind=DP), ALLOCATABLE :: aux(:,:)
  REAL (DP)                 :: rtmp(2)
  REAL (DP),    ALLOCATABLE :: hr(:,:,:), sr(:,:)
  REAL (DP),    ALLOCATABLE :: en(:),ctmp(:)
  REAL(kind=DP) :: rr
  REAL(kind=DP), ALLOCATABLE :: ovec2(:,:)
 !
  ! ... external functions
  !
  REAL (DP), EXTERNAL :: DDOT

  !
  !
  CALL start_clock( 'diago_cg' )
  !
 
  !
  !
  ALLOCATE( spsi( ndim ) )
  ALLOCATE( scg(  ndim ) )
  ALLOCATE( hpsi( ndim ) )
  ALLOCATE( g(    ndim ) )
  ALLOCATE( cg(   ndim ) )
  ALLOCATE( g0(   ndim ) )
  ALLOCATE( ppsi( ndim ) )
  !    
  ALLOCATE( lagrange( max_state) )
  !
  avg_iter = 0.D0
  notconv  = 0
  moved    = 0
  l_all_ok=.true.
  l_first_out=.false.
  !
  ! ... every eigenfunction is calculated separately
  !
  write(stdout,*) 'ATTENZIONE1'
  FLUSH(stdout)

  l_blk= (ndim)/nproc
  if(l_blk*nproc < ndim) l_blk = l_blk+1
  nbegin=mpime*l_blk+1
  nend=nbegin+l_blk-1
  if(nend > ndim) nend=ndim
  nsize=nend-nbegin+1!WARNING it could be < 1


!initialization
  DO iw = 1, max_state
     DO ig = 1, ndim
         rr  = randy()!rndm()
         ovec(ig,iw)=rr
      END DO
   END DO
   allocate(aux(ndim,2))
   ALLOCATE( hr( max_state, max_state, 2 ) )
   ALLOCATE( sr( max_state, max_state ) )
   ALLOCATE( en( max_state) ,ctmp(max_state))
   DO m = 1, max_state
     call gradient(ovec(1:ndim,m),aux(1:ndim,1))
     aux(:,2)=ovec(:,m)

     if(nsize > 0 )then
        !CALL DGEMV( 'T', nsize, 2, 1.D0, aux(nbegin:nend,1:2), nsize, ovec(nbegin:nend,m), 1, 0.D0, rtmp, 1 )
        CALL DGEMV( 'T', nsize, 2, 1.D0, aux(nbegin,1), ndim, ovec(nbegin,m), 1, 0.D0, rtmp, 1 )
     else
        rtmp(1:2)=0.d0
     endif
     call mp_sum(rtmp(1:2),world_comm)

     hr(m,m,1) = rtmp(1)
     sr(m,m)   = rtmp(2)

     DO j = m + 1, max_state
        if(nsize>0) then
           !CALL DGEMV( 'T', nsize, 2, 1.D0, aux(nbegin:nend,1:2), nsize, ovec(nbegin:nend,j), 1, 0.D0, rtmp, 1 )
           CALL DGEMV( 'T', nsize, 2, 1.D0, aux(nbegin,1), ndim, ovec(nbegin,j), 1, 0.D0, rtmp, 1 )
        else
           rtmp(1:2)=0.d0
        endif
     
        hr(j,m,1) = rtmp(1)
        sr(j,m)   = rtmp(2)

        hr(m,j,1) = rtmp(1)
        sr(m,j)   = rtmp(2)

     END DO

  END DO

  write(stdout,*) 'ATTENZIONE2'
  FLUSH(stdout)

  call mp_sum(hr(:,:,1),world_comm)
  call mp_sum(sr(:,:),world_comm)
  write(stdout,*) 'Call rdiaghg'
  FLUSH(stdout)

  CALL rdiaghg( max_state, max_state, hr, sr, max_state, en, hr(1,1,2) )
  write(stdout,*) 'Done'
  FLUSH(stdout)

  e(1:max_state) = en(1:max_state)

 
!  DO i = 1,ndim

!     DO m = 1, max_state

!        ctmp(m) = SUM( hr(:,m,2) * ovec(i,:) )

!     END DO

!     ovec(i,1:max_state) = ctmp(1:max_state)

!  END DO

  allocate(ovec2(ndim,max_state))
  ovec2(:,:)=ovec(:,:)
  ovec(:,:)=0.d0
  if(nsize > 0) then
     call dgemm('N','N',nsize,max_state,max_state,1.d0,ovec2(nbegin:nend,1:max_state),&
   &nsize,hr(1:max_state,1:max_state,2),max_state,0.d0,ovec(nbegin:nend,1:max_state),nsize)
  endif
  call mp_sum(ovec(:,:),world_comm)
 

  deallocate(ovec2)
  deallocate(aux)
  deallocate(hr,sr)
  deallocate(en,ctmp)



  write(stdout,*) 'ATTENZIONE3'
  FLUSH(stdout)







states:  DO m = 1, max_state
    
   write(stdout,*) 'ATTENZIONE4',m
  FLUSH(stdout)

    
     !
     ! ... calculate S|psi>
     !
     !CALL s_1psi( ndmx, ndim, psi(1,m), spsi )
     spsi(:)=ovec(:,m)
     !
     ! ... orthogonalize starting eigenfunction to those already calculated
     !
     if(nsize>0) then
        !CALL DGEMV( 'T', nsize, m, 1.D0, ovec(nbegin:nend,1:m), nsize, spsi(nbegin:nend), 1, 0.D0, lagrange, 1 )
        CALL DGEMV( 'T', nsize, m, 1.D0, ovec(nbegin,1), SIZE(ovec,1), spsi(nbegin), 1, 0.D0, lagrange, 1 )
     else
        lagrange(:)=0.d0
     endif
     !
     call mp_sum(lagrange(1:m),world_comm)
    
        !
   
     psi_norm = lagrange(m)
     !
     DO j = 1, m - 1
        !
        ovec(:,m)  = ovec(:,m) - lagrange(j) * ovec(:,j)
        !
        psi_norm = psi_norm - lagrange(j)**2
        !
     END DO
     !
     psi_norm = SQRT( psi_norm )
    
          !
     ovec(:,m) = ovec(:,m) / psi_norm
     !
     ! ... calculate starting gradient (|hpsi> = H|psi>) ...
     !
     call gradient(ovec(1:ndim,m),hpsi)
  
    
     spsi(1:ndim)=ovec(1:ndim,m)
     !
     ! ... and starting eigenvalue (e = <y|PHP|y> = <psi|H|psi>)
     !
     ! ... NB:  DDOT(2*ndim,a,1,b,1) = DBLE( ZDOTC(ndim,a,1,b,1) )
     !
     if(nsize>0) then
        e(m) = DDOT( nsize, ovec(nbegin:nend,m), 1, hpsi(nbegin:nend), 1 )
     else
        e(m)=0.d0
     endif
        !
     call mp_sum(e(m),world_comm)
     !
         !
     ! ... start iteration for this band
     ! 
     
     
     iterate: DO iter = 1, maxter
        !
        ! ... calculate  P (PHP)|y>
        ! ... ( P = preconditioning matrix, assumed diagonal )
        !
        g(1:ndim)    = hpsi(1:ndim)! / precondition(:)
        ppsi(1:ndim) = spsi(1:ndim)! / precondition(:)
        !
        ! ... ppsi is now S P(P^2)|y> = S P^2|psi>)
        !
        if(nsize>0) then
           es(1) = DDOT( nsize, spsi(nbegin:nend), 1, g(nbegin:nend), 1 )
           es(2) = DDOT( nsize, spsi(nbegin:nend), 1, ppsi(nbegin:nend), 1 )
        else
           es(1:2)=0.d0
        endif
        call mp_sum(es(1:2),world_comm)
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
        !CALL s_1psi( ndmx, ndim, g(1), scg(1) )
        scg(1:ndim)=g(1:ndim)
        !
        
        if(nsize> 0) then
           !CALL DGEMV( 'T', nsize, ( m - 1 ), 1.D0, &
           !         ovec(nbegin:nend,1:m-1), nsize, scg(nbegin:nend), 1, 0.D0, lagrange, 1 )
           CALL DGEMV( 'T', nsize, ( m - 1 ), 1.D0, ovec(nbegin,1), SIZE(ovec,1), scg(nbegin), 1, 0.D0, lagrange, 1 )
        else
           lagrange(1:m-1)=0.d0
        endif
           !
        call mp_sum(lagrange(1:m-1),world_comm)
        !
        !
        DO j = 1, ( m - 1 )
           !
           g(:)   = g(:)   - lagrange(j) * ovec(:,j)
           scg(:) = scg(:) - lagrange(j) * ovec(:,j)
           !
        END DO
        !
        IF ( iter /= 1 ) THEN
           !
           ! ... gg1 is <g(n+1)|S|g(n)> (used in Polak-Ribiere formula)
           !
           if(nsize>0) then
              gg1 = DDOT( nsize, g(nbegin:nend), 1, g0(nbegin:nend), 1 )
           else
              gg1=0.d0
           endif
          !
           call mp_sum(gg1,world_comm)
           !
           !
        END IF
        !
        ! ... gg is <g(n+1)|S|g(n+1)>
        !
        g0(:) = scg(:)
        !
        g0(1:ndim) = g0(1:ndim)! * precondition(:)
        !
        if(nsize>0) then
           gg = DDOT( nsize, g(nbegin:nend), 1, g0(nbegin:nend), 1 )
        else
           gg=0.d0
        endif
        !
        call mp_sum(gg,world_comm)
        !
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
           cg(:) = cg(:) - psi_norm * ovec(:,m)
           !
        END IF
        !
        ! ... |cg> contains now the conjugate gradient
        !
        ! ... |scg> is S|cg>
        !
         call gradient(cg,ppsi)

        
        scg(1:ndim)=cg(1:ndim)
        !
        if(nsize>0) then
           cg0 = DDOT( nsize, cg(nbegin:nend), 1, scg(nbegin:nend), 1 )
        else
           cg0=0.d0
        endif
        !
        call mp_sum(cg0,world_comm)
        !
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
        if(nsize>0) then
           a0 = 2.D0 * DDOT( nsize, ovec(nbegin:nend,m), 1, ppsi(nbegin:nend), 1 )
        else
           a0=0.d0
        endif
        !
        !
        a0 = a0 / cg0
        !
        call mp_sum(a0,world_comm)
        !
        if(nsize>0) then
           b0 = DDOT( nsize, cg(nbegin:nend), 1, ppsi(nbegin:nend), 1 )
        else
           b0=0.d0
        endif
        !
        !
        b0 = b0 / cg0**2
        !
        call mp_sum(b0,world_comm)

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
        ovec(:,m) = cost * ovec(:,m) + sint / cg0 * cg(:)
        !
        ! ... here one could test convergence on the energy
        !
        IF ( ABS( e(m) - e0 ) < ethr ) THEN
           write(stdout,*) 'State:',m,'Iterations:',iter,e(m)
           FLUSH(stdout)
          
           EXIT iterate
        ELSE
           l_all_ok=.false.
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
           write(stdout,*) 'DO REORDER:',m
           FLUSH(stdout)
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
           ppsi(:) = ovec(:,m)
           !
           DO j = m, i + 1, - 1
              !
              e(j) = e(j-1)
              !
              ovec(:,j) = ovec(:,j-1)
              !
           END DO
           !
           e(i) = e0
           !
           ovec(:,i) = ppsi(:)
           !
           ! ... this procedure should be good if only a few inversions occur,
           ! ... extremely inefficient if eigenvectors are often in bad order
           ! ... ( but this should not happen )
           !
        END IF
        !
     END IF
     if(abs(e(m))<cutoff)    EXIT
     
        !
  END DO states
  found_state=m
  if(found_state>max_state) found_state=max_state
  !
  avg_iter = avg_iter / DBLE( found_state )
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
  CALL stop_clock( 'diago_cg' )

  RETURN

  CONTAINS
    SUBROUTINE gradient(vec,grad)
!apply gradient
      implicit none
      
      REAL(kind=DP), INTENT(in)  :: vec(ndim)
      REAL(kind=DP), INTENT(out) :: grad(ndim)

      grad(:)=0.d0
      if(nsize>0) then
         if(.not.l_para) then
            call dgemm('T','N',nsize,1,ndim,-1.d0,omat(1:ndim,nbegin:nend),ndim,vec,ndim,0.d0,grad(nbegin:nend),nsize)
         else
            call dgemm('T','N',nsize,1,ndim,-1.d0,omat(1:ndim,1:nsize),ndim,vec,ndim,0.d0,grad(nbegin:nend),nsize)
         endif
      endif
      call mp_sum(grad(1:ndim),world_comm)

      return

    END SUBROUTINE gradient
  !
END SUBROUTINE diago_cg
