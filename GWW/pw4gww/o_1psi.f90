!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!

SUBROUTINE o_rcgdiagg( npwx, npw, nbnd, psi, e, precondition, &
                     ethr, maxter, reorder, notconv, avg_iter ,&
                      numv, v_states,hdiag,ptype,fcw_number,fcw_state,fcw_mat)
  !----------------------------------------------------------------------------
  !
  ! ... "poor man" iterative diagonalization of a complex hermitian matrix
  ! ... through preconditioned conjugate gradient algorithm
  ! ... Band-by-band algorithm with minimal use of memory

  !
  USE constants, ONLY : pi
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : gstart
  USE mp_global, ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE mp_world,  ONLY : world_comm
  USE fft_base,  ONLY : dffts
  USE io_global, ONLY :stdout
   
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,      INTENT(IN)    :: npwx, npw, nbnd, maxter
   REAL (DP),    INTENT(IN)    :: precondition(npw), ethr
  COMPLEX (DP), INTENT(INOUT) :: psi(npwx,nbnd)
  REAL (DP),    INTENT(INOUT) :: e(nbnd)
  INTEGER,      INTENT(OUT)   :: notconv
  REAL (DP),    INTENT(OUT)   :: avg_iter

  INTEGER, INTENT(in) :: numv!number of valence states                                                               
  REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,numv)!valence states in real space
  REAL(kind=DP), INTENT(in) :: hdiag(npw)!inverse of estimation of diagonal part of hamiltonian
  INTEGER, INTENT(in) :: ptype!type of approximation for O operator
  INTEGER, INTENT(in) :: fcw_number!number of "fake conduction" states for O matrix method 
  COMPLEX(kind=DP) :: fcw_state(npw,fcw_number)! "fake conduction" states for O matrix method
  REAL(kind=DP) :: fcw_mat(fcw_number,fcw_number)! "fake conduction" matrix                                                            


  !
  ! ... local variables
  !
  INTEGER                   :: i, j, m, iter, moved,ig
  REAL (DP),    ALLOCATABLE :: lagrange(:)
  COMPLEX (DP), ALLOCATABLE :: hpsi(:), spsi(:), g(:), cg(:), &
                               scg(:), ppsi(:), g0(:)  
  REAL (DP)                 :: psi_norm, a0, b0, gg0, gamma, gg, gg1, &
                               cg0, e0, es(2)
  REAL (DP)                 :: theta, cost, sint, cos2t, sin2t
  LOGICAL                   :: reorder
  INTEGER                   :: npw2, npwx2
  REAL (DP)                 :: empty_ethr,sca
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
     spsi(1:npw)=psi(1:npw,m)
    
     !
     ! ... orthogonalize starting eigenfunction to those already calculated
     !
     CALL DGEMV( 'T', npw2, m, 2.D0, psi, npwx2, spsi, 1, 0.D0, lagrange, 1 )
     !
     IF ( gstart == 2 ) lagrange(1:m) = lagrange(1:m) - psi(1,1:m) * spsi(1)
     !
     CALL mp_sum( lagrange( 1:m ), world_comm)
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
     !CALL h_1psi( npwx, npw, psi(1,m), hpsi, spsi )
     call o_1psi_gamma( numv, v_states, psi(1,m), hpsi,.false.,hdiag, ptype,fcw_number,fcw_state,fcw_mat,ethr)
     spsi(1:npw)=psi(1:npw,m)
     

    !
     ! ... and starting eigenvalue (e = <y|PHP|y> = <psi|H|psi>)
     !
     ! ... NB:  ddot(2*npw,a,1,b,1) = DBLE( zdotc(npw,a,1,b,1) )
     !
     e(m) = 2.D0 * ddot( npw2, psi(1,m), 1, hpsi, 1 )
     !
     IF ( gstart == 2 ) e(m) = e(m) - psi(1,m) * hpsi(1)
     !
     CALL mp_sum( e(m), world_comm)
     !
     ! ... start iteration for this band
     !
     iterate: DO iter = 1, maxter
        !
        ! ... calculate  P (PHP)|y>
        ! ... ( P = preconditioning matrix, assumed diagonal )
        !
        g(1:npw)    = hpsi(1:npw)! / precondition(:)
        ppsi(1:npw) = spsi(1:npw)! / precondition(:)
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
        CALL mp_sum(  es, world_comm )
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
        !CALL s_1psi( npwx, npw, g(1), scg(1) )
        scg(1:npw)=g(1:npw)
        !
        CALL DGEMV( 'T', npw2, ( m - 1 ), 2.D0, &
                    psi, npwx2, scg, 1, 0.D0, lagrange, 1 )
        !
        IF ( gstart == 2 ) &
           lagrange(1:m-1) = lagrange(1:m-1) - psi(1,1:m-1) * scg(1)
        !
        CALL mp_sum( lagrange( 1 : m-1 ), world_comm)
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
           CALL mp_sum(  gg1, world_comm )
           !
        END IF
        !
        ! ... gg is <g(n+1)|S|g(n+1)>
        !
        g0(:) = scg(:)
        !
        g0(1:npw) = g0(1:npw)! * precondition(:)
        !
        gg = 2.D0 * ddot( npw2, g(1), 1, g0(1), 1 )
        !
        IF ( gstart == 2 ) gg = gg - g(1) * g0(1)
        !
        CALL mp_sum(  gg, world_comm )
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
        !CALL h_1psi( npwx, npw, cg(1), ppsi(1), scg(1) )
        call o_1psi_gamma( numv, v_states, cg, ppsi,.false.,hdiag, ptype,fcw_number,fcw_state,fcw_mat,ethr)
        sca=0.d0
        do ig=1,npw
           sca=sca+2.d0*dble(conjg(cg(ig))*ppsi(ig))
        enddo
        if(gstart==2) sca=sca-dble(conjg(cg(1))*ppsi(1))
        call mp_sum(sca, world_comm)
       
     
        scg(1:npw)=cg(1:npw)
  !
        cg0 = 2.D0 * ddot( npw2, cg(1), 1, scg(1), 1 )
        !
        IF ( gstart == 2 ) cg0 = cg0 - cg(1) * scg(1)
        !
        CALL mp_sum(  cg0, world_comm )
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
        CALL mp_sum(  a0, world_comm )
        !
        b0 = 2.D0 * ddot( npw2, cg(1), 1, ppsi(1), 1 )
        !
        IF ( gstart == 2 ) b0 = b0 - cg(1) * ppsi(1)
        !
        b0 = b0 / cg0**2
        !
        CALL mp_sum(  b0, world_comm )
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
       
           !
        IF ( ABS( e(m) - e0 ) < ethr ) EXIT iterate
           !
        
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
END SUBROUTINE o_rcgdiagg



!
!----------------------------------------------------------------------------
SUBROUTINE o_1psi_gamma( numv, v_states, psi, opsi,l_freq,hdiag, ptype,fcw_number,fcw_state,fcw_mat,ethr)
  !----------------------------------------------------------------------------
  !
  !
  !this subroutines applies the O oprator to a state psi
  !IT WORKS ONLY FOR NORMCONSERVING PSEUDOPOTENTIALS
  !the valence states in G space must be in evc  
  ! Gamma point version

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE kinds,    ONLY : DP
   USE wannier_gw
   USE gvect
   USE constants, ONLY : e2, pi, tpi, fpi
   USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
   USE wvfct,    ONLY : g2kin, npwx, npw, nbnd, et
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : world_comm, mpime, nproc
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE kinds, ONLY : DP
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE becmod,           ONLY : becp,allocate_bec_type,deallocate_bec_type
   USE uspp,                 ONLY : vkb, nkb, okvan
   USE g_psi_mod,            ONLY : h_diag, s_diag
   USE klist,                ONLY : xk,igk_k
   
  !
  IMPLICIT NONE

  INTEGER, INTENT(in) :: numv!number of valence states
  REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,numv)!valence states in real space
  COMPLEX(kind=DP), INTENT(in) :: psi(npw)!input wavefunction
  COMPLEX(kind=DP), INTENT(out) :: opsi(npw)!O|\psi>
  LOGICAL, INTENT(in) :: l_freq!if true estimates the operator a 0 frequency
  REAL(kind=DP), INTENT(in) :: hdiag(npw)!inverse of estimation of diagonal part of hamiltonian
  INTEGER, INTENT(in) :: ptype!type of approximation for O operator
  INTEGER, INTENT(in) :: fcw_number!number of "fake conduction" states for O matrix method
  COMPLEX(kind=DP) :: fcw_state(npw,fcw_number)! "fake conduction" states for O matrix method
  REAL(kind=DP) :: fcw_mat(fcw_number,fcw_number)! "fake conduction" matrix
  REAL(kind=DP), INTENT(in) :: ethr!threshold on (H-e) inversion
  

  REAL(kind=DP), ALLOCATABLE :: psi_r(:,:), psi_v(:)
  COMPLEX(kind=DP), ALLOCATABLE :: psi_g(:,:), h_psi_g(:,:),s_psi_g(:,:)
  REAL(kind=DP) :: ec
  INTEGER :: iv,ii,jj,ig

  REAL(kind=DP),ALLOCATABLE  :: p_terms(:),s_terms(:)
  INTEGER :: l_blk,nbegin,nend,nsize
  !REAL(kind=DP), ALLOCATABLE :: h_diag (:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: psi_g2(:,:)
  INTEGER :: kter
  LOGICAL :: lconv_root,lfirst
  REAL(kind=DP) :: anorm
  EXTERNAL :: hpsi_pw4gww,cg_psi_pw4gww
  COMPLEX(kind=DP), POINTER, SAVE :: tmp_psi(:,:)


  allocate(psi_r(dffts%nnr,2),psi_v(dffts%nnr))
  if(pmat_type/=0) then
     allocate(h_psi_g(npw,2),s_psi_g(npw,2),psi_g(npw,2))
  else
     allocate(h_psi_g(npw,numv),s_psi_g(npw,numv),psi_g(npw,numv))
  endif
  
!fourier transform psi to R space

  opsi(1:npw)=(0.d0,0.d0)


  if(pmat_type==0 .or. pmat_type==1 .or. pmat_type==2) then
     psic(:)=(0.d0,0.d0)
     psic(nls(1:npw))  = psi(1:npw)
     psic(nlsm(1:npw)) = CONJG( psi(1:npw) )
     CALL invfft ('Wave', psic, dffts)
     psi_v(:)= DBLE(psic(:))
  endif
  
  if(pmat_type==0) then
     call start_clock('opsi_total')
     call allocate_bec_type ( nkb, numv, becp)
     IF ( nkb > 0 )  CALL init_us_2( npw, igk_k(1,1), xk(1,1), vkb )
     if(.not.associated(tmp_psi)) then
        allocate( tmp_psi(npw,num_nbndv(1)))
        lfirst=.true.
     else
        lfirst=.false.
     endif
     allocate (h_diag(npw, numv),s_diag(npw,numv))
     allocate(psi_g2(npw,numv))
     !      
     ! compute the kinetic energy
     !                                                                                                                                                   
        do ig = 1, npw
           g2kin (ig) = ( g (1,ig)**2 + g (2,ig)**2 + g (3,ig)**2 ) * tpiba2
        enddo
        h_diag=0.d0
        do iv = 1, numv
           do ig = 1, npw
              !h_diag(ig,ibnd)=1.d0/max(1.0d0,g2kin(ig)/eprec(ibnd,ik))
              !h_diag(ig,iv) = 1.D0 + g2kin(ig) + &
              !     SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
              h_diag(ig,iv)=g2kin(ig)
              !h_diag(ig,iv)=1.d0/max(1.0d0,g2kin(ig)/8.d0)
           enddo
        enddo
      do iv=1,numv,2
!!product with psi_v
         if(iv/=numv) then
           psi_r(1:dffts%nnr,1)=psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv)
           psi_r(1:dffts%nnr,2)=psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv+1)
        else
           psi_r(1:dffts%nnr,1)=psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv)
        endif
!!fourier transfrm to G
        if(iv/=numv) then
           psic(1:dffts%nnr)=cmplx(psi_r(1:dffts%nnr,1),psi_r(1:dffts%nnr,2))
        else
           psic(1:dffts%nnr)=cmplx(psi_r(1:dffts%nnr,1),0.d0)
        endif
        call start_clock('opsi_fft')
        CALL fwfft ('Wave', psic, dffts)
        call stop_clock('opsi_fft')
        if(iv/=numv) then
           psi_g(1:npw,iv)=0.5d0*(psic(nls(1:npw))+conjg(psic(nlsm(1:npw))))
           psi_g(1:npw,iv+1)=(0.d0,-0.5d0)*(psic(nls(1:npw))-conjg(psic(nlsm(1:npw))))
           if(gstart==2) psi_g(1,iv)=dble(psi_g(1,iv))
           if(gstart==2) psi_g(1,iv+1)=dble(psi_g(1,iv+1))
        else
           psi_g(1:npw,iv)=psic(nls(1:npw))
           if(gstart==2) psi_g(1,iv)=dble(psi_g(1,iv))
        endif

!!project on conduction manifold                                
        call start_clock('opsi_pc')
        if(iv/=numv) then
           call pc_operator(psi_g(:,iv),1,.false.)
           call pc_operator(psi_g(:,iv+1),1,.false.)
        else
           call pc_operator(psi_g(:,iv),1,.false.)
        endif
        call stop_clock('opsi_pc')
     enddo

     write(stdout,*) 'DEBUG1'
     FLUSH(stdout)
!call (H-e)^-1 solver
     if(.true.) then
        psi_g2(1:npw,1:numv)=psi_g(1:npw,1:numv)
     else
        psi_g2(1:npw,1:numv)=tmp_psi(1:npw,1:numv)
     endif
     write(stdout,*) 'DEBUG1.5'
     FLUSH(stdout)
     call cgsolve_all_gamma (hpsi_pw4gww,cg_psi_pw4gww,et(1,1),psi_g,psi_g2, &
              h_diag,npw,npw,ethr,1,kter,lconv_root,anorm,numv,1)

     tmp_psi(1:npw,1:numv)=psi_g2(1:npw,1:numv)
     write(stdout,*) 'DEBUG2',kter,lconv_root,anorm
     FLUSH(stdout)


     do iv=1,numv,2
        !!project on conduction manifold                                                                                                                            
        call start_clock('opsi_pc')
        if(iv/=numv) then
           call pc_operator(psi_g2(:,iv),1,.false.)
           call pc_operator(psi_g2(:,iv+1),1,.false.)
        else
           call pc_operator(psi_g2(:,iv),1,.false.)
        endif
        call stop_clock('opsi_pc')


        !!fourier transform to R space 
        psic(:)=(0.d0,0.d0)
        if(iv/=numv) then
           psic(nls(1:npw))  = psi_g2(1:npw,iv)+(0.d0,1.d0)*psi_g2(1:npw,iv+1)
           psic(nlsm(1:npw)) = CONJG( psi_g2(1:npw,iv) )+(0.d0,1.d0)*conjg(psi_g2(1:npw,iv+1))
        else
           psic(nls(1:npw))  = psi_g2(1:npw,iv)
           psic(nlsm(1:npw)) = CONJG( psi_g2(1:npw,iv) )
        endif
        call start_clock('opsi_fft')
        CALL invfft ('Wave', psic, dffts)
        call stop_clock('opsi_fft')
        if(iv/=numv) then
           psi_r(:,1)= DBLE(psic(:))
           psi_r(:,2)= dimag(psic(:))
        else
           psi_r(:,1)= DBLE(psic(:))
        endif
!!product with psi_v          
        if(iv/=numv) then
           psi_r(1:dffts%nnr,1)=psi_r(1:dffts%nnr,1)*v_states(1:dffts%nnr,iv)
           psi_r(1:dffts%nnr,2)=psi_r(1:dffts%nnr,2)*v_states(1:dffts%nnr,iv+1)
        else
           psi_r(1:dffts%nnr,1)=psi_r(1:dffts%nnr,1)*v_states(1:dffts%nnr,iv)
        endif
!!fourier transform in G space!! sum up results
!!TAKE CARE OF SIGN 
        if(iv/=numv) then
           psic(:)=cmplx(psi_r(:,1),psi_r(:,2))
        else
           psic(:)=cmplx(psi_r(:,1),0.d0)
        endif
        CALL fwfft ('Wave', psic, dffts)
        if(iv/=numv) then
           opsi(1:npw)=opsi(1:npw)-0.5d0*(psic(nls(1:npw))+conjg(psic(nlsm(1:npw))))
           opsi(1:npw)=opsi(1:npw)-(0.d0,-0.5d0)*(psic(nls(1:npw))-conjg(psic(nlsm(1:npw))))
        else
           opsi(1:npw)=opsi(1:npw)-psic(nls(igk_k(1:npw,1)))
        endif
     enddo
     deallocate(h_diag,s_diag)
     deallocate(psi_g2)
     call deallocate_bec_type(becp)
     call stop_clock('opsi_total')
    ! call print_clock('opsi_total')
    ! call print_clock('opsi_fft')
    ! call print_clock('opsi_pc')                   
    
  else if(pmat_type==1) then
!loop on v
     call start_clock('opsi_total')
     do iv=1,numv,2
!!product with psi_v
        if(iv/=numv) then
           psi_r(1:dffts%nnr,1)=psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv)
           psi_r(1:dffts%nnr,2)=psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv+1)
        else
           psi_r(1:dffts%nnr,1)=psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv)
        endif
!!fourier transfrm to G
        if(iv/=numv) then
           psic(1:dffts%nnr)=cmplx(psi_r(1:dffts%nnr,1),psi_r(1:dffts%nnr,2))
        else
           psic(1:dffts%nnr)=cmplx(psi_r(1:dffts%nnr,1),0.d0)
        endif
        call start_clock('opsi_fft')
        CALL fwfft ('Wave', psic, dffts)
        call stop_clock('opsi_fft')
        if(iv/=numv) then
           psi_g(1:npw,1)=0.5d0*(psic(nls(1:npw))+conjg(psic(nlsm(1:npw))))
           psi_g(1:npw,2)=(0.d0,-0.5d0)*(psic(nls(1:npw))-conjg(psic(nlsm(1:npw))))
           if(gstart==2) psi_g(1,1)=dble(psi_g(1,1))
           if(gstart==2) psi_g(1,2)=dble(psi_g(1,2))
        else
           psi_g(1:npw,1)=psic(nls(1:npw))
           if(gstart==2) psi_g(1,1)=dble(psi_g(1,1))
        endif

!!project on conduction manifold
        call start_clock('opsi_pc')
        if(iv/=numv) then
           call pc_operator(psi_g(:,1),1,.false.)
           call pc_operator(psi_g(:,2),1,.false.)
        else
           call pc_operator(psi_g(:,1),1,.false.)
        endif
     !call pc_operator_test(psi_g)
        if(l_freq) then
          
          !    call cgsolve_all_gamma (h_psi, cg_psi, e, d0psi, dpsi, h_diag, &
          !         ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd, npol)
          
          
        endif
        call stop_clock('opsi_pc')
!!fourier transform to R space
        psic(:)=(0.d0,0.d0)
        if(iv/=numv) then
           psic(nls(1:npw))  = psi_g(1:npw,1)+(0.d0,1.d0)*psi_g(1:npw,2)
           psic(nlsm(1:npw)) = CONJG( psi_g(1:npw,1) )+(0.d0,1.d0)*conjg(psi_g(1:npw,2))
        else
           psic(nls(1:npw))  = psi_g(1:npw,1)
           psic(nlsm(1:npw)) = CONJG( psi_g(1:npw,1) )
        endif
        call start_clock('opsi_fft')
        CALL invfft ('Wave', psic, dffts)
        call stop_clock('opsi_fft')
        if(iv/=numv) then
           psi_r(:,1)= DBLE(psic(:))
           psi_r(:,2)= dimag(psic(:))
        else
           psi_r(:,1)= DBLE(psic(:))
        endif
!!product with psi_v

        if(iv/=numv) then
           psi_r(1:dffts%nnr,1)=psi_r(1:dffts%nnr,1)*v_states(1:dffts%nnr,iv)
           psi_r(1:dffts%nnr,2)=psi_r(1:dffts%nnr,2)*v_states(1:dffts%nnr,iv+1)
        else
           psi_r(1:dffts%nnr,1)=psi_r(1:dffts%nnr,1)*v_states(1:dffts%nnr,iv)
        endif
!!fourier transform in G space
!! sum up results
!!TAKE CARE OF SIGN

        if(iv/=numv) then
           psic(:)=cmplx(psi_r(:,1),psi_r(:,2))
        else
           psic(:)=cmplx(psi_r(:,1),0.d0)
        endif
        CALL fwfft ('Wave', psic, dffts)
        if(iv/=numv) then
           opsi(1:npw)=opsi(1:npw)-0.5d0*(psic(nls(1:npw))+conjg(psic(nlsm(1:npw))))
           opsi(1:npw)=opsi(1:npw)-(0.d0,-0.5d0)*(psic(nls(1:npw))-conjg(psic(nlsm(1:npw))))
        else
           opsi(1:npw)=opsi(1:npw)-psic(nls(igk_k(1:npw,1)))
        endif

     enddo
     call stop_clock('opsi_total')
    ! call print_clock('opsi_total')
    ! call print_clock('opsi_fft')
    ! call print_clock('opsi_pc')
  else if(pmat_type==2) then
     psi_r(:,1)=0.d0
     do iv=1,numv
        psi_r(1:dffts%nnr,1)=psi_r(1:dffts%nnr,1)+psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv)
     enddo
     
     psic(:)=cmplx(psi_r(:,1),0.d0)
     CALL fwfft ('Wave', psic, dffts)
     psi_g(1:npw,1)=psic(nls(igk_k(1:npw,1)))
     if(gstart==2) psi_g(1,1)=dble(psi_g(1,1))
     
     call pc_operator(psi_g(:,1),1,.false.)
     psi_g(1:npw,1)=psi_g(1:npw,1)*hdiag(1:npw)
     call pc_operator(psi_g(:,1),1,.false.)

     psic(nls(igk_k(1:npw,1)))  = psi_g(1:npw,1)
     psic(nlsm(igk_k(1:npw,1))) = CONJG( psi_g(1:npw,1) )
     CALL invfft ('Wave', psic, dffts)
     psi_r(:,1)=dble(psic(:))
  
     psi_v(:)= psi_r(:,1)
     psi_r(:,1)=0.d0
     do iv=1,numv
        psi_r(1:dffts%nnr,1)=psi_r(1:dffts%nnr,1)+psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv)
     enddo

     psic(:)=cmplx(psi_r(:,1),0.d0)
     CALL fwfft ('Wave', psic, dffts)
     opsi(1:npw)=-psic(nls(igk_k(1:npw,1)))
    
  else!cases 3,4
!form scalar products
     allocate(p_terms(fcw_number),s_terms(fcw_number))
     call dgemm('T','N',fcw_number,1,2*npw,2.d0,fcw_state,2*npw,psi,2*npw,0.d0,p_terms,fcw_number)
     if(gstart==2) then
        do ii=1,fcw_number
           p_terms(ii)=p_terms(ii)-dble(conjg(fcw_state(1,ii))*psi(1))
        enddo
     endif
     call mp_sum(p_terms(:),world_comm)
     
!multiply to D matrix     
     l_blk= (fcw_number)/nproc
     if(l_blk*nproc < (fcw_number)) l_blk = l_blk+1
     nbegin=mpime*l_blk+1
     nend=nbegin+l_blk-1
     if(nend > fcw_number) nend=fcw_number
     nsize=nend-nbegin+1
   

     s_terms(:)=0.d0
     if(nsize>0) then
        call dgemm('T','N',nsize,1,fcw_number,1.d0,fcw_mat,fcw_number,p_terms,fcw_number,0.d0,&
          &s_terms(nbegin:nend),nsize)
     endif
  
!collect from processors
     call mp_sum(s_terms,world_comm)

!multiply with gamma functions
     call dgemm('N','N',2*npw,1,fcw_number,-1.d0,fcw_state,2*npw,s_terms,fcw_number,0.d0,opsi,2*npw)
     
   
     deallocate(p_terms,s_terms)
  endif

  if(gstart==2) opsi(1)=(0.d0,0.d0)
  !
  deallocate(psi_r, psi_g,psi_v)
  deallocate(h_psi_g,s_psi_g)
  !
  RETURN
  !
END SUBROUTINE o_1psi_gamma


SUBROUTINE evc_to_real(numv, v_states)
!this subroutine fourier transform states from evc
!to real space

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE kinds,    ONLY : DP
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE wavefunctions_module, ONLY : evc, psic
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE klist,  ONLY : igk_k

   implicit none

  INTEGER, INTENT(in) :: numv!number of states to be transformed
  REAL(kind=DP), INTENT(out) :: v_states(dffts%nnr,numv)!target arrsy
  
  INTEGER :: iv

  do iv=1,numv,2
     psic(:)=(0.d0,0.d0)
     if(iv < numv) then
        psic(nls(igk_k(1:npw,1)))  = evc(1:npw,iv)+(0.d0,1.d0)*evc(1:npw,iv+1)
        psic(nlsm(igk_k(1:npw,1))) = CONJG( evc(1:npw,iv) )+(0.d0,1.d0)*CONJG(evc(1:npw,iv+1))
     else
        psic(nls(igk_k(1:npw,1)))  = evc(1:npw,iv)
        psic(nlsm(igk_k(1:npw,1))) = CONJG( evc(1:npw,iv) )
     endif
     CALL invfft ('Wave', psic, dffts)
     if(iv<numv) then
        v_states(1:dffts%nnr,iv)=dble(psic(1:dffts%nnr))
        v_states(1:dffts%nnr,iv+1)=dimag(psic(1:dffts%nnr))
     else
        v_states(1:dffts%nnr,iv)=dble(psic(1:dffts%nnr))
     endif

  enddo
return

END SUBROUTINE evc_to_real


SUBROUTINE o_basis_init(numpw,o_basis,numv,v_states,cutoff, ptype,fcw_number,fcw_state,fcw_mat,ethr)
!this subroutine initializes the polarization basis to
!random values and diagonalize with respect to the O matrix


  USE gvect, ONLY : g
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE kinds,    ONLY : DP
  USE wvfct,    ONLY : g2kin, npwx, npw, nbnd
  USE wavefunctions_module, ONLY : evc, psic
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  USE constants, ONLY : tpi
  USE random_numbers, ONLY : randy
  USE klist,                ONLY : xk,igk_k
  USE cell_base,            ONLY : tpiba2
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  

   implicit none

   INTEGER, INTENT(in) :: numv!number of states to be transformed
   REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,numv)!valence states in real space
   INTEGER, INTENT(in) :: numpw!dimesion of the polarization basis
   COMPLEX(kind=DP), INTENT(out) :: o_basis(npw,numpw)!polarization basis
   REAL(kind=DP), INTENT(in) :: cutoff!cutoff for planewaves
   INTEGER, INTENT(in) :: ptype!type of approximation for O operator
   INTEGER, INTENT(in) :: fcw_number!number of "fake conduction" states for O matrix method
   COMPLEX(kind=DP) :: fcw_state(npw,fcw_number)! "fake conduction" states for O matrix method
   REAL(kind=DP) :: fcw_mat(fcw_number,fcw_number)! "fake conduction" matrix
   REAL(kind=DP), INTENT(in) :: ethr!threshold o_1psi_gamma


   INTEGER :: iw,ig
   REAL(kind=DP) :: rr, arg
   REAL(kind=DP), ALLOCATABLE :: e(:)
   REAL(kind=DP), ALLOCATABLE :: hdiag(:)


  
   allocate(hdiag(npw))

   g2kin(1:npw) = ( (g(1,igk_k(1:npw,1)) )**2 + &
        ( g(2,igk_k(1:npw,1)) )**2 + &
        ( g(3,igk_k(1:npw,1)) )**2 ) * tpiba2

  
   do ig=1,npw
      !if(g2kin(ig) <= 3.d0) then
      if(g2kin(ig) <= cutoff) then
         hdiag(ig)=1.d0!/(1.d0+g2kin(ig))
      else
        hdiag(ig)=0.d0
     endif
  enddo
  hdiag(:)=1.d0
  

   allocate( e(numpw))

   DO iw = 1, numpw
      !
      DO ig = 1, npw
                   !
         rr  = randy()
         
         arg = tpi * randy()
                   !
         o_basis(ig,iw) = CMPLX( rr*COS( arg ), rr*SIN( arg ) ) / &
              ( g(1,igk_k(ig,1))**2 +  g(2,igk_k(ig,1))**2 + g(3,igk_k(ig,1))**2 + 1.D0) 
                   !
      END DO
                !
   END DO

   !CALL rinitcgg( npwx, npw, n_starting_wfc, &
   !                         nbnd, wfcatom, wfcatom, etatom )

  

   call o_rinitcgg( npwx, npw, numpw, numpw, o_basis, o_basis, e, numv, v_states,hdiag,ptype,fcw_number,fcw_state,fcw_mat,ethr)
   do iw=1,numpw
      write(stdout,*) 'E', iw, e(iw)
   enddo

   deallocate(e)
   deallocate(hdiag)
   return

 END SUBROUTINE o_basis_init


 SUBROUTINE o_basis_write(numpw, o_basis,lcutoff,ecutoff, l_pbc)

!this subroutines writes the basis of the polarization on file
   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE kinds,    ONLY : DP
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE wavefunctions_module, ONLY : evc, psic
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE wavefunctions_module, ONLY : psic
   USE io_files,  ONLY : diropn, prefix
   USE gvect,     ONLY : ngm, gg,gstart
   USE cell_base, ONLY: tpiba2
   USE wannier_gw, ONLY : max_ngm
   USE gvect, ONLY : nl
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE klist,   ONLY : igk_k

   implicit none 

   INTEGER, EXTERNAL :: find_free_unit

   INTEGER, INTENT(inout) :: numpw!dimension of the polarization basis
   COMPLEX(kind=DP), INTENT(in) :: o_basis(npw,numpw)
   REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum  
   LOGICAL, INTENT(in) :: lcutoff !if true uses cutoff on G defined by ecutoff
   LOGICAL, INTENT(in) :: l_pbc!if true PBC are assumed and element G=0 is added as the last one
   
   INTEGER :: iw, iungprod, ig,ngm_max
   LOGICAL :: exst
   COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:)

   allocate(tmp_g(max_ngm))

   if(lcutoff) then
      ngm_max=0
      do ig=1,ngm
         if(gg(ig)*tpiba2 >= ecutoff) exit
         ngm_max=ngm_max+1
      enddo
   else
      ngm_max=ngm
   endif

   write(stdout,*) 'NGM MAX:', ngm_max, ngm

   iungprod = find_free_unit()
   CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )



   do iw=1,numpw

      psic(:)=(0.d0,0.d0)
      psic(nls(igk_k(1:npw,1)))  = o_basis(1:npw,iw)
      psic(nlsm(igk_k(1:npw,1))) = CONJG( o_basis(1:npw,iw))
      tmp_g(1:max_ngm)=psic(nl(1:max_ngm))
      if(gstart==2) tmp_g(1)=(0.d0,0.d0)
      CALL davcio(tmp_g, max_ngm*2,iungprod,iw,1)
   enddo

   if(l_pbc) then
      numpw=numpw+1
      tmp_g(1:max_ngm)=(0.d0,0.d0)
      if(gstart==2) tmp_g(1)=(1.d0,0.d0)
      CALL davcio(tmp_g, max_ngm*2,iungprod,numpw,1)
   endif
   close(iungprod)
   deallocate(tmp_g)

   return

 END SUBROUTINE o_basis_write


!----------------------------------------------------------------------------
SUBROUTINE o_1psi_gamma_real( numv, v_states, psi, opsi)
  !----------------------------------------------------------------------------
  !
  !
  !this subroutines applies the O oprator to a state psi
  !IT WORKS ONLY FOR NORMCONSERVING PSEUDOPOTENTIALS
  !the valence states in G space must be in evc
  ! Gamma point version in real space

   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE kinds,    ONLY : DP
   USE wannier_gw
   USE gvect
   USE constants, ONLY : e2, pi, tpi, fpi
   USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE wavefunctions_module, ONLY : evc, psic
   USE mp, ONLY : mp_sum, mp_barrier, mp_bcast
   USE mp_world, ONLY : mpime, world_comm
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE klist,  ONLY : igk_k

  USE kinds, ONLY : DP
  !
  IMPLICIT NONE

  INTEGER, INTENT(in) :: numv!number of valence states
  REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,numv)!valence states in real space
  COMPLEX(kind=DP), INTENT(in) :: psi(npw)!input wavefunction
  COMPLEX(kind=DP), INTENT(out) :: opsi(npw)!O|\psi>

  REAL(kind=DP), ALLOCATABLE :: psi_r(:), psi_v(:), psi_w(:)
  COMPLEX(kind=DP), ALLOCATABLE :: psi_g(:)
  REAL(kind=DP), ALLOCATABLE :: prod(:)
  INTEGER :: iv

  allocate(psi_r(dffts%nnr),psi_g(npw),psi_v(dffts%nnr))
  allocate(prod(numv),psi_w(dffts%nnr))

!fourier transform psi to R space

  opsi(1:npw)=(0.d0,0.d0)

  psic(:)=(0.d0,0.d0)
  psic(nls(igk_k(1:npw,1)))  = psi(1:npw)
  psic(nlsm(igk_k(1:npw,1))) = CONJG( psi(1:npw) )
  CALL invfft ('Wave', psic, dffts)
  psi_v(1:dffts%nnr)= DBLE(psic(1:dffts%nnr))
  psi_w(:)=0.d0
!loop on v
  do iv=1,numv

!!product with psi_v
     psi_r(1:dffts%nnr)=psi_v(1:dffts%nnr)*v_states(1:dffts%nnr, iv)
!!project on conduction manifold
     call dgemm('T','N',numv,1,dffts%nnr,1.d0,v_states,dffts%nnr,psi_r,dffts%nnr,0.d0,prod,numv)
     call mp_sum(prod(1:numv), world_comm)
     prod(:)=prod(:)/dble(dffts%nr1*dffts%nr2*dffts%nr3)
     call dgemm('N','N',dffts%nnr,1,numv,-1.d0,v_states,dffts%nnr,prod,numv,1.d0, psi_r,dffts%nnr)
!


     psi_r(1:dffts%nnr)=psi_r(1:dffts%nnr)*v_states(1:dffts%nnr,iv)

!add up result (with sign)
     psi_w(:)=psi_w(:)-psi_r(:)



  enddo
!!fourier transform in G space


  psic(:)=cmplx(psi_w(:),0.d0)
  CALL fwfft ('Wave', psic, dffts)
  opsi(1:npw)=psic(nls(igk_k(1:npw,1)))


  if(gstart==2) opsi(1)=dble(opsi(1))
  !
  deallocate(psi_r, psi_g,psi_v)
  deallocate(prod,psi_w)
  !
  RETURN
  !
END SUBROUTINE o_1psi_gamma_real



 SUBROUTINE o_basis_test(numv,v_states,numpw, lcutoff,ecutoff)

!this subroutines writes the basis of the polarization on file      
   USE io_global,            ONLY : stdout, ionode, ionode_id
   USE kinds,    ONLY : DP
   USE wvfct,    ONLY : npwx, npw, nbnd
   USE wavefunctions_module, ONLY : evc, psic
   USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE wavefunctions_module, ONLY : psic
   USE io_files,  ONLY : prefix, tmp_dir, diropn
   USE gvect,     ONLY : ngm, gg,gstart
   USE cell_base, ONLY: tpiba2
   USE wannier_gw, ONLY : max_ngm
   USE gvect, ONLY : nl
   USE mp, ONLY : mp_sum
   USE mp_world, ONLY : world_comm
   USE fft_base,             ONLY : dfftp, dffts
   USE fft_interfaces,       ONLY : fwfft, invfft
   USE klist, ONLY : igk_k



   implicit none

   INTEGER, EXTERNAL :: find_free_unit
   INTEGER, INTENT(in) :: numv!number of valence states
   REAL(kind=DP), INTENT(in) :: v_states(dffts%nnr,numv)!valence states in real space             
   INTEGER, INTENT(in) :: numpw!dimension of the polarization basis
   REAL(kind=DP), INTENT(in) :: ecutoff!cutoff in Rydberg for g sum
   LOGICAL, INTENT(in) :: lcutoff !if true uses cutoff on G defined by ecutoff    

   INTEGER :: iw, iungprod, ig,ngm_max
   LOGICAL :: exst
   COMPLEX(kind=DP), ALLOCATABLE :: tmp_g(:), psi1(:),psi2(:)
   REAL(kind=DP) :: sca

   allocate(tmp_g(max_ngm),psi1(npw),psi2(npw))
   
   if(lcutoff) then
      ngm_max=0
      do ig=1,ngm
         if(gg(ig)*tpiba2 >= ecutoff) exit
         ngm_max=ngm_max+1
      enddo
   else
      ngm_max=ngm
   endif

   write(stdout,*) 'NGM MAX:', ngm_max, ngm

   iungprod = find_free_unit()
   CALL diropn( iungprod, 'wiwjwfc_red', max_ngm*2, exst )

   do iw=1,numpw
      CALL davcio(tmp_g, max_ngm*2,iungprod,iw,-1)
      psic(:)=(0.d0,0.d0)
      do ig=1,max_ngm
         psic(nls(ig))=tmp_g(ig)
         psic(nlsm(ig))=conjg(tmp_g(ig))
      enddo
      do ig=1,npw
         psi1(ig)=psic(nls(igk_k(ig,1)))
      enddo

      call o_1psi_gamma( numv, v_states, psi1, psi2)
      sca=0.d0
      do ig=1,npw
         sca=sca+2.d0*dble(conjg(psi2(ig))*psi2(ig))
      enddo
      if(gstart==2) sca=sca-dble(conjg(psi2(1))*psi2(1))
      call mp_sum(sca,world_comm)
      sca=dsqrt(sca)
      psi2(:)=psi2(:)/sca
      sca=0.d0
      do ig=1,npw
         sca=sca+2.d0*dble(conjg(psi1(ig))*psi2(ig))
      enddo
      if(gstart==2) sca=sca-dble(conjg(psi1(1))*psi2(1))
      call mp_sum(sca,world_comm)
      write(stdout,*) 'o basis test:',iw,sca

   enddo


   close(iungprod)
   deallocate(tmp_g,psi1,psi2)

   return

 END SUBROUTINE o_basis_test

