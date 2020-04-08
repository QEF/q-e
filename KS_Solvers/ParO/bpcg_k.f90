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
! The following file is for solving step 2 of the parallel orbital updating algorithm.
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE bpcg_k( hs_psi, g_1psi, psi0, spsi0, npw, npwx, nbnd, npol, nvec, psi, hpsi, spsi, ethr, e, nhpsi )
  !----------------------------------------------------------------------------
  !
  ! Block Preconditioned Conjugate Gradient solution of the linear system
  !
  !      [ H - e S ]|\tilde\psi> = Pc [ e S - H ] |psi>
  !
  ! the search targets the space orthogonal to the current best wfcs (psi0);
  ! the solution is sought until the residual norm is a fixed fraction of the RHS norm
  ! in this way the more accurate is the original problem the more accuratly the correction is computed
  !
  ! in order to avoid un-necessary HSpsi evaluations this version assumes psi,hpsi and spsi are all
  ! provided in input and return their estimate for further use
  !
  USE util_param,    ONLY : DP, stdout
  USE mp_bands_util, ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! Following varibales are temporary

  COMPLEX(DP),INTENT(IN) :: psi0(npwx*npol,nbnd)  ! psi0  needed to compute the Pv projection
  COMPLEX(DP),INTENT(IN) :: spsi0(npwx*npol,nbnd) ! Spsi0  needed to compute the Pv projection
  INTEGER,  INTENT(IN)   :: npw, npwx, nbnd, npol, nvec ! input dimensions 
  REAL(DP), INTENT(IN)   :: ethr                  ! threshold for convergence.
  REAL(DP), INTENT(INOUT)   :: e(nvec)            ! current estimate of the target eigenvalues
  COMPLEX(DP),INTENT(INOUT) :: psi(npwx*npol,nvec),hpsi(npwx*npol,nvec),spsi(npwx*npol,nvec) ! 
                                                  ! input: the current estimate of the wfcs
                                                  ! output: the estimated correction vectors
  INTEGER, INTENT(INOUT) :: nhpsi                 ! (updated) number of Hpsi evaluations
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 5 ! maximum number of CG iterations 
  !
  COMPLEX(DP), ALLOCATABLE ::  b(:,:),                        & ! RHS for testing
                               p(:,:), hp(:,:), sp(:,:), z(:,:) ! additional working vetors

  COMPLEX(DP), ALLOCATABLE ::  spsi0vec (:,:) ! the product of spsi0 and a group of vectors

  REAL(DP), ALLOCATABLE :: g0(:), g1(:), g2(:), alpha(:), gamma(:), ethr_cg(:), ff(:), ff0(:)
  INTEGER, ALLOCATABLE :: cg_iter(:)
  REAL(DP) :: beta, ee
  INTEGER  :: kdim, kdmx, i, l, block_size, done, nactive, nnew, newdone
  !
  REAL(DP), EXTERNAL :: DDOT
  COMPLEX(DP), EXTERNAL :: ZDOTC

  EXTERNAL  hs_psi, g_1psi
  ! hs_1psi( npwx, npw, psi, hpsi, spsi )
  ! hs_psi( npwx, npw, nvec, psi, hpsi, spsi )
  !
  CALL start_clock( 'pcg' ); !write (6,*) ' enter pcg' , e(1:2), 'npol = ', npol ; FLUSH(6)
  !
  kdim  = npwx*(npol-1) + npw
  kdmx  = npwx* npol
  block_size = min(nvec,64) 
  !
  ALLOCATE( g0( block_size ), g1( block_size ), g2( block_size ), alpha( block_size ), gamma( block_size ) )
  ALLOCATE( ethr_cg( block_size ), ff( block_size ), ff0( block_size ), cg_iter( block_size ) )
  ALLOCATE( z( kdmx, block_size ), b( kdmx, block_size ) )
  ALLOCATE( p(kdmx,block_size), hp(kdmx,block_size), sp(kdmx,block_size) )
  ALLOCATE( spsi0vec(nbnd, block_size) )
  !
  done    = 0  ! the number of correction vectors already solved
  nactive = 0  ! the number of correction vectors currently being updated
  cg_iter = 0  ! how many iteration each active vector has completed (<= maxter)

  MAIN_LOOP: & ! This is a continuous loop. It terminates only when nactive vanishes
  DO
     nnew = min(done+block_size,nvec)-(done+nactive) ! number of new corrections to be added to the seach

     if ( nnew > 0 ) then    ! add nnew vectors to the active list
        !write(6,*) ' nnew =', nnew
        do l=nactive+1,nactive+nnew; i=l+done
           !write(6,*) ' l =',l,' i =',i
           !write (6,*) ' enter pcg' , e(i), 'npol = ', npol ; FLUSH(6)
           b(:,l) = e(i) * spsi(:,i) - hpsi(:,i)               ! initial gradient and saved RHS for later
           z(:,l) = b(:,l); call g_1psi(npwx,npw,z(:,l),e(i)) ! initial preconditioned gradient
        end do
     !- project on conduction bands
        CALL start_clock( 'pcg:ortho' )
        CALL ZGEMM( 'C', 'N', nbnd, nnew, kdim, ONE, spsi0, kdmx, z(:,nactive+1), kdmx, ZERO, spsi0vec, nbnd )
        CALL mp_sum( spsi0vec, intra_bgrp_comm )
        CALL ZGEMM( 'N', 'N', kdim, nnew, nbnd, (-1.D0,0.D0), psi0, kdmx, spsi0vec, nbnd, ONE, z(:,nactive+1), kdmx )
        CALL stop_clock( 'pcg:ortho' )
     !-
        do l=nactive+1,nactive+nnew; i=l+done
           g0(l) = DDOT( 2*kdim, z(:,l), 1, b(:,l), 1 )
        end do
        CALL mp_sum( g0(nactive+1:nactive+nnew), intra_bgrp_comm ) ! g0 = < initial z | initial gradient b >
     
        do l=nactive+1,nactive+nnew; i=l+done
           !write(6,*) ' l =',l,' i =',i
           ff(l) = 0.d0 ; ff0(l) = ff(l)
           !write (6,*) 0, g0(l), ff(l)

           ! ethr_cg = ethr  ! CG convergence threshold could be set from input but it is not ...
           ethr_cg(l) = 1.0D-2  ! it makes more sense to fix the convergence of the CG solution to a 
                                ! fixed function of the RHS (see ethr_cg update later).
           ethr_cg(l) = max ( 0.01*ethr, ethr_cg(l) * g0(l) ) ! here we set the convergence of the correction
           !write(6,*) 'ethr_cg :', ethr_cg(l)

           ! zero the trial solution
           psi(:,i) = ZERO ; hpsi(:,i) = ZERO ; spsi(:,i) = ZERO
           ! initial search direction
           p(:,l) = z(:,l)

           cg_iter(l) = 0 ! this is a new correction vector, reset its interation count
        end do
        nactive = nactive + nnew
     end if
     !write(6,*) ' done  =',done, ' nactive  =',nactive

!  iterate: !  DO cg_iter = 1, maxter  ! THIS IS THE ENTRY POINT OF THE PCG LOOP

     if ( nactive == 0 ) EXIT MAIN_LOOP               ! this is the only MAIN_LOOP EXIT condition

     cg_iter(1:nactive) = cg_iter(1:nactive) + 1      ! update interation counters

     CALL start_clock( 'pcg:hs_1psi' )
!     do l = 1, nactive  ! THIS COULD/SHOULD BE A GLOBAL CALL (ONLY WITHIN ONE BGRP THOUGH)
!        CALL hs_1psi( npwx, npw, p(:,l), hp(:,l), sp(:,l) ) ! apply H to a single wavefunction (no bgrp parallelization here!)
!     end do
     CALL hs_psi( npwx, npw, nactive, p, hp, sp ) ! apply H to a single wavefunction (no bgrp parallelization here!)
     CALL stop_clock( 'pcg:hs_1psi' )

     do l = 1, nactive; i=l+done
        gamma(l) = DDOT( 2*kdim, p(:,l), 1, hp(:,l), 1 ) - e(i) * DDOT( 2*kdim, p(:,l), 1, sp(:,l), 1 ) 
     end do
     CALL mp_sum( gamma(1:nactive), intra_bgrp_comm ) ! gamma = < p | hp - e sp >

     do l = 1, nactive; i=l+done
        !write(6,*) ' l =',l,' i =',i

        alpha(l) = g0(l)/gamma(l)
        !write(6,*) 'g0, gamma, alpha :', g0(l), gamma(l), alpha(l)

        psi(:,i)  = psi(:,i)  + alpha(l) * p(:,l)     ! updated solution
        hpsi(:,i) = hpsi(:,i) + alpha(l) * hp(:,l)    ! updated solution
        spsi(:,i) = spsi(:,i) + alpha(l) * sp(:,l)    ! updated solution

        g2(l) = DDOT(2*kdim,z(:,l),1,b(:,l),1) + e(i) * DDOT(2*kdim,z(:,l),1,spsi(:,i),1) - DDOT(2*kdim,z(:,l),1,hpsi(:,i),1)
     end do
     CALL mp_sum( g2(1:nactive), intra_bgrp_comm )    ! g2 = < old z | new gradient b + e spsi - hpsi >
     do l = 1, nactive; i=l+done                      ! update the preconditioned gradient
        z(:,l) = b(:,l) + e(i) * spsi(:,i) - hpsi(:,i); call g_1psi(npwx,npw,z(:,l),e(i))
     end do
  !- project on conduction bands
     CALL start_clock( 'pcg:ortho' )
     CALL ZGEMM( 'C', 'N', nbnd, nactive, kdim, ONE, spsi0, kdmx, z, kdmx, ZERO, spsi0vec, nbnd )
     CALL mp_sum( spsi0vec, intra_bgrp_comm )
     CALL ZGEMM( 'N', 'N', kdim, nactive, nbnd, (-1.D0,0.D0), psi0, kdmx, spsi0vec, nbnd, ONE, z, kdmx )
     CALL stop_clock( 'pcg:ortho' )
  !-
     do l = 1, nactive; i=l+done
        g1(l) = DDOT(2*kdim,z(:,l),1,b(:,l),1) + e(i) * DDOT(2*kdim,z(:,l),1,spsi(:,i),1) - DDOT(2*kdim,z(:,l),1,hpsi(:,i),1)
     end do
     CALL mp_sum( g1(1:nactive), intra_bgrp_comm )   ! g1 = < new z | new gradient b + e spsi - hpsi >

     do l = 1, nactive; i = l + done                 ! evaluate the function ff
        ff(l) = -0.5_DP * ( e(i)*DDOT(2*kdim,psi(:,i),1,spsi(:,i),1)-DDOT(2*kdim,psi(:,i),1,hpsi(:,i),1) ) &
                - DDOT(2*kdim,psi(:,i),1,b(:,l),1)
     end do
     CALL mp_sum( ff(1:nactive), intra_bgrp_comm )   ! function minimum -0.5 < psi | e spsi - hpsi > - < psi | b >

     newdone = 0  ! number of correction vectors that converge (or are done) at this iteration
     do l = 1, nactive; i = l + done
        !write (6,*) cg_iter(l), g1(l), ff(l),  gamma(l)

        IF ( ff(l) > ff0(l) .AND. ff0(l) < 0.d0 ) THEN
           psi(:,i)  = psi(:,i)  - alpha(l) * p(:,l) ! fallback solution: if last iter failed to improve ff0
           hpsi(:,i) = hpsi(:,i) - alpha(l) * hp(:,l)! exit whitout updating and ...
           spsi(:,i) = spsi(:,i) - alpha(l) * sp(:,l)! hope next time it'll be better
        END IF

        !write(6,*) 'g0, g1, g2 :', g0(l), g1(l), g2(l)
        !write(6,*) 'ff0, ff : ', ff0(l), ff(l)
        IF ( ABS ( g1(l) ) < ethr_cg(l) .OR. ( ff(l) > ff0(l) ) .OR. cg_iter(l) == maxter) THEN ! EXIT iterate
           !write (6,*) ' exit  pcg loop'
           !write(6,*) ' l =',l,' i =',i
           !if ( cg_iter(l) == maxter.and. ABS(g1(l)) > ethr_cg(l))  write (6,*) 'CG not converged maxter exceeded', cg_iter(l), g1(l), g0(l), ethr_cg(l)
           !IF ( ABS ( g1(l) ) < ethr_cg(l)) write (6,*) 'CG correction converged ', cg_iter(l), g1(l), ethr_cg(l)
           !IF ( ABS ( g1(l) ) > g0(l)     ) write (6,*) 'CG not converged ', cg_iter(l), g1(l), g0(l), ethr_cg(l)

           nhpsi = nhpsi + cg_iter(l)    ! update nhpsi count
           IF (.NOT. (ABS(g1(l))< ethr_cg(l) .OR. (ff(l)>ff0(l)) ) .AND. cg_iter(l)==maxter) nhpsi = nhpsi + 1 ! because this would be the count

           newdone = newdone + 1         ! one more solution found (or no more active anyway)
           !write(6,*) ' newdone = ', newdone
       
           CALL start_clock( 'pcg:move' )
           !write(6,*) ' swapping converged psi/hpsi/spsi i = ',i, " with i' = ",done+newdone
           ! swap the terminated vector with the first in the list of the active ones
           p (:,l) = psi (:,done+newdone) ; psi (:,done+newdone) = psi (:,i) ; psi (:,i) = p (:,l)
           hp(:,l) = hpsi(:,done+newdone) ; hpsi(:,done+newdone) = hpsi(:,i) ; hpsi(:,i) = hp(:,l)
           sp(:,l) = spsi(:,done+newdone) ; spsi(:,done+newdone) = spsi(:,i) ; spsi(:,i) = sp(:,l)
           ee      = e(done+newdone)      ; e(done+newdone)      = e(i)      ; e(i)      = ee

           !write(6,*) ' overwrite converged p/hp/etc l = ',l, ' with newdone = ',newdone
           ! move information of the swapped active vector in the right place to keep going
           p(:,l) = p(:,newdone)          ; hp(:,l) = p(:,newdone)           ; sp(:,l) = sp(:,newdone)
           b(:,l) = b(:,newdone) ; z(:,l) = z(:,newdone) ; ff0(l) = ff0(newdone) ; ff(l) = ff(newdone)
           alpha(l) = alpha(newdone) ; g0(l) = g0(newdone) ; g1(l) = g1(newdone) ; g2(l) = g2(newdone)
           cg_iter(l) = cg_iter(newdone) ; ethr_cg(l) = ethr_cg(newdone)
           CALL stop_clock( 'pcg:move' )

        ELSE

           !write(6,*) ' l =',l,' i =',i
           beta   = (g1(l)-g2(l))/g0(l)         ! Polak - Ribiere style update
           g0(l)  = g1(l)                       ! < new z | new gradient >  ->  < old z | old gradient >
           p(:,l) = z(:,l) + beta * p(:,l)      ! updated search direction
           !write(6,*) 'beta :', beta

           ff0(l) = ff(l)                       ! updated minimum value reached by the function

        END IF
     end do
    
     IF ( newdone > 0 ) THEN

        done    = done    + newdone
        nactive = nactive - newdone

        !write(6,*) ' there have been ', newdone, ' new converged solution'
        !write(6,*) ' done = ', done, ' nactive =', nactive

        CALL start_clock( 'pcg:move' )
        do l=1, nactive

           !write(6,*) ' l+newdone =',l+newdone,'  ->   l =',l
           p (:,l) = p (:,l+newdone) ; hp(:,l) = hp(:,l+newdone) ; sp(:,l) = sp(:,l+newdone)
           b(:,l) = b(:,l+newdone) ; z(:,l) = z(:,l+newdone) ; ff0(l) = ff0(l+newdone) ; ff(l) = ff(l+newdone)
           g0(l) = g0(l+newdone) ; g1(l) = g1(l+newdone) ; g2(l) = g2(l+newdone)
           cg_iter(l) = cg_iter(l+newdone) ; ethr_cg(l) = ethr_cg(l+newdone)

        end do
        CALL stop_clock( 'pcg:move' )

     END IF

!  END DO iterate      Here is where the pcg loop would terminate

  END DO  MAIN_LOOP
  !write (6,*) ' exit  pcg loop'

  DEALLOCATE( spsi0vec )
  DEALLOCATE( b, p, hp, sp, z )
  DEALLOCATE( ethr_cg, ff, ff0, cg_iter )
  DEALLOCATE( g0, g1, g2, alpha, gamma )
  !
  CALL stop_clock( 'pcg' )
  !
  RETURN
  !
END SUBROUTINE bpcg_k
