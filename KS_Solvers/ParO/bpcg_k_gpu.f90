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
! * GPU version Ivan Carnimeo
!
! The following file is for solving step 2 of the parallel orbital updating algorithm.
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE bpcg_k_gpu( hs_psi_gpu, g_1psi_gpu, psi0_d, spsi0_d, npw, npwx, nbnd, npol, nvec, &
                                              psi_d, hpsi_d, spsi_d, ethr, e_d, nhpsi )
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
#if defined (__CUDA)
  USE cudafor
#endif
  USE util_param,    ONLY : DP, stdout
  USE mp_bands_util, ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! Following varibales are temporary

  COMPLEX(DP),INTENT(IN) :: psi0_d(npwx*npol,nbnd)  ! psi0  needed to compute the Pv projection
  COMPLEX(DP),INTENT(IN) :: spsi0_d(npwx*npol,nbnd) ! Spsi0  needed to compute the Pv projection
  INTEGER,  INTENT(IN)   :: npw, npwx, nbnd, npol, nvec ! input dimensions 
  REAL(DP), INTENT(IN)   :: ethr                  ! threshold for convergence.
  REAL(DP), INTENT(INOUT)   :: e_d(nvec)            ! current estimate of the target eigenvalues
  COMPLEX(DP),INTENT(INOUT) :: psi_d(npwx*npol,nvec),hpsi_d(npwx*npol,nvec),spsi_d(npwx*npol,nvec) ! 
                                                  ! input: the current estimate of the wfcs
                                                  ! output: the estimated correction vectors
  INTEGER, INTENT(INOUT) :: nhpsi                 ! (updated) number of Hpsi evaluations
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 5 ! maximum number of CG iterations 
  !
  COMPLEX(DP), ALLOCATABLE ::  b_d(:,:),                        & ! RHS for testing
                               p_d(:,:), hp_d(:,:), sp_d(:,:), z_d(:,:) ! additional working vetors

  COMPLEX(DP), ALLOCATABLE ::  spsi0vec_d (:,:) ! the product of spsi0 and a group of vectors

  REAL(DP), ALLOCATABLE :: g0(:), g1(:), g2(:), gamma(:), ethr_cg(:), ff(:), ff0(:)
  REAL(DP), ALLOCATABLE :: alpha(:)
  INTEGER, ALLOCATABLE :: cg_iter(:)
  REAL(DP) :: beta, ee
  INTEGER  :: kdim, kdmx, i, l, block_size, done, nactive, nnew, newdone
  !
  REAL(DP), EXTERNAL :: gpu_DDOT

  EXTERNAL  g_1psi_gpu, hs_psi_gpu
  ! hs_1psi( npwx, npw, psi, hpsi, spsi )
  ! hs_psi( npwx, npw, nvec, psi, hpsi, spsi )
  !
  !
  INTEGER :: ii, jj  ! cuf kernel indeces
  REAL(DP) :: tmp  
#if defined (__CUDA)
  attributes(device) :: psi_d, hpsi_d, spsi_d, psi0_d, spsi0_d, spsi0vec_d
  attributes(device) :: e_d 
  attributes(device) :: b_d, p_d, hp_d, sp_d, z_d
#endif 
  !
  CALL start_clock( 'pcg' ); !write (6,*) ' enter pcg' , e(1:2), 'npol = ', npol ; FLUSH(6)
  !
  kdim  = npwx*(npol-1) + npw
  kdmx  = npwx* npol
  block_size = min(nvec,64) 
  !
  ALLOCATE( g0( block_size ), g1( block_size ), g2( block_size ), gamma( block_size ) )
  ALLOCATE( ethr_cg( block_size ), ff( block_size ), ff0( block_size ), cg_iter( block_size ) )
  ALLOCATE( z_d( kdmx, block_size ), b_d( kdmx, block_size ) )
  ALLOCATE( p_d(kdmx,block_size), hp_d(kdmx,block_size), sp_d(kdmx,block_size) )
  ALLOCATE( spsi0vec_d(nbnd, block_size) )
  ALLOCATE( alpha( block_size ) )
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
!$cuf kernel do(1)
           DO ii = 1, kdmx
             b_d(ii,l) = e_d(i) * spsi_d(ii,i) - hpsi_d(ii,i)               ! initial gradient and saved RHS for later
           END DO 
!$cuf kernel do(1)
           DO ii = 1, kdmx
             z_d(ii,l) = b_d(ii,l)
           END DO 
           call g_1psi_gpu(npwx,npw,z_d(:,l),e_d(i)) ! initial preconditioned gradient
        end do
     !- project on conduction bands
        CALL start_clock( 'pcg:ortho' )
        CALL gpu_ZGEMM( 'C', 'N', nbnd, nnew, kdim, ONE, spsi0_d, kdmx, z_d(:,nactive+1), kdmx, ZERO, &
                                                                                        spsi0vec_d, nbnd )
        CALL mp_sum( spsi0vec_d, intra_bgrp_comm )
        CALL gpu_ZGEMM( 'N', 'N', kdim, nnew, nbnd, (-1.D0,0.D0), psi0_d, kdmx, spsi0vec_d, nbnd, ONE, &
                                                                                      z_d(:,nactive+1), kdmx )
        CALL stop_clock( 'pcg:ortho' )
     !-
        do l=nactive+1,nactive+nnew; i=l+done
           g0(l) = gpu_DDOT( 2*kdim, z_d(:,l), 1, b_d(:,l), 1 )
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
!$cuf kernel do(1)
           DO ii = 1, kdmx
             psi_d(ii,i) = ZERO 
             hpsi_d(ii,i) = ZERO 
             spsi_d(ii,i) = ZERO
           END DO 
           ! initial search direction
!$cuf kernel do(1)
           DO ii = 1, kdmx
             p_d(ii,l) = z_d(ii,l)
           END DO 
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
     CALL hs_psi_gpu( npwx, npw, nactive, p_d, hp_d, sp_d ) ! apply H to a single wavefunction (no bgrp parallelization here!)
     CALL stop_clock( 'pcg:hs_1psi' )
     do l = 1, nactive; i=l+done
        gamma(l) = gpu_DDOT( 2*kdim, p_d(:,l), 1, hp_d(:,l), 1 )
        gamma(l) = gamma(l) - e_d(i) * gpu_DDOT( 2*kdim, p_d(:,l), 1, sp_d(:,l), 1 ) 
     end do
     CALL mp_sum( gamma(1:nactive), intra_bgrp_comm ) ! gamma = < p | hp - e sp >

     do l = 1, nactive; i=l+done
        !write(6,*) ' l =',l,' i =',i

        alpha(l) = g0(l)/gamma(l)
        !write(6,*) 'g0, gamma, alpha :', g0(l), gamma(l), alpha(l)
        tmp = alpha(l)
!$cuf kernel do(1)
        DO ii = 1, kdmx
          psi_d(ii,i)  = psi_d(ii,i)  + tmp * p_d(ii,l)     ! updated solution
          hpsi_d(ii,i) = hpsi_d(ii,i) + tmp * hp_d(ii,l)    ! updated solution
          spsi_d(ii,i) = spsi_d(ii,i) + tmp * sp_d(ii,l)    ! updated solution
        END DO 

        g2(l) = gpu_DDOT(2*kdim,z_d(:,l),1,b_d(:,l),1) &
                + e_d(i) * gpu_DDOT(2*kdim,z_d(:,l),1,spsi_d(:,i),1) &
                - gpu_DDOT(2*kdim,z_d(:,l),1,hpsi_d(:,i),1)
     end do
     CALL mp_sum( g2(1:nactive), intra_bgrp_comm )    ! g2 = < old z | new gradient b + e spsi - hpsi >
     do l = 1, nactive; i=l+done                      ! update the preconditioned gradient
!$cuf kernel do(1)
        DO ii = 1, kdmx
          z_d(ii,l) = b_d(ii,l) + e_d(i) * spsi_d(ii,i) - hpsi_d(ii,i) 
        END DO 
        call g_1psi_gpu(npwx,npw,z_d(:,l),e_d(i))
     end do
  !- project on conduction bands
     CALL start_clock( 'pcg:ortho' )
     CALL gpu_ZGEMM( 'C', 'N', nbnd, nactive, kdim, ONE, spsi0_d, kdmx, z_d, kdmx, ZERO, spsi0vec_d, nbnd )
     CALL mp_sum( spsi0vec_d, intra_bgrp_comm )
     CALL gpu_ZGEMM( 'N', 'N', kdim, nactive, nbnd, (-1.D0,0.D0), psi0_d, kdmx, spsi0vec_d, nbnd, ONE, z_d, kdmx )
     CALL stop_clock( 'pcg:ortho' )
  !-
     do l = 1, nactive; i=l+done
        g1(l) = gpu_DDOT(2*kdim,z_d(:,l),1,b_d(:,l),1) &
                + e_d(i) * gpu_DDOT(2*kdim,z_d(:,l),1,spsi_d(:,i),1) &
                - gpu_DDOT(2*kdim,z_d(:,l),1,hpsi_d(:,i),1)
     end do
     CALL mp_sum( g1(1:nactive), intra_bgrp_comm )   ! g1 = < new z | new gradient b + e spsi - hpsi >
     do l = 1, nactive; i = l + done                 ! evaluate the function ff
        ff(l) = -0.5_DP * ( e_d(i)*gpu_DDOT(2*kdim,psi_d(:,i),1,spsi_d(:,i),1)&
                                                -gpu_DDOT(2*kdim,psi_d(:,i),1,hpsi_d(:,i),1) )
        ff(l) = ff(l) - gpu_DDOT(2*kdim,psi_d(:,i),1,b_d(:,l),1)
     end do
     CALL mp_sum( ff(1:nactive), intra_bgrp_comm )   ! function minimum -0.5 < psi | e spsi - hpsi > - < psi | b >

     newdone = 0  ! number of correction vectors that converge (or are done) at this iteration
     do l = 1, nactive; i = l + done
        !write (6,*) cg_iter(l), g1(l), ff(l),  gamma(l)

        IF ( ff(l) > ff0(l) .AND. ff0(l) < 0.d0 ) THEN
           tmp = alpha(l)
!$cuf kernel do(1)
           DO ii = 1, kdmx
             psi_d(ii,i)  = psi_d(ii,i)  - tmp * p_d(ii,l) ! fallback solution: if last iter failed to improve ff0
             hpsi_d(ii,i) = hpsi_d(ii,i) - tmp * hp_d(ii,l)! exit whitout updating and ...
             spsi_d(ii,i) = spsi_d(ii,i) - tmp * sp_d(ii,l)! hope next time it'll be better
           END DO 
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
!$cuf kernel do(1)
           DO ii = 1, kdmx
             p_d (ii,l) = psi_d (ii,done+newdone) 
             psi_d (ii,done+newdone) = psi_d (ii,i) 
             psi_d (ii,i) = p_d (ii,l)
             hp_d(ii,l) = hpsi_d(ii,done+newdone) 
             hpsi_d(ii,done+newdone) = hpsi_d(ii,i) 
             hpsi_d(ii,i) = hp_d(ii,l)
             sp_d(ii,l) = spsi_d(ii,done+newdone) 
             spsi_d(ii,done+newdone) = spsi_d(ii,i)  
             spsi_d(ii,i) = sp_d(ii,l)
           END DO 

           ee = e_d(done+newdone)       
!$cuf kernel do(1)
           DO ii = 1, 1
             e_d(done+newdone) = e_d(i)      
           END DO 
           e_d(i) = ee

           !write(6,*) ' overwrite converged p/hp/etc l = ',l, ' with newdone = ',newdone
           ! move information of the swapped active vector in the right place to keep going
!$cuf kernel do(1)
           DO ii = 1, kdmx
             p_d(ii,l) = p_d(ii,newdone)          
             hp_d(ii,l) = p_d(ii,newdone)           
             sp_d(ii,l) = sp_d(ii,newdone)
             b_d(ii,l) = b_d(ii,newdone) 
             z_d(ii,l) = z_d(ii,newdone) 
           END DO 

           alpha(l) = alpha(newdone) 

           ff0(l) = ff0(newdone) 
           ff(l) = ff(newdone)
           g0(l) = g0(newdone) 
           g1(l) = g1(newdone) 
           g2(l) = g2(newdone)
           cg_iter(l) = cg_iter(newdone) 
           ethr_cg(l) = ethr_cg(newdone)
           CALL stop_clock( 'pcg:move' )

        ELSE

           !write(6,*) ' l =',l,' i =',i
           beta   = (g1(l)-g2(l))/g0(l)         ! Polak - Ribiere style update
           g0(l)  = g1(l)                       ! < new z | new gradient >  ->  < old z | old gradient >
!$cuf kernel do(1)
           DO ii = 1, kdmx
             p_d(ii,l) = z_d(ii,l) + beta * p_d(ii,l)      ! updated search direction
           END DO 
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
!$cuf kernel do(1)
           DO ii = 1, kdmx
             p_d (ii,l) = p_d (ii,l+newdone) 
             hp_d(ii,l) = hp_d(ii,l+newdone) 
             sp_d(ii,l) = sp_d(ii,l+newdone)
             b_d(ii,l) = b_d(ii,l+newdone)  
             z_d(ii,l) = z_d(ii,l+newdone) 
           END DO 

           ff0(l) = ff0(l+newdone) ; ff(l) = ff(l+newdone)
           g0(l) = g0(l+newdone) ; g1(l) = g1(l+newdone) ; g2(l) = g2(l+newdone)
           cg_iter(l) = cg_iter(l+newdone) ; ethr_cg(l) = ethr_cg(l+newdone)

        end do
        CALL stop_clock( 'pcg:move' )

     END IF

!  END DO iterate      Here is where the pcg loop would terminate
    
  END DO  MAIN_LOOP
  !write (6,*) ' exit  pcg loop'

  DEALLOCATE( spsi0vec_d )
  DEALLOCATE( b_d, p_d, hp_d, sp_d, z_d )
  DEALLOCATE( ethr_cg, ff, ff0, cg_iter )
  DEALLOCATE( g0, g1, g2, gamma )
  DEALLOCATE( alpha )
  !
  CALL stop_clock( 'pcg' )
  !
  RETURN
  !
END SUBROUTINE bpcg_k_gpu
