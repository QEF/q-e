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
! * GPU porting Ivan Carnimeo
!
! The following file is for solving step 2 of the parallel orbital updating algorithm.
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE bpcg_gamma( hs_psi_ptr, g_psi_ptr, psi0, spsi0, npw, npwx, nbnd, nvec, psi, hpsi, spsi, ethr, e, nhpsi )
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
  USE mp_bands_util, ONLY : intra_bgrp_comm, gstart
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! Following variables are temporary

  COMPLEX(DP),INTENT(IN) :: psi0(npwx,nbnd)  ! psi0  needed to compute the Pv projection
  COMPLEX(DP),INTENT(IN) :: spsi0(npwx,nbnd) ! Spsi0  needed to compute the Pv projection
  INTEGER,  INTENT(IN)   :: npw, npwx, nbnd, nvec ! input dimensions 
  REAL(DP), INTENT(IN)   :: ethr                  ! threshold for convergence.
  REAL(DP), INTENT(INOUT)   :: e(nvec)            ! current estimate of the target eigenvalues
  COMPLEX(DP),INTENT(INOUT) :: psi(npwx,nvec),hpsi(npwx,nvec),spsi(npwx,nvec) ! 
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

  REAL(DP), ALLOCATABLE    ::  spsi0vec (:,:) ! the product of spsi0 and a group of vectors
  !$acc declare device_resident(spsi0vec)

  REAL(DP), ALLOCATABLE :: g0(:), g1(:), g2(:), alpha(:), gamma(:), ethr_cg(:), ff(:), ff0(:)
  !$acc declare device_resident(g0, g1, g2, alpha, gamma, ethr_cg, ff, ff0)
  INTEGER, ALLOCATABLE :: cg_iter(:)
  !$acc declare device_resident(cg_iter)
  INTEGER :: cg_iter_l ! cg_iter(l) (useful for some GPU optimization)
  !
  LOGICAL :: ff_check, ff0_check, g1_check, iter_check  ! exit iteration condition checks
  !
  REAL(DP) :: beta, ee
  INTEGER  :: npw2, npwx2, i, l, block_size, done, nactive, nnew, newdone
  !
  REAL(DP), EXTERNAL :: MYDDOT_VECTOR_GPU
  !$acc routine(MYDDOT_VECTOR_GPU) vector
  !

  EXTERNAL  hs_psi_ptr, g_psi_ptr
  ! hs_1psi_ptr( npwx, npw, psi, hpsi, spsi )
  ! hs_psi_ptr( npwx, npw, nvec, psi, hpsi, spsi )
  !
  INTEGER :: ii ! kernel indeces
  !
  CALL start_clock( 'pcg' ); !write (6,*) ' enter pcg' , e(1:2) ; FLUSH(6)
  !
  npw2  = 2*npw
  npwx2 = 2*npwx
  block_size = min(nvec,64) 
  !
  ALLOCATE( g0( block_size ), g1( block_size ), g2( block_size ), alpha( block_size ), gamma( block_size ) )
  ALLOCATE( ethr_cg( block_size ), ff( block_size ), ff0( block_size ), cg_iter( block_size ) )
  ALLOCATE( z( npwx, block_size ), b( npwx, block_size ) )
  ALLOCATE( p(npwx,block_size), hp(npwx,block_size), sp(npwx,block_size) )
  ALLOCATE( spsi0vec(nbnd, block_size) )
  !$acc enter data create(p, hp, sp, z, b)
  !
  done    = 0  ! the number of correction vectors already solved
  nactive = 0  ! the number of correction vectors currently being updated

  !$acc kernels 
  cg_iter = 0  ! how many iteration each active vector has completed (<= maxter)
  !$acc end kernels

  MAIN_LOOP: & ! This is a continuous loop. It terminates only when nactive vanishes
  DO
     nnew = min(done+block_size,nvec)-(done+nactive) ! number of new corrections to be added to the search
     if ( nnew > 0 ) then    ! add nnew vectors to the active list
        !write(6,*) ' nnew =', nnew
        !$acc parallel  
        !$acc loop gang private(i)
        do l=nactive+1,nactive+nnew; i=l+done
           !write(6,*) ' l =',l,' i =',i
           !write (6,*) ' enter pcg' , e(i) ; FLUSH(6)
           !$acc loop vector
           DO ii = 1, npwx
             b(ii,l) = e(i) * spsi(ii,i) - hpsi(ii,i)               ! initial gradient and saved RHS for later
           END DO 
           !$acc loop vector
           DO ii = 1, npwx
             z(ii,l) = b(ii,l)
           END DO 
        end do
        !$acc end parallel

        CALL g_psi_ptr( npwx, npw, nnew, 1, z(1,nactive+1), e(nactive+done+1) )

     !- project on conduction bands
        CALL start_clock( 'pcg:ortho' )
        !$acc host_data use_device(spsi0, psi0, z, spsi0vec)
        CALL MYDGEMM( 'T','N', nbnd,nnew,npw2, 2.D0, spsi0, npwx2, z(:,nactive+1), npwx2, 0.D0, spsi0vec, nbnd )
        IF ( gstart == 2 ) CALL MYDGER( nbnd, nnew, -1.D0, spsi0, npwx2, z(:,nactive+1), npwx2, spsi0vec, nbnd )
        CALL mp_sum( spsi0vec, intra_bgrp_comm )
        CALL MYDGEMM( 'N','N', npw2,nnew,nbnd,-1.D0, psi0, npwx2, spsi0vec, nbnd, 1.D0, z(:,nactive+1), npwx2 )
        !$acc end host_data
        CALL stop_clock( 'pcg:ortho' )
     !-

        !$acc parallel loop 
        do l=nactive+1,nactive+nnew
           g0(l) = 2.D0*MYDDOT_VECTOR_GPU(npw2,z(:,l),b(:,l))
           IF (gstart==2) g0(l)=g0(l)-CONJG(z(1,l))*b(1,l)
        end do

        !$acc host_data use_device(g0)
        CALL mp_sum( g0(nactive+1:nactive+nnew), intra_bgrp_comm ) ! g0 = < initial z | initial gradient b >
        !$acc end host_data 

        !$acc parallel 
        !$acc loop gang private(i)
        do l=nactive+1,nactive+nnew; i=l+done
           !write(6,*) ' l =',l,' i =',i
           ff(l) = 0.d0 ; ff0(l) = ff(l)
           !write (6,*) 0, g0(l), ff(l)

           ! ethr_cg = ethr  ! CG convergence threshold could be set from input but it is not ...
           ethr_cg(l) = 1.0D-2  ! it makes more sense to fix the convergence of the CG solution to a 
                                ! fixed function of the RHS (see ethr_cg update later).
           ethr_cg(l) = max ( 0.01*ethr, ethr_cg(l) * g0(l) ) ! here we set the convergence of the correction
           !write(6,*) 'ethr_cg :', ethr_cg(l)

           !$acc loop vector 
           do ii = 1, npwx
             ! zero the trial solution
             psi(ii,i) = ZERO 
             hpsi(ii,i) = ZERO 
             spsi(ii,i) = ZERO
             ! initial search direction
             p(ii,l) = z(ii,l)
           end do 
           cg_iter(l) = 0 ! this is a new correction vector, reset its interation count
        end do
        !$acc end parallel

        nactive = nactive + nnew
     end if
     !write(6,*) ' done  =',done, ' nactive  =',nactive

!  iterate: !  DO cg_iter = 1, maxter  ! THIS IS THE ENTRY POINT OF THE PCG LOOP

     if ( nactive == 0 ) EXIT MAIN_LOOP               ! this is the only MAIN_LOOP EXIT condition

     !$acc kernels 
     cg_iter(1:nactive) = cg_iter(1:nactive) + 1      ! update interation counters
     !$acc end kernels

     CALL start_clock( 'pcg:hs_1psi' )
!     do l = 1, nactive  ! THIS COULD/SHOULD BE A GLOBAL CALL (ONLY WITHIN ONE BGRP THOUGH)
!        CALL hs_1psi( npwx, npw, p(:,l), hp(:,l), sp(:,l) ) ! apply H to a single wavefunction (no bgrp parallelization here!)
!     end do
     CALL hs_psi_ptr( npwx, npw, nactive, p, hp, sp ) ! apply H to a single wavefunction (no bgrp parallelization here!)
     CALL stop_clock( 'pcg:hs_1psi' )

     !$acc parallel loop private(i)     
     do l = 1, nactive; i=l+done
        gamma(l) = 2.D0*MYDDOT_VECTOR_GPU(npw2,p(:,l),hp(:,l)) &
                       - e(i) * 2.D0*MYDDOT_VECTOR_GPU(npw2,p(:,l),sp(:,l)) 
        IF (gstart==2) gamma(l) = gamma(l) - CONJG(p(1,l))*hp(1,l) + e(i) * CONJG(p(1,l))*sp(1,l)
     end do

     !$acc host_data use_device(gamma)
     CALL mp_sum( gamma(1:nactive), intra_bgrp_comm ) ! gamma = < p | hp - e sp >
     !$acc end host_data

     !$acc parallel 
     !$acc loop gang private(i)
     do l = 1, nactive; i=l+done
        !write(6,*) ' l =',l,' i =',i

        alpha(l) = g0(l)/gamma(l)
        !write(6,*) 'g0, gamma, alpha :', g0(l), gamma(l), alpha(l)
        !$acc loop vector 
        DO ii = 1, npwx
          psi(ii,i)  = psi(ii,i)  + alpha(l) * p(ii,l)     ! updated solution
          hpsi(ii,i) = hpsi(ii,i) + alpha(l) * hp(ii,l)    ! updated solution
          spsi(ii,i) = spsi(ii,i) + alpha(l) * sp(ii,l)    ! updated solution
        END DO 

        g2(l) = 2.D0 * ( MYDDOT_VECTOR_GPU(npw2,z(:,l),b(:,l)) &
                + e(i) * MYDDOT_VECTOR_GPU(npw2,z(:,l),spsi(:,i)) &
                - MYDDOT_VECTOR_GPU(npw2,z(:,l),hpsi(:,i)) )
        IF (gstart==2) g2(l) = g2(l) - CONJG(z(1,l))*b(1,l) - e(i)*CONJG(z(1,l))*spsi(1,i) + CONJG(z(1,l))*hpsi(1,i)
     end do
     !$acc end parallel

     !$acc host_data use_device(g2)
     CALL mp_sum( g2(1:nactive), intra_bgrp_comm )    ! g2 = < old z | new gradient b + e spsi - hpsi >
     !$acc end host_data

     !$acc parallel loop collapse(2) private(i)
     do l = 1, nactive                                ! update the preconditioned gradient
        DO ii = 1, npwx
          i=l+done
          z(ii,l) = b(ii,l) + e(i) * spsi(ii,i) - hpsi(ii,i) 
        END DO 
     end do

     CALL g_psi_ptr( npwx, npw, nactive, 1, z, e(done+1) )

  !- project on conduction bands
     CALL start_clock( 'pcg:ortho' )
     !$acc host_data use_device(spsi0, psi0, z, spsi0vec)
     CALL MYDGEMM( 'T','N', nbnd,nactive,npw2, 2.D0, spsi0, npwx2, z, npwx2, 0.D0, spsi0vec, nbnd )
     IF ( gstart == 2 ) CALL MYDGER( nbnd, nactive, -1.D0, spsi0, npwx2, z, npwx2, spsi0vec, nbnd )
     CALL mp_sum( spsi0vec, intra_bgrp_comm )
     CALL MYDGEMM( 'N','N', npw2,nactive,nbnd,-1.D0, psi0, npwx2, spsi0vec, nbnd, 1.D0, z, npwx2 )
     !$acc end host_data
     CALL stop_clock( 'pcg:ortho' )
  !-
     !$acc parallel loop private(i) 
     do l = 1, nactive; i=l+done
        g1(l) = 2.D0 * ( MYDDOT_VECTOR_GPU(npw2,z(:,l),b(:,l)) &
                + e(i) * MYDDOT_VECTOR_GPU(npw2,z(:,l),spsi(:,i)) &
                       - MYDDOT_VECTOR_GPU(npw2,z(:,l),hpsi(:,i)) )
        IF (gstart==2) g1(l) = g1(l) - CONJG(z(1,l)) * ( b(1,l) + e(i) * spsi(1,i) - hpsi(1,i) )
     end do

     !$acc host_data use_device(g1)
     CALL mp_sum( g1(1:nactive), intra_bgrp_comm )   ! g1 = < new z | new gradient b + e spsi - hpsi >
     !$acc end host_data 

     !$acc parallel loop private(i) 
     do l = 1, nactive; i = l + done                 ! evaluate the function ff
        ff(l) = - ( e(i)*MYDDOT_VECTOR_GPU(npw2,psi(:,i),spsi(:,i)) &
                        -MYDDOT_VECTOR_GPU(npw2,psi(:,i),hpsi(:,i)) ) &
                - 2.D0 * MYDDOT_VECTOR_GPU(npw2,psi(:,i),b(:,l))
        if (gstart==2) ff(l) = ff(l) + 0.5D0 * CONJG(psi(1,i))*( e(i)*spsi(1,i) - hpsi(1,i) + 2.D0 * b(1,l) )
     end do

     !$acc host_data use_device(ff)
     CALL mp_sum( ff(1:nactive), intra_bgrp_comm )   ! function minimum -0.5 < psi | e spsi - hpsi > - < psi | b >
     !$acc end host_data

     newdone = 0  ! number of correction vectors that converge (or are done) at this iteration
     CALL start_clock( 'pcg:move1' )
     do l = 1, nactive; i = l + done
        !write (6,*) cg_iter(l), g1(l), ff(l),  gamma(l)

        !$acc serial copyout(ff_check, ff0_check, g1_check, cg_iter_l, iter_check) copyin(maxter)
        ff_check   = ff(l) > ff0(l)  
        ff0_check  = ff0(l) < 0.d0 
        g1_check   = ABS ( g1(l) ) < ethr_cg(l)  
        cg_iter_l  = cg_iter(l)
        iter_check = cg_iter_l == maxter 
        !$acc end serial 

        !write (6,*) cg_iter(l), g1(l), ff(l),  gamma(l)

        IF ( ff_check .AND. ff0_check ) THEN
           !$acc parallel loop 
           DO ii = 1, npwx
             psi(ii,i)  = psi(ii,i)  - alpha(l) * p(ii,l) ! fallback solution: if last iter failed to improve ff0
             hpsi(ii,i) = hpsi(ii,i) - alpha(l) * hp(ii,l)! exit whitout updating and ...
             spsi(ii,i) = spsi(ii,i) - alpha(l) * sp(ii,l)! hope next time it'll be better
           END DO 
        END IF

        !write(6,*) 'g0, g1, g2 :', g0(l), g1(l), g2(l)
        !write(6,*) 'ff0, ff : ', ff0(l), ff(l)
        IF ( g1_check .OR. ff_check .OR. iter_check ) THEN ! EXIT iterate
           !write (6,*) ' exit  pcg loop'
           !write(6,*) ' l =',l,' i =',i
           !if ( cg_iter(l) == maxter.and. ABS(g1(l)) > ethr_cg(l))  write (6,*) 'CG not converged maxter exceeded', cg_iter(l), g1(l), g0(l), ethr_cg(l)
           !IF ( ABS ( g1(l) ) < ethr_cg(l)) write (6,*) 'CG correction converged ', cg_iter(l), g1(l), ethr_cg(l)
           !IF ( ABS ( g1(l) ) > g0(l)     ) write (6,*) 'CG not converged ', cg_iter(l), g1(l), g0(l), ethr_cg(l)

           nhpsi = nhpsi + cg_iter_l    ! update nhpsi count
           IF (.NOT. ( g1_check .OR. ff_check ) .AND. iter_check ) nhpsi = nhpsi + 1 ! because this would be the count

           newdone = newdone + 1         ! one more solution found (or no more active anyway)

           !write(6,*) ' newdone = ', newdone

           !write(6,*) ' swapping converged psi/hpsi/spsi i = ',i, " with i' = ",done+newdone
           ! swap the terminated vector with the first in the list of the active ones
           !$acc kernels copyin(done, newdone)
           !$acc loop independent
           DO ii = 1, npwx
             p (ii,l) = psi (ii,done+newdone) ; psi (ii,done+newdone) = psi (ii,i) ; psi (ii,i) = p (ii,l)
             hp(ii,l) = hpsi(ii,done+newdone) ; hpsi(ii,done+newdone) = hpsi(ii,i) ; hpsi(ii,i) = hp(ii,l)
             sp(ii,l) = spsi(ii,done+newdone) ; spsi(ii,done+newdone) = spsi(ii,i) ; spsi(ii,i) = sp(ii,l)
           END DO 

           ee = e(done+newdone)       
           e(done+newdone) = e(i)      
           e(i) = ee

           !write(6,*) ' overwrite converged p/hp/etc l = ',l, ' with newdone = ',newdone
           ! move information of the swapped active vector in the right place to keep going
           !$acc loop independent
           DO ii = 1, npwx
             p(ii,l) = p(ii,newdone) ; hp(ii,l) = p(ii,newdone)  ; sp(ii,l) = sp(ii,newdone)
             b(ii,l) = b(ii,newdone) ;  z(ii,l) = z(ii,newdone)
           END DO

           ff0(l) = ff0(newdone) ; ff(l) = ff(newdone)
           alpha(l) = alpha(newdone) ; g0(l) = g0(newdone) ; g1(l) = g1(newdone) ; g2(l) = g2(newdone)
           cg_iter(l) = cg_iter(newdone) ; ethr_cg(l) = ethr_cg(newdone)
           !$acc end kernels

        ELSE

           !write(6,*) ' l =',l,' i =',i
           !$acc kernels
           beta   = (g1(l)-g2(l))/g0(l)         ! Polak - Ribiere style update
           g0(l)  = g1(l)                       ! < new z | new gradient >  ->  < old z | old gradient >
           !$acc loop independent
           DO ii = 1, npwx
             p(ii,l) = z(ii,l) + beta * p(ii,l)      ! updated search direction
           END DO 
           !write(6,*) 'beta :', beta

           ff0(l) = ff(l)                       ! updated minimum value reached by the function
           !$acc end kernels

        END IF
     end do
     CALL stop_clock( 'pcg:move1' )

     IF ( newdone > 0 ) THEN

        done    = done    + newdone
        nactive = nactive - newdone

        !write(6,*) ' there have been ', newdone, ' new converged solution'
        !write(6,*) ' done = ', done, ' nactive =', nactive

        CALL start_clock( 'pcg:move2' )
        do l=1, nactive

           !write(6,*) ' l+newdone =',l+newdone,'  ->   l =',l
           !$acc parallel loop 
           DO ii = 1, npwx
             p(ii,l) = p(ii,l+newdone) ; hp(ii,l) = hp(ii,l+newdone) ; sp(ii,l) = sp(ii,l+newdone)
             b(ii,l) = b(ii,l+newdone) ;  z(ii,l) =  z(ii,l+newdone) 
           END DO 
           !$acc kernels
           ff0(l) = ff0(l+newdone) ; ff(l) = ff(l+newdone)
           g0(l) = g0(l+newdone) ; g1(l) = g1(l+newdone) ; g2(l) = g2(l+newdone)
           cg_iter(l) = cg_iter(l+newdone) ; ethr_cg(l) = ethr_cg(l+newdone)
           !$acc end kernels
        end do
        CALL stop_clock( 'pcg:move2' )

     END IF

!  END DO iterate      Here is where the pcg loop would terminate
    
  END DO  MAIN_LOOP
  !write (6,*) ' exit  pcg loop'

  !$acc exit data delete(p, hp, sp, z, b)
  DEALLOCATE( spsi0vec )
  DEALLOCATE( b, p, hp, sp, z )
  DEALLOCATE( ethr_cg, ff, ff0, cg_iter )
  DEALLOCATE( g0, g1, g2, alpha, gamma )
  !
  CALL stop_clock( 'pcg' )
  !
  RETURN
  !
END SUBROUTINE bpcg_gamma
