!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#define DIIS_DEBUG
#undef  DIIS_DEBUG
!
!#define CG_STEP
#undef CG_STEP
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
MODULE diis_module
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  !     ( H - e S ) |psi> = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an overlap matrix, |psi> is a complex vector.  
  !
  ! ... a variant of the DIIS Residual Minimization Method is used
  !
  ! ... written by Carlo Sbraccia ( 25/04/2004 )
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  ! ... public methods
  !
  PUBLIC :: cdiisg, rdiisg
  !
  ! ... internal DIIS variables
  !
  INTEGER        :: maxter      = 20       ! maximum number of iterations
  INTEGER        :: sd_maxter_0 = 2        ! initial number of trial iterations 
  REAL (KIND=DP) :: lambda_0    = 1.0D0    ! initial trial steps
  !
  ! ... module procedures :
  !  
  CONTAINS
    !
    ! ... real routines
    !
    !------------------------------------------------------------------------
    SUBROUTINE rdiisg( ndim, ndmx, nbands, diis_ndim, psi, &
                       e, ethr, btype, notcnv, diis_iter, iter )
      !------------------------------------------------------------------------
      !
      ! ... gamma version of the diis algorithm :
      ! ...    real wavefunctions with only half plane waves stored  
      !
      USE gvect, ONLY : gstart
      !
      IMPLICIT NONE
      !
      ! ... on INPUT
      !
      INTEGER :: ndim, ndmx, nbands, diis_ndim, btype(nbands), iter
        ! dimension of the matrix to be diagonalized
        ! leading dimension of matrix psi, as declared in the calling pgm unit
        ! integer number of searched low-lying roots
        ! maximum dimension of the reduced basis set :
        !   the basis set is refreshed when its dimension would exceed diis_ndim
        ! band type ( 1 = occupied, 0 = empty )
        ! scf iteration
      REAL (KIND=DP) :: ethr
        ! energy threshold for convergence :
        !   eigenvector improvement is stopped, when two consecutive estimates 
        !   of the eigenvalue differ by less than ethr.
      !
      ! ... on OUTPUT
      !
      COMPLEX (KIND=DP) :: psi(ndmx,nbands)
        !  psi contains the refined estimates of the eigenvectors
      REAL (KIND=DP) :: e(nbands)
        ! contains the estimated eigenvalues
      INTEGER :: diis_iter  
        ! average number of iterations performed per band
      INTEGER :: notcnv
        ! number of unconverged bands
      !
      ! ... LOCAL variables
      !
      INTEGER :: kter, ib, cnv, n, m
        ! counter on iterations
        ! counter on bands
        ! number of converged bands
        ! do-loop counters
      REAL :: empty_ethr
        ! threshold for empty bands    
      REAL (KIND=DP), ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:)
        ! H matrix on the reduced basis
        ! S matrix on the reduced basis
        ! the eigenvectors of the Hamiltonian
      COMPLEX (KIND=DP), ALLOCATABLE :: hpsi(:,:), spsi(:,:), aux(:,:)
        ! the product of H and psi
        ! the product of S and psi
        ! work space, contains the residual vector
      COMPLEX (KIND=DP), ALLOCATABLE :: psi_old(:,:,:), hpsi_old(:,:,:), &
                                        spsi_old(:,:,:)
        ! diis-workspace: old eigenvectors
        ! diis-workspace: old product of H and psi
        ! diis-workspace: old product of S and psi
      REAL (KIND=DP), ALLOCATABLE :: e_old(:,:)
        ! eigenvalues of the previous iteration
      REAL (KIND=DP), ALLOCATABLE :: lambda(:)
        ! trial step length 
      INTEGER, ALLOCATABLE :: nbase(:), sd_maxter(:)
        ! counter on the reduced basis vectors
        ! maximum number of trial steps (different for each band)
      LOGICAL, ALLOCATABLE :: conv(:), frozen(:)
        ! .TRUE. if the band is converged
        ! .TRUE. if the band was converged at the previous step
      !
      !
      CALL start_clock( 'diis' )
      !
      ! ... work-space allocation
      !
      ALLOCATE( psi_old(  ndmx, diis_ndim , nbands ) )
      ALLOCATE( hpsi_old( ndmx, diis_ndim , nbands ) )
      ALLOCATE( spsi_old( ndmx, diis_ndim , nbands ) )
      ALLOCATE( e_old( diis_ndim, nbands ) )
      ALLOCATE( hpsi( ndmx, nbands ) )
      ALLOCATE( spsi( ndmx, nbands ) )
      ALLOCATE( aux(  ndmx, nbands ) )
      ALLOCATE( hr(   nbands, nbands ) )
      ALLOCATE( sr(   nbands, nbands ) )
      ALLOCATE( vr(   nbands, nbands ) )
      ALLOCATE( lambda( nbands ) )
      ALLOCATE( nbase(     nbands ) )
      ALLOCATE( sd_maxter( nbands ) )
      ALLOCATE( conv(      nbands ) )
      ALLOCATE( frozen(    nbands ) )
      !
      ! ... initialization of control variables
      !
      kter      = 0
      conv      = .FALSE.
      cnv       = 0
      notcnv    = nbands
      nbase     = 0 
      lambda    = lambda_0
      sd_maxter = sd_maxter_0
      !
      ! ... threshold for empty bands
      !
      empty_ethr = MAX( ( ethr * 50.D0 ), 1.D-3 )
      !
      IF ( iter > 1 ) THEN
         !
         ! ... after the first scf iteration one SD step is supposed to
         ! ... be sufficient. 
         !
         sd_maxter = 1
         !
      END IF
      !
      ! ... initialization of the work-space
      !
      e_old      = 0.D0
      e_old(1,:) = e(:)
      psi_old    = ZERO
      hpsi_old   = ZERO
      spsi_old   = ZERO
      !
      ! ... diagonalization loop
      !
      iterate: DO
         !
         kter = kter + 1
         !
         IF ( kter > maxter ) EXIT iterate
         !
#if defined (DIIS_DEBUG)
         PRINT *, ""
         PRINT *, "kter = ", kter
         PRINT *, ""
#endif
         !
         ! ... the size of the history-subspace is increased ( different for 
         ! ... for each band )
         !
         WHERE ( .NOT. conv(:) )
            !
            nbase(:) = MIN( ( nbase(:) + 1 ), diis_ndim )
            !
         END WHERE        
         !
         ! ... bands are reordered so that converged bands come first
         !
         CALL reorder_bands()
         !
         ! ... here we compute H|psi> and S|psi> for not converged bands
         !
         CALL h_psi( ndmx, ndim, notcnv, psi(:,cnv+1), hpsi(:,cnv+1) )
         CALL s_psi( ndmx, ndim, notcnv, psi(:,cnv+1), spsi(:,cnv+1) )
         !
         ! ... here we set up the hamiltonian and overlap matrix 
         ! ... on the subspace :
         !
         ! ...    hc_ij = <psi_i|H|psi_j>  and  sc_ij = <psi_i|S|psi_j>
         !
         CALL DGEMM( 'T', 'N', nbands, nbands, 2*ndim, 2.D0, &
                     psi, 2*ndmx, hpsi, 2*ndmx, 0.D0, hr, nbands )
         CALL DGEMM( 'T', 'N', nbands, nbands, 2*ndim, 2.D0, &
                     psi, 2*ndmx, spsi, 2*ndmx, 0.D0, sr, nbands )
         !
         IF ( gstart == 2 ) THEN
            !
            CALL DGER( nbands, nbands, -1.D0, psi, &
                       2*ndmx, hpsi, 2*ndmx, hr, nbands )
            CALL DGER( nbands, nbands, -1.D0, psi, &
                       2*ndmx, spsi, 2*ndmx, sr, nbands )
            !
         END IF              
         !
         CALL reduce( nbands * nbands, hr )
         CALL reduce( nbands * nbands, sr )
         !
         ! ... the subspace rotation is performed by solving the eigenvalue 
         ! ... problem on the subspace :
         !     
         ! ...    ( hr - e * sr ) |phi> = 0
         !
         CALL rdiaghg( nbands, nbands, hr, sr, nbands, e, vr )
         !
         ! ... convergence is checked here
         !
         frozen = conv
         !
         WHERE( btype(:) == 1 )
            !
            conv(:) = ( ( kter > 1 ) .AND. &
                        ( ABS( e(:) - e_old(1,:) ) < ethr ) )
            !
         ELSEWHERE
            !
            conv(:) = ( ( kter > 1 ) .AND. &
                        ( ABS( e(:) - e_old(1,:) ) < empty_ethr ) )
            !
         END WHERE
         !
         ! ... |psi>,  H|psi>  and  S|psi>  are rotated
         !
         CALL DGEMM( 'N', 'N', 2*ndim, nbands, nbands, 1.D0, &
                     psi, 2*ndmx, vr, nbands, 0.D0, aux, 2*ndmx ) 
         !     
         psi(:,:) = aux(:,:)
         !
         CALL DGEMM( 'N', 'N', 2*ndim, nbands, nbands, 1.D0, &
                     hpsi, 2*ndmx, vr, nbands, 0.D0, aux, 2*ndmx ) 
         !     
         hpsi(:,:) = aux(:,:)
         !
         CALL DGEMM( 'N', 'N', 2*ndim, nbands, nbands, 1.D0, &
                     spsi, 2*ndmx, vr, nbands, 0.D0, aux, 2*ndmx ) 
         !     
         spsi(:,:) = aux(:,:)
         !
#if defined (DIIS_DEBUG)
         PRINT *, "eigenvalues :"
         PRINT *, e(:)
         PRINT *, "variation :"
         PRINT *, ABS( e(:) - e_old(1,:) )
         PRINT *, conv(:)
#endif
         !
         ! ... exit if all bands are converged
         !
         IF ( ALL( conv(:) ) ) EXIT iterate
         !
         ! ... a new step is performed
         !
         bands_loop: DO ib = 1, nbands
            !
            ! ... the two highest bands are always optimized 
            ! ... ( even if converged )
            !
            IF ( conv(ib) .AND. ( ib <= ( nbands - 2 ) ) ) CYCLE bands_loop
            !
            ! ... holes-sniffer
            !
            IF ( ( ib <= ( nbands - 2 ) ) .AND.  &
                 frozen(ib) .AND. ( kter > sd_maxter(ib) ) ) THEN
               !
               ! ... an hole has been detected :  the work-space is cleaned
               !
#if defined (DIIS_DEBUG)               
               PRINT *, "HOLES-SNIFFER: ", &
                        "an hole has been detected near band ", ib
#endif               
               !
               nbase(ib) = 1
               !
               e_old(:,ib)      = 0.D0
               psi_old(:,:,ib)  = ZERO
               hpsi_old(:,:,ib) = ZERO
               spsi_old(:,:,ib) = ZERO
               !
               sd_maxter(ib) = kter + sd_maxter_0               
               !
               maxter = MAX( maxter, sd_maxter(ib) )
               !
            END IF
            !
            ! ... the two highest bands are never optimized with the 
            ! ... DIIS procedure
            !
            IF ( ( kter <= sd_maxter(ib) ) .OR. ( ib > ( nbands - 2 ) ) ) THEN
               !
               ! ... trial step :
               !     
               IF ( ( e(ib) < e_old(1,ib) ) .OR. ( kter == 1 ) ) THEN
                  !
                  IF ( kter == 1 ) THEN
                     !
                     ! ... first step :  lambda is kept fixed
                     !
                  ELSE
                     !
                     ! ... step accepted :  lambda is increased
                     !
                     lambda(ib) = lambda(ib) * 1.25D0
                     !
                  END IF   
                  !
               ELSE
                  !
                  ! ... step rejected :  lambda is decreased
                  !
                  lambda(ib) = 0.75D0 * lambda(ib)
                  !
               END IF
               !
               ! ... here we compute the residual vector ( stored in aux ) :
               !
               ! ...    |res> = ( H - e * S ) |psi>           
               !
               aux(:,ib) = hpsi(:,ib) - e(ib) * spsi(:,ib)
               !
            ELSE   
               !
               ! ... DIIS step :  the best eigenvector and residual vector 
               ! ...              are computed
               !
               CALL diis_step( ib )           
               !
            END IF
            !
            ! ... DIIS work-space is refreshed ( for not converged bands only )
            !
            DO n = ( diis_ndim - 1 ), 2, -1
               !
               e_old(n,ib)      = e_old(n-1,ib)
               psi_old(:,n,ib)  = psi_old(:,n-1,ib)
               hpsi_old(:,n,ib) = hpsi_old(:,n-1,ib)                 
               spsi_old(:,n,ib) = spsi_old(:,n-1,ib)
               !
            END DO
            !
            e_old(1,ib)      = e(ib)
            psi_old(:,1,ib)  = psi(:,ib)   
            hpsi_old(:,1,ib) = hpsi(:,ib)
            spsi_old(:,1,ib) = spsi(:,ib)
            !
         END DO bands_loop
         !
#if defined (CG_STEP)
         !
         CALL cg_step()
         !
#else         
         !
         ! ... preconditioning of all residual vecotrs
         !
         CALL g_psi( ndmx, ndim, nbands, aux, e )           
         !
#endif         
         !           
         ! ... here we compute the new eigenvectors for non converged bands 
         ! ... only 
         ! ... ( the highest two bands are always optimized, even if converged )
         !
         FORALL( ib = 1: nbands, .NOT. conv(ib) .OR. ib > ( nbands - 2 ) )
            !
            psi(:,ib) = psi(:,ib) - lambda(ib) * aux(:,ib)
            !
         END FORALL   
         !
      END DO iterate
      !
      ! ... this is an overestimate of the real number of iterations
      !
      diis_iter = MIN( kter, maxter )
      !
      ! ... the number of not-converged bands is computed
      !
      notcnv = nbands
      !
      DO ib = 1, nbands
         !
         IF ( conv(ib) ) notcnv = notcnv - 1
         !
      END DO
      !
      ! ... work-space deallocation
      !
      DEALLOCATE( frozen )
      DEALLOCATE( conv )
      DEALLOCATE( sd_maxter )
      DEALLOCATE( nbase )
      DEALLOCATE( lambda )
      DEALLOCATE( vr )
      DEALLOCATE( sr )
      DEALLOCATE( hr )
      DEALLOCATE( aux )
      DEALLOCATE( hpsi )
      DEALLOCATE( spsi )
      DEALLOCATE( e_old )
      DEALLOCATE( psi_old )
      DEALLOCATE( hpsi_old )
      DEALLOCATE( spsi_old )
      !
      CALL stop_clock( 'diis' )
      !
      RETURN
      !
      CONTAINS
        !
        ! ... internal procedures
        !
        !--------------------------------------------------------------------
        SUBROUTINE reorder_bands()
          !--------------------------------------------------------------------
          !
          ! ... this routine is used to reorder the bands :
          ! ... converged bands come first.
          ! ... for this pourpose an auxiliary vector is used
          !
          IMPLICIT NONE
          !
          !
          cnv    = 0
          notcnv = nbands
          !
          IF ( ( kter <= 2 ) .OR. ( .NOT. ANY( conv(:) ) ) ) RETURN
          !
          DO ib = 1, nbands
             !
             IF ( conv(ib) .AND. ( ib <= ( nbands - 2 ) ) ) THEN
                !
                cnv = cnv + 1
                !
                aux(:,cnv) = psi(:,ib)
                !
                hpsi(:,cnv) = hpsi(:,ib)
                spsi(:,cnv) = spsi(:,ib)              
                !
             ELSE
                !
                aux(:,notcnv) = psi(:,ib)
                !
                notcnv = notcnv - 1
                !
             END IF
             !
          END DO
          !
          notcnv = nbands - notcnv
          !
          psi(:,:) = aux(:,:)
          !
          RETURN
          !
        END SUBROUTINE reorder_bands
        !      
#if defined (CG_STEP)
        !
        !--------------------------------------------------------------------
        SUBROUTINE cg_step()
          !--------------------------------------------------------------------
          !
          IMPLICIT NONE
          !
          REAL (KIND=DP)              :: gamma, num, den
          REAL (KIND=DP), ALLOCATABLE :: res_old(:), pres(:), pres_old(:)
          !
          REAL (KIND=DP), EXTERNAL :: DDOT      
          !
          !
          IF ( kter >= 2 ) THEN
             !
             ! ... preconditioned conjugate gradients step
             !
             ALLOCATE( res_old(  ndmx ) )
             ALLOCATE( pres(     ndmx ) )
             ALLOCATE( pres_old( ndmx ) )
             !
             bands_loop: DO ib = 1, nbands
                !
                IF ( conv(ib) ) CYCLE bands_loop
                !
                pres    = ( hpsi_old(:,1,ib) - e_old(1,ib) * spsi_old(:,1,ib) )
                res_old = ( hpsi_old(:,2,ib) - e_old(2,ib) * spsi_old(:,2,ib) )
                !
                pres_old = res_old
                !
                CALL g_psi( ndmx, ndim, 1, pres,     e_old(1,ib) )
                CALL g_psi( ndmx, ndim, 1, pres_old, e_old(2,ib) )
                !
                num = 2.D0 * DDOT( 2*ndim, pres, 1, aux(:,ib), 1 )
                !
                IF ( gstart == 2 ) num = num - pres(1) * aux(1,ib)
                !
                den = 2.D0 * DDOT( 2*ndim, pres_old, 1, res_old, 1 )
                !
                IF ( gstart == 2 ) den = den - pres_old(1) * res_old(1)
                !
                gamma = num / den
                !
                aux(:,ib) = pres - gamma * ( psi_old(:,1,ib) - psi_old(:,2,ib) )
                !
             END DO bands_loop 
             !
             DEALLOCATE( res_old )
             DEALLOCATE( pres )
             DEALLOCATE( pres_old )             
             !
          ELSE
             !
             ! ... standard preconditioned steepest descent step
             !
             CALL g_psi( ndmx, ndim, nbands, aux, e )
             !
          END IF
          !
          RETURN
          !
        END SUBROUTINE cg_step
        !
#endif
        !          
        !--------------------------------------------------------------------
        SUBROUTINE diis_step( ib )
          !--------------------------------------------------------------------
          !
          IMPLICIT NONE
          !
          INTEGER, INTENT(IN)            :: ib
          INTEGER                        :: dim, n
          REAL (KIND=DP)                 :: psiSpsi
          REAL (KIND=DP), ALLOCATABLE    :: e_small(:)
          REAL (KIND=DP), ALLOCATABLE    :: rr_small(:,:), &
                                            sr_small(:,:), &
                                            vr_small(:,:)
          COMPLEX (KIND=DP), ALLOCATABLE :: all_psi(:,:),  &
                                            all_hpsi(:,:), &
                                            all_spsi(:,:), &
                                            all_res(:,:)
          !
          REAL (KIND=DP), EXTERNAL :: DDOT
          !
          !
          dim = nbase(ib)
          !
          ! ... internal work-space allocation
          !
          ALLOCATE( e_small( dim ) )
          ALLOCATE( rr_small( dim, dim ) )
          ALLOCATE( sr_small( dim, dim ) )
          ALLOCATE( vr_small( dim, dim ) )
          ALLOCATE( all_psi(  ndmx, dim ) )
          ALLOCATE( all_hpsi( ndmx, dim ) )
          ALLOCATE( all_spsi( ndmx, dim ) )
          ALLOCATE( all_res(  ndmx, dim ) )
          !
          ! ... the history of this band is recostructed
          !
          !
          ! ... the history of this band is reconstructed
          !
          all_psi(:,1)  = psi(:,ib)
          all_hpsi(:,1) = hpsi(:,ib)
          all_spsi(:,1) = spsi(:,ib)
          e_small(1)    = e(ib)
          !
          all_psi(:,2:dim)  = psi_old(:,:,ib)
          all_hpsi(:,2:dim) = hpsi_old(:,:,ib)
          all_spsi(:,2:dim) = spsi_old(:,:,ib)
          e_small(2:dim)    = e_old(:,ib)
          !
          ! ... orthogonalization
          !
          CALL cgramg1( ndmx, dim, ndim, 1, dim, all_psi, all_spsi, all_hpsi )
          !
          FORALL( n = 1: dim ) &
             all_res(:,n) = ( all_hpsi(:,n) - e_small(n) * all_spsi(:,n) )
          !   
          ! ... here we construct the matrices :
          ! ...    rr_ij = <res_i|res_j>  and  sr_ij = <psi_i|S|psi_j>
          !
          CALL DGEMM( 'T', 'N', dim, dim, 2*ndim, 2.D0, all_res(:,:), &
                      2*ndmx, all_res(:,:), 2*ndmx, 0.D0, rr_small, dim )
          !
          IF ( gstart == 2 ) &
             CALL DGER( dim, dim, -1.D0, all_res(:,:), &
                        2*ndmx, all_res(:,:), 2*ndmx, rr_small, dim )
          !
          CALL reduce( dim * dim, rr_small )
          !
          sr_small(n,n) = 0.D0
          !
          FORALL( n = 1: dim ) sr_small(n,n) = 1.D0
          ! 
          ! ... diagonalize the reduced hamiltonian
          !
          CALL rdiaghg( dim, 1, rr_small, sr_small, dim, e_small, vr_small )
          !
          ! ... here we compute the best estimate of the |psi>, H|psi>, S|psi>
          !
          CALL DGEMM( 'N', 'N', 2*ndim, 1, dim, 1.D0, all_psi(:,:), &
                      2*ndmx, vr_small(:,1), dim, 0.D0, psi(:,ib), 2*ndmx )
          !
          CALL DGEMM( 'N', 'N', 2*ndim, 1, dim, 1.D0, all_hpsi(:,:), &
                      2*ndmx, vr_small(:,1), dim, 0.D0, hpsi(:,ib), 2*ndmx )
          !
          CALL DGEMM( 'N', 'N', 2*ndim, 1, dim, 1.D0, all_spsi(:,:), &
                      2*ndmx, vr_small(:,1), dim, 0.D0, spsi(:,ib), 2*ndmx )
          !
          psiSpsi = 2.D0 * DDOT( 2*ndim, psi(:,ib), 1, spsi(:,ib), 1 )
          !
          IF ( gstart == 2 ) psiSpsi = psiSpsi - psi(1,ib) * spsi(1,ib)
          !
          CALL reduce( 1, psiSpsi )
          !          
          e(ib) = 2.D0 * DDOT( 2*ndim, psi(:,ib), 1, hpsi(:,ib), 1 )
          !
          IF ( gstart == 2 ) e(ib) = e(ib) - psi(1,ib) * hpsi(1,ib)
          !
          CALL reduce( 1, e(ib) )          
          !
          e(ib) = e(ib) / psiSpsi
          !
          ! ... here we compute the best estimate of the residual vector
          !
          aux(:,ib) = ( hpsi(:,ib) - e(ib) * spsi(:,ib) ) / psiSpsi
          !
          ! ... internal work-space deallocation
          !
          DEALLOCATE( e_small )
          DEALLOCATE( rr_small )
          DEALLOCATE( sr_small )
          DEALLOCATE( vr_small )
          DEALLOCATE( all_psi )
          DEALLOCATE( all_hpsi )
          DEALLOCATE( all_spsi )
          DEALLOCATE( all_res )
          !   
          RETURN
          !
        END SUBROUTINE diis_step
        !
    END SUBROUTINE rdiisg    
    !
    ! ... complex routines
    !
    !------------------------------------------------------------------------
    SUBROUTINE cdiisg( ndim, ndmx, nbands, diis_ndim, psi, &
                       e, ethr, btype, notcnv, diis_iter, iter )
      !------------------------------------------------------------------------
      !
      ! ... k-points version of the diis algorithm
      !
      IMPLICIT NONE
      !
      ! ... on INPUT
      !
      INTEGER :: ndim, ndmx, nbands, diis_ndim, btype(nbands), iter
        ! dimension of the matrix to be diagonalized
        ! leading dimension of matrix psi, as declared in the calling pgm unit
        ! integer number of searched low-lying roots
        ! maximum dimension of the reduced basis set :
        !   the basis set is refreshed when its dimension would exceed diis_ndim
        ! band type ( 1 = occupied, 0 = empty )
        ! scf iteration
      REAL (KIND=DP) :: ethr
        ! energy threshold for convergence :
        !   eigenvector improvement is stopped, when two consecutive estimates 
        !   of the eigenvalue differ by less than ethr.
      !
      ! ... on OUTPUT
      !
      COMPLEX (KIND=DP) :: psi(ndmx,nbands)
        !  psi contains the refined estimates of the eigenvectors
      REAL (KIND=DP) :: e(nbands)
        ! contains the estimated eigenvalues
      INTEGER :: diis_iter  
        ! average number of iterations performed per band
      INTEGER :: notcnv
        ! number of unconverged bands
      !
      ! ... LOCAL variables
      !
      INTEGER :: kter, ib, cnv, n, m
        ! counter on iterations
        ! counter on bands
        ! number of converged bands
        ! do-loop counters
      REAL :: empty_ethr
        ! threshold for empty bands  
      COMPLEX (KIND=DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
        ! H matrix on the reduced basis
        ! S matrix on the reduced basis
        ! the eigenvectors of the Hamiltonian
      COMPLEX (KIND=DP), ALLOCATABLE :: hpsi(:,:), spsi(:,:), aux(:,:)
        ! the product of H and psi
        ! the product of S and psi
        ! work space, contains the residual vector
      COMPLEX (KIND=DP), ALLOCATABLE :: psi_old(:,:,:), hpsi_old(:,:,:), &
                                        spsi_old(:,:,:)
        ! diis-workspace: old eigenvectors
        ! diis-workspace: old product of H and psi
        ! diis-workspace: old product of S and psi
      REAL (KIND=DP), ALLOCATABLE :: e_old(:,:)
        ! eigenvalues of the previous iteration
      REAL (KIND=DP), ALLOCATABLE :: lambda(:)      
        ! trial step length 
      INTEGER, ALLOCATABLE :: nbase(:), sd_maxter(:)
        ! counter on the reduced basis vectors
        ! maximum number of trial steps (different for each band)
      LOGICAL, ALLOCATABLE :: conv(:), frozen(:)
        ! .TRUE. if the band is converged
        ! .TRUE. if the band was converged at the previous step
      !
      !
      CALL start_clock( 'diis' )
      !
      ! ... work-space allocation
      !
      ALLOCATE( psi_old(  ndmx, ( diis_ndim - 1 ) , nbands ) )
      ALLOCATE( hpsi_old( ndmx, ( diis_ndim - 1 ) , nbands ) )
      ALLOCATE( spsi_old( ndmx, ( diis_ndim - 1 ) , nbands ) )
      ALLOCATE( e_old( ( diis_ndim - 1 ), nbands ) )
      ALLOCATE( hpsi( ndmx, nbands ) )
      ALLOCATE( spsi( ndmx, nbands ) )
      ALLOCATE( aux(  ndmx, nbands ) )
      ALLOCATE( hc(   nbands, nbands ) )
      ALLOCATE( sc(   nbands, nbands ) )
      ALLOCATE( vc(   nbands, nbands ) )
      ALLOCATE( lambda( nbands ) )
      ALLOCATE( nbase(     nbands ) )
      ALLOCATE( sd_maxter( nbands ) )
      ALLOCATE( conv(      nbands ) )
      ALLOCATE( frozen(    nbands ) )
      !
      ! ... initialization of control variables
      !
      kter      = 0
      conv      = .FALSE.
      cnv       = 0
      notcnv    = nbands
      nbase     = 0 
      lambda    = lambda_0
      sd_maxter = sd_maxter_0
      !
      ! ... threshold for empty bands
      !
      empty_ethr = MAX( ( ethr * 50.D0 ), 1.D-3 )
      !
      IF ( iter > 1 ) THEN
         !
         ! ... after the first scf iteration one trial step is supposed to
         ! ... be sufficient. 
         !
         sd_maxter = 1
         !
      END IF
      !
      ! ... initialization of the work-space
      !
      e_old      = 0.D0
      e_old(1,:) = e(:)
      psi_old    = ZERO
      hpsi_old   = ZERO
      spsi_old   = ZERO
      !
      ! ... diagonalization loop
      !
      iterate: DO
         !
         kter = kter + 1
         !
         IF ( kter > maxter ) EXIT iterate
         !
#if defined (DIIS_DEBUG)
         PRINT *, ""
         PRINT *, "kter = ", kter
         PRINT *, ""
#endif
         !
         ! ... the size of the history-subspace is increased ( different for 
         ! ... for each band )
         !
         WHERE ( .NOT. conv(:) )
            !
            nbase(:) = MIN( ( nbase(:) + 1 ), diis_ndim )
            !
         END WHERE        
         !
         ! ... bands are reordered so that converged bands come first
         !
         CALL reorder_bands()
         !
         ! ... here we compute H|psi> and S|psi> for not converged bands
         !
         CALL h_psi( ndmx, ndim, notcnv, psi(:,cnv+1), hpsi(:,cnv+1) )
         CALL s_psi( ndmx, ndim, notcnv, psi(:,cnv+1), spsi(:,cnv+1) )
         !
         ! ... here we set up the hamiltonian and overlap matrix 
         ! ... on the subspace :
         !
         ! ...    hc_ij = <psi_i|H|psi_j>  and  sc_ij = <psi_i|S|psi_j>
         !
         CALL ZGEMM( 'C', 'N', nbands, nbands, ndim, ONE, &
                     psi, ndmx, hpsi, ndmx, ZERO, hc, nbands )
         !
         CALL ZGEMM( 'C', 'N', nbands, nbands, ndim, ONE, &
                     psi, ndmx, spsi, ndmx, ZERO, sc, nbands )
         !
         CALL reduce( 2 * nbands * nbands, hc )
         CALL reduce( 2 * nbands * nbands, sc )
         !
         ! ... the subspace rotation is performed by solving the eigenvalue 
         ! ... problem on the subspace :
         !     
         ! ...    ( hc - e * sc ) |phi> = 0
         !
         CALL cdiaghg( nbands, nbands, hc, sc, nbands, e, vc )
         !
         ! ... convergence is checked here
         !
         frozen = conv
         !
         WHERE( btype(:) == 1 )
            !
            conv(:) = ( ( kter > 1 ) .AND. &
                        ( ABS( e(:) - e_old(1,:) ) < ethr ) )
            !
         ELSEWHERE
            !
            conv(:) = ( ( kter > 1 ) .AND. &
                        ( ABS( e(:) - e_old(1,:) ) < empty_ethr ) )
            !
         END WHERE
         !
         ! ... |psi>,  H|psi>  and  S|psi>  are rotated
         !
         CALL ZGEMM( 'N', 'N', ndim, nbands, nbands, ONE, &
                     psi, ndmx, vc, nbands, ZERO, aux, ndmx ) 
         !     
         psi(:,:) = aux(:,:)
         !
         CALL ZGEMM( 'N', 'N', ndim, nbands, nbands, ONE, &
                     hpsi, ndmx, vc, nbands, ZERO, aux, ndmx ) 
         !     
         hpsi(:,:) = aux(:,:)
         !
         CALL ZGEMM( 'N', 'N', ndim, nbands, nbands, ONE, &
                     spsi, ndmx, vc, nbands, ZERO, aux, ndmx ) 
         !     
         spsi(:,:) = aux(:,:)
         !
#if defined (DIIS_DEBUG)
         PRINT *, "eigenvalues :"
         PRINT *, e(:)
         PRINT *, "variation :"
         PRINT *, ABS( e(:) - e_old(1,:) )
         PRINT *, conv(:)
#endif
         !
         ! ... exit if all bands are converged
         !
         IF ( ALL( conv(:) ) ) EXIT iterate
         !
         ! ... a new step is performed
         !
         bands_loop: DO ib = 1, nbands
            !
            ! ... the two highest bands are always optimized 
            ! ... ( even if converged )
            !
            IF ( conv(ib) .AND. ( ib <= ( nbands - 2 ) ) ) CYCLE bands_loop
            !
            ! ... holes-sniffer
            !
            IF ( ( ib <= ( nbands - 2 ) ) .AND.  &
                 frozen(ib) .AND. ( kter > sd_maxter(ib) ) ) THEN
               !
               ! ... an hole has been detected :  the work-space is cleaned
               !
#if defined (DIIS_DEBUG)               
               PRINT *, "HOLES-SNIFFER: ", &
                        "an hole has been detected near band ", ib
#endif               
               !
               nbase(ib) = 1
               !
               e_old(:,ib)      = 0.D0
               psi_old(:,:,ib)  = ZERO
               hpsi_old(:,:,ib) = ZERO
               spsi_old(:,:,ib) = ZERO
               !
               sd_maxter(ib) = kter + sd_maxter_0               
               !
               maxter = MAX( maxter, sd_maxter(ib) )
               !
            END IF
            !
            ! ... the two highest bands are never optimized with the 
            ! ... DIIS procedure
            !
            IF ( ( kter <= sd_maxter(ib) ) .OR. ( ib > ( nbands - 2 ) ) ) THEN
               !
               ! ... trial step :
               !     
               IF ( ( e(ib) < e_old(1,ib) ) .OR. ( kter == 1 ) ) THEN
                  !
                  IF ( kter == 1 ) THEN
                     !
                     ! ... first step :  lambda is kept fixed
                     !
                  ELSE
                     !
                     ! ... step accepted :  lambda is increased
                     !
                     lambda(ib) = lambda(ib) * 1.25D0
                     !
                  END IF   
                  !
               ELSE
                  !
                  ! ... step rejected :  lambda is decreased
                  !
                  lambda(ib) = 0.75D0 * lambda(ib)
                  !
               END IF
               !
               ! ... here we compute the residual vector ( stored in aux ) :
               !
               ! ...    |res> = ( H - e * S ) |psi>           
               !
               aux(:,ib) = hpsi(:,ib) - e(ib) * spsi(:,ib)
               !
            ELSE   
               !
               ! ... DIIS step :  the best eigenvector and residual vector 
               ! ...              are computed
               !
               CALL diis_step( ib )           
               !
            END IF
            !
            ! ... DIIS work-space is refreshed ( for not converged bands only )
            !
            DO n = ( diis_ndim - 1 ), 2, -1
               !
               e_old(n,ib)      = e_old(n-1,ib)
               psi_old(:,n,ib)  = psi_old(:,n-1,ib)
               hpsi_old(:,n,ib) = hpsi_old(:,n-1,ib)                 
               spsi_old(:,n,ib) = spsi_old(:,n-1,ib)
               !
            END DO
            !
            e_old(1,ib)      = e(ib)
            psi_old(:,1,ib)  = psi(:,ib)   
            hpsi_old(:,1,ib) = hpsi(:,ib)
            spsi_old(:,1,ib) = spsi(:,ib)
            !
         END DO bands_loop
         !
#if defined (CG_STEP)
         !
         CALL cg_step()
         !
#else         
         !
         ! ... preconditioning of all residual vecotrs
         !
         CALL g_psi( ndmx, ndim, nbands, aux, e )           
         !
#endif         
         !           
         ! ... here we compute the new eigenvectors for "non converged"
         ! ... bands only 
         ! ... ( the highest two bands are always optimized, even if converged )
         !
         FORALL( ib = 1: nbands, &
                 ( .NOT. conv(ib) ) .OR. ( ib > ( nbands - 2 ) ) )
            !
            psi(:,ib) = psi(:,ib) - lambda(ib) * aux(:,ib)
            !
         END FORALL   
         !
      END DO iterate
      !
      ! ... this is an overestimate of the real number of iterations
      !
      diis_iter = MIN( kter, maxter )
      !
      ! ... the number of not-converged bands is computed
      !
      notcnv = nbands
      !
      DO ib = 1, nbands
         !
         IF ( conv(ib) ) notcnv = notcnv - 1
         !
      END DO
      !
      ! ... work-space deallocation
      !
      DEALLOCATE( frozen )
      DEALLOCATE( conv )
      DEALLOCATE( sd_maxter )
      DEALLOCATE( nbase )
      DEALLOCATE( lambda )
      DEALLOCATE( vc )
      DEALLOCATE( sc )
      DEALLOCATE( hc )
      DEALLOCATE( aux )
      DEALLOCATE( hpsi )
      DEALLOCATE( spsi )
      DEALLOCATE( e_old )
      DEALLOCATE( psi_old )
      DEALLOCATE( hpsi_old )
      DEALLOCATE( spsi_old )
      !
      CALL stop_clock( 'diis' )
      !
      RETURN
      !
      CONTAINS
        !
        ! ... internal procedures
        !
        !--------------------------------------------------------------------
        SUBROUTINE reorder_bands()
          !--------------------------------------------------------------------
          !
          ! ... this routine is used to reorder the bands :
          ! ... converged bands come first.
          ! ... for this pourpose an auxiliary vector is used
          !
          IMPLICIT NONE
          !
          !
          cnv    = 0
          notcnv = nbands
          !
          IF ( ( kter <= 2 ) .OR. ( .NOT. ANY( conv(:) ) ) ) RETURN
          !
          DO ib = 1, nbands
             !
             IF ( conv(ib) .AND. ( ib <= ( nbands - 2 ) ) ) THEN
                !
                cnv = cnv + 1
                !
                aux(:,cnv) = psi(:,ib)
                !
                hpsi(:,cnv) = hpsi(:,ib)
                spsi(:,cnv) = spsi(:,ib)              
                !
             ELSE
                !
                aux(:,notcnv) = psi(:,ib)
                !
                notcnv = notcnv - 1
                !
             END IF
             !
          END DO
          !
          notcnv = nbands - notcnv
          !
          psi(:,:) = aux(:,:)
          !
          RETURN
          !
        END SUBROUTINE reorder_bands
        !
#if defined (CG_STEP)
        !
        !--------------------------------------------------------------------
        SUBROUTINE cg_step()
          !--------------------------------------------------------------------
          !
          IMPLICIT NONE
          !
          COMPLEX (KIND=DP)              :: gamma, num, den
          COMPLEX (KIND=DP), ALLOCATABLE :: res_old(:), pres(:), pres_old(:)
          !
          COMPLEX (KIND=DP), EXTERNAL :: ZDOTC      
          !
          !
          ALLOCATE( res_old(  ndmx ) )
          ALLOCATE( pres(     ndmx ) )
          ALLOCATE( pres_old( ndmx ) )          
          !
          bands_loop: DO ib = 1, nbands
             !
             IF ( nbase(ib) >= 2 ) THEN
                !
                ! ... preconditioned conjugate gradients step
                !
                IF ( conv(ib) ) CYCLE bands_loop
                !
                pres    = ( hpsi_old(:,1,ib) - e_old(1,ib) * spsi_old(:,1,ib) )
                res_old = ( hpsi_old(:,2,ib) - e_old(2,ib) * spsi_old(:,2,ib) )
                !
                pres_old = res_old
                !
                CALL g_psi( ndmx, ndim, 1, pres,     e_old(1,ib) )
                CALL g_psi( ndmx, ndim, 1, pres_old, e_old(2,ib) )
                !
                num = ZDOTC( ndim, pres, 1, aux(:,ib), 1 )
                !
                CALL reduce( 2, num )
                !
                den = ZDOTC( ndim, pres_old, 1, res_old, 1 )
                !
                CALL reduce( 2, den )
                !
                gamma = num / den
                !
                aux(:,ib) = pres - gamma * ( psi_old(:,1,ib) - psi_old(:,2,ib) )
                !
             ELSE
                !
                ! ... standard preconditioned steepest descent step
                !
                CALL g_psi( ndmx, ndim, 1, aux(:,ib), e(ib) )
                !
             END IF
             !
          END DO bands_loop 
          !
          DEALLOCATE( res_old )
          DEALLOCATE( pres )
          DEALLOCATE( pres_old )
          !
          RETURN
          !
        END SUBROUTINE cg_step
        !
#endif        
        !
        !--------------------------------------------------------------------
        SUBROUTINE diis_step( ib )
          !--------------------------------------------------------------------
          !
          IMPLICIT NONE
          !
          INTEGER, INTENT(IN)            :: ib
          INTEGER                        :: dim, n, m
          REAL (KIND=DP)                 :: psiSpsi
          REAL (KIND=DP),    ALLOCATABLE :: e_small(:)
          COMPLEX (KIND=DP), ALLOCATABLE :: rc_small(:,:), &
                                            vc_small(:,:)
          COMPLEX (KIND=DP), ALLOCATABLE :: all_psi(:,:),  &
                                            all_hpsi(:,:), &
                                            all_spsi(:,:), &
                                            all_res(:,:)
          !
          COMPLEX(KIND=DP), EXTERNAL :: ZDOTC
          !
          !
          dim = nbase(ib)
          !
          ! ... internal work-space allocation
          !
          ALLOCATE( e_small( dim ) )
          ALLOCATE( rc_small( dim, dim ) )
          ALLOCATE( vc_small( dim, dim ) )
          ALLOCATE( all_psi(  ndmx, dim ) )
          ALLOCATE( all_hpsi( ndmx, dim ) )
          ALLOCATE( all_spsi( ndmx, dim ) )
          ALLOCATE( all_res(  ndmx, dim ) )
          !
          ! ... the history of this band is reconstructed
          !
          all_psi(:,1)  = psi(:,ib)
          all_hpsi(:,1) = hpsi(:,ib)
          all_spsi(:,1) = spsi(:,ib)
          e_small(1)    = e(ib)
          !
          all_psi(:,2:dim)  = psi_old(:,:,ib)
          all_hpsi(:,2:dim) = hpsi_old(:,:,ib)
          all_spsi(:,2:dim) = spsi_old(:,:,ib)
          e_small(2:dim)    = e_old(:,ib)
          !
          ! ... orthogonalization
          !
          CALL cgramg1( ndmx, dim, ndim, 1, dim, all_psi, all_spsi, all_hpsi )
          !
          FORALL( n = 1: dim ) &
             all_res(:,n) = ( all_hpsi(:,n) - e_small(n) * all_spsi(:,n) )
          !   
          ! ... here we construct the matrix :  rc_ij = <res_i|res_j>
          !
          CALL ZGEMM( 'C', 'N', dim, dim, ndim, ONE, all_res(:,:), &
                      ndmx, all_res(:,:), ndmx, ZERO, rc_small, dim )
          !
          CALL reduce( 2 * dim * dim, rc_small )
          ! 
          ! ... diagonalize the reduced hamiltonian
          !
          CALL cdiagh( dim, rc_small, dim, e_small, vc_small )
          !
          ! ... here we compute the best estimate of the |psi>, H|psi>, S|psi>
          !
          CALL ZGEMM( 'N', 'N', ndim, 1, dim, ONE, all_psi(:,:), &
                      ndmx, vc_small(:,1), dim, ZERO, psi(:,ib), ndmx )
          !
          CALL ZGEMM( 'N', 'N', ndim, 1, dim, ONE, all_hpsi(:,:), &
                      ndmx, vc_small(:,1), dim, ZERO, hpsi(:,ib), ndmx )
          !
          CALL ZGEMM( 'N', 'N', ndim, 1, dim, ONE, all_spsi(:,:), &
                      ndmx, vc_small(:,1), dim, ZERO, spsi(:,ib), ndmx )
          !
          psiSpsi = REAL( ZDOTC( ndim, psi(:,ib), 1, spsi(:,ib), 1 ) )
          !
          CALL reduce( 1, psiSpsi )
          !
          e(ib) = REAL( ZDOTC( ndim, psi(:,ib), 1, hpsi(:,ib), 1 ) ) / psiSpsi
          !
          CALL reduce( 1, e(ib) )
          !
          ! ... here we compute the best estimate of the residual vector
          !
          aux(:,ib) = hpsi(:,ib) - e(ib) * spsi(:,ib)
          !
          ! ... internal work-space deallocation
          !
          DEALLOCATE( e_small )
          DEALLOCATE( rc_small )
          DEALLOCATE( vc_small )
          DEALLOCATE( all_psi )
          DEALLOCATE( all_hpsi )
          DEALLOCATE( all_spsi )
          DEALLOCATE( all_res )
          !   
          RETURN
          !
        END SUBROUTINE diis_step
        !
    END SUBROUTINE cdiisg
    !
END MODULE diis_module
