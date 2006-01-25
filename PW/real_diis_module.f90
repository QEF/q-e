!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#define DIIS_DEBUG
!
!#define WINDOW_ORTHO
!#define SHOW_OVERLAP
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
#include "f_defs.h"
! 
! ... Description of the DIIS algorithm in diis_base.f90
!
! ... written by Carlo Sbraccia ( 08/06/2004 )
!  
!
!----------------------------------------------------------------------------
MODULE real_diis_module
  !----------------------------------------------------------------------------
  !
  ! ... real hamiltonian case
  !
  USE gvect, ONLY : gstart
  USE diis_base
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  ! ... public methods
  !
  PUBLIC :: rdiisg
  !
  ! ... specific real DIIS variables :
  !
  SAVE
  !
  REAL (DP), ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:)
    ! H matrix on the reduced basis
    ! S matrix on the reduced basis
    ! the eigenvectors of the Hamiltonian
  !
  ! ... external functions
  !
  REAL (DP), EXTERNAL :: DDOT
  !
  ! ... module procedures :
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE rdiisg( ndim, ndmx, nbnd, psi, e, &
                       btype, notcnv, diis_iter, iter )
      !------------------------------------------------------------------------
      !
      ! ... gamma_only version of the DIIS algorithm
      !
      USE gvect,         ONLY : gstart
      !
      USE g_psi_mod,     ONLY : h_diag
      USE control_flags, ONLY : diis_ndim, ethr, istep
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      INTEGER, INTENT(IN) :: ndim, ndmx, nbnd, iter
        ! dimension of the matrix to be diagonalized
        ! leading dimension of matrix psi, as declared in the calling pgm unit
        ! integer number of searched low-lying roots
        ! scf iteration
      INTEGER, INTENT(INOUT) :: btype(nbnd)
        ! band type ( 1 = occupied, 0 = empty )
      COMPLEX (DP), INTENT(INOUT) :: psi(ndmx,nbnd)
        !  psi contains the refined estimates of the eigenvectors
      REAL (DP), INTENT(INOUT) :: e(nbnd)
        ! contains the estimated eigenvalues
      INTEGER, INTENT(OUT) :: diis_iter
        ! average number of iterations performed per band
      INTEGER , INTENT(OUT):: notcnv
        ! number of unconverged bands
      !
      ! ... local variables
      !
      INTEGER :: ib
      !
      !
      CALL start_clock( 'diis' )
      !
      ! ... initialization of control variables
      !
      nbase      = 0
      diis_iter  = 0
      ndim2      = 2 * ndim
      ndmx2      = 2 * ndmx
      nbnd_diis  = ( nbnd - ncgbnd )
      diis_ndim1 = ( diis_ndim - 1 )
      !
      ! ... work-space allocation
      !
      CALL allocate_base( ndmx, nbnd )
      !
      ALLOCATE( hr( nbnd, nbnd ) )
      ALLOCATE( sr( nbnd, nbnd ) )
      ALLOCATE( vr( nbnd, nbnd ) )
      !
      ! ... again control arrays
      !
      conv = .FALSE.
      !
      ! ... initialization of the work-space
      !
      e_old    = 0.D0
      psi_old  = ZERO
      hpsi_old = ZERO
      spsi_old = ZERO
      !
      ! ... threshold for empty bands
      !
      empty_ethr = MAX( ( ethr * 5.D0 ), empty_bands_ethr_min )
      !
#if defined (DIIS_DEBUG)
      PRINT *, "input eigenvalues :"
      PRINT *, e(:)
#endif            
      !
      ! ... initialization
      !
      IF ( ( iter == 1 ) .AND. ( istep == 1 ) ) THEN
         !
         ! ... some sweeps over bands are performed at the first scf
         ! ... iteration of the first ionic step :  
         ! ... this to reduce the amount of holes in the energy spectrum
         !
         CALL init_steps( ndim, ndmx, nbnd, psi, e, btype, diis_iter )
         !
         ! ... first check on convergence
         !
         IF ( ALL( conv(:) ) ) THEN
            !
            notcnv = 0
            !
            CALL terminate()
            !
            RETURN
            !
         END IF
         !
      ELSE
         !
         ! ... initialization step for DIIS
         !
         ! ... here we compute H|psi> and S|psi> for all the bands
         !
         CALL h_psi( ndmx, ndim, nbnd, psi(1,1), hpsi(1,1) )
         CALL s_psi( ndmx, ndim, nbnd, psi(1,1), spsi(1,1) )
         !
         ! ... here we set up the hamiltonian and overlap matrix
         ! ... on the subspace :
         !
         ! ...    hr_ij = <psi_i|H|psi_j>  and  sr_ij = <psi_i|S|psi_j>
         !
         CALL DGEMM( 'T', 'N', nbnd, nbnd, ndim2, 2.D0, &
                     psi, ndmx2, hpsi, ndmx2, 0.D0, hr, nbnd )
         CALL DGEMM( 'T', 'N', nbnd, nbnd, ndim2, 2.D0, &
                     psi, ndmx2, spsi, ndmx2, 0.D0, sr, nbnd )
         !
         IF ( gstart == 2 ) THEN
            !
            CALL DGER( nbnd, nbnd, -1.D0, psi, ndmx2, hpsi, ndmx2, hr, nbnd )
            CALL DGER( nbnd, nbnd, -1.D0, psi, ndmx2, spsi, ndmx2, sr, nbnd )
            !
         END IF              
         !
         CALL reduce( nbnd * nbnd, hr )
         CALL reduce( nbnd * nbnd, sr )
         !
         ! ... the subspace rotation is performed by solving the eigenvalue
         ! ... problem on the subspace :
         !
         ! ...    ( hr - e * sr ) |phi> = 0
         !
         CALL rdiaghg( nbnd, nbnd, hr, sr, nbnd, e, vr )
         !
         ! ... |psi>,  H|psi>  and  S|psi>  are rotated
         !
         CALL DGEMM( 'N', 'N', ndim2, nbnd, nbnd, 1.D0, &
                     psi, ndmx2, vr, nbnd, 0.D0, aux, ndmx2 )
         !
         psi(:,:) = aux(:,:)
         !
         CALL DGEMM( 'N', 'N', ndim2, nbnd, nbnd, 1.D0, &
                     hpsi, ndmx2, vr, nbnd, 0.D0, aux, ndmx2 )
         !
         hpsi(:,:) = aux(:,:)
         !
         CALL DGEMM( 'N', 'N', ndim2, nbnd, nbnd, 1.D0, &
                     spsi, ndmx2, vr, nbnd, 0.D0, aux, ndmx2 )
         !
         spsi(:,:) = aux(:,:)
         !
      END IF
      !
      ! ... then DIIS-diagonalization is performed only on the lowest 
      ! ... "nbnd_diis" bands
      !
      CALL diis_with_ortho( ndim, ndmx, psi, e, btype, diis_iter )
      !
      ! ... finally, all eventual holes left by the DIIS-diagonalization are
      ! ... filled with a CG-based diagonalization of the topmost "ncgbnd" 
      ! ... bands
      !
      holes_sniffer: DO
         !
         CALL holes_filler( ndmx, ndim, nbnd, psi, e, h_diag, diis_iter )
         !
         IF ( no_holes( ndmx, nbnd, btype, psi, e ) ) EXIT holes_sniffer
         ! 
      END DO  holes_sniffer
      !
      ! ... the number of not-converged bands is computed
      !
      notcnv = COUNT( .NOT. conv(:) )
      !
      CALL terminate()
      !
      RETURN
      !
      CONTAINS
        !
        !--------------------------------------------------------------------
        SUBROUTINE terminate()
          !--------------------------------------------------------------------
          !
          IMPLICIT NONE
          !
          ! ... work-space deallocation
          !
          CALL deallocate_base()
          !
          DEALLOCATE( vr )
          DEALLOCATE( sr )
          DEALLOCATE( hr )
          !
          CALL stop_clock( 'diis' )
          !
#if defined (DIIS_DEBUG)
          PRINT *, "final eigenvalues :"
          PRINT *, e(:)
          PRINT *, "variation :"
          PRINT *, ABS( e(:) - e_ref(:) )
          PRINT *, conv(:)
#endif          
          !
          RETURN
          !
        END SUBROUTINE terminate
        !
    END SUBROUTINE rdiisg
    !
    !------------------------------------------------------------------------
    SUBROUTINE init_steps( ndim, ndmx, nbnd, psi, e, btype, diis_iter )
      !------------------------------------------------------------------------
      !
      USE control_flags, ONLY : ethr
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      INTEGER, INTENT(IN) :: ndim, ndmx, nbnd
        ! dimension of the matrix to be diagonalized
        ! leading dimension of matrix psi, as declared in the calling pgm unit
        ! integer number of searched low-lying roots
      INTEGER, INTENT(INOUT) :: btype(nbnd)
        ! band type ( 1 = occupied, 0 = empty )
      COMPLEX (DP), INTENT(INOUT) :: psi(ndmx,nbnd)
        !  psi contains the refined estimates of the eigenvectors
      REAL (DP), INTENT(INOUT) :: e(nbnd)
        ! contains the estimated eigenvalues
      INTEGER, INTENT(INOUT) :: diis_iter  
        ! average number of iterations performed per band         
      !
      ! ... local variables
      !
      INTEGER        :: sweep, trial_step
      INTEGER        :: ib
      INTEGER        :: notcnv
      REAL (DP) :: psiSpsi
      !
      !
      sweep = 0
      !
      CALL h_psi( ndmx, ndim, nbnd, psi(1,1), hpsi(1,1) )
      CALL s_psi( ndmx, ndim, nbnd, psi(1,1), spsi(1,1) )
      !
      iterate: DO
         !
         sweep = sweep + 1
         !
         IF ( sweep > max_sweeps ) EXIT iterate
         !
#if defined (DIIS_DEBUG)
         WRITE( *, '(/"sweep = ",I3/)') sweep
#endif
         !
         ! ... old eigenvalues are saved
         !
         e_ref(:) = e(:)
         !
         DO trial_step = 1, max_trial_steps
            !
            ! ... residual vectors
            !
            DO ib = 1, nbnd
               !
               IF ( conv(ib) ) CYCLE
               !
               psiSpsi = 2.D0 * DDOT( ndim2, psi(1,ib), 1, spsi(1,ib), 1 )
               !
               IF ( gstart == 2 ) psiSpsi = psiSpsi - psi(1,ib) * spsi(1,ib)
               !
               CALL reduce( 1, psiSpsi )
               !
               psiSpsi = 1.D0 / psiSpsi
               !
               psi(1,ib)  = psi(1,ib)  * psiSpsi
               hpsi(1,ib) = hpsi(1,ib) * psiSpsi
               spsi(1,ib) = spsi(1,ib) * psiSpsi
               !
               e(ib) = 2.D0 * DDOT( ndim2, psi(1,ib), 1, hpsi(1,ib), 1 )
               !
               IF ( gstart == 2 ) e(ib) = e(ib) - psi(1,ib) * hpsi(1,ib)
               !
               CALL reduce( 1, e(ib) )
               !
               aux(:,ib) = hpsi(:,ib) - e(ib) * spsi(:,ib)
               !
               CALL g_psi( ndmx, ndim, 1, aux(1,ib), e(ib) )               
               !
            END DO
            !           
            ! ... trial step
            !
            FORALL( ib = 1: nbnd, .NOT. conv(ib) )
               !
               psi(:,ib) = psi(:,ib) - sweeps_lambda * aux(:,ib)
               !
            END FORALL
            !
            ! ... bands are reordered so that converged bands come first
            !
            CALL reorder_bands( psi, notcnv, nbnd, +1 )
            !
            ! ... here we compute H|psi> and S|psi> for not converged bands
            !
            CALL h_psi( ndmx, ndim, notcnv, psi(1,cnv+1), hpsi(1,cnv+1) )
            CALL s_psi( ndmx, ndim, notcnv, psi(1,cnv+1), spsi(1,cnv+1) ) 
            !
            ! ... bands are back-ordered
            !
            CALL reorder_bands( psi, notcnv, nbnd, -1 )                       
            !
         END DO
         !
         ! ... here we set up the hamiltonian and overlap matrix
         ! ... on the subspace :
         !
         ! ...    hr_ij = <psi_i|H|psi_j>  and  sr_ij = <psi_i|S|psi_j>
         !
         CALL DGEMM( 'T', 'N', nbnd, nbnd, ndim2, 2.D0, &
                     psi, ndmx2, hpsi, ndmx2, 0.D0, hr, nbnd )
         CALL DGEMM( 'T', 'N', nbnd, nbnd, ndim2, 2.D0, &
                     psi, ndmx2, spsi, ndmx2, 0.D0, sr, nbnd )
         !
         IF ( gstart == 2 ) THEN
            !
            CALL DGER( nbnd, nbnd, -1.D0, psi, ndmx2, hpsi, ndmx2, hr, nbnd )
            CALL DGER( nbnd, nbnd, -1.D0, psi, ndmx2, spsi, ndmx2, sr, nbnd )
            !
         END IF              
         !
         CALL reduce( nbnd * nbnd, hr )
         CALL reduce( nbnd * nbnd, sr )
         !
         ! ... the subspace rotation is performed by solving the eigenvalue
         ! ... problem on the subspace :
         !
         ! ...    ( hr - e * sr ) |phi> = 0
         !
         CALL rdiaghg( nbnd, nbnd, hr, sr, nbnd, e, vr )         
         !
         ! ... |psi>  are rotated
         !
         CALL DGEMM( 'N', 'N', ndim2, nbnd, nbnd, 1.D0, &
                     psi, ndmx2, vr, nbnd, 0.D0, aux, ndmx2 )
         !
         psi(:,:) = aux(:,:)
         !
         ! ... convergence is checked here
         !
         WHERE( btype(:) == 1 )
            !
            conv(:) = ( ( ABS( e(:) - e_ref(:) ) < ethr ) )
            !
         ELSEWHERE
            !
            conv(:) = ( ( ABS( e(:) - e_ref(:) ) < empty_ethr ) )
            !
         END WHERE
         !
         ! ... exit if all bands are converged
         !
         IF ( ALL( conv(:) ) ) EXIT iterate
         !
         ! ... H|psi>  and  S|psi>  are rotated
         !
         CALL DGEMM( 'N', 'N', ndim2, nbnd, nbnd, 1.D0, &
                     hpsi, ndmx2, vr, nbnd, 0.D0, aux, ndmx2 )
         !     
         hpsi(:,:) = aux(:,:)
         !
         CALL DGEMM( 'N', 'N', ndim2, nbnd, nbnd, 1.D0, &
                     spsi, ndmx2, vr, nbnd, 0.D0, aux, ndmx2 )
         !     
         spsi(:,:) = aux(:,:)
         !
#if defined (DIIS_DEBUG)
         PRINT *, "eigenvalues :"
         PRINT *, e(:)
         PRINT *, "variation :"
         PRINT *, ABS( e(:) - e_ref(:) )
         PRINT *, conv(:)
#endif
         !
      END DO iterate
      !
      ! ... if some band is not yet coverged a final trial-step is done
      !
      IF ( ANY( .NOT. conv(:) ) ) THEN
         !
         nbase = 1
         !
         DO ib = 1, nbnd
            !
            IF ( conv(ib) ) CYCLE
            !
            IF ( ib <= nbnd_diis ) THEN
               !
               ! ... the first step in the DIIS history is saved
               !
               e_old(1,ib)      = e(ib)
               psi_old(:,1,ib)  = psi(:,ib)
               hpsi_old(:,1,ib) = hpsi(:,ib)
               spsi_old(:,1,ib) = spsi(:,ib)
               !
            END IF   
            !
            aux(:,ib) = hpsi(:,ib) - e(ib) * spsi(:,ib)
            !
            CALL g_psi( ndmx, ndim, 1, aux(1,ib), e(ib) )
            !
            psi(:,ib) = psi(:,ib) - diis_lambda * aux(:,ib)
            !
         END DO   
         !
      END IF
      !
      ! ... this is an overestimate of the real number of iterations
      !
      diis_iter = diis_iter + MIN( sweep, max_sweeps ) * 2
      !
      RETURN
      !
    END SUBROUTINE init_steps
    !
    !------------------------------------------------------------------------
    SUBROUTINE diis_with_ortho( ndim, ndmx, psi, e, btype, diis_iter )
      !------------------------------------------------------------------------
      !
      USE control_flags, ONLY : diis_ndim, ethr
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      INTEGER, INTENT(IN) :: ndim, ndmx
        ! dimension of the matrix to be diagonalized
        ! leading dimension of matrix psi, as declared in the calling pgm unit
      INTEGER, INTENT(INOUT) :: btype(nbnd_diis)
        ! band type ( 1 = occupied, 0 = empty )
      COMPLEX (DP), INTENT(INOUT) :: psi(ndmx,nbnd_diis)
        !  psi contains the refined estimates of the eigenvectors
      REAL (DP), INTENT(INOUT) :: e(nbnd_diis)
        ! contains the estimated eigenvalues
      INTEGER, INTENT(INOUT) :: diis_iter  
        ! average number of iterations performed per band                 
      !
      ! ... local variables
      !
      INTEGER        :: kter, ib, jb
      INTEGER        :: notcnv
      REAL (DP) :: psi_norm, overlap
      !
      !
#if defined (WINDOW_ORTHO)
      !
      ortho_win = MAX( ortho_win_min, &
                       0.05D0 * ( MAXVAL( e(:) - MINVAL( e(:) ) ) ) )
      !
#endif
      !
      ! ... the first step in the DIIS history is saved
      !
      FORALL( ib = 1: nbnd_diis )
         !
         e_old(1,ib)      = e(ib)
         psi_old(:,1,ib)  = psi(:,ib)
         hpsi_old(:,1,ib) = hpsi(:,ib)
         spsi_old(:,1,ib) = spsi(:,ib)
         !
      END FORALL
      !
      nbase = 1
      !
      ! ... new residual vectors
      !
      FORALL( ib = 1 : nbnd_diis )
         !
         aux(:,ib) = hpsi(:,ib) - e(ib) * spsi(:,ib)
         !
      END FORALL
      !
      ! ... preconditioning of all the residual vecotrs
      !
      CALL g_psi( ndmx, ndim, nbnd_diis, aux, e )
      !
      ! ... new trial wavefunctions
      !
      FORALL( ib = 1 : nbnd_diis )
         !
         psi(:,ib) = psi(:,ib) - diis_lambda * aux(:,ib)
         !
      END FORALL
      !
      kter = 1
      !
      ! ... DIIS loop
      !
      iterate: DO
         !
         kter = kter + 1
         !
         IF ( kter > maxter ) EXIT iterate
         !
#if defined (DIIS_DEBUG)
         WRITE( *, '(/"kter = ",I3/)') kter
#endif
         !
         ! ... reference eigenvalues are saved
         !
         e_ref(1:nbnd_diis) = e(1:nbnd_diis)
         !
         ! ... the size of the history-subspace is increased
         !
         nbase = MIN( ( nbase + 1 ), diis_ndim )
         !
         ! ... bands are reordered so that converged bands come first
         !
         CALL reorder_bands( psi, notcnv, nbnd_diis, +1 )
         !
         ! ... here we compute H|psi> and S|psi> for not converged bands
         !
         CALL h_psi( ndmx, ndim, notcnv, psi(1,cnv+1), hpsi(1,cnv+1) )
         CALL s_psi( ndmx, ndim, notcnv, psi(1,cnv+1), spsi(1,cnv+1) )
         !
         ! ... bands are back-ordered
         !
         CALL reorder_bands( psi, notcnv, nbnd_diis, -1 )
         !
         ! ... DIIS step : the best |psi>, H|psi> and S|spi> are computed here
         !
         CALL diis_step()
         !
         ! ... orthogonalization step of all the best eigenvectors ( an energy 
         ! ... window "ortho_win" is used here )
         !
#if defined (WINDOW_ORTHO)         
         !
         CALL start_clock( 'cgramg1' )
         !
         DO ib = 1, nbnd_diis
            !
#if defined (SHOW_OVERLAP)
            WRITE( 999, '(//"VECTOR = ",I3)' ) ib
#endif
            !          
            DO jb = 1, ( ib - 1 )
               !
               IF ( ( e(ib) - e(jb) ) > ortho_win ) CYCLE 
               !
               overlap = 2.D0 * DDOT( ndim2, psi(1,jb), 1, spsi(1,ib), 1 )
               !
               IF ( gstart == 2 ) overlap = overlap - psi(1,jb) * spsi(1,ib)
               !
               CALL reduce( 1, overlap )
               !
#if defined (SHOW_OVERLAP)
               WRITE( 999, '("OVERLAP(",I3,",",I3,") = ",2(2X,F14.10))' ) &
                   jb, ib, overlap, ( e(ib) - e(jb) )
#endif
               !
               psi(:,ib)  = psi(:,ib)  - overlap * psi(:,jb)
               hpsi(:,ib) = hpsi(:,ib) - overlap * hpsi(:,jb)
               spsi(:,ib) = spsi(:,ib) - overlap * spsi(:,jb)
               !
            END DO
            !
            psi_norm = 2.D0 * DDOT( ndim2, psi(1,ib), 1, spsi(1,ib), 1 )
            !
            IF ( gstart == 2 ) psi_norm = psi_norm - psi(1,ib) * spsi(1,ib)
            !
            CALL reduce( 1, psi_norm )
            !
            IF ( psi_norm < eps32 ) THEN
               !
               PRINT *, psi_norm
               !
               CALL errore( 'diis_step', ' negative norm in S ', 1 )
               !
            END IF
            !
            psi_norm = 1.D0 / SQRT( psi_norm )
            !
            psi(:,ib)  = psi_norm * psi(:,ib)
            hpsi(:,ib) = psi_norm * hpsi(:,ib)
            spsi(:,ib) = psi_norm * spsi(:,ib)
            !
         END DO
         !
         CALL stop_clock( 'cgramg1' )
         !
#else
         !
         CALL cgramg1( ndmx, nbnd_diis, ndim, 1, nbnd_diis, psi, spsi, hpsi )
         !
#endif
         !         
         ! ... convergence is checked here
         !
         WHERE( btype(1:nbnd_diis) == 1 )
            !
            conv(1:nbnd_diis) = ( ABS( e(1:nbnd_diis) - e_ref(1:nbnd_diis) ) < ethr )
            !
         ELSEWHERE
            !
            conv(1:nbnd_diis) = ( ABS( e(1:nbnd_diis) - e_ref(1:nbnd_diis) ) < empty_ethr )
            !
         END WHERE
         !
#if defined (DIIS_DEBUG)
         PRINT *, "eigenvalues :"
         PRINT *, e(1:nbnd_diis)
         PRINT *, "variation :"
         PRINT *, ABS( e(1:nbnd_diis) - e_ref(1:nbnd_diis) )
         PRINT *, conv(1:nbnd_diis)
#endif         
         !
         ! ... exit if all bands are converged
         !
         IF ( ALL( conv(1:nbnd_diis) ) ) EXIT iterate
         !
         ! ... new preconditioned residual vectors
         !
         DO ib = 1, nbnd_diis
            !
            IF ( conv(ib) ) CYCLE
            !
            aux(:,ib) = hpsi(:,ib) - e(ib) * spsi(:,ib)
            !
            CALL g_psi( ndmx, ndim, 1, aux(1,ib), e(ib) )
            !
         END DO
         !           
         ! ... here we compute the new eigenvectors for not converged
         ! ... bands only 
         !
         FORALL( ib = 1: nbnd_diis, .NOT. conv(ib) )
            !
            psi(:,ib) = psi(:,ib) - diis_lambda * aux(:,ib)
            !
         END FORALL
         !
      END DO iterate
      !
#if defined (SHOW_OVERLAP)
      !
      WRITE( 999, '(//"FINAL CHECK ON ORTHOGONALIZATION ")' )
      !
      DO ib = 1, nbnd_diis
         !
         WRITE( 999, '(//"VECTOR = ",I3)' ) ib
         !          
         DO jb = 1, ib
            !
            overlap = 2.D0 * DDOT( ndim2, psi(1,jb), 1, spsi(1,ib), 1 )
            !
            IF ( gstart == 2 ) overlap = overlap - psi(1,jb) * spsi(1,ib)
            !
            CALL reduce( 1, overlap )
            !
            WRITE( 999, '("OVERLAP(",I3,",",I3,") = ",3(2X,F14.10))' ) &
                jb, ib, overlap, ( e(ib) - e(jb) )
            !
         END DO
         !
      END DO
#endif
      !
      ! ... this is an overestimate of the real number of iterations
      !
      diis_iter = diis_iter + MIN( kter, maxter )
      !
      RETURN
      !
      CONTAINS
        !
        !--------------------------------------------------------------------
        SUBROUTINE diis_step()
          !--------------------------------------------------------------------
          !
          IMPLICIT NONE
          !
          INTEGER                        :: ib, n, dim, dim_new
          REAL (DP)                 :: psiSpsi
          REAL (DP),    ALLOCATABLE :: e_small(:)
          REAL (DP),    ALLOCATABLE :: rr_small(:,:), &
                                            sr_small(:,:), &
                                            vr_small(:,:)
          COMPLEX (DP), ALLOCATABLE :: all_psi(:,:),  &
                                            all_hpsi(:,:), &
                                            all_spsi(:,:), &
                                            all_res(:,:)
          REAL (DP),    ALLOCATABLE :: ps(:)
          !
          !
          dim     = MIN( ( nbase - 1 ), diis_ndim1 )
          dim_new = MIN( nbase, diis_ndim1 )
          !
          ! ... internal work-space allocation
          !
          ALLOCATE( e_small( nbase ) )
          ALLOCATE( rr_small( nbase, nbase ) )
          ALLOCATE( sr_small( nbase, nbase ) )
          ALLOCATE( vr_small( nbase, nbase ) )
          ALLOCATE( all_psi(  ndmx, nbase ) )
          ALLOCATE( all_hpsi( ndmx, nbase ) )
          ALLOCATE( all_spsi( ndmx, nbase ) )
          ALLOCATE( all_res(  ndmx, nbase ) )
          ALLOCATE( ps( nbase ) )
          !
          ! ... initialization
          !
          e_small  = 0.D0
          rr_small = 0.D0
          sr_small = 0.D0
          vr_small = 0.D0
          all_psi  = ZERO
          all_hpsi = ZERO
          all_spsi = ZERO
          all_res  = ZERO
          !
          bands_loop: DO ib = 1, nbnd_diis
             !
             IF ( conv(ib) ) CYCLE bands_loop
             !
             ! ... the history of this band is reconstructed
             !
             all_psi(:,1)  = psi(:,ib)
             all_hpsi(:,1) = hpsi(:,ib)
             all_spsi(:,1) = spsi(:,ib)
             !
             all_psi(:,2:nbase)  = psi_old(:,1:dim,ib)
             all_hpsi(:,2:nbase) = hpsi_old(:,1:dim,ib)
             all_spsi(:,2:nbase) = spsi_old(:,1:dim,ib)
             !
             ! ... orthogonalization of the new wavefunction to the history
             !
             DO n = 1, nbase
                ! 
                ps(n) = 2.D0 * DDOT( ndim2, all_psi(1,1), 1, all_spsi(1,n), 1 )
                !
                IF ( gstart == 2 ) ps(n) = ps(n) - all_psi(1,1) * all_spsi(1,n)
                !
             END DO
             !
             CALL reduce( nbase, ps )
             !
             DO n = 2, nbase
                !
                all_psi(:,1)  = all_psi(:,1)  - ps(n) * all_psi(:,n)
                all_hpsi(:,1) = all_hpsi(:,1) - ps(n) * all_hpsi(:,n)
                all_spsi(:,1) = all_spsi(:,1) - ps(n) * all_spsi(:,n)
                !
             END DO
             !
             psiSpsi = 2.D0 * DDOT( ndim2, all_psi(1,1), 1, all_spsi(1,1), 1 )
             !
             IF ( gstart == 2 ) psiSpsi = psiSpsi - all_psi(1,1) * all_spsi(1,1)
             !
             CALL reduce( 1, psiSpsi )
             !
             IF ( psiSpsi < eps32 ) THEN
                !
                PRINT *, psiSpsi
                !
                CALL errore( 'diis_step', ' negative norm in S ', 1 )
                !
             END IF
             !
             psiSpsi = 1.D0 / SQRT( psiSpsi )
             ! 
             all_psi(:,1)  = psiSpsi * all_psi(:,1)
             all_hpsi(:,1) = psiSpsi * all_hpsi(:,1)
             all_spsi(:,1) = psiSpsi * all_spsi(:,1)
             !
             ! ... the enrgy of the new wavefunction is computed here
             !
             e(ib) = 2.D0 * DDOT( ndim2, psi(1,ib), 1, hpsi(1,ib), 1 )
             !
             IF ( gstart == 2 ) e(ib) = e(ib) - psi(1,ib) * hpsi(1,ib)
             !
             CALL reduce( 1, e(ib) )
             !
             ! ... the "energy" history of this band is reconstructed
             !
             e_small(1)       = e(ib)
             e_small(2:nbase) = e_old(1:dim,ib)
             !
             ! ... DIIS work-space is refreshed
             !
             psi_old(:,1:dim_new,ib)  = all_psi(:,1:dim_new)
             hpsi_old(:,1:dim_new,ib) = all_hpsi(:,1:dim_new)
             spsi_old(:,1:dim_new,ib) = all_spsi(:,1:dim_new)
             e_old(1:dim_new,ib)      = e_small(1:dim_new)
             !
             ! ... residual vectors are finally computed
             !
             FORALL( n = 1: nbase ) 
                !
                all_res(:,n) = ( all_hpsi(:,n) - e_small(n) * all_spsi(:,n) )
                !
             END FORALL             
             !
             ! ... here we construct the matrices :
             ! ...    rr_ij = <res_i|res_j>  and  sr_ij = <psi_i|S|psi_j>
             !
             CALL DGEMM( 'T', 'N', nbase, nbase, ndim2, 2.D0, all_res, &
                         ndmx2, all_res, ndmx2, 0.D0, rr_small, nbase )
             !
             IF ( gstart == 2 ) &
                CALL DGER( nbase, nbase, -1.D0, all_res, &
                           ndmx2, all_res, ndmx2, rr_small, nbase )
             !
             CALL reduce( nbase * nbase, rr_small )
             !
             sr_small(:,:) = 0.D0
             !
             FORALL( n = 1: nbase ) sr_small(n,n) = 1.D0
             ! 
             ! ... diagonalize the reduced hamiltonian
             !
             CALL rdiaghg( nbase, 1, rr_small, &
                           sr_small, nbase, e_small, vr_small )
             !
             ! ... here we compute the best estimate of the |psi>, H|psi>
             ! ... and S|psi>
             !
             CALL DGEMM( 'N', 'N', ndim2, 1, nbase, 1.D0, all_psi(1,1), &
                         ndmx2, vr_small(1,1), nbase, 0.D0, psi(1,ib), ndmx2 )
             !
             CALL DGEMM( 'N', 'N', ndim2, 1, nbase, 1.D0, all_hpsi(1,1), &
                         ndmx2, vr_small(1,1), nbase, 0.D0, hpsi(1,ib), ndmx2 )
             !
             CALL DGEMM( 'N', 'N', ndim2, 1, nbase, 1.D0, all_spsi(1,1), &
                         ndmx2, vr_small(1,1), nbase, 0.D0, spsi(1,ib), ndmx2 )
             !
             psiSpsi = 2.D0 * DDOT( ndim2, psi(1,ib), 1, spsi(1,ib), 1 )
             !
             IF ( gstart == 2 ) psiSpsi = psiSpsi - psi(1,ib) * spsi(1,ib)
             !
             CALL reduce( 1, psiSpsi )
             !
             ! ... |psi>,  H|psi>  and  S|psi> are normalized
             !
             psiSpsi = 1.D0 / SQRT( psiSpsi )
             !
             psi(:,ib)  = psi(:,ib)  * psiSpsi
             hpsi(:,ib) = hpsi(:,ib) * psiSpsi
             spsi(:,ib) = spsi(:,ib) * psiSpsi
             !             
             e(ib) = 2.D0 * DDOT( ndim2, psi(1,ib), 1, hpsi(1,ib), 1 )
             !
             IF ( gstart == 2 ) e(ib) = e(ib) - psi(1,ib) * hpsi(1,ib)
             !
             CALL reduce( 1, e(ib) )
             !
          END DO bands_loop             
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
          DEALLOCATE( ps )
          !   
          RETURN
          !
        END SUBROUTINE diis_step
        !
    END SUBROUTINE diis_with_ortho
    !
    !------------------------------------------------------------------------
    SUBROUTINE holes_filler( ndmx, ndim, nbnd, psi, e, precondition, diis_iter )
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : pi, tpi
      USE random_numbers, ONLY : rndm
      !
      IMPLICIT NONE
      !
      ! ... I/O variables 
      !
      INTEGER,           INTENT(IN)    :: ndmx, ndim, nbnd
      REAL (DP),    INTENT(IN)    :: precondition(ndim)
      COMPLEX (DP), INTENT(INOUT) :: psi(ndmx,nbnd)
      REAL (DP),    INTENT(INOUT) :: e(nbnd)
      INTEGER,           INTENT(INOUT) :: diis_iter
      !
      ! ... local variables
      !
      INTEGER                        :: i, j, ib
      INTEGER                        :: cgter, holes_filler_iter
      REAL (DP),    ALLOCATABLE :: lagrange(:)
      COMPLEX (DP), ALLOCATABLE :: g(:), cg(:), scg(:), ppsi(:), g0(:)
      REAL (DP)                 :: psi_norm, a0, b0, gg0, gamma, gg, &
                                        gg1, theta, cg0, e0, es(2)
      REAL (DP)                 :: rr, arg
      !
      !
      ALLOCATE( scg(  ndmx ) )    
      ALLOCATE( g(    ndmx ) )    
      ALLOCATE( cg(   ndmx ) )    
      ALLOCATE( g0(   ndmx ) )    
      ALLOCATE( ppsi( ndmx ) )    
      ALLOCATE( lagrange( nbnd ) )    
      !
      holes_filler_iter = 0
      !
      ! ... every eigenfunction is calculated separately
      !
      DO ib = ( nbnd_diis + 1 ), nbnd
         !
#if defined (DIIS_DEBUG)         
         PRINT *, "CG: BAND = ", ib
#endif
         !
         ! ... add some noise to |psi>
         !
         rr = DDOT( ndim2, aux(1,ib), 1, aux(1,ib), 1 )
         !
         rr = SQRT( rr )
         !
         DO i = 1, ndim
            !
            arg = tpi * rndm()
            !
            psi(i,ib) = psi(i,ib) + rr * CMPLX( COS( arg ), SIN( arg ) )
            !
         END DO
         !
         IF ( gstart == 2 ) psi(1,ib) = CMPLX( DBLE( psi(1,ib) ), 0.D0 )
         !
         ! ... calculate S|psi>
         !        
         CALL s_1psi( ndmx, ndim, psi(1,ib), spsi(1,ib) )
         !
         ! ... orthogonalize starting eigenfunction to those already calculated
         !
         DO j = 1, ( ib - 1 )
            !
            lagrange(j) = 2.D0 * DDOT( ndim2, psi(1,j), 1, spsi(1,ib), 1 )
            !
            IF ( gstart == 2 ) lagrange(j) = lagrange(j) - psi(1,j) * spsi(1,ib)
            !
         END DO
         !
         CALL reduce( ( ib - 1 ), lagrange )
         !
         DO j = 1, ( ib - 1 )
            !
            psi(:,ib)  = psi(:,ib)  - lagrange(j) * psi(:,j)
            spsi(:,ib) = spsi(:,ib) - lagrange(j) * spsi(:,j)
            !
         END DO
         !
         psi_norm = 2.D0 * DDOT( ndim2, psi(1,ib), 1, spsi(1,ib), 1 )
         !
         IF ( gstart == 2 ) psi_norm = psi_norm - psi(1,ib) * spsi(1,ib)
         !
         CALL reduce( 1, psi_norm )
         !
         IF ( psi_norm < eps32 ) THEN
            !
            PRINT *, psi_norm
            !
            CALL errore( 'holes_filler', ' negative norm in S ', 1 )
            !
         END IF
         !
         psi_norm = 1.D0 / SQRT( psi_norm )
         !
         psi(:,ib)  = psi(:,ib)  * psi_norm
         spsi(:,ib) = spsi(:,ib) * psi_norm
         !
         ! ... calculate starting gradient ( |hpsi> = H|psi> ) ...
         !
         CALL h_psi( ndmx, ndim, 1, psi(1,ib), hpsi(1,ib) )
         !
         ! ... and starting eigenvalue ( e = <y|PHP|y> = <psi|H|psi> )
         !
         e(ib) = 2.D0 * DDOT( ndim2, psi(1,ib), 1, hpsi(1,ib), 1 )
         !
         IF ( gstart == 2 ) e(ib) = e(ib) - psi(1,ib) * hpsi(1,ib)
         !
         CALL reduce( 1, e(ib) )
         !
#if defined (DIIS_DEBUG)                     
         PRINT *, "INITIAL EIGENVALUE = ", e(ib)
#endif
         !
         ! ... start iteration for this band
         !
         iterate: DO cgter = 1, cg_maxter
            !
#if defined (DIIS_DEBUG)                     
            PRINT *, "ITERATION = ", cgter
#endif            
            !
            ! ... calculate  P (PHP)|y>
            ! ... ( P = preconditioning matrix, assumed diagonal )
            !
            g(1:ndim)    = hpsi(1:ndim,ib) / precondition(1:ndim)
            ppsi(1:ndim) = spsi(1:ndim,ib) / precondition(1:ndim)
            !
            ! ... ppsi is now S P(P^2)|y> = S P^2|psi>)
            !
            es(1) = 2.D0 * DDOT( ndim2, spsi(1,ib), 1, g(1), 1 )
            es(2) = 2.D0 * DDOT( ndim2, spsi(1,ib), 1, ppsi(1), 1 )
            !
            IF ( gstart == 2 ) THEN
               !
               es(1) = es(1) - spsi(1,ib) * g(1)
               es(2) = es(2) - spsi(1,ib) * ppsi(1)
               !
            END IF
            !
            CALL reduce( 2, es )
            !
            es(1) = es(1) / es(2)
            !
            g(:) = g(:) - es(1) * ppsi(:)
            !
            ! ... e1 = <y| S P^2 PHP|y> / <y| S S P^2|y> ensures that 
            ! ... <g| S P^2|y> = 0
            ! ... orthogonalize to lowest eigenfunctions (already calculated)
            !
            ! ... scg is used as workspace
            !
            CALL s_1psi( ndmx, ndim, g(1), scg(1) )
            !
            DO j = 1, ( ib - 1 )
               !
               lagrange(j) = 2.D0 * DDOT( ndim2, psi(1,j), 1, scg(1), 1 )
               !
               IF ( gstart == 2 ) lagrange(j) = lagrange(j) - psi(1,j) * scg(1)
               !
            END DO
            !
            CALL reduce( ( ib - 1 ), lagrange )
            !
            DO j = 1, ( ib - 1 )
               !
               g(:)   = g(:)   - lagrange(j) * psi(:,j)
               scg(:) = scg(:) - lagrange(j) * psi(:,j)
               !
            END DO
            !
            IF ( cgter /= 1 ) THEN
               !
               ! ... gg1 is <g(n+1)|S|g(n)> (used in Polak-Ribiere formula)
               !
               gg1 = 2.D0 * DDOT( ndim2, g(1), 1, g0(1), 1 )
               !
               IF ( gstart == 2 ) gg1 = gg1 - g(1) * g0(1)
               !
               CALL reduce( 1, gg1 )
               !
            END IF
            !
            ! ... gg is <g(n+1)|S|g(n+1)>
            !
            g0(:) = scg(:)
            !
            g0(1:ndim) = g0(1:ndim) * precondition(1:ndim)
            !
            gg = 2.D0 * DDOT( ndim2, g(1), 1, g0(1), 1 )
            !
            IF ( gstart == 2 ) gg = gg - g(1) * g0(1)
            !
            CALL reduce( 1, gg )
            !
            IF ( cgter == 1 ) THEN
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
               psi_norm = gamma * cg0 * SIN( theta )
               !
               cg(:) = cg(:) - psi_norm * psi(:,ib)
               !
            END IF
            !
            ! ... |cg> contains now the conjugate gradient
            !
            ! ... |scg> is S|cg>
            !
            CALL h_1psi( ndmx, ndim, cg(1), ppsi(1), scg(1) )
            !
            cg0 = 2.D0 * DDOT( ndim2, cg(1), 1, scg(1), 1 )
            !
            IF ( gstart == 2 ) cg0 = cg0 - cg(1) * scg(1)
            !
            CALL reduce( 1, cg0 )
            !
            cg0 = SQRT( cg0 )
            !
            ! ... |ppsi> contains now HP|cg>
            ! ... minimize <y(t)|PHP|y(t)> , where :
            ! ...                          |y(t)> = cos(t)|y> + sin(t)/cg0 |cg>
            ! ... Note that  <y|P^2S|y> = 1, <y|P^2S|cg> = 0 ,
            ! ...           <cg|P^2S|cg> = cg0^2
            ! ... so that the result is correctly normalized :
            ! ...                           <y(t)|P^2S|y(t)> = 1
            !
            a0 = 4.D0 * DDOT( ndim2, psi(1,ib), 1, ppsi(1), 1 )
            !
            IF ( gstart == 2 ) a0 = a0 - psi(1,ib) * ppsi(1)
            !
            a0 = a0 / cg0
            !
            CALL reduce( 1, a0 )
            !
            b0 = 2.D0 * DDOT( ndim2, cg(1), 1, ppsi(1), 1 )
            !
            IF ( gstart == 2 ) b0 = b0 - cg(1) * ppsi(1)
            !
            b0 = b0 / cg0**2
            !
            CALL reduce( 1, b0 )
            !
            e0 = e(ib)
            !
            theta = 0.5D0 * ATAN( a0 / ( e0 - b0 ) )
            !
            es(1) = 0.5D0 * ( ( e0 - b0 ) * COS( 2.D0 * theta ) + &
                              a0 * SIN( 2.D0 * theta ) + e0 + b0 )       
            es(2) = 0.5D0 * ( - ( e0 - b0 ) * COS( 2.D0 * theta ) - &
                              a0 * SIN( 2.D0 * theta ) + e0 + b0 )
            !
            ! ... there are two possible solutions, choose the minimum
            !
            IF ( es(2) < es(1) ) theta = theta + 0.5D0 * pi
            !
            ! ... new estimate of the eigenvalue
            !
            e(ib) = MIN( es(1), es(2) )
            !
            ! ... upgrade |psi>, H|psi> and S|psi>
            !
            psi(:,ib) = COS( theta ) * psi(:,ib) + &
                        SIN( theta ) / cg0 * cg(:)
            !
            spsi(:,ib) = COS( theta ) * spsi(:,ib) + &
                         SIN( theta ) / cg0 * scg(:)
            !
            hpsi(:,ib) = COS( theta ) * hpsi(:,ib) + &
                         SIN( theta ) / cg0 * ppsi(:)
            !
#if defined (DIIS_DEBUG)            
            PRINT *, e(ib), ABS( e(ib) - e0 )
#endif
            !
            ! ... here one could test convergence on the energy
            !
            IF ( ABS( e(ib) - e0 ) < holes_filler_ethr ) THEN
               !
               conv(ib) = .TRUE.
               !
               EXIT iterate
               !
            END IF   
            !
         END DO iterate
         !
         holes_filler_iter = holes_filler_iter + MIN( cgter, maxter )
         !
      END DO
      !
      diis_iter = diis_iter + ANINT( DBLE( holes_filler_iter ) / DBLE( nbnd ) )
      !
      DEALLOCATE( lagrange )
      DEALLOCATE( ppsi )
      DEALLOCATE( g0 )
      DEALLOCATE( cg )
      DEALLOCATE( g )
      DEALLOCATE( scg )
      !
      RETURN
      !
    END SUBROUTINE holes_filler
    !
END MODULE real_diis_module
