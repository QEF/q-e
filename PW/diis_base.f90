!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!#define DIIS_DEBUG
!#define HOLES_FILLER_DEBUG
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
#include "f_defs.h"
! 
! ... Iterative solution of the eigenvalue problem:
!
!     ( H - e S ) |psi> = 0
!
! ... where H is an symmetric operator, e is a real scalar,
! ... S is an overlap matrix, |psi> is a complex vector.
!
! ... The following variant of the DIIS Residual Minimization Method, based
! ... on three steps, is used :
!
! ... 1) At the first "scf" step of the first "ionic" step, the diagonalization
! ...    starts with an initialization procedure. Two (or more) sweeps on all
! ...    the bands are performed. A sweep consists in two (or more if required)
! ...    unconstrained trial steps (steepest descent steps) followed by a
! ...    subspace rotation. If convergence is not achieved in this phase a
! ...    final trial step is done.
! ...    In all the other cases the initialization procedure consists in a 
! ...    subspace rotation followed by a trial step.
!
! ... 2) Diagonalization based on the DIIS algorithm is performed on the lowest
! ...    nbnd_diis bands. Orthogonalization of the eigenvectors is done at
! ...    each step. The possibility of orthogonalizing a given band only to
! ...    those inside an energy window is also implemented.
! 
! ... 3) The topmost "ncgbnd" bands are converged using a standard 
! ...    conjugate-gradient procedure. This ensures that eventual holes
! ...    left by the DIIS algorithm can be identified and filled.
! ...    If "ncgbnd" holes have been found this procedure is repeated on the 
! ...    topmost "ncgbnd" bands (previously optimized with the DIIS algorithm).
!
! ... written by Carlo Sbraccia ( 08/06/2004 )
!  
!
!----------------------------------------------------------------------------
MODULE diis_base
  !----------------------------------------------------------------------------
  !
  USE constants, ONLY : eps32
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE  
  !
  SAVE
  !
  ! ... control parameters
  !
  INTEGER,        PARAMETER :: maxter = 15
    ! maximum number of iterations
  INTEGER,        PARAMETER :: cg_maxter = 100
    ! maximum number of iterations in the cg procedure
  INTEGER,        PARAMETER :: max_sweeps = 3
    ! max number of sweeps
  INTEGER,        PARAMETER :: max_trial_steps = 2
    ! number of trial steps per sweep
  INTEGER,        PARAMETER :: ncgbnd = 2
    ! number of empty bands used by the holes filler
  REAL (DP), PARAMETER :: sweeps_lambda = 2.D0
    ! sweeps step-length  
  REAL (DP), PARAMETER :: diis_lambda = 2.5D0
    ! diis step-length  
  REAL (DP), PARAMETER :: ortho_win_min = 0.03674930873264048097D0 
    ! minimum energy window (Ry) for orthogonalization
  REAL (DP), PARAMETER :: holes_filler_ethr = 1.D-6
    ! convergence threshold for the holes filler
  REAL (DP), PARAMETER :: empty_bands_ethr_min = 1.D-5
    ! minimum convergence threshold on empty bands
  !
  ! ... module variables
  !
  REAL (DP) :: ortho_win
    ! energy window (Ry) for orthogonalization
  REAL (DP) :: empty_ethr 
    ! threshold for empty bands
  COMPLEX (DP), ALLOCATABLE :: hpsi(:,:), spsi(:,:), aux(:,:)
    ! the product of H and psi
    ! the product of S and psi
    ! auxiliary work-space
  COMPLEX (DP), ALLOCATABLE :: psi_old(:,:,:), hpsi_old(:,:,:), &
                                    spsi_old(:,:,:)
    ! DIIS-workspace: old eigenvectors
    ! DIIS-workspace: old product of H and psi
    ! DIIS-workspace: old product of S and psi
  REAL (DP), ALLOCATABLE :: e_old(:,:)
    ! DIIS-workspace: old eigenvalues
  REAL (DP), ALLOCATABLE :: e_ref(:)
    ! eigenvalues of the previous iteration  
  INTEGER :: nbase
    ! counter on the reduced basis vectors
  LOGICAL, ALLOCATABLE :: conv(:)
    ! .TRUE. if the band is converged  
  INTEGER, ALLOCATABLE :: mask(:)
    ! array used to reorder bands
  INTEGER :: cnv, npw2, npwx2, nbnd_diis, diis_npw1
    ! number of converged bands
    ! 2 * npw
    ! 2 * npwx
    ! nbnd - ncgbnd
    ! diis_npw - 1
  !
  CONTAINS
    !
    ! ... allocation/deallocation procedures
    !
    !-----------------------------------------------------------------------
    SUBROUTINE allocate_base( npwx, nbnd )
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: npwx, nbnd
      !
      ALLOCATE( psi_old(  npwx, diis_npw1, nbnd_diis ) )
      ALLOCATE( hpsi_old( npwx, diis_npw1, nbnd_diis ) )
      ALLOCATE( spsi_old( npwx, diis_npw1, nbnd_diis ) )
      ALLOCATE( e_old( diis_npw1, nbnd_diis ) )
      ALLOCATE( e_ref( nbnd ) )
      ALLOCATE( hpsi( npwx, nbnd ) )
      ALLOCATE( spsi( npwx, nbnd ) )
      ALLOCATE( aux(  npwx, nbnd ) )
      ALLOCATE( conv( nbnd ) )
      ALLOCATE( mask( nbnd ) )       
      !
    END SUBROUTINE allocate_base
    !
    !-----------------------------------------------------------------------
    SUBROUTINE deallocate_base()
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      IF ( ALLOCATED( mask ) )     DEALLOCATE( mask )
      IF ( ALLOCATED( conv ) )     DEALLOCATE( conv )
      IF ( ALLOCATED( aux ) )      DEALLOCATE( aux )
      IF ( ALLOCATED( hpsi ) )     DEALLOCATE( hpsi )
      IF ( ALLOCATED( spsi ) )     DEALLOCATE( spsi )
      IF ( ALLOCATED( e_ref ) )    DEALLOCATE( e_ref )
      IF ( ALLOCATED( e_old ) )    DEALLOCATE( e_old )
      IF ( ALLOCATED( psi_old ) )  DEALLOCATE( psi_old )
      IF ( ALLOCATED( hpsi_old ) ) DEALLOCATE( hpsi_old )
      IF ( ALLOCATED( spsi_old ) ) DEALLOCATE( spsi_old )
      !
    END SUBROUTINE deallocate_base
    !
    ! ... other procedures
    !
    !------------------------------------------------------------------------
    FUNCTION no_holes( npwx, nbnd, btype, psi, e )
      !------------------------------------------------------------------------
      !
      USE control_flags, ONLY : ethr
      !
      IMPLICIT NONE
      !
      LOGICAL                          :: no_holes
      INTEGER,           INTENT(IN)    :: npwx, nbnd
      INTEGER,           INTENT(INOUT) :: btype(nbnd)
      REAL (DP),    INTENT(INOUT) :: e(nbnd)
      COMPLEX (DP), INTENT(INOUT) :: psi(npwx,nbnd)
      !
      INTEGER                        :: ib, i, j, moved, btype_tmp
      REAL (DP)                 :: e_tmp, local_ethr
      COMPLEX (DP), ALLOCATABLE :: psi_tmp(:), hpsi_tmp(:), spsi_tmp(:)
      !
      !
      ALLOCATE( psi_tmp(  npwx ) )
      ALLOCATE( hpsi_tmp( npwx ) )
      ALLOCATE( spsi_tmp( npwx ) )
      !
      no_holes = .TRUE.
      !
      moved = 0
      !
#if defined (HOLES_FILLER_DEBUG)
      PRINT *, "OLD SPECTRUM:"
      PRINT *, E(:)
#endif      
      !
      DO ib = 2, nbnd
         !
         ! ... reorder eigenvalues if they are not in the right order
         ! ... ( this can and will happen in not-so-special cases )
         ! ... this procedure should be good if only a few inversions 
         ! ... occur, extremely inefficient if eigenvectors are often 
         ! ... in bad order (but this should not happen)
         !
         IF ( btype(ib) == 1 ) THEN
            !
            local_ethr = 0.5D0 * ethr
            !
         ELSE
            !
            local_ethr = 0.5D0 * empty_ethr
            !
         END IF
         !
         IF ( e(ib) < ( e(ib-1) - local_ethr ) ) THEN
            !
            ! ... if the last calculated eigenvalue is not the largest...
            !
            DO i = ( ib - 2 ), 1, - 1
               !
               IF ( e(ib) > ( e(i) - local_ethr ) ) EXIT
               !
            END DO
            !
            i = i + 1
            !
            moved = moved + 1
            !
            ! ... last calculated eigenvalue should be 
            ! ... in the i-th position: reorder
            !
            e_tmp       = e(ib)
            btype_tmp   = btype(ib) 
            psi_tmp(:)  = psi(:,ib) 
            hpsi_tmp(:) = hpsi(:,ib)
            spsi_tmp(:) = spsi(:,ib) 
            !
            DO j = ib, ( i + 1 ), - 1
               !
               e(j)      = e(j-1)
               btype(j)  = btype(j-1)
               psi(:,j)  = psi(:,j-1)
               hpsi(:,j) = hpsi(:,j-1)
               spsi(:,j) = spsi(:,j-1)
               !
            END DO
            !
            e(i)      = e_tmp
            btype(i)  = btype_tmp
            psi(:,i)  = psi_tmp(:)
            hpsi(:,i) = hpsi_tmp(:)
            spsi(:,i) = spsi_tmp(:)
            !
         END IF
         !
      END DO
      !
      IF ( moved >= ncgbnd ) no_holes = .FALSE.
      !
#if defined (HOLES_FILLER_DEBUG)
      !
      PRINT *, "NEW SPECTRUM:"
      PRINT *, E(:)
      !
      IF ( moved > 0 .AND. moved < ncgbnd ) THEN
         !
         PRINT *, moved, "  HOLES FOUND AND FILLED"
         !
      ELSE IF ( moved == ncgbnd ) THEN
         !
         PRINT *, moved, "  HOLES FOUND AND FILLED"
         PRINT *, "SEARCHING FOR OTHERS..."
         !
      END IF
      !
      PRINT *, "no_holes = ", no_holes
      !
#endif
      !
      DEALLOCATE( psi_tmp )
      DEALLOCATE( hpsi_tmp )
      DEALLOCATE( spsi_tmp )
      !
      RETURN
      !
    END FUNCTION no_holes 
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE reorder_bands( psi, notcnv, nbnd, order )
      !------------------------------------------------------------------------
      !
      ! ... this routine is used to reorder the bands :
      ! ... converged bands come first.
      ! ... for this pourpose an auxiliary vector (aux) is used.
      ! ... the same routine is used to back-order the bands.
      !
      IMPLICIT NONE
      !
      INTEGER,           INTENT(IN)    :: order  
      INTEGER,           INTENT(IN)    :: nbnd
      COMPLEX (DP), INTENT(INOUT) :: psi(:,:)
      INTEGER,           INTENT(INOUT) :: notcnv
      !
      INTEGER :: ib
      !
      !
      IF ( order > 0 ) THEN
         !
         cnv    = 0
         notcnv = nbnd
         !
         ! ... return if no eigenvector is converged
         !
         IF ( .NOT. ANY( conv(1:nbnd) ) ) RETURN
         !
         DO ib = 1, nbnd
            !
            IF ( conv(ib) ) THEN
               !
               cnv = cnv + 1
               !
               mask(ib) = cnv
               !
               aux(:,cnv) = psi(:,ib)
               !
               IF ( cnv == ib ) CYCLE
               !
               hpsi(:,cnv) = hpsi(:,ib)
               spsi(:,cnv) = spsi(:,ib)              
               !
            ELSE
               ! 
               mask(ib) = notcnv
               !
               aux(:,notcnv) = psi(:,ib)
               !
               ! ... notcnv is used as a counter
               !
               notcnv = notcnv - 1
               !
            END IF
            !
         END DO
         !
         psi(:,1:nbnd) = aux(:,1:nbnd)
         !
         ! ... the correct number of not-converged bands is computed here
         !
         notcnv = nbnd - notcnv
         !
      ELSE
         !
         ! ... return if no eigenvector is converged
         !
         IF ( .NOT. ANY( conv(1:nbnd) ) ) RETURN
         !
         aux(:,1:nbnd) = psi(:,mask(1:nbnd))
         !
         psi(:,1:nbnd) = aux(:,1:nbnd)
         !
         aux(:,1:nbnd) = hpsi(:,mask(1:nbnd))
         !
         hpsi(:,1:nbnd) = aux(:,1:nbnd)
         !
         aux(:,1:nbnd) = spsi(:,mask(1:nbnd))
         !
         spsi(:,1:nbnd) = aux(:,1:nbnd)
         !
      END IF
      !
      RETURN      
      !
    END SUBROUTINE reorder_bands         
    !
END MODULE diis_base
