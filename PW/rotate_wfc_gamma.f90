!
! Copyright (C) 2003-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#include "f_defs.h"

!----------------------------------------------------------------------------
SUBROUTINE rotate_wfc_gamma( npwx, npw, nstart, gstart, nbnd, &
                             psi, overlap, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_wfc for Gamma-only calculations
  ! ... This version assumes real wavefunctions (k=0) with only
  ! ... half plane waves stored: psi(-G)=psi*(G), except G=0
  !
  USE kinds,         ONLY : DP
  USE control_flags, ONLY : gamma_only 
  USE mp_global,     ONLY : intra_pool_comm
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER :: npw, npwx, nstart, nbnd, gstart
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! first G with nonzero norm
  LOGICAL :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP) :: psi(npwx,nstart), evc(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e(nbnd)
    ! eigenvalues
  !
  ! ... auxiliary variables:
  !
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  REAL(DP),    ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:), en(:)
  !
  ALLOCATE( aux(  npwx, nstart ) )    
  ALLOCATE( hr( nstart, nstart ) )    
  ALLOCATE( sr( nstart, nstart ) )    
  ALLOCATE( vr( nstart, nstart ) )    
  ALLOCATE( en( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  CALL h_psi( npwx, npw, nstart, psi, aux )
  !
  CALL DGEMM( 'T', 'N', nstart, nstart, 2 * npw, 2.D0 , psi, &
              2 * npwx, aux, 2 * npwx, 0.D0, hr, nstart )
  !  
  IF ( gstart == 2 ) &
     call DGER( nstart, nstart, -1.D0, psi, 2 * npwx, aux, &
                2 * npwx, hr, nstart )
  !     
#if defined (__PARA)
  CALL mp_sum(  hr , intra_pool_comm )
#endif
  !     
  IF ( overlap ) THEN
     !
     CALL s_psi( npwx, npw, nstart, psi, aux )
     !
     CALL DGEMM( 'T', 'N', nstart, nstart, 2 * npw, 2.D0 , psi, &
                 2 * npwx, aux, 2 * npwx, 0.D0, sr, nstart )
     !            
     IF ( gstart == 2 ) &
        CALL DGER( nstart, nstart, -1.D0, psi, 2 * npwx, &
                   aux, 2 * npwx, sr, nstart )
     !              
  ELSE
     !
     CALL DGEMM( 'T', 'N', nstart, nstart, 2 * npw, 2.D0, psi, &
                 2 * npwx, psi, 2 * npwx, 0.D0, sr, nstart )
     !
     IF ( gstart == 2 ) &
        CALL DGER( nstart, nstart, -1.D0, psi, 2 * npwx, &
                   psi, 2 * npwx, sr, nstart )
     !
  END IF
  !
#if defined (__PARA)
  CALL mp_sum(  sr , intra_pool_comm )
#endif
  !
  ! ... Diagonalize
  !
  CALL rdiaghg( nstart, nbnd, hr, sr, nstart, en, vr )
  !
  e(:) = en(1:nbnd)
  !
  ! ... update the basis set
  !
  CALL DGEMM( 'N', 'N', 2 * npw, nbnd, nstart, 1.D0, psi, 2 * npwx, &
              vr, nstart, 0.D0, aux, 2 * npwx )
  !   
  evc(:,:) = aux(:,1:nbnd)
  !
  DEALLOCATE( en )
  DEALLOCATE( vr )
  DEALLOCATE( sr )
  DEALLOCATE( hr )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE rotate_wfc_gamma
!
!
!----------------------------------------------------------------------------
SUBROUTINE protate_wfc_gamma( npwx, npw, nstart, gstart, nbnd, psi, overlap, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Parallel version of rotate_wfc for Gamma-only calculations
  ! ... Subroutine with distributed matrices, written by Carlo Cavazzoni
  ! ... This version assumes real wavefunctions (k=0) with only
  ! ... half plane waves stored: psi(-G)=psi*(G), except G=0
  !
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : gamma_only 
  USE mp_global,        ONLY : npool, nproc_pool, me_pool, root_pool, &
                               intra_pool_comm, init_ortho_group, me_image, &
                               ortho_comm, np_ortho, me_ortho, ortho_comm_id,&
                               leg_ortho
  USE descriptors,      ONLY : descla_siz_ , descla_init , lambda_node_ , &
                               la_nx_ , la_nrl_ , la_n_ , &
                               ilac_ , ilar_ , nlar_ , nlac_ , la_npc_ , &
                               la_npr_ , la_me_ , la_comm_ , &
                               la_myr_ , la_myc_ , nlax_
  USE parallel_toolkit, ONLY : dsqmred, dsqmdst, dsqmsym
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier, mp_end

  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER :: npw, npwx, nstart, nbnd, gstart
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! first G with nonzero norm
  LOGICAL :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP) :: psi(npwx,nstart), evc(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e(nbnd)
    ! eigenvalues
  !
  ! ... auxiliary variables:
  !
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  REAL(DP),    ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:), en(:)
  !
  INTEGER :: desc( descla_siz_ )
    ! matrix distribution descriptors
  INTEGER :: nx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  INTEGER, ALLOCATABLE :: desc_ip( :, :, : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
  !
  ALLOCATE( desc_ip( descla_siz_ , np_ortho(1), np_ortho(2) ) )
  ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ) )
  !
  CALL desc_init( nstart, desc, desc_ip )
  !
  ALLOCATE( aux(  npwx, nstart ) )    
  ALLOCATE( hr( nx, nx ) )    
  ALLOCATE( sr( nx, nx ) )    
  ALLOCATE( vr( nx, nx ) )    
  ALLOCATE( en( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  CALL h_psi( npwx, npw, nstart, psi, aux )
  !
  CALL compute_distmat( hr, psi, aux )
  !
  IF ( overlap ) THEN
     !
     CALL s_psi( npwx, npw, nstart, psi, aux )
     !
     CALL compute_distmat( sr, psi, aux )
     !              
  ELSE
     !
     CALL compute_distmat( sr, psi, psi )
     !
  END IF
  !
  ! ... Diagonalize
  !
  CALL prdiaghg( nstart, hr, sr, nx, en, vr, desc )
  !
  e(:) = en(1:nbnd)
  !
  ! ... update the basis set
  !
  CALL refresh_evc( )
  !   
  evc(:,:) = aux(:,1:nbnd)
  !
  DEALLOCATE( desc_ip )
  DEALLOCATE( rank_ip )
  DEALLOCATE( en )
  DEALLOCATE( vr )
  DEALLOCATE( sr )
  DEALLOCATE( hr )
  DEALLOCATE( aux )
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE desc_init( nsiz, desc, desc_ip )
     !
     INTEGER, INTENT(IN)  :: nsiz
     INTEGER, INTENT(OUT) :: desc(:)
     INTEGER, INTENT(OUT) :: desc_ip(:,:,:)
     INTEGER :: i, j, rank
     INTEGER :: coor_ip( 2 )
     ! 
     CALL descla_init( desc, nsiz, nsiz, np_ortho, me_ortho, ortho_comm, ortho_comm_id )
     ! 
     nx = desc( nlax_ )
     !
     DO j = 0, desc( la_npc_ ) - 1
        DO i = 0, desc( la_npr_ ) - 1
           coor_ip( 1 ) = i
           coor_ip( 2 ) = j
           CALL descla_init( desc_ip(:,i+1,j+1), desc( la_n_ ), desc( la_nx_ ), &
                             np_ortho, coor_ip, ortho_comm, 1 )
           CALL GRID2D_RANK( 'R', desc( la_npr_ ), desc( la_npc_ ), i, j, rank )
           rank_ip( i+1, j+1 ) = rank * leg_ortho
        END DO
     END DO
     !
     la_proc = .FALSE.
     IF( desc( lambda_node_ ) > 0 ) la_proc = .TRUE.
     !
     RETURN
  END SUBROUTINE desc_init
  !
  !
  SUBROUTINE compute_distmat( dm, v, w )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm 
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     REAL(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = 0.0d0
     !
     DO ipc = 1, desc( la_npc_ ) !  loop on column procs 
        !
        nc = desc_ip( nlac_ , 1, ipc )
        ic = desc_ip( ilac_ , 1, ipc )
        !
        DO ipr = 1, ipc ! use symmetry for the loop on row procs
           !
           nr = desc_ip( nlar_ , ipr, ipc )
           ir = desc_ip( ilar_ , ipr, ipc )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL DGEMM( 'T', 'N', nr, nc, 2*npw, 2.D0 , &
                       v(1,ir), 2*npwx, w(1,ic), 2*npwx, 0.D0, work, nx )

           IF ( gstart == 2 ) &
              CALL DGER( nr, nc, -1.D0, v(1,ir), 2*npwx, w(1,ic), 2*npwx, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, intra_pool_comm )

        END DO
        !
     END DO
     !
     CALL dsqmsym( nstart, dm, nx, desc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat
  !
  !
  SUBROUTINE refresh_evc( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     REAL(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc( la_npc_ )
        !
        nc = desc_ip( nlac_ , 1, ipc )
        ic = desc_ip( ilac_ , 1, ipc )
        !
        IF( ic <= nbnd ) THEN
           !
           nc = min( nc, nbnd - ic + 1 )
           !
           beta = 0.0d0

           DO ipr = 1, desc( la_npr_ )
              !
              nr = desc_ip( nlar_ , ipr, ipc )
              ir = desc_ip( ilar_ , ipr, ipc )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == desc( la_myr_ ) .AND. ipc-1 == desc( la_myc_ ) .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vr(:,1:nc), root, intra_pool_comm )
                 CALL DGEMM( 'N', 'N', 2*npw, nc, nr, 1.D0, &
                          psi(1,ir), 2*npwx, vr, nx, beta, aux(1,ic), 2*npwx )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, intra_pool_comm )
                 CALL DGEMM( 'N', 'N', 2*npw, nc, nr, 1.D0, &
                          psi(1,ir), 2*npwx, vtmp, nx, beta, aux(1,ic), 2*npwx )
              END IF
              ! 

              beta = 1.0d0

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )

     RETURN
  END SUBROUTINE refresh_evc
  !
END SUBROUTINE protate_wfc_gamma
