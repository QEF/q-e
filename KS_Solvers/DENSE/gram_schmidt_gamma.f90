!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0._DP, 0._DP )
!
!--------------------------------------------------------------------------
SUBROUTINE gram_schmidt_gamma( npwx, npw, nbnd, psi, hpsi, spsi, e, &
                               uspp, eigen, reorder, nbsize )
  !--------------------------------------------------------------------------
  !
  ! ... Gram-Schmidt orthogonalization, for Gamma-only calculations.
  ! ... blocking algorithm is used.
  !
  USE util_param,     ONLY : DP, eps16
  USE mp,            ONLY : mp_sum, mp_max, mp_bcast
  USE mp_bands_util, ONLY : gstart, inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npw, npwx, nbnd
  COMPLEX(DP), INTENT(INOUT) :: psi (npwx,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: hpsi(npwx,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: spsi(npwx,nbnd)
  REAL(DP),    INTENT(OUT)   :: e(nbnd)
  LOGICAL,     INTENT(IN)    :: uspp
  LOGICAL,     INTENT(IN)    :: eigen
  LOGICAL,     INTENT(IN)    :: reorder
  INTEGER,     INTENT(IN)    :: nbsize
  !
  ! ... local variables
  !
  LOGICAL                  :: eigen_
  INTEGER                  :: npw2, npwx2
  INTEGER                  :: iblock, nblock
  INTEGER                  :: iblock_start, iblock_end
  INTEGER                  :: jblock_start, jblock_end
  INTEGER                  :: ibnd_start, ibnd_end
  INTEGER                  :: jbnd_start, jbnd_end
  COMPLEX(DP), ALLOCATABLE :: phi(:,:), hphi(:,:), sphi(:,:)
  INTEGER,     ALLOCATABLE :: owner_bgrp_id(:)
  !
  COMPLEX(DP), ALLOCATABLE :: sr(:), sr2(:,:)
  !
  CALL start_clock( 'gsorth' )
  !
  eigen_ = eigen
  !
  IF ( reorder ) THEN
     !
     eigen_ = .TRUE.
     !
  END IF
  !
  npw2  = 2 * npw
  npwx2 = 2 * npwx
  !
  nblock = nbnd / nbsize
  IF ( MOD( nbnd, nbsize ) /= 0 ) nblock = nblock + 1
  !
  CALL divide( inter_bgrp_comm, nblock, iblock_start, iblock_end )
  !
  IF ( my_bgrp_id >= nblock ) THEN
     !
     iblock_start = nblock + 1
     iblock_end   = nblock
     !
  END IF
  !
  ALLOCATE( phi ( npwx, nbnd ) )
  IF ( eigen_ ) ALLOCATE( hphi( npwx, nbnd ) )
  IF ( uspp )   ALLOCATE( sphi( npwx, nbnd ) )
  !$acc enter data create( phi, sphi, hphi )
  !
  !$acc kernels
  phi = ZERO
  !
  IF ( eigen_ ) hphi = ZERO
  !
  IF ( uspp )   sphi = ZERO
  !$acc end kernels
  !
  ! ... Set owners of blocks
  !
  ALLOCATE( owner_bgrp_id(nblock)) 
  owner_bgrp_id = 0
  !
  DO iblock = 1, nblock
     !
     IF ( iblock_start <= iblock .AND. iblock <= iblock_end ) &
     owner_bgrp_id(iblock) = my_bgrp_id
     !
  END DO
  !
  CALL mp_max( owner_bgrp_id, inter_bgrp_comm )
  !
  ! ... Set Im[ psi(G=0) ] - needed for numerical stability
  !
  IF ( gstart == 2 ) THEN
    !$acc kernels
    psi(1,1:nbnd) = CMPLX( DBLE( psi(1,1:nbnd) ), 0._DP, kind=DP )
    !$acc end kernels
  END IF
  !
  ! ... Set initial : |phi_j> = |psi_j>
  !
  !$acc host_data use_device(psi, phi)
  CALL MYDCOPY( npwx2 * nbnd, psi(1,1), 1, phi(1,1), 1 )
  !$acc end host_data
  !
  ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
  IF ( gstart == 2 ) THEN
    !$acc kernels
    phi(1,1:nbnd) = CMPLX( DBLE( phi(1,1:nbnd) ), 0._DP, kind=DP )
    !$acc end kernels
  END IF
  !
  IF ( eigen_ ) THEN
     !
     !$acc host_data use_device(hpsi, hphi)
     CALL MYDCOPY( npwx2 * nbnd, hpsi(1,1), 1, hphi(1,1), 1 )
     !$acc end host_data
     !
     ! NOTE: set Im[ H*phi(G=0) ] - needed for numerical stability
     IF ( gstart == 2 ) THEN
       !$acc kernels
       hphi(1,1:nbnd) = CMPLX( DBLE( hphi(1,1:nbnd) ), 0._DP, kind=DP )
       !$acc end kernels
     END IF
     !
  END IF
  !
  IF ( uspp ) THEN
     !
     !$acc host_data use_device( spsi, sphi )
     CALL MYDCOPY( npwx2 * nbnd, spsi(1,1), 1, sphi(1,1), 1 )
     !$acc end host_data
     !
     ! NOTE: set Im[ S*phi(G=0) ] - needed for numerical stability
     IF ( gstart == 2 ) THEN 
       !$acc kernels
       sphi(1,1:nbnd) = CMPLX( DBLE( sphi(1,1:nbnd) ), 0._DP, kind=DP )
       !$acc end kernels
     END IF
     !
  END IF
  !
  ! ... Buffer allocation 
  !
  ALLOCATE( sr(nbsize))
  IF (nbnd .GT. nbsize) ALLOCATE( sr2(nbsize, nbnd - nbsize))
  !$acc enter data create( sr, sr2 )
  !
  ! ... Blocking loop
  !
  DO iblock = 1, nblock
     !
     ! ... Orthogonalize diagonal block by standard Gram-Schmidt
     !
     ibnd_start = ( iblock - 1 ) * nbsize + 1
     ibnd_end   = MIN( iblock * nbsize, nbnd )
     !
     IF ( owner_bgrp_id(iblock) == my_bgrp_id ) &
     CALL gram_schmidt_diag( ibnd_start, ibnd_end )
     !
     ! ... Bcast diagonal block
     !
     !$acc host_data use_device( phi, hphi, sphi )
     CALL mp_bcast( phi(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
     !
     IF ( eigen_ ) &
     CALL mp_bcast( hphi(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
     !
     IF ( uspp ) &
     CALL mp_bcast( sphi(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
     !$acc end host_data
     !
     ! ... Project off-diagonal block outside of diagonal block
     !
     jblock_start = MAX( iblock_start, iblock + 1 )
     jblock_end   = iblock_end
     !
     jbnd_start = ( jblock_start - 1 ) * nbsize + 1
     jbnd_end   = MIN( jblock_end * nbsize, nbnd )
     !
     IF ( jblock_start <= jblock_end .AND. jbnd_start <= jbnd_end ) &
     CALL project_offdiag( ibnd_start, ibnd_end, jbnd_start, jbnd_end )
     !
  END DO
  !
  ! ... Buffer Realese
  !
  !$acc exit data delete( sr, sr2 )
  DEALLOCATE (sr) 
  IF (nbnd .GT. nbsize) DEALLOCATE(sr2) 
  !
  ! ... Copy psi <- phi
  !
  !$acc host_data use_device( phi, psi, sphi, spsi, hphi, hpsi )
  CALL MYDCOPY( npwx2 * nbnd, phi(1,1), 1, psi(1,1), 1 )
  !
  IF ( eigen_ ) &
  CALL MYDCOPY( npwx2 * nbnd, hphi(1,1), 1, hpsi(1,1), 1 )
  !
  IF ( uspp ) &
  CALL MYDCOPY( npwx2 * nbnd, sphi(1,1), 1, spsi(1,1), 1 )
  !$acc end host_data
  !
  ! ... Calculate energy eigenvalues
  !
  IF ( eigen_ ) CALL energyeigen( )
  !
  ! ... Sort wave functions
  !
  IF ( reorder ) CALL sort_vectors( )
  !
  DEALLOCATE( owner_bgrp_id )
  !
  !$acc exit data delete( phi, hphi, sphi)
  DEALLOCATE( phi )
  IF ( eigen_ ) DEALLOCATE( hphi )
  IF ( uspp   ) DEALLOCATE( sphi )
  !
  CALL stop_clock( 'gsorth' )
  !
  RETURN
  !
  !
CONTAINS
  !
  !
  SUBROUTINE gram_schmidt_diag( ibnd_start, ibnd_end )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: ibnd_start, ibnd_end
    !
    INTEGER               :: ibnd
    REAL(DP)              :: norm
    REAL(DP)              :: psi_ibnd
    REAL(DP), EXTERNAL    :: MYDDOT
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( ibnd > ibnd_start ) THEN
          !
          ! ... <phi_j| S |psi_i>
          !
          IF ( uspp ) THEN
             !
             !$acc host_data use_device( phi, spsi, sr )
             CALL MYDGEMV( 'T', npw2, ibnd - ibnd_start, 2._DP, phi(1,ibnd_start), npwx2, &
                         spsi(1,ibnd), 1, 0._DP, sr(1), 1 )
             !$acc end host_data
             !
             IF ( gstart == 2 ) THEN
                !$acc kernels copyout( psi_ibnd )
                psi_ibnd = -spsi(1,ibnd)
                !$acc end kernels
                !$acc host_data use_device( phi, sr )
                CALL MYDAXPY( ibnd - ibnd_start, psi_ibnd , phi(1,ibnd_start), npwx2, &
                         sr(1), 1 )
                !$acc end host_data
             END IF
             !
          ELSE
             !
             !$acc host_data use_device( phi, psi, sr )
             CALL MYDGEMV( 'T', npw2, ibnd - ibnd_start, 2._DP, phi(1,ibnd_start), npwx2, &
                         psi(1,ibnd), 1, 0._DP, sr(1), 1 )
             !$acc end host_data
             !
             IF ( gstart == 2 ) THEN
                !$acc kernels copyout(psi_ibnd)
                psi_ibnd = -psi(1,ibnd)
                !$acc end kernels
                !$acc host_data use_device( phi, sr )
                CALL MYDAXPY( ibnd - ibnd_start, psi_ibnd, phi(1,ibnd_start), npwx2, &
                            sr(1), 1 )
                !$acc end host_data
             END IF
             !
          END IF
          !
          !$acc host_data use_device(sr)
          CALL mp_sum( sr, intra_bgrp_comm )
          !$acc end host_data
          !
          ! ... phi_i = phi_i - phi_j * <phi_j| S |psi_i>
          !
          !$acc host_data use_device(phi, sr)
          CALL MYDGEMV( 'N', npw2, ibnd - ibnd_start, -1._DP, phi(1,ibnd_start), npwx2, &
                      sr(1), 1, 1._DP, phi(1,ibnd), 1 )
          !$acc end host_data
          !
          ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
          IF ( gstart == 2 ) THEN
            !$acc kernels
            phi(1,ibnd) = CMPLX( DBLE( phi(1,ibnd) ), 0._DP, kind=DP )
            !$acc end kernels
          END IF
          !
          IF ( eigen_ ) THEN
             !
             !$acc host_data use_device(hphi, sr)
             CALL MYDGEMV( 'N', npw2, ibnd - ibnd_start, -1._DP, hphi(1,ibnd_start), npwx2, &
                         sr(1), 1, 1._DP, hphi(1,ibnd), 1 )
             !$acc end host_data
             !
             ! NOTE: set Im[ H*phi(G=0) ] - needed for numerical stability
             IF ( gstart == 2 ) THEN
                !$acc kernels
                hphi(1,ibnd) = CMPLX( DBLE( hphi(1,ibnd) ), 0._DP, kind=DP )
                !$acc end kernels
             END IF
             !
          END IF
          !
          IF ( uspp ) THEN
             !
             !$acc host_data use_device(sphi, sr)
             CALL MYDGEMV( 'N', npw2, ibnd - ibnd_start, -1._DP, sphi(1,ibnd_start), npwx2, &
                         sr(1), 1, 1._DP, sphi(1,ibnd), 1 )
             !$acc end host_data
             !
             ! NOTE: set Im[ S*phi(G=0) ] - needed for numerical stability
             IF ( gstart == 2 ) THEN
                !$acc kernels
                sphi(1,ibnd) = CMPLX( DBLE( sphi(1,ibnd) ), 0._DP, kind=DP )
                !$acc end kernels
             END IF
             !
          END IF
          !
       END IF
       !
       ! ... Normalize : phi_i = phi_i / SQRT(<phi_i| S |phi_i>)
       !
       IF ( uspp ) THEN
          !
          !$acc host_data use_device(phi, sphi)
          norm = 2._DP * MYDDOT( npw2, phi(1,ibnd), 1, sphi(1,ibnd), 1 )
          !$acc end host_data
          !
          IF ( gstart == 2 ) THEN
             !$acc kernels copy(norm)
             norm = norm - DBLE( phi(1,ibnd) ) * DBLE ( sphi(1,ibnd) )
             !$acc end kernels 
          END IF
          !
       ELSE
          !
          !$acc host_data use_device(phi)
          norm = 2._DP * MYDDOT( npw2, phi(1,ibnd), 1, phi(1,ibnd), 1 )
          !$acc end host_data
          !
          IF ( gstart == 2 ) THEN
             !$acc kernels copy(norm) 
             norm = norm - DBLE( phi(1,ibnd) ) * DBLE ( phi(1,ibnd) )
             !$acc end kernels
          END IF
          !
       END IF
       !
       CALL mp_sum( norm, intra_bgrp_comm )
       !
       norm = SQRT( MAX( norm, 0._DP ) )
       !
       IF ( norm < eps16 ) &
       CALL errore( ' gram_schmidt_gamma ', ' vectors are linear dependent ', 1 )
       !
       !$acc host_data use_device( phi)
       CALL MYDSCAL( npw2, 1._DP / norm, phi(1,ibnd), 1 )
       !$acc end host_data
       !
       ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
       IF ( gstart == 2 ) THEN
             !$acc kernels
             phi(1,ibnd) = CMPLX( DBLE( phi(1,ibnd) ), 0._DP, kind=DP )
             !$acc end kernels
       END IF
       !
       IF ( eigen_ ) THEN
          !
          !$acc host_data use_device(hphi)
          CALL MYDSCAL( npw2, 1._DP / norm, hphi(1,ibnd), 1 )
          !$acc end host_data
          !
          ! NOTE: set Im[ H*phi(G=0) ] - needed for numerical stability
          IF ( gstart == 2 ) THEN
              !$acc kernels
              hphi(1,ibnd) = CMPLX( DBLE( hphi(1,ibnd) ), 0._DP, kind=DP )
              !$acc end kernels
          END IF
          !
       END IF
       !
       IF ( uspp ) THEN
          !
          !$acc host_data use_device(sphi)
          CALL MYDSCAL( npw2, 1._DP / norm, sphi(1,ibnd), 1 )
          !$acc end host_data
          !
          ! NOTE: set Im[ S*phi(G=0) ] - needed for numerical stability
          IF ( gstart == 2 ) THEN
              !$acc kernels
              sphi(1,ibnd) = CMPLX( DBLE( sphi(1,ibnd) ), 0._DP, kind=DP )
              !$acc end kernels
          END IF
          !
       END IF
       !
    END DO
    !
    RETURN
    !
  END SUBROUTINE gram_schmidt_diag
  !
  !
  SUBROUTINE project_offdiag( ibnd_start, ibnd_end, jbnd_start, jbnd_end )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: ibnd_start, ibnd_end
    INTEGER, INTENT(IN)  :: jbnd_start, jbnd_end
    !
    INTEGER               :: ibnd_size
    INTEGER               :: jbnd_size
    !
    ibnd_size = ibnd_end - ibnd_start + 1
    jbnd_size = jbnd_end - jbnd_start + 1
    !
    ! ... <phi_i| S |psi_j>
    !
    !$acc host_data use_device(phi, spsi, psi, sr2)
    IF ( uspp ) THEN
       !
       CALL MYDGEMM( 'T', 'N', ibnd_size, jbnd_size, npw2, 2._DP, phi(1,ibnd_start), npwx2, &
                   spsi(1,jbnd_start), npwx2, 0._DP, sr2(1,1), nbsize )
       !
       IF ( gstart == 2 ) &
       CALL MYDGER( ibnd_size, jbnd_size, -1._DP, psi(1,ibnd_start), npwx2, &
                  spsi(1,jbnd_start), npwx2, sr2(1,1), nbsize )
       !
    ELSE
       !
       CALL MYDGEMM( 'T', 'N', ibnd_size, jbnd_size, npw2, 2._DP, phi(1,ibnd_start), npwx2, &
                   psi(1,jbnd_start), npwx2, 0._DP, sr2(1,1), nbsize )
       !
       IF ( gstart == 2 ) &
       CALL MYDGER( ibnd_size, jbnd_size, -1._DP, psi(1,ibnd_start), npwx2, &
                  psi(1,jbnd_start), npwx2, sr2(1,1), nbsize )
       !
    END IF
    !
    CALL mp_sum( sr2, intra_bgrp_comm )
    !
    ! ... phi_j = phi_j - phi_i * <phi_i| S |psi_j>
    !
    CALL MYDGEMM( 'N', 'N', npw2, jbnd_size, ibnd_size, -1._DP, phi(1,ibnd_start), npwx2, &
                sr2(1,1), nbsize, 1._DP, phi(1,jbnd_start), npwx2 )
    !$acc end host_data
    !
    ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
    IF ( gstart == 2 ) THEN
              !$acc kernels
              phi(1, jbnd_start:jbnd_end) = &
                        CMPLX( DBLE( phi(1, jbnd_start:jbnd_end ) ), 0._DP, kind=DP )
              !$acc end kernels
    END IF
    !
    IF ( eigen_ ) THEN
       !
       !$acc host_data use_device(hphi, sr2)
       CALL MYDGEMM( 'N', 'N', npw2, jbnd_size, ibnd_size, -1._DP, hphi(1,ibnd_start), npwx2, &
                   sr2(1,1), nbsize, 1._DP, hphi(1,jbnd_start), npwx2 )
       !$acc end host_data
       !
       ! NOTE: set Im[ H*phi(G=0) ] - needed for numerical stability
       IF ( gstart == 2 ) THEN
              !$acc kernels
              hphi(1, jbnd_start:jbnd_end) = &
                          CMPLX( DBLE( hphi(1, jbnd_start:jbnd_end) ), 0._DP, kind=DP )
              !$acc end kernels
       END IF
       !
    END IF
    !
    IF ( uspp ) THEN
       !
       !$acc host_data use_device(sphi, sr2)
       CALL MYDGEMM( 'N', 'N', npw2, jbnd_size, ibnd_size, -1._DP, sphi(1,ibnd_start), npwx2, &
                   sr2(1,1), nbsize, 1._DP, sphi(1,jbnd_start), npwx2 )
       !$acc end host_data
       !
       ! NOTE: set Im[ S*phi(G=0) ] - needed for numerical stability
       IF ( gstart == 2 ) THEN
              !$acc kernels
              sphi(1, jbnd_start:jbnd_end) = &
                          CMPLX( DBLE( sphi(1,jbnd_start:jbnd_end) ), 0._DP, kind=DP )
              !$acc end kernels
       END IF
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE project_offdiag
  !
  !
  SUBROUTINE energyeigen( )
    !
    IMPLICIT NONE
    !
    INTEGER :: ibnd, ibnd_start, ibnd_end
    !
    REAL(DP), EXTERNAL :: MYDDOT_VECTOR_GPU
    !$acc routine( MYDDOT_VECTOR_GPU ) vector
    !
    ! ... <psi_i| H |psi_i>
    !
    !$acc kernels
    e(:) = 0._DP
    !$acc end kernels
    !
    CALL divide( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
    !
    !$acc parallel copyin(npw2, ibnd_start, ibnd_end, gstart) 
    !$acc loop gang 
    DO ibnd = ibnd_start, ibnd_end
       !
       e(ibnd) = 2._DP * MYDDOT_VECTOR_GPU( npw2, psi(1,ibnd), hpsi(1,ibnd) )
       !
       IF ( gstart == 2 ) e(ibnd) = e(ibnd) - DBLE( psi(1,ibnd) ) * DBLE ( hpsi(1,ibnd) )
       !
    END DO
    !$acc end parallel
    !
    !$acc host_data use_device(e)
    CALL mp_sum( e(ibnd_start:ibnd_end), intra_bgrp_comm )
    CALL mp_sum( e, inter_bgrp_comm )
    !$acc end host_data
    !
    RETURN
    !
  END SUBROUTINE energyeigen
  !
  !
  SUBROUTINE sort_vectors( )
    !
    IMPLICIT NONE
    !
    INTEGER  :: ibnd
    INTEGER  :: nswap
    REAL(DP) :: e0
    EXTERNAL :: MYDSWAP_VECTOR_GPU  
    !$acc routine(MYDSWAP_VECTOR_GPU) vector
    !
10  nswap = 0
    !
    !$acc parallel copy(nswap) copyin(npw2)
    !$acc loop gang reduction(+:nswap) private(e0)
    DO ibnd = 2, nbnd
       !
       IF ( e(ibnd) < e(ibnd-1) ) THEN
          !
          nswap = nswap + 1
          !
          e0        = e(ibnd)
          e(ibnd)   = e(ibnd-1)
          e(ibnd-1) = e0
          !
          CALL MYDSWAP_VECTOR_GPU( npw2, psi(1,ibnd),  psi(1,ibnd-1) )
          !
          IF ( eigen_ ) &
          CALL MYDSWAP_VECTOR_GPU( npw2, hpsi(1,ibnd), hpsi(1,ibnd-1))
          !
          IF ( uspp ) &
          CALL MYDSWAP_VECTOR_GPU( npw2, spsi(1,ibnd), spsi(1,ibnd-1))
          !
       END IF
       !
    END DO
    !$acc end parallel
    !
    IF ( nswap > 0 ) GOTO 10
    !
    RETURN
    !
  END SUBROUTINE sort_vectors
  !
  !
END SUBROUTINE gram_schmidt_gamma
