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
SUBROUTINE gram_schmidt_gamma_gpu( npwx, npw, nbnd, psi_d, hpsi_d, spsi_d, e, &
                               uspp, eigen, reorder, nbsize )
  !--------------------------------------------------------------------------
  !
  ! ... Gram-Schmidt orthogonalization, for Gamma-only calculations.
  ! ... blocking algorithm is used.
  !
  USE util_param,     ONLY : DP, eps16
  USE mp,            ONLY : mp_sum, mp_max, mp_bcast
  USE mp_bands_util, ONLY : gstart, inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id
  USE device_memcpy_m,        ONLY : dev_memcpy, dev_memset
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npw, npwx, nbnd
  COMPLEX(DP), INTENT(INOUT) :: psi_d (npwx,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: hpsi_d(npwx,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: spsi_d(npwx,nbnd)
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
  ! ... device variables
  !
  INTEGER     :: ii, buf_start, buf_end, buf_size, info
  COMPLEX(DP),ALLOCATABLE  :: phi_d (:,:)
  COMPLEX(DP),ALLOCATABLE  :: hphi_d(:,:)
  COMPLEX(DP),ALLOCATABLE  :: sphi_d(:,:)
  !
#if defined (__CUDA)
  attributes(device) :: psi_d, hpsi_d, spsi_d
  attributes(device) :: phi_d, hphi_d, sphi_d
#endif 
  !
  COMPLEX(DP), ALLOCATABLE :: sr_d(:), sr2_d(:,:)
#if defined (__CUDA)
  attributes(device) :: sr_d, sr2_d
#endif 
  !
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
  ALLOCATE( phi_d ( npwx, nbnd ) )
  IF ( eigen_ ) ALLOCATE( hphi_d( npwx, nbnd ) )
  IF ( uspp )   ALLOCATE( sphi_d( npwx, nbnd ) )
  !
  phi_d = ZERO
  !
  IF ( eigen_ ) hphi_d = ZERO
  !
  IF ( uspp )   sphi_d = ZERO
!
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
  !
  IF ( gstart == 2 ) THEN
!$cuf kernel do(1)
     DO ii =1,nbnd
        psi_d(1,ii) = CMPLX( DBLE( psi_d(1,ii) ), 0._DP, kind=DP )
     END DO
  END IF
  !
  ! ... Set initial : |phi_j> = |psi_j>
  !
  CALL DCOPY_gpu( npwx2 * nbnd, psi_d(1,1), 1, phi_d(1,1), 1 )
  !
  ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
  !
  IF ( gstart == 2 ) THEN 
!$cuf kernel do(1)
     DO ii =1,nbnd
        phi_d(1,ii) = CMPLX( DBLE( phi_d(1,ii) ), 0._DP, kind=DP )
     END DO
  END IF
  !
  !
  IF ( eigen_ ) THEN
     !
     CALL DCOPY_gpu( npwx2 * nbnd, hpsi_d(1,1), 1, hphi_d(1,1), 1 )
     !
     ! NOTE: set Im[ H*phi(G=0) ] - needed for numerical stability
     IF ( gstart == 2 ) THEN
!$cuf kernel do(1)
        DO ii =1,nbnd
           hphi_d(1,ii) = CMPLX( DBLE( hphi_d(1,ii) ), 0._DP, kind=DP )
        END DO
     END IF
     !
  END IF
  !
  IF ( uspp ) THEN
     !
     CALL DCOPY_gpu( npwx2 * nbnd, spsi_d(1,1), 1, sphi_d(1,1), 1 )
     !
     ! NOTE: set Im[ S*phi(G=0) ] - needed for numerical stability
     IF ( gstart == 2 ) THEN 
!$cuf kernel do(1)
        DO ii =1,nbnd
           sphi_d(1,ii) = CMPLX( DBLE( sphi_d(1,ii) ), 0._DP, kind=DP )
        END DO
     END IF
     !
  END IF
  !
  ! ... Buffer allocation 
  !
  !
  ALLOCATE( sr_d(nbsize))
  IF (nbnd .GT. nbsize) ALLOCATE( sr2_d(nbsize, nbnd - nbsize))
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
     CALL gram_schmidt_diag_gpu( ibnd_start, ibnd_end )
     !
     ! ... Bcast diagonal block
     !
     CALL mp_bcast( phi_d(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
     !
     IF ( eigen_ ) &
     CALL mp_bcast( hphi_d(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
     !
     IF ( uspp ) &
     CALL mp_bcast( sphi_d(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
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
     CALL project_offdiag_gpu( ibnd_start, ibnd_end, jbnd_start, jbnd_end )
     !
  END DO
  !
  ! ... Buffer Realese
  !
  DEALLOCATE (sr_d) 
  IF (nbnd .GT. nbsize) DEALLOCATE(sr2_d) 
  !
  !
  ! ... Copy psi <- phi
  !
  !CALL DCOPY( npwx2 * nbnd, phi(1,1), 1, psi(1,1), 1 )
  CALL DCOPY_gpu( npwx2 * nbnd, phi_d(1,1), 1, psi_d(1,1), 1 )
  !
  IF ( eigen_ ) &
  CALL DCOPY_gpu( npwx2 * nbnd, hphi_d(1,1), 1, hpsi_d(1,1), 1 )
  !CALL DCOPY( npwx2 * nbnd, hphi(1,1), 1, hpsi(1,1), 1 )
  !
  IF ( uspp ) &
  CALL DCOPY_gpu( npwx2 * nbnd, sphi_d(1,1), 1, spsi_d(1,1), 1 )
  !CALL DCOPY( npwx2 * nbnd, sphi(1,1), 1, spsi(1,1), 1 )
  !
  ! ... Calculate energy eigenvalues
  !
 IF ( eigen_ ) CALL energyeigen_gpu( )
  !
  ! ... Sort wave functions
  !
  IF ( reorder ) CALL sort_vectors_gpu( )
  !
  DEALLOCATE( owner_bgrp_id )
  !
  DEALLOCATE( phi_d )
  IF ( eigen_ ) DEALLOCATE( hphi_d )
  IF ( uspp   ) DEALLOCATE( sphi_d )
  !
  RETURN
  !
  !
  CALL stop_clock( 'gsorth' )
  !
  !
CONTAINS
  !
  !
  SUBROUTINE gram_schmidt_diag_gpu( ibnd_start, ibnd_end )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: ibnd_start, ibnd_end
    !
    INTEGER               :: ibnd
    REAL(DP)              :: norm
    REAL(DP)              :: psi_ibnd
    REAL(DP), EXTERNAL    :: gpu_DDOT 
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( ibnd > ibnd_start ) THEN
          !
          ! ... <phi_j| S |psi_i>
          !
          IF ( uspp ) THEN
             !
             CALL DGEMV_gpu( 'T', npw2, ibnd - ibnd_start, 2._DP, phi_d(1,ibnd_start), npwx2, &
                         spsi_d(1,ibnd), 1, 0._DP, sr_d(1), 1 )
             !
             IF ( gstart == 2 ) THEN
                psi_ibnd = -spsi_d(1,ibnd)
                CALL DAXPY_gpu( ibnd - ibnd_start, psi_ibnd , phi_d(1,ibnd_start), npwx2, &
                         sr_d(1), 1 )
             END IF
             !
          ELSE
             !
             CALL DGEMV_gpu( 'T', npw2, ibnd - ibnd_start, 2._DP, phi_d(1,ibnd_start), npwx2, &
                         psi_d(1,ibnd), 1, 0._DP, sr_d(1), 1 )
             !
             IF ( gstart == 2 ) THEN

                psi_ibnd = -psi_d(1,ibnd)
                CALL DAXPY_gpu( ibnd - ibnd_start, psi_ibnd, phi_d(1,ibnd_start), npwx2, &
                            sr_d(1), 1 )
             END IF
             !
          END IF
          !
          CALL mp_sum( sr_d, intra_bgrp_comm )
          !
          ! ... phi_i = phi_i - phi_j * <phi_j| S |psi_i>
          !
          CALL DGEMV_gpu( 'N', npw2, ibnd - ibnd_start, -1._DP, phi_d(1,ibnd_start), npwx2, &
                      sr_d(1), 1, 1._DP, phi_d(1,ibnd), 1 )
          !
          ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
          IF ( gstart == 2 ) THEN
!$cuf kernel do(1)
            DO ii=1,1
               phi_d(1,ibnd) = CMPLX( DBLE( phi_d(1,ibnd) ), 0._DP, kind=DP )
             END DO
          END IF
          !
          IF ( eigen_ ) THEN
             !
             CALL DGEMV_gpu( 'N', npw2, ibnd - ibnd_start, -1._DP, hphi_d(1,ibnd_start), npwx2, &
                         sr_d(1), 1, 1._DP, hphi_d(1,ibnd), 1 )
             !
             ! NOTE: set Im[ H*phi(G=0) ] - needed for numerical stability
             IF ( gstart == 2 ) THEN                
!$cuf kernel do(1)
                DO ii=1,1
                   hphi_d(1,ibnd) = CMPLX( DBLE( hphi_d(1,ibnd) ), 0._DP, kind=DP )
                END DO
             END IF
             !
          END IF
          !
          IF ( uspp ) THEN
             !
             CALL DGEMV_gpu( 'N', npw2, ibnd - ibnd_start, -1._DP, sphi_d(1,ibnd_start), npwx2, &
                         sr_d(1), 1, 1._DP, sphi_d(1,ibnd), 1 )
             !
             ! NOTE: set Im[ S*phi(G=0) ] - needed for numerical stability
             IF ( gstart == 2 ) THEN
!$cuf kernel do(1)
                DO ii=1,1
                   sphi_d(1,ibnd) = CMPLX( DBLE( sphi_d(1,ibnd) ), 0._DP, kind=DP )
                END DO
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
          norm = 2._DP * gpu_DDOT( npw2, phi_d(1,ibnd), 1, sphi_d(1,ibnd), 1 )
          !
          IF ( gstart == 2 ) THEN
!$cuf kernel do(1)
             DO ii=1,1          
                norm = norm - DBLE( phi_d(1,ibnd) ) * DBLE ( sphi_d(1,ibnd) )
             END DO
          END IF
          !
       ELSE
          !
          norm = 2._DP * gpu_DDOT( npw2, phi_d(1,ibnd), 1, phi_d(1,ibnd), 1 )
          !
          IF ( gstart == 2 ) THEN
!$cuf kernel do(1)
             DO ii=1,1
                norm = norm - DBLE( phi_d(1,ibnd) ) * DBLE ( phi_d(1,ibnd) )
             END DO
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
       CALL DSCAL_gpu( npw2, 1._DP / norm, phi_d(1,ibnd), 1 )
       !
       ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
       IF ( gstart == 2 ) THEN
!$cuf kernel do(1)
             DO ii=1,1          
                phi_d(1,ibnd) = CMPLX( DBLE( phi_d(1,ibnd) ), 0._DP, kind=DP )
             END DO
       END IF
       !
       IF ( eigen_ ) THEN
          !
          CALL DSCAL_gpu( npw2, 1._DP / norm, hphi_d(1,ibnd), 1 )
          !
          ! NOTE: set Im[ H*phi(G=0) ] - needed for numerical stability
          IF ( gstart == 2 ) THEN
!$cuf kernel do(1)
              DO ii=1,1          
                 hphi_d(1,ibnd) = CMPLX( DBLE( hphi_d(1,ibnd) ), 0._DP, kind=DP )
              END DO
          END IF
          !
       END IF
       !
       IF ( uspp ) THEN
          !
          CALL DSCAL_gpu( npw2, 1._DP / norm, sphi_d(1,ibnd), 1 )
          !
          ! NOTE: set Im[ S*phi(G=0) ] - needed for numerical stability
          IF ( gstart == 2 ) THEN
!$cuf kernel do(1)
              DO ii=1,1          
                sphi_d(1,ibnd) = CMPLX( DBLE( sphi_d(1,ibnd) ), 0._DP, kind=DP )
              END DO
          END IF
          !
       END IF
       !
    END DO
    !
    RETURN
    !
  END SUBROUTINE gram_schmidt_diag_gpu
  !
  !
  SUBROUTINE project_offdiag_gpu( ibnd_start, ibnd_end, jbnd_start, jbnd_end )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: ibnd_start, ibnd_end
    INTEGER, INTENT(IN)  :: jbnd_start, jbnd_end
    !
    INTEGER               :: ibnd_size
    INTEGER               :: jbnd_size
    !
    !
    ibnd_size = ibnd_end - ibnd_start + 1
    jbnd_size = jbnd_end - jbnd_start + 1
    !
    ! ... <phi_i| S |psi_j>
    !
    IF ( uspp ) THEN
       !
       CALL gpu_DGEMM( 'T', 'N', ibnd_size, jbnd_size, npw2, 2._DP, phi_d(1,ibnd_start), npwx2, &
                   spsi_d(1,jbnd_start), npwx2, 0._DP, sr2_d(1,1), nbsize )
       !
       IF ( gstart == 2 ) &
       CALL gpu_DGER( ibnd_size, jbnd_size, -1._DP, psi_d(1,ibnd_start), npwx2, &
                  spsi_d(1,jbnd_start), npwx2, sr2_d(1,1), nbsize )
       !
    ELSE
       !
       CALL gpu_DGEMM( 'T', 'N', ibnd_size, jbnd_size, npw2, 2._DP, phi_d(1,ibnd_start), npwx2, &
                   psi_d(1,jbnd_start), npwx2, 0._DP, sr2_d(1,1), nbsize )
       !
       IF ( gstart == 2 ) &
       CALL gpu_DGER( ibnd_size, jbnd_size, -1._DP, psi_d(1,ibnd_start), npwx2, &
                  psi_d(1,jbnd_start), npwx2, sr2_d(1,1), nbsize )
       !
    END IF
    !
    CALL mp_sum( sr2_d, intra_bgrp_comm )
    !
    ! ... phi_j = phi_j - phi_i * <phi_i| S |psi_j>
    !
    CALL gpu_DGEMM( 'N', 'N', npw2, jbnd_size, ibnd_size, -1._DP, phi_d(1,ibnd_start), npwx2, &
                sr2_d(1,1), nbsize, 1._DP, phi_d(1,jbnd_start), npwx2 )
    !
    ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
    IF ( gstart == 2 ) THEN
!$cuf kernel do(1)
              DO ii=jbnd_start, jbnd_end          
              phi_d(1, ii) = &
                        CMPLX( DBLE( phi_d(1, ii) ), 0._DP, kind=DP )
              END DO
    END IF
    !
    IF ( eigen_ ) THEN
       !
       CALL gpu_DGEMM( 'N', 'N', npw2, jbnd_size, ibnd_size, -1._DP, hphi_d(1,ibnd_start), npwx2, &
                   sr2_d(1,1), nbsize, 1._DP, hphi_d(1,jbnd_start), npwx2 )
       !
       ! NOTE: set Im[ H*phi(G=0) ] - needed for numerical stability
       IF ( gstart == 2 ) THEN 
!$cuf kernel do(1)
              DO ii= jbnd_start, jbnd_end         
                 hphi_d(1, ii) = &
                          CMPLX( DBLE( hphi_d(1, ii) ), 0._DP, kind=DP )
              END DO
       END IF
       !
    END IF
    !
    IF ( uspp ) THEN
       !
       CALL gpu_DGEMM( 'N', 'N', npw2, jbnd_size, ibnd_size, -1._DP, sphi_d(1,ibnd_start), npwx2, &
                   sr2_d(1,1), nbsize, 1._DP, sphi_d(1,jbnd_start), npwx2 )
       !
       ! NOTE: set Im[ S*phi(G=0) ] - needed for numerical stability
       IF ( gstart == 2 ) THEN 
!$cuf kernel do(1)
              DO ii=jbnd_start, jbnd_end
                  sphi_d(1, ii ) = &
                          CMPLX( DBLE( sphi_d(1, ii) ), 0._DP, kind=DP )
              END DO
       END IF
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE project_offdiag_gpu
  !
  !
  SUBROUTINE energyeigen_gpu( )
    !
    IMPLICIT NONE
    !
    INTEGER :: ibnd, ibnd_start, ibnd_end
    REAL(DP) :: e_ibnd
    !
    REAL(DP), EXTERNAL :: gpu_DDOT
    !
    ! ... <psi_i| H |psi_i>
    !
    e(:) = 0._DP
    !
    CALL divide( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       e(ibnd) = 2._DP * gpu_DDOT( npw2, psi_d(1,ibnd), 1, hpsi_d(1,ibnd), 1 )
       !
       IF ( gstart == 2 ) THEN
!$cuf kernel do(1)
          DO ii=1,1
             e_ibnd = DBLE( psi_d(1,ibnd) ) * DBLE ( hpsi_d(1,ibnd) )
          END DO
          e(ibnd) = e(ibnd) - e_ibnd
       END IF
       !
    END DO
    !
    CALL mp_sum( e(ibnd_start:ibnd_end), intra_bgrp_comm )
    CALL mp_sum( e, inter_bgrp_comm )
    !
    RETURN
    !
  END SUBROUTINE energyeigen_gpu
  !
  !
  SUBROUTINE sort_vectors_gpu( )
    !
    IMPLICIT NONE
    !
    INTEGER  :: ibnd
    INTEGER  :: nswap
    REAL(DP) :: e0
    !
10  nswap = 0
    !
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
          CALL DSWAP_gpu( npw2, psi_d(1,ibnd), 1, psi_d(1,ibnd-1), 1 )
          !
          IF ( eigen_ ) &
          CALL DSWAP_gpu( npw2, hpsi_d(1,ibnd), 1, hpsi_d(1,ibnd-1), 1 )
          !
          IF ( uspp ) &
          CALL DSWAP_gpu( npw2, spsi_d(1,ibnd), 1, spsi_d(1,ibnd-1), 1 )
          !
       END IF
       !
    END DO
    !
    IF ( nswap > 0 ) GOTO 10
    !
    RETURN
    !
  END SUBROUTINE sort_vectors_gpu
  !
  !
END SUBROUTINE gram_schmidt_gamma_gpu
