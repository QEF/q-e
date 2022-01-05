!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0._DP, 0._DP )
#define ONE  ( 1._DP, 0._DP )
#define MONE (-1._DP, 0._DP )
!
!--------------------------------------------------------------------------
SUBROUTINE gram_schmidt_k_gpu( npwx, npw, nbnd, npol, psi_d, hpsi_d, spsi_d, e, &
                           uspp, eigen, reorder, nbsize )
  !--------------------------------------------------------------------------
  !
  ! ... Gram-Schmidt orthogonalization, for k-point calculations.
  ! ... blocking algorithm is used.
  !
  USE util_param,     ONLY : DP, eps16
  USE mp,            ONLY : mp_sum, mp_max, mp_bcast
  USE mp_bands_util, ONLY : inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id
  USE device_fbuff_m,         ONLY : buffer => dev_buf
  USE device_memcpy_m,        ONLY : dev_memcpy, dev_memset
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npw, npwx, nbnd, npol
  COMPLEX(DP), INTENT(INOUT) :: psi_d (npwx*npol,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: hpsi_d(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: spsi_d(npwx*npol,nbnd)
  REAL(DP),    INTENT(INOUT) :: e(nbnd)
  LOGICAL,     INTENT(IN)    :: uspp
  LOGICAL,     INTENT(IN)    :: eigen
  LOGICAL,     INTENT(IN)    :: reorder
  INTEGER,     INTENT(IN)    :: nbsize
  !
  ! ... local variables
  !
  LOGICAL                  :: eigen_
  INTEGER                  :: kdim, kdmx
  INTEGER                  :: iblock, nblock
  INTEGER                  :: iblock_start, iblock_end
  INTEGER                  :: jblock_start, jblock_end
  INTEGER                  :: ibnd_start, ibnd_end
  INTEGER                  :: jbnd_start, jbnd_end
  INTEGER,     ALLOCATABLE :: owner_bgrp_id(:)
  INTEGER                  :: buf_start, buf_end, buf_size
  !
  ! ... device variables 
  !
  INTEGER :: ii, jj, kk, info
  COMPLEX(DP), ALLOCATABLE :: phi_d(:,:), hphi_d(:,:), sphi_d(:,:)
#if defined (__CUDA)
  attributes(device) :: psi_d, hpsi_d, spsi_d
  attributes(device) :: phi_d, hphi_d, sphi_d
#endif
  !
  COMPLEX(DP), ALLOCATABLE :: sc_d(:), sc2_d(:,:)
#if defined (__CUDA)
  attributes(device) :: sc_d, sc2_d
#endif 
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
  !
  eigen_ = eigen
  !
  IF ( reorder ) THEN
     !
     eigen_ = .TRUE.
     !
  END IF
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
  ALLOCATE( phi_d ( kdmx, nbnd ) )
  IF ( eigen_ ) ALLOCATE( hphi_d( kdmx, nbnd ) )
  IF ( uspp )   ALLOCATE( sphi_d( kdmx, nbnd ) )
  !
  ALLOCATE( owner_bgrp_id( nblock ) )
  !
!$cuf kernel do(2)
  DO ii = 1, kdmx  
    DO jj = 1, nbnd
      phi_d(ii, jj) = ZERO 
    END DO 
  END DO   
  !
  IF ( eigen_ ) THEN 
!$cuf kernel do(2)
    DO ii = 1, kdmx  
      DO jj = 1, nbnd
        hphi_d(ii, jj) = ZERO 
      END DO 
    END DO   
  END IF  
  !
  IF ( uspp )   THEN 
!$cuf kernel do(2)
    DO ii = 1, kdmx  
      DO jj = 1, nbnd
        sphi_d(ii, jj) = ZERO 
      END DO 
    END DO   
  END IF  
  !
  ! ... Set owers of blocks
  !
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
  ! ... Set initial : |phi_j> = |psi_j>
  !
  CALL ZCOPY_gpu( kdmx * nbnd, psi_d(1,1), 1, phi_d(1,1), 1 )
  !
  IF ( eigen_ ) &
  CALL ZCOPY_gpu( kdmx * nbnd, hpsi_d(1,1), 1, hphi_d(1,1), 1 )
  !
  IF ( uspp ) &
  CALL ZCOPY_gpu( kdmx * nbnd, spsi_d(1,1), 1, sphi_d(1,1), 1 )
  !
  !
  ! ... Allocate buffers 
  !
  ALLOCATE (sc_d(nbsize)) 
  IF ( nbnd .GT. nbsize)  ALLOCATE (sc2_d(nbsize, nbnd - nbsize)) 
  !
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
     !
     !
     CALL project_offdiag_gpu( ibnd_start, ibnd_end, jbnd_start, jbnd_end )
     !
     !
  END DO
  !
  ! ... Buffer Realese
  !
  DEALLOCATE (sc_d) 
  IF (nbnd .GT. nbsize) DEALLOCATE (sc2_d) 
  !
  !
  ! ... Copy psi <- phi
  !
  CALL ZCOPY_gpu( kdmx * nbnd, phi_d(1,1), 1, psi_d(1,1), 1 )
  !
  IF ( eigen_ ) &
  CALL ZCOPY_gpu( kdmx * nbnd, hphi_d(1,1), 1, hpsi_d(1,1), 1 )
  !
  IF ( uspp ) &
  CALL ZCOPY_gpu( kdmx * nbnd, sphi_d(1,1), 1, spsi_d(1,1), 1 )
  !
  ! ... Calculate energy eigenvalues
  !
  IF ( eigen_ ) CALL energyeigen_gpu( )
  !
  ! ... Sort wave functions
  !
  IF ( reorder ) CALL sort_vectors_gpu( )
  !
  DEALLOCATE( phi_d )
  IF ( eigen_ ) DEALLOCATE( hphi_d )
  IF ( uspp )   DEALLOCATE( sphi_d )
  !
  DEALLOCATE( owner_bgrp_id )
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
    INTEGER                  :: ibnd
    REAL(DP)                 :: norm
    COMPLEX(DP), EXTERNAL    :: ZDOTC_gpu
    !
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( ibnd > ibnd_start ) THEN
          !
          ! ... <phi_j| S |psi_i>
          !
          IF ( uspp ) THEN
             !
             CALL ZGEMV_gpu( 'C', kdim, ibnd - ibnd_start, ONE, phi_d(1,ibnd_start), kdmx, &
                         spsi_d(1,ibnd), 1, ZERO, sc_d(1), 1 )
             !
          ELSE
             !
             CALL ZGEMV_gpu( 'C', kdim, ibnd - ibnd_start, ONE, phi_d(1,ibnd_start), kdmx, &
                         psi_d(1,ibnd), 1, ZERO, sc_d(1), 1 )
             !
          END IF
          !
          !
          CALL mp_sum( sc_d, intra_bgrp_comm )
          !
          ! ... phi_i = phi_i - phi_j * <phi_j| S |psi_i>
          !
          CALL ZGEMV_gpu( 'N', kdim, ibnd - ibnd_start, MONE, phi_d(1,ibnd_start), kdmx, &
                      sc_d(1), 1, ONE, phi_d(1,ibnd), 1 )
          !
          !
          IF ( eigen_ ) &
          CALL ZGEMV_gpu( 'N', kdim, ibnd - ibnd_start, MONE, hphi_d(1,ibnd_start), kdmx, &
                      sc_d(1), 1, ONE, hphi_d(1,ibnd), 1 )
          !
          IF ( uspp ) &
          CALL ZGEMV_gpu( 'N', kdim, ibnd - ibnd_start, MONE, sphi_d(1,ibnd_start), kdmx, &
                      sc_d(1), 1, ONE, sphi_d(1,ibnd), 1 )
          !
       END IF
       !
       ! ... Normalize : phi_i = phi_i / SQRT(<phi_i| S |phi_i>)
       !
       IF ( uspp ) THEN
          !
          norm = DBLE( ZDOTC_gpu( kdim, phi_d(1,ibnd), 1, sphi_d(1,ibnd), 1 ) )
          !
       ELSE
          !
          norm = DBLE( ZDOTC_gpu( kdim, phi_d(1,ibnd), 1, phi_d(1,ibnd), 1 ) )
          !
       END IF
       !
       CALL mp_sum( norm, intra_bgrp_comm )
       !
       norm = SQRT( MAX( norm, 0._DP ) )
       !
       IF ( norm < eps16 ) &
       CALL errore( ' gram_schmidt_k ', ' vectors are linear dependent ', 1 )
       !
       CALL ZDSCAL_gpu( kdim, 1._DP / norm, phi_d(1,ibnd), 1 )
       !
       IF ( eigen_ ) &
       CALL ZDSCAL_gpu( kdim, 1._DP / norm, hphi_d(1,ibnd), 1 )
       !
       IF ( uspp ) &
       CALL ZDSCAL_gpu( kdim, 1._DP / norm, sphi_d(1,ibnd), 1 )
       !
    END DO
    !
    !
    RETURN
    !
  END SUBROUTINE gram_schmidt_diag
  !
  !
  SUBROUTINE project_offdiag_gpu( ibnd_start, ibnd_end, jbnd_start, jbnd_end )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: ibnd_start, ibnd_end
    INTEGER, INTENT(IN)  :: jbnd_start, jbnd_end
    !
    INTEGER                  :: ibnd_size
    INTEGER                  :: jbnd_size 
    !
    ibnd_size = ibnd_end - ibnd_start + 1
    jbnd_size = jbnd_end - jbnd_start + 1  
    !
    ! ... <phi_i| S |psi_j>
    !
    IF ( uspp ) THEN
       !
       CALL gpu_ZGEMM( 'C', 'N', ibnd_size, jbnd_size, kdim, ONE, phi_d(1,ibnd_start), kdmx, &
                   spsi_d(1,jbnd_start), kdmx, ZERO, sc2_d(1,1), nbsize )
       !
    ELSE
       !
       CALL gpu_ZGEMM( 'C', 'N', ibnd_size, jbnd_size, kdim, ONE, phi_d(1,ibnd_start), kdmx, &
                   psi_d(1,jbnd_start), kdmx, ZERO, sc2_d(1,1), nbsize )
       !
    END IF
    !
    CALL mp_sum( sc2_d, intra_bgrp_comm )
    !
    ! ... phi_j = phi_j - phi_i * <phi_i| S |psi_j>
    !
    CALL gpu_ZGEMM( 'N', 'N', kdim, jbnd_size, ibnd_size, MONE, phi_d(1,ibnd_start), kdmx, &
                sc2_d(1,1), nbsize, ONE, phi_d(1,jbnd_start), kdmx )
    !
    IF ( eigen_ ) &
    CALL gpu_ZGEMM( 'N', 'N', kdim, jbnd_size, ibnd_size, MONE, hphi_d(1,ibnd_start), kdmx, &
                sc2_d(1,1), nbsize, ONE, hphi_d(1,jbnd_start), kdmx )
    !
    IF ( uspp ) &
    CALL gpu_ZGEMM( 'N', 'N', kdim, jbnd_size, ibnd_size, MONE, sphi_d(1,ibnd_start), kdmx, &
                sc2_d(1,1), nbsize, ONE, sphi_d(1,jbnd_start), kdmx )
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
    !
    COMPLEX(DP), EXTERNAL :: ZDOTC_gpu
    !
    ! ... <psi_i| H |psi_i>
    !
    e(:) = 0._DP
    !
    CALL divide( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       e(ibnd) = DBLE( ZDOTC_gpu( kdim, psi_d(1,ibnd), 1, hpsi_d(1,ibnd), 1 ) )
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
          CALL ZSWAP_gpu( kdim, psi_d(1,ibnd), 1, psi_d(1,ibnd-1), 1 )
          !
          IF ( eigen_ ) &
          CALL ZSWAP_gpu( kdim, hpsi_d(1,ibnd), 1, hpsi_d(1,ibnd-1), 1 )
          !
          IF ( uspp ) &
          CALL ZSWAP_gpu( kdim, spsi_d(1,ibnd), 1, spsi_d(1,ibnd-1), 1 )
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
END SUBROUTINE gram_schmidt_k_gpu
