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
SUBROUTINE gram_schmidt_k( npwx, npw, nbnd, npol, psi, hpsi, spsi, e, &
                           uspp, eigen, reorder, nbsize )
  !--------------------------------------------------------------------------
  !
  ! ... Gram-Schmidt orthogonalization, for k-point calculations.
  ! ... blocking algorithm is used.
  !
  USE util_param,     ONLY : DP, eps16
  USE mp,            ONLY : mp_sum, mp_max, mp_bcast
  USE mp_bands_util, ONLY : inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npw, npwx, nbnd, npol
  COMPLEX(DP), INTENT(INOUT) :: psi (npwx*npol,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: hpsi(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: spsi(npwx*npol,nbnd)
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
  COMPLEX(DP), ALLOCATABLE :: phi(:,:), hphi(:,:), sphi(:,:)
  INTEGER,     ALLOCATABLE :: owner_bgrp_id(:)
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
  ALLOCATE( phi ( kdmx, nbnd ) )
  IF ( eigen_ ) ALLOCATE( hphi( kdmx, nbnd ) )
  IF ( uspp )   ALLOCATE( sphi( kdmx, nbnd ) )
  ALLOCATE( owner_bgrp_id( nblock ) )
  !
  phi = ZERO
  !
  IF ( eigen_ ) hphi = ZERO
  !
  IF ( uspp )   sphi = ZERO
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
  CALL ZCOPY( kdmx * nbnd, psi(1,1), 1, phi(1,1), 1 )
  !
  IF ( eigen_ ) &
  CALL ZCOPY( kdmx * nbnd, hpsi(1,1), 1, hphi(1,1), 1 )
  !
  IF ( uspp ) &
  CALL ZCOPY( kdmx * nbnd, spsi(1,1), 1, sphi(1,1), 1 )
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
     CALL mp_bcast( phi(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
     !
     IF ( eigen_ ) &
     CALL mp_bcast( hphi(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
     !
     IF ( uspp ) &
     CALL mp_bcast( sphi(:,ibnd_start:ibnd_end), owner_bgrp_id(iblock), inter_bgrp_comm )
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
  ! ... Copy psi <- phi
  !
  CALL ZCOPY( kdmx * nbnd, phi(1,1), 1, psi(1,1), 1 )
  !
  IF ( eigen_ ) &
  CALL ZCOPY( kdmx * nbnd, hphi(1,1), 1, hpsi(1,1), 1 )
  !
  IF ( uspp ) &
  CALL ZCOPY( kdmx * nbnd, sphi(1,1), 1, spsi(1,1), 1 )
  !
  ! ... Calculate energy eigenvalues
  !
  IF ( eigen_ ) CALL energyeigen( )
  !
  ! ... Sort wave functions
  !
  IF ( reorder ) CALL sort_vectors( )
  !
  DEALLOCATE( phi )
  IF ( eigen_ ) DEALLOCATE( hphi )
  IF ( uspp )   DEALLOCATE( sphi )
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
    COMPLEX(DP), ALLOCATABLE :: sc(:)
    REAL(DP)                 :: norm
    REAL(DP), EXTERNAL       :: DDOT
    !
    ALLOCATE( sc( ibnd_start:ibnd_end ) )
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( ibnd > ibnd_start ) THEN
          !
          ! ... <phi_j| S |psi_i>
          !
          IF ( uspp ) THEN
             !
             CALL ZGEMV( 'C', kdim, ibnd - ibnd_start, ONE, phi(1,ibnd_start), kdmx, &
                         spsi(1,ibnd), 1, ZERO, sc(ibnd_start), 1 )
             !
          ELSE
             !
             CALL ZGEMV( 'C', kdim, ibnd - ibnd_start, ONE, phi(1,ibnd_start), kdmx, &
                         psi(1,ibnd), 1, ZERO, sc(ibnd_start), 1 )
             !
          END IF
          !
          CALL mp_sum( sc, intra_bgrp_comm )
          !
          ! ... phi_i = phi_i - phi_j * <phi_j| S |psi_i>
          !
          CALL ZGEMV( 'N', kdim, ibnd - ibnd_start, MONE, phi(1,ibnd_start), kdmx, &
                      sc(ibnd_start), 1, ONE, phi(1,ibnd), 1 )
          !
          IF ( eigen_ ) &
          CALL ZGEMV( 'N', kdim, ibnd - ibnd_start, MONE, hphi(1,ibnd_start), kdmx, &
                      sc(ibnd_start), 1, ONE, hphi(1,ibnd), 1 )
          !
          IF ( uspp ) &
          CALL ZGEMV( 'N', kdim, ibnd - ibnd_start, MONE, sphi(1,ibnd_start), kdmx, &
                      sc(ibnd_start), 1, ONE, sphi(1,ibnd), 1 )
          !
       END IF
       !
       ! ... Normalize : phi_i = phi_i / SQRT(<phi_i| S |phi_i>)
       !
       IF ( uspp ) THEN
          !
          norm = DDOT( 2*kdim, phi(1,ibnd), 1, sphi(1,ibnd), 1 )
          !
       ELSE
          !
          norm = DDOT ( 2*kdim, phi(1,ibnd), 1, phi(1,ibnd), 1 ) 
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
       CALL ZDSCAL( kdim, 1._DP / norm, phi(1,ibnd), 1 )
       !
       IF ( eigen_ ) &
       CALL ZDSCAL( kdim, 1._DP / norm, hphi(1,ibnd), 1 )
       !
       IF ( uspp ) &
       CALL ZDSCAL( kdim, 1._DP / norm, sphi(1,ibnd), 1 )
       !
    END DO
    !
    DEALLOCATE( sc )
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
    INTEGER                  :: ibnd_size
    INTEGER                  :: jbnd_size
    COMPLEX(DP), ALLOCATABLE :: sc(:,:)
    !
    ibnd_size = ibnd_end - ibnd_start + 1
    jbnd_size = jbnd_end - jbnd_start + 1
    !
    ALLOCATE( sc( ibnd_start:ibnd_end, jbnd_start:jbnd_end ) )
    !
    ! ... <phi_i| S |psi_j>
    !
    IF ( uspp ) THEN
       !
       CALL ZGEMM( 'C', 'N', ibnd_size, jbnd_size, kdim, ONE, phi(1,ibnd_start), kdmx, &
                   spsi(1,jbnd_start), kdmx, ZERO, sc(ibnd_start,jbnd_start), ibnd_size )
       !
    ELSE
       !
       CALL ZGEMM( 'C', 'N', ibnd_size, jbnd_size, kdim, ONE, phi(1,ibnd_start), kdmx, &
                   psi(1,jbnd_start), kdmx, ZERO, sc(ibnd_start,jbnd_start), ibnd_size )
       !
    END IF
    !
    CALL mp_sum( sc, intra_bgrp_comm )
    !
    ! ... phi_j = phi_j - phi_i * <phi_i| S |psi_j>
    !
    CALL ZGEMM( 'N', 'N', kdim, jbnd_size, ibnd_size, MONE, phi(1,ibnd_start), kdmx, &
                sc(ibnd_start,jbnd_start), ibnd_size, ONE, phi(1,jbnd_start), kdmx )
    !
    IF ( eigen_ ) &
    CALL ZGEMM( 'N', 'N', kdim, jbnd_size, ibnd_size, MONE, hphi(1,ibnd_start), kdmx, &
                sc(ibnd_start,jbnd_start), ibnd_size, ONE, hphi(1,jbnd_start), kdmx )
    !
    IF ( uspp ) &
    CALL ZGEMM( 'N', 'N', kdim, jbnd_size, ibnd_size, MONE, sphi(1,ibnd_start), kdmx, &
                sc(ibnd_start,jbnd_start), ibnd_size, ONE, sphi(1,jbnd_start), kdmx )
    !
    DEALLOCATE( sc )
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
    REAL(DP), EXTERNAL :: DDOT
    !
    ! ... <psi_i| H |psi_i>
    !
    e(:) = 0._DP
    !
    CALL divide( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       e(ibnd) = DDOT( 2*kdim, psi(1,ibnd), 1, hpsi(1,ibnd), 1 ) 
       !
    END DO
    !
    CALL mp_sum( e(ibnd_start:ibnd_end), intra_bgrp_comm )
    CALL mp_sum( e, inter_bgrp_comm )
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
          CALL ZSWAP( kdim, psi(1,ibnd), 1, psi(1,ibnd-1), 1 )
          !
          IF ( eigen_ ) &
          CALL ZSWAP( kdim, hpsi(1,ibnd), 1, hpsi(1,ibnd-1), 1 )
          !
          IF ( uspp ) &
          CALL ZSWAP( kdim, spsi(1,ibnd), 1, spsi(1,ibnd-1), 1 )
          !
       END IF
       !
    END DO
    !
    IF ( nswap > 0 ) GOTO 10
    !
    RETURN
    !
  END SUBROUTINE sort_vectors
  !
  !
END SUBROUTINE gram_schmidt_k
