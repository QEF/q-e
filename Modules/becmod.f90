!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
MODULE becmod
  !
  ! ... *bec* contain <beta|psi> - used in h_psi, s_psi, many other places
  ! ... calbec( npw, beta, psi, betapsi [, nbnd ] ) is an interface calculating
  ! ...    betapsi(i,j)  = <beta(i)|psi(j)>   (the sum is over npw components)
  ! ... or betapsi(i,s,j)= <beta(i)|psi(s,j)> (s=polarization index)
  !
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : gamma_only, smallmem
  USE gvect,            ONLY : gstart
  USE noncollin_module, ONLY : noncolin, npol
  !
  SAVE
  !
  TYPE bec_type
     REAL(DP),   ALLOCATABLE :: r(:,:)    ! appropriate for gammaonly
     COMPLEX(DP),ALLOCATABLE :: k(:,:)    ! appropriate for generic k
     COMPLEX(DP),ALLOCATABLE :: nc(:,:,:)   ! appropriate for noncolin
     INTEGER :: comm
     INTEGER :: nbnd
     INTEGER :: nproc
     INTEGER :: mype
     INTEGER :: nbnd_loc
     INTEGER :: ibnd_begin
  END TYPE bec_type
  !
  TYPE (bec_type) :: becp  ! <beta|psi>

  PRIVATE
  !
  INTERFACE calbec
     !
     MODULE PROCEDURE calbec_k, calbec_gamma, calbec_gamma_nocomm, calbec_nc, calbec_bec_type
     !
  END INTERFACE

  INTERFACE becscal
     !
     MODULE PROCEDURE becscal_nck, becscal_gamma
     !
  END INTERFACE
  !
  PUBLIC :: bec_type, becp, allocate_bec_type, deallocate_bec_type, calbec, &
            beccopy, becscal, is_allocated_bec_type
  !
CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_bec_type ( npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !_
    USE mp_bands, ONLY: intra_bgrp_comm
    USE mp,       ONLY: mp_get_comm_null
    !
    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    TYPE (bec_type), INTENT (inout) :: betapsi ! NB: must be INOUT otherwise
                                               !  the allocatd array is lost
    INTEGER, INTENT (in) :: npw
    INTEGER, OPTIONAL :: nbnd
    !
    INTEGER :: local_nbnd
    INTEGER, EXTERNAL :: ldim_block, gind_block
    INTEGER :: m_loc, m_begin, ip
    REAL(DP), ALLOCATABLE :: dtmp(:,:)
    !
    IF ( present (nbnd) ) THEN
        local_nbnd = nbnd
    ELSE
        local_nbnd = size ( psi, 2)
    ENDIF

    IF ( gamma_only ) THEN
       !
       IF( betapsi%comm == mp_get_comm_null() ) THEN
          !
          CALL calbec_gamma ( npw, beta, psi, betapsi%r, local_nbnd, intra_bgrp_comm )
          !
       ELSE
          !
          ALLOCATE( dtmp( SIZE( betapsi%r, 1 ), SIZE( betapsi%r, 2 ) ) )
          !
          DO ip = 0, betapsi%nproc - 1
             m_loc   = ldim_block( betapsi%nbnd , betapsi%nproc, ip )
             m_begin = gind_block( 1,  betapsi%nbnd, betapsi%nproc, ip )
             IF( ( m_begin + m_loc - 1 ) > local_nbnd ) m_loc = local_nbnd - m_begin + 1
             IF( m_loc > 0 ) THEN
                CALL calbec_gamma ( npw, beta, psi(:,m_begin:m_begin+m_loc-1), dtmp, m_loc, betapsi%comm )
                IF( ip == betapsi%mype ) THEN
                   betapsi%r(:,1:m_loc) = dtmp(:,1:m_loc)
                END IF
             END IF
          END DO

          DEALLOCATE( dtmp )
          !
       END IF
       !
    ELSEIF ( noncolin) THEN
       !
       CALL  calbec_nc ( npw, beta, psi, betapsi%nc, local_nbnd )
       !
    ELSE
       !
       CALL  calbec_k ( npw, beta, psi, betapsi%k, local_nbnd )
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE calbec_bec_type
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_gamma_nocomm ( npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    USE mp_bands, ONLY: intra_bgrp_comm
    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    REAL (DP), INTENT (out) :: betapsi(:,:)
    INTEGER, INTENT (in) :: npw
    INTEGER, OPTIONAL :: nbnd
    INTEGER :: m
    IF ( present (nbnd) ) THEN
        m = nbnd
    ELSE
        m = size ( psi, 2)
    ENDIF
    CALL calbec_gamma ( npw, beta, psi, betapsi, m, intra_bgrp_comm )
    RETURN
    !
  END SUBROUTINE calbec_gamma_nocomm
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_gamma ( npw, beta, psi, betapsi, nbnd, comm )
    !-----------------------------------------------------------------------
    !
    ! ... matrix times matrix with summation index (k=1,npw) running on
    ! ... half of the G-vectors or PWs - assuming k=0 is the G=0 component:
    ! ... betapsi(i,j) = 2Re(\sum_k beta^*(i,k)psi(k,j)) + beta^*(i,0)psi(0,j)
    !
    USE mp,        ONLY : mp_sum

    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    REAL (DP), INTENT (out) :: betapsi(:,:)
    INTEGER, INTENT (in) :: npw
    INTEGER, INTENT (in) :: nbnd
    INTEGER, INTENT (in) :: comm 
    !
    INTEGER :: nkb, npwx, m
    !
    m = nbnd
    !
    nkb = size (beta, 2)
    IF ( nkb == 0 ) RETURN
    !
    CALL start_clock( 'calbec' )
    IF ( npw == 0 ) betapsi(:,:)=0.0_DP
    npwx= size (beta, 1)
    IF ( npwx /= size (psi, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
#ifdef DEBUG
    WRITE (*,*) 'calbec gamma'
    WRITE (*,*)  nkb,  size (betapsi,1) , m , size (betapsi, 2)
#endif
    IF ( nkb /= size (betapsi,1) .or. m > size (betapsi, 2) ) &
      CALL errore ('calbec', 'size mismatch', 3)
    !
    IF ( m == 1 ) THEN
        !
        CALL DGEMV( 'C', 2*npw, nkb, 2.0_DP, beta, 2*npwx, psi, 1, 0.0_DP, &
                     betapsi, 1 )
        IF ( gstart == 2 ) betapsi(:,1) = betapsi(:,1) - beta(1,:)*psi(1,1)
        !
    ELSE
        !
        CALL DGEMM( 'C', 'N', nkb, m, 2*npw, 2.0_DP, beta, 2*npwx, psi, &
                    2*npwx, 0.0_DP, betapsi, nkb )
        IF ( gstart == 2 ) &
           CALL DGER( nkb, m, -1.0_DP, beta, 2*npwx, psi, 2*npwx, betapsi, nkb )
        !
    ENDIF
    !
    CALL mp_sum( betapsi( :, 1:m ), comm )
    !
    CALL stop_clock( 'calbec' )
    !
    RETURN
    !
  END SUBROUTINE calbec_gamma
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_k ( npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !
    ! ... matrix times matrix with summation index (k=1,npw) running on
    ! ... G-vectors or PWs : betapsi(i,j) = \sum_k beta^*(i,k) psi(k,j)
    !
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp,       ONLY : mp_sum

    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    COMPLEX (DP), INTENT (out) :: betapsi(:,:)
    INTEGER, INTENT (in) :: npw
    INTEGER, OPTIONAL :: nbnd
    !
    INTEGER :: nkb, npwx, m
    !
    nkb = size (beta, 2)
    IF ( nkb == 0 ) RETURN
    !
    CALL start_clock( 'calbec' )
    IF ( npw == 0 ) betapsi(:,:)=(0.0_DP,0.0_DP)
    npwx= size (beta, 1)
    IF ( npwx /= size (psi, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
    IF ( present (nbnd) ) THEN
        m = nbnd
    ELSE
        m = size ( psi, 2)
    ENDIF
#ifdef DEBUG
    WRITE (*,*) 'calbec k'
    WRITE (*,*)  nkb,  size (betapsi,1) , m , size (betapsi, 2)
#endif
    IF ( nkb /= size (betapsi,1) .or. m > size (betapsi, 2) ) &
      CALL errore ('calbec', 'size mismatch', 3)
    !
    IF ( m == 1 ) THEN
       !
       CALL ZGEMV( 'C', npw, nkb, (1.0_DP,0.0_DP), beta, npwx, psi, 1, &
                   (0.0_DP, 0.0_DP), betapsi, 1 )
       !
    ELSE
       !
       CALL ZGEMM( 'C', 'N', nkb, m, npw, (1.0_DP,0.0_DP), &
                 beta, npwx, psi, npwx, (0.0_DP,0.0_DP), betapsi, nkb )
       !
    ENDIF
    !
    CALL mp_sum( betapsi( :, 1:m ), intra_bgrp_comm )
    !
    CALL stop_clock( 'calbec' )
    !
    RETURN
    !
  END SUBROUTINE calbec_k
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_nc ( npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !
    ! ... matrix times matrix with summation index (k below) running on
    ! ... G-vectors or PWs corresponding to two different polarizations:
    ! ... betapsi(i,1,j) = \sum_k=1,npw beta^*(i,k) psi(k,j)
    ! ... betapsi(i,2,j) = \sum_k=1,npw beta^*(i,k) psi(k+npwx,j)
    !
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp,       ONLY : mp_sum

    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    COMPLEX (DP), INTENT (out) :: betapsi(:,:,:)
    INTEGER, INTENT (in) :: npw
    INTEGER, OPTIONAL :: nbnd
    !
    INTEGER :: nkb, npwx, npol, m
    !
    nkb = size (beta, 2)
    IF ( nkb == 0 ) RETURN
    !
    CALL start_clock ('calbec')
    IF ( npw == 0 ) betapsi(:,:,:)=(0.0_DP,0.0_DP)
    npwx= size (beta, 1)
    IF ( 2*npwx /= size (psi, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
    IF ( present (nbnd) ) THEN
        m = nbnd
    ELSE
        m = size ( psi, 2)
    ENDIF
    npol= size (betapsi, 2)
#ifdef DEBUG
    WRITE (*,*) 'calbec nc'
    WRITE (*,*)  nkb,  size (betapsi,1) , m , size (betapsi, 3)
#endif
    IF ( nkb /= size (betapsi,1) .or. m > size (betapsi, 3) ) &
      CALL errore ('calbec', 'size mismatch', 3)
    !
    CALL ZGEMM ('C', 'N', nkb, m*npol, npw, (1.0_DP, 0.0_DP), beta, &
              npwx, psi, npwx, (0.0_DP, 0.0_DP),  betapsi, nkb)
    !
    CALL mp_sum( betapsi( :, :, 1:m ), intra_bgrp_comm )
    !
    CALL stop_clock( 'calbec' )
    !
    RETURN
    !
  END SUBROUTINE calbec_nc
  !
  !
  !-----------------------------------------------------------------------
  FUNCTION is_allocated_bec_type (bec) RESULT (isalloc)
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    TYPE (bec_type) :: bec
    LOGICAL :: isalloc
    isalloc = (allocated(bec%r) .or. allocated(bec%nc) .or. allocated(bec%k))
    RETURN
    !
    !-----------------------------------------------------------------------
  END FUNCTION is_allocated_bec_type
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE allocate_bec_type ( nkb, nbnd, bec, comm )
    !-----------------------------------------------------------------------
    USE mp, ONLY: mp_size, mp_rank, mp_get_comm_null
    IMPLICIT NONE
    TYPE (bec_type) :: bec
    INTEGER, INTENT (in) :: nkb, nbnd
    INTEGER, INTENT (in), OPTIONAL :: comm
    INTEGER :: ierr, nbnd_siz
    INTEGER, EXTERNAL :: ldim_block, gind_block
    !
    nbnd_siz = nbnd
    bec%comm = mp_get_comm_null()
    bec%nbnd = nbnd
    bec%mype = 0
    bec%nproc = 1
    bec%nbnd_loc = nbnd
    bec%ibnd_begin = 1
    !
    IF( PRESENT( comm ) .AND. gamma_only .AND. smallmem ) THEN
       bec%comm = comm
       bec%nproc = mp_size( comm )
       IF( bec%nproc > 1 ) THEN
          nbnd_siz   = nbnd / bec%nproc
          IF( MOD( nbnd, bec%nproc ) /= 0 ) nbnd_siz = nbnd_siz + 1
          bec%mype  = mp_rank( bec%comm )
          bec%nbnd_loc   = ldim_block( becp%nbnd , bec%nproc, bec%mype )
          bec%ibnd_begin = gind_block( 1,  becp%nbnd, bec%nproc, bec%mype )
       END IF
    END IF
    !
    IF ( gamma_only ) THEN
       !
       ALLOCATE( bec%r( nkb, nbnd_siz ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' allocate_bec_type ', ' cannot allocate bec%r ', ABS(ierr) )
       !
       bec%r(:,:)=0.0D0
       !
    ELSEIF ( noncolin) THEN
       !
       ALLOCATE( bec%nc( nkb, npol, nbnd_siz ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' allocate_bec_type ', ' cannot allocate bec%nc ', ABS(ierr) )
       !
       bec%nc(:,:,:)=(0.0D0,0.0D0)
       !
    ELSE
       !
       ALLOCATE( bec%k( nkb, nbnd_siz ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' allocate_bec_type ', ' cannot allocate bec%k ', ABS(ierr) )
       !
       bec%k(:,:)=(0.0D0,0.0D0)
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE allocate_bec_type
  !
  !-----------------------------------------------------------------------
  SUBROUTINE deallocate_bec_type (bec)
    !-----------------------------------------------------------------------
    !
    USE mp, ONLY: mp_get_comm_null
    IMPLICIT NONE
    TYPE (bec_type) :: bec
    !
    bec%comm = mp_get_comm_null()
    bec%nbnd = 0
    !
    IF (allocated(bec%r))  DEALLOCATE(bec%r)
    IF (allocated(bec%nc)) DEALLOCATE(bec%nc)
    IF (allocated(bec%k))  DEALLOCATE(bec%k)
    !
    RETURN
    !
  END SUBROUTINE deallocate_bec_type

  SUBROUTINE beccopy(bec, bec1, nkb, nbnd)
    IMPLICIT NONE
    TYPE(bec_type), INTENT(in) :: bec
    TYPE(bec_type)  :: bec1
    INTEGER, INTENT(in) :: nkb, nbnd

    IF (gamma_only) THEN
       CALL dcopy(nkb*nbnd, bec%r, 1, bec1%r, 1)
    ELSEIF (noncolin) THEN
       CALL zcopy(nkb*npol*nbnd, bec%nc, 1, bec1%nc,  1)
    ELSE
       CALL zcopy(nkb*nbnd, bec%k, 1, bec1%k, 1)
    ENDIF

    RETURN
  END SUBROUTINE beccopy

  SUBROUTINE becscal_nck(alpha, bec, nkb, nbnd)
    IMPLICIT NONE
    TYPE(bec_type), INTENT(INOUT) :: bec
    COMPLEX(DP), INTENT(IN) :: alpha
    INTEGER, INTENT(IN) :: nkb, nbnd

    IF (gamma_only) THEN
       CALL errore('becscal_nck','called in the wrong case',1)
    ELSEIF (noncolin) THEN
       CALL zscal(nkb*npol*nbnd, alpha, bec%nc, 1)
    ELSE
       CALL zscal(nkb*nbnd, alpha, bec%k, 1)
    ENDIF

    RETURN
  END SUBROUTINE becscal_nck

  SUBROUTINE becscal_gamma(alpha, bec, nkb, nbnd)
    IMPLICIT NONE
    TYPE(bec_type), INTENT(INOUT) :: bec
    REAL(DP), INTENT(IN) :: alpha
    INTEGER, INTENT(IN) :: nkb, nbnd

    IF (gamma_only) THEN
       CALL dscal(nkb*nbnd, alpha, bec%r, 1)
    ELSE
       CALL errore('becscal_gamma','called in the wrong case',1)
    ENDIF

    RETURN
  END SUBROUTINE becscal_gamma

END MODULE becmod
