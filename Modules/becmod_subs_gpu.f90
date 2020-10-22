!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
MODULE becmod_subs_gpum
  
  ! NOTA BENE : THE SUBROUTINES IN THIS FILE ARE ONLY PARTIALLY TESTED!
   
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
  USE becmod_gpum,      ONLY : bec_type_d
  !
  SAVE
  !
  PRIVATE
  !
  INTERFACE calbec_gpu
     !
     MODULE PROCEDURE calbec_k_gpu, calbec_gamma_gpu, calbec_gamma_nocomm_gpu, calbec_nc_gpu, calbec_bec_type_gpu
     !
  END INTERFACE

  INTERFACE becscal_gpu
     !
     MODULE PROCEDURE becscal_nck_gpu, becscal_gamma_gpu
     !
  END INTERFACE
  !
  PUBLIC :: allocate_bec_type_gpu, deallocate_bec_type_gpu, calbec_gpu, &
            beccopy_gpu, becscal_gpu, is_allocated_bec_type_gpu, &
            synchronize_bec_type_gpu, &
            using_becp_auto, using_becp_d_auto
  !
CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_bec_type_gpu ( npw, beta_d, psi_d, betapsi_d, nbnd )
    !-----------------------------------------------------------------------
    !_
    USE mp_bands, ONLY: intra_bgrp_comm
    USE mp,       ONLY: mp_get_comm_null
    !
    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta_d(:,:), psi_d(:,:)
    TYPE (bec_type_d), TARGET, INTENT (inout) :: betapsi_d ! NB: must be INOUT otherwise
                                               !  the allocatd array is lost
    INTEGER, INTENT (in) :: npw
    INTEGER, OPTIONAL :: nbnd
#if defined(__CUDA)
    attributes(DEVICE) :: beta_d, psi_d
#endif
    !
    INTEGER :: local_nbnd
    INTEGER, EXTERNAL :: ldim_block, gind_block
    INTEGER :: m_loc, m_begin, ip
    REAL(DP), ALLOCATABLE :: dtmp_d(:,:)   ! replace this with buffers !
    INTEGER :: i, j, nkb
    REAL(DP), POINTER :: betapsi_d_r_d(:,:)
#if defined(__CUDA)
    attributes(DEVICE) :: dtmp_d, betapsi_d_r_d
#endif
    !
    IF ( present (nbnd) ) THEN
        local_nbnd = nbnd
    ELSE
        local_nbnd = size ( psi_d, 2)
    ENDIF

    IF ( gamma_only ) THEN
       !
       IF( betapsi_d%comm == mp_get_comm_null() ) THEN
          !
          CALL calbec_gamma_gpu ( npw, beta_d, psi_d, betapsi_d%r_d, local_nbnd, intra_bgrp_comm )
          !
       ELSE
          !
          ALLOCATE( dtmp_d( SIZE( betapsi_d%r_d, 1 ), SIZE( betapsi_d%r_d, 2 ) ) )
          !
          DO ip = 0, betapsi_d%nproc - 1
             m_loc   = ldim_block( betapsi_d%nbnd , betapsi_d%nproc, ip )
             m_begin = gind_block( 1,  betapsi_d%nbnd, betapsi_d%nproc, ip )
             IF( ( m_begin + m_loc - 1 ) > local_nbnd ) m_loc = local_nbnd - m_begin + 1
             IF( m_loc > 0 ) THEN
                CALL calbec_gamma_gpu ( npw, beta_d, psi_d(:,m_begin:m_begin+m_loc-1), dtmp_d, m_loc, betapsi_d%comm )
                IF( ip == betapsi_d%mype ) THEN
                   nkb = SIZE( betapsi_d%r_d, 1 )
                   betapsi_d_r_d => betapsi_d%r_d
!$cuf kernel do(2) <<<*,*>>>
                   DO j=1,m_loc
                      DO i=1, nkb
                         betapsi_d_r_d(i,j) = dtmp_d(i,j)
                      END DO
                   END DO
                END IF
             END IF
          END DO

          DEALLOCATE( dtmp_d )
          !
       END IF
       !
    ELSEIF ( noncolin) THEN
       !
       CALL  calbec_nc_gpu ( npw, beta_d, psi_d, betapsi_d%nc_d, local_nbnd )
       !
    ELSE
       !
       CALL  calbec_k_gpu ( npw, beta_d, psi_d, betapsi_d%k_d, local_nbnd )
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE calbec_bec_type_gpu
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_gamma_nocomm_gpu ( npw, beta_d, psi_d, betapsi_d, nbnd )
    !-----------------------------------------------------------------------
    USE mp_bands, ONLY: intra_bgrp_comm
    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta_d(:,:), psi_d(:,:)
    REAL (DP), INTENT (out) :: betapsi_d(:,:)
    INTEGER, INTENT (in) :: npw
    INTEGER, OPTIONAL :: nbnd
    INTEGER :: m
#if defined(__CUDA)
    attributes(DEVICE) :: beta_d, psi_d, betapsi_d
#endif
    IF ( present (nbnd) ) THEN
        m = nbnd
    ELSE
        m = size ( psi_d, 2)
    ENDIF
    CALL calbec_gamma_gpu ( npw, beta_d, psi_d, betapsi_d, m, intra_bgrp_comm )
    RETURN
    !
  END SUBROUTINE calbec_gamma_nocomm_gpu
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_gamma_gpu ( npw, beta_d, psi_d, betapsi_d, nbnd, comm )
    !-----------------------------------------------------------------------
    !
    ! ... matrix times matrix with summation index (k=1,npw) running on
    ! ... half of the G-vectors or PWs - assuming k=0 is the G=0 component:
    ! ... betapsi(i,j) = 2Re(\sum_k beta^*(i,k)psi(k,j)) + beta^*(i,0)psi(0,j)
    !
    USE mp,        ONLY : mp_sum, mp_size
#if defined(__CUDA)
    USE cudafor
    USE cublas
#endif
    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta_d(:,:), psi_d(:,:)
    REAL (DP), INTENT (out) :: betapsi_d(:,:)
    INTEGER, INTENT (in) :: npw
    INTEGER, INTENT (in) :: nbnd
    INTEGER, INTENT (in) :: comm 
    !
#if defined(__CUDA)
    attributes(DEVICE) :: beta_d, psi_d, betapsi_d
#endif
    INTEGER :: nkb, npwx, m
    INTEGER :: i,j
    !
    m = nbnd
    !
    nkb = size (beta_d, 2)
    IF ( nkb == 0 ) RETURN
    !
    CALL start_clock( 'calbec' )
    IF ( npw == 0 ) betapsi_d(:,:)=0.0_DP
    npwx= size (beta_d, 1)
    IF ( npwx /= size (psi_d, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
#if defined(DEBUG)
    WRITE (*,*) 'calbec gamma'
    WRITE (*,*)  nkb,  size (betapsi_d,1) , m , size (betapsi_d, 2)
#endif
    IF ( nkb /= size (betapsi_d,1) .or. m > size (betapsi_d, 2) ) &
      CALL errore ('calbec', 'size mismatch', 3)
    !
    IF ( m == 1 ) THEN
        !
        CALL cudaDGEMV( 'C', 2*npw, nkb, 2.0_DP, beta_d, 2*npwx, psi_d, 1, 0.0_DP, &
                     betapsi_d, 1 )
        IF ( gstart == 2 ) THEN
           !betapsi_d(:,1) = betapsi_d(:,1) - beta_d(1,:)*psi_d(1,1)
           !$cuf kernel do(1) <<<*,*>>>
           DO i=1, nkb
              betapsi_d(i,1) = betapsi_d(i,1) - DBLE(beta_d(1,i)*psi_d(1,1))
           END DO
        END IF
        !
    ELSE
        !
        CALL DGEMM( 'C', 'N', nkb, m, 2*npw, 2.0_DP, beta_d, 2*npwx, psi_d, &
                    2*npwx, 0.0_DP, betapsi_d, nkb )
        IF ( gstart == 2 ) &
           CALL cudaDGER( nkb, m, -1.0_DP, beta_d, 2*npwx, psi_d, 2*npwx, betapsi_d, nkb )
        !
    ENDIF
    !
    IF (mp_size(comm) > 1) CALL mp_sum( betapsi_d( :, 1:m ), comm )
    !
    CALL stop_clock( 'calbec' )
    !
    RETURN
    !
  END SUBROUTINE calbec_gamma_gpu
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_k_gpu ( npw, beta_d, psi_d, betapsi_d, nbnd )
    !-----------------------------------------------------------------------
    !
    ! ... matrix times matrix with summation index (k=1,npw) running on
    ! ... G-vectors or PWs : betapsi(i,j) = \sum_k beta^*(i,k) psi(k,j)
    !
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp,       ONLY : mp_sum, mp_size
#if defined(__CUDA)
    USE cudafor
    USE cublas
#endif
    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta_d(:,:), psi_d(:,:)
    COMPLEX (DP), INTENT (out) :: betapsi_d(:,:)
    INTEGER, INTENT (in) :: npw
    INTEGER, OPTIONAL :: nbnd
    !
    INTEGER :: nkb, npwx, m
    !
#if defined(__CUDA)
    attributes(device) :: beta_d, psi_d, betapsi_d
#endif
    nkb = size (beta_d, 2)
    IF ( nkb == 0 ) RETURN
    !
    CALL start_clock( 'calbec' )
    IF ( npw == 0 ) betapsi_d(:,:)=(0.0_DP,0.0_DP)
    npwx= size (beta_d, 1)
    IF ( npwx /= size (psi_d, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
    IF ( present (nbnd) ) THEN
        m = nbnd
    ELSE
        m = size ( psi_d, 2)
    ENDIF
#if defined(DEBUG)
    WRITE (*,*) 'calbec k'
    WRITE (*,*)  nkb,  size (betapsi_d,1) , m , size (betapsi_d, 2)
#endif
    IF ( nkb /= size (betapsi_d,1) .or. m > size (betapsi_d, 2) ) &
      CALL errore ('calbec', 'size mismatch', 3)
    !
    IF ( m == 1 ) THEN
       !
       CALL ZGEMV( 'C', npw, nkb, (1.0_DP,0.0_DP), beta_d, npwx, psi_d, 1, &
                   (0.0_DP, 0.0_DP), betapsi_d, 1 )
       !
    ELSE
       !
       CALL ZGEMM( 'C', 'N', nkb, m, npw, (1.0_DP,0.0_DP), &
                 beta_d, npwx, psi_d, npwx, (0.0_DP,0.0_DP), betapsi_d, nkb )
       !
    ENDIF
    !
    IF (mp_size(intra_bgrp_comm) > 1) CALL mp_sum( betapsi_d( :, 1:m ), intra_bgrp_comm )
    !
    CALL stop_clock( 'calbec' )
    !
    RETURN
    !
  END SUBROUTINE calbec_k_gpu
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_nc_gpu ( npw, beta_d, psi_d, betapsi_d, nbnd )
    !-----------------------------------------------------------------------
    !
    ! ... matrix times matrix with summation index (k below) running on
    ! ... G-vectors or PWs corresponding to two different polarizations:
    ! ... betapsi(i,1,j) = \sum_k=1,npw beta^*(i,k) psi(k,j)
    ! ... betapsi(i,2,j) = \sum_k=1,npw beta^*(i,k) psi(k+npwx,j)
    !
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp,       ONLY : mp_sum, mp_size
#if defined(__CUDA)
    USE cudafor
    USE cublas
#endif
    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta_d(:,:), psi_d(:,:)
    COMPLEX (DP), INTENT (out) :: betapsi_d(:,:,:)
    INTEGER, INTENT (in) :: npw
    INTEGER, OPTIONAL :: nbnd
    !
    INTEGER :: nkb, npwx, npol, m
    !
#if defined(__CUDA)
    attributes(device) :: beta_d, psi_d, betapsi_d
#endif
    nkb = size (beta_d, 2)
    IF ( nkb == 0 ) RETURN
    !
    CALL start_clock ('calbec')
    IF ( npw == 0 ) betapsi_d(:,:,:)=(0.0_DP,0.0_DP)
    npwx= size (beta_d, 1)
    IF ( 2*npwx /= size (psi_d, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
    IF ( present (nbnd) ) THEN
        m = nbnd
    ELSE
        m = size ( psi_d, 2)
    ENDIF
    npol= size (betapsi_d, 2)
#if defined(DEBUG)
    WRITE (*,*) 'calbec nc'
    WRITE (*,*)  nkb,  size (betapsi_d,1) , m , size (betapsi_d, 3)
#endif
    IF ( nkb /= size (betapsi_d,1) .or. m > size (betapsi_d, 3) ) &
      CALL errore ('calbec', 'size mismatch', 3)
    !
    CALL ZGEMM ('C', 'N', nkb, m*npol, npw, (1.0_DP, 0.0_DP), beta_d, &
              npwx, psi_d, npwx, (0.0_DP, 0.0_DP),  betapsi_d, nkb)
    !
    IF (mp_size(intra_bgrp_comm) > 1)  CALL mp_sum( betapsi_d( :, :, 1:m ), intra_bgrp_comm )
    !
    CALL stop_clock( 'calbec' )
    !
    RETURN
    !
  END SUBROUTINE calbec_nc_gpu
  !
  !
  !-----------------------------------------------------------------------
  FUNCTION is_allocated_bec_type_gpu (bec_d) RESULT (isalloc)
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    TYPE (bec_type_d) :: bec_d
    LOGICAL :: isalloc
    isalloc = (allocated(bec_d%r_d) .or. allocated(bec_d%nc_d) .or. allocated(bec_d%k_d))
    RETURN
    !
    !-----------------------------------------------------------------------
  END FUNCTION is_allocated_bec_type_gpu
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE allocate_bec_type_gpu ( nkb, nbnd, bec_d, comm )
    !-----------------------------------------------------------------------
    USE mp, ONLY: mp_size, mp_rank, mp_get_comm_null
    USE device_memcpy_m, ONLY : dev_memset
    IMPLICIT NONE
    TYPE (bec_type_d) :: bec_d
    INTEGER, INTENT (in) :: nkb, nbnd
    INTEGER, INTENT (in), OPTIONAL :: comm
    INTEGER :: ierr, nbnd_siz
    INTEGER, EXTERNAL :: ldim_block, gind_block
    !
    nbnd_siz = nbnd
    bec_d%comm = mp_get_comm_null()
    bec_d%nbnd = nbnd
    bec_d%mype = 0
    bec_d%nproc = 1
    bec_d%nbnd_loc = nbnd
    bec_d%ibnd_begin = 1
    !
    IF( PRESENT( comm ) .AND. gamma_only .AND. smallmem ) THEN
       bec_d%comm = comm
       bec_d%nproc = mp_size( comm )
       IF( bec_d%nproc > 1 ) THEN
          nbnd_siz   = nbnd / bec_d%nproc
          IF( MOD( nbnd, bec_d%nproc ) /= 0 ) nbnd_siz = nbnd_siz + 1
          bec_d%mype  = mp_rank( bec_d%comm )
          bec_d%nbnd_loc   = ldim_block( bec_d%nbnd , bec_d%nproc, bec_d%mype )
          bec_d%ibnd_begin = gind_block( 1,  bec_d%nbnd, bec_d%nproc, bec_d%mype )
       END IF
    END IF
    !
    IF ( gamma_only ) THEN
       !
       ALLOCATE( bec_d%r_d( nkb, nbnd_siz ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' allocate_bec_type ', ' cannot allocate bec_d%r ', ABS(ierr) )
       !
       CALL dev_memset(bec_d%r_d, 0.0D0, (/1,nkb/), 1, (/1, nbnd_siz/), 1)
       !
    ELSEIF ( noncolin) THEN
       !
       ALLOCATE( bec_d%nc_d( nkb, npol, nbnd_siz ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' allocate_bec_type ', ' cannot allocate bec_d%nc ', ABS(ierr) )
       !
       CALL dev_memset(bec_d%nc_d, (0.0D0,0.0D0), (/1, nkb/), 1, (/1, npol/), 1, (/1, nbnd_siz/), 1)
       !
    ELSE
       !
       ALLOCATE( bec_d%k_d( nkb, nbnd_siz ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' allocate_bec_type ', ' cannot allocate bec_d%k ', ABS(ierr) )
       !
       CALL dev_memset(bec_d%k_d, (0.0D0,0.0D0), (/1, nkb/), 1, (/1, npol/), 1)
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE allocate_bec_type_gpu
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE deallocate_bec_type_gpu (bec_d)
    !-----------------------------------------------------------------------
    !
    USE mp, ONLY: mp_get_comm_null
    IMPLICIT NONE
    TYPE (bec_type_d) :: bec_d
    !
    bec_d%comm = mp_get_comm_null()
    bec_d%nbnd = 0
    !
    IF (allocated(bec_d%r_d))  DEALLOCATE(bec_d%r_d)
    IF (allocated(bec_d%nc_d)) DEALLOCATE(bec_d%nc_d)
    IF (allocated(bec_d%k_d))  DEALLOCATE(bec_d%k_d)
    !
    RETURN
    !
  END SUBROUTINE deallocate_bec_type_gpu
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE 	synchronize_bec_type_gpu (bec_d, bec, what)
    !-----------------------------------------------------------------------
    !
    ! ... Updates a device or host version of a bec_type variable.
    ! ... Direction 'h' updates host version, direction 'd' the device
    ! ... version.
    !
    USE mp, ONLY: mp_get_comm_null
    USE becmod, ONLY : bec_type
    IMPLICIT NONE
    TYPE (bec_type_d) :: bec_d
    TYPE (bec_type) :: bec
    CHARACTER, INTENT(IN) :: what
    !
    IF ( gamma_only ) THEN
       !
       IF (.not. (allocated(bec_d%r_d) .and. allocated(bec%r))) &
          CALL errore('becmod_gpu', 'Unallocated array',1)
       SELECT CASE(what)
        CASE('d')
          bec_d%r_d = bec%r
        CASE('h')
          bec%r     = bec_d%r_d
        CASE DEFAULT
          CALL errore('becmod_gpu', 'Invalid command',2)
       END SELECT
       !
    ELSEIF ( noncolin) THEN
       !
       IF (.not. (allocated(bec_d%nc_d) .and. allocated(bec%nc))) &
          CALL errore('becmod_gpu', 'Unallocated array',3)
       SELECT CASE(what)
        CASE('d')
          bec_d%nc_d = bec%nc
        CASE('h')
          bec%nc     = bec_d%nc_d
        CASE DEFAULT
          CALL errore('becmod_gpu', 'Invalid command',4)
       END SELECT
       !
    ELSE
       !
       IF (.not. (allocated(bec_d%k_d) .and. allocated(bec%k))) &
          CALL errore('becmod_gpu', 'Unallocated array',5)
       SELECT CASE(what)
        CASE('d')
          bec_d%k_d = bec%k
        CASE('h')
          bec%k     = bec_d%k_d
        CASE DEFAULT
          CALL errore('becmod_gpu', 'Invalid command',6)
       END SELECT
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE 	synchronize_bec_type_gpu
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE beccopy_gpu(bec, bec1, nkb, nbnd)
#if defined(__CUDA)
    USE cudafor
    USE cublas
#endif
    IMPLICIT NONE
    TYPE(bec_type_d), INTENT(in) :: bec
    TYPE(bec_type_d)  :: bec1
    INTEGER, INTENT(in) :: nkb, nbnd

    IF (gamma_only) THEN
       CALL dcopy(nkb*nbnd, bec%r_d, 1, bec1%r_d, 1)
    ELSEIF (noncolin) THEN
       CALL zcopy(nkb*npol*nbnd, bec%nc_d, 1, bec1%nc_d,  1)
    ELSE
       CALL zcopy(nkb*nbnd, bec%k_d, 1, bec1%k_d, 1)
    ENDIF

    RETURN
  END SUBROUTINE beccopy_gpu

  SUBROUTINE becscal_nck_gpu(alpha, bec_d, nkb, nbnd)
#if defined(__CUDA)
    USE cudafor
    USE cublas
#endif
    IMPLICIT NONE
    TYPE(bec_type_d), INTENT(INOUT) :: bec_d
    COMPLEX(DP), INTENT(IN) :: alpha
    INTEGER, INTENT(IN) :: nkb, nbnd
    IF (gamma_only) THEN
       CALL errore('becscal_nck','called in the wrong case',1)
    ELSEIF (noncolin) THEN
       CALL zscal(nkb*npol*nbnd, alpha, bec_d%nc_d, 1)
    ELSE
       CALL zscal(nkb*nbnd, alpha, bec_d%k_d, 1)
    ENDIF

    RETURN
  END SUBROUTINE becscal_nck_gpu

  SUBROUTINE becscal_gamma_gpu(alpha, bec_d, nkb, nbnd)
#if defined(__CUDA)
    USE cudafor
    USE cublas
#endif
    IMPLICIT NONE
    TYPE(bec_type_d), INTENT(INOUT) :: bec_d
    REAL(DP), INTENT(IN) :: alpha
    INTEGER, INTENT(IN) :: nkb, nbnd

    IF (gamma_only) THEN
       CALL dscal(nkb*nbnd, alpha, bec_d%r_d, 1)
    ELSE
       CALL errore('becscal_gamma','called in the wrong case',1)
    ENDIF

    RETURN
  END SUBROUTINE becscal_gamma_gpu
  !
  SUBROUTINE using_becp_auto(intento)
    USE becmod_gpum, ONLY : using_becp_r
    USE becmod_gpum, ONLY : using_becp_k
    USE becmod_gpum, ONLY : using_becp_nc
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: intento
    !
    !
    IF ( gamma_only ) THEN
       !
       CALL using_becp_r(intento)
       !
    ELSEIF ( noncolin) THEN
       !
       CALL using_becp_nc(intento)
       !
    ELSE
       !
       CALL using_becp_k(intento)
       !
    ENDIF
  END SUBROUTINE using_becp_auto
  !
  SUBROUTINE using_becp_d_auto(intento)
    USE becmod_gpum, ONLY : using_becp_r_d
    USE becmod_gpum, ONLY : using_becp_k_d
    USE becmod_gpum, ONLY : using_becp_nc_d
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: intento
    !
    !
    IF ( gamma_only ) THEN
       !
       CALL using_becp_r_d(intento)
       !
    ELSEIF ( noncolin) THEN
       !
       CALL using_becp_nc_d(intento)
       !
    ELSE
       !
       CALL using_becp_k_d(intento)
       !
    ENDIF
  END SUBROUTINE using_becp_d_auto


END MODULE becmod_subs_gpum
