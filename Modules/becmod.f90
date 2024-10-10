!
! Copyright (C) 2001-2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE becmod
  !---------------------------------------------------------------------------
  !! This module contains \(\langle\text{beta}|\text{psi}\rangle\) - used
  !! in \(\texttt{h_psi}\), \(\texttt{s_psi}\) and many other places.  
  !! \(\texttt{calbec}\) is an interface calculating \(\text{betapsi}(i,j) =
  !! \langle \text{beta}(i)|\text{psi}(j)\rangle \) (the sum is over \(\text{npw}\) 
  !! components) or \( \text{betapsi}(i,s,j)= \langle\text{beta}(i)|\text{psi}(s,j)
  !! \rangle\) (s=polarization index).
  !
  USE kinds,            ONLY : DP, offload_kind_acc, offload_kind_cpu, offload_kind_omp
  USE control_flags,    ONLY : gamma_only, smallmem, offload_acc, offload_cpu, offload_omp
  USE gvect,            ONLY : gstart
  USE noncollin_module, ONLY : noncolin, npol
  !
  SAVE
  !
  TYPE bec_type
     REAL(DP),   ALLOCATABLE :: r(:,:)
     !! appropriate for gammaonly
     COMPLEX(DP),ALLOCATABLE :: k(:,:)
     !! appropriate for generic k
     COMPLEX(DP),ALLOCATABLE :: nc(:,:,:)
     !! appropriate for noncolin
     INTEGER :: nbnd
  END TYPE bec_type
  !
  TYPE (bec_type) :: becp
  !! \(\langle\text{beta}|\text{psi}\rangle\)
  !
  PRIVATE
  !
  INTERFACE calbec
     !
     MODULE PROCEDURE &
!                     usage: call calbec( offload_type, beta, psi, betapsi )
!                             - beta, psi, betapsi can be CPU, OpenACC (or OpenMP5), 
!                             - CPU, OpenACC (and OpenMP5) cases are distinguished by type(offload_type)
                      calbec_k_acc, calbec_gamma_acc, calbec_nc_acc, calbec_bec_type_acc, &
                      calbec_k_cpu, calbec_gamma_cpu, calbec_nc_cpu, calbec_bec_type_cpu, &
!                     usage: call calbec( beta, psi, betapsi ) ("old" way to call calbec on CPU)
!                             - beta, psi, betapsi are CPU-only 
!                             - this allows to keep unchanged ALL the "old" calbec calls in the code
                      calbec_k,     calbec_gamma,     calbec_nc,     calbec_bec_type
     !
  END INTERFACE
  !
  INTERFACE becscal
     !
     MODULE PROCEDURE becscal_nck, becscal_gamma
     !
  END INTERFACE
  !
  INTERFACE becupdate
     !
     MODULE PROCEDURE becupdate_0D_cpu, becupdate_1D_cpu, becupdate_2D_cpu, &
                      becupdate_0D_acc, becupdate_1D_acc, becupdate_2D_acc 
     !
  END INTERFACE
  !
  PUBLIC :: bec_type, becp, calbec, beccopy, becscal, becupdate, is_allocated_bec_type, &
            allocate_bec_type,     deallocate_bec_type, &
            allocate_bec_type_acc, deallocate_bec_type_acc
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_bec_type_acc ( offload, npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !
    ! beta, psi, betapsi, are assumed OpenACC data on GPU
    !
    IMPLICIT NONE
    TYPE(offload_kind_acc), INTENT(IN) :: offload
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    TYPE (bec_type), INTENT (inout) :: betapsi ! NB: must be INOUT otherwise
                                               !  the allocatd array is lost
    INTEGER, INTENT (in) :: npw
    INTEGER, OPTIONAL :: nbnd
    !
    INTEGER :: local_nbnd
    !
    IF ( present (nbnd) ) THEN
        local_nbnd = nbnd
    ELSE
        local_nbnd = size ( psi, 2)
    ENDIF

    IF ( gamma_only ) THEN
       !
       CALL calbec_gamma_acc ( offload_acc, npw, beta, psi, betapsi%r, local_nbnd )
       !
    ELSEIF ( noncolin) THEN
       !
       CALL  calbec_nc_acc ( offload_acc, npw, beta, psi, betapsi%nc, local_nbnd )
       !
    ELSE
       !
       CALL  calbec_k_acc ( offload_acc, npw, beta, psi, betapsi%k, local_nbnd )
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE calbec_bec_type_acc
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_bec_type_cpu ( offload, npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(offload_kind_cpu), INTENT(IN) :: offload
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    TYPE (bec_type), INTENT (inout) :: betapsi ! NB: must be INOUT otherwise
                                               !  the allocatd array is lost
    INTEGER, INTENT (in) :: npw
    INTEGER, OPTIONAL :: nbnd
    INTEGER :: m
    !
    IF ( present (nbnd) ) THEN
        m = nbnd
    ELSE
        m = size ( psi, 2)
    ENDIF
    Call calbec_bec_type ( npw, beta, psi, betapsi, m )
    !
    RETURN
    !
  END SUBROUTINE calbec_bec_type_cpu
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_bec_type ( npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    TYPE (bec_type), INTENT (inout) :: betapsi ! NB: must be INOUT otherwise
                                               !  the allocatd array is lost
    INTEGER, INTENT (in) :: npw
    INTEGER, OPTIONAL :: nbnd
    !
    INTEGER :: local_nbnd
    !
    IF ( present (nbnd) ) THEN
        local_nbnd = nbnd
    ELSE
        local_nbnd = size ( psi, 2)
    ENDIF

    IF ( gamma_only ) THEN
       !
       CALL calbec_gamma ( npw, beta, psi, betapsi%r, local_nbnd )
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
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_gamma_acc ( offload, npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !! matrix times matrix with summation index (k=1,npw) running on
    !! half of the G-vectors or PWs - assuming k=0 is the G=0 component:
    !
    !! $$ betapsi(i,j) = 2Re(\sum_k beta^*(i,k)psi(k,j)) + beta^*(i,0)psi(0,j) $$
    !
    ! beta, psi, betapsi, are assumed OpenACC data on GPU
    !
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp,        ONLY : mp_sum, mp_size
    !
    IMPLICIT NONE
    TYPE(offload_kind_acc), INTENT(IN) :: offload
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    REAL (DP), INTENT (out) :: betapsi(:,:)
    INTEGER, INTENT (in) :: npw
    INTEGER, INTENT (in), OPTIONAL :: nbnd
    !
    INTEGER :: nkb, npwx, m
    !
    IF ( present (nbnd) ) THEN
        m = nbnd
    ELSE
        m = size ( psi, 2)
    ENDIF
    !
    nkb = size (beta, 2)
    IF ( nkb == 0 ) RETURN
    !
    CALL start_clock( 'calbec' )
    IF ( npw == 0 ) THEN
      !$acc kernels
      betapsi(:,:)=0.0_DP
      !$acc end kernels
    END IF
    npwx= size (beta, 1)
    IF ( npwx /= size (psi, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
#if defined(DEBUG)
    WRITE (*,*) 'calbec gamma'
    WRITE (*,*)  nkb,  size (betapsi,1) , m , size (betapsi, 2)
#endif
    IF ( nkb /= size (betapsi,1) .or. m > size (betapsi, 2) ) &
      CALL errore ('calbec', 'size mismatch', 3)
    !
    IF ( m == 1 ) THEN
        !
        !$acc host_data use_device(beta,psi,betapsi)
        CALL MYDGEMV( 'C', 2*npw, nkb, 2.0_DP, beta, 2*npwx, psi, 1, 0.0_DP, &
                     betapsi, 1 )
        !$acc end host_data
        IF ( gstart == 2 ) THEN
          !$acc kernels
          betapsi(:,1) = betapsi(:,1) - beta(1,:)*psi(1,1)
          !$acc end kernels
        END IF
        !
    ELSE
        !
        !$acc host_data use_device(beta,psi,betapsi)
        CALL MYDGEMM( 'C', 'N', nkb, m, 2*npw, 2.0_DP, beta, 2*npwx, psi, &
                    2*npwx, 0.0_DP, betapsi, nkb )
        !$acc end host_data
        IF ( gstart == 2 ) THEN
          !$acc host_data use_device(beta,psi,betapsi)
          CALL MYDGER( nkb, m, -1.0_DP, beta, 2*npwx, psi, 2*npwx, betapsi, nkb )
          !$acc end host_data
        END IF
        !
    ENDIF
    !
    IF (mp_size(intra_bgrp_comm) > 1) THEN
      !$acc host_data use_device(betapsi)
      CALL mp_sum( betapsi( :, 1:m ), intra_bgrp_comm )
      !$acc end host_data
    END IF
    !
    CALL stop_clock( 'calbec' )
    !
    RETURN
    !
  END SUBROUTINE calbec_gamma_acc
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_gamma_cpu ( offload, npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(offload_kind_cpu), INTENT(IN) :: offload
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    REAL (DP), INTENT (out) :: betapsi(:,:)
    INTEGER, INTENT (in) :: npw
    INTEGER, INTENT (in), OPTIONAL :: nbnd
    !
    Call calbec_gamma ( npw, beta, psi, betapsi, nbnd )
    !
  END SUBROUTINE calbec_gamma_cpu
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_gamma ( npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !! matrix times matrix with summation index (k=1,npw) running on
    !! half of the G-vectors or PWs - assuming k=0 is the G=0 component:
    !
    !! $$ betapsi(i,j) = 2Re(\sum_k beta^*(i,k)psi(k,j)) + beta^*(i,0)psi(0,j) $$
    !
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp,        ONLY : mp_sum
    !
    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    REAL (DP), INTENT (out) :: betapsi(:,:)
    INTEGER, INTENT (in) :: npw
    INTEGER, INTENT (in), optional :: nbnd
    !
    INTEGER :: nkb, npwx, m
    !
    IF ( present (nbnd) ) THEN
        m = nbnd
    ELSE
        m = size ( psi, 2)
    ENDIF
    !
    nkb = size (beta, 2)
    IF ( nkb == 0 ) RETURN
    !
    CALL start_clock( 'calbec' )
    IF ( npw == 0 ) betapsi(:,:)=0.0_DP
    npwx= size (beta, 1)
    IF ( npwx /= size (psi, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
#if defined(DEBUG)
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
    CALL mp_sum( betapsi( :, 1:m ), intra_bgrp_comm )
    !
    CALL stop_clock( 'calbec' )
    !
    RETURN
    !
  END SUBROUTINE calbec_gamma
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_k_acc ( offload, npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !! Matrix times matrix with summation index (k=1,npw) running on
    !! G-vectors or PWs:
    !! $$ betapsi(i,j) = \sum_k beta^*(i,k) psi(k,j)$$
    !
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp,       ONLY : mp_sum, mp_size

    IMPLICIT NONE
    TYPE(offload_kind_acc), INTENT(IN) :: offload
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
    IF ( npw == 0 ) THEN
      !$acc kernels
      betapsi(:,:)=(0.0_DP,0.0_DP)
      !$acc end kernels
    END IF
    npwx= size (beta, 1)
    IF ( npwx /= size (psi, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
    IF ( present (nbnd) ) THEN
        m = nbnd
    ELSE
        m = size ( psi, 2)
    ENDIF
#if defined(DEBUG)
    WRITE (*,*) 'calbec k'
    WRITE (*,*)  nkb,  size (betapsi,1) , m , size (betapsi, 2)
#endif
    IF ( nkb /= size (betapsi,1) .or. m > size (betapsi, 2) ) &
      CALL errore ('calbec', 'size mismatch', 3)
    !
    IF ( m == 1 ) THEN
       !
       !$acc host_data use_device(beta, psi, betapsi)
       CALL MYZGEMV( 'C', npw, nkb, (1.0_DP,0.0_DP), beta, npwx, psi, 1, &
                   (0.0_DP, 0.0_DP), betapsi, 1 )
       !$acc end host_data
       !
    ELSE
       !
       !$acc host_data use_device(beta, psi, betapsi)
       CALL MYZGEMM( 'C', 'N', nkb, m, npw, (1.0_DP,0.0_DP), &
                 beta, npwx, psi, npwx, (0.0_DP,0.0_DP), betapsi, nkb )
       !$acc end host_data
       !
    ENDIF
    !
    IF (mp_size(intra_bgrp_comm) > 1) THEN
      !$acc host_data use_device(betapsi)
      CALL mp_sum( betapsi( :, 1:m ), intra_bgrp_comm )
      !$acc end host_data
    END IF
    !
    CALL stop_clock( 'calbec' )
    !
    RETURN
    !
  END SUBROUTINE calbec_k_acc
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_k_cpu ( offload, npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(offload_kind_cpu), INTENT(IN) :: offload
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    COMPLEX (DP), INTENT (out) :: betapsi(:,:)
    INTEGER, INTENT (in) :: npw
    INTEGER, OPTIONAL :: nbnd
    INTEGER :: m
    !
    IF ( present (nbnd) ) THEN
        m = nbnd
    ELSE
        m = size ( psi, 2)
    ENDIF
    Call calbec_k ( npw, beta, psi, betapsi, m )
    !
    RETURN
    !
  END SUBROUTINE calbec_k_cpu
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_k ( npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !! Matrix times matrix with summation index (k=1,npw) running on
    !! G-vectors or PWs:
    !! $$ betapsi(i,j) = \sum_k beta^*(i,k) psi(k,j)$$
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
#if defined(DEBUG)
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
  SUBROUTINE calbec_nc_acc ( offload, npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !! Matrix times matrix with summation index (k below) running on
    !! G-vectors or PWs corresponding to two different polarizations:
    !
    !! * \(betapsi(i,1,j) = \sum_k=1,npw beta^*(i,k) psi(k,j)\)
    !! * \(betapsi(i,2,j) = \sum_k=1,npw beta^*(i,k) psi(k+npwx,j)\)
    !
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp,       ONLY : mp_sum, mp_size

    IMPLICIT NONE
    TYPE(offload_kind_acc), INTENT(IN) :: offload
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
    IF ( npw == 0 ) THEN
      !$acc kernels
      betapsi(:,:,:)=(0.0_DP,0.0_DP)
      !$acc end kernels
    END IF
    npwx= size (beta, 1)
    IF ( 2*npwx /= size (psi, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
    IF ( present (nbnd) ) THEN
        m = nbnd
    ELSE
        m = size ( psi, 2)
    ENDIF
    npol= size (betapsi, 2)
#if defined(DEBUG)
    WRITE (*,*) 'calbec nc'
    WRITE (*,*)  nkb,  size (betapsi,1) , m , size (betapsi, 3)
#endif
    IF ( nkb /= size (betapsi,1) .or. m > size (betapsi, 3) ) &
      CALL errore ('calbec', 'size mismatch', 3)
    !
    !$acc host_data use_device(beta, psi, betapsi)
    CALL MYZGEMM ('C', 'N', nkb, m*npol, npw, (1.0_DP, 0.0_DP), beta, &
              npwx, psi, npwx, (0.0_DP, 0.0_DP),  betapsi, nkb)
    !$acc end host_data
    !
    IF (mp_size(intra_bgrp_comm) > 1) THEN
      !$acc host_data use_device(betapsi)
      CALL mp_sum( betapsi( :, :, 1:m ), intra_bgrp_comm )
      !$acc end host_data
    END IF
    !
    CALL stop_clock( 'calbec' )
    !
    RETURN
    !
  END SUBROUTINE calbec_nc_acc
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_nc_cpu ( offload, npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE(offload_kind_cpu), INTENT(IN) :: offload
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    COMPLEX (DP), INTENT (out) :: betapsi(:,:,:)
    INTEGER, INTENT (in) :: npw
    INTEGER, OPTIONAL :: nbnd
    INTEGER :: m
    !
    IF ( present (nbnd) ) THEN
        m = nbnd
    ELSE
        m = size ( psi, 2)
    ENDIF
    Call calbec_nc ( npw, beta, psi, betapsi, m )
    !
    RETURN
    !
  END SUBROUTINE calbec_nc_cpu
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_nc ( npw, beta, psi, betapsi, nbnd )
    !-----------------------------------------------------------------------
    !! Matrix times matrix with summation index (k below) running on
    !! G-vectors or PWs corresponding to two different polarizations:
    !
    !! * \(betapsi(i,1,j) = \sum_k=1,npw beta^*(i,k) psi(k,j)\)
    !! * \(betapsi(i,2,j) = \sum_k=1,npw beta^*(i,k) psi(k+npwx,j)\)
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
#if defined(DEBUG)
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
  SUBROUTINE allocate_bec_type_acc ( nkb, nbnd, bec, comm )
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    TYPE (bec_type) :: bec
    INTEGER, INTENT (in) :: nkb, nbnd
    INTEGER, INTENT (in), OPTIONAL :: comm
    INTEGER :: ierr
    !
    bec%nbnd = nbnd
    !
    IF( PRESENT( comm ) .AND. gamma_only .AND. smallmem ) THEN
       CALL errore ('allocate_bec_type_acc', 'discontinued feature', 1)
    END IF
    !
    !$acc enter data copyin(bec)
    !
    IF ( gamma_only ) THEN
       !
       ALLOCATE( bec%r( nkb, nbnd ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' allocate_bec_type_acc ', ' cannot allocate bec%r ', ABS(ierr) )
       !
       bec%r(:,:)=0.0D0
       !$acc enter data copyin(bec%r)
       !
    ELSEIF ( noncolin) THEN
       !
       ALLOCATE( bec%nc( nkb, npol, nbnd ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' allocate_bec_type_acc ', ' cannot allocate bec%nc ', ABS(ierr) )
       !
       bec%nc(:,:,:)=(0.0D0,0.0D0)
       !$acc enter data copyin(bec%nc)
       !
    ELSE
       !
       ALLOCATE( bec%k( nkb, nbnd ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' allocate_bec_type_acc ', ' cannot allocate bec%k ', ABS(ierr) )
       !
       bec%k(:,:)=(0.0D0,0.0D0)
       !$acc enter data copyin(bec%k)
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE allocate_bec_type_acc
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE allocate_bec_type ( nkb, nbnd, bec, comm )
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    TYPE (bec_type) :: bec
    INTEGER, INTENT (in) :: nkb, nbnd
    INTEGER, INTENT (in), OPTIONAL :: comm
    INTEGER :: ierr
    !
    bec%nbnd = nbnd
    !
    IF( PRESENT( comm ) .AND. gamma_only .AND. smallmem ) THEN
       CALL errore ('allocate_bec_type', 'discontinued feature', 1)
    END IF
    !
    IF ( gamma_only ) THEN
       !
       ALLOCATE( bec%r( nkb, nbnd ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' allocate_bec_type ', ' cannot allocate bec%r ', ABS(ierr) )
       !
       bec%r(:,:)=0.0D0
       !
    ELSEIF ( noncolin) THEN
       !
       ALLOCATE( bec%nc( nkb, npol, nbnd ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' allocate_bec_type ', ' cannot allocate bec%nc ', ABS(ierr) )
       !
       bec%nc(:,:,:)=(0.0D0,0.0D0)
       !
    ELSE
       !
       ALLOCATE( bec%k( nkb, nbnd ), STAT=ierr )
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
  SUBROUTINE deallocate_bec_type_acc (bec)
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE (bec_type) :: bec
    !
    bec%nbnd = 0
    !
    IF (allocated(bec%r))  THEN
      !$acc exit data delete(bec%r)
      DEALLOCATE(bec%r)
    END IF
    IF (allocated(bec%nc)) THEN
       !$acc exit data delete(bec%nc)
       DEALLOCATE(bec%nc)
    END IF
    IF (allocated(bec%k))  THEN
      !$acc exit data delete(bec%k)
      DEALLOCATE(bec%k)
    END IF
    !
    !$acc exit data delete(bec)
    !
    RETURN
    !
  END SUBROUTINE deallocate_bec_type_acc

  !
  !-----------------------------------------------------------------------
  SUBROUTINE deallocate_bec_type (bec)
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    TYPE (bec_type) :: bec
    !
    bec%nbnd = 0
    !
    IF (allocated(bec%r))  THEN
      DEALLOCATE(bec%r)
    END IF
    IF (allocated(bec%nc)) THEN
       DEALLOCATE(bec%nc)
    END IF
    IF (allocated(bec%k))  THEN
      DEALLOCATE(bec%k)
    END IF
    !
    RETURN
    !
  END SUBROUTINE deallocate_bec_type

  SUBROUTINE beccopy(bec, bec1, nkb, nbnd, comm)
    USE mp, ONLY: mp_size, mp_sum
    IMPLICIT NONE
    TYPE(bec_type), INTENT(in) :: bec
    TYPE(bec_type)  :: bec1
    INTEGER, INTENT(in) :: nkb, nbnd
    INTEGER, INTENT (in), OPTIONAL :: comm

    INTEGER :: nbgrp, ib_start, ib_end, this_bgrp_nbnd

    nbgrp = 1; ib_start = 1; ib_end = nbnd ; this_bgrp_nbnd = nbnd
    IF( PRESENT( comm ) ) THEN
       nbgrp = mp_size( comm )
       call divide( comm, nbnd, ib_start, ib_end) ; this_bgrp_nbnd = ib_end - ib_start + 1
    END IF

    IF (gamma_only) THEN
       if(nbgrp>1) bec1%r = 0.d0
       CALL dcopy(nkb*this_bgrp_nbnd, bec%r, 1, bec1%r(1,ib_start), 1)
       if (nbgrp > 1) CALL mp_sum( bec1%r, comm )
    ELSEIF (noncolin) THEN
       if(nbgrp>1) bec1%nc = ( 0.d0, 0.d0 )
       CALL zcopy(nkb*npol*this_bgrp_nbnd, bec%nc, 1, bec1%nc(1,1,ib_start),  1)
       if (nbgrp > 1) CALL mp_sum( bec1%nc, comm )
    ELSE
       if(nbgrp>1) bec1%k = ( 0.d0, 0.d0 )
       CALL zcopy(nkb*this_bgrp_nbnd, bec%k, 1, bec1%k(1,ib_start), 1)
       if (nbgrp > 1) CALL mp_sum( bec1%k, comm )
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

  SUBROUTINE becupdate_0D_cpu( offload, bec_host, bec_device )
  IMPLICIT NONE
  TYPE(offload_kind_cpu), INTENT(IN)  :: offload 
  TYPE(bec_type),         INTENT(IN)  :: bec_device
  TYPE(bec_type),         INTENT(OUT) :: bec_host
  !
  RETURN
  !
  END SUBROUTINE becupdate_0D_cpu

  SUBROUTINE becupdate_1D_cpu( offload, bec_host, idx, ni, bec_device )
  IMPLICIT NONE
  TYPE(offload_kind_cpu), INTENT(IN)  :: offload
  INTEGER,                INTENT(IN)  :: idx, ni  
  TYPE(bec_type),         INTENT(IN)  :: bec_device
  TYPE(bec_type),         INTENT(OUT) :: bec_host(ni)
  !
  RETURN
  !
  END SUBROUTINE becupdate_1D_cpu

  SUBROUTINE becupdate_2D_cpu( offload, bec_host, idx, ni, jdx, nj, bec_device )
  IMPLICIT NONE
  TYPE(offload_kind_cpu), INTENT(IN)  :: offload
  INTEGER,                INTENT(IN)  :: idx, ni, jdx, nj
  TYPE(bec_type),         INTENT(IN)  :: bec_device
  TYPE(bec_type),         INTENT(OUT) :: bec_host(ni,nj)
  !
  RETURN
  !
  END SUBROUTINE becupdate_2D_cpu

  SUBROUTINE becupdate_0D_acc( offload, bec_host, bec_device )
  IMPLICIT NONE
  TYPE(offload_kind_acc), INTENT(IN)  :: offload 
  TYPE(bec_type),         INTENT(IN)  :: bec_device
  TYPE(bec_type),         INTENT(OUT) :: bec_host
  !
  !$acc data copy(bec_host) present(bec_device)
  !
  if(allocated(bec_device%r)) then
    !$acc kernels copyout(bec_host%r)
    bec_host%r(:,:) = bec_device%r(:,:)
    !$acc end kernels
  elseif(allocated(bec_device%k)) then
    !$acc kernels copyout(bec_host%k)
    bec_host%k(:,:) = bec_device%k(:,:)
    !$acc end kernels
  elseif(allocated(bec_device%nc)) then
    !$acc kernels copyout(bec_host%nc)
    bec_host%nc(:,:,:) = bec_device%nc(:,:,:)
    !$acc end kernels
  endif  
  !
  !$acc end data
  !
  RETURN
  !
  END SUBROUTINE becupdate_0D_acc

  SUBROUTINE becupdate_1D_acc( offload, bec_host, idx, ni, bec_device )
  IMPLICIT NONE
  TYPE(offload_kind_acc), INTENT(IN)  :: offload
  INTEGER,                INTENT(IN)  :: idx, ni  
  TYPE(bec_type),         INTENT(IN)  :: bec_device
  TYPE(bec_type),         INTENT(OUT) :: bec_host(ni)
  !
  !$acc data copy(bec_host(idx)) present(bec_device)
  !
  if(allocated(bec_device%r)) then
    !$acc kernels copyout(bec_host(idx)%r)
    bec_host(idx)%r(:,:) = bec_device%r(:,:)
    !$acc end kernels
  elseif(allocated(bec_device%k)) then
    !$acc kernels copyout(bec_host(idx)%k)
    bec_host(idx)%k(:,:) = bec_device%k(:,:)
    !$acc end kernels
  elseif(allocated(bec_device%nc)) then
    !$acc kernels copyout(bec_host(idx)%nc)
    bec_host(idx)%nc(:,:,:) = bec_device%nc(:,:,:)
    !$acc end kernels
  endif  
  !
  !$acc end data
  !
  RETURN
  !
  END SUBROUTINE becupdate_1D_acc

  SUBROUTINE becupdate_2D_acc( offload, bec_host, idx, ni, jdx, nj, bec_device )
  IMPLICIT NONE
  TYPE(offload_kind_acc), INTENT(IN)  :: offload
  INTEGER,                INTENT(IN)  :: idx, ni, jdx, nj
  TYPE(bec_type),         INTENT(IN)  :: bec_device
  TYPE(bec_type),         INTENT(OUT) :: bec_host(ni,nj)
  !
  !$acc data copy(bec_host(idx,jdx)) present(bec_device)
  !
  if(allocated(bec_device%r)) then
    !$acc kernels copyout(bec_host(idx,jdx)%r)
    bec_host(idx,jdx)%r(:,:) = bec_device%r(:,:)
    !$acc end kernels
  elseif(allocated(bec_device%k)) then
    !$acc kernels copyout(bec_host(idx,jdx)%k)
    bec_host(idx,jdx)%k(:,:) = bec_device%k(:,:)
    !$acc end kernels
  elseif(allocated(bec_device%nc)) then
    !$acc kernels copyout(bec_host(idx,jdx)%nc)
    bec_host(idx,jdx)%nc(:,:,:) = bec_device%nc(:,:,:)
    !$acc end kernels
  endif  
  !
  !$acc end data
  !
  RETURN
  !
  END SUBROUTINE becupdate_2D_acc

END MODULE becmod
