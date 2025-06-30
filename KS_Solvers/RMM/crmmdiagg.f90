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
!
!----------------------------------------------------------------------------
SUBROUTINE crmmdiagg( h_psi_ptr, s_psi_ptr, npwx, npw, nbnd, npol, psi, hpsi, spsi, e, &
                      g2kin, btype, ethr, ndiis, uspp, do_hpsi, is_exx, notconv, rmm_iter )
  !----------------------------------------------------------------------------
  !
  ! ... Iterative diagonalization of a complex hermitian matrix
  ! ... through preconditioned RMM-DIIS algorithm.
  !
  USE util_param,     ONLY : DP, eps14, eps16
  USE mp,            ONLY : mp_sum, mp_bcast
  USE mp_bands_util, ONLY : inter_bgrp_comm, intra_bgrp_comm, me_bgrp, root_bgrp, &
                            root_bgrp_id, use_bgrp_in_hpsi
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npwx, npw, nbnd, npol
  COMPLEX(DP), INTENT(INOUT) :: psi (npwx*npol,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: hpsi(npwx*npol,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: spsi(npwx*npol,nbnd)
  REAL(DP),    INTENT(INOUT) :: e(nbnd)
  REAL(DP),    INTENT(IN)    :: g2kin(npwx)
  INTEGER,     INTENT(IN)    :: btype(nbnd)
  REAL(DP),    INTENT(IN)    :: ethr
  INTEGER,     INTENT(IN)    :: ndiis
  LOGICAL,     INTENT(IN)    :: uspp
  LOGICAL,     INTENT(IN)    :: do_hpsi
  LOGICAL,     INTENT(IN)    :: is_exx
  INTEGER,     INTENT(OUT)   :: notconv
  REAL(DP),    INTENT(OUT)   :: rmm_iter
  !
  ! ... local variables
  !
  INTEGER                  :: ierr
  INTEGER                  :: idiis
  INTEGER                  :: motconv
  INTEGER                  :: kdim, kdmx
  INTEGER                  :: ibnd, ibnd_start, ibnd_end, ibnd_size
  INTEGER,     ALLOCATABLE :: ibnd_index(:)
  INTEGER,     ALLOCATABLE :: jbnd_index(:)
  REAL(DP)                 :: empty_ethr
  COMPLEX(DP), ALLOCATABLE :: phi(:,:,:), hphi(:,:,:), sphi(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: kpsi(:,:), hkpsi(:,:), skpsi(:,:)
  COMPLEX(DP), ALLOCATABLE :: hc(:,:,:), sc(:,:,:)
  REAL(DP),    ALLOCATABLE :: php(:,:), psp(:,:)
  REAL(DP),    ALLOCATABLE :: ew(:), hw(:), sw(:)
  LOGICAL,     ALLOCATABLE :: conv(:)
  !
  REAL(DP),    PARAMETER   :: SREF = 0.50_DP
  REAL(DP),    PARAMETER   :: SMIN = 0.05_DP
  REAL(DP),    PARAMETER   :: SMAX = 1.00_DP
  !
  EXTERNAL :: h_psi_ptr, s_psi_ptr
    ! h_psi_ptr(npwx,npw,nbnd,psi,hpsi)
    !     calculates H|psi>
    ! s_psi_ptr(npwx,npw,nbnd,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nbnd)
  !
  CALL start_clock( 'crmmdiagg' )
  !
  empty_ethr = MAX( ( ethr * 5._DP ), 1.E-5_DP )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx * npol
     kdmx = npwx * npol
     !
  END IF
  !
  CALL divide( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
  !
  ibnd_size = MAX( ibnd_end - ibnd_start + 1, 0 )
  !
  IF( ibnd_size == 0 ) CALL errore( ' crmmdiagg ', ' ibnd_size == 0 ', 1 )
  !
  ALLOCATE( phi( kdmx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate phi ', ABS(ierr) )
  !
  ALLOCATE( hphi( kdmx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate hphi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     !
     ALLOCATE( sphi( kdmx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate sphi ', ABS(ierr) )
     !
  END IF
  !$acc enter data create(phi, sphi, hphi)
  !
  ALLOCATE( kpsi( kdmx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate kpsi ', ABS(ierr) )
  !
  ALLOCATE( hkpsi( kdmx, nbnd ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate hkpsi ', ABS(ierr) )
  !$acc enter data create(kpsi, hkpsi)
  !
  IF ( uspp ) THEN
     !
     ALLOCATE( skpsi( kdmx, nbnd ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate skpsi ', ABS(ierr) )
     !$acc enter data create(skpsi)
     !
  END IF
  !
  ALLOCATE( hc( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate hc ', ABS(ierr) )
  !
  ALLOCATE( sc( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate sc ', ABS(ierr) )
  !
  ALLOCATE( php( ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate php ', ABS(ierr) )
  !
  ALLOCATE( psp( ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  IF( ierr /= 0 ) CALL errore( ' crmmdiagg ', ' cannot allocate psp ', ABS(ierr) )
  !
  ALLOCATE( ew( nbnd ) )
  ALLOCATE( hw( nbnd ) )
  ALLOCATE( sw( nbnd ) )
  !
  ALLOCATE( conv( nbnd ) )
  ALLOCATE( ibnd_index( nbnd ) )
  ALLOCATE( jbnd_index( ibnd_start:ibnd_end ) )
  !
  !$acc kernels
  phi  = ZERO
  hphi = ZERO
  IF ( uspp ) sphi = ZERO
  !
  kpsi  = ZERO
  hkpsi = ZERO
  IF ( uspp ) skpsi = ZERO
  !$acc end kernels
  !
  hc = ZERO
  sc = ZERO
  !
  php = 0._DP
  psp = 0._DP
  !
  ew = e
  hw = e
  sw = 1._DP
  !
  conv = .FALSE.
  ibnd_index = 0
  jbnd_index = 0
  !
  FORALL( ibnd = 1:nbnd )              ibnd_index(ibnd) = ibnd
  FORALL( ibnd = ibnd_start:ibnd_end ) jbnd_index(ibnd) = ibnd - ibnd_start + 1
  !
  rmm_iter = 0._DP
  notconv  = nbnd
  motconv  = ibnd_size
  !
  ! ... Calculate H |psi> and S |psi>, if required
  !
  IF ( do_hpsi ) THEN
     !
     CALL calc_hpsi( )
     !
  END IF
  !
  ! ... RMM-DIIS's loop
  !
  DO idiis = 1, ndiis
     !
     rmm_iter = rmm_iter + DBLE( notconv )
     !
     ! ... Perform DIIS
     !
     CALL do_diis( idiis )
     !
     ! ... Line searching
     !
     CALL line_search( )
     !
     ! ... Calculate eigenvalues and check convergence
     !
     CALL eigenvalues( )
     !
     IF ( notconv == 0 ) EXIT
     !
  END DO
  !
  rmm_iter = rmm_iter / DBLE( nbnd )
  !
  ! ... Merge wave functions
  !
  IF ( ibnd_start > 1 ) THEN
     !
     !$acc kernels
     psi (:,1:(ibnd_start-1)) = ZERO
     hpsi(:,1:(ibnd_start-1)) = ZERO
     IF ( uspp ) &
     spsi(:,1:(ibnd_start-1)) = ZERO
     !$acc end kernels
     !
  END IF
  !
  IF ( ibnd_end < nbnd ) THEN
     !
     !$acc kernels
     psi (:,(ibnd_end+1):nbnd) = ZERO
     hpsi(:,(ibnd_end+1):nbnd) = ZERO
     IF ( uspp ) &
     spsi(:,(ibnd_end+1):nbnd) = ZERO
     !$acc end kernels
     !
  END IF
  !
  !$acc host_data use_device(psi, hpsi, spsi)
  CALL mp_sum( psi,  inter_bgrp_comm )
  CALL mp_sum( hpsi, inter_bgrp_comm )
  IF ( uspp ) &
  CALL mp_sum( spsi, inter_bgrp_comm )
  !$acc end host_data
  !
  !$acc exit data delete(phi, sphi, hphi)
  DEALLOCATE( phi )
  DEALLOCATE( hphi )
  IF ( uspp ) DEALLOCATE( sphi )
  !$acc exit data delete(kpsi, hkpsi, skpsi)
  DEALLOCATE( kpsi )
  DEALLOCATE( hkpsi )
  IF ( uspp ) DEALLOCATE( skpsi )
  DEALLOCATE( hc )
  DEALLOCATE( sc )
  DEALLOCATE( php )
  DEALLOCATE( psp )
  DEALLOCATE( ew )
  DEALLOCATE( hw )
  DEALLOCATE( sw )
  DEALLOCATE( conv )
  DEALLOCATE( ibnd_index )
  DEALLOCATE( jbnd_index )
  !
  CALL stop_clock( 'crmmdiagg' )
  !
  RETURN
  !
  !
CONTAINS
  !
  !
  SUBROUTINE calc_hpsi( )
    !
    IMPLICIT NONE
    !
    INTEGER :: ibnd
    !
    REAL(DP), EXTERNAL :: MYDDOT
    !
    ! ... Operate the Hamiltonian : H |psi>
    !
    !$acc kernels
    hpsi = ZERO
    !$acc end kernels
    !
    CALL h_psi_ptr( npwx, npw, nbnd, psi, hpsi )
    !
    ! ... Operate the Overlap : S |psi>
    !
    IF ( uspp ) THEN
       !
       !$acc kernels
       spsi = ZERO
       !$acc end kernels
       !
       CALL s_psi_ptr( npwx, npw, nbnd, psi, spsi )
       !
    END IF
    !
    ! ... Matrix element : <psi| H |psi>
    !
    !$acc host_data use_device(psi, hpsi)
    DO ibnd = ibnd_start, ibnd_end
       !
       hw(ibnd) = MYDDOT( 2*kdim, psi(1,ibnd), 1, hpsi(1,ibnd), 1 )
       !
    END DO
    !$acc end host_data
    !
    CALL mp_sum( hw(ibnd_start:ibnd_end), intra_bgrp_comm )
    !
    ! ... Matrix element : <psi| S |psi>
    !
    !$acc host_data use_device(psi, spsi)
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( uspp ) THEN
          !
          sw(ibnd) = MYDDOT( 2*kdim, psi(1,ibnd), 1, spsi(1,ibnd), 1 )
          !
       ELSE
          !
          sw(ibnd) = MYDDOT( 2*kdim, psi(1,ibnd), 1, psi(1,ibnd), 1 )
          !
       END IF
       !
    END DO
    !$acc end host_data
    !
    CALL mp_sum( sw(ibnd_start:ibnd_end), intra_bgrp_comm )
    !
    ! ... Energy eigenvalues
    !
    IF( ANY( sw(ibnd_start:ibnd_end) <= eps16 ) ) &
    CALL errore( ' crmmdiagg ', ' sw <= 0 ', 1 )
    !
    ew(1:nbnd) = 0._DP
    ew(ibnd_start:ibnd_end) = hw(ibnd_start:ibnd_end) / sw(ibnd_start:ibnd_end)
    !
    CALL mp_sum( ew, inter_bgrp_comm )
    !
    e(1:nbnd) = ew(1:nbnd)
    !
    RETURN
    !
  END SUBROUTINE calc_hpsi
  !
  !
  SUBROUTINE do_diis( idiis )
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: idiis
    !
    INTEGER                  :: ibnd, jbnd, kbnd
    INTEGER                  :: kdiis
    REAL(DP)                 :: norm
    COMPLEX(DP)              :: ec
    COMPLEX(DP), ALLOCATABLE :: vec1(:)
    COMPLEX(DP), ALLOCATABLE :: vec2(:,:)
    COMPLEX(DP), ALLOCATABLE :: vc(:)
    COMPLEX(DP) :: kvc  ! vc(kdiis) for cuf kernel
    COMPLEX(DP), ALLOCATABLE :: tc(:,:)
    !
    ALLOCATE( vec1( kdmx ) )
    ALLOCATE( vec2( kdmx, idiis ) )
    IF ( idiis > 1 )   ALLOCATE( vc( idiis ) )
    IF ( motconv > 0 ) ALLOCATE( tc( idiis, motconv ) )
    !$acc enter data create( vec1, vec2, tc )
    !
    ! ... Save current wave functions and matrix elements
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       !$acc host_data use_device(psi, hpsi, spsi, phi, hphi, sphi )
       CALL MYZCOPY( kdim, psi (1,ibnd), 1, phi (1,ibnd,idiis), 1 )
       CALL MYZCOPY( kdim, hpsi(1,ibnd), 1, hphi(1,ibnd,idiis), 1 )
       IF ( uspp ) &
       CALL MYZCOPY( kdim, spsi(1,ibnd), 1, sphi(1,ibnd,idiis), 1 )
       !$acc end host_data
       !
       php(ibnd,idiis) = hw(ibnd)
       psp(ibnd,idiis) = sw(ibnd)
       !
    END DO
    !
    ! ... <R_i|R_j>
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       ! ... Residual vectors : |R> = (H - e S) |psi>
       !
       DO kdiis = 1, idiis
          !
          ec = CMPLX( php(ibnd,kdiis), 0._DP, kind=DP )
          !
          !$acc host_data use_device(hphi, sphi, phi, vec2)
          CALL MYZCOPY( kdim, hphi(1,ibnd,kdiis), 1, vec2(1,kdiis), 1 )
          !
          IF ( uspp ) THEN
             !
             CALL MYZAXPY( kdim, -ec, sphi(1,ibnd,kdiis), 1, vec2(1,kdiis), 1 )
             !
          ELSE
             !
             CALL MYZAXPY( kdim, -ec, phi(1,ibnd,kdiis), 1, vec2(1,kdiis), 1 )
             !
          END IF
          !$acc end host_data
          !
       END DO
       !
       ec = CMPLX( php(ibnd,idiis), 0._DP, kind=DP )
       !
       !$acc host_data use_device(hphi, sphi, phi, vec1)
       CALL MYZCOPY( kdim, hphi(1,ibnd,idiis), 1, vec1(1), 1 )
       !
       IF ( uspp ) THEN
          !
          CALL MYZAXPY( kdim, -ec, sphi(1,ibnd,idiis), 1, vec1(1), 1 )
          !
       ELSE
          !
          CALL MYZAXPY( kdim, -ec, phi(1,ibnd,idiis), 1, vec1(1), 1 )
          !
       END IF
       !$acc end host_data
       !
       !$acc host_data use_device(vec1, vec2, tc)
       CALL MYZGEMV( 'C', kdim, idiis, ONE, vec2(1,1), kdmx, &
                   vec1(1), 1, ZERO, tc(1,jbnd), 1 )
       !$acc end host_data
       !
    END DO
    !
    IF ( motconv > 0 ) THEN
       !
       !$acc host_data use_device(tc)
       CALL mp_sum( tc, intra_bgrp_comm )
       !$acc end host_data
       !
    END IF
    !
    !$acc update self(tc) 
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       hc(1:idiis,idiis,ibnd) = tc(1:idiis,jbnd)
       hc(idiis,1:idiis,ibnd) = CONJG( tc(1:idiis,jbnd) )
       hc(idiis,idiis,ibnd)   = CMPLX( DBLE( tc(idiis,jbnd) ), 0._DP, kind=DP )
       !
    END DO
    !
    ! ... <phi_i| S |phi_j>
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       !$acc host_data use_device(phi, sphi, hphi, vec1, vec2)
       DO kdiis = 1, idiis
          !
          CALL MYZCOPY( kdim, phi(1,ibnd,kdiis), 1, vec2(1,kdiis), 1 )
          !
       END DO
       !
       IF ( uspp ) THEN
          !
          CALL MYZCOPY( kdim, sphi(1,ibnd,idiis), 1, vec1(1), 1 )
          !
       ELSE
          !
          CALL MYZCOPY( kdim, phi(1,ibnd,idiis), 1, vec1(1), 1 )
          !
       END IF
       !$acc end host_data
       !
       !$acc host_data use_device(vec1, vec2, tc)
       CALL MYZGEMV( 'C', kdim, idiis, ONE, vec2(1,1), kdmx, &
                   vec1(1), 1, ZERO, tc(1,jbnd), 1 )
       !$acc end host_data
       !
    END DO
    !
    IF ( motconv > 0 ) THEN
       !
       !$acc host_data use_device(tc)
       CALL mp_sum( tc, intra_bgrp_comm )
       !$acc end host_data
       !
    END IF
    !
    !$acc update self(tc) 
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       sc(1:idiis,idiis,ibnd) = tc(1:idiis,jbnd)
       sc(idiis,1:idiis,ibnd) = CONJG( tc(1:idiis,jbnd) )
       sc(idiis,idiis,ibnd)   = CMPLX( DBLE( tc(idiis,jbnd) ), 0._DP, kind=DP )
       !
    END DO
    !
    ! ... Update current wave functions and residual vectors
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       kbnd = ibnd_index(ibnd)
       !
       IF ( idiis > 1 ) THEN
          !
          ! ... solve Rv = eSv
          !
          IF ( me_bgrp == root_bgrp ) CALL diag_diis( ibnd, idiis, vc(:) )
          CALL mp_bcast( vc, root_bgrp, intra_bgrp_comm )
          !
          !$acc kernels
          psi (:,ibnd) = ZERO
          hpsi(:,ibnd) = ZERO
          IF ( uspp ) spsi(:,ibnd) = ZERO
          kpsi(:,kbnd) = ZERO
          !$acc end kernels
          !
          DO kdiis = 1, idiis
             !
             ! ... Wave functions
             !
             kvc = vc(kdiis) 
             !
             !$acc host_data use_device(psi, hpsi, spsi, phi, hphi, sphi )
             CALL MYZAXPY( kdim, kvc, phi (1,ibnd,kdiis), 1, psi (1,ibnd), 1 )
             CALL MYZAXPY( kdim, kvc, hphi(1,ibnd,kdiis), 1, hpsi(1,ibnd), 1 )
             IF (uspp) CALL MYZAXPY( kdim, kvc, sphi(1,ibnd,kdiis), 1, spsi(1,ibnd), 1 )
             !$acc end host_data
             !
             ! ... Residual vectors
             !
             ec = CMPLX( php(ibnd,kdiis), 0._DP, kind=DP )
             !
             !$acc host_data use_device(phi, hphi, sphi, vec1 )
             CALL MYZCOPY( kdim, hphi(1,ibnd,kdiis), 1, vec1(1), 1 )
             !
             IF ( uspp ) THEN
                !
                CALL MYZAXPY( kdim, -ec, sphi(1,ibnd,kdiis), 1, vec1(1), 1 )
                !
             ELSE
                !
                CALL MYZAXPY( kdim, -ec, phi(1,ibnd,kdiis), 1, vec1(1), 1 )
                !
             END IF
             !$acc end host_data
             !
             !$acc host_data use_device(kpsi, vec1)
             CALL MYZAXPY( kdim, kvc, vec1(1), 1, kpsi(1,kbnd), 1 )
             !$acc end host_data
             !
          END DO
          !
       ELSE
          !
          ! ... Wave functions
          !
          norm = SQRT( sw(ibnd) )
          !$acc host_data use_device(psi, hpsi, spsi)
          CALL MYZDSCAL( kdim, 1._DP / norm, psi (1,ibnd), 1 )
          CALL MYZDSCAL( kdim, 1._DP / norm, hpsi(1,ibnd), 1 )
          IF ( uspp ) &
          CALL MYZDSCAL( kdim, 1._DP / norm, spsi(1,ibnd), 1 )
          !$acc end host_data
          !
          ! ... Residual vectors
          !
          ec = CMPLX( hw(ibnd), 0._DP, kind=DP )
          !
          !$acc host_data use_device(psi, spsi, hpsi, kpsi)
          CALL MYZCOPY( kdim, hpsi(1,ibnd), 1, kpsi(1,kbnd), 1 )
          !
          IF ( uspp ) THEN
             !
             CALL MYZAXPY( kdim, -ec, spsi(1,ibnd), 1, kpsi(1,kbnd), 1 )
             !
          ELSE
             !
             CALL MYZAXPY( kdim, -ec, psi(1,ibnd), 1, kpsi(1,kbnd), 1 )
             !
          END IF
          !$acc end host_data
          !
       END IF
       !
    END DO
    !
    !$acc exit data delete( vec1, vec2, tc)
    DEALLOCATE( vec1 )
    DEALLOCATE( vec2 )
    IF ( idiis > 1 )   DEALLOCATE( vc )
    IF ( motconv > 0 ) DEALLOCATE( tc )
    !
    RETURN
    !
  END SUBROUTINE do_diis
  !
  !
  SUBROUTINE diag_diis( ibnd, idiis, vc )
    !
    IMPLICIT NONE
    !
    INTEGER,     INTENT(IN)  :: ibnd
    INTEGER,     INTENT(IN)  :: idiis
    COMPLEX(DP), INTENT(OUT) :: vc(idiis)
    !
    INTEGER                  :: info
    INTEGER                  :: ndim, kdim
    INTEGER                  :: i, imin
    REAL(DP)                 :: emin
    REAL(DP)                 :: vnrm
    COMPLEX(DP), ALLOCATABLE :: h1(:,:)
    COMPLEX(DP), ALLOCATABLE :: h2(:,:)
    COMPLEX(DP), ALLOCATABLE :: h3(:,:)
    COMPLEX(DP), ALLOCATABLE :: s1(:,:)
    COMPLEX(DP), ALLOCATABLE :: x1(:,:)
    COMPLEX(DP), ALLOCATABLE :: u1(:)
    REAL(DP),    ALLOCATABLE :: e1(:)
    INTEGER                  :: nwork
    COMPLEX(DP), ALLOCATABLE :: work(:)
    REAL(DP),    ALLOCATABLE :: rwork(:)
    !
    REAL(DP), EXTERNAL    :: DDOT
    !
    ndim  = idiis
    nwork = 3 * ndim
    !
    ALLOCATE( h1( ndim, ndim ) )
    ALLOCATE( h2( ndim, ndim ) )
    ALLOCATE( h3( ndim, ndim ) )
    ALLOCATE( s1( ndim, ndim ) )
    ALLOCATE( x1( ndim, ndim ) )
    ALLOCATE( u1( ndim ) )
    ALLOCATE( e1( ndim ) )
    ALLOCATE( work( nwork ) )
    ALLOCATE( rwork( 3 * ndim - 2 ) )
    !
    h1(1:ndim,1:ndim) = hc(1:ndim,1:ndim,ibnd)
    s1(1:ndim,1:ndim) = sc(1:ndim,1:ndim,ibnd)
    !
    CALL ZHEEV( 'V', 'U', ndim, s1, ndim, e1, work, nwork, rwork, info )
    !
    IF( info /= 0 ) CALL errore( ' crmmdiagg ', ' cannot solve diis ', ABS(info) )
    !
    kdim = 0
    !
    x1 = ZERO
    !
    DO i = 1, ndim
       !
       IF ( e1(i) > eps14 ) THEN
          !
          kdim = kdim + 1
          !
          x1(:,kdim) = s1(:,i) / SQRT(e1(i))
          !
       END IF
       !
    END DO
    !
    IF ( kdim <= 1 ) THEN
       !
       vc        = ZERO
       vc(idiis) = ONE
       !
       GOTO 10
       !
    END IF
    !
    h2 = ZERO
    !
    CALL ZGEMM( 'N', 'N', ndim, kdim, ndim, ONE, h1, ndim, x1, ndim, ZERO, h2, ndim )
    !
    h3 = ZERO
    !
    CALL ZGEMM( 'C', 'N', kdim, kdim, ndim, ONE, x1, ndim, h2, ndim, ZERO, h3, ndim )
    !
    e1 = 0._DP
    !
    CALL ZHEEV( 'V', 'U', kdim, h3, ndim, e1, work, nwork, rwork, info )
    !
    IF( info /= 0 ) CALL errore( ' crmmdiagg ', ' cannot solve diis ', ABS(info) )
    !
    imin = 1
    emin = e1(1)
    !
    DO i = 2, kdim
       !
       IF ( ABS( e1(i) ) < ABS( emin ) ) THEN
          !
          imin = i
          emin = e1(i)
          !
       END IF
       !
    END DO
    !
    CALL ZGEMV( 'N', ndim, kdim, ONE, x1, ndim, h3(:,imin), 1, ZERO, vc, 1 )
    !
    s1(1:ndim,1:ndim) = sc(1:ndim,1:ndim,ibnd)
    !
    CALL ZGEMV( 'N', ndim, ndim, ONE, s1, ndim, vc, 1, ZERO, u1, 1 )
    !
    vnrm = SQRT( DDOT( 2*ndim, vc, 1, u1, 1 ) ) 
    !
    vc = vc / vnrm
    !
10  DEALLOCATE( h1 )
    DEALLOCATE( h2 )
    DEALLOCATE( h3 )
    DEALLOCATE( s1 )
    DEALLOCATE( x1 )
    DEALLOCATE( u1 )
    DEALLOCATE( e1 )
    DEALLOCATE( work )
    DEALLOCATE( rwork )
    !
    RETURN
    !
  END SUBROUTINE diag_diis
  !
  !
  SUBROUTINE line_search( )
    !
    IMPLICIT NONE
    !
    INTEGER               :: ig, ipol
    INTEGER               :: ibnd, jbnd, kbnd
    LOGICAL               :: para_hpsi
    REAL(DP)              :: psir, psii, psi2
    REAL(DP)              :: kdiag, k1, k2
    REAL(DP)              :: x, x2, x3, x4
    REAL(DP), ALLOCATABLE :: ekin(:)
    REAL(DP)              :: a, b
    REAL(DP)              :: ene0, ene1
    REAL(DP)              :: step, norm
    REAL(DP)              :: php, khp, khk
    REAL(DP)              :: psp, ksp, ksk
    REAL(DP), ALLOCATABLE :: hmat(:,:), smat(:,:)
    REAL(DP), ALLOCATABLE :: heig(:), seig(:)
    REAL(DP)              :: c1, c2
    COMPLEX(DP)           :: z1, z2
    REAL(DP), ALLOCATABLE :: coef(:,:)
    !
    REAL(DP), EXTERNAL :: MYDDOT
    !
    REAL(DP) :: ekinj
    INTEGER :: idx
    !
    IF ( motconv > 0 ) THEN
       !
       ALLOCATE( ekin( motconv ) )
       ALLOCATE( hmat( 3, motconv ) )
       ALLOCATE( smat( 3, motconv ) )
       ALLOCATE( heig( motconv ) )
       ALLOCATE( seig( motconv ) )
       ALLOCATE( coef( 2, motconv ) )
       !
    END IF
    !
    ! ... Kinetic energy
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       ekinj = 0._DP
       !
       !$acc kernels copy(ekinj)
       DO ipol = 1, npol
          !
          DO ig = 1, npw
             !
             psir = DBLE ( psi(ig+(ipol-1)*npwx,ibnd) )
             psii = AIMAG( psi(ig+(ipol-1)*npwx,ibnd) )
             psi2 = psir * psir + psii * psii
             ekinj = ekinj + g2kin(ig) * psi2
             !
          END DO
          !
       END DO
       !$acc end kernels
       !
       ekin(jbnd) = ekinj
       !
    END DO
    !
    IF ( motconv > 0 ) THEN
       !
       CALL mp_sum( ekin, intra_bgrp_comm )
       !
    END IF
    !
    ! ... Preconditioning vectors : K (H - e S) |psi>
    !
    ! ... G.Kresse and J.Furthmuller, PRB 54, 11169 (1996)
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       ekinj = ekin(jbnd) 
       !
       kbnd = ibnd_index(ibnd)
       !
       !$acc kernels
       DO ipol = 1, npol
          !
          DO ig = 1, npw
             !
             x  = g2kin(ig) / ( 1.5_DP * ekinj )
             x2 = x * x
             x3 = x * x2
             x4 = x * x3
             !
             k1 = 27._DP + 18._DP * x + 12._DP * x2 + 8._DP * x3
             k2 = k1 + 16._DP * x4
             kdiag = ( -4._DP / 3._DP / ekinj ) * k1 / k2
             !
             kpsi(ig+(ipol-1)*npwx,kbnd) = kdiag * kpsi(ig+(ipol-1)*npwx,kbnd)
             !
          END DO
          !
       END DO
       !$acc end kernels
       !
    END DO
    !
    ! ... Share kpsi for all band-groups
    !
    IF ( use_bgrp_in_hpsi .AND. ( .NOT. is_exx ) .AND. ( notconv > 1 ) ) THEN
       !
       para_hpsi = .TRUE.
       !
    ELSE
       !
       para_hpsi = .FALSE.
       !
    END IF
    !
    IF ( ( .NOT. para_hpsi ) .OR. ( notconv /= nbnd ) ) THEN
       !
       DO ibnd = 1, ( ibnd_start - 1)
          !
          idx = ibnd_index(ibnd) 
          !
          IF ( .NOT. conv(ibnd) ) THEN 
            !$acc kernels
            kpsi(:,idx) = ZERO
            !$acc end kernels
          END IF 
          !
       END DO
       !
       DO ibnd = ( ibnd_end + 1 ), nbnd
          !
          idx = ibnd_index(ibnd) 
          !
          IF ( .NOT. conv(ibnd) ) THEN
            !$acc kernels
            kpsi(:,idx) = ZERO
            !$acc end kernels
          END IF 
          !
       END DO
       !
       !$acc host_data use_device(kpsi)
       CALL mp_sum( kpsi(:,1:notconv), inter_bgrp_comm )
       !$acc end host_data
       !
    END IF
    !
    ! ... Operate the Hamiltonian : H K (H - eS) |psi>
    !
    CALL h_psi_ptr( npwx, npw, notconv, kpsi, hkpsi )
    !
    ! ... Operate the Overlap : S K (H - eS) |psi>
    !
    IF ( uspp ) CALL s_psi_ptr( npwx, npw, notconv, kpsi, skpsi )
    !
    ! ... Create 2 x 2 matrix
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       kbnd = ibnd_index(ibnd)
       !
       !$acc host_data use_device(psi, hpsi, spsi, kpsi, hkpsi, skpsi)
       php = MYDDOT( 2*kdim, psi (1,ibnd), 1, hpsi (1,ibnd), 1 ) 
       khp = MYDDOT( 2*kdim, kpsi(1,kbnd), 1, hpsi (1,ibnd), 1 ) 
       khk = MYDDOT( 2*kdim, kpsi(1,kbnd), 1, hkpsi(1,kbnd), 1 ) 
       !
       IF ( uspp ) THEN
          !
          psp = MYDDOT(2*kdim, psi (1,ibnd), 1, spsi (1,ibnd), 1 ) 
          ksp = MYDDOT(2*kdim, kpsi(1,kbnd), 1, spsi (1,ibnd), 1 ) 
          ksk = MYDDOT(2*kdim, kpsi(1,kbnd), 1, skpsi(1,kbnd), 1 ) 
          !
       ELSE
          !
          psp = MYDDOT( 2*kdim, psi (1,ibnd), 1, psi (1,ibnd), 1 ) 
          ksp = MYDDOT( 2*kdim, kpsi(1,kbnd), 1, psi (1,ibnd), 1 ) 
          ksk = MYDDOT( 2*kdim, kpsi(1,kbnd), 1, kpsi(1,kbnd), 1 ) 
          !
       END IF
       !$acc end host_data
       !
       hmat(1,jbnd) = php
       hmat(2,jbnd) = khp
       hmat(3,jbnd) = khk
       !
       smat(1,jbnd) = psp
       smat(2,jbnd) = ksp
       smat(3,jbnd) = ksk
       !
    END DO
    !
    IF ( motconv > 0 ) THEN
       !
       CALL mp_sum( hmat, intra_bgrp_comm )
       CALL mp_sum( smat, intra_bgrp_comm )
       !
    END IF
    !
    ! ... Line searching for each band
    !
    IF ( me_bgrp == root_bgrp ) THEN
       !
       DO ibnd = ibnd_start, ibnd_end
          !
          IF ( conv(ibnd) ) CYCLE
          !
          jbnd = jbnd_index(ibnd)
          !
          php = hmat(1,jbnd)
          khp = hmat(2,jbnd)
          khk = hmat(3,jbnd)
          !
          psp = smat(1,jbnd)
          ksp = smat(2,jbnd)
          ksk = smat(3,jbnd)
          IF( psp <= eps16 ) CALL errore( ' crmmdiagg ', ' psp <= 0 ', 1 )
          !
          norm = psp + 2._DP * ksp * SREF + ksk * SREF * SREF
          IF( norm <= eps16 ) CALL errore( ' crmmdiagg ', ' norm <= 0 ', 1 )
          !
          ene0 = php / psp
          ene1 = ( php + 2._DP * khp * SREF + khk * SREF * SREF ) / norm
          !
          a = 2._DP * ( khp * psp - php * ksp ) / psp / psp
          b = ( ene1 - ene0 - a * SREF ) / SREF / SREF
          !
          IF( ABS( b ) > eps16 ) THEN
             step  = -0.5_DP * a / b
          ELSE
             IF ( a < 0._DP ) THEN
                step = SMAX
             ELSE
                step = SMIN
             END IF
          END IF
          !
          step  = MAX( SMIN, step )
          step  = MIN( SMAX, step )
          !
          norm  = psp + 2._DP * ksp * step + ksk * step * step
          IF( norm <= eps16 ) CALL errore( ' crmmdiagg ', ' norm <= 0 ', 2 )
          norm  = SQRT( norm )
          !
          coef(1,jbnd) = 1._DP / norm
          coef(2,jbnd) = step  / norm
          !
          ! ... Update current matrix elements
          !
          c1 = coef(1,jbnd)
          c2 = coef(2,jbnd)
          !
          heig(jbnd) = php * c1 * c1 + 2._DP * khp * c1 * c2 + khk * c2 * c2
          seig(jbnd) = psp * c1 * c1 + 2._DP * ksp * c1 * c2 + ksk * c2 * c2
          !
       END DO
       !
    END IF
    !
    IF ( motconv > 0 ) THEN
       !
       CALL mp_bcast( coef, root_bgrp, intra_bgrp_comm )
       CALL mp_bcast( heig, root_bgrp, intra_bgrp_comm )
       CALL mp_bcast( seig, root_bgrp, intra_bgrp_comm )
       !
    END IF
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       hw(ibnd) = heig(jbnd)
       sw(ibnd) = seig(jbnd)
       !
    END DO
    !
    ! ... Update current wave functions
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       kbnd = ibnd_index(ibnd)
       !
       z1 = CMPLX( coef(1,jbnd), 0._DP, kind=DP )
       z2 = CMPLX( coef(2,jbnd), 0._DP, kind=DP )
       !
       !$acc host_data use_device(psi, hpsi, spsi, kpsi, hkpsi, skpsi)
       CALL MYZSCAL( kdim, z1, psi (1,ibnd), 1 )
       CALL MYZAXPY( kdim, z2, kpsi(1,kbnd), 1, psi(1,ibnd), 1 )
       !
       CALL MYZSCAL( kdim, z1, hpsi (1,ibnd), 1 )
       CALL MYZAXPY( kdim, z2, hkpsi(1,kbnd), 1, hpsi(1,ibnd), 1 )
       !
       IF ( uspp ) THEN
          !
          CALL MYZSCAL( kdim, z1, spsi (1,ibnd), 1 )
          CALL MYZAXPY( kdim, z2, skpsi(1,kbnd), 1, spsi(1,ibnd), 1 )
          !
       END IF
       !$acc end host_data
       !
    END DO
    !
    IF ( motconv > 0 ) THEN
       !
       DEALLOCATE( ekin )
       DEALLOCATE( hmat )
       DEALLOCATE( smat )
       DEALLOCATE( heig )
       DEALLOCATE( seig )
       DEALLOCATE( coef )
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE line_search
  !
  !
  SUBROUTINE eigenvalues( )
    !
    IMPLICIT NONE
    !
    INTEGER :: ibnd
    !
    ! ... Energy eigenvalues
    !
    IF( ANY( sw(ibnd_start:ibnd_end) <= eps16 ) ) &
    CALL errore( ' crmmdiagg ', ' sw <= 0 ', 1 )
    !
    ew(1:nbnd) = 0._DP
    ew(ibnd_start:ibnd_end) = hw(ibnd_start:ibnd_end) / sw(ibnd_start:ibnd_end)
    !
    CALL mp_sum( ew, inter_bgrp_comm )
    !
    ! ... Check convergence
    !
    WHERE( btype(1:nbnd) == 1 )
       !
       conv(1:nbnd) = conv(1:nbnd) .OR. ( ( ABS( ew(1:nbnd) - e(1:nbnd) ) < ethr ) )
       !
    ELSEWHERE
       !
       conv(1:nbnd) = conv(1:nbnd) .OR. ( ( ABS( ew(1:nbnd) - e(1:nbnd) ) < empty_ethr ) )
       !
    END WHERE
    !
    CALL mp_bcast( conv, root_bgrp_id, inter_bgrp_comm )
    !
    ! ... Count not converged bands
    !
    notconv = 0
    !
    DO ibnd = 1, nbnd
       !
       IF ( conv(ibnd) ) THEN
          !
          ibnd_index(ibnd) = 0
          !
       ELSE
          !
          notconv = notconv + 1
          !
          ibnd_index(ibnd) = notconv
          !
       END IF
       !
    END DO
    !
    motconv = 0
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) THEN
          !
          jbnd_index(ibnd) = 0
          !
       ELSE
          !
          motconv = motconv + 1
          !
          jbnd_index(ibnd) = motconv
          !
       END IF
       !
    END DO
    !
    ! ... Save current eigenvalues
    !
    e(1:nbnd) = ew(1:nbnd)
    !
    RETURN
    !
  END SUBROUTINE eigenvalues
  !
  !
END SUBROUTINE crmmdiagg
