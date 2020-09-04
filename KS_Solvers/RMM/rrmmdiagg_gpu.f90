!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0._DP, 0._DP )
#if !defined(__CUDA)
#define cublasDgemm dgemm
#endif
!
!----------------------------------------------------------------------------
SUBROUTINE rrmmdiagg_gpu( h_psi_gpu, s_psi_gpu, npwx, npw, nbnd, psi, hpsi, spsi, & 
                          e, g2kin, btype, ethr, ndiis, uspp, do_hpsi, is_exx, notconv, rmm_iter )
  !----------------------------------------------------------------------------
  !
  ! ... Iterative diagonalization of a complex hermitian matrix
  ! ... through preconditioned RMM-DIIS algorithm.
  !
#if defined(__CUDA)
  use cudafor
!  use cublas
#endif
  USE util_param,    ONLY : DP, stdout, eps14, eps16
  USE mp,            ONLY : mp_sum, mp_bcast
  USE mp_bands_util, ONLY : gstart, inter_bgrp_comm, intra_bgrp_comm, me_bgrp, root_bgrp, &
                            root_bgrp_id, use_bgrp_in_hpsi
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER,     INTENT(IN)    :: npwx, npw, nbnd
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
  COMPLEX(DP), INTENT(INOUT) :: psi (npwx,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: hpsi(npwx,nbnd)
  COMPLEX(DP), INTENT(INOUT) :: spsi(npwx,nbnd)
  REAL(DP),    INTENT(INOUT) :: e(nbnd)
  REAL(DP),    INTENT(IN)    :: g2kin(npwx)
!****************GPU-VARIABLES***************************************
  COMPLEX(DP), ALLOCATABLE :: psi_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: hpsi_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: spsi_d(:,:)
  REAL(DP),    ALLOCATABLE :: e_d(:)
  REAL(DP),    ALLOCATABLE :: g2kin_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: psi_d, hpsi_d, spsi_d, e_d, g2kin_d
#endif
!********************************************************************
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
  INTEGER                  :: npw2, npwx2
  INTEGER                  :: ibnd, ibnd_start, ibnd_end, ibnd_size
  INTEGER,     ALLOCATABLE :: ibnd_index(:)
  INTEGER,     ALLOCATABLE :: jbnd_index(:)
  LOGICAL,     ALLOCATABLE :: conv(:)
  REAL(DP)                 :: empty_ethr
  !
  COMPLEX(DP), ALLOCATABLE :: phi(:,:,:), hphi(:,:,:), sphi(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: kpsi(:,:), hkpsi(:,:), skpsi(:,:)
  REAL(DP),    ALLOCATABLE :: hr(:,:,:), sr(:,:,:)
  REAL(DP),    ALLOCATABLE :: php(:,:), psp(:,:)
  REAL(DP),    ALLOCATABLE :: ew(:), hw(:), sw(:)
!****************GPU-VARIABLES***************************************
  COMPLEX(DP), ALLOCATABLE :: phi_d(:,:,:), hphi_d(:,:,:), sphi_d(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: kpsi_d(:,:), hkpsi_d(:,:), skpsi_d(:,:)
  REAL(DP),    ALLOCATABLE :: hr_d(:,:,:), sr_d(:,:,:)
  REAL(DP),    ALLOCATABLE :: php_d(:,:), psp_d(:,:)
  REAL(DP),    ALLOCATABLE :: ew_d(:), hw_d(:), sw_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: phi_d, hphi_d, sphi_d, kpsi_d, hkpsi_d, skpsi_d
  attributes(DEVICE) :: hr_d, sr_d, php_d, psp_d
  attributes(DEVICE) :: ew_d, hw_d, sw_d
#endif
!********************************************************************
  INTEGER :: i, j, k
  !
  REAL(DP),    PARAMETER   :: SREF = 0.50_DP
  REAL(DP),    PARAMETER   :: SMIN = 0.05_DP
  REAL(DP),    PARAMETER   :: SMAX = 1.00_DP
  !
  EXTERNAL :: h_psi_gpu, s_psi_gpu
  !EXTERNAL :: h_psi, s_psi
    ! h_psi(npwx,npw,nbnd,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nbnd,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nbnd)
  !
  CALL start_clock( 'rrmmdiagg' )
  !
  IF ( gstart == -1 ) CALL errore( ' rrmmdiagg ', 'gstart variable not initialized', 1 )
  !
  empty_ethr = MAX( ( ethr * 5._DP ), 1.E-5_DP )
  !
  npw2  = 2 * npw
  npwx2 = 2 * npwx
  
  CALL divide( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
  !
  ibnd_size = MAX( ibnd_end - ibnd_start + 1, 0 )
  !
  IF( ibnd_size == 0 ) CALL errore( ' rrmmdiagg ', ' ibnd_size == 0 ', 1 )
  !
  !***GPU***
  ALLOCATE (psi_d(npwx,nbnd)) 
  ALLOCATE (hpsi_d(npwx,nbnd)) 
  ALLOCATE (spsi_d(npwx,nbnd))
  ALLOCATE (e_d(nbnd)) 
  ALLOCATE (g2kin_d(npwx)) 
  hpsi_d = hpsi
  psi_d  = psi
  IF ( uspp ) &
  spsi_d = spsi
  e_d = e
  !g2kin_d = g2kin
  !***GPU***
  !
  ALLOCATE( phi( npwx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  ALLOCATE( phi_d( npwx, ibnd_start:ibnd_end, ndiis ), STAT=ierr ) !***GPU***
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate phi ', ABS(ierr) )
  !
  ALLOCATE( hphi( npwx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  ALLOCATE( hphi_d( npwx, ibnd_start:ibnd_end, ndiis ), STAT=ierr ) !***GPU***
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate hphi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     !
     ALLOCATE( sphi( npwx, ibnd_start:ibnd_end, ndiis ), STAT=ierr )
     ALLOCATE( sphi_d( npwx, ibnd_start:ibnd_end, ndiis ), STAT=ierr ) !***GPU***
     IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate sphi ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( kpsi( npwx, nbnd ), STAT=ierr )
  ALLOCATE( kpsi_d( npwx, nbnd ), STAT=ierr ) !***GPU***
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate kpsi ', ABS(ierr) )
  !
  ALLOCATE( hkpsi( npwx, nbnd ), STAT=ierr )
  ALLOCATE( hkpsi_d( npwx, nbnd ), STAT=ierr ) !***GPU***
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate hkpsi ', ABS(ierr) )
  
  IF ( uspp ) THEN
     !
     ALLOCATE( skpsi( npwx, nbnd ), STAT=ierr )
     ALLOCATE( skpsi_d( npwx, nbnd ), STAT=ierr ) !***GPU***
     IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate skpsi ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( hr( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )
  ALLOCATE( hr_d( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )  !***GPU***
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate hr ', ABS(ierr) )
  !
  ALLOCATE( sr( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )
  ALLOCATE( sr_d( ndiis, ndiis, ibnd_start:ibnd_end ), STAT=ierr )  !***GPU***
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate sr ', ABS(ierr) )
  !
  ALLOCATE( php( ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  ALLOCATE( php_d( ibnd_start:ibnd_end, ndiis ), STAT=ierr )  !***GPU***
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate php ', ABS(ierr) )
  !
  ALLOCATE( psp( ibnd_start:ibnd_end, ndiis ), STAT=ierr )
  ALLOCATE( psp_d( ibnd_start:ibnd_end, ndiis ), STAT=ierr )  !***GPU***
  IF( ierr /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot allocate psp ', ABS(ierr) )
  !
  ALLOCATE( ew( nbnd ) )
  ALLOCATE( ew_d( nbnd ) )  !***GPU***
  ALLOCATE( hw( nbnd ) )
  ALLOCATE( hw_d( nbnd ) )  !***GPU***
  ALLOCATE( sw( nbnd ) )
  ALLOCATE( sw_d( nbnd ) )  !***GPU***
  !
  ALLOCATE( conv( nbnd ) )
  ALLOCATE( ibnd_index( nbnd ) )
  ALLOCATE( jbnd_index( ibnd_start:ibnd_end ) )
  !
  phi  = ZERO
  phi_d  = ZERO  !***GPU***
  hphi = ZERO
  hphi_d = ZERO  !***GPU***
  IF ( uspp ) sphi = ZERO
  IF ( uspp ) sphi_d = ZERO !***GPU***
  !
  kpsi  = ZERO
  kpsi_d  = ZERO  !***GPU***
  hkpsi = ZERO
  hkpsi_d = ZERO  !***GPU***
  IF ( uspp ) skpsi = ZERO
  IF ( uspp ) skpsi_d = ZERO !***GPU***
  !
  hr = 0._DP
  hr_d = 0._DP  !***GPU***
  sr = 0._DP
  sr_d = 0._DP  !***GPU***
  !
  php = 0._DP
  php_d = 0._DP  !***GPU***
  psp = 0._DP
  psp_d = 0._DP   !***GPU***
  !
  ew = e
  hw = e
!$cuf kernel do(1) <<<*,*>>>
  DO i=1, nbnd
     ew_d(i) = e_d(i)   !***GPU***
     hw_d(i) = e_d(i)   !***GPU***
  END DO
  sw = 1._DP
  sw_d = 1._DP  !***GPU***
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
     !CALL calc_hpsi( )
     CALL calc_hpsi_gpu( )
     hpsi   = hpsi_d 
     psi    =  psi_d
     IF ( uspp ) &
     spsi(1:npwx,1:nbnd) = spsi_d(1:npwx,1:nbnd)
     hw(1:nbnd) = hw_d(1:nbnd)
     sw(1:nbnd) = sw_d(1:nbnd)
     ew(1:nbnd) = ew_d(1:nbnd)
     e(1:nbnd)  = e_d(1:nbnd)
     !
  END IF
  !
  ! ... Set Im[ psi(G=0) ] - needed for numerical stability
  !
  IF ( gstart == 2 ) THEN
     !
     psi (1,1:nbnd) = CMPLX( DBLE( psi (1,1:nbnd) ), 0._DP, kind=DP )
     hpsi(1,1:nbnd) = CMPLX( DBLE( hpsi(1,1:nbnd) ), 0._DP, kind=DP )
     IF ( uspp ) &
     spsi(1,1:nbnd) = CMPLX( DBLE( spsi(1,1:nbnd) ), 0._DP, kind=DP )
!$cuf kernel do(1) <<<*,*>>>
     DO i=1, nbnd
        psi_d (1,i) = CMPLX( DBLE( psi_d (1,i) ), 0._DP, kind=DP )!***GPU***
        hpsi_d(1,i) = CMPLX( DBLE( hpsi_d(1,i) ), 0._DP, kind=DP )!***GPU***
        IF ( uspp ) &
        spsi_d(1,i) = CMPLX( DBLE( spsi_d(1,i) ), 0._DP, kind=DP )!***GPU***
     END DO
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
     !CALL do_diis_gpu( idiis )
     WRITE(stdout,*) '******do do_diis*****'
     !
     !psi(1:npwx,1:nbnd)  = psi_d(1:npwx,1:nbnd)
     !hpsi(1:npwx,1:nbnd) = hpsi_d(1:npwx,1:nbnd)
     !IF ( uspp ) &
     !spsi(1:npwx,1:nbnd) = spsi_d(1:npwx,1:nbnd)
     !
     ! ... Line searching
     !
     CALL line_search( )
     WRITE(stdout,*) '******do line_search*****'
     !
     ! ... Calculate eigenvalues and check convergence
     !
     CALL eigenvalues( )
     WRITE(stdout,*) '******do eigenvalues*****'
     !
     IF ( notconv == 0 ) EXIT
     !
  END DO
  WRITE(stdout,*) '******out from rrm do*****'
  !
  rmm_iter = rmm_iter / DBLE( nbnd )
  !
  ! ... Merge wave functions
  !
  IF ( ibnd_start > 1 ) THEN
     !
     psi (:,1:(ibnd_start-1)) = ZERO
     hpsi(:,1:(ibnd_start-1)) = ZERO
     IF ( uspp ) &
     spsi(:,1:(ibnd_start-1)) = ZERO
!$cuf kernel do(2) <<<*,*>>>
     DO i=1,npwx
        DO j=1, ibnd_start-1
           psi_d (i,j) = ZERO
           hpsi_d(i,j) = ZERO
           IF ( uspp ) &
           spsi_d(i,j) = ZERO
        END DO
     END DO
     !
  END IF
  !
  IF ( ibnd_end < nbnd ) THEN
     !
     psi (:,(ibnd_end+1):nbnd) = ZERO
     hpsi(:,(ibnd_end+1):nbnd) = ZERO
     IF ( uspp ) &
     spsi(:,(ibnd_end+1):nbnd) = ZERO
!$cuf kernel do(2) <<<*,*>>>
     DO i=1,npwx
        DO j= ibnd_end+1, nbnd
           psi_d (i,j) = ZERO
           hpsi_d(i,j) = ZERO
           IF ( uspp ) &
           spsi_d(i,j) = ZERO
        END DO
     END DO
     !
  END IF
  !
  CALL mp_sum( psi,  inter_bgrp_comm )
  CALL mp_sum( hpsi, inter_bgrp_comm )
  IF ( uspp ) &
  CALL mp_sum( spsi, inter_bgrp_comm )
  !
  DEALLOCATE( phi )
  DEALLOCATE( hphi )
  IF ( uspp ) DEALLOCATE( sphi )
  DEALLOCATE( kpsi )
  DEALLOCATE( hkpsi )
  IF ( uspp ) DEALLOCATE( skpsi )
  DEALLOCATE( hr )
  DEALLOCATE( sr )
  DEALLOCATE( php )
  DEALLOCATE( psp )
  DEALLOCATE( ew )
  DEALLOCATE( hw )
  DEALLOCATE( sw )
  !***GPU******************
  DEALLOCATE (psi_d) 
  DEALLOCATE (hpsi_d) 
  DEALLOCATE (spsi_d)
  DEALLOCATE (e_d)
  DEALLOCATE (g2kin_d) 
  ! 
  DEALLOCATE( phi_d )
  DEALLOCATE( hphi_d )
  IF ( uspp ) DEALLOCATE( sphi_d )
  DEALLOCATE( kpsi_d )
  DEALLOCATE( hkpsi_d )
  IF ( uspp ) DEALLOCATE( skpsi_d )
  DEALLOCATE( hr_d )
  DEALLOCATE( sr_d )
  DEALLOCATE( php_d )
  DEALLOCATE( psp_d )
  DEALLOCATE( ew_d )
  DEALLOCATE( hw_d )
  DEALLOCATE( sw_d )
  !************************
  DEALLOCATE( conv )
  DEALLOCATE( ibnd_index )
  DEALLOCATE( jbnd_index )
  !
  CALL stop_clock( 'rrmmdiagg' )
  !
  WRITE(stdout,*) '******1diag done*****'
  RETURN
  !
  !
CONTAINS
  !
  SUBROUTINE calc_hpsi_gpu( )
    !
    IMPLICIT NONE
    !
    INTEGER :: ibnd
    !
!    REAL(DP), EXTERNAL :: DDOT
    REAL(DP), EXTERNAL :: gpu_DDOT
    REAL(DP)           :: tmp 
    !
    ! ... Operate the Hamiltonian : H |psi>
    !
    hpsi_d = ZERO
    !
    CALL h_psi_gpu( npwx, npw, nbnd, psi_d, hpsi_d )
    !
    ! ... Operate the Overlap : S |psi>
    !
    IF ( uspp ) THEN
       !
       spsi_d = ZERO
       !
       CALL s_psi_gpu( npwx, npw, nbnd, psi_d, spsi_d )
       !
    END IF
    !
    ! ... Matrix element : <psi| H |psi>
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       tmp = 2._DP * gpu_DDOT( npw2, psi_d(1,ibnd), 1, hpsi_d(1,ibnd), 1 )
       !
       IF ( gstart == 2 ) THEN
          !
          hw_d(ibnd )= tmp - gpu_DDOT( 1, psi_d(1,ibnd), 1, hpsi_d(1,ibnd),1)
        ELSE
          ! 
          hw_d(ibnd) = tmp
          !
       END IF
       !
    END DO
    !
    CALL mp_sum( hw_d(ibnd_start:ibnd_end), intra_bgrp_comm )
    !
    ! ... Matrix element : <psi| S |psi>
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( uspp ) THEN
          !
          tmp = 2._DP * gpu_DDOT( npw2, psi_d(1,ibnd), 1, spsi_d(1,ibnd), 1 )
          !
          IF ( gstart == 2 ) THEN
             !
             sw_d(ibnd )= tmp - gpu_DDOT( 1, psi_d(1,ibnd), 1, spsi_d(1,ibnd),1)
             !
          ELSE
             sw_d(ibnd) = tmp
          END IF
          !
       ELSE
          !
          tmp = 2._DP * gpu_DDOT( npw2, psi_d(1,ibnd), 1, psi_d(1,ibnd), 1 )
          !
          IF ( gstart == 2 ) THEN
             !
             sw_d(ibnd )= tmp - gpu_DDOT( 1, psi_d(1,ibnd), 1, psi_d(1,ibnd),1)
             !
          ELSE
             sw_d(ibnd) = tmp
          END IF
          !
       END IF
       !
    END DO
    !
    CALL mp_sum( sw_d(ibnd_start:ibnd_end), intra_bgrp_comm )
    !
    ! ... Energy eigenvalues
    !
    DO i=ibnd_start, ibnd_end
       tmp= sw_d(i)
       IF(tmp <= eps16 ) CALL errore( ' rrmmdiagg ', ' sw <= 0 ', 1 )
    END DO
    !
!$cuf kernel do(1) <<<*,*>>>
    DO i=1, nbnd
       ew_d(i) = 0._DP
    END DO
    !
!$cuf kernel do(1) <<<*,*>>>
    DO i=ibnd_start,ibnd_end
       ew_d(i) = hw_d(i) / sw_d(i)
    END DO
    !
    CALL mp_sum( ew_d, inter_bgrp_comm )
    !
!$cuf kernel do(1) <<<*,*>>>
    DO i=1, nbnd
       e_d(i) = ew_d(i)
    END DO
    !
    RETURN
    !
  END SUBROUTINE calc_hpsi_gpu
  !
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
    REAL(DP)                 :: er
    COMPLEX(DP), ALLOCATABLE :: vec1(:)
    COMPLEX(DP), ALLOCATABLE :: vec2(:,:)
    REAL(DP),    ALLOCATABLE :: vr(:)
    REAL(DP),    ALLOCATABLE :: tr(:,:)
    !
    !GPU!
    COMPLEX(DP), ALLOCATABLE :: vec1_d(:)
    COMPLEX(DP), ALLOCATABLE :: vec2_d(:,:)
    REAL(DP),    ALLOCATABLE :: vr_d(:)
    REAL(DP),    ALLOCATABLE :: tr_d(:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: vec1_d, vec2_d, vr_d, tr_d
#endif
    !
    ALLOCATE( vec1( npwx ) )
    ALLOCATE( vec2( npwx, idiis ) )
    IF ( idiis > 1 )   ALLOCATE( vr( idiis ) )
    IF ( motconv > 0 ) ALLOCATE( tr( idiis, motconv ) )
    !GPU!
    ALLOCATE( vec1_d( npwx ) )
    ALLOCATE( vec2_d( npwx, idiis ) )
    IF ( idiis > 1 )   ALLOCATE( vr_d( idiis ) )
    IF ( motconv > 0 ) ALLOCATE( tr_d( idiis, motconv ) )
    !
    !
    ! ... Save current wave functions and matrix elements
    !
    !DO ibnd = ibnd_start, ibnd_end
    !   !
    !   IF ( conv(ibnd) ) CYCLE
    !   !
    !   CALL DCOPY( npw2, psi (1,ibnd), 1, phi (1,ibnd,idiis), 1 )
    !   CALL DCOPY( npw2, hpsi(1,ibnd), 1, hphi(1,ibnd,idiis), 1 )
    !   IF ( uspp ) &
    !   CALL DCOPY( npw2, spsi(1,ibnd), 1, sphi(1,ibnd,idiis), 1 )
    !   !
    !   php(ibnd,idiis) = hw(ibnd)
    !   psp(ibnd,idiis) = sw(ibnd)
    !   !
    !END DO
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       CALL gpu_DCOPY( npw2, psi_d (1,ibnd), 1, phi_d (1,ibnd,idiis), 1 )
       CALL gpu_DCOPY( npw2, hpsi_d(1,ibnd), 1, hphi_d(1,ibnd,idiis), 1 )
       IF ( uspp ) &
       CALL gpu_DCOPY( npw2, spsi_d(1,ibnd), 1, sphi_d(1,ibnd,idiis), 1 )
       !
   END DO
   !
!$cuf kernel do(1) <<<*,*>>>
   DO ibnd = ibnd_start, ibnd_end
       php_d(ibnd,idiis) = hw_d(ibnd)
       psp_d(ibnd,idiis) = sw_d(ibnd)
   END DO
   phi  = phi_d
   hphi = hphi_d
   sphi = sphi_d
   php = php_d
   psp = psp_d
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
          er = php_d(ibnd,kdiis)
          !
          CALL gpu_DCOPY( npw2, hphi_d(1,ibnd,kdiis), 1, vec2_d(1,kdiis), 1 )
          !
          IF ( uspp ) THEN
             !
             CALL gpu_DAXPY( npw2, -er, sphi_d(1,ibnd,kdiis), 1, vec2_d(1,kdiis), 1 )
             !
          ELSE
             !
             CALL gpu_DAXPY( npw2, -er, phi_d(1,ibnd,kdiis), 1, vec2_d(1,kdiis), 1 )
             !
          END IF
          !
       END DO
       !
       vec2 = vec2_d
       er = php(ibnd,idiis)
       !
       CALL gpu_DCOPY( npw2, hphi_d(1,ibnd,idiis), 1, vec1_d(1), 1 )
       !
       IF ( uspp ) THEN
          !
          CALL gpu_DAXPY( npw2, -er, sphi_d(1,ibnd,idiis), 1, vec1_d(1), 1 )
          !
       ELSE
          !
          CALL gpu_DAXPY( npw2, -er, phi_d(1,ibnd,idiis), 1, vec1_d(1), 1 )
          !
       END IF
       !
       vec1 = vec1_d
       !
       CALL gpu_DGEMV( 'T', npw2, idiis, 2._DP, vec2_d(1,1), npwx, &
                   vec1_d(1), 1, 0._DP, tr_d(1,jbnd), 1 )
       !
       tr = tr_d
       !
       IF ( gstart == 2 ) &
       tr(1:idiis,jbnd) = tr(1:idiis,jbnd) - DBLE( vec2(1,1:idiis) ) * DBLE( vec1(1) )
       !
    END DO
    !
    IF ( motconv > 0 ) THEN
       !
       CALL mp_sum( tr, intra_bgrp_comm )
       !
    END IF
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       hr(1:idiis,idiis,ibnd) = tr(1:idiis,jbnd)
       hr(idiis,1:idiis,ibnd) = tr(1:idiis,jbnd)
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
       DO kdiis = 1, idiis
          !
          CALL DCOPY( npw2, phi(1,ibnd,kdiis), 1, vec2(1,kdiis), 1 )
          !
       END DO
       !
       IF ( uspp ) THEN
          !
          CALL DCOPY( npw2, sphi(1,ibnd,idiis), 1, vec1(1), 1 )
          !
       ELSE
          !
          CALL DCOPY( npw2, phi(1,ibnd,idiis), 1, vec1(1), 1 )
          !
       END IF
       !
       CALL DGEMV( 'T', npw2, idiis, 2._DP, vec2(1,1), npwx2, &
                   vec1(1), 1, 0._DP, tr(1,jbnd), 1 )
       !
       IF ( gstart == 2 ) &
       tr(1:idiis,jbnd) = tr(1:idiis,jbnd) - DBLE( vec2(1,1:idiis) ) * DBLE( vec1(1) )
       !
    END DO
    !
    IF ( motconv > 0 ) THEN
       !
       CALL mp_sum( tr, intra_bgrp_comm )
       !
    END IF
    !
    DO ibnd = ibnd_start, ibnd_end
       !
       IF ( conv(ibnd) ) CYCLE
       !
       jbnd = jbnd_index(ibnd)
       !
       sr(1:idiis,idiis,ibnd) = tr(1:idiis,jbnd)
       sr(idiis,1:idiis,ibnd) = tr(1:idiis,jbnd)
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
          IF ( me_bgrp == root_bgrp ) CALL diag_diis( ibnd, idiis, vr(:) )
          CALL mp_bcast( vr, root_bgrp, intra_bgrp_comm )
          !
          psi (:,ibnd) = ZERO
          hpsi(:,ibnd) = ZERO
          IF ( uspp ) spsi(:,ibnd) = ZERO
          kpsi(:,kbnd) = ZERO
          !
          DO kdiis = 1, idiis
             !
             ! ... Wave functions
             !
             CALL DAXPY( npw2, vr(kdiis), phi (1,ibnd,kdiis), 1, psi (1,ibnd), 1 )
             CALL DAXPY( npw2, vr(kdiis), hphi(1,ibnd,kdiis), 1, hpsi(1,ibnd), 1 )
             IF ( uspp ) &
             CALL DAXPY( npw2, vr(kdiis), sphi(1,ibnd,kdiis), 1, spsi(1,ibnd), 1 )
             !
             ! ... Residual vectors
             !
             er = php(ibnd,kdiis)
             !
             CALL DCOPY( npw2, hphi(1,ibnd,kdiis), 1, vec1(1), 1 )
             !
             IF ( uspp ) THEN
                !
                CALL DAXPY( npw2, -er, sphi(1,ibnd,kdiis), 1, vec1(1), 1 )
                !
             ELSE
                !
                CALL DAXPY( npw2, -er, phi(1,ibnd,kdiis), 1, vec1(1), 1 )
                !
             END IF
             !
             CALL DAXPY( npw2, vr(kdiis), vec1(1), 1, kpsi(1,kbnd), 1 )
             !
          END DO
          !
       ELSE
          !
          ! ... Wave functions
          !
          norm = SQRT( sw(ibnd) )
          CALL DSCAL( npw2, 1._DP / norm, psi (1,ibnd), 1 )
          CALL DSCAL( npw2, 1._DP / norm, hpsi(1,ibnd), 1 )
          IF ( uspp ) &
          CALL DSCAL( npw2, 1._DP / norm, spsi(1,ibnd), 1 )
          !
          ! ... Residual vectors
          !
          er = hw(ibnd)
          !
          CALL DCOPY( npw2, hpsi(1,ibnd), 1, kpsi(1,kbnd), 1 )
          !
          IF ( uspp ) THEN
             !
             CALL DAXPY( npw2, -er, spsi(1,ibnd), 1, kpsi(1,kbnd), 1 )
             !
          ELSE
             !
             CALL DAXPY( npw2, -er, spsi(1,ibnd), 1, kpsi(1,kbnd), 1 )
             !
          END IF
          !
       END IF
       !
       ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
       IF ( gstart == 2 ) THEN
          !
          psi (1,ibnd) = CMPLX( DBLE( psi (1,ibnd) ), 0._DP, kind=DP )
          hpsi(1,ibnd) = CMPLX( DBLE( hpsi(1,ibnd) ), 0._DP, kind=DP )
          IF ( uspp ) &
          spsi(1,ibnd) = CMPLX( DBLE( spsi(1,ibnd) ), 0._DP, kind=DP )
          kpsi(1,kbnd) = CMPLX( DBLE( kpsi(1,kbnd) ), 0._DP, kind=DP )
          !
       END IF
       !
    END DO
    !
    DEALLOCATE( vec1 )
    DEALLOCATE( vec2 )
    IF ( idiis > 1 )   DEALLOCATE( vr )
    IF ( motconv > 0 ) DEALLOCATE( tr )
    !
    RETURN
    !
  END SUBROUTINE do_diis
  !
  !
  !
  SUBROUTINE diag_diis( ibnd, idiis, vr )
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN)  :: ibnd
    INTEGER,  INTENT(IN)  :: idiis
    REAL(DP), INTENT(OUT) :: vr(idiis)
    !
    INTEGER               :: info
    INTEGER               :: ndim, kdim
    INTEGER               :: i, imin
    REAL(DP)              :: emin
    REAL(DP)              :: vnrm
    REAL(DP), ALLOCATABLE :: h1(:,:)
    REAL(DP), ALLOCATABLE :: h2(:,:)
    REAL(DP), ALLOCATABLE :: h3(:,:)
    REAL(DP), ALLOCATABLE :: s1(:,:)
    REAL(DP), ALLOCATABLE :: x1(:,:)
    REAL(DP), ALLOCATABLE :: u1(:)
    REAL(DP), ALLOCATABLE :: e1(:)
    INTEGER               :: nwork
    REAL(DP), ALLOCATABLE :: work(:)
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
    !
    h1(1:ndim,1:ndim) = hr(1:ndim,1:ndim,ibnd)
    s1(1:ndim,1:ndim) = sr(1:ndim,1:ndim,ibnd)
    !
    CALL DSYEV( 'V', 'U', ndim, s1, ndim, e1, work, nwork, info )
    !
    IF( info /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot solve diis s1', ABS(info) )
    !
    kdim = 0
    !
    x1 = 0._DP
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
       vr        = 0._DP
       vr(idiis) = 1._DP
       !
       GOTO 10
       !
    END IF
    !
    h2 = 0._DP
    !
    CALL DGEMM( 'N', 'N', ndim, kdim, ndim, 1._DP, h1, ndim, x1, ndim, 0._DP, h2, ndim )
    !
    h3 = 0._DP
    !
    CALL DGEMM( 'T', 'N', kdim, kdim, ndim, 1._DP, x1, ndim, h2, ndim, 0._DP, h3, ndim )
    !
    e1 = 0._DP
    !
    CALL DSYEV( 'V', 'U', kdim, h3, ndim, e1, work, nwork, info )
    !
    IF( info /= 0 ) CALL errore( ' rrmmdiagg ', ' cannot solve diis h3', ABS(info) )
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
    CALL DGEMV( 'N', ndim, kdim, 1._DP, x1, ndim, h3(:,imin), 1, 0._DP, vr, 1 )
    !
    s1(1:ndim,1:ndim) = sr(1:ndim,1:ndim,ibnd)
    !
    CALL DGEMV( 'N', ndim, ndim, 1._DP, s1, ndim, vr, 1, 0._DP, u1, 1 )
    !
    vnrm = SQRT( DDOT( ndim, vr, 1, u1, 1 ) )
    !
    vr = vr / vnrm
    !
10  DEALLOCATE( h1 )
    DEALLOCATE( h2 )
    DEALLOCATE( h3 )
    DEALLOCATE( s1 )
    DEALLOCATE( x1 )
    DEALLOCATE( u1 )
    DEALLOCATE( e1 )
    DEALLOCATE( work )
    !
    RETURN
    !
  END SUBROUTINE diag_diis
  !
  !
  !
  SUBROUTINE line_search( )
    !
    IMPLICIT NONE
    !
    INTEGER               :: ig
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
    REAL(DP), ALLOCATABLE :: coef(:,:)
    !
    REAL(DP), EXTERNAL    :: DDOT
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
       ekin(jbnd) = 0._DP
       !
       DO ig = gstart, npw
          !
          psir = DBLE ( psi(ig,ibnd) )
          psii = AIMAG( psi(ig,ibnd) )
          psi2 = psir * psir + psii * psii
          ekin(jbnd) = ekin(jbnd) + 2._DP * g2kin(ig) * psi2
          !
       END DO
       !
       IF ( gstart == 2 ) THEN
          !
          psir = DBLE ( psi(ig,ibnd) )
          psi2 = psir * psir
          ekin(jbnd) = ekin(jbnd) + g2kin(1) * psi2
          !
       END IF
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
       kbnd = ibnd_index(ibnd)
       !
       DO ig = 1, npw
          !
          x  = g2kin(ig) / ( 1.5_DP * ekin(jbnd) )
          x2 = x * x
          x3 = x * x2
          x4 = x * x3
          !
          k1 = 27._DP + 18._DP * x + 12._DP * x2 + 8._DP * x3
          k2 = k1 + 16._DP * x4
          kdiag = ( -4._DP / 3._DP / ekin(jbnd) ) * k1 / k2
          !
          kpsi(ig,kbnd) = kdiag * kpsi(ig,kbnd)
          !
       END DO
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
          IF ( .NOT. conv(ibnd) ) &
          kpsi(:,ibnd_index(ibnd)) = ZERO
          !
       END DO
       !
       DO ibnd = ( ibnd_end + 1 ), nbnd
          !
          IF ( .NOT. conv(ibnd) ) &
          kpsi(:,ibnd_index(ibnd)) = ZERO
          !
       END DO
       !
       CALL mp_sum( kpsi(:,1:notconv), inter_bgrp_comm )
       !
    END IF
    !
    ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
    IF ( gstart == 2 ) &
    kpsi (1,1:notconv) = CMPLX( DBLE( kpsi (1,1:notconv) ), 0._DP, kind=DP )
    !
    ! ... Operate the Hamiltonian : H K (H - eS) |psi>
    !
    !kpsi_d(1:npwx,1:nbnd) = kpsi(1:npwx,1:nbnd)
    !hkpsi_d(1:npwx,1:nbnd) = hkpsi(1:npwx,1:nbnd)
    kpsi_d = kpsi
    hkpsi_d = hkpsi
    CALL h_psi_gpu( npwx, npw, notconv, kpsi_d, hkpsi_d )
    hkpsi = hkpsi_d
    kpsi = kpsi_d
    !hkpsi(1:npwx,1:nbnd) = hkpsi_d(1:npwx,1:nbnd)
    !kpsi(1:npwx,1:nbnd) = kpsi_d(1:npwx,1:nbnd)
    !
    ! ... Operate the Overlap : S K (H - eS) |psi>
    !
    IF ( uspp ) THEN
       kpsi_d  =  kpsi
       skpsi_d =  skpsi
       CALL s_psi_gpu( npwx, npw, notconv, kpsi, skpsi )
       skpsi = skpsi_d
    END IF
    !
    ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
    IF ( gstart == 2 ) THEN
       !
       hkpsi(1,1:nbnd) = CMPLX( DBLE( hkpsi(1,1:nbnd) ), 0._DP, kind=DP )
       IF ( uspp ) &
       skpsi(1,1:nbnd) = CMPLX( DBLE( skpsi(1,1:nbnd) ), 0._DP, kind=DP )
       !
    END IF
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
       php = 2._DP * DDOT( npw2, psi (1,ibnd), 1, hpsi (1,ibnd), 1 )
       khp = 2._DP * DDOT( npw2, kpsi(1,kbnd), 1, hpsi (1,ibnd), 1 )
       khk = 2._DP * DDOT( npw2, kpsi(1,kbnd), 1, hkpsi(1,kbnd), 1 )
       !
       IF ( gstart == 2 ) THEN
          !
          php = php - DBLE( psi (1,ibnd) ) * DBLE ( hpsi (1,ibnd) )
          khp = khp - DBLE( kpsi(1,kbnd) ) * DBLE ( hpsi (1,ibnd) )
          khk = khk - DBLE( kpsi(1,kbnd) ) * DBLE ( hkpsi(1,kbnd) )
          !
       END IF
       !
       IF ( uspp ) THEN
          !
          psp = 2._DP * DDOT( npw2, psi (1,ibnd), 1, spsi (1,ibnd), 1 )
          ksp = 2._DP * DDOT( npw2, kpsi(1,kbnd), 1, spsi (1,ibnd), 1 )
          ksk = 2._DP * DDOT( npw2, kpsi(1,kbnd), 1, skpsi(1,kbnd), 1 )
          !
          IF ( gstart == 2 ) THEN
             !
             psp = psp - DBLE( psi (1,ibnd) ) * DBLE ( spsi (1,ibnd) )
             ksp = ksp - DBLE( kpsi(1,kbnd) ) * DBLE ( spsi (1,ibnd) )
             ksk = ksk - DBLE( kpsi(1,kbnd) ) * DBLE ( skpsi(1,kbnd) )
             !
          END IF
          !
       ELSE
          !
          psp = 2._DP * DDOT( npw2, psi (1,ibnd), 1, psi (1,ibnd), 1 )
          ksp = 2._DP * DDOT( npw2, kpsi(1,kbnd), 1, psi (1,ibnd), 1 )
          ksk = 2._DP * DDOT( npw2, kpsi(1,kbnd), 1, kpsi(1,kbnd), 1 )
          !
          IF ( gstart == 2 ) THEN
             !
             psp = psp - DBLE( psi (1,ibnd) ) * DBLE ( psi (1,ibnd) )
             ksp = ksp - DBLE( kpsi(1,kbnd) ) * DBLE ( psi (1,ibnd) )
             ksk = ksk - DBLE( kpsi(1,kbnd) ) * DBLE ( kpsi(1,kbnd) )
             !
          END IF
          !
       END IF
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
          IF( psp <= eps16 ) CALL errore( ' rrmmdiagg ', ' psp <= 0 ', 1 )
          !
          norm = psp + 2._DP * ksp * SREF + ksk * SREF * SREF
          IF( norm <= eps16 ) CALL errore( ' rrmmdiagg ', ' norm <= 0 ', 1 )
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
          IF( norm <= eps16 ) CALL errore( ' rrmmdiagg ', ' norm <= 0 ', 2 )
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
       c1 = coef(1,jbnd)
       c2 = coef(2,jbnd)
       !
       CALL DSCAL( npw2, c1, psi (1,ibnd), 1 )
       CALL DAXPY( npw2, c2, kpsi(1,kbnd), 1, psi(1,ibnd), 1 )
       !
       CALL DSCAL( npw2, c1, hpsi (1,ibnd), 1 )
       CALL DAXPY( npw2, c2, hkpsi(1,kbnd), 1, hpsi(1,ibnd), 1 )
       !
       IF ( uspp ) THEN
          !
          CALL DSCAL( npw2, c1, spsi (1,ibnd), 1 )
          CALL DAXPY( npw2, c2, skpsi(1,kbnd), 1, spsi(1,ibnd), 1 )
          !
       END IF
       !
       ! NOTE: set Im[ phi(G=0) ] - needed for numerical stability
       IF ( gstart == 2 ) THEN
          !
          psi (1,ibnd) = CMPLX( DBLE( psi (1,ibnd) ), 0._DP, kind=DP )
          hpsi(1,ibnd) = CMPLX( DBLE( hpsi(1,ibnd) ), 0._DP, kind=DP )
          IF ( uspp ) &
          spsi(1,ibnd) = CMPLX( DBLE( spsi(1,ibnd) ), 0._DP, kind=DP )
          !
       END IF
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
    CALL errore( ' rrmmdiagg ', ' sw <= 0 ', 1 )
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
END SUBROUTINE rrmmdiagg_gpu
