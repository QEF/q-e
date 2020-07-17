!
! Copyright (C) 2004-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE xc_lda_lsda
!
USE kinds,     ONLY: DP
USE funct,     ONLY: get_iexch, get_icorr, is_libxc,  &
                     exx_is_active, get_exx_fraction, &
                     get_finite_size_cell_volume
!
IMPLICIT NONE
!
PRIVATE
SAVE
!
!  LDA and LSDA exchange-correlation drivers
PUBLIC :: xc, xc_lda, xc_lsda
PRIVATE :: xc_lda_gpu, xc_lsda_gpu
PUBLIC :: change_threshold_lda
!
!  density threshold (set to default value)
REAL(DP) :: rho_threshold = 1.E-10_DP
!
 CONTAINS
!
!
!-----------------------------------------------------------------------
SUBROUTINE change_threshold_lda( rho_thr_in )
  !--------------------------------------------------------------------
  !! Change rho threshold.
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rho_thr_in
  !
  rho_threshold = rho_thr_in
  !
  RETURN
  !
END SUBROUTINE
!
!
!---------------------------------------------------------------------------
SUBROUTINE xc( length, sr_d, sv_d, rho_in, ex_out, ec_out, vx_out, vc_out )
  !-------------------------------------------------------------------------
  !! Wrapper routine. Calls xc-driver routines from internal libraries
  !! of q-e or from the external libxc, depending on the input choice.
  !
  !! NOTE: look at 'PP/src/benchmark_libxc.f90' to see the differences
  !!       between q-e and libxc.
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
#endif
  USE control_flags, ONLY: use_gpu
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: sr_d
  !! spin dimension of rho
  INTEGER, INTENT(IN) :: sv_d
  !! spin dimension of v
  REAL(DP), INTENT(IN) :: rho_in(length,sr_d)
  !! Charge density
  REAL(DP), INTENT(OUT) :: ex_out(length)
  !! \(\epsilon_x(rho)\) ( NOT \(E_x(\text{rho})\) )
  REAL(DP), INTENT(OUT) :: vx_out(length,sv_d)
  !! \(dE_x(\text{rho})/d\text{rho}  ( NOT d\epsilon_x(\text{rho})/d\text{rho} )
  REAL(DP), INTENT(OUT) :: ec_out(length)
  !! \(\epsilon_c(rho)\) ( NOT \(E_c(\text{rho})\) )
  REAL(DP), INTENT(OUT) :: vc_out(length,sv_d)
  !! \(dE_c(\text{rho})/d\text{rho}  ( NOT d\epsilon_c(\text{rho})/d\text{rho} )
  !
  ! ... local variables
  !
#if defined(__LIBXC)
  TYPE(xc_f03_func_t) :: xc_func
  TYPE(xc_f03_func_info_t) :: xc_info1, xc_info2
  INTEGER :: fkind_x
  REAL(DP) :: amag
  REAL(DP), ALLOCATABLE :: rho_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_lxc(:), vc_lxc(:)
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
#endif
  !
  REAL(DP), ALLOCATABLE :: arho(:), zeta(:)
  !
  INTEGER :: ir, iexch, icorr
  !
  iexch = get_iexch()
  icorr = get_icorr()
  !
  ex_out = 0.0_DP ; vx_out = 0.0_DP
  ec_out = 0.0_DP ; vc_out = 0.0_DP
  !
#if defined(__LIBXC)
  !
  fkind_x = -1
  lengthxc = length
  !
  IF ( ANY(is_libxc(1:2)) ) THEN
    !
    ALLOCATE( rho_lxc(length*sv_d) )
    ALLOCATE( vx_lxc(length*sv_d), vc_lxc(length*sv_d) )
    !
    ! ... set libxc input
    SELECT CASE( sr_d )
    CASE( 1 )
       !
       rho_lxc(:) = ABS(rho_in(:,1))
       !
    CASE( 2 )
       !
       DO ir = 1, length
          rho_lxc(2*ir-1) = (rho_in(ir,1) + rho_in(ir,2)) * 0.5_DP
          rho_lxc(2*ir)   = (rho_in(ir,1) - rho_in(ir,2)) * 0.5_DP
       ENDDO
       !
    CASE( 4 )
       !
       DO ir = 1, length
          amag = SQRT( SUM(rho_in(ir,2:4)**2) )
          rho_lxc(2*ir-1) = (rho_in(ir,1) + amag) * 0.5_DP
          rho_lxc(2*ir)   = (rho_in(ir,1) - amag) * 0.5_DP
       ENDDO
       !
    CASE DEFAULT
       !
       CALL errore( 'xc_LDA', 'Wrong number of spin dimensions', 1 )
       !
    END SELECT
    !
  ENDIF
  !
  !
  ! ... EXCHANGE
  IF ( is_libxc(1) ) THEN
     CALL xc_f03_func_init( xc_func, iexch, sv_d )
       xc_info1 = xc_f03_func_get_info( xc_func )
       CALL xc_f03_func_set_dens_threshold( xc_func, rho_threshold )
       fkind_x  = xc_f03_func_info_get_kind( xc_info1 )
       CALL xc_f03_lda_exc_vxc( xc_func, lengthxc, rho_lxc(1), ex_out(1), vx_lxc(1) )
     CALL xc_f03_func_end( xc_func )
  ENDIF
  !
  ! ... CORRELATION
  IF ( is_libxc(2) ) THEN
     CALL xc_f03_func_init( xc_func, icorr, sv_d )
      xc_info2 = xc_f03_func_get_info( xc_func )
      CALL xc_f03_func_set_dens_threshold( xc_func, rho_threshold )
      CALL xc_f03_lda_exc_vxc( xc_func, lengthxc, rho_lxc(1), ec_out(1), vc_lxc(1) )
     CALL xc_f03_func_end( xc_func )
  ENDIF
  !
  IF ( ((.NOT.is_libxc(1)) .OR. (.NOT.is_libxc(2))) &
        .AND. fkind_x/=XC_EXCHANGE_CORRELATION ) THEN
     !
     SELECT CASE( sr_d )
     CASE( 1 )
        !
        CALL xc_lda( length, ABS(rho_in(:,1)), ex_out, ec_out, vx_out(:,1), vc_out(:,1) )
        !
     CASE( 2 )
        !
        ALLOCATE( arho(length), zeta(length) )
        arho = ABS(rho_in(:,1))
        WHERE (arho > rho_threshold) zeta(:) = rho_in(:,2) / arho(:)
        CALL xc_lsda( length, arho, zeta, ex_out, ec_out, vx_out, vc_out )
        DEALLOCATE( arho, zeta )
        !
     CASE( 4 )
        !
        ALLOCATE( arho(length), zeta(length) )
        arho = ABS( rho_in(:,1) )
        WHERE (arho > rho_threshold) zeta(:) = SQRT( rho_in(:,2)**2 + rho_in(:,3)**2 + &
                                             rho_in(:,4)**2 ) / arho(:) ! amag/arho
        CALL xc_lsda( length, arho, zeta, ex_out, ec_out, vx_out, vc_out )
        DEALLOCATE( arho, zeta )
        !
     CASE DEFAULT
        !
        CALL errore( 'xc_LDA', 'Wrong ns input', 2 )
        !
     END SELECT
     !
  ENDIF
  !
  !  ... fill output arrays
  !  
  IF (sv_d == 1) THEN
     IF (is_libxc(1)) vx_out(:,1) = vx_lxc(:)
     IF (is_libxc(2)) vc_out(:,1) = vc_lxc(:)
  ELSE
     IF (is_libxc(1)) THEN
        DO ir = 1, length
           vx_out(ir,1) = vx_lxc(2*ir-1)
           vx_out(ir,2) = vx_lxc(2*ir)
        ENDDO
     ENDIF
     IF (is_libxc(2)) THEN
        DO ir = 1, length
           vc_out(ir,1) = vc_lxc(2*ir-1)
           vc_out(ir,2) = vc_lxc(2*ir)
        ENDDO
     ENDIF
  ENDIF
  !
  IF (ANY(is_libxc(1:2))) THEN
     DEALLOCATE( rho_lxc )
     DEALLOCATE( vx_lxc, vc_lxc )
  ENDIF
  !
#else
  !
  SELECT CASE( sr_d )
  CASE( 1 )
     !
     IF (use_gpu) THEN
        CALL xc_lda_gpu( length, ABS(rho_in(:,1)), ex_out, ec_out, vx_out(:,1), vc_out(:,1) )
     ELSE
        CALL xc_lda( length, ABS(rho_in(:,1)), ex_out, ec_out, vx_out(:,1), vc_out(:,1) )
     END IF
     !
  CASE( 2 )
     !
     ALLOCATE( arho(length), zeta(length) )
     !
     arho = ABS(rho_in(:,1))
     WHERE (arho > rho_threshold) zeta(:) = rho_in(:,2) / arho(:)
     !
     IF (use_gpu) THEN
        CALL xc_lsda_gpu( length, arho, zeta, ex_out, ec_out, vx_out, vc_out )
     ELSE
        CALL xc_lsda( length, arho, zeta, ex_out, ec_out, vx_out, vc_out )
     END IF
     !
     DEALLOCATE( arho, zeta )
     ! 
  CASE( 4 )
     !
     ALLOCATE( arho(length), zeta(length) )
     !
     arho = ABS( rho_in(:,1) )
     WHERE (arho > rho_threshold) zeta(:) = SQRT( rho_in(:,2)**2 + rho_in(:,3)**2 + &
                                          rho_in(:,4)**2 ) / arho(:) ! amag/arho
     !
     IF (use_gpu) THEN
        CALL xc_lsda_gpu( length, arho, zeta, ex_out, ec_out, vx_out, vc_out )
     ELSE
        CALL xc_lsda( length, arho, zeta, ex_out, ec_out, vx_out, vc_out )
     END IF
     !
     DEALLOCATE( arho, zeta )
     !
  CASE DEFAULT
     !
     CALL errore( 'xc_LDA', 'Wrong ns input', 2 )
     !
  END SELECT
  !
#endif
  !
  !
  RETURN
  !
END SUBROUTINE xc
!
!----------------------------------------------------------------------------
!-------  LDA-LSDA DRIVERS --------------------------------------------------
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
SUBROUTINE xc_lda( length, rho_in, ex_out, ec_out, vx_out, vc_out )
  !--------------------------------------------------------------------------
  !! LDA exchange and correlation functionals - Hartree a.u.
  !
  !! * Exchange:
  !!    * Slater;
  !!    * relativistic Slater.
  !! * Correlation:
  !!    * Ceperley-Alder (Perdew-Zunger parameters);
  !!    * Vosko-Wilk-Nusair;
  !!    * Lee-Yang-Parr;
  !!    * Perdew-Wang;
  !!    * Wigner;
  !!    * Hedin-Lundqvist;
  !!    * Ortiz-Ballone (Perdew-Zunger formula);
  !!    * Ortiz-Ballone (Perdew-Wang formula);
  !!    * Gunnarsson-Lundqvist.
  !
  !! NOTE:
  !! $$ E_x = \int E_x(\text{rho}) dr, E_x(\text{rho}) = 
  !!               \text{rho}\epsilon_c(\text{rho})\ . $$
  !! Same for correlation.
  !
  USE exch_lda
  USE corr_lda
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_in
  !! Charge density
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ex_out
  !! \(\epsilon_x(rho)\) ( NOT \(E_x(\text{rho})\) )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vx_out
  !! \(dE_x(\text{rho})/d\text{rho}\)  ( NOT \(d\epsilon_x(\text{rho})/d\text{rho}\) )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec_out
  !! \(\epsilon_c(rho)\) ( NOT \(E_c(\text{rho})\) )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vc_out
  !! \(dE_c(\text{rho})/d\text{rho}\)  ( NOT \(d\epsilon_c(\text{rho})/d\text{rho}\) )
  !
  ! ... local variables
  !
  INTEGER  :: ir, iexch, icorr
  REAL(DP) :: rho, rs
  REAL(DP) :: ex, ec, ec_
  REAL(DP) :: vx, vc, vc_
  REAL(DP) :: exx_fraction
  REAL(DP) :: finite_size_cell_volume
  LOGICAL :: exx_started, is_there_finite_size_corr
  REAL(DP), PARAMETER :: third = 1.0_DP/3.0_DP, &
                         pi34 = 0.6203504908994_DP, e2 = 2.0_DP
  !                      pi34 = (3/4pi)^(1/3)
  !
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
  !
  iexch = get_iexch()
  icorr = get_icorr()
  exx_started = exx_is_active()
  exx_fraction = get_exx_fraction()
  IF (iexch==8 .OR. icorr==10) THEN
    CALL get_finite_size_cell_volume( is_there_finite_size_corr, &
                                      finite_size_cell_volume )
    !
    IF (.NOT. is_there_finite_size_corr) CALL errore( 'XC',&
        'finite size corrected exchange used w/o initialization', 1 )
  ENDIF
  !
!$omp parallel if(ntids==1)
!$omp do private( rho, rs, ex, ec, ec_, vx, vc, vc_ )
  DO ir = 1, length
     !
     rho = ABS(rho_in(ir))
     !
     ! ... RHO THRESHOLD
     !
     IF ( rho > rho_threshold ) THEN
        rs = pi34 / rho**third
     ELSE
        ex_out(ir) = 0.0_DP  ;  ec_out(ir) = 0.0_DP
        vx_out(ir) = 0.0_DP  ;  vc_out(ir) = 0.0_DP
        CYCLE
     ENDIF
     !
     ! ... EXCHANGE
     !
     SELECT CASE( iexch )
     CASE( 1 )                      ! 'sla'
        !
        CALL slater( rs, ex, vx )
        !
     CASE( 2 )                      ! 'sl1'
        !
        CALL slater1( rs, ex, vx )
        !
     CASE( 3 )                      ! 'rxc'
        !
        CALL slater_rxc( rs, ex, vx )
        !
     CASE( 4, 5 )                   ! 'oep','hf'
        !
        IF ( exx_started ) THEN
           ex = 0.0_DP
           vx = 0.0_DP
        ELSE
           CALL slater( rs, ex, vx )
        ENDIF
        !
     CASE( 6, 7 )                   ! 'pb0x' or 'DF-cx-0', or 'DF2-0',
        !                           ! 'B3LYP'
        CALL slater( rs, ex, vx )
        IF ( exx_started ) THEN
           ex = (1.0_DP - exx_fraction) * ex
           vx = (1.0_DP - exx_fraction) * vx
        ENDIF
        !
     CASE( 8 )                      ! 'sla+kzk'
        !
        CALL slaterKZK( rs, ex, vx, finite_size_cell_volume )
        !
     CASE( 9 )                      ! 'X3LYP'
        !
        CALL slater( rs, ex, vx )
        IF ( exx_started ) THEN
           ex = (1.0_DP - exx_fraction) * ex
           vx = (1.0_DP - exx_fraction) * vx
        ENDIF
        !
     CASE DEFAULT
        !
        ex = 0.0_DP
        vx = 0.0_DP
        !
     END SELECT
     !
     !
     ! ... CORRELATION
     !
     SELECT CASE( icorr )
     CASE( 1 )
        !
        CALL pz( rs, 1, ec, vc )
        !
     CASE( 2 )
        !
        CALL vwn( rs, ec, vc )
        !
     CASE( 3 )
        !
        CALL lyp( rs, ec, vc )
        !
     CASE( 4 )
        !
        CALL pw( rs, 1, ec, vc )
        !
     CASE( 5 )
        !
        CALL wignerc( rs, ec, vc )
        !
     CASE( 6 )
        !
        CALL hl( rs, ec, vc )
        !
     CASE( 7 )
        !
        CALL pz( rs, 2, ec, vc )
        ! 
     CASE( 8 )
        !
        CALL pw( rs, 2, ec, vc )
        !
     CASE( 9 )
        !
        CALL gl( rs, ec, vc )
        !
     CASE( 10 )
        !
        CALL pzKZK( rs, ec, vc, finite_size_cell_volume )
        !
     CASE( 11 )
        !
        CALL vwn1_rpa( rs, ec, vc )
        !
     CASE( 12 )                ! 'B3LYP'
        !
        CALL vwn( rs, ec, vc )
        ec = 0.19_DP * ec
        vc = 0.19_DP * vc
        !
        CALL lyp( rs, ec_, vc_ )
        ec = ec + 0.81_DP * ec_
        vc = vc + 0.81_DP * vc_
        !
     CASE( 13 )                ! 'B3LYP-V1R'
        !
        CALL vwn1_rpa ( rs, ec, vc )
        ec = 0.19_DP * ec
        vc = 0.19_DP * vc
        !
        CALL lyp( rs, ec_, vc_ )
        ec = ec + 0.81_DP * ec_
        vc = vc + 0.81_DP * vc_
        !
     CASE( 14 )                ! 'X3LYP'
        !
        CALL vwn1_rpa( rs, ec, vc )
        ec = 0.129_DP * ec
        vc = 0.129_DP * vc
        !
        CALL lyp( rs, ec_, vc_ )
        ec = ec + 0.871_DP * ec_
        vc = vc + 0.871_DP * vc_
        !
     CASE DEFAULT
        !
        ec = 0.0_DP
        vc = 0.0_DP
        !
     END SELECT
     !
     ex_out(ir) = ex  ;  ec_out(ir) = ec
     vx_out(ir) = vx  ;  vc_out(ir) = vc
     !
  ENDDO
!$omp end do
!$omp end parallel
  !
  !
  RETURN
  !
END SUBROUTINE xc_lda
!
!
!
!----------------------------------------------------------------------------
SUBROUTINE xc_lda_gpu( length, rho_in, ex_out, ec_out, vx_out, vc_out )
  !--------------------------------------------------------------------------
  !! gpu version of xc_lda -- highly experimental --
  !
  USE exch_lda_gpu
  USE corr_lda_gpu
  !
  IMPLICIT NONE
  !
  
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  REAL(DP), INTENT(IN) :: rho_in(length)
  !! Charge density
  REAL(DP), INTENT(OUT) :: ex_out(length)
  !! \(\epsilon_x(rho)\) ( NOT \(E_x(\text{rho})\) )
  REAL(DP), INTENT(OUT) :: vx_out(length)
  !! \(dE_x(\text{rho})/d\text{rho}\)  ( NOT \(d\epsilon_x(\text{rho})/d\text{rho}\) )
  REAL(DP), INTENT(OUT) :: ec_out(length)
  !! \(\epsilon_c(rho)\) ( NOT \(E_c(\text{rho})\) )
  REAL(DP), INTENT(OUT) :: vc_out(length)
  !! \(dE_c(\text{rho})/d\text{rho}\)  ( NOT \(d\epsilon_c(\text{rho})/d\text{rho}\) )
  !
  ! ... local variables
  !
  INTEGER :: ir
  INTEGER :: iexch, icorr
  REAL(DP) :: exx_fraction
  REAL(DP) :: finite_size_cell_volume
  LOGICAL  :: exx_started, is_there_finite_size_corr
  REAL(DP), PARAMETER :: third = 1.0_DP/3.0_DP, &
                         pi34 = 0.6203504908994_DP, e2 = 2.0_DP
                               
  !
  !
  ! .... GPU DEVICE VARIABLES ..........
  !
  REAL(DP), ALLOCATABLE :: rs_d(:)
  REAL(DP), ALLOCATABLE :: ex_d(:)
  REAL(DP), ALLOCATABLE :: vx_d(:)
  REAL(DP), ALLOCATABLE :: ec_d(:), ec_d_(:)
  REAL(DP), ALLOCATABLE :: vc_d(:), vc_d_(:)
#if defined(__CUDA)
  attributes(DEVICE) :: rs_d, ex_d, vx_d, ec_d, ec_d_, vc_d, vc_d_
#endif
  !
  INTEGER :: iflag_ext 
  REAL(DP) :: rho_threshold_d
  !
  !
  ALLOCATE( rs_d(length) )
  ALLOCATE( ex_d(length) )
  ALLOCATE( vx_d(length) )
  ALLOCATE( ec_d(length) )
  ALLOCATE( vc_d(length) )
  !
  iexch = get_iexch()
  icorr = get_icorr()
  !
  IF ( icorr==12 .OR. icorr==13 .OR. icorr==14 ) THEN
    ALLOCATE( ec_d_(length) )
    ALLOCATE( vc_d_(length) )
  ENDIF
  !
  exx_started = exx_is_active()
  exx_fraction = get_exx_fraction()
  IF (iexch==8 .OR. icorr==10) THEN
    CALL get_finite_size_cell_volume( is_there_finite_size_corr, &
                                      finite_size_cell_volume )
    !
    IF (.NOT. is_there_finite_size_corr) CALL errore( 'XC',&
        'finite size corrected exchange used w/o initialization', 1 )
  ENDIF
  !
  rs_d = rho_in
  !
  rho_threshold_d = rho_threshold
  !
!$cuf kernel do(1) <<<*,*>>>
  DO ir = 1, length
     !
     !rho_d = ABS(rho_in_d(ir))
     !
     ! ... RHO THRESHOLD
     !
     IF ( rs_d(ir) > rho_threshold_d ) THEN
        rs_d(ir) = pi34 / rs_d(ir)**third
     ELSE
        ex_d(ir) = 0.0d0  ;  ec_d(ir) = 0.0d0
        vx_d(ir) = 0.0d0  ;  vc_d(ir) = 0.0d0
        CYCLE
     ENDIF
     !
     ! ... EXCHANGE
     !
     SELECT CASE( iexch )
     CASE( 1 )                      ! 'sla'
        !
        CALL slater_d( rs_d(ir), ex_d(ir), vx_d(ir) )
        !
     CASE( 2 )                      ! 'sl1'
        !
        CALL slater1_d( rs_d(ir), ex_d(ir), vx_d(ir) )
        !
     CASE( 3 )                      ! 'rxc'
        !
        CALL slater_rxc_d( rs_d(ir), ex_d(ir), vx_d(ir) )
        !
     CASE( 4, 5 )                   ! 'oep','hf'
        !
        IF ( exx_started ) THEN
           ex_d(ir) = 0.0d0
           vx_d(ir) = 0.0d0
        ELSE
           CALL slater_d( rs_d(ir), ex_d(ir), vx_d(ir) )
        ENDIF
        !
     CASE( 6, 7 )                   ! 'pb0x' or 'DF-cx-0', or 'DF2-0',
        !                           ! 'B3LYP'
        CALL slater_d( rs_d(ir), ex_d(ir), vx_d(ir) )
        IF ( exx_started ) THEN
           ex_d(ir) = (1.0d0 - exx_fraction) * ex_d(ir)
           vx_d(ir) = (1.0d0 - exx_fraction) * vx_d(ir)
        ENDIF
        !
     CASE( 8 )                      ! 'sla+kzk'
        !
        CALL slaterKZK_d( rs_d(ir), ex_d(ir), vx_d(ir), finite_size_cell_volume )
        !
     CASE( 9 )                      ! 'X3LYP'
        !
        CALL slater_d( rs_d(ir), ex_d(ir), vx_d(ir) )
        IF ( exx_started ) THEN
           ex_d(ir) = (1.0d0 - exx_fraction) * ex_d(ir)
           vx_d(ir) = (1.0d0 - exx_fraction) * vx_d(ir)
        ENDIF
        !
     CASE DEFAULT
        !
        ex_d(ir) = 0.0d0
        vx_d(ir) = 0.0d0
        !
     END SELECT
     !
     !
     ! ... CORRELATION
     !
     SELECT CASE( icorr )
     CASE( 1 )
        !
        iflag_ext=1
        CALL pz_d( rs_d(ir), iflag_ext, ec_d(ir), vc_d(ir) )
        !
     CASE( 2 )
        !
        CALL vwn_d( rs_d(ir), ec_d(ir), vc_d(ir) )
        !
     CASE( 3 )
        !
        CALL lyp_d( rs_d(ir), ec_d(ir), vc_d(ir) )
        !
     CASE( 4 )
        !
        iflag_ext=1
        CALL pw_d( rs_d(ir), iflag_ext, ec_d(ir), vc_d(ir) )
        !
     CASE( 5 )
        !
        CALL wignerc_d( rs_d(ir), ec_d(ir), vc_d(ir) )
        !
     CASE( 6 )
        !
        CALL hl_d( rs_d(ir), ec_d(ir), vc_d(ir) )
        !
     CASE( 7 )
        !
        iflag_ext=2
        CALL pz_d( rs_d(ir), iflag_ext, ec_d(ir), vc_d(ir) )
        ! 
     CASE( 8 )
        !
        iflag_ext=2
        CALL pw_d( rs_d(ir), iflag_ext, ec_d(ir), vc_d(ir) )
        !
     CASE( 9 )
        !
        CALL gl_d( rs_d(ir), ec_d(ir), vc_d(ir) )
        !
     CASE( 10 )
        !
        CALL pzKZK_d( rs_d(ir), ec_d(ir), vc_d(ir), finite_size_cell_volume )
        !
     CASE( 11 )
        !
        CALL vwn1_rpa_d( rs_d(ir), ec_d(ir), vc_d(ir) )
        !
     CASE( 12 )                ! 'B3LYP'
        !
        CALL vwn_d( rs_d(ir), ec_d(ir), vc_d(ir) )
        ec_d(ir) = 0.19d0 * ec_d(ir)
        vc_d(ir) = 0.19d0 * vc_d(ir)
        !
        CALL lyp_d( rs_d(ir), ec_d_(ir), vc_d_(ir) )
        ec_d(ir) = ec_d(ir) + 0.81d0 * ec_d_(ir)
        vc_d(ir) = vc_d(ir) + 0.81d0 * vc_d_(ir)
        !
     CASE( 13 )                ! 'B3LYP-V1R'
        !
        CALL vwn1_rpa_d( rs_d(ir), ec_d(ir), vc_d(ir) )
        ec_d(ir) = 0.19d0 * ec_d(ir)
        vc_d(ir) = 0.19d0 * vc_d(ir)
        !
        CALL lyp_d( rs_d(ir), ec_d_(ir), vc_d_(ir) )
        ec_d(ir) = ec_d(ir) + 0.81d0 * ec_d_(ir)
        vc_d(ir) = vc_d(ir) + 0.81d0 * vc_d_(ir)
        !
     CASE( 14 )                ! 'X3LYP'
        !
        CALL vwn1_rpa_d( rs_d(ir), ec_d(ir), vc_d(ir) )
        ec_d(ir) = 0.129d0 * ec_d(ir)
        vc_d(ir) = 0.129d0 * vc_d(ir)
        !
        CALL lyp_d( rs_d(ir), ec_d_(ir), vc_d_(ir) )
        ec_d(ir) = ec_d(ir) + 0.871d0 * ec_d_(ir)
        vc_d(ir) = vc_d(ir) + 0.871d0 * vc_d_(ir)
        !
     CASE DEFAULT
        !
        ec_d(ir) = 0.0d0
        vc_d(ir) = 0.0d0
        !
     END SELECT
     !
  ENDDO
  !
  ex_out = ex_d
  vx_out = vx_d
  ec_out = ec_d
  vc_out = vc_d
  !
  !
  DEALLOCATE( rs_d )
  DEALLOCATE( ex_d )
  DEALLOCATE( vx_d )
  DEALLOCATE( ec_d )
  DEALLOCATE( vc_d )
  IF ( icorr==12 .OR. icorr==13 .OR. icorr==14 ) THEN
    DEALLOCATE( ec_d_ )
    DEALLOCATE( vc_d_ )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE xc_lda_gpu
!

!-----------------------------------------------------------------------------
SUBROUTINE xc_lsda( length, rho_in, zeta_in, ex_out, ec_out, vx_out, vc_out )
  !-----------------------------------------------------------------------------
  !! LSD exchange and correlation functionals - Hartree a.u.
  !
  !! * Exchange:
  !!    * Slater (alpha=2/3).
  !! * Correlation:
  !!    * Ceperley & Alder (Perdew-Zunger parameters);
  !!    * Perdew & Wang.
  !
  USE exch_lda
  USE corr_lda
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_in
  !! Total charge density
  REAL(DP), INTENT(IN),  DIMENSION(length) :: zeta_in
  !! zeta = mag / rho_tot
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ex_out
  !! \(\epsilon_x(rho)\) ( NOT \(E_x(\text{rho})\) )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec_out
  !! \(\epsilon_c(rho)\) ( NOT \(E_c(\text{rho})\) )
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vx_out
  !! \(dE_x(\text{rho})/d\text{rho}\)  ( NOT \(d\epsilon_x(\text{rho})/d\text{rho}\) )
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vc_out
  !! \(dE_c(\text{rho})/d\text{rho}\)  ( NOT \(d\epsilon_c(\text{rho})/d\text{rho}\) )
  !
  ! ...  local variables
  !
  INTEGER  :: ir, iexch, icorr
  REAL(DP) :: rho, rs, zeta
  REAL(DP) :: ex, ec, ec_
  REAL(DP) :: vx(2), vc(2), vc_(2)
  REAL(DP) :: exx_fraction
  LOGICAL :: exx_started
  !
  REAL(DP), PARAMETER :: third = 1.0_DP/3.0_DP, &
                         pi34 = 0.6203504908994_DP
  !                      pi34 = (3/4pi)^(1/3)
  !
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  ntids = omp_get_num_threads()
#endif
  !
  iexch = get_iexch()
  icorr = get_icorr()
  exx_started = exx_is_active()
  exx_fraction = get_exx_fraction()
  !
!$omp parallel if(ntids==1)
!$omp do private( rho, rs, zeta, ex, ec, ec_, vx, vc, vc_ )
  DO ir = 1, length
     !
     zeta = zeta_in(ir)
     IF (ABS(zeta) > 1.D0) zeta = SIGN( 1.D0, zeta )
     !
     rho = ABS(rho_in(ir))
     !
     IF ( rho > rho_threshold ) THEN
        rs = pi34 / rho**third
     ELSE
        ex_out(ir) = 0.0_DP  ;  vx_out(ir,:) = 0.0_DP
        ec_out(ir) = 0.0_DP  ;  vc_out(ir,:) = 0.0_DP
        CYCLE
     ENDIF
     !
     !
     ! ... EXCHANGE
     !
     SELECT CASE( iexch )
     CASE( 1 )                                      ! 'sla'
        !
        CALL slater_spin( rho, zeta, ex, vx(1), vx(2) )
        !
     CASE( 2 )                                      ! 'sl1'
        !
        CALL slater1_spin( rho, zeta, ex, vx(1), vx(2) )
        !
     CASE( 3 )                                      ! 'rxc'
        !
        CALL slater_rxc_spin( rho, zeta, ex, vx(1), vx(2) )
        !
     CASE( 4, 5 )                                   ! 'oep','hf'
        !
        IF ( exx_started ) THEN
           ex = 0.0_DP
           vx = 0.0_DP
        ELSE
           CALL slater_spin( rho, zeta, ex, vx(1), vx(2) )
        ENDIF
        !
     CASE( 6 )                                      ! 'pb0x'
        !
        CALL slater_spin( rho, zeta, ex, vx(1), vx(2) )
        IF ( exx_started ) THEN
           ex = (1.0_DP - exx_fraction) * ex
           vx = (1.0_DP - exx_fraction) * vx
        ENDIF
        !
     CASE( 7 )                                      ! 'B3LYP'
        !
        CALL slater_spin( rho, zeta, ex, vx(1), vx(2) )
        IF ( exx_started ) THEN
           ex = (1.0_DP - exx_fraction) * ex
           vx = (1.0_DP - exx_fraction) * vx
        ENDIF
        !
     CASE( 9 )                                      ! 'X3LYP'
        !
        CALL slater_spin( rho, zeta, ex, vx(1), vx(2) )
        IF ( exx_started ) THEN
           ex = (1.0_DP - exx_fraction) * ex
           vx = (1.0_DP - exx_fraction) * vx
        ENDIF
        !
     CASE DEFAULT
        !
        ex = 0.0_DP
        vx = 0.0_DP
        !
     END SELECT
     !
     !
     ! ... CORRELATION
     !
     SELECT CASE( icorr )
     CASE( 0 )
        !
        ec = 0.0_DP
        vc = 0.0_DP
        !
     CASE( 1 )
        !
        CALL pz_spin( rs, zeta, ec, vc(1), vc(2) )
        !
     CASE( 2 )
        !
        CALL vwn_spin( rs, zeta, ec, vc(1), vc(2) )
        !
     CASE( 3 )
        !
        CALL lsd_lyp( rho, zeta, ec, vc(1), vc(2) )       ! from CP/FPMD (more_functionals)
        !
     CASE( 4 )
        !
        CALL pw_spin( rs, zeta, ec, vc(1), vc(2) )
        !
     CASE( 12 )                                           ! 'B3LYP'
        !
        CALL vwn_spin( rs, zeta, ec, vc(1), vc(2) )
        ec = 0.19_DP * ec
        vc = 0.19_DP * vc
        !
        CALL lsd_lyp( rho, zeta, ec_, vc_(1), vc_(2) )    ! from CP/FPMD (more_functionals)
        ec = ec + 0.81_DP * ec_
        vc = vc + 0.81_DP * vc_
        !     
     CASE( 13 )                                           ! 'B3LYP-V1R'
        !
        CALL vwn1_rpa_spin( rs, zeta, ec, vc(1), vc(2) )
        ec = 0.19_DP * ec
        vc = 0.19_DP * vc
        !
        CALL lsd_lyp( rho, zeta, ec_, vc_(1), vc_(2) )    ! from CP/FPMD (more_functionals)
        ec = ec + 0.81_DP * ec_
        vc = vc + 0.81_DP * vc_
        !
     CASE( 14 )                                           ! 'X3LYP
        !
        CALL vwn1_rpa_spin( rs, zeta, ec, vc(1), vc(2) )
        ec = 0.129_DP * ec
        vc = 0.129_DP * vc
        !
        CALL lsd_lyp( rho, zeta, ec_, vc_(1), vc_(2) )    ! from CP/FPMD (more_functionals)
        ec = ec + 0.871_DP * ec_
        vc = vc + 0.871_DP * vc_
        !
     CASE DEFAULT
        !
        ec = 0.0_DP
        vc = 0.0_DP
        !
     END SELECT
     !
     ex_out(ir) = ex  ;  vx_out(ir,:) = vx(:)
     ec_out(ir) = ec  ;  vc_out(ir,:) = vc(:)
     !
  ENDDO
!$omp end do
!$omp end parallel
  !
  !
  RETURN
  !
END SUBROUTINE xc_lsda
!
!
!-----------------------------------------------------------------------------
SUBROUTINE xc_lsda_gpu( length, rho_in, zeta_in, ex_out, ec_out, vx_out, vc_out )
  !-----------------------------------------------------------------------------
  !! gpu version of xc_lsda -- highly experimental --
  !
  USE exch_lda_gpu
  USE corr_lda_gpu
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_in
  !! Total charge density
  REAL(DP), INTENT(IN),  DIMENSION(length) :: zeta_in
  !! zeta = mag / rho_tot
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ex_out
  !! \(\epsilon_x(rho)\) ( NOT \(E_x(\text{rho})\) )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: ec_out
  !! \(\epsilon_c(rho)\) ( NOT \(E_c(\text{rho})\) )
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vx_out
  !! \(dE_x(\text{rho})/d\text{rho}\)  ( NOT \(d\epsilon_x(\text{rho})/d\text{rho}\) )
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vc_out
  !! \(dE_c(\text{rho})/d\text{rho}\)  ( NOT \(d\epsilon_c(\text{rho})/d\text{rho}\) )
  !
  ! ...  local variables
  !
  INTEGER :: ir
  !
  INTEGER :: iexch, icorr
  REAL(8) :: exx_fraction
  LOGICAL :: exx_started
  REAL(8), PARAMETER :: third = 1.0d0/3.0d0, &
                        pi34 = 0.6203504908994d0
  !                     pi34 = (3/4pi)^(1/3)
  !
  !
  ! .... GPU DEVICE VARIABLES ..........
  !
  REAL(8), ALLOCATABLE :: rs_d(:)
  REAL(8), ALLOCATABLE :: rho_d(:)
  REAL(8), ALLOCATABLE :: zeta_d(:)
  REAL(8), ALLOCATABLE :: ex_d(:)
  REAL(8), ALLOCATABLE :: vx_d(:,:)
  REAL(8), ALLOCATABLE :: ec_d(:),   ec_d_(:)
  REAL(8), ALLOCATABLE :: vc_d(:,:), vc_d_(:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: rs_d, rho_d, zeta_d, ex_d, vx_d, ec_d, ec_d_, vc_d, vc_d_
#endif
  !
  REAL(8) :: rho_threshold_d
#if defined(__CUDA)
  attributes(DEVICE) :: rho_threshold_d
#endif
  !
  !
  ALLOCATE( rho_d(length) )
  ALLOCATE( rs_d(length)  )
  ALLOCATE( zeta_d(length) )
  ALLOCATE( ex_d(length) )
  ALLOCATE( vx_d(length,2) )
  ALLOCATE( ec_d(length) )
  ALLOCATE( vc_d(length,2) )
  !
  iexch = get_iexch()
  icorr = get_icorr()
  !
  IF ( icorr==12 .OR. icorr==13 .OR. icorr==14 ) THEN
    ALLOCATE( ec_d_(length) )
    ALLOCATE( vc_d_(length,2) )
  ENDIF
  !
  exx_started = exx_is_active()
  exx_fraction = get_exx_fraction()
  !
  rho_threshold_d = rho_threshold
  !
  rho_d = rho_in
  rs_d  = rho_in
  zeta_d = zeta_in
  !
  !
!$cuf kernel do(1) <<<*,*>>>
  DO ir = 1, length
     !
     IF (ABS(zeta_d(ir)) > 1.d0) zeta_d(ir) = SIGN( 1.d0, zeta_d(ir) )
     !
     IF ( rho_d(ir) > rho_threshold_d ) THEN
        rs_d(ir) = pi34 / rho_d(ir)**third
     ELSE
        ex_d(ir)=0.0d0 ; vx_d(ir,:)=0.0d0 !; vx_d(ir,2)=0.0d0
        ec_d(ir)=0.0d0 ; vc_d(ir,1)=0.0d0 ; vc_d(ir,2)=0.0d0
        CYCLE
     ENDIF
     !
     !
     ! ... EXCHANGE
     !
     SELECT CASE( iexch )
     CASE( 1 )                                      ! 'sla'
        !
        CALL slater_spin_d( rho_d(ir), zeta_d(ir), ex_d(ir), vx_d(ir,1), vx_d(ir,2)  )
        !
     CASE( 2 )                                      ! 'sl1'
        !
        CALL slater1_spin_d( rho_d(ir), zeta_d(ir), ex_d(ir), vx_d(ir,1), vx_d(ir,2) )
        !
     CASE( 3 )                                      ! 'rxc'
        !
        CALL slater_rxc_spin_d( rho_d(ir), zeta_d(ir), ex_d(ir), vx_d(ir,1), vx_d(ir,2) )
        !
     CASE( 4, 5 )                                   ! 'oep','hf'
        !
        IF ( exx_started ) THEN
           ex_d(ir) = 0.0d0
           vx_d(ir,1) = 0.0d0
           vx_d(ir,2) = 0.0d0
        ELSE
           CALL slater_spin_d( rho_d(ir), zeta_d(ir), ex_d(ir), vx_d(ir,1), vx_d(ir,2) )
        ENDIF
        !
     CASE( 6 )                                      ! 'pb0x'
        !
        CALL slater_spin_d( rho_d(ir), zeta_d(ir), ex_d(ir), vx_d(ir,1), vx_d(ir,2) )
        IF ( exx_started ) THEN
           ex_d(ir) = (1.0d0 - exx_fraction) * ex_d(ir)
           vx_d(ir,1) = (1.0d0 - exx_fraction) * vx_d(ir,1)
           vx_d(ir,2) = (1.0d0 - exx_fraction) * vx_d(ir,2)
        ENDIF
        !
     CASE( 7 )                                      ! 'B3LYP'
        !
        CALL slater_spin_d( rho_d(ir), zeta_d(ir), ex_d(ir), vx_d(ir,1), vx_d(ir,2) )
        IF ( exx_started ) THEN
           ex_d(ir) = (1.0d0 - exx_fraction) * ex_d(ir)
           vx_d(ir,1) = (1.0d0 - exx_fraction) * vx_d(ir,1)
           vx_d(ir,2) = (1.0d0 - exx_fraction) * vx_d(ir,2)
        ENDIF
        !
     CASE( 9 )                                      ! 'X3LYP'
        !
        CALL slater_spin_d( rho_d(ir), zeta_d(ir), ex_d(ir), vx_d(ir,1), vx_d(ir,2) )
        IF ( exx_started ) THEN
           ex_d(ir) = (1.0d0 - exx_fraction) * ex_d(ir)
           vx_d(ir,1) = (1.0d0 - exx_fraction) * vx_d(ir,1)
           vx_d(ir,2) = (1.0d0 - exx_fraction) * vx_d(ir,2)
        ENDIF
        !
     CASE DEFAULT
        !
        ex_d(ir) = 0.0d0
        vx_d(ir,1) = 0.0d0
        vx_d(ir,2) = 0.0d0
        !
     END SELECT
     !
     !
     ! ... CORRELATION
     !
     SELECT CASE( icorr )
     CASE( 0 )
        !
        ec_d(ir) = 0.0d0
        vc_d(ir,1) = 0.0d0
        vc_d(ir,2) = 0.0d0
        !
     CASE( 1 )
        !
        CALL pz_spin_d( rs_d(ir), zeta_d(ir), ec_d(ir), vc_d(ir,1), vc_d(ir,2) )
        !
     CASE( 2 )
        !
        CALL vwn_spin_d( rs_d(ir), zeta_d(ir), ec_d(ir), vc_d(ir,1), vc_d(ir,2) )
        !
     CASE( 3 )
        !
        CALL lsd_lyp_d( rho_d(ir), zeta_d(ir), ec_d(ir), vc_d(ir,1), vc_d(ir,2) )     ! from CP/FPMD (more_functionals)
        !
     CASE( 4 )
        !
        CALL pw_spin_d( rs_d(ir), zeta_d(ir), ec_d(ir), vc_d(ir,1), vc_d(ir,2) )
        !
     CASE( 12 )                                           ! 'B3LYP'
        !
        CALL vwn_spin_d( rs_d(ir), zeta_d(ir), ec_d(ir), vc_d(ir,1), vc_d(ir,2) )
        ec_d(ir) = 0.19d0 * ec_d(ir)
        vc_d(ir,1) = 0.19d0 * vc_d(ir,1)
        vc_d(ir,2) = 0.19d0 * vc_d(ir,2)
        !
        CALL lsd_lyp_d( rho_d(ir), zeta_d(ir), ec_d_(ir), vc_d_(ir,1), vc_d_(ir,2) )   ! from CP/FPMD (more_functionals)
        ec_d(ir) = ec_d(ir) + 0.81d0 * ec_d_(ir)
        vc_d(ir,1) = vc_d(ir,1) + 0.81d0 * vc_d_(ir,1)
        vc_d(ir,2) = vc_d(ir,2) + 0.81d0 * vc_d_(ir,2)
        !     
     CASE( 13 )                                           ! 'B3LYP-V1R'
        !
        CALL vwn1_rpa_spin_d( rs_d(ir), zeta_d(ir), ec_d(ir), vc_d(ir,1), vc_d(ir,2) )
        ec_d(ir) = 0.19d0 * ec_d(ir)
        vc_d(ir,1) = 0.19d0 * vc_d(ir,1)
        vc_d(ir,2) = 0.19d0 * vc_d(ir,2)
        !
        CALL lsd_lyp_d( rho_d(ir), zeta_d(ir), ec_d_(ir), vc_d_(ir,1), vc_d_(ir,2) )  ! from CP/FPMD (more_functionals)
        ec_d(ir) = ec_d(ir) + 0.81d0 * ec_d_(ir)
        vc_d(ir,1) = vc_d(ir,1) + 0.81d0 * vc_d_(ir,1)
        vc_d(ir,2) = vc_d(ir,2) + 0.81d0 * vc_d_(ir,2)
        !
     CASE( 14 )                                           ! 'X3LYP
        !
        CALL vwn1_rpa_spin_d( rs_d(ir), zeta_d(ir), ec_d(ir), vc_d(ir,1), vc_d(ir,2) )
        ec_d(ir) = 0.129d0 * ec_d(ir)
        vc_d(ir,1) = 0.129d0 * vc_d(ir,1)
        vc_d(ir,2) = 0.129d0 * vc_d(ir,2)
        !
        CALL lsd_lyp_d( rho_d(ir), zeta_d(ir), ec_d_(ir), vc_d_(ir,1), vc_d_(ir,2) )   ! from CP/FPMD (more_functionals)
        ec_d(ir) = ec_d(ir) + 0.871d0 * ec_d_(ir)
        vc_d(ir,1) = vc_d(ir,1) + 0.871d0 * vc_d_(ir,1)
        vc_d(ir,2) = vc_d(ir,2) + 0.871d0 * vc_d_(ir,2)
        !
     CASE DEFAULT
        !
        !CALL errore( 'xc_lda_lsda_drivers (xc_lsda)', 'not implemented', icorr )
        !
     END SELECT
     !
     !
  ENDDO
  !
  ex_out = ex_d  ;  vx_out = vx_d
  ec_out = ec_d  ;  vc_out = vc_d
  !
  DEALLOCATE( rho_d )
  DEALLOCATE( rs_d  )
  DEALLOCATE( zeta_d )
  DEALLOCATE( ex_d )
  DEALLOCATE( vx_d )
  DEALLOCATE( ec_d )
  DEALLOCATE( vc_d )
  IF ( icorr==12 .OR. icorr==13 .OR. icorr==14 ) THEN
     DEALLOCATE( ec_d_ )
     DEALLOCATE( vc_d_ )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE xc_lsda_gpu
!
!
END MODULE xc_lda_lsda
