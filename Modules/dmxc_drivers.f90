!
! Copyright (C) 2004-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
SUBROUTINE dmxc( length, sr_d, rho_in, dmuxc )
  !---------------------------------------------------------------------
  !! Wrapper routine. Calls dmxc-driver routines from internal libraries
  !! or from the external one 'libxc', depending on the input choice.
  !
  ! Only two possibilities in the present version (LDA only):
  ! 1) iexch libxc + icorr libxc
  ! 2) iexch qe    + icorr qe
  !
  USE kinds,            ONLY: DP
  USE funct,            ONLY: get_iexch, get_icorr, is_libxc
  USE xc_lda_lsda,      ONLY: xc_lda, xc_lsda
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
#endif
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: sr_d
  !! number of spin components
  REAL(DP), INTENT(IN) :: rho_in(length,sr_d)
  !! charge density
  REAL(DP), INTENT(OUT) :: dmuxc(length,sr_d,sr_d)
  !! the derivative of the xc potential
  !
  ! ... local variables
  !
#if defined(__LIBXC)
  TYPE(xc_f03_func_t) :: xc_func
  TYPE(xc_f03_func_info_t) :: xc_info1, xc_info2
  INTEGER :: pol_unpol
  REAL(DP), ALLOCATABLE :: rho_lxc(:)
  REAL(DP), ALLOCATABLE :: dmxc_lxc(:), dmex_lxc(:), dmcr_lxc(:)
  LOGICAL :: exch_lxc_avail, corr_lxc_avail
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
#endif
  !
  INTEGER :: iexch, icorr
  INTEGER :: ir, length_lxc, length_dlxc
  REAL(DP), PARAMETER :: small = 1.E-10_DP, rho_trash = 0.5_DP
  !
  iexch = get_iexch()
  icorr = get_icorr()
  !
#if defined(__LIBXC)
  !
  lengthxc = length
  !
  IF ( (is_libxc(1) .OR. iexch==0) .AND. (is_libxc(2) .OR. icorr==0)) THEN
    !
    length_lxc = length*sr_d
    !
    ! ... set libxc input
    SELECT CASE( sr_d )
    CASE( 1 )
      !
      ALLOCATE( rho_lxc(length_lxc) )
      pol_unpol = 1
      rho_lxc = rho_in(:,1) 
      !
    CASE( 2 )
      !
      ALLOCATE( rho_lxc(length_lxc) )
      pol_unpol = 2
      DO ir = 1, length
        rho_lxc(2*ir-1) = rho_in(ir,1)
        rho_lxc(2*ir)   = rho_in(ir,2)
      ENDDO
      !
    CASE( 4 )
      !
      CALL errore( 'dmxc', 'The derivative of the xc potential with libxc &
                           &is not available for noncollinear case', 1 )
      !
    CASE DEFAULT
      !
      CALL errore( 'dmxc', 'Wrong number of spin dimensions', 2 )
      !
    END SELECT
    !
    length_dlxc = length
    IF (pol_unpol == 2) length_dlxc = length*3
    !
    !
    ALLOCATE( dmex_lxc(length_dlxc), dmcr_lxc(length_dlxc), &
              dmxc_lxc(length_dlxc) )
    !
    ! ... DERIVATIVE FOR EXCHANGE
    dmex_lxc(:) = 0.0_DP
    IF (iexch /= 0) THEN    
       CALL xc_f03_func_init( xc_func, iexch, pol_unpol )
        xc_info1 = xc_f03_func_get_info( xc_func )
        CALL xc_f03_lda_fxc( xc_func, lengthxc, rho_lxc(1), dmex_lxc(1) )
       CALL xc_f03_func_end( xc_func )
    ENDIF    
    !
    ! ... DERIVATIVE FOR CORRELATION
    dmcr_lxc(:) = 0.0_DP
    IF (icorr /= 0) THEN
       CALL xc_f03_func_init( xc_func, icorr, pol_unpol )
        xc_info2 = xc_f03_func_get_info( xc_func )
        CALL xc_f03_lda_fxc( xc_func, lengthxc, rho_lxc(1), dmcr_lxc(1) )
       CALL xc_f03_func_end( xc_func )
    ENDIF
    !
    dmxc_lxc = (dmex_lxc + dmcr_lxc)*2.0_DP
    !
    IF (sr_d == 1) THEN
      dmuxc(:,1,1) = dmxc_lxc(:)
    ELSEIF (sr_d == 2) THEN
      DO ir = 1, length
        dmuxc(ir,1,1) = dmxc_lxc(3*ir-2)
        dmuxc(ir,1,2) = dmxc_lxc(3*ir-1)
        dmuxc(ir,2,1) = dmxc_lxc(3*ir-1)
        dmuxc(ir,2,2) = dmxc_lxc(3*ir)
      ENDDO
    ENDIF
    !
    DEALLOCATE( dmex_lxc, dmcr_lxc, dmxc_lxc )
    DEALLOCATE( rho_lxc )
    !
  ELSEIF ((.NOT.is_libxc(1)) .AND. (.NOT.is_libxc(2)) ) THEN
    !
    IF ( sr_d == 1 ) CALL dmxc_lda( length, rho_in(:,1), dmuxc(:,1,1) )
    IF ( sr_d == 2 ) CALL dmxc_lsda( length, rho_in, dmuxc )
    !
  ELSE
    !
    CALL errore( 'dmxc', 'Derivatives of exchange and correlation terms, &
                        & at present, must be both qe or both libxc.', 3 )
    !
  ENDIF
  !
#else
  !
  SELECT CASE( sr_d )
  CASE( 1 )
     !
     CALL dmxc_lda( length, rho_in(:,1), dmuxc(:,1,1) )
     !
  CASE( 2 )
     !
     CALL dmxc_lsda( length, rho_in, dmuxc )
     ! 
  CASE( 4 )
     !
     CALL dmxc_nc( length, rho_in(:,1), rho_in(:,2:4), dmuxc )
     !
  CASE DEFAULT
     !
     CALL errore( 'xc_LDA', 'Wrong ns input', 4 )
     !
  END SELECT
  !
#endif
  !
  !
  RETURN
  !
END SUBROUTINE
!
!
!-----------------------------------------------------------------------
SUBROUTINE dmxc_lda( length, rho_in, dmuxc )
  !---------------------------------------------------------------------
  !! Computes the derivative of the xc potential with respect to the 
  !! local density.
  !
  USE xc_lda_lsda,  ONLY: xc_lda
  USE exch_lda,     ONLY: slater
  USE funct,        ONLY: get_iexch, get_icorr
  USE kinds,        ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the input/output arrays
  REAL(DP), INTENT(IN),  DIMENSION(length) :: rho_in
  !! the charge density ( positive )
  REAL(DP), INTENT(OUT), DIMENSION(length) :: dmuxc
  !! the derivative of the xc potential
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: ex, vx
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: arho, rhoaux, dr
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: ec, vc
  !
  REAL(DP) :: rho, rs, ex_s, vx_s
  REAL(DP) :: dpz
  INTEGER  :: iexch, icorr
  INTEGER  :: iflg, ir, i1, i2, f1, f2
  !
  REAL(DP), PARAMETER :: small = 1.E-30_DP, e2 = 2.0_DP,        &
                         pi34 = 0.75_DP/3.141592653589793_DP,   &
                         third = 1.0_DP/3.0_DP, rho_trash = 0.5_DP, &
                         rs_trash = 1.0_DP
#if defined(_OPENMP)
  INTEGER :: ntids
  INTEGER, EXTERNAL :: omp_get_num_threads
  !
  !
  ntids = omp_get_num_threads()
#endif  
  !
  iexch = get_iexch()
  icorr = get_icorr()
  !
  dmuxc = 0.0_DP
  !
  ! ... first case: analytical derivatives available
  !
  IF (iexch == 1 .AND. icorr == 1) THEN
  !
!$omp parallel if(ntids==1)
!$omp do private( rs, rho, ex_s, vx_s , iflg)
     DO ir = 1, length
        !
        rho = rho_in(ir)
        IF ( rho < -small ) rho = -rho_in(ir)
        !
        IF ( rho > small ) THEN
           rs = (pi34 / rho)**third
        ELSE
           dmuxc(ir) = 0.0_DP
           CYCLE
        ENDIF
        !
        CALL slater( rs, ex_s, vx_s )
        dmuxc(ir) = vx_s / (3.0_DP * rho)
        !
        iflg = 2
        IF (rs < 1.0_DP) iflg = 1
        dmuxc(ir) = dmuxc(ir) + dpz( rs, iflg )
        dmuxc(ir) = dmuxc(ir) * SIGN(1.0_DP,rho_in(ir))
        !
     ENDDO
!$omp end do
!$omp end parallel
     !
  ELSE
     !
     ! ... second case: numerical derivatives
     !
     ALLOCATE( ex(2*length), vx(2*length)  )
     ALLOCATE( ec(2*length), vc(2*length)  )
     ALLOCATE( arho(length), dr(length), rhoaux(2*length) )
     !
     i1 = 1         ;  f1 = length             !two blocks:  [ rho+dr ]
     i2 = length+1  ;  f2 = 2*length           !             [ rho-dr ]              
     !
     arho = ABS(rho_in)
     dr = 0.0_DP
     WHERE ( arho > small ) dr = MIN( 1.E-6_DP, 1.E-4_DP * rho_in )
     !
     rhoaux(i1:f1) = arho+dr
     rhoaux(i2:f2) = arho-dr
     !
     CALL xc_lda( length*2, rhoaux, ex, ec, vx, vc )
     !
     WHERE ( arho < small ) dr = 1.0_DP ! ... to avoid NaN in the next operation
     !
     dmuxc(:) = (vx(i1:f1) + vc(i1:f1) - vx(i2:f2) - vc(i2:f2)) / &
                (2.0_DP * dr(:))
     !
     DEALLOCATE( ex, vx  )
     DEALLOCATE( ec, vc  )
     DEALLOCATE( dr, rhoaux )
     !
     WHERE ( arho < small ) dmuxc = 0.0_DP
     ! however a higher threshold is already present in xc_lda()
     dmuxc(:) = dmuxc(:) * SIGN(1.0_DP,rho_in(:))
     !
     DEALLOCATE( arho )
     !
  ENDIF
  !
  ! bring to rydberg units
  !
  dmuxc = e2 * dmuxc
  !
  RETURN
  !
END SUBROUTINE dmxc_lda
!
!
!-----------------------------------------------------------------------
SUBROUTINE dmxc_lsda( length, rho_in, dmuxc )
!-----------------------------------------------------------------------
  !! Computes the derivative of the xc potential with respect to the 
  !! local density in the spin-polarized case.
  !
  USE kinds,          ONLY: DP
  USE funct,          ONLY: get_iexch, get_icorr
  USE xc_lda_lsda,    ONLY: xc_lsda
  USE exch_lda,       ONLY: slater
  USE corr_lda,       ONLY: pz, pz_polarized
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the input/output arrays
  REAL(DP), INTENT(IN), DIMENSION(length,2) :: rho_in
  !! spin-up and spin-down charge density
  REAL(DP), INTENT(OUT), DIMENSION(length,2,2) :: dmuxc
  !! u-u, u-d, d-u, d-d derivatives of the XC functional
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: rhotot(:), zeta(:), zeta_eff(:)
  !
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: aux1, aux2, dr, dz
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: rhoaux, zetaux
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: vx, vc, vxc
  REAL(DP) :: ecu, ecp, ex_s
  REAL(DP) :: vcu, vcp, vx_s
  !
  REAL(DP) :: fz, fz1, fz2, dmcu, dmcp, aa, bb, cc
  REAL(DP) :: rs, zeta_s
  !
  REAL(DP) :: dpz, dpz_polarized
  !
  INTEGER :: iexch, icorr
  INTEGER :: ir, is, iflg
  INTEGER :: i1, i2, i3, i4
  INTEGER :: f1, f2, f3, f4
  !
  REAL(DP), PARAMETER :: small = 1.E-30_DP, e2 = 2.0_DP,      &
                         pi34 = 0.75_DP/3.141592653589793_DP, &
                         third = 1.0_DP/3.0_DP, p43 = 4.0_DP/3.0_DP, &
                         p49 = 4.0_DP/9.0_DP, m23 = -2.0_DP/3.0_DP
  !
  iexch = get_iexch()
  icorr = get_icorr()
  !
  dmuxc = 0.0_DP
  ALLOCATE(rhotot(length)) 
  rhotot(:) = rho_in(:,1) + rho_in(:,2)
  !
  IF (iexch == 1 .AND. icorr == 1) THEN
     !
     ! ... first case: analytical derivative available
     !
     !$omp parallel do default(private) shared(length,rhotot, rho_in, dmuxc )   
     DO ir = 1, length
        !
        IF (rhotot(ir) < small) CYCLE
        zeta_s = (rho_in(ir,1) - rho_in(ir,2)) / rhotot(ir)
        IF (ABS(zeta_s) > 1.0_DP) CYCLE
        !
        ! ... exchange
        !
        rs = ( pi34 / (2.0_DP * rho_in(ir,1)) )**third
        CALL slater( rs, ex_s, vx_s )
        !
        dmuxc(ir,1,1) = vx_s / (3.0_DP * rho_in(ir,1))
        !
        rs = ( pi34 / (2.0_DP * rho_in(ir,2)) )**third
        CALL slater( rs, ex_s, vx_s )
        !
        dmuxc(ir,2,2) = vx_s / (3.0_DP * rho_in(ir,2))
        !
        ! ... correlation
        !
        rs = (pi34 / rhotot(ir))**third
        !
        CALL pz( rs, 1, ecu, vcu )
        CALL pz_polarized( rs, ecp, vcp )
        !
        fz  = ( (1.0_DP + zeta_s)**p43 + (1.0_DP - zeta_s)**p43 - 2.0_DP ) &
                  / (2.0_DP**p43 - 2.0_DP)
        fz1 = p43 * ( (1.0_DP + zeta_s)**third - (1.0_DP - zeta_s)**third) &
                  / (2.0_DP**p43 - 2.0_DP)
        fz2 = p49 * ( (1.0_DP + zeta_s)**m23   + (1.0_DP - zeta_s)**m23)   &
                  / (2.0_DP**p43 - 2.0_DP)
        !
        iflg = 2
        IF (rs < 1.0_DP) iflg = 1
        !
        dmcu = dpz( rs, iflg )
        dmcp = dpz_polarized( rs, iflg )
        !
        aa = dmcu + fz * (dmcp - dmcu)
        bb = 2.0_DP * fz1 * (vcp - vcu - (ecp - ecu) ) / rhotot(ir)
        cc = fz2 * (ecp - ecu) / rhotot(ir)
        !
        dmuxc(ir,1,1) = dmuxc(ir,1,1) + aa + (1.0_DP - zeta_s) * bb +  &
                                             (1.0_DP - zeta_s)**2 * cc
        dmuxc(ir,2,1) = dmuxc(ir,2,1) + aa + (-zeta_s) * bb +          &
                                             (zeta_s**2 - 1.0_DP) * cc
        dmuxc(ir,1,2) = dmuxc(ir,2,1)
        dmuxc(ir,2,2) = dmuxc(ir,2,2) + aa - (1.0_DP + zeta_s) * bb +  &
                                             (1.0_DP + zeta_s)**2 * cc
     ENDDO
     !
  ELSE
     !
     !
     ALLOCATE( vx(4*length,2) , vc(4*length,2), vxc(2*length,2) )
     ALLOCATE( rhoaux(4*length), zetaux(4*length) )
     ALLOCATE( aux1(4*length) , aux2(4*length) )
     ALLOCATE( dr(length), dz(length) )
     ALLOCATE( zeta(length), zeta_eff(length)) 
     !
     i1 = 1     ;   f1 = length          !  four blocks:  [ rho+dr , zeta    ]
     i2 = f1+1  ;   f2 = 2*length        !                [ rho-dr , zeta    ]
     i3 = f2+1  ;   f3 = 3*length        !                [ rho    , zeta+dz ]
     i4 = f3+1  ;   f4 = 4*length        !                [ rho    , zeta-dz ]
     !
     !
     dz(:) = 1.E-6_DP  ! dz(:) = MIN( 1.d-6, 1.d-4*ABS(zeta(:)) )
     !
     ! ... THRESHOLD STUFF AND dr(:)
     dr(:) = 0.0_DP
     zeta(:) = 0.0_dp
     zeta_eff(:) = 0.0_dp
     DO ir = 1, length
        IF (rhotot(ir) > small) THEN
           zeta_s = (rho_in(ir,1) - rho_in(ir,2)) / rhotot(ir)
           zeta(ir) = zeta_s
           ! ... If zeta is too close to +-1, the derivative is computed at a slightly
           ! smaller zeta
           zeta_eff(ir) = SIGN( MIN( ABS(zeta_s), (1.0_DP-2.0_DP*dz(ir)) ), zeta_s )
           dr(ir) = MIN( 1.E-6_DP, 1.E-4_DP * rhotot(ir) )
           IF (ABS(zeta_s) > 1.0_DP) THEN  
             rhotot(ir) = 0.d0 ;  dr(ir) = 0.d0 ! e.g. vx=vc=0.0
           ENDIF
        ENDIF
     ENDDO
     !
     rhoaux(i1:f1) = rhotot + dr    ;   zetaux(i1:f1) = zeta
     rhoaux(i2:f2) = rhotot - dr    ;   zetaux(i2:f2) = zeta
     rhoaux(i3:f3) = rhotot         ;   zetaux(i3:f3) = zeta_eff + dz
     rhoaux(i4:f4) = rhotot         ;   zetaux(i4:f4) = zeta_eff - dz
     !
     CALL xc_lsda( length*4, rhoaux, zetaux, aux1, aux2, vx, vc )
     !
     WHERE (rhotot <= small)  ! ... to avoid NaN in the next operations
        dr=1.0_DP ; rhotot=0.5d0
     END WHERE
     !
     dmuxc(:,1,1) = ( vx(i1:f1,1) + vc(i1:f1,1) - vx(i2:f2,1) - vc(i2:f2,1) ) / (2.0_DP*dr)
     dmuxc(:,2,2) = ( vx(i1:f1,2) + vc(i1:f1,2) - vx(i2:f2,2) - vc(i2:f2,2) ) / (2.0_DP*dr)
     !
     aux1(i1:f1) = 1.0_DP / rhotot(:) / (2.0_DP*dz(:))
     aux1(i2:f2) = aux1(i1:f1)
     !
     vxc(i1:f2,1) = ( vx(i3:f4,1) + vc(i3:f4,1) ) * aux1(i1:f2)
     vxc(i1:f2,2) = ( vx(i3:f4,2) + vc(i3:f4,2) ) * aux1(i1:f2)
     !
     dmuxc(:,2,1) = dmuxc(:,1,1) - (vxc(i1:f1,1) - vxc(i2:f2,1)) * (1.0_DP+zeta)
     dmuxc(:,1,2) = dmuxc(:,2,2) + (vxc(i1:f1,2) - vxc(i2:f2,2)) * (1.0_DP-zeta)
     dmuxc(:,1,1) = dmuxc(:,1,1) + (vxc(i1:f1,1) - vxc(i2:f2,1)) * (1.0_DP-zeta)
     dmuxc(:,2,2) = dmuxc(:,2,2) - (vxc(i1:f1,2) - vxc(i2:f2,2)) * (1.0_DP+zeta)
     !
     DEALLOCATE( vx, vc, vxc )
     DEALLOCATE( rhoaux, zetaux )
     DEALLOCATE( aux1, aux2 )
     DEALLOCATE( dr, dz )
     !
  ENDIF
  !
  ! ... bring to Rydberg units
  !
  dmuxc = e2 * dmuxc
  !
  RETURN
  !
END SUBROUTINE dmxc_lsda
!
!
!-----------------------------------------------------------------------
SUBROUTINE dmxc_nc( length, rho_in, m, dmuxc )
!-----------------------------------------------------------------------
  !! Computes the derivative of the xc potential with respect to the 
  !! local density and magnetization in the non-collinear case.
  !
  USE xc_lda_lsda,  ONLY: xc_lsda
  USE kinds,        ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the input/output arrays
  REAL(DP), INTENT(IN), DIMENSION(length) :: rho_in
  !! total charge density
  REAL(DP), INTENT(IN), DIMENSION(length,3) :: m
  !! magnetization vector
  REAL(DP), INTENT(OUT), DIMENSION(length,4,4) :: dmuxc
  !! derivative of XC functional
  !
  ! ... local variables
  !
  REAL(DP), DIMENSION(length) :: rhotot, amag, zeta, zeta_eff, dr, dz
  REAL(DP), DIMENSION(length) :: vs
  LOGICAL,  DIMENSION(length) :: is_null
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: rhoaux, zetaux
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: aux1, aux2
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: vx, vc
  REAL(DP), DIMENSION(length) :: dvxc_rho, dbx_rho, dby_rho, dbz_rho
  !
  REAL(DP) :: dvxc_mx, dvxc_my, dvxc_mz, &
              dbx_mx, dbx_my, dbx_mz,    &
              dby_mx, dby_my, dby_mz,    &
              dbz_mx, dbz_my, dbz_mz
  REAL(DP) :: zeta_s
  !
  INTEGER :: i1, i2, i3, i4, i5, i
  INTEGER :: f1, f2, f3, f4, f5
  !
  REAL(DP), PARAMETER :: small = 1.E-30_DP, e2 = 2.0_DP, &
                         rho_trash = 0.5_DP, zeta_trash = 0.5_DP, &
                         amag_trash= 0.025_DP
  !
  dmuxc = 0.0_DP
  !
  ALLOCATE( rhoaux(length*5), zetaux(length*5) )
  ALLOCATE( aux1(length*5), aux2(length*5) )
  ALLOCATE( vx(length*5,2), vc(length*5,2) )
  !
  rhotot = rho_in
  zeta   = zeta_trash
  amag   = amag_trash
  is_null = .FALSE.
  !
  i1 = 1     ;   f1 = length    !           five blocks:  [ rho    , zeta    ]   
  i2 = f1+1  ;   f2 = 2*length  !                         [ rho+dr , zeta    ]   
  i3 = f2+1  ;   f3 = 3*length  !                         [ rho-dr , zeta    ]   
  i4 = f3+1  ;   f4 = 4*length  !                         [ rho    , zeta+dz ]   
  i5 = f4+1  ;   f5 = 5*length  !                         [ rho    , zeta-dz ]   
  !
  dz = 1.0E-6_DP     !dz = MIN( 1.d-6, 1.d-4*ABS(zeta) ) 
  !
  DO i = 1, length
     zeta_s = zeta_trash
     IF (rhotot(i) <= small) THEN
        rhotot(i) = rho_trash
        is_null(i) = .TRUE.
     ENDIF
     amag(i) = SQRT( m(i,1)**2 + m(i,2)**2 + m(i,3)**2 )
     IF (rhotot(i) > small) zeta_s = amag(i) / rhotot(i)
     zeta(i) = zeta_s
     zeta_eff(i) = SIGN( MIN( ABS(zeta_s), (1.0_DP-2.0_DP*dz(i)) ), zeta_s )
     IF (ABS(zeta_s) > 1.0_DP) is_null(i) = .TRUE.
  ENDDO
  !
  dr = MIN( 1.E-6_DP, 1.E-4_DP * rhotot )
  !
  rhoaux(i1:f1) = rhotot         ;   zetaux(i1:f1) = zeta
  rhoaux(i2:f2) = rhotot + dr    ;   zetaux(i2:f2) = zeta
  rhoaux(i3:f3) = rhotot - dr    ;   zetaux(i3:f3) = zeta
  rhoaux(i4:f4) = rhotot         ;   zetaux(i4:f4) = zeta_eff + dz
  rhoaux(i5:f5) = rhotot         ;   zetaux(i5:f5) = zeta_eff - dz
  !
  !
  CALL xc_lsda( length*5, rhoaux, zetaux, aux1, aux2, vx, vc )
  !
  !
  vs(:) = 0.5_DP*( vx(i1:f1,1)+vc(i1:f1,1)-vx(i1:f1,2)-vc(i1:f1,2) )
  !
  dvxc_rho(:) = ((vx(i2:f2,1) + vc(i2:f2,1) - vx(i3:f3,1) - vc(i3:f3,1)) + &
                (vx(i2:f2,2) + vc(i2:f2,2) - vx(i3:f3,2) - vc(i3:f3,2))) / (4.0_DP*dr)
  !
  aux2(1:length) =  vx(i2:f2,1) + vc(i2:f2,1) - vx(i3:f3,1) - vc(i3:f3,1) - &
                  ( vx(i2:f2,2) + vc(i2:f2,2) - vx(i3:f3,2) - vc(i3:f3,2) )
  !
  WHERE (amag > 1.E-10_DP)
    dbx_rho(:) = aux2(1:length) * m(:,1) / (4.0_DP*dr*amag)
    dby_rho(:) = aux2(1:length) * m(:,2) / (4.0_DP*dr*amag)
    dbz_rho(:) = aux2(1:length) * m(:,3) / (4.0_DP*dr*amag)
  END WHERE  
  !
  aux1(1:length) =  vx(i4:f4,1) + vc(i4:f4,1) - vx(i5:f5,1) - vc(i5:f5,1) + &
                    vx(i4:f4,2) + vc(i4:f4,2) - vx(i5:f5,2) - vc(i5:f5,2)
  aux2(1:length) =  vx(i4:f4,1) + vc(i4:f4,1) - vx(i5:f5,1) - vc(i5:f5,1) - &
                  ( vx(i4:f4,2) + vc(i4:f4,2) - vx(i5:f5,2) - vc(i5:f5,2) )
  !
  DO i = 1, length
     !
     IF ( is_null(i) ) THEN
        dmuxc(i,:,:) = 0.0_DP
        CYCLE
     ENDIF
     !
     IF (amag(i) <= 1.E-10_DP) THEN
        dmuxc(i,1,1) = dvxc_rho(i)
        CYCLE
     ENDIF
     !
     dvxc_rho(i) = dvxc_rho(i) - aux1(i) * zeta(i)/rhotot(i) / (4.0_DP*dz(i))
     dbx_rho(i)  = dbx_rho(i)  - aux2(i) * m(i,1) * zeta(i)/rhotot(i) / (4.0_DP*dz(i)*amag(i))
     dby_rho(i)  = dby_rho(i)  - aux2(i) * m(i,2) * zeta(i)/rhotot(i) / (4.0_DP*dz(i)*amag(i))
     dbz_rho(i)  = dbz_rho(i)  - aux2(i) * m(i,3) * zeta(i)/rhotot(i) / (4.0_DP*dz(i)*amag(i))
     !
     dmuxc(i,1,1) = dvxc_rho(i)
     dmuxc(i,2,1) = dbx_rho(i)
     dmuxc(i,3,1) = dby_rho(i)
     dmuxc(i,4,1) = dbz_rho(i)
     !
     ! ... Here the derivatives with respect to m
     !
     dvxc_mx = aux1(i) * m(i,1) / rhotot(i) / (4.0_DP*dz(i)*amag(i))
     dvxc_my = aux1(i) * m(i,2) / rhotot(i) / (4.0_DP*dz(i)*amag(i))
     dvxc_mz = aux1(i) * m(i,3) / rhotot(i) / (4.0_DP*dz(i)*amag(i))
     !
     dbx_mx  = (aux2(i) * m(i,1) * m(i,1) * amag(i)/rhotot(i) / (4.0_DP*dz(i)) + &
                vs(i) * (m(i,2)**2+m(i,3)**2)) / amag(i)**3
     dbx_my  = (aux2(i) * m(i,1) * m(i,2) * amag(i)/rhotot(i) / (4.0_DP*dz(i)) - &
                vs(i) * m(i,1) * m(i,2) ) / amag(i)**3
     dbx_mz  = (aux2(i) * m(i,1) * m(i,3) * amag(i)/rhotot(i) / (4.0_DP*dz(i)) - &
                vs(i) * m(i,1) * m(i,3) ) / amag(i)**3
     !
     dby_mx  = dbx_my
     dby_my  = (aux2(i) * m(i,2) * m(i,2) * amag(i)/rhotot(i) / (4.0_DP*dz(i)) + &
                vs(i) * (m(i,1)**2 + m(i,3)**2)) / amag(i)**3
     dby_mz  = (aux2(i) * m(i,2) * m(i,3) * amag(i)/rhotot(i) / (4.0_DP*dz(i)) - &
                vs(i) * m(i,2) * m(i,3)) / amag(i)**3
     !
     dbz_mx  = dbx_mz
     dbz_my  = dby_mz
     dbz_mz  = (aux2(i) * m(i,3) * m(i,3) * amag(i)/rhotot(i) / (4.0_DP*dz(i)) + &
                vs(i)*(m(i,1)**2 + m(i,2)**2)) / amag(i)**3
     !
     ! ... assigns values to dmuxc and sets to zero trash points
     !
     dmuxc(i,1,2) = dvxc_mx 
     dmuxc(i,1,3) = dvxc_my  
     dmuxc(i,1,4) = dvxc_mz 
     !
     dmuxc(i,2,2) = dbx_mx 
     dmuxc(i,2,3) = dbx_my 
     dmuxc(i,2,4) = dbx_mz 
     !
     dmuxc(i,3,2) = dby_mx 
     dmuxc(i,3,3) = dby_my 
     dmuxc(i,3,4) = dby_mz 
     !
     dmuxc(i,4,2) = dbz_mx 
     dmuxc(i,4,3) = dbz_my 
     dmuxc(i,4,4) = dbz_mz 
     !
  ENDDO
  !
  ! ... brings to rydberg units
  !
  dmuxc = e2 * dmuxc
  !
  DEALLOCATE( rhoaux, zetaux)
  DEALLOCATE( aux1, aux2 )
  DEALLOCATE( vx, vc )
  !
  RETURN
  !
END SUBROUTINE dmxc_nc
!
!
!-----------------------------------------------------------------------
FUNCTION dpz( rs, iflg )
  !-----------------------------------------------------------------------
  !!  Derivative of the correlation potential with respect to local density
  !!  Perdew and Zunger parameterization of the Ceperley-Alder functional.
  !
  USE kinds,     ONLY: DP
  USE constants, ONLY: pi, fpi
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rs
  INTEGER,  INTENT(IN) :: iflg
  REAL(DP) :: dpz
  !
  !  ... local variables
  !  a,b,c,d,gc,b1,b2 are the parameters defining the functional
  !
  REAL(DP), PARAMETER :: a = 0.0311d0, b = -0.048d0, c = 0.0020d0, &
       d = -0.0116d0, gc = -0.1423d0, b1 = 1.0529d0, b2 = 0.3334d0,&
       a1 = 7.0d0 * b1 / 6.d0, a2 = 4.d0 * b2 / 3.d0
  REAL(DP) :: x, den, dmx, dmrs
  !
  IF (iflg == 1) THEN
     dmrs = a / rs + 2.d0 / 3.d0 * c * (LOG(rs) + 1.d0) + &
          (2.d0 * d-c) / 3.d0
  ELSE
     x = SQRT(rs)
     den = 1.d0 + x * (b1 + x * b2)
     dmx = gc * ( (a1 + 2.d0 * a2 * x) * den - 2.d0 * (b1 + 2.d0 * &
           b2 * x) * (1.d0 + x * (a1 + x * a2) ) ) / den**3
     dmrs = 0.5d0 * dmx / x
  ENDIF
  !
  dpz = - fpi * rs**4.d0 / 9.d0 * dmrs
  !
  RETURN
  !
END FUNCTION dpz
!
!
!-----------------------------------------------------------------------
FUNCTION dpz_polarized( rs, iflg )
  !-----------------------------------------------------------------------
  !!  Derivative of the correlation potential with respect to local density
  !!  Perdew and Zunger parameterization of the Ceperley-Alder functional.  |
  !!  Spin-polarized case.
  !
  USE kinds,     ONLY: DP
  USE constants, ONLY: pi, fpi
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rs
  INTEGER,  INTENT(IN) :: iflg
  REAL(DP) :: dpz_polarized
  !
  ! ... local variables
  !
  !  a,b,c,d,gc,b1,b2 are the parameters defining the functional
  !
  REAL(DP), PARAMETER :: a=0.01555_DP, b=-0.0269_DP, c=0.0007_DP,  &
                         d=-0.0048_DP, gc=-0.0843_DP, b1=1.3981_DP,&
                         b2=0.2611_DP, a1=7.0_DP*b1/6._DP, a2=4._DP*b2/3._DP
  REAL(DP) :: x, den, dmx, dmrs
  !
  !
  IF (iflg == 1) THEN
     dmrs = a/rs + 2._DP/3._DP * c * (LOG(rs) + 1._DP) + &
            (2._DP*d - c)/3._DP
  ELSE
     x = SQRT(rs)
     den = 1._DP + x * (b1 + x*b2)
     dmx = gc * ( (a1 + 2._DP * a2 * x) * den - 2._DP * (b1 + 2._DP * &
           b2 * x) * (1._DP + x * (a1 + x*a2) ) ) / den**3
     dmrs = 0.5d0 * dmx/x
  ENDIF
  !
  dpz_polarized = - fpi * rs**4._DP / 9._DP * dmrs
  !
  !
  RETURN
  !
END FUNCTION dpz_polarized
