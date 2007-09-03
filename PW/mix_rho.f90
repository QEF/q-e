!
! Copyright (C) 2002-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#define ZERO ( 0.D0, 0.D0 )
#define ONLY_SMOOTH_G
!
!----------------------------------------------------------------------------
SUBROUTINE mix_rho( rhocout, rhocin, taukout, taukin, nsout, nsin, alphamix, &
                    dr2, tr2_min, iter, n_iter, conv )
  !----------------------------------------------------------------------------
  !
  ! ... Modified Broyden's method for charge density mixing
  ! ...         D.D. Johnson PRB 38, 12807 (1988)
  !
  ! ... On output: the mixed density is in rhocin
  !
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat
  USE gvect,          ONLY : ngm
  USE gsmooth,        ONLY : ngms
  USE ldaU,           ONLY : lda_plus_u, Hubbard_lmax
  USE funct,          ONLY : dft_is_meta
  USE lsda_mod,       ONLY : nspin
  USE control_flags,  ONLY : imix, ngm0, tr2, io_level
  USE io_files,       ONLY : find_free_unit
  USE cell_base,      ONLY : omega
  !
  IMPLICIT NONE
  !
  ! ... First the I/O variable
  !
  INTEGER :: &
    iter,                  &!  (in)  counter of the number of iterations
    n_iter                  !  (in)  numb. of iterations used in mixing
  COMPLEX(DP) :: &
    rhocin(ngm,nspin), rhocout(ngm,nspin)
  COMPLEX(DP) :: &
    taukin(ngm,nspin), taukout(ngm,nspin)
  REAL(DP) :: &
    nsout(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat), &!
    nsin(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat),  &!
    alphamix,              &! (in) mixing factor
    dr2                     ! (out) the estimated errr on the energy
  REAL(DP) :: &
    tr2_min       ! estimated error from diagonalization. If the estimated
                  ! scf error is smaller than this, exit: a more accurate 
                  ! diagonalization is needed
  LOGICAL :: &
    conv          ! (out) if true the convergence has been reached
  INTEGER, PARAMETER :: &
    maxmix = 25             ! max number of iterations for charge mixing
  !
  ! ... Here the local variables
  !
  LOGICAL :: &
    tmeta         ! set true if dft_is_meta
  INTEGER ::    &
    iunmix,        &! I/O unit number of charge density file
    iunmix2,       &! I/O unit number of tauk file
    iunmix3,       &! I/O unit number of ns file
    iter_used,     &! actual number of iterations used
    ipos,          &! index of the present iteration
    inext,         &! index of the next iteration
    i, j,          &! counters on number of iterations
    iwork(maxmix), &! dummy array used as output by libr. routines
    info,          &! flag saying if the exec. of libr. routines was ok
    ldim            ! 2 * Hubbard_lmax + 1
  COMPLEX(DP), ALLOCATABLE :: &
    rhoinsave(:,:),     &! rhoinsave(ngm0,nspin): work space
    rhoutsave(:,:),     &! rhoutsave(ngm0,nspin): work space
    taukinsave(:,:),    &! 
    taukoutsave(:,:),   &! 
    nsinsave(:,:,:,:),  &!
    nsoutsave(:,:,:,:)   !
  REAL(DP) :: &
    betamix(maxmix,maxmix), &
    gamma0,                 &
    work(maxmix)
  LOGICAL :: &
    savetofile,  &! save intermediate steps on file $prefix."mix","mix.tauk"...
    exst          ! if true the file exists
  !
  ! ... saved variables and arrays
  !
  INTEGER, SAVE :: &
    mixrho_iter = 0    ! history of mixing
  COMPLEX(DP), ALLOCATABLE, SAVE :: &
    df(:,:,:),        &! information from preceding iterations
    dv(:,:,:)          !     "  "       "     "        "  "
  COMPLEX(DP), ALLOCATABLE, SAVE :: &
    df_k(:,:,:),      &! idem
    dv_k(:,:,:)        ! idem
  REAL(DP), ALLOCATABLE, SAVE :: &
    df_ns(:,:,:,:,:), &! idem 
    dv_ns(:,:,:,:,:)   ! idem
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: rho_ddot, ns_ddot, tauk_ddot
  !
  !
  CALL start_clock( 'mix_rho' )
  !
#if defined (ONLY_SMOOTH_G)
  ngm0 = ngms
#else
  ngm0 = ngm
#endif
  !
  mixrho_iter = iter
  !
  IF ( n_iter > maxmix ) CALL errore( 'mix_rho', 'n_iter too big', 1 )
  !
  IF ( lda_plus_u ) ldim = 2 * Hubbard_lmax + 1
  !
  savetofile = (io_level > 1)
  !
  rhocout(:,:) = rhocout(:,:) - rhocin(:,:)
  tmeta = dft_is_meta()
  IF (tmeta) taukout(:,:) = taukout(:,:) - taukin(:,:)
  IF ( lda_plus_u ) nsout(:,:,:,:) = nsout(:,:,:,:) - nsin(:,:,:,:)
  !
  dr2 = rho_ddot( rhocout, rhocout, ngm, ngm, nspin, ngm )
  IF ( tmeta ) dr2 = dr2 + tauk_ddot( taukout, taukout, ngm, ngm, nspin, ngm )
  IF ( lda_plus_u ) dr2 = dr2 + ns_ddot( nsout, nsout, nspin )
  !
  conv = ( dr2 < tr2 )
  !
  IF ( conv .OR. dr2 < tr2_min ) THEN
     !
     ! ... if convergence is achieved or if the self-consistency error (dr2) is
     ! ... smaller than the estimated error due to diagonalization (tr2_min),
     ! ... exit and leave rhocin and rhocout unchanged
     !
     IF ( ALLOCATED( df ) )    DEALLOCATE( df )
     IF ( ALLOCATED( dv ) )    DEALLOCATE( dv )
     IF ( ALLOCATED( df_k ) )  DEALLOCATE( df_k )
     IF ( ALLOCATED( dv_k ) )  DEALLOCATE( dv_k )
     IF ( ALLOCATED( df_ns ) ) DEALLOCATE( df_ns )
     IF ( ALLOCATED( dv_ns ) ) DEALLOCATE( dv_ns )
     !
     rhocout(:,:) = rhocout(:,:) + rhocin(:,:)
     IF ( tmeta ) taukout(:,:) = taukout(:,:) + taukin(:,:)
     IF ( lda_plus_u ) nsout(:,:,:,:) = nsout(:,:,:,:) + nsin(:,:,:,:)
     !
     CALL stop_clock( 'mix_rho' )
     !
     RETURN
     !
  END IF
  !
  IF ( savetofile ) THEN
     !
     iunmix = find_free_unit()
     CALL diropn( iunmix, 'mix', 2*ngm0*nspin, exst )
     !
     IF ( tmeta ) THEN
        iunmix2 = find_free_unit()
        CALL diropn( iunmix2, 'mix.tauk', 2*ngm0*nspin, exst )
     END IF
     !
     IF ( lda_plus_u ) THEN
        iunmix3 = find_free_unit()
        CALL diropn( iunmix3, 'mix.ns', ldim*ldim*nspin*nat, exst )
     END IF
     !
     IF ( mixrho_iter > 1 .AND. .NOT. exst ) THEN
        !
        CALL infomsg( 'mix_rho', 'file not found, restarting' )
        mixrho_iter = 1
        !
     END IF
     !
  END IF
  !
  IF ( savetofile .OR. mixrho_iter == 1 ) THEN
     !
     IF ( .NOT. ALLOCATED( df ) ) ALLOCATE( df( ngm0, nspin, n_iter ) )
     IF ( .NOT. ALLOCATED( dv ) ) ALLOCATE( dv( ngm0, nspin, n_iter ) )
     !
     IF ( tmeta ) THEN
        IF ( .NOT. ALLOCATED( df_k ) ) ALLOCATE( df_k( ngm0, nspin, n_iter ) )
        IF ( .NOT. ALLOCATED( dv_k ) ) ALLOCATE( dv_k( ngm0, nspin, n_iter ) )
     END IF
     !
     IF ( lda_plus_u ) THEN
        IF ( .NOT. ALLOCATED( df_ns ) ) &
           ALLOCATE( df_ns( ldim, ldim, nspin, nat, n_iter ) )
        IF ( .NOT. ALLOCATED( dv_ns ) ) &
           ALLOCATE( dv_ns( ldim, ldim, nspin, nat, n_iter ) )
     END IF
     !
  END IF
  !
  ! ... iter_used = mixrho_iter-1  if  mixrho_iter <= n_iter
  ! ... iter_used = n_iter         if  mixrho_iter >  n_iter
  !
  iter_used = MIN( ( mixrho_iter - 1 ), n_iter )
  !
  ! ... ipos is the position in which results from the present iteration
  ! ... are stored. ipos=mixrho_iter-1 until ipos=n_iter, then back to 1,2,...
  !
  ipos = mixrho_iter - 1 - ( ( mixrho_iter - 2 ) / n_iter ) * n_iter
  !
  IF ( mixrho_iter > 1 ) THEN
     !
     IF ( savetofile ) THEN
        !
        CALL davcio( df(1,1,ipos), 2*ngm0*nspin, iunmix, 1, -1 )
        CALL davcio( dv(1,1,ipos), 2*ngm0*nspin, iunmix, 2, -1 )
        !
        IF ( tmeta ) THEN
           !
           CALL davcio( df_k(1,1,ipos), 2*ngm0*nspin, iunmix2, 1, -1 )
           CALL davcio( dv_k(1,1,ipos), 2*ngm0*nspin, iunmix2, 2, -1 )
           !
        END IF
        IF ( lda_plus_u ) THEN
           !
           CALL davcio( df_ns(1,1,1,1,ipos),ldim*ldim*nspin*nat,iunmix3,1,-1 )
           CALL davcio( dv_ns(1,1,1,1,ipos),ldim*ldim*nspin*nat,iunmix3,2,-1 )
           !
        END IF
        !
     END IF
     !
     df(:,:,ipos) = df(:,:,ipos) - rhocout(1:ngm0,:)
     dv(:,:,ipos) = dv(:,:,ipos) - rhocin (1:ngm0,:)
     !
     IF ( tmeta) THEN
        df_k(:,:,ipos) = df_k(:,:,ipos) - taukout(1:ngm0,:)
        dv_k(:,:,ipos) = dv_k(:,:,ipos) - taukin (1:ngm0,:)
     END IF
     !
     IF ( lda_plus_u ) THEN
        !
        df_ns(:,:,:,:,ipos) = df_ns(:,:,:,:,ipos) - nsout
        dv_ns(:,:,:,:,ipos) = dv_ns(:,:,:,:,ipos) - nsin
        !
     END IF
     !
  END IF
  !
  IF ( savetofile ) THEN
     !
     DO i = 1, iter_used
        !
        IF ( i /= ipos ) THEN
           !
           CALL davcio( df(1,1,i), 2*ngm0*nspin, iunmix, 2*i+1, -1 )
           CALL davcio( dv(1,1,i), 2*ngm0*nspin, iunmix, 2*i+2, -1 )
           !
           IF ( tmeta ) THEN
              CALL davcio( df_k(1,1,i), 2*ngm0*nspin, iunmix2, 2*i+1, -1 )
              CALL davcio( dv_k(1,1,i), 2*ngm0*nspin, iunmix2, 2*i+2, -1 )
           END IF
           !
           IF ( lda_plus_u ) THEN
              CALL davcio(df_ns(1,1,1,1,i),ldim*ldim*nspin*nat,iunmix3,2*i+1,-1)
              CALL davcio(dv_ns(1,1,1,1,i),ldim*ldim*nspin*nat,iunmix3,2*i+2,-1)
           END IF
           !
        END IF
        !
     END DO
     ! not good if ngm0 != ngm 
     !! CALL davcio( rhocout, 2*ngm0*nspin, iunmix, 1, 1 )
     !! CALL davcio( rhocin,  2*ngm0*nspin, iunmix, 2, 1 )
     ! use instead:
     ALLOCATE( rhoinsave( ngm0, nspin ) )
     rhoinsave(:,:) = rhocout(1:ngm0,:)
     CALL davcio( rhoinsave, 2*ngm0*nspin, iunmix, 1, 1 )
     rhoinsave(:,:) = rhocin(1:ngm0,:)
     CALL davcio( rhoinsave,  2*ngm0*nspin, iunmix, 2, 1 )
     !
     IF ( mixrho_iter > 1 ) THEN
        !
        CALL davcio( df(1,1,ipos), 2*ngm0*nspin, iunmix, 2*ipos+1, 1 )
        CALL davcio( dv(1,1,ipos), 2*ngm0*nspin, iunmix, 2*ipos+2, 1 )
        !
     END IF
     !
     IF ( tmeta ) THEN
        !
        !! CALL davcio( taukout, 2*ngm0*nspin, iunmix2, 1, 1 )
        !! CALL davcio( taukin,  2*ngm0*nspin, iunmix2, 2, 1 )
        rhoinsave(:,:) = taukout(1:ngm0,:)
        CALL davcio( rhoinsave, 2*ngm0*nspin, iunmix2, 1, 1 )
        rhoinsave(:,:) = taukin(1:ngm0,:)
        CALL davcio( rhoinsave,  2*ngm0*nspin, iunmix2, 2, 1 )
        !
        IF ( mixrho_iter > 1 ) THEN
           !
           CALL davcio( df_k(1,1,ipos), 2*ngm0*nspin, iunmix2, 2*ipos+1, 1 )
           CALL davcio( dv_k(1,1,ipos), 2*ngm0*nspin, iunmix2, 2*ipos+2, 1 )
           !
        END IF
        !
     END IF
     !
     DEALLOCATE( rhoinsave )
     !
     IF ( lda_plus_u ) THEN
        !
        CALL davcio( nsout, ldim*ldim*nspin*nat, iunmix3, 1, 1 )
        CALL davcio( nsin , ldim*ldim*nspin*nat, iunmix3, 2, 1 )
        !
        IF ( mixrho_iter > 1 ) THEN
           !
           CALL davcio( df_ns(1,1,1,1,ipos), &
                        ldim*ldim*nspin*nat, iunmix3, 2*ipos+1, 1 )
           CALL davcio( dv_ns(1,1,1,1,ipos), &
                        ldim*ldim*nspin*nat, iunmix3, 2*ipos+2, 1 )
        END IF
        !
     END IF
     !
  ELSE
     !
     ALLOCATE( rhoinsave( ngm0, nspin ), rhoutsave( ngm0, nspin ) )
     !
     rhoinsave = rhocin(1:ngm0,:)
     rhoutsave = rhocout(1:ngm0,:)
     !
     IF ( tmeta ) THEN
        !
        ALLOCATE( taukinsave( ngm0, nspin ), taukoutsave( ngm0, nspin ) )
        !
        taukinsave = taukin(1:ngm0,:)
        taukoutsave = taukout(1:ngm0,:)
        !
     END IF
     !
     IF ( lda_plus_u ) THEN
        !
        ALLOCATE( nsinsave(  ldim, ldim, nspin, nat ), &
                  nsoutsave( ldim, ldim, nspin, nat ) )
        nsinsave  = nsin
        nsoutsave = nsout
        !
     END IF
     !
  END IF
  !
  DO i = 1, iter_used
     !
     DO j = i, iter_used
        !
        betamix(i,j) = rho_ddot( df(1,1,j), df(1,1,i), ngm0, ngm0, nspin, ngm0 )
        !
        IF ( tmeta ) betamix(i,j) = betamix(i,j) + &
                     tauk_ddot( df_k(1,1,j), df_k(1,1,i), ngm0,ngm0,nspin,ngm0 )
        !
        IF ( lda_plus_u ) betamix(i,j) = betamix(i,j) + &
                          ns_ddot( df_ns(1,1,1,1,j), df_ns(1,1,1,1,i), nspin )
        !
        betamix(j,i) = betamix(i,j)
        !
     END DO
     !
  END DO
  !
  CALL DSYTRF( 'U', iter_used, betamix, maxmix, iwork, work, maxmix, info )
  CALL errore( 'broyden', 'factorization', info )
  !
  CALL DSYTRI( 'U', iter_used, betamix, maxmix, iwork, work, info )
  CALL errore( 'broyden', 'DSYTRI', info )
  !
  FORALL( i = 1:iter_used, &
          j = 1:iter_used, j > i ) betamix(j,i) = betamix(i,j)
  !
  DO i = 1, iter_used
     !
     work(i) = rho_ddot( df(1,1,i), rhocout, ngm0, ngm, nspin, ngm0 )
     !
     IF ( tmeta ) &
        work(i) = tauk_ddot( df_k(1,1,i), taukout, ngm0, ngm, nspin, ngm0 )
     !
     IF ( lda_plus_u ) &
        work(i) = work(i) + ns_ddot( df_ns(1,1,1,1,i), nsout, nspin )
     !
  END DO
  !
  DO i = 1, iter_used
     !
     gamma0 = SUM( betamix(1:iter_used,i)*work(1:iter_used) )
     !
     rhocin(1:ngm0,:)  = rhocin(1:ngm0,:)  - gamma0*dv(:,:,i)
     rhocout(1:ngm0,:) = rhocout(1:ngm0,:) - gamma0*df(:,:,i)
     !
     IF ( tmeta ) THEN
        taukin(1:ngm0,:)  = taukin(1:ngm0,:)  - gamma0*dv_k(:,:,i)
        taukout(1:ngm0,:) = taukout(1:ngm0,:) - gamma0*df_k(:,:,i)
     END IF
     !
     IF ( lda_plus_u ) THEN
        !
        nsin  = nsin  - gamma0*dv_ns(:,:,:,:,i)
        nsout = nsout - gamma0*df_ns(:,:,:,:,i)
        !
     END IF
     !
  END DO
  !
  ! ... auxiliary vectors dv and df not needed anymore
  !
  IF ( savetofile ) THEN
     !
     IF ( lda_plus_u ) THEN
        !
        CLOSE( iunmix3, STATUS = 'KEEP' )
        !
        DEALLOCATE( df_ns, dv_ns )
        !
     END IF
     !
     IF ( tmeta ) THEN
        !
        CLOSE( iunmix2, STATUS = 'KEEP' )
        !
        DEALLOCATE( df_k, dv_k )
        !
     END IF
     !
     CLOSE( iunmix, STATUS = 'KEEP' )
     !
     DEALLOCATE( df, dv )
     !
  ELSE
     !
     inext = mixrho_iter - ( ( mixrho_iter - 1 ) / n_iter ) * n_iter
     !
     IF ( lda_plus_u ) THEN
        !
        df_ns(:,:,:,:,inext) = nsoutsave
        dv_ns(:,:,:,:,inext) = nsinsave
        !
        DEALLOCATE( nsinsave, nsoutsave )
        !
     END IF
     !
     IF ( tmeta ) THEN
        !
        df_k(:,:,inext) = taukoutsave(:,:)
        dv_k(:,:,inext) = taukinsave(:,:)
        !
        DEALLOCATE( taukinsave, taukoutsave )
        !
     END IF
     !
     df(:,:,inext) = rhoutsave(:,:)
     dv(:,:,inext) = rhoinsave(:,:)
     !
     DEALLOCATE( rhoinsave, rhoutsave )
     !
  END IF
  !
  ! ... preconditioning the new search direction
  !
  IF ( imix == 1 ) THEN
     !
     CALL approx_screening( rhocout )
     !
  ELSE IF ( imix == 2 ) THEN
     !
     CALL approx_screening2( rhocout, rhocin )
     !
  END IF
  !
  ! ... set new trial density
  !
  rhocin = rhocin + alphamix * rhocout
  !
  IF ( tmeta ) taukin = taukin + alphamix * taukout
  !
  IF ( lda_plus_u ) nsin = nsin + alphamix * nsout
  !
  CALL stop_clock( 'mix_rho' )
  !
  RETURN
  !
END SUBROUTINE mix_rho
!
!----------------------------------------------------------------------------
FUNCTION rho_ddot( rho1, rho2, ngm1, ngm2, nspin, gf )
  !----------------------------------------------------------------------------
  !
  ! ... calculates 4pi/G^2*rho1(-G)*rho2(G) = V1_Hartree(-G)*rho2(G)
  ! ... used as an estimate of the self-consistency error on the energy
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, tpi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gg, gstart
  USE wvfct,         ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN) :: ngm1, ngm2, nspin, gf
  COMPLEX(DP), INTENT(IN) :: rho1(ngm1,nspin), rho2(ngm2,nspin)
  REAL(DP)                :: rho_ddot
  !
  REAL(DP) :: fac
  INTEGER  :: ig
  !
  !
  fac = e2 * fpi / tpiba2
  !
  rho_ddot = 0.D0
  !
  IF ( nspin == 1 ) THEN
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + REAL( CONJG( rho1(ig,1) )*rho2(ig,1) ) / gg(ig)
        !
     END DO
     !
     rho_ddot = fac*rho_ddot
     !
     IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     ! ... first the charge
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + REAL( CONJG( rho1(ig,1)+rho1(ig,2) ) * &
                                         ( rho2(ig,1)+rho2(ig,2) ) ) / gg(ig)
        !
     END DO
     !
     rho_ddot = fac*rho_ddot
     !
     IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
     !
     ! ... then the magnetization
     !
     fac = e2 * fpi / tpi**2  ! lambda = 1 a.u.
     !
     ! ... G=0 term
     !
     IF ( gstart == 2 ) THEN
        !
        rho_ddot = rho_ddot + fac * REAL( CONJG( rho1(1,1) - rho1(1,2) ) * &
                                               ( rho2(1,1) - rho2(1,2) ) )
        !
     END IF
     !
     IF ( gamma_only ) fac = 2.D0 * fac
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + fac * REAL( CONJG( rho1(ig,1) - rho1(ig,2) ) * &
                                               ( rho2(ig,1) - rho2(ig,2) ) )
        !
     END DO
     !
  ELSE IF ( nspin == 4 ) THEN
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + REAL( CONJG( rho1(ig,1) )*rho2(ig,1) ) / gg(ig)
        !
     END DO
     !
     rho_ddot = fac*rho_ddot
     !
     IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
     !
     fac = e2*fpi / (tpi**2)  ! lambda=1 a.u.
     !
     IF ( gstart == 2 ) THEN
        !
        rho_ddot = rho_ddot + &
                   fac * ( REAL( CONJG( rho1(1,2))*(rho2(1,2) ) ) + &
                           REAL( CONJG( rho1(1,3))*(rho2(1,3) ) ) + &
                           REAL( CONJG( rho1(1,4))*(rho2(1,4) ) ) )
        !
     END IF
     !
     IF ( gamma_only ) fac = 2.D0 * fac
     !
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + &
                   fac * ( REAL( CONJG( rho1(ig,2))*(rho2(ig,2) ) ) + &
                           REAL( CONJG( rho1(ig,3))*(rho2(ig,3) ) ) + &
                           REAL( CONJG( rho1(ig,4))*(rho2(ig,4) ) ) )
        !
     END DO
     !
  END IF
  !
  rho_ddot = rho_ddot * omega * 0.5D0
  !
  CALL reduce( 1, rho_ddot )
  !
  RETURN
  !
END FUNCTION rho_ddot
!
!----------------------------------------------------------------------------
FUNCTION tauk_ddot( tauk1, tauk2, ngm1, ngm2, nspin, gf )
  !----------------------------------------------------------------------------
  !
  ! ... calculates 4pi/G^2*rho1(-G)*rho2(G) = V1_Hartree(-G)*rho2(G)
  ! ... used as an estimate of the self-consistency error on the energy
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, tpi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gg, gstart
  USE wvfct,         ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN) :: ngm1, ngm2, nspin, gf
  COMPLEX(DP), INTENT(IN) :: tauk1(ngm1,nspin), tauk2(ngm2,nspin)
  REAL(DP)                :: tauk_ddot
  !
  REAL(DP) :: fac
  INTEGER  :: ig
  !
  tauk_ddot = 0.D0
  !
!  if (.true. ) return
  IF ( nspin == 1 ) THEN
     !
     DO ig = gstart, gf
        tauk_ddot = tauk_ddot + REAL( CONJG( tauk1(ig,1) )*tauk2(ig,1) ) 
     END DO
     !
     IF ( gamma_only ) tauk_ddot = 2.D0 * tauk_ddot
     !
     ! ... G=0 term
     !
     IF ( gstart == 2 ) THEN
        !
        tauk_ddot = tauk_ddot + REAL( CONJG( tauk1(1,1) ) * tauk2(1,1) )
        !
     END IF
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     DO ig = gstart, gf
        !
        tauk_ddot = tauk_ddot + REAL( CONJG( tauk1(ig,1) + tauk1(ig,2) ) * &
                                           ( tauk2(ig,1) + tauk2(ig,2) ) ) 
        tauk_ddot = tauk_ddot + REAL( CONJG( tauk1(ig,1) - tauk1(ig,2) ) * &
                                           ( tauk2(ig,1) - tauk2(ig,2) ) )
        !
     END DO
     !
     IF ( gamma_only ) tauk_ddot = 2.D0 * tauk_ddot
     !
     ! ... G=0 term
     !
     IF ( gstart == 2 ) THEN
        !
        tauk_ddot = tauk_ddot + REAL( CONJG( tauk1(1,1) + tauk1(1,2) ) * &
                                           ( tauk2(1,1) + tauk2(1,2) ) )
        tauk_ddot = tauk_ddot + REAL( CONJG( tauk1(1,1) - tauk1(1,2) ) * &
                                           ( tauk2(1,1) - tauk2(1,2) ) )
        !
     END IF
     !
  ELSE IF ( nspin == 4 ) THEN
     !
     DO ig = gstart, gf
        !
        tauk_ddot = tauk_ddot + &
                         ( REAL( CONJG( tauk1(ig,1))*(tauk2(ig,1) ) ) + &
                           REAL( CONJG( tauk1(ig,2))*(tauk2(ig,2) ) ) + &
                           REAL( CONJG( tauk1(ig,3))*(tauk2(ig,3) ) ) + &
                           REAL( CONJG( tauk1(ig,4))*(tauk2(ig,4) ) ) )
        !
     END DO
     !
     IF ( gamma_only ) tauk_ddot = 2.D0 * tauk_ddot
     !
     IF ( gstart == 2 ) THEN
        !
        tauk_ddot = tauk_ddot + &
                         ( REAL( CONJG( tauk1(1,1))*(tauk2(1,1) ) ) + &
                           REAL( CONJG( tauk1(1,2))*(tauk2(1,2) ) ) + &
                           REAL( CONJG( tauk1(1,3))*(tauk2(1,3) ) ) + &
                           REAL( CONJG( tauk1(1,4))*(tauk2(1,4) ) ) )
        !
     END IF
     !
  END IF
  !
  fac = e2 * fpi / tpi**2  ! lambda = 1 a.u.
  !
  tauk_ddot = fac * tauk_ddot * omega * 0.5D0
  !
  CALL reduce( 1, tauk_ddot )
  !
  RETURN
  !
END FUNCTION tauk_ddot
!
!----------------------------------------------------------------------------
FUNCTION ns_ddot( ns1, ns2, nspin )
  !----------------------------------------------------------------------------
  !
  ! ... calculates U/2 \sum_i ns1(i)*ns2(i)
  ! ... used as an estimate of the self-consistency error on the 
  ! ... LDA+U correction to the energy
  !
  USE kinds,     ONLY : DP
  USE ldaU,      ONLY : lda_plus_u, Hubbard_lmax, Hubbard_l, Hubbard_U, &
                        Hubbard_alpha
  USE ions_base, ONLY : nat, ityp
  !
  IMPLICIT NONE  
  !
  INTEGER,  INTENT(IN) :: nspin
  REAL(DP), INTENT(IN) :: ns1(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat), &
                          ns2(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  REAL(DP)             :: ns_ddot
  !
  INTEGER :: na, nt, m1, m2
  !
  !
  ns_ddot = 0.D0
  !
  IF ( .NOT. lda_plus_u ) RETURN
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     IF ( Hubbard_U(nt) /= 0.D0 .OR. Hubbard_alpha(nt) /= 0.D0 ) THEN
        !
        m1 = 2 * Hubbard_l(nt) + 1
        m2 = 2 * Hubbard_l(nt) + 1
        !
        ns_ddot = ns_ddot + 0.5D0 * Hubbard_U(nt) * &
                            SUM( ns1(:m1,:m2,:nspin,na)*ns2(:m1,:m2,:nspin,na) )
        !
     END IF
     !
  END DO
  !
  IF ( nspin == 1 ) ns_ddot = 2.D0*ns_ddot
  !
  RETURN
  !
END FUNCTION ns_ddot
!
!----------------------------------------------------------------------------
SUBROUTINE approx_screening( drho )
  !----------------------------------------------------------------------------
  !
  ! ... apply an average TF preconditioning to drho
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, pi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gg, ngm
  USE klist,         ONLY : nelec
  USE lsda_mod,      ONLY : nspin
  USE control_flags, ONLY : ngm0
  !
  IMPLICIT NONE  
  !
  COMPLEX(DP), INTENT(INOUT) :: drho(ngm,nspin) ! (in/out)
  !
  REAL(DP) :: rrho, rmag, rs, agg0
  INTEGER  :: ig
  !
  !
  rs = ( 3.D0 * omega / fpi / nelec )**( 1.D0 / 3.D0 )
  !
  agg0 = ( 12.D0 / pi )**( 2.D0 / 3.D0 ) / tpiba2 / rs
  !
  IF ( nspin == 1 .OR. nspin == 4 ) THEN
     !
     drho(:ngm0,1) =  drho(:ngm0,1) * gg(:ngm0) / ( gg(:ngm0) + agg0 )
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     DO ig = 1, ngm0
        !
        rrho = ( drho(ig,1) + drho(ig,2) ) * gg(ig) / ( gg(ig) + agg0 )
        rmag = ( drho(ig,1) - drho(ig,2) )
        !
        drho(ig,1) =  0.5D0*( rrho + rmag )
        drho(ig,2) =  0.5D0*( rrho - rmag )
        !
     END DO
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE approx_screening
!
#if defined (ONLY_SMOOTH_G)
!----------------------------------------------------------------------------
SUBROUTINE approx_screening2( drho, rhobest )
  !----------------------------------------------------------------------------
  !
  ! ... apply a local-density dependent TF preconditioning to drho
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : e2, pi, tpi, fpi, eps8, eps32
  USE cell_base,            ONLY : omega, tpiba2
  USE gsmooth,              ONLY : nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, &
                                   nrxxs, nls, nlsm
  USE gvect,                ONLY : gg, ngm
  USE wavefunctions_module, ONLY : psic
  USE klist,                ONLY : nelec
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : gamma_only
  USE control_flags,        ONLY : ngm0
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: drho(ngm,nspin)
  COMPLEX(DP), INTENT(IN)    :: rhobest(ngm,nspin)
  !
  INTEGER, PARAMETER :: mmx = 12
  !
  INTEGER :: &
    iwork(mmx), i, j, m, info
  REAL(DP) :: &
    rs, min_rs, max_rs, avg_rsm1, target, dr2_best
  REAL(DP) :: &
    aa(mmx,mmx), invaa(mmx,mmx), bb(mmx), work(mmx), vec(mmx), agg0
  COMPLEX(DP), ALLOCATABLE :: &
    v(:,:),     &! v(ngm0,mmx)
    w(:,:),     &! w(ngm0,mmx)
    dv(:),      &! dv(ngm0)
    vbest(:),   &! vbest(ngm0)
    wbest(:)     ! wbest(ngm0)
  REAL(DP), ALLOCATABLE :: &
    alpha(:)     ! alpha(nrxxs)
  !
  COMPLEX(DP)         :: rrho, rmag
  INTEGER             :: ir, ig
  REAL(DP), PARAMETER :: one_third = 1.D0 / 3.D0
  REAL(DP), EXTERNAL  :: rho_ddot
  !
  !
  IF ( nspin == 2 ) THEN
     !
     DO ig = 1, ngm0
        !
        rrho = drho(ig,1) + drho(ig,2)
        rmag = drho(ig,1) - drho(ig,2)
        !        
        drho(ig,1) = rrho
        drho(ig,2) = rmag
        !
     END DO
     !
  END IF
  !
  target = 0.D0
  !
  IF ( gg(1) < eps8 ) drho(1,1) = ZERO
  !
  ALLOCATE( alpha( nrxxs ) )
  ALLOCATE( v( ngm0, mmx ), &
            w( ngm0, mmx ), dv( ngm0 ), vbest( ngm0 ), wbest( ngm0 ) )
  !
  v(:,:)   = ZERO
  w(:,:)   = ZERO
  dv(:)    = ZERO
  vbest(:) = ZERO
  wbest(:) = ZERO
  !
  ! ... calculate alpha from density
  !
  psic(:) = ZERO
  !
  IF ( nspin == 2 ) THEN
     !
     psic(nls(:ngm0)) = ( rhobest(:ngm0,1) + rhobest(:ngm0,2) )
     !
  ELSE
     !
     psic(nls(:ngm0)) = rhobest(:ngm0,1)
     !
  END IF
  !
  IF ( gamma_only ) psic(nlsm(:ngm0)) = CONJG( psic(nls(:ngm0)) )
  !
  CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 1 )
  !
  alpha(:) = REAL( psic(1:nrxxs) )
  !
  min_rs   = ( 3.D0 * omega / fpi / nelec )**one_third
  max_rs   = min_rs
  avg_rsm1 = 0.D0
  !
  DO ir = 1, nrxxs
     !
     alpha(ir) = ABS( alpha(ir) )
     !
     IF ( alpha(ir) > eps32 ) THEN
        !
        rs        = ( 3.D0 / fpi / alpha(ir) )**one_third
        min_rs    = MIN( min_rs, rs )
        avg_rsm1  = avg_rsm1 + 1.D0 / rs
        max_rs    = MAX( max_rs, rs )
        alpha(ir) = rs
        !
     END IF   
     !
  END DO
  !
  CALL reduce( 1, avg_rsm1 )
  !
  CALL extreme( min_rs, -1 )
  CALL extreme( max_rs, +1 )
  !
  alpha = 3.D0 * ( tpi / 3.D0 )**( 5.D0 / 3.D0 ) * alpha
  !
  avg_rsm1 = ( nr1s*nr2s*nr3s ) / avg_rsm1
  rs       = ( 3.D0 * omega / fpi / nelec )**one_third
  agg0     = ( 12.D0 / pi )**( 2.D0 / 3.D0 ) / tpiba2 / avg_rsm1
  !
  ! ... calculate deltaV and the first correction vector
  !
  psic(:) = ZERO
  !
  psic(nls(:ngm0)) = drho(:ngm0,1)
  !
  IF ( gamma_only ) psic(nlsm(:ngm0)) = CONJG( psic(nls(:ngm0)) )
  !
  CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +1 )
  !
  psic(:nrxxs) = psic(:nrxxs) * alpha(:)
  !
  CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1 )
  !
  dv(:) = psic(nls(:ngm0)) * gg(:ngm0) * tpiba2
  v(:,1)= psic(nls(:ngm0)) * gg(:ngm0) / ( gg(:ngm0) + agg0 )
  !
  m       = 1
  aa(:,:) = 0.D0
  bb(:)   = 0.D0
  !
  repeat_loop: DO
     !
     ! ... generate the vector w
     !     
     w(:,m) = fpi * e2 * v(:,m)
     !
     psic(:) = ZERO
     !
     psic(nls(:ngm0)) = v(:,m)
     !
     IF ( gamma_only ) psic(nlsm(:ngm0)) = CONJG( psic(nls(:ngm0)) )
     !
     CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +1 )
     !
     psic(:nrxxs) = psic(:nrxxs) * alpha(:)
     !
     CALL cft3s( psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1 )
     !
     w(:,m) = w(:,m) + gg(:ngm0) * tpiba2 * psic(nls(:ngm0))
     !
     ! ... build the linear system
     !
     DO i = 1, m
        !
        aa(i,m) = rho_ddot( w(1,i), w(1,m), ngm0, ngm0, 1, ngm0)
        !
        aa(m,i) = aa(i,m)
        !
     END DO
     !
     bb(m) = rho_ddot( w(1,m), dv, ngm0, ngm0, 1, ngm0 )
     !
     ! ... solve it -> vec
     !
     invaa = aa
     !
     CALL DSYTRF( 'U', m, invaa, mmx, iwork, work, mmx, info )
     CALL errore( 'broyden', 'factorization', info )
     !
     CALL DSYTRI( 'U', m, invaa, mmx, iwork, work, info )
     CALL errore( 'broyden', 'DSYTRI', info )
     !     
     FORALL( i = 1:m, j = 1:m, j > i ) invaa(j,i) = invaa(i,j)
     !
     FORALL( i = 1:m ) vec(i) = SUM( invaa(i,:)*bb(:) )
     !
     vbest(:) = ZERO
     wbest(:) = dv(:)
     !
     DO i = 1, m
        !
        vbest = vbest + vec(i) * v(:,i)
        wbest = wbest - vec(i) * w(:,i)
        !
     END DO
     !
     dr2_best = rho_ddot( wbest, wbest, ngm0, ngm0, 1, ngm0 )
     !
     IF ( target == 0.D0 ) target = 1.D-6 * dr2_best
     !
     IF ( dr2_best < target ) THEN
        !
        drho(:ngm0,1) = vbest(:)
        !
        IF ( nspin == 2 ) THEN
           !
           DO ig = 1, ngm0
              !
              rrho = drho(ig,1)
              rmag = drho(ig,2)
              !
              drho(ig,1) = 0.5D0 * ( rrho + rmag )
              drho(ig,2) = 0.5D0 * ( rrho - rmag )
              !
           END DO
           !
        END IF
        !
        DEALLOCATE( alpha, v, w, dv, vbest, wbest )
        !
        EXIT repeat_loop
        !
     ELSE IF ( m >= mmx ) THEN
        !
        m = 1
        !
        v(:,m)  = vbest(:)
        aa(:,:) = 0.D0
        bb(:)   = 0.D0
        !
        CYCLE repeat_loop
        !
     END IF
     !
     m = m + 1
     !
     v(:,m) = wbest(:) / ( gg(:ngm0) + agg0 )
     !
  END DO repeat_loop
  !
  RETURN
  !
END SUBROUTINE approx_screening2
#else
!----------------------------------------------------------------------------
SUBROUTINE approx_screening2( drho, rhobest )
  !----------------------------------------------------------------------------
  !
  ! ... apply a local-density dependent TF preconditioning to drho
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : e2, pi, tpi, fpi, eps8, eps32
  USE cell_base,            ONLY : omega, tpiba2
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   nl, nlm, gg, ngm
  USE klist,                ONLY : nelec
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  USE control_flags,        ONLY : ngm0
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: drho(ngm,nspin)
  COMPLEX(DP), INTENT(IN)    :: rhobest(ngm,nspin)
  !
  INTEGER, PARAMETER :: mmx = 12
  !
  INTEGER :: &
    iwork(mmx), i, j, m, info
  REAL(DP) :: &
    rs, min_rs, max_rs, avg_rsm1, target, dr2_best
  REAL(DP) :: &
    aa(mmx,mmx), invaa(mmx,mmx), bb(mmx), work(mmx), vec(mmx), agg0
  COMPLEX(DP), ALLOCATABLE :: &
    v(:,:),     &! v(ngm0,mmx)
    w(:,:),     &! w(ngm0,mmx)
    dv(:),      &! dv(ngm0)
    vbest(:),   &! vbest(ngm0)
    wbest(:)     ! wbest(ngm0)
  REAL(DP), ALLOCATABLE :: &
    alpha(:)     ! alpha(nrxx)
  !
  COMPLEX(DP)         :: rrho, rmag
  INTEGER             :: ir, ig
  REAL(DP), PARAMETER :: one_third = 1.D0 / 3.D0
  REAL(DP), EXTERNAL  :: rho_ddot
  !
  !
  IF ( nspin == 2 ) THEN
     !
     DO ig = 1, ngm0
        !
        rrho = drho(ig,1) + drho(ig,2)
        rmag = drho(ig,1) - drho(ig,2)
        !        
        drho(ig,1) = rrho
        drho(ig,2) = rmag
        !
     END DO
     !
  END IF
  !
  target = 0.D0
  !
  IF ( gg(1) < eps8 ) drho(1,1) = ZERO
  !
  ALLOCATE( alpha( nrxx ), v( ngm0, mmx ), &
            w( ngm0, mmx ), dv( ngm0 ), vbest( ngm0 ), wbest( ngm0 ) )
  !
  v(:,:)   = ZERO
  w(:,:)   = ZERO
  dv(:)    = ZERO
  vbest(:) = ZERO
  wbest(:) = ZERO
  !
  ! ... calculate alpha from density
  !
  psic(:) = ZERO
  !
  IF ( nspin == 2 ) THEN
     !
     psic(nl(:ngm0)) = ( rhobest(:ngm0,1) + rhobest(:ngm0,2) )
     !
  ELSE
     !
     psic(nl(:ngm0)) = rhobest(:ngm0,1)
     !
  END IF
  !
  IF ( gamma_only ) psic(nlm(:ngm0)) = CONJG( psic(nl(:ngm0)) )
  !
  CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
  !
  alpha(:) = REAL( psic(:) )
  !
  min_rs   = ( 3.D0 * omega / fpi / nelec )**one_third
  max_rs   = min_rs
  avg_rsm1 = 0.D0
  !
  DO ir = 1, nrxx
     !
     alpha(ir) = ABS( alpha(ir) )
     !
     IF ( alpha(ir) > eps32 ) THEN
        !
        rs        = ( 3.D0 / fpi / alpha(ir) )**one_third
        min_rs    = MIN( min_rs, rs )
        avg_rsm1  = avg_rsm1 + 1.D0 / rs
        max_rs    = MAX( max_rs, rs )
        alpha(ir) = rs
        !
     END IF   
     !
  END DO
  !
  CALL reduce( 1, avg_rsm1 )
  !
  CALL extreme( min_rs, -1 )
  CALL extreme( max_rs, +1 )
  !
  alpha = 3.D0 * ( tpi / 3.D0 )**( 5.D0 / 3.D0 ) * alpha
  !
  avg_rsm1 = ( nr1*nr2*nr3 ) / avg_rsm1
  rs       = ( 3.D0 * omega / fpi / nelec )**one_third
  agg0     = ( 12.D0 / pi )**( 2.D0 / 3.D0 ) / tpiba2 / avg_rsm1
  !
  ! ... calculate deltaV and the first correction vector
  !
  psic(:) = ZERO
  !
  psic(nl(:ngm0)) = drho(:,1)
  !
  IF ( gamma_only ) psic(nlm(:ngm0)) = CONJG( psic(nl(:ngm0)) )
  !
  CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, +1 )
  !
  psic(:) = psic(:) * alpha(:)
  !
  CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
  !
  dv(:) = psic(nl(:ngm0)) * gg(:ngm0) * tpiba2
  v(:,1)= psic(nl(:ngm0)) * gg(:ngm0) / ( gg(:ngm0) + agg0 )
  !
  m       = 1
  aa(:,:) = 0.D0
  bb(:)   = 0.D0
  !
  repeat_loop: DO
     !
     ! ... generate the vector w
     !     
     w(:,m) = fpi * e2 * v(:,m)
     !
     psic(:) = ZERO
     !
     psic(nl(:ngm0)) = v(:,m)
     !
     IF ( gamma_only ) psic(nlm(:ngm0)) = CONJG( psic(nl(:ngm0)) )
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, +1 )
     !
     psic(:) = psic(:) * alpha(:)
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
     !
     w(:,m) = w(:,m) + gg(:ngm0) * tpiba2 * psic(nl(:ngm0))
     !
     ! ... build the linear system
     !
     DO i = 1, m
        !
        aa(i,m) = rho_ddot( w(1,i), w(1,m), ngm0, ngm0, 1, ngm0)
        !
        aa(m,i) = aa(i,m)
        !
     END DO
     !
     bb(m) = rho_ddot( w(1,m), dv, ngm0, ngm0, 1, ngm0 )
     !
     ! ... solve it -> vec
     !
     invaa = aa
     !
     CALL DSYTRF( 'U', m, invaa, mmx, iwork, work, mmx, info )
     CALL errore( 'broyden', 'factorization', info )
     !
     CALL DSYTRI( 'U', m, invaa, mmx, iwork, work, info )
     CALL errore( 'broyden', 'DSYTRI', info )
     !     
     FORALL( i = 1:m, j = 1:m, j > i ) invaa(j,i) = invaa(i,j)
     !
     FORALL( i = 1:m ) vec(i) = SUM( invaa(i,:)*bb(:) )
     !
     vbest(:) = ZERO
     wbest(:) = dv(:)
     !
     DO i = 1, m
        !
        vbest = vbest + vec(i) * v(:,i)
        wbest = wbest - vec(i) * w(:,i)
        !
     END DO
     !
     dr2_best = rho_ddot( wbest, wbest, ngm0, ngm0, 1, ngm0 )
     !
     IF ( target == 0.D0 ) target = 1.D-6 * dr2_best
     !
     IF ( dr2_best < target ) THEN
        !
        drho(:ngm0,1) = vbest(:)
        !
        IF ( nspin == 2 ) THEN
           !
           DO ig = 1, ngm0
              !
              rrho = drho(ig,1)
              rmag = drho(ig,2)
              !
              drho(ig,1) = 0.5D0 * ( rrho + rmag )
              drho(ig,2) = 0.5D0 * ( rrho - rmag )
              !
           END DO
           !
        END IF
        !
        DEALLOCATE( alpha, v, w, dv, vbest, wbest )
        !
        EXIT repeat_loop
        !
     ELSE IF ( m >= mmx ) THEN
        !
        m = 1
        !
        v(:,m)  = vbest(:)
        aa(:,:) = 0.D0
        bb(:)   = 0.D0
        !
        CYCLE repeat_loop
        !
     END IF
     !
     m = m + 1
     !
     v(:,m) = wbest(:) / ( gg(:ngm0) + agg0 )
     !
  END DO repeat_loop
  !
  RETURN
  !
END SUBROUTINE approx_screening2
#endif
