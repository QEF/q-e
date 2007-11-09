!
! Copyright (C) 2002-2006 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#define ZERO ( 0._dp, 0._dp )
!
!----------------------------------------------------------------------------
SUBROUTINE mix_rho( input_rhout, rhoin, &
                    input_becout, becin, input_nsout, nsin,       &
                    alphamix, dr2, tr2_min, iter, n_iter, conv )
  !----------------------------------------------------------------------------
  !
  ! ... Modified Broyden's method for charge density mixing
  ! ...         D.D. Johnson PRB 38, 12807 (1988)
  !
  ! ... On output: the mixed density is in rhoin, mixed augmentation
  ! ...            channel occ. is in becin
  !                input_rhocout, input_becout etc are unchanged
  !
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat
  USE gvect,          ONLY : nrxx, ngm
  USE gsmooth,        ONLY : ngms
  USE ldaU,           ONLY : lda_plus_u, Hubbard_lmax
  USE funct,          ONLY : dft_is_meta
  USE lsda_mod,       ONLY : nspin
  USE control_flags,  ONLY : imix, ngm0, tr2, io_level
  USE io_files,       ONLY : find_free_unit
  USE cell_base,      ONLY : omega
  ! ... for PAW:
  USE uspp_param,           ONLY : nhm
  USE grid_paw_variables,   ONLY : okpaw
  USE rad_paw_routines,     ONLY : paw_ddot
  USE scf,            ONLY : scf_type, create_scf_type, destroy_scf_type, &
                             mix_type, create_mix_type, destroy_mix_type, &
                             assign_scf_to_mix_type, assign_mix_to_scf_type, &
                             mix_type_AXPY, diropn_mix_file, close_mix_file, &
                             davcio_mix_type, rho_ddot, high_frequency_mixing
   USE io_global,     ONLY : stdout
  !
  IMPLICIT NONE
  integer :: kilobytes
  !
  ! ... First the I/O variable
  !
  INTEGER :: &
    iter,        &!  (in) counter of the number of iterations
    n_iter        !  (in) numb. of iterations used in mixing
  REAL(DP) :: &
    alphamix ,   &! (in) mixing factor
    dr2           ! (out) the estimated errr on the energy
  REAL(DP) :: &
    tr2_min       ! estimated error in diagonalization. If the estimated
                  ! scf error is smaller than this, exit: a more accurate 
                  ! diagonalization is needed
  LOGICAL :: &
    conv          ! (out) .true. if the convergence has been reached

  type(scf_type), intent(in)    :: input_rhout
  REAL(DP),    intent(in)    :: input_becout(nhm*(nhm+1)/2,nat,nspin), &! PAW
             input_nsout(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)   ! LDA+U
  type(scf_type), intent(inout) :: rhoin
  REAL(DP),    intent(inout) :: becin (nhm*(nhm+1)/2,nat,nspin),       &!PAW
                   nsin(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)    ! LDA+U
  !
  ! ... Here the local variables
  !
  type(mix_type) :: rhout, rhoin_m
  REAL(DP),    allocatable   :: becout(:,:,:), &! PAW
                                nsout(:,:,:,:)  ! LDA+U
  INTEGER, PARAMETER :: &
    maxmix = 25             ! max number of iterations for charge mixing
  INTEGER ::    &
    iunmix,        &! I/O unit number of charge density file in G-space
    iunmix_ns,     &! I/O unit number of ns file
    iunmix_paw,    &! I/O unit number of PAW file
    iter_used,     &! actual number of iterations used
    ipos,          &! index of the present iteration
    inext,         &! index of the next iteration
    i, j,          &! counters on number of iterations
    iwork(maxmix), &! dummy array used as output by libr. routines
    info,          &! flag saying if the exec. of libr. routines was ok
    ldim            ! 2 * Hubbard_lmax + 1
  type(mix_type) :: rhoin_save, rhout_save
  COMPLEX(DP), ALLOCATABLE :: &
    nsinsave(:,:,:,:),  &!
    nsoutsave(:,:,:,:)   !
  REAL(DP), ALLOCATABLE :: &
    becinsave(:,:,:),   &
    becoutsave(:,:,:)
  REAL(DP) :: &
    betamix(maxmix,maxmix), &
    gamma0,                 &
    work(maxmix)
  LOGICAL :: &
    savetofile,  &! save intermediate steps on file $prefix."mix",...
    exst          ! if true the file exists
  !
  ! ... saved variables and arrays
  !
  INTEGER, SAVE :: &
    mixrho_iter = 0    ! history of mixing
  
  TYPE(mix_type), ALLOCATABLE, SAVE :: &
    df(:),        &! information from preceding iterations
    dv(:)          !     "  "       "     "        "  "
  REAL(DP), ALLOCATABLE, SAVE :: &
    df_ns(:,:,:,:,:), &! idem 
    dv_ns(:,:,:,:,:)   ! idem
  REAL(DP), ALLOCATABLE, SAVE :: &
    df_bec(:,:,:,:), &! idem !PAW
    dv_bec(:,:,:,:)   ! idem !PAW
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: ns_ddot
#ifdef __GRID_PAW
  REAL(DP), EXTERNAL :: rho1_ddot
#endif
  !
  CALL start_clock( 'mix_rho' )
  !
  ngm0 = ngms
  !
  mixrho_iter = iter
  !
  IF ( n_iter > maxmix ) CALL errore( 'mix_rho', 'n_iter too big', 1 )
  !
  IF ( lda_plus_u ) ldim = 2 * Hubbard_lmax + 1
  !
  savetofile = (io_level > 1)
  !
  ! define rhocout variables and copy input_rhocout in there
  !
  call create_mix_type(rhout)
  call assign_scf_to_mix_type(input_rhout, rhout)
  if (lda_plus_u) allocate (nsout(ldim,ldim,nspin,nat))
  if (okpaw) allocate (becout(nhm*(nhm+1)/2,nat,nspin))
  !
  call create_mix_type(rhoin_m)
  call assign_scf_to_mix_type(rhoin, rhoin_m)
  call mix_type_AXPY ( -1.d0, rhoin_m, rhout )
  IF ( lda_plus_u ) nsout(:,:,:,:) = input_nsout(:,:,:,:) - nsin(:,:,:,:)
  IF ( okpaw )      becout(:,:,:)  = input_becout(:,:,:)  - becin(:,:,:) !PAW
 
  !
  dr2 = rho_ddot( rhout, rhout, ngms )  !!!! this used to be ngm NOT ngms
  IF ( lda_plus_u ) dr2 = dr2 + ns_ddot( nsout, nsout, nspin )
  IF ( okpaw )      dr2 = dr2 + PAW_ddot ( becout, becout ) !PAW
  !
  conv = ( dr2 < tr2 )
  !
  IF ( conv .OR. dr2 < tr2_min ) THEN
     !
     ! ... if convergence is achieved or if the self-consistency error (dr2) is
     ! ... smaller than the estimated error due to diagonalization (tr2_min),
     ! ... exit and leave rhoin and rhocout unchanged
     !
     IF ( ALLOCATED( df ) ) THEN
         DO i=1, n_iter
            call destroy_mix_type(df(i))
         END DO
         DEALLOCATE( df )
     END IF
     IF ( ALLOCATED( dv ) ) THEN
         DO i=1, n_iter
            call destroy_mix_type(dv(i))
         END DO
         DEALLOCATE( dv )
     END IF
     IF ( ALLOCATED( df_ns ) )  DEALLOCATE( df_ns )
     IF ( ALLOCATED( dv_ns ) )  DEALLOCATE( dv_ns )
     IF ( ALLOCATED( df_bec ) ) DEALLOCATE ( df_bec ) !PAW
     IF ( ALLOCATED( dv_bec ) ) DEALLOCATE ( dv_bec ) !PAW
     !
     call destroy_mix_type(rhoin_m)
     call destroy_mix_type(rhout)
     if (lda_plus_u) deallocate (nsout)
     if (okpaw) deallocate (becout)
 
     CALL stop_clock( 'mix_rho' )
     !
     RETURN
     !
  END IF
  !
  IF ( savetofile ) THEN
     !
     iunmix = find_free_unit()
     CALL diropn_mix_file( iunmix, 'mix', exst )
     !
     IF ( lda_plus_u ) THEN
        iunmix_ns = find_free_unit()
        CALL diropn( iunmix_ns, 'mix.ns', ldim*ldim*nspin*nat, exst )
     END IF
     !
     IF ( okpaw ) then
        iunmix_paw = find_free_unit()
        CALL diropn( iunmix_paw, 'mix.bec', &
                     ( (nhm * (nhm + 1)/2) * nat * nspin ), exst )
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
     IF ( .NOT. ALLOCATED( df ) ) THEN
        ALLOCATE( df( n_iter ) )
        DO i=1,n_iter
           CALL create_mix_type( df(i) )
        END DO
     END IF
     IF ( .NOT. ALLOCATED( dv ) ) THEN
        ALLOCATE( dv( n_iter ) )
        DO i=1,n_iter
           CALL create_mix_type( dv(i) )
        END DO
     END IF
     !
     IF ( lda_plus_u ) THEN
        IF ( .NOT. ALLOCATED( df_ns ) ) &
           ALLOCATE( df_ns( ldim, ldim, nspin, nat, n_iter ) )
        IF ( .NOT. ALLOCATED( dv_ns ) ) &
           ALLOCATE( dv_ns( ldim, ldim, nspin, nat, n_iter ) )
     END IF
     IF ( okpaw ) THEN
        IF ( .NOT. ALLOCATED( df_bec ) ) &
             ALLOCATE( df_bec( nhm * (nhm + 1)/2, nat, nspin, n_iter ) )
        IF ( .NOT. ALLOCATED( dv_bec ) ) &
             ALLOCATE( dv_bec( nhm * (nhm + 1)/2, nat, nspin, n_iter ) )
        !
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
        CALL davcio_mix_type( df(ipos), iunmix, 1, -1 )
        CALL davcio_mix_type( dv(ipos), iunmix, 2, -1 )
        IF ( lda_plus_u ) THEN
           CALL davcio( df_ns(1,1,1,1,ipos),ldim*ldim*nspin*nat,iunmix_ns,1,-1 )
           CALL davcio( dv_ns(1,1,1,1,ipos),ldim*ldim*nspin*nat,iunmix_ns,2,-1 )
        END IF
        IF ( okpaw ) THEN
           CALL davcio( df_bec(1,1,1,ipos),(nhm*(nhm+1)/2)*nat*nspin,iunmix_paw,1,-1 )
           CALL davcio( dv_bec(1,1,1,ipos),(nhm*(nhm+1)/2)*nat*nspin,iunmix_paw,2,-1 )
        END IF
        !
     END IF
     !
     call mix_type_AXPY ( -1.d0, rhout, df(ipos) )
     call mix_type_AXPY ( -1.d0, rhoin_m, dv(ipos) )
     IF ( lda_plus_u ) THEN
        df_ns(:,:,:,:,ipos) = df_ns(:,:,:,:,ipos) - nsout
        dv_ns(:,:,:,:,ipos) = dv_ns(:,:,:,:,ipos) - nsin
     END IF
     IF ( okpaw ) THEN
        df_bec(:,:,:,ipos) = df_bec(:,:,:,ipos) - becout
        dv_bec(:,:,:,ipos) = dv_bec(:,:,:,ipos) - becin
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
           CALL davcio_mix_type( df(i), iunmix, 2*i+1, -1 )
           CALL davcio_mix_type( dv(i), iunmix, 2*i+2, -1 )
           IF ( lda_plus_u ) THEN
              CALL davcio(df_ns(1,1,1,1,i),ldim*ldim*nspin*nat,iunmix_ns,2*i+1,-1)
              CALL davcio(dv_ns(1,1,1,1,i),ldim*ldim*nspin*nat,iunmix_ns,2*i+2,-1)
           END IF
           IF ( okpaw ) THEN
              CALL davcio(df_bec(1,1,1,i),(nhm*(nhm+1)/2)*nat*nspin,iunmix_paw,2*i+1,-1)
              CALL davcio(dv_bec(1,1,1,i),(nhm*(nhm+1)/2)*nat*nspin,iunmix_paw,2*i+2,-1)
           END IF
           !
        END IF
        !
     END DO
     CALL davcio_mix_type( rhout,   iunmix, 1, 1 )
     CALL davcio_mix_type( rhoin_m, iunmix, 2, 1 )
     !
     IF ( mixrho_iter > 1 ) THEN
        CALL davcio_mix_type( df(ipos), iunmix, 2*ipos+1, 1 )
        CALL davcio_mix_type( dv(ipos), iunmix, 2*ipos+2, 1 )
     END IF
     !
     IF ( lda_plus_u ) THEN
        !
        CALL davcio( nsout, ldim*ldim*nspin*nat, iunmix_ns, 1, 1 )
        CALL davcio( nsin , ldim*ldim*nspin*nat, iunmix_ns, 2, 1 )
        !
        IF ( mixrho_iter > 1 ) THEN
           !
           CALL davcio( df_ns(1,1,1,1,ipos), &
                        ldim*ldim*nspin*nat, iunmix_ns, 2*ipos+1, 1 )
           CALL davcio( dv_ns(1,1,1,1,ipos), &
                        ldim*ldim*nspin*nat, iunmix_ns, 2*ipos+2, 1 )
        END IF
        !
     END IF
     !
     IF ( okpaw ) THEN
        !
        CALL davcio( becout, (nhm*(nhm+1)/2)*nat*nspin, iunmix_paw, 1, 1 )
        CALL davcio( becin , (nhm*(nhm+1)/2)*nat*nspin, iunmix_paw, 2, 1 )
        !
        IF ( mixrho_iter > 1 ) THEN
           !
           CALL davcio( df_bec(1,1,1,ipos), (nhm*(nhm+1)/2)*nat*nspin, &
                        iunmix_paw, 2*ipos+1, 1 )
           CALL davcio( dv_bec(1,1,1,ipos), (nhm*(nhm+1)/2)*nat*nspin, &
                        iunmix_paw, 2*ipos+2, 1 )
        END IF
        !
     END IF
     !
  ELSE
     !
     call create_mix_type (rhoin_save)
     call create_mix_type (rhout_save)
     !
     rhoin_save = rhoin_m
     rhout_save = rhout
     !
     IF ( lda_plus_u ) THEN
        ALLOCATE( nsinsave(  ldim, ldim, nspin, nat ), &
                  nsoutsave( ldim, ldim, nspin, nat ) )
        nsinsave  = nsin
        nsoutsave = nsout
     END IF
     !
     IF ( okpaw ) THEN
        ALLOCATE( becinsave (nhm*(nhm+1)/2,nat,nspin), &
                  becoutsave(nhm*(nhm+1)/2,nat,nspin) )
        becinsave  = becin
        becoutsave = becout
     END IF
     !
  END IF
  !
  DO i = 1, iter_used
     !
     DO j = i, iter_used
        !
        betamix(i,j) = rho_ddot( df(j), df(i), ngm0 )
        !
        IF ( lda_plus_u ) betamix(i,j) = betamix(i,j) + &
                          ns_ddot( df_ns(1,1,1,1,j), df_ns(1,1,1,1,i), nspin )
        !
        IF ( okpaw ) &
           betamix(i,j) = betamix(i,j) + &
                          PAW_ddot( df_bec(1,1,1,j), df_bec(1,1,1,i) )
        !
        betamix(j,i) = betamix(i,j) !symmetrize
        !
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
     work(i) = rho_ddot( df(i), rhout, ngm0 )
     !
     IF ( lda_plus_u ) &
        work(i) = work(i) + ns_ddot( df_ns(1,1,1,1,i), nsout, nspin )
     !
     IF ( okpaw ) &
        work(i) = work(i) + PAW_ddot( df_bec(1,1,1,i), becout )
     !
  END DO
  !
  DO i = 1, iter_used
     !
     gamma0 = SUM( betamix(1:iter_used,i)*work(1:iter_used) )
     !
     call mix_type_AXPY ( -gamma0, dv(i), rhoin_m )
     call mix_type_AXPY ( -gamma0, df(i), rhout )
     !
     IF ( lda_plus_u ) THEN
        nsin  = nsin  - gamma0*dv_ns(:,:,:,:,i)
        nsout = nsout - gamma0*df_ns(:,:,:,:,i)
     END IF
     !
     IF ( okpaw ) THEN
        becin  = becin  - gamma0 * dv_bec(:,:,:,i)
        becout = becout - gamma0 * df_bec(:,:,:,i)
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
        CLOSE( iunmix_ns, STATUS = 'KEEP' )
        !
        DEALLOCATE( df_ns, dv_ns )
        !
     END IF
     !
     IF ( okpaw ) THEN
        !
        CLOSE( iunmix_paw, STATUS = 'KEEP' )
        !
        DEALLOCATE( df_bec, dv_bec )
        !
     END IF
     !
     call close_mix_file( iunmix )
     !
     IF ( ALLOCATED( df ) ) THEN
         DO i=1, n_iter
            call destroy_mix_type(df(i))
         END DO
         DEALLOCATE( df )
     END IF
     IF ( ALLOCATED( dv ) ) THEN
         DO i=1, n_iter
            call destroy_mix_type(dv(i))
         END DO
         DEALLOCATE( dv )
     END IF
     !
  ELSE
     !
     inext = mixrho_iter - ( ( mixrho_iter - 1 ) / n_iter ) * n_iter
     !
     IF ( lda_plus_u ) THEN
        df_ns(:,:,:,:,inext) = nsoutsave
        dv_ns(:,:,:,:,inext) = nsinsave
        DEALLOCATE( nsinsave, nsoutsave )
     END IF
     !
     IF ( okpaw ) THEN
        df_bec(:,:,:,inext) = becoutsave
        dv_bec(:,:,:,inext) = becinsave
        DEALLOCATE( becinsave, becoutsave )
     END IF
     !
     df(inext) = rhout_save
     dv(inext) = rhoin_save
     !
     call destroy_mix_type( rhoin_save )
     call destroy_mix_type( rhout_save )
     !
  END IF
  !
  ! ... preconditioning the new search direction
  !
  IF ( imix == 1 ) THEN
     !
     CALL approx_screening( rhout )
     !
  ELSE IF ( imix == 2 ) THEN
     !
     CALL approx_screening2( rhout, rhoin_m )
     !
  END IF
  !
  ! ... set new trial density
  !
  call mix_type_AXPY ( alphamix, rhout, rhoin_m )
  IF ( lda_plus_u ) nsin = nsin + alphamix * nsout
  IF ( okpaw ) becin = becin + alphamix * becout
  !
  ! ... simple mixing for high_frequencies (and set to zero the smooth ones)
  call high_frequency_mixing ( rhoin, input_rhout, alphamix )
  ! ... add the mixed rho for the smooth frequencies
  call assign_mix_to_scf_type(rhoin_m,rhoin)
  !
  call destroy_mix_type(rhout)
  call destroy_mix_type(rhoin_m)
  if (lda_plus_u) deallocate (nsout)
  if (okpaw) deallocate (becout)

  CALL stop_clock( 'mix_rho' )
  !
  RETURN
  !
END SUBROUTINE mix_rho
!
#ifdef __GRID_PAW
!----------------------------------------------------------------------------
FUNCTION rho1_ddot( bec1, bec2 )
  !----------------------------------------------------------------------------
  !
  ! ... calculates 4pi/G^2*rho1(-G)*rho2(G) = V1_Hartree(-G)*rho2(G),
  ! ... where rho1 and rho2 are 1-center charges (AE and PS)
  ! ... input variables are the augmentation channel occupations
  ! ... used as an estimate of the self-consistency error on the energy
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, tpi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : ngm, nl, nlm, gg, g, gstart
  USE lsda_mod,      ONLY : nspin
  USE wvfct,         ONLY : gamma_only
  !  
  USE grid_paw_variables, ONLY : pp, okpaw, prodp, prodpt, prod0p, prod0pt
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,         ONLY : nh, nhm 
  USE wvfct,              ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  REAL(DP), INTENT(IN) :: &
     bec1(nhm*(nhm+1)/2,nat,nspin), &
     bec2(nhm*(nhm+1)/2,nat,nspin)
  !
  REAL(DP)                :: rho1_ddot
  !
  ! ... and the local variables
  !
  REAL(DP) :: fac   ! a multiplicative factor
  !
  INTEGER :: gi, ig, na, nt, ih, jh, ijh, ijh2, is
  ! counters
  !
  COMPLEX(DP), POINTER :: prodp_(:,:,:), prod0p_(:,:,:)
  INTEGER :: i_what
  REAL(DP):: i_sign
  !
  rho1_ddot = 0._dp
  !
  IF ( .NOT. okpaw ) RETURN
  !
  gi = gstart
  !
  fac = e2 * fpi / tpiba2
  !  
  whattodo: DO i_what=1, 2
     !
     NULLIFY(prodp_,prod0p_)
     IF (i_what==1) THEN
        prodp_ => prodp
        prod0p_ => prod0p
     ELSE IF (i_what==2) THEN
        prodp_ => prodpt
        prod0p_ => prod0pt
     END IF
     i_sign = DBLE(1-2*(i_what-1)) ! = +1 if i_what==1, -1 if i_what==2
     !
     DO ijh = 1, nhm*(nhm+1)/2
        !
        DO ijh2 = 1, nhm*(nhm+1)/2
           !
           DO na = 1, nat
              !
              nt = ityp (na)
              IF ( nspin == 1 ) THEN
                 !
                 rho1_ddot = rho1_ddot + i_sign * fac * &
                 bec1(ijh,na,1) * prodp_(ijh, ijh2, nt) * bec2(ijh2,na,1)
                 !
!!$              gamma_only case not yet implemented
!!$              IF ( gamma_only ) rho1_ddot = 2.D0 * rho1_ddot
                 !
              ELSE IF ( nspin == 2 ) THEN
                 !
                 ! ... first the charge
                 !
                 rho1_ddot = rho1_ddot + i_sign * fac * &
                  (bec1(ijh,na,1)+bec1(ijh,na,2)) * prodp_(ijh,ijh2,nt) * &
                  (bec2(ijh2,na,1)+bec2(ijh2,na,2))
                 ! 
!!$              IF ( gamma_only ) rho1_ddot = 2.D0 * rho1_ddot
                 !
                 ! ... then the magnetization
                 !
                 fac = e2 * fpi / tpi**2  ! lambda = 1 a.u.
                 !
                 ! ... G=0 term
                 !
                 IF ( gstart == 2 ) THEN
                    !
                    rho1_ddot = rho1_ddot + i_sign * fac * &
                    (bec1(ijh,na,1)-bec1(ijh,na,2)) * prod0p_(ijh,ijh2,nt) * &
                    (bec2(ijh2,na,1)-bec2(ijh2,na,2))
                    !
                 END IF
                 !
!!$              IF ( gamma_only ) fac = 2.D0 * fac
                 !
                 rho1_ddot = rho1_ddot + i_sign * fac * &
                   (bec1(ijh,na,1)-bec1(ijh,na,2)) * prodp_(ijh,ijh2,nt) * &
                   (bec2(ijh2,na,1)-bec2(ijh2,na,2))
                 !
!!$           non-collinear case not yet implemented
!!$           ELSE IF ( nspin == 4 ) THEN
!!$              !
!!$              rho1_ddot = rho1_ddot + fac * DBLE((-1)**(i_what-1)) * 
!!$              bec1(ijh,na,1) * prodp_(ijh, ijh2, nt) * bec2(ijh2,na,1)
!!$              !
!!$              IF ( gamma_only ) rho1_ddot = 2.D0 * rho1_ddot
!!$              !
!!$              fac = e2*fpi / (tpi**2)  ! lambda=1 a.u.
!!$              !
!!$              IF ( gstart == 2 ) THEN
!!$                 !
!!$                 rho1_ddot = rho1_ddot + DBLE((-1)**(i_what-1)) * fac * &
!!$                 bec1(ijh,na,2) * prod0p_(ijh,ijh2,nt) * bec2(ijh2,na,2)+ &
!!$                 bec1(ijh,na,3) * prod0p_(ijh,ijh2,nt) * bec2(ijh2,na,3)+ &
!!$                 bec1(ijh,na,4) * prod0p_(ijh,ijh2,nt) * bec2(ijh2,na,4)
!!$                 !
!!$              END IF
!!$              !
!!$              IF ( gamma_only ) fac = 2.D0 * fac
!!$              !
!!$              rho1_ddot = rho1_ddot + DBLE((-1)**(i_what-1)) * fac * &
!!$              ( bec1(ijh,na,2) * prodp_(ijh,ijh2,nt) * bec2(ijh2,na,2)+ &
!!$                bec1(ijh,na,3) * prodp_(ijh,ijh2,nt) * bec2(ijh2,na,3)+ &
!!$                bec1(ijh,na,4) * prodp_(ijh,ijh2,nt) * bec2(ijh2,na,4) )
!!$                 !
              END IF
              !
           END DO
           !
        END DO
        !
     END DO
     !
  END DO whattodo
  !
  rho1_ddot = rho1_ddot * omega * 0.5D0
  !
  CALL reduce( 1, rho1_ddot )
  !
  RETURN
  !
END FUNCTION rho1_ddot
#endif
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
  USE gvect,         ONLY : gg, ngm, &
                            nr1, nr2, nr3, nrx1, nrx2, nrx3, nl, nlm
  USE klist,         ONLY : nelec
  USE lsda_mod,      ONLY : nspin
  USE control_flags, ONLY : ngm0
  USE scf,           ONLY : mix_type
  USE wvfct,         ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  !
  IMPLICIT NONE  
  !
  type (mix_type), intent(INOUT) :: drho ! (in/out)
  !
  REAL(DP) :: rrho, rmag, rs, agg0
  INTEGER  :: ig, is
  !
  rs = ( 3.D0 * omega / fpi / nelec )**( 1.D0 / 3.D0 )
  !
  agg0 = ( 12.D0 / pi )**( 2.D0 / 3.D0 ) / tpiba2 / rs
  !
  IF ( nspin == 1 .OR. nspin == 4 ) THEN
     !
     drho%of_g(:ngm0,1) =  drho%of_g(:ngm0,1) * gg(:ngm0) / (gg(:ngm0)+agg0)
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     DO ig = 1, ngm0
        !
        rrho = ( drho%of_g(ig,1) + drho%of_g(ig,2) ) * gg(ig) / (gg(ig)+agg0)
        rmag = ( drho%of_g(ig,1) - drho%of_g(ig,2) )
        !
        drho%of_g(ig,1) =  0.5D0*( rrho + rmag )
        drho%of_g(ig,2) =  0.5D0*( rrho - rmag )
        !
     END DO
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE approx_screening
!
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
  USE gvect,                ONLY : gg, ngm, &
                                   nr1, nr2, nr3, nrx1, nrx2, nrx3, nl, nlm
  USE wavefunctions_module, ONLY : psic
  USE klist,                ONLY : nelec
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : gamma_only
  USE control_flags,        ONLY : ngm0
  USE scf,                  ONLY : mix_type, local_tf_ddot
  !
  IMPLICIT NONE
  !
  type(mix_type), intent(inout) :: drho
  type(mix_type), intent(in) :: rhobest
  !
  INTEGER, PARAMETER :: mmx = 12
  !
  INTEGER :: &
    iwork(mmx), i, j, m, info, is
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
  !
  !
  IF ( nspin == 2 ) THEN
     !
     DO ig = 1, ngm0
        !
        rrho = drho%of_g(ig,1) + drho%of_g(ig,2)
        rmag = drho%of_g(ig,1) - drho%of_g(ig,2)
        !        
        drho%of_g(ig,1) = rrho
        drho%of_g(ig,2) = rmag
        !
     END DO
     !
  END IF
  !
  target = 0.D0
  !
  IF ( gg(1) < eps8 ) drho%of_g(1,1) = ZERO
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
     psic(nls(:ngm0)) = ( rhobest%of_g(:ngm0,1) + rhobest%of_g(:ngm0,2) )
     !
  ELSE
     !
     psic(nls(:ngm0)) = rhobest%of_g(:ngm0,1)
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
  psic(nls(:ngm0)) = drho%of_g(:ngm0,1)
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
        aa(i,m) = local_tf_ddot( w(1,i), w(1,m), ngm0)
        !
        aa(m,i) = aa(i,m)
        !
     END DO
     !
     bb(m) = local_tf_ddot( w(1,m), dv, ngm0)
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
     dr2_best = local_tf_ddot( wbest, wbest, ngm0 )
     !
     IF ( target == 0.D0 ) target = 1.D-6 * dr2_best
     !
     IF ( dr2_best < target ) THEN
        !
        drho%of_g(:ngm0,1) = vbest(:)
        !
        IF ( nspin == 2 ) THEN
           !
           DO ig = 1, ngm0
              !
              rrho = drho%of_g(ig,1)
              rmag = drho%of_g(ig,2)
              !
              drho%of_g(ig,1) = 0.5D0 * ( rrho + rmag )
              drho%of_g(ig,2) = 0.5D0 * ( rrho - rmag )
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
