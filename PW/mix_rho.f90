!
! Copyright (C) 2002-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
#undef DEBUG
!
!----------------------------------------------------------------------------
SUBROUTINE mix_rho( rhout, rhoin, nsout, nsin, alphamix, dr2, ethr, ethr_min, &
                    iter, n_iter, filename, conv )
  !----------------------------------------------------------------------------
  !
  ! ... Modified Broyden's method for charge density mixing
  ! ...         d.d. johnson prb 38, 12807 (1988)
  ! ... On output: the mixed density is in rhoin, rhout is UNCHANGED
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   nl, nlm
  USE ldaU,                 ONLY : lda_plus_u, Hubbard_lmax
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : imix, ngm0, tr2
  USE wvfct,                ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  !
  IMPLICIT NONE
  !
  ! ... First the I/O variable
  !
  CHARACTER(LEN=42) :: &
    filename                !  (in) I/O filename for mixing history
                            !  if absent everything is kept in memory
  INTEGER ::    &
    iter,                  &!  (in)  counter of the number of iterations
    n_iter                  !  (in)  numb. of iterations used in mixing
  REAL(KIND=DP) :: &
    rhout(nrxx,nspin),     &! (in) the "out" density; (out) rhout-rhoin
    rhoin(nrxx,nspin),     &! (in) the "in" density; (out) the new dens.
    nsout(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat), &!
    nsin(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat),  &!
    alphamix,              &! (in) mixing factor
    dr2                     ! (out) the estimated errr on the energy
  REAL (KIND=DP) :: &
    ethr,        &! actual threshold for diagonalization
    ethr_min      ! minimal threshold for diagonalization
                  ! if in output ethr >= ethr_min a more accurate
                  ! diagonalization is needed
  LOGICAL :: &
    conv                    ! (out) if true the convergence has been reached
  INTEGER, PARAMETER :: &
    maxmix = 25             ! max number of iterations for charge mixing
  !
  ! ... Here the local variables
  !
  INTEGER ::    &
    iunmix,        &! I/O unit number of charge density file
    iunmix2,       &! I/O unit number of ns file
    iunit,         &! counter on I/O unit numbers
    iter_used,     &! actual number of iterations used
    ipos,          &! index of the present iteration
    inext,         &! index of the next iteration
    i, j,          &! counters on number of iterations
    is,            &! counter on spin component
    ig,            &! counter on G-vectors
    iwork(maxmix), &! dummy array used as output by libr. routines
    info            ! flag saying if the exec. of libr. routines was ok
  COMPLEX(KIND=DP), ALLOCATABLE :: &
    rhocin(:,:),        &! rhocin(ngm0,nspin)
    rhocout(:,:),       &! rhocout(ngm0,nspin)
    rhoinsave(:),       &! rhoinsave(ngm0*nspin): work space
    rhoutsave(:),       &! rhoutsave(ngm0*nspin): work space
    nsinsave(:,:,:,:),  &!
    nsoutsave(:,:,:,:)   !
  REAL(KIND=DP), ALLOCATABLE, SAVE :: &
    df_ns(:,:,:,:,:), &! df_ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,n_iter):idem
    dv_ns(:,:,:,:,:)   ! dv_ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat,n_iter):idem            
  COMPLEX(KIND=DP), ALLOCATABLE, SAVE :: &
    df(:,:),   &! df(ngm0*nspin,n_iter): information from preceding iterations
    dv(:,:)     ! dv(ngm0*nspin,n_iter):    "  "       "     "        "  "
  INTEGER :: &
    ldim
  REAL(KIND=DP) :: &
    betamix(maxmix,maxmix), &
    gamma0,                 &
    work(maxmix),           &
    dehar
  LOGICAL :: &
    saveonfile,   &! save intermediate steps on file "filename"
    opnd,         &! if true the file is already opened
    exst           ! if true the file exists
  REAL(KIND=DP), EXTERNAL :: rho_dot_product, ns_dot_product, fn_dehar
  !
  !
  CALL start_clock( 'mix_rho' )
  !
  IF ( iter < 1 ) CALL errore( 'mix_rho', 'iter is wrong', 1 )
  !
  IF ( n_iter > maxmix ) CALL errore( 'mix_rho', 'n_iter too big', 1 )
  !
  IF ( lda_plus_u ) ldim = 2 * Hubbard_lmax + 1
  !
  saveonfile = ( filename /= ' ' )
  !
  !CALL DAXPY( nrxx*nspin, -1.D0, rhoin, 1, rhout, 1 )
  !
  ALLOCATE( rhocin(ngm0,nspin), rhocout(ngm0,nspin) )
  !
  ! ... psic is used as work space - must be already allocated !
  !
  DO is = 1, nspin
     !
     psic(:) = CMPLX( rhoin(:,is), 0.D0 )
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1 )
     !
     rhocin(:,is) = psic(nl(:))
     !
     psic(:) = CMPLX( rhout(:,is), 0.D0 )
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1 )
     !
     rhocout(:,is) = psic(nl(:)) - rhocin(:,is)
     !
  END DO
  !
  IF ( lda_plus_u ) nsout(:,:,:,:) = nsout(:,:,:,:) - nsin(:,:,:,:)
  !
  dr2 = rho_dot_product( rhocout, rhocout ) + ns_dot_product( nsout, nsout )
  !
  ! ... ethr_min is dr2 * 1.D0 / nelec or - 1.0 if no check is required
  !
  IF ( ethr_min > 0.D0 ) THEN
     !
     ethr_min = dr2 * ethr_min
     !
     IF ( ethr > ethr_min ) THEN
        !
        ethr = ethr_min
        !
        ! ... rhoin and rhout are unchanged
        !
        DEALLOCATE( rhocin, rhocout )
        !
        CALL stop_clock( 'mix_rho' )
        !
        RETURN
        !
     END IF
     !   
  END IF
  ! 
  conv = ( dr2 < tr2 )
  !
  dehar = fn_dehar( rhocout )
  !
#if defined (DEBUG)
  !IF ( lda_plus_u ) WRITE( stdout,*) ' ns_dr2 =', ns_dot_product(nsout,nsout)
  IF ( conv ) THEN
     !
     WRITE( stdout, 100 ) &
         dr2, rho_dot_product(rhocout,rhocout) + ns_dot_product(nsout,nsout)
     WRITE( stdout, '(" dehar =",F15.8)' ) dehar
     !
  END IF
#endif
  !
  IF ( saveonfile ) THEN
     !
     DO iunit = 99, 1, - 1
        !
        INQUIRE( UNIT = iunit, OPENED = opnd )
        iunmix = iunit
        IF ( .NOT. opnd ) GO TO 10
        !
     END DO
     !
     CALL errore( 'mix_rho', 'free unit not found?!?', 1 )
     !
10   CONTINUE
     !
     IF ( lda_plus_u ) THEN
        !
        DO iunit = ( iunmix - 1 ), 1, - 1
           !
           INQUIRE( UNIT = iunit, OPENED = opnd )
           iunmix2 = iunit
           IF ( .NOT. opnd ) GO TO 20
           !
        END DO
        CALL errore( 'mix_rho', 'second free unit not found?!?', 1 )
        !
20      CONTINUE
        !
     END IF
     !
     IF ( conv ) THEN
        !
        CALL diropn( iunmix, filename, ( 2 * ngm0 * nspin ), exst )
        CLOSE( UNIT = iunmix, STATUS = 'DELETE' )
        !
        IF ( lda_plus_u ) THEN
           !
           CALL diropn( iunmix2, TRIM( filename ) // '.ns', &
                        ( ldim * ldim * nspin * nat ), exst )
           CLOSE( UNIT = iunmix2, STATUS = 'DELETE' )
           !
        END IF
        !
        DEALLOCATE( rhocin, rhocout )
        !
        CALL stop_clock( 'mix_rho' )
        !
        RETURN
        !
     END IF
     !
     CALL diropn(iunmix,filename,2*ngm0*nspin,exst)
     !
     IF ( lda_plus_u ) &
        CALL diropn( iunmix2, TRIM( filename ) // '.ns', &
                     ( ldim * ldim * nspin * nat ), exst )
     !
     IF ( iter > 1 .AND. .NOT. exst ) THEN
        !
        CALL errore( 'mix_rho','file not found, restarting', - 1 )
        iter = 1
        !
     END IF
     !
     ALLOCATE( df(ngm0*nspin,n_iter), dv(ngm0*nspin,n_iter) )
     !
     IF ( lda_plus_u ) &
        ALLOCATE( df_ns(ldim,ldim,nspin,nat,n_iter), &
                  dv_ns(ldim,ldim,nspin,nat,n_iter) )
     !
  ELSE
     IF ( iter == 1 ) THEN
        !
        ALLOCATE( df(ngm0*nspin,n_iter), dv(ngm0*nspin,n_iter) )
        !
        IF ( lda_plus_u ) &
           ALLOCATE( df_ns(ldim,ldim,nspin,nat,n_iter),&
                     dv_ns(ldim,ldim,nspin,nat,n_iter) )
        !
     END IF
     !
     IF ( conv ) THEN
        !
        IF ( lda_plus_u ) DEALLOCATE( df_ns, dv_ns )
        !
        DEALLOCATE( df, dv )
        DEALLOCATE( rhocin, rhocout )
        !
        CALL stop_clock( 'mix_rho' )
        !
        RETURN
        !
     END IF
     !
     ALLOCATE( rhoinsave(ngm0*nspin), rhoutsave(ngm0*nspin) )
     !
     IF ( lda_plus_u ) &
        ALLOCATE( nsinsave (ldim,ldim,nspin,nat), &
                  nsoutsave(ldim,ldim,nspin,nat) )
     !
  END IF
  !
  ! ... copy only the high frequency Fourier component into rhoin
  ! ...                                            (NB: rhout=rhout-rhoin)
  !
  CALL DCOPY( nrxx * nspin, rhout, 1, rhoin, 1 )
  !
  DO is = 1, nspin
     !
     psic(:) = ( 0.D0, 0.D0 )
     !
     psic(nl(:)) = rhocin(:,is)+rhocout(:,is)
     !
     IF ( gamma_only ) &
        psic(nlm(:)) = CONJG( psic(nl(:)) )
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
     CALL DAXPY( nrxx, -1.D0, psic, 2, rhoin(1,is), 1 )
     !
  END DO
  !
  ! ... iter_used = iter-1  if iter <= n_iter
  ! ... iter_used = n_iter  if iter >  n_iter
  !
  iter_used = MIN( ( iter - 1 ), n_iter )
  !
  ! ... ipos is the position in which results from the present iteration
  ! ... are stored. ipos=iter-1 until ipos=n_iter, then back to 1,2,...
  !
  ipos = iter - 1 - ( ( iter - 2 ) / n_iter ) * n_iter
  !
  IF ( iter > 1 ) THEN
     !
     IF ( saveonfile ) THEN
        !
        CALL davcio( df(1,ipos), 2*ngm0*nspin, iunmix, 1, -1 )
        CALL davcio( dv(1,ipos), 2*ngm0*nspin, iunmix, 2, -1 )
        !
        IF ( lda_plus_u ) THEN
           !
           CALL davcio(df_ns(1,1,1,1,ipos),ldim*ldim*nspin*nat,iunmix2,1,-1)
           CALL davcio(dv_ns(1,1,1,1,ipos),ldim*ldim*nspin*nat,iunmix2,2,-1)
           !
        END IF
        !
     END IF
     !
     CALL DAXPY( 2*ngm0*nspin, -1.D0, rhocout, 1, df(1,ipos), 1 )
     CALL DAXPY( 2*ngm0*nspin, -1.D0, rhocin , 1, dv(1,ipos), 1 )
     !norm = sqrt(rho_dot_product(df(1,ipos),df(1,ipos)) + &
     !            ns_dot_product(df_ns(1,1,1,1,ipos),df_ns(1,1,1,1,ipos)) )
     !call DSCAL (2*ngm0*nspin,-1.d0/norm,df(1,ipos),1)
     !call DSCAL (2*ngm0*nspin,-1.d0/norm,dv(1,ipos),1)
     !
     IF ( lda_plus_u ) THEN
        !
        CALL DAXPY(ldim*ldim*nspin*nat,-1.D0,nsout,1,df_ns(1,1,1,1,ipos),1)
        CALL DAXPY(ldim*ldim*nspin*nat,-1.D0,nsin ,1,dv_ns(1,1,1,1,ipos),1)
        !
     END IF
     !
  END IF
  !
  IF ( saveonfile ) THEN
     !
     DO i = 1, iter_used
        !
        IF ( i /= ipos ) THEN
           !
           CALL davcio( df(1,i), 2*ngm0*nspin, iunmix, 2*i+1, -1 )
           CALL davcio( dv(1,i), 2*ngm0*nspin, iunmix, 2*i+2, -1 )
           !
           IF ( lda_plus_u ) THEN
              !
              CALL davcio(df_ns(1,1,1,1,i),ldim*ldim*nspin*nat,iunmix2,2*i+1,-1)
              CALL davcio(dv_ns(1,1,1,1,i),ldim*ldim*nspin*nat,iunmix2,2*i+2,-1)
              !
           END IF
           !
        END IF
        !
     END DO
     !
     CALL davcio( rhocout, 2*ngm0*nspin, iunmix, 1, 1 )
     CALL davcio( rhocin , 2*ngm0*nspin, iunmix, 2, 1 )
     !
     IF ( iter > 1 ) THEN
        !
        CALL davcio( df(1,ipos), 2*ngm0*nspin, iunmix, 2*ipos+1, 1 )
        CALL davcio( dv(1,ipos), 2*ngm0*nspin, iunmix, 2*ipos+2, 1 )
        !
     END IF
     !
     IF ( lda_plus_u ) THEN
        !
        CALL davcio( nsout, ldim*ldim*nspin*nat, iunmix2, 1, 1 )
        CALL davcio( nsin , ldim*ldim*nspin*nat, iunmix2, 2, 1 )
        !
        IF ( iter > 1 ) THEN
           !
           CALL davcio( df_ns(1,1,1,1,ipos), ldim*ldim*nspin*nat, &
                        iunmix2, 2*ipos+1, 1 )
           CALL davcio( dv_ns(1,1,1,1,ipos), ldim*ldim*nspin*nat, &
                        iunmix2, 2*ipos+2, 1 )
        END IF
        !
     END IF
     !
  ELSE
     !
     CALL DCOPY( 2*ngm0*nspin, rhocin , 1, rhoinsave, 1 )
     CALL DCOPY( 2*ngm0*nspin, rhocout, 1, rhoutsave, 1 )
     !
     IF ( lda_plus_u ) THEN
        !
        CALL DCOPY( ldim*ldim*nspin*nat, nsin , 1, nsinsave , 1 )
        CALL DCOPY( ldim*ldim*nspin*nat, nsout, 1, nsoutsave, 1 )
        !
     END IF
     !
  END IF
  !
  DO i = 1, iter_used
     DO j = i, iter_used
        !
        betamix(i,j) = rho_dot_product( df(1,j), df(1,i) ) 
        !
        IF ( lda_plus_u ) &
           betamix(i,j) = betamix(i,j) + &
                          ns_dot_product( df_ns(1,1,1,1,j), df_ns(1,1,1,1,i) )
        !
     END DO
  END DO
  !
  CALL DSYTRF( 'U', iter_used, betamix, maxmix, iwork, work, maxmix, info )
  CALL errore( 'broyden', 'factorization', info )
  CALL DSYTRI( 'U', iter_used, betamix, maxmix, iwork, work, info )
  CALL errore( 'broyden', 'DSYTRI', info )
  !
  DO i = 1, iter_used
     DO j = ( i + 1 ), iter_used
        !
        betamix(j,i) = betamix(i,j)
        !
     END DO
  END DO
  !
  DO i = 1, iter_used
     !
     work(i) = rho_dot_product( df(1,i), rhocout )
     !
     IF ( lda_plus_u ) &
        work(i) = work(i) + ns_dot_product( df_ns(1,1,1,1,i), nsout )
     !
  END DO
  !
  DO i = 1, iter_used
     !
     gamma0 = 0.D0
     !
     DO j = 1, iter_used
        gamma0 = gamma0 + betamix(j,i) * work(j)
     END DO
     !
     CALL DAXPY( 2*ngm0*nspin, -gamma0, dv(1,i), 1, rhocin, 1 )
     CALL DAXPY( 2*ngm0*nspin, -gamma0, df(1,i), 1, rhocout, 1 )
     !
     IF ( lda_plus_u ) THEN
        CALL DAXPY(ldim*ldim*nspin*nat,-gamma0,dv_ns(1,1,1,1,i),1,nsin(1,1,1,1) ,1)
        CALL DAXPY(ldim*ldim*nspin*nat,-gamma0,df_ns(1,1,1,1,i),1,nsout(1,1,1,1),1)
     END IF
     !
  END DO
  !
#if defined (DEBUG)
  WRITE( stdout, 100 ) &
      dr2, rho_dot_product( rhocout, rhocout ) + ns_dot_product( nsout, nsout )
  WRITE( stdout, '(" dehar =",F15.8)' ) dehar
#endif
  !
  ! ... auxiliary vectors dv and df not needed anymore
  !
  IF ( saveonfile ) THEN
     !
     IF ( lda_plus_u ) THEN
        CLOSE( iunmix2, STATUS = 'KEEP' )
        DEALLOCATE( df_ns, dv_ns )
     END IF
     CLOSE( iunmix, STATUS = 'KEEP' )
     DEALLOCATE( df, dv )
     !
  ELSE
     !
     inext = iter - ( ( iter - 1 ) / n_iter ) * n_iter
     !
     IF ( lda_plus_u ) THEN
        CALL DCOPY(ldim*ldim*nspin*nat,nsoutsave,1,df_ns(1,1,1,1,inext),1)
        CALL DCOPY(ldim*ldim*nspin*nat,nsinsave ,1,dv_ns(1,1,1,1,inext),1)
        DEALLOCATE( nsinsave, nsoutsave )
     END IF
     !
     CALL DCOPY( 2*ngm0*nspin, rhoutsave, 1, df(1,inext), 1 )
     CALL DCOPY( 2*ngm0*nspin, rhoinsave, 1, dv(1,inext), 1 )
     DEALLOCATE( rhoinsave, rhoutsave )
     !
  END IF
  !
  ! ... preconditioning the new search direction (if imix.gt.0)
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
  CALL DAXPY( 2*ngm0*nspin, alphamix, rhocout, 1, rhocin, 1 )
  !
  DO is = 1, nspin
     !
     psic(:) = ( 0.D0, 0.D0 )
     !
     psic(nl(:)) = rhocin(:,is)
     !
     IF ( gamma_only ) &
        psic(nlm(:)) = conjg ( psic(nl(:)) )
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
     CALL DAXPY( nrxx, 1.D0, psic, 2, rhoin(1,is), 1 )
     !
  END DO
  !
  IF ( lda_plus_u ) & 
     CALL DAXPY( ldim*ldim*nspin*nat, alphamix, nsout, 1, nsin, 1 )
  !
  ! ... clean up
  !
  DEALLOCATE( rhocout )
  DEALLOCATE( rhocin )
  !
  CALL stop_clock( 'mix_rho' )
  !
  RETURN
  !
100  FORMAT( ' dr2 =', 1PE15.1, ' internal_best_dr2= ', 1PE15.1 )
  !
END SUBROUTINE mix_rho
!
!
!----------------------------------------------------------------------------
FUNCTION rho_dot_product( rho1, rho2 )
  !----------------------------------------------------------------------------
  !
  ! ... this function evaluates the dot product between two input densities
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, tpi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gstart, gg
  USE lsda_mod,      ONLY : nspin
  USE control_flags, ONLY : ngm0
  USE wvfct,         ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  REAL(KIND=DP) :: &
    rho_dot_product                    ! (out) the function value
  COMPLEX(KIND=DP) :: &
    rho1(ngm0,nspin), rho2(ngm0,nspin) ! (in) the two densities
  !
  ! ... and the local variables
  !
  REAL(KIND=DP) :: &
    fac   ! a multiplicative factors
  INTEGER  :: &
    is, ig
  !
  !
  rho_dot_product = 0.D0
  !
  IF ( nspin == 1 ) THEN
     !
     is = 1
     !
     DO ig = gstart, ngm0
        !
        fac = e2 * fpi / ( tpiba2 * gg(ig) )
        !
        rho_dot_product = rho_dot_product + &
                          fac * REAL( CONJG( rho1(ig,is) ) * rho2(ig,is) )
        !
     END DO
     !
     IF ( gamma_only ) rho_dot_product = 2.D0 * rho_dot_product
     !
  ELSE
     !
     DO ig = gstart, ngm0
        !
        fac = e2 * fpi / ( tpiba2 * gg(ig) )
        !
        rho_dot_product = rho_dot_product + &
                          fac * REAL( CONJG( rho1(ig,1) + rho1(ig,2) ) * &
                                           ( rho2(ig,1) + rho2(ig,2) ) )
        !
     END DO
     !
     IF ( gamma_only ) rho_dot_product = 2.D0 * rho_dot_product
     !
     fac = e2 * fpi / tpi**2  ! lambda = 1 a.u.
     !
     ! ... G=0 term
     !
     IF ( gstart == 2 ) THEN
        !
        rho_dot_product = rho_dot_product + &
                          fac * REAL( CONJG( rho1(1,1) - rho1(1,2) ) * &
                                           ( rho2(1,1) - rho2(1,2) ) )
        !
     END IF
     !
     IF ( gamma_only ) fac = 2.D0 * fac
     !
     DO ig = gstart, ngm0
        !
        rho_dot_product = rho_dot_product + &
                          fac * REAL( CONJG( rho1(ig,1) - rho1(ig,2) ) * &
                                           ( rho2(ig,1) - rho2(ig,2) ) )
        !
     END DO
     !
  END IF
  !
  rho_dot_product = rho_dot_product * omega / 2.D0
  !
#if defined (__PARA)
  CALL reduce( 1, rho_dot_product )
#endif
  !
  RETURN
  !
END FUNCTION rho_dot_product
!
!
!----------------------------------------------------------------------------
FUNCTION ns_dot_product( ns1, ns2 )
  !----------------------------------------------------------------------------
  !
  ! ... this function evaluates the dot product between two input densities
  !
  USE kinds,      ONLY : DP
  USE ldaU,       ONLY : lda_plus_u, Hubbard_lmax, Hubbard_l, Hubbard_U, &
                         Hubbard_alpha
  USE ions_base,  ONLY : nat, ityp
  USE lsda_mod,   ONLY : nspin
  !
  IMPLICIT NONE  
  !
  ! ... I/O variables
  !
  REAL(KIND=DP) :: &
    ns_dot_product                                    ! (out) the function value
  REAL(KIND=DP) :: &
    ns1(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat), &
    ns2(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)  ! (in) the two ns 
  !
  ! ... and the local variables
  !
  REAL(KIND=DP) :: &
    sum
  INTEGER :: &
    na, nt, is, m1, m2
  !
  !
  ns_dot_product = 0.D0
  !
  IF ( .NOT. lda_plus_u ) RETURN
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     IF ( Hubbard_U(nt) /= 0.D0 .OR. Hubbard_alpha(nt) /= 0.D0 ) THEN
        !
        sum = 0.D0
        !
        DO is = 1, nspin
           DO m1 = 1, ( 2 * Hubbard_l(nt) + 1 )
              DO m2 = m1, ( 2 * Hubbard_l(nt) + 1 )
                 !
                 sum = sum + ns1(m1,m2,is,na) * ns2(m2,m1,is,na)
                 !
              END DO
           END DO
        END DO
        !
        ns_dot_product = ns_dot_product + 0.5D0 * Hubbard_U(nt) * sum
        !
     END IF
     !
  END DO
  !
  IF (nspin == 1) ns_dot_product = 2.D0 * ns_dot_product
  !
  RETURN
  !
END FUNCTION ns_dot_product

!----------------------------------------------------------------------------
FUNCTION fn_dehar( drho )
  !----------------------------------------------------------------------------
  !
  ! ... this function evaluates the residual hartree energy of drho
  !
  USE kinds,          ONLY : DP
  USE constants,      ONLY : e2, fpi
  USE cell_base,      ONLY : omega, tpiba2
  USE gvect,          ONLY : gstart, gg
  USE lsda_mod,       ONLY : nspin
  USE control_flags,  ONLY : ngm0
  USE wvfct,          ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  REAL(KIND=DP) :: &
    fn_dehar         ! (out) the function value
  COMPLEX(KIND=DP) :: &
    drho(ngm0,nspin) ! (in) the density difference
  !
  ! ... and the local variables
  !
  REAL(KIND=DP) :: &
    fac   ! a multiplicative factors
  INTEGER :: &
    is, ig
  !
  !
  fn_dehar = 0.D0
  !
  IF ( nspin == 1 ) THEN
     !
     is=1
     !
     DO ig = gstart, ngm0
        !
        fac = e2 * fpi / ( tpiba2 * gg(ig) )
        !
        fn_dehar = fn_dehar +  fac * ABS( drho(ig,is) )**2
        !
     END DO
     !
  ELSE
     !
     DO ig = gstart, ngm0
        !
        fac = e2 * fpi / ( tpiba2 * gg(ig) )
        !
        fn_dehar = fn_dehar +  fac * ABS( drho(ig,1) + drho(ig,2) )**2
        !
     END DO
     !
  END IF
  !
  IF ( gamma_only ) THEN
     !
     fn_dehar = fn_dehar * omega
     !
  ELSE
     !
     fn_dehar = fn_dehar * omega / 2.D0
     !
  END IF
  !
#if defined (__PARA)
  CALL reduce( 1, fn_dehar )
#endif
  !
  RETURN
  !
END FUNCTION fn_dehar
!
!
!----------------------------------------------------------------------------
SUBROUTINE approx_screening( drho )
  !----------------------------------------------------------------------------
  !
  ! ... apply an average TF preconditioning to drho
  !
  USE kinds,          ONLY : DP
  USE constants,      ONLY : e2, pi, fpi
  USE cell_base,      ONLY : omega, tpiba2
  USE gvect,          ONLY : gstart, gg
  USE klist,          ONLY : nelec
  USE lsda_mod,       ONLY : nspin
  USE control_flags,  ONLY : ngm0
  !
  IMPLICIT NONE  
  !
  ! ... I/O variables
  !
  COMPLEX(KIND=DP) :: &
    drho(ngm0,nspin) ! (in/out)
  !
  ! ... and the local variables
  !
  REAL(KIND=DP) :: &
    rrho, rmag, rs, agg0
  INTEGER :: &
    is, ig
  !
  !
  rs = ( 3.D0 * omega / fpi / nelec )**( 1.D0 / 3.D0 )
  !
  agg0 = ( 12.D0 / pi )**( 2.D0 / 3.D0 ) / tpiba2 / rs
  !
#if defined (DEBUG)
  WRITE( stdout,'(A,F12.6,A,F12.6)') ' avg rs  =', rs, ' avg rho =', nelec/omega
#endif
  !
  IF ( nspin == 1 ) THEN
     !
     drho(:,1) =  drho(:,1) * gg(:) / ( gg(:) + agg0 )
     !
  ELSE
     !
     DO ig = 1, ngm0
        !
        rrho = ( drho(ig,1) + drho(ig,2) ) * gg(ig) / ( gg(ig) + agg0 )
        rmag = ( drho(ig,1) - drho(ig,2) )
        drho(ig,1) =  0.5D0 * ( rrho + rmag )
        drho(ig,2) =  0.5D0 * ( rrho - rmag )
        !
     END DO
     !
  END IF
  !
  RETURN
  !
END SUBROUTINE approx_screening
!
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
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   nl, nlm, gg
  USE klist,                ONLY : nelec
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : ngm0
  USE wvfct,                ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  COMPLEX(KIND=DP) :: &
    drho(ngm0,nspin), rhobest(ngm0,nspin)
  !
  ! ... and the local variables
  !
  INTEGER, PARAMETER :: & 
    mmx = 12
  INTEGER :: &
    iwork(mmx), i, j, m, info, nspin_save
  REAL(KIND=DP) :: &
    rs, min_rs, max_rs, avg_rsm1, target, dr2_best, ccc, cbest, l2smooth
  REAL(KIND=DP) :: &
    aa(mmx,mmx), invaa(mmx,mmx), bb(mmx), work(mmx), vec(mmx), agg0
  COMPLEX(KIND=DP) :: &
    rrho, rmag
  COMPLEX(KIND=DP), ALLOCATABLE :: &
    v(:,:),     &! v(ngm0,mmx)
    w(:,:),     &! w(ngm0,mmx)
    dv(:),      &! dv(ngm0)
    vbest(:),   &! vbest(ngm0)
    wbest(:)     ! wbest(ngm0)
  REAL(KIND=DP), ALLOCATABLE :: &
    alpha(:)     ! alpha(nrxx)
  INTEGER :: &
    is, ir, ig
  REAL(KIND=DP), EXTERNAL :: rho_dot_product
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
  nspin_save = nspin
  nspin      = 1
  is         = 1
  target     = 0.D0
  !
  IF ( gg(1) < eps8 ) drho(1,is) = ( 0.D0, 0.D0 )
  !
  ALLOCATE( alpha(nrxx), v(ngm0,mmx), w(ngm0,mmx), &
            dv(ngm0), vbest(ngm0), wbest(ngm0) )
  !
  v(:,:)   = ( 0.D0, 0.D0 )
  w(:,:)   = ( 0.D0, 0.D0 )
  dv(:)    = ( 0.D0, 0.D0 )
  vbest(:) = ( 0.D0, 0.D0 )
  wbest(:) = ( 0.D0, 0.D0 )
  !
  ! ... calculate alpha from density smoothed with a lambda=0 a.u.
  !
  l2smooth = 0.D0
  psic(:)  = ( 0.D0, 0.D0 )
  !
  IF ( nspin == 1 ) THEN
     !
     psic(nl(:)) = rhobest(:,1) * EXP( - 0.5D0 * l2smooth * tpiba2 * gg(:) )
     !
  ELSE
     !
     psic(nl(:)) = ( rhobest(:,1) + rhobest(:,2) ) * &
                   EXP( - 0.5D0 * l2smooth * tpiba2 * gg(:) )
     !
  END IF
  !
  IF ( gamma_only ) THEN
     !
     psic(nlm(:)) = CONJG( psic(nl(:)) )
     !
  END IF
  !
  CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
  !
  alpha(:) = REAL( psic(:) )
  !
  min_rs   = ( 3.D0 * omega / fpi / nelec )**( 1.D0 / 3.D0 )
  max_rs   = min_rs
  avg_rsm1 = 0.D0
  !
  DO ir = 1, nrxx
     !
     alpha(ir) = ABS( alpha(ir) )
     !
     IF ( alpha(ir) > eps32 ) THEN
        !
        rs        = ( 3.D0 / fpi / alpha(ir) )**( 1.D0 /3.D0 )
        min_rs    = MIN( min_rs, rs )
        avg_rsm1  = avg_rsm1 + 1.D0 / rs
        max_rs    = MAX( max_rs, rs )
        alpha(ir) = rs
        !
     END IF   
     !
  END DO
  !
#if defined (__PARA)
  CALL reduce( 1, avg_rsm1 )
  CALL extreme( min_rs, -1 )
  CALL extreme( max_rs, +1 )
#endif
  !
  CALL DSCAL( nrxx, 3.D0*( tpi / 3.D0 )**( 5.D0 / 3.D0 ), alpha, 1 )
  !
  avg_rsm1 = ( nr1 * nr2 * nr3 ) / avg_rsm1
  rs       = ( 3.D0 * omega / fpi / nelec )**( 1.D0 / 3.D0 )
  agg0     = ( 12.D0 / pi )**( 2.D0 / 3.D0 ) / tpiba2 / avg_rsm1
  !
#if defined (DEBUG)
  WRITE( stdout,'(A,5F12.6)') ' min/avgm1/max rs  =', min_rs,avg_rsm1,max_rs,rs
#endif
  !
  ! ... calculate deltaV and the first correction vector
  !
  psic(:) = ( 0.D0, 0.D0 )
  !
  psic(nl(:)) = drho(:,is)
  !
  IF ( gamma_only ) &
     psic(nlm(:)) = CONJG( psic(nl(:)) )
  !
  CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
  !
  psic(:) = psic(:) * alpha(:)
  !
  CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1 )
  !
  dv(:) = psic(nl(:)) * gg(:) * tpiba2
  v(:,1)= psic(nl(:)) * gg(:) / ( gg(:) + agg0 )
  !
  m       = 1
  ccc     = rho_dot_product( dv, dv )
  aa(:,:) = 0.D0
  bb(:)   = 0.D0
  !
  repeat_loop: DO
     !
     ! ... generate the vector w
     !
     w(:,m) = gg(:) * tpiba2 * v(:,m)
     !
     psic(:) = ( 0.D0, 0.D0 )
     !
     psic(nl(:)) = v(:,m)
     !
     IF ( gamma_only ) &
        psic(nlm(:)) = conjg( psic(nl(:)) )
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
     !
     WHERE( alpha(:) > eps32 )
        !
        psic(:) = psic(:) * fpi * e2 / alpha(:)
        !
     END WHERE
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1 )
     !
     w(:,m) = w(:,m) + psic(nl(:))
     !
     ! ... build the linear system
     !
     DO i = 1, m
        aa(i,m) = rho_dot_product( w(1,i), w(1,m) )
        aa(m,i) = aa(i,m)
     END DO
     !
     bb(m) = rho_dot_product( w(1,m), dv )
     !
     ! ... solve it -> vec
     !
     CALL DCOPY( mmx * mmx, aa, 1, invaa, 1 )
     CALL DSYTRF( 'U', m, invaa, mmx, iwork, work, mmx, info )
     CALL errore( 'broyden', 'factorization', info )
     CALL DSYTRI( 'U', m, invaa, mmx, iwork, work, info )
     CALL errore( 'broyden', 'DSYTRI', info )
     !
     DO i = 1, m
        DO j = ( i + 1 ), m
           !
           invaa(j,i) = invaa(i,j)
           !
        END DO
     END DO
     !
     DO i = 1, m
        !
        vec(i) = 0.D0
        !
        DO j = 1, m
           vec(i) = vec(i) + invaa(i,j) * bb(j)
        END DO
        !
     END DO
     !
     vbest(:) = ( 0.D0, 0.D0 )
     wbest(:) = dv(:)
     !
     DO i = 1, m
        CALL DAXPY( 2 * ngm0,  vec(i), v(1,i), 1, vbest, 1 )
        CALL DAXPY( 2 * ngm0, -vec(i), w(1,i), 1, wbest, 1 )
     END DO
     !
     cbest = ccc
     !
     DO i = 1, m
        cbest = cbest - bb(i) * vec(i)
     END DO
     !
     dr2_best = rho_dot_product( wbest, wbest )
     !
     IF ( target == 0.D0 ) target = 1.D-6 * dr2_best
     !
     IF ( dr2_best < target ) THEN
        !
        psic(:) = ( 0.D0, 0.D0 )
        !
        psic(nl(:)) = vbest(:)
        !
        IF ( gamma_only ) &
           psic(nlm(:)) = CONJG( psic(nl(:)) )
        !
        CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1 )
        !
        WHERE( alpha(:) > eps32 )
           !
           psic(:) = psic(:) / alpha(:)
           !
        END WHERE
        !
        CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1 )
        !
        drho(:,is) = psic(nl(:))
        !
        nspin = nspin_save
        !
        IF ( nspin == 2 ) THEN
           DO ig = 1, ngm0
              rrho = drho(ig,1)
              rmag = drho(ig,2)
              drho(ig,1) = 0.5D0 * ( rrho + rmag )
              drho(ig,2) = 0.5D0 * ( rrho - rmag )
           END DO
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
     v(:,m) = wbest(:) / ( gg(:) + agg0 )
     !
  END DO repeat_loop
  !
  RETURN
  !
END SUBROUTINE approx_screening2
