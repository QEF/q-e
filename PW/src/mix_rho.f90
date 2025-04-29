!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#define ZERO ( 0._dp, 0._dp )
!
! This macro force the normalization of betamix matrix, usually not necessary
!#define __NORMALIZE_BETAMIX
!
!----------------------------------------------------------------------------
SUBROUTINE mix_rho( input_rhout, rhoin, alphamix, dr2, tr2_min, iter, n_iter,&
                    iunmix, conv )
  !----------------------------------------------------------------------------
  !! * Modified Broyden's method for charge density mixing: D.D. Johnson,
  !!   PRB 38, 12807 (1988) ;
  !! * Thomas-Fermi preconditioning described in: Raczkowski, Canning, Wang,
  !!   PRB 64,121101 (2001) ;
  !! * Extended to mix also quantities needed for PAW, meta-GGA, DFT+U(+V) ;
  !! * Electric field (all these are included into \(\text{mix_type}\)) ;
  !! * On output: the mixed density is in \(\text{rhoin}\), 
  !!   \(\text{input_rhout}\) is unchanged.
  !
  USE kinds,          ONLY : DP
  USE ions_base,      ONLY : nat, ityp, ntyp => nsp
  USE gvect,          ONLY : ngm
  USE gvecs,          ONLY : ngms
  USE lsda_mod,       ONLY : nspin
  USE control_flags,  ONLY : imix, ngm0, tr2, io_level
  ! ... for PAW:
  USE uspp_param,     ONLY : nhm
  USE ener,           ONLY : ef
  USE gcscf_module,   ONLY : lgcscf, gcscf_gh, gcscf_mu, gcscf_eps
  USE scf,            ONLY : scf_type, create_scf_type, destroy_scf_type, &
                             mix_type, create_mix_type, destroy_mix_type, &
                             assign_scf_to_mix_type, assign_mix_to_scf_type, &
                             mix_type_AXPY, davcio_mix_type, rho_ddot, &
                             high_frequency_mixing, nsg_ddot, &
                             mix_type_COPY, mix_type_SCAL
  USE io_global,      ONLY : stdout
  USE gcscf_module,   ONLY : lgcscf, gcscf_gh, gcscf_mu, gcscf_eps
  USE ldaU,           ONLY : lda_plus_u, lda_plus_u_kind, ldim_u, nsnew, neighood, &
                             max_num_neighbors, nsg, nsgnew, Hubbard_l, Hubbard_lmax
  USE buffers,        ONLY : open_buffer, close_buffer, get_buffer, save_buffer
#if defined (__OSCDFT)
  USE plugin_flags,     ONLY : use_oscdft
  USE oscdft_base,      ONLY : oscdft_ctx
  USE oscdft_functions, ONLY : oscdft_mix_rho
  USE ldaU,             ONLY : ldmx
  USE oscdft_functions, ONLY : oscdft_constrain_ns
#endif
  !
  IMPLICIT NONE
  !
  ! ... First the I/O variable
  !
  INTEGER, INTENT(IN) :: iter
  !! counter of the number of iterations
  INTEGER, INTENT(IN) :: n_iter
  !! number of iterations used in mixing
  INTEGER, INTENT(IN) :: iunmix
  !! I/O unit where data from previous iterations is stored
  REAL(DP), INTENT(IN) :: alphamix
  !! mixing factor
  REAL(DP), INTENT(IN) :: tr2_min
  !! estimated error in diagonalization. If the estimated
  !! scf error is smaller than this, exit: a more accurate 
  !! diagonalization is needed
  REAL(DP), INTENT(OUT) :: dr2
  !! the estimated error on the energy
  LOGICAL, INTENT(OUT) :: conv
  !! .TRUE. if the convergence has been reached
  !
  TYPE(scf_type), INTENT(INOUT) :: input_rhout
  TYPE(scf_type), INTENT(INOUT) :: rhoin
  !
  ! ... local variables
  !
  TYPE(mix_type) :: rhout_m, rhoin_m
  INTEGER, PARAMETER :: &
    maxmix = 25     ! max number of iterations for charge mixing
  INTEGER ::       &
    iter_used,     &! actual number of iterations used
    ipos,          &! index of the present iteration
    inext,         &! index of the next iteration
    i, j,          &! counters on number of iterations
    info,          &! flag saying if the exec. of libr. routines was ok
    ldim,          &! 2 * Hubbard_lmax + 1
    iunmix_nsg,    &! the unit for Hubbard mixing within DFT+U+V
    nt,            &! index of the atomic type
    nword           ! size the DFT+U+V-related arrays
  REAL(DP),ALLOCATABLE :: betamix(:,:), work(:)
  INTEGER, ALLOCATABLE :: iwork(:)
  COMPLEX(DP), ALLOCATABLE :: nsginsave(:,:,:,:,:),  nsgoutsave(:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: deltansg(:,:,:,:,:)
  LOGICAL :: exst, exst_mem, exst_file
  REAL(DP) :: gamma0
#if defined(__NORMALIZE_BETAMIX)
  REAL(DP) :: norm2, obn
#endif
  !
  ! ... saved variables and arrays
  !
  INTEGER, SAVE :: &
    mixrho_iter = 0    ! history of mixing
  TYPE(mix_type), ALLOCATABLE, SAVE :: &
    df(:),        &! information from preceding iterations
    dv(:)          !     "  "       "     "        "  "
  REAL(DP) :: norm
  INTEGER, PARAMETER :: read_ = -1, write_ = +1
  !
  ! ... external functions
  !
  INTEGER, EXTERNAL :: find_free_unit
  COMPLEX(DP), ALLOCATABLE :: df_nsg(:,:,:,:,:,:), dv_nsg(:,:,:,:,:,:)
  !
#if defined (__OSCDFT)
  INTEGER :: iunmix_cons, & ! the unit for the constraint mixing
             nword_cons
  LOGICAL :: exst_mem_cons, exst_file_cons
  COMPLEX(DP), ALLOCATABLE :: delta_cons(:,:,:,:), t_cons1(:,:,:, :),   &
                              t_cons2(:,:,:, :), df_cons(:,:,:,:,:),    &
                              dv_cons(:,:,:,:,:), cons_insave(:,:,:,:), & 
                              cons_outsave(:,:,:,:)
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
  ! define mix_type variables and copy scf_type variables there
  !
  call create_mix_type(rhout_m)
  call create_mix_type(rhoin_m)
  !
  ! DFT+U+V case
  !
  IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) THEN
     ldim = 0
     DO nt = 1, ntyp
        ldim = MAX(ldim, ldim_u(nt))
     ENDDO
     ALLOCATE ( deltansg (ldim, ldim, max_num_neighbors, nat, nspin) )
     deltansg(:,:,:,:,:) = nsgnew(:,:,:,:,:) - nsg(:,:,:,:,:)
  ENDIF
  !
  call assign_scf_to_mix_type(rhoin, rhoin_m)
  call assign_scf_to_mix_type(input_rhout, rhout_m)
  call mix_type_AXPY ( -1.d0, rhoin_m, rhout_m )
  !
  IF ( lgcscf ) THEN
     !
     dr2 = rho_ddot( rhout_m, rhout_m, ngms, gcscf_gh )
     !
  ELSE
     !
     dr2 = rho_ddot( rhout_m, rhout_m, ngms )  !!!! this used to be ngm NOT ngms
     !
  END IF
  !
  IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) &
      dr2 = dr2 + nsg_ddot( deltansg, deltansg, nspin )
  !
  IF (dr2 < 0.0_DP) CALL errore('mix_rho','negative dr2',1)
  !
  conv = ( dr2 < tr2 )
  !
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     IF (conv .AND. .NOT.oscdft_ctx%conv) conv = .FALSE.
  ENDIF
#endif
  !
  IF ( lgcscf ) THEN
     !
     conv = conv .AND. ( ABS( ef - gcscf_mu ) < gcscf_eps )
     !
  END IF
#if defined(__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) CALL oscdft_mix_rho(oscdft_ctx, conv)
#endif
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
     !
     call destroy_mix_type(rhoin_m)
     call destroy_mix_type(rhout_m)
     !
     IF ( ALLOCATED( dv_nsg ) ) DEALLOCATE( dv_nsg )
     IF ( ALLOCATED( df_nsg ) ) DEALLOCATE( df_nsg )
     IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) THEN
        nsgnew(:,:,:,:,:) = deltansg(:,:,:,:,:) + nsg(:,:,:,:,:)
        DEALLOCATE(deltansg)
     ENDIF
     !
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     IF (oscdft_ctx%is_constraint .AND. oscdft_ctx%conv) THEN
        IF (ALLOCATED(dv_cons) ) DEALLOCATE( dv_cons )
        IF (ALLOCATED(df_cons) ) DEALLOCATE( df_cons )
        IF (ALLOCATED(delta_cons) ) DEALLOCATE( delta_cons )
     ENDIF
  ENDIF
#endif   
     ! 
     CALL stop_clock( 'mix_rho' )
     !
     RETURN
     !
  END IF
  !
  IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) THEN
     iunmix_nsg = find_free_unit()
     nword = ldim * ldim * max_num_neighbors * nat * nspin * n_iter
     CALL open_buffer( iunmix_nsg, 'mix.nsg', nword, io_level, exst_mem, exst_file)
     ALLOCATE( df_nsg(ldim,ldim,max_num_neighbors,nat,nspin,n_iter) )
     ALLOCATE( dv_nsg(ldim,ldim,max_num_neighbors,nat,nspin,n_iter) )
  ENDIF
  !
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     ALLOCATE (t_cons1(ldmx, ldmx, nspin, nat))
     ALLOCATE (t_cons2(ldmx, ldmx, nspin, nat))
     ALLOCATE (delta_cons(ldmx, ldmx, nspin, nat))
     t_cons1(:,:,:,:) = oscdft_ctx%constraint(:,:,:,:)
     IF (lda_plus_u) THEN
        IF (lda_plus_u_kind==0) THEN
           CALL oscdft_constrain_ns (oscdft_ctx, Hubbard_lmax, Hubbard_l, &
                          input_rhout%ns, oscdft_ctx%conv, dr2, tr2, iter)
        ELSEIF (lda_plus_u_kind==2) THEN
           ALLOCATE (nsnew(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat))
           CALL oscdft_nsg(3)
           CALL oscdft_constrain_ns (oscdft_ctx, Hubbard_lmax, Hubbard_l, &
                       nsnew, oscdft_ctx%conv, dr2, tr2, iter)
           DEALLOCATE (nsnew)
        ENDIF
     ENDIF
     IF (.NOT. oscdft_ctx%conv) THEN
        delta_cons(:,:,:,:) = oscdft_ctx%constraint(:,:,:,:) - t_cons1(:,:,:,:)
        t_cons2(:,:,:,:) = oscdft_ctx%constraint(:,:,:,:)
     ENDIF
     IF (.NOT. oscdft_ctx%conv) THEN
        iunmix_cons = find_free_unit() + 1
        nword_cons = ldmx * ldmx * nat * nspin * n_iter
        CALL open_buffer( iunmix_cons, 'mix.hubcons', nword_cons, io_level, exst_mem_cons, exst_file_cons)
        ALLOCATE( df_cons(ldmx,ldmx,nspin, nat, n_iter) )
        ALLOCATE( dv_cons(ldmx,ldmx,nspin, nat, n_iter) )
     ENDIF
  ENDIF
#endif
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
     CALL davcio_mix_type( df(ipos), iunmix, 1, read_ )
     CALL davcio_mix_type( dv(ipos), iunmix, 2, read_ )
     !
     call mix_type_AXPY ( -1.d0, rhout_m, df(ipos) )
     call mix_type_AXPY ( -1.d0, rhoin_m, dv(ipos) )
     !
     IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) THEN
        CALL get_buffer ( df_nsg, nword, iunmix_nsg, 1 )
        CALL get_buffer ( dv_nsg, nword, iunmix_nsg, 2 )
        df_nsg(:,:,:,:,:,ipos) = df_nsg(:,:,:,:,:,ipos) - deltansg
        dv_nsg(:,:,:,:,:,ipos) = dv_nsg(:,:,:,:,:,ipos) - nsg
     ENDIF
     !
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     IF (oscdft_ctx%is_constraint .AND. .NOT.oscdft_ctx%conv) THEN
        CALL get_buffer ( df_cons, nword_cons, iunmix_cons, 1 )
        CALL get_buffer ( dv_cons, nword_cons, iunmix_cons, 2 )
        df_cons(:,:,:,:,ipos) = df_cons(:,:,:,:,ipos) - delta_cons
        dv_cons(:,:,:,:,ipos) = dv_cons(:,:,:,:,ipos) - t_cons2
     ENDIF
  ENDIF
#endif
     !
#if defined(__NORMALIZE_BETAMIX)
     ! NORMALIZE
     ! TODO: need to check compatibility gcscf and lda_plus_u
     IF ( lgcscf ) THEN
        norm2 = rho_ddot( df(ipos), df(ipos), ngm0, gcscf_gh )
     ELSE
        norm2 = rho_ddot( df(ipos), df(ipos), ngm0 )
     END IF
     IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) THEN
        norm2 = norm2 + nsg_ddot( df_nsg(1,1,1,1,1,ipos), df_nsg(1,1,1,1,1,ipos), nspin )
     ENDIF
     obn = 1.d0/sqrt(norm2)
     call mix_type_SCAL (obn,df(ipos))
     call mix_type_SCAL (obn,dv(ipos))
     IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) THEN
        df_nsg(:,:,:,:,:,ipos) = df_nsg(:,:,:,:,:,ipos) * obn
        dv_nsg(:,:,:,:,:,ipos) = dv_nsg(:,:,:,:,:,ipos) * obn
     ENDIF
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     IF (oscdft_ctx%is_constraint .AND. .NOT.oscdft_ctx%conv) THEN
        df_cons(:,:,:,:,ipos) = df_cons(:,:,:,:,ipos) * obn
        dv_cons(:,:,:,:,ipos) = dv_cons(:,:,:,:,ipos) * obn
     ENDIF
  ENDIF
#endif
#endif
     !
  END IF
  !
  DO i = 1, iter_used
     !
     IF ( i /= ipos ) THEN
        !
        CALL davcio_mix_type( df(i), iunmix, 2*i+1, read_ )
        CALL davcio_mix_type( dv(i), iunmix, 2*i+2, read_ )
     END IF
     !
  END DO
  !
  CALL davcio_mix_type( rhout_m, iunmix, 1, write_ )
  CALL davcio_mix_type( rhoin_m, iunmix, 2, write_ )
  !
  IF ( mixrho_iter > 1 ) THEN
     CALL davcio_mix_type( df(ipos), iunmix, 2*ipos+1, write_ )
     CALL davcio_mix_type( dv(ipos), iunmix, 2*ipos+2, write_ )
  END IF
  !
  IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) THEN 
     !
     ALLOCATE( nsginsave(  ldim,ldim,max_num_neighbors,nat,nspin ), &
               nsgoutsave( ldim,ldim,max_num_neighbors,nat,nspin ) )
     nsginsave  = (0.d0,0.d0)
     nsgoutsave = (0.d0,0.d0)
     nsginsave  = nsg
     nsgoutsave = deltansg !nsgnew
     !
  ENDIF
  !
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     IF (oscdft_ctx%is_constraint .AND. .NOT.oscdft_ctx%conv) THEN
        ALLOCATE( cons_insave(  ldmx,ldmx,nspin, nat) )
        ALLOCATE( cons_outsave( ldmx,ldmx,nspin, nat ) )
        cons_insave  = (0.d0, 0.d0)
        cons_outsave = (0.d0, 0.d0)
        cons_insave(:, :, :, :)  = t_cons2(:, :, :, :)
        cons_outsave(:, :, :, :) = delta_cons(:, :, :, :)
     ENDIF
  ENDIF
#endif
  !
  ! Nothing else to do on first iteration
  skip_on_first: &
  IF (iter_used > 0) THEN
    !
    ALLOCATE(betamix(iter_used, iter_used)) !iter_used))
    betamix = 0._dp
    !
    DO i = 1, iter_used
        !
        DO j = i, iter_used
            !
            IF ( lgcscf ) THEN
               !
               betamix(i,j) = rho_ddot( df(j), df(i), ngm0, gcscf_gh )
               !
            ELSE
               !
               betamix(i,j) = rho_ddot( df(j), df(i), ngm0 )
               !
            END IF
            !
            IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) &
               betamix(i,j) = betamix(i,j) + &
                  nsg_ddot( df_nsg(1,1,1,1,1,j), df_nsg(1,1,1,1,1,i), nspin )
            !
            betamix(j,i) = betamix(i,j)
            !
        END DO
        !
    END DO
    !
    allocate(work(iter_used), iwork(iter_used))
    CALL DSYTRF( 'U', iter_used, betamix, iter_used, iwork, work, iter_used, info )
    CALL errore( 'broyden', 'factorization', abs(info) )
    !
    CALL DSYTRI( 'U', iter_used, betamix, iter_used, iwork, work, info )
    CALL errore( 'broyden', 'DSYTRI', abs(info) )    !
    deallocate(iwork)
    !
    FORALL( i = 1:iter_used, &
            j = 1:iter_used, j > i ) betamix(j,i) = betamix(i,j)
    !
    DO i = 1, iter_used
       !
       IF ( lgcscf ) THEN
          !
          work(i) = rho_ddot( df(i), rhout_m, ngm0, gcscf_gh )
          !
       ELSE
          !
          work(i) = rho_ddot( df(i), rhout_m, ngm0 )
          !
       END IF
       !
        IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) &
           work(i) = work(i) + nsg_ddot( df_nsg(1,1,1,1,1,i), deltansg, nspin )
    END DO
    !
    DO i = 1, iter_used
        !
        gamma0 = DOT_PRODUCT( betamix(1:iter_used,i), work(1:iter_used) )
        !
        call mix_type_AXPY ( -gamma0, dv(i), rhoin_m )
        call mix_type_AXPY ( -gamma0, df(i), rhout_m )
        !
        IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) THEN
           nsg(:,:,:,:,:) = nsg(:,:,:,:,:) - gamma0*dv_nsg(:,:,:,:,:,i)
           deltansg(:,:,:,:,:) = deltansg(:,:,:,:,:) - gamma0*df_nsg(:,:,:,:,:,i)
        ENDIF
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     IF (oscdft_ctx%is_constraint .AND. .NOT.oscdft_ctx%conv) THEN
        t_cons2(:, :, :, :) = t_cons2(:, :, :, :) - gamma0 * dv_cons(:, :, :, :, i)
        delta_cons(:, :, :, :) = delta_cons(:, :, :, :) - gamma0 * df_cons(:, :, :, :, i)
     ENDIF
  ENDIF
#endif
        !
    END DO
    DEALLOCATE(betamix, work)
    !
    ! ... auxiliary vectors dv and df not needed anymore
    !
  ENDIF skip_on_first
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
  IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) THEN
     inext = mixrho_iter - ( ( mixrho_iter - 1 ) / n_iter ) * n_iter
     df_nsg(:,:,:,:,:,inext) = nsgoutsave
     dv_nsg(:,:,:,:,:,inext) = nsginsave
     IF (ALLOCATED(nsginsave))  DEALLOCATE(nsginsave)
     IF (ALLOCATED(nsgoutsave)) DEALLOCATE(nsgoutsave)
  ENDIF
  !
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     IF (oscdft_ctx%is_constraint .AND. .NOT.oscdft_ctx%conv) THEN
        inext = mixrho_iter - ( ( mixrho_iter - 1 ) / n_iter ) * n_iter
        df_cons(:,:,:,:,inext) = cons_outsave(:, :, :, :)
        dv_cons(:,:,:,:,inext) = cons_insave(:, :, :, :)
     ENDIF
     IF (ALLOCATED(cons_insave))  DEALLOCATE(cons_insave)
     IF (ALLOCATED(cons_outsave)) DEALLOCATE(cons_outsave)
  ENDIF
#endif
  !
  ! ... preconditioning the new search direction
  !
  IF ( imix == 1 ) THEN
     !
     CALL approx_screening( rhout_m )
     !
  ELSE IF ( imix == 2 ) THEN
     !
     CALL approx_screening2( rhout_m, rhoin_m )
     !
  END IF
  !
  ! ... set new trial density
  !
  call mix_type_AXPY ( alphamix, rhout_m, rhoin_m )
  IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) THEN
     nsg = nsg + alphamix * deltansg
     IF (ALLOCATED(deltansg)) DEALLOCATE(deltansg)
  ENDIF
  !
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     IF (oscdft_ctx%is_constraint) THEN
        IF (.NOT.oscdft_ctx%conv) THEN
           oscdft_ctx%constraint(:, :, :, :) = t_cons2(:, :, :, :) + &
                                     alphamix * delta_cons(:, :, :, :)
        ELSE
           oscdft_ctx%constraint(:, :, :, :) = 0.8 * t_cons1(:, :, :, :) 
        ENDIF
     ENDIF
  ENDIF
#endif
  !
  ! ... simple mixing for high_frequencies (and set to zero the smooth ones)
  call high_frequency_mixing ( rhoin, input_rhout, alphamix )
  ! ... add the mixed rho for the smooth frequencies
  call assign_mix_to_scf_type(rhoin_m,rhoin)
  !
  call destroy_mix_type(rhout_m)
  call destroy_mix_type(rhoin_m)
  !
  IF (lda_plus_u .AND. lda_plus_u_kind.EQ.2) THEN
     CALL save_buffer ( df_nsg, nword, iunmix_nsg, 1 )
     CALL save_buffer ( dv_nsg, nword, iunmix_nsg, 2 )
     DEALLOCATE( dv_nsg )
     DEALLOCATE( df_nsg )
     CALL close_buffer(iunmix_nsg, 'keep')
  ENDIF
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     IF (oscdft_ctx%is_constraint .AND. .NOT.oscdft_ctx%conv) THEN
        CALL save_buffer ( df_cons, nword_cons, iunmix_cons, 1 )
        CALL save_buffer ( dv_cons, nword_cons, iunmix_cons, 2 )
        DEALLOCATE( dv_cons )
        DEALLOCATE( df_cons )
        CALL close_buffer(iunmix_cons, 'keep')
     ENDIF
     IF (ALLOCATED(t_cons1)) DEALLOCATE( t_cons1 )
     IF (ALLOCATED(t_cons2)) DEALLOCATE( t_cons2 )
  ENDIF
#endif
  !
  CALL stop_clock( 'mix_rho' )
  !
  RETURN
  !
END SUBROUTINE mix_rho
!
!----------------------------------------------------------------------------
SUBROUTINE approx_screening( drho )
  !----------------------------------------------------------------------------
  !! Apply an average TF preconditioning to \(\text{drho}\).
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, pi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gg, ngm
  USE klist,         ONLY : nelec
  USE control_flags, ONLY : ngm0
  USE scf,           ONLY : mix_type
  USE gcscf_module,  ONLY : lgcscf, gcscf_gk
  !
  IMPLICIT NONE
  !
  type (mix_type), intent(INOUT) :: drho ! (in/out)
  !$acc declare present(drho, drho%of_g, gg) 
  !
  REAL(DP) :: rs, agg0, bgg0
  !
  rs = ( 3.D0 * omega / fpi / nelec )**( 1.D0 / 3.D0 )
  !
  agg0 = ( 12.D0 / pi )**( 2.D0 / 3.D0 ) / tpiba2 / rs
  !
  IF ( lgcscf ) THEN
     !
     bgg0 = gcscf_gk * gcscf_gk / tpiba2
     !$acc kernels
     drho%of_g(:ngm0,1) =  drho%of_g(:ngm0,1) * (gg(:ngm0)+bgg0) &
                        / (gg(:ngm0)+agg0+bgg0)
     !$acc end kernels
     !
  ELSE
     !
     !$acc kernels
     drho%of_g(:ngm0,1) =  drho%of_g(:ngm0,1) * gg(:ngm0) / (gg(:ngm0)+agg0)
     !$acc end kernels
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
  !! Apply a local-density dependent TF preconditioning to \(\text{drho}\).
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : e2, pi, tpi, fpi, eps8, eps32
  USE cell_base,            ONLY : omega, tpiba2
  USE gvect,                ONLY : gg, ngm
  USE klist,                ONLY : nelec
  USE control_flags,        ONLY : ngm0, gamma_only
  USE scf,                  ONLY : mix_type, local_tf_ddot
  USE mp,                   ONLY : mp_sum
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE fft_base,             ONLY : dffts
  USE fft_rho,              ONLY : rho_r2g, rho_g2r
  USE gcscf_module,         ONLY : lgcscf, gcscf_gk, gcscf_gh
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
    avg_rsm1, target, dr2_best
  REAL(DP) :: &
    aa(mmx,mmx), invaa(mmx,mmx), bb(mmx), work(mmx), vec(mmx), agg0, bgg0
  COMPLEX(DP), ALLOCATABLE :: &
    v(:,:),     &! v(ngm0,mmx)
    w(:,:),     &! w(ngm0,mmx)
    dv(:),      &! dv(ngm0)
    vbest(:),   &! vbest(ngm0)
    wbest(:),   &! wbest(ngm0)
    auxg(:,:)    ! auxg(dffts%nnr,1) 
  REAL(DP), ALLOCATABLE :: &
    alpha(:),   &! alpha(dffts%nnr)
    auxr(:)      ! auxr(dffts%nnr)
  !
  INTEGER             :: ir, ig
  REAL(DP), PARAMETER :: one_third = 1.D0 / 3.D0
  INTEGER :: dffts_nnr  
  INTEGER :: mmx_refreshed
  INTEGER :: MAX_MMX_REFRESHES = 4
  ! parameter setting how many times the solver's iterative space
  ! is refreshed before quitting
  !$acc data present(drho, drho%of_g, rhobest, rhobest%of_g) 
  !

  dffts_nnr = dffts%nnr
  target = 0.D0
  !
  IF ( (.NOT. lgcscf) .AND. gg(1) < eps8 ) THEN 
    !$acc kernels
    drho%of_g(1,1) = ZERO
    !$acc end kernels
  END IF 
  !
  ALLOCATE( auxr(dffts_nnr), auxg(dffts_nnr,1) )
  ALLOCATE( alpha( dffts_nnr ) )
  ALLOCATE( v( ngm0, mmx ), &
            w( ngm0, mmx ), dv( ngm0 ), vbest( ngm0 ), wbest( ngm0 ) )
  !
  !$acc enter data create(v, w, dv, vbest, wbest, auxg, auxr, alpha) 
  CALL rho_g2r( dffts, rhobest%of_g(:,1), auxr )
  !
  avg_rsm1 = 0.D0
  !
#if defined(_OPENACC)
  !$acc parallel loop reduction(+:avg_rsm1)
#else
  !$omp parallel do reduction(+:avg_rsm1)
#endif
  DO ir = 1, dffts_nnr
     alpha(ir) = ABS( auxr(ir) )
     !
     IF ( alpha(ir) > eps32 ) THEN
        !
        alpha(ir) = ( 3.D0 / fpi / alpha(ir) )**one_third
        avg_rsm1  = avg_rsm1 + 1.D0 / alpha(ir)
        !
     END IF
     !
     alpha(ir) = 3.D0 * ( tpi / 3.D0 )**( 5.D0 / 3.D0 ) * alpha(ir)
     !
  END DO
#if !defined(_OPENACC)
  !$omp end parallel do
#endif
  !
  CALL mp_sum( avg_rsm1 , intra_bgrp_comm )
  avg_rsm1 = ( dffts%nr1*dffts%nr2*dffts%nr3 ) / avg_rsm1
  agg0     = ( 12.D0 / pi )**( 2.D0 / 3.D0 ) / tpiba2 / avg_rsm1
  IF ( lgcscf ) bgg0 = gcscf_gk * gcscf_gk / tpiba2
  !
  IF ( lgcscf ) THEN
     !
     bgg0 = gcscf_gk * gcscf_gk / tpiba2
     !
  END IF
  !
  ! ... calculate deltaV and the first correction vector
  !
  CALL rho_g2r( dffts, drho%of_g(:,1), auxr )
  !
#if defined(_OPENACC)
  !$acc parallel loop 
#else
  !$omp parallel do
#endif 
  DO ir = 1, dffts_nnr
     auxr(ir) = auxr(ir) * alpha(ir)
  ENDDO
#if !defined(_OPENACC)
  !$omp end parallel do
#endif
  !
  CALL rho_r2g( dffts, auxr, auxg )
  !
  IF ( lgcscf ) THEN
     !
#if defined (_OPENACC) 
     !$acc parallel loop 
#else
     !$omp parallel do
#endif 
     DO ig = 1, ngm0
        dv(ig) = auxg(ig,1) * ( gg(ig) + bgg0 ) * tpiba2
        v(ig,1)= auxg(ig,1) * ( gg(ig) + bgg0 ) / ( gg(ig) + agg0 + bgg0 )
     ENDDO
#if !defined(_OPENACC) 
     !$omp end parallel do
#endif
     !
  ELSE
     !
#if defined(_OPENACC) 
     !$acc parallel loop
#else
     !$omp parallel do
#endif 
     DO ig = 1, ngm0
        dv(ig) = auxg(ig,1) * gg(ig) * tpiba2
        v(ig,1)= auxg(ig,1) * gg(ig) / ( gg(ig) + agg0 )
     ENDDO
#if !defined(_OPENACC)
     !$omp end parallel do
#endif
     !
  END IF
  !
  m       = 1
  mmx_refreshed = 0 
  aa(:,:) = 0.D0
  bb(:)   = 0.D0
  !
  repeat_loop: DO
     !
     ! ... generate the vector w
     !
#if defined (_OPENACC) 
     !$acc parallel loop 
#else     
     !$omp parallel
        !$omp do
#endif
        DO ig = 1, ngm0
           !
           w(ig,m) = fpi * e2 * v(ig,m)
           !
        ENDDO
#if !defined(_OPENACC) 
        !$omp end do nowait
     !$omp end parallel
#endif
     !
     CALL rho_g2r( dffts, v(:,m), auxr )
     !
#if defined(_OPENACC) 
     !$acc parallel loop
#else
     !$omp parallel do
#endif
     DO ir = 1, dffts_nnr
        auxr(ir) = auxr(ir) * alpha(ir)
     ENDDO
#if !defined(_OPENACC)
     !$omp end parallel do
#endif
     !
     CALL rho_r2g( dffts, auxr, auxg )
     !
     IF ( lgcscf ) THEN
        !
#if defined (_OPENACC)
        !$acc parallel loop 
#else
        !$omp parallel do
#endif
        DO ig = 1, ngm0
           w(ig,m) = w(ig,m) + ( gg(ig) + bgg0 ) * tpiba2 * auxg(ig,1)
        ENDDO
#if !defined(_OPENACC)
        !$omp end parallel do
#endif
        !
     ELSE
        !
#if defined (_OPENACC)
        !$acc parallel loop
#else 
        !$omp parallel do
#endif
        DO ig = 1, ngm0
           w(ig,m) = w(ig,m) + gg(ig) * tpiba2 * auxg(ig,1)
        ENDDO
#if !defined(_OPENACC)
        !$omp end parallel do
#endif
        !
     END IF
     !
     ! ... build the linear system
     !

     DO i = 1, m
        !
        IF ( lgcscf ) THEN
           !
           aa(i,m) = local_tf_ddot( w(1,i), w(1,m), ngm0, gcscf_gh)
           !
        ELSE
           !
           aa(i,m) = local_tf_ddot( w(1,i), w(1,m), ngm0)
           !
        END IF
        !
        aa(m,i) = aa(i,m)
        !
     END DO
     !
     IF ( lgcscf ) THEN
        !
        bb(m) = local_tf_ddot( w(1,m), dv, ngm0, gcscf_gh)
        !
     ELSE
        !
        bb(m) = local_tf_ddot( w(1,m), dv, ngm0)
        !
     END IF
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
#if defined(_OPENACC) 
     !$acc parallel loop
     do ig = 1, ngm0
       vbest(ig) = ZERO 
       wbest(ig) = dv(ig)  
       do i =1, m 
         vbest(ig) = vbest(ig) + vec(i) * v(ig,i) 
         wbest(ig) = wbest(ig) - vec(i) * w(ig,i)
       end do 
     end do 
#else        
     !$omp parallel
        !$omp do
        DO ig = 1, ngm0
           vbest(ig) = ZERO
           wbest(ig) = dv(ig)
        ENDDO
        !$omp end do nowait
        !
        DO i = 1, m
           !$omp do
           DO ig = 1, ngm0
              vbest(ig) = vbest(ig) + vec(i) * v(ig,i)
              wbest(ig) = wbest(ig) - vec(i) * w(ig,i)
           ENDDO
           !$omp end do nowait
        END DO
     !$omp end parallel
#endif 
     !
     IF ( lgcscf ) THEN
        !
        dr2_best = local_tf_ddot( wbest, wbest, ngm0, gcscf_gh )
        !
     ELSE
        ! 
        dr2_best = local_tf_ddot( wbest, wbest, ngm0 )
        !
     END IF
     !
     IF ( target == 0.D0 ) target = MAX( 1.D-12, 1.D-6*dr2_best )
     !
     IF ( dr2_best < target .OR. (& 
          m >=mmx .AND. mmx_refreshed >= MAX_MMX_REFRESHES) & 
        ) THEN
        ! exit if converged or after the solver has been restarted 
        !MAX_MMX_REFRESHES times, avoiding a possible infinite loop
#if defined (_OPENACC) 
        !$acc parallel loop 
#else 
        !$omp parallel
           !$omp do
#endif
           DO ig = 1, ngm0
              drho%of_g(ig,1) = vbest(ig)
           ENDDO
#if !defined(_OPENACC) 
           !$omp end do nowait
           !
        !$omp end parallel
#endif
        !
        !$acc exit data finalize delete(auxr, auxg, alpha, v, w, dv, vbest, wbest) 
        DEALLOCATE( auxr, auxg )
        DEALLOCATE( alpha, v, w, dv, vbest, wbest )
        !
        EXIT repeat_loop
        !
     ELSE IF ( m >= mmx ) THEN
        !
        m = 1
        mmx_refreshed = mmx_refreshed + 1 
        !
#if defined(_OPENACC) 
        !$acc parallel loop
#else
        !$omp parallel do
#endif 
        DO ig = 1, ngm0
           v(ig,m)  = vbest(ig)
        ENDDO
#if !defined(_OPENACC) 
        !$omp end parallel do
#endif 
        aa(:,:) = 0.D0
        bb(:)   = 0.D0
        !
        CYCLE repeat_loop 
        !
     END IF
     !
     m = m + 1
     !
     IF ( lgcscf ) THEN
        !
#if defined(_OPENACC) 
        !$acc parallel loop
#else
        !$omp parallel do
#endif
        DO ig = 1, ngm0
           v(ig,m) = wbest(ig) / ( gg(ig) + agg0 + bgg0 )
        ENDDO
#if !defined(_OPENACC) 
        !$omp end parallel do
#endif
        !
     ELSE
        !
#if defined(_OPENACC)
        !$acc parallel loop
#else 
        !$omp parallel do
#endif 
        DO ig = 1, ngm0
           v(ig,m) = wbest(ig) / ( gg(ig) + agg0 )
        ENDDO
#if !defined(_OPENACC) 
        !$omp end parallel do
#endif 
        !
     END IF
     !
  END DO repeat_loop
  !
  !$acc end data
  RETURN
  !
END SUBROUTINE approx_screening2
