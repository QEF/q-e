!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE mix
  !--------------------------------------------------------------------------
  !! This module contains variables and auxiliary routines needed for
  !! the self-consistent cycle.  
  !
  USE kinds,           ONLY : DP
  USE lsda_mod,        ONLY : nspin
  USE ldaU,            ONLY : lda_plus_u, Hubbard_lmax, lda_plus_u_kind, ldmx, &
                              ldmx_b, ldmx_tot, max_num_neighbors, &
                              is_hubbard_back, orbital_resolved
  USE ions_base,       ONLY : nat
  USE buffers,         ONLY : open_buffer, close_buffer, get_buffer, save_buffer
  USE xc_lib,          ONLY : xclib_dft_is
  USE fft_base,        ONLY : dfftp
  USE fft_rho,         ONLY : rho_g2r
  USE gvect,           ONLY : ngm
  USE gvecs,           ONLY : ngms
  USE ions_base,       ONLY : ntyp => nsp
  USE paw_variables,   ONLY : okpaw
  USE uspp_param,      ONLY : nhm
  USE control_flags,   ONLY : lxdm, sic
  USE scf,             ONLY : scf_type
  !
  SAVE
  !
  TYPE mix_type
     !! For quantities directly involved in the mixing
     COMPLEX(DP), ALLOCATABLE :: of_g(:,:)
     !! the charge density in G-space
     COMPLEX(DP), ALLOCATABLE :: kin_g(:,:)
     !! the kinetic energy density in G-space
     REAL(DP),    ALLOCATABLE :: ns(:,:,:,:)
     !! the DFT+U occupation matrix
     REAL(DP),    ALLOCATABLE :: nsb(:,:,:,:)
     !! the DFT+U occupation matrix (background states)
     COMPLEX(DP), ALLOCATABLE :: ns_nc(:,:,:,:)
     !! the DFT+U occupation matrix noncollinear case 
     COMPLEX(DP), ALLOCATABLE :: nsg(:,:,:,:,:)
     !! the DFT+U+V generalized occupation matrix
     REAL(DP),    ALLOCATABLE :: bec(:,:,:)
     !! PAW corrections to hamiltonian
     REAL(DP) :: el_dipole
     !! electronic dipole, if a dipole field is present
     COMPLEX(DP), ALLOCATABLE :: pol_g(:,:)  
     !! polaron density in G-space
  END TYPE mix_type
  !
  !! DFT+U, colinear and noncolinear cases
  !! These variables are set every time create_scf_type is called
  LOGICAL :: lda_plus_u_co
  !! true if LDA+U, collinear case
  LOGICAL :: lda_plus_u_cob
  !! true if LDA+U, collinear case (background states)
  LOGICAL :: lda_plus_u_nc
  !! true if LDA+U, noncollinear case
  LOGICAL :: lda_plus_u_v
  !! true if LDA+U+V case
  LOGICAL :: need_ked
  !! true if kinetic energy density is present
  INTEGER, PRIVATE  :: record_length, &
                       rlen_rho=0,  rlen_kin=0,  rlen_ldaU=0,  rlen_bec=0,&
                       rlen_dip=0, rlen_ldaUb=0, rlen_pol=0, &
                       start_rho=0, start_kin=0, start_ldaU=0, start_bec=0, &
                       start_dipole=0, start_ldaUb=0, start_pol=0
  COMPLEX(DP), ALLOCATABLE:: io_buffer(:)
  !! buffer for storing mix_type variables
  !
  PRIVATE
  PUBLIC :: mix_rho, open_mix_file, close_mix_file
  !
CONTAINS
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
  USE ldaU,           ONLY : lda_plus_u, lda_plus_u_kind, ldim_u, neighood, &
                             max_num_neighbors, Hubbard_l, Hubbard_lmax
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
    nt,            &! index of the atomic type
    nword           ! size the DFT+U+V-related arrays
  REAL(DP), ALLOCATABLE :: betamix(:,:), work(:)
  REAL(DP), ALLOCATABLE :: nsnew(:,:,:,:)
  INTEGER, ALLOCATABLE :: iwork(:)
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
  IF (dr2 < 0.0_DP) CALL errore('mix_rho','negative dr2',1)
  !
  conv = ( dr2 < tr2 )
  !
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     conv = conv .AND. oscdft_ctx%conv
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
           CALL oscdft_nsg3 (input_rhout%nsg,nsnew)
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
     obn = 1.d0/sqrt(norm2)
     call mix_type_SCAL (obn,df(ipos))
     call mix_type_SCAL (obn,dv(ipos))
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
    END DO
    !
    DO i = 1, iter_used
        !
        gamma0 = DOT_PRODUCT( betamix(1:iter_used,i), work(1:iter_used) )
        !
        call mix_type_AXPY ( -gamma0, dv(i), rhoin_m )
        call mix_type_AXPY ( -gamma0, df(i), rhout_m )
        !
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
 !----------------------------------------------------
 SUBROUTINE create_mix_type( rho )
   !--------------------------------------------------
   !
   IMPLICIT NONE
   !
   TYPE(mix_type) :: rho
   !
   ALLOCATE( rho%of_g(ngms,nspin) )
  !$acc enter data copyin(rho) create(rho%of_g(1:ngms, 1:nspin))
   !
  !$acc kernels 
   rho%of_g = 0._dp
  !$acc end kernels
   !
   IF (need_ked) THEN
      ALLOCATE( rho%kin_g(ngms,nspin) )
     !$acc enter data create(rho%kin_g(1:ngms, 1:nspin))
     !$acc kernels
      rho%kin_g = 0._dp
     !$acc end kernels
   ENDIF
   !
   IF (lda_plus_u_nc) THEN
      ALLOCATE( rho%ns_nc(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) )
      rho%ns_nc = 0._dp
   ENDIF
   IF (lda_plus_u_co) THEN
      ALLOCATE( rho%ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat) )
      rho%ns = 0._dp
   ENDIF
   IF (lda_plus_u_cob) THEN
      ALLOCATE( rho%nsb(ldmx_b,ldmx_b,nspin,nat) )
      rho%nsb = 0._dp
   ENDIF
   IF (lda_plus_u_v ) THEN
      ALLOCATE ( rho%nsg(ldmx_tot,ldmx_tot,max_num_neighbors,nat,nspin) )
      rho%nsg = 0.0_dp
   END IF
   !
   IF (okpaw) THEN
      ALLOCATE( rho%bec(nhm*(nhm+1)/2,nat,nspin) )
      rho%bec = 0._dp
   ENDIF
   !
   rho%el_dipole = 0._dp
   !
   IF (sic) THEN
      ALLOCATE(rho%pol_g(ngms,nspin))
      rho%pol_g = 0._dp
   END IF
   !
   RETURN
   !
 END SUBROUTINE create_mix_type
 !
 !
 !------------------------------------------------------
 SUBROUTINE destroy_mix_type( rho )
   !----------------------------------------------------
   !! Deallocates a \(\text{mix_type}\) object.
   !
   IMPLICIT NONE
   !
   TYPE(mix_type) :: rho
   !
   
   IF (ALLOCATED(rho%of_g) )  THEN
    !$acc exit data finalize delete(rho%of_g) 
     DEALLOCATE( rho%of_g  )
   END IF 
   IF (ALLOCATED(rho%kin_g))  THEN
    !$acc exit data finalize delete(rho%kin_g) 
     DEALLOCATE( rho%kin_g )
   END IF
  !$acc exit data finalize delete(rho)  
   IF (ALLOCATED(rho%ns)   )  DEALLOCATE( rho%ns    )
   IF (ALLOCATED(rho%nsb)  )  DEALLOCATE( rho%nsb   )
   IF (ALLOCATED(rho%ns_nc))  DEALLOCATE( rho%ns_nc )
   IF (ALLOCATED(rho%nsg)  )  DEALLOCATE( rho%nsg   )
   IF (ALLOCATED(rho%bec)  )  DEALLOCATE( rho%bec   )
   !
   RETURN
   !
 END SUBROUTINE destroy_mix_type
 !
 !
 !-----------------------------------------------------
 SUBROUTINE assign_scf_to_mix_type( rho_s, rho_m )
   !----------------------------------------------------
   !! It fills a \(\text{mix_type}\) object starting from a
   !! \(\text{scf_type}\) one.
   !
   IMPLICIT NONE
   !
   TYPE(scf_type), INTENT(IN)  :: rho_s
   TYPE(mix_type), INTENT(INOUT) :: rho_m
   !
  !$acc enter data present_or_copyin(rho_s, rho_s%of_g) 
  !$acc kernels present(rho_m, rho_m%of_g, rho_s%of_g) 
   rho_m%of_g(1:ngms,1:nspin) = rho_s%of_g(1:ngms,1:nspin)
  !$acc end kernels 
   IF (sic) rho_m%pol_g(1:ngms,:) = rho_s%pol_g(1:ngms,:)
   !
   IF (need_ked) THEN
    !$acc enter data present_or_copyin(rho_s%kin_g)
    !$acc kernels present(rho_m%kin_g, rho_s%kin_g) 
     rho_m%kin_g(1:ngms,:) = rho_s%kin_g(1:ngms,:)
    !$acc end kernels
    !$acc exit data delete(rho_s%kin_g)
   END IF 
  !$acc exit data delete(rho_s, rho_s%of_g)  
   IF (lda_plus_u_nc)  rho_m%ns_nc  = rho_s%ns_nc
   IF (lda_plus_u_co)  rho_m%ns     = rho_s%ns
   IF (lda_plus_u_cob) rho_m%nsb    = rho_s%nsb
   IF (lda_plus_u_v)   rho_m%nsg    = rho_s%nsg
   IF (okpaw)          rho_m%bec    = rho_s%bec
   !
   rho_m%el_dipole = rho_s%el_dipole
   !
   RETURN
   !
 END SUBROUTINE assign_scf_to_mix_type
 !
 !
 !-----------------------------------------------------------------
 SUBROUTINE assign_mix_to_scf_type( rho_m, rho_s )
   !----------------------------------------------------------------
   !! It fills a \(\text{scf_type}\) object starting from a 
   !! \(\text{mix_type}\) one.
   !
   IMPLICIT NONE
   !
   TYPE(mix_type), INTENT(IN) :: rho_m
   TYPE(scf_type), INTENT(INOUT) :: rho_s
   !
   INTEGER :: is
   !   
  !$acc enter data present_or_copyin(rho_s) present_or_copyin(rho_s%of_g, rho_s%of_r)    
  !$acc kernels 
   rho_s%of_g(1:ngms,:) = rho_m%of_g(1:ngms,:)
  !$acc end kernels 
   CALL rho_g2r( dfftp, rho_s%of_g, rho_s%of_r )
  !$acc exit data copyout(rho_s%of_r, rho_s%of_g)  
   !
   IF (sic) THEN
      rho_s%pol_g(1:ngms,:) = rho_m%pol_g(1:ngms,:)
      CALL rho_g2r( dfftp, rho_s%pol_g, rho_s%pol_r )
   END IF
   !
   IF ( need_ked ) THEN
     !$acc enter data present_or_copyin(rho_s%kin_g, rho_s%kin_r)
     !$acc kernels
      rho_s%kin_g(1:ngms,:) = rho_m%kin_g(:,:)
     !$acc end kernels
      CALL rho_g2r( dfftp, rho_s%kin_g, rho_s%kin_r )
     !$acc exit data copyout(rho_s%kin_r, rho_s%kin_g) 
   ENDIF
   !
  !$acc exit data delete(rho_s) 
   IF (lda_plus_u_nc)  rho_s%ns_nc(:,:,:,:) = rho_m%ns_nc(:,:,:,:)
   IF (lda_plus_u_co)  rho_s%ns(:,:,:,:)    = rho_m%ns(:,:,:,:)
   IF (lda_plus_u_cob) rho_s%nsb(:,:,:,:)   = rho_m%nsb(:,:,:,:)
   IF (lda_plus_u_v)   rho_s%nsg(:,:,:,:,:)   = rho_m%nsg(:,:,:,:,:)
   IF (okpaw)          rho_s%bec(:,:,:)     = rho_m%bec(:,:,:)
   !
   rho_s%el_dipole = rho_m%el_dipole
   !
   RETURN
   !
 END SUBROUTINE assign_mix_to_scf_type
 !
 !
 !----------------------------------------------------------------------------
 SUBROUTINE mix_type_AXPY( A, X, Y )
  !----------------------------------------------------------------------------
  !! Works like daxpy for \(\text{scf_type}\) variables: \(Y = A\cdot X + Y\)
  ! NB: A is a REAL(DP) number
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP) :: A
  TYPE(mix_type), INTENT(IN)    :: X
  TYPE(mix_type), INTENT(INOUT) :: Y
  !
  integer :: calls = 0 
  calls = calls + 1 
 !$acc data  present(X,Y)
 !$acc kernels present(X%of_g, Y%of_g) 
  Y%of_g = Y%of_g  + A * X%of_g
 !$acc end kernels 
  !
  IF (need_ked) THEN 
   !$acc kernels present(X%kin_g, Y%kin_g)
    Y%kin_g     = Y%kin_g     + A * X%kin_g
   !$acc end kernels
  END IF 
  IF (lda_plus_u_nc)           Y%ns_nc     = Y%ns_nc     + A * X%ns_nc
  IF (lda_plus_u_co)           Y%ns        = Y%ns        + A * X%ns
  IF (lda_plus_u_cob)          Y%nsb       = Y%nsb       + A * X%nsb
  IF (lda_plus_u_v)            Y%nsg       = Y%nsg       + A * X%nsg
  IF (okpaw)                   Y%bec       = Y%bec       + A * X%bec
  IF (sic)                     Y%pol_g     = Y%pol_g     + A * X%pol_g
  ! No need to spare an operation on a single number
  ! IF (dipfield)                Y%el_dipole = Y%el_dipole + A * X%el_dipole
  Y%el_dipole = Y%el_dipole + A * X%el_dipole
  !
 !$acc end data
  RETURN
  !
 END SUBROUTINE mix_type_AXPY
 !
 !
 !----------------------------------------------------------------------------
 SUBROUTINE mix_type_COPY( X, Y )
  !----------------------------------------------------------------------------
  !! Works like DCOPY for \(\text{mix_type}\) copy variables: \(Y = X\).
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  TYPE(mix_type), INTENT(IN)    :: X
  TYPE(mix_type), INTENT(INOUT) :: Y
  !
 !$acc data present_or_copyin(Y,X)
 !$acc kernels  present_or_copyin(X%of_g, Y%of_g) 
  Y%of_g  = X%of_g
 !$acc end kernels
  !
  IF (need_ked) THEN
   !$acc kernels present_or_copyin(X%kin_g, Y%kin_g) 
    Y%kin_g     = X%kin_g
   !$acc end kernels
  END IF
  IF (lda_plus_u_nc)           Y%ns_nc     = X%ns_nc
  IF (lda_plus_u_co)           Y%ns        = X%ns
  IF (lda_plus_u_cob)          Y%nsb       = X%nsb
  IF (lda_plus_u_v)            Y%nsg       = X%nsg
  IF (okpaw)                   Y%bec       = X%bec
  IF (sic)                     Y%pol_g     = X%pol_g
  Y%el_dipole = X%el_dipole
  !
 !$acc end data
  RETURN
  !
 END SUBROUTINE mix_type_COPY
 !
 !
 !----------------------------------------------------------------------------
 SUBROUTINE mix_type_SCAL( A, X )
  !----------------------------------------------------------------------------
  !! Works like DSCAL for \(\text{mix_type}\) copy variables: \(X = A \cdot X\)  
  !! NB: A is a REAL(DP) number
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  !
  REAL(DP),       INTENT(IN)    :: A
  TYPE(mix_type), INTENT(INOUT) :: X
  !
  !
 !$acc data present_or_copyin(X)
 !$acc kernels present_or_copyin(X%of_g) 
  X%of_g(:,:) = A * X%of_g(:,:)
 !$acc end kernels
  !
  IF (need_ked) THEN
   !$acc kernels present_or_copyin(X%kin_g)
    X%kin_g     = A * X%kin_g
   !$acc end kernels
  END IF 
  IF (lda_plus_u_nc)           X%ns_nc     = A * X%ns_nc
  IF (lda_plus_u_co)           X%ns        = A * X%ns
  IF (lda_plus_u_cob)          X%nsb       = A * X%nsb
  IF (lda_plus_u_v)            X%nsg       = A * X%nsg
  IF (okpaw)                   X%bec       = A * X%bec
  IF (sic)                     X%pol_g     = A * X%pol_g
  X%el_dipole = A * X%el_dipole
  !
 !$acc end data
  RETURN
  !
 END SUBROUTINE mix_type_SCAL
 !
 !
 !---------------------------------------------------------------------
 SUBROUTINE high_frequency_mixing( rhoin, input_rhout, alphamix )
   !-------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   TYPE (scf_type), INTENT(INOUT) :: rhoin
   TYPE (scf_type), INTENT(IN) :: input_rhout
   REAL(DP), INTENT(IN) :: alphamix
   !
   ! ... local variable
   !
   INTEGER :: is

   call start_clock('high_freq_mix') 
   !
   !$acc data present_or_copyin(rhoin, rhoin%of_g, rhoin%of_r) 
   IF (ngms < ngm ) THEN
      !
      rhoin%of_g = rhoin%of_g + alphamix * (input_rhout%of_g-rhoin%of_g)
      rhoin%of_g(1:ngms,1:nspin) = (0.d0,0.d0)
      CALL rho_g2r( dfftp, rhoin%of_g, rhoin%of_r )
      !
      IF (need_ked) THEN
         rhoin%kin_g = rhoin%kin_g + alphamix * ( input_rhout%kin_g-rhoin%kin_g)
         rhoin%kin_g(1:ngms,1:nspin) = (0.d0,0.d0)
         CALL rho_g2r( dfftp, rhoin%kin_g, rhoin%kin_r )
      ENDIF
      !
      IF(sic) THEN
         rhoin%pol_g = rhoin%pol_g + alphamix * (input_rhout%pol_g-rhoin%pol_g)
         rhoin%pol_g(1:ngms,1:nspin) = (0.d0,0.d0)
         CALL rho_g2r( dfftp, rhoin%pol_g, rhoin%pol_r )
      END IF
      !
   ELSE
      !
      rhoin%of_g(:,:)= (0.d0,0.d0)
      rhoin%of_r(:,:)= 0.d0
      IF (need_ked) THEN
         rhoin%kin_g(:,:)= (0.d0,0.d0)
         rhoin%kin_r(:,:)= 0.d0
      ENDIF
      IF(sic) then
         rhoin%pol_g(:,:)= (0.d0,0.d0)
         rhoin%pol_r(:,:)= 0.d0
      END IF
      !
   ENDIF
   !
   IF (lda_plus_u_nc)  rhoin%ns_nc(:,:,:,:) = 0.d0
   IF (lda_plus_u_co)  rhoin%ns(:,:,:,:)    = 0.d0
   IF (lda_plus_u_cob) rhoin%nsb(:,:,:,:)   = 0.d0
   IF (lda_plus_u_v)   rhoin%nsg(:,:,:,:,:) = 0.d0
   !
   !$acc end data 
   call stop_clock('high_freq_mix') 
   RETURN
   !
 END SUBROUTINE high_frequency_mixing 
 !
 !
 !------------------------------------------------------------------------
 SUBROUTINE open_mix_file( iunit, extension, exst )
   !------------------------------------------------------------------------
   !
   !! open_mix_files performs some initializations, must be called first
   !
   USE control_flags,  ONLY : io_level
   !
   IMPLICIT NONE
   !
   CHARACTER(LEN=*), INTENT(IN) :: extension
   INTEGER, INTENT(IN) :: iunit
   LOGICAL :: exst
   !
   lda_plus_u_co = lda_plus_u .AND. .NOT. (nspin == 4 ) .AND. .NOT. ( lda_plus_u_kind == 2)
   lda_plus_u_nc = lda_plus_u .AND.       (nspin == 4 ) .AND. .NOT. ( lda_plus_u_kind == 2)
   lda_plus_u_cob = lda_plus_u_co .AND. ANY( is_hubbard_back(1:ntyp) )
   lda_plus_u_v   =  ( lda_plus_u_kind == 2 )
   !
   need_ked = xclib_dft_is('meta') .OR. lxdm
   ! define lengths (in real numbers) of different record chunks
   !
   rlen_rho = 2 * ngms * nspin
   IF (need_ked)                rlen_kin  = 2 * ngms * nspin 
   IF (lda_plus_u_co)           rlen_ldaU = (2*Hubbard_lmax+1)**2 *nspin*nat
   IF (lda_plus_u_cob)          rlen_ldaUb = (ldmx_b)**2 *nspin*nat
   IF (lda_plus_u_nc)           rlen_ldaU = 2 * (2*Hubbard_lmax+1)**2 *nspin*nat
   IF (lda_plus_u_v)            rlen_ldaU = 2 * (ldmx_tot**2*max_num_neighbors*nat*nspin)
   IF (okpaw)                   rlen_bec  = (nhm*(nhm+1)/2) * nat * nspin
   IF (sic)                     rlen_pol  = 2*ngms*nspin
   rlen_dip  = 1
   !
   ! define the starting point of the different chunks. Beware: each starting point
   ! is the index of a COMPLEX array. When real arrays with odd dimension are copied
   ! to/from the complex array io_buffer, the last complex number will be half-filled
   ! but must still be counted as one!
   start_rho    = 1
   start_kin    = start_rho  + rlen_rho / 2
   start_ldaU   = start_kin  + rlen_kin / 2
   IF (lda_plus_u_cob) THEN
      start_ldaUb = start_ldaU + ( rlen_ldaU + 1 ) / 2
      start_bec = start_ldaUb + ( rlen_ldaUb + 1 ) / 2
   ELSE
      ! FIXME: not always present
      start_bec = start_ldaU + ( rlen_ldaU + 1 ) / 2
   ENDIF
   start_dipole = start_bec  + ( rlen_bec + 1 ) / 2
   start_pol    = start_dipole + ( rlen_dip + 1 ) / 2
   !
   ! define total record length, in complex numbers
   record_length = start_pol + rlen_pol - 1
   !
   ! open file and allocate io_buffer
   CALL open_buffer( iunit, extension, record_length, io_level, exst )
   !
   ALLOCATE( io_buffer(record_length) )
   ! setting to zero -prevents trouble with "holes" due to odd dimensions of real
   ! arrays
   io_buffer(:) = (0.0_dp, 0.0_dp)
   !
   RETURN
   !
 END SUBROUTINE open_mix_file
 !
 !
 !------------------------------------------------------------------------
 SUBROUTINE close_mix_file( iunit, stat )
   !---------------------------------------------------------------------
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: iunit
   CHARACTER(LEN=*), INTENT(IN) :: stat
   !
   DEALLOCATE( io_buffer )
   !
   CALL close_buffer( iunit, TRIM(stat) ) 
   !
   RETURN
   !
 END SUBROUTINE close_mix_file
 !
 !
 !------------------------------------------------------------
 SUBROUTINE davcio_mix_type( rho, iunit, record, iflag )
   !----------------------------------------------------------
   !
   IMPLICIT NONE
   !
   TYPE(mix_type) :: rho
   INTEGER, INTENT(IN) :: iunit, record, iflag
   !
   IF (iflag > 0) THEN
      !
     !$acc update self(rho%of_g) 
      CALL DCOPY(rlen_rho,rho%of_g,1,io_buffer(start_rho),1)
      !
      IF (need_ked) THEN
       !$acc update self(rho%kin_g) 
        CALL DCOPY(rlen_kin, rho%kin_g,1,io_buffer(start_kin), 1)
      END IF 
      IF (lda_plus_u_nc)           CALL DCOPY(rlen_ldaU,rho%ns_nc,1,io_buffer(start_ldaU),1)
      IF (lda_plus_u_co)           CALL DCOPY(rlen_ldaU,rho%ns,   1,io_buffer(start_ldaU),1)
      IF (lda_plus_u_v)            CALL DCOPY(rlen_ldaU,rho%nsg,1,io_buffer(start_ldaU),1)
      IF (lda_plus_u_cob)          CALL DCOPY(rlen_ldaUb,rho%nsb, 1,io_buffer(start_ldaUb),1)
      IF (okpaw)                   CALL DCOPY(rlen_bec, rho%bec,  1,io_buffer(start_bec), 1)
      !
      io_buffer(start_dipole) = CMPLX( rho%el_dipole, 0.0_dp, KIND=DP )
      IF (sic)                     CALL DCOPY(rlen_pol, rho%pol_g, 1,io_buffer(start_pol),1)
      !
      CALL save_buffer( io_buffer, record_length, iunit, record )   
      !
   ELSEIF (iflag < 0 ) THEN
      !
      CALL get_buffer( io_buffer, record_length, iunit, record )
      !
      CALL DCOPY(rlen_rho,io_buffer(start_rho),1,rho%of_g,1)
     !$acc update device(rho%of_g)
      !
      IF (need_ked) THEN 
        CALL DCOPY(rlen_kin, io_buffer(start_kin), 1,rho%kin_g,1)
       !$acc update device(rho%kin_g) 
      END IF
      IF (lda_plus_u_co)           CALL DCOPY(rlen_ldaU,io_buffer(start_ldaU),1,rho%ns,   1)
      IF (lda_plus_u_cob)          CALL DCOPY(rlen_ldaUb,io_buffer(start_ldaUb),1,rho%nsb,1)
      IF (lda_plus_u_nc)           CALL DCOPY(rlen_ldaU,io_buffer(start_ldaU),1,rho%ns_nc,1)
      IF (lda_plus_u_v)            CALL DCOPY(rlen_ldaU,io_buffer(start_ldaU),1,rho%nsg,1)
      IF (okpaw)                   CALL DCOPY(rlen_bec, io_buffer(start_bec), 1,rho%bec,  1)
      !
      rho%el_dipole = DBLE( io_buffer(start_dipole) )
      IF (sic)                     CALL DCOPY(rlen_pol, io_buffer(start_pol), 1,rho%pol_g,1)
      !
   ENDIF
   !
 END SUBROUTINE davcio_mix_type
 !
 !
 !-----------------------------------------------------------------------------------
FUNCTION rho_ddot( rho1, rho2, gf, g0 )
  !----------------------------------------------------------------------------------
  !! Calculates \(4\pi/G^2\ \rho_1(-G)\ \rho_2(G) = V1_\text{Hartree}(-G)\ \rho_2(G)\)
  !! used as an estimate of the self-consistency error on the energy.
  !
  USE kinds,           ONLY : DP
  USE constants,       ONLY : e2, tpi, fpi
  USE cell_base,       ONLY : omega, tpiba2
  USE gvect,           ONLY : gg, gstart
  USE control_flags,   ONLY : gamma_only
  USE noncollin_module,ONLY : noncolin
  USE paw_onecenter,   ONLY : paw_ddot
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp,              ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE(mix_type), INTENT(IN) :: rho1
  !! first density matrix
  TYPE(mix_type), INTENT(IN) :: rho2
  !! second density matrix
  INTEGER, INTENT(IN) :: gf
  !! points delimiter
  REAL(DP), OPTIONAL, INTENT(IN) :: g0
  !! factorized G-vector norm of G=0 used in GC-SCF
  REAL(DP) :: rho_ddot
  !! output: see function comments
  !
 !$acc declare present(rho1, rho2)
  ! ... local variables
  !
  REAL(DP) :: fac
  REAL(DP) :: gg0
  REAL(DP) :: rho0
  INTEGER  :: ig
  !
  fac = e2 * fpi / tpiba2
  !
  rho_ddot = 0.D0
  !
  IF ( PRESENT(g0) ) THEN
     !
     gg0 = g0 * g0 / tpiba2
     !
  ELSE
     !
     gg0 = -1.0_DP
     !
  END IF
  !
  IF ( gg0 > 0.0_DP ) THEN
     !
    !$acc parallel loop reduction(+:rho_ddot)
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + &
                   REAL( CONJG( rho1%of_g(ig,1) )*rho2%of_g(ig,1), DP ) / ( gg(ig) + gg0 )
        !
     END DO
    !$acc end parallel loop 
     !
     IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
     !
     IF ( gstart == 2 ) THEN
        !
       !$acc update host(rho1%of_g(1,1), rho2%of_g(1,1))
        rho_ddot = rho_ddot + &
                   REAL( CONJG( rho1%of_g(1,1) )*rho2%of_g(1,1), DP ) / ( gg(1) + gg0 )
        !
     END IF
     !
  ELSE
     !
    !$acc parallel loop reduction(+:rho_ddot)
     DO ig = gstart, gf
        !
        rho_ddot = rho_ddot + &
                   REAL( CONJG( rho1%of_g(ig,1) )*rho2%of_g(ig,1), DP ) / gg(ig)
        !
     END DO
    !$acc end parallel loop
     !
     IF ( gamma_only ) rho_ddot = 2.D0 * rho_ddot
     !
  END IF
  !
  rho_ddot = fac*rho_ddot
  !
  IF ( nspin >= 2 )  THEN
     fac = e2*fpi / tpi**2  ! lambda=1 a.u.
     IF ( gstart == 2 ) THEN
        !$acc update host(rho1%of_g(1,2:nspin), rho2%of_g(1,2:nspin))
        rho_ddot = rho_ddot + &
                fac * SUM(REAL(CONJG( rho1%of_g(1,2:nspin))*(rho2%of_g(1,2:nspin) ), DP))
     ENDIF
     !
     IF ( gamma_only ) fac = 2.D0 * fac
     !
    !$acc parallel loop reduction(+:rho_ddot)
     DO ig = gstart, gf
        rho_ddot = rho_ddot + &
              fac * SUM(REAL(CONJG( rho1%of_g(ig,2:nspin))*(rho2%of_g(ig,2:nspin) ), DP))
     ENDDO
    !$acc end parallel do
  ENDIF
  !
  rho_ddot = rho_ddot * omega * 0.5D0
  !
  CALL mp_sum( rho_ddot, intra_bgrp_comm )
  !
  IF (xclib_dft_is('meta')) rho_ddot = rho_ddot + tauk_ddot( rho1, rho2, gf )
  !
  IF (lda_plus_u ) THEN 
      IF ( orbital_resolved ) THEN
         IF ( noncolin ) THEN
            rho_ddot = rho_ddot + ns_ddot_um_nc( rho1, rho2 )
         ELSE
            rho_ddot = rho_ddot + ns_ddot_um( rho1, rho2 )
         ENDIF
      ELSE
         rho_ddot = rho_ddot + ns_ddot( rho1, rho2 )
      ENDIF
  ENDIF
  ! 
  ! Beware: paw_ddot has a hidden parallelization on all processors
  !         it must be called on all processors or else it will hang
  ! Beware: commented out because it yields too often negative values
  ! IF (okpaw)         rho_ddot = rho_ddot + paw_ddot(rho1%bec, rho2%bec)
  !
  ! Dipole is zero if not set
  ! IF (dipfield)
  rho_ddot = rho_ddot + (e2/2.0_DP)*(rho1%el_dipole * rho2%el_dipole)*omega/fpi
  !
  RETURN
  !
END FUNCTION rho_ddot
!
!
!----------------------------------------------------------------------------
FUNCTION tauk_ddot( rho1, rho2, gf )
  !----------------------------------------------------------------------------
  !! Calculates \(4\pi/G^2\ \rho_1(-G)\ \rho_2(G) = V1_\text{Hartree}(-G)\ \rho_2(G)\)
  !! used as an estimate of the self-consistency error on the energy - kinetic density
  !! version.
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : e2, tpi, fpi
  USE cell_base,     ONLY : omega, tpiba2
  USE gvect,         ONLY : gg, gstart
  USE control_flags, ONLY : gamma_only
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  TYPE(mix_type), INTENT(IN) :: rho1
  !! first kinetic density
 !$acc declare present(rho1)
  TYPE(mix_type), INTENT(IN) :: rho2
  !! second kinetic density
 !$acc declare present(rho2)
  INTEGER, INTENT(IN) :: gf
  !! point delimiter
  REAL(DP) :: tauk_ddot
  !! output: see function comments
  !
  ! ... local variables
  !
  REAL(DP) :: fac
  INTEGER :: ig
  !
  tauk_ddot = 0.D0
  !
  !  write (*,*) rho1%kin_g(1:4,1)
  !  if (.true. ) stop
  !
 !$acc parallel loop reduction(+:tauk_ddot)
  DO ig = gstart, gf
     tauk_ddot = tauk_ddot + DBLE( CONJG( rho1%kin_g(ig,1) )*rho2%kin_g(ig,1) ) 
  ENDDO
  !
  IF ( nspin==1 .AND. gamma_only ) tauk_ddot = 2.D0 * tauk_ddot
  !
  ! ... G=0 term
  !
  IF ( gstart == 2 ) THEN
    !$acc update host(rho1%kin_g(1,1:nspin), rho2%kin_g(1,1:nspin))
     tauk_ddot = tauk_ddot + DBLE( CONJG( rho1%kin_g(1,1) ) * rho2%kin_g(1,1) )
  ENDIF
  !
  IF ( nspin >= 2 ) THEN
    !$acc parallel loop reduction (+: tauk_ddot)
     DO ig = gstart, gf
        tauk_ddot = tauk_ddot + &
          SUM( REAL( CONJG( rho1%kin_g(ig,2:nspin))*(rho2%kin_g(ig,2:nspin) ), DP))
     ENDDO
    !$acc end parallel loop
     !
     IF ( gamma_only ) tauk_ddot = 2.D0 * tauk_ddot
     !
     ! ... G=0 term
     IF ( gstart == 2 ) THEN
        tauk_ddot = tauk_ddot + &
          SUM(REAL(CONJG( rho1%kin_g(1,2:nspin))*(rho2%kin_g(1,2:nspin) ), DP))
     ENDIF
     !
     IF ( nspin == 2 ) tauk_ddot = 0.5D0 *  tauk_ddot 
  ENDIF
  !
  fac = e2 * fpi / tpi**2  ! lambda = 1 a.u.
  !
  tauk_ddot = fac * tauk_ddot * omega * 0.5D0
  !
  CALL mp_sum( tauk_ddot, intra_bgrp_comm )
  !
  RETURN
  !
END FUNCTION tauk_ddot
!
!
!----------------------------------------------------------------------------
FUNCTION ns_ddot( rho1, rho2 )
  !---------------------------------------------------------------------------
  !! Calculates \(U/2 \sum_i \text{ns1}(i)\ \text{ns2}(i)\) used as an estimate
  !! of the self-consistency error on the DFT+U correction to the energy.
  !
  USE kinds,     ONLY : DP
  USE ldaU,      ONLY : Hubbard_l, Hubbard_U, Hubbard_U2, ldim_back, &
                        lda_plus_u_kind, is_hubbard, is_hubbard_back
  USE ions_base, ONLY : nat, ityp
  !
  IMPLICIT NONE  
  !
  TYPE(mix_type), INTENT(IN) :: rho1
  !! first Hubbard ns
  TYPE(mix_type), INTENT(IN) :: rho2
  !! second Hubbard ns
  REAL(DP) :: ns_ddot
  !! output: see function comments
  !
  ! ... local variables
  !
  INTEGER :: na, nt, m1, m2
  !
  ns_ddot = 0.D0
  !
  IF (lda_plus_u_kind.EQ.2) THEN
     ns_ddot = nsg_ddot( rho1%nsg, rho2%nsg, nspin )
     RETURN
  END IF
  !
  DO na = 1, nat
     nt = ityp(na)
     IF ( is_hubbard(nt) ) THEN
        !
        m1 = 2 * Hubbard_l(nt) + 1
        m2 = 2 * Hubbard_l(nt) + 1
        !
        IF (nspin == 4) THEN
          ns_ddot = ns_ddot + 0.5D0 * Hubbard_U(nt) * &
                  SUM( CONJG(rho1%ns_nc(:m1,:m2,:nspin,na))*rho2%ns_nc(:m1,:m2,:nspin,na) )
        ELSE
          ns_ddot = ns_ddot + 0.5D0 * Hubbard_U(nt) * &
                  SUM( rho1%ns(:m1,:m2,:nspin,na)*rho2%ns(:m1,:m2,:nspin,na) )
        ENDIF
        !
     ENDIF
     !
     ! Background part
     !
     IF ( is_hubbard_back(nt) .AND. lda_plus_u_kind.EQ.0 ) THEN
        !
        m1 = ldim_back(nt)
        m2 = ldim_back(nt)
        !
        ns_ddot = ns_ddot + 0.5D0 * Hubbard_U2(nt) * &
                SUM( rho1%nsb(:m1,:m2,:nspin,na)*rho2%nsb(:m1,:m2,:nspin,na) )
        !
     ENDIF
     !
  ENDDO
  !
  IF ( nspin == 1 ) ns_ddot = 2.D0*ns_ddot
  !
  RETURN
  !
END FUNCTION ns_ddot
!
!----------------------------------------------------------------------------
FUNCTION nsg_ddot( nsg1, nsg2, nspin )
  !----------------------------------------------------------------------------
  !! Calculates \(U/2 \sum_i \text{nsg1}(i) \text{nsg2}(i)\)
  !! used as an estimate of the self-consistency error on the 
  !! DFT+U+V correction to the energy
  !
  USE kinds,     ONLY : DP
  USE ldaU,      ONLY : Hubbard_V, max_num_neighbors, is_hubbard, &
                        is_hubbard_back, ldim_u, neighood, ldmx,  &
                        ldmx_tot, at_sc
  USE ions_base, ONLY : nat, ityp
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: nsg1(ldmx_tot,ldmx_tot,max_num_neighbors,nat,nspin), &
                             nsg2(ldmx_tot,ldmx_tot,max_num_neighbors,nat,nspin)
  INTEGER,  INTENT(IN) :: nspin
  REAL(DP)             :: nsg_ddot
  !
  INTEGER :: na1, nt1, m1, m2, na2, nt2, viz, equiv_na2, i_type
  INTEGER, EXTERNAL :: find_viz, type_interaction
  !
  nsg_ddot = 0.D0
  !
  DO na1 = 1, nat
     nt1 = ityp(na1)
     IF ( (is_hubbard(nt1) .OR. is_hubbard_back(nt1)) .AND.ldim_u(nt1).GT.0 ) THEN
        DO viz = 1, neighood(na1)%num_neigh
           na2 = neighood(na1)%neigh(viz)
           equiv_na2 = at_sc(na2)%at
           nt2 = ityp(equiv_na2)
           IF ((Hubbard_V(na1,na2,1) .NE. 0.d0) .OR. &
               (Hubbard_V(na1,na2,2) .NE. 0.d0) .OR. &
               (Hubbard_V(na1,na2,3) .NE. 0.d0) .OR. &
               (Hubbard_V(na1,na2,4) .NE. 0.d0) ) THEN
              DO m1=1,ldim_u(nt1)
                 DO m2=1,ldim_u(nt2)
                    i_type = type_interaction(na1,m1,equiv_na2,m2)
                    nsg_ddot = nsg_ddot + 0.5D0 * ABS(Hubbard_V(na1,na2,i_type)) * &
                      REAL(SUM( nsg1(m2,m1,viz, na1,:nspin)*  &
                                CONJG(nsg2(m2,m1,viz, na1,:nspin) ) ) )
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDIF
  ENDDO
  !
  IF ( nspin == 1 ) nsg_ddot = 2.D0*nsg_ddot
  !
  RETURN
  !
END FUNCTION nsg_ddot
!
!
!----------------------------------------------------------------------------
FUNCTION ns_ddot_um( rho1, rho2 )
  !---------------------------------------------------------------------------
  !! Calculates \(U/2 \sum_i \text{ns1}(i)\ \text{ns2}(i)\) used as an estimate
  !! of the self-consistency error on the orbital-resolved DFT+U correction to the energy.
  !
  USE kinds,     ONLY : DP
  USE ldaU,      ONLY : Hubbard_l, Hubbard_U, Hubbard_U2, ldim_back, &
                        lda_plus_u_kind, is_hubbard, eigenvecs_ref, &
                        Hubbard_lmax, Hubbard_Um, apply_U
  USE ions_base, ONLY : nat, ityp
  USE constants, ONLY : eps16, RYTOEV
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE  
  !
  TYPE(mix_type), INTENT(IN) :: rho1
  !! first Hubbard ns
  TYPE(mix_type), INTENT(IN) :: rho2
  !! second Hubbard ns
  REAL(DP) :: ns_ddot_um
  !! output: see function comments
  !
  ! ... local variables
  !
  COMPLEX(DP)  :: vet1(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin)
  COMPLEX(DP)  :: vet2(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin)
  INTEGER      :: order1(2*Hubbard_lmax+1), order2(2*Hubbard_lmax+1)
  INTEGER      :: na, nt, ldim, is, m, index1, index2
  REAL(DP)     :: lambda1(2*Hubbard_lmax+1,nspin), lambda2(2*Hubbard_lmax+1,nspin)
  !
  ns_ddot_um = 0.D0
  !
  IF (.NOT. apply_U) RETURN
  ! if apply_um is still .FALSE.
  ! do not (yet) apply Hubbard U corrections.
  !
  DO na = 1, nat
    nt = ityp(na)
    IF ( is_hubbard(nt) ) THEN
      !
      ldim = 2 * Hubbard_l(nt) + 1
      !
      vet1(:,:,:) = CMPLX(0.D0,0.D0, kind=dp)
      vet2(:,:,:) = CMPLX(0.D0,0.D0, kind=dp)
      lambda1(:,:) = 0.D0
      lambda2(:,:) = 0.D0
      !
      ! diagonalize old- and new occupation matrix
      CALL diag_ns( ldim, rho1%ns(1:ldim,1:ldim,:,na), lambda1(1:ldim,:), vet1(1:ldim,1:ldim,:) )
      CALL diag_ns( ldim, rho2%ns(1:ldim,1:ldim,:,na), lambda2(1:ldim,:), vet2(1:ldim,1:ldim,:) )
      ! 
      DO is = 1, nspin
         ! order eigenvectors
         order1(:) = 0
         order2(:) = 0
         CALL order_eigenvecs( order1(1:ldim), vet1(1:ldim,1:ldim,is), &
                                 eigenvecs_ref(1:ldim,1:ldim,is,na), ldim )
         CALL order_eigenvecs( order2(1:ldim), vet2(1:ldim,1:ldim,is), &
                                 eigenvecs_ref(1:ldim,1:ldim,is,na), ldim )
         !
         IF ( ANY(ABS(Hubbard_Um(:,is,nt)) .GT. eps16) ) THEN
            ! compute U(m)/2*SUM(ns1*ns2)
            DO m = 1, ldim
               !
               ! find the index where the order vector is
               ! equal to m to apply the same Hubbard_Um
               ! to the same eigenstates 
               index1 = FINDLOC(order1,m,dim=1)
               index2 = FINDLOC(order2,m,dim=1)
               !
               ns_ddot_um = ns_ddot_um + 0.5D0 * Hubbard_Um(m,is,nt) * &
                      lambda1(index1,is) * lambda2(index2,is)
               !
               ! This can be removed once the merge request is approved
#if defined(__DEBUG)
               WRITE(stdout,'(5X,"m: ", i1,", is: ", i1, ", index1:", i1, " ,&
                      & index2:", i1)') m, is,index1,index2
               WRITE(stdout,'(5X,"U: ", f5.3,", lambda1: ", f7.4,", lambda2: ", &
                      & f7.4,", ns_ddot_um: ", f7.4)') &
                      Hubbard_Um(m,is,nt)*RYTOEV,lambda1(index1,is),lambda2(index2,is),ns_ddot_um
#endif
            !
            ENDDO
            !
         ENDIF
         !
      ENDDO
      !
     ENDIF
     !
  ENDDO
  !
  IF ( nspin == 1 ) ns_ddot_um = 2.D0*ns_ddot_um
  !
  RETURN
  !
END FUNCTION ns_ddot_um
!
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
FUNCTION ns_ddot_um_nc( rho1, rho2 )
   !---------------------------------------------------------------------------
   !! Calculates \(U/2 \sum_i \text{ns1}(i)\ \text{ns2}(i)\) used as an estimate
   !! of the self-consistency error on the orbital-resolved DFT+U correction to the energy.
   !
   USE kinds,     ONLY : DP
   USE ldaU,      ONLY : Hubbard_l, ldim_back, &
                         lda_plus_u_kind, is_hubbard, eigenvecs_ref, &
                         Hubbard_lmax, Hubbard_Um_nc, apply_U
   USE ions_base, ONLY : nat, ityp
   USE constants, ONLY : eps16, RYTOEV
   USE io_global, ONLY : stdout
   !
   IMPLICIT NONE  
   !
   TYPE(mix_type), INTENT(IN) :: rho1
   !! first Hubbard ns
   TYPE(mix_type), INTENT(IN) :: rho2
   !! second Hubbard ns
   REAL(DP) :: ns_ddot_um_nc
   !! output: see function comments
   !
   ! ... local variables
   !
   ! For NC case we allocate arrays as 2*(2l+1)
   COMPLEX(DP)  :: vet1(4*Hubbard_lmax+2,4*Hubbard_lmax+2)
   COMPLEX(DP)  :: vet2(4*Hubbard_lmax+2,4*Hubbard_lmax+2)
   REAL(DP)     :: lambda1(4*Hubbard_lmax+2), lambda2(4*Hubbard_lmax+2)
   INTEGER      :: order1(4*Hubbard_lmax+2), order2(4*Hubbard_lmax+2)
   INTEGER      :: na, nt, ldim, is, m, index1, index2

   !
   ns_ddot_um_nc = 0.D0
   !
   IF (.NOT. apply_U) RETURN
   ! if apply_um is still .FALSE.
   ! do not (yet) apply Hubbard U corrections.
   !
   DO na = 1, nat
     nt = ityp(na)
     IF ( is_hubbard(nt) ) THEN
       !
       ldim = 2 * Hubbard_l(nt) + 1
       !
       vet1(:,:) = CMPLX(0.D0,0.D0, kind=dp)
       vet2(:,:) = CMPLX(0.D0,0.D0, kind=dp)
       lambda1(:) = 0.D0
       lambda2(:) = 0.D0
       !
       ! diagonalize old- and new occupation matrix
       CALL diag_ns_nc( ldim, rho1%ns(1:ldim,1:ldim,:,na), lambda1(1:2*ldim), vet1(1:2*ldim,1:2*ldim) )
       CALL diag_ns_nc( ldim, rho2%ns(1:ldim,1:ldim,:,na), lambda2(1:2*ldim), vet2(1:2*ldim,1:2*ldim) )
       ! 
       ! order eigenvectors
       order1(:) = 0
       order2(:) = 0
       CALL order_eigenvecs( order1(1:2*ldim), vet1(1:2*ldim,1:2*ldim), &
                               eigenvecs_ref(1:2*ldim,1:2*ldim,1,na), 2*ldim )
       CALL order_eigenvecs( order2(1:2*ldim), vet2(1:2*ldim,1:2*ldim), &
                               eigenvecs_ref(1:2*ldim,1:2*ldim,1,na), 2*ldim )
       !
       IF ( ANY(ABS(Hubbard_Um_nc(:,nt)) .GT. eps16) ) THEN
          ! compute U(m)/2*SUM(ns1*ns2)
          DO m = 1, 2*ldim
             !
             ! find the index where the order vector is
             ! equal to m to apply the same Hubbard_Um
             ! to the same eigenstates 
             index1 = FINDLOC(order1,m,dim=1)
             index2 = FINDLOC(order2,m,dim=1)
             !
             ns_ddot_um_nc = ns_ddot_um_nc + 0.5D0 * Hubbard_Um_nc(m,nt) * &
                           lambda1(index1) * lambda2(index2)
             !
             ! This can be removed once the merge request is approved
#if defined(__DEBUG)
             WRITE(stdout,'(5X,"m: ", i1,", is: ", i1, ", index1:", i1, " , &
                   & index2:", i1)') m, is,index1,index2
             WRITE(stdout,'(5X,"U: ", f5.3,", lambda1: ", f7.4,", lambda2: ",&
                   & f7.4,", ns_ddot_um_nc: ", f7.4)') &
                   Hubbard_Um_nc(m,nt)*RYTOEV,lambda1(index1),lambda2(index2),ns_ddot_um_nc
#endif
          !
          ENDDO
       !
       ENDIF
      !
      ENDIF
      !
   ENDDO
   !
RETURN
   !
 END FUNCTION ns_ddot_um_nc
 !
 !----------------------------------------------------------------------------
FUNCTION local_tf_ddot( rho1, rho2, ngm0, g0 )
  !----------------------------------------------------------------------------
  !! Calculates \(4\pi/G^2\ \rho_1(-G)\ \rho_2(G) = V1_\text{Hartree}(-G)\ \rho_2(G)\)
  !! used as an estimate of the self-consistency error on the energy - version 
  !! for the case with local-density dependent TF preconditioning to drho.
  !
  USE kinds,           ONLY : DP
  USE constants,       ONLY : e2, fpi
  USE cell_base,       ONLY : omega, tpiba2
  USE gvect,           ONLY : gg, gstart
  USE control_flags,   ONLY : gamma_only
  USE mp_bands,        ONLY : intra_bgrp_comm
  USE mp,              ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ngm0
  !! input length (local number of smooth vectors)
  COMPLEX(DP), INTENT(IN) :: rho1(ngm0)
  !! see main comment
  COMPLEX(DP), INTENT(IN) :: rho2(ngm0)
  !! see main comment
  REAL(DP), OPTIONAL, INTENT(IN) :: g0
  !! factorized G-vector norm of G=0 used in GC-SCF
  REAL(DP) :: local_tf_ddot
  !! see main comment
  !
  ! ... local variables
  !
  REAL(DP) :: fac
  REAL(DP) :: gg0
  INTEGER  :: ig
  ! 
  !
  fac = e2 * fpi / tpiba2
  !
  IF ( PRESENT(g0) ) THEN
     !
     gg0 = g0 * g0 / tpiba2
     !
  ELSE
     !
     gg0 = 0.0_DP
     !
  END IF
  !
  local_tf_ddot = 0.D0
  !$acc data present_or_copyin(rho1, rho2)
#if defined(_OPENACC)
  !$acc parallel loop reduction(+:local_tf_ddot) 
#else
  !$omp parallel do reduction(+:local_tf_ddot)
#endif
  DO ig = gstart, ngm0
     local_tf_ddot = local_tf_ddot + DBLE( CONJG(rho1(ig))*rho2(ig) ) / ( gg(ig) + gg0 )
  END DO
#if !defined(_OPENACC) 
  !$omp end parallel do
#endif
  !$acc end data
  !
  IF ( gamma_only ) local_tf_ddot = 2.D0 * local_tf_ddot
  IF ( gstart == 2 .AND. gg0 > 0.0_dp ) THEN
     ! This is the G=0 term, that for gamma_only must not be counted twice
     local_tf_ddot = local_tf_ddot + DBLE( CONJG(rho1(1))*rho2(1) ) / ( gg(1) + gg0 )
  END IF
  local_tf_ddot = fac * local_tf_ddot * omega * 0.5D0
  CALL mp_sum( local_tf_ddot, intra_bgrp_comm )
  !
  RETURN
  !
END FUNCTION local_tf_ddot
!
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
!
END MODULE mix
