!
! Copyright (C) 2001-2014 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
PROGRAM benchmark_libxc
  !
  !------------------------------------------------------------------------------------!
  !  REMEMBER to comment eventual libxc blocks in the functional routines in 'Modules' !
  !  folder in order to run consistent tests (in the present version they should be    !
  !  already absent, however).                                                         !
  !------------------------------------------------------------------------------------!
  !
#if defined(__LIBXC)
  !
  USE xc_f90_types_m
  USE xc_f90_lib_m
  !
  USE xc_lda_lsda
  !
  IMPLICIT NONE
  !
  !-------- Common vars ----------------------
  INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(14,200)
  INTEGER, PARAMETER :: nnr = 6
  CHARACTER(LEN=120) :: e_q, f_q
  INTEGER :: ii, ns, np, quit, i_sub
  LOGICAL :: LDA, POLARIZED, ENERGY_ONLY, DF_OK
  REAL(DP), PARAMETER :: null = 0.0_DP, pi34 = 0.6203504908994_DP
  !
  !----------QE vars --------------------------
  INTEGER :: iexch_qe, icorr_qe
  REAL(DP) :: rs(nnr)
  REAL(DP), ALLOCATABLE :: rho_qe(:,:)
  REAL(DP), ALLOCATABLE :: rho_tot(:), zeta(:)
  REAL(DP), ALLOCATABLE :: ex_qe(:), ec_qe(:)
  REAL(DP), ALLOCATABLE :: vx_qe(:,:), vc_qe(:,:)
  REAL(DP), ALLOCATABLE :: dmuxc(:,:,:)
  !
  !--------- LIBXC vars -----------------------
  TYPE(xc_f90_pointer_t) :: xc_func
  TYPE(xc_f90_pointer_t) :: xc_info1, xc_info2
  CHARACTER(LEN=120) :: name1, name2
  INTEGER :: iexch_lxc, icorr_lxc
  INTEGER :: pol_unpol
  REAL(DP), ALLOCATABLE :: rho_lxc(:)
  REAL(DP), ALLOCATABLE :: ex_lxc(:), ec_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_lxc(:), vc_lxc(:)
  REAL(DP), ALLOCATABLE :: dex_lxc(:), dcr_lxc(:), df_lxc(:)
  !
  !
  !
  ! **************************************************************************
  ! *------------------------------------------------------------------------*
  ! * libxc funct. indexes: http://bigdft.org/Wiki/index.php?title=XC_codes  *
  ! * qe      "       "   : see comments in Modules/funct.f90                *
  ! *------------------------------------------------------------------------*
  ! *                                                                        *
  ! *  ... some examples for LDA case:                                       *
  ! *                                                                        *
  ! *                                                                        *
  ! *                         |   q-e     |   libxc  |                       *
  ! *                         |___________|__________|                       *
  ! *             slater (x)  |    1      |     1    |                       *
  ! *             pz (c)      |    1      |     9    |                       *
  ! *                         |           |          |                       *
  ! *             wigner (c)  |    5      |     2    |                       *
  ! *             vwn (c)     |    2      |     7    |                       *
  ! *             pw (c)      |    4      |    12    |                       *
  ! *             ...         |   ...     |    ...   |                       *
  ! *                                                                        *
  ! **************************************************************************
  !
  !
  !
  WRITE (*,'(/,1x,a)', ADVANCE='no') "Derivative of xc?(y/n) "
  READ(*,*) f_q
  DF_OK = .FALSE.
  IF ( TRIM(f_q) == 'y' ) DF_OK = .TRUE.
  IF ( TRIM(f_q) /= 'y' .AND. TRIM(f_q) /= 'n' ) THEN
     PRINT *, CHAR(10)//"ERROR: it is yes (y) or no (n)"//CHAR(10)
     RETURN
  ENDIF
  !
  IF ( .NOT.DF_OK ) THEN
     WRITE (*,'(/,1x,a)', ADVANCE='no') "Energy only (y/n)? "
     READ(*,*) e_q
     ENERGY_ONLY = .FALSE.
     IF ( TRIM(e_q) == 'y' ) ENERGY_ONLY = .TRUE.
     IF ( TRIM(e_q) /= 'y' .AND. TRIM(e_q) /= 'n' ) THEN
        PRINT *, CHAR(10)//"ERROR: it is yes (y) or no (n)"//CHAR(10)
        RETURN
     ENDIF
  ENDIF
  !
  WRITE (*,'(/,1x,a)', ADVANCE='no') "Polarization switch (1 unpolarized,  & 
                                                         & 2 polarized):  "
  READ(*,*) ns
  IF ( ns/=1 .AND. ns/=2 ) THEN
     PRINT *, CHAR(10)//"ERROR: you can only choose 1 or 2"//CHAR(10)
     RETURN
  ENDIF
  WRITE (*,'(/,1x,a)') "-- Functional indexes "
  WRITE (*,'(/,1x,a)', ADVANCE='no') "iexch_libxc  icorr_libxc: "
  READ(*,*) iexch_lxc, icorr_lxc
  WRITE (*,'(/,1x,a)', ADVANCE='no') "iexch_qe  icorr_qe: "
  READ(*,*) iexch_qe, icorr_qe
  IF (ns == 2 .AND. icorr_qe/=0 .AND. icorr_qe/=1 .AND. icorr_qe/=2 .AND. &
                    icorr_qe/=4 .AND. icorr_qe/=8 .AND. icorr_qe/=3 .AND. &
                    icorr_qe/=7 .AND. icorr_qe/=13) THEN
     PRINT *, CHAR(10)//" ERROR: icorr_qe not available at these conditions"//CHAR(10)
     RETURN
  ENDIF
  !
  !
  POLARIZED = .FALSE.
  IF (ns == 2) THEN
     POLARIZED = .TRUE.
  ENDIF
  !
  pol_unpol = XC_UNPOLARIZED
  np = 1
  !
  IF ( ns == 2 ) THEN
     pol_unpol = XC_POLARIZED
     np = 3
  ENDIF
  !
  !
  ! *******************************************************************************
  ! *                     Polarized case:                                         *
  ! *                                                                             *
  ! *  qe =>  rho_qe(:,1) = up    |      libxc     =>        rho_lxc(2n+1) = up   *
  ! *         rho_qe(:,2) = down  |            (dim=2*nnr)   rho_lxc(2n+2) = down *
  ! *                             |                                               *
  ! *         dmxc(:,1,1) = uu    |            (dim=3*nnr)   df_lxc(3n+1) = uu    *
  ! *         dmxc(:,1,2) = ud    |                          df_lxc(3n+2) = ud    *
  ! *         dmxc(:,2,1) = du    |                          df_lxc(3n+2) = du    *
  ! *         dmxc(:,2,2) = dd    |                          df_lxc(3n+3) = dd    *
  ! *                                                                             *
  ! *******************************************************************************
  !
  !
  ! ------ Allocations --------
  !
  ! ... qe
  !
  ALLOCATE( rho_qe(nnr,ns) )
  IF ( POLARIZED ) ALLOCATE( rho_tot(nnr), zeta(nnr) )
  ALLOCATE( ex_qe(nnr), ec_qe(nnr) )
  ALLOCATE( vx_qe(nnr,ns), vc_qe(nnr,ns) )
  IF ( DF_OK ) THEN
    IF ( .NOT.POLARIZED ) ALLOCATE( dmuxc(nnr,1,1) )
    IF ( POLARIZED ) ALLOCATE( dmuxc(nnr,2,2) )
  ENDIF
  !
  ! ... libxc
  !
  ALLOCATE( rho_lxc(nnr*ns) )
  ALLOCATE( ex_lxc(nnr), ec_lxc(nnr) )
  ALLOCATE( vx_lxc(nnr*ns), vc_lxc(nnr*ns) )
  IF ( DF_OK ) THEN
     IF ( .NOT.POLARIZED ) ALLOCATE( dex_lxc(nnr), dcr_lxc(nnr), df_lxc(nnr) )
     IF ( POLARIZED ) ALLOCATE( dex_lxc(nnr*3), dcr_lxc(nnr*3), df_lxc(nnr*3) )
  ENDIF
  !
  !
  !----- Initializations -----
  !
  ! ... qe
  !
  rho_qe = 0.0_DP
  IF ( POLARIZED ) THEN
     rho_tot = 0.0_DP
     zeta = 0.0_DP
  ENDIF
  !
  ! ... libcx
  !
  rho_lxc = 0.0_DP
  !
  !
  ! -------- Setting up an arbitrary input for both qe and libxc -----
  !
  ! ... qe
  !
  DO ii = 1, nnr
     !
     rho_qe(ii,1) = DBLE(ii)/DBLE(nnr+2)
     !
     IF ( POLARIZED ) THEN
        !
        rho_qe(ii,2) = (1.0_DP - rho_qe(ii,1))*0.7_DP
        rho_tot(ii) = rho_qe(ii,1) + rho_qe(ii,2)
        zeta(ii) = (rho_qe(ii,1) - rho_qe(ii,2)) / rho_tot(ii)
        !
     ENDIF
     !
  ENDDO
  !
  ! ... libxc
  !
  DO ii = 1, nnr
     !
     IF ( .NOT. POLARIZED ) THEN
        !
        rho_lxc(ii) = rho_qe(ii,1)
        !
     ELSE
        !
        rho_lxc(2*ii-1) = rho_qe(ii,1)
        rho_lxc(2*ii) = rho_qe(ii,2)
        !
     ENDIF
     !
  ENDDO
  !
  !
  !-------- Calculation of energies and potential -------------------
  !
  !------ LIBXC ------
  !   
  CALL xc_f90_func_init( xc_func, xc_info1, iexch_lxc, pol_unpol )         ! ... EXCHANGE   
  CALL xc_f90_lda_exc_vxc( xc_func, nnr, rho_lxc(1), ex_lxc(1), vx_lxc(1) )   
  !   
  IF ( DF_OK ) CALL xc_f90_lda_fxc( xc_func, nnr, rho_lxc(1), dex_lxc(1) )   
  CALL xc_f90_func_end( xc_func )   
  !    
  CALL xc_f90_func_init( xc_func, xc_info2, icorr_lxc, pol_unpol )         ! ... CORRELATION   
  CALL xc_f90_lda_exc_vxc( xc_func, nnr, rho_lxc(1), ec_lxc(1), vc_lxc(1) )   
  !   
  IF ( DF_OK ) THEN   
     CALL xc_f90_lda_fxc( xc_func, nnr, rho_lxc(1), dcr_lxc(1) )   
     df_lxc = (dex_lxc + dcr_lxc) * 2.0_DP   
  ENDIF   
  CALL xc_f90_func_end( xc_func )   
  !   
  !----- QE ----------   
  !   
  CALL select_lda_functionals( iexch_qe, icorr_qe )   ! ... EXCHANGE and CORRELATION   
  !   
  IF ( DF_OK ) THEN   
     IF ( .NOT.POLARIZED ) THEN   
        CALL dmxc_lda( nnr, rho_qe(:,1), dmuxc(:,1,1) )   
     ELSE   
        CALL dmxc_lsda( nnr, rho_qe, dmuxc )   
     ENDIF   
  ENDIF   
  !   
  IF ( .NOT. POLARIZED ) THEN   
     CALL xc_lda( nnr, rho_qe(:,1), ex_qe, ec_qe, vx_qe(:,1), vc_qe(:,1) )   
  ELSE   
     CALL xc_lsda( nnr, rho_tot, zeta, ex_qe, ec_qe, vx_qe, vc_qe )   
  ENDIF   
  !
  !------------------
  !
  CALL xc_f90_info_name( xc_info1, name1 )
  CALL xc_f90_info_name( xc_info2, name2 )
  !
  PRINT *, "=================================== "//CHAR(10)//" "
  PRINT *, "libxc functionals:"//CHAR(10)//" "
  PRINT *, "Exchange: ", TRIM(name1)
  PRINT *, "Correlation: ", TRIM(name2)
  PRINT *, " "
  !
  DO ii = 1, nnr, nnr-1   
     WRITE(*,*) ' '   
     WRITE(*,*) ' '   
     WRITE(*,909) ii, nnr   
     IF ( .NOT. POLARIZED ) THEN   
        WRITE (*, 401 ) rho_qe(ii,1)   
     ELSE   
        WRITE (*, 402 ) rho_qe(ii,1), rho_qe(ii,2)   
     ENDIF   
     PRINT *, " "   
     IF (.NOT. DF_OK) THEN   
        PRINT *, "=== Exchange and correlation energies: ==="   
        WRITE (*,102) ex_qe(ii),  ec_qe(ii)   
        WRITE (*,202) ex_lxc(ii), ec_lxc(ii)   
        PRINT *, " --- "   
        WRITE (*,302) ex_qe(ii)-ex_lxc(ii), ec_qe(ii)-ec_lxc(ii)   
        !   
        IF (.NOT. ENERGY_ONLY) THEN   
           PRINT *, " "   
           PRINT *, "=== Exchange potential ==="   
           IF ( .NOT. POLARIZED ) THEN   
              WRITE (*,101) vx_qe(ii,1)   
              WRITE (*,201) vx_lxc(ii)   
              PRINT *, " --- "   
              WRITE (*,301) vx_qe(ii,1)-vx_lxc(ii)   
           ELSEIF ( POLARIZED ) THEN   
              WRITE (*,102) vx_qe(ii,1), vx_qe(ii,2)   
              WRITE (*,202) vx_lxc(2*ii-1), vx_lxc(2*ii)   
              PRINT *, " --- "   
              WRITE (*,302) vx_qe(ii,1)-vx_lxc(2*ii-1), vx_qe(ii,2)-vx_lxc(2*ii)   
           ENDIF   
           PRINT *, " "   
           PRINT *, "=== Correlation potential ==="   
           IF ( .NOT. POLARIZED ) THEN   
              WRITE (*,101) vc_qe(ii,1)   
              WRITE (*,201) vc_lxc(ii)   
              PRINT *, " --- "   
              WRITE (*,301) vc_qe(ii,1)-vc_lxc(ii)   
           ELSEIF ( POLARIZED ) THEN   
              WRITE (*,102) vc_qe(ii,1), vc_qe(ii,2)   
              WRITE (*,202) vc_lxc(2*ii-1), vc_lxc(2*ii)   
              PRINT *, " --- "   
              WRITE (*,302) vc_qe(ii,1)-vc_lxc(2*ii-1), vc_qe(ii,2)-vc_lxc(2*ii)   
           ENDIF   
           PRINT *, "--- ---"   
        ENDIF   
     ELSE   
        PRINT *, " "   
        PRINT *, "=== First derivative of xc functional: ==="   
        IF ( .NOT. POLARIZED ) THEN   
           WRITE (*,101) dmuxc(ii,1,1)   
           WRITE (*,201) df_lxc(ii)   
           PRINT *, " --- "   
           WRITE (*,301) dmuxc(ii,1,1)-df_lxc(ii)   
        ELSE   
           WRITE (*,104) dmuxc(ii,1,1), dmuxc(ii,2,1), dmuxc(ii,2,2), dmuxc(ii,1,2)   
           WRITE (*,203) df_lxc(3*ii-2), df_lxc(3*ii-1), df_lxc(3*ii)   
           PRINT *, " --- "   
           WRITE (*,303) dmuxc(ii,1,1)-df_lxc(3*ii-2), dmuxc(ii,2,1)-df_lxc(3*ii-1), &   
                         dmuxc(ii,2,2)-df_lxc(3*ii)   
        ENDIF   
     ENDIF   
  ENDDO   
  !
  101 FORMAT('qe: ',3x,F17.14)
  102 FORMAT('qe: ',3x,F17.14,4x,F17.14)
  103 FORMAT('qe: ',3x,F17.14,4x,F17.14,4x,F17.14)
  104 FORMAT('qe: ',3x,F17.14,4x,F17.14,4x,F17.14,4x,F17.14)
  !
  201 FORMAT('libxc: ',F17.14)
  202 FORMAT('libxc: ',F17.14,4x,F17.14)
  203 FORMAT('libxc: ',F17.14,4x,F17.14,4x,F17.14)
  !
  301 FORMAT('diffs: ',F17.14)
  302 FORMAT('diffs: ',F17.14,4x,F17.14)
  303 FORMAT('diffs: ',F17.14,4x,F17.14,4x,F17.14)
  !
  401 FORMAT('rho: ',F17.14)
  402 FORMAT('rho(up,down): ',F17.14,4x,F17.14)
  !
  909 FORMAT('grid-point: ',I4,' of',I4)
  !
  DEALLOCATE( rho_qe )
  IF ( POLARIZED ) DEALLOCATE( rho_tot, zeta )
  DEALLOCATE( ex_qe, ec_qe )
  DEALLOCATE( vx_qe, vc_qe )
  IF (DF_OK) DEALLOCATE( dmuxc )
  !
  DEALLOCATE( rho_lxc )
  DEALLOCATE( ex_lxc, ec_lxc )
  DEALLOCATE( vx_lxc, vc_lxc )
  IF (DF_OK) DEALLOCATE( dex_lxc, dcr_lxc, df_lxc )
  !
  PRINT *, " "
  !
  RETURN
  !
#else
  !
  PRINT *, "ERROR: library libxc not included."
  RETURN
  !
#endif
  !
END PROGRAM benchmark_libxc
