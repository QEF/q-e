!
! Copyright (C) 2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!============================================================================
!============================================================================
PROGRAM xclib_test
  !==========================================================================
  !! Testing program for xc\_lib library in QE. Different options:
  !
  !! * dft-info: provides infos on the input DFT (both QE and Libxc);
  !! * xc-benchmark: difference with respect to a given set of benchmark data
  !!   (on file);
  !! * gen-benchmark: generates set of benchmark data on file;
  !! * dft-comparison: shows difference between two DFTs (E and V differences).
  !
  !! Available cases:
  !
  !! * LDA;
  !! * derivative of LDA (dmxc);
  !! * GGA;
  !! * derivative of GGA (dgcxc);
  !! * metaGGA.
  !
  USE kind_l,      ONLY: DP
  USE constants_l, ONLY: pi
  USE xc_lib,      ONLY: xclib_set_dft_from_name, xclib_set_exx_fraction, &
                         xclib_get_ID, xclib_reset_dft, xc_gcx,           &
                         xclib_dft_is_libxc, xclib_init_libxc,            &
                         xclib_finalize_libxc
  USE xclib_parallel_include
#if defined(__LIBXC)
  USE xc_f03_lib_m
  USE dft_par_mod, ONLY: xc_func, xc_info
#endif
  !
  IMPLICIT NONE
  !
#if defined(__MPI)
  INTEGER    STATUS(MPI_STATUS_SIZE)
#else
#define MPI_MAX_PROCESSOR_NAME 64
#endif
  !
#if defined(__LIBXC)
  CHARACTER(LEN=120) :: lxc_kind, lxc_family
  INTEGER :: n_ext, id(6)
#endif
  !
  INTEGER :: mype, npes, comm, ntgs, root
  LOGICAL :: iope
  INTEGER :: i, ierr, ierrm
  INTEGER :: nnodes, nlen
  !
  INTEGER, PARAMETER :: stdin  = 5
  INTEGER, PARAMETER :: stdout = 6
  !
  !-------- Grid dim vars --------------------
  INTEGER, PARAMETER :: npoints = 90000
  INTEGER :: nr, nnr, nrpe, nnr_b, nnr_int, nnrb, nnrt, nnrbt
  INTEGER :: nnrit, nnrbit, nskip
  !
  !-------- Input vars -----------------------
  CHARACTER(LEN=30) :: test, family
  CHARACTER(LEN=30) :: dft1, dft2
  INTEGER :: nspin
  LOGICAL :: DF_OK
  !
  !---------- DFT infos -------------------------
  INTEGER :: iexch1, icorr1, igcx1, igcc1, imeta1, imetac1
  INTEGER :: iexch2, icorr2, igcx2, igcc2, imeta2, imetac2
  LOGICAL :: LDA, GGA, MGGA, POLARIZED, is_libxc(6)
  CHARACTER(LEN=120) :: name1, name2
  !
  !-------- Various params -------------------
  REAL(DP), PARAMETER :: null=0.0_DP, pi34=0.6203504908994_DP
  REAL(DP), PARAMETER :: thresh_lda  = 0.d0, & !1.E-6_DP, &
                         thresh_gga  = 0.d0, & !1.E-6_DP, &
                         thresh_mgga = 0.d0    !1.E-6_DP
  REAL(DP), PARAMETER :: diff_thr_e_lda  = 1.0E-6_DP,  &
                         diff_thr_v_lda  = 1.0E-6_DP,  &
                         diff_thr_e_gga  = 1.0E-12_DP, &
                         diff_thr_vgga   = 1.0E-12_DP, &
                         diff_thr_e_mgga = 1.0E-12_DP, &
                         diff_thr_vmgga  = 1.0E-12_DP, &
                         diff_thr_dmuxc  = 1.0E-6_DP,  &
                         diff_thr_dv     = 1.0E-16_DP
  REAL(DP) :: fact, exx_frctn
  !
  !---------- Indexes ---------------------------
  INTEGER :: ii, ns, np, ipol, ithr, nthr, iip, iout
  !
  !---------- XClib input vars ------------------
  REAL(DP), ALLOCATABLE :: rho(:,:), rho_tz(:,:), rho_b(:,:),rhotz_b(:,:)
  REAL(DP), ALLOCATABLE :: grho(:,:,:), grh(:,:,:), grho_b(:,:,:), &
                           grh_b(:,:,:)
  REAL(DP), ALLOCATABLE :: tau(:,:), tau_b(:,:)                  
  REAL(DP) :: grho2(2), grho_ud
  REAL(DP) :: rhoi(2), grhoi(3,2), taui(2)
  !
  !--------- dft1 vars --------------------------
  REAL(DP), ALLOCATABLE :: ex1(:), ec1(:)
  REAL(DP), ALLOCATABLE :: exg1(:), ecg1(:)
  REAL(DP), ALLOCATABLE :: vx1(:,:), vc1(:,:)
  REAL(DP), ALLOCATABLE :: dmuxc1(:,:,:)
  REAL(DP), ALLOCATABLE :: dvxcrr1(:,:,:), dvxcsr1(:,:,:), &
                           dvxcss1(:,:,:)
  REAL(DP), ALLOCATABLE :: v1x1(:,:), v2x1(:,:), v3x1(:,:)
  REAL(DP), ALLOCATABLE :: v1c1(:,:), v2c_ud1(:)
  REAL(DP), ALLOCATABLE :: v2c1(:,:), v3c1(:,:)
  !
  !--------- dft2 vars ---------------------------
  REAL(DP), ALLOCATABLE :: ex2(:), ec2(:)
  REAL(DP), ALLOCATABLE :: exg2(:), ecg2(:)
  REAL(DP), ALLOCATABLE :: vx2(:,:), vc2(:,:)
  REAL(DP), ALLOCATABLE :: dmuxc2(:,:,:)
  REAL(DP), ALLOCATABLE :: dvxcrr2(:,:,:), dvxcsr2(:,:,:), &
                           dvxcss2(:,:,:)
  REAL(DP), ALLOCATABLE :: v1x2(:,:), v2x2(:,:), v3x2(:,:)
  REAL(DP), ALLOCATABLE :: v1c2(:,:), v2c2(:,:), v2c_ud2(:)
  REAL(DP), ALLOCATABLE :: v2cm1(:,:,:), v2cm2(:,:,:), v3c2(:,:)
  !
  !----------diff vars ---------------------------
  LOGICAL :: ex_is_out, ec_is_out, vx_is_out, vc_is_out, &
             something_out, dmuxc_is_out
  LOGICAL :: v1x_is_out, v2x_is_out, v1c_is_out, v2c_is_out, &
             v3x_is_out, v3c_is_out, &
             dvxcrr_is_out, dvxcsr_is_out, dvxcss_is_out, dvgga_is_out
  ! ... LDA aver
  REAL(DP) :: ex_aver_b(2),   ec_aver_b(2),   &
              vx_aver_b(1,2), vc_aver_b(1,2), &
              dv_aver_b(3)
  ! ... GGA/MGGA aver
  REAL(DP) :: v1x_aver_b(1,2), v1c_aver_b(1,2), &
              v2x_aver_b(1,2), v2c_aver_b(1,3), &
              v2c_ud1_aver(2), v2c_ud1_min(2), v2c_ud1_max(2), &
              dvrr_aver_b(1,3), dvsr_aver_b(1,3), &
              dvss_aver_b(1,3)
  ! ... MGGA aver
  REAL(DP) :: v3x_aver_b(1,2), v3c_aver_b(1,2)
  !
  REAL(DP) :: aver_sndu, aver_recu
  REAL(DP) :: vaver(2), vmax(2), vmin(2)
  !
  CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME), ALLOCATABLE :: proc_name(:)
  CHARACTER(LEN=MPI_MAX_PROCESSOR_NAME), ALLOCATABLE :: node_name(:)
  INTEGER, ALLOCATABLE :: proc2node(:)
#if defined(_OPENMP)
  INTEGER, EXTERNAL :: omp_get_max_threads
#endif
  !
#if defined(_OPENMP)
  INTEGER :: PROVIDED
#endif
  !
  NAMELIST/input_namelist/ test, nspin, family, DF_OK, dft1, dft2
  !
  NAMELIST/lda_benchmark_data /ex_aver_b, ec_aver_b, vx_aver_b, vc_aver_b,    &
                               ex2, ec2, vx2, vc2
  NAMELIST/dlda_benchmark_data/dv_aver_b, dmuxc2
  NAMELIST/gga_benchmark_data /ex_aver_b, ec_aver_b, v1x_aver_b, v1c_aver_b,  &
                               v2x_aver_b, v2c_aver_b,v2c_ud1_aver, ex2, ec2, &
                               v1x2, v1c2, v2x2, v2c2, v2c_ud2
  NAMELIST/dgga_benchmark_data/dvrr_aver_b, dvsr_aver_b, dvss_aver_b,dvxcrr2, &
                               dvxcsr2, dvxcss2
  NAMELIST/mgga_benchmark_data/ex_aver_b, ec_aver_b, ex_aver_b, ec_aver_b,    &
                               v1x_aver_b, v1c_aver_b, v2x_aver_b, v2c_aver_b,&
                               v3x_aver_b, v3c_aver_b, ex2, ec2, v1x2, v1c2,  &
                               v2x2, v2c2, v3x2, v3c2
#if defined(__MPI)
  !
#if defined(_OPENMP)
  CALL MPI_Init_thread(MPI_THREAD_FUNNELED, PROVIDED, ierr)
#else
  CALL MPI_Init(ierr)
#endif
  CALL mpi_comm_rank(MPI_COMM_WORLD,mype,ierr)
  CALL mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
  comm = MPI_COMM_WORLD
  ntgs = 1
  root = 0
  IF(mype==root) THEN
     iope = .TRUE.
  ELSE
     iope = .FALSE.
  ENDIF
#else
  mype = 0
  npes = 1
  comm = 0
  ntgs = 1
  root = 0
  iope = .TRUE.
#endif
  !
  ! ... init
  test = 'none'
  dft1 = 'none'
  dft2 = 'none'
  DF_OK = .FALSE.
  nspin = 1
  !
  !==========================================================================
  ! GET INPUT FROM FILE
  !==========================================================================
  !
  IF (mype==root) THEN
    READ( stdin, input_namelist )
    IF ( test(1:4)=='gen-' ) THEN
      test = 'exe-benchmark'
      WRITE( stdout, input_namelist )
      test = 'gen-benchmark'
    ENDIF  
  ENDIF
  !
  !==========================================================================
  ! MPI MANAGEMENT
  !==========================================================================
  !
#if defined(__MPI)
  CALL MPI_BCAST( test,   30, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST( family, 30, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST( dft1,   30, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST( dft2,   30, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST( nspin,   1, MPI_INT,       0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST( DF_OK,   1, MPI_LOGICAL,   0, MPI_COMM_WORLD, ierr )
#endif
  !
  ALLOCATE( proc_name( npes ) )
  ALLOCATE( node_name( npes ) )
  ALLOCATE( proc2node( npes ) )
  nnodes = 0
  !
  DO i = 1, npes
#if defined(__MPI)
    IF ( mype == i-1 ) THEN
       CALL MPI_Get_processor_name( proc_name(i), nlen, ierr )
    ENDIF
    CALL MPI_BCAST( nlen, 1, MPI_INT, i-1, MPI_COMM_WORLD, ierr )
    CALL MPI_BCAST( proc_name(i), MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER,&
                    i-1, MPI_COMM_WORLD, ierr )
#else
     proc_name(i) = 'localhost'
#endif
!      IF( mype == root ) THEN
!        write(6,310)  i, proc_name(i)
!      ENDIF
310 FORMAT(' pe = ',I5,' name = ', A20) 
    DO ii = 1, nnodes
       IF ( proc_name(i) == node_name(ii) ) THEN
          EXIT
       ENDIF
    ENDDO
    !
    IF ( ii > nnodes ) THEN
       nnodes = nnodes + 1
       node_name(nnodes) = proc_name(i)
    ENDIF
    proc2node(i) = ii
  ENDDO
  !
  IF ( mype == root .AND. test/='gen-benchmark' ) THEN
    WRITE(stdout,*) ; WRITE(stdout,*) ' --- XC_LIB TESTING PROGRAM --- '
    WRITE(stdout,*)
    WRITE(stdout,*) 'Node list:'
    WRITE(stdout,*) 
    DO ii = 1, nnodes
       WRITE(stdout,310)  ii, node_name(ii)
    ENDDO
    WRITE(stdout,*) ; WRITE(stdout,*) 'Test type: ', test
    WRITE(stdout,*)
  ENDIF
  !
#if defined(__MPI)
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
#endif
  !
  !==========================================================================
  ! PRINT DFT INFOS
  !==========================================================================
  !
  IF ( TRIM(test)=='dft-info' .AND. mype==root ) THEN
    !
    CALL xclib_set_dft_from_name( dft1 )
    !
    iexch1 = xclib_get_ID('LDA','EXCH')
    is_libxc(1) = xclib_dft_is_libxc('LDA','EXCH')
    icorr1 = xclib_get_ID('LDA','CORR')
    is_libxc(2) = xclib_dft_is_libxc('LDA','CORR')
    IF (iexch1+icorr1/=0)  LDA = .TRUE.
    igcx1 = xclib_get_ID('GGA','EXCH')
    is_libxc(3) = xclib_dft_is_libxc('GGA','EXCH')
    igcc1 = xclib_get_ID('GGA','CORR')
    is_libxc(4) = xclib_dft_is_libxc('GGA','CORR')
    IF (igcx1+igcc1/=0)    GGA = .TRUE.
    imeta1  = xclib_get_ID('MGGA','EXCH')
    is_libxc(5) = xclib_dft_is_libxc('MGGA','EXCH')
    imetac1 = xclib_get_ID('MGGA','CORR')
    is_libxc(6) = xclib_dft_is_libxc('MGGA','CORR')
    IF (imeta1+imetac1/=0) MGGA = .TRUE.
    !
    WRITE(stdout,*) " "
    WRITE(stdout,*) "=================================== "//CHAR(10)//" "
    WRITE(stdout,*) "XC functional IDs:"
    WRITE(stdout,*) CHAR(10)//"LDA IDs"
    WRITE(stdout,121) iexch1, is_libxc(1), icorr1, is_libxc(2)
    WRITE(stdout,*) CHAR(10)//"GGA IDs"
    WRITE(stdout,121) igcx1, is_libxc(3), igcc1, is_libxc(4)
    WRITE(stdout,*) CHAR(10)//"MGGA IDs"
    WRITE(stdout,121) imeta1, is_libxc(5), imetac1, is_libxc(6)
    !
    IF (ANY(.NOT.is_libxc(:))) THEN
      WRITE(stdout,*) CHAR(10)//"References for QE functionals are temporarily&
                                & listed in Modules/funct.f90"
    ENDIF
    !
#if defined(__LIBXC)
    !
    IF (xclib_dft_is_libxc('ANY')) CALL xclib_init_libxc( 1 )
    !
    WRITE(stdout,*) CHAR(10)//"LIBXC functional infos:"
    !
    id(1) = iexch1 ; id(2) = icorr1
    id(3) = igcx1  ; id(4) = igcc1
    id(5) = imeta1 ; id(6) = imetac1
    !
    DO i = 1, 6
      IF (is_libxc(i)) THEN
        WRITE(stdout,*) CHAR(10)//"Functional with ID:", id(i)
        !
        SELECT CASE( xc_f03_func_info_get_kind(xc_info(i)) )
        CASE( XC_EXCHANGE )
          WRITE(lxc_kind, '(a)') 'Exchange functional'
        CASE( XC_CORRELATION )
          WRITE(lxc_kind, '(a)') 'Correlation functional'
        CASE( XC_EXCHANGE_CORRELATION )
          WRITE(lxc_kind, '(a)') 'Exchange+Correlation functional'
        CASE( XC_KINETIC )
          WRITE(lxc_kind, '(a)') 'Kinetic energy functional - not implemented&
                                 &in QE.'
        CASE DEFAULT
          WRITE(lxc_kind, '(a)') 'Unknown kind'
        END SELECT
        !
        SELECT CASE( xc_f03_func_info_get_family(xc_info(i)) )
        CASE( XC_FAMILY_LDA )
          WRITE(lxc_family,'(a)') "LDA"
        CASE( XC_FAMILY_GGA )
          WRITE(lxc_family,'(a)') "GGA"
        CASE( XC_FAMILY_HYB_GGA )
          WRITE(lxc_family,'(a)') "Hybrid GGA"
        CASE( XC_FAMILY_MGGA )
          WRITE(lxc_family,'(a)') "MGGA"
        CASE( XC_FAMILY_HYB_MGGA )
          WRITE(lxc_family,'(a)') "Hybrid MGGA"
        CASE DEFAULT
          WRITE(lxc_family,'(a)') "unknown"
        END SELECT
        !
        WRITE(*,'("The functional ''", a, "'' is a ", a, ", it belongs to &
               &the ''", a, "'' family and is defined in the reference(s): &
               &")') TRIM(xc_f03_func_info_get_name(xc_info(i))), TRIM(lxc_kind)&
               ,TRIM(lxc_family)
        ii = 0
        DO WHILE( ii >= 0 )
         WRITE(*,'(a,i1,2a)') '[',ii+1,'] ',TRIM(xc_f03_func_reference_get_ref( &
                                   xc_f03_func_info_get_references(xc_info(i), ii)))
        ENDDO
        !
        WRITE(stdout,*)
        n_ext = xc_f03_func_info_get_n_ext_params( xc_info(i) )
        WRITE(stdout,*) 'Number of external parameters: ', n_ext
        !
        IF ( n_ext/=0 ) THEN
          DO ii = 0, n_ext-1
            WRITE(stdout,*) &
              TRIM(xc_f03_func_info_get_ext_params_description(xc_info(i), ii))
            WRITE(stdout,*) 'Default value: ', &
                   xc_f03_func_info_get_ext_params_default_value(xc_info(i), ii)
          ENDDO
        ENDIF
        !
      ENDIF
    ENDDO
    !
    IF (xclib_dft_is_libxc('ANY')) CALL xclib_finalize_libxc()
    !
#endif
    !
    121 FORMAT('Exch: ',I3,' is libxc: ',L1,';  Corr: ',I3,' is libxc: ',L1 )
    !
    GOTO 10
    !
  ENDIF !dft-info
  !
  !
  ! ... point distribution over CPUs
  !
  nr = 0
  IF ( MOD(npoints,npes)<=(mype+1) ) nr = 1
  nnr = npoints/npes + ABS(nr-1)
  IF (nr==1) nrpe = mype*nnr
  IF (nr==0) nrpe = npoints-(npoints/npes)*(npes-mype-1)
  !
  ns = nspin
  !
  ! ... initialization of averages
  !
  ex_aver_b  = 0._DP ; ec_aver_b  = 0._DP
  vx_aver_b  = 0._DP ; vc_aver_b  = 0._DP
  v1x_aver_b = 0._DP ; v2x_aver_b = 0._DP
  v3x_aver_b = 0._DP ; v1c_aver_b = 0._DP
  v2c_aver_b = 0._DP ; v3c_aver_b = 0._DP
  dv_aver_b  = 0._DP
  dvrr_aver_b = 0._DP
  dvsr_aver_b = 0._DP
  dvss_aver_b = 0._DP
  !
  ! ... initialize first DFT
  !
  CALL xclib_set_dft_from_name( dft1 )
  !
  LDA = .FALSE.
  GGA = .FALSE.
  MGGA= .FALSE.
  !
  iexch1 = xclib_get_ID('LDA','EXCH')
  icorr1 = xclib_get_ID('LDA','CORR')
  IF (iexch1+icorr1/=0)  LDA = .TRUE.
  igcx1 = xclib_get_ID('GGA','EXCH')
  igcc1 = xclib_get_ID('GGA','CORR')
  IF (igcx1+igcc1/=0)    GGA = .TRUE.
  imeta1  = xclib_get_ID('MGGA','EXCH')
  imetac1 = xclib_get_ID('MGGA','CORR')
  IF (imeta1+imetac1/=0) MGGA = .TRUE.
  !
  IF ( MGGA .AND. DF_OK ) THEN
    WRITE(stdout,*) " "
    WRITE(stdout,*) "ERROR: derivative not available with MGGA."
    GO TO 10
  ENDIF
  !
  IF (ns == 2 .AND. icorr1/=0 .AND. icorr1/=1 .AND. icorr1/=2 .AND. &
                    icorr1/=4 .AND. icorr1/=8 .AND. icorr1/=3 .AND. &
                    icorr1/=7 .AND. icorr1/=13) THEN
     WRITE(stdout,*) CHAR(10)//" ERROR: icorr1 not available at these &
                               &conditions"//CHAR(10)
     GO TO 10
  ENDIF
  !
  !
  ! ... index stuff
  !
  POLARIZED = .FALSE.
  IF (ns == 2) THEN
     POLARIZED = .TRUE.
  ENDIF
  !
  IF (.NOT. POLARIZED) THEN
    nnr_b = 1
  ELSE
    IF (LDA ) nnr_b = 2
    IF (GGA ) nnr_b = 4
    IF (MGGA) nnr_b = 5
  ENDIF  
  !
  IF (LDA ) nthr = 1
  IF (GGA ) nthr = 2
  IF (MGGA) nthr = 3
  !
  nnrt = nnr+nthr
  nnrb = nnr
  IF (test(5:13)=='benchmark') nnrb = nnr_b
  nnrbt = nnrb+nthr
  ! 
  np = 1
  IF (ns==2) np = 3
  !
  IF (xclib_dft_is_libxc('ANY')) CALL xclib_init_libxc( ns )
  !
  !==========================================================================
  ! ALLOCATIONS OF XC I/O ARRAYS
  !==========================================================================
  ! ... input
  !
  ALLOCATE( rho(nnrt,ns), rho_tz(nnrt,ns) )
  IF ( GGA .OR. MGGA ) ALLOCATE( grho(3,nnrt,ns) )
  IF ( MGGA ) ALLOCATE( tau(nnrt,ns) )
  !
  IF ( test(5:13)=='benchmark' ) THEN
    ALLOCATE( rhotz_b(nnrbt,ns), rho_b(nnr_b+nthr,ns) )
    IF ( GGA .OR. MGGA ) ALLOCATE( grho_b(3,nnr_b+nthr,ns) )
    IF ( MGGA ) ALLOCATE( tau_b(nnr_b+nthr,ns) )
  ENDIF
  !
  ! ... dft1 output arrays
  !
  IF (.NOT. DF_OK) ALLOCATE( ex1(nnrt), ec1(nnrt) )
  !
  IF ( LDA .OR. GGA ) THEN
     IF ( LDA ) THEN
       IF ( .NOT. DF_OK ) ALLOCATE( vx1(nnrt,ns), vc1(nnrt,ns) )
       IF ( DF_OK ) ALLOCATE( dmuxc1(nnrt,ns,ns) )
     ENDIF
     !
     IF ( GGA ) THEN
       IF ( .NOT. DF_OK ) THEN
         ALLOCATE( exg1(nnrt), ecg1(nnrt) )
         ALLOCATE( v1x1(nnrt,ns), v2x1(nnrt,ns) )
         ALLOCATE( v1c1(nnrt,ns), v2c1(nnrt,ns), v2c_ud1(nnrt) )
       ELSE
         ALLOCATE( grh(nnrt,3,ns) )
         IF ( test(5:13)=='benchmark' ) ALLOCATE( grh_b(nnrt,3,ns) )
         ALLOCATE( dvxcrr1(nnrt,ns,ns), dvxcsr1(nnrt,ns,ns), &
                   dvxcss1(nnrt,ns,ns) )
       ENDIF
     ENDIF
  ELSEIF ( MGGA ) THEN
     ALLOCATE( v1x1(nnrt,ns), v2x1(nnrt,ns), v3x1(nnrt,ns) )
     ALLOCATE( v1c1(nnrt,ns), v2c1(nnrt,ns), v3c1(nnrt,ns) )
     ALLOCATE( v2cm1(np,nnrt,ns) )
  ENDIF
  !
  ! ... dft2 output / benchmark data arrays
  !
  IF (.NOT. DF_OK) ALLOCATE( ex2(nnrbt), ec2(nnrbt) )
  !
  IF ( LDA .OR. GGA ) THEN
     IF ( LDA ) THEN
       IF ( .NOT. DF_OK ) ALLOCATE( vx2(nnrbt,ns), vc2(nnrbt,ns) )
       IF ( DF_OK ) ALLOCATE( dmuxc2(nnrbt,ns,ns) )
     ENDIF
     !
     IF ( GGA ) THEN
       IF ( .NOT. DF_OK ) THEN
         ALLOCATE( exg2(nnrbt), ecg2(nnrbt) )
         ALLOCATE( v1x2(nnrbt,ns), v2x2(nnrbt,ns) )
         ALLOCATE( v1c2(nnrbt,ns), v2c2(nnrbt,ns), v2c_ud2(nnrbt) )
       ELSE
         ALLOCATE( dvxcrr2(nnrbt,ns,ns), dvxcsr2(nnrbt,ns,ns), &
                   dvxcss2(nnrbt,ns,ns) )
       ENDIF
     ENDIF
  ELSEIF ( MGGA ) THEN
     ALLOCATE( v1x2(nnrbt,ns), v2x2(nnrbt,ns), v3x2(nnrbt,ns) )
     ALLOCATE( v1c2(nnrbt,ns), v2c2(nnrbt,ns), v3c2(nnrbt,ns) )
     ALLOCATE( v2cm2(np,nnrbt,ns) )
  ENDIF
  !
  !==========================================================================
  ! SET (PROVISIONAL) INPUT GRID FOR BENCHMARK TEST:
  ! =========================================================================
  ! LDA
  ! rho unpol (1p): 0.6
  ! 
  ! rho pol (2p):   0.6  0.1
  !                 0.1  0.6
  ! --------------------------
  ! GGA
  ! rho unpol (1p): 0.6       grho: 0.1 0.2 0.3
  ! 
  ! rho pol (4p):   0.6  0.1  grho: 0.1 0.2 0.3  0.4 0.3 0.2
  !                 0.1  0.6        0.1 0.2 0.3  0.4 0.3 0.2 
  !                 0.6  0.1        0.4 0.3 0.2  0.1 0.2 0.3
  !                 0.1  0.6        0.4 0.3 0.2  0.3 0.2 0.1
  ! --------------------------
  ! MGGA
  ! rho unpol (1p): 0.6       grho: 0.1 0.2 0.3               tau: 0.1
  ! 
  ! rho pol (5p):   0.6  0.1  grho: 0.1 0.2 0.3  0.4 0.3 0.2  tau: 0.1  0.2
  !                 0.1  0.6        0.1 0.2 0.3  0.4 0.3 0.2       0.1  0.2
  !                 0.6  0.1        0.4 0.3 0.2  0.1 0.2 0.3       0.1  0.2
  !                 0.1  0.6        0.4 0.3 0.2  0.3 0.2 0.1       0.1  0.2
  !                 0.1  0.6        0.4 0.3 0.2  0.3 0.2 0.1       0.2  0.1
  ! --------------------------  
  !
  IF (test(5:13)=='benchmark' .AND. mype==root) THEN
    IF (nspin == 1) THEN
      rho_b(1,1) = 0.6_DP ; rhotz_b(1,1) = 0.6_DP
      IF (family/='LDA')  grho_b(:,1,1) = (/ 0.1_DP, 0.2_DP, 0.3_DP /)
      IF (family=='MGGA') tau_b(1,1) = 0.1_DP
    ELSE
      DO i = 1, nnr_b
        IF (MOD(i,2)==0) rho_b(i,:) = (/ 0.6_DP, 0.1_DP /) 
        IF (MOD(i,2)/=0) rho_b(i,:) = (/ 0.1_DP, 0.6_DP /)
        IF (MOD(i,2)==0) rhotz_b(i,:) = (/ 0.7_DP, 0.5_DP /) 
        IF (MOD(i,2)/=0) rhotz_b(i,:) = (/ 0.7_DP, -0.5_DP /)
        IF ( family/='LDA') THEN
          IF (i<=2) THEN
            grho_b(:,i,1) = (/ 0.1_DP, 0.2_DP, 0.3_DP /)
            grho_b(:,i,2) = (/ 0.4_DP, 0.3_DP, 0.2_DP /)
          ELSE
            grho_b(:,i,1) = (/ 0.4_DP, 0.3_DP, 0.2_DP /)
            grho_b(:,i,2) = (/ 0.1_DP, 0.2_DP, 0.3_DP /)
          ENDIF
        ENDIF
      ENDDO
      !
      IF (family=='MGGA') THEN
        rho_b(5,:) = (/ 0.1_DP, 0.6_DP /)
        grho_b(:,5,1) = (/ 0.4_DP, 0.3_DP, 0.2_DP /)
        grho_b(:,5,2) = (/ 0.1_DP, 0.2_DP, 0.3_DP /)
        DO i = 1, 4
          tau_b(i,:) = (/ 0.1_DP, 0.2_DP/)
        ENDDO  
        tau_b(5,:) = (/ 0.2_DP, 0.1_DP /)
      ENDIF
    ENDIF
  ENDIF
  !
  !==========================================================================
  ! READ BENCHMARK DATA FROM FILE
  !==========================================================================
  !
  IF (test=='exe-benchmark' .AND. mype==root) THEN
    IF (.NOT. DF_OK) THEN
      IF (family=='LDA' ) READ(stdin, lda_benchmark_data)
      IF (family=='GGA' ) READ(stdin, gga_benchmark_data)
      IF (family=='MGGA') READ(stdin, mgga_benchmark_data)
    ELSE
      IF (family=='LDA')  READ(stdin, dlda_benchmark_data)
      IF (family=='GGA')  READ(stdin, dgga_benchmark_data)
    ENDIF
  ENDIF
  !
  !==========================================================================
  ! BUILD ARBITRARY INPUT FOR A LARGE GRID (for Etot and Vtot test)
  !==========================================================================
  !
  fact = (3.d0/5.d0)*(3.d0*pi*pi)**(2.0/3.0)
  !
  rho  = 0.0_DP
  IF ( GGA .OR. MGGA ) grho = 0.0_DP
  IF ( MGGA ) tau = 0.0_DP
  !
  DO ii = 1, nnr
     !
     iip = nrpe+ii
     !
     rho(ii,1) = DBLE(iip)/DBLE(npoints+2)
     IF (.NOT. POLARIZED) rho_tz(ii,1) = rho(ii,1)
     !  
     IF ( GGA .OR. MGGA ) THEN
       grho(1,ii,1) = ABS( 0.05_DP + 0.8_DP*SIN(DBLE(iip)) )
       grho(2,ii,1) = ABS( 0.05_DP + 0.7_DP*SIN(DBLE(iip)) )
       grho(3,ii,1) = ABS( 0.05_DP + 0.6_DP*SIN(DBLE(iip)) )
     ENDIF  
     !
     IF ( MGGA ) tau(ii,1) = fact*ABS(rho(ii,1)*ns)**(5._DP/3._DP)/ns
     !  
     IF ( POLARIZED ) THEN  
        !  
        rho(ii,2) = (1.0_DP - rho(ii,1))*0.7_DP  
        rho_tz(ii,1) = rho(ii,1) + rho(ii,2)  
        rho_tz(ii,2) = rho(ii,1) - rho(ii,2)  
        !  
        IF ( GGA .OR. MGGA ) THEN  
           grho(1,ii,2) = ABS( (1.0_DP - grho(1,ii,1))*0.7_DP )  
           grho(2,ii,2) = ABS( (1.0_DP - grho(2,ii,1))*0.6_DP )
           grho(3,ii,2) = ABS( (1.0_DP - grho(3,ii,1))*0.5_DP )
        ENDIF
        !  
        IF ( MGGA ) tau(ii,2) = fact*ABS(rho(ii,2)*ns)**(5._DP/3._DP)/ns
        !  
     ENDIF
     !
  ENDDO
  !
  !--- THRESHOLD POINTS ---
  !
  rho(nnr+1,1) = thresh_lda/3.0_DP
  IF (.NOT. POLARIZED) rho_tz(nnr+1,1) = rho(nnr+1,1)
  IF ( POLARIZED ) THEN
     rho(nnr+1,2) = rho(nnr,1)
     rho_tz(nnr+1,1) = rho(nnr+1,1) + rho(nnr+1,2)
     rho_tz(nnr+1,2) = rho(nnr+1,1) - rho(nnr+1,2)
  ENDIF
  !
  IF ( GGA .OR. MGGA ) THEN
    grho(:,nnr+1,1) = grho(:,nnr,1)
    !
    rho(nnr+2,:) = rho(nnr,:)
    IF (.NOT. POLARIZED) rho_tz(nnr+2,1) = rho(nnr+2,1)
    !
    grho(:,nnr+2,1) = thresh_gga/10
    IF ( POLARIZED ) THEN
      rho_tz(nnr+2,1) = rho(nnr+2,1) + rho(nnr+2,2)
      rho_tz(nnr+2,2) = rho(nnr+2,1) - rho(nnr+2,2)
      grho(:,nnr+2,2) = grho(:,nnr,2)
    ENDIF
    !
    IF ( MGGA ) THEN  
      tau(nnr+1,:) = tau(nnr,:)
      tau(nnr+2,:) = tau(nnr,:)
      rho(nnr+3,:) = rho(nnr,:)
      grho(:,nnr+3,:) = grho(:,nnr,:)
      tau(nnr+3,1) = thresh_mgga/10
      IF ( POLARIZED ) tau(nnr+3,2) = tau(nnr,2)
    ENDIF
    !
  ENDIF
  !
  IF (mype==root .AND. test(5:13)=='benchmark') THEN
    rho_b(nnrb+1:nnrbt,:) = rho(nnr+1:nnr+nthr,:)
    rhotz_b(nnrb+1:nnrbt,:) = rho_tz(nnr+1:nnr+nthr,:)
    IF (family/='LDA') grho_b(:,nnrb+1:nnrbt,:) = grho(:,nnr+1:nnr+nthr,:)
    IF (family=='MGGA') tau_b(nnrb+1:nnrbt,:) = tau(nnr+1:nnr+nthr,:)
    !
    rho(1:nnrbt,1:ns) = rho_b(1:nnrbt,1:ns)
    rho_tz(1:nnrbt,1:ns) = rhotz_b(1:nnrbt,1:ns)
    IF (family/='LDA')  grho(:,1:nnrbt,1:ns) = grho_b(:,1:nnrbt,1:ns)
    IF (family=='MGGA') tau(1:nnrbt,1:ns) = tau_b(1:nnrbt,1:ns)
  ENDIF
  !
  !==========================================================================
  ! CALCULATION OF ENERGIES AND POTENTIAL ARRAYS
  !==========================================================================
  !
  ! ... xclib calls for DFT1
  !
  IF ( LDA ) THEN
     !
     IF (.NOT. DF_OK ) CALL xc( nnrt, ns, ns, rho_tz, ex1, ec1, vx1, vc1 )
     IF ( DF_OK ) CALL dmxc( nnrt, ns, rho, dmuxc1 )
     !
  ENDIF   
  !
  IF ( GGA ) THEN
    !
    IF ( .NOT. DF_OK ) THEN
      !
      IF ( .NOT. LDA ) THEN
        ex1 = 0.d0  ;  ec1 = 0.d0
        vx1 = 0.d0  ;  vc1 = 0.d0
      ENDIF
      !
      CALL xc_gcx( nnrt, ns, rho, grho, exg1, ecg1, v1x1, v2x1, v1c1, &
                   v2c1, v2c_ud1 )
      !
      ex1 = ex1*rho_tz(:,1) + exg1
      ec1 = ec1*rho_tz(:,1) + ecg1
      v1x1 = v1x1 + vx1
      v1c1 = v1c1 + vc1
      !
    ELSE
      !
      DO ii = 1, nnrt
        grh(ii,1:3,1) = grho(1:3,ii,1)
        IF (ns==2) grh(ii,1:3,2) = grho(1:3,ii,2)
      ENDDO
      !
      CALL dgcxc( nnrt, ns, rho, grh, dvxcrr1, dvxcsr1, dvxcss1 )
      !
      dvxcrr1 = dvxcrr1 + dmuxc1
      !
    ENDIF
     !
  ENDIF
  !
  IF ( MGGA ) THEN
     CALL xc_metagcx( nnrt, ns, np, rho, grho, tau, ex1, ec1, v1x1, &
                      v2x1, v3x1, v1c1, v2cm1, v3c1 )
     v2c1 = v2cm1(1,:,:)
  ENDIF
  !
  ! ... xclib calls for DFT2
  !
  IF (test == 'dft-comparison'.OR. test(1:4)=='gen-') THEN
    !
    IF (test == 'dft-comparison') THEN
      CALL xclib_reset_dft()
      CALL xclib_set_dft_from_name( dft2 )
      IF (xclib_dft_is_libxc('ANY')) CALL xclib_init_libxc( ns )
    ENDIF
    IF (test(1:4)=='gen-') dft2 = dft1
    !
    !
    LDA = .FALSE.
    GGA = .FALSE.
    MGGA= .FALSE.
    iexch2 = xclib_get_ID('LDA','EXCH')
    icorr2 = xclib_get_ID('LDA','CORR')
    IF (iexch2+icorr2/=0)  LDA = .TRUE.
    igcx2 = xclib_get_ID('GGA','EXCH')
    igcc2 = xclib_get_ID('GGA','CORR')
    IF (igcx2+igcc2/=0)    GGA = .TRUE.
    imeta2  = xclib_get_ID('MGGA','EXCH')
    imetac2 = xclib_get_ID('MGGA','CORR')
    IF (imeta2+imetac2/=0) MGGA = .TRUE.
    !
    IF ( LDA ) THEN
      IF ( .NOT. DF_OK ) CALL xc( nnrbt, ns, ns, rho_tz(1:nnrbt,:), &
                                  ex2, ec2, vx2, vc2 )
      IF ( DF_OK ) CALL dmxc( nnrbt, ns, rho(1:nnrbt,:), dmuxc2 )
    ENDIF
    !
    !
    IF ( GGA ) THEN
      IF ( .NOT. DF_OK ) THEN
        !
        IF ( .NOT. LDA ) THEN
          ex2 = 0.d0  ;  ec2 = 0.d0
          vx2 = 0.d0  ;  vc2 = 0.d0
        ENDIF
        !
        CALL xc_gcx( nnrbt, ns, rho(1:nnrbt,:), grho(:,1:nnrbt,:), exg2, &
                     ecg2, v1x2, v2x2, v1c2, v2c2, v2c_ud2 )
        !
        ex2 = ex2*rho_tz(1:nnrbt,1) + exg2
        ec2 = ec2*rho_tz(1:nnrbt,1) + ecg2
        v1x2 = v1x2 + vx2
        v1c2 = v1c2 + vc2
        !
      ELSE
        !
        DO ii = 1, nnrbt
          grh(ii,1:3,1) = grho(1:3,ii,1)
          IF (ns==2) grh(ii,1:3,2) = grho(1:3,ii,2)
        ENDDO
        !
        CALL dgcxc( nnrbt, ns, rho(1:nnrbt,:), grh(1:nnrbt,:,:), dvxcrr2, &
                    dvxcsr2, dvxcss2 )
        !
        IF ( LDA ) dvxcrr2 = dvxcrr2 + dmuxc2
        !
      ENDIF
    ENDIF
    !
    IF ( MGGA ) THEN
      CALL xc_metagcx( nnrbt, ns, np, rho(1:nnrbt,:), grho(:,1:nnrbt,:), &
                       tau(1:nnrbt,:), ex2, ec2, v1x2, v2x2, v3x2, v1c2, &
                       v2cm2, v3c2 )
      v2c2 = v2cm2(1,:,:)
    ENDIF    
    !
  ENDIF
  !
  !
  !==========================================================================
  ! COMPUTE AND PRINT TEST RESULTS
  !==========================================================================
  !
  IF (mype==root .AND. test(1:4)/='gen-') THEN
    WRITE(stdout,*) ' '
    WRITE(stdout,911) npoints
  ENDIF  
  !
  !... calculate statistics of E over a large number of points (npoints)
  !
  IF ( .NOT. DF_OK ) THEN
    CALL evxc_stats( 'Ex', ex1, ex2, ex_aver_b )
    CALL evxc_stats( 'Ec', ec1, ec2, ec_aver_b )
  ENDIF
  !
  !
  IF ( LDA .AND. .NOT. GGA ) THEN
     !
     ! ... calculate statistics of V over a large number of points (npoints)
     !
     IF ( .NOT. DF_OK ) THEN
       CALL evxc_stats( 'Vx', vx1, vx2, vx_aver_b(1,:) )
       CALL evxc_stats( 'Vc', vc1, vc2, vc_aver_b(1,:) )
     ELSE
       CALL derivxc_stats( 'dmuxc', dmuxc1, dmuxc2, dv_aver_b )
     ENDIF
     !
     !
     IF (test(5:13)=='benchmark') THEN
       ex_is_out = .TRUE.     ; ec_is_out = .TRUE.
       vx_is_out = .TRUE.     ; vc_is_out = .TRUE. 
       something_out = .TRUE. ; dmuxc_is_out = .TRUE.
     ENDIF  
     !
     ! 
     iout = 0
     !
     DO ii = 1, nnrb
        !
        IF ( mype/=root ) CYCLE
        !
        IF ( test == 'dft-comparison' ) THEN
          !
          IF ( .NOT. DF_OK ) THEN
            ex_is_out = is_it_out( diff_thr_e_lda, 1, ex1(ii:ii), ex2(ii:ii) )
            ec_is_out = is_it_out( diff_thr_e_lda, 1, ec1(ii:ii), ec2(ii:ii) )
            vx_is_out = is_it_out( diff_thr_v_lda,ns, vx1(ii,:), vx2(ii,:) )
            vc_is_out = is_it_out( diff_thr_v_lda,ns, vc1(ii,:), vc2(ii,:) )
            something_out = ANY((/ex_is_out, ec_is_out, vx_is_out, vc_is_out/))
          ELSE
            dmuxc_is_out = is_dit_out( diff_thr_dmuxc, dmuxc1(ii,:,:), &
                                       dmuxc2(ii,:,:) )
          ENDIF
          !
        ENDIF
        !        
        IF ( something_out .OR. dmuxc_is_out ) THEN
          !
          iout = iout + 1
          !
          IF (iout<=10) THEN
            IF (test(1:4)/='gen-') WRITE(stdout,*) " "
            !
            IF ( test=='exe-benchmark' ) THEN
              rhoi(1:ns)=rho_b(ii,1:ns)
              WRITE(stdout,909) nrpe+ii, nnr_b
            ELSEIF ( test=='dft-comparison' ) THEN
              rhoi(1:ns) = rho(ii,1:ns)
              WRITE(stdout,909) nrpe+ii, npoints
            ENDIF
            !
            IF (test/='gen-benchmark') THEN
              IF ( .NOT. POLARIZED ) WRITE(stdout, 401 ) rhoi(1)
              IF (       POLARIZED ) WRITE(stdout, 402 ) rhoi(1), rhoi(2)
              WRITE(stdout,*) " "
            ENDIF
            !
            IF (.NOT. DF_OK) THEN
              IF ( ex_is_out ) CALL print_diff( 'Ex', ex1(ii:ii), ex2(ii:ii) )
              IF ( ec_is_out ) CALL print_diff( 'Ec', ec1(ii:ii), ec2(ii:ii) )
              IF ( vx_is_out ) CALL print_diff( 'Vx', vx1(ii,:), vx2(ii,:) )
              IF ( vc_is_out ) CALL print_diff( 'Vc', vc1(ii,:), vc2(ii,:) )
            ELSE  
              CALL print_diff2( 'dmuxc', dmuxc1(ii,:,:), dmuxc2(ii,:,:) )
            ENDIF !df_ok
            !
          ENDIF
          !
        ENDIF
        !
     ENDDO
     !
     ! ... THRESHOLD TEST 
     !
     IF ( mype == root ) THEN
       IF ( test/='gen-benchmark' ) THEN
         WRITE(stdout,*) " "
         WRITE(stdout,*) "--- INPUT THRESHOLD CHECK ---"
         WRITE(stdout,*) " "
       ENDIF  
       !
       DO ithr = 1, nthr
         !
         nnrit  = nnr + ithr
         nnrbit = nnrb + ithr
         ! 
         IF ( test/='gen-benchmark' ) THEN
           WRITE(stdout,*) " "
           WRITE(stdout,910) ithr, nthr
         ENDIF
         !
         rhoi(1:ns) = rho(nnrb+ithr,1:ns)
         IF ( test(5:13)=='benchmark' ) rhoi(1:ns)=rho_b(nnrbit,1:ns)
         !
         IF (test/='gen-benchmark') THEN
           IF ( .NOT. POLARIZED ) WRITE(stdout, 401 ) rhoi(1)
           IF (       POLARIZED ) WRITE(stdout, 402 ) rhoi(1), rhoi(2)
           WRITE(stdout,*) " "
         ENDIF
         !
         IF (.NOT. DF_OK) THEN
           CALL print_diff('Ex',ex1(nnrit:nnrit),ex2(nnrbit:nnrbit))
           CALL print_diff('Ec',ec1(nnrit:nnrit),ec2(nnrbit:nnrbit))
           CALL print_diff( 'Vx', vx1(nnrit,:), vx2(nnrbit,:) )
           CALL print_diff( 'Vc', vc1(nnrit,:), vc2(nnrbit,:) )
         ELSE
           CALL print_diff2('dmuxc',dmuxc1(nnrit,:,:),dmuxc2(nnrbit,:,:))
         ENDIF !df_ok
       ENDDO
     ENDIF
     !
     !
  ELSEIF ( GGA ) THEN
     !
     ! ... calculate statistics over a large number of points (npoints)
     !
     IF ( .NOT. DF_OK ) THEN
       CALL evxc_stats( 'V1x', v1x1, v1x2, v1x_aver_b )
       CALL evxc_stats( 'V2x', v2x1, v2x2, v2x_aver_b )
       CALL evxc_stats( 'V1c', v1c1, v1c2, v1c_aver_b )
       CALL evxc_stats( 'V2c', v2c1, v2c2, v2c_aver_b )
     ELSE
       CALL derivxc_stats( 'dvxcrr', dvxcrr1, dvxcrr2, dvrr_aver_b )
       CALL derivxc_stats( 'dvxcsr', dvxcsr1, dvxcsr2, dvsr_aver_b )
       CALL derivxc_stats( 'dvxcss', dvxcss1, dvxcss2, dvss_aver_b )
     ENDIF
     !
     !
     IF (test(5:13)=='benchmark')  THEN
       ex_is_out = .TRUE.  ; ec_is_out = .TRUE.
       v1x_is_out = .TRUE. ; v2x_is_out = .TRUE.
       v1c_is_out = .TRUE. ; v2c_is_out = .TRUE.
       dvxcrr_is_out = .TRUE. ; dvxcsr_is_out = .TRUE.
       dvxcss_is_out = .TRUE.
       something_out = .TRUE. ; dvgga_is_out = .TRUE.
     ENDIF
     !
     iout = 0
     !
     DO ii = 1, nnrb
       !
       IF (mype/=root) CYCLE
       !
       IF (test=='dft-comparison') THEN
         IF ( .NOT. DF_OK ) THEN
           ex_is_out = is_it_out( diff_thr_e_gga, 1, ex1(ii:ii), ex2(ii:ii) )
           ec_is_out = is_it_out( diff_thr_e_gga, 1, ec1(ii:ii), ec2(ii:ii) )
           v1x_is_out= is_it_out( diff_thr_vgga,ns, v1x1(ii,:),v1x2(ii,:) )
           v2x_is_out= is_it_out( diff_thr_vgga,ns, v2x1(ii,:),v2x2(ii,:) )
           v1c_is_out= is_it_out( diff_thr_vgga,ns, v1c1(ii,:),v1c2(ii,:) )
           v2c_is_out= is_it_out( diff_thr_vgga,np, v2c1(ii,:),v2c2(ii,:), &
                                                 v2c_ud1(ii), v2c_ud2(ii) )
           something_out=ANY((/ex_is_out, ec_is_out, v1x_is_out, v2x_is_out, &
                               v1c_is_out, v2c_is_out/) )
         ELSE           
           dvxcrr_is_out = is_dit_out(diff_thr_dv,dvxcrr1(ii,:,:),dvxcrr2(ii,:,:))
           dvxcsr_is_out = is_dit_out(diff_thr_dv,dvxcsr1(ii,:,:),dvxcsr2(ii,:,:))
           dvxcss_is_out = is_dit_out(diff_thr_dv,dvxcss1(ii,:,:),dvxcss2(ii,:,:))
           dvgga_is_out = ANY((/dvxcrr_is_out, dvxcsr_is_out, dvxcss_is_out/))
         ENDIF
       ENDIF
       !
       !
       IF ( something_out .OR. dvgga_is_out ) THEN
         !
         iout = iout + 1 
         ! 
         IF (iout<=10) THEN 
           !
           IF ( test=='exe-benchmark' ) THEN 
             WRITE(stdout,*) " "
             rhoi(1:ns)=rho_b(ii,1:ns) ; grhoi(:,1:ns) = grho_b(:,ii,1:ns)
             WRITE(stdout,909) nrpe+ii, nnr_b
           ELSEIF ( test=='dft-comparison' ) THEN 
             WRITE(stdout,*) " "
             rhoi(1:ns) = rho(ii,1:ns) ; grhoi(:,1:ns) = grho(:,ii,1:ns) 
             WRITE(stdout,909) nrpe+ii, npoints
           ENDIF
           !
           IF ( test/='gen-benchmark' ) THEN 
             IF (.NOT. POLARIZED ) THEN 
               WRITE(stdout,401) rhoi(1) 
               grho2(1) = grhoi(1,1)**2 + grhoi(2,1)**2 + grhoi(3,1)**2 
               WRITE(stdout,501) grho2(1) 
             ELSE 
               WRITE(stdout,402) rhoi(1), rhoi(2) 
               grho2(1) = grhoi(1,1)**2 + grhoi(2,1)**2 + grhoi(3,1)**2 
               grho2(2) = grhoi(1,2)**2 + grhoi(2,2)**2 + grhoi(3,2)**2 
               grho_ud  = grhoi(1,1) * grhoi(1,2) + & 
                          grhoi(2,1) * grhoi(2,2) + & 
                          grhoi(3,1) * grhoi(3,2) 
               WRITE(stdout,503) grho2(1), grho_ud, grho2(2) 
             ENDIF 
             WRITE(stdout,*) " "
           ENDIF  
           ! 
           ! 
           IF (.NOT. DF_OK) THEN 
             ! 
             IF (ex_is_out) CALL print_diff( 'Ex', ex1(ii:ii), ex2(ii:ii) ) 
             IF (ec_is_out) CALL print_diff( 'Ec', ec1(ii:ii), ec2(ii:ii) ) 
             IF (v1x_is_out) CALL print_diff( 'V1x',v1x1(ii,:), v1x2(ii,:) ) 
             IF (v2x_is_out) CALL print_diff( 'V2x',v2x1(ii,:), v2x2(ii,:) ) 
             IF (v1c_is_out) CALL print_diff( 'V1c',v1c1(ii,:), v1c2(ii,:) ) 
             IF (v2c_is_out) CALL print_diff( 'V2c',v2c1(ii,:), v2c2(ii,:), & 
                                               v2c_ud1(ii), v2c_ud2(ii) )
           ELSE
             ! 
             IF (test/='gen-benchmark') WRITE(stdout,*) " " 
             ! 
             IF (dvxcrr_is_out) CALL print_diff2( 'dvxcrr',dvxcrr1(ii,:,:), &
                                                           dvxcrr2(ii,:,:) ) 
             IF (dvxcsr_is_out) CALL print_diff2( 'dvxcsr',dvxcsr1(ii,:,:), &
                                                           dvxcsr2(ii,:,:) ) 
             IF (dvxcss_is_out) CALL print_diff2( 'dvxcss',dvxcss1(ii,:,:), &
                                                           dvxcss2(ii,:,:) ) 
           ENDIF 
           ! 
         ENDIF !iout 
         ! 
       ENDIF 
       ! 
     ENDDO 
     ! 
     ! 
     ! ... THRESHOLD TEST 
     ! 
     IF (mype==root) THEN
       IF (test(1:4)/='gen-' ) THEN
         WRITE(stdout,*) " " 
         WRITE(stdout,*) "--- INPUT THRESHOLD CHECK ---" 
         WRITE(stdout,*) " " 
       ENDIF  
       ! 
       DO ithr = 1, nthr 
         !
         nnrbit = nnrb + ithr
         !
         rhoi(1:ns) = rho(nnrbit,1:ns)
         grhoi(:,1:ns) = grho(:,nnrbit,1:ns) 
         IF ( test=='exe-benchmark' ) THEN 
           rhoi(1:ns)=rho_b(nnrbit,1:ns)
           grhoi(:,1:ns) = grho_b(:,nnrbit,1:ns) 
         ENDIF   
         ! 
         IF (test/='gen-benchmark') THEN
           IF (.NOT. POLARIZED ) THEN 
             WRITE(stdout,*) " "
             WRITE(stdout,401) rhoi(1) 
             grho2(1) = grhoi(1,1)**2 + grhoi(2,1)**2 + grhoi(3,1)**2 
             WRITE(stdout,501) grho2(1) 
           ELSE 
             WRITE(stdout,*) " "
             WRITE(stdout,402) rhoi(1), rhoi(2) 
             grho2(1) = grhoi(1,1)**2 + grhoi(2,1)**2 + grhoi(3,1)**2 
             grho2(2) = grhoi(1,2)**2 + grhoi(2,2)**2 + grhoi(3,2)**2 
             grho_ud  = grhoi(1,1) * grhoi(1,2) + & 
                        grhoi(2,1) * grhoi(2,2) + & 
                        grhoi(3,1) * grhoi(3,2) 
             WRITE(stdout,503) grho2(1), grho_ud, grho2(2) 
           ENDIF 
           WRITE(stdout,*) " " 
         ENDIF  
         ! 
         IF (.NOT. DF_OK) THEN 
           CALL print_diff( 'Ex', ex1(nnrbit:nnrbit), ex2(nnrbit:nnrbit) ) 
           CALL print_diff( 'Ec', ec1(nnrbit:nnrbit), ec2(nnrbit:nnrbit) ) 
           CALL print_diff( 'V1x', v1x1(nnrbit,:), v1x2(nnrbit,:) ) 
           CALL print_diff( 'V2x', v2x1(nnrbit,:), v2x2(nnrbit,:) ) 
           CALL print_diff( 'V1c', v1c1(nnrbit,:), v1c2(nnrbit,:) ) 
           CALL print_diff( 'V2c', v2c1(nnrbit,:), v2c2(nnrbit,:),  & 
                                   v2c_ud1(nnrbit), v2c_ud2(nnrbit) )
         ELSE   
           CALL print_diff2('dvxcrr',dvxcrr1(nnrbit,:,:), dvxcrr2(nnrbit,:,:)) 
           CALL print_diff2('dvxcsr',dvxcsr1(nnrbit,:,:), dvxcsr2(nnrbit,:,:)) 
           CALL print_diff2('dvxcss',dvxcss1(nnrbit,:,:), dvxcss2(nnrbit,:,:)) 
         ENDIF 
         ! 
       ENDDO !df_ok 
       ! 
     ENDIF 
     ! 
     !
  ELSEIF ( MGGA ) THEN 
     !
     ! ... calculate statistics over a large number of points (npoints) 
     ! 
     CALL evxc_stats( 'V1x', v1x1, v1x2, v1x_aver_b(1,:) )
     CALL evxc_stats( 'V2x', v2x1, v2x2, v2x_aver_b(1,:) )
     CALL evxc_stats( 'V1c', v1c1, v1c2, v1c_aver_b(1,:) )
     CALL evxc_stats( 'V2c', v2c1, v2c2, v2c_aver_b(1,:) )
     CALL evxc_stats( 'V3x', v3x1, v3x2, v3x_aver_b(1,:) )
     CALL evxc_stats( 'V3c', v3c1, v3c2, v3c_aver_b(1,:) )
     ! 
     ! ... calculate values over a few benchmark points (nnr_b) 
     ! 
     IF (mype == root) THEN 
       !
       IF (test(5:13)=='benchmark') THEN
         ex_is_out = .TRUE.  ; ec_is_out = .TRUE.
         v1x_is_out = .TRUE. ; v2x_is_out = .TRUE. 
         v3x_is_out = .TRUE. ; v1c_is_out = .TRUE. 
         v2c_is_out = .TRUE. ; v3c_is_out = .TRUE. 
         something_out = .TRUE. 
       ENDIF 
       ! 
       iout = 0 
       !
       DO ii = 1, nnrb 
         !
         IF (test=='dft-comparison') THEN 
           !
           ex_is_out = is_it_out( diff_thr_e_mgga, 1, ex1(ii:ii), ex2(ii:ii) ) 
           ec_is_out = is_it_out( diff_thr_e_mgga, 1, ec1(ii:ii), ec2(ii:ii) ) 
           v1x_is_out = is_it_out( diff_thr_vmgga,ns,v1x1(ii,:),v1x2(ii,:) )
           v2x_is_out = is_it_out( diff_thr_vmgga,ns,v2x1(ii,:),v2x2(ii,:) )
           v3x_is_out = is_it_out( diff_thr_vmgga,ns,v3x1(ii,:),v3x2(ii,:) )
           v1c_is_out = is_it_out( diff_thr_vmgga,ns,v1c1(ii,:),v1c2(ii,:) )
           v2c_is_out = is_it_out( diff_thr_vmgga,ns,v2c1(ii,:),v2c2(ii,:) )
           v3c_is_out = is_it_out( diff_thr_vmgga,ns,v3c1(ii,:),v3c2(ii,:) )
           something_out = ANY((/ex_is_out,ec_is_out, v1x_is_out, v2x_is_out, &
                                 v3x_is_out, v1c_is_out, v2c_is_out, &
                                 v3c_is_out/))
         ENDIF
         !
         IF ( something_out ) THEN 
            ! 
            iout = iout + 1 
            ! 
            IF (iout<=10) THEN 
              IF (test(1:4)/='gen-') THEN
                WRITE(stdout,*) " " 
                IF ( test=='exe-benchmark' ) THEN
                  WRITE(stdout,909) ii, nnr_b
                  rhoi(1:ns) = rho_b(ii,1:ns)
                  grhoi(:,1:ns) = grho_b(:,ii,1:ns)
                  taui(1:ns) = tau_b(ii,1:ns) 
                ELSE 
                  WRITE(stdout,909) ii, npoints
                  rhoi(1:ns) = rho(ii,1:ns)
                  grhoi(:,1:ns) = grho(:,ii,1:ns) 
                  taui(1:ns) = tau(ii,1:ns)
                ENDIF
              ENDIF  
              ! 
              IF (test/='gen-benchmark') THEN
                IF (.NOT. POLARIZED ) THEN 
                  WRITE(stdout,401) rhoi(1) 
                  grho2(1) = grhoi(1,1)**2 + grhoi(2,1)**2 + grhoi(3,1)**2 
                  WRITE(stdout,501) grho2(1) 
                  WRITE(stdout,601) taui(1) 
                ELSE 
                  WRITE(stdout,402) rhoi(1), rhoi(2) 
                  grho2(1) = grhoi(1,1)**2 + grhoi(2,1)**2 + grhoi(3,1)**2 
                  grho2(2) = grhoi(1,2)**2 + grhoi(2,2)**2 + grhoi(3,2)**2 
                  grho_ud  = grhoi(1,1) * grhoi(1,2) + & 
                             grhoi(2,1) * grhoi(2,2) + & 
                             grhoi(3,1) * grhoi(3,2) 
                  WRITE(stdout,503) grho2(1), grho_ud, grho2(2) 
                  WRITE(stdout,602) taui(1), taui(2) 
                ENDIF 
              ENDIF
              ! 
              IF (ex_is_out ) CALL print_diff( 'Ex', ex1(ii:ii), ex2(ii:ii) )
              IF (ec_is_out ) CALL print_diff( 'Ec', ec1(ii:ii), ec2(ii:ii) )
              IF (v1x_is_out) CALL print_diff( 'V1x',v1x1(ii,:), v1x2(ii,:) )
              IF (v2x_is_out) CALL print_diff( 'V2x',v2x1(ii,:), v2x2(ii,:) )
              IF (v1c_is_out) CALL print_diff( 'V1c',v1c1(ii,:), v1c2(ii,:) )
              IF (v2c_is_out) CALL print_diff( 'V2c',v2c1(ii,:), v2c2(ii,:) )
              IF (v3x_is_out) CALL print_diff( 'V3x',v3x1(ii,:), v3x2(ii,:) )
              IF (v3c_is_out) CALL print_diff( 'V3c',v3c1(ii,:), v3c2(ii,:) )
              ! 
           ENDIF  
         ENDIF 
         ! 
       ENDDO 
     ENDIF   
     ! 
     ! 
     ! ... THRESHOLD TEST 
     ! 
     IF ( mype==root ) THEN 
       !
       IF (test/='gen-benchmark') THEN
         WRITE(stdout,*) " "  
         WRITE(stdout,*) "--- INPUT THRESHOLD CHECK ---"  
         WRITE(stdout,*) " "  
       ENDIF  
       !  
       DO ithr = 1, nthr  
         !
         nnrbit = nnrb + ithr
         !
         rhoi(1:ns) = rho(nnrbit,1:ns)
         grhoi(:,1:ns) = grho(:,nnr+ithr,1:ns) 
         taui(1:ns) = tau(nnrbit,1:ns) 
         IF ( test=='exe-benchmark' ) THEN 
           rhoi(1:ns) = rho_b(nnrbit,1:ns)
           grhoi(:,1:ns) = grho_b(:,nnrbit,1:ns) 
           taui(1:ns) = tau_b(nnrbit,1:ns) 
         ENDIF 
         !
         IF (test/='gen-benchmark') THEN
           WRITE(stdout,*) " " 
           WRITE(stdout,910) ithr, nthr 
           IF (.NOT. POLARIZED ) THEN    
             WRITE(stdout,401) rhoi(1) 
             grho2(1) = grhoi(1,1)**2 + grhoi(2,1)**2 + grhoi(3,1)**2    
             WRITE(stdout,501) grho2(1)  
             WRITE(stdout,601) taui(1) 
           ELSE    
             WRITE(stdout,402) rhoi(1), rhoi(2)    
             grho2(1) = grhoi(1,1)**2 + grhoi(2,1)**2 + grhoi(3,1)**2    
             grho2(2) = grhoi(1,2)**2 + grhoi(2,2)**2 + grhoi(3,2)**2    
             grho_ud  = grhoi(1,1) * grhoi(1,2) + &
                        grhoi(2,1) * grhoi(2,2) + &
                        grhoi(3,1) * grhoi(3,2)    
             WRITE(stdout,503) grho2(1), grho_ud, grho2(2) 
             WRITE(stdout,602) taui(1), taui(2) 
           ENDIF
           WRITE(stdout,*) " " 
         ENDIF
         !    
         CALL print_diff( 'Ex', ex1(nnrbit:nnrbit), ex2(nnrbit:nnrbit) )    
         CALL print_diff( 'Ec', ec1(nnrbit:nnrbit), ec2(nnrbit:nnrbit) )
         CALL print_diff( 'V1x', v1x1(nnrbit,:), v1x2(nnrbit,:) )
         CALL print_diff( 'V2x', v2x1(nnrbit,:), v2x2(nnrbit,:) )
         CALL print_diff( 'V1c', v1c1(nnrbit,:), v1c2(nnrbit,:) )
         CALL print_diff( 'V2c', v2c1(nnrbit,:), v2c2(nnrbit,:) )
         CALL print_diff( 'V3x', v3x1(nnrbit,:), v3x2(nnrbit,:) )
         CALL print_diff( 'V3c', v3c1(nnrbit,:), v3c2(nnrbit,:) )
         !  
       ENDDO
       ! 
     ENDIF  !mype
     ! 
  ENDIF
  !
  !
  IF (test=='gen-benchmark' .AND. mype==root) THEN
    IF (.NOT. DF_OK) THEN
      IF ( family=='LDA' ) WRITE(stdout,lda_benchmark_data)
      IF ( family=='GGA' ) WRITE(stdout,gga_benchmark_data)
      IF ( family=='MGGA') WRITE(stdout,mgga_benchmark_data)
    ELSE
      IF ( family=='LDA' ) WRITE(stdout,dlda_benchmark_data)
      IF ( family=='GGA' ) WRITE(stdout,dgga_benchmark_data)
    ENDIF
  ENDIF  
  ! 
  !
  401 FORMAT('rho: ',F17.14)
  402 FORMAT('rho(up,down): ',F17.14,4x,F17.14)
  !
  501 FORMAT('grho2: ',F17.14)
  502 FORMAT('grho2(uu,dd): ',F17.14,4x,F17.14)
  503 FORMAT('grho2(uu,ud,dd): ',F17.14,4x,F17.14,4x,F17.14)
  !
  601 FORMAT('tau: ',F17.14)
  602 FORMAT('tau(up,down): ',F17.14,4x,F17.14)
  !
  909 FORMAT('grid-point: ',I5,' of ',I5)
  910 FORMAT('threshold-point: ',I4,' of ',I4)
  911 FORMAT(' TOTAL VALUES OVER ',I5,' POINTS')
  !
  !
  !==========================================================================
  ! FINALIZE
  !==========================================================================
  !
  IF (xclib_dft_is_libxc('ANY')) CALL xclib_finalize_libxc()
  !
  DEALLOCATE( rho, rho_tz )
  !
  IF ( GGA .OR. MGGA ) DEALLOCATE( grho )
  IF ( MGGA ) DEALLOCATE( tau )
  !
  IF ( test(5:13)=='benchmark' ) THEN
    DEALLOCATE( rho_b, rhotz_b )
    IF ( GGA .OR. MGGA ) DEALLOCATE( grho_b )
    IF ( MGGA ) DEALLOCATE( tau_b )
  ENDIF
  !
  ! ... dft1 output arrays
  !
  IF (.NOT. DF_OK) DEALLOCATE( ex1, ec1 )
  !
  IF ( LDA .OR. GGA ) THEN
     IF ( LDA ) THEN
       IF ( .NOT. DF_OK ) DEALLOCATE( vx1, vc1 )
       IF ( DF_OK ) DEALLOCATE( dmuxc1 )
     ENDIF
     !
     IF ( GGA ) THEN
       IF ( .NOT. DF_OK ) THEN
         DEALLOCATE( exg1, ecg1 )
         DEALLOCATE( v1x1, v2x1 )
         DEALLOCATE( v1c1, v2c1, v2c_ud1 )
       ELSE
         DEALLOCATE( grh )
         IF ( test(5:13)=='benchmark' ) DEALLOCATE( grh_b )
         DEALLOCATE( dvxcrr1, dvxcsr1, dvxcss1 )
       ENDIF
     ENDIF
  ELSEIF ( MGGA ) THEN
     DEALLOCATE( v1x1, v2x1, v3x1 )
     DEALLOCATE( v1c1, v2c1, v3c1 )
     DEALLOCATE( v2cm1 )
  ENDIF
  !
  ! ... dft2 output / benchmark data arrays
  !
  IF (.NOT. DF_OK) DEALLOCATE( ex2, ec2 )
  !
  IF ( LDA .OR. GGA ) THEN
     IF ( LDA ) THEN
       IF ( .NOT. DF_OK ) DEALLOCATE( vx2, vc2 )
       IF ( DF_OK ) DEALLOCATE( dmuxc2 )
     ENDIF
     !
     IF ( GGA ) THEN
       IF ( .NOT. DF_OK ) THEN
         DEALLOCATE( exg2, ecg2 )
         DEALLOCATE( v1x2, v2x2 )
         DEALLOCATE( v1c2, v2c2, v2c_ud2 )
       ELSE
         DEALLOCATE( dvxcrr2, dvxcsr2, dvxcss2 )
       ENDIF
     ENDIF
  ELSEIF ( MGGA ) THEN
     DEALLOCATE( v1x2, v2x2, v3x2 )
     DEALLOCATE( v1c2, v2c2, v3c2 )
     DEALLOCATE( v2cm2 )
  ENDIF
  !
  DEALLOCATE( proc_name )
  DEALLOCATE( node_name )
  DEALLOCATE( proc2node )
  !
  10 CONTINUE
  !
#if defined(__MPI)
  CALL mpi_finalize( ierr )
#endif
  !
  WRITE(stdout,*) " "
  !
  STOP
  !
  !
 CONTAINS
 !
 !
 !------------------------------------------------------------------------
 SUBROUTINE diff_average( thr, x_dft1, x_dft2, aver_abs_perc, nnr_nt )
  !----------------------------------------------------------------------
  !! Calculates average difference (both absolute and percentage) between
  !! dft1 and dft2 quantities.
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: thr
  REAL(DP), INTENT(IN)  :: x_dft1(nnr), x_dft2(nnr)
  REAL(DP), INTENT(OUT) :: aver_abs_perc(2)
  !! 1: absolute difference;  2: percentage difference
  INTEGER, INTENT(OUT) :: nnr_nt
  !
  INTEGER  :: i
  !REAL(DP) :: calc_perc_diff
  REAL(DP) :: abs_diff, perc_diff
  !
  nnr_nt = 0
  aver_abs_perc = 0._DP
  !
  DO i = 1, nnr
    abs_diff  = ABS(x_dft1(i) - x_dft2(i))
    perc_diff = calc_perc_diff( thr, x_dft1(i), x_dft2(i) )
    aver_abs_perc(1) = aver_abs_perc(1) + abs_diff
    IF ( perc_diff < 0._DP ) CYCLE
    nnr_nt = nnr_nt+1
    aver_abs_perc(2) = aver_abs_perc(2) + perc_diff
  ENDDO
  !
  aver_abs_perc(1) = aver_abs_perc(1) / DBLE(npoints)
  aver_abs_perc(2) = aver_abs_perc(2) / DBLE(npoints)
  !
  RETURN
  !
 END SUBROUTINE diff_average
 !
 !
 !------------------------------------------------------------------------
 SUBROUTINE diff_max( thr, x_dft1, x_dft2, max_abs_perc )
  !-----------------------------------------------------------------------
  !! Finds the max difference (both absolute and percentage) between dft1
  !! and dft2 quantities.
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: thr
  REAL(DP), INTENT(IN)  :: x_dft1(nnr), x_dft2(nnr)
  REAL(DP), INTENT(OUT) :: max_abs_perc(2)
  !
  INTEGER :: i
  !REAL(DP) :: calc_perc_diff
  REAL(DP) :: abs_diff, perc_diff
  REAL(DP) :: abs_diff_prev, perc_diff_prev
  !
  abs_diff_prev  = 0.0_DP
  perc_diff_prev = 0.0_DP
  !
  DO i = 1, nnr
    abs_diff  = ABS(x_dft1(i) - x_dft2(i))
    perc_diff = calc_perc_diff( thr, x_dft1(i), x_dft2(i) )
    IF ( abs_diff > abs_diff_prev ) THEN
      max_abs_perc(1) = abs_diff
      abs_diff_prev = abs_diff
    ENDIF
    IF ( perc_diff > perc_diff_prev ) THEN
      max_abs_perc(2) = perc_diff
      perc_diff_prev = perc_diff
    ENDIF
  ENDDO
  !
  RETURN
  !
 END SUBROUTINE diff_max
 !
 !
 !-------------------------------------------------------------------------
 SUBROUTINE diff_min( thr, x_dft1, x_dft2, min_abs_perc )
  !------------------------------------------------------------------------
  !! Finds the max difference (both absoulte and percentage) between dft1
  !! and dft2 quantities.
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: thr
  REAL(DP), INTENT(IN)  :: x_dft1(nnr), x_dft2(nnr)
  REAL(DP), INTENT(OUT) :: min_abs_perc(2)
  !
  INTEGER  :: i
  !REAL(DP) :: calc_perc_diff
  REAL(DP) :: abs_diff, perc_diff
  REAL(DP) :: perc_diff_prev, abs_diff_prev
  !
  abs_diff_prev  = 1000._DP
  perc_diff_prev = 1000._DP
  !
  DO i = 1, nnr
    abs_diff  = ABS(x_dft1(i) - x_dft2(i))
    perc_diff = calc_perc_diff( thr, x_dft1(i), x_dft2(i) )
    IF ( abs_diff < abs_diff_prev .and. abs_diff>0._DP ) THEN
      min_abs_perc(1) = abs_diff
      abs_diff_prev = abs_diff
    ENDIF
    IF ( perc_diff < 0._DP ) CYCLE
    IF ( perc_diff < perc_diff_prev ) THEN
      min_abs_perc(2) = perc_diff
      perc_diff_prev = perc_diff
    ENDIF
  ENDDO
  !
  RETURN
  !
 END SUBROUTINE diff_min
 !
 !--------------------------------------------------------------------
 SUBROUTINE print_stat( what, vaver, vmax, vmin, averref )
  !------------------------------------------------------------------
  !! Prints average, max and min differences between XC arrays
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(IN) :: what
  REAL(DP), INTENT(IN) :: vaver(2), vmax(2), vmin(2)
  REAL(DP), OPTIONAL :: averref
  !
  IF (test=='dft-comparison') THEN
    WRITE(stdout,*) " "
    WRITE(stdout,*) " ", TRIM(what)
    WRITE(stdout,*) "AVR abs: ", vaver(1), "   AVR %: ", vaver(2)
    WRITE(stdout,*) "MAX abs: ", vmax(1),  "   MAX %: ", vmax(2)
    WRITE(stdout,*) "MIN abs: ", vmin(1),  "   MIN %: ", vmin(2)
  ELSEIF (test=='exe-benchmark') THEN
    WRITE(stdout,*) " "
    WRITE(stdout,*) " ", TRIM(what)
    WRITE(stdout,*) "AVR test: ", vaver(1)
    WRITE(stdout,*) "AVR ref : ", averref
    WRITE(stdout,*) "diff    : ", vaver(1)-averref
  ENDIF
  !
 END SUBROUTINE print_stat
 !
 !------------------------------------------------------------------
 SUBROUTINE print_diff( what, x_dft1, x_dft2, x_ud1, x_ud2 )
  !-----------------------------------------------------------------
  !! Prints difference between dft1 and dft2 output arrays.
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: what
  REAL(DP), INTENT(IN) :: x_dft1(ns), x_dft2(ns)
  REAL(DP), INTENT(IN), OPTIONAL :: x_ud1, x_ud2
  !
  IF (test=='gen-benchmark') RETURN
  !
  WRITE(stdout,*) " "
  WRITE(stdout,*) what
  !
  IF ( .NOT. POLARIZED .OR. what(1:1)=='E' ) THEN
    WRITE(stdout,101) x_dft1(1)
    WRITE(stdout,201) x_dft2(1)
    WRITE(stdout,*) " --- "
    WRITE(stdout,301) ABS(x_dft1(1)-x_dft2(1))
  ELSEIF ( POLARIZED ) THEN
    IF ( .NOT. PRESENT(x_ud1) ) THEN
      WRITE(stdout,102) x_dft1(1), x_dft1(2)
      WRITE(stdout,202) x_dft2(1), x_dft2(2)
      WRITE(stdout,*) " --- "
      WRITE(stdout,302) ABS(x_dft1(1)-x_dft2(1)), ABS(x_dft1(2)-x_dft2(2))
    ELSE
      WRITE(stdout,103) x_dft1(1), x_ud1, x_dft1(2)
      WRITE(stdout,203) x_dft2(1), x_ud2, x_dft2(2)
      WRITE(stdout,*) " --- "
      WRITE(stdout,303) ABS(x_dft1(1)-x_dft2(1)), ABS(x_ud1-x_ud2), &
                        ABS(x_dft1(2)-x_dft2(2))
    ENDIF
  ENDIF
  !
  101 FORMAT('dft1/test: ',F17.14)
  102 FORMAT('dft1/test: ',F17.14,4x,F17.14)
  103 FORMAT('dft1/test: ',F17.14,4x,F17.14,4x,F17.14)
  201 FORMAT('dft2/ref:  ',F17.14)
  202 FORMAT('dft2/ref:  ',F17.14,4x,F17.14)
  203 FORMAT('dft2/ref:  ',F17.14,4x,F17.14,4x,F17.14)
  301 FORMAT('diff: ',5x,F17.14)
  302 FORMAT('diff: ',5x,F17.14,4x,F17.14)
  303 FORMAT('diff: ',5x,F17.14,4x,F17.14,4x,F17.14)
  !
 END SUBROUTINE print_diff
 !
 !----------------------------------------------------------------------
 SUBROUTINE print_diff2( what, dxc1, dxc2 )
  !---------------------------------------------------------------------
  !! Same as \(\texttt{print_diff}\), but for derivatives of Vxc.
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: what
  REAL(DP), INTENT(IN) :: dxc1(ns,ns), dxc2(ns,ns)
  !
  IF (test=='gen-benchmark') RETURN
  !
  WRITE(stdout,*) " "
  WRITE(stdout,*) what
  !
  IF ( .NOT. POLARIZED ) THEN
    WRITE(stdout,101) dxc1(1,1)
    WRITE(stdout,201) dxc2(1,1)
    WRITE(stdout,*) " --- "
    WRITE(stdout,301) dxc1(1,1)-dxc2(1,1)
  ELSE
    WRITE(stdout,103) dxc1(1,1), dxc1(2,1), dxc1(2,2)
    WRITE(stdout,203) dxc2(1,1), dxc2(2,1), dxc2(2,2)
    WRITE(stdout,*) " --- "
    WRITE(stdout,303) dxc1(1,1)-dxc2(1,1), dxc1(2,1)-dxc2(2,1), &
                      dxc1(2,2)-dxc2(2,2)
  ENDIF
  !
  101 FORMAT('dft1/test: ',F17.14)
  103 FORMAT('dft1/test: ',F17.14,4x,F17.14,4x,F17.14)
  201 FORMAT('dft2/ref:  ',F17.14)
  203 FORMAT('dft2/ref:  ',F17.14,4x,F17.14,4x,F17.14)
  301 FORMAT('diff: ',4x,F17.14)
  303 FORMAT('diff: ',4x,F17.14,4x,F17.14,4x,F17.14)
  !
 END SUBROUTINE print_diff2
 !
 !-------------------------------------------------------------------------
 SUBROUTINE evxc_stats( what, xc_1, xc_2, aver )
  !------------------------------------------------------------------------
  !! If test=dft-comparison calculates average, max and min difference 
  !! between output arrays of dft1 and dft2.  
  !! If test=exe-benchmark calculates difference between total energy
  !! and potential calculated over npoints k-points and values taken from
  !! benchmark data file.  
  !! If test=gen-benchmark calculates the total energy and potential 
  !! over npoints k-points.
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: what
  REAL(DP), INTENT(IN) :: xc_1(nnrt,ns)
  REAL(DP), INTENT(IN) :: xc_2(nnrt,ns)
  REAL(DP), INTENT(INOUT) :: aver(2)
  !
  REAL(DP) :: xc_aver(2,2), xc_max(2,2), xc_min(2,2)
  REAL(DP) :: thr, aver_snd(1), aver_rec(1)
  INTEGER :: ierr2
  !
  IF (LDA .AND. .NOT.GGA) THEN
    IF (what(1:1)=='E') thr = diff_thr_e_lda
    IF (what(1:1)=='V') thr = diff_thr_v_lda
  ELSEIF ( GGA ) THEN
    IF (what(1:1)=='E') thr = diff_thr_e_lda
    IF (what(1:1)=='V') thr = diff_thr_v_lda
  ELSEIF ( MGGA ) THEN
    IF (what(1:1)=='E') thr = diff_thr_e_mgga
    IF (what(1:1)=='V') thr = diff_thr_vmgga
  ENDIF
  !
  IF (test/='gen-benchmark') THEN
    WRITE(stdout,*) " "
    IF ( POLARIZED .AND. what(1:1)/='E' ) WRITE(stdout,*) TRIM(what)
  ENDIF
  !
  xc_aver=0._DP ; xc_max=0._DP ; xc_min=0._DP
  !
  IF ( test=='dft-comparison' ) THEN
    CALL diff_average( thr, xc_1(1:nnr,1), xc_2(1:nnr,1), xc_aver(:,1),nnr_int)
    CALL diff_max( thr, xc_1(1:nnr,1), xc_2(1:nnr,1), xc_max(:,1) )
    CALL diff_min( thr, xc_1(1:nnr,1), xc_2(1:nnr,1), xc_min(:,1) )
  ELSE
    xc_aver(1,1) = SUM(xc_1(1:nnr,1))/DBLE(npoints)
    xc_max(1,1) = MAXVAL( xc_1(1:nnr,1) )
    xc_min(1,1) = MINVAL( xc_1(1:nnr,1) )
    !
#if defined(__MPI)
    aver_snd = xc_aver(1,1)
    CALL MPI_REDUCE( aver_snd, aver_rec, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, &
                     comm, ierr2 )
    xc_aver(1:1,1) = aver_rec
#endif
  ENDIF
  !
  IF ( .NOT. POLARIZED .OR. what(1:1)=='E' ) THEN
     IF (mype==root) THEN
       IF (test=='dft-comparison')  THEN
         CALL print_stat( what, xc_aver(:,1), xc_max(:,1), xc_min(:,1) )
       ELSE  
         CALL print_stat( what, xc_aver(:,1),xc_max(:,1),xc_min(:,1),aver(1) )
       ENDIF
     ENDIF
  ELSE
    IF ( test=='dft-comparison' ) THEN
      CALL diff_average(thr,xc_1(1:nnr,2),xc_2(1:nnr,2),xc_aver(:,2),nnr_int)
      CALL diff_max( thr, xc_1(1:nnr,2), xc_2(1:nnr,2), xc_max(:,2) )
      CALL diff_min( thr, xc_1(1:nnr,2), xc_2(1:nnr,2), xc_min(:,2) )
      !
      IF (TRIM(what)=='V2c' .AND. GGA ) THEN
        CALL diff_average( diff_thr_vgga, v2c_ud1, v2c_ud2, vaver, nnr_int )
        CALL diff_max( diff_thr_vgga, v2c_ud1, v2c_ud2, vmax )
        CALL diff_min( diff_thr_vgga, v2c_ud1, v2c_ud2, vmin )
        !
        IF (mype==root) CALL print_stat( 'cross', vaver, vmax, vmin )
      ENDIF
      !      
    ELSE
      xc_aver(1,2) = SUM(xc_1(1:nnr,2))/npoints
      xc_max(1,2) = MAXVAL( xc_1(1:nnr,2) )
      xc_min(1,2) = MINVAL( xc_1(1:nnr,2) )
      !
#if defined(__MPI)
      aver_snd = xc_aver(1,2)
      CALL MPI_REDUCE( aver_snd, aver_rec, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       0, comm, ierr2 )
      xc_aver(1:1,2) = aver_rec
#endif
      !
      IF (TRIM(what)=='V2c' .AND. GGA ) THEN
        v2c_ud1_aver(1) = SUM(v2c_ud1(1:nnr))/npoints
        v2c_ud1_max(1) = MAXVAL( v2c_ud1(1:nnr) )
        v2c_ud1_min(1) = MINVAL( v2c_ud1(1:nnr) )
        ! 
#if defined(__MPI)
        aver_sndu = v2c_ud1_aver(1)
        CALL MPI_REDUCE( aver_sndu, aver_recu, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, 0, comm, ierrm )
        v2c_ud1_aver(1) = aver_recu
#endif
        !
        IF (mype==root) CALL print_stat( 'cross', v2c_ud1_aver, v2c_ud1_max, &
                                         v2c_ud1_min, v2c_aver_b(1,3) )
      ENDIF
      !
    ENDIF
    !
    IF (mype==root) THEN
      CALL print_stat( 'up',  xc_aver(:,1),xc_max(:,1),xc_min(:,1), aver(1) )
      CALL print_stat( 'down',xc_aver(:,2),xc_max(:,2),xc_min(:,2), aver(2) )
    ENDIF 
    !
  ENDIF
  !
  IF (test=='gen-benchmark') THEN
     aver = xc_aver(1,:)
     IF (TRIM(what)=='V2c'.AND.GGA.AND.ns==2) v2c_aver_b(1,3)=v2c_ud1_aver(1)
  ENDIF
  !
  RETURN
  !
 END SUBROUTINE evxc_stats
 !
 !
 !---------------------------------------------------------------------
 SUBROUTINE derivxc_stats( what, dxc_qe, dxc_lxc, aver )
  !--------------------------------------------------------------------
  !! Same as \(\texttt{evxc_stats}\), but for derivatives of Vxc.
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: what
  REAL(DP), INTENT(IN) :: dxc_qe(nnrt,ns,ns)
  REAL(DP), INTENT(IN) :: dxc_lxc(nnrt,ns,ns)
  REAL(DP), INTENT(INOUT) :: aver(3)
  !
  REAL(DP) :: dxc_aver(2,np), dxc_max(2,np), dxc_min(2,np)
  REAL(DP) :: thr, aver_snd(1), aver_rec(1)
  INTEGER :: ierr2
  !
  IF (LDA .AND. .NOT.GGA) thr = diff_thr_dmuxc
  IF ( GGA ) thr = diff_thr_dv
  !
  IF (test/='gen-benchmark') THEN
    WRITE(stdout,*) " "
    WRITE(stdout,*) what
  ENDIF
  !
  dxc_aver=0._DP ; dxc_min=0._DP ; dxc_max=0._DP
  !
  IF ( test=='dft-comparison' ) THEN
    CALL diff_average( thr, dxc_qe(1:nnr,1,1), dxc_lxc(1:nnr,1,1), &
                       dxc_aver(1:nnr,1), nnr_int )
    CALL diff_max( thr, dxc_qe(1:nnr,1,1), dxc_lxc(1:nnr,1,1),     &
                   dxc_max(1:nnr,1) )
    CALL diff_min( thr, dxc_qe(1:nnr,1,1), dxc_lxc(1:nnr,1,1),     &
                   dxc_min(1:nnr,1) )
  ELSE
    dxc_aver(1,1) = SUM(dxc_qe(1:nnr,1,1))/DBLE(npoints)
    dxc_max(1,1) = MAXVAL( dxc_qe(1:nnr,1,1) )
    dxc_min(1,1) = MINVAL( dxc_qe(1:nnr,1,1) )
    !
#if defined(__MPI)
    aver_snd = dxc_aver(1,1)
    CALL MPI_REDUCE( aver_snd, aver_rec, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                     0, comm, ierr2 )
    dxc_aver(1:1,1) = aver_rec
#endif
  ENDIF
  !
  IF ( .NOT. POLARIZED ) THEN
    IF (mype==root) THEN
      IF (test=='dft-comparison') THEN
        CALL print_stat( what, dxc_aver(1:nnr,1), dxc_max(1:nnr,1),   &
                         dxc_min(1:nnr,1) )
      ELSE
        CALL print_stat( what, dxc_aver(1:nnrb,1), dxc_max(1:nnrb,1), &
                         dxc_min(1:nnrb,1), aver(1) )
      ENDIF  
    ENDIF  
  ELSE
    IF ( test=='dft-comparison' ) THEN
      CALL diff_average( thr, dxc_qe(1:nnr,1,2), dxc_lxc(1:nnr,1,2), &
                         dxc_aver(1:nnr,2), nnr_int )
      CALL diff_max( thr, dxc_qe(1:nnr,1,2), dxc_lxc(1:nnr,1,2),     &
                     dxc_max(1:nnr,2) )
      CALL diff_min( thr, dxc_qe(1:nnr,1,2), dxc_lxc(1:nnr,1,2),     &
                     dxc_min(1:nnr,2) )
      !
      CALL diff_average( thr, dxc_qe(1:nnr,2,2), dxc_lxc(1:nnr,2,2), &
                         dxc_aver(1:nnr,3), nnr_int )
      CALL diff_max( thr, dxc_qe(1:nnr,2,2), dxc_lxc(1:nnr,2,2),     &
                     dxc_max(1:nnr,3) )
      CALL diff_min( thr, dxc_qe(1:nnr,2,2), dxc_lxc(1:nnr,2,2),     &
                     dxc_min(1:nnr,3) )
    ELSE
      dxc_aver(1,2) = SUM(dxc_qe(1:nnr,1,2))/DBLE(npoints)
      dxc_max(1,2) = MAXVAL( dxc_qe(1:nnr,1,2) )
      dxc_min(1,2) = MINVAL( dxc_qe(1:nnr,1,2) )
      !
#if defined(__MPI)
      aver_snd = dxc_aver(1,2)
      CALL MPI_REDUCE( aver_snd, aver_rec, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, 0, comm, ierr2 )
      dxc_aver(1:1,2) = aver_rec
#endif
      dxc_aver(1,3) = SUM(dxc_qe(1:nnr,2,2))/DBLE(npoints)
      dxc_max(1,3) = MAXVAL( dxc_qe(1:nnr,2,2) )
      dxc_min(1,3) = MINVAL( dxc_qe(1:nnr,2,2) )
      !
#if defined(__MPI)
      aver_snd = dxc_aver(1,3)
      CALL MPI_REDUCE( aver_snd, aver_rec, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, 0, comm, ierr2 )
      dxc_aver(1:1,3) = aver_rec
#endif
    ENDIF
    !
    IF (mype==root) THEN
      CALL print_stat( 'up-up', dxc_aver(:,1), dxc_max(:,1), dxc_min(:,1),  &
                       aver(1) )
      CALL print_stat( 'up-down', dxc_aver(:,2), dxc_max(:,2), dxc_min(:,2),&
                       aver(2) )
      CALL print_stat( 'down-down',dxc_aver(:,3), dxc_max(:,3),dxc_min(:,3),&
                       aver(3) )
    ENDIF  
    !
  ENDIF
  !
  IF (test=='gen-benchmark') aver(1:np) = dxc_aver(1,:)
  !
  RETURN
  !
 END SUBROUTINE derivxc_stats
 !
 !
 !------------------------------------------------------------------------
 FUNCTION is_it_out( diff_thr, dm, x_dft1, x_dft2, x_ud_1, x_ud_2 )
  !----------------------------------------------------------------------
  !! TRUE if the difference between \(\text{x_df1}\) and \(\text{x_df2}\)
  !! is bigger than the difference threshold.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: dm
  REAL(DP), INTENT(IN) :: diff_thr, x_dft1(dm), x_dft2(dm)
  REAL(DP), INTENT(IN), OPTIONAL :: x_ud_1, x_ud_2
  LOGICAL :: is_it_out, is_it_out_ud
  !
  is_it_out = ANY(ABS(x_dft1(1:dm)-x_dft2(1:dm)) > diff_thr)
  !
  IF (PRESENT(x_ud_1)) THEN
    is_it_out_ud =  ABS(x_ud_1-x_ud_2) > diff_thr
    is_it_out = ANY( (/ is_it_out, is_it_out_ud /) )
  ENDIF
  !
 END FUNCTION
 !
 !--------------------------------------------------------------------------
 FUNCTION is_dit_out( diff_thr, dx_dft1, dx_dft2 )
  !------------------------------------------------------------------------
  !! TRUE if the difference between \(\text{dx_df1}\) and \(\text{dx_df2}\)
  !! is bigger than the difference threshold.
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: diff_thr, dx_dft1(ns,ns), dx_dft2(ns,ns)
  REAL(DP) :: dxc_diff(np)
  LOGICAL :: is_dit_out
  !
  dxc_diff(1) = ABS(dx_dft1(1,1)-dx_dft2(1,1))
  IF ( POLARIZED ) THEN
    dxc_diff(2) = ABS(dx_dft1(2,1)-dx_dft2(2,1))
    dxc_diff(3) = ABS(dx_dft1(2,2)-dx_dft2(2,2))
  ENDIF
  !
  is_dit_out = ANY(dxc_diff(:) > diff_thr)
  !
 END FUNCTION
 !
 !------------------------------------------------------------------------
 FUNCTION calc_perc_diff( thr, x_qe, x_lxc )
  !----------------------------------------------------------------------
  !! Calculates difference between qe and libxc quantities in percentage.
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: thr
  REAL(DP), INTENT(IN) :: x_qe, x_lxc
  REAL(DP) :: calc_perc_diff
  !
  REAL(DP) :: perc_diff
  !
  perc_diff = -1.d0
  IF ( ABS(x_qe)<10.d0*thr .AND. ABS(x_qe-x_lxc)<10.d0*thr ) RETURN
  IF ( ABS(x_qe)==0.d0 .AND. ABS(x_qe-x_lxc)>thr ) calc_perc_diff = 100.d0
  IF ( ABS(x_qe)>thr ) calc_perc_diff = ABS( (x_qe-x_lxc)/x_qe )*100.d0
  !
  RETURN
  !
 END FUNCTION calc_perc_diff
 !
END PROGRAM xclib_test
