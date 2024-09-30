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
  !! Testing program for xc\_lib library in QE:
  !
  !! * generate: stores a set of dft benchmark data on xml file;
  !! * execute:  runs XClib to calculate data then compares them 
  !!             to the benchmark set in the xml file.
  !
  !! Choices:
  !
  !! * ALL_TERMS: run test over all available terms for each family
  !!              and kind;
  !! * ALL_SHORT: run test over all available full dfts in QE internal
  !!              dft list;
  !! * ALL_LIBXC: run test over all Libxc dfts usable in QE.
  !
  !! Other specific options:
  !
  !! * family (LDA, GGA, MGGA);
  !! * polarization;
  !! * derivative of xc potential.
  !
  !! See README.TEST file for more details.
  !
  USE kind_l,         ONLY: DP
  USE constants_l,    ONLY: pi
  USE beef_interface, ONLY: beefsetmode
  USE xc_lib,         ONLY: xclib_set_dft_from_name, xclib_set_exx_fraction, &
                            xclib_get_ID, xclib_reset_dft, xc, xc_gcx,       &
                            xc_metagcx, xclib_dft_is_libxc, xclib_init_libxc,&
                            xclib_finalize_libxc, xclib_set_finite_size_volume,&
                            xclib_set_auxiliary_flags, xclib_dft_is, start_exx,&
                            set_libxc_ext_param, dmxc
  USE xclib_utils_and_para
  !--xml
  USE xmltools,       ONLY: xml_open_file, xml_closefile,xmlr_readtag,  &
                            xmlw_writetag, xmlw_opentag, xmlw_closetag, &
                            xmlr_opentag, xmlr_closetag, get_attr, add_attr
#if defined(__LIBXC)
  USE xc_f03_lib_m
  USE dft_setting_params, ONLY: xc_func, xc_info, xc_kind_error, libxc_flags
#endif
  USE dft_setting_params,   ONLY: is_libxc
  USE dft_setting_routines, ONLY: capital
  !
  USE qe_dft_list, ONLY: nxc, ncc, ngcx, ngcc, nmeta, n_dft, &
                         dft_LDAx_name, dft_LDAc_name, dft_GGAx_name, &
                         dft_GGAc_name, dft_MGGA_name, dft_full
  USE qe_dft_refs
  !
  IMPLICIT NONE
  !
#include "qe_version.h"  
  !
#if defined(__MPI)
  INTEGER    STATUS(MPI_STATUS_SIZE)
#else
#define MPI_MAX_PROCESSOR_NAME 64
#endif
  !
#if defined(__LIBXC)
  INTEGER :: major, minor, micro, ifamily, fkind
  TYPE(xc_f03_func_t) :: xc_func0
  TYPE(xc_f03_func_info_t) :: xc_info0
#endif
  !
  INTEGER :: mype, npes, comm, ntgs, root
  LOGICAL :: iope
  INTEGER :: i, ierr, ierrm
  INTEGER :: nnodes, nlen
  !
  INTEGER, PARAMETER :: stdin=5
  REAL(DP) :: time_tot1, time_tot2, time(6)=0.d0
  !
  !-------- Grid dim vars --------------------
  INTEGER, PARAMETER :: npoints=90000
  INTEGER :: nr, nnr, nrpe, nnr_b, nnrbt
  INTEGER :: is, is_min, is_max
  !
  !-------- Input vars -----------------------
  CHARACTER(LEN=15) :: test, family, fam_init
  CHARACTER(LEN=32) :: dft, dft_init, dft_lxc, xc_kind, xmldft, xmlfamily
  CHARACTER(LEN=30) :: polarization, xmlpolarization
  CHARACTER(LEN=15) :: input_err_where=''
  CHARACTER(LEN=40) :: input_err=''
  LOGICAL :: xc_derivative, xmlxc_derivative, show_time
  !
  !---------- DFT infos -------------------------
  INTEGER :: iexch1, icorr1, igcx1, igcc1, imeta1, imetac1
  INTEGER :: id_vec(6), n_qe_func, naver
  LOGICAL :: LDA, GGA, MGGA, POLARIZED
  !
  !-------- Various params -------------------
  REAL, PARAMETER :: volume=0.1d0
  REAL(DP), PARAMETER :: null=0.0_DP, pi34=0.6203504908994_DP
  REAL(DP), PARAMETER :: thresh_lda  = 0.d0, & !1.E-6_DP, &
                         thresh_gga  = 0.d0, & !1.E-6_DP, &
                         thresh_mgga = 0.d0    !1.E-6_DP
  REAL(DP), PARAMETER :: diff_thr_e_lda  = 1.0E-8_DP,  &
                         diff_thr_v_lda  = 1.0E-8_DP,  &
                         diff_thr_e_gga  = 1.0E-10_DP, &
                         diff_thr_vgga   = 1.0E-10_DP, &
                         diff_thr_e_mgga = 1.0E-10_DP, &
                         diff_thr_vmgga  = 1.0E-10_DP, &
                         diff_thr_dmuxc  = 1.0E-6_DP,  &
                         diff_thr_dv     = 1.0E-10_DP
  REAL(DP) :: fact
  REAL(DP) :: aver_thresh = 10.E-8        !^^--- to optimize
  !
  !---------- Loop indexes ---------------------------
  INTEGER :: id, ii, ns, np, nthr, iip, iout, iaverout, iavernull, &
             l, idterm
  !
  !---------- XC-input vars ------------------
  REAL(DP), ALLOCATABLE :: rho(:,:), rho_tz(:,:)
  REAL(DP), ALLOCATABLE :: grho(:,:,:), grh(:,:,:)
  REAL(DP), ALLOCATABLE :: tau(:,:)
  REAL(DP) :: grho2(2), grho_ud
  REAL(DP) :: rhoi(2), grhoi(3,2), taui(2)
  !
  !--------- set1 vars --------------------------
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
  !--------- set2 vars ---------------------------
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
  !----------match vars ---------------------------
  LOGICAL :: ex_is_out, ec_is_out, vx_is_out, vc_is_out, &
             something_out, dmuxc_is_out
  LOGICAL :: v1x_is_out, v2x_is_out, v1c_is_out, v2c_is_out, &
             v3x_is_out, v3c_is_out, &
             dvxcrr_is_out, dvxcsr_is_out, dvxcss_is_out, dvgga_is_out
  ! ... LDA aver
  REAL(DP) :: ex_aver(2),   ec_aver(2),   &
              vx_aver(1,2), vc_aver(1,2), &
              dv_aver(3)
  ! ... GGA/MGGA aver
  REAL(DP) :: v1x_aver(1,2), v1c_aver(1,2), &
              v2x_aver(1,2), v2c_aver(1,3), &
              v2c_ud1_aver(2), &
              dvrr_aver(1,3), dvsr_aver(1,3), &
              dvss_aver(1,3)
  ! ... MGGA aver
  REAL(DP) :: v3x_aver(1,2), v3c_aver(1,2)
  REAL(DP) :: aver_sndu, aver_recu
  !
  ! ... xml
  INTEGER :: tag_err
  LOGICAL :: not_found_skip
  CHARACTER(LEN=1) :: dummy
  CHARACTER(LEN=30) :: filename_xml=""
  CHARACTER(LEN=48) :: xc_data="XC_DATA__________"
  ! ... output
  INTEGER :: iunpun, iun, nlen1, nlen2
  LOGICAL :: found, exc_term=.TRUE., cor_term=.TRUE.
  CHARACTER(LEN=40), PARAMETER :: failed='**FAILED**', &
                                  skipped='**skipped - by default**', &
                                  skipped2='**skipped - not found in xml**',&
                                  skipped3='**skipped - needs Libxc**',&
                                  skipped4='**skipped - Libxc dft not usable in QE**'
  CHARACTER(LEN=18), PARAMETER :: passed='match', stored='stored', &
                                  passed0='match (but null!)'
  CHARACTER(LEN=6) :: gen_version = ''
  CHARACTER(LEN=5) :: libxc_version='none', libxc_gen_version = ''
  !
  ! ... MPI/OpenMP stuff
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
  NAMELIST/input_namelist/ test, filename_xml, polarization, family, &
                           xc_derivative, dft, show_time
  !
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
  dft = 'none'
  xc_derivative = .FALSE.
  polarization = 'none'
  !
  nowarning = .TRUE.
  !
  !==========================================================================
  ! GET INPUT FROM FILE
  !==========================================================================
  !
  IF (mype==root) THEN
    !
    READ( stdin, input_namelist )
    !
    DO i = 1, LEN_TRIM(test)
      test(i:i) = capital( test(i:i) )
    ENDDO
    DO i = 1, LEN_TRIM(dft)
      dft(i:i) = capital( dft(i:i) )
    ENDDO
    DO i = 1, LEN_TRIM(family)
      family(i:i) = capital( family(i:i) )
    ENDDO
    DO i = 1, LEN_TRIM(polarization)
      polarization(i:i) = capital( polarization(i:i) )
    ENDDO
    IF ( dft(1:4)=='ALL_' ) family = 'ALL'
    !
    IF ( test=='GENERATE' ) THEN
      !
      iunpun = xml_open_file( "./"//TRIM(filename_xml) )
      IF ( iunpun == -1 ) THEN
        CALL xclib_error( 'xclib_test', 'GENERATE test cannot open xml a data file', 1 )
      ENDIF
      !
      CALL xmlw_opentag( "XCTEST-DATA-SET" )
      CALL add_attr( "DFT", dft )
      CALL add_attr( "FAMILY", family )
      CALL add_attr( "VXC_DERIVATIVE", xc_derivative )
      CALL add_attr( "SPIN", polarization )
      CALL xmlw_writetag( "HEADER", "" )
      !
    ELSEIF ( TRIM(test)=='EXECUTE' ) THEN
      !
      INQUIRE( FILE = filename_xml, exist=found )
      IF (.NOT. found ) THEN
        CALL xclib_error( 'xclib_test', 'xml data file not found', 1 )
      ENDIF
      !
      iun = xml_open_file( filename_xml )
      IF ( iun==-1 ) THEN
        CALL xclib_error( 'xclib_test', 'xml data file not readable', 2 )
      ENDIF
      !
      CALL xmlr_opentag( "XCTEST-DATA-SET" )
      !
      CALL xmlr_readtag( "HEADER", dummy )
      CALL get_attr( "DFT", xmldft )
      CALL get_attr( "FAMILY", xmlfamily )
      CALL get_attr( "VXC_DERIVATIVE", xmlxc_derivative )
      CALL get_attr( "SPIN", xmlpolarization )
      !
      IF ( dft(1:4)=='ALL_' .AND. TRIM(dft)/=TRIM(xmldft) ) input_err_where = 'dft'
      IF ( TRIM(xmlpolarization)/='BOTH' .AND. polarization/=xmlpolarization ) &
                                                   input_err_where = 'polarization'
      IF ( xmlxc_derivative .NEQV. xc_derivative ) input_err_where = 'xc_derivative'
      IF ( TRIM(input_err_where) /= '' ) THEN
        input_err = 'Input mismatch in: '
        WRITE(input_err(20:35), '(a)') input_err_where
        CALL xclib_error( 'xclib_test', input_err, 1 )
      ENDIF
      IF ( TRIM(family) /= TRIM(xmlfamily) ) THEN
        WRITE(stdout,*) 'WARNING: input mismatch in FAMILY. The one in the xml &
                        &file will be considered only.'
        family = xmlfamily
      ENDIF
      !
    ELSE
      !
      CALL xclib_error( 'xclib_test', 'Wrong input test.', 2 )
      !
    ENDIF
  ENDIF
  !
  !==========================================================================
  ! MPI MANAGEMENT
  !==========================================================================
  !
#if defined(__MPI)
  CALL MPI_BCAST( test,         30, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST( family,       30, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST( dft,          30, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST( polarization, 30, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST( xc_derivative, 1, MPI_LOGICAL,   0, MPI_COMM_WORLD, ierr )
#endif
  !
  ALLOCATE( proc_name(npes) )
  ALLOCATE( node_name(npes) )
  ALLOCATE( proc2node(npes) )
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
  IF ( mype == root .AND. test/='GENERATE' ) THEN
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
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
#endif
  !
  ! ... point distribution over CPUs
  !
  nr = 0
  IF ( MOD(npoints,npes)<=(mype+1) ) nr = 1
  nnr = npoints/npes + ABS(nr-1)
  IF (nr==1) nrpe = mype*nnr
  IF (nr==0) nrpe = npoints-(npoints/npes)*(npes-mype-1)
  !
  ! ... openacc init (otherwise it offsets the wall time of the first test)
  !
#if defined(_OPENACC)
  !$acc data create( time )
  !$acc end data
#endif
  !
  ! ... capitalize input vars
  !
  IF (TRIM(polarization)=='UNPOLARIZED' ) THEN
    is_min = 1  ;  is_max = 1
  ELSEIF (TRIM(polarization)=='POLARIZED' ) THEN
    is_min = 2  ;  is_max = 2
  ELSEIF (TRIM(polarization)=='BOTH' ) THEN
    is_min = 1  ;  is_max = 2
  ELSE
    CALL xclib_error( 'xclib_test', 'Wrong input polarization.', 3 )
  ENDIF
  !
  fam_init = family
  dft_init = dft
  !
  n_qe_func = 1
  IF (TRIM(dft)=='ALL_TERMS') n_qe_func = nxc+ncc+ngcx+ngcc+nmeta+5
  IF (TRIM(dft)=='ALL_SHORT') n_qe_func = n_dft
  IF (TRIM(dft)=='ALL_LIBXC') n_qe_func = 999
  !
  ! ... check QE and Libxc versions used for storing and the current ones
  !
  IF ( mype==root ) THEN
#if defined(__LIBXC)
    CALL xc_f03_version( major, minor, micro )
    libxc_version = ' . . '
    WRITE(libxc_version(1:1), '(i1)') major
    WRITE(libxc_version(3:3), '(i1)') minor
    WRITE(libxc_version(5:5), '(i1)') micro
#else
    IF (dft=='ALL_LIBXC') CALL xclib_error( 'xclib_test','Libxc library not linked.',4 )
#endif
    IF ( TRIM(test)=='GENERATE' ) THEN
      CALL xmlw_opentag( "QE_Libxc_VERSION_TEST" )
      CALL xmlw_writetag( "QE_version", version_number )
      CALL xmlw_writetag( "Libxc_version", libxc_version )
      CALL xmlw_closetag()
    ELSEIF ( TRIM(test)=='EXECUTE' ) THEN
      CALL xmlr_opentag( "QE_Libxc_VERSION_TEST" )
      CALL xmlr_readtag( "QE_version", gen_version )
      CALL xmlr_readtag( "Libxc_version", libxc_gen_version )
      CALL xmlr_closetag()
      !
      WRITE(stdout,*) 'QE version used for data storing: ', gen_version
      WRITE(stdout,*) 'QE version currently in use: ', version_number
      WRITE(stdout,*) ''
      WRITE(stdout,*) 'Libxc version used for data storing: ', libxc_gen_version
      WRITE(stdout,*) 'Libxc version currently in use: ', libxc_version
      WRITE(stdout,*) ''
    ENDIF
  ENDIF
  !
  ! ... main loop over functionals to test
  !
  DO id = 1, n_qe_func
   DO is = is_min, is_max
    !
    time = 0.d0
    ns = is
    !
    CALL xclib_reset_dft()
    !
    IF ( dft_init=='ALL_TERMS' ) THEN
      IF (id<=nxc+1) THEN
         idterm = id-1
         dft = dft_LDAx_name(idterm)
         xc_kind = 'X'
      ELSEIF (id>=nxc+2 .AND. id<=nxc+ncc+2) THEN
         idterm = id-nxc-2
         dft = dft_LDAc_name(idterm)
         xc_kind = 'C'
      ELSEIF (id>=nxc+ncc+3 .AND. id<=nxc+ncc+ngcx+3) THEN
         idterm = id-nxc-ncc-3
         dft = dft_GGAx_name(idterm)
         xc_kind = 'X'
      ELSEIF (id>=nxc+ncc+ngcx+4 .AND. id<=nxc+ncc+ngcx+ngcc+4) THEN
         idterm = id-nxc-ncc-ngcx-4
         dft = dft_GGAc_name(idterm)
         xc_kind = 'C'
      ELSEIF (id>=nxc+ncc+ngcx+ngcc+5 .AND. id<=nxc+ncc+ngcx+ngcc+nmeta+5) THEN
         idterm = id-nxc-ncc-ngcx-ngcc-5
         dft = dft_MGGA_name(idterm)
         xc_kind = 'XC'
      ENDIF
    ELSEIF ( dft_init=='ALL_SHORT' ) THEN
      dft = dft_full(id)%name
    ENDIF
    !
    ! ... initialization of averages
    !
    ex_aver   = 0._DP  ; ec_aver   = 0._DP
    vx_aver   = 0._DP  ; vc_aver   = 0._DP
    v1x_aver  = 0._DP  ; v2x_aver  = 0._DP
    v3x_aver  = 0._DP  ; v1c_aver  = 0._DP
    v2c_aver  = 0._DP  ; v3c_aver  = 0._DP
    v2c_ud1_aver = 0.0_DP
    dv_aver   = 0._DP  ; dvrr_aver = 0._DP
    dvsr_aver = 0._DP  ; dvss_aver = 0._DP
    !
    ! ... initialize first DFT
    !
    xc_data="XC_DATA__________"
    !
    ! ... overlaps between full dfts and shortnames
    !
    IF ( dft_init=='ALL_TERMS' ) THEN
      IF ( TRIM(dft)=='BLYP' ) dft='+BLYP'
      IF ( TRIM(dft)=='PZ'   ) dft='+PZ'
    ENDIF 
    !
    ! ... skipped cases (some need further checks)
    !
    IF ( TRIM(dft)=='xxxx' .OR. TRIM(dft)=='NONE' ) THEN
      id_vec(6)=idterm
      IF (mype==root) CALL print_test_status( skipped )
      CYCLE
    ENDIF
#if !defined(__LIBXC)
    IF ( TRIM(dft)=='TB09'  .OR. TRIM(dft)=='SCAN' .OR. &
         TRIM(dft)=='SCAN0' .OR. TRIM(dft)=='SCA0' .OR. &
         TRIM(dft)=='R2SCAN'.OR. TRIM(dft)=='RSCAN') THEN
      id_vec(6)=idterm
      IF (mype==root) CALL print_test_status( skipped3 )
      CYCLE
    ENDIF
#else
    ! libxc id=576 -> Libxc5.1.5 bug (same name of id=575)
    IF ( dft_init=='ALL_LIBXC' ) THEN
      IF ( id==576 ) THEN
        IF (mype==root) CALL print_test_status( skipped )
        CYCLE
      ENDIF
      !
      dft_lxc = xc_f03_functional_get_name( id )
      IF ( TRIM(dft_lxc) == '' ) CYCLE
      !
      fkind=-100 ; ifamily=-100
      CALL xc_f03_func_init( xc_func0, id, 1 )
      xc_info0 = xc_f03_func_get_info( xc_func0 )
      fkind = xc_f03_func_info_get_kind( xc_info0 )
      ifamily = xc_f03_func_info_get_family( xc_info0 )
      CALL xc_f03_func_end( xc_func0 )
      !
      dft = 'XC-000I-000I-000I-000I-000I-000I'
      IF ( ifamily==XC_FAMILY_LDA .AND. (fkind==XC_EXCHANGE .OR. fkind==XC_KINETIC) ) THEN
        WRITE( dft(4:6),   '(i3.3)' ) id
        WRITE( dft(7:7),   '(a)' ) 'L'
      ELSEIF ( ifamily==XC_FAMILY_LDA .AND. (fkind==XC_CORRELATION .OR. &
                                        fkind==XC_EXCHANGE_CORRELATION) ) THEN
        WRITE( dft(9:11),  '(i3.3)' ) id
        WRITE( dft(12:12), '(a)' ) 'L'
      ELSEIF ( (ifamily==XC_FAMILY_GGA .OR. ifamily==XC_FAMILY_HYB_GGA) .AND. &
               (fkind==XC_EXCHANGE .OR. fkind==XC_KINETIC) ) THEN
        WRITE( dft(14:16), '(i3.3)' ) id
        WRITE( dft(17:17), '(a)' ) 'L'
      ELSEIF ( (ifamily==XC_FAMILY_GGA  .OR. ifamily==XC_FAMILY_HYB_GGA) .AND. &
               (fkind==XC_CORRELATION .OR. fkind==XC_EXCHANGE_CORRELATION) ) THEN
        WRITE( dft(19:21), '(i3.3)' ) id
        WRITE( dft(22:22), '(a)' ) 'L'
      ELSEIF ( (ifamily==XC_FAMILY_MGGA .OR. ifamily==XC_FAMILY_HYB_MGGA) .AND. &
               (fkind==XC_EXCHANGE .OR. fkind==XC_KINETIC) ) THEN
        WRITE( dft(24:26), '(i3.3)' ) id
        WRITE( dft(27:27), '(a)' ) 'L'
      ELSEIF ( (ifamily==XC_FAMILY_MGGA .OR. ifamily==XC_FAMILY_HYB_MGGA) .AND. &
               (fkind==XC_CORRELATION .OR. fkind==XC_EXCHANGE_CORRELATION) ) THEN
        WRITE( dft(29:31), '(i3.3)' ) id
        WRITE( dft(32:32), '(a)' ) 'L'
      ENDIF
    ENDIF
#endif
    !
    CALL xclib_set_dft_from_name( dft )
    !
    LDA = .FALSE. ; GGA = .FALSE. ; MGGA = .FALSE.
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
    id_vec(1) = iexch1  ;  id_vec(2) = icorr1
    id_vec(3) = igcx1   ;  id_vec(4) = igcc1
    id_vec(5) = imeta1  ;  id_vec(6) = imetac1
    !
#if defined(__LIBXC)
    IF (xclib_dft_is_libxc( 'ANY' )) CALL xclib_init_libxc( ns, .FALSE. )
    !
    IF ( igcc1==428 .AND. is_libxc(4) ) THEN
!    IF ( igcx1==426 .AND. is_libxc(3) ) THEN
      ! Example of how to change an external parameter in a Libxc
      ! functional (HYB_GGA_XC_HSE06).
      ! Arguments:
      ! 1- family-kind index (1:LDAx, 2:LDAc, 3:GGAx, ...);
      ! 2- parameter index (you find it with xc_infos.x);
      ! 3- new value of the parameter.
      CALL set_libxc_ext_param( 4, 0, 0.25d0  )
      CALL set_libxc_ext_param( 4, 1, 0.106d0 )
      CALL set_libxc_ext_param( 4, 2, 0.106d0 )
      !
    ENDIF
    !
    IF ( xc_kind_error ) THEN
      CALL print_test_status( skipped4 )
      CALL xclib_finalize_libxc()
      CYCLE
    ENDIF
    IF ( xc_derivative .AND. ANY(libxc_flags(:,2)==0) ) THEN
      CALL print_test_status( skipped )
      CALL xclib_finalize_libxc()
      CYCLE
    ENDIF
    !
    IF (dft_init=='ALL_LIBXC') THEN
      fkind=-10
      DO l = 1, 6
        IF (id_vec(l)/=0.AND.is_libxc(l)) fkind=xc_f03_func_info_get_kind( xc_info(l) )
      ENDDO
      !
      IF (fkind==XC_EXCHANGE) THEN
        xc_kind = 'X'
      ELSEIF (fkind==XC_CORRELATION) THEN
        xc_kind = 'C'
      ELSEIF (fkind==XC_EXCHANGE_CORRELATION) THEN
        xc_kind = 'C'
      ENDIF
    ENDIF
#else
    IF (xclib_dft_is_libxc( 'ANY' )) THEN
      CALL print_test_status( skipped3 )
      CYCLE
    ENDIF
#endif
    !
    CALL xclib_set_auxiliary_flags( .FALSE. )
    IF ( xclib_dft_is('hybrid') ) CALL start_exx
    !
    IF ( fam_init(1:3)=='ALL' .AND. LDA  )  family='LDA'
    IF ( fam_init(1:3)=='ALL' .AND. GGA  )  family='GGA'
    IF ( fam_init(1:3)=='ALL' .AND. MGGA )  family='MGGA'
    !
    IF ( iexch1+icorr1+igcx1+igcc1+imeta1+imetac1==0 ) CYCLE
    !
    IF ( dft_init=='ALL_TERMS' ) THEN
      IF ( LDA .AND. GGA ) THEN
        IF ( id<=nxc+ncc+2 ) THEN
          GGA=.FALSE.
          IF (iexch1/=0 .AND. icorr1/=0 ) THEN
            IF ( TRIM(xc_kind)=='C') CYCLE
            IF ( TRIM(xc_kind)=='X') xc_kind='XC'
          ENDIF
          family='LDA'
        ELSE
          LDA=.FALSE.
          IF (igcx1/=0 .AND. igcc1/=0 ) THEN
            IF ( TRIM(xc_kind)=='C') CYCLE
            IF ( TRIM(xc_kind)=='X') xc_kind='XC'
          ENDIF
        ENDIF
      ENDIF
      IF ( LDA .AND. .NOT. GGA ) THEN
        IF ( id<=nxc+ncc+2 ) THEN
          IF (iexch1/=0 .AND. icorr1/=0 ) THEN
            IF ( TRIM(xc_kind)=='C') CYCLE
            IF ( TRIM(xc_kind)=='X') xc_kind='XC'
          ENDIF
          family='LDA'
        ENDIF
      ENDIF
      exc_term = xc_kind/='C'
      cor_term = xc_kind/='X'
    ELSE
      exc_term = (iexch1+igcx1+imeta1  /= 0)
      cor_term = (icorr1+igcc1+imetac1 /= 0)
      IF ( exc_term ) xc_kind='X'
      IF ( cor_term ) xc_kind='C'
      IF ( exc_term .AND. cor_term ) xc_kind='XC'
    ENDIF
    !
    IF (MGGA) THEN
      IF (.NOT. xc_derivative ) THEN
        xc_kind ='XC'
        exc_term = .TRUE.
        cor_term = .TRUE.
      ELSE
        CALL xclib_finalize_libxc()
        CYCLE
      ENDIF
    ENDIF
    !
    IF (ns == 2 .AND. (.NOT.is_libxc(2)) .AND. icorr1/=0 .AND. &
               icorr1/=1 .AND. icorr1/=2 .AND. icorr1/=4 .AND. &
               icorr1/=8 .AND. icorr1/=3 .AND. icorr1/=7 .AND. &
               icorr1/=13) THEN
      IF (mype==root) CALL print_test_status( skipped )         
      CYCLE
    ENDIF
    !
    IF (dft_init=='ALL_TERMS') THEN
      IF ( (LDA  .AND. exc_term .AND. dft_LDAx(iexch1)%wrn(1:12)=='never called') .OR. &
           (LDA  .AND. cor_term .AND. dft_LDAx(icorr1)%wrn(1:12)=='never called') .OR. &
           (GGA  .AND. exc_term .AND. dft_GGAx(igcx1)%wrn(1:12) =='never called') .OR. &
           (GGA  .AND. cor_term .AND. dft_GGAx(igcc1)%wrn(1:12) =='never called') .OR. &
           (MGGA .AND. exc_term .AND. dft_GGAx(imeta1)%wrn(1:12)=='never called') ) THEN
        !
        IF (mype==root) CALL print_test_status( skipped )
        CYCLE
        !
      ENDIF
    ENDIF
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
    nnrbt = nnr_b+nthr
    !
    np = 1
    IF (ns==2) np = 3
    !
    ! ... number of averages to calculate
    IF (.NOT.xc_derivative) THEN
      IF ( LDA ) naver = ns+1
      IF ( LDA .AND. TRIM(xc_kind)=='XC') naver = 2*ns+2
      IF ( GGA .AND. TRIM(xc_kind)=='X' ) naver = 2*ns+1
      IF ( GGA .AND. TRIM(xc_kind)=='C' ) naver = 3*ns+np
      IF ( GGA .AND. TRIM(xc_kind)=='XC') naver = 4*ns+np
      IF ( MGGA ) naver = 6*ns+2
    ELSE
      IF ( LDA ) naver = np
      IF ( GGA ) naver = 3*np
    ENDIF
    !
    !==========================================================================
    ! ALLOCATIONS OF XC I/O ARRAYS
    !==========================================================================
    !
    ! ... input
    !
    ALLOCATE( rho(nnr,ns), rho_tz(nnr,ns) )
    IF ( GGA .OR. MGGA ) ALLOCATE( grho(3,nnr,ns) )
    IF ( MGGA ) ALLOCATE( tau(nnr,ns) )
    !
    ! ... dft1 output arrays
    !
    IF (.NOT. xc_derivative) ALLOCATE( ex1(nnr), ec1(nnr) )
    !
    IF ( LDA .OR. GGA ) THEN
       IF ( LDA ) THEN
         IF ( .NOT. xc_derivative ) ALLOCATE( vx1(nnr,ns), vc1(nnr,ns) )
         IF ( xc_derivative ) ALLOCATE( dmuxc1(nnr,ns,ns) )
       ENDIF
       !
       IF ( GGA ) THEN
         IF ( .NOT. xc_derivative ) THEN
           ALLOCATE( exg1(nnr), ecg1(nnr) )
           ALLOCATE( v1x1(nnr,ns), v2x1(nnr,ns) )
           ALLOCATE( v1c1(nnr,ns), v2c1(nnr,ns), v2c_ud1(nnr) )
           v2c_ud1 = 0.d0
         ELSE
           ALLOCATE( grh(nnr,3,ns) )
           ALLOCATE( dvxcrr1(nnr,ns,ns), dvxcsr1(nnr,ns,ns), &
                     dvxcss1(nnr,ns,ns) )
         ENDIF
       ENDIF
    ELSEIF ( MGGA ) THEN
       ALLOCATE( v1x1(nnr,ns), v2x1(nnr,ns), v3x1(nnr,ns) )
       ALLOCATE( v1c1(nnr,ns), v2c1(nnr,ns), v3c1(nnr,ns) )
       ALLOCATE( v2cm1(np,nnr,ns) )
    ENDIF
    !
    ! ... dft2 output / benchmark data arrays
    !
    IF ( TRIM(test)=='EXECUTE' ) THEN
      IF (.NOT. xc_derivative) ALLOCATE( ex2(nnrbt), ec2(nnrbt) )
      !
      IF ( LDA .OR. GGA ) THEN
         IF ( LDA ) THEN
           IF ( .NOT. xc_derivative ) ALLOCATE( vx2(nnrbt,ns), vc2(nnrbt,ns) )
           IF ( xc_derivative ) ALLOCATE( dmuxc2(nnrbt,ns,ns) )
         ENDIF
         !
         IF ( GGA ) THEN
           IF ( .NOT. xc_derivative ) THEN
             ALLOCATE( exg2(nnrbt), ecg2(nnrbt) )
             ALLOCATE( v1x2(nnrbt,ns), v2x2(nnrbt,ns) )
             ALLOCATE( v1c2(nnrbt,ns), v2c2(nnrbt,ns), v2c_ud2(nnrbt) )
             v2c_ud1 = 0.d0
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
    ENDIF
    !
    !==========================================================================
    ! READ BENCHMARK DATA FROM FILE
    !==========================================================================
    !
    nlen1 = LEN(TRIM(dft))
    nlen2 = LEN(TRIM(family))
    IF ( dft_init=='ALL_TERMS' ) THEN
      WRITE(xc_data(9:8+nlen1),'(a)') dft(1:nlen1)
      WRITE(xc_data(14:13+nlen2),'(a)') family(1:nlen2)
      IF (is==1) WRITE(xc_data(18:30),'(a)') 'UNP'
      IF (is==2) WRITE(xc_data(18:30),'(a)') 'POL'
    ELSEIF ( dft_init=='ALL_SHORT' ) THEN
      WRITE(xc_data(9:8+nlen1),'(a)') dft(1:nlen1)
      IF (is==1) WRITE(xc_data(9+nlen1:),'(a)') 'UNP'
      IF (is==2) WRITE(xc_data(9+nlen1:),'(a)') 'POL'
    ELSEIF ( dft_init=='ALL_LIBXC' .OR. dft_init(1:4)/='ALL_' ) THEN
      xc_data="XC_DATA_______________________________________"
      WRITE(xc_data(9:8+nlen1),'(a)') dft(1:nlen1)
      WRITE(xc_data(42:41+nlen2),'(a)') family(1:nlen2)
      IF (is==1) WRITE(xc_data(42+nlen2:),'(a)') 'UNP'
      IF (is==2) WRITE(xc_data(42+nlen2:),'(a)') 'POL'
    ENDIF
    !
    ! ... read data set from xml file
    !
    not_found_skip = .FALSE.
    IF (test=='EXECUTE' .AND. mype==root) THEN
      CALL xmlr_opentag( TRIM(xc_data), tag_err )
      IF (tag_err==0) THEN
        CALL xmlr_readtag( "time_tot", time_tot2 )
        IF (.NOT. xc_derivative) THEN
          IF ( exc_term ) THEN
            CALL xmlr_readtag( "EX_AVER", ex_aver(:) )
            CALL xmlr_readtag( "EX", ex2(:) )
            IF ( family=='LDA' ) THEN
              CALL xmlr_readtag( "VX_AVER", vx_aver(:,:) )
              CALL xmlr_readtag( "VX", vx2(:,:) )
            ELSE
              CALL xmlr_readtag( "V1X_AVER", v1x_aver(:,:) )
              CALL xmlr_readtag( "V2X_AVER", v2x_aver(:,:) )
              CALL xmlr_readtag( "V1X", v1x2(:,:) )
              CALL xmlr_readtag( "V2X", v2x2(:,:) )
              IF ( family=='MGGA' ) THEN
                CALL xmlr_readtag( "V3X_AVER", v3x_aver(:,:) )
                CALL xmlr_readtag( "V3X", v3x2(:,:) )
              ENDIF
            ENDIF
          ENDIF
          IF ( cor_term ) THEN
            CALL xmlr_readtag( "EC_AVER", ec_aver(:) )
            CALL xmlr_readtag( "EC", ec2(:) )
            IF ( family=='LDA' ) THEN
              CALL xmlr_readtag( "VC_AVER", vc_aver(:,:) )
              CALL xmlr_readtag( "VC", vc2(:,:) )
            ELSE
              CALL xmlr_readtag( "V1C_AVER", v1c_aver(:,:) )
              CALL xmlr_readtag( "V2C_AVER", v2c_aver(:,:) )
              IF ( family=='GGA' .AND. is==2 ) CALL xmlr_readtag( "V2Cud_AVER", v2c_ud1_aver(:) )
              IF ( family=='MGGA' ) CALL xmlr_readtag( "V3C_AVER", v3c_aver(:,:) )
              CALL xmlr_readtag( "V1C", v1c2(:,:) )
              CALL xmlr_readtag( "V2C", v2c2(:,:) )
              IF ( family=='GGA' .AND. is==2 ) CALL xmlr_readtag( "V2Cud", v2c_ud2(:) )
              IF ( family=='MGGA' ) CALL xmlr_readtag( "V3C", v3c2(:,:) )
            ENDIF
          ENDIF
        ELSE !xc_derivative
          IF (family=='LDA') THEN 
            CALL xmlr_readtag( "dV_AVER", dv_aver(:) )
            CALL xmlr_readtag( "dV", dmuxc2(:,:,:) )
          ELSE
            CALL xmlr_readtag( "dVrr_AVER", dvrr_aver(:,:) )
            CALL xmlr_readtag( "dVsr_AVER", dvsr_aver(:,:) )
            CALL xmlr_readtag( "dVss_AVER", dvss_aver(:,:) )
            CALL xmlr_readtag( "dVXCrr", dvxcrr2(:,:,:) )
            CALL xmlr_readtag( "dVXCsr", dvxcsr2(:,:,:) )
            CALL xmlr_readtag( "dVXCss", dvxcss2(:,:,:) )
          ENDIF
        ENDIF
        CALL xmlr_closetag()
      ELSE
        CALL print_test_status( skipped2 )
        not_found_skip = .TRUE.
      ENDIF
    ENDIF
    !
#if defined(__MPI)
    CALL MPI_BCAST( not_found_skip, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr )
#endif
    IF (not_found_skip) GOTO 10
    !
    !==========================================================================
    ! BUILD ARBITRARY LARGE GRID FOR AVERAGE CALCULATIONS
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
       IF ( MGGA ) THEN
          !tau(ii,1) = fact*ABS(rho(ii,1)*ns)**(5._DP/3._DP)/ns
          tau(ii,1) = SQRT( ABS(SUM(grho(:,ii,1)**2/(rho(ii,1)*3._DP*SIN(DBLE(iip))))) )
       ENDIF
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
          IF ( MGGA ) THEN
            !tau(ii,2) = fact*ABS(rho(ii,2)*ns)**(5._DP/3._DP)/ns
            tau(ii,2) = SQRT( ABS(SUM(grho(:,ii,2)**2/(rho(ii,2)*3._DP*SIN(DBLE(iip))))) )
          ENDIF  
          !  
       ENDIF
       !
    ENDDO
    !
    !==========================================================================
    ! BUILD A SHORT INPUT GRID FOR DIRECT COMPARISON
    !==========================================================================
    !
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
    IF (mype==root) THEN
      IF (ns == 1) THEN
        rho(1,1) = 0.6_DP ; rho_tz(1,1) = 0.6_DP
        IF (family/='LDA')  grho(:,1,1) = (/ 0.1_DP, 0.2_DP, 0.3_DP /)
        IF (family=='MGGA') tau(1,1) = 0.1_DP
      ELSE
        DO i = 1, nnr_b
          IF (MOD(i,2)==0) rho(i,:) = (/ 0.6_DP, 0.1_DP /)
          IF (MOD(i,2)/=0) rho(i,:) = (/ 0.1_DP, 0.6_DP /)
          IF (MOD(i,2)==0) rho_tz(i,:) = (/ 0.7_DP, 0.5_DP /)
          IF (MOD(i,2)/=0) rho_tz(i,:) = (/ 0.7_DP, -0.5_DP /)
          IF ( family/='LDA') THEN
            IF (i<=2) THEN
              grho(:,i,1) = (/ 0.1_DP, 0.2_DP, 0.3_DP /)
              grho(:,i,2) = (/ 0.4_DP, 0.3_DP, 0.2_DP /)
            ELSE
              grho(:,i,1) = (/ 0.4_DP, 0.3_DP, 0.2_DP /)
              grho(:,i,2) = (/ 0.1_DP, 0.2_DP, 0.3_DP /)
            ENDIF
          ENDIF
        ENDDO
        !
        IF (family=='MGGA') THEN
          rho(5,:) = (/ 0.1_DP, 0.6_DP /)
          grho(:,5,1) = (/ 0.4_DP, 0.3_DP, 0.2_DP /)
          grho(:,5,2) = (/ 0.1_DP, 0.2_DP, 0.3_DP /)
          DO i = 1, 4
            tau(i,:) = (/ 0.1_DP, 0.2_DP/)
          ENDDO
          tau(5,:) = (/ 0.2_DP, 0.1_DP /)
        ENDIF
      ENDIF
    ENDIF
    !
    ! --- THRESHOLD POINTS ---
    !
    IF (mype==root) THEN
      rho(nnr_b+1,1) = thresh_lda/3.0_DP
      IF (.NOT. POLARIZED) rho_tz(nnr_b+1,1) = rho(nnr_b+1,1)
      IF ( POLARIZED ) THEN
         rho(nnr_b+1,2) = rho(nnr_b,1)
         rho_tz(nnr_b+1,1) = rho(nnr_b+1,1) + rho(nnr_b+1,2)
         rho_tz(nnr_b+1,2) = rho(nnr_b+1,1) - rho(nnr_b+1,2)
      ENDIF
      !
      IF ( GGA .OR. MGGA ) THEN
        grho(:,nnr_b+1,1) = grho(:,nnr_b,1)
        !
        rho(nnr_b+2,:) = rho(nnr_b,:)
        IF (.NOT. POLARIZED) rho_tz(nnr_b+2,1) = rho(nnr_b+2,1)
        !
        grho(:,nnr_b+2,1) = thresh_gga/10._DP
        IF ( POLARIZED ) THEN
          rho_tz(nnr_b+2,1) = rho(nnr_b+2,1) + rho(nnr_b+2,2)
          rho_tz(nnr_b+2,2) = rho(nnr_b+2,1) - rho(nnr_b+2,2)
          grho(:,nnr_b+2,2) = grho(:,nnr_b,2)
        ENDIF
        !
        IF ( MGGA ) THEN
          tau(nnr_b+1,:) = tau(nnr_b,:)
          tau(nnr_b+2,:) = tau(nnr_b,:)
          rho(nnr_b+3,:) = rho(nnr_b,:)
          grho(:,nnr_b+3,:) = grho(:,nnr_b,:)
          tau(nnr_b+3,1) = thresh_mgga/10._DP
          IF ( POLARIZED ) tau(nnr_b+3,2) = tau(nnr_b,2)
        ENDIF
        !
      ENDIF
    ENDIF
    !
    !==========================================================================
    ! CALCULATION OF ENERGIES AND POTENTIAL ARRAYS
    !==========================================================================
    !
    ! ... xclib calls for DFT
    !
    IF ( LDA ) THEN
       !
       IF ((iexch1==8 .AND. .NOT.is_libxc(1)) .OR. (icorr1==10 .AND. &
            .NOT. is_libxc(2))) CALL xclib_set_finite_size_volume( volume )
       !
       IF (mype==root) time(1) = get_wall_time()
       IF (.NOT. xc_derivative ) CALL xc( nnr, ns, ns, rho_tz, ex1, ec1, vx1, vc1 )
       IF (      xc_derivative ) CALL dmxc( nnr, ns, rho, dmuxc1 )
       IF (mype==root) time(2) = get_wall_time()
       !
    ENDIF
    !
    IF ( GGA ) THEN
      !
      IF ( .NOT. xc_derivative ) THEN
        !
        IF ( .NOT. LDA ) THEN
          ex1 = 0.d0  ;  ec1 = 0.d0
        ENDIF
        !
        IF ( dft(1:3)=='BEE' ) CALL beefsetmode(-1) !** beeforder can be manually 
        !                                           !   changed here for other checks
        !
        IF (mype==root) time(3) = get_wall_time()
        CALL xc_gcx( nnr, ns, rho, grho, exg1, ecg1, v1x1, v2x1, v1c1, &
                     v2c1, v2c_ud1 )
        IF (mype==root) time(4) = get_wall_time()
        !
        ex1 = ex1*rho_tz(:,1) + exg1
        ec1 = ec1*rho_tz(:,1) + ecg1
        IF ( LDA ) THEN
          v1x1 = v1x1 + vx1
          v1c1 = v1c1 + vc1
        ENDIF
        !
      ELSE
        !
        DO ii = 1, nnr
          grh(ii,1:3,1) = grho(1:3,ii,1)
          IF (ns==2) grh(ii,1:3,2) = grho(1:3,ii,2)
        ENDDO
        !
        IF (mype==root) time(3) = get_wall_time()
        CALL dgcxc( nnr, ns, rho, grh, dvxcrr1, dvxcsr1, dvxcss1 )
        IF (mype==root) time(4) = get_wall_time()
        !
        IF ( LDA ) dvxcrr1 = dvxcrr1 + dmuxc1
        !
      ENDIF
      !
    ENDIF
    !
    IF ( MGGA ) THEN
       IF (mype==root) time(5) = get_wall_time()
       CALL xc_metagcx( nnr, ns, np, rho, grho, tau, ex1, ec1, v1x1, &
                        v2x1, v3x1, v1c1, v2cm1, v3c1 )
       IF (mype==root) time(6) = get_wall_time()
       v2c1 = v2cm1(1,:,:)
    ENDIF
    !
    !==========================================================================
    ! COMPUTE AND PRINT TEST RESULTS
    !==========================================================================
    !
    time_tot1 = (time(6)-time(5))+(time(4)-time(3))+(time(2)-time(1))
    !
    iaverout = 0
    iavernull = 0
    !
    IF ( .NOT. xc_derivative ) THEN
      IF (exc_term) CALL evxc_stats( 'Ex', ex1, ex_aver )
      IF (cor_term) CALL evxc_stats( 'Ec', ec1, ec_aver )
    ENDIF
    !
    IF ( LDA .AND. .NOT. GGA ) THEN
       IF ( .NOT. xc_derivative ) THEN
         IF (exc_term) CALL evxc_stats( 'Vx', vx1, vx_aver(1,:) )
         IF (cor_term) CALL evxc_stats( 'Vc', vc1, vc_aver(1,:) )
       ELSE
         CALL derivxc_stats( 'dmuxc', dmuxc1, dv_aver )
       ENDIF
    ELSEIF ( GGA ) THEN
       IF ( .NOT. xc_derivative ) THEN
         IF (exc_term) CALL evxc_stats( 'V1x', v1x1, v1x_aver )
         IF (exc_term) CALL evxc_stats( 'V2x', v2x1, v2x_aver )
         IF (cor_term) CALL evxc_stats( 'V1c', v1c1, v1c_aver )
         IF (cor_term) CALL evxc_stats( 'V2c', v2c1, v2c_aver ) 
       ELSE
         CALL derivxc_stats( 'dvxcrr', dvxcrr1, dvrr_aver )
         CALL derivxc_stats( 'dvxcsr', dvxcsr1, dvsr_aver )
         CALL derivxc_stats( 'dvxcss', dvxcss1, dvss_aver )
       ENDIF
    ELSEIF ( MGGA ) THEN 
       CALL evxc_stats( 'V1x', v1x1, v1x_aver(1,:) )
       CALL evxc_stats( 'V2x', v2x1, v2x_aver(1,:) )
       CALL evxc_stats( 'V1c', v1c1, v1c_aver(1,:) )
       CALL evxc_stats( 'V2c', v2c1, v2c_aver(1,:) )
       CALL evxc_stats( 'V3x', v3x1, v3x_aver(1,:) )
       CALL evxc_stats( 'V3c', v3c1, v3c_aver(1,:) )
    ENDIF
    !
    ! ... store data set in xml file
    !
    IF (TRIM(test)=='GENERATE' .AND. mype==root) THEN
      CALL xmlw_opentag( TRIM(xc_data) )
      CALL xmlw_writetag( "time_tot", time_tot1 )
      IF (.NOT. xc_derivative) THEN
        IF ( exc_term ) THEN
          CALL xmlw_writetag( "EX_AVER", ex_aver(:) )
          CALL xmlw_writetag( "EX", ex1(1:nnrbt) )
          IF ( family=='LDA' ) THEN
            CALL xmlw_writetag( "VX_AVER", vx_aver(:,:) )
            CALL xmlw_writetag( "VX", vx1(1:nnrbt,:) )
          ELSE
            CALL xmlw_writetag( "V1X_AVER", v1x_aver(:,:) )
            CALL xmlw_writetag( "V2X_AVER", v2x_aver(:,:) )
            CALL xmlw_writetag( "V1X", v1x1(1:nnrbt,:) )
            CALL xmlw_writetag( "V2X", v2x1(1:nnrbt,:) )
            IF ( family=='MGGA' ) THEN
              CALL xmlw_writetag( "V3X_AVER", v3x_aver(:,:) )
              CALL xmlw_writetag( "V3X", v3x1(1:nnrbt,:) )
            ENDIF
          ENDIF
        ENDIF
        IF ( cor_term ) THEN
          CALL xmlw_writetag( "EC_AVER", ec_aver(:) )
          CALL xmlw_writetag( "EC", ec1(1:nnrbt) )
          IF ( family=='LDA' ) THEN
            CALL xmlw_writetag( "VC_AVER", vc_aver(:,:) )
            CALL xmlw_writetag( "VC", vc1(1:nnrbt,:) )
          ELSE
            CALL xmlw_writetag( "V1C_AVER", v1c_aver(:,:) )
            CALL xmlw_writetag( "V2C_AVER", v2c_aver(:,:) )
            IF ( family=='GGA' .AND. is==2 ) CALL xmlw_writetag( "V2Cud_AVER", v2c_ud1_aver(:) )
            IF ( family=='MGGA' ) CALL xmlw_writetag( "V3C_AVER", v3c_aver(:,:) )
            CALL xmlw_writetag( "V1C", v1c1(1:nnrbt,:) )
            CALL xmlw_writetag( "V2C", v2c1(1:nnrbt,:) )
            IF ( family=='GGA' .AND. is==2 ) CALL xmlw_writetag( "V2Cud", v2c_ud1(1:nnrbt) )
            IF ( family=='MGGA' ) CALL xmlw_writetag( "V3C", v3c1(1:nnrbt,:) )
          ENDIF
        ENDIF
      ELSE !xc_derivative
        IF (family=='LDA') THEN
          CALL xmlw_writetag( "dV_AVER", dv_aver(:) )
          CALL xmlw_writetag( "dV", dmuxc1(1:nnrbt,:,:) )
        ELSE
          CALL xmlw_writetag( "dVrr_AVER", dvrr_aver(:,:) )
          CALL xmlw_writetag( "dVsr_AVER", dvsr_aver(:,:) )
          CALL xmlw_writetag( "dVss_AVER", dvss_aver(:,:) )
          CALL xmlw_writetag( "dVXCrr", dvxcrr1(1:nnrbt,:,:) )
          CALL xmlw_writetag( "dVXCsr", dvxcsr1(1:nnrbt,:,:) )
          CALL xmlw_writetag( "dVXCss", dvxcss1(1:nnrbt,:,:) )
        ENDIF
      ENDIF
      CALL xmlw_closetag()
      CALL print_test_status( stored )
      GOTO 10
    ENDIF
    !
    !
    IF (mype == root) THEN
       !
       IF ( LDA .AND. .NOT. GGA ) THEN
         !
         iout = 0
         !
         DO ii = 1, nnrbt
            !
            IF ( .NOT. xc_derivative ) THEN
              ex_is_out = exc_term .AND. is_it_out( diff_thr_e_lda, 1, ex1(ii:ii), ex2(ii:ii) )
              ec_is_out = cor_term .AND. is_it_out( diff_thr_e_lda, 1, ec1(ii:ii), ec2(ii:ii) )
              vx_is_out = exc_term .AND. is_it_out( diff_thr_v_lda,ns, vx1(ii,:), vx2(ii,:) )
              vc_is_out = cor_term .AND. is_it_out( diff_thr_v_lda,ns, vc1(ii,:), vc2(ii,:) )
              something_out = ANY((/ex_is_out, ec_is_out, vx_is_out, vc_is_out/))
            ELSE
              dmuxc_is_out = is_dit_out( diff_thr_dmuxc, dmuxc1(ii,:,:), &
                                         dmuxc2(ii,:,:) )
            ENDIF
            !
            IF ( (.NOT.xc_derivative.AND.something_out) .OR. (xc_derivative.AND.dmuxc_is_out) ) THEN
              !
              iout = iout + 1
              !
              WRITE(stdout,*) " "
              IF ( ii > nnr_b ) THEN
                WRITE(stdout,*) "--- threshold points ---"
                WRITE(stdout,*) " "
              ENDIF
              !
              rhoi(1:ns) = rho(ii,1:ns)
              WRITE(stdout,909) nrpe+ii, nnr_b
              !
              IF ( .NOT. POLARIZED ) WRITE(stdout, 401 ) rhoi(1)
              IF (       POLARIZED ) WRITE(stdout, 402 ) rhoi(1), rhoi(2)
              WRITE(stdout,*) " "
              !
              IF (.NOT. xc_derivative) THEN
                IF ( exc_term .AND. ex_is_out ) CALL print_diff( 'Ex', ex1(ii:ii), ex2(ii:ii) )
                IF ( cor_term .AND. ec_is_out ) CALL print_diff( 'Ec', ec1(ii:ii), ec2(ii:ii) )
                IF ( exc_term .AND. vx_is_out ) CALL print_diff( 'Vx', vx1(ii,:),  vx2(ii,:)  )
                IF ( cor_term .AND. vc_is_out ) CALL print_diff( 'Vc', vc1(ii,:),  vc2(ii,:)  )
              ELSE
                CALL print_diff2( 'dmuxc', dmuxc1(ii,:,:), dmuxc2(ii,:,:) )
              ENDIF !df_ok
              !
            ENDIF
            !
         ENDDO
         !
      ELSEIF ( GGA ) THEN
         !
         iout = 0
         !
         DO ii = 1, nnrbt
           !
           IF ( .NOT. xc_derivative ) THEN
             ex_is_out = exc_term .AND. is_it_out( diff_thr_e_gga, 1, ex1(ii:ii), ex2(ii:ii) )
             ec_is_out = cor_term .AND. is_it_out( diff_thr_e_gga, 1, ec1(ii:ii), ec2(ii:ii) )
             v1x_is_out= exc_term .AND. is_it_out( diff_thr_vgga, ns, v1x1(ii,:), v1x2(ii,:) )
             v2x_is_out= exc_term .AND. is_it_out( diff_thr_vgga, ns, v2x1(ii,:), v2x2(ii,:) )
             v1c_is_out= cor_term .AND. is_it_out( diff_thr_vgga, ns, v1c1(ii,:), v1c2(ii,:) )
             IF ( .NOT. POLARIZED ) THEN
               v2c_is_out= cor_term .AND. is_it_out( diff_thr_vgga, ns, v2c1(ii,:), v2c2(ii,:) )
             ELSE
               v2c_is_out= cor_term .AND. is_it_out( diff_thr_vgga, ns, v2c1(ii,:), v2c2(ii,:), &
                                                     v2c_ud1(ii), v2c_ud2(ii) )
             ENDIF
             something_out=ANY((/ex_is_out, ec_is_out, v1x_is_out, v2x_is_out, &
                                 v1c_is_out, v2c_is_out/) )
           ELSE
             dvxcrr_is_out = is_dit_out(diff_thr_dv,dvxcrr1(ii,:,:),dvxcrr2(ii,:,:))
             dvxcsr_is_out = is_dit_out(diff_thr_dv,dvxcsr1(ii,:,:),dvxcsr2(ii,:,:))
             dvxcss_is_out = is_dit_out(diff_thr_dv,dvxcss1(ii,:,:),dvxcss2(ii,:,:))
             dvgga_is_out = ANY((/dvxcrr_is_out, dvxcsr_is_out, dvxcss_is_out/))
           ENDIF
           !
           IF ( (.NOT.xc_derivative.AND.something_out) .OR. (xc_derivative.AND.dvgga_is_out) ) THEN
             !
             iout = iout + 1
             !
             IF ( ii>nnr_b ) THEN
                WRITE(stdout,*) "--- threshold points ---"
                WRITE(stdout,*) " "
             ENDIF 
             !
             WRITE(stdout,*) " "
             rhoi(1:ns)=rho(ii,1:ns) ; grhoi(:,1:ns) = grho(:,ii,1:ns)
             WRITE(stdout,909) nrpe+ii, nnrbt
             !
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
             ! 
             IF (.NOT. xc_derivative) THEN 
               ! 
               IF (exc_term .AND. ex_is_out)  CALL print_diff( 'Ex', ex1(ii:ii), ex2(ii:ii) )
               IF (cor_term .AND. ec_is_out)  CALL print_diff( 'Ec', ec1(ii:ii), ec2(ii:ii) )
               IF (exc_term .AND. v1x_is_out) CALL print_diff( 'V1x',v1x1(ii,:), v1x2(ii,:) )
               IF (exc_term .AND. v2x_is_out) CALL print_diff( 'V2x',v2x1(ii,:), v2x2(ii,:) )
               IF (cor_term .AND. v1c_is_out) CALL print_diff( 'V1c',v1c1(ii,:), v1c2(ii,:) )
               IF ( .NOT. POLARIZED ) THEN
                 IF (cor_term .AND. v2c_is_out) CALL print_diff( 'V2c',v2c1(ii,:), v2c2(ii,:))
               ELSE
                 IF (cor_term .AND. v2c_is_out) CALL print_diff( 'V2c',v2c1(ii,:), v2c2(ii,:),&
                                                                 v2c_ud1(ii), v2c_ud2(ii) )
               ENDIF
             ELSE
               ! 
               !WRITE(stdout,*) " " 
               ! 
               IF (dvxcrr_is_out) CALL print_diff2( 'dvxcrr',dvxcrr1(ii,:,:), &
                                                             dvxcrr2(ii,:,:) )
               IF (dvxcsr_is_out) CALL print_diff2( 'dvxcsr',dvxcsr1(ii,:,:), &
                                                             dvxcsr2(ii,:,:) )
               IF (dvxcss_is_out) CALL print_diff2( 'dvxcss',dvxcss1(ii,:,:), &
                                                             dvxcss2(ii,:,:) )
             ENDIF
             !
           ENDIF  
           ! 
         ENDDO 
         ! 
         !
      ELSEIF ( MGGA ) THEN 
         !
         ! ... calculate values over a few benchmark points (nnr_b) 
         ! 
         iout = 0 
         !
         DO ii = 1, nnrbt
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
           !
           IF ( something_out ) THEN 
              ! 
              iout = iout + 1 
              
              IF ( ii>nnr_b ) THEN
                WRITE(stdout,*) "--- threshold points ---"
                WRITE(stdout,*) " "
              ENDIF
              !
              WRITE(stdout,*) " "   
              WRITE(stdout,909) ii, nnr_b  
              rhoi(1:ns) = rho(ii,1:ns)  
              grhoi(:,1:ns) = grho(:,ii,1:ns)  
              taui(1:ns) = tau(ii,1:ns)      
              !   
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
           !
         ENDDO
      ENDIF
      !
      IF ( TRIM(test)=='EXECUTE' ) THEN
        IF (iout+iaverout/=0) CALL print_test_status( failed )
        IF (iout+iaverout==0) THEN
          IF (iavernull/=naver) CALL print_test_status( passed )
          IF (iavernull==naver) CALL print_test_status( passed0 )
        ENDIF
      ENDIF  
      !
    ENDIF
    !
    !
    !==========================================================================
    ! FINALIZE
    !==========================================================================
    !
    10 CONTINUE
    !
    IF (xclib_dft_is_libxc('ANY')) CALL xclib_finalize_libxc()
    !
    DEALLOCATE( rho, rho_tz )
    !
    IF ( GGA .OR. MGGA ) DEALLOCATE( grho )
    IF ( MGGA ) DEALLOCATE( tau )
    !
    ! ... set1 output arrays
    !
    IF (.NOT. xc_derivative) DEALLOCATE( ex1, ec1 )
    !
    IF ( LDA .OR. GGA ) THEN
       IF ( LDA ) THEN
         IF ( .NOT. xc_derivative ) DEALLOCATE( vx1, vc1 )
         IF ( xc_derivative ) DEALLOCATE( dmuxc1 )
       ENDIF
       !
       IF ( GGA ) THEN
         IF ( .NOT. xc_derivative ) THEN
           DEALLOCATE( exg1, ecg1 )
           DEALLOCATE( v1x1, v2x1 )
           DEALLOCATE( v1c1, v2c1, v2c_ud1 )
         ELSE
           DEALLOCATE( grh )
           DEALLOCATE( dvxcrr1, dvxcsr1, dvxcss1 )
         ENDIF
       ENDIF
    ELSEIF ( MGGA ) THEN
       DEALLOCATE( v1x1, v2x1, v3x1 )
       DEALLOCATE( v1c1, v2c1, v3c1 )
       DEALLOCATE( v2cm1 )
    ENDIF
    !
    ! ... set2 output / benchmark data arrays
    !
    IF ( TRIM(test)=='EXECUTE' ) THEN
      IF (.NOT. xc_derivative) DEALLOCATE( ex2, ec2 )
      !
      IF ( LDA .OR. GGA ) THEN
         IF ( LDA ) THEN
           IF ( .NOT. xc_derivative ) DEALLOCATE( vx2, vc2 )
           IF ( xc_derivative ) DEALLOCATE( dmuxc2 )
         ENDIF
         !
         IF ( GGA ) THEN
           IF ( .NOT. xc_derivative ) THEN
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
    ENDIF
    !
   ENDDO 
  ENDDO ! end of main loop over dfts
  !
  IF (mype==root) THEN
    IF (TRIM(test)=='EXECUTE' ) CALL xmlr_closetag()
    IF (TRIM(test)=='GENERATE') CALL xmlw_closetag()
    !
    CALL xml_closefile( )
  ENDIF
  !
  401 FORMAT('rho: ',F17.14)
  402 FORMAT('rho(up,down): ',F17.14,4x,F17.14)
  !
  501 FORMAT('grho2: ',F17.14)
  503 FORMAT('grho2(uu,ud,dd): ',F17.14,4x,F17.14,4x,F17.14)
  !
  601 FORMAT('tau: ',F17.14)
  602 FORMAT('tau(up,down): ',F17.14,4x,F17.14)
  !
  909 FORMAT('grid-point: ',I5,' of ',I5)
  !   
  DEALLOCATE( proc_name )
  DEALLOCATE( node_name )
  DEALLOCATE( proc2node )
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
 !--------------------------------------------------------------------
 SUBROUTINE print_aver( what, vaver, averref, what2 )
  !------------------------------------------------------------------
  !! Prints average, max and min differences between XC arrays
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(IN) :: what
  CHARACTER(len=*), INTENT(IN), OPTIONAL :: what2
  REAL(DP), INTENT(IN) :: vaver(2)
  REAL(DP), OPTIONAL :: averref
  !
  IF (ABS(vaver(1)-averref)>aver_thresh) THEN
    WRITE(stdout,*) " "
    IF (PRESENT(what2)) THEN
      WRITE(stdout,*) " ", TRIM(what),TRIM(what2)
    ELSE
      WRITE(stdout,*) " ", TRIM(what)
    ENDIF
    WRITE(stdout,*) "AVR test: ", vaver(1)
    WRITE(stdout,*) "AVR ref : ", averref
    WRITE(stdout,*) "diff    : ", vaver(1)-averref
    iaverout=iaverout+1
  ENDIF
  !
  IF (vaver(1)==0.d0 .AND. averref==0.d0) iavernull=iavernull+1
  !
 END SUBROUTINE print_aver
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
  101 FORMAT('test: ',F17.14)
  102 FORMAT('test: ',F17.14,4x,F17.14)
  103 FORMAT('test: ',F17.14,4x,F17.14,4x,F17.14)
  201 FORMAT('ref:  ',F17.14)
  202 FORMAT('ref:  ',F17.14,4x,F17.14)
  203 FORMAT('ref:  ',F17.14,4x,F17.14,4x,F17.14)
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
  101 FORMAT('test: ',F17.14)
  103 FORMAT('test: ',F17.14,4x,F17.14,4x,F17.14)
  201 FORMAT('ref:  ',F17.14)
  203 FORMAT('ref:  ',F17.14,4x,F17.14,4x,F17.14)
  301 FORMAT('diff: ',4x,F17.14)
  303 FORMAT('diff: ',4x,F17.14,4x,F17.14,4x,F17.14)
  !
 END SUBROUTINE print_diff2
 !
 !-------------------------------------------------------------------------
 SUBROUTINE evxc_stats( what, xc_1, aver )
  !------------------------------------------------------------------------
  !! If test=EXECUTE calculates difference between total energy
  !! and potential calculated over npoints k-points and values taken from
  !! benchmark data file.  
  !! If test=GENERATE calculates the total energy and potential 
  !! over npoints k-points.
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: what
  REAL(DP), INTENT(IN) :: xc_1(nnr,ns)
  REAL(DP), INTENT(INOUT) :: aver(2)
  !
  REAL(DP) :: xc_aver(2,2)
  REAL(DP) :: thr, aver_snd(1), aver_rec(1)
  INTEGER :: ierr2
  !
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
  !IF (mype==root .AND. TRIM(test)=='EXECUTE') THEN
  !  WRITE(stdout,*) " "
  !  IF ( POLARIZED .AND. what(1:1)/='E' ) WRITE(stdout,*) TRIM(what)
  !ENDIF
  !
  xc_aver=0._DP
  !
  xc_aver(1,1) = SUM(xc_1(1:nnr,1))/DBLE(npoints)
  !
#if defined(__MPI)
  aver_snd = xc_aver(1,1)
  CALL MPI_REDUCE( aver_snd, aver_rec, 1, MPI_DOUBLE_PRECISION, MPI_SUM,0, &
                   comm, ierr2 )
  xc_aver(1:1,1) = aver_rec
#endif
  !
  IF ( .NOT. POLARIZED .OR. what(1:1)=='E' ) THEN
    IF (mype==root .AND. TRIM(test)=='EXECUTE') CALL print_aver( what, xc_aver(:,1), aver(1) )
  ELSE
    xc_aver(1,2) = SUM(xc_1(1:nnr,2))/DBLE(npoints)
    !
#if defined(__MPI)
    aver_snd = xc_aver(1,2)
    CALL MPI_REDUCE( aver_snd, aver_rec, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                     0, comm, ierr2 )
    xc_aver(1:1,2) = aver_rec
#endif
    !
    IF (TRIM(what)=='V2c' .AND. GGA ) THEN
      v2c_ud1_aver(1) = SUM(v2c_ud1(1:nnr))/DBLE(npoints)
      ! 
#if defined(__MPI)
      aver_sndu = v2c_ud1_aver(1)
      CALL MPI_REDUCE( aver_sndu, aver_recu, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, 0, comm, ierrm )
      v2c_ud1_aver(1) = aver_recu
#endif
      !
      IF (mype==root .AND. TRIM(test)=='EXECUTE') CALL print_aver( what, &
                                 v2c_ud1_aver, v2c_aver(1,3), 'cross' )
    ENDIF
    !
    IF (mype==root .AND. TRIM(test)=='EXECUTE') THEN
      CALL print_aver( what, xc_aver(:,1), aver(1), ' up' )
      CALL print_aver( what, xc_aver(:,2), aver(2), ' down' )
    ENDIF 
    !
  ENDIF
  !
  IF (TRIM(test)=='GENERATE') THEN
     aver = xc_aver(1,:)
     IF (TRIM(what)=='V2c'.AND.GGA.AND.ns==2) v2c_aver(1,3)=v2c_ud1_aver(1)
  ENDIF
  !
  RETURN
  !
 END SUBROUTINE evxc_stats
 !
 !
 !---------------------------------------------------------------------
 SUBROUTINE derivxc_stats( what, dxc_qe, aver )
  !--------------------------------------------------------------------
  !! Same as \(\texttt{evxc_stats}\), but for derivatives of Vxc.
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: what
  REAL(DP), INTENT(IN) :: dxc_qe(nnr,ns,ns)
  REAL(DP), INTENT(INOUT) :: aver(3)
  !
  REAL(DP) :: dxc_aver(2,np)
  REAL(DP) :: thr, aver_snd(1), aver_rec(1)
  INTEGER :: ierr2
  !
  IF (LDA .AND. .NOT.GGA) thr = diff_thr_dmuxc
  IF ( GGA ) thr = diff_thr_dv
  !
  dxc_aver=0._DP
  !
  dxc_aver(1,1) = SUM(dxc_qe(1:nnr,1,1))/DBLE(npoints)
  !
#if defined(__MPI)
  aver_snd = dxc_aver(1,1)
  CALL MPI_REDUCE( aver_snd, aver_rec, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                   0, comm, ierr2 )
  dxc_aver(1:1,1) = aver_rec
#endif
  !
  IF ( .NOT. POLARIZED ) THEN
    IF (mype==root .AND. TRIM(test)=='EXECUTE') CALL print_aver( what, &
                                          dxc_aver(1:nnr_b,1), aver(1) )
  ELSE
    dxc_aver(1,2) = SUM(dxc_qe(1:nnr,1,2))/DBLE(npoints)
    !
#if defined(__MPI)
    aver_snd = dxc_aver(1,2)
    CALL MPI_REDUCE( aver_snd, aver_rec, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, 0, comm, ierr2 )
    dxc_aver(1:1,2) = aver_rec
#endif
    dxc_aver(1,3) = SUM(dxc_qe(1:nnr,2,2))/DBLE(npoints)
    !
#if defined(__MPI)
    aver_snd = dxc_aver(1,3)
    CALL MPI_REDUCE( aver_snd, aver_rec, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, 0, comm, ierr2 )
    dxc_aver(1:1,3) = aver_rec
#endif
    !
    IF (mype==root .AND. TRIM(test)=='EXECUTE') THEN
      CALL print_aver( what, dxc_aver(:,1), aver(1), ' up-up' )
      CALL print_aver( what, dxc_aver(:,2), aver(2), ' up-down' )
      CALL print_aver( what, dxc_aver(:,3), aver(3), ' down-down' )
    ENDIF
    !
  ENDIF
  !
  IF (TRIM(test)=='GENERATE') aver(1:np) = dxc_aver(1,:)
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
  INTEGER :: j
  !
  is_it_out = ANY(ABS(x_dft1(1:dm)-x_dft2(1:dm)) > diff_thr)
  !
  IF (PRESENT(x_ud_1)) THEN
    is_it_out_ud =  ABS(x_ud_1-x_ud_2) > diff_thr
    is_it_out = ANY( (/ is_it_out, is_it_out_ud /) )
  ENDIF
  !
  ! ... to flush out NaN if any
  DO j = 1, dm
    IF ( .NOT.x_dft1(j)>0.111_DP .AND. .NOT.x_dft1(j)<=0.111_DP .OR. &
        (.NOT.x_dft2(j)>0.111_DP .AND. .NOT.x_dft2(j)<=0.111_DP) ) is_it_out = .TRUE.
  ENDDO
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
  INTEGER :: j
  !
  dxc_diff(1) = ABS(dx_dft1(1,1)-dx_dft2(1,1))
  !
  IF ( POLARIZED ) THEN
    dxc_diff(2) = ABS(dx_dft1(2,1)-dx_dft2(2,1))
    dxc_diff(3) = ABS(dx_dft1(2,2)-dx_dft2(2,2))
  ENDIF
  !
  is_dit_out = ANY(dxc_diff(1:np) > diff_thr)
  !
  ! ... to flush out NaN if any
  DO j = 1, np
    IF ( .NOT.dxc_diff(j)>0.111_DP .AND. .NOT.dxc_diff(j)<=0.111_DP ) is_dit_out = .TRUE.
  ENDDO
  !
 END FUNCTION
 !
 !--------------------------------------------------------------------------
 SUBROUTINE print_test_status( status )
  !------------------------------------------------------------------------
  !! Print test status on screen.
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: status
  CHARACTER(LEN=30)  :: dft_out
  CHARACTER(LEN=100) :: test_output_gen
  CHARACTER(LEN=115) :: test_output_exe
  INTEGER :: j, id_term
  !
  dft_out = dft
  IF (dft_init=='ALL_LIBXC') dft_out = dft_lxc
  !
  IF (TRIM(test)=='GENERATE') THEN
    test_output_gen = ''
    IF (dft_init=='ALL_TERMS') THEN
      DO j = 1, 6
        IF (id_vec(j)/=0) id_term = id_vec(j)
      ENDDO
      IF (is==is_min) WRITE(test_output_gen(1:3), '(i3)') id_term
      WRITE(test_output_gen(5:8),  '(a)') TRIM(family)
      WRITE(test_output_gen(10:11),'(a)') TRIM(xc_kind)
    ELSE
      IF (is==is_min) WRITE(test_output_gen(1:3), '(i3)') id
    ENDIF
    IF (is==1) WRITE(test_output_gen(13:17), '(a)') 'UNPOL'
    IF (is==2) WRITE(test_output_gen(13:15), '(a)') 'POL'
    WRITE(test_output_gen(19:54), '(a)') TRIM(dft_out)
    WRITE(test_output_gen(56:),'(a)') TRIM(status)
    WRITE(stdout,*) test_output_gen
  ELSEIF (TRIM(test)=='EXECUTE') THEN
    test_output_exe = ''
    IF (dft_init=='ALL_TERMS') THEN
      DO j = 1, 6
        IF (id_vec(j)/=0) id_term = id_vec(j)
      ENDDO
      IF (is==is_min) WRITE(test_output_exe(1:3), '(i3)') id_term
      WRITE(test_output_exe(5:8), '(a)') TRIM(family)
      WRITE(test_output_exe(10:11),'(a)') TRIM(xc_kind)
    ELSE
      IF (is==is_min) WRITE(test_output_exe(1:3), '(i3)') id
    ENDIF
    IF (is==1) WRITE(test_output_exe(13:17), '(a)') 'UNPOL'
    IF (is==2) WRITE(test_output_exe(13:15), '(a)') 'POL'
    WRITE(test_output_exe(19:54), '(a)') TRIM(dft_out)
    WRITE(test_output_exe(56:96), '(a)') TRIM(status)
    IF (status==passed .AND. show_time) THEN
      WRITE(test_output_exe(80:85), '(a)') 'time:'
      WRITE(test_output_exe(86:92), '(F6.3)') time_tot1
      WRITE(test_output_exe(93:100), '(a)') 's  incr:'
      WRITE(test_output_exe(101:110), '(F8.3)') time_tot1/time_tot2*100.d0-100.d0
      WRITE(test_output_exe(111:111), '(a)') '%'
    ENDIF
    WRITE(stdout,*) test_output_exe
  ENDIF
  !
  RETURN
  !
 END SUBROUTINE  
 !
 !---------------------------------------------------------------------------
 FUNCTION get_wall_time()
    !------------------------------------------------------------------------
    !! Get wall time (copied from LAXlib example).
    REAL(DP) :: get_wall_time
#if defined(__MPI)
    get_wall_time = MPI_WTIME()
#else
    INTEGER :: cr, nc
    REAL(DP), SAVE :: t0 = -1.0
    !
    CALL system_clock(count_rate=cr)
    CALL system_clock(count=nc)
    IF ( t0 < 0.0 ) t0 = DBLE(nc)/cr
    get_wall_time = DBLE(nc)/cr - t0
#endif
  END FUNCTION get_wall_time
  ! 
END PROGRAM xclib_test
