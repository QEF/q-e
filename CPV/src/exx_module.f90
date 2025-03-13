MODULE exx_module
  !----------------------------------------------------------------------------------------------------------------
  !! Exact exchange calculation using maximally localized Wannier functions 
  !! The references for this algorithm are:  
  !! (i)  theory: X. Wu , A. Selloni, and R. Car, Phys. Rev. B 79, 085102 (2009).  
  !! (ii) implementation: H.-Y. Ko, J. Jia, B. Santra, X. Wu, R. Car, and 
  !!      R. A. DiStasio Jr., J. Chem. Theory Comput. 16, 3757 (2020).
  !
  !! Contributors:  
  !! Cornell University   - Hsin-Yu Ko, Junteng Jia, Robert A. DiStasio  
  !! Princeton University - Marcos F. Calegari Andrade, Lingzhu Kong, Zhaofeng Li, Roberto Car  
  !! Temple University    - Biswajit Santra, Xifan Wu, Charles W. Swartz
  !
  !! Code Version 1.1 (November 2020)
  !
  !----------------------------------------------------------------------------------------------------------------
  ! For the developers: 
  ! the related files are exx_cg.f90 exx_es.f90 exx_gs.f90 exx_module.f90 exx_pair.f90 exx_psi.f90 exx_vofr.f90
  ! to find related changes made in other files look for the string 'exx_wf related'
  !----------------------------------------------------------------------------------------------------------------
  !
  USE cell_base,          ONLY: h                  !cell at time t (current time step), also used in r = h s
  USE cell_base,          ONLY: ainv               !h^-1 matrix for converting between r and s coordinates via s = h^-1 r)
  USE cell_base,          ONLY: omega              !cell volume (in au^3)
  USE cell_base,          ONLY: isotropic          !if .TRUE. then "volume" option is chosen for cell_dofree
  USE cell_base,          ONLY: ibrav              !ibrav determines the symmetry of the simulation cell. See manual for acceptable values ... 
  USE cell_base,          ONLY: ref_at             !reference cell parameters in ref_alat unit
  USE command_line_options, ONLY : ndiag_          !-ndiag indicator; =0 -ndiag not specified; = N if -ndiag N specified  
  USE constants,          ONLY: pi                 !pi in double-precision
  USE constants,          ONLY: fpi                !4.0*pi in double-precision
  USE control_flags,      ONLY: lwfnscf            !
  USE control_flags,      ONLY: lwfpbe0nscf        !non selfconsitent pbe0 calculation for empty bands .. 
  USE control_flags,      ONLY: nbeg               !nbeg<0 in from_scratch calculations ...
  USE control_flags,      ONLY: thdyn              !if .TRUE. then variable cell calculation is turned on ..   
  USE cp_main_variables,  ONLY: idesc              !descriptor type
  USE electrons_base,     ONLY: nbsp               !number of electronic bands/states ...
  USE electrons_base,     ONLY: nspin              !spin unpolarized (npsin=1) vs. spin polarized (nspin=2) specification
  USE electrons_base,     ONLY: nupdwn             !number of states with up and down spin 
  USE fft_base,           ONLY: dffts              !FFT derived data type
  USE fft_base,           ONLY: dfftp              !FFT derived data type 
  USE xc_lib,             ONLY: xclib_get_exx_fraction, &
                                stop_exx, start_exx
  USE input_parameters,   ONLY: ref_alat           !alat of reference cell ..
  USE input_parameters,   ONLY: ref_cell           !.true. if reference cell parameters are in use, .false. otherwise ... 
  USE io_global,          ONLY: stdout             !print/write argument for standard output (to output file)
  USE kinds,              ONLY: DP                 !double-precision kind (selected_real_kind(14,200))
  USE mp,                 ONLY: mp_sum             !MPI collection with sum
  USE mp_global,          ONLY: nproc_image        !number of processors
  USE mp_global,          ONLY: me_image           !processor number (0,1,...,nproc_image-1)
  USE mp_global,          ONLY: me_bgrp            !processor number (0,1,...,nproc_image-1)
  USE mp_global,          ONLY: intra_image_comm   !standard MPI communicator
  USE parallel_include
  USE wannier_base,       ONLY: neigh         ! maximum number of neighbors set in the calculation to initialize many arrays  ...
  USE wannier_base,       ONLY: dis_cutoff    ! radius to obtain WF pairs 
  USE wannier_base,       ONLY: exx_ps_rcut_s ! radius of the poisson sphere for self orbital
  USE wannier_base,       ONLY: exx_me_rcut_s ! radius of the ME sphere for self orbital
  USE wannier_base,       ONLY: exx_ps_rcut_p ! radius of the poisson sphere for pair orbital
  USE wannier_base,       ONLY: exx_me_rcut_p ! radius of the ME sphere for pair orbital
  USE wannier_base,       ONLY: vnbsp 
  USE wannier_base,       ONLY: texx_cube     ! flag to toggle the cube domain
  USE fft_helper_subroutines
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! PUBLIC variables 
  !
  INTEGER, PARAMETER, PUBLIC          :: lmax=6
  !! maximum angular momentum
  INTEGER, PARAMETER, PUBLIC          :: nord1=3
  !! order of expansion for first derivatives (points on one side)
  INTEGER, PARAMETER, PUBLIC          :: nord2=3
  !! order of expansion for second derivatives (points on one side)
  !
  INTEGER, PUBLIC                     :: n_exx
  !! index of exx steps
  INTEGER, PUBLIC                     :: sc_fac
  !! scaling factor to parallelize on (nproc_image=integer multiple of nbsp)
  ! sc_fac=INT(nproc_image/nbsp); the parallelization can be generalized
  ! for any number of procesors by using an array of sc_fac(1:nbsp)
  INTEGER, PUBLIC                     :: np_in_sp_s
  !! number of grid points in the PS sphere for self orbital in Full Grid 
  INTEGER, PUBLIC                     :: np_in_sp_me_s
  !! number of grid points in the ME sphere for self orbital in Full Grid 
  INTEGER, PUBLIC                     :: np_in_sp_p
  !! number of grid points in the PS sphere for pair orbital in Full Grid
  INTEGER, PUBLIC                     :: np_in_sp_me_p
  !! number of grid points in the ME sphere for pair orbital in Full Grid
  !
  INTEGER, PUBLIC                     :: my_nbspx
  !! parallelization/distribution of orbitals over processors 
  INTEGER,  ALLOCATABLE, PUBLIC       :: my_nbsp(:)
  !! parallelization/distribution of orbitals over processors 
  INTEGER,  ALLOCATABLE, PUBLIC       :: my_nxyz(:)
  !! parallelization/distribution of orbitals over processors 
  INTEGER,  ALLOCATABLE, PUBLIC       :: index_my_nbsp(:,:)
  !! parallelization/distribution of orbitals over processors 
  INTEGER,  ALLOCATABLE, PUBLIC       :: rk_of_obtl(:)
  !! parallelization/distribution of orbitals over processors 
  INTEGER,  ALLOCATABLE, PUBLIC       :: lindex_of_obtl(:)
  !! parallelization/distribution of orbitals over processors 
  INTEGER,  ALLOCATABLE, PUBLIC       :: pair_label(:,:)
  !! the orbital label of previous pair potential/density
  INTEGER,  ALLOCATABLE, PUBLIC       :: pair_step(:,:)
  !! the last step that we use previous pair potential/density
  INTEGER,  ALLOCATABLE, PUBLIC       :: pair_status(:,:)
  !! the status of this pair for guess
  !
  ! conversion between 3D index (i,j,k) and 1D index 
  ! odthothd_in_sp(3, 1:np_in_sp_p) is for inner sphere (1st shell)
  ! odthothd_in_sp(3, np_in_sp_p+1:np_tmp_1) is for 2nd shell
  ! odthothd_in_sp(3, np_tmp_1+1:np_tmp_2) is for 3rd shell
  ! odthothd_in_sp(3, np_tmp_2:np_in_sp_me_s) is for 4th shell
  INTEGER,  ALLOCATABLE, PUBLIC       :: odtothd_in_sp(:,:)
  !! 1D to 3D index converter in local grid within largest ME sphere
  INTEGER,  ALLOCATABLE, PUBLIC       :: thdtood_in_sp(:,:,:)
  !! 3D to 1D index converter in local grid within largest ME sphere
  INTEGER,  ALLOCATABLE, PUBLIC       :: thdtood(:,:,:)
  !! 3D to 1D index converter in global grid (nr1,nr2,nr3 => nr1*nr2*nr3) 
  !
  REAL(DP), ALLOCATABLE, PUBLIC       :: rhopr(:,:)
  !
  REAL(DP), ALLOCATABLE, PUBLIC       :: xx_in_sp(:)
  !! distance between centre of the simulation box and grid points along X direction .. 
  REAL(DP), ALLOCATABLE, PUBLIC       :: yy_in_sp(:)
  !! distance between centre of the simulation box and grid points along Y direction .. 
  REAL(DP), ALLOCATABLE, PUBLIC       :: zz_in_sp(:)
  !! distance between centre of the simulation box and grid points along Z direction .. 
  !
  REAL(DP), ALLOCATABLE, PUBLIC       :: sc_xx_in_sp(:)
  !! scaled distance between centre of the simulation box and grid points along X direction .. 
  REAL(DP), ALLOCATABLE, PUBLIC       :: sc_yy_in_sp(:)
  !! scaled distance between centre of the simulation box and grid points along Y direction .. 
  REAL(DP), ALLOCATABLE, PUBLIC       :: sc_zz_in_sp(:)
  !! scaled distance between centre of the simulation box and grid points along Z direction .. 
  !
  REAL(DP), ALLOCATABLE, PUBLIC       :: selfv(:,:,:)
  !! self potential stored in Poisson sphere as guess potential ...
  REAL(DP), ALLOCATABLE, PUBLIC       :: pair_dist(:,:,:)
  !! pair distance stored in order to extrapolate guess potential ...    
  REAL(DP), ALLOCATABLE, PUBLIC       :: pairv(:,:,:,:)
  !! pair potential stored in Poisson sphere as guess potential ...
  !
  REAL(DP), ALLOCATABLE, PUBLIC       :: exx_potential(:,:)
  !! EXX potential which is passed/added to the DFT-GGA potential ... 
  !
  REAL(DP), ALLOCATABLE, PUBLIC       :: clm(:,:)
  REAL(DP), ALLOCATABLE, PUBLIC       :: coeke(:,:,:)           ! coeke(neighbor, d/di, d/dj)
  REAL(DP), ALLOCATABLE, PUBLIC       :: coe_1st_derv(:,:)      ! coe_1st_derv(neighbor, d/di)
  REAL(DP), ALLOCATABLE, PUBLIC       :: vwc(:,:)
  REAL(DP), PUBLIC                    :: h_init(3,3)
  REAL(DP), PUBLIC                    :: dexx_dh(3,3)           ! dexx/dhab for vofrho.f90
  REAL(DP), PUBLIC                    :: exxalfa                ! fraction of exx mixing (locally used in CP)
  real(dp), allocatable, public  :: psime_pair_recv(:,:,:) !! recieving buffer holding MLWF on local subdomains
  !! the first dimension is the (flattened) 1d local grid index; the second is the jth neighbor orbital index;
  !! the third is the ith local orbital index of the current MPI process
  real(dp), allocatable, public  :: psime_pair_send(:,:,:) !! sending buffer holding MLWF on local subdomains
  !! the first dimension is the (flattened) 1d local grid index; the second is the jth neighbor orbital index;
  !! the third is the ith local orbital index of the current MPI process
  !
#if defined(_OPENMP)
  INTEGER, EXTERNAL                   :: omp_get_max_threads
#endif
#if defined __CUDA
  REAL(DP), ALLOCATABLE, PUBLIC, DEVICE       :: coe_1st_derv_d(:,:)      ! coe_1st_derv(neighbor, d/di)
  REAL(DP), ALLOCATABLE, PUBLIC, DEVICE       :: coeke_d(:,:,:)           ! coeke(neighbor, d/di, d/dj)
  REAL(DP), ALLOCATABLE, PUBLIC, DEVICE       :: coemicf_d(:,:,:)         ! coefficient for preconditioner
#endif
  ! cubic domain related variables
  INTEGER, PARAMETER, PUBLIC          :: lm_mx=((lmax+1)*(lmax+2))/2  ! size of flattened 1d lm angular momentum ... 
  !----------------------------------------------------------------------------------------------------------------
  ! JJ: new variables
  !----------------------------------------------------------------------------------------------------------------
  INTEGER, PUBLIC                     :: s_me_r(6), n_s_me      ! start/ending indexes for 3-dimension in self multipole cube
  INTEGER, PUBLIC                     :: s_ps_r(6), n_s_ps      ! start/ending indexes for 3-dimension in self poisson cube
  INTEGER, PUBLIC                     :: p_me_r(6), n_p_me      ! start/ending indexes for 3-dimension in pair multipole cube
  INTEGER, PUBLIC                     :: p_ps_r(6), n_p_ps      ! start/ending indexes for 3-dimension in pair poisson cube 
  !----------------------------------------------------------------------------------------------------------------
  INTEGER, PUBLIC                     :: PScubeSL_s
  INTEGER, PUBLIC                     :: PScubeSL_p
  INTEGER, PUBLIC                     :: MEcubeSL_s
  INTEGER, PUBLIC                     :: MEcubeSL_p
  INTEGER, PUBLIC                     :: nrg(3)
  INTEGER, PUBLIC                     :: nrgr(3)
  !----------------------------------------------------------------------------------------------------------------
  REAL(DP), ALLOCATABLE, PUBLIC       :: me_cs(:,:,:,:)
  !! coordinate  of every point in ME cube
  REAL(DP), ALLOCATABLE, PUBLIC       :: me_rs(:,:,:,:)
  !! distance^n  of every point in ME cube
  REAL(DP), ALLOCATABLE, PUBLIC       :: me_ri(:,:,:,:)
  !! distance^-n of every point in ME cube
  COMPLEX(DP), ALLOCATABLE, PUBLIC    :: me_rc(:,:,:,:)
  !! exp(i*m*phi_j) of everg point in PS cube
#ifdef __CUDA
  REAL(DP), ALLOCATABLE, PUBLIC   , DEVICE :: me_cs_d(:,:,:,:)       ! coordinate  of every point in ME cube
  REAL(DP), ALLOCATABLE, PUBLIC   , DEVICE :: me_rs_d(:,:,:,:)       ! distance^n  of every point in ME cube
  REAL(DP), ALLOCATABLE, PUBLIC   , DEVICE :: me_ri_d(:,:,:,:)       ! distance^-n of every point in ME cube
  COMPLEX(DP), ALLOCATABLE, PUBLIC, DEVICE :: me_rc_d(:,:,:,:)       ! exp(i*m*phi_j) of everg point in PS cube
#endif
  REAL(DP), ALLOCATABLE, PUBLIC       :: selfrho(:,:,:)
  !! self density stored in Poisson sphere for guess potential ...
  REAL(DP), ALLOCATABLE, PUBLIC       :: pairrho(:,:,:,:)
  !! pair density stored in Poisson sphere for guess potential ...
  REAL(DP), PUBLIC                    :: fbsscale
  !! coefficient for preconditioner
  REAL(DP), ALLOCATABLE, PUBLIC       :: coemicf(:,:,:)
  !! coefficient for preconditioner
  real(dp), allocatable, public  :: rho_ps(:)
  real(dp), allocatable, public  :: pot_ps(:)
  INTEGER ::  i, iobtl, gindex_of_iobtl, irank, proc, tmp_iobtl, ndiag_n, ndiag_nx, ndiag_i
#ifdef __CUDA
  attributes(device) :: rho_ps, pot_ps
  attributes(pinned) :: psime_pair_send, psime_pair_recv
  real(dp), allocatable, device :: psime_pair_recv_d(:,:,:),psime_pair_send_d(:,:,:)
  real(dp), allocatable, device :: psi_d(:,:)!,rhops_d(:)
  real(dp), allocatable, device :: vpsil_d(:,:)
#endif
  !==========================================================================
  !
  ! PRIVATE variables 
  !
  INTEGER, PRIVATE :: nr1,nr2,nr3                           !real space grid dimensions (global first, second, and third dimensions of the 3D grid)
  INTEGER, PRIVATE :: nr1r,nr2r,nr3r                        !reduced real space grid dimensions (global first, second, and third dimensions of the 3D grid)
  INTEGER, PRIVATE :: ierr                                  !error
  REAL(DP), PRIVATE :: h_in(3,3)                            !cell parameter to initialize ... 
  !
  ! PUBLIC subroutines
  !
  PUBLIC :: exx_initialize
  PUBLIC :: exx_finalize
  PUBLIC :: getnpinsp 
  PUBLIC :: exx_setup 
  PUBLIC :: exx_setup_nscf
  PUBLIC :: fornberg
  PUBLIC :: exx_energy_cell_derivative
  !
  ! PRIVATE subroutines
  !
  PRIVATE :: setclm 
  PRIVATE :: exx_pot_derivative
  !
  !
CONTAINS
  !
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE exx_initialize()
      !! EXX initialization called in \(\texttt{init_run}\).
      !
      IMPLICIT NONE
      !
      include 'laxlib.fh'
      !
      INTEGER ::  i, iobtl, gindex_of_iobtl, irank, proc, tmp_iobtl, ndiag_n, ndiag_nx, ndiag_i
      REAL(DP) :: hx,hy,hz
      CHARACTER (len=300) :: print_str
      ndiag_i = ndiag_
      !
      ! Start of calculation banner...
      !
      WRITE(stdout,'(/,3X,"----------------------------------------------------")')
      WRITE(stdout,'(3X,"Exact Exchange Using Wannier Function Initialization")')
      WRITE(stdout,'(3X,"----------------------------------------------------")')
      !WRITE(stdout,'(/,3X,"The references for this algorithm are:",/ &
      !    &,5X,"(i)  theory: X. Wu , A. Selloni, and R. Car, Phys. Rev. B 79, 085102 (2009).",/ &
      !    &,5X,"(ii) implementation: H.-Y. Ko, B. Santra, R. A. DiStasio, L. Kong, Z. Li, X. Wu, and R. Car, arxiv.")')
      !
      WRITE(stdout,'(/,3X,"Parallelization info :")')
      WRITE(stdout,'(5X,"electronic states   ",3X,I7)') nbsp
      WRITE(stdout,'(5X,"MPI tasks           ",3X,I7)') nproc_image 
#if defined(_OPENMP)
      WRITE(stdout,'(5X,"OpenMP threads/MPI task",3X,I4)') omp_get_max_threads() 
#endif
      WRITE(stdout,'(5X,"Taskgroups          ",3X,I7)') fftx_ntgrp(dffts)
      !
      ! the fraction of exact exchange is stored here
      !
      exxalfa = xclib_get_exx_fraction()
      !
      IF(nbeg.LT.0) THEN
        !
        WRITE(stdout,'(/,3X,"This is an exact exchange calculation from scratch.",/ &
            &,3X,"To generate Wannier functions, exact exchange will",/ & 
            &,3X,"be turned off in the beginning ...  ")')
        !
        CALL stop_exx()  ! exx_is_active = .FALSE.
        !
      ELSE 
        !
        CALL start_exx() ! exx_is_active = .TRUE.
        !
      END IF
      !
      ! print few key input parameters used in exx calculations ...  
      !
      !WRITE(stdout,'(/,3X,"Below are the user defined parameters controlling the accuracy of the exact exchange energy")')
      !
      WRITE(stdout,'(/,3X,"parameters used in exact exchange calculation",/ &
          &,5X,"radius to compute pairs:",F6.1,1X,"A.U.",3X,"maximum number of neighbors:",I5)')dis_cutoff,neigh
      !
      IF (texx_cube) THEN
        WRITE(stdout,'(/,3X,"parameters used to solve Poisson equation",/ &
            &,5X,"cube side length for self potential:",F6.1,1X,"A.U."   &
            &,3X,"cube side length for pair potential:",F6.1,1X,"A.U.",/ &
            &,5X,"Poisson solver discretized using",I4,3X,"points in each dimension")')2*exx_ps_rcut_s,2*exx_ps_rcut_p,2*nord2+1
        !
        WRITE(stdout,'(/,3X,"parameters used for multipole expansion",/   &
            &,5X,"cube side length for self potential:",F6.1,1X,"A.U."   &
            &,3X,"cube side length for pair potential:",F6.1,1X,"A.U.",/ &
            &,5X,"maximum angular momentum:",I4)')2*exx_me_rcut_s,2*exx_me_rcut_p,lmax
      ELSE
        WRITE(stdout,'(/,3X,"parameters used to solve Poisson equation",/ &
            &,5X,"radius for self potential:",F6.1,1X,"A.U."   &
            &,3X,"radius for pair potential:",F6.1,1X,"A.U.",/ &
            &,5X,"Poisson solver discretized using",I4,3X,"points in each dimension")')exx_ps_rcut_s,exx_ps_rcut_p,2*nord2+1
        !
        WRITE(stdout,'(/,3X,"parameters used for multipole expansion",/   &
            &,5X,"radius for self potential:",F6.1,1X,"A.U."   &
            &,3X,"radius for pair potential:",F6.1,1X,"A.U.",/ &
            &,5X,"maximum angular momentum:",I4)')exx_me_rcut_s,exx_me_rcut_p,lmax
      END IF
      !
      ! Error messages for inconsistencies with current version of code...
      !
      IF(exx_ps_rcut_p.GT.exx_ps_rcut_s) CALL errore('exx_module','EXX calculation error :  &
          & The exx_ps_rcut_pair should be set smaller than or equal to the exx_ps_rcut_self',1)
      !
      IF(exx_ps_rcut_p.GE.exx_me_rcut_p) CALL errore('exx_module','EXX calculation error :  &
          & The exx_ps_rcut_pair should be set smaller than the exx_me_rcut_pair',1)
      !
      IF(exx_ps_rcut_p.GE.exx_me_rcut_s) CALL errore('exx_module','EXX calculation error :  &
          & The exx_ps_rcut_pair should be set smaller than the exx_me_rcut_self',1)
      !
      IF(exx_me_rcut_p.GT.exx_me_rcut_s) CALL errore('exx_module','EXX calculation error :  &
          & The exx_me_rcut_pair should be set smaller than or equal to the exx_me_rcut_self',1)
      !
      IF(exx_ps_rcut_s.GE.exx_me_rcut_s) CALL errore('exx_module','EXX calculation error :  &
          & The exx_ps_rcut_self should be set smaller than the exx_me_rcut_self',1)
#ifdef __CUDA
      IF(.not. texx_cube) CALL errore('exx_module','EXX calculation error :  &
          & Only cubic subdomain implemented in Exx',1)
#endif

      !
      hx=DSQRT(h(1,1)*h(1,1)+h(2,1)*h(2,1)+h(3,1)*h(3,1))
      hy=DSQRT(h(1,2)*h(1,2)+h(2,2)*h(2,2)+h(3,2)*h(3,2))
      hz=DSQRT(h(1,3)*h(1,3)+h(2,3)*h(2,3)+h(3,3)*h(3,3))
      IF(2*exx_me_rcut_s.GT.MIN(hx,hy,hz)) CALL errore('exx_module','EXX calculation error :  &
          & The exx_me_rcut_self should be set smaller than half the minimum cell length',1)
      !
      IF(fftx_ntgrp(dffts).GT.1) CALL errore('exx_module','EXX calculation error : &
          & taskgroup (-ntg) > 1 needed for zeta>1 calculations currently broken and will&
          & be fixed in an up-coming major update. Please contact Robert A. DiStasio Jr.&
          & (distasio@cornell.edu) if you should need assistance reverting to an earlier&
          & version with working taskgroup supoprt.',1)
      !
      hx=DSQRT(h(1,1)*h(1,1)+h(2,1)*h(2,1)+h(3,1)*h(3,1))
      hy=DSQRT(h(1,2)*h(1,2)+h(2,2)*h(2,2)+h(3,2)*h(3,2))
      hz=DSQRT(h(1,3)*h(1,3)+h(2,3)*h(2,3)+h(3,3)*h(3,3))
      IF(2*exx_me_rcut_s.GT.MIN(hx,hy,hz)) CALL errore('exx_module','EXX calculation error :  &
          & The exx_me_rcut_self should be set smaller than half the minimum cell length',1)
      !
      IF(fftx_ntgrp(dffts).GT.1) CALL errore('exx_module','EXX calculation error : &
          & taskgroup (-ntg) > 1 needed for zeta>1 calculations currently unavailable and will&
          & become available in the next major update. Please contact Robert A. DiStasio Jr.&
          & (distasio@cornell.edu) if you should need assistance reverting to an earlier&
          & version with working taskgroup support.',1)
      !
      IF(nproc_image.GE.nbsp) THEN
        !
        IF(MOD(nproc_image,nbsp).NE.0) CALL errore('exx_module','EXX calculation error :  &
            & number of MPI tasks has to be integer multiple of number of electronic states ',1)
        !
      END IF
      !
      IF(nproc_image.GT.nbsp.AND.nproc_image.LT.dfftp%nr3) CALL errore('exx_module',&
          'EXX calculation error :  for systems with few electronic states (less than nr3x)&
          & the number of MPI tasks has to be less or equal to the number of electronic states ',1)
      !
      IF(nproc_image.GT.nbsp.AND.MOD(nbsp,2).NE.0) CALL errore('exx_module','EXX calculation error :  &
          & for odd number of electronic states the number of MPI tasks has to be less or&
          & equal to the number of electronic states ',1)
      !
      IF(nspin.EQ.1) THEN 
        !
        ndiag_n=(INT(DSQRT(DBLE(idesc(LAX_DESC_N,1)))))**2
        !
        IF(nproc_image.GT.idesc(LAX_DESC_N,1).AND.ndiag_i.EQ.0) THEN
          !
          WRITE(print_str,'(3X,"EXX calculation error : use -ndiag N option in the execution of cp.x. &
              & Set N to any perfect square number which is equal to or less than the number of electronic states &
              & ***** Use -ndiag",1X,I7)') ndiag_n
          !
          CALL errore('exx_module',print_str,1)
          !
        END IF
        !
      ELSE IF(nspin.EQ.2) THEN
        !
        ! ndiag_n: suggested optimal ndiag value
        !
        ndiag_nx=MIN(idesc(LAX_DESC_N,1),idesc(LAX_DESC_N,2))
        ndiag_n=(INT(DSQRT(DBLE(ndiag_nx))))**2
        !
        IF(nproc_image.GT.ndiag_nx.AND.ndiag_i.EQ.0) THEN
          !
          WRITE(print_str,'(3X,"EXX calculation error : use -ndiag N option in the execution of cp.x. &
          & Set N to any perfect square number which is equal to or less than the number of electronic states &
          & ***** Use -ndiag",1X,I7)') ndiag_n
          !
          CALL errore('exx_module',print_str,1)
          !
        END IF
        !
      END IF      
      !
      IF(fftx_ntgrp(dffts).GT.1) CALL errore('exx_module','EXX calculation error : taskgroup no longer supported for exx.',1)
      IF((nproc_image.LE.nbsp).AND.(fftx_ntgrp(dffts).GT.1)) CALL errore('exx_module','EXX calculation error :  &
          & use taskgroup (-ntg) = 1 when number of MPI tasks is less or equal to the number of electronic states',1)
      !
      ! to fix this issue. see file exx_psi.f90, exx_gs.f90
      IF(nproc_image.GT.nbsp.AND.MOD(dffts%nnr,fftx_ntgrp(dffts)).NE.0) CALL errore('exx_module','EXX calculation error : &
          & (nr1x * nr2x) is not integer multiple of the number of task groups. Change task groups such that &
          & (nr1x * nr2x) becomes integer multiple of the number of task groups. Otherwise restrict number of MPI tasks &
          & up to the number of electronic states.',1)
      !
      ! to fix this issue. see file exx_psi.f90, exx_gs.f90
      IF((nproc_image.GT.nbsp).AND.(MOD(dffts%nnr,2).NE.0)) CALL errore('exx_module','EXX calculation error : &
          & (nr1x * nr2x) is an odd number. The calculation will work when the number of MPI tasks is smaller than&
          & or equal to the electronic bands. Otherwise, change ecutwfc to make (nr1x * nr2x) an even number.',1)
      !
      ! to fix this issue. see file exx_psi.f90, exx_gs.f90
      IF((nproc_image.GT.nbsp).AND.MOD(nbsp,2*fftx_ntgrp(dffts)).NE.0) CALL errore('exx_module','EXX calculation error : &
          & number of electronic states is not integer multiple of two times the number of task groups. &
          & Either change the number of taskgroups or restrict number of MPI tasks up to the number of electronic states.',1)
      !
      ! to fix this issue. see file exx_psi.f90, exx_gs.f90
      IF(nproc_image.GT.nbsp) THEN
        !
        IF (NINT(2**(LOG(DBLE(INT(nproc_image / dfftp%nr3))) / LOG(2.0))).EQ.1) THEN
          !
          ! NINT(2**(LOG(DBLE(INT(nproc_image / dfftp%nr3))) / LOG(2.0))) is the largest possible task group that one may use in this implementation
          !
          write(stdout,*) 
          write(stdout,*) "**********************************************************************************************"
          write(stdout,*) "*****************************   EXX PARALLELIZATION SUGGESTION   *****************************"
          write(stdout,*) "**********************************************************************************************"
          write(stdout,*) "You may want to use number of MPI tasks = ", CEILING(DBLE(2.0*dfftp%nr3)/DBLE(nbsp))*nbsp,& 
            & "(combined with -ntg 2)"
          !
        ELSE IF (NINT(2**(LOG(DBLE(INT(nproc_image / dfftp%nr3))) / LOG(2.0))).GT.fftx_ntgrp(dffts)) THEN
          !
          write(stdout,*) 
          write(stdout,*) "**********************************************************************************************"
          write(stdout,*) "*****************************   EXX PARALLELIZATION SUGGESTION   *****************************"
          write(stdout,*) "**********************************************************************************************"
          write(stdout,*) "You may want to change number of taskgroups (-ntg) to ", &
            NINT(2**(LOG(DBLE(INT(nproc_image / dfftp%nr3))) / LOG(2.0)))
        END IF
        !
        IF(fftx_ntgrp(dffts).EQ.1) THEN
          !
          write(stdout,*) 
          write(stdout,*) "**********************************************************************************************"
          write(stdout,*) "*****************************         Possible Solutions         *****************************"
          write(stdout,*) "**********************************************************************************************"
          write(stdout,*) "For the following error printing, we offer the following solutions:"
          write(stdout,*) " 1. (generic) Use numer of MPI tasks equal to the number of electronic states and -ntg 1"
          write(stdout,*) " 2. (advance) Tune numer of MPI tasks as integral multiple of the number of electronic states "
          write(stdout,*) "              such that this number is greater than or equal to the nr3x * (Taskgroups number)"
          write(stdout,*) "              Notice that the (Taskgroups number) should be at least 2."
          write(stdout,*) "              (See EXX PARALLELIZATION SUGGESTION above)"
          !
          CALL errore('exx_module','EXX calculation error : &
            & One needs number of task groups =  2^n where n is a positive integer when number of MPI tasks is greater than &
            & the number of electronic states. See above for Possible Solutions',1)
          !
        ELSE IF (NINT(2**(LOG(DBLE(fftx_ntgrp(dffts))) / LOG(2.0))).NE.fftx_ntgrp(dffts)) THEN
          !
          ! NINT(2**(LOG(DBLE(fftx_ntgrp(dffts))) / LOG(2.0))) is the largest power of 2 that is smaller or equal to dffts%nogrp
          !
          CALL errore('exx_module','EXX calculation error : &
            & One needs number of task groups =  2^n where n is a positive integer when number of MPI tasks is greater than &
            & the number of electronic states.',1)
        END IF
        !
      END IF
      !
      IF((fftx_ntgrp(dffts).GT.1).AND.(dfftp%nr3*fftx_ntgrp(dffts).GT.nproc_image)) &
          & CALL errore('exx_module','EXX calculation error : &
          & (nr3x * number of taskgroups) is greater than the number of MPI tasks. Change the number of MPI tasks or the number &
          & of taskgroups or both. To estimate ntg, find the value of nr3x in the output and compute (MPI task/nr3x) and take &
          & the integer value.',1)
      ! 
      IF ( lwfpbe0nscf ) THEN
        CALL errore('exx_module','Non self-consistent EXX calculation is possibly not working',1)
      END IF
      !
      !Index of EXX step is initialized ..
      ! 
      n_exx=0
      !
      !Grid sizes are initialized ..
      !
      nr1=dfftp%nr1; nr2=dfftp%nr2; nr3=dfftp%nr3
      nr1r=nr1/2; nr2r=nr2/2; nr3r=nr3/2
      IF(MOD(nr1,2).EQ.1) nr1r=(nr1+1)/2
      IF(MOD(nr2,2).EQ.1) nr2r=(nr2+1)/2
      IF(MOD(nr3,2).EQ.1) nr3r=(nr3+1)/2
      !
      ! calculates clm = (l-m)!/(l+m)!
      !
      ALLOCATE( clm(0:lmax, 0:lmax) )
      CALL setclm(lmax, clm)
      !
      ! calculate the finite difference coefficients for 1st and 2nd derivatives
      !   coe_1st_derv and coeke are computed by subroutine fornberg
      !   which is called in exx_gs.f90
      !
      ALLOCATE( coe_1st_derv(-nord1:nord1, 3)) ! coe_1st_derv(neighbor, d/di)
      ALLOCATE( coeke(-nord2:nord2, 3,3))      ! coeke(neighbor, d/di, d/dj)
#ifdef __CUDA
      ALLOCATE(coe_1st_derv_d, source=coe_1st_derv)
      ALLOCATE(coeke_d,   source=coeke)
#endif
      !
      IF (texx_cube) then
        CALL exx_initialize_cube
      ELSE
        CALL exx_initialize_sphere
      END IF
      !
      CALL exx_initialize_common
      !
      RETURN
      !
  END SUBROUTINE exx_initialize
  !--------------------------------------------------------------------------------------------------------------

  SUBROUTINE  exx_initialize_cube()
    IMPLICIT NONE
    REAL(DP) :: hx, hy, hz                             !grid spacing along lattice directions
    ALLOCATE( coemicf(-nord2:nord2, 3,3))    ! coeke(neighbor, d/di, d/dj)
#ifdef __CUDA
    ALLOCATE(coemicf_d, source=coemicf)
#endif
    nrg(1)=nr1;   nrg(2)=nr2;   nrg(3)=nr3
    nrgr(1)=nr1r; nrgr(2)=nr2r; nrgr(3)=nr3r
    !
    ! ALLOCATE exx_potential
    !
    IF (nproc_image .LT. nbsp) THEN
      !
      ALLOCATE( exx_potential(dffts%nr1*dffts%nr2*dffts%my_nr3p,nbsp) )
      !
    ELSE
      !
      IF ( dffts%has_task_groups ) THEN
        !
        ALLOCATE( exx_potential(dffts%nnr,nproc_image/fftx_ntgrp(dffts)) )
        !
        IF(MOD(nproc_image,fftx_ntgrp(dffts)).NE.0) CALL errore &
          & ('exx_module','EXX calculation is not working when &
          & number of MPI tasks (nproc_image) is not integer multiple of number of taskgroups',1)
        !
      ELSE
        !
        ALLOCATE( exx_potential(dffts%nr1x*dffts%nr2x*dffts%my_nr3p,nproc_image) ) !
        !
      END IF
      !
    END IF
    !
    exx_potential=0.0_DP
    !
    ! Compute number of grid points in Poisson and multipole spheres and store
    ! information to interchange between global and local (in sphere) grid indices ..
    !
    IF(thdyn) THEN
      !
      ! the number of grid points in Poisson and multipole spheres are computed
      ! using the simulation cell given in the input to keep those numbers
      ! constant throughout a variable cell dynamics ...
      !
      h_in=h_init
      !
      WRITE(stdout,'(/,3X,"This is a variable cell calculation. In this implementation",/ &
        &,3X,"the number of grid points used in Poisson cube and",/ &
        &,3X,"multipole expansion cube are kept constant as computed",/ &
        &,3X,"from the simulation cell given in the input...")')
      !
      WRITE(stdout,'(/3X,"cell from input:")')
      WRITE(stdout,'(3X,3F12.6)') h_in
      !
      ! Alternatively, reference cell could be used ... 
      !
      !IF(ref_cell.EQV..TRUE.) h_in=ref_at*ref_alat 
      !
    ELSE
      !
      ! for constant volume calculation ... 
      !
      h_in=h
      !
    END IF
    !
    !--------------------------------------------------------------------------
    ! calculate grid spaing along each direction
    !--------------------------------------------------------------------------
    hx=DSQRT(h(1,1)*h(1,1)+h(2,1)*h(2,1)+h(3,1)*h(3,1))/nr1
    hy=DSQRT(h(1,2)*h(1,2)+h(2,2)*h(2,2)+h(3,2)*h(3,2))/nr2
    hz=DSQRT(h(1,3)*h(1,3)+h(2,3)*h(2,3)+h(3,3)*h(3,3))/nr3
    !--------------------------------------------------------------------------
    ! Number of grid points inside Poisson (PS) and multipole expansion (ME)
    ! subdomains. _s: self; _p: pair
    ! We make the side length of the cube as the diameter of the spherical
    ! subdomain
    !--------------------------------------------------------------------------
    PScubeSL_s = IDNINT((2*exx_ps_rcut_s/MIN(hx,hy,hz))/8.0D0)*8-1
    PScubeSL_p = IDNINT((2*exx_ps_rcut_p/MIN(hx,hy,hz))/8.0D0)*8-1
    MEcubeSL_s = IDNINT((2*exx_me_rcut_s/MIN(hx,hy,hz))/8.0D0)*8-1
    MEcubeSL_p = IDNINT((2*exx_me_rcut_p/MIN(hx,hy,hz))/8.0D0)*8-1
    !--------------------------------------------------------------------------
    s_ps_r(1) = nr1r - (PScubeSL_s-1)/2;
    s_ps_r(2) = nr2r - (PScubeSL_s-1)/2;
    s_ps_r(3) = nr3r - (PScubeSL_s-1)/2;
    s_ps_r(4) = nr1r - (PScubeSL_s-1)/2 + (PScubeSL_s-1);
    s_ps_r(5) = nr2r - (PScubeSL_s-1)/2 + (PScubeSL_s-1);
    s_ps_r(6) = nr3r - (PScubeSL_s-1)/2 + (PScubeSL_s-1);
    !--------------------------------------------------------------------------
    p_ps_r(1) = nr1r - (PScubeSL_p-1)/2;
    p_ps_r(2) = nr2r - (PScubeSL_p-1)/2;
    p_ps_r(3) = nr3r - (PScubeSL_p-1)/2;
    p_ps_r(4) = nr1r - (PScubeSL_p-1)/2 + (PScubeSL_p-1);
    p_ps_r(5) = nr2r - (PScubeSL_p-1)/2 + (PScubeSL_p-1);
    p_ps_r(6) = nr3r - (PScubeSL_p-1)/2 + (PScubeSL_p-1);
    !--------------------------------------------------------------------------
    s_me_r(1) = nr1r - (MEcubeSL_s-1)/2;
    s_me_r(2) = nr2r - (MEcubeSL_s-1)/2;
    s_me_r(3) = nr3r - (MEcubeSL_s-1)/2;
    s_me_r(4) = nr1r - (MEcubeSL_s-1)/2 + (MEcubeSL_s-1);
    s_me_r(5) = nr2r - (MEcubeSL_s-1)/2 + (MEcubeSL_s-1);
    s_me_r(6) = nr3r - (MEcubeSL_s-1)/2 + (MEcubeSL_s-1);
    !--------------------------------------------------------------------------
    p_me_r(1) = nr1r - (MEcubeSL_p-1)/2;
    p_me_r(2) = nr2r - (MEcubeSL_p-1)/2;
    p_me_r(3) = nr3r - (MEcubeSL_p-1)/2;
    p_me_r(4) = nr1r - (MEcubeSL_p-1)/2 + (MEcubeSL_p-1);
    p_me_r(5) = nr2r - (MEcubeSL_p-1)/2 + (MEcubeSL_p-1);
    p_me_r(6) = nr3r - (MEcubeSL_p-1)/2 + (MEcubeSL_p-1);
    !--------------------------------------------------------------------------
    n_s_ps = PScubeSL_s**3
    n_p_ps = PScubeSL_p**3
    n_s_me = MEcubeSL_s**3
    n_p_me = MEcubeSL_p**3
    !--------------------------------------------------------------------------
    !
    ALLOCATE( me_cs(3,       s_me_r(1):s_me_r(4),s_me_r(2):s_me_r(5),s_me_r(3):s_me_r(6)), stat=ierr)
    ALLOCATE( me_rs(0:lmax,  s_me_r(1):s_me_r(4),s_me_r(2):s_me_r(5),s_me_r(3):s_me_r(6)), stat=ierr)
    ALLOCATE( me_ri(0:lmax+1,s_me_r(1):s_me_r(4),s_me_r(2):s_me_r(5),s_me_r(3):s_me_r(6)), stat=ierr)
    ALLOCATE( me_rc(0:lmax,  s_me_r(1):s_me_r(4),s_me_r(2):s_me_r(5),s_me_r(3):s_me_r(6)), stat=ierr)
    me_cs=0.0_DP; me_rs=0.0_DP; me_ri=0.0_DP; me_rc=0.0_DP
#ifdef __CUDA
      ALLOCATE(me_cs_d, source=me_cs)
      ALLOCATE(me_rs_d, source=me_rs)
      ALLOCATE(me_ri_d, source=me_ri)
      ALLOCATE(me_rc_d, source=me_rc)
#endif
    !
    WRITE(stdout,'(/,3X,"number of grid points in Poisson cube ",/ &
      &,5X,"self potential:",I8,3X,"pair potential:",I8)') n_s_ps, n_p_ps
    WRITE(stdout,'(/,3X,"number of grid points in multipole expansion cube: ",/ &
      &,5X,"self potential:",I8,3X,"pair potential:",I8)') n_s_me, n_p_me
    RETURN
  END SUBROUTINE exx_initialize_cube

  SUBROUTINE  exx_initialize_sphere()
    IMPLICIT NONE
    !
    ! ALLOCATE exx_potential
    !
    IF (nproc_image .LT. nbsp) THEN
      !
      ALLOCATE( exx_potential(dffts%nr1*dffts%nr2*dffts%my_nr3p,nbsp) )
      !
    ELSE
      !
      IF ( dffts%has_task_groups ) THEN
        !
        ALLOCATE( exx_potential(dffts%nnr,nproc_image/fftx_ntgrp(dffts)) )
        !
        IF(MOD(nproc_image,fftx_ntgrp(dffts)).NE.0) CALL errore &
          & ('exx_module','EXX calculation is not working when &
          & number of MPI tasks (nproc_image) is not integer multiple of number of taskgroups',1)
        !
      ELSE
        !
        ALLOCATE( exx_potential(dffts%nr1x*dffts%nr2x*dffts%my_nr3p,nproc_image) ) !
        !
      END IF
      !
    END IF
    !
    exx_potential=0.0_DP
    !
    ! Compute number of grid points in Poisson and multipole spheres and store
    ! information to interchange between global and local (in sphere) grid indices ..
    !
    IF(thdyn) THEN
      !
      ! the number of grid points in Poisson and multipole spheres are computed
      ! using the simulation cell given in the input to keep those numbers
      ! constant throughout a variable cell dynamics ...
      !
      h_in=h_init
      !
      WRITE(stdout,'(/,3X,"This is a variable cell calculation. In this implementation",/ &
        &,3X,"the number of grid points used in Poisson sphere and",/ &
        &,3X,"multipole expansion sphere are kept constant as computed",/ &
        &,3X,"from the simulation cell given in the input...")')
      !
      WRITE(stdout,'(/3X,"cell from input:")')
      WRITE(stdout,'(3X,3F12.6)') h_in
      !
      ! Alternatively, reference cell could be used ... 
      !
      !IF(ref_cell.EQV..TRUE.) h_in=ref_at*ref_alat 
      !
    ELSE
      !
      ! for constant volume calculation ... 
      !
      h_in=h
      !
    END IF
    !
    ! get number of points inside the Poisson and (Poisson + Multipole expansion) spheres
    ! for self and pair spheres (using different cutoff radii)
    !
    CALL getnpinsp(exx_ps_rcut_s, exx_me_rcut_s, np_in_sp_s, np_in_sp_me_s )
    CALL getnpinsp(exx_ps_rcut_p, exx_me_rcut_p, np_in_sp_p, np_in_sp_me_p )
    !
    ! exx_setup   orders the local index in the order such that
    !                      1 ... np_in_sp_s      is in self Poisson sphere
    !             np_in_sp_s ... np_in_sp_me_s   is in self Multipole sphere
    !                      1 ... np_in_sp_p      is in pair Poisson sphere
    !             np_in_sp_p ... np_in_sp_me_p   is in self Multipole sphere
    !
    !             and computes the x y z values for multipole expansion
    !          ( notice that it is still xyz even in non-orthogonal cells)
    !
    CALL exx_setup( )
    !
    WRITE(stdout,'(/,3X,"number of grid points in Poisson spehere ",/ &
      &,5X,"self potential:",I8,3X,"pair potential:",I8)')np_in_sp_s,np_in_sp_p
    WRITE(stdout,'(/,3X,"number of grid points in multipole expansion spehere: ",/ &
      &,5X,"self potential:",I8,3X,"pair potential:",I8)')np_in_sp_me_s,np_in_sp_me_p

    RETURN
  END SUBROUTINE exx_initialize_sphere

  SUBROUTINE exx_initialize_common()
    IMPLICIT NONE
    !
    ! Variables for parallelization are initialized ....
    !
    ! parallelization scaling factors
    !
    IF(nproc_image.LE.nbsp) THEN
      sc_fac = 1
    ELSE
      sc_fac = nproc_image/nbsp
    END IF
    !
    ! my_nbspx is the maxval of my_nbsp(:), ie the maximum number of bands a processor (mpi task) can have
    !
    IF(nproc_image .LE. nbsp) THEN
      my_nbspx   = nbsp / nproc_image
      IF( MOD(nbsp, nproc_image) /= 0)THEN
        my_nbspx = my_nbspx + 1
      ENDIF
    ELSE
      my_nbspx   = (nbsp / nproc_image) + 1
    END IF
    !
    !print *, 'my_nbspx =', my_nbspx
    !
    ALLOCATE( my_nxyz ( nproc_image ) )
    ALLOCATE( my_nbsp ( nproc_image ) )
    !
    IF(nproc_image .LE. nbsp) THEN
      my_nbsp(:) = nbsp/nproc_image
      IF( MOD(nbsp, nproc_image) /= 0)THEN
        DO i = 1, nproc_image
          IF( (i-1) < MOD(nbsp, nproc_image) ) my_nbsp(i) = my_nbsp(i)+1
        END DO
      END IF
    ELSE
      my_nbsp(:) = 1
    END IF
    !
    ! ** Note that in this case .. 
    ! this is incorrect:   my_nxyz(:) = nr1*nr2*dffts%npp(me_image+1)
    my_nxyz(:) = nr1*nr2*dffts%nr3p
    !
    !DEBUG
    !WRITE(stdout,'("my_nbsp")')
    !WRITE(stdout,'(20I5)')my_nbsp
    !WRITE(stdout,'("my_nxyz")')
    !WRITE(stdout,'(20I7)')my_nxyz
    !DEBUG
    !
    ALLOCATE( index_my_nbsp (my_nbspx, nproc_image) )
    ALLOCATE( rk_of_obtl ( nbsp ) )
    ALLOCATE( lindex_of_obtl( nbsp ) )
    !
    IF(nproc_image .LE. nbsp) THEN
      index_my_nbsp(:, :) = nbsp + 1 ! set orbital index to beyond nbsp to indicate this entry is not used...
      DO irank = 1, nproc_image
        DO iobtl = 1, my_nbsp(irank)
          gindex_of_iobtl = iobtl
          DO proc = 1, irank - 1, 1
            gindex_of_iobtl = gindex_of_iobtl + my_nbsp(proc)
          END DO
          IF( gindex_of_iobtl <= nbsp) THEN
            index_my_nbsp(iobtl, irank) = gindex_of_iobtl
          END IF
        END DO
      END DO
    ELSE
      DO proc = 1, nproc_image
        index_my_nbsp(1, proc) = proc
      END DO
    END IF
    !
    !DEBUG
    !WRITE(stdout,'("index_my_nbsp")')
    !WRITE(stdout,'(20I5)')index_my_nbsp(1,:)
    !WRITE(stdout,'(20I5)')index_my_nbsp(2,:)
    !DEBUG
    !
    DO iobtl = 1, nbsp
      rk_of_obtl(iobtl) = 0
      tmp_iobtl = iobtl
      DO proc = 1, nproc_image
        tmp_iobtl = tmp_iobtl - my_nbsp(proc)
        IF (tmp_iobtl <= 0) THEN
          rk_of_obtl(iobtl) = proc - 1
          !print *, 'lrk_of_iobtl=', proc-1, rk_of_obtl(iobtl) 
          EXIT
        END IF
      END DO
    END DO
    !
    !DEBUG
    !WRITE(stdout,'("rk_of_obtl")')
    !WRITE(stdout,'(20I5)')rk_of_obtl
    !DEBUG
    !
    DO iobtl = 1, nbsp
      lindex_of_obtl(iobtl) = iobtl
      DO proc = 1, nproc_image
        IF (lindex_of_obtl(iobtl) <= my_nbsp(proc)) EXIT
        lindex_of_obtl(iobtl) = lindex_of_obtl(iobtl) - my_nbsp(proc)
      END DO
    END DO
    !
    !DEBUG
    !WRITE(stdout,'("lindex_of_obtl")')
    !WRITE(stdout,'(20I5)')lindex_of_obtl
    !DEBUG
    !
    ! Allocate other variables ....
    !
    IF ( lwfpbe0nscf .or. lwfnscf ) ALLOCATE( rhopr( dfftp%nnr, nspin ) )
    !
    IF ( lwfpbe0nscf ) THEN
      ALLOCATE( vwc(3, vnbsp) )
      ALLOCATE( pair_dist( 3, neigh, my_nbspx), stat=ierr )
      pair_dist (:,:,:) = 0.0_DP
      ALLOCATE( pairv( np_in_sp_p, 3, neigh, my_nbspx), stat=ierr )
      pairv (:,:,:,:) = 0.0_DP
    END IF
    !
    WRITE(stdout,'(/,3X,"----------------------------------------------------")')
    return
  end subroutine exx_initialize_common

  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE exx_finalize()
      !
      IF( ALLOCATED( clm )     )        DEALLOCATE( clm)
      IF( ALLOCATED( coeke)    )        DEALLOCATE( coeke)
      IF( ALLOCATED( coe_1st_derv)   )  DEALLOCATE( coe_1st_derv)
#ifdef __CUDA
      IF( ALLOCATED( coe_1st_derv_d)    )      DEALLOCATE( coe_1st_derv_d)
      IF( ALLOCATED( coeke_d)    )      DEALLOCATE( coeke_d)
      IF( ALLOCATED( coemicf_d)    )    DEALLOCATE( coemicf_d)
#endif
      IF( ALLOCATED( exx_potential ) )  DEALLOCATE( exx_potential )
      IF( ALLOCATED( rhopr ) )          DEALLOCATE( rhopr )
      IF( ALLOCATED( vwc)    )          DEALLOCATE( vwc )
      IF( ALLOCATED( selfv ) )          DEALLOCATE( selfv )
      IF( ALLOCATED( pairv ) )          DEALLOCATE( pairv )
      IF( ALLOCATED( pair_dist ) )      DEALLOCATE( pair_dist )
      IF( ALLOCATED( pair_label ) )     DEALLOCATE( pair_label )
      IF( ALLOCATED( pair_step ) )      DEALLOCATE( pair_step )
      IF( ALLOCATED( pair_status ) )    DEALLOCATE( pair_status )
      IF( ALLOCATED( my_nxyz ) )        DEALLOCATE( my_nxyz )
      IF( ALLOCATED( my_nbsp ) )        DEALLOCATE( my_nbsp )
      IF( ALLOCATED( index_my_nbsp ) )  DEALLOCATE( index_my_nbsp)
      IF( ALLOCATED( rk_of_obtl ) )     DEALLOCATE( rk_of_obtl)
      IF( ALLOCATED( lindex_of_obtl ) ) DEALLOCATE( lindex_of_obtl )
      !
      IF( ALLOCATED( odtothd_in_sp ) )  DEALLOCATE( odtothd_in_sp )
      IF( ALLOCATED( thdtood_in_sp ) )  DEALLOCATE( thdtood_in_sp )
      IF( ALLOCATED( thdtood  ))        DEALLOCATE( thdtood)
      IF( ALLOCATED( xx_in_sp ))        DEALLOCATE( xx_in_sp )
      IF( ALLOCATED( yy_in_sp ))        DEALLOCATE( yy_in_sp )
      IF( ALLOCATED( zz_in_sp ))        DEALLOCATE( zz_in_sp )
      IF( ALLOCATED( sc_xx_in_sp ))     DEALLOCATE( sc_xx_in_sp )
      IF( ALLOCATED( sc_yy_in_sp ))     DEALLOCATE( sc_yy_in_sp )
      IF( ALLOCATED( sc_zz_in_sp ))     DEALLOCATE( sc_zz_in_sp )
      IF (ALLOCATED(psime_pair_send))  DEALLOCATE(psime_pair_send)
      IF (ALLOCATED(psime_pair_recv))  DEALLOCATE(psime_pair_recv)
      IF (ALLOCATED( rho_ps ))          DEALLOCATE( rho_ps )
      IF (ALLOCATED( pot_ps ))          DEALLOCATE( pot_ps )
      IF( ALLOCATED( me_cs ) )          DEALLOCATE( me_cs )
      IF( ALLOCATED( me_rs ) )          DEALLOCATE( me_rs )
      IF( ALLOCATED( me_ri ) )          DEALLOCATE( me_ri )
      IF( ALLOCATED( me_rc ) )          DEALLOCATE( me_rc )
      IF( ALLOCATED( selfrho ) )        DEALLOCATE( selfrho )
      IF( ALLOCATED( pairrho ) )        DEALLOCATE( pairrho )
      IF( ALLOCATED( coemicf)    )      DEALLOCATE( coemicf)
#ifdef __CUDA
      IF( ALLOCATED( me_cs_d ) )          DEALLOCATE( me_cs_d )
      IF( ALLOCATED( me_rs_d ) )          DEALLOCATE( me_rs_d )
      IF( ALLOCATED( me_ri_d ) )          DEALLOCATE( me_ri_d )
      IF( ALLOCATED( me_rc_d ) )          DEALLOCATE( me_rc_d )
      IF (ALLOCATED(psime_pair_send_d))  DEALLOCATE(psime_pair_send_d)
      IF (ALLOCATED(psime_pair_recv_d))  DEALLOCATE(psime_pair_recv_d)
      IF(ALLOCATED(psi_d   ))  DEALLOCATE(psi_d   )
      IF(ALLOCATED(vpsil_d))  DEALLOCATE(vpsil_d)
#endif
      !
      RETURN
      !
  END SUBROUTINE exx_finalize
  !--------------------------------------------------------------------------------------------------------------

  !======================================================================
  ! Set up the real-space grid for exact exchange.
  ! Define two spheres, possion solver is called for potential inside the inner one 
  ! Multipole expansion is used for points inside the outer one and all other points
  ! are set to zero.
  !
  ! Adapted from PARSEC by Lingzhu Kong. http://parsec.ices.utexas.edu/
  !========================================================================
  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE getnpinsp( exx_ps_rcut, exx_me_rcut, np_in_sp, np_in_sp_me )
      !
      IMPLICIT NONE
      !
      REAL(DP) :: exx_ps_rcut, exx_me_rcut
      INTEGER  :: np_in_sp, np_in_sp_me
      !
      REAL(DP) :: x,y,z,dqs(3),dq(3),dist
      INTEGER  :: i,j,k,np_in_sp2
      ! --------------------------------------------------------------------
      np_in_sp=0;np_in_sp2=0;np_in_sp_me= 0
      !
      DO k = 1,nr3
        DO j = 1, nr2
          DO i =1, nr1
            !
            ! distances between Grid points and center of the simulation cell in S space
            ! center of the box is set to grid point at int(nr1/2), int(nr2/2), int(nr3/2) for every cell
            ! 
            dqs(1) = (DBLE(i)/DBLE(nr1)) - DBLE(INT(nr1/2))/DBLE(nr1)
            dqs(2) = (DBLE(j)/DBLE(nr2)) - DBLE(INT(nr2/2))/DBLE(nr2)
            dqs(3) = (DBLE(k)/DBLE(nr3)) - DBLE(INT(nr3/2))/DBLE(nr3)
            !
            ! Here we are computing distances between Grid points and center of the simulation cell, so no MIC is needed ...
            ! Compute distance between grid point and the center of the simulation cell in R space 
            !
            dq(1)=h_in(1,1)*dqs(1)+h_in(1,2)*dqs(2)+h_in(1,3)*dqs(3)   !r_i = h s_i
            dq(2)=h_in(2,1)*dqs(1)+h_in(2,2)*dqs(2)+h_in(2,3)*dqs(3)   !r_i = h s_i
            dq(3)=h_in(3,1)*dqs(1)+h_in(3,2)*dqs(2)+h_in(3,3)*dqs(3)   !r_i = h s_i
            !
            dist = DSQRT(dq(1)*dq(1)+dq(2)*dq(2)+dq(3)*dq(3))
            !
            IF (dist .LE. exx_ps_rcut) THEN
              !
              np_in_sp = np_in_sp + 1
              !
            ELSE IF (dist .LE. exx_me_rcut) THEN
              !
              np_in_sp2 = np_in_sp2 + 1
              !
            END IF
            !
            np_in_sp_me = np_in_sp+np_in_sp2
            !
          END DO !i
        END DO !j
      END DO !k
      !
      RETURN
  END SUBROUTINE getnpinsp
  !--------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE exx_setup( )
      !
      IMPLICIT NONE
      !
      ! ====================================================================
      ! INPUT VARIABLES
      !
      ! LOCAL VARIABLES
      INTEGER  :: np, npsp1, npsp2, npsp3, npsp4, np_tmp_1, np_tmp_2 
      REAL(DP) :: x,y,z, dq(3),dqs(3),dist, rcut_sp2, rcut_sp3
      INTEGER  :: i,j,k, ierr, tmp
      ! --------------------------------------------------------------------
      !
      !! This part would be necessary if number of points in PS and ME sphere
      !! change every step
      !IF(n_exx.GT.0) THEN 
      !    IF( ALLOCATED( odtothd_in_sp ) )  DEALLOCATE(odtothd_in_sp )
      !    IF( ALLOCATED( thdtood_in_sp ) )  DEALLOCATE(thdtood_in_sp )
      !    IF( ALLOCATED( thdtood  ))        DEALLOCATE(thdtood)
      !    IF( ALLOCATED( xx_in_sp ))        DEALLOCATE(xx_in_sp )
      !    IF( ALLOCATED( yy_in_sp ))        DEALLOCATE(yy_in_sp )
      !    IF( ALLOCATED( zz_in_sp ))        DEALLOCATE(zz_in_sp )
      !END IF
      !
      ALLOCATE( odtothd_in_sp(3, np_in_sp_me_s ), stat=ierr )
      ALLOCATE( thdtood_in_sp(nr1, nr2, nr3), stat=ierr )
      ALLOCATE( thdtood(nr1, nr2, nr3), stat=ierr )
      ALLOCATE( xx_in_sp(1:np_in_sp_me_s), stat=ierr )
      ALLOCATE( yy_in_sp(1:np_in_sp_me_s), stat=ierr )
      ALLOCATE( zz_in_sp(1:np_in_sp_me_s), stat=ierr )
      ALLOCATE( sc_xx_in_sp(1:np_in_sp_me_s), stat=ierr )
      ALLOCATE( sc_yy_in_sp(1:np_in_sp_me_s), stat=ierr )
      ALLOCATE( sc_zz_in_sp(1:np_in_sp_me_s), stat=ierr )
      !
      xx_in_sp=0.0_DP; sc_xx_in_sp=0.0_DP
      yy_in_sp=0.0_DP; sc_yy_in_sp=0.0_DP
      zz_in_sp=0.0_DP; sc_zz_in_sp=0.0_DP
      !
      thdtood_in_sp = 0; odtothd_in_sp = 0; thdtood = 0
      !
      !this is the cutoff acording to the size of sphere
      !
      rcut_sp2=MIN(exx_me_rcut_p,exx_ps_rcut_s)
      rcut_sp3=MAX(exx_me_rcut_p,exx_ps_rcut_s)
      np_tmp_1=MIN(np_in_sp_s,np_in_sp_me_p)
      np_tmp_2=MAX(np_in_sp_s,np_in_sp_me_p)
      np    = 0
      npsp1 = 0; npsp2 = 0; npsp3 = 0; npsp4 = 0
      !
      DO k = 1,nr3
        DO j = 1, nr2
          DO i =1, nr1
            !
            np = np + 1
            !
            thdtood(i,j,k) = np
            !
            ! distances between Grid points and center of the simulation cell in S space
            ! center of the box is set to grid point at int(nr1/2), int(nr2/2), int(nr3/2) for every cell
            ! 
            dqs(1) = (DBLE(i)/DBLE(nr1)) - DBLE(INT(nr1/2))/DBLE(nr1)
            dqs(2) = (DBLE(j)/DBLE(nr2)) - DBLE(INT(nr2/2))/DBLE(nr2)
            dqs(3) = (DBLE(k)/DBLE(nr3)) - DBLE(INT(nr3/2))/DBLE(nr3)
            !
            ! Here we are computing distances between Grid points and center of the simulation cell, so no MIC is needed ...
            ! Compute distance between grid point and the center of the simulation cell in R space 
            !
            dq(1)=h_in(1,1)*dqs(1)+h_in(1,2)*dqs(2)+h_in(1,3)*dqs(3)   !r_i = h s_i
            dq(2)=h_in(2,1)*dqs(1)+h_in(2,2)*dqs(2)+h_in(2,3)*dqs(3)   !r_i = h s_i
            dq(3)=h_in(3,1)*dqs(1)+h_in(3,2)*dqs(2)+h_in(3,3)*dqs(3)   !r_i = h s_i
            !
            dist = DSQRT(dq(1)*dq(1)+dq(2)*dq(2)+dq(3)*dq(3))
            !
            IF (dist .LE. exx_ps_rcut_p) THEN
              npsp1 = npsp1 + 1
              thdtood_in_sp(i,j,k)  = npsp1
              odtothd_in_sp(1,npsp1) = i
              odtothd_in_sp(2,npsp1) = j
              odtothd_in_sp(3,npsp1) = k
              !
              xx_in_sp(npsp1) = dq(1); sc_xx_in_sp(npsp1) = dqs(1)
              yy_in_sp(npsp1) = dq(2); sc_yy_in_sp(npsp1) = dqs(2)
              zz_in_sp(npsp1) = dq(3); sc_zz_in_sp(npsp1) = dqs(3)
              !
            ELSE IF (dist .LE. rcut_sp2) THEN
              !
              npsp2 = npsp2 + 1
              tmp = npsp2 + np_in_sp_p
              thdtood_in_sp(i,j,k)  = tmp
              odtothd_in_sp(1, tmp) = i
              odtothd_in_sp(2, tmp) = j
              odtothd_in_sp(3, tmp) = k
              !
              xx_in_sp(tmp) = dq(1); sc_xx_in_sp(tmp) = dqs(1)
              yy_in_sp(tmp) = dq(2); sc_yy_in_sp(tmp) = dqs(2)
              zz_in_sp(tmp) = dq(3); sc_zz_in_sp(tmp) = dqs(3)
              !
            ELSE IF (dist .LE. rcut_sp3) THEN
              !
              npsp3 = npsp3 + 1
              tmp = npsp3 + np_tmp_1
              thdtood_in_sp(i,j,k)  = tmp
              odtothd_in_sp(1, tmp) = i
              odtothd_in_sp(2, tmp) = j
              odtothd_in_sp(3, tmp) = k
              !
              xx_in_sp(tmp) = dq(1); sc_xx_in_sp(tmp) = dqs(1)
              yy_in_sp(tmp) = dq(2); sc_yy_in_sp(tmp) = dqs(2)
              zz_in_sp(tmp) = dq(3); sc_zz_in_sp(tmp) = dqs(3)
              !
            ELSE IF (dist .LE. exx_me_rcut_s) THEN
              !
              npsp4 = npsp4 + 1
              tmp = npsp4 + np_tmp_2
              thdtood_in_sp(i,j,k)  = tmp
              odtothd_in_sp(1, tmp) = i
              odtothd_in_sp(2, tmp) = j
              odtothd_in_sp(3, tmp) = k
              !
              xx_in_sp(tmp) = dq(1); sc_xx_in_sp(tmp) = dqs(1)
              yy_in_sp(tmp) = dq(2); sc_yy_in_sp(tmp) = dqs(2)
              zz_in_sp(tmp) = dq(3); sc_zz_in_sp(tmp) = dqs(3)
              !
            END IF
          END DO
        END DO
      END DO
      !
      IF ( npsp1 .NE. np_in_sp_p) THEN
        WRITE(stdout, *)&
            'number of points in the 1st shell does not match', npsp1, np_in_sp_p
        WRITE(stdout, *)'STOP in exx_setup'
        RETURN
      END IF
      !
      IF ( npsp2 .NE. np_tmp_1-np_in_sp_p) THEN
        WRITE(stdout,*)&
            'number of points in the 2nd shell does not match', npsp2, np_tmp_1-np_in_sp_p
        WRITE(stdout, *)'STOP in exx_setup'
        RETURN
      END IF
      !
      IF ( npsp3 .NE. np_tmp_2-np_tmp_1) THEN
        WRITE(stdout,*)&
            'number of points in the 3rd shell does not match', npsp3, np_tmp_2-np_tmp_1
        WRITE(stdout, *)'STOP in exx_setup'
        RETURN
      END IF
      !
      IF ( npsp4 .NE. np_in_sp_me_s-np_tmp_2) THEN
        WRITE(stdout,*)&
            'number of points in the 4th shell does not match', npsp4, np_in_sp_me_s-np_tmp_2
        WRITE(stdout, *)'STOP in exx_setup'
        RETURN
      END IF
      !
      RETURN
      !
  END SUBROUTINE exx_setup
  !--------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE exx_setup_nscf( nnrtot, lpole, clm, factor, wc, vwc, nbsp, vnbsp )
      ! HK: Broken at this moment need to be fixed
      !
      IMPLICIT NONE
      !
      ! ====================================================================
      ! INPUT VARIABLES
      INTEGER  nnrtot, lpole, nbsp, vnbsp
      REAL(DP) factor, wc(3, nbsp), vwc(3, vnbsp)
      REAL(DP) clm(0:lpole, 0:lpole)
      !
      ! ====================================================================
      !      integer  odtothd_in_sp(3,np_in_sp), odtothd_in_sp2(3,np_in_sp2)
      !      integer  thdtood_in_sp(nr1, nr2, nr3),thdtood(nr1,nr2,nr3)
      ! ====================================================================
      ! LOCAL VARIABLES
      INTEGER   np, npsp, npsp2
      REAL(DP)  x,y,z, dist, exx_ps_rcut, exx_me_rcut
      INTEGER   i,j,k, ierr, tmp,np_in_sp
      REAL(DP) :: lenA,lenB,lenC                       !length of lattice  A, B, C in a.u.
      REAL(DP) :: centerx,centery,centerz              !coordinate of the center of the simulation box
      REAL(DP) :: hx,hy,hz                             !grid spacing along lattice directions
      REAL(DP) :: Jim(3,3)                             ! [d/d x]        [d/d a]
      !                                      jacobian    |d/d y| = [J]  |d/d b|
      !                                                  [d/d z]        [d/d c]
      ! --------------------------------------------------------------------
      !
      exx_ps_rcut=exx_ps_rcut_s
      exx_me_rcut=exx_me_rcut_s
      np_in_sp=np_in_sp_s
      !
      ! should work for all type of cells
      !
      lenA=DSQRT(h(1,1)*h(1,1)+h(2,1)*h(2,1)+h(3,1)*h(3,1))
      lenB=DSQRT(h(1,2)*h(1,2)+h(2,2)*h(2,2)+h(3,2)*h(3,2))
      lenC=DSQRT(h(1,3)*h(1,3)+h(2,3)*h(2,3)+h(3,3)*h(3,3))
      ! 
      hx = lenA / nr1  !grid spacing in Lattice 1 direction
      hy = lenB / nr2  !grid spacing in Lattice 2 direction
      hz = lenC / nr3  !grid spacing in Lattice 3 direction
      !
      centerx = 0.5_DP * lenA 
      centery = 0.5_DP * lenB
      centerz = 0.5_DP * lenC
      !
      ! For points in the 1st sphere, one needs to know if its finite-difference neighbors
      ! are inside or outside. We set the thdtood_in_sp to be np_in_sp + 1 for outside neighbors
      ALLOCATE( odtothd_in_sp(3, np_in_sp_me_s ), stat=ierr )
      ALLOCATE( thdtood_in_sp(nr1, nr2, nr3), stat=ierr )
      ALLOCATE( thdtood(nr1, nr2, nr3), stat=ierr )
      ALLOCATE( xx_in_sp(1:np_in_sp_me_s), stat=ierr )
      ALLOCATE( yy_in_sp(1:np_in_sp_me_s), stat=ierr )
      ALLOCATE( zz_in_sp(1:np_in_sp_me_s), stat=ierr )
      !
      xx_in_sp=0.0_DP
      yy_in_sp=0.0_DP
      zz_in_sp=0.0_DP
      !
      thdtood_in_sp = 0
      odtothd_in_sp = 0
      thdtood = 0
      !
      np    = 0
      npsp  = 0
      npsp2 = 0
      do k = 1,nr3
        do j = 1, nr2
          do i =1, nr1
            np = np + 1
            thdtood(i,j,k) = np
            !
            x = i * hx -centerx
            y = j * hy -centery
            z = k * hz -centerz
            dist = sqrt(x*x + y*y + z*z)
            !
            if(dist .le. exx_ps_rcut)then
              npsp = npsp + 1
              thdtood_in_sp(i,j,k)  = npsp
              odtothd_in_sp(1,npsp) = i
              odtothd_in_sp(2,npsp) = j
              odtothd_in_sp(3,npsp) = k
              !
              xx_in_sp(npsp) = x
              yy_in_sp(npsp) = y
              zz_in_sp(npsp) = z
            elseif(dist .le. exx_me_rcut)then
              npsp2 = npsp2 + 1
              tmp = npsp2 + np_in_sp
              thdtood_in_sp(i,j,k)  = tmp
              odtothd_in_sp(1, tmp) = i
              odtothd_in_sp(2, tmp) = j
              odtothd_in_sp(3, tmp) = k
              !
              xx_in_sp(tmp) = x
              yy_in_sp(tmp) = y
              zz_in_sp(tmp) = z
              !
            endif
          enddo
        enddo
      enddo
      !
      write(6,*)' npsp in exx_setup =', npsp, npsp2
      if( npsp .ne. np_in_sp)then
        write(6, *)'number of points in the 1st sphere does not match', npsp, np_in_sp
        write(6, *)'STOP in exx_setup'
        return
      endif
      !
      if( npsp2+npsp .ne. np_in_sp_me_s)then
        write(6,*)'number of points in the 2nd sphere does not match', npsp2+npsp, np_in_sp_me_s
        write(6, *)'STOP in exx_setup'
        return
      endif
      if (texx_cube) then
        ALLOCATE( me_cs(3,       s_me_r(1):s_me_r(4),s_me_r(2):s_me_r(5),s_me_r(3):s_me_r(6)), stat=ierr)
        ALLOCATE( me_rs(0:lmax,  s_me_r(1):s_me_r(4),s_me_r(2):s_me_r(5),s_me_r(3):s_me_r(6)), stat=ierr)
        ALLOCATE( me_ri(0:lmax+1,s_me_r(1):s_me_r(4),s_me_r(2):s_me_r(5),s_me_r(3):s_me_r(6)), stat=ierr)
        ALLOCATE( me_rc(0:lmax,  s_me_r(1):s_me_r(4),s_me_r(2):s_me_r(5),s_me_r(3):s_me_r(6)), stat=ierr)
        me_cs=0.0_DP; me_rs=0.0_DP; me_ri=0.0_DP; me_rc=0.0_DP
      end if
      !
      !========================================================================
      !
      ! get hx*d/dx, hx^2*d^2/dx^2 stencil and cross coefficients
      !
      CALL fornberg(nord1,nord2,coe_1st_derv(:,1),coeke(:,1,1),coeke(:,1,2),ierr)
      !
      if (ierr .ne. 0) then
        write(6,*) ' ERROR: Wrong parameter in CALL of Fornberg'
        write(6,*) ' STOP in init_var'
        return
      endif
      !
      !      renormalize coekes with respect to the grid spacing
      !
      ! first derivative coefficients
      !
      coe_1st_derv(:,3) = coe_1st_derv(:,1)/hz    ! d/dz stencil
      coe_1st_derv(:,2) = coe_1st_derv(:,1)/hy    ! d/dy stencil
      coe_1st_derv(:,1) = coe_1st_derv(:,1)/hx    ! d/dx stencil
      !
      ! axial derivatives
      !
      coeke(:,3,3) = -coeke(:,1,1)/(hz*hz*factor) ! -d^2/dz^2/factor stencil
      coeke(:,2,2) = -coeke(:,1,1)/(hy*hy*factor) ! -d^2/dy^2/factor stencil
      coeke(:,1,1) = -coeke(:,1,1)/(hx*hx*factor) ! -d^2/dx^2/factor stencil
      !
      ! cross derivatives
      !
      coeke(:,2,3) = -coeke(:,1,2)/(hy*hz*factor) ! -d^2/dydz/factor stencil
      coeke(:,1,3) = -coeke(:,1,2)/(hx*hz*factor) ! -d^2/dydz/factor stencil
      coeke(:,1,2) = -coeke(:,1,2)/(hx*hy*factor) ! -d^2/dydz/factor stencil
      !
      ! J = transpose(ainv).(diag(a))
      !
      Jim(:,1) = ainv(1,:)*lenA
      Jim(:,2) = ainv(2,:)*lenB ! i={xyz}, m={abc}
      Jim(:,3) = ainv(3,:)*lenC
      !
      ! weigh coeke with the Jacobian --
      !
      ! axial derivatives
      !
      coeke(:,3,3) = (Jim(1,3)**2+Jim(2,3)**2+Jim(3,3)**2)*coeke(:,3,3)
      coeke(:,2,2) = (Jim(1,2)**2+Jim(2,2)**2+Jim(3,2)**2)*coeke(:,2,2)
      coeke(:,1,1) = (Jim(1,1)**2+Jim(2,1)**2+Jim(3,1)**2)*coeke(:,1,1)
      !
      ! cross derivatives
      !
      coeke(:,2,3) = 2.0_DP*(Jim(1,2)*Jim(1,3)+Jim(2,2)*Jim(2,3)+Jim(3,2)*Jim(3,3))*coeke(:,2,3)
      coeke(:,1,3) = 2.0_DP*(Jim(1,1)*Jim(1,3)+Jim(2,1)*Jim(2,3)+Jim(3,1)*Jim(3,3))*coeke(:,1,3)
      coeke(:,1,2) = 2.0_DP*(Jim(1,1)*Jim(1,2)+Jim(2,1)*Jim(2,2)+Jim(3,1)*Jim(3,2))*coeke(:,1,2)
      coeke(:,3,2) = coeke(:,2,3) ! symmetry of coeke
      coeke(:,2,1) = coeke(:,1,2) ! symmetry of coeke
      coeke(:,3,1) = coeke(:,1,3) ! symmetry of coeke
      !==========================================================================
      !
      CALL setclm(lpole, clm)
      !
      CALL getwc(wc, vwc, vnbsp, nbsp, h(1,1), h(2,2), h(3,3))
      !
      RETURN
      !
  END SUBROUTINE exx_setup_nscf
  !--------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE getwc(wc, vwc, vnbsp, nbsp, a1, a2, a3)
      ! it is called by the exx_setup_nscf (Broken)
      !
      IMPLICIT NONE
      !
      INTEGER     vnbsp, nbsp, ir
      REAl(DP)    wc(3, nbsp), vwc(3, vnbsp), a1, a2, a3
      !
      do ir = 1, nbsp
        read(407, *) wc(1,ir), wc(2,ir), wc(3,ir)
      end do
      !
      do ir = 1, vnbsp
        read(408,*)vwc(1,ir), vwc(2,ir), vwc(3,ir)
      enddo
      !
      do ir = 1, vnbsp
        if (vwc(1, ir) < 0) then
          vwc(1,ir) = vwc(1,ir) + a1
        end if
        if (vwc(2, ir) < 0) then
          vwc(2,ir) = vwc(2,ir) + a2
        end if
        if (vwc(3, ir) < 0) then
          vwc(3,ir) = vwc(3,ir) + a3
        end if
      end do
      !
      do ir = 1, nbsp
        if (wc(1, ir) < 0) then
          wc(1,ir) = wc(1,ir) + a1
        end if
        if (wc(2, ir) < 0) then
          wc(2,ir) = wc(2,ir) + a2
        end if
        if (wc(3, ir) < 0) then
          wc(3,ir) = wc(3,ir) + a3
        end if
      end do
      !
      RETURN
  END SUBROUTINE getwc
  !--------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------
  ! ==================================================================
  ! SUBROUTINE to set the various clm coefficients. Separated so as to
  ! not clutter up the above code. calculates clm = (l-m)!/(l+m)! for 
  ! l = 0,lpole, m=0,l. 
  ! evaluating these to 20 decimal places is certainly as accurate 
  ! as doing the calculations in the code but less work.
  ! ==================================================================
  SUBROUTINE setclm(lpole, clm)
      !
      IMPLICIT NONE
      ! INPUT: order of mutipole expansion
      INTEGER lpole
      ! OUTPUT: clm coefficients
      REAL*8 clm(0:lpole, 0:lpole)
      !
      clm(0,0) = 1.00000000000000000000e+00
      IF (lpole .GE. 1) THEN
        clm(1,0) = 1.00000000000000000000e+00
        clm(1,1) = 5.00000000000000000000e-01
      END IF
      IF (lpole .GE. 2) THEN
        clm(2,0) = 1.00000000000000000000e+00
        clm(2,1) = 1.66666666666666666670e-01
        clm(2,2) = 4.16666666666666666670e-02
      END IF
      IF (lpole .GE. 3) THEN
        clm(3,0) = 1.00000000000000000000e+00
        clm(3,1) = 8.33333333333333333330e-02
        clm(3,2) = 8.33333333333333333330e-03
        clm(3,3) = 1.38888888888888888890e-03
      END IF
      IF (lpole .GE. 4) THEN
        clm(4,0) = 1.00000000000000000000e+00
        clm(4,1) = 5.00000000000000000000e-02
        clm(4,2) = 2.77777777777777777780e-03
        clm(4,3) = 1.98412698412698412700e-04
        clm(4,4) = 2.48015873015873015870e-05
      END IF
      IF (lpole .GE. 5) THEN
        clm(5,0) = 1.00000000000000000000e+00
        clm(5,1) = 3.33333333333333333330e-02
        clm(5,2) = 1.19047619047619047620e-03
        clm(5,3) = 4.96031746031746031750e-05
        clm(5,4) = 2.75573192239858906530e-06
        clm(5,5) = 2.75573192239858906530e-07
      END IF
      IF (lpole .GE. 6) THEN
        clm(6,0) = 1.00000000000000000000e+00
        clm(6,1) = 2.38095238095238095240e-02
        clm(6,2) = 5.95238095238095238100e-04
        clm(6,3) = 1.65343915343915343920e-05
        clm(6,4) = 5.51146384479717813050e-07
        clm(6,5) = 2.50521083854417187750e-08
        clm(6,6) = 2.08767569878680989790e-09
      END IF
      IF (lpole .GE. 7) THEN
        clm(7,0) = 1.00000000000000000000e+00
        clm(7,1) = 1.78571428571428571430e-02
        clm(7,2) = 3.30687830687830687830e-04
        clm(7,3) = 6.61375661375661375660e-06
        clm(7,4) = 1.50312650312650312650e-07
        clm(7,5) = 4.17535139757361979580e-09
        clm(7,6) = 1.60590438368216145990e-10
        clm(7,7) = 1.14707455977297247140e-11
      END IF
      IF (lpole .GE. 8) THEN
        clm(8,0) = 1.00000000000000000000e+00
        clm(8,1) = 1.38888888888888888890e-02
        clm(8,2) = 1.98412698412698412700e-04
        clm(8,3) = 3.00625300625300625300e-06
        clm(8,4) = 5.01042167708834375500e-08
        clm(8,5) = 9.63542630209296875960e-10
        clm(8,6) = 2.29414911954594494280e-11
        clm(8,7) = 7.64716373181981647590e-13
        clm(8,8) = 4.77947733238738529740e-14
      END IF
      IF (lpole .GE. 9) THEN
        clm(9,0) = 1.00000000000000000000e+00
        clm(9,1) = 1.11111111111111111110e-02
        clm(9,2) = 1.26262626262626262630e-04
        clm(9,3) = 1.50312650312650312650e-06
        clm(9,4) = 1.92708526041859375190e-08
        clm(9,5) = 2.75297894345513393130e-10
        clm(9,6) = 4.58829823909188988550e-12
        clm(9,7) = 9.55895466477477059490e-14
        clm(9,8) = 2.81145725434552076320e-15
        clm(9,9) = 1.56192069685862264620e-16
      END IF
      !
      RETURN
  END SUBROUTINE setclm  
  !--------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------
  !===============================================================
  !
  !Coefficients for the first & second order numerical derivative
  !under the centered finite difference scheme.
  !Bengt Fornberg,  Exxon Res. & Eng. Co., NJ 08801
  !'bfornbe@erenj.com'
  !David M. Sloan,  Dept. of Mathematics, U. of Strathclyde,
  !Glasgow G1 1HX, Scotland,  'd.sloan@strath.ac.uk'
  !Acta Numerica 94,  Cambridge Univ. Press (1994)
  !
  !===============================================================
  SUBROUTINE fornberg(norder1,norder2,coe1,coe,coe_cross,ierr)
      !
      !     add functionality for the general cell
      !     we follow B. Fornberg in  Math. Comp. 51 (1988), 699-706
      !
      USE kinds, ONLY : DP
      !
      IMPLICIT NONE
      !
      !     Input/Output variables:
      !
      !     order of expansion of derivative. 
      !     it is the number of neighbors used ON ONE SIDE.
      !     the maximum order implemented is 20.
      INTEGER, INTENT(IN) :: norder1, norder2
      !
      !  coe1 - coefficients for the first derivative
      REAL(DP), INTENT(OUT) :: coe1(-norder1:norder1)
      !
      !  coe  - coefficients for the axial derivative (2 nd)
      REAL(DP), INTENT(OUT) :: coe(-norder2:norder2)
      !
      !  coe_cross - coefficients for the cross derivative (2 nd)
      REAL(DP), INTENT(OUT) :: coe_cross(-norder2:norder2)
      !
      !     ierr - error flag
      INTEGER, INTENT(OUT) :: ierr
      !     
      !     Work variables:
      !
      !     counters 
      INTEGER i
      !     ---------------------------------------------------------------
      !
      !     First order derivative
      !
      ierr = 0
      !
      SELECT CASE (norder1)
      CASE (1)
        coe1(1) =  0.50000000000000D+00
      CASE (2)
        coe1(1) =  2.d0/3.d0
        coe1(2) = -1.d0/12.d0
      CASE (3)
        coe1(1) =  3.d0/4.d0
        coe1(2) = -3.d0/20.d0
        coe1(3) =  1.d0/60.d0
      CASE (4)
        coe1(1) =  4.d0/5.d0
        coe1(2) = -1.d0/5.d0
        coe1(3) =  4.d0/105.d0
        coe1(4) = -1.d0/280.d0
      CASE (5)
        coe1(1) =  0.8333333333D+00
        coe1(2) = -0.2380952381D+00
        coe1(3) =  0.5952380952D-01
        coe1(4) = -0.9920634921D-02
        coe1(5) =  0.7936507937D-03
      CASE (6)
        coe1(1) =  0.8571428571D+00
        coe1(2) = -0.2678571429D+00
        coe1(3) =  0.7936507937D-01
        coe1(4) = -0.1785714286D-01
        coe1(5) =  0.2597402597D-02
        coe1(6) = -0.1803751804D-03
      CASE (7)
        coe1(1) =  0.8750000000D+00
        coe1(2) = -0.2916666667D+00
        coe1(3) =  0.9722222222D-01
        coe1(4) = -0.2651515152D-01
        coe1(5) =  0.5303030303D-02
        coe1(6) = -0.6798756799D-03
        coe1(7) =  0.4162504163D-04
      CASE (8)
        coe1(1) =  0.8888888889D+00
        coe1(2) = -0.3111111111D+00
        coe1(3) =  0.1131313131D+00
        coe1(4) = -0.3535353535D-01
        coe1(5) =  0.8702408702D-02
        coe1(6) = -0.1554001554D-02
        coe1(7) =  0.1776001776D-03
        coe1(8) = -0.9712509713D-05
      CASE (9)
        coe1(1) =  0.9000000000D+00
        coe1(2) = -0.3272727273D+00
        coe1(3) =  0.1272727273D+00
        coe1(4) = -0.4405594406D-01
        coe1(5) =  0.1258741259D-01
        coe1(6) = -0.2797202797D-02
        coe1(7) =  0.4495504496D-03
        coe1(8) = -0.4627725216D-04
        coe1(9) =  0.2285296403D-05
      CASE (10)
        coe1(1) =  0.9090909091D+00
        coe1(2) = -0.3409090909D+00
        coe1(3) =  0.1398601399D+00
        coe1(4) = -0.5244755245D-01
        coe1(5) =  0.1678321678D-01
        coe1(6) = -0.4370629371D-02
        coe1(7) =  0.8814714697D-03
        coe1(8) = -0.1285479227D-03
        coe1(9) =  0.1202787580D-04
        coe1(10)= -0.5412544112D-06
      END SELECT
      !
      coe1(0) = 0.0_DP
      DO i = 1,norder1
        coe1(-i) = -coe1(i)
      END DO
      !
      !     Second order derivative
      !
      coe_cross(0) = 0.0_DP
      !
      SELECT CASE (norder2)
      CASE (1)
        coe(0) = -0.20000000000000D+01
        coe(1) =  0.10000000000000D+01
        ! cross derivative
        coe_cross(1) =  0.25_DP
      CASE (2)
        coe(0) = -0.25000000000000D+01
        coe(1) =  0.13333333333333D+01
        coe(2) = -0.83333333333333D-01
        ! cross derivative
        coe_cross(1) =  1.0_DP/3.0_DP
        coe_cross(2) = -1.0_DP/48.0_DP
      CASE (3)
        coe(0) = -0.27222222222222D+01
        coe(1) =  0.15000000000000D+01
        coe(2) = -0.15000000000000D+00
        coe(3) =  0.11111111111111D-01
        ! cross derivative
        coe_cross(1) =  3.0_DP/8.0_DP
        coe_cross(2) = -3.0_DP/80.0_DP
        coe_cross(3) =  1.0_DP/360.0_DP
      CASE (4)
        coe(0) = -0.28472222222222D+01
        coe(1) =  0.16000000000000D+01
        coe(2) = -0.20000000000000D+00
        coe(3) =  0.25396825396825D-01
        coe(4) = -0.17857142857143D-02
        ! cross derivative
        CALL errore('fornberg','cross derivative not yet &
            &implemented for this nord2', nord2)
      CASE (5)
        coe(0) = -0.29272222222222D+01
        coe(1) =  0.16666666666667D+01
        coe(2) = -0.23809523809524D+00
        coe(3) =  0.39682539682540D-01
        coe(4) = -0.49603174603175D-02
        coe(5) =  0.31746031746032D-03
        ! cross derivative
        CALL errore('fornberg','cross derivative not yet &
            &implemented for this nord2', nord2)
      CASE (6)
        coe(0) = -0.29827777777778D+01
        coe(1) =  0.17142857142857D+01
        coe(2) = -0.26785714285714D+00
        coe(3) =  0.52910052910053D-01
        coe(4) = -0.89285714285714D-02
        coe(5) =  0.10389610389610D-02
        coe(6) = -0.60125060125060D-04
        ! cross derivative
        CALL errore('fornberg','cross derivative not yet &
            &implemented for this nord2', nord2)
      CASE (7)
        coe(0) = -0.30235941043084D+01
        coe(1) =  0.17500000000000D+01
        coe(2) = -0.29166666666667D+00
        coe(3) =  0.64814814814815D-01
        coe(4) = -0.13257575757576D-01
        coe(5) =  0.21212121212121D-02
        coe(6) = -0.22662522662523D-03
        coe(7) =  0.11892869035726D-04
        ! cross derivative
        CALL errore('fornberg','cross derivative not yet &
            &implemented for this nord2', nord2)
      CASE (8)
        coe(0) = -0.30548441043084D+01
        coe(1) =  0.17777777777778D+01
        coe(2) = -0.31111111111111D+00
        coe(3) =  0.75420875420875D-01
        coe(4) = -0.17676767676768D-01
        coe(5) =  0.34809634809635D-02
        coe(6) = -0.51800051800052D-03
        coe(7) =  0.50742907885765D-04
        coe(8) = -0.24281274281274D-05
        ! cross derivative
        CALL errore('fornberg','cross derivative not yet &
            &implemented for this nord2', nord2)
      CASE (9)
        coe(0) = -0.30795354623331D+01
        coe(1) =  0.18000000000000D+01
        coe(2) = -0.32727272727273D+00
        coe(3) =  0.84848484848485D-01
        coe(4) = -0.22027972027972D-01
        coe(5) =  0.50349650349650D-02
        coe(6) = -0.93240093240093D-03
        coe(7) =  0.12844298558584D-03
        coe(8) = -0.11569313039901D-04
        coe(9) =  0.50784364509855D-06
        ! cross derivative
        CALL errore('fornberg','cross derivative not yet &
            &implemented for this nord2', nord2)
      CASE (10)
        coe(0) = -0.30995354623331D+01
        coe(1) =  0.18181818181818D+01
        coe(2) = -0.34090909090909D+00
        coe(3) =  0.93240093240093D-01
        coe(4) = -0.26223776223776D-01
        coe(5) =  0.67132867132867D-02
        coe(6) = -0.14568764568765D-02
        coe(7) =  0.25184899134479D-03
        coe(8) = -0.32136980666392D-04
        coe(9) =  0.26728612899924D-05
        coe(10)= -0.10825088224469D-06
        ! cross derivative
        CALL errore('fornberg','cross derivative not yet &
            &implemented for this nord2', nord2)
      END SELECT
      !
      DO i = 1,norder2
        coe(-i) = coe(i)
        coe_cross(-i) = coe_cross(i)
      END DO
      !
      RETURN
  END SUBROUTINE fornberg
  !--------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE exx_pot_derivative (np_in_sp_me, n, v, dvdr)
      !
      !  To calculate the cell derivatives for the EXX calculation, one needs
      !  the derivatives of the exchange potential.
      !
      IMPLICIT NONE
      ! io
      INTEGER , INTENT(IN)  :: np_in_sp_me, n
      REAL(DP), INTENT(IN)  :: v(np_in_sp_me)
      REAL(DP), INTENT(OUT) :: dvdr(n,3)
      ! loc
      INTEGER               :: i, ish, ii, jj, kk
      REAL(DP)              :: v1, v2, v3
      !
      !  --------------------------------------------------------
      !  initialize
      dvdr(:,:) = 0.0_DP
      !
      !  derivatives
      !
      !$omp parallel do private(ii,jj,kk,v1,v2,v3)
      DO i = 1, n ! running inside PS sphere
        !
        ii = odtothd_in_sp(1,i)
        jj = odtothd_in_sp(2,i)
        kk = odtothd_in_sp(3,i)
        !
        DO ish = 1, nord1
          ! since the derivative coeff is odd in parity,
          !     we combine neighbors
          v1 = v( thdtood_in_sp( ii+ish, jj,     kk))    - &
              v( thdtood_in_sp( ii-ish, jj,     kk))
          !
          v2 = v( thdtood_in_sp( ii,     jj+ish, kk))    - &
              v( thdtood_in_sp( ii,     jj-ish, kk))
          !
          v3 = v( thdtood_in_sp( ii,     jj,     kk+ish))- &
              v( thdtood_in_sp( ii,     jj,     kk-ish))
          !
          dvdr(i,1) = dvdr(i,1) + coe_1st_derv(ish,1)*v1
          dvdr(i,2) = dvdr(i,2) + coe_1st_derv(ish,2)*v2
          dvdr(i,3) = dvdr(i,3) + coe_1st_derv(ish,3)*v3
          !
        END DO
        !
      END DO
      !$omp end parallel do
      !
      RETURN
      !
  END SUBROUTINE exx_pot_derivative
  !--------------------------------------------------------------------------------------------------------------

  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE exx_energy_cell_derivative (np_in_sp_me, n, tran, &
          vl, ha, hb, hc, rhol, dexx_dhab)
      !
      IMPLICIT NONE
      ! io
      INTEGER , INTENT(IN)  :: np_in_sp_me, n, tran(3)
      REAL(DP), INTENT(IN)  :: vl(np_in_sp_me), rhol(np_in_sp_me) ! rhol pass in side n is good enough
      REAL(DP), INTENT(IN)  :: ha, hb, hc
      REAL(DP), INTENT(OUT) :: dexx_dhab(3,3)
      ! loc
      INTEGER               :: i, ish, ii, jj, kk, alpha, beta
      REAL(DP)              :: dvdr(n,3), Imn(3,3), r_alpha(3)
      REAL(DP)              :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
      !
      ! initialize
      !
      tmp1 = 0.0_DP; tmp2 = 0.0_DP; tmp3 = 0.0_DP
      tmp4 = 0.0_DP; tmp5 = 0.0_DP; tmp6 = 0.0_DP
      !
      Imn(:,:) = 0.0_DP
      !
      ! Calculate the potential derivative
      !
      CALL exx_pot_derivative (np_in_sp_me, n, vl, dvdr)
      !
      ! Integration:
      !             I_{mn} = \int \dd r \rho(r) r_m v_n (r)
      !
      !$omp parallel do private(i,ii,jj,kk,r_alpha) reduction(+:tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
      DO i = 1, n ! running inside PS sphere
        !
        ! Calculate r_m on the fly
        !
        ii = odtothd_in_sp(1,i) - tran(1) ! include the l2goff (1d->3d)
        jj = odtothd_in_sp(2,i) - tran(2) ! include the l2goff (1d->3d)
        kk = odtothd_in_sp(3,i) - tran(3) ! include the l2goff (1d->3d)
        !
        r_alpha(1) = DBLE(ii)*ha
        r_alpha(2) = DBLE(jj)*hb
        r_alpha(3) = DBLE(kk)*hc
        !
        ! Combining potential derivatives with cell parameters
        !
        tmp1 = tmp1 + rhol(i)*r_alpha(1)*dvdr(i,1)
        tmp2 = tmp2 + rhol(i)*r_alpha(1)*dvdr(i,2)
        tmp3 = tmp3 + rhol(i)*r_alpha(1)*dvdr(i,3)
        tmp4 = tmp4 + rhol(i)*r_alpha(2)*dvdr(i,2)
        tmp5 = tmp5 + rhol(i)*r_alpha(2)*dvdr(i,3)
        tmp6 = tmp6 + rhol(i)*r_alpha(3)*dvdr(i,3)
        !
      END DO
      !$omp end parallel do
      !
      Imn(1,1) = tmp1; Imn(1,2) = tmp2; Imn(1,3) = tmp3
      Imn(2,1) = tmp2; Imn(2,2) = tmp4; Imn(2,3) = tmp5
      Imn(3,1) = tmp3; Imn(3,2) = tmp5; Imn(3,3) = tmp6
      !
      ! Calculate the stress tensor from Imn
      !   (d exx / d hab )_{ij} =  2 (sum_k  Imn_{ik}(ainv)_{jk})
      DO alpha = 1, 3
        DO beta = 1, 3
          dexx_dhab(alpha,beta) = Imn(alpha,1)*ainv(beta,1)+&
              Imn(alpha,2)*ainv(beta,2)+Imn(alpha,3)*ainv(beta,3)
        END DO
      END DO
      !
      dexx_dhab = 2.0_DP*dexx_dhab
      !
      RETURN
  END SUBROUTINE exx_energy_cell_derivative

  !--------------------------------------------------------------------------------------------------------------
  SUBROUTINE exx_energy_cell_derivative_cube(me_r, ps_r, tran, rho, pot, &
                                      ha_proj, hb_proj, hc_proj, Jim, dexx_dhab)
    !
    IMPLICIT NONE
    !----------------------------------------------------------------------------
    ! pass in variable
    !----------------------------------------------------------------------------
    INTEGER , INTENT(IN)  :: me_r(6)
    INTEGER , INTENT(IN)  :: ps_r(6)
    INTEGER , INTENT(IN)  :: tran(3)
    REAL(DP), INTENT(IN)  :: rho(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP), INTENT(IN)  :: pot(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP), INTENT(IN)  :: ha_proj(3)
    REAL(DP), INTENT(IN)  :: hb_proj(3)
    REAL(DP), INTENT(IN)  :: hc_proj(3)
    REAL(DP), INTENT(IN)  :: Jim(3,3)  ! jacobian [d/d x]        [d/d a]
    !                                             |d/d y| = [J]  |d/d b|
    !                                             [d/d z]        [d/d c]
    REAL(DP), INTENT(OUT) :: dexx_dhab(3,3)
    !----------------------------------------------------------------------------
    ! local variable
    !----------------------------------------------------------------------------
    INTEGER               :: i,j,k,ish
    INTEGER               :: ii,jj,kk
    INTEGER               :: alpha, beta
    REAL(DP)              :: Imn(3,3), r_alpha(3)
    REAL(DP)              :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
    !
    REAL(DP) :: dvdri1, dvdri2, dvdri3 
    REAL(DP) :: dvdr1, dvdr2, dvdr3 
    integer              :: tran1, tran2, tran3
    REAL(DP)              :: r_alpha1, r_alpha2, r_alpha3
    !REAL(DP), ALLOCATABLE :: dvdr(:,:,:,:)
#ifdef __CUDA
    REAL(DP), device  :: ha_proj_d(3)
    REAL(DP), device  :: hb_proj_d(3)
    REAL(DP), device  :: hc_proj_d(3)
    REAL(DP), device  :: Jim_d(3,3)
    attributes(device) :: pot, rho
#endif
    !----------------------------------------------------------------------------------
    
    !----------------------------------------------------------------------------------
    ! initialize
    !----------------------------------------------------------------------------------
    !ALLOCATE(dvdr(3,ps_r(1):ps_r(4),ps_r(2):ps_r(5),ps_r(3):ps_r(6))); dvdr=0.0D0
    !----------------------------------------------------------------------------------
    tmp1 = 0.0_DP; tmp2 = 0.0_DP; tmp3 = 0.0_DP
    tmp4 = 0.0_DP; tmp5 = 0.0_DP; tmp6 = 0.0_DP
    !----------------------------------------------------------------------------------
    Imn(:,:) = 0.0_DP
#ifdef __CUDA
    ha_proj_d = ha_proj
    hb_proj_d = hb_proj
    hc_proj_d = hc_proj
    Jim_d = Jim
#endif
    !----------------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------------
    ! Calculate the potential derivative
    !----------------------------------------------------------------------------------
    !CALL exx_pot_derivative(me_r, ps_r, pot, Jim, dvdr)
    !----------------------------------------------------------------------------------

    tran1 = tran(1)
    tran2 = tran(2)
    tran3 = tran(3)
    
    !----------------------------------------------------------------------------------
    ! Integration: I_{mn} = \int \dd r \rho(r) r_m dv_n (r)
    !----------------------------------------------------------------------------------
#ifdef __CUDA
    associate (ha_proj=>ha_proj_d, hb_proj=>hb_proj_d, hc_proj=>hc_proj_d, Jim=>Jim_d, &
    &      coe_1st_derv=>coe_1st_derv_d )
    !$cuf kernel do (3)
#else
   !$omp parallel do collapse(3) reduction(+:tmp1,tmp2,tmp3,tmp4,tmp5,tmp6) &
   !$omp private(i,j,k,ii,jj,kk,ish) &
   !$omp private(r_alpha1, r_alpha2, r_alpha3) &
   !$omp private(dvdri1, dvdri2, dvdri3, dvdr1, dvdr2, dvdr3) &
   !$omp firstprivate(ha_proj, hb_proj, hc_proj) &
   !$omp firstprivate(Jim,tran1,tran2,tran3)
#endif
    DO k = ps_r(3),ps_r(6)
      DO j = ps_r(2),ps_r(5)
        DO i = ps_r(1),ps_r(4)
          !----------------------------------------------------
          ! Calculate r_m on the fly
          !----------------------------------------------------
          ii = i - tran1 ! include the l2goff (1d->3d)
          jj = j - tran2 ! include the l2goff (1d->3d)
          kk = k - tran3 ! include the l2goff (1d->3d)
          !----------------------------------------------------
          r_alpha1 = DBLE(ii)*ha_proj(1)+DBLE(jj)*hb_proj(1)+DBLE(kk)*hc_proj(1)
          r_alpha2 = DBLE(ii)*ha_proj(2)+DBLE(jj)*hb_proj(2)+DBLE(kk)*hc_proj(2)
          r_alpha3 = DBLE(ii)*ha_proj(3)+DBLE(jj)*hb_proj(3)+DBLE(kk)*hc_proj(3)
          !----------------------------------------------------
          dvdri1 = 0.d0; dvdri2 = 0.d0; dvdri3 = 0.d0
          DO ish = 1, nord1
            ! since the derivative coeff is odd in parity, we combine neighbors
            dvdri1 = dvdri1 + coe_1st_derv(ish,1)*(pot(i+ish,j,k)-pot(i-ish,j,k))
            dvdri2 = dvdri2 + coe_1st_derv(ish,2)*(pot(i,j+ish,k)-pot(i,j-ish,k))
            dvdri3 = dvdri3 + coe_1st_derv(ish,3)*(pot(i,j,k+ish)-pot(i,j,k-ish))
          END DO
          !----------------------------------------------------
          dvdr1 = Jim(1,1)*dvdri1 + Jim(1,2)*dvdri2 + Jim(1,3)*dvdri3
          dvdr2 = Jim(2,1)*dvdri1 + Jim(2,2)*dvdri2 + Jim(2,3)*dvdri3
          dvdr3 = Jim(3,1)*dvdri1 + Jim(3,2)*dvdri2 + Jim(3,3)*dvdri3
          !----------------------------------------------------
          !----------------------------------------------------
          ! Combining potential derivatives with cell parameters
          !----------------------------------------------------
          tmp1 = tmp1 + rho(i,j,k)*r_alpha1*dvdr1
          tmp2 = tmp2 + rho(i,j,k)*r_alpha1*dvdr2
          tmp3 = tmp3 + rho(i,j,k)*r_alpha1*dvdr3
          tmp4 = tmp4 + rho(i,j,k)*r_alpha2*dvdr2
          tmp5 = tmp5 + rho(i,j,k)*r_alpha2*dvdr3
          tmp6 = tmp6 + rho(i,j,k)*r_alpha3*dvdr3
          !----------------------------------------------------
        END DO
      END DO
    END DO
#ifdef __CUDA
    end associate
#else
    !$omp end parallel do
#endif
    !----------------------------------------------------------------------------------
    Imn(1,1) = tmp1; Imn(1,2) = tmp2; Imn(1,3) = tmp3
    Imn(2,1) = tmp2; Imn(2,2) = tmp4; Imn(2,3) = tmp5
    Imn(3,1) = tmp3; Imn(3,2) = tmp5; Imn(3,3) = tmp6
    !----------------------------------------------------------------------------------
    ! Calculate the stress tensor from Imn: (d exx / d hab )_{ij} =  2 (sum_k  Imn_{ik}(ainv)_{jk})
    !----------------------------------------------------------------------------------
    DO alpha = 1, 3
      DO beta = 1, 3
        dexx_dhab(alpha,beta) = Imn(alpha,1)*ainv(beta,1)+Imn(alpha,2)*ainv(beta,2)+Imn(alpha,3)*ainv(beta,3)
      END DO
    END DO
    !----------------------------------------------------------------------------------
    dexx_dhab = 2.0_DP*dexx_dhab
    !----------------------------------------------------------------------------------
    !
    !DEALLOCATE(dvdr)
    !
    RETURN
  END SUBROUTINE exx_energy_cell_derivative_cube
  !--------------------------------------------------------------------------------------------------------------

#ifdef __CUDA
  attributes(host,device) &
#endif
  integer function l2gcb(n,l,t)
    !! This function is the cube analogue of the \(\texttt{l2goff}\) function in
    !! \(\texttt{exx_gs.f90}\).  
    !! These functions provides a local to global grid transformation to allow
    !! sparse matrix/vector operations.
    implicit none
    integer, value :: n, l, t
    l2gcb = MOD(l-t-1+n, n)+1
  end function l2gcb
END MODULE exx_module
