!
! Copyright (C) 2002-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(HAVE_GITREV)
#include "git-rev.h"
#endif
!
!==-----------------------------------------------------------------------==!
MODULE environment
  !==-----------------------------------------------------------------------==!
  !! Environment management.
  USE kinds, ONLY: DP
  USE io_files, ONLY: crash_file, nd_nmbr
  USE io_global, ONLY: stdout, meta_ionode
  USE mp_world,  ONLY: nproc, nnode
  USE mp_images, ONLY: me_image, my_image_id, root_image, nimage, &
      nproc_image
  USE mp_pools,  ONLY: npool
  USE mp_bands,  ONLY: ntask_groups, nproc_bgrp, nbgrp, nyfft
  USE global_version, ONLY: version_number
  USE fox_init_module, ONLY: fox_init
  USE command_line_options, ONLY : nmany_
  USE clib_wrappers, ONLY : get_mem_avail
#if defined(__HDF5)
  USE qeh5_base_module,   ONLY: initialize_hdf5, finalize_hdf5
#endif

  IMPLICIT NONE

  ! ... code name
  CHARACTER(LEN=20) :: code = 'notset'

  SAVE

  PRIVATE

  PUBLIC :: environment_start
  PUBLIC :: environment_end
  PUBLIC :: opening_message
  PUBLIC :: parallel_info
  PUBLIC :: print_cuda_info

  !==-----------------------------------------------------------------------==!
CONTAINS
  !==-----------------------------------------------------------------------==!

  SUBROUTINE environment_start( code_in )

    CHARACTER(LEN=*), INTENT(IN) :: code_in

    LOGICAL           :: exst, debug = .false.
    CHARACTER(LEN=80) :: code_version, uname
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    CHARACTER(LEN=3)           :: env_maxdepth
    INTEGER :: ios, crashunit, max_depth 

#if defined(__TRACE)
    CALL get_environment_variable('ESPRESSO_MAX_DEPTH', env_maxdepth)
    IF (env_maxdepth .NE. ' ') THEN 
      READ(env_maxdepth,'(I3)',iostat=ios) max_depth
      IF (ios == 0 ) CALL set_trace_max_depth( max_depth )
    END IF
#endif
    ! ... use ".FALSE." to disable all clocks except the total cpu time clock
    ! ... use ".TRUE."  to enable clocks
    CALL init_clocks(.TRUE.) 
    code = code_in
    CALL start_clock( TRIM(code) )

    code_version = TRIM (code) // " v." // TRIM (version_number)

    ! ... for compatibility with PWSCF
#if defined(__MPI)
    nd_nmbr = TRIM ( int_to_char( me_image+1 ))
#else
    nd_nmbr = ' '
#endif

    IF( meta_ionode ) THEN

       ! ...  search for file CRASH and delete it

       INQUIRE( FILE=TRIM(crash_file), EXIST=exst )
       IF( exst ) THEN
          OPEN( NEWUNIT=crashunit, FILE=TRIM(crash_file), STATUS='OLD',IOSTAT=ios )
          IF (ios==0) THEN
             CLOSE( UNIT=crashunit, STATUS='DELETE', IOSTAT=ios )
          ELSE
             WRITE(stdout,'(5x,"Remark: CRASH file could not be deleted")')
          END IF
       END IF

    ELSE
       ! ... one processor per image (other than meta_ionode)
       ! ... or, for debugging purposes, all processors,
       ! ... open their own standard output file
!#define DEBUG
#if defined(DEBUG)
       debug = .true.
#endif
       IF (me_image == root_image .OR. debug ) THEN
          uname = 'out.' // trim(int_to_char( my_image_id )) // '_' // &
               trim(int_to_char( me_image))
          OPEN ( unit = stdout, file = TRIM(uname),status='unknown')
       ELSE
#if defined(_WIN32)
          OPEN ( unit = stdout, file='NUL:', status='unknown' )
#else
          OPEN ( unit = stdout, file='/dev/null', status='unknown' )
#endif
       END IF

    END IF
    !
    CALL opening_message( code_version )
#if defined(__MPI)
    CALL parallel_info ( )
#else
    CALL serial_info()
#endif
    CALL fox_init()
#if defined(__HDF5)
    CALL initialize_hdf5()
#endif
    !
    WRITE(stdout,'(5x, I0, A, A)') get_mem_avail()/1024, &
                &" MiB available memory on the printing compute node ", &
                &"when the environment starts"
    WRITE(stdout, *)
  END SUBROUTINE environment_start

  !==-----------------------------------------------------------------------==!

  SUBROUTINE environment_end( code_in )
    ! next line for back-compatibility
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: code_in
#if defined(_HDF5)
    CALL finalize_hdf5()
#endif
    IF ( meta_ionode ) WRITE( stdout, * )

    IF ( PRESENT(code_in) ) THEN
       code = code_in
    END IF
    IF ( code == 'notset' ) THEN
       WRITE( stdout,'(5X,"WARNING: environment_end needs a call to environment_start at the beginning of the run")')
    ELSE
       CALL stop_clock(  TRIM(code) )
       CALL print_clock( TRIM(code) )
    END IF
    CALL closing_message( )

    IF( meta_ionode ) THEN
       WRITE( stdout,'(A)')      '   JOB DONE.'
       WRITE( stdout,3335)
    END IF
3335 FORMAT('=',78('-'),'=')
    FLUSH(stdout)

    RETURN
  END SUBROUTINE environment_end

  !==-----------------------------------------------------------------------==!

  SUBROUTINE opening_message( code_version )

    CHARACTER(LEN=*), INTENT(IN) :: code_version
    CHARACTER(LEN=9)  :: cdate, ctime

    CALL date_and_tim( cdate, ctime )
    !
    WRITE( stdout, '(/5X,"Program ",A," starts on ",A9," at ",A9)' ) &
         TRIM(code_version), cdate, ctime
#if defined(HAVE_GITREV)
    WRITE( stdout, '(8X, "Git branch: ", A)' ) &
      GIT_BRANCH_RAW
    WRITE( stdout, '(8X, "Last git commit: ", A)' ) &
      GIT_HASH_RAW
    WRITE( stdout, '(8X, "Last git commit date: ", A)' ) & 
      GIT_COMMIT_LAST_CHANGED_RAW
    WRITE( stdout, '(8X, "Last git commit subject: ", A)' ) & 
      GIT_COMMIT_SUBJECT_RAW
#endif
    !
    WRITE( stdout, '(/5X,"This program is part of the open-source Quantum ",&
         &    "ESPRESSO suite", &
         &/5X,"for quantum simulation of materials; please cite",   &
         &/9X,"""P. Giannozzi et al., J. Phys.:Condens. Matter 21 ",&
         &    "395502 (2009);", &
         &/9X,"""P. Giannozzi et al., J. Phys.:Condens. Matter 29 ",&
         &    "465901 (2017);", &
         &/9X,"""P. Giannozzi et al., J. Chem. Phys. 152 ",&
         &    "154105 (2020);", &
         &/9X," URL http://www.quantum-espresso.org"", ", &
         &/5X,"in publications or presentations arising from this work. More details at",&
         &/5x,"http://www.quantum-espresso.org/quote")' )
    RETURN
  END SUBROUTINE opening_message

  !==-----------------------------------------------------------------------==!

  SUBROUTINE closing_message( )

    CHARACTER(LEN=9)  :: cdate, ctime
    CHARACTER(LEN=80) :: time_str

    CALL date_and_tim( cdate, ctime )

    time_str = 'This run was terminated on:  ' // ctime // ' ' // cdate

    IF( meta_ionode ) THEN
       WRITE( stdout,*)
       WRITE( stdout,3334) time_str
       WRITE( stdout,3335)
    END IF

3334 FORMAT(3X,A60,/)
3335 FORMAT('=',78('-'),'=')

    RETURN
  END SUBROUTINE closing_message

  !==-----------------------------------------------------------------------==!
  SUBROUTINE parallel_info ( )
    !
#if defined(_OPENMP)
    INTEGER, EXTERNAL :: omp_get_max_threads
    !
    WRITE( stdout, '(/5X,"Parallel version (MPI & OpenMP), running on ",&
         &I7," processor cores")' ) nproc * omp_get_max_threads()
    !
    WRITE( stdout, '(5X,"Number of MPI processes:           ",I7)' ) nproc
    !
    WRITE( stdout, '(5X,"Threads/MPI process:               ",I7)' ) &
         omp_get_max_threads()
#else
    WRITE( stdout, '(/5X,"Parallel version (MPI), running on ",&
         &I5," processors")' ) nproc
#endif
    !
    WRITE( stdout, '(/5X,"MPI processes distributed on ",&
         &I5," nodes")' ) nnode
    IF ( nimage > 1 ) WRITE( stdout, &
         '(5X,"path-images division:  nimage    = ",I7)' ) nimage
    IF ( npool > 1 ) WRITE( stdout, &
         '(5X,"K-points division:     npool     = ",I7)' ) npool
    IF ( nbgrp > 1 ) WRITE( stdout, &
         '(5X,"band groups division:  nbgrp     = ",I7)' ) nbgrp
    IF ( nproc_bgrp > 1 ) WRITE( stdout, &
         '(5X,"R & G space division:  proc/nbgrp/npool/nimage = ",I7)' ) nproc_bgrp
    IF ( nyfft > 1 ) WRITE( stdout, &
         '(5X,"wavefunctions fft division:  Y-proc x Z-proc = ",2I7)' ) &
         nyfft, nproc_bgrp / nyfft
    IF ( ntask_groups > 1 ) WRITE( stdout, &
         '(5X,"wavefunctions fft division:  task group distribution",/,34X,"#TG    x Z-proc = ",2I7)' ) &
         ntask_groups, nproc_bgrp / ntask_groups
    IF ( nmany_ > 1) WRITE( stdout, '(5X,"FFT bands division:     nmany     = ",I7)' ) nmany_
    !
  END SUBROUTINE parallel_info

  !==-----------------------------------------------------------------------==!
  SUBROUTINE serial_info ( )
    !
#if defined(_OPENMP)
    INTEGER, EXTERNAL :: omp_get_max_threads
#endif
    !
#if defined(_OPENMP)
    WRITE( stdout, '(/5X,"Serial multi-threaded version, running on ",&
         &I4," processor cores")' ) omp_get_max_threads()
    !
#else
    WRITE( stdout, '(/5X,"Serial version")' )
#endif
    !
  END SUBROUTINE serial_info
!
!-----------------------------------------------------------------------
SUBROUTINE print_cuda_info(check_use_gpu) 
  !-----------------------------------------------------------------------
  !
  USE io_global,       ONLY : stdout
  USE control_flags,   ONLY : use_gpu_=> use_gpu, iverbosity
  USE mp_world,        ONLY : nnode, nproc
  USE mp,              ONLY : mp_sum, mp_max
#if defined(__CUDA)
  USE cudafor
#endif 
  !
  IMPLICIT NONE
  !
  LOGICAL, OPTIONAL,INTENT(IN)  :: check_use_gpu 
  !! if present and trues the internal variable use_gpu is checked
#if defined (__CUDA) 
  INTEGER :: idev, ndev, ierr
  TYPE (cudaDeviceProp) :: prop
  LOGICAL               :: use_gpu = .TRUE. 
  !
  IF ( PRESENT(check_use_gpu) ) THEN 
    IF (check_use_gpu) use_gpu = use_gpu_ 
  END IF 
  !
  ierr = cudaGetDevice( idev )
  IF (ierr /= 0) CALL errore('summary', 'cannot get device id', ierr)
  ierr = cudaGetDeviceCount( ndev )
  IF (ierr /= 0) CALL errore('summary', 'cannot get device count', ierr)
  !
  IF (use_gpu) THEN
     WRITE( stdout, '(/,5X,"GPU acceleration is ACTIVE. ",i2," visible GPUs per MPI rank")' ) ndev
#if defined(__GPU_MPI)
     WRITE( stdout, '(5x, "GPU-aware MPI enabled")')
#endif
     WRITE( stdout, '()' )
  ELSE
     WRITE( stdout, '(/,5X,"GPU acceleration is NOT ACTIVE.",/)' )
  END IF
  !
  ! User friendly, approximated warning. FIXME: does it really work?
  ! In order to get this done right, one needs an intra_node communicator
  !
  IF (nproc > ndev * nnode * 2) &
     CALL infomsg('print_cuda_info', &
      'High GPU oversubscription detected. Are you sure this is what you want?')
  !
  ! Verbose information for advanced users
  IF (iverbosity > 0) THEN
     WRITE( stdout, '(/,5X,"GPU used by master process:",/)' )
     ! Device info taken from
     ! https://devblogs.nvidia.com/how-query-device-properties-and-handle-errors-cuda-fortran/
     ierr = cudaGetDeviceProperties(prop, idev)
     WRITE(stdout,"(5X,'   Device Number: ',i0)") idev
     WRITE(stdout,"(5X,'   Device name: ',a)") trim(prop%name)
     WRITE(stdout,"(5X,'   Compute capability: ',i0, i0)") prop%major, prop%minor
     WRITE(stdout,"(5X,'   Ratio of single to double precision performance: ',i0)") prop%singleToDoublePrecisionPerfRatio
     WRITE(stdout,"(5X,'   Memory Clock Rate (KHz): ', i0)") &
       prop%memoryClockRate
     WRITE(stdout,"(5X,'   Memory Bus Width (bits): ', i0)") &
       prop%memoryBusWidth
     WRITE(stdout,"(5X,'   Peak Memory Bandwidth (GB/s): ', f8.2)") &
       2.0d0*prop%memoryClockRate*(prop%memoryBusWidth/8.0d0)/10.0**6
     WRITE(stdout,"(5X,'   Total Global Memory (MiB): ',i0)") prop%totalGlobalMem/2**20
     ! Memory reported in MiB, not in MB to match what "nvidia-smi" reports.
     ! 1 MiB = 2^20 bytes = 1048576 bytes
     ! 1 MB  = 10^6 bytes = 1000000 bytes
     WRITE(stdout,'()')
  END IF
  !
#endif
  !
END SUBROUTINE print_cuda_info

  !==-----------------------------------------------------------------------==!
END MODULE environment
!==-----------------------------------------------------------------------==!
