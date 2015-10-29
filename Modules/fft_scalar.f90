!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------!
! FFT scalar drivers Module - contains machine-dependent routines for      !
! FFTW, FFTW3, ESSL (both 3d for serial execution and 1d+2d FFTs for       !
! parallel execution; NEC ASL libraries (3d only, no parallel execution)   !
! Written by Carlo Cavazzoni, modified by P. Giannozzi, contributions      !
! by Martin Hilgemans, Guido Roma, Pascal Thibaudeau, Stephane Lefranc,    !
! Nicolas Lacorne, Filippo Spiga, Nicola Varini - Last update Jul 2015     !
!--------------------------------------------------------------------------!

#include "fft_defs.h"

#if defined __DFTI
#include "mkl_dfti.f90"
#endif

!=----------------------------------------------------------------------=!
   MODULE fft_scalar
!=----------------------------------------------------------------------=!
       USE kinds

       USE, intrinsic ::  iso_c_binding
       
#if defined __DFTI
       USE MKL_DFTI ! -- this can be found in the MKL include directory
#endif
       IMPLICIT NONE
        SAVE

        PRIVATE
        PUBLIC :: cft_1z, cft_2xy, cft_b, cfft3d, cfft3ds
        PUBLIC :: good_fft_dimension, allowed, good_fft_order
        PUBLIC :: cft_b_omp_init, cft_b_omp

! ...   Local Parameter

        !   ndims   Number of different FFT tables that the module
        !           could keep into memory without reinitialization
        !   nfftx   Max allowed fft dimension

        INTEGER, PARAMETER :: ndims = 3, nfftx = 2049

        !   Workspace that is statically allocated is defined here
        !   in order to avoid multiple copies of the same workspace
        !   lwork:   Dimension of the work space array (if any)

#if defined __FFTW 

        INTEGER   :: cft_b_dims( 4 )
!$omp threadprivate (cft_b_dims)
        C_POINTER :: cft_b_bw_planz = 0
!$omp threadprivate (cft_b_bw_planz)
        C_POINTER :: cft_b_bw_planx = 0
!$omp threadprivate (cft_b_bw_planx)
        C_POINTER :: cft_b_bw_plany = 0
!$omp threadprivate (cft_b_bw_plany)

#elif defined __FFTW3

#if defined __OPENMP
#include "fftw3.f03"
#else
#include "fftw3.f"
#endif

#elif defined __DFTI
        TYPE dfti_descriptor_array
           TYPE(DFTI_DESCRIPTOR), POINTER :: desc
        END TYPE

#elif ( defined __ESSL || defined __LINUX_ESSL )

        !   ESSL IBM library: see the ESSL manual for DCFT

        INTEGER, PARAMETER :: lwork = 20000 + ( 2*nfftx + 256 ) * 64 + 3*nfftx
        REAL (DP) :: work( lwork )

#endif

!=----------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------=!

!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "z"
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cft_1z(c, nsl, nz, ldz, isign, cout)

!     driver routine for nsl 1d complex fft's of length nz
!     ldz >= nz is the distance between sequences to be transformed
!     (ldz>nz is used on some architectures to reduce memory conflicts)
!     input  :  c(ldz*nsl)   (complex)
!     output : cout(ldz*nsl) (complex - NOTA BENE: transform is not in-place!)
!     isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nz, nsl, ldz) are stored and re-used if available

     INTEGER, INTENT(IN) :: isign
     INTEGER, INTENT(IN) :: nsl, nz, ldz

     COMPLEX (DP) :: c(:), cout(:)

     REAL (DP)  :: tscale
     INTEGER    :: i, err, idir, ip, void
     INTEGER, SAVE :: zdims( 3, ndims ) = -1
     INTEGER, SAVE :: icurrent = 1
     LOGICAL :: done

#if defined __HPM
     INTEGER :: OMP_GET_THREAD_NUM
#endif
     INTEGER :: tid

     ! ...   Machine-Dependent parameters, work arrays and tables of factors

     !   ltabl   Dimension of the tables of factors calculated at the
     !           initialization stage

#if defined __OPENMP
     INTEGER :: offset, ldz_t
     INTEGER :: omp_get_max_threads
     EXTERNAL :: omp_get_max_threads
#endif


#if defined __FFTW || defined __FFTW3

     !   Pointers to the "C" structures containing FFT factors ( PLAN )
     !   C_POINTER is defined in include/fft_defs.h
     !   for 32bit executables, C_POINTER is integer(4)
     !   for 64bit executables, C_POINTER is integer(8)

     C_POINTER, SAVE :: fw_planz( ndims ) = 0
     C_POINTER, SAVE :: bw_planz( ndims ) = 0

#elif defined __DFTI

     !   Intel MKL native FFT driver

     TYPE(DFTI_DESCRIPTOR_ARRAY), SAVE :: hand( ndims )
     LOGICAL, SAVE :: dfti_first = .TRUE.
     INTEGER :: dfti_status = 0
     !
     IF( dfti_first .EQ. .TRUE. ) THEN
        DO ip = 1, ndims
           hand(ip)%desc => NULL()
        END DO
        dfti_first = .FALSE.
     END IF

#elif defined __ESSL || defined __LINUX_ESSL

     !   ESSL IBM library: see the ESSL manual for DCFT

     INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx
     REAL (DP), SAVE :: fw_tablez( ltabl, ndims )
     REAL (DP), SAVE :: bw_tablez( ltabl, ndims )

#elif defined __SX6

     !   NEC MathKeisan

     INTEGER, PARAMETER :: ltabl = 2 * nfftx + 64
     REAL (DP), SAVE :: tablez (ltabl, ndims)
     REAL (DP)       :: work(4*nz*nsl)
     COMPLEX (DP)    :: DUMMY
     INTEGER, SAVE :: isys = 1

#endif

     IF( nsl < 0 ) THEN
       CALL errore(" fft_scalar: cft_1z ", " nsl out of range ", nsl)
     END IF

     !
     !   Here initialize table only if necessary
     !

     DO ip = 1, ndims

        !   first check if there is already a table initialized
        !   for this combination of parameters

        done = ( nz == zdims(1,ip) )
#if defined __ESSL || defined __LINUX_ESSL || defined __FFTW3 || __DFTI

        !   The initialization in ESSL and FFTW v.3 depends on all three parameters

        done = done .AND. ( nsl == zdims(2,ip) ) .AND. ( ldz == zdims(3,ip) )
#endif
        IF (done) EXIT
     END DO

     IF( .NOT. done ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       ! WRITE( stdout, fmt="('DEBUG cft_1z, reinitializing tables ', I3)" ) icurrent
#if defined __FFTW

       IF( fw_planz( icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_planz( icurrent) )
       IF( bw_planz( icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_planz( icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_planz( icurrent), nz, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_planz( icurrent), nz, idir)

#elif defined __FFTW3

#if defined __OPENMP
       CALL dfftw_cleanup_threads() 
       void = fftw_init_threads()
       CALL dfftw_plan_with_nthreads(omp_get_max_threads())      
#endif

       IF( fw_planz( icurrent) /= 0 ) CALL dfftw_destroy_plan( fw_planz( icurrent) )
       IF( bw_planz( icurrent) /= 0 ) CALL dfftw_destroy_plan( bw_planz( icurrent) )
       idir = -1
       CALL dfftw_plan_many_dft( fw_planz( icurrent), 1, nz, nsl, c, &
            (/SIZE(c)/), 1, ldz, cout, (/SIZE(cout)/), 1, ldz, idir, FFTW_ESTIMATE)
       idir = 1
       CALL dfftw_plan_many_dft( bw_planz( icurrent), 1, nz, nsl, c, &
            (/SIZE(c)/), 1, ldz, cout, (/SIZE(cout)/), 1, ldz, idir, FFTW_ESTIMATE)

#elif defined __DFTI

       CALL qe_dfti_init()

#elif defined __ESSL || defined __LINUX_ESSL

       tscale = 1.0_DP / nz

       CALL DCFT ( 1, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl,  1, &
          tscale, fw_tablez(1, icurrent), ltabl, work(1), lwork)
       CALL DCFT ( 1, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, -1, &
          1.0_DP, bw_tablez(1, icurrent), ltabl, work(1), lwork)


#elif defined __SX6

       CALL ZZFFTM (0, nz, 1, 1.0_DP, DUMMY, ldz, DUMMY, ldz, &
                    tablez (1, icurrent), work, isys)
#else

       CALL errore(' cft_1z ',' no scalar fft driver specified ', 1)

#endif

       zdims(1,icurrent) = nz; zdims(2,icurrent) = nsl; zdims(3,icurrent) = ldz;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined __FFT_CLOCKS
     CALL start_clock( 'cft_1z' )
#endif

#if defined __FFTW

#if defined __OPENMP

     ldz_t = ldz
     !
     IF (isign < 0) THEN
!$omp parallel default(none) private(tid,offset,i,tscale) shared(c,isign,nsl,fw_planz,ip,nz,cout,ldz) &
!$omp &        firstprivate(ldz_t)
!$omp do
       DO i=1, nsl
          offset = 1 + ((i-1)*ldz_t)
          CALL FFT_Z_STICK_SINGLE(fw_planz( ip), c(offset), ldz_t)
       END DO
!$omp end do
!$omp end parallel
       tscale = 1.0_DP / nz
       cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl ) * tscale
     ELSE IF (isign > 0) THEN
!$omp parallel default(none) private(tid,offset,i) shared(c,isign,nsl,bw_planz,ip,cout,ldz) &
!$omp &        firstprivate(ldz_t)
!$omp do
       DO i=1, nsl
          offset =  1 + ((i-1)* ldz_t)
          CALL FFT_Z_STICK_SINGLE(bw_planz( ip), c(offset), ldz_t)
       END DO
!$omp end do
!$omp workshare
       cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl )
!$omp end workshare
!$omp end parallel
     END IF


#else

     IF (isign < 0) THEN
        CALL FFT_Z_STICK(fw_planz( ip), c(1), ldz, nsl)
        tscale = 1.0_DP / nz
        cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl ) * tscale
     ELSE IF (isign > 0) THEN
        CALL FFT_Z_STICK(bw_planz( ip), c(1), ldz, nsl)
        cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl )
     END IF

#endif

#elif defined __FFTW3

     IF (isign < 0) THEN
        CALL dfftw_execute_dft( fw_planz( ip), c, cout)
        tscale = 1.0_DP / nz
        cout( 1 : ldz * nsl ) = cout( 1 : ldz * nsl ) * tscale
     ELSE IF (isign > 0) THEN
        CALL dfftw_execute_dft( bw_planz( ip), c, cout)
     END IF

#elif defined __DFTI

     IF (isign < 0) THEN
        dfti_status = DftiComputeForward(hand(ip)%desc, c )
        IF(dfti_status /= 0) &
           CALL errore(' cft_1z ',' stopped in DftiComputeForward ', dfti_status )
        cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl )
     ELSE IF (isign > 0) THEN
        dfti_status = DftiComputeBackward(hand(ip)%desc, c )
        IF(dfti_status /= 0) &
           CALL errore(' cft_1z ',' stopped in DftiComputeBackward ', dfti_status )
        cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl )
     END IF

#elif defined __SX6

     IF ( isign < 0 ) THEN
        idir   = -1
        tscale = 1.0_DP / nz
     ELSE IF ( isign > 0 ) THEN
        idir   = 1
        tscale = 1.0_DP
     END IF
     IF (isign /= 0) CALL ZZFFTM (idir, nz, nsl, tscale, c(1), ldz, &
          cout(1), ldz, tablez (1, ip), work, isys)

#elif defined __ESSL || defined __LINUX_ESSL

     ! essl uses a different convention for forward/backward transforms
     ! wrt most other implementations: notice the sign of "idir"

     IF( isign < 0 ) THEN
        idir   =+1
        tscale = 1.0_DP / nz
        CALL DCFT (0, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, idir, &
             tscale, fw_tablez(1, ip), ltabl, work, lwork)
     ELSE IF( isign > 0 ) THEN
        idir   =-1
        tscale = 1.0_DP
        CALL DCFT (0, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, idir, &
             tscale, bw_tablez(1, ip), ltabl, work, lwork)
     END IF

#else

    CALL errore(' cft_1z ',' no scalar fft driver specified ', 1)

#endif

#if defined __FFT_CLOCKS
     CALL stop_clock( 'cft_1z' )
#endif

     RETURN


   CONTAINS

     SUBROUTINE dummy_init()
        RETURN
     END SUBROUTINE dummy_init

#ifdef __DFTI
     SUBROUTINE qe_dfti_init()

       if( ASSOCIATED( hand( icurrent )%desc ) ) THEN
          dfti_status = DftiFreeDescriptor( hand( icurrent )%desc )
          IF( dfti_status /= 0) THEN
             WRITE(*,*) "stopped in DftiFreeDescriptor", dfti_status
             STOP
          ENDIF
       END IF

     dfti_status = DftiCreateDescriptor(hand( icurrent )%desc, DFTI_DOUBLE, DFTI_COMPLEX, 1,nz)
     IF(dfti_status /= 0)  &
       CALL errore(' cft_1z ',' stopped in DftiCreateDescriptor ', dfti_status )

     dfti_status = DftiSetValue(hand( icurrent )%desc, DFTI_NUMBER_OF_TRANSFORMS,nsl)
     IF(dfti_status /= 0) &
       CALL errore(' cft_1z ',' stopped in DFTI_NUMBER_OF_TRANSFORMS ', dfti_status )

     dfti_status = DftiSetValue(hand( icurrent )%desc,DFTI_INPUT_DISTANCE, ldz )
     IF(dfti_status /= 0) &
       CALL errore(' cft_1z ',' stopped in DFTI_INPUT_DISTANCE ', dfti_status )

     dfti_status = DftiSetValue(hand( icurrent )%desc, DFTI_PLACEMENT, DFTI_INPLACE)
     IF(dfti_status /= 0) &
       CALL errore(' cft_1z ',' stopped in DFTI_PLACEMENT ', dfti_status )

     tscale = 1.0_DP/nz
     dfti_status = DftiSetValue( hand( icurrent )%desc, DFTI_FORWARD_SCALE, tscale);
     IF(dfti_status /= 0) &
       CALL errore(' cft_1z ',' stopped in DFTI_FORWARD_SCALE ', dfti_status )

     dfti_status = DftiSetValue( hand( icurrent )%desc, DFTI_BACKWARD_SCALE, DBLE(1) );
     IF(dfti_status /= 0) &
       CALL errore(' cft_1z ',' stopped in DFTI_BACKWARD_SCALE ', dfti_status )

     dfti_status = DftiCommitDescriptor(hand( icurrent )%desc)
     IF(dfti_status /= 0) &
       CALL errore(' cft_1z ',' stopped in DftiCommitDescriptor ', dfti_status )

     END SUBROUTINE qe_dfti_init
#endif


   END SUBROUTINE cft_1z

!
!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "x" and "y" direction
!
!
!
!=----------------------------------------------------------------------=!
!
!

   SUBROUTINE cft_2xy(r, nzl, nx, ny, ldx, ldy, isign, pl2ix)

!     driver routine for nzl 2d complex fft's of lengths nx and ny
!     input : r(ldx*ldy)  complex, transform is in-place
!     ldx >= nx, ldy >= ny are the physical dimensions of the equivalent
!     2d array: r2d(ldx, ldy) (x first dimension, y second dimension)
!     (ldx>nx, ldy>ny used on some architectures to reduce memory conflicts)
!     pl2ix(nx) (optional) is 1 for columns along y to be transformed
!     isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nx,ny,nzl,ldx) are stored and re-used if available

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: isign, ldx, ldy, nx, ny, nzl
     INTEGER, OPTIONAL, INTENT(IN) :: pl2ix(:)
     COMPLEX (DP) :: r( : )
     INTEGER :: i, k, j, err, idir, ip, kk, void
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims( 4, ndims) = -1
     LOGICAL :: dofft( nfftx ), done
     INTEGER, PARAMETER  :: stdout = 6

#if defined __HPM
     INTEGER :: OMP_GET_THREAD_NUM
#endif
#if defined __OPENMP
     INTEGER :: offset
     INTEGER :: nx_t, ny_t, nzl_t, ldx_t, ldy_t
     INTEGER  :: itid, mytid, ntids
     INTEGER  :: omp_get_thread_num, omp_get_num_threads,omp_get_max_threads
     EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

#if defined __DFTI

     TYPE(DFTI_DESCRIPTOR_ARRAY), SAVE :: hand( ndims )
     LOGICAL, SAVE :: dfti_first = .TRUE.
     INTEGER :: dfti_status = 0

#elif defined __FFTW || defined __FFTW3

#if defined __FFTW && __FFTW_ALL_XY_PLANES
     C_POINTER, SAVE :: fw_plan_2d( ndims ) = 0
     C_POINTER, SAVE :: bw_plan_2d( ndims ) = 0
#else
     C_POINTER, SAVE :: fw_plan( 2, ndims ) = 0
     C_POINTER, SAVE :: bw_plan( 2, ndims ) = 0
#endif


#elif defined __ESSL || defined __LINUX_ESSL

     INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx
     REAL (DP), SAVE :: fw_tablex( ltabl, ndims ), fw_tabley( ltabl, ndims )
     REAL (DP), SAVE :: bw_tablex( ltabl, ndims ), bw_tabley( ltabl, ndims )

#elif defined __SX6

     INTEGER, PARAMETER :: ltabl = 2*nfftx + 64
     REAL (DP), SAVE :: tablex(ltabl, ndims), tabley(ltabl, ndims)
     REAL (DP)       :: work(4*nx*ny)
     COMPLEX (DP) :: XY(ldx*ny)
     COMPLEX (DP) :: DUMMY
     INTEGER, SAVE :: isys = 1

#endif


     dofft( 1 : nx ) = .TRUE.
     IF( PRESENT( pl2ix ) ) THEN
       IF( SIZE( pl2ix ) < nx ) &
         CALL errore( ' cft_2xy ', ' wrong dimension for arg no. 8 ', 1 )
       DO i = 1, nx
         IF( pl2ix(i) < 1 ) dofft( i ) = .FALSE.
       END DO
     END IF

     ! WRITE( stdout,*) 'DEBUG: ', COUNT( dofft )

     !
     !   Here initialize table only if necessary
     !
#if defined __DFTI
     IF( dfti_first .EQ. .TRUE. ) THEN
        DO ip = 1, ndims
           hand(ip)%desc => NULL()
        END DO
        dfti_first = .FALSE.
     END IF
#endif

     DO ip = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       done = ( ny == dims(1,ip) ) .AND. ( nx == dims(3,ip) )
       done = done .AND. ( ldx == dims(2,ip) ) .AND.  ( nzl == dims(4,ip) )
       IF (done) EXIT

     END DO

     IF( .NOT. done ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       ! WRITE( stdout, fmt="('DEBUG cft_2xy, reinitializing tables ', I3)" ) icurrent

#if defined __DFTI

       if( ASSOCIATED( hand( icurrent )%desc ) ) THEN
          dfti_status = DftiFreeDescriptor( hand( icurrent )%desc )
          IF( dfti_status /= 0) THEN
             WRITE(*,*) "stopped in DftiFreeDescriptor", dfti_status
             STOP
          ENDIF
       END IF

       dfti_status = DftiCreateDescriptor(hand( icurrent )%desc, DFTI_DOUBLE, DFTI_COMPLEX, 2,(/nx,ny/))
       IF(dfti_status /= 0) THEN
          WRITE(*,*) "stopped in DftiCreateDescriptor", dfti_status
          STOP
       ENDIF
       dfti_status = DftiSetValue(hand( icurrent )%desc, DFTI_NUMBER_OF_TRANSFORMS,nzl)
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DFTI_NUMBER_OF_TRANSFORMS", dfti_status
          STOP
       ENDIF
       dfti_status = DftiSetValue(hand( icurrent )%desc,DFTI_INPUT_DISTANCE, ldx*ldy )
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DFTI_INPUT_DISTANCE", dfti_status
          STOP
       ENDIF
       dfti_status = DftiSetValue(hand( icurrent )%desc, DFTI_PLACEMENT, DFTI_INPLACE)
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DFTI_PLACEMENT", dfti_status
          STOP
       ENDIF
       tscale = 1.0_DP/ (nx * ny )
       dfti_status = DftiSetValue( hand( icurrent )%desc, DFTI_FORWARD_SCALE, tscale);
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DFTI_FORWARD_SCALE", dfti_status
          STOP
       ENDIF
       dfti_status = DftiSetValue( hand( icurrent )%desc, DFTI_BACKWARD_SCALE, DBLE(1) );
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DFTI_BACKWARD_SCALE", dfti_status
          STOP
       ENDIF
       dfti_status = DftiCommitDescriptor(hand( icurrent )%desc)
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DftiCommitDescriptor", dfti_status
          STOP
       ENDIF


#elif defined __FFTW

#if defined __FFTW_ALL_XY_PLANES
       IF( fw_plan_2d( icurrent) /= 0 )  CALL DESTROY_PLAN_2D(fw_plan_2d(icurrent) )
       IF( bw_plan_2d( icurrent) /= 0 )  CALL DESTROY_PLAN_2D(bw_plan_2d(icurrent) )
       idir = -1; CALL CREATE_PLAN_2D( fw_plan_2d(icurrent), nx, ny, idir)
       idir =  1; CALL CREATE_PLAN_2D( bw_plan_2d(icurrent), nx, ny, idir)
#else
       IF( fw_plan( 2,icurrent) /= 0 )   CALL DESTROY_PLAN_1D( fw_plan( 2,icurrent) )
       IF( bw_plan( 2,icurrent) /= 0 )   CALL DESTROY_PLAN_1D( bw_plan( 2,icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 2,icurrent), ny, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 2,icurrent), ny, idir)

       IF( fw_plan( 1,icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_plan( 1,icurrent) )
       IF( bw_plan( 1,icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_plan( 1,icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 1,icurrent), nx, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 1,icurrent), nx, idir)
#endif

#elif defined __FFTW3

#if defined __OPENMP
       CALL dfftw_cleanup_threads() 
       void = fftw_init_threads()
       CALL dfftw_plan_with_nthreads(omp_get_max_threads())      
#endif

       IF ( ldx /= nx .OR. ldy /= ny ) THEN
          IF( fw_plan(2,icurrent) /= 0 )  CALL dfftw_destroy_plan( fw_plan(2,icurrent) )
          IF( bw_plan(2,icurrent) /= 0 )  CALL dfftw_destroy_plan( bw_plan(2,icurrent) )
          idir = -1
          CALL dfftw_plan_many_dft( fw_plan(2,icurrent), 1, ny, 1, r(1:), &
               (/ldx*ldy/), ldx, 1, r(1:), (/ldx*ldy/), ldx, 1, idir, &
               FFTW_ESTIMATE)
          idir =  1
          CALL dfftw_plan_many_dft( bw_plan(2,icurrent), 1, ny, 1, r(1:), &
               (/ldx*ldy/), ldx, 1, r(1:), (/ldx*ldy/), ldx, 1, idir, &
               FFTW_ESTIMATE)

          IF( fw_plan(1,icurrent) /= 0 ) CALL dfftw_destroy_plan( fw_plan(1,icurrent) )
          IF( bw_plan(1,icurrent) /= 0 ) CALL dfftw_destroy_plan( bw_plan(1,icurrent) )
          idir = -1
          CALL dfftw_plan_many_dft( fw_plan(1,icurrent), 1, nx, ny, r(1:), &
               (/ldx*ldy/), 1, ldx, r(1:), (/ldx*ldy/), 1, ldx, idir, &
               FFTW_ESTIMATE)
          idir =  1
          CALL dfftw_plan_many_dft( bw_plan(1,icurrent), 1, nx, ny, r(1:), &
               (/ldx*ldy/), 1, ldx, r(1:), (/ldx*ldy/), 1, ldx, idir, &
               FFTW_ESTIMATE)
       ELSE
          IF( fw_plan( 1, icurrent) /= 0 ) CALL dfftw_destroy_plan( fw_plan( 1, icurrent) )
          IF( bw_plan( 1, icurrent) /= 0 ) CALL dfftw_destroy_plan( bw_plan( 1, icurrent) )
          idir = -1
          CALL dfftw_plan_many_dft( fw_plan( 1, icurrent), 2, (/nx, ny/), nzl,&
               r(1:), (/nx, ny/), 1, nx*ny, r(1:), (/nx, ny/), 1, nx*ny, idir,&
               FFTW_ESTIMATE)
          idir = 1
          CALL dfftw_plan_many_dft( bw_plan( 1, icurrent), 2, (/nx, ny/), nzl,&
               r(1:), (/nx, ny/), 1, nx*ny, r(1:), (/nx, ny/), 1, nx*ny, idir,&
               FFTW_ESTIMATE)
       END IF

#elif defined __ESSL || defined __LINUX_ESSL

#if defined __OPENMP

       tscale = 1.0_DP / ( nx * ny )
       CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, nx,  1, 1.0_DP, &
          fw_tabley( 1, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, nx, -1, 1.0_DP, &
          bw_tabley(1, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny,  1, &
          tscale, fw_tablex( 1, icurrent), ltabl, work(1), lwork)
       CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny, -1, &
          1.0_DP, bw_tablex(1, icurrent), ltabl, work(1), lwork)

#else

       tscale = 1.0_DP / ( nx * ny )
       CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, 1,  1, 1.0_DP, &
          fw_tabley( 1, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, 1, -1, 1.0_DP, &
          bw_tabley(1, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny,  1, &
          tscale, fw_tablex( 1, icurrent), ltabl, work(1), lwork)
       CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny, -1, &
          1.0_DP, bw_tablex(1, icurrent), ltabl, work(1), lwork)

#endif

#elif defined __SX6


       CALL ZZFFT(0, ny, 1.0_DP, DUMMY, DUMMY,              &
                  tabley (1, icurrent), work, isys)
       CALL ZZFFTM  (0, nx, 1, 1.0_DP, DUMMY, ldx, DUMMY, ldx,           &
                     tablex(1, icurrent), work, isys)
#else

       CALL errore(' cft_2xy ',' no scalar fft driver specified ', 1)

#endif

       dims(1,icurrent) = ny; dims(2,icurrent) = ldx;
       dims(3,icurrent) = nx; dims(4,icurrent) = nzl;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined __FFT_CLOCKS
     CALL start_clock( 'cft_2xy' )
#endif


#if defined __DFTI

     IF( isign < 0 ) THEN
        !
        dfti_status = DftiComputeForward(hand(ip)%desc, r(:))
        IF(dfti_status /= 0)THEN
           WRITE(*,*) "stopped in DftiComputeForward", dfti_status
           STOP
        ENDIF
        !
     ELSE IF( isign > 0 ) THEN
        !
        dfti_status = DftiComputeBackward(hand(ip)%desc, r(:))
        IF(dfti_status /= 0)THEN
           WRITE(*,*) "stopped in DftiComputeBackward", dfti_status
           STOP
        ENDIF
        !
     END IF


#elif defined __FFTW

#if defined __FFTW_ALL_XY_PLANES

     IF( isign < 0 ) THEN
        !
        tscale = 1.0_DP / ( nx * ny )
        !
        CALL fftw_inplace_drv_2d( fw_plan_2d(ip), nzl, r(1), 1, ldx*ldy )
        CALL ZDSCAL( ldx * ldy * nzl, tscale, r(1), 1)
        !
     ELSE IF( isign > 0 ) THEN
        !
        CALL fftw_inplace_drv_2d( bw_plan_2d(ip), nzl, r(1), 1, ldx*ldy )
        !
     END IF

#elif defined __OPENMP

     nx_t  = nx
     ny_t  = ny
     nzl_t = nzl
     ldx_t = ldx 
     ldy_t = ldy
     !
     IF( isign < 0 ) THEN
        !
        tscale = 1.0_DP / ( nx * ny )
        !
!$omp parallel default(none) private(offset,itid,mytid,ntids,k,j,i)  shared(r,dofft,ip,fw_plan,nzl,nx,ny,ldx,ldy,tscale)  &
!$omp & firstprivate(nx_t, ny_t, nzl_t, ldx_t, ldy_t)

!$omp do
        DO i=1,nzl
           offset = 1+ ((i-1)*(ldx_t*ldy_t))
           CALL FFT_X_STICK_SINGLE( fw_plan(1,ip), r(offset), nx_t, ny_t, nzl_t, ldx_t, ldy_t )
        END DO
!$omp end do

        mytid = omp_get_thread_num()  ! take the thread ID
        ntids = omp_get_num_threads() ! take the number of threads
        itid  = 0

        do i = 1, nx
          do k = 1, nzl
            IF( dofft( i ) ) THEN
              IF( itid == mytid ) THEN
                j = i + ldx_t*ldy_t * ( k - 1 )
                call FFT_Y_STICK(fw_plan(2,ip), r(j), ny_t, ldx_t)
              END IF
              itid = MOD( itid + 1, ntids )
            END IF
          end do
        end do

!$omp barrier
 
!$omp workshare
        r = r * tscale
!$omp end workshare

!$omp end parallel

        ! CALL ZDSCAL( ldx * ldy * nzl, tscale, r(1), 1)
        !
     ELSE IF( isign > 0 ) THEN
        !
!$omp parallel default(none) private(offset,itid,mytid,ntids,k,j,i) shared(r,nx,nzl,dofft,ip,bw_plan) &
!$omp & firstprivate(nx_t, ny_t, nzl_t, ldx_t, ldy_t)

        mytid = omp_get_thread_num()  ! take the thread ID
        ntids = omp_get_num_threads() ! take the number of threads
        itid  = 0

        do i = 1, nx
          do k = 1, nzl
            IF( dofft( i ) ) THEN
              IF( itid == mytid ) THEN
                j = i + ldx_t*ldy_t * ( k - 1 )
                call FFT_Y_STICK( bw_plan(2,ip), r(j), ny_t, ldx_t)
              END IF
              itid = MOD( itid + 1, ntids )
            END IF
          end do
        end do

!$omp barrier

!$omp do
        DO i=1,nzl
           offset = 1+ ((i-1)*(ldx_t*ldy_t))
           CALL FFT_X_STICK_SINGLE( bw_plan(1,ip), r(offset), nx_t, ny_t, nzl_t, ldx_t, ldy_t )
        END DO
!$omp end do
!$omp end parallel
        !
     END IF

#else

     IF( isign < 0 ) THEN

       CALL FFT_X_STICK( fw_plan(1,ip), r(1), nx, ny, nzl, ldx, ldy )

       do i = 1, nx
         do k = 1, nzl
           IF( dofft( i ) ) THEN
             j = i + ldx*ldy * ( k - 1 )
             call FFT_Y_STICK(fw_plan(2,ip), r(j), ny, ldx)
           END IF
         end do
       end do
       tscale = 1.0_DP / ( nx * ny )
       CALL ZDSCAL( ldx * ldy * nzl, tscale, r(1), 1)

     ELSE IF( isign > 0 ) THEN

       do i = 1, nx
         do k = 1, nzl
           IF( dofft( i ) ) THEN
             j = i + ldx*ldy * ( k - 1 )
             call FFT_Y_STICK( bw_plan(2,ip), r(j), ny, ldx)
           END IF
         end do
       end do

       CALL FFT_X_STICK( bw_plan(1,ip), r(1), nx, ny, nzl, ldx, ldy )

    END IF

#endif

#elif defined __FFTW3

     IF ( ldx /= nx .OR. ldy /= ny ) THEN
        IF( isign < 0 ) THEN
           do j = 0, nzl-1
              CALL dfftw_execute_dft( fw_plan (1, ip), &
                   r(1+j*ldx*ldy:), r(1+j*ldx*ldy:))
           end do
           do i = 1, nx
              do k = 1, nzl
                 IF( dofft( i ) ) THEN
                    j = i + ldx*ldy * ( k - 1 )
                    call dfftw_execute_dft( fw_plan ( 2, ip), r(j:), r(j:))
                 END IF
              end do
           end do
           tscale = 1.0_DP / ( nx * ny )
           CALL ZDSCAL( ldx * ldy * nzl, tscale, r(1), 1)
        ELSE IF( isign > 0 ) THEN
           do i = 1, nx
              do k = 1, nzl
                 IF( dofft( i ) ) THEN
                    j = i + ldx*ldy * ( k - 1 )
                    call dfftw_execute_dft( bw_plan ( 2, ip), r(j:), r(j:))
                 END IF
              end do
           end do
           do j = 0, nzl-1
              CALL dfftw_execute_dft( bw_plan( 1, ip), &
                   r(1+j*ldx*ldy:), r(1+j*ldx*ldy:))
           end do
        END IF
     ELSE
        IF( isign < 0 ) THEN
           call dfftw_execute_dft( fw_plan( 1, ip), r(1:), r(1:))
           tscale = 1.0_DP / ( nx * ny )
           CALL ZDSCAL( ldx * ldy * nzl, tscale, r(1), 1)
        ELSE IF( isign > 0 ) THEN
           call dfftw_execute_dft( bw_plan( 1, ip), r(1:), r(1:))
        END IF
     END IF

#elif defined __ESSL || defined __LINUX_ESSL

#if defined __OPENMP

   IF( isign < 0 ) THEN
      tscale = 1.0_DP / ( nx * ny )
      do k = 1, nzl
         kk = 1 + ( k - 1 ) * ldx * ldy
         CALL DCFT ( 0, r( kk ), 1, ldx, r( kk ), 1, ldx, nx, ny, &
              1, tscale, fw_tablex( 1, ip ), ltabl, work( 1 ), lwork)
         CALL DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, nx, &
              1, 1.0_DP, fw_tabley(1, ip), ltabl, work( 1 ), lwork)
      end do
   ELSE IF( isign > 0 ) THEN
      DO k = 1, nzl
         kk = 1 + ( k - 1 ) * ldx * ldy
         CALL DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, nx, &
                   -1, 1.0_DP, bw_tabley(1, ip), ltabl, work( 1 ), lwork)
         CALL DCFT ( 0, r( kk ), 1, ldx, r( kk ), 1, ldx, nx, ny, &
                   -1, 1.0_DP, bw_tablex(1, ip), ltabl, work( 1 ), lwork)
      END DO
   END IF

#else

   IF( isign < 0 ) THEN
      idir = 1
      tscale = 1.0_DP / ( nx * ny )
      do k = 1, nzl
         kk = 1 + ( k - 1 ) * ldx * ldy
         CALL DCFT ( 0, r(kk), 1, ldx, r(kk), 1, ldx, nx, ny, idir, &
              tscale, fw_tablex( 1, ip ), ltabl, work( 1 ), lwork)
         do i = 1, nx
            IF( dofft( i ) ) THEN
               kk = i + ( k - 1 ) * ldx * ldy
               call DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, 1, &
                    idir, 1.0_DP, fw_tabley(1, ip), ltabl, work( 1 ), lwork)
            END IF
         end do
      end do
   ELSE IF( isign > 0 ) THEN
      idir = -1
      DO k = 1, nzl
         do i = 1, nx
            IF( dofft( i ) ) THEN
               kk = i + ( k - 1 ) * ldx * ldy
               call DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, 1, &
                    idir, 1.0_DP, bw_tabley(1, ip), ltabl, work( 1 ), lwork)
            END IF
         end do
         kk = 1 + ( k - 1 ) * ldx * ldy
         CALL DCFT ( 0, r( kk ), 1, ldx, r( kk ), 1, ldx, nx, ny, idir, &
              1.0_DP, bw_tablex(1, ip), ltabl, work( 1 ), lwork)
      END DO
   END IF
#endif

#elif defined __SX6

      IF( isign < 0 ) THEN

       idir = -1
       tscale = 1.0_DP / (nx * ny)
       DO k = 0, nzl-1
          kk = k * ldx * ldy
! FORWARD: ny FFTs in the X direction
          CALL ZZFFTM ( idir, nx, ny, tscale, r(kk+1), ldx, r(kk+1), ldx,   &
                        tablex (1, ip), work(1), isys )
! FORWARD: nx FFTs in the Y direction
          DO i = 1, nx
             IF ( dofft(i) ) THEN
                DO j = 0, ny-1
                   XY(j+1) = r(i + (j) * ldx + kk)
                END DO
                CALL ZZFFT(idir, ny, 1.0_DP, XY, XY, tabley (1, ip),      &
                           work(1), isys)
                DO j = 0, ny-1
                   r(i + (j) * ldx + kk) = XY(j+1)
                END DO
             END IF
          END DO
       END DO

     ELSE IF ( isign > 0 ) THEN

       idir = 1
       tscale = 1.0_DP
       DO k = 0, nzl-1
! BACKWARD: nx FFTs in the Y direction
          kk = (k) * ldx * ldy
          DO i = 1, nx
             IF ( dofft(i) ) THEN
                DO j = 0, ny-1
                   XY(j+1) = r(i + (j) * ldx + kk)
                END DO
                CALL ZZFFT(idir, ny, 1.0_DP, XY, XY, tabley (1, ip),      &
                           work(1), isys)
                DO j = 0, ny-1
                   r(i + (j) * ldx + kk) = XY(j+1)
                END DO
             END IF
          END DO
! BACKWARD: ny FFTs in the X direction
          CALL ZZFFTM ( idir, nx, ny, tscale, r(kk+1), ldx, r(kk+1), ldx,   &
                        tablex (1, ip), work(1), isys )
       END DO

     END IF

#else

     CALL errore(' cft_2xy ',' no scalar fft driver specified ', 1)

#endif

#if defined __FFT_CLOCKS
     CALL stop_clock( 'cft_2xy' )
#endif

     RETURN

   END SUBROUTINE cft_2xy


!
!=----------------------------------------------------------------------=!
!
!
!
!         3D scalar FFTs
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cfft3d( f, nx, ny, nz, ldx, ldy, ldz, isign )

  !     driver routine for 3d complex fft of lengths nx, ny, nz
  !     input  :  f(ldx*ldy*ldz)  complex, transform is in-place
  !     ldx >= nx, ldy >= ny, ldz >= nz are the physical dimensions
  !     of the equivalent 3d array: f3d(ldx,ldy,ldz)
  !     (ldx>nx, ldy>ny, ldz>nz may be used on some architectures
  !      to reduce memory conflicts - not implemented for FFTW)
  !     isign > 0 : f(G) => f(R)   ; isign < 0 : f(R) => f(G)
  !
  !     Up to "ndims" initializations (for different combinations of input
  !     parameters nx,ny,nz) are stored and re-used if available

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nx, ny, nz, ldx, ldy, ldz, isign
     COMPLEX (DP) :: f(:)
     INTEGER :: i, k, j, err, idir, ip
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims(3,ndims) = -1

#if defined __FFTW || defined __FFTW3

     C_POINTER, save :: fw_plan(ndims) = 0
     C_POINTER, save :: bw_plan(ndims) = 0

#elif defined __DFTI

     !   Intel MKL native FFT driver

     TYPE(DFTI_DESCRIPTOR_ARRAY), SAVE :: hand(ndims)
     LOGICAL, SAVE :: dfti_first = .TRUE.
     INTEGER :: dfti_status = 0
     !
     IF( dfti_first .EQ. .TRUE. ) THEN
        DO ip = 1, ndims
           hand(ip)%desc => NULL()
        END DO
 
        dfti_first = .FALSE.
     END IF

#elif defined __SX6

     INTEGER, PARAMETER :: ltabl = 60
     INTEGER, PARAMETER :: lwork = 195+6*nfftx
     INTEGER, SAVE  :: iw0(ltabl, ndims)
     INTEGER :: k_off, kj_offset
     REAL (DP), SAVE :: auxp (lwork, ndims)
     ! not sure whether auxp is work space or not
     COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: cw2
     COMPLEX (DP) :: f_out(size(f))

#  if defined ASL && defined MICRO
     INTEGER :: nbtasks
     COMMON/NEC_ASL_PARA/nbtasks
#  endif

#endif

     IF ( nx < 1 ) &
         call errore('cfft3d',' nx is less than 1 ', 1)
     IF ( ny < 1 ) &
         call errore('cfft3d',' ny is less than 1 ', 1)
     IF ( nz < 1 ) &
         call errore('cfft3',' nz is less than 1 ', 1)

#if defined __SX6
#  if defined ASL
       ALLOCATE (cw2(ldx*ldy*ldz))
       CALL zfc3cl (f(1), nx, ny, nz, ldx, ldy, ldz, err)
#  else
       ALLOCATE (cw2(6*ldx*ldy*ldz))
#  endif
#endif
     !
     !   Here initialize table only if necessary
     !
     ip = -1
     DO i = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       IF ( ( nx == dims(1,i) ) .and. &
            ( ny == dims(2,i) ) .and. &
            ( nz == dims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF
     END DO

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one


#if defined __FFTW
       IF ( nx /= ldx .or. ny /= ldy .or. nz /= ldz ) &
         call errore('cfft3','not implemented',1)

       IF( fw_plan(icurrent) /= 0 ) CALL DESTROY_PLAN_3D( fw_plan(icurrent) )
       IF( bw_plan(icurrent) /= 0 ) CALL DESTROY_PLAN_3D( bw_plan(icurrent) )
       idir = -1; CALL CREATE_PLAN_3D( fw_plan(icurrent), nx, ny, nz, idir)
       idir =  1; CALL CREATE_PLAN_3D( bw_plan(icurrent), nx, ny, nz, idir)

#elif defined __FFTW3

       IF ( nx /= ldx .or. ny /= ldy .or. nz /= ldz ) &
            call errore('cfft3','not implemented',3)
       IF( fw_plan(icurrent) /= 0 ) CALL dfftw_destroy_plan( fw_plan(icurrent) )
       IF( bw_plan(icurrent) /= 0 ) CALL dfftw_destroy_plan( bw_plan(icurrent) )
       idir = -1
       CALL dfftw_plan_dft_3d ( fw_plan(icurrent), nx, ny, nz, f(1:), &
            f(1:), idir, FFTW_ESTIMATE)
       idir =  1
       CALL dfftw_plan_dft_3d ( bw_plan(icurrent), nx, ny, nz, f(1:), &
            f(1:), idir, FFTW_ESTIMATE)

#elif defined __DFTI

      if( ASSOCIATED( hand(icurrent)%desc ) ) THEN
          dfti_status = DftiFreeDescriptor( hand(icurrent)%desc )
          IF( dfti_status /= 0) THEN
             WRITE(*,*) "stopped in DftiFreeDescriptor", dfti_status
             STOP
          ENDIF
       END IF

       dfti_status = DftiCreateDescriptor(hand(icurrent)%desc, DFTI_DOUBLE, DFTI_COMPLEX, 3,(/nx,ny,nz/))
       IF(dfti_status /= 0) THEN
          WRITE(*,*) "stopped in DftiCreateDescriptor", dfti_status
          STOP
       ENDIF
       dfti_status = DftiSetValue(hand(icurrent)%desc, DFTI_NUMBER_OF_TRANSFORMS,1)
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DFTI_NUMBER_OF_TRANSFORMS", dfti_status
          STOP
       ENDIF
       dfti_status = DftiSetValue(hand(icurrent)%desc, DFTI_PLACEMENT, DFTI_INPLACE)
       IF(dfti_status /= 0)THEN
         WRITE(*,*) "stopped in DFTI_PLACEMENT", dfti_status
         STOP
      ENDIF
       tscale = 1.0_DP/ (nx * ny * nz)
       dfti_status = DftiSetValue( hand(icurrent)%desc, DFTI_FORWARD_SCALE, tscale);
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DFTI_FORWARD_SCALE", dfti_status
          STOP
       ENDIF
       tscale = 1.0_DP
       dfti_status = DftiSetValue( hand(icurrent)%desc, DFTI_BACKWARD_SCALE, tscale );
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DFTI_BACKWARD_SCALE", dfti_status
          STOP
       ENDIF

       dfti_status = DftiCommitDescriptor(hand(icurrent)%desc)
       IF(dfti_status /= 0) THEN
          WRITE(*,*) "stopped in DftiCreateDescriptor", dfti_status
          STOP
       ENDIF




#elif defined __ESSL || defined __LINUX_ESSL

       ! no initialization for 3d FFT's from ESSL

#elif defined __SX6

#  if defined ASL
#    if defined MICRO
       CALL hfc3fb (nx,ny,nz, f(1) , ldx, ldy, ldz, 0, &
            iw0(1,icurrent), auxp(1,icurrent), cw2(1), nbtasks, err)
#    else
       CALL zfc3fb (nx,ny,nz, f(1), ldx, ldy, ldz, 0, &
             iw0(1,icurrent), auxp(1,icurrent), cw2(1), err)
#    endif
#  else
       ! for some reason the error variable is not set by this driver on NEC SX machines
       err = 0 
       CALL ZZFFT3D (0, nx,ny,nz, 1.0_DP, f(1), ldx, ldy, &
          &             f(1), ldx, ldy, auxp(1,icurrent), cw2(1), err)
#  endif

       IF (err /= 0) CALL errore('cfft3d','FFT init returned an error ', err)

#else

       CALL errore(' cfft3d ',' no scalar fft driver specified ', 1)

#endif

       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF

     !
     !   Now perform the 3D FFT using the machine specific driver
     !

#if defined __FFTW
     IF( isign < 0 ) THEN
       call FFTW_INPLACE_DRV_3D( fw_plan(ip), 1, f(1), 1, 1 )
       tscale = 1.0_DP / DBLE( nx * ny * nz )
       call ZDSCAL( nx * ny * nz, tscale, f(1), 1)

     ELSE IF( isign > 0 ) THEN
       call FFTW_INPLACE_DRV_3D( bw_plan(ip), 1, f(1), 1, 1 )
     END IF

#elif defined __FFTW3

   IF( isign < 0 ) THEN
      call dfftw_execute_dft( fw_plan(ip), f(1:), f(1:))
      tscale = 1.0_DP / DBLE( nx * ny * nz )
      call ZDSCAL( nx * ny * nz, tscale, f(1), 1)

   ELSE IF( isign > 0 ) THEN

      call dfftw_execute_dft( bw_plan(ip), f(1:), f(1:))

   END IF

#elif defined __DFTI

     IF( isign < 0 ) THEN
        !
        dfti_status = DftiComputeForward(hand(ip)%desc, f(1:))
        IF(dfti_status /= 0)THEN
           WRITE(*,*) "stopped in DftiComputeForward", dfti_status
           STOP
        ENDIF
        !
     ELSE IF( isign > 0 ) THEN
        !
        dfti_status = DftiComputeBackward(hand(ip)%desc, f(1:))
        IF(dfti_status /= 0)THEN
           WRITE(*,*) "stopped in DftiComputeBackward", dfti_status
           STOP
        ENDIF
        !
     END IF




#elif defined __ESSL || defined __LINUX_ESSL

     IF ( isign < 0 ) THEN
       tscale = 1.0_DP / ( nx * ny * nz )
       idir = +1
     ELSE IF( isign > 0 ) THEN
       tscale = 1.0_DP
       idir = -1
     END IF

     IF( isign /= 0 ) CALL dcft3( f(1), ldx,ldx*ldy, f(1), ldx,ldx*ldy, &
          nx,ny,nz, idir, tscale, work(1), lwork)

#elif defined __SX6

#  if defined ASL
#    if defined MICRO
     CALL hfc3bf (nx,ny,nz, f(1), ldx,ldy, ldz, &
          -isign, iw0(1,ip), auxp(1,ip), cw2(1), nbtasks, err)
#    else
     CALL zfc3bf (nx,ny,nz, f(1), ldx,ldy, ldz, &
          -isign, iw0(1,ip), auxp(1,ip), cw2(1), err)
#    endif
     IF ( isign < 0) THEN
        tscale = 1.0_DP / DBLE( nx * ny * nz )
        call ZDSCAL( ldx * ldy * ldz, tscale, f(1), 1)
     END IF
#  else
     ! for some reason the error variable is not set by this driver on NEC SX machines
     err = 0 
     tscale = 1.0_DP
     IF ( isign < 0) THEN
        tscale = tscale / DBLE( nx * ny * nz )
     END IF
     CALL ZZFFT3D (isign, nx,ny,nz, tscale, f(1), ldx,ldy, &
          f_out(1), ldx,ldy, auxp(1,ip), cw2(1), err)
!$omp parallel do private(j,i,k_off,kj_offset)
     do k=1,nz
        k_off = (k-1)*ldx*ldy
        do j=1,ny
           kj_offset = (j-1)*ldx + k_off
           do i=1,nx
              f(i+kj_offset) = f_out(i+kj_offset)
           end do
        end do
     end do
!$omp end parallel do
#   endif
     IF (err /= 0) CALL errore('cfft3d','FFT returned an error ', err)
     DEALLOCATE(cw2)

#endif

     RETURN
   END SUBROUTINE cfft3d

!
!=----------------------------------------------------------------------=!
!
!
!
!         3D scalar FFTs,  but using sticks!
!
!
!
!=----------------------------------------------------------------------=!
!

SUBROUTINE cfft3ds (f, nx, ny, nz, ldx, ldy, ldz, isign, &
     do_fft_x, do_fft_y)
  !
  !     driver routine for 3d complex "reduced" fft - see cfft3d
  !     The 3D fft are computed only on lines and planes which have
  !     non zero elements. These lines and planes are defined by
  !     the two integer vectors do_fft_x(ldy*nz) and do_fft_y(nz)
  !     (1 = perform fft, 0 = do not perform fft)
  !     This routine is implemented only for fftw, essl, acml
  !     If not implemented, cfft3d is called instead
  !
  !----------------------------------------------------------------------
  !
  implicit none

  integer :: nx, ny, nz, ldx, ldy, ldz, isign
  !
  !   logical dimensions of the fft
  !   physical dimensions of the f array
  !   sign of the transformation

  complex(DP) :: f ( ldx * ldy * ldz )
  integer :: do_fft_x(:), do_fft_y(:)
  !
  integer :: m, incx1, incx2
  INTEGER :: i, k, j, err, idir, ip,  ii, jj
  REAL(DP) :: tscale
  INTEGER, SAVE :: icurrent = 1
  INTEGER, SAVE :: dims(3,ndims) = -1

#if defined __FFTW || __FFTW3

  C_POINTER, SAVE :: fw_plan ( 3, ndims ) = 0
  C_POINTER, SAVE :: bw_plan ( 3, ndims ) = 0

#elif defined __ESSL || defined __LINUX_ESSL

  INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx
  REAL (DP), SAVE :: fw_table( ltabl, 3, ndims )
  REAL (DP), SAVE :: bw_table( ltabl, 3, ndims )

#else

  CALL cfft3d (f, nx, ny, nz, ldx, ldy, ldz, isign)
  RETURN

#endif

  tscale = 1.0_DP

  ! WRITE( stdout, fmt="('DEBUG cfft3ds :',6I6)") nx, ny, nz, ldx, ldy, ldz
  ! WRITE( stdout, fmt="('DEBUG cfft3ds :',24I2)") do_fft_x
  ! WRITE( stdout, fmt="('DEBUG cfft3ds :',24I2)") do_fft_y


  IF( ny /= ldy ) &
    CALL errore(' cfft3ds ', ' wrong dimensions: ny /= ldy ', 1 )

     ip = -1
     DO i = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       IF( ( nx == dims(1,i) ) .and. ( ny == dims(2,i) ) .and. &
           ( nz == dims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF

     END DO

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

#if defined __FFTW

       IF( fw_plan( 1, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_plan( 1, icurrent) )
       IF( bw_plan( 1, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_plan( 1, icurrent) )
       IF( fw_plan( 2, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_plan( 2, icurrent) )
       IF( bw_plan( 2, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_plan( 2, icurrent) )
       IF( fw_plan( 3, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_plan( 3, icurrent) )
       IF( bw_plan( 3, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_plan( 3, icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 1, icurrent), nx, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 1, icurrent), nx, idir)
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 2, icurrent), ny, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 2, icurrent), ny, idir)
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 3, icurrent), nz, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 3, icurrent), nz, idir)

#elif defined __FFTW3
       IF( fw_plan( 1, icurrent) /= 0 ) &
            CALL dfftw_destroy_plan( fw_plan( 1, icurrent) )
       IF( bw_plan( 1, icurrent) /= 0 ) &
            CALL dfftw_destroy_plan( bw_plan( 1, icurrent) )
       IF( fw_plan( 2, icurrent) /= 0 ) &
            CALL dfftw_destroy_plan( fw_plan( 2, icurrent) )
       IF( bw_plan( 2, icurrent) /= 0 ) &
            CALL dfftw_destroy_plan( bw_plan( 2, icurrent) )
       IF( fw_plan( 3, icurrent) /= 0 ) &
            CALL dfftw_destroy_plan( fw_plan( 3, icurrent) )
       IF( bw_plan( 3, icurrent) /= 0 ) &
            CALL dfftw_destroy_plan( bw_plan( 3, icurrent) )
       idir = -1
       CALL dfftw_plan_many_dft( fw_plan( 1, icurrent), &
            1, nx, 1, f(1:), (/ldx, ldy, ldz/), 1, ldx, &
            f(1:), (/ldx, ldy, ldz/), 1, ldx, idir, FFTW_ESTIMATE)
       idir = 1
       CALL dfftw_plan_many_dft( bw_plan( 1, icurrent), &
            1, nx, 1, f(1:), (/ldx, ldy, ldz/), 1, ldx, &
            f(1:), (/ldx, ldy, ldz/), 1, ldx, idir, FFTW_ESTIMATE)
       idir = -1
       CALL dfftw_plan_many_dft( fw_plan( 2, icurrent), &
            1, ny, nx, f(1:), (/ldx, ldy, ldz/), ldx, 1, &
            f(1:), (/ldx, ldy, ldz/), ldx, 1, idir, FFTW_ESTIMATE)
       idir = 1
       CALL dfftw_plan_many_dft( bw_plan( 2, icurrent), &
            1, ny, nx, f(1:), (/ldx, ldy, ldz/), ldx, 1, &
            f(1:), (/ldx, ldy, ldz/), ldx, 1, idir, FFTW_ESTIMATE)
       idir = -1
       CALL dfftw_plan_many_dft( fw_plan( 3, icurrent), &
            1, nz, nx*ny, f(1:), (/ldx, ldy, ldz/), ldx*ldy, 1, &
            f(1:), (/ldx, ldy, ldz/), ldx*ldy, 1, idir, FFTW_ESTIMATE)
       idir = 1
       CALL dfftw_plan_many_dft( bw_plan( 3, icurrent), &
            1, nz, nx*ny, f(1:), (/ldx, ldy, ldz/), ldx*ldy, 1, &
            f(1:), (/ldx, ldy, ldz/), ldx*ldy, 1, idir, FFTW_ESTIMATE)

#elif defined __ESSL || defined __LINUX_ESSL
       !
       ! ESSL sign convention for fft's is the opposite of the "usual" one
       !
       tscale = 1.0_DP
       !  x - direction
       incx1 = 1; incx2 = ldx; m = 1
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nx, m,  1, 1.0_DP, &
          fw_table( 1, 1, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nx, m, -1, 1.0_DP, &
          bw_table(1, 1, icurrent), ltabl, work(1), lwork )
       !  y - direction
       incx1 = ldx; incx2 = 1; m = nx;
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, ny, m,  1, 1.0_DP, &
          fw_table( 1, 2, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, ny, m, -1, 1.0_DP, &
          bw_table(1, 2, icurrent), ltabl, work(1), lwork )
       !  z - direction
       incx1 = ldx * ldy; incx2 = 1; m = ldx * ny
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nz, m,  1, 1.0_DP, &
          fw_table(1, 3, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nz, m, -1, 1.0_DP, &
          bw_table(1, 3, icurrent), ltabl, work(1), lwork )

#else

       CALL errore(' cfft3ds ',' no scalar fft driver specified ', 1)

#endif

       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF


     IF ( isign > 0 ) THEN

        !
        !  i - direction ...
        !

        incx1 = 1;  incx2 = ldx;  m = 1

        do k = 1, nz
           do j = 1, ny
              jj = j + ( k - 1 ) * ldy
              ii = 1 + ldx * ( jj - 1 )
              if ( do_fft_x( jj ) == 1 ) THEN
#if defined __FFTW
                call FFTW_INPLACE_DRV_1D( bw_plan( 1, ip), m, f( ii ), incx1, incx2 )
#elif defined __FFTW3
                call dfftw_execute_dft( bw_plan( 1, ip), f( ii: ), f( ii: ) )
#elif defined __ESSL || defined __LINUX_ESSL
                call dcft (0, f (ii), incx1,incx2, f (ii), incx1,incx2, nx, m, &
                -isign, 1.0_DP, bw_table ( 1, 1,  ip ), ltabl, work( 1 ), lwork)
#else
                call errore(' cfft3ds ',' no scalar fft driver specified ', 2)
#endif
              endif
           enddo
        enddo

        !
        !  ... j-direction ...
        !

        incx1 = ldx;  incx2 = 1;  m = nx

        do k = 1, nz
           ii = 1 + ldx * ldy * ( k - 1 )
           if ( do_fft_y( k ) == 1 ) then
#if defined __FFTW
             call FFTW_INPLACE_DRV_1D( bw_plan( 2, ip), m, f( ii ), incx1, incx2 )
#elif defined __FFTW3
             call dfftw_execute_dft( bw_plan( 2, ip), f( ii: ), f( ii: ) )
#elif defined __ESSL || defined __LINUX_ESSL
             call dcft (0, f (ii), incx1, incx2, f (ii), incx1, incx2, nx, m, &
               -isign, 1.0_DP, bw_table ( 1, 2,  ip ), ltabl, work( 1 ), lwork)
#else
             call errore(' cfft3ds ',' no scalar fft driver specified ', 3)
#endif
           endif
        enddo

        !
        !     ... k-direction
        !

        incx1 = ldx * ldy;  incx2 = 1;  m = ldx * ny

#if defined __FFTW
        call FFTW_INPLACE_DRV_1D( bw_plan( 3, ip), m, f( 1 ), incx1, incx2 )
#elif defined __FFTW3
        call dfftw_execute_dft( bw_plan( 3, ip), f(1:), f(1:) )
#elif defined __ESSL || defined __LINUX_ESSL
        call dcft (0, f( 1 ), incx1, incx2, f( 1 ), incx1, incx2, nz, m, &
          -isign, 1.0_DP, bw_table ( 1, 3, ip ), ltabl, work( 1 ), lwork)
#endif

     ELSE

        !
        !     ... k-direction
        !

        incx1 = ldx * ny;  incx2 = 1;  m = ldx * ny

#if defined __FFTW
        call FFTW_INPLACE_DRV_1D( fw_plan( 3, ip), m, f( 1 ), incx1, incx2 )
#elif defined __FFTW3
        call dfftw_execute_dft( fw_plan( 3, ip), f(1:), f(1:) )
#elif defined __ESSL || defined __LINUX_ESSL
         call dcft (0, f( 1 ), incx1, incx2, f( 1 ), incx1, incx2, nz, m, &
          -isign, 1.0_DP, fw_table ( 1, 3, ip ), ltabl, work( 1 ), lwork)

#endif

        !
        !     ... j-direction ...
        !

        incx1 = ldx;  incx2 = 1;  m = nx

        do k = 1, nz
           ii = 1 + ldx * ldy * ( k - 1 )
           if ( do_fft_y ( k ) == 1 ) then
#if defined __FFTW
             call FFTW_INPLACE_DRV_1D( fw_plan( 2, ip), m, f( ii ), incx1, incx2 )
#elif defined __FFTW3
             call dfftw_execute_dft( fw_plan( 2, ip), f( ii: ), f( ii: ) )
#elif defined __ESSL || defined __LINUX_ESSL
             call dcft (0, f (ii), incx1, incx2, f (ii), incx1, incx2, ny, m, &
               -isign, 1.0_DP, fw_table ( 1, 2, ip ), ltabl, work( 1 ), lwork)
#else
             call errore(' cfft3ds ',' no scalar fft driver specified ', 4)
#endif
           endif
        enddo

        !
        !     i - direction ...
        !

        incx1 = 1;  incx2 = ldx;  m = 1

        do k = 1, nz
           do j = 1, ny
              jj = j + ( k - 1 ) * ldy
              ii = 1 + ldx * ( jj - 1 )
              if ( do_fft_x( jj ) == 1 ) then
#if defined __FFTW
                call FFTW_INPLACE_DRV_1D( fw_plan( 1, ip), m, f( ii ), incx1, incx2 )
#elif defined __FFTW3
                call dfftw_execute_dft( fw_plan( 1, ip), f( ii: ), f( ii: ) )
#elif defined __ESSL || defined __LINUX_ESSL
                call dcft (0, f (ii), incx1,incx2, f (ii), incx1,incx2, nx, m, &
                 -isign, 1.0_DP, fw_table ( 1, 1, ip ), ltabl, work( 1 ), lwork)
#else
                call errore(' cfft3ds ',' no scalar fft driver specified ', 5)
#endif
              endif
           enddo
        enddo

        call DSCAL (2 * ldx * ldy * nz, 1.0_DP/(nx * ny * nz), f(1), 1)

     END IF
     RETURN
   END SUBROUTINE cfft3ds

!
!=----------------------------------------------------------------------=!
!
!
!
!         3D parallel FFT on sub-grids
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cft_b ( f, nx, ny, nz, ldx, ldy, ldz, imin3, imax3, sgn )

!     driver routine for 3d complex fft's on box grid, parallel case
!     fft along xy is done only on planes that correspond to dense grid
!     planes on the current processor, i.e. planes with imin3 <= nz <= imax3
!     implemented for essl, fftw, scsl, complib, only for sgn=1 (f(R) => f(G))
!     (beware: here the "essl" convention for the sign of the fft is used!)
!
      implicit none
      integer nx,ny,nz,ldx,ldy,ldz,imin3,imax3,sgn
      complex(dp) :: f(:)

      integer isign, naux, ibid, nplanes, nstart, k
      real(DP) :: tscale

      integer :: ip, i
      integer, save :: icurrent = 1
      integer, save :: dims( 4, ndims ) = -1

#if defined __FFTW || __FFTW3 || __DFTI

      C_POINTER, save :: bw_planz(  ndims ) = 0
      C_POINTER, save :: bw_planx(  ndims ) = 0
      C_POINTER, save :: bw_plany(  ndims ) = 0
      C_POINTER, save :: bw_planxy( ndims ) = 0

#elif defined __ESSL || defined __LINUX_ESSL

      INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx
      real(dp), save :: aux3( ltabl, ndims )
      real(dp), save :: aux2( ltabl, ndims )
      real(dp), save :: aux1( ltabl, ndims )

#endif

      isign = -sgn
      tscale = 1.0_DP

      if ( isign > 0 ) then
         call errore('cft_b','not implemented',isign)
      end if
!
! 2d fft on xy planes - only needed planes are transformed
! note that all others are left in an unusable state
!
      nplanes = imax3 - imin3 + 1
      nstart  = ( imin3 - 1 ) * ldx * ldy + 1

      !
      !   Here initialize table only if necessary
      !

      ip = -1
      DO i = 1, ndims

        !   first check if there is already a table initialized
        !   for this combination of parameters

        IF ( ( nx == dims(1,i) ) .and. ( ny == dims(2,i) ) .and. &
             ( nz == dims(3,i) ) .and. ( nplanes == dims(4,i) ) ) THEN
           ip = i
           EXIT
        END IF

      END DO

      IF( ip == -1 ) THEN

        !   no table exist for these parameters
        !   initialize a new one

#if defined __FFTW

        if ( bw_planz(icurrent) /= 0 ) &
             call DESTROY_PLAN_1D( bw_planz(icurrent) )
        call CREATE_PLAN_1D( bw_planz(icurrent), nz, 1 )

        if ( bw_planx(icurrent) /= 0 ) &
             call DESTROY_PLAN_1D( bw_planx(icurrent) )
        call CREATE_PLAN_1D( bw_planx(icurrent), nx, 1 )

        if ( bw_plany(icurrent) /= 0 ) &
             call DESTROY_PLAN_1D( bw_plany(icurrent) )
        call CREATE_PLAN_1D( bw_plany(icurrent), ny, 1 )

        if ( bw_planxy(icurrent) /= 0 ) &
             call DESTROY_PLAN_2D( bw_planxy(icurrent) )
        call CREATE_PLAN_2D( bw_planxy(icurrent), nx, ny, 1 )
!
#elif defined __FFTW3

        if ( bw_planz(icurrent) /= 0 ) &
             call dfftw_destroy_plan(bw_planz(icurrent))
        call dfftw_plan_many_dft( bw_planz(icurrent), 1, nz, ldx*ldy, &
             f(1:), (/SIZE(f)/), ldx*ldy, 1, f(1:), (/SIZE(f)/), ldx*ldy, 1, &
             1, FFTW_ESTIMATE )

        if ( bw_planxy(icurrent) /= 0 ) &
             call dfftw_destroy_plan(bw_planxy(icurrent))
        call dfftw_plan_many_dft( bw_planxy(icurrent), 2, (/nx, ny/), nplanes,&
             f(nstart:),  (/ldx, ldy/), 1, ldx*ldy, f(nstart:), (/ldx, ldy/), &
             1, ldx*ldy, 1, FFTW_ESTIMATE )

#elif defined __ESSL || defined __LINUX_ESSL

         if( nz /= dims(3,icurrent) ) then
           call dcft( 1, f(1), ldx*ldy, 1, f(1), ldx*ldy, 1, nz, ldx*ldy, &
              isign, tscale, aux3(1,icurrent), ltabl, work(1), lwork)
         end if
         call dcft( 1, f(1), 1, ldx, f(1), 1, ldx, nx, ldy*nplanes, isign, &
              tscale, aux1(1,icurrent), ltabl, work(1), lwork)
         if( ny /= dims(2,icurrent) ) then
           call dcft( 1, f(1), ldx, 1, f(1), ldx, 1, ny, ldx, isign,    &
              tscale, aux2(1,icurrent), ltabl, work(1), lwork)
         end if

#else

        CALL errore(' cft_b ',' no scalar fft driver specified ', 1)

#endif

        dims(1,icurrent) = nx; dims(2,icurrent) = ny
        dims(3,icurrent) = nz; dims(4,icurrent) = nplanes
        ip = icurrent
        icurrent = MOD( icurrent, ndims ) + 1

      END IF


#if defined __FFTW

      !
      !  fft along Z
      !
      call FFTW_INPLACE_DRV_1D( bw_planz(ip), ldx*ldy, f(1), ldx*ldy, 1 )
      !
      !  fft along Y
      !  fft along X
      !
      do k = imin3, imax3
        call FFTW_INPLACE_DRV_1D( bw_plany(ip), nx, f((k-1)*ldx*ldy + 1), ldx, 1 )
        call FFTW_INPLACE_DRV_1D( bw_planx(ip), ny, f((k-1)*ldx*ldy + 1), 1, ldx )
      end do   

#elif defined __FFTW3 || defined __DFTI

      call dfftw_execute_dft(bw_planz(ip), f(1:), f(1:))
      call dfftw_execute_dft(bw_planxy(ip), f(nstart:), f(nstart:))

#elif defined __ESSL || defined __LINUX_ESSL

      !   fft in the z-direction...

      call dcft( 0, f(1), ldx*ldy, 1, f(1), ldx*ldy, 1, nz, ldx*ldy, isign, &
           tscale, aux3(1,ip), ltabl, work(1), lwork)

      !   x-direction

      call dcft( 0, f(nstart), 1, ldx, f(nstart), 1, ldx, nx, ldy*nplanes, &
           isign, tscale, aux1(1,ip), ltabl, work(1), lwork)

      !   y-direction

      DO K = imin3, imax3
        nstart = ( k - 1 ) * ldx * ldy + 1
        call dcft( 0, f(nstart), ldx, 1, f(nstart), ldx, 1, ny, ldx, isign, &
             tscale, aux2(1,ip), ltabl, work(1), lwork)
      END DO

#endif

      RETURN
   END SUBROUTINE cft_b

!
!=----------------------------------------------------------------------=!
!
!
!
!   3D parallel FFT on sub-grids, to be called inside OpenMP region
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cft_b_omp_init ( nx, ny, nz )

!     driver routine for 3d complex fft's on box grid, init subroutine
!
      implicit none
      integer, INTENT(IN) :: nx,ny,nz
      !
      !   Here initialize table 
      !
#if defined __FFTW

!$omp parallel

      IF( cft_b_bw_planz  == 0 ) THEN
         CALL CREATE_PLAN_1D( cft_b_bw_planz, nz, 1 )
         cft_b_dims(3) = nz
      END IF
      IF( cft_b_bw_planx  == 0 ) THEN
         CALL CREATE_PLAN_1D( cft_b_bw_planx, nx, 1 )
         cft_b_dims(1) = nx
      END IF
      IF( cft_b_bw_plany  == 0 ) THEN
         CALL CREATE_PLAN_1D( cft_b_bw_plany, ny, 1 )
         cft_b_dims(2) = ny
      END IF

!$omp end parallel

#else

      CALL errore(' cft_b_omp_init ',' no scalar fft driver specified ', 1)

#endif

     RETURN
   END SUBROUTINE cft_b_omp_init


   SUBROUTINE cft_b_omp ( f, nx, ny, nz, ldx, ldy, ldz, imin3, imax3, sgn )

!     driver routine for 3d complex fft's on box grid, parallel (MPI+OpenMP) case
!     fft along xy is done only on planes that correspond to dense grid
!     planes on the current processor, i.e. planes with imin3 <= nz <= imax3
!     implemented ONLY for internal fftw, and only for sgn=1 (f(R) => f(G))
!     (beware: here the "essl" convention for the sign of the fft is used!)
!
!     This driver is meant for calls inside parallel OpenMP sections
!
      implicit none
      integer, INTENT(IN) :: nx,ny,nz,ldx,ldy,ldz,imin3,imax3,sgn
      complex(dp) :: f(:)

      INTEGER, SAVE :: k
!$omp threadprivate (k)

      if ( -sgn > 0 ) then
         CALL errore('cft_b_omp','forward transform not implemented',1)
      end if

#if defined __FFTW

      IF ( ( cft_b_bw_planz == 0 ) .or. ( cft_b_bw_planx == 0 ) .or. ( cft_b_bw_plany == 0 ) ) THEN
         CALL errore('cft_b_omp','plan not initialized',1)
      END IF

      !  consistency check

      IF ( ( nx /= cft_b_dims(1) ) .or. ( ny /= cft_b_dims(2) ) .or. ( nz /= cft_b_dims(3) ) ) THEN
         CALL errore('cft_b_omp', 'dimensions are inconsistent with the existing plan',1) 
      END IF

      !  fft along Z
      !
      call FFTW_INPLACE_DRV_1D( cft_b_bw_planz, ldx*ldy, f(1), ldx*ldy, 1 )
      !
      !  fft along Y
      !  fft along X
      !
      do k = imin3, imax3
        call FFTW_INPLACE_DRV_1D( cft_b_bw_plany, nx, f((k-1)*ldx*ldy + 1), ldx, 1 )
        call FFTW_INPLACE_DRV_1D( cft_b_bw_planx, ny, f((k-1)*ldx*ldy + 1), 1, ldx )
      end do   

#else

      CALL errore(' cft_b_omp ',' no scalar fft driver specified ', 1)

#endif

     RETURN
   END SUBROUTINE cft_b_omp


!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT support Functions/Subroutines
!
!
!
!=----------------------------------------------------------------------=!
!

!
integer function good_fft_dimension (n)
  !
  ! Determines the optimal maximum dimensions of fft arrays
  ! Useful on some machines to avoid memory conflicts
  !
  USE kinds, only : DP
  IMPLICIT NONE
  INTEGER :: n, nx
  REAL(DP) :: log2n
  !
  ! this is the default: max dimension = fft dimension
  nx = n
  !
#if defined(__ESSL) || defined(__LINUX_ESSL)
  log2n = LOG ( dble (n) ) / LOG ( 2.0_DP )
  ! log2n is the logarithm of n in base 2
  IF ( ABS (NINT(log2n) - log2n) < 1.0d-8 ) nx = n + 1
  ! if n is a power of 2 (log2n is integer) increase dimension by 1
#elif defined(__SX6)
  !
  if (mod (n, 2) ==0) nx = n + 1
  ! for nec vector machines: if n is even increase dimension by 1
  !
#endif
  !
  good_fft_dimension = nx
  return
end function good_fft_dimension


!=----------------------------------------------------------------------=!

function allowed (nr)


  ! find if the fft dimension is a good one
  ! a "bad one" is either not implemented (as on IBM with ESSL)
  ! or implemented but with awful performances (most other cases)

  USE kinds

  implicit none
  integer :: nr

  logical :: allowed
  integer :: pwr (5)
  integer :: mr, i, fac, p, maxpwr
  integer :: factors( 5 ) = (/ 2, 3, 5, 7, 11 /)

  ! find the factors of the fft dimension

  mr  = nr
  pwr = 0
  factors_loop: do i = 1, 5
     fac = factors (i)
     maxpwr = NINT ( LOG( DBLE (mr) ) / LOG( DBLE (fac) ) ) + 1
     do p = 1, maxpwr
        if ( mr == 1 ) EXIT factors_loop
        if ( MOD (mr, fac) == 0 ) then
           mr = mr / fac
           pwr (i) = pwr (i) + 1
        endif
     enddo
  end do factors_loop

  IF ( nr /= ( mr * 2**pwr (1) * 3**pwr (2) * 5**pwr (3) * 7**pwr (4) * 11**pwr (5) ) ) &
     CALL errore (' allowed ', ' what ?!? ', 1 )

  if ( mr /= 1 ) then

     ! fft dimension contains factors > 11 : no good in any case

     allowed = .false.

  else

#if defined __ESSL || defined __LINUX_ESSL

     ! IBM machines with essl libraries

     allowed =  ( pwr(1) >= 1 ) .and. ( pwr(2) <= 2 ) .and. ( pwr(3) <= 1 ) .and. &
                ( pwr(4) <= 1 ) .and. ( pwr(5) <= 1 ) .and. &
                ( ( (pwr(2) == 0 ) .and. ( pwr(3) + pwr(4) + pwr(5) ) <= 2 ) .or. &
                  ( (pwr(2) /= 0 ) .and. ( pwr(3) + pwr(4) + pwr(5) ) <= 1 ) )
#else

     ! fftw and all other cases: no factors 7 and 11

     allowed = ( ( pwr(4) == 0 ) .and. ( pwr(5) == 0 ) )

#endif

  endif

  return
end function allowed

!=----------------------------------------------------------------------=!

   INTEGER FUNCTION good_fft_order( nr, np )

!
!    This function find a "good" fft order value greater or equal to "nr"
!
!    nr  (input) tentative order n of a fft
!
!    np  (optional input) if present restrict the search of the order
!        in the ensemble of multiples of np
!
!    Output: the same if n is a good number
!         the closest higher number that is good
!         an fft order is not good if not implemented (as on IBM with ESSL)
!         or implemented but with awful performances (most other cases)
!

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: nr
     INTEGER, OPTIONAL, INTENT(IN) :: np
     INTEGER :: new

     new = nr
     IF( PRESENT( np ) ) THEN
       DO WHILE( ( ( .NOT. allowed( new ) ) .OR. ( MOD( new, np ) /= 0 ) ) .AND. ( new <= nfftx ) )
         new = new + 1
       END DO
     ELSE
       DO WHILE( ( .NOT. allowed( new ) ) .AND. ( new <= nfftx ) )
         new = new + 1
       END DO
     END IF

     IF( new > nfftx ) &
       CALL errore( ' good_fft_order ', ' fft order too large ', new )

     good_fft_order = new

     RETURN
   END FUNCTION good_fft_order


!=----------------------------------------------------------------------=!
   END MODULE fft_scalar
!=----------------------------------------------------------------------=!
