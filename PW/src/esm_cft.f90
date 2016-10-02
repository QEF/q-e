!
! Copyright (C) 2007-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------

#include "../../include/fft_defs.h"

MODULE esm_cft
  !--------------------------------------------------------------------------
  !
  ! ... this module contains the variables and SUBROUTINEs needed for the 
  ! ... 1-Dimentinal FFT called from ESM routines
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: esm_cft_1z_init, esm_cft_1z
  !
  LOGICAL :: cft_initialized = .false.
  
  !   nfftx   Max allowed fft dimension

  INTEGER, PARAMETER :: nfftx = 2049

#if defined __FFTW || defined __FFTW3

  !   Pointers to the "C" structures containing FFT factors ( PLAN )
  !   C_POINTER is defined in include/fft_defs.h
  !   for 32bit executables, C_POINTER is integer(4)
  !   for 64bit executables, C_POINTER is integer(8)

  C_POINTER :: fw_planz = 0
  C_POINTER :: bw_planz = 0

#elif defined(__LINUX_ESSL)

  !   ESSL IBM library: see the ESSL manual for DCFT

  !   Workspace that is dynamically allocated in initialization is defined here
  !   lwork:   Dimension of the work space array (if any)

  !   ltabl   Dimension of the tables of factors calculated at the
  !           initialization stage
  
  INTEGER, PARAMETER :: lwork = 20000 + ( 2*nfftx + 256 ) * 64 + 3*nfftx
  REAL (DP), ALLOCATABLE :: work( :, : )

  INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx
  REAL (DP) :: fw_tablez( ltabl )
  REAL (DP) :: bw_tablez( ltabl )

#elif defined(__SCSL)

  !   SGI scientific library scsl

  INTEGER, PARAMETER :: lwork = 2 * nfftx
  COMPLEX (DP), ALLOCATABLE :: work( :, : )

  INTEGER, PARAMETER :: ltabl = 2 * nfftx + 256
  REAL (DP) :: tablez (ltabl)
  INTEGER :: isys(0:1) = (/ 1, 1 /)

#elif defined(__SX6)

  !   NEC MathKeisan

  INTEGER :: lwork
  REAL (DP) :: work( :, : )

  INTEGER, PARAMETER :: ltabl = 2 * nfftx + 64
  REAL (DP) :: tablez (ltabl)
  INTEGER :: isys = 1

#elif defined(__SUNPERF)
  
  !   SUN sunperf library

  INTEGER, PARAMETER :: ltabl = 4 * nfftx + 15
  REAL (DP) :: tablez (ltabl)

#endif

#if defined __FFTW3

  !  Only FFTW_ESTIMATE is actually used

#define  FFTW_MEASURE  0
#define  FFTW_ESTIMATE 64

#endif

CONTAINS

SUBROUTINE esm_cft_1z_init(nsl, nz, ldz)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nsl, nz, ldz
  REAL (DP)  :: tscale
  INTEGER    :: idir, nth
#if defined __FFTW3 || defined __LINUX_ESSL
  COMPLEX(DP), ALLOCATABLE :: c(:), cout(:)
#endif
#if defined __OPENMP
  INTEGER, EXTERNAL   :: OMP_GET_MAX_THREADS
#endif

#if defined __SCSL
  
  !   SGI scientific library scsl
  
  REAL (DP)       :: DUMMY
  
#elif defined(__SX6)
  
  !   NEC MathKeisan
  
  COMPLEX (DP)    :: DUMMY
  
#endif

  !   initialize a new one

  nth = 1
#if defined(__OPENMP)
  nth = OMP_GET_MAX_THREADS()
#endif

#if defined __FFTW3 || defined __LINUX_ESSL
  ALLOCATE ( c( ldz*nsl ) )
  ALLOCATE ( cout( ldz*nsl ) )
#endif

#if defined __FFTW

  IF( fw_planz /= 0 ) CALL DESTROY_PLAN_1D( fw_planz )
  IF( bw_planz /= 0 ) CALL DESTROY_PLAN_1D( bw_planz )
  idir = -1; CALL CREATE_PLAN_1D( fw_planz, nz, idir)
  idir =  1; CALL CREATE_PLAN_1D( bw_planz, nz, idir)

#elif defined(__FFTW3)

  IF( fw_planz /= 0 ) CALL dfftw_destroy_plan( fw_planz )
  IF( bw_planz /= 0 ) CALL dfftw_destroy_plan( bw_planz )
  idir = -1
  CALL dfftw_plan_many_dft( fw_planz, 1, nz, nsl, c, &
    (/ldz*nsl/), 1, ldz, cout, (/ldz*nsl/), 1, ldz, idir, FFTW_ESTIMATE)
  idir = 1
  CALL dfftw_plan_many_dft( bw_planz, 1, nz, nsl, c, &
    (/ldz*nsl/), 1, ldz, cout, (/ldz*nsl/), 1, ldz, idir, FFTW_ESTIMATE)
  DEALLOCATE ( c )
  DEALLOCATE ( cout )

#elif defined(__LINUX_ESSL)

  IF( ALLOCATED( work ) ) DEALLOCATE( work )
  ALLOCATE( work( lwork, 0:nth-1 ) )
  tscale = 1.0_DP / nz
  CALL DCFT ( 1, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl,  1, &
    tscale, fw_tablez, ltabl, work(1,1), lwork)
  CALL DCFT ( 1, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, -1, &
    1.0_DP, bw_tablez, ltabl, work(1,1), lwork)

#elif defined(__SCSL)

  IF( ALLOCATED( work ) ) DEALLOCATE( work )
  ALLOCATE( work( lwork, 0:nth-1 ) )
  CALL ZZFFTM (0, nz, 0, 0.0_DP, DUMMY, 1, DUMMY, 1, &
    tablez, DUMMY, isys)

#elif defined(__SX6)
  
  lwork = 4 * nz * nsl
  IF( ALLOCATED( work ) ) DEALLOCATE( work )
  ALLOCATE( work( lwork, 0:nth-1 ) )
  CALL ZZFFTM (0, nz, 1, 1.0_DP, DUMMY, ldz, DUMMY, ldz, &
    tablez, work(1,1), isys)

#elif defined(__SUNPERF)

  CALL zffti (nz, tablez )
  
#else

#if defined __FFTW3 || defined __LINUX_ESSL
  IF( ALLOCATED( c ) ) DEALLOCATE( c )
  IF( ALLOCATED( cout ) ) DEALLOCATE( cout )
#endif

  CALL errore(' esm_cft_1z_init ',' no scalar fft driver specified ', 1)

#endif

  cft_initialized = .true.

  RETURN
END SUBROUTINE esm_cft_1z_init


! esm_cft_1z is modified version of cft_1z in fft_scalar.f90
! esm_cft_1z is supposed to be called from OpenMP parallel directive

SUBROUTINE esm_cft_1z(c, nsl, nz, ldz, isign, cout)

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
  INTEGER    :: i, idir, ith

#if defined __OPENMP
  INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM
#endif

  IF( nsl < 0 ) THEN
    CALL errore(" esm: esm_cft_1z ", " nsl out of range ", nsl)
  END IF

  IF( .NOT. cft_initialized ) THEN
    CALL errore(" esm: esm_cft_1z ", " not initialized ", nsl)
  END IF
  
  ith = 0
#if defined __OPENMP
  ith = OMP_GET_THREAD_NUM()
#endif

  !
  !   Now perform the FFTs using machine specific drivers
  !
  
#if defined __FFT_CLOCKS
  CALL start_clock( 'esm_cft_1z' )
#endif
  
#if defined __FFTW
  
  IF (isign < 0) THEN
    CALL FFT_Z_STICK(fw_planz, c(1), ldz, nsl)
    tscale = 1.0_DP / nz
    cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl ) * tscale
  ELSE IF (isign > 0) THEN
    CALL FFT_Z_STICK(bw_planz, c(1), ldz, nsl)
    cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl )
  END IF
  
#elif defined(__FFTW3)

  IF (isign < 0) THEN
    CALL dfftw_execute_dft( fw_planz, c, cout)
    tscale = 1.0_DP / nz
    cout( 1 : ldz * nsl ) = cout( 1 : ldz * nsl ) * tscale
  ELSE IF (isign > 0) THEN
    CALL dfftw_execute_dft( bw_planz, c, cout)
  END IF
  
#elif defined(__SCSL)
  
  IF ( isign < 0 ) THEN
    idir   = -1
    tscale = 1.0_DP / nz
  ELSE IF ( isign > 0 ) THEN
    idir   = 1
    tscale = 1.0_DP
  END IF
  IF (isign /= 0) CALL ZZFFTM (idir, nz, nsl, tscale, c(1), ldz, &
    cout(1), ldz, tablez, work(1,ith), isys)
  
#elif defined(__SX6)
  IF ( isign < 0 ) THEN
    idir   = -1
    tscale = 1.0_DP / nz
  ELSE IF ( isign > 0 ) THEN
    idir   = 1
    tscale = 1.0_DP
  END IF
  IF (isign /= 0) CALL ZZFFTM (idir, nz, nsl, tscale, c(1), ldz, &
    cout(1), ldz, tablez, work(1,ith), isys)
  
#elif defined(__LINUX_ESSL)
  
  ! essl uses a different convention for forward/backward transforms
  ! wrt most other implementations: notice the sign of "idir"
  
  IF( isign < 0 ) THEN
    idir   =+1
    tscale = 1.0_DP / nz
    CALL DCFT (0, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, idir, &
      tscale, fw_tablez, ltabl, work(1,ith), lwork)
  ELSE IF( isign > 0 ) THEN
    idir   =-1
    tscale = 1.0_DP
    CALL DCFT (0, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, idir, &
      tscale, bw_tablez, ltabl, work(1,ith), lwork)
  END IF
  
#elif defined(__SUNPERF)
  
  IF ( isign < 0) THEN
    DO i = 1, nsl
      CALL zfftf ( nz, c(1+(i-1)*ldz), tablez )
    END DO
    cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl ) / nz
  ELSE IF( isign > 0 ) THEN
    DO i = 1, nsl
      CALL zfftb ( nz, c(1+(i-1)*ldz), tablez )
    enddo
    cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl )
  END IF
  
#else
  
  CALL errore(' esm_cft_1z ',' no scalar fft driver specified ', 1)
  
#endif
  
#if defined __FFT_CLOCKS
  CALL stop_clock( 'esm_cft_1z' )
#endif
  
  RETURN
END SUBROUTINE esm_cft_1z

END MODULE esm_cft



