
!------------------------------------------------------------------------------!
! CONFIGURATION FILE
!------------------------------------------------------------------------------!
! Here you can choose one of the preset configurations
! Commenting all these lines, you can drive the configuration using -D options.

#ifdef __AIX

@PROCESS OPTIMIZE(0)

#define __IOTK_BINARY_FORMAT "IBM-SP/XLF"
#define __IOTK_LOGICAL1 1
#define __IOTK_LOGICAL2 2
#define __IOTK_LOGICAL3 4
#define __IOTK_LOGICAL4 8
#define __IOTK_INTEGER1 4
#define __IOTK_INTEGER2 2
#define __IOTK_INTEGER3 1
#define __IOTK_INTEGER4 8
#define __IOTK_REAL1    8
#define __IOTK_REAL2    4
#undef  __IOTK_REAL3 
#undef  __IOTK_REAL4

#define __IOTK_MAXRANK 5

#else

#define __IOTK_UNKNOWN_PROCESS
#define __IOTK_BINARY_FORMAT "UNKNOWN"
#define __IOTK_LOGICAL1 1
#define __IOTK_LOGICAL2 2
#define __IOTK_LOGICAL3 4
#define __IOTK_LOGICAL4 8
#define __IOTK_INTEGER1 4
#define __IOTK_INTEGER2 2
#define __IOTK_INTEGER3 1
#define __IOTK_INTEGER4 8
#define __IOTK_REAL1    8
#define __IOTK_REAL2    4
#undef  __IOTK_REAL3 
#undef  __IOTK_REAL4

#define __IOTK_MAXRANK 5

#endif

#ifdef __MPI
#  define __IOTK_MPI_ABORT
#else
#  undef  __IOTK_MPI_ABORT
#endif

!------------------------------------------------------------------------------!

