#if defined __LINUX
#  define FPMD_MESSAGE_MAX_SIZE         262144
#  define FPMD_INTEGER_MESSAGE_MAX_SIZE  65536
#  define FPMD_REAL_MESSAGE_MAX_SIZE     32768
#  define FPMD_COMPLEX_MESSAGE_MAX_SIZE  16384
#else
#  define FPMD_MESSAGE_MAX_SIZE          20000000
#  define FPMD_INTEGER_MESSAGE_MAX_SIZE   5000000
#  define FPMD_REAL_MESSAGE_MAX_SIZE      2500000
#  define FPMD_COMPLEX_MESSAGE_MAX_SIZE   1250000
#endif

#if defined __CRAY 

#  define FPMD_MPI_COMPLEX MPI_COMPLEX
#  define FPMD_MPI_REAL    MPI_REAL

#  define FPMD_DGEMTX   SGEMTX
#  define FPMD_DGEMX   SGEMX
#  define FPMD_ZGEFA   CGEFA
#  define FPMD_ZGEDI   CGEDI
#  define FPMD_DGESUB   SGESUB
#  define FPMD_DGER     SGER
#  define FPMD_DGEMV    SGEMV
#  define FPMD_ZGEMV    CGEMV
#  define FPMD_DCOPY    SCOPY
#  define FPMD_ZCOPY    CCOPY
#  define FPMD_DGEMM    SGEMM
#  define FPMD_ZGEMM    CGEMM
#  define FPMD_DGEMUL    SGEMUL
#  define FPMD_IDAMAX  ISAMAX
#  define FPMD_ZSCAL   CSCAL
#  define FPMD_DSCAL   SSCAL
#  define FPMD_ZDSCAL   CSSCAL
#  define FPMD_DSWAP   SSWAP
#  define FPMD_ZSWAP   CSWAP
#  define FPMD_DYAX    SYAX
#  define FPMD_DNRM2   SNRM2
#  define FPMD_DAXPY   SAXPY
#  define FPMD_ZAXPY   CAXPY
#  define FPMD_ZDOTU   CDOTU
#  define FPMD_ZDOTC   CDOTC
#  define FPMD_DDOT   SDOT

#else

#  define FPMD_MPI_COMPLEX MPI_DOUBLE_COMPLEX
#  define FPMD_MPI_REAL    MPI_DOUBLE_PRECISION

#  define FPMD_DGEMTX   DGEMTX
#  define FPMD_DGEMX   DGEMX
#  define FPMD_ZGEFA   ZGEFA
#  define FPMD_ZGEDI   ZGEDI
#  define FPMD_DGESUB   DGESUB
#  define FPMD_DGER     DGER
#  define FPMD_DGEMV    DGEMV
#  define FPMD_ZGEMV    ZGEMV
#  define FPMD_DCOPY    DCOPY
#  define FPMD_ZCOPY    ZCOPY
#  define FPMD_DGEMM    DGEMM
#  define FPMD_ZGEMM    ZGEMM
#  define FPMD_DGEMUL    DGEMUL
#  define FPMD_IDAMAX  IDAMAX
#  define FPMD_ZSCAL   ZSCAL
#  define FPMD_DSCAL   DSCAL
#  define FPMD_ZDSCAL   ZDSCAL
#  define FPMD_DSWAP   DSWAP
#  define FPMD_ZSWAP   ZSWAP
#  define FPMD_DYAX    DYAX
#  define FPMD_DNRM2   DNRM2
#  define FPMD_DAXPY   DAXPY
#  define FPMD_ZAXPY   ZAXPY
#  define FPMD_ZDOTU   ZDOTU
#  define FPMD_ZDOTC   ZDOTC
#  define FPMD_DDOT   DDOT

#endif

