
#if defined (__SX4) || defined (ORIGIN) || defined (T3D) || defined (FUJ64) || defined (HITACHI)
#  define FLUSH
#endif

#if defined(CRAYY) || defined(SGI) || defined(PC) || defined(AIX) || defined(T3D) ||defined(HITACHI) || defined(SUN)
#  define C_POINTER  integer
#endif

#if defined(ORIGIN) || defined (FUJ64)|| defined (DEC) || defined (__SX6)
#  define C_POINTER  integer*8
#endif

#if defined(__SX4)
#  define C_POINTER  real*8
#endif

#if defined(__SX4)
#  define DIRECT_IO_FACTOR 1
#elif defined(DEC)
#  define DIRECT_IO_FACTOR 2
#else
#  define DIRECT_IO_FACTOR 8 
#endif

#if defined(CRAYY) || defined (__SX4) || defined (T3D)

#  define DREAL       real
#  define DIMAG       aimag
#  define DCMPLX      cmplx

#  define DAXPY       saxpy
#  define DCOPY       scopy
#  define DDOT        sdot
#  define DGETRF      sgetrf
#  define DGETRI      sgetri
#  define DGEMV       sgemv
#  define DGEMM       sgemm
#  define DNRM2       snrm2
#  define DSCAL       sscal
#  define DSYTRF      ssytrf
#  define DSYTRI      ssytri
#  define DSYGV       ssygv
#  define DSYGVX      ssygvx
#  define DSWAP       sswap
#  define ILAENV      ilaenv
#  define ZAXPY       caxpy
#  define ZCOPY       ccopy
#  define ZDOTC       cdotc
#  define ZDOTU       cdotu
#  define ZGEMM       cgemm
#  define ZGEMV       cgemv
#  define ZHEEV       cheev
#  define ZHEGV       chegv
#  define ZHEGVX      chegvx
#  define ZHETRD      CHETRD
#  define ZSCAL       cscal

#  define IZAMAX   icamax
#  define DYAX     syax
#  define ZSWAP    cswap
#  define ZDSCAL   csscal
#  define IDAMAX   isamax
#  define DGEMUL   sgemul
#  define DGESUB   sgesub
#  define DGER     sger
#  define ZGEFA    cgefa
#  define ZGEDI    cgedi
#  define DGEMTX   sgemtx
#  define DGEMX    sgemx


#else

#  define DREAL       dreal
#  define DCMPLX      dcmplx
#  define DIMAG       dimag

#  if defined(ADD_BLAS_TWO_UNDERSCORES)

#    define DAXPY       daxpy__
#    define DCOPY       dcopy__
#    define DDOT        ddot__
#    define DGETRF      dgetrf__
#    define DGETRI      dgetri__
#    define DGEMV       dgemv__
#    define DGEMM       dgemm__
#    define DNRM2       dnrm2__
#    define DSCAL       dscal__
#    define DSYTRF      dsytrf__
#    define DSYTRI      dsytri__
#    define DSYGV       dsygv__
#    define DSYGVX      dsygvx__
#    define DSWAP       dswap__
#    define ILAENV      ilaenv__
#    define ZAXPY       zaxpy__
#    define ZCOPY       zcopy__
#    define ZDOTC       zdotc__
#    define ZDOTU       zdotu__
#    define ZGEMM       zgemm__
#    define ZGEMV       zgemv__
#    define ZHEEV       zheev__
#    define ZHEGV       zhegv__
#    define ZHEGVX      zhegvx__
#    define ZSCAL       zscal__

#    define IZAMAX   izamax__
#    define DYAX     dyax__
#    define ZSWAP    zswap__
#    define ZDSCAL   zdscal__
#    define IDAMAX   idamax__
#    define DGEMUL   dgemul__
#    define DGESUB   dgesub__
#    define DGER     dger__
#    define ZGEFA    zgefa__
#    define ZGEDI    zgedi__
#    define DGEMTX   dgemtx__
#    define DGEMX    dgemx__

#  elif defined(ADD_BLAS_ONE_UNDERSCORE)

#    define DAXPY       daxpy_
#    define DCOPY       dcopy_
#    define DDOT        ddot_
#    define DGETRF      dgetrf_
#    define DGETRI      dgetri_
#    define DGEMV       dgemv_
#    define DGEMM       dgemm_
#    define DNRM2       dnrm2_
#    define DSCAL       dscal_
#    define DSYTRF      dsytrf_
#    define DSYTRI      dsytri_
#    define DSYGV       dsygv_
#    define DSYGVX      dsygvx_
#    define DSWAP       dswap_
#    define ILAENV      ilaenv_
#    define ZAXPY       zaxpy_
#    define ZCOPY       zcopy_
#    define ZDOTC       zdotc_
#    define ZDOTU       zdotu_
#    define ZGEMM       zgemm_
#    define ZGEMV       zgemv_
#    define ZHEEV       zheev_
#    define ZHEGV       zhegv_
#    define ZHEGVX      zhegvx_
#    define ZSCAL       zscal_

#    define IZAMAX   izamax_
#    define DYAX     dyax_
#    define ZSWAP    zswap_
#    define ZDSCAL   zdscal_
#    define IDAMAX   idamax_
#    define DGEMUL   dgemul_
#    define DGESUB   dgesub_
#    define DGER     dger_
#    define ZGEFA    zgefa_
#    define ZGEDI    zgedi_
#    define DGEMTX   dgemtx_
#    define DGEMX    dgemx_

#  else

#    define DAXPY       daxpy
#    define DCOPY       dcopy
#    define DDOT        ddot
#    define DGETRF      dgetrf
#    define DGETRI      dgetri
#    define DGEMV       dgemv
#    define DGEMM       dgemm
#    define DNRM2       dnrm2
#    define DSCAL       dscal
#    define DSYTRF      dsytrf
#    define DSYTRI      dsytri
#    define DSYGV       dsygv
#    define DSYGVX      dsygvx
#    define DSWAP       dswap
#    define ILAENV      ilaenv
#    define ZAXPY       zaxpy
#    define ZCOPY       zcopy
#    define ZDOTC       zdotc
#    define ZDOTU       zdotu
#    define ZGEMM       zgemm
#    define ZGEMV       zgemv
#    define ZHEEV       zheev
#    define ZHEGV       zhegv
#    define ZHEGVX      zhegvx
#    define ZSCAL       zscal

#    define IZAMAX   izamax
#    define DYAX     dyax
#    define ZSWAP    zswap
#    define ZDSCAL   zdscal
#    define IDAMAX   idamax
#    define DGEMUL   dgemul
#    define DGESUB   dgesub
#    define DGER     dger
#    define ZGEFA    zgefa
#    define ZGEDI    zgedi
#    define DGEMTX   dgemtx
#    define DGEMX    dgemx

#  endif

#  define ZHETRD      ZHETRD

#endif

#if defined(ABSOFT)
#  define getenv getenv_
#  define getarg getarg_
#  define iargc  iargc_
#endif

#if defined(T3D_BENCHLIB)
#  define DCOPY       scopy_t3e
#  define DGEMM       s_gemm
#  define ZGEMM       c_gemm
#endif

