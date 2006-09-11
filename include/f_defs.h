!
!  Copyright (C) 2002-2004 PWSCF,FPMD,CPV groups
!  This file is distributed under the terms of the
!  GNU General Public License. See the file `License'
!  in the root directory of the present distribution,
!  or http://www.gnu.org/copyleft/gpl.txt .
!

#if defined(__AIX) || defined(__MAC64) || defined(FUJ64)|| defined(__ALPHA) || defined(__SX6) || defined(__LINUX64) || defined(__HP64) || defined(__ALTIX)|| defined(__ORIGIN) 
#  define C_POINTER  integer(kind=8)
#else
#  define C_POINTER  integer(kind=4)
#endif

#if defined(__ALPHA)
#  define DIRECT_IO_FACTOR 2
#else
#  define DIRECT_IO_FACTOR 8 
#endif

#if defined (__XLF) || defined(__ABSOFT)
#  define flush  flush_
#endif

#if defined(__ABSOFT)
#  define getenv getenv_
#  define getarg getarg_
#  define iargc  iargc_
#endif

#if defined(__LAM) && ( defined (__LINUX) || defined (__LINUX64) )
#  define MPI_REAL8 MPI_DOUBLE_PRECISION
#endif

#  define DREAL       @@_use_DBLE_instead@@
#  define DIMAG       @@_use_AIMAG_instead@@
#  define DCMPLX      @@_use_CMPLX_instead@@
#  define CMPLX(a,b)  cmplx(a,b,kind=DP)

#if defined(ADD_BLAS_TWO_UNDERSCORES)

#    define DAXPY    daxpy__
#    define DCOPY    dcopy__
#    define DDOT     ddot__
#    define DGETRF   dgetrf__
#    define DGETRI   dgetri__
#    define DGEMV    dgemv__
#    define DGEMM    dgemm__
#    define DGER     dger__
#    define DNRM2    dnrm2__
#    define DPOTRF   dpotrf__
#    define DPOTRS   dpotrs__
#    define DSCAL    dscal__
#    define DSPEV    dspev__
#    define DSYTRF   dsytrf__
#    define DSYTRI   dsytri__
#    define DSYEV    dsyev__
#    define DSYGV    dsygv__
#    define DSYGVX   dsygvx__
#    define DSWAP    dswap__
#    define ILAENV   ilaenv__
#    define IDAMAX   idamax__
#    define IZAMAX   izamax__
#    define ZAXPY    zaxpy__
#    define ZCOPY    zcopy__
#    define ZDOTC    zdotc__
#    define ZDOTU    zdotu__
#    define ZGEMM    zgemm__
#    define ZGEMV    zgemv__
#    define ZGESV    zgesv__
#    define ZGESVD   zgesvd__
#    define ZGGEV    zggev__
#    define ZHEEV    zheev__
#    define ZHEEVX   zheevx__
#    define ZHEGV    zhegv__
#    define ZHEGVX   zhegvx__
#    define ZHPEV    zhpev__
#    define ZSCAL    zscal__
#    define ZSWAP    zswap__
#    define ZDSCAL   zdscal__

#  elif defined(ADD_BLAS_ONE_UNDERSCORE)

#    define DAXPY    daxpy_
#    define DCOPY    dcopy_
#    define DDOT     ddot_
#    define DGETRF   dgetrf_
#    define DGETRI   dgetri_
#    define DGEMV    dgemv_
#    define DGEMM    dgemm_
#    define DGER     dger_
#    define DNRM2    dnrm2_
#    define DPOTRF   dpotrf_
#    define DPOTRS   dpotrs_
#    define DSCAL    dscal_
#    define DSPEV    dspev_
#    define DSYTRF   dsytrf_
#    define DSYTRI   dsytri_
#    define DSYEV    dsyev_
#    define DSYGV    dsygv_
#    define DSYGVX   dsygvx_
#    define DSWAP    dswap_
#    define ILAENV   ilaenv_
#    define IDAMAX   idamax_
#    define IZAMAX   izamax_
#    define ZAXPY    zaxpy_
#    define ZCOPY    zcopy_
#    define ZDOTC    zdotc_
#    define ZDOTU    zdotu_
#    define ZGEMM    zgemm_
#    define ZGEMV    zgemv_
#    define ZGESV    zgesv_
#    define ZGESVD   zgesvd_
#    define ZGGEV    zggev_
#    define ZHEEV    zheev_
#    define ZHEEVX   zheevx_
#    define ZHEGV    zhegv_
#    define ZHEGVX   zhegvx_
#    define ZHPEV    zhpev_
#    define ZSCAL    zscal_
#    define ZSWAP    zswap_
#    define ZDSCAL   zdscal_

#  else

#    define DAXPY    daxpy
#    define DCOPY    dcopy
#    define DDOT     ddot
#    define DGETRF   dgetrf
#    define DGETRI   dgetri
#    define DGEMV    dgemv
#    define DGEMM    dgemm
#    define DGER     dger
#    define DNRM2    dnrm2
#    define DPOTRF   dpotrf
#    define DPOTRS   dpotrs
#    define DSCAL    dscal
#    define DSPEV    dspev
#    define DSYTRF   dsytrf
#    define DSYTRI   dsytri
#    define DSYEV    dsyev
#    define DSYGV    dsygv
#    define DSYGVX   dsygvx
#    define DSWAP    dswap
#    define ILAENV   ilaenv
#    define IDAMAX   idamax
#    define IZAMAX   izamax
#    define ZAXPY    zaxpy
#    define ZCOPY    zcopy
#    define ZDOTC    zdotc
#    define ZDOTU    zdotu
#    define ZGEMM    zgemm
#    define ZGEMV    zgemv
#    define ZGESV    zgesv
#    define ZGESVD   zgesvd
#    define ZGGEV    zggev
#    define ZHEEV    zheev
#    define ZHEEVX   zheevx
#    define ZHEGV    zhegv
#    define ZHEGVX   zhegvx
#    define ZHPEV    zhpev
#    define ZSCAL    zscal
#    define ZSWAP    zswap
#    define ZDSCAL   zdscal

#endif

