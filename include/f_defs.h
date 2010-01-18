!
!  Copyright (C) 2002-2006 Quantum ESPRESSO group
!  This file is distributed under the terms of the
!  GNU General Public License. See the file `License'
!  in the root directory of the present distribution,
!  or http://www.gnu.org/copyleft/gpl.txt .
!
#  define DREAL       @@_use_DBLE_instead@@
#  define DIMAG       @@_use_AIMAG_instead@@
#  define DCMPLX      @@_use_CMPLX_instead@@
#  define DFLOAT      @@_use_DBLE_instead@@
#  define CMPLX(a,b)  cmplx(a,b,kind=DP)

#if defined(ADD_BLAS_TWO_UNDERSCORES)
#    define C_NAME(name) name ## __
#elif defined(ADD_BLAS_ONE_UNDERSCORE)
#    define C_NAME(name) name ## _
#else
#    define C_NAME(name) name
#endif

#define DAXPY    C_NAME(daxpy)
#define DCOPY    C_NAME(dcopy)
#define DDOT     C_NAME(ddot)
#define DGETRF   C_NAME(dgetrf)
#define DGETRI   C_NAME(dgetri)
#define DGEMV    C_NAME(dgemv)
#define DGEMM    C_NAME(dgemm)
#define DGER     C_NAME(dger)
#define DNRM2    C_NAME(dnrm2)
#define DPOTRF   C_NAME(dpotrf)
#define DPOTRS   C_NAME(dpotrs)
#define DSCAL    C_NAME(dscal)
#define DSPEV    C_NAME(dspev)
#define DSYTRF   C_NAME(dsytrf)
#define DSYTRI   C_NAME(dsytri)
#define DSYEV    C_NAME(dsyev)
#define DSYGV    C_NAME(dsygv)
#define DSYGVX   C_NAME(dsygvx)
#define DSWAP    C_NAME(dswap)
#define ILAENV   C_NAME(ilaenv)
#define IDAMAX   C_NAME(idamax)
#define IZAMAX   C_NAME(izamax)
#define ZAXPY    C_NAME(zaxpy)
#define ZCOPY    C_NAME(zcopy)
#define ZDOTC    C_NAME(zdotc)
#define ZDOTU    C_NAME(zdotu)
#define ZGEMM    C_NAME(zgemm)
#define ZGEMV    C_NAME(zgemv)
#define ZGESV    C_NAME(zgesv)
#define ZGESVD   C_NAME(zgesvd)
#define ZGGEV    C_NAME(zggev)
#define ZHEEV    C_NAME(zheev)
#define ZHEEVX   C_NAME(zheevx)
#define ZHEGV    C_NAME(zhegv)
#define ZHEGVX   C_NAME(zhegvx)
#define ZHPEV    C_NAME(zhpev)
#define ZSCAL    C_NAME(zscal)
#define ZSWAP    C_NAME(zswap)
#define ZDSCAL   C_NAME(zdscal)

