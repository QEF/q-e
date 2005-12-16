
! these four routines fix a bug of the Accelerate.framework implementation 
! of BLAS on Mac OS X . Copied from:
! http://developer.apple.com/hardware/ve/errata.html#fortran_conventions
! by Stefano Baroni, December 10, 2005

      double complex function zdotc(n, zx, incx, zy, incy)
      double complex zx(*), zy(*), z
      integer n, incx, incy
      
      call cblas_zdotc_sub(%val(n), zx, %val(incx), zy, %val(incy), z)
      
      zdotc = z
      return
      end
      
      double complex function zdotu(n, zx, incx, zy, incy)
      double complex zx(*), zy(*), z
      integer n, incx, incy
      
      call cblas_zdotu_sub(%val(n), zx, %val(incx), zy, %val(incy), z)
      
      zdotu = z
      return
      end
      
      complex function cdotc(n, cx, incx, cy, incy)
      complex cx(*), cy(*), c
      integer n, incx, incy
      
      call cblas_cdotc_sub(%val(n), cx, %val(incx), cy, %val(incy), c)
      
      cdotc = c
      return
      end

      complex function cdotu(n, cx, incx, cy, incy)
      complex cx(*), cy(*), c
      integer n, incx, incy
      
      call cblas_cdotu_sub(%val(n), cx, %val(incx), cy, %val(incy), c)
      
      cdotu = c
      return
      end
