MODULE hip_kernels
  use omp_lib
  USE, INTRINSIC :: iso_c_binding

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: scalar_init, scalar_multiply, scalar_multiply_3D, loop2d_scatter_hip

#if defined(__HIP)
  INTERFACE
    SUBROUTINE scalar_init_(s, dev_ptr,val, stream) &
        & BIND(C, name="c_scalar_init_")
      USE iso_c_binding
      INTEGER(C_INT), INTENT(IN), VALUE   :: s
      TYPE(C_PTR), VALUE                  :: dev_ptr
      REAL(C_DOUBLE), VALUE               :: val
      TYPE(C_PTR),    VALUE               :: stream
    END SUBROUTINE scalar_init_
  END INTERFACE

  INTERFACE
    SUBROUTINE scalar_multiply_(s, dev_ptr,val, stream) &
        & BIND(C, name="c_scalar_multiply_")
      USE iso_c_binding
      INTEGER(C_INT), INTENT(IN), VALUE   :: s
      TYPE(C_PTR), VALUE                  :: dev_ptr
      REAL(C_DOUBLE), VALUE               :: val
     TYPE(C_PTR),    VALUE               :: stream
    END SUBROUTINE scalar_multiply_
  END INTERFACE

  INTERFACE
    SUBROUTINE loop2d_scatter_hip_( drz, f_in, f_out, dft_ismap, nppx, nnp, of1, of2, npp, nswip, stream ) &
        & BIND(C, name="loop2d_scatter_hip_")
      USE iso_c_binding
      TYPE(C_PTR), VALUE  :: f_in, dft_ismap
      TYPE(C_PTR), VALUE :: f_out
      INTEGER(C_INT), INTENT(in), VALUE :: drz, nppx, nnp, npp, nswip, of1, of2
      TYPE(C_PTR), VALUE :: stream
    END SUBROUTINE loop2d_scatter_hip_
  END INTERFACE
#endif

  CONTAINS

  SUBROUTINE scalar_init(a,val,s,stream)
    COMPLEX(8), INTENT(inout)  :: a(:)
    REAL(8), INTENT(in)        :: val
    INTEGER(C_INT), INTENT(in) :: s
    TYPE(C_PTR)                :: stream
#if defined(__HIP)
    !$omp target data use_device_addr(a)
    CALL scalar_init_(s,c_loc(a),val,stream)
    !$omp end target data
#endif
  END SUBROUTINE scalar_init

  SUBROUTINE scalar_multiply(a,val,s,stream)
    COMPLEX(8), INTENT(inout)  :: a(:)
    REAL(8), INTENT(in)        :: val
    INTEGER(C_INT), INTENT(in) :: s
    TYPE(C_PTR)                :: stream
#if defined(__HIP)
    !$omp target data use_device_addr(a)
    CALL scalar_multiply_(s,c_loc(a),val,stream)
    !$omp end target data
#endif
  END SUBROUTINE scalar_multiply

  SUBROUTINE scalar_multiply_3D(a,val,s,stream)
    COMPLEX(8), INTENT(inout)  :: a(:,:,:)
    REAL(8), INTENT(in)        :: val
    INTEGER(C_INT), INTENT(in) :: s
    TYPE(C_PTR)                :: stream
#if defined(__HIP)
    !$omp target data use_device_addr(a)
    CALL scalar_multiply_(s,c_loc(a),val,stream)
    !$omp end target data
#endif
  END SUBROUTINE scalar_multiply_3D

  SUBROUTINE loop2d_scatter_hip( drz, f_in, f_out, dft_ismap, nppx, nnp, of1, of2, npp, nswip, stream )
    COMPLEX(8), INTENT(in)  :: f_in(:)
    COMPLEX(8), INTENT(inout) :: f_out(:)
    INTEGER, INTENT(in) :: dft_ismap(:)
    INTEGER(C_INT), INTENT(in) :: drz, nppx, nnp, npp, nswip, of1, of2
    TYPE(C_PTR) :: stream
#if defined(__HIP)
    !$omp target data use_device_addr(f_in, f_out, dft_ismap)
    CALL loop2d_scatter_hip_( drz, c_loc(f_in), c_loc(f_out), c_loc(dft_ismap), nppx, nnp, of1, of2, &
                              npp, nswip, stream )
    !$omp end target data
#endif
  END SUBROUTINE loop2d_scatter_hip
  
END MODULE hip_kernels






