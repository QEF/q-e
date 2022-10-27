MODULE fftw_interfaces
  
  INTERFACE
     SUBROUTINE create_plan_1d( plan, nz, i ) BIND(C,name="create_plan_1d")
       USE iso_c_binding
       TYPE(C_PTR) :: plan
       INTEGER(KIND=C_INT) :: nz, i
     END SUBROUTINE create_plan_1d
     SUBROUTINE destroy_plan_1d( plan ) BIND(C,name="destroy_plan_1d") 
       USE iso_c_binding
       TYPE(C_PTR) :: plan
     END SUBROUTINE destroy_plan_1d
     
     SUBROUTINE create_plan_2d( plan, nx, ny, i ) BIND(C,name="create_plan_2d")
       USE iso_c_binding
       TYPE(C_PTR) :: plan
       INTEGER(KIND=C_INT) :: nx, ny, i
     END SUBROUTINE create_plan_2d
     SUBROUTINE destroy_plan_2d( plan ) BIND(C,name="destroy_plan_2d")
       USE iso_c_binding
       TYPE(C_PTR) :: plan
     END SUBROUTINE destroy_plan_2d
     
     SUBROUTINE create_plan_3d( plan, nx, ny, nz, i ) BIND(C,name="create_plan_3d")
       USE iso_c_binding
       TYPE(C_PTR) :: plan
       INTEGER(KIND=C_INT) :: nx, ny, nz, i
     END SUBROUTINE create_plan_3d
     SUBROUTINE destroy_plan_3d( plan ) BIND(C,name="destroy_plan_3d")
       USE iso_c_binding
       TYPE(C_PTR) :: plan
     END SUBROUTINE destroy_plan_3d
     
     SUBROUTINE fft_x_stick( plan, xy, nx, ny, nz, ldx, ldy ) &
          BIND(C,name="fft_x_stick")
       USE iso_c_binding
       TYPE(C_PTR) :: plan
       INTEGER(KIND=C_INT) :: nx, ny, nz, ldx, ldy
       COMPLEX(KIND=C_DOUBLE) :: xy
     END SUBROUTINE fft_x_stick

     SUBROUTINE fft_x_stick_single( plan, xy, nx, ny, nz, ldx, ldy ) &
          BIND(C,name="fft_x_stick_single")
       USE iso_c_binding
       TYPE(C_PTR) :: plan
       INTEGER(KIND=C_INT) :: nx, ny, nz, ldx, ldy
       COMPLEX(KIND=C_DOUBLE) :: xy
     END SUBROUTINE fft_x_stick_single
     
     SUBROUTINE fft_y_stick( plan, xy, ny, ldx ) BIND(C,name="fft_y_stick")
       USE iso_c_binding
       TYPE(C_PTR) :: plan
       INTEGER(KIND=C_INT) :: ny, ldx
       COMPLEX(KIND=C_DOUBLE) :: xy
     END SUBROUTINE fft_y_stick
     
     SUBROUTINE fft_z_stick( plan, z, ldz, nzl ) BIND(C,name="fft_z_stick")
       USE iso_c_binding
       TYPE(C_PTR) :: plan
       INTEGER(KIND=C_INT) :: ldz, nzl
       COMPLEX(KIND=C_DOUBLE) :: z
     END SUBROUTINE fft_z_stick
     
     SUBROUTINE fft_z_stick_single( plan, z, ldz ) BIND(C,name="fft_z_stick_single")
       USE iso_c_binding
       TYPE(C_PTR) :: plan
       INTEGER(KIND=C_INT) :: ldz
       COMPLEX(KIND=C_DOUBLE) :: z
     END SUBROUTINE fft_z_stick_single
     
     SUBROUTINE fftw_inplace_drv_1d( plan, m, z, inc1, inc2 ) &
          BIND(C,name="fftw_inplace_drv_1d")
       USE iso_c_binding
       TYPE(C_PTR) :: plan
       INTEGER(KIND=C_INT) :: m, inc1, inc2
       COMPLEX(KIND=C_DOUBLE) :: z
     END SUBROUTINE fftw_inplace_drv_1d
     SUBROUTINE fftw_inplace_drv_3d( plan, m, z, inc1, inc2 ) &
          BIND(C,name="fftw_inplace_drv_3d")
       USE iso_c_binding
       TYPE(C_PTR) :: plan
       INTEGER(KIND=C_INT) :: m, inc1, inc2
       COMPLEX(KIND=C_DOUBLE) :: z
     END SUBROUTINE fftw_inplace_drv_3d
  END INTERFACE

END MODULE fftw_interfaces
