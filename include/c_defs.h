/*
  Copyright (C) 2002-2006 Quantum-ESPRESSO group
  This file is distributed under the terms of the
  GNU General Public License. See the file `License'
  in the root directory of the present distribution,
  or http://www.gnu.org/copyleft/gpl.txt .
*/

/* Machine-dependent Fortran to C calling convention
   redefine C symbols so that Fortran finds them
   Absoft:
      capital letters, no added underscores (leave as is)
   XLF (Aix, Mac OS-X), HP-UX:
      lowercase, with no added underscores 
   G95, EKOPath, Alpha Linux:
      lowercase, with one added underscore if the name does
      not contain underscores, with two if it does!
   Most other cases: 
      lowercase, with one added underscore 
*/

#if !defined (__ABSOFT)

/* Absoft: do nothing */

#  if defined (__XLF) || defined (__HP)

/* convert to lowercase */

#  define FFTW_INPLACE_DRV_1D fftw_inplace_drv_1d
#  define FFTW_INPLACE_DRV_2D fftw_inplace_drv_2d
#  define FFTW_INPLACE_DRV_3D fftw_inplace_drv_3d
#  define CREATE_PLAN_2D create_plan_2d
#  define CREATE_PLAN_3D create_plan_3d
#  define CREATE_PLAN_1D create_plan_1d
#  define DESTROY_PLAN_1D destroy_plan_1d
#  define DESTROY_PLAN_2D destroy_plan_2d
#  define DESTROY_PLAN_3D destroy_plan_3d
#  define FFT_X_STICK fft_x_stick
#  define FFT_Y_STICK fft_y_stick
#  define FFT_Z_STICK fft_z_stick
#  define LN_ALLOC ln_alloc 
#  define LN_DEALLOC ln_dealloc
#  define LN_SET ln_set
#  define LN_ACTIVATE ln_activate
#  define LN_IND ln_ind
#  define CCLOCK cclock
#  define ELAPSED_SECONDS elapsed_seconds
#  define C_MKDIR c_mkdir
#  define MEMSTAT memstat

#  elif defined (__G95) || defined (__EKO) || (defined __ALPHA && defined __LINUX64)

/* convert to lowercase, add one or two underscores */

#  define FFTW_INPLACE_DRV_1D fftw_inplace_drv_1d__
#  define FFTW_INPLACE_DRV_2D fftw_inplace_drv_2d__
#  define FFTW_INPLACE_DRV_3D fftw_inplace_drv_3d__
#  define CREATE_PLAN_1D create_plan_1d__
#  define CREATE_PLAN_2D create_plan_2d__
#  define CREATE_PLAN_3D create_plan_3d__
#  define DESTROY_PLAN_1D destroy_plan_1d__
#  define DESTROY_PLAN_2D destroy_plan_2d__
#  define DESTROY_PLAN_3D destroy_plan_3d__
#  define FFT_X_STICK fft_x_stick__
#  define FFT_Y_STICK fft_y_stick__
#  define FFT_Z_STICK fft_z_stick__
#  define LN_ALLOC ln_alloc__
#  define LN_DEALLOC ln_dealloc__
#  define LN_SET ln_set__
#  define LN_ACTIVATE ln_activate__
#  define LN_IND ln_ind__
#  define CCLOCK cclock_
#  define ELAPSED_SECONDS elapsed_seconds__
#  define C_MKDIR c_mkdir__
#  define MEMSTAT memstat_

#  else

/* convert to lowercase, add one underscore */

#  define FFTW_INPLACE_DRV_1D fftw_inplace_drv_1d_
#  define FFTW_INPLACE_DRV_2D fftw_inplace_drv_2d_
#  define FFTW_INPLACE_DRV_3D fftw_inplace_drv_3d_
#  define CREATE_PLAN_1D create_plan_1d_
#  define CREATE_PLAN_2D create_plan_2d_
#  define CREATE_PLAN_3D create_plan_3d_
#  define DESTROY_PLAN_1D destroy_plan_1d_
#  define DESTROY_PLAN_2D destroy_plan_2d_
#  define DESTROY_PLAN_3D destroy_plan_3d_
#  define FFT_X_STICK fft_x_stick_
#  define FFT_Y_STICK fft_y_stick_
#  define FFT_Z_STICK fft_z_stick_
#  define LN_ALLOC ln_alloc_ 
#  define LN_DEALLOC ln_dealloc_
#  define LN_SET ln_set_
#  define LN_ACTIVATE ln_activate_
#  define LN_IND ln_ind_
#  define CCLOCK cclock_
#  define ELAPSED_SECONDS elapsed_seconds_
#  define C_MKDIR c_mkdir_
#  define MEMSTAT memstat_

#  endif

#endif
