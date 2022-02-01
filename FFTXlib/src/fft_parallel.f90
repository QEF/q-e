!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------==!
!
!     Parallel 3D FFT high level Driver
!     ( Charge density and Wave Functions )
!
!     Written and maintained by Carlo Cavazzoni !     Last update Apr. 2009
!     Rewritten  by Stefano de Gironcoli        !     Sep-Nov 2016
!
!!=---------------------------------------------------------------------==!
!
MODULE fft_parallel
!
   USE fft_param
   IMPLICIT NONE
   SAVE
!
CONTAINS
!
!  General purpose driver
!
!----------------------------------------------------------------------------
SUBROUTINE tg_cft3s( f, dfft, isgn )
  !----------------------------------------------------------------------------
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +-3 : parallel 3d fft for wavefunctions with task group
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (cft_1z)
  !! ...              transpose across nodes           (fft_scatter_yz)
  !! ...              fft along y using pencils        (cft_1y)
  !! ...              transpose across nodes           (fft_scatter_xy)
  !! ...              fft along x using pencils        (cft_1x)
  !
  !! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  !! ...              fft along x using pencils        (cft_1x)
  !! ...              transpose across nodes           (fft_scatter_xy)
  !! ...              fft along y using pencils        (cft_1y)
  !! ...              transpose across nodes           (fft_scatter_yz)
  !! ...              fft along z using pencils        (cft_1z)
  !
  ! If task_group_fft_is_active the FFT acts on a number of wfcs equal to
  ! dfft%nproc2, the number of Y-sections in which a plane is divided.
  ! Data are reshuffled by the fft_scatter_tg routine so that each of the
  ! dfft%nproc2 subgroups (made by dfft%nproc3 procs) deals with whole planes
  ! of a single wavefunciton.
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  !
  USE fft_scalar, ONLY : cft_1z
  USE fft_scatter,ONLY : fft_scatter_xy, fft_scatter_yz, &
          fft_scatter_tg, fft_scatter_tg_opt
  USE fft_types,  ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfft   ! descriptor of fft data layout
  COMPLEX(DP), INTENT(inout)             :: f( : ) ! array containing data to be transformed
  INTEGER, INTENT(in)                    :: isgn   ! fft direction (potential: +/-1, wave: +/-2, wave_tg: +/-3)
  !
  INTEGER                  :: n1, n2, n3, nx1, nx2, nx3
  INTEGER                  :: nnr_
  INTEGER                  :: nsticks_x, nsticks_y, nsticks_z
  COMPLEX(DP), ALLOCATABLE :: aux (:)
  INTEGER                  :: i
  !
  !write (6,*) 'enter tg_cft3s ',isgn ; write(6,*) ; FLUSH(6)
  n1  = dfft%nr1  ; n2  = dfft%nr2  ; n3  = dfft%nr3
  nx1 = dfft%nr1x ; nx2 = dfft%nr2x ; nx3 = dfft%nr3x
  !
  if (abs(isgn) == 1 ) then       ! potential fft
     nnr_ = dfft%nnr
     nsticks_x = dfft%my_nr2p * dfft%my_nr3p
     nsticks_y = dfft%nr1p(dfft%mype2+1) * dfft%my_nr3p
     nsticks_z = dfft%nsp(dfft%mype+1)
  else if (abs(isgn) == 2 ) then  ! wave func fft
     nnr_ = dfft%nnr
     nsticks_x = dfft%my_nr2p * dfft%my_nr3p
     nsticks_y = dfft%nr1w(dfft%mype2+1) * dfft%my_nr3p
     nsticks_z = dfft%nsw(dfft%mype+1)
  else if (abs(isgn) == 3 ) then  ! wave func fft with task groups
     nnr_ = dfft%nnr_tg
     nsticks_x = dfft%nr2 * dfft%my_nr3p
     nsticks_y = dfft%nr1w_tg * dfft%my_nr3p
     nsticks_z = dfft%nsw_tg(dfft%mype+1)
  else
     CALL fftx_error__( ' tg_cft3s', ' wrong value of isgn ', 10+abs(isgn) )
  end if
  ALLOCATE( aux( nnr_ ) )
  !
  IF ( isgn > 0 ) THEN  ! G -> R
     if (isgn==+3) then
        call fft_scatter_tg_opt ( dfft, f, aux, nnr_, isgn)
     else
!$omp parallel do
        do i=1, nsticks_z*nx3
           aux(i) = f(i)
        enddo
!$omp end parallel do
     endif
     ! Gz, Gy, Gx
     CALL cft_1z( aux, nsticks_z, n3, nx3, isgn, f )
     ! Rz, Gy, Gx
     CALL fft_scatter_yz ( dfft, f, aux, nnr_, isgn )
     ! Gy, Gx, Rz
     CALL cft_1z( aux, nsticks_y, n2, nx2, isgn, f )
     ! Ry, Gx, Rz
     CALL fft_scatter_xy ( dfft, f, aux, nnr_, isgn )
     ! Gx, Ry, Rz
     CALL cft_1z( aux, nsticks_x, n1, nx1, isgn, f )
     ! Rx, Ry, Rz
     ! clean garbage beyond the intended dimension.
     if (nsticks_x*nx1 < nnr_) f(nsticks_x*nx1+1:nnr_) = (0.0_DP,0.0_DP)
     !
  ELSE                  ! R -> G
     !
     ! Rx, Ry, Rz
     CALL cft_1z( f, nsticks_x, n1, nx1, isgn, aux )
     ! Gx, Ry, Rz
     CALL fft_scatter_xy ( dfft, f, aux, nnr_, isgn )
     ! Ry, Gx, Rz
     CALL cft_1z( f, nsticks_y, n2, nx2, isgn, aux )
     ! Gy, Gx, Rz
     CALL fft_scatter_yz ( dfft, f, aux, nnr_, isgn )
     ! Rz, Gy, Gx
     CALL cft_1z( f, nsticks_z, n3, nx3, isgn, aux )
     ! Gz, Gy, Gx
     if (isgn==-3) then
        call fft_scatter_tg_opt ( dfft, aux, f, nnr_, isgn)
     else
!$omp parallel do
        do i=1, nsticks_z*nx3
           f(i) = aux(i)
        enddo
!$omp end parallel do
     endif
  ENDIF
  !write (6,99) f(1:400); write(6,*); FLUSH(6)
  !
  DEALLOCATE( aux )
  !
  !if (.true.) stop
  RETURN
99 format ( 20 ('(',2f12.9,')') )

  !
END SUBROUTINE tg_cft3s

!
!  Specific driver for the new 'many' call
!
!----------------------------------------------------------------------------
SUBROUTINE many_cft3s( f, dfft, isgn, howmany )
  !----------------------------------------------------------------------------
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +-3 : parallel 3d fft for wavefunctions with task group
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (cft_1z)
  !! ...              transpose across nodes           (fft_scatter_yz)
  !! ...              fft along y using pencils        (cft_1y)
  !! ...              transpose across nodes           (fft_scatter_xy)
  !! ...              fft along x using pencils        (cft_1x)
  !
  !! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  !! ...              fft along x using pencils        (cft_1x)
  !! ...              transpose across nodes           (fft_scatter_xy)
  !! ...              fft along y using pencils        (cft_1y)
  !! ...              transpose across nodes           (fft_scatter_yz)
  !! ...              fft along z using pencils        (cft_1z)
  !
  ! If task_group_fft_is_active the FFT acts on a number of wfcs equal to
  ! dfft%nproc2, the number of Y-sections in which a plane is divided.
  ! Data are reshuffled by the fft_scatter_tg routine so that each of the
  ! dfft%nproc2 subgroups (made by dfft%nproc3 procs) deals with whole planes
  ! of a single wavefunciton.
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  !
  USE fft_scalar,  ONLY : cft_1z
  USE fft_scatter, ONLY : fft_scatter_xy, fft_scatter_yz, &
         fft_scatter_tg_opt, fft_scatter_many_xy, fft_scatter_many_yz
  USE fft_types,   ONLY : fft_type_descriptor
  USE omp_lib
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(inout), TARGET :: dfft     ! descriptor of fft data layout
  COMPLEX(DP), INTENT(inout)                        :: f( : )   ! array containing data to be transformed
  INTEGER, INTENT(in)                               :: isgn     ! fft direction (potential: +/-1, wave: +/-2, wave_tg: +/-3)
  INTEGER, INTENT(in)                               :: howmany  ! number of FFTs grouped together
  !
  INTEGER                          :: n1, n2, n3, nx1, nx2, nx3
  INTEGER                          :: nnr_
  INTEGER                          :: nsticks_x, nsticks_y, nsticks_z
  INTEGER                          :: nsticks_yx, nsticks_zx
  INTEGER                          :: i, j
  COMPLEX(DP), POINTER             :: aux(:)
  !
  !write (6,*) 'enter tg_cft3s ',isgn ; write(6,*) ; FLUSH(6)
  n1  = dfft%nr1  ; n2  = dfft%nr2  ; n3  = dfft%nr3
  nx1 = dfft%nr1x ; nx2 = dfft%nr2x ; nx3 = dfft%nr3x
  !
  nnr_ = dfft%nnr
  nsticks_x  = dfft%my_nr2p  * dfft%my_nr3p
  if (abs(isgn) == 1 ) then       ! potential fft
     nsticks_y  = dfft%nr1p(dfft%mype2+1) * dfft%my_nr3p
     nsticks_yx = MAXVAL(dfft%nr1p) * MAXVAL(dfft%nr3p)
     nsticks_z  = dfft%nsp(dfft%mype+1)
     nsticks_zx = MAXVAL(dfft%nsp)
  else if (abs(isgn) == 2 ) then  ! wave func fft
     nsticks_y  = dfft%nr1w(dfft%mype2+1) * dfft%my_nr3p
     nsticks_yx = MAXVAL(dfft%nr1w) * MAXVAL(dfft%nr3p)
     nsticks_z  = dfft%nsw(dfft%mype+1)
     nsticks_zx = MAXVAL(dfft%nsw)
  else if (abs(isgn) == 3 ) then  ! wave func fft with task groups
     CALL fftx_error__( ' many_cft3s', ' Taskgroup and many not supported ', 10+abs(isgn) )
  else
     CALL fftx_error__( ' many_cft3s', ' wrong value of isgn ', 10+abs(isgn) )
  end if
  aux => dfft%aux
  !
#if defined(__FFT_OPENMP_TASKS)
  call execute_using_tasks()
#else
  call execute_using_threads()
#endif
  !write (6,99) f_d(1:400); write(6,*); FLUSH(6)
  !
  RETURN
99 format ( 20 ('(',2f12.9,')') )
  !
CONTAINS
  !
#if defined(__FFT_OPENMP_TASKS)
  !
  SUBROUTINE execute_using_tasks()
    !
    ! This implementation requires thread-safe cft_1z and multi-threaded MPI
    !
    IMPLICIT NONE
    !
    INTEGER :: comm2s
    !
#if defined(_OPENMP) && !defined(__FFT_SCALAR_THREAD_SAFE)
    CALL fftx_error__('many_cft3s::execute_using_tasks', &
    'needs thread-safe fft_scalar backend selected at compile time.',1)
#endif
    !
    IF ( isgn > 0 ) THEN  ! G -> R
!$omp parallel default(none)                        &
!$omp          private(i, j, comm2s)                &
!$omp          shared(howmany, f, nnr_, isgn, dfft) &
!$omp          shared(nsticks_z, n3, nx3)           &
!$omp          shared(nsticks_y, n2, nx2)           &
!$omp          shared(nsticks_x, n1, nx1)           &
!$omp          shared(nsticks_zx, nsticks_yx, aux)
       !
!$omp do
       DO i = 0, howmany-1
          CALL cft_1z( f(i*nnr_+1:), nsticks_z, n3, nx3, isgn, aux(nx3*nsticks_zx*i+1:) )
       ENDDO
!$omp end do
       !
!$omp single
       CALL fft_scatter_many_yz( dfft, aux, f, isgn, howmany )
       !
       DO i = 0, howmany-1
!$omp task depend(out:f(i*nnr_+1:(i+1)*nnr_))
          CALL cft_1z( f(i*nnr_+1:), nsticks_y, n2, nx2, isgn, aux(i*nnr_+1:), in_place=.true. )
!$omp end task
!$omp task depend(in:f(i*nnr_+1:(i+1)*nnr_)) depend(out:aux(i*nnr_+1:(i+1)*nnr_))
          comm2s = dfft%comm2s( i + 1 )
          CALL fft_scatter_xy( dfft, f(i*nnr_+1:), aux(i*nnr_+1:), nnr_, isgn, comm2s )
!$omp end task
!$omp task depend(in:aux(i*nnr_+1:(i+1)*nnr_))
          CALL cft_1z( aux(i*nnr_+1:), nsticks_x, n1, nx1, isgn, f(i*nnr_+1:) )
!$omp end task
       ENDDO
!$omp end single
       !
       if (nsticks_x*nx1 < nnr_) then
!$omp do collapse(2)
          DO i = 0, howmany-1
             do j=nsticks_x*nx1+1, nnr_
                f(j+i*nnr_) = (0.0_DP,0.0_DP)
             end do
          END DO
!$omp end do
       endif
!$omp end parallel
    ELSE                  ! R -> G
!$omp parallel default(none)                        &
!$omp          private(i, j, comm2s)                &
!$omp          shared(howmany, f, isgn, nnr_, dfft) &
!$omp          shared(nsticks_z, n3, nx3)           &
!$omp          shared(nsticks_y, n2, nx2)           &
!$omp          shared(nsticks_x, n1, nx1)           &
!$omp          shared(nsticks_zx, nsticks_yx, aux)
       !
!$omp single
!$omp taskgroup
       DO i = 0, howmany-1
!$omp task depend(out:aux(i*nnr_+1:(i+1)*nnr_))
          CALL cft_1z( f(i*nnr_+1:), nsticks_x, n1, nx1, isgn, aux(i*nnr_+1:) )
!$omp end task
!$omp task depend(in:aux(i*nnr_+1:(i+1)*nnr_)) depend(out:f(i*nnr_+1:(i+1)*nnr_))
          comm2s = dfft%comm2s( i + 1 )
          CALL fft_scatter_xy( dfft, f(i*nnr_+1:), aux(i*nnr_+1:), nnr_, isgn, comm2s )
!$omp end task
!$omp task depend(in:f(i*nnr_+1:(i+1)*nnr_))
          CALL cft_1z( f(i*nnr_+1:), nsticks_y, n2, nx2, isgn, aux(i*nnr_+1:), in_place=.true. )
!$omp end task
       ENDDO
!$omp end taskgroup
       !
       CALL fft_scatter_many_yz( dfft, aux, f, isgn, howmany )
!$omp end single
       !
!$omp do
       DO i = 0, howmany-1
          CALL cft_1z( aux(nx3*nsticks_zx*i+1:), nsticks_z, n3, nx3, isgn, f(i*nnr_+1:) )
       ENDDO
!$omp end do
       !
!$omp do collapse(2)
       DO i = 0, howmany-1
          DO j=nsticks_z*nx3+1, nnr_
             f(j+i*nnr_) = (0.0_DP,0.0_DP)
          ENDDO
       ENDDO
!$omp end do
!$omp end parallel
    ENDIF
  ENDSUBROUTINE execute_using_tasks
#endif
  !
  SUBROUTINE execute_using_threads()
    !
    ! This implementation requires thread-safe cft_1z
    !
    IMPLICIT NONE
    !
#if defined(_OPENMP) && !defined(__FFT_SCALAR_THREAD_SAFE)
    CALL fftx_error__('many_cft3s::execute_using_threads', &
    'needs thread-safe fft_scalar backend selected at compile time.',1)
#endif
    !
    IF ( isgn > 0 ) THEN  ! G -> R
!$omp parallel default(none)                                        &
!$omp          private(i, j)                                        &
!$omp          shared(howmany, f, nnr_, isgn, dfft, aux)            &
!$omp          shared(nsticks_z, n3, nx3)                           &
!$omp          shared(nsticks_y, n2, nx2)                           &
!$omp          shared(nsticks_x, n1, nx1)                           &
!$omp          shared(nsticks_zx, nsticks_yx)
    !
!$omp do
       DO i = 0, howmany-1
          DO j=1, nsticks_z*nx3
             aux(j+i*nnr_) = f(j+i*nnr_)
          ENDDO
       ENDDO
!$omp end do
       !
!$omp do
       DO i = 0, howmany-1
          CALL cft_1z( aux(i*nnr_+1:), nsticks_z, n3, nx3, isgn, f(nx3*nsticks_zx*i+1:) )
       ENDDO
!$omp end do
       !
!$omp single
       CALL fft_scatter_many_yz( dfft, f, dfft%aux, isgn, howmany )
!$omp end single
       !
!$omp do
       DO i = 0, howmany-1
          CALL cft_1z( aux(i*nnr_+1:), nsticks_y, n2, nx2, isgn, f(nx2*nsticks_yx*i+1:) )
       ENDDO
!$omp end do
       !
!$omp single
       CALL fft_scatter_many_xy ( dfft, f, aux, isgn, howmany )
!$omp end single
       !
!$omp do
       DO i = 0, howmany-1
          CALL cft_1z( aux(i*nnr_+1:), nsticks_x, n1, nx1, isgn, f(i*nnr_+1:) )
       ENDDO
!$omp end do
       !
!$omp do
       DO i = 0, howmany-1
          if (nsticks_x*nx1 < nnr_) then
             do j=nsticks_x*nx1+1, nnr_
                 f(j+i*nnr_) = (0.0_DP,0.0_DP)
             end do
          endif
       END DO
!$omp end do
!$omp end parallel
    ELSE                  ! R -> G
!$omp parallel default(none)                                        &
!$omp          private(i)                                           &
!$omp          shared(howmany, f, isgn, nnr_, dfft, aux)            &
!$omp          shared(nsticks_z, n3, nx3)                           &
!$omp          shared(nsticks_y, n2, nx2)                           &
!$omp          shared(nsticks_x, n1, nx1)                           &
!$omp          shared(nsticks_zx, nsticks_yx)
       !
!$omp do
       DO i = 0, howmany-1
          CALL cft_1z( f(i*nnr_+1:), nsticks_x, n1, nx1, isgn, aux(i*nnr_+1:) )
       ENDDO
!$omp end do
       !
!$omp single
       CALL fft_scatter_many_xy ( dfft, f, aux, isgn, howmany )
!$omp end single
       !
!$omp do
       DO i = 0, howmany-1
          CALL cft_1z( f(nx2*nsticks_yx*i+1:), nsticks_y, n2, nx2, isgn, aux(i*nnr_+1:))
       ENDDO
!$omp end do
       !
!$omp single
       CALL fft_scatter_many_yz( dfft, f, dfft%aux, isgn, howmany )
!$omp end single
       !
!$omp do
       DO i = 0, howmany-1
          CALL cft_1z( f(nx3*nsticks_zx*i+1:), nsticks_z, n3, nx3, isgn, aux(i*nnr_+1:) )
       ENDDO
!$omp end do
       !
!$omp do
       DO i = 0, howmany-1
          DO j=0, nsticks_z-1
             f(i*nnr_+j*nx3+1:i*nnr_+j*nx3+n3) = aux(i*nnr_+j*nx3+1:i*nnr_+j*nx3+n3)
          ENDDO
       ENDDO
!$omp end do
!$omp end parallel
    ENDIF
  ENDSUBROUTINE execute_using_threads
ENDSUBROUTINE many_cft3s


#if defined(__CUDA)
!
!  General purpose driver, GPU version
!
!----------------------------------------------------------------------------
SUBROUTINE tg_cft3s_gpu( f_d, dfft, isgn )
  !----------------------------------------------------------------------------
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +-3 : parallel 3d fft for wavefunctions with task group
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (cft_1z)
  !! ...              transpose across nodes           (fft_scatter_yz)
  !! ...              fft along y using pencils        (cft_1y)
  !! ...              transpose across nodes           (fft_scatter_xy)
  !! ...              fft along x using pencils        (cft_1x)
  !
  !! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  !! ...              fft along x using pencils        (cft_1x)
  !! ...              transpose across nodes           (fft_scatter_xy)
  !! ...              fft along y using pencils        (cft_1y)
  !! ...              transpose across nodes           (fft_scatter_yz)
  !! ...              fft along z using pencils        (cft_1z)
  !
  ! If task_group_fft_is_active the FFT acts on a number of wfcs equal to 
  ! dfft%nproc2, the number of Y-sections in which a plane is divided. 
  ! Data are reshuffled by the fft_scatter_tg routine so that each of the 
  ! dfft%nproc2 subgroups (made by dfft%nproc3 procs) deals with whole planes 
  ! of a single wavefunciton.
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  !
  USE cudafor
  !USE nvtx_fft
  USE fft_buffers,    ONLY : check_buffers_size, aux_d => dev_space_fftparallel
  USE fft_scalar,     ONLY : cft_1z_gpu
  USE fft_scatter_gpu,ONLY : fft_scatter_xy_gpu, fft_scatter_yz_gpu, fft_scatter_tg_gpu
  USE fft_scatter_gpu,ONLY : fft_scatter_tg_opt_gpu
  USE fft_types,  ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfft   ! descriptor of fft data layout
  COMPLEX(DP), DEVICE, INTENT(inout)     :: f_d( : ) ! array containing data to be transformed
  INTEGER, INTENT(in)                    :: isgn   ! fft direction (potential: +/-1, wave: +/-2, wave_tg: +/-3)
  !
  INTEGER                          :: n1, n2, n3, nx1, nx2, nx3
  INTEGER                          :: nnr_
  INTEGER                          :: nsticks_x, nsticks_y, nsticks_z
  !COMPLEX(DP), POINTER, DEVICE     :: aux_d (:)
  INTEGER                          :: ierr, i
  INTEGER(kind = cuda_stream_kind) :: stream  = 0
  !
  !write (6,*) 'enter tg_cft3s ',isgn ; write(6,*) ; FLUSH(6)
  n1  = dfft%nr1  ; n2  = dfft%nr2  ; n3  = dfft%nr3
  nx1 = dfft%nr1x ; nx2 = dfft%nr2x ; nx3 = dfft%nr3x
  !
  if (abs(isgn) == 1 ) then       ! potential fft
     nnr_ = dfft%nnr
     nsticks_x = dfft%my_nr2p * dfft%my_nr3p
     nsticks_y = dfft%nr1p(dfft%mype2+1) * dfft%my_nr3p
     nsticks_z = dfft%nsp(dfft%mype+1)
  else if (abs(isgn) == 2 ) then  ! wave func fft
     nnr_ = dfft%nnr
     nsticks_x = dfft%my_nr2p * dfft%my_nr3p
     nsticks_y = dfft%nr1w(dfft%mype2+1) * dfft%my_nr3p
     nsticks_z = dfft%nsw(dfft%mype+1)
  else if (abs(isgn) == 3 ) then  ! wave func fft with task groups
     nnr_ = dfft%nnr_tg
     nsticks_x = dfft%nr2 * dfft%my_nr3p
     nsticks_y = dfft%nr1w_tg * dfft%my_nr3p
     nsticks_z = dfft%nsw_tg(dfft%mype+1)
  else
     CALL fftx_error__( ' tg_cft3s', ' wrong value of isgn ', 10+abs(isgn) )
  end if
  !
  CALL check_buffers_size( dfft )
  !
  IF ( isgn > 0 ) THEN  ! G -> R
     !CALL nvtxStartRangeAsync("tg_cft3s_gpu G->R", 1)
     if (isgn==+3) then 
        call fft_scatter_tg_opt_gpu ( dfft, f_d, aux_d, nnr_, isgn, stream)
        CALL cft_1z_gpu( aux_d, nsticks_z, n3, nx3, isgn, f_d, stream )
     else
        !aux_d(1:nnr_)=f_d(1:nnr_) ! not limiting the range to dfft%nnr may crash when size(f_d)>size(aux_d)
        !ierr = cudaMemcpy( aux_d(1), f_d(1), nnr_, cudaMemcpyDeviceToDevice )
        CALL cft_1z_gpu( f_d, nsticks_z, n3, nx3, isgn, aux_d, stream, in_place=.true. )
     endif
     CALL fft_scatter_yz_gpu ( dfft, f_d, aux_d, nnr_, isgn )
     CALL cft_1z_gpu( aux_d, nsticks_y, n2, nx2, isgn, f_d, stream )
     CALL fft_scatter_xy_gpu ( dfft, f_d, aux_d, nnr_, isgn, stream )
     CALL cft_1z_gpu( aux_d, nsticks_x, n1, nx1, isgn, f_d, stream )
     ! clean garbage beyond the intended dimension. should not be needed but apparently it is !
     if (nsticks_x*nx1 < nnr_) then
        !$cuf kernel do(1)<<<*,*,0,stream>>>
        do i=nsticks_x*nx1+1, nnr_
            f_d(i) = (0.0_DP,0.0_DP)
        end do
     endif
     !CALL nvtxEndRangeAsync()
     !
  ELSE                  ! R -> G
     !
     !CALL nvtxStartRangeAsync("tg_cft3s_gpu R->G", 2)
     CALL cft_1z_gpu( f_d, nsticks_x, n1, nx1, isgn, aux_d, stream )
     CALL fft_scatter_xy_gpu ( dfft, f_d, aux_d, nnr_, isgn, stream )
     CALL cft_1z_gpu( f_d, nsticks_y, n2, nx2, isgn, aux_d, stream )
     CALL fft_scatter_yz_gpu ( dfft, f_d, aux_d, nnr_, isgn )

     if (isgn==-3) then
        CALL cft_1z_gpu( f_d, nsticks_z, n3, nx3, isgn, aux_d, stream )
        ! clean garbage beyond the intended dimension. should not be needed but apparently it is !
        if (nsticks_z*nx3 < nnr_) then
           !$cuf kernel do(1)<<<*,*,0,stream>>>
           do i=nsticks_z*nx3+1, nnr_
               aux_d(i) = (0.0_DP,0.0_DP)
           end do
        endif
        call fft_scatter_tg_opt_gpu ( dfft, aux_d, f_d, nnr_, isgn, stream)
     else
        !f_d(1:nnr_)=aux_d(1:nnr_) ! not limiting the range to dfft%nnr may crash when size(f_d)>size(aux_d)
        !ierr = cudaMemcpy( f_d(1), aux_d(1), nnr_, cudaMemcpyDeviceToDevice )
        CALL cft_1z_gpu( f_d, nsticks_z, n3, nx3, isgn, aux_d, stream, in_place=.true. )
        ! clean garbage beyond the intended dimension. should not be needed but apparently it is !
        if (nsticks_z*nx3 < nnr_) then
           !$cuf kernel do(1)<<<*,*,0,stream>>>
           do i=nsticks_z*nx3+1, nnr_
               f_d(i) = (0.0_DP,0.0_DP)
           end do
        endif
     endif
     !CALL nvtxEndRangeAsync()
  ENDIF
  !write (6,99) f_d(1:400); write(6,*); FLUSH(6)
  !
  !if (.true.) stop
  RETURN
99 format ( 20 ('(',2f12.9,')') )
  !
END SUBROUTINE tg_cft3s_gpu

!
!  Specific driver for the new 'many' call, GPU version
!
!----------------------------------------------------------------------------
SUBROUTINE many_cft3s_gpu( f_d, dfft, isgn, howmany )
  !----------------------------------------------------------------------------
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +-3 : parallel 3d fft for wavefunctions with task group
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (cft_1z)
  !! ...              transpose across nodes           (fft_scatter_yz)
  !! ...              fft along y using pencils        (cft_1y)
  !! ...              transpose across nodes           (fft_scatter_xy)
  !! ...              fft along x using pencils        (cft_1x)
  !
  !! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  !! ...              fft along x using pencils        (cft_1x)
  !! ...              transpose across nodes           (fft_scatter_xy)
  !! ...              fft along y using pencils        (cft_1y)
  !! ...              transpose across nodes           (fft_scatter_yz)
  !! ...              fft along z using pencils        (cft_1z)
  !
  ! If task_group_fft_is_active the FFT acts on a number of wfcs equal to 
  ! dfft%nproc2, the number of Y-sections in which a plane is divided. 
  ! Data are reshuffled by the fft_scatter_tg routine so that each of the 
  ! dfft%nproc2 subgroups (made by dfft%nproc3 procs) deals with whole planes 
  ! of a single wavefunciton.
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  !
  USE cudafor
  !USE nvtx_fft
  USE fft_buffers,    ONLY : check_buffers_size, aux_d => dev_space_fftparallel
  USE fft_scalar,     ONLY : cft_1z_gpu
  USE fft_scatter_gpu,ONLY : fft_scatter_xy_gpu, fft_scatter_yz_gpu, &
                              & fft_scatter_tg_gpu, fft_scatter_many_yz_gpu, &
                              & fft_scatter_tg_opt_gpu
  USE fft_types,  ONLY : fft_type_descriptor
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfft     ! descriptor of fft data layout
  COMPLEX(DP), DEVICE, INTENT(inout)     :: f_d( : ) ! array containing data to be transformed
  INTEGER, INTENT(in)                    :: isgn     ! fft direction (potential: +/-1, wave: +/-2, wave_tg: +/-3)
  INTEGER, INTENT(in)                    :: howmany  ! number of FFTs grouped together
  !
  INTEGER                          :: n1, n2, n3, nx1, nx2, nx3
  INTEGER                          :: nnr_
  INTEGER                          :: nsticks_x, nsticks_y, nsticks_z
  INTEGER                          :: nsticks_zx
  !COMPLEX(DP), POINTER, DEVICE     :: aux_d (:)
  INTEGER                          :: ierr, i, j
  INTEGER                          :: nstreams
  INTEGER(kind = cuda_stream_kind) :: stream  = 0 ! cuda_default_stream
  !
  !write (6,*) 'enter tg_cft3s ',isgn ; write(6,*) ; FLUSH(6)
  n1  = dfft%nr1  ; n2  = dfft%nr2  ; n3  = dfft%nr3
  nx1 = dfft%nr1x ; nx2 = dfft%nr2x ; nx3 = dfft%nr3x
  !
  if (abs(isgn) == 1 ) then       ! potential fft
     nnr_ = dfft%nnr
     nsticks_x = dfft%my_nr2p * dfft%my_nr3p
     nsticks_y = dfft%nr1p(dfft%mype2+1) * dfft%my_nr3p
     nsticks_z = dfft%nsp(dfft%mype+1)
     nsticks_zx = MAXVAL(dfft%nsp)
  else if (abs(isgn) == 2 ) then  ! wave func fft
     nnr_ = dfft%nnr
     nsticks_x = dfft%my_nr2p * dfft%my_nr3p
     nsticks_y = dfft%nr1w(dfft%mype2+1) * dfft%my_nr3p
     nsticks_z = dfft%nsw(dfft%mype+1)
     nsticks_zx = MAXVAL(dfft%nsw)
  else if (abs(isgn) == 3 ) then  ! wave func fft with task groups
     CALL fftx_error__( ' many_cft3s', ' Taskgroup and many not supported ', 10+abs(isgn) )
  else
     CALL fftx_error__( ' many_cft3s', ' wrong value of isgn ', 10+abs(isgn) )
  end if
  !
  CALL check_buffers_size( dfft, howmany )
  !
  nstreams = dfft%nstream_many
  ierr = cudaDeviceSynchronize()
  !
  IF ( isgn > 0 ) THEN  ! G -> R
     !CALL nvtxStartRangeAsync("many_cft3s_gpu G->R", 1)
     DO i = 0, howmany-1
        CALL cft_1z_gpu( f_d(i*nnr_+1:), nsticks_z, n3, nx3, isgn, aux_d(nx3*nsticks_zx*i+1:), &
                         stream ) !dfft%stream_many((mod(i,nstreams)+1)) )
     END DO
     !
     ! this brings back data from packed to unpacked
     CALL fft_scatter_many_yz_gpu ( dfft, aux_d(1), f_d(1), howmany*nnr_, isgn, howmany )
     !
     DO i = 0, howmany-1
        CALL cft_1z_gpu( f_d(i*nnr_+1:), nsticks_y, n2, nx2, isgn, aux_d(i*nnr_+1:), &
                         stream, in_place=.true. ) ! dfft%stream_many((mod(i,nstreams)+1)), in_place=.true. )
     END DO
     DO i = 0, howmany-1
        CALL fft_scatter_xy_gpu ( dfft, f_d(i*nnr_+1:), aux_d(i*nnr_+1:), nnr_, isgn, &
                                  stream ) ! dfft%stream_many((mod(i,nstreams)+1)) )
     END DO
     !
     DO i = 0, howmany-1
        CALL cft_1z_gpu( aux_d(i*nnr_+1:), nsticks_x, n1, nx1, isgn, f_d(i*nnr_+1:), &
                         stream ) ! dfft%stream_many((mod(i,nstreams)+1)) )
     END DO
     !
     ! clean garbage beyond the intended dimension. should not be needed but apparently it is !
     IF (nsticks_x*nx1 < nnr_) THEN
!$cuf kernel do(2)<<<*,*>>>
        DO i = 0, howmany-1
           DO j=nsticks_x*nx1+1, nnr_
               f_d(j+i*nnr_) = (0.0_DP,0.0_DP)
           END DO
        END DO
     ENDIF
     !
     !CALL nvtxEndRangeAsync()
     !
  ELSE                  ! R -> G
     !
     !CALL nvtxStartRangeAsync("many_cft3s_gpu R->G", 2)
     DO i = 0, howmany-1
        CALL cft_1z_gpu( f_d(i*nnr_+1:), nsticks_x, n1, nx1, isgn, aux_d(i*nnr_+1:), &
                         stream ) ! dfft%stream_many((mod(i,nstreams)+1)) )
     END DO
     DO i = 0, howmany-1
        CALL fft_scatter_xy_gpu ( dfft, f_d(i*nnr_+1), aux_d(i*nnr_+1), nnr_, isgn, &
                                  stream ) ! dfft%stream_many((mod(i,nstreams)+1)) )
     END DO
     !
     DO i = 0, howmany-1
        CALL cft_1z_gpu( f_d(i*nnr_+1:), nsticks_y, n2, nx2, isgn, aux_d(i*nnr_+1:), &
                         stream, in_place=.true. ) ! dfft%stream_many((mod(i,nstreams)+1)), in_place=.true. )
     END DO
     !
     CALL fft_scatter_many_yz_gpu ( dfft, aux_d, f_d, howmany*nnr_, isgn, howmany )
     !
     DO i = 0, howmany-1
        CALL cft_1z_gpu( aux_d(nx3*nsticks_zx*i+1:), nsticks_z, n3, nx3, isgn, f_d(i*nnr_+1:), &
                         stream ) ! dfft%stream_many((mod(i,nstreams)+1)) )
     END DO
     !
     ! clean garbage beyond the intended dimension. should not be needed but apparently it is !
     IF (nsticks_z*nx3 < nnr_) THEN
!$cuf kernel do(2)<<<*,*>>>
        DO i = 0, howmany-1
           DO j=nsticks_z*nx3+1, nnr_
                f_d(j+i*nnr_) = (0.0_DP,0.0_DP)
           END DO
        END DO
     ENDIF
     !CALL nvtxEndRangeAsync()
  ENDIF
  !write (6,99) f_d(1:400); write(6,*); FLUSH(6)
  !
  !DO i=1,nstreams
  !    ierr = cudaStreamSynchronize( dfft%stream_many(i) )
  !END DO
  !
  RETURN
99 format ( 20 ('(',2f12.9,')') )
  !
END SUBROUTINE many_cft3s_gpu

#endif
!--------------------------------------------------------------------------------
!   Auxiliary routines to read/write from/to a distributed array
!   NOT optimized for efficiency .... just to show how one can access the data
!--------------------------------------------------------------------------------
!
COMPLEX (DP) FUNCTION get_f_of_R (i,j,k,f,dfft)
!------  read from a distributed complex array f(:) in direct space
!
  USE fft_types,  ONLY : fft_type_descriptor
  IMPLICIT NONE
  TYPE (fft_type_descriptor), INTENT(IN) :: dfft
  INTEGER, INTENT (IN) :: i,j,k
  COMPLEX(DP), INTENT (IN) :: f(:)
  INTEGER :: kk, ii, jj, ip, ierr
  COMPLEX(DP) :: f_aux

  IF ( i <= 0 .OR. i > dfft%nr1 ) CALL fftx_error__( ' get_f_of_R', ' first  index out of range ', 1 )
  IF ( j <= 0 .OR. j > dfft%nr2 ) CALL fftx_error__( ' get_f_of_R', ' second index out of range ', 2 )
  IF ( k <= 0 .OR. k > dfft%nr3 ) CALL fftx_error__( ' get_f_of_R', ' third  index out of range ', 3 )

#if defined(__MPI)
  do ip = 1, dfft%nproc3
     if ( dfft%i0r3p(ip) < k ) kk = ip
  end do
  do ip = 1, dfft%nproc2
     if ( dfft%i0r2p(ip) < j ) jj = ip
  end do
  ii  = i + dfft%nr1x * ( j - dfft%i0r2p(jj) - 1 ) + dfft%nr1x * dfft%nr2p(jj) * ( k - dfft%i0r3p(kk) - 1 )
  f_aux = (0.d0,0.d0)
  if ( (jj == (dfft%mype2 + 1)) .and. (kk == (dfft%mype3 +1)) ) f_aux = f(ii)
  CALL MPI_ALLREDUCE( f_aux, get_f_of_R,   2, MPI_DOUBLE_PRECISION, MPI_SUM, dfft%comm, ierr )
#else
  ii = i + dfft%nr1x * (j-1) + dfft%nr1x*dfft%nr2x * (k-1)
  get_f_of_R = f(ii)
#endif
END FUNCTION get_f_of_R

SUBROUTINE put_f_of_R (f_in,i,j,k,f,dfft)
!------  write on a distributed complex array f(:) in direct space
!
  USE fft_types,  ONLY : fft_type_descriptor
  IMPLICIT NONE
  TYPE (fft_type_descriptor), INTENT(IN) :: dfft
  INTEGER, INTENT (IN) :: i,j,k
  COMPLEX(DP), INTENT (IN) :: f_in
  COMPLEX(DP), INTENT (INOUT) :: f(:)
  INTEGER :: kk, ii, jj, ip, ierr

  IF ( i <= 0 .OR. i > dfft%nr1 ) CALL fftx_error__( ' put_f_of_R', ' first  index out of range ', 1 )
  IF ( j <= 0 .OR. j > dfft%nr2 ) CALL fftx_error__( ' put_f_of_R', ' second index out of range ', 2 )
  IF ( k <= 0 .OR. k > dfft%nr3 ) CALL fftx_error__( ' put_f_of_R', ' third  index out of range ', 3 )

#if defined(__MPI)
  do ip = 1, dfft%nproc3
     if ( dfft%i0r3p(ip) < k ) kk = ip
  end do
  do ip = 1, dfft%nproc2
     if ( dfft%i0r2p(ip) < j ) jj = ip
  end do
  ii  = i + dfft%nr1x * ( j - dfft%i0r2p(jj) - 1 ) + dfft%nr1x * dfft%nr2p(jj) * ( k - dfft%i0r3p(kk) - 1 )
  if ( (jj == (dfft%mype2 + 1)) .and. (kk == (dfft%mype3 +1)) ) f(ii) = f_in
#else
  ii = i + dfft%nr1x * (j-1) + dfft%nr1x*dfft%nr2x * (k-1)
  f(ii) = f_in
#endif

END SUBROUTINE put_f_of_R

COMPLEX (DP) FUNCTION get_f_of_G (i,j,k,f,dfft)
!------  read from a distributed complex array f(:) in reciprocal space
!
  USE fft_types,  ONLY : fft_type_descriptor
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: i,j,k
  COMPLEX(DP), INTENT (IN) :: f(:)
  TYPE (fft_type_descriptor), INTENT(IN) :: dfft
  INTEGER :: ii, jj, ip, ierr
  COMPLEX(DP) :: f_aux

  IF ( i <= 0 .OR. i > dfft%nr1 ) CALL fftx_error__( ' get_f_of_G', ' first  index out of range ', 1 )
  IF ( j <= 0 .OR. j > dfft%nr2 ) CALL fftx_error__( ' get_f_of_G', ' second index out of range ', 2 )
  IF ( k <= 0 .OR. k > dfft%nr3 ) CALL fftx_error__( ' get_f_of_G', ' third  index out of range ', 3 )

#if defined(__MPI)
  ii = i + dfft%nr1x * (j -1)
  jj = dfft%isind(ii)  ! if jj is zero this G vector does not belong to this processor
  f_aux = (0.d0,0.d0)
  if ( jj > 0 ) f_aux = f( k + dfft%nr3x * (jj -1))
  CALL MPI_ALLREDUCE( f_aux, get_f_of_G,   2, MPI_DOUBLE_PRECISION, MPI_SUM, dfft%comm, ierr )
#else
  ii = i + dfft%nr1 * (j-1) + dfft%nr1*dfft%nr2 * (k-1)
  get_f_of_G = f(ii)
#endif
END FUNCTION get_f_of_G

SUBROUTINE put_f_of_G (f_in,i,j,k,f,dfft)
!------  write on a distributed complex array f(:) in reciprocal space
!
  USE fft_types,  ONLY : fft_type_descriptor
  IMPLICIT NONE
  COMPLEX(DP), INTENT (IN) :: f_in
  INTEGER, INTENT (IN) :: i,j,k
  COMPLEX(DP), INTENT (INOUT) :: f(:)
  TYPE (fft_type_descriptor), INTENT(IN) :: dfft
  INTEGER :: ii, jj

  IF ( i <= 0 .OR. i > dfft%nr1 ) CALL fftx_error__( ' put_f_of_G', ' first  index out of range ', 1 )
  IF ( j <= 0 .OR. j > dfft%nr2 ) CALL fftx_error__( ' put_f_of_G', ' second index out of range ', 2 )
  IF ( k <= 0 .OR. k > dfft%nr3 ) CALL fftx_error__( ' put_f_of_G', ' third  index out of range ', 3 )

#if defined(__MPI)
  ii = i + dfft%nr1x * (j -1)
  jj = dfft%isind(ii)  ! if jj is zero this G vector does not belong to this processor
  if ( jj > 0 )   f( k + dfft%nr3x * (jj -1)) = f_in
#else
  ii = i + dfft%nr1 * (j-1) + dfft%nr1*dfft%nr2 * (k-1)
  f(ii) = f_in
#endif
END SUBROUTINE put_f_of_G

END MODULE fft_parallel
