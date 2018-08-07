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
!
!     Parallel 3D FFT high level Driver
!     ( Charge density and Wave Functions )
!
!     Written and maintained by Carlo Cavazzoni
!     Last update Apr. 2009
!
!!=---------------------------------------------------------------------==!
!
MODULE fft_parallel_2d
!
#ifdef __CUDA
   USE cudafor
#endif

   USE fft_param
   IMPLICIT NONE
   SAVE

#ifdef __CUDA
  COMPLEX(DP), PINNED, ALLOCATABLE :: f_h(:), aux_h(:), aux2_h(:)
  COMPLEX(DP), DEVICE, ALLOCATABLE :: aux_d(:), aux2_d(:)
#endif



!
CONTAINS
!
!  General purpose driver, including Task groups parallelization
!
!----------------------------------------------------------------------------
SUBROUTINE tg_cft3s( f, dfft, isgn )
  !----------------------------------------------------------------------------
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !                  NOT IMPLEMENTED WITH TASK GROUPS
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (cft_1z)
  !! ...              transpose across nodes           (fft_scatter)
  !! ...                 and reorder
  ! ...              fft along y (using planes) and x (cft_2xy)
  ! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  ! ...              fft along x and y(using planes)  (cft_2xy)
  ! ...              transpose across nodes           (fft_scatter)
  ! ...                 and reorder
  ! ...              fft along z using pencils        (cft_1z)
  !
  ! ...  The array "planes" signals whether a fft is needed along y :
  ! ...    planes(i)=0 : column f(i,*,*) empty , don't do fft along y
  ! ...    planes(i)=1 : column f(i,*,*) filled, fft along y needed
  ! ...  "empty" = no active components are present in f(i,*,*)
  ! ...            after (isgn>0) or before (isgn<0) the fft on z direction
  !
  ! ...  Note that if isgn=+/-1 (fft on rho and pot.) all fft's are needed
  ! ...  and all planes(i) are set to 1
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  ! Task Group added by Costas Bekas, Oct. 2005, adapted from the CPMD code
  ! (Alessandro Curioni) and revised by Carlo Cavazzoni 2007.
  !
  USE fft_scalar, ONLY : cft_1z, cft_2xy
  USE scatter_mod_2d,   ONLY : fft_scatter
  USE fft_types,  ONLY : fft_type_descriptor

  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout)    :: f( : )  ! array containing data to be transformed
  TYPE (fft_type_descriptor), INTENT(in) :: dfft
                                           ! descriptor of fft data layout
  INTEGER, INTENT(in)           :: isgn    ! fft direction
  !
  ! the following ifdef prevents usage of directive in older ifort versions
  !
#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
  ! the following is a workaround for Intel 12.1 bug
#if __INTEL_COMPILER  < 9999
!dir$ attributes align: 4096 :: aux
#endif
#endif
#endif
  INTEGER                    :: me_p
  INTEGER                    :: n1, n2, n3, nx1, nx2, nx3
  COMPLEX(DP), ALLOCATABLE   :: aux (:)
  INTEGER                    :: planes( dfft%nr1x )
  !LOGICAL                    :: use_tg
  !
  !
  !use_tg = dfft%has_task_groups
  IF (dfft%has_task_groups) CALL fftx_error__( ' tg_cft3s ', ' task groups on large mesh not implemented ', 1 )
  !
  n1  = dfft%nr1
  n2  = dfft%nr2
  n3  = dfft%nr3
  nx1 = dfft%nr1x
  nx2 = dfft%nr2x
  nx3 = dfft%nr3x
  !
  ALLOCATE( aux( dfft%nnr ) )
  !
  me_p = dfft%mype + 1
  !
  IF ( isgn > 0 ) THEN
     !
     IF ( isgn /= 2 ) THEN
        !
        CALL cft_1z( f, dfft%nsp( me_p ), n3, nx3, isgn, aux )
        !
        planes = dfft%iplp
        !
     ELSE
        !
        CALL cft_1z( f, dfft%nsw( me_p ), n3, nx3, isgn, aux )
        !
        planes = dfft%iplw
        !
     ENDIF
     !
     CALL fw_scatter( isgn ) ! forward scatter from stick to planes
     !
     CALL cft_2xy( f, dfft%nr3p( me_p ), n1, n2, nx1, nx2, isgn, planes )
     !
  ELSE
     !
     IF ( isgn == -1 ) THEN
        !
        planes = dfft%iplp
        !
     ELSE IF ( isgn == -2 ) THEN
        !
        planes = dfft%iplw
        !
     ENDIF
     !
     CALL cft_2xy( f, dfft%nr3p( me_p ), n1, n2, nx1, nx2, isgn, planes)
     !
     CALL bw_scatter( isgn )
     !
     IF ( isgn /= -2 ) THEN
        !
        CALL cft_1z( aux, dfft%nsp( me_p ), n3, nx3, isgn, f )
        !
     ELSE
        !
        CALL cft_1z( aux, dfft%nsw( me_p ), n3, nx3, isgn, f )
        !
     ENDIF
     !
  ENDIF
  !
  DEALLOCATE( aux )
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE fw_scatter( iopt )

        !Transpose data for the 2-D FFT on the x-y plane
        !
        !NOGRP*dfft%nnr: The length of aux and f
        !nr3x: The length of each Z-stick
        !aux: input - output
        !f: working space
        !isgn: type of scatter
        !dfft%nsw(me) holds the number of Z-sticks proc. me has.
        !dfft%nr3p: number of planes per processor
        !
     !
     USE scatter_mod_2d, ONLY: fft_scatter
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == 2 ) THEN
        !
        CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsw, dfft%nr3p, iopt )
        !
     ELSEIF( iopt == 1 ) THEN
        !
        CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsp, dfft%nr3p, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE fw_scatter

  !

  SUBROUTINE bw_scatter( iopt )
     !
     USE scatter_mod_2d, ONLY: fft_scatter
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == -2 ) THEN
        !
        CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsw, dfft%nr3p, iopt )
        !
     ELSEIF( iopt == -1 ) THEN
        !
        CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsp, dfft%nr3p, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE bw_scatter
  !
END SUBROUTINE tg_cft3s
!
!
!

#if defined(__CUDA)
!----------------------------------------------------------------------------
SUBROUTINE tg_cft3s_gpu( f_d, dfft, isgn )
  !----------------------------------------------------------------------------
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !                  NOT IMPLEMENTED WITH TASK GROUPS
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (cft_1z)
  !! ...              transpose across nodes           (fft_scatter)
  !! ...                 and reorder
  ! ...              fft along y (using planes) and x (cft_2xy)
  ! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  ! ...              fft along x and y(using planes)  (cft_2xy)
  ! ...              transpose across nodes           (fft_scatter)
  ! ...                 and reorder
  ! ...              fft along z using pencils        (cft_1z)
  !
  ! ...  The array "planes" signals whether a fft is needed along y :
  ! ...    planes(i)=0 : column f(i,*,*) empty , don't do fft along y
  ! ...    planes(i)=1 : column f(i,*,*) filled, fft along y needed
  ! ...  "empty" = no active components are present in f(i,*,*)
  ! ...            after (isgn>0) or before (isgn<0) the fft on z direction
  !
  ! ...  Note that if isgn=+/-1 (fft on rho and pot.) all fft's are needed
  ! ...  and all planes(i) are set to 1
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  ! Task Group added by Costas Bekas, Oct. 2005, adapted from the CPMD code
  ! (Alessandro Curioni) and revised by Carlo Cavazzoni 2007.
  !
  USE fft_scalar, ONLY : cft_1z_gpu, cft_2xy_gpu
  USE scatter_mod_2d_gpu,   ONLY : fft_scatter_gpu
  USE fft_types,  ONLY : fft_type_descriptor

  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfft
  COMPLEX(DP), DEVICE, INTENT(inout)    :: f_d( dfft%nnr ) ! array containing data to be transformed
                                           ! descriptor of fft data layout
  INTEGER, INTENT(in)           :: isgn    ! fft direction
  !
  ! the following ifdef prevents usage of directive in older ifort versions
  !
#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
  ! the following is a workaround for Intel 12.1 bug
#if __INTEL_COMPILER  < 9999
!dir$ attributes align: 4096 :: yf, aux
#endif
#endif
#endif
  INTEGER                    :: me_p, istat
  INTEGER                    :: n1, n2, n3, nx1, nx2, nx3
  COMPLEX(DP), ALLOCATABLE   :: yf(:)
  INTEGER                    :: planes( dfft%nr1x )
  INTEGER(kind = cuda_stream_kind) :: stream  = 0
  !
  !
  n1  = dfft%nr1
  n2  = dfft%nr2
  n3  = dfft%nr3
  nx1 = dfft%nr1x
  nx2 = dfft%nr2x
  nx3 = dfft%nr3x
  !
  IF( dfft%has_task_groups ) CALL fftx_error__( ' tg_cft3s ', ' task groups in 2D + 1D decomposition not implemented ', 1 )

  IF( ALLOCATED( f_h   ) .and. SIZE( f_h   ) < dfft%nnr ) DEALLOCATE( f_h   )
  IF( ALLOCATED( aux_h ) .and. SIZE( aux_h ) < dfft%nnr ) DEALLOCATE( aux_h )
  IF( ALLOCATED( aux_d ) .and. SIZE( aux_d ) < dfft%nnr ) DEALLOCATE( aux_d )

  IF( .not. ALLOCATED( f_h   ) ) ALLOCATE( f_h  ( dfft%nnr ) )
  IF( .not. ALLOCATED( aux_h ) ) ALLOCATE( aux_h( dfft%nnr ) )
  IF( .not. ALLOCATED( aux_d ) ) ALLOCATE( aux_d( dfft%nnr ) )

  !
  me_p = dfft%mype + 1
  !
  IF ( isgn > 0 ) THEN
     !
     IF ( isgn /= 2 ) THEN
        !
        CALL cft_1z_gpu( f_d, dfft%nsp( me_p ), n3, nx3, isgn, aux_d, stream )
        !
        planes = dfft%iplp
        !
     ELSE
        !
        CALL cft_1z_gpu( f_d, dfft%nsw( me_p ), n3, nx3, isgn, aux_d, stream )
        !
        planes = dfft%iplw
        !
     ENDIF
     !
     CALL fw_scatter_gpu( isgn ) ! forward scatter from stick to planes
     !
     CALL cft_2xy_gpu( f_d, aux_d, dfft%nr3p( me_p ), n1, n2, nx1, nx2, isgn, stream, planes )
     !
  ELSE
     !
     IF ( isgn /= -2 ) THEN
        !
        planes = dfft%iplp
        !
     ELSE
        !
        planes = dfft%iplw
        !
     ENDIF

     CALL cft_2xy_gpu( f_d, aux_d, dfft%nr3p( me_p ), n1, n2, nx1, nx2, isgn, stream, planes)
     !
     CALL bw_scatter_gpu( isgn )
     !
     !f_d = (0.d0, 0.d0)
     !
     IF ( isgn /= -2 ) THEN
        !
        CALL cft_1z_gpu( aux_d, dfft%nsp( me_p ), n3, nx3, isgn, f_d, stream )
         !
     ELSE
        !
        CALL cft_1z_gpu( aux_d, dfft%nsw( me_p ), n3, nx3, isgn, f_d, stream )
        !
     ENDIF
     !
  ENDIF
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE fw_scatter_gpu( iopt )

        !Transpose data for the 2-D FFT on the x-y plane
        !
        !NOGRP*dfft%nnr: The length of aux and f
        !nr3x: The length of each Z-stick
        !aux: input - output
        !f: working space
        !isgn: type of scatter
        !dfft%nsw(me) holds the number of Z-sticks proc. me has.
        !dfft%nr3p: number of planes per processor
        !
     !
     USE scatter_mod_2d_gpu, ONLY: fft_scatter_gpu
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == 2 ) THEN
        !
        CALL fft_scatter_gpu( dfft, aux_d, aux_h, nx3, dfft%nnr, f_d, f_h, dfft%nsw, dfft%nr3p, iopt )
        !
     ELSEIF( iopt == 1 ) THEN
        !
        CALL fft_scatter_gpu( dfft, aux_d, aux_h, nx3, dfft%nnr, f_d, f_h, dfft%nsp, dfft%nr3p, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE fw_scatter_gpu

  !

  SUBROUTINE bw_scatter_gpu( iopt )
     !
     USE scatter_mod_2d_gpu, ONLY: fft_scatter_gpu
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == -2 ) THEN
        !
        CALL fft_scatter_gpu( dfft, aux_d, aux_h, nx3, dfft%nnr, f_d, f_h, dfft%nsw, dfft%nr3p, iopt )
        !
     ELSEIF( iopt == -1 ) THEN
        !
        CALL fft_scatter_gpu( dfft, aux_d, aux_h, nx3, dfft%nnr, f_d, f_h, dfft%nsp, dfft%nr3p, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE bw_scatter_gpu
  !
END SUBROUTINE tg_cft3s_gpu

SUBROUTINE tg_cft3s_batch_gpu( f_d, dfft, isgn, batchsize )
  !----------------------------------------------------------------------------
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !                  NOT IMPLEMENTED WITH TASK GROUPS
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (cft_1z)
  !! ...              transpose across nodes           (fft_scatter)
  !! ...                 and reorder
  ! ...              fft along y (using planes) and x (cft_2xy)
  ! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  ! ...              fft along x and y(using planes)  (cft_2xy)
  ! ...              transpose across nodes           (fft_scatter)
  ! ...                 and reorder
  ! ...              fft along z using pencils        (cft_1z)
  !
  ! ...  The array "planes" signals whether a fft is needed along y :
  ! ...    planes(i)=0 : column f(i,*,*) empty , don't do fft along y
  ! ...    planes(i)=1 : column f(i,*,*) filled, fft along y needed
  ! ...  "empty" = no active components are present in f(i,*,*)
  ! ...            after (isgn>0) or before (isgn<0) the fft on z direction
  !
  ! ...  Note that if isgn=+/-1 (fft on rho and pot.) all fft's are needed
  ! ...  and all planes(i) are set to 1
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  ! Task Group added by Costas Bekas, Oct. 2005, adapted from the CPMD code
  ! (Alessandro Curioni) and revised by Carlo Cavazzoni 2007.
  !
  USE fft_scalar, ONLY : cft_1z_gpu, cft_2xy_gpu
  USE scatter_mod_2d_gpu,   ONLY : fft_scatter_gpu, fft_scatter_batch_a_gpu, fft_scatter_batch_b_gpu
  USE fft_types,  ONLY : fft_type_descriptor

  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfft
                                           ! descriptor of fft data layout
  INTEGER, INTENT(in)           :: isgn    ! fft direction
  INTEGER, INTENT(in)           :: batchsize   
  COMPLEX(DP), DEVICE, INTENT(inout)    :: f_d( batchsize * dfft%nnr ) ! array containing data to be transformed
  !
  INTEGER                    :: me_p, istat, i, j, currsize
  INTEGER                    :: n1, n2, n3, nx1, nx2, nx3, ncpx, nppx, proc
  COMPLEX(DP), ALLOCATABLE   :: yf(:)
  INTEGER                    :: planes( dfft%nr1x )
  INTEGER(kind = cuda_stream_kind) :: stream  = 0
  !
  !
  n1  = dfft%nr1
  n2  = dfft%nr2
  n3  = dfft%nr3
  nx1 = dfft%nr1x
  nx2 = dfft%nr2x
  nx3 = dfft%nr3x
  !
  IF( ALLOCATED( f_h   ) .and. SIZE( f_h   ) < batchsize * dfft%nnr ) DEALLOCATE( f_h   )
  IF( ALLOCATED( aux_h ) .and. SIZE( aux_h ) < batchsize * dfft%nnr ) DEALLOCATE( aux_h )
  IF( ALLOCATED( aux_d ) .and. SIZE( aux_d ) < batchsize * dfft%nnr ) DEALLOCATE( aux_d )
  IF( ALLOCATED( aux2_h ) .and. SIZE( aux2_h ) < batchsize * dfft%nnr ) DEALLOCATE( aux2_h )
  IF( ALLOCATED( aux2_d ) .and. SIZE( aux2_d ) < batchsize * dfft%nnr ) DEALLOCATE( aux2_d )

  IF( .not. ALLOCATED( f_h   ) ) ALLOCATE( f_h  ( batchsize * dfft%nnr ) )
  IF( .not. ALLOCATED( aux_h ) ) ALLOCATE( aux_h( batchsize * dfft%nnr ) )
  IF( .not. ALLOCATED( aux_d ) ) ALLOCATE( aux_d( batchsize * dfft%nnr ) )
  IF( .not. ALLOCATED( aux2_h ) ) ALLOCATE( aux2_h( batchsize * dfft%nnr ) )
  IF( .not. ALLOCATED( aux2_d ) ) ALLOCATE( aux2_d( batchsize * dfft%nnr ) )
  !
  me_p = dfft%mype + 1

  ncpx = 0
  nppx = 0
  DO proc = 1, dfft%nproc
     ncpx = max( ncpx, dfft%nsw ( proc ) )
     nppx = max( nppx, dfft%nr3p ( proc ) )
  ENDDO
  !
  IF ( isgn > 0 ) THEN
     DO j = 0, batchsize-1, dfft%subbatchsize
       currsize = min(dfft%subbatchsize, batchsize - j)
       !
       IF ( isgn /= 2 ) THEN
          !
          CALL fftx_error__( ' tg_cft3s_gpu_batch ', ' isgn /= 2 not implemented ', 1 )
          !
          CALL cft_1z_gpu( f_d, dfft%nsp( me_p ), n3, nx3, isgn, aux_d, stream )
          !
          planes = dfft%iplp
          !
       ELSE
          !

          DO i = 0, currsize - 1
            CALL cft_1z_gpu( f_d((j+i)*dfft%nnr + 1:), dfft%nsw( me_p ), n3, nx3, isgn, aux_d(j*dfft%nnr + i*ncpx*nx3 +1:), dfft%a2a_comp )
          ENDDO

          i = cudaEventRecord(dfft%bevents(j/dfft%subbatchsize + 1), dfft%a2a_comp)
          i = cudaStreamWaitEvent(dfft%bstreams(j/dfft%subbatchsize + 1), dfft%bevents(j/dfft%subbatchsize + 1), 0)

          !
          planes = dfft%iplw
          !
       ENDIF

       IF (j > 0) i = cudaStreamWaitEvent(dfft%bstreams(j/dfft%subbatchsize + 1), dfft%bevents(j/dfft%subbatchsize), 0)

       IF (dfft%nproc .ne. 1) THEN
         CALL fft_scatter_batch_a_gpu( dfft, aux_d(j*dfft%nnr + 1:), aux_h(j*dfft%nnr + 1:), nx3, dfft%nnr, f_d(j*dfft%nnr + 1:), &
         f_h(j*dfft%nnr + 1:), aux2_d(j*dfft%nnr + 1:), aux2_h(j*dfft%nnr + 1:), dfft%nsw, dfft%nr3p, isgn, currsize, j/dfft%subbatchsize + 1 )
       ELSE
         CALL fft_scatter_batch_a_gpu( dfft, aux2_d(j*dfft%nnr + 1:), aux2_h(j*dfft%nnr + 1:), nx3, dfft%nnr, f_d(j*dfft%nnr + 1:), &
         f_h(j*dfft%nnr + 1:), aux_d(j*dfft%nnr + 1:), aux_h(j*dfft%nnr + 1:), dfft%nsw, dfft%nr3p, isgn, currsize, j/dfft%subbatchsize + 1 )
       ENDIF

     ENDDO

     DO j = 0, batchsize-1, dfft%subbatchsize
       currsize = min(dfft%subbatchsize, batchsize - j)

       IF (dfft%nproc .ne. 1) THEN
         CALL fft_scatter_batch_b_gpu( dfft, aux_d(j*dfft%nnr + 1:), aux_h(j*dfft%nnr + 1:), nx3, dfft%nnr, f_d(j*dfft%nnr + 1:), &
         f_h(j*dfft%nnr + 1:), aux2_d(j*dfft%nnr + 1:), aux2_h(j*dfft%nnr + 1:), dfft%nsw, dfft%nr3p, isgn, currsize, j/dfft%subbatchsize + 1 )
       ELSE
         CALL fft_scatter_batch_b_gpu( dfft, aux2_d(j*dfft%nnr + 1:), aux2_h(j*dfft%nnr + 1:), nx3, dfft%nnr, f_d(j*dfft%nnr + 1:), &
         f_h(j*dfft%nnr + 1:), aux_d(j*dfft%nnr + 1:), aux_h(j*dfft%nnr + 1:), dfft%nsw, dfft%nr3p, isgn, currsize, j/dfft%subbatchsize + 1 )
       ENDIF

       IF (currsize == dfft%subbatchsize) THEN
         CALL cft_2xy_gpu( f_d(j*dfft%nnr + 1:), aux_d(j*dfft%nnr + 1:), currsize * nppx, n1, n2, nx1, nx2, isgn, dfft%a2a_comp, planes )
       ELSE
         DO i = 0, currsize - 1
           CALL cft_2xy_gpu( f_d((j+i)*dfft%nnr + 1:), aux_d((j+i)*dfft%nnr + 1:), dfft%nr3p( me_p ), n1, n2, nx1, nx2, isgn,  &
           dfft%a2a_comp, planes )
         ENDDO
       ENDIF

     ENDDO

!     i = cudaDeviceSynchronize()

     !
  ELSE
!     i = cudaDeviceSynchronize()

     DO j = 0, batchsize-1, dfft%subbatchsize
       currsize = min(dfft%subbatchsize, batchsize - j)
       !
       IF ( isgn /= -2 ) THEN
          !
          CALL fftx_error__( ' tg_cft3s_gpu_batch ', ' isgn /= -2 not implemented ', 1 )
          !
          planes = dfft%iplp
          !
       ELSE
          !
          planes = dfft%iplw
          !
       ENDIF

       IF (currsize == dfft%subbatchsize) THEN
         CALL cft_2xy_gpu( f_d(j*dfft%nnr + 1:), aux_d(j*dfft%nnr + 1:), currsize * nppx, n1, n2, nx1, nx2, isgn, dfft%a2a_comp, planes )
       ELSE
         DO i = 0, currsize - 1
           CALL cft_2xy_gpu( f_d((j+i)*dfft%nnr + 1:), aux_d((j+i)*dfft%nnr + 1:), dfft%nr3p( me_p ), n1, n2, nx1, nx2, isgn, dfft%a2a_comp, planes )
         ENDDO
       ENDIF

       if (j > 0) i = cudaStreamWaitEvent(dfft%bstreams(j/dfft%subbatchsize + 1), dfft%bevents(j/dfft%subbatchsize), 0)

       IF (dfft%nproc .ne. 1) THEN
         CALL fft_scatter_batch_a_gpu( dfft, aux_d(j*dfft%nnr + 1:), aux_h(j*dfft%nnr + 1:), nx3, dfft%nnr, f_d(j*dfft%nnr + 1:), &
         f_h(j*dfft%nnr + 1:), aux2_d(j*dfft%nnr + 1:), aux2_h(j*dfft%nnr + 1:), dfft%nsw, dfft%nr3p, isgn, currsize, j/dfft%subbatchsize + 1 )
       ELSE
         CALL fft_scatter_batch_a_gpu( dfft, aux2_d(j*dfft%nnr + 1:), aux2_h(j*dfft%nnr + 1:), nx3, dfft%nnr, f_d(j*dfft%nnr + 1:), &
         f_h(j*dfft%nnr + 1:), aux_d(j*dfft%nnr + 1:), aux_h(j*dfft%nnr + 1:), dfft%nsw, dfft%nr3p, isgn, currsize, j/dfft%subbatchsize + 1 )
       ENDIF

     ENDDO

     DO j = 0, batchsize-1, dfft%subbatchsize
       currsize = min(dfft%subbatchsize, batchsize - j)

       IF (dfft%nproc .ne. 1) THEN
         CALL fft_scatter_batch_b_gpu( dfft, aux_d(j*dfft%nnr + 1:), aux_h(j*dfft%nnr + 1:), nx3, dfft%nnr, f_d(j*dfft%nnr + 1:), &
         f_h(j*dfft%nnr + 1:), aux2_d(j*dfft%nnr + 1:), aux2_h(j*dfft%nnr + 1:), dfft%nsw, dfft%nr3p, isgn, currsize, j/dfft%subbatchsize + 1 )
       ELSE
         CALL fft_scatter_batch_b_gpu( dfft, aux2_d(j*dfft%nnr + 1:), aux2_h(j*dfft%nnr + 1:), nx3, dfft%nnr, f_d(j*dfft%nnr + 1:), &
         f_h(j*dfft%nnr + 1:), aux_d(j*dfft%nnr + 1:), aux_h(j*dfft%nnr + 1:), dfft%nsw, dfft%nr3p, isgn, currsize, j/dfft%subbatchsize + 1 )
       ENDIF

       i = cudaEventRecord(dfft%bevents(j/dfft%subbatchsize + 1), dfft%bstreams(j/dfft%subbatchsize + 1))
       i = cudaStreamWaitEvent(dfft%a2a_comp, dfft%bevents(j/dfft%subbatchsize + 1), 0)


       DO i = 0, currsize - 1
         CALL cft_1z_gpu( aux_d(j*dfft%nnr + i*ncpx*nx3 + 1:), dfft%nsw( me_p ), n3, nx3, isgn, f_d((j+i)*dfft%nnr + 1:), dfft%a2a_comp )
       ENDDO

     ENDDO

!     i = cudaDeviceSynchronize()

  ENDIF
  !
  !
  RETURN
  !
  !
END SUBROUTINE tg_cft3s_batch_gpu
#endif
!


END MODULE fft_parallel_2d
