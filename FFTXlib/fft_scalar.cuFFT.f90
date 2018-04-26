!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!=----------------------------------------------------------------------=!
   MODULE fft_scalar_cuFFT
!=----------------------------------------------------------------------=!
#ifdef __CUDA
       USE fft_param
!! iso_c_binding provides C_PTR, C_NULL_PTR, C_ASSOCIATED
       USE iso_c_binding

       USE cudafor
       USE cufft
      
       IMPLICIT NONE
       SAVE
       PRIVATE
       PUBLIC :: cft_1z_gpu, cft_2xy_gpu, cfft3d_gpu, cfft3ds_gpu

!=----------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------=!


   SUBROUTINE cft_1z_gpu(c_d, nsl, nz, ldz, isign, cout_d, stream, in_place)

!     driver routine for nsl 1d complex fft's of length nz
!     ldz >= nz is the distance between sequences to be transformed
!     (ldz>nz is used on some architectures to reduce memory conflicts)
!     input  :  c_d(ldz*nsl)   (complex)
!     ### GPU VERION IN PLACE!!! #### output : cout_d(ldz*nsl) (complex - NOTA BENE: transform is not in-place!)
!     isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nz, nsl, ldz) are stored and re-used if available
#ifdef TRACK_FLOPS
     USE flops_tracker, ONLY : fft_ops
#endif
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: isign
     INTEGER, INTENT(IN) :: nsl, nz, ldz
     INTEGER(kind = cuda_stream_kind) :: stream
     LOGICAL, INTENT(IN), optional :: in_place

     COMPLEX (DP), DEVICE :: c_d(:), cout_d(:)

     REAL (DP)  :: tscale
     INTEGER    :: i, err, idir, ip, void, istat
#ifdef TRACK_FLOPS
     REAL (DP), SAVE :: zflops( ndims ) = 0.d0
#endif
     INTEGER, SAVE :: zdims( 3, ndims ) = -1
     INTEGER, SAVE :: icurrent = 1
     LOGICAL :: found
     LOGICAL :: is_inplace

     INTEGER :: tid

     ! Contains FFT factors ( PLAN )
     INTEGER, SAVE :: cufft_planz( ndims ) = 0

     IF( nsl < 0 ) THEN
       CALL fftx_error__(" fft_scalar: cft_1z ", " nsl out of range ", nsl)
     END IF

     !
     !   Here initialize table only if necessary
     !
     CALL lookup()

     IF( .NOT. found ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       CALL init_plan()

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !
     IF ( present( in_place ) ) THEN
       is_inplace = in_place
     ELSE
       is_inplace = .false.
     END IF
     !
     istat = cufftSetStream(cufft_planz(ip), stream)

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'GPU_cft_1z' )
#endif

     IF (isign < 0) THEN
        !print *,"exec cufft FWD",nz,ldz,nsl
        !call flush(6)
        istat = cufftExecZ2Z( cufft_planz( ip), c_d(1), c_d(1), CUFFT_FORWARD )
        tscale = 1.0_DP / nz
        IF (is_inplace) THEN
!$cuf kernel do(1) <<<*,*,0,stream>>>
           DO i=1, ldz * nsl
              c_d( i ) = c_d( i ) * tscale
           END DO
        ELSE
!$cuf kernel do(1) <<<*,*,0,stream>>>
           DO i=1, ldz * nsl
              cout_d( i ) = c_d( i ) * tscale
           END DO
        END IF
     ELSE IF (isign > 0) THEN
        !print *,"exec cufft INV",nz,ldz,nsl
        !call flush(6)
        IF (is_inplace) THEN
           istat = cufftExecZ2Z( cufft_planz( ip), c_d(1), c_d(1), CUFFT_INVERSE ) !CUFFT_FORWARD )
        ELSE
           istat = cufftExecZ2Z( cufft_planz( ip), c_d(1), cout_d(1), CUFFT_INVERSE ) !CUFFT_FORWARD )
        END IF
     END IF

#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'GPU_cft_1z' )
#endif

#ifdef TRACK_FLOPS
     fft_ops = fft_ops + zflops( ip )
#endif

     RETURN

     CONTAINS 

     SUBROUTINE lookup()
     DO ip = 1, ndims
        !   first check if there is already a table initialized
        !   for this combination of parameters
        found = ( nz == zdims(1,ip) ) .AND. ( nsl == zdims(2,ip) ) .AND. ( ldz == zdims(3,ip) )
        IF (found) EXIT
     END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()
       IMPLICIT NONE
       INTEGER, PARAMETER :: RANK=1
       INTEGER :: FFT_DIM(RANK), DATA_DIM(RANK)
       INTEGER :: STRIDE, DIST, BATCH

        FFT_DIM(1) = nz
       DATA_DIM(1) = ldz
            STRIDE = 1
              DIST = ldz
             BATCH = nsl

       IF( cufft_planz( icurrent) /= 0 ) istat = cufftDestroy( cufft_planz( icurrent) )

       istat = cufftPlanMany( cufft_planz( icurrent), RANK, FFT_DIM, & 
                              DATA_DIM, STRIDE, DIST, &
                              DATA_DIM, STRIDE, DIST, &
                              CUFFT_Z2Z, BATCH )

#if defined(__CUDA_DEBUG)
       print *,"INIT CUFFT Z PLAN: ",nz,"x",nsl,"x",ldz
#endif

#ifdef TRACK_FLOPS
       zflops( icurrent ) = 5.0d0 * REAL( nz ) * log( REAL( nz ) )/log( 2.d0 )
#endif

       zdims(1,icurrent) = nz; zdims(2,icurrent) = nsl; zdims(3,icurrent) = ldz;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

   END SUBROUTINE cft_1z_gpu

   SUBROUTINE cft_2xy_gpu(r_d, temp_d, nzl, nx, ny, ldx, ldy, isign, stream, pl2ix)

!     driver routine for nzl 2d complex fft's of lengths nx and ny
!     input : r_d(ldx*ldy)  complex, transform is in-place
!     ldx >= nx, ldy >= ny are the physical dimensions of the equivalent
!     2d array: r2d(ldx, ldy) (x first dimension, y second dimension)
!     (ldx>nx, ldy>ny used on some architectures to reduce memory conflicts)
!     pl2ix(nx) (optional) is 1 for columns along y to be transformed
!     isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nx,ny,nzl,ldx) are stored and re-used if available
!#ifdef TRACK_FLOPS
!     USE flops_tracker, ONLY : fft_ops
!#endif
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: isign, ldx, ldy, nx, ny, nzl
     INTEGER(kind = cuda_stream_kind) :: stream
     INTEGER, OPTIONAL, INTENT(IN) :: pl2ix(:)
!pgi$ ignore_tkr r_d, temp_d
     COMPLEX (DP), DEVICE :: r_d(ldx,ldy,nzl), temp_d(ldy,nzl,ldx)
     INTEGER :: i, k, j, err, idir, ip, kk, void, istat
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims( 6, ndims) = -1
     ! dims(5,:) = batch_1
     ! dims(6,:) = batch_2
     LOGICAL :: dofft( nfftx ), found
     INTEGER, PARAMETER  :: stdout = 6
#ifdef TRACK_FLOPS
     REAL (DP), SAVE :: xyflops( ndims ) = 0.d0
#endif

#if defined(__FFTW_ALL_XY_PLANES)
     INTEGER, SAVE :: cufft_plan_2d( ndims ) = 0
#else
     INTEGER, SAVE :: cufft_plan_x( ndims ) = 0
     INTEGER, SAVE :: cufft_plan_y( 2, ndims ) = 0
#endif
     INTEGER :: batch_1, batch_2
     !C_POINTER, SAVE :: fw_plan_2d( ndims ) = 0
     !C_POINTER, SAVE :: bw_plan_2d( ndims ) = 0

     dofft( 1 : nx ) = .TRUE.
     batch_1 = nx
     batch_2 = 0
     IF( PRESENT( pl2ix ) ) THEN
       IF( SIZE( pl2ix ) < nx ) &
         CALL fftx_error__( ' cft_2xy ', ' wrong dimension for arg no. 8 ', 1 )
       DO i = 1, nx
         IF( pl2ix(i) < 1 ) dofft( i ) = .FALSE.
       END DO

       i=1
       do while(pl2ix(i) >= 1 .and. i<=nx); i=i+1; END DO
       batch_1 = i-1
       do while(pl2ix(i) < 1 .and. i<=nx); i=i+1; END DO
       batch_2 = nx-i+1
#if 0
       !batch_2_start = i
       !do while(pl2ix(i) >= 1 .and. i<nx); i=i+1; END DO
      
       do while( i<=nx ) 

         do while(pl2ix(i) < 1 .and. i<=nx); i=i+1; END DO
         batch_start = i
         do while(pl2ix(i) >= 1 .and. i<=nx); i=i+1; END DO
         batch_end = i-1
         batch_count = batch_end - batch_start + 1
!         print *,"batch: ",batch_start,batch_end,batch_count
       enddo
#endif
     END IF

     !
     !   Here initialize table only if necessary
     !

     CALL lookup()

     IF( .NOT. found ) THEN

       !   no table exist for these parameters
       !   initialize a new one
       CALL init_plan()

     END IF

#if defined(__FFTW_ALL_XY_PLANES)
     istat = cufftSetStream(cufft_plan_2d(ip), stream)
#else
     istat = cufftSetStream(cufft_plan_x(ip), stream)
     istat = cufftSetStream(cufft_plan_y(1,ip), stream)
     istat = cufftSetStream(cufft_plan_y(2,ip), stream)
#endif

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'GPU_cft_2xy' )
#endif

     IF( isign < 0 ) THEN
        !
        tscale = 1.0_DP / ( nx * ny )
        !
#if defined(__FFTW_ALL_XY_PLANES)
        istat = cufftExecZ2Z( cufft_plan_2d(ip), r_d(1,1,1), r_d(1,1,1), CUFFT_FORWARD )
#else
        istat = cufftExecZ2Z( cufft_plan_x(ip), r_d(1,1,1), r_d(1,1,1), CUFFT_FORWARD )
        if(istat) print *,"error in fftxy fftx istat = ",istat

!$cuf kernel do(3) <<<*,(16,16,1), 0, stream>>>
        DO k=1, nzl
           DO i=1, ldx 
              DO j=1, ldy
                temp_d(j,k,i) = r_d(i,j,k)
              END DO
           END DO
        END DO
        
        
        if(batch_1>0) then
           istat = cufftExecZ2Z( cufft_plan_y(1,ip), temp_d(1,1,1), temp_d(1,1,1), CUFFT_FORWARD )
           if(istat) print *,"error in fftxy ffty batch_1 istat = ",istat
        end if

        if(batch_2>0) then     
           istat = cufftExecZ2Z( cufft_plan_y(2,ip), temp_d(1,1,nx-batch_2+1), temp_d(1,1,nx-batch_2+1), CUFFT_FORWARD )
           if(istat) print *,"error in fftxy ffty batch_2 istat = ",istat
        end if

!$cuf kernel do(3) <<<*,(16,16,1), 0, stream>>>
        DO k=1, nzl
           DO j=1, ldy
             DO i=1, ldx
                r_d(i,j,k) = temp_d(j,k,i) * tscale
              END DO
           END DO
        END DO
#endif

        !CALL ZDSCAL( ldx * ldy * nzl, tscale, r_d(1), 1)
        !
     ELSE IF( isign > 0 ) THEN
        !
        !print *,"exec cufft INV",nx,ny,ldx,ldy,nzl
#if defined(__FFTW_ALL_XY_PLANES)
        istat = cufftExecZ2Z( cufft_plan_2d(ip), r_d(1,1,1), r_d(1,1,1), CUFFT_INVERSE )
#else
!$cuf kernel do(3) <<<*,(16,16,1), 0, stream>>>
        DO k=1, nzl
           DO i=1, ldx
              DO j=1, ldy
                temp_d(j,k,i) = r_d(i,j,k)
              END DO
           END DO
        END DO

        if(batch_1>0) then
           istat = cufftExecZ2Z( cufft_plan_y(1,ip), temp_d(1,1,1), temp_d(1,1,1), CUFFT_INVERSE )
           if(istat) print *,"error in fftxy ffty batch_1 istat = ",istat
        end if

        if(batch_2>0) then
           istat = cufftExecZ2Z( cufft_plan_y(2,ip), temp_d(1,1,nx-batch_2+1), temp_d(1,1,nx-batch_2+1), CUFFT_INVERSE )
           if(istat) print *,"error in fftxy ffty batch_2 istat = ",istat
        end if

!$cuf kernel do(3) <<<*,(16,16,1), 0, stream>>>
        DO k=1, nzl
           DO j=1, ldy
             DO i=1, ldx
                r_d(i,j,k) = temp_d(j,k,i)
              END DO
           END DO
        END DO

!        do i = 1, nx
!           IF( dofft( i ) ) THEN
!             istat = cufftExecZ2Z( cufft_plan_y(ip), r_d(i), r_d(i), CUFFT_INVERSE )
!             if(istat) print *,"error in fftxy ffty istat = ",istat,i
!           END IF
!        end do

        istat = cufftExecZ2Z( cufft_plan_x(ip), r_d(1,1,1), r_d(1,1,1), CUFFT_INVERSE )
        if(istat) print *,"error in fftxy fftx istat = ",istat

#endif
        !
     END IF


#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'GPU_cft_2xy' )
#endif

#ifdef TRACK_FLOPS
     fft_ops = fft_ops + xyflops( ip )
#endif

     RETURN

   CONTAINS

     SUBROUTINE lookup()
     DO ip = 1, ndims
       !   first check if there is already a table initialized
       !   for this combination of parameters
       found = ( ny == dims(1,ip) ) .AND. ( nx == dims(3,ip) )
       found = found .AND. ( ldx == dims(2,ip) ) .AND.  ( nzl == dims(4,ip) )
       found = found .AND. ( batch_1 == dims(5,ip) ) .AND. (batch_2 == dims(6,ip) )
       IF (found) EXIT
     END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()
       IMPLICIT NONE
#if defined(__FFTW_ALL_XY_PLANES)
       INTEGER, PARAMETER :: RANK=2
       INTEGER :: FFT_DIM(RANK), DATA_DIM(RANK)
       INTEGER :: STRIDE, DIST, BATCH

        FFT_DIM(1) = ny
        FFT_DIM(2) = nx
       DATA_DIM(1) = ldy
       DATA_DIM(2) = ldx
            STRIDE = 1
              DIST = ldx*ldy
             BATCH = nzl

       IF( cufft_plan_2d( icurrent) /= 0 )  istat = cufftDestroy( cufft_plan_2d(icurrent) )

       istat = cufftPlanMany( cufft_plan_2d( icurrent), RANK, FFT_DIM, &
                              DATA_DIM, STRIDE, DIST, &
                              DATA_DIM, STRIDE, DIST, &
                              CUFFT_Z2Z, BATCH )

#if defined(__CUDA_DEBUG)
       print *,"INIT CUFFT ALL_XY PLAN: ",nx,"x",ny,"x",nzl,"ldx:",ldx,"batch:",batch_1,batch_2
#endif

#else
       INTEGER, PARAMETER :: RANK=1
       INTEGER :: FFT_DIM_X(RANK), DATA_DIM_X(RANK), FFT_DIM_Y(RANK), DATA_DIM_Y(RANK)
       INTEGER :: STRIDE_X, STRIDE_Y, DIST_X, DIST_Y, BATCH_X, BATCH_Y1, BATCH_Y2

        FFT_DIM_X(1) = nx
       DATA_DIM_X(1) = ldx
            STRIDE_X = 1
              DIST_X = ldx
             BATCH_X = ny*nzl

        FFT_DIM_Y(1) = ny
       DATA_DIM_Y(1) = ldy
            STRIDE_Y = 1
              DIST_Y = ldy
            BATCH_Y1 = nzl*BATCH_1
            BATCH_Y2 = nzl*BATCH_2


       IF( cufft_plan_x( icurrent) /= 0 )  istat = cufftDestroy( cufft_plan_x(icurrent) )
       IF( cufft_plan_y( 1, icurrent) /= 0 )  istat = cufftDestroy( cufft_plan_y(1,icurrent) )
       IF( cufft_plan_y( 2, icurrent) /= 0 )  istat = cufftDestroy( cufft_plan_y(2,icurrent) )

#if defined(__CUDA_DEBUG)
       print *,"INIT CUFFT XY PLAN: ",nx,"x",ny,"x",nzl,"ldx:",ldx,"batch:",batch_1,batch_2
#endif

       istat = cufftPlanMany( cufft_plan_x( icurrent), RANK, FFT_DIM_X, &
                              DATA_DIM_X, STRIDE_X, DIST_X, &
                              DATA_DIM_X, STRIDE_X, DIST_X, &
                              CUFFT_Z2Z, BATCH_X )

       istat = cufftPlanMany( cufft_plan_y( 1, icurrent), RANK, FFT_DIM_Y, &
                              DATA_DIM_Y, STRIDE_Y, DIST_Y, &
                              DATA_DIM_Y, STRIDE_Y, DIST_Y, &
                              CUFFT_Z2Z, BATCH_Y1 )

       istat = cufftPlanMany( cufft_plan_y( 2, icurrent), RANK, FFT_DIM_Y, &
                              DATA_DIM_Y, STRIDE_Y, DIST_Y, &
                              DATA_DIM_Y, STRIDE_Y, DIST_Y, &
                              CUFFT_Z2Z, BATCH_Y2 )


#endif

#ifdef TRACK_FLOPS
       xyflops( icurrent ) = REAL( ny*nzl )                    * 5.0d0 * REAL( nx ) * log( REAL( nx )  )/log( 2.d0 ) &
                           + REAL( nzl*BATCH_1 + nzl*BATCH_2 ) * 5.0d0 * REAL( ny ) * log( REAL( ny )  )/log( 2.d0 )

#endif

       dims(1,icurrent) = ny; dims(2,icurrent) = ldx;
       dims(3,icurrent) = nx; dims(4,icurrent) = nzl;
       dims(5,icurrent) = BATCH_1; dims(6,icurrent) = BATCH_2;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

   END SUBROUTINE cft_2xy_gpu
!
!=----------------------------------------------------------------------=!
!
!
!
!         3D scalar FFTs
!
!
!
!=----------------------------------------------------------------------=!
!


   SUBROUTINE cfft3d_gpu( f_d, nx, ny, nz, ldx, ldy, ldz, howmany, isign, stream )

  !     driver routine for 3d complex fft of lengths nx, ny, nz
  !     input  :  f_d(ldx*ldy*ldz)  complex, transform is in-place
  !     ldx >= nx, ldy >= ny, ldz >= nz are the physical dimensions
  !     of the equivalent 3d array: f3d(ldx,ldy,ldz)
  !     (ldx>nx, ldy>ny, ldz>nz may be used on some architectures
  !      to reduce memory conflicts - not implemented for FFTW)
  !     isign > 0 : f(G) => f(R)   ; isign < 0 : f(R) => f(G)
  !
  !     Up to "ndims" initializations (for different combinations of input
  !     parameters nx,ny,nz) are stored and re-used if available

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
     COMPLEX (DP), device :: f_d(:)
     INTEGER(kind = cuda_stream_kind) :: stream
     INTEGER :: i, k, j, err, idir, ip, istat
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims(4,ndims) = -1

!     C_POINTER, save :: fw_plan(ndims) = 0
!     C_POINTER, save :: bw_plan(ndims) = 0
     INTEGER, SAVE :: cufft_plan_3d( ndims ) = 0


     IF ( nx < 1 ) &
         call fftx_error__('cfft3d',' nx is less than 1 ', 1)
     IF ( ny < 1 ) &
         call fftx_error__('cfft3d',' ny is less than 1 ', 1)
     IF ( nz < 1 ) &
         call fftx_error__('cfft3',' nz is less than 1 ', 1)

     !
     !   Here initialize table only if necessary
     !
     CALL lookup()

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       CALL init_plan()

     END IF

     !
     !   Now perform the 3D FFT using the machine specific driver
     !
     istat = cufftSetStream(cufft_plan_3d(ip), stream)
     IF( isign < 0 ) THEN

        istat = cufftExecZ2Z( cufft_plan_3d(ip), f_d(1), f_d(1), CUFFT_FORWARD )

       tscale = 1.0_DP / DBLE( nx * ny * nz )
!$cuf kernel do(1) <<<*,*,0,stream>>>
        DO i=1, nx*ny*nz*howmany
           f_d( i ) = f_d( i ) * tscale
        END DO
!       call ZDSCAL( nx * ny * nz, tscale, f_d(1), 1)

     ELSE IF( isign > 0 ) THEN

!       call FFTW_INPLACE_DRV_3D( bw_plan(ip), 1, f_d(1), 1, 1 )
        istat = cufftExecZ2Z( cufft_plan_3d(ip), f_d(1), f_d(1), CUFFT_INVERSE )

     END IF

     RETURN

   CONTAINS 

     SUBROUTINE lookup()
     ip = -1
     DO i = 1, ndims
       !   first check if there is already a table initialized
       !   for this combination of parameters
       IF ( ( nx == dims(1,i) ) .and. &
            ( ny == dims(2,i) ) .and. &
            ( nz == dims(3,i) ) .and. &
            ( howmany == dims(4,i) ) ) THEN
         ip = i
         EXIT
       END IF
     END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()
       INTEGER, PARAMETER :: RANK=3
       INTEGER :: FFT_DIM(RANK), DATA_DIM(RANK)
       INTEGER :: STRIDE, DIST, BATCH

        FFT_DIM(1) = nz
        FFT_DIM(2) = ny
        FFT_DIM(3) = nx
       DATA_DIM(1) = ldz
       DATA_DIM(2) = ldy
       DATA_DIM(3) = ldx
            STRIDE = 1
              DIST = ldx*ldy*ldz
             BATCH = howmany

       IF( cufft_plan_3d( icurrent) /= 0 )  istat = cufftDestroy( cufft_plan_3d(icurrent) )

       istat = cufftPlanMany( cufft_plan_3d( icurrent), RANK, FFT_DIM, &
                              DATA_DIM, STRIDE, DIST, &
                              DATA_DIM, STRIDE, DIST, &
                              CUFFT_Z2Z, BATCH )

       !IF ( nx /= ldx .or. ny /= ldy .or. nz /= ldz ) &
       !  call fftx_error__('cfft3','not implemented',1)
       !IF( fw_plan(icurrent) /= 0 ) CALL DESTROY_PLAN_3D( fw_plan(icurrent) )
       !IF( bw_plan(icurrent) /= 0 ) CALL DESTROY_PLAN_3D( bw_plan(icurrent) )
       !idir = -1; CALL CREATE_PLAN_3D( fw_plan(icurrent), nx, ny, nz, idir)
       !idir =  1; CALL CREATE_PLAN_3D( bw_plan(icurrent), nx, ny, nz, idir)
       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
       dims(4,icurrent) = howmany
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

   END SUBROUTINE cfft3d_gpu
   
   SUBROUTINE cfft3ds_gpu (f_d, nx, ny, nz, ldx, ldy, ldz, howmany, isign, &
     do_fft_z, do_fft_y, stream)
     !
     !     driver routine for 3d complex "reduced" fft - see cfft3d
     !     The 3D fft are computed only on lines and planes which have
     !     non zero elements. These lines and planes are defined by
     !     the two integer vectors do_fft_y(nx) and do_fft_z(ldx*ldy)
     !     (1 = perform fft, 0 = do not perform fft)
     !     This routine is implemented only for fftw, essl, acml
     !     If not implemented, cfft3d is called instead
     !
     !----------------------------------------------------------------------
     !
     implicit none
         INTEGER, PARAMETER  :: stdout = 6
     
     integer :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
     !
     !   logical dimensions of the fft
     !   physical dimensions of the f_d array
     !   sign of the transformation
     
     complex(DP),device :: f_d ( ldx * ldy * ldz )
     integer :: do_fft_y(:), do_fft_z(:)
     integer(kind = cuda_stream_kind) :: stream
     !
     integer :: m, incx1, incx2
     INTEGER :: i, k, j, err, idir, ip,  ii, jj, istat
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims(3,ndims) = -1
     
     !TYPE(C_PTR), SAVE :: fw_plan ( 3, ndims ) = C_NULL_PTR
     !TYPE(C_PTR), SAVE :: bw_plan ( 3, ndims ) = C_NULL_PTR
     INTEGER, SAVE :: cufft_plan_1d( 3, ndims ) = 0
     !
     ! The current version of this function is massively outperformed by
     !  cfft3d_gpu. Leaving the call to full 3D FFT for the time being.
     ! 
     CALL cfft3d_gpu (f_d, nx, ny, nz, ldx, ldy, ldz, howmany, isign, stream)
     return
     
     tscale = 1.0_DP
     
     IF( ny /= ldy ) &
       CALL fftx_error__(' cfft3ds ', ' wrong dimensions: ny /= ldy ', 1 )
     IF( howmany < 1 ) &
       CALL fftx_error__(' cfft3ds ', ' howmany less than one ', 1 )
     IF( howmany /= 1 ) &
       CALL fftx_error__(' cfft3ds ', ' howmany different from 1, not yet implemented for cuFFT ', 1 )

     CALL lookup()
     
     IF( ip == -1 ) THEN
     
       !   no table exist for these parameters
       !   initialize a new one
     
       CALL init_plan()
     
     END IF
     
     istat = cufftSetStream(cufft_plan_1d(3, ip), stream)
     istat = cufftSetStream(cufft_plan_1d(2, ip), stream)
     istat = cufftSetStream(cufft_plan_1d(1, ip), stream)
     IF ( isign > 0 ) THEN
     
        !
        !  k-direction ...
        !
     
        !incx1 = ldx * ldy;  incx2 = 1;  m = 1
     
        do i =1, nx
           do j =1, ny
              ii = i + ldx * (j-1)
              if ( do_fft_z(ii) > 0) then
                 ! call dfftw_execute_dft( bw_plan( 3, ip), f_d( ii:), f_d( ii:) )
                 istat = cufftExecZ2Z( cufft_plan_1d(3, ip), f_d( ii:), f_d( ii:), CUFFT_INVERSE )
              end if
           end do
        end do
     
        !
        !  ... j-direction ...
        !
     
        !incx1 = ldx;  incx2 = ldx*ldy;  m = nz
     
        do i = 1, nx
           if ( do_fft_y( i ) == 1 ) then
             !call dfftw_execute_dft( bw_plan( 2, ip), f_d( i: ), f_d( i: ) )
             istat = cufftExecZ2Z( cufft_plan_1d(2, ip), f_d(i:), f_d( i:), CUFFT_INVERSE )
           endif
        enddo
     
        !
        !  ... i - direction
        !
     
        !incx1 = 1;  incx2 = ldx;  m = ldy*nz
     
        !call dfftw_execute_dft( bw_plan( 1, ip), f_d( 1: ), f_d( 1: ) )
        istat = cufftExecZ2Z( cufft_plan_1d(1, ip), f_d( 1:), f_d( 1:), CUFFT_INVERSE )
     
     ELSE
     
        !
        !  i - direction ...
        !
     
        !incx1 = 1;  incx2 = ldx;  m = ldy*nz
     
        !call dfftw_execute_dft( fw_plan( 1, ip), f_d( 1: ), f_d( 1: ) )
        istat = cufftExecZ2Z( cufft_plan_1d(1, ip), f_d( 1:), f_d( 1:), CUFFT_FORWARD )
     
        !
        !  ... j-direction ...
        !
     
        !incx1 = ldx;  incx2 = ldx*ldy;  m = nz
     
        do i = 1, nx
           if ( do_fft_y ( i ) == 1 ) then
             !call dfftw_execute_dft( fw_plan( 2, ip), f_d( i: ), f_d( i: ) )
             istat = cufftExecZ2Z( cufft_plan_1d(2, ip), f_d( i:), f_d( i:), CUFFT_FORWARD )
           endif
        enddo
     
        !
        !  ... k-direction
        !
     
        !incx1 = ldx * ny;  incx2 = 1;  m = 1
     
        do i = 1, nx
           do j = 1, ny
              ii = i + ldx * (j-1)
              if ( do_fft_z ( ii) > 0) then
                 !call dfftw_execute_dft( fw_plan( 3, ip), f_d(ii:), f_d(ii:) )
                 istat = cufftExecZ2Z( cufft_plan_1d(3, ip), f_d( ii:), f_d( ii:), CUFFT_FORWARD )
              end if
           end do
        end do
     
        !call DSCAL (2 * ldx * ldy * nz, 1.0_DP/(nx * ny * nz), f_d(1), 1)
        tscale = 1.0_DP / DBLE( nx * ny * nz )
        !$cuf kernel do(1) <<<*,*,0,stream>>>
        DO i=1, nx*ny*nz
           f_d( i ) = f_d( i ) * tscale
        END DO
     
     END IF
     RETURN
     
      CONTAINS
     
        SUBROUTINE lookup()
        ip = -1
        DO i = 1, ndims
          !   first check if there is already a table initialized
          !   for this combination of parameters
          IF( ( nx == dims(1,i) ) .and. ( ny == dims(2,i) ) .and. &
              ( nz == dims(3,i) ) ) THEN
            ip = i
            EXIT
          END IF
        END DO
        END SUBROUTINE lookup
     
        SUBROUTINE init_plan()
           INTEGER, PARAMETER :: RANK=3
           INTEGER :: FFT_DIM(RANK), HOW_MANY(RANK), DATA_DIM(RANK)
           INTEGER :: STRIDE(RANK), DIST(RANK)
           INTEGER :: i
           
            FFT_DIM(1) = nx
            FFT_DIM(2) = ny
            FFT_DIM(3) = nz
           HOW_MANY(1) = ny*nz
           HOW_MANY(2) = nz
           HOW_MANY(3) = 1
           DATA_DIM(1) = ldx
           DATA_DIM(2) = ldy
           DATA_DIM(3) = ldz
             STRIDE(1) = 1
             STRIDE(2) = ldx
             STRIDE(3) = ldx*ldy
               DIST(1) = ldx
               DIST(2) = ldx*ldy
               DIST(3) = 1
        
          IF( cufft_plan_1d( 1, icurrent) /= 0 )  istat = cufftDestroy( cufft_plan_1d(1, icurrent) )
          IF( cufft_plan_1d( 2, icurrent) /= 0 )  istat = cufftDestroy( cufft_plan_1d(2, icurrent) )
          IF( cufft_plan_1d( 3, icurrent) /= 0 )  istat = cufftDestroy( cufft_plan_1d(3, icurrent) )

       
       !CALL dfftw_plan_many_dft( fw_plan( 1, icurrent), &
       !     1, nx, ny*nz, f(1:), (/ldx, ldy, ldz/), 1, ldx, &
       !                   f(1:), (/ldx, ldy, ldz/), 1, ldx, idir, FFTW_ESTIMATE)
       !CALL dfftw_plan_many_dft( fw_plan( 2, icurrent), &
       !     1, ny, nz, f(1:), (/ldx, ldy, ldz/), ldx, ldx*ldy, &
       !                f(1:), (/ldx, ldy, ldz/), ldx, ldx*ldy, idir, FFTW_ESTIMATE)
       !CALL dfftw_plan_many_dft( fw_plan( 3, icurrent), &
       !     1, nz, 1, f(1:), (/ldx, ldy, ldz/), ldx*ldy, 1, &
       !               f(1:), (/ldx, ldy, ldz/), ldx*ldy, 1, idir, FFTW_ESTIMATE)
          
          
          DO i = 1,3 
            istat = cufftPlanMany( cufft_plan_1d( i, icurrent), 1, FFT_DIM(i), &
                              DATA_DIM(i), STRIDE(i), DIST(i), &
                              DATA_DIM(i), STRIDE(i), DIST(i), &
                              CUFFT_Z2Z, HOW_MANY(i) )
            ! set current dimension as a lookup parameter along the three directions
            dims(1,icurrent) = FFT_DIM(i)
          END DO
          ! dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
          
          ip = icurrent
          icurrent = MOD( icurrent, ndims ) + 1
        END SUBROUTINE init_plan
     
     END SUBROUTINE cfft3ds_gpu
#endif
!=----------------------------------------------------------------------=!
 END MODULE fft_scalar_cuFFT
!=----------------------------------------------------------------------=!


