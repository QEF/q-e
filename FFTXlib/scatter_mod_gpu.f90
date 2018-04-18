!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
! FFT base Module.
! Written by Carlo Cavazzoni, modified by Paolo Giannozzi
! Rewritten by Stefano de Gironcoli, ported to GPU by Pietro Bonfa'
!----------------------------------------------------------------------
!
!=----------------------------------------------------------------------=!
   MODULE scatter_mod_gpu
!=----------------------------------------------------------------------=!
#if defined(__CUDA)
        USE fft_types, ONLY: fft_type_descriptor
        USE fft_param

        IMPLICIT NONE


        SAVE

        PRIVATE

        PUBLIC :: fft_type_descriptor
        PUBLIC :: fft_scatter_xy_gpu, fft_scatter_yz_gpu, fft_scatter_tg_gpu, &
                  fft_scatter_tg_opt_gpu

!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!
!
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_xy_gpu ( desc, f_in_d, f_aux_d, nxx_, isgn )
  !-----------------------------------------------------------------------
  !
  ! transpose of the fft xy planes across the desc%comm2 communicator
  !
  ! a) From Y-oriented columns to X-oriented partial slices (isgn > 0)
  !    Active columns along the Y direction corresponding to a subset of the 
  !    active X values and a range of Z values (in this order) are stored 
  !    consecutively for each processor and are such that the subgroup owns
  !    all data for a range of Z values.
  !
  !    The Y pencil -> X-oriented partial slices transposition is performed 
  !    in the subgroup of processors (desc%comm2) owning this range of Z values.
  !
  !    The transpose takes place in two steps:
  !    1) on each processor the columns are sliced into sections along Y
  !       that are stored one after the other. On each processor, slices for
  !       processor "iproc2" are desc%nr2p(iproc2)*desc%nr1p(me2)*desc%my_nr3p big.
  !    2) all processors communicate to exchange slices (all sectin of columns with 
  !       Y in the slice belonging to "me" must be received, all the others 
  !       must be sent to "iproc2")
  !
  !    Finally one gets the "partial slice" representation: each processor has 
  !    all the X values of desc%my_nr2p Y and desc%my_nr3p Z values. 
  !    Data are organized with the X index running fastest, then Y, then Z.
  !
  !    f_in  contains the input Y columns, is destroyed on output
  !    f_aux contains the output X-oriented partial slices.
  !
  !  b) From planes to columns (isgn < 0)
  !
  !    Quite the same in the opposite direction
  !    f_aux contains the input X-oriented partial slices, is destroyed on output
  !    f_in  contains the output Y columns.
  !
  USE cudafor
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in) :: desc
  INTEGER, INTENT(in)                    :: nxx_, isgn
  COMPLEX (DP), DEVICE, INTENT(inout)    :: f_in_d (nxx_), f_aux_d (nxx_)

#if defined(__MPI)
  !
  COMPLEX (DP), ALLOCATABLE    :: f_in(:), f_aux(:)
  INTEGER :: ierr, me2, nproc2, iproc2, ncpx, my_nr2p, nr2px, ip, ip0
  INTEGER :: i, it, j, k, kfrom, kdest, offset, ioff, mc, m1, m3, i1, icompact, sendsize
  INTEGER, ALLOCATABLE :: ncp_(:)
  INTEGER, ALLOCATABLE, DEVICE :: nr1p__d(:), indx_d(:,:)
  !
#if defined(__NON_BLOCKING_SCATTER)
  INTEGER :: sh(desc%nproc2), rh(desc%nproc2)
#endif
  me2    = desc%mype2 + 1
  nproc2 = desc%nproc2 ; if ( abs(isgn) == 3 ) nproc2 = 1 

  ! allocate auxiliary array for columns distribution
  ALLOCATE ( ncp_(nproc2), nr1p__d(nproc2), indx_d(desc%nr1x,nproc2) )
  if ( abs (isgn) == 1 ) then          ! It's a potential FFT
     ncp_ = desc%nr1p * desc%my_nr3p
     nr1p__d= desc%nr1p
     indx_d = desc%indp
     my_nr2p=desc%my_nr2p
  else if ( abs (isgn) == 2 ) then     ! It's a wavefunction FFT
     ncp_ = desc%nr1w * desc%my_nr3p
     nr1p__d= desc%nr1w
     indx_d = desc%indw
     my_nr2p=desc%my_nr2p
  else if ( abs (isgn) == 3 ) then     ! It's a wavefunction FFT with task group
     ncp_ = desc%nr1w_tg * desc%my_nr3p! 
     nr1p__d= desc%nr1w_tg               !
     indx_d(1:desc%nr1x,1) = desc%indw_tg(1:desc%nr1x)          ! 
     my_nr2p=desc%nr2x                 ! in task group FFTs whole Y colums are distributed
  end if
  !
  CALL start_clock ('fft_scatt_xy')
  !
  ! calculate the message size
  !
  if ( abs (isgn) == 3 ) then
     nr2px = desc%nr2x ! if it's a task group FFT whole planes are distributed
  else
     nr2px = MAXVAL ( desc%nr2p )  ! maximum number of Y values to be disributed
  end if

  ncpx  = MAXVAL ( ncp_ )       ! maximum number of Y columns to be disributed
  
  sendsize = ncpx * nr2px       ! dimension of the scattered chunks (safe value)
  !
  ALLOCATE(f_in(nxx_), f_aux(nxx_))         ! allocate host copy of f_in and f_aux
  !
  ierr = 0
  IF (isgn.gt.0) THEN

     IF (nproc2==1) GO TO 10
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
     offset = 0
     DO iproc2 = 1, nproc2
        kdest = ( iproc2 - 1 ) * sendsize
        kfrom = offset
        !DO k = 1, ncp_(me2)
        !  !DO i = 1, desc%nr2p( iproc2 )
        !      !f_aux ( kdest + i ) =  f_in ( kfrom + i )
        !  !ENDDO
        !   
        !  !kdest = kdest + nr2px
        !  !kfrom = kfrom + desc%nr2x
        !ENDDO
        
        ierr = cudaMemcpy2D( f_aux(kdest + 1), nr2px, f_in_d(kfrom + 1 ), desc%nr2x, desc%nr2p( iproc2 ), ncp_(me2), cudaMemcpyDeviceToHost )
        offset = offset + desc%nr2p( iproc2 )
#if defined(__NON_BLOCKING_SCATTER)
        CALL mpi_isend( f_aux( (iproc2-1)*sendsize + 1 ), sendsize, &
                        MPI_DOUBLE_COMPLEX, iproc2-1, me2, desc%comm2, &
                        sh( iproc2 ), ierr )
#endif
     ENDDO
     !
     ! step two: communication  across the    nproc3    group
     !
#if defined(__NON_BLOCKING_SCATTER)
     DO iproc2 = 1, nproc2
        !
        ! now post the receive
        !
        CALL mpi_irecv( f_in( (iproc2-1)*sendsize + 1 ), sendsize, &
                        MPI_DOUBLE_COMPLEX, iproc2-1, MPI_ANY_TAG, &
                        desc%comm2, rh( iproc2 ), ierr )
        !IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', ' forward receive info<>0', abs(ierr) )
     ENDDO
     
     call mpi_waitall( nproc2, sh, MPI_STATUSES_IGNORE, ierr )
     !
#else
     CALL mpi_alltoall (f_aux(1), sendsize, MPI_DOUBLE_COMPLEX, f_in(1), &
                        sendsize, MPI_DOUBLE_COMPLEX, desc%comm2, ierr)
     !
     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
#endif
     f_in_d(1:nxx_) = f_in(1:nxx_)
     !
10   CONTINUE
     !
     !f_aux = (0.0_DP, 0.0_DP)
     !$cuf kernel do (1) <<<*,*>>>
     do i = lbound(f_aux_d,1), ubound(f_aux_d,1)
       f_aux_d(i) = (0.d0, 0.d0)
     end do
     !
#if defined(__NON_BLOCKING_SCATTER)
     if (nproc2 > 1) call mpi_waitall( nproc2, rh, MPI_STATUSES_IGNORE, ierr )
#endif
     !
     DO iproc2 = 1, nproc2
        
        !$cuf kernel do(2) <<<*,*>>>
        DO i = 1, ncp_( iproc2 )
           DO j = 1, my_nr2p
              it = ( iproc2 - 1 ) * sendsize + (i-1)*nr2px
              m3 = (i-1)/nr1p__d(iproc2)+1
              i1 = mod(i-1,nr1p__d(iproc2))+1
              m1 = indx_d(i1,iproc2)
              icompact = m1 + (m3-1)*desc%nr1x*my_nr2p + (j-1)*desc%nr1x
           ! old do loop started here
              !f_aux( m1 + (j-1)*desc%nr1x + (m3-1)*desc%nr1x*my_nr2p ) = f_in( j + it )
              f_aux_d( icompact ) = f_in_d( j + it )
              !icompact = icompact + desc%nr1x
           ENDDO
           !it = it + nr2px
        ENDDO
     ENDDO

  ELSE
     !
     !  "backward" scatter from planes to columns
     !
     DO iproc2 = 1, nproc2
        
        !$cuf kernel do (1) <<<*,*>>>
        DO i = 1, ncp_( iproc2 )
           DO j = 1, my_nr2p
              it = ( iproc2 - 1 ) * sendsize + (i-1)*nr2px
              m3 = (i-1)/nr1p__d(iproc2)+1 ; i1  = mod(i-1,nr1p__d(iproc2))+1 ;  m1 = indx_d(i1,iproc2)
              icompact = m1 + (m3-1)*desc%nr1x*my_nr2p + (j-1)*desc%nr1x
           !DO j = 1, my_nr2p
              !f_in( j + it ) = f_aux( m1 + (j-1)*desc%nr1x + (m3-1)*desc%nr1x*my_nr2p )
              f_in_d( j + it ) = f_aux_d( icompact )
           !   icompact = icompact + desc%nr1x
           ENDDO
           !it = it + nr2px
        ENDDO
        ! Async copy here?
        
#if defined(__NON_BLOCKING_SCATTER)
#error Not with GPUs now
        IF( nproc2 > 1 ) CALL mpi_isend( f_in( ( iproc2 - 1 ) * sendsize + 1 ), &
                              sendsize, MPI_DOUBLE_COMPLEX, iproc2-1, me2, &
                              desc%comm2, sh( iproc2 ), ierr )
#endif
     ENDDO
     IF (nproc2==1) GO TO 20
     !
     !  step two: communication
     !
     f_in(1:nxx_) = f_in_d(1:nxx_)
#if ! defined(__NON_BLOCKING_SCATTER)
     CALL mpi_alltoall (f_in(1), sendsize, MPI_DOUBLE_COMPLEX, f_aux(1), &
                        sendsize, MPI_DOUBLE_COMPLEX, desc%comm2, ierr)

     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
#else
     DO iproc2 = 1, nproc2
        CALL mpi_irecv( f_aux( (iproc2-1)*sendsize + 1 ), sendsize, &
                        MPI_DOUBLE_COMPLEX, iproc2-1, MPI_ANY_TAG, &
                        desc%comm2, rh(iproc2), ierr )
        ! IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', ' backward receive info<>0', abs(ierr) )
     ENDDO
        
     call mpi_waitall( nproc2, sh, MPI_STATUSES_IGNORE, ierr )
     call mpi_waitall( nproc2, rh, MPI_STATUSES_IGNORE, ierr )
#endif
     !
     !  step one: store contiguously the columns
     !
     ! ensures that no garbage is present in the output 
     ! not useless ... clean the array to be returned from the garbage of previous A2A step
     !f_in = (0.0_DP, 0.0_DP) !
     !$cuf kernel do (1) <<<*,*>>>
     do i = lbound(f_in_d,1), ubound(f_in_d,1)
       f_in_d(i) = (0.d0, 0.d0)
     end do
     !
     offset = 0
     DO iproc2 = 1, nproc2
        kdest = ( iproc2 - 1 ) * sendsize
        kfrom = offset 
        !DO k = 1, ncp_(me2)
        !   DO i = 1, desc%nr2p( iproc2 )
        !      f_in ( kfrom + i ) = f_aux ( kdest + i )
        !   ENDDO
        !   kdest = kdest + nr2px
        !   kfrom = kfrom + desc%nr2x
        !ENDDO
        ierr = cudaMemcpy2D( f_in_d(kfrom +1 ), desc%nr2x, f_aux(kdest + 1), nr2px, desc%nr2p( iproc2 ), ncp_(me2), cudaMemcpyHostToDevice )
        offset = offset + desc%nr2p( iproc2 )
     ENDDO

20   CONTINUE


  ENDIF

  DEALLOCATE ( ncp_ , nr1p__d, indx_d , f_in, f_aux)
  CALL stop_clock ('fft_scatt_xy')

#endif

  RETURN
99 format ( 20 ('(',2f12.9,')') )

END SUBROUTINE fft_scatter_xy_gpu
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_yz_gpu ( desc, f_in_d, f_aux_d, nxx_, isgn )
  !-----------------------------------------------------------------------
  !
  ! transpose of the fft yz planes across the desc%comm3 communicator
  !
  ! a) From Z-oriented columns to Y-oriented colums (isgn > 0)
  !    Active columns (or sticks or pencils) along the Z direction for each 
  !    processor are stored consecutively and are such that they correspond 
  !    to a subset of the active X values.
  !
  !    The pencil -> slices transposition is performed in the subgroup 
  !    of processors (desc%comm3) owning these X values.
  !
  !    The transpose takes place in two steps:
  !    1) on each processor the columns are sliced into sections along Z
  !       that are stored one after the other. On each processor, slices for
  !       processor "iproc3" are desc%nr3p(iproc3)*desc%nsw/nsp(me) big.
  !    2) all processors communicate to exchange slices (all columns with 
  !       Z in the slice belonging to "me" must be received, all the others 
  !       must be sent to "iproc3")
  !
  !    Finally one gets the "slice" representation: each processor has 
  !    desc%nr3p(mype3) Z values of all the active pencils along Y for the
  !    X values of the current group. Data are organized with the Y index
  !    running fastest, then the reordered X values, then Z.
  !
  !    f_in  contains the input Z columns, is destroyed on output
  !    f_aux contains the output Y colums.
  !
  !  b) From planes to columns (isgn < 0)
  !
  !    Quite the same in the opposite direction
  !    f_aux contains the input Y columns, is destroyed on output
  !    f_in  contains the output Z columns.
  !
  USE cudafor
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in) :: desc
  INTEGER, INTENT(in)                    :: nxx_, isgn
  COMPLEX (DP), DEVICE, INTENT(inout)    :: f_in_d (nxx_), f_aux_d (nxx_)

#if defined(__MPI)
  COMPLEX (DP),    ALLOCATABLE :: f_in(:), f_aux(:)
  INTEGER, DEVICE, POINTER, CONTIGUOUS :: desc_ismap_d(:)
  !
  INTEGER :: ierr, me, me2, me2_start, me2_end, me3, nproc3, iproc3, ncpx, nr3px, ip, ip0
  INTEGER :: i, it, it0, k, kfrom, kdest, offset, ioff, mc, m1, m2, i1,  sendsize
  INTEGER, ALLOCATABLE :: ncp_(:)
  INTEGER, DEVICE, POINTER, CONTIGUOUS :: ir1p__d(:)
  INTEGER :: my_nr1p_
  !
#if defined(__NON_BLOCKING_SCATTER)
  INTEGER :: sh(desc%nproc3), rh(desc%nproc3)
#endif
  !
  me     = desc%mype  + 1
  me2    = desc%mype2 + 1
  me3    = desc%mype3 + 1
  nproc3 = desc%nproc3

  ! allocate auxiliary array for columns distribution
  ALLOCATE ( ncp_( desc%nproc) )

  me2_start = me2 ; me2_end = me2
  if ( abs (isgn) == 1 ) then      ! It's a potential FFT
     ncp_ = desc%nsp 
     my_nr1p_ = count (desc%ir1p   > 0)
     ir1p__d => desc%ir1p_d
  else if ( abs (isgn) == 2 ) then ! It's a wavefunction FFT
     ncp_ = desc%nsw
     my_nr1p_ = count (desc%ir1w   > 0)
     ir1p__d => desc%ir1w_d
  else if ( abs (isgn) == 3 ) then ! It's a wavefunction FFT with task group
     ncp_ = desc%nsw
     my_nr1p_ = count (desc%ir1w_tg > 0)
     ir1p__d => desc%ir1w_tg_d
     me2_start = 1 ; me2_end = desc%nproc2
  end if
  !
  CALL start_clock ('fft_scatt_yz')
  !
  ! calculate the message size
  !
  nr3px = MAXVAL ( desc%nr3p )  ! maximum number of Z values to be exchanged
  ncpx  = MAXVAL ( ncp_ )       ! maximum number of Z columns to be exchanged
  if (abs(isgn)==3) then
     ncpx = ncpx * desc%nproc2  ! if it's a task group FFT groups of columns are exchanged
  end if
  
  sendsize = ncpx * nr3px       ! dimension of the scattered chunks
  
  ALLOCATE(f_in(nxx_), f_aux(nxx_))         ! allocate host copy of f_in and f_aux
  desc_ismap_d => desc%ismap_d

  ierr = 0
  IF (isgn.gt.0) THEN

     IF (nproc3==1) GO TO 10
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
     offset = 0
     DO iproc3 = 1, nproc3
        kdest = ( iproc3 - 1 ) * sendsize
        kfrom = offset
        DO me2 = me2_start, me2_end
           ip = desc%iproc(me2,me3)
        !   DO k = 1, ncp_ (ip)   ! was ncp_(me3)
        !      DO i = 1, desc%nr3p( iproc3 )
        !         f_aux ( kdest + i ) =  f_in ( kfrom + i )
        !      ENDDO
        !      kdest = kdest + nr3px
        !      kfrom = kfrom + desc%nr3x
        !   ENDDO

        !      integer function cudaMemcpy2D(dst, dpitch, src, spitch, width, height, kdir)
           ierr = cudaMemcpy2D( f_aux(kdest + 1), nr3px, f_in_d(kfrom + 1 ), desc%nr3x, desc%nr3p( iproc3 ), ncp_(ip), cudaMemcpyDeviceToHost )
           kdest = kdest + nr3px*ncp_ (ip)
           kfrom = kfrom + desc%nr3x*ncp_ (ip)
        ENDDO
        offset = offset + desc%nr3p( iproc3 )
#if defined(__NON_BLOCKING_SCATTER)
        CALL mpi_isend( f_aux( ( iproc3 - 1 ) * sendsize + 1 ), sendsize, &
                        MPI_DOUBLE_COMPLEX, iproc3-1, me3, desc%comm3, &
                        sh( iproc3 ), ierr )
#endif
     ENDDO
     !
     ! ensures that no garbage is present in the output 
     ! useless; the later accessed elements are overwritten by the A2A step
     !
     ! step two: communication  across the    nproc3    group
     !
#if defined(__NON_BLOCKING_SCATTER)
     DO iproc3 = 1, nproc3
        !
        ! now post the receive
        !
        CALL mpi_irecv( f_in( (iproc3-1)*sendsize + 1 ), sendsize, &
                        MPI_DOUBLE_COMPLEX, iproc3-1, MPI_ANY_TAG, &
                        desc%comm3, rh( iproc3 ), ierr )
        !IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', ' forward receive info<>0', abs(ierr) )
        !
        !
     ENDDO
     
     call mpi_waitall( nproc3, sh, MPI_STATUSES_IGNORE, ierr )
     !
#else
     CALL mpi_alltoall (f_aux(1), sendsize, MPI_DOUBLE_COMPLEX, f_in(1), &
                        sendsize, MPI_DOUBLE_COMPLEX, desc%comm3, ierr)

     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
#endif     
     !
     f_in_d(1:nxx_) = f_in(1:nxx_)
     !
10   CONTINUE
     !
     !f_aux = (0.0_DP, 0.0_DP) !
     !$cuf kernel do (1) <<<*,*>>>
     do i = lbound(f_aux_d,1), ubound(f_aux_d,1)
       f_aux_d(i) = (0.d0, 0.d0)
     end do
     
#if defined(__NON_BLOCKING_SCATTER)
     if (nproc3 > 1) call mpi_waitall( nproc3, rh, MPI_STATUSES_IGNORE, ierr )
#endif
     !
     DO iproc3 = 1, desc%nproc3
        it0 = ( iproc3 - 1 ) * sendsize
        DO me2 = me2_start, me2_end
           ip = desc%iproc( me2, iproc3)
           ioff = desc%iss(ip)
!$cuf kernel do(2) <<<*,*>>>
           DO i = 1, ncp_( ip ) ! was ncp_(iproc3)
              DO k = 1, desc%my_nr3p
                 it = it0 + (i-1)*nr3px
                 mc = desc_ismap_d( i + ioff ) ! this is  m1+(m2-1)*nr1x  of the  current pencil
                 m1 = mod (mc-1,desc%nr1x) + 1 ; m2 = (mc-1)/desc%nr1x + 1 
                 i1 = m2 + ( ir1p__d(m1) - 1 ) * desc%nr2x + (k-1)*desc%nr2x*my_nr1p_
              
                 f_aux_d( i1 ) = f_in_d( k + it )
                 !i1 = i1 + desc%nr2x*my_nr1p_
              ENDDO
           ENDDO
           it0 = it0 + ncp_( ip )*nr3px
        ENDDO
     ENDDO

  ELSE
     !
     !  "backward" scatter from planes to columns
     !
     DO iproc3 = 1, desc%nproc3
        it0 = ( iproc3 - 1 ) * sendsize
        DO me2 = me2_start, me2_end
           ip = desc%iproc(me2, iproc3)
           ioff = desc%iss(ip)
!$cuf kernel do(2) <<<*,*>>>
           DO i = 1, ncp_( ip )
              DO k = 1, desc%my_nr3p
                 it = it0 + (i-1)*nr3px
                 mc = desc_ismap_d( i + ioff ) ! this is  m1+(m2-1)*nr1x  of the  current pencil
                 m1 = mod (mc-1,desc%nr1x) + 1 ; m2 = (mc-1)/desc%nr1x + 1 
                 i1 = m2 + ( ir1p__d(m1) - 1 ) * desc%nr2x + (k-1)*(desc%nr2x * my_nr1p_)
              
                 f_in_d( k + it ) = f_aux_d( i1 )
                 !i1 = i1 + desc%nr2x * my_nr1p_
              ENDDO
              !it = it + nr3px
           ENDDO
           it0 = it0 + ncp_( ip )*nr3px
        ENDDO
#if defined(__NON_BLOCKING_SCATTER)
#error Non blocking with GPU not yet available
        IF( nproc3 > 1 ) CALL mpi_isend( f_in( ( iproc3 - 1 ) * sendsize + 1 ), &
                                        sendsize, MPI_DOUBLE_COMPLEX, iproc3-1, &
                                        me3, desc%comm3, sh( iproc3 ), ierr )
#endif
     ENDDO
     
     IF( nproc3 == 1 ) GO TO 20
     !
     f_in(1:nxx_) = f_in_d(1:nxx_)
     !
     !  step two: communication
     !
#if ! defined(__NON_BLOCKING_SCATTER)
     CALL mpi_alltoall (f_in(1), sendsize, MPI_DOUBLE_COMPLEX, f_aux(1), &
                        sendsize, MPI_DOUBLE_COMPLEX, desc%comm3, ierr)

     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
#else
     DO iproc3 = 1, desc%nproc3
        CALL mpi_irecv( f_aux( (iproc3-1)*sendsize + 1 ), sendsize, &
                        MPI_DOUBLE_COMPLEX, iproc3-1, MPI_ANY_TAG, &
                        desc%comm3, rh(iproc3), ierr )
        ! IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', ' backward receive info<>0', abs(ierr) )
     ENDDO
        
     call mpi_waitall( desc%nproc3, sh, MPI_STATUSES_IGNORE, ierr )
     call mpi_waitall( desc%nproc3, rh, MPI_STATUSES_IGNORE, ierr )
#endif
     !
     !  step one: store contiguously the columns
     !
     offset = 0
     DO iproc3 = 1, nproc3
        kdest = ( iproc3 - 1 ) * sendsize
        kfrom = offset
        DO me2 = me2_start, me2_end
           ip = desc%iproc(me2,me3)
        !   DO k = 1, ncp_ (ip) 
        !      DO i = 1, desc%nr3p( iproc3 )
        !         f_in ( kfrom + i ) = f_aux ( kdest + i )
        !      ENDDO
        !      kdest = kdest + nr3px
        !      kfrom = kfrom + desc%nr3x
        !   ENDDO
          ierr = cudaMemcpy2D( f_in_d(kfrom +1 ), desc%nr3x, f_aux(kdest + 1), nr3px, desc%nr3p( iproc3 ), ncp_ (ip), cudaMemcpyHostToDevice )
          kdest = kdest + nr3px*ncp_ (ip)
          kfrom = kfrom + desc%nr3x*ncp_ (ip)
        ENDDO
        offset = offset + desc%nr3p( iproc3 )
     ENDDO

     ! clean extra array elements in each stick

     IF( desc%nr3x /= desc%nr3 ) THEN
!$cuf kernel do(2) <<<*,*>>>
       DO k = 1, ncp_ ( desc%mype+1 ) 
          DO i = desc%nr3, desc%nr3x
             f_in_d( (k-1)*desc%nr3x + i ) = 0.0d0 
          END DO
       END DO
     END IF


20   CONTINUE

  ENDIF

  DEALLOCATE ( ncp_ , f_aux, f_in)
  CALL stop_clock ('fft_scatt_yz')

#endif

  RETURN
98 format ( 10 ('(',2f12.9,')') )
99 format ( 20 ('(',2f12.9,')') )

END SUBROUTINE fft_scatter_yz_gpu
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_tg_gpu ( desc, f_in_d, f_aux_d, nxx_, isgn )
  !-----------------------------------------------------------------------
  !
  ! task group wavefunction redistribution
  !
  ! a) (isgn >0 ) From many-wfc partial-plane arrangement to single-wfc whole-plane one
  !
  ! b) (isgn <0 ) From single-wfc whole-plane arrangement to many-wfc partial-plane one
  !
  ! in both cases:
  !    f_in  contains the input data, is overwritten with the desired output
  !    f_aux is used as working array, may contain garbage in output
  !
  USE cudafor
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in) :: desc
  INTEGER, INTENT(in)                    :: nxx_, isgn
  COMPLEX (DP), DEVICE, INTENT(inout)    :: f_in_d (nxx_), f_aux_d (nxx_)
  COMPLEX (DP), ALLOCATABLE              :: f_in(:), f_aux(:)

  INTEGER :: ierr

  CALL start_clock ('fft_scatt_tg')

  if ( abs (isgn) /= 3 ) call fftx_error__ ('fft_scatter_tg', 'wrong call', 1 )
  ALLOCATE(f_in(nxx_), f_aux(nxx_))
#if defined(__MPI)
  !
  f_aux(1:nxx_) = f_in_d(1:nxx_)
  if ( isgn > 0 ) then

     CALL MPI_ALLTOALLV( f_aux,  desc%tg_snd, desc%tg_sdsp, MPI_DOUBLE_COMPLEX, &
                         f_in, desc%tg_rcv, desc%tg_rdsp, MPI_DOUBLE_COMPLEX, desc%comm2, ierr)
     IF( ierr /= 0 ) CALL fftx_error__( 'fft_scatter_tg', ' alltoall error 1 ', abs(ierr) )

  else

     CALL MPI_ALLTOALLV( f_aux,  desc%tg_rcv, desc%tg_rdsp, MPI_DOUBLE_COMPLEX, &
                         f_in, desc%tg_snd, desc%tg_sdsp, MPI_DOUBLE_COMPLEX, desc%comm2, ierr)
     IF( ierr /= 0 ) CALL fftx_error__( 'fft_scatter_tg', ' alltoall error 2 ', abs(ierr) )
  end if
  f_in_d(1:nxx_) = f_in(1:nxx_)
#endif
  DEALLOCATE(f_in, f_aux)
  CALL stop_clock ('fft_scatt_tg')
  RETURN
99 format ( 20 ('(',2f12.9,')') )

END SUBROUTINE fft_scatter_tg_gpu
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_tg_opt_gpu ( desc, f_in_d, f_out_d, nxx_, isgn )
  !-----------------------------------------------------------------------
  !
  ! task group wavefunction redistribution
  !
  ! a) (isgn >0 ) From many-wfc partial-plane arrangement to single-wfc whole-plane one
  !
  ! b) (isgn <0 ) From single-wfc whole-plane arrangement to many-wfc partial-plane one
  !
  ! in both cases:
  !    f_in  contains the input data
  !    f_out contains the output data
  !
  USE cudafor
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in) :: desc
  INTEGER, INTENT(in)                    :: nxx_, isgn
  COMPLEX (DP), DEVICE, INTENT(inout)    :: f_in_d (nxx_), f_out_d (nxx_)
  COMPLEX (DP), ALLOCATABLE              :: f_in(:), f_out(:)

  INTEGER :: ierr

  CALL start_clock ('fft_scatt_tg')

  if ( abs (isgn) /= 3 ) call fftx_error__ ('fft_scatter_tg', 'wrong call', 1 )
  !
  ALLOCATE(f_in(nxx_), f_out(nxx_))
  f_in(1:nxx_) = f_in_d(1:nxx_)
#if defined(__MPI)
  !
  if ( isgn > 0 ) then

     CALL MPI_ALLTOALLV( f_in,  desc%tg_snd, desc%tg_sdsp, MPI_DOUBLE_COMPLEX, &
                         f_out, desc%tg_rcv, desc%tg_rdsp, MPI_DOUBLE_COMPLEX, desc%comm2, ierr)
     IF( ierr /= 0 ) CALL fftx_error__( 'fft_scatter_tg', ' alltoall error 1 ', abs(ierr) )

  else

     CALL MPI_ALLTOALLV( f_in,  desc%tg_rcv, desc%tg_rdsp, MPI_DOUBLE_COMPLEX, &
                         f_out, desc%tg_snd, desc%tg_sdsp, MPI_DOUBLE_COMPLEX, desc%comm2, ierr)
     IF( ierr /= 0 ) CALL fftx_error__( 'fft_scatter_tg', ' alltoall error 2 ', abs(ierr) )
  end if

#endif
  f_out_d(1:nxx_) = f_out(1:nxx_)
  DEALLOCATE(f_in, f_out)
  CALL stop_clock ('fft_scatt_tg')

  RETURN
99 format ( 20 ('(',2f12.9,')') )

END SUBROUTINE fft_scatter_tg_opt_gpu

#endif
!=----------------------------------------------------------------------=!
END MODULE scatter_mod_gpu
!=----------------------------------------------------------------------=!
