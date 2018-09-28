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
!----------------------------------------------------------------------
!
#ifdef __CUDA
!=----------------------------------------------------------------------=!
   MODULE scatter_mod_2d_gpu
!=----------------------------------------------------------------------=!

        USE fft_types, ONLY: fft_type_descriptor
        USE fft_param


        USE cudafor

        IMPLICIT NONE

        SAVE

        PRIVATE

        PUBLIC :: fft_scatter_gpu, fft_scatter_gpu_batch, fft_scatter_batch_a_gpu, fft_scatter_batch_b_gpu

!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!

!----------------------------------------------------------------------------
SUBROUTINE fft_scatter_gpu ( dfft, f_in_d, f_in, nr3x, nxx_, f_aux_d, f_aux, ncp_, npp_, isgn )
  !
  USE cudafor
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), TARGET, INTENT(in) :: dfft
  INTEGER, INTENT(in)           :: nr3x, nxx_, isgn, ncp_ (:), npp_ (:)
  COMPLEX (DP), DEVICE, INTENT(inout)   :: f_in_d (nxx_), f_aux_d (nxx_)
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)
  INTEGER :: cuf_i, cuf_j, nswip
  INTEGER :: istat
  INTEGER, POINTER, DEVICE :: p_ismap_d(:)
#if defined(__MPI)

  INTEGER :: srh(2*dfft%nproc)
  INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
  INTEGER :: me_p, nppx, mc, j, npp, nnp, ii, it, ip, ioff, sendsiz, ncpx, ipp, nblk, nsiz
  !
  INTEGER, ALLOCATABLE, DIMENSION(:) :: offset_proc
  INTEGER :: iter, dest, sorc
  INTEGER :: istatus(MPI_STATUS_SIZE)


  p_ismap_d => dfft%ismap_d
  !
  me     = dfft%mype + 1
  !
  nprocp = dfft%nproc
  !
  istat = cudaDeviceSynchronize()
  !
  ncpx = maxval(ncp_)
  nppx = maxval(npp_)

  ! This should never happend and should go away
  IF ( dfft%nproc == 1 ) THEN
     nppx = dfft%nr3x
  END IF
  sendsiz = ncpx * nppx
  !

  ierr = 0
  IF (isgn.gt.0) THEN

     IF (nprocp==1) GO TO 10
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
     offset = 0
     !f_aux = (0.d0, 0.d0)

     DO gproc = 1, nprocp
        !
        kdest = ( gproc - 1 ) * sendsiz
        kfrom = offset
        !
#ifdef __GPU_MPI

!$cuf kernel do(2) <<<*,*>>>
        DO k = 1, ncp_ (me)
           DO i = 1, npp_ ( gproc )
             f_aux_d( kdest + i + (k-1)*nppx ) = f_in_d( kfrom + i + (k-1)*nr3x )
           END DO
        END DO

#else
        istat = cudaMemcpy2D( f_aux(kdest + 1), nppx, f_in_d(kfrom + 1 ), nr3x, npp_(gproc), ncp_(me), cudaMemcpyDeviceToHost )
        if( istat ) CALL fftx_error__("fft_scatter", "ERROR cudaMemcpy2D failed : ", istat)
#endif

        offset = offset + npp_ ( gproc )
     ENDDO
     !
     ! maybe useless; ensures that no garbage is present in the output
     !
     !! f_in = 0.0_DP
     !
     ! step two: communication
     !
     gcomm = dfft%comm


     CALL start_clock ('a2a_fw')
#ifdef __GPU_MPI

     istat = cudaDeviceSynchronize()
     DO iter = 2, nprocp
        IF(IAND(nprocp, nprocp-1) == 0) THEN
          sorc = IEOR( me-1, iter-1 )
        ELSE
          sorc = MOD(me-1 - (iter-1) + nprocp, nprocp)
        ENDIF

        call MPI_IRECV( f_in_d((sorc)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, srh(iter-1), ierr )

     ENDDO

     DO iter = 2, nprocp
        IF(IAND(nprocp, nprocp-1) == 0) THEN
          dest = IEOR( me-1, iter-1 )
        ELSE
          dest = MOD(me-1 + (iter-1), nprocp)
        ENDIF

        call MPI_ISEND( f_aux_d((dest)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, dest, 0, gcomm, srh(iter+nprocp-2), ierr )

     ENDDO

     istat = cudaMemcpyAsync( f_in_d( (me-1)*sendsiz + 1), f_aux_d((me-1)*sendsiz + 1), sendsiz, stream=dfft%a2a_comp )

     call MPI_WAITALL(2*nprocp-2, srh, MPI_STATUSES_IGNORE, ierr)
     istat = cudaDeviceSynchronize()
#else
     CALL mpi_alltoall (f_aux(1), sendsiz, MPI_DOUBLE_COMPLEX, f_in(1), sendsiz, MPI_DOUBLE_COMPLEX, gcomm, ierr)
#endif
     CALL stop_clock ('a2a_fw')

     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )

#ifndef __GPU_MPI
     f_in_d(1:sendsiz*dfft%nproc) = f_in(1:sendsiz*dfft%nproc)
#endif

     !
10   CONTINUE
     
     !f_aux_d = (0.d0, 0.d0)
     !$cuf kernel do (1) <<<*,*>>>
     do i = lbound(f_aux_d,1), ubound(f_aux_d,1)
       f_aux_d(i) = (0.d0, 0.d0)
     end do

     IF( isgn == 1 ) THEN

        npp = dfft%nr3p( me )
        nnp = dfft%nnp

        DO ip = 1, dfft%nproc
           ioff = dfft%iss( ip )
           nswip = dfft%nsp( ip )
!$cuf kernel do(2) <<<*,*>>>
           DO cuf_j = 1, npp
              DO cuf_i = 1, nswip
                 it = ( ip - 1 ) * sendsiz + (cuf_i-1)*nppx
                 mc = p_ismap_d( cuf_i + ioff )
                 f_aux_d( mc + ( cuf_j - 1 ) * nnp ) = f_in_d( cuf_j + it )
              ENDDO
           ENDDO
        ENDDO
     ELSE

        npp  = dfft%nr3p( me )
        nnp  = dfft%nnp
        !
        ip = 1
        !
        DO gproc = 1, dfft%nproc
           !
           ioff = dfft%iss( ip )
           nswip =  dfft%nsw( ip )
           !
!$cuf kernel do(2) <<<*,*>>>
           DO cuf_j = 1, npp
              DO cuf_i = 1, nswip
                 !
                 mc = p_ismap_d( cuf_i + ioff )
                 !
                 it = (cuf_i-1) * nppx + ( gproc - 1 ) * sendsiz
                 !
                 f_aux_d( mc + ( cuf_j - 1 ) * nnp ) = f_in_d( cuf_j + it )
                 !
              ENDDO
              !
           ENDDO
           !
           ip = ip + 1
           !
        ENDDO
     END IF
  ELSE
     !
     !  "backward" scatter from planes to columns
     !
     IF( isgn == -1 ) THEN

        npp = dfft%nr3p( me )
        nnp = dfft%nnp

        DO ip = 1, dfft%nproc
           ioff = dfft%iss( ip )
           nswip = dfft%nsp( ip )
!$cuf kernel do(2) <<<*,*>>>
           DO cuf_j = 1, npp
              DO cuf_i = 1, nswip
                 mc = p_ismap_d( cuf_i + ioff )
                 it = ( ip - 1 ) * sendsiz + (cuf_i-1)*nppx
                 f_in_d( cuf_j + it ) = f_aux_d( mc + ( cuf_j - 1 ) * nnp )
              ENDDO
           ENDDO

        ENDDO
     ELSE

        npp  = dfft%nr3p( me )
        nnp  = dfft%nnp
        !
        ip = 1
        !
        DO gproc = 1, dfft%nproc
           !
           ioff = dfft%iss( ip )
           !
           nswip = dfft%nsw( ip )
!$cuf kernel do(2) <<<*,*>>>
           DO cuf_j = 1, npp
              DO cuf_i = 1, nswip
                 !
                 mc = p_ismap_d( cuf_i + ioff )
                 !
                 it = (cuf_i-1) * nppx + ( gproc - 1 ) * sendsiz
                 !
                 f_in_d( cuf_j + it ) = f_aux_d( mc + ( cuf_j - 1 ) * nnp )
              ENDDO
              !
           ENDDO
           !
           ip = ip + 1
           !
        ENDDO
     END IF
     !
     IF( nprocp == 1 ) GO TO 20
     !
     !  step two: communication
     !
     gcomm = dfft%comm
     !
#ifndef __GPU_MPI
     f_in(1:sendsiz*dfft%nproc) = f_in_d(1:sendsiz*dfft%nproc)
#endif
     !
     ! CALL mpi_barrier (gcomm, ierr)  ! why barrier? for buggy openmpi over ib
     CALL start_clock ('a2a_bw')
#ifdef __GPU_MPI

     istat = cudaDeviceSynchronize()
     DO iter = 2, nprocp
        IF(IAND(nprocp, nprocp-1) == 0) THEN
          sorc = IEOR( me-1, iter-1 )
        ELSE
          sorc = MOD(me-1 - (iter-1) + nprocp, nprocp)
        ENDIF

        call MPI_IRECV( f_aux_d((sorc)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, srh(iter-1), ierr )

     ENDDO

     DO iter = 2, nprocp
        IF(IAND(nprocp, nprocp-1) == 0) THEN
          dest = IEOR( me-1, iter-1 )
        ELSE
          dest = MOD(me-1 + (iter-1), nprocp)
        ENDIF

        call MPI_ISEND( f_in_d((dest)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, dest, 0, gcomm, srh(iter+nprocp-2), ierr )

     ENDDO

     istat = cudaMemcpyAsync( f_aux_d( (me-1)*sendsiz + 1), f_in_d((me-1)*sendsiz + 1), sendsiz, stream=dfft%a2a_comp )

     call MPI_WAITALL(2*nprocp-2, srh, MPI_STATUSES_IGNORE, ierr)
     istat = cudaDeviceSynchronize()
#else
     CALL mpi_alltoall (f_in(1), sendsiz, MPI_DOUBLE_COMPLEX, f_aux(1), sendsiz, MPI_DOUBLE_COMPLEX, gcomm, ierr)
#endif
     CALL stop_clock ('a2a_bw')
     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
     !
     !  step one: store contiguously the columns
     !
     !! f_in = 0.0_DP
     !
     offset = 0

     DO gproc = 1, nprocp
        !
        kdest = ( gproc - 1 ) * sendsiz
        kfrom = offset
        !
#ifdef __GPU_MPI

!$cuf kernel do(2) <<<*,*>>>
        DO k = 1, ncp_ (me)
           DO i = 1, npp_ ( gproc )
             f_in_d( kfrom + i + (k-1)*nr3x ) = f_aux_d( kdest + i + (k-1)*nppx )
           END DO
        END DO

#else
        istat = cudaMemcpy2D( f_in_d(kfrom +1 ), nr3x, f_aux(kdest + 1), nppx, npp_(gproc), ncp_(me), cudaMemcpyHostToDevice )
#endif
        offset = offset + npp_ ( gproc )
     ENDDO

20   CONTINUE

  ENDIF

  istat = cudaDeviceSynchronize()

#endif

  RETURN

END SUBROUTINE fft_scatter_gpu

SUBROUTINE fft_scatter_gpu_batch ( dfft, f_in_d, f_in, nr3x, nxx_, f_aux_d, f_aux, ncp_, npp_, isgn, batchsize, srh )
  !
  USE cudafor
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), TARGET, INTENT(in) :: dfft
  INTEGER, INTENT(in)           :: nr3x, nxx_, isgn, ncp_ (:), npp_ (:)
  COMPLEX (DP), DEVICE, INTENT(inout)   :: f_in_d (batchsize * nxx_), f_aux_d (batchsize * nxx_)
  COMPLEX (DP), INTENT(inout)   :: f_in (batchsize * nxx_), f_aux (batchsize * nxx_)
  INTEGER, INTENT(IN) :: batchsize
  INTEGER, INTENT(INOUT) :: srh(2*dfft%nproc)
  INTEGER :: cuf_i, cuf_j, nswip
  INTEGER :: istat
  INTEGER, POINTER, DEVICE :: p_ismap_d(:)
#if defined(__MPI)

  INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
  INTEGER :: me_p, nppx, mc, j, npp, nnp, nnr, ii, it, ip, ioff, sendsiz, ncpx, ipp, nblk, nsiz
  !

  INTEGER, ALLOCATABLE, DIMENSION(:) :: offset_proc, kdest_proc, kfrom_proc
  INTEGER :: iter, dest, sorc
  INTEGER :: istatus(MPI_STATUS_SIZE)


  p_ismap_d => dfft%ismap_d

  me     = dfft%mype + 1
  !
  nprocp = dfft%nproc
  !
  istat = cudaDeviceSynchronize()
  !
  !ALLOCATE( offset_proc( nprocp ), kdest_proc( nprocp ), kfrom_proc( nprocp ) )
  ncpx = maxval(ncp_)
  nppx = maxval(npp_)
  !
  IF ( dfft%nproc == 1 ) THEN
     nppx = dfft%nr3x
  END IF
  !
  sendsiz = batchsize * ncpx * nppx
  nnr     = dfft%nnr
  !
  ierr = 0
  IF (isgn.gt.0) THEN

     IF (nprocp==1) GO TO 10
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
     offset = 0

     DO gproc = 1, nprocp
        !
        kdest = ( gproc - 1 ) * sendsiz
        kfrom = offset
        !
#ifdef __GPU_MPI
        !
!$cuf kernel do(2) <<<*,*>>>
        DO k = 1, batchsize * ncpx
           DO i = 1, npp_ ( gproc )
             f_aux_d( kdest + i + (k-1)*nppx ) = f_in_d( kfrom + i + (k-1)*nr3x )
           END DO
        END DO
        !
#else
        istat = cudaMemcpy2D( f_aux(kdest + 1), nppx, f_in_d(kfrom + 1 ), nr3x, npp_(gproc), batchsize * ncpx, cudaMemcpyDeviceToHost )
        if( istat ) CALL fftx_error__ ("ERROR cudaMemcpy2D failed : ",istat)
#endif
        !
        offset = offset + npp_ ( gproc )
     ENDDO
     !
     ! maybe useless; ensures that no garbage is present in the output
     !
     !! f_in = 0.0_DP
     !
     ! step two: communication
     !
     gcomm = dfft%comm

     CALL start_clock ('a2a_fw')

     istat = cudaDeviceSynchronize()
     DO iter = 2, nprocp
        IF(IAND(nprocp, nprocp-1) == 0) THEN
          sorc = IEOR( me-1, iter-1 )
        ELSE
          sorc = MOD(me-1 - (iter-1) + nprocp, nprocp)
        ENDIF

#ifdef __GPU_MPI
        call MPI_IRECV( f_in_d((sorc)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, srh(iter-1), ierr )
#else
        call MPI_IRECV( f_in((sorc)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, srh(iter-1), ierr )
#endif

     ENDDO

     DO iter = 2, nprocp
        IF(IAND(nprocp, nprocp-1) == 0) THEN
          dest = IEOR( me-1, iter-1 )
        ELSE
          dest = MOD(me-1 + (iter-1), nprocp)
        ENDIF

#ifdef __GPU_MPI
        call MPI_ISEND( f_aux_d((dest)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, dest, 0, gcomm, srh(iter+nprocp-2), ierr )
#else
        call MPI_ISEND( f_aux((dest)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, dest, 0, gcomm, srh(iter+nprocp-2), ierr )
#endif

     ENDDO

#ifdef __GPU_MPI
     istat = cudaMemcpyAsync( f_in_d( (me-1)*sendsiz + 1), f_aux_d((me-1)*sendsiz + 1), sendsiz, stream=dfft%a2a_comp )
#else
     f_in((me-1)*sendsiz + 1 : me*sendsiz) = f_aux((me-1)*sendsiz + 1 : me*sendsiz)
#endif

     call MPI_WAITALL(2*nprocp-2, srh, MPI_STATUSES_IGNORE, ierr)
     istat = cudaDeviceSynchronize()

     CALL stop_clock ('a2a_fw')

     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )

#ifndef __GPU_MPI
     f_in_d(1:sendsiz*dfft%nproc) = f_in(1:sendsiz*dfft%nproc)
#endif

     !
10   CONTINUE
     
     ! Zero out f_aux_d
     !$cuf kernel do (1) <<<*,*>>>
     do i = lbound(f_aux_d,1), ubound(f_aux_d,1)
       f_aux_d(i) = (0.d0, 0.d0)
     end do

     IF( isgn == 1 ) THEN

        npp = dfft%nr3p( me )
        nnp = dfft%nnp

        DO ip = 1, dfft%nproc
           ioff = dfft%iss( ip )
           nswip = dfft%nsp( ip )
!$cuf kernel do(2) <<<*,*>>>
           DO cuf_j = 1, npp
              DO cuf_i = 1, nswip
                 it = ( ip - 1 ) * sendsiz + (cuf_i-1)*nppx
                 mc = p_ismap_d( cuf_i + ioff )
                 f_aux_d( mc + ( cuf_j - 1 ) * nnp ) = f_in_d( cuf_j + it )
              ENDDO
           ENDDO
        ENDDO
     ELSE

        npp  = dfft%nr3p( me )
        nnp  = dfft%nnp
        !
        ip = 1
        !
        DO gproc = 1,  dfft%nproc
           !
           ioff = dfft%iss( ip )
           nswip =  dfft%nsw( ip )
           !
!$cuf kernel do(3) <<<*,*>>>
           DO i = 0, batchsize-1
              DO cuf_j = 1, npp
                 DO cuf_i = 1, nswip
                    !
                    mc = p_ismap_d( cuf_i + ioff )
                    !
                    it = (cuf_i-1) * nppx + ( gproc - 1 ) * sendsiz + i*nppx*ncpx
                    !
                    f_aux_d( mc + ( cuf_j - 1 ) * nnp + i*nnr ) = f_in_d( cuf_j + it )
                 ENDDO
                 !
              ENDDO
           ENDDO
           !
           ip = ip + 1
           !
           !
        ENDDO
     END IF
  ELSE
     !
     !  "backward" scatter from planes to columns
     !
     IF( isgn == -1 ) THEN

        npp = dfft%nr3p( me )
        nnp = dfft%nnp

        DO ip = 1, dfft%nproc
           ioff = dfft%iss( ip )
           nswip = dfft%nsp( ip )
!$cuf kernel do(2) <<<*,*>>>
           DO cuf_j = 1, npp
              DO cuf_i = 1, nswip
                 mc = p_ismap_d( cuf_i + ioff )
                 it = ( ip - 1 ) * sendsiz + (cuf_i-1)*nppx
                    f_in_d( cuf_j + it ) = f_aux_d( mc + ( cuf_j - 1 ) * nnp )
              ENDDO
           ENDDO

        ENDDO
     ELSE

        npp  = dfft%nr3p( me )
        nnp  = dfft%nnp
        !
        DO ip = 1, dfft%nproc
           !
           !
           ioff = dfft%iss( ip )
           !
           nswip = dfft%nsw( ip )
!$cuf kernel do(3) <<<*,*>>>
           DO i = 0, batchsize-1
              DO cuf_j = 1, npp
                 DO cuf_i = 1, nswip
                    !
                    mc = p_ismap_d( cuf_i + ioff )
                    !
                    it = (cuf_i-1) * nppx + ( ip - 1 ) * sendsiz + i*nppx*ncpx
                    !
                    f_in_d( cuf_j + it ) = f_aux_d( mc + ( cuf_j - 1 ) * nnp + i*nnr )
                    !
                 ENDDO
                 !
              ENDDO
           ENDDO
           !
        ENDDO
     END IF
     !
     IF( nprocp == 1 ) GO TO 20
     !
     !  step two: communication
     !
     gcomm = dfft%comm

#ifndef __GPU_MPI
     f_in(1:sendsiz*dfft%nproc) = f_in_d(1:sendsiz*dfft%nproc)
#endif

     ! CALL mpi_barrier (gcomm, ierr)  ! why barrier? for buggy openmpi over ib
     CALL start_clock ('a2a_bw')


     istat = cudaDeviceSynchronize()
     DO iter = 2, nprocp
        IF(IAND(nprocp, nprocp-1) == 0) THEN
          sorc = IEOR( me-1, iter-1 )
        ELSE
          sorc = MOD(me-1 - (iter-1) + nprocp, nprocp)
        ENDIF

#ifdef __GPU_MPI
        call MPI_IRECV( f_aux_d((sorc)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, srh(iter-1), ierr )
#else
        call MPI_IRECV( f_aux((sorc)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, srh(iter-1), ierr )
#endif

     ENDDO

     DO iter = 2, nprocp
        IF(IAND(nprocp, nprocp-1) == 0) THEN
          dest = IEOR( me-1, iter-1 )
        ELSE
          dest = MOD(me-1 + (iter-1), nprocp)
        ENDIF

#ifdef __GPU_MPI
        call MPI_ISEND( f_in_d((dest)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, dest, 0, gcomm, srh(iter+nprocp-2), ierr )
#else
        call MPI_ISEND( f_in((dest)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, dest, 0, gcomm, srh(iter+nprocp-2), ierr )
#endif

     ENDDO

#ifdef __GPU_MPI
     istat = cudaMemcpyAsync( f_aux_d( (me-1)*sendsiz + 1), f_in_d((me-1)*sendsiz + 1), sendsiz, stream=dfft%a2a_comp )
#else
     f_aux( (me-1)*sendsiz + 1:me*sendsiz) = f_in((me-1)*sendsiz + 1:me*sendsiz)
#endif

     call MPI_WAITALL(2*nprocp-2, srh, MPI_STATUSES_IGNORE, ierr)
     istat = cudaDeviceSynchronize()

     CALL stop_clock ('a2a_bw')
     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
     !
     !  step one: store contiguously the columns
     !
     !! f_in = 0.0_DP
     !
     offset = 0

     DO gproc = 1, nprocp
        !
        kdest = ( gproc - 1 ) * sendsiz
        kfrom = offset
        !
#ifdef __GPU_MPI

!$cuf kernel do(2) <<<*,*>>>
        DO k = 1, batchsize * ncpx
           DO i = 1, npp_ ( gproc )
             f_in_d( kfrom + i + (k-1)*nr3x ) = f_aux_d( kdest + i + (k-1)*nppx )
           END DO
        END DO

#else
        istat = cudaMemcpy2D( f_in_d(kfrom +1 ), nr3x, f_aux(kdest + 1), nppx, npp_(gproc), batchsize * ncpx, cudaMemcpyHostToDevice )
#endif
        offset = offset + npp_ ( gproc )
     ENDDO

20   CONTINUE

  ENDIF

  istat = cudaDeviceSynchronize()

#endif

  RETURN

END SUBROUTINE fft_scatter_gpu_batch

SUBROUTINE fft_scatter_batch_a_gpu ( dfft, f_in_d, f_in, nr3x, nxx_, f_aux_d, f_aux, f_aux2_d, f_aux2, ncp_, npp_, isgn, batchsize, batch_id )
  !
  USE cudafor
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), TARGET, INTENT(in) :: dfft
  INTEGER, INTENT(in)           :: nr3x, nxx_, isgn, ncp_ (:), npp_ (:)
  COMPLEX (DP), DEVICE, INTENT(inout)   :: f_in_d (batchsize * nxx_), f_aux_d (batchsize * nxx_), f_aux2_d(batchsize * nxx_)
  COMPLEX (DP), INTENT(inout)   :: f_in (batchsize * nxx_), f_aux (batchsize * nxx_), f_aux2(batchsize * nxx_)
  INTEGER, INTENT(IN) :: batchsize, batch_id
  INTEGER :: cuf_i, cuf_j, nswip
  INTEGER :: istat
  INTEGER, POINTER, DEVICE :: p_ismap_d(:)
#if defined(__MPI)
  INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
  INTEGER :: me_p, nppx, mc, j, npp, nnp, nnr, ii, it, ip, ioff, sendsiz, ncpx, ipp, nblk, nsiz
  !
  LOGICAL :: use_tg
  INTEGER, ALLOCATABLE, DIMENSION(:) :: offset_proc !, kdest_proc, kfrom_proc_
  INTEGER :: iter, dest, sorc
  INTEGER :: istatus(MPI_STATUS_SIZE)

  p_ismap_d => dfft%ismap_d
  me     = dfft%mype + 1
  !
  nprocp = dfft%nproc
  !
  !istat = cudaDeviceSynchronize()
  !

#ifdef __IPC
#ifndef __GPU_MPI
  call get_ipc_peers( dfft%IPC_PEER )
#endif
#endif
  !
  ncpx = maxval(ncp_)
  nppx = maxval(npp_)
  !
  IF ( dfft%nproc == 1 ) THEN
     nppx = dfft%nr3x
  END IF
  !
  sendsiz = batchsize * ncpx * nppx
  nnr     = dfft%nnr
  !
  ierr = 0

  IF (isgn.gt.0) THEN

     IF (nprocp==1) GO TO 10
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
     ALLOCATE( offset_proc( nprocp ) )
     offset = 0
     DO proc = 1, nprocp
        offset_proc( proc ) = offset
        offset = offset + npp_ ( proc )
     ENDDO

     DO iter = 2, nprocp
        IF(IAND(nprocp, nprocp-1) == 0) THEN
          dest = IEOR( me-1, iter-1 )
        ELSE
          dest = MOD(me-1 + (iter-1), nprocp)
        ENDIF
        proc = dest + 1

        !
        kdest = ( proc - 1 ) * sendsiz
        kfrom = offset_proc( proc )
        !
#ifdef __GPU_MPI
        istat = cudaMemcpy2DAsync( f_aux_d(kdest + 1), nppx, f_in_d(kfrom + 1 ), nr3x, npp_(proc), batchsize * ncpx,cudaMemcpyDeviceToDevice, dfft%bstreams(batch_id) )
        if( istat ) CALL fftx_error__ ("ERROR cudaMemcpy2D failed : ",istat)

#else
#ifdef __IPC
        IF(dfft%IPC_PEER( dest + 1 ) .eq. 1) THEN
           istat = cudaMemcpy2DAsync( f_aux_d(kdest + 1), nppx, f_in_d(kfrom + 1 ), nr3x, npp_(proc), batchsize * ncpx,cudaMemcpyDeviceToDevice, dfft%bstreams(batch_id) )
           if( istat ) CALL fftx_error__ ("ERROR cudaMemcpy2D failed : ",istat)
        ELSE
           istat = cudaMemcpy2DAsync( f_aux(kdest + 1), nppx, f_in_d(kfrom + 1 ), nr3x, npp_(proc), batchsize * ncpx,cudaMemcpyDeviceToHost, dfft%bstreams(batch_id) )
           if( istat ) CALL fftx_error__ ("ERROR cudaMemcpy2D failed : ",istat)
        ENDIF
#else
        istat = cudaMemcpy2DAsync( f_aux(kdest + 1), nppx, f_in_d(kfrom + 1 ), nr3x, npp_(proc), batchsize * ncpx,cudaMemcpyDeviceToHost, dfft%bstreams(batch_id) )
        if( istat ) CALL fftx_error__ ("ERROR cudaMemcpy2D failed : ",istat)
#endif
#endif
     ENDDO

     istat = cudaEventRecord( dfft%bevents(batch_id), dfft%bstreams(batch_id) )
     DEALLOCATE( offset_proc )
     !
10   CONTINUE
     
  ELSE
     !
     !  "backward" scatter from planes to columns
     !
     IF( isgn == -1 ) THEN

        npp = dfft%nr3p( me )
        nnp = dfft%nnp

        DO iter = 1, nprocp
           IF(IAND(nprocp, nprocp-1) == 0) THEN
              dest = IEOR( me-1, iter-1 )
           ELSE
              dest = MOD(me-1 + (iter-1), nprocp)
           ENDIF

           ip = dest + 1
           ioff = dfft%iss( ip )
           nswip = dfft%nsp( ip )
!$cuf kernel do(2) <<<*,*,0,dfft%a2a_comp>>>
           DO cuf_j = 1, npp
              DO cuf_i = 1, nswip
                 mc = p_ismap_d( cuf_i + ioff )
                 it = ( ip - 1 ) * sendsiz + (cuf_i-1)*nppx
                 f_aux2_d( cuf_j + it ) = f_aux_d( mc + ( cuf_j - 1 ) * nnp )
              ENDDO
           ENDDO

        ENDDO

     ELSE

        npp  = dfft%nr3p( me )
        nnp  = dfft%nnp
        !
        DO iter = 1, nprocp
           IF(IAND(nprocp, nprocp-1) == 0) THEN
              dest = IEOR( me-1, iter-1 )
           ELSE
              dest = MOD(me-1 + (iter-1), nprocp)
           ENDIF
           gproc = dest + 1
           !
           ioff = dfft%iss( gproc )
           !
           nswip = dfft%nsw( gproc )
!$cuf kernel do(3) <<<*,*, 0, dfft%a2a_comp>>>
           DO i = 0, batchsize-1
              DO cuf_j = 1, npp
                 DO cuf_i = 1, nswip
                   !
                   mc = p_ismap_d( cuf_i + ioff )
                   !
                   it = (cuf_i-1) * nppx + ( gproc - 1 ) * sendsiz + i*nppx*ncpx
                   !
                   f_aux2_d( cuf_j + it ) = f_aux_d( mc + ( cuf_j - 1 ) * nnp + i*nnr )
                 ENDDO
                 !
              ENDDO
           ENDDO
           !
        ENDDO
     END IF

#ifndef __GPU_MPI
     i = cudaEventRecord(dfft%bevents(batch_id), dfft%a2a_comp)
     i = cudaStreamWaitEvent(dfft%bstreams(batch_id), dfft%bevents(batch_id), 0)

     DO proc = 1, nprocp
        IF (proc .ne. me) THEN
#ifdef __IPC
           IF(dfft%IPC_PEER( proc ) .eq. 0) THEN
              kdest = ( proc - 1 ) * sendsiz
              istat = cudaMemcpyAsync( f_aux2(kdest+1), f_aux2_d(kdest+1), sendsiz, stream=dfft%bstreams(batch_id) )
           ENDIF
#else
           kdest = ( proc - 1 ) * sendsiz
           istat = cudaMemcpyAsync( f_aux2(kdest+1), f_aux2_d(kdest+1), sendsiz, stream=dfft%bstreams(batch_id) )
#endif
        ENDIF
     ENDDO
#endif


#ifdef __GPU_MPI
     istat = cudaEventRecord( dfft%bevents(batch_id), dfft%a2a_comp )
#else
     istat = cudaEventRecord( dfft%bevents(batch_id), dfft%bstreams(batch_id) )
#endif
  ENDIF

  !istat = cudaDeviceSynchronize()

#endif

  RETURN

END SUBROUTINE fft_scatter_batch_a_gpu


SUBROUTINE fft_scatter_batch_b_gpu ( dfft, f_in_d, f_in, nr3x, nxx_, f_aux_d, f_aux, f_aux2_d, f_aux2, ncp_, npp_, isgn, batchsize, batch_id )
  !
  USE cudafor
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), TARGET, INTENT(in) :: dfft
  INTEGER, INTENT(in)           :: nr3x, nxx_, isgn, ncp_ (:), npp_ (:)
  COMPLEX (DP), DEVICE, INTENT(inout)   :: f_in_d (batchsize * nxx_), f_aux_d (batchsize * nxx_), f_aux2_d (batchsize * nxx_)
  COMPLEX (DP), INTENT(inout)   :: f_in (batchsize * nxx_), f_aux (batchsize * nxx_), f_aux2(batchsize * nxx_)
  INTEGER, INTENT(IN) :: batchsize, batch_id
  INTEGER :: cuf_i, cuf_j, nswip
  INTEGER :: istat
  INTEGER, POINTER, DEVICE :: p_ismap_d(:)
#if defined(__MPI)

  INTEGER :: k, offset, proc, ierr, me, nprocp, gproc, gcomm, i, kdest, kfrom
  INTEGER :: me_p, nppx, mc, j, npp, nnp, nnr, ii, it, ip, ioff, sendsiz, ncpx, ipp, nblk, nsiz

  INTEGER :: iter, dest, sorc, req_cnt
  INTEGER :: istatus(MPI_STATUS_SIZE)

  p_ismap_d => dfft%ismap_d

  me     = dfft%mype + 1
  !
  nprocp = dfft%nproc
  !
  !istat = cudaDeviceSynchronize()
  !
  ncpx = maxval(ncp_)
  nppx = maxval(npp_)
  !
  IF ( dfft%nproc == 1 ) THEN
     nppx = dfft%nr3x
  END IF
  !
  sendsiz = batchsize * ncpx * nppx
  nnr     = dfft%nnr
#ifdef __IPC
  call get_ipc_peers( dfft%IPC_PEER )
#endif
  !

  ierr = 0
  IF (isgn.gt.0) THEN

     IF (nprocp==1) GO TO 10
     ! step two: communication
     !
     gcomm = dfft%comm

     ! JR Note: Holding off staging receives until buffer is packed.
     istat = cudaEventSynchronize( dfft%bevents(batch_id) ) 
     CALL start_clock ('A2A')
#ifdef __IPC
     !TODO: possibly remove this barrier by ensuring recv buffer is not used by previous operation
     call MPI_Barrier( gcomm, ierr )
#endif
     req_cnt = 0

     DO iter = 2, nprocp
        IF(IAND(nprocp, nprocp-1) == 0) THEN
          sorc = IEOR( me-1, iter-1 )
        ELSE
          sorc = MOD(me-1 - (iter-1) + nprocp, nprocp)
        ENDIF
#ifdef __IPC
        IF(dfft%IPC_PEER( sorc + 1 ) .eq. 0) THEN
#endif
#ifdef __GPU_MPI
           call MPI_IRECV( f_aux2_d((sorc)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, dfft%srh(req_cnt+1, batch_id), ierr )
#else
           call MPI_IRECV( f_aux2((sorc)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, dfft%srh(req_cnt+1, batch_id), ierr )
#endif
           req_cnt = req_cnt + 1
#ifdef __IPC
        ENDIF
#endif

     ENDDO

     DO iter = 2, nprocp
        IF(IAND(nprocp, nprocp-1) == 0) THEN
          dest = IEOR( me-1, iter-1 )
        ELSE
          dest = MOD(me-1 + (iter-1), nprocp)
        ENDIF
#ifdef __IPC
        IF(dfft%IPC_PEER( dest + 1 ) .eq. 1) THEN
           call ipc_send( f_aux_d((dest)*sendsiz + 1), sendsiz, f_aux2_d((me-1)*sendsiz + 1), 1, dest, gcomm, ierr )
        ELSE
#endif
#ifdef __GPU_MPI
           call MPI_ISEND( f_aux_d((dest)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, dest, 0, gcomm, dfft%srh(req_cnt+1, batch_id), ierr )
#else
           call MPI_ISEND( f_aux((dest)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, dest, 0, gcomm, dfft%srh(req_cnt+1, batch_id), ierr )
#endif
           req_cnt = req_cnt + 1
#ifdef __IPC
        ENDIF
#endif
     ENDDO

     offset = 0
     DO proc = 1, me-1
        offset = offset + npp_ ( proc )
     ENDDO
     istat = cudaMemcpy2DAsync( f_aux2_d((me-1)*sendsiz + 1), nppx, f_in_d(offset + 1 ), nr3x, npp_(me), batchsize * ncpx,cudaMemcpyDeviceToDevice, dfft%bstreams(batch_id) )

     if(req_cnt .gt. 0) then
        call MPI_WAITALL(req_cnt, dfft%srh(1:req_cnt, batch_id), MPI_STATUSES_IGNORE, ierr)
     endif

#ifdef __IPC
     call sync_ipc_sends( gcomm )
     call MPI_Barrier( gcomm, ierr )
#endif
     CALL stop_clock ('A2A')

     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )

#ifndef __GPU_MPI
     DO proc = 1, nprocp
        IF (proc .ne. me) THEN
#ifdef __IPC
           IF(dfft%IPC_PEER( proc ) .eq. 0) THEN
              kdest = ( proc - 1 ) * sendsiz
              istat = cudaMemcpyAsync( f_aux2_d(kdest+1), f_aux2(kdest+1), sendsiz, stream=dfft%bstreams(batch_id) )
           ENDIF
#else
           kdest = ( proc - 1 ) * sendsiz
           istat = cudaMemcpyAsync( f_aux2_d(kdest+1), f_aux2(kdest+1), sendsiz, stream=dfft%bstreams(batch_id) )
#endif
        ENDIF
     ENDDO
#endif

    i = cudaEventRecord(dfft%bevents(batch_id), dfft%bstreams(batch_id))
    i = cudaStreamWaitEvent(dfft%a2a_comp, dfft%bevents(batch_id), 0)


     !
10   CONTINUE
     
     ! Zero out f_aux_d
     !$cuf kernel do (1) <<<*,*,0,dfft%a2a_comp>>>
     do i = lbound(f_aux_d,1), ubound(f_aux_d,1)
       f_aux_d(i) = (0.d0, 0.d0)
     end do

     IF( isgn == 1 ) THEN

        npp = dfft%nr3p( me )
        nnp = dfft%nnp

        DO ip = 1, nprocp
           ioff = dfft%iss( ip )
           nswip = dfft%nsp( ip )
!$cuf kernel do(2) <<<*,*,0,dfft%a2a_comp>>>
           DO cuf_j = 1, npp
              DO cuf_i = 1, nswip
                 it = ( ip - 1 ) * sendsiz + (cuf_i-1)*nppx
                 mc = p_ismap_d( cuf_i + ioff )
                 f_aux_d( mc + ( cuf_j - 1 ) * nnp ) = f_aux2_d( cuf_j + it )
              ENDDO
           ENDDO
        ENDDO
     ELSE
        !
        npp  = dfft%nr3p( me )
        nnp  = dfft%nnp
        !
        ip = 1
        !
        DO gproc = 1, nprocp
           !
           ioff = dfft%iss( ip )
           nswip =  dfft%nsw( ip )
           !
!$cuf kernel do(3) <<<*,*,0,dfft%a2a_comp>>>
           DO i = 0, batchsize-1
              DO cuf_j = 1, npp
                DO cuf_i = 1, nswip
                   !
                   mc = p_ismap_d( cuf_i + ioff )
                   !
                   it = (cuf_i-1) * nppx + ( gproc - 1 ) * sendsiz + i*nppx*ncpx
                   !
                   f_aux_d( mc + ( cuf_j - 1 ) * nnp + i*nnr ) = f_aux2_d( cuf_j + it )
                 ENDDO
                   !
              ENDDO
           ENDDO
           !
           ip = ip + 1
           !
        ENDDO
     END IF
  ELSE
     !
     !  "backward" scatter from planes to columns
     !
     IF( nprocp == 1 ) GO TO 20
     !
     !  step two: communication
     !
     gcomm = dfft%comm
     !
     ! JR Note: Holding off staging receives until buffer is packed.
     istat = cudaEventSynchronize( dfft%bevents(batch_id) ) 
     CALL start_clock ('A2A')
#ifdef __IPC
     ! TODO: possibly remove this barrier
     call MPI_Barrier( gcomm, ierr )
#endif
     req_cnt = 0

     DO iter = 2, nprocp
        IF(IAND(nprocp, nprocp-1) == 0) THEN
          sorc = IEOR( me-1, iter-1 )
        ELSE
          sorc = MOD(me-1 - (iter-1) + nprocp, nprocp)
        ENDIF
#ifdef __IPC
        IF(dfft%IPC_PEER( sorc + 1 ) .eq. 0) THEN
#endif
#ifdef __GPU_MPI
           call MPI_IRECV( f_aux_d((sorc)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, dfft%srh(req_cnt+1, batch_id), ierr )
#else
           call MPI_IRECV( f_aux((sorc)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, sorc, 0, gcomm, dfft%srh(req_cnt+1, batch_id), ierr )
#endif
           req_cnt = req_cnt + 1
#ifdef __IPC
        ENDIF
#endif

     ENDDO

     DO iter = 2, nprocp
        IF(IAND(nprocp, nprocp-1) == 0) THEN
          dest = IEOR( me-1, iter-1 )
        ELSE
          dest = MOD(me-1 + (iter-1), nprocp)
        ENDIF
#ifdef __IPC
        IF(dfft%IPC_PEER( dest + 1 ) .eq. 1) THEN
           call ipc_send( f_aux2_d((dest)*sendsiz + 1), sendsiz, f_aux_d((me-1)*sendsiz + 1), 0, dest, gcomm, ierr )
        ELSE
#endif
#ifdef __GPU_MPI
           call MPI_ISEND( f_aux2_d((dest)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, dest, 0, gcomm, dfft%srh(req_cnt+1, batch_id), ierr )
#else
           call MPI_ISEND( f_aux2((dest)*sendsiz + 1), sendsiz, MPI_DOUBLE_COMPLEX, dest, 0, gcomm, dfft%srh(req_cnt+1, batch_id), ierr )
#endif
           req_cnt = req_cnt + 1
#ifdef __IPC
        ENDIF
#endif
     ENDDO

     offset = 0
     DO proc = 1, me-1
        offset = offset + npp_ ( proc )
     ENDDO
     istat = cudaMemcpy2DAsync( f_in_d(offset + 1), nr3x, f_aux2_d((me-1)*sendsiz + 1), nppx, npp_(me), batchsize * ncpx, &
     cudaMemcpyDeviceToDevice, dfft%bstreams(batch_id) )

     if(req_cnt .gt. 0) then
        call MPI_WAITALL(req_cnt, dfft%srh(1:req_cnt, batch_id), MPI_STATUSES_IGNORE, ierr)
     endif
#ifdef __IPC
     call sync_ipc_sends( gcomm )
     call MPI_Barrier( gcomm, ierr )
#endif
     CALL stop_clock ('A2A')

     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
     !
     !  step one: store contiguously the columns
     !
     !! f_in = 0.0_DP
     !
     offset = 0

     DO gproc = 1, nprocp
        !
        kdest = ( gproc - 1 ) * sendsiz
        kfrom = offset
        !
        IF (gproc .ne. me) THEN
#ifdef __GPU_MPI

!         This commented code is left here for helping understand the following calls to CUDA APIs
!
!!$cuf kernel do(2) <<<*,*, 0, dfft%bstreams(batch_id)>>>
!         !DO k = 1, ncp_ (me)
!         DO k = 1, batchsize * ncpx
!            DO i = 1, npp_ ( gproc )
!              f_in_d( kfrom + i + (k-1)*nr3x ) = f_aux_d( kdest + i + (k-1)*nppx )
!            END DO
!         END DO
          istat = cudaMemcpy2DAsync( f_in_d(kfrom +1 ), nr3x, f_aux_d(kdest + 1), nppx, npp_(gproc), batchsize * ncpx, &
          cudaMemcpyDeviceToDevice, dfft%bstreams(batch_id) )

#else
#ifdef __IPC
          IF(dfft%IPC_PEER( gproc ) .eq. 1) THEN
               istat = cudaMemcpy2DAsync( f_in_d(kfrom +1 ), nr3x, f_aux_d(kdest + 1), nppx, npp_(gproc), batchsize * ncpx, &
               cudaMemcpyDeviceToDevice, dfft%bstreams(batch_id) )
          ELSE
               istat = cudaMemcpy2DAsync( f_in_d(kfrom +1 ), nr3x, f_aux(kdest + 1), nppx, npp_(gproc), batchsize * ncpx, &
               cudaMemcpyHostToDevice, dfft%bstreams(batch_id) )
          ENDIF
#else

          istat = cudaMemcpy2DAsync( f_in_d(kfrom +1 ), nr3x, f_aux(kdest + 1), nppx, npp_(gproc), batchsize * ncpx, &
          cudaMemcpyHostToDevice, dfft%bstreams(batch_id) )
#endif
#endif
        ENDIF
        offset = offset + npp_ ( gproc )
     ENDDO

20   CONTINUE

  ENDIF

  !istat = cudaDeviceSynchronize()

#endif

  RETURN

END SUBROUTINE fft_scatter_batch_b_gpu

!
!=----------------------------------------------------------------------=!
END MODULE scatter_mod_2d_gpu
!=----------------------------------------------------------------------=!
#endif
