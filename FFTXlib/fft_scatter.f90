!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
! Basic scatter operations needed by parallel FFT
! Written by Carlo Cavazzoni and Stefano de Gironcoli, modified by PG
!
!----------------------------------------------------------------------
!
!=----------------------------------------------------------------------=!
   MODULE fft_scatter
!=----------------------------------------------------------------------=!

        USE fft_types, ONLY: fft_type_descriptor
        USE fft_param

        IMPLICIT NONE

        SAVE

        PRIVATE

        PUBLIC :: fft_scatter_xy, fft_scatter_yz, fft_scatter_many_xy, fft_scatter_many_yz
        PUBLIC :: fft_scatter_tg, fft_scatter_tg_opt

!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!
!
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_xy ( desc, f_in, f_aux, nxx_, isgn, comm )
  !-----------------------------------------------------------------------
  !
  ! Transpose of the fft xy planes across the comm communicator.
  ! If the optional comm is not provided as input, the transpose is made
  ! across desc%comm2 communicator.
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
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in) :: desc
  INTEGER, INTENT(in)           :: nxx_, isgn
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)
  INTEGER, OPTIONAL, INTENT(in) :: comm
  INTEGER :: nr1_temp(1), comm_

#if defined(__MPI)
!$omp master
  CALL start_clock ('fft_scatt_xy')
!$omp end master

  IF (PRESENT(comm)) THEN
     comm_ = comm
  ELSE
     comm_ = desc%comm2
  ENDIF
  !
  if ( abs (isgn) == 1 ) then          ! It's a potential FFT
     CALL impl_xy( MAXVAL ( desc%nr2p ), desc%nproc2, desc%my_nr2p, desc%nr1p, desc%indp, desc%iplp, comm_)
  else if ( abs (isgn) == 2 ) then     ! It's a wavefunction FFT
     CALL impl_xy( MAXVAL ( desc%nr2p ), desc%nproc2, desc%my_nr2p, desc%nr1w, desc%indw, desc%iplw, comm_)
  else if ( abs (isgn) == 3 ) then     ! It's a wavefunction FFT with task group
     ! in task group FFTs whole Y colums are distributed
     nr1_temp = desc%nr1w_tg
     CALL impl_xy( desc%nr2x, 1, desc%nr2x, nr1_temp, desc%indw_tg, desc%iplw, comm_)
  end if
  !
!$omp master
  CALL stop_clock ('fft_scatt_xy')
!$omp end master

  RETURN

  CONTAINS

  SUBROUTINE impl_xy(nr2px, nproc2, my_nr2p, nr1p_, indx, iplx, mpi_comm)
  IMPLICIT NONE
  !
  INTEGER, INTENT(in):: nr2px, nproc2, my_nr2p
  INTEGER, INTENT(in):: nr1p_(nproc2), indx(desc%nr1x,nproc2), iplx(desc%nr1x), mpi_comm
  !
  INTEGER :: ierr, me2, iproc2, ncpx
  INTEGER :: i, it, j, k, kfrom, kdest, mc, m1, m3, i1, icompact, sendsize
  !
#if defined(__NON_BLOCKING_SCATTER)
  INTEGER :: sh(desc%nproc2), rh(desc%nproc2)
#endif

  me2    = desc%mype2 + 1
  ncpx = MAXVAL(nr1p_) * desc%my_nr3p       ! maximum number of Y columns to be disributed
  ! calculate the message size
  sendsize = ncpx * nr2px       ! dimension of the scattered chunks (safe value)

  ierr = 0
  IF (isgn.gt.0) THEN
     !
     IF (nproc2==1) GO TO 10
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
!$omp parallel do collapse(2) private(kdest,kfrom)
     DO iproc2 = 1, nproc2
        DO k = 0, nr1p_(me2)*desc%my_nr3p-1
           kdest = ( iproc2 - 1 ) * sendsize + nr2px * k
           kfrom = desc%nr2p_offset(iproc2) + desc%nr2x *k
           DO i = 1, desc%nr2p( iproc2 )
              f_aux ( kdest + i ) =  f_in ( kfrom + i )
           ENDDO
        ENDDO
     ENDDO
!$omp end parallel do
#if defined(__NON_BLOCKING_SCATTER)
     DO iproc2 = 1, nproc2
        CALL mpi_isend( f_aux( (iproc2-1)*sendsize + 1 ), sendsize, &
                        MPI_DOUBLE_COMPLEX, iproc2-1, me2, mpi_comm, &
                        sh( iproc2 ), ierr )
     ENDDO
#endif
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
                        mpi_comm, rh( iproc2 ), ierr )
        !IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', ' forward receive info<>0', abs(ierr) )
     ENDDO

     call mpi_waitall( nproc2, sh, MPI_STATUSES_IGNORE, ierr )
     !
#else
     CALL mpi_alltoall (f_aux(1), sendsize, MPI_DOUBLE_COMPLEX, f_in(1), &
                        sendsize, MPI_DOUBLE_COMPLEX, mpi_comm, ierr)
     !
     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
#endif
     !
10   CONTINUE
     !
#if defined(__NON_BLOCKING_SCATTER)
     if (nproc2 > 1) call mpi_waitall( nproc2, rh, MPI_STATUSES_IGNORE, ierr )
#endif
     !
!$omp parallel
!$omp do collapse(2) private(it,m3,i1,m1,icompact)
     DO iproc2 = 1, nproc2
        DO i = 0, ncpx-1
           IF(i>=nr1p_(iproc2)*desc%my_nr3p) CYCLE ! control i from 0 to nr1p_(iproc2)*desc%my_nr3p-1
           it = ( iproc2 - 1 ) * sendsize + nr2px * i
           m3 = i/nr1p_(iproc2)+1
           i1 = mod(i,nr1p_(iproc2))+1
           m1 = indx(i1,iproc2)
           icompact = m1 + (m3-1)*desc%nr1x*my_nr2p
           DO j = 1, my_nr2p
              !f_aux( m1 + (j-1)*desc%nr1x + (m3-1)*desc%nr1x*my_nr2p ) = f_in( j + it )
              f_aux( icompact ) = f_in( j + it )
              icompact = icompact + desc%nr1x
           ENDDO
        ENDDO
     ENDDO
!$omp end do nowait
     !
     ! clean extra array elements in each stick
     !
!$omp do
     DO k = 1, my_nr2p*desc%my_nr3p
        DO i1 = 1, desc%nr1x
           IF(iplx(i1)==0) f_aux(desc%nr1x*(k-1)+i1) = (0.0_DP, 0.0_DP)
        ENDDO
     ENDDO
!$omp end do nowait
!$omp end parallel
     !
  ELSE
     !
     !  "backward" scatter from planes to columns
     !
!$omp parallel do collapse(2) private(it,m3,i1,m1,icompact)
     DO iproc2 = 1, nproc2
        DO i = 0, ncpx-1
           IF(i>=nr1p_(iproc2)*desc%my_nr3p) CYCLE ! control i from 0 to nr1p_(iproc2)*desc%my_nr3p-1
           it = ( iproc2 - 1 ) * sendsize + nr2px * i
           m3 = i/nr1p_(iproc2)+1
           i1 = mod(i,nr1p_(iproc2))+1
           m1 = indx(i1,iproc2)
           icompact = m1 + (m3-1)*desc%nr1x*my_nr2p
           DO j = 1, my_nr2p
              !f_in( j + it ) = f_aux( m1 + (j-1)*desc%nr1x + (m3-1)*desc%nr1x*my_nr2p )
              f_in( j + it ) = f_aux( icompact )
              icompact = icompact + desc%nr1x
           ENDDO
        ENDDO
     ENDDO
!$omp end parallel do
#if defined(__NON_BLOCKING_SCATTER)
     DO iproc2 = 1, nproc2
        IF( nproc2 > 1 ) CALL mpi_isend( f_in( ( iproc2 - 1 ) * sendsize + 1 ), &
                              sendsize, MPI_DOUBLE_COMPLEX, iproc2-1, me2, &
                              mpi_comm, sh( iproc2 ), ierr )
     ENDDO
#endif
     IF (nproc2==1) GO TO 20

     !
     !  step two: communication
     !
#if ! defined(__NON_BLOCKING_SCATTER)
     CALL mpi_alltoall (f_in(1), sendsize, MPI_DOUBLE_COMPLEX, f_aux(1), &
                        sendsize, MPI_DOUBLE_COMPLEX, mpi_comm, ierr)

     IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
#else
     DO iproc2 = 1, nproc2
        CALL mpi_irecv( f_aux( (iproc2-1)*sendsize + 1 ), sendsize, &
                        MPI_DOUBLE_COMPLEX, iproc2-1, MPI_ANY_TAG, &
                        mpi_comm, rh(iproc2), ierr )
        ! IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', ' backward receive info<>0', abs(ierr) )
     ENDDO

     call mpi_waitall( nproc2, sh, MPI_STATUSES_IGNORE, ierr )
     call mpi_waitall( nproc2, rh, MPI_STATUSES_IGNORE, ierr )
#endif
     !
     !  step one: store contiguously the columns
     !
!$omp parallel do collapse(2) private(kdest,kfrom)
     DO iproc2 = 1, nproc2
        DO k = 0, nr1p_(me2)*desc%my_nr3p-1
           kdest = ( iproc2 - 1 ) * sendsize + nr2px * k
           kfrom = desc%nr2p_offset(iproc2) + desc%nr2x * k
           DO i = 1, desc%nr2p( iproc2 )
              f_in ( kfrom + i ) = f_aux ( kdest + i )
           ENDDO
        ENDDO
     ENDDO
!$omp end parallel do
     !
     ! clean extra array elements in each stick
     !
     IF( desc%nr2x /= desc%nr2 ) THEN
        DO k = 1, nr1p_(me2)*desc%my_nr3p
           f_in(desc%nr2x*(k-1)+desc%nr2+1:desc%nr2x*k) = (0.0_DP, 0.0_DP)
        ENDDO
     ENDIF

20   CONTINUE

  ENDIF

  END SUBROUTINE impl_xy

#endif

END SUBROUTINE fft_scatter_xy
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_many_xy ( desc, f_in, f_aux, isgn, howmany)
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
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in), TARGET :: desc
  INTEGER, INTENT(in)                            :: isgn
  INTEGER, INTENT(in)                            :: howmany

  COMPLEX (DP), INTENT(inout)   :: f_in (:), f_aux (:)

#if defined(__MPI)
  !
  INTEGER :: ierr, me2, nproc2, iproc2, ncpx, nr2px, ip, icompact
  INTEGER :: i, j, it, it0, k, kfrom, kdest, offset, m1, m2, i1, sendsize
  INTEGER, POINTER :: iplx(:), nr1p_(:), indx(:,:)
  !
  me2    = desc%mype2 + 1
  nproc2 = desc%nproc2
  !
  if ( abs (isgn) == 1 ) then      ! It's a potential FFT
     iplx     => desc%iplp
     nr1p_    => desc%nr1p
     indx     => desc%indp
  else if ( abs (isgn) == 2 ) then ! It's a wavefunction FFT
     iplx     => desc%iplw
     nr1p_    => desc%nr1w
     indx     => desc%indw
  else if ( abs (isgn) == 3 ) then ! It's a wavefunction FFT with task group
     print *, "ERRORE, this should never happen!"
  end if
  !
  CALL start_clock ('fft_scatt_many_xy')
  !
  ! calculate the message size
  !
  nr2px = MAXVAL ( desc%nr2p )                  ! maximum number of Y values to be exchanged
  ncpx  = MAXVAL ( nr1p_ ) * MAXVAL(desc%nr3p)  ! maximum number of Y columns to be exchanged

  sendsize = howmany * ncpx * nr2px             ! dimension of the scattered chunks
  !
  ierr = 0
  !
  IF (isgn.gt.0) THEN
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
     offset = 0

     IF (nproc2==1) THEN
!$omp parallel default(none)                                   &
!$omp          private(i, j, k, ip, it, m1, m2, i1, icompact)  &
!$omp          shared(howmany, ncpx, nr2px, nr1p_, iplx)       &
!$omp&         shared(desc, f_aux, f_in, indx)
        ip = nr1p_(1) * desc%my_nr3p + 1
!$omp do
        DO k=0, howmany-1
           DO j = 1, ncpx
              IF ( j >= ip ) CYCLE
              it = (j - 1) * nr2px + k * ncpx * nr2px
              m2 = (j - 1) / nr1p_(1) + 1
              i1 = mod((j - 1),nr1p_(1)) + 1
              m1 = indx(i1,1)
              icompact = m1 + (m2 - 1) * desc%nr1x * desc%my_nr2p
              DO i = 1, desc%my_nr2p
                 f_aux( icompact + k * desc%nnr ) = f_in( i + it )
                 icompact = icompact + desc%nr1x
              ENDDO
           ENDDO
        ENDDO
!$omp end do
!$omp do
        DO k=0, howmany-1
           DO j = 1, desc%my_nr2p*desc%my_nr3p
              DO i1 = 1, desc%nr1x
                 IF(iplx(i1)==0) f_aux(desc%nr1x*(j-1)+i1+k*desc%nnr) = (0.0_DP, 0.0_DP)
              ENDDO
           ENDDO
        ENDDO
!$omp end do
!$omp end parallel
     ELSE
!$omp parallel default(none)                                  &
!$omp          private(iproc2, i, j, k, kdest, kfrom)         &
!$omp          shared(nproc2, howmany, ncpx, sendsize, nr1p_) &
!$omp&         shared(nr2px, desc, f_aux, f_in, me2)
!$omp do
        DO k = 0, howmany-1
           DO iproc2 = 1, nproc2
              DO j = 0, nr1p_(me2) * desc%my_nr3p - 1
                 kdest = ( iproc2 - 1 ) * sendsize + nr2px * (j + k*ncpx)
                 kfrom = desc%nr2p_offset(iproc2)  + desc%nr2x * (j + k*ncpx)
                    DO i = 1, desc%nr2p( iproc2 )
                       f_aux ( kdest + i ) =  f_in ( kfrom + i )
                    ENDDO
              ENDDO
           ENDDO
        ENDDO
!$omp end do
!$omp end parallel

        CALL mpi_alltoall (f_aux(1), sendsize, MPI_DOUBLE_COMPLEX, f_in(1), &
                           sendsize, MPI_DOUBLE_COMPLEX, desc%comm2, ierr)

        IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )

!$omp parallel default(none)                                               &
!$omp          private(iproc2, i, j, k, ip, it, m1, m2, icompact, it0, i1) &
!$omp          shared(nproc2, howmany, ncpx, sendsize, nr2px, indx)        &
!$omp&         shared(desc, f_aux, f_in, iplx, nr1p_)
!$omp do
        DO k=0, howmany-1
           DO iproc2 = 1, nproc2
              ip = nr1p_(iproc2) * desc%my_nr3p + 1
              it0 = ( iproc2 - 1 ) * sendsize
              DO j = 1, ncpx
                 IF ( j >= ip ) CYCLE
                 it = it0 + (j - 1) * nr2px + k * ncpx * nr2px
                 m2 = (j - 1) / nr1p_(iproc2) + 1
                 i1 = mod((j - 1),nr1p_(iproc2)) + 1
                 m1 = indx(i1,iproc2)
                 icompact = m1 + (m2 - 1) * desc%nr1x * desc%my_nr2p
                 DO i = 1, desc%my_nr2p
                    f_aux( icompact + k*desc%nnr ) = f_in( i + it )
                    icompact = icompact + desc%nr1x
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
!$omp end do
!$omp do
        DO k=0, howmany-1
           DO j = 1, desc%my_nr2p*desc%my_nr3p
              DO i1 = 1, desc%nr1x
                 IF(iplx(i1)==0) f_aux(desc%nr1x*(j-1)+i1+k*desc%nnr) = (0.0_DP, 0.0_DP)
              ENDDO
           ENDDO
        ENDDO
!$omp end do
!$omp end parallel
     ENDIF
     !
  ELSE
     !
     !  "backward" scatter from planes to columns
     !
     IF( nproc2 == 1 ) THEN
!$omp parallel default(none)                                   &
!$omp          private(i, j, k, ip, it, m1, m2, i1, icompact)  &
!$omp          shared(howmany, ncpx, nr2px, nr1p_)             &
!$omp&         shared(desc, f_aux, f_in, indx)
        ip = nr1p_(1) * desc%my_nr3p + 1
!$omp do
        DO k = 0, howmany - 1
           DO j = 1, ncpx
              IF ( j >= ip ) CYCLE
              it = (j - 1) * nr2px + k * ncpx * nr2px
              m2 = (j - 1) / nr1p_(1) + 1
              i1 = mod((j - 1),nr1p_(1)) + 1
              m1 = indx(i1,1)
              icompact = m1 + (m2 - 1) * desc%nr1x * desc%my_nr2p
              DO i = 1, desc%my_nr2p
                 f_in( i + it ) = f_aux( icompact + k * desc%nnr )
                 icompact = icompact + desc%nr1x
              ENDDO
           ENDDO
        ENDDO
!$omp end do
!$omp end parallel
     ELSE
!$omp parallel default(none)                                               &
!$omp          private(iproc2, i, j, k, ip, it, m1, m2, icompact, it0, i1) &
!$omp          shared(nproc2, howmany, ncpx, sendsize, nr2px, indx)        &
!$omp&         shared(desc, f_aux, f_in, nr1p_)
!$omp do
        DO k=0, howmany-1
           DO iproc2 = 1, nproc2
              ip = nr1p_(iproc2) * desc%my_nr3p + 1
              it0 = ( iproc2 - 1 ) * sendsize
              DO j = 1, ncpx
                 IF ( j >= ip ) CYCLE
                 it = it0 + (j - 1) * nr2px + k * ncpx * nr2px
                 m2 = (j - 1) / nr1p_(iproc2) + 1
                 i1 = mod((j - 1),nr1p_(iproc2))+1
                 m1 = indx(i1,iproc2)
                 icompact = m1 + (m2 - 1) * desc%nr1x * desc%my_nr2p
                 DO i = 1, desc%my_nr2p
                    f_in( i + it ) = f_aux( icompact + k*desc%nnr )
                    icompact = icompact + desc%nr1x
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
!$omp end do
!$omp end parallel
        !
        !
        !  step two: communication
        !
        CALL mpi_alltoall (f_in(1), sendsize, MPI_DOUBLE_COMPLEX, f_aux(1), &
                           sendsize, MPI_DOUBLE_COMPLEX, desc%comm2, ierr)

        IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
        !
        !  step one: store contiguously the columns
        !
!$omp parallel default(none)                                  &
!$omp          private(iproc2, i, j, k, kdest, kfrom)         &
!$omp          shared(nproc2, howmany, ncpx, sendsize, nr2px) &
!$omp&         shared(desc, f_aux, f_in, nr1p_, me2)
!$omp do
        DO k = 0, howmany-1
           DO iproc2 = 1, nproc2
              DO j = 0, nr1p_(me2) * desc%my_nr3p - 1
                 kdest = ( iproc2 - 1 ) * sendsize + nr2px * (j + k*ncpx)
                 kfrom = desc%nr2p_offset(iproc2)  + desc%nr2x * (j + k*ncpx)
                    DO i = 1, desc%nr2p( iproc2 )
                       f_in ( kfrom + i ) =  f_aux ( kdest + i )
                    ENDDO
              ENDDO
           ENDDO
        ENDDO
!$omp end do
           ! clean extra array elements in each stick
        IF( desc%nr2x /= desc%nr2 ) THEN
!$omp do
           DO k = 0, howmany-1
              DO j = 1, nr1p_(me2)*desc%my_nr3p
                 DO i = desc%nr2+1, desc%nr2x
                    f_in( k*ncpx*desc%nr2x + (j-1)*desc%nr2x + i) = 0.0d0
                 END DO
              END DO
           ENDDO
!$omp end do
        ENDIF
!$omp end parallel
     ENDIF
  ENDIF

  CALL stop_clock ('fft_scatt_many_xy')

#endif

  RETURN
98 format ( 10 ('(',2f12.9,')') )
99 format ( 20 ('(',2f12.9,')') )

END SUBROUTINE fft_scatter_many_xy
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_yz ( desc, f_in, f_aux, nxx_, isgn )
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
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in) :: desc
  INTEGER, INTENT(in)           :: nxx_, isgn
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)

#if defined(__MPI)
  !
  CALL start_clock ('fft_scatt_yz')

  if ( abs (isgn) == 1 ) then      ! It's a potential FFT
     CALL impl_yz(desc%mype2+1, desc%mype2+1, desc%nsp, desc%ir1p, desc%nsp_offset)
  else if ( abs (isgn) == 2 ) then ! It's a wavefunction FFT
     CALL impl_yz(desc%mype2+1, desc%mype2+1, desc%nsw, desc%ir1w, desc%nsw_offset)
  else if ( abs (isgn) == 3 ) then ! It's a wavefunction FFT with task group
     CALL impl_yz(1, desc%nproc2, desc%nsw, desc%ir1w_tg, desc%nsw_offset)
  end if

  CALL stop_clock ('fft_scatt_yz')

  RETURN

  CONTAINS

  SUBROUTINE impl_yz(me2_start, me2_end, ncp_, ir1p_, me2_iproc3_offset)
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: me2_start, me2_end
  INTEGER, INTENT(in) :: ncp_(desc%nproc), ir1p_(desc%nr1x), me2_iproc3_offset(desc%nproc2, desc%nproc3)
  !
  INTEGER :: ierr, me2, me3, nproc3, iproc3, ncpx, nr3px, ip
  INTEGER :: i, it, k, kfrom, kdest, mc, m1, m2, i1, sendsize
  INTEGER :: my_nr1p_
  !
#if defined(__NON_BLOCKING_SCATTER)
  INTEGER :: sh(desc%nproc3), rh(desc%nproc3)
#endif
  !
  me3    = desc%mype3 + 1
  nproc3 = desc%nproc3
  !
  my_nr1p_ = count (ir1p_ > 0)
  !
  ! calculate the message size
  !
  nr3px = MAXVAL ( desc%nr3p )  ! maximum number of Z values to be exchanged
  ncpx  = MAXVAL ( ncp_ )       ! maximum number of Z columns to be exchanged
  sendsize = ncpx * nr3px * ( me2_end - me2_start + 1 )  ! dimension of the scattered chunks
  !
  ierr = 0
  IF (isgn.gt.0) THEN
     !
     IF (nproc3==1) GO TO 10
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
!$omp parallel do collapse(3) private(kdest,kfrom)
     DO iproc3 = 1, nproc3
        DO me2 = me2_start, me2_end
           DO k = 0, ncpx-1   ! was ncp_(me3)
              IF (k>=ncp_(desc%iproc(me2,me3))) CYCLE ! control k from 0 to ncp_(desc%iproc(me2,me3))-1
              kdest = ( iproc3 - 1 ) * sendsize + nr3px * ( me2_iproc3_offset( me2 - me2_start + 1, me3 ) + k )
              kfrom = desc%nr3p_offset(iproc3) + desc%nr3x * ( me2_iproc3_offset( me2 - me2_start + 1, me3 ) + k )
              DO i = 1, desc%nr3p( iproc3 )
                 f_aux ( kdest + i ) =  f_in ( kfrom + i )
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!$omp end parallel do
#if defined(__NON_BLOCKING_SCATTER)
     DO iproc3 = 1, nproc3
        CALL mpi_isend( f_aux( ( iproc3 - 1 ) * sendsize + 1 ), sendsize, &
                        MPI_DOUBLE_COMPLEX, iproc3-1, me3, desc%comm3, &
                        sh( iproc3 ), ierr )
     ENDDO
#endif
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
10   CONTINUE
     !
#if defined(__NON_BLOCKING_SCATTER)
     if (nproc3 > 1) call mpi_waitall( nproc3, rh, MPI_STATUSES_IGNORE, ierr )
#endif
     !
!$omp parallel
     ! ensures that no garbage is present in the output
     !
!$omp do
     DO k = 1, desc%my_nr3p*my_nr1p_*desc%nr2x
        f_aux(k) = (0.0_DP, 0.0_DP)
     ENDDO
!$omp end do
!$omp do collapse(3) private(ip,it,mc,m1,m2,i1)
     DO iproc3 = 1, desc%nproc3
        DO me2 = me2_start, me2_end
           DO i = 1, ncpx ! was ncp_(iproc3)
              ip = desc%iproc( me2, iproc3)
              IF ( i>ncp_(ip) ) CYCLE
              it = ( iproc3 - 1 ) * sendsize + nr3px * ( me2_iproc3_offset( me2 - me2_start + 1, iproc3 ) + i - 1 )
              mc = desc%ismap( i + desc%iss(ip) ) ! this is  m1+(m2-1)*nr1x  of the  current pencil
              m1 = mod (mc-1,desc%nr1x) + 1
              m2 = (mc-1)/desc%nr1x + 1
              i1 = m2 + ( ir1p_(m1) - 1 ) * desc%nr2x
              DO k = 1, desc%my_nr3p
                 f_aux( i1 ) = f_in( k + it )
                 i1 = i1 + desc%nr2x*my_nr1p_
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!$omp end do nowait
!$omp end parallel
     !
  ELSE
     !
     !  "backward" scatter from planes to columns
     !
!$omp parallel do collapse(3) private(ip,it,mc,m1,m2,i1)
     DO iproc3 = 1, desc%nproc3
        DO me2 = me2_start, me2_end
           DO i = 1, ncpx
              ip = desc%iproc( me2, iproc3)
              IF ( i>ncp_(ip) ) CYCLE
              it = ( iproc3 - 1 ) * sendsize + nr3px * ( me2_iproc3_offset( me2 - me2_start + 1, iproc3 ) + i - 1 )
              mc = desc%ismap( i + desc%iss(ip) ) ! this is  m1+(m2-1)*nr1x  of the  current pencil
              m1 = mod (mc-1,desc%nr1x) + 1
              m2 = (mc-1)/desc%nr1x + 1
              i1 = m2 + ( ir1p_(m1) - 1 ) * desc%nr2x
              DO k = 1, desc%my_nr3p
                 f_in( k + it ) = f_aux( i1 )
                 i1 = i1 + desc%nr2x * my_nr1p_
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!$omp end parallel do
#if defined(__NON_BLOCKING_SCATTER)
     DO iproc3 = 1, desc%nproc3
        IF( nproc3 > 1 ) CALL mpi_isend( f_in( ( iproc3 - 1 ) * sendsize + 1 ), &
                                        sendsize, MPI_DOUBLE_COMPLEX, iproc3-1, &
                                        me3, desc%comm3, sh( iproc3 ), ierr )
     ENDDO
#endif

     IF( nproc3 == 1 ) GO TO 20
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
!$omp parallel do collapse(3) private(kdest,kfrom)
     DO iproc3 = 1, nproc3
        DO me2 = me2_start, me2_end
           DO k = 0, ncpx-1   ! was ncp_(me3)
              IF (k>=ncp_(desc%iproc(me2,me3))) CYCLE ! control k from 0 to ncp_(desc%iproc(me2,me3))-1
              kdest = ( iproc3 - 1 ) * sendsize + nr3px * ( me2_iproc3_offset( me2 - me2_start + 1, me3 ) + k )
              kfrom = desc%nr3p_offset(iproc3) + desc%nr3x * ( me2_iproc3_offset( me2 - me2_start + 1, me3 ) + k )
              DO i = 1, desc%nr3p( iproc3 )
                 f_in ( kfrom + i ) = f_aux ( kdest + i )
              ENDDO
           ENDDO
        ENDDO
     ENDDO
!$omp end parallel do

     ! clean extra array elements in each stick

     IF( desc%nr3x /= desc%nr3 ) THEN
        DO k = 1, ncp_(desc%mype+1)
           f_in(desc%nr3x*(k-1)+desc%nr3+1:desc%nr3x*k) = (0.0_DP, 0.0_DP)
        END DO
     END IF

20   CONTINUE

  ENDIF

  END SUBROUTINE impl_yz

#endif

END SUBROUTINE fft_scatter_yz
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_many_yz ( desc, f_in, f_aux, isgn, howmany )
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
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in), TARGET :: desc
  COMPLEX (DP), INTENT(inout)                    :: f_in(:), f_aux(:)
  INTEGER, INTENT(in)                            :: isgn, howmany

#if defined(__MPI)
  !
  INTEGER :: ierr, me, me2, me3, nproc3, iproc3, ncpx, nr3px, ip
  INTEGER :: i, j, it, it0, k, kfrom, kdest, mc, m1, m2, i1, sendsize
  INTEGER, POINTER :: ncp_(:), ir1p_(:)
  INTEGER :: my_nr1p_
  !
  me     = desc%mype  + 1
  me2    = desc%mype2 + 1
  me3    = desc%mype3 + 1
  nproc3 = desc%nproc3
  !
  if ( abs (isgn) == 1 ) then      ! It's a potential FFT
     ncp_     => desc%nsp
     ir1p_    => desc%ir1p
     my_nr1p_ = count (desc%ir1p   > 0)
  else if ( abs (isgn) == 2 ) then ! It's a wavefunction FFT
     ncp_     => desc%nsw
     ir1p_    => desc%ir1w
     my_nr1p_ = count (desc%ir1w   > 0)
  else if ( abs (isgn) == 3 ) then ! It's a wavefunction FFT with task group
     print *, "ERRORE, this should never happen!"
  end if
  !
  CALL start_clock ('fft_scatt_many_yz')
  !
  ! calculate the message size
  !
  nr3px = MAXVAL ( desc%nr3p )  ! maximum number of Z values to be exchanged
  ncpx  = MAXVAL ( ncp_ )       ! maximum number of Z columns to be exchanged

  sendsize = howmany * ncpx * nr3px       ! dimension of the scattered chunks
  !
  ierr = 0
  IF (isgn.gt.0) THEN
     !
     ! "forward" scatter from columns to planes
     !
     ! step one: store contiguously the slices
     !
     IF (nproc3==1) THEN
!$omp parallel default(none)                                   &
!$omp          private(i, j, k, it, mc, m1, m2, i1)            &
!$omp          shared(howmany, ncpx, nr3px, desc, f_aux, f_in) &
!$omp&         shared(ncp_, me2, ir1p_, my_nr1p_)
!$omp do collapse(2)
        DO k=0, howmany-1
           DO j = 1, desc%my_nr3p*desc%nr2x*my_nr1p_
              f_aux(j+k*desc%nnr) = (0.0_DP, 0.0_DP)
           ENDDO
        ENDDO
!$omp end do
!$omp do collapse(2)
        DO k=0, howmany-1
           DO j = 1, ncp_(desc%iproc( me2, 1))
              it = (j-1)*nr3px + k*ncpx*nr3px
              mc = desc%ismap( j + desc%iss(desc%iproc( me2, 1)) )
              m1 = mod (mc-1,desc%nr1x) + 1
              m2 = (mc-1)/desc%nr1x + 1
              i1 = m2 + ( ir1p_(m1) - 1 ) * desc%nr2x
              DO i = 1, desc%my_nr3p
                 f_aux( i1 + k*desc%nnr ) = f_in( it + i)
                 i1 = i1 + desc%nr2x * my_nr1p_
              ENDDO
           ENDDO
        ENDDO
!$omp end do
!$omp end parallel
     ELSE
!$omp parallel default(none)                                   &
!$omp          private(i, j, k, kdest, kfrom, iproc3)          &
!$omp&         private(it, mc, m1, m2, i1)                     &
!$omp&         shared(nproc3, howmany, ncpx, sendsize, nr3px)  &
!$omp&         shared(ncp_, ir1p_, my_nr1p_, ierr)             &
!$omp&         shared(desc, f_aux, f_in, me2, me3)
!$omp do collapse(3)
        DO k = 0, howmany-1
           DO iproc3 = 1, nproc3
              DO j = 0, ncp_(desc%iproc(me2,me3))-1
                 kdest = ( iproc3 - 1 ) * sendsize + nr3px * (j + k*ncpx)
                 kfrom = desc%nr3p_offset(iproc3)  + desc%nr3x * (j + k*ncpx)
                 DO i = 1, desc%nr3p( iproc3 )
                    f_aux ( kdest + i ) =  f_in ( kfrom + i )
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
!$omp end do
        !
        ! ensures that no garbage is present in the output
        ! useless; the later accessed elements are overwritten by the A2A step
        !
        ! step two: communication  across the    nproc3    group
        !
!$omp single
        CALL mpi_alltoall (f_aux(1), sendsize, MPI_DOUBLE_COMPLEX, f_in(1), &
                           sendsize, MPI_DOUBLE_COMPLEX, desc%comm3, ierr)

        IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
!$omp end single
        !
!$omp do collapse(2)
        DO k=0, howmany-1
           DO j = 1, desc%my_nr3p*desc%nr2x*my_nr1p_
              f_aux(j+k*desc%nnr) = (0.0_DP, 0.0_DP)
           ENDDO
        ENDDO
!$omp end do
!$omp do collapse(3)
        DO k=0, howmany-1
           DO iproc3 = 1, nproc3
              DO j = 1, ncpx
                 IF ( j>ncp_(desc%iproc( me2, iproc3)) ) CYCLE
                 it = (iproc3 - 1) * sendsize + (j-1)*nr3px + k*ncpx*nr3px
                 mc = desc%ismap( j + desc%iss(desc%iproc( me2, iproc3)) )
                 m1 = mod (mc-1,desc%nr1x) + 1
                 m2 = (mc-1)/desc%nr1x + 1
                 i1 = m2 + ( ir1p_(m1) - 1 ) * desc%nr2x
                 DO i = 1, desc%my_nr3p
                    f_aux( i1 + k*desc%nnr ) = f_in( it + i)
                    i1 = i1 + desc%nr2x * my_nr1p_
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
!$omp end do
!$omp end parallel
     ENDIF
     !
  ELSE
     !
     !  "backward" scatter from planes to columns
     !
     IF( nproc3 == 1 ) THEN
!$omp parallel default(none)                                   &
!$omp          private(i, j, k, it, mc, m1, m2, i1)            &
!$omp          shared(howmany, ncpx, nr3px, desc, f_aux, f_in) &
!$omp&         shared(ncp_, me2, ir1p_, my_nr1p_)
!$omp do collapse(2)
        DO k = 0, howmany - 1
           DO j = 1, ncp_(desc%iproc(me2,1))
              it = (j-1)*nr3px + k*ncpx*nr3px
              mc = desc%ismap( j + desc%iss(desc%iproc(me2,1)))
              m1 = mod (mc-1,desc%nr1x) + 1
              m2 = (mc-1)/desc%nr1x + 1
              i1 = m2 + ( ir1p_(m1) - 1 ) * desc%nr2x
              DO i = 1, desc%my_nr3p
                 f_in( i + it ) = f_aux( i1  + k*desc%nnr )
                 i1 = i1 + desc%nr2x*my_nr1p_
              ENDDO
           ENDDO
        ENDDO
!$omp end do
!$omp end parallel
     ELSE
     !
!$omp parallel default(none)                                             &
!$omp          private(iproc3, i, j, k, ip, it, mc, m1, m2, i1, it0,kdest,kfrom)     &
!$omp          shared(nproc3, howmany, ncpx, sendsize, nr3px, desc, me3) &
!$omp&         shared(f_aux, f_in, ncp_, me2, ir1p_, my_nr1p_, ierr)
!$omp do collapse(3)
        DO k = 0, howmany - 1
           DO iproc3 = 1, nproc3
              DO j = 1, ncpx
                 IF ( j>ncp_(desc%iproc( me2, iproc3)) ) CYCLE
                 it = (iproc3-1) * sendsize + (j-1)*nr3px + k*ncpx*nr3px
                 mc = desc%ismap( j + desc%iss(desc%iproc( me2, iproc3)) )
                 m1 = mod (mc-1,desc%nr1x) + 1
                 m2 = (mc-1)/desc%nr1x + 1
                 i1 = m2 + ( ir1p_(m1) - 1 ) * desc%nr2x
                 DO i = 1, desc%my_nr3p
                    f_in( i + it ) = f_aux( i1  + k*desc%nnr )
                    i1 = i1 + desc%nr2x*my_nr1p_
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
!$omp end do
        !
!$omp single
        !
        !  step two: communication
        !
        CALL mpi_alltoall (f_in(1), sendsize, MPI_DOUBLE_COMPLEX, f_aux(1), &
                           sendsize, MPI_DOUBLE_COMPLEX, desc%comm3, ierr)

        IF( abs(ierr) /= 0 ) CALL fftx_error__ ('fft_scatter', 'info<>0', abs(ierr) )
!$omp end single
        !
        !  step one: store contiguously the columns
        !
!$omp do collapse(3)
        DO k = 0, howmany-1
           DO iproc3 = 1, nproc3
              DO j = 0, ncp_(desc%iproc(me2,me3))-1
                 kdest = ( iproc3 - 1 ) * sendsize + nr3px * (j + k*ncpx)
                 kfrom = desc%nr3p_offset(iproc3)  + desc%nr3x * (j + k*ncpx)
                 DO i = 1, desc%nr3p( iproc3 )
                    f_in( kfrom + i ) = f_aux( kdest + i )
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
!$omp end do
        ! clean extra array elements in each stick
        IF( desc%nr3x /= desc%nr3 ) THEN
!$omp do collapse(2)
           DO k = 0, howmany-1
              DO j = 1, ncp_ ( desc%mype+1 )
                 DO i = desc%nr3+1, desc%nr3x
                    f_in( k*ncpx*desc%nr3x + (j-1)*desc%nr3x + i) = 0.0d0
                 END DO
              END DO
           END DO
!$omp end do
        ENDIF
!$omp end parallel
     ENDIF
     !
  ENDIF

  CALL stop_clock ('fft_scatt_many_yz')

#endif

  RETURN
98 format ( 10 ('(',2f12.9,')') )
99 format ( 20 ('(',2f12.9,')') )

END SUBROUTINE fft_scatter_many_yz

!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_tg ( desc, f_in, f_aux, nxx_, isgn )
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
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in) :: desc
  INTEGER, INTENT(in)           :: nxx_, isgn
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)

  INTEGER :: ierr

  CALL start_clock ('fft_scatt_tg')

  if ( abs (isgn) /= 3 ) call fftx_error__ ('fft_scatter_tg', 'wrong call', 1 )

#if defined(__MPI)
  !
  f_aux = f_in
  if ( isgn > 0 ) then

     CALL MPI_ALLTOALLV( f_aux,  desc%tg_snd, desc%tg_sdsp, MPI_DOUBLE_COMPLEX, &
                         f_in, desc%tg_rcv, desc%tg_rdsp, MPI_DOUBLE_COMPLEX, desc%comm2, ierr)
     IF( ierr /= 0 ) CALL fftx_error__( 'fft_scatter_tg', ' alltoall error 1 ', abs(ierr) )

  else

     CALL MPI_ALLTOALLV( f_aux,  desc%tg_rcv, desc%tg_rdsp, MPI_DOUBLE_COMPLEX, &
                         f_in, desc%tg_snd, desc%tg_sdsp, MPI_DOUBLE_COMPLEX, desc%comm2, ierr)
     IF( ierr /= 0 ) CALL fftx_error__( 'fft_scatter_tg', ' alltoall error 2 ', abs(ierr) )
   end if

#endif
  CALL stop_clock ('fft_scatt_tg')

  RETURN

END SUBROUTINE fft_scatter_tg
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_tg_opt ( desc, f_in, f_out, nxx_, isgn )
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
  IMPLICIT NONE

  TYPE (fft_type_descriptor), INTENT(in) :: desc
  INTEGER, INTENT(in)           :: nxx_, isgn
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_out (nxx_)

  INTEGER :: ierr

  CALL start_clock ('fft_scatt_tg')

  if ( abs (isgn) /= 3 ) call fftx_error__ ('fft_scatter_tg', 'wrong call', 1 )

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
  CALL stop_clock ('fft_scatt_tg')

  RETURN

END SUBROUTINE fft_scatter_tg_opt
!
!=----------------------------------------------------------------------=!
   END MODULE fft_scatter
!=----------------------------------------------------------------------=!
