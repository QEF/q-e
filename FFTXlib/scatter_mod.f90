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
! Rewritten by Stefano de Gironcoli
!----------------------------------------------------------------------
!
!=----------------------------------------------------------------------=!
   MODULE scatter_mod
!=----------------------------------------------------------------------=!

        USE fft_types, ONLY: fft_type_descriptor
        USE fft_param

        IMPLICIT NONE

        INTERFACE gather_grid
           MODULE PROCEDURE gather_real_grid, gather_complex_grid
        END INTERFACE

        INTERFACE scatter_grid
           MODULE PROCEDURE scatter_real_grid, scatter_complex_grid
        END INTERFACE

        SAVE

        PRIVATE

        PUBLIC :: fft_type_descriptor
        PUBLIC :: gather_grid, scatter_grid
        PUBLIC :: fft_scatter_xy, fft_scatter_yz, fft_scatter_tg, fft_scatter_tg_opt
        PUBLIC :: cgather_sym, cgather_sym_many, cscatter_sym_many

!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!
!
!
!-----------------------------------------------------------------------
SUBROUTINE fft_scatter_xy ( desc, f_in, f_aux, nxx_, isgn )
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

  TYPE (fft_type_descriptor), INTENT(in) :: desc
  INTEGER, INTENT(in)           :: nxx_, isgn
  COMPLEX (DP), INTENT(inout)   :: f_in (nxx_), f_aux (nxx_)

  INTEGER :: nr1_temp(1)

#if defined(__MPI)
  CALL start_clock ('fft_scatt_xy')
  !
  if ( abs (isgn) == 1 ) then          ! It's a potential FFT
     CALL impl_xy( MAXVAL ( desc%nr2p ), desc%nproc2, desc%my_nr2p, desc%nr1p, desc%indp, desc%iplp)
  else if ( abs (isgn) == 2 ) then     ! It's a wavefunction FFT
     CALL impl_xy( MAXVAL ( desc%nr2p ), desc%nproc2, desc%my_nr2p, desc%nr1w, desc%indw, desc%iplw)
  else if ( abs (isgn) == 3 ) then     ! It's a wavefunction FFT with task group
     ! in task group FFTs whole Y colums are distributed
     nr1_temp = desc%nr1w_tg
     CALL impl_xy( desc%nr2x, 1, desc%nr2x, nr1_temp, desc%indw_tg, desc%iplw)
  end if
  !
  CALL stop_clock ('fft_scatt_xy')

  RETURN

  CONTAINS

  SUBROUTINE impl_xy(nr2px, nproc2, my_nr2p, nr1p_, indx, iplx)
  IMPLICIT NONE
  !
  INTEGER, INTENT(in):: nr2px, nproc2, my_nr2p
  INTEGER, INTENT(in):: nr1p_(nproc2), indx(desc%nr1x,nproc2), iplx(desc%nr1x)
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
                        MPI_DOUBLE_COMPLEX, iproc2-1, me2, desc%comm2, &
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
                              desc%comm2, sh( iproc2 ), ierr )
     ENDDO
#endif
     IF (nproc2==1) GO TO 20

     !
     !  step two: communication
     !
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
     CALL impl_yz(desc%mype2+1, desc%mype2+1, desc%nsp, desc%ir1p)
  else if ( abs (isgn) == 2 ) then ! It's a wavefunction FFT
     CALL impl_yz(desc%mype2+1, desc%mype2+1, desc%nsw, desc%ir1w)
  else if ( abs (isgn) == 3 ) then ! It's a wavefunction FFT with task group
     CALL impl_yz(1, desc%nproc2, desc%nsw, desc%ir1w_tg)
  end if

  CALL stop_clock ('fft_scatt_yz')

  RETURN

  CONTAINS

  SUBROUTINE impl_yz(me2_start, me2_end, ncp_, ir1p_)
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: me2_start, me2_end
  INTEGER, INTENT(in) :: ncp_(desc%nproc), ir1p_(desc%nr1x)
  !
  INTEGER :: ierr, me2, me3, nproc3, iproc3, ncpx, nr3px, ip
  INTEGER :: i, it, k, kfrom, kdest, mc, m1, m2, i1, sendsize
  INTEGER, ALLOCATABLE :: me2_offset(:), me2_iproc3_offset(:,:)
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
  ALLOCATE ( me2_offset( me2_end - me2_start + 1 ) )
  ALLOCATE ( me2_iproc3_offset( me2_end - me2_start + 1, nproc3 ) )
  !
  me2_offset(1) = 0
  me2_iproc3_offset(1,1:nproc3) = 0
  DO me2 = me2_start, me2_end-1
    me2_offset(me2+1) = me2_offset(me2) + ncp_(desc%iproc(me2,me3))
    me2_iproc3_offset(me2+1,1:nproc3) = me2_iproc3_offset(me2,1:nproc3) + ncp_(desc%iproc(me2,1:nproc3))
  ENDDO
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
  !
  DEALLOCATE ( me2_offset , me2_iproc3_offset )

  END SUBROUTINE impl_yz

#endif

END SUBROUTINE fft_scatter_yz
!
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
!----------------------------------------------------------------------------
SUBROUTINE gather_real_grid ( dfft, f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers a distributed real-space FFT grid to dfft%root, that is,
  ! ... the first processor of input descriptor dfft - version for real arrays
  !
  ! ... REAL*8  f_in  = distributed variable (dfft%nnr)
  ! ... REAL*8  f_out = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(in) :: f_in (:)
  REAL(DP), INTENT(inout):: f_out(:)
  TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info, offset_in, offset_aux, ir3
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfft%nproc-1), recvcount(0:dfft%nproc-1)
  REAL(DP), ALLOCATABLE ::  f_aux(:) 
  !
  IF( size( f_in ) < dfft%nnr ) &
     CALL fftx_error__( ' gather_real_grid ', ' f_in too small ', dfft%nnr-size( f_in ) )
  !
  CALL start_clock( 'rgather_grid' )

  !write (6,*) 'gather grid ok 0 ', dfft%nproc, dfft%nproc2, dfft%nproc3
  ALLOCATE ( f_aux(dfft%nr1x * dfft%nr2x * dfft%my_nr3p) )
  ! 1) gather within the comm2 communicator
  displs = 0
  DO proc =0, ( dfft%nproc2 -1 )
     recvcount(proc) = dfft%nr1x * dfft%nr2p(proc+1) 
     if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
  end do
  offset_in = 1; offset_aux = 1
  do ir3 = 1, dfft%my_nr3p
    !write (6,*) 'gather grid ok 1 ir3=', ir3
     info = 0
     CALL MPI_GATHERV( f_in(offset_in) , recvcount(dfft%mype2), MPI_DOUBLE_PRECISION, &
                       f_aux(offset_aux), recvcount, displs, MPI_DOUBLE_PRECISION, dfft%root,      &
                       dfft%comm2, info )
     CALL fftx_error__( 'gather_real_grid', 'info<>0', info )
    !write (6,*) 'gather grid ok 2 ir3=', ir3
     offset_in  = offset_in + dfft%nr1x * dfft%my_nr2p
     offset_aux = offset_aux + dfft%nr1x * dfft%nr2
  end do

  ! 2) gather within the comm3 communicator 
  displs = 0
  DO proc = 0, ( dfft%nproc3 - 1 )
     recvcount(proc) = dfft%nr1x * dfft%nr2x * dfft%nr3p(proc+1) 
     if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
  ENDDO
  info = 0
  !write (6,*) 'gather grid ok 3'
  CALL MPI_GATHERV( f_aux, recvcount(dfft%mype3), MPI_DOUBLE_PRECISION, &
                    f_out, recvcount, displs, MPI_DOUBLE_PRECISION, dfft%root,      &
                    dfft%comm3, info )
  !write (6,*) 'gather grid ok 4'
  CALL fftx_error__( ' gather_real_grid', 'info<>0', info )
  !
  ! ... the following check should be performed only on processor dfft%root
  ! ... otherwise f_out must be allocated on all processors even if not used
  !
  info = size( f_out ) - displs( dfft%nproc3-1 ) - recvcount( dfft%nproc3-1 )
  IF( info < 0 ) &
     CALL fftx_error__( ' gather_real_grid ', ' f_out too small ', -info )
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'rgather_grid' )
  !
#else
  CALL fftx_error__(' gather_real_grid', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE gather_real_grid

!----------------------------------------------------------------------------
SUBROUTINE gather_complex_grid ( dfft, f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... gathers a distributed real-space FFT grid to dfft%root, that is,
  ! ... the first processor of input descriptor dfft - complex arrays
  !
  ! ... COMPLEX*16  f_in  = distributed variable (dfft%nnr)
  ! ... COMPLEX*16  f_out = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(in) :: f_in (:)
  COMPLEX(DP), INTENT(inout):: f_out(:)
  TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:) 
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info, offset_in, offset_aux, ir3
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfft%nproc-1), recvcount(0:dfft%nproc-1)
  !
  
  CALL start_clock( 'cgather_grid' )
  !write (*,*) 'gcgather_grid size(f_in),dfft%nnr',size(f_in), dfft%nnr ; FLUSH(6)
  IF( 2*size( f_in ) < dfft%nnr ) &
     CALL fftx_error__( ' gather_complex_grid ', ' f_in too small ', dfft%nnr-size( f_in ) )
  !
  !write (6,*) 'gather grid ok 0 ', dfft%nproc, dfft%nproc2, dfft%nproc3
  ALLOCATE ( f_aux(dfft%nr1x * dfft%nr2x * dfft%my_nr3p ) )
  ! 1) gather within the comm2 communicator
  displs = 0
  DO proc =0, ( dfft%nproc2 -1 )
     recvcount(proc) = 2 * dfft%nr1x * dfft%nr2p(proc+1) 
     if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
  end do
  offset_in = 1; offset_aux = 1
  do ir3 = 1, dfft%my_nr3p
     info = 0
     !write (6,*) 'gather grid ok 1 ir3=', ir3
     CALL MPI_GATHERV( f_in(offset_in) , recvcount(dfft%mype2), MPI_DOUBLE_PRECISION, &
                       f_aux(offset_aux), recvcount, displs, MPI_DOUBLE_PRECISION, dfft%root,      &
                       dfft%comm2, info )
     CALL fftx_error__( 'gather_complex_grid', 'info<>0', info )
     !write (6,*) 'gather grid ok 2 ir3=', ir3
     offset_in  = offset_in + dfft%nr1x * dfft%my_nr2p
     offset_aux = offset_aux + dfft%nr1x * dfft%nr2
  end do

  ! 2) gather within the comm3 communicator 
  displs = 0
  DO proc = 0, ( dfft%nproc3 - 1 )
     recvcount(proc) = 2 * dfft%nr1x * dfft%nr2x * dfft%nr3p(proc+1) 
     if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
  ENDDO
  !
  ! ... the following check should be performed only on processor dfft%root
  ! ... otherwise f_out must be allocated on all processors even if not used
  !
  !write (*,*) 'gcgather_grid 2*size(f_out)',2*size(f_out) ; FLUSH(6)
  !write (*,*) 'gcgather_grid displ+recv',dfft%nproc3, displs(dfft%nproc3-1) + recvcount(dfft%nproc3-1); FLUSH(6)
  info = 2*size( f_out ) - displs( dfft%nproc3 - 1 ) - recvcount( dfft%nproc3-1 ) ; FLUSH(6)
  IF( info < 0 ) CALL fftx_error__( ' gather_complex_grid ', ' f_out too small ', -info )

  info = 0
  !write (6,*) 'gather grid ok 3'
  CALL MPI_GATHERV( f_aux, recvcount(dfft%mype3), MPI_DOUBLE_PRECISION, &
                    f_out, recvcount, displs, MPI_DOUBLE_PRECISION, dfft%root,      &
                    dfft%comm3, info )
  !write (6,*) 'gather grid ok 4'
  CALL fftx_error__( 'gather_complex_grid', 'info<>0', info )
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'cgather_grid' )
  !
#else
  CALL fftx_error__('gather_complex_grid', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE gather_complex_grid

!----------------------------------------------------------------------------
SUBROUTINE scatter_real_grid ( dfft, f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters a real-space FFT grid from dfft%root, first processor of
  ! ... input descriptor dfft, to all others - opposite of "gather_grid"
  !
  ! ... REAL*8  f_in  = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
  ! ... REAL*8  f_out = distributed variable (dfft%nnr)
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(in) :: f_in (:)
  REAL(DP), INTENT(inout):: f_out(:)
  TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
  REAL(DP), ALLOCATABLE ::  f_aux(:) 
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info, offset_in, offset_aux, ir3
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfft%nproc-1), sendcount(0:dfft%nproc-1)
  !
  !
  CALL start_clock( 'rscatter_grid' )
  !
  !write (6,*) 'scatter grid ok 0'
  ALLOCATE ( f_aux(dfft%nr1x * dfft%nr2x * dfft%my_nr3p) )
  ! 1) scatter within the comm3 communicator 
  displs = 0
  DO proc = 0, ( dfft%nproc3 - 1 )
     sendcount(proc) = dfft%nr1x * dfft%nr2x * dfft%nr3p(proc+1) 
     if (proc > 0) displs(proc) = displs(proc-1) + sendcount(proc-1)
  ENDDO
  info = size( f_in ) - displs( dfft%nproc3 - 1 ) - sendcount( dfft%nproc3 - 1 )
  IF( info < 0 ) CALL fftx_error__( ' scatter_real_grid ', ' f_in too small ', -info )
  info = 0
  !write (6,*) 'scatter grid ok 1'
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_aux, sendcount(dfft%mype3), MPI_DOUBLE_PRECISION, &
                     dfft%root, dfft%comm3, info )
  !write (6,*) 'scatter grid ok 2'
  CALL fftx_error__( 'scatter_real_grid', 'info<>0', info )

  ! 2) scatter within the comm2 communicator
  IF( size( f_out ) < dfft%nnr ) &
     CALL fftx_error__( ' scatter_real_grid ', ' f_out too small ', dfft%nnr-size( f_out ) )
  !
  displs = 0 ; f_out = 0.0D0
  DO proc =0, ( dfft%nproc2 -1 )
     sendcount(proc) = dfft%nr1x * dfft%nr2p(proc+1) 
     if (proc > 0) displs(proc) = displs(proc-1) + sendcount(proc-1)
  end do
  offset_in = 1 ; offset_aux = 1
  do ir3 = 1, dfft%my_nr3p
     info = 0
  !write (6,*) 'scatter grid ok 3, ir3=', ir3
     CALL MPI_SCATTERV( f_aux(offset_aux), sendcount, displs, MPI_DOUBLE_PRECISION,   &
                        f_out(offset_in), sendcount(dfft%mype2), MPI_DOUBLE_PRECISION, &
                        dfft%root, dfft%comm2, info )
  !write (6,*) 'scatter grid ok 4, ir3=', ir3
     CALL fftx_error__( 'scatter_real_grid', 'info<>0', info )
     offset_in  = offset_in + dfft%nr1x * dfft%my_nr2p
     offset_aux = offset_aux + dfft%nr1x * dfft%nr2
  end do
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'rscatter_grid' )
  !
#else
  CALL fftx_error__('scatter_real_grid', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE scatter_real_grid
!----------------------------------------------------------------------------
SUBROUTINE scatter_complex_grid ( dfft, f_in, f_out )
  !----------------------------------------------------------------------------
  !
  ! ... scatters a real-space FFT grid from dfft%root, first processor of
  ! ... input descriptor dfft, to all others - opposite of "gather_grid"
  !
  ! ... COMPLEX*16  f_in  = gathered variable (dfft%nr1x*dfft%nr2x*dfft%nr3x)
  ! ... COMPLEX*16  f_out = distributed variable (dfft%nnr)
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(in) :: f_in (:)
  COMPLEX(DP), INTENT(inout):: f_out(:)
  TYPE ( fft_type_descriptor ), INTENT(IN) :: dfft
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:) 
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info, offset_in, offset_aux, ir3
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfft%nproc-1), sendcount(0:dfft%nproc-1)
  !
  CALL start_clock( 'cscatter_grid' )
  !
  !write (6,*) 'scatter grid ok 0'
  ALLOCATE ( f_aux(dfft%nr1x * dfft%nr2x * dfft%my_nr3p ) )
  ! 1) scatter within the comm3 communicator 
  displs = 0
  DO proc = 0, ( dfft%nproc3 - 1 )
     sendcount(proc) = 2 * dfft%nr1x * dfft%nr2x * dfft%nr3p(proc+1) 
     if (proc > 0) displs(proc) = displs(proc-1) + sendcount(proc-1)
  ENDDO
  !
  !write(*,*) 'cscatter_grid 2*size(f_in) ', 2*size(f_in); FLUSH(6)
  !write(*,*) 'cscatter_grid displ+send ', dfft%nproc3, displs(dfft%nproc3-1) + sendcount(dfft%nproc3-1); FLUSH(6)
  info = 2*size( f_in ) - displs( dfft%nproc3 - 1 ) - sendcount( dfft%nproc3 - 1 )
  IF( info < 0 ) &
     CALL fftx_error__( ' scatter_complex_grid ', ' f_in too small ', -info )
  !
  info = 0
  !write (6,*) 'scatter grid ok 1'
  CALL MPI_SCATTERV( f_in, sendcount, displs, MPI_DOUBLE_PRECISION,   &
                     f_aux, sendcount(dfft%mype3), MPI_DOUBLE_PRECISION, &
                     dfft%root, dfft%comm3, info )
  !write (6,*) 'scatter grid ok 2'
  CALL fftx_error__( ' scatter_complex_grid', 'info<>0', info )

  ! 2) scatter within the comm2 communicator
  !write(*,*) 'cscatter_grid size(f_out), dfft%nnr ', size(f_out),dfft%nnr; FLUSH(6)
  IF( size( f_out ) < dfft%nnr ) &
     CALL fftx_error__( ' scatter_complex_grid ', ' f_out too small ', dfft%nnr-size( f_out ) )
  !
  displs = 0 ; f_out = 0.0D0
  DO proc =0, ( dfft%nproc2 -1 )
     sendcount(proc) = 2 * dfft%nr1x * dfft%nr2p(proc+1) 
     if (proc > 0) displs(proc) = displs(proc-1) + sendcount(proc-1)
  end do
  offset_in = 1 ; offset_aux = 1
  do ir3 = 1, dfft%my_nr3p
     info = 0
  !write (6,*) 'scatter grid ok 3, ir3=', ir3
     CALL MPI_SCATTERV( f_aux(offset_aux), sendcount, displs, MPI_DOUBLE_PRECISION,   &
                        f_out(offset_in), sendcount(dfft%mype2), MPI_DOUBLE_PRECISION, &
                        dfft%root, dfft%comm2, info )
  !write (6,*) 'scatter grid ok 4, ir3=', ir3
     CALL fftx_error__( 'scatter_complex_grid', 'info<>0', info )
     offset_in  = offset_in + dfft%nr1x * dfft%my_nr2p
     offset_aux = offset_aux + dfft%nr1x * dfft%nr2x
  end do
  !
  ! ... the following check should be performed only on processor dfft%root
  ! ... otherwise f_in must be allocated on all processors even if not used
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'cscatter_grid' )
  !
#else
  CALL fftx_error__('scatter_complex_grid', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE scatter_complex_grid
!
! ... "gather"-like subroutines
!
!-----------------------------------------------------------------------
SUBROUTINE cgather_sym( dfftp, f_in, f_out )
  !-----------------------------------------------------------------------
  !
  ! ... gather complex data for symmetrization (used in phonon code)
  ! ... Differs from gather_grid because mpi_allgatherv is used instead
  ! ... of mpi_gatherv - all data is gathered on ALL processors
  ! ... COMPLEX*16  f_in  = distributed variable (nrxx)
  ! ... COMPLEX*16  f_out = gathered variable (nr1x*nr2x*nr3x)
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfftp
  COMPLEX(DP) :: f_in( : ), f_out(:)
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:) 
  !
#if defined(__MPI)
  !
  INTEGER :: proc, info, offset_in, offset_aux, ir3
  ! ... the following are automatic arrays
  INTEGER :: displs(0:dfftp%nproc-1), recvcount(0:dfftp%nproc-1)
  !
  CALL start_clock( 'cgather' )

  ALLOCATE ( f_aux(dfftp%nr1x * dfftp%nr2x * dfftp%my_nr3p ) )
  !
  CALL MPI_BARRIER( dfftp%comm, info )
  !
  ! 1) gather within the comm2 communicator
  displs = 0
  DO proc =0, ( dfftp%nproc2 -1 )
     recvcount(proc) = 2 * dfftp%nr1x * dfftp%nr2p(proc+1) 
     if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
  end do
  offset_in = 1; offset_aux = 1
  do ir3 = 1, dfftp%my_nr3p
     info = 0
     CALL MPI_ALLGATHERV( f_in(offset_in), recvcount(dfftp%mype2), MPI_DOUBLE_PRECISION, &
                          f_aux(offset_aux),recvcount, displs, MPI_DOUBLE_PRECISION, &
                          dfftp%comm2, info )
     CALL fftx_error__( 'cgather_sym', 'info<>0', info )
     offset_in  = offset_in  + dfftp%nr1x * dfftp%my_nr2p
     offset_aux = offset_aux + dfftp%nr1x * dfftp%nr2x
  end do

  ! 2) gather within the comm3 communicator 
  displs = 0
  DO proc = 0, ( dfftp%nproc3 - 1 )
     recvcount(proc) = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3p(proc+1) 
     if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
  ENDDO
  info = 0
  CALL MPI_ALLGATHERV( f_aux, recvcount(dfftp%mype3), MPI_DOUBLE_PRECISION, &
                       f_out, recvcount, displs, MPI_DOUBLE_PRECISION, &
                       dfftp%comm3, info )
  CALL fftx_error__( 'cgather_sym', 'info<>0', info )
  !
  ! ... the following check should be performed only on processor dfft%root
  ! ... otherwise f_out must be allocated on all processors even if not used
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'cgather' )
  !
#else
  CALL fftx_error__('cgather_sym', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE cgather_sym
!
!
!-----------------------------------------------------------------------
SUBROUTINE cgather_sym_many( dfftp, f_in, f_out, nbnd, nbnd_proc, start_nbnd_proc )
  !-----------------------------------------------------------------------
  !
  ! ... Written by A. Dal Corso
  !
  ! ... This routine generalizes cgather_sym, receiveng nbnd complex 
  ! ... distributed functions and collecting nbnd_proc(dfftp%mype+1) 
  ! ... functions in each processor.
  ! ... start_nbnd_proc(dfftp%mype+1), says where the data for each processor
  ! ... start in the distributed variable
  ! ... COMPLEX*16  f_in  = distributed variable (nrxx,nbnd)
  ! ... COMPLEX*16  f_out = gathered variable (nr1x*nr2x*nr3x,nbnd_proc(dfftp%mype+1))
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfftp
  INTEGER :: nbnd, nbnd_proc(dfftp%nproc), start_nbnd_proc(dfftp%nproc)
  COMPLEX(DP) :: f_in(dfftp%nnr,nbnd)
  COMPLEX(DP) :: f_out(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,nbnd_proc(dfftp%mype+1))
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:) 
  !
#if defined(__MPI)
  !
  INTEGER :: proc_, proc2_, proc3_, proc, info, offset_in, offset_aux, ir3, nr2px
  INTEGER :: ibnd, jbnd, iii
  INTEGER :: displs(0:dfftp%nproc-1), recvcount(0:dfftp%nproc-1)
  !
  !
  CALL start_clock( 'cgather' )
  !
  ALLOCATE ( f_aux(dfftp%nr1x * dfftp%nr2x * dfftp%my_nr3p ) )
  !
  CALL MPI_BARRIER( dfftp%comm, info )
  !
  f_out = (0.d0,0.d0)

  !write (*,*) 'enter cgather ', nbnd
  !write (*,*) nbnd_proc
  !write (*,*) start_nbnd_proc
  !write (*,*) 'dfftp%nproc',dfftp%nproc 
  !do ibnd =1,3
  !write (*,*) 'evc = ',ibnd
  !write (*,*) f_in(1:3,ibnd)
  !end do


  nr2px = MAXVAL ( dfftp%nr2p )  ! maximum number of Y values to be disributed

  DO proc_ = 0, dfftp%nproc - 1
     ! define the processor index of the two sub-communicators
     proc2_ = dfftp%iproc2(proc_ + 1) -1 ; proc3_ = dfftp%iproc3(proc_ + 1) -1
     !write (*,*) ' proc_, proc2_, proc3_, dfftp%mype2', proc_, proc2_, proc3_, dfftp%mype2
     DO ibnd = 1, nbnd_proc(proc_ +1)
        jbnd = start_nbnd_proc(proc_ +1) + ibnd - 1
        !write (*,*) ' proc_, ibnd, jbnd ', proc_, ibnd, jbnd
        f_aux(:) = (0.d0,0.d0)
        ! 1) gather within the comm2 communicator
        displs = 0
        DO proc =0, ( dfftp%nproc2 -1 )
           recvcount(proc) = 2 * dfftp%nr1x * dfftp%nr2p(proc+1)
           if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
        end do
        !write (*,*) 'dfftp%nr1x ', dfftp%nr1x
        !write (*,*) 'dfftp%nr2x ', dfftp%nr2x
        !write (*,*) 'dfftp%nr3x ', dfftp%nr3x
        !write (*,*) 'dfftp%nr2p ', dfftp%nr2p, ' nnr3px ', nr2px
        !write (*,*) 'recvcount ', recvcount(0:dfftp%nproc2 -1)
        !write (*,*) 'displs ', displs(0:dfftp%nproc2 -1)
        offset_in = 1; offset_aux = 1
        do ir3 = 1, dfftp%my_nr3p
           info = 0
           CALL MPI_GATHERV( f_in(offset_in,jbnd), recvcount(dfftp%mype2), MPI_DOUBLE_PRECISION, &
                             f_aux(offset_aux), recvcount, displs, MPI_DOUBLE_PRECISION, &
                             proc2_, dfftp%comm2, info )
           CALL fftx_error__( 'cgather_sym_many', 'info<>0', info )
           offset_in  = offset_in  + dfftp%nr1x * dfftp%my_nr2p
           offset_aux = offset_aux + dfftp%nr1x * dfftp%nr2x
        end do
        !write(*,*) ' -> f_in(...+1:...+3,jbnd) ',jbnd
        !offset_in = 0
        !do iii =1,2*dfftp%nr2x
        !   write(*,'(i4,3("(",2f10.7,") "))') iii, f_in(offset_in+1:offset_in+3, jbnd)
        !   offset_in  = offset_in  + dfftp%nr1x 
        !enddo
        IF (dfftp%mype2==proc2_) THEN
        !   write(*,*) ' -> f_aux(...+1:...+3)'
        !   offset_aux = 0
        !   do iii =1,4*dfftp%nr2x
        !      write(*,'(i4,3("(",2f10.7,") "))') iii, f_aux(offset_aux+1:offset_aux+3)
        !      offset_aux = offset_aux + dfftp%nr1x
        !   enddo
        !   write(*,*) ' -> F_AUX(...+1:...+3)'
        !   offset_aux = 0
        !   do iii =1,dfftp%my_nr3p
        !      write(*,'(i4,3("(",2f10.7,") "))') iii, f_aux(offset_aux+1:offset_aux+3)
        !      offset_aux  = offset_aux  + dfftp%nr1x*dfftp%nr2x
        !   enddo
           ! 2) gather within the comm3 communicator 
           displs = 0
           DO proc = 0, ( dfftp%nproc3 - 1 )
              recvcount(proc) = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3p(proc+1) 
              if (proc > 0) displs(proc) = displs(proc-1) + recvcount(proc-1)
           ENDDO
           info = 0
           CALL MPI_GATHERV( f_aux, recvcount(dfftp%mype3), MPI_DOUBLE_PRECISION, &
                             f_out(1,ibnd), recvcount, displs, MPI_DOUBLE_PRECISION, &
                             proc3_, dfftp%comm3, info )
           CALL fftx_error__( 'cgather_sym_many', 'info<>0', info )
   
        !   IF (dfftp%mype3==proc3_) THEN
        !      write(*,*) ' -> f_out(...+1:...+3,ibnd) ',ibnd
        !      offset_in  = 0
        !      do iii =1,dfftp%nr3x
        !         write(*,'(i4,3("(",2f10.7,") "))') iii, f_out(offset_in+1:offset_in+3, ibnd)
        !         offset_in  = offset_in  + dfftp%nr1x*dfftp%nr2x
        !      enddo
        !   END IF
        END IF

     END DO
  END DO
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'cgather' )
  !
#else
  CALL fftx_error__('cgather_sym_many', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE cgather_sym_many
!
!----------------------------------------------------------------------------
SUBROUTINE cscatter_sym_many( dfftp, f_in, f_out, target_ibnd, nbnd, nbnd_proc, &
                              start_nbnd_proc   )
  !----------------------------------------------------------------------------
  !
  ! ... Written by A. Dal Corso
  !
  ! ... generalizes cscatter_sym. It assumes that each processor has
  ! ... a certain number of bands (nbnd_proc(dfftp%mype+1)). The processor 
  ! ... that has target_ibnd scatters it to all the other processors 
  ! ... that receive a distributed part of the target function. 
  ! ... start_nbnd_proc(dfftp%mype+1) is used to identify the processor
  ! ... that has the required band
  !
  ! ... COMPLEX*16  f_in  = gathered variable (nr1x*nr2x*nr3x, nbnd_proc(dfftp%mype+1) )
  ! ... COMPLEX*16  f_out = distributed variable (nrxx)
  !
  IMPLICIT NONE
  !
  TYPE (fft_type_descriptor), INTENT(in) :: dfftp
  INTEGER :: nbnd, nbnd_proc(dfftp%nproc), start_nbnd_proc(dfftp%nproc)
  COMPLEX(DP) :: f_in(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x,nbnd_proc(dfftp%mype+1))
  COMPLEX(DP) :: f_out(dfftp%nnr)
  COMPLEX(DP), ALLOCATABLE ::  f_aux(:) 
  INTEGER :: target_ibnd
  !
#if defined(__MPI)
  !
  INTEGER :: proc_, proc2_, proc3_, proc, info, offset_out, offset_aux, ir3
  INTEGER :: displs(0:dfftp%nproc-1), sendcount(0:dfftp%nproc-1)
  INTEGER :: ibnd, jbnd
  !
  CALL start_clock( 'cscatter_sym' )
  !
  ALLOCATE ( f_aux(dfftp%nr1x * dfftp%nr2x * dfftp%my_nr3p ) )
  !
  f_out = (0.0_DP, 0.0_DP)
  !
  CALL MPI_BARRIER( dfftp%comm, info )
  !
  DO proc_ = 0, dfftp%nproc - 1
     ! define the processor index of the two sub-communicators
     proc2_ = dfftp%iproc2(proc_ + 1) -1 ; proc3_ = dfftp%iproc3(proc_ + 1) -1
     DO ibnd = 1, nbnd_proc(proc_+1)
        jbnd = start_nbnd_proc(proc_+1) + ibnd - 1
        IF (jbnd/=target_ibnd) CYCLE
        IF (dfftp%mype2==proc2_) THEN
           ! 1) scatter within the comm3 communicator 
           displs = 0
           DO proc = 0, ( dfftp%nproc3 - 1 )
              sendcount(proc) = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3p(proc+1) 
              if (proc > 0) displs(proc) = displs(proc-1) + sendcount(proc-1)
           ENDDO
           info = 0
           CALL MPI_SCATTERV( f_in(1,ibnd), sendcount, displs, MPI_DOUBLE_PRECISION,   &
                              f_aux, sendcount(dfftp%mype3), MPI_DOUBLE_PRECISION, &
                              proc3_, dfftp%comm3, info )
        ELSE
           f_aux=(0.d0,0.d0)
        END IF
        ! 2) scatter within the comm2 communicator
        displs = 0 ; f_out = 0.0D0
        DO proc =0, ( dfftp%nproc2 -1 )
           sendcount(proc) = 2 * dfftp%nr1x * dfftp%nr2p(proc+1) 
           if (proc > 0) displs(proc) = displs(proc-1) + sendcount(proc-1)
        end do
        offset_out = 1 ; offset_aux = 1
        do ir3 = 1, dfftp%my_nr3p
           info = 0
           CALL MPI_SCATTERV( f_aux(offset_aux), sendcount, displs, MPI_DOUBLE_PRECISION,   &
                              f_out(offset_out), sendcount(dfftp%mype2), MPI_DOUBLE_PRECISION, &
                              proc2_, dfftp%comm2, info )
           CALL fftx_error__( 'gather_grid', 'info<>0', info )
           offset_out = offset_out + dfftp%nr1x * dfftp%my_nr2p
           offset_aux = offset_aux + dfftp%nr1x * dfftp%nr2x
        end do
     ENDDO
  ENDDO
  !
  DEALLOCATE ( f_aux )
  !
  CALL stop_clock( 'cscatter_sym' )
  !
#else
  CALL fftx_error__('cscatter_sym_many', 'do not use in serial execution', 1)
#endif
  !
  RETURN
  !
END SUBROUTINE cscatter_sym_many


!=----------------------------------------------------------------------=!
   END MODULE scatter_mod
!=----------------------------------------------------------------------=!
!
!
!---------------------------------------------------------------------
subroutine fftsort (n, ia)  
  !---------------------------------------------------------------------
  ! sort an integer array ia(1:n) into ascending order using heapsort algorithm.
  ! n is input, ia is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ia).
  ! in case of equal values in the data array (ia) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  implicit none  
  !-input/output variables
  integer :: n  
  integer :: ia (2,n)  
  !-local variables
  integer :: i, ir, j, l
  integer :: iia(2)  
  ! nothing to order
  if (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  
  ir = n  
10 continue  
  ! still in hiring phase
  if (l.gt.1) then  
     l = l - 1  
     iia(:) = ia (:,l)  
     ! in retirement-promotion phase.
  else  
     ! clear a space at the end of the array
     iia(:) = ia (:,ir)  
     !
     ! retire the top of the heap into it
     ia (:,ir) = ia (:,1)  
     !
     ! decrease the size of the corporation
     ir = ir - 1  
     ! done with the last promotion
     if (ir.eq.1) then  
        ! the least competent worker at all !
        ia (:,1) = iia(:)  
        !
        return  
     endif
  endif
  ! wheter in hiring or promotion phase, we
  i = l  
  ! set up to place iia in its proper level
  j = l + l  
  !
  do while (j.le.ir)  
     if (j.lt.ir) then  
        if (ia (1,j) .lt. ia (1,j + 1) ) then  
           j = j + 1  
        endif
     endif
     ! demote iia
     if (iia(1).lt.ia (1,j) ) then  
        ia (:,i) = ia (:,j)  
        i = j  
        j = j + j  
     else  
        ! set j to terminate do-while loop
        j = ir + 1  
     endif
  enddo
  ia (:,i) = iia(:)  
  goto 10  
  !
end subroutine fftsort

