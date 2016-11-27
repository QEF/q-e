!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if defined(__DFTI)
#include "mkl_dfti.f90"
!=----------------------------------------------------------------------=!
   MODULE fft_scalar_dfti
!=----------------------------------------------------------------------=!

       USE MKL_DFTI ! -- this can be found in the MKL include directory
       USE fft_param

       IMPLICIT NONE
        SAVE

        PRIVATE
        PUBLIC :: cft_1z, cft_2xy, cfft3d, cfft3ds

        TYPE dfti_descriptor_array
           TYPE(DFTI_DESCRIPTOR), POINTER :: desc
        END TYPE

!=----------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------=!

!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "z"
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cft_1z(c, nsl, nz, ldz, isign, cout)

!     driver routine for nsl 1d complex fft's of length nz
!     ldz >= nz is the distance between sequences to be transformed
!     (ldz>nz is used on some architectures to reduce memory conflicts)
!     input  :  c(ldz*nsl)   (complex)
!     output : cout(ldz*nsl) (complex - NOTA BENE: transform is not in-place!)
!     isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nz, nsl, ldz) are stored and re-used if available

     INTEGER, INTENT(IN) :: isign
     INTEGER, INTENT(IN) :: nsl, nz, ldz

     COMPLEX (DP) :: c(:), cout(:)

     REAL (DP)  :: tscale
     INTEGER    :: i, err, idir, ip, void
     INTEGER, SAVE :: zdims( 3, ndims ) = -1
     INTEGER, SAVE :: icurrent = 1
     LOGICAL :: found

     INTEGER :: tid

#if defined(__OPENMP)
     INTEGER :: offset, ldz_t
     INTEGER :: omp_get_max_threads
     EXTERNAL :: omp_get_max_threads
#endif

     !   Intel MKL native FFT driver

     TYPE(DFTI_DESCRIPTOR_ARRAY), SAVE :: hand( ndims )
     LOGICAL, SAVE :: dfti_first = .TRUE.
     INTEGER :: dfti_status = 0
     !
     CALL check_dims()
     !
     !   Here initialize table only if necessary
     !
     CALL lookup()

     IF( .NOT. found ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       CALL init_dfti()

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'cft_1z' )
#endif

     IF (isign < 0) THEN
        dfti_status = DftiComputeForward(hand(ip)%desc, c, cout )
        IF(dfti_status /= 0) &
           CALL fftx_error__(' cft_1z ',' stopped in DftiComputeForward ', dfti_status )
     ELSE IF (isign > 0) THEN
        dfti_status = DftiComputeBackward(hand(ip)%desc, c, cout )
        IF(dfti_status /= 0) &
           CALL fftx_error__(' cft_1z ',' stopped in DftiComputeBackward ', dfti_status )
     END IF

#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'cft_1z' )
#endif

     RETURN

   CONTAINS !=------------------------------------------------=!

     SUBROUTINE check_dims()
     IF( nsl < 0 ) THEN
       CALL fftx_error__(" fft_scalar: cft_1z ", " nsl out of range ", nsl)
     END IF
     END SUBROUTINE check_dims

     SUBROUTINE lookup()
     IF( dfti_first ) THEN
        DO ip = 1, ndims
           hand(ip)%desc => NULL()
        END DO
        dfti_first = .FALSE.
     END IF
     DO ip = 1, ndims
        !   first check if there is already a table initialized
        !   for this combination of parameters
        found = ( nz == zdims(1,ip) )
        !   The initialization in ESSL and FFTW v.3 depends on all three parameters
        found = found .AND. ( nsl == zdims(2,ip) ) .AND. ( ldz == zdims(3,ip) )
        IF (found) EXIT
     END DO
     END SUBROUTINE lookup

     SUBROUTINE init_dfti()

       if( ASSOCIATED( hand( icurrent )%desc ) ) THEN
          dfti_status = DftiFreeDescriptor( hand( icurrent )%desc )
          IF( dfti_status /= 0) THEN
             WRITE(*,*) "stopped in DftiFreeDescriptor", dfti_status
             STOP
          ENDIF
       END IF

       dfti_status = DftiCreateDescriptor(hand( icurrent )%desc, DFTI_DOUBLE, DFTI_COMPLEX, 1,nz)
       IF(dfti_status /= 0)  &
         CALL fftx_error__(' cft_1z ',' stopped in DftiCreateDescriptor ', dfti_status )

       dfti_status = DftiSetValue(hand( icurrent )%desc, DFTI_NUMBER_OF_TRANSFORMS,nsl)
       IF(dfti_status /= 0) &
         CALL fftx_error__(' cft_1z ',' stopped in DFTI_NUMBER_OF_TRANSFORMS ', dfti_status )

       dfti_status = DftiSetValue(hand( icurrent )%desc,DFTI_INPUT_DISTANCE, ldz )
       IF(dfti_status /= 0) &
         CALL fftx_error__(' cft_1z ',' stopped in DFTI_INPUT_DISTANCE ', dfti_status )

       dfti_status = DftiSetValue(hand( icurrent )%desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
       IF(dfti_status /= 0) &
         CALL fftx_error__(' cft_1z ',' stopped in DFTI_PLACEMENT ', dfti_status )

       dfti_status = DftiSetValue(hand( icurrent )%desc,DFTI_OUTPUT_DISTANCE, ldz )
       IF(dfti_status /= 0) &
         CALL fftx_error__(' cft_1z ',' stopped in DFTI_OUTPUT_DISTANCE ', dfti_status )

       tscale = 1.0_DP/nz
       dfti_status = DftiSetValue( hand( icurrent )%desc, DFTI_FORWARD_SCALE, tscale);
       IF(dfti_status /= 0) &
         CALL fftx_error__(' cft_1z ',' stopped in DFTI_FORWARD_SCALE ', dfti_status )

       dfti_status = DftiSetValue( hand( icurrent )%desc, DFTI_BACKWARD_SCALE, DBLE(1) );
       IF(dfti_status /= 0) &
         CALL fftx_error__(' cft_1z ',' stopped in DFTI_BACKWARD_SCALE ', dfti_status )

       dfti_status = DftiCommitDescriptor(hand( icurrent )%desc)
       IF(dfti_status /= 0) &
         CALL fftx_error__(' cft_1z ',' stopped in DftiCommitDescriptor ', dfti_status )

       zdims(1,icurrent) = nz; zdims(2,icurrent) = nsl; zdims(3,icurrent) = ldz;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END SUBROUTINE init_dfti


   END SUBROUTINE cft_1z

!
!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "x" and "y" direction
!
!
!
!=----------------------------------------------------------------------=!
!
!

   SUBROUTINE cft_2xy(r, nzl, nx, ny, ldx, ldy, isign, pl2ix)

!     driver routine for nzl 2d complex fft's of lengths nx and ny
!     input : r(ldx*ldy)  complex, transform is in-place
!     ldx >= nx, ldy >= ny are the physical dimensions of the equivalent
!     2d array: r2d(ldx, ldy) (x first dimension, y second dimension)
!     (ldx>nx, ldy>ny used on some architectures to reduce memory conflicts)
!     pl2ix(nx) (optional) is 1 for columns along y to be transformed
!     isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nx,ny,nzl,ldx) are stored and re-used if available

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: isign, ldx, ldy, nx, ny, nzl
     INTEGER, OPTIONAL, INTENT(IN) :: pl2ix(:)
     COMPLEX (DP) :: r( : )
     INTEGER :: i, k, j, err, idir, ip, kk, void
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims( 4, ndims) = -1
     LOGICAL :: dofft( nfftx ), found
     INTEGER, PARAMETER  :: stdout = 6

#if defined(__OPENMP)
     INTEGER :: offset
     INTEGER :: nx_t, ny_t, nzl_t, ldx_t, ldy_t
     INTEGER  :: itid, mytid, ntids
     INTEGER  :: omp_get_thread_num, omp_get_num_threads,omp_get_max_threads
     EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

     TYPE(DFTI_DESCRIPTOR_ARRAY), SAVE :: hand( ndims )
     LOGICAL, SAVE :: dfti_first = .TRUE.
     INTEGER :: dfti_status = 0

     dofft( 1 : nx ) = .TRUE.
     IF( PRESENT( pl2ix ) ) THEN
       IF( SIZE( pl2ix ) < nx ) &
         CALL fftx_error__( ' cft_2xy ', ' wrong dimension for arg no. 8 ', 1 )
       DO i = 1, nx
         IF( pl2ix(i) < 1 ) dofft( i ) = .FALSE.
       END DO
     END IF

     !
     !   Here initialize table only if necessary
     !

     CALL lookup()

     IF( .NOT. found ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       CALL init_dfti()

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'cft_2xy' )
#endif

     IF( isign < 0 ) THEN
        !
        dfti_status = DftiComputeForward(hand(ip)%desc, r(:))
        IF(dfti_status /= 0)THEN
           WRITE(*,*) "stopped in DftiComputeForward", dfti_status
           STOP
        ENDIF
        !
     ELSE IF( isign > 0 ) THEN
        !
        dfti_status = DftiComputeBackward(hand(ip)%desc, r(:))
        IF(dfti_status /= 0)THEN
           WRITE(*,*) "stopped in DftiComputeBackward", dfti_status
           STOP
        ENDIF
        !
     END IF


#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'cft_2xy' )
#endif

     RETURN

   CONTAINS !=------------------------------------------------=!

     SUBROUTINE check_dims()
     END SUBROUTINE check_dims

     SUBROUTINE lookup()
     IF( dfti_first ) THEN
        DO ip = 1, ndims
           hand(ip)%desc => NULL()
        END DO
        dfti_first = .FALSE.
     END IF
     DO ip = 1, ndims
       !   first check if there is already a table initialized
       !   for this combination of parameters
       found = ( ny == dims(1,ip) ) .AND. ( nx == dims(3,ip) )
       found = found .AND. ( ldx == dims(2,ip) ) .AND.  ( nzl == dims(4,ip) )
       IF (found) EXIT
     END DO
     END SUBROUTINE lookup

     SUBROUTINE init_dfti()

       if( ASSOCIATED( hand( icurrent )%desc ) ) THEN
          dfti_status = DftiFreeDescriptor( hand( icurrent )%desc )
          IF( dfti_status /= 0) THEN
             WRITE(*,*) "stopped in DftiFreeDescriptor", dfti_status
             STOP
          ENDIF
       END IF

       dfti_status = DftiCreateDescriptor(hand( icurrent )%desc, DFTI_DOUBLE, DFTI_COMPLEX, 2,(/nx,ny/))
       IF(dfti_status /= 0) THEN
          WRITE(*,*) "stopped in DftiCreateDescriptor", dfti_status
          STOP
       ENDIF
       dfti_status = DftiSetValue(hand( icurrent )%desc, DFTI_NUMBER_OF_TRANSFORMS,nzl)
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DFTI_NUMBER_OF_TRANSFORMS", dfti_status
          STOP
       ENDIF
       dfti_status = DftiSetValue(hand( icurrent )%desc,DFTI_INPUT_DISTANCE, ldx*ldy )
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DFTI_INPUT_DISTANCE", dfti_status
          STOP
       ENDIF
       dfti_status = DftiSetValue(hand( icurrent )%desc, DFTI_PLACEMENT, DFTI_INPLACE)
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DFTI_PLACEMENT", dfti_status
          STOP
       ENDIF
       tscale = 1.0_DP/ (nx * ny )
       dfti_status = DftiSetValue( hand( icurrent )%desc, DFTI_FORWARD_SCALE, tscale);
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DFTI_FORWARD_SCALE", dfti_status
          STOP
       ENDIF
       dfti_status = DftiSetValue( hand( icurrent )%desc, DFTI_BACKWARD_SCALE, DBLE(1) );
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DFTI_BACKWARD_SCALE", dfti_status
          STOP
       ENDIF
       dfti_status = DftiCommitDescriptor(hand( icurrent )%desc)
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in DftiCommitDescriptor", dfti_status
          STOP
       ENDIF

       dims(1,icurrent) = ny; dims(2,icurrent) = ldx;
       dims(3,icurrent) = nx; dims(4,icurrent) = nzl;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_dfti

   END SUBROUTINE cft_2xy


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
   SUBROUTINE cfft3d( f, nx, ny, nz, ldx, ldy, ldz, howmany, isign )

  !     driver routine for 3d complex fft of lengths nx, ny, nz
  !     input  :  f(ldx*ldy*ldz)  complex, transform is in-place
  !     ldx >= nx, ldy >= ny, ldz >= nz are the physical dimensions
  !     of the equivalent 3d array: f3d(ldx,ldy,ldz)
  !     (ldx>nx, ldy>ny, ldz>nz may be used on some architectures
  !      to reduce memory conflicts - not implemented for FFTW)
  !     isign > 0 : f(G) => f(R)   ; isign < 0 : f(R) => f(G)
  !
  !     howmany: perform this many ffts, separated by ldx*ldy*ldz in memory
  !     Up to "ndims" initializations (for different combinations of input
  !     parameters nx,ny,nz) are stored and re-used if available

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
     COMPLEX (DP) :: f(:)
     INTEGER :: i, k, j, err, idir, ip
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims(4,ndims) = -1

     !   Intel MKL native FFT driver

     TYPE(DFTI_DESCRIPTOR_ARRAY), SAVE :: hand(ndims)
     LOGICAL, SAVE :: dfti_first = .TRUE.
     INTEGER :: dfti_status = 0
     !

     CALL check_dims()

     !
     !   Here initialize table only if necessary
     !

     CALL lookup()

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       CALL init_dfti()

     END IF

     !
     !   Now perform the 3D FFT using the machine specific driver
     !

     IF( isign < 0 ) THEN
        !
        dfti_status = DftiComputeForward(hand(ip)%desc, f(1:))
        IF(dfti_status /= 0)THEN
           WRITE(*,*) "stopped in cfft3d, DftiComputeForward", dfti_status
           STOP
        ENDIF
        !
     ELSE IF( isign > 0 ) THEN
        !
        dfti_status = DftiComputeBackward(hand(ip)%desc, f(1:))
        IF(dfti_status /= 0)THEN
           WRITE(*,*) "stopped in cfft3d, DftiComputeBackward", dfti_status
           STOP
        ENDIF
        !
     END IF

     RETURN

   CONTAINS !=------------------------------------------------=!

     SUBROUTINE check_dims()
     IF ( nx < 1 ) &
         call fftx_error__('cfft3d',' nx is less than 1 ', 1)
     IF ( ny < 1 ) &
         call fftx_error__('cfft3d',' ny is less than 1 ', 1)
     IF ( nz < 1 ) &
         call fftx_error__('cfft3d',' nz is less than 1 ', 1)
     IF ( howmany < 1 ) &
         call fftx_error__('cfft3d',' howmany is less than 1 ', 1)
     END SUBROUTINE check_dims

     SUBROUTINE lookup()
     IF( dfti_first ) THEN
        DO ip = 1, ndims
           hand(ip)%desc => NULL()
        END DO
        dfti_first = .FALSE.
     END IF
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

     SUBROUTINE init_dfti()
      if( ASSOCIATED( hand(icurrent)%desc ) ) THEN
          dfti_status = DftiFreeDescriptor( hand(icurrent)%desc )
          IF( dfti_status /= 0) THEN
             WRITE(*,*) "stopped in cfft3d, DftiFreeDescriptor", dfti_status
             STOP
          ENDIF
       END IF

       dfti_status = DftiCreateDescriptor(hand(icurrent)%desc, DFTI_DOUBLE, DFTI_COMPLEX, 3,(/nx,ny,nz/))
       IF(dfti_status /= 0) THEN
          WRITE(*,*) "stopped in cfft3d, DftiCreateDescriptor", dfti_status
          STOP
       ENDIF
       dfti_status = DftiSetValue(hand(icurrent)%desc, DFTI_NUMBER_OF_TRANSFORMS,howmany)
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in cfft3d, DFTI_NUMBER_OF_TRANSFORMS", dfti_status
          STOP
       ENDIF
       dfti_status = DftiSetValue(hand(icurrent)%desc, DFTI_INPUT_DISTANCE, ldx*ldy*ldz)
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in cfft3dm, DFTI_INPUT_DISTANCE", dfti_status
          STOP
       ENDIF
       dfti_status = DftiSetValue(hand(icurrent)%desc, DFTI_PLACEMENT, DFTI_INPLACE)
       IF(dfti_status /= 0)THEN
         WRITE(*,*) "stopped in cfft3d, DFTI_PLACEMENT", dfti_status
         STOP
       ENDIF
       tscale = 1.0_DP/ (nx * ny * nz)
       dfti_status = DftiSetValue( hand(icurrent)%desc, DFTI_FORWARD_SCALE, tscale);
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in cfft3d, DFTI_FORWARD_SCALE", dfti_status
          STOP
       ENDIF
       tscale = 1.0_DP
       dfti_status = DftiSetValue( hand(icurrent)%desc, DFTI_BACKWARD_SCALE, tscale );
       IF(dfti_status /= 0)THEN
          WRITE(*,*) "stopped in cfft3d, DFTI_BACKWARD_SCALE", dfti_status
          STOP
       ENDIF

       dfti_status = DftiCommitDescriptor(hand(icurrent)%desc)
       IF(dfti_status /= 0) THEN
          WRITE(*,*) "stopped in cfft3d, DftiCreateDescriptor", dfti_status
          STOP
       ENDIF

       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz; dims(4,icurrent) = howmany
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_dfti

   END SUBROUTINE cfft3d

!
!=----------------------------------------------------------------------=!
!
!
!
!         3D scalar FFTs,  but using sticks!
!
!
!
!=----------------------------------------------------------------------=!
!

SUBROUTINE cfft3ds (f, nx, ny, nz, ldx, ldy, ldz, howmany, isign, do_fft_z, do_fft_y)
  !
  implicit none

  integer :: nx, ny, nz, ldx, ldy, ldz, isign, howmany
  !
  complex(DP) :: f ( ldx * ldy * ldz )
  integer :: do_fft_y(:), do_fft_z(:)
  !
  CALL cfft3d (f, nx, ny, nz, ldx, ldy, ldz, howmany, isign)

END SUBROUTINE cfft3ds
#else
   MODULE fft_scalar_dfti
#endif
!=----------------------------------------------------------------------=!
END MODULE fft_scalar_dfti
!=----------------------------------------------------------------------=!
