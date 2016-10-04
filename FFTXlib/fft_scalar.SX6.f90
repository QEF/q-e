!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------!
! FFT scalar drivers Module - contains machine-dependent routines for      !
! FFTW, FFTW3, ESSL (both 3d for serial execution and 1d+2d FFTs for       !
! parallel execution; NEC ASL libraries (3d only, no parallel execution)   !
! Written by Carlo Cavazzoni, modified by P. Giannozzi, contributions      !
! by Martin Hilgemans, Guido Roma, Pascal Thibaudeau, Stephane Lefranc,    !
! Nicolas Lacorne, Filippo Spiga, Nicola Varini - Last update Jul 2015     !
!--------------------------------------------------------------------------!

!=----------------------------------------------------------------------=!
   MODULE fft_scalar
!=----------------------------------------------------------------------=!

       USE, intrinsic ::  iso_c_binding
       
       IMPLICIT NONE
        SAVE

        PRIVATE
        PUBLIC :: cft_1z, cft_2xy, cfft3d, cfft3ds

! ...   Local Parameter

#include "fft_param.f90"

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
     LOGICAL :: done

     INTEGER :: tid

     ! ...   Machine-Dependent parameters, work arrays and tables of factors

     !   ltabl   Dimension of the tables of factors calculated at the
     !           initialization stage

#if defined(__OPENMP)
     INTEGER :: offset, ldz_t
     INTEGER :: omp_get_max_threads
     EXTERNAL :: omp_get_max_threads
#endif

     !   NEC MathKeisan

     INTEGER, PARAMETER :: ltabl = 2 * nfftx + 64
     REAL (DP), SAVE :: tablez (ltabl, ndims)
     REAL (DP)       :: work(4*nz*nsl)
     COMPLEX (DP)    :: DUMMY
     INTEGER, SAVE :: isys = 1

     IF( nsl < 0 ) THEN
       CALL fftx_error__(" fft_scalar: cft_1z ", " nsl out of range ", nsl)
     END IF

     !
     !   Here initialize table only if necessary
     !

     DO ip = 1, ndims

        !   first check if there is already a table initialized
        !   for this combination of parameters

        done = ( nz == zdims(1,ip) )
        IF (done) EXIT
     END DO

     IF( .NOT. done ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       ! WRITE( stdout, fmt="('DEBUG cft_1z, reinitializing tables ', I3)" ) icurrent

       CALL ZZFFTM (0, nz, 1, 1.0_DP, DUMMY, ldz, DUMMY, ldz, &
                    tablez (1, icurrent), work, isys)

       zdims(1,icurrent) = nz; zdims(2,icurrent) = nsl; zdims(3,icurrent) = ldz;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'cft_1z' )
#endif


     IF ( isign < 0 ) THEN
        idir   = -1
        tscale = 1.0_DP / nz
     ELSE IF ( isign > 0 ) THEN
        idir   = 1
        tscale = 1.0_DP
     END IF
     IF (isign /= 0) CALL ZZFFTM (idir, nz, nsl, tscale, c(1), ldz, &
          cout(1), ldz, tablez (1, ip), work, isys)

#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'cft_1z' )
#endif

     RETURN


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
     LOGICAL :: dofft( nfftx ), done
     INTEGER, PARAMETER  :: stdout = 6

#if defined(__OPENMP)
     INTEGER :: offset
     INTEGER :: nx_t, ny_t, nzl_t, ldx_t, ldy_t
     INTEGER  :: itid, mytid, ntids
     INTEGER  :: omp_get_thread_num, omp_get_num_threads,omp_get_max_threads
     EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif


     INTEGER, PARAMETER :: ltabl = 2*nfftx + 64
     REAL (DP), SAVE :: tablex(ltabl, ndims), tabley(ltabl, ndims)
     REAL (DP)       :: work(4*nx*ny)
     COMPLEX (DP) :: XY(ldx*ny)
     COMPLEX (DP) :: DUMMY
     INTEGER, SAVE :: isys = 1


     dofft( 1 : nx ) = .TRUE.
     IF( PRESENT( pl2ix ) ) THEN
       IF( SIZE( pl2ix ) < nx ) &
         CALL fftx_error__( ' cft_2xy ', ' wrong dimension for arg no. 8 ', 1 )
       DO i = 1, nx
         IF( pl2ix(i) < 1 ) dofft( i ) = .FALSE.
       END DO
     END IF

     ! WRITE( stdout,*) 'DEBUG: ', COUNT( dofft )

     !
     !   Here initialize table only if necessary
     !

     DO ip = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       done = ( ny == dims(1,ip) ) .AND. ( nx == dims(3,ip) )
       done = done .AND. ( ldx == dims(2,ip) ) .AND.  ( nzl == dims(4,ip) )
       IF (done) EXIT

     END DO

     IF( .NOT. done ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       ! WRITE( stdout, fmt="('DEBUG cft_2xy, reinitializing tables ', I3)" ) icurrent

       CALL ZZFFT(0, ny, 1.0_DP, DUMMY, DUMMY,              &
                  tabley (1, icurrent), work, isys)
       CALL ZZFFTM  (0, nx, 1, 1.0_DP, DUMMY, ldx, DUMMY, ldx,           &
                     tablex(1, icurrent), work, isys)

       dims(1,icurrent) = ny; dims(2,icurrent) = ldx;
       dims(3,icurrent) = nx; dims(4,icurrent) = nzl;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'cft_2xy' )
#endif


      IF( isign < 0 ) THEN

       idir = -1
       tscale = 1.0_DP / (nx * ny)
       DO k = 0, nzl-1
          kk = k * ldx * ldy
! FORWARD: ny FFTs in the X direction
          CALL ZZFFTM ( idir, nx, ny, tscale, r(kk+1), ldx, r(kk+1), ldx,   &
                        tablex (1, ip), work(1), isys )
! FORWARD: nx FFTs in the Y direction
          DO i = 1, nx
             IF ( dofft(i) ) THEN
                DO j = 0, ny-1
                   XY(j+1) = r(i + (j) * ldx + kk)
                END DO
                CALL ZZFFT(idir, ny, 1.0_DP, XY, XY, tabley (1, ip),      &
                           work(1), isys)
                DO j = 0, ny-1
                   r(i + (j) * ldx + kk) = XY(j+1)
                END DO
             END IF
          END DO
       END DO

     ELSE IF ( isign > 0 ) THEN

       idir = 1
       tscale = 1.0_DP
       DO k = 0, nzl-1
! BACKWARD: nx FFTs in the Y direction
          kk = (k) * ldx * ldy
          DO i = 1, nx
             IF ( dofft(i) ) THEN
                DO j = 0, ny-1
                   XY(j+1) = r(i + (j) * ldx + kk)
                END DO
                CALL ZZFFT(idir, ny, 1.0_DP, XY, XY, tabley (1, ip),      &
                           work(1), isys)
                DO j = 0, ny-1
                   r(i + (j) * ldx + kk) = XY(j+1)
                END DO
             END IF
          END DO
! BACKWARD: ny FFTs in the X direction
          CALL ZZFFTM ( idir, nx, ny, tscale, r(kk+1), ldx, r(kk+1), ldx,   &
                        tablex (1, ip), work(1), isys )
       END DO

     END IF

#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'cft_2xy' )
#endif

     RETURN

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
  !     Up to "ndims" initializations (for different combinations of input
  !     parameters nx,ny,nz) are stored and re-used if available

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
     COMPLEX (DP) :: f(:)
     INTEGER :: i, k, j, err, idir, ip
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims(3,ndims) = -1

     INTEGER, PARAMETER :: ltabl = 60
     INTEGER, PARAMETER :: lwork = 195+6*nfftx
     INTEGER, SAVE  :: iw0(ltabl, ndims)
     INTEGER :: k_off, kj_offset
     REAL (DP), SAVE :: auxp (lwork, ndims)
     ! not sure whether auxp is work space or not
     COMPLEX(DP), DIMENSION(:), ALLOCATABLE :: cw2
     COMPLEX (DP) :: f_out(size(f))

#if defined(ASL) && defined(MICRO)
     INTEGER :: nbtasks
     COMMON/NEC_ASL_PARA/nbtasks
#endif


     IF ( nx < 1 ) &
         call fftx_error__('cfft3d',' nx is less than 1 ', 1)
     IF ( ny < 1 ) &
         call fftx_error__('cfft3d',' ny is less than 1 ', 1)
     IF ( nz < 1 ) &
         call fftx_error__('cfft3d',' nz is less than 1 ', 1)
     IF ( howmany /= 1 ) &
         call fftx_error__('cfft3d',' homany different from 1, not yet implemented for SX6 ', 1)

#if defined(ASL)
       ALLOCATE (cw2(ldx*ldy*ldz))
       CALL zfc3cl (f(1), nx, ny, nz, ldx, ldy, ldz, err)
#else
       ALLOCATE (cw2(6*ldx*ldy*ldz))
#endif
     !
     !   Here initialize table only if necessary
     !
     ip = -1
     DO i = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       IF ( ( nx == dims(1,i) ) .and. &
            ( ny == dims(2,i) ) .and. &
            ( nz == dims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF
     END DO

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

#if defined(ASL)
#if defined(MICRO)
       CALL hfc3fb (nx,ny,nz, f(1) , ldx, ldy, ldz, 0, &
            iw0(1,icurrent), auxp(1,icurrent), cw2(1), nbtasks, err)
#else
       CALL zfc3fb (nx,ny,nz, f(1), ldx, ldy, ldz, 0, &
             iw0(1,icurrent), auxp(1,icurrent), cw2(1), err)
#endif
#else
       ! for some reason the error variable is not set by this driver on NEC SX machines
       err = 0 
       CALL ZZFFT3D (0, nx,ny,nz, 1.0_DP, f(1), ldx, ldy, &
          &             f(1), ldx, ldy, auxp(1,icurrent), cw2(1), err)
#endif

       IF (err /= 0) CALL fftx_error__('cfft3d','FFT init returned an error ', err)

       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF

     !
     !   Now perform the 3D FFT using the machine specific driver
     !

#if defined(ASL)
#if defined(MICRO)
     CALL hfc3bf (nx,ny,nz, f(1), ldx,ldy, ldz, &
          -isign, iw0(1,ip), auxp(1,ip), cw2(1), nbtasks, err)
#else
     CALL zfc3bf (nx,ny,nz, f(1), ldx,ldy, ldz, &
          -isign, iw0(1,ip), auxp(1,ip), cw2(1), err)
#endif
     IF ( isign < 0) THEN
        tscale = 1.0_DP / DBLE( nx * ny * nz )
        call ZDSCAL( ldx * ldy * ldz, tscale, f(1), 1)
     END IF
#else
     ! for some reason the error variable is not set by this driver on NEC SX machines
     err = 0 
     tscale = 1.0_DP
     IF ( isign < 0) THEN
        tscale = tscale / DBLE( nx * ny * nz )
     END IF
     CALL ZZFFT3D (isign, nx,ny,nz, tscale, f(1), ldx,ldy, &
          f_out(1), ldx,ldy, auxp(1,ip), cw2(1), err)
!$omp parallel do private(j,i,k_off,kj_offset)
     do k=1,nz
        k_off = (k-1)*ldx*ldy
        do j=1,ny
           kj_offset = (j-1)*ldx + k_off
           do i=1,nx
              f(i+kj_offset) = f_out(i+kj_offset)
           end do
        end do
     end do
!$omp end parallel do
#endif

     IF (err /= 0) CALL fftx_error__('cfft3d','FFT returned an error ', err)
     DEALLOCATE(cw2)


     RETURN
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

SUBROUTINE cfft3ds (f, nx, ny, nz, ldx, ldy, ldz, howmany, isign, &
     do_fft_z, do_fft_y)
  !
  !     driver routine for 3d complex "reduced" fft - see cfft3d
  !     The 3D fft are computed only on lines and planes which have
  !     non zero elements. These lines and planes are defined by
  !     the two integer vectors do_fft_y(nx) and do_fft_z(ldx*ny)
  !     (1 = perform fft, 0 = do not perform fft)
  !     This routine is implemented only for fftw, essl, acml
  !     If not implemented, cfft3d is called instead
  !
  !----------------------------------------------------------------------
  !
  implicit none

  integer :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
  !
  !   logical dimensions of the fft
  !   physical dimensions of the f array
  !   sign of the transformation

  complex(DP) :: f ( ldx * ldy * ldz )
  integer :: do_fft_y(:), do_fft_z(:)
  !
  integer :: m, incx1, incx2
  INTEGER :: i, k, j, err, idir, ip,  ii, jj
  REAL(DP) :: tscale
  INTEGER, SAVE :: icurrent = 1
  INTEGER, SAVE :: dims(3,ndims) = -1

  CALL cfft3d (f, nx, ny, nz, ldx, ldy, ldz, howmany, isign)
  RETURN
END SUBROUTINE cfft3ds

!=----------------------------------------------------------------------=!
   END MODULE fft_scalar
!=----------------------------------------------------------------------=!
