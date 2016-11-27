!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

#include "fft_defs.h"

!=----------------------------------------------------------------------=!
   MODULE fft_scalar_FFTW
!=----------------------------------------------------------------------=!

       USE fft_param
       
       IMPLICIT NONE
       SAVE

       PRIVATE
       PUBLIC :: cft_1z, cft_2xy, cfft3d, cfft3ds

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

     !   Pointers to the "C" structures containing FFT factors ( PLAN )
     !   C_POINTER is defined in include/fft_defs.h
     !   for 32bit executables, C_POINTER is integer(4)
     !   for 64bit executables, C_POINTER is integer(8)

     C_POINTER, SAVE :: fw_planz( ndims ) = 0
     C_POINTER, SAVE :: bw_planz( ndims ) = 0

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

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'cft_1z' )
#endif


#if defined(__OPENMP)

     ldz_t = ldz
     !
     IF (isign < 0) THEN
!$omp parallel default(none) private(tid,offset,i,tscale) shared(c,isign,nsl,fw_planz,ip,nz,cout,ldz) &
!$omp &        firstprivate(ldz_t)
!$omp do
       DO i=1, nsl
          offset = 1 + ((i-1)*ldz_t)
          CALL FFT_Z_STICK_SINGLE(fw_planz( ip), c(offset), ldz_t)
       END DO
!$omp end do
!$omp end parallel
       tscale = 1.0_DP / nz
       cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl ) * tscale
     ELSE IF (isign > 0) THEN
!$omp parallel default(none) private(tid,offset,i) shared(c,isign,nsl,bw_planz,ip,cout,ldz) &
!$omp &        firstprivate(ldz_t)
!$omp do
       DO i=1, nsl
          offset =  1 + ((i-1)* ldz_t)
          CALL FFT_Z_STICK_SINGLE(bw_planz( ip), c(offset), ldz_t)
       END DO
!$omp end do
!$omp workshare
       cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl )
!$omp end workshare
!$omp end parallel
     END IF


#else

     IF (isign < 0) THEN
        CALL FFT_Z_STICK(fw_planz( ip), c(1), ldz, nsl)
        tscale = 1.0_DP / nz
        cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl ) * tscale
     ELSE IF (isign > 0) THEN
        CALL FFT_Z_STICK(bw_planz( ip), c(1), ldz, nsl)
        cout( 1 : ldz * nsl ) = c( 1 : ldz * nsl )
     END IF

#endif

#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'cft_1z' )
#endif

     RETURN

   CONTAINS

     SUBROUTINE lookup()
     DO ip = 1, ndims
        !   first check if there is already a table initialized
        !   for this combination of parameters
        found = ( nz == zdims(1,ip) )
        IF (found) EXIT
     END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()
       IF( fw_planz( icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_planz( icurrent) )
       IF( bw_planz( icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_planz( icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_planz( icurrent), nz, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_planz( icurrent), nz, idir)
       zdims(1,icurrent) = nz; zdims(2,icurrent) = nsl; zdims(3,icurrent) = ldz;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

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


#if defined(__FFTW_ALL_XY_PLANES)
     C_POINTER, SAVE :: fw_plan_2d( ndims ) = 0
     C_POINTER, SAVE :: bw_plan_2d( ndims ) = 0
#else
     C_POINTER, SAVE :: fw_plan( 2, ndims ) = 0
     C_POINTER, SAVE :: bw_plan( 2, ndims ) = 0
#endif


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

       CALL init_plan()

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'cft_2xy' )
#endif


#if defined(__FFTW_ALL_XY_PLANES)

     IF( isign < 0 ) THEN
        !
        tscale = 1.0_DP / ( nx * ny )
        !
        CALL fftw_inplace_drv_2d( fw_plan_2d(ip), nzl, r(1), 1, ldx*ldy )
        CALL ZDSCAL( ldx * ldy * nzl, tscale, r(1), 1)
        !
     ELSE IF( isign > 0 ) THEN
        !
        CALL fftw_inplace_drv_2d( bw_plan_2d(ip), nzl, r(1), 1, ldx*ldy )
        !
     END IF

#elif defined(__OPENMP)

     nx_t  = nx
     ny_t  = ny
     nzl_t = nzl
     ldx_t = ldx 
     ldy_t = ldy
     !
     IF( isign < 0 ) THEN
        !
        tscale = 1.0_DP / ( nx * ny )
        !
!$omp parallel default(none) private(offset,itid,mytid,ntids,k,j,i)  shared(r,dofft,ip,fw_plan,nzl,nx,ny,ldx,ldy,tscale)  &
!$omp & firstprivate(nx_t, ny_t, nzl_t, ldx_t, ldy_t)

!$omp do
        DO i=1,nzl
           offset = 1+ ((i-1)*(ldx_t*ldy_t))
           CALL FFT_X_STICK_SINGLE( fw_plan(1,ip), r(offset), nx_t, ny_t, nzl_t, ldx_t, ldy_t )
        END DO
!$omp end do

        mytid = omp_get_thread_num()  ! take the thread ID
        ntids = omp_get_num_threads() ! take the number of threads
        itid  = 0

        do i = 1, nx
          do k = 1, nzl
            IF( dofft( i ) ) THEN
              IF( itid == mytid ) THEN
                j = i + ldx_t*ldy_t * ( k - 1 )
                call FFT_Y_STICK(fw_plan(2,ip), r(j), ny_t, ldx_t)
              END IF
              itid = MOD( itid + 1, ntids )
            END IF
          end do
        end do

!$omp barrier
 
!$omp workshare
        r = r * tscale
!$omp end workshare

!$omp end parallel
        !
     ELSE IF( isign > 0 ) THEN
        !
!$omp parallel default(none) private(offset,itid,mytid,ntids,k,j,i) shared(r,nx,nzl,dofft,ip,bw_plan) &
!$omp & firstprivate(nx_t, ny_t, nzl_t, ldx_t, ldy_t)

        mytid = omp_get_thread_num()  ! take the thread ID
        ntids = omp_get_num_threads() ! take the number of threads
        itid  = 0

        do i = 1, nx
          do k = 1, nzl
            IF( dofft( i ) ) THEN
              IF( itid == mytid ) THEN
                j = i + ldx_t*ldy_t * ( k - 1 )
                call FFT_Y_STICK( bw_plan(2,ip), r(j), ny_t, ldx_t)
              END IF
              itid = MOD( itid + 1, ntids )
            END IF
          end do
        end do

!$omp barrier

!$omp do
        DO i=1,nzl
           offset = 1+ ((i-1)*(ldx_t*ldy_t))
           CALL FFT_X_STICK_SINGLE( bw_plan(1,ip), r(offset), nx_t, ny_t, nzl_t, ldx_t, ldy_t )
        END DO
!$omp end do
!$omp end parallel
        !
     END IF

#else

     IF( isign < 0 ) THEN

       CALL FFT_X_STICK( fw_plan(1,ip), r(1), nx, ny, nzl, ldx, ldy )

       do i = 1, nx
         do k = 1, nzl
           IF( dofft( i ) ) THEN
             j = i + ldx*ldy * ( k - 1 )
             call FFT_Y_STICK(fw_plan(2,ip), r(j), ny, ldx)
           END IF
         end do
       end do

       tscale = 1.0_DP / ( nx * ny )
       CALL ZDSCAL( ldx * ldy * nzl, tscale, r(1), 1)

     ELSE IF( isign > 0 ) THEN

       do i = 1, nx
         do k = 1, nzl
           IF( dofft( i ) ) THEN
             j = i + ldx*ldy * ( k - 1 )
             call FFT_Y_STICK( bw_plan(2,ip), r(j), ny, ldx)
           END IF
         end do
       end do

       CALL FFT_X_STICK( bw_plan(1,ip), r(1), nx, ny, nzl, ldx, ldy )

    END IF

#endif


#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'cft_2xy' )
#endif

     RETURN

   CONTAINS

     SUBROUTINE lookup()
     DO ip = 1, ndims
       !   first check if there is already a table initialized
       !   for this combination of parameters
       found = ( ny == dims(1,ip) ) .AND. ( nx == dims(3,ip) )
       found = found .AND. ( ldx == dims(2,ip) ) .AND.  ( nzl == dims(4,ip) )
       IF (found) EXIT
     END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()
#if defined __FFTW_ALL_XY_PLANES
       IF( fw_plan_2d( icurrent) /= 0 )  CALL DESTROY_PLAN_2D(fw_plan_2d(icurrent) )
       IF( bw_plan_2d( icurrent) /= 0 )  CALL DESTROY_PLAN_2D(bw_plan_2d(icurrent) )
       idir = -1; CALL CREATE_PLAN_2D( fw_plan_2d(icurrent), nx, ny, idir)
       idir =  1; CALL CREATE_PLAN_2D( bw_plan_2d(icurrent), nx, ny, idir)
#else
       IF( fw_plan( 2,icurrent) /= 0 )   CALL DESTROY_PLAN_1D( fw_plan( 2,icurrent) )
       IF( bw_plan( 2,icurrent) /= 0 )   CALL DESTROY_PLAN_1D( bw_plan( 2,icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 2,icurrent), ny, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 2,icurrent), ny, idir)

       IF( fw_plan( 1,icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_plan( 1,icurrent) )
       IF( bw_plan( 1,icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_plan( 1,icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 1,icurrent), nx, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 1,icurrent), nx, idir)
#endif
       dims(1,icurrent) = ny; dims(2,icurrent) = ldx;
       dims(3,icurrent) = nx; dims(4,icurrent) = nzl;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

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

     C_POINTER, save :: fw_plan(ndims) = 0
     C_POINTER, save :: bw_plan(ndims) = 0

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

     IF( isign < 0 ) THEN

       call FFTW_INPLACE_DRV_3D( fw_plan(ip), 1, f(1), 1, 1 )

       tscale = 1.0_DP / DBLE( nx * ny * nz )
       call ZDSCAL( nx * ny * nz, tscale, f(1), 1)

     ELSE IF( isign > 0 ) THEN

       call FFTW_INPLACE_DRV_3D( bw_plan(ip), 1, f(1), 1, 1 )

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
            ( nz == dims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF
     END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()
       IF ( nx /= ldx .or. ny /= ldy .or. nz /= ldz ) &
         call fftx_error__('cfft3','not implemented',1)
       IF( fw_plan(icurrent) /= 0 ) CALL DESTROY_PLAN_3D( fw_plan(icurrent) )
       IF( bw_plan(icurrent) /= 0 ) CALL DESTROY_PLAN_3D( bw_plan(icurrent) )
       idir = -1; CALL CREATE_PLAN_3D( fw_plan(icurrent), nx, ny, nz, idir)
       idir =  1; CALL CREATE_PLAN_3D( bw_plan(icurrent), nx, ny, nz, idir)
       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

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
  !     the two integer vectors do_fft_y(nx) and do_fft_z(ldx*ldy)
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

  complex(DP) :: f ( ldx * ldy * ldz * howmany )
  integer :: do_fft_z(:), do_fft_y(:)
  !
  integer :: m, incx1, incx2
  INTEGER :: i, k, j, err, idir, ip,  ii, jj, h, ldh
  REAL(DP) :: tscale
  INTEGER, SAVE :: icurrent = 1
  INTEGER, SAVE :: dims(3,ndims) = -1

  C_POINTER, SAVE :: fw_plan ( 3, ndims ) = 0
  C_POINTER, SAVE :: bw_plan ( 3, ndims ) = 0

  tscale = 1.0_DP
  ldh = ldx * ldy * ldz

  IF( ny /= ldy ) &
    CALL fftx_error__(' cfft3ds ', ' wrong dimensions: ny /= ldy ', 1 )
  IF( howmany < 1 ) &
    CALL fftx_error__(' cfft3ds ', ' howmany less than one ', 1 )

     CALL lookup()

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       CALL init_plan()

     END IF


     IF ( isign > 0 ) THEN

        DO h = 0, howmany - 1
           !
           !  k-direction ...
           !

           incx1 = ldx * ldy;  incx2 = 1;  m = 1

           do i =1, nx
              do j = 1, ny
                 ii = i + ldx * (j-1)
                 if ( do_fft_z(ii) > 0 ) then
                    call FFTW_INPLACE_DRV_1D( bw_plan( 3, ip), m, f( ii + h*ldh ), incx1, incx2 )
                 end if
              end do
           end do

           !
           !  ... j-direction ...
           !

           incx1 = ldx;  incx2 = ldx*ldy;  m = nz
   
           do i = 1, nx
              if ( do_fft_y( i ) == 1 ) then
                call FFTW_INPLACE_DRV_1D( bw_plan( 2, ip), m, f( i + h*ldh ), incx1, incx2 )
              endif
           enddo

           !
           !  ... i - direction
           !

           incx1 = 1;  incx2 = ldx;  m = ldy*nz

           call FFTW_INPLACE_DRV_1D( bw_plan( 1, ip), m, f( 1 + h*ldh ), incx1, incx2 )

        END DO

     ELSE

        DO h = 0, howmany - 1
           !
           !  i - direction ...
           !

           incx1 = 1;  incx2 = ldx;  m = ldy*nz

           call FFTW_INPLACE_DRV_1D( fw_plan( 1, ip), m, f( 1 + h*howmany ), incx1, incx2 )

           !
           !  ... j-direction ...
           !

           incx1 = ldx;  incx2 = ldx*ldy;  m = nz

           do i = 1, nx
              if ( do_fft_y ( i ) == 1 ) then
                call FFTW_INPLACE_DRV_1D( fw_plan( 2, ip), m, f( i + h*howmany ), incx1, incx2 )
              endif
           enddo

           !
           !  ... k-direction
           !

           incx1 = ldx * ny;  incx2 = 1;  m = 1
 
           do i = 1, nx
              do j = 1, ny
                 ii = i + ldx * (j -1)
                 if ( do_fft_z ( ii ) > 0 ) then
                    call FFTW_INPLACE_DRV_1D( fw_plan( 3, ip), m, f( ii + h*howmany ), incx1, incx2 )
                 end if
              end do
           end do

           call DSCAL (2 * ldx * ldy * nz, 1.0_DP/(nx * ny * nz), f(1+ h*howmany ), 1)
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
       IF( fw_plan( 1, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_plan( 1, icurrent) )
       IF( bw_plan( 1, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_plan( 1, icurrent) )
       IF( fw_plan( 2, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_plan( 2, icurrent) )
       IF( bw_plan( 2, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_plan( 2, icurrent) )
       IF( fw_plan( 3, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_plan( 3, icurrent) )
       IF( bw_plan( 3, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_plan( 3, icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 1, icurrent), nx, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 1, icurrent), nx, idir)
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 2, icurrent), ny, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 2, icurrent), ny, idir)
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 3, icurrent), nz, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 3, icurrent), nz, idir)

       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

   END SUBROUTINE cfft3ds

!=----------------------------------------------------------------------=!
 END MODULE fft_scalar_FFTW
!=----------------------------------------------------------------------=!

