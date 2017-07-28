!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------=!
   MODULE fft_smallbox
!=----------------------------------------------------------------------=!

!! iso_c_binding provides C_PTR, C_NULL_PTR, C_ASSOCIATED
       USE iso_c_binding
       IMPLICIT NONE
       SAVE

        PRIVATE
        PUBLIC :: cft_b, cft_b_omp_init, cft_b_omp

! ...   Local Parameter

        INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

        !   ndims   Number of different FFT tables that the module
        !           could keep into memory without reinitialization

        INTEGER, PARAMETER :: ndims = 3

        !   Workspace that is statically allocated is defined here
        !   in order to avoid multiple copies of the same workspace
        !   lwork:   Dimension of the work space array (if any)

        INTEGER   :: cft_b_dims( 4 )
!$omp threadprivate (cft_b_dims)
        TYPE(C_PTR) :: cft_b_bw_planz = C_NULL_PTR
!$omp threadprivate (cft_b_bw_planz)
        TYPE(C_PTR) :: cft_b_bw_planx = C_NULL_PTR
!$omp threadprivate (cft_b_bw_planx)
        TYPE(C_PTR) :: cft_b_bw_plany = C_NULL_PTR
!$omp threadprivate (cft_b_bw_plany)

!=----------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------=!

!
!=----------------------------------------------------------------------=!
!
!
!
!         3D parallel FFT on sub-grids
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cft_b ( f, nx, ny, nz, ldx, ldy, ldz, imin3, imax3, sgn )

!     driver routine for 3d complex fft's on box grid, parallel case
!     fft along xy is done only on planes that correspond to dense grid
!     planes on the current processor, i.e. planes with imin3 <= nz <= imax3
!     implemented for FFTW, only for sgn=1 (f(R) => f(G))
!     (beware: here the "essl" convention for the sign of the fft is used!)
!
      USE fftw_interfaces
      implicit none
      integer nx,ny,nz,ldx,ldy,ldz,imin3,imax3,sgn
      complex(dp) :: f(:)

      integer isign, naux, ibid, nplanes, nstart, k
      real(DP) :: tscale

      integer :: ip, i
      integer, save :: icurrent = 1
      integer, save :: dims( 4, ndims ) = -1

      TYPE(C_PTR), save :: bw_planz(  ndims ) = C_NULL_PTR
      TYPE(C_PTR), save :: bw_planx(  ndims ) = C_NULL_PTR
      TYPE(C_PTR), save :: bw_plany(  ndims ) = C_NULL_PTR
      TYPE(C_PTR), save :: bw_planxy( ndims ) = C_NULL_PTR

      isign = -sgn
      tscale = 1.0_DP

      if ( isign > 0 ) then
         call fftx_error__('cft_b','not implemented',isign)
      end if
!
! 2d fft on xy planes - only needed planes are transformed
! note that all others are left in an unusable state
!
      nplanes = imax3 - imin3 + 1
      nstart  = ( imin3 - 1 ) * ldx * ldy + 1

      !
      !   Here initialize table only if necessary
      !

      ip = -1
      DO i = 1, ndims

        !   first check if there is already a table initialized
        !   for this combination of parameters

        IF ( ( nx == dims(1,i) ) .and. ( ny == dims(2,i) ) .and. &
             ( nz == dims(3,i) ) .and. ( nplanes == dims(4,i) ) ) THEN
           ip = i
           EXIT
        END IF

      END DO

      IF( ip == -1 ) THEN

        !   no table exist for these parameters
        !   initialize a new one

        if ( C_ASSOCIATED(bw_planz(icurrent)) ) &
             call DESTROY_PLAN_1D( bw_planz(icurrent) )
        call CREATE_PLAN_1D( bw_planz(icurrent), nz, 1 )

        if ( C_ASSOCIATED(bw_planx(icurrent)) ) &
             call DESTROY_PLAN_1D( bw_planx(icurrent) )
        call CREATE_PLAN_1D( bw_planx(icurrent), nx, 1 )

        if ( C_ASSOCIATED(bw_plany(icurrent)) ) &
             call DESTROY_PLAN_1D( bw_plany(icurrent) )
        call CREATE_PLAN_1D( bw_plany(icurrent), ny, 1 )

        if ( C_ASSOCIATED(bw_planxy(icurrent)) ) &
             call DESTROY_PLAN_2D( bw_planxy(icurrent) )
        call CREATE_PLAN_2D( bw_planxy(icurrent), nx, ny, 1 )
!
        dims(1,icurrent) = nx; dims(2,icurrent) = ny
        dims(3,icurrent) = nz; dims(4,icurrent) = nplanes
        ip = icurrent
        icurrent = MOD( icurrent, ndims ) + 1

      END IF

      !
      !  fft along Z
      !
      call FFTW_INPLACE_DRV_1D( bw_planz(ip), ldx*ldy, f(1), ldx*ldy, 1 )
      !
      !  fft along Y
      !  fft along X
      !
      do k = imin3, imax3
        call FFTW_INPLACE_DRV_1D( bw_plany(ip), nx, f((k-1)*ldx*ldy + 1), ldx, 1 )
        call FFTW_INPLACE_DRV_1D( bw_planx(ip), ny, f((k-1)*ldx*ldy + 1), 1, ldx )
      end do   

      RETURN
   END SUBROUTINE cft_b

!
!=----------------------------------------------------------------------=!
!
!
!
!   3D parallel FFT on sub-grids, to be called inside OpenMP region
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cft_b_omp_init ( nx, ny, nz )

!     driver routine for 3d complex fft's on box grid, init subroutine
!
      USE fftw_interfaces
      implicit none
      integer, INTENT(IN) :: nx,ny,nz
      !
      !   Here initialize table 
      !
!$omp parallel

      IF( .NOT. C_ASSOCIATED(cft_b_bw_planz) ) THEN
         CALL CREATE_PLAN_1D( cft_b_bw_planz, nz, 1 )
         cft_b_dims(3) = nz
      END IF
      IF( .NOT. C_ASSOCIATED(cft_b_bw_planx) ) THEN
         CALL CREATE_PLAN_1D( cft_b_bw_planx, nx, 1 )
         cft_b_dims(1) = nx
      END IF
      IF( .NOT. C_ASSOCIATED(cft_b_bw_plany) ) THEN
         CALL CREATE_PLAN_1D( cft_b_bw_plany, ny, 1 )
         cft_b_dims(2) = ny
      END IF

!$omp end parallel

     RETURN
   END SUBROUTINE cft_b_omp_init


   SUBROUTINE cft_b_omp ( f, nx, ny, nz, ldx, ldy, ldz, imin3, imax3, sgn )

!     driver routine for 3d complex fft's on box grid, parallel (MPI+OpenMP) case
!     fft along xy is done only on planes that correspond to dense grid
!     planes on the current processor, i.e. planes with imin3 <= nz <= imax3
!     implemented ONLY for internal fftw, and only for sgn=1 (f(R) => f(G))
!     (beware: here the "essl" convention for the sign of the fft is used!)
!
!     This driver is meant for calls inside parallel OpenMP sections
!
      USE fftw_interfaces
      implicit none
      integer, INTENT(IN) :: nx,ny,nz,ldx,ldy,ldz,imin3,imax3,sgn
      complex(dp) :: f(:)

      INTEGER, SAVE :: k
!$omp threadprivate (k)

      if ( -sgn > 0 ) then
         CALL fftx_error__('cft_b_omp','forward transform not implemented',1)
      end if

      IF ( .NOT. C_ASSOCIATED(cft_b_bw_planz) .or. &
           .NOT. C_ASSOCIATED(cft_b_bw_planx) .or. &
           .NOT. C_ASSOCIATED(cft_b_bw_plany) ) THEN
         CALL fftx_error__('cft_b_omp','plan not initialized',1)
      END IF

      !  consistency check

      IF ( ( nx /= cft_b_dims(1) ) .or. ( ny /= cft_b_dims(2) ) .or. ( nz /= cft_b_dims(3) ) ) THEN
         CALL fftx_error__('cft_b_omp', 'dimensions are inconsistent with the existing plan',1) 
      END IF

      !  fft along Z
      !
      call FFTW_INPLACE_DRV_1D( cft_b_bw_planz, ldx*ldy, f(1), ldx*ldy, 1 )
      !
      !  fft along Y
      !  fft along X
      !
      do k = imin3, imax3
        call FFTW_INPLACE_DRV_1D( cft_b_bw_plany, nx, f((k-1)*ldx*ldy + 1), ldx, 1 )
        call FFTW_INPLACE_DRV_1D( cft_b_bw_planx, ny, f((k-1)*ldx*ldy + 1), 1, ldx )
      end do   

     RETURN
   END SUBROUTINE cft_b_omp


!=----------------------------------------------------------------------=!
   END MODULE fft_smallbox
!=----------------------------------------------------------------------=!

