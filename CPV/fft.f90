!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  ----------------------------------------------
!  This Module written by Carlo Cavazzoni 
!  Last modified April 2003
!  ----------------------------------------------

#include "f_defs.h"

!=---------------------------------------------------------------------=!
   MODULE fft
!=---------------------------------------------------------------------=!

     USE kinds
     USE cp_types
     USE fft_types, ONLY: fft_dlay_descriptor
     USE stick, ONLY: dfftp
     USE fft_scalar

     IMPLICIT NONE
     SAVE

     PRIVATE

     INTEGER, PARAMETER :: FFT_MODE_WAVE = 1
     INTEGER, PARAMETER :: FFT_MODE_POTE = 2

     REAL(dbl) :: tims(4,2) = 0.0d0

     COMPLEX(dbl), ALLOCATABLE :: zwrk(:,:)
     INTEGER, ALLOCATABLE :: ig2st(:,:)
     INTEGER :: ngw, ng

     INTERFACE pfwfft
       MODULE PROCEDURE pfwfft1, pfwfft2
     END INTERFACE
     INTERFACE pinvfft
       MODULE PROCEDURE pinvfft1, pinvfft2
     END INTERFACE

     REAL(dbl)  :: total_fft_wf_time = 0.0d0
     INTEGER    :: number_of_fft_wf  = 0
     REAL(dbl)  :: total_fft_time = 0.0_dbl
     INTEGER    :: number_of_fft  = 0
     LOGICAL    :: first          = .true.
     LOGICAL    :: tk             = .false.

     COMPLEX(dbl), PARAMETER :: czero = (0.0_dbl,0.0_dbl)
     COMPLEX(dbl), PARAMETER :: uimag = (0.0_dbl,1.0_dbl)

     REAL(dbl) :: cclock
     EXTERNAL  :: cclock

     PUBLIC :: pfwfft, pinvfft, fft_initialize, fft_closeup, fft_setup
     PUBLIC :: pw_fwfft, pw_invfft, fft_wf_initialize, fft_time_stat
     PUBLIC :: FFT_MODE_WAVE, FFT_MODE_POTE, tims

!----------------------------------------------------------------------
   CONTAINS
!----------------------------------------------------------------------

   SUBROUTINE fft_allocate_workspace( nr3x, nsp, nsw )

     IMPLICIT NONE
     INTEGER :: nr3x, nsp(:), nsw(:)
     INTEGER :: ierr, nsx

     nsx = MAX( MAXVAL( nsp ), MAXVAL( nsw ) )

     IF( ALLOCATED( zwrk ) ) THEN
       IF( SIZE( zwrk, 1 ) /= nr3x .OR. SIZE( zwrk, 2 ) /= nsx ) THEN
         DEALLOCATE(zwrk, STAT=ierr )
         IF( ierr /= 0 ) &
           CALL errore(' fft_allocate_workspace ', ' deallocation of zwrk failed ', ierr )
       END IF
     END IF

     IF( .NOT. ALLOCATED( zwrk ) ) THEN
       ALLOCATE( zwrk( nr3x, nsx ), STAT=ierr )
       IF( ierr /= 0 ) &
         CALL errore(' fft_allocate_workspace ', ' allocation of zwrk failed ', ierr )
       zwrk = 0.0d0
     END IF

     RETURN
   END SUBROUTINE

!----------------------------------------------------------------------

   SUBROUTINE fft_initialize

     COMPLEX(dbl) :: dum(1,1,1)

     CALL fft_allocate_workspace( dfftp%nr3x, dfftp%nsp, dfftp%nsw )
     CALL pc3fft_drv( dum, zwrk, 0, dfftp, FFT_MODE_POTE )

     RETURN
   END SUBROUTINE fft_initialize

!----------------------------------------------------------------------

   SUBROUTINE fft_wf_initialize

     COMPLEX(dbl) :: dum(1,1,1)

     CALL fft_allocate_workspace( dfftp%nr3x, dfftp%nsp, dfftp%nsw )
     CALL pc3fft_drv(dum, zwrk, 0, dfftp, FFT_MODE_WAVE)

     RETURN
   END SUBROUTINE fft_wf_initialize

!----------------------------------------------------------------------

   SUBROUTINE fft_closeup

     IMPLICIT NONE

     INTEGER :: ierr

     IF( ALLOCATED( zwrk ) ) THEN
       DEALLOCATE(zwrk, STAT=ierr)
       IF( ierr /= 0 ) CALL errore(' fft_closeup ', ' deallocation of zwrk failed ', ierr )
     END IF

     IF( ALLOCATED( ig2st ) ) THEN
       DEALLOCATE(ig2st, STAT=ierr)
       IF( ierr /= 0 ) CALL errore(' fft_closeup ', ' deallocation of ig2st failed ', ierr )
     END IF

     first = .TRUE.

     RETURN
   END SUBROUTINE fft_closeup

!----------------------------------------------------------------------
!
   SUBROUTINE fft_setup(gv, kp)

     USE mp_global, ONLY: nproc
     USE brillouin, ONLY: kpoints

     TYPE (recvecs), INTENT(IN) :: gv
     TYPE (kpoints), INTENT(IN) :: kp

     INTEGER :: ierr
     INTEGER :: i, ig

     !
     !  Body
     !

     tk = .NOT. ( kp%scheme == 'gamma' )

     IF( gv%ng_l < gv%ngw_l ) &
       CALL errore( ' fft_setup ', ' inconsistend number of plane waves ', 1 )

     IF( ALLOCATED( ig2st ) ) DEALLOCATE( ig2st )
     ALLOCATE( ig2st( SIZE( gv%ig2st, 1 ), gv%ng_l ) )

     DO ig = 1, gv%ng_l
       ig2st(:,ig) = gv%ig2st(:,ig)
     END DO

     ng  = gv%ng_l
     ngw = gv%ngw_l

     tims = 0.0d0

     first = .false.

     RETURN
   END SUBROUTINE fft_setup

!
!----------------------------------------------------------------------
!

!  BEGIN manual -------------------------------------------------------------

   SUBROUTINE fft_time_stat( nstep )

!  Print some statistics about time wasted by fft routines
!  --------------------------------------------------------------------------
!  END manual ---------------------------------------------------------------

     USE io_global, ONLY: stdout

     INTEGER, INTENT(IN) :: nstep
     REAL(dbl)           :: tloop, tav, ttot
     INTEGER :: i, j

     ttot   = total_fft_wf_time
     tav    = ttot / number_of_fft_wf

     WRITE( stdout,*)
     WRITE( stdout,*)   '  Wave functions FFT execution time statistics'
     WRITE( stdout,500) '  total number of ffts ..', number_of_fft_wf
     WRITE( stdout,501) '  average time per fft ..', tav
     WRITE( stdout,501) '  total fft time    .....', ttot
     WRITE( stdout,*)

     ttot   = total_fft_time
     tav    = ttot / number_of_fft

     WRITE( stdout,*)
     WRITE( stdout,*)   '  Charge density and potential FFT execution statistics'
     WRITE( stdout,500) '  total number of ffts ..', number_of_fft
     WRITE( stdout,501) '  average time per fft ..', tav
     WRITE( stdout,501) '  total fft time    .....', ttot
     WRITE( stdout,*)

     WRITE( stdout, fmt ="(/,3X,'PC3FFT TIMINGS')" )
     WRITE( stdout,910)
     WRITE( stdout,999) ( ( tims(i,j), i = 1, 4), j = 1, 2 )

910  FORMAT('    FFTXW    FFTYW    FFTZW    TRASW    FFTXP    FFTYP    FFTZP    TRASP')
999  FORMAT(8(F9.3))
500  FORMAT(1X,A,I10)
501  FORMAT(1X,A,F10.5)

     RETURN
   END SUBROUTINE fft_time_stat

!
!=---------------------------------------------------------------------==!
!
!
!     Charge Density and Potentials FFTs high level Drivers
!
!
!=---------------------------------------------------------------------==!
!
   SUBROUTINE pfwfft2( c, cpsi )
!
     !   This subroutine COMPUTE on the charge density grid :
     !      C = FWFFT( PSI )
!
      IMPLICIT NONE

      COMPLEX(dbl), INTENT(INOUT) :: cpsi(:,:,:)
      COMPLEX(dbl), INTENT(OUT) :: C(:)
      REAL(dbl) :: t1
      INTEGER :: j, ierr, k, is

      t1 = cclock()

      IF ( first ) &
        CALL errore( ' pfwfft 2 ', ' fft not initialized ', 1 )

      IF ( SIZE( cpsi, 1 ) /= dfftp%nr1x ) &
        CALL errore( ' pfwfft 2 ', ' inconsistent array dimensions ', 1 )
      IF ( SIZE( cpsi, 2 ) /= dfftp%nr2x ) &
        CALL errore( ' pfwfft 2 ', ' inconsistent array dimensions ', 2 )
      IF ( SIZE( cpsi, 3 ) /= dfftp%npl ) &
        CALL errore( ' pfwfft 2 ', ' inconsistent array dimensions ', 3 )
      IF ( SIZE( zwrk, 1 ) /= dfftp%nr3x ) &
        CALL errore( ' pfwfft 2 ', ' inconsistent array dimensions ', 4 )

      CALL pc3fft_drv(cpsi(1,1,1), zwrk(1,1), -1, dfftp, FFT_MODE_POTE)
      DO j = 1, ng
        k  = ig2st( 1, j ); is = ig2st( 2, j )
        C( j ) = zwrk( k, is )
      END DO

      total_fft_time =  total_fft_time + ( cclock() - T1 )
      number_of_fft  =  number_of_fft  + 1

     RETURN
   END SUBROUTINE pfwfft2

!----------------------------------------------------------------------

   SUBROUTINE pfwfft1( c, a )
!
     !   This subroutine COMPUTE on the charge density grid :
     !     C = FWFFT( A )
     !       C is a COMPLEX 1D array (reciprocal space)
     !       A is a REAL 3D array (real space)
!
      IMPLICIT NONE

      REAL(dbl),  INTENT(IN) :: A(:,:,:)
      COMPLEX(dbl) :: C(:)

      COMPLEX(dbl), allocatable :: psi(:,:,:)
      REAL(dbl) :: t1
      INTEGER :: is, i, j, k, nxl, nyl, nzl, ierr

      T1 = cclock()

      IF ( first ) &
        CALL errore( ' pfwfft 1 ', ' fft not initialized ', 1 )

      IF ( SIZE( A, 1 ) /= dfftp%nr1x ) &
        CALL errore( ' pfwfft 1 ', ' inconsistent array dimensions ', 1 )
      IF ( SIZE( A, 2 ) /= dfftp%nr2x ) &
        CALL errore( ' pfwfft 1 ', ' inconsistent array dimensions ', 2 )
      IF ( SIZE( A, 3 ) /= dfftp%npl ) &
        CALL errore( ' pfwfft 1 ', ' inconsistent array dimensions ', 3 )
      IF ( SIZE( zwrk, 1 ) /= dfftp%nr3x ) &
        CALL errore( ' pfwfft 1 ', ' inconsistent array dimensions ', 4 )

      nxl = SIZE(A,1)
      nyl = SIZE(A,2)
      nzl = SIZE(A,3)

      ALLOCATE( psi(nxl,nyl,nzl), STAT=ierr)
      IF( ierr /= 0 )  call errore(' PFWFFT ', ' allocation of psi failed ' ,0)

      psi = CMPLX( A )

      !CALL pc3fft_stick(psi, zwrk, -1, dfftp, FFT_MODE_POTE)
      CALL pc3fft_drv(psi(1,1,1), zwrk(1,1), -1, dfftp, FFT_MODE_POTE)

      DO j = 1, ng
        k  = ig2st( 1, j ); is = ig2st( 2, j )
        C( j ) = zwrk( k, is )
      END DO

      DEALLOCATE( psi, STAT=ierr )
      IF( ierr /= 0 )  call errore(' PFWFFT ', ' deallocation of psi failed ' ,0)

      total_fft_time =  total_fft_time + ( cclock() - T1 )
      number_of_fft  =  number_of_fft  + 1

      RETURN
   END SUBROUTINE pfwfft1

!----------------------------------------------------------------------

   SUBROUTINE pinvfft2(c, a, alpha)
!
     !   This subroutine COMPUTE on the charge density grid :
     !     C = ALPHA * C + INVFFT( A )   if alpha is present
     !     C =     INVFFT( A )           otherwise
!

      IMPLICIT NONE

      REAL(dbl), INTENT(INOUT)  :: C(:,:,:)
      REAL(dbl), INTENT(IN), OPTIONAL :: ALPHA
      COMPLEX(dbl), INTENT(IN) :: A(:)

      INTEGER :: i, j, k, nxl, nyl, nzl, is, is_i, k_i, ierr
      COMPLEX(dbl), ALLOCATABLE :: psi(:,:,:)
      REAL(dbl) t1
!
      T1 = cclock()

      IF ( first ) &
        CALL errore(' pinvfft ',' fft not initialized ', 0 )

      IF ( SIZE( c, 1 ) /= dfftp%nr1x ) &
        CALL errore( ' pinvfft 2 ', ' inconsistent array dimensions ', 1 )
      IF ( SIZE( c, 2 ) /= dfftp%nr2x ) &
        CALL errore( ' pinvfft 2 ', ' inconsistent array dimensions ', 2 )
      IF ( SIZE( c, 3 ) /= dfftp%npl ) &
        CALL errore( ' pinvfft 2 ', ' inconsistent array dimensions ', 3 )
      IF ( SIZE( zwrk, 1 ) /= dfftp%nr3x ) &
        CALL errore( ' pinvfft 2 ', ' inconsistent array dimensions ', 4 )

      nxl = SIZE( c, 1 )
      nyl = SIZE( c, 2 )
      nzl = SIZE( c, 3 )

      ALLOCATE( psi( nxl, nyl, nzl ), STAT=ierr)
      IF( ierr /= 0 )  call errore(' PFWFFT ', ' allocation of psi failed ' ,0)

      zwrk = 0.0d0

      IF(.NOT. tk) THEN
        DO j = 1, ng
          k   = ig2st(1,j); is   = ig2st(2,j)
          k_i = ig2st(3,j); is_i = ig2st(4,j)
          zwrk(k,     is) = A(j)
          zwrk(k_i, is_i) = CONJG( A(j) )
        END DO
      ELSE
        DO j = 1, ng
          k = ig2st(1,j); is = ig2st(2,j)
          zwrk(k,is) = A(j)
        END DO
      END IF

      !CALL pc3fft_stick(psi(:,:,:), zwrk(:,:), +1, dfftp, FFT_MODE_POTE)
      CALL pc3fft_drv(psi(1,1,1), zwrk(1,1), +1, dfftp, FFT_MODE_POTE)

      IF( .NOT. PRESENT( alpha ) ) THEN
        c = REAL(psi)
      ELSE
        c = c + alpha * REAL(psi)
      END IF

      DEALLOCATE(psi, STAT=ierr)
      IF( ierr /= 0 )  call errore(' PFWFFT ', ' deallocation of psi failed ' ,0)

      total_fft_time =  total_fft_time + ( cclock() - t1 )
      number_of_fft  =  number_of_fft  + 1

      RETURN
   END SUBROUTINE pinvfft2

!----------------------------------------------------------------------

   SUBROUTINE pinvfft1(gamma, c, alpha, a)

     !   This subroutine COMPUTE on the charge density grid :
     !     C = gamma * C + alpha * INVFFT( A )

      IMPLICIT NONE

      REAL(dbl), INTENT(in) :: alpha, gamma
      REAL(dbl)    :: C(:,:,:)
      COMPLEX(dbl) :: A(:)

      INTEGER :: i, j, k, k_i, is, is_i, ierr
      COMPLEX(dbl), ALLOCATABLE :: psi(:,:,:)
      REAL(dbl) :: t1
!
      T1 = cclock()

      IF ( first ) &
        CALL errore(' pinvfft 1 ',' fft not initialized ', 1 )

      IF ( SIZE( c, 1 ) /= dfftp%nr1x ) &
        CALL errore( ' pinvfft 1 ', ' inconsistent array dimensions ', 1 )
      IF ( SIZE( c, 2 ) /= dfftp%nr2x ) &
        CALL errore( ' pinvfft 1 ', ' inconsistent array dimensions ', 2 )
      IF ( SIZE( c, 3 ) /= dfftp%npl ) &
        CALL errore( ' pinvfft 1 ', ' inconsistent array dimensions ', 3 )
      IF ( SIZE( zwrk, 1 ) /= dfftp%nr3x ) &
        CALL errore( ' pinvfft 1 ', ' inconsistent array dimensions ', 4 )

      ALLOCATE( psi( SIZE(c,1), SIZE(c,2), SIZE(c,3) ), STAT=ierr)
      IF( ierr /= 0 )  call errore(' PINVFFT ', ' allocation of psi failed ', 1 )

      zwrk = 0.0d0

      IF(.NOT. tk) THEN
        DO j = 1, ng
          k   = ig2st( 1, j ); is   = ig2st( 2, j )
          k_i = ig2st( 3, j ); is_i = ig2st( 4, j )
          zwrk(k,     is) = A(j)
          zwrk(k_i, is_i) = CONJG( A(j) )
        END DO
      ELSE
        DO j = 1, ng
          k   = ig2st(1,j);   is   = ig2st(2,j)
          zwrk(k,     is) = A(j)
        END DO
      END IF

      !CALL pc3fft_stick(psi, zwrk, +1, dfftp, FFT_MODE_POTE)
      CALL pc3fft_drv(psi(1,1,1), zwrk(1,1), +1, dfftp, FFT_MODE_POTE)

      c = gamma * c + alpha * REAL( psi )

      DEALLOCATE(psi, STAT=ierr)
      IF( ierr /= 0 )  call errore(' PINVFFT ', ' deallocation of psi failed ' ,0)

      total_fft_time =  total_fft_time + ( cclock() - t1 )
      number_of_fft  =  number_of_fft  + 1

      RETURN
   END SUBROUTINE pinvfft1

!=---------------------------------------------------------------------==!
!
!
!     Wave Function FFTs high level Drivers
!
!
!=---------------------------------------------------------------------==!

   SUBROUTINE pw_fwfft( psi, c, ca )

     !   This subroutine COMPUTE on the wave function sub-grid :

     USE stick, ONLY: dfftp

     COMPLEX(dbl) :: C(:)
     COMPLEX(dbl), OPTIONAL :: CA(:)
     COMPLEX(dbl) :: psi(:,:,:)
     REAL(dbl)  :: T1
     INTEGER :: is, is_i, k, k_i, j

     t1 = cclock()

     IF ( first ) &
       CALL errore(' pw_fwfft 1 ',' fft not initialized ', 1 )

     IF ( SIZE( psi, 1 ) /= dfftp%nr1x ) &
       CALL errore( ' pw_fwfft 1 ', ' inconsistent array dimensions ', 1 )
     IF ( SIZE( psi, 2 ) /= dfftp%nr2x ) &
       CALL errore( ' pw_fwfft 1 ', ' inconsistent array dimensions ', 2 )
     IF ( SIZE( psi, 3 ) /= dfftp%npl ) &
       CALL errore( ' pw_fwfft 1 ', ' inconsistent array dimensions ', 3 )
     IF ( SIZE( zwrk, 1 ) /= dfftp%nr3x ) &
       CALL errore( ' pw_fwfft 1 ', ' inconsistent array dimensions ', 4 )

     !CALL pc3fft_stick(psi, zwrk, -1, dfftp, FFT_MODE_WAVE)
     CALL pc3fft_drv(psi(1,1,1), zwrk(1,1), -1, dfftp, FFT_MODE_WAVE)

     IF( PRESENT( ca ) ) THEN
       DO j = 1, ngw
         k   = ig2st(1,j); is   = ig2st(2,j)
         k_i = ig2st(3,j); is_i = ig2st(4,j)
         ca(j) = zwrk(k_i, is_i)
         c (j) = zwrk(k  , is  )
       END DO
     ELSE
       DO j = 1, ngw
         k = ig2st(1,j); is = ig2st(2,j)
         c(j) = zwrk(k,is)
       END DO
     END IF

     total_fft_wf_time =  total_fft_wf_time + ( cclock() - t1 )
     number_of_fft_wf  =  number_of_fft_wf  + 1

     RETURN
   END SUBROUTINE pw_fwfft

!
!==---------------------------------------------------------------------==!
!

   SUBROUTINE pw_invfft( psi, c, ca )

     !   This subroutine COMPUTE on the wave function sub-grid :

     USE stick, ONLY: dfftp

     COMPLEX(dbl), INTENT(IN) :: C(:)
     COMPLEX(dbl), INTENT(IN), OPTIONAL :: CA(:)
     COMPLEX(dbl) :: psi(:,:,:)

     COMPLEX(dbl) :: cca1, cca2
     REAL(dbl) :: T1, T2
     INTEGER :: is, is_i, k, k_i, j

     IF ( first ) &
       CALL errore(' pw_invfft 1 ',' fft not initialized ', 1 )

     IF ( SIZE( psi, 1 ) /= dfftp%nr1x ) &
       CALL errore( ' pw_invfft 1 ', ' inconsistent array dimensions ', 1 )
     IF ( SIZE( psi, 2 ) /= dfftp%nr2x ) &
       CALL errore( ' pw_invfft 1 ', ' inconsistent array dimensions ', 2 )
     IF ( SIZE( psi, 3 ) /= dfftp%npl ) &
       CALL errore( ' pw_invfft 1 ', ' inconsistent array dimensions ', 3 )
     IF ( SIZE( zwrk, 1 ) /= dfftp%nr3x ) &
       CALL errore( ' pw_invfft 1 ', ' inconsistent array dimensions ', 4 )

     T1 = cclock()

     zwrk = czero

     ! WRITE(6,*) 'DEBUG pw_invfft ', ngw
     IF( PRESENT( ca ) ) THEN
       DO j = 1, ngw
         cca1 =       C(j)  + uimag *       CA(j)
         cca2 = CONJG(C(j)) + uimag * CONJG(CA(j))
         k   = ig2st(1,j); is   = ig2st(2,j)
         k_i = ig2st(3,j); is_i = ig2st(4,j)
         zwrk(k  , is  ) = cca1
         zwrk(k_i, is_i) = cca2
         ! WRITE( *, * ) k, is, k_i, is_i
         ! WRITE( *, fmt='(2D24.16)' ) zwrk(k  , is  )
         ! WRITE( *, fmt='(2D24.16)' ) zwrk(k_i, is_i)
       END DO
     ELSE
       DO j = 1, ngw
         k = ig2st(1,j); is = ig2st(2,j)
         zwrk(k,is) = c(j)
       END DO
     END IF

     CALL pc3fft_drv(psi(1,1,1), zwrk(1,1), +1, dfftp, FFT_MODE_WAVE)

     total_fft_wf_time =  total_fft_wf_time + ( cclock() - T1 )
     number_of_fft_wf  =  number_of_fft_wf  + 1

     RETURN
   END SUBROUTINE pw_invfft



!==---------------------------------------------------------------------==!
   END MODULE fft
!==---------------------------------------------------------------------==!
