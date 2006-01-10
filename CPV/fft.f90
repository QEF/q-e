!
! Copyright (C) 2002-2005 FPMD-CPV groups
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

     IMPLICIT NONE
     SAVE

     PRIVATE

     INTEGER, PARAMETER :: FFT_MODE_WAVE = 2
     INTEGER, PARAMETER :: FFT_MODE_POTE = 1

     INTEGER :: ngw, ng, ngs

     INTERFACE pfwfft
       MODULE PROCEDURE pfwfft1, pfwfft2
     END INTERFACE
     INTERFACE pinvfft
       MODULE PROCEDURE pinvfft2
     END INTERFACE

     REAL(DP)  :: total_fft_wf_time = 0.0d0
     INTEGER    :: number_of_fft_wf  = 0
     REAL(DP)  :: total_fft_time = 0.0d0
     INTEGER    :: number_of_fft  = 0
     LOGICAL    :: first          = .true.
     LOGICAL    :: tk             = .false.

     REAL(DP), EXTERNAL :: cclock

     PUBLIC :: pfwfft, pinvfft, fft_closeup, fft_setup
     PUBLIC :: pw_fwfft, pw_invfft, fft_time_stat

!----------------------------------------------------------------------
   CONTAINS
!----------------------------------------------------------------------

   SUBROUTINE fft_closeup

     IMPLICIT NONE

     ng  = 0
     ngw = 0
     ngs = 0 
     first = .TRUE.

     RETURN
   END SUBROUTINE fft_closeup

!----------------------------------------------------------------------
!
   SUBROUTINE fft_setup( lgamma_ , ng_ , ngs_ , ngw_ )

     INTEGER, INTENT(IN) :: ng_ , ngs_ , ngw_
     LOGICAL, INTENT(IN) :: lgamma_

     tk = .NOT. lgamma_

     IF( ng_ < ngw_ .OR. ngs_ < ngw_ ) &
       CALL errore( ' fft_setup ', ' inconsistend number of plane waves ', 1 )

     ng  = ng_
     ngw = ngw_
     ngs = ngs_

     first = .false.

     RETURN
   END SUBROUTINE fft_setup

!
!----------------------------------------------------------------------
!

   SUBROUTINE fft_time_stat( nstep )

     !  Print some statistics about time wasted by fft routines

     USE io_global, ONLY: stdout
     USE fft_base, ONLY: fft_timing

     INTEGER, INTENT(IN) :: nstep
     REAL(DP)           :: tloop, tav, ttot
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
     WRITE( stdout,999) ( ( fft_timing(i,j), i = 1, 4), j = 1, 2 )

910  FORMAT('    FFTXP    FFTYP    FFTZP    TRASP    FFTXW    FFTYW    FFTZW    TRASW')
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

     USE fft_types, ONLY: fft_dlay_descriptor
     USE fft_base, ONLY: dfftp
     USE mp_global, ONLY: mpime
     USE io_global, ONLY: stdout

      IMPLICIT NONE

      COMPLEX(DP), INTENT(INOUT) :: cpsi(:)
      COMPLEX(DP), INTENT(OUT) :: C(:)
      COMPLEX(DP), ALLOCATABLE :: psi(:)
      COMPLEX(DP), ALLOCATABLE :: zwrk(:) 
      REAL(DP) :: t1
      INTEGER :: ierr

      t1 = cclock()

      IF ( first ) &
        CALL errore( ' pfwfft 2 ', ' fft not initialized ', 1 )

      IF ( SIZE( cpsi ) /= dfftp%nnr ) THEN
        WRITE( stdout, * ) 'Values = ', SIZE( cpsi ), dfftp%nnr
        CALL errore( ' pfwfft 2 ', ' inconsistent array dimensions ', 1 )
      END IF

#if defined __PARA

      ALLOCATE( zwrk( dfftp%nsp( mpime + 1 ) * dfftp%nr3x ) )

      CALL pc3fft_drv(cpsi(1), zwrk, -1, dfftp, FFT_MODE_POTE)

      CALL psi2c( zwrk, SIZE( zwrk ), c(1), c(1), ng, 10 )

      DEALLOCATE( zwrk )

#else

      ALLOCATE( psi( SIZE( cpsi ) ) )

      psi = cpsi

      CALL fwfft( psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )

      CALL psi2c( psi, dfftp%nnr, c(1), c(1), ng, 10 )

      DEALLOCATE( psi )

#endif

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

     USE fft_types, ONLY: fft_dlay_descriptor
     USE fft_base, ONLY: dfftp
     USE mp_global, ONLY: mpime
     use recvecs_indexes, only: nm, np

      IMPLICIT NONE

      REAL(DP),  INTENT(IN) :: A(:)
      COMPLEX(DP) :: C(:)

      COMPLEX(DP), allocatable :: psi(:)
      COMPLEX(DP), ALLOCATABLE :: zwrk(:) 
      REAL(DP) :: t1
      INTEGER :: ierr, ig, k, is

      T1 = cclock()

      IF ( first ) &
        CALL errore( ' pfwfft 1 ', ' fft not initialized ', 1 )

      IF ( SIZE( A ) /= dfftp%nnr ) &
        CALL errore( ' pfwfft 1 ', ' inconsistent array dimensions ', 1 )

      ALLOCATE( psi( SIZE(A) ), STAT=ierr)
      IF( ierr /= 0 )  call errore(' PFWFFT ', ' allocation of psi failed ' ,0)

      psi = CMPLX( A, 0.d0 )

#if defined __PARA

      ALLOCATE( zwrk( dfftp%nsp( mpime + 1 ) * dfftp%nr3x ) )

      CALL pc3fft_drv(psi(1), zwrk, -1, dfftp, FFT_MODE_POTE)

      CALL psi2c( zwrk(1), SIZE( zwrk ), c(1), c(1), ng, 10 )

      DEALLOCATE( zwrk )


#else

      CALL fwfft( psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )
      CALL psi2c( psi, dfftp%nnr, c(1), c(1), ng, 10 )

#endif

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
     USE fft_types, ONLY: fft_dlay_descriptor
     USE fft_base, ONLY: dfftp
     USE mp_global, ONLY: mpime

      IMPLICIT NONE

      REAL(DP), INTENT(INOUT)  :: C(:)
      REAL(DP), INTENT(IN), OPTIONAL :: ALPHA
      COMPLEX(DP), INTENT(IN) :: A(:)

      INTEGER :: ierr
      COMPLEX(DP), ALLOCATABLE :: psi(:)
      COMPLEX(DP), ALLOCATABLE :: zwrk(:) 
      REAL(DP) t1
!
      T1 = cclock()

      IF ( first ) &
        CALL errore(' pinvfft ',' fft not initialized ', 0 )

      IF ( SIZE( c ) /= dfftp%nnr ) &
        CALL errore( ' pinvfft 2 ', ' inconsistent array dimensions ', 1 )

      ALLOCATE( psi( SIZE( c ) ), STAT=ierr)
      IF( ierr /= 0 )  call errore(' PFWFFT ', ' allocation of psi failed ' ,0)

#if defined __PARA

      ALLOCATE( zwrk( dfftp%nsp( mpime + 1 ) * dfftp%nr3x ) )

      IF( tk ) THEN
        CALL c2psi( zwrk, SIZE( zwrk ), a(1), a(1), ng, 10 )
      ELSE
        CALL c2psi( zwrk, SIZE( zwrk ), a(1), a(1), ng, 11 )
      END IF

      CALL pc3fft_drv(psi(1), zwrk, +1, dfftp, FFT_MODE_POTE)

      DEALLOCATE( zwrk )

#else

      IF( tk ) THEN
        CALL c2psi( psi, dfftp%nnr, a(1), a(1), ng, 10 )
      ELSE
        CALL c2psi( psi, dfftp%nnr, a(1), a(1), ng, 11 )
      END IF

      CALL invfft( psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )

#endif

      IF( .NOT. PRESENT( alpha ) ) THEN
        c = DBLE( psi )
      ELSE
        c = c + alpha * DBLE( psi )
      END IF

      DEALLOCATE(psi, STAT=ierr)
      IF( ierr /= 0 )  call errore(' PFWFFT ', ' deallocation of psi failed ' ,0)

      total_fft_time =  total_fft_time + ( cclock() - t1 )
      number_of_fft  =  number_of_fft  + 1

      RETURN
   END SUBROUTINE pinvfft2

!=---------------------------------------------------------------------==!
!
!
!     Wave Function FFTs high level Drivers
!
!
!=---------------------------------------------------------------------==!

   SUBROUTINE pw_fwfft( psi, c, ca )

     !   This subroutine COMPUTE on the wave function sub-grid :

     USE fft_types, ONLY: fft_dlay_descriptor
     USE fft_base, ONLY: dffts
     USE mp_global, ONLY: mpime

     COMPLEX(DP) :: C(:)
     COMPLEX(DP), OPTIONAL :: CA(:)
     COMPLEX(DP) :: psi(:)
     REAL(DP)  :: T1
     INTEGER :: ierr
     COMPLEX(DP), ALLOCATABLE :: psitmp(:)
     COMPLEX(DP), ALLOCATABLE :: zwrk(:) 

     t1 = cclock()

     IF ( first ) &
       CALL errore(' pw_fwfft 1 ',' fft not initialized ', 1 )

     IF ( SIZE( psi ) /= dffts%nnr ) &
       CALL errore( ' pw_fwfft 1 ', ' inconsistent array dimensions ', 1 )

#if defined __PARA


     ALLOCATE( zwrk( dffts%nsp( mpime + 1 ) * dffts%nr3x ) )

     CALL pc3fft_drv(psi(1), zwrk, -1, dffts, FFT_MODE_WAVE)

     IF( PRESENT( ca ) ) THEN
       CALL psi2c( zwrk, SIZE( zwrk ), c(1), ca(1), ngw, 2 )
     ELSE
       CALL psi2c( zwrk, SIZE( zwrk ), c(1), c(1), ngw, 0 )
     END IF

     DEALLOCATE( zwrk )

#else

     ALLOCATE( psitmp( SIZE( psi ) ) )
 
     psitmp = psi

     CALL fwfftw( psitmp, dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x )

     IF( PRESENT( ca ) ) THEN
       CALL psi2c( psitmp, dffts%nnr, c(1), ca(1), ngw, 2 )
     ELSE
       CALL psi2c( psitmp, dffts%nnr, c(1), c(1), ngw, 0 )
     END IF

     DEALLOCATE( psitmp )

#endif

     total_fft_wf_time =  total_fft_wf_time + ( cclock() - t1 )
     number_of_fft_wf  =  number_of_fft_wf  + 1

     RETURN
   END SUBROUTINE pw_fwfft

!
!==---------------------------------------------------------------------==!
!

   SUBROUTINE pw_invfft( psi, c, ca )

     !   This subroutine COMPUTE on the wave function sub-grid :

     USE fft_types, ONLY: fft_dlay_descriptor
     USE fft_base, ONLY: dffts
     USE mp_global, ONLY: mpime

     COMPLEX(DP), INTENT(IN) :: C(:)
     COMPLEX(DP), INTENT(IN), OPTIONAL :: CA(:)
     COMPLEX(DP) :: psi(:)
     COMPLEX(DP), ALLOCATABLE :: zwrk(:) 
     REAL(DP) :: T1

     IF ( first ) &
       CALL errore(' pw_invfft 1 ',' fft not initialized ', 1 )

     T1 = cclock()


     IF ( SIZE( psi ) /= dffts%nnr ) &
       CALL errore( ' pw_invfft 1 ', ' inconsistent array dimensions ', 1 )

#if defined __PARA

     ALLOCATE( zwrk( dffts%nsp( mpime + 1 ) * dffts%nr3x ) )

     IF( PRESENT( ca ) ) THEN
       CALL c2psi( zwrk, SIZE( zwrk ), c(1), ca(1), ngw, 2 )
     ELSE
       IF( tk ) THEN
         CALL c2psi( zwrk, SIZE( zwrk ), c(1), c(1), ngw, 0 )
       ELSE
         CALL c2psi( zwrk, SIZE( zwrk ), c(1), c(1), ngw, 1 )
       END IF
     END IF

     CALL pc3fft_drv(psi(1), zwrk, +1, dffts, FFT_MODE_WAVE)

     DEALLOCATE( zwrk )

#else
 
     IF( PRESENT( ca ) ) THEN
       CALL c2psi( psi, dffts%nnr, c(1), ca(1), ngw, 2 )
     ELSE
       IF( tk ) THEN
         CALL c2psi( psi, dffts%nnr, c(1), c(1), ngw, 0 )
       ELSE
         CALL c2psi( psi, dffts%nnr, c(1), c(1), ngw, 1 )
       END IF
     END IF

     CALL ivfftw( psi, dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x )

#endif

     total_fft_wf_time =  total_fft_wf_time + ( cclock() - T1 )
     number_of_fft_wf  =  number_of_fft_wf  + 1

     RETURN
   END SUBROUTINE pw_invfft



!==---------------------------------------------------------------------==!
   END MODULE fft
!==---------------------------------------------------------------------==!


!-----------------------------------------------------------------------
      subroutine invfft(f,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!-----------------------------------------------------------------------
! inverse fourier transform of potentials and charge density
! on the dense grid . On output, f is overwritten
!
      use fft_cp, only: cfft_cp
      use para_mod, only: dfftp
      use fft_scalar, only: cfft3d

      complex(8) f(nr1x*nr2x*nr3x)
      integer nr1,nr2,nr3,nr1x,nr2x,nr3x
      call start_clock( 'fft' )
#ifdef __PARA
      call cfft_cp(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,1,dfftp)
#else
# if defined __AIX || __FFTW || __COMPLIB || __SCSL
      call cfft3d(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,1)
# else
      call cfft3(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,1)
# endif
#endif
      call stop_clock( 'fft' )
!
      return
      end subroutine invfft

!-----------------------------------------------------------------------
      subroutine ivfftb(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3)
!-----------------------------------------------------------------------
! inverse fourier transform of Q functions (Vanderbilt pseudopotentials)
! on the  box grid . On output, f is overwritten
!
      use fft_scalar, only: cfft3d
      integer nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3
      complex(8) f(nr1bx*nr2bx*nr3bx)

!     in a parallel execution, not all processors calls this routine
!     then we should avoid clocks, otherwise the program hangs in print_clock
#if ! defined __PARA
      call start_clock( 'fftb' )
#endif

#ifdef __PARA
      call cfftpb(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3,1)
#else
# if defined __AIX || __FFTW || __COMPLIB || __SCSL
      call cfft3d(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,1)
# else
      call cfft3b(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,1)
# endif
#endif

#if ! defined __PARA
      call stop_clock( 'fftb' )
#endif
!
      return
      end subroutine ivfftb


!-----------------------------------------------------------------------
      subroutine ivffts(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!-----------------------------------------------------------------------
! inverse fourier transform of  potentials and charge density
! on the smooth grid . On output, f is overwritten
!
      use fft_cp, only: cfft_cp
      use para_mod, only: dffts
      use fft_scalar, only: cfft3d
      complex(8) f(nr1sx*nr2sx*nr3sx)
      integer nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx
      call start_clock( 'ffts' )
#ifdef __PARA
      call cfft_cp(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,1,dffts)
#else
# if defined __AIX || __FFTW || __COMPLIB || __SCSL
      call cfft3d(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,1)
# else
      call cfft3s(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,1)
# endif
#endif
      call stop_clock( 'ffts' )
!
      return
      end subroutine ivffts


!-----------------------------------------------------------------------
      subroutine ivfftw(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!-----------------------------------------------------------------------
! inverse fourier transform of wavefunctions
! on the smooth grid . On output, f is overwritten
!
      use fft_cp, only: cfft_cp
      use para_mod, only: dffts
      use fft_scalar, only: cfft3d, cfft3ds
      complex(8) f(nr1sx*nr2sx*nr3sx)
      integer nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx
      call start_clock('fftw')
#ifdef __PARA
      call cfft_cp(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,2,dffts)
#else
# if defined __AIX || __FFTW
      call cfft3ds(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,1,dffts%isind, dffts%iplw)
# elif defined __COMPLIB || __SCSL
      call cfft3d(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,1)
# else
      call cfft3s(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,1)
# endif
#endif
      call stop_clock('fftw')
!
      return
      end subroutine ivfftw


!-----------------------------------------------------------------------
      subroutine fwfft(f,nr1,nr2,nr3,nr1x,nr2x,nr3x)
!-----------------------------------------------------------------------
! forward fourier transform of potentials and charge density 
! on the dense grid . On output, f is overwritten
! 
      use fft_cp, only: cfft_cp
      use para_mod, only: dfftp
      use fft_scalar, only: cfft3d
      complex(8) f(nr1x*nr2x*nr3x)
      integer nr1,nr2,nr3,nr1x,nr2x,nr3x
      call start_clock( 'fft' )
#ifdef __PARA
      call cfft_cp(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,-1,dfftp)
#else 
# if defined __AIX || __FFTW || __COMPLIB || __SCSL
      call cfft3d(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,-1)
# else
      call cfft3(f,nr1,nr2,nr3,nr1x,nr2x,nr3x,-1)
# endif
#endif
      call stop_clock( 'fft' )
      return
      end subroutine fwfft


!-----------------------------------------------------------------------
      subroutine fwffts(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!-----------------------------------------------------------------------
! forward fourier transform of potentials and charge density
! on the smooth grid . On output, f is overwritten
!
      use fft_cp, only: cfft_cp
      use para_mod, only: dffts
      use fft_scalar, only: cfft3d
      complex(8) f(nr1sx*nr2sx*nr3sx)
      integer nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx
      call start_clock( 'ffts' )
#ifdef __PARA
      call cfft_cp(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,-1,dffts)
#else
# if defined __AIX || __FFTW || __COMPLIB || __SCSL
      call cfft3d(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,-1)
# else
      call cfft3s(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,-1)
# endif
#endif
      call stop_clock( 'ffts' )
      return
      end subroutine fwffts


!-----------------------------------------------------------------------
      subroutine fwfftw(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx)
!-----------------------------------------------------------------------
! forward fourier transform of potentials and charge density
! on the smooth grid . On output, f is overwritten
!
      use fft_cp, only: cfft_cp
      use para_mod, only: dffts
      use fft_scalar, only: cfft3d, cfft3ds
      complex(8) f(nr1sx*nr2sx*nr3sx)
      integer nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx
      call start_clock( 'fftw' )
#ifdef __PARA
      call cfft_cp(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,-2,dffts)
#else
# if defined __AIX || __FFTW 
      call cfft3ds(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,-1,dffts%isind, dffts%iplw)
# elif defined __COMPLIB || __SCSL
      call cfft3d(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,-1)
# else
      call cfft3s(f,nr1s,nr2s,nr3s,nr1sx,nr2sx,nr3sx,-1)
# endif
#endif
      call stop_clock( 'fftw' )
      return
      end subroutine fwfftw


!-----------------------------------------------------------------------
      subroutine ivfftbold(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3)
!-----------------------------------------------------------------------
! inverse fourier transform of Q functions (Vanderbilt pseudopotentials)
! on the  box grid . On output, f is overwritten
!
      use fft_scalar, only: cfft3d
      integer nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3
      complex(8) f(nr1bx*nr2bx*nr3bx)
      call start_clock(' ivfftbold ' )

#ifdef __PARA
      ! call cfft3d(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,1)  ! DEBUG
      call cfftpb(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,irb3,1)
#else
# if defined __AIX || __FFTW || __COMPLIB || __SCSL
      call cfft3d(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,1)
# else
      call cfft3b(f,nr1b,nr2b,nr3b,nr1bx,nr2bx,nr3bx,1)
# endif
#endif

      call stop_clock(' ivfftbold ' )
!
      return
      end subroutine ivfftbold

!-----------------------------------------------------------------------

     SUBROUTINE c2psi( psi, nnr, c, ca, ng, iflg )
       !
       use recvecs_indexes, only: nm, np
       use gvecs, only: nms, nps
       use kinds, only: DP

       implicit none

       complex(DP) :: psi(*), c(*), ca(*)
       integer, intent(in) :: nnr, ng, iflg

       complex(DP), parameter :: ci=(0.0d0,1.0d0)
       integer :: ig
       
         psi( 1 : nnr ) = 0.0d0

         !
         !  iflg "cases"
         !
         !  0, 10   Do not use gamma symmetry
         !
         !  1, 11   set psi using a wf with Gamma symmetry
         !
         !  2, 12   set psi combining two wf with Gamma symmetry
         !

         SELECT CASE ( iflg )
           !
           !  Case 0, 1 and 2  SMOOTH MESH
           !
           CASE ( 0 )
             !
             do ig = 1, ng
               psi( nps( ig ) ) = c( ig )
             end do
             !
           CASE ( 1 )
             !
             do ig = 1, ng
               psi( nms( ig ) ) = CONJG( c( ig ) )
               psi( nps( ig ) ) = c( ig )
             end do
             !
           CASE ( 2 )
             !
             do ig = 1, ng
               psi( nms( ig ) ) = CONJG( c( ig ) ) + ci * conjg( ca( ig ) )
               psi( nps( ig ) ) = c( ig ) + ci * ca( ig )
             end do

           !
           !  Case 10, 11 and 12  DENSE MESH
           !
           CASE ( 10 )
             !
             do ig = 1, ng
               psi( np( ig ) ) = c( ig )
             end do
             !
           CASE ( 11 )
             !
             do ig = 1, ng
               psi( nm( ig ) ) = CONJG( c( ig ) )
               psi( np( ig ) ) = c( ig )
             end do
             !
           CASE ( 12 )
             !
             do ig = 1, ng
               psi( nm( ig ) ) = CONJG( c( ig ) ) + ci * conjg( ca( ig ) )
               psi( np( ig ) ) = c( ig ) + ci * ca( ig )
             end do
             !
           CASE DEFAULT
             !
             CALL errore(" c2psi "," wrong value for iflg ", ABS( iflg ) )

         END SELECT

       return
     END SUBROUTINE c2psi

!-----------------------------------------------------------------------

     SUBROUTINE psi2c( psi, nnr, c, ca, ng, iflg )

       use recvecs_indexes, only: nm, np
       use gvecs, only: nms, nps
       use kinds, only: DP

       implicit none

       complex(DP) :: psi(*), c(*), ca(*)
       integer, intent(in) :: nnr, ng, iflg

       complex(DP), parameter :: ci=(0.0d0,1.0d0)
       integer :: ig

         !
         !  iflg "cases"
         !
         !  0, 10   Do not use gamma symmetry
         !
         !  1, 11   set psi using a wf with Gamma symmetry
         !
         !  2, 12   set psi combining two wf with Gamma symmetry
         !
       
         SELECT CASE ( iflg )

           !
           !  Case 0, 1 and 2  SMOOTH MESH
           !
           CASE ( 0 )
             !
             do ig = 1, ng
               c( ig ) = psi( nps( ig ) )
             end do
             !
           CASE ( 1 )
             !
             CALL errore(" psi2c "," wrong value for iflg ", 11 )
             !
           CASE ( 2 )
             !
             DO ig = 1, ng
               ca(ig) = psi( nms( ig ) )
               c (ig) = psi( nps( ig ) )
             END DO

           !
           !  Case 10, 11 and 12  DENSE MESH
           !
           CASE ( 10 )
             !
             do ig = 1, ng
               c( ig ) = psi( np( ig ) )
             end do
             !
           CASE ( 11 )
             !
             CALL errore(" psi2c "," wrong value for iflg ", 1 )
             !
           CASE ( 12 )
             !
             DO ig = 1, ng
               ca(ig) = psi( nm( ig ) )
               c (ig) = psi( np( ig ) )
             END DO

           CASE DEFAULT
             !
             CALL errore(" psi2c "," wrong value for iflg ", ABS( iflg ) )

         END SELECT
       
       return
     END SUBROUTINE psi2c
