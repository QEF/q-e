!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------==!
!
!
!     Parallel 3D FFT high level Driver
!     ( Charge density and Wave Functions )
!
!     Written and maintained by Carlo Cavazzoni
!     Last update Apr. 2009
!
!!=---------------------------------------------------------------------==!
!
MODULE fft_parallel
!
   IMPLICIT NONE
   SAVE

   INTEGER, PARAMETER :: DP = selected_real_kind(14,200)

   PRIVATE :: DP
!
CONTAINS
!
!  General purpose driver, including Task groups parallelization
!
!----------------------------------------------------------------------------
SUBROUTINE tg_cft3s( f, dfft, isgn, use_task_groups )
  !----------------------------------------------------------------------------
  !
  !! ... isgn = +-1 : parallel 3d fft for rho and for the potential
  !                  NOT IMPLEMENTED WITH TASK GROUPS
  !! ... isgn = +-2 : parallel 3d fft for wavefunctions
  !
  !! ... isgn = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  !! ...              fft along z using pencils        (cft_1z)
  !! ...              transpose across nodes           (fft_scatter)
  !! ...                 and reorder
  ! ...              fft along y (using planes) and x (cft_2xy)
  ! ... isgn = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  ! ...              fft along x and y(using planes)  (cft_2xy)
  ! ...              transpose across nodes           (fft_scatter)
  ! ...                 and reorder
  ! ...              fft along z using pencils        (cft_1z)
  !
  ! ...  The array "planes" signals whether a fft is needed along y :
  ! ...    planes(i)=0 : column f(i,*,*) empty , don't do fft along y
  ! ...    planes(i)=1 : column f(i,*,*) filled, fft along y needed
  ! ...  "empty" = no active components are present in f(i,*,*)
  ! ...            after (isgn>0) or before (isgn<0) the fft on z direction
  !
  ! ...  Note that if isgn=+/-1 (fft on rho and pot.) all fft's are needed
  ! ...  and all planes(i) are set to 1
  !
  ! This driver is based on code written by Stefano de Gironcoli for PWSCF.
  ! Task Group added by Costas Bekas, Oct. 2005, adapted from the CPMD code
  ! (Alessandro Curioni) and revised by Carlo Cavazzoni 2007.
  !
  USE fft_scalar, ONLY : cft_1z, cft_2xy
  USE scatter_mod,   ONLY : fft_scatter
  USE fft_types,  ONLY : fft_dlay_descriptor

  !
  IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  !
  COMPLEX(DP), INTENT(inout)    :: f( : )  ! array containing data to be transformed
  TYPE (fft_dlay_descriptor), INTENT(inout) :: dfft
                                           ! descriptor of fft data layout
  INTEGER, INTENT(in)           :: isgn    ! fft direction
  LOGICAL, OPTIONAL, INTENT(in) :: use_task_groups
                                           ! specify if you want to use task groups parallelization
  !
  ! the following ifdef prevents usage of directive in older ifort versions
  !
#if defined(__INTEL_COMPILER)
#if __INTEL_COMPILER  >= 1300
  ! the following is a workaround for Intel 12.1 bug
#if __INTEL_COMPILER  < 9999
!dir$ attributes align: 4096 :: yf, aux
#endif
#endif
#endif
  INTEGER                    :: me_p
  INTEGER                    :: n1, n2, n3, nx1, nx2, nx3
  COMPLEX(DP), ALLOCATABLE   :: yf(:), aux (:)
  INTEGER                    :: planes( dfft%nr1x )
  LOGICAL                    :: use_tg
  !
  !
  IF( present( use_task_groups ) ) THEN
     use_tg = use_task_groups
  ELSE
     use_tg = .false.
  ENDIF
  !
  IF( use_tg .and. .not. dfft%have_task_groups ) &
     CALL fftx_error__( ' tg_cft3s ', ' call requiring task groups for a descriptor without task groups ', 1 )
  !
  n1  = dfft%nr1
  n2  = dfft%nr2
  n3  = dfft%nr3
  nx1 = dfft%nr1x
  nx2 = dfft%nr2x
  nx3 = dfft%nr3x
  !
  IF( use_tg ) THEN
     ALLOCATE( aux( dfft%nogrp * dfft%tg_nnr ) )
     ALLOCATE( YF ( dfft%nogrp * dfft%tg_nnr ) )
  ELSE
     ALLOCATE( aux( dfft%nnr ) )
  ENDIF
  !
  me_p = dfft%mype + 1
  !
  IF ( isgn > 0 ) THEN
     !
     IF ( isgn /= 2 ) THEN
        !
        IF( use_tg ) &
           CALL fftx_error__( ' tg_cft3s ', ' task groups on large mesh not implemented ', 1 )
        !
        CALL cft_1z( f, dfft%nsp( me_p ), n3, nx3, isgn, aux )
        !
        planes = dfft%iplp
        !
     ELSE
        !
        IF( use_tg ) THEN
           CALL pack_group_sticks( f, yf, dfft )
           CALL cft_1z( yf, dfft%tg_nsw( me_p ), n3, nx3, isgn, aux )
        ELSE
           CALL cft_1z( f, dfft%nsw( me_p ), n3, nx3, isgn, aux )
        ENDIF
        !
        planes = dfft%iplw
        !
     ENDIF
     !
     CALL fw_scatter( isgn ) ! forward scatter from stick to planes
     !
     IF( use_tg ) THEN
        CALL cft_2xy( f, dfft%tg_npp( me_p ), n1, n2, nx1, nx2, isgn, planes )
     ELSE
        CALL cft_2xy( f, dfft%npp( me_p ), n1, n2, nx1, nx2, isgn, planes )
     ENDIF
     !
  ELSE
     !
     IF ( isgn /= -2 ) THEN
        !
        IF( use_tg ) &
           CALL fftx_error__( ' tg_cft3s ', ' task groups on large mesh not implemented ', 1 )
        !
        planes = dfft%iplp
        !
     ELSE
        !
        planes = dfft%iplw
        !
     ENDIF

     IF( use_tg ) THEN
        CALL cft_2xy( f, dfft%tg_npp( me_p ), n1, n2, nx1, nx2, isgn, planes )
     ELSE
        CALL cft_2xy( f, dfft%npp( me_p ), n1, n2, nx1, nx2, isgn, planes)
     ENDIF
     !
     CALL bw_scatter( isgn )
     !
     IF ( isgn /= -2 ) THEN
        !
        CALL cft_1z( aux, dfft%nsp( me_p ), n3, nx3, isgn, f )
         !
     ELSE
        !
        IF( use_tg ) THEN
           CALL cft_1z( aux, dfft%tg_nsw( me_p ), n3, nx3, isgn, yf )
           CALL unpack_group_sticks( yf, f, dfft )
        ELSE
           CALL cft_1z( aux, dfft%nsw( me_p ), n3, nx3, isgn, f )
        ENDIF
        !
     ENDIF
     !
  ENDIF
  !
  IF( use_tg ) THEN
     DEALLOCATE( yf )
  ENDIF

  DEALLOCATE( aux )
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE fw_scatter( iopt )

        !Transpose data for the 2-D FFT on the x-y plane
        !
        !NOGRP*dfft%nnr: The length of aux and f
        !nr3x: The length of each Z-stick
        !aux: input - output
        !f: working space
        !isgn: type of scatter
        !dfft%nsw(me) holds the number of Z-sticks proc. me has.
        !dfft%npp: number of planes per processor
        !
     !
     USE scatter_mod, ONLY: fft_scatter
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == 2 ) THEN
        !
        IF( use_tg ) THEN
           !
           CALL fft_scatter( dfft, aux, nx3, dfft%nogrp*dfft%tg_nnr, f, dfft%tg_nsw, dfft%tg_npp, iopt, use_tg )
           !
        ELSE
           !
           CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsw, dfft%npp, iopt )
           !
        ENDIF
        !
     ELSEIF( iopt == 1 ) THEN
        !
        CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsp, dfft%npp, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE fw_scatter

  !

  SUBROUTINE bw_scatter( iopt )
     !
     USE scatter_mod, ONLY: fft_scatter
     !
     INTEGER, INTENT(in) :: iopt
     !
     IF( iopt == -2 ) THEN
        !
        IF( use_tg ) THEN
           !
           CALL fft_scatter( dfft, aux, nx3, dfft%nogrp*dfft%tg_nnr, f, dfft%tg_nsw, dfft%tg_npp, iopt, use_tg )
           !
        ELSE
           !
           CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsw, dfft%npp, iopt )
           !
        ENDIF
        !
     ELSEIF( iopt == -1 ) THEN
        !
        CALL fft_scatter( dfft, aux, nx3, dfft%nnr, f, dfft%nsp, dfft%npp, iopt )
        !
     ENDIF
     !
     RETURN
  END SUBROUTINE bw_scatter
  !
END SUBROUTINE tg_cft3s
!
!
!
!----------------------------------------------------------------------------
SUBROUTINE fw_tg_cft3_z( f_in, dfft, f_out )
  !----------------------------------------------------------------------------
  !
  USE fft_scalar, ONLY : cft_1z
  USE fft_types,  ONLY : fft_dlay_descriptor
  !
  IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  !
  COMPLEX(DP), INTENT(inout)    :: f_in( : )  ! INPUT array containing data to be transformed
  COMPLEX(DP), INTENT(inout)   :: f_out (:)  ! OUTPUT
  TYPE (fft_dlay_descriptor), INTENT(inout) :: dfft ! descriptor of fft data layout
  !
  CALL cft_1z( f_in, dfft%tg_nsw( dfft%mype + 1 ), dfft%nr3, dfft%nr3x, 2, f_out )
  !
END SUBROUTINE fw_tg_cft3_z
!
!----------------------------------------------------------------------------
SUBROUTINE bw_tg_cft3_z( f_out, dfft, f_in )
  !----------------------------------------------------------------------------
  !
  USE fft_scalar, ONLY : cft_1z
  USE fft_types,  ONLY : fft_dlay_descriptor
  !
  IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  !
  COMPLEX(DP), INTENT(inout)    :: f_out( : ) ! OUTPUT
  COMPLEX(DP), INTENT(inout)   :: f_in (:) ! INPUT array containing data to be transformed
  TYPE (fft_dlay_descriptor), INTENT(inout) :: dfft ! descriptor of fft data layout
  !
  CALL cft_1z( f_in, dfft%tg_nsw( dfft%mype + 1 ), dfft%nr3, dfft%nr3x, -2, f_out )
  !
END SUBROUTINE bw_tg_cft3_z
!
!----------------------------------------------------------------------------
SUBROUTINE fw_tg_cft3_scatter( f, dfft, aux )
  !----------------------------------------------------------------------------
  !
  USE scatter_mod,   ONLY : fft_scatter
  USE fft_types,  ONLY : fft_dlay_descriptor
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout)    :: f( : ), aux( : )  ! array containing data to be transformed
  TYPE (fft_dlay_descriptor), INTENT(inout) :: dfft     ! descriptor of fft data layout
  !
  CALL fft_scatter( dfft, aux, dfft%nr3x, dfft%nogrp*dfft%tg_nnr, f, dfft%tg_nsw, dfft%tg_npp, 2, .true. )
  !
END SUBROUTINE fw_tg_cft3_scatter
!
!----------------------------------------------------------------------------
SUBROUTINE bw_tg_cft3_scatter( f, dfft, aux )
  !----------------------------------------------------------------------------
  !
  USE scatter_mod,   ONLY : fft_scatter
  USE fft_types,  ONLY : fft_dlay_descriptor
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout)    :: f( : ), aux( : )  ! array containing data to be transformed
  TYPE (fft_dlay_descriptor), INTENT(inout) :: dfft     ! descriptor of fft data layout
  !
  CALL fft_scatter( dfft, aux, dfft%nr3x, dfft%nogrp*dfft%tg_nnr, f, dfft%tg_nsw, dfft%tg_npp, -2, .true. )
  !
END SUBROUTINE bw_tg_cft3_scatter
!
!----------------------------------------------------------------------------
SUBROUTINE fw_tg_cft3_xy( f, dfft )
  !----------------------------------------------------------------------------
  !
  USE fft_scalar, ONLY : cft_2xy
  USE fft_types,  ONLY : fft_dlay_descriptor
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout)    :: f( : ) ! INPUT/OUTPUT array containing data to be transformed
  TYPE (fft_dlay_descriptor), INTENT(inout) :: dfft ! descriptor of fft data layout
  INTEGER                    :: planes( dfft%nr1x )
  !
  planes = dfft%iplw
  CALL cft_2xy( f, dfft%tg_npp( dfft%mype + 1 ), dfft%nr1, dfft%nr2, dfft%nr1x, dfft%nr2x, 2, planes )
  !
END SUBROUTINE fw_tg_cft3_xy
!
!----------------------------------------------------------------------------
SUBROUTINE bw_tg_cft3_xy( f, dfft )
  !----------------------------------------------------------------------------
  !
  USE fft_scalar, ONLY : cft_2xy
  USE fft_types,  ONLY : fft_dlay_descriptor
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout)    :: f( : ) ! INPUT/OUTPUT  array containing data to be transformed
  TYPE (fft_dlay_descriptor), INTENT(inout) :: dfft ! descriptor of fft data layout
  INTEGER                    :: planes( dfft%nr1x )
  !
  planes = dfft%iplw
  CALL cft_2xy( f, dfft%tg_npp( dfft%mype + 1 ), dfft%nr1, dfft%nr2, dfft%nr1x, dfft%nr2x, -2, planes )
  !
END SUBROUTINE bw_tg_cft3_xy


!----------------------------------------------------------------------------
  SUBROUTINE pack_group_sticks( f, yf, dfft )

     USE fft_types,  ONLY : fft_dlay_descriptor

     IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif

     COMPLEX(DP), INTENT(in)    :: f( : )  ! array containing all bands, and gvecs distributed across processors
     COMPLEX(DP), INTENT(out)    :: yf( : )  ! array containing bands collected into task groups
     TYPE (fft_dlay_descriptor), INTENT(inout) :: dfft
     INTEGER                     :: ierr
     !
     IF( dfft%tg_rdsp(dfft%nogrp) + dfft%tg_rcv(dfft%nogrp) > size( yf ) ) THEN
        CALL fftx_error__( 'pack_group_sticks' , ' inconsistent size ', 1 )
     ENDIF
     IF( dfft%tg_psdsp(dfft%nogrp) + dfft%tg_snd(dfft%nogrp) > size( f ) ) THEN
        CALL fftx_error__( 'pack_group_sticks', ' inconsistent size ', 2 )
     ENDIF

     CALL start_clock( 'ALLTOALL' )
     !
     !  Collect all the sticks of the different states,
     !  in "yf" processors will have all the sticks of the OGRP

#if defined(__MPI)

     CALL MPI_ALLTOALLV( f(1), dfft%tg_snd, dfft%tg_psdsp, MPI_DOUBLE_COMPLEX, yf(1), dfft%tg_rcv, &
      &                     dfft%tg_rdsp, MPI_DOUBLE_COMPLEX, dfft%ogrp_comm, IERR)
     IF( ierr /= 0 ) THEN
        CALL fftx_error__( 'pack_group_sticks', ' alltoall error 1 ', abs(ierr) )
     ENDIF

#else

     IF( dfft%tg_rcv(dfft%nogrp) /= dfft%tg_snd(dfft%nogrp) ) THEN
        CALL fftx_error__( 'pack_group_sticks', ' inconsistent size ', 3 )
     ENDIF

     yf( 1 : dfft%tg_rcv(dfft%nogrp) ) =  f( 1 : dfft%tg_snd(dfft%nogrp) )

#endif

     CALL stop_clock( 'ALLTOALL' )
     !
     !YF Contains all ( ~ NOGRP*dfft%nsw(me) ) Z-sticks
     !
     RETURN
  END SUBROUTINE pack_group_sticks

  SUBROUTINE pack_group_sticks_i( f, yf, dfft, req )

     USE fft_types,  ONLY : fft_dlay_descriptor

     IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif

     COMPLEX(DP), INTENT(in)    :: f( : )  ! array containing all bands, and gvecs distributed across processors
     COMPLEX(DP), INTENT(out)    :: yf( : )  ! array containing bands collected into task groups
     TYPE (fft_dlay_descriptor), INTENT(inout) :: dfft
     INTEGER, INTENT(OUT) :: req
     INTEGER :: ierr
     !
     IF( dfft%tg_rdsp(dfft%nogrp) + dfft%tg_rcv(dfft%nogrp) > size( yf ) ) THEN
        CALL fftx_error__( 'pack_group_sticks' , ' inconsistent size ', 1 )
     ENDIF
     IF( dfft%tg_psdsp(dfft%nogrp) + dfft%tg_snd(dfft%nogrp) > size( f ) ) THEN
        CALL fftx_error__( 'pack_group_sticks', ' inconsistent size ', 2 )
     ENDIF

     CALL start_clock( 'ALLTOALL' )
     !
     !  Collect all the sticks of the different states,
     !  in "yf" processors will have all the sticks of the OGRP

#if defined(__MPI)

     CALL MPI_IALLTOALLV( f(1), dfft%tg_snd, dfft%tg_psdsp, MPI_DOUBLE_COMPLEX, yf(1), dfft%tg_rcv, &
      &                     dfft%tg_rdsp, MPI_DOUBLE_COMPLEX, dfft%ogrp_comm, req, IERR)
     IF( ierr /= 0 ) THEN
        CALL fftx_error__( 'pack_group_sticks', ' alltoall error 1 ', abs(ierr) )
     ENDIF

#else

     IF( dfft%tg_rcv(dfft%nogrp) /= dfft%tg_snd(dfft%nogrp) ) THEN
        CALL fftx_error__( 'pack_group_sticks', ' inconsistent size ', 3 )
     ENDIF

     yf( 1 : dfft%tg_rcv(dfft%nogrp) ) =  f( 1 : dfft%tg_snd(dfft%nogrp) )

#endif

     CALL stop_clock( 'ALLTOALL' )
     !
     !YF Contains all ( ~ NOGRP*dfft%nsw(me) ) Z-sticks
     !
     RETURN
  END SUBROUTINE pack_group_sticks_i

  !
  SUBROUTINE unpack_group_sticks( yf, f, dfft )

     USE fft_types,  ONLY : fft_dlay_descriptor

     IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif

     COMPLEX(DP), INTENT(out)    :: f( : )  ! array containing all bands, and gvecs distributed across processors
     COMPLEX(DP), INTENT(in)    :: yf( : )  ! array containing bands collected into task groups
     TYPE (fft_dlay_descriptor), INTENT(inout) :: dfft
     !
     !  Bring pencils back to their original distribution
     !
     INTEGER                     :: ierr
     !
     IF( dfft%tg_usdsp(dfft%nogrp) + dfft%tg_snd(dfft%nogrp) > size( f ) ) THEN
        CALL fftx_error__( 'unpack_group_sticks', ' inconsistent size ', 3 )
     ENDIF
     IF( dfft%tg_rdsp(dfft%nogrp) + dfft%tg_rcv(dfft%nogrp) > size( yf ) ) THEN
        CALL fftx_error__( 'unpack_group_sticks', ' inconsistent size ', 4 )
     ENDIF

     CALL start_clock( 'ALLTOALL' )

#if defined(__MPI)
     CALL MPI_Alltoallv( yf(1), &
          dfft%tg_rcv, dfft%tg_rdsp, MPI_DOUBLE_COMPLEX, f(1), &
          dfft%tg_snd, dfft%tg_usdsp, MPI_DOUBLE_COMPLEX, dfft%ogrp_comm, IERR)
     IF( ierr /= 0 ) THEN
        CALL fftx_error__( 'unpack_group_sticks', ' alltoall error 2 ', abs(ierr) )
     ENDIF
#endif

     CALL stop_clock( 'ALLTOALL' )

     RETURN
  END SUBROUTINE unpack_group_sticks

  SUBROUTINE unpack_group_sticks_i( yf, f, dfft, req )

     USE fft_types,  ONLY : fft_dlay_descriptor

     IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
     INTEGER :: req
#endif

     COMPLEX(DP), INTENT(out)    :: f( : )  ! array containing all bands, and gvecs distributed across processors
     COMPLEX(DP), INTENT(in)    :: yf( : )  ! array containing bands collected into task groups
     TYPE (fft_dlay_descriptor), INTENT(inout) :: dfft
     !
     !  Bring pencils back to their original distribution
     !
     INTEGER                     :: ierr
     !
     IF( dfft%tg_usdsp(dfft%nogrp) + dfft%tg_snd(dfft%nogrp) > size( f ) ) THEN
        CALL fftx_error__( 'unpack_group_sticks', ' inconsistent size ', 3 )
     ENDIF
     IF( dfft%tg_rdsp(dfft%nogrp) + dfft%tg_rcv(dfft%nogrp) > size( yf ) ) THEN
        CALL fftx_error__( 'unpack_group_sticks', ' inconsistent size ', 4 )
     ENDIF

     CALL start_clock( 'ALLTOALL' )

#if defined(__MPI)
     CALL MPI_IAlltoallv( yf(1), &
          dfft%tg_rcv, dfft%tg_rdsp, MPI_DOUBLE_COMPLEX, f(1), &
          dfft%tg_snd, dfft%tg_usdsp, MPI_DOUBLE_COMPLEX, dfft%ogrp_comm, req, IERR)
     IF( ierr /= 0 ) THEN
        CALL fftx_error__( 'unpack_group_sticks', ' alltoall error 2 ', abs(ierr) )
     ENDIF
#endif

     CALL stop_clock( 'ALLTOALL' )

     RETURN
  END SUBROUTINE unpack_group_sticks_i

SUBROUTINE tg_gather( dffts, v, tg_v )
   !
   USE fft_types,      ONLY : fft_dlay_descriptor

   ! T.G.
   ! NOGRP:      Number of processors per orbital task group

   IMPLICIT NONE
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif

   TYPE(fft_dlay_descriptor), INTENT(inout) :: dffts

   REAL(DP) :: v(:)
   REAL(DP) :: tg_v(:)

   INTEGER :: nsiz, i, ierr, nsiz_tg
   INTEGER :: recv_cnt( dffts%nogrp ), recv_displ( dffts%nogrp )

   nsiz_tg = dffts%tg_nnr * dffts%nogrp

   IF( size( tg_v ) < nsiz_tg ) &
      CALL fftx_error__( ' tg_gather ', ' tg_v too small ', ( nsiz_tg - size( tg_v ) ) )

   nsiz = dffts%npp( dffts%mype+1 ) * dffts%nr1x * dffts%nr2x

   IF( size( v ) < nsiz ) &
      CALL fftx_error__( ' tg_gather ', ' v too small ',  ( nsiz - size( v ) ) )

   !
   !  The potential in v is distributed across all processors
   !  We need to redistribute it so that it is completely contained in the
   !  processors of an orbital TASK-GROUP
   !
   recv_cnt(1)   = dffts%npp( dffts%nolist(1) + 1 ) * dffts%nr1x * dffts%nr2x
   recv_displ(1) = 0
   DO i = 2, dffts%nogrp
      recv_cnt(i) = dffts%npp( dffts%nolist(i) + 1 ) * dffts%nr1x * dffts%nr2x
      recv_displ(i) = recv_displ(i-1) + recv_cnt(i-1)
   ENDDO

   ! clean only elements that will not be overwritten
   !
   DO i = recv_displ(dffts%nogrp) + recv_cnt( dffts%nogrp ) + 1, size( tg_v )
      tg_v( i ) = 0.0d0
   ENDDO

#if defined(__MPI)

   CALL MPI_Allgatherv( v(1), nsiz, MPI_DOUBLE_PRECISION, &
        tg_v(1), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, dffts%ogrp_comm, IERR)

   IF( ierr /= 0 ) &
      CALL fftx_error__( ' tg_gather ', ' MPI_Allgatherv ', abs( ierr ) )

#endif

END SUBROUTINE tg_gather
!
!  Complex version of previous routine
!
SUBROUTINE tg_cgather( dffts, v, tg_v )
   !
   USE fft_types,      ONLY : fft_dlay_descriptor

   ! T.G.
   ! NOGRP:      Number of processors per orbital task group

   IMPLICIT NONE
#if defined(__MPI)
   INCLUDE 'mpif.h'
#endif

   TYPE(fft_dlay_descriptor), INTENT(inout) :: dffts

   COMPLEX(DP) :: v(:)
   COMPLEX(DP) :: tg_v(:)

   INTEGER :: nsiz, i, ierr, nsiz_tg
   INTEGER :: recv_cnt( dffts%nogrp ), recv_displ( dffts%nogrp )

   nsiz_tg = dffts%tg_nnr * dffts%nogrp

   IF( size( tg_v ) < nsiz_tg ) &
      CALL fftx_error__( ' tg_gather ', ' tg_v too small ', ( nsiz_tg - size( tg_v ) ) )

   nsiz = dffts%npp( dffts%mype+1 ) * dffts%nr1x * dffts%nr2x

   IF( size( v ) < nsiz ) &
      CALL fftx_error__( ' tg_gather ', ' v too small ',  ( nsiz - size( v ) ) )

   !
   !  The potential in v is distributed across all processors
   !  We need to redistribute it so that it is completely contained in the
   !  processors of an orbital TASK-GROUP
   !
   recv_cnt(1)   = dffts%npp( dffts%nolist(1) + 1 ) * dffts%nr1x * dffts%nr2x
   recv_displ(1) = 0
   DO i = 2, dffts%nogrp
      recv_cnt(i) = dffts%npp( dffts%nolist(i) + 1 ) * dffts%nr1x * dffts%nr2x
      recv_displ(i) = recv_displ(i-1) + recv_cnt(i-1)
   ENDDO

   ! clean only elements that will not be overwritten
   !
   DO i = recv_displ(dffts%nogrp) + recv_cnt( dffts%nogrp ) + 1, size( tg_v )
      tg_v( i ) = (0.0d0,0.0d0)
   ENDDO
   !
   ! The quantities are complex, multiply the cunters by 2 and gather
   ! real numbers
   !
   nsiz = 2 * nsiz
   recv_cnt = 2 * recv_cnt
   recv_displ = 2 * recv_displ
#if defined(__MPI)

   CALL MPI_Allgatherv( v(1), nsiz, MPI_DOUBLE_PRECISION, &
        tg_v(1), recv_cnt, recv_displ, MPI_DOUBLE_PRECISION, dffts%ogrp_comm, IERR)

   IF( ierr /= 0 ) &
      CALL fftx_error__( ' tg_cgather ', ' MPI_Allgatherv ', abs( ierr ) )

#endif

END SUBROUTINE tg_cgather
END MODULE fft_parallel
